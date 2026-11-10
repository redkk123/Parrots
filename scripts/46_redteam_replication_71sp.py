#!/usr/bin/env python3
"""
Red-Team Audit Question 8: Holdout Replication on the 13 New Species (71-Species Extension)
and Out-of-Sample Holdout Evaluation of Key Phase 3 Discovery Genes (PTEN, MAD2L2, HSPA8, etc.).
"""
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from lib.phylo import Tree  # noqa: E402

OUT = ROOT / "results/redteam_audit"
G71 = ROOT / "results/phase2_genes_71sp"
OUT.mkdir(parents=True, exist_ok=True)
G71.mkdir(parents=True, exist_ok=True)

CODE = {
    a + b + c: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"[16 * i + 4 * j + k]
    for i, a in enumerate("TCAG") for j, b in enumerate("TCAG") for k, c in enumerate("TCAG")
}


def translate_codon(c):
    return CODE.get(c.upper(), "X") if "-" not in c and "N" not in c.upper() else "-"


def translate_seq(s):
    return "".join(translate_codon(s[i:i + 3]) for i in range(0, len(s) - 2, 3))


def parse_fasta(path):
    seqs = {}
    name, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)
            name, buf = line[1:].strip().split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def main():
    print("[RedTeam Q8] Preparing 71-species traits, master tree, and AA alignments...", flush=True)
    spec0b = importlib.util.spec_from_file_location("p0b", ROOT / "scripts/05_phase0b_analysis.py")
    p0b = importlib.util.module_from_spec(spec0b)
    spec0b.loader.exec_module(p0b)

    spec_sel = importlib.util.spec_from_file_location("sel", ROOT / "scripts/26_select_genomes.py")
    sel = importlib.util.module_from_spec(spec_sel)
    spec_sel.loader.exec_module(sel)

    full = Tree.from_newick((ROOT / "results/phase0/supertree_burgio2019_full.nwk").read_text())
    ph = pd.read_csv(ROOT / "results/phase0/phenotypes.tsv", sep="\t")
    ph = ph[ph.tree_tip.isin(set(full.tip_labels()))]
    b = ph.dropna(subset=["log_brain_mass", "log_body_mass"])
    _, r_b, _ = p0b.subset_fit(full, b, "log_brain_mass", ["log_body_mass"])

    m = pd.read_csv(ROOT / "results/phase0b/species_data_matrix_with_residuals.tsv", sep="\t")
    t = m[["species", "tree_tip", "rel_log_life_expectancy"]].copy()
    t["rel_log_brain"] = t.tree_tip.map(r_b)
    both = t.dropna(subset=["rel_log_brain", "rel_log_life_expectancy"])
    _, r_c, _ = p0b.subset_fit(full, both, "rel_log_brain", ["rel_log_life_expectancy"])
    t["rel_brain_given_le"] = t.tree_tip.map(r_c)

    sg71 = pd.read_csv(ROOT / "results/phase1_audit/species_genomes.tsv", sep="\t")
    keep71 = set(sg71.species[sg71.status.str.startswith("included")])
    sg58 = pd.read_csv(ROOT / "results/phase1_audit/analysis_species.tsv", sep="\t")
    keep58 = set(sg58.species[sg58.status.str.startswith("included")])

    t71 = t[t.species.isin(keep71)].copy()
    t71["species_tip"] = t71.species.str.replace(" ", "_")
    t71["clade"] = ["Psittaculidae" if sel.group(s.split()[0]) != "other_Psittaciformes" else "rest" for s in t71.species]
    t71["cohort"] = ["discovery_58" if s in keep58 else "holdout_13" for s in t71.species]
    t71.to_csv(G71 / "rer_traits_71sp.tsv", sep="\t", index=False, float_format="%.5g")
    t71[t71["cohort"] == "holdout_13"].to_csv(OUT / "8_holdout_13_species_list.tsv", sep="\t", index=False, float_format="%.5g")

    tip_of = dict(zip(m.species, m.tree_tip))
    tips71 = {tip_of[s]: s.replace(" ", "_") for s in t71.species if s in tip_of}
    pr71 = full.prune(set(tips71))
    pr71.label = [tips71.get(x, x) for x in pr71.label]
    (G71 / "master_burgio2019_71sp.nwk").write_text(pr71.to_newick())

    # Translate all 97 codon alignments in results/phase2_genes_71sp/alignments/ to AA
    aa_dir = G71 / "aa"
    aa_dir.mkdir(parents=True, exist_ok=True)
    for cf in sorted((G71 / "alignments").glob("*.codon.fasta")):
        g = cf.name.split(".")[0]
        seqs = parse_fasta(cf)
        with open(aa_dir / f"{g}.fa", "w") as f:
            for sp, cds in seqs.items():
                f.write(f">{sp}\n{translate_seq(cds)}\n")

    # Write and run R script for 71-species panel RERconverge + 58-species 45/13 holdout cross-validation
    r_code = """
suppressPackageStartupMessages({
  library(Rcpp)
  library(RERconverge)
  library(TreeTools)
  library(impute)
  library(ape)
  library(phangorn)
  library(parallel)
})

out_dir <- "results/redteam_audit"
g71_dir <- "results/phase2_genes_71sp"
ncores <- 4

mtree71 <- ape::multi2di(ape::read.tree(file.path(g71_dir, "master_burgio2019_71sp.nwk")), random = FALSE)
gt71_file <- file.path(g71_dir, "genetrees_burgio2019_71sp.txt")

if (!file.exists(gt71_file)) {
  alnfiles <- sort(list.files(path = file.path(g71_dir, "aa"), pattern = "\\\\.fa$", full.names = TRUE))
  cat("[RedTeam Q8 R] Estimating LG branch lengths for", length(alnfiles), "panel genes across 71 species...\\n")
  est_one <- function(f) {
    g <- sub("\\\\.fa$", "", basename(f))
    tryCatch({
      aln <- phangorn::read.phyDat(f, type = "AA", format = "fasta")
      inboth <- intersect(names(aln), mtree71$tip.label)
      if (length(inboth) < 10) return(NULL)
      gt <- ape::unroot(ape::multi2di(ape::keep.tip(mtree71, inboth), random = FALSE))
      aln <- subset(aln, subset = inboth)
      gt$edge.length <- rep(0.01, nrow(gt$edge))
      fit <- phangorn::optim.pml(
        phangorn::pml(gt, aln, model = "LG", k = 1),
        optEdge = TRUE, rearrangement = "none", model = "LG",
        control = phangorn::pml.control(trace = 0, epsilon = 1e-4, maxit = 5)
      )
      paste(g, ape::write.tree(fit$tree), sep = "\\t")
    }, error = function(e) NULL)
  }
  lines <- unlist(parallel::mclapply(alnfiles, est_one, mc.cores = ncores))
  writeLines(lines, gt71_file)
}

trees71 <- readTrees(gt71_file)
t71 <- read.delim(file.path(g71_dir, "rer_traits_71sp.tsv"), stringsAsFactors = FALSE)

tv_le_71 <- setNames(t71$rel_log_life_expectancy, t71$species_tip); tv_le_71 <- tv_le_71[!is.na(tv_le_71)]
tv_br_71 <- setNames(t71$rel_log_brain, t71$species_tip); tv_br_71 <- tv_br_71[!is.na(tv_br_71)]

sp58 <- t71$species_tip[t71$cohort == "discovery_58"]
sp13 <- t71$species_tip[t71$cohort == "holdout_13"]

pdf(NULL)
rers_71_le <- getAllResiduals(trees71, useSpecies = names(tv_le_71), transform = "sqrt", use.weights = TRUE, norm = "scale", min.sp = 10, min.valid = 20)
rers_71_br <- getAllResiduals(trees71, useSpecies = names(tv_br_71), transform = "sqrt", use.weights = TRUE, norm = "scale", min.sp = 10, min.valid = 20)

res_le_71 <- correlateWithContinuousPhenotype(rers_71_le, char2Paths(tv_le_71, trees71), min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)
res_br_71 <- correlateWithContinuousPhenotype(rers_71_br, char2Paths(tv_br_71, trees71), min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)

# Discovery 58 subset on the same panel trees
tv_le_58 <- tv_le_71[names(tv_le_71) %in% sp58]
tv_br_58 <- tv_br_71[names(tv_br_71) %in% sp58]
res_le_58 <- correlateWithContinuousPhenotype(rers_71_le, char2Paths(tv_le_58, trees71), min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)
res_br_58 <- correlateWithContinuousPhenotype(rers_71_br, char2Paths(tv_br_58, trees71), min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)

# Holdout 13 subset alone on the same panel trees (min.sp = 5 since holdout has 13 species total)
tv_le_13 <- tv_le_71[names(tv_le_71) %in% sp13]
tv_br_13 <- tv_br_71[names(tv_br_71) %in% sp13]
res_le_13 <- correlateWithContinuousPhenotype(rers_71_le, char2Paths(tv_le_13, trees71), min.sp = 5, winsorizeRER = 3, winsorizetrait = 3)
res_br_13 <- correlateWithContinuousPhenotype(rers_71_br, char2Paths(tv_br_13, trees71), min.sp = 5, winsorizeRER = 3, winsorizetrait = 3)

panel_rep <- data.frame(
  gene = rownames(res_le_71),
  Rho_LE_58sp = res_le_58[rownames(res_le_71), "Rho"],
  P_LE_58sp = res_le_58[rownames(res_le_71), "P"],
  Rho_LE_holdout13 = res_le_13[rownames(res_le_71), "Rho"],
  N_branches_LE_holdout13 = res_le_13[rownames(res_le_71), "N"],
  P_LE_holdout13 = res_le_13[rownames(res_le_71), "P"],
  Rho_LE_71sp = res_le_71[, "Rho"],
  N_branches_LE_71sp = res_le_71[, "N"],
  P_LE_71sp = res_le_71[, "P"],
  Rho_Brain_58sp = res_br_58[rownames(res_le_71), "Rho"],
  P_Brain_58sp = res_br_58[rownames(res_le_71), "P"],
  Rho_Brain_holdout13 = res_br_13[rownames(res_le_71), "Rho"],
  N_branches_Brain_holdout13 = res_br_13[rownames(res_le_71), "N"],
  P_Brain_holdout13 = res_br_13[rownames(res_le_71), "P"],
  Rho_Brain_71sp = res_br_71[rownames(res_le_71), "Rho"],
  P_Brain_71sp = res_br_71[rownames(res_le_71), "P"]
)
panel_rep <- panel_rep[order(panel_rep$P_LE_71sp), ]
write.table(panel_rep, file.path(out_dir, "8_holdout13_and_71sp_panel_replication.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE)

# Now run 30 clade-stratified 45-train / 13-test holdout splits on the 7,756-gene matrix to evaluate
# out-of-sample replication for PTEN, MAD2L2, HSPA8, CTNNB1, PSMA4, SSU72, TOMM20, ANAPC11, etc.
cat("[RedTeam Q8 R] Running 30 clade-stratified 45-train / 13-test holdout splits on 7,756 genes...\\n")
trees58 <- readTrees("results/phase3_genome/genetrees_burgio2019.txt")
traits58 <- read.delim("results/phase3_genome/rer_traits.tsv", stringsAsFactors = FALSE)
tv58_le <- setNames(traits58$rel_log_life_expectancy, traits58$species_tip); tv58_le <- tv58_le[!is.na(tv58_le)]
tv58_br <- setNames(traits58$rel_log_brain, traits58$species_tip); tv58_br <- tv58_br[!is.na(tv58_br)]
rers58 <- getAllResiduals(trees58, useSpecies = names(tv58_le), transform = "sqrt", use.weights = TRUE, norm = "scale", min.sp = 10, min.valid = 20)

key_genes <- c("PTEN", "MAD2L2", "HSPA8", "CTNNB1", "CTNNBIP1", "PSMA4", "SSU72", "TOMM20", "ANAPC11", "RTN1", "SOX14", "NCALD", "CBLN4")
set.seed(20261110)
psl_sp <- traits58$species_tip[traits58$clade == "Psittaculidae"]
rest_sp <- traits58$species_tip[traits58$clade == "rest"]

n_splits <- 30
train_le_mat <- matrix(NA_real_, nrow = length(key_genes), ncol = n_splits, dimnames = list(key_genes, NULL))
test_le_mat <- matrix(NA_real_, nrow = length(key_genes), ncol = n_splits, dimnames = list(key_genes, NULL))
test_br_mat <- matrix(NA_real_, nrow = length(key_genes), ncol = n_splits, dimnames = list(key_genes, NULL))
coupling_test13 <- numeric(n_splits)

for (s in 1:n_splits) {
  test_sp <- c(sample(psl_sp, 3), sample(rest_sp, 10))
  train_sp <- setdiff(traits58$species_tip, test_sp)
  p_tr_le <- suppressMessages(char2Paths(tv58_le[names(tv58_le) %in% train_sp], trees58))
  p_te_le <- suppressMessages(char2Paths(tv58_le[names(tv58_le) %in% test_sp], trees58))
  p_te_br <- suppressMessages(char2Paths(tv58_br[names(tv58_br) %in% test_sp], trees58))
  r_tr_le <- correlateWithContinuousPhenotype(rers58, p_tr_le, min.sp = 8, winsorizeRER = 3, winsorizetrait = 3)
  r_te_le <- correlateWithContinuousPhenotype(rers58, p_te_le, min.sp = 5, winsorizeRER = 3, winsorizetrait = 3)
  r_te_br <- correlateWithContinuousPhenotype(rers58, p_te_br, min.sp = 5, winsorizeRER = 3, winsorizetrait = 3)
  train_le_mat[, s] <- r_tr_le[key_genes, "Rho"]
  test_le_mat[, s] <- r_te_le[key_genes, "Rho"]
  test_br_mat[, s] <- r_te_br[key_genes, "Rho"]
  coupling_test13[s] <- cor(r_te_le$Rho, r_te_br$Rho, use = "complete.obs")
}

orig_le <- read.delim("results/phase3_genome/rer_burgio2019_rel_log_life_expectancy.tsv", stringsAsFactors = FALSE)
rownames(orig_le) <- orig_le$gene
orig_br <- read.delim("results/phase3_genome/rer_burgio2019_rel_log_brain.tsv", stringsAsFactors = FALSE)
rownames(orig_br) <- orig_br$gene

cv_df <- data.frame(
  gene = key_genes,
  full58_Rho_LE = orig_le[key_genes, "Rho"],
  full58_P_LE = orig_le[key_genes, "P"],
  mean_train45_Rho_LE = rowMeans(train_le_mat, na.rm = TRUE),
  mean_holdout13_Rho_LE = rowMeans(test_le_mat, na.rm = TRUE),
  sd_holdout13_Rho_LE = apply(test_le_mat, 1, sd, na.rm = TRUE),
  sign_concordance_holdout13_LE = rowMeans(sign(test_le_mat) == sign(orig_le[key_genes, "Rho"]), na.rm = TRUE),
  full58_Rho_Brain = orig_br[key_genes, "Rho"],
  mean_holdout13_Rho_Brain = rowMeans(test_br_mat, na.rm = TRUE),
  sign_concordance_holdout13_Brain = rowMeans(sign(test_br_mat) == sign(orig_br[key_genes, "Rho"]), na.rm = TRUE)
)
write.table(cv_df, file.path(out_dir, "8_holdout_cv_key_genes.tsv"), sep = "\\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[RedTeam Q8 R] Done! Mean coupling in 13-sp holdout sets = %.4f (SD=%.4f); Cor(Rho_LE_58, Rho_LE_71) across 86 panel genes = %.4f\\n",
            mean(coupling_test13, na.rm=TRUE), sd(coupling_test13, na.rm=TRUE),
            cor(panel_rep$Rho_LE_58sp, panel_rep$Rho_LE_71sp, use="complete.obs")))
"""
    r_script_path = ROOT / "scripts/R/redteam_replication_71sp.R"
    r_script_path.write_text(r_code)
    env = os.environ.copy()
    env["PATH"] = f"{ROOT / '.renv45/bin'}:{env.get('PATH', '')}"
    subprocess.run([str(ROOT / ".renv45/bin/Rscript"), str(r_script_path)], env=env, check=True)


if __name__ == "__main__":
    main()
