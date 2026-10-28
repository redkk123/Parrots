# RERconverge for H6 (genome-wide: longevity, brain, brain beyond longevity) and H8 (panel genes within Psittaculidae
# and within the rest), as fixed in results/phase0/analysis_plan.md (2026-10-28).
# NOT YET RUN: written against the RERconverge vignette API (v0.3) before the package could be installed; check
# function signatures against the installed version on the first run.
# Inputs: data/genome_alignments/aa/*.fa (scripts/34_align_genome.py), results/phase3_genome/{rer_traits.tsv,
# master_<tree>.nwk} (scripts/35_prepare_rer_inputs.py), optional GO gene sets data/external/go/go_bp.gmt (human
# symbols). Usage: Rscript scripts/R/rerconverge_h6_h8.R [burgio2019|smith2024]   (default burgio2019)
# Outputs: results/phase3_genome/rer_<tree>_<trait>.tsv, rer_<tree>_h8_<clade>.tsv, rer_<tree>_enrichment_<trait>.tsv
suppressPackageStartupMessages(library(RERconverge))
args <- commandArgs(trailingOnly = TRUE)
tree_name <- if (length(args)) args[1] else "burgio2019"
out <- "results/phase3_genome"
nperm <- as.integer(Sys.getenv("RER_PERMS", "1000"))
set.seed(20261028)

master <- file.path(out, paste0("master_", tree_name, ".nwk"))
genetrees <- file.path(out, paste0("genetrees_", tree_name, ".txt"))
if (!file.exists(genetrees)) {
  # branch lengths per gene on the fixed topology (phangorn, LG)
  estimatePhangornTreeAll(alndir = "data/genome_alignments/aa", pattern = "\\.fa$", treefile = master,
                          output.file = genetrees, submodel = "LG", type = "AA", format = "fasta")
}
trees <- readTrees(genetrees, minSpecs = 10)
traits <- read.delim(file.path(out, "rer_traits.tsv"))
mtree <- ape::read.tree(master)

run_trait <- function(col, species = NULL, label = col) {
  tv <- setNames(traits[[col]], traits$species_tip)
  tv <- tv[!is.na(tv)]
  if (!is.null(species)) tv <- tv[names(tv) %in% species]
  rers <- getAllResiduals(trees, useSpecies = names(tv), transform = "sqrt", weighted = TRUE, scale = TRUE,
                          plot = FALSE)
  cp <- char2Paths(tv, trees)
  res <- correlateWithContinuousPhenotype(rers, cp, min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)
  perms <- getPermsContinuous(nperm, tv, rers, annotlist = NULL, trees = trees,
                              mastertree = ape::keep.tip(mtree, names(tv)), calculateenrich = FALSE)
  res$permP <- permpvalcor(res, perms)
  res$permQ <- p.adjust(res$permP, "BH")
  res$gene <- rownames(res)
  write.table(res, file.path(out, sprintf("rer_%s_%s.tsv", tree_name, label)), sep = "\t", quote = FALSE,
              row.names = FALSE)
  gmt <- "data/external/go/go_bp.gmt"
  if (file.exists(gmt) && is.null(species)) {
    annots <- list(GO_BP = read.gmt(gmt))
    enr <- fastwilcoxGMTall(getStat(res), annots, outputGeneVals = FALSE)
    write.table(enr$GO_BP, file.path(out, sprintf("rer_%s_enrichment_%s.tsv", tree_name, label)), sep = "\t",
                quote = FALSE)
  }
  invisible(res)
}

# H6: genome-wide, each phenotype separately; brain beyond longevity is the mandatory analysis for shared genes
for (col in c("rel_log_life_expectancy", "rel_log_brain", "rel_brain_given_le")) run_trait(col)

# H8 item 2: panel genes within Psittaculidae and within the rest (longevity)
panel <- read.delim("data/pilot_genes/canonical_genes.tsv")$gene_symbol
for (cl in c("Psittaculidae", "rest")) {
  sp <- traits$species_tip[traits$clade == cl]
  res <- run_trait("rel_log_life_expectancy", species = sp, label = paste0("h8_", cl))
  keep <- res[res$gene %in% panel, ]
  keep$panelQ <- p.adjust(keep$permP, "BH")
  write.table(keep, file.path(out, sprintf("rer_%s_h8_%s_panel.tsv", tree_name, cl)), sep = "\t", quote = FALSE,
              row.names = FALSE)
}
cat("RERconverge done:", tree_name, "\n")
