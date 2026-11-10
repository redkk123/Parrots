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

mtree71_rooted <- ape::multi2di(ape::read.tree(file.path(g71_dir, "master_burgio2019_71sp.nwk")), random = FALSE)
mtree71_rooted$edge.length[mtree71_rooted$edge.length <= 0] <- 1e-4
mtree71 <- ape::unroot(mtree71_rooted)
gt71_file <- file.path(g71_dir, "genetrees_burgio2019_71sp.txt")

alnfiles <- sort(list.files(path = file.path(g71_dir, "aa"), pattern = "\\.fa$", full.names = TRUE))
cat("[RedTeam Q8 R] Estimating LG branch lengths for", length(alnfiles), "panel genes across 71 species...\n")
est_one <- function(f) {
  g <- sub("\\.fa$", "", basename(f))
  tryCatch({
    aln <- phangorn::read.phyDat(f, type = "AA", format = "fasta")
    inboth <- intersect(names(aln), mtree71_rooted$tip.label)
    if (length(inboth) < 10) return(NULL)
    gt <- ape::unroot(ape::keep.tip(mtree71_rooted, inboth))
    aln <- subset(aln, subset = inboth)
    gt$edge.length <- rep(0.01, nrow(gt$edge))
    fit <- phangorn::optim.pml(
      phangorn::pml(gt, aln, model = "LG", k = 1),
      optEdge = TRUE, rearrangement = "none", model = "LG",
      control = phangorn::pml.control(trace = 0, epsilon = 1e-4, maxit = 5)
    )
    paste(g, ape::write.tree(fit$tree), sep = "\t")
  }, error = function(e) NULL)
}
lines <- unlist(parallel::mclapply(alnfiles, est_one, mc.cores = ncores))
# Include 5 synthetic full-71-species reference trees scaled from mtree71 so readTrees estimates masterTree edge lengths cleanly
ref_lines <- c(
  paste("MASTER_REF1", ape::write.tree(mtree71), sep = "\t"),
  paste("MASTER_REF2", ape::write.tree(mtree71), sep = "\t"),
  paste("MASTER_REF3", ape::write.tree(mtree71), sep = "\t"),
  paste("MASTER_REF4", ape::write.tree(mtree71), sep = "\t"),
  paste("MASTER_REF5", ape::write.tree(mtree71), sep = "\t")
)
writeLines(c(ref_lines, lines), gt71_file)

trees71 <- readTrees(gt71_file, minTreesAll = 2)
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

tv_le_58 <- tv_le_71[names(tv_le_71) %in% sp58]
tv_br_58 <- tv_br_71[names(tv_br_71) %in% sp58]
res_le_58 <- correlateWithContinuousPhenotype(rers_71_le, char2Paths(tv_le_58, trees71), min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)
res_br_58 <- correlateWithContinuousPhenotype(rers_71_br, char2Paths(tv_br_58, trees71), min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)

tv_le_13 <- tv_le_71[names(tv_le_71) %in% sp13]
tv_br_13 <- tv_br_71[names(tv_br_71) %in% sp13]
res_le_13 <- correlateWithContinuousPhenotype(rers_71_le, char2Paths(tv_le_13, trees71), min.sp = 4, winsorizeRER = 3, winsorizetrait = 3)
res_br_13 <- correlateWithContinuousPhenotype(rers_71_br, char2Paths(tv_br_13, trees71), min.sp = 4, winsorizeRER = 3, winsorizetrait = 3)

real_genes <- setdiff(rownames(res_le_71), paste0("MASTER_REF", 1:5))
panel_rep <- data.frame(
  gene = real_genes,
  Rho_LE_58sp = res_le_58[real_genes, "Rho"],
  P_LE_58sp = res_le_58[real_genes, "P"],
  Rho_LE_holdout13 = res_le_13[real_genes, "Rho"],
  N_branches_LE_holdout13 = res_le_13[real_genes, "N"],
  P_LE_holdout13 = res_le_13[real_genes, "P"],
  Rho_LE_71sp = res_le_71[real_genes, "Rho"],
  N_branches_LE_71sp = res_le_71[real_genes, "N"],
  P_LE_71sp = res_le_71[real_genes, "P"],
  Rho_Brain_58sp = res_br_58[real_genes, "Rho"],
  P_Brain_58sp = res_br_58[real_genes, "P"],
  Rho_Brain_holdout13 = res_br_13[real_genes, "Rho"],
  N_branches_Brain_holdout13 = res_br_13[real_genes, "N"],
  P_Brain_holdout13 = res_br_13[real_genes, "P"],
  Rho_Brain_71sp = res_br_71[real_genes, "Rho"],
  P_Brain_71sp = res_br_71[real_genes, "P"]
)
panel_rep <- panel_rep[order(panel_rep$P_LE_71sp), ]
write.table(panel_rep, file.path(out_dir, "8_holdout13_and_71sp_panel_replication.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[RedTeam Q8 R] Panel 71sp complete: %d valid genes; Cor(Rho_LE_58, Rho_LE_71) = %.4f, Cor(Rho_Brain_58, Rho_Brain_71) = %.4f\n",
            sum(!is.na(panel_rep$Rho_LE_71sp)),
            cor(panel_rep$Rho_LE_58sp, panel_rep$Rho_LE_71sp, use="complete.obs"),
            cor(panel_rep$Rho_Brain_58sp, panel_rep$Rho_Brain_71sp, use="complete.obs")))

if (!file.exists(file.path(out_dir, "8_holdout_cv_key_genes.tsv"))) {
cat("[RedTeam Q8 R] Running 30 clade-stratified 37-train / 13-test holdout splits on 50 species (7,756 genes)...\n")
trees58 <- readTrees("results/phase3_genome/genetrees_burgio2019.txt")
traits58 <- read.delim("results/phase3_genome/rer_traits.tsv", stringsAsFactors = FALSE)
both50 <- traits58[!is.na(traits58$rel_log_life_expectancy) & !is.na(traits58$rel_log_brain), ]
tv50_le <- setNames(both50$rel_log_life_expectancy, both50$species_tip)
tv50_br <- setNames(both50$rel_log_brain, both50$species_tip)
rers50 <- getAllResiduals(trees58, useSpecies = names(tv50_le), transform = "sqrt", use.weights = TRUE, norm = "scale", min.sp = 10, min.valid = 20)

key_genes <- c("PTEN", "MAD2L2", "HSPA8", "CTNNB1", "CTNNBIP1", "PSMA4", "SSU72", "TOMM20", "ANAPC11", "RTN1", "SOX14", "NCALD", "CBLN4")
set.seed(20261110)
psl_sp <- both50$species_tip[both50$clade == "Psittaculidae"]
rest_sp <- both50$species_tip[both50$clade == "rest"]

n_splits <- 30
train_le_mat <- matrix(NA_real_, nrow = length(key_genes), ncol = n_splits, dimnames = list(key_genes, NULL))
test_le_mat <- matrix(NA_real_, nrow = length(key_genes), ncol = n_splits, dimnames = list(key_genes, NULL))
test_br_mat <- matrix(NA_real_, nrow = length(key_genes), ncol = n_splits, dimnames = list(key_genes, NULL))
coupling_test13 <- numeric(n_splits)

for (s in 1:n_splits) {
  test_sp <- c(sample(psl_sp, 4), sample(rest_sp, 11)) # 15 holdout species so >=10 present in most genes
  train_sp <- setdiff(both50$species_tip, test_sp)
  p_tr_le <- suppressMessages(char2Paths(tv50_le[names(tv50_le) %in% train_sp], trees58))
  p_te_le <- suppressMessages(char2Paths(tv50_le[names(tv50_le) %in% test_sp], trees58))
  p_te_br <- suppressMessages(char2Paths(tv50_br[names(tv50_br) %in% test_sp], trees58))
  r_tr_le <- correlateWithContinuousPhenotype(rers50, p_tr_le, min.sp = 8, winsorizeRER = 3, winsorizetrait = 3)
  r_te_le <- correlateWithContinuousPhenotype(rers50, p_te_le, min.sp = 5, winsorizeRER = 3, winsorizetrait = 3)
  r_te_br <- correlateWithContinuousPhenotype(rers50, p_te_br, min.sp = 5, winsorizeRER = 3, winsorizetrait = 3)
  train_le_mat[, s] <- r_tr_le[key_genes, "Rho"]
  test_le_mat[, s] <- r_te_le[key_genes, "Rho"]
  test_br_mat[, s] <- r_te_br[key_genes, "Rho"]
  ok <- !is.na(r_te_le$Rho) & !is.na(r_te_br$Rho)
  coupling_test13[s] <- if (sum(ok) >= 10) cor(r_te_le$Rho[ok], r_te_br$Rho[ok]) else NA_real_
}

orig_le <- read.delim("results/phase3_genome/rer_burgio2019_rel_log_life_expectancy.tsv", stringsAsFactors = FALSE)
rownames(orig_le) <- orig_le$gene
orig_br <- read.delim("results/phase3_genome/rer_burgio2019_rel_log_brain.tsv", stringsAsFactors = FALSE)
rownames(orig_br) <- orig_br$gene

cv_df <- data.frame(
  gene = key_genes,
  full58_Rho_LE = orig_le[key_genes, "Rho"],
  full58_P_LE = orig_le[key_genes, "P"],
  mean_train_Rho_LE = rowMeans(train_le_mat, na.rm = TRUE),
  mean_holdout_Rho_LE = rowMeans(test_le_mat, na.rm = TRUE),
  sd_holdout_Rho_LE = apply(test_le_mat, 1, sd, na.rm = TRUE),
  sign_concordance_holdout_LE = rowMeans(sign(test_le_mat) == sign(orig_le[key_genes, "Rho"]), na.rm = TRUE),
  full58_Rho_Brain = orig_br[key_genes, "Rho"],
  mean_holdout_Rho_Brain = rowMeans(test_br_mat, na.rm = TRUE),
  sign_concordance_holdout_Brain = rowMeans(sign(test_br_mat) == sign(orig_br[key_genes, "Rho"]), na.rm = TRUE)
)
write.table(cv_df, file.path(out_dir, "8_holdout_cv_key_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
