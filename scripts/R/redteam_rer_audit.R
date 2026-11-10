#!/usr/bin/env Rscript
# Red-Team Audit of Phase 3 RERconverge Results (Tests 1, 2, 4, 5)
suppressPackageStartupMessages({
  library(Rcpp)
  library(RERconverge)
  library(TreeTools)
  library(impute)
  library(geiger)
  library(parallel)
  library(ape)
  library(phytools)
})

# Exact C++ implementation from scripts/R/rerconverge_h6_h8.R
cppFunction("
NumericMatrix fastCorPerms(NumericMatrix RER, NumericMatrix Pmat, int min_sp, int winN) {
  int ng = RER.nrow(), nb = RER.ncol(), np = Pmat.ncol();
  NumericMatrix out(ng, np);
  std::fill(out.begin(), out.end(), NA_REAL);
  for (int i = 0; i < ng; ++i) {
    std::vector<int> idx;
    for (int j = 0; j < nb; ++j) {
      if (!NumericMatrix::is_na(RER(i, j)) && !NumericMatrix::is_na(Pmat(j, 0))) idx.push_back(j);
    }
    int n = idx.size();
    if (n < min_sp) continue;
    std::vector<double> x(n), xs(n);
    for (int k = 0; k < n; ++k) x[k] = xs[k] = RER(i, idx[k]);
    std::sort(xs.begin(), xs.end(), std::greater<double>());
    if (n > 2 * winN) {
      double hi = xs[winN], lo = xs[n - 1 - winN];
      for (int k = 0; k < n; ++k) { if (x[k] > hi) x[k] = hi; if (x[k] < lo) x[k] = lo; }
    }
    double mx = 0; for (int k = 0; k < n; ++k) mx += x[k]; mx /= n;
    double sx = 0; for (int k = 0; k < n; ++k) { x[k] -= mx; sx += x[k]*x[k]; }
    if (sx <= 0) continue;
    sx = std::sqrt(sx);
    std::vector<double> y(n), ys(n);
    for (int p = 0; p < np; ++p) {
      for (int k = 0; k < n; ++k) y[k] = ys[k] = Pmat(idx[k], p);
      std::sort(ys.begin(), ys.end(), std::greater<double>());
      if (n > 2 * winN) {
        double hi = ys[winN], lo = ys[n - 1 - winN];
        for (int k = 0; k < n; ++k) { if (y[k] > hi) y[k] = hi; if (y[k] < lo) y[k] = lo; }
      }
      double my = 0; for (int k = 0; k < n; ++k) my += y[k]; my /= n;
      double sy = 0, dot = 0;
      for (int k = 0; k < n; ++k) { double dy = y[k] - my; sy += dy*dy; dot += x[k]*dy; }
      if (sy > 0) out(i, p) = dot / (sx * std::sqrt(sy));
    }
  }
  return out;
}
")

out_dir <- "results/redteam_audit"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ncores <- as.integer(Sys.getenv("RER_CORES", "8"))

cat("[RedTeam] Loading trees from genetrees_burgio2019.txt...\n")
trees <- readTrees("results/phase3_genome/genetrees_burgio2019.txt")

traits <- read.delim("results/phase3_genome/rer_traits.tsv", stringsAsFactors = FALSE)
ph_full <- read.delim("results/phase0/phenotypes.tsv", stringsAsFactors = FALSE)
ph_full$species_tip <- gsub(" ", "_", ph_full$species)
traits <- merge(traits, ph_full[, c("species_tip", "log_body_mass", "clutch_size", "incubation_days", "fledging_age_days")],
                by = "species_tip", all.x = TRUE)

pdf(NULL)
tv_le <- setNames(traits$rel_log_life_expectancy, traits$species_tip); tv_le <- tv_le[!is.na(tv_le)]
tv_br <- setNames(traits$rel_log_brain, traits$species_tip); tv_br <- tv_br[!is.na(tv_br)]

cat("[RedTeam] Computing RER matrices (58 sp for LE, 50 sp for Brain)...\n")
rers_58 <- getAllResiduals(trees, useSpecies = names(tv_le), transform = "sqrt",
                           use.weights = TRUE, norm = "scale", min.sp = 10, min.valid = 20)
rers_50 <- getAllResiduals(trees, useSpecies = names(tv_br), transform = "sqrt",
                           use.weights = TRUE, norm = "scale", min.sp = 10, min.valid = 20)

run_with_seed <- function(tv, rers_mat, seed_base, nperm = 1000) {
  min_sp <- 10
  cp <- char2Paths(tv, trees)
  res <- correlateWithContinuousPhenotype(rers_mat, cp, min.sp = min_sp, winsorizeRER = 3, winsorizetrait = 3)
  mt_root <- ape::multi2di(phangorn::midpoint(trees$masterTree), random = FALSE)
  mt_sub <- ape::multi2di(ape::keep.tip(mt_root, names(tv)), random = FALSE)
  mt_sub$edge.length[mt_sub$edge.length <= 0] <- 1e-6
  Pmat <- do.call(cbind, parallel::mclapply(1:nperm, function(k) {
    set.seed(seed_base + k)
    suppressMessages(char2Paths(simpermvec(tv, mt_sub), trees))
  }, mc.cores = ncores))
  corRho <- fastCorPerms(rers_mat, Pmat, min_sp, 3)
  rownames(corRho) <- rownames(rers_mat)
  perms <- list(corRho = corRho)
  pv <- permpvalcor(res, perms)
  res$permP <- pv$permpval
  res$permStat <- pv$permstats
  res$permQ <- p.adjust(res$permP, "BH")
  res$gene <- rownames(res)
  list(res = res, cp = cp, Pmat = Pmat, corRho = corRho)
}

# ==============================================================================
# TEST 1: Reproducibility & Random Seed Sensitivity
# ==============================================================================
cat("[RedTeam] Running Test 1 (Reproducibility & Seed Sensitivity)...\n")
orig_le <- read.delim("results/phase3_genome/rer_burgio2019_rel_log_life_expectancy.tsv", stringsAsFactors = FALSE)
rownames(orig_le) <- orig_le$gene
orig_br <- read.delim("results/phase3_genome/rer_burgio2019_rel_log_brain.tsv", stringsAsFactors = FALSE)
rownames(orig_br) <- orig_br$gene
orig_br_given_le <- read.delim("results/phase3_genome/rer_burgio2019_rel_brain_given_le.tsv", stringsAsFactors = FALSE)
rownames(orig_br_given_le) <- orig_br_given_le$gene

seed2_le_out <- run_with_seed(tv_le, rers_58, 999999, 200)
seed2_br_out <- run_with_seed(tv_br, rers_50, 999999, 200)
seed2_le <- seed2_le_out$res[orig_le$gene, ]
seed2_br <- seed2_br_out$res[orig_br$gene, ]

# Also compute Rho_LE on the exact 50-species RER matrix (rers_50) so we compare 58-sp vs 50-sp RER matrix
cp_le_50 <- char2Paths(tv_le[names(tv_br)], trees)
res_le_50 <- correlateWithContinuousPhenotype(rers_50, cp_le_50, min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)

if (!file.exists(file.path(out_dir, "1_reproducibility.tsv"))) {
  repro_df <- data.frame(
    metric = c(
      "r(Rho_LE, Rho_Brain) [seed=20261028, saved files]",
      "r(Rho_LE, Rho_Brain) [seed=999999, rerun]",
      "r(Rho_LE_50sp, Rho_Brain_50sp) [both on identical 50-species RER matrix]",
      "max_abs_diff_Rho_LE (saved vs seed 999999 rerun)",
      "max_abs_diff_Rho_Brain (saved vs seed 999999 rerun)",
      "cor(permP_LE_seed20261028, permP_LE_seed999999)",
      "cor(permP_Brain_seed20261028, permP_Brain_seed999999)",
      "median_abs_diff_permP_LE",
      "median_abs_diff_permP_Brain",
      "n_shared_sig_genes_seed20261028 (P<0.05 & permP<0.05 both traits, same sign)",
      "n_shared_sig_genes_seed999999 (P<0.05 & permP<0.05 both traits, same sign)",
      "n_shared_sig_genes_overlap_both_seeds"
    ),
    value = c(
      cor(orig_le$Rho, orig_br$Rho, use = "complete.obs"),
      cor(seed2_le$Rho, seed2_br$Rho, use = "complete.obs"),
      cor(res_le_50$Rho, orig_br$Rho, use = "complete.obs"),
      max(abs(orig_le$Rho - seed2_le$Rho), na.rm = TRUE),
      max(abs(orig_br$Rho - seed2_br$Rho), na.rm = TRUE),
      cor(orig_le$permP, seed2_le$permP, use = "complete.obs"),
      cor(orig_br$permP, seed2_br$permP, use = "complete.obs"),
      median(abs(orig_le$permP - seed2_le$permP), na.rm = TRUE),
      median(abs(orig_br$permP - seed2_br$permP), na.rm = TRUE),
      sum(orig_le$P < 0.05 & orig_le$permP < 0.05 & orig_br$P < 0.05 & orig_br$permP < 0.05 & sign(orig_le$Rho) == sign(orig_br$Rho), na.rm = TRUE),
      sum(seed2_le$P < 0.05 & seed2_le$permP < 0.05 & seed2_br$P < 0.05 & seed2_br$permP < 0.05 & sign(seed2_le$Rho) == sign(seed2_br$Rho), na.rm = TRUE),
      sum((orig_le$P < 0.05 & orig_le$permP < 0.05 & orig_br$P < 0.05 & orig_br$permP < 0.05 & sign(orig_le$Rho) == sign(orig_br$Rho)) &
          (seed2_le$P < 0.05 & seed2_le$permP < 0.05 & seed2_br$P < 0.05 & seed2_br$permP < 0.05 & sign(seed2_le$Rho) == sign(seed2_br$Rho)), na.rm = TRUE)
    )
  )
  write.table(repro_df, file.path(out_dir, "1_reproducibility.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

# ==============================================================================
# TEST 2: Noise Floor of Coupling, Unrelated Phenotypes & Joint Permulations
# ==============================================================================
cat("[RedTeam] Running Test 2 (Noise floor, unrelated traits, independent & coupled permulations)...\n")
paths_le <- seed2_le_out$cp
paths_br <- seed2_br_out$cp

# Unrelated traits (computed against rers_58 on shared species)
tv_mass <- setNames(traits$log_body_mass, traits$species_tip); tv_mass <- tv_mass[!is.na(tv_mass)]
paths_mass <- char2Paths(tv_mass, trees)
res_mass <- correlateWithContinuousPhenotype(rers_58, paths_mass, min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)

sub_cs <- traits[!is.na(traits$clutch_size) & !is.na(traits$log_body_mass), ]
tv_clutch <- setNames(sub_cs$clutch_size, sub_cs$species_tip)
paths_clutch <- char2Paths(tv_clutch, trees)
res_clutch <- correlateWithContinuousPhenotype(rers_58, paths_clutch, min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)

tv_rel_clutch <- setNames(resid(lm(clutch_size ~ log_body_mass, data = sub_cs)), sub_cs$species_tip)
paths_rel_clutch <- char2Paths(tv_rel_clutch, trees)
res_rel_clutch <- correlateWithContinuousPhenotype(rers_58, paths_rel_clutch, min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)

sub_inc <- traits[!is.na(traits$incubation_days) & !is.na(traits$log_body_mass), ]
tv_rel_inc <- setNames(resid(lm(incubation_days ~ log_body_mass, data = sub_inc)), sub_inc$species_tip)
paths_rel_inc <- char2Paths(tv_rel_inc, trees)
res_rel_inc <- correlateWithContinuousPhenotype(rers_58, paths_rel_inc, min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)

# Null 2A: 200 independent permulated phenotype pairs (seed 999999 for LE vs seed 555555 for Brain)
seed3_br_out <- run_with_seed(tv_br, rers_50, 555555, 200)
corMat_le <- seed2_le_out$corRho
corMat_br_indep <- seed3_br_out$corRho
r_null_indep <- sapply(1:200, function(b) cor(corMat_le[, b], corMat_br_indep[, b], use = "complete.obs"))
path_cor_indep <- sapply(1:200, function(b) cor(seed2_le_out$Pmat[, b], seed3_br_out$Pmat[, b], use = "complete.obs"))

# Null 2B: 200 phenotypically coupled Brownian motion pairs (simulated jointly with 2x2 ratematrix on the 50 species)
both50 <- traits[!is.na(traits$rel_log_life_expectancy) & !is.na(traits$rel_log_brain), ]
mt_root <- ape::multi2di(phangorn::midpoint(trees$masterTree), random = FALSE)
mt50 <- ape::multi2di(ape::keep.tip(mt_root, both50$species_tip), random = FALSE)
mt50$edge.length[mt50$edge.length <= 0] <- 1e-6
trait_mat50 <- cbind(LE = both50$rel_log_life_expectancy, Brain = both50$rel_log_brain)
rownames(trait_mat50) <- both50$species_tip
v_2x2 <- geiger::ratematrix(mt50, trait_mat50)
set.seed(20261110)
sims_2x2 <- geiger::sim.char(mt50, v_2x2, nsim = 200, model = "BM", root = 0)

P1_coupled <- do.call(cbind, parallel::mclapply(1:200, function(b) {
  suppressMessages(char2Paths(sims_2x2[, 1, b], trees))
}, mc.cores = ncores))
P2_coupled <- do.call(cbind, parallel::mclapply(1:200, function(b) {
  suppressMessages(char2Paths(sims_2x2[, 2, b], trees))
}, mc.cores = ncores))

cor1_coupled <- fastCorPerms(rers_50, P1_coupled, 10, 3)
cor2_coupled <- fastCorPerms(rers_50, P2_coupled, 10, 3)
r_null_coupled <- sapply(1:200, function(b) cor(cor1_coupled[, b], cor2_coupled[, b], use = "complete.obs"))
path_cor_coupled <- sapply(1:200, function(b) cor(P1_coupled[, b], P2_coupled[, b], use = "complete.obs"))

# Null 2C: Observed LE fixed vs 200 permulations of Brain
r_obsLE_vs_permBr <- sapply(1:200, function(b) cor(orig_le$Rho, corMat_br_indep[, b], use = "complete.obs"))

coupling_df <- data.frame(
  comparison = c(
    "Observed: rel_log_life_expectancy vs rel_log_brain (saved files)",
    "Observed (50sp matched RER): rel_log_life_expectancy vs rel_log_brain",
    "Control 1: rel_log_life_expectancy vs log_body_mass (n=58)",
    "Control 2: rel_log_life_expectancy vs clutch_size raw (n=51)",
    "Control 3: rel_log_life_expectancy vs rel_clutch_size (residual on mass, n=51)",
    "Control 4: rel_log_life_expectancy vs rel_incubation_days (residual on mass, n=50)",
    "Null 2A: Independent permulated traits (200 pairs)",
    "Null 2B: Jointly simulated coupled BM traits (200 pairs, empirical 2x2 ratematrix)",
    "Null 2C: Observed LE fixed vs 200 permulated Brain traits"
  ),
  tip_cor_pearson = c(
    cor(both50$rel_log_life_expectancy, both50$rel_log_brain),
    cor(both50$rel_log_life_expectancy, both50$rel_log_brain),
    cor(traits$rel_log_life_expectancy, traits$log_body_mass, use = "complete.obs"),
    cor(sub_cs$rel_log_life_expectancy, sub_cs$clutch_size, use = "complete.obs"),
    cor(sub_cs$rel_log_life_expectancy, tv_rel_clutch),
    cor(sub_inc$rel_log_life_expectancy, tv_rel_inc),
    0.0,
    mean(sapply(1:200, function(b) cor(sims_2x2[,1,b], sims_2x2[,2,b]))),
    0.0
  ),
  branch_path_cor = c(
    cor(paths_le, paths_br, use = "complete.obs"),
    cor(cp_le_50, paths_br, use = "complete.obs"),
    cor(paths_le, paths_mass, use = "complete.obs"),
    cor(paths_le, paths_clutch, use = "complete.obs"),
    cor(paths_le, paths_rel_clutch, use = "complete.obs"),
    cor(paths_le, paths_rel_inc, use = "complete.obs"),
    mean(path_cor_indep),
    mean(path_cor_coupled),
    mean(sapply(1:200, function(b) cor(paths_le, seed3_br_out$Pmat[, b], use = "complete.obs")))
  ),
  genome_Rho_cor = c(
    cor(orig_le$Rho, orig_br$Rho, use = "complete.obs"),
    cor(res_le_50$Rho, orig_br$Rho, use = "complete.obs"),
    cor(orig_le$Rho, res_mass$Rho, use = "complete.obs"),
    cor(orig_le$Rho, res_clutch$Rho, use = "complete.obs"),
    cor(orig_le$Rho, res_rel_clutch$Rho, use = "complete.obs"),
    cor(orig_le$Rho, res_rel_inc$Rho, use = "complete.obs"),
    mean(r_null_indep),
    mean(r_null_coupled),
    mean(r_obsLE_vs_permBr)
  ),
  notes = c(
    sprintf("r_Rho = %.4f; branch_path_cor = %.4f", cor(orig_le$Rho, orig_br$Rho, use="complete.obs"), cor(paths_le, paths_br, use="complete.obs")),
    sprintf("r_Rho = %.4f; branch_path_cor = %.4f", cor(res_le_50$Rho, orig_br$Rho, use="complete.obs"), cor(cp_le_50, paths_br, use="complete.obs")),
    sprintf("r_Rho = %.4f; branch_path_cor = %.4f", cor(orig_le$Rho, res_mass$Rho, use="complete.obs"), cor(paths_le, paths_mass, use="complete.obs")),
    sprintf("r_Rho = %.4f; branch_path_cor = %.4f", cor(orig_le$Rho, res_clutch$Rho, use="complete.obs"), cor(paths_le, paths_clutch, use="complete.obs")),
    sprintf("r_Rho = %.4f; branch_path_cor = %.4f", cor(orig_le$Rho, res_rel_clutch$Rho, use="complete.obs"), cor(paths_le, paths_rel_clutch, use="complete.obs")),
    sprintf("r_Rho = %.4f; branch_path_cor = %.4f", cor(orig_le$Rho, res_rel_inc$Rho, use="complete.obs"), cor(paths_le, paths_rel_inc, use="complete.obs")),
    sprintf("mean=%.4f, SD=%.4f, mean|r|=%.4f, 95%%=[%.4f, %.4f], max=%.4f, p(r>=0.666)=%.4g; cor(r_Rho, path_cor)=%.4f",
            mean(r_null_indep), sd(r_null_indep), mean(abs(r_null_indep)), quantile(r_null_indep, 0.025), quantile(r_null_indep, 0.975), max(r_null_indep), (sum(r_null_indep >= 0.6658)+1)/201, cor(r_null_indep, path_cor_indep)),
    sprintf("mean=%.4f, SD=%.4f, 95%%=[%.4f, %.4f], max=%.4f, p(r>=0.666)=%.4g; cor(r_Rho, path_cor)=%.4f",
            mean(r_null_coupled), sd(r_null_coupled), quantile(r_null_coupled, 0.025), quantile(r_null_coupled, 0.975), max(r_null_coupled), (sum(r_null_coupled >= 0.6658)+1)/201, cor(r_null_coupled, path_cor_coupled)),
    sprintf("mean=%.4f, SD=%.4f, mean|r|=%.4f, 95%%=[%.4f, %.4f], max=%.4f, p(r>=0.666)=%.4g",
            mean(r_obsLE_vs_permBr), sd(r_obsLE_vs_permBr), mean(abs(r_obsLE_vs_permBr)), quantile(r_obsLE_vs_permBr, 0.025), quantile(r_obsLE_vs_permBr, 0.975), max(r_obsLE_vs_permBr), (sum(r_obsLE_vs_permBr >= 0.6658)+1)/201)
  )
)
write.table(coupling_df, file.path(out_dir, "2_coupling_null_and_controls.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ==============================================================================
# TEST 4: Leave-One-Clade-Out Sensitivity
# ==============================================================================
cat("[RedTeam] Running Test 4 (Leave-One-Clade-Out)...\n")
clade_defs <- list(
  "All_50_species_with_brain" = character(0),
  "Drop_Strigopidae (Strigops + Nestor, -2 sp)" = c("Strigops_habroptila", "Nestor_notabilis"),
  "Drop_Cacatuidae (Cacatua, Calyptorhynchus, Eolophus, Probosciger, -6 sp)" = both50$species_tip[grep("^(Cacatua|Calyptorhynchus|Eolophus|Probosciger)_", both50$species_tip)],
  "Drop_Loriini (Chalcopsitta, Eos, Trichoglossus, -5 sp)" = both50$species_tip[grep("^(Chalcopsitta|Eos|Trichoglossus|Lorius|Vini|Parvipsitta)_", both50$species_tip)],
  "Drop_Melopsittacus_undulatus (-1 sp)" = c("Melopsittacus_undulatus"),
  "Drop_All_Psittaculidae (-13 sp, leave 37 non-Psittaculidae)" = both50$species_tip[both50$clade == "Psittaculidae"],
  "Drop_Arini_Macaws_Conures (Ara, Anodorhynchus, Pyrrhura, etc., -12 sp)" = both50$species_tip[grep("^(Ara|Anodorhynchus|Primolius|Orthopsittaca|Diopsittaca|Guaruba|Aratinga|Eupsittula|Psittacara|Pyrrhura|Enicognathus|Cyanoliseus|Nandayus|Brotogeris|Pionites|Deroptyus)_", both50$species_tip)],
  "Drop_Androglossini (Amazona, Pionus, Pyrilia, Triclaria, -9 sp)" = both50$species_tip[grep("^(Amazona|Pionus|Pyrilia|Graydidascalus|Triclaria|Alipiopsitta)_", both50$species_tip)]
)

loco_rows <- list()
for (cname in names(clade_defs)) {
  drop_sp <- clade_defs[[cname]]
  keep_df <- both50[!both50$species_tip %in% drop_sp, ]
  t_le <- setNames(keep_df$rel_log_life_expectancy, keep_df$species_tip)
  t_br <- setNames(keep_df$rel_log_brain, keep_df$species_tip)
  p_le <- suppressMessages(char2Paths(t_le, trees))
  p_br <- suppressMessages(char2Paths(t_br, trees))
  r_le <- correlateWithContinuousPhenotype(rers_50, p_le, min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)
  r_br <- correlateWithContinuousPhenotype(rers_50, p_br, min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)
  loco_rows[[length(loco_rows) + 1]] <- data.frame(
    subset = cname,
    n_species = nrow(keep_df),
    tip_cor_pearson = cor(t_le, t_br),
    tip_cor_spearman = cor(t_le, t_br, method = "spearman"),
    branch_path_cor = cor(p_le, p_br, use = "complete.obs"),
    genome_Rho_cor = cor(r_le$Rho, r_br$Rho, use = "complete.obs"),
    PTEN_Rho_LE = r_le["PTEN", "Rho"],
    PTEN_Rho_Brain = r_br["PTEN", "Rho"],
    MAD2L2_Rho_LE = r_le["MAD2L2", "Rho"],
    MAD2L2_Rho_Brain = r_br["MAD2L2", "Rho"]
  )
}
loco_df <- do.call(rbind, loco_rows)
write.table(loco_df, file.path(out_dir, "4_leave_one_clade_out.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ==============================================================================
# TEST 5: True Residual Brain & True Partial Correlation
# ==============================================================================
cat("[RedTeam] Running Test 5 (True Residual Brain & Partial Correlation)...\n")
# 5A. Tip-level OLS residual on the 50 species: rel_log_brain ~ rel_log_life_expectancy
tv_br_ols <- setNames(resid(lm(rel_log_brain ~ rel_log_life_expectancy, data = both50)), both50$species_tip)
res_br_ols_out <- run_with_seed(tv_br_ols, rers_50, 20261028, 200)
res_br_ols <- res_br_ols_out$res[orig_le$gene, ]

# 5B. Branch-path-level orthogonalization: P_Brain orthogonalized against P_LE on the 50-species branches
paths_le50 <- cp_le_50
paths_br50 <- paths_br
ok_b <- !is.na(paths_le50) & !is.na(paths_br50)
paths_br_orth_path <- rep(NA_real_, length(paths_br50))
paths_br_orth_path[ok_b] <- resid(lm(paths_br50[ok_b] ~ paths_le50[ok_b]))
paths_le_orth_path <- rep(NA_real_, length(paths_le50))
paths_le_orth_path[ok_b] <- resid(lm(paths_le50[ok_b] ~ paths_br50[ok_b]))

res_br_orth_path <- correlateWithContinuousPhenotype(rers_50, paths_br_orth_path, min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)
res_le_orth_path <- correlateWithContinuousPhenotype(rers_50, paths_le_orth_path, min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)

# 5C. Exact per-gene partial correlation: r(RER_g, P_Brain | P_LE) and r(RER_g, P_LE | P_Brain)
r_le50_vec <- res_le_50[orig_le$gene, "Rho"]
r_br50_vec <- orig_br[orig_le$gene, "Rho"]
c_path <- cor(paths_le50, paths_br50, use = "complete.obs")
partial_br_given_le <- (r_br50_vec - r_le50_vec * c_path) / sqrt(pmax(1e-12, (1 - r_le50_vec^2) * (1 - c_path^2)))
partial_le_given_br <- (r_le50_vec - r_br50_vec * c_path) / sqrt(pmax(1e-12, (1 - r_br50_vec^2) * (1 - c_path^2)))

key_genes <- c("PTEN", "CTNNB1", "CTNNBIP1", "PSMA4", "SSU72", "TOMM20", "ANAPC11", "MAD2L2", "HSPA8", "RTN1", "SOX14", "NCALD", "CBLN4")
key_df <- data.frame(
  gene = key_genes,
  Rho_LE_58sp = orig_le[key_genes, "Rho"],
  P_LE_58sp = orig_le[key_genes, "P"],
  Rho_Brain_raw = orig_br[key_genes, "Rho"],
  P_Brain_raw = orig_br[key_genes, "P"],
  Rho_Brain_PGLS_beta004 = orig_br_given_le[key_genes, "Rho"],
  P_Brain_PGLS_beta004 = orig_br_given_le[key_genes, "P"],
  Rho_Brain_OLS_tip_orth = res_br_ols[key_genes, "Rho"],
  P_Brain_OLS_tip_orth = res_br_ols[key_genes, "P"],
  permP_Brain_OLS_tip_orth = res_br_ols[key_genes, "permP"],
  Rho_Brain_branch_orth = res_br_orth_path[key_genes, "Rho"],
  P_Brain_branch_orth = res_br_orth_path[key_genes, "P"],
  Rho_LE_branch_orth_Brain = res_le_orth_path[key_genes, "Rho"],
  P_LE_branch_orth_Brain = res_le_orth_path[key_genes, "P"]
)
write.table(key_df, file.path(out_dir, "5_residual_brain_key_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summary_5 <- data.frame(
  method = c(
    "1. Unadjusted rel_log_brain (50 sp)",
    "2. PGLS residual rel_brain_given_le (scripts/35_prepare_rer_inputs.py, beta_GLS=0.040, lambda=0.90)",
    "3. Tip-OLS orthogonalized rel_brain_given_le_OLS (beta_OLS=0.258 on 50 tips, tip_cor=0.000)",
    "4. Branch-path orthogonalized paths_brain_orth_le (branch_path_cor=0.000)",
    "5. Per-gene analytical partial correlation r(RER_g, P_Brain | P_LE)"
  ),
  tip_cor_with_LE = c(
    cor(both50$rel_log_life_expectancy, both50$rel_log_brain),
    cor(both50$rel_log_life_expectancy, both50$rel_brain_given_le),
    cor(both50$rel_log_life_expectancy, tv_br_ols),
    NA_real_,
    NA_real_
  ),
  branch_path_cor_with_LE = c(
    cor(paths_le50, paths_br50, use = "complete.obs"),
    cor(paths_le50, char2Paths(setNames(both50$rel_brain_given_le, both50$species_tip), trees), use = "complete.obs"),
    cor(paths_le50, res_br_ols_out$cp, use = "complete.obs"),
    cor(paths_le50, paths_br_orth_path, use = "complete.obs"),
    0.0
  ),
  genome_Rho_cor_with_LE_58sp = c(
    cor(orig_le$Rho, orig_br$Rho, use = "complete.obs"),
    cor(orig_le$Rho, orig_br_given_le$Rho, use = "complete.obs"),
    cor(orig_le$Rho, res_br_ols$Rho, use = "complete.obs"),
    cor(orig_le$Rho, res_br_orth_path[orig_le$gene, "Rho"], use = "complete.obs"),
    cor(orig_le$Rho, partial_br_given_le, use = "complete.obs")
  ),
  genome_Rho_cor_with_LE_50sp = c(
    cor(r_le50_vec, orig_br$Rho, use = "complete.obs"),
    cor(r_le50_vec, orig_br_given_le$Rho, use = "complete.obs"),
    cor(r_le50_vec, res_br_ols$Rho, use = "complete.obs"),
    cor(r_le50_vec, res_br_orth_path[orig_le$gene, "Rho"], use = "complete.obs"),
    cor(r_le50_vec, partial_br_given_le, use = "complete.obs")
  ),
  genome_Rho_cor_with_Brain_raw = c(
    1.0,
    cor(orig_br$Rho, orig_br_given_le$Rho, use = "complete.obs"),
    cor(orig_br$Rho, res_br_ols$Rho, use = "complete.obs"),
    cor(orig_br$Rho, res_br_orth_path[orig_le$gene, "Rho"], use = "complete.obs"),
    cor(orig_br$Rho, partial_br_given_le, use = "complete.obs")
  )
)
write.table(summary_5, file.path(out_dir, "5_residual_brain_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("[RedTeam] R audit script completed successfully.\n")
