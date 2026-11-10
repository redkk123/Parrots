# RERconverge for H6 (genome-wide: longevity, brain, brain beyond longevity) and H8 (panel genes within Psittaculidae
# and within the rest), as fixed in results/phase0/analysis_plan.md (2026-10-28).
# Inputs: data/genome_alignments/aa/*.fa (scripts/34_align_genome.py), results/phase3_genome/{rer_traits.tsv,
# master_<tree>.nwk} (scripts/35_prepare_rer_inputs.py), optional GO gene sets data/external/go/go_bp.gmt (human
# symbols). Usage: Rscript scripts/R/rerconverge_h6_h8.R [burgio2019|smith2024]   (default burgio2019)
# Outputs: results/phase3_genome/rer_<tree>_<trait>.tsv, rer_<tree>_h8_<clade>_panel.tsv, rer_<tree>_enrichment_<trait>.tsv
suppressPackageStartupMessages({
  library(Rcpp)
  library(RERconverge)
  library(TreeTools)
  library(impute)
  library(geiger)
  library(parallel)
})

# Exact C++ implementation of RERconverge's winsorized Pearson correlation across N permulated path vectors
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

args <- commandArgs(trailingOnly = TRUE)
tree_name <- if (length(args)) args[1] else "burgio2019"
out <- "results/phase3_genome"
nperm <- as.integer(Sys.getenv("RER_PERMS", "1000"))
ncores <- as.integer(Sys.getenv("RER_CORES", "8"))
set.seed(20261028)

master <- file.path(out, paste0("master_", tree_name, ".nwk"))
genetrees <- file.path(out, paste0("genetrees_", tree_name, ".txt"))
mtree <- ape::multi2di(ape::read.tree(master), random = FALSE)

if (!file.exists(genetrees)) {
  alnfiles <- sort(list.files(path = "data/genome_alignments/aa", pattern = "\\.fa$", full.names = TRUE))
  cat("Estimating LG branch lengths for", length(alnfiles), "genes on", tree_name, "using", ncores, "cores...\n")
  est_one <- function(f) {
    g <- sub("\\.fa$", "", basename(f))
    tryCatch({
      aln <- phangorn::read.phyDat(f, type = "AA", format = "fasta")
      inboth <- intersect(names(aln), mtree$tip.label)
      if (length(inboth) < 10) return(NULL)
      gt <- ape::unroot(ape::multi2di(ape::keep.tip(mtree, inboth), random = FALSE))
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
  writeLines(lines, genetrees)
}

need_h6 <- any(!file.exists(file.path(out, sprintf("rer_%s_%s.tsv", tree_name,
                                                   c("rel_log_life_expectancy", "rel_log_brain", "rel_brain_given_le")))))
trees <- if (need_h6) readTrees(genetrees) else NULL
traits <- read.delim(file.path(out, "rer_traits.tsv"))
gmt <- "data/external/go/go_bp.gmt"
annots <- if (need_h6 && file.exists(gmt)) list(GO_BP = read.gmt(gmt)) else NULL

run_trait <- function(col, species = NULL, label = col) {
  tv <- setNames(traits[[col]], traits$species_tip)
  tv <- tv[!is.na(tv)]
  if (!is.null(species)) tv <- tv[names(tv) %in% species]
  t_use <- if (is.null(species)) trees else readTrees(genetrees, useSpecies = names(tv))
  min_sp <- if (length(tv) < 20) 8 else 10
  min_val <- min(20, max(8, floor(0.6 * length(tv))))
  pdf(NULL)
  rers <- getAllResiduals(t_use, useSpecies = names(tv), transform = "sqrt",
                          use.weights = TRUE, norm = "scale", min.sp = min_sp, min.valid = min_val)
  cp <- char2Paths(tv, t_use)
  res <- correlateWithContinuousPhenotype(rers, cp, min.sp = min_sp, winsorizeRER = 3, winsorizetrait = 3)
  mt_root <- ape::multi2di(phangorn::midpoint(t_use$masterTree), random = FALSE)
  mt_sub <- ape::multi2di(ape::keep.tip(mt_root, names(tv)), random = FALSE)
  mt_sub$edge.length[mt_sub$edge.length <= 0] <- 1e-6
  Pmat <- do.call(cbind, parallel::mclapply(1:nperm, function(k) {
    set.seed(20261028 + k)
    suppressMessages(char2Paths(simpermvec(tv, mt_sub), t_use))
  }, mc.cores = ncores))
  corRho <- fastCorPerms(rers, Pmat, min_sp, 3)
  rownames(corRho) <- rownames(rers)
  perms <- list(corRho = corRho)
  pv <- permpvalcor(res, perms)
  res$permP <- pv$permpval
  res$permStat <- pv$permstats
  res$permQ <- p.adjust(res$permP, "BH")
  res$gene <- rownames(res)
  res <- res[, c("gene", "Rho", "N", "P", "p.adj", "permP", "permQ", "permStat")]
  write.table(res, file.path(out, sprintf("rer_%s_%s.tsv", tree_name, label)), sep = "\t", quote = FALSE,
              row.names = FALSE)
  if (!is.null(annots) && is.null(species)) {
    enr <- fastwilcoxGMTall(getStat(res), annots, outputGeneVals = FALSE)
    enr_df <- enr$GO_BP
    enr_df$term <- rownames(enr_df)
    enr_df <- enr_df[, c("term", setdiff(colnames(enr_df), "term"))]
    write.table(enr_df, file.path(out, sprintf("rer_%s_enrichment_%s.tsv", tree_name, label)), sep = "\t",
                quote = FALSE, row.names = FALSE)
  }
  invisible(res)
}

# H6: genome-wide, each phenotype separately; brain beyond longevity is the mandatory analysis for shared genes
for (col in c("rel_log_life_expectancy", "rel_log_brain", "rel_brain_given_le")) {
  f_out <- file.path(out, sprintf("rer_%s_%s.tsv", tree_name, col))
  if (!file.exists(f_out)) {
    cat("Running H6 trait:", col, "on", tree_name, "...\n")
    run_trait(col)
  }
}

# H8 item 2: panel genes within Psittaculidae and within the rest (longevity)
panel <- read.delim("data/pilot_genes/canonical_genes.tsv")$gene_symbol
for (cl in c("Psittaculidae", "rest")) {
  cat("Running H8 within-clade:", cl, "on", tree_name, "...\n")
  sp <- traits$species_tip[traits$clade == cl]
  res <- run_trait("rel_log_life_expectancy", species = sp, label = paste0("h8_", cl))
  keep <- res[res$gene %in% panel, ]
  keep$panelQ <- p.adjust(keep$permP, "BH")
  keep$panel_padj_analytic <- p.adjust(keep$P, "BH")
  write.table(keep, file.path(out, sprintf("rer_%s_h8_%s_panel.tsv", tree_name, cl)), sep = "\t", quote = FALSE,
              row.names = FALSE)
}
cat("RERconverge done:", tree_name, "\n")
