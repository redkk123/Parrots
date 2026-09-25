# Detect shifts in the optimum of relative longevity along the tree (scOU, PhylogeneticEM;
# Bastide et al. 2017 JRSS-B, 2018 Syst Biol). Exploratory: residuals are treated as data,
# and tip measurement error is not modelled.
# Run from the repository root with an R that has PhylogeneticEM:
#   Rscript scripts/R/regime_shifts.R
suppressPackageStartupMessages({library(ape); library(phytools); library(PhylogeneticEM); library(jsonlite)})
full <- read.tree("results/phase0/supertree_burgio2019_full.nwk")
m <- read.delim("results/phase0b/species_data_matrix_with_residuals.tsv")
outcomes <- c(rel_LE = "rel_log_life_expectancy", rel_maxL_ssadj = "rel_log_max_longevity_ssadj")
set.seed(20260926)
summary <- list(PhylogeneticEM = as.character(packageVersion("PhylogeneticEM")))

clade_label <- function(tr, edge) {
  node <- tr$edge[edge, 2]
  tips <- if (node <= Ntip(tr)) tr$tip.label[node] else extract.clade(tr, node)$tip.label
  gen <- sort(unique(sub("_.*", "", tips)))
  list(n_tips = length(tips), genera = paste(if (length(gen) > 6) c(gen[1:6], "...") else gen, collapse = ";"),
       tips = paste(sort(tips), collapse = ";"))
}

for (nm in names(outcomes)) {
  col <- outcomes[[nm]]
  d <- m[!is.na(m[[col]]), ]
  tr <- keep.tip(full, d$tree_tip)
  tr <- multi2di(tr, random = FALSE)                      # polytomies -> zero-length edges
  tr$edge.length[tr$edge.length <= 0] <- 1e-6
  tr <- force.ultrametric(tr, method = "extend", message = FALSE)
  h <- max(node.depth.edgelength(tr)); tr$edge.length <- tr$edge.length / h   # unit height
  y <- matrix(d[[col]][match(tr$tip.label, d$tree_tip)], nrow = 1, dimnames = list(nm, tr$tip.label))
  t0 <- Sys.time()
  fit <- PhyloEM(phylo = tr, Y_data = y, process = "scOU", random.root = TRUE, stationary.root = TRUE,
                 K_max = 12, parallel_alpha = TRUE, Ncores = 4, method.selection = c("LINselect", "DDSE", "Djump"),
                 progress.bar = FALSE)
  el <- difftime(Sys.time(), t0, units = "mins")
  saveRDS(fit, sprintf("data/processed/phyloem_%s.rds", nm))
  sel <- list()
  # univariate names: LINselect -> "BGHuni", DDSE -> "DDSE_BM1", Djump -> "Djump_BM1"
  for (meth in c("BGHuni", "DDSE_BM1", "Djump_BM1")) {
    pr <- tryCatch(params_process(fit, method.selection = meth), error = function(e) NULL)
    if (is.null(pr)) { sel[[meth]] <- list(error = "criterion not available"); next }
    sh <- pr$shifts
    rows <- if (length(sh$edges)) lapply(seq_along(sh$edges), function(i) c(
      list(edge = sh$edges[i], shift_value = sh$values[i]), clade_label(tr, sh$edges[i]))) else list()
    sel[[meth]] <- list(K = length(sh$edges), alpha_per_tree_height = unname(pr$selection.strength),
                        half_life_myr = log(2) / unname(pr$selection.strength) * h,
                        root_optimum = unname(pr$optimal.value), shifts = rows)
  }
  # likelihood profile over K (at the selected alpha)
  prof <- data.frame(K = 0:12, loglik = sapply(0:12, function(k) {
    v <- tryCatch(fit$alpha_max$results_summary$log_likelihood[fit$alpha_max$results_summary$K_try == k],
                  error = function(e) NA); if (length(v)) v[1] else NA }))
  summary[[nm]] <- list(n_species = Ntip(tr), minutes = as.numeric(el), selection = sel,
                        loglik_by_K = prof)
  png(sprintf("results/phase0b/figures/fig4_regime_shifts_%s.png", nm), width = 1800, height = 3600, res = 200)
  plot(fit, method.selection = "BGHuni", show.tip.label = TRUE, label_cex = 0.35, value_in_box = TRUE,
       shifts_bg = "white", show_axis_traits = TRUE)
  dev.off()
  cat(nm, ": K(BGHuni) =", sel$BGHuni$K, " K(DDSE) =", sel$DDSE_BM1$K, " minutes =", round(el, 1), "\n")
}
write(toJSON(summary, auto_unbox = TRUE, pretty = TRUE, digits = 5), "results/phase0b/regime_shifts.json")
