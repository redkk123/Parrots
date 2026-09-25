# Cross-check the Python phylogenetic code (scripts/lib/phylo.py) against ape / nlme / phylolm / phytools.
# Run from the repository root: Rscript scripts/R/validate_phase0.R
suppressPackageStartupMessages({library(ape); library(nlme); library(phylolm); library(phytools); library(jsonlite)})
tree <- read.tree("results/phase0/species_tree.nwk")
d <- read.delim("results/phase0/phenotypes.tsv")
d <- d[match(tree$tip.label, d$tree_tip), ]
rownames(d) <- d$tree_tip
py <- fromJSON("results/phase0/phase0_summary.json")
out <- list(r_version = R.version.string,
            packages = sapply(c("ape", "nlme", "phylolm", "phytools"), function(p) as.character(packageVersion(p))))

# 1. Pagel's lambda PGLS
g <- gls(log_life_expectancy ~ log_body_mass, data = d, method = "ML",
         correlation = corPagel(0.5, tree, form = ~tree_tip, fixed = FALSE))
p_lam <- phylolm(log_life_expectancy ~ log_body_mass, data = d, phy = tree, model = "lambda")
p_bm <- phylolm(log_life_expectancy ~ log_body_mass, data = d, phy = tree, model = "BM")
p_ou <- phylolm(log_life_expectancy ~ log_body_mass, data = d, phy = tree, model = "OUfixedRoot")
out$pgls <- data.frame(
  implementation = c("python lib/phylo.py", "nlme::gls corPagel (ML)", "phylolm lambda"),
  lambda = c(py$allometry$lambda, coef(g$modelStruct$corStruct, unconstrained = FALSE), p_lam$optpar),
  slope = c(py$allometry$slope, coef(g)[2], coef(p_lam)[2]),
  slope_se = c(py$allometry$slope_se, sqrt(diag(vcov(g)))[2], sqrt(diag(p_lam$vcov))[2]))
out$phylolm_aic <- c(BM = AIC(p_bm), lambda = AIC(p_lam), OU = AIC(p_ou))
out$phylolm_ou_alpha_per_tree_height <- p_ou$optpar * max(node.depth.edgelength(tree))

# 2. Ancestral states of relative life expectancy: phytools::fastAnc vs Python BM estimates
pa <- read.delim("results/phase0/ancestral_states_rel_logLE.tsv")
tips <- pa[pa$is_tip %in% c("True", TRUE), ]
x <- setNames(tips$est_rel_logLE, tips$label)[tree$tip.label]
fa <- fastAnc(tree, x, vars = TRUE)
ntip <- Ntip(tree)
keys <- sapply(ntip + seq_len(tree$Nnode), function(nd) {
  lab <- tree$tip.label[phangorn_free_desc <- (function(nd) { # descendant tips without phangorn
    kids <- nd; tipset <- integer(0)
    while (length(kids)) { k <- kids[1]; kids <- kids[-1]
      if (k <= ntip) tipset <- c(tipset, k) else kids <- c(kids, tree$edge[tree$edge[, 1] == k, 2]) }
    tipset })(nd)]
  paste0(length(lab), "|", min(lab), "|", max(lab))
})
stopifnot(!anyDuplicated(keys))
cmp <- merge(data.frame(clade_key = keys, r_est = fa$ace, r_var = fa$var),
             pa[, c("clade_key", "est_rel_logLE", "se")], by = "clade_key")
out$ancestral <- list(n_nodes_matched = nrow(cmp), n_nodes = tree$Nnode,
                      max_abs_diff_estimate = max(abs(cmp$r_est - cmp$est_rel_logLE)),
                      cor_estimate = cor(cmp$r_est, cmp$est_rel_logLE),
                      cor_se = cor(sqrt(cmp$r_var), cmp$se))
dir.create("results/phase0b", showWarnings = FALSE)
write(toJSON(out, auto_unbox = TRUE, pretty = TRUE, digits = 6), "results/phase0b/validation_python_vs_R.json")
cat(toJSON(out, auto_unbox = TRUE, pretty = TRUE, digits = 5), "\n")
