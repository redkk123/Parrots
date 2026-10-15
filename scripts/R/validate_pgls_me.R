# Reference values for scripts/lib/phylo.py::fit_pgls_me (PGLS with known per-tip measurement error).
# Reads data/external/me_validation/{tree.nwk,data.tsv} written by scripts/19_measurement_error.py and writes
# data/external/me_validation/r_results.json.
#   (a) BM + ME regression, ML:  phytools::pgls.SEy(y ~ x, corClass = corBrownian, se = se, method = "ML")
#   (b) lambda + ME, intercept only:  phytools::phylosig(method = "lambda", se = se)
suppressPackageStartupMessages({library(ape); library(phytools); library(nlme); library(jsonlite)})
d <- read.delim("data/external/me_validation/data.tsv")
tr <- read.tree("data/external/me_validation/tree.nwk")
tr$edge.length <- tr$edge.length / max(node.depth.edgelength(tr))   # unit height, as in phylo.py
rownames(d) <- d$tip
se <- setNames(d$se, d$tip)
fa <- pgls.SEy(y ~ x, data = d, corClass = corBrownian, tree = tr, se = se, method = "ML")
y <- setNames(d$y, d$tip)
fb <- phylosig(tr, y, method = "lambda", se = se)
# (c) direct ML, written independently of phylo.py: V = sig2 * C_lambda + diag(se^2), beta by GLS.
C <- vcv(tr)[d$tip, d$tip]; X <- cbind(1, d$x); E <- diag(d$se^2)
nll <- function(p, lam_free) {
  s2 <- exp(p[1]); lam <- if (lam_free) p[2] else 1
  Cl <- C * lam; diag(Cl) <- diag(C)
  V <- s2 * Cl + E; Vi <- solve(V)
  b <- solve(t(X) %*% Vi %*% X, t(X) %*% Vi %*% d$y); r <- d$y - X %*% b
  0.5 * (nrow(X) * log(2 * pi) + as.numeric(determinant(V)$modulus) + as.numeric(t(r) %*% Vi %*% r))
}
fit_direct <- function(lam_free) {
  o <- if (lam_free) optim(c(log(0.1), 0.5), nll, lam_free = TRUE, method = "L-BFGS-B",
                           lower = c(-12, 0), upper = c(5, 1)) else
         optim(log(0.1), nll, lam_free = FALSE, method = "Brent", lower = -12, upper = 5)
  s2 <- exp(o$par[1]); lam <- if (lam_free) o$par[2] else 1
  Cl <- C * lam; diag(Cl) <- diag(C); Vi <- solve(s2 * Cl + E)
  A <- solve(t(X) %*% Vi %*% X); b <- A %*% t(X) %*% Vi %*% d$y
  list(beta = as.numeric(b), se = sqrt(diag(A)), lambda = lam, sig2 = s2, loglik = -o$value)
}
out <- list(
  # pgls.SEy re-estimates a free residual scale in gls(), i.e. it multiplies the known se^2 by an estimated
  # factor: a different model, kept for the record only.
  pgls_SEy_different_model = list(beta = unname(coef(fa)), se = unname(sqrt(diag(vcov(fa)))),
                                  loglik = as.numeric(logLik(fa))),
  bm_me = fit_direct(FALSE), lambda_me_regression = fit_direct(TRUE),
  lambda_me = list(lambda = fb$lambda, sig2 = fb$sig2, loglik = fb$logL))
write(toJSON(out, auto_unbox = TRUE, digits = 10, pretty = TRUE), "data/external/me_validation/r_results.json")
print(out)
