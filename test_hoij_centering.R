# Focused regression test: for a saturated covariance model, the exact
# weighted ML covariance is quadratic in fixed-total observation weights.
# HOIJ-2 must therefore reproduce it up to numerical derivative error.
# Run from the repository root with Rscript test_hoij_centering.R.

suppressPackageStartupMessages(library(lavaan))
source("hoij_core.R")
source("hoij_lavaan.R")

set.seed(2718)
dat <- data.frame(x = rnorm(80), y = rnorm(80))
dat$y <- 0.7 * dat$x + dat$y
model <- 'x ~~ x + y
          y ~~ y'

for (with_means in c(FALSE, TRUE)) {
  fit <- sem(model, data = dat, estimator = "ML",
             meanstructure = with_means)
  ans <- hoij_lavaan(fit, B = 100L, order = 2L, seed = 19,
                     details = TRUE)
  X <- as.matrix(lavInspect(fit, "data"))
  n <- nrow(X)
  exact <- matrix(NA_real_, nrow(ans$weights), ncol(ans$replicates),
                   dimnames = dimnames(ans$replicates))
  for (b in seq_len(nrow(ans$weights))) {
    w <- ans$weights[b, ]
    mu <- colSums(X * w) / n
    Xc <- sweep(X, 2L, mu, "-")
    S <- crossprod(Xc * w, Xc) / n
    exact[b, "x~~x"] <- S["x", "x"]
    exact[b, "x~~y"] <- S["x", "y"]
    exact[b, "y~~y"] <- S["y", "y"]
    if (with_means) {
      exact[b, "x~1"] <- mu["x"]
      exact[b, "y~1"] <- mu["y"]
    }
  }
  err <- max(abs(ans$replicates - exact))
  cat(sprintf("meanstructure = %s, maximum absolute error = %.3g\n",
              with_means, err))
  stopifnot(is.finite(err), err < 5e-5)
}
