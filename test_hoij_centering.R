# Focused regression test: for a saturated covariance model, the exact
# weighted ML covariance is quadratic in fixed-total observation weights.
# HOIJ-2 must therefore reproduce it up to numerical derivative error.
# Run from the repository root with Rscript test_hoij_centering.R.

suppressPackageStartupMessages(library(lavaan))
source("hoij_core.R")
source("hoij_lavaan.R")

set.seed(2718)
dat <- data.frame(x = 3 + rnorm(80), y = -2 + rnorm(80))
dat$y <- 0.7 * dat$x + dat$y
model <- 'x ~~ x + y
          y ~~ y'

answers <- list()
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
  answers[[as.character(with_means)]] <- ans
}

## Both routes use identical weights and must agree on the covariances.
cov_names <- colnames(answers$`FALSE`$replicates)
stopifnot(identical(answers$`FALSE`$weights, answers$`TRUE`$weights),
          max(abs(answers$`FALSE`$replicates -
                  answers$`TRUE`$replicates[, cov_names])) < 5e-5)

## Exercise the actual simulation bootstrap helper without running the
## simulation or loading its optional VITA dependencies.
helper_env <- new.env(parent = globalenv())
for (expr in parse("03_simulation_study.R")) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.symbol(expr[[2]]) && as.character(expr[[2]]) %in%
      c("make_boot_partable", "fit_boot_cov")) eval(expr, helper_env)
}
PT <- helper_env$make_boot_partable(fit)
for (b in 1:5) {
  w <- answers$`TRUE`$weights[b, ]
  dat_b <- dat[rep.int(seq_len(nrow(dat)), w), , drop = FALSE]
  from_moments <- helper_env$fit_boot_cov(
    PT, cov(dat_b), nrow(dat_b), mu_b = colMeans(dat_b))
  from_data <- sem(model, data = dat_b, estimator = "ML",
                   meanstructure = TRUE, se = "none")
  stopifnot(!is.null(from_moments), lavInspect(from_data, "converged"),
            identical(names(coef(from_moments)), names(coef(from_data))),
            max(abs(coef(from_moments) - coef(from_data))) < 1e-5)
}
cat("Joint and profiled routes agree; simulation moment refits match raw data.\n")
