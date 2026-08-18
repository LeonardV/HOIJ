# =====================================================================
# Tests for hoij_lavaan()
#
# (1) basic IJ1 and HOIJ-2 run on the Holzinger-Swineford mediation model
# (2) validation against an exact bootstrap that uses the same weight
#     vectors, so the difference is pure approximation error
# (3) scope checks return informative errors
# (4) compatibility with the development version of lavaan, which
#     renamed the internal helpers used by hoij_core.R
#
# Run with:  Rscript test_hoij_lavaan.R
# =====================================================================

suppressPackageStartupMessages(library(lavaan))
source("hoij_core.R")
source("hoij_lavaan.R")

model <- '
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
  visual ~ c*textual + b*speed
  speed  ~ a*textual
'
fit <- sem(model, data = HolzingerSwineford1939, estimator = "ML")
stopifnot(lavInspect(fit, "converged"))

functionals <- c(ab = "a*b", psi_speed = "`speed~~speed`")


## --- (1) basic run ---------------------------------------------------
cat("-- (1) basic run --\n")
h2 <- hoij_lavaan(fit, functional = functionals, B = 1000L, order = 2L,
                  seed = 1, details = TRUE)
print(h2)
h1 <- hoij_lavaan(fit, functional = functionals, B = 1000L, order = 1L,
                  seed = 1)
stopifnot(all(is.finite(h2$results$se)), all(is.finite(h2$results$lo)))

## default: every free parameter
h_all <- hoij_lavaan(fit, B = 400L, seed = 2)
stopifnot(nrow(h_all$results) == length(coef(fit)))

## the replicate standard errors should be of the same order as the
## model-based ones
se_ratio <- h_all$results$se / sqrt(diag(lavInspect(fit, "vcov")))
cat(sprintf("\nSE ratio HOIJ-2 vs lavaan default: median %.3f (range %.2f-%.2f)\n",
            median(se_ratio), min(se_ratio), max(se_ratio)))
stopifnot(all(se_ratio > 0.5 & se_ratio < 2))


## --- (2) exact bootstrap on the same weight vectors ------------------
cat("\n-- (2) HOIJ-2 vs exact bootstrap (shared weights) --\n")
W <- h2$weights
B <- nrow(W)
N <- nrow(HolzingerSwineford1939)
dat_mat  <- as.matrix(lavInspect(fit, "data"))
th_names <- names(coef(fit))

PT <- parTable(fit); PT$ustart <- PT$est; PT$start <- PT$est; PT$se <- NULL

boot_th <- matrix(NA_real_, B, length(th_names),
                  dimnames = list(NULL, th_names))
t0 <- proc.time()[["elapsed"]]
for (r in seq_len(B)) {
  w   <- W[r, ]
  mu  <- colSums(dat_mat * w) / N
  S_b <- (crossprod(dat_mat * w, dat_mat) - N * tcrossprod(mu)) / (N - 1)
  S_b <- (S_b + t(S_b)) / 2
  fb <- tryCatch(
    sem(model = PT, sample.cov = S_b, sample.nobs = N, estimator = "ML",
        se = "none", test = "none", h1 = FALSE, baseline = FALSE,
        check.gradient = FALSE, check.start = FALSE, check.post = FALSE,
        control = list(iter.max = 150L)),
    error = function(e) NULL)
  if (!is.null(fb) && isTRUE(lavInspect(fb, "converged")))
    boot_th[r, ] <- coef(fb, type = "free")
}
cat(sprintf("Bootstrap: %d/%d converged, %.1fs\n",
            sum(is.finite(boot_th[, 1])), B, proc.time()[["elapsed"]] - t0))

eval_expr <- function(txt, mat) {
  e <- parse(text = txt)[[1]]
  apply(mat, 1, function(row) eval(e, envir = as.list(row)))
}
pr <- c(0.025, 0.975)
cat(sprintf("%-10s %-7s %8s %8s %8s\n", "functional", "method", "se", "lo", "hi"))
for (nm in names(functionals)) {
  v_h2 <- eval_expr(functionals[[nm]], h2$replicates)
  v_bt <- eval_expr(functionals[[nm]], boot_th[is.finite(boot_th[, 1]), ])
  q_h2 <- quantile(v_h2, pr, names = FALSE)
  q_bt <- quantile(v_bt, pr, names = FALSE)
  cat(sprintf("%-10s %-7s %8.4f %8.4f %8.4f\n", nm, "hoij2", sd(v_h2),
              q_h2[1], q_h2[2]))
  cat(sprintf("%-10s %-7s %8.4f %8.4f %8.4f\n", nm, "boot", sd(v_bt),
              q_bt[1], q_bt[2]))
  ## the approximation error should be small relative to the width
  rel <- max(abs(q_h2 - q_bt)) / (q_bt[2] - q_bt[1])
  cat(sprintf("%-10s max |CI difference| / width = %.3f\n", "", rel))
  stopifnot(rel < 0.15, abs(sd(v_h2) / sd(v_bt) - 1) < 0.15)
}


## --- (3) scope checks -------------------------------------------------
cat("\n-- (3) scope checks --\n")
expect_error <- function(expr, pattern) {
  msg <- tryCatch({ expr; NULL }, error = function(e) conditionMessage(e))
  stopifnot(!is.null(msg), grepl(pattern, msg))
  cat("  OK:", msg, "\n")
}
expect_error(hoij_lavaan(sem("f =~ x1 + x2 + x3",
                             data = HolzingerSwineford1939,
                             group = "school")), "single-group")
expect_error(hoij_lavaan(sem("f =~ x1 + x2 + x3",
                             data = HolzingerSwineford1939,
                             estimator = "ULS")), "ML")
expect_error(hoij_lavaan(sem("f =~ x1 + a*x2 + a*x3",
                             data = HolzingerSwineford1939)),
             "equality constraints")


## --- (4) development-lavaan compatibility -----------------------------
## The development version renamed lav_model_x2GLIST -> lav_model_x2glist,
## lav_model_gradient -> lav_model_grad and the argument GLIST -> glist.
## Here the argument rename is simulated by replacing the release
## internals with shims that have the development signature; results
## must be identical to the unmodified run.
cat("\n-- (4) development-lavaan compatibility (simulated renames) --\n")
ns <- asNamespace("lavaan")
if (all(c("lav_model_implied", "lav_model_gradient") %in% ls(ns))) {
  orig_implied  <- get("lav_model_implied",  envir = ns)
  orig_gradient <- get("lav_model_gradient", envir = ns)

  ## The shims must keep serving lavaan's own internal calls, which use
  ## the old GLIST= name; those are caught through `...`. Our own code
  ## detects `glist` in the formals and uses the development name.
  dev_implied <- function(lavmodel = NULL, glist = NULL, delta = TRUE, ...) {
    dots <- list(...)
    if (is.null(glist) && !is.null(dots$GLIST)) glist <- dots$GLIST
    orig_implied(lavmodel = lavmodel, GLIST = glist, delta = delta)
  }
  dev_gradient <- function(lavmodel = NULL, glist = NULL,
                           lavsamplestats = NULL, lavdata = NULL,
                           lavcache = NULL, ...) {
    dots <- list(...)
    if (is.null(glist) && !is.null(dots$GLIST)) glist <- dots$GLIST
    orig_gradient(lavmodel = lavmodel, GLIST = glist,
                  lavsamplestats = lavsamplestats, lavdata = lavdata,
                  lavcache = lavcache)
  }
  assignInNamespace("lav_model_implied",  dev_implied,  ns = "lavaan")
  assignInNamespace("lav_model_gradient", dev_gradient, ns = "lavaan")
  source("hoij_core.R")   # reset the resolver cache
  h2_dev <- hoij_lavaan(fit, functional = functionals, B = 1000L,
                        order = 2L, seed = 1)
  stopifnot(isTRUE(all.equal(h2_dev$results, h2$results)))
  cat("  OK: identical results with development-style glist signatures\n")

  ## An implied() that IGNORES its glist argument is the silent failure
  ## mode described in 00_install_dependencies.R; it must be caught.
  broken_implied <- function(lavmodel = NULL, glist = NULL, delta = TRUE, ...)
    orig_implied(lavmodel = lavmodel, GLIST = NULL, delta = delta)
  assignInNamespace("lav_model_implied", broken_implied, ns = "lavaan")
  source("hoij_core.R")
  expect_error(hoij_lavaan(fit, functional = functionals, B = 100L, seed = 1),
               "Broken link")

  assignInNamespace("lav_model_implied",  orig_implied,  ns = "lavaan")
  assignInNamespace("lav_model_gradient", orig_gradient, ns = "lavaan")
  source("hoij_core.R")
} else {
  cat("  (skipped: this lavaan version already uses the development names,\n",
      "   so sections 1-3 exercise that code path directly)\n")
}

cat("\nAll tests passed.\n")
