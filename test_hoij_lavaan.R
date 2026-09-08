# =====================================================================
# Tests for hoij_lavaan()
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
fit <- sem(model, data = HolzingerSwineford1939, estimator = "ML",
           meanstructure = TRUE)
stopifnot(length(coef(fit)) == 30L)

functionals <- c(ab = "a*b", psi_visual = "`visual~~visual`", psi_speed = "`speed~~speed`")


## --- (1) basic run ---------------------------------------------------
h2 <- hoij_lavaan(fit, functional = functionals, B = 1000L, order = 2L,
                  seed = 42, details = TRUE)
print(h2)
h1 <- hoij_lavaan(fit, functional = functionals, B = 1000L, order = 1L,
                  seed = 42)
print(h1)

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
        sample.mean = mu, meanstructure = TRUE, sample.cov.rescale = TRUE,
        se = "none", test = "none", h1 = FALSE, baseline = FALSE,
        check.gradient = FALSE, check.start = FALSE, check.post = FALSE,
        control = list(iter.max = 150L)),
    error = function(e) NULL)
  if (!is.null(fb) && isTRUE(lavInspect(fb, "converged")))
    boot_th[r, ] <- coef(fb, type = "free")
}
cat(sprintf("Bootstrap: %d/%d converged, %.1fs\n",
            sum(is.finite(boot_th[, 1])), B, proc.time()[["elapsed"]] - t0))

## every free parameter, against the bootstrap on the same weights
se_h2 <- apply(h2$replicates, 2, sd)
se_bt <- apply(boot_th[is.finite(boot_th[, 1]), , drop = FALSE], 2, sd)
cat(sprintf("SE ratio HOIJ-2 vs exact bootstrap, all %d parameters: %.3f-%.3f\n",
            length(se_h2), min(se_h2 / se_bt), max(se_h2 / se_bt)))
stopifnot(all(se_h2 / se_bt > 0.75), all(se_h2 / se_bt < 1.30))

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
  rel <- max(abs(q_h2 - q_bt)) / (q_bt[2] - q_bt[1])
  cat(sprintf("%-10s max |CI difference| / width = %.3f\n", "", rel))
  stopifnot(rel < 0.25, abs(sd(v_h2) / sd(v_bt) - 1) < 0.20)
}
