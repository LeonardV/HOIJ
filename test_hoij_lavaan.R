# ============================================================
# Tests voor hoij_lavaan()
#
# (1) basisrun IJ1 + HOIJ-2 op het HS-mediatiemodel;
# (2) validatie tegen een EXACTE bootstrap met dezelfde
#     gewichtsvectoren (verschil = pure approximatiefout);
# (3) scope-checks geven nette fouten.
#
# Draaien:  Rscript test_hoij_lavaan.R
# ============================================================

suppressPackageStartupMessages(library(lavaan))
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

## ── (1) basisrun ──────────────────────────────────────────────
cat("── (1) Basisrun ──\n")
h2 <- hoij_lavaan(fit, functional = functionals,
                  B = 1000L, order = 2L, seed = 1, details = TRUE)
print(h2)
h1 <- hoij_lavaan(fit, functional = functionals,
                  B = 1000L, order = 1L, seed = 1)
stopifnot(all(is.finite(h2$results$se)), all(is.finite(h2$results$lo)))

## default: alle vrije parameters
h_all <- hoij_lavaan(fit, B = 400L, seed = 2)
stopifnot(nrow(h_all$results) == length(coef(fit)))
## IJ/bootstrap-SE's moeten dezelfde orde van grootte hebben als de
## klassieke SE's (ratio ruwweg in [0.5, 2])
se_ratio <- h_all$results$se / sqrt(diag(lavInspect(fit, "vcov")))
cat(sprintf("\nSE-ratio HOIJ-2 vs. lavaan-standaard: mediaan %.3f (bereik %.2f–%.2f)\n",
            median(se_ratio), min(se_ratio), max(se_ratio)))
stopifnot(all(se_ratio > 0.5 & se_ratio < 2))

## ── (2) exacte bootstrap met DEZELFDE gewichten ───────────────
cat("\n── (2) HOIJ-2 vs. exacte bootstrap (gedeelde gewichten) ──\n")
W  <- h2$weights
B  <- nrow(W)
N  <- nrow(HolzingerSwineford1939)
dat_mat <- as.matrix(lavInspect(fit, "data"))
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
cat(sprintf("Bootstrap: %d/%d geconvergeerd, %.1fs\n",
            sum(is.finite(boot_th[, 1])), B, proc.time()[["elapsed"]] - t0))

eval_expr <- function(txt, mat) {
  e <- parse(text = txt)[[1]]
  apply(mat, 1, function(row) eval(e, envir = as.list(row)))
}
pr <- c(0.025, 0.975)
cat(sprintf("%-10s %-7s %8s %8s %8s\n", "functional", "methode", "se", "lo", "hi"))
for (nm in names(functionals)) {
  v_h2 <- eval_expr(functionals[[nm]], h2$replicates)
  v_bt <- eval_expr(functionals[[nm]], boot_th[is.finite(boot_th[, 1]), ])
  q_h2 <- quantile(v_h2, pr, names = FALSE); q_bt <- quantile(v_bt, pr, names = FALSE)
  cat(sprintf("%-10s %-7s %8.4f %8.4f %8.4f\n", nm, "hoij2", sd(v_h2), q_h2[1], q_h2[2]))
  cat(sprintf("%-10s %-7s %8.4f %8.4f %8.4f\n", nm, "boot",  sd(v_bt), q_bt[1], q_bt[2]))
  ## approximatiefout moet klein zijn t.o.v. de intervalbreedte
  rel <- max(abs(q_h2 - q_bt)) / (q_bt[2] - q_bt[1])
  cat(sprintf("%-10s max |CI-verschil| / breedte = %.3f\n", "", rel))
  stopifnot(rel < 0.15,
            abs(sd(v_h2) / sd(v_bt) - 1) < 0.15)
}

## ── (3) scope-checks ──────────────────────────────────────────
cat("\n── (3) Scope-checks ──\n")
expect_error <- function(expr, patroon) {
  msg <- tryCatch({ expr; NULL }, error = function(e) conditionMessage(e))
  stopifnot(!is.null(msg), grepl(patroon, msg))
  cat("  OK:", msg, "\n")
}
fit_2g  <- sem("f =~ x1 + x2 + x3", data = HolzingerSwineford1939,
               group = "school")
expect_error(hoij_lavaan(fit_2g), "1 groep")
fit_uls <- sem("f =~ x1 + x2 + x3", data = HolzingerSwineford1939,
               estimator = "ULS")
expect_error(hoij_lavaan(fit_uls), "ML")
fit_eq  <- sem("f =~ x1 + a*x2 + a*x3", data = HolzingerSwineford1939)
expect_error(hoij_lavaan(fit_eq), "restricties")

cat("\nAlle tests geslaagd.\n")
