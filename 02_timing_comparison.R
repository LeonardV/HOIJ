# =====================================================================
# 02_timing_comparison.R
#
# Companion code for:
#   Vanbrabant, L., & Rosseel, Y. Approximating percentile bootstrap
#   confidence intervals in SEM without repeated refitting: A tutorial
#   on the second-order infinitesimal jackknife.
#
# Reproduces Table 3: timing decomposition of the exact and the
# approximate bootstrap for two models fitted to N = 500 observations,
#   (i)  the D = 21 latent mediation model of Section 3, and
#   (ii) an orthogonal bifactor model for the same nine indicators
#        (D = 27), which is deliberately harder to refit.
#
# The comparison separates the one-time derivative setup from the
# per-replicate cost. Per-replicate costs are measured on N_TIMING = 200
# shared weight vectors and extrapolated as
#   approximate: setup + B * cost per replicate
#   exact      : B * mean refit time
# which is exactly the arithmetic reported in the table. Wall-clock
# numbers are hardware specific; the scaling pattern (one fit plus
# derivatives, versus B full refits) is the point.
#
# Run 00_install_dependencies.R once before this script.
# =====================================================================

source("hoij_core.R")
suppressPackageStartupMessages(library(lavaan))

## --- settings --------------------------------------------------------
N_EX         <- 500     # sample size of the example data set
N_TIMING     <- 200     # replicates used for the per-replicate timings
N_FIT_TIMING <- 30      # repeated fits for a stable single-fit timing
B_TARGETS    <- c(1000, 5000)
KAPPA        <- 0.5     # trust-region damping, as in the other scripts
SEED_MED     <- 20260706
SEED_BIF     <- 20260810

out_dir <- "hoij_timing_output"
if (!dir.exists(out_dir)) dir.create(out_dir)

cat("lavaan", as.character(packageVersion("lavaan")), "\n")
hoij_selftest()


# ---------------------------------------------------------------------
# Models
#
# The bifactor model is not specified in the article text. We use the
# standard orthogonal bifactor structure on the same nine indicators: a
# general factor g plus the three original group factors, all latent
# covariances fixed to zero and all factor variances fixed to one
# (std.lv = TRUE). That gives 18 loadings + 9 residual variances = 27
# free parameters.
# ---------------------------------------------------------------------
model_med <- '
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9

  speed   ~ a*textual
  visual  ~ b*speed + c*textual
'

model_bifactor <- '
  g       =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9

  g ~~ 0*visual
  g ~~ 0*textual
  g ~~ 0*speed
  visual  ~~ 0*textual
  visual  ~~ 0*speed
  textual ~~ 0*speed
'


# ---------------------------------------------------------------------
# Median wall-clock time of a single warm-start ML fit
# ---------------------------------------------------------------------
time_repeated_fit <- function(model_syntax, dat, start_fit, std_lv, n_fit) {
  times <- rep(NA_real_, n_fit)
  for (i in seq_len(n_fit)) {
    gc(FALSE)
    tt <- system.time({
      f <- tryCatch(sem(model_syntax, data = dat, se = "none",
                        estimator = "ML", std.lv = std_lv, start = start_fit),
                    error = function(e) NULL)
    })["elapsed"]
    if (!is.null(f) && lavInspect(f, "converged")) times[i] <- as.numeric(tt)
  }
  out <- c(median = median(times, na.rm = TRUE),
           mean   = mean(times, na.rm = TRUE),
           min    = min(times, na.rm = TRUE),
           max    = max(times, na.rm = TRUE),
           n_fail = sum(is.na(times)))
  if (!is.finite(out["median"])) stop("no converged repeated fits")
  out
}


# ---------------------------------------------------------------------
# Full timing decomposition for one model
# ---------------------------------------------------------------------
time_one_model <- function(model_syntax, D_expected, label, seed_data,
                           std_lv = FALSE) {
  cat(sprintf("\n===== %s (D = %d) =====\n", label, D_expected))

  ## Population = Holzinger-Swineford ML estimates. No effect override:
  ## this example is about computational cost, not about skewness.
  fit_pop <- sem(model_syntax, data = HolzingerSwineford1939, se = "none",
                 std.lv = std_lv)
  stopifnot(lavInspect(fit_pop, "converged"))

  set.seed(seed_data)
  dat <- simulateData(parTable(fit_pop), sample.nobs = N_EX)

  fit <- sem(model_syntax, data = dat, se = "none", estimator = "ML",
             std.lv = std_lv)
  stopifnot(lavInspect(fit, "converged"))

  theta0 <- coef(fit, type = "free")
  D <- length(theta0); N <- nrow(dat)
  stopifnot(D == D_expected)

  fit_timing <- time_repeated_fit(model_syntax, dat, fit, std_lv,
                                  N_FIT_TIMING)
  t_fit <- unname(fit_timing["median"])
  cat(sprintf("N = %d; median warm-start fit = %.4f s (%d fits, %d failed)\n",
              N, t_fit, N_FIT_TIMING, as.integer(fit_timing["n_fail"])))

  ## --- setup, timed in the three blocks reported in the table --------
  t_scores <- system.time({
    Scores <- lavScores(fit, scaling = TRUE)
    H.inv  <- lavTech(fit, "inverted.information.observed")
    H_obs  <- lavTech(fit, "information.observed")
  })["elapsed"]

  grad_F <- make_grad_F(fit)
  t_cal <- system.time({
    cal <- calibrate_alpha(grad_F, theta0, H_obs)
  })["elapsed"]
  if (!is.finite(cal$alpha) || cal$spread > 0.1)
    stop(sprintf("[%s] alpha calibration failed (spread = %.3g)",
                 label, cal$spread))

  t_J <- system.time({
    J_all <- compute_all_J(fit, theta0)
  })["elapsed"]

  ## compute_T_tensor_grad() uses exactly 2 * D^2 gradient evaluations
  t_T <- system.time({
    T_arr <- compute_T_tensor_grad(grad_F, theta0, cal$alpha)
  })["elapsed"]

  ## The scale calibration is a fixed cost of obtaining Jhat^-1 on the
  ## log-likelihood scale, so it is reported with the scores block.
  t_scores_hinv <- as.numeric(t_scores + t_cal)
  t_setup <- t_scores_hinv + as.numeric(t_J) + as.numeric(t_T)
  cat(sprintf("setup: scores+Jhat^-1 = %.4f s | J_i = %.4f s | Khat = %.4f s | total = %.4f s\n",
              t_scores_hinv, as.numeric(t_J), as.numeric(t_T), t_setup))

  ## --- per-replicate cost: approximate bootstrap ---------------------
  set.seed(seed_data + 1)
  W_counts <- t(rmultinom(N_TIMING, size = N, prob = rep(1 / N, N)))
  dW <- W_counts - 1L

  t_loop <- system.time({
    ij1 <- ij1_replicates(theta0, Scores, H.inv, dW)
    hoij2_replicates(theta0, ij1$C, dW, H.inv, J_all, T_arr, kappa = KAPPA)
  })["elapsed"]
  t_rep_approx <- as.numeric(t_loop) / N_TIMING
  cat(sprintf("approximate bootstrap: %.6f s per replicate (n = %d)\n",
              t_rep_approx, N_TIMING))

  ## --- per-replicate cost: exact bootstrap, same weight vectors ------
  t_refit <- rep(NA_real_, N_TIMING)
  n_fail <- 0L
  for (bb in seq_len(N_TIMING)) {
    dat_b <- dat[rep.int(seq_len(N), W_counts[bb, ]), , drop = FALSE]
    tt <- system.time({
      fit_b <- tryCatch(sem(model_syntax, data = dat_b, se = "none",
                            estimator = "ML", std.lv = std_lv, start = fit),
                        error = function(e) NULL)
    })["elapsed"]
    if (!is.null(fit_b) && lavInspect(fit_b, "converged")) {
      t_refit[bb] <- as.numeric(tt)
    } else {
      n_fail <- n_fail + 1L
    }
  }
  if (sum(!is.na(t_refit)) < 0.5 * N_TIMING)
    stop(sprintf("[%s] too many failed refits (%d/%d) for a reliable timing",
                 label, n_fail, N_TIMING))

  ## Failed refits are excluded from the mean, so the reported exact
  ## bootstrap cost is conservative.
  t_rep_exact <- mean(t_refit, na.rm = TRUE)
  cat(sprintf("exact bootstrap: %.6f s per replicate (%d/%d converged)\n",
              t_rep_exact, N_TIMING - n_fail, N_TIMING))

  ## --- totals, speed-up and break-even -------------------------------
  tot_approx <- t_setup + B_TARGETS * t_rep_approx
  tot_exact  <- B_TARGETS * t_rep_exact
  Bstar <- if (t_rep_exact > t_rep_approx)
    ceiling(t_setup / (t_rep_exact - t_rep_approx)) else Inf

  list(label = label, D = D, fit_timing = fit_timing, t_fit = t_fit,
       t_scores_hinv = t_scores_hinv, t_J = as.numeric(t_J),
       t_T = as.numeric(t_T), t_setup = t_setup,
       t_rep_approx = t_rep_approx,
       tot_approx_1000 = tot_approx[1], tot_approx_5000 = tot_approx[2],
       t_rep_exact = t_rep_exact,
       tot_exact_1000 = tot_exact[1], tot_exact_5000 = tot_exact[2],
       speedup_1000 = tot_exact[1] / tot_approx[1], Bstar = Bstar,
       n_fail_exact = n_fail)
}


# ---------------------------------------------------------------------
# Run both models
# ---------------------------------------------------------------------
res_med <- time_one_model(model_med, D_expected = 21, label = "Mediation",
                          seed_data = SEED_MED, std_lv = FALSE)
res_bif <- time_one_model(model_bifactor, D_expected = 27, label = "Bifactor",
                          seed_data = SEED_BIF, std_lv = TRUE)


# ---------------------------------------------------------------------
# Table 3
# ---------------------------------------------------------------------
rows <- list(
  c("Median warm-start ML fit",            "t_fit",           4),
  c("casewise scores and Jhat^-1",         "t_scores_hinv",   4),
  c("casewise curvatures J_i",             "t_J",             4),
  c("third-derivative array Khat",         "t_T",             4),
  c("setup total",                         "t_setup",         4),
  c("per replication (approximate)",       "t_rep_approx",    5),
  c("total, B = 1,000 (approximate)",      "tot_approx_1000", 3),
  c("total, B = 5,000 (approximate)",      "tot_approx_5000", 3),
  c("per replication (mean refit)",        "t_rep_exact",     5),
  c("total, B = 1,000 (exact)",            "tot_exact_1000",  3),
  c("total, B = 5,000 (exact)",            "tot_exact_5000",  3),
  c("Speed-up factor at B = 1,000",        "speedup_1000",    1),
  c("Break-even B*",                       "Bstar",           0))

tab_timing <- data.frame(
  Row           = vapply(rows, `[`, character(1), 1),
  Mediation_D21 = vapply(rows, function(r) as.numeric(res_med[[r[2]]]),
                         numeric(1)),
  Bifactor_D27  = vapply(rows, function(r) as.numeric(res_bif[[r[2]]]),
                         numeric(1)))

cat("\n===== Table 3 =====\n")
print(tab_timing, row.names = FALSE, digits = 4)

fmt <- function(x, d) formatC(as.numeric(x), digits = d, format = "f")
tex <- file.path(out_dir, "tab_timing.tex")
con <- file(tex, "w")
writeLines(c(
  sprintf("\\emph{Approximate bootstrap: setup (once)} & & \\\\"),
  sprintf("\\quad casewise scores and $\\Jhat^{-1}$ & %s & %s \\\\",
          fmt(res_med$t_scores_hinv, 4), fmt(res_bif$t_scores_hinv, 4)),
  sprintf("\\quad casewise curvatures $J_i$ & %s & %s \\\\",
          fmt(res_med$t_J, 4), fmt(res_bif$t_J, 4)),
  sprintf("\\quad third-derivative array $\\Khat$ & %s & %s \\\\",
          fmt(res_med$t_T, 4), fmt(res_bif$t_T, 4)),
  sprintf("\\quad setup total & %s & %s \\\\",
          fmt(res_med$t_setup, 4), fmt(res_bif$t_setup, 4)),
  "\\addlinespace",
  sprintf("\\emph{Approximate bootstrap: replication loop} & & \\\\"),
  sprintf("\\quad per replication & %s & %s \\\\",
          fmt(res_med$t_rep_approx, 5), fmt(res_bif$t_rep_approx, 5)),
  sprintf("\\quad total, $B=1{,}000$ & %s & %s \\\\",
          fmt(res_med$tot_approx_1000, 3), fmt(res_bif$tot_approx_1000, 3)),
  sprintf("\\quad total, $B=5{,}000$ & %s & %s \\\\",
          fmt(res_med$tot_approx_5000, 3), fmt(res_bif$tot_approx_5000, 3)),
  "\\addlinespace",
  sprintf("\\emph{Exact bootstrap} & & \\\\"),
  sprintf("\\quad per replication (mean refit) & %s & %s \\\\",
          fmt(res_med$t_rep_exact, 5), fmt(res_bif$t_rep_exact, 5)),
  sprintf("\\quad total, $B=1{,}000$ & %s & %s \\\\",
          fmt(res_med$tot_exact_1000, 3), fmt(res_bif$tot_exact_1000, 3)),
  sprintf("\\quad total, $B=5{,}000$ & %s & %s \\\\",
          fmt(res_med$tot_exact_5000, 3), fmt(res_bif$tot_exact_5000, 3)),
  "\\addlinespace",
  sprintf("Speed-up factor at $B=1{,}000$ & %s & %s \\\\",
          fmt(res_med$speedup_1000, 1), fmt(res_bif$speedup_1000, 1)),
  sprintf("Break-even $B^{*}$ & %s & %s \\\\",
          format(res_med$Bstar), format(res_bif$Bstar))), con)
close(con)
cat("  written:", tex, "\n")

cat(sprintf("\nFailed exact refits during timing: mediation %d/%d, bifactor %d/%d\n",
            res_med$n_fail_exact, N_TIMING, res_bif$n_fail_exact, N_TIMING))

stamp <- format(Sys.time(), "%Y%m%d_%H%M")
write.csv(tab_timing, file.path(out_dir, sprintf("tab_timing_%s.csv", stamp)),
          row.names = FALSE)
saveRDS(list(mediation = res_med, bifactor = res_bif, N = N_EX,
             N_TIMING = N_TIMING, N_FIT_TIMING = N_FIT_TIMING,
             B_TARGETS = B_TARGETS, kappa = KAPPA,
             sessionInfo = sessionInfo()),
        file.path(out_dir, sprintf("timing_raw_%s.rds", stamp)))
cat(sprintf("Done. Output in %s\n", out_dir))
