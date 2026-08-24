# =====================================================================
# 03_simulation_study.R
#
# Design
#   population   latent mediation model of Eq. (2); population values are
#                the ML estimates on the Holzinger-Swineford data
#                (D = 21 free parameters), no effect override
#   functionals  primary   ab, psi_speed = speed~~speed
#                secondary r2_speed, omega_speed
#   factors      N    in {100, 200, 500}
#                dist in {normal, nonnormal}   (VITA/copula, lognormal
#                     margins, skewness ~ 2, excess kurtosis ~ 7.86)
#                spec in {correct, misspec}    (population residual
#                     covariance x4 ~~ x7, calibrated to a population
#                     RMSEA of .04 for the analysis model)
#   constants    S = 1,000 data sets per cell; B = 1,000 multinomial
#                weight vectors per data set, drawn once and shared by
#                the exact bootstrap, IJ1 and HOIJ-2
#   methods      wald_expected, wald_hw, mc_hw, ij1, hoij2, boot
#                (plus hoij2_sens and boot_bca as sensitivity checks;
#                 neither is reported in the article)
#   estimand     phi(theta*). For the normal conditions theta* is
#                obtained by fitting the analysis model to the
#                population covariance matrix; for the VITA conditions
#                from a large simulated sample, because the copula
#                calibration matches the target covariance only
#                numerically
#   outcomes     coverage, left/right tail error, median width, median
#                computation time per data set, failure rates
#
# Cost of a full run: 12 cells x S x (1 + B) = about 12 million lavaan
# fits. Run with SMOKE_TEST = TRUE first.
#
# Run 00_install_dependencies.R once before this script.
# =====================================================================

source("hoij_core.R")
HOIJ_CORE <- normalizePath("hoij_core.R")

for (pkg in c("lavaan", "MASS", "covsim", "rvinecopulib")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' is required; see 00_install_dependencies.R.")
}
suppressPackageStartupMessages({
  library(lavaan); library(parallel); library(covsim); library(rvinecopulib)
})
cat("lavaan", as.character(packageVersion("lavaan")), "\n")
cat("cores detected:", detectCores(), "\n")


# ---------------------------------------------------------------------
# 1. Settings
# ---------------------------------------------------------------------
SMOKE_TEST <- FALSE      # TRUE: pipeline test (S = 10, B = 100)

S <- if (SMOKE_TEST)  10L else 1000L   # data sets per cell
B <- if (SMOKE_TEST) 100L else 1000L   # shared weight vectors per data set
R_MC  <- B                             # Monte Carlo draws, parity with B
ALPHA <- 0.05                          # nominal 95% intervals

SPREAD_TOL  <- 0.1        # tolerance on the gradient-Hessian check
MAX_REDRAW  <- 50L        # cap on redraws per data-set slot
INCLUDE_BCA <- TRUE       # BCa sensitivity check; not reported

MAIN_FNS      <- c("ab", "psi_speed")
SECONDARY_FNS <- c("r2_speed", "omega_speed")
TABLE_FNS     <- c(MAIN_FNS, SECONDARY_FNS)

MISSPEC_PAIR <- c("x4", "x7")   # omitted residual covariance
RMSEA_TARGET <- 0.04            # population RMSEA of the analysis model
N_PSEUDO     <- 1e6             # sample.nobs when fitting to Sigma
N_TRUTH   <- if (SMOKE_TEST)  20000L else  200000L  # VITA truth sample
VITA_NMAX <- if (SMOKE_TEST)  50000L else 1000000L  # VITA calibration size

## Lognormal margins with univariate skewness 2; the implied excess
## kurtosis is about 7.86.
VITA_SDLOG <- uniroot(function(s) {
  v <- exp(s^2); (v + 2) * sqrt(v - 1) - 2
}, interval = c(0.05, 2))$root
VITA_SKEW_TARGET   <- 2
VITA_EXKURT_TARGET <- { v <- exp(VITA_SDLOG^2); v^4 + 2*v^3 + 3*v^2 - 6 }

SEED_BASE <- 2026L

## lavaan fits are CPU bound, so physical cores determine throughput.
ncores <- max(1L, detectCores(logical = FALSE) - 1L)
if (is.na(ncores) || ncores < 1L) ncores <- max(1L, detectCores() - 1L)

out_dir <- file.path(getwd(), "hoij_sim_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

design <- expand.grid(N = c(100L, 200L, 500L),
                      dist = c("normal", "nonnormal"),
                      spec = c("correct", "misspec"),
                      stringsAsFactors = FALSE)
design$cell <- seq_len(nrow(design))

cat(sprintf("Design: %d cells | S = %d | B = %d | workers = %d\n",
            nrow(design), S, B, ncores))
cat(sprintf("Estimated number of lavaan fits: %.1f million\n",
            nrow(design) * S * (1 + B) / 1e6))
if (SMOKE_TEST) cat(">> SMOKE_TEST = TRUE: pipeline test, not the reported run\n")


# ---------------------------------------------------------------------
# 2. Model, population and misspecification
# ---------------------------------------------------------------------
model_analysis <- '
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
  visual ~ c*textual + b*speed
  speed  ~ a*textual
'

fit_hs <- sem(model_analysis, data = HolzingerSwineford1939, estimator = "ML")
stopifnot(lavInspect(fit_hs, "converged"))
theta_pop <- coef(fit_hs, type = "free")
stopifnot(length(theta_pop) == 21L)
cat(sprintf("Population fitted on Holzinger-Swineford data: D = %d\n",
            length(theta_pop)))

## The population is written as a syntax string with fixed values, so
## that simulateData() cannot fall back on starting values.
build_pop_syntax <- function(fit, extra = NULL) {
  pt <- parTable(fit)
  pt <- pt[pt$op %in% c("=~", "~", "~~"), , drop = FALSE]
  paste(c(sprintf("%s %s %.12g*%s", pt$lhs, pt$op, pt$est, pt$rhs), extra),
        collapse = "\n")
}

Sigma_pop <- lavInspect(fit_hs, "implied")$cov[, ]
ov_names  <- rownames(Sigma_pop)
stopifnot(all(MISSPEC_PAIR %in% ov_names))

## Adding a residual covariance delta to Sigma = Lambda Psi Lambda' +
## Theta changes exactly one off-diagonal cell.
make_Sigma_mis <- function(delta) {
  Sm <- Sigma_pop
  Sm[MISSPEC_PAIR[1], MISSPEC_PAIR[2]] <-
    Sm[MISSPEC_PAIR[1], MISSPEC_PAIR[2]] + delta
  Sm[MISSPEC_PAIR[2], MISSPEC_PAIR[1]] <-
    Sm[MISSPEC_PAIR[1], MISSPEC_PAIR[2]]
  Sm
}

## Population RMSEA of the analysis model at a given delta:
## fit to Sigma with a large N, F0 = 2 * fmin, RMSEA = sqrt(F0 / df).
pop_rmsea <- function(delta) {
  Sm <- make_Sigma_mis(delta)
  if (!tryCatch({ chol(Sm); TRUE }, error = function(e) FALSE)) return(10)
  f <- tryCatch(sem(model_analysis, sample.cov = Sm, sample.nobs = N_PSEUDO,
                    estimator = "ML", se = "none"), error = function(e) NULL)
  if (is.null(f) || !lavInspect(f, "converged")) return(10)
  sqrt(max(2 * unname(fitMeasures(f, "fmin")), 0) / unname(fitMeasures(f, "df")))
}

cat(sprintf("\nCalibrating misspecification (%s ~~ %s) to RMSEA = %.3f ...\n",
            MISSPEC_PAIR[1], MISSPEC_PAIR[2], RMSEA_TARGET))
delta_star <- uniroot(function(d) pop_rmsea(d) - RMSEA_TARGET,
                      interval = c(0.005, 0.5), tol = 1e-5)$root
rmsea_achieved <- pop_rmsea(delta_star)
cat(sprintf("  delta* = %.4f  ->  population RMSEA = %.4f\n",
            delta_star, rmsea_achieved))

pop_syntax <- list(
  correct = build_pop_syntax(fit_hs),
  misspec = build_pop_syntax(fit_hs, extra = sprintf(
    "%s ~~ %.12g*%s", MISSPEC_PAIR[1], delta_star, MISSPEC_PAIR[2])))
Sigma_by_spec <- list(correct = Sigma_pop, misspec = make_Sigma_mis(delta_star))


# ---------------------------------------------------------------------
# 3. Non-normal data generation (VITA/covsim)
# ---------------------------------------------------------------------
lnorm_meanlog_for_variance <- function(target_var, sdlog = VITA_SDLOG) {
  ## Var[lognormal] = (exp(sdlog^2) - 1) * exp(2 * meanlog + sdlog^2)
  0.5 * (log(target_var / (exp(sdlog^2) - 1)) - sdlog^2)
}

calibrate_vita <- function(Sigma_target, label, seed_offset = 0L) {
  cat(sprintf("\nCalibrating VITA for %s (Nmax = %d, skew ~ %.2f, excess kurtosis ~ %.2f) ...\n",
              label, VITA_NMAX, VITA_SKEW_TARGET, VITA_EXKURT_TARGET))
  set.seed(SEED_BASE + 7000L + seed_offset)
  margins <- lapply(diag(Sigma_target), function(v)
    list(distr = "lnorm", meanlog = lnorm_meanlog_for_variance(v),
         sdlog = VITA_SDLOG))
  covsim::vita(margins = margins, sigma.target = Sigma_target,
               family_set = c("clayton", "gauss", "joe", "gumbel", "frank"),
               Nmax = VITA_NMAX, verbose = FALSE, cores = 1L)
}

vita_by_spec <- list(correct = calibrate_vita(Sigma_by_spec$correct, "correct", 1L),
                     misspec = calibrate_vita(Sigma_by_spec$misspec, "misspec", 2L))

simulate_vita_data <- function(N, spec, center = TRUE) {
  X <- as.data.frame(rvinecopulib::rvine(N, vita_by_spec[[spec]]))
  names(X) <- ov_names
  ## ML without a mean structure uses the covariance matrix only;
  ## centring stabilises the numerically computed sample covariance.
  if (center) X[] <- lapply(X, function(z) as.numeric(z - mean(z)))
  X
}

moment_diag <- function(dat) data.frame(
  skew   = sapply(dat, function(x) mean((x - mean(x))^3) / sd(x)^3),
  exkurt = sapply(dat, function(x) mean((x - mean(x))^4) / sd(x)^4 - 3))


# ---------------------------------------------------------------------
# 4. Functionals
#
# speed is endogenous, so Var(speed) = a^2 Var(textual) + psi_speed.
# Each functional exists in a scalar and a vectorised form; the latter
# is applied to a matrix of replicates in one call.
# ---------------------------------------------------------------------
fn_names <- TABLE_FNS

var_speed_scalar <- function(th)
  as.numeric(th["a"])^2 * as.numeric(th["textual~~textual"]) +
  as.numeric(th["speed~~speed"])

var_speed_vec <- function(mat, nms)
  mat[, nms == "a"]^2 * mat[, nms == "textual~~textual"] +
  mat[, nms == "speed~~speed"]

functionals_scalar <- list(
  ab        = function(th) as.numeric(th["a"] * th["b"]),
  psi_speed = function(th) as.numeric(th["speed~~speed"]),
  r2_speed  = function(th) {
    var_s <- var_speed_scalar(th)
    if (!is.finite(var_s) || var_s < 1e-12) return(NA_real_)
    1 - as.numeric(th["speed~~speed"]) / var_s
  },
  omega_speed = function(th) {
    ls    <- 1 + as.numeric(th["speed=~x8"]) + as.numeric(th["speed=~x9"])
    num   <- ls^2 * var_speed_scalar(th)
    den   <- num + as.numeric(th["x7~~x7"]) + as.numeric(th["x8~~x8"]) +
      as.numeric(th["x9~~x9"])
    if (!is.finite(den) || abs(den) < 1e-12) return(NA_real_)
    num / den
  })

functionals_vec <- list(
  ab        = function(mat, nms) mat[, nms == "a"] * mat[, nms == "b"],
  psi_speed = function(mat, nms) mat[, nms == "speed~~speed"],
  r2_speed  = function(mat, nms) {
    var_s <- var_speed_vec(mat, nms)
    res <- 1 - mat[, nms == "speed~~speed"] / var_s
    res[!is.finite(var_s) | var_s < 1e-12] <- NA_real_
    res
  },
  omega_speed = function(mat, nms) {
    ls  <- 1 + mat[, nms == "speed=~x8"] + mat[, nms == "speed=~x9"]
    num <- ls^2 * var_speed_vec(mat, nms)
    den <- num + mat[, nms == "x7~~x7"] + mat[, nms == "x8~~x8"] +
      mat[, nms == "x9~~x9"]
    res <- num / den
    res[!is.finite(den) | abs(den) < 1e-12] <- NA_real_
    res
  })

stopifnot(identical(fn_names, c("ab", "psi_speed", "r2_speed", "omega_speed")),
          all(fn_names %in% names(functionals_scalar)),
          all(fn_names %in% names(functionals_vec)))

## The closed forms above are checked against lavaan's own output; an
## error here would contaminate every result, so this is a hard stop.
local({
  th   <- coef(fit_hs, type = "free")
  estm <- lavInspect(fit_hs, "est")
  clv  <- lavInspect(fit_hs, "cov.lv")
  ind  <- paste0("x", 7:9)
  lsum <- sum(estm$lambda[ind, "speed"])
  omega_ref <- lsum^2 * clv["speed", "speed"] /
    (lsum^2 * clv["speed", "speed"] + sum(diag(estm$theta)[ind]))
  stopifnot(
    abs(functionals_scalar$omega_speed(th) - omega_ref) < 1e-10,
    abs(functionals_scalar$psi_speed(th) - th["speed~~speed"]) < 1e-12,
    abs(functionals_scalar$r2_speed(th) -
          (1 - th["speed~~speed"] / clv["speed", "speed"])) < 1e-10)
  cat(sprintf("Functional check OK: psi_speed = %.4f, r2_speed = %.4f, omega_speed = %.4f\n",
              functionals_scalar$psi_speed(th), functionals_scalar$r2_speed(th),
              omega_ref))
})


# ---------------------------------------------------------------------
# 5. Pseudo-true values (coverage targets)
# ---------------------------------------------------------------------
pseudo_truth <- list()
for (sp in names(Sigma_by_spec)) {
  ## normal: exact, by fitting the analysis model to the population Sigma
  f <- sem(model_analysis, sample.cov = Sigma_by_spec[[sp]],
           sample.nobs = N_PSEUDO, estimator = "ML", se = "none")
  stopifnot(lavInspect(f, "converged"))
  th <- coef(f, type = "free")
  pseudo_truth[[paste("normal", sp, sep = "_")]] <-
    sapply(fn_names, function(fn) functionals_scalar[[fn]](th))

  ## non-normal: from a large VITA sample, because the copula matches
  ## the target covariance only numerically
  set.seed(SEED_BASE + 8000L + match(sp, names(Sigma_by_spec)))
  dat_big <- simulate_vita_data(N_TRUTH, sp, center = TRUE)
  f <- sem(model_analysis, sample.cov = cov(dat_big), sample.nobs = N_TRUTH,
           estimator = "ML", se = "none")
  if (!lavInspect(f, "converged"))
    stop("pseudo-truth fit did not converge for nonnormal/", sp)
  th <- coef(f, type = "free")
  pseudo_truth[[paste("nonnormal", sp, sep = "_")]] <-
    sapply(fn_names, function(fn) functionals_scalar[[fn]](th))
}

## reproduction check: for normal/correct, theta* must return phi(theta_pop)
stopifnot(max(abs(pseudo_truth$normal_correct -
                    sapply(fn_names, function(fn)
                      functionals_scalar[[fn]](theta_pop)))) < 1e-6)

cat("\nPseudo-true values (coverage target per dist x spec):\n")
print(round(do.call(rbind, pseudo_truth), 4))


# ---------------------------------------------------------------------
# 6. Pre-flight checks
# ---------------------------------------------------------------------
local({
  set.seed(SEED_BASE)
  md <- moment_diag(simulate_vita_data(50000L, "correct", center = TRUE))
  ok <- abs(mean(md$skew) - VITA_SKEW_TARGET) < 0.35 &&
    abs(mean(md$exkurt) - VITA_EXKURT_TARGET) < 1.75
  cat(sprintf("VITA check: skew = %.2f (target %.2f), excess kurtosis = %.2f (target %.2f)  %s\n",
              mean(md$skew), VITA_SKEW_TARGET, mean(md$exkurt),
              VITA_EXKURT_TARGET, if (ok) "OK" else "FAIL"))
  if (!ok) warning("VITA moments deviate from target", call. = FALSE)
})


# ---------------------------------------------------------------------
# 7. Helpers for the simulation
# ---------------------------------------------------------------------
elapsed_sec <- function(t0) as.numeric(proc.time()[["elapsed"]] - t0)

## BCa interval with IJ1-based acceleration (sensitivity check only)
bca_ci <- function(boot_vals, theta_hat, acc_vals, alpha = ALPHA,
                   min_n = 40L) {
  boot_vals <- boot_vals[is.finite(boot_vals)]
  if (length(boot_vals) < min_n || !is.finite(theta_hat))
    return(c(lo = NA_real_, hi = NA_real_))
  p0 <- mean(boot_vals < theta_hat)
  p0 <- min(max(p0, 1 / (length(boot_vals) + 1)),
            length(boot_vals) / (length(boot_vals) + 1))
  z0 <- qnorm(p0); z_a <- qnorm(1 - alpha / 2)
  acc <- 0
  if (!is.null(acc_vals)) {
    av <- acc_vals[is.finite(acc_vals)]
    if (length(av) >= 10) {
      d <- mean(av) - av
      acc <- sum(d^3) / (6 * (sum(d^2))^1.5)
      if (!is.finite(acc)) acc <- 0
    }
  }
  q <- quantile(boot_vals,
                c(pnorm(z0 + (z0 - z_a) / (1 - acc * (z0 - z_a))),
                  pnorm(z0 + (z0 + z_a) / (1 - acc * (z0 + z_a)))),
                names = FALSE, type = 7)
  c(lo = q[1], hi = q[2])
}

## Bootstrap refits reuse a pre-parsed parameter table with ustart set
## to the ML estimates: lavaan then parses the syntax once per data set
## instead of B times, and every replicate starts warm at theta-hat.
make_boot_partable <- function(fit_k) {
  PT <- parTable(fit_k); PT$ustart <- PT$est; PT$start <- PT$est
  PT$se <- NULL
  PT
}

## Fit to a weighted covariance matrix with all side work switched off.
## check.post = FALSE only disables the check inside the fitter; the
## explicit post.check below still counts inadmissible replicates.
fit_boot_cov <- function(PT_boot, S_b, N_sim, iter.max = 150L) {
  fb <- tryCatch(sem(model = PT_boot, sample.cov = S_b, sample.nobs = N_sim,
                     estimator = "ML", se = "none", test = "none",
                     h1 = FALSE, baseline = FALSE, check.gradient = FALSE,
                     check.start = FALSE, check.post = FALSE,
                     control = list(iter.max = iter.max)),
                 error = function(e) NULL)
  if (!is.null(fb) &&
      tryCatch(isTRUE(lavInspect(fb, "converged")), error = function(e) FALSE))
    return(fb)
  NULL
}

## One-off check that the parameter-table route reproduces a strictly
## converged fit and keeps the coefficient names of the primary fit.
# check_bootstrap_refit_route <- function(tol = 1e-4) {
#   set.seed(SEED_BASE + 909L)
#   dat0 <- simulateData(pop_syntax$correct, sample.nobs = 100L)
#   fit0 <- sem(model_analysis, data = dat0, estimator = "ML",
#               se = "robust.huber.white")
#   stopifnot(lavInspect(fit0, "converged"))
#   S_b <- cov(as.matrix(dat0)[sample.int(nrow(dat0), nrow(dat0), TRUE), ])
# 
#   PT0   <- make_boot_partable(fit0)
#   fb_pt <- fit_boot_cov(PT0, S_b, nrow(dat0))
#   if (is.null(fb_pt)) stop("parameter-table refit route does not converge")
#   th_pt <- coef(fb_pt, type = "free")
#   stopifnot(identical(names(th_pt), names(coef(fit0, type = "free"))))
# 
#   arb <- sem(model = PT0, sample.cov = S_b, sample.nobs = nrow(dat0),
#              estimator = "ML", se = "none", test = "none",
#              control = list(iter.max = 2000L, eval.max = 4000L,
#                             rel.tol = 1e-10))
#   d <- max(abs(th_pt - coef(arb, type = "free")[names(th_pt)]))
#   cat(sprintf("Bootstrap refit route: max |theta(iter.max=150) - theta(strict)| = %.2e (tol %.0e)  %s\n",
#               d, tol, if (d < tol) "OK" else "FAIL"))
#   if (d >= tol) stop("refit route deviates too much from a strict fit")
#   invisible(TRUE)
# }
# check_bootstrap_refit_route()

safe_aggregate <- function(formula, data, FUN, ..., label = deparse(formula)) {
  vars <- all.vars(formula)
  miss <- setdiff(vars, names(data))
  if (length(miss))
    stop(sprintf("aggregation '%s': missing variable(s) %s", label,
                 paste(miss, collapse = ", ")))
  d <- data[complete.cases(data[vars]), , drop = FALSE]
  if (!nrow(d)) stop(sprintf("aggregation '%s': no valid rows", label))
  aggregate(formula, data = d, FUN = FUN, ...)
}


# ---------------------------------------------------------------------
# 8. One data set: six interval methods on shared weight vectors
# ---------------------------------------------------------------------
methods_all <- c("wald_expected", "wald_hw", "mc_hw", "ij1", "hoij2", "boot",
                 "hoij2_sens", if (INCLUDE_BCA) "boot_bca")

run_dataset <- function(cell_row, s) {
  N_sim <- cell_row$N
  set.seed(SEED_BASE + 100000L * cell_row$cell + s)   # seed per task
  syntax <- pop_syntax[[cell_row$spec]]
  truth  <- pseudo_truth[[paste(cell_row$dist, cell_row$spec, sep = "_")]]

  ## --- primary fit; redraw only on non-convergence ------------------
  t0 <- proc.time()[["elapsed"]]
  n_redraw <- 0L; fit_k <- NULL; dat_k <- NULL
  repeat {
    dat_try <- tryCatch(
      if (cell_row$dist == "nonnormal")
        simulate_vita_data(N_sim, cell_row$spec, center = TRUE)
      else simulateData(syntax, sample.nobs = N_sim),
      error = function(e) NULL)
    f_try <- if (is.null(dat_try)) NULL else
      tryCatch(sem(model_analysis, data = dat_try, estimator = "ML",
                   se = "robust.huber.white"), error = function(e) NULL)
    if (tryCatch(!is.null(f_try) && isTRUE(lavInspect(f_try, "converged")),
                 error = function(e) FALSE)) {
      fit_k <- f_try; dat_k <- dat_try; break
    }
    n_redraw <- n_redraw + 1L
    if (n_redraw >= MAX_REDRAW) break
  }
  if (is.null(fit_k))
    return(list(res = NULL,
                diag = data.frame(cell_row, s = s, n_redraw = n_redraw,
                                  converged = FALSE, inadmissible_primary = NA,
                                  boot_fail = NA_real_, boot_inadmiss = NA_real_,
                                  hoij_ok = FALSE, fallb_Hobs = 0L,
                                  fallb_spread = 0L, fallb_deriv = 0L,
                                  grad_spread = NA_real_,
                                  hoij_frac_inadmiss = NA_real_,
                                  error = NA_character_)))
  t_fit <- elapsed_sec(t0)

  ## An inadmissible but converged primary fit is kept: excluding it
  ## would condition coverage on admissibility, while the pseudo-true
  ## value is unconditional, and would selectively remove the small-N
  ## data sets where estimator skewness matters most.
  admiss_primary <- tryCatch(
    isTRUE(suppressWarnings(lavInspect(fit_k, "post.check"))),
    error = function(e) TRUE)

  theta0   <- coef(fit_k, type = "free")
  th_names <- names(theta0)
  D        <- length(theta0)
  spl <- strsplit(th_names, "~~", fixed = TRUE)
  var_idx <- which(vapply(spl, function(z) length(z) == 2 && z[1] == z[2],
                          logical(1)))

  V_hw  <- tryCatch(lavInspect(fit_k, "vcov")[th_names, th_names, drop = FALSE],
                    error = function(e) NULL)
  V_expected <- tryCatch(
    lavTech(fit_k, "inverted.information.expected") / N_sim,
    error = function(e) NULL)

  ## --- Wald intervals ------------------------------------------------
  t0 <- proc.time()[["elapsed"]]
  ci_wald_expected <- lapply(fn_names, function(fn) if (is.null(V_expected))
    c(lo = NA_real_, hi = NA_real_) else
      wald_ci(functionals_scalar[[fn]], theta0, V_expected, alpha = ALPHA))
  t_wexp <- elapsed_sec(t0)
  t0 <- proc.time()[["elapsed"]]
  ci_wald_hw <- lapply(fn_names, function(fn) if (is.null(V_hw))
    c(lo = NA_real_, hi = NA_real_) else
      wald_ci(functionals_scalar[[fn]], theta0, V_hw, alpha = ALPHA))
  t_whw <- elapsed_sec(t0)
  names(ci_wald_expected) <- names(ci_wald_hw) <- fn_names

  ## --- Monte Carlo (HW) ---------------------------------------------
  t0 <- proc.time()[["elapsed"]]
  draws_mc <- NULL
  if (!is.null(V_hw)) {
    L <- tryCatch(t(chol(V_hw + diag(1e-10, D))), error = function(e) NULL)
    if (!is.null(L)) {
      draws_mc <- sweep(matrix(rnorm(R_MC * D), R_MC, D) %*% t(L), 2,
                        theta0, "+")
      colnames(draws_mc) <- th_names
    }
  }
  t_mc <- elapsed_sec(t0)

  ## --- one weight matrix, shared by bootstrap, IJ1 and HOIJ-2 -------
  t0 <- proc.time()[["elapsed"]]
  W_counts <- t(rmultinom(B, size = N_sim, prob = rep(1, N_sim)))   # B x N
  dW <- W_counts - 1L
  t_w <- elapsed_sec(t0)

  ## --- exact bootstrap on the weighted covariance matrix ------------
  ## Without a mean structure the sample covariance matrix is the full
  ## sufficient statistic, so fitting the p x p weighted covariance is
  ## equivalent to fitting the N x p resampled data, but faster.
  t0 <- proc.time()[["elapsed"]]
  boot_th <- matrix(NA_real_, B, D, dimnames = list(NULL, th_names))
  n_boot_fail <- 0L; n_boot_inadmiss <- 0L
  PT_boot <- make_boot_partable(fit_k)
  dat_mat <- as.matrix(dat_k)
  mu_all  <- W_counts %*% dat_mat / N_sim
  for (r in 1:B) {
    w   <- W_counts[r, ]
    S_b <- (crossprod(dat_mat * w, dat_mat) -
              N_sim * tcrossprod(mu_all[r, ])) / (N_sim - 1)
    S_b <- (S_b + t(S_b)) / 2
    fit_b <- fit_boot_cov(PT_boot, S_b, N_sim)
    if (is.null(fit_b)) { n_boot_fail <- n_boot_fail + 1L; next }
    adm <- tryCatch(isTRUE(suppressWarnings(lavInspect(fit_b, "post.check"))),
                    error = function(e) TRUE)
    ## Non-converged or inadmissible refits are dropped from the
    ## bootstrap percentiles, but their weight vectors still feed IJ1
    ## and HOIJ-2, which cannot fail in this way.
    if (!adm) { n_boot_inadmiss <- n_boot_inadmiss + 1L; next }
    th_b <- tryCatch(coef(fit_b, type = "free"), error = function(e) NULL)
    if (!is.null(th_b) && length(th_b) == D) boot_th[r, ] <- th_b
  }
  t_boot <- elapsed_sec(t0)

  ## --- IJ1, Eq. (7) --------------------------------------------------
  t0 <- proc.time()[["elapsed"]]
  ij_th <- NULL; C_mat <- NULL
  Scores <- tryCatch(lavScores(fit_k, scaling = TRUE), error = function(e) NULL)
  H.inv  <- tryCatch(lavTech(fit_k, "inverted.information.observed"),
                     error = function(e) NULL)
  if (!is.null(Scores) && !is.null(H.inv)) {
    ij1 <- ij1_replicates(theta0, Scores, H.inv, dW)
    ij_th <- ij1$theta; C_mat <- ij1$C
  }
  t_scores <- elapsed_sec(t0)

  ## --- HOIJ-2, Eq. (8): derivative setup and replication loop -------
  hoij_th <- NULL; hoij_inadmiss <- NULL
  fallb_Hobs <- 0L; fallb_spread <- 0L; fallb_deriv <- 0L
  grad_spread <- NA_real_
  t_hsetup <- 0; t_hloop <- 0
  if (!is.null(C_mat)) {
    t0 <- proc.time()[["elapsed"]]
    H_obs <- tryCatch(lavTech(fit_k, "information.observed"),
                      error = function(e) NULL)
    if (is.null(H_obs)) {
      fallb_Hobs <- 1L
    } else {
      grad_F <- make_grad_F(fit_k)
      chk <- tryCatch(check_gradient_hessian(grad_F, theta0, H_obs),
                      error = function(e) NULL)
      if (!is.null(chk) && is.finite(chk$spread)) grad_spread <- chk$spread
      ## No silent fallback: if the second-order step cannot be computed
      ## reliably, hoij2 is recorded as NA and the reason is counted.
      if (is.null(chk) || !is.finite(chk$spread) ||
          chk$spread > SPREAD_TOL) {
        fallb_spread <- 1L
      } else {
        T_arr <- tryCatch(compute_T_tensor_grad(grad_F, theta0),
                          error = function(e) NULL)
        H_all <- tryCatch(compute_all_H(fit_k, theta0), error = function(e) NULL)
        if (is.null(T_arr) || is.null(H_all)) {
          fallb_deriv <- 1L
        } else {
          t_hsetup <- elapsed_sec(t0)
          t0 <- proc.time()[["elapsed"]]
          hoij_th <- hoij2_replicates(theta0, C_mat, dW, H.inv, H_all, T_arr)
          hoij_inadmiss <- apply(hoij_th[, var_idx, drop = FALSE] < 0, 1, any)
          t_hloop <- elapsed_sec(t0)
        }
      }
    }
    if (is.null(hoij_th)) t_hsetup <- elapsed_sec(t0)
  }

  ## IJ acceleration for the optional BCa interval (no delete-1 refits)
  acc_th <- if (INCLUDE_BCA && !is.null(Scores) && !is.null(H.inv))
    sweep(Scores %*% H.inv, 2, theta0, "+") else NULL
  if (!is.null(acc_th)) colnames(acc_th) <- th_names

  ## --- intervals and bookkeeping per functional x method ------------
  min_n <- max(40L, ceiling(0.5 * B))
  out <- vector("list", length(fn_names) * length(methods_all)); oi <- 0L
  for (fn in fn_names) {
    f_v <- functionals_vec[[fn]]
    fn_hat <- tryCatch(functionals_scalar[[fn]](theta0),
                       error = function(e) NA_real_)
    tv <- truth[[fn]]

    mc_vals   <- if (!is.null(draws_mc)) f_v(draws_mc, th_names) else NULL
    boot_vals <- f_v(boot_th, th_names)
    ij_vals   <- if (!is.null(ij_th))   f_v(ij_th, th_names)   else NULL
    hoij_vals <- if (!is.null(hoij_th)) f_v(hoij_th, th_names) else NULL

    cis <- list(
      wald_expected = ci_wald_expected[[fn]][c("lo", "hi")],
      wald_hw       = ci_wald_hw[[fn]][c("lo", "hi")],
      mc_hw    = if (is.null(mc_vals)) c(NA_real_, NA_real_) else
        percentile_ci(mc_vals, ALPHA, max(40L, ceiling(0.5 * R_MC))),
      ij1      = if (is.null(ij_vals)) c(NA_real_, NA_real_) else
        percentile_ci(ij_vals, ALPHA, min_n),
      ## primary convention: inadmissible approximate values are kept
      hoij2    = if (is.null(hoij_vals)) c(NA_real_, NA_real_) else
        percentile_ci(hoij_vals, ALPHA, min_n),
      boot     = percentile_ci(boot_vals, ALPHA, min_n),
      ## sensitivity check: drop replicates with a negative variance
      hoij2_sens = if (is.null(hoij_vals)) c(NA_real_, NA_real_) else
        percentile_ci(hoij_vals[!hoij_inadmiss], ALPHA, min_n))
    if (INCLUDE_BCA)
      cis$boot_bca <- bca_ci(boot_vals, fn_hat,
                             if (is.null(acc_th)) NULL else
                               f_v(acc_th, th_names), min_n = min_n)

    ## HOIJ-2 timing includes the primary fit, the derivative setup and
    ## all B weight evaluations.
    tms <- c(wald_expected = t_fit + t_wexp, wald_hw = t_fit + t_whw,
             mc_hw = t_fit + t_mc, ij1 = t_fit + t_w + t_scores,
             hoij2 = t_fit + t_w + t_scores + t_hsetup + t_hloop,
             boot = t_fit + t_w + t_boot,
             hoij2_sens = t_fit + t_w + t_scores + t_hsetup + t_hloop,
             boot_bca = t_fit + t_w + t_boot)

    boot_ok <- is.finite(boot_th[, 1])
    ffin <- c(
      wald_expected = NA_real_, wald_hw = NA_real_,
      mc_hw = if (is.null(mc_vals)) NA_real_ else mean(is.finite(mc_vals)),
      ij1   = if (is.null(ij_vals)) NA_real_ else mean(is.finite(ij_vals)),
      hoij2 = if (is.null(hoij_vals)) NA_real_ else mean(is.finite(hoij_vals)),
      boot  = mean(boot_ok),
      hoij2_sens = if (is.null(hoij_vals)) NA_real_ else
        mean(is.finite(hoij_vals) & !hoij_inadmiss),
      boot_bca = mean(boot_ok))

    for (m in methods_all) {
      lo <- unname(cis[[m]][1]); hi <- unname(cis[[m]][2])
      oi <- oi + 1L
      out[[oi]] <- data.frame(
        cell = cell_row$cell, N = N_sim, dist = cell_row$dist,
        spec = cell_row$spec, s = s, functional = fn, method = m,
        lo = lo, hi = hi, est = fn_hat, truth = tv,
        covered = if (is.finite(lo) && is.finite(hi))
          as.integer(tv >= lo & tv <= hi) else NA_integer_,
        miss_left  = if (is.finite(lo)) as.integer(tv < lo) else NA_integer_,
        miss_right = if (is.finite(hi)) as.integer(tv > hi) else NA_integer_,
        width = if (is.finite(lo) && is.finite(hi)) hi - lo else NA_real_,
        time_s = as.numeric(tms[[m]]), frac_finite = as.numeric(ffin[[m]]),
        stringsAsFactors = FALSE)
    }
  }

  list(res = do.call(rbind, out),
       diag = data.frame(cell_row, s = s, n_redraw = n_redraw,
                         converged = TRUE,
                         inadmissible_primary = !admiss_primary,
                         boot_fail = n_boot_fail / B,
                         boot_inadmiss = n_boot_inadmiss / B,
                         hoij_ok = !is.null(hoij_th), fallb_Hobs = fallb_Hobs,
                         fallb_spread = fallb_spread, fallb_deriv = fallb_deriv,
                         grad_spread = grad_spread,
                         hoij_frac_inadmiss = if (is.null(hoij_inadmiss))
                           NA_real_ else mean(hoij_inadmiss),
                         error = NA_character_))
}


# ---------------------------------------------------------------------
# 9. Parallel execution
#
# The seed is set per (cell, data set) task inside run_dataset(), so the
# result is exactly reproducible and independent of the scheduling.
# ---------------------------------------------------------------------
tasks <- do.call(rbind, lapply(seq_len(nrow(design)), function(ci)
  data.frame(cell = ci, s = seq_len(S))))
cat(sprintf("\n%d tasks (cell x data set)\n", nrow(tasks)))

## A task that dies takes the whole parLapplyLB() call with it, so an
## unexpected error is caught, recorded in 'diags$error' and skipped;
## it is never silent, because the run prints the offending tasks below.
run_task <- function(ti) {
  tk <- tasks[ti, ]
  cell_row <- design[design$cell == tk$cell, ]
  tryCatch(run_dataset(cell_row, tk$s), error = function(e)
    list(res = NULL,
         diag = data.frame(cell_row, s = tk$s, n_redraw = NA_integer_,
                           converged = NA, inadmissible_primary = NA,
                           boot_fail = NA_real_, boot_inadmiss = NA_real_,
                           hoij_ok = NA, fallb_Hobs = NA_integer_,
                           fallb_spread = NA_integer_, fallb_deriv = NA_integer_,
                           grad_spread = NA_real_,
                           hoij_frac_inadmiss = NA_real_,
                           error = conditionMessage(e),
                           stringsAsFactors = FALSE)))
}

cl <- makeCluster(ncores, type = "PSOCK")
clusterExport(cl, "HOIJ_CORE")
invisible(clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(lavaan); library(covsim); library(rvinecopulib)
  })
  source(HOIJ_CORE)
  ## guard against BLAS oversubscription (workers x BLAS threads)
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1)
  }
  NULL
}))
clusterExport(cl, c("design", "tasks", "run_dataset", "pop_syntax",
                    "Sigma_by_spec", "vita_by_spec", "ov_names",
                    "pseudo_truth", "model_analysis", "fn_names",
                    "functionals_scalar", "functionals_vec", "methods_all",
                    "S", "B", "R_MC", "ALPHA", "SPREAD_TOL",
                    "MAX_REDRAW", "INCLUDE_BCA",
                    "SEED_BASE", "make_boot_partable", "fit_boot_cov",
                    "simulate_vita_data", "elapsed_sec", "moment_diag",
                    "var_speed_scalar", "var_speed_vec", "bca_ci"))

## Load-balanced chunks; task order is shuffled deterministically so
## that slow (N = 500) and fast (N = 100) tasks are mixed within chunks.
set.seed(SEED_BASE)
task_order <- sample(seq_len(nrow(tasks)))
CHUNK  <- max(ncores * 20L, 48L)
chunks <- split(task_order, ceiling(seq_along(task_order) / CHUNK))
res_list <- vector("list", length(chunks))
t_start <- proc.time()[3]
for (ch in seq_along(chunks)) {
  res_list[[ch]] <- parLapplyLB(cl, chunks[[ch]], run_task)
  el <- proc.time()[3] - t_start
  done <- sum(lengths(chunks[seq_len(ch)]))
  cat(sprintf("  chunk %d/%d | tasks %d/%d | %.1f min elapsed | ETA %.1f min\n",
              ch, length(chunks), done, nrow(tasks), el / 60,
              el / 60 * (nrow(tasks) / done - 1)))
  ## checkpoint per chunk; to resume, read the partial files back in
  saveRDS(res_list[[ch]],
          #file.path(out_dir, sprintf("hoij_sim_N100_partial_chunk%03d.rds", ch)))
          file.path(out_dir, sprintf("hoij_sim_partial_chunk%03d.rds", ch)))
}
stopCluster(cl)

flat    <- unlist(res_list, recursive = FALSE)
results <- do.call(rbind, lapply(flat, `[[`, "res"))
diags   <- do.call(rbind, lapply(flat, `[[`, "diag"))
if (is.null(results) || nrow(results) == 0L)
  stop("no valid result rows were produced; inspect 'diags'")

n_err <- sum(!is.na(diags$error))
if (n_err) {
  cat(sprintf("\n%d of %d tasks failed with an error:\n", n_err, nrow(tasks)))
  print(head(unique(diags$error[!is.na(diags$error)]), 10L))
} else cat("\nNo task failed with an error.\n")

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
#write.csv(results, file.path(out_dir, sprintf("hoij_sim_results_S%d_%s.csv",
#                                              S, stamp)), row.names = FALSE)
saveRDS(list(results = results, diags = diags, design = design,
             pseudo_truth = pseudo_truth, delta_star = delta_star,
             rmsea_achieved = rmsea_achieved,
             config = list(S = S, B = B, R_MC = R_MC, ALPHA = ALPHA,
                           INCLUDE_BCA = INCLUDE_BCA,
                           N_TRUTH = N_TRUTH, VITA_NMAX = VITA_NMAX,
                           VITA_SDLOG = VITA_SDLOG, SEED_BASE = SEED_BASE,
                           SMOKE_TEST = SMOKE_TEST),
             sessionInfo = sessionInfo()),
        #file.path(out_dir, sprintf("hoij_sim_N100_S%d_%s.rds", S, stamp)))
        file.path(out_dir, sprintf("hoij_sim_full_S%d_%s.rds", S, stamp)))
cat(sprintf("\nRaw results written (%d rows)\n", nrow(results)))


# ---------------------------------------------------------------------
# 10. Aggregation per cell x functional x method
# ---------------------------------------------------------------------
summ <- Reduce(function(a, b)
  merge(a, b, by = c("N", "dist", "spec", "functional", "method")),
  list(
    safe_aggregate(cbind(covered, miss_left, miss_right) ~
                     N + dist + spec + functional + method, results,
                   function(x) mean(x, na.rm = TRUE), label = "coverage"),
    safe_aggregate(width ~ N + dist + spec + functional + method, results,
                   median, na.rm = TRUE, label = "median width"),
    safe_aggregate(time_s ~ N + dist + spec + functional + method, results,
                   median, na.rm = TRUE, label = "median time"),
    setNames(safe_aggregate(covered ~ N + dist + spec + functional + method,
                            results, function(x) sum(is.finite(x)),
                            label = "n valid"),
             c("N", "dist", "spec", "functional", "method", "n_valid"))))

summ$mcse_cov <- sqrt(0.95 * 0.05 / pmax(summ$n_valid, 1))
mcse_ref <- sqrt(0.95 * 0.05 / S)
cat(sprintf("\nMonte Carlo error on coverage at full validity (S = %d): +/- %.2f percentage points\n",
            S, 100 * mcse_ref))

fail_tab <- aggregate(cbind(n_redraw, inadmissible_primary, boot_fail,
                            boot_inadmiss, hoij_frac_inadmiss) ~ N + dist + spec,
                      data = diags, FUN = function(x) mean(x, na.rm = TRUE))
cat("\nFailure and HOIJ diagnostics per cell (means):\n")
print(fail_tab, row.names = FALSE, digits = 3)


# ---------------------------------------------------------------------
# 11. Table bodies and Figure 2
# ---------------------------------------------------------------------
method_order <- c("wald_expected", "wald_hw", "mc_hw", "ij1", "hoij2", "boot")
method_label <- c(wald_expected = "Wald (Expected)", wald_hw = "Wald (HW)",
                  mc_hw = "Monte Carlo (HW)", ij1 = "IJ1 percentile",
                  hoij2 = "HOIJ-2 percentile", boot = "Bootstrap percentile")
functional_label <- list(ab = expression(italic(ab)),
                         psi_speed = expression(psi[speed]),
                         r2_speed = expression(R^2 * "" [speed]),
                         omega_speed = expression(omega[speed]))
cond_order <- data.frame(
  dist = c("normal", "normal", "nonnormal", "nonnormal"),
  spec = c("correct", "misspec", "correct", "misspec"),
  lab  = c("Normal data, correct model", "Normal data, misspecified model",
           "Non-normal data, correct model",
           "Non-normal data, misspecified model"))

fmt <- function(x, d = 3) ifelse(is.finite(x),
                                 formatC(x, digits = d, format = "f"), "--")

emit_sim_table <- function(fn, N_val, file) {
  con <- file(file, "w"); on.exit(close(con))
  for (ci in seq_len(nrow(cond_order))) {
    writeLines(sprintf("\\multicolumn{5}{l}{\\emph{%s}} \\\\",
                       cond_order$lab[ci]), con)
    for (m in method_order) {
      r <- summ[summ$functional == fn & summ$N == N_val &
                  summ$dist == cond_order$dist[ci] &
                  summ$spec == cond_order$spec[ci] & summ$method == m, ]
      if (!nrow(r)) r <- data.frame(covered = NA, miss_left = NA,
                                    miss_right = NA, width = NA)
      writeLines(sprintf("\\quad %-22s & %s & %s & %s & %s \\\\",
                         method_label[m], fmt(100 * r$covered[1], 1),
                         fmt(100 * r$miss_left[1], 1),
                         fmt(100 * r$miss_right[1], 1), fmt(r$width[1], 3)),
                 con)
    }
    if (ci < nrow(cond_order)) writeLines("\\midrule", con)
  }
  cat("  written:", file, "\n")
}

for (N_val in unique(design$N))
  for (fn in TABLE_FNS)
    emit_sim_table(fn, N_val,
                   file.path(out_dir, sprintf("tab_sim_%s_N%d.tex", fn, N_val)))

## median time per data set and failure rates, pooled over cells
tt <- safe_aggregate(time_s ~ method, results, median, na.rm = TRUE,
                     label = "time per method")
fail_by_m <- c(
  wald_expected = 0, wald_hw = 0,
  mc_hw = 100 * (1 - mean(results$frac_finite[results$method == "mc_hw"],
                          na.rm = TRUE)),
  ij1   = 100 * (1 - mean(results$frac_finite[results$method == "ij1"],
                          na.rm = TRUE)),
  hoij2 = 100 * (1 - mean(results$frac_finite[results$method == "hoij2"],
                          na.rm = TRUE)),
  boot  = 100 * (mean(diags$boot_fail, na.rm = TRUE) +
                   mean(diags$boot_inadmiss, na.rm = TRUE)))
t_boot_med <- tt$time_s[tt$method == "boot"]
con <- file(file.path(out_dir, "tab_sim_time.tex"), "w")
for (m in method_order)
  writeLines(sprintf("%-22s & %s & %s & %s \\\\", method_label[m],
                     fmt(tt$time_s[tt$method == m], 2), fmt(fail_by_m[m], 2),
                     fmt(tt$time_s[tt$method == m] / t_boot_med, 3)), con)
close(con)
cat("  written:", file.path(out_dir, "tab_sim_time.tex"), "\n")

pdf(file.path(out_dir, "fig-sim-coverage.pdf"), width = 11, height = 6.8)
n_col <- length(unique(design$N))
layout(rbind(matrix(seq_len(length(MAIN_FNS) * n_col), nrow = length(MAIN_FNS),
                    byrow = TRUE),
             rep(length(MAIN_FNS) * n_col + 1, n_col)),
       heights = c(rep(1, length(MAIN_FNS)), 0.24))
op <- par(mar = c(4, 10.5, 2.5, 1), mgp = c(2.2, 0.7, 0))

pch_map <- c(wald_expected = 1, wald_hw = 2, mc_hw = 15, ij1 = 16, hoij2 = 17,
             boot = 18)
col_map <- setNames(c("black", "grey45", "firebrick", "steelblue"),
                    paste(cond_order$dist, cond_order$spec))

for (fn in MAIN_FNS) {
  for (N_val in sort(unique(design$N))) {
    plot(NA, xlim = c(0.80, 1.00), ylim = c(0.5, length(method_order) + 0.5),
         yaxt = "n", xlab = "Empirical coverage", ylab = "",
         main = bquote(.(functional_label[[fn]]) * "," ~ N == .(N_val)))
    rect(0.95 - 2 * mcse_ref, 0, 0.95 + 2 * mcse_ref, length(method_order) + 1,
         col = "grey90", border = NA)
    abline(v = 0.95, lty = 2, lwd = 1.2)
    axis(2, at = seq_along(method_order), las = 1, cex.axis = 1.0,
         labels = gsub("--", "-", method_label[method_order]))
    for (ci in seq_len(nrow(cond_order))) {
      key <- paste(cond_order$dist[ci], cond_order$spec[ci])
      for (mi in seq_along(method_order)) {
        r <- summ[summ$functional == fn & summ$N == N_val &
                    summ$dist == cond_order$dist[ci] &
                    summ$spec == cond_order$spec[ci] &
                    summ$method == method_order[mi], ]
        if (nrow(r) && is.finite(r$covered[1])) {
          y_i <- mi + (ci - 2.5) * 0.13
          if (is.finite(r$mcse_cov[1]))
            segments(r$covered[1] - 2 * r$mcse_cov[1], y_i,
                     r$covered[1] + 2 * r$mcse_cov[1], y_i, lwd = 2,
                     col = adjustcolor(col_map[key], alpha.f = 0.45))
          points(r$covered[1], y_i, pch = pch_map[method_order[mi]],
                 col = col_map[key], cex = 1.1)
        }
      }
    }
  }
}
par(mar = c(0, 0, 0, 0)); plot.new()
legend("top", legend = cond_order$lab, col = col_map, pch = 15, bty = "n",
       cex = 1.05, pt.cex = 1.5, x.intersp = 0.8, horiz = TRUE,
       title = "Condition", title.font = 2)
legend("bottom", fill = "grey90", border = NA, bty = "n", cex = 1.0,
       legend = expression(paste("± 2 Monte Carlo SEs around ", .95)))
par(op); layout(1)
dev.off()
cat("  written:", file.path(out_dir, "fig-sim-coverage.pdf"), "\n")

cat(sprintf("\nDone. %s\n", if (SMOKE_TEST)
  "This was a SMOKE_TEST; set SMOKE_TEST <- FALSE for the reported run." else
    sprintf("Full run completed (S = %d, B = %d).", S, B)))



# =====================================================================
# 04_regenerate_fig_sim_coverage.R
#
# =====================================================================

out_dir <- "hoij_sim_output"   # <-- pas aan naar jouw pad indien nodig

# ---------------------------------------------------------------------
# 1. Chunkbestanden herladen en samenvoegen
# ---------------------------------------------------------------------
chunk_files <- sort(list.files(out_dir,
                               pattern = "^hoij_sim_partial_chunk[0-9]+\\.rds$",
                               full.names = TRUE))
if (!length(chunk_files))
  stop("Geen chunkbestanden gevonden in '", out_dir, "'")
cat(sprintf("Gevonden chunkbestanden: %d\n", length(chunk_files)))

flat <- unlist(lapply(chunk_files, readRDS), recursive = FALSE)
cat(sprintf("Totaal aantal samengevoegde taken (cel x dataset): %d\n", length(flat)))

results <- do.call(rbind, lapply(flat, `[[`, "res"))
diags   <- do.call(rbind, lapply(flat, `[[`, "diag"))

if (is.null(results) || !nrow(results))
  stop("Geen geldige resultaatrijen gevonden; controleer 'diags'.")

n_err <- sum(!is.na(diags$error))
if (n_err) {
  cat(sprintf("%d taken faalden met een fout (eerste unieke meldingen):\n", n_err))
  print(head(unique(diags$error[!is.na(diags$error)]), 10L))
} else cat("Geen taken met een fout.\n")

cat("\nAantal taken per cel (diags$cell):\n")
print(table(diags$cell))

# ---------------------------------------------------------------------
# 2. Design en constantes reconstrueren
# ---------------------------------------------------------------------
design <- expand.grid(N = c(100L, 200L, 500L),
                      dist = c("normal", "nonnormal"),
                      spec = c("correct", "misspec"),
                      stringsAsFactors = FALSE)
design$cell <- seq_len(nrow(design))

S <- length(unique(results$s))
mcse_ref <- sqrt(0.95 * 0.05 / S)
cat(sprintf("\nAfgeleide S (datasets per cel) = %d  ->  MCSE-referentie = %.4f\n",
            S, mcse_ref))

safe_aggregate <- function(formula, data, FUN, ..., label = deparse(formula)) {
  vars <- all.vars(formula)
  miss <- setdiff(vars, names(data))
  if (length(miss))
    stop(sprintf("aggregation '%s': missing variable(s) %s", label,
                 paste(miss, collapse = ", ")))
  d <- data[complete.cases(data[vars]), , drop = FALSE]
  if (!nrow(d)) stop(sprintf("aggregation '%s': no valid rows", label))
  aggregate(formula, data = d, FUN = FUN, ...)
}

# ---------------------------------------------------------------------
# 3. Aggregatie per cel x functional x methode (identiek aan sectie 10
#    van 03_simulation_study.R)
# ---------------------------------------------------------------------
summ <- Reduce(function(a, b)
  merge(a, b, by = c("N", "dist", "spec", "functional", "method")),
  list(
    safe_aggregate(cbind(covered, miss_left, miss_right) ~
                     N + dist + spec + functional + method, results,
                   function(x) mean(x, na.rm = TRUE), label = "coverage"),
    safe_aggregate(width ~ N + dist + spec + functional + method, results,
                   median, na.rm = TRUE, label = "median width"),
    safe_aggregate(time_s ~ N + dist + spec + functional + method, results,
                   median, na.rm = TRUE, label = "median time"),
    setNames(safe_aggregate(covered ~ N + dist + spec + functional + method,
                            results, function(x) sum(is.finite(x)),
                            label = "n valid"),
             c("N", "dist", "spec", "functional", "method", "n_valid"))))

summ$mcse_cov <- sqrt(0.95 * 0.05 / pmax(summ$n_valid, 1))

# ---------------------------------------------------------------------
# 4. Figuur: 3 rijen (ab, psi_speed, omega_speed) x 3 N-kolommen
# ---------------------------------------------------------------------
MAIN_FNS <- c("ab", "psi_speed", "omega_speed")   # <-- derde rij toegevoegd
method_order <- c("wald_expected", "wald_hw", "mc_hw", "ij1", "hoij2", "boot")
method_label <- c(wald_expected = "Wald (Expected)", wald_hw = "Wald (HW)",
                  mc_hw = "Monte Carlo (HW)", ij1 = "IJ1 percentile",
                  hoij2 = "HOIJ-2 percentile", boot = "Bootstrap percentile")
functional_label <- list(ab = quote(italic(ab)),
                         psi_speed = quote(psi[speed]),
                         r2_speed = quote(R^2 * "" [speed]),
                         omega_speed = quote(omega[speed]))
cond_order <- data.frame(
  dist = c("normal", "normal", "nonnormal", "nonnormal"),
  spec = c("correct", "misspec", "correct", "misspec"),
  lab  = c("Normal data, correct model", "Normal data, misspecified model",
           "Non-normal data, correct model",
           "Non-normal data, misspecified model"))

fig_path <- file.path(out_dir, "fig-sim-coverage.pdf")
pdf(fig_path, width = 11, height = 9.7)   # hoogte omhoog voor de 3e rij
n_col <- length(unique(design$N))
layout(rbind(matrix(seq_len(length(MAIN_FNS) * n_col), nrow = length(MAIN_FNS),
                    byrow = TRUE),
             rep(length(MAIN_FNS) * n_col + 1, n_col)),
       heights = c(rep(1, length(MAIN_FNS)), 0.2))
op <- par(mar = c(4, 10.5, 2.5, 1), mgp = c(2.2, 0.7, 0))

pch_map <- c(wald_expected = 1, wald_hw = 2, mc_hw = 15, ij1 = 16, hoij2 = 17,
             boot = 18)
col_map <- setNames(c("black", "grey45", "firebrick", "steelblue"),
                    paste(cond_order$dist, cond_order$spec))

for (fn in MAIN_FNS) {
  for (N_val in sort(unique(design$N))) {
    plot(NA, xlim = c(0.80, 1.00), ylim = c(0.5, length(method_order) + 0.5),
         yaxt = "n", xlab = "Empirical coverage", ylab = "",
         main = bquote(.(functional_label[[fn]]) * "," ~ N == .(N_val)))
    rect(0.95 - 2 * mcse_ref, 0, 0.95 + 2 * mcse_ref, length(method_order) + 1,
         col = "grey90", border = NA)
    abline(v = 0.95, lty = 2, lwd = 1.2)
    axis(2, at = seq_along(method_order), las = 1, cex.axis = 1.0,
         labels = gsub("--", "-", method_label[method_order]))
    for (ci in seq_len(nrow(cond_order))) {
      key <- paste(cond_order$dist[ci], cond_order$spec[ci])
      for (mi in seq_along(method_order)) {
        r <- summ[summ$functional == fn & summ$N == N_val &
                    summ$dist == cond_order$dist[ci] &
                    summ$spec == cond_order$spec[ci] &
                    summ$method == method_order[mi], ]
        if (nrow(r) && is.finite(r$covered[1])) {
          y_i <- mi + (ci - 2.5) * 0.13
          if (is.finite(r$mcse_cov[1]))
            segments(r$covered[1] - 2 * r$mcse_cov[1], y_i,
                     r$covered[1] + 2 * r$mcse_cov[1], y_i, lwd = 2,
                     col = adjustcolor(col_map[key], alpha.f = 0.45))
          points(r$covered[1], y_i, pch = pch_map[method_order[mi]],
                 col = col_map[key], cex = 1.1)
        }
      }
    }
  }
}
par(mar = c(0, 0, 0, 0)); plot.new()
legend("top", legend = cond_order$lab, col = col_map, pch = 15, bty = "n",
       cex = 1.05, pt.cex = 1.5, x.intersp = 0.8, horiz = TRUE,
       title = "Condition", title.font = 2)
legend("bottom", fill = "grey90", border = NA, bty = "n", cex = 1.0,
       legend = expression(paste("± 2 Monte Carlo SEs around ", .95)))
par(op); layout(1)
dev.off()
cat("  written:", fig_path, "\n")