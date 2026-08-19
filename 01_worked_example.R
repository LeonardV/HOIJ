# =====================================================================
# 01_worked_example.R
#
# Companion code for:
#   Vanbrabant, L., & Rosseel, Y. Approximating percentile bootstrap
#   confidence intervals in SEM without repeated refitting: A tutorial
#   on the second-order infinitesimal jackknife.
#
# Reproduces Section 3 (tutorial example):
#   Table 2  distributional SD, skewness and 95% limits for ab,
#            psi_speed and omega_speed
#   Table 5  the same summaries for |ab|, theta_11, R2_speed and P_M
#   Figure 1 fig-shape.pdf, replicate densities for the three targets
#
# One data set (N = 500) from a weak-effect latent mediation population,
# and B = 5,000 multinomial weight vectors that are shared by IJ1,
# HOIJ-2 and the exact bootstrap, so that any difference between the
# approximate and the exact distribution is approximation error rather
# than a different set of bootstrap draws.
#
# Run 00_install_dependencies.R once before this script.
# =====================================================================

source("hoij_core.R")
suppressPackageStartupMessages(library(lavaan))

## --- settings --------------------------------------------------------
SEED_DATA    <- 20260706
SEED_WEIGHTS <- 20260707
SEED_MC      <- 20260708

N_EX     <- 500      # sample size of the single example data set
B        <- 5000     # weight vectors shared by IJ1, HOIJ-2 and bootstrap
R_MC     <- 5000     # Monte Carlo draws
KAPPA    <- 0.5      # trust-region damping of the second-order step
ALPHA_CI <- 0.05     # nominal 95% intervals

## Weak-effect population: only the two mediation paths are overridden,
## all remaining population values are the Holzinger-Swineford estimates.
## These two values must match the values reported in Section 3.1.
EFFECT_PARS <- c(a = 0.10, b = 0.10)

out_dir <- "hoij_worked_example_output"
if (!dir.exists(out_dir)) dir.create(out_dir)

cat("lavaan", as.character(packageVersion("lavaan")), "\n")
hoij_selftest()


# ---------------------------------------------------------------------
# 1. Model, population and data
# ---------------------------------------------------------------------
model_med <- '
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9

  speed   ~ a*textual
  visual  ~ b*speed + c*textual
'

fit_hs <- sem(model_med, data = HolzingerSwineford1939, se = "none")
stopifnot(lavInspect(fit_hs, "converged"))

pt_pop <- parTable(fit_hs)
for (par in names(EFFECT_PARS)) {
  pt_pop$ustart[pt_pop$label == par] <- EFFECT_PARS[par]
  pt_pop$est[pt_pop$label == par]    <- EFFECT_PARS[par]
}

set.seed(SEED_DATA)
dat <- simulateData(pt_pop, sample.nobs = N_EX)

fit <- sem(model_med, data = dat, se = "none", estimator = "ML")
stopifnot(lavInspect(fit, "converged"))

theta0   <- coef(fit, type = "free")
th_names <- names(theta0)
D        <- length(theta0)
N        <- nrow(dat)
stopifnot(D == 21, N == N_EX)

cat(sprintf("Analysis model: D = %d, N = %d\n", D, N))
cat(sprintf("a = %.4f, b = %.4f, ab = %.4f, psi_speed = %.4f\n",
            theta0["a"], theta0["b"], theta0["a"] * theta0["b"],
            theta0["speed~~speed"]))


# ---------------------------------------------------------------------
# 2. Functionals
#
# All functionals are written for a matrix of parameter vectors (one
# replicate per row), so that they can be applied to the Monte Carlo,
# IJ1, HOIJ-2 and bootstrap replicates in one call. The functional is
# always evaluated AFTER the replicate has been constructed.
# ---------------------------------------------------------------------
par_cols <- function(mat, nms) {
  idx <- match(nms, colnames(mat))
  if (anyNA(idx)) stop("free parameter(s) not found: ",
                       paste(nms[is.na(idx)], collapse = ", "))
  mat[, idx, drop = FALSE]
}

## model-implied latent variances
var_textual <- function(mat) as.numeric(par_cols(mat, "textual~~textual"))
var_speed <- function(mat) {
  p <- par_cols(mat, c("a", "textual~~textual", "speed~~speed"))
  p[, 1]^2 * p[, 2] + p[, 3]
}
var_visual <- function(mat) {
  p <- par_cols(mat, c("a", "b", "c", "textual~~textual",
                       "visual~~visual"))
  p[, 2]^2 * var_speed(mat) + p[, 3]^2 * p[, 4] +
    2 * p[, 2] * p[, 3] * p[, 1] * p[, 4] + p[, 5]
}

## McDonald's omega for a factor with a unit-loading marker item
omega_factor <- function(mat, factor, items, var_fun) {
  lambda_sum <- 1 + rowSums(par_cols(mat, paste0(factor, "=~", items[-1])))
  theta_sum  <- rowSums(par_cols(mat, paste0(items, "~~", items)))
  num   <- lambda_sum^2 * var_fun(mat)
  denom <- num + theta_sum
  out <- num / denom
  out[!is.finite(denom) | abs(denom) < 1e-10] <- NA_real_
  as.numeric(out)
}

functionals_all <- list(
  ## Table 2: main teaching targets
  ab          = function(mat) par_cols(mat, "a")[, 1] * par_cols(mat, "b")[, 1],
  psi_speed   = function(mat) as.numeric(par_cols(mat, "speed~~speed")),
  omega_speed = function(mat) omega_factor(mat, "speed", paste0("x", 7:9),
                                           var_speed),
  ## Table 5: supporting functionals
  abs_ab      = function(mat) abs(par_cols(mat, "a")[, 1] *
                                    par_cols(mat, "b")[, 1]),
  theta11     = function(mat) as.numeric(par_cols(mat, "x1~~x1")),
  r2_speed    = function(mat) {
    psi <- as.numeric(par_cols(mat, "speed~~speed"))
    vs  <- var_speed(mat)
    out <- 1 - psi / vs
    out[!is.finite(vs) | vs < 1e-10] <- NA_real_
    out
  },
  pm          = function(mat) {
    ab_val <- par_cols(mat, "a")[, 1] * par_cols(mat, "b")[, 1]
    denom  <- ab_val + par_cols(mat, "c")[, 1]
    out <- ab_val / denom
    out[!is.finite(denom) | abs(denom) < 1e-6] <- NA_real_
    out
  },
  ## reported as a diagnostic only (Section 3 uses the speed factor)
  omega_vis   = function(mat) omega_factor(mat, "visual", paste0("x", 1:3),
                                           var_visual),
  omega_text  = function(mat) omega_factor(mat, "textual", paste0("x", 4:6),
                                           var_textual)
)

FN_MAIN      <- c("ab", "psi_speed", "omega_speed")          # Table 2
FN_SECONDARY <- c("abs_ab", "theta11", "r2_speed", "pm")     # Table 5
fn_names     <- c(FN_MAIN, FN_SECONDARY)

## a single parameter vector as a one-row matrix
as_row <- function(th) matrix(th, nrow = 1,
                              dimnames = list(NULL, names(th)))
fn_hat <- vapply(functionals_all, function(f) f(as_row(theta0)), numeric(1))

## the two identity-type functionals must reproduce lavaan exactly
stopifnot(
  abs(fn_hat["psi_speed"] - theta0["speed~~speed"]) < 1e-12,
  abs(fn_hat["r2_speed"] -
        (1 - theta0["speed~~speed"] /
           lavInspect(fit, "cov.lv")["speed", "speed"])) < 1e-10,
  abs(fn_hat["omega_speed"] - {
    est <- lavInspect(fit, "est")
    ls  <- sum(est$lambda[paste0("x", 7:9), "speed"])
    num <- ls^2 * lavInspect(fit, "cov.lv")["speed", "speed"]
    num / (num + sum(diag(est$theta)[paste0("x", 7:9)]))
  }) < 1e-10
)

cat("\nPoint estimates:\n"); print(round(fn_hat, 4))


# ---------------------------------------------------------------------
# 3. Covariance matrices for the two Wald comparators
# ---------------------------------------------------------------------
fit_inf <- sem(model_med, data = dat, estimator = "ML",
               se = "standard", information = "expected")
fit_hw  <- sem(model_med, data = dat, estimator = "ML",
               se = "robust.huber.white")
stopifnot(lavInspect(fit_inf, "converged"), lavInspect(fit_hw, "converged"))

V_inf <- lavInspect(fit_inf, "vcov")[th_names, th_names, drop = FALSE]
V_hw  <- lavInspect(fit_hw,  "vcov")[th_names, th_names, drop = FALSE]
stopifnot(max(abs(coef(fit_hw, type = "free")[th_names] - theta0)) < 1e-6)


# ---------------------------------------------------------------------
# 4. One-time HOIJ setup: scores, curvature and third derivatives
# ---------------------------------------------------------------------
Scores <- lavScores(fit, scaling = TRUE)                  # N x D
H.inv  <- lavTech(fit, "inverted.information.observed")   # Jhat^-1
H_obs  <- lavTech(fit, "information.observed")            # Jhat
dimnames(H.inv) <- list(th_names, th_names)

grad_F <- make_grad_F(fit)
chk <- check_gradient_hessian(grad_F, theta0, H_obs)
if (!is.finite(chk$spread) || chk$spread > 0.1)
  stop(sprintf("derivative check failed (spread = %.3g)", chk$spread))

J_all <- compute_all_J(fit, theta0)                # N x D x D
T_arr <- compute_T_tensor_grad(grad_F, theta0)     # D x D x D, equals -Khat


# ---------------------------------------------------------------------
# 5. Shared weights, exact bootstrap, IJ1, HOIJ-2 and Monte Carlo
# ---------------------------------------------------------------------
set.seed(SEED_WEIGHTS)
W_counts <- t(rmultinom(B, size = N, prob = rep(1 / N, N)))   # B x N
dW <- W_counts - 1L

cat(sprintf("\nExact bootstrap: %d refits ...\n", B))
boot_th <- matrix(NA_real_, B, D, dimnames = list(NULL, th_names))
for (bb in seq_len(B)) {
  dat_b <- dat[rep.int(seq_len(N), W_counts[bb, ]), , drop = FALSE]
  fit_b <- tryCatch(sem(model_med, data = dat_b, se = "none",
                        estimator = "ML", start = fit),
                    error = function(e) NULL)
  if (!is.null(fit_b) && lavInspect(fit_b, "converged"))
    boot_th[bb, ] <- coef(fit_b, type = "free")
  if (bb %% 1000 == 0) cat(sprintf("  %d / %d\n", bb, B))
}
valid <- complete.cases(boot_th)
cat(sprintf("  %d converged, %d failed (%.2f%%)\n",
            sum(valid), sum(!valid), 100 * mean(!valid)))

ij1   <- ij1_replicates(theta0, Scores, H.inv, dW)               # Eq. (7)
hoij2 <- hoij2_replicates(theta0, ij1$C, dW, H.inv, J_all, T_arr,
                          kappa = KAPPA)                          # Eq. (8)
cat(sprintf("HOIJ-2 damping: fraction s < 1 = %.3f, mean s = %.3f\n",
            mean(hoij2$s < 1), mean(hoij2$s)))

set.seed(SEED_MC)
L_mc  <- t(chol(V_hw + diag(1e-10, D)))
mc_th <- sweep(matrix(rnorm(R_MC * D), R_MC, D) %*% t(L_mc), 2, theta0, "+")
colnames(mc_th) <- th_names

## The approximate replicates are restricted to the weight vectors on
## which the exact bootstrap converged, so all four distributions are
## based on the same weights.
replicates <- list(mc_hw = mc_th,
                   ij1   = ij1$theta[valid, , drop = FALSE],
                   hoij2 = hoij2$theta[valid, , drop = FALSE],
                   boot  = boot_th[valid, , drop = FALSE])


# ---------------------------------------------------------------------
# 6. Tables 2 and 5
# ---------------------------------------------------------------------
method_label <- c(wald_inf = "Wald--delta (Inf)",
                  wald_hw  = "Wald--delta (HW)",
                  mc_hw    = "Monte Carlo (HW)",
                  ij1      = "IJ1 percentile",
                  hoij2    = "HOIJ-2 percentile",
                  boot     = "Bootstrap percentile")

summarise_functional <- function(fn) {
  f_mat <- functionals_all[[fn]]
  f_sca <- function(th) f_mat(as_row(th))
  est   <- fn_hat[fn]

  wald_row <- function(method, V) {
    ci <- wald_ci(f_sca, theta0, V, alpha = ALPHA_CI)
    data.frame(functional = fn, method = method, SE = unname(ci["se"]),
               skewness = 0, lo = unname(ci["lo"]), hi = unname(ci["hi"]))
  }
  perc_row <- function(method) {
    v  <- f_mat(replicates[[method]])
    vf <- v[is.finite(v)]
    ci <- percentile_ci(vf, alpha = ALPHA_CI)
    data.frame(functional = fn, method = method,
               SE = if (length(vf) > 1) sd(vf) else NA_real_,
               skewness = skewness(vf),
               lo = unname(ci["lo"]), hi = unname(ci["hi"]))
  }
  rbind(wald_row("wald_inf", V_inf), wald_row("wald_hw", V_hw),
        do.call(rbind, lapply(names(replicates), perc_row)))
}

tab_intervals <- do.call(rbind, lapply(c(fn_names, "omega_vis", "omega_text"),
                                       summarise_functional))
tab_intervals$method_label <- method_label[tab_intervals$method]

cat("\n===== interval summaries =====\n")
print(tab_intervals[, c("functional", "method_label", "SE", "skewness",
                        "lo", "hi")], row.names = FALSE, digits = 4)

## LaTeX table bodies, in the column order of the manuscript
functional_tex <- c(ab = "$ab$", psi_speed = "$\\psi_{\\mathrm{speed}}$",
                    omega_speed = "$\\omega_{\\mathrm{speed}}$",
                    abs_ab = "$|ab|$", theta11 = "$\\theta_{11}$",
                    r2_speed = "$R^2_{\\mathrm{speed}}$", pm = "$P_M$")
fmt <- function(x, d = 3) if (is.na(x)) "--" else
  formatC(x, digits = d, format = "f")

emit_tex <- function(fns, file) {
  con <- file(file, "w"); on.exit(close(con))
  for (fn in fns) {
    sub <- tab_intervals[tab_intervals$functional == fn, ]
    first <- TRUE
    for (m in names(method_label)) {
      r <- sub[sub$method == m, ]
      writeLines(sprintf("%s & %s & %s & %s & %s & %s \\\\",
                         if (first) functional_tex[fn] else "",
                         method_label[m], fmt(r$SE), fmt(r$skewness),
                         fmt(r$lo), fmt(r$hi)), con)
      first <- FALSE
    }
    if (fn != tail(fns, 1)) writeLines("\\midrule", con)
  }
  cat("  written:", file, "\n")
}
emit_tex(FN_MAIN,      file.path(out_dir, "tab_intervals.tex"))
emit_tex(FN_SECONDARY, file.path(out_dir, "tab_intervals_secondary.tex"))


# ---------------------------------------------------------------------
# 7. Figure 1
# ---------------------------------------------------------------------
plot_panel <- function(fn, xlab) {
  f <- functionals_all[[fn]]
  dens <- lapply(replicates[c("boot", "hoij2", "ij1", "mc_hw")], function(m) {
    v <- f(m); density(v[is.finite(v)])
  })
  plot(NA, bty = "l", ylab = "Density", xlab = xlab, main = "",
       xlim = range(vapply(dens, function(d) range(d$x), numeric(2))),
       ylim = c(0, max(vapply(dens, function(d) max(d$y), numeric(1)))))
  polygon(dens$boot$x, dens$boot$y,
          col = adjustcolor("grey60", alpha.f = 0.5), border = NA)
  lines(dens$hoij2$x, dens$hoij2$y, lwd = 2, lty = 1)
  lines(dens$ij1$x,   dens$ij1$y,   lwd = 2, lty = 2)
  lines(dens$mc_hw$x, dens$mc_hw$y, lwd = 2, lty = 3)
  abline(v = fn_hat[fn], lwd = 1)
}

pdf(file.path(out_dir, "fig-shape.pdf"), width = 13, height = 4.2)
par(mfrow = c(1, 3), mar = c(4, 4, 1, 1))
plot_panel("ab", expression(italic(ab)))
legend("topright", bty = "n", cex = 0.8, merge = FALSE,
       legend = c("Bootstrap", "HOIJ-2", "IJ1", "Monte Carlo"),
       fill = c(adjustcolor("grey60", alpha.f = 0.5), NA, NA, NA),
       border = c("grey60", NA, NA, NA),
       lty = c(NA, 1, 2, 3), lwd = c(NA, 2, 2, 2))
## psi_speed is an identity functional, so Monte Carlo coincides with
## Wald--delta (HW) by construction: any difference from the bootstrap
## in this panel comes from the reweighted estimator itself.
plot_panel("psi_speed", expression(psi[speed]))
plot_panel("omega_speed", expression(omega[speed]))
dev.off()
cat("  written:", file.path(out_dir, "fig-shape.pdf"), "\n")


# ---------------------------------------------------------------------
# 8. Replicate-level agreement between HOIJ-2 and the exact bootstrap
# ---------------------------------------------------------------------
tab_pairwise <- do.call(rbind, lapply(fn_names, function(fn) {
  f <- functionals_all[[fn]]
  v_boot  <- f(replicates$boot)
  v_hoij2 <- f(replicates$hoij2)
  both <- is.finite(v_boot) & is.finite(v_hoij2)
  data.frame(functional = fn, estimate = fn_hat[fn],
             boot_skew = skewness(v_boot[both]),
             hoij2_skew = skewness(v_hoij2[both]),
             cor_boot_hoij2 = cor(v_boot[both], v_hoij2[both]),
             n_pair = sum(both))
}))
cat("\n===== HOIJ-2 vs exact bootstrap, same weight vectors =====\n")
print(tab_pairwise, row.names = FALSE, digits = 4)


# ---------------------------------------------------------------------
# 9. Save
# ---------------------------------------------------------------------
stamp <- format(Sys.time(), "%Y%m%d_%H%M")
write.csv(tab_intervals,
          file.path(out_dir, sprintf("tab_intervals_%s.csv", stamp)),
          row.names = FALSE)
saveRDS(list(seeds = c(data = SEED_DATA, weights = SEED_WEIGHTS,
                       mc = SEED_MC),
             effect_pars = EFFECT_PARS, N = N, B = B, R_MC = R_MC, D = D,
             kappa = KAPPA, theta0 = theta0, fn_hat = fn_hat,
             valid = valid, replicates = replicates, W_counts = W_counts,
             damping_s = hoij2$s, V_inf = V_inf, V_hw = V_hw,
             tab_intervals = tab_intervals, tab_pairwise = tab_pairwise,
             sessionInfo = sessionInfo()),
        file.path(out_dir, sprintf("worked_example_%s.rds", stamp)))
cat(sprintf("\nDone. Output in %s\n", out_dir))
