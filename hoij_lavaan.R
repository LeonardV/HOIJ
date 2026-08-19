# =====================================================================
# hoij_lavaan(): HOIJ-2 standard errors and confidence intervals for a
# fitted lavaan model
#
# Companion code for:
#   Vanbrabant, L., & Rosseel, Y. Approximating percentile bootstrap
#   confidence intervals in SEM without repeated refitting: A tutorial
#   on the second-order infinitesimal jackknife.
#
# This is the reusable version of the machinery used in the analysis
# scripts: it approximates the bootstrap distribution of (functions of)
# the free parameters from a single fitted model, without refitting.
#
# For each multinomial weight vector w (as a bootstrap resample would
# produce), with C = (w - 1)' S Jhat^-1 (see hoij_core.R),
#   IJ1    : theta(w) ~ theta-hat - C                         Eq. (7)
#   HOIJ-2 : theta(w) ~ theta-hat - C + s * d2                Eq. (8)
#            d2 = Jhat^-1 J_delta C - 1/2 Jhat^-1 T(C, C)
#            s  = min(1, kappa * ||C|| / ||d2||)
# The functional is evaluated on each replicate; the standard error is
# the standard deviation and the interval the percentile interval of
# those values.
#
# Scope (enforced by .hoij_check_fit):
#   single group, complete data, estimator ML, continuous indicators,
#   no equality constraints, no multilevel or sampling weights
#
# Usage:
#   source("hoij_core.R"); source("hoij_lavaan.R")
#   fit <- sem(model, data = dat)
#   hoij_lavaan(fit, functional = c("a*b", "speed~~speed"))
# =====================================================================

if (!exists("compute_all_J", mode = "function")) source("hoij_core.R")


# ---------------------------------------------------------------------
# Scope checks
# ---------------------------------------------------------------------
.hoij_check_fit <- function(fit) {
  if (!inherits(fit, "lavaan"))
    stop("'fit' must be a fitted lavaan object.", call. = FALSE)
  if (!isTRUE(lavaan::lavInspect(fit, "converged")))
    stop("The lavaan model did not converge.", call. = FALSE)
  if (lavaan::lavInspect(fit, "ngroups") != 1L)
    stop("hoij_lavaan() currently supports single-group models only.",
         call. = FALSE)
  if (lavaan::lavInspect(fit, "nlevels") > 1L)
    stop("hoij_lavaan() does not support multilevel models.", call. = FALSE)
  opt <- lavaan::lavInspect(fit, "options")
  if (!opt$estimator %in% "ML")
    stop("hoij_lavaan() requires estimator = \"ML\" (found: ",
         opt$estimator, ").", call. = FALSE)
  if (opt$missing %in% c("ml", "fiml", "ml.x", "two.stage",
                         "robust.two.stage"))
    stop("hoij_lavaan() requires complete data (missing = \"", opt$missing,
         "\" is not supported).", call. = FALSE)
  if (isTRUE(lavaan::lavInspect(fit, "categorical")))
    stop("hoij_lavaan() does not support categorical indicators.",
         call. = FALSE)
  if (fit@Model@eq.constraints ||
      (!is.null(fit@Model@ceq.function) &&
       !identical(body(fit@Model@ceq.function), quote(NULL)) &&
       length(fit@Model@ceq.linear.idx) +
       length(fit@Model@ceq.nonlinear.idx) > 0))
    stop("hoij_lavaan() does not support equality constraints; defined ",
         "parameters can be passed through 'functional' instead.",
         call. = FALSE)
  invisible(TRUE)
}


# ---------------------------------------------------------------------
# Functional interface
#   NULL       -> every free parameter
#   character  -> an expression in parameter names, e.g. "a*b" or
#                 "1 - `speed~~speed` / `visual~~visual`"
#                 (names that are not syntactic R names need backticks)
#   function   -> phi(theta), with theta the named free-parameter vector
#   a (named) list of expressions and/or functions is also allowed
# ---------------------------------------------------------------------
.hoij_make_functionals <- function(functional, th_names) {
  as_fun <- function(x, label) {
    if (is.function(x)) return(x)
    if (is.character(x) && length(x) == 1L) {
      expr <- parse(text = x)[[1]]
      return(function(theta) eval(expr, envir = as.list(theta)))
    }
    stop("Functional '", label, "' must be a function or a character ",
         "expression.", call. = FALSE)
  }
  if (is.null(functional)) {
    fns <- lapply(th_names, function(nm) { force(nm)
      function(theta) unname(theta[nm]) })
    names(fns) <- th_names
    return(fns)
  }
  if (is.function(functional)) functional <- list(functional)
  if (is.character(functional)) functional <- as.list(functional)
  if (!is.list(functional))
    stop("'functional' must be NULL, a character vector, a function, or a ",
         "list of these.", call. = FALSE)

  nms <- names(functional)
  if (is.null(nms)) nms <- rep("", length(functional))
  labels <- vapply(seq_along(functional), function(i)
    if (nzchar(nms[i])) nms[i] else
      if (is.character(functional[[i]])) functional[[i]] else
        sprintf("phi%d", i), character(1))

  fns <- lapply(seq_along(functional), function(i)
    as_fun(functional[[i]], labels[i]))
  names(fns) <- labels
  fns
}


# ---------------------------------------------------------------------
# Main function
#
# @param fit         converged lavaan object (single group, ML, complete data)
# @param functional  NULL, character expression(s), function(s), or a list
# @param B           number of multinomial weight vectors
# @param order       1 = IJ1 (linear), 2 = HOIJ-2 (default)
# @param kappa       trust-region damping of the second-order step: its
#                    norm is capped at kappa times the norm of the
#                    first-order step. The bound is active for a
#                    sizeable minority of weight vectors, so report
#                    diagnostics$frac_damped; kappa = Inf gives Eq. (8)
#                    unmodified
# @param level       confidence level of the percentile interval
# @param admissibility "keep" (default: all replicates count) or "drop"
#                    (replicates with a negative variance parameter are
#                    removed before the SE and interval are computed)
# @param alpha_spread_tol tolerance on the scale calibration
# @param seed        optional seed for the weight draws
# @param details     TRUE: also return the replicates and the weight matrix
# @return object of class "hoij_lavaan" with $results (est, se, lo, hi),
#         $diagnostics and optionally $replicates / $weights
# ---------------------------------------------------------------------
hoij_lavaan <- function(fit, functional = NULL, B = 1000L, order = 2L,
                        kappa = 0.5, level = 0.95,
                        admissibility = c("keep", "drop"),
                        alpha_spread_tol = 0.1, seed = NULL,
                        details = FALSE) {

  admissibility <- match.arg(admissibility)
  stopifnot(order %in% c(1L, 2L), B >= 40L, level > 0, level < 1, kappa > 0)
  .hoij_check_fit(fit)
  if (!is.null(seed)) set.seed(seed)

  theta0   <- lavaan::coef(fit, type = "free")
  th_names <- names(theta0)
  D        <- length(theta0)
  N        <- nrow(fit@Data@X[[1]])
  fns      <- .hoij_make_functionals(functional, th_names)

  ## Guard against a silently broken link to lavaan's internals: a
  ## perturbation of theta must change the casewise log-likelihood.
  th_pert <- theta0; th_pert[1] <- th_pert[1] + 1e-3
  if (identical(compute_loglik_casewise(fit, theta0),
                compute_loglik_casewise(fit, th_pert)))
    stop("Broken link to the lavaan internals: perturbing theta does not ",
         "change the casewise log-likelihood. See 00_install_dependencies.R.",
         call. = FALSE)

  ## variance parameters, used for the admissibility of replicates
  spl <- strsplit(th_names, "~~", fixed = TRUE)
  var_idx <- which(vapply(spl, function(z) length(z) == 2 && z[1] == z[2],
                          logical(1)))

  ## --- one-time setup ------------------------------------------------
  t0 <- proc.time()[["elapsed"]]
  Scores <- lavaan::lavScores(fit, scaling = TRUE)              # N x D
  H.inv  <- lavaan::lavTech(fit, "inverted.information.observed")
  if (is.null(Scores) || is.null(H.inv))
    stop("Scores or observed information are not available for this model.",
         call. = FALSE)

  alpha <- NA_real_; alpha_spread <- NA_real_; J_all <- NULL; T_arr <- NULL
  if (order == 2L) {
    H_obs  <- lavaan::lavTech(fit, "information.observed")
    grad_F <- make_grad_F(fit)
    cal <- calibrate_alpha(grad_F, theta0, H_obs)
    alpha <- cal$alpha; alpha_spread <- cal$spread
    if (!is.finite(alpha) || cal$spread > alpha_spread_tol)
      stop(sprintf(paste0("Scale calibration is inconsistent (spread = %.3g ",
                          "> %.3g): the second-order step cannot be computed ",
                          "reliably for this model. Use order = 1 (IJ1) or a ",
                          "bootstrap."), cal$spread, alpha_spread_tol),
           call. = FALSE)
    T_arr <- compute_T_tensor_grad(grad_F, theta0, alpha)
    J_all <- compute_all_J(fit, theta0)
  }
  t_setup <- proc.time()[["elapsed"]] - t0

  ## --- weights and replicates ----------------------------------------
  t0 <- proc.time()[["elapsed"]]
  W  <- t(rmultinom(B, size = N, prob = rep(1, N)))             # B x N
  dW <- W - 1L

  ij1 <- ij1_replicates(theta0, Scores, H.inv, dW)
  if (order == 2L) {
    h2 <- hoij2_replicates(theta0, ij1$C, dW, H.inv, J_all, T_arr,
                           kappa = kappa)
    theta_rep <- h2$theta; s_vec <- h2$s
  } else {
    theta_rep <- ij1$theta; s_vec <- rep(1, B)
  }
  inadmiss <- if (length(var_idx))
    apply(theta_rep[, var_idx, drop = FALSE] < 0, 1, any) else rep(FALSE, B)
  t_rep <- proc.time()[["elapsed"]] - t0

  ## --- SE and percentile interval per functional ----------------------
  keep <- if (admissibility == "drop") !inadmiss else rep(TRUE, B)
  if (sum(keep) < max(40L, ceiling(0.5 * B)))
    warning("Fewer than half of the replicates are admissible; the ",
            "standard errors and intervals may be unreliable.", call. = FALSE)
  pr <- c((1 - level) / 2, 1 - (1 - level) / 2)

  res <- do.call(rbind, lapply(names(fns), function(nm) {
    phi  <- fns[[nm]]
    est  <- tryCatch(as.numeric(phi(theta0)), error = function(e) NA_real_)
    vals <- vapply(seq_len(B), function(i)
      tryCatch(as.numeric(phi(theta_rep[i, ])), error = function(e) NA_real_),
      numeric(1))
    v <- vals[keep & is.finite(vals)]
    if (length(v) >= 40L) {
      q <- quantile(v, pr, names = FALSE, type = 7); se <- sd(v)
    } else {
      q <- c(NA_real_, NA_real_); se <- NA_real_
    }
    data.frame(functional = nm, est = est, se = se, lo = q[1], hi = q[2],
               n_used = length(v), stringsAsFactors = FALSE)
  }))
  rownames(res) <- NULL

  out <- list(
    results = res,
    diagnostics = list(
      order = order, B = B, kappa = kappa, level = level,
      admissibility = admissibility, alpha = alpha,
      alpha_spread = alpha_spread,
      frac_damped = if (order == 2L) mean(s_vec < 1) else NA_real_,
      mean_s = if (order == 2L) mean(s_vec) else NA_real_,
      frac_inadmissible = mean(inadmiss),
      time_setup_s = t_setup, time_replicates_s = t_rep, N = N, D = D),
    call = match.call())
  if (details) { out$replicates <- theta_rep; out$weights <- W }
  class(out) <- "hoij_lavaan"
  out
}


# ---------------------------------------------------------------------
# Print method
# ---------------------------------------------------------------------
print.hoij_lavaan <- function(x, digits = 3, ...) {
  d <- x$diagnostics
  cat(if (d$order == 2L) "HOIJ-2 (second-order infinitesimal jackknife)\n"
      else "IJ1 (first-order infinitesimal jackknife)\n")
  cat(sprintf("B = %d weight vectors | %d%% percentile CI | N = %d, D = %d\n",
              d$B, round(100 * d$level), d$N, d$D))
  cat(sprintf("setup %.2fs + replicates %.2fs\n",
              d$time_setup_s, d$time_replicates_s))
  if (d$order == 2L)
    cat(sprintf("alpha = %.2f (spread %.1e) | damped %.1f%% (mean s = %.2f)\n",
                d$alpha, d$alpha_spread, 100 * d$frac_damped, d$mean_s))
  if (d$frac_inadmissible > 0)
    cat(sprintf("inadmissible replicates: %.1f%% (%s)\n",
                100 * d$frac_inadmissible,
                if (d$admissibility == "keep") "kept" else "dropped"))
  cat("\n")
  r <- x$results
  r[c("est", "se", "lo", "hi")] <- lapply(r[c("est", "se", "lo", "hi")],
                                          round, digits = digits)
  print(r, row.names = FALSE)
  invisible(x)
}
