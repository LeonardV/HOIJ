# =====================================================================
# hoij_core.R -- computational kernel 
#
# Companion code for:
#   Vanbrabant, L., & Rosseel, Y. Approximating percentile bootstrap
#   confidence intervals in SEM without repeated refitting: A tutorial
#   on the second-order infinitesimal jackknife.
#
# Notation follows the article:
#   s_i(theta)   casewise score                                 Eq. (2)
#   H_i(theta)   casewise observed information                  Eq. (3)
#   Hhat         mean casewise information at theta-hat
#   g_delta      weight perturbation of the score               Eq. (6)
#   H_delta      weight perturbation of the information         Eq. (9)
#   Khat(u, v)   third-derivative contraction                   Eq. (12)
#   IJ1          theta-hat + Hhat^-1 g_delta                    Eq. (7)
#   HOIJ-2       IJ1 - Hhat^-1 H_delta d + 1/2 Hhat^-1 Khat(d, d), Eq. (8)
#
# lavaan conventions relied upon:
#   lavScores(fit, scaling = TRUE)         = -s_i(theta-hat) / N
#   lavTech(fit, "information.observed")   = sum_i H_i(theta-hat) / N = Hhat
#   compute_T_tensor_grad()                = -Khat
# =====================================================================

# lavaan internals

.hoij_internals <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    find_fun <- function(nms) {
      for (nm in nms) {
        f <- tryCatch(get(nm, envir = asNamespace("lavaan")),
                      error = function(e) NULL)
        if (is.function(f)) return(f)
      }
      stop("lavaan internal not found (tried: ", paste(nms, collapse = ", "),
           "); this lavaan version (",
           as.character(utils::packageVersion("lavaan")),
           ") is not supported.", call. = FALSE)
    }
    glist_arg <- function(f) {
      fa <- names(formals(f))
      if ("glist" %in% fa) "glist" else if ("GLIST" %in% fa) "GLIST" else
        stop("lav_model_implied() no longer has a glist/GLIST argument; ",
             "hoij_core.R needs to be updated.", call. = FALSE)
    }
    implied  <- find_fun("lav_model_implied")
    gradient <- find_fun(c("lav_model_gradient", "lav_model_grad"))
    cache <<- list(
      x2glist       = find_fun(c("lav_model_x2glist", "lav_model_x2GLIST")),
      implied       = implied,
      gradient      = gradient,
      implied_glist = glist_arg(implied),
      grad_glist    = glist_arg(gradient))
    cache
  }
})



# Casewise log-likelihood, Eq. (1)
compute_loglik_casewise <- function(fit, theta) {
  X <- fit@Data@X[[1]]
  N <- nrow(X); p <- ncol(X)

  ints  <- .hoij_internals()
  GLIST <- ints$x2glist(fit@Model, x = theta)
  args  <- list(fit@Model); args[[ints$implied_glist]] <- GLIST
  implied <- do.call(ints$implied, args)

  Sigma <- implied$cov[[1]]
  mu    <- implied$mean[[1]]
  if (is.null(mu) || length(mu) == 0) mu <- colMeans(X)

  Sigma_inv <- tryCatch(solve(Sigma), error = function(e) MASS::ginv(Sigma))
  log_det   <- determinant(Sigma, logarithm = TRUE)$modulus[1]

  X_centered <- sweep(X, 2, mu, "-")
  quad_form  <- rowSums((X_centered %*% Sigma_inv) * X_centered)

  as.numeric(-0.5 * (p * log(2 * pi) + log_det) - 0.5 * quad_form)
}


# Casewise observed information H_i, Eq. (3)
compute_all_H <- function(fit, theta0, delta = 1e-5) {
  D <- length(theta0); N <- nrow(fit@Data@X[[1]])
  H_array <- array(0, dim = c(N, D, D))
  ll_0 <- compute_loglik_casewise(fit, theta0)

  bump <- function(idx, sgn) {
    th <- theta0; th[idx] <- th[idx] + sgn * delta; th
  }

  for (k in seq_len(D)) {
    for (l in k:D) {
      if (k == l) {
        ll_p <- compute_loglik_casewise(fit, bump(k,  1))
        ll_m <- compute_loglik_casewise(fit, bump(k, -1))
        H_array[, k, k] <- -(ll_p - 2 * ll_0 + ll_m) / delta^2
      } else {
        th_pp <- bump(k, 1);  th_pp[l] <- th_pp[l] + delta
        th_pm <- bump(k, 1);  th_pm[l] <- th_pm[l] - delta
        th_mp <- bump(k, -1); th_mp[l] <- th_mp[l] + delta
        th_mm <- bump(k, -1); th_mm[l] <- th_mm[l] - delta
        H_array[, k, l] <- -(compute_loglik_casewise(fit, th_pp) -
                             compute_loglik_casewise(fit, th_pm) -
                             compute_loglik_casewise(fit, th_mp) +
                             compute_loglik_casewise(fit, th_mm)) / (4 * delta^2)
        H_array[, l, k] <- H_array[, k, l]
      }
    }
  }
  H_array
}


# Analytic gradient of lavaan's fit function F, as a function of theta
make_grad_F <- function(fit) {
  ints <- .hoij_internals()
  lavmodel <- fit@Model; lavsamplestats <- fit@SampleStats
  lavdata  <- fit@Data;  lavcache       <- fit@Cache

  function(theta) {
    args <- list(lavmodel = lavmodel, lavsamplestats = lavsamplestats,
                 lavdata = lavdata, lavcache = lavcache)
    args[[ints$grad_glist]] <- ints$x2glist(lavmodel, x = theta)
    as.numeric(do.call(ints$gradient, args))
  }
}


# Consistency check on the finite-difference route
check_gradient_hessian <- function(grad_F, theta0, H_observed, h = 1e-5) {
  D <- length(theta0)
  H_grad <- matrix(NA_real_, D, D)
  for (k in seq_len(D)) {
    tp <- theta0; tp[k] <- tp[k] + h
    tm <- theta0; tm[k] <- tm[k] - h
    H_grad[, k] <- (grad_F(tp) - grad_F(tm)) / (2 * h)
  }
  H_grad <- (H_grad + t(H_grad)) / 2

  idx <- abs(H_grad) > 1e-6 * max(abs(H_grad))
  ratio <- as.numeric(H_observed)[idx] / as.numeric(H_grad)[idx]

  list(spread = max(abs(ratio / median(ratio) - 1)),
       ratio = median(ratio))
}


# Third-derivative array, Eq. (12)
compute_T_tensor_grad <- function(grad_F, theta, h = 1e-4) {
  D <- length(theta)
  T_arr <- array(0, dim = c(D, D, D))
  g0 <- grad_F(theta)

  for (l in seq_len(D)) {
    for (m in l:D) {
      if (l == m) {
        tp <- theta; tp[l] <- tp[l] + h
        tm <- theta; tm[l] <- tm[l] - h
        col <- (grad_F(tp) - 2 * g0 + grad_F(tm)) / h^2
      } else {
        t_pp <- theta; t_pp[l] <- t_pp[l] + h; t_pp[m] <- t_pp[m] + h
        t_pm <- theta; t_pm[l] <- t_pm[l] + h; t_pm[m] <- t_pm[m] - h
        t_mp <- theta; t_mp[l] <- t_mp[l] - h; t_mp[m] <- t_mp[m] + h
        t_mm <- theta; t_mm[l] <- t_mm[l] - h; t_mm[m] <- t_mm[m] - h
        col <- (grad_F(t_pp) - grad_F(t_pm) -
                grad_F(t_mp) + grad_F(t_mm)) / (4 * h^2)
      }
      T_arr[, l, m] <- col
      T_arr[, m, l] <- col
    }
  }

  (T_arr + aperm(T_arr, c(2, 1, 3)) + aperm(T_arr, c(3, 2, 1)) +
    aperm(T_arr, c(1, 3, 2)) + aperm(T_arr, c(2, 3, 1)) +
    aperm(T_arr, c(3, 1, 2))) / 6
}


# First-order replicates, Eq. (7)
ij1_replicates <- function(theta0, Scores, H.inv, dW) {
  C_mat <- (dW %*% Scores) %*% H.inv                       # B x D
  theta_rep <- sweep(-C_mat, 2, theta0, "+")
  colnames(theta_rep) <- names(theta0)
  list(theta = theta_rep, C = C_mat)
}


# Second-order replicates, Eq. (8)
hoij2_replicates <- function(theta0, C_mat, dW, H.inv, H_all, T_arr) {
  D <- length(theta0); B <- nrow(C_mat); N <- dim(H_all)[1]

  Tmat  <- matrix(T_arr, nrow = D)                         # D x D^2
  H_2d  <- matrix(H_all, nrow = N, ncol = D * D)           # N x D^2
  HW_2d <- (dW %*% H_2d) / N                               # B x D^2, Eq. (9)
  HT    <- H.inv %*% Tmat

  theta_rep <- matrix(NA_real_, B, D, dimnames = list(NULL, names(theta0)))
  for (i in seq_len(B)) {
    c_vec  <- C_mat[i, ]
    H_dw_i <- matrix(HW_2d[i, ], D, D)

    Bc <- drop(H.inv %*% H_dw_i %*% c_vec)
    Ac <- 0.5 * drop(HT %*% as.vector(tcrossprod(c_vec)))

    theta_rep[i, ] <- theta0 - c_vec + (Bc - Ac)
  }
  theta_rep
}


## Utilities functions

# numerical gradient of a scalar functional phi(theta)
num_grad_f <- function(f, theta, delta = 1e-5) {
  vapply(seq_along(theta), function(k) {
    tp <- theta; tp[k] <- tp[k] + delta
    tm <- theta; tm[k] <- tm[k] - delta
    tryCatch((f(tp) - f(tm)) / (2 * delta), error = function(e) NA_real_)
  }, numeric(1))
}

# symmetric delta-method interval
wald_ci <- function(f, theta0, vcov_mat, alpha = 0.05) {
  grad <- num_grad_f(f, theta0)
  se   <- sqrt(max(0, as.numeric(t(grad) %*% vcov_mat %*% grad)))
  est  <- unname(tryCatch(as.numeric(f(theta0)), error = function(e) NA_real_))
  c(lo = est - qnorm(1 - alpha / 2) * se,
    hi = est + qnorm(1 - alpha / 2) * se, se = se)
}

# equal-tailed percentile interval over replicate functional values
percentile_ci <- function(vals, alpha = 0.05, min_n = 40L) {
  vals <- vals[is.finite(vals)]
  if (length(vals) < min_n) return(c(lo = NA_real_, hi = NA_real_))
  q <- quantile(vals, c(alpha / 2, 1 - alpha / 2), names = FALSE, type = 7)
  c(lo = q[1], hi = q[2])
}

skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L) return(NA_real_)
  m <- mean(x); v <- mean((x - m)^2)
  if (!is.finite(v) || v <= 0) return(NA_real_)
  mean((x - m)^3) / v^1.5
}


# ---------------------------------------------------------------------
# Self-test of the lavaan conventions listed in the header
#
# (a) lavScores(scaling = TRUE) = -s_i(theta-hat) / N
# (b) information.observed      = sum_i H_i(theta-hat) / N
# (c) the finite-difference gradient reproduces that information
# (d) compute_T_tensor_grad() matches a direct third derivative
# (e) inverted.information.expected / N is lavaan's standard vcov,
#     the covariance matrix behind the Wald (Expected) comparator
# ---------------------------------------------------------------------
hoij_selftest <- function(tol_rel = 0.01, verbose = TRUE) {
  say <- function(...) if (verbose) cat(sprintf(...))
  say("-- hoij_core self-test (lavaan %s) --\n",
      as.character(packageVersion("lavaan")))

  fit <- lavaan::sem("f =~ x1 + x2 + x3",
                     data = lavaan::HolzingerSwineford1939,
                     estimator = "ML", se = "robust.huber.white")
  th0 <- lavaan::coef(fit, type = "free")
  D <- length(th0); N <- nrow(fit@Data@X[[1]]); h <- 1e-6

  ## (a) lavScores(scaling = TRUE) = -s_i / N
  S_num <- vapply(seq_len(D), function(k) {
    tp <- th0; tp[k] <- tp[k] + h
    tm <- th0; tm[k] <- tm[k] - h
    (compute_loglik_casewise(fit, tp) -
       compute_loglik_casewise(fit, tm)) / (2 * h)
  }, numeric(N))
  r_sc <- median(as.numeric(lavaan::lavScores(fit, scaling = TRUE)) /
                   as.numeric(S_num)) * N
  ok_a <- is.finite(r_sc) && abs(r_sc + 1) < tol_rel
  say("  (a) lavScores scale    : N * ratio = %+.6f (expect -1)  %s\n",
      r_sc, if (ok_a) "OK" else "FAIL")

  ## (b) information.observed = sum_i H_i / N  (also catches H_i == 0)
  H_sum <- apply(compute_all_H(fit, th0), c(2, 3), sum)
  H_obs <- lavaan::lavTech(fit, "information.observed")
  r_H <- median(as.numeric(H_obs) / as.numeric(H_sum)) * N
  ok_b <- is.finite(r_H) && abs(r_H - 1) < tol_rel
  say("  (b) observed info      : N * ratio = %+.6f (expect +1)  %s\n",
      r_H, if (ok_b) "OK" else "FAIL")

  ## (c) the finite differences reproduce lavaan's observed information
  grad_F <- make_grad_F(fit)
  chk <- tryCatch(check_gradient_hessian(grad_F, th0, H_obs),
                  error = function(e) NULL)
  ok_c <- !is.null(chk) && is.finite(chk$spread) && chk$spread < 0.01 &&
    abs(chk$ratio - 1) < tol_rel
  say("  (c) gradient Hessian   : ratio = %+.6f, spread = %.2e  %s\n",
      if (is.null(chk)) NA else chk$ratio,
      if (is.null(chk)) NA else chk$spread, if (ok_c) "OK" else "FAIL")

  ## (d) T array against a direct third derivative of -mean log-likelihood
  ok_d <- FALSE
  if (ok_c) {
    T_arr <- compute_T_tensor_grad(grad_F, th0)
    f_tot <- function(th) -sum(compute_loglik_casewise(fit, th)) / N
    hh <- 1e-3
    k <- which.max(abs(T_arr[cbind(1:D, 1:D, 1:D)]))
    pert <- function(sgn) { th <- th0; th[k] <- th[k] + sgn * hh; th }
    t_dir <- (f_tot(pert(2)) - 2 * f_tot(pert(1)) +
                2 * f_tot(pert(-1)) - f_tot(pert(-2))) / (2 * hh^3)
    r_T <- T_arr[k, k, k] / t_dir
    ok_d <- is.finite(r_T) && abs(r_T - 1) < 0.05
    say("  (d) third derivatives  : ratio = %+.6f (expect +1)   %s\n",
        r_T, if (ok_d) "OK" else "FAIL")
  } else {
    say("  (d) third derivatives  : skipped (check (c) failed)\n")
  }

  ## (e) expected-information vcov used by the Wald (Expected) comparator
  fit_std <- lavaan::sem("f =~ x1 + x2 + x3",
                         data = lavaan::HolzingerSwineford1939,
                         estimator = "ML", se = "standard")
  r_e <- max(abs(lavaan::lavTech(fit, "inverted.information.expected") / N /
                   lavaan::lavInspect(fit_std, "vcov") - 1))
  ok_e <- is.finite(r_e) && r_e < 1e-6
  say("  (e) expected info scale: max rel. deviation = %.2e  %s\n",
      r_e, if (ok_e) "OK" else "FAIL")

  ok <- ok_a && ok_b && ok_c && ok_d && ok_e
  if (ok) say("  self-test PASSED\n\n") else
    warning("hoij_core self-test FAILED; fix this before interpreting any ",
            "IJ1/HOIJ-2 output.", call. = FALSE)
  invisible(ok)
}
