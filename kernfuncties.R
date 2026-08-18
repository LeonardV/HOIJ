compute_loglik_casewise <- function(fit, theta) {
  X <- fit@Data@X[[1]]
  N <- nrow(X); p <- ncol(X)
  GLIST <- lavaan:::lav_model_x2glist(fit@Model, x = theta)
  implied <- lavaan:::lav_model_implied(fit@Model, GLIST = GLIST)
  Sigma <- implied$cov[[1]]
  mu <- implied$mean[[1]]
  if (is.null(mu) || length(mu) == 0) mu <- colMeans(X)
  Sigma_inv <- tryCatch(solve(Sigma), error = function(e) MASS::ginv(Sigma))
  log_det <- determinant(Sigma, logarithm = TRUE)$modulus[1]
  const <- -0.5 * (p * log(2 * pi) + log_det)
  X_centered <- sweep(X, 2, mu, "-")
  quad_form <- rowSums((X_centered %*% Sigma_inv) * X_centered)
  as.numeric(const - 0.5 * quad_form)
}

compute_all_J <- function(fit, theta0, delta = 1e-5) {
  D <- length(theta0); N <- nrow(fit@Data@X[[1]])
  J_array <- array(0, dim = c(N, D, D))
  ll_0 <- compute_loglik_casewise(fit, theta0)
  for (k in 1:D) {
    for (l in k:D) {
      if (k == l) {
        tp <- theta0; tp[k] <- tp[k] + delta
        tm <- theta0; tm[k] <- tm[k] - delta
        ll_p <- compute_loglik_casewise(fit, tp)
        ll_m <- compute_loglik_casewise(fit, tm)
        J_array[, k, k] <- -((ll_p - 2*ll_0 + ll_m) / (delta^2))
      } else {
        tpp <- theta0; tpp[k] <- tpp[k]+delta; tpp[l] <- tpp[l]+delta
        tpm <- theta0; tpm[k] <- tpm[k]+delta; tpm[l] <- tpm[l]-delta
        tmp_ <- theta0; tmp_[k] <- tmp_[k]-delta; tmp_[l] <- tmp_[l]+delta
        tmm <- theta0; tmm[k] <- tmm[k]-delta; tmm[l] <- tmm[l]-delta
        J_array[, k, l] <- -((compute_loglik_casewise(fit, tpp) -
                                compute_loglik_casewise(fit, tpm) -
                                compute_loglik_casewise(fit, tmp_) +
                                compute_loglik_casewise(fit, tmm)) / (4*delta^2))
        J_array[, l, k] <- J_array[, k, l]
      }
    }
  }
  J_array
}

make_grad_F <- function(fit) {
  lavmodel       <- fit@Model
  lavsamplestats <- fit@SampleStats
  lavdata        <- fit@Data
  lavcache       <- fit@Cache
  function(theta) {
    GLIST <- lavaan:::lav_model_x2glist(lavmodel, x = theta)
    as.numeric(lavaan:::lav_model_gradient(
      lavmodel       = lavmodel,
      GLIST          = GLIST,
      lavsamplestats = lavsamplestats,
      lavdata        = lavdata,
      lavcache       = lavcache))
  }
}

calibrate_alpha <- function(grad_F, theta0, H_observed, h = 1e-5) {
  D_loc <- length(theta0)
  H_grad <- matrix(NA_real_, D_loc, D_loc)
  for (k in 1:D_loc) {
    tp <- theta0; tp[k] <- tp[k] + h
    tm <- theta0; tm[k] <- tm[k] - h
    H_grad[, k] <- (grad_F(tp) - grad_F(tm)) / (2 * h)
  }
  H_grad <- (H_grad + t(H_grad)) / 2
  idx    <- abs(H_grad) > 1e-6 * max(abs(H_grad))
  ratio  <- as.numeric(H_observed)[idx] / as.numeric(H_grad)[idx]
  alpha  <- median(ratio)
  spread <- max(abs(ratio / alpha - 1))
  list(alpha = alpha, spread = spread)
}

compute_T_tensor_grad <- function(grad_F, theta, alpha, h = 1e-4) {
  D_loc <- length(theta)
  T_arr <- array(0, dim = c(D_loc, D_loc, D_loc))
  g0    <- grad_F(theta)
  for (l in 1:D_loc) {
    for (m in l:D_loc) {
      if (l == m) {
        tp <- theta; tp[l] <- tp[l] + h
        tm <- theta; tm[l] <- tm[l] - h
        col <- (grad_F(tp) - 2 * g0 + grad_F(tm)) / h^2
      } else {
        tpp <- theta; tpp[l] <- tpp[l] + h; tpp[m] <- tpp[m] + h
        tpm <- theta; tpm[l] <- tpm[l] + h; tpm[m] <- tpm[m] - h
        tmp_ <- theta; tmp_[l] <- tmp_[l] - h; tmp_[m] <- tmp_[m] + h
        tmm <- theta; tmm[l] <- tmm[l] - h; tmm[m] <- tmm[m] - h
        col <- (grad_F(tpp) - grad_F(tpm) - grad_F(tmp_) + grad_F(tmm)) / (4*h^2)
      }
      T_arr[, l, m] <- alpha * col
      T_arr[, m, l] <- alpha * col
    }
  }
  T_arr <- (T_arr +
              aperm(T_arr, c(2, 1, 3)) + aperm(T_arr, c(3, 2, 1)) +
              aperm(T_arr, c(1, 3, 2)) + aperm(T_arr, c(2, 3, 1)) +
              aperm(T_arr, c(3, 1, 2))) / 6
  T_arr
}
