# ============================================================
# HOIJ Worked Example - Table tab:intervals + Figure fig:shape
#
# Zes intervalmethoden x zeven functionalen, op dezelfde zwak-effect-
# populatie en dataset als HOIJ_worked_example_mediation_v2.R.
#
# Belangrijk:
# - abs_ab is opgenomen als niet-glad functionaal.
# - omega_vis, omega_text en omega_speed worden alle drie berekend.
# - Voor fig:shape wordt omega_speed gebruikt.
# - De tabel tab:intervals gebruikt omega_speed.
# - psi_speed (= speed~~speed, de ruwe storingsvariantie van de
#   speed-vergelijking) is toegevoegd als contrastcasus: g = identiteit,
#   dus alle scheefheid is parameter-eigen (kromming van de
#   score-vergelijking) in plaats van transformatie-geinduceerd zoals
#   bij ab/omega.
# - r2_speed vervangt r2_vis. Dit is het verklaarde-variantiefunctionaal
#   voor de speed-vergelijking: R^2_speed = 1 - psi_speed / Var(speed).
# ============================================================


# ─────────────────────────────────────────────────────────────
# 0. INSTELLINGEN
# ─────────────────────────────────────────────────────────────

SEED_DATA        <- 20260706
SEED_WEIGHTS     <- 20260707
SEED_MC          <- 20260708
N_EX             <- 500
B                <- 5000
R_MC             <- 5000
KAPPA_DAMP       <- 0.5
ALPHA_CI         <- 0.05
EFFECT_PARS_ZWAK <- c(a = 0.05, b = 0.10)

out_dir <- "hoij_we1_intervals_output"
if (!dir.exists(out_dir)) dir.create(out_dir)


# ─────────────────────────────────────────────────────────────
# 1. PACKAGES
# ─────────────────────────────────────────────────────────────

if (!requireNamespace("lavaan", quietly = TRUE)) install.packages("lavaan")
library(lavaan)

cat("lavaan", as.character(packageVersion("lavaan")), "geladen\n")


# ─────────────────────────────────────────────────────────────
# 2. HOIJ-KERNFUNCTIES
# ─────────────────────────────────────────────────────────────

compute_loglik_casewise <- function(fit, theta) {
  X <- fit@Data@X[[1]]
  N <- nrow(X)
  p <- ncol(X)

  glist <- lavaan:::lav_model_x2glist(fit@Model, x = theta)
  implied <- lavaan:::lav_model_implied(fit@Model, glist = glist)

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
  D <- length(theta0)
  N <- nrow(fit@Data@X[[1]])
  J_array <- array(0, dim = c(N, D, D))

  ll_0 <- compute_loglik_casewise(fit, theta0)

  for (k in seq_len(D)) {
    for (l in k:D) {
      if (k == l) {
        tp <- theta0
        tm <- theta0
        tp[k] <- tp[k] + delta
        tm[k] <- tm[k] - delta

        ll_p <- compute_loglik_casewise(fit, tp)
        ll_m <- compute_loglik_casewise(fit, tm)

        J_array[, k, k] <- -((ll_p - 2 * ll_0 + ll_m) / delta^2)
      } else {
        tpp <- theta0
        tpm <- theta0
        tmp_ <- theta0
        tmm <- theta0

        tpp[k] <- tpp[k] + delta
        tpp[l] <- tpp[l] + delta

        tpm[k] <- tpm[k] + delta
        tpm[l] <- tpm[l] - delta

        tmp_[k] <- tmp_[k] - delta
        tmp_[l] <- tmp_[l] + delta

        tmm[k] <- tmm[k] - delta
        tmm[l] <- tmm[l] - delta

        J_array[, k, l] <- -(
          compute_loglik_casewise(fit, tpp) -
            compute_loglik_casewise(fit, tpm) -
            compute_loglik_casewise(fit, tmp_) +
            compute_loglik_casewise(fit, tmm)
        ) / (4 * delta^2)

        J_array[, l, k] <- J_array[, k, l]
      }
    }
  }

  J_array
}


make_grad_F <- function(fit) {
  lavmodel        <- fit@Model
  lavsamplestats  <- fit@SampleStats
  lavdata         <- fit@Data
  lavcache        <- fit@Cache

  function(theta) {
    glist <- lavaan:::lav_model_x2glist(lavmodel, x = theta)

    as.numeric(
      lavaan:::lav_model_grad(
        lavmodel       = lavmodel,
        glist          = glist,
        lavsamplestats = lavsamplestats,
        lavdata        = lavdata,
        lavcache       = lavcache
      )
    )
  }
}


calibrate_alpha <- function(grad_F, theta0, H_observed, h = 1e-5) {
  D_loc <- length(theta0)
  H_grad <- matrix(NA_real_, D_loc, D_loc)

  for (k in seq_len(D_loc)) {
    tp <- theta0
    tm <- theta0
    tp[k] <- tp[k] + h
    tm[k] <- tm[k] - h

    H_grad[, k] <- (grad_F(tp) - grad_F(tm)) / (2 * h)
  }

  H_grad <- (H_grad + t(H_grad)) / 2

  idx <- abs(H_grad) > 1e-6 * max(abs(H_grad))
  ratio <- as.numeric(H_observed)[idx] / as.numeric(H_grad)[idx]

  alpha <- median(ratio)
  spread <- max(abs(ratio / alpha - 1))

  list(alpha = alpha, spread = spread)
}


compute_T_tensor_grad <- function(grad_F, theta, alpha, h = 1e-4) {
  D_loc <- length(theta)
  T_arr <- array(0, dim = c(D_loc, D_loc, D_loc))
  g0 <- grad_F(theta)

  for (l in seq_len(D_loc)) {
    for (m in l:D_loc) {
      if (l == m) {
        tp <- theta
        tm <- theta
        tp[l] <- tp[l] + h
        tm[l] <- tm[l] - h

        col <- (grad_F(tp) - 2 * g0 + grad_F(tm)) / h^2
      } else {
        tpp <- theta
        tpm <- theta
        tmp_ <- theta
        tmm <- theta

        tpp[l] <- tpp[l] + h
        tpp[m] <- tpp[m] + h

        tpm[l] <- tpm[l] + h
        tpm[m] <- tpm[m] - h

        tmp_[l] <- tmp_[l] - h
        tmp_[m] <- tmp_[m] + h

        tmm[l] <- tmm[l] - h
        tmm[m] <- tmm[m] - h

        col <- (
          grad_F(tpp) -
            grad_F(tpm) -
            grad_F(tmp_) +
            grad_F(tmm)
        ) / (4 * h^2)
      }

      T_arr[, l, m] <- alpha * col
      T_arr[, m, l] <- alpha * col
    }
  }

  T_arr <- (
    T_arr +
      aperm(T_arr, c(2, 1, 3)) +
      aperm(T_arr, c(3, 2, 1)) +
      aperm(T_arr, c(1, 3, 2)) +
      aperm(T_arr, c(2, 3, 1)) +
      aperm(T_arr, c(3, 1, 2))
  ) / 6

  T_arr
}


num_grad_f <- function(f, theta, delta = 1e-5) {
  sapply(seq_along(theta), function(k) {
    tp <- theta
    tm <- theta
    tp[k] <- tp[k] + delta
    tm[k] <- tm[k] - delta

    tryCatch(
      (f(tp) - f(tm)) / (2 * delta),
      error = function(e) NA_real_
    )
  })
}


skewness <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) < 3L) return(NA_real_)

  m <- mean(x)
  v <- mean((x - m)^2)

  if (!is.finite(v) || v <= 0) return(NA_real_)

  mean((x - m)^3) / v^1.5
}


# ─────────────────────────────────────────────────────────────
# 3. MODEL, POPULATIE EN DATA
# ─────────────────────────────────────────────────────────────

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

# Zwak-effect-override: alleen a en b worden overschreven.
pt_pop$ustart[pt_pop$label == "a"] <- EFFECT_PARS_ZWAK["a"]
pt_pop$ustart[pt_pop$label == "b"] <- EFFECT_PARS_ZWAK["b"]
pt_pop$est[pt_pop$label == "a"]    <- EFFECT_PARS_ZWAK["a"]
pt_pop$est[pt_pop$label == "b"]    <- EFFECT_PARS_ZWAK["b"]

set.seed(SEED_DATA)
dat <- simulateData(pt_pop, sample.nobs = N_EX)

fit <- sem(model_med, data = dat, se = "none", estimator = "ML")
stopifnot(lavInspect(fit, "converged"))

theta0   <- coef(fit, type = "free")
th_names <- names(theta0)
D        <- length(theta0)
N        <- nrow(dat)

stopifnot(D == 21, N == N_EX)

required_parameters <- c(
  "a", "b", "c",
  "textual~~textual",
  "speed~~speed",
  "visual~~visual",
  "x1~~x1"
)
missing_parameters <- setdiff(required_parameters, th_names)
if (length(missing_parameters) > 0L) {
  stop(sprintf(
    "De volgende vereiste vrije parameters ontbreken: %s",
    paste(missing_parameters, collapse = ", ")
  ))
}

cat(sprintf(
  "Analysemodel gefit: D = %d, N = %d\n",
  D, N
))
cat(sprintf(
  "a_hat = %.4f, b_hat = %.4f, ab_hat = %.4f, psi_speed_hat = %.4f\n",
  theta0["a"], theta0["b"], theta0["a"] * theta0["b"], theta0["speed~~speed"]
))


# ─────────────────────────────────────────────────────────────
# 4. FUNCTIONALEN: BASIS + DRIE OMEGA'S
# ─────────────────────────────────────────────────────────────

var_textual_vec <- function(mat, nms) {
  mat[, nms == "textual~~textual"]
}

var_speed_vec <- function(mat, nms) {
  a     <- mat[, nms == "a"]
  phi_t <- mat[, nms == "textual~~textual"]
  psi1  <- mat[, nms == "speed~~speed"]

  a^2 * phi_t + psi1
}

var_visual_vec <- function(mat, nms) {
  a     <- mat[, nms == "a"]
  b     <- mat[, nms == "b"]
  cc    <- mat[, nms == "c"]
  phi_t <- mat[, nms == "textual~~textual"]
  psi1  <- mat[, nms == "speed~~speed"]
  psi2  <- mat[, nms == "visual~~visual"]

  var_speed <- a^2 * phi_t + psi1
  b^2 * var_speed + cc^2 * phi_t + 2 * b * cc * a * phi_t + psi2
}

get_par_vec <- function(mat, nms, nm) {
  idx <- which(nms == nm)
  if (length(idx) != 1L) {
    stop(sprintf(
      "Parameter %s moet exact 1 keer voorkomen, maar komt %d keer voor.",
      nm, length(idx)
    ))
  }

  as.numeric(mat[, idx, drop = TRUE])
}

get_par_mat <- function(mat, nms, nm_vec) {
  idx <- vapply(nm_vec, function(nm) {
    ii <- which(nms == nm)
    if (length(ii) != 1L) {
      stop(sprintf(
        "Parameter %s moet exact 1 keer voorkomen, maar komt %d keer voor.",
        nm, length(ii)
      ))
    }
    ii
  }, integer(1))

  out <- mat[, idx, drop = FALSE]
  storage.mode(out) <- "double"
  out
}

omega_factor_vec <- function(mat, nms, factor, items, var_fun) {
  loading_names <- paste0(factor, "=~", items[-1])
  resid_names   <- paste0(items, "~~", items)

  # Gebruik expliciet drop = FALSE. Bij nrow(mat) == 1 maakt sapply()
  # anders van twee vrije loadings een vector van lengte 2, waarna
  # rowSums() ten onrechte twee rijen ziet.
  lambda_free <- get_par_mat(mat, nms, loading_names)
  theta_free  <- get_par_mat(mat, nms, resid_names)

  lambda_sum <- 1 + rowSums(lambda_free)
  theta_sum  <- rowSums(theta_free)
  vv         <- as.numeric(var_fun(mat, nms))

  if (length(vv) != nrow(mat)) {
    stop(sprintf(
      "var_fun voor %s geeft lengte %d terug, maar nrow(mat) = %d.",
      factor, length(vv), nrow(mat)
    ))
  }

  denom <- lambda_sum^2 * vv + theta_sum
  out <- lambda_sum^2 * vv / denom

  out[is.na(denom) | denom < 1e-10 | is.na(vv) | vv < 1e-10] <- NA_real_
  as.numeric(out)
}

functionals_all <- list(
  ab = function(mat, nms) {
    mat[, nms == "a"] * mat[, nms == "b"]
  },

  abs_ab = function(mat, nms) {
    abs(mat[, nms == "a"] * mat[, nms == "b"])
  },

  theta11 = function(mat, nms) {
    mat[, nms == "x1~~x1"]
  },

  # Rauwe storingsvariantie van de speed-vergelijking (g = identiteit).
  psi_speed = function(mat, nms) {
    mat[, nms == "speed~~speed"]
  },

  # Verklaarde variantie in de speed-vergelijking.
  # In dit model: Var(speed) = a^2 Var(textual) + psi_speed.
  r2_speed = function(mat, nms) {
    psi1 <- mat[, nms == "speed~~speed"]
    vv   <- var_speed_vec(mat, nms)

    r2 <- 1 - psi1 / vv
    r2[is.na(vv) | vv < 1e-10] <- NA_real_

    r2
  },

  pm = function(mat, nms) {
    ab_val <- mat[, nms == "a"] * mat[, nms == "b"]
    denom  <- ab_val + mat[, nms == "c"]

    res <- ab_val / denom
    res[is.na(denom) | abs(denom) < 1e-6] <- NA_real_

    res
  },

  omega_vis = function(mat, nms) {
    omega_factor_vec(
      mat = mat,
      nms = nms,
      factor = "visual",
      items = c("x1", "x2", "x3"),
      var_fun = var_visual_vec
    )
  },

  omega_text = function(mat, nms) {
    omega_factor_vec(
      mat = mat,
      nms = nms,
      factor = "textual",
      items = c("x4", "x5", "x6"),
      var_fun = var_textual_vec
    )
  },

  omega_speed = function(mat, nms) {
    omega_factor_vec(
      mat = mat,
      nms = nms,
      factor = "speed",
      items = c("x7", "x8", "x9"),
      var_fun = var_speed_vec
    )
  }
)

fn_all_names <- names(functionals_all)

functionals_all_scalar <- lapply(functionals_all, function(fv) {
  function(th) {
    val <- as.numeric(
      fv(
        matrix(th, nrow = 1, dimnames = list(NULL, names(th))),
        names(th)
      )
    )

    if (length(val) != 1L) {
      stop(sprintf(
        "Scalar functional gaf lengte %d terug in plaats van lengte 1.",
        length(val)
      ))
    }

    val
  }
})

fn_all_hat <- vapply(fn_all_names, function(fn_nm) {
  functionals_all_scalar[[fn_nm]](theta0)
}, numeric(1))

stopifnot(is.numeric(fn_all_hat), length(fn_all_hat) == length(fn_all_names))

cat("\nPuntschattingen alle functionalen:\n")
print(round(fn_all_hat, 4))


# ─────────────────────────────────────────────────────────────
# 5. V_inf EN V_hw VIA lavaan
# ─────────────────────────────────────────────────────────────

fit_inf <- tryCatch(
  sem(
    model_med,
    data        = dat,
    estimator   = "ML",
    se          = "standard",
    information = "expected"
  ),
  error = function(e) NULL
)

if (is.null(fit_inf) || !lavInspect(fit_inf, "converged")) {
  stop("Wald-delta (Inf)-fit met expected information mislukt.")
}

V_inf <- lavInspect(fit_inf, "vcov")[th_names, th_names, drop = FALSE]

fit_hw <- tryCatch(
  sem(
    model_med,
    data      = dat,
    estimator = "ML",
    se        = "robust.huber.white"
  ),
  error = function(e) NULL
)

if (is.null(fit_hw) || !lavInspect(fit_hw, "converged")) {
  stop("Wald-delta (HW)-fit met Huber-White sandwich mislukt.")
}

V_hw <- lavInspect(fit_hw, "vcov")[th_names, th_names, drop = FALSE]

stopifnot(max(abs(coef(fit_inf, type = "free")[th_names] - theta0[th_names])) < 1e-6)
stopifnot(max(abs(coef(fit_hw,  type = "free")[th_names] - theta0[th_names])) < 1e-6)

z_crit <- qnorm(1 - ALPHA_CI / 2)


# ─────────────────────────────────────────────────────────────
# 6. HOIJ-BASISGROOTHEDEN
# ─────────────────────────────────────────────────────────────

Scores <- lavScores(fit, scaling = TRUE)
H.inv  <- lavTech(fit, "inverted.information.observed")
H_obs  <- lavTech(fit, "information.observed")

if (is.null(Scores) || is.null(H.inv) || is.null(H_obs)) {
  stop("Scores of observed information niet beschikbaar op deze lavaan-versie.")
}

dimnames(H.inv) <- list(th_names, th_names)

grad_F <- make_grad_F(fit)
cal <- calibrate_alpha(grad_F, theta0, H_obs)

if (!is.finite(cal$alpha) || cal$spread > 0.1) {
  stop(sprintf("Alpha-kalibratie mislukt (spread = %.3g).", cal$spread))
}

J_all <- compute_all_J(fit, theta0, delta = 1e-5)
T_arr <- compute_T_tensor_grad(grad_F, theta0, cal$alpha)


# ─────────────────────────────────────────────────────────────
# 7. GEDEELDE GEWICHTEN + EXACTE BOOTSTRAP
# ─────────────────────────────────────────────────────────────

set.seed(SEED_WEIGHTS)
W_counts <- t(rmultinom(B, size = N, prob = rep(1 / N, N)))
DW <- W_counts - 1L

cat(sprintf("\nExacte bootstrap: %d herfits...\n", B))

boot_th <- matrix(NA_real_, B, D, dimnames = list(NULL, th_names))

for (bb in seq_len(B)) {
  idx <- rep.int(seq_len(N), W_counts[bb, ])
  dat_b <- dat[idx, , drop = FALSE]

  fit_b <- tryCatch(
    sem(
      model_med,
      data      = dat_b,
      se        = "none",
      estimator = "ML",
      start     = fit
    ),
    error = function(e) NULL
  )

  if (!is.null(fit_b) && lavInspect(fit_b, "converged")) {
    th_b <- coef(fit_b, type = "free")
    if (length(th_b) == D) boot_th[bb, ] <- th_b
  }

  if (bb %% 1000 == 0) {
    cat(sprintf("  %d / %d\n", bb, B))
  }
}

valid <- stats::complete.cases(boot_th)

cat(sprintf(
  "Herfits voltooid: %d geconvergeerd, %d gefaald (%.2f%%)\n",
  sum(valid), sum(!valid), 100 * mean(!valid)
))


# ─────────────────────────────────────────────────────────────
# 8. IJ1 EN HOIJ-2
# ─────────────────────────────────────────────────────────────

G_mat <- DW %*% Scores
C_mat <- G_mat %*% H.inv

ij_th <- sweep(-C_mat, 2, theta0, "+")
colnames(ij_th) <- th_names

Tmat     <- matrix(T_arr, nrow = D)
J_all_2d <- matrix(J_all, nrow = N, ncol = D * D)
JW_2d    <- (DW %*% J_all_2d) / N
HT       <- H.inv %*% Tmat

hoij_th <- matrix(NA_real_, B, D, dimnames = list(NULL, th_names))
s_vec   <- rep(NA_real_, B)

for (i in seq_len(B)) {
  c_vec  <- C_mat[i, ]
  J_dw_i <- matrix(JW_2d[i, ], D, D)

  Bc <- drop(H.inv %*% J_dw_i %*% c_vec)

  kron_cc <- as.vector(tcrossprod(c_vec))
  Ac <- 0.5 * drop(HT %*% kron_cc)

  d1 <- -c_vec
  d2 <- Bc - Ac

  n1 <- sqrt(sum(d1^2))
  n2 <- sqrt(sum(d2^2))

  s <- if (n2 > 0) min(1, KAPPA_DAMP * n1 / n2) else 1

  s_vec[i] <- s
  hoij_th[i, ] <- theta0 + d1 + s * d2
}

cat(sprintf(
  "Demping HOIJ-2: fractie s < 1 = %.3f, gemiddelde s = %.3f\n",
  mean(s_vec < 1), mean(s_vec)
))


# ─────────────────────────────────────────────────────────────
# 9. MONTE CARLO (HW)
# ─────────────────────────────────────────────────────────────

set.seed(SEED_MC)

L_mc <- t(chol(V_hw + diag(1e-10, D)))
Z_mc <- matrix(rnorm(R_MC * D), R_MC, D)
mc_th <- sweep(Z_mc %*% t(L_mc), 2, theta0, "+")

colnames(mc_th) <- th_names


# ─────────────────────────────────────────────────────────────
# 10. OMEGA-DIAGNOSTIEK EN KEUZE VOOR FIGUUR/TABEL
# ─────────────────────────────────────────────────────────────

omega_candidates <- c("omega_vis", "omega_text", "omega_speed")

omega_shape_diag <- do.call(rbind, lapply(omega_candidates, function(fn_nm) {
  f <- functionals_all[[fn_nm]]

  v_boot <- f(boot_th[valid, , drop = FALSE], th_names)
  v_mc   <- f(mc_th, th_names)

  v_boot <- v_boot[is.finite(v_boot)]
  v_mc   <- v_mc[is.finite(v_mc)]

  skew_boot <- skewness(v_boot)
  skew_mc   <- skewness(v_mc)

  data.frame(
    functional       = fn_nm,
    skew_boot        = skew_boot,
    skew_mc          = skew_mc,
    abs_skew_boot    = abs(skew_boot),
    abs_diff_boot_mc = abs(skew_boot - skew_mc),
    n_boot           = length(v_boot),
    n_mc             = length(v_mc),
    stringsAsFactors = FALSE
  )
}))

omega_shape_diag <- omega_shape_diag[
  order(-omega_shape_diag$abs_diff_boot_mc, -omega_shape_diag$abs_skew_boot),
]

cat("\nOmega shape diagnostics:\n")
print(omega_shape_diag, row.names = FALSE, digits = 4)

# Handmatige keuze voor manuscriptfiguur en tabellen:
# omega_speed wordt gekozen omdat deze de speed-factor betreft en in de
# huidige run ook de grootste discrepantie tussen bootstrap- en MC-shape liet zien.
omega_for_shape <- "omega_speed"

omega_selected_row <- omega_shape_diag[
  omega_shape_diag$functional == omega_for_shape,
  ,
  drop = FALSE
]

cat(sprintf(
  "\nGekozen omega voor figuur en tabel: %s (bootstrap-skewness = %.3f; |boot - MC skew| = %.3f)\n",
  omega_for_shape,
  omega_selected_row$skew_boot,
  omega_selected_row$abs_diff_boot_mc
))

omega_label_expr <- expression(omega[speed])
omega_tex_label  <- "$\\omega_{\\mathrm{speed}}$"


# Zeven functionalen voor de manuscript-tabel:
# ab, abs_ab, theta11, psi_speed, r2_speed, pm, en omega_speed.
# psi_speed is g = identiteit. r2_speed is de verklaarde variantie in
# dezelfde speed-vergelijking.
fn_names <- c("ab", "abs_ab", "theta11", "psi_speed", "r2_speed", "pm", omega_for_shape)

functionals_vec <- functionals_all[fn_names]
functionals_scalar <- functionals_all_scalar[fn_names]
fn_hat <- fn_all_hat[fn_names]

stopifnot(
  length(fn_names) == 7L,
  anyDuplicated(fn_names) == 0L,
  "psi_speed" %in% fn_names,
  "r2_speed" %in% fn_names,
  !"r2_vis" %in% fn_names,
  all(fn_names %in% names(functionals_all))
)

psi_speed_check <- functionals_vec$psi_speed(
  matrix(theta0, nrow = 1, dimnames = list(NULL, th_names)),
  th_names
)
if (!isTRUE(all.equal(
  as.numeric(psi_speed_check),
  as.numeric(theta0["speed~~speed"]),
  tolerance = 1e-12
))) {
  stop("psi_speed is niet gelijk aan de vrije parameter speed~~speed.")
}

r2_speed_check <- functionals_vec$r2_speed(
  matrix(theta0, nrow = 1, dimnames = list(NULL, th_names)),
  th_names
)
r2_speed_manual <- 1 - as.numeric(theta0["speed~~speed"]) /
  as.numeric(var_speed_vec(matrix(theta0, nrow = 1, dimnames = list(NULL, th_names)), th_names))
if (!isTRUE(all.equal(
  as.numeric(r2_speed_check),
  as.numeric(r2_speed_manual),
  tolerance = 1e-12
))) {
  stop("r2_speed komt niet overeen met 1 - speed~~speed / Var(speed).")
}


# ─────────────────────────────────────────────────────────────
# 11. TABEL tab:intervals VULLEN
# ─────────────────────────────────────────────────────────────

qlo <- ALPHA_CI / 2
qhi <- 1 - ALPHA_CI / 2

method_order <- c("wald_inf", "wald_hw", "mc_hw", "ij1", "hoij2", "boot")

method_labels <- c(
  wald_inf = "Wald--delta (Inf)",
  wald_hw  = "Wald--delta (HW)",
  mc_hw    = "Monte Carlo (HW)",
  ij1      = "IJ1 percentile",
  hoij2    = "HOIJ-2 percentile",
  boot     = "Bootstrap percentile"
)

rows <- list()

for (fn_nm in fn_names) {
  f <- functionals_vec[[fn_nm]]
  fscalar <- functionals_scalar[[fn_nm]]
  est <- fn_hat[fn_nm]

  grad_inf <- num_grad_f(fscalar, theta0)
  se_inf <- sqrt(max(0, as.numeric(t(grad_inf) %*% V_inf %*% grad_inf)))

  rows[[length(rows) + 1]] <- data.frame(
    functional = fn_nm,
    method = "wald_inf",
    SE = se_inf,
    skewness = 0,
    lo = est - z_crit * se_inf,
    hi = est + z_crit * se_inf,
    stringsAsFactors = FALSE
  )

  grad_hw <- num_grad_f(fscalar, theta0)
  se_hw <- sqrt(max(0, as.numeric(t(grad_hw) %*% V_hw %*% grad_hw)))

  rows[[length(rows) + 1]] <- data.frame(
    functional = fn_nm,
    method = "wald_hw",
    SE = se_hw,
    skewness = 0,
    lo = est - z_crit * se_hw,
    hi = est + z_crit * se_hw,
    stringsAsFactors = FALSE
  )

  draws_list <- list(
    mc_hw = mc_th,
    ij1   = ij_th[valid, , drop = FALSE],
    hoij2 = hoij_th[valid, , drop = FALSE],
    boot  = boot_th[valid, , drop = FALSE]
  )

  for (m in names(draws_list)) {
    v <- f(draws_list[[m]], th_names)
    vf <- v[is.finite(v)]

    rows[[length(rows) + 1]] <- data.frame(
      functional = fn_nm,
      method = m,
      SE = if (length(vf) > 1) sd(vf) else NA_real_,
      skewness = if (length(vf) >= 40) skewness(vf) else NA_real_,
      lo = if (length(vf) >= 40) quantile(vf, qlo, names = FALSE) else NA_real_,
      hi = if (length(vf) >= 40) quantile(vf, qhi, names = FALSE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
}

tab_intervals <- do.call(rbind, rows)
tab_intervals$method_label <- method_labels[tab_intervals$method]

tab_intervals <- tab_intervals[
  order(
    match(tab_intervals$functional, fn_names),
    match(tab_intervals$method, method_order)
  ),
]

cat("\n══════════ tab:intervals ══════════\n")
print(
  tab_intervals[, c("functional", "method_label", "SE", "skewness", "lo", "hi")],
  row.names = FALSE,
  digits = 4
)


# ─────────────────────────────────────────────────────────────
# 12. LATEX-REGELS GENEREREN
# ─────────────────────────────────────────────────────────────

fmt <- function(x, d = 3) {
  if (is.na(x)) "--" else formatC(x, digits = d, format = "f")
}

functional_tex <- c(
  ab          = "$ab$",
  abs_ab      = "$|ab|$",
  theta11     = "$\\theta_{11}$",
  psi_speed   = "$\\psi_{\\mathrm{speed}}$",
  r2_speed    = "$R^2_{\\mathrm{speed}}$",
  pm          = "$P_M$",
  omega_vis   = "$\\omega_{\\mathrm{vis}}$",
  omega_text  = "$\\omega_{\\mathrm{text}}$",
  omega_speed = "$\\omega_{\\mathrm{speed}}$"
)

cat("\n══════════ LaTeX-regels tab:intervals ══════════\n")

for (fn_nm in fn_names) {
  cat(sprintf("%% %s\n", fn_nm))

  sub <- tab_intervals[tab_intervals$functional == fn_nm, ]

  first <- TRUE
  for (m in method_order) {
    r <- sub[sub$method == m, ]

    fn_print <- if (first) functional_tex[fn_nm] else ""
    first <- FALSE

    cat(sprintf(
      "%s & %s & %s & %s & %s & %s \\\\\n",
      fn_print,
      method_labels[m],
      fmt(r$SE, 3),
      fmt(r$skewness, 3),
      fmt(r$lo, 3),
      fmt(r$hi, 3)
    ))
  }

  if (fn_nm != tail(fn_names, 1)) {
    cat("\\midrule\n")
  }
}


# ─────────────────────────────────────────────────────────────
# 13. FIGUUR fig:shape
# ─────────────────────────────────────────────────────────────

plot_panel <- function(fn_nm, xlab) {
  f <- functionals_all[[fn_nm]]

  v_boot  <- f(boot_th[valid, , drop = FALSE], th_names)
  v_hoij2 <- f(hoij_th[valid, , drop = FALSE], th_names)
  v_ij1   <- f(ij_th[valid, , drop = FALSE], th_names)
  v_mc    <- f(mc_th, th_names)
  est     <- fn_all_hat[fn_nm]

  d_boot  <- density(v_boot[is.finite(v_boot)])
  d_hoij2 <- density(v_hoij2[is.finite(v_hoij2)])
  d_ij1   <- density(v_ij1[is.finite(v_ij1)])
  d_mc    <- density(v_mc[is.finite(v_mc)])

  xr <- range(d_boot$x, d_hoij2$x, d_ij1$x, d_mc$x)
  yr <- range(0, d_boot$y, d_hoij2$y, d_ij1$y, d_mc$y)

  plot(
    NA,
    xlim = xr,
    ylim = yr,
    xlab = xlab,
    ylab = "Density",
    main = "",
    bty = "l"
  )

  polygon(
    d_boot$x,
    d_boot$y,
    col = adjustcolor("grey60", alpha.f = 0.5),
    border = NA
  )

  lines(d_hoij2$x, d_hoij2$y, lwd = 2, lty = 1)
  lines(d_ij1$x, d_ij1$y, lwd = 2, lty = 2)
  lines(d_mc$x, d_mc$y, lwd = 2, lty = 3)

  abline(v = est, lwd = 1, col = "black")
}

pdf(file.path(out_dir, "fig-shape.pdf"), width = 13, height = 4.2)
par(mfrow = c(1, 3), mar = c(4, 4, 1, 1))

plot_panel("ab", expression(italic(ab)))

legend(
  "topright",
  legend = c("Bootstrap", "HOIJ-2", "IJ1", "Monte Carlo"),
  fill = c(adjustcolor("grey60", alpha.f = 0.5), NA, NA, NA),
  border = c("grey60", NA, NA, NA),
  lty = c(NA, 1, 2, 3),
  lwd = c(NA, 2, 2, 2),
  merge = FALSE,
  bty = "n",
  cex = 0.8
)

# psi_speed: g = identiteit, dus MC (HW) valt hier per constructie
# samen met Wald--delta (HW). Elk verschil met bootstrap/HOIJ-2 in dit
# paneel komt dus uit de reweighted estimator.
plot_panel("psi_speed", expression(psi[speed]))

plot_panel(omega_for_shape, omega_label_expr)

dev.off()

cat(sprintf(
  "\nFiguur opgeslagen: %s\n",
  file.path(out_dir, "fig-shape.pdf")
))


# ─────────────────────────────────────────────────────────────
# 14. PAIRWISE AGREEMENT EXACT BOOTSTRAP VS HOIJ-2
# ─────────────────────────────────────────────────────────────

pairwise_rows <- list()

for (fn_nm in fn_names) {
  f <- functionals_vec[[fn_nm]]

  v_boot  <- f(boot_th[valid, , drop = FALSE], th_names)
  v_hoij2 <- f(hoij_th[valid, , drop = FALSE], th_names)

  both <- is.finite(v_boot) & is.finite(v_hoij2)

  pairwise_rows[[length(pairwise_rows) + 1]] <- data.frame(
    functional = fn_nm,
    estimate = fn_hat[fn_nm],
    boot_lo = unname(quantile(v_boot[both], qlo)),
    boot_hi = unname(quantile(v_boot[both], qhi)),
    boot_skew = skewness(v_boot[both]),
    hoij2_lo = unname(quantile(v_hoij2[both], qlo)),
    hoij2_hi = unname(quantile(v_hoij2[both], qhi)),
    hoij2_skew = skewness(v_hoij2[both]),
    cor_boot_hoij2 = cor(v_boot[both], v_hoij2[both]),
    n_pair = sum(both),
    stringsAsFactors = FALSE
  )
}

tab_pairwise <- do.call(rbind, pairwise_rows)

cat("\n══════════ tab:boot-vs-hoij waarden ══════════\n")
print(tab_pairwise, row.names = FALSE, digits = 4)


# ─────────────────────────────────────────────────────────────
# 15. WEGSCHRIJVEN
# ─────────────────────────────────────────────────────────────

ts <- format(Sys.time(), "%Y%m%d_%H%M")

out_csv_intervals <- file.path(out_dir, sprintf("tab_intervals_%s.csv", ts))
out_csv_pairwise  <- file.path(out_dir, sprintf("tab_pairwise_boot_hoij_%s.csv", ts))
out_csv_omega     <- file.path(out_dir, sprintf("omega_shape_diag_%s.csv", ts))
out_rds           <- file.path(out_dir, sprintf("we1_intervals_draws_%s.rds", ts))

write.csv(tab_intervals, out_csv_intervals, row.names = FALSE)
write.csv(tab_pairwise, out_csv_pairwise, row.names = FALSE)
write.csv(omega_shape_diag, out_csv_omega, row.names = FALSE)

saveRDS(
  list(
    seed_data = SEED_DATA,
    seed_weights = SEED_WEIGHTS,
    seed_mc = SEED_MC,
    effect_pars = EFFECT_PARS_ZWAK,
    N = N,
    B = B,
    R_mc = R_MC,
    D = D,
    theta0 = theta0,
    th_names = th_names,
    fn_all_hat = fn_all_hat,
    fn_hat = fn_hat,
    fn_names = fn_names,
    omega_shape_diag = omega_shape_diag,
    omega_for_shape = omega_for_shape,
    valid = valid,
    boot_th = boot_th,
    hoij_th = hoij_th,
    ij_th = ij_th,
    mc_th = mc_th,
    V_inf = V_inf,
    V_hw = V_hw,
    tab_intervals = tab_intervals,
    tab_pairwise = tab_pairwise,
    W_counts = W_counts
  ),
  out_rds
)

cat(sprintf(
  "\nOpgeslagen:\n  %s\n  %s\n  %s\n  %s\n",
  out_csv_intervals,
  out_csv_pairwise,
  out_csv_omega,
  out_rds
))
