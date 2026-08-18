# ============================================================
# HOIJ Timing - tab:timing (mediation D=21 vs. bifactor D=27)
#
# [W19] Bifactor-model expliciet gekozen (niet in de manuscripttekst
#       gespecificeerd): orthogonaal bifactor-model op dezelfde 9
#       Holzinger-Swineford-items -- algemene factor g plus de drie
#       oorspronkelijke specifieke factoren (visual/textual/speed),
#       alle onderlinge covarianties op 0 (standaard bifactor-
#       identificatie). Telling: g-blok 8 vrije ladingen + 1 variantie
#       = 9; elk specifiek-factor-blok 2 vrije ladingen + 1 variantie
#       = 3 x 3 = 9; 9 residuele varianties. Totaal D = 9+9+9 = 27,
#       zoals in de tabel vereist. x1 dient als marker voor zowel g
#       als visual -- gebruikelijk maar kan Heywood-gevoeligheid
#       verhogen; vandaar de tryCatch/telling bij de herfits.
# [W20] Populatie = HS-MLE-schattingen, GEEN effect-override: dit
#       voorbeeld gaat over rekenkosten, niet over scheefheid, dus
#       geen reden om van de standaardpopulatie af te wijken.
# [W21] Timing per replicatie wordt gemeten op een representatieve
#       steekproef van N_TIMING replicaties (default 200), niet op de
#       volle B, en daarna geextrapoleerd naar B=1.000/5.000 via
#       vermenigvuldiging. Dat is precies de rekenregel die de tabel
#       zelf hanteert (setup + B x per-replicatiekosten), dus geen
#       verlies aan geldigheid -- alleen nodig omdat 5.000 volledige
#       herfits van het D=27-model onnodig lang zouden duren voor een
#       timingschatting.
# [W22] compute_T_tensor_grad doet exact 2*D^2 gradient-evaluaties
#       (D diagonaaltermen x 2 evals, D*(D-1)/2 buiten-diagonaaltermen
#       x 4 evals => 2D + 2D(D-1) = 2D^2), overeenkomstig de claim
#       "on the order of 2D^2 gradient evaluations" in de tekst.
#
# Wall-clock is hardware-afhankelijk; de complexiteitsclaim (1 fit +
# afgeleiden vs. B volledige herfits) blijft primair, zoals elders in
# het manuscript benadrukt.
# ============================================================


# ─────────────────────────────────────────────────────────────
# 0. INSTELLINGEN
# ─────────────────────────────────────────────────────────────

N_EX         <- 500
N_TIMING     <- 200          # aantal herfits/replicaties voor de timingschatting
N_FIT_TIMING <- 30           # aantal herhaalde fits voor stabiele single-fit timing
B_TARGETS    <- c(1000, 5000)
KAPPA_DAMP   <- 0.5
SEED_MED     <- 20260706
SEED_BIF     <- 20260810
out_dir      <- "hoij_timing_output"

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
        
        J_array[, k, k] <- -((ll_p - 2 * ll_0 + ll_m) / (delta^2))
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
          (
            compute_loglik_casewise(fit, tpp) -
              compute_loglik_casewise(fit, tpm) -
              compute_loglik_casewise(fit, tmp_) +
              compute_loglik_casewise(fit, tmm)
          ) / (4 * delta^2)
        )
        
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
    
    as.numeric(
      lavaan:::lav_model_gradient(
        lavmodel       = lavmodel,
        GLIST          = GLIST,
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


# Herhaalde fit-timing: stabieler dan één system.time()-meting.
# We gebruiken warm starts, omdat de exacte bootstrap-refits in dit script
# ook start = fit gebruiken. De rij in de tabel heet daarom expliciet
# "Median warm-start ML fit".
time_fit_repeated <- function(model_syntax, dat, start_fit,
                              std_lv = FALSE,
                              n_fit_timing = 30) {
  times <- rep(NA_real_, n_fit_timing)
  
  for (ii in seq_len(n_fit_timing)) {
    gc(FALSE)
    
    tt <- system.time({
      fit_tmp <- tryCatch(
        sem(
          model_syntax,
          data      = dat,
          se        = "none",
          estimator = "ML",
          std.lv    = std_lv,
          start     = start_fit
        ),
        error = function(e) NULL
      )
    })["elapsed"]
    
    if (!is.null(fit_tmp) && lavInspect(fit_tmp, "converged")) {
      times[ii] <- as.numeric(tt)
    }
  }
  
  out <- c(
    median = median(times, na.rm = TRUE),
    mean   = mean(times, na.rm = TRUE),
    min    = min(times, na.rm = TRUE),
    max    = max(times, na.rm = TRUE),
    n_fail = sum(is.na(times))
  )
  
  if (!is.finite(out["median"])) {
    stop("Repeated ML-fit timing failed: no converged repeated fits.")
  }
  
  out
}


# ─────────────────────────────────────────────────────────────
# 3. MODELDEFINITIES
# ─────────────────────────────────────────────────────────────

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
  visual ~~ 0*textual
  visual ~~ 0*speed
  textual ~~ 0*speed
'


# ─────────────────────────────────────────────────────────────
# 4. GENERIEKE TIMINGFUNCTIE VOOR EEN MODEL
# ─────────────────────────────────────────────────────────────

time_one_model <- function(model_syntax, D_expected, label,
                           seed_data, std_lv = FALSE) {
  cat(sprintf("\n========== %s (D=%d) ==========\n", label, D_expected))
  
  fit_pop <- sem(
    model_syntax,
    data   = HolzingerSwineford1939,
    se     = "none",
    std.lv = std_lv
  )
  
  if (!lavInspect(fit_pop, "converged")) {
    stop(sprintf("[%s] populatiefit niet geconvergeerd.", label))
  }
  
  pt_pop <- parTable(fit_pop)
  
  set.seed(seed_data)
  dat <- simulateData(pt_pop, sample.nobs = N_EX)
  
  # Analysefit: deze fit is nodig als referentiepunt voor de HOIJ-setup.
  # De tijd ervan wordt niet meer als één losse system.time()-meting
  # gerapporteerd; daarvoor gebruiken we hieronder repeated warm-start timing.
  fit <- sem(
    model_syntax,
    data      = dat,
    se        = "none",
    estimator = "ML",
    std.lv    = std_lv
  )
  
  if (!lavInspect(fit, "converged")) {
    stop(sprintf("[%s] analysefit op de voorbeelddataset niet geconvergeerd.", label))
  }
  
  theta0   <- coef(fit, type = "free")
  th_names <- names(theta0)
  D        <- length(theta0)
  N        <- nrow(dat)
  
  stopifnot(D == D_expected)
  
  # Stabielere single-fit timing: mediaan over meerdere warm-start fits.
  fit_timing <- time_fit_repeated(
    model_syntax  = model_syntax,
    dat           = dat,
    start_fit     = fit,
    std_lv        = std_lv,
    n_fit_timing  = N_FIT_TIMING
  )
  
  t_fit <- unname(fit_timing["median"])
  
  cat(sprintf(
    "D = %d bevestigd; N = %d; median warm-start fit = %.4f s (n=%d, fail=%d)\n",
    D, N, t_fit, N_FIT_TIMING, as.integer(fit_timing["n_fail"])
  ))
  
  # -- setup: scores + Jhat^-1 --
  t_sc <- system.time({
    Scores <- lavScores(fit, scaling = TRUE)
    H.inv  <- lavTech(fit, "inverted.information.observed")
    H_obs  <- lavTech(fit, "information.observed")
  })["elapsed"]
  
  if (is.null(Scores) || is.null(H.inv) || is.null(H_obs)) {
    stop(sprintf("[%s] Scores/observed information niet beschikbaar.", label))
  }
  
  dimnames(H.inv) <- list(th_names, th_names)
  
  grad_F <- make_grad_F(fit)
  
  t_cal <- system.time({
    cal <- calibrate_alpha(grad_F, theta0, H_obs)
  })["elapsed"]
  
  if (!is.finite(cal$alpha) || cal$spread > 0.1) {
    stop(sprintf("[%s] alpha-kalibratie mislukt (spread = %.3g).",
                 label, cal$spread))
  }
  
  # -- setup: casewise curvatures J_i --
  t_J <- system.time({
    J_all <- compute_all_J(fit, theta0, delta = 1e-5)
  })["elapsed"]
  
  # -- setup: derde-afgeleide-array Khat --
  t_T <- system.time({
    T_arr <- compute_T_tensor_grad(grad_F, theta0, cal$alpha)
  })["elapsed"]
  
  t_setup <- as.numeric(t_sc + t_cal + t_J + t_T)
  
  cat(sprintf(
    "setup: scores+Jhat^-1=%.4f s | J_i=%.4f s | Khat=%.4f s | totaal=%.4f s\n",
    as.numeric(t_sc + t_cal), as.numeric(t_J), as.numeric(t_T), t_setup
  ))
  
  # -- per-replicatiekosten: approximate bootstrap (HOIJ-2-lus) --
  set.seed(seed_data + 1)
  W_counts <- t(rmultinom(N_TIMING, size = N, prob = rep(1 / N, N)))
  DW <- W_counts - 1L
  
  t_loop <- system.time({
    G_mat <- DW %*% Scores
    C_mat <- G_mat %*% H.inv
    
    Tmat <- matrix(T_arr, nrow = D)
    J_all_2d <- matrix(J_all, nrow = N, ncol = D * D)
    JW_2d <- (DW %*% J_all_2d) / N
    HT <- H.inv %*% Tmat
    
    for (i in seq_len(N_TIMING)) {
      c_vec <- C_mat[i, ]
      J_dw_i <- matrix(JW_2d[i, ], D, D)
      
      Bc <- drop(H.inv %*% J_dw_i %*% c_vec)
      
      kron_cc <- as.vector(tcrossprod(c_vec))
      Ac <- 0.5 * drop(HT %*% kron_cc)
      
      d1 <- -c_vec
      d2 <- Bc - Ac
      
      n1 <- sqrt(sum(d1^2))
      n2 <- sqrt(sum(d2^2))
      
      s <- if (n2 > 0) min(1, KAPPA_DAMP * n1 / n2) else 1
      
      invisible(theta0 + d1 + s * d2)
    }
  })["elapsed"]
  
  t_rep_approx <- as.numeric(t_loop) / N_TIMING
  
  cat(sprintf(
    "approximate bootstrap: %.6f s/replicatie (op %d replicaties)\n",
    t_rep_approx, N_TIMING
  ))
  
  # -- per-replicatiekosten: exacte bootstrap (herfits, warme start) --
  t_refit_i <- rep(NA_real_, N_TIMING)
  n_fail <- 0L
  
  for (bb in seq_len(N_TIMING)) {
    idx <- rep.int(seq_len(N), W_counts[bb, ])
    dat_b <- dat[idx, , drop = FALSE]
    
    tt <- system.time({
      fit_b <- tryCatch(
        sem(
          model_syntax,
          data      = dat_b,
          se        = "none",
          estimator = "ML",
          start     = fit,
          std.lv    = std_lv
        ),
        error = function(e) NULL
      )
    })["elapsed"]
    
    if (!is.null(fit_b) && lavInspect(fit_b, "converged")) {
      t_refit_i[bb] <- as.numeric(tt)
    } else {
      n_fail <- n_fail + 1L
    }
  }
  
  if (sum(!is.na(t_refit_i)) < 0.5 * N_TIMING) {
    stop(sprintf(
      "[%s] te veel niet-geconvergeerde herfits (%d/%d) voor een betrouwbare timingschatting.",
      label, n_fail, N_TIMING
    ))
  }
  
  t_refit_mean <- mean(t_refit_i, na.rm = TRUE)
  
  cat(sprintf(
    "exacte bootstrap: %.6f s/replicatie (%d/%d geconvergeerd)\n",
    t_refit_mean, N_TIMING - n_fail, N_TIMING
  ))
  
  # -- totalen, speed-up, break-even --
  tot_approx <- setNames(
    t_setup + B_TARGETS * t_rep_approx,
    paste0("B", B_TARGETS)
  )
  
  tot_exact <- setNames(
    B_TARGETS * t_refit_mean,
    paste0("B", B_TARGETS)
  )
  
  speedup_1000 <- tot_exact["B1000"] / tot_approx["B1000"]
  
  Bstar <- if (t_refit_mean > t_rep_approx) {
    ceiling(t_setup / (t_refit_mean - t_rep_approx))
  } else {
    Inf
  }
  
  list(
    label = label,
    D = D,
    
    t_fit = t_fit,
    fit_timing = fit_timing,
    
    t_scores_hinv = as.numeric(t_sc + t_cal),
    t_J = as.numeric(t_J),
    t_T = as.numeric(t_T),
    t_setup = t_setup,
    
    t_rep_approx = t_rep_approx,
    tot_approx_1000 = tot_approx["B1000"],
    tot_approx_5000 = tot_approx["B5000"],
    
    t_rep_exact = t_refit_mean,
    tot_exact_1000 = tot_exact["B1000"],
    tot_exact_5000 = tot_exact["B5000"],
    
    speedup_1000 = speedup_1000,
    Bstar = Bstar,
    
    n_fail_exact = n_fail
  )
}


# ─────────────────────────────────────────────────────────────
# 5. UITVOEREN VOOR BEIDE MODELLEN
# ─────────────────────────────────────────────────────────────

res_med <- time_one_model(
  model_syntax = model_med,
  D_expected   = 21,
  label        = "Mediation",
  seed_data    = SEED_MED,
  std_lv       = FALSE
)

res_bif <- time_one_model(
  model_syntax = model_bifactor,
  D_expected   = 27,
  label        = "Bifactor",
  seed_data    = SEED_BIF,
  std_lv       = TRUE
)


# ─────────────────────────────────────────────────────────────
# 6. TABEL tab:timing SAMENSTELLEN EN LATEX GENEREREN
# ─────────────────────────────────────────────────────────────

fmt <- function(x, d = 3) {
  formatC(as.numeric(x), digits = d, format = "f")
}

cat("\n══════════ tab:timing ══════════\n")

row_labels <- c(
  "Median warm-start ML fit",
  "casewise scores and Jhat^-1",
  "casewise curvatures J_i",
  "third-derivative array Khat",
  "setup total",
  "per replication (approx)",
  "total, B=1000 (approx)",
  "total, B=5000 (approx)",
  "per replication (mean refit)",
  "total, B=1000 (exact)",
  "total, B=5000 (exact)",
  "Speed-up factor at B=1000",
  "Break-even B*"
)

get_vals <- function(r) {
  c(
    r$t_fit,
    r$t_scores_hinv,
    r$t_J,
    r$t_T,
    r$t_setup,
    r$t_rep_approx,
    r$tot_approx_1000,
    r$tot_approx_5000,
    r$t_rep_exact,
    r$tot_exact_1000,
    r$tot_exact_5000,
    r$speedup_1000,
    r$Bstar
  )
}

tab_timing <- data.frame(
  Row = row_labels,
  Mediation_D21 = get_vals(res_med),
  Bifactor_D27  = get_vals(res_bif)
)

print(tab_timing, row.names = FALSE, digits = 4)

cat("\n══════════ Fit timing diagnostics ══════════\n")
cat("Mediation repeated warm-start fit timing:\n")
print(res_med$fit_timing, digits = 4)
cat("Bifactor repeated warm-start fit timing:\n")
print(res_bif$fit_timing, digits = 4)

cat("\n══════════ LaTeX-regels (tab:timing) ══════════\n")

cat(sprintf(
  "Median warm-start ML fit & %s & %s \\\\\n",
  fmt(res_med$t_fit, 4), fmt(res_bif$t_fit, 4)
))

cat("\\addlinespace\n")

cat(sprintf(
  "\\quad casewise scores and $\\Jhat^{-1}$ & %s & %s \\\\\n",
  fmt(res_med$t_scores_hinv, 4), fmt(res_bif$t_scores_hinv, 4)
))

cat(sprintf(
  "\\quad casewise curvatures $J_i$ & %s & %s \\\\\n",
  fmt(res_med$t_J, 4), fmt(res_bif$t_J, 4)
))

cat(sprintf(
  "\\quad third-derivative array $\\Khat$ & %s & %s \\\\\n",
  fmt(res_med$t_T, 4), fmt(res_bif$t_T, 4)
))

cat(sprintf(
  "\\quad setup total & %s & %s \\\\\n",
  fmt(res_med$t_setup, 4), fmt(res_bif$t_setup, 4)
))

cat(sprintf(
  "\\quad per replication & %s & %s \\\\\n",
  fmt(res_med$t_rep_approx, 5), fmt(res_bif$t_rep_approx, 5)
))

cat(sprintf(
  "\\quad total, $B=1{,}000$ & %s & %s \\\\\n",
  fmt(res_med$tot_approx_1000, 3), fmt(res_bif$tot_approx_1000, 3)
))

cat(sprintf(
  "\\quad total, $B=5{,}000$ & %s & %s \\\\\n",
  fmt(res_med$tot_approx_5000, 3), fmt(res_bif$tot_approx_5000, 3)
))

cat("\\addlinespace\n")

cat(sprintf(
  "\\quad per replication (mean refit) & %s & %s \\\\\n",
  fmt(res_med$t_rep_exact, 5), fmt(res_bif$t_rep_exact, 5)
))

cat(sprintf(
  "\\quad total, $B=1{,}000$ & %s & %s \\\\\n",
  fmt(res_med$tot_exact_1000, 3), fmt(res_bif$tot_exact_1000, 3)
))

cat(sprintf(
  "\\quad total, $B=5{,}000$ & %s & %s \\\\\n",
  fmt(res_med$tot_exact_5000, 3), fmt(res_bif$tot_exact_5000, 3)
))

cat("\\addlinespace\n")

cat(sprintf(
  "Speed-up factor at $B=1{,}000$ & %s & %s \\\\\n",
  fmt(res_med$speedup_1000, 1), fmt(res_bif$speedup_1000, 1)
))

cat(sprintf(
  "Break-even $B^{*}$ & %s & %s \\\\\n",
  format(res_med$Bstar), format(res_bif$Bstar)
))

cat(sprintf(
  "\n[W20-check] mislukte exacte herfits tijdens timing: mediation=%d/%d, bifactor=%d/%d\n",
  res_med$n_fail_exact, N_TIMING, res_bif$n_fail_exact, N_TIMING
))

cat(
  "(Hoger dan verwacht bij bifactor kan wijzen op Heywood-gevoeligheid\n",
  "van het dubbele-marker-item x1; zie [W19].)\n",
  sep = ""
)


# ─────────────────────────────────────────────────────────────
# 7. WEGSCHRIJVEN
# ─────────────────────────────────────────────────────────────

ts <- format(Sys.time(), "%Y%m%d_%H%M")

out_csv <- file.path(out_dir, sprintf("tab_timing_%s.csv", ts))
out_rds <- file.path(out_dir, sprintf("timing_raw_%s.rds", ts))

write.csv(tab_timing, out_csv, row.names = FALSE)

saveRDS(
  list(
    res_med      = res_med,
    res_bif      = res_bif,
    N_TIMING     = N_TIMING,
    N_FIT_TIMING = N_FIT_TIMING,
    B_TARGETS    = B_TARGETS
  ),
  out_rds
)

cat(sprintf("\nOpgeslagen:\n  %s\n  %s\n", out_csv, out_rds))