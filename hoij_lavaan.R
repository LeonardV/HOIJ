# ============================================================
# hoij_lavaan(): Higher-Order Infinitesimal Jackknife voor lavaan
#
# Berekent standaardfouten en percentiel-betrouwbaarheidsintervallen
# voor (functies van) de vrije parameters van een gefit lavaan-model,
# ZONDER bootstrap-herfits. De bootstrapverdeling wordt benaderd met
# een eerste- (IJ1) of tweede-orde (HOIJ-2) Taylor-expansie van de
# schatter in de observatiegewichten, zoals in HOIJ_simulatiestudie §8.
#
# Per gewichtsvector w (multinomiaal, zoals een bootstrap-resample):
#   c    = H^{-1} S' (w - 1)                        (invloedsterm)
#   IJ1  : theta(w) ~ theta_hat - c
#   HOIJ2: d2 = H^{-1} J(dw) c - 1/2 H^{-1} T (c (x) c)
#          theta(w) ~ theta_hat - c + s*d2,
#          s = min(1, kappa*||c||/||d2||)           (trust-region-demping)
#
# met S = casewise scores (lavScores), H = geobserveerde informatie,
# J = casewise geobserveerde informatie (N x D x D) en T = derde-orde
# tensor van de loglikelihood (D x D x D).
#
# SCOPE (hard afgedwongen, zie .hoij_check_fit):
#   - 1 groep, complete data (listwise), estimator ML
#   - geen categorische variabelen, geen gelijkheidsrestricties
#   - geen multilevel / sampling weights
#
# Afhankelijkheden: lavaan (lavScores, lavTech, lavInspect), MASS (ginv-
# fallback). De rekenkern gebruikt drie niet-geexporteerde lavaan-
# functies (lav_model_x2glist, lav_model_implied, lav_model_gradient);
# zie ONDERZOEK_HOIJ_lavaan.md §2 voor het migratiepad naar de
# geexporteerde API.
# ============================================================


# ─────────────────────────────────────────────────────────────
# Interne rekenkern (identiek aan `kernfuncties`, hernoemd .hoij_*)
# ─────────────────────────────────────────────────────────────

## Versie-robuuste resolver voor niet-geexporteerde lavaan-internals:
## lavaan hernoemde rond 0.6-18 o.a. lav_model_x2GLIST -> lav_model_x2glist.
.hoij_lav_internal <- function(...) {
  for (nm in c(...)) {
    f <- tryCatch(get(nm, envir = asNamespace("lavaan")),
                  error = function(e) NULL)
    if (is.function(f)) return(f)
  }
  stop("lavaan-internal niet gevonden (geprobeerd: ",
       paste(c(...), collapse = ", "), "); deze lavaan-versie (",
       as.character(utils::packageVersion("lavaan")),
       ") wordt niet ondersteund.", call. = FALSE)
}

.hoij_x2glist  <- function() .hoij_lav_internal("lav_model_x2glist",
                                                "lav_model_x2GLIST")
.hoij_implied  <- function() .hoij_lav_internal("lav_model_implied")
.hoij_gradient <- function() .hoij_lav_internal("lav_model_gradient")

.hoij_loglik_casewise <- function(fit, theta) {
  X <- fit@Data@X[[1]]
  N <- nrow(X); p <- ncol(X)
  GLIST <- .hoij_x2glist()(fit@Model, x = theta)
  implied <- .hoij_implied()(fit@Model, GLIST = GLIST)
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

.hoij_all_J <- function(fit, theta0, delta = 1e-5) {
  D <- length(theta0); N <- nrow(fit@Data@X[[1]])
  J_array <- array(0, dim = c(N, D, D))
  ll_0 <- .hoij_loglik_casewise(fit, theta0)
  for (k in 1:D) {
    for (l in k:D) {
      if (k == l) {
        tp <- theta0; tp[k] <- tp[k] + delta
        tm <- theta0; tm[k] <- tm[k] - delta
        ll_p <- .hoij_loglik_casewise(fit, tp)
        ll_m <- .hoij_loglik_casewise(fit, tm)
        J_array[, k, k] <- -((ll_p - 2*ll_0 + ll_m) / (delta^2))
      } else {
        tpp <- theta0; tpp[k] <- tpp[k]+delta; tpp[l] <- tpp[l]+delta
        tpm <- theta0; tpm[k] <- tpm[k]+delta; tpm[l] <- tpm[l]-delta
        tmp_ <- theta0; tmp_[k] <- tmp_[k]-delta; tmp_[l] <- tmp_[l]+delta
        tmm <- theta0; tmm[k] <- tmm[k]-delta; tmm[l] <- tmm[l]-delta
        J_array[, k, l] <- -((.hoij_loglik_casewise(fit, tpp) -
                                .hoij_loglik_casewise(fit, tpm) -
                                .hoij_loglik_casewise(fit, tmp_) +
                                .hoij_loglik_casewise(fit, tmm)) / (4*delta^2))
        J_array[, l, k] <- J_array[, k, l]
      }
    }
  }
  J_array
}

.hoij_grad_F <- function(fit) {
  lavmodel       <- fit@Model
  lavsamplestats <- fit@SampleStats
  lavdata        <- fit@Data
  lavcache       <- fit@Cache
  function(theta) {
    GLIST <- .hoij_x2glist()(lavmodel, x = theta)
    as.numeric(.hoij_gradient()(
      lavmodel       = lavmodel,
      GLIST          = GLIST,
      lavsamplestats = lavsamplestats,
      lavdata        = lavdata,
      lavcache       = lavcache))
  }
}

.hoij_calibrate_alpha <- function(grad_F, theta0, H_observed, h = 1e-5) {
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

.hoij_T_tensor <- function(grad_F, theta, alpha, h = 1e-4) {
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


# ─────────────────────────────────────────────────────────────
# Scope-validatie: geen stille degradatie, informatieve fouten
# ─────────────────────────────────────────────────────────────

.hoij_check_fit <- function(fit) {
  if (!inherits(fit, "lavaan"))
    stop("'fit' moet een gefit lavaan-object zijn.", call. = FALSE)
  if (!isTRUE(lavaan::lavInspect(fit, "converged")))
    stop("Het lavaan-model is niet geconvergeerd.", call. = FALSE)
  if (lavaan::lavInspect(fit, "ngroups") != 1L)
    stop("hoij_lavaan() ondersteunt vooralsnog alleen 1 groep.", call. = FALSE)
  if (lavaan::lavInspect(fit, "nlevels") > 1L)
    stop("hoij_lavaan() ondersteunt geen multilevel-modellen.", call. = FALSE)
  opt <- lavaan::lavInspect(fit, "options")
  if (!opt$estimator %in% "ML")
    stop("hoij_lavaan() vereist estimator = \"ML\" (nu: ",
         opt$estimator, ").", call. = FALSE)
  if (opt$missing %in% c("ml", "fiml", "ml.x", "two.stage", "robust.two.stage"))
    stop("hoij_lavaan() ondersteunt geen missing-data-methoden (missing = \"",
         opt$missing, "\"); gebruik complete data.", call. = FALSE)
  if (isTRUE(lavaan::lavInspect(fit, "categorical")))
    stop("hoij_lavaan() ondersteunt geen categorische (ordinale) variabelen.",
         call. = FALSE)
  if (fit@Model@eq.constraints ||
      (!is.null(fit@Model@ceq.function) &&
       !identical(body(fit@Model@ceq.function), quote(NULL)) &&
       length(fit@Model@ceq.linear.idx) + length(fit@Model@ceq.nonlinear.idx) > 0))
    stop("hoij_lavaan() ondersteunt geen gelijkheidsrestricties; ",
         "gedefinieerde parameters kunnen wel via 'functional'.", call. = FALSE)
  invisible(TRUE)
}


# ─────────────────────────────────────────────────────────────
# Functionalen-interface
#   NULL           -> alle vrije parameters
#   character      -> expressies in parameternamen, bv. "a*b" of
#                     "1 - `speed~~speed` / `visual~~visual`"
#                     (niet-syntactische namen tussen backticks)
#   function       -> phi(theta) met theta = benoemde vrije-parametervector
#   (benoemde) lijst van functies en/of expressies mag ook
# ─────────────────────────────────────────────────────────────

.hoij_make_functionals <- function(functional, th_names) {
  as_fun <- function(x, label) {
    if (is.function(x)) return(x)
    if (is.character(x) && length(x) == 1L) {
      ## backtick-vrije parameternamen die geen geldige R-namen zijn
      ## (bv. speed~~speed) kunnen in de expressie met backticks worden
      ## aangeduid; evaluatie gebeurt in een omgeving met alle namen.
      expr <- parse(text = x)[[1]]
      return(function(theta) eval(expr, envir = as.list(theta)))
    }
    stop("Functionaal '", label, "' moet een functie of een ",
         "character-expressie zijn.", call. = FALSE)
  }
  if (is.null(functional)) {
    fns <- lapply(th_names, function(nm) {
      force(nm); function(theta) unname(theta[nm])
    })
    names(fns) <- th_names
    return(fns)
  }
  if (is.function(functional)) functional <- list(functional)
  if (is.character(functional)) functional <- as.list(functional)
  if (!is.list(functional))
    stop("'functional' moet NULL, een character-vector, een functie of ",
         "een lijst daarvan zijn.", call. = FALSE)
  nms <- names(functional)
  if (is.null(nms)) nms <- rep("", length(functional))
  auto <- vapply(seq_along(functional), function(i) {
    if (nzchar(nms[i])) nms[i]
    else if (is.character(functional[[i]])) functional[[i]]
    else sprintf("phi%d", i)
  }, character(1))
  fns <- lapply(seq_along(functional), function(i)
    as_fun(functional[[i]], auto[i]))
  names(fns) <- auto
  fns
}


# ─────────────────────────────────────────────────────────────
# Hoofdfunctie
# ─────────────────────────────────────────────────────────────

#' Higher-Order Infinitesimal Jackknife SE's en CI's voor lavaan
#'
#' @param fit        geconvergeerd lavaan-object (1 groep, ML, complete data)
#' @param functional NULL (alle vrije parameters), character-expressie(s) in
#'                   parameternamen (bv. "a*b"), functie(s) phi(theta), of
#'                   een (benoemde) lijst daarvan
#' @param B          aantal multinomiale gewichtsvectoren (pseudo-replicaten)
#' @param order      1 = IJ1 (lineair), 2 = HOIJ-2 (default)
#' @param kappa      trust-region-demping van de tweede-orde stap
#' @param level      betrouwbaarheidsniveau van het percentielinterval
#' @param admissibility "keep" (default: alle replicaten tellen mee) of
#'                   "drop" (replicaten met negatieve variantieparameters
#'                   worden voor SE/CI geschrapt; sensitiviteitsvariant)
#' @param alpha_spread_tol tolerantie op de alpha-kalibratiespread
#' @param seed       optionele seed voor de gewichtstrekking
#' @param details    TRUE: bewaar ook theta-replicaten en gewichtsmatrix W
#'                   (voor benchmarking tegen een bootstrap met dezelfde W)
#' @return object van klasse "hoij_lavaan" met $results (est, se, lo, hi),
#'         $diagnostics en optioneel $replicates/$weights
hoij_lavaan <- function(fit,
                        functional = NULL,
                        B = 1000L,
                        order = 2L,
                        kappa = 0.5,
                        level = 0.95,
                        admissibility = c("keep", "drop"),
                        alpha_spread_tol = 0.1,
                        seed = NULL,
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

  ## indices variantieparameters (toelaatbaarheid replicaten)
  spl <- strsplit(th_names, "~~", fixed = TRUE)
  var_idx <- which(vapply(spl, function(z)
    length(z) == 2 && z[1] == z[2], logical(1)))

  ## ── eenmalige setup ──
  t0 <- proc.time()[["elapsed"]]
  Scores <- lavaan::lavScores(fit, scaling = TRUE)              # N x D
  H.inv  <- lavaan::lavTech(fit, "inverted.information.observed")
  if (is.null(Scores) || is.null(H.inv))
    stop("Scores of geobserveerde informatie niet beschikbaar voor dit model.",
         call. = FALSE)

  alpha <- NA_real_; alpha_spread <- NA_real_
  Tmat <- NULL; J_2d <- NULL
  if (order == 2L) {
    H_obs  <- lavaan::lavTech(fit, "information.observed")
    grad_F <- .hoij_grad_F(fit)
    cal <- .hoij_calibrate_alpha(grad_F, theta0, H_obs)
    alpha <- cal$alpha; alpha_spread <- cal$spread
    if (!is.finite(alpha) || cal$spread > alpha_spread_tol)
      stop(sprintf(paste0("Alpha-kalibratie inconsistent (spread = %.3g > %.3g): ",
                          "de tweede-orde stap is voor dit model niet ",
                          "betrouwbaar berekenbaar. Gebruik order = 1 (IJ1) ",
                          "of een bootstrap."), cal$spread, alpha_spread_tol),
           call. = FALSE)
    T_arr <- .hoij_T_tensor(grad_F, theta0, alpha)
    J_all <- .hoij_all_J(fit, theta0)
    Tmat  <- matrix(T_arr, nrow = D)                            # D x D^2
    J_2d  <- matrix(J_all, nrow = N, ncol = D * D)              # N x D^2
  }
  t_setup <- proc.time()[["elapsed"]] - t0

  ## ── gewichten en replicatielus (gevectoriseerd waar mogelijk) ──
  t0 <- proc.time()[["elapsed"]]
  W  <- t(rmultinom(B, size = N, prob = rep(1, N)))             # B x N
  dW <- W - 1L
  C_mat <- (dW %*% Scores) %*% H.inv                            # B x D

  theta_rep <- sweep(-C_mat, 2, theta0, "+")                    # IJ1
  s_vec <- rep(1, B)
  if (order == 2L) {
    JW <- (dW %*% J_2d) / N                                     # B x D^2
    HT <- H.inv %*% Tmat                                        # D x D^2
    for (i in seq_len(B)) {
      c_vec <- C_mat[i, ]
      J_dw  <- matrix(JW[i, ], D, D)
      Bc    <- drop(H.inv %*% (J_dw %*% c_vec))
      Ac    <- 0.5 * drop(HT %*% as.vector(tcrossprod(c_vec)))
      d2    <- Bc - Ac
      n1 <- sqrt(sum(c_vec^2)); n2 <- sqrt(sum(d2^2))
      s  <- if (n2 > 0) min(1, kappa * n1 / n2) else 1
      s_vec[i] <- s
      theta_rep[i, ] <- theta0 - c_vec + s * d2
    }
  }
  colnames(theta_rep) <- th_names
  inadmiss <- if (length(var_idx))
    apply(theta_rep[, var_idx, drop = FALSE] < 0, 1, any)
  else rep(FALSE, B)
  t_rep <- proc.time()[["elapsed"]] - t0

  ## ── SE + percentiel-CI per functionaal ──
  keep <- if (admissibility == "drop") !inadmiss else rep(TRUE, B)
  if (sum(keep) < max(40L, ceiling(0.5 * B)))
    warning("Minder dan de helft van de replicaten toelaatbaar; ",
            "SE/CI zijn mogelijk onbetrouwbaar.", call. = FALSE)
  pr <- c((1 - level) / 2, 1 - (1 - level) / 2)

  res <- do.call(rbind, lapply(names(fns), function(nm) {
    phi  <- fns[[nm]]
    est  <- tryCatch(as.numeric(phi(theta0)), error = function(e) NA_real_)
    vals <- vapply(seq_len(B), function(i)
      tryCatch(as.numeric(phi(theta_rep[i, ])),
               error = function(e) NA_real_), numeric(1))
    v <- vals[keep & is.finite(vals)]
    if (length(v) >= 40L) {
      q  <- quantile(v, pr, names = FALSE, type = 7)
      se <- sd(v)
    } else { q <- c(NA_real_, NA_real_); se <- NA_real_ }
    data.frame(functional = nm, est = est, se = se,
               lo = q[1], hi = q[2],
               n_used = length(v), stringsAsFactors = FALSE)
  }))
  rownames(res) <- NULL

  out <- list(
    results     = res,
    diagnostics = list(
      order = order, B = B, kappa = kappa, level = level,
      admissibility = admissibility,
      alpha = alpha, alpha_spread = alpha_spread,
      frac_damped = if (order == 2L) mean(s_vec < 1) else NA_real_,
      mean_s      = if (order == 2L) mean(s_vec) else NA_real_,
      frac_inadmissible = mean(inadmiss),
      time_setup_s = t_setup, time_replicates_s = t_rep,
      N = N, D = D),
    call = match.call())
  if (details) { out$replicates <- theta_rep; out$weights <- W }
  class(out) <- "hoij_lavaan"
  out
}


# ─────────────────────────────────────────────────────────────
# Print-methode
# ─────────────────────────────────────────────────────────────

print.hoij_lavaan <- function(x, digits = 3, ...) {
  d <- x$diagnostics
  cat(sprintf("%s (B = %d gewichtsvectoren, %d%% percentiel-CI)\n",
              if (d$order == 2L) "HOIJ-2 (higher-order infinitesimal jackknife)"
              else "IJ1 (eerste-orde infinitesimal jackknife)",
              d$B, round(100 * d$level)))
  cat(sprintf("N = %d, D = %d vrije parameters | setup %.2fs + replicaten %.2fs\n",
              d$N, d$D, d$time_setup_s, d$time_replicates_s))
  if (d$order == 2L)
    cat(sprintf("alpha = %.2f (spread %.2e) | gedempt: %.1f%% (mean s = %.3f)\n",
                d$alpha, d$alpha_spread, 100 * d$frac_damped, d$mean_s))
  if (d$frac_inadmissible > 0)
    cat(sprintf("Niet-toelaatbare replicaten (negatieve variantie): %.1f%% (%s)\n",
                100 * d$frac_inadmissible,
                if (d$admissibility == "keep") "meegenomen" else "geschrapt"))
  cat("\n")
  r <- x$results
  r[c("est", "se", "lo", "hi")] <- lapply(r[c("est", "se", "lo", "hi")],
                                          round, digits = digits)
  print(r, row.names = FALSE)
  invisible(x)
}
