# =====================================================================
# 00_install_dependencies.R
#
# Companion code for:
#   Vanbrabant, L., & Rosseel, Y. Approximating percentile bootstrap
#   confidence intervals in SEM without repeated refitting: A tutorial
#   on the second-order infinitesimal jackknife.
#
# Run this script once before 01_worked_example.R,
# 02_timing_comparison.R and 03_simulation_study.R.
#
# WHY THE DEVELOPMENT VERSION OF lavaan?
# The casewise curvature J_i and the third-derivative array Khat are not
# part of lavaan's public API. They are obtained by re-evaluating the
# model-implied moments at perturbed parameter values, which uses three
# internal helpers. Those helpers were renamed between the CRAN release
# and the development version:
#
#   CRAN (<= 0.6.17)            development
#   ------------------------    ------------------------------
#   lav_model_x2GLIST()         lav_model_x2glist()
#   lav_model_gradient()        lav_model_grad()
#   lav_model_implied(GLIST=)   lav_model_implied(glist=)
#
# The argument rename is the dangerous one: lav_model_implied() accepts
# `...`, so a call using the wrong spelling is silently ignored and the
# implied moments are returned at the *fitted* parameters. The casewise
# log-likelihood is then constant in theta and J_i collapses to zero
# without any error being raised.
#
# hoij_core.R resolves both spellings at run time, so the scripts work
# with either version, and hoij_selftest() verifies the scaling
# conventions before anything is computed. The development version is
# nevertheless what the reported results were produced with, and is what
# we recommend installing.
# =====================================================================

## 1. Development version of lavaan --------------------------------------
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}

remotes::install_github("yrosseel/lavaan", upgrade = "never")

## 2. Remaining packages -------------------------------------------------
##    MASS            generalised inverse fallback (all scripts)
##    covsim          VITA/copula non-normal data generation (script 03)
##    rvinecopulib    vine copula objects returned by covsim::vita()
pkgs <- c("MASS", "covsim", "rvinecopulib")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1),
                             quietly = TRUE)]
if (length(missing_pkgs)) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

## 3. Check --------------------------------------------------------------
suppressPackageStartupMessages(library(lavaan))
cat("lavaan", as.character(packageVersion("lavaan")), "installed\n")

internals_ok <- any(vapply(
  c("lav_model_x2glist", "lav_model_x2GLIST"),
  function(nm) exists(nm, envir = asNamespace("lavaan")), logical(1)))
if (!internals_ok) {
  stop("This lavaan version exposes neither lav_model_x2glist() nor ",
       "lav_model_x2GLIST(); hoij_core.R needs to be updated.")
}
cat("lavaan internals required by hoij_core.R are available\n")
