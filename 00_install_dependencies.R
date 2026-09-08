# =====================================================================
# 00_install_dependencies.R
#
# Run this script once before 01_worked_example.R,
# 02_timing_comparison.R and 03_simulation_study.R.
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

## Check the joint mean/covariance derivative conventions before analysis.
source("hoij_core.R")
if (!hoij_selftest(meanstructure = TRUE))
  stop("Joint mean/covariance derivative self-test failed.")
