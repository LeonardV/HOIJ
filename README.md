# HOIJ-2 for SEM: replication code

Replication code for:

> Vanbrabant, L., & Rosseel, Y. *Approximating percentile bootstrap
> confidence intervals in SEM without repeated refitting: The second-order infinitesimal jackknife.*

The second-order infinitesimal jackknife (HOIJ-2) approximates the
nonparametric bootstrap distribution of (functions of) the parameters of
a fitted SEM from a single model fit plus a one-time derivative setup,
instead of refitting the model for every bootstrap resample.

## Files

| File | What it does |
|---|---|
| `00_install_dependencies.R` | Installs the development version of **lavaan** and the remaining packages. Run once. |
| `hoij_core.R` | Computational kernel: casewise scores and curvature, third derivatives, the IJ1 and HOIJ-2 replicates. |
| `01_worked_example.R` | Section 3. |
| `02_timing_comparison.R` | Subsection 3.3: timing decomposition for the mediation (D = 21) and bifactor (D = 27) models. |
| `03_simulation_study.R` | Section 4. |
| `hoij_lavaan.R` | `hoij_lavaan()`, a reusable function that returns HOIJ-2 standard errors and percentile intervals for your own fitted lavaan model. |
| `test_hoij_lavaan.R` | Tests for `hoij_lavaan()`, including a comparison against an exact bootstrap on the same weight vectors. |

Scripts assume the repository root as the working directory and write
their output to a script-specific subdirectory.

```r
source("00_install_dependencies.R")   # once
source("01_worked_example.R")
```

`03_simulation_study.R` is the expensive one (12 cells x 1,000 data sets
x 1,001 fits). Run it with `SMOKE_TEST <- TRUE` first.

## Using HOIJ-2 on your own model

```r
source("hoij_core.R")
source("hoij_lavaan.R")

fit <- lavaan::sem(model, data = mydata, estimator = "ML")

hoij_lavaan(fit, functional = c(ab = "a*b", psi = "`speed~~speed`"),
            B = 1000, order = 2, seed = 1)
```

`functional` accepts expressions in the free-parameter names, functions
of the parameter vector, or `NULL` for every free parameter. The current
scope is single-group ML with complete data and continuous indicators,
without equality constraints.

## An implementation notes

The third derivatives are second differences of lavaan's analytic
gradient, taken without any rescaling: for normal-theory ML the function
lavaan optimises is the negative mean log-likelihood up to a constant,
which self-test (b) verifies. `hoij_lavaan()` therefore refuses
`likelihood = "wishart"`, which rescales that objective by roughly
N/(N-1). `check_gradient_hessian()` confirms that the finite differences
reproduce lavaan's analytic observed information before they are used.