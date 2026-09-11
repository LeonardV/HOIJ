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
| `02_timing_comparison.R` | Subsection 3.3: timing decomposition for the mediation (D = 30) and bifactor (D = 36) models, including nine free intercepts each. |
| `03_simulation_study.R` | Section 4. |
| `hoij_lavaan.R` | `hoij_lavaan()`, a reusable function that returns HOIJ-2 standard errors and percentile intervals for your own fitted lavaan model. |

Scripts assume the repository root as the working directory and write
their output to a script-specific subdirectory ending in `_means`.
These directories separate the joint mean/covariance runs from earlier
outputs and simulation checkpoints. Both output paths in the simulation
script, including the figure-regeneration section, use the new directory.

```r
source("00_install_dependencies.R")   # once
source("01_worked_example.R")
```

`03_simulation_study.R` is the expensive one (12 cells x 1,000 data sets
x 1,001 fits). Run it with `SMOKE_TEST <- TRUE` first.

## Using HOIJ-2 on your own model

```r
library(lavaan)

source("hoij_core.R") 
source("hoij_lavaan.R")

model <- '
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
  visual ~ c*textual + b*speed
  speed  ~ a*textual
'
fit <- sem(model, data = HolzingerSwineford1939, estimator = "ML",
           meanstructure = TRUE)

hoij_lavaan(fit, functional = c(ab = "a*b"), B = 5000, order = 2, seed = 42)

HOIJ-2 (second-order infinitesimal jackknife)
B = 5000 weight vectors | 95% percentile CI | N = 301, D = 30
setup 5.39s + replicates 4.11s
derivative check 9.7e-06
inadmissible replicates: 1.3% (kept)

 functional   est    se    lo    hi n_used
         ab 0.095 0.053 0.017 0.228   5000
```

`functional` accepts expressions in the free-parameter names, functions
of the parameter vector, or `NULL` for every free parameter. The current
scope is single-group ML with complete data and continuous indicators,
without equality constraints.
