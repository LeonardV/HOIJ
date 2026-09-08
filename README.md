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
| `test_hoij_lavaan.R` | Tests for `hoij_lavaan()`, including a comparison against an exact bootstrap on the same weight vectors. |
| `test_hoij_centering.R` | Checks exact weighted means/covariances, agreement of joint and profiled routes, and the simulation bootstrap helper against raw-data refits. |

Scripts assume the repository root as the working directory and write
their output to a script-specific subdirectory ending in `_means`.
These directories separate the joint mean/covariance runs from earlier
outputs and simulation checkpoints. Both output paths in the simulation
script, including the figure-regeneration section, use the new directory.

```r
source("00_install_dependencies.R")   # once
source("test_hoij_centering.R")
source("test_hoij_lavaan.R")
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

```

`functional` accepts expressions in the free-parameter names, functions
of the parameter vector, or `NULL` for every free parameter. The current
scope is single-group ML with complete data and continuous indicators,
without equality constraints.

## Implementation notes

All analysis and bootstrap fits now use `meanstructure = TRUE`. The nine
observed intercepts are free and the latent means remain fixed at zero.
Scores, observed information, casewise curvature and third derivatives
are computed for the full parameter vector, including those intercepts.
This accounts for changing sample means within the joint HOIJ-2 expansion.
The simulated populations retain zero means, as in the original design.

The simulation bootstrap supplies both `sample.mean` and `sample.cov`.
Bootstrap covariances use divisor N - 1 and `sample.cov.rescale = TRUE`,
so lavaan converts them to the normal-theory ML covariance with divisor N.
Population covariance matrices instead use `sample.cov.rescale = FALSE`.
The Wald and Monte Carlo comparisons use the same full parameter vector.

Rerun the worked example, timing comparison and simulation before updating
the manuscript. Earlier numerical examples and the existing supplementary
PDFs have not been regenerated for this route. The extra parameters change
the derivative setup cost, so earlier timings do not describe these scripts.

The third derivatives are second differences of lavaan's analytic
gradient, taken without any rescaling: for normal-theory ML the function
lavaan optimises is the negative mean log-likelihood up to a constant,
which self-test (b) verifies. `hoij_lavaan()` therefore refuses
`likelihood = "wishart"`, which rescales that objective by roughly
N/(N-1). `check_gradient_hessian()` confirms that the finite differences
reproduce lavaan's analytic observed information before they are used.

For compatibility, covariance-only ML (`meanstructure = FALSE`) remains
supported by the reusable function. This route
profiles out the observed means. Recentring a reweighted sample changes
its ML covariance by a term that is quadratic in the mean shift, which
the fixed-centre casewise derivatives do not capture. `hoij2_replicates()`
therefore takes the fitted object as its `fit` argument and adds this
profiled-mean correction; with an explicit mean structure the joint
derivatives already contain it and no correction is applied.
`test_hoij_centering.R` checks both cases against the exact weighted
covariance and verifies the estimated means on the joint route.
