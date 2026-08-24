# drmeta

**Design-Indexed Location-Scale Meta-Analysis**

`drmeta` fits meta-analytic location-scale models in which residual
between-study heterogeneity is an exponential function of a prespecified
design-robustness index:

```
y_i ~ N( x_i'beta ,  v_i + tau0^2 * exp(-gamma * DR_i) )
```

The constrained form imposes `gamma >= 0`, encoding the directional hypothesis
that unexplained heterogeneity does not increase as design robustness
improves. Setting `gamma = 0` recovers the conventional random-effects model
exactly.

## What this model does and does not do

A scale model changes how precision is allocated across studies. It does
**not** identify or remove a systematic design-linked shift in the conditional
mean. If the mean varies with design robustness, the scale-only pooled
estimate is a model-dependent weighted average, not the effect of a
hypothetical perfectly designed study. Supply location moderators via `mods`
when design-linked mean differences are plausible, and do not describe scale
reweighting as confounding adjustment.

## Installation

```r
# install.packages("remotes")
remotes::install_github("causalfragility-lab/drmeta")
```

## Quick start

```r
library(drmeta)

bcg <- utils::read.csv(
  system.file("extdata", "bcg_design_robustness.csv", package = "drmeta")
)

# Constrained fit
fit <- drmeta(yi = bcg[["yi"]], vi = bcg[["vi"]], dr = bcg[["dr"]])
summary(fit)

# Minimum comparison set
re  <- drmeta(bcg[["yi"]], bcg[["vi"]], bcg[["dr"]], gamma_fixed = 0)  # random effects
ls  <- drmeta(bcg[["yi"]], bcg[["vi"]], bcg[["dr"]], constrained = FALSE)
jls <- drmeta(bcg[["yi"]], bcg[["vi"]], bcg[["dr"]], mods = 1 - bcg[["dr"]])

# Interpretation over observed support, not extrapolated to [0,1]
dr_scale_attenuation(fit)
dr_plot_vfun(fit)

# Boundary-aware inference and influence
drmeta_bootstrap_gamma(fit, B = 999, seed = 1)
dr_loo(fit)
```

## Reporting

Report `gamma` together with the observed range of the design index, the
convergence and boundary status, and the fitted variance ratio or attenuation
over a prespecified contrast **within** observed support. The magnitude of
`gamma` alone depends on the scaling of the index and is not an invariant
measure of design sensitivity.

For a primary analysis, construct the design index without using the realized
effect estimate, its standard error, its p-value, its confidence interval, or
any outcome-dependent diagnostic. Scores that use such information define an
exploratory analysis.

## Status

Current version 0.2.2. `NEWS.md` documents the breaking changes from
0.1.0, including the functions removed and the withdrawal of the
design-explained variance decomposition, and the 0.2.1 change to
`drmeta_bootstrap_gamma()`, which now short-circuits when the constrained
scale gradient is estimated at the boundary.

## License

MIT
