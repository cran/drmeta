# drmeta 0.2.2

* Documentation update and archival release. No changes to package code.


# drmeta 0.2.1

## Behavior changes

* `drmeta_bootstrap_gamma()` now short-circuits when the constrained estimate
  already lies at `gamma = 0`. In that situation the null and alternative fits
  coincide, so the likelihood-ratio statistic is zero by construction and the
  simulated null distribution has a large point mass at zero. Bootstrapping
  ties at zero returns a value that reads like a p-value but only records the
  proportion of replicates that also reached the boundary. The function now
  returns `statistic = 0`, `p.value = 1`, and `B_used = 0` without running any
  replicates.

* The returned object gains a logical `boundary` component, `TRUE` when the
  short-circuit applied and `FALSE` otherwise. `print()` explains in the
  boundary case why no bootstrap was run and points to `constrained = FALSE`
  as the appropriate directional follow-up.

## New features

* `drmeta_bootstrap_gamma()` gains a `tol` argument controlling the threshold
  below which the observed likelihood-ratio statistic is treated as exactly
  zero. Default `1e-8`.

## Documentation

* The `Description` field cites Self and Liang (1987) for the boundary problem
  and Viechtbauer and Lopez-Lopez (2022) for the general location-scale parent
  model, and states that a scale model reweights studies rather than adjusting
  the mean for design-linked bias.

# drmeta 0.2.0

## Breaking changes

* The package is repositioned as a constrained specialization of the
  meta-analytic location-scale model rather than a standalone
  "variance-function framework for causal credibility". Documentation no
  longer describes the design index as a measure of causal credibility.

* `dr_heterogeneity()` no longer returns the design-residual variance
  decomposition or `R2_DR`. That decomposition was not identified by the
  model and is superseded by `dr_scale_attenuation()`, which is explicitly a
  description of the fitted scale function and not a proportion of
  heterogeneity explained.

* The following functions from 0.1.0 have been removed: `dr_variance()`,
  `dr_weights()`, `dr_forest()`, `dr_funnel()`, `dr_plot()`, and
  `dr_pub_bias()`. `dr_variance()` is superseded by `dr_scale_predict()`;
  study weights are available as the `weights` component of a fit. The
  publication-bias and plotting functions were removed pending revision
  against the corrected specification.

* `drmeta()` gains `mods`, `constrained`, `gamma_max`, and `gamma_fixed`
  arguments. The pooled effect is now returned as `beta`, a named vector of
  location coefficients, rather than a scalar `mu`.

## New features

* Location moderators and the joint design-indexed location-scale model.

* Exact estimation at the boundary `gamma = 0`, via direct optimization over
  a finite interval plus an explicit boundary comparison.

* Unrestricted scale fits (`constrained = FALSE`) as a directional diagnostic;
  a negative estimated gradient contradicts the substantive constraint rather
  than merely failing to support it.

* `drmeta_bootstrap_gamma()`, a parametric-bootstrap test of the boundary null.

* `dr_scale_predict()` and `dr_scale_attenuation()` for observed-support
  summaries of the fitted scale function.

* `dr_loo()` for leave-one-out influence on both the location and scale
  components.

* `dr_plot_vfun()` for the fitted variance function, drawn over observed
  support by default.

* `confint()` and `logLik()` methods. `logLik()` reports the criterion
  actually used and refuses to relabel a REML fit as ML.

* Warnings when the design index has fewer than three distinct values, when
  its observed range is narrow, and when the estimated gradient sits at the
  optimization bound.

## Bug fixes

* The bundled BCG example data used a `DR` column alongside `dr_weight`, so
  the documented `bcg$dr` silently partial-matched onto `dr_weight` and fed
  fitted weights into the model as the design index. The column is now named
  `dr`, the stale weight columns are removed, and examples index with `[[`.

* Roxygen blocks in the diagnostics file were separated from their functions
  by a blank line, so regenerating the documentation would have dropped
  `dr_scale_predict()`, `dr_scale_attenuation()`, and `dr_heterogeneity()`
  from the namespace.

* All exported functions now document a return value.

* The default `gamma_max` is 8, matching the optimization range reported in
  the accompanying manuscript. It was 20.
