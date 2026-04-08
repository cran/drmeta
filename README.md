# drmeta <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/causalfragility-lab/drmeta/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/causalfragility-lab/drmeta/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

## Design-Robust Meta-Analysis: A Variance-Function Framework for Causal Credibility

**drmeta** implements the DR-Meta model of Hait (2026) — a
variance-function random-effects framework that embeds causal design
robustness directly into the heterogeneity structure of meta-analysis.

### The core idea

Standard random-effects meta-analysis absorbs design-quality differences into
a single constant $\tau^2$. This allows large-sample but weakly-identified
studies to dominate pooled estimates. DR-Meta instead models between-study
variance as a monotone-decreasing function of a design robustness index
$\text{DR}_i \in [0,1]$:

$$y_i \sim N\!\bigl(\mu,\; \underbrace{v_i}_{\text{sampling}} +
\underbrace{\tau^2(\text{DR}_i;\,\psi)}_{\text{design-dependent}}\bigr)$$

Study weights $w_i = 1/\sigma_i^2$ therefore automatically down-weight weaker
designs — not via ad-hoc quality multipliers, but through a coherent
likelihood-based mechanism.

### Special cases

| Condition | DR-Meta reduces to |
|---|---|
| $\text{DR}_i = c$ (constant) | Classical random-effects (Proposition 1) |
| $\tau_0^2 = 0$ | Fixed-effects meta-analysis (Proposition 2) |

## Installation
```r
# Install the development version from GitHub:
devtools::install_github("causalfragility-lab/drmeta")
```

## Quick start
```r
library(drmeta)

# 1. Build the DR index
dr <- dr_from_design(c("RCT", "DiD", "OLS", "IV", "matching"))

# or from continuous sub-scores:
dr <- dr_score(
  balance = c(0.9, 0.6, 0.3, 0.8, 0.5),
  overlap  = c(0.85, 0.7, 0.4, 0.75, 0.6),
  design   = dr_from_design(c("RCT", "DiD", "OLS", "IV", "matching"))
)

# 2. Fit DR-Meta
fit <- drmeta(yi = my_effects, vi = my_variances, dr = dr)
summary(fit)

# 3. Visualise
dr_forest(fit)       # forest plot ordered by DR
dr_plot(fit)         # weights vs DR (Lemma 3)
dr_plot_vfun(fit)    # estimated variance function

# 4. Diagnostics
dr_heterogeneity(fit)  # Q, I², Proposition 6 decomposition
dr_loo(fit)            # leave-one-out influence

# 5. Publication bias
dr_pub_bias(fit)  # PET / PEESE / Egger with DR-Meta weights
dr_funnel(fit)    # funnel plot
```

## Key functions

| Function | Description |
|---|---|
| `drmeta()` | Fit DR-Meta via profiled REML (core estimation) |
| `dr_score()` | Build DR index from continuous sub-scores |
| `dr_from_design()` | Build DR index from design type labels |
| `dr_weights()` | Compute and inspect DR-Meta weights |
| `dr_variance()` | Evaluate the fitted variance function |
| `dr_heterogeneity()` | Cochran's Q, I², H², Proposition 6 decomposition |
| `dr_loo()` | Leave-one-out influence diagnostics |
| `dr_forest()` | Forest plot ordered by DR |
| `dr_plot()` | Weight vs DR diagnostic plot |
| `dr_plot_vfun()` | Variance function curve |
| `dr_pub_bias()` | PET / PEESE / Egger publication bias |
| `dr_funnel()` | Funnel plot with DR-coloured points |
| `normalize_01()` | Rescale any vector to [0,1] |

## Connection to the literature

DR-Meta is a special case of the **meta-analytic location-scale model**
(Viechtbauer & Lopez-Lopez, 2022) in which the scale component is explicitly
motivated by causal design theory. It is complementary to bias-model
frameworks (Turner et al., 2009; Rhodes et al., 2015; Mathur & VanderWeele,
2020), which adjust the *mean* for confounding while DR-Meta adjusts the
*variance*.

## References

- Hait, S. (2026). *Design-Robust Meta-Analysis: A Variance-Function
  Framework for Causal Credibility*. Unpublished manuscript.
- Viechtbauer, W., & Lopez-Lopez, J. A. (2022). *Research Synthesis Methods,
  13*, 697–715. <https://doi.org/10.1002/jrsm.1562>
- Turner, R. M., et al. (2009). *JRSS-A, 172*, 21–47.
- Rhodes, K. M., et al. (2015). *JRSS-A, 178*, 641–666.
- Mathur, M. B., & VanderWeele, T. J. (2020). *JASA, 115*, 1190–1204.
