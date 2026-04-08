## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 7,
  fig.height = 5
)

## ----dr-from-design-----------------------------------------------------------
library(drmeta)

designs <- c("RCT", "IV", "DiD", "matching", "OLS", "OLS", "RCT", "DiD")
dr <- dr_from_design(designs)
data.frame(Design = designs, DR = dr)

## ----dr-from-subscores--------------------------------------------------------
k <- 8
balance  <- c(0.95, 0.80, 0.60, 0.70, 0.30, 0.25, 0.90, 0.65)
overlap  <- c(0.90, 0.75, 0.55, 0.80, 0.40, 0.35, 0.85, 0.70)
design_s <- dr_from_design(designs)   # reuse labels from above

dr_composite <- dr_score(
  balance  = balance,
  overlap  = overlap,
  design   = design_s,
  weights  = c(2, 2, 3)   # design quality weighted most
)
round(dr_composite, 3)

## ----fit-drmeta---------------------------------------------------------------
set.seed(42)
k  <- 20
dr <- runif(k, 0.1, 0.95)
vi <- runif(k, 0.008, 0.055)

# True data-generating model (Sections 3–4, Hait 2025):
#   tau0sq = 0.04, gamma = 2, mu = 0.30
tau2_true <- 0.04 * exp(-2 * dr)
yi <- rnorm(k, mean = 0.30, sd = sqrt(vi + tau2_true))

# Fit DR-Meta (exponential variance function, REML)
fit <- drmeta(yi = yi, vi = vi, dr = dr)
print(fit)

## ----summary-drmeta-----------------------------------------------------------
summary(fit)

## ----compare-vfun-------------------------------------------------------------
fit_lin <- drmeta(yi, vi, dr, vfun = "linear")
cat("Exponential AIC:", round(AIC(fit), 2), "\n")
cat("Linear      AIC:", round(AIC(fit_lin), 2), "\n")

## ----forest-plot, fig.height=7------------------------------------------------
dr_forest(fit, main = "DR-Meta Forest Plot — Simulated Data")

## ----vfun-plot----------------------------------------------------------------
dr_plot_vfun(fit, main = "Estimated Variance Function")

## ----weight-plot--------------------------------------------------------------
dr_plot(fit, main = "DR-Meta Weights vs Design Robustness")

## ----heterogeneity------------------------------------------------------------
het <- dr_heterogeneity(fit)

cat("--- Heterogeneity Summary ---\n")
print(het$summary)

cat("\n--- Proposition 6 Decomposition ---\n")
print(het$decomposition)
cat(sprintf("\n%.1f%% of total heterogeneity is design-explained (R2_DR).\n",
            het$decomposition$R2_DR * 100))

cat("\n--- Per-Study Q Contributions (top 5) ---\n")
top5 <- head(het$contributions[order(-het$contributions$pct_Q), ], 5)
print(top5)

## ----loo----------------------------------------------------------------------
loo <- dr_loo(fit)
cat("Influential studies:\n")
print(subset(loo$summary, influential == TRUE))
cat("\nFull LOO summary (first 5 rows):\n")
print(head(loo$summary, 5))

## ----pub-bias-----------------------------------------------------------------
pb <- dr_pub_bias(fit)
cat("--- PET ---\n");   print(pb$PET)
cat("\n--- PEESE ---\n"); print(pb$PEESE)
cat("\n--- Recommendation ---\n"); cat(pb$recommendation, "\n")

## ----funnel-------------------------------------------------------------------
dr_funnel(fit, main = "DR-Meta Funnel Plot")

