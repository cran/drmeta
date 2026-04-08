# =============================================================================
# test-drmeta.R  --  testthat tests for the drmeta package
#
# Run with: devtools::test()  or  testthat::test_file("tests/testthat/test-drmeta.R")
# =============================================================================

library(testthat)
library(drmeta)

# Shared fixture: 20-study dataset used across tests
.make_data <- function(k = 20, seed = 42) {
  set.seed(seed)
  dr        <- runif(k, 0.1, 0.9)
  vi        <- runif(k, 0.01, 0.05)
  tau2_true <- 0.04 * exp(-2 * dr)
  yi        <- rnorm(k, 0.3, sqrt(vi + tau2_true))
  list(yi = yi, vi = vi, dr = dr)
}

d <- .make_data()


# =============================================================================
# 1. drmeta() -- core estimation
# =============================================================================

test_that("drmeta() returns a drmeta object with expected components", {
  fit <- drmeta(d$yi, d$vi, d$dr)

  expect_s3_class(fit, "drmeta")

  # Required list components
  expected <- c("mu", "se", "ci.lb", "ci.ub", "zval", "pval",
                "tau0sq", "gamma", "tau2_i", "sigma2_i", "weights",
                "loglik", "reml_loglik", "AIC", "BIC",
                "k", "yi", "vi", "dr", "slab", "vfun", "method",
                "converged", "optim_out", "call")
  expect_true(all(expected %in% names(fit)))
})

test_that("drmeta() pooled estimate is a finite scalar", {
  fit <- drmeta(d$yi, d$vi, d$dr)

  expect_length(fit$mu, 1)
  expect_true(is.finite(fit$mu))
  expect_true(fit$se > 0)
  expect_true(fit$ci.lb < fit$mu)
  expect_true(fit$mu  < fit$ci.ub)
})

test_that("drmeta() variance-function parameters are non-negative", {
  fit <- drmeta(d$yi, d$vi, d$dr)

  expect_gte(fit$tau0sq, 0)
  expect_gte(fit$gamma,  0)
  expect_true(all(fit$tau2_i   >= 0))
  expect_true(all(fit$sigma2_i >  0))
})

test_that("drmeta() weights sum to 1", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  expect_equal(sum(fit$weights), 1, tolerance = 1e-8)
})

test_that("drmeta() works with linear variance function", {
  fit <- drmeta(d$yi, d$vi, d$dr, vfun = "linear")
  expect_s3_class(fit, "drmeta")
  expect_equal(fit$vfun, "linear")
  expect_true(is.finite(fit$mu))
})

test_that("drmeta() works with ML estimation", {
  fit <- drmeta(d$yi, d$vi, d$dr, method = "ML")
  expect_s3_class(fit, "drmeta")
  expect_equal(fit$method, "ML")
  expect_true(is.finite(fit$mu))
})

test_that("drmeta() AIC < BIC is not guaranteed, but both are finite", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  expect_true(is.finite(fit$AIC))
  expect_true(is.finite(fit$BIC))
})

test_that("drmeta() accepts custom study labels", {
  labs <- paste0("S", seq_along(d$yi))
  fit  <- drmeta(d$yi, d$vi, d$dr, slab = labs)
  expect_equal(fit$slab, labs)
})

test_that("drmeta() warns and defaults when dr is NULL", {
  expect_warning(
    fit <- drmeta(d$yi, d$vi, dr = NULL),
    regexp = "DR_i = 0.5"
  )
  expect_true(all(fit$dr == 0.5))
})


# =============================================================================
# 2. S3 methods
# =============================================================================

test_that("print.drmeta() runs without error", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  expect_output(print(fit), regexp = "DR-Meta")
})

test_that("summary.drmeta() runs without error and prints all sections", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  out <- capture.output(summary(fit))
  expect_true(any(grepl("Pooled Effect", out)))
  expect_true(any(grepl("Variance-Function", out)))
  expect_true(any(grepl("Model Fit", out)))
})

test_that("coef.drmeta() returns named vector of length 3", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  co  <- coef(fit)
  expect_length(co, 3)
  expect_named(co, c("mu", "tau0sq", "gamma"))
})

test_that("confint.drmeta() returns a data frame with correct columns", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  ci  <- confint(fit)
  expect_s3_class(ci, "data.frame")
  expect_true(all(c("lower", "upper", "estimate") %in% names(ci)))
  expect_true(ci$lower < ci$estimate)
  expect_true(ci$estimate < ci$upper)
})

test_that("fitted.drmeta() returns vector of length k all equal to mu", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  f   <- fitted(fit)
  expect_length(f, fit$k)
  expect_true(all(f == fit$mu))
})

test_that("residuals.drmeta() raw and standardised work", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  r   <- residuals(fit)
  rs  <- residuals(fit, type = "standardised")
  expect_length(r,  fit$k)
  expect_length(rs, fit$k)
  expect_false(isTRUE(all.equal(r, rs)))   # they differ
})

test_that("logLik.drmeta() returns logLik-class object", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  ll  <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_true(is.finite(as.numeric(ll)))
})

test_that("AIC() and BIC() dispatch correctly", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  expect_true(is.finite(AIC(fit)))
  expect_true(is.finite(BIC(fit)))
})


# =============================================================================
# 3. Edge cases
# =============================================================================

test_that("drmeta() errors with fewer than 3 studies", {
  expect_error(
    drmeta(c(0.1, 0.2), c(0.01, 0.02), c(0.5, 0.8)),
    regexp = "At least 3"
  )
})

test_that("drmeta() errors with non-positive vi", {
  yi <- d$yi[1:5]; vi <- d$vi[1:5]; dr <- d$dr[1:5]
  vi[2] <- 0
  expect_error(drmeta(yi, vi, dr), regexp = "strictly positive")
})

test_that("drmeta() errors with dr outside [0,1]", {
  yi <- d$yi[1:5]; vi <- d$vi[1:5]; dr <- d$dr[1:5]
  dr[1] <- -0.1
  expect_error(drmeta(yi, vi, dr), regexp = "\\[0, 1\\]")
})

test_that("drmeta() handles NA removal with a warning", {
  yi <- d$yi; vi <- d$vi; dr <- d$dr
  yi[3] <- NA
  expect_warning(fit <- drmeta(yi, vi, dr), regexp = "missing")
  expect_equal(fit$k, length(d$yi) - 1L)
})

test_that("drmeta() works with k = 3 (minimum)", {
  set.seed(1)
  # k=3 may not converge -- that is acceptable; we just check it runs
  suppressWarnings(
    fit <- drmeta(rnorm(3, 0.3), runif(3, 0.01, 0.05), c(0.2, 0.5, 0.9))
  )
  expect_s3_class(fit, "drmeta")
  expect_true(is.finite(fit$mu))
})

test_that("drmeta() with constant DR reduces to finite result (Proposition 1)", {
  set.seed(7)
  k  <- 15
  yi <- rnorm(k, 0.3, 0.15)
  vi <- runif(k, 0.01, 0.04)
  dr <- rep(0.5, k)   # constant DR
  fit <- drmeta(yi, vi, dr)
  expect_s3_class(fit, "drmeta")
  expect_true(is.finite(fit$mu))
  # All tau2_i should be identical when DR is constant
  expect_true(diff(range(fit$tau2_i)) < 1e-10)
})

test_that("drmeta() convergence flag is logical", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  expect_type(fit$converged, "logical")
})


# =============================================================================
# 4. dr_score() and dr_from_design()
# =============================================================================

test_that("dr_score() returns values in [0,1]", {
  k       <- 10
  balance <- runif(k)
  overlap <- runif(k)
  dr      <- dr_score(balance = balance, overlap = overlap)
  expect_length(dr, k)
  expect_true(all(dr >= 0 & dr <= 1))
})

test_that("dr_score() equal weights gives simple mean", {
  b <- c(0.2, 0.8)
  o <- c(0.6, 0.4)
  dr <- dr_score(balance = b, overlap = o)
  expect_equal(as.numeric(dr), (b + o) / 2, tolerance = 1e-10)
})

test_that("dr_score() respects custom weights", {
  b  <- c(1, 0)
  o  <- c(0, 1)
  dr <- dr_score(balance = b, overlap = o, weights = c(1, 0))
  expect_equal(as.numeric(dr), b, tolerance = 1e-10)
})

test_that("dr_score() clips values outside [0,1] with warning", {
  expect_warning(
    dr <- dr_score(x = c(-0.1, 0.5, 1.2)),
    regexp = "clipped"
  )
  expect_true(all(dr >= 0 & dr <= 1))
})

test_that("dr_score() errors with zero sub-scores", {
  expect_error(dr_score(), regexp = "at least one")
})

test_that("dr_score() carries subscores attribute", {
  dr <- dr_score(a = c(0.5, 0.7), b = c(0.3, 0.9))
  expect_true(!is.null(attr(dr, "subscores")))
  expect_s3_class(attr(dr, "subscores"), "data.frame")
})

test_that("dr_from_design() maps known labels correctly", {
  scores <- dr_from_design(c("rct", "ols", "iv", "did"))
  expect_equal(scores, c(1.00, 0.25, 0.75, 0.60))
})

test_that("dr_from_design() is case-insensitive", {
  expect_equal(
    dr_from_design(c("RCT", "OLS")),
    dr_from_design(c("rct", "ols"))
  )
})

test_that("dr_from_design() warns on unknown labels", {
  expect_warning(dr_from_design("unknown_xyz"), regexp = "Unrecognised")
})

test_that("dr_from_design() uses default_score for unknown labels", {
  suppressWarnings(
    s <- dr_from_design("unknown_xyz", default_score = 0.33)
  )
  expect_equal(s, 0.33)
})

test_that("dr_from_design() respects custom_map", {
  s <- dr_from_design("mydesign", custom_map = c(mydesign = 0.65),
                      warn_unknown = FALSE)
  expect_equal(s, 0.65)
})

test_that("normalize_01() rescales to [0,1]", {
  x  <- c(2, 5, 8)
  nx <- normalize_01(x)
  expect_equal(nx, c(0, 0.5, 1), tolerance = 1e-10)
})

test_that("normalize_01() returns zeros for constant input", {
  expect_equal(normalize_01(c(3, 3, 3)), c(0, 0, 0))
})


# =============================================================================
# 5. Heterogeneity diagnostics
# =============================================================================

test_that("dr_heterogeneity() returns list with summary, decomposition, contributions", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  het <- dr_heterogeneity(fit)

  expect_type(het, "list")
  expect_named(het, c("summary", "decomposition", "contributions"))
  expect_s3_class(het$summary,       "data.frame")
  expect_s3_class(het$decomposition, "data.frame")
  expect_s3_class(het$contributions, "data.frame")
})

test_that("dr_heterogeneity() summary has expected columns", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  het <- dr_heterogeneity(fit)
  expect_true(all(c("k", "Q", "df", "pval", "I2", "H2") %in% names(het$summary)))
})

test_that("dr_heterogeneity() I2 is in [0, 100]", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  I2  <- dr_heterogeneity(fit)$summary$I2
  expect_gte(I2, 0)
  expect_lte(I2, 100)
})

test_that("dr_heterogeneity() pct_Q contributions sum to ~100", {
  fit   <- drmeta(d$yi, d$vi, d$dr)
  pct   <- dr_heterogeneity(fit)$contributions$pct_Q
  expect_equal(sum(pct), 100, tolerance = 0.01)
})

test_that("dr_heterogeneity() errors on non-drmeta input", {
  expect_error(dr_heterogeneity(list(a = 1)), regexp = "drmeta")
})


# =============================================================================
# 6. Leave-one-out diagnostics
# =============================================================================

test_that("dr_loo() returns list with summary and full", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  loo <- dr_loo(fit)

  expect_type(loo, "list")
  expect_named(loo, c("summary", "full"))
  expect_s3_class(loo$summary, "data.frame")
  expect_s3_class(loo$full,    "drmeta")
})

test_that("dr_loo() summary has k rows (one per study)", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  loo <- dr_loo(fit)
  expect_equal(nrow(loo$summary), fit$k)
})

test_that("dr_loo() delta_mu is finite for all studies", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  loo <- dr_loo(fit)
  expect_true(all(is.finite(loo$summary$delta_mu)))
})

test_that("dr_loo() influential column is logical", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  loo <- dr_loo(fit)
  expect_type(loo$summary$influential, "logical")
})

test_that("dr_loo() errors on non-drmeta input", {
  expect_error(dr_loo(list()), regexp = "drmeta")
})


# =============================================================================
# 7. Publication bias
# =============================================================================

test_that("dr_pub_bias() returns list with PET, PEESE, Egger, recommendation", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  pb  <- dr_pub_bias(fit)

  expect_type(pb, "list")
  expect_true(all(c("PET", "PEESE", "Egger", "recommendation") %in% names(pb)))
})

test_that("dr_pub_bias() PET is a data frame with expected rows", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  pb  <- dr_pub_bias(fit)
  expect_s3_class(pb$PET, "data.frame")
  expect_true("intercept" %in% rownames(pb$PET))
  expect_true("sei"       %in% rownames(pb$PET))
})

test_that("dr_pub_bias() recommendation is a character string", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  pb  <- dr_pub_bias(fit)
  expect_type(pb$recommendation, "character")
  expect_gt(nchar(pb$recommendation), 0)
})

test_that("dr_pub_bias() can run subset of tests", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  pb  <- dr_pub_bias(fit, test = "PET")
  expect_true("PET"  %in% names(pb))
  expect_false("PEESE" %in% names(pb))
})

test_that("dr_pub_bias() errors on non-drmeta input", {
  expect_error(dr_pub_bias(42), regexp = "drmeta")
})


# =============================================================================
# 8. Plot functions (smoke tests -- just check they run without error)
# =============================================================================

test_that("dr_forest() runs without error", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  pdf(NULL)
  expect_no_error(dr_forest(fit))
  dev.off()
})

test_that("dr_forest() returns invisible data frame", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  pdf(NULL)
  result <- dr_forest(fit)
  dev.off()
  expect_s3_class(result, "data.frame")
})

test_that("dr_plot() runs without error", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  pdf(NULL)
  expect_no_error(dr_plot(fit))
  dev.off()
})

test_that("dr_plot_vfun() runs without error", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  pdf(NULL)
  expect_no_error(dr_plot_vfun(fit))
  dev.off()
})

test_that("dr_funnel() runs without error", {
  fit <- drmeta(d$yi, d$vi, d$dr)
  pdf(NULL)
  expect_no_error(dr_funnel(fit))
  dev.off()
})


# =============================================================================
# 9. Weight and variance utilities
# =============================================================================

test_that("dr_weights() returns positive weights of correct length", {
  vi <- runif(10, 0.01, 0.05)
  dr <- runif(10)
  w  <- dr_weights(vi, dr, tau0sq = 0.04, gamma = 1.5)
  expect_length(w, 10)
  expect_true(all(w > 0))
})

test_that("dr_weights() normalise = 'sum1' gives weights summing to 1", {
  vi <- runif(10, 0.01, 0.05)
  dr <- runif(10)
  w  <- dr_weights(vi, dr, tau0sq = 0.04, gamma = 1.5, normalise = "sum1")
  expect_equal(sum(w), 1, tolerance = 1e-10)
})

test_that("dr_weights() normalise = 'pct' gives weights summing to 100", {
  vi <- runif(10, 0.01, 0.05)
  dr <- runif(10)
  w  <- dr_weights(vi, dr, tau0sq = 0.04, gamma = 1.5, normalise = "pct")
  expect_equal(sum(w), 100, tolerance = 1e-8)
})

test_that("dr_variance() returns non-negative values", {
  grid <- seq(0, 1, by = 0.1)
  tau2 <- dr_variance(grid, tau0sq = 0.04, gamma = 1.5)
  expect_true(all(tau2 >= 0))
  expect_length(tau2, length(grid))
})

test_that("dr_variance() is monotone decreasing for exponential vfun", {
  grid <- seq(0, 1, by = 0.05)
  tau2 <- dr_variance(grid, tau0sq = 0.04, gamma = 1.5)
  expect_true(all(diff(tau2) <= 0))
})

test_that("dr_variance() accepts a fitted drmeta object", {
  fit  <- drmeta(d$yi, d$vi, d$dr)
  tau2 <- dr_variance(fit)
  expect_length(tau2, fit$k)
  expect_true(all(tau2 >= 0))
})
