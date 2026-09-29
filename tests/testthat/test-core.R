test_that("gamma zero nests random effects", {
  yi <- c(-.2, .1, .3, .5, -.1, .2)
  vi <- rep(.04, length(yi)); dr <- seq(0, 1, length.out = length(yi))
  f <- drmeta(yi, vi, dr, gamma_fixed = 0)
  expect_equal(f$gamma, 0)
  expect_true(is.finite(f$tau0sq))
  expect_true(all(abs(diff(f$tau2i)) < 1e-12))
})

test_that("higher DR gives higher weight when vi is equal", {
  yi <- c(0, .1, .2, .3, .4, .5); vi <- rep(.03, 6); dr <- seq(0, 1, length.out = 6)
  f <- drmeta(yi, vi, dr, gamma_fixed = 2)
  expect_true(all(diff(f$weights) > 0))
})

test_that("weights do not exceed fixed-effect weights", {
  yi <- c(0, .1, .2, .3, .4, .5); vi <- seq(.02, .07, length.out = 6)
  dr <- seq(0, 1, length.out = 6)
  f <- drmeta(yi, vi, dr, gamma_fixed = 2)
  expect_true(all(1 / f$sigma2 <= 1 / vi + 1e-12))
})

test_that("attenuation is bounded for constrained fits", {
  yi <- c(0, .1, .2, .3, .4, .5); vi <- rep(.03, 6); dr <- seq(0, 1, length.out = 6)
  f <- drmeta(yi, vi, dr, gamma_fixed = 2)
  a <- dr_scale_attenuation(f)
  expect_gte(a, 0); expect_lte(a, 1)
})

test_that("the constrained fit can attain the exact boundary", {
  set.seed(11)
  k <- 40; dr <- runif(k); vi <- runif(k, .005, .03)
  yi <- rnorm(k, .3, sqrt(vi + .10 * exp(3 * dr)))   # wrong-direction gradient
  expect_equal(drmeta(yi, vi, dr, .quiet = TRUE)$gamma, 0)
  expect_lt(drmeta(yi, vi, dr, constrained = FALSE, .quiet = TRUE)$gamma, 0)
})

test_that("the bundled example data uses an unambiguous dr column", {
  path <- system.file("extdata", "bcg_design_robustness.csv", package = "drmeta")
  bcg <- utils::read.csv(path)
  expect_true("dr" %in% names(bcg))
  # Guards against partial matching of bcg$dr onto some other dr-prefixed column.
  expect_identical(bcg$dr, bcg[["dr"]])
  expect_true(all(bcg[["dr"]] >= 0 & bcg[["dr"]] <= 1))
})

test_that("location moderators are estimated and returned", {
  set.seed(3)
  k <- 40; dr <- runif(k); vi <- runif(k, .005, .03)
  yi <- rnorm(k, .3 + .2 * (1 - dr), sqrt(vi + .10 * exp(-3 * dr)))
  f <- drmeta(yi, vi, dr, mods = 1 - dr, .quiet = TRUE)
  expect_length(coef(f), 2L)
  expect_equal(dim(confint(f)), c(2L, 2L))
})

test_that("logLik carries df and nobs and refuses relabelling", {
  set.seed(5)
  k <- 30; dr <- runif(k); vi <- runif(k, .005, .03)
  yi <- rnorm(k, .3, sqrt(vi + .05 * exp(-2 * dr)))
  f <- drmeta(yi, vi, dr, method = "REML", .quiet = TRUE)
  ll <- logLik(f)
  expect_equal(attr(ll, "df"), 3L)
  expect_equal(attr(ll, "nobs"), k - 1L)
  expect_error(logLik(f, REML = FALSE), "refit")
})

test_that("leave-one-out returns one row per study", {
  set.seed(7)
  k <- 20; dr <- runif(k); vi <- runif(k, .005, .03)
  yi <- rnorm(k, .3, sqrt(vi + .10 * exp(-3 * dr)))
  labels <- paste("Study", seq_len(k))
  f <- drmeta(yi, vi, dr, slab = labels, .quiet = TRUE)
  loo <- dr_loo(f)
  expect_equal(nrow(loo), k)
  expect_identical(loo$study, labels)
  expect_true(all(c("delta_est", "delta_gamma", "influential") %in% names(loo)))
  override <- paste("Trial", seq_len(k))
  expect_identical(dr_loo(f, slab = override)$study, override)
  expect_error(dr_loo(f, slab = override[-1]), "one non-missing label")
})

test_that("dr_score normalises weights and clips to the unit interval", {
  s <- suppressWarnings(dr_score(a = c(0, .5, 1.4), b = c(1, .5, 0), weights = c(3, 1)))
  expect_true(all(s >= 0 & s <= 1))
  expect_equal(unname(s[3]), .75 * 1 + .25 * 0)
  expect_error(dr_score(a = c(0, .5, 1), b = c(1, .5, 0),
                        weights = c(0, 0)), "nonnegative")
})

test_that("input validation rejects malformed data", {
  expect_error(drmeta(1:5, rep(-1, 5), seq(0, 1, length.out = 5)), "positive")
  expect_error(drmeta(1:5, rep(.1, 5), seq(0, 2, length.out = 5)), "\\[0,1\\]")
  expect_error(drmeta(1:2, rep(.1, 2), c(0, 1)), "at least three")
  expect_error(drmeta(1:5, rep(.1, 5), seq(0, 1, length.out = 5),
                      slab = letters[1:4]), "one non-missing label")
})

test_that("bootstrap short-circuits at the boundary", {
  path <- system.file("extdata", "bcg_design_robustness.csv", package = "drmeta")
  bcg <- utils::read.csv(path)
  fit <- drmeta(bcg[["yi"]], bcg[["vi"]], bcg[["dr"]], .quiet = TRUE)
  expect_equal(fit$gamma, 0)
  b <- drmeta_bootstrap_gamma(fit, B = 99, seed = 1)
  expect_true(b$boundary)
  expect_equal(b$statistic, 0)
  expect_equal(b$p.value, 1)
  expect_equal(b$B_used, 0L)
})

test_that("bootstrap runs for a clearly positive scale gradient", {
  dr <- rep(c(0, .5, 1), each = 12)
  amp <- rep(c(.8, .25, .04), each = 12)
  yi <- .3 + rep(c(-1, 1), 18) * amp
  vi <- rep(.002, length(yi))
  fit <- suppressWarnings(drmeta(yi, vi, dr, method = "ML", .quiet = TRUE))
  expect_gt(fit$gamma, 0)
  b <- suppressWarnings(drmeta_bootstrap_gamma(fit, B = 9, seed = 9))
  expect_false(b$boundary)
  expect_gt(b$B_used, 0L)
  expect_true(b$p.value > 0 && b$p.value <= 1)
})
