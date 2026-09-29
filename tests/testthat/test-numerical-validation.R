test_that("ML unrestricted estimates match metafor location-scale estimates", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("metadat")

  dat <- metadat::dat.landenberger2005
  lab <- tolower(trimws(as.character(dat$design)))
  dat$dr <- c(nonequiv = 0, match = .5, rct = 1)[lab]
  expect_false(anyNA(dat$dr))
  dat <- metafor::escalc("OR", ai = n.cbt.non, bi = n.cbt.rec,
                         ci = n.ctrl.non, di = n.ctrl.rec, data = dat)

  mf <- metafor::rma(yi, vi, scale = ~ dr, data = dat, method = "ML")
  dm <- drmeta(dat$yi, dat$vi, dat$dr, method = "ML",
               constrained = FALSE, slab = dat$study)
  expect_lt(abs(unname(dm$beta[1] - mf$beta[1])), 1e-5)
  expect_lt(abs(unname(dm$tau0sq - exp(mf$alpha[1]))), 1e-5)
  expect_lt(abs(unname(dm$gamma + mf$alpha[2])), 1e-5)
  expect_lt(abs(as.numeric(logLik(dm)) - as.numeric(logLik(mf))), 1e-5)
})

test_that("the gamma zero joint fit matches ordinary meta-regression", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("metadat")

  dat <- metadat::dat.landenberger2005
  lab <- tolower(trimws(as.character(dat$design)))
  dat$dr <- c(nonequiv = 0, match = .5, rct = 1)[lab]
  dat <- metafor::escalc("OR", ai = n.cbt.non, bi = n.cbt.rec,
                         ci = n.ctrl.non, di = n.ctrl.rec, data = dat)
  dm <- drmeta(dat$yi, dat$vi, dat$dr, mods = 1 - dat$dr,
               method = "REML", slab = dat$study)
  mf <- metafor::rma(yi, vi, mods = ~ I(1 - dr), data = dat,
                     method = "REML")
  expect_equal(dm$gamma, 0)
  expect_lt(max(abs(unname(dm$beta) - unname(c(stats::coef(mf))))), 1e-5)
  expect_lt(abs(dm$tau0sq - mf$tau2), 1e-4)
})

test_that("informative starts recover the Anderson-deletion interior optimum", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("metadat")

  dat <- metadat::dat.landenberger2005
  lab <- tolower(trimws(as.character(dat$design)))
  dat$dr <- c(nonequiv = 0, match = .5, rct = 1)[lab]
  dat <- metafor::escalc("OR", ai = n.cbt.non, bi = n.cbt.rec,
                         ci = n.ctrl.non, di = n.ctrl.rec, data = dat)
  dat <- dat[dat$study != "Anderson (2002)", ]
  bounds <- list(alpha.min = c(-Inf, -8), alpha.max = c(Inf, 0),
                 alpha.init = c(log(.04), -.4))
  mf <- suppressWarnings(metafor::rma(yi, vi, scale = ~ dr, data = dat,
                                      method = "ML", control = bounds))
  dm <- drmeta(dat$yi, dat$vi, dat$dr, method = "ML", slab = dat$study)
  expect_lt(abs(unname(dm$gamma + mf$alpha[2])), 1e-4)
  expect_lt(abs(unname(dm$tau0sq - exp(mf$alpha[1]))), 1e-5)
  expect_lt(abs(as.numeric(logLik(dm)) - as.numeric(logLik(mf))), 1e-5)
})

test_that("grouped shape check detects the full-data interior peak", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("metadat")

  dat <- metadat::dat.landenberger2005
  lab <- tolower(trimws(as.character(dat$design)))
  dat$dr <- c(nonequiv = 0, match = .5, rct = 1)[lab]
  dat <- metafor::escalc("OR", ai = n.cbt.non, bi = n.cbt.rec,
                         ci = n.ctrl.non, di = n.ctrl.rec, data = dat)
  fit <- drmeta(dat$yi, dat$vi, dat$dr, method = "ML", slab = dat$study)
  groups <- factor(lab, levels = c("nonequiv", "match", "rct"))
  chk <- dr_shape_check(fit, groups)
  expect_identical(chk$pattern, "interior peak")
  expect_lt(abs(chk$statistic - 5.9915), 1e-3)
  expect_lt(abs(chk$p.value - .05), 2e-3)
  expect_lt(abs(chk$estimates$tau2[2] - .157), 5e-3)
})
