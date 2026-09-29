## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.width = 6, fig.height = 4)
library(drmeta)
library(metafor)
library(metadat)

## ----data---------------------------------------------------------------------
dat <- metadat::dat.landenberger2005
design <- tolower(trimws(as.character(dat$design)))
dat$dr <- c(nonequiv = 0, match = .5, rct = 1)[design]
stopifnot(!anyNA(dat$dr))
dat <- metafor::escalc("OR", ai = n.cbt.non, bi = n.cbt.rec,
                       ci = n.ctrl.non, di = n.ctrl.rec, data = dat)
table(dat$dr)

## ----four-fits----------------------------------------------------------------
constant <- drmeta(dat$yi, dat$vi, dat$dr, gamma_fixed = 0,
                   slab = dat$study)
directional <- drmeta(dat$yi, dat$vi, dat$dr, slab = dat$study)
unrestricted <- drmeta(dat$yi, dat$vi, dat$dr, constrained = FALSE,
                       slab = dat$study)
joint <- drmeta(dat$yi, dat$vi, dat$dr, mods = 1 - dat$dr,
                slab = dat$study)

data.frame(
  model = c("constant", "directional", "unrestricted", "joint"),
  beta0 = vapply(list(constant, directional, unrestricted, joint),
                 function(x) unname(x$beta[1]), numeric(1)),
  tau0sq = vapply(list(constant, directional, unrestricted, joint),
                  function(x) x$tau0sq, numeric(1)),
  gamma = vapply(list(constant, directional, unrestricted, joint),
                 function(x) x$gamma, numeric(1))
)

## ----shape-check--------------------------------------------------------------
ordered_design <- factor(design,
                         levels = c("nonequiv", "match", "rct"))
shape <- dr_shape_check(directional, ordered_design)
shape

## ----boundary-----------------------------------------------------------------
drmeta_bootstrap_gamma(directional, B = 999, seed = 20260928)

## ----loo----------------------------------------------------------------------
loo <- dr_loo(directional)
anderson <- loo[loo$study == "Anderson (2002)", ]
anderson[, c("study", "est_loo", "tau0sq_loo", "gamma_loo")]
sum(loo$gamma_loo > 0, na.rm = TRUE)

## ----anderson-bootstrap-------------------------------------------------------
keep <- dat$study != "Anderson (2002)"
without_anderson <- drmeta(dat$yi[keep], dat$vi[keep], dat$dr[keep],
                           slab = dat$study[keep])
anderson_test <- drmeta_bootstrap_gamma(without_anderson, B = 999,
                                        seed = 20260928)
c(LR = anderson_test$statistic, p = anderson_test$p.value)

## ----score-sensitivity--------------------------------------------------------
randomized_vs_rest <- as.numeric(dat$dr == 1)
binary_fit <- drmeta(dat$yi, dat$vi, randomized_vs_rest,
                     slab = dat$study)
binary_fit$gamma

matched_at_075 <- ifelse(dat$dr == .5, .75, dat$dr)
spacing_fit <- drmeta(dat$yi, dat$vi, matched_at_075,
                      slab = dat$study)
c(primary = directional$gamma,
  matched_at_075 = spacing_fit$gamma,
  randomized_vs_rest = binary_fit$gamma)

## ----numerical-validation-----------------------------------------------------
mf <- metafor::rma(yi, vi, scale = ~ dr, data = dat, method = "ML")
dm <- drmeta(dat$yi, dat$vi, dat$dr, constrained = FALSE,
             method = "ML", slab = dat$study)
stopifnot(
  abs(unname(dm$beta[1] - mf$beta[1])) < 1e-5,
  abs(unname(dm$tau0sq - exp(mf$alpha[1]))) < 1e-5,
  abs(unname(dm$gamma + mf$alpha[2])) < 1e-5,
  abs(as.numeric(logLik(dm)) - as.numeric(logLik(mf))) < 1e-5
)

## ----session-info-------------------------------------------------------------
sessionInfo()

