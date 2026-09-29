# Independent numerical validation of the drmeta likelihood and diagnostics
# Initially audited under R 4.3.3 and rerun through the 0.2.3 finalization gate
# on 2026-09-28 with metafor 4.4-0 and metadat 1.2-0.
# Mapping: metafor scale = ~ dr, log link  =>  alpha[1] = log(tau0sq), alpha[2] = -gamma.
# metafor bounds must be finite on BOTH sides of a constrained coefficient, so the
# directional fit uses alpha[2] in [-8, 0], matching drmeta's default gamma_max = 8.

suppressMessages({library(metafor); library(drmeta)})

dat <- metadat::dat.landenberger2005
lab <- tolower(trimws(as.character(dat$design)))
dat$dr <- c(nonequiv = 0, match = 0.5, rct = 1)[lab]
stopifnot(!anyNA(dat$dr))
dat <- escalc("OR", ai = n.cbt.non, bi = n.cbt.rec,
              ci = n.ctrl.non, di = n.ctrl.rec, data = dat)

bnd <- list(alpha.min = c(-Inf, -8), alpha.max = c(Inf, 0))
mf  <- function(m) c(beta0 = unname(m$beta[1]), tau0sq = exp(m$alpha[1]),
                     gamma = -m$alpha[2], logLik = as.numeric(logLik(m)))
dm  <- function(m) c(beta0 = unname(m$beta[1]), tau0sq = m$tau0sq,
                     gamma = m$gamma, logLik = as.numeric(logLik(m)))

compare <- function(d, method) {
  suppressWarnings(rbind(
    metafor_unrestricted = mf(rma(yi, vi, scale = ~ dr, data = d, method = method)),
    drmeta_unrestricted  = dm(drmeta(d$yi, d$vi, d$dr, constrained = FALSE, method = method)),
    metafor_bounded      = mf(rma(yi, vi, scale = ~ dr, data = d, method = method, control = bnd)),
    drmeta_directional   = dm(drmeta(d$yi, d$vi, d$dr, method = method))))
}

# 1. Full data. ML log-likelihoods agree exactly; REML criteria differ by an additive
#    constant (metafor includes the log|X'X| term), so compare REML point estimates only.
print(compare(dat, "ML"),   digits = 6)
print(compare(dat, "REML"), digits = 6)

# 2. Joint model at the gamma = 0 boundary must equal ordinary meta-regression.
jd <- drmeta(dat$yi, dat$vi, dat$dr, mods = 1 - dat$dr)
jr <- rma(yi, vi, mods = ~ I(1 - dr), data = dat)
stopifnot(jd$gamma == 0,
          isTRUE(all.equal(unname(jd$beta), unname(c(coef(jr))), tolerance = 1e-5)),
          isTRUE(all.equal(jd$tau0sq, jr$tau2, tolerance = 1e-4)))

# 3. Without Anderson (2002): interior optimum. metafor's bounded fit from default
#    start values stops at alpha[2] = 0 with a lower likelihood; supplying a start
#    value recovers the drmeta solution.
d2 <- dat[dat$study != "Anderson (2002)", ]
print(compare(d2, "ML"), digits = 6)
fixed_start <- rma(yi, vi, scale = ~ dr, data = d2, method = "ML",
                   control = c(bnd, list(alpha.init = c(log(0.04), -0.4))))
print(mf(fixed_start), digits = 6)

# 4. Is heterogeneity design-dependent at all? A categorical scale model does not
#    impose monotonicity. Full data: matched groups carry most residual variance.
cat_full <- rma(yi, vi, scale = ~ factor(design), data = dat, method = "ML")
cat_null <- rma(yi, vi, scale = ~ 1, data = dat, method = "ML")
print(anova(cat_full, cat_null))
cat_d2 <- rma(yi, vi, scale = ~ factor(design), data = d2, method = "ML")
print(anova(cat_d2, rma(yi, vi, scale = ~ 1, data = d2, method = "ML")))

# 5. Group-specific fitted variances in score order.
design_order <- factor(lab, levels = c("nonequiv", "match", "rct"))
cat_groups <- rma(yi, vi, scale = ~ design_order,
                  data = dat, method = "ML")
design_grid <- data.frame(
  design_order = factor(levels(design_order), levels = levels(design_order))
)
scale_matrix <- model.matrix(~ design_order, data = design_grid)
tau2_by_design <- setNames(
  exp(as.vector(scale_matrix %*% as.numeric(cat_groups$alpha))),
  levels(design_order)
)
print(tau2_by_design, digits = 6)

# 6. Applied diagnostics reported in the paper and vignette.
directional_reml <- drmeta(dat$yi, dat$vi, dat$dr, slab = dat$study)
loo <- dr_loo(directional_reml)
print(table(positive_gamma = loo$gamma_loo > 0, useNA = "ifany"))
randomized_vs_rest <- as.numeric(dat$dr == 1)
binary_fit <- drmeta(dat$yi, dat$vi, randomized_vs_rest, slab = dat$study)
print(c(randomized_vs_rest_gamma = binary_fit$gamma))

without_anderson_reml <- drmeta(d2$yi, d2$vi, d2$dr, slab = d2$study)
deletion_boot <- drmeta_bootstrap_gamma(without_anderson_reml, B = 999,
                                        seed = 20260928)
print(c(LR = deletion_boot$statistic, p = deletion_boot$p.value))

# 7. Numerical invariants used by the package regression tests.
full_mf <- rma(yi, vi, scale = ~ dr, data = dat, method = "ML")
full_dm <- drmeta(dat$yi, dat$vi, dat$dr, constrained = FALSE, method = "ML")
stopifnot(
  abs(unname(full_dm$beta[1] - full_mf$beta[1])) < 1e-5,
  abs(unname(full_dm$tau0sq - exp(full_mf$alpha[1]))) < 1e-5,
  abs(unname(full_dm$gamma + full_mf$alpha[2])) < 1e-5,
  abs(as.numeric(logLik(full_dm)) - as.numeric(logLik(full_mf))) < 1e-5,
  abs(unname(without_anderson_reml$gamma) - 0.457) < 2e-3,
  sum(loo$gamma_loo > 0, na.rm = TRUE) == 23L,
  abs(unname(binary_fit$gamma) - 1.78) < .02
)
