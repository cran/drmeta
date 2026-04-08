# =============================================================================
# drmeta.R  —  Core estimation for the DR-Meta model
#
# Fits a random-effects meta-analysis in which between-study heterogeneity
# is a monotone-decreasing function of each study's design robustness index
# DR_i in [0,1].  The variance function is:
#
#   Exponential (default):  tau^2(DR_i) = tau0sq * exp(-gamma * DR_i)
#   Linear:                 tau^2(DR_i) = tau0sq * max(1 - gamma * DR_i, 0)
#
# Parameters (tau0sq >= 0, gamma >= 0) are estimated by REML or ML via
# numerical optimisation; the pooled estimate mu is then the
# inverse-total-variance weighted mean.
#
# References: Hait (2025), "Design-Robust Meta-Analysis: A Variance-Function
# Framework for Causal Credibility".
# =============================================================================


# =============================================================================
# Internal helpers
# =============================================================================

#' @keywords internal
.make_tau2_fn <- function(vfun) {
  if (vfun == "exponential") {
    function(dr, tau0sq, gamma) tau0sq * exp(-gamma * dr)
  } else {
    # Linear: tau0sq * max(1 - gamma * dr, 0)
    function(dr, tau0sq, gamma) tau0sq * pmax(1 - gamma * dr, 0)
  }
}

# Log-likelihood (and REML log-likelihood) for DR-Meta
# theta = c(log_tau0sq, log_gamma) to enforce positivity (or 0 on boundary)
.drmeta_loglik <- function(theta, yi, vi, dr, tau2_fn, drmeta_method) {
  tau0sq <- exp(theta[1])
  gamma  <- exp(theta[2])

  tau2_i   <- tau2_fn(dr, tau0sq, gamma)
  sigma2_i <- vi + tau2_i

  if (any(sigma2_i <= 0)) return(1e10)

  w    <- 1 / sigma2_i
  mu   <- sum(w * yi) / sum(w)
  resid <- yi - mu

  ll <- -0.5 * (sum(log(sigma2_i)) + sum(resid^2 / sigma2_i))

  if (drmeta_method == "REML") {
    ll <- ll - 0.5 * log(sum(w))
  }

  -ll  # return negative log-likelihood for minimisation
}


# =============================================================================
# Main estimation function
# =============================================================================

#' Fit a Design-Robust Meta-Analysis (DR-Meta) Model
#'
#' Fits a random-effects meta-analysis model in which between-study
#' heterogeneity is a monotone-decreasing function of each study's design
#' robustness index \eqn{\mathrm{DR}_i \in [0,1]}.  Studies with higher design
#' robustness receive less heterogeneity weight, implementing Proposition 1 of
#' Hait (2025).
#'
#' @param yi Numeric vector of \eqn{k} effect-size estimates.
#' @param vi Numeric vector of \eqn{k} sampling variances (all must be > 0).
#' @param dr Numeric vector of \eqn{k} design robustness indices in the
#'   interval \eqn{[0, 1]}.  If \code{NULL}, a warning is issued and all
#'   studies are assigned \eqn{\mathrm{DR}_i = 0.5}.
#' @param vfun Variance function: \code{"exponential"} (default) or
#'   \code{"linear"}.
#' @param method Estimation method: \code{"REML"} (default) or \code{"ML"}.
#' @param slab Optional character vector of study labels.
#' @param control List of control arguments passed to \code{stats::optim}.
#'
#' @return An object of class \code{"drmeta"} (a named list). Key components:
#'   \code{mu} (pooled estimate), \code{se}, \code{ci.lb}, \code{ci.ub},
#'   \code{zval}, \code{pval}, \code{tau0sq}, \code{gamma}, \code{tau2_i},
#'   \code{sigma2_i}, \code{weights}, \code{loglik}, \code{reml_loglik},
#'   \code{AIC}, \code{BIC}, \code{k}, \code{yi}, \code{vi}, \code{dr},
#'   \code{slab}, \code{vfun}, \code{method}, \code{converged},
#'   \code{optim_out}, \code{call}.
#'
#' @references
#' Hait, S. (2025). Design-Robust Meta-Analysis: A Variance-Function
#' Framework for Causal Credibility.
#'
#' @examples
#' set.seed(42)
#' k  <- 20
#' dr <- runif(k, 0.1, 0.9)
#' vi <- runif(k, 0.01, 0.05)
#' tau2_true <- 0.04 * exp(-2 * dr)
#' yi <- rnorm(k, 0.3, sqrt(vi + tau2_true))
#'
#' fit <- drmeta(yi, vi, dr)
#' print(fit)
#' summary(fit)
#'
#' @seealso \code{\link{dr_heterogeneity}}, \code{\link{dr_loo}},
#'   \code{\link{dr_pub_bias}}, \code{\link{dr_forest}},
#'   \code{\link{dr_score}}, \code{\link{dr_from_design}}
#' @export
drmeta <- function(yi, vi, dr = NULL,
                   vfun    = c("exponential", "linear"),
                   method  = c("REML", "ML"),
                   slab    = NULL,
                   control = list()) {

  cl     <- match.call()
  vfun   <- match.arg(vfun)
  method <- match.arg(method)

  # ---- DR default ----------------------------------------------------------
  if (is.null(dr)) {
    warning("No DR_i supplied; setting DR_i = 0.5 for all studies.",
            call. = FALSE)
    dr <- rep(0.5, length(yi))
  }

  # ---- Input validation ----------------------------------------------------
  k_raw <- length(yi)

  if (!is.numeric(yi))
    stop("`yi` must be a numeric vector.", call. = FALSE)
  if (!is.numeric(vi))
    stop("`vi` must be a numeric vector.", call. = FALSE)
  if (!is.numeric(dr))
    stop("`dr` must be a numeric vector.", call. = FALSE)
  if (length(vi) != k_raw || length(dr) != k_raw)
    stop("`yi`, `vi`, and `dr` must all have the same length.", call. = FALSE)
  if (any(vi <= 0, na.rm = TRUE))
    stop("`vi` must be strictly positive.", call. = FALSE)
  if (any(dr < 0 | dr > 1, na.rm = TRUE))
    stop("`dr` values must lie in [0, 1].", call. = FALSE)

  # ---- NA removal ----------------------------------------------------------
  keep <- complete.cases(yi, vi, dr)
  if (any(!keep)) {
    warning(sprintf(
      "%d stud%s with missing values removed.",
      sum(!keep), if (sum(!keep) == 1) "y" else "ies"
    ), call. = FALSE)
    yi <- yi[keep]; vi <- vi[keep]; dr <- dr[keep]
    if (!is.null(slab)) slab <- slab[keep]
  }

  k <- length(yi)
  if (k < 3L)
    stop("At least 3 studies are required.", call. = FALSE)

  # ---- Study labels --------------------------------------------------------
  if (is.null(slab)) slab <- paste0("Study", seq_len(k))

  # ---- Variance function ---------------------------------------------------
  tau2_fn <- .make_tau2_fn(vfun)

  # ---- Starting values: classical DerSimonian-Laird tau^2 -----------------
  w0      <- 1 / vi
  mu0     <- sum(w0 * yi) / sum(w0)
  Q0      <- sum(w0 * (yi - mu0)^2)
  C0      <- sum(w0) - sum(w0^2) / sum(w0)
  tau2_DL <- max((Q0 - (k - 1)) / C0, 1e-6)

  # Start: tau0sq ≈ tau2_DL, gamma ≈ 1
  theta0 <- c(log(tau2_DL), log(1))

  # ---- Optimisation --------------------------------------------------------
  ctrl_def <- list(maxit = 2000, factr = 1e7)
  ctrl     <- modifyList(ctrl_def, control)

  opt <- tryCatch(
    stats::optim(
      par     = theta0,
      fn      = .drmeta_loglik,
      yi      = yi, vi = vi, dr = dr,
      tau2_fn = tau2_fn,
      drmeta_method = method,
      method  = "L-BFGS-B",
      lower   = c(-20, -10),
      upper   = c(10,  10),
      control = ctrl
    ),
    error = function(e) {
      stats::optim(
        par     = theta0,
        fn      = .drmeta_loglik,
        yi      = yi, vi = vi, dr = dr,
        tau2_fn = tau2_fn,
        drmeta_method = method,
        method  = "Nelder-Mead",
        control = ctrl
      )
    }
  )

  tau0sq <- exp(opt$par[1])
  gamma  <- exp(opt$par[2])

  # ---- Derived quantities --------------------------------------------------
  tau2_i   <- tau2_fn(dr, tau0sq, gamma)
  sigma2_i <- vi + tau2_i

  w        <- 1 / sigma2_i
  mu       <- sum(w * yi) / sum(w)
  se       <- sqrt(1 / sum(w))
  ci.lb    <- mu - qnorm(0.975) * se
  ci.ub    <- mu + qnorm(0.975) * se
  zval     <- mu / se
  pval     <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  weights  <- w / sum(w)  # normalised to sum 1

  # ---- Log-likelihoods and information criteria ----------------------------
  resid    <- yi - mu
  ll_ml    <- -0.5 * (sum(log(sigma2_i)) + sum(resid^2 / sigma2_i))
  ll_reml  <- ll_ml - 0.5 * log(sum(w))

  np  <- 3L  # mu, tau0sq, gamma (ML parameters; for AIC/BIC use ML)
  AIC <- -2 * ll_ml + 2 * np
  BIC <- -2 * ll_ml + log(k) * np

  # ---- Return --------------------------------------------------------------
  out <- list(
    mu          = mu,
    se          = se,
    ci.lb       = ci.lb,
    ci.ub       = ci.ub,
    zval        = zval,
    pval        = pval,
    tau0sq      = tau0sq,
    gamma       = gamma,
    tau2_i      = tau2_i,
    sigma2_i    = sigma2_i,
    weights     = weights,
    loglik      = ll_ml,
    reml_loglik = ll_reml,
    AIC         = AIC,
    BIC         = BIC,
    k           = k,
    yi          = yi,
    vi          = vi,
    dr          = dr,
    slab        = slab,
    vfun        = vfun,
    method      = method,
    converged   = (opt$convergence == 0L),
    optim_out   = opt,
    call        = cl
  )
  class(out) <- "drmeta"
  out
}


# =============================================================================
# S3 methods
# =============================================================================

#' Print Method for drmeta Objects
#'
#' @param x A fitted `"drmeta"` object.
#' @param digits Number of significant digits.  Default 4.
#' @param ... Ignored.
#' @return Invisibly returns the original \code{drmeta} object \code{x},
#'   unchanged. This function is called for its side effect of printing a
#'   formatted summary of the fitted DR-Meta model to the console.
#' @export
print.drmeta <- function(x, digits = 4, ...) {
  cat("\n--- DR-Meta: Design-Robust Random-Effects Model ---\n\n")
  cat(sprintf("Variance function : %s\n", x$vfun))
  cat(sprintf("Estimation method : %s\n", x$method))
  cat(sprintf("Studies (k)       : %d\n\n", x$k))

  cat(sprintf("Pooled estimate   : %.*f\n", digits, x$mu))
  cat(sprintf("Std. error        : %.*f\n", digits, x$se))
  cat(sprintf("95%% CI            : [%.*f, %.*f]\n",
              digits, x$ci.lb, digits, x$ci.ub))
  cat(sprintf("z = %.*f,  p = %.*f\n\n", digits, x$zval, digits, x$pval))

  cat(sprintf("tau0^2 (DR=0)     : %.*f\n", digits, x$tau0sq))
  cat(sprintf("gamma             : %.*f\n", digits, x$gamma))
  if (!x$converged)
    cat("WARNING: optimiser did not converge.\n")
  invisible(x)
}


#' Summary Method for drmeta Objects
#'
#' @param object A fitted `"drmeta"` object.
#' @param digits Number of significant digits.  Default 4.
#' @param ... Ignored.
#' @return Invisibly returns the fitted \code{drmeta} object \code{object},
#'   unchanged. Called for its side effect of printing a detailed formatted
#'   summary — including the pooled estimate, confidence interval, z-test,
#'   variance-function parameters, and model fit statistics — to the console.
#' @export
summary.drmeta <- function(object, digits = 4, ...) {
  cat("\n=== DR-Meta Summary ===\n\n")

  cat("Call:\n  ")
  cat(deparse(object$call), "\n\n")

  cat("--- Pooled Effect ---\n")
  cat(sprintf("  mu   = %.*f  (SE = %.*f)\n",
              digits, object$mu, digits, object$se))
  cat(sprintf("  95%% CI: [%.*f, %.*f]\n",
              digits, object$ci.lb, digits, object$ci.ub))
  cat(sprintf("  z = %.*f,  p = %.*f\n\n",
              digits, object$zval, digits, object$pval))

  cat("--- Variance-Function Parameters ---\n")
  cat(sprintf("  tau0^2 = %.*f  (heterogeneity at DR=0)\n",
              digits, object$tau0sq))
  cat(sprintf("  gamma  = %.*f  (decay rate)\n",
              digits, object$gamma))
  cat(sprintf("  vfun   = \"%s\",  method = \"%s\"\n\n",
              object$vfun, object$method))

  cat("--- Model Fit ---\n")
  cat(sprintf("  logLik (ML)   = %.*f\n", digits, object$loglik))
  cat(sprintf("  logLik (REML) = %.*f\n", digits, object$reml_loglik))
  cat(sprintf("  AIC = %.*f,  BIC = %.*f\n\n",
              digits, object$AIC, digits, object$BIC))

  cat(sprintf("Converged: %s\n", object$converged))
  invisible(object)
}


#' Extract Coefficients from a drmeta Object
#'
#' Returns a named numeric vector of the three model parameters:
#' `mu`, `tau0sq`, and `gamma`.
#'
#' @param object A fitted `"drmeta"` object.
#' @param ... Ignored.
#' @return A named numeric vector of length 3 with the estimated model
#'   parameters: \code{mu} (pooled effect estimate), \code{tau0sq}
#'   (baseline between-study variance at DR = 0), and \code{gamma}
#'   (variance-function decay rate).
#' @export
coef.drmeta <- function(object, ...) {
  c(mu = object$mu, tau0sq = object$tau0sq, gamma = object$gamma)
}


#' Confidence Interval for a drmeta Object
#'
#' Returns a data frame with the estimate and 95% confidence interval for
#' the pooled effect \eqn{\hat\mu}.
#'
#' @param object A fitted `"drmeta"` object.
#' @param parm Ignored (only `mu` is returned).
#' @param level Confidence level.  Default 0.95.
#' @param ... Ignored.
#' @return A data frame with one row (\code{mu}) and three columns:
#'   \code{estimate} (the pooled effect \eqn{\hat\mu}), \code{lower},
#'   and \code{upper} (confidence interval bounds at the requested
#'   \code{level}, default 95\%).
#' @export
confint.drmeta <- function(object, parm = NULL, level = 0.95, ...) {
  alpha <- 1 - level
  z     <- qnorm(1 - alpha / 2)
  data.frame(
    estimate = object$mu,
    lower    = object$mu - z * object$se,
    upper    = object$mu + z * object$se,
    row.names = "mu"
  )
}


#' Fitted Values for a drmeta Object
#'
#' Returns a vector of length \eqn{k} where every element equals the pooled
#' estimate \eqn{\hat\mu} (the model has a single intercept, so all fitted
#' values are identical).
#'
#' @param object A fitted `"drmeta"` object.
#' @param ... Ignored.
#' @return A numeric vector of length \eqn{k} where every element equals
#'   the pooled estimate \eqn{\hat\mu}. Because DR-Meta has a single
#'   intercept, all studies share the same fitted value.
#' @export
fitted.drmeta <- function(object, ...) {
  rep(object$mu, object$k)
}


#' Residuals for a drmeta Object
#'
#' @param object A fitted `"drmeta"` object.
#' @param type `"raw"` (default) or `"standardised"`.
#' @param ... Ignored.
#' @return A numeric vector of length \eqn{k} of residuals. When
#'   \code{type = "raw"} (default), returns observed minus fitted values
#'   (\eqn{y_i - \hat\mu}). When \code{type = "standardised"}, each
#'   residual is divided by \eqn{\sqrt{\hat\sigma^2_i}} (the square root
#'   of the total study variance under the fitted model).
#' @export
residuals.drmeta <- function(object, type = c("raw", "standardised"), ...) {
  type  <- match.arg(type)
  resid <- object$yi - object$mu
  if (type == "standardised") resid <- resid / sqrt(object$sigma2_i)
  resid
}


#' Log-Likelihood for a drmeta Object
#'
#' @param object A fitted `"drmeta"` object.
#' @param REML Logical.  If `TRUE`, returns the REML log-likelihood.
#'   Default `FALSE` (ML).
#' @param ... Ignored.
#' @return An object of class \code{"logLik"}. The numeric value is the
#'   maximised log-likelihood (ML or REML, depending on \code{REML}).
#'   The object carries two attributes: \code{df} (number of parameters,
#'   always 3: \code{mu}, \code{tau0sq}, \code{gamma}) and \code{nobs}
#'   (number of studies \eqn{k}).
#' @export
logLik.drmeta <- function(object, REML = FALSE, ...) {
  val <- if (REML) object$reml_loglik else object$loglik
  attr(val, "df")    <- 3L
  attr(val, "nobs")  <- object$k
  class(val) <- "logLik"
  val
}
