#' Fit a design-indexed meta-analytic location-scale model
#'
#' @param yi Numeric vector of study effect estimates.
#' @param vi Numeric vector of sampling variances.
#' @param dr Numeric design-robustness index in \eqn{[0, 1]}.
#' @param mods Optional numeric vector or matrix of location moderators. An
#'   intercept is added automatically.
#' @param slab Optional study labels. Defaults to the row numbers and is used
#'   by study-level diagnostics such as \code{\link{dr_loo}}.
#' @param method Either "REML" or "ML".
#' @param constrained If TRUE, impose gamma >= 0. If FALSE, allow either sign.
#' @param gamma_max Finite optimization bound for abs(gamma). The default,
#'   8, matches the range used in the simulation study reported in the
#'   accompanying manuscript. Refit under a wider bound whenever the estimate
#'   piles up near it.
#' @param gamma_fixed Optional fixed value for gamma, used for null fitting.
#' @param control Optional list passed to optim().
#' @param .quiet If TRUE, suppress the design-diagnostic warnings about the
#'   distribution of \code{dr}. Intended for internal refitting, where the
#'   design index is identical across calls and the diagnostics would repeat
#'   without conveying new information. Errors and input validation are
#'   unaffected.
#' @return An object of class \code{drmeta}: a list whose components include
#'   \code{beta} (location coefficients), \code{vcov}, \code{ci_beta},
#'   \code{tau0sq}, \code{gamma}, \code{tau2i}, \code{sigma2},
#'   \code{weights}, \code{residuals}, \code{slab}, \code{logLik}, \code{method},
#'   \code{constrained}, \code{convergence}, and the supplied data.
#' @examples
#' path <- system.file("extdata", "bcg_design_robustness.csv", package = "drmeta")
#' bcg <- utils::read.csv(path)
#' fit <- drmeta(yi = bcg[["yi"]], vi = bcg[["vi"]], dr = bcg[["dr"]])
#' summary(fit)
#' @export
drmeta <- function(yi, vi, dr, mods = NULL, slab = NULL,
                   method = c("REML", "ML"),
                   constrained = TRUE, gamma_max = 8,
                   gamma_fixed = NULL, control = list(),
                   .quiet = FALSE) {
  method <- match.arg(method)
  yi <- as.numeric(yi); vi <- as.numeric(vi); dr <- as.numeric(dr)
  k <- length(yi)
  if (length(vi) != k || length(dr) != k) stop("yi, vi, and dr must have equal length.")
  if (k < 3L || any(!is.finite(yi)) || any(!is.finite(vi)) || any(!is.finite(dr)))
    stop("Inputs must be finite and contain at least three studies.")
  if (any(vi <= 0)) stop("All sampling variances must be positive.")
  if (any(dr < 0 | dr > 1)) stop("dr must lie in [0,1].")
  if (is.null(slab)) slab <- as.character(seq_len(k))
  if (length(slab) != k || anyNA(slab))
    stop("slab must contain one non-missing label per study.")
  slab <- as.character(slab)
  if (!is.null(mods)) {
    M <- as.matrix(mods)
    if (nrow(M) != k) stop("mods must have one row per study.")
    X <- cbind(`(Intercept)` = 1, M)
  } else X <- matrix(1, nrow = k, ncol = 1, dimnames = list(NULL, "(Intercept)"))
  p <- ncol(X)
  if (k <= p + 1L) stop("Too few studies for the requested location model.")
  if (qr(X)$rank < p) stop("Location-model matrix is rank deficient.")
  
  gls_at <- function(log_tau0sq, gamma) {
    tau0sq <- exp(log_tau0sq)
    tau2i <- tau0sq * exp(-gamma * dr)
    sigma2 <- vi + tau2i
    w <- 1 / sigma2
    XtWX <- crossprod(X, w * X)
    beta <- tryCatch(solve(XtWX, crossprod(X, w * yi)), error = function(e) NULL)
    if (is.null(beta)) return(NULL)
    res <- as.vector(yi - X %*% beta)
    list(beta = as.vector(beta), vcov = solve(XtWX), residuals = res,
         tau0sq = tau0sq, gamma = gamma, tau2i = tau2i,
         sigma2 = sigma2, weights = w, XtWX = XtWX)
  }
  
  objective <- function(par) {
    log_tau <- par[1]
    gamma <- if (is.null(gamma_fixed)) par[2] else gamma_fixed
    z <- gls_at(log_tau, gamma)
    if (is.null(z)) return(.Machine$double.xmax / 100)
    quad <- sum(z$residuals^2 / z$sigma2)
    if (method == "ML") {
      0.5 * (k * log(2*pi) + sum(log(z$sigma2)) + quad)
    } else {
      ld <- as.numeric(determinant(z$XtWX, logarithm = TRUE)$modulus)
      0.5 * ((k-p) * log(2*pi) + sum(log(z$sigma2)) + ld + quad)
    }
  }
  
  # Stable moment-based start. The exact value is not inferentially important.
  w_fe <- 1 / vi
  beta_fe <- solve(crossprod(X, w_fe * X), crossprod(X, w_fe * yi))
  r_fe <- as.vector(yi - X %*% beta_fe)
  tau_start <- max(1e-8, stats::var(r_fe) - mean(vi))
  log_start <- log(tau_start)
  ctl <- utils::modifyList(list(maxit = 2000, factr = 1e7), control)
  
  if (is.null(gamma_fixed)) {
    lower_g <- if (constrained) 0 else -abs(gamma_max)
    upper_g <- abs(gamma_max)
    opt <- stats::optim(c(log_start, 0), objective, method = "L-BFGS-B",
                        lower = c(log(.Machine$double.eps), lower_g),
                        upper = c(log(max(1, stats::var(yi) * 1e4)), upper_g),
                        control = ctl)
    gamma_hat <- opt$par[2]
    # Explicit boundary comparison avoids numerical failure to select gamma=0.
    if (constrained) {
      opt0 <- stats::optimize(function(x) objective(c(x, 0)),
                              interval = c(log(.Machine$double.eps), log(max(1, stats::var(yi) * 1e4))))
      if (opt0$objective <= opt$value + 1e-10) {
        opt$par <- c(opt0$minimum, 0); opt$value <- opt0$objective; gamma_hat <- 0
      }
    }
  } else {
    if (!is.finite(gamma_fixed)) stop("gamma_fixed must be finite.")
    if (constrained && gamma_fixed < 0) stop("gamma_fixed must be nonnegative for a constrained fit.")
    opt0 <- stats::optimize(function(x) objective(c(x, gamma_fixed)),
                            interval = c(log(.Machine$double.eps), log(max(1, stats::var(yi) * 1e4))))
    opt <- list(par = c(opt0$minimum, gamma_fixed), value = opt0$objective,
                convergence = 0, message = NULL)
    gamma_hat <- gamma_fixed
  }
  
  fit <- gls_at(opt$par[1], gamma_hat)
  df <- k - p
  crit <- stats::qt(.975, df = df)
  ci_beta <- cbind(fit$beta - crit * sqrt(diag(fit$vcov)),
                   fit$beta + crit * sqrt(diag(fit$vcov)))
  colnames(ci_beta) <- c("ci.lb", "ci.ub"); rownames(ci_beta) <- colnames(X)
  names(fit$beta) <- colnames(X); dimnames(fit$vcov) <- list(colnames(X), colnames(X))
  
  out <- c(fit, list(yi = yi, vi = vi, dr = dr, slab = slab,
                     X = X, k = k, p = p, df = df,
                     method = method, constrained = constrained, gamma_max = gamma_max,
                     gamma_fixed = gamma_fixed, logLik = -opt$value,
                     convergence = opt$convergence, message = opt$message,
                     ci_beta = ci_beta, call = match.call()))
  class(out) <- "drmeta"
  
  # Design diagnostics describe the dr vector only. They are invariant across
  # refits that hold dr fixed, so internal callers suppress them via .quiet.
  if (!.quiet) {
    if (length(unique(dr)) < 3L)
      warning("The design-robustness index has fewer than three distinct values; gamma may be weakly identified.")
    if (diff(range(dr)) < .20)
      warning("The observed DR range is narrow; interpret the scale gradient cautiously.")
    if (is.null(gamma_fixed) && abs(abs(gamma_hat) - abs(gamma_max)) < 1e-3)
      warning("The scale gradient is at the optimization bound; refit under a wider gamma_max before interpreting it.")
  }
  out
}

#' @rdname drmeta-methods
#' @export
coef.drmeta <- function(object, ...) object$beta

#' @rdname drmeta-methods
#' @export
vcov.drmeta <- function(object, ...) object$vcov

#' @rdname drmeta-methods
#' @export
fitted.drmeta <- function(object, ...) as.vector(object$X %*% object$beta)

#' @rdname drmeta-methods
#' @export
residuals.drmeta <- function(object, ...) object$residuals

#' @rdname drmeta-methods
#' @export
print.drmeta <- function(x, ...) {
  cat("Design-indexed location-scale meta-analysis\n")
  cat("Method:", x$method, " | k =", x$k, " | gamma constraint:",
      if (x$constrained) ">= 0" else "unrestricted", "\n")
  print(round(x$beta, 5))
  cat("tau0^2 =", format(x$tau0sq, digits = 5),
      " gamma =", format(x$gamma, digits = 5), "\n")
  invisible(x)
}

#' @rdname drmeta-methods
#' @export
summary.drmeta <- function(object, ...) {
  se <- sqrt(diag(object$vcov)); tval <- object$beta / se
  tab <- cbind(estimate = object$beta, se = se, t = tval,
               p = 2 * stats::pt(abs(tval), df = object$df, lower.tail = FALSE),
               object$ci_beta)
  out <- list(coefficients = tab, tau0sq = object$tau0sq, gamma = object$gamma,
              attenuation = dr_scale_attenuation(object), logLik = object$logLik,
              method = object$method, k = object$k, df = object$df,
              convergence = object$convergence, constrained = object$constrained)
  class(out) <- "summary.drmeta"; out
}

#' @rdname drmeta-methods
#' @export
print.summary.drmeta <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Design-indexed location-scale meta-analysis\n\n")
  
  # printCoefmat() expects the p-value column to be the final column.
  # Confidence limits are therefore printed separately.
  test_columns <- c("estimate", "se", "t", "p")
  ci_columns <- c("ci.lb", "ci.ub")
  
  stats::printCoefmat(
    x$coefficients[, test_columns, drop = FALSE],
    digits = digits,
    P.values = TRUE,
    has.Pvalue = TRUE
  )
  
  cat("\n", round(100 * 0.95), "% confidence intervals:\n", sep = "")
  print(
    round(x$coefficients[, ci_columns, drop = FALSE], digits = digits),
    quote = FALSE
  )
  
  cat(
    "\ntau0^2:", format(x$tau0sq, digits = digits),
    " gamma:", format(x$gamma, digits = digits),
    " attenuation over observed DR range:",
    format(x$attenuation, digits = digits), "\n"
  )
  cat(
    "Criterion:", x$method,
    " logLik =", format(x$logLik, digits = digits),
    " convergence =", x$convergence, "\n"
  )
  invisible(x)
}

#' S3 methods for drmeta objects
#'
#' @param object,x A fitted \code{drmeta} object.
#' @param parm Character or numeric selection of location parameters.
#' @param level Confidence level.
#' @param REML If TRUE, report the restricted log-likelihood.
#' @param digits Number of significant digits.
#' @param ... Ignored.
#' @return \code{coef} and \code{fitted} and \code{residuals} return numeric
#'   vectors; \code{vcov} returns the location covariance matrix;
#'   \code{confint} returns a matrix of confidence limits; \code{logLik}
#'   returns an object of class \code{logLik} carrying \code{df} and
#'   \code{nobs} attributes; \code{summary} returns an object of class
#'   \code{summary.drmeta}; the print methods return their argument invisibly.
#' @name drmeta-methods
NULL

#' @rdname drmeta-methods
#' @export
confint.drmeta <- function(object, parm = NULL, level = 0.95, ...) {
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    stop("level must be a single number in (0, 1).")
  se <- sqrt(diag(object$vcov))
  crit <- stats::qt(1 - (1 - level) / 2, df = object$df)
  out <- cbind(object$beta - crit * se, object$beta + crit * se)
  colnames(out) <- paste0(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2),
                                 trim = TRUE), " %")
  rownames(out) <- names(object$beta)
  if (!is.null(parm)) out <- out[parm, , drop = FALSE]
  out
}

#' @rdname drmeta-methods
#' @export
logLik.drmeta <- function(object, REML = NULL, ...) {
  is_reml <- if (is.null(REML)) identical(object$method, "REML") else isTRUE(REML)
  if (!is.null(REML) && is_reml != identical(object$method, "REML"))
    stop("The fit was obtained under ", object$method,
         "; refit with method = ", if (is_reml) "\"REML\"" else "\"ML\"",
         " rather than relabelling the criterion.")
  val <- object$logLik
  # Scale parameters: tau0^2 always, plus gamma when it was estimated.
  n_scale <- 1L + as.integer(is.null(object$gamma_fixed))
  attr(val, "df") <- object$p + n_scale
  attr(val, "nobs") <- if (is_reml) object$k - object$p else object$k
  class(val) <- "logLik"
  val
}
