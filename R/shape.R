#' Check whether a monotone scale model misses a grouped pattern
#'
#' Fits a categorical scale model and a constant-scale model with
#' \code{metafor::rma()}, using the location model already stored in a fitted
#' \code{drmeta} object. This diagnostic is useful when the design score has a
#' small number of ordered categories: an interior peak or trough cannot be
#' represented by the package's monotone exponential variance function.
#'
#' @param object A fitted \code{drmeta} object.
#' @param groups A vector defining scale groups. By default the observed
#'   design-robustness values are used. For an ordered interpretation, supply
#'   a factor whose levels are in the intended order.
#' @param method Either \code{"ML"} or \code{"REML"}. ML is the default for
#'   the likelihood-ratio comparison.
#' @return An object of class \code{dr_shape_check}. It contains grouped
#'   variance estimates, the likelihood-ratio statistic and p-value, the
#'   inferred qualitative pattern, and the two fitted \code{metafor} models.
#' @details The likelihood-ratio p-value is a descriptive model check. It uses
#'   the usual chi-square reference with one fewer degree of freedom than the
#'   number of groups. Estimates at or near zero can make that reference
#'   approximate. The diagnostic does not select a new primary score or alter
#'   the fitted \code{drmeta} model.
#' @examples
#' \dontrun{
#' path <- system.file("extdata", "bcg_design_robustness.csv", package = "drmeta")
#' bcg <- utils::read.csv(path)
#' fit <- drmeta(bcg[["yi"]], bcg[["vi"]], bcg[["dr"]])
#' dr_shape_check(fit)
#' }
#' @export
dr_shape_check <- function(object, groups = object$dr,
                           method = c("ML", "REML")) {
  if (!inherits(object, "drmeta")) stop("object must be a drmeta fit.")
  if (!requireNamespace("metafor", quietly = TRUE))
    stop("dr_shape_check() requires the suggested package 'metafor'.")
  method <- match.arg(method)
  if (length(groups) != object$k || anyNA(groups))
    stop("groups must contain one non-missing value per study.")

  if (is.factor(groups)) {
    group <- droplevels(groups)
  } else if (is.numeric(groups)) {
    group <- factor(groups, levels = sort(unique(groups)))
  } else {
    group <- factor(groups, levels = unique(groups))
  }
  if (nlevels(group) < 2L)
    stop("groups must contain at least two distinct values.")

  dat <- data.frame(yi = object$yi, vi = object$vi, group = group)
  mods <- if (object$p > 1L) object$X[, -1, drop = FALSE] else NULL
  if (is.null(mods)) {
    full <- metafor::rma(yi = dat$yi, vi = dat$vi, scale = ~ group,
                         data = dat, method = method)
    null <- metafor::rma(yi = dat$yi, vi = dat$vi, scale = ~ 1,
                         data = dat, method = method)
  } else {
    full <- metafor::rma(yi = dat$yi, vi = dat$vi, mods = mods,
                         scale = ~ group, data = dat, method = method)
    null <- metafor::rma(yi = dat$yi, vi = dat$vi, mods = mods,
                         scale = ~ 1, data = dat, method = method)
  }

  lev <- levels(group)
  group_grid <- data.frame(group = factor(lev, levels = lev))
  scale_matrix <- stats::model.matrix(~ group, data = group_grid)
  tau2 <- exp(as.vector(scale_matrix %*% as.numeric(full$alpha)))
  estimates <- data.frame(
    group = lev,
    n = as.integer(table(group)[lev]),
    tau2 = tau2,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  delta <- diff(tau2)
  tol <- sqrt(.Machine$double.eps) * max(1, tau2)
  pattern <- if (all(delta <= tol)) {
    "nonincreasing"
  } else if (all(delta >= -tol)) {
    "nondecreasing"
  } else if (length(tau2) >= 3L &&
             which.max(tau2) %in% seq.int(2L, length(tau2) - 1L)) {
    "interior peak"
  } else if (length(tau2) >= 3L &&
             which.min(tau2) %in% seq.int(2L, length(tau2) - 1L)) {
    "interior trough"
  } else {
    "nonmonotone"
  }

  lr <- max(0, 2 * (as.numeric(stats::logLik(full)) -
                    as.numeric(stats::logLik(null))))
  df <- nlevels(group) - 1L
  structure(
    list(estimates = estimates, statistic = lr, df = df,
         p.value = stats::pchisq(lr, df = df, lower.tail = FALSE),
         pattern = pattern, method = method, full = full, null = null,
         call = match.call()),
    class = "dr_shape_check"
  )
}

#' Print a grouped scale-shape diagnostic
#'
#' @param x A \code{dr_shape_check} object.
#' @param digits Number of significant digits.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}.
#' @export
print.dr_shape_check <- function(x,
                                 digits = max(3L, getOption("digits") - 3L),
                                 ...) {
  cat("Grouped scale-shape diagnostic\n\n")
  print(x$estimates, digits = digits, row.names = FALSE)
  cat("\nPattern:", x$pattern, "\n")
  cat("LR =", format(x$statistic, digits = digits),
      " df =", x$df,
      " p =", format(x$p.value, digits = digits), "\n")
  invisible(x)
}
