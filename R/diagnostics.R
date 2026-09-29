#' Predict the fitted between-study variance
#'
#' Evaluates the fitted variance function \eqn{\tau^2(d) = \tau_0^2
#' \exp(-\gamma d)} over a grid of design-robustness values. The default grid
#' spans the observed range of the design index rather than the full interval
#' \eqn{[0, 1]}, because values outside the observed support are
#' extrapolations of the fitted scale model.
#'
#' @param object A fitted \code{drmeta} object.
#' @param dr New design-robustness values in \eqn{[0, 1]}. Defaults to a grid
#'   spanning the observed range.
#' @return A data frame with columns \code{dr} and \code{tau2}, giving the
#'   fitted residual between-study variance at each supplied value.
#' @examples
#' path <- system.file("extdata", "bcg_design_robustness.csv", package = "drmeta")
#' bcg <- utils::read.csv(path)
#' fit <- drmeta(yi = bcg[["yi"]], vi = bcg[["vi"]], dr = bcg[["dr"]])
#' head(dr_scale_predict(fit))
#' @export
dr_scale_predict <- function(object, dr = seq(min(object$dr), max(object$dr), length.out = 100)) {
  if (!inherits(object, "drmeta")) stop("object must be a drmeta fit.")
  dr <- as.numeric(dr)
  if (any(!is.finite(dr))) stop("dr must be finite.")
  if (any(dr < 0 | dr > 1)) stop("dr must lie in [0,1].")
  data.frame(dr = dr, tau2 = object$tau0sq * exp(-object$gamma * dr))
}

#' Proportional attenuation in fitted heterogeneity
#'
#' Returns the proportional reduction in fitted residual between-study
#' heterogeneity between two design-robustness values,
#' \eqn{A_\tau = 1 - \exp(-\gamma (d_H - d_L))}. This is a summary of the
#' fitted scale function only. It is not a causal \eqn{R^2}, not a proportion
#' of heterogeneity explained, and not evidence that design quality caused a
#' reduction in heterogeneity.
#'
#' @param object A fitted \code{drmeta} object.
#' @param d_low Lower DR value; defaults to the observed minimum.
#' @param d_high Higher DR value; defaults to the observed maximum.
#' @return A single numeric value. For a constrained fit it lies in
#'   \eqn{[0, 1)}; for an unrestricted fit with a negative gradient it is
#'   negative, indicating heterogeneity that increases with the design index.
#' @examples
#' path <- system.file("extdata", "bcg_design_robustness.csv", package = "drmeta")
#' bcg <- utils::read.csv(path)
#' fit <- drmeta(yi = bcg[["yi"]], vi = bcg[["vi"]], dr = bcg[["dr"]])
#' dr_scale_attenuation(fit)
#' dr_scale_attenuation(fit, d_low = 0.25, d_high = 0.75)
#' @export
dr_scale_attenuation <- function(object, d_low = min(object$dr), d_high = max(object$dr)) {
  if (!inherits(object, "drmeta")) stop("object must be a drmeta fit.")
  if (!is.numeric(d_low) || !is.numeric(d_high) ||
      length(d_low) != 1L || length(d_high) != 1L ||
      !is.finite(d_low) || !is.finite(d_high))
    stop("d_low and d_high must each be a single finite number.")
  if (d_high < d_low) stop("d_high must be at least d_low.")
  1 - exp(-object$gamma * (d_high - d_low))
}

#' Classical heterogeneity diagnostics and design-indexed summaries
#'
#' Reports the classical fixed-effect residual \eqn{Q} statistic and \eqn{I^2}
#' alongside the fitted design-indexed scale parameters. Unlike version 0.1.0,
#' this function does not return a design-explained variance decomposition;
#' that quantity was not identified by the model and has been replaced by
#' \code{\link{dr_scale_attenuation}}.
#'
#' @param object A fitted \code{drmeta} object.
#' @return A list with components \code{Q}, \code{df}, \code{p}, \code{I2},
#'   \code{tau0sq}, \code{gamma}, \code{attenuation_observed}, and
#'   \code{tau2_by_study}.
#' @examples
#' path <- system.file("extdata", "bcg_design_robustness.csv", package = "drmeta")
#' bcg <- utils::read.csv(path)
#' fit <- drmeta(yi = bcg[["yi"]], vi = bcg[["vi"]], dr = bcg[["dr"]])
#' dr_heterogeneity(fit)
#' @export
dr_heterogeneity <- function(object) {
  if (!inherits(object, "drmeta")) stop("object must be a drmeta fit.")
  yi <- object$yi; vi <- object$vi; X <- object$X
  w_fe <- 1 / vi
  b_fe <- solve(crossprod(X, w_fe * X), crossprod(X, w_fe * yi))
  r_fe <- as.vector(yi - X %*% b_fe)
  Q <- sum(w_fe * r_fe^2)
  df <- object$k - object$p
  I2 <- if (Q > 0) max(0, (Q - df) / Q) else 0
  list(Q = Q, df = df, p = stats::pchisq(Q, df, lower.tail = FALSE),
       I2 = I2, tau0sq = object$tau0sq, gamma = object$gamma,
       attenuation_observed = dr_scale_attenuation(object),
       tau2_by_study = object$tau2i)
}

#' Leave-one-out influence diagnostics
#'
#' Refits the model with each study omitted in turn and records how the
#' location and scale estimates move. Both components are reported, because a
#' study can be uninfluential for the pooled mean while dominating the scale
#' gradient.
#'
#' @param object A fitted \code{drmeta} object.
#' @param slab Optional study labels. By default, labels stored in
#'   \code{object$slab} are used; row-number labels are used for older objects
#'   without that component.
#' @return A data frame with one row per omitted study, containing the study
#'   label, its design-robustness value, the leave-one-out estimate of the
#'   first location coefficient and its change from the full fit, the
#'   leave-one-out \code{tau0sq} and \code{gamma} with their changes, a
#'   \code{converged} flag, and a logical \code{influential} column flagging
#'   changes in the first location coefficient exceeding twice its
#'   full-model standard error.
#' @examples
#' path <- system.file("extdata", "bcg_design_robustness.csv", package = "drmeta")
#' bcg <- utils::read.csv(path)
#' fit <- drmeta(yi = bcg[["yi"]], vi = bcg[["vi"]], dr = bcg[["dr"]])
#' dr_loo(fit)
#' @export
dr_loo <- function(object, slab = object$slab) {
  if (!inherits(object, "drmeta")) stop("object must be a drmeta fit.")
  k <- object$k
  if (k < 4L) stop("Leave-one-out requires at least four studies.")
  if (is.null(slab)) slab <- as.character(seq_len(k))
  if (length(slab) != k || anyNA(slab))
    stop("slab must contain one non-missing label per study.")
  slab <- as.character(slab)
  mods_full <- if (object$p > 1L) object$X[, -1, drop = FALSE] else NULL
  se_full <- sqrt(object$vcov[1, 1])
  res <- lapply(seq_len(k), function(i) {
    fit <- tryCatch(
      drmeta(object$yi[-i], object$vi[-i], object$dr[-i],
             mods = if (is.null(mods_full)) NULL else mods_full[-i, , drop = FALSE],
             slab = slab[-i],
             method = object$method, constrained = object$constrained,
             gamma_max = object$gamma_max, gamma_fixed = object$gamma_fixed,
             .quiet = TRUE),
      error = function(e) NULL)
    if (is.null(fit))
      return(data.frame(est_loo = NA_real_, tau0sq_loo = NA_real_,
                        gamma_loo = NA_real_, converged = FALSE))
    data.frame(est_loo = unname(fit$beta[1]), tau0sq_loo = fit$tau0sq,
               gamma_loo = fit$gamma, converged = identical(as.numeric(fit$convergence), 0))
  })
  out <- do.call(rbind, res)
  out <- cbind(study = slab, dr = object$dr, out,
               stringsAsFactors = FALSE)
  out$delta_est <- out$est_loo - unname(object$beta[1])
  out$delta_tau0sq <- out$tau0sq_loo - object$tau0sq
  out$delta_gamma <- out$gamma_loo - object$gamma
  out$influential <- is.finite(out$delta_est) & abs(out$delta_est) > 2 * se_full
  rownames(out) <- NULL
  out[, c("study", "dr", "est_loo", "delta_est", "tau0sq_loo", "delta_tau0sq",
          "gamma_loo", "delta_gamma", "converged", "influential")]
}
