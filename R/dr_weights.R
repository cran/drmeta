# =============================================================================
# dr_weights.R  —  Weight and variance-function utilities
# =============================================================================

#' Compute DR-Meta Study Weights
#'
#' Returns the DR-Meta inverse-total-variance weights
#' \deqn{w_i = \frac{1}{v_i + \hat\tau^2(\mathrm{DR}_i)},}
#' given estimated variance-function parameters.  Optionally normalises weights
#' to sum to 1 or to 100.
#'
#' This function is primarily a utility for diagnostics and visualisation;
#' weights are also returned as part of the \code{"drmeta"} object produced by
#' \code{\link{drmeta}}.
#'
#' @param vi Numeric vector of sampling variances (\eqn{v_i > 0}).
#' @param dr Numeric vector of design robustness indices in \eqn{[0,1]}.
#' @param tau0sq Non-negative scalar: estimated baseline heterogeneity
#'   \eqn{\hat\tau_0^2}.
#' @param gamma Non-negative scalar: estimated variance-function slope
#'   \eqn{\hat\gamma}.
#' @param vfun Variance function: `"exponential"` (default) or `"linear"`.
#' @param normalise Character: `"none"` (raw weights, default), `"sum1"`
#'   (divide by sum so weights sum to 1), or `"pct"` (multiply by 100 after
#'   normalising to sum to 1).
#'
#' @return A numeric vector of weights, the same length as `vi`.
#'
#' @examples
#' vi  <- c(0.02, 0.03, 0.015, 0.025, 0.01)
#' dr  <- c(0.9,  0.4,  0.7,   0.2,   1.0)
#' dr_weights(vi, dr, tau0sq = 0.04, gamma = 1.5)
#' dr_weights(vi, dr, tau0sq = 0.04, gamma = 1.5, normalise = "pct")
#'
#' @seealso \code{\link{drmeta}}, \code{\link{dr_variance}}
#' @export
dr_weights <- function(vi, dr,
                       tau0sq,
                       gamma,
                       vfun      = c("exponential", "linear"),
                       normalise = c("none", "sum1", "pct")) {

  vfun      <- match.arg(vfun)
  normalise <- match.arg(normalise)

  if (!is.numeric(vi) || any(vi <= 0))
    stop("`vi` must be a numeric vector of positive values.", call. = FALSE)
  if (!is.numeric(dr) || any(dr < 0 | dr > 1))
    stop("`dr` must be numeric with values in [0,1].", call. = FALSE)
  if (length(vi) != length(dr))
    stop("`vi` and `dr` must have the same length.", call. = FALSE)
  if (!is.numeric(tau0sq) || length(tau0sq) != 1L || tau0sq < 0)
    stop("`tau0sq` must be a non-negative scalar.", call. = FALSE)
  if (!is.numeric(gamma) || length(gamma) != 1L || gamma < 0)
    stop("`gamma` must be a non-negative scalar.", call. = FALSE)

  tau2_fn  <- .make_tau2_fn(vfun)
  tau2_i   <- tau2_fn(dr, tau0sq, gamma)
  sigma2_i <- vi + tau2_i
  w        <- 1 / sigma2_i

  if (normalise == "sum1") w <- w / sum(w)
  if (normalise == "pct")  w <- 100 * w / sum(w)
  w
}


#' Evaluate the DR-Meta Variance Function
#'
#' Evaluates \eqn{\tau^2(\mathrm{DR};\,\psi)} for a grid of DR values or for
#' the study-level DR indices from a fitted \code{"drmeta"} model.
#' Useful for visualising how heterogeneity varies with design robustness.
#'
#' @param dr Numeric vector of design robustness values in \eqn{[0,1]}, or
#'   a fitted \code{"drmeta"} object (in which case `tau0sq`, `gamma`, and
#'   `vfun` are extracted automatically and `dr` from the model is used).
#' @param tau0sq Scalar \eqn{\hat\tau_0^2}.  Ignored if `dr` is a
#'   \code{"drmeta"} object.
#' @param gamma Scalar \eqn{\hat\gamma}.  Ignored if `dr` is a
#'   \code{"drmeta"} object.
#' @param vfun Variance function: `"exponential"` or `"linear"`.  Ignored if
#'   `dr` is a \code{"drmeta"} object.
#'
#' @return A numeric vector of \eqn{\tau^2(\mathrm{DR}_i)} values.
#'
#' @examples
#' # Evaluate on a grid
#' grid <- seq(0, 1, by = 0.1)
#' tau2 <- dr_variance(grid, tau0sq = 0.04, gamma = 1.5)
#' plot(grid, tau2, type = "l", xlab = "DR", ylab = expression(tau^2))
#'
#' # Extract from a fitted model
#' set.seed(1)
#' fit <- drmeta(yi = rnorm(10, 0.3), vi = runif(10, 0.01, 0.05),
#'               dr = runif(10))
#' dr_variance(fit)
#'
#' @seealso \code{\link{drmeta}}, \code{\link{dr_weights}}
#' @export
dr_variance <- function(dr, tau0sq = NULL, gamma = NULL,
                        vfun = c("exponential", "linear")) {

  if (inherits(dr, "drmeta")) {
    obj    <- dr
    tau0sq <- obj$tau0sq
    gamma  <- obj$gamma
    vfun   <- obj$vfun
    dr     <- obj$dr
  }

  vfun <- match.arg(vfun)

  if (is.null(tau0sq) || is.null(gamma))
    stop("Supply `tau0sq` and `gamma`, or pass a fitted `drmeta` object.",
         call. = FALSE)

  tau2_fn <- .make_tau2_fn(vfun)
  tau2_fn(dr, tau0sq, gamma)
}
