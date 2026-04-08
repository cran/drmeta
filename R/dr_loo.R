# =============================================================================
# dr_loo.R  —  Leave-one-out influence diagnostics for DR-Meta
# =============================================================================

#' Leave-One-Out Influence Diagnostics for DR-Meta
#'
#' For each study, re-fits the DR-Meta model after excluding that study and
#' records how the pooled estimate, confidence interval, variance-function
#' parameters, and heterogeneity change.  Studies with large absolute
#' \eqn{\Delta\hat\mu} or large shifts in \eqn{\hat\tau_0^2}/\eqn{\hat\gamma}
#' are considered influential.
#'
#' @param object A fitted \code{"drmeta"} object from \code{\link{drmeta}}.
#' @param parallel Logical.  If `TRUE` and the \pkg{parallel} package is
#'   available, uses \code{parallel::mclapply} for the LOO loop (Unix/macOS
#'   only).  Default `FALSE`.
#' @param mc.cores Integer.  Number of cores for parallel execution.  Default
#'   is \code{parallel::detectCores() - 1}.
#'
#' @return A list with:
#' \describe{
#'   \item{`summary`}{A data frame with one row per study and columns:
#'     `study`, `DR`, `est_loo` (LOO pooled estimate), `ci.lb_loo`,
#'     `ci.ub_loo`, `tau0sq_loo`, `gamma_loo`, `delta_mu` (change in
#'     estimate), `delta_tau0sq`, `delta_gamma`, `influential` (logical:
#'     |delta_mu| > 2 * SE of full model).}
#'   \item{`full`}{The original full-model \code{"drmeta"} object.}
#' }
#'
#' @examples
#' set.seed(7)
#' k  <- 12
#' dr <- runif(k, 0.1, 0.9)
#' vi <- runif(k, 0.01, 0.05)
#' tau2_true <- 0.04 * exp(-2 * dr)
#' yi <- rnorm(k, 0.3, sqrt(vi + tau2_true))
#'
#' fit <- drmeta(yi, vi, dr)
#' loo <- dr_loo(fit)
#' loo$summary
#'
#' @export
dr_loo <- function(object, parallel = FALSE, mc.cores = NULL) {

  if (!inherits(object, "drmeta"))
    stop("`object` must be a fitted `drmeta` model.", call. = FALSE)

  yi     <- object$yi
  vi     <- object$vi
  dr     <- object$dr
  slab   <- object$slab
  k      <- object$k
  vfun   <- object$vfun
  method <- object$method

  full_mu     <- object$mu
  full_se     <- object$se
  full_tau0sq <- object$tau0sq
  full_gamma  <- object$gamma

  # Influence threshold: |delta_mu| > 2 * SE_full
  thresh <- 2 * full_se

  # ---- LOO loop -----------------------------------------------------------
  .fit_loo <- function(i) {
    yi_i  <- yi[-i];  vi_i <- vi[-i];  dr_i <- dr[-i]
    fit_i <- tryCatch(
      drmeta(yi = yi_i, vi = vi_i, dr = dr_i,
             vfun = vfun, method = method),
      error = function(e) NULL
    )
    if (is.null(fit_i)) {
      return(data.frame(
        study       = slab[i],
        DR          = dr[i],
        est_loo     = NA_real_,
        ci.lb_loo   = NA_real_,
        ci.ub_loo   = NA_real_,
        tau0sq_loo  = NA_real_,
        gamma_loo   = NA_real_,
        delta_mu    = NA_real_,
        delta_tau0sq = NA_real_,
        delta_gamma = NA_real_,
        influential = NA
      ))
    }
    data.frame(
      study       = slab[i],
      DR          = dr[i],
      est_loo     = fit_i$mu,
      ci.lb_loo   = fit_i$ci.lb,
      ci.ub_loo   = fit_i$ci.ub,
      tau0sq_loo  = fit_i$tau0sq,
      gamma_loo   = fit_i$gamma,
      delta_mu    = fit_i$mu - full_mu,
      delta_tau0sq = fit_i$tau0sq - full_tau0sq,
      delta_gamma = fit_i$gamma - full_gamma,
      influential = abs(fit_i$mu - full_mu) > thresh
    )
  }

  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    cores <- mc.cores %||% max(1L, parallel::detectCores() - 1L)
    rows  <- parallel::mclapply(seq_len(k), .fit_loo, mc.cores = cores)
  } else {
    rows <- lapply(seq_len(k), .fit_loo)
  }

  summary_df <- do.call(rbind, rows)
  rownames(summary_df) <- NULL

  attr(summary_df, "full_mu")     <- full_mu
  attr(summary_df, "full_se")     <- full_se
  attr(summary_df, "full_tau0sq") <- full_tau0sq
  attr(summary_df, "full_gamma")  <- full_gamma
  attr(summary_df, "threshold")   <- thresh

  list(summary = summary_df, full = object)
}

# internal null-coalescing helper (also used in other files)
`%||%` <- function(x, y) if (is.null(x)) y else x
