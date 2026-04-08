# =============================================================================
# dr_heterogeneity.R  —  Heterogeneity diagnostics for DR-Meta
#
# Implements the design-residual decomposition (Proposition 6, Hait 2025):
#
#   tau^2_total = E[tau^2(DR_i)]  +  Var(u_i | DR_i)
#
# Also provides Cochran's Q, I^2, H^2, and per-study Q contributions,
# computed using DR-Meta weights (not standard inverse-variance weights).
# =============================================================================

#' Heterogeneity Diagnostics for DR-Meta
#'
#' Computes a suite of heterogeneity statistics for a fitted \code{"drmeta"}
#' model, including Cochran's Q (with DR-Meta weights), \eqn{I^2}, \eqn{H^2},
#' the design-residual variance decomposition of Proposition 6 (Hait, 2025),
#' and per-study contributions to \eqn{Q}.
#'
#' @section Design-Residual Decomposition (Proposition 6):
#' The total between-study heterogeneity can be decomposed as:
#' \deqn{\tau^2_{\text{total}} = \mathbb{E}[\tau^2(\mathrm{DR}_i)]
#'   + \mathrm{Var}(u_i \mid \mathrm{DR}_i),}
#' where the first term is the **design-explained** heterogeneity
#' (captured by the variance function) and the second is the
#' **design-residual** heterogeneity (unexplained by DR).  This decomposition
#' is analogous to R-squared in meta-regression.
#'
#' The design-explained proportion is
#' \deqn{R^2_{\mathrm{DR}} =
#'   \frac{\mathbb{E}[\tau^2(\mathrm{DR}_i)]}{\tau^2_{\mathrm{total}}}.}
#'
#' @param object A fitted \code{"drmeta"} object from \code{\link{drmeta}}.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{`summary`}{A one-row data frame with: `k`, `Q`, `df`, `pval`,
#'     `tau2_mean` (mean design-specific heterogeneity), `I2`, `H2`.}
#'   \item{`decomposition`}{A one-row data frame with the Proposition 6
#'     decomposition: `tau2_design_explained`, `tau2_residual`,
#'     `tau2_total`, `R2_DR` (proportion explained by design).}
#'   \item{`contributions`}{A data frame with per-study Q contributions:
#'     `study`, `DR`, `tau2_i`, `q_i`, `pct_Q`.}
#' }
#'
#' @references
#' Hait, S. (2025). *Design-Robust Meta-Analysis: A Variance-Function
#' Framework for Causal Credibility*. Proposition 6.
#'
#' @examples
#' set.seed(42)
#' k <- 15
#' dr <- runif(k, 0.1, 0.9)
#' vi <- runif(k, 0.01, 0.05)
#' tau2_true <- 0.04 * exp(-2 * dr)
#' yi <- rnorm(k, 0.3, sqrt(vi + tau2_true))
#'
#' fit <- drmeta(yi, vi, dr)
#' het <- dr_heterogeneity(fit)
#' het$summary
#' het$decomposition
#' het$contributions
#'
#' @export
dr_heterogeneity <- function(object) {

  if (!inherits(object, "drmeta"))
    stop("`object` must be a fitted `drmeta` model.", call. = FALSE)

  yi       <- object$yi
  vi       <- object$vi
  dr       <- object$dr
  tau2_i   <- object$tau2_i
  sigma2_i <- object$sigma2_i
  mu       <- object$mu
  k        <- object$k
  slab     <- object$slab

  # ---- Cochran's Q (with DR-Meta weights) --------------------------------
  # Q = sum( w_i * (y_i - mu)^2 ) where w_i = 1 / sigma2_i  (DR weights)
  w     <- 1 / sigma2_i
  resid <- yi - mu
  q_i   <- w * resid^2
  Q     <- sum(q_i)
  df_Q  <- k - 1L
  pval  <- stats::pchisq(Q, df = df_Q, lower.tail = FALSE)

  # ---- I^2 and H^2 -------------------------------------------------------
  # Use the typical heterogeneity formula adapted to DR-Meta sigma2
  # I^2 = (Q - df) / Q * 100
  I2 <- max(100 * (Q - df_Q) / Q, 0)
  H2 <- Q / df_Q

  # Mean design-specific tau^2 across studies
  tau2_mean <- mean(tau2_i)

  # ---- Proposition 6 decomposition ---------------------------------------
  # Design-explained: E[tau^2(DR_i)] = mean(tau2_i)
  tau2_design <- tau2_mean

  # Design-residual: estimated as the excess heterogeneity not captured by
  # the variance function.  We estimate this as the REML tau^2 from a
  # classical RE model on the DR-adjusted residuals.
  # Simpler approximation: total tau^2 from Q-based estimator on
  # residuals from mu_DR.
  C        <- sum(w) - sum(w^2) / sum(w)
  tau2_total <- max((Q - df_Q) / C, 0)
  tau2_resid <- max(tau2_total - tau2_design, 0)
  R2_DR      <- if (tau2_total > 0) tau2_design / tau2_total else NA_real_

  # ---- Summaries ---------------------------------------------------------
  summary_df <- data.frame(
    k         = k,
    Q         = round(Q, 4),
    df        = df_Q,
    pval      = round(pval, 4),
    tau2_mean = round(tau2_mean, 4),
    I2        = round(I2, 2),
    H2        = round(H2, 4),
    row.names  = NULL,
    check.names = FALSE
  )

  decomp_df <- data.frame(
    tau2_design_explained = round(tau2_design, 4),
    tau2_residual         = round(tau2_resid, 4),
    tau2_total            = round(tau2_total, 4),
    R2_DR                 = round(R2_DR, 4),
    row.names  = NULL,
    check.names = FALSE
  )

  contrib_df <- data.frame(
    study  = slab,
    DR     = round(dr, 3),
    tau2_i = round(tau2_i, 4),
    q_i    = round(q_i, 4),
    pct_Q  = round(100 * q_i / Q, 2),
    row.names  = NULL,
    check.names = FALSE
  )

  list(
    summary      = summary_df,
    decomposition = decomp_df,
    contributions = contrib_df
  )
}
