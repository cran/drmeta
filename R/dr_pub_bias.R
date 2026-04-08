# =============================================================================
# dr_pub_bias.R  —  Publication bias tools adapted for DR-Meta weights
#
# PET and PEESE use DR-Meta weights (1/sigma2_i) as the regression weights,
# so small-study effects are assessed conditional on design-adjusted precision.
# =============================================================================

#' Publication Bias Assessment for DR-Meta (PET/PEESE/Egger/Funnel)
#'
#' Performs a suite of publication-bias and small-study-effects tests adapted
#' for DR-Meta weights.  PET and PEESE regressions use
#' \eqn{w_i = 1/\hat\sigma_i^2 = 1/(v_i + \hat\tau^2(\mathrm{DR}_i))} as
#' regression weights, so precision is design-adjusted.
#'
#' @details
#' **PET** (Precision Effect Test): regresses \eqn{y_i} on \eqn{se_i =
#' \sqrt{v_i}}, with DR-Meta weights.  A significant slope implies small-study
#' bias; the intercept estimates the publication-bias-corrected effect.
#'
#' **PEESE** (Precision Effect Estimate with Standard Error): regresses
#' \eqn{y_i} on \eqn{v_i}, with DR-Meta weights.  Generally preferred when the
#' true effect is non-zero.
#'
#' **Egger test**: Egger-type regression using the standardised effect
#' \eqn{y_i / se_i} on \eqn{1 / se_i}, with DR-Meta weights.
#'
#' **Funnel asymmetry**: classic funnel plot with DR-Meta summary.
#'
#' @param object A fitted \code{"drmeta"} object from \code{\link{drmeta}}.
#' @param test Character vector of tests to run; any subset of
#'   `c("PET", "PEESE", "Egger")`.  Default: all three.
#' @param alpha Significance level for PET/PEESE intercept test.  Default 0.05.
#'
#' @return A list with elements:
#' \describe{
#'   \item{`PET`}{Summary of PET regression (data frame with intercept,
#'     slope, SE, z, p).}
#'   \item{`PEESE`}{Summary of PEESE regression.}
#'   \item{`Egger`}{Summary of Egger regression.}
#'   \item{`recommendation`}{Character string: use PET intercept if PET slope
#'     is significant and effect is small; use PEESE intercept otherwise.}
#' }
#'
#' @examples
#' set.seed(99)
#' k  <- 20
#' dr <- runif(k, 0.1, 0.9)
#' vi <- runif(k, 0.005, 0.08)
#' # Introduce small-study effect: smaller studies overestimate
#' yi <- rnorm(k, 0.3 + 0.5 * sqrt(vi), sqrt(vi + 0.03 * exp(-1.5 * dr)))
#' fit <- drmeta(yi, vi, dr)
#' pb  <- dr_pub_bias(fit)
#' pb$PET
#' pb$PEESE
#' pb$recommendation
#'
#' @export
dr_pub_bias <- function(object,
                        test  = c("PET", "PEESE", "Egger"),
                        alpha = 0.05) {

  if (!inherits(object, "drmeta"))
    stop("`object` must be a fitted `drmeta` model.", call. = FALSE)

  test <- match.arg(test, several.ok = TRUE)

  yi       <- object$yi
  vi       <- object$vi
  sei      <- sqrt(vi)
  w        <- 1 / object$sigma2_i   # DR-Meta weights
  k        <- object$k

  results <- list()

  # ---- PET: y_i ~ intercept + sei, weighted by w_i -----------------------
  if ("PET" %in% test) {
    pet <- .wls_summary(yi, cbind(1, sei), w, names = c("intercept", "sei"))
    results$PET <- pet
  }

  # ---- PEESE: y_i ~ intercept + vi, weighted by w_i ----------------------
  if ("PEESE" %in% test) {
    peese <- .wls_summary(yi, cbind(1, vi), w, names = c("intercept", "vi"))
    results$PEESE <- peese
  }

  # ---- Egger: (y_i/sei) ~ (1/sei) + intercept, weighted by w_i -----------
  if ("Egger" %in% test) {
    precision <- 1 / sei
    stand_yi  <- yi / sei
    egger     <- .wls_summary(stand_yi, cbind(precision, 1), w,
                               names = c("precision", "intercept"))
    results$Egger <- egger
  }

  # ---- Recommendation (PET-PEESE rule) ------------------------------------
  rec <- "Unable to evaluate: PET not run."
  if ("PET" %in% test && "PEESE" %in% test) {
    pet_slope_p <- results$PET["sei", "pval"]
    pet_intcpt  <- results$PET["intercept", "estimate"]
    if (!is.na(pet_slope_p) && pet_slope_p < alpha) {
      rec <- sprintf(
        paste0("PET slope is significant (p = %.4f): small-study bias likely.\n",
               "  Use PET intercept as bias-corrected estimate: %.4f.\n",
               "  If true effect is likely non-zero, prefer PEESE intercept: %.4f."),
        pet_slope_p,
        pet_intcpt,
        results$PEESE["intercept", "estimate"]
      )
    } else {
      rec <- sprintf(
        paste0("PET slope is not significant (p = %.4f): no strong evidence ",
               "of small-study bias.\n",
               "  Use main DR-Meta pooled estimate: %.4f."),
        if (is.na(pet_slope_p)) 1 else pet_slope_p,
        object$mu
      )
    }
  }
  results$recommendation <- rec

  results
}


#' Funnel Plot for DR-Meta
#'
#' Creates a funnel plot for a fitted \code{"drmeta"} model.  The horizontal
#' axis shows effect-size estimates; the vertical axis shows standard errors.
#' Point sizes are proportional to DR-Meta weights and point colour encodes
#' design robustness.
#'
#' @param object A fitted \code{"drmeta"} object.
#' @param contours Logical.  If `TRUE` (default), adds 95\% and 99\% funnel
#'   contours around the pooled estimate.
#' @param xlab X-axis label.
#' @param ylab Y-axis label (default: reversed SE axis).
#' @param main Plot title.
#' @param col_low Colour for low-DR studies.  Default `"#D6604D"`.
#' @param col_high Colour for high-DR studies.  Default `"#2166AC"`.
#'
#' @return Invisibly returns `NULL`.
#'
#' @examples
#' set.seed(42)
#' k <- 15
#' dr <- runif(k)
#' vi <- runif(k, 0.01, 0.06)
#' yi <- rnorm(k, 0.3, sqrt(vi + 0.04 * exp(-1.5 * dr)))
#' fit <- drmeta(yi, vi, dr)
#' dr_funnel(fit)
#'
#' @export
dr_funnel <- function(object,
                      contours = TRUE,
                      xlab     = "Effect size",
                      ylab     = "Standard error",
                      main     = "DR-Meta Funnel Plot",
                      col_low  = "#D6604D",
                      col_high = "#2166AC") {

  if (!inherits(object, "drmeta"))
    stop("`object` must be a fitted `drmeta` model.", call. = FALSE)

  yi   <- object$yi
  sei  <- sqrt(object$vi)
  dr   <- object$dr
  w    <- object$weights
  mu   <- object$mu

  # Colour by DR (gradient low→high)
  dr_norm <- (dr - min(dr)) / max(diff(range(dr)), 1e-6)
  cols    <- .blend_colours(col_low, col_high, dr_norm)

  # Point size proportional to weight
  cex_pts <- 0.6 + 2.0 * (w / max(w))

  sei_max <- max(sei) * 1.15
  ylim    <- c(sei_max, 0)  # y-axis reversed (larger SE at bottom)
  xlim    <- range(c(yi - 2 * sei, yi + 2 * sei, mu + c(-1, 1) * 0.1))

  plot(yi, sei,
       pch  = 21,
       bg   = cols,
       col  = "grey40",
       cex  = cex_pts,
       xlim = xlim,
       ylim = ylim,
       xlab = xlab,
       ylab = ylab,
       main = main)

  abline(v = mu, lty = 2, col = "grey50")

  if (contours) {
    sei_seq <- seq(0, sei_max, length.out = 200)
    lines(mu + qnorm(0.975) * sei_seq, sei_seq, col = "grey60", lty = 3)
    lines(mu - qnorm(0.975) * sei_seq, sei_seq, col = "grey60", lty = 3)
    lines(mu + qnorm(0.995) * sei_seq, sei_seq, col = "grey80", lty = 3)
    lines(mu - qnorm(0.995) * sei_seq, sei_seq, col = "grey80", lty = 3)
  }

  # Colour legend for DR
  legend("topright",
         legend = c("High DR", "Low DR"),
         pch    = 21,
         pt.bg  = c(col_high, col_low),
         col    = "grey40",
         bty    = "n", cex = 0.85)

  invisible(NULL)
}


# =============================================================================
# Internal WLS helper
# =============================================================================

.wls_summary <- function(y, X, w, names = NULL) {
  # Weighted least squares: solve (X'WX)^{-1} X'Wy
  W     <- diag(w)
  XtW   <- t(X) %*% W
  XtWX  <- XtW %*% X
  XtWy  <- XtW %*% y
  beta  <- tryCatch(
    solve(XtWX, XtWy),
    error = function(e) rep(NA_real_, ncol(X))
  )
  resid <- y - X %*% beta
  # Residual variance (heteroskedasticity-consistent approximation)
  sigma2_hat <- sum(w * resid^2) / (length(y) - ncol(X))
  vcov_b <- tryCatch(
    solve(XtWX) * sigma2_hat,
    error = function(e) matrix(NA_real_, ncol(X), ncol(X))
  )
  se   <- sqrt(diag(vcov_b))
  zval <- beta / se
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)

  out <- data.frame(
    estimate = as.numeric(beta),
    se       = se,
    zval     = zval,
    pval     = pval,
    row.names = if (is.null(names)) paste0("b", seq_along(beta)) else names
  )
  out
}

# Blend two hex colours by fraction t in [0,1]
.blend_colours <- function(col0, col1, t) {
  rgb0 <- col2rgb(col0) / 255
  rgb1 <- col2rgb(col1) / 255
  r <- outer(1 - t, rgb0[1, ]) + outer(t, rgb1[1, ])
  g <- outer(1 - t, rgb0[2, ]) + outer(t, rgb1[2, ])
  b <- outer(1 - t, rgb0[3, ]) + outer(t, rgb1[3, ])
  rgb(r, g, b)
}
