# =============================================================================
# dr_plot.R  —  Visualisation functions for DR-Meta
# =============================================================================

#' Forest Plot for DR-Meta
#'
#' Draws a forest plot for a fitted \code{"drmeta"} model using base graphics.
#' Studies are ordered by design robustness (largest DR at top by default).
#' Point sizes are proportional to DR-Meta weights; a vertical reference line
#' and summary diamond are included.
#'
#' @param object A fitted \code{"drmeta"} object.
#' @param order_by Character: `"dr"` (default, sort by DR_i descending),
#'   `"yi"` (sort by effect size), or `"none"` (original order).
#' @param xlab X-axis label.  Default `"Effect size"`.
#' @param main Plot title.
#' @param col_point Colour for study-level estimate points.  Default `"#2166AC"`.
#' @param col_diamond Colour for the summary diamond.  Default `"#D6604D"`.
#' @param col_dr Colour for the DR bar on the left.  Default `"#4DAC26"`.
#' @param show_dr Logical.  If `TRUE` (default), displays a coloured DR bar
#'   indicating design robustness strength.
#' @param xlim Numeric vector of length 2 for x-axis limits.  If `NULL`,
#'   computed automatically.
#' @param cex_study Scaling factor for study labels.  Default 0.8.
#' @param ... Further graphical arguments (ignored).
#'
#' @return Invisibly returns the data frame used for plotting (ordered studies).
#'
#' @examples
#' set.seed(42)
#' k <- 10
#' dr <- runif(k, 0.1, 0.9)
#' vi <- runif(k, 0.01, 0.05)
#' yi <- rnorm(k, 0.3, sqrt(vi + 0.04 * exp(-2 * dr)))
#' fit <- drmeta(yi, vi, dr)
#' dr_forest(fit)
#'
#' @export
dr_forest <- function(object,
                      order_by   = c("dr", "yi", "none"),
                      xlab       = "Effect size",
                      main       = "DR-Meta Forest Plot",
                      col_point  = "#2166AC",
                      col_diamond = "#D6604D",
                      col_dr     = "#4DAC26",
                      show_dr    = TRUE,
                      xlim       = NULL,
                      cex_study  = 0.8,
                      ...) {

  if (!inherits(object, "drmeta"))
    stop("`object` must be a fitted `drmeta` model.", call. = FALSE)

  order_by <- match.arg(order_by)

  # ---- Build plotting data frame -----------------------------------------
  ci.lb_i <- object$yi - qnorm(0.975) * sqrt(object$vi)
  ci.ub_i <- object$yi + qnorm(0.975) * sqrt(object$vi)

  df <- data.frame(
    study  = object$slab,
    yi     = object$yi,
    vi     = object$vi,
    dr     = object$dr,
    ci.lb  = ci.lb_i,
    ci.ub  = ci.ub_i,
    weight = object$weights * 100,
    stringsAsFactors = FALSE
  )

  if (order_by == "dr")   df <- df[order(df$dr,  decreasing = FALSE), ]
  if (order_by == "yi")   df <- df[order(df$yi,  decreasing = FALSE), ]
  k <- nrow(df)

  if (is.null(xlim)) {
    pad  <- 0.1 * diff(range(c(df$ci.lb, df$ci.ub, object$ci.lb, object$ci.ub)))
    xlim <- range(c(df$ci.lb, df$ci.ub, object$ci.lb, object$ci.ub)) + c(-pad, pad)
  }

  # ---- Layout ------------------------------------------------------------
  left_margin <- if (show_dr) 7.5 else 6
  op <- par(mar = c(4, left_margin, 3, 2), no.readonly = TRUE)
  on.exit(par(op))

  y_pos <- seq_len(k) + 2  # leave 2 rows at bottom for diamond

  plot(NA, xlim = xlim, ylim = c(0, max(y_pos) + 1),
       xlab = xlab, ylab = "", main = main,
       yaxt = "n", bty = "l")

  abline(v = 0, lty = 2, col = "grey60")

  # ---- Study-level rows --------------------------------------------------
  cex_pts <- 0.4 + 1.6 * (df$weight / max(df$weight))

  for (i in seq_len(k)) {
    y <- y_pos[i]
    # CI line
    lines(c(df$ci.lb[i], df$ci.ub[i]), c(y, y), col = "grey40", lwd = 1.2)
    # Point (size proportional to weight)
    points(df$yi[i], y, pch = 15, col = col_point, cex = cex_pts[i])
    # Study label
    mtext(df$study[i], side = 2, at = y, las = 1, cex = cex_study,
          adj = 1, line = 0.5)

    # DR bar
    if (show_dr) {
      dr_col <- adjustcolor(col_dr, alpha.f = 0.3 + 0.7 * df$dr[i])
      rect(xlim[1] - 0.02 * diff(xlim), y - 0.35,
           xlim[1] - 0.002 * diff(xlim), y + 0.35,
           col = dr_col, border = NA)
    }
  }

  # ---- Summary diamond ---------------------------------------------------
  y_diam <- 1.2
  mu     <- object$mu
  ci_lb  <- object$ci.lb
  ci_ub  <- object$ci.ub

  polygon(
    x = c(ci_lb, mu, ci_ub, mu),
    y = c(y_diam, y_diam + 0.4, y_diam, y_diam - 0.4),
    col = col_diamond, border = col_diamond
  )

  # Summary label
  mtext(sprintf("Pooled (DR-Meta): %.3f [%.3f, %.3f]",
                mu, ci_lb, ci_ub),
        side = 1, line = 2.5, cex = 0.85, adj = 0)

  # Separator line above diamond
  abline(h = y_pos[1] - 0.5, col = "grey60", lwd = 0.8)

  if (show_dr) {
    mtext("DR", side = 2, at = max(y_pos) + 0.7, las = 1,
          cex = 0.7, col = col_dr, line = 0.5)
  }

  invisible(df)
}


#' Weight Diagnostic Plot for DR-Meta
#'
#' Plots DR-Meta study weights against design robustness (\eqn{\mathrm{DR}_i}),
#' illustrating the monotone ordering of Lemma 3 (Hait, 2025).  An overlaid
#' curve shows the theoretical weight function holding sampling variance at
#' its median value.
#'
#' @param object A fitted \code{"drmeta"} object.
#' @param col_pts Colour for individual study points.  Default `"#2166AC"`.
#' @param col_curve Colour for the theoretical weight curve.  Default `"#D6604D"`.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param main Plot title.
#' @param show_labels Logical.  If `TRUE` (default), labels outlier points.
#' @param ... Additional graphical arguments passed to `plot`.
#'
#' @return Invisibly returns a data frame with study-level weight information.
#'
#' @examples
#' set.seed(42)
#' k <- 12
#' dr <- runif(k, 0.1, 0.9)
#' vi <- runif(k, 0.01, 0.05)
#' yi <- rnorm(k, 0.3, sqrt(vi + 0.04 * exp(-2 * dr)))
#' fit <- drmeta(yi, vi, dr)
#' dr_plot(fit)
#'
#' @export
dr_plot <- function(object,
                    col_pts    = "#2166AC",
                    col_curve  = "#D6604D",
                    xlab       = expression(paste("Design robustness (", DR[i], ")")),
                    ylab       = "DR-Meta weight (normalised, %)",
                    main       = "DR-Meta Weight vs Design Robustness",
                    show_labels = TRUE,
                    ...) {

  if (!inherits(object, "drmeta"))
    stop("`object` must be a fitted `drmeta` model.", call. = FALSE)

  w_pct <- object$weights * 100
  dr    <- object$dr
  vi    <- object$vi
  slab  <- object$slab

  tau2_fn  <- .make_tau2_fn(object$vfun)
  tau0sq   <- object$tau0sq
  gamma    <- object$gamma
  vi_med   <- stats::median(vi)

  # Theoretical curve at median vi
  dr_grid   <- seq(0, 1, length.out = 200)
  tau2_grid <- tau2_fn(dr_grid, tau0sq, gamma)
  w_grid    <- 1 / (vi_med + tau2_grid)
  # Normalise curve relative to study weights for comparability
  w_grid_n  <- w_grid / sum(1 / (vi + tau2_fn(dr, tau0sq, gamma))) * 100

  ylim <- c(0, max(w_pct, w_grid_n) * 1.1)

  plot(dr, w_pct,
       pch = 16, col = col_pts,
       xlab = xlab, ylab = ylab, main = main,
       xlim = c(0, 1), ylim = ylim, ...)

  lines(dr_grid, w_grid_n, col = col_curve, lwd = 2, lty = 2)

  legend("topleft",
         legend = c("Study weights", "Theoretical curve (median vi)"),
         pch = c(16, NA), lty = c(NA, 2),
         col = c(col_pts, col_curve), bty = "n", cex = 0.85)

  if (show_labels) {
    thresh <- mean(w_pct) + 1.5 * stats::sd(w_pct)
    idx    <- which(w_pct > thresh)
    if (length(idx) > 0)
      text(dr[idx], w_pct[idx], labels = slab[idx],
           pos = 3, cex = 0.7, col = col_pts)
  }

  invisible(data.frame(study = slab, dr = dr, weight_pct = w_pct))
}


#' Variance Function Plot for DR-Meta
#'
#' Plots the fitted variance function \eqn{\tau^2(\mathrm{DR};\,\hat\psi)}
#' as a curve from DR = 0 to DR = 1, with study-level \eqn{\hat\tau^2_i}
#' overlaid as points.
#'
#' @param object A fitted \code{"drmeta"} object.
#' @param col_curve Line colour for variance function.  Default `"#D6604D"`.
#' @param col_pts Point colour for study-level \eqn{\hat\tau^2_i}.
#'   Default `"#2166AC"`.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param main Plot title.
#'
#' @return Invisibly returns a data frame with the plotting grid.
#'
#' @examples
#' set.seed(1)
#' fit <- drmeta(rnorm(12, 0.3), runif(12, 0.01, 0.04), runif(12))
#' dr_plot_vfun(fit)
#'
#' @export
dr_plot_vfun <- function(object,
                         col_curve = "#D6604D",
                         col_pts   = "#2166AC",
                         xlab = expression(paste("Design robustness (", DR[i], ")")),
                         ylab = expression(paste("Between-study variance ", tau^2)),
                         main = "DR-Meta Variance Function") {

  if (!inherits(object, "drmeta"))
    stop("`object` must be a fitted `drmeta` model.", call. = FALSE)

  tau2_fn <- .make_tau2_fn(object$vfun)
  tau0sq  <- object$tau0sq
  gamma   <- object$gamma

  dr_grid   <- seq(0, 1, length.out = 300)
  tau2_grid <- tau2_fn(dr_grid, tau0sq, gamma)

  ylim <- c(0, max(tau2_grid, object$tau2_i) * 1.15)

  plot(dr_grid, tau2_grid,
       type = "l", col = col_curve, lwd = 2.5,
       xlab = xlab, ylab = ylab, main = main,
       ylim = ylim, xlim = c(0, 1))

  points(object$dr, object$tau2_i, pch = 16, col = col_pts, cex = 1.1)

  legend("topright",
         legend = c(
           sprintf("Fitted: %s  tau0^2=%.4f  gamma=%.3f", object$vfun, tau0sq, gamma),
           expression(paste("Study-level ", hat(tau)[i]^2))
         ),
         lty = c(1, NA), pch = c(NA, 16),
         col = c(col_curve, col_pts), bty = "n", cex = 0.85)

  abline(h = tau0sq, lty = 3, col = "grey50")
  text(0.02, tau0sq * 1.04,
       labels = expression(paste(hat(tau)[0]^2, " (DR=0)")),
       cex = 0.75, col = "grey40", adj = 0)

  invisible(data.frame(dr = dr_grid, tau2 = tau2_grid))
}
