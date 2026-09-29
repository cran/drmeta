#' Plot the fitted variance function
#'
#' Draws the fitted residual between-study variance against the design index,
#' with study-level fitted variances overlaid. By default the curve is drawn
#' only over the observed range of the design index, since the fitted function
#' outside that range is an extrapolation.
#'
#' @param object A fitted \code{drmeta} object.
#' @param extrapolate If TRUE, draw the curve over the whole interval
#'   \eqn{[0, 1]} and shade the region beyond the observed support. Default
#'   FALSE.
#' @param col_curve Colour of the fitted curve.
#' @param col_pts Colour of the study-level points.
#' @param xlab,ylab,main Axis labels and title.
#' @param ... Further arguments passed to \code{\link[graphics]{plot}}.
#' @return Invisibly, the data frame of plotted grid values.
#' @examples
#' path <- system.file("extdata", "bcg_design_robustness.csv", package = "drmeta")
#' bcg <- utils::read.csv(path)
#' fit <- drmeta(yi = bcg[["yi"]], vi = bcg[["vi"]], dr = bcg[["dr"]])
#' dr_plot_vfun(fit)
#' @export
dr_plot_vfun <- function(object, extrapolate = FALSE,
                         col_curve = "#D6604D", col_pts = "#2166AC",
                         xlab = "Design robustness", ylab = "Fitted tau^2",
                         main = "Fitted variance function", ...) {
  if (!inherits(object, "drmeta")) stop("object must be a drmeta fit.")
  obs <- range(object$dr)
  grid_range <- if (isTRUE(extrapolate)) c(0, 1) else obs
  grid <- dr_scale_predict(object, seq(grid_range[1], grid_range[2], length.out = 200))
  ylim <- range(c(grid$tau2, object$tau2i), finite = TRUE)
  graphics::plot(grid$dr, grid$tau2, type = "n", xlab = xlab, ylab = ylab,
                 main = main, ylim = ylim, ...)
  if (isTRUE(extrapolate)) {
    usr <- graphics::par("usr")
    if (obs[1] > 0)
      graphics::rect(usr[1], usr[3], obs[1], usr[4],
                     col = grDevices::adjustcolor("grey70", alpha.f = 0.25), border = NA)
    if (obs[2] < 1)
      graphics::rect(obs[2], usr[3], usr[2], usr[4],
                     col = grDevices::adjustcolor("grey70", alpha.f = 0.25), border = NA)
  }
  graphics::lines(grid$dr, grid$tau2, col = col_curve, lwd = 2)
  graphics::points(object$dr, object$tau2i, col = col_pts, pch = 19)
  graphics::box()
  invisible(grid)
}
