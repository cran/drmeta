#' drmeta: Design-Robust Meta-Analysis
#'
#' Fits the DR-Meta variance-function random-effects model (Hait, 2025),
#' where between-study heterogeneity is a monotone-decreasing function of
#' each study's design robustness index.
#'
#' @importFrom grDevices adjustcolor col2rgb rgb
#' @importFrom graphics abline legend lines mtext par points polygon rect text
#' @importFrom stats pnorm qnorm complete.cases median
#' @importFrom utils modifyList
#' @docType package
#' @name drmeta-package
"_PACKAGE"
