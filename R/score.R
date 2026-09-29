#' Rescale a numeric vector to the unit interval
#'
#' @param x A numeric vector. Missing values are ignored during rescaling and
#'   preserved in the result.
#' @return A numeric vector of the same length rescaled to \eqn{[0, 1]}. If all
#'   non-missing values are equal, a zero vector is returned.
#' @examples
#' normalize_01(c(2, 5, 8))
#' normalize_01(c(1, 1, 1))
#' @export
normalize_01 <- function(x) {
  x <- as.numeric(x)
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng))) stop("x must contain at least one finite value.")
  if (diff(rng) == 0) return(ifelse(is.na(x), NA_real_, 0))
  (x - rng[1]) / diff(rng)
}

#' Build a design-robustness index from sub-scores
#'
#' Combines user-supplied sub-score dimensions into a single index in
#' \eqn{[0, 1]} as a normalised weighted average. The components, their
#' coding rules, and their weights should be prespecified and reported.
#'
#' For confirmatory use, none of the sub-scores may be a function of the
#' realized effect estimate, its standard error, its p-value, its confidence
#' interval, or any outcome-dependent diagnostic. A score that uses such
#' information defines an exploratory analysis, not a primary scale moderator.
#' This function cannot verify that condition; it is the analyst's
#' responsibility.
#'
#' @param ... Named numeric vectors of equal length, each a sub-score in
#'   \eqn{[0, 1]}.
#' @param weights Optional relative weights, one per sub-score. Defaults to
#'   equal weighting.
#' @param warn_range If TRUE, warn when any value falls outside \eqn{[0, 1]}
#'   before clipping.
#' @return A numeric vector in \eqn{[0, 1]}, carrying an attribute
#'   \code{"subscores"} that holds the clipped inputs and the resulting index.
#' @examples
#' dr_score(balance = c(0.9, 0.6, 0.4),
#'          overlap = c(0.8, 0.7, 0.5),
#'          weights = c(2, 1))
#' @export
dr_score <- function(..., weights = NULL, warn_range = TRUE) {
  parts <- list(...)
  if (!length(parts)) stop("Supply at least one sub-score.")
  if (is.null(names(parts)) || any(!nzchar(names(parts))))
    stop("Every sub-score must be named.")
  lens <- vapply(parts, length, integer(1))
  if (length(unique(lens)) != 1L) stop("All sub-scores must have equal length.")
  mat <- vapply(parts, as.numeric, numeric(lens[1]))
  if (!is.matrix(mat)) mat <- matrix(mat, nrow = lens[1],
                                     dimnames = list(NULL, names(parts)))
  if (any(!is.finite(mat))) stop("Sub-scores must be finite.")
  if (warn_range && any(mat < 0 | mat > 1))
    warning("Some sub-scores fall outside [0,1] and were clipped.")
  mat <- pmin(pmax(mat, 0), 1)
  w <- if (is.null(weights)) rep(1, ncol(mat)) else as.numeric(weights)
  if (length(w) != ncol(mat)) stop("weights must have one entry per sub-score.")
  if (any(!is.finite(w)) || any(w < 0) || sum(w) <= 0)
    stop("weights must be nonnegative, finite, and not all zero.")
  w <- w / sum(w)
  dr <- as.vector(mat %*% w)
  attr(dr, "subscores") <- data.frame(mat, dr = dr)
  dr
}

#' Map study design labels to illustrative design-robustness scores
#'
#' Provides a convenience mapping from design-type labels to numeric scores.
#'
#' This mapping is a demonstration default, not a validated measurement model
#' and not a universal hierarchy of study quality. A design label does not
#' determine the magnitude or direction of bias in a particular evidence base.
#' For substantive work the index should be constructed for the specific
#' synthesis, outcome, and identification problem at hand, ordinarily with
#' \code{\link{dr_score}}, and reported under alternative defensible codings.
#'
#' @param design Character vector of design labels (case-insensitive).
#' @param custom_map Optional named numeric vector overriding or extending the
#'   defaults.
#' @param default_score Score assigned to unrecognised labels. Default 0.25.
#' @param warn_unknown If TRUE, warn about unrecognised labels.
#' @return A numeric vector of scores in \eqn{[0, 1]}, the same length as
#'   \code{design}.
#' @seealso \code{\link{dr_score}}
#' @examples
#' dr_from_design(c("RCT", "DiD", "OLS", "IV", "matching"))
#' dr_from_design("house_method", custom_map = c(house_method = 0.65))
#' @export
dr_from_design <- function(design, custom_map = NULL, default_score = 0.25,
                           warn_unknown = TRUE) {
  base_map <- c(rct = 1.00, rd = 0.75, rdd = 0.75, iv = 0.75,
                did = 0.60, diff_in_diff = 0.60, matching = 0.50, psm = 0.50,
                ipw = 0.45, propensity = 0.45, regression = 0.25, ols = 0.25,
                cross_section = 0.20, descriptive = 0.10)
  if (!is.null(custom_map)) {
    if (is.null(names(custom_map)) || any(!nzchar(names(custom_map))))
      stop("custom_map must be a named numeric vector.")
    cm <- as.numeric(custom_map)
    if (any(!is.finite(cm)) || any(cm < 0 | cm > 1))
      stop("custom_map values must be finite and lie in [0,1].")
    names(cm) <- tolower(names(custom_map))
    base_map[names(cm)] <- cm
  }
  key <- tolower(trimws(as.character(design)))
  out <- unname(base_map[key])
  unknown <- is.na(out)
  if (any(unknown)) {
    if (warn_unknown)
      warning("Unrecognised design labels assigned the default score: ",
              paste(unique(key[unknown]), collapse = ", "))
    out[unknown] <- default_score
  }
  out
}
