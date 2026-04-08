# =============================================================================
# dr_score.R  —  Construct the Design Robustness Index DR_i
#
# The DR_i index summarises a study's causal identification strength on [0,1].
# Larger values indicate stronger design robustness and lower susceptibility
# to bias (Assumption A2; Section 3.1 of Hait, 2025).
#
# Two approaches are provided:
#   1. dr_score()   — weighted composite from continuous sub-scores
#   2. dr_from_design() — hierarchy-based look-up from design type labels
# =============================================================================

#' Construct a Design Robustness Index (DR_i)
#'
#' Computes a scalar design robustness index \eqn{\mathrm{DR}_i \in [0,1]}
#' for each study by forming a weighted composite of user-supplied sub-scores.
#' This is the recommended way to operationalise the DR index described in
#' Section 3.1 of Hait (2025).
#'
#' @details
#' Sub-scores are first clipped to \eqn{[0,1]} and then combined as a
#' normalised weighted average:
#' \deqn{\mathrm{DR}_i = \sum_j \tilde{w}_j\, s_{ij},
#'   \qquad \tilde{w}_j = w_j / \sum_j w_j.}
#' The result is therefore guaranteed to lie in \eqn{[0,1]}.
#'
#' Typical sub-score dimensions for quasi-experimental studies include:
#' \itemize{
#'   \item **Balance**: covariate balance between treatment and control (e.g.,
#'     standardised mean difference < 0.1 scores 1.0).
#'   \item **Overlap**: common-support / propensity-score overlap.
#'   \item **Design**: study design type — see \code{\link{dr_from_design}}.
#'   \item **Transparency**: pre-registration, data/code availability.
#' }
#'
#' @param ... Named numeric vectors, each of length \eqn{k} (number of
#'   studies), representing individual sub-score dimensions.  Each element
#'   must lie in \eqn{[0, 1]}.  Names are used in the returned data frame.
#' @param weights Optional numeric vector of the same length as the number of
#'   sub-score arguments, giving the relative importance of each dimension.
#'   Defaults to equal weighting.
#' @param warn_range Logical.  If `TRUE` (default), warns when any input value
#'   lies outside \eqn{[0,1]} before clipping.
#'
#' @return A numeric vector of length \eqn{k} with values in \eqn{[0,1]}.
#'   The vector carries an attribute `"subscores"` containing a data frame of
#'   the clipped sub-scores and the final DR index.
#'
#' @examples
#' k <- 5
#' balance  <- c(0.9, 0.6, 0.4, 0.8, 0.3)
#' overlap  <- c(0.8, 0.7, 0.5, 0.9, 0.4)
#' design   <- c(1.0, 0.5, 0.5, 0.75, 0.25)
#'
#' # Equal weights
#' dr <- dr_score(balance = balance, overlap = overlap, design = design)
#' dr
#'
#' # Down-weight transparency
#' transp <- c(1, 0, 0, 1, 0)
#' dr_w <- dr_score(balance = balance, overlap = overlap,
#'                  design = design, transparency = transp,
#'                  weights = c(2, 2, 3, 1))
#' dr_w
#'
#' @seealso \code{\link{dr_from_design}}, \code{\link{drmeta}}
#' @export
dr_score <- function(..., weights = NULL, warn_range = TRUE) {

  subs <- list(...)

  if (length(subs) == 0L)
    stop("Supply at least one sub-score vector.", call. = FALSE)

  # Names
  nms <- names(subs)
  if (is.null(nms) || any(nms == ""))
    nms <- paste0("score", seq_along(subs))

  # Check all vectors are numeric and same length
  if (!all(vapply(subs, is.numeric, logical(1L))))
    stop("All sub-score arguments must be numeric vectors.", call. = FALSE)

  lens <- vapply(subs, length, integer(1L))
  if (length(unique(lens)) > 1L)
    stop("All sub-score vectors must have the same length (k).", call. = FALSE)
  k <- lens[1L]

  # Warn and clip to [0,1]
  subs_mat <- do.call(cbind, subs)
  if (warn_range && any(subs_mat < 0 | subs_mat > 1, na.rm = TRUE))
    warning("Some sub-scores lie outside [0,1] and will be clipped.",
            call. = FALSE)
  subs_mat <- pmax(pmin(subs_mat, 1), 0)

  # Weights
  p <- length(subs)
  if (is.null(weights)) {
    weights <- rep(1, p)
  } else {
    if (!is.numeric(weights) || length(weights) != p)
      stop("`weights` must be a numeric vector of length equal to the ",
           "number of sub-scores (", p, ").", call. = FALSE)
    if (any(weights < 0))
      stop("`weights` must be non-negative.", call. = FALSE)
  }
  w_norm <- weights / sum(weights)

  dr <- as.numeric(subs_mat %*% w_norm)

  # Attach sub-score data frame as attribute
  ss_df <- as.data.frame(subs_mat)
  colnames(ss_df) <- nms
  ss_df$DR <- dr
  attr(dr, "subscores") <- ss_df

  dr
}


#' Design Robustness from Study Design Type
#'
#' Maps a vector of study design type labels to a numeric design robustness
#' score in \eqn{[0,1]}, using a pre-specified hierarchy of causal credibility.
#' This is a convenient starting point for operationalising the DR index when
#' only design type is available.
#'
#' @details
#' The default hierarchy follows the causal inference literature (Rubin, 2008;
#' Rosenbaum, 2010; Imbens & Rubin, 2015):
#'
#' | Design type label | Default DR score |
#' |---|---|
#' | `"rct"` | 1.00 |
#' | `"rd"`, `"rdd"` | 0.75 |
#' | `"iv"` | 0.75 |
#' | `"did"`, `"diff_in_diff"` | 0.60 |
#' | `"matching"`, `"psm"` | 0.50 |
#' | `"ipw"`, `"propensity"` | 0.45 |
#' | `"regression"`, `"ols"` | 0.25 |
#' | `"cross_section"` | 0.20 |
#' | `"descriptive"` | 0.10 |
#'
#' Users can override or extend this table via the `custom_map` argument.
#'
#' @param design Character vector of design type labels (case-insensitive).
#' @param custom_map Optional named numeric vector to override or add design
#'   types, e.g. `c(my_design = 0.65)`.
#' @param default_score Numeric score assigned to unrecognised design labels.
#'   Default is 0.25 (conservative).
#' @param warn_unknown Logical.  If `TRUE` (default), warns about unrecognised
#'   labels.
#'
#' @return A numeric vector of design robustness scores in \eqn{[0,1]},
#'   the same length as `design`.
#'
#' @examples
#' designs <- c("RCT", "DiD", "OLS", "IV", "matching", "unknown_design")
#' dr_from_design(designs)
#'
#' # Custom override
#' dr_from_design(designs, custom_map = c(unknown_design = 0.35))
#'
#' @seealso \code{\link{dr_score}}, \code{\link{drmeta}}
#' @export
dr_from_design <- function(design,
                           custom_map    = NULL,
                           default_score = 0.25,
                           warn_unknown  = TRUE) {

  if (!is.character(design))
    stop("`design` must be a character vector.", call. = FALSE)

  # Default hierarchy
  base_map <- c(
    rct          = 1.00,
    rd           = 0.75,
    rdd          = 0.75,
    iv           = 0.75,
    did          = 0.60,
    diff_in_diff = 0.60,
    matching     = 0.50,
    psm          = 0.50,
    ipw          = 0.45,
    propensity   = 0.45,
    regression   = 0.25,
    ols          = 0.25,
    cross_section = 0.20,
    descriptive  = 0.10
  )

  if (!is.null(custom_map)) {
    if (!is.numeric(custom_map) || is.null(names(custom_map)))
      stop("`custom_map` must be a named numeric vector.", call. = FALSE)
    names(custom_map) <- tolower(trimws(names(custom_map)))
    base_map[names(custom_map)] <- custom_map
  }

  design_lc <- tolower(trimws(design))
  known     <- design_lc %in% names(base_map)

  if (warn_unknown && any(!known)) {
    unk <- unique(design[!known])
    warning("Unrecognised design types set to default_score = ", default_score,
            ": ", paste(unk, collapse = ", "), call. = FALSE)
  }

  scores <- ifelse(known, base_map[design_lc], default_score)
  unname(scores)
}


#' Normalise a Numeric Vector to the \eqn{[0, 1]} Interval
#'
#' Linearly rescales a numeric vector to the \eqn{[0,1]} interval. Useful for
#' standardising individual sub-score components before aggregation with
#' \code{\link{dr_score}}.
#'
#' @param x A numeric vector.  `NA` values are ignored during rescaling.
#' @return A numeric vector rescaled to \eqn{[0,1]}.  If all non-missing
#'   values are equal, returns a zero vector (to avoid division by zero).
#'
#' @examples
#' normalize_01(c(2, 5, 8))   # returns c(0, 0.5, 1)
#' normalize_01(c(1, 1, 1))   # returns c(0, 0, 0)
#'
#' @export
normalize_01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}
