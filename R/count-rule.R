# R/count-rule.R
# -----------------------------------------------------------------------------
# The record of how a cell count was set
# -----------------------------------------------------------------------------
# select_resolution() writes one list, `count_rule`, when the count is set on
# the path that names an estimate: the profile carries one, or the criterion
# is "size" or "loo".  It says which rule set the count, for which estimate,
# and what the profile scored at that count.  A selection from a profile
# without an estimate by "cp", "reliability", "elbow" or "moran_z" carries
# none, and neither does a count given as a number.


#' The rows a size divides, read off a profile's bounds
#'
#' The rows with a response when the profile had one and some row carried
#' it, every point otherwise: the number \code{resolution_profile()} divides
#' by \code{obs_per_cell}.
#'
#' @param bounds \code{attr(<profile>, "bounds")}, or \code{NULL} for an
#'   object that has lost it (a column subset).
#' @return One integer; \code{NA} when \code{bounds} holds no count.
#' @keywords internal
#' @noRd
.count_rule_rows <- function(bounds) {
  if (!is.list(bounds)) return(NA_integer_)
  one <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v)
  if (one(bounds$n_response) && bounds$n_response > 0) return(as.integer(bounds$n_response))
  if (one(bounds$n) && bounds$n > 0) return(as.integer(bounds$n))
  NA_integer_
}


#' A profile's value of one column at one level
#'
#' @param profile A \code{resolution_profile}, or a subset of one.
#' @param level The level (cell count).
#' @param column Name of the column.
#' @return One number; \code{NA} when the profile has no such column or no
#'   such level.
#' @keywords internal
#' @noRd
.profile_at <- function(profile, level, column) {
  if (!all(c("levels", column) %in% names(profile))) return(NA_real_)
  at <- match(as.integer(level), as.integer(profile[["levels"]]))
  if (is.na(at)) return(NA_real_)
  as.numeric(profile[[column]])[at]
}


#' Whether a level's cells beat one mean, for the profile's estimate
#'
#' The rule of \code{.no_structure()} on the numbers a profile with an
#' estimate keeps: the leave-one-out error at the level, that of one mean
#' for every scored row, and the F test of the cells against one mean.  The
#' comparison is made only while the \code{loo} column is the one the
#' references were recorded with (\code{attr(, "estimate")$loo_levels}); a
#' profile whose column was edited has no reference.
#'
#' @param profile A \code{resolution_profile}.
#' @param level The level (cell count).
#' @param tol The margin by which the cells have to beat one mean.
#' @param n_searched The number of levels the level was chosen among.
#' @return A list: \code{structured}, \code{structure_p} and \code{beats}
#'   as \code{.no_structure()} returns them, \code{loo} (the error at the
#'   level) and \code{one_cell} (that of one mean).  Every element is
#'   \code{NA} for a profile without an estimate, without the column, or
#'   whose column is not the one recorded.
#' @keywords internal
#' @noRd
.estimate_structure <- function(profile, level, tol, n_searched = 1L) {
  none <- list(structured = NA, structure_p = NA_real_, beats = NA,
               loo = NA_real_, one_cell = NA_real_)
  est <- attr(profile, "estimate", exact = TRUE)
  if (!is.list(est) || !all(c("levels", "loo") %in% names(profile))) return(none)
  lv  <- profile[["levels"]]
  ref <- est$loo_levels
  one <- est$loo_one_cell %||% NA_real_
  same <- length(one) == 1L && is.finite(one) && !is.null(ref) &&
    all(as.character(lv) %in% names(ref)) &&
    isTRUE(all.equal(unname(ref[as.character(lv)]), as.numeric(profile[["loo"]])))
  if (!same) return(none)
  loo <- .profile_at(profile, level, "loo")
  p   <- if (is.null(est$f_p)) NA_real_ else unname(est$f_p[as.character(level)])[1L]
  r   <- .no_structure(loo, one, p, tol, n_searched = n_searched)
  list(structured = r$structured,
       structure_p = if (is.finite(r$structure_p)) r$structure_p else NA_real_,
       beats = r$beats, loo = loo, one_cell = one)
}


#' The record of how a count was set
#'
#' Built by \code{select_resolution()}, which puts it in the selection as
#' \code{count_rule} and, marked as sitting beside the scored centres, on
#' the selection's seeds.
#'
#' @param profile The \code{resolution_profile} the count was read off.
#' @param rule The criterion applied: \code{"size"}, \code{"loo"},
#'   \code{"cp"}, \code{"reliability"}, \code{"elbow"} or \code{"moran_z"}.
#' @param cells The count chosen.
#' @param obs_per_cell_asked The size asked for; \code{NA} unless
#'   \code{rule} is \code{"size"}.
#' @param structured Whether the cells at that count beat one mean for the
#'   profile's estimate, for the \code{scored} list.
#' @param scored_cells Whether the centres the record sits beside are the
#'   ones the profile scored: \code{NA} for a record that stands for a count
#'   alone.
#' @return A list of class \code{spatialkit_count_rule} with the fields
#'   \code{version}, \code{rule}, \code{estimate}, \code{cells},
#'   \code{obs_per_cell}, \code{obs_per_cell_asked}, \code{response},
#'   \code{n}, \code{n_points}, \code{n_fit}, \code{seed}, \code{features},
#'   \code{scored}, \code{partition}, \code{crs} and \code{scored_cells}.
#'   \code{estimate} is the estimate the rule was computed for: the
#'   profile's for \code{"size"} and \code{"loo"}, \code{"mean"} for
#'   \code{"cp"}, which scores plain means, and \code{NA} for the other
#'   rules.  \code{scored} is \code{NULL} for a profile without an estimate.
#' @keywords internal
#' @noRd
.count_rule_new <- function(profile, rule, cells, obs_per_cell_asked = NA_real_,
                            structured = NA, scored_cells = NA) {
  cells  <- as.integer(cells)
  bounds <- attr(profile, "bounds", exact = TRUE)
  est    <- attr(profile, "estimate", exact = TRUE)
  has_est <- is.list(est) && !is.null(est$value)
  int_of <- function(v)
    if (is.numeric(v) && length(v) == 1L && is.finite(v)) as.integer(v) else NA_integer_
  n <- .count_rule_rows(bounds)
  scored <- NULL
  if (has_est) {
    occupied <- if (is.null(est$cells)) NA_integer_
                else int_of(unname(est$cells[as.character(cells)])[1L])
    scored <- list(
      estimate     = as.character(est$value),
      n            = int_of(est$n),
      cells        = occupied,
      loo          = .profile_at(profile, cells, "loo"),
      loo_se       = .profile_at(profile, cells, "loo_se"),
      loo_one_cell = as.numeric(est$loo_one_cell %||% NA_real_)[1L],
      share_alone  = .profile_at(profile, cells, "share_alone"),
      own_weight   = .profile_at(profile, cells, "own_weight"),
      structured   = as.logical(structured)[1L],
      rows         = as.character(est$rows %||% NA_character_)[1L])
  }
  part   <- .partition_record(profile)
  nstart <- attr(profile, "nstart", exact = TRUE)
  crs    <- attr(profile, "centres_crs", exact = TRUE)
  structure(list(
    version  = 1L,
    rule     = rule,
    estimate = if (rule %in% c("size", "loo"))
                 (if (has_est) as.character(est$value) else NA_character_)
               else if (identical(rule, "cp")) "mean"
               else NA_character_,
    cells    = cells,
    obs_per_cell       = if (is.na(n)) NA_real_ else n / cells,
    obs_per_cell_asked = as.numeric(obs_per_cell_asked)[1L],
    response = if (has_est && !is.null(est$response)) as.character(est$response)[1L]
               else NA_character_,
    n        = n,
    n_points = if (is.list(bounds)) int_of(bounds$n) else NA_integer_,
    n_fit    = if (is.list(bounds)) int_of(bounds$n_sample) else NA_integer_,
    # A profile built with `seed = NULL` keeps none, and neither does one
    # without an estimate.
    seed     = if (has_est && !is.null(est$seed)) est$seed else NA_integer_,
    features = if (has_est && !is.null(est$features)) as.character(est$features)[1L]
               else NA_character_,
    scored   = scored,
    partition = c(part[.partition_fields()],
                  list(nstart = if (is.null(nstart)) NA_integer_ else nstart)),
    crs      = if (inherits(crs, "crs")) crs else sf::NA_crs_,
    scored_cells = as.logical(scored_cells)[1L]),
    class = "spatialkit_count_rule")
}
