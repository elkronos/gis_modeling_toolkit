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
#
# The record travels as the attribute "count_rule": from the selection's
# seeds to the seeds a seeding function returns, to the cells built on them
# (and `params$count_rule` of the tessellation), to the layer assigned to the
# cells, to the summary.  summarize_by_cell() reads it and checks two
# things: that it reports the estimate the count was chosen for, and that
# the cells and rows it is given are the ones the profile scored.


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


#' The centres a profile scored at one level
#'
#' The matrix \code{resolution_profile()} keeps in \code{attr(, "centres")}
#' for the level, which its estimate gave each scored row to by
#' \code{.nearest_centre()}, as it is: the coordinates in the profile's CRS
#' and the centres in the order they were numbered in.
#'
#' @param profile A \code{resolution_profile}.
#' @param level The level (cell count).
#' @return A two-column numeric matrix without dimnames, one row per centre;
#'   \code{NULL} when the profile holds no centres at that level (a column
#'   subset, which has lost them).
#' @keywords internal
#' @noRd
.profile_centres <- function(profile, level) {
  cen <- attr(profile, "centres", exact = TRUE)
  xy  <- if (is.list(cen)) cen[[as.character(level)]] else NULL
  if (!is.matrix(xy) || !is.numeric(xy) || nrow(xy) < 1L || ncol(xy) < 2L) return(NULL)
  xy <- unname(xy[, 1:2, drop = FALSE])
  storage.mode(xy) <- "double"
  xy
}


#' The centres a record sits beside
#'
#' @param rec A record of a count (\code{.count_rule_get()}), or \code{NULL}.
#' @return \code{rec$scored$centres} when the record says it sits beside the
#'   centres the profile scored (\code{scored_cells} is \code{TRUE}) and holds
#'   one finite centre per cell of the count, two at least, with a CRS object
#'   (which may be \code{NA}); \code{NULL} otherwise, and then no point is
#'   placed by the profile's rule.
#' @keywords internal
#' @noRd
.scored_centres <- function(rec) {
  if (is.null(rec) || !isTRUE(rec$scored_cells) || !is.list(rec$scored) ||
      !inherits(rec$crs, "crs"))
    return(NULL)
  cen <- rec$scored$centres
  if (!is.matrix(cen) || !is.numeric(cen) || ncol(cen) != 2L || nrow(cen) < 2L ||
      !identical(nrow(cen), as.integer(rec$cells)) || !all(is.finite(cen)))
    return(NULL)
  cen
}


#' The CRS the cells of a record's centres are drawn in
#'
#' @param rec A record of a count, or \code{NULL}.
#' @return The record's CRS when it sits beside the centres the profile
#'   scored (\code{.scored_centres()}) and that CRS is a projected one;
#'   \code{NULL} otherwise.
#' @keywords internal
#' @noRd
.scored_crs <- function(rec) {
  if (is.null(.scored_centres(rec))) return(NULL)
  crs <- rec$crs
  if (is.na(crs) || isTRUE(sf::st_is_longlat(crs))) return(NULL)
  crs
}


#' The nearest and the next nearest centre of every point
#'
#' @param x,y Coordinates of the points.
#' @param cx,cy Coordinates of the centres, two at least.
#' @param use_fnn Whether to use \pkg{FNN}.  Without it the search is
#'   exact, over every centre.
#' @return A list: \code{first}, the nearest centre of each point, and
#'   \code{d1} and \code{d2}, the squared distances to it and to the next
#'   nearest, both computed in R from the coordinates.
#' @keywords internal
#' @noRd
.two_nearest <- function(x, y, cx, cy, use_fnn = .has_backend("FNN")) {
  if (isTRUE(use_fnn) && length(cx) > 16L) {
    nn <- matrix(FNN::get.knnx(cbind(cx, cy), cbind(x, y), k = 2L)$nn.index, length(x), 2L)
    a <- nn[, 1L]; b <- nn[, 2L]
    da <- (x - cx[a])^2 + (y - cy[a])^2
    db <- (x - cx[b])^2 + (y - cy[b])^2
    swap <- db < da
    return(list(first = ifelse(swap, b, a), d1 = pmin(da, db), d2 = pmax(da, db)))
  }
  d1 <- (x - cx[1L])^2 + (y - cy[1L])^2
  d2 <- rep(Inf, length(x))
  first <- rep.int(1L, length(x))
  for (j in seq_along(cx)[-1L]) {
    d <- (x - cx[j])^2 + (y - cy[j])^2
    nearer <- d < d1
    d2 <- ifelse(nearer, d1, pmin(d2, d))
    d1[nearer] <- d[nearer]
    first[nearer] <- j
  }
  list(first = first, d1 = d1, d2 = d2)
}


#' Which centre each seed is
#'
#' Seeds are the centres a record holds when they are as many, and each lies
#' at a different centre: within one part in a thousand of its distance from
#' the next nearest centre, or within \code{.edge_tol()} of the centres'
#' coordinates when that is more.  The seeds a selection carries, and those a
#' seeding function returns for one, are the centres in the profile's order,
#' in the CRS of the points they were placed for; put back in the record's
#' CRS they have moved by the rounding of the two transformations, a few
#' millimetres at most (Lambert azimuthal equal-area across a continent), and
#' a thousandth of the distance between two centres is far more than that
#' for centres more than a few metres apart.  Seeds in another order, sorted
#' by their cell for example, are matched all the same.  Seeds placed by a
#' fit of their own lie metres from the centres and are not.
#'
#' @param xy Two-column matrix of the seeds' coordinates, in the record's
#'   CRS.
#' @param cen The record's centres (\code{.scored_centres()}).
#' @param use_fnn Whether to use \pkg{FNN}.
#' @return Integer vector, for each centre the row of \code{xy} that is it;
#'   \code{NULL} when the seeds are not the centres.
#' @keywords internal
#' @noRd
.seeds_as_centres <- function(xy, cen, use_fnn = .has_backend("FNN")) {
  if (is.null(cen) || !is.matrix(xy) || ncol(xy) < 2L) return(NULL)
  K <- nrow(cen)
  if (nrow(xy) != K || K < 2L || !all(is.finite(xy[, 1:2]))) return(NULL)
  two <- .two_nearest(xy[, 1L], xy[, 2L], cen[, 1L], cen[, 2L], use_fnn)
  if (anyDuplicated(two$first)) return(NULL)
  tol <- .edge_tol(c(range(cen[, 1L]), range(cen[, 2L])))
  if (!all(sqrt(two$d1) <= pmax(1e-3 * sqrt(two$d2), tol))) return(NULL)
  of <- integer(K)
  of[two$first] <- seq_len(K)
  of
}


#' How far from an edge a point can lie and still be on it
#'
#' The corners of Voronoi cells are computed in floating point, so a point
#' the same distance from two centres can lie a little inside either cell.
#' Rounding grows with the size of the coordinates, so the margin is one
#' part in 1e10 of the largest coordinate of the layer: 0.44 mm at a UTM
#' northing of 4,400 km, far more than the rounding of the corners and far
#' less than the spacing of gridded locations.  The cells
#' \code{create_voronoi_polygons()} draws on seeds that are the centres are
#' drawn around the centres as the record holds them
#' (\code{.seeds_as_centres()}), so rounding the seeds took on through
#' another CRS does not move their edges.
#'
#' @param bb A bounding box (\code{sf::st_bbox()}), in the CRS the distance
#'   is taken in.
#' @return One positive number, in the units of that CRS.
#' @keywords internal
#' @noRd
.edge_tol <- function(bb) {
  m <- suppressWarnings(max(abs(as.numeric(bb)), na.rm = TRUE))
  1e-10 * max(1, if (is.finite(m)) m else 1)
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
#'   Its last element, \code{centres}, is the matrix of the centres the
#'   profile scored at the count, in \code{crs} and in the profile's order
#'   (\code{NULL} when the profile no longer holds them): each scored row was
#'   given to the nearest of them, a tie to the lower-numbered one, and
#'   \code{assign_features_to_polygons()} and \code{cell_values_at()} read
#'   them to place a point on the edge between two cells the same way.
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
      rows         = as.character(est$rows %||% NA_character_)[1L],
      centres      = .profile_centres(profile, cells))
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


#' The record of a count a layer carries
#'
#' A record is read only when it is one this code can read: a list of class
#' \code{spatialkit_count_rule} whose \code{version} is \code{1L} and
#' which names one rule, one count and one estimate (which may be
#' \code{NA}).  Anything else is treated as no record, so a record of
#' another version, or one a caller cut down, is not checked against and not
#' handed on.
#'
#' @param x Any object.
#' @return \code{attr(x, "count_rule")} when it is such a record;
#'   \code{NULL} otherwise.
#' @keywords internal
#' @noRd
.count_rule_get <- function(x) {
  rec <- attr(x, "count_rule", exact = TRUE)
  if (!is.list(rec) || !inherits(rec, "spatialkit_count_rule")) return(NULL)
  one <- function(v) length(v) == 1L
  readable <- identical(rec$version, 1L) &&
    is.character(rec$rule) && one(rec$rule) && !is.na(rec$rule) &&
    is.numeric(rec$cells) && one(rec$cells) && !is.na(rec$cells) &&
    one(rec$estimate) && (is.character(rec$estimate) || is.na(rec$estimate))
  if (readable) rec else NULL
}


#' Put the record of a count on a layer
#'
#' @param x The layer: seeds, cells or assigned rows.
#' @param rec The record, or \code{NULL} for none.
#' @param scored_cells \code{NULL} to leave the record's \code{scored_cells}
#'   as it is; otherwise its new value: \code{TRUE} beside the centres the
#'   profile scored, \code{FALSE} beside other centres placed at that count.
#' @return \code{x}, with \code{attr(, "count_rule")} set when \code{rec} is
#'   a record and left as it was otherwise.
#' @keywords internal
#' @noRd
.count_rule_put <- function(x, rec, scored_cells = NULL) {
  if (is.null(rec)) return(x)
  if (!is.null(scored_cells)) rec$scored_cells <- as.logical(scored_cells)[1L]
  attr(x, "count_rule") <- rec
  x
}


# The fields in which two records of one path agree, in the order they are
# compared: the first that differs is the one a warning names.  The first
# five say which selection the cells were built from; `scored_cells` says
# whether they are that selection's scored cells or other cells built from
# its count (a lattice, seeds a function placed itself), so rows assigned to
# the one and summarised with the other differ in it and in nothing else.
.count_rule_same_fields <- c("rule", "estimate", "cells", "crs", "scored$loo",
                             "scored_cells")

#' The first field in which two records differ
#'
#' A layer of assigned rows and the cells it is summarised with each carry
#' the record of the selection their cells were built from.  When the two
#' selections are not the same one, the records differ in the rule, the
#' estimate, the count, the CRS or the leave-one-out error the profile
#' scored.  When they are the same selection and the cells are not the same
#' cells (its scored cells beside a lattice sized by it), the records differ
#' in \code{scored_cells}.
#'
#' @param a,b Two records.
#' @return \code{NULL} when the two agree in all six; otherwise a list with
#'   \code{field}, the name of the first that differs, and \code{a} and
#'   \code{b}, its two values as text.
#' @keywords internal
#' @noRd
.count_rule_differs <- function(a, b) {
  value <- function(rec, field) switch(field,
    "crs"        = rec$crs,
    "scored$loo" = if (is.list(rec$scored)) rec$scored$loo else NULL,
    rec[[field]])
  text <- function(v, field) {
    if (is.null(v)) return("none")
    if (identical(field, "crs"))
      return(if (inherits(v, "crs") && !is.na(v)) as.character(v$input)[1L] else "NA")
    if (is.numeric(v) && !is.integer(v)) return(format(v, digits = 8))
    v <- as.character(v)[1L]
    if (is.na(v)) "NA" else v
  }
  same <- function(x, y, field) {
    if (is.null(x) || is.null(y)) return(is.null(x) && is.null(y))
    if (identical(field, "crs")) {
      if (!inherits(x, "crs") || !inherits(y, "crs")) return(identical(x, y))
      if (is.na(x) || is.na(y)) return(is.na(x) && is.na(y))
      return(isTRUE(x == y))
    }
    if (identical(field, "scored$loo")) return(.same_loo(x, y))
    identical(x, y)
  }
  for (field in .count_rule_same_fields) {
    x <- value(a, field); y <- value(b, field)
    if (!same(x, y, field))
      return(list(field = field, a = text(x, field), b = text(y, field)))
  }
  NULL
}


#' Whether two leave-one-out errors are one number
#'
#' @param a,b Two errors.
#' @return \code{TRUE} when both are missing, or when they agree to a
#'   relative tolerance of \code{.RECORD_TOL}; \code{FALSE} for a missing
#'   error against a number.
#' @keywords internal
#' @noRd
.same_loo <- function(a, b) {
  a <- as.numeric(a)[1L]; b <- as.numeric(b)[1L]
  if (is.na(a) || is.na(b)) return(is.na(a) && is.na(b))
  abs(a - b) <= .RECORD_TOL * max(abs(a), abs(b))
}


#' Check a summary against the record of its cells' count
#'
#' Three things, in this order, each with a warning of its own where it
#' fails.  The warnings are raised here; the summary is returned as asked.
#' \enumerate{
#'   \item Two records.  When the cells handed to the summary carry another
#'     record than the assigned rows do, the rows were assigned to other
#'     cells than those: class \code{spatialkit_cells_differ} with
#'     \code{what = "record"}.  The checks below run against the record of
#'     the assigned rows.
#'   \item The estimate.  When the count was chosen for plain means or for
#'     pooled values and the summary reports the other: class
#'     \code{spatialkit_estimate_mismatch}.
#'   \item The cells and the rows.  The engine is run on the response with
#'     the estimate the profile scored, and its scored rows, its occupied
#'     cells and its leave-one-out error are compared with the record's, the
#'     last to a relative tolerance of \code{.RECORD_TOL}.  The first that
#'     differs is named: class \code{spatialkit_cells_differ} with
#'     \code{what} one of \code{"n"}, \code{"cells"} and \code{"loo"}.
#' }
#' The third is skipped, and the reason kept, where the record's numbers do
#' not describe the cells in hand or cannot be computed again: the profile
#' scored no estimate, or read the selection half; other cells were built
#' from the count, or the record does not say whose centres these are; the
#' response is another column, or was not summarised; the input has no
#' coordinates; the features summarised are not points; the fit was stopped
#' before an exact search above the limit (\code{.work_guard()}).
#'
#' Features that are not points are assigned by largest overlap, where the
#' profile scored each at one point in the cell of the nearest centre, so a
#' difference there says nothing about rows dropped.  It is the layer
#' summarised that decides: the points of such features
#' (\code{coerce_to_points()}), assigned as points, are the rows the profile
#' scored, and are checked like any points, whatever the record's
#' \code{features} says of the layer profiled.
#'
#' @param rec The record the summary reads: that of the assigned rows, or of
#'   the cells when the rows carry none.
#' @param other The record of the cells when the assigned rows carry one
#'   too; \code{NULL} otherwise.
#' @param estimate The estimate the summary reports.
#' @param response_var The call's \code{response_var}, or \code{NULL}.
#' @param summarised The response columns the call summarised.
#' @param layer The layer summarised, as the caller passed it.
#' @param found A function of one estimate, \code{"mean"} or
#'   \code{"pooled"}, that returns what the engine finds for the response on
#'   the rows summarised: a list with \code{n} (scored rows), \code{cells}
#'   (the cells that hold one), \code{loo} (the leave-one-out error of that
#'   estimate) and \code{status} (\code{"work limit"} when the fit was
#'   stopped).  It is called only when the third check runs.
#' @return \code{rec} with one more element, \code{check}: a list with
#'   \code{estimate} (\code{"match"}, \code{"mismatch"}, or \code{NA} for a
#'   count chosen for neither estimate) and \code{cells} (\code{"match"},
#'   \code{"differ"}, or \code{"not checked: <reason>"}).  \code{cells} is
#'   \code{"differ"} also when the two records differ.
#' @keywords internal
#' @noRd
.count_rule_check <- function(rec, other = NULL, estimate, response_var, summarised,
                              layer, found) {
  value_name <- function(est) if (identical(est, "pooled")) "pooled value" else "plain mean"
  check <- list(estimate = NA_character_, cells = NA_character_)

  # ---- 1. Two records ---------------------------------------------------------
  records_differ <- FALSE
  if (!is.null(other)) {
    d <- .count_rule_differs(rec, other)
    if (!is.null(d)) {
      records_differ <- TRUE
      .warn_classed("spatialkit_cells_differ",
                    paste0("summarize_by_cell(): `assigned_points_sf` and `cells_sf` carry ",
                           "different count records (%s: %s against %s), so the rows were ",
                           "assigned to other cells than the `cells_sf` given. The record of ",
                           "the assigned layer is used."),
                    d$field, d$a, d$b,
                    data = list(what = "record", field = d$field, scored = d$a, found = d$b))
    }
  }

  # ---- 2. The estimate --------------------------------------------------------
  chosen <- rec$estimate
  if (is.character(chosen) && length(chosen) == 1L && chosen %in% c("mean", "pooled")) {
    check$estimate <- if (identical(chosen, estimate)) "match" else "mismatch"
    if (identical(check$estimate, "mismatch"))
      .warn_classed("spatialkit_estimate_mismatch",
                    paste0("summarize_by_cell(): the number of cells (%d, set by %s) was ",
                           "chosen for the %s, and this summary reports the %s. %s"),
                    as.integer(rec$cells), rec$rule, value_name(chosen), value_name(estimate),
                    if (identical(chosen, "pooled"))
                      paste0("At this cell size plain means of small cells are noisy: pass ",
                             "estimate = \"pooled\", or choose the count for plain means ",
                             "(resolution_profile(estimate = \"mean\")).")
                    else
                      paste0("Pooled values are valid on these cells; the count that suits ",
                             "them is set by a size (criterion = \"size\")."),
                    data = list(chosen_for = chosen, reported = estimate, rule = rec$rule,
                                cells = as.integer(rec$cells)))
  }

  # ---- 3. The cells and the rows ----------------------------------------------
  sc <- rec$scored
  why <- if (!is.list(sc) || !isTRUE(sc$estimate %in% c("mean", "pooled")))
           "the profile scored no estimate"
         else if (!identical(sc$rows, "all")) "the profile read the selection half"
         else if (isFALSE(rec$scored_cells)) "other cells were built from this count"
         else if (!isTRUE(rec$scored_cells)) "the record does not say whose centres these are"
         else if (is.null(response_var)) "the response was not summarised"
         else if (!identical(response_var, rec$response)) "the response is another column"
         else if (!response_var %in% summarised) "the response was not summarised"
         else if (!inherits(layer, "sf")) "the input has no coordinates"
         else if (!all(sf::st_geometry_type(layer, by_geometry = TRUE) == "POINT"))
           "the features are not points"
         else NULL
  if (is.null(why)) {
    got <- found(sc$estimate)
    if (identical(got$status, "work limit")) why <- "work limit"
  }
  if (!is.null(why)) {
    check$cells <- if (records_differ) "differ" else paste("not checked:", why)
    rec$check <- check
    return(rec)
  }

  int_of <- function(v) if (is.numeric(v) && length(v) == 1L) as.integer(v) else NA_integer_
  num_of <- function(v) if (is.numeric(v) && length(v) == 1L) as.numeric(v) else NA_real_
  a <- list(n = int_of(sc$n), cells = int_of(sc$cells), loo = num_of(sc$loo))
  b <- list(n = int_of(got$n), cells = int_of(got$cells), loo = num_of(got$loo))
  what <- if (!identical(a$n, b$n)) "n"
          else if (!identical(a$cells, b$cells)) "cells"
          else if (!.same_loo(a$loo, b$loo)) "loo"
          else NULL
  check$cells <- if (is.null(what) && !records_differ) "match" else "differ"
  if (!is.null(what)) {
    num <- function(v) if (is.na(v)) "NA" else format(v, digits = 6)
    .warn_classed("spatialkit_cells_differ",
                  paste0("summarize_by_cell(): these are not the cells and rows the count ",
                         "was scored on (%s differs). The profile scored %d rows of '%s' in ",
                         "%d cells with a leave-one-out error of %s; here %d rows in %d ",
                         "cells give %s. %s"),
                  c(n = "the number of rows", cells = "the number of cells",
                    loo = "the leave-one-out error")[[what]],
                  a$n, rec$response, a$cells, num(a$loo), b$n, b$cells, num(b$loo),
                  if (identical(what, "loo"))
                    paste0("The same number of rows and cells was found, so rows sit in ",
                           "other cells than the profile put them in, or values changed: ",
                           "rows exactly on a boundary between two cells (gridded ",
                           "locations) can be assigned to the other cell, and coordinates ",
                           "that reached this CRS by another transformation can settle a ",
                           "tie between neighbouring cells the other way.")
                  else
                    paste0("Rows were dropped or moved between cells: points outside every ",
                           "cell, a boundary that clips the cells, a subset of the rows, or ",
                           "another layer."),
                  data = list(what = what, scored = a, found = b))
  }
  rec$check <- check
  rec
}


#' @export
print.spatialkit_count_rule <- function(x, ...) {
  label  <- function(name) sprintf("  %-12s: ", name)
  indent <- strrep(" ", 16L)
  one    <- function(v) length(v) == 1L && !is.na(v)
  # An error is written as print.spatialkit_pooled() writes it, so the two
  # records of one summary show one number one way.
  num    <- function(v) if (is.numeric(v) && one(v)) sprintf("%.2e", v) else "NA"
  count  <- function(k, what)
    sprintf("%s %s%s", format(k), what, if (isTRUE(k == 1L)) "" else "s")

  # ---- the count and what set it ----------------------------------------------
  asked <- x$obs_per_cell_asked
  how <- c(if (identical(x$rule, "size") && is.numeric(asked) && one(asked))
             paste(.per_cell_text(asked), "asked"),
           if (is.numeric(x$n) && one(x$n))
             sprintf("%s / %s = %.2f", count(x$n, "row"), count(x$cells, "cell"),
                     x$n / x$cells))
  cat(sprintf("Count: %s, set by %s%s\n", count(x$cells, "cell"), x$rule,
              if (length(how)) sprintf(" (%s)", paste(how, collapse = "; ")) else ""))

  # ---- the estimate the rule was computed for ---------------------------------
  cat(label("chosen for"),
      if (identical(x$estimate, "pooled")) "the pooled value"
      else if (identical(x$estimate, "mean")) "the plain mean"
      # A size and the leave-one-out error score the profile's estimate, so
      # without one the object the count was read off named none (a profile
      # built without an estimate, or a column subset, which loses it).
      else if (isTRUE(x$rule %in% c("size", "loo")))
        "no cell value: the profile names no estimate"
      else sprintf("no cell value: %s scores none", x$rule),
      "\n", sep = "")

  # ---- the cells --------------------------------------------------------------
  part   <- x$partition
  nstart <- part$nstart
  starts <- .partition_label(part, if (is.numeric(nstart) && one(nstart)) nstart else NULL)
  # The seed with every digit: 100000, not 1e+05.
  seed <- if (is.numeric(x$seed) && one(x$seed)) format(x$seed, scientific = FALSE)
  if (!is.null(seed))
    starts <- if (endsWith(starts, ")"))
      paste0(substr(starts, 1L, nchar(starts) - 1L), "; seed ", seed, ")")
    else sprintf("%s (seed %s)", starts, seed)
  fitted <- if (!one(x$n_fit) || !one(x$n_points)) NULL
            else if (x$n_fit < x$n_points)
              sprintf("fitted to a subsample of %d of the %d points",
                      as.integer(x$n_fit), as.integer(x$n_points))
            else "fitted to every row"
  half  <- is.list(x$scored) && identical(x$scored$rows, "selection half")
  whose <- if (isTRUE(x$scored_cells))
             paste0("these are the cells the profile scored", if (half) ", on half the rows")
           else if (isFALSE(x$scored_cells))
             "other cells were built from this count: not the cells the profile scored"
           else "the record stands for the count and sits beside no centres"
  # What built the cells, wrapped where the density kind's settings make it
  # long; whose cells they are starts a line of its own.
  built <- strwrap(paste0(paste(c(if (identical(part$kind, "density")) "density"
                                  else if (identical(part$kind, "kmeans")) "k-means",
                                  starts, fitted), collapse = ", "), ";"),
                   width = 64)
  cat(label("cells"), paste(built, collapse = paste0("\n", indent)), "\n",
      indent, whose, "\n", sep = "")

  # ---- what the profile scored at the count -----------------------------------
  sc <- x$scored
  if (is.list(sc))
    cat(label("scored"),
        sprintf("loo %s on %s in %s (one cell: %s)%s\n", num(sc$loo),
                if (half) sprintf("the %s of the selection half", count(sc$n, "row"))
                else count(sc$n, "row"),
                count(sc$cells, "cell"), num(sc$loo_one_cell),
                if (is.numeric(sc$share_alone) && one(sc$share_alone))
                  sprintf("; %.1f%% of rows alone", 100 * sc$share_alone)
                else ""),
        sep = "")

  # ---- the summary's checks, on a summary only --------------------------------
  ck <- x$check
  if (is.list(ck)) {
    skipped <- is.character(ck$cells) && startsWith(ck$cells, "not checked: ")
    cat(label(if (skipped) "not checked" else "checked"),
        if (identical(ck$cells, "match")) "the cells and rows summarised are the ones scored"
        else if (identical(ck$cells, "differ"))
          "the cells and rows summarised are not the ones scored"
        else if (skipped) substring(ck$cells, nchar("not checked: ") + 1L)
        else ck$cells,
        "\n", sep = "")
    if (identical(ck$estimate, "mismatch"))
      cat(indent,
          if (identical(x$estimate, "pooled"))
            "the summary reports the plain mean, not the pooled value"
          else "the summary reports the pooled value, not the plain mean",
          "\n", sep = "")
  }
  invisible(x)
}
