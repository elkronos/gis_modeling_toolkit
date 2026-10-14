# R/cell-values-at.R
# -----------------------------------------------------------------------------
# The value of a cell at other locations
# -----------------------------------------------------------------------------
# cell_values_at() gives new locations the value and interval of their cell,
# by the rule the held-out rows of an evaluation are scored with: the cell
# of the nearest seed (or the polygon that holds the location), and for a
# location in no cell, or in a cell with neither a scored row nor a value,
# the cell of the nearest scored row.  Beside each value it returns the
# location's distance from the nearest observation.


#' The columns of a summary that hold a column's cell value
#'
#' The response's columns are read before a predictor's, and within either
#' the pooled value before the plain mean.
#'
#' @param nms The names of the summary's columns.
#' @param v The summarised column.
#' @return \code{NULL} when the summary holds neither value of \code{v};
#'   otherwise a list: \code{value} (the column's name), \code{estimate}
#'   (\code{"pooled"} or \code{"mean"}) and \code{lo}, \code{hi} (the names
#'   of the interval's columns, or \code{NULL} when the summary does not hold
#'   both beside a pooled value).
#' @keywords internal
#' @noRd
.values_at_columns <- function(nms, v) {
  for (p in c("resp_", "pred_")) {
    pooled <- paste0("..pooled_", p, v)
    plain  <- paste0(p, "mean_", v)
    if (pooled %in% nms) {
      lo <- paste0("..obs_lo_", p, v)
      hi <- paste0("..obs_hi_", p, v)
      both <- all(c(lo, hi) %in% nms)
      return(list(value = pooled, estimate = "pooled",
                  lo = if (both) lo, hi = if (both) hi))
    }
    if (plain %in% nms)
      return(list(value = plain, estimate = "mean", lo = NULL, hi = NULL))
  }
  NULL
}


#' The geometry of a layer, one point per feature
#'
#' @param layer An sf layer.
#' @return An sf layer with the geometry column and no other, row for row
#'   with \code{layer}.  A feature that is not a point is reduced to one by
#'   \code{coerce_to_points(, "auto")}, whose conditions are passed on; an
#'   empty geometry stays an empty point.
#' @keywords internal
#' @noRd
.values_at_points <- function(layer) {
  g <- sf::st_sf(geometry = sf::st_geometry(layer))
  if (nrow(g) > 0L && !all(sf::st_geometry_type(g, by_geometry = TRUE) == "POINT"))
    g <- coerce_to_points(g, "auto")
  g
}


#' A point layer where the scored rows are read
#'
#' @param pts A layer of points, as \code{.values_at_points()} returns it.
#' @param arg The argument the layer came from, for the messages.
#' @param what What the layer holds, in the words of the messages.
#' @param crs The working CRS; \code{NA} when the coordinates of
#'   \code{points} are used as they are.
#' @param points_crs The CRS of \code{points} as it was passed.
#' @param as_lonlat Whether \code{points} carries no CRS and its coordinates
#'   were taken for longitude and latitude.
#' @return \code{pts} in the CRS the scored rows are measured in.  A layer
#'   without a CRS beside \code{points} without one is read as
#'   \code{points} was read: as EPSG:4326 when \code{as_lonlat}, and as it is
#'   otherwise.
#' @keywords internal
#' @noRd
.values_at_place <- function(pts, arg, what, crs, points_crs, as_lonlat) {
  own <- sf::st_crs(pts)
  if (is.na(own)) {
    if (!is.na(points_crs))
      stop(sprintf(paste0("cell_values_at(): `%s` has no CRS and `points` has one; set it ",
                          "with sf::st_set_crs() so the %s can be placed."), arg, what),
           call. = FALSE)
    if (as_lonlat) pts <- sf::st_transform(sf::st_set_crs(pts, 4326), crs)
    return(pts)
  }
  if (is.na(crs))
    stop(sprintf(paste0("cell_values_at(): `%s` has a CRS and `points` has none; set the ",
                        "CRS of `points` with sf::st_set_crs() so the %s can be placed."),
                 arg, what), call. = FALSE)
  ensure_projected(pts, target_crs = crs)
}


#' The coordinates of a point layer
#'
#' @param pts A layer of points.
#' @return A two-column matrix, one row per point; \code{NA} for an empty
#'   point and for one with a coordinate beyond \code{.POOL_XY_MAX}, from
#'   which no distance can be taken.
#' @keywords internal
#' @noRd
.values_at_xy <- function(pts) {
  if (nrow(pts) == 0L) return(matrix(numeric(0), 0L, 2L))
  xy  <- unname(sf::st_coordinates(pts)[, 1:2, drop = FALSE])
  far <- abs(xy[, 1L]) > .POOL_XY_MAX | abs(xy[, 2L]) > .POOL_XY_MAX
  xy[which(far | !is.finite(xy[, 1L]) | !is.finite(xy[, 2L])), ] <- NA_real_
  xy
}


#' The nearest centre of every point, where many points can tie
#'
#' The answer of \code{.nearest_centre()}: the nearest centre by squared
#' distance in R's arithmetic, a tie to the lower centre number.  That
#' function searches a point whose two nearest centres tie against every
#' centre, which is the right trade where the centres are a partition's and
#' few points tie.  Here the centres can be every scored location of a layer
#' and the points a grid laid over them, so that most points tie; searched
#' that way the time grows with the points times the centres.  A tied point
#' is given a longer shortlist from the kd-tree instead, twice as long each
#' round, until its furthest member is further than its nearest by more than
#' one part in a thousand million: every centre as near as the nearest is
#' then on the list, and the lowest-numbered of them is taken.  A point whose
#' shortlist would hold every centre is searched exactly.
#'
#' @param x,y Coordinates of the points.
#' @param cx,cy Coordinates of the centres.
#' @param use_fnn Whether to use \pkg{FNN}.  Without it, and for up to 16
#'   centres, the search is \code{.nearest_centre()}'s exact one.
#' @return Integer vector: the centre of each point.
#' @keywords internal
#' @noRd
.values_at_nearest <- function(x, y, cx, cy, use_fnn) {
  K <- length(cx)
  if (!isTRUE(use_fnn) || K <= 16L) return(.nearest_centre(x, y, cx, cy, use_fnn = FALSE))
  cen  <- cbind(cx, cy)
  out  <- integer(length(x))
  todo <- seq_along(x)
  k <- 2L
  while (length(todo)) {
    if (k >= K) {
      out[todo] <- .nearest_exact(x[todo], y[todo], cx, cy)
      break
    }
    qx <- x[todo]; qy <- y[todo]
    nn <- matrix(FNN::get.knnx(cen, cbind(qx, qy), k = k)$nn.index, length(todo), k)
    d2 <- matrix((qx - cx[nn])^2 + (qy - cy[nn])^2, length(todo), k)
    best <- do.call(pmin, lapply(seq_len(k), function(j) d2[, j]))
    far  <- do.call(pmax, lapply(seq_len(k), function(j) d2[, j]))
    # The lowest-numbered of the centres at the least distance.
    nn[d2 > best] <- K + 1L
    pick <- do.call(pmin, lapply(seq_len(k), function(j) nn[, j]))
    done <- far > best * (1 + 1e-9)
    out[todo[done]] <- pick[done]
    todo <- todo[!done]
    k <- 2L * k
  }
  out
}


#' The centres the profile scored, when the seeds are those centres
#'
#' The seeds of a \code{select_resolution()} result, and the seeds
#' \code{get_voronoi_seeds()} hands on from one, are the centres the profile
#' scored, and come back in the CRS of the points they were placed for.  Put
#' back in the CRS the profile measured in, they differ from the centres by
#' the rounding of the two transformations, which can settle a tie between
#' two seeds the other way, and a caller may hand them in another order
#' (sorted by their cell, say), which settles a tie by that order.  The
#' centres as the record holds them, in its order, settle it as the profile
#' did.
#'
#' @param rule The summary's record of the count, or \code{NULL}.
#' @param seed_xy The seeds' coordinates in the working CRS, one row per row
#'   of \code{seeds}.
#' @param crs The working CRS.
#' @param use_fnn Whether \pkg{FNN} is used, as the call's other searches
#'   are told.
#' @return When the record sits beside the centres
#'   (\code{.scored_centres()}), the working CRS is the record's (or both
#'   have none), and the seeds are the centres (\code{.seeds_as_centres()}):
#'   a list of \code{centres}, the record's, and \code{seed}, for each centre
#'   the row of \code{seeds} that is it.  \code{NULL} otherwise, and the
#'   seeds are searched as they are.
#' @keywords internal
#' @noRd
.values_at_centres <- function(rule, seed_xy, crs, use_fnn = .has_backend("FNN")) {
  cen <- .scored_centres(rule)
  if (is.null(cen)) return(NULL)
  if (!inherits(crs, "crs") || is.na(crs) != is.na(rule$crs)) return(NULL)
  if (!is.na(crs) && !isTRUE(crs == rule$crs)) return(NULL)
  of <- .seeds_as_centres(seed_xy, cen, use_fnn)
  if (is.null(of)) return(NULL)
  list(centres = cen, seed = of)
}


#' Cell IDs written in the class of a summary's ID column
#'
#' IDs are matched as text.  An ID that matches a row of the summary is that
#' row's.  One that matches none (a cell the summary has no row for) is
#' converted from its text: to a number for a numeric column, to a level
#' added to the levels for a factor.
#'
#' @param text Character vector of IDs, \code{NA} for none.
#' @param ids The summary's ID column.
#' @return A vector of the class of \code{ids}, as long as \code{text}, when
#'   every ID can be written in that class.  When one that matches no row
#'   cannot (a string that is no number or a fraction beside an integer
#'   column, a string that is not the text of a number beside a double one,
#'   any ID beside a column of another class than character, factor, integer
#'   and double), the result is \code{text}: an ID is never turned into
#'   another ID or into \code{NA}.
#' @keywords internal
#' @noRd
.values_at_ids <- function(text, ids) {
  at   <- match(text, .id_as_character(ids), incomparables = NA)
  miss <- which(!is.na(text) & is.na(at))
  if (is.character(ids)) return(text)
  if (is.factor(ids)) {
    lev <- c(levels(ids), setdiff(unique(text[miss]), levels(ids)))
    return(factor(text, levels = lev, ordered = is.ordered(ids)))
  }
  out <- ids[at]
  if (!length(miss)) return(out)
  read <- if (is.integer(ids) && !is.object(ids)) suppressWarnings(as.integer(text[miss]))
          else if (is.double(ids) && !is.object(ids)) suppressWarnings(as.numeric(text[miss]))
          else NULL
  # Written back, the number has to be the ID it was read from.
  if (is.null(read) || !identical(.id_as_character(read), text[miss])) return(text)
  out[miss] <- read
  out
}


#' Cell values at other locations
#'
#' Gives each location of `newdata` the value of the cell it lies in, read
#' from a [summarize_by_cell()] result: the pooled value when the summary
#' holds one, the plain mean otherwise, with the interval for one
#' observation when the summary holds it. Beside each value it returns the
#' location's distance from the nearest observation, since a cell's value
#' describes the part of the cell where its observations lie.
#'
#' The function works in six steps.
#'
#' 1. *The working CRS.* The CRS the summary measured in: the one in
#'    `attr(values, "count_rule")` when that record names a projected CRS,
#'    else the one in `attr(values, "pooled")`, else the CRS of `points`
#'    when it is projected and the one [ensure_projected()] picks for it
#'    when it is not. `points`, `newdata` and `seeds` are transformed to it.
#'    When neither `points` nor `newdata` has a CRS the coordinates are used
#'    as they are, unless those of `points` look like longitude and latitude:
#'    both are then read as EPSG:4326, as [summarize_by_cell()] reads
#'    `points`, with [ensure_projected()]'s warning.
#' 2. *The scored rows.* The rows of `points` with a finite `response_var`,
#'    a cell ID and finite coordinates. Their distinct locations are the
#'    observations every distance below is taken from.
#' 3. *The cell of a location.* With `seeds`, the cell of the nearest seed.
#'    That is the Voronoi cell before its clip to the cover, so every
#'    location has one, unless its nearest seed has an `NA` ID
#'    ([build_tessellation()] gives `NA` in `index` to a seed whose cell the
#'    clip removed). Without `seeds`, the polygon of `values` that holds the
#'    location, by [assign_features_to_polygons()] with its default
#'    `tie_break` and the summary's record of the count, so a location on
#'    the edge between two cells the profile scored is placed as that
#'    function places such a row (see "Ties"); a location in no polygon has
#'    no cell. Polygons without a
#'    CRS are read as that function read them for the rows of `points`: the
#'    locations are handed to it in the CRS of `points`, whatever CRS
#'    `newdata` arrived in.
#' 4. *The value.* The cell's pooled value of `response_var`,
#'    `..pooled_<prefix><response_var>` with the prefix `resp_` or `pred_`,
#'    when `values` has that column; otherwise its mean,
#'    `<prefix>mean_<response_var>`. The column `estimate` says which. For a
#'    column summarised as response and as predictor the response's columns
#'    are read. The interval is the cell's `..obs_lo_<prefix><response_var>`
#'    and `..obs_hi_<prefix><response_var>` when `values` has them.
#' 5. *The fill.* A location with no cell, or whose cell holds no scored row
#'    and has no value in `values`, takes the value and interval of the cell
#'    of its nearest scored row. `filled` is then `TRUE` and `filled_from` is
#'    that cell's ID. The ID column still holds the location's own cell, so
#'    the two columns together say where a location lies and where its value
#'    came from. A cell that holds a scored row is not filled. Where `values`
#'    has no value for such a cell (the summary could not compute the pooled
#'    columns, or has no row for the cell), its locations have `value = NA`
#'    and `filled = FALSE`; a filled value is `NA` in the same case, when
#'    the cell it is taken from has none.
#' 6. *The distance.* `dist_nearest` is the distance from the location to
#'    the nearest location of a scored row, in the units of the working CRS.
#'    `dist_ratio` is that distance over the layer's spacing: the median,
#'    over the distinct locations of the scored rows, of the distance to the
#'    nearest other one.
#'
#' Cell IDs are matched as text. The ID column of `values` is character when
#' the summary joined a `cells_sf` whose IDs were of another class than the
#' assigned layer's, while `points` and `seeds` keep the class they had; the
#' IDs of all three are written out as text before any match, and both ID
#' columns returned have the class of the ID column of `values` wherever
#' that class can hold every ID (see Value).
#'
#' @section Ties:
#' A location the same distance from two seeds takes the one that comes
#' first in `seeds`. For the seeds of a [select_resolution()] result that is
#' the lower-numbered centre, the rule [resolution_profile()] scores with.
#' When the summary's record of the count holds the centres the profile
#' scored and `seeds` are those centres (the seeds a selection carries, or
#' those a seeding function returns for it), in any order and any CRS, the
#' distances are taken to the centres as the record holds them, a tie goes
#' to the lower-numbered centre, and each centre gives the ID of the seed
#' that stands at it. Seeds count as the centres when they are as many and
#' each lies at a different centre, within a thousandth of its distance from
#' the next nearest one, so a seed that came back from another CRS with its
#' coordinates rounded, or seeds sorted by their cell, do not settle a tie
#' another way. A location on the edge between two such cells then takes
#' the cell [assign_features_to_polygons()] gives a row there, with `seeds`
#' and without. Other seeds, among them a subset of the centres, are
#' searched as given. A
#' location the same distance from two locations of scored rows takes the
#' one with the lower x coordinate, then the lower y; among scored rows at
#' one location that lie in different cells, the cell that comes first in
#' the order of the IDs (ascending for numbers, the C-locale order for
#' strings, the order of the levels for a factor). Neither rule depends on
#' the order of the rows of `points`, or on whether the \pkg{FNN} package is
#' installed. A tie is a tie between distances computed from the coordinates
#' in the working CRS, so on gridded locations the same layer in another CRS
#' or unit can put a location on a boundary in the other cell.
#'
#' @section Away from the observations:
#' The value and the interval describe locations among the observations.
#' The interval's two ends are quantiles of leave-one-out errors, each the
#' error of predicting an observation from the others around it, so it says
#' nothing of a location far from every observation: a gap in the layer, or
#' ground beyond its edge. Away from the observations no interval here is
#' calibrated. `dist_nearest` and `dist_ratio` say how far each location is
#' from the data. No threshold is applied to either and no location is
#' flagged or left out. [area_of_applicability()] and [fold_separation()]
#' are the tools for deciding where a prediction applies.
#'
#' A simulation shows how the interval's coverage falls with `dist_ratio`
#' (`dev/sim-pooled-values.R` in the package's repository; "Pooled values" in
#' [summarize_by_cell()] describes it). Each of 200 layers has 1,500
#' clustered observations in 333 cells, 300 further observations held out
#' at random among them and 300 held out as a block, which leaves a gap in
#' the layer. Over the held-out observations of all layers, the interval at
#' `obs_level = 0.90` held 91.5 percent of those with a `dist_ratio` below
#' 1, 89.2 percent from 1 to 2, 83.9 percent from 2 to 4, 77.9 percent from
#' 4 to 8 and 69.5 percent above 8. It held 89.6 percent of the observations
#' held out at random and 75.1 percent of those in the block.
#'
#' @section Without the FNN package:
#' The function makes three nearest-neighbour searches: each location's
#' nearest scored location, each location's nearest seed, and each scored
#' location's nearest other one (for the spacing). With \pkg{FNN} each uses
#' a kd-tree. Without it each compares every pair, and a search that would
#' compare more than `getOption("spatialkit.exact_pairs_max", 2e9)` pairs
#' stops the call first, with an error of class `spatialkit_work_limit`
#' that names the search (the fields `work`, `limit` and `what`). The search
#' for the spacing does not depend on the number of locations asked for:
#' above about 44,700 distinct scored locations the call stops without
#' \pkg{FNN} however few they are. The result is the same with and without
#' the package.
#'
#' @param newdata The locations: an sf layer. A feature that is not a point
#'   is reduced to one with [`coerce_to_points(, "auto")`][coerce_to_points].
#'   No column of it is read. With no CRS while `points` has one, the call
#'   stops.
#' @param values The [summarize_by_cell()] result for `points`. It must hold
#'   the ID column and a cell value of `response_var`: the pooled value
#'   (`estimate = "pooled"`) or the mean (the default `agg_funs`). With
#'   `seeds = NULL` it must also carry the cell polygons, which it does when
#'   it was made with `cells_sf`.
#' @param points The assigned layer `values` was computed from: an sf layer
#'   with the ID column and `response_var`, and at least one scored row.
#' @param response_var Name of the summarised column wanted: the response,
#'   or a predictor.
#' @param id_col Name of the ID column, found as [summarize_by_cell()] finds
#'   it: when `points` has no column of that name, the first of `"poly_id"`,
#'   `"polygon_id"` and `"cell_id"` it has. It cannot have the name of
#'   another column of the result (`value`, `estimate`, `obs_lo`, `obs_hi`,
#'   `filled`, `filled_from`, `dist_nearest`, `dist_ratio`).
#' @param seeds `NULL`, or the seeds the cells were built on: an sf POINT
#'   layer with a column named as the ID column that gives each seed's cell.
#'   For cells from `tess <- build_tessellation(seeds)` that is
#'   `seeds$cell_id <- tess$index` when the ID column is `"cell_id"`. Two
#'   seeds at different locations with one ID are an error; seeds that share
#'   a location share its cell, as `index` gives them. `NA` marks a seed
#'   without a cell.
#'
#' @return A data frame with one row per row of `newdata`, in its order:
#'   \describe{
#'     \item{the ID column}{The cell the location lies in, under the name of
#'       the ID column; `NA` for a location in no cell.}
#'     \item{`value`}{The cell's value, or for a filled row the value of the
#'       cell in `filled_from`.}
#'     \item{`estimate`}{`"pooled"` or `"mean"`: which value of the summary
#'       `value` is.}
#'     \item{`obs_lo`, `obs_hi`}{The interval for one observation, when
#'       `values` holds it.}
#'     \item{`filled`}{`TRUE` when the location has no cell, or its cell
#'       neither a scored row nor a value, and the value is that of the cell
#'       of the nearest scored row.}
#'     \item{`filled_from`}{For a filled row, the ID of the cell the value
#'       came from; `NA` otherwise.}
#'     \item{`dist_nearest`}{Distance to the nearest location of a scored
#'       row, in the units of the working CRS.}
#'     \item{`dist_ratio`}{`dist_nearest` over the spacing; `NA` when the
#'       scored rows sit on one location.}
#'   }
#'   A row of `newdata` with an empty geometry, or with a coordinate that is
#'   not finite or lies beyond 1e150, has `NA` in every column. The two ID
#'   columns have the class of the ID column of `values`; when a cell that
#'   `values` has no row for has an ID that class cannot hold (a string
#'   beside numbers), both are character.
#'
#'   What produced the frame is in `attr(, "values_at")`, a list of class
#'   `"spatialkit_values_at"` with a print method: `assignment` (`"seeds"`
#'   or `"polygons"`), `crs` (the working CRS), `column` (`response_var`),
#'   `estimate`, `obs_level` (the level of the interval columns; `NULL`
#'   without them, and `NA` when `values` holds them and no longer says at
#'   what level), `n` (the scored rows), `locations` (their distinct
#'   locations), `spacing` (what `dist_ratio` divides by; `NA` with one
#'   location),
#'   `filled` (how many rows of `newdata` were filled) and `count_rule` (the
#'   summary's record of how the cell count was set, or `NULL`). It describes
#'   the call: a subset of the frame's rows keeps it as it is.
#'
#'   The function raises no classed warning of its own. It passes on what
#'   [coerce_to_points()], a CRS transformation and
#'   [assign_features_to_polygons()] raise, but for one warning: the one
#'   that function raises when no location at all falls in a polygon is
#'   neither raised nor logged, since every such row has `filled = TRUE`,
#'   which says the same row by row.
#' @examples
#' library(sf)
#' set.seed(1)
#' n <- 400
#' east  <- runif(n, 0, 1e4)
#' north <- runif(n, 0, 1e4)
#' all_pts <- st_as_sf(
#'   data.frame(x = 5e5 + east, y = 5e6 + north,
#'              income = 30 + 5 * sin(east / 2000) + 5 * cos(north / 1500) +
#'                rnorm(n)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' pts    <- all_pts[1:350, ]
#' newpts <- all_pts[351:400, ]
#'
#' # The count from a size, the cells, the observations in them, the values.
#' prof  <- resolution_profile(pts, "income", estimate = "pooled",
#'                             obs_per_cell = 4.5, sac = "auto")
#' sel   <- select_resolution(prof, criterion = "size")
#' seeds <- get_voronoi_seeds(method = "kmeans", n = sel, sample_points = pts)
#' tess  <- build_tessellation(seeds, method = "voronoi")
#' asg   <- assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id")
#' out   <- summarize_by_cell(asg, "income", id_col = "cell_id",
#'                            cells_sf = tess$cells, estimate = "pooled",
#'                            obs_level = 0.90)
#'
#' # The value at other locations, each placed in the cell of its nearest
#' # seed.  `tess$index` gives the cell of each seed.
#' seeds$cell_id <- tess$index
#' at <- cell_values_at(newpts, out, asg, "income", id_col = "cell_id",
#'                      seeds = seeds)
#' head(at)
#' attr(at, "values_at")
#'
#' # Without the seeds a location is placed in the polygon that holds it.  One
#' # beyond the layer's edge lies in no polygon: it is filled from the cell of
#' # the nearest observation, and `dist_ratio` says how far that is.
#' beyond <- st_sf(geometry = st_sfc(st_point(c(5e5 + 2e4, 5e6 + 5e3)), crs = 32632))
#' cell_values_at(rbind(st_sf(geometry = st_geometry(newpts[1:2, ])), beyond),
#'                out, asg, "income", id_col = "cell_id")
#' @family aggregation
#' @seealso [summarize_by_cell()], which computes `values`;
#'   [build_tessellation()], whose `index` gives each seed's cell;
#'   [area_of_applicability()] and [fold_separation()] for where a
#'   prediction applies.
#' @export
cell_values_at <- function(newdata, values, points, response_var, id_col = "poly_id",
                           seeds = NULL) {
  # ---- the arguments ---------------------------------------------------------
  if (!inherits(newdata, "sf"))
    stop("cell_values_at(): `newdata` must be an sf object.", call. = FALSE)
  if (!inherits(points, "sf"))
    stop("cell_values_at(): `points` must be an sf object.", call. = FALSE)
  if (!is.data.frame(values))
    stop("cell_values_at(): `values` must be the summarize_by_cell() result for ",
         "`points` (a data frame).", call. = FALSE)
  if (!is.character(response_var) || length(response_var) != 1L ||
      is.na(response_var) || !nzchar(response_var))
    stop("cell_values_at(): `response_var` must be a single column name (a character ",
         "string).", call. = FALSE)
  if (!is.character(id_col) || length(id_col) != 1L || is.na(id_col) || !nzchar(id_col))
    stop("cell_values_at(): `id_col` must be a single column name (a character string).",
         call. = FALSE)

  # The ID column, found in `points` as summarize_by_cell() finds it; the
  # summary's own ID column carries that name.
  rows <- sf::st_drop_geometry(points)
  id_candidates <- unique(c(id_col, "poly_id", "polygon_id", "cell_id"))
  id_found <- id_candidates[id_candidates %in% names(rows)]
  if (!length(id_found))
    stop("cell_values_at(): could not find an ID column in `points`. Looked for: ",
         paste(id_candidates, collapse = ", "), ".", call. = FALSE)
  id_col <- id_found[[1L]]
  if (!(id_col %in% names(values)))
    stop(sprintf(paste0("cell_values_at(): `values` has no ID column '%s', the ID column ",
                        "of `points`; it must be the summarize_by_cell() result for ",
                        "`points`."), id_col), call. = FALSE)
  # The frame returned names its first column after the ID column, so that
  # name cannot be another column's.
  reserved <- c("value", "estimate", "obs_lo", "obs_hi", "filled", "filled_from",
                "dist_nearest", "dist_ratio")
  if (id_col %in% reserved)
    stop(sprintf(paste0("cell_values_at(): the ID column is named '%s', which is a column ",
                        "of the result (%s). Rename it in `points`, `values` and `seeds`."),
                 id_col, paste(reserved, collapse = ", ")), call. = FALSE)
  if (!(response_var %in% names(rows)))
    stop(sprintf("cell_values_at(): `points` has no column '%s'.", response_var),
         call. = FALSE)

  cols <- .values_at_columns(names(values), response_var)
  if (is.null(cols)) {
    # The names of the role the column has in `values`: a predictor's when the
    # summary holds another statistic of it as a predictor and none as the
    # response, the response's otherwise.
    of_it <- endsWith(names(values), paste0("_", response_var))
    p <- if (any(of_it & startsWith(names(values), "pred_")) &&
             !any(of_it & startsWith(names(values), "resp_"))) "pred_" else "resp_"
    stop(sprintf(paste0("cell_values_at(): `values` has no cell value of '%s' (neither ",
                        "`..pooled_%s` nor `%s`). Summarise it with summarize_by_cell() ",
                        "first, with the mean in `agg_funs` or with estimate = \"pooled\"."),
                 response_var, paste0(p, response_var),
                 paste0(p, "mean_", response_var)), call. = FALSE)
  }

  by_seeds <- !is.null(seeds)
  if (!by_seeds) {
    polygons <- inherits(values, "sf") && nrow(values) > 0L &&
      all(sf::st_geometry_type(values, by_geometry = TRUE) %in% c("POLYGON", "MULTIPOLYGON"))
    if (!polygons)
      stop(sprintf(paste0("cell_values_at(): `values` carries no cell polygons, so a ",
                          "location cannot be placed in a cell. Build `values` with ",
                          "`cells_sf`, or pass the seeds the cells were built on as ",
                          "`seeds` (seeds$%s <- tess$index)."), id_col), call. = FALSE)
  } else {
    if (!inherits(seeds, "sf") || !(id_col %in% names(sf::st_drop_geometry(seeds))) ||
        !all(sf::st_geometry_type(seeds, by_geometry = TRUE) == "POINT"))
      stop(sprintf(paste0("cell_values_at(): `seeds` must be an sf POINT layer with a ",
                          "column '%s' that gives each seed's cell; for cells from ",
                          "build_tessellation(seeds), set seeds$%s <- tess$index."),
                   id_col, id_col), call. = FALSE)
    seed_ids <- .id_as_character(seeds[[id_col]])
    # Seeds at different locations with one ID are refused, and counted: every
    # seed whose ID another seed elsewhere carries.  Seeds on one location
    # are not: build_tessellation() gives them one cell, and `index` repeats
    # its ID for each of them.
    shared <- unique(seed_ids[duplicated(seed_ids) & !is.na(seed_ids)])
    if (length(shared)) {
      at_seed <- sf::st_coordinates(seeds)[, 1:2, drop = FALSE]
      among   <- which(seed_ids %in% shared)
      first   <- among[match(seed_ids[among], seed_ids[among])]
      apart   <- !((at_seed[among, 1L] == at_seed[first, 1L] &
                    at_seed[among, 2L] == at_seed[first, 2L]) %in% TRUE)
      twice   <- sum(seed_ids[among] %in% seed_ids[among][apart])
      if (twice > 0L)
        stop(sprintf(paste0("cell_values_at(): %d seeds have the cell of another in ",
                            "`seeds$%s`; each seed stands for one cell. Set the column ",
                            "from the tessellation built on these seeds ",
                            "(seeds$%s <- tess$index)."),
                     twice, id_col, id_col), call. = FALSE)
    }
  }

  points_crs <- sf::st_crs(points)
  if (is.na(sf::st_crs(newdata)) && !is.na(points_crs))
    stop("cell_values_at(): `newdata` has no CRS and `points` has one; set it with ",
         "sf::st_set_crs() so the locations can be placed.", call. = FALSE)

  # ---- 1. the working CRS ----------------------------------------------------
  # The record's, else the pooled record's, each when it is a projected CRS;
  # else what the summary's own rule gives for `points`.
  rule   <- .count_rule_get(values)
  pooled <- attr(values, "pooled")
  measured <- function(crs)
    inherits(crs, "crs") && !is.na(crs) && !isTRUE(sf::st_is_longlat(crs))
  crs_given <- if (!is.null(rule) && measured(rule$crs)) rule$crs
               else if (is.list(pooled) && measured(pooled$crs)) pooled$crs
               else NULL
  at_rows <- .pooled_xy(points, crs_given)
  crs     <- at_rows$crs
  # Points without a CRS whose coordinates the step above took for lon/lat
  # have moved; a layer without a CRS beside them is read the same way.
  as_lonlat <- FALSE
  if (is.na(points_crs) && !is.na(crs) && nrow(points) > 0L) {
    raw <- .values_at_xy(suppressWarnings(.values_at_points(points)))
    as_lonlat <- isTRUE(any(at_rows$xy != raw, na.rm = TRUE))
  }

  # ---- 2. the scored rows and their distinct locations -----------------------
  y      <- rows[[response_var]]
  row_id <- rows[[id_col]]
  xy     <- at_rows$xy
  scored <- is.numeric(y) & is.finite(y) & !is.na(row_id) &
    is.finite(xy[, 1L]) & is.finite(xy[, 2L])
  n_scored <- sum(scored)
  if (n_scored == 0L)
    stop(sprintf(paste0("cell_values_at(): `points` has no row with a value of '%s', a ",
                        "cell and a location."), response_var), call. = FALSE)
  sx  <- xy[scored, 1L]
  sy  <- xy[scored, 2L]
  sid <- row_id[scored]
  # The cells in the order of their IDs, and the rows by location and then by
  # cell: the first row of a location is the one in the first of those cells.
  cells <- unique(sid)
  cells <- cells[order(cells, method = "radix")]
  gi    <- match(sid, cells)
  o     <- order(sx, sy, gi)
  sx <- sx[o]; sy <- sy[o]
  first <- c(TRUE, sx[-1L] != sx[-n_scored] | sy[-1L] != sy[-n_scored])
  ux <- sx[first]
  uy <- sy[first]
  U  <- length(ux)
  u_cell <- .id_as_character(cells)[gi[o][first]]

  # ---- the locations asked for -----------------------------------------------
  m     <- nrow(newdata)
  where <- .values_at_points(newdata)
  placed <- .values_at_place(where, "newdata", "locations", crs, points_crs, as_lonlat)
  q     <- .values_at_xy(placed)
  live  <- which(is.finite(q[, 1L]) & is.finite(q[, 2L]))
  qx <- q[live, 1L]
  qy <- q[live, 2L]

  if (by_seeds) {
    seed_xy <- .values_at_xy(.values_at_place(.values_at_points(seeds), "seeds", "seeds",
                                              crs, points_crs, as_lonlat))
    seed_at <- which(is.finite(seed_xy[, 1L]) & is.finite(seed_xy[, 2L]))
  }

  # ---- the bound on the exact searches ----------------------------------------
  # All three before any of them starts.  The third does not depend on the
  # locations asked for.
  use_fnn <- .has_backend("FNN")
  if (!use_fnn) {
    limit <- .exact_pairs_max()
    n_loc <- function(k) sprintf("%s location%s", format(k), if (k == 1L) "" else "s")
    .work_guard(as.numeric(length(live)) * U, limit, what = "nearest scored location",
                caller = "cell_values_at",
                search = sprintf("finding the nearest of %d scored location%s for %s", U,
                                 if (U == 1L) "" else "s", n_loc(length(live))))
    if (by_seeds)
      .work_guard(as.numeric(length(live)) * length(seed_at), limit, what = "nearest seed",
                  caller = "cell_values_at",
                  search = sprintf("finding the nearest of %d seed%s for %s",
                                   length(seed_at), if (length(seed_at) == 1L) "" else "s",
                                   n_loc(length(live))))
    if (U > 1L)
      .work_guard(as.numeric(U) * U, limit, what = "spacing", caller = "cell_values_at",
                  search = sprintf(paste0("the spacing of the scored rows (the nearest other ",
                                          "of each of %d scored locations)"), U))
  }

  # ---- 3. the cell of a location ----------------------------------------------
  cell <- rep(NA_character_, length(live))
  if (by_seeds) {
    if (length(live) && length(seed_at)) {
      # Seeds that are the centres the profile scored, in any order, are
      # searched at the centres as the record holds them and in its order, so
      # a location the same distance from two goes where the profile and the
      # assignment put it, whatever rounding a change of CRS brought to the
      # seeds.  Each centre takes the ID of the seed that is it.
      as_cen <- .values_at_centres(rule, seed_xy, crs, use_fnn)
      if (is.null(as_cen)) {
        if (!is.null(.scored_centres(rule)))
          .log_info(paste0("cell_values_at(): the %d seeds are not the %d centres the ",
                           "summary's record holds, so they are searched as given, and a ",
                           "location the same distance from two seeds takes the one that ",
                           "comes first in `seeds`."),
                    nrow(seed_xy), nrow(.scored_centres(rule)))
        j <- .values_at_nearest(qx, qy, seed_xy[seed_at, 1L], seed_xy[seed_at, 2L], use_fnn)
        cell <- seed_ids[seed_at][j]
      } else {
        cen <- as_cen$centres
        j <- .values_at_nearest(qx, qy, cen[, 1L], cen[, 2L], use_fnn)
        cell <- seed_ids[as_cen$seed[j]]
      }
    }
  } else if (length(live) && !all(is.na(values[[id_col]]))) {
    # The assignment is given the points and no column of `newdata`, so none
    # can collide with the ID column, and the cells' polygons with their IDs
    # and no other column or class of `values`.  A row of `values` without
    # an ID is no cell.
    feat <- where[live, , drop = FALSE]
    has_id   <- !is.na(values[[id_col]])
    cells_at <- sf::st_sf(
      stats::setNames(data.frame(values[[id_col]][has_id], stringsAsFactors = FALSE), id_col),
      geometry = sf::st_geometry(values)[has_id])
    # With the summary's record of the count, which the rows were assigned
    # with: a location on the edge between two cells the profile scored goes
    # where the profile put such a row.
    cells_at <- .count_rule_put(cells_at, rule)
    # `frame` is the CRS the coordinates of `points` are read in, which is
    # where the rows were assigned from.  Points without a CRS are given it
    # beside polygons that carry one, and not left for the assignment to
    # stamp with the polygons' CRS.  Polygons without a CRS take the CRS of
    # the features handed to that function, so the points go there in
    # `frame`, as the rows did, whatever CRS `newdata` arrived in.
    frame <- if (!is.na(points_crs)) points_crs
             else if (as_lonlat) sf::st_crs(4326) else crs
    if (is.na(sf::st_crs(feat))) {
      if (!is.na(sf::st_crs(cells_at)) && !is.na(frame)) feat <- sf::st_set_crs(feat, frame)
    } else if (is.na(sf::st_crs(cells_at)) && !is.na(frame) && sf::st_crs(feat) != frame) {
      feat <- sf::st_transform(feat, frame)
    }
    # That function's warning that no feature falls in any polygon is neither
    # raised nor logged: every such row is returned with filled = TRUE.  Its
    # other warnings pass.
    answered <- .answer_for_classed("spatialkit_unassigned")
    held <- tryCatch(
      assign_features_to_polygons(feat, cells_at, polygon_id_col = id_col,
                                  keep_unassigned = TRUE, quiet = TRUE),
      finally = answered())
    if (nrow(held) != length(live))
      stop("cell_values_at(): the assignment returned ", nrow(held), " rows for ",
           length(live), " locations.", call. = FALSE)
    cell <- .id_as_character(held[[id_col]])
  }

  # ---- 5 and 6. the nearest scored row, and how far it is ----------------------
  near <- if (length(live)) .values_at_nearest(qx, qy, ux, uy, use_fnn) else integer(0)
  dist <- sqrt((qx - ux[near])^2 + (qy - uy[near])^2)
  spacing <- NA_real_
  if (U > 1L) {
    other   <- .nearest_location_search(cbind(ux, uy), seq_len(U), use_fnn,
                                        caller = "cell_values_at")
    spacing <- stats::median(sqrt((ux - ux[other])^2 + (uy - uy[other])^2))
  }

  # ---- 4 and 5. the value, the cell's own or the nearest scored row's ----------
  value_ids <- .id_as_character(values[[id_col]])
  val <- as.numeric(values[[cols$value]])
  row_of <- function(id) match(id, value_ids, incomparables = NA)
  # Filled: no cell, or a cell with no value in `values` and no scored row.
  # A cell that holds a scored row keeps what `values` has for it, be it NA.
  filled <- is.na(val[row_of(cell)]) & !(cell %in% .id_as_character(cells))
  from   <- ifelse(filled, u_cell[near], NA_character_)
  taken  <- row_of(ifelse(filled, from, cell))

  # ---- the frame: one row per row of `newdata` ---------------------------------
  whole <- function(x, na) { out <- rep(na, m); out[live] <- x; out }
  ids   <- .values_at_ids(c(whole(cell, NA_character_), whole(from, NA_character_)),
                          values[[id_col]])
  has_interval <- !is.null(cols$lo)
  out <- data.frame(id = ids[seq_len(m)],
                    value = whole(val[taken], NA_real_),
                    estimate = whole(rep(cols$estimate, length(live)), NA_character_),
                    stringsAsFactors = FALSE)
  if (has_interval) {
    out$obs_lo <- whole(as.numeric(values[[cols$lo]])[taken], NA_real_)
    out$obs_hi <- whole(as.numeric(values[[cols$hi]])[taken], NA_real_)
  }
  out$filled       <- whole(filled, NA)
  out$filled_from  <- ids[m + seq_len(m)]
  out$dist_nearest <- whole(dist, NA_real_)
  out$dist_ratio   <- whole(dist / spacing, NA_real_)
  names(out)[1L] <- id_col

  obs_level <- if (!has_interval) NULL
               else if (is.list(pooled) && is.numeric(pooled$obs_level) &&
                        length(pooled$obs_level) == 1L) pooled$obs_level
               else NA_real_
  attr(out, "values_at") <- structure(
    list(assignment = if (by_seeds) "seeds" else "polygons",
         crs        = crs,
         column     = response_var,
         estimate   = cols$estimate,
         obs_level  = obs_level,
         n          = n_scored,
         locations  = U,
         spacing    = spacing,
         filled     = sum(filled),
         count_rule = rule),
    class = "spatialkit_values_at")
  out
}


#' @export
print.spatialkit_values_at <- function(x, ...) {
  label  <- function(name) sprintf("  %-12s: ", name)
  indent <- strrep(" ", 16L)
  count  <- function(k, what)
    sprintf("%s %s%s", format(k), what, if (isTRUE(k == 1L)) "" else "s")
  cat(sprintf("Cell values at other locations: %s, the %s\n", x$column,
              if (identical(x$estimate, "pooled")) "pooled value" else "plain mean"))
  cat(label("cell"),
      if (identical(x$assignment, "seeds"))
        "of the nearest seed (the Voronoi cell before its clip to the cover)"
      else "the polygon of `values` that holds the location",
      "\n", sep = "")
  if (!is.null(x$obs_level))
    cat(label("interval"),
        if (is.numeric(x$obs_level) && length(x$obs_level) == 1L && !is.na(x$obs_level))
          paste0(.level_text(x$obs_level), "% ")
        else "",
        "for ONE OBSERVATION (obs_lo, obs_hi), centred on the value\n", sep = "")
  crs <- .fold_crs_label(x$crs)
  cat(label("scored rows"), count(x$n, "row"), " at ", count(x$locations, "location"),
      if (is.na(crs)) ", coordinates without a CRS" else paste0(", measured in ", crs),
      "\n", sep = "")
  cat(label("spacing"),
      if (is.na(x$spacing))
        "not available: the scored rows sit on one location, so dist_ratio is NA"
      else paste0(format(signif(x$spacing, 4)), ", the median distance from a scored ",
                  "location to the\n", indent, "nearest other; dist_ratio is ",
                  "dist_nearest over it"),
      "\n", sep = "")
  cat(label("filled"),
      if (isTRUE(x$filled == 0L)) "no location"
      else paste0(count(x$filled, "location"), if (isTRUE(x$filled == 1L)) " has" else " have",
                  " no cell with a value and ",
                  if (isTRUE(x$filled == 1L)) "takes" else "take", " that of\n", indent,
                  "the cell of the nearest scored row (filled_from)"),
      "\n", sep = "")
  rule <- x$count_rule
  if (is.list(rule) && length(rule$cells) == 1L && length(rule$rule) == 1L)
    cat(label("count"), count(rule$cells, "cell"), ", set by ", rule$rule,
        " (the summary's record, in `count_rule`)\n", sep = "")
  if (!is.null(x$obs_level))
    cat("Away from the observations no interval here is calibrated.\n")
  invisible(x)
}
