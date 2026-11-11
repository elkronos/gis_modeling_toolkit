# =============================================================================
# Prediction to a surface
# =============================================================================

#' Build a regular point grid covering a bounding box
#'
#' Cell centres, not corners, so every point lies inside the cell it represents.
#'
#' @param bb An \code{sf} bbox.
#' @param crs Target CRS for the result.
#' @param cell_size Edge length in CRS units, or NULL to derive from n_cells.
#' @param n_cells Approximate number of cells when \code{cell_size} is NULL.
#' @return An sf POINT layer with columns \code{..grid_x}, \code{..grid_y}.
#' @keywords internal
#' @noRd
# Ceiling on the number of prediction cells.  A surface of more than a few
# million points is not something anyone plots or writes out; a request for one
# is a units mistake, and refusing it in milliseconds beats failing to allocate
# after minutes.
.surface_max_cells <- 5e6

.make_prediction_grid <- function(bb, crs, cell_size = NULL, n_cells = 10000L) {
  w <- as.numeric(bb[["xmax"]] - bb[["xmin"]])
  h <- as.numeric(bb[["ymax"]] - bb[["ymin"]])
  if (!is.finite(w) || !is.finite(h) || w <= 0 || h <= 0)
    stop("predict_surface(): the fitted data has no finite extent.", call. = FALSE)

  if (is.null(cell_size)) {
    # as.integer() silently returns NA above 2^31, and the error below then
    # blamed `cell_size` -- an argument the caller did not pass.  Validate the
    # argument the caller actually gave.
    if (!is.numeric(n_cells) || length(n_cells) != 1L || !is.finite(n_cells) ||
        n_cells < 1)
      stop("predict_surface(): `n_cells` must be a single positive finite number.",
           call. = FALSE)
    if (n_cells > .surface_max_cells)
      stop(sprintf(paste0("predict_surface(): `n_cells` = %s is above the %s ",
                          "cells this function will build."),
                   format(n_cells, big.mark = ",", scientific = FALSE),
                   format(.surface_max_cells, big.mark = ",",
                          scientific = FALSE)), call. = FALSE)
    cell_size <- sqrt((w * h) / n_cells)
  }
  if (!is.numeric(cell_size) || length(cell_size) != 1L ||
      !is.finite(cell_size) || cell_size <= 0)
    stop("predict_surface(): `cell_size` must be a positive finite number.",
         call. = FALSE)

  # Refuse an absurd grid before expand.grid() builds it.  The classic units
  # mistake -- a cell_size in kilometres on a metre CRS -- asked for 3.7e10
  # cells and R aborted trying to allocate 139.6 GB (or swapped the machine on
  # a bigger one).  create_grid_polygons() guards the identical mistake and
  # names the CRS units in its message; this builder is the sibling that did
  # not.
  n_est <- ceiling(w / cell_size) * ceiling(h / cell_size)
  if (is.finite(n_est) && n_est > .surface_max_cells) {
    unit_lbl <- tryCatch({
      u <- sf::st_crs(crs)$units_gdal
      if (is.null(u) || is.na(u) || !nzchar(u)) "CRS units" else u
    }, error = function(e) "CRS units")
    stop(sprintf(paste0("predict_surface(): a cell size of %s on an extent of ",
                        "%s x %s would produce about %s cells, above the %s ",
                        "this function will build. Check that `cell_size` is in ",
                        "the CRS units of the fitted data (%s), or raise ",
                        "`n_cells` instead."),
                 format(signif(cell_size, 6)), format(signif(w, 4)),
                 format(signif(h, 4)),
                 format(n_est, big.mark = ",", scientific = FALSE),
                 format(.surface_max_cells, big.mark = ",", scientific = FALSE),
                 unit_lbl),
         call. = FALSE)
  }

  # seq(from, to, by) ERRORS with "wrong sign in 'by' argument" when from > to,
  # rather than returning length 0 -- so a cell_size wider than the extent used
  # to abort with an internal-looking message, and the length-0 guards below
  # were unreachable.  Build the axis explicitly instead.
  .axis <- function(lo, hi) {
    lo <- as.numeric(lo); hi <- as.numeric(hi)
    # Enough cells to cover the extent, centred on it.  floor() cells anchored
    # at the lower bound left the remainder -- up to a cell wide -- uncovered
    # on the east and north: cell_size = 100 on a 980 x 956 extent covered
    # 900 x 900, and 14 of 120 training points lay in no cell.  The grid now
    # overhangs the box by less than one cell, split evenly on both sides, so
    # every centre stays inside it.  The relative tolerance keeps an exact
    # multiple exact: 0.3 / 0.1 is 2.9999999999999996, and a ratio a hair
    # above an integer must not gain a column either.
    r <- (hi - lo) / cell_size
    n <- ceiling(r - sqrt(.Machine$double.eps) * max(1, r))
    if (!is.finite(n) || n < 1L) n <- 1L     # a cell wider than the extent
    offset <- (hi - lo - n * cell_size) / 2
    lo + offset + cell_size / 2 + seq.int(0L, n - 1L) * cell_size
  }
  xs <- .axis(bb[["xmin"]], bb[["xmax"]])
  ys <- .axis(bb[["ymin"]], bb[["ymax"]])

  g <- expand.grid(..grid_x = xs, ..grid_y = ys, KEEP.OUT.ATTRS = FALSE)
  out <- sf::st_as_sf(g, coords = c("..grid_x", "..grid_y"), crs = crs,
                      remove = FALSE)
  attr(out, "cell_size") <- cell_size
  out
}


#' The covariate row nearest each grid point
#'
#' A grid point between covariate points of a lattice (the grid of
#' \code{predict_surface()} laid over the training points, say) is the same
#' distance from several of them.  \code{sf::st_nearest_feature()} gave it
#' whichever one GEOS met first, and the distances differ in their last bits
#' once the unit or the origin of the coordinates changes, so the surface
#' took another covariate row in another frame.  For point covariates the
#' search is that of \code{cell_values_at()} (\code{.values_at_nearest()}):
#' squared distances within \code{.TIE_TOL} (one part in ten million) of the
#' least, plus the margin \code{.centres_tol()} of the points on the
#' distance (\code{.tie_bound()}), are a tie, and the tie goes to the point
#' in the lowest column
#' (\code{.near_columns()}: first coordinates within \code{.TIE_TOL} of the
#' extent are one column), then in the lowest row by the same rule read on
#' y, then the lower row of the layer.  The extent is the diagonal of the
#' points' box, or the margin of an edge at their largest coordinate
#' (\code{.edge_tol()}) over \code{.TIE_TOL} where that is more
#' (\code{.near_points_extent()}): two
#' records of one station that differ in their last bits, alone in the
#' layer, have a box whose diagonal is last-bits small, and a raw y read
#' after the column let those bits choose, which move with the unit and the
#' origin.  A layer of lines or polygons (or one that mixes them with
#' points) is searched by \code{.surface_nearest_shapes()}, the same rule
#' read on distances to the features and on the centres of their bounding
#' boxes.  A layer in degrees is searched by
#' \code{sf::st_nearest_feature()}.  A grid point without a location gets \code{NA} whatever the
#' layer, where \code{st_nearest_feature()} gave a line, polygon or point in
#' degrees a row.
#'
#' @param pts The grid's points, an sfc in the CRS of \code{covariates}.
#' @param covariates The covariate layer.
#' @return Integer vector, one row of \code{covariates} per point; \code{NA}
#'   for a point without a location and where no covariate row has one.
#' @keywords internal
#' @noRd
.surface_nearest_rows <- function(pts, covariates) {
  # Points in an sfc_GEOMETRY column are read as the same points in
  # sfc_POINT (.typed_points()): sf::st_coordinates() refuses that class.
  pts <- .typed_points(pts)
  cg  <- .typed_points(sf::st_geometry(covariates))
  if (!isTRUE(sf::st_is_longlat(cg)) &&
      !all(sf::st_geometry_type(cg, by_geometry = TRUE) == "POINT"))
    return(.surface_nearest_shapes(pts, cg))
  out  <- rep(NA_integer_, length(pts))
  gxy  <- .values_at_xy(sf::st_sf(geometry = pts))
  q    <- which(!is.na(gxy[, 1L]))
  if (!length(q)) return(out)
  if (isTRUE(sf::st_is_longlat(cg))) {
    out[q] <- sf::st_nearest_feature(pts[q], cg)
    return(out)
  }
  cxy  <- .values_at_xy(sf::st_sf(geometry = cg))
  live <- which(!is.na(cxy[, 1L]))
  if (!length(live)) return(out)
  cx <- cxy[live, 1L]; cy <- cxy[live, 2L]
  # By column, then by row read the same way on y, then row of the layer:
  # order() keeps the row order among ties.  The step is never less than the
  # margin of an edge (.near_points_extent(), .edge_tol()): points
  # that all lie within the last bits of each other have a last-bits small
  # box.
  ex <- .near_points_extent(cx, cy)
  o  <- order(.near_columns(cx, cy, ex), .near_columns(cy, cx, ex))
  j <- .values_at_nearest(gxy[q, 1L], gxy[q, 2L], cx[o], cy[o],
                          use_fnn = .has_backend("FNN"))
  out[q] <- live[o[j]]
  out
}


#' The line or polygon nearest each grid point
#'
#' A grid point on the edge two covariate polygons share (or at the shared
#' end of two lines) is at distance 0 from both, and
#' \code{sf::st_nearest_feature()} gave it whichever one GEOS met first,
#' which moved with the last bits of the coordinates and with the row order
#' of the layer.  Here the least distance \code{d} is
#' \code{st_nearest_feature()}'s, and the features within
#' \code{d * (1 + .TIE_TOL) + .edge_tol()} of the point (the margin of the
#' layer's bounding box, \code{.edge_join()}'s) are a tie.  A point with one
#' such feature takes it.  Among several, the one whose
#' bounding-box centre is in the lowest column (\code{.near_columns()} over
#' the centres of the layer's features that have a location, with centre x
#' values within \code{.TIE_TOL} of the diagonal of the layer's bounding box
#' one column), then in the lowest row by the same rule read on y
#' (\code{.near_columns(y, x)}), then in the lower row of the layer, is
#' taken.  The extent is the layer's, not the centres' own: where every
#' centre lies on one point (a polygon with a hole and the island that fills
#' it) the centres' box is last-bits small, and the tie would go by those
#' bits.  The features near a point are shortlisted by the square of
#' half-side a hair above that reach around it
#' (\code{st_intersects()} uses GEOS's index, where
#' \code{st_is_within_distance()} takes a distance for every pair of a
#' projected layer), then by the distance from the point to their bounding
#' box, which is never more than the distance to the feature (a point far
#' from the layer has a large square), and kept by their distance.  A point
#' without a location (an empty point, or a coordinate that is not finite or
#' lies beyond \code{.POOL_XY_MAX}) gets \code{NA}.
#'
#' @param pts The grid's points, an sfc in the CRS of \code{cg}.
#' @param cg The covariates' geometry, not in degrees.
#' @return Integer vector, one feature of \code{cg} per point; \code{NA}
#'   for a point without a location and where \code{st_nearest_feature()}
#'   gives none.
#' @keywords internal
#' @noRd
.surface_nearest_shapes <- function(pts, cg) {
  nn <- sf::st_nearest_feature(pts, cg)
  pxy <- .values_at_xy(sf::st_sf(geometry = pts))
  nn[is.na(pxy[, 1L])] <- NA_integer_
  q <- which(!is.na(nn))
  if (!length(q)) return(nn)
  d   <- as.numeric(sf::st_distance(pts[q], cg[nn[q]], by_element = TRUE))
  tol <- .edge_tol(sf::st_bbox(cg))
  r   <- d * (1 + .TIE_TOL) + tol
  # The square holds every feature within r of the point; its corners are
  # computed in floating point, so it is taken a hair wider.
  sq   <- sf::st_buffer(pts[q], r * (1 + 1e-6) + tol, endCapStyle = "SQUARE")
  hits <- sf::st_intersects(sq, cg)
  many <- which(lengths(hits) > 1L)
  if (!length(many)) return(nn)

  # The bounding box of each feature (NA for an empty one), for its centre
  # and for the box distance below.
  live <- which(!sf::st_is_empty(cg))
  bbm  <- matrix(NA_real_, length(cg), 4L)
  bbm[live, ] <- matrix(vapply(cg[live], function(g) as.numeric(sf::st_bbox(g)),
                               numeric(4)), ncol = 4L, byrow = TRUE)
  # The rank of each feature by its box centre: column, then row by the
  # same rule on y, then row of the layer (order() keeps it among ties).
  # The step is a part of the layer's diagonal: the centres' own box is
  # last-bits small where every centre lies on one point (a polygon with a
  # hole and the island that fills it).
  cx <- (bbm[live, 1L] + bbm[live, 3L]) / 2
  cy <- (bbm[live, 2L] + bbm[live, 4L]) / 2
  bb <- as.numeric(sf::st_bbox(cg))
  ex <- sqrt((bb[3L] - bb[1L])^2 + (bb[4L] - bb[2L])^2)
  rank <- rep(NA_integer_, length(cg))
  rank[live[order(.near_columns(cx, cy, ex), .near_columns(cy, cx, ex))]] <-
    seq_along(live)

  pk <- rep(many, lengths(hits[many]))
  fj <- unlist(hits[many], use.names = FALSE)
  # The distance to a feature's box is never more than the distance to the
  # feature, so a pair whose box lies beyond the reach is not a tie.  Far
  # from the layer the square holds many features that are not near.
  px <- pxy[q[pk], 1L]; py <- pxy[q[pk], 2L]
  bx <- pmax(bbm[fj, 1L] - px, 0, px - bbm[fj, 3L])
  by <- pmax(bbm[fj, 2L] - py, 0, py - bbm[fj, 4L])
  near <- sqrt(bx^2 + by^2) <= r[pk] * (1 + 1e-6) + tol
  pk <- pk[near]; fj <- fj[near]
  two <- pk %in% pk[duplicated(pk)]
  pk <- pk[two]; fj <- fj[two]
  if (!length(pk)) return(nn)
  dj <- as.numeric(sf::st_distance(pts[q[pk]], cg[fj], by_element = TRUE))
  ok <- dj <= r[pk]
  pk <- pk[ok]; fj <- fj[ok]
  # The lowest rank among the features each point ties with.
  o  <- order(pk, rank[fj])
  first <- !duplicated(pk[o])
  nn[q[pk[o][first]]] <- fj[o][first]
  nn
}


#' Predict a fitted model onto a regular grid
#'
#' Builds a prediction surface over the extent of the training data (or over a
#' grid you supply), predicts in chunks, and returns an \code{sf} layer.
#'
#' \code{predict_surface()} builds the grid, joins the covariates, aligns the
#' CRSs and calls \code{predict()} in chunks of \code{chunk_size} rows.  Rows
#' do not interact, so chunking changes memory use and not the predictions.
#' For a \code{bayesian_fit} the draw matrix is \code{n_draws} by the number of
#' rows, which on a fine grid outgrows memory long before the fit does.
#'
#' @section Coordinate reference systems and geometry:
#' \code{grid}, \code{boundary} and \code{covariates} are brought into the
#' fit's CRS first.  A layer without a CRS is read the way the training data
#' were (the assumption recorded on the fit), with a warning, and then
#' reprojected; when that reading is longitude and latitude, a layer with an
#' x outside -180 to 360 or a y beyond 90 north or south cannot be in it, and
#' the call stops.  Without this a CRS-less grid could land thousands of
#' kilometres from the covariates, and every cell would take the same nearest
#' feature.  A layer still without a CRS after that is taken as EPSG:4326 and
#' reprojected when its coordinates look like longitude and latitude, and
#' otherwise stamped with the fit's CRS, with a warning either way; a fit's
#' CRS in longitude and latitude is not stamped on coordinates outside
#' longitude -180 to 360 and latitude -90 to 90, which is refused.  See
#' \link{spatialkit-geometry} for these rules in general.
#'
#' A grid of polygons (\code{\link{create_grid_polygons}()} output, say) is
#' reduced to one representative point per cell, as
#' \code{\link{coerce_to_points}()} does, so covariates are taken at the
#' location predicted for; \code{boundary} then keeps the cells whose point
#' falls inside it.  All three layers are read in x and y: one whose
#' geometries mix XY with XYZ or XYM, or carry M (XYM or XYZM), is read
#' without its Z and M, with no warning for a boundary or covariates.  A
#' POINT grid whose coordinates are not moved (one in
#' the fit's CRS, or one without a CRS that is stamped with it) is returned
#' with its rows as given, Z and M included; one that mixes dimensions and is
#' reprojected into the fit's CRS comes back in XY with a warning of class
#' \code{"spatialkit_zm_dropped"}.  The points a grid of polygons is reduced
#' to are in XY when its geometries mix dimensions or carry M.
#'
#' A grid row without a location (an empty geometry, or a coordinate that is
#' not finite or lies beyond 1e150) gets an \code{NA} prediction; in a grid
#' that is reprojected, or reduced from polygons, it comes back as an empty
#' point.  A covariate row without a location is never the nearest feature.
#'
#' @section Ties between nearest covariate features:
#' A grid point the same distance from several features of
#' \code{covariates} takes one of them by a fixed rule, so a grid laid
#' between the points of a lattice takes the same covariates in another unit
#' or with another origin, and with or without the \pkg{FNN} package.
#' \itemize{
#'   \item Points.  Distances count as equal when the squared distance is no
#'     more than the square of the least distance lengthened by a margin,
#'     plus one part in ten million of the least squared distance.  Among
#'     them the point lowest in x is taken, then the one lowest in y, then
#'     the first in row order.  x values count as one when they lie within
#'     one part in ten million of the diagonal of the points' bounding box,
#'     or within the margin where that is more; y values are read the same
#'     way.  A grid point near two points that lie within the margin of each
#'     other (two records of one station that differ only in the last bits
#'     of their coordinates) takes the same one however near it lies.
#'   \item Lines and polygons, or a layer that mixes them with points (on the
#'     edge two polygons share, say).  Distances within one part in ten
#'     million of the least, plus the margin, count as the same.  Among them
#'     the feature whose bounding-box centre is lowest in x is taken, then
#'     lowest in y, x values and y values within one part in ten million of
#'     the diagonal of the bounding box of \code{covariates} counting as
#'     one, then the first in row order.
#'   \item A layer in degrees.  The nearest feature is
#'     \code{sf::st_nearest_feature()}'s, which settles a tie in the order
#'     its search meets the features (GEOS, or s2 while
#'     \code{sf::sf_use_s2()} is on).
#'   \item The margin.  It is one part in 1e10 of the largest coordinate of
#'     the features that have a location, at most a millionth of the
#'     diagonal of their bounding box and at least four times the spacing of
#'     doubles at that coordinate.  It grows with the largest coordinate:
#'     about 0.46 mm at a UTM northing near 4,550 km, about half a
#'     micrometre once the origin lies next to data 5 km across.  So a grid
#'     point whose distances to two features differ by less than the margin
#'     with one origin and by more with another takes the one the tie rule
#'     picks with the first and the nearer with the second.
#' }
#'
#' @param object A \code{spatial_fit} (from \code{\link{fit_rf_model}()},
#'   \code{\link{fit_gwr_model}()} or
#'   \code{\link{fit_bayesian_spatial_model}()}).
#' @param grid Optional \code{sf} POINT layer to predict onto, with at least
#'   one row.  When \code{NULL}, a regular grid is built over the training
#'   extent.  A grid of polygons is reduced to points; see "Coordinate
#'   reference systems and geometry".
#' @param cell_size Grid resolution in CRS units.  Ignored when \code{grid} is
#'   supplied; when \code{NULL}, derived from \code{n_cells}.  A value that
#'   would produce more than 5,000,000 cells is refused, naming the implied
#'   count and the CRS units.  The usual cause is a value in the wrong unit.
#'   The grid is centred on the training bounding box and covers it: when the
#'   extent is not a whole number of cells, it overhangs the box by less than
#'   one cell, split evenly between the two sides.  A \code{cell_size} wider
#'   than the extent yields a single centred cell.
#' @param n_cells Approximate cell count used to derive \code{cell_size}.
#'   Default 10000.  Must be a single positive finite number and at most
#'   5,000,000; anything else is an error.  Also ignored when \code{grid} is
#'   supplied.  The grid you pass is used verbatim.
#' @param boundary Optional polygonal \code{sf}/\code{sfc}; grid points
#'   outside it are dropped.  Aligned to the fit's CRS as \code{grid} is.
#' @param covariates Optional \code{sf} layer carrying the model's
#'   predictors.  Required when the model has predictors and \code{grid} does
#'   not already contain them.  Values are taken from the nearest feature
#'   (see "Ties between nearest covariate features").  Aligned to the fit's
#'   CRS as \code{grid} is, with the same warning when it has no CRS.
#' @param chunk_size Rows per prediction call.  Default 5000.  A
#'   performance setting only: rows do not interact, and for a \code{bayesian_fit}
#'   the GP boundary is held at its fitted value whatever the chunk holds
#'   (see \link{spatial_fit-predict}).  The two warnings in which
#'   \code{predict()} on a GWR fit counts locations (classes
#'   \code{"spatialkit_gwr_extrapolation"} and
#'   \code{"spatialkit_gwr_not_estimable"}) are raised once for the grid,
#'   with its totals, whatever the number of chunks.
#' @param se Logical; also return a standard error (or posterior SD) column
#'   where the backend supports it.  Default FALSE.  For a \code{bayesian_fit} this is
#'   the SD of the posterior draws \code{predict()} returns, and those are of
#'   the expected value by default (\code{type = "epred"}): the uncertainty
#'   of the mean surface, not of a new observation, which also carries the
#'   observation noise.  For the predictive SD, the one that goes with
#'   prediction intervals and \code{cv_bayes()}'s calibration, pass
#'   \code{type = "predict"} as well.
#' @param ... Passed to \code{predict()}, e.g. \code{type = "predict"} for a
#'   \code{bayesian_fit}.  Not \code{draws}, which this function sets itself
#'   and refuses here.
#' @return An \code{sf} POINT layer with a \code{.pred} column (and
#'   \code{.pred_se} when \code{se = TRUE} and available; one a supplied
#'   \code{grid} already carried, from an earlier surface, is removed
#'   otherwise).  For an
#'   auto-generated grid the resolution is attached as attribute
#'   \code{"cell_size"}.  For a user-supplied \code{grid} it is only whatever
#'   \code{"cell_size"} attribute that object already carried.  That is usually
#'   \code{NULL}, and \code{NULL} for certain if the grid had to be
#'   re-projected, since \code{st_transform()} does not preserve custom
#'   attributes.  The resolution of a grid you built is not this function's to
#'   infer.
#' @seealso \code{\link{area_of_applicability}()} to mark where the surface
#'   applies, and \link{spatial_fit-predict} for the predictions of each
#'   backend.
#' @family prediction
#' @examplesIf requireNamespace("ranger", quietly = TRUE)
#' # Any spatial_fit works here; a forest keeps the example free of the
#' # optional GWR and Stan backends.
#' library(sf)
#' set.seed(1)
#' n <- 120
#' pts <- st_as_sf(
#'   data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
#'              elev = rnorm(n)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' pts$price <- 10 + 0.01 * (st_coordinates(pts)[, 1] - 5e5) +
#'   2 * pts$elev + rnorm(n)
#' fit  <- fit_rf_model(pts, "price", "elev", num_trees = 100, seed = 1)
#' surf <- predict_surface(fit, n_cells = 500, covariates = pts)
#' surf[".pred"]               # one prediction per grid cell, as an sf layer
#' # Check where the surface extrapolates before mapping it.  The grid took
#' # its covariates from the nearest observation, so here nothing is
#' # outside; covariates from another source, such as a raster converted to
#' # sf points, are where this bites.
#' area_of_applicability(surf, model = fit)
#' @export
predict_surface <- function(object, grid = NULL, cell_size = NULL,
                            n_cells = 10000L, boundary = NULL,
                            covariates = NULL, chunk_size = 5000L,
                            se = FALSE, ...) {
  if (!inherits(object, "spatial_fit"))
    stop("predict_surface(): `object` must be a spatial_fit.", call. = FALSE)

  # `draws` is this function's to set: it asks for the draw matrix itself when
  # se = TRUE.  Passed through `...` it reached predict() too, and a backend
  # that honours it returned an n_draws x n matrix that as.numeric() flattened
  # into .pred column by column -- cell 1's draws, then cell 2's -- with only
  # a cryptic length warning; with se = TRUE the duplicated argument failed
  # inside try() and was reported as a backend without draws.  A prefix
  # counts, since R's argument matching would complete it.
  dot_nms <- names(list(...))
  if (!is.null(dot_nms) && any(nzchar(dot_nms) & startsWith("draws", dot_nms)))
    stop("predict_surface(): `draws` cannot be passed through `...`; ",
         "predict_surface() requests the posterior draws itself when se = TRUE ",
         "and returns their SD as .pred_se. For the draw matrix, call ",
         "predict(object, newdata = grid, draws = TRUE).", call. = FALSE)

  train <- object$data_sf
  if (!inherits(train, "sf"))
    stop("predict_surface(): the fit carries no training geometry.", call. = FALSE)
  target_crs <- sf::st_crs(train)

  # A `grid` or `covariates` layer still without a CRS once the fit's own
  # assumption has been replayed is aligned to the fit's CRS as `boundary`
  # is, and as the other functions align such a layer: reprojected from
  # EPSG:4326 when it looks like lon/lat, otherwise stamped, either way with
  # an R warning naming this function and the argument.  ensure_projected()
  # stamped it with a log line naming neither, which knitr,
  # spatialkit_quiet() and tryCatch() never show -- and a wrongly stamped
  # layer puts every covariate lookup in the wrong place.  A layer the
  # replay marked as belonging to a CRS-less fit's own space is left there.
  .align_to_fit <- function(x, what) {
    if (is.na(sf::st_crs(x)) && !is.null(.crs_or_null(target_crs)) &&
        !identical(attr(x, "crs_assumed"), "none"))
      .transform_or_stamp(x, target_crs, what = what, caller = "predict_surface")
    else
      ensure_projected(x, target_crs = .crs_or_null(target_crs))
  }

  # ---- grid ----------------------------------------------------------------
  if (is.null(grid)) {
    grid <- .make_prediction_grid(sf::st_bbox(train), target_crs,
                                  cell_size = cell_size, n_cells = n_cells)
  } else {
    if (!inherits(grid, "sf"))
      stop("predict_surface(): `grid` must be an sf object.", call. = FALSE)
    # The boundary-clip path below guards its own empty result; without this
    # an empty `grid` reaches seq(1L, 0L, by = chunk_size) and aborts with
    # "wrong sign in 'by' argument", which names nothing the caller passed.
    if (nrow(grid) == 0L)
      stop("predict_surface(): `grid` has no rows; there is nothing to ",
           "predict onto.", call. = FALSE)
    # The same replay `covariates` gets, and for exactly the same reason: a
    # CRS-less grid whose OWN bounding box does not trip the lon/lat heuristic
    # had the fit's projected CRS stamped onto raw degrees, landing 5,300 km
    # from the correctly-reprojected covariates.  st_nearest_feature() then
    # handed every cell the same covariate row and the whole surface collapsed
    # to one constant -- silently, with no error anywhere.
    # A grid that mixes XY and XYZ (or XYM) geometries is read in x and y:
    # sf reads a layer in the dimension of its first geometry, so a move
    # into the fit's CRS stopped with "OGR error", and GEOS stopped on its
    # polygons ("ParseException", "Points of LinearRing do not form a closed
    # linestring") or on any that carry M.  A POINT grid is returned as the
    # surface: one whose coordinates are moved into the fit's CRS below
    # comes back in XY with the warning of class spatialkit_zm_dropped, as
    # ensure_projected() returns such a layer, under this function's name.
    # One whose coordinates stay as they are (already in that CRS, or
    # without a CRS and stamped with it or left in a CRS-less fit's space)
    # keeps its rows as given, as area_of_applicability() keeps stamped
    # `newdata`, and the joins below read a copy of it in XY.  The points a
    # grid of polygons is reduced to are this function's own, so that grid,
    # and one that carries M, are read in XY without a word.
    grid_moves <- function(g) {
      gcrs <- sf::st_crs(g)
      if (!is.na(gcrs)) {
        if (identical(gcrs, target_crs)) return(FALSE)
        # ensure_projected() with no target moves only a lon/lat layer.
        if (is.na(target_crs)) return(isTRUE(sf::st_is_longlat(g)))
        return(TRUE)
      }
      if (is.na(target_crs)) return(FALSE)
      # The replay below gives a CRS-less grid the fit's assumption, which
      # ensure_projected() then moves into the fit's CRS; without one,
      # .align_to_fit() reprojects it when it looks like lon/lat and stamps
      # it otherwise.
      assumed <- attr(train, "crs_assumed")
      if (!is.null(assumed) && !identical(assumed, "none"))
        return(!identical(sf::st_crs(assumed), target_crs))
      isTRUE(.looks_like_lonlat(g)$lonlat)
    }
    grid_points <- all(sf::st_geometry_type(grid, by_geometry = TRUE) == "POINT")
    moves <- grid_points && grid_moves(grid)
    if (!grid_points)
      grid <- .xy_if_mixed_or_m(grid)
    else if (moves)
      grid <- .xy_if_mixed_warn(
        grid, "predict_surface", "grid",
        why = paste0("so it is read in x and y before it is brought into ",
                     "the fit's CRS: Z and M are dropped and the surface is ",
                     "returned in XY"))
    # A row with a coordinate that is not finite or lies beyond .POOL_XY_MAX
    # has no location.  It is emptied as given, before the move into the
    # fit's CRS, which can bring it back finite (a northing of 1e200 in
    # EPSG:32616 lands at a finite point in EPSG:32617), and before a grid of
    # polygons is reduced to points; such a row was predicted at the point
    # the move gave it.  It comes back empty with an NA prediction, as the
    # same row emptied.  A POINT grid that stays where it is keeps the row as
    # given, and predict() gives it NA.
    if (!grid_points || moves) grid <- .empty_unlocated(grid)
    grid <- .replay_crs_assumption(grid, train, "predict_surface", "grid")
    grid <- .align_to_fit(grid, "grid")
    # A polygon grid -- create_grid_polygons() output, say -- was used as it
    # was.  st_nearest_feature() then gave each cell whichever covariate point
    # inside it the spatial index returned first, not the one at its centre,
    # while predict() pointized the cell by itself, so the covariates and the
    # location predicted at no longer matched: predictions off by up to 2.6 on
    # a 0-30 response, and a row shuffle of `covariates` moved them by up to
    # 4.7.  Reduce it to representative points first, as
    # area_of_applicability() does, so the surface is the POINT layer the
    # manual promises.  A cell whose centroid cannot be taken is refused in
    # this function's words, naming `grid` (.pointize_auto()).
    if (!grid_points)
      grid <- .xy_if_mixed(.pointize_auto(grid, "predict_surface", "`grid`"))
  }
  res <- attr(grid, "cell_size")
  # The joins below read the grid's points in x and y; a grid returned with
  # its rows as given can still mix dimensions or carry M.
  grid_xy <- function() .xy_if_mixed_or_m(sf::st_geometry(grid))

  # ---- clip ----------------------------------------------------------------
  if (!is.null(boundary)) {
    # Read only for its x and y: a boundary that mixes XY and XYZ (or XYM)
    # geometries, or whose geometries carry M, is read in XY before it is
    # moved (sf::st_transform() stopped with "OGR error" on an XYZ first
    # row) and before GEOS reads it.
    bnd <- .replay_crs_assumption(.xy_if_mixed_or_m(sf::st_geometry(boundary)),
                                  train, "predict_surface", "boundary")
    # A boundary still without a CRS is aligned to the fit's as other
    # functions align one: an R warning naming this function and argument.
    # ensure_projected()'s stamp was a log line naming neither.
    bnd <- if (is.null(.crs_or_null(target_crs)))
      ensure_projected(bnd)
    else .transform_or_stamp(bnd, target_crs, what = "boundary",
                             caller = "predict_surface")
    keep <- lengths(sf::st_intersects(grid_xy(), bnd)) > 0L
    grid <- grid[keep, , drop = FALSE]
    if (nrow(grid) == 0L)
      stop("predict_surface(): no grid points fall inside `boundary`.",
           call. = FALSE)
  }

  # ---- covariates ----------------------------------------------------------
  preds <- object$predictor_vars
  missing_preds <- setdiff(preds, names(grid))
  if (length(missing_preds) > 0L) {
    if (is.null(covariates))
      stop("predict_surface(): the model uses predictor(s) ",
           paste(sQuote(missing_preds), collapse = ", "),
           " which are absent from the grid. Supply `covariates`, or a `grid` ",
           "that already carries them.", call. = FALSE)
    if (!inherits(covariates, "sf"))
      stop("predict_surface(): `covariates` must be an sf object.", call. = FALSE)
    # st_nearest_feature() against an empty layer fails inside sf with a
    # message that names neither argument.
    if (nrow(covariates) == 0L)
      stop("predict_surface(): `covariates` has no rows, so there is no ",
           "nearest feature to take predictor values from.", call. = FALSE)

    cov_missing <- setdiff(missing_preds, names(covariates))
    if (length(cov_missing) > 0L)
      stop("predict_surface(): `covariates` lacks column(s) ",
           paste(sQuote(cov_missing), collapse = ", "), ".", call. = FALSE)

    # Read only for its x and y and its columns: a layer that mixes XY and
    # XYZ (or XYM) geometries, or whose geometries carry M, is read in XY
    # before it is moved, so ensure_projected() has no copy of it to warn
    # about, and before GEOS reads it.
    # A row of `covariates` without a location (a coordinate that is not
    # finite or lies beyond .POOL_XY_MAX) is emptied as given, before the
    # move, which can bring it back finite, where it was taken as the nearest
    # feature of the grid points around it.  An empty row is never the
    # nearest feature.
    covariates <- .replay_crs_assumption(
      .empty_unlocated(.xy_if_mixed_or_m(covariates)), train,
      "predict_surface", "covariates")
    covariates <- .align_to_fit(covariates, "covariates")
    nn  <- .surface_nearest_rows(grid_xy(), covariates)
    cdf <- sf::st_drop_geometry(covariates)[nn, missing_preds, drop = FALSE]
    for (cn in missing_preds) grid[[cn]] <- cdf[[cn]]
  }

  # ---- chunked prediction --------------------------------------------------
  n <- nrow(grid)
  # Inf is the natural way to ask for "one chunk, do not split", and it is
  # exactly what as.integer() turns into NA -- after which seq(by = NA) fails
  # with "invalid '(to - from)/by'", naming nothing the caller passed.
  .check_scalar(chunk_size, "chunk_size", "predict_surface", min = 1,
                max = .Machine$integer.max, what = "a single positive number")
  chunk_size <- max(1L, as.integer(chunk_size))
  starts <- seq(1L, n, by = chunk_size)

  preds_vec <- rep(NA_real_, n)
  se_vec    <- if (isTRUE(se)) rep(NA_real_, n) else NULL
  se_ok     <- isTRUE(se)

  # With se = TRUE the point predictions can be read off the same draw matrix
  # instead of drawing the posterior a second time per chunk -- but only when
  # the caller left the summary at its default, since predict() returns
  # colMeans(draws) for summary = "mean" and the column medians otherwise.
  # Deriving unconditionally would silently change .pred for anyone passing
  # summary = "median".
  #
  # Test a PREFIX match, not `.dots$summary`.  `$` on a list partial-matches
  # only in the other direction, so list(summ = "median")$summary is NULL --
  # but R's argument matching does partial-match a supplied name onto a
  # formal, so predict_surface(fit, se = TRUE, summ = "median") is legal and
  # arrives at predict() as summary = "median".  The old test missed it and
  # .pred silently became colMeans(d): exactly what the paragraph above says
  # must not happen.
  .dots <- list(...)
  dot_names <- names(.dots)
  supplied_summary <- !is.null(dot_names) &&
    any(nzchar(dot_names) & startsWith("summary", dot_names))
  derive_pred <- se_ok && !supplied_summary

  # A GWR's predict() counts, per call, the locations it cannot estimate and
  # those beyond the reach of the fit, and warns with that call's counts.
  # Predicting in chunks would pass on one warning per chunk; collected here,
  # each class is raised once with the totals over the grid.  Other backends
  # signal neither class, so for them this only evaluates the loop.
  .gwr_tally_predict_warnings({
    for (s in starts) {
      e   <- min(s + chunk_size - 1L, n)
      idx <- s:e
      blk <- grid[idx, , drop = FALSE]

      got_pred <- FALSE
      if (se_ok) {
        d <- try(stats::predict(object, newdata = blk, draws = TRUE, ...),
                 silent = TRUE)
        if (inherits(d, "try-error") || !is.matrix(d)) {
          se_ok <- FALSE          # backend has no draw interface; report once
        } else {
          se_vec[idx] <- apply(d, 2L, stats::sd)
          if (derive_pred) {
            preds_vec[idx] <- colMeans(d)
            got_pred <- TRUE
          }
        }
      }

      if (!got_pred) {
        p <- try(stats::predict(object, newdata = blk, ...), silent = TRUE)
        if (inherits(p, "try-error"))
          stop("predict_surface(): prediction failed on rows ", s, "-", e, ": ",
               as.character(p), call. = FALSE)
        # One value per row, or the assignment below recycles or truncates
        # whatever came back into .pred without a word that means anything.
        if (length(p) != length(idx))
          stop(sprintf(paste0("predict_surface(): predict() returned %d value(s) ",
                              "for the %d rows %d-%d; expected one per row. Check ",
                              "what the backend's predict() returns for the ",
                              "arguments passed through `...`."),
                       length(p), length(idx), s, e), call. = FALSE)
        preds_vec[idx] <- as.numeric(p)
      }
    }
  }, caller = "predict_surface", n_total = n)

  if (isTRUE(se) && !se_ok)
    .log_warn(paste0("predict_surface(): `se = TRUE` but this backend does not ",
                     "expose posterior draws; returning predictions only."))

  grid$.pred <- preds_vec
  # A grid that is an earlier surface carries that model's .pred_se.  It was
  # kept whenever this call did not replace it -- beside the new .pred, and
  # even after the log said "returning predictions only".
  grid$.pred_se <- if (isTRUE(se) && se_ok) se_vec else NULL

  attr(grid, "cell_size") <- res
  grid
}
