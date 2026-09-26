# -----------------------------------------------------------------------------
# CRS handling shared by the tessellation builders
# -----------------------------------------------------------------------------

#' Give a CRS-less points layer or boundary the other one's CRS
#'
#' When exactly one of the two has a CRS, the other is interpreted from it
#' the way harmonize_crs() does, warning either way: coordinates that look
#' like lon/lat are taken as EPSG:4326 and reprojected, anything else is
#' stamped.  Before this, the builders handed sf two layers in different CRSs
#' and every method stopped with sf's bare "st_crs(x) == st_crs(y) is not
#' TRUE": UTM points read from a CSV with a UTM boundary, or projected points
#' with a boundary read from a file that lost its .prj.  A CRS-less points
#' layer that does not look like lon/lat cannot be put in a GEOGRAPHIC
#' boundary's CRS -- stamping degrees onto UTM numbers is wrong -- so that is
#' refused with a message that says what to do.
#'
#' @param points,boundary sf layers; \code{boundary} may be \code{NULL}.
#' @param caller Function name for the messages.
#' @return A list with \code{points} and \code{boundary}.
#' @keywords internal
#' @noRd
.resolve_crsless_pair <- function(points, boundary, caller) {
  if (is.null(boundary)) return(list(points = points, boundary = boundary))
  pcrs <- sf::st_crs(points)
  bcrs <- sf::st_crs(boundary)
  if (is.na(bcrs) && !is.na(pcrs)) {
    boundary <- .transform_or_stamp(boundary, pcrs, "boundary", caller)
  } else if (is.na(pcrs) && !is.na(bcrs)) {
    if (isTRUE(sf::st_is_longlat(bcrs)) && !isTRUE(.looks_like_lonlat(points)$lonlat))
      stop(sprintf(paste0(
        "%s(): `points_sf` has no CRS and its coordinates do not look like ",
        "lon/lat, but `boundary` is in a geographic CRS (%s), so the two ",
        "cannot be placed in one space. Set the CRS of `points_sf` with ",
        "sf::st_crs(), or pass a boundary in the points' own projected CRS."),
        caller, .fold_crs_label(bcrs)), call. = FALSE)
    points <- .transform_or_stamp(points, bcrs, "points_sf", caller)
    attr(points, "crs_assumed") <- NULL
  }
  list(points = points, boundary = boundary)
}


#' Is `crs` a geographic (lon/lat) CRS?
#' @keywords internal
#' @noRd
.is_geographic_crs <- function(crs) {
  cc <- tryCatch(sf::st_crs(crs), error = function(e) sf::NA_crs_)
  !is.na(cc) && isTRUE(sf::st_is_longlat(cc))
}


#' Return a layer built in a projected CRS in the geographic CRS asked for
#'
#' A geographic \code{crs} names the CRS a tessellation is returned in; the
#' cells are built in metres (see build_tessellation()).  Long edges are
#' densified first, to a hundredth of the layer's extent, so a straight edge
#' in the working projection keeps its course in lon/lat instead of being
#' read as a great circle between its two ends.
#'
#' @param x sf layer (or \code{NULL}).
#' @param crs_out The geographic \code{sf::crs}, or \code{NULL} for no change.
#' @keywords internal
#' @noRd
.to_output_crs <- function(x, crs_out) {
  if (is.null(x) || is.null(crs_out)) return(x)
  bb  <- sf::st_bbox(x)
  ext <- max(as.numeric(bb["xmax"] - bb["xmin"]), as.numeric(bb["ymax"] - bb["ymin"]))
  if (is.finite(ext) && ext > 0) x <- sf::st_segmentize(x, dfMaxLength = ext / 100)
  sf::st_transform(x, crs_out)
}


#' Project a lon/lat boundary for a grid of equal-area cells
#'
#' The CRS ensure_projected() picks for distances, unless that CRS distorts
#' areas across the boundary by more than \code{.area_error_tol}, in which
#' case the equal-area choice (\code{purpose = "area"}) is used.  A UTM zone
#' on a local extent is within a quarter of a percent and is kept, so local
#' grids are unchanged; Web Mercator over near-global extents (whole cells
#' differing five-fold in area) and a zone stretched past its width are not.
#' Projected input is returned as it is.  Used by create_grid_polygons() and
#' create_grid_polygons_cached(), so the two lay a grid in the same CRS.
#'
#' @param boundary sf polygon layer.
#' @param caller Function name for the log line.
#' @keywords internal
#' @noRd
.project_for_grid <- function(boundary, caller = "create_grid_polygons") {
  crs0 <- sf::st_crs(boundary)
  proj <- ensure_projected(boundary)
  lonlat <- .is_longlat(boundary) ||
    (is.na(crs0) && identical(attr(proj, "crs_assumed"), "EPSG:4326"))
  if (!lonlat) return(proj)
  err <- .crs_area_error(proj, grid = TRUE)
  if (!is.finite(err) || err <= .area_error_tol) return(proj)
  src <- if (is.na(crs0)) sf::st_set_crs(boundary, 4326) else boundary
  eq  <- ensure_projected(src, purpose = "area")
  .log_warn(paste0("%s(): %s distorts areas across this boundary by up to %.1f%%, ",
                   "so its cells would not be equal-area; laying the grid in %s ",
                   "(ensure_projected(purpose = \"area\")) instead. Pass `crs` ",
                   "to choose the CRS yourself."),
            caller, .fold_crs_label(proj), 100 * err, .fold_crs_label(eq))
  if (is.na(crs0)) attr(eq, "crs_assumed") <- "EPSG:4326"
  eq
}

# -----------------------------------------------------------------------------
# Clip Target
# -----------------------------------------------------------------------------

#' Build a polygonal clip target from points and/or a boundary
#'
#' Resolves a single polygon to tessellate within. With a `boundary` it is
#' that boundary (optionally buffered by `expand`); without one it is the
#' axis-aligned bounding box of `points_sf` (the rectangle in the working
#' CRS), again optionally buffered, or a small buffer around the points when
#' they all share one x or one y. Reach for it to build the `boundary` that
#' `method = "hex"` and `"square"` require, to check that a study-area
#' polygon actually contains the observations before tessellating, or to
#' pass the same envelope to [create_voronoi_polygons()] and
#' [create_grid_polygons()] so that two tessellations of one dataset cover
#' identical ground.
#'
#' It is not the target [build_tessellation()] derives on its own:
#' `method = "voronoi"` without a boundary clips to the convex hull of the
#' points buffered by 2 percent of its diagonal, and there `expand` is always
#' a distance. A bounding box over a non-rectangular point cloud includes
#' corners with no data, so a hex or square grid laid over it has cells that
#' hold no points; pass the study-area polygon when there is one.
#'
#' @param points_sf An sf object with POINT/MULTIPOINT geometry.
#' @param boundary Optional polygonal sf object. One with no CRS, given with
#'   points that have one, is interpreted in the points' CRS as
#'   [harmonize_crs()] does, with a warning: coordinates that look like
#'   lon/lat are taken as EPSG:4326 and reprojected, anything else is stamped.
#' @param expand Numeric expansion distance or fraction (0–1 = fraction of
#'   extent). Absolute values are expressed in the units of the CRS the clip
#'   target is built in. Because [sf::st_buffer()] interprets `dist` as
#'   **metres** for lon/lat data while the fraction-of-extent form is derived
#'   from a bounding box measured in **degrees**, lon/lat input is projected to
#'   a local projected CRS first (see `@return`), so both forms agree.
#' @param quiet Logical; suppress this function's progress \code{message()}s.
#'   It does not silence R warnings, nor the package's console log echo
#'   (see \code{\link{spatialkit_quiet}} for that). Default \code{FALSE}.
#' @return An sf polygon layer representing the clip target. For lon/lat input
#'   the layer is returned in the automatically selected local projected CRS,
#'   not the input CRS; a message reports this unless `quiet = TRUE`.
#' @family spatial data preparation
#' @examples
#' library(sf)
#' set.seed(1)
#' pts <- st_as_sf(
#'   data.frame(x = 5e5 + runif(30, 0, 100), y = 5e6 + runif(30, 0, 100)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' # No boundary: the bounding box, expanded by 10% of the extent
#' box <- clip_target_for(pts, expand = 0.1, quiet = TRUE)
#' st_area(box)
#' @export
clip_target_for <- function(points_sf, boundary = NULL, expand = 0, quiet = FALSE) {
  .msg <- function(...) if (!quiet) message(...)
  .assert_sf(points_sf, c("POINT", "MULTIPOINT"), "points_sf")
  # .expand_distance() below TOLERATES a malformed `expand` by returning 0,
  # which turned `expand = c(0.05, 0.05)` into a silent no-op -- the returned
  # bbox was byte-identical to expand = 0, with no condition raised -- and a
  # degenerate bbox reached `if (expand == 0)` in the caller before that
  # tolerance applied, aborting on NA or on a length-2 vector.  Check once,
  # here, so both paths get the same answer.
  .check_scalar(expand, "expand", "clip_target_for", min = 0,
                what = "a single non-negative number (a fraction of the extent below 1, otherwise a distance)")

  .expand_distance <- function(ref_geom, expand) {
    if (!is.numeric(expand) || length(expand) != 1 || is.na(expand) || expand == 0) return(0)
    bb <- sf::st_bbox(ref_geom)
    dx <- as.numeric(bb$xmax - bb$xmin)
    dy <- as.numeric(bb$ymax - bb$ymin)
    if (expand > 0 && expand < 1) return(max(dx, dy) * expand)
    expand
  }

  # The fraction-of-extent form of `expand` is derived from the bounding box,
  # which for lon/lat data is measured in DEGREES -- but with s2 enabled
  # st_buffer() reads `dist` as METRES, so the two disagree by five orders of
  # magnitude.  Project to a local projected CRS first so the distance that is
  # computed and the distance that is buffered share the same units.  The
  # boundary is aligned to the projected points immediately below.
  if (.is_longlat(points_sf)) {
    .msg("clip_target_for(): input is lon/lat; projecting to a local projected ",
         "CRS so `expand` is measured in projected units. The returned clip ",
         "target uses that CRS, not the input CRS.")
    points_sf <- ensure_projected(points_sf)
  }

  crs_pts <- sf::st_crs(points_sf)

  if (!is.null(boundary)) {
    # .align_crs() leaves a CRS-less boundary as it is, so the target came
    # back with no CRS -- and for lon/lat points, in degrees, with `expand =
    # 20` buffering by 20 degrees rather than 20 metres.  Interpret it in the
    # points' (projected) CRS instead, as build_tessellation() does.
    if (is.na(sf::st_crs(boundary)) && !is.na(crs_pts))
      boundary <- .transform_or_stamp(boundary, crs_pts, "boundary", "clip_target_for")
    boundary <- .align_crs(boundary, points_sf)
    if (!any(sf::st_geometry_type(boundary) %in% c("POLYGON", "MULTIPOLYGON")))
      stop("clip_target_for(): `boundary` must be polygonal.")
    boundary <- .safe_make_valid(boundary)
    dist <- .expand_distance(boundary, expand)
    tgt <- if (dist != 0) suppressWarnings(sf::st_buffer(boundary, dist = dist)) else boundary
    return(sf::st_sf(geometry = sf::st_geometry(tgt)))
  }

  pts_geom <- sf::st_geometry(points_sf)
  if (length(pts_geom) == 0) stop("clip_target_for(): `points_sf` is empty.")

  bb <- sf::st_bbox(pts_geom)
  zero_w <- isTRUE(all.equal(as.numeric(bb$xmin), as.numeric(bb$xmax)))
  zero_h <- isTRUE(all.equal(as.numeric(bb$ymin), as.numeric(bb$ymax)))

  if (zero_w || zero_h) {
    .msg("clip_target_for(): degenerate bbox; using small buffer around points.")
    dist_use <- .expand_distance(points_sf, if (expand == 0) 0.05 else expand)
    if (dist_use <= 0) dist_use <- 1
    tgt <- suppressWarnings(sf::st_buffer(sf::st_union(pts_geom), dist = dist_use))
    tgt <- .safe_make_valid(tgt)
    sf::st_crs(tgt) <- crs_pts
    return(sf::st_sf(geometry = sf::st_geometry(tgt)))
  }

  tgt <- sf::st_as_sfc(bb, crs = crs_pts)
  dist <- .expand_distance(tgt, expand)
  if (dist != 0) tgt <- suppressWarnings(sf::st_buffer(tgt, dist = dist))
  tgt <- .safe_make_valid(tgt)
  sf::st_sf(geometry = tgt)
}

# -----------------------------------------------------------------------------
# Point-to-cell spatial index
# -----------------------------------------------------------------------------

#' Assign each point to its containing cell, returning cell_id values
#'
#' Uses st_intersects then falls back to st_nearest_feature for unmatched
#' points. The returned vector uses the cell_id column of `cells_sf`, not
#' raw row indices.
#'
#' @param pts An sf POINT object.
#' @param cells_sf An sf polygon object with a `cell_id` column.
#' @return Integer vector of cell_id values, one per row of `pts`, with
#'   \code{NA} for points that fall outside every cell, carrying the
#'   nearest-cell repairs as \code{attr(, "snapped")} (\code{n}, \code{which},
#'   \code{distance}); the callers move that record into \code{params}.
#' @keywords internal
#' @noRd
.build_point_cell_index <- function(pts, cells_sf) {
  n <- nrow(pts)
  if (n == 0L || nrow(cells_sf) == 0L)
    return(structure(rep(NA_integer_, n),
                     snapped = list(n = 0L, which = integer(0), distance = numeric(0))))

  cell_ids <- cells_sf$cell_id
  hits <- sf::st_intersects(pts, cells_sf)
  
  index <- vapply(hits, function(row_idx) {
    row_idx <- as.integer(row_idx)
    if (length(row_idx) == 0L)     return(NA_integer_)
    if (length(row_idx) == 1L)     return(cell_ids[row_idx])
    # Ambiguous: pick the cell with the smallest cell_id
    min(cell_ids[row_idx])
  }, integer(1))

  # Points that fall in no cell.  For a Voronoi tessellation, or a grid clipped
  # to a boundary, "outside every cell" means outside the study area, and
  # snapping such a point to the NEAREST cell silently counts it as if it were
  # inside: a summary built from this index counted all 40 points of a layer
  # whose study area held 10, while assign_features_to_polygons() on the same
  # cells correctly reported 30 misses.  Nearest-cell assignment is still the
  # right repair for a point a hair outside a cell edge (floating-point, a
  # clipped sliver), so keep it -- but bound it by the cell size and say how
  # many points were snapped, and leave the genuinely-outside ones NA.
  unmatched <- which(is.na(index))
  snapped   <- list(n = 0L, which = integer(0), distance = numeric(0))
  if (length(unmatched) > 0L) {
    near  <- suppressWarnings(sf::st_nearest_feature(pts[unmatched, ], cells_sf))
    valid <- is.finite(near)
    if (any(valid)) {
      d <- suppressWarnings(as.numeric(sf::st_distance(
        sf::st_geometry(pts[unmatched[valid], ]),
        sf::st_geometry(cells_sf)[near[valid]], by_element = TRUE)))
      # A tolerance of one thousandth of the median cell width: wide enough for
      # any rounding at a shared edge, far too narrow to swallow a point that
      # is genuinely outside the study area.
      cw  <- suppressWarnings(stats::median(vapply(
        sf::st_geometry(cells_sf), function(g) {
          bb <- sf::st_bbox(g)
          as.numeric(bb[["xmax"]] - bb[["xmin"]])
        }, numeric(1)), na.rm = TRUE))
      tol <- if (is.finite(cw) && cw > 0) cw / 1000 else 0
      snap <- valid
      snap[valid] <- is.finite(d) & d <= tol
      index[unmatched[snap]] <- cell_ids[near[snap]]
      # Which points were repaired, and how far outside each sat: the comment
      # above promised to say how many, and nothing did.
      d_all <- rep(NA_real_, length(unmatched)); d_all[valid] <- d
      snapped <- list(n = sum(snap), which = as.integer(unmatched[snap]),
                      distance = as.numeric(d_all[snap]))
      if (snapped$n > 0L)
        .log_info(paste0("build_tessellation(): %d point(s) sitting just outside ",
                         "every cell (at most %.3g units, within a thousandth of ",
                         "the median cell width) were assigned to the nearest cell."),
                  snapped$n, max(snapped$distance))
    }
  }
  n_out <- sum(is.na(index))
  if (n_out > 0L)
    .log_info(paste0("build_tessellation(): %d of %d point(s) fall outside every ",
                     "cell and are recorded as NA in `index`."),
              n_out, length(index))

  attr(index, "snapped") <- snapped
  index
}

# -----------------------------------------------------------------------------
# Voronoi Tessellation
# -----------------------------------------------------------------------------

#' Create Voronoi polygons from points with CRS handling and optional clipping
#'
#' Assigns every location in the study area to its nearest input point, giving
#' one cell per point. This is the tessellation to reach for when the
#' observations themselves define the regions of interest (sampling sites,
#' monitoring stations, service points), because cell size then adapts to
#' sampling density instead of being imposed by a fixed grid: dense areas get
#' small cells and sparse areas large ones. Prefer [create_grid_polygons()]
#' instead when you need equal-area cells or a resolution independent of where
#' the data happen to be.
#'
#' The heavy lifting is [sf::st_voronoi()].  What this adds is the surrounding
#' bookkeeping: projecting lon/lat input, building and buffering an envelope
#' so edge cells are bounded, clipping to `boundary`, restoring the
#' point-to-cell correspondence that `st_voronoi()` scrambles, and stamping
#' stable `cell_id` values.
#'
#' The generators are the points' vertices, not the features. A MULTIPOINT
#' feature with several vertices therefore gets one cell per vertex, and its
#' `index` entry is the smallest `cell_id` among the cells it touches; the
#' others are referenced by no feature. Other functions in the package
#' ([prep_model_data()], [make_folds()]) reduce such a feature to its
#' centroid instead, so cast to POINT, or take centroids, first if one cell
#' per feature is what you want.
#'
#' @param points_sf An sf object with POINT/MULTIPOINT geometries.
#' @param boundary Optional polygonal sf object. When exactly one of
#'   `points_sf` and `boundary` has a CRS, the other is interpreted in it as
#'   [harmonize_crs()] does, with a warning (lon/lat-looking coordinates are
#'   reprojected from EPSG:4326, others are stamped); CRS-less points that do
#'   not look like lon/lat cannot take a geographic boundary's CRS, and are
#'   refused with an error.
#' @param expand Numeric; absolute distance, in the working CRS's units, by
#'   which the boundary (or the hull derived from the points) is grown before
#'   the diagram is built. With `clip = TRUE` the cells are clipped to the
#'   grown boundary, so they reach `expand` beyond the study area, and a point
#'   up to `expand` outside it gets a cell and an `index` value. The grown
#'   boundary is the one returned as `boundary`.
#' @param clip Logical; intersect cells with boundary.
#' @param keep_duplicates Logical. Has no effect on the result: coincident
#'   points are merged before the diagram is built either way, so they share
#'   one cell and all of them are indexed to it.
#' @param crs Optional target CRS. A projected CRS is the working CRS. A
#'   geographic one (EPSG:4326, say) is the CRS the result is returned in:
#'   the cells are built in the local projected CRS [ensure_projected()]
#'   picks for the points, so they are nearest-point cells on the ground, and
#'   are then transformed, with long edges densified.
#' @param quiet Logical; suppress this function's progress \code{message()}s.
#'   It does not silence R warnings, nor the package's console log echo
#'   (see \code{\link{spatialkit_quiet}} for that). Default \code{FALSE}.
#' @return A list with \code{cells}, \code{index}, \code{boundary},
#'   \code{method} and \code{params}.  \code{index} holds one \code{cell_id}
#'   per row of \code{points_sf}, and \code{NA} for a point that falls outside
#'   every cell, which means outside the study area (grown by \code{expand}
#'   when it is positive), so a summary built from it counts only the points
#'   the tessellation actually covers.  \code{boundary} is the boundary the
#'   cells were built in: the one supplied or derived, grown by
#'   \code{expand}.
#' @family tessellation
#' @examples
#' library(sf)
#' set.seed(1)
#' pts <- st_as_sf(
#'   data.frame(x = 5e5 + runif(15, 0, 100), y = 5e6 + runif(15, 0, 100)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' res <- create_voronoi_polygons(pts, quiet = TRUE)
#' res$cells   # one polygon per unique point, with stable cell_id
#' res$index   # cell_id assignment for each input point
#' @export
create_voronoi_polygons <- function(
    points_sf, boundary = NULL, expand = 0, clip = TRUE,
    keep_duplicates = FALSE, crs = NULL, quiet = FALSE
) {
  .assert_sf(points_sf, c("POINT", "MULTIPOINT"), "points_sf")
  if (nrow(points_sf) < 1) stop("create_voronoi_polygons(): `points_sf` has no rows.")
  .msg <- function(...) if (!quiet) message(...)

  pts <- points_sf
  crs_out <- NULL
  if (!is.null(crs)) {
    pts <- .transform_or_stamp(pts, crs, "points_sf", "create_voronoi_polygons")
    if (!is.null(boundary))
      boundary <- .transform_or_stamp(boundary, crs, "boundary", "create_voronoi_polygons")
    # A geographic `crs` is where the cells are RETURNED, not where they are
    # built: st_voronoi() on degrees is not a nearest-point partition (at
    # 55N, 21% of sampled locations sat in another point's cell) and s2 read
    # the degree-sized hull buffer below as metres.
    if (.is_geographic_crs(crs)) {
      crs_out <- sf::st_crs(pts)
      pts <- ensure_projected(pts)
      if (!is.null(boundary)) boundary <- .align_crs(boundary, pts)
    }
  } else {
    if (.is_longlat(pts)) pts <- ensure_projected(pts)
    if (!is.null(boundary)) {
      # One side with no CRS takes the other's (.align_crs() leaves it as it
      # is, and sf then refused the pair with its bare CRS-mismatch error).
      pair <- .resolve_crsless_pair(pts, boundary, "create_voronoi_polygons")
      pts <- pair$points; boundary <- pair$boundary
      # CRS-less lon/lat points that just took a geographic boundary's CRS.
      if (.is_longlat(pts)) pts <- ensure_projected(pts)
      boundary <- .align_crs(boundary, pts)
    }
  }

  if (is.null(boundary)) {
    .msg("create_voronoi_polygons(): deriving boundary via convex hull of points.")
    hull <- sf::st_convex_hull(sf::st_union(sf::st_geometry(pts)))
    # Buffer the hull slightly so edge Voronoi cells are not clipped to
    # degenerate slivers along the convex hull boundary.
    bb   <- sf::st_bbox(hull)
    diag <- sqrt((bb$xmax - bb$xmin)^2 + (bb$ymax - bb$ymin)^2)
    hull <- suppressWarnings(sf::st_buffer(hull, dist = 0.02 * max(diag, 1)))
    hull_sfc <- sf::st_sfc(hull, crs = sf::st_crs(pts))
    boundary <- sf::st_sf(geometry = hull_sfc)
  } else {
    if (!any(sf::st_geometry_type(boundary) %in% c("POLYGON", "MULTIPOLYGON")))
      stop("create_voronoi_polygons(): `boundary` must be polygonal.")
  }
  boundary <- .safe_make_valid(boundary)

  # `isTRUE(is.numeric(expand)) && expand > 0` short-circuits to FALSE for a
  # character `expand`, dropping the requested expansion with no condition at
  # all, and raises R's own "'length = 2' in coercion to 'logical(1)'" for a
  # length-2 one -- an error for input that build_tessellation(method = "hex")
  # accepts without complaint.
  .check_scalar(expand, "expand", "create_voronoi_polygons", min = 0,
                what = "a single non-negative buffer distance")

  boundary_expanded <- if (expand > 0) {
    suppressWarnings(sf::st_buffer(boundary, dist = expand))
  } else boundary

  geom_pts <- sf::st_geometry(pts)
  pts_for_graph <- if (isTRUE(keep_duplicates)) geom_pts else .dedup_points(geom_pts)
  mp  <- sf::st_union(pts_for_graph)
  env <- sf::st_union(sf::st_geometry(boundary_expanded))
  vor <- suppressWarnings(sf::st_voronoi(mp, envelope = env))
  cells <- sf::st_collection_extract(vor, "POLYGON", warn = FALSE)
  cells_sfc <- sf::st_sfc(cells, crs = sf::st_crs(pts))
  cells <- sf::st_sf(geometry = cells_sfc)
  cells <- .safe_make_valid(cells)

  if (clip) {
    clip_to <- if (expand > 0) boundary_expanded else boundary
    # Union the boundary first (as the Voronoi envelope above already does).
    # Intersecting against a multi-feature boundary splits every straddling
    # cell into one row per boundary feature and grafts the boundary's
    # attribute columns onto the result, breaking the documented
    # "one polygon per unique point" contract.
    clip_to <- sf::st_union(sf::st_geometry(clip_to))
    cells <- suppressWarnings(sf::st_intersection(cells, clip_to))
    # st_intersection can produce non-polygon slivers (LINESTRING, POINT,
    # GEOMETRYCOLLECTION); keep only POLYGON/MULTIPOLYGON and non-empty rows
    gtypes <- as.character(sf::st_geometry_type(cells, by_geometry = TRUE))
    keep <- gtypes %in% c("POLYGON", "MULTIPOLYGON") & !sf::st_is_empty(cells)
    cells <- cells[keep, , drop = FALSE]
  }

  if (nrow(cells) > 0) {
    cells <- ensure_stable_poly_id(cells, id_col = "cell_id")
  } else {
    cells$cell_id <- integer(0)
  }

  # Build point → cell index using cell_id values.  The snapping record rides
  # on `params`, not on `index`, so `index` stays a plain integer vector.
  index   <- .build_point_cell_index(pts, cells)
  snapped <- attr(index, "snapped")
  attr(index, "snapped") <- NULL

  list(
    cells    = .to_output_crs(cells, crs_out),
    index    = index,
    # The boundary the cells were built in and clipped to.  With expand > 0
    # that is the grown one: returning the ungrown input described a study
    # area the cells overhang by `expand` (1.93 km^2 of cells against a
    # returned 1 km^2) and in which indexed points lay outside.
    boundary = .to_output_crs(boundary_expanded, crs_out),
    method   = "voronoi",
    params   = list(clip = clip, expand = expand, keep_duplicates = keep_duplicates,
                    snapped = snapped)
  )
}

# -----------------------------------------------------------------------------
# Grid Tessellations (Hex / Square)
# -----------------------------------------------------------------------------

#' Create square or hexagonal grid polygons over a boundary
#'
#' Lays a regular grid of equal-area cells over `boundary` and clips it to that
#' boundary. Reach for this rather than [create_voronoi_polygons()] when cell
#' size should be a decision you make, instead of one dictated by where the
#' observations happen to be. That is the case when you need per-cell rates
#' comparable across the map, or a resolution that stays fixed as the sample
#' grows. Hexagons (`type = "hex"`) avoid the axis-aligned artefacts of squares
#' and give every cell the same distance to all six neighbours, which matters
#' for anything that reads neighbourhoods.
#'
#' Size the grid with exactly one of `target_cells` (roughly how many cells you
#' want, the package derives the rest), `cellsize` (a fixed edge length in CRS
#' units) or `n` (a fixed number of columns and rows). See `@param cellsize`
#' for what happens when more than one is given.
#'
#' @param boundary Polygonal sf or sfc object.
#' @param target_cells Optional approximate desired number of cells.  The cell
#'   \emph{size} is derived from it as \code{sqrt(area / target_cells)}, where
#'   \code{area} is that of the boundary's bounding box, so square grids get
#'   square cells; for hex grids the count is adjusted for hexagonal packing
#'   density and the size rounded so that a whole number of hexagons spans
#'   the longer side of the box.  Neither depends on which way the boundary
#'   lies.  The word "approximate" is load bearing: a
#'   grid of square cells over an elongated bounding box needs more of them
#'   than a grid of rectangles would (a 1000 x 1 strip at
#'   \code{target_cells = 9} yields cells of side 10.5 and about 95 of them),
#'   and clipping to an irregular boundary moves the count again.  Pass
#'   \code{cellsize} when the count matters more than the shape.
#' @param type Grid type: `"square"` (the default) or `"hex"`.
#' @param cellsize Optional numeric cell size (length 1 or 2), in the units of
#'   the working CRS. Takes precedence over `n`: if both are supplied,
#'   `cellsize` is used, `n` is ignored and a warning is logged. Supply exactly
#'   one of `target_cells`, `cellsize` and `n`.  For `type = "hex"` a hexagon
#'   is defined by a single edge-to-edge distance, so only `cellsize[1]` is
#'   used and a differing `cellsize[2]` is ignored with a logged warning.
#' @param n Optional grid resolution (integer, length 1 or 2) giving the number
#'   of columns and rows to divide the boundary's bounding box into; the cell
#'   size is derived from it. [sf::st_make_grid()] derives hexagon placement
#'   from `cellsize` alone, so for `type = "hex"` `n` does not set the number of
#'   cells, although it does change the grid, because the `cellsize` derived
#'   from it is what the hexagons are built with. Ignored (with a logged
#'   warning) when `cellsize` is also supplied. Passing both would otherwise
#'   truncate the grid to `n[1]` x `n[2]` cells anchored at the bounding-box
#'   corner, covering only part of the boundary.
#' @param clip Logical; clip grid to boundary.
#' @param crs Optional target CRS. When `NULL` (default) a lon/lat boundary is
#'   projected with [ensure_projected()], which changes the CRS of the returned
#'   grid; a message reports this unless `quiet = TRUE`. When that CRS would
#'   distort cell areas across the boundary by more than 1 percent (Web
#'   Mercator over a near-global extent, a UTM zone stretched well past its
#'   width), `ensure_projected(purpose = "area")` is used instead, with a
#'   logged warning, so the cells stay equal-area; a local extent keeps its
#'   UTM zone. A geographic `crs` (EPSG:4326, say) is the CRS the grid is
#'   returned in: a grid sized by `target_cells` or `n` is laid in that
#'   projected CRS and then transformed, with long edges densified. A
#'   `cellsize` is in the units of `crs`, so with a geographic `crs` it is in
#'   degrees and the grid is laid in degrees, as asked; such cells are not
#'   equal-area.
#' @param quiet Logical; suppress this function's progress \code{message()}s.
#'   It does not silence R warnings, nor the package's console log echo
#'   (see \code{\link{spatialkit_quiet}} for that). Default \code{FALSE}.
#' @param max_cells Upper bound on the number of cells the grid may have,
#'   estimated from the boundary's bounding box before anything is built.
#'   Default \code{1e6}. A \code{cellsize} in the wrong units (metres on a
#'   boundary in kilometres, say) asks for a grid that cannot be built, and
#'   this refuses it with a message before memory is exhausted. Set to
#'   \code{Inf} to disable.
#' @return An sf polygon layer with poly_id column.
#' @family tessellation
#' @examples
#' library(sf)
#' bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
#'   c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
#' ))), crs = 32632))
#' grid_sq  <- create_grid_polygons(bnd, target_cells = 100, type = "square")
#' grid_hex <- create_grid_polygons(bnd, target_cells = 100, type = "hex")
#' nrow(grid_sq)
#' # Hex counts run above target because clipping keeps every hexagon that
#' # merely overhangs the boundary; the inflation is proportionally larger
#' # at small target_cells.
#' nrow(grid_hex)
#' @export
create_grid_polygons <- function(
    boundary, target_cells = NULL, type = c("square", "hex"),
    cellsize = NULL, n = NULL, clip = TRUE, crs = NULL, quiet = FALSE,
    max_cells = 1e6
) {
  type <- match.arg(type)
  .msg <- function(...) if (!quiet) message(...)

  .as_sf <- function(x) {
    if (inherits(x, "sf")) return(x)
    if (inherits(x, "sfc")) return(sf::st_sf(geometry = x))
    stop("create_grid_polygons(): 'boundary' must be an sf or sfc object.")
  }

  boundary <- .as_sf(boundary)
  if (!all(as.character(sf::st_geometry_type(boundary, by_geometry = TRUE)) %in%
           c("POLYGON", "MULTIPOLYGON")))
    stop("create_grid_polygons(): 'boundary' must be polygonal (POLYGON/MULTIPOLYGON).")

  crs_out <- NULL
  if (!is.null(crs)) {
    boundary <- .transform_or_stamp(boundary, crs, "boundary", "create_grid_polygons")
    # A geographic `crs` is where the grid is RETURNED.  Laid in degrees the
    # cells were neither square nor equal-area, and clipping them under s2
    # often stopped with "Edge 0 is degenerate" or left points inside the
    # boundary with no cell.  An explicit `cellsize` is in the units of
    # `crs`, degrees, so a grid sized by it is still laid in degrees.
    if (.is_geographic_crs(crs) && is.null(cellsize)) {
      crs_out  <- sf::st_crs(boundary)
      boundary <- .project_for_grid(boundary)
    }
  } else {
    crs_before <- sf::st_crs(boundary)
    boundary <- .project_for_grid(boundary)
    if (!identical(crs_before, sf::st_crs(boundary)))
      .msg("create_grid_polygons(): projecting `boundary` to a local projected ",
           "CRS; the returned grid uses that CRS. Pass `crs` to control it.")
  }
  boundary <- .safe_make_valid(boundary)

  bb <- sf::st_bbox(boundary)
  w  <- as.numeric(bb["xmax"] - bb["xmin"])
  h  <- as.numeric(bb["ymax"] - bb["ymin"])
  if (!(is.finite(w) && is.finite(h) && w > 0 && h > 0))
    stop("create_grid_polygons(): boundary bbox has non-positive extent.")
  env <- sf::st_as_sfc(bb, crs = sf::st_crs(boundary))

  # Parse and validate `n` once, up front, so an invalid value is rejected the
  # same way whether or not `cellsize` was also supplied (it used to be
  # silently coerced to NULL in the cellsize branch and to error otherwise).
  if (!is.null(n)) {
    if (length(n) == 1L) n <- rep(n, 2L)
    n <- suppressWarnings(as.integer(n[1:2]))
    if (any(is.na(n)) || any(n < 1))
      stop("create_grid_polygons(): 'n' must be integer >= 1.", call. = FALSE)
  }

  # Derive n and/or cellsize
  #
  # `cellsize_supplied` records whether the CALLER fixed the cell size, as
  # opposed to the package deriving it from `n` or `target_cells`.  It decides
  # whether `n` is forwarded to st_make_grid() below; see the note there.
  cellsize_supplied <- !is.null(cellsize)
  if (!is.null(cellsize)) {
    if (length(cellsize) == 1L) cellsize <- rep(cellsize, 2L)
    if (length(cellsize) != 2L || any(!is.finite(cellsize)) || any(cellsize <= 0))
      stop("create_grid_polygons(): 'cellsize' must be positive numeric (length 1 or 2).")
    if (!is.null(n)) {
      .log_warn(paste0("create_grid_polygons(): both `cellsize` and `n` were ",
                       "supplied; `cellsize` wins and `n` (%s) is ignored. ",
                       "Pass one or the other -- `n` would cap the grid at ",
                       "%d x %d cells anchored at the bbox corner, leaving ",
                       "most of the boundary uncovered."),
                paste(n, collapse = " x "), n[1L], n[2L])
      n <- NULL
    }
    # Same rule for `target_cells`, which reaches here as build_tessellation()'s
    # `approx_n_cells`.  Without this, build_tessellation(method = "hex",
    # approx_n_cells = 25, cellsize = 10) returned however many cells a
    # 10-unit lattice holds and said nothing about the 25.
    if (!is.null(target_cells)) {
      .log_warn(paste0("create_grid_polygons(): both `cellsize` and ",
                       "`target_cells` were supplied; `cellsize` wins and ",
                       "`target_cells` (%s) is ignored. Pass one or the other."),
                format(target_cells))
      target_cells <- NULL
    }
  } else if (!is.null(n)) {
    cellsize <- c(w / n[1], h / n[2])
  } else {
    if (is.null(target_cells) || !is.finite(target_cells) || target_cells < 1)
      stop("create_grid_polygons(): supply either 'cellsize', 'n', or a positive 'target_cells'.")
    
    # For st_make_grid(), hex `cellsize` is the distance between opposite
    # edges, so a hexagon's area is (sqrt(3)/2) * cellsize^2 — smaller than
    # a square of the same cellsize.  A cellsize derived for N squares
    # therefore yields ~N / (sqrt(3)/2) hexagons.  Shrink the effective
    # target by that factor so the derived cellsize produces roughly
    # `target_cells` hexes (previously this divided instead of multiplying,
    # overshooting the requested count by ~33%).
    effective_target <- if (identical(type, "hex")) {
      target_cells * (sqrt(3) / 2)
    } else {
      target_cells
    }
    nx <- max(1L, round(sqrt(effective_target * (w / h))))
    ny <- max(1L, round(ceiling(effective_target / nx)))
    n  <- c(nx, ny)
    # cellsize = c(w/nx, h/ny) forces an exact tiling of the bbox, which means
    # the "square" grid is only square when w/nx happens to equal h/ny: on a
    # 2:1 boundary target_cells = 9 gave 50 x 33.3 cells, and on a 1000:1 strip
    # 10.5:1 rectangles and 95 cells for a requested 9.  Derive ONE edge from
    # the area per cell so the cells are square, as the type name and the
    # "equal-area cells" documentation both say, and let the grid overhang the
    # bbox (clip = TRUE trims it, and `n` is dropped below so st_make_grid()
    # covers the whole boundary).
    if (!identical(type, "hex")) {
      side <- sqrt((w * h) / effective_target)
      if (is.finite(side) && side > 0) {
        cellsize <- c(side, side)
        n <- NULL
      } else {
        cellsize <- c(w / nx, h / ny)
      }
    } else {
      # st_make_grid() builds hexagons from cellsize[1] alone, and that was
      # w / nx: the WIDTH of the box over a column count rounded, and floored
      # at 1, from the aspect ratio.  A tall narrow boundary therefore got
      # hexagons about as wide as the whole box -- a 1 x 1000 strip at target
      # 9 gave 1734 of them where the same strip lying flat gave 89.  Count
      # along the LONGER side instead.  For a boundary at least as wide as it
      # is tall that is exactly w / nx, so those grids are unchanged; a tall
      # one now gets the grid its lying-down twin gets.
      long <- max(w, h)
      side <- long / max(1L, round(sqrt(effective_target * long / min(w, h))))
      cellsize <- c(side, side)
    }
  }

  # st_make_grid(square = FALSE) uses only cellsize[1] for a hexagon -- the
  # distance between opposite edges -- so a length-2 `cellsize` had its second
  # component silently ignored while the max_cells estimate below used BOTH,
  # putting the estimate out by cellsize[2]/cellsize[1]: a grid refused at
  # max_cells = 200 via c(10, 10) was built (1,319 cells) via c(10, 100).
  # Collapse it here so the guard measures the grid that will actually be
  # built.
  if (identical(type, "hex") && length(cellsize) == 2L &&
      !isTRUE(all.equal(cellsize[1L], cellsize[2L]))) {
    if (cellsize_supplied)
      .log_warn(paste0("create_grid_polygons(): hexagonal cells are defined by a ",
                       "single edge-to-edge distance; using cellsize[1] = %s and ",
                       "ignoring cellsize[2] = %s."),
                format(cellsize[1L]), format(cellsize[2L]))
    cellsize[2L] <- cellsize[1L]
  }

  # Refuse an absurd grid before building it.  A `cellsize` in the wrong
  # units -- 10 on a boundary measured in metres that spans 100 km -- asks for
  # 1e8 polygons, which is ~80 GB and hours of st_make_grid() with nothing to
  # say until R itself fails to allocate.  Estimate the count from the bbox
  # first (hex cells are ~15% smaller, so the estimate is inflated by that).
  n_est <- ceiling(w / cellsize[1L]) * ceiling(h / cellsize[2L])
  if (identical(type, "hex")) n_est <- n_est / (sqrt(3) / 2)
  if (is.finite(max_cells) && n_est > max_cells)
    stop(sprintf(paste0("create_grid_polygons(): a cell size of %s x %s on a ",
                        "boundary of %s x %s would produce about %s cells, above ",
                        "`max_cells` = %s. Check that `cellsize` is in the ",
                        "boundary's CRS units (%s), or raise `max_cells` if the ",
                        "count is intended."),
                 format(cellsize[1L], digits = 4), format(cellsize[2L], digits = 4),
                 format(w, digits = 4), format(h, digits = 4),
                 format(n_est, big.mark = ",", scientific = FALSE, digits = 3),
                 format(max_cells, big.mark = ",", scientific = FALSE),
                 sf::st_crs(boundary)$units_gdal %||% "unknown"),
         call. = FALSE)

  grid_args <- list(x = env, what = "polygons",
                    square = identical(type, "square"))
  # REGRESSION NOTE -- do not "simplify" this back to passing both whenever
  # both are non-NULL.
  #
  # st_make_grid() does NOT ignore `n` when `cellsize` is given: for square
  # grids it uses `cellsize` for the cell dimensions AND `nx = n[1]`,
  # `ny = n[2]` for the counts, anchored at the bbox corner.  That is exactly
  # what we want when the PACKAGE derived `cellsize` from `n` (the `n` and
  # `target_cells` branches above): omitting `n` there makes sf recompute
  # nx = ceiling(w / cellsize[1]), which floating-point division pushes one
  # past the intended count (e.g. w = 100, n = 9 gives 100/(100/9) = 9.0000...4
  # -> 10 columns).
  #
  # It is exactly what we do NOT want when the CALLER supplied `cellsize`: an
  # unrelated `n` then truncates the grid to n[1] x n[2] cells in one corner
  # of the bbox and silently leaves the rest of the boundary uncovered
  # (cellsize = 25 with n = 2 on a 100x100 boundary covered 2500 of 10000 --
  # and clip = TRUE discards nothing, so it looks like an ordinary grid).
  # `n` is dropped with a warning in that branch above, so it is NULL here.
  #
  # For hex grids sf short-circuits to make_hex_grid() and reads `cellsize`
  # only, so the extra argument is inert there.
  if (!is.null(cellsize)) grid_args$cellsize <- cellsize
  if (!is.null(n) && !cellsize_supplied) grid_args$n <- n
  grid_sfc <- do.call(sf::st_make_grid, grid_args)
  if (length(grid_sfc) == 0L)
    stop("create_grid_polygons(): st_make_grid() produced zero cells.")

  grid_sfc <- sf::st_sfc(grid_sfc, crs = sf::st_crs(boundary))
  grid_sf <- sf::st_sf(poly_id = seq_along(grid_sfc), geometry = grid_sfc)

  if (isTRUE(clip)) {
    # Union first: a multi-feature boundary would otherwise split straddling
    # cells into one row per boundary feature and graft the boundary's
    # attribute columns onto the grid.
    clip_to <- sf::st_union(sf::st_geometry(boundary))
    grid_sf <- suppressWarnings(sf::st_intersection(.safe_make_valid(grid_sf), clip_to))
    # st_intersection can produce non-polygon slivers; keep only valid polygons
    gtypes <- as.character(sf::st_geometry_type(grid_sf, by_geometry = TRUE))
    keep <- gtypes %in% c("POLYGON", "MULTIPOLYGON") & !sf::st_is_empty(grid_sf)
    grid_sf <- grid_sf[keep, , drop = FALSE]
    grid_sf$poly_id <- seq_len(nrow(grid_sf))
  }
  .to_output_crs(grid_sf, crs_out)
}

# -----------------------------------------------------------------------------
# Unified Tessellation Builder
# -----------------------------------------------------------------------------

#' Build a tessellation (Voronoi, Delaunay triangles, hex grid, or square grid)
#'
#' The single entry point for turning a point pattern into analysis regions, and
#' the first step of the package's pipeline.  It wraps the four tessellation
#' methods behind one interface that handles CRS projection, clipping and stable
#' cell identifiers consistently, and returns the cell layer together with the
#' point-to-cell index that \code{\link{assign_features_to_polygons}()} and
#' \code{\link{summarize_by_cell}()} consume.  Prefer it to the individual
#' constructors whenever you might want to compare methods: the return shape
#' does not change with \code{method}, so swapping \code{"voronoi"} for
#' \code{"hex"} costs one argument.
#'
#' Which method to reach for.  \code{"voronoi"} gives one cell per point, so
#' resolution follows sampling density.  That is the choice when the
#' observations themselves define the regions.  \code{"hex"} and
#' \code{"square"} give equal-area cells on a fixed grid, so cell size is a
#' decision you make, and the data does not make it for you; hexagons avoid the
#' axis-aligned artefacts of squares and have uniform neighbour distances.
#' \code{"triangles"} returns the Delaunay triangulation, useful for
#' interpolation and adjacency work; it is not meant as an aggregation unit.
#' \code{\link{determine_optimal_levels}()} will suggest a cell count from the
#' spatial structure of the data.
#'
#' \code{"voronoi"} and \code{"triangles"} are built on the points' vertices:
#' a MULTIPOINT feature with several vertices gets one cell (or triangle
#' corner) per vertex, and its \code{index} entry is the smallest
#' \code{cell_id} among the cells it touches.  See
#' \code{\link{create_voronoi_polygons}()}.
#'
#' @param points_sf An sf object with POINT/MULTIPOINT geometry.
#' @param boundary Polygonal sf/sfc study area. **Required** for
#'   `method = "hex"` and `method = "square"`, which have no extent of their
#'   own to lay a grid over and error without it; supply the study-area polygon,
#'   or build one from the points with [clip_target_for()]. **Optional** for
#'   `method = "voronoi"` and `method = "triangles"`, which derive their extent
#'   from the points themselves and use `boundary` only to clip the result when
#'   `clip = TRUE`. When exactly one of `points_sf` and `boundary` has a CRS,
#'   the other is interpreted in it as [harmonize_crs()] does, with a warning;
#'   CRS-less points that do not look like lon/lat cannot take a geographic
#'   boundary's CRS and are refused with an error. When neither has one, both
#'   stay in the same unnamed planar space.
#' @param method One of "voronoi", "triangles", "hex", "square".
#' @param approx_n_cells Approximate number of cells.  Read by
#'   \code{method = "hex"} and \code{"square"} only: \code{"voronoi"} grows
#'   one cell per input point and \code{"triangles"} one triangle per
#'   neighbouring triple, so neither has a count to set, and both warn that the
#'   argument was ignored.  For a Voronoi cell count, place the seeds with
#'   \code{\link{get_voronoi_seeds}()} and tessellate those.  For hex
#'   grids the target is adjusted for packing density; the actual count after
#'   clipping to an irregular boundary may differ noticeably.  Besides a
#'   number, this accepts what the level-selection step returned: the integer
#'   vector of ranked candidates from \code{\link{determine_optimal_levels}()}
#'   (its first element is used), a \code{\link{select_resolution}()} result
#'   (its \code{$best}), or a \code{\link{resolution_profile}()} (read with
#'   \code{select_resolution()} at its default criterion).  The count used
#'   is returned as \code{params$approx_n_cells} and where it came from as
#'   \code{params$approx_n_cells_from} (\code{NULL} for a plain number).
#'   A count read off a profile or selection is a number of k-means cells:
#'   every one occupied, and small where the points are dense.  A lattice
#'   lays that many equal cells over the whole boundary, so on clustered
#'   points many of them hold no point (about half, on six clusters in a
#'   square); on evenly spread points it matches.  \code{params$cells_occupied}
#'   and \code{params$cells_empty} report how the points filled the grid, and
#'   a count that came from a profile or selection warns when fewer than
#'   three quarters of it are occupied.  For cells that follow the points,
#'   seed a Voronoi tessellation with
#'   \code{get_voronoi_seeds(method = "kmeans", n = <the count>, sample_points = <the points>)}.
#' @param cellsize Numeric cell size, in the units of the working CRS.  Read by
#'   \code{method = "hex"} and \code{"square"} only; the other two methods warn
#'   that it was ignored.  When both \code{cellsize} and \code{approx_n_cells}
#'   are given, \code{cellsize} wins and \code{approx_n_cells} is ignored with
#'   a logged warning; supply one or the other.  With a geographic \code{crs}
#'   it is in that CRS's degrees, and the grid is laid in degrees.
#' @param expand Buffer distance, in the working CRS's units, by which the
#'   Voronoi boundary is grown before the diagram is built. Applied by
#'   `method = "voronoi"` only; the `"hex"`, `"square"` and `"triangles"`
#'   methods ignore it (the value you passed is still echoed back in
#'   `params$expand`). With `clip = TRUE` the cells are clipped to the grown
#'   boundary, which is the one returned as `boundary`: cells reach `expand`
#'   beyond the study area, and a point up to `expand` outside it is indexed.
#' @param clip Logical; clip to boundary.
#' @param keep_duplicates Logical. Has no effect on the cells or the index:
#'   coincident points are merged before a Voronoi diagram or a Delaunay
#'   triangulation is built either way, and every one of them is indexed to
#'   the cell they share.
#' @param crs Optional target CRS. A projected CRS is the working CRS. A
#'   geographic one (EPSG:4326, say) is the CRS the result is returned in: the
#'   cells are built in the local projected CRS [ensure_projected()] picks for
#'   the points, indexed there, and then transformed with long edges
#'   densified, so Voronoi cells are nearest-point cells on the ground and
#'   grid cells are laid in metres rather than degrees. The exception is a hex
#'   or square grid sized by `cellsize`, which is in degrees and so is laid in
#'   degrees.
#' @param quiet Logical; suppress this function's progress \code{message()}s.
#'   It does not silence R warnings, nor the package's console log echo
#'   (see \code{\link{spatialkit_quiet}} for that). Default \code{FALSE}.
#' @return A list with components:
#'   \describe{
#'     \item{`cells`}{An sf polygon layer, one row per cell. It always carries
#'       a `cell_id` column; the `"hex"` and `"square"` methods additionally
#'       carry `poly_id`, which holds the same values.}
#'     \item{`index`}{Integer vector of `cell_id` values, one per row of
#'       `points_sf`, and `NA` for a point that falls inside no cell, that is,
#'       one outside the study area. Only a point within a thousandth of the
#'       median cell width of a cell is snapped to it. That covers points
#'       sitting exactly on a shared edge, and leaves points outside the study
#'       area as `NA`. A summary built from `index` therefore counts only the
#'       points the tessellation actually covers.}
#'     \item{`boundary`}{The boundary used (possibly derived and/or
#'       reprojected, and for `"voronoi"` grown by `expand`).}
#'     \item{`method`}{The method actually used.}
#'     \item{`params`}{The parameters the tessellation was built with, plus
#'       `snapped`, the record of that nearest-cell repair: a list with `n`,
#'       `which` (row positions in `points_sf`) and `distance` (how far
#'       outside every cell each sat, in CRS units). For `"hex"` and
#'       `"square"` also `cells_occupied` and `cells_empty`, the number of
#'       cells that hold at least one point and that hold none.}
#'   }
#' @family tessellation
#' @examples
#' library(sf)
#' set.seed(1)
#' pts <- st_as_sf(
#'   data.frame(x = 5e5 + runif(20, 0, 100), y = 5e6 + runif(20, 0, 100)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' tess <- build_tessellation(pts, method = "voronoi", quiet = TRUE)
#' tess$cells
#' @export
build_tessellation <- function(
    points_sf, boundary = NULL,
    method = c("voronoi", "triangles", "hex", "square"),
    approx_n_cells = NULL, cellsize = NULL, expand = 0,
    clip = TRUE, keep_duplicates = FALSE, crs = NULL, quiet = FALSE
) {
  .msg <- function(...) if (!quiet) message(...)
  method <- match.arg(method)
  .assert_sf(points_sf, c("POINT", "MULTIPOINT"), "points_sf")

  # The level-selection step's own answer is accepted here, so the count
  # need not be carried between the two calls by hand.
  n_cells <- .resolve_cell_count(approx_n_cells, "approx_n_cells", "build_tessellation")
  approx_n_cells      <- if (is.null(n_cells)) NULL else n_cells$n
  approx_n_cells_from <- if (is.null(n_cells)) NULL else n_cells$from

  # `approx_n_cells` and `cellsize` size the hex and square lattices and
  # nothing else: Voronoi grows one cell per input point and Delaunay one
  # triangle per neighbouring triple, so neither has a count to set.  Both used
  # to be dropped in silence, so `build_tessellation(pts, method = "voronoi",
  # approx_n_cells = 25)` returned one cell per observation -- the degenerate
  # case -- with nothing to say 25 had been asked for.  Warn rather than stop:
  # the call still produces a valid tessellation, just not the one intended.
  #
  # The `params` clause is voronoi-only on purpose.  That branch returns
  # create_voronoi_polygons()'s own list, which has no slot for either
  # argument; the triangles branch does echo `approx_n_cells` back (but not
  # `cellsize`), so claiming otherwise for it would be false.
  if (!method %in% c("hex", "square")) {
    ignored <- c(if (!is.null(approx_n_cells)) "approx_n_cells",
                 if (!is.null(cellsize)) "cellsize")
    if (length(ignored) > 0L)
      .warn_and_log(
        "build_tessellation(method = \"%s\") ignores the grid-sizing %s %s. %s",
        method,
        if (length(ignored) > 1L) "arguments" else "argument",
        paste(sprintf("`%s`", ignored), collapse = " and "),
        if (identical(method, "voronoi"))
          paste("`params` does not record the request either. Voronoi grows",
                "one cell per input point: to control the cell count, place",
                "seeds with get_voronoi_seeds() and tessellate those, or use",
                "method = \"hex\" or \"square\".")
        else
          paste("Delaunay produces one triangle per neighbouring triple, so",
                "the count follows from the points; use method = \"hex\" or",
                "\"square\" to set it."))
  }

  # --- CRS handling ---
  crs_out <- NULL
  if (!is.null(crs)) {
    points_sf <- .transform_or_stamp(points_sf, crs, "points_sf", "build_tessellation")
    if (!is.null(boundary))
      boundary <- .transform_or_stamp(boundary, crs, "boundary", "build_tessellation")
    # A geographic `crs` names the CRS the result is RETURNED in; the cells
    # are built in metres.  st_voronoi(), st_make_grid() and the hull buffer
    # all work on raw coordinates, so in degrees Voronoi cells stopped being
    # a nearest-point partition (21% of sampled locations at 55N sat in
    # another point's cell), grid cells were neither square nor equal-area,
    # and clipping them under s2 often stopped with "Edge 0 is degenerate"
    # or gave points inside the boundary an NA index.  Build in the local
    # projected CRS, index there, and transform at the end (finish() below).
    # An explicit `cellsize` is in the units of `crs`, degrees, so a lattice
    # sized by it is still laid in degrees, as asked.
    if (.is_geographic_crs(crs) &&
        !(method %in% c("hex", "square") && !is.null(cellsize))) {
      crs_out   <- sf::st_crs(points_sf)
      points_sf <- ensure_projected(points_sf)
      if (!is.null(boundary)) boundary <- .align_crs(boundary, points_sf)
    }
  } else {
    # ONE decision for points and boundary together.  A CRS-less layer goes
    # through the same lon/lat heuristic as everywhere else; if it is taken as
    # lon/lat, a CRS-less boundary is given the same interpretation before
    # being aligned, and if it is not, the boundary stays in the same unnamed
    # space (crs_arg = NA below tells create_grid_polygons() not to project it
    # on its own).  Deciding separately for the two -- points left as they
    # were, boundary projected inside create_grid_polygons() -- put grid and
    # points in different CRSs and hex/square died in st_intersects() on
    # input that voronoi/triangles accepted.
    if (.is_longlat(points_sf)) {
      .msg("build_tessellation(): projecting points to a local UTM CRS.")
      points_sf <- ensure_projected(points_sf)
    } else if (is.na(sf::st_crs(points_sf))) {
      points_sf <- ensure_projected(points_sf)
    }
    if (!is.null(boundary)) {
      # Only a POSITIVE assumption is a CRS.  ensure_projected() records
      # "none" for CRS-less points it left planar, and st_crs("none") is an
      # error ("invalid crs: none"), so every method failed on CRS-less
      # planar points with a CRS-less boundary -- including the documented
      # boundary = clip_target_for(pts) -- and such data could not be gridded
      # at all.  Left alone, the two stay in the same unnamed space.
      assumed <- attr(points_sf, "crs_assumed")
      if (identical(assumed, "EPSG:4326") && is.na(sf::st_crs(boundary)))
        boundary <- sf::st_set_crs(boundary, sf::st_crs(assumed))
      # One side with no CRS takes the other's; .align_crs() leaves a
      # CRS-less side as it is, and sf then stopped every method with its
      # bare "st_crs(x) == st_crs(y) is not TRUE".
      pair <- .resolve_crsless_pair(points_sf, boundary, "build_tessellation")
      points_sf <- pair$points; boundary <- pair$boundary
      boundary <- .align_crs(boundary, points_sf)
    }
  }
  finish <- function(res) {
    if (is.null(crs_out)) return(res)
    res$cells    <- .to_output_crs(res$cells, crs_out)
    res$boundary <- .to_output_crs(res$boundary, crs_out)
    res
  }
  if (!is.null(boundary)) {
    if (!any(sf::st_geometry_type(boundary) %in% c("POLYGON", "MULTIPOLYGON")))
      stop("build_tessellation(): `boundary` must be polygonal.")
    boundary <- .safe_make_valid(boundary)
  }

  # A CRS-less `points_sf` yields NA_crs_, which is a list rather than NULL and
  # so is not treated as "no CRS supplied" downstream -- create_voronoi_polygons()
  # would call st_transform() on a CRS-less object and fail.  Normalise to NULL.
  pts_crs <- sf::st_crs(points_sf)
  crs_arg <- if (is.na(pts_crs)) NULL else pts_crs

  # ---- Voronoi ----
  if (identical(method, "voronoi")) {
    return(finish(create_voronoi_polygons(
      points_sf = points_sf, boundary = boundary, expand = expand,
      clip = clip, keep_duplicates = keep_duplicates,
      crs = crs_arg, quiet = quiet
    )))
  }

  # ---- Hex / Square ----
  if (method %in% c("hex", "square")) {
    if (is.null(boundary))
      stop("build_tessellation(): `boundary` is required for hex/square grids.")
    # Pass the points' CRS through: without it create_grid_polygons() would
    # re-project the boundary on its own, leaving grid and points in different
    # CRSs and breaking the st_intersects() below.
    grid <- create_grid_polygons(
      boundary = boundary, target_cells = approx_n_cells,
      type = method, cellsize = cellsize, clip = clip,
      # NA (not NULL) when the points are CRS-less: keep the boundary in the
      # points' unnamed space instead of letting create_grid_polygons() run
      # its own projection heuristic on it.
      crs = if (is.null(crs_arg)) sf::NA_crs_ else crs_arg,
      quiet = quiet
    )
    # Build point-to-cell index for grid tessellations
    id_col <- if ("poly_id" %in% names(grid)) "poly_id" else "cell_id"
    if (!id_col %in% names(grid)) {
      grid$cell_id <- seq_len(nrow(grid))
      id_col <- "cell_id"
    }
    # `cell_id` is the id column every other method returns, so keep it on the
    # grid too rather than dropping it after indexing: downstream helpers (and
    # plot_tessellation_map(fill_col = "cell_id")) can then treat all four
    # methods alike.  `poly_id` is retained for backward compatibility.
    if (id_col != "cell_id") {
      grid$cell_id <- grid[[id_col]]
    }
    index   <- .build_point_cell_index(points_sf, grid)
    snapped <- attr(index, "snapped")
    attr(index, "snapped") <- NULL

    # How the points filled the lattice.  A count read off resolution_profile()
    # or select_resolution() is a number of k-means cells: all occupied, small
    # where the points are dense.  The same number of EQUAL cells over the
    # boundary leaves the gaps between clusters empty -- on six clusters in a
    # square, 33 of 77 hexagons held a point for a profile count of 56 -- and
    # nothing said so.  Report occupancy always, and warn when a count that
    # came from a profile or selection leaves under three quarters occupied
    # (evenly spread points fill 1.2 times the count).
    occupied <- length(unique(index[!is.na(index)]))
    if (!is.null(approx_n_cells_from) && is.null(cellsize) &&
        occupied < 0.75 * approx_n_cells)
      .warn_and_log(paste0(
        "build_tessellation(method = \"%s\"): %d of the %d cells hold a point, ",
        "against %s cells from %s. That count is of k-means cells, every one ",
        "occupied and dense where the points are; a lattice of equal cells ",
        "over clustered points leaves many empty. For cells that follow the ",
        "points, place seeds with get_voronoi_seeds(method = \"kmeans\", n = ",
        "<the count>, sample_points = <the points>) and use method = ",
        "\"voronoi\"."),
        method, occupied, nrow(grid), format(approx_n_cells), approx_n_cells_from)

    return(finish(list(
      cells = grid, index = index, boundary = boundary, method = method,
      params = list(approx_n_cells = approx_n_cells,
                    approx_n_cells_from = approx_n_cells_from,
                    cellsize = cellsize,
                    clip = clip, keep_duplicates = keep_duplicates,
                    expand = expand, snapped = snapped,
                    cells_occupied = occupied,
                    cells_empty = nrow(grid) - occupied)
    )))
  }

  # ---- Delaunay triangles ----
  if (identical(method, "triangles")) {
    pts <- if (isTRUE(keep_duplicates)) points_sf else .dedup_points(points_sf)
    if (nrow(pts) < 3L) stop("build_tessellation(triangles): need at least 3 unique points.")
    coords <- sf::st_coordinates(pts)[, 1:2, drop = FALSE]
    # Points on one line have no triangulation.  qhull returned a 0 x 3
    # matrix without an error, the fallback below then logged that
    # delaunayn() had failed, which it had not, and the call returned no cells
    # and an index of NAs (a single transect, 10 points on y = 0.5x).  Say so,
    # as the three-point rule above does.
    if (qr(scale(coords, scale = FALSE))$rank < 2L)
      stop("build_tessellation(triangles): the points are collinear, so their ",
           "Delaunay triangulation is empty. Use method = \"voronoi\", which ",
           "handles points on a line.", call. = FALSE)

    tri_sfc <- NULL
    why_geos <- "package 'geometry' is not installed"
    if (requireNamespace("geometry", quietly = TRUE)) {
      # Centre the points before qhull sees them.  It triangulates by lifting
      # each point onto x^2 + y^2, and at projected magnitudes (a UTM
      # northing near 5e6, or 9e6 south of the equator) the lift has no
      # precision left to separate points a few metres apart: they were
      # dropped as "coplanar" and never became vertices, with nothing said.
      # 200 points over 100 m at (5e5, 5e6) gave 26 triangles instead of 386,
      # yet every point still fell in one.  Translation does not change the
      # Delaunay triangulation.  The shift is the bbox midpoint, not the
      # mean: min and max do not depend on row order, so a permuted layer
      # hands qhull bit-identical coordinates and gets the same triangles
      # and cell_ids, as before.  Rings below use the original coordinates,
      # so the output coordinates are untouched.
      ctr <- (apply(coords, 2L, min) + apply(coords, 2L, max)) / 2
      tri_idx <- try(geometry::delaunayn(sweep(coords, 2L, ctr)), silent = TRUE)
      why_geos <- if (inherits(tri_idx, "try-error"))
        sprintf("geometry::delaunayn() failed (%s)",
                trimws(conditionMessage(attr(tri_idx, "condition"))))
      else "geometry::delaunayn() returned no triangles"
      if (!inherits(tri_idx, "try-error") && length(tri_idx)) {
        polys <- vector("list", nrow(tri_idx))
        for (i in seq_len(nrow(tri_idx))) {
          idx <- tri_idx[i, ]
          ring <- rbind(coords[idx, , drop = FALSE], coords[idx[1], , drop = FALSE])
          
          signed_area <- sum(
            ring[-nrow(ring), 1] * ring[-1, 2] -
            ring[-1, 1] * ring[-nrow(ring), 2]
          ) / 2
          if (signed_area < 0) ring <- ring[rev(seq_len(nrow(ring))), ]

          polys[[i]] <- sf::st_polygon(list(ring))
        }
        tri_sfc <- sf::st_sfc(polys, crs = sf::st_crs(pts))
      }
    }

    if (is.null(tri_sfc)) {
      # Name the reason that applies: the message used to blame a missing
      # package or a failure whichever had happened.
      .log_warn(paste0(
        "build_tessellation(triangles): %s. Falling back to GEOS via ",
        "sf::st_triangulate() instead of qhull; the result is still the ",
        "Delaunay triangulation of the input points, but degenerate (e.g. ",
        "co-circular) configurations may be resolved differently."), why_geos)
      # st_triangulate() accepts a MULTIPOINT and returns the true Delaunay
      # triangulation of it.  Triangulating the convex-hull POLYGON instead
      # (as this used to) discards every interior point.
      tri_sfc <- sf::st_triangulate(sf::st_union(sf::st_geometry(pts)))
      tri_sfc <- sf::st_collection_extract(tri_sfc, "POLYGON", warn = FALSE)
      tri_sfc <- sf::st_sfc(tri_sfc, crs = sf::st_crs(pts))
      if (length(tri_sfc) == 0L)
        stop("build_tessellation(triangles): the Delaunay triangulation of ",
             "these points is empty (they are collinear to within rounding). ",
             "Use method = \"voronoi\".", call. = FALSE)
    }

    tri_sf <- sf::st_sf(geometry = .safe_make_valid(tri_sfc))
    if (!is.null(boundary) && isTRUE(clip)) {
      # Union first: a multi-feature boundary would otherwise split straddling
      # triangles into one row per boundary feature and graft the boundary's
      # attribute columns onto the result.
      clip_to <- sf::st_union(sf::st_geometry(boundary))
      tri_sf <- suppressWarnings(sf::st_intersection(tri_sf, clip_to))
      # Keep only POLYGON/MULTIPOLYGON (drop slivers) and non-empty
      gtypes <- as.character(sf::st_geometry_type(tri_sf, by_geometry = TRUE))
      keep <- gtypes %in% c("POLYGON", "MULTIPOLYGON") & !sf::st_is_empty(tri_sf)
      tri_sf <- tri_sf[keep, , drop = FALSE]
    }
    tri_sf$cell_id <- seq_len(nrow(tri_sf))

    # Build point-to-cell index for triangles
    index   <- .build_point_cell_index(points_sf, tri_sf)
    snapped <- attr(index, "snapped")
    attr(index, "snapped") <- NULL

    return(finish(list(
      cells = tri_sf, index = index, boundary = boundary, method = "triangles",
      params = list(clip = clip, approx_n_cells = approx_n_cells,
                    approx_n_cells_from = approx_n_cells_from,
                    keep_duplicates = keep_duplicates, expand = expand,
                    snapped = snapped)
    )))
  }

  stop("build_tessellation(): unknown method.")
}
