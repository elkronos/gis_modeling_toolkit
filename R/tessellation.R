# -----------------------------------------------------------------------------
# Clip Target
# -----------------------------------------------------------------------------

#' Build a polygonal clip target from points and/or a boundary
#'
#' Resolves the single polygon that every tessellation method clips against.
#' With a `boundary` it is that boundary (optionally buffered by `expand`);
#' without one it is the convex hull of `points_sf`, again optionally buffered.
#' Reach for it when you want to see or reuse the exact clip target
#' [build_tessellation()] will apply — for instance to check that a study-area
#' polygon actually contains the observations before tessellating, or to pass
#' the same envelope to [create_voronoi_polygons()] and
#' [create_grid_polygons()] so that two tessellations of one dataset cover
#' identical ground.
#'
#' @param points_sf An sf object with POINT/MULTIPOINT geometry.
#' @param boundary Optional polygonal sf object.
#' @param expand Numeric expansion distance or fraction (0–1 = fraction of
#'   extent). Absolute values are expressed in the units of the CRS the clip
#'   target is built in. Because [sf::st_buffer()] interprets `dist` as
#'   **metres** for lon/lat data while the fraction-of-extent form is derived
#'   from a bounding box measured in **degrees**, lon/lat input is projected to
#'   a local projected CRS first (see `@return`), so both forms agree.
#' @param quiet Logical; suppress messages.
#' @return An sf polygon layer representing the clip target. For lon/lat input
#'   the layer is returned in the automatically selected local projected CRS,
#'   not the input CRS; a message reports this unless `quiet = TRUE`.
#' @export
clip_target_for <- function(points_sf, boundary = NULL, expand = 0, quiet = FALSE) {
  .msg <- function(...) if (!quiet) message(...)
  .assert_sf(points_sf, c("POINT", "MULTIPOINT"), "points_sf")

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
#' @return Integer vector of cell_id values, one per row of `pts`.
#' @keywords internal
#' @noRd
.build_point_cell_index <- function(pts, cells_sf) {
  n <- nrow(pts)
  if (n == 0L || nrow(cells_sf) == 0L) return(rep(NA_integer_, n))

  cell_ids <- cells_sf$cell_id
  hits <- sf::st_intersects(pts, cells_sf)
  
  index <- vapply(hits, function(row_idx) {
    row_idx <- as.integer(row_idx)
    if (length(row_idx) == 0L)     return(NA_integer_)
    if (length(row_idx) == 1L)     return(cell_ids[row_idx])
    # Ambiguous: pick the cell with the smallest cell_id
    min(cell_ids[row_idx])
  }, integer(1))

  # Fallback for unmatched points: assign to nearest cell
  unmatched <- which(is.na(index))
  if (length(unmatched) > 0L) {
    near <- suppressWarnings(sf::st_nearest_feature(pts[unmatched, ], cells_sf))
    valid <- is.finite(near)
    index[unmatched[valid]] <- cell_ids[near[valid]]
  }

  index
}

# -----------------------------------------------------------------------------
# Voronoi Tessellation
# -----------------------------------------------------------------------------

#' Create Voronoi polygons from points with robust CRS and optional clipping
#'
#' Assigns every location in the study area to its nearest input point, giving
#' one cell per point. This is the tessellation to reach for when the
#' observations themselves define the regions of interest — sampling sites,
#' monitoring stations, service points — because cell size then adapts to
#' sampling density instead of being imposed by a fixed grid: dense areas get
#' small cells and sparse areas large ones. Prefer [create_grid_polygons()]
#' instead when you need equal-area cells or a resolution independent of where
#' the data happen to be.
#'
#' The heavy lifting is [sf::st_voronoi()]; what this adds is the surrounding
#' bookkeeping — projecting lon/lat input, building and buffering an envelope
#' so edge cells are bounded, clipping to `boundary`, restoring the
#' point-to-cell correspondence that `st_voronoi()` scrambles, and stamping
#' stable `cell_id` values.
#'
#' @param points_sf An sf object with POINT/MULTIPOINT geometries.
#' @param boundary Optional polygonal sf object.
#' @param expand Numeric; absolute buffer distance for the envelope.
#' @param clip Logical; intersect cells with boundary.
#' @param keep_duplicates Logical; keep coincident points for graph construction.
#' @param crs Optional target CRS.
#' @param quiet Logical; suppress messages.
#' @return A list with cells, index, boundary, method, params.
#' @family tessellation
#' @examples
#' library(sf)
#' set.seed(1)
#' pts <- st_as_sf(
#'   data.frame(x = runif(15, 0, 100), y = runif(15, 0, 100)),
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
  if (!is.null(crs)) {
    pts <- .transform_or_stamp(pts, crs, "points_sf", "create_voronoi_polygons")
    if (!is.null(boundary))
      boundary <- .transform_or_stamp(boundary, crs, "boundary", "create_voronoi_polygons")
  } else {
    if (.is_longlat(pts)) pts <- ensure_projected(pts)
    if (!is.null(boundary)) boundary <- .align_crs(boundary, pts)
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

  boundary_expanded <- if (isTRUE(is.numeric(expand)) && expand > 0) {
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
    clip_to <- if (isTRUE(is.numeric(expand)) && expand > 0) boundary_expanded else boundary
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

  # Build point → cell index using cell_id values
  index <- .build_point_cell_index(pts, cells)

  list(
    cells    = cells,
    index    = index,
    boundary = boundary,
    method   = "voronoi",
    params   = list(clip = clip, expand = expand, keep_duplicates = keep_duplicates)
  )
}

# -----------------------------------------------------------------------------
# Grid Tessellations (Hex / Square)
# -----------------------------------------------------------------------------

#' Create square or hexagonal grid polygons over a boundary
#'
#' Lays a regular grid of equal-area cells over `boundary` and clips it to that
#' boundary. Reach for this rather than [create_voronoi_polygons()] when cell
#' size should be a decision you make — because you need per-cell rates
#' comparable across the map, or a resolution that stays fixed as the sample
#' grows — instead of one dictated by where the observations happen to be.
#' Hexagons (`type = "hex"`) avoid the axis-aligned artefacts of squares and
#' give every cell the same distance to all six neighbours, which matters for
#' anything that reads neighbourhoods.
#'
#' Size the grid with exactly one of `target_cells` (roughly how many cells you
#' want, the package derives the rest), `cellsize` (a fixed edge length in CRS
#' units) or `n` (a fixed number of columns and rows). See `@param cellsize`
#' for what happens when more than one is given.
#'
#' @param boundary Polygonal sf or sfc object.
#' @param target_cells Optional approximate desired number of cells. For hex
#'   grids the count is adjusted for hexagonal packing density, but the final
#'   cell count after clipping to an irregular boundary may deviate
#'   substantially from the requested value.
#' @param type Grid type: `"square"` (the default) or `"hex"`.
#' @param cellsize Optional numeric cell size (length 1 or 2), in the units of
#'   the working CRS. Takes precedence over `n`: if both are supplied,
#'   `cellsize` is used, `n` is ignored and a warning is logged. Supply exactly
#'   one of `target_cells`, `cellsize` and `n`.
#' @param n Optional grid resolution (integer, length 1 or 2) giving the number
#'   of columns and rows to divide the boundary's bounding box into; the cell
#'   size is derived from it. Applies to square grids only; [sf::st_make_grid()]
#'   derives hexagon placement from `cellsize` alone. Ignored (with a logged
#'   warning) when `cellsize` is also supplied — passing both would otherwise
#'   truncate the grid to `n[1]` x `n[2]` cells anchored at the bounding-box
#'   corner, covering only part of the boundary.
#' @param clip Logical; clip grid to boundary.
#' @param crs Optional target CRS. When `NULL` (default) the boundary is
#'   projected with [ensure_projected()], which changes the CRS of the returned
#'   grid; a message reports this unless `quiet = TRUE`.
#' @param quiet Logical; suppress messages.
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
    cellsize = NULL, n = NULL, clip = TRUE, crs = NULL, quiet = FALSE
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

  if (!is.null(crs)) {
    boundary <- .transform_or_stamp(boundary, crs, "boundary", "create_grid_polygons")
  } else {
    crs_before <- sf::st_crs(boundary)
    boundary <- ensure_projected(boundary)
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
    cellsize <- c(w / nx, h / ny)
  }

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
  grid_sf
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
#' \code{\link{summarize_by_cell}()} consume.  Use it rather than the
#' individual constructors whenever you might want to compare methods: the
#' return shape does not change with \code{method}, so swapping
#' \code{"voronoi"} for \code{"hex"} costs one argument.
#'
#' Which method to reach for.  \code{"voronoi"} gives one cell per point, so
#' resolution follows sampling density -- the choice when the observations
#' themselves define the regions.  \code{"hex"} and \code{"square"} give
#' equal-area cells on a fixed grid, so cell size is a decision you make rather
#' than one the data makes for you; hexagons avoid the axis-aligned artefacts of
#' squares and have uniform neighbour distances.  \code{"triangles"} returns
#' the Delaunay triangulation, useful for interpolation and adjacency work
#' rather than as an aggregation unit.  \code{\link{determine_optimal_levels}()}
#' will suggest a cell count from the spatial structure of the data.
#'
#' @param points_sf An sf object with POINT/MULTIPOINT geometry.
#' @param boundary Polygonal sf/sfc study area. **Required** for
#'   `method = "hex"` and `method = "square"`, which have no extent of their
#'   own to lay a grid over and error without it; supply the study-area polygon,
#'   or build one from the points with [clip_target_for()]. **Optional** for
#'   `method = "voronoi"` and `method = "triangles"`, which derive their extent
#'   from the points themselves and use `boundary` only to clip the result when
#'   `clip = TRUE`.
#' @param method One of "voronoi", "triangles", "hex", "square".
#' @param approx_n_cells Approximate number of cells (grid methods). For hex
#'   grids the target is adjusted for packing density; the actual count after
#'   clipping to an irregular boundary may differ noticeably.
#' @param cellsize Numeric cell size (grid methods).
#' @param expand Buffer distance for the Voronoi envelope. Applied by
#'   `method = "voronoi"` only; the `"hex"`, `"square"` and `"triangles"`
#'   methods ignore it (the value you passed is still echoed back in
#'   `params$expand`).
#' @param clip Logical; clip to boundary.
#' @param keep_duplicates Logical; keep duplicate points.
#' @param crs Optional target CRS.
#' @param quiet Logical; suppress messages.
#' @return A list with components:
#'   \describe{
#'     \item{`cells`}{An sf polygon layer, one row per cell. It always carries
#'       a `cell_id` column; the `"hex"` and `"square"` methods additionally
#'       carry `poly_id`, which holds the same values.}
#'     \item{`index`}{Integer vector of `cell_id` values, one per row of
#'       `points_sf`.}
#'     \item{`boundary`}{The boundary used (possibly derived and/or reprojected).}
#'     \item{`method`}{The method actually used.}
#'     \item{`params`}{The parameters the tessellation was built with.}
#'   }
#' @family tessellation
#' @examples
#' library(sf)
#' set.seed(1)
#' pts <- st_as_sf(
#'   data.frame(x = runif(20, 0, 100), y = runif(20, 0, 100)),
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

  # --- CRS handling ---
  if (!is.null(crs)) {
    points_sf <- .transform_or_stamp(points_sf, crs, "points_sf", "build_tessellation")
    if (!is.null(boundary))
      boundary <- .transform_or_stamp(boundary, crs, "boundary", "build_tessellation")
  } else {
    if (.is_longlat(points_sf)) {
      .msg("build_tessellation(): projecting points to a local UTM CRS.")
      points_sf <- ensure_projected(points_sf)
    }
    if (!is.null(boundary)) boundary <- .align_crs(boundary, points_sf)
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
    return(create_voronoi_polygons(
      points_sf = points_sf, boundary = boundary, expand = expand,
      clip = clip, keep_duplicates = keep_duplicates,
      crs = crs_arg, quiet = quiet
    ))
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
      type = method, cellsize = cellsize, clip = clip, crs = crs_arg,
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
    index <- .build_point_cell_index(points_sf, grid)

    return(list(
      cells = grid, index = index, boundary = boundary, method = method,
      params = list(approx_n_cells = approx_n_cells, cellsize = cellsize,
                    clip = clip, keep_duplicates = keep_duplicates,
                    expand = expand)
    ))
  }

  # ---- Delaunay triangles ----
  if (identical(method, "triangles")) {
    pts <- if (isTRUE(keep_duplicates)) points_sf else .dedup_points(points_sf)
    if (nrow(pts) < 3L) stop("build_tessellation(triangles): need at least 3 unique points.")
    coords <- sf::st_coordinates(pts)[, 1:2, drop = FALSE]

    tri_sfc <- NULL
    if (requireNamespace("geometry", quietly = TRUE)) {
      tri_idx <- try(geometry::delaunayn(coords), silent = TRUE)
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
      .log_warn(
        "build_tessellation(triangles): package 'geometry' unavailable or delaunayn() failed. Falling back to GEOS via sf::st_triangulate() instead of qhull; the result is still the Delaunay triangulation of the input points, but degenerate (e.g. co-circular) configurations may be resolved differently."
      )
      # st_triangulate() accepts a MULTIPOINT and returns the true Delaunay
      # triangulation of it.  Triangulating the convex-hull POLYGON instead
      # (as this used to) discards every interior point.
      tri_sfc <- sf::st_triangulate(sf::st_union(sf::st_geometry(pts)))
      tri_sfc <- sf::st_collection_extract(tri_sfc, "POLYGON", warn = FALSE)
      tri_sfc <- sf::st_sfc(tri_sfc, crs = sf::st_crs(pts))
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
    index <- .build_point_cell_index(points_sf, tri_sf)

    return(list(
      cells = tri_sf, index = index, boundary = boundary, method = "triangles",
      params = list(clip = clip, approx_n_cells = approx_n_cells,
                    keep_duplicates = keep_duplicates, expand = expand)
    ))
  }

  stop("build_tessellation(): unknown method.")
}
