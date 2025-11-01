# =============================================================================
# Geospatial Modeling Toolkit (Robust & Modular)
# - Tessellation: Voronoi / Hex / Square
# - Assignment: map features to cells (robust CRS handling)
# - Seeding: k-means / random / provided points
# - Models: GWR (spgwr) and Bayesian spatial (spBayes; exp/spherical/Matérn)
# - Level selection: elbow heuristic on WSS
# - Utilities: pointization for non-point geometries, summaries, plotting
#
# Works with arbitrary user datasets (not tied to any specific region).
#
# Dependencies: logger, sf, sp, spgwr, spBayes, deldir, ggplot2, dplyr, tidyr, mvtnorm
# =============================================================================

suppressPackageStartupMessages({
  library(logger)
  library(sf)
  library(sp)
  library(spgwr)
  library(spBayes)
  library(deldir)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(mvtnorm)
})

# -----------------------------------------------------------------------------
# Logging
# -----------------------------------------------------------------------------
log_appender(appender_file("model_log.log"))
log_threshold(INFO)

# -----------------------------------------------------------------------------
# Small utilities
# -----------------------------------------------------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

.is_longlat <- function(x) {
  cr <- sf::st_crs(x)
  if (is.na(cr)) return(FALSE)
  out <- try(sf::st_is_longlat(cr), silent = TRUE)
  isTRUE(out)
}

# -----------------------------------------------------------------------------
# CRS & Geometry Utilities
# -----------------------------------------------------------------------------

#' Choose a locally appropriate projected CRS (UTM by centroid, hemisphere-aware)
#'
#' Computes a “best guess” projected coordinate reference system (CRS) for an
#' `sf`/`sfc` object. If the input is already projected (i.e., not long/lat),
#' its CRS is returned unchanged. If the input is in geographic coordinates
#' (lon/lat), a UTM CRS is chosen from the centroid longitude and hemisphere:
#' EPSG `326##` for the northern hemisphere and `327##` for the southern,
#' where `##` is the UTM zone. When a lon/lat CRS cannot be confirmed, the
#' function falls back to Web Mercator (EPSG:3857). If `x` is not an `sf`/`sfc`
#' object, `sf::NA_crs_` is returned.
#'
#' @details
#' The decision flow is:
#'
#' 1. If `x` is not `sf`/`sfc`, return `sf::NA_crs_`.
#' 2. If `x` has a non-long/lat CRS already, return that CRS.
#' 3. Otherwise, try to get lon/lat coordinates by transforming to EPSG:4326.
#' 4. If lon/lat cannot be confirmed, return EPSG:3857 (fallback).
#' 5. Compute the centroid of the unioned geometry in lon/lat, derive the UTM
#'    zone from longitude, and select:
#'    - EPSG: `32600 + zone` (north),
#'    - EPSG: `32700 + zone` (south).
#'
#' Fallbacks to EPSG:3857 also occur if the centroid cannot be computed or does
#' not yield valid numeric coordinates.
#'
#' @note
#' - UTM is optimal for local analyses but less suitable for very large extents,
#'   multi-zone datasets, areas near the poles, or geometries crossing the
#'   antimeridian. In such cases the 3857 fallback or a hand-picked local CRS
#'   may be preferable.
#' - This function **returns** a CRS object; it does not modify `x`. Apply the
#'   result with `st_transform(x, .pick_local_projected_crs(x))` (or use a
#'   wrapper like `ensure_projected()` if available).
#'
#' @param x An `sf` or `sfc` object whose geometry extent is used to infer a
#'   suitable local projected CRS.
#'
#' @return An `sf` CRS object (class `crs`) to use for local distance/area
#'   computations. Possible values include the existing projected CRS of `x`,
#'   a UTM CRS (EPSG 326##/327##), EPSG:3857 as a fallback, or `sf::NA_crs_`
#'   if `x` is not `sf`/`sfc`.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' crs_guess <- .pick_local_projected_crs(nc)
#' nc_proj <- st_transform(nc, crs_guess)
#' st_crs(nc_proj)
#' }
#'
#' @seealso
#' \code{\link[sf]{st_crs}}, \code{\link[sf]{st_transform}},
#' \code{\link[sf]{st_is_longlat}}
#'
#' @keywords internal
.pick_local_projected_crs <- function(x) {
  if (!inherits(x, c("sf", "sfc"))) return(sf::NA_crs_)
  crs <- sf::st_crs(x)
  if (!is.na(crs) && !.is_longlat(x)) return(crs)
  
  x_ll <- tryCatch({
    if (is.na(crs)) stop("No CRS set.")
    sf::st_transform(x, 4326)
  }, error = function(e) x)
  
  if (is.na(sf::st_crs(x_ll)) || !.is_longlat(x_ll)) return(sf::st_crs(3857))
  
  ctr <- sf::st_coordinates(sf::st_centroid(sf::st_union(sf::st_geometry(x_ll))))
  if (!is.numeric(ctr) || length(ctr) < 2) return(sf::st_crs(3857))
  lon <- ctr[1]; lat <- ctr[2]
  
  zone <- floor((lon + 180) / 6) + 1
  zone <- max(1, min(60, zone))
  epsg <- if (!is.na(lat) && lat < 0) 32700 + zone else 32600 + zone
  sf::st_crs(epsg)
}

#' Ensure an object uses a projected CRS (optionally match a target CRS)
#'
#' Guarantees that an `sf`/`sfc` object is in a projected (planar) CRS so that
#' distances/areas behave sensibly. If `target_crs` is provided, the object is
#' set/transformed to that CRS. Otherwise:
#'
#' - If `x` already has a non-long/lat CRS, it is returned unchanged.
#' - If `x` has a long/lat CRS, a locally appropriate projected CRS is chosen
#'   via [`.pick_local_projected_crs()`] and `x` is transformed to it.
#' - If `x` has **no CRS** and its bounding box appears to be lon/lat
#'   (within `[-180,180] × [-90,90]`), the CRS is assumed to be WGS84
#'   (EPSG:4326) and then projected using [`.pick_local_projected_crs()`].
#' - If `x` has **no CRS** and does **not** look like lon/lat, it is returned
#'   unchanged (we avoid guessing a projected CRS blindly).
#'
#' @details
#' When `target_crs` is supplied, it may be anything accepted by
#' `sf::st_crs()` (e.g. EPSG code, proj4string, WKT, or another `sf`/`sfc`
#' object). If `x` has no CRS, the function **sets** `target_crs` on `x`
#' (no reprojection is possible). If `x` has a CRS different from
#' `target_crs`, `x` is **transformed** to `target_crs`.
#'
#' The heuristic for missing-CRS lon/lat is intentionally conservative:
#' only if all bbox limits fall within plausible geographic degrees will
#' WGS84 be assumed; otherwise the object is left untouched.
#'
#' @param x An `sf` or `sfc` object to check/transform. If not `sf`/`sfc`,
#'   the object is returned as-is.
#' @param target_crs Optional. A CRS specification to match, accepted by
#'   `sf::st_crs()` (e.g., `4326`, `"EPSG:3857"`, a WKT string, or an object
#'   carrying a CRS like an `sf`/`sfc`).
#'
#' @return The same class as `x` (`sf` or `sfc`), with a projected CRS ensured.
#'   Depending on inputs, this may be:
#'   - `x` transformed to `target_crs`,
#'   - `x` with `target_crs` set (if `x` had no CRS),
#'   - `x` transformed to a locally chosen projected CRS, or
#'   - `x` unchanged (if already projected, or missing CRS that doesn't look
#'     like lon/lat).
#'
#' @examples
#' \dontrun{
#' library(sf)
#' # Example data
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#'
#' # 1) Already projected? Return unchanged.
#' nc_proj <- st_transform(nc, 32119)  # NAD83 / NC (ftUS)
#' st_crs(ensure_projected(nc_proj))
#'
#' # 2) Long/lat -> auto-pick local projected CRS
#' nc_ll <- st_transform(nc, 4326)
#' nc_auto <- ensure_projected(nc_ll)
#' st_is_longlat(nc_auto)   # FALSE
#'
#' # 3) Missing CRS but bbox looks like lon/lat -> assume 4326 then project
#' nc_drop <- st_set_crs(nc_ll, NA)
#' nc_fixed <- ensure_projected(nc_drop)
#' st_is_longlat(nc_fixed)  # FALSE
#'
#' # 4) Force a particular target CRS
#' nc_web <- ensure_projected(nc_ll, target_crs = 3857)
#' st_crs(nc_web)$epsg      # 3857
#'
#' # 5) Match a boundary's CRS
#' boundary <- st_union(nc)
#' nc_matched <- ensure_projected(nc_ll, boundary)
#' identical(st_crs(nc_matched), st_crs(boundary))  # TRUE
#' }
#'
#' @seealso
#' [sf::st_crs()], [sf::st_transform()], [sf::st_is_longlat()],
#' [`.pick_local_projected_crs()`]
#'
#' @family CRS helpers
ensure_projected <- function(x, target_crs = NULL) {
  if (!(inherits(x, "sf") || inherits(x, "sfc"))) return(x)
  
  if (!is.null(target_crs)) {
    tcrs <- sf::st_crs(target_crs)
    if (is.na(tcrs)) return(x)
    xcrs <- sf::st_crs(x)
    if (is.na(xcrs)) {
      sf::st_crs(x) <- tcrs
      return(x)
    }
    if (!identical(xcrs, tcrs)) x <- sf::st_transform(x, tcrs)
    return(x)
  }
  
  xcrs <- sf::st_crs(x)
  if (is.na(xcrs)) {
    bb <- try(sf::st_bbox(x), silent = TRUE)
    if (!inherits(bb, "try-error") &&
        all(is.finite(bb)) &&
        bb["xmin"] >= -180 && bb["xmax"] <= 180 &&
        bb["ymin"] >= -90  && bb["ymax"] <= 90) {
      sf::st_crs(x) <- sf::st_crs(4326)
      tr <- .pick_local_projected_crs(x)
      x <- sf::st_transform(x, tr)
    } else {
      # Likely projected but unknown; leave as-is
    }
    return(x)
  }
  
  if (.is_longlat(x)) {
    tr <- .pick_local_projected_crs(x)
    x <- sf::st_transform(x, tr)
  }
  x
}

#' Harmonize (align) coordinate reference systems (CRS) of two `sf`/`sfc` objects
#'
#' Ensures two spatial objects use a common CRS with minimal guessing:
#'
#' - If one object has a missing CRS and the other has a defined CRS, the
#'   missing one is **assigned** the other's CRS (no reprojection involved).
#' - If both objects have defined but different CRSs, `b` is **transformed**
#'   to the CRS of `a` using `sf::st_transform()`.
#' - If both CRSs are defined and identical, both inputs are returned unchanged.
#' - If both CRSs are missing, both are returned unchanged (no heuristics applied).
#'
#' This helper is useful before spatial joins, intersections, plotting overlays,
#' or any operation that requires a shared CRS.
#'
#' @param a An `sf` or `sfc` object whose CRS will be treated as the reference
#'   when both inputs have defined CRSs.
#' @param b An `sf` or `sfc` object to align with `a`'s CRS as needed.
#'
#' @return A named `list` with two elements:
#' \describe{
#'   \item{`a`}{The (possibly CRS-assigned) first input, unchanged in geometry.}
#'   \item{`b`}{The second input, either CRS-assigned to match `a` or transformed
#'               to `a`'s CRS when both CRSs are defined and differ.}
#' }
#'
#' @details
#' When only one object has a defined CRS, that CRS is copied to the other
#' (i.e., **setting** the CRS, not reprojecting). This mirrors common practice
#' when two datasets are known to be in the same spatial frame but one has lost
#' its CRS metadata.
#'
#' If both CRSs are defined and differ, only `b` is reprojected to `a`'s CRS,
#' making `a` the reference. Choose which argument you pass as `a` accordingly.
#'
#' The function does not attempt to infer CRSs when both are missing; use
#' [`ensure_projected()`] or your own domain knowledge to set CRSs first.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' a <- st_transform(nc, 32119)     # NAD83 / NC
#' b <- st_transform(nc, 4326)      # WGS84 lon/lat
#'
#' out <- harmonize_crs(a, b)
#' st_crs(out$a) == st_crs(out$b)   # TRUE; b was transformed to match a
#'
#' # If b had missing CRS but is actually in the same frame as a:
#' b_missing <- st_set_crs(b, NA)
#' out2 <- harmonize_crs(a, b_missing)
#' st_crs(out2$b)                   # now set to a's CRS (no transform)
#' }
#'
#' @seealso
#' [sf::st_crs()], [sf::st_transform()], [ensure_projected()]
#'
#' @family CRS helpers
harmonize_crs <- function(a, b) {
  if (!(inherits(a, c("sf","sfc")) && inherits(b, c("sf","sfc")))) {
    stop("harmonize_crs() expects sf/sfc for both arguments.")
  }
  if (is.na(sf::st_crs(a)) && !is.na(sf::st_crs(b))) sf::st_crs(a) <- sf::st_crs(b)
  if (!is.na(sf::st_crs(a)) && is.na(sf::st_crs(b))) sf::st_crs(b) <- sf::st_crs(a)
  if (!is.na(sf::st_crs(a)) && !identical(sf::st_crs(a), sf::st_crs(b))) {
    b <- sf::st_transform(b, sf::st_crs(a))
  }
  list(a = a, b = b)
}

#' Coerce arbitrary geometries to representative POINTs
#'
#' Converts any `sf` object to POINT geometries for downstream tasks like
#' seeding tessellations, spatial assignment, and plotting. Supports several
#' strategies: `"surface"` (`sf::st_point_on_surface()`), `"centroid"`
#' (`sf::st_centroid()`), and `"line_midpoint"` (a 50% sample along lines).
#' The default `"auto"` chooses a sensible strategy based on the input geometry
#' types.
#'
#' @param x An `sf` object (may contain POINT, LINESTRING, POLYGON, and multi-*
#'   variants). All non-geometry columns are preserved.
#' @param strategy One of `c("auto", "surface", "centroid", "line_midpoint")`.
#'   \describe{
#'     \item{`"auto"`}{Returns `x` unchanged if already POINT/MULTIPOINT.
#'       Uses point-on-surface for (multi)polygons and midpoints for (multi)lines.
#'       For mixed/other types, falls back to point-on-surface.}
#'     \item{`"surface"`}{Uses `sf::st_point_on_surface()`; label point lies
#'       within polygon interiors if geometries are valid.}
#'     \item{`"centroid"`}{Uses `sf::st_centroid()`; may lie outside concave
#'       polygons unless projected to a suitable CRS first.}
#'     \item{`"line_midpoint"`}{Uses `sf::st_line_sample(sample = 0.5)` to
#'       produce midpoints; collections are extracted and cast to POINT.}
#'   }
#'
#' @return An `sf` with POINT geometry:
#' \itemize{
#'   \item Same row count as `x`. If an intermediate operation would change the
#'         length (e.g., collections), a safe fallback re-computes as centroids
#'         and logs a warning.
#'   \item All original attributes retained.
#'   \item CRS preserved from the input; internal projections are temporary.
#' }
#'
#' @details
#' For numerically stable geometric operations, if `x` is in lon/lat
#' (`sf::st_is_longlat()`), the function temporarily projects the geometry
#' to a locally appropriate projected CRS via a helper (see
#' `.pick_local_projected_crs()`), performs the computation, then transforms
#' the resulting POINTs back to the original CRS. Projection failures fall back
#' to computing in the original CRS with a logged warning.
#'
#' When using `"line_midpoint"`, midpoints are computed with
#' `sf::st_line_sample(sample = 0.5)`. If the result contains MULTIPOINTs or
#' collections, they are extracted (`sf::st_collection_extract("POINT")`) and
#' cast to POINTs. Should a transform back to the original CRS fail, the
#' function falls back to centroids.
#'
#' Mixed-geometry inputs are supported. In `"auto"` mode, POINT/MULTIPOINT rows
#' are returned as-is; polygonal rows use point-on-surface; linear rows use
#' midpoints; remaining types fall back to point-on-surface.
#'
#' Diagnostic messages are emitted via `logger::log_info()` / `log_warn()` when
#' projections are applied or when fallbacks occur.
#'
#' @section Robustness notes:
#' * The function never guesses or assigns a new persistent CRS to `x`; any
#'   projection is temporary and result POINTs are returned in the original CRS.
#' * If invalid geometries cause failures in geometric operations, the function
#'   attempts a reasonable fallback (e.g., centroids) and logs a warning.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#'
#' # Polygon centroids
#' pts1 <- coerce_to_points(nc, strategy = "centroid")
#'
#' # Safer interior labels for polygons
#' pts2 <- coerce_to_points(nc, strategy = "surface")
#'
#' # Midpoints for road lines (example)
#' # roads <- st_read("roads.gpkg")
#' # pts3 <- coerce_to_points(roads, strategy = "line_midpoint")
#'
#' # Mixed input, choose automatically
#' # mix <- rbind(nc[1:3,], st_as_sf(roads[1:3,]))
#' # pts_auto <- coerce_to_points(mix, strategy = "auto")
#' }
#'
#' @seealso
#' [sf::st_point_on_surface()], [sf::st_centroid()], [sf::st_line_sample()],
#' [ensure_projected()], [harmonize_crs()]
#'
#' @family geometry conversion helpers
#' @export
coerce_to_points <- function(x, strategy = c("auto","surface","centroid","line_midpoint")) {
  strategy <- match.arg(strategy)
  if (!inherits(x, "sf")) stop("coerce_to_points() expects an sf object.")
  
  g <- sf::st_geometry(x)
  
  .project_for_calc <- function(g, parent_sf, op_name) {
    if (.is_longlat(parent_sf)) {
      trg <- .pick_local_projected_crs(parent_sf)
      log_info("coerce_to_points(): projecting long/lat to {as.character(trg)} for {op_name}")
      gp <- try(sf::st_transform(g, trg), silent = TRUE)
      if (inherits(gp, "try-error")) {
        log_warn("coerce_to_points(): projection failed for {op_name}; using original CRS")
        gp <- g
      }
      return(list(geom = gp, crs_back = sf::st_crs(g)))
    }
    list(geom = g, crs_back = NULL)
  }
  
  .on_surface <- function(g, parent_sf) {
    ctx <- .project_for_calc(g, parent_sf, "point_on_surface")
    lbl <- sf::st_point_on_surface(ctx$geom)           # compute on sfc to avoid attribute warning
    if (!is.null(ctx$crs_back) && !identical(sf::st_crs(ctx$geom), ctx$crs_back)) {
      lbl <- try(sf::st_transform(lbl, ctx$crs_back), silent = TRUE)
      if (inherits(lbl, "try-error")) lbl <- sf::st_point_on_surface(g)
    }
    lbl
  }
  
  .centroid <- function(g, parent_sf) {
    ctx <- .project_for_calc(g, parent_sf, "centroid")
    cen <- sf::st_centroid(ctx$geom)                   # compute on sfc
    if (!is.null(ctx$crs_back) && !identical(sf::st_crs(ctx$geom), ctx$crs_back)) {
      cen <- try(sf::st_transform(cen, ctx$crs_back), silent = TRUE)
      if (inherits(cen, "try-error")) cen <- sf::st_centroid(g)
    }
    cen
  }
  
  .line_midpoint <- function(g, parent_sf) {
    ctx <- .project_for_calc(g, parent_sf, "line_midpoint")
    pts <- sf::st_line_sample(ctx$geom, sample = 0.5)  # sfc POINT or MULTIPOINT
    # Avoid “already POINT” warning by checking geometry type first
    gtypes <- unique(sf::st_geometry_type(pts, by_geometry = TRUE))
    if (!all(gtypes %in% "POINT")) {
      pts <- suppressWarnings(sf::st_collection_extract(pts, "POINT"))
    }
    pts <- sf::st_cast(pts, "POINT")
    if (!is.null(ctx$crs_back) && !identical(sf::st_crs(ctx$geom), ctx$crs_back)) {
      pts <- try(sf::st_transform(pts, ctx$crs_back), silent = TRUE)
      if (inherits(pts, "try-error")) pts <- sf::st_cast(sf::st_geometry(sf::st_centroid(g)), "POINT")
    }
    pts
  }
  
  safe_bind <- function(geom) {
    if (length(geom) != nrow(x)) {
      log_warn("coerce_to_points(): geometry count mismatch ({length(geom)} vs {nrow(x)}); falling back to centroids")
      geom <- .centroid(g, x)
    }
    sf::st_sf(sf::st_drop_geometry(x), geometry = geom)
  }
  
  gcls <- unique(sf::st_geometry_type(g, by_geometry = TRUE))
  if (strategy == "auto") {
    if (all(gcls %in% c("POINT","MULTIPOINT"))) return(x)
    if (all(gcls %in% c("POLYGON","MULTIPOLYGON"))) return(safe_bind(.on_surface(g, x)))
    if (all(gcls %in% c("LINESTRING","MULTILINESTRING"))) return(safe_bind(.line_midpoint(g, x)))
    return(safe_bind(.on_surface(g, x)))
  }
  if (strategy == "surface")  return(safe_bind(.on_surface(g, x)))
  if (strategy == "centroid") return(safe_bind(.centroid(g, x)))
  if (strategy == "line_midpoint") return(safe_bind(.line_midpoint(g, x)))
  x
}


# -----------------------------------------------------------------------------
# Voronoi Tessellation
# -----------------------------------------------------------------------------

#' Create bounded Voronoi polygons from seed points
#'
#' Builds a Voronoi (Dirichlet) tessellation from seed locations and (optionally)
#' clips the result to a supplied boundary. Duplicate or coincident seeds are
#' handled deterministically via a tiny radial jitter to guarantee valid tiles.
#'
#' @param points An `sf` object whose geometry column is `POINT`/`MULTIPOINT`.
#'   Must contain at least **two distinct** coordinates after de-duplication.
#' @param clip_with Optional `sf`/`sfc` polygonal layer used to bound (clip)
#'   the Voronoi diagram. If supplied, CRS is harmonized to `points` internally.
#'
#' @return An `sf` object of polygon tiles with:
#' \itemize{
#'   \item `geometry`: POLYGON/MULTIPOLYGON Voronoi cells (clipped as needed)
#'   \item `poly_id`: sequential integer ID for each cell (1..n)
#' }
#'
#' @details
#' \strong{Preconditions}
#' * `points` must be `sf` with `POINT`/`MULTIPOINT` geometry.
#' * At least two distinct seed locations are required; otherwise an error
#'   is thrown.
#'
#' \strong{Duplicate handling}
#' * Exact duplicate coordinates are detected (at ~1e-12 precision) and are
#'   spread deterministically on a tiny circle whose radius scales with the
#'   data extent. This avoids degenerate tiles while keeping seed locations
#'   effectively unchanged at plotting scales.
#'
#' \strong{Bounding and clipping}
#' * A bounded diagram is constructed by passing a rectangle (`rw`) to
#'   \pkg{deldir}. If `clip_with` is supplied, its bbox defines the bound and
#'   the final tiles are intersected with the union of `clip_with`. Otherwise,
#'   the tiles are clipped to the bbox of `points`.
#'
#' \strong{CRS}
#' * If `clip_with` is provided, its CRS is aligned to `points` using
#'   [harmonize_crs()]. The output inherits the CRS of `points`.
#'
#' \strong{Validity}
#' * Tiles are passed through `sf::st_make_valid()` before clipping to mitigate
#'   minor topology issues from triangulation/numerics.
#'
#' @section Algorithm:
#' 1. Validate inputs and (if needed) harmonize CRS with `clip_with`.
#' 2. Extract seed coordinates and deterministically jitter duplicates.
#' 3. Compute Delaunay/Voronoi via \code{deldir::deldir()} and
#'    \code{deldir::tile.list()}.
#' 4. Build `sfc` polygons, make valid, and clip either to `clip_with`
#'    (unioned) or the points' bbox.
#' 5. Return as `sf` with a `poly_id` column.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' set.seed(1)
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' nc_state <- st_union(nc)
#' # Sample seeds inside the state boundary
#' seeds <- st_sample(nc_state, size = 40, type = "random", exact = TRUE) |> st_sf()
#' # Voronoi clipped to the state
#' vor_nc <- create_voronoi_polygons(seeds, clip_with = nc_state)
#' plot(st_geometry(vor_nc), col = NA, border = 'red')
#' plot(st_geometry(nc), add = TRUE, border = 'black')
#' }
#'
#' @seealso [harmonize_crs()], \code{\link[deldir]{deldir}},
#'   \code{\link[deldir]{tile.list}}, \code{\link[sf]{st_make_valid}},
#'   \code{\link[sf]{st_intersection}}
#'
#' @importFrom deldir deldir tile.list
#' @importFrom sf st_geometry_type st_coordinates st_bbox st_union st_as_sfc
#'   st_make_valid st_sfc st_crs st_as_sf st_set_crs st_intersection st_polygon
#'
#' @family tessellation helpers
#' @export
create_voronoi_polygons <- function(points, clip_with = NULL) {
  log_info("Creating Voronoi polygons")
  if (!inherits(points, "sf")) stop("`points` must be an sf object.")
  gt <- unique(sf::st_geometry_type(points, by_geometry = TRUE))
  if (!all(gt %in% c("POINT","MULTIPOINT"))) stop("`points` must have POINT/MULTIPOINT geometry.")
  if (nrow(points) <= 1) stop("Insufficient points for Voronoi (n <= 1).")
  
  if (!is.null(clip_with)) {
    hw <- harmonize_crs(points, clip_with)
    points <- hw$a
    clip_with <- hw$b
  }
  
  coords <- sf::st_coordinates(points)[, 1:2, drop = FALSE]
  if (nrow(unique(data.frame(x = coords[,1], y = coords[,2]))) < 2) {
    stop("Voronoi requires at least two distinct seed locations.")
  }
  
  # deterministic micro-offset for exact duplicates
  key <- interaction(round(coords[,1], 12), round(coords[,2], 12), drop = TRUE)
  tab <- table(key)
  if (any(tab > 1)) {
    rng <- max(diff(range(coords[,1])), diff(range(coords[,2])))
    if (!is.finite(rng) || rng == 0) rng <- 1
    r <- 1e-6 * rng
    dups <- names(tab)[tab > 1]
    for (gk in dups) {
      idx <- which(key == gk)
      m <- length(idx)
      ang <- seq(0, 2*pi, length.out = m + 1)[1:m]
      coords[idx,1] <- coords[idx,1] + r * cos(ang)
      coords[idx,2] <- coords[idx,2] + r * sin(ang)
    }
  }
  
  # bounded Voronoi
  bb_src <- if (!is.null(clip_with)) sf::st_bbox(clip_with) else sf::st_bbox(points)
  rw <- c(bb_src["xmin"], bb_src["xmax"], bb_src["ymin"], bb_src["ymax"])
  
  vor <- deldir::deldir(coords[,1], coords[,2], rw = rw)
  tiles <- deldir::tile.list(vor)
  
  polygons <- lapply(tiles, function(tile) {
    xy <- cbind(tile$x, tile$y)
    if (nrow(xy) < 3) return(NULL)
    if (!all(xy[1, ] == xy[nrow(xy), ])) xy <- rbind(xy, xy[1, ])
    sf::st_polygon(list(xy))
  })
  polygons <- Filter(Negate(is.null), polygons)
  
  sfc_polys <- sf::st_sfc(polygons, crs = sf::st_crs(points))
  sfc_polys <- suppressWarnings(sf::st_make_valid(sfc_polys))
  
  # clip result; intersections on sfc → sfc, so convert to sf afterwards
  sfc_clipped <- if (!is.null(clip_with)) {
    suppressWarnings(sf::st_intersection(sfc_polys, sf::st_union(clip_with)))
  } else {
    bb <- sf::st_as_sfc(sf::st_bbox(points))
    if (is.na(sf::st_crs(bb))) bb <- sf::st_set_crs(bb, sf::st_crs(points))
    suppressWarnings(sf::st_intersection(sfc_polys, bb))
  }
  
  polys_sf <- sf::st_sf(geometry = sfc_clipped)
  polys_sf$poly_id <- seq_len(nrow(polys_sf))
  polys_sf
}

# -----------------------------------------------------------------------------
# Grid Tessellations (Hex / Square)
# -----------------------------------------------------------------------------

#' Create a clipped hex or square grid targeting a cell count
#'
#' Generates a regular grid (hexagonal or square) over a polygonal boundary and
#' clips the cells to the boundary. The grid resolution is adjusted iteratively
#' so the number of output cells is close to `target_cells`.
#'
#' @param boundary An `sf`/`sfc` polygonal layer defining the area to cover.
#'   It is dissolved via `sf::st_union()` and validated with
#'   `sf::st_make_valid()` before grid generation. \strong{Use a projected CRS}
#'   (e.g., UTM or equal-area) so distances/areas are meaningful.
#' @param target_cells Approximate number of cells desired (integer > 0).
#'   The algorithm iteratively adapts the cell size until the realized count is
#'   within a tolerance band of ~15% (0.85–1.15×).
#' @param type Grid geometry to create: `"hex"` (hexagonal) or `"square"`.
#' @param max_iter Maximum number of refinement iterations for cell-size tuning.
#'
#' @return An `sf` object of grid polygons clipped to `boundary` with:
#' \itemize{
#'   \item `geometry`: POLYGON/MULTIPOLYGON cells
#'   \item `poly_id`: sequential integer identifier (1..n)
#' }
#'
#' @details
#' Cell size is initialized from \eqn{\sqrt{\mathrm{area}/\mathrm{target\_cells}}}
#' using the boundary area, then refined by comparing the realized number of
#' cells after clipping against `target_cells`. If there are too many cells,
#' the cell size is increased; if too few, decreased. Refinement stops when the
#' count falls within the tolerance or when `max_iter` is reached.
#'
#' \strong{CRS}: The function does not reproject inputs. For consistent cell
#' sizing, provide `boundary` in a suitable projected CRS. Consider calling
#' `ensure_projected()` beforehand.
#'
#' @section Algorithm:
#' 1. Dissolve and validate `boundary`; compute its area.
#' 2. Initialize cell size from area/target.
#' 3. Loop up to `max_iter`:
#'    * Create grid via `sf::st_make_grid()` (hex or square).
#'    * Clip to `boundary` (`sf::st_intersection()`), count cells.
#'    * Adjust cell size up/down based on the count ratio.
#' 4. Return clipped cells as `sf` with `poly_id`.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(dplyr)
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' state <- st_union(nc)
#' # Ensure projected CRS for sensible areas (example: Web Mercator here, but
#' # prefer a local UTM / equal-area CRS in real analyses)
#' state <- st_transform(state, 3857)
#' hex_g <- create_grid_polygons(state, target_cells = 200, type = "hex")
#' sq_g  <- create_grid_polygons(state, target_cells = 200, type = "square")
#' plot(st_geometry(hex_g), border = "red"); plot(st_geometry(state), add = TRUE)
#' }
#'
#' @seealso [create_voronoi_polygons()], [build_tessellation()], [ensure_projected()]
#'
#' @importFrom sf st_union st_make_valid st_area st_make_grid st_intersection st_as_sf
#' @importFrom dplyr mutate row_number
#'
#' @family tessellation helpers
#' @export
create_grid_polygons <- function(boundary, target_cells, type = c("hex","square"), max_iter = 8) {
  type <- match.arg(type)
  if (is.null(boundary)) stop("`boundary` must be provided for grid tessellation.")
  
  boundary <- sf::st_union(boundary) |> sf::st_make_valid()
  A <- as.numeric(sf::st_area(boundary))
  if (!is.finite(A) || A <= 0) stop("Boundary has invalid area for grid generation.")
  
  cell <- sqrt(A / max(1, target_cells))
  grid <- NULL
  for (it in seq_len(max_iter)) {
    grid <- sf::st_make_grid(boundary, cellsize = c(cell, cell), what = "polygons", square = (type == "square"))
    grid <- suppressWarnings(sf::st_intersection(grid, boundary))
    n_cells <- length(grid)
    if (n_cells == 0) { cell <- cell * 0.7; next }
    ratio <- n_cells / target_cells
    if (ratio > 1.15) {
      cell <- cell * sqrt(ratio)     # too many cells → enlarge
    } else if (ratio < 0.85) {
      cell <- cell / sqrt(1/ratio)   # too few cells → shrink
    } else break
  }
  sf::st_as_sf(grid) |>
    dplyr::mutate(poly_id = dplyr::row_number())
}

# -----------------------------------------------------------------------------
# Seeding Strategies for Voronoi
# -----------------------------------------------------------------------------

#' Generate k-means seed points from feature coordinates
#'
#' Runs k-means on the input feature coordinates and returns the cluster
#' centers as seed points (typically for Voronoi tessellation). The requested
#' `k` is clamped to the valid range \eqn{[1, n-1]} to avoid degenerate
#' solutions when there are few input features.
#'
#' @param points_sf An `sf` object whose geometry represents feature locations
#'   (preferably POINT/MULTIPOINT). Distances are computed on the raw
#'   coordinates; for sensible clustering, provide a projected CRS (e.g., UTM).
#'   Use [ensure_projected()] or [coerce_to_points()] as needed.
#' @param k Desired number of clusters (seeds). Internally constrained to
#'   `1L <= k <= (nrow(points_sf) - 1L)`.
#' @param set_seed Integer seed for reproducible k-means initialization.
#'
#' @return An `sf` POINT layer of cluster centers in the same CRS as
#'   `points_sf`. The returned object contains a single `geometry` column
#'   (no additional attributes).
#'
#' @details
#' The function extracts coordinates via `sf::st_coordinates()` and fits
#' `stats::kmeans` with `iter.max = 50` and `nstart = 10`. The resulting
#' centers are converted to an `sf` object using the CRS from `points_sf`.
#'
#' \strong{CRS}: No reprojection is performed. For meaningful spatial clusters,
#' supply `points_sf` in a projected CRS where Euclidean distance is appropriate.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' pts <- st_as_sf(data.frame(x = rnorm(200), y = rnorm(200)),
#'                 coords = c("x","y"), crs = 3857)
#' seeds <- voronoi_seeds_kmeans(pts, k = 12)
#' plot(st_geometry(pts), pch = 16, cex = 0.5)
#' plot(st_geometry(seeds), col = "red", pch = 3, cex = 1.2, add = TRUE)
#' }
#'
#' @seealso [voronoi_seeds_random()], [create_voronoi_polygons()],
#'   [build_tessellation()], [coerce_to_points()], [ensure_projected()]
#'
#' @importFrom sf st_coordinates st_as_sf st_crs
#' @importFrom stats kmeans
#'
#' @family tessellation seeding helpers
#' @export
voronoi_seeds_kmeans <- function(points_sf, k, set_seed = 456) {
  coords <- sf::st_coordinates(points_sf)
  n <- nrow(coords)
  k <- max(1L, min(k, n - 1L))
  set.seed(set_seed)
  km <- stats::kmeans(coords, centers = k, iter.max = 50, nstart = 10)
  cent <- as.data.frame(km$centers); names(cent) <- c("x","y")
  sf::st_as_sf(cent, coords = c("x","y"), crs = sf::st_crs(points_sf))
}

#' Generate random seed points within a boundary
#'
#' Draws `k` random points uniformly within the provided `boundary`
#' and returns them as an `sf` POINT layer. The boundary is dissolved
#' via `sf::st_union()` so sampling occurs across the whole area.
#'
#' @param boundary An `sf` or `sfc` object defining the sampling region.
#'   For areal sampling, supply POLYGON/MULTIPOLYGON geometry. If LINESTRING
#'   or POINT geometries are provided, sampling follows those geometries
#'   (per `sf::st_sample()` semantics).
#' @param k Integer; number of seed points to generate.
#' @param set_seed Integer seed for reproducible sampling.
#'
#' @return An `sf` POINT object of length `k` in the same CRS as `boundary`,
#'   with a single `geometry` column.
#'
#' @details
#' Sampling is performed with `sf::st_sample(type = "random", exact = TRUE)`
#' on `sf::st_union(boundary)`. For spatially uniform sampling by area,
#' pass `boundary` in a projected CRS (e.g., UTM). If `boundary` is in
#' geographic coordinates (lon/lat), the result is uniform in degrees,
#' not in meters. Use [ensure_projected()] if needed.
#'
#' The function sets the CRS on the returned points to `sf::st_crs(boundary)`.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' # Example boundary (North Carolina counties dissolved)
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' bnd <- st_union(nc)
#' # Project for area-correct sampling (optional but recommended)
#' bnd <- st_transform(bnd, 32119)  # NAD83 / NC
#'
#' seeds <- voronoi_seeds_random(bnd, k = 20, set_seed = 42)
#' plot(st_geometry(bnd), col = NA, border = 'grey40')
#' plot(st_geometry(seeds), pch = 16, col = 'red', add = TRUE)
#' }
#'
#' @seealso [voronoi_seeds_kmeans()], [create_voronoi_polygons()],
#'   [build_tessellation()], [ensure_projected()]
#'
#' @importFrom sf st_sample st_union st_set_crs st_crs st_sf
#'
#' @family tessellation seeding helpers
#' @export
voronoi_seeds_random <- function(boundary, k, set_seed = 456) {
  set.seed(set_seed)
  pts <- sf::st_sample(sf::st_union(boundary), size = k, type = "random", exact = TRUE)
  sf::st_sf(geometry = pts) |> sf::st_set_crs(sf::st_crs(boundary))
}

# -----------------------------------------------------------------------------
# Assignment: map features to tessellation polygons
# -----------------------------------------------------------------------------

#' Assign features to polygon cells
#'
#' Tags each feature in an `sf` layer with the index of the first polygon
#' it intersects from a polygon layer, returning only the features that were
#' successfully assigned. A new integer column `polygon_id` is added to the
#' returned `sf` object.
#'
#' @param features An `sf` object containing the geometries to assign
#'   (POINT/LINESTRING/POLYGON, etc.).
#' @param polygons_sf An `sf` or `sfc` object representing the target cells,
#'   typically POLYGON/MULTIPOLYGON. If an `sfc` vector is supplied it is
#'   converted to `sf` internally.
#'
#' @return An `sf` object consisting of the subset of `features` that intersect
#'   at least one polygon, with an added `polygon_id` column (integer) giving
#'   the **row index** of the first intersecting polygon in `polygons_sf`.
#'
#' @details
#' - Coordinate reference systems of the two inputs are harmonized via
#'   [harmonize_crs()] before computing intersections.
#' - Intersections are computed with `sf::st_intersects()` and the **first**
#'   hit (if multiple) is recorded. Features with no hits are dropped.
#' - `polygon_id` is based on the row position in `polygons_sf` (1..N). If your
#'   polygon layer has its own identifier column (e.g., `poly_id`), you can map
#'   indices to that identifier after the call using a join.
#' - Because `st_intersects()` is used, points lying exactly on shared
#'   boundaries may intersect multiple cells; only the first is retained. If
#'   strict containment is desired, consider replacing the logic with
#'   `sf::st_within()` in your own workflow.
#'
#' @section Performance:
#' For large layers, ensure `polygons_sf` is polygonal and in a projected CRS
#' appropriate for your region. `sf` will typically use prepared geometries to
#' accelerate `st_intersects()`, but performance still depends on geometry
#' complexity and count.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' # Polygons: dissolve NC counties to a single boundary and make a grid
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' bnd <- st_union(nc)
#' bnd <- st_transform(bnd, 32119)  # NAD83 / NC
#' grid <- st_make_grid(bnd, n = c(8, 6), what = "polygons")
#' grid <- st_intersection(grid, bnd)
#' grid_sf <- st_as_sf(grid)
#'
#' # Sample points inside the boundary
#' set.seed(42)
#' pts <- st_sample(bnd, size = 500, type = "random", exact = TRUE) |> st_sf()
#'
#' # Assign points to grid cells
#' assigned <- assign_features_to_polygons(pts, grid_sf)
#' head(assigned$polygon_id)
#'
#' # If grid has an ID column and you want that instead of row index:
#' grid_sf$cell_id <- seq_len(nrow(grid_sf))
#' assigned$cell_id <- grid_sf$cell_id[assigned$polygon_id]
#' }
#'
#' @seealso [harmonize_crs()], [coerce_to_points()],
#'   [create_voronoi_polygons()], [create_grid_polygons()], [build_tessellation()]
#'
#' @importFrom sf st_as_sf st_geometry st_intersects
#'
#' @family tessellation utilities
#' @export
assign_features_to_polygons <- function(features, polygons_sf) {
  if (!inherits(features, "sf")) stop("`features` must be an sf object.")
  if (!(inherits(polygons_sf, "sf") || inherits(polygons_sf, "sfc"))) {
    stop("`polygons_sf` must be sf/sfc POLYGON/MULTIPOLYGON.")
  }
  if (inherits(polygons_sf, "sfc")) polygons_sf <- sf::st_as_sf(polygons_sf)
  
  hh <- harmonize_crs(features, polygons_sf)
  features <- hh$a
  polygons_sf <- hh$b
  
  hits <- sf::st_intersects(features, sf::st_geometry(polygons_sf))
  features$polygon_id <- vapply(hits, function(ix) if (length(ix)) ix[1] else NA_integer_, 1L)
  features <- features[!is.na(features$polygon_id), , drop = FALSE]
  features
}

# -----------------------------------------------------------------------------
# Level Selection (cluster counts / target cells)
# -----------------------------------------------------------------------------

#' Choose reasonable tessellation levels via an elbow heuristic
#'
#' Estimates a small set of candidate cluster counts (number of polygons/cells)
#' by running k-means over feature coordinates and detecting the "elbow" in the
#' total within-cluster sum of squares (WSS) curve. Returns up to `top_n`
#' levels centered on the detected elbow (elbow − 1, elbow, elbow + 1).
#'
#' @param points_sf An `sf` object with one coordinate per row (typically POINT
#'   geometries). If your data are lines/polygons or otherwise not one point per
#'   row, first convert with [coerce_to_points()].
#' @param max_levels Integer, maximum number of clusters to consider. The search
#'   actually evaluates `k = 1..max_k` where `max_k = min(max_levels, n - 1)`,
#'   with `n` the number of rows in `points_sf`.
#' @param top_n Integer, number of levels to return (up to three are available
#'   around the elbow). Defaults to `3`.
#'
#' @return An integer vector of candidate levels (cluster counts), ordered in
#'   descending order, of length `<= top_n`. If fewer than three levels are
#'   available (e.g., very small `n`), fewer values are returned. If `n < 3`,
#'   the function returns `1L`.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Extracts coordinates with `sf::st_coordinates(points_sf)`.
#'   \item Computes k-means (Euclidean) for `k = 1..max_k` and records WSS.
#'   \item Uses a discrete second-derivative heuristic (`which.min(diff(diff(WSS))) + 2`)
#'         to pick the elbow.
#'   \item Returns up to three levels: elbow − 1, elbow, elbow + 1 (bounded to
#'         `[1, max_k]`), then keeps the top `top_n`.
#' }
#'
#' \strong{CRS note:} k-means uses Euclidean distances in the native coordinate
#' units. For geographic (lon/lat) coordinates, consider projecting to a local
#' planar CRS (see [ensure_projected()]) to avoid large distortion.
#'
#' @section Robustness:
#' The routine guards against tiny samples and caps `k` by `n - 1` to avoid
#' degenerate partitions.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' # Make 200 random points inside North Carolina (sf example data)
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' bnd <- st_union(nc)
#' set.seed(1)
#' pts <- st_sample(bnd, size = 200, type = "random", exact = TRUE) |> st_sf()
#' pts <- st_transform(pts, 32119)  # project for sensible Euclidean distances
#'
#' determine_optimal_levels(pts, max_levels = 12, top_n = 3)
#' }
#'
#' @seealso [coerce_to_points()], [build_tessellation()], [ensure_projected()]
#'
#' @importFrom sf st_coordinates
#' @importFrom stats kmeans
#'
#' @family tessellation utilities
#' @export
determine_optimal_levels <- function(points_sf, max_levels = 10, top_n = 3) {
  coords <- sf::st_coordinates(points_sf)
  n <- nrow(coords)
  if (n < 3) return(1L)
  max_k <- max(2L, min(max_levels, n - 1L))
  wss <- numeric(max_k)
  for (k in 1:max_k) {
    wss[k] <- stats::kmeans(coords, centers = k, iter.max = 20, nstart = 5)$tot.withinss
  }
  elbow_k <- if (length(wss) >= 4) which.min(diff(diff(wss))) + 2L else which.min(wss)
  picks <- sort(unique(pmax(1L, pmin(max_k, c(elbow_k - 1L, elbow_k, elbow_k + 1L)))), decreasing = TRUE)
  head(picks, min(top_n, length(picks)))
}

# -----------------------------------------------------------------------------
# Modeling: GWR and Bayesian Spatial (spBayes)
# -----------------------------------------------------------------------------

#' Validate presence of modeling variables in an input dataset
#'
#' Confirms that the response and all predictor variables exist as columns in
#' `data_sf`. If any are missing, an informative error is thrown listing the
#' absent variable names. When all variables are present, the function returns
#' (invisibly) `TRUE`.
#'
#' @param data_sf An `sf` object (or `data.frame`) containing the modeling
#'   variables as columns.
#' @param response_var A single string naming the response variable column.
#' @param predictor_vars A character vector of predictor (feature) column names.
#'
#' @return Invisibly returns `TRUE` when all requested variables are present.
#'   Otherwise, the function stops with an error describing which variables are
#'   missing.
#'
#' @details This helper performs a lightweight schema check only; it does not
#'   modify data or validate types. It is intended to be called by higher-level
#'   modeling routines prior to fitting (e.g., to guard early and produce clear
#'   error messages).
#'
#' @examples
#' \dontrun{
#' library(sf)
#' d <- data.frame(y = rnorm(10), x1 = runif(10), x2 = runif(10))
#' d <- st_as_sf(d, coords = c("x1", "x2"), remove = FALSE)
#' .validate_model_inputs(d, "y", c("x1", "x2"))  # returns TRUE (invisibly)
#'
#' # This will error because x3 is missing:
#' .validate_model_inputs(d, "y", c("x1", "x3"))
#' }
#'
#' @seealso [fit_gwr_model()], [fit_bayesian_spatial_model()]
#'
#' @family modeling utilities
#' @keywords internal
#' @noRd
.validate_model_inputs <- function(data_sf, response_var, predictor_vars) {
  miss <- setdiff(c(response_var, predictor_vars), names(data_sf))
  if (length(miss)) stop(sprintf("Missing variables in data: %s", paste(miss, collapse = ", ")))
  invisible(TRUE)
}

#' Fit a Geographically Weighted Regression (GWR) with adaptive bandwidth
#'
#' Wraps \pkg{spgwr} to fit a GWR model using an adaptively selected kernel
#' bandwidth (as a proportion of the sample) via [spgwr::gwr.sel()]. The
#' function safeguards non-finite selections by defaulting to `0.75` and caps
#' values at `< 1` (set to `0.99`) to avoid the degenerate full-neighborhood
#' case. The fitted object includes a hat matrix so that AICc is available.
#'
#' @param data_sf An `sf` object containing the response and predictor columns
#'   and a valid geometry column. For stable distance calculations, prefer a
#'   projected CRS in linear units (e.g., meters). Upstream helpers such as
#'   [ensure_projected()] can be used to enforce this.
#' @param response_var A single string naming the response variable column.
#' @param predictor_vars A character vector of predictor (feature) column names.
#'
#' @return A named list with components:
#' \describe{
#'   \item{`model`}{The fitted `spgwr::gwr` object.}
#'   \item{`bandwidth`}{Numeric scalar in `(0, 1]` giving the adaptive
#'         bandwidth proportion actually used.}
#'   \item{`AICc`}{Numeric scalar: corrected Akaike Information Criterion
#'         extracted from the fit (requires `hatmatrix = TRUE`).}
#' }
#'
#' @details
#' Internally, `data_sf` is coerced to `sp` classes (`Spatial*`) because
#' \pkg{spgwr} operates on `sp` rather than `sf`. The adaptive bandwidth is
#' chosen by [spgwr::gwr.sel()] with `adapt = TRUE`. If the selector returns a
#' non-finite result, a fallback of `0.75` is used. If it equals or exceeds `1`,
#' the value is reduced to `0.99` to avoid numerically unstable behavior.
#'
#' Distances are interpreted in the units of the data’s CRS. Using an
#' appropriate projected CRS (e.g., UTM) is recommended; see
#' [ensure_projected()] and [harmonize_crs()] in this toolkit.
#'
#' @section Logging:
#' Uses \pkg{logger} to emit progress messages (e.g., "Fitting GWR model").
#'
#' @examples
#' \dontrun{
#' library(sf)
#' # Toy data
#' set.seed(1)
#' df <- data.frame(y = rnorm(50), x1 = runif(50), x2 = rnorm(50),
#'                  lon = runif(50, -96, -95), lat = runif(50, 29.5, 30.2))
#' sf_pts <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
#' sf_pts <- ensure_projected(sf_pts)  # enforce projected CRS for distance
#'
#' res <- fit_gwr_model(sf_pts, response_var = "y", predictor_vars = c("x1","x2"))
#' res$bandwidth
#' res$AICc
#' }
#'
#' @seealso [fit_bayesian_spatial_model()], [determine_optimal_levels()],
#'   [ensure_projected()], [spgwr::gwr()], [spgwr::gwr.sel()]
#'
#' @importFrom spgwr gwr gwr.sel
#' @family modeling helpers
#' @export
fit_gwr_model <- function(data_sf, response_var, predictor_vars) {
  log_info("Fitting GWR model")
  .validate_model_inputs(data_sf, response_var, predictor_vars)
  
  data_sp <- as(data_sf, "Spatial")
  fml <- as.formula(paste(response_var, "~", paste(predictor_vars, collapse = " + ")))
  bw <- suppressWarnings(spgwr::gwr.sel(fml, data = data_sp, adapt = TRUE))
  # bw is a proportion (0,1] for adaptive kernels
  if (!is.finite(bw) || is.na(bw)) bw <- 0.75
  if (isTRUE(all.equal(bw, 1)) || bw >= 1) bw <- 0.99
  
  # Use adaptive bandwidth by passing it as the 'adapt' proportion
  fit <- spgwr::gwr(fml, data = data_sp, adapt = bw, hatmatrix = TRUE)
  list(model = fit, bandwidth = bw, AICc = fit$results$AICc)
}

#' Heuristic prior bounds for spatial decay parameter \eqn{\phi}
#'
#' Computes lower/upper bounds for the spatial decay/range parameter
#' \eqn{\phi} used in covariance models (e.g., exponential, Matérn) based on
#' pairwise interpoint distances. The rule-of-thumb links \eqn{\phi} to a
#' “practical range” via \eqn{\phi \approx 3 / r}. Here we set
#' \code{lower = 3 / max(distance)} and \code{upper = max(lower * 1.2, 3 / dq)},
#' where \code{dq} is a lower quantile (default 10\%) of nonzero distances.
#' Falls back to \code{c(0.001, 1)} if distances are degenerate.
#'
#' @param coords_xy A numeric matrix/data frame of coordinates \code{(n x d)}
#'   in **linear units** (e.g., meters). Use a projected CRS upstream so that
#'   pairwise distances are meaningful.
#' @param q_small Numeric scalar in \code{(0, 1)} giving the small distance
#'   quantile used for the upper bound (default \code{0.1} for the 10th
#'   percentile of positive distances).
#'
#' @return A numeric length-2 vector \code{c(lower, upper)} providing
#'   heuristic bounds for \eqn{\phi}.
#'
#' @details
#' The bounds aim to be weakly informative yet scale-aware:
#' \itemize{
#'   \item \strong{Lower}: prevents excessively long ranges
#'         (\eqn{\phi} too small) by tying to the domain diameter.
#'   \item \strong{Upper}: prevents overly short ranges (\eqn{\phi} too large)
#'         by tying to a small—but nonzero—interpoint spacing (quantile).
#' }
#' If all pairwise distances are zero or missing, the function returns
#' \code{c(0.001, 1)} as a conservative fallback.
#'
#' @section Units:
#' Distances inherit the units of \code{coords_xy}. Always project geographic
#' coordinates before calling (see \code{ensure_projected()}).
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' xy <- cbind(runif(100, 0, 10000), runif(100, 0, 8000))  # meters
#' .phi_prior_bounds(xy)          # default (10% quantile)
#' .phi_prior_bounds(xy, 0.05)    # tighter upper bound via 5% quantile
#' }
#'
#' @seealso \code{fit_bayesian_spatial_model()}
#'
#' @importFrom stats dist quantile
#' @keywords internal
.phi_prior_bounds <- function(coords_xy, q_small = 0.1) {
  d <- as.numeric(dist(coords_xy))
  d <- d[d > 0]
  if (!length(d)) return(c(0.001, 1))
  dmax <- max(d)
  dq <- stats::quantile(d, probs = q_small, names = FALSE, type = 7)
  lower <- 3 / dmax
  upper <- max(lower * 1.2, 3 / dq)
  c(as.numeric(lower), as.numeric(upper))
}

#' Matérn correlation function \eqn{\rho(h;\,\phi,\nu)}
#'
#' Computes the Matérn correlation for pairwise distances \code{h} with
#' decay (range) parameter \eqn{\phi>0} and smoothness \eqn{\nu>0}:
#' \deqn{\rho(h) \;=\; \frac{2^{\,1-\nu}}{\Gamma(\nu)} \, (\phi h)^{\nu}\, K_{\nu}(\phi h),}
#' where \eqn{K_{\nu}} is the modified Bessel function of the second kind.
#' The diagonal (zero distance) is set to 1 explicitly, and tiny negative
#' numerical artifacts are truncated at 0.
#'
#' @param h A numeric vector/matrix of distances (usually an \eqn{n\times n}
#'   pairwise distance matrix) in **linear units** (e.g., meters). If a vector
#'   is supplied it is coerced to a 1-row matrix.
#' @param phi Positive decay/range parameter. Larger \eqn{\phi} yields faster
#'   correlation decay; a common rule-of-thumb is a practical range
#'   \eqn{r \approx 3/\phi}.
#' @param nu Positive smoothness parameter controlling differentiability of the
#'   underlying process (\eqn{\nu = 0.5} gives the exponential model; larger
#'   \eqn{\nu} approaches Gaussian-like smoothness).
#'
#' @return A numeric matrix the same dimension as \code{as.matrix(h)} with
#'   entries in \eqn{[0,1]} and ones on the diagonal.
#'
#' @details
#' For numerical stability, entries with \eqn{h=0} are set to 1 and any
#' negative values from floating-point error are clamped to 0. Ensure
#' geographic coordinates are projected upstream so \code{h} represents
#' meaningful Euclidean distances (see \code{ensure_projected()}).
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' xy <- cbind(runif(20, 0, 10000), runif(20, 0, 8000))  # meters
#' D  <- as.matrix(dist(xy))
#' R  <- .matern_rho(D, phi = 0.0005, nu = 1.5)
#' range(R)  # correlations in [0,1]
#' }
#'
#' @references
#' Stein, M. L. (1999) *Interpolation of Spatial Data: Some Theory for Kriging*.
#' Rasmussen, C. E., & Williams, C. K. I. (2006) *Gaussian Processes for Machine Learning*.
#'
#' @seealso \code{fit_bayesian_spatial_model()}, \code{.phi_prior_bounds()}
#' @keywords internal
.matern_rho <- function(h, phi, nu) {
  h <- as.matrix(h)
  z <- phi * h
  rho <- matrix(0, nrow(h), ncol(h))
  diag(rho) <- 1
  mask <- z > 0
  if (any(mask)) {
    zpos <- z[mask]
    fac <- (2^(1 - nu)) / gamma(nu)
    rho[mask] <- fac * (zpos^nu) * besselK(zpos, nu)
    rho[mask] <- pmax(rho[mask], 0)
  }
  rho
}

#' Stable multivariate normal log-likelihood with diagonal jitter
#'
#' Evaluates \eqn{\log p(y \mid \mu, V)} under a multivariate Normal model
#' using \code{mvtnorm::dmvnorm}. If the evaluation fails (typically because
#' \code{V} is nearly singular or not numerically positive definite), a tiny
#' diagonal jitter \eqn{\epsilon I} with
#' \eqn{\epsilon = 10^{-8}\,\mathrm{mean}(\mathrm{diag}(V))} is added and the
#' likelihood is recomputed.
#'
#' @param y Numeric vector of observations of length \eqn{n}.
#' @param mu Numeric mean vector of length \eqn{n}.
#' @param V Numeric \eqn{n \times n} covariance matrix (symmetric, ideally
#'   positive definite) corresponding to \code{y}.
#'
#' @return A single numeric value: the log-likelihood \eqn{\log p(y \mid \mu, V)}.
#'
#' @details
#' This helper is designed for robustness in MCMC or model-selection loops where
#' covariance matrices may become ill-conditioned due to rounding, duplicated
#' locations, or extreme hyperparameters. The jitter is applied only if the
#' first evaluation errors; it is not a substitute for proper model
#' regularization. If jittering is triggered frequently, revisit the covariance
#' construction (e.g., add a nugget or reproject coordinates).
#'
#' The input \code{V} is modified only locally (a copy with jitter is used if
#' needed). Ensure that \code{y}, \code{mu}, and \code{V} are dimensionally
#' consistent; no recycling is performed.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n  <- 3
#' y  <- c(0.2, -0.1, 0.3)
#' mu <- rep(0, n)
#' # Nearly singular covariance (high correlation)
#' V  <- matrix(0.99, n, n); diag(V) <- 1
#' # Evaluate robustly:
#' ll <- .safe_ll(y, mu, V)
#' ll
#' }
#'
#' @seealso \code{mvtnorm::dmvnorm}, \code{fit_bayesian_spatial_model()}
#' @keywords internal
.safe_ll <- function(y, mu, V) {
  out <- try(mvtnorm::dmvnorm(y, mean = mu, sigma = V, log = TRUE), silent = TRUE)
  if (inherits(out, "try-error")) {
    diag(V) <- diag(V) + 1e-8 * mean(diag(V))
    out <- mvtnorm::dmvnorm(y, mean = mu, sigma = V, log = TRUE)
  }
  out
}

#' Bayesian spatial regression via spBayes with DIC-like proxy
#'
#' Fits a Bayesian spatial linear model using \pkg{spBayes} (\code{spLM}) on
#' point-referenced data stored in an \pkg{sf} object and returns the fitted
#' model, recovered posterior samples, and a lightweight DIC-style score
#' computed from posterior predictive log-likelihoods.
#'
#' @param data_sf An \code{sf} object containing the response, predictors, and
#'   point geometries. Coordinates are extracted from \code{st_geometry(data_sf)}.
#'   For meaningful spatial scales, prefer a projected CRS (e.g., meters).
#' @param response_var String; the name of the response variable column in
#'   \code{data_sf}.
#' @param predictor_vars Character vector; names of predictor columns in
#'   \code{data_sf}. An intercept is included automatically via the model
#'   formula \code{response_var ~ predictor_vars}.
#' @param n.samples Integer; number of MCMC samples to draw in \code{spBayes::spLM}.
#'   The function uses a fixed 50\% burn-in when calling \code{spBayes::spRecover}.
#'   Default is \code{5000}.
#' @param cov_model One of \code{"exponential"}, \code{"spherical"}, or
#'   \code{"matern"}; selects the spatial covariance family used by \code{spLM}.
#'
#' @details
#' The design matrix \eqn{X} is built from \code{predictor_vars}; coordinates are
#' taken from the geometry of \code{data_sf}. Prior bounds for the decay
#' parameter \eqn{\phi} are derived from interpoint distances using
#' \code{.phi_prior_bounds()}, which sets a weakly-informative uniform prior over a
#' scale consistent with the geometry units. For the Matérn model, a uniform
#' prior for smoothness \eqn{\nu \sim \mathcal{U}(0.5, 2.5)} is used with an initial
#' value of 1.5.
#'
#' Starting values and tuning parameters are modestly informative:
#' \itemize{
#'   \item \code{starting = list(beta = 0, phi = mean(phi_bounds), sigma.sq = 1, tau.sq = 1 [, nu = 1.5])}
#'   \item \code{tuning = list(phi = 0.1, sigma.sq = 0.1, tau.sq = 0.1 [, nu = 0.05])}
#' }
#' After running \code{spLM}, the function calls \code{spRecover} with a
#' burn-in of \eqn{\lfloor 0.5\,n.samples \rfloor}.
#'
#' A DIC-like proxy is computed to aid model/level comparison:
#' from a subset of posterior parameter draws, covariance matrices \eqn{V} are
#' formed, the posterior mean \eqn{\bar\beta} is used to obtain \eqn{\mu = X\bar\beta},
#' and log-likelihoods \eqn{\log p(y \mid \mu, V)} are evaluated via
#' \code{mvtnorm::dmvnorm}. If evaluation fails due to near-singular \eqn{V},
#' \code{.safe_ll()} adds a tiny diagonal jitter. The proxy returned is
#' \deqn{\mathrm{DIC}_{\mathrm{proxy}} = -2\,\mathrm{mean}(\ell) + \mathrm{var}(\ell).}
#' This is a heuristic (not Spiegelhalter's formal DIC) intended for relative
#' model selection across tessellation levels/covariance families.
#'
#' \strong{Units and CRS:} The decay parameter \eqn{\phi} and distance matrix
#' \eqn{D} depend on coordinate units. If your data are in longitude/latitude,
#' reproject to a suitable local projected CRS (e.g., UTM) beforehand to ensure
#' interpretable length scales and stable numerics.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{model}: the \code{spBayes::spLM} fitted model object.
#'   \item \code{samples}: the object returned by \code{spBayes::spRecover}
#'     (posterior samples of regression coefficients and, when applicable,
#'     covariance parameters).
#'   \item \code{DIC}: numeric scalar, the DIC-like proxy described above
#'     (smaller is better, for relative comparison).
#'   \item \code{phi_prior}: length-2 numeric vector with the \eqn{\phi} prior
#'     bounds used (\code{c(lower, upper)}).
#' }
#'
#' @section Errors and diagnostics:
#' The function validates that \code{response_var} and \code{predictor_vars} exist
#' in \code{data_sf}. If \code{data_sf} lacks geometry or contains duplicate/very
#' close points leading to ill-conditioning, the likelihood computation may rely
#' on jitter via \code{.safe_ll()}. Frequent jittering suggests reprojecting,
#' deduplicating, adding a nugget, or revisiting priors.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' set.seed(1)
#' pts <- st_as_sf(data.frame(x = runif(80), y = runif(80)),
#'                 coords = c("x","y"), crs = 3857)
#' pts$X1 <- rnorm(nrow(pts))
#' pts$X2 <- rnorm(nrow(pts))
#' pts$Y  <- 1 + 2*pts$X1 - 1*pts$X2 + rnorm(nrow(pts), sd = 0.3)
#'
#' fit <- fit_bayesian_spatial_model(
#'   data_sf = pts,
#'   response_var = "Y",
#'   predictor_vars = c("X1","X2"),
#'   n.samples = 3000,
#'   cov_model = "matern"
#' )
#' fit$DIC
#' }
#'
#' @seealso
#' \code{\link[spBayes]{spLM}}, \code{\link[spBayes]{spRecover}},
#' \code{\link[mvtnorm]{dmvnorm}}, \code{.phi_prior_bounds()}, \code{.matern_rho()},
#' \code{.safe_ll()}
#'
#' @export
fit_bayesian_spatial_model <- function(data_sf, response_var, predictor_vars, n.samples = 5000,
                                       cov_model = c("exponential","spherical","matern")) {
  cov_model <- match.arg(cov_model)
  log_info(sprintf("Fitting Bayesian spatial model (%s); n.samples=%d", cov_model, n.samples))
  .validate_model_inputs(data_sf, response_var, predictor_vars)
  
  df <- sf::st_drop_geometry(data_sf)
  xy <- sf::st_coordinates(sf::st_geometry(data_sf))
  coords_xy <- as.matrix(xy)
  
  fml <- as.formula(paste(response_var, "~", paste(predictor_vars, collapse = " + ")))
  X <- model.matrix(fml, df); y <- df[[response_var]]; p <- ncol(X)
  
  phi_bounds <- .phi_prior_bounds(coords_xy)
  priors <- list(
    beta.Norm   = list(rep(0, p), diag(100, p)),
    phi.Unif    = phi_bounds,
    sigma.sq.IG = c(2, 1),
    tau.sq.IG   = c(2, 1)
  )
  starting <- list(beta = rep(0, p), phi = mean(phi_bounds), sigma.sq = 1, tau.sq = 1)
  tuning   <- list(phi = 0.1, sigma.sq = 0.1, tau.sq = 0.1)
  
  if (cov_model == "matern") {
    priors$nu.Unif <- c(0.5, 2.5)
    starting$nu <- 1.5
    tuning$nu <- 0.05
  }
  
  sm <- spBayes::spLM(
    fml, data = df, coords = coords_xy,
    starting = starting, tuning = tuning, priors = priors,
    cov.model = cov_model, n.samples = n.samples, verbose = TRUE
  )
  
  burn <- floor(0.5 * n.samples)
  rec  <- spBayes::spRecover(sm, start = burn)
  
  theta <- as.matrix(sm$p.theta.samples)                 # sigma.sq, tau.sq, phi, (nu if matern)
  betaS <- as.matrix(rec$p.beta.recover.samples)
  beta_bar <- colMeans(betaS)
  
  D <- as.matrix(stats::dist(coords_xy))
  mu <- as.vector(X %*% beta_bar)
  
  .cov_from_theta <- function(ssq, tsq, phi, D) {
    ssq * exp(-phi * D) + diag(tsq, nrow(D))
  }
  if (cov_model == "spherical") {
    .cov_from_theta <- function(ssq, tsq, phi, D) {
      a <- 1 / max(phi, .Machine$double.eps)
      h <- D
      rho <- matrix(0, nrow(h), ncol(h))
      mask <- h < a + .Machine$double.eps
      ha <- (h[mask] / a)
      rho[mask] <- 1 - 1.5*ha + 0.5*ha^3
      diag(rho) <- 1
      ssq * rho + diag(tsq, nrow(h))
    }
  }
  if (cov_model == "matern") {
    .cov_from_theta <- function(ssq, tsq, phi, D, nu) {
      rho <- .matern_rho(D, phi = phi, nu = nu)
      ssq * rho + diag(tsq, nrow(D))
    }
  }
  
  idx <- seq_len(nrow(theta))
  if (length(idx) > 300) idx <- unique(round(seq(1, length(idx), length.out = 300)))
  
  ll_vals <- vapply(idx, function(i) {
    ssq <- theta[i, "sigma.sq"]; tsq <- theta[i, "tau.sq"]; phi <- theta[i, "phi"]
    if (cov_model == "matern") {
      this_nu <- if ("nu" %in% colnames(theta)) theta[i, "nu"] else 1.5
      V <- .cov_from_theta(ssq, tsq, phi, D, this_nu)
    } else {
      V <- .cov_from_theta(ssq, tsq, phi, D)
    }
    .safe_ll(y, mu, V)
  }, numeric(1))
  
  dic_proxy <- -2 * mean(ll_vals) + stats::var(ll_vals)
  list(model = sm, samples = rec, DIC = dic_proxy, phi_prior = phi_bounds)
}

# -----------------------------------------------------------------------------
# Tessellation Pipeline (Voronoi / Hex / Square)
# -----------------------------------------------------------------------------

#' Build tessellations and assign features to cells
#'
#' Constructs a tessellation over a study area using one of three methods
#' (Voronoi, hexagonal grid, or square grid) at one or more target levels, and
#' assigns input features to the resulting polygons. Robustly harmonizes CRS,
#' pointizes non-point geometries, and (optionally) clips to a supplied boundary.
#'
#' @param features_sf An \code{sf} object containing the features to be assigned
#'   to tessellation cells. May contain POINT/LINE/POLYGON geometries; non-point
#'   geometries are converted to representative points via \code{coerce_to_points()}.
#' @param levels Integer or numeric vector; target numbers of cells (for grids)
#'   or seeds (for Voronoi). Each value is processed separately; the return is a
#'   named list keyed by \code{as.character(level)}.
#' @param method One of \code{"voronoi"}, \code{"hex"}, or \code{"square"}.
#'   \itemize{
#'     \item \strong{voronoi}: Seeds are generated and Voronoi tiles are clipped
#'       to \code{boundary} (if provided) or the features' bounding box.
#'     \item \strong{hex}/\strong{square}: A grid is generated over \code{boundary}
#'       (required for good behavior; if omitted, the features' bounding box is used).
#'   }
#' @param boundary Optional \code{sf} or \code{sfc} POLYGON/MULTIPOLYGON used to
#'   constrain/clip the tessellation. If provided, CRSs are harmonized and
#'   \code{features_sf} is projected to match. If omitted:
#'   \itemize{
#'     \item Voronoi: tiles are clipped to the features' bounding box.
#'     \item Grids: a bounding-box polygon is derived from the features.
#'   }
#' @param seeds Seeding strategy for Voronoi; one of \code{"kmeans"},
#'   \code{"random"}, or \code{"provided"}. Ignored for hex/square grids.
#'   \itemize{
#'     \item \code{"kmeans"}: cluster the (pointized) features into \code{k = level}
#'       groups and use cluster centroids as seeds.
#'     \item \code{"random"}: sample \code{k = level} points uniformly within
#'       \code{boundary} (or the features' bbox if boundary is missing).
#'     \item \code{"provided"}: use \code{provided_seed_points} as-is (see below).
#'   }
#' @param provided_seed_points Optional \code{sf} POINT/MULTIPOINT object used
#'   when \code{seeds = "provided"}. If the count differs from \code{level},
#'   the function proceeds with the provided count and logs a warning.
#' @param pointize Strategy for converting non-point geometries to points prior
#'   to seeding/assignment; one of \code{"auto"}, \code{"surface"},
#'   \code{"centroid"}, or \code{"line_midpoint"}. See
#'   \code{coerce_to_points()} for details. Default \code{"auto"} chooses a
#'   sensible method by geometry type.
#' @param seed_kmeans_from Deprecated compatibility parameter; currently treated
#'   as \code{"features"} and ignored.
#'
#' @details
#' \strong{CRS handling:} All inputs are passed through \code{ensure_projected()}
#' to avoid long/lat distance artifacts. When \code{boundary} is supplied,
#' \code{features_sf} is projected to match \code{boundary}; otherwise a local
#' projected CRS is chosen heuristically.
#'
#' \strong{Voronoi path:} Seeds are prepared per \code{seeds}. Voronoi polygons
#' are built via \code{create_voronoi_polygons()} and clipped to \code{boundary}
#' (if provided) or a bounding box. Exact duplicate seeds receive a deterministic
#' micro-offset to ensure valid tessellation. If fewer than two seeds remain, a
#' single polygon equal to the clip target is returned.
#'
#' \strong{Grid path:} Hex/square grids are generated with
#' \code{create_grid_polygons()}, adjusting cell size iteratively to approach
#' \code{target_cells = level}.
#'
#' \strong{Assignment:} Features (pointized if needed) are assigned to polygons
#' using \code{st_intersects()} via \code{assign_features_to_polygons()}, adding a
#' \code{polygon_id} column. Unassigned features are dropped.
#'
#' @return A named \code{list} where each element corresponds to a value in
#'   \code{levels} (name \code{as.character(level)}). Each element is a \code{list}
#'   with components:
#'   \itemize{
#'     \item \code{polygons}: \code{sf} POLYGON/MULTIPOLYGON with column
#'       \code{poly_id} (1..n).
#'     \item \code{data}: \code{sf} POINT layer of assigned features with
#'       \code{polygon_id} referencing \code{poly_id}.
#'     \item \code{method}: the tessellation method used.
#'     \item \code{level}: the level (target K / target cells).
#'   }
#'
#' @section Logs & warnings:
#' The function emits informative log messages (via \pkg{logger}) about CRS
#' choices, boundary coercion, seed counts, and grid sizing. If
#' \code{seeds = "provided"} and \code{nrow(provided_seed_points) != level}, a
#' warning is logged and the provided seeds are used.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' set.seed(123)
#' # Synthetic points within a square AOI
#' aoi <- st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 1000, ymax = 1000)), crs = 3857)
#' pts <- st_sample(aoi, 300) |> st_sf() |> st_set_crs(3857)
#'
#' # Voronoi (kmeans seeds), k = 8, clipped to AOI
#' vt <- build_tessellation(pts, levels = 8, method = "voronoi", boundary = st_sf(geometry = aoi))
#' plot(st_geometry(vt[["8"]]$polygons)); points(st_coordinates(vt[["8"]]$data))
#'
#' # Hex grid with ~20 cells
#' ht <- build_tessellation(pts, levels = 20, method = "hex", boundary = st_sf(geometry = aoi))
#' }
#'
#' @seealso
#' \code{\link{coerce_to_points}}, \code{\link{ensure_projected}},
#' \code{\link{harmonize_crs}}, \code{\link{create_voronoi_polygons}},
#' \code{\link{create_grid_polygons}}, \code{\link{assign_features_to_polygons}},
#' \code{\link{voronoi_seeds_kmeans}}, \code{\link{voronoi_seeds_random}}
#'
#' @export
build_tessellation <- function(
    features_sf,
    levels,
    method = c("voronoi","hex","square"),
    boundary = NULL,
    seeds = c("kmeans","random","provided"),
    provided_seed_points = NULL,
    pointize = c("auto","surface","centroid","line_midpoint"),
    seed_kmeans_from = c("features","assigned_points")  # kept for API compatibility; acts like "features"
) {
  method <- match.arg(method)
  seeds  <- match.arg(seeds)
  pointize <- match.arg(pointize)
  seed_kmeans_from <- match.arg(seed_kmeans_from)
  
  # Prepare working CRS and points for seeding/assignment
  if (!is.null(boundary)) {
    boundary <- ensure_projected(boundary)
    # DEBUG: boundary class before coercion
    log_info("build_tessellation(): initial boundary class: {paste(class(boundary), collapse=',')}")
    # Ensure boundary is sf (not just sfc) for dplyr ops
    if (inherits(boundary, "sfc")) boundary <- sf::st_sf(geometry = boundary)
    features_sf <- ensure_projected(features_sf, boundary)
  } else {
    features_sf <- ensure_projected(features_sf)
  }
  pts_for_seeds <- coerce_to_points(features_sf, strategy = pointize)
  
  # Validate provided seeds if used
  if (seeds == "provided") {
    if (is.null(provided_seed_points)) stop("Seeds set to 'provided' but `provided_seed_points` is NULL.")
    if (!inherits(provided_seed_points, "sf")) stop("`provided_seed_points` must be an sf object.")
    gtp <- unique(sf::st_geometry_type(provided_seed_points, by_geometry = TRUE))
    if (!all(gtp %in% c("POINT","MULTIPOINT"))) stop("`provided_seed_points` must have POINT/MULTIPOINT geometry.")
  }
  
  res <- list()
  for (lvl in levels) {
    log_info("build_tessellation(): method={method}, level={lvl}, seeds={seeds}")
    
    if (method == "voronoi") {
      if (seeds == "provided") {
        sp <- provided_seed_points
        hh <- harmonize_crs(pts_for_seeds, sp)
        pts_for_seeds <- hh$a; sp <- hh$b
        if (nrow(sp) != lvl) {
          log_warn("Provided {nrow(sp)} seeds but level={lvl}; proceeding with provided seed count.")
        }
      } else if (seeds == "kmeans") {
        sp <- voronoi_seeds_kmeans(pts_for_seeds, k = lvl)
      } else { # random
        if (is.null(boundary)) {
          b <- sf::st_as_sfc(sf::st_bbox(pts_for_seeds))
          b <- sf::st_set_crs(b, sf::st_crs(pts_for_seeds))
          sp_geom <- sf::st_sample(b, size = lvl, type = "random", exact = TRUE)
          sp <- sf::st_sf(geometry = sp_geom) |> sf::st_set_crs(sf::st_crs(pts_for_seeds))
        } else {
          sp <- voronoi_seeds_random(boundary, k = lvl)
          hw <- harmonize_crs(pts_for_seeds, sp); pts_for_seeds <- hw$a; sp <- hw$b
        }
      }
      
      # Clip target: always coerce to sf (not sfc)
      clip_target <- if (!is.null(boundary)) boundary else {
        tmp <- sf::st_as_sfc(sf::st_bbox(pts_for_seeds))
        tmp <- sf::st_set_crs(tmp, sf::st_crs(pts_for_seeds))
        sf::st_sf(geometry = tmp)
      }
      log_info("build_tessellation(): clip_target class: {paste(class(clip_target), collapse=',')}")
      
      if (nrow(sp) < 2) {
        log_info("build_tessellation(): Voronoi with <2 seeds; using clip_target as single polygon")
        polys <- dplyr::mutate(clip_target, poly_id = 1L)
      } else {
        polys <- create_voronoi_polygons(sp, clip_with = clip_target)
      }
      
    } else {
      # hex or square grids
      if (is.null(boundary)) {
        tmp <- sf::st_as_sfc(sf::st_bbox(pts_for_seeds))
        tmp <- sf::st_set_crs(tmp, sf::st_crs(pts_for_seeds))
        boundary <- sf::st_sf(geometry = tmp)
      }
      log_info("build_tessellation(): grid path, boundary class pre-grid: {paste(class(boundary), collapse=',')}")
      polys <- create_grid_polygons(boundary, target_cells = lvl, type = if (method == "hex") "hex" else "square")
      hh <- harmonize_crs(polys, pts_for_seeds); polys <- hh$a; pts_for_seeds <- hh$b
    }
    
    # Safety: ensure polygons are sf before assignment
    if (inherits(polys, "sfc")) {
      log_warn("build_tessellation(): 'polys' is sfc; converting to sf before assignment")
      polys <- sf::st_sf(geometry = polys)
    }
    
    assigned <- assign_features_to_polygons(pts_for_seeds, polys)
    res[[as.character(lvl)]] <- list(polygons = polys, data = assigned, method = method, level = lvl)
  }
  res
}

# -----------------------------------------------------------------------------
# Modeling Pipeline Wrapper
# -----------------------------------------------------------------------------

#' Evaluate spatial models across tessellations and levels
#'
#' Builds one or more tessellations (Voronoi / hex / square) at specified
#' \code{levels}, assigns features to cells, and fits requested models
#' (GWR and/or Bayesian) to compare goodness-of-fit metrics across designs.
#' Handles CRS projection robustly and converts non-point geometries to points
#' for seeding/assignment.
#'
#' @param data_sf \code{sf} object containing the analysis data, including
#'   \code{response_var}, \code{predictor_vars}, and a geometry column. May be
#'   POINT/LINE/POLYGON; non-point geometries are pointized internally.
#' @param response_var Character scalar; name of the response variable column in
#'   \code{data_sf}.
#' @param predictor_vars Character vector; names of predictor columns to include
#'   in the model formula.
#' @param levels Integer or numeric vector of target cell counts / seed counts
#'   to evaluate. If \code{NULL}, an automatic selection is derived via
#'   \code{\link{determine_optimal_levels}} using \code{max_levels} and
#'   \code{top_n}.
#' @param tessellation Character vector of tessellation methods to evaluate;
#'   any of \code{"voronoi"}, \code{"hex"}, \code{"square"}. Multiple values are
#'   allowed; each is evaluated at all \code{levels}.
#' @param boundary Optional \code{sf}/\code{sfc} POLYGON/MULTIPOLYGON to clip
#'   tessellations and to define the working CRS. If \code{NULL}, a bbox derived
#'   from \code{data_sf} is used as the clip area.
#' @param seeds Seeding strategy for Voronoi tessellations; one of
#'   \code{"kmeans"}, \code{"random"}, or \code{"provided"}. Ignored for grids.
#' @param provided_seed_points Optional \code{sf} POINT/MULTIPOINT to use when
#'   \code{seeds = "provided"}. If count differs from the level, evaluation
#'   proceeds with the provided count (a warning is logged).
#' @param pointize Strategy for converting non-point geometries prior to
#'   assignment/seeding; one of \code{"auto"}, \code{"surface"},
#'   \code{"centroid"}, or \code{"line_midpoint"}. See
#'   \code{\link{coerce_to_points}}.
#' @param seed_kmeans_from Deprecated compatibility parameter; currently treated
#'   as \code{"features"} and passed through to \code{\link{build_tessellation}}.
#' @param models Character vector indicating which models to fit; any of
#'   \code{"GWR"}, \code{"Bayesian"}.
#' @param n.samples Integer; number of MCMC samples for the Bayesian model
#'   (passed to \code{\link{fit_bayesian_spatial_model}}).
#' @param cov_model Covariance model for the Bayesian fit; one of
#'   \code{"exponential"}, \code{"spherical"}, \code{"matern"}.
#' @param max_levels Integer; maximum \emph{k} to consider when selecting levels
#'   automatically.
#' @param top_n Integer; number of candidate levels to return from the elbow
#'   heuristic when \code{levels = NULL}.
#'
#' @details
#' \strong{CRS handling:} If \code{boundary} is provided, both \code{boundary}
#' and \code{data_sf} are projected to a consistent projected CRS via
#' \code{\link{ensure_projected}}. Otherwise, a locally appropriate projected
#' CRS is chosen heuristically.
#'
#' \strong{Tessellations:} For each method in \code{tessellation} and each
#' \code{level}, a tessellation is produced with
#' \code{\link{build_tessellation}}; features are assigned to cells using
#' \code{\link{assign_features_to_polygons}} (unassigned features are dropped).
#'
#' \strong{Models and metrics:}
#' \itemize{
#'   \item \emph{GWR}: fitted via \pkg{spgwr}; the reported \code{Metric} is
#'   AICc from the fit.
#'   \item \emph{Bayesian}: fitted via \pkg{spBayes}; the reported
#'   \code{Metric} is a DIC-like proxy
#'   \eqn{-2\,\mathrm{mean}(\ell) + \mathrm{var}(\ell)} computed from the
#'   log-likelihood over posterior draws.
#' }
#' Lower \code{Metric} indicates better fit for both models.
#'
#' @return A \code{list} with components:
#' \describe{
#'   \item{\code{results}}{A tibble with columns
#'     \code{Tessellation} (method name),
#'     \code{Level} (integer),
#'     \code{Model} (\code{"GWR"} or \code{"Bayesian"}),
#'     \code{Metric} (numeric score; lower is better).}
#'   \item{\code{tess_data}}{A nested list of tessellation outputs keyed by
#'     method then level; each entry contains \code{polygons} (sf with
#'     \code{poly_id}) and \code{data} (assigned points with \code{polygon_id}).}
#'   \item{\code{levels}}{The \code{levels} actually evaluated (useful if
#'     auto-selected).}
#'   \item{\code{params}}{A list echoing key parameters used for the run.}
#' }
#'
#' @section Logging:
#' Informative messages (via \pkg{logger}) are emitted about CRS choices,
#' missing assignments, and per-level/model progress.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' set.seed(1)
#' # Build synthetic data
#' bb <- st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 10000, ymax = 8000)), crs = 3857)
#' pts <- st_sample(bb, 400) |> st_sf() |> st_set_crs(3857)
#' pts$yvar <- rnorm(nrow(pts))
#' pts$x1 <- rnorm(nrow(pts))
#' pts$x2 <- rnorm(nrow(pts))
#'
#' # Compare Voronoi vs Hex at two levels for GWR only
#' out <- evaluate_models(
#'   data_sf = pts |> dplyr::rename(resp = yvar),
#'   response_var = "resp",
#'   predictor_vars = c("x1","x2"),
#'   levels = c(6, 12),
#'   tessellation = c("voronoi","hex"),
#'   boundary = st_sf(geometry = bb),
#'   models = "GWR",
#'   n.samples = 1000,
#'   cov_model = "exponential"
#' )
#' dplyr::arrange(out$results, Metric)
#' }
#'
#' @seealso
#' \code{\link{build_tessellation}},
#' \code{\link{coerce_to_points}},
#' \code{\link{determine_optimal_levels}},
#' \code{\link{fit_gwr_model}},
#' \code{\link{fit_bayesian_spatial_model}},
#' \code{\link{plot_tessellation_map}},
#' \code{\link{summarize_by_cell}}
#'
#' @export
evaluate_models <- function(
    data_sf,
    response_var,
    predictor_vars,
    levels = NULL,
    tessellation = c("voronoi","hex","square"),
    boundary = NULL,
    seeds = c("kmeans","random","provided"),
    provided_seed_points = NULL,
    pointize = c("auto","surface","centroid","line_midpoint"),
    seed_kmeans_from = c("features","assigned_points"),
    models = c("GWR","Bayesian"),
    n.samples = 2500,
    cov_model = c("exponential","spherical","matern"),
    max_levels = 10,
    top_n = 3
) {
  tessellation <- unique(match.arg(tessellation, several.ok = TRUE))
  seeds <- match.arg(seeds)
  pointize <- match.arg(pointize)
  seed_kmeans_from <- match.arg(seed_kmeans_from)
  cov_model <- match.arg(cov_model)
  
  # Prepare data and boundary CRS
  if (!is.null(boundary)) {
    boundary <- ensure_projected(boundary)
    data_sf  <- ensure_projected(data_sf, boundary)
  } else {
    data_sf  <- ensure_projected(data_sf)
  }
  
  # Coerce to representative POINTS for assignment & seeding
  data_pts <- coerce_to_points(data_sf, strategy = pointize)
  
  # Select default levels if not provided (via elbow on WSS)
  if (is.null(levels)) {
    levels <- determine_optimal_levels(data_pts, max_levels = max_levels, top_n = top_n)
  }
  
  # Build tessellations
  tess_data <- list()
  for (tess in tessellation) {
    tess_data[[tess]] <- build_tessellation(
      features_sf = data_pts, levels = levels, method = tess,
      boundary = boundary, seeds = seeds,
      provided_seed_points = provided_seed_points,
      pointize = pointize, seed_kmeans_from = seed_kmeans_from
    )
  }
  
  # Fit models for each tessellation & level
  results <- tibble::tibble(
    Tessellation = character(),
    Level        = integer(),
    Model        = character(),
    Metric       = numeric()
  )
  
  for (tess in tessellation) {
    for (lvl in levels) {
      obj <- tess_data[[tess]][[as.character(lvl)]]
      if (is.null(obj) || nrow(obj$data) == 0) {
        log_warn(sprintf("No assigned data for %s level %s", tess, lvl))
        next
      }
      df <- obj$data
      
      if ("GWR" %in% models) {
        gwr_res <- fit_gwr_model(df, response_var, predictor_vars)
        results <- dplyr::bind_rows(results, tibble::tibble(
          Tessellation = tess, Level = lvl, Model = "GWR", Metric = gwr_res$AICc
        ))
      }
      if ("Bayesian" %in% models) {
        bayes_res <- fit_bayesian_spatial_model(df, response_var, predictor_vars, n.samples = n.samples, cov_model = cov_model)
        results <- dplyr::bind_rows(results, tibble::tibble(
          Tessellation = tess, Level = lvl, Model = "Bayesian", Metric = bayes_res$DIC
        ))
      }
    }
  }
  
  list(
    results   = results,
    tess_data = tess_data,
    levels    = levels,
    params    = list(
      response_var = response_var,
      predictor_vars = predictor_vars,
      tessellation = tessellation,
      seeds = seeds,
      pointize = pointize,
      cov_model = cov_model,
      n.samples = n.samples
    )
  )
}

# -----------------------------------------------------------------------------
# Optional: Aggregations & Plotting
# -----------------------------------------------------------------------------

#' Summarize assigned points by tessellation cell
#'
#' Aggregates feature-level data by \code{polygon_id} and returns counts and
#' (optionally) mean summaries of the response and selected predictors.
#' Designed to consume the output of \code{\link{assign_features_to_polygons}}.
#'
#' @param assigned_points_sf An \code{sf} object of assigned features that
#'   includes a \code{polygon_id} column (as produced by
#'   \code{\link{assign_features_to_polygons}}). Additional numeric columns may
#'   be present for response/predictors.
#' @param response_var Optional character scalar naming a numeric response
#'   column to summarize with a mean. If \code{NULL} or not present, the
#'   response mean is omitted.
#' @param predictor_vars Optional character vector of predictor column names to
#'   summarize with means. Names not found in \code{assigned_points_sf} are
#'   silently ignored.
#'
#' @details
#' The function drops geometry and computes:
#' \itemize{
#'   \item \code{n}: count of assigned features per \code{polygon_id}.
#'   \item \code{mean_response}: mean of \code{response_var} (if provided and present).
#'   \item \code{mean_<pred>}: means for each predictor in \code{predictor_vars}
#'     that exists in the data.
#' }
#' All means are computed with \code{na.rm = TRUE}. The result is a plain
#' data frame/tibble (not \code{sf}).
#'
#' @return A tibble/data frame with one row per \code{polygon_id}, always
#'   including \code{polygon_id} and \code{n}, and conditionally including
#'   \code{mean_response} and \code{mean_<predictor>} columns.
#'
#' @examples
#' \dontrun{
#' # Suppose `asgn` is the `data` element returned by build_tessellation()[[k]]
#' # and contains columns: polygon_id, resp, x1, x2, geometry
#' tbl <- summarize_by_cell(asgn, response_var = "resp", predictor_vars = c("x1","x2"))
#' head(tbl)
#' }
#'
#' @seealso
#' \code{\link{assign_features_to_polygons}},
#' \code{\link{build_tessellation}},
#' \code{\link{plot_tessellation_map}}
#'
#' @export
summarize_by_cell <- function(assigned_points_sf, response_var = NULL, predictor_vars = NULL) {
  df <- sf::st_drop_geometry(assigned_points_sf)
  out <- df |>
    dplyr::group_by(polygon_id) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  
  if (!is.null(response_var) && response_var %in% names(df)) {
    resp <- df |>
      dplyr::group_by(polygon_id) |>
      dplyr::summarise(mean_response = mean(.data[[response_var]], na.rm = TRUE), .groups = "drop")
    out <- dplyr::left_join(out, resp, by = "polygon_id")
  }
  
  if (!is.null(predictor_vars)) {
    keep <- predictor_vars[predictor_vars %in% names(df)]
    if (length(keep)) {
      preds <- df |>
        dplyr::group_by(polygon_id) |>
        dplyr::summarise(dplyr::across(dplyr::all_of(keep), ~mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
                         .groups = "drop")
      out <- dplyr::left_join(out, preds, by = "polygon_id")
    }
  }
  out
}

#' Plot tessellation polygons and assigned points
#'
#' Draws tessellation \code{polygons_sf} and the assigned \code{points_sf} on a
#' clean map, optionally clipped/outlined by a \code{boundary}, over a vector
#' \code{basemap}, or on top of OSM tiles. Labels can show each polygon's ID or
#' the count of assigned points per polygon.
#'
#' @param polygons_sf An \code{sf} or \code{sfc} object of tessellation polygons.
#'   If a \code{poly_id} column is absent it will be created as
#'   \code{seq_len(nrow(polygons_sf))}.
#' @param points_sf An \code{sf} or \code{sfc} object of features assigned to
#'   polygons; must include a \code{polygon_id} column (as produced by
#'   \code{\link{assign_features_to_polygons}}).
#' @param title Optional character; plot title.
#' @param boundary Optional \code{sf}/\code{sfc} geometry used as an outline
#'   (drawn with no fill) to provide geographic context or clipping reference.
#' @param basemap Optional \code{sf}/\code{sfc} layer (e.g., counties) drawn
#'   beneath polygons with a light fill to give visual context.
#' @param use_osm_tiles Logical; if \code{TRUE} and \pkg{ggspatial} is available,
#'   renders web tiles beneath the data and reprojects layers to EPSG:3857.
#'   If \pkg{ggspatial} is not installed, falls back to the plain \code{sf}
#'   rendering path.
#' @param osm_type Character tile provider passed to
#'   \code{ggspatial::annotation_map_tile()} (e.g., \code{"osm"}).
#' @param osm_zoom Optional integer for tile zoom delta (passed to
#'   \code{zoomin} in \pkg{ggspatial}); \code{NULL} lets \pkg{ggspatial} choose.
#' @param basemap_alpha Numeric in \code{[0,1]} controlling the transparency of
#'   \code{basemap}.
#' @param label One of \code{"poly_id"}, \code{"count"}, or \code{"none"}.
#'   Controls the text drawn at each polygon's interior point: the polygon ID,
#'   the count of assigned points (computed on the fly), or no labels.
#' @param show_counts Deprecated logical/NULL for backward compatibility.
#'   \code{TRUE} is equivalent to \code{label = "count"}; \code{FALSE} to
#'   \code{label = "none"}. If not \code{NULL}, overrides \code{label}.
#'
#' @details
#' Inputs supplied as \code{sfc} are coerced to \code{sf}. Coordinate reference
#' systems are harmonized across \code{polygons_sf}, \code{points_sf},
#' \code{boundary}, and \code{basemap} via \code{\link{harmonize_crs}}.
#' Label anchor points are computed with \code{sf::st_point_on_surface()}.
#' To avoid warnings and deformed label placement on geographic (lon/lat) CRS,
#' the function temporarily projects polygons to a suitable local projected CRS,
#' computes interior points, then transforms labels back.
#'
#' When \code{use_osm_tiles = TRUE} and \pkg{ggspatial} is available, all layers
#' are transformed to Web Mercator (EPSG:3857) and tiles are drawn with
#' \code{ggspatial::annotation_map_tile(type = osm_type, zoomin = osm_zoom)}.
#' Otherwise, a plain \code{ggplot2} rendering path is used with optional
#' \code{basemap} and \code{boundary}.
#'
#' Fill colors for polygons use \code{ggplot2::scale_fill_viridis_d()} (guide
#' suppressed), point colors use \code{ggplot2::scale_color_viridis_d(name =
#' "Polygon ID")}, and a minimalist theme is applied. The legend shows the
#' mapping of \code{polygon_id} for points.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' # Suppose bt <- build_tessellation(...); choose a level:
#' obj <- bt[["12"]]
#' p <- plot_tessellation_map(
#'   polygons_sf = obj$polygons,
#'   points_sf   = obj$data,
#'   title       = "Voronoi (k = 12)",
#'   boundary    = some_boundary_sf,
#'   basemap     = some_admin_units_sf,
#'   label       = "poly_id"
#' )
#' print(p)
#'
#' # With OSM tiles (requires ggspatial):
#' p2 <- plot_tessellation_map(obj$polygons, obj$data,
#'   title = "Hex grid with tiles", boundary = some_boundary_sf,
#'   use_osm_tiles = TRUE, osm_type = "osm", label = "count")
#' }
#'
#' @seealso
#' \code{\link{build_tessellation}},
#' \code{\link{assign_features_to_polygons}},
#' \code{\link{summarize_by_cell}}
#'
#' @importFrom ggplot2 ggplot geom_sf geom_sf_text labs coord_sf theme_minimal
#'   theme scale_fill_viridis_d scale_color_viridis_d guides guide_legend
#'   element_blank element_line element_text
#' @export
plot_tessellation_map <- function(
    polygons_sf,
    points_sf,
    title = NULL,
    boundary = NULL,
    basemap = NULL,
    use_osm_tiles = FALSE,
    osm_type = "osm",
    osm_zoom = NULL,
    basemap_alpha = 0.35,
    label = c("poly_id","count","none"),
    show_counts = NULL   # back-compat; TRUE -> "count", FALSE -> "none"
) {
  label <- match.arg(label)
  
  # back-compat switch
  if (!is.null(show_counts)) {
    label <- if (isTRUE(show_counts)) "count" else "none"
  }
  
  # Coerce sfc -> sf
  if (inherits(polygons_sf, "sfc")) polygons_sf <- sf::st_sf(geometry = polygons_sf)
  if (inherits(points_sf,   "sfc")) points_sf   <- sf::st_sf(geometry = points_sf)
  if (!is.null(boundary) && inherits(boundary, "sfc")) boundary <- sf::st_sf(geometry = boundary)
  if (!is.null(basemap)  && inherits(basemap,  "sfc")) basemap  <- sf::st_sf(geometry = basemap)
  
  # Ensure poly_id exists
  if (!("poly_id" %in% names(polygons_sf))) {
    polygons_sf$poly_id <- seq_len(nrow(polygons_sf))
  }
  
  # Harmonize CRS
  if (!is.null(boundary)) {
    hb <- harmonize_crs(polygons_sf, boundary); polygons_sf <- hb$a; boundary <- hb$b
  }
  hb2 <- harmonize_crs(polygons_sf, points_sf); polygons_sf <- hb2$a; points_sf <- hb2$b
  if (!is.null(basemap)) {
    hb3 <- harmonize_crs(polygons_sf, basemap); polygons_sf <- hb3$a; basemap <- hb3$b
  }
  
  # Add counts only when needed
  if (identical(label, "count")) {
    counts <- points_sf |>
      sf::st_drop_geometry() |>
      dplyr::count(polygon_id, name = "n_pts")
    polygons_sf <- polygons_sf |>
      dplyr::left_join(counts, by = c("poly_id" = "polygon_id")) |>
      dplyr::mutate(n_pts = dplyr::coalesce(n_pts, 0L))
  }
  
  # Label geometry helper (avoid long/lat warnings)
  make_labels <- function(poly_sf) {
    g <- sf::st_geometry(poly_sf)
    if (.is_longlat(poly_sf)) {
      g_proj <- sf::st_geometry(ensure_projected(poly_sf))
      lbl    <- sf::st_point_on_surface(g_proj)
      sf::st_transform(lbl, sf::st_crs(poly_sf))
    } else {
      sf::st_point_on_surface(g)
    }
  }
  
  # OSM tiles path (EPSG:3857) or plain sf path
  if (isTRUE(use_osm_tiles) && requireNamespace("ggspatial", quietly = TRUE)) {
    target_crs <- sf::st_crs(3857)
    polygons_d <- sf::st_transform(polygons_sf, target_crs)
    points_d   <- sf::st_transform(points_sf,   target_crs)
    boundary_d <- if (!is.null(boundary)) sf::st_transform(boundary, target_crs) else NULL
    basemap_d  <- if (!is.null(basemap))  sf::st_transform(basemap,  target_crs) else NULL
    
    lab_geom <- make_labels(polygons_d)
    lab_vals <- if (identical(label, "count")) polygons_d$n_pts else polygons_d$poly_id
    label_pts <- sf::st_sf(.lab = lab_vals, geometry = lab_geom)
    
    ggplot2::ggplot() +
      ggspatial::annotation_map_tile(type = osm_type, zoomin = osm_zoom) +
      { if (!is.null(basemap_d))  ggplot2::geom_sf(data = basemap_d,  fill = "grey85", color = "grey40", alpha = basemap_alpha, linewidth = 0.2) else NULL } +
      { if (!is.null(boundary_d)) ggplot2::geom_sf(data = boundary_d, fill = NA, color = "grey20",  linewidth = 0.5) else NULL } +
      ggplot2::geom_sf(data = polygons_d, ggplot2::aes(fill = factor(poly_id)), color = "white", alpha = 0.28, linewidth = 0.4) +
      ggplot2::geom_sf(data = points_d,   ggplot2::aes(color = factor(polygon_id)), size = 1.6, alpha = 0.9) +
      { if (!identical(label, "none")) ggplot2::geom_sf_text(data = label_pts, ggplot2::aes(label = .lab), size = 3.4, fontface = "bold") else NULL } +
      ggplot2::scale_fill_viridis_d(guide = "none") +
      ggplot2::scale_color_viridis_d(name = "Polygon ID") +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3))) +
      ggplot2::labs(title = title %||% "Tessellation Map") +
      ggplot2::coord_sf(crs = target_crs, expand = FALSE) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position  = "right",
        plot.title = ggplot2::element_text(face = "bold")
      )
    
  } else {
    lab_geom <- make_labels(polygons_sf)
    lab_vals <- if (identical(label, "count")) polygons_sf$n_pts else polygons_sf$poly_id
    label_pts <- sf::st_sf(.lab = lab_vals, geometry = lab_geom)
    
    ggplot2::ggplot() +
      { if (!is.null(basemap))  ggplot2::geom_sf(data = basemap,  fill = "grey85", color = "grey40", alpha = basemap_alpha, linewidth = 0.2) else NULL } +
      { if (!is.null(boundary)) ggplot2::geom_sf(data = boundary, fill = NA, color = "grey20",  linewidth = 0.5) else NULL } +
      ggplot2::geom_sf(data = polygons_sf, ggplot2::aes(fill = factor(poly_id)), color = "white", alpha = 0.28, linewidth = 0.4) +
      ggplot2::geom_sf(data = points_sf,   ggplot2::aes(color = factor(polygon_id)), size = 1.6, alpha = 0.9) +
      { if (!identical(label, "none")) ggplot2::geom_sf_text(data = label_pts, ggplot2::aes(label = .lab), size = 3.4, fontface = "bold") else NULL } +
      ggplot2::scale_fill_viridis_d(guide = "none") +
      ggplot2::scale_color_viridis_d(name = "Polygon ID") +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3))) +
      ggplot2::labs(title = title %||% "Tessellation Map") +
      ggplot2::coord_sf(expand = FALSE) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.2),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position  = "right",
        plot.title = ggplot2::element_text(face = "bold")
      )
  }
}
