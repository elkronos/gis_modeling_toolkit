# =============================================================================
# Geospatial Modeling Toolkit (Robust & Modular)
#
# Core capabilities
# - Tessellation builders:
#     • Voronoi (bounded / expanded / clipped)
#     • Delaunay triangles (geometry::delaunayn when available; sf fallback)
#     • Hexagonal grids (inside boundary)
#     • Square grids (inside boundary)
# - Assignment:
#     • Robustly map features to cells with careful CRS alignment & validity fixes
# - Seeding:
#     • K-means, random, or user-provided point seeds
# - Models:
#     • Geographically Weighted Regression (spgwr)
#     • Bayesian spatial regression (spBayes) with exponential / spherical / Matérn
# - Level selection:
#     • Elbow heuristic on within-cluster sum of squares (WSS)
# - Utilities:
#     • Pointization for non-point geometries
#     • Safe polygon validity repair (sf/lwgeom or buffered fallback)
#     • Automatic projection to a suitable projected CRS when inputs are lon/lat
#     • Summaries & plotting helpers (ggplot2)
#
# Design notes
# - Works with arbitrary datasets; not tied to any region.
# - Fails safe: geometry validity, duplicate handling, degenerate bbox handling.
# - CRS-aware at every step; avoids accidental lon/lat geometry ops.
#
# Required packages
#   logger, sf, sp, spgwr, spBayes, deldir, ggplot2, dplyr, tidyr, mvtnorm
#
# Optional (recommended) packages used when available
#   lwgeom   : faster/safer st_make_valid where supported
#   geometry : fast Delaunay triangulation (delaunayn); falls back to sf otherwise
#
# =============================================================================

suppressPackageStartupMessages({
  library(logger)    # structured messages
  library(sf)        # modern vector GIS
  library(sp)        # legacy classes (spgwr/spBayes interop)
  library(spgwr)     # geographically weighted regression
  library(spBayes)   # Bayesian spatial models
  library(deldir)    # Voronoi/Dirichlet + Delaunay helpers used by Voronoi workflows
  library(ggplot2)   # plotting utilities
  library(dplyr)     # data wrangling
  library(tidyr)     # data reshaping
  library(mvtnorm)   # multivariate normals in Bayesian pieces
})

# Note: lwgeom and geometry are used opportunistically via requireNamespace()
#       inside helper functions; they are not hard dependencies.

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

#' Pick a sensible local projected CRS for an `sf`/`sfc` object
#'
#' Chooses an appropriate projected coordinate reference system for spatial
#' data, favoring a UTM zone based on the dataset's geographic centroid when
#' possible, and falling back to Web Mercator (EPSG:3857) when longitude/latitude
#' information cannot be reliably determined.
#'
#' @details
#' The selection proceeds as follows:
#' \enumerate{
#'   \item If `x` is not an `sf` or `sfc` object, return \code{sf::NA_crs_}.
#'   \item If `x` already has a non-long/lat CRS, return that CRS unchanged.
#'   \item Otherwise attempt to ensure the data are in WGS84 (\code{EPSG:4326})
#'         to read longitude/latitude:
#'         \itemize{
#'           \item If a CRS is set on `x`, transform to 4326.
#'           \item If no CRS is set or the transform fails, proceed with `x`
#'                 as-is (which will be treated as unknown CRS in the next step).
#'         }
#'   \item If `x` (or the transformed copy) is not in long/lat, return
#'         Web Mercator (\code{EPSG:3857}) as a safe fallback.
#'   \item Compute the centroid of the unioned geometry in long/lat and derive
#'         the UTM zone from the centroid longitude:
#'         \code{zone = floor((lon + 180) / 6) + 1}, clamped to 1–60.
#'         Use the southern-hemisphere series \code{EPSG:32700 + zone} if
#'         the centroid latitude is < 0, otherwise use the northern series
#'         \code{EPSG:32600 + zone}.
#' }
#'
#' This helper relies on an internal predicate \code{.is_longlat()} to detect
#' geographic CRSs. It never reprojects your data; it only returns an
#' \code{sf::crs} object you can use for projection elsewhere (e.g.,
#' \code{sf::st_transform(x, .pick_local_projected_crs(x))}).
#'
#' @param x An \code{sf} or \code{sfc} object.
#'
#' @return An \code{sf::crs} object. Specifically:
#' \itemize{
#'   \item The existing projected CRS of `x` (if `x` is already projected);
#'   \item A UTM CRS (\code{EPSG:32601–32660} or \code{EPSG:32701–32760})
#'         chosen from the centroid of `x` in lon/lat; or
#'   \item \code{sf::st_crs(3857)} if lon/lat cannot be determined.
#'   \item \code{sf::NA_crs_} if `x` is not an \code{sf}/\code{sfc} object.
#' }
#'
#' @seealso \code{\link[sf]{st_crs}}, \code{\link[sf]{st_transform}}
#'
#' @examples
#' \donttest{
#' library(sf)
#' pts <- st_as_sf(
#'   data.frame(lon = c(12.5, 12.6), lat = c(41.9, 42.0)),
#'   coords = c("lon", "lat"), crs = 4326
#' )
#' crs_choice <- .pick_local_projected_crs(pts)
#' crs_choice$epsg  # likely a UTM 33N code (e.g., 32633) for Rome area
#'
#' # Already projected: original CRS is returned unchanged
#' pts_proj <- st_transform(pts, 32633)
#' identical(st_crs(pts_proj), .pick_local_projected_crs(pts_proj))
#' }
#'
#' @keywords internal
#' @noRd
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

#' Ensure an object has a projected CRS (with sensible defaults)
#'
#' Coerces spatial objects to a projected coordinate reference system (CRS)
#' suitable for distance/area calculations. If `target_crs` is supplied, the
#' object is aligned to that CRS. Otherwise, the function attempts to detect
#' whether the data are in geographic lon/lat and, if so, transforms to a
#' locally appropriate projected CRS chosen by \code{.pick_local_projected_crs()}.
#'
#' @details
#' Behavior by case:
#' \itemize{
#'   \item \strong{Non-spatial input} (\code{!inherits(x, "sf")} and
#'         \code{!inherits(x, "sfc")}): returned unchanged.
#'
#'   \item \strong{`target_crs` provided}: Parsed with \code{sf::st_crs()}.
#'         \itemize{
#'           \item If \code{target_crs} cannot be parsed (\code{NA} CRS), return
#'                 \code{x} unchanged.
#'           \item If \code{x} has \code{NA} CRS, set its CRS to
#'                 \code{target_crs} \emph{without} reprojecting (assign only),
#'                 then return.
#'           \item If \code{x} has a CRS different from \code{target_crs},
#'                 transform via \code{sf::st_transform(x, target_crs)} and
#'                 return. If already identical, return unchanged.
#'         }
#'
#'   \item \strong{`target_crs` not provided}:
#'         \itemize{
#'           \item If \code{x} has \code{NA} CRS, a heuristic checks whether the
#'                 bounding box looks like degrees: \code{xmin/xmax} within
#'                 \code{[-180, 180]} and \code{ymin/ymax} within \code{[-90, 90]}.
#'                 If so, assign \code{EPSG:4326} (WGS84), then transform to a
#'                 locally appropriate projected CRS using
#'                 \code{.pick_local_projected_crs(x)}.
#'                 Otherwise (likely already projected but unknown), leave as-is.
#'           \item If \code{x} has a CRS and it is geographic
#'                 (\code{.is_longlat(x) == TRUE}), transform to the local
#'                 projected CRS returned by \code{.pick_local_projected_crs(x)}.
#'                 If \code{x} is already projected, return unchanged.
#'         }
#' }
#'
#' Notes:
#' \itemize{
#'   \item This helper never forces a transform if the input is not an
#'         \code{sf}/\code{sfc} object.
#'   \item When \code{x} has \code{NA} CRS and \code{target_crs} is given, the
#'         CRS is \emph{assigned} (no coordinate change), which is the safest
#'         behavior when the original CRS is unknown.
#'   \item The lon/lat heuristic relies on \code{sf::st_bbox(x)} being finite.
#'         If the bbox cannot be computed, or values are out of range, the
#'         object is returned unchanged.
#'   \item Transformations may throw if the CRS definitions are invalid or the
#'         pipeline cannot be built by PROJ; such errors are not caught here.
#' }
#'
#' @param x An \code{sf} or \code{sfc} object (other objects are returned
#'   unchanged).
#' @param target_crs Optional target CRS. Anything accepted by
#'   \code{sf::st_crs()}, e.g., an EPSG integer (e.g., \code{32633}),
#'   an \code{sf::crs} object, a WKT string, or another object with a valid CRS.
#'
#' @return \code{x}, potentially with a new projected CRS. The return type/class
#'   matches the input. In summary:
#'   \itemize{
#'     \item If \code{target_crs} is provided, an object in that CRS (or CRS
#'           assigned if \code{x} had \code{NA} CRS).
#'     \item If \code{target_crs} is not provided and \code{x} is lon/lat, an
#'           object transformed to a local projected CRS (typically UTM).
#'     \item Otherwise, the original object unchanged.
#'   }
#'
#' @seealso \code{\link[sf]{st_crs}}, \code{\link[sf]{st_transform}},
#'   \code{\link[sf]{st_bbox}}, and the internal helpers
#'   \code{.pick_local_projected_crs()} and \code{.is_longlat()}.
#'
#' @examples
#' \donttest{
#' library(sf)
#' # Example lon/lat points near Rome
#' pts_ll <- st_as_sf(data.frame(lon = c(12.5, 12.6), lat = c(41.9, 42.0)),
#'                    coords = c("lon", "lat"), crs = 4326)
#'
#' # Auto-pick a local projected CRS (likely UTM 33N) and transform
#' pts_proj <- ensure_projected(pts_ll)
#' st_is_longlat(pts_proj) # FALSE
#'
#' # Force a specific target CRS
#' pts_32633 <- ensure_projected(pts_ll, target_crs = 32633)
#'
#' # If CRS is missing but bbox looks like degrees, assume WGS84 then project
#' st_crs(pts_ll) <- NA_crs_
#' pts_guess <- ensure_projected(pts_ll)
#'
#' # Non-sf objects are returned unchanged
#' ensure_projected(data.frame(x = 1:3))
#' }
#'
#' @keywords internal
#' @noRd
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

#' Harmonize CRS between two spatial objects
#'
#' Aligns two \pkg{sf} objects to a common coordinate reference system (CRS).
#' When \code{target_crs} is provided, both inputs are aligned to that CRS.
#' Otherwise, if their CRSs differ, one side is transformed (or, if needed,
#' assigned) to match the other according to \code{prefer}.
#'
#' @details
#' Logic by case:
#' \itemize{
#'   \item \strong{Explicit \code{target_crs}}: Both \code{a} and \code{b} are
#'         aligned to \code{target_crs}. If an input has an \code{NA} CRS, it is
#'         \emph{assigned} that CRS via \code{sf::st_set_crs()} (no reprojection).
#'         If it has a valid CRS, the function attempts \code{sf::st_transform()};
#'         on failure, it falls back to assigning the CRS and emits a warning.
#'
#'   \item \strong{Both CRSs \code{NA}}: Returned unchanged.
#'
#'   \item \strong{Exactly one CRS \code{NA}}: The \code{NA} side is assigned
#'         the other's CRS (no reprojection).
#'
#'   \item \strong{Both CRSs set and identical}: Returned unchanged.
#'
#'   \item \strong{Both CRSs set and different}: Transform the non-preferred side
#'         to the preferred side's CRS using \code{sf::st_transform()}; on error,
#'         fall back to assigning the CRS and warn.
#' }
#'
#' \strong{Important}: Assigning a CRS with \code{sf::st_set_crs()} does \emph{not}
#' change coordinates. Only use assignment when you are certain the existing
#' coordinates are already expressed in that CRS.
#'
#' @param a,b Objects of class \code{sf} or \code{sfc}.
#' @param prefer Which object's CRS to keep when both CRSs are set but different
#'   and \code{target_crs} is \code{NULL}. One of \code{"a"} or \code{"b"}.
#' @param target_crs Optional target CRS to apply to both objects. Anything
#'   accepted by \code{sf::st_crs()} (e.g., EPSG integer like \code{3857},
#'   WKT string, \code{sf::crs} object).
#'
#' @return A named list with components:
#' \itemize{
#'   \item \code{a}: The input \code{a} aligned to the chosen CRS.
#'   \item \code{b}: The input \code{b} aligned to the chosen CRS.
#' }
#'
#' @examples
#' \donttest{
#' library(sf)
#'
#' # Two small point layers in different CRSs
#' a <- st_as_sf(data.frame(x = 0:1, y = 0:1), coords = c("x", "y"), crs = 4326)
#' b <- st_transform(a, 3857)
#'
#' # Prefer 'a' -> transform b to EPSG:4326
#' out1 <- harmonize_crs(a, b, prefer = "a")
#' st_crs(out1$a) == st_crs(out1$b)  # TRUE
#'
#' # Prefer 'b' -> transform a to EPSG:3857
#' out2 <- harmonize_crs(a, b, prefer = "b")
#'
#' # Force a specific target CRS for both
#' out3 <- harmonize_crs(a, b, target_crs = 3395)
#'
#' # One side has NA CRS -> it is assigned the other's CRS (no reprojection)
#' st_crs(b) <- NA_crs_
#' out4 <- harmonize_crs(a, b)  # b gets assigned EPSG:4326
#' }
#'
#' @seealso \code{\link[sf]{st_crs}}, \code{\link[sf]{st_transform}},
#'   \code{\link[sf]{st_set_crs}}, and \code{\link{ensure_projected}}.
#'
#' @export
harmonize_crs <- function(a, b, prefer = c("a", "b"), target_crs = NULL) {
  if (!inherits(a, c("sf", "sfc"))) stop("harmonize_crs(): `a` must be sf or sfc.")
  if (!inherits(b, c("sf", "sfc"))) stop("harmonize_crs(): `b` must be sf or sfc.")
  prefer <- match.arg(prefer)
  
  crs_a <- sf::st_crs(a)
  crs_b <- sf::st_crs(b)
  
  .safe_transform <- function(x, to) {
    tryCatch(
      sf::st_transform(x, to),
      error = function(e) {
        warning(sprintf("harmonize_crs(): st_transform() failed (%s); setting CRS instead.", conditionMessage(e)),
                call. = FALSE)
        sf::st_set_crs(x, to)
      }
    )
  }
  
  # If an explicit target CRS is provided, align both to it.
  if (!is.null(target_crs)) {
    a <- if (!is.na(crs_a)) .safe_transform(a, target_crs) else sf::st_set_crs(a, target_crs)
    b <- if (!is.na(crs_b)) .safe_transform(b, target_crs) else sf::st_set_crs(b, target_crs)
    return(list(a = a, b = b))
  }
  
  # Both NA -> unchanged
  if (is.na(crs_a) && is.na(crs_b)) {
    return(list(a = a, b = b))
  }
  
  # One NA -> set that CRS to the other's (no transform)
  if (is.na(crs_a) && !is.na(crs_b)) {
    a <- sf::st_set_crs(a, crs_b)
    return(list(a = a, b = b))
  }
  if (!is.na(crs_a) && is.na(crs_b)) {
    b <- sf::st_set_crs(b, crs_a)
    return(list(a = a, b = b))
  }
  
  # Both set: equal -> unchanged; different -> transform one side.
  if (identical(sf::st_crs(a), sf::st_crs(b))) {
    return(list(a = a, b = b))
  }
  
  if (prefer == "a") {
    b <- .safe_transform(b, crs_a)
  } else {
    a <- .safe_transform(a, crs_b)
  }
  
  list(a = a, b = b)
}

#' Fast center points of per-feature bounding boxes
#'
#' Computes the midpoint of the axis-aligned bounding box for each feature in
#' an \pkg{sf} object and returns those midpoints as an \code{sfc} point vector.
#' This is faster and more general than \code{sf::st_centroid()} for cases where
#' a simple center-of-extent is sufficient (e.g., labeling), and it works for
#' any geometry type (POINT/LINE/POLYGON, including multi- geometries).
#'
#' @details
#' For each geometry \eqn{g}, the function evaluates
#' \eqn{((xmin+xmax)/2,\ (ymin+ymax)/2)} from \code{sf::st_bbox(g)} and wraps the
#' result in \code{sf::st_point()}, assembling the results into an
#' \code{sfc} with the same CRS as \code{x}.
#'
#' \strong{Notes and caveats}
#' \itemize{
#'   \item The returned points are centers of the \emph{bounding boxes}, not
#'         true geometric centroids. They may lie outside highly concave or
#'         donut-shaped polygons.
#'   \item Coordinates are computed in the input CRS. If \code{x} is geographic
#'         (lon/lat), values are simple arithmetic midpoints and do not account
#'         for wrap-around at the antimeridian or high-latitude distortions.
#'         Consider projecting first via \code{sf::st_transform()}.
#'   \item Empty geometries may cause errors in \code{sf::st_bbox()}; ensure
#'         inputs are non-empty and valid as needed.
#' }
#'
#' @param x An object of class \code{sf}. Only the geometry column is used.
#'
#' @return An \code{sfc} (POINT) vector with one point per feature in \code{x},
#'   carrying the same CRS as \code{x}.
#'
#' @examples
#' \donttest{
#' library(sf)
#' # Simple polygon and line
#' p <- st_as_sf(data.frame(
#'   id = 1,
#'   wkt = "POLYGON((0 0, 2 0, 2 1, 0 1, 0 0))"
#' ), wkt = "wkt", crs = 4326)
#'
#' l <- st_as_sf(data.frame(
#'   id = 1,
#'   wkt = "LINESTRING(10 0, 12 2, 11 4)"
#' ), wkt = "wkt", crs = 4326)
#'
#' # Bounding-box centers
#' c_poly <- .bbox_center_sfc(p)
#' c_line <- .bbox_center_sfc(l)
#'
#' st_coordinates(c_poly)
#' st_coordinates(c_line)
#' }
#'
#' @seealso \code{\link[sf]{st_bbox}}, \code{\link[sf]{st_centroid}},
#'   \code{\link[sf]{st_point}}, \code{\link[sf]{st_transform}}
#'
#' @keywords internal
.bbox_center_sfc <- function(x) {
  stopifnot(inherits(x, "sf"))
  pts <- lapply(sf::st_geometry(x), function(g) {
    bb <- sf::st_bbox(g)
    sf::st_point(c((bb["xmin"] + bb["xmax"])/2, (bb["ymin"] + bb["ymax"])/2))
  })
  sf::st_sfc(pts, crs = sf::st_crs(x))
}

#' Coerce arbitrary geometries to representative points
#'
#' Converts the geometry column of an \pkg{sf} object to POINTs using one of
#' several strategies. This is useful for labeling, tessellation seeding, or
#' reducing non-point features to a single representative location while
#' preserving all attributes.
#'
#' @details
#' The \code{mode} controls how points are derived:
#'
#' \describe{
#'   \item{\code{"auto"} (default)}{
#'     Chooses a sensible strategy by geometry type:
#'     \itemize{
#'       \item \strong{POINT}: passed through unchanged.
#'       \item \strong{MULTIPOINT}: first coordinate is used.
#'       \item \strong{LINESTRING}: midpoint along the line (see
#'             \emph{Projection behavior} below).
#'       \item \strong{MULTILINESTRING}: midpoint of the longest constituent
#'             linestring (with projection behavior as for LINESTRING); if no
#'             parts are present, falls back to centroid.
#'       \item \strong{POLYGON}/\strong{MULTIPOLYGON}: robust
#'             \code{sf::st_point_on_surface()}.
#'       \item other types: \code{sf::st_centroid()} as a fallback.
#'     }
#'   }
#'   \item{\code{"centroid"}}{Geometric centroid via \code{sf::st_centroid()}.}
#'   \item{\code{"point_on_surface"} or \code{"surface"}}{
#'     Point guaranteed to lie on the surface via \code{sf::st_point_on_surface()}.
#'     The value \code{"surface"} is an alias.
#'   }
#'   \item{\code{"line_midpoint"}}{
#'     For LINESTRING geometries, uses \code{sf::st_line_sample(..., sample = 0.5)}
#'     to place a point halfway along the line length. For non-line features,
#'     falls back to centroid. If any feature is a MULTILINESTRING or a
#'     GEOMETRYCOLLECTION, the function errors (cast to LINESTRING first).
#'   }
#'   \item{\code{"bbox_center"}}{
#'     Numeric center of each feature's axis-aligned bounding box (see
#'     \code{.bbox_center_sfc()}); faster than centroids but may lie outside
#'     concave or donut polygons.
#'   }
#' }
#'
#' \strong{Projection behavior.}
#' Some operations (e.g., \code{"line_midpoint"} and the LINESTRING /
#' MULTILINESTRING branches of \code{"auto"}) are most meaningful in a projected
#' CRS. When \code{tmp_project = TRUE} and the input is geographic (lon/lat),
#' the function temporarily projects to a reasonable local projected CRS via
#' \code{ensure_projected()} to compute the midpoint, then transforms back to
#' the original CRS. When \code{tmp_project = FALSE} in lon/lat, a spherical-safe
#' fallback using \code{sf::st_centroid()} is used for lines.
#'
#' \strong{Notes and caveats.}
#' \itemize{
#'   \item \code{"bbox_center"} returns centers of extents, not true centroids.
#'   \item MULTIPOINT reduction selects the first point deterministically.
#'   \item Empty or invalid geometries may cause downstream errors; validate
#'         inputs as needed.
#' }
#'
#' @param x An \code{sf} object. All attributes are preserved; only the geometry
#'   column is replaced with POINTs.
#' @param mode Character scalar, one of \code{"auto"}, \code{"centroid"},
#'   \code{"point_on_surface"}, \code{"surface"}, \code{"line_midpoint"},
#'   or \code{"bbox_center"}. See \emph{Details}.
#' @param tmp_project Logical; if \code{TRUE} (default), temporarily project
#'   from geographic CRS to a suitable local projected CRS for line-based
#'   midpoints, returning results in the original CRS. Ignored when the input
#'   is already projected or when the selected \code{mode} does not require it.
#'
#' @return An \code{sf} object of the same row count and attributes as \code{x},
#'   with geometry coerced to POINTs and CRS preserved.
#'
#' @examples
#' \donttest{
#' library(sf)
#'
#' # Polygons -> points on surface
#' poly <- st_as_sf(data.frame(
#'   id = 1,
#'   wkt = "POLYGON((0 0, 2 0, 2 1, 0 1, 0 0))"
#' ), wkt = "wkt", crs = 4326)
#' coerce_to_points(poly, mode = "point_on_surface")
#'
#' # Line midpoint (auto chooses midpoint for LINESTRING)
#' ln <- st_as_sf(data.frame(
#'   id = 1,
#'   wkt = "LINESTRING(0 0, 1 1, 2 0)"
#' ), wkt = "wkt", crs = 4326)
#' coerce_to_points(ln, mode = "auto")
#' coerce_to_points(ln, mode = "line_midpoint", tmp_project = TRUE)
#'
#' # Bounding-box centers (fast labels)
#' coerce_to_points(poly, mode = "bbox_center")
#' }
#'
#' @seealso
#' \code{\link[sf]{st_centroid}}, \code{\link[sf]{st_point_on_surface}},
#' \code{\link[sf]{st_line_sample}}, \code{\link[sf]{st_transform}},
#' \code{\link[sf]{st_is_longlat}}, \code{ensure_projected()},
#' \code{.bbox_center_sfc()}
#'
#' @export
coerce_to_points <- function(
    x,
    mode = c("auto","centroid","point_on_surface","surface","line_midpoint","bbox_center"),
    tmp_project = TRUE
) {
  stopifnot(inherits(x, "sf"))
  mode <- match.arg(mode)
  if (identical(mode, "surface")) mode <- "point_on_surface"
  if (nrow(x) == 0L) return(x)
  
  g   <- sf::st_geometry(x)
  crs <- sf::st_crs(x)
  
  # -- bbox_center: numeric center of bbox (no polygon creation) ---------------
  if (mode == "bbox_center") {
    centers <- .bbox_center_sfc(x)  # defined elsewhere
    return(sf::st_set_geometry(x, centers))
  }
  
  # -- direct spherical-safe ops ----------------------------------------------
  if (mode == "centroid") {
    return(sf::st_set_geometry(x, sf::st_centroid(g)))
  }
  if (mode == "point_on_surface") {
    return(sf::st_set_geometry(x, sf::st_point_on_surface(g)))
  }
  
  # Helper: is this object in lon/lat?
  is_ll <- .is_longlat(x)  # defined elsewhere
  
  # Internal: coerce a single geometry (sfg or length-1 sfc) to a 1-row sf with CRS
  .as_sf_single <- function(geom, crs) {
    sfc1 <- if (inherits(geom, "sfc")) {
      if (is.na(sf::st_crs(geom)) && !is.na(crs)) sf::st_crs(geom) <- crs
      geom
    } else {
      sf::st_sfc(geom, crs = crs)
    }
    sf::st_sf(geometry = sfc1)
  }
  
  # Helper: midpoint for a single LINESTRING in the CRS of `x`, projecting if needed
  midpoint_linestring <- function(geom) {
    if (is_ll && !tmp_project) {
      # safe spherical fallback
      return(sf::st_centroid(geom))
    } else {
      x1    <- .as_sf_single(geom, crs)
      x_proj <- if (tmp_project) ensure_projected(x1) else x1
      midp   <- sf::st_line_sample(sf::st_geometry(x_proj), sample = 0.5)
      midp   <- sf::st_cast(midp, "POINT")
      # transform back to original CRS if needed
      if (!identical(sf::st_crs(x_proj), crs)) midp <- sf::st_transform(midp, crs)
      return(sf::st_geometry(midp))
    }
  }
  
  if (mode == "line_midpoint") {
    # LINESTRING only; MULTILINESTRING must be cast by the caller (UAT expects an error)
    gtypes <- as.character(sf::st_geometry_type(g, by_geometry = TRUE))
    if (any(gtypes %in% c("MULTILINESTRING","GEOMETRYCOLLECTION"))) {
      stop("line_midpoint only supports LINESTRING; cast MULTILINESTRING to LINESTRING first.")
    }
    out <- vector("list", length(g))
    for (i in seq_along(g)) {
      gt <- as.character(sf::st_geometry_type(g[i], by_geometry = TRUE))
      if (gt == "LINESTRING") {
        out[[i]] <- midpoint_linestring(g[i])[[1]]
      } else {
        # non-lines: use centroid as a reasonable fallback
        out[[i]] <- sf::st_centroid(g[i])[[1]]
      }
    }
    return(sf::st_set_geometry(x, sf::st_sfc(out, crs = crs)))
  }
  
  # ---------------------- mode == "auto" --------------------------------------
  out <- vector("list", length(g))
  for (i in seq_along(g)) {
    gt <- as.character(sf::st_geometry_type(g[i], by_geometry = TRUE))
    
    if (gt == "POINT") {
      out[[i]] <- g[[i]]  # pass-through
      
    } else if (gt == "MULTIPOINT") {
      # choose first coordinate
      coords <- sf::st_coordinates(sf::st_cast(g[i], "POINT"))
      out[[i]] <- sf::st_point(coords[1, 1:2])
      
    } else if (gt == "LINESTRING") {
      # midpoint (planar) or spherical centroid fallback
      out[[i]] <- midpoint_linestring(g[i])[[1]]
      
    } else if (gt == "MULTILINESTRING") {
      if (is_ll && !tmp_project) {
        # spherical-safe fallback
        out[[i]] <- sf::st_centroid(g[i])[[1]]
      } else {
        # pick the longest segment, then midpoint
        x1     <- .as_sf_single(g[i], crs)
        x_proj <- if (tmp_project) ensure_projected(x1) else x1
        parts  <- suppressWarnings(sf::st_cast(sf::st_geometry(x_proj), "LINESTRING"))
        if (length(parts) == 0L) {
          out[[i]] <- sf::st_centroid(g[i])[[1]]
        } else {
          lens <- as.numeric(sf::st_length(parts))
          j    <- if (length(lens)) which.max(lens) else 1L
          mp   <- sf::st_line_sample(parts[j], sample = 0.5)
          mp   <- sf::st_cast(mp, "POINT")
          # back to original CRS if needed
          if (!identical(sf::st_crs(x_proj), crs)) mp <- sf::st_transform(mp, crs)
          out[[i]] <- mp[[1]]
        }
      }
      
    } else if (gt %in% c("POLYGON","MULTIPOLYGON")) {
      # robust point guaranteed on the surface
      out[[i]] <- sf::st_point_on_surface(g[i])[[1]]
      
    } else {
      # fallback for other types: centroid
      out[[i]] <- sf::st_centroid(g[i])[[1]]
    }
  }
  
  sf::st_set_geometry(x, sf::st_sfc(out, crs = crs))
}

#' Generate seed points for Voronoi tessellation
#'
#' Creates an \pkg{sf} POINT layer of "seed" locations to drive a Voronoi
#' tessellation (or any cell-based partition). Multiple strategies are
#' supported: user-provided points, uniform random sampling within a boundary,
#' or k-means clustering of a sampling cloud.
#'
#' @details
#' The \code{method} determines how seeds are produced:
#'
#' \describe{
#'   \item{\code{"provided"}}{
#'     Use the supplied \code{seeds} (must be an \code{sf} POINT layer).
#'     The argument \code{n} is ignored. If \code{boundary} is supplied and
#'     both data have defined (and different) CRSs, \code{seeds} are transformed
#'     into \code{boundary}'s CRS.
#'   }
#'   \item{\code{"random"}}{
#'     Draw \code{n} points uniformly at random inside \code{boundary}.
#'     The boundary is first validated and unioned for robust sampling.
#'     When available, \code{sf::st_sample(..., exact = TRUE)} is used to obtain
#'     exactly \code{n} points; otherwise a fallback truncates/resamples to match
#'     \code{n}. Requires a polygonal \code{boundary}.
#'   }
#'   \item{\code{"kmeans"}}{
#'     Cluster a cloud of points with k-means and use the \code{n} cluster
#'     centers as seeds. If \code{sample_points} is provided (POINT \code{sf}),
#'     it is used as the cloud; otherwise a uniform cloud is pre-sampled inside
#'     \code{boundary} (\code{max(2000, 50 * n)} points). Requires a polygonal
#'     \code{boundary} unless \code{sample_points} is given. The clustering
#'     is performed on raw coordinates (Euclidean k-means); for meaningful
#'     distances, supply inputs in an appropriate projected CRS.
#'   }
#' }
#'
#' \strong{CRS handling.}
#' If \code{boundary} is supplied and has a defined CRS, outputs are returned
#' in the boundary's CRS. Inputs with defined but different CRSs are transformed
#' to match the boundary before processing. If \code{boundary} is \code{NULL},
#' outputs retain the CRS of the provided inputs (\code{seeds} or
#' \code{sample_points}/pre-sample cloud).
#'
#' \strong{Reproducibility.}
#' If \code{set_seed} is provided, \code{set.seed()} is called for reproducible
#' random sampling and k-means initialization; the prior RNG state is restored
#' on exit.
#'
#' @param boundary Optional polygonal \code{sf} object (\code{POLYGON} or
#'   \code{MULTIPOLYGON}) that defines the sampling/working area. Required for
#'   \code{method = "random"} and for \code{method = "kmeans"} unless
#'   \code{sample_points} is supplied.
#' @param method Character; one of \code{"kmeans"}, \code{"random"}, or
#'   \code{"provided"}. See \emph{Details}.
#' @param n Integer; number of seeds to return. Required for
#'   \code{method = "random"} and \code{method = "kmeans"}. Ignored for
#'   \code{method = "provided"}.
#' @param seeds \code{sf} POINT object of user-provided seeds; required when
#'   \code{method = "provided"}. All rows are used, and a sequential
#'   \code{seed_id} is added.
#' @param sample_points Optional \code{sf} POINT object providing the cloud to
#'   be clustered when \code{method = "kmeans"}. If provided together with a
#'   polygonal \code{boundary} and both have defined (different) CRSs,
#'   \code{sample_points} are transformed to the boundary CRS.
#' @param kmeans_nstart Integer; passed to \code{stats::kmeans()} \code{nstart}.
#'   Default \code{10}.
#' @param kmeans_iter Integer; passed to \code{stats::kmeans()} \code{iter.max}.
#'   Default \code{100}.
#' @param set_seed Optional integer to seed R's RNG for reproducibility.
#'
#' @return An \code{sf} POINT object with:
#' \itemize{
#'   \item \code{seed_id}: integer sequence \code{1..nrow(output)}.
#'   \item \code{method}: character indicating the generation method
#'         (\code{"provided"}, \code{"random"}, or \code{"kmeans"}).
#'   \item \code{geometry}: seed points; in \code{boundary}'s CRS when provided,
#'         otherwise in the CRS of the input used to generate seeds.
#' }
#'
#' @section Notes:
#' \itemize{
#'   \item For \code{"kmeans"}, clustering occurs in raw coordinate space;
#'         supply projected coordinates to avoid distortion in geographic CRSs.
#'   \item \code{boundary} is validated (\code{sf::st_make_valid()}) and unioned
#'         (\code{sf::st_union()}) prior to sampling for robustness.
#'   \item The function does not space seeds deterministically; use k-means for
#'         more even coverage relative to \code{sample_points} or the uniform
#'         cloud.
#' }
#'
#' @examples
#' \donttest{
#' library(sf)
#'
#' # Make a simple rectangular boundary
#' b <- st_as_sf(data.frame(
#'   id = 1,
#'   wkt = "POLYGON((0 0, 10 0, 10 5, 0 5, 0 0))"
#' ), wkt = "wkt", crs = 3857)
#'
#' # Random seeds
#' s_rand <- get_voronoi_seeds(boundary = b, method = "random", n = 50, set_seed = 1)
#'
#' # K-means seeds (auto-sampled cloud inside boundary)
#' s_km <- get_voronoi_seeds(boundary = b, method = "kmeans", n = 20, set_seed = 1)
#'
#' # Provided seeds (ignore n)
#' pts <- st_as_sf(data.frame(x = c(2, 5, 8), y = c(1, 3, 4)), coords = c("x","y"), crs = 3857)
#' s_prov <- get_voronoi_seeds(method = "provided", seeds = pts)
#'
#' # K-means using an external sample cloud (e.g., observations)
#' cloud <- st_as_sf(data.frame(
#'   x = runif(200, 0, 10),
#'   y = runif(200, 0, 5)
#' ), coords = c("x","y"), crs = 3857)
#' s_km2 <- get_voronoi_seeds(method = "kmeans", n = 10, sample_points = cloud, set_seed = 42)
#' }
#'
#' @seealso
#' \code{\link{create_voronoi_polygons}}, \code{\link{build_tessellation}},
#' \code{\link[sf]{st_sample}}, \code{\link[stats]{kmeans}}
#'
#' @importFrom sf st_make_valid st_union st_sample st_transform st_crs st_sf
#' @importFrom sf st_sfc st_point st_coordinates
#' @importFrom stats kmeans
#' @export
get_voronoi_seeds <- function(boundary = NULL,
                              method = c("kmeans", "random", "provided"),
                              n = NULL,
                              seeds = NULL,
                              sample_points = NULL,
                              kmeans_nstart = 10,
                              kmeans_iter = 100,
                              set_seed = NULL) {
  method <- match.arg(method)
  
  # -- small helpers -----------------------------------------------------------
  assert_sf <- function(x, what = c("POINT", "POLYGON", "MULTIPOLYGON")) {
    if (!inherits(x, "sf")) stop("Expected an 'sf' object.", call. = FALSE)
    gcls <- unique(sf::st_geometry_type(x, by_geometry = TRUE))
    if (!any(gcls %in% what)) {
      stop("Geometry must be one of: ", paste(what, collapse = ", "), call. = FALSE)
    }
  }
  
  same_or_transform_to_boundary <- function(x) {
    # If boundary is present and both CRSs are defined & different, transform x -> boundary CRS
    if (!is.null(boundary) &&
        !is.na(sf::st_crs(boundary)) &&
        !is.na(sf::st_crs(x)) &&
        !(sf::st_crs(boundary) == sf::st_crs(x))) {
      x <- sf::st_transform(x, sf::st_crs(boundary))
    }
    x
  }
  
  boundary_union <- function(b) {
    # Valid + union for sampling robustness
    b <- sf::st_make_valid(b)
    # suppress warnings about attributes during union
    sf::st_union(b)
  }
  
  # -- validate inputs ---------------------------------------------------------
  if (!is.null(boundary)) assert_sf(boundary, c("POLYGON", "MULTIPOLYGON"))
  
  if (method %in% c("random", "kmeans") && is.null(n)) {
    stop("Argument 'n' is required for method = '", method, "'.", call. = FALSE)
  }
  
  if (method == "provided") {
    if (is.null(seeds)) stop("Provide 'seeds' (sf POINT) for method = 'provided'.", call. = FALSE)
    assert_sf(seeds, "POINT")
  }
  
  if (!is.null(sample_points)) {
    assert_sf(sample_points, "POINT")
  }
  
  if (!is.null(set_seed)) {
    old_seed <- .Random.seed
    on.exit({ if (!is.null(old_seed)) .Random.seed <<- old_seed }, add = TRUE)
    set.seed(set_seed)
  }
  
  # -- main branches -----------------------------------------------------------
  out <- switch(
    method,
    
    "provided" = {
      s <- seeds
      # Transform seeds to boundary CRS when boundary is supplied
      s <- same_or_transform_to_boundary(s)
      s$seed_id <- seq_len(nrow(s))
      s$method  <- "provided"
      s
    },
    
    "random" = {
      if (is.null(boundary)) stop("Random seeds require 'boundary'.", call. = FALSE)
      
      b <- boundary_union(boundary)
      # st_sample exact=TRUE (if available) ensures length n
      pts <- try(sf::st_sample(b, size = n, type = "random", exact = TRUE), silent = TRUE)
      if (inherits(pts, "try-error")) {
        # Fallback for older sf without exact=
        pts <- sf::st_sample(b, size = n, type = "random")
        # If sample overshoots/undershoots, truncate or resample lightly
        if (length(pts) != n) {
          pts <- pts[seq_len(min(length(pts), n))]
          while (length(pts) < n) {
            pts_more <- sf::st_sample(b, size = n - length(pts), type = "random")
            pts <- c(pts, pts_more)
          }
        }
      }
      
      s <- sf::st_sf(seed_id = seq_len(n), method = "random", geometry = pts, crs = sf::st_crs(boundary))
      s
    },
    
    "kmeans" = {
      # Determine the cloud of points to cluster
      if (!is.null(sample_points)) {
        # Use provided sample points (transform to boundary if needed)
        if (!is.null(boundary)) sample_points <- same_or_transform_to_boundary(sample_points)
        cloud <- sample_points
      } else {
        if (is.null(boundary)) {
          stop("K-means seeds require 'boundary' (or provide 'sample_points').", call. = FALSE)
        }
        # Pre-sample a reasonably dense uniform cloud in boundary, then cluster
        b <- boundary_union(boundary)
        cloud_n <- max(2000L, 50L * as.integer(n))
        cloud_geom <- try(sf::st_sample(b, size = cloud_n, type = "random", exact = TRUE), silent = TRUE)
        if (inherits(cloud_geom, "try-error")) {
          cloud_geom <- sf::st_sample(b, size = cloud_n, type = "random")
        }
        cloud <- sf::st_sf(geometry = cloud_geom, crs = sf::st_crs(boundary))
      }
      
      xy <- sf::st_coordinates(cloud)
      km <- stats::kmeans(x = xy, centers = n, iter.max = kmeans_iter, nstart = kmeans_nstart)
      centers <- km$centers
      centers_sfc <- sf::st_sfc(
        lapply(seq_len(nrow(centers)), function(i) sf::st_point(centers[i, ])),
        crs = sf::st_crs(cloud)
      )
      
      # Make sure seeds are in boundary CRS if boundary provided
      s <- sf::st_sf(seed_id = seq_len(n), method = "kmeans", geometry = centers_sfc)
      s <- same_or_transform_to_boundary(s)
      s
    }
  )
  
  # Ensure exactly POINT geometries and desired CRS if boundary exists
  if (!is.null(boundary) &&
      !is.na(sf::st_crs(boundary)) &&
      !(sf::st_crs(out) == sf::st_crs(boundary))) {
    out <- sf::st_transform(out, sf::st_crs(boundary))
  }
  
  out
}

#' Build a polygonal clip target from points and/or a boundary
#'
#' Constructs a polygonal \code{sf} layer suitable for clipping tessellations or
#' other geometry operations. If a \code{boundary} is provided, it is validated,
#' optionally expanded (buffered), and returned (in the points' CRS). If no
#' \code{boundary} is provided, the function creates an expanded bounding box
#' around the input points; degenerate cases (zero-width/height bbox) are
#' handled by buffering the point union.
#'
#' @details
#' The resulting polygon is intended as a robust "target extent" for clipping.
#'
#' \strong{CRS handling}
#' \itemize{
#'   \item If \code{boundary} is supplied and both \code{points_sf} and
#'         \code{boundary} have defined CRSs that differ, \code{boundary} is
#'         transformed to the CRS of \code{points_sf} before further steps.
#'   \item The returned object always uses the CRS of \code{points_sf}.
#' }
#'
#' \strong{Geometry validation}
#' \itemize{
#'   \item \code{boundary} (when supplied) is validated using
#'         \code{sf::st_make_valid()} if available (or a safe fallback).
#'   \item The constructed target (buffered boundary or bbox/buffer) is also
#'         validated before return.
#' }
#'
#' \strong{Expansion distance}
#' \itemize{
#'   \item If \code{expand} is strictly between \code{0} and \code{1}, it is
#'         interpreted as a fraction of the maximum side of the reference
#'         geometry's bounding box (the boundary if present, otherwise the
#'         points' bbox). The resulting absolute distance is used to buffer.
#'   \item If \code{expand} is \code{0}, no expansion occurs.
#'   \item If \code{expand} is \code{>= 1} or negative, it is treated as an
#'         absolute buffer distance in the CRS units (negative values shrink).
#'   \item When working in geographic (lon/lat) CRSs, absolute distances are in
#'         degrees; for distance-aware buffering, supply/project to a suitable
#'         planar CRS beforehand.
#' }
#'
#' \strong{Degenerate bounding boxes}
#' \itemize{
#'   \item If no \code{boundary} is given and the points' bbox has zero width or
#'         height (e.g., identical or collinear points), a small buffer around
#'         the union of points is used to form a polygonal target. A message is
#'         printed unless \code{quiet = TRUE}.
#' }
#'
#' @param points_sf An \code{sf} object whose geometry must be \code{POINT} or
#'   \code{MULTIPOINT}. Used for CRS reference and, when \code{boundary} is
#'   \code{NULL}, to derive the bounding box for the clip target.
#' @param boundary Optional polygonal \code{sf} object (\code{POLYGON} or
#'   \code{MULTIPOLYGON}). If supplied, it defines the clip target (after
#'   validation and optional expansion). Must be polygonal if not \code{NULL}.
#' @param expand Numeric. Expansion (buffer) applied to the reference geometry:
#'   \itemize{
#'     \item \code{0 < expand < 1}: fraction of max(bbox width, bbox height).
#'     \item \code{expand == 0}: no expansion.
#'     \item \code{expand >= 1} or negative: absolute distance in CRS units
#'           (negative shrinks).
#'   }
#' @param quiet Logical; if \code{FALSE} (default), emits informational messages
#'   for degenerate-bbox handling and other notable branches.
#'
#' @return An \code{sf} polygon layer (single-column \code{geometry}) representing
#'   the clip target. The CRS matches that of \code{points_sf}. When
#'   \code{boundary} is provided, the geometry is the (optionally buffered)
#'   boundary in the points' CRS; otherwise it is an expanded bbox (or a small
#'   buffer around the points if the bbox is degenerate).
#'
#' @examples
#' \donttest{
#' library(sf)
#'
#' # Example points in a projected CRS
#' pts <- st_as_sf(
#'   data.frame(x = c(100, 200, 250), y = c(50, 120, 90)),
#'   coords = c("x","y"), crs = 3857
#' )
#'
#' # 1) No boundary: build an expanded bbox target (10% of bbox max side)
#' tgt1 <- clip_target_for(points_sf = pts, boundary = NULL, expand = 0.1)
#'
#' # 2) With a polygonal boundary and absolute expansion of 500 meters (CRS units)
#' poly <- st_as_sf(
#'   data.frame(
#'     wkt = "POLYGON((0 0, 400 0, 400 200, 0 200, 0 0))"
#'   ), wkt = "wkt", crs = 3857
#' )
#' tgt2 <- clip_target_for(points_sf = pts, boundary = poly, expand = 500)
#'
#' # 3) Degenerate bbox (identical points): small buffer is used
#' pts_same <- st_as_sf(data.frame(x = 0, y = 0), coords = c("x","y"), crs = 3857)
#' tgt3 <- clip_target_for(points_sf = pts_same, expand = 0.05)
#' }
#'
#' @seealso
#' \code{\link{build_tessellation}},
#' \code{\link[sf]{st_make_valid}},
#' \code{\link[sf]{st_buffer}},
#' \code{\link[sf]{st_as_sfc}},
#' \code{\link[sf]{st_bbox}}
#'
#' @importFrom sf st_geometry_type st_crs st_transform st_make_valid
#' @importFrom sf st_buffer st_bbox st_as_sfc st_union st_geometry st_sf
#' @export
clip_target_for <- function(
    points_sf,
    boundary = NULL,
    expand   = 0,
    quiet    = FALSE
) {
  # ---- debug helper ----
  .msg <- function(...) if (!quiet) message(...)
  
  # ---- safety ----
  if (!inherits(points_sf, "sf")) {
    stop("clip_target_for(): `points_sf` must be an sf object.")
  }
  gtypes <- as.character(sf::st_geometry_type(points_sf, by_geometry = TRUE))
  if (!all(gtypes %in% c("POINT", "MULTIPOINT"))) {
    stop("clip_target_for(): `points_sf` geometry must be POINT/MULTIPOINT.")
  }
  
  # ---- helpers ----
  .safe_make_valid <- function(g) {
    if ("st_make_valid" %in% getNamespaceExports("sf")) {
      return(suppressWarnings(sf::st_make_valid(g)))
    }
    if (requireNamespace("lwgeom", quietly = TRUE) &&
        "st_make_valid" %in% getNamespaceExports("lwgeom")) {
      return(suppressWarnings(lwgeom::st_make_valid(g)))
    }
    suppressWarnings(sf::st_buffer(g, 0))
  }
  
  .align_crs <- function(a, b) {
    if (is.null(a) || is.null(b)) return(a)
    if (is.na(sf::st_crs(a)) || is.na(sf::st_crs(b))) return(a)
    if (sf::st_crs(a) == sf::st_crs(b)) return(a)
    sf::st_transform(a, sf::st_crs(b))
  }
  
  .expand_distance <- function(ref_geom, expand) {
    # If 0 < expand < 1, interpret as fraction of max bbox side.
    # Otherwise, treat as absolute distance in CRS units.
    if (!is.numeric(expand) || length(expand) != 1 || is.na(expand) || expand == 0) return(0)
    bb <- sf::st_bbox(ref_geom)
    dx <- as.numeric(bb$xmax - bb$xmin)
    dy <- as.numeric(bb$ymax - bb$ymin)
    if (expand > 0 && expand < 1) {
      return(max(dx, dy) * expand)
    }
    expand
  }
  
  # ---- build target ----
  crs_pts <- sf::st_crs(points_sf)
  
  if (!is.null(boundary)) {
    # Use provided boundary (aligned, valid), with optional expansion
    boundary <- .align_crs(boundary, points_sf)
    
    if (!any(sf::st_geometry_type(boundary) %in% c("POLYGON", "MULTIPOLYGON"))) {
      stop("clip_target_for(): `boundary` must be polygonal.")
    }
    
    boundary <- .safe_make_valid(boundary)
    dist <- .expand_distance(boundary, expand)
    tgt <- if (dist != 0) {
      suppressWarnings(sf::st_buffer(boundary, dist = dist))
    } else {
      boundary
    }
    return(sf::st_sf(geometry = sf::st_geometry(tgt)))
  }
  
  # No boundary: use an expanded bbox (or a small buffer if degenerate)
  pts_geom <- sf::st_geometry(points_sf)
  if (length(pts_geom) == 0) stop("clip_target_for(): `points_sf` is empty.")
  
  bb <- sf::st_bbox(pts_geom)
  # Handle degenerate bbox (all points identical or collinear with zero width/height)
  zero_w <- isTRUE(all.equal(as.numeric(bb$xmin), as.numeric(bb$xmax)))
  zero_h <- isTRUE(all.equal(as.numeric(bb$ymin), as.numeric(bb$ymax)))
  
  if (zero_w || zero_h) {
    # Create a minimal envelope around the points using a tiny buffer,
    # scaled by a reasonable fallback distance.
    .msg("clip_target_for(): degenerate bbox; using small buffer around points.")
    # Fallback distance: 1% of a nominal unit or 1 if uncomputable
    fallback <- 1
    dist0 <- .expand_distance(points_sf, if (expand == 0) 0.05 else expand)
    dist_use <- if (dist0 > 0) dist0 else fallback
    tgt <- suppressWarnings(sf::st_buffer(sf::st_union(pts_geom), dist = dist_use))
    tgt <- .safe_make_valid(tgt)
    sf::st_crs(tgt) <- crs_pts
    return(sf::st_sf(geometry = sf::st_geometry(tgt)))
  }
  
  # Regular bbox
  tgt <- sf::st_as_sfc(bb, crs = crs_pts)
  dist <- .expand_distance(tgt, expand)
  if (dist != 0) {
    tgt <- suppressWarnings(sf::st_buffer(tgt, dist = dist))
  }
  tgt <- .safe_make_valid(tgt)
  sf::st_sf(geometry = tgt)
}

# -----------------------------------------------------------------------------
# Voronoi Tessellation
# -----------------------------------------------------------------------------

#' Create Voronoi polygons from points with robust CRS and optional clipping
#'
#' Builds a planar Voronoi tessellation from input point geometries and
#' optionally clips the result to a polygonal boundary. Handles geographic
#' (lon/lat) inputs by projecting to a suitable planar CRS, validates polygon
#' geometries, and (optionally) deduplicates coincident points for graph
#' construction without losing the original point–cell mapping.
#'
#' @details
#' \strong{Inputs and geometry types}
#' \itemize{
#'   \item \code{points_sf} must be an \code{sf} object whose geometries are
#'         \code{POINT} or \code{MULTIPOINT}. At least one row is required.
#'   \item \code{boundary}, when supplied, must be polygonal
#'         (\code{POLYGON}/\code{MULTIPOLYGON}). If omitted, a convex hull of
#'         \code{points_sf} is used as the working boundary (a message is
#'         printed unless \code{quiet = TRUE}).
#' }
#'
#' \strong{CRS handling}
#' \itemize{
#'   \item If \code{crs} is provided, both \code{points_sf} and \code{boundary}
#'         (when present) are transformed to that CRS before tessellation.
#'   \item If \code{crs} is \code{NULL} and \code{points_sf} is in geographic
#'         coordinates, the function projects to a local planar CRS
#'         (prefers user-defined \code{ensure_projected()}, otherwise uses a UTM
#'         guess). The \code{boundary} (if present) is aligned to the points'
#'         CRS.
#'   \item The returned \code{cells} always use the CRS of the transformed
#'         \code{points_sf}.
#' }
#'
#' \strong{Boundary preparation and expansion}
#' \itemize{
#'   \item Polygon inputs are validated using \code{lwgeom::st_make_valid()}
#'         when available; otherwise a \code{st_buffer(., 0)} fallback is used.
#'   \item \code{expand} is interpreted as an absolute buffer distance in CRS
#'         units. When \code{expand > 0}, the boundary is buffered outward by
#'         \code{expand} and that expanded boundary is used as the Voronoi
#'         envelope (and for clipping when \code{clip = TRUE}).
#' }
#'
#' \strong{Duplicates and stable mapping}
#' \itemize{
#'   \item When \code{keep_duplicates = FALSE} (default), coincident points are
#'         removed \emph{only} for the Voronoi graph construction to avoid
#'         degeneracies. The final point–cell index is still computed for every
#'         original input row.
#'   \item Each input point row is assigned to exactly one cell via spatial
#'         intersection; ties (points lying on shared boundaries) are resolved
#'         deterministically by choosing the smallest \code{cell_id}.
#'         If a point does not intersect any cell after clipping, the nearest
#'         cell is used.
#' }
#'
#' \strong{Cell identifiers}
#' \itemize{
#'   \item If a user helper \code{ensure_stable_poly_id()} is available, it is
#'         used to generate a stable \code{cell_id}; otherwise sequential
#'         identifiers \code{1..n} are assigned.
#' }
#'
#' @param points_sf An \code{sf} object with \code{POINT}/\code{MULTIPOINT}
#'   geometries from which the Voronoi diagram is constructed.
#' @param boundary Optional polygonal \code{sf} object used to constrain the
#'   tessellation. If \code{NULL}, a convex hull of \code{points_sf} is used.
#' @param expand Numeric; absolute buffer distance (in CRS units) applied to the
#'   \code{boundary} to enlarge the Voronoi envelope (and, if \code{clip = TRUE},
#'   the clipping geometry). Use \code{0} for no expansion.
#' @param clip Logical; if \code{TRUE} (default), intersect Voronoi cells with
#'   the (expanded) \code{boundary}. If \code{FALSE}, the un-clipped Voronoi
#'   polygons within the envelope are returned.
#' @param keep_duplicates Logical; if \code{FALSE} (default), coincident points
#'   are deduplicated for the Voronoi graph construction only; the returned
#'   index still has one entry per input row.
#' @param crs Optional target CRS (any form accepted by \code{sf::st_crs()}) to
#'   which inputs are transformed prior to tessellation. If \code{NULL}, the
#'   current CRS is used; geographic inputs are projected automatically.
#' @param quiet Logical; suppresses informational messages when \code{TRUE}.
#'
#' @return A named \code{list} with components:
#' \describe{
#'   \item{\code{cells}}{An \code{sf} polygon layer of Voronoi cells with a
#'     \code{cell_id} column. CRS matches the (possibly transformed) input
#'     points.}
#'   \item{\code{index}}{An integer vector of length \code{nrow(points_sf)} that
#'     maps each input point row to the row number of its associated cell in
#'     \code{cells}.}
#'   \item{\code{boundary}}{The polygonal boundary actually used (possibly the
#'     convex hull of points and/or buffered), in the cells' CRS.}
#'   \item{\code{method}}{The string \code{"voronoi"}.}
#'   \item{\code{params}}{A list of parameters used: \code{clip},
#'     \code{expand}, \code{keep_duplicates}.}
#' }
#'
#' @examples
#' \donttest{
#' library(sf)
#'
#' # Example points (Web Mercator)
#' pts <- st_as_sf(data.frame(x = c(0, 100, 200, 150),
#'                            y = c(0,  50, 100, 150)),
#'                 coords = c("x","y"), crs = 3857)
#'
#' # Simple Voronoi clipped to a rectangular boundary, expanded by 50 units
#' b <- st_as_sf(
#'   data.frame(wkt = "POLYGON(( -50 -50, 300 -50, 300 250, -50 250, -50 -50 ))"),
#'   wkt = "wkt", crs = 3857
#' )
#' res <- create_voronoi_polygons(pts, boundary = b, expand = 50, clip = TRUE)
#' plot(st_geometry(res$boundary), col = NA, border = 'grey')
#' plot(st_geometry(res$cells), add = TRUE)
#' points(st_coordinates(pts))
#'
#' # Geographic input: automatic projection and convex-hull boundary
#' pts_ll <- st_transform(pts, 4326)
#' res_ll <- create_voronoi_polygons(pts_ll, boundary = NULL, clip = TRUE)
#' }
#'
#' @seealso
#' \code{\link{build_tessellation}},
#' \code{\link{clip_target_for}},
#' \code{\link[sf]{st_voronoi}},
#' \code{\link[sf]{st_intersection}},
#' \code{\link[sf]{st_make_valid}},
#' \code{\link[lwgeom]{st_make_valid}}
#'
#' @importFrom sf st_geometry_type st_crs st_transform st_union st_convex_hull
#' @importFrom sf st_voronoi st_collection_extract st_sf st_buffer st_intersection
#' @importFrom sf st_is_empty st_nearest_feature st_geometry
#' @importFrom lwgeom st_make_valid
#' @export
create_voronoi_polygons <- function(
    points_sf,
    boundary        = NULL,
    expand          = 0,
    clip            = TRUE,
    keep_duplicates = FALSE,
    crs             = NULL,
    quiet           = FALSE
) {
  if (!inherits(points_sf, "sf")) {
    stop("create_voronoi_polygons(): `points_sf` must be an sf object.")
  }
  gtypes <- as.character(sf::st_geometry_type(points_sf, by_geometry = TRUE))
  if (!all(gtypes %in% c("POINT", "MULTIPOINT"))) {
    stop("create_voronoi_polygons(): `points_sf` geometry must be POINT/MULTIPOINT.")
  }
  if (nrow(points_sf) < 1) {
    stop("create_voronoi_polygons(): `points_sf` has no rows.")
  }
  
  # ---- helpers ---------------------------------------------------------------
  .safe_make_valid <- function(g) {
    if (requireNamespace("lwgeom", quietly = TRUE) &&
        "st_make_valid" %in% getNamespaceExports("lwgeom")) {
      return(suppressWarnings(lwgeom::st_make_valid(g)))
    }
    suppressWarnings(sf::st_buffer(g, 0))
  }
  .is_longlat <- function(x) {
    out <- tryCatch(sf::st_is_longlat(x), error = function(e) NA)
    isTRUE(out)
  }
  .guess_utm <- function(x) {
    ctr <- try(sf::st_coordinates(sf::st_centroid(sf::st_union(sf::st_geometry(x)))), silent = TRUE)
    if (inherits(ctr, "try-error") || length(ctr) < 2) return(4326)
    lon <- as.numeric(ctr[1]); lat <- as.numeric(ctr[2])
    zone <- floor((lon + 180) / 6) + 1
    if (lat >= 0) 32600 + zone else 32700 + zone
  }
  .ensure_projected <- function(x) {
    if (exists("ensure_projected", mode = "function")) return(ensure_projected(x))
    if (.is_longlat(x)) {
      epsg <- .guess_utm(x)
      if (!quiet) message(sprintf("create_voronoi_polygons(): projecting to EPSG:%s.", epsg))
      return(sf::st_transform(x, epsg))
    }
    x
  }
  .align_crs <- function(a, b) {
    if (is.null(b)) return(a)
    if (is.na(sf::st_crs(a)) || is.na(sf::st_crs(b))) return(a)
    if (sf::st_crs(a) == sf::st_crs(b)) return(a)
    sf::st_transform(a, sf::st_crs(b))
  }
  .dedupe_points <- function(g) {
    m <- sf::st_coordinates(g)
    key <- paste0(round(m[,1], 10), "_", round(m[,2], 10))
    keep <- !duplicated(key)
    structure(list(geom = g[keep, , drop = FALSE], keep_idx = keep, key = key), class = "dedup_pack")
  }
  
  # ---- CRS handling ----------------------------------------------------------
  pts <- points_sf
  if (!is.null(crs)) {
    pts <- sf::st_transform(pts, crs)
    if (!is.null(boundary)) boundary <- sf::st_transform(boundary, crs)
  } else {
    if (.is_longlat(pts)) pts <- .ensure_projected(pts)
    if (!is.null(boundary)) boundary <- .align_crs(boundary, pts)
  }
  
  # ---- boundary (derive if missing) ------------------------------------------
  if (is.null(boundary)) {
    if (!quiet) message("create_voronoi_polygons(): deriving boundary via convex hull of points.")
    boundary <- sf::st_convex_hull(sf::st_union(sf::st_geometry(pts)))
    boundary <- sf::st_sf(geometry = boundary, crs = sf::st_crs(pts))
  } else {
    if (!any(sf::st_geometry_type(boundary) %in% c("POLYGON", "MULTIPOLYGON"))) {
      stop("create_voronoi_polygons(): `boundary` must be polygonal.")
    }
  }
  boundary <- .safe_make_valid(boundary)
  
  # expanded boundary for envelope (optional)
  boundary_expanded <- if (isTRUE(is.numeric(expand)) && expand > 0) {
    suppressWarnings(sf::st_buffer(boundary, dist = expand))
  } else boundary
  
  # ---- Voronoi construction --------------------------------------------------
  geom_pts <- sf::st_geometry(pts)
  
  # optional dedupe for voronoi graph – doesn't affect mapping (we map from originals)
  pts_for_graph <- if (isTRUE(keep_duplicates)) geom_pts else .dedupe_points(geom_pts)$geom
  
  # st_voronoi expects a single MULTIPOINT / GEOMETRYCOLLECTION
  mp  <- sf::st_union(pts_for_graph)
  env <- sf::st_geometry(boundary_expanded)
  vor <- suppressWarnings(sf::st_voronoi(mp, envelope = env))
  cells <- sf::st_collection_extract(vor, "POLYGON", warn = FALSE)
  cells <- sf::st_sf(geometry = cells, crs = sf::st_crs(pts))
  cells <- .safe_make_valid(cells)
  
  if (clip) {
    clip_to <- if (isTRUE(is.numeric(expand)) && expand > 0) boundary_expanded else boundary
    cells <- suppressWarnings(sf::st_intersection(cells, clip_to))
    empt <- sf::st_is_empty(cells)
    if (any(empt, na.rm = TRUE)) cells <- cells[!empt, , drop = FALSE]
  }
  
  # stable cell ids
  if (nrow(cells) > 0) {
    if (exists("ensure_stable_poly_id", mode = "function")) {
      cells <- ensure_stable_poly_id(cells, id_col = "cell_id")
    } else {
      cells$cell_id <- seq_len(nrow(cells))
    }
  } else {
    cells$cell_id <- integer(0)
  }
  
  # ---- map points -> cells (duplicate-safe, deterministic) -------------------
  index <- rep(NA_integer_, nrow(pts))
  if (nrow(cells) > 0) {
    hits <- sf::st_intersects(pts, cells)
    # pre-extract the vector of cell_ids to allow deterministic min() tie-break
    cell_ids <- cells$cell_id
    for (i in seq_along(hits)) {
      ids <- as.integer(hits[[i]])
      if (length(ids) == 1L) {
        index[i] <- ids
      } else if (length(ids) >= 2L) {
        # deterministic choice if a point lies exactly on a boundary (multiple cells)
        # choose the minimum cell_id so duplicates pick the same one
        index[i] <- ids[which.min(cell_ids[ids])]
      } else {
        # rare: no intersection (e.g., exact-after-clip); fallback to nearest cell
        near <- suppressWarnings(sf::st_nearest_feature(pts[i, ], cells))
        if (is.finite(near)) index[i] <- as.integer(near)
      }
    }
  }
  
  # ---- return ---------------------------------------------------------------
  list(
    cells    = cells,
    index    = index,
    boundary = boundary,
    method   = "voronoi",
    params   = list(
      clip            = clip,
      expand          = expand,
      keep_duplicates = keep_duplicates
    )
  )
}

# -----------------------------------------------------------------------------
# Grid Tessellations (Hex / Square)
# -----------------------------------------------------------------------------

#' Create square or hexagonal grid polygons over a boundary (with robust CRS handling)
#'
#' Builds a planar grid of polygons (square or hexagonal) over the bounding box
#' of a polygonal \code{boundary}, and optionally clips the grid to that
#' boundary (respecting interior holes). Cell size and/or grid resolution can be
#' specified directly, or derived from a target number of cells. Inputs in
#' geographic coordinates are projected to a suitable local planar CRS to ensure
#' stable distances and areas.
#'
#' @details
#' \strong{Inputs and geometry types}
#' \itemize{
#'   \item \code{boundary} must be an \code{sf} or \code{sfc} object with
#'         polygonal geometry (\code{POLYGON}/\code{MULTIPOLYGON}). A positive
#'         bounding-box extent is required.
#' }
#'
#' \strong{CRS handling}
#' \itemize{
#'   \item If \code{crs} is supplied, \code{boundary} is transformed to that CRS
#'         before grid construction.
#'   \item If \code{crs} is \code{NULL} and \code{boundary} is geographic
#'         (lon/lat), the function projects to a local planar CRS (prefers a
#'         package helper \code{ensure_projected()}, if available; otherwise a
#'         UTM zone is guessed from the boundary centroid).
#'   \item The returned grid is in the CRS used to build it (either the
#'         user-specified \code{crs}, a projected CRS chosen automatically, or
#'         the original CRS if it was already projected).
#' }
#'
#' \strong{Grid sizing}
#' \itemize{
#'   \item Provide \code{cellsize} (scalar or length-2 \code{c(dx, dy)} in CRS
#'         units) to control polygon dimensions directly.
#'   \item Provide \code{n} (scalar or length-2 \code{c(nx, ny)}) to request a
#'         specific number of cells along the x and y directions; \code{cellsize}
#'         is then derived from the boundary bbox.
#'   \item If neither \code{cellsize} nor \code{n} is given, \code{target_cells}
#'         is used to derive a reasonable \code{n} and \code{cellsize} from the
#'         bbox aspect ratio.
#'   \item If both \code{cellsize} and \code{n} are provided, both are forwarded
#'         to \code{sf::st_make_grid()} (they should be consistent). If \code{n}
#'         is invalid (e.g., non-positive), it is ignored in favor of
#'         \code{cellsize}.
#' }
#'
#' \strong{Grid type and clipping}
#' \itemize{
#'   \item \code{type = "square"} produces a square grid; \code{type = "hex"}
#'         produces a hexagonal grid (via \code{square = FALSE} in
#'         \code{sf::st_make_grid()}).
#'   \item When \code{clip = TRUE} (default), the raw grid is intersected with
#'         the (validated) boundary, preserving interior holes. When
#'         \code{clip = FALSE}, the full grid covering the boundary's bbox is
#'         returned.
#' }
#'
#' \strong{Validity}
#' \itemize{
#'   \item The boundary is validated prior to clipping using
#'         \code{sf::st_make_valid()} when available (otherwise, a
#'         \code{st_buffer(., 0)} fallback is used) to reduce topology errors.
#' }
#'
#' @param boundary Polygonal \code{sf} or \code{sfc} object defining the area of
#'   interest; must be \code{POLYGON}/\code{MULTIPOLYGON}.
#' @param target_cells Optional positive number giving the approximate desired
#'   total number of cells when neither \code{cellsize} nor \code{n} is
#'   supplied. Used to derive a grid resolution consistent with the boundary
#'   bbox aspect.
#' @param type Grid type: \code{"square"} or \code{"hex"}.
#' @param cellsize Optional numeric cell size in CRS units; either a scalar
#'   (applied to both x and y) or a length-2 vector \code{c(dx, dy)}. Must be
#'   strictly positive.
#' @param n Optional grid resolution; an integer scalar (applied to both axes)
#'   or length-2 integer vector \code{c(nx, ny)} with values \eqn{\ge 1}.
#' @param clip Logical; if \code{TRUE} (default), intersect the grid with the
#'   boundary to keep only areas inside (holes are respected). If \code{FALSE},
#'   return the full un-clipped grid covering the boundary's bbox.
#' @param crs Optional target CRS (any form accepted by \code{sf::st_crs()}) to
#'   which the boundary is transformed before tessellation. If \code{NULL}, a
#'   local planar CRS is chosen automatically for geographic inputs.
#' @param quiet Logical; suppresses informational messages when \code{TRUE}.
#'
#' @return An \code{sf} polygon layer containing the grid cells, with a
#'   \code{poly_id} column (sequential ids). The layer’s CRS is the CRS used for
#'   grid construction.
#'
#' @examples
#' \donttest{
#' library(sf)
#'
#' # Simple rectangular boundary (Web Mercator)
#' b <- st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 1000, ymax = 600), crs = st_crs(3857)))
#' b <- st_sf(geometry = b)
#'
#' # Square grid from target number of cells
#' g_sq <- create_grid_polygons(boundary = b, target_cells = 30, type = "square")
#'
#' # Hex grid with explicit cellsize
#' g_hex <- create_grid_polygons(boundary = b, cellsize = 100, type = "hex")
#'
#' # Grid with explicit resolution (nx, ny), no clipping
#' g_nc <- create_grid_polygons(boundary = b, n = c(5, 3), type = "square", clip = FALSE)
#'
#' plot(st_geometry(g_sq), main = "Square grid (clipped)"); plot(st_geometry(b), add = TRUE)
#' }
#'
#' @seealso
#' \code{\link{build_tessellation}},
#' \code{\link{create_voronoi_polygons}},
#' \code{\link[sf]{st_make_grid}},
#' \code{\link[sf]{st_intersection}},
#' \code{\link[sf]{st_make_valid}}
#'
#' @importFrom sf st_geometry_type st_crs st_transform st_bbox st_as_sfc
#' @importFrom sf st_make_grid st_sf st_buffer st_intersection st_make_valid
#' @export
create_grid_polygons <- function(
    boundary,
    target_cells = NULL,         # approximate desired number of cells
    type         = c("square","hex"),
    cellsize     = NULL,         # numeric scalar or length-2 (x,y)
    n            = NULL,         # integer scalar or length-2 (nx,ny)
    clip         = TRUE,
    crs          = NULL,
    quiet        = FALSE
) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("create_grid_polygons(): package 'sf' is required.")
  }
  type <- match.arg(type)
  
  # ---- helpers ---------------------------------------------------------------
  .msg <- function(...) if (!quiet) message(...)
  .as_sf <- function(x) {
    if (inherits(x, "sf")) return(x)
    if (inherits(x, "sfc")) return(sf::st_sf(geometry = x))
    stop("create_grid_polygons(): 'boundary' must be an sf or sfc object.")
  }
  .is_poly <- function(x) {
    gtypes <- as.character(sf::st_geometry_type(x, by_geometry = TRUE))
    all(gtypes %in% c("POLYGON","MULTIPOLYGON"))
  }
  .make_valid <- function(g) {
    # Prefer sf::st_make_valid if exported; otherwise buffer(.,0) fallback
    if ("st_make_valid" %in% getNamespaceExports("sf")) {
      return(suppressWarnings(sf::st_make_valid(g)))
    }
    return(suppressWarnings(sf::st_buffer(g, 0)))
  }
  .ensure_projected_local <- function(x, tgt = NULL) {
    # Use your package helper if present
    if (exists("ensure_projected", mode = "function")) {
      return(ensure_projected(x, tgt))
    }
    # Minimal fallback: if long/lat -> guess UTM zone; else return as-is
    is_ll <- tryCatch(sf::st_is_longlat(x), error = function(e) NA)
    if (isTRUE(is_ll)) {
      ctr  <- sf::st_coordinates(sf::st_centroid(sf::st_union(sf::st_geometry(x))))
      lon  <- as.numeric(ctr[1]); lat <- as.numeric(ctr[2])
      zone <- floor((lon + 180)/6) + 1
      epsg <- if (lat >= 0) 32600 + zone else 32700 + zone
      return(sf::st_transform(x, epsg))
    }
    x
  }
  
  # ---- boundary prep ---------------------------------------------------------
  boundary <- .as_sf(boundary)
  .msg(sprintf("create_grid_polygons(): boundary class = [%s]",
               paste(class(boundary), collapse = ", ")))
  .msg(sprintf("create_grid_polygons(): boundary geom types = [%s]",
               paste(unique(as.character(sf::st_geometry_type(boundary))), collapse = ", ")))
  
  if (!.is_poly(boundary)) {
    stop("create_grid_polygons(): 'boundary' must be polygonal (POLYGON/MULTIPOLYGON).")
  }
  
  if (!is.null(crs)) {
    boundary <- sf::st_transform(boundary, crs)
  } else {
    boundary <- .ensure_projected_local(boundary)
  }
  boundary <- .make_valid(boundary)
  
  bb <- sf::st_bbox(boundary)
  w  <- as.numeric(bb["xmax"] - bb["xmin"])
  h  <- as.numeric(bb["ymax"] - bb["ymin"])
  if (!(is.finite(w) && is.finite(h) && w > 0 && h > 0)) {
    stop("create_grid_polygons(): boundary bbox has non-positive extent.")
  }
  .msg(sprintf("create_grid_polygons(): bbox width=%.6g height=%.6g (CRS units)", w, h))
  env <- sf::st_as_sfc(bb, crs = sf::st_crs(boundary))
  
  # ---- derive n and/or cellsize ----------------------------------------------
  if (!is.null(cellsize)) {
    if (length(cellsize) == 1L) cellsize <- rep(cellsize, 2L)
    if (length(cellsize) != 2L || any(!is.finite(cellsize)) || any(cellsize <= 0)) {
      stop("create_grid_polygons(): 'cellsize' must be positive numeric (length 1 or 2).")
    }
    if (!is.null(n)) {
      if (length(n) == 1L) n <- rep(as.integer(n), 2L)
      n <- as.integer(n[1:2])
      if (any(is.na(n)) || any(n < 1)) {
        .msg("create_grid_polygons(): ignoring invalid 'n' and using 'cellsize' only.")
        n <- NULL
      }
    }
  } else if (!is.null(n)) {
    if (length(n) == 1L) n <- rep(as.integer(n), 2L)
    n <- as.integer(n[1:2])
    if (any(is.na(n)) || any(n < 1)) {
      stop("create_grid_polygons(): 'n' must be integer (length 1 or 2) and >= 1.")
    }
    cellsize <- c(w / n[1], h / n[2])
  } else {
    if (is.null(target_cells) || !is.finite(target_cells) || target_cells < 1) {
      stop("create_grid_polygons(): supply either 'cellsize', 'n', or a positive 'target_cells'.")
    }
    nx <- max(1L, round(sqrt(target_cells * (w / h))))
    ny <- max(1L, round(ceiling(target_cells / nx)))
    n  <- c(nx, ny)
    cellsize <- c(w / nx, h / ny)
    .msg(sprintf("create_grid_polygons(): derived n=c(%d,%d); cellsize=%.6g x %.6g (CRS units)",
                 nx, ny, cellsize[1], cellsize[2]))
  }
  .msg(sprintf("create_grid_polygons(): type=%s; using cellsize=%.6g x %.6g; n=%s",
               type, cellsize[1], cellsize[2],
               if (is.null(n)) "NULL" else paste0("c(", paste(n, collapse=","), ")")))
  
  # ---- build grid ------------------------------------------------------------
  grid_sfc <- sf::st_make_grid(
    env,
    cellsize = cellsize,
    n        = n,                 # may be NULL; that's fine when cellsize is given
    what     = "polygons",
    square   = identical(type, "square")
  )
  raw_n <- length(grid_sfc)
  .msg(sprintf("create_grid_polygons(): raw grid cells = %d", raw_n))
  if (raw_n == 0L) {
    stop("create_grid_polygons(): st_make_grid() produced zero cells — check inputs/parameters.")
  }
  
  grid_sf <- sf::st_sf(poly_id = seq_along(grid_sfc), geometry = grid_sfc, crs = sf::st_crs(boundary))
  
  # ---- clip (respect holes) --------------------------------------------------
  if (isTRUE(clip)) {
    grid_sf <- suppressWarnings(sf::st_intersection(.make_valid(grid_sf), boundary))
    empt <- sf::st_is_empty(grid_sf)
    if (any(empt, na.rm = TRUE)) grid_sf <- grid_sf[!empt, , drop = FALSE]
    grid_sf$poly_id <- seq_len(nrow(grid_sf))  # re-pack ids
    .msg(sprintf("create_grid_polygons(): clipped grid cells = %d", nrow(grid_sf)))
  }
  
  grid_sf
}

# -----------------------------------------------------------------------------
# Seeding Strategies for Voronoi
# -----------------------------------------------------------------------------

#' K-means seed generation from point coordinates (for Voronoi or grid seeding)
#'
#' Runs \code{\link[stats]{kmeans}} on the coordinates of an \code{sf} point
#' layer and returns the resulting cluster centroids as an \code{sf} POINT
#' layer. The number of centers is clamped to the interval
#' \eqn{[1, n - 1]} (where \eqn{n} is the number of input points) to avoid
#' degenerate cases; at least one center is always returned.
#'
#' @details
#' This function treats the input coordinates as a numeric matrix and applies
#' Euclidean \code{k}-means. If your data are in geographic (lon/lat) degrees,
#' the clustering is performed in degree units; for distance-aware results,
#' project your data to an appropriate planar CRS before calling this function
#' (e.g., via a helper like \code{ensure_projected()} in this package or
#' \code{\link[sf]{st_transform}}).
#'
#' The output CRS matches the CRS of \code{points_sf}. The centroids are the
#' \emph{cluster centers} reported by \code{\link[stats]{kmeans}} (which are not
#' necessarily original observations).
#'
#' @param points_sf An \code{sf} object with POINT or MULTIPOINT geometries to
#'   be clustered.
#' @param k Integer; requested number of clusters. Internally clamped to
#'   \eqn{\max(1, \min(k, n-1))}, where \eqn{n} is the number of input points.
#' @param set_seed Optional integer passed to \code{set.seed()} to make the
#'   \code{kmeans} initialization reproducible. Default is \code{456}.
#'
#' @return An \code{sf} object of \eqn{k^\ast} POINT geometries (where
#'   \eqn{k^\ast = \max(1, \min(k, n-1))}) in the same CRS as \code{points_sf}.
#'   No additional attributes are included.
#'
#' @examples
#' \donttest{
#' library(sf)
#' set.seed(1)
#' pts <- data.frame(x = runif(200, 0, 1000), y = runif(200, 0, 800))
#' pts <- st_as_sf(pts, coords = c("x","y"), crs = 3857)
#'
#' # Request 10 clusters (will clamp if necessary)
#' seeds <- voronoi_seeds_kmeans(pts, k = 10, set_seed = 42)
#'
#' plot(st_geometry(pts), pch = 16, cex = 0.4, main = "K-means seeds")
#' plot(st_geometry(seeds), pch = 3, cex = 1.2, add = TRUE)
#' }
#'
#' @seealso
#' \code{\link{get_voronoi_seeds}},
#' \code{\link{create_voronoi_polygons}},
#' \code{\link{build_tessellation}},
#' \code{\link[sf]{st_transform}},
#' \code{\link[stats]{kmeans}}
#'
#' @importFrom sf st_coordinates st_crs st_as_sf
#' @importFrom stats kmeans
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

#' Random seed generation within a polygonal boundary (Voronoi/grid seeding)
#'
#' Draws \eqn{k} points uniformly at random from inside a polygonal
#' \code{boundary} using \code{\link[sf]{st_sample}} and returns them as an
#' \code{sf} POINT layer. Holes in the boundary are respected (no points are
#' drawn in holes), and multipart boundaries are handled via
#' \code{\link[sf]{st_union}}.
#'
#' @details
#' Sampling is performed in the coordinate space of \code{boundary}. If
#' \code{boundary} is in geographic coordinates (longitude/latitude), the draw
#' is uniform in degrees, which is not area-uniform. For area-uniform sampling,
#' project \code{boundary} to a suitable planar (preferably equal-area) CRS
#' before calling this function (e.g., via \code{\link[sf]{st_transform}} or a
#' helper like \code{ensure_projected()}).
#'
#' The function uses \code{st_sample(..., exact = TRUE)} to guarantee exactly
#' \eqn{k} points (requires an \pkg{sf} version that supports the \code{exact}
#' argument). The returned object’s CRS is set to \code{sf::st_crs(boundary)}.
#'
#' @param boundary An \code{sf} or \code{sfc} polygonal object (POLYGON or
#'   MULTIPOLYGON) defining the region to sample within.
#' @param k Integer; the number of random seed points to generate.
#' @param set_seed Integer passed to \code{\link[base]{set.seed}} for
#'   reproducibility. Default is \code{456}.
#'
#' @return An \code{sf} object of \eqn{k} POINT geometries in the same CRS as
#'   \code{boundary}. No additional attributes are included.
#'
#' @examples
#' \donttest{
#' library(sf)
#'
#' # Simple 1 km square boundary in a projected CRS
#' bb  <- st_bbox(c(xmin = 0, ymin = 0, xmax = 1000, ymax = 1000), crs = 3857)
#' bnd <- st_as_sfc(bb)
#'
#' seeds <- voronoi_seeds_random(bnd, k = 25, set_seed = 123)
#'
#' plot(st_geometry(bnd), border = 'grey40', main = "Random seeds")
#' plot(st_geometry(seeds), pch = 16, cex = 0.8, add = TRUE)
#' }
#'
#' @seealso
#' \code{\link{get_voronoi_seeds}},
#' \code{\link{voronoi_seeds_kmeans}},
#' \code{\link{create_voronoi_polygons}},
#' \code{\link{build_tessellation}},
#' \code{\link[sf]{st_sample}}, \code{\link[sf]{st_transform}}
#'
#' @importFrom sf st_sample st_union st_sf st_set_crs st_crs
#' @export
voronoi_seeds_random <- function(boundary, k, set_seed = 456) {
  set.seed(set_seed)
  pts <- sf::st_sample(sf::st_union(boundary), size = k, type = "random", exact = TRUE)
  sf::st_sf(geometry = pts) |> sf::st_set_crs(sf::st_crs(boundary))
}

# -----------------------------------------------------------------------------
# Assignment: map features to tessellation polygons
# -----------------------------------------------------------------------------

#' Assign features to polygons and attach a polygon ID
#'
#' Joins an \code{sf} layer of input features (points, lines, or polygons) to a
#' polygon layer and returns the features with a polygon identifier column
#' attached. Coordinate reference systems (CRS) are harmonized safely using
#' \code{\link{harmonize_crs}} before the spatial join. If multiple polygons
#' match a single feature, at most one match is retained (via
#' \code{largest = TRUE}) so each feature receives a single polygon ID.
#'
#' @details
#' The polygon identifier column is controlled by \code{polygon_id_col}. If that
#' column does not exist on \code{polygons_sf}, a new integer ID is created in
#' the returned data (1..nrow(polygons)). If the polygons already contain a
#' likely identifier, one of the following candidates is used automatically:
#' \code{polygon_id_col}, \code{"poly_id"}, \code{"polygon_id"}, \code{"id"},
#' \code{"cell_id"}, \code{"grid_id"}.
#'
#' The spatial relationship used for the join is configurable through
#' \code{predicate}. By default, \code{\link[sf]{st_within}} assigns a feature
#' to the polygon that contains it. Alternative predicates such as
#' \code{\link[sf]{st_intersects}} or \code{\link[sf]{st_contains}} may be
#' provided (e.g., \code{predicate = sf::st_intersects}).
#'
#' When a feature matches multiple polygons (e.g., overlapping polygons),
#' \code{\link[sf]{st_join}} is called with \code{largest = TRUE} to keep at most
#' one polygon per feature. If this tie-breaking is not appropriate for your
#' data, resolve overlaps beforehand (e.g., by dissolving polygons) or perform a
#' custom join.
#'
#' The result is returned in the original CRS of \code{features_sf} when that CRS
#' is available.
#'
#' @param features_sf An \code{sf} object containing the features to assign
#'   (geometry may be POINT/LINESTRING/POLYGON, etc.).
#' @param polygons_sf An \code{sf} (or \code{sfc}) polygonal layer providing the
#'   target polygons for assignment.
#' @param polygon_id_col Character string naming the polygon identifier column to
#'   add to the features. If this column does not exist on
#'   \code{polygons_sf}, a new sequential ID is created. Default: \code{"poly_id"}.
#' @param keep_unassigned Logical; if \code{FALSE} (default), features that do
#'   not match any polygon are dropped. If \code{TRUE}, unmatched features are
#'   retained with \code{NA} in \code{polygon_id_col}.
#' @param predicate A binary spatial predicate function from \pkg{sf} used by
#'   \code{\link[sf]{st_join}} (e.g., \code{sf::st_within},
#'   \code{sf::st_intersects}, \code{sf::st_contains}). Default:
#'   \code{sf::st_within}.
#'
#' @return An \code{sf} object with the same rows (subject to
#'   \code{keep_unassigned}) and geometry as \code{features_sf}, plus a single
#'   polygon identifier column named \code{polygon_id_col}.
#'
#' @examples
#' \donttest{
#' library(sf)
#'
#' # Two simple polygons
#' polys <- st_as_sf(data.frame(
#'   id = c(10L, 20L),
#'   wkt = c("POLYGON((0 0, 2 0, 2 2, 0 2, 0 0))",
#'           "POLYGON((2 0, 4 0, 4 2, 2 2, 2 0))")
#' ), wkt = "wkt", crs = 3857)
#'
#' # Some points
#' pts <- st_as_sf(data.frame(
#'   x = c(0.5, 1.5, 2.5, 3.5, 10),
#'   y = c(0.5, 1.5, 1.0, 1.0, 10)
#' ), coords = c("x","y"), crs = 3857)
#'
#' out1 <- assign_features_to_polygons(
#'   features_sf   = pts,
#'   polygons_sf   = polys,
#'   polygon_id_col = "id",
#'   keep_unassigned = TRUE
#' )
#'
#' # Using a different predicate (intersects)
#' out2 <- assign_features_to_polygons(
#'   pts, polys, polygon_id_col = "id", predicate = sf::st_intersects
#' )
#' }
#'
#' @seealso
#' \code{\link{harmonize_crs}},
#' \code{\link[sf]{st_join}},
#' \code{\link[sf]{st_within}},
#' \code{\link[sf]{st_intersects}}
#'
#' @importFrom sf st_join st_transform st_as_sf st_crs
#' @export
assign_features_to_polygons <- function(
    features_sf,
    polygons_sf,
    polygon_id_col = "poly_id",
    keep_unassigned = FALSE,
    predicate = sf::st_within
) {
  if (!inherits(features_sf, "sf"))
    stop("assign_features_to_polygons(): `features_sf` must be an sf object.")
  if (!inherits(polygons_sf, "sf")) {
    if (inherits(polygons_sf, "sfc")) polygons_sf <- sf::st_as_sf(polygons_sf) else
      stop("assign_features_to_polygons(): `polygons_sf` must be sf/sfc.")
  }
  if (nrow(polygons_sf) == 0L)
    stop("assign_features_to_polygons(): `polygons_sf` has zero rows.")
  
  # Harmonize CRS without assuming geometry column name
  orig_crs <- sf::st_crs(features_sf)
  hh <- harmonize_crs(features_sf, polygons_sf)
  f <- hh$a
  p <- hh$b
  
  # Ensure there is an ID column on polygons; create one if needed
  id_candidates <- c(polygon_id_col, "poly_id", "polygon_id", "id", "cell_id", "grid_id")
  id_col <- id_candidates[id_candidates %in% names(p)][1]
  if (length(id_col) == 0L || is.na(id_col)) {
    id_col <- polygon_id_col
    p[[id_col]] <- seq_len(nrow(p))
  }
  
  # Keep only the ID (geometry stays attached automatically)
  p_sel <- p[, id_col, drop = FALSE]
  
  # Spatial join: one polygon per feature; largest handles rare multi-matches
  joined <- sf::st_join(f, p_sel, join = predicate, left = TRUE, largest = TRUE)
  
  # Normalize ID column name to `polygon_id_col`
  if (!identical(id_col, polygon_id_col)) {
    names(joined)[names(joined) == id_col] <- polygon_id_col
  }
  
  # Filter out unassigned unless requested otherwise
  if (!keep_unassigned) {
    joined <- joined[!is.na(joined[[polygon_id_col]]), , drop = FALSE]
  }
  
  # Return in the original CRS (if it existed)
  if (!is.na(orig_crs)) {
    joined <- sf::st_transform(joined, orig_crs)
  }
  joined
}

# -----------------------------------------------------------------------------
# Level Selection (cluster counts / target cells)
# -----------------------------------------------------------------------------

#' Select an elbow (knee) from a WSS curve
#'
#' Heuristically selects the "elbow" (a.k.a. knee) from a vector of
#' within-cluster sum of squares (WSS) values as a function of the number of
#' clusters \eqn{k}. The algorithm looks for the index \eqn{k} in the range
#' \code{min_k:max_k} that maximizes the discrete second difference of the WSS
#' curve. For very short sequences (\eqn{< 3} points), it falls back to the
#' mid-point of the inspected range. Optionally, neighboring \eqn{k} values are
#' returned as candidate alternatives.
#'
#' @details
#' Let \eqn{w_k} denote WSS at \eqn{k}. Over the restricted index set
#' \code{k_idx <- min_k:max_k}, the first and second discrete differences are
#' computed as:
#' \deqn{d1_i = w_{k_{i+1}} - w_{k_i}, \quad d2_i = d1_{i+1} - d1_i.}
#' The elbow \eqn{k^\*} is chosen as the \eqn{k} corresponding to the maximum
#' \eqn{d2} (with an offset so that \eqn{k^\*} is aligned to the original
#' \eqn{k}-indexing of \code{wss}). Any non-finite entries in \code{d2} are
#' treated as \code{-Inf}. If all \code{d2} are non-finite or the inspected
#' sequence has fewer than three points, the elbow defaults to the mid-point
#' \code{floor((min_k + max_k) / 2)}.
#'
#' If \code{return_neighbors = TRUE}, the function also returns the set
#' \code{c(k* - 1, k*, k* + 1)} clipped to \code{[min_k, max_k]} and with
#' duplicates removed.
#'
#' The vector \code{wss} is assumed to be indexed such that \code{wss[k]}
#' corresponds to WSS at \eqn{k} clusters (typically monotonically decreasing).
#'
#' @param wss Numeric vector of within-cluster sum of squares indexed by
#'   \eqn{k}, i.e., \code{wss[k]} corresponds to the WSS at \eqn{k} clusters.
#'   Must have length \eqn{\ge 2}.
#' @param max_k Integer upper bound on \eqn{k} to inspect (inclusive). Defaults
#'   to \code{length(wss)}; values above \code{length(wss)} are truncated.
#' @param min_k Integer lower bound on \eqn{k} to inspect (inclusive). Defaults
#'   to \code{1}; values less than 1 are raised to 1.
#' @param return_neighbors Logical; if \code{TRUE} (default), also return the
#'   neighboring \eqn{k} values \code{k* - 1} and \code{k* + 1} (clipped to the
#'   inspected interval) as \code{candidates}. If \code{FALSE}, only \code{k*}
#'   is returned.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{knee_k} Integer; the selected elbow \eqn{k^\*}.
#'   \item \code{candidates} Integer vector of candidate \eqn{k} values. If
#'         \code{return_neighbors = TRUE}, this is the clipped set
#'         \code{c(k* - 1, k*, k* + 1)} with duplicates removed; otherwise it is
#'         just \code{knee_k}.
#'   \item \code{diagnostics} A list with:
#'     \itemize{
#'       \item \code{wss} The \code{wss} slice actually inspected
#'             (\code{wss[min_k:max_k]}).
#'       \item \code{d1} First differences of that slice (\code{diff(wss)}).
#'       \item \code{d2} Second differences (\code{diff(d1)}), possibly empty if
#'             fewer than three points were inspected.
#'     }
#' }
#'
#' @references
#' Thorndike, R. L. (1953). Who belongs in the family? \emph{Psychometrika},
#' 18(4), 267–276. (The original "elbow" method.)
#'
#' @examples
#' # A simple monotone-decreasing WSS curve for k = 1..10
#' wss <- c(100, 70, 55, 47, 42, 39, 37, 36, 35, 34)
#'
#' # Default search over the full range
#' res <- .elbow_from_wss(wss)
#' res$knee_k
#' res$candidates
#'
#' # Restrict the search range
#' .elbow_from_wss(wss, min_k = 2, max_k = 8)
#'
#' # If the range is very short, falls back to the mid-point
#' .elbow_from_wss(wss, min_k = 4, max_k = 5)
#'
#' @seealso \code{\link[stats]{kmeans}} for producing WSS curves to which this
#'   heuristic is often applied.
#'
#' @keywords internal
#' @noRd
.elbow_from_wss <- function(wss,
                            max_k = length(wss),
                            min_k = 1L,
                            return_neighbors = TRUE) {
  if (!is.numeric(wss) || length(wss) < 2L)
    stop(".elbow_from_wss(): `wss` must be numeric length >= 2.")
  min_k <- as.integer(min_k); max_k <- as.integer(max_k)
  min_k <- max(1L, min_k); max_k <- min(length(wss), max_k)
  if (min_k >= max_k) stop(".elbow_from_wss(): need at least two k values.")
  
  k_idx <- seq.int(min_k, max_k)
  wss_k <- as.numeric(wss[k_idx])
  
  if (length(wss_k) < 3L) {
    knee_k <- floor((min_k + max_k) / 2)
    candidates <- if (return_neighbors) {
      sort(unique(pmin(max_k, pmax(min_k, c(knee_k - 1L, knee_k, knee_k + 1L)))))
    } else knee_k
    return(list(
      knee_k = knee_k,
      candidates = candidates,
      diagnostics = list(wss = wss_k, d1 = diff(wss_k), d2 = numeric(0))
    ))
  }
  
  d1 <- diff(wss_k)
  d2 <- diff(d1)
  d2[!is.finite(d2)] <- -Inf
  
  knee_k <- if (all(!is.finite(d2))) {
    floor((min_k + max_k) / 2)
  } else {
    (min_k - 1L) + (which.max(d2) + 1L)
  }
  
  candidates <- if (return_neighbors) {
    sort(unique(pmin(max_k, pmax(min_k, c(knee_k - 1L, knee_k, knee_k + 1L)))))
  } else knee_k
  
  list(knee_k = knee_k, candidates = candidates, diagnostics = list(wss = wss_k, d1 = d1, d2 = d2))
}

#' Determine an optimal number of spatial levels (clusters) via an elbow heuristic
#'
#' Computes a simple within-cluster sum-of-squares (WSS) curve over \eqn{k=1,\dots,K_{\max}}
#' using k-means on projected feature coordinates and selects candidate \eqn{k}
#' values around the elbow (knee) of that curve. This is useful for choosing
#' the number of tessellation cells or cluster “levels” for spatial models.
#'
#' @details
#' The input \code{data_sf} can contain any geometry type. Non-point geometries
#' are first converted to representative points with
#' \code{\link{coerce_to_points}(mode = "auto")}. If the data are in a geographic
#' CRS (lon/lat), they are projected to a local planar CRS via
#' \code{\link{ensure_projected}} so that Euclidean distances are meaningful for
#' k-means and WSS calculations.
#'
#' For scalability, if there are more than \code{sample_n} features, a random
#' subset of size \code{sample_n} is used to compute the WSS curve (controlled
#' by \code{set_seed} for reproducibility). The WSS for \eqn{k=1} is computed
#' analytically from the overall centroid; for \eqn{k \ge 2} it is taken from
#' \code{\link[stats]{kmeans}} with modest defaults (\code{iter.max = 50},
#' \code{nstart = 5}). If k-means fails at some \eqn{k}, the previous WSS value
#' is reused as a conservative fallback.
#'
#' The elbow is detected by \code{\link{.elbow_from_wss}}, which chooses the
#' \eqn{k} that maximizes the discrete second difference of the WSS curve (with
#' a short-series fallback to the mid-point). Up to \code{top_n} neighboring
#' candidates around the elbow are returned.
#'
#' @param data_sf An \pkg{sf} object containing the input features. Any geometry
#'   type is accepted; non-point geometries are converted to points internally.
#' @param max_levels Integer. Upper bound on the number of levels (clusters)
#'   to consider. The actual maximum evaluated is
#'   \code{min(max_levels, N - 1)}, where \code{N} is the number of points after
#'   any sub-sampling (and at least 2).
#' @param top_n Integer. How many candidate \eqn{k} values to return around the
#'   detected elbow. Defaults to \code{3}.
#' @param sample_n Integer. If \code{nrow(data_sf) > sample_n}, randomly sample
#'   \code{sample_n} features to build the WSS curve for speed. Set larger for
#'   higher fidelity at the cost of compute.
#' @param set_seed Integer seed used for both sub-sampling (when applicable)
#'   and k-means initialization to make results reproducible.
#'
#' @return An integer vector of candidate numbers of levels (clusters), in
#'   ascending order, clipped to \code{[1, k_max]} and with duplicates removed.
#'   If there are fewer than 3 points available, the function returns \code{1L}.
#'
#' @section Algorithm (summary):
#' \enumerate{
#'   \item Convert geometries to points (\code{coerce_to_points(mode = "auto")}).
#'   \item Project to a local planar CRS if needed (\code{ensure_projected()}).
#'   \item Optionally sub-sample to \code{sample_n} points.
#'   \item Compute WSS for \eqn{k = 1, \dots, k_{\max}} using k-means.
#'   \item Detect the elbow with \code{.elbow_from_wss()} and return up to
#'         \code{top_n} neighboring \eqn{k} values.
#' }
#'
#' @note The returned candidates are heuristics; domain knowledge and downstream
#'   validation should guide the final choice of \eqn{k}.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' set.seed(1)
#' pts <- data.frame(
#'   x = rnorm(200, sd = 3),
#'   y = rnorm(200, sd = 1.5)
#' )
#' sf_pts <- st_as_sf(pts, coords = c("x","y"), crs = 4326)
#'
#' # Suggest cluster counts (levels) around the elbow
#' determine_optimal_levels(sf_pts, max_levels = 10, top_n = 3, sample_n = 150, set_seed = 42)
#' }
#'
#' @seealso \code{\link[stats]{kmeans}}, \code{\link{.elbow_from_wss}},
#'   \code{\link{coerce_to_points}}, \code{\link{ensure_projected}}
#'
#' @export
determine_optimal_levels <- function(data_sf,
                                     max_levels = 12L,
                                     top_n = 3L,
                                     sample_n = 1500L,
                                     set_seed = 123L) {
  if (!inherits(data_sf, "sf")) stop("determine_optimal_levels(): `data_sf` must be an sf object.")
  
  # points + projected coords
  if (!all(sf::st_geometry_type(data_sf, by_geometry = TRUE) %in% c("POINT","MULTIPOINT"))) {
    data_sf <- coerce_to_points(data_sf, "auto")
  }
  data_sf <- ensure_projected(data_sf)
  
  xy <- sf::st_coordinates(data_sf)[, 1:2, drop = FALSE]
  n  <- nrow(xy)
  if (n < 3L) return(as.integer(1L))
  
  # optional sub-sample for speed
  set.seed(set_seed)
  if (n > sample_n) {
    xy <- xy[sample(seq_len(n), sample_n), , drop = FALSE]
  }
  
  k_max <- max(2L, min(as.integer(max_levels), nrow(xy) - 1L))
  wss <- numeric(k_max)
  
  for (k in seq_len(k_max)) {
    if (k == 1L) {
      ctr <- colMeans(xy)
      wss[k] <- sum(rowSums((xy - matrix(ctr, nrow(xy), 2, byrow = TRUE))^2))
    } else {
      km <- try(stats::kmeans(xy, centers = k, iter.max = 50, nstart = 5), silent = TRUE)
      if (inherits(km, "try-error")) {
        # fallback: reuse previous WSS or compute via simple cluster assignment to random centers
        wss[k] <- wss[k - 1L]
      } else {
        wss[k] <- km$tot.withinss
      }
    }
  }
  
  elbow <- .elbow_from_wss(wss, max_k = k_max, min_k = 1L, return_neighbors = TRUE)
  out <- as.integer(head(elbow$candidates, max(1L, as.integer(top_n))))
  out[out < 1L]   <- 1L
  out[out > k_max] <- k_max
  unique(out)
}

# -----------------------------------------------------------------------------
# Modeling: GWR and Bayesian Spatial (spBayes)
# -----------------------------------------------------------------------------

#' Validate presence of model variables in a data frame / sf object
#'
#' Checks that the specified response and predictor variables exist as columns
#' of \code{data_sf}. Intended as a lightweight guard prior to fitting models
#' (e.g., GWR or Bayesian spatial models).
#'
#' @param data_sf A \code{data.frame} or \pkg{sf} object containing the data.
#'   Only the column names are inspected; geometry (if present) is ignored.
#' @param response_var A single character string naming the response variable
#'   column that must be present in \code{data_sf}.
#' @param predictor_vars A character vector of predictor (feature) column names
#'   that must be present in \code{data_sf}. May be \code{NULL} or length 0 if
#'   no predictors are required.
#'
#' @return Invisibly returns \code{TRUE} if all requested columns are present.
#'   Otherwise, the function throws an error listing the missing variables.
#'
#' @section Errors:
#' If any of \code{response_var} or \code{predictor_vars} are not found in
#' \code{names(data_sf)}, an error is raised of the form
#' \dQuote{Missing variables in data: var1, var2, ...}.
#'
#' @examples
#' df <- data.frame(y = 1:5, x1 = rnorm(5), x2 = runif(5))
#' .validate_model_inputs(df, response_var = "y", predictor_vars = c("x1","x2"))
#' # returns (invisibly) TRUE
#'
#' \dontrun{
#' .validate_model_inputs(df, response_var = "y", predictor_vars = c("x1","x3"))
#' # Error: Missing variables in data: x3
#' }
#'
#' @keywords internal
#' @noRd
.validate_model_inputs <- function(data_sf, response_var, predictor_vars) {
  miss <- setdiff(c(response_var, predictor_vars), names(data_sf))
  if (length(miss)) stop(sprintf("Missing variables in data: %s", paste(miss, collapse = ", ")))
  invisible(TRUE)
}

#' Prepare and sanitize an \pkg{sf} dataset for spatial modeling
#'
#' Ensures the input features are point geometries in a projected CRS and
#' removes rows with missing or non-finite values in the modeling columns.
#' Optionally aligns the CRS to a provided \code{boundary}. This is a common
#' pre-processing step before fitting spatial models (e.g., GWR or Bayesian
#' spatial models) where planar distances/areas are expected.
#'
#' @param data_sf An \pkg{sf} object containing the response/predictor columns
#'   and geometries. If geometries are not points, they are converted to points
#'   according to \code{pointize}.
#' @param response_var A length-1 character string naming the response variable
#'   column in \code{data_sf}.
#' @param predictor_vars A non-empty character vector naming predictor (feature)
#'   columns in \code{data_sf}.
#' @param boundary Optional \pkg{sf} or \pkg{sfc} object used only to align CRS.
#'   When supplied, \code{data_sf} is projected to the CRS of \code{boundary}
#'   (via \code{\link{ensure_projected}}). The geometry type of
#'   \code{boundary} is not enforced; it is not used for clipping here.
#' @param pointize Strategy for converting non-point geometries to points,
#'   one of:
#'   \itemize{
#'     \item \code{"auto"} (default) — sensible type-specific choice:
#'       centroids or midpoints for lines; \code{st_point_on_surface()} for
#'       polygons; first point for multipoints; pass-through for points.
#'     \item \code{"surface"} — alias of \code{"point_on_surface"}.
#'     \item \code{"centroid"} — \code{sf::st_centroid()}.
#'     \item \code{"point_on_surface"} — \code{sf::st_point_on_surface()} (always
#'       inside polygon where possible).
#'     \item \code{"line_midpoint"} — midpoint along each \strong{LINESTRING}
#'       (projects temporarily when needed).
#'     \item \code{"bbox_center"} — numeric center of each feature's bounding
#'       box (no guarantee of lying on/in the geometry).
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates presence of \code{response_var} and \code{predictor_vars}
#'         in \code{data_sf}.
#'   \item Coerces geometries to points using \code{\link{coerce_to_points}}
#'         if needed, per \code{pointize}.
#'   \item Projects to a suitable planar CRS using \code{\link{ensure_projected}}:
#'         to \code{boundary}'s CRS if provided; otherwise chooses a local
#'         projected CRS when input is in lon/lat.
#'   \item Drops rows with missing values across the modeling columns and, for
#'         numeric modeling columns only, drops rows with non-finite values
#'         (e.g., \code{NaN}, \code{Inf}). A warning is emitted indicating how
#'         many rows were removed; if the \pkg{logger} package is available,
#'         \code{logger::log_warn()} is used, otherwise a base warning.
#' }
#'
#' @return An \pkg{sf} object (points) in a projected CRS (aligned to
#'   \code{boundary} when supplied), containing only rows with complete and
#'   finite modeling values. Column names and all non-geometry attributes are
#'   preserved for retained rows.
#'
#' @section Errors:
#' \itemize{
#'   \item If \code{data_sf} is not an \pkg{sf} object.
#'   \item If \code{response_var} is not a single character string.
#'   \item If \code{predictor_vars} is not a non-empty character vector.
#'   \item If any required columns are missing from \code{data_sf}.
#' }
#'
#' @seealso
#' \code{\link{coerce_to_points}} for geometry pointization,
#' \code{\link{ensure_projected}} for robust CRS projection/alignment.
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' # Example data: polygons with attributes
#' polys <- st_as_sf(
#'   data.frame(id = 1:3, y = rnorm(3), x1 = runif(3), x2 = rnorm(3)),
#'   wkt = c("POLYGON((0 0,1 0,1 1,0 1,0 0))",
#'           "POLYGON((1 0,2 0,2 1,1 1,1 0))",
#'           "POLYGON((0 1,1 1,1 2,0 2,0 1))"),
#'   crs = 4326
#' )
#'
#' # Optional boundary used only for CRS alignment
#' bnd <- st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 2, ymax = 2), crs = st_crs(4326)))
#'
#' # Prepare for modeling: polygon -> point, projected CRS, drop incomplete rows
#' data_pts <- prep_model_data(
#'   data_sf       = polys,
#'   response_var  = "y",
#'   predictor_vars= c("x1","x2"),
#'   boundary      = bnd,
#'   pointize      = "point_on_surface"
#' )
#' }
#'
#' @export
prep_model_data <- function(data_sf,
                            response_var,
                            predictor_vars,
                            boundary = NULL,
                            pointize = c("auto","surface","centroid","line_midpoint","bbox_center")) {
  # ---- Validate inputs -------------------------------------------------------
  if (!inherits(data_sf, "sf")) stop("prep_model_data(): 'data_sf' must be an sf object.")
  pointize <- match.arg(pointize)
  
  if (!is.character(response_var) || length(response_var) != 1L)
    stop("prep_model_data(): 'response_var' must be a single column name.")
  
  if (!is.character(predictor_vars) || length(predictor_vars) < 1L)
    stop("prep_model_data(): 'predictor_vars' must be a non-empty character vector.")
  
  req_cols <- c(response_var, predictor_vars)
  miss <- setdiff(req_cols, names(data_sf))
  if (length(miss)) {
    stop("prep_model_data(): missing required column(s): ", paste(miss, collapse = ", "))
  }
  
  # ---- Geometry coercion (points) -------------------------------------------
  if (!all(sf::st_geometry_type(data_sf, by_geometry = TRUE) %in% c("POINT","MULTIPOINT"))) {
    data_sf <- coerce_to_points(data_sf, pointize)
  }
  
  # ---- CRS alignment/projection ---------------------------------------------
  if (!is.null(boundary)) {
    bnd <- if (inherits(boundary, "sfc")) sf::st_as_sf(boundary) else boundary
    if (!inherits(bnd, "sf")) stop("prep_model_data(): 'boundary' must be sf/sfc when supplied.")
    bnd <- ensure_projected(bnd)
    data_sf <- ensure_projected(data_sf, bnd)
  } else {
    data_sf <- ensure_projected(data_sf)
  }
  
  # ---- Row filtering: complete cases & finite numeric values ----------------
  df <- sf::st_drop_geometry(data_sf)[, req_cols, drop = FALSE]
  
  # complete cases across all requested columns
  ok_cc <- stats::complete.cases(df)
  
  # finite check only on numeric columns (avoid errors on factors/characters)
  num_mask <- vapply(df, is.numeric, logical(1))
  ok_fin <- if (any(num_mask)) {
    apply(as.matrix(df[, num_mask, drop = FALSE]), 1L, function(r) all(is.finite(r)))
  } else {
    rep(TRUE, nrow(df))
  }
  
  keep <- ok_cc & ok_fin
  dropped <- sum(!keep)
  if (dropped > 0) {
    if (requireNamespace("logger", quietly = TRUE)) {
      logger::log_warn("prep_model_data(): dropping %d row(s) with non-finite or missing values in modeling columns.", dropped)
    } else {
      warning(sprintf("prep_model_data(): dropping %d row(s) with non-finite or missing values in modeling columns.", dropped))
    }
  }
  
  data_sf[keep, , drop = FALSE]
}

# --- small internal helper ----------------------------------------------------

#' Compute basic regression error metrics (RMSE, MAE, R²)
#'
#' Calculates common accuracy metrics comparing observed values (\code{y})
#' to predictions (\code{yhat}). Non-finite pairs (where either \code{y} or
#' \code{yhat} is not finite) are removed prior to calculation.
#'
#' @param y Numeric vector of observed/true values. Should be the same length
#'   as \code{yhat} and correspond element-wise.
#' @param yhat Numeric vector of predicted/estimated values. Should be the same
#'   length as \code{y} and correspond element-wise.
#'
#' @details
#' After removing non-finite pairs, the following are computed:
#' \itemize{
#'   \item \strong{n}: number of finite pairs used.
#'   \item \strong{RMSE}: root mean squared error,
#'     \deqn{\mathrm{RMSE} = \sqrt{\frac{1}{n}\sum_i (y_i - \hat{y}_i)^2}}
#'   \item \strong{MAE}: mean absolute error,
#'     \deqn{\mathrm{MAE} = \frac{1}{n}\sum_i |y_i - \hat{y}_i|}
#'   \item \strong{R2}: unadjusted coefficient of determination,
#'     \deqn{R^2 = 1 - \frac{\sum_i (y_i - \hat{y}_i)^2}{\sum_i (y_i - \bar{y})^2}}
#'     computed relative to the mean of \code{y}. If the total sum of squares
#'     is zero (e.g., \code{y} is constant), \code{R2} is returned as \code{NA}.
#' }
#'
#' @return A base \code{data.frame} with columns:
#' \itemize{
#'   \item \code{n} (integer): number of finite pairs.
#'   \item \code{RMSE} (numeric): root mean squared error (same units as \code{y}).
#'   \item \code{MAE} (numeric): mean absolute error (same units as \code{y}).
#'   \item \code{R2} (numeric): coefficient of determination in \eqn{[0, 1]} when defined,
#'         otherwise \code{NA_real_}.
#' }
#'
#' @note If all pairs are filtered out (i.e., \code{n == 0}), \code{RMSE} and
#'   \code{MAE} will be \code{NaN} and \code{R2} will be \code{NA}.
#'
#' @examples
#' y    <- c(1.0, 2.0, 3.0, 4.0, 5.0)
#' yhat <- c(0.9, 2.1, 2.8, 4.2, 5.1)
#' .compute_reg_metrics(y, yhat)
#'
#' # With non-finite values (filtered pairwise)
#' y2    <- c(1, 2, NA, 4, 5)
#' yhat2 <- c(1, NaN, 3, 4, Inf)
#' .compute_reg_metrics(y2, yhat2)
#'
#' @keywords internal
#' @noRd
.compute_reg_metrics <- function(y, yhat) {
  ok <- is.finite(y) & is.finite(yhat)
  y    <- y[ok]
  yhat <- yhat[ok]
  n <- length(y)
  rss <- sum((y - yhat)^2)
  tss <- sum( (y - mean(y))^2 )
  rmse <- sqrt(rss / n)
  mae  <- mean(abs(y - yhat))
  r2   <- if (tss > 0) 1 - rss / tss else NA_real_
  data.frame(n = n, RMSE = rmse, MAE = mae, R2 = r2)
}

#' Fit a Geographically Weighted Regression (GWR) with adaptive bandwidth
#'
#' Fits an adaptive-kernel GWR using \pkg{spgwr} on an \code{sf} dataset.
#' Non-point geometries are converted to representative points (via
#' \code{\link{coerce_to_points}}), rows with missing or non-finite modeling
#' values are dropped (with a warning), and geographic (lon/lat) data are
#' projected to a local metric CRS (via \code{\link{ensure_projected}}) to
#' ensure stable distance calculations.
#'
#' @param data_sf An \code{sf} object containing the response, predictors,
#'   and geometries. If geometries are not \code{POINT}/\code{MULTIPOINT},
#'   they are coerced to points automatically.
#' @param response_var A length-1 character string naming the response column
#'   in \code{data_sf}.
#' @param predictor_vars A non-empty character vector naming the predictor
#'   columns in \code{data_sf}.
#'
#' @details
#' \strong{Preprocessing}
#' \itemize{
#'   \item If \code{data_sf} geometries are not \code{POINT}/\code{MULTIPOINT},
#'         they are converted to points using \code{coerce_to_points("auto")}.
#'   \item Rows are filtered to those that are complete across \code{response_var}
#'         and \code{predictor_vars}; non-finite numeric values are also removed.
#'         The number of dropped rows is reported via \pkg{logger} if available
#'         (otherwise a warning).
#'   \item If the data are in longitude/latitude, they are projected to a local,
#'         metric CRS using \code{ensure_projected()} before distance-based
#'         computations. If \code{ensure_projected()} is not available, an error
#'         is thrown when data are long/lat.
#' }
#'
#' \strong{Bandwidth selection and fitting}
#' \itemize{
#'   \item Bandwidth is selected adaptively via \code{spgwr::gwr.sel(..., adapt = TRUE)}.
#'         If selection fails or returns a non-positive/NA value, a fallback value
#'         of \code{0.75} is used. If the selected value is \eqn{\ge 1}, it is
#'         clamped to \code{0.99}.
#'   \item The model is then fit with \code{spgwr::gwr(..., adapt = bw, hatmatrix = TRUE)}.
#' }
#'
#' \strong{AICc extraction}
#' The function extracts AICc from common \pkg{spgwr} output locations
#' (\code{$GW.diagnostic$AICc}, \code{$results$AICc}, or a \code{$results} data.frame
#' column). If not found, \code{AICc} is returned as \code{NA_real_}.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{model}}{The fitted \code{spgwr::gwr} object.}
#'   \item{\code{bandwidth}}{Numeric scalar, the adaptive bandwidth used (in \eqn{(0,1)}).}
#'   \item{\code{adaptive}}{Logical flag (always \code{TRUE}).}
#'   \item{\code{AICc}}{Numeric AICc value if available, otherwise \code{NA_real_}.}
#' }
#'
#' @section Dependencies:
#' Requires \pkg{spgwr} and \pkg{sp}. The function also relies on \pkg{sf}
#' and optionally \pkg{logger} for warnings. Helper functions
#' \code{\link{coerce_to_points}} and \code{\link{ensure_projected}} are used if present.
#'
#' @seealso \code{\link{prep_model_data}}, \code{\link{coerce_to_points}},
#'   \code{\link{ensure_projected}}, \code{spgwr::gwr}, \code{spgwr::gwr.sel}
#'
#' @examples
#' \dontrun{
#' library(sf)
#' # Example data: random points with a simple linear relation
#' set.seed(1)
#' pts <- st_as_sf(
#'   data.frame(x = runif(200, -1, 1), y = runif(200, -1, 1)),
#'   coords = c("x","y"), crs = 4326
#' )
#' pts$resp <- rnorm(200)
#' pts$pred1 <- rnorm(200)
#' pts$pred2 <- rnorm(200)
#'
#' fit <- fit_gwr_model(
#'   data_sf       = pts,
#'   response_var  = "resp",
#'   predictor_vars = c("pred1","pred2")
#' )
#' fit$bandwidth
#' fit$AICc
#' }
#'
#' @export
fit_gwr_model <- function(data_sf, response_var, predictor_vars) {
  if (!inherits(data_sf, "sf")) {
    stop("fit_gwr_model(): `data_sf` must be an sf object.")
  }
  if (!requireNamespace("spgwr", quietly = TRUE)) {
    stop("fit_gwr_model(): package 'spgwr' is required.")
  }
  if (!requireNamespace("sp", quietly = TRUE)) {
    stop("fit_gwr_model(): package 'sp' is required.")
  }
  
  # Ensure point geometry
  gtypes <- as.character(sf::st_geometry_type(data_sf, by_geometry = TRUE))
  if (!all(gtypes %in% c("POINT", "MULTIPOINT"))) {
    data_sf <- coerce_to_points(data_sf, "auto")
  }
  
  # Keep only rows that are complete AND finite on modeling columns
  keep_cols <- unique(c(response_var, predictor_vars))
  miss <- setdiff(keep_cols, names(data_sf))
  if (length(miss)) {
    stop(sprintf("fit_gwr_model(): missing columns: %s", paste(miss, collapse = ", ")))
  }
  
  vals   <- sf::st_drop_geometry(data_sf)[, keep_cols, drop = FALSE]
  cc_na  <- stats::complete.cases(vals)
  numcol <- vapply(vals, is.numeric, logical(1))
  finite_ok <- if (any(numcol)) {
    apply(as.matrix(vals[, numcol, drop = FALSE]), 1L, function(x) all(is.finite(x)))
  } else {
    rep(TRUE, nrow(vals))
  }
  cc <- cc_na & finite_ok
  
  if (!all(cc)) {
    n_drop <- sum(!cc)
    if (requireNamespace("logger", quietly = TRUE)) {
      logger::log_warn("fit_gwr_model(): dropping %d rows with non-finite/NA in modeling columns.", n_drop)
    } else {
      warning(sprintf("fit_gwr_model(): dropping %d rows with non-finite/NA in modeling columns.", n_drop))
    }
  }
  dat <- data_sf[cc, , drop = FALSE]
  
  # ---- CRS guard: project lon/lat to a metric CRS for stable distances ----
  is_ll <- tryCatch({
    if (isTRUE(get0("st_is_longlat", envir = asNamespace("sf"), inherits = FALSE))) {
      sf::st_is_longlat(dat)
    } else if (exists(".is_longlat", mode = "function")) {
      .is_longlat(dat)
    } else {
      FALSE
    }
  }, error = function(e) FALSE)
  
  if (isTRUE(is_ll)) {
    if (!exists("ensure_projected", mode = "function")) {
      stop("fit_gwr_model(): data are long/lat but `ensure_projected()` is not available.")
    }
    dat <- ensure_projected(dat)
  }
  
  # Build formula (handles non-syntactic names safely)
  fml <- stats::reformulate(termlabels = predictor_vars, response = response_var)
  
  # Convert to sp for spgwr
  sp_dat <- methods::as(dat, "Spatial")
  
  # Helper: fixed-distance fallback (retained for future use)
  .fixed_bw_fallback <- function(spobj) {
    coords <- sp::coordinates(spobj)
    n <- nrow(coords)
    s <- if (n > 1000) sample.int(n, 1000) else seq_len(n)
    d <- as.matrix(stats::dist(coords[s, , drop = FALSE]))
    stats::median(d[upper.tri(d)], na.rm = TRUE)
  }
  
  # ---- Adaptive bandwidth selection ----
  bw <- try(suppressWarnings(spgwr::gwr.sel(fml, data = sp_dat, adapt = TRUE)), silent = TRUE)
  if (inherits(bw, "try-error") || !is.finite(bw) || is.na(bw) || bw <= 0) {
    if (requireNamespace("logger", quietly = TRUE)) {
      logger::log_warn("fit_gwr_model(): adaptive bandwidth selection failed; using 0.75.")
    } else {
      warning("fit_gwr_model(): adaptive bandwidth selection failed; using 0.75.")
    }
    bw <- 0.75
  } else if (bw >= 1) {
    if (requireNamespace("logger", quietly = TRUE)) {
      logger::log_warn("fit_gwr_model(): adaptive bandwidth >= 1; clamping to 0.99.")
    } else {
      warning("fit_gwr_model(): adaptive bandwidth >= 1; clamping to 0.99.")
    }
    bw <- 0.99
  }
  
  # ---- Fit GWR (adaptive) ----
  fit <- try(spgwr::gwr(fml, data = sp_dat, adapt = bw, hatmatrix = TRUE), silent = TRUE)
  if (inherits(fit, "try-error")) {
    stop(sprintf("fit_gwr_model(): GWR fit failed: %s", as.character(fit)))
  }
  
  # ---- AICc extraction (defensive across spgwr variants) ----
  AICc_val <- NA_real_
  if (!is.null(fit$GW.diagnostic) && !is.null(fit$GW.diagnostic$AICc)) {
    AICc_val <- suppressWarnings(as.numeric(fit$GW.diagnostic$AICc))
  } else if (!is.null(fit$results) && is.list(fit$results) && !is.null(fit$results$AICc)) {
    AICc_val <- suppressWarnings(as.numeric(fit$results$AICc))
  } else if (!is.null(fit$results) && is.data.frame(fit$results) && "AICc" %in% names(fit$results)) {
    AICc_val <- suppressWarnings(as.numeric(fit$results$AICc[1]))
  }
  
  list(
    model     = fit,
    bandwidth = as.numeric(bw),
    adaptive  = TRUE,
    AICc      = AICc_val
  )
}

#' Heuristic bounds for the spatial decay parameter \eqn{\phi}
#'
#' Computes rule-of-thumb lower/upper bounds for the spatial decay (range)
#' parameter \eqn{\phi} used in common covariance models (e.g., exponential,
#' spherical, Matérn as parameterized in \pkg{spBayes}, where correlation
#' typically decays like \eqn{\exp(-\phi d)}). Bounds are derived from the
#' empirical pairwise distances among sampling locations.
#'
#' @param coords_xy A numeric matrix/data frame of coordinates with two
#'   columns \code{(x, y)} (projected/metric units are strongly recommended).
#'   Distances are computed via \code{stats::dist()}; zero distances (duplicate
#'   locations) are removed before computing quantiles.
#' @param q_small Numeric scalar in \eqn{(0, 1)} giving the lower distance
#'   quantile used to define the upper bound (default \code{0.1} for the
#'   10th percentile of the nonzero pairwise distances).
#'
#' @details
#' Let \eqn{d_{\max}} denote the maximum nonzero pairwise distance and
#' \eqn{d_q} the \code{q_small} quantile of nonzero pairwise distances.
#' Using the “effective range” heuristic \eqn{d_{\text{eff}} \approx 3/\phi}
#' for exponential-type models:
#' \itemize{
#'   \item The lower bound is set to \eqn{\phi_L = 3 / d_{\max}}, allowing
#'         very long spatial dependence across the entire domain.
#'   \item The upper bound is set to
#'         \eqn{\phi_U = \max(1.2\,\phi_L,\ 3 / d_q)}, allowing short-range
#'         dependence at about the \code{q_small} inter-site distance and
#'         ensuring \eqn{\phi_U > \phi_L}.
#' }
#'
#' If all pairwise distances are zero (e.g., a single unique location),
#' the function returns the conservative fallback \code{c(0.001, 1)}.
#'
#' @return A numeric vector of length 2: \code{c(lower, upper)} for \eqn{\phi}.
#'   These are suitable, for example, as \code{phi.Unif = c(lower, upper)}
#'   when specifying priors in \pkg{spBayes}.
#'
#' @section Units and CRS:
#' Distances are computed in the units of \code{coords_xy}. For meaningful
#' bounds, supply projected (metric) coordinates (e.g., UTM). If your data are
#' in longitude/latitude degrees, project first (e.g., with \pkg{sf}:
#' \code{ensure_projected()} in this package).
#'
#' @examples
#' set.seed(1)
#' xy <- cbind(runif(200, 0, 10000), runif(200, 0, 8000))  # meters
#' .phi_prior_bounds(xy)
#'
#' # Tighter upper bound by using a larger lower-distance quantile:
#' .phi_prior_bounds(xy, q_small = 0.25)
#'
#' @seealso \code{\link[spBayes]{spLM}}, \code{\link[spBayes]{spMvLM}}
#'
#' @export
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

#' Matérn correlation function \eqn{\rho(h \,|\, \phi, \nu)}
#'
#' Evaluates the Matérn correlation for a set of pairwise distances \code{h}
#' using decay/range parameter \code{phi} and smoothness \code{nu}:
#' \deqn{\rho(h) \;=\; \frac{2^{\,1-\nu}}{\Gamma(\nu)} \, (\phi h)^{\nu} \, K_{\nu}(\phi h),}
#' where \eqn{K_{\nu}(\cdot)} is the modified Bessel function of the second kind and
#' \eqn{\rho(0)=1}.
#'
#' @param h A numeric object coercible to a matrix via \code{as.matrix(h)}.
#'   Commonly a numeric vector of distances or a symmetric pairwise-distance
#'   matrix (e.g., \code{as.matrix(dist(coords))}). Units must match the
#'   inverse units of \code{phi}.
#' @param phi Positive numeric scalar; decay/range parameter (units: \eqn{1/\text{distance}}).
#'   Larger \code{phi} yields faster correlation decay. For the exponential case
#'   (\eqn{\nu=0.5}), a common “effective range” heuristic is \eqn{\approx 3/\phi}.
#' @param nu Positive numeric scalar; smoothness parameter. Typical choices are
#'   \code{0.5} (exponential), \code{1.5}, \code{2.5}; as \code{nu -> Inf} the kernel
#'   approaches Gaussian.
#'
#' @details
#' The implementation:
#' \itemize{
#'   \item Coerces \code{h} to a matrix and computes \eqn{z=\phi h}.
#'   \item Initializes a zero matrix, then sets the diagonal to 1 explicitly to
#'         enforce \eqn{\rho(0)=1}.
#'   \item For strictly positive distances (\eqn{z>0}), evaluates the closed-form
#'         Matérn correlation using \code{\link[base]{besselK}} and clips tiny
#'         negative numerical artifacts to 0.
#' }
#'
#' \strong{Note on duplicate locations:} Off-diagonal entries with \code{h == 0}
#' (exact duplicates) are \emph{not} set to 1 by this function and will remain 0.
#' In typical workflows duplicates are removed upstream; if you need off-diagonal
#' zeros replaced by 1 for duplicate pairs, handle that before/after calling this
#' helper.
#'
#' The function assumes \code{phi > 0} and \code{nu > 0}; supplying non-positive
#' values will yield undefined or \code{NaN} results.
#'
#' @return A numeric matrix the same dimension as \code{as.matrix(h)} containing
#'   Matérn correlations. The diagonal entries are 1.
#'
#' @references
#' Matérn, B. (1960). \emph{Spatial Variation}. Springer.
#'
#' Rasmussen, C. E., & Williams, C. K. I. (2006). \emph{Gaussian Processes for
#' Machine Learning}. MIT Press.
#'
#' @seealso \code{\link[base]{besselK}}, \code{\link[base]{gamma}}
#'
#' @examples
#' # Pairwise distances for random points in a square (projected units)
#' set.seed(42)
#' coords <- cbind(runif(20, 0, 10000), runif(20, 0, 10000))
#' H <- as.matrix(dist(coords))
#'
#' # Matérn correlation with nu = 1.5 (once-differentiable) and phi = 1/3000
#' R <- .matern_rho(H, phi = 1/3000, nu = 1.5)
#' range(R)
#'
#' # Correlation as a function of distance for multiple nu values
#' d <- seq(0, 10000, length.out = 100)
#' r_exp  <- .matern_rho(d, phi = 1/3000, nu = 0.5)  # exponential
#' r_m32  <- .matern_rho(d, phi = 1/3000, nu = 1.5)
#' r_m52  <- .matern_rho(d, phi = 1/3000, nu = 2.5)
#' # r_exp, r_m32, r_m52 are column vectors (matrices with 100 x 1)
#'
#' @keywords internal
#' @noRd
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

#' Safe multivariate normal log-likelihood with jitter fallback
#'
#' Computes the log-density of a multivariate normal \eqn{y \sim \mathcal{N}(\mu, V)}
#' using \code{mvtnorm::dmvnorm}. If the evaluation fails (e.g., because \code{V}
#' is near-singular or not numerically positive definite), a small diagonal
#' jitter is added and the evaluation is retried.
#'
#' @param y Numeric vector of observations.
#' @param mu Numeric vector of means, same length as \code{y}.
#' @param V Numeric covariance matrix of dimension \eqn{length(y) \times length(y)}.
#'   Should be symmetric and (numerically) positive semi-definite.
#'
#' @details
#' The function first attempts \code{mvtnorm::dmvnorm(y, mean = mu, sigma = V, log = TRUE)}.
#' On error, it perturbs the covariance by adding
#' \deqn{\epsilon I,\quad \epsilon = 10^{-8}\,\big(\mathrm{mean}(\mathrm{diag}(V)) + \epsilon_{\mathrm{mach}}\big),}
#' and retries once. The perturbation is applied to a local copy of \code{V}.
#' If the second attempt fails, the underlying error is propagated.
#'
#' @return A single numeric value: the log-density \eqn{\log p(y \mid \mu, V)}.
#'
#' @section Dependencies:
#' Requires the \pkg{mvtnorm} package.
#'
#' @examples
#' set.seed(1)
#' mu <- c(0, 0)
#' V  <- matrix(c(1, 0.9, 0.9, 1), 2, 2)
#' y  <- c(0.2, -0.1)
#' ll <- .safe_ll(y, mu, V)
#'
#' # Near-singular covariance (will trigger jitter fallback)
#' V_sing <- matrix(c(1, 1, 1, 1), 2, 2)
#' ll2 <- .safe_ll(y, mu, V_sing)
#'
#' @keywords internal
#' @noRd
.safe_ll <- function(y, mu, V) {
  out <- try(mvtnorm::dmvnorm(y, mean = mu, sigma = V, log = TRUE), silent = TRUE)
  if (inherits(out, "try-error")) {
    eps <- 1e-8 * (mean(diag(V)) + .Machine$double.eps)
    diag(V) <- diag(V) + eps
    out <- mvtnorm::dmvnorm(y, mean = mu, sigma = V, log = TRUE)
  }
  out
}

#' Fit a Bayesian spatial regression with a 2D Gaussian Process (via {brms})
#'
#' Fits a Bayesian regression model with a spatial Gaussian Process (GP) term
#' using \pkg{brms}. Geometry is coerced to points (if needed), data are projected
#' to a metric CRS for stable distances, and a GP smoother
#' \code{gp(..x, ..y, k = gp_k, c = gp_c)} is added to the right-hand side of the
#' model formula. Optionally computes PSIS-LOO diagnostics.
#'
#' @param data_sf An \pkg{sf} object containing the response, predictors, and
#'   geometries. Non-point geometries are coerced to points according to
#'   \code{pointize}.
#' @param response_var Single string; the name of the response column in \code{data_sf}.
#' @param predictor_vars Character vector of predictor column names in \code{data_sf}.
#'   Use \code{character(0)} to fit an intercept-only model plus the GP term.
#' @param family A \pkg{brms} family object (e.g., \code{gaussian()}, \code{poisson()},
#'   \code{bernoulli()}, ...). Default is \code{gaussian()}.
#' @param gp_k Positive integer passed to \code{brms::gp()} as \code{k}. Controls the
#'   rank/complexity of the GP approximation (larger \code{k} = more flexible, slower).
#' @param gp_c Positive numeric passed to \code{brms::gp()} as \code{c}. See
#'   \code{?brms::gp} for details; typically influences the GP scale.
#' @param prior Optional \pkg{brms} prior specification created with
#'   \code{brms::prior()}, \code{brms::set_prior()}, or a list thereof. Default \code{NULL}
#'   lets \pkg{brms} use its defaults.
#' @param chains Number of MCMC chains. Passed to \code{brms::brm()}.
#' @param iter Total iterations per chain. Passed to \code{brms::brm()}.
#' @param warmup Warmup iterations per chain (default \code{floor(iter/2)}). Passed to
#'   \code{brms::brm()}.
#' @param cores Number of parallel cores for sampling. Default uses all but one.
#' @param seed Integer seed for reproducibility. Passed to \code{brms::brm()}.
#' @param backend Sampling backend: \code{"auto"} (default; prefers \pkg{cmdstanr} if
#'   available, otherwise \pkg{rstan}), or explicitly \code{"cmdstanr"} or \code{"rstan"}.
#' @param control Named list of sampler controls passed to \code{brms::brm()} (e.g.,
#'   \code{list(adapt_delta = 0.9, max_treedepth = 12)}).
#' @param compute_loo Logical; if \code{TRUE} (default), compute PSIS-LOO via
#'   \code{brms::loo()}. When available, LOOIC is returned; otherwise it is derived as
#'   \code{-2 * elpd\_loo}.
#' @param pointize One of \code{"auto"}, \code{"surface"}, \code{"centroid"},
#'   \code{"line_midpoint"}, or \code{"bbox_center"}. Controls how non-point geometries
#'   are converted to points using \code{coerce_to_points()} before modeling.
#'   Default \code{"auto"}.
#' @param boundary Optional polygonal \pkg{sf}/\pkg{sfc} object used only to harmonize
#'   CRS/projection with \code{data_sf}. It is not used as a mask in the model.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates inputs and selects a sampling backend.
#'   \item Calls \code{prep_model_data()} to (a) coerce geometries to points according to
#'         \code{pointize}, (b) project data to a metric CRS (using
#'         \code{ensure_projected()}), and (c) drop rows with missing/non-finite modeling
#'         values.
#'   \item Extracts projected coordinates into auxiliary columns \code{..x}, \code{..y}.
#'   \item Builds a formula of the form
#'         \code{response ~ predictors + gp(..x, ..y, k = gp_k, c = gp_c)} and fits with
#'         \code{brms::brm()}.
#'   \item Optionally computes PSIS-LOO via \code{brms::loo()} and returns LOO objects and
#'         LOOIC (if available).
#' }
#'
#' Models are fit in projected (non-long/lat) coordinates to ensure that spatial
#' distances are in consistent metric units. If \code{data_sf} is long/lat, it will be
#' projected via \code{ensure_projected()} (typically to a local UTM).
#'
#' @return A named \code{list} with components:
#' \describe{
#'   \item{\code{model}}{A \code{brmsfit} object.}
#'   \item{\code{loo}}{A PSIS-LOO object from \code{brms::loo()}, or \code{NULL} if not computed.}
#'   \item{\code{looic}}{Numeric LOOIC if available, otherwise \code{NA_real_}.}
#'   \item{\code{coords}}{Character vector \code{c("..x","..y")} indicating the coordinate columns used.}
#'   \item{\code{n}}{Integer number of observations used in the fit.}
#' }
#'
#' @section Dependencies:
#' Requires \pkg{brms}. Depending on \code{backend}, also requires \pkg{cmdstanr}
#' (with a working CmdStan installation) or \pkg{rstan} with a functioning C++ toolchain.
#' Relies on helper utilities in this package: \code{prep_model_data()},
#' \code{coerce_to_points()}, and \code{ensure_projected()}.
#'
#' @note
#' This function does not mask/clip data to \code{boundary}; \code{boundary} is used only
#' for CRS harmonization. For reproducibility, set \code{seed}. Large \code{gp_k} increases
#' fidelity but can substantially slow down sampling.
#'
#' @seealso \code{\link{prep_model_data}}, \code{\link{coerce_to_points}},
#'   \code{\link{ensure_projected}}, \code{\link[brms]{gp}}, \code{\link[brms]{brm}},
#'   \code{\link[brms]{loo}}
#'
#' @examples
#' \dontrun{
#' # Synthetic example (requires brms + working backend)
#' library(sf)
#' set.seed(1)
#' pts  <- st_as_sf(data.frame(x = runif(200), y = runif(200)),
#'                  coords = c("x","y"), crs = 4326)
#' pts  <- st_transform(pts, 3857)  # project to meters
#' XY   <- st_coordinates(pts)
#' eta  <- sin(XY[,1] / 5e4) + cos(XY[,2] / 5e4)
#' y    <- eta + rnorm(nrow(pts), sd = 0.3)
#' dat  <- cbind(pts, y = y, z1 = rnorm(nrow(pts)))
#'
#' fit <- fit_bayesian_spatial_model(
#'   data_sf        = dat,
#'   response_var   = "y",
#'   predictor_vars = "z1",
#'   gp_k           = 5,
#'   gp_c           = 1.5,
#'   backend        = "cmdstanr",
#'   chains         = 2, iter = 1000, warmup = 500, cores = 2
#' )
#' print(fit$model)
#' if (!is.null(fit$loo)) print(fit$loo)
#' }
#'
#' @export
fit_bayesian_spatial_model <- function(
    data_sf,
    response_var,
    predictor_vars,
    family      = gaussian(),
    gp_k        = 5,
    gp_c        = 1.5,
    prior       = NULL,
    chains      = 4,
    iter        = 2000,
    warmup      = floor(iter / 2),
    cores       = max(1L, parallel::detectCores() - 1L),
    seed        = 123,
    backend     = c("auto", "cmdstanr", "rstan"),
    control     = list(adapt_delta = 0.9, max_treedepth = 12),
    compute_loo = TRUE,
    pointize    = "auto",
    boundary    = NULL
) {
  if (!inherits(data_sf, "sf")) {
    stop("fit_bayesian_spatial_model(): `data_sf` must be an sf object.")
  }
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("fit_bayesian_spatial_model(): package 'brms' is required.")
  }
  
  backend <- match.arg(backend)
  if (identical(backend, "auto")) {
    backend <- if (requireNamespace("cmdstanr", quietly = TRUE)) "cmdstanr" else "rstan"
  }
  
  # ---- Centralized prep ------------------------------------------------------
  if (!exists("prep_model_data", mode = "function")) {
    stop("fit_bayesian_spatial_model(): `prep_model_data()` not found.")
  }
  
  dat_sf <- prep_model_data(
    data_sf        = data_sf,
    response_var   = response_var,
    predictor_vars = predictor_vars,
    boundary       = boundary,
    pointize       = pointize
  )
  
  # Must be POINT/MULTIPOINT and projected
  gtypes <- as.character(sf::st_geometry_type(dat_sf, by_geometry = TRUE))
  if (!all(gtypes %in% c("POINT", "MULTIPOINT"))) {
    stop("fit_bayesian_spatial_model(): geometry must be POINT/MULTIPOINT after prep_model_data().")
  }
  is_ll <- tryCatch(sf::st_is_longlat(dat_sf), error = function(e) FALSE)
  if (isTRUE(is_ll)) {
    if (!exists("ensure_projected", mode = "function")) {
      stop("fit_bayesian_spatial_model(): data are long/lat but `ensure_projected()` is not available.")
    }
    dat_sf <- ensure_projected(dat_sf)
  }
  
  # ---- Build data frame with coords -----------------------------------------
  coords <- sf::st_coordinates(dat_sf)
  if (!all(c("X", "Y") %in% colnames(coords))) {
    colnames(coords)[1:2] <- c("X", "Y")
  }
  dat_df <- cbind(
    sf::st_drop_geometry(dat_sf),
    `..x` = coords[, "X"],
    `..y` = coords[, "Y"]
  )
  
  # ---- Build formula: predictors (+ intercept) + gp(..x, ..y) ---------------
  if (length(predictor_vars) == 0L) {
    rhs_terms <- "1"
  } else {
    base_fml  <- stats::reformulate(termlabels = predictor_vars)
    rhs_terms <- as.character(base_fml)[2L]  # safe with non-syntactic names
  }
  gp_term <- sprintf("gp(..x, ..y, k = %s, c = %s)", gp_k, gp_c)
  fml_txt <- sprintf("%s ~ %s + %s", response_var, rhs_terms, gp_term)
  fml     <- stats::as.formula(fml_txt)
  
  # ---- Fit with brms (use do.call to set backend cleanly) --------------------
  brm_args <- list(
    formula = fml,
    data    = dat_df,
    family  = family,
    prior   = prior,
    chains  = chains,
    iter    = iter,
    warmup  = warmup,
    cores   = cores,
    seed    = seed,
    control = control
  )
  if (backend %in% c("cmdstanr", "rstan")) {
    brm_args$backend <- backend
  }
  
  fit <- try(do.call(brms::brm, brm_args), silent = TRUE)
  if (inherits(fit, "try-error")) {
    stop(sprintf("fit_bayesian_spatial_model(): brms fit failed: %s", as.character(fit)))
  }
  
  # ---- PSIS-LOO (optional) ---------------------------------------------------
  loo_obj <- NULL
  looic   <- NA_real_
  if (isTRUE(compute_loo)) {
    loo_try <- try(brms::loo(fit), silent = TRUE)
    if (!inherits(loo_try, "try-error")) {
      loo_obj <- loo_try
      est <- try(loo_obj$estimates, silent = TRUE)
      if (!inherits(est, "try-error") && is.matrix(est)) {
        rn <- rownames(est)
        if ("looic" %in% rn) {
          looic <- suppressWarnings(as.numeric(est["looic", "Estimate"]))
        } else if ("elpd_loo" %in% rn) {
          # fallback: LOOIC = -2 * elpd_loo
          looic <- suppressWarnings(as.numeric(-2 * est["elpd_loo", "Estimate"]))
        }
      }
    } else {
      if (requireNamespace("logger", quietly = TRUE)) {
        logger::log_warn("fit_bayesian_spatial_model(): LOO computation failed.")
      } else {
        warning("fit_bayesian_spatial_model(): LOO computation failed.", call. = FALSE)
      }
    }
  }
  
  list(
    model  = fit,
    loo    = if (!is.null(loo_obj)) loo_obj else NULL,
    looic  = looic,
    coords = c("..x", "..y"),
    n      = nrow(dat_df)
  )
}

# -----------------------------------------------------------------------------
# Modeling Pipeline Wrapper
# -----------------------------------------------------------------------------

#' Evaluate spatial models (GWR and/or Bayesian GP) with optional cross-validation
#'
#' Compares two spatial modeling strategies—Geographically Weighted Regression
#' (GWR, via \pkg{spgwr}) and a Bayesian regression with a 2D Gaussian Process
#' (via \pkg{brms})—either using \emph{cross-validation} or \emph{in-sample}
#' fits. Inputs are harmonized (geometry coerced to points if needed and projected
#' to a metric CRS) and model performance metrics are reported.
#'
#' @param data_sf An \pkg{sf} object containing the response, predictors, and
#'   geometries. Non-point geometries are allowed and will be converted to points
#'   according to \code{pointize}.
#' @param response_var Single string; name of the response column in \code{data_sf}.
#' @param predictor_vars Character vector of predictor column names in
#'   \code{data_sf}. May be \code{character(0)} to fit an intercept-only model.
#' @param do_cv Logical; if \code{TRUE} (default) evaluates models using
#'   cross-validation (via \code{cv_gwr()} and/or \code{cv_bayes()} if available).
#'   If \code{FALSE}, fits models in-sample and computes metrics on the training
#'   data.
#' @param folds Optional fold assignment object accepted by your
#'   \code{cv_gwr()} / \code{cv_bayes()} implementations. If \code{NULL}, those
#'   functions are expected to generate \code{k}-fold splits using \code{k} and
#'   \code{seed}.
#' @param k Integer; number of folds when \code{folds} is \code{NULL}. Default \code{5}.
#' @param seed Integer seed for reproducible CV splits (and any internal sampling).
#' @param boundary Optional polygonal \pkg{sf}/\pkg{sfc} object used solely to
#'   harmonize CRS with \code{data_sf} (via \code{prep_model_data()}); it is not
#'   used to mask training data.
#' @param pointize One of \code{"auto"}, \code{"surface"}, \code{"centroid"},
#'   \code{"line_midpoint"}, or \code{"bbox_center"}. Controls how non-point
#'   geometries are converted to points before modeling; passed to
#'   \code{coerce_to_points()} inside \code{prep_model_data()}.
#' @param gwr_args Named list of extra arguments for \code{fit_gwr_model()}.
#'   Used only in the in-sample path (\code{do_cv = FALSE}); not forwarded to
#'   \code{cv_gwr()} unless your implementation accepts them explicitly.
#' @param bayes_args Named list of extra arguments for
#'   \code{fit_bayesian_spatial_model()}. In the CV path, these are passed as
#'   \code{fit_args} to \code{cv_bayes()} \emph{only if} that function has a
#'   \code{fit_args} formal; otherwise they are ignored. In the in-sample path,
#'   they are forwarded directly to \code{fit_bayesian_spatial_model()}.
#' @param summary Character; when summarizing Bayesian posterior predictions in
#'   the in-sample path, use \code{"mean"} or \code{"median"} (default \code{"mean"}).
#' @param models Character vector choosing which models to evaluate. Allowed
#'   values are \code{"GWR"} and \code{"Bayesian"}. If a requested model's
#'   dependencies/functions are unavailable, it is silently dropped with a
#'   message; if none remain, the function errors.
#' @param quiet Logical; if \code{TRUE}, suppress progress messages.
#'
#' @details
#' \strong{Cross-validation path (\code{do_cv = TRUE}, default):}
#' \itemize{
#'   \item Requires helper functions \code{cv_gwr()} and/or \code{cv_bayes()} to be
#'         available in the namespace. If \pkg{spgwr} + \code{fit_gwr_model()} +
#'         \code{cv_gwr()} are not available, GWR is dropped; similarly for
#'         \pkg{brms} + \code{fit_bayesian_spatial_model()} + \code{cv_bayes()}.
#'   \item Each CV function is called with \code{data_sf}, \code{response_var},
#'         \code{predictor_vars}, \code{folds}/\code{k}/\code{seed}, and
#'         CRS/pointization options. If \code{cv_bayes()} exposes a \code{fit_args}
#'         argument, \code{bayes_args} is supplied there.
#'   \item The returned objects from the CV helpers are included verbatim, and an
#'         overall comparison table (\code{comparison}) is assembled with columns
#'         \code{model}, \code{RMSE}, \code{MAE}, \code{R2}, and
#'         (if available) \code{n_pred}.
#' }
#'
#' \strong{In-sample path (\code{do_cv = FALSE}):}
#' \itemize{
#'   \item Calls \code{prep_model_data()} to coerce geometries to points (per
#'         \code{pointize}), project to a metric CRS, and drop rows with non-finite
#'         values in modeling columns.
#'   \item Fits GWR via \code{fit_gwr_model()} (if selected) and extracts fitted
#'         values from the \code{$SDF} slot (or falls back to
#'         \code{spgwr::gwr.predict()} if needed).
#'   \item Fits the Bayesian model via \code{fit_bayesian_spatial_model()} (if selected)
#'         and computes fitted responses using posterior expectations
#'         (\code{posterior_epred}) summarized by \code{summary}.
#'   \item For each fitted model, returns RMSE, MAE, and \(R^2\) computed against
#'         the in-sample truth.
#' }
#'
#' @return
#' If \code{do_cv = TRUE}, a list containing some of:
#' \describe{
#'   \item{\code{gwr_cv}}{The object returned by \code{cv_gwr()} (if GWR evaluated).}
#'   \item{\code{bayes_cv}}{The object returned by \code{cv_bayes()} (if Bayesian evaluated).}
#'   \item{\code{comparison}}{A \code{data.frame}/tibble with columns such as
#'         \code{model}, \code{RMSE}, \code{MAE}, \code{R2}, and optionally \code{n_pred}.}
#' }
#'
#' If \code{do_cv = FALSE}, a list with components:
#' \describe{
#'   \item{\code{formula}}{Model formula used (as a string).}
#'   \item{\code{data}}{The processed \pkg{sf} data used for fitting.}
#'   \item{\code{gwr}}{List with GWR fit information (if fitted): \code{model},
#'         \code{bandwidth}, \code{adaptive = TRUE}, and \code{AICc}.}
#'   \item{\code{bayes}}{A \code{brmsfit} (or related) object for the Bayesian model
#'         (if fitted).}
#'   \item{\code{metrics}}{A \code{data.frame} summarizing in-sample RMSE, MAE, and \(R^2\)
#'         for each fitted model.}
#' }
#'
#' @section Dependencies:
#' Optional runtime dependencies include \pkg{spgwr} for GWR and \pkg{brms}
#' (plus a working backend such as \pkg{cmdstanr} or \pkg{rstan}) for the Bayesian
#' model. Cross-validation also requires package/helper functions
#' \code{cv_gwr()} and/or \code{cv_bayes()} to be present. Data preparation relies
#' on internal helpers: \code{prep_model_data()}, \code{coerce_to_points()},
#' and \code{ensure_projected()}.
#'
#' @seealso
#' \code{\link{prep_model_data}}, \code{\link{fit_gwr_model}},
#' \code{\link{fit_bayesian_spatial_model}}, \code{\link{cv_gwr}},
#' \code{\link{cv_bayes}}
#'
#' @examples
#' \dontrun{
#' library(sf)
#' set.seed(42)
#' n <- 250
#' pts <- st_as_sf(
#'   data.frame(x = runif(n), y = runif(n), z1 = rnorm(n), z2 = rnorm(n)),
#'   coords = c("x","y"), crs = 4326
#' )
#' pts <- st_transform(pts, 3857)
#' XY  <- st_coordinates(pts)
#' y   <- 1 + 0.6*pts$z1 - 0.4*pts$z2 + sin(XY[,1]/5e4) + rnorm(n, sd = 0.3)
#' dat <- cbind(pts, y = y)
#'
#' # Cross-validated comparison (requires cv_gwr()/cv_bayes())
#' eval_cv <- evaluate_models(
#'   data_sf        = dat,
#'   response_var   = "y",
#'   predictor_vars = c("z1","z2"),
#'   do_cv          = TRUE,
#'   k              = 5,
#'   models         = c("GWR","Bayesian")
#' )
#' eval_cv$comparison
#'
#' # In-sample fits (no CV)
#' eval_fit <- evaluate_models(
#'   data_sf        = dat,
#'   response_var   = "y",
#'   predictor_vars = c("z1","z2"),
#'   do_cv          = FALSE,
#'   models         = c("GWR","Bayesian"),
#'   bayes_args     = list(chains = 2, iter = 1000, warmup = 500, backend = "cmdstanr")
#' )
#' eval_fit$metrics
#' }
#'
#' @export
evaluate_models <- function(
    data_sf,
    response_var,
    predictor_vars,
    do_cv      = TRUE,
    folds      = NULL,
    k          = 5,
    seed       = 123,
    boundary   = NULL,
    pointize   = "auto",
    gwr_args   = list(),
    bayes_args = list(),
    summary    = c("mean", "median"),
    models     = c("GWR", "Bayesian"),
    quiet      = FALSE
) {
  summary <- match.arg(summary)
  .msg <- function(...) if (!quiet) message(...)
  .has <- function(pkg) requireNamespace(pkg, quietly = TRUE)
  .filter_args <- function(fun, args_list) {
    fml <- try(formals(fun), silent = TRUE)
    if (inherits(fml, "try-error")) return(list())
    ok_names <- names(fml)
    args_list[names(args_list) %in% ok_names]
  }
  
  # ---- guards ----
  if (!inherits(data_sf, "sf")) stop("evaluate_models(): `data_sf` must be an sf object.")
  if (!is.character(response_var) || length(response_var) != 1L)
    stop("evaluate_models(): `response_var` must be a single column name.")
  if (!all(predictor_vars %in% names(sf::st_drop_geometry(data_sf))))
    stop("evaluate_models(): some `predictor_vars` are missing in `data_sf`.")
  
  models <- unique(intersect(models, c("GWR","Bayesian")))
  if (length(models) == 0L) models <- "GWR"
  
  has_spgwr <- .has("spgwr")
  has_brms  <- .has("brms")
  
  has_fit_gwr   <- exists("fit_gwr_model", mode = "function")
  has_fit_bayes <- exists("fit_bayesian_spatial_model", mode = "function")
  has_cv_gwr    <- exists("cv_gwr", mode = "function")
  has_cv_bayes  <- exists("cv_bayes", mode = "function")
  
  if ("GWR" %in% models && !(has_spgwr && has_fit_gwr && has_cv_gwr)) {
    .msg("evaluate_models(): dropping GWR (needs spgwr + fit_gwr_model() + cv_gwr()).")
    models <- setdiff(models, "GWR")
  }
  if ("Bayesian" %in% models && !(has_brms && has_fit_bayes && has_cv_bayes)) {
    .msg("evaluate_models(): dropping Bayesian (needs brms + fit_bayesian_spatial_model() + cv_bayes()).")
    models <- setdiff(models, "Bayesian")
  }
  if (length(models) == 0L) stop("evaluate_models(): no viable models after availability checks.")
  .msg("evaluate_models(): using models = [", paste(models, collapse = ", "), "]")
  
  # Ensure an internal row id (useful for CV impls that expect original indexing)
  if (!("..row_id" %in% names(data_sf))) data_sf$`..row_id` <- seq_len(nrow(data_sf))
  
  # ---------------- CV path ----------------
  if (isTRUE(do_cv)) {
    if (!is.null(seed)) set.seed(seed)
    
    comparison_rows <- list()
    out <- list()
    
    if ("GWR" %in% models) {
      .msg("evaluate_models(): running CV for GWR ...")
      base_args <- list(
        data_sf        = data_sf,
        response_var   = response_var,
        predictor_vars = predictor_vars,
        folds          = folds,
        k              = k,
        seed           = seed,
        boundary       = boundary,
        pointize       = pointize
        # NOTE: do NOT pass fit_args here unless cv_gwr() supports it
      )
      args_gwr <- .filter_args(cv_gwr, base_args)
      gwr_cv <- do.call(cv_gwr, args_gwr)
      out$gwr_cv <- gwr_cv
      
      ov <- try(as.data.frame(gwr_cv$overall), silent = TRUE)
      if (inherits(ov, "try-error") || nrow(ov) == 0L)
        ov <- data.frame(RMSE=NA_real_, MAE=NA_real_, R2=NA_real_)
      if (!"n_pred" %in% names(ov)) ov$n_pred <- length(predictor_vars)
      ov$model <- "GWR"
      comparison_rows[["GWR"]] <- ov[, intersect(c("model","RMSE","MAE","R2","n_pred"), names(ov)), drop = FALSE]
    }
    
    if ("Bayesian" %in% models) {
      .msg("evaluate_models(): running CV for Bayesian ...")
      base_args <- list(
        data_sf        = data_sf,
        response_var   = response_var,
        predictor_vars = predictor_vars,
        folds          = folds,
        k              = k,
        seed           = seed,
        boundary       = boundary,
        pointize       = pointize,
        summary        = summary
        # Only pass fit_args if cv_bayes() supports it
      )
      # try to add fit_args if accepted
      if ("fit_args" %in% names(formals(cv_bayes))) base_args$fit_args <- bayes_args
      args_bayes <- .filter_args(cv_bayes, base_args)
      bayes_cv <- do.call(cv_bayes, args_bayes)
      out$bayes_cv <- bayes_cv
      
      ov <- try(as.data.frame(bayes_cv$overall), silent = TRUE)
      if (inherits(ov, "try-error") || nrow(ov) == 0L)
        ov <- data.frame(RMSE=NA_real_, MAE=NA_real_, R2=NA_real_)
      if (!"n_pred" %in% names(ov)) ov$n_pred <- length(predictor_vars)
      ov$model <- "Bayesian"
      comparison_rows[["Bayesian"]] <- ov[, intersect(c("model","RMSE","MAE","R2","n_pred"), names(ov)), drop = FALSE]
    }
    
    comparison <- dplyr::bind_rows(comparison_rows)
    return(c(out, list(comparison = comparison)))
  }
  
  # ---------------- In-sample path ----------------
  if (!exists("prep_model_data", mode = "function"))
    stop("evaluate_models(): `prep_model_data()` not found (needed for in-sample fits).")
  
  dat_sf <- prep_model_data(
    data_sf        = data_sf,
    response_var   = response_var,
    predictor_vars = predictor_vars,
    boundary       = boundary,
    pointize       = pointize
  )
  
  gtypes <- as.character(sf::st_geometry_type(dat_sf, by_geometry = TRUE))
  if (!all(gtypes %in% c("POINT","MULTIPOINT")))
    stop("evaluate_models(): geometry must be POINT/MULTIPOINT after prep_model_data().")
  if (isTRUE(tryCatch(sf::st_is_longlat(dat_sf), error = function(e) FALSE)))
    stop("evaluate_models(): data are still lon/lat after prep; ensure `prep_model_data()` projects.")
  
  # Add x/y if absent
  if (!all(c("x","y") %in% names(dat_sf))) {
    crd <- sf::st_coordinates(dat_sf)
    dat_sf$x <- crd[,1]; dat_sf$y <- crd[,2]
  }
  
  fml <- if (length(predictor_vars) == 0L) {
    stats::as.formula(paste(response_var, "~ 1"))
  } else {
    stats::as.formula(paste(response_var, "~", paste(predictor_vars, collapse = " + ")))
  }
  
  y_true <- sf::st_drop_geometry(dat_sf)[[response_var]]
  .rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))
  .mae  <- function(y, yhat) mean(abs(y - yhat), na.rm = TRUE)
  .r2   <- function(y, yhat) {
    ym <- mean(y, na.rm = TRUE)
    1 - sum((y - yhat)^2, na.rm = TRUE)/sum((y - ym)^2, na.rm = TRUE)
  }
  
  ret <- list(formula = deparse(fml), data = dat_sf)
  metrics_rows <- list()
  
  # ---- GWR (optional) ----
  if ("GWR" %in% models) {
    .msg("evaluate_models(): fitting in-sample GWR ...")
    gwr_fit <- do.call(
      fit_gwr_model,
      c(list(
        data_sf        = dat_sf,
        response_var   = response_var,
        predictor_vars = predictor_vars
      ), gwr_args)
    )
    ret$gwr <- list(
      model     = gwr_fit$model,
      bandwidth = gwr_fit$bandwidth,
      adaptive  = TRUE,
      AICc      = if (!is.null(gwr_fit$AICc)) gwr_fit$AICc else NA_real_
    )
    
    yhat_gwr <- tryCatch({
      if (!is.null(gwr_fit$model$SDF)) {
        sdf <- gwr_fit$model$SDF
        if ("pred"    %in% names(sdf)) as.numeric(sdf$pred)
        else if ("yhat"   %in% names(sdf)) as.numeric(sdf$yhat)
        else if ("fitted" %in% names(sdf)) as.numeric(sdf$fitted)
        else stop("No fitted column in SDF")
      } else stop("No SDF slot on GWR fit")
    }, error = function(e) rep(NA_real_, length(y_true)))
    
    if (!any(is.finite(yhat_gwr))) {
      .msg("  GWR: SDF predictions not found; trying spgwr::gwr.predict() fallback ...")
      sp_dat <- methods::as(dat_sf, "Spatial")
      bw <- gwr_fit$bandwidth
      pred_try <- try(
        spgwr::gwr.predict(
          formula     = fml,
          data        = sp_dat,
          predictdata = sp_dat,
          adapt       = bw,
          se.fit      = FALSE
        ), silent = TRUE
      )
      if (!inherits(pred_try, "try-error") && !is.null(pred_try$SDF)) {
        sdfp <- pred_try$SDF
        if ("pred" %in% names(sdfp))      yhat_gwr <- as.numeric(sdfp$pred)
        else if ("yhat" %in% names(sdfp)) yhat_gwr <- as.numeric(sdfp$yhat)
      }
    }
    
    met <- if (!any(is.finite(yhat_gwr))) {
      data.frame(model="GWR", RMSE=NA_real_, MAE=NA_real_, R2=NA_real_, n=0L)
    } else {
      data.frame(
        model="GWR",
        RMSE=.rmse(y_true, yhat_gwr),
        MAE =.mae (y_true, yhat_gwr),
        R2  =.r2  (y_true, yhat_gwr),
        n   =sum(is.finite(y_true) & is.finite(yhat_gwr))
      )
    }
    metrics_rows[["GWR"]] <- met
  }
  
  # ---- Bayesian (optional) ----
  if ("Bayesian" %in% models) {
    .msg("evaluate_models(): fitting in-sample Bayesian ...")
    bayes_fit <- do.call(
      fit_bayesian_spatial_model,
      c(list(
        data_sf        = dat_sf,
        response_var   = response_var,
        predictor_vars = predictor_vars
      ), bayes_args)
    )
    bayes_obj <- if (!is.null(bayes_fit$model)) bayes_fit$model else bayes_fit
    ret$bayes <- bayes_obj
    
    newdata_df <- sf::st_drop_geometry(dat_sf)
    crd <- sf::st_coordinates(dat_sf)
    newdata_df$`..x` <- crd[,1]; newdata_df$`..y` <- crd[,2]
    if (!all(c("x","y") %in% names(newdata_df))) { newdata_df$x <- crd[,1]; newdata_df$y <- crd[,2] }
    
    yhat_bayes <- rep(NA_real_, length(y_true)); pred_ok <- FALSE
    if (inherits(bayes_obj, "brmsfit") && has_brms) {
      draws <- try(brms::posterior_epred(bayes_obj, newdata = newdata_df), silent = TRUE)
      if (!inherits(draws, "try-error") && is.matrix(draws)) {
        yhat_bayes <- if (summary == "mean") colMeans(draws) else apply(draws, 2L, stats::median)
        pred_ok <- TRUE
      }
    }
    if (!pred_ok && inherits(bayes_obj, "stanreg") && .has("rstanarm")) {
      draws <- try(rstanarm::posterior_linpred(bayes_obj, newdata = newdata_df, transform = TRUE), silent = TRUE)
      if (!inherits(draws, "try-error") && is.matrix(draws)) {
        yhat_bayes <- if (summary == "mean") colMeans(draws) else apply(draws, 2L, stats::median)
        pred_ok <- TRUE
      }
    }
    if (!pred_ok) {
      p1 <- try(stats::predict(bayes_obj, newdata = newdata_df, type = "response"), silent = TRUE)
      if (!inherits(p1, "try-error") && is.numeric(p1) && length(p1) == length(y_true)) {
        yhat_bayes <- as.numeric(p1); pred_ok <- TRUE
      }
    }
    if (!pred_ok) {
      p2 <- try(stats::predict(bayes_obj, newdata = newdata_df), silent = TRUE)
      if (!inherits(p2, "try-error") && is.numeric(p2) && length(p2) == length(y_true)) {
        yhat_bayes <- as.numeric(p2); pred_ok <- TRUE
      }
    }
    
    met <- if (!pred_ok) {
      data.frame(model="Bayesian", RMSE=NA_real_, MAE=NA_real_, R2=NA_real_, n=0L)
    } else {
      data.frame(
        model="Bayesian",
        RMSE=.rmse(y_true, yhat_bayes),
        MAE =.mae (y_true, yhat_bayes),
        R2  =.r2  (y_true, yhat_bayes),
        n   =sum(is.finite(y_true) & is.finite(yhat_bayes))
      )
    }
    metrics_rows[["Bayesian"]] <- met
  }
  
  ret$metrics <- dplyr::bind_rows(metrics_rows)
  ret
}

# -----------------------------------------------------------------------------
# Optional: Aggregations & Plotting
# -----------------------------------------------------------------------------

#' Summarize features by polygon/cell ID (counts and means)
#'
#' Aggregates an \pkg{sf} point dataset—where each row has a polygon/cell
#' identifier—into one row per cell, returning the per-cell count and optional
#' means for a response and/or selected predictor variables.
#'
#' Geometry is dropped prior to summarization. The cell identifier column is
#' resolved robustly by checking, in order, \code{id_col}, \code{"poly_id"},
#' and \code{"polygon_id"}. If none is found, the function errors.
#'
#' @param assigned_points_sf An \pkg{sf} object of features that have already
#'   been assigned to polygons/cells (e.g., via
#'   \code{\link{assign_features_to_polygons}}). Must contain a cell identifier
#'   column (see \code{id_col}).
#' @param response_var Optional single string naming a numeric response column
#'   for which the per-cell mean (\code{mean_response}) will be computed.
#'   If the column is absent or non-numeric, it is skipped (with a message
#'   when \code{quiet = FALSE}).
#' @param predictor_vars Optional character vector of predictor column names.
#'   Only predictors that are present and numeric are summarized; each retained
#'   predictor contributes a per-cell mean column named \code{mean_<predictor>}.
#'   Missing or non-numeric predictors are skipped (with a message when
#'   \code{quiet = FALSE}).
#' @param id_col String with the preferred name of the polygon/cell ID column.
#'   If this column is not found, the function falls back to \code{"poly_id"}
#'   and then \code{"polygon_id"}.
#' @param quiet Logical; if \code{FALSE}, prints informative messages about
#'   ID resolution and skipped variables. Default \code{TRUE}.
#'
#' @details
#' The function performs:
#' \enumerate{
#'   \item \code{sf::st_drop_geometry()} to remove geometry.
#'   \item Resolution of the ID column among \code{id_col}, \code{"poly_id"},
#'         and \code{"polygon_id"}.
#'   \item A base per-cell count (\code{n}).
#'   \item Optional per-cell mean for \code{response_var} (if present and numeric),
#'         stored in \code{mean_response}.
#'   \item Optional per-cell means for numeric \code{predictor_vars}, each named
#'         \code{mean_<predictor>}.
#' }
#' Joins are performed solely by the resolved ID column name; the ID column
#' name is preserved in the output.
#'
#' @return A tibble/data.frame with one row per polygon/cell ID containing:
#' \itemize{
#'   \item the ID column (name as resolved),
#'   \item \code{n}: the number of assigned features in the cell,
#'   \item \code{mean_response}: optional mean of \code{response_var},
#'   \item \code{mean_<predictor>}: optional means for each numeric predictor
#'         found in \code{predictor_vars}.
#' }
#'
#' @seealso \code{\link{assign_features_to_polygons}}
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(dplyr)
#'
#' # Example points with a polygon/cell id
#' set.seed(1)
#' pts <- st_as_sf(
#'   data.frame(
#'     x = runif(100), y = runif(100),
#'     poly_id = sample(1:5, 100, replace = TRUE),
#'     y_resp = rnorm(100),
#'     x1 = rnorm(100), x2 = rnorm(100)
#'   ),
#'   coords = c("x","y"), crs = 4326
#' )
#'
#' # Summarize counts and means by cell
#' summary_tbl <- summarize_by_cell(
#'   assigned_points_sf = pts,
#'   response_var       = "y_resp",
#'   predictor_vars     = c("x1","x2"),
#'   id_col             = "poly_id",
#'   quiet              = FALSE
#' )
#' summary_tbl
#' }
#'
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr group_by summarise left_join across all_of n
#' @importFrom rlang .data
#' @export
summarize_by_cell <- function(assigned_points_sf,
                              response_var  = NULL,
                              predictor_vars = NULL,
                              id_col        = "poly_id",
                              quiet         = TRUE) {
  .msg <- function(...) if (!quiet) message(...)
  # Drop geometry
  df <- sf::st_drop_geometry(assigned_points_sf)
  
  # Resolve ID column robustly
  id_candidates <- unique(c(id_col, "poly_id", "polygon_id"))
  id_found <- id_candidates[id_candidates %in% names(df)]
  if (length(id_found) == 0L) {
    stop(
      "summarize_by_cell(): could not find an ID column. ",
      "Looked for: ", paste(id_candidates, collapse = ", "), "."
    )
  }
  id_col <- id_found[[1]]
  .msg(sprintf("summarize_by_cell(): using id_col = '%s'", id_col))
  
  # Base count per cell
  out <- df |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  
  # Optional response mean
  if (!is.null(response_var)) {
    if (response_var %in% names(df)) {
      if (!is.numeric(df[[response_var]])) {
        .msg(sprintf("summarize_by_cell(): response_var '%s' is not numeric; skipping mean.", response_var))
      } else {
        resp <- df |>
          dplyr::group_by(.data[[id_col]]) |>
          dplyr::summarise(mean_response = mean(.data[[response_var]], na.rm = TRUE),
                           .groups = "drop")
        # Correct join: just join by the id column name
        out <- dplyr::left_join(out, resp, by = id_col)
      }
    } else {
      .msg(sprintf("summarize_by_cell(): response_var '%s' not found; skipping.", response_var))
    }
  }
  
  # Optional predictor means (only keep numeric predictors that are present)
  if (!is.null(predictor_vars)) {
    keep <- predictor_vars[predictor_vars %in% names(df)]
    if (length(keep)) {
      is_num <- vapply(keep, function(nm) is.numeric(df[[nm]]), logical(1))
      keep_num <- keep[is_num]
      if (!all(is_num)) {
        .msg(sprintf(
          "summarize_by_cell(): non-numeric predictors skipped: %s",
          paste(keep[!is_num], collapse = ", ")
        ))
      }
      if (length(keep_num)) {
        preds <- df |>
          dplyr::group_by(.data[[id_col]]) |>
          dplyr::summarise(
            dplyr::across(
              dplyr::all_of(keep_num),
              ~ mean(.x, na.rm = TRUE),
              .names = "mean_{.col}"
            ),
            .groups = "drop"
          )
        out <- dplyr::left_join(out, preds, by = id_col)
      }
    } else {
      .msg("summarize_by_cell(): none of the requested predictor_vars are present; skipping.")
    }
  }
  
  # Ensure the ID column name is preserved (dplyr already does this)
  names(out)[names(out) == id_col] <- id_col
  out
}

#' Plot a tessellation map with optional boundary, seeds, and features
#'
#' Builds a layered \pkg{ggplot2} map of polygon tessellations and optional
#' overlays for a study boundary, seed points, and additional features.
#' Optionally maps a fill aesthetic from a tessellation attribute and draws
#' per-cell text labels.
#'
#' @param tessellation_sf An \pkg{sf} object of \code{POLYGON}/\code{MULTIPOLYGON}
#'   geometries representing the tessellation to plot. Required.
#' @param boundary Optional \pkg{sf}/\pkg{sfc} polygon layer to draw as an
#'   outline over the tessellation.
#' @param seeds_sf Optional \pkg{sf}/\pkg{sfc} point layer of seed locations to
#'   overlay.
#' @param features_sf Optional \pkg{sf}/\pkg{sfc} layer of additional features
#'   (points/lines/polygons) to overlay beneath the boundary outline.
#' @param fill_col Optional string giving the name of a column in
#'   \code{tessellation_sf} to map to polygon fill. If \code{NULL} (default) or
#'   the column is absent, tessellation polygons are drawn with outlines only.
#' @param palette Name of the \pkg{viridis} palette to use when \code{fill_col}
#'   is provided. Passed to \code{ggplot2::scale_fill_viridis_c()} for numeric
#'   data and \code{ggplot2::scale_fill_viridis_d()} for non-numeric data.
#'   Default \code{"viridis"}.
#' @param na_fill Fill color used for \code{NA} values in the fill scale.
#'   Default \code{"grey90"}.
#' @param tile_alpha Alpha transparency for tessellation polygons when
#'   \code{fill_col} is mapped. Default \code{0.9}.
#' @param outline_col Outline color for tessellation polygons. Default \code{"white"}.
#' @param outline_size Outline linewidth for tessellation polygons. Default \code{0.2}.
#' @param features_col Color for \code{features_sf}. For points this controls
#'   point color; for lines/polygons it controls stroke color. Default \code{"#333333"}.
#' @param features_size Size/linewidth for \code{features_sf}. Default \code{0.5}.
#' @param seeds_col Color for \code{seeds_sf} points. Default \code{"#1f77b4"}.
#' @param seeds_size Point size for \code{seeds_sf}. Default \code{1.5}.
#' @param boundary_col Outline color for \code{boundary}. Default \code{"#111111"}.
#' @param boundary_size Outline linewidth for \code{boundary}. Default \code{0.6}.
#' @param labels Logical; if \code{TRUE}, draws per-cell text labels at
#'   \code{sf::st_point_on_surface()} of each tessellation polygon. Default \code{FALSE}.
#' @param label_col Name of the tessellation column to use for labels
#'   (only used when \code{labels = TRUE}). Default \code{"grid_id"}.
#' @param label_size Text size for labels. Default \code{2.7}.
#' @param legend Logical; show the fill legend when \code{fill_col} is mapped.
#'   Default \code{TRUE}.
#' @param legend_title Optional title for the fill legend. If \code{NULL},
#'   the legend is titled with \code{fill_col}.
#' @param theme A \pkg{ggplot2} theme to apply. Default \code{ggplot2::theme_void()}.
#' @param target_crs Optional CRS to use for plotting. If supplied, all layers
#'   with a defined CRS are transformed to this CRS. If \code{NULL}, the CRS of
#'   \code{tessellation_sf} (if defined) is used. If no plotting CRS can be
#'   determined, layers are drawn "as-is" (a warning is issued).
#' @param title,subtitle,caption Optional plot annotations passed to
#'   \code{ggplot2::labs()}.
#'
#' @details
#' \strong{CRS handling:}
#' \itemize{
#'   \item If \code{target_crs} is provided, layers with a defined CRS are transformed
#'         to it.
#'   \item Otherwise, if \code{tessellation_sf} has a defined CRS, it is used for plotting
#'         and other layers are transformed to match.
#'   \item Layers lacking a defined CRS are left untransformed (a warning is issued).
#' }
#'
#' \strong{Fill mapping:}
#' When \code{fill_col} is provided and exists in \code{tessellation_sf}, the
#' column is mapped to \code{fill}. Numeric columns use
#' \code{scale\_fill\_viridis\_c(option = palette, na.value = na\_fill)},
#' while non-numeric columns use \code{scale\_fill\_viridis\_d(option = palette,
#' na.value = na\_fill)}. If \code{fill_col} is \code{NULL} or missing in the
#' data, polygons are drawn with outlines only.
#'
#' \strong{Layer order:}
#' Tessellation polygons are drawn first, then \code{features_sf}, then
#' \code{seeds_sf}, and finally the \code{boundary} outline. Labels (if enabled)
#' are drawn on top using \code{geom_sf_text()}.
#'
#' @return A \pkg{ggplot2} object representing the map. The object can be further
#'   modified with additional \pkg{ggplot2} layers or scales.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(ggplot2)
#'
#' # Simple demo tessellation (grid) and boundary
#' b <- st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 10, ymax = 6)), crs = 3857)
#' grid <- st_make_grid(b, n = c(5,3), what = "polygons")
#' tess <- st_sf(grid_id = seq_along(grid), value = rnorm(length(grid)),
#'               geometry = grid)
#' boundary <- st_as_sf(b)
#'
#' # Random seed and feature points
#' set.seed(1)
#' seeds   <- st_as_sf(data.frame(x = runif(10, 0, 10), y = runif(10, 0, 6)),
#'                     coords = c("x","y"), crs = 3857)
#' feats   <- st_as_sf(data.frame(x = runif(50, 0, 10), y = runif(50, 0, 6)),
#'                     coords = c("x","y"), crs = 3857)
#'
#' # Plot with a mapped fill and labels
#' plot_tessellation_map(
#'   tessellation_sf = tess,
#'   boundary        = boundary,
#'   seeds_sf        = seeds,
#'   features_sf     = feats,
#'   fill_col        = "value",
#'   labels          = TRUE,
#'   label_col       = "grid_id",
#'   title           = "Tessellation map",
#'   subtitle        = "with features and seeds"
#' )
#' }
#'
#' @seealso \code{\link{create_voronoi_polygons}}, \code{\link{create_grid_polygons}},
#'   \code{\link{assign_features_to_polygons}}
#'
#' @importFrom sf st_geometry_type st_crs st_transform st_as_sf st_is_empty st_point_on_surface
#' @importFrom ggplot2 ggplot geom_sf geom_sf_text aes labs coord_sf
#' @importFrom ggplot2 scale_fill_viridis_c scale_fill_viridis_d theme_void
#' @export
plot_tessellation_map <- function(tessellation_sf,
                                  boundary = NULL,
                                  seeds_sf = NULL,
                                  features_sf = NULL,
                                  fill_col = NULL,
                                  palette = "viridis",
                                  na_fill = "grey90",
                                  tile_alpha = 0.9,
                                  outline_col = "white",
                                  outline_size = 0.2,
                                  features_col = "#333333",
                                  features_size = 0.5,
                                  seeds_col = "#1f77b4",
                                  seeds_size = 1.5,
                                  boundary_col = "#111111",
                                  boundary_size = 0.6,
                                  labels = FALSE,
                                  label_col = "grid_id",
                                  label_size = 2.7,
                                  legend = TRUE,
                                  legend_title = NULL,
                                  theme = ggplot2::theme_void(),
                                  target_crs = NULL,
                                  title = NULL,
                                  subtitle = NULL,
                                  caption = NULL) {
  
  # --- checks -----------------------------------------------------------------
  if (missing(tessellation_sf) || is.null(tessellation_sf)) {
    stop("plot_tessellation_map(): 'tessellation_sf' is required.")
  }
  if (!inherits(tessellation_sf, "sf")) {
    stop("plot_tessellation_map(): 'tessellation_sf' must be an sf object.")
  }
  gtypes <- as.character(sf::st_geometry_type(tessellation_sf, by_geometry = TRUE))
  if (!all(gtypes %in% c("POLYGON","MULTIPOLYGON"))) {
    stop("plot_tessellation_map(): 'tessellation_sf' must be POLYGON/MULTIPOLYGON.")
  }
  
  .warn <- function(fmt, ...) {
    msg <- sprintf(fmt, ...)
    if (requireNamespace("logger", quietly = TRUE)) logger::log_warn(msg) else warning(msg, call. = FALSE)
  }
  
  # --- pick plot CRS ----------------------------------------------------------
  pick_crs <- function() {
    if (!is.null(target_crs)) return(sf::st_crs(target_crs))
    crs_tess <- sf::st_crs(tessellation_sf)
    if (!is.na(crs_tess)) return(crs_tess)
    sf::st_crs(NA)
  }
  plot_crs <- pick_crs()
  if (is.na(plot_crs)) {
    .warn("plot_tessellation_map(): tessellation CRS is undefined and no 'target_crs' provided; drawing layers as-is.")
  }
  
  # --- coerce & transform helpers --------------------------------------------
  to_sf <- function(x) if (inherits(x, "sfc")) sf::st_as_sf(x) else x
  
  maybe_transform <- function(x) {
    if (is.null(x)) return(NULL)
    x <- to_sf(x)
    if (!inherits(x, "sf") && !inherits(x, "sfc")) {
      stop("plot_tessellation_map(): all layers must be sf/sfc when provided.")
    }
    if (is.na(plot_crs)) return(x)     # nowhere to transform
    x_crs <- sf::st_crs(x)
    if (is.na(x_crs)) {
      .warn("plot_tessellation_map(): a layer has undefined CRS; leaving it untransformed.")
      return(x)
    }
    if (x_crs == plot_crs) return(x)
    sf::st_transform(x, plot_crs)
  }
  
  tess <- maybe_transform(tessellation_sf)
  bnd  <- maybe_transform(boundary)
  sds  <- maybe_transform(seeds_sf)
  fea  <- maybe_transform(features_sf)
  
  # --- build ggplot -----------------------------------------------------------
  p <- ggplot2::ggplot()
  
  # fill mapping (optional)
  has_fill <- !is.null(fill_col) && fill_col %in% names(tess)
  if (has_fill) {
    tess$`..__fill__` <- tess[[fill_col]]
    # tiles with mapped fill: DO NOT set fill=NA here; let aes handle it
    p <- p + ggplot2::geom_sf(
      data = tess,
      mapping = ggplot2::aes(fill = `..__fill__`),
      color = outline_col,
      linewidth = outline_size,
      alpha = tile_alpha
    )
  } else {
    # outlines only
    p <- p + ggplot2::geom_sf(
      data = tess,
      color = outline_col,
      linewidth = outline_size,
      fill = NA
    )
  }
  
  # features below boundary
  if (!is.null(fea) && nrow(fea) > 0L && !all(sf::st_is_empty(fea))) {
    p <- p + ggplot2::geom_sf(
      data = fea,
      color = features_col,
      linewidth = features_size, # used by LINESTRING/POLYGON
      size = features_size       # used by POINTs
    )
  }
  
  # seeds
  if (!is.null(sds) && nrow(sds) > 0L && !all(sf::st_is_empty(sds))) {
    p <- p + ggplot2::geom_sf(
      data = sds,
      color = seeds_col,
      size  = seeds_size
    )
  }
  
  # boundary outline
  if (!is.null(bnd) && nrow(bnd) > 0L && !all(sf::st_is_empty(bnd))) {
    p <- p + ggplot2::geom_sf(
      data = bnd,
      fill = NA,
      color = boundary_col,
      linewidth = boundary_size
    )
  }
  
  # labels
  if (isTRUE(labels)) {
    if (!label_col %in% names(tess)) {
      .warn("plot_tessellation_map(): 'labels=TRUE' but label_col '%s' not found; skipping labels.", label_col)
    } else if (nrow(tess) > 0L && !all(sf::st_is_empty(tess))) {
      centers <- suppressWarnings(sf::st_point_on_surface(tess))
      centers$`..__lab__` <- tess[[label_col]]
      p <- p + ggplot2::geom_sf_text(
        data = centers,
        ggplot2::aes(label = `..__lab__`),
        size = label_size
      )
    }
  }
  
  # fill scale (if mapped)
  if (has_fill) {
    label <- if (!is.null(legend_title)) legend_title else fill_col
    is_cont <- is.numeric(tess$`..__fill__`)
    if (is_cont) {
      p <- p + ggplot2::scale_fill_viridis_c(
        option = palette, na.value = na_fill,
        name = label, guide = if (legend) "colourbar" else "none"
      )
    } else {
      p <- p + ggplot2::scale_fill_viridis_d(
        option = palette, na.value = na_fill,
        name = label, guide = if (legend) "legend" else "none", drop = FALSE
      )
    }
  }
  
  # coords, theme, labels
  p <- p +
    ggplot2::coord_sf(crs = if (!is.na(plot_crs)) plot_crs else NULL, expand = FALSE) +
    theme +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption)
  
  p
}

#' Create spatial cross-validation folds (random K-fold, block K-fold, buffered LOO)
#'
#' Builds train/test splits for point data using one of three strategies:
#' random \emph{K}-fold, spatial block \emph{K}-fold, or buffered leave-one-out
#' (LOO). The function is CRS-aware, coerces non-point geometries to points when
#' needed, and uses a projected CRS for distance-based operations.
#'
#' @param points_sf An \pkg{sf} object. If geometries are not
#'   \code{POINT}/\code{MULTIPOINT}, they are coerced to points via
#'   \code{\link{coerce_to_points}(x, "auto")} before fold construction. Rows are
#'   not reordered.
#' @param k Integer. Number of folds for \code{method = "random_kfold"} and
#'   \code{method = "block_kfold"}. If \code{k > nrow(points\_sf)}, it is reduced
#'   to \code{nrow(points\_sf)}. Ignored for \code{method = "buffered_loo"}
#'   (which produces \eqn{n} folds).
#' @param method Character; one of:
#'   \itemize{
#'     \item \code{"random_kfold"} — random partition into \code{k} folds of
#'       approximately equal size.
#'     \item \code{"block_kfold"} — build a grid of spatial blocks over a
#'       region and assign entire blocks to folds (spatial CV).
#'     \item \code{"buffered_loo"} — leave-one-out where each test point is
#'       buffered and all training points within \code{buffer} of the test point
#'       (including itself) are excluded.
#'   }
#' @param seed Optional integer RNG seed. If supplied, the seed is set for the
#'   duration of the function and the previous RNG state is restored on exit.
#' @param block_nx,block_ny Optional integers giving the number of grid cells in
#'   the \emph{x} and \emph{y} directions for \code{method = "block_kfold"}.
#'   When \code{NULL}, these are derived automatically from the region aspect
#'   ratio and \code{block_multiplier * k}.
#' @param block_multiplier Numeric scalar (default \code{3}). Only used when
#'   \code{block_nx/block_ny} are not provided. The target number of blocks is
#'   approximately \code{block_multiplier * k}, and \code{block_nx, block_ny}
#'   are chosen to respect the region's width/height ratio.
#' @param boundary Optional \pkg{sf}/\pkg{sfc} polygon layer defining the region
#'   for \code{method = "block_kfold"}. If provided, it is transformed to the
#'   (projected) CRS of the points, validated, unioned, and used to clip the
#'   grid. If \code{NULL}, the grid is built over the points' bounding box. If
#'   some points fall outside the supplied boundary, the function unions the
#'   boundary with the points' bbox to ensure all points can be assigned.
#' @param buffer Positive numeric distance (in CRS units) used by
#'   \code{method = "buffered_loo"} to exclude neighbors around each test point.
#'   Required when \code{method = "buffered_loo"}.
#' @param drop_empty_blocks Logical (default \code{TRUE}). For
#'   \code{method = "block_kfold"}, drop grid cells that contain no points and
#'   remap block IDs before allocating blocks to folds.
#'
#' @details
#' \strong{CRS and geometry handling}
#' \itemize{
#'   \item If \code{points_sf} is not \code{POINT}/\code{MULTIPOINT}, it is
#'         converted to points with \code{\link{coerce_to_points}(x, "auto")}.
#'   \item Distance-based operations (\code{"block_kfold"} grid sizing/centroids
#'         and \code{"buffered_loo"} neighbor search) are performed in a
#'         projected CRS: points are transformed with
#'         \code{\link{ensure_projected}()} (and \code{boundary} is transformed
#'         to match when provided).
#' }
#'
#' \strong{Method specifics}
#' \describe{
#'   \item{\code{random_kfold}}{Shuffles row indices and splits them into
#'     \code{k} approximately equal folds.}
#'   \item{\code{block_kfold}}{Builds a rectangular grid over the \code{boundary}
#'     (or points' bbox), clips it to the region, assigns each point to the
#'     first intersecting grid cell, optionally drops empty cells, and greedily
#'     allocates whole blocks to folds to balance point counts. If the number of
#'     non-empty blocks is less than \code{k}, \code{k} is reduced accordingly.
#'     Points not intersecting any block are assigned to the nearest grid
#'     centroid.}
#'   \item{\code{buffered_loo}}{For each observation \eqn{i}, the test set is
#'     \eqn{\{i\}} and the training set excludes all observations within
#'     \code{buffer} of \eqn{i} (including \eqn{i}). Produces \eqn{n} folds.}
#' }
#'
#' \strong{Reproducibility} — If \code{seed} is supplied, the RNG state is set
#' within the function and restored on exit.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{method}}{The method used (character).}
#'   \item{\code{k}}{Number of folds actually produced (integer). For
#'     \code{"buffered_loo"} this equals \eqn{n}.}
#'   \item{\code{folds}}{A list of length \code{k}; each element is a list with
#'     integer vectors \code{$train} and \code{$test} giving row indices (1-based)
#'     into the input \code{points_sf} after internal coercion/projection (row
#'     order is preserved).}
#'   \item{\code{assignment}}{A tibble with columns
#'     \code{row_id} (original input row index) and \code{fold} (the assigned
#'     fold number for that row).}
#'   \item{\code{params}}{A list of method-specific settings used (e.g., RNG
#'     \code{seed}, derived grid dimensions, buffer distance, whether a
#'     \code{boundary} was supplied).}
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' # Make some points in a rectangle (projected CRS)
#' set.seed(1)
#' pts <- st_as_sf(
#'   data.frame(x = runif(100, 0, 1000), y = runif(100, 0, 600)),
#'   coords = c("x","y"), crs = 3857
#' )
#'
#' # 1) Random K-fold
#' folds_rnd <- make_folds(pts, k = 5, method = "random_kfold", seed = 42)
#' str(folds_rnd$folds[[1]])
#'
#' # 2) Block K-fold over a boundary polygon
#' bnd <- st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 1000, ymax = 600)), crs = 3857)
#' bnd <- st_sf(geometry = bnd)
#' folds_blk <- make_folds(
#'   pts, k = 5, method = "block_kfold",
#'   boundary = bnd, block_multiplier = 3, seed = 42
#' )
#'
#' # 3) Buffered leave-one-out (buffer in CRS units)
#' folds_bloo <- make_folds(pts, k = 0, method = "buffered_loo", buffer = 100)
#' }
#'
#' @seealso \code{\link{coerce_to_points}}, \code{\link{ensure_projected}}
#'
#' @importFrom sf st_geometry_type st_make_valid st_union st_intersects
#' @importFrom sf st_as_sfc st_bbox st_set_crs st_sf st_make_grid st_intersection
#' @importFrom sf st_centroid st_distance st_is_within_distance
#' @importFrom tibble tibble
#' @export
make_folds <- function(points_sf,
                       k,
                       method = c("random_kfold","block_kfold","buffered_loo"),
                       seed = NULL,
                       block_nx = NULL,
                       block_ny = NULL,
                       block_multiplier = 3,
                       boundary = NULL,
                       buffer = NULL,
                       drop_empty_blocks = TRUE) {
  method <- match.arg(method)
  
  # ----- lightweight logging wrappers (no hard {logger} dependency) ----------
  .log_warn <- function(fmt, ...) {
    msg <- sprintf(fmt, ...)
    if (requireNamespace("logger", quietly = TRUE)) logger::log_warn(msg) else warning(msg, call. = FALSE)
  }
  .log_info <- function(fmt, ...) {
    msg <- sprintf(fmt, ...)
    if (requireNamespace("logger", quietly = TRUE)) logger::log_info(msg) else message(msg)
  }
  
  # ----- seed hygiene: set locally and restore on exit ------------------------
  if (!is.null(seed)) {
    old_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (old_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    on.exit({
      if (old_exists) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
  }
  
  if (!inherits(points_sf, "sf")) stop("make_folds(): `points_sf` must be an sf object.")
  # Ensure POINTs (preserve attributes) and a projected CRS where needed
  if (!all(sf::st_geometry_type(points_sf, by_geometry = TRUE) %in% c("POINT","MULTIPOINT"))) {
    points_sf <- coerce_to_points(points_sf, "auto")
  }
  # Keep an index to refer back to original row order
  points_sf$..row_id <- seq_len(nrow(points_sf))
  
  # Helper to build return object
  .ret <- function(method, k, folds, assignment, params) {
    list(
      method = method,
      k = k,
      folds = folds,
      assignment = assignment,
      params = params
    )
  }
  
  # ------------------------------
  # RANDOM K-FOLD
  # ------------------------------
  if (method == "random_kfold") {
    n <- nrow(points_sf)
    if (k < 2) k <- 2L
    if (k > n) {
      .log_warn("make_folds(random_kfold): k (%d) > n (%d); reducing k to n.", k, n)
      k <- n
    }
    idx <- sample.int(n, n)  # shuffled indices
    # Split approximately evenly
    sizes <- rep(floor(n / k), k)
    remainder <- n - sum(sizes)
    if (remainder > 0) sizes[seq_len(remainder)] <- sizes[seq_len(remainder)] + 1L
    
    splits <- vector("list", k)
    start <- 1L
    assign_vec <- integer(n)
    for (j in seq_len(k)) {
      stop <- start + sizes[j] - 1L
      test_idx <- idx[start:stop]
      train_idx <- setdiff(idx, test_idx)
      splits[[j]] <- list(train = train_idx, test = test_idx)
      assign_vec[test_idx] <- j
      start <- stop + 1L
    }
    assignment <- tibble::tibble(row_id = points_sf$..row_id, fold = assign_vec)
    return(.ret(
      method = method,
      k = k,
      folds = splits,
      assignment = assignment,
      params = list(seed = seed)
    ))
  }
  
  # ------------------------------
  # BLOCK K-FOLD (spatial blocks)
  # ------------------------------
  if (method == "block_kfold") {
    if (k < 2) k <- 2L
    
    # Work in a projected CRS
    pts <- ensure_projected(points_sf)
    reg <- if (!is.null(boundary)) {
      b <- ensure_projected(boundary, sf::st_crs(pts))
      if (inherits(b, "sfc")) b <- sf::st_sf(geometry = b)
      b <- sf::st_make_valid(sf::st_union(b))
      # Robust inside check for multi-part boundaries:
      mat <- sf::st_intersects(pts, b, sparse = FALSE)
      inside_any <- apply(mat, 1L, any)
      if (!all(inside_any)) {
        bb_pts <- sf::st_as_sfc(sf::st_bbox(pts)) |> sf::st_set_crs(sf::st_crs(pts))
        b <- suppressWarnings(sf::st_union(sf::st_make_valid(sf::st_sf(geometry = c(b, bb_pts)))))
      }
      b
    } else {
      sf::st_as_sfc(sf::st_bbox(pts)) |> sf::st_set_crs(sf::st_crs(pts)) |> sf::st_sf()
    }
    
    # Derive grid dims if not provided
    if (is.null(block_nx) || is.null(block_ny)) {
      bb <- sf::st_bbox(reg)
      w  <- as.numeric(bb["xmax"] - bb["xmin"])
      h  <- as.numeric(bb["ymax"] - bb["ymin"])
      ratio <- if (h > 0) w / h else 1
      target_blocks <- max(1L, round(block_multiplier * k))
      nx <- max(1L, round(sqrt(target_blocks * ratio)))
      ny <- max(1L, round(max(1, target_blocks / nx)))
    } else {
      nx <- as.integer(block_nx); ny <- as.integer(block_ny)
      if (nx < 1 || ny < 1) stop("make_folds(block_kfold): block_nx and block_ny must be >= 1.")
    }
    
    .log_info("make_folds(block_kfold): making grid %d x %d (≈ %d blocks) over region.", nx, ny, nx*ny)
    
    grid <- sf::st_make_grid(reg, n = c(nx, ny), what = "polygons", square = TRUE)
    # Clip to region outline; drop empties
    grid <- suppressWarnings(sf::st_intersection(grid, sf::st_union(reg)))
    grid_sf <- sf::st_as_sf(grid)
    if (nrow(grid_sf) == 0) stop("make_folds(block_kfold): generated grid has zero cells after clipping.")
    
    # Assign each point to the first intersecting grid cell
    hits <- sf::st_intersects(pts, grid_sf)
    block_id <- vapply(hits, function(ix) if (length(ix)) ix[1] else NA_integer_, 1L)
    pts$..block_id <- block_id
    
    if (drop_empty_blocks) {
      # Keep only blocks that actually contain points
      used_blocks <- sort(unique(pts$..block_id[!is.na(pts$..block_id)]))
      grid_sf <- grid_sf[used_blocks, , drop = FALSE]
      # Re-map block ids to 1..B
      remap <- match(pts$..block_id, used_blocks)
      pts$..block_id <- remap
    }
    
    # Handle points outside all blocks (rare, but be safe)
    if (anyNA(pts$..block_id)) {
      n_na <- sum(is.na(pts$..block_id))
      .log_warn("make_folds(block_kfold): %d points not assigned to any block; assigning to nearest grid centroid.", n_na)
      cent <- sf::st_centroid(sf::st_geometry(grid_sf))
      na_idx <- which(is.na(pts$..block_id))
      dmat <- sf::st_distance(sf::st_geometry(pts[na_idx, ]), cent)
      pts$..block_id[na_idx] <- apply(as.matrix(dmat), 1, which.min)
    }
    
    B <- max(pts$..block_id, na.rm = TRUE)
    if (!is.finite(B) || B < 1) stop("make_folds(block_kfold): no usable blocks found.")
    if (B < k) {
      .log_warn("make_folds(block_kfold): number of non-empty blocks (%d) < k (%d); reducing k to %d.", B, k, B)
      k <- B
    }
    
    # Count points per block and sort descending
    blk_sizes <- as.integer(table(factor(pts$..block_id, levels = seq_len(B))))
    order_blk <- order(blk_sizes, decreasing = TRUE)
    blk_ids_sorted <- seq_len(B)[order_blk]
    blk_sizes_sorted <- blk_sizes[order_blk]
    
    # Greedy allocation of blocks to folds to balance point counts
    fold_loads <- integer(k)
    fold_blocks <- vector("list", k)
    for (i in seq_along(blk_ids_sorted)) {
      # Choose fold with minimal current load; tie-break with RNG (seeded above if provided)
      j <- which(fold_loads == min(fold_loads))
      if (length(j) > 1) j <- sample(j, 1L)
      fold_blocks[[j]] <- c(fold_blocks[[j]], blk_ids_sorted[i])
      fold_loads[j] <- fold_loads[j] + blk_sizes_sorted[i]
    }
    
    # Build train/test index sets per fold
    splits <- vector("list", k)
    assign_vec <- integer(nrow(pts))
    for (j in seq_len(k)) {
      test_idx <- which(pts$..block_id %in% fold_blocks[[j]])
      train_idx <- setdiff(seq_len(nrow(pts)), test_idx)
      splits[[j]] <- list(train = train_idx, test = test_idx)
      assign_vec[test_idx] <- j
    }
    
    assignment <- tibble::tibble(row_id = pts$..row_id, fold = assign_vec)
    
    return(.ret(
      method = method,
      k = k,
      folds = splits,
      assignment = assignment,
      params = list(
        seed = seed,
        grid_nx = nx,
        grid_ny = ny,
        blocks_used = B,
        block_multiplier = block_multiplier,
        boundary_supplied = !is.null(boundary)
      )
    ))
  }
  
  # ------------------------------
  # BUFFERED LEAVE-ONE-OUT
  # ------------------------------
  if (method == "buffered_loo") {
    if (is.null(buffer) || !is.numeric(buffer) || buffer <= 0) {
      stop("make_folds(buffered_loo): `buffer` (positive numeric, in CRS units) is required.")
    }
    # Work in projected CRS for reliable metric distances
    pts <- ensure_projected(points_sf)
    n <- nrow(pts)
    
    .log_info("make_folds(buffered_loo): building neighbor lists within buffer = %g (CRS units).", buffer)
    # For each i, list of j within 'buffer' (includes i)
    nb <- sf::st_is_within_distance(sf::st_geometry(pts), sf::st_geometry(pts), dist = buffer)
    
    splits <- vector("list", n)
    assignment <- tibble::tibble(row_id = pts$..row_id, fold = seq_len(n))
    for (i in seq_len(n)) {
      excl <- sort(unique(nb[[i]]))
      test_idx  <- i
      train_idx <- setdiff(seq_len(n), excl)
      splits[[i]] <- list(train = train_idx, test = test_idx)
    }
    
    return(.ret(
      method = method,
      k = n,
      folds = splits,
      assignment = assignment,
      params = list(buffer = buffer)
    ))
  }
  
  stop("make_folds(): unsupported method.")
}

#' K-fold or user-supplied cross-validation for Geographically Weighted Regression (GWR)
#'
#' Runs out-of-sample evaluation for a GWR model using either internally created
#' random \emph{K}-folds or a user-supplied fold definition. Data are prepared
#' via \code{\link{prep_model_data}()}, coerced to points when needed, and must
#' be in a projected CRS (metric units) before fitting. Bandwidth can be fixed
#' (adaptive fraction) or selected per fold via \code{spgwr::gwr.sel()}.
#'
#' @param data_sf An \pkg{sf} object containing the features to model. If the
#'   geometry is not \code{POINT}/\code{MULTIPOINT}, it is coerced to points
#'   according to \code{pointize} by \code{\link{prep_model_data}()}.
#' @param response_var Character scalar. Name of the response variable column in
#'   \code{data_sf}.
#' @param predictor_vars Character vector. Names of predictor columns in
#'   \code{data_sf}. Use an empty vector for an intercept-only model.
#' @param folds Optional list of folds. Each element must be a list with integer
#'   vectors \code{$train} and \code{$test} giving \emph{original} row indices
#'   (before preprocessing). If omitted, random \emph{K}-folds are created on the
#'   kept rows after preprocessing using \code{k} and \code{seed}. You can
#'   obtain a compatible structure from \code{\link{make_folds}()} (use the
#'   \code{$folds} element or the \code{$assignment} to build one).
#' @param k Integer. Number of folds to create when \code{folds} is \code{NULL}.
#'   Ignored when \code{folds} is supplied.
#' @param seed Integer RNG seed used when creating random folds. Ignored when
#'   \code{folds} is supplied.
#' @param bandwidth Optional numeric in \eqn{(0,1)}: adaptive bandwidth fraction
#'   (proportion of neighbors) for GWR. If \code{NULL} (default), the bandwidth
#'   is selected on the training subset of each fold via
#'   \code{spgwr::gwr.sel(adapt = TRUE)}. Values \eqn{\le 0} trigger per-fold
#'   selection; values \eqn{\ge 1} are clamped to \code{0.99}.
#' @param boundary Optional \pkg{sf}/\pkg{sfc} polygon used only to help
#'   \code{\link{prep_model_data}()} choose/align a projected CRS. It is not
#'   used directly by \code{cv_gwr()}.
#' @param pointize Character; how to convert non-point geometries if present in
#'   \code{data_sf}. Passed to \code{\link{prep_model_data}()} (see that
#'   function for supported modes such as \code{"auto"}, \code{"surface"},
#'   \code{"centroid"}, etc.).
#'
#' @details
#' \strong{Preprocessing:} The function requires \code{\link{prep_model_data}()}
#' to be available. That helper (a) validates columns, (b) coerces geometry to
#' points when necessary, (c) ensures a projected CRS, and (d) drops rows with
#' missing/non-finite modeling values. A column \code{..row_id} (the original
#' row index) is preserved and used to remap user-supplied folds to the kept
#' rows.
#'
#' \strong{CRS:} GWR is distance-based. After preprocessing, the input must be in
#' a projected CRS; the function stops if data remain in lon/lat.
#'
#' \strong{Bandwidth:} When \code{bandwidth} is \code{NULL}, an adaptive
#' fraction is selected on each training fold via \code{spgwr::gwr.sel()}.
#' Otherwise, the supplied adaptive fraction is used for fitting and predicting.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{overall}}{Data frame with aggregate RMSE, MAE, R2 across all
#'   test predictions (\code{n\_pred} = number of finite predictions).}
#'   \item{\code{fold_metrics}}{Data frame of per-fold metrics: \code{fold},
#'   \code{n\_test}, \code{n\_pred}, \code{RMSE}, \code{MAE}, \code{R2},
#'   \code{bandwidth} (adaptive fraction used), and training-only \code{train\_AICc}
#'   when available.}
#'   \item{\code{predictions}}{Data frame with columns \code{..row_id} (original
#'   row index), \code{fold}, observed \code{y}, and predicted \code{yhat}.}
#'   \item{\code{folds}}{The fold definition actually used (list of \code{$train}
#'   $test} index vectors) after remapping to kept rows.}
#'   \item{\code{formula}}{The model formula used.}
#'   \item{\code{adaptive}}{Logical flag (always \code{TRUE} here; GWR is fit
#'   with adaptive bandwidth).}
#' }
#'
#' @section Failure handling:
#' Bandwidth selection or prediction failures on a fold are handled
#' defensively: the bandwidth falls back to \code{0.75} (adaptive fraction), and
#' predictions revert to computing \eqn{X \beta} from fold-specific local
#' coefficients when \code{spgwr::gwr.predict()} does not provide a prediction
#' column. Folds with no finite predictions are skipped (and reported in
#' \code{fold_metrics} only when metrics can be computed).
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' # Toy projected data
#' set.seed(1)
#' pts <- st_as_sf(
#'   data.frame(
#'     x = runif(200, 0, 2000),
#'     y = runif(200, 0, 1200),
#'     y_resp = rnorm(200),
#'     x1 = rnorm(200),
#'     x2 = rnorm(200)
#'   ),
#'   coords = c("x","y"), crs = 3857
#' )
#'
#' # Random 5-fold CV (bandwidth selected per fold)
#' res_cv <- cv_gwr(
#'   data_sf = pts,
#'   response_var = "y_resp",
#'   predictor_vars = c("x1","x2"),
#'   k = 5, seed = 42
#' )
#' res_cv$overall
#'
#' # With user folds from make_folds()
#' fdef <- make_folds(pts, k = 5, method = "random_kfold", seed = 42)
#' res_cv2 <- cv_gwr(
#'   data_sf = pts,
#'   response_var = "y_resp",
#'   predictor_vars = c("x1","x2"),
#'   folds = fdef$folds
#' )
#' }
#'
#' @seealso \code{\link{fit_gwr_model}}, \code{\link{make_folds}},
#'   \code{\link{prep_model_data}}, \code{\link{evaluate_models}}
#'
#' @importFrom sf st_geometry_type st_is_longlat st_drop_geometry
#' @importFrom stats as.formula model.matrix na.omit
#' @importFrom methods as
#' @importFrom spgwr gwr.sel gwr.predict gwr
#' @export
cv_gwr <- function(
    data_sf,
    response_var,
    predictor_vars,
    folds     = NULL,
    k         = 5,
    seed      = 123,
    bandwidth = NULL,
    boundary  = NULL,
    pointize  = "auto"
) {
  if (!inherits(data_sf, "sf")) stop("cv_gwr(): `data_sf` must be an sf object.")
  if (!requireNamespace("spgwr", quietly = TRUE)) stop("cv_gwr(): package 'spgwr' is required.")
  if (!requireNamespace("sp", quietly = TRUE))    stop("cv_gwr(): package 'sp' is required.")
  
  # Preserve original row ids for remapping and run centralized prep
  if (!("..row_id" %in% names(data_sf))) data_sf$`..row_id` <- seq_len(nrow(data_sf))
  if (!exists("prep_model_data", mode = "function")) {
    stop("cv_gwr(): `prep_model_data()` not found. Please define it and try again.")
  }
  
  dat_sf <- prep_model_data(
    data_sf        = data_sf,
    response_var   = response_var,
    predictor_vars = predictor_vars,
    boundary       = boundary,
    pointize       = pointize
  )
  
  if (!("..row_id" %in% names(dat_sf))) {
    stop("cv_gwr(): `prep_model_data()` must preserve a `..row_id` column.")
  }
  
  # Ensure POINT/MULTIPOINT and projected (prep is expected to do this)
  gtypes <- as.character(sf::st_geometry_type(dat_sf, by_geometry = TRUE))
  if (!all(gtypes %in% c("POINT", "MULTIPOINT"))) {
    stop("cv_gwr(): geometry must be POINT/MULTIPOINT after prep_model_data().")
  }
  is_ll <- tryCatch(sf::st_is_longlat(dat_sf), error = function(e) FALSE)
  if (isTRUE(is_ll)) {
    stop("cv_gwr(): data must be in a projected CRS after prep_model_data().")
  }
  
  keep_idx <- dat_sf$`..row_id`
  
  # Build formula
  fml <- if (length(predictor_vars) == 0L) {
    stats::as.formula(paste(response_var, "~ 1"))
  } else {
    stats::as.formula(paste(response_var, "~", paste(predictor_vars, collapse = " + ")))
  }
  
  # Folds: use provided (remap) or create at random on kept rows
  remapped_folds <- list()
  if (is.null(folds)) {
    set.seed(seed)
    ids <- keep_idx
    fold_id <- sample(rep(seq_len(k), length.out = length(ids)))
    for (i in seq_len(k)) {
      remapped_folds[[i]] <- list(train = ids[fold_id != i], test = ids[fold_id == i])
    }
  } else {
    for (i in seq_along(folds)) {
      tr_pos <- match(folds[[i]]$train, keep_idx)
      te_pos <- match(folds[[i]]$test,  keep_idx)
      tr_ids <- keep_idx[stats::na.omit(tr_pos)]
      te_ids <- keep_idx[stats::na.omit(te_pos)]
      remapped_folds[[i]] <- list(train = tr_ids, test = te_ids)
    }
  }
  
  # Logger helper
  .log_warn <- function(...) {
    msg <- sprintf(...)
    if (requireNamespace("logger", quietly = TRUE)) logger::log_warn(msg) else warning(msg)
  }
  
  # Validate / clamp user-supplied adaptive fraction
  use_fixed_bw <- !is.null(bandwidth)
  if (use_fixed_bw) {
    if (!is.finite(bandwidth) || is.na(bandwidth) || bandwidth <= 0) {
      .log_warn("cv_gwr(): invalid `bandwidth`; falling back to per-fold selection.")
      use_fixed_bw <- FALSE
    } else if (bandwidth >= 1) {
      .log_warn("cv_gwr(): `bandwidth` >= 1; clamping to 0.99.")
      bandwidth <- 0.99
    }
  }
  
  # Storage
  pred_rows  <- list()
  fold_stats <- list()
  
  for (i in seq_along(remapped_folds)) {
    tr_ids <- remapped_folds[[i]]$train
    te_ids <- remapped_folds[[i]]$test
    
    tr_pos <- stats::na.omit(match(tr_ids, keep_idx))
    te_pos <- stats::na.omit(match(te_ids, keep_idx))
    
    if (length(tr_pos) < 2L || length(te_pos) < 1L) {
      .log_warn("cv_gwr(): fold %d skipped (train n=%d, test n=%d).", i, length(tr_pos), length(te_pos))
      next
    }
    
    train_sf <- dat_sf[tr_pos, , drop = FALSE]
    test_sf  <- dat_sf[te_pos, , drop = FALSE]
    
    sp_train <- methods::as(train_sf, "Spatial")
    sp_test  <- methods::as(test_sf,  "Spatial")
    
    # Bandwidth selection on training set (adaptive)
    bw_i <- bandwidth
    if (!use_fixed_bw) {
      bw_try <- try(suppressWarnings(spgwr::gwr.sel(fml, data = sp_train, adapt = TRUE)), silent = TRUE)
      if (inherits(bw_try, "try-error") || !is.finite(bw_try) || is.na(bw_try) || bw_try <= 0) {
        .log_warn("cv_gwr(): fold %d bandwidth selection failed; using 0.75.", i)
        bw_i <- 0.75
      } else if (bw_try >= 1) {
        .log_warn("cv_gwr(): fold %d bandwidth >= 1; clamping to 0.99.", i)
        bw_i <- 0.99
      } else {
        bw_i <- as.numeric(bw_try)
      }
    } else {
      bw_i <- as.numeric(bandwidth)
    }
    
    # Training-only fit to capture AICc (optional but useful)
    train_AICc <- NA_real_
    fit_train <- try(
      suppressWarnings(spgwr::gwr(fml, data = sp_train, adapt = bw_i, hatmatrix = TRUE, se.fit = FALSE)),
      silent = TRUE
    )
    if (!inherits(fit_train, "try-error") && !is.null(fit_train$GW.diagnostic)) {
      if (!is.null(fit_train$GW.diagnostic$AICc)) {
        train_AICc <- as.numeric(fit_train$GW.diagnostic$AICc)
      }
    }
    
    # Predict on test locations (primary path)
    y_true <- sf::st_drop_geometry(test_sf)[[response_var]]
    y_hat  <- rep(NA_real_, length(te_pos))
    
    pred_obj <- try(
      suppressWarnings(
        spgwr::gwr.predict(
          formula     = fml,
          data        = sp_train,
          predictdata = sp_test,
          adapt       = bw_i,
          se.fit      = FALSE
        )
      ),
      silent = TRUE
    )
    
    if (!inherits(pred_obj, "try-error") && is.list(pred_obj) && !is.null(pred_obj$SDF)) {
      sdf <- pred_obj$SDF
      pred_df <- if (methods::is(sdf, "Spatial")) sdf@data else as.data.frame(sdf)
      nm  <- names(pred_df)
      hit <- nm[grepl("^(pred|prediction|yhat|fit)", nm, ignore.case = TRUE)]
      if (length(hit)) {
        y_hat <- as.numeric(pred_df[[hit[1]]])
      } else {
        # Fallback: recompute as X * local betas at fit points (if present)
        coef_df <- pred_df
        mm      <- stats::model.matrix(fml, data = sf::st_drop_geometry(test_sf))
        cn_mm <- make.names(colnames(mm))
        cn_cf <- make.names(colnames(coef_df))
        keep  <- intersect(cn_mm, cn_cf)
        if (length(keep)) {
          mm2 <- mm[, match(keep, cn_mm), drop = FALSE]
          cf2 <- as.matrix(coef_df[, match(keep, cn_cf), drop = FALSE])
          if (nrow(mm2) == nrow(cf2)) y_hat <- rowSums(mm2 * cf2)
        }
      }
    } else {
      # Robust backup: fit at test points, then X * betas
      fit_obj <- try(
        suppressWarnings(spgwr::gwr(fml, data = sp_train, adapt = bw_i,
                                    fit.points = sp_test, hatmatrix = FALSE)),
        silent = TRUE
      )
      if (!inherits(fit_obj, "try-error") && !is.null(fit_obj$SDF)) {
        coef_df <- fit_obj$SDF@data
        mm      <- stats::model.matrix(fml, data = sf::st_drop_geometry(test_sf))
        cn_mm <- make.names(colnames(mm))
        cn_cf <- make.names(colnames(coef_df))
        keep  <- intersect(cn_mm, cn_cf)
        if (length(keep)) {
          mm2 <- mm[, match(keep, cn_mm), drop = FALSE]
          cf2 <- as.matrix(coef_df[, match(keep, cn_cf), drop = FALSE])
          if (nrow(mm2) == nrow(cf2)) y_hat <- rowSums(mm2 * cf2)
        }
      } else {
        .log_warn("cv_gwr(): fold %d prediction failed; yhat set to NA.", i)
      }
    }
    
    # Metrics for this fold
    ok <- is.finite(y_true) & is.finite(y_hat)
    n_ok <- sum(ok)
    if (n_ok == 0L) {
      .log_warn("cv_gwr(): fold %d produced no finite predictions.", i)
      next
    }
    
    resid <- y_true[ok] - y_hat[ok]
    rmse  <- sqrt(mean(resid^2))
    mae   <- mean(abs(resid))
    ssr   <- sum(resid^2)
    sst   <- sum((y_true[ok] - mean(y_true[ok]))^2)
    r2    <- if (sst > 0) 1 - ssr / sst else NA_real_
    
    fold_stats[[length(fold_stats) + 1L]] <- data.frame(
      fold = i, n_test = length(y_true), n_pred = n_ok,
      RMSE = rmse, MAE = mae, R2 = r2, bandwidth = bw_i, train_AICc = train_AICc,
      stringsAsFactors = FALSE
    )
    
    pred_rows[[length(pred_rows) + 1L]] <- data.frame(
      `..row_id` = test_sf$`..row_id`,
      fold       = i,
      y          = as.numeric(y_true),
      yhat       = as.numeric(y_hat),
      stringsAsFactors = FALSE
    )
  }
  
  # Bind results
  preds <- if (length(pred_rows)) do.call(rbind, pred_rows) else
    data.frame(`..row_id` = integer(), fold = integer(), y = numeric(), yhat = numeric())
  
  folds_df <- if (length(fold_stats)) do.call(rbind, fold_stats) else
    data.frame(fold = integer(), n_test = integer(), n_pred = integer(),
               RMSE = numeric(), MAE = numeric(), R2 = numeric(),
               bandwidth = numeric(), train_AICc = numeric())
  
  # Overall metrics across all folds
  ok_all <- is.finite(preds$y) & is.finite(preds$yhat)
  if (any(ok_all)) {
    resid_all <- preds$y[ok_all] - preds$yhat[ok_all]
    rmse_all  <- sqrt(mean(resid_all^2))
    mae_all   <- mean(abs(resid_all))
    ssr_all   <- sum(resid_all^2)
    sst_all   <- sum((preds$y[ok_all] - mean(preds$y[ok_all]))^2)
    r2_all    <- if (sst_all > 0) 1 - ssr_all / sst_all else NA_real_
    overall   <- data.frame(RMSE = rmse_all, MAE = mae_all, R2 = r2_all,
                            n_pred = sum(ok_all), stringsAsFactors = FALSE)
  } else {
    overall <- data.frame(RMSE = NA_real_, MAE = NA_real_, R2 = NA_real_,
                          n_pred = 0L, stringsAsFactors = FALSE)
  }
  
  list(
    overall      = overall,
    fold_metrics = folds_df,
    predictions  = preds,
    folds        = remapped_folds,
    formula      = deparse(fml),
    adaptive     = TRUE
  )
}

#' K-fold cross-validation for the Bayesian spatial model (brms / rstanarm backends)
#'
#' Performs out-of-sample evaluation for a Bayesian spatial regression model
#' fitted by \code{\link{fit_bayesian_spatial_model}()} across either user-supplied
#' folds or internally generated random \emph{K}-folds. Geometry is validated and
#' prepared via \code{\link{prep_model_data}()} and predictions are produced from
#' posterior draws (preferred) or generic \code{predict()} as a fallback.
#'
#' @param data_sf An \pkg{sf} object with observations and geometry. If geometry
#'   is not \code{POINT}/\code{MULTIPOINT}, it will be coerced to points by
#'   \code{\link{prep_model_data}()} according to \code{pointize}.
#' @param response_var Character scalar. Name of the response column in \code{data_sf}.
#' @param predictor_vars Character vector of predictor column names in \code{data_sf}.
#'   Use an empty vector for an intercept-only model.
#' @param folds Optional list of folds. Each element must be a list with integer
#'   vectors \code{$train} and \code{$test} that refer to the \emph{original}
#'   row indices in \code{data_sf} (before preprocessing). If omitted, random
#'   \emph{K}-folds are created on the kept rows after preprocessing using \code{k}
#'   and \code{seed}. You can obtain a compatible structure from
#'   \code{\link{make_folds}()} (use its \code{$folds} element).
#' @param k Integer \eqn{\ge 2}. Number of folds to create when \code{folds} is
#'   \code{NULL}. Ignored otherwise.
#' @param seed Integer RNG seed used when creating random folds. Ignored when
#'   \code{folds} is supplied.
#' @param boundary Optional \pkg{sf}/\pkg{sfc} polygon used only by
#'   \code{\link{prep_model_data}()} to help choose/align a projected CRS. Not
#'   used directly by \code{cv_bayes()}.
#' @param pointize Character; how to convert non-point geometries if present in
#'   \code{data_sf}. Passed to \code{\link{prep_model_data}()}. Common values:
#'   \code{"auto"} (default), \code{"surface"}, \code{"centroid"},
#'   \code{"line_midpoint"}, \code{"bbox_center"}.
#' @param fit_args Named list of additional arguments forwarded to
#'   \code{\link{fit_bayesian_spatial_model}()} (e.g., \code{family}, \code{gp_k},
#'   \code{gp_c}, \code{prior}, \code{backend}, \code{chains}, \code{iter}, etc.).
#' @param summary Character; how to summarize posterior predictions when available:
#'   \code{"mean"} (default) or \code{"median"}.
#'
#' @details
#' \strong{Preprocessing.} The function requires \code{\link{prep_model_data}()}
#' and \code{\link{fit_bayesian_spatial_model}()}. Preprocessing validates
#' modeling columns, coerces geometry to points if needed, ensures a projected
#' CRS (metric units), and drops rows with missing/non-finite values in the
#' modeling fields. It also preserves an \code{..row_id} column (original row
#' index) that this function uses to remap user-supplied folds to the kept rows.
#'
#' \strong{CRS.} Bayesian spatial models here use a Gaussian Process term over
#' coordinates. If data remain in geographic coordinates after preprocessing,
#' they are projected with \code{ensure_projected()} when available; otherwise an
#' error is raised. Coordinates used for prediction include both \code{x}/\code{y}
#' and \code{..x}/\code{..y} columns derived from geometry (the latter match the
#' defaults used by \code{\link{fit_bayesian_spatial_model}()} for \code{brms::gp()}).
#'
#' \strong{Prediction.} For \pkg{brms} fits, predictions use
#' \code{brms::posterior_epred()} and are summarized by \code{summary}; their
#' posterior standard deviation is returned as \code{yhat_sd}. For \pkg{rstanarm}
#' fits, \code{rstanarm::posterior_linpred(transform = TRUE)} is used. If neither
#' path is available, \code{predict(type = "response", se.fit = TRUE)} is tried,
#' falling back to plain \code{predict()}.
#'
#' \strong{Fold-level information criteria.} When the fitted object is a
#' \code{brmsfit}, the function attempts to compute LOOIC (\code{brms::loo()})
#' and WAIC (\code{brms::waic()}) for the training fit of each fold; failures are
#' tolerated and result in \code{NA}.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{overall}}{Data frame of aggregate metrics across all held-out
#'     predictions: \code{RMSE}, \code{MAE}, \code{R2}, and \code{n_pred}
#'     (number of finite predictions).}
#'   \item{\code{fold_metrics}}{Data frame with one row per evaluated fold,
#'     containing \code{fold}, \code{n_test}, \code{n_pred}, \code{RMSE},
#'     \code{MAE}, \code{R2}, \code{n_draws} (number of posterior draws used
#'     when applicable), and optionally \code{LOOIC}, \code{WAIC} (for
#'     \pkg{brms} models).}
#'   \item{\code{predictions}}{Data frame with columns \code{..row_id} (original
#'     row index), \code{fold}, observed \code{y}, predicted \code{yhat}, and
#'     \code{yhat_sd} (posterior SD or prediction SE when available).}
#'   \item{\code{folds}}{The fold definition actually used (list of
#'     \code{$train}/\code{$test} index vectors) after remapping to the kept rows.}
#'   \item{\code{formula}}{The bare modeling formula (without the spatial GP term),
#'     for reference/logging.}
#' }
#'
#' @section Failure handling:
#' Folds with too few training or testing observations are skipped. Model fitting
#' or prediction failures on a fold are reported via warnings (optionally through
#' \pkg{logger} if installed) and the fold contributes no metrics/predictions.
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' # Simulated projected point data
#' set.seed(1)
#' n <- 150
#' pts <- st_as_sf(
#'   data.frame(
#'     x = runif(n, 0, 2000),
#'     y = runif(n, 0, 2000),
#'     y_resp = rnorm(n),
#'     x1 = rnorm(n),
#'     x2 = rnorm(n)
#'   ),
#'   coords = c("x","y"), crs = 3857
#' )
#'
#' # 5-fold CV using defaults (requires brms)
#' res <- cv_bayes(
#'   data_sf = pts,
#'   response_var = "y_resp",
#'   predictor_vars = c("x1","x2"),
#'   k = 5, seed = 42,
#'   fit_args = list(iter = 1000, chains = 2)  # forwarded to fit_bayesian_spatial_model()
#' )
#' res$overall
#'
#' # With user-supplied folds
#' fdef <- make_folds(pts, k = 5, method = "random_kfold", seed = 42)
#' res2 <- cv_bayes(
#'   data_sf = pts,
#'   response_var = "y_resp",
#'   predictor_vars = c("x1","x2"),
#'   folds = fdef$folds
#' )
#' }
#'
#' @seealso
#' \code{\link{fit_bayesian_spatial_model}}, \code{\link{prep_model_data}},
#' \code{\link{make_folds}}, \code{\link{evaluate_models}}
#'
#' @importFrom sf st_geometry_type st_is_longlat st_coordinates st_drop_geometry
#' @importFrom stats as.formula predict na.omit
#' @export
cv_bayes <- function(
    data_sf,
    response_var,
    predictor_vars,
    folds     = NULL,
    k         = 5,
    seed      = 123,
    boundary  = NULL,
    pointize  = "auto",
    fit_args  = list(),
    summary   = c("mean", "median")
) {
  summary <- match.arg(summary)
  
  if (!inherits(data_sf, "sf")) stop("cv_bayes(): `data_sf` must be an sf object.")
  if (!exists("prep_model_data", mode = "function")) {
    stop("cv_bayes(): `prep_model_data()` not found.")
  }
  if (!exists("fit_bayesian_spatial_model", mode = "function")) {
    stop("cv_bayes(): `fit_bayesian_spatial_model()` not found.")
  }
  
  # Keep original ids for remapping; run centralized prep
  if (!("..row_id" %in% names(data_sf))) data_sf$`..row_id` <- seq_len(nrow(data_sf))
  
  dat_sf <- prep_model_data(
    data_sf        = data_sf,
    response_var   = response_var,
    predictor_vars = predictor_vars,
    boundary       = boundary,
    pointize       = pointize
  )
  
  if (!("..row_id" %in% names(dat_sf))) {
    stop("cv_bayes(): `prep_model_data()` must preserve a `..row_id` column.")
  }
  
  # Ensure POINT/MULTIPOINT and projected
  gtypes <- as.character(sf::st_geometry_type(dat_sf, by_geometry = TRUE))
  if (!all(gtypes %in% c("POINT", "MULTIPOINT"))) {
    stop("cv_bayes(): geometry must be POINT/MULTIPOINT after prep_model_data().")
  }
  is_ll <- tryCatch(sf::st_is_longlat(dat_sf), error = function(e) FALSE)
  if (isTRUE(is_ll)) {
    if (!exists("ensure_projected", mode = "function")) {
      stop("cv_bayes(): data are long/lat but `ensure_projected()` is unavailable.")
    }
    dat_sf <- ensure_projected(dat_sf)
  }
  
  # Coordinates (convenience; some models may use x/y columns)
  if (!all(c("x", "y") %in% names(dat_sf))) {
    crd <- sf::st_coordinates(dat_sf)
    dat_sf$x <- crd[, 1]
    dat_sf$y <- crd[, 2]
  }
  
  keep_idx <- dat_sf$`..row_id`
  
  # Bare formula for reference/logging
  fml <- if (length(predictor_vars) == 0L) {
    stats::as.formula(paste(response_var, "~ 1"))
  } else {
    stats::as.formula(paste(response_var, "~", paste(predictor_vars, collapse = " + ")))
  }
  
  # Folds: remap or create
  remapped_folds <- list()
  if (is.null(folds)) {
    set.seed(seed)
    ids <- keep_idx
    fold_id <- sample(rep(seq_len(k), length.out = length(ids)))
    for (i in seq_len(k)) {
      remapped_folds[[i]] <- list(train = ids[fold_id != i], test = ids[fold_id == i])
    }
  } else {
    for (i in seq_along(folds)) {
      tr_pos <- match(folds[[i]]$train, keep_idx)
      te_pos <- match(folds[[i]]$test,  keep_idx)
      remapped_folds[[i]] <- list(
        train = keep_idx[stats::na.omit(tr_pos)],
        test  = keep_idx[stats::na.omit(te_pos)]
      )
    }
  }
  
  .log_warn <- function(...) {
    msg <- sprintf(...)
    if (requireNamespace("logger", quietly = TRUE)) logger::log_warn(msg) else warning(msg)
  }
  
  pred_rows  <- list()
  fold_stats <- list()
  
  # ----------------------------- CV loop --------------------------------------
  for (i in seq_along(remapped_folds)) {
    tr_ids <- remapped_folds[[i]]$train
    te_ids <- remapped_folds[[i]]$test
    tr_pos <- stats::na.omit(match(tr_ids, keep_idx))
    te_pos <- stats::na.omit(match(te_ids, keep_idx))
    
    if (length(tr_pos) < 2L || length(te_pos) < 1L) {
      .log_warn("cv_bayes(): fold %d skipped (train n=%d, test n=%d).", i, length(tr_pos), length(te_pos))
      next
    }
    
    train_sf <- dat_sf[tr_pos, , drop = FALSE]
    test_sf  <- dat_sf[te_pos, , drop = FALSE]
    
    # ----------------------------- Fit ----------------------------------------
    fit <- try(
      do.call(
        fit_bayesian_spatial_model,
        c(list(
          data_sf        = train_sf,
          response_var   = response_var,
          predictor_vars = predictor_vars
        ), fit_args)
      ),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) {
      .log_warn("cv_bayes(): fold %d model fit failed; skipping.", i)
      next
    }
    model_obj <- if (!is.null(fit$model)) fit$model else fit
    
    # -------------------------- Predict (FIX A) -------------------------------
    # Ensure test data include the coordinates columns expected by spatial terms,
    # e.g., brms::gp(..x, ..y). We add ..x/..y from geometry and keep x/y too.
    test_df <- sf::st_drop_geometry(test_sf)
    crd_te  <- sf::st_coordinates(test_sf)
    test_df$`..x` <- crd_te[, 1]
    test_df$`..y` <- crd_te[, 2]
    if (!all(c("x", "y") %in% names(test_df))) {
      test_df$x <- crd_te[, 1]
      test_df$y <- crd_te[, 2]
    }
    
    y_true  <- test_df[[response_var]]
    y_hat   <- rep(NA_real_, length(y_true))
    y_sd    <- rep(NA_real_, length(y_true))
    n_draw  <- NA_integer_
    
    pred_ok <- FALSE
    
    # brms
    if (inherits(model_obj, "brmsfit") && requireNamespace("brms", quietly = TRUE)) {
      draws <- try(brms::posterior_epred(model_obj, newdata = test_df), silent = TRUE)
      if (!inherits(draws, "try-error") && is.matrix(draws)) {
        n_draw <- nrow(draws)
        y_hat  <- if (summary == "mean") colMeans(draws) else apply(draws, 2L, stats::median)
        y_sd   <- apply(draws, 2L, stats::sd)
        pred_ok <- TRUE
      }
    }
    
    # rstanarm
    if (!pred_ok && inherits(model_obj, "stanreg") && requireNamespace("rstanarm", quietly = TRUE)) {
      draws <- try(rstanarm::posterior_linpred(model_obj, newdata = test_df, transform = TRUE), silent = TRUE)
      if (!inherits(draws, "try-error") && is.matrix(draws)) {
        n_draw <- nrow(draws)
        y_hat  <- if (summary == "mean") colMeans(draws) else apply(draws, 2L, stats::median)
        y_sd   <- apply(draws, 2L, stats::sd)
        pred_ok <- TRUE
      }
    }
    
    # Generic predict()
    if (!pred_ok) {
      p1 <- try(stats::predict(model_obj, newdata = test_df, type = "response", se.fit = TRUE), silent = TRUE)
      if (!inherits(p1, "try-error")) {
        if (is.list(p1) && !is.null(p1$fit)) {
          y_hat <- as.numeric(p1$fit)
          if (!is.null(p1$se.fit)) y_sd <- as.numeric(p1$se.fit)
          pred_ok <- TRUE
        } else if (is.numeric(p1) && length(p1) == length(y_true)) {
          y_hat <- as.numeric(p1)
          pred_ok <- TRUE
        }
      }
    }
    if (!pred_ok) {
      p2 <- try(stats::predict(model_obj, newdata = test_df), silent = TRUE)
      if (!inherits(p2, "try-error") && is.numeric(p2) && length(p2) == length(y_true)) {
        y_hat <- as.numeric(p2)
        pred_ok <- TRUE
      }
    }
    if (!pred_ok) .log_warn("cv_bayes(): fold %d prediction failed; yhat set to NA.", i)
    
    # ------------------------- Per-fold metrics -------------------------------
    ok   <- is.finite(y_true) & is.finite(y_hat)
    n_ok <- sum(ok)
    if (n_ok == 0L) {
      .log_warn("cv_bayes(): fold %d produced no finite predictions.", i)
      next
    }
    resid <- y_true[ok] - y_hat[ok]
    rmse  <- sqrt(mean(resid^2))
    mae   <- mean(abs(resid))
    ssr   <- sum(resid^2)
    sst   <- sum((y_true[ok] - mean(y_true[ok]))^2)
    r2    <- if (sst > 0) 1 - ssr / sst else NA_real_
    
    # ---------------- Info criteria (brms) — FIX B ----------------------------
    looic <- NA_real_; waic <- NA_real_
    if (inherits(model_obj, "brmsfit") && requireNamespace("loo", quietly = TRUE)) {
      # LOO
      loo_res <- try(brms::loo(model_obj), silent = TRUE)
      if (!inherits(loo_res, "try-error") && !is.null(loo_res$estimates)) {
        est <- loo_res$estimates
        if ("looic" %in% rownames(est)) {
          looic <- as.numeric(est["looic", "Estimate"])
        } else if ("elpd_loo" %in% rownames(est)) {
          looic <- -2 * as.numeric(est["elpd_loo", "Estimate"])
        }
      }
      # WAIC
      waic_res <- try(brms::waic(model_obj), silent = TRUE)
      if (!inherits(waic_res, "try-error") && !is.null(waic_res$estimates)) {
        est <- waic_res$estimates
        if ("waic" %in% rownames(est)) {
          waic <- as.numeric(est["waic", "Estimate"])
        }
      }
    }
    
    fold_stats[[length(fold_stats) + 1L]] <- data.frame(
      fold    = i,
      n_test  = length(y_true),
      n_pred  = n_ok,
      RMSE    = rmse,
      MAE     = mae,
      R2      = r2,
      n_draws = n_draw,
      LOOIC   = looic,
      WAIC    = waic,
      stringsAsFactors = FALSE
    )
    
    pred_rows[[length(pred_rows) + 1L]] <- data.frame(
      `..row_id` = test_sf$`..row_id`,
      fold       = i,
      y          = as.numeric(y_true),
      yhat       = as.numeric(y_hat),
      yhat_sd    = ifelse(is.finite(y_sd), as.numeric(y_sd), NA_real_),
      stringsAsFactors = FALSE
    )
  }
  
  # ------------------------------ Bind results --------------------------------
  preds <- if (length(pred_rows)) do.call(rbind, pred_rows) else
    data.frame(`..row_id` = integer(), fold = integer(), y = numeric(), yhat = numeric(), yhat_sd = numeric())
  
  folds_df <- if (length(fold_stats)) do.call(rbind, fold_stats) else
    data.frame(fold = integer(), n_test = integer(), n_pred = integer(),
               RMSE = numeric(), MAE = numeric(), R2 = numeric(),
               n_draws = integer(), LOOIC = numeric(), WAIC = numeric())
  
  # Overall metrics
  ok_all <- is.finite(preds$y) & is.finite(preds$yhat)
  if (any(ok_all)) {
    resid_all <- preds$y[ok_all] - preds$yhat[ok_all]
    rmse_all  <- sqrt(mean(resid_all^2))
    mae_all   <- mean(abs(resid_all))
    ssr_all   <- sum(resid_all^2)
    sst_all   <- sum((preds$y[ok_all] - mean(preds$y[ok_all]))^2)
    r2_all    <- if (sst_all > 0) 1 - ssr_all / sst_all else NA_real_
    overall   <- data.frame(RMSE = rmse_all, MAE = mae_all, R2 = r2_all,
                            n_pred = sum(ok_all), stringsAsFactors = FALSE)
  } else {
    overall <- data.frame(RMSE = NA_real_, MAE = NA_real_, R2 = NA_real_,
                          n_pred = 0L, stringsAsFactors = FALSE)
  }
  
  list(
    overall      = overall,
    fold_metrics = folds_df,
    predictions  = preds,
    folds        = remapped_folds,
    formula      = deparse(fml)
  )
}

#' Build a tessellation (Voronoi, Delaunay triangles, hex grid, or square grid)
#'
#' Creates polygonal tessellations from point data using one of several methods:
#' Voronoi polygons, Delaunay triangles, hexagonal grids, or square grids. The
#' function normalizes coordinate reference systems (CRS), optionally clips to a
#' boundary, and returns a structured list suitable for downstream analysis.
#'
#' @param points_sf An \pkg{sf} object with \code{POINT}/\code{MULTIPOINT}
#'   geometry. Required for all methods; used as the seed set for Voronoi and
#'   triangles, and only to infer/align CRS for grid methods.
#' @param boundary Optional \pkg{sf}/\pkg{sfc} polygon or multipolygon to limit
#'   the tessellation. Required for grid methods (\code{hex}, \code{square}).
#'   For \code{voronoi}, if omitted a convex hull of the points is used as the
#'   envelope; for \code{triangles}, clipping to \code{boundary} is optional
#'   (controlled by \code{clip}).
#' @param method Character; one of \code{"voronoi"}, \code{"triangles"},
#'   \code{"hex"}, \code{"square"}.
#' @param approx_n_cells Approximate number of cells to target (used only when
#'   \code{method %in% c("hex","square")} and \code{cellsize} is not supplied).
#'   The actual number will depend on the boundary extent and clipping.
#' @param cellsize Numeric length 1 or 2 specifying grid cell size in CRS units
#'   (x[, y]) for grid methods. If length 1, the same value is used for both
#'   axes. Ignored for \code{voronoi} and \code{triangles}.
#' @param expand Numeric; buffer distance in CRS units applied to the Voronoi
#'   envelope before constructing cells (only used for \code{method = "voronoi"}).
#'   Set to \code{0} for no expansion. Ignored for \code{triangles} and grid
#'   methods.
#' @param clip Logical; if \code{TRUE}, clip output polygons to \code{boundary}
#'   (when supplied). For \code{voronoi}, clipping happens to the (possibly
#'   expanded) boundary; for \code{triangles}, triangles are intersected with the
#'   boundary; for grids, polygons are generated inside the boundary so clipping
#'   is effectively always \code{TRUE}.
#' @param keep_duplicates Logical; when \code{FALSE} (default), duplicate point
#'   coordinates are removed before constructing the Voronoi graph or triangle
#'   mesh. For Voronoi, the returned index maps all original (including duplicate)
#'   points back to cells deterministically.
#' @param crs Optional target CRS (anything accepted by \code{sf::st_crs()}). If
#'   supplied, \code{points_sf} and \code{boundary} are transformed to this CRS
#'   before tessellation. If omitted and \code{points_sf} is in geographic
#'   coordinates, a local projected CRS is chosen via \code{ensure_projected()}
#'   (if available); otherwise \code{boundary} is aligned to the points' CRS.
#' @param quiet Logical; suppress informative messages when \code{TRUE}.
#'
#' @details
#' \strong{CRS handling:}
#' \itemize{
#'   \item If \code{crs} is provided, both inputs are transformed to it.
#'   \item If \code{points_sf} is long/lat and \code{crs} is \code{NULL}, the
#'         function projects points to a local projected CRS (via
#'         \code{ensure_projected()} if available). The \code{boundary} is then
#'         transformed to match.
#'   \item Validity of polygonal \code{boundary} is enforced with
#'         \code{sf::st_make_valid()} (or a safe fallback).
#' }
#'
#' \strong{Methods:}
#' \describe{
#'   \item{\code{voronoi}}{Builds Voronoi cells using \code{sf::st_voronoi()}
#'     with an envelope derived from \code{boundary} (buffered by
#'     \code{expand} if \code{expand > 0}) or, if \code{boundary} is missing,
#'     from the convex hull of \code{points_sf}. Cells are optionally clipped to
#'     the envelope. A point-to-cell index is returned.}
#'   \item{\code{triangles}}{Builds a Delaunay triangulation over the unique
#'     points using \pkg{geometry} (\code{geometry::delaunayn}) when available,
#'     otherwise falls back to \code{sf::st_triangulate()} on the convex hull.
#'     Triangles may be clipped to \code{boundary} if \code{clip = TRUE}.}
#'   \item{\code{hex}, \code{square}}{Creates a hexagonal or square grid within
#'     \code{boundary} using \code{\link{create_grid_polygons}()} with either
#'     \code{cellsize} or \code{approx_n_cells}. \code{boundary} is required.}
#' }
#'
#' \strong{Dependencies & fallbacks:}
#' The function requires \pkg{sf}. The \pkg{geometry} package is optional but
#' preferred for \code{triangles}; if unavailable or failing, an \pkg{sf}
#' triangulation fallback is used. Polygon validity is ensured via
#' \code{sf::st_make_valid()} when exported (or a \code{buffer(., 0)} fallback).
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{cells}}{An \pkg{sf} polygon layer of tessellation cells.}
#'   \item{\code{index}}{For \code{method = "voronoi"}, an integer vector of
#'     length \code{nrow(points_sf)} mapping each input point to its cell row
#'     index (resolving ties deterministically). \code{NULL} for other methods.}
#'   \item{\code{boundary}}{The polygon boundary used (supplied or derived).}
#'   \item{\code{method}}{The method used: \code{"voronoi"}, \code{"triangles"},
#'     \code{"hex"}, or \code{"square"}.}
#'   \item{\code{params}}{A list echoing key parameters used (e.g.,
#'     \code{clip}, \code{expand}, \code{keep_duplicates}, \code{cellsize},
#'     \code{approx_n_cells}).}
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' # Sample points (projected)
#' set.seed(1)
#' pts <- st_as_sf(
#'   data.frame(x = runif(200, 0, 2000), y = runif(200, 0, 2000)),
#'   coords = c("x","y"), crs = 3857
#' )
#' bnd <- st_as_sfc(st_bbox(pts)) |> st_as_sf()
#'
#' # Voronoi (clipped to boundary, slight expansion of envelope)
#' v <- build_tessellation(pts, boundary = bnd, method = "voronoi", expand = 50)
#' plot(st_geometry(v$cells))
#'
#' # Delaunay triangles (clipped)
#' t <- build_tessellation(pts, boundary = bnd, method = "triangles", clip = TRUE)
#' plot(st_geometry(t$cells))
#'
#' # Hex grid targeting ~150 cells
#' h <- build_tessellation(pts, boundary = bnd, method = "hex", approx_n_cells = 150)
#' plot(st_geometry(h$cells))
#'
#' # Square grid with fixed cell size (CRS units)
#' s <- build_tessellation(pts, boundary = bnd, method = "square", cellsize = 250)
#' plot(st_geometry(s$cells))
#' }
#'
#' @seealso
#' \code{\link{create_voronoi_polygons}}, \code{\link{create_grid_polygons}},
#' \code{\link{clip_target_for}}
#'
#' @importFrom sf st_geometry_type st_crs st_transform st_is_longlat st_coordinates
#'   st_convex_hull st_union st_triangulate st_collection_extract st_sfc st_buffer
#'   st_intersection st_is_empty st_sf
#' @export
build_tessellation <- function(
    points_sf,
    boundary        = NULL,
    method          = c("voronoi", "triangles", "hex", "square"),
    approx_n_cells  = NULL,
    cellsize        = NULL,
    expand          = 0,
    clip            = TRUE,
    keep_duplicates = FALSE,
    crs             = NULL,
    quiet           = FALSE
) {
  # ---- debug helper ----
  .msg <- function(...) if (!quiet) message(...)
  method <- match.arg(method)
  
  # ---- safety & helpers ----
  if (!inherits(points_sf, "sf")) stop("build_tessellation(): `points_sf` must be an sf object.")
  gtypes <- as.character(sf::st_geometry_type(points_sf, by_geometry = TRUE))
  if (!all(gtypes %in% c("POINT", "MULTIPOINT"))) {
    stop("build_tessellation(): `points_sf` geometry must be POINT/MULTIPOINT.")
  }
  
  .safe_make_valid <- function(g) {
    # Prefer sf::st_make_valid if available
    if ("st_make_valid" %in% getNamespaceExports("sf")) {
      return(suppressWarnings(sf::st_make_valid(g)))
    }
    # Try lwgeom only if the symbol actually exists
    if (requireNamespace("lwgeom", quietly = TRUE) &&
        "st_make_valid" %in% getNamespaceExports("lwgeom")) {
      return(suppressWarnings(lwgeom::st_make_valid(g)))
    }
    # Fallback
    suppressWarnings(sf::st_buffer(g, 0))
  }
  
  .dedupe_points <- function(g) {
    m <- sf::st_coordinates(g)
    if (nrow(m) == 0) return(g[FALSE, , drop = FALSE])
    key <- paste0(round(m[, 1], 10), "_", round(m[, 2], 10))
    keep <- !duplicated(key)
    g[keep, , drop = FALSE]
  }
  
  .align_crs <- function(a, b) {
    if (is.null(b)) return(a)
    if (is.na(sf::st_crs(a)) || is.na(sf::st_crs(b))) return(a)
    if (sf::st_crs(a) == sf::st_crs(b)) return(a)
    sf::st_transform(a, sf::st_crs(b))
  }
  
  .is_longlat <- function(x) {
    out <- tryCatch(sf::st_is_longlat(x), error = function(e) NA)
    isTRUE(out)
  }
  
  # ---- CRS handling ----
  if (!is.null(crs)) {
    points_sf <- sf::st_transform(points_sf, crs)
    if (!is.null(boundary)) boundary <- sf::st_transform(boundary, crs)
  } else {
    if (.is_longlat(points_sf)) {
      .msg("build_tessellation(): projecting points to a local UTM CRS.")
      points_sf <- ensure_projected(points_sf)
    }
    if (!is.null(boundary)) boundary <- .align_crs(boundary, points_sf)
  }
  
  if (!is.null(boundary)) {
    if (!any(sf::st_geometry_type(boundary) %in% c("POLYGON", "MULTIPOLYGON"))) {
      stop("build_tessellation(): `boundary` must be polygonal.")
    }
    boundary <- .safe_make_valid(boundary)
  }
  
  # ---- Branch by method ----
  if (identical(method, "voronoi")) {
    res <- create_voronoi_polygons(
      points_sf       = points_sf,
      boundary        = boundary,
      expand          = expand,
      clip            = clip,
      keep_duplicates = keep_duplicates,
      crs             = sf::st_crs(points_sf),
      quiet           = quiet
    )
    # normalize the return shape (already a list from create_voronoi_polygons)
    return(res)
  }
  
  if (identical(method, "hex") || identical(method, "square")) {
    if (is.null(boundary)) {
      stop("build_tessellation(): `boundary` is required for hex/square grids.")
    }
    type <- if (identical(method, "hex")) "hex" else "square"
    # Delegate to the grid builder (it already logs sizing)
    grid <- create_grid_polygons(
      boundary      = boundary,
      target_cells  = approx_n_cells,
      type          = type,
      cellsize      = cellsize,
      quiet         = quiet
    )
    return(list(
      cells    = grid,
      index    = NULL,
      boundary = boundary,
      method   = method,
      params   = list(
        approx_n_cells  = approx_n_cells,
        cellsize        = cellsize,
        clip            = TRUE,   # grids are generated inside boundary already
        keep_duplicates = keep_duplicates,
        expand          = 0
      )
    ))
  }
  
  if (identical(method, "triangles")) {
    # Triangles: prefer Delaunay on the input points via geometry::delaunayn
    pts <- if (isTRUE(keep_duplicates)) points_sf else .dedupe_points(points_sf)
    if (nrow(pts) < 3L) stop("build_tessellation(triangles): need at least 3 unique points.")
    
    coords <- sf::st_coordinates(pts)[, 1:2, drop = FALSE]
    
    tri_sfc <- NULL
    used <- "geometry::delaunayn"
    if (requireNamespace("geometry", quietly = TRUE)) {
      tri_idx <- try(geometry::delaunayn(coords), silent = TRUE)
      if (!inherits(tri_idx, "try-error") && length(tri_idx)) {
        # Build polygons from simplices
        polys <- vector("list", nrow(tri_idx))
        for (i in seq_len(nrow(tri_idx))) {
          idx <- tri_idx[i, ]
          ring <- rbind(coords[idx, , drop = FALSE], coords[idx[1], , drop = FALSE])
          polys[[i]] <- sf::st_polygon(list(ring))
        }
        tri_sfc <- sf::st_sfc(polys, crs = sf::st_crs(pts))
      } else {
        used <- "sf::st_triangulate fallback"
      }
    } else {
      used <- "sf::st_triangulate fallback"
    }
    
    if (is.null(tri_sfc)) {
      # Fallback: triangulate the convex hull surface
      hull <- sf::st_convex_hull(sf::st_union(sf::st_geometry(pts)))
      tri_sfc <- sf::st_triangulate(hull)
      tri_sfc <- sf::st_collection_extract(tri_sfc, "POLYGON", warn = FALSE)
      tri_sfc <- sf::st_sfc(tri_sfc, crs = sf::st_crs(pts))
    }
    
    .msg(sprintf("build_tessellation(triangles): built %d triangles via %s.",
                 length(tri_sfc), used))
    
    tri_sf <- sf::st_sf(geometry = .safe_make_valid(tri_sfc))
    if (!is.null(boundary) && isTRUE(clip)) {
      tri_sf <- suppressWarnings(sf::st_intersection(tri_sf, boundary))
      empt <- sf::st_is_empty(tri_sf)
      if (any(empt, na.rm = TRUE)) tri_sf <- tri_sf[!empt, , drop = FALSE]
      .msg(sprintf("build_tessellation(triangles): clipped triangles = %d.", nrow(tri_sf)))
    }
    
    tri_sf$cell_id <- seq_len(nrow(tri_sf))
    
    return(list(
      cells    = tri_sf,
      index    = NULL,
      boundary = boundary,
      method   = "triangles",
      params   = list(
        clip            = clip,
        approx_n_cells  = approx_n_cells,
        keep_duplicates = keep_duplicates,
        expand          = 0
      )
    ))
  }
  
  stop("build_tessellation(): unknown method.")
}

#' Create deterministic, stable polygon IDs based on spatial sort keys
#'
#' Ensures that a polygon layer has a reproducible, deterministic identifier
#' column by sorting features using representative point coordinates (and
#' secondary tie-breakers) and then assigning sequential IDs. This is useful for
#' keeping cell/feature IDs stable across runs when geometry input order may
#' vary (e.g., after unions, intersections, or grid generation).
#'
#' @param polygons_sf An \pkg{sf} or \pkg{sfc} object containing only
#'   polygonal features (\code{POLYGON} or \code{MULTIPOLYGON}). If an \code{sfc}
#'   is supplied it is coerced to \code{sf}. Non-polygon rows are dropped
#'   silently; if no polygon rows remain the function errors.
#' @param id_col Character scalar. Name of the identifier column to write on the
#'   returned object. If a column with this name already exists, it will be
#'   overwritten.
#' @param method One of \code{"centroid"}, \code{"surface_point"}, or
#'   \code{"bbox_center"} specifying how to compute the representative point
#'   used to derive the sorting key:
#'   \itemize{
#'     \item \code{"centroid"} — geometric centroid (may fall outside concave polygons).
#'     \item \code{"surface_point"} — a point guaranteed to lie on/inside each polygon
#'           (\code{sf::st_point_on_surface()}).
#'     \item \code{"bbox_center"} — numeric center of each feature's bounding box
#'           (fastest; purely geometric, not spherical).
#'   }
#' @param make_valid Logical; if \code{TRUE} (default), applies
#'   \code{sf::st_make_valid()} to input polygons before computing sort keys.
#' @param transform_for_sort A CRS understood by \code{sf::st_crs()} (e.g., an
#'   EPSG code like \code{4326}, a PROJ string, or an \code{sf} CRS object)
#'   used \emph{only} to compute the sort key coordinates. If \code{polygons_sf}
#'   has an undefined CRS or the transform fails, the existing CRS is used
#'   instead. Set to \code{NULL} to disable any transform for sorting.
#'
#' @details
#' The function (i) optionally validates polygons, (ii) optionally transforms a
#' copy to \code{transform_for_sort} to make coordinates comparable, (iii)
#' computes a representative point per polygon, (iv) builds a total order using
#' the representative point's \code{x} then \code{y} coordinates, with ties
#' broken by polygon \emph{area} (ascending) and finally the original row index,
#' and (v) assigns sequential IDs (\code{1..n}) into \code{id_col}. Multi-part
#' features (\code{MULTIPOLYGON}) are \emph{not} exploded; IDs are assigned per
#' input row.
#'
#' The output is re-ordered by the computed sort key. All original attributes
#' are preserved; only the row order changes and \code{id_col} is (over)written.
#'
#' @return An \pkg{sf} polygon layer containing the same polygon rows (possibly
#'   re-ordered) and an added/overwritten identifier column \code{id_col} with
#'   values \code{1..n}.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' # Example polygons (unordered)
#' p <- st_as_sf(
#'   data.frame(id = 101:103),
#'   wkt = c("POLYGON((0 0,2 0,2 1,0 1,0 0))",
#'           "POLYGON((3 0,4 0,4 2,3 2,3 0))",
#'           "POLYGON((0 2,1 2,1 3,0 3,0 2))"),
#'   crs = 4326
#' )
#'
#' # Stable IDs by centroid in EPSG:4326
#' p1 <- ensure_stable_poly_id(p, id_col = "poly_id", method = "centroid", transform_for_sort = 4326)
#'
#' # Stable IDs with inside point and without transforming CRS for sorting
#' p2 <- ensure_stable_poly_id(p, id_col = "cell_id", method = "surface_point", transform_for_sort = NULL)
#' }
#'
#' @seealso \code{\link[sf]{st_point_on_surface}}, \code{\link[sf]{st_centroid}},
#'   \code{\link[sf]{st_make_valid}}, \code{\link[sf]{st_transform}}
#'
#' @importFrom sf st_as_sf st_geometry_type st_make_valid st_transform
#'   st_point_on_surface st_centroid st_bbox st_point st_sfc st_crs
#'   st_geometry st_coordinates st_area st_sf
#' @export
ensure_stable_poly_id <- function(polygons_sf,
                                  id_col = "poly_id",
                                  method = c("centroid", "surface_point", "bbox_center"),
                                  make_valid = TRUE,
                                  transform_for_sort = 4326) {
  # Normalize to sf
  if (inherits(polygons_sf, "sfc")) polygons_sf <- sf::st_as_sf(polygons_sf)
  if (!inherits(polygons_sf, "sf")) {
    stop("ensure_stable_poly_id(): `polygons_sf` must be an sf/sfc object.")
  }
  # Keep only POLYGON/MULTIPOLYGON rows; DO NOT explode MULTIPOLYGONs
  gtypes <- as.character(sf::st_geometry_type(polygons_sf, by_geometry = TRUE))
  keep   <- gtypes %in% c("POLYGON", "MULTIPOLYGON")
  polygons_sf <- polygons_sf[keep, , drop = FALSE]
  if (nrow(polygons_sf) == 0L) stop("ensure_stable_poly_id(): no polygon rows found.")
  
  if (isTRUE(make_valid)) {
    polygons_sf <- sf::st_make_valid(polygons_sf)
  }
  
  method <- match.arg(method)
  
  # Work on a transformed copy for the sort key (if transform makes sense)
  sort_sf <- polygons_sf
  if (!is.null(transform_for_sort) && !is.na(sf::st_crs(sort_sf))) {
    sort_sf <- tryCatch(sf::st_transform(sort_sf, transform_for_sort), error = function(e) sort_sf)
  }
  
  # Representative points
  rep_pts <- switch(
    method,
    centroid      = sf::st_centroid(sort_sf),
    surface_point = sf::st_point_on_surface(sort_sf),
    bbox_center   = {
      geoms <- sf::st_geometry(sort_sf)
      centers <- lapply(geoms, function(g) {
        bb <- sf::st_bbox(g)
        sf::st_point(c((bb["xmin"] + bb["xmax"]) / 2, (bb["ymin"] + bb["ymax"]) / 2))
      })
      sf::st_sf(geometry = sf::st_sfc(centers, crs = sf::st_crs(sort_sf)))
    }
  )
  
  # Coordinates + tie-breakers
  xy <- suppressWarnings(sf::st_coordinates(sf::st_geometry(rep_pts)))
  if (!is.matrix(xy) || nrow(xy) != nrow(polygons_sf)) {
    stop("ensure_stable_poly_id(): failed to compute representative coordinates.")
  }
  area <- suppressWarnings(as.numeric(sf::st_area(sort_sf)))
  area[!is.finite(area)] <- 0
  idx0 <- seq_len(nrow(polygons_sf))
  
  # Build total order: x, y, area (asc), original row
  ord <- do.call(order, list(xy[, 1], xy[, 2], area, idx0))
  
  out <- polygons_sf[ord, , drop = FALSE]
  out[[id_col]] <- seq_len(nrow(out))
  out
}

#' Internal package cache environment
#'
#' A singleton \code{environment} used to store ephemeral, in-session caches
#' (e.g., memoized results or temporary objects) shared across functions in this
#' package. The environment is created with \code{emptyenv()} as its parent to
#' avoid accidental name lookups/inheritance from the search path.
#'
#' @details
#' This object is initialized at load time as:
#' \code{.gmt_cache <- new.env(parent = emptyenv())}.
#' It is not exported and is intended for internal use only. If you need to
#' clear it, you can do so with:
#' \preformatted{
#'   rm(list = ls(envir = .gmt_cache, all.names = TRUE), envir = .gmt_cache)
#' }
#'
#' @format An environment.
#'
#' @name .gmt_cache
#' @keywords internal
#' @noRd
.gmt_cache <- new.env(parent = emptyenv())

#' Build a deterministic cache key from geometry, CRS, and parameters
#'
#' Constructs a stable, hash-based key that uniquely identifies a tessellation/
#' grid request given a boundary geometry, the requested grid/tile \code{type},
#' an approximate target cell count, and any additional parameters supplied in
#' \code{...}. The key is deterministic across sessions and insensitive to
#' multipart polygon ordering and the order of named arguments in \code{...}.
#'
#' @details
#' The key is composed of three parts joined by \code{"::"}:
#' \enumerate{
#'   \item \strong{type}: included verbatim.
#'   \item \strong{target cell count}: coerced with \code{as.integer()} (may be
#'         \code{NA} if not provided).
#'   \item \strong{hash}: \code{digest::digest()} of a list containing:
#'     \itemize{
#'       \item \code{geom_wkt}: \code{sf::st_as_text(sf::st_union(boundary))}.
#'             If the union fails, falls back to \code{st_as_text(boundary)}.
#'             Using the union reduces sensitivity to the ordering of multipart
#'             geometries.
#'       \item \code{crs}: a token derived from \code{sf::st_crs(boundary)}:
#'             the \code{$input} string if available, otherwise the EPSG code
#'             (\code{$epsg}); if neither is available, the literal
#'             \code{"NA_CRS"}.
#'       \item \code{args}: the \code{...} arguments, normalized by sorting
#'             lexicographically by (possibly missing) names. Unnamed arguments
#'             are treated as empty-name entries so that different call orders
#'             produce the same key.
#'     }
#' }
#'
#' This helper is intended for internal caching (e.g., using \code{.gmt_cache})
#' and does not modify inputs. It requires \pkg{sf} and \pkg{digest}.
#'
#' @param boundary An \code{sf} or \code{sfc} object (typically
#'   \code{POLYGON}/\code{MULTIPOLYGON}) used to derive a geometry fingerprint
#'   and CRS token. Only geometry and CRS are read; attributes are ignored.
#' @param type Character scalar describing the tessellation/grid type (e.g.,
#'   \code{"hex"}, \code{"square"}, \code{"voronoi"}, etc.). Included verbatim
#'   in the key prefix.
#' @param target_cells Optional numeric/integer giving an approximate number of
#'   desired cells. Coerced with \code{as.integer()} and included in the key
#'   (may be \code{NA}).
#' @param ... Additional parameters that materially affect the tessellation or
#'   grid (e.g., \code{cellsize}, \code{expand}). Values are captured and
#'   included in the hash after name-based normalization. The function does not
#'   interpret these arguments; it only fingerprints them.
#'
#' @return A length-1 character vector: \code{"<type>::<target>::<hash>"}.
#'
#' @section Notes:
#' Different CRSs for the same coordinates will yield different keys because the
#' CRS token is included. If \code{boundary} has no CRS, the token becomes
#' \code{"NA_CRS"} and keys will match only when the geometry and parameters are
#' otherwise identical.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' b <- st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 1000, ymax = 800)), crs = 3857)
#' .cache_key(b, type = "hex", target_cells = 250, cellsize = 100)
#' }
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom sf st_crs st_union st_as_text
#' @importFrom digest digest
.cache_key <- function(boundary, type, target_cells, ...) {
  # Normalize and capture CRS token
  crs_obj <- sf::st_crs(boundary)
  crs_token <- if (!is.null(crs_obj) && !is.na(crs_obj)) {
    inp <- crs_obj$input
    eps <- crs_obj$epsg
    if (!is.null(inp) && !is.na(inp) && nzchar(as.character(inp))) {
      as.character(inp)
    } else if (!is.null(eps) && !is.na(eps)) {
      as.character(eps)
    } else {
      "NA_CRS"
    }
  } else {
    "NA_CRS"
  }
  
  # Hash geometry + params; union to reduce multipart ordering effects
  geom_wkt <- tryCatch(
    sf::st_as_text(sf::st_union(boundary)),
    error = function(e) sf::st_as_text(boundary)
  )
  
  # Capture dots in a stable way (sort by (possibly missing) names using base R)
  dots <- list(...)
  if (length(dots)) {
    nms <- names(dots)
    if (!is.null(nms)) {
      nms_ord <- nms
      nms_ord[is.na(nms_ord) | !nzchar(nms_ord)] <- ""
      ord <- order(nms_ord, na.last = TRUE)
      dots <- dots[ord]
    }
  }
  
  paste0(
    type, "::", as.integer(target_cells), "::",
    digest::digest(list(geom_wkt = geom_wkt, crs = crs_token, args = dots))
  )
}

#' Create and cache grid polygons over a boundary
#'
#' Builds a square or hexagonal grid clipped to a polygonal \code{boundary}
#' (via [create_grid_polygons()]) and memoizes the result in an environment so
#' repeated calls with the same inputs return instantly. The cache key is
#' deterministic and derived from the boundary geometry/CRS and all
#' grid-affecting parameters (via [\.cache_key()]). Returned polygons are
#' assigned a stable identifier column using [ensure_stable_poly_id()].
#'
#' @details
#' This wrapper:
#' \enumerate{
#'   \item Normalizes \code{boundary} to an \code{sf} object and projects it to
#'         a local metric CRS with [ensure_projected()] so cell sizing is in
#'         linear units.
#'   \item Computes a cache key using [\.cache_key()] that fingerprints the
#'         (projected) geometry, its CRS token, the grid \code{type}, the
#'         \code{target_cells} value, and every argument passed in \code{...}.
#'   \item On a cache hit, returns the cached \code{sf} grid immediately.
#'   \item On a miss, calls [create_grid_polygons()] with the supplied arguments,
#'         normalizes IDs with [ensure_stable_poly_id()], stores the result in
#'         \code{cache_env}, and returns it.
#' }
#'
#' The cache is an R environment (default: [\.gmt_cache]) and persists for the
#' life of the R session. The cache key includes a CRS token; changing the CRS
#' (or the boundary geometry) produces a different key.
#'
#' @param boundary An \code{sf} or \code{sfc} polygonal object
#'   (\code{POLYGON}/\code{MULTIPOLYGON}) defining the region of interest. It is
#'   coerced to \code{sf} if needed and projected to a local metric CRS for grid
#'   construction.
#' @param target_cells Approximate desired number of cells (integer/numeric).
#'   Used by [create_grid_polygons()] when \code{cellsize} and \code{n} are not
#'   supplied.
#' @param type Grid type to generate; one of \code{"hex"} or \code{"square"}.
#' @param ... Additional arguments forwarded to [create_grid_polygons()], such as
#'   \code{cellsize}, \code{n}, \code{clip}, \code{crs}, and \code{quiet}.
#' @param cache_env Environment used to store memoized grids. Defaults to
#'   [\.gmt_cache].
#'
#' @return An \code{sf} data frame of grid polygons (clipped to the boundary if
#'   \code{clip = TRUE}), with a stable identifier column (default:
#'   \code{poly_id}) created by [ensure_stable_poly_id()].
#'
#' @examples
#' \dontrun{
#' library(sf)
#' # Simple boundary (projected)
#' b <- st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 2000, ymax = 1200)), crs = 3857)
#'
#' # First call builds and caches
#' g1 <- create_grid_polygons_cached(b, target_cells = 250, type = "hex", clip = TRUE)
#'
#' # Second call (same args) returns instantly from cache
#' g2 <- create_grid_polygons_cached(b, target_cells = 250, type = "hex", clip = TRUE)
#' identical(st_geometry(g1), st_geometry(g2))
#' }
#'
#' @seealso [create_grid_polygons()], [ensure_stable_poly_id()], [\.cache_key()],
#'   [ensure_projected()]
#'
#' @export
#' @importFrom sf st_as_sf
create_grid_polygons_cached <- function(boundary,
                                        target_cells,
                                        type = c("hex","square"),
                                        ...,
                                        cache_env = .gmt_cache) {
  type <- match.arg(type)
  
  # Normalize boundary to sf and ensure projected CRS for robust grid sizing
  bnd <- if (inherits(boundary, "sfc")) sf::st_as_sf(boundary) else boundary
  if (!inherits(bnd, "sf")) stop("create_grid_polygons_cached(): 'boundary' must be sf/sfc POLYGON/MULTIPOLYGON.")
  bnd <- ensure_projected(bnd)
  
  # Compose cache key (include all grid-affecting extras)
  key <- .cache_key(bnd, type, target_cells, ...)
  
  # Cache hit
  if (exists(key, envir = cache_env, inherits = FALSE)) {
    return(get(key, envir = cache_env, inherits = FALSE))
  }
  
  # Build, normalize IDs, then cache
  out <- create_grid_polygons(bnd, target_cells = target_cells, type = type, ...)
  out <- ensure_stable_poly_id(out)
  
  assign(key, out, envir = cache_env)
  out
}

#' Cross-validated comparison of spatial models with optional tessellation output
#'
#' Runs k-fold (or user-supplied) cross-validation for one or both of:
#' \itemize{
#'   \item Geographically Weighted Regression (GWR) via \code{\link{cv_gwr}}
#'   \item A Bayesian spatial model via \code{\link{cv_bayes}}
#' }
#' and returns fold-level metrics, overall metrics, and (optionally) a
#' lightweight tessellation object that can be used for mapping/diagnostics.
#' The tessellation is \emph{not} used in the CV itself; it is provided to help
#' visualize the study region and any gridding you asked for.
#'
#' @details
#' This function is a convenience orchestrator:
#' \enumerate{
#'   \item Validates inputs and ensures an internal row id (\code{..row_id}).
#'   \item Optionally builds a tessellation over a projected boundary:
#'         \itemize{
#'           \item For \code{tess_method \%in\% c("grid","hex","square")}, calls
#'                 \code{\link{create_grid_polygons}} (default
#'                 \code{target_cells = 30}) and passes any extra arguments from
#'                 \code{tess_args}.
#'           \item For \code{tess_method \%in\% c("voronoi","triangles")}, calls
#'                 \code{\link{build_tessellation}} with the input points (coerced
#'                 to POINTs if needed) and any \code{tess_args}; the returned
#'                 \code{$cells} are included in the output.
#'         }
#'   \item Runs \code{\link{cv_gwr}} and/or \code{\link{cv_bayes}} depending on
#'         \code{models}, forwarding only arguments accepted by those functions.
#' }
#'
#' The \code{boundary} (if supplied) is projected to a metric CRS via
#' \code{\link{ensure_projected}} for stable geometry operations. If no
#' \code{boundary} is provided, a slightly expanded convex hull of the (projected)
#' input points is used for tessellation.
#'
#' @param data_sf An \code{sf} object containing the response, predictors, and
#'   geometries. Non-point geometries are allowed; they will be coerced to
#'   points according to \code{pointize} when building the tessellation and
#'   inside the CV helpers.
#' @param response_var Character scalar; name of the response variable column.
#' @param predictor_vars Character vector; names of predictor columns to include.
#' @param k Integer; number of folds when \code{folds} is \code{NULL}. Default \code{5}.
#' @param seed Optional integer seed for reproducibility of random fold
#'   assignment (when \code{folds} is \code{NULL}) and any stochastic components.
#' @param folds Optional list of precomputed fold splits of the form
#'   \code{list(list(train = int, test = int), ...)} using original row indices.
#'   If provided, \code{k} and \code{seed} are ignored for fold creation.
#' @param boundary Optional polygon \code{sf}/\code{sfc} used to derive the
#'   tessellation extent. If \code{NULL}, a buffered convex hull of the input
#'   points (in a projected CRS) is used.
#' @param pointize Strategy used to coerce non-point geometries to points when
#'   needed; passed to \code{\link{coerce_to_points}}. One of
#'   \code{"auto"}, \code{"surface"}, \code{"centroid"}, \code{"line_midpoint"},
#'   or \code{"bbox_center"}. Default \code{"auto"}.
#' @param tess_method Tessellation type to create for diagnostics. One of
#'   \code{"grid"}, \code{"hex"}, \code{"square"}, \code{"voronoi"}, or
#'   \code{"triangles"}. For \code{"grid"} the grid is square; for
#'   \code{"hex"} and \code{"square"} the corresponding grid is generated.
#' @param tess_args Named list of additional arguments forwarded to the
#'   tessellation builders (\code{\link{create_grid_polygons}} or
#'   \code{\link{build_tessellation}}), e.g. \code{target_cells}, \code{cellsize},
#'   \code{expand}, \code{clip}, etc.
#' @param summary Character; how to summarize posterior predictions in the
#'   Bayesian CV path. One of \code{"mean"} or \code{"median"} (forwarded to
#'   \code{\link{cv_bayes}}). Default \code{"mean"}.
#' @param models Character vector selecting which models to run. Any subset of
#'   \code{c("GWR","Bayesian")}. Models are skipped automatically if their
#'   required packages/functions are unavailable.
#' @param quiet Logical; if \code{TRUE}, suppresses informational messages.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{\code{overall}}{A data frame with one row per fitted model containing
#'         RMSE, MAE, R2, and \code{n_pred} aggregated over folds.}
#'   \item{\code{by_fold}}{A data frame of fold-level metrics for each model
#'         (e.g., RMSE, MAE, R2, bandwidth/AICc for GWR; LOOIC/WAIC for Bayesian
#'         when available).}
#'   \item{\code{tessellation}}{A list describing the diagnostic tessellation:
#'         \code{list(method = <character>, args = <list>, cells = <sf POLYGON>)}.}
#' }
#'
#' @section Dependencies:
#' \itemize{
#'   \item GWR path requires package \pkg{spgwr} and function \code{\link{cv_gwr}}.
#'   \item Bayesian path requires package \pkg{brms} and function \code{\link{cv_bayes}}.
#'   \item Tessellation helpers use \pkg{sf} and functions
#'         \code{\link{create_grid_polygons}} and/or \code{\link{build_tessellation}}.
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' # Example points (projected)
#' pts <- st_as_sf(data.frame(x = runif(300, 0, 10000),
#'                            y = runif(300, 0, 8000),
#'                            yvar = rnorm(300),
#'                            x1 = rnorm(300), x2 = rnorm(300)),
#'                 coords = c("x","y"), crs = 3857)
#'
#' res <- evaluate_models_cv(
#'   data_sf       = pts,
#'   response_var  = "yvar",
#'   predictor_vars= c("x1","x2"),
#'   k             = 5,
#'   tess_method   = "hex",
#'   tess_args     = list(target_cells = 40),
#'   models        = c("GWR","Bayesian")
#' )
#'
#' res$overall     # model comparison
#' res$by_fold     # per-fold metrics
#' res$tessellation$cells  # diagnostic grid polygons
#' }
#'
#' @seealso \code{\link{cv_gwr}}, \code{\link{cv_bayes}},
#'   \code{\link{create_grid_polygons}}, \code{\link{build_tessellation}},
#'   \code{\link{ensure_projected}}, \code{\link{coerce_to_points}}
#'
#' @export
evaluate_models_cv <- function(
    data_sf,
    response_var,
    predictor_vars,
    k            = 5,
    seed         = 123,
    folds        = NULL,
    boundary     = NULL,
    pointize     = "auto",
    tess_method  = c("grid", "hex", "square", "voronoi", "triangles"),
    tess_args    = list(),
    summary      = c("mean","median"),
    models       = c("GWR","Bayesian"),
    quiet        = FALSE
) {
  summary     <- match.arg(summary)
  tess_method <- match.arg(tess_method)
  .msg <- function(...) if (!quiet) message(...)
  .has <- function(pkg) requireNamespace(pkg, quietly = TRUE)
  .filter_args <- function(fun, args_list) {
    fml <- try(formals(fun), silent = TRUE)
    if (inherits(fml, "try-error")) return(list())
    ok <- names(fml)
    args_list[names(args_list) %in% ok]
  }
  
  # ---- guards ----
  if (!inherits(data_sf, "sf")) stop("evaluate_models_cv(): `data_sf` must be an sf object.")
  if (!is.character(response_var) || length(response_var) != 1L)
    stop("evaluate_models_cv(): `response_var` must be a single column name.")
  if (!all(predictor_vars %in% names(sf::st_drop_geometry(data_sf))))
    stop("evaluate_models_cv(): some `predictor_vars` are missing in `data_sf`.")
  
  models <- unique(intersect(models, c("GWR","Bayesian")))
  if (length(models) == 0L) models <- "GWR"
  
  has_spgwr <- .has("spgwr")
  has_brms  <- .has("brms")
  
  has_cv_gwr   <- exists("cv_gwr",   mode = "function")
  has_cv_bayes <- exists("cv_bayes", mode = "function")
  
  if ("GWR" %in% models && !(has_spgwr && has_cv_gwr)) {
    .msg("evaluate_models_cv(): dropping GWR (needs spgwr + cv_gwr()).")
    models <- setdiff(models, "GWR")
  }
  if ("Bayesian" %in% models && !(has_brms && has_cv_bayes)) {
    .msg("evaluate_models_cv(): dropping Bayesian (needs brms + cv_bayes()).")
    models <- setdiff(models, "Bayesian")
  }
  if (length(models) == 0L) stop("evaluate_models_cv(): no viable models after availability checks.")
  .msg("evaluate_models_cv(): using models = [", paste(models, collapse = ", "), "]")
  
  # Ensure internal row id for CV implementations that rely on original indexing
  if (!("..row_id" %in% names(data_sf))) data_sf$`..row_id` <- seq_len(nrow(data_sf))
  
  # ---------------- Tessellation (diagnostic/return) ----------------
  # Build/normalize boundary in a projected CRS
  proj_pts <- try({
    pts_tmp <- if (inherits(sf::st_geometry_type(data_sf), "sfc_POINT")) data_sf else coerce_to_points(data_sf, pointize)
    ensure_projected(pts_tmp)
  }, silent = TRUE)
  if (inherits(proj_pts, "try-error")) proj_pts <- ensure_projected(data_sf)
  
  boundary_proj <- if (!is.null(boundary)) ensure_projected(boundary, proj_pts) else {
    g   <- sf::st_geometry(proj_pts)
    hull <- sf::st_convex_hull(sf::st_union(g))
    # expand buffer ~1% of bbox diagonal (safe for meters in UTM-ish CRS)
    bb   <- sf::st_bbox(hull)
    diag <- sqrt((bb$xmax - bb$xmin)^2 + (bb$ymax - bb$ymin)^2)
    sf::st_buffer(hull, dist = 0.01 * diag)
  }
  
  cells <- NULL
  if (tess_method %in% c("grid","hex","square")) {
    if (!exists("create_grid_polygons", mode = "function"))
      stop("evaluate_models_cv(): create_grid_polygons() not found.")
    args <- modifyList(list(boundary = boundary_proj,
                            target_cells = 30,
                            type = if (tess_method %in% c("hex","square")) tess_method else "square"),
                       tess_args, keep.null = TRUE)
    cells <- do.call(create_grid_polygons, args)
  } else if (tess_method %in% c("voronoi","triangles")) {
    if (!exists("build_tessellation", mode = "function"))
      stop("evaluate_models_cv(): build_tessellation() not found.")
    method_map <- if (tess_method == "voronoi") "voronoi" else "triangles"
    args <- modifyList(list(points = coerce_to_points(proj_pts, "auto"),
                            boundary = boundary_proj,
                            method = method_map),
                       tess_args, keep.null = TRUE)
    bt <- do.call(build_tessellation, args)
    cells <- bt$cells
  } else {
    stop("evaluate_models_cv(): unsupported tessellation method: ", tess_method)
  }
  
  tessellation <- list(method = tess_method, args = tess_args, cells = cells)
  
  # ---------------- Run CV models ----------------
  if (!is.null(seed)) set.seed(seed)
  
  comparison_rows <- list()
  by_fold_rows    <- list()
  
  if ("GWR" %in% models) {
    .msg("evaluate_models_cv(): running CV for GWR ...")
    base_args <- list(
      data_sf        = data_sf,
      response_var   = response_var,
      predictor_vars = predictor_vars,
      folds          = folds,
      k              = k,
      seed           = seed,
      boundary       = boundary,
      pointize       = pointize
      # DO NOT pass fit_args unless cv_gwr() supports it
    )
    args_gwr <- .filter_args(cv_gwr, base_args)
    gwr_cv <- do.call(cv_gwr, args_gwr)
    
    ov <- try(as.data.frame(gwr_cv$overall), silent = TRUE)
    if (inherits(ov, "try-error") || nrow(ov) == 0L)
      ov <- data.frame(RMSE=NA_real_, MAE=NA_real_, R2=NA_real_)
    if (!"n_pred" %in% names(ov)) ov$n_pred <- length(predictor_vars)
    ov$model <- "GWR"
    comparison_rows[["GWR"]] <- ov[, intersect(c("model","RMSE","MAE","R2","n_pred"), names(ov)), drop = FALSE]
    
    bf <- try(as.data.frame(gwr_cv$fold_metrics), silent = TRUE)
    if (!inherits(bf, "try-error") && nrow(bf)) {
      if (!"fold" %in% names(bf)) bf$fold <- seq_len(nrow(bf))
      bf$model <- "GWR"
      by_fold_rows[["GWR"]] <- bf
    }
  }
  
  if ("Bayesian" %in% models) {
    .msg("evaluate_models_cv(): running CV for Bayesian ...")
    base_args <- list(
      data_sf        = data_sf,
      response_var   = response_var,
      predictor_vars = predictor_vars,
      folds          = folds,
      k              = k,
      seed           = seed,
      boundary       = boundary,
      pointize       = pointize,
      summary        = summary
    )
    if ("fit_args" %in% names(formals(cv_bayes))) base_args$fit_args <- list()
    args_bayes <- .filter_args(cv_bayes, base_args)
    bayes_cv <- do.call(cv_bayes, args_bayes)
    
    ov <- try(as.data.frame(bayes_cv$overall), silent = TRUE)
    if (inherits(ov, "try-error") || nrow(ov) == 0L)
      ov <- data.frame(RMSE=NA_real_, MAE=NA_real_, R2=NA_real_)
    if (!"n_pred" %in% names(ov)) ov$n_pred <- length(predictor_vars)
    ov$model <- "Bayesian"
    comparison_rows[["Bayesian"]] <- ov[, intersect(c("model","RMSE","MAE","R2","n_pred"), names(ov)), drop = FALSE]
    
    bf <- try(as.data.frame(bayes_cv$fold_metrics), silent = TRUE)
    if (!inherits(bf, "try-error") && nrow(bf)) {
      if (!"fold" %in% names(bf)) bf$fold <- seq_len(nrow(bf))
      bf$model <- "Bayesian"
      by_fold_rows[["Bayesian"]] <- bf
    }
  }
  
  overall <- dplyr::bind_rows(comparison_rows)
  by_fold <- dplyr::bind_rows(by_fold_rows)
  
  list(
    overall      = overall,
    by_fold      = by_fold,
    tessellation = tessellation
  )
}
