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

#' Ensure an object is projected (planar) in a sensible CRS
#'
#' If `x` is lon/lat, reprojects to a local metric CRS suitable for planar ops.
#' If `target` is supplied, transforms to `st_crs(target)`. If `x` has `NA`
#' CRS **and** its coordinates do not look like lon/lat, the CRS is left as
#' `NA` (no guessing).
#'
#' @param x An `sf`/`sfc`.
#' @param target Optional `sf`/`sfc` used to copy CRS.
#'
#' @return `x` in a projected CRS (or unchanged if appropriate).
#' @export
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
#' Returns both inputs in the same CRS. If one has a CRS and the other does not,
#' the one with `NA` is transformed/assigned to match. If both are `NA`, both
#' are returned unchanged.
#'
#' @param a,b `sf`/`sfc` objects.
#' @return A list with components `a` and `b` in the same CRS.
#' @export
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

#' Numeric center of feature bounding boxes
#'
#' Computes the center of each feature's bbox without creating polygons,
#' avoiding s2 "degenerate ring" issues for zero-area bboxes (e.g., points).
#'
#' @param x An `sf` object.
#' @return An `sfc_POINT` with `st_crs(x)`.
#' @keywords internal
#' @noRd
.bbox_center_sfc <- function(x) {
  stopifnot(inherits(x, "sf"))
  pts <- lapply(sf::st_geometry(x), function(g) {
    bb <- sf::st_bbox(g)
    sf::st_point(c((bb["xmin"] + bb["xmax"])/2, (bb["ymin"] + bb["ymax"])/2))
  })
  sf::st_sfc(pts, crs = sf::st_crs(x))
}

#' Coerce geometries to representative points (robust, CRS-safe)
#'
#' Converts any `sf` geometry (points, lines, polygons, multiparts) to a
#' `POINT` geometry using one of several strategies.
#'
#' @param x An `sf` object.
#' @param mode Character; one of
#'   `c("auto","centroid","point_on_surface","surface","line_midpoint","bbox_center")`.
#'   (Alias: `"surface"` -> `"point_on_surface"`.)
#' @param tmp_project Logical; when `TRUE` and `x` is longitude/latitude,
#'   temporarily projects to a local metric CRS for planar operations
#'   (centroid/line midpoint). When `FALSE` on lon/lat **lines**,
#'   `"auto"` and `"line_midpoint"` avoid `st_line_sample()` and fall
#'   back to centroids to prevent s2 errors.
#'
#' @return An `sf` with the same rows/attributes as `x`, but `POINT`
#'   geometries in the same CRS as `x`.
#'
#' @details
#' * `"auto"`:
#'   - POINT/MULTIPOINT → representative point (first coordinate).
#'   - LINESTRING → midpoint (planar); if lon/lat and `tmp_project = FALSE`, use centroid fallback.
#'   - MULTILINESTRING → longest segment’s midpoint (planar) when `tmp_project = TRUE`,
#'     otherwise centroid fallback.
#'   - POLYGON/MULTIPOLYGON → `st_point_on_surface()`.
#'   - Other types → centroid fallback.
#' * `"centroid"`: `st_centroid()` (spherical OK).
#' * `"point_on_surface"` (alias `"surface"`): robust point guaranteed on the surface.
#' * `"line_midpoint"`: midpoint for LINESTRING only (errors for MULTILINESTRING).
#' * `"bbox_center"`: numeric center of each feature’s bbox using an internal
#'   numeric approach (no temporary polygons), so it’s safe for zero-area bboxes.
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
    centers <- .bbox_center_sfc(x)  # uses numeric bbox centers; defined elsewhere in this file
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
  is_ll <- .is_longlat(x)
  
  # Helper: midpoint for a single LINESTRING in the CRS of `x`, projecting if needed
  midpoint_linestring <- function(geom) {
    if (is_ll && !tmp_project) {
      # safe spherical fallback
      return(sf::st_centroid(geom))
    } else {
      x_proj <- if (tmp_project) ensure_projected(sf::st_as_sf(geom, crs = crs)) else sf::st_as_sf(geom, crs = crs)
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
        g_i   <- if (tmp_project) sf::st_geometry(ensure_projected(sf::st_as_sf(g[i], crs = crs))) else g[i]
        parts <- suppressWarnings(sf::st_cast(g_i, "LINESTRING"))
        if (length(parts) == 0L) {
          out[[i]] <- sf::st_centroid(g[i])[[1]]
        } else {
          lens <- as.numeric(sf::st_length(parts))
          j    <- if (length(lens)) which.max(lens) else 1L
          mp   <- sf::st_line_sample(parts[j], sample = 0.5)
          mp   <- sf::st_cast(mp, "POINT")
          # back to original CRS if needed
          mp   <- if (!identical(sf::st_crs(parts), crs)) sf::st_transform(mp, crs) else mp
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

# -----------------------------------------------------------------------------
# Voronoi Tessellation
# -----------------------------------------------------------------------------

#' Voronoi polygons from points (duplicates replicated to one tile per input seed)
#'
#' Builds a Voronoi tessellation from input points and (optionally) clips it
#' to a supplied polygon. Exact coordinate duplicates are consolidated for the
#' tessellation step, but the output is then **replicated** so there is **one
#' polygon per input row**; duplicate seeds will share the same polygon geometry.
#'
#' @param points   An `sf` with POINT geometry (other geometries are coerced
#'                 via `coerce_to_points(points, "auto")`).
#' @param clip_with Optional `sf`/`sfc` polygon to clip against. If `NULL`,
#'                 a slightly expanded bbox of the seeds is used to bound the tessellation.
#'
#' @return An `sf` POLYGON layer with one row **per input seed** and a column
#'         `seed_id` giving the original row index (1..nrow(points)).
#' @export
create_voronoi_polygons <- function(points, clip_with = NULL) {
  if (!inherits(points, "sf")) stop("create_voronoi_polygons(): `points` must be an sf object.")
  # Ensure POINTs
  if (!all(sf::st_geometry_type(points, by_geometry = TRUE) %in% c("POINT","MULTIPOINT"))) {
    points <- coerce_to_points(points, "auto")
  }
  
  # Keep original order/index for replication later
  points$..seed_id <- seq_len(nrow(points))
  
  # Work in a projected CRS for robust geometry ops
  crs_in <- sf::st_crs(points)
  pts    <- ensure_projected(points)
  
  # ---- Build clipping polygon (projected) ----
  clip_poly <- NULL
  if (!is.null(clip_with)) {
    clip_poly <- ensure_projected(clip_with, sf::st_crs(pts))
    if (inherits(clip_poly, "sfc")) clip_poly <- sf::st_sf(geometry = clip_poly)
    clip_poly <- sf::st_make_valid(sf::st_union(clip_poly))
  } else {
    bb <- sf::st_bbox(pts)
    w <- as.numeric(bb["xmax"] - bb["xmin"]); if (!is.finite(w)) w <- 1
    h <- as.numeric(bb["ymax"] - bb["ymin"]); if (!is.finite(h)) h <- 1
    bb_exp <- c(
      xmin = bb["xmin"] - 0.01 * w,
      ymin = bb["ymin"] - 0.01 * h,
      xmax = bb["xmax"] + 0.01 * w,
      ymax = bb["ymax"] + 0.01 * h
    )
    clip_poly <- sf::st_as_sfc(bb_exp) |>
      sf::st_set_crs(sf::st_crs(pts)) |>
      sf::st_make_valid()
  }
  
  # ---- Consolidate exact duplicate coordinates (for tessellation only) ----
  C     <- sf::st_coordinates(pts)
  key   <- paste(C[, 1], C[, 2], sep = "|")
  uniq  <- !duplicated(key)
  seeds <- pts[uniq, , drop = FALSE]
  
  if (nrow(seeds) == 0L) stop("create_voronoi_polygons(): no non-duplicate seeds remain.")
  
  # Tiny deterministic jitter for numerical stability (no RNG used)
  bb_s   <- sf::st_bbox(seeds)
  span   <- max(as.numeric(bb_s["xmax"] - bb_s["xmin"]),
                as.numeric(bb_s["ymax"] - bb_s["ymin"]), 1)
  eps    <- span * 1e-9
  sgeom  <- sf::st_geometry(seeds)
  for (i in seq_len(nrow(seeds))) {
    dx <- ((i %% 7) - 3) / 3 * eps
    dy <- (((i * 3) %% 7) - 3) / 3 * eps
    xy <- sf::st_coordinates(sgeom[i]) + c(dx, dy)
    sgeom[i] <- sf::st_point(xy)
  }
  sf::st_geometry(seeds) <- sgeom
  
  # ---- Voronoi tessellation (projected), then clip & keep POLYGON ----
  mp      <- sf::st_union(seeds)  # MULTIPOINT
  v_raw   <- suppressWarnings(sf::st_voronoi(mp, envelope = clip_poly))
  v_poly  <- sf::st_collection_extract(v_raw, "POLYGON")
  v_sf    <- sf::st_as_sf(v_poly)
  
  v_clip  <- suppressWarnings(sf::st_intersection(v_sf, clip_poly))
  v_clip  <- sf::st_collection_extract(v_clip, "POLYGON")
  v_clip  <- sf::st_make_valid(sf::st_as_sf(v_clip))
  
  # ---- Empty fallback: return clip target polygons replicated per seed ----
  if (nrow(v_clip) == 0L) {
    clip_sf <- if (inherits(clip_poly, "sf")) clip_poly else sf::st_as_sf(clip_poly)
    cont <- sf::st_within(pts, clip_sf)
    idx  <- vapply(cont, function(ix) if (length(ix)) ix[1] else NA_integer_, 1L)
    if (any(is.na(idx))) {
      near <- sf::st_nearest_feature(pts[is.na(idx), ], clip_sf)
      idx[is.na(idx)] <- near
    }
    geom_rep <- sf::st_geometry(clip_sf)[idx]
    out <- sf::st_sf(seed_id = points$..seed_id, geometry = geom_rep)
    return(sf::st_transform(out, crs_in))
  }
  
  # ---- Dissolve multi-part clips into one polygon per UNIQUE seed -----------
  # Tag each clipped piece with its nearest unique seed, then dissolve.
  cents      <- sf::st_centroid(v_clip)
  owner_idx  <- sf::st_nearest_feature(cents, seeds)  # 1..n_unique
  v_clip$..owner <- as.integer(owner_idx)
  
  # IMPORTANT FIX: let sf::summarise do the union implicitly (avoid naming geometry)
  v_cells <- suppressWarnings(
    v_clip |>
      dplyr::group_by(..owner) |>
      dplyr::summarise(.groups = "drop")
  )
  v_cells <- sf::st_make_valid(v_cells)
  
  # Sanity: ensure one row per unique seed; pad any missings
  if (nrow(v_cells) != nrow(seeds)) {
    full <- data.frame(..owner = seq_len(nrow(seeds)))
    v_cells <- full |>
      dplyr::left_join(v_cells, by = "..owner")
    
    miss <- which(is.na(sf::st_is_empty(v_cells$geometry)) | sf::st_is_empty(v_cells$geometry))
    if (length(miss)) {
      have <- setdiff(seq_len(nrow(v_cells)), miss)
      if (length(have)) {
        cent_have <- sf::st_centroid(v_cells[have, ])
        near_idx  <- sf::st_nearest_feature(v_cells[miss, ], cent_have)
        sf::st_geometry(v_cells)[miss] <- sf::st_geometry(v_cells[have, ])[near_idx]
      } else {
        env_sf <- if (inherits(clip_poly, "sf")) clip_poly else sf::st_as_sf(clip_poly)
        sf::st_geometry(v_cells)[miss] <- sf::st_geometry(env_sf)[1]
      }
    }
  }
  
  # ---- Replicate to one polygon per ORIGINAL seed (duplicates share geometry) ----
  uniq_keys  <- key[uniq]
  grp_of_row <- match(key, uniq_keys)  # length = nrow(points), values in 1..n_unique
  
  geom_rep <- sf::st_geometry(v_cells)[grp_of_row]
  
  out <- sf::st_sf(
    seed_id = points$..seed_id,
    geometry = sf::st_sfc(geom_rep, crs = sf::st_crs(pts))
  )
  
  sf::st_transform(out, crs_in)
}


# -----------------------------------------------------------------------------
# Grid Tessellations (Hex / Square)
# -----------------------------------------------------------------------------

#' Grid polygons over a boundary with iterative cell-size tuning
#'
#' Generates either a hexagonal or square grid over `boundary`, clips to the
#' boundary, and keeps **only POLYGON** geometries. Cell size is tuned
#' iteratively to approach `target_cells`. Guards against zero/invalid outputs;
#' if clipping yields zero polygons at the chosen size, it expands the search
#' or falls back to the boundary polygon(s).
#'
#' @param boundary `sf`/`sfc` polygon to cover.
#' @param target_cells Integer target number of output cells inside `boundary`.
#' @param type `"hex"` (default) or `"square"`.
#' @param ... Passed to `sf::st_make_grid()` (e.g., `offset`, etc.).
#' @param max_iter Integer, tuning iterations (default 25).
#' @param widen_factor Numeric, how wide the initial search bracket is around
#'   the heuristic cell size (default 8).
#' @return An `sf` POLYGON layer of grid cells intersected with the boundary.
#' @export
create_grid_polygons <- function(boundary,
                                 target_cells,
                                 type = c("hex","square"),
                                 ...,
                                 max_iter = 25,
                                 widen_factor = 8) {
  type <- match.arg(type)
  if (missing(target_cells) || !is.finite(target_cells) || target_cells < 1) {
    target_cells <- 1L
  }
  
  if (inherits(boundary, "sfc")) boundary <- sf::st_sf(geometry = boundary)
  if (!inherits(boundary, "sf")) stop("create_grid_polygons(): `boundary` must be sf/sfc.")
  # Project for stable areas/lengths, unify/make valid
  bnd <- ensure_projected(boundary)
  bpoly <- sf::st_make_valid(sf::st_union(bnd))
  crs0 <- sf::st_crs(boundary)
  
  # Helper: build & count cells for a given cellsize
  build_grid <- function(cs) {
    square_flag <- identical(type, "square")
    g <- try(suppressWarnings(
      sf::st_make_grid(bpoly, cellsize = cs, what = "polygons", square = square_flag, ...)
    ), silent = TRUE)
    if (inherits(g, "try-error") || length(g) == 0) {
      return(list(sf = NULL, n = 0L))
    }
    gsf <- sf::st_as_sf(g)
    # Clip and keep only polygons
    clipped <- try(suppressWarnings(sf::st_intersection(gsf, bpoly)), silent = TRUE)
    if (inherits(clipped, "try-error") || nrow(clipped) == 0) {
      return(list(sf = NULL, n = 0L))
    }
    clipped <- sf::st_collection_extract(clipped, "POLYGON")
    clipped <- sf::st_make_valid(sf::st_as_sf(clipped))
    list(sf = clipped, n = nrow(clipped))
  }
  
  # Heuristic initial size from area / target
  A <- as.numeric(sf::st_area(bpoly))
  if (!is.finite(A) || A <= 0) stop("create_grid_polygons(): boundary area is non-positive.")
  cs0 <- sqrt(A / max(1, target_cells))
  lo <- cs0 / widen_factor
  hi <- cs0 * widen_factor
  
  best <- NULL
  best_err <- Inf
  
  # Ensure bracket produces at least something; expand if needed
  for (j in 1:4) {
    test_lo <- build_grid(lo)
    test_hi <- build_grid(hi)
    if (test_lo$n > 0 || test_hi$n > 0) break
    lo <- lo / 2
    hi <- hi * 2
  }
  
  # Iterative tuning (binary search on cellsize)
  for (iter in seq_len(max_iter)) {
    mid <- sqrt(lo * hi) # geometric mean for scale invariance
    res <- build_grid(mid)
    if (!is.null(res$sf)) {
      err <- abs(res$n - target_cells)
      if (err < best_err) {
        best <- res$sf
        best_err <- err
        # early exit if exact hit
        if (best_err == 0) break
      }
      # Too many cells → cellsize too small → increase it
      if (res$n > target_cells) {
        lo <- mid
      } else {
        hi <- mid
      }
    } else {
      # invalid build → increase cell size to simplify geometry
      lo <- mid
    }
  }
  
  # If we still failed to get anything, try a last-resort coarse grid
  if (is.null(best) || nrow(best) == 0) {
    logger::log_warn("create_grid_polygons(): grid build failed after tuning; returning boundary polygon(s).")
    out <- sf::st_as_sf(bpoly)
    return(sf::st_transform(out, crs0))
  }
  
  sf::st_transform(best, crs0)
}

# Small `%||%` helper used above
`%||%` <- function(a, b) if (!is.null(a)) a else b

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


#' Assign features to polygons (geometry/CRS safe)
#'
#' @title Assign features to polygons
#' @description
#' Spatially assigns each feature in `features_sf` to a cell in `polygons_sf`,
#' returning the input features with an added polygon identifier column.
#' The function is geometry-column agnostic and CRS-safe: it harmonizes CRSs
#' between inputs and restores the original CRS of `features_sf` on output.
#'
#' @details
#' * **Polygon ID handling:** If `polygons_sf` already contains an identifier
#'   column, it is used. Candidate names searched (in order) are
#'   `c(polygon_id_col, "polygon_id", "id", "cell_id", "grid_id", "poly_id")`.
#'   If none exist, a new sequential integer ID is created.
#' * **Join predicate:** The spatial predicate can be customized via `predicate`.
#'   The default `sf::st_within` excludes boundary-only touches; use
#'   `sf::st_intersects` to treat boundary touches as inside.
#' * **Multiple matches / ties:** If a feature matches multiple polygons
#'   (common with inclusive predicates or non-point features), the join uses
#'   `largest = TRUE` to keep the single “largest” match as determined by *sf*.
#' * **Unassigned features:** Set `keep_unassigned = TRUE` to retain features
#'   that do not match any polygon (their polygon ID will be `NA`).
#' * **Input types:** `polygons_sf` may be an `sf` or `sfc` object; `sfc`
#'   inputs are converted to `sf` internally.
#'
#' @param features_sf An `sf` object containing the features to assign.
#' @param polygons_sf An `sf` or `sfc` polygon layer (e.g., a grid/tessellation)
#'   used for assignment.
#' @param polygon_id_col Character string giving the desired name of the ID
#'   column in the output (default: `"polygon_id"`). If a column with this name
#'   exists in `polygons_sf`, it is used; otherwise a new sequential ID is
#'   created under this name.
#' @param keep_unassigned Logical; if `FALSE` (default), rows that do not match
#'   any polygon are dropped. If `TRUE`, such rows are kept with `NA` in the
#'   ID column.
#' @param predicate A binary spatial predicate function from *sf* used to
#'   determine membership (default: `sf::st_within`). Common alternatives:
#'   `sf::st_intersects`, `sf::st_contains`, `sf::st_covered_by`, etc.
#'
#' @return
#' An `sf` object containing the original feature attributes and geometry,
#' plus the polygon ID column named `polygon_id_col`. The CRS is restored to
#' the original CRS of `features_sf` if it was defined.
#'
#' @section Column & geometry behavior:
#' The function never assumes a specific geometry column name; it relies on
#' `sf`’s internal geometry column tracking. Only the polygon identifier column
#' is selected from `polygons_sf` and merged onto the features.
#'
#' @section CRS handling:
#' Input CRSs are aligned via [`harmonize_crs()`], and the output is transformed
#' back to the original CRS of `features_sf` when present.
#'
#' @section Ties & overlaps:
#' When a feature relates to multiple polygons, the join uses
#' `largest = TRUE` to keep a single best match. For points this typically
#' resolves boundary ambiguities; for lines/polygons it selects the largest
#' spatial relation as computed by *sf*.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' set.seed(1)
#'
#' # Random points (WGS84)
#' pts <- st_as_sf(
#'   data.frame(x = runif(50, -95.5, -95.1), y = runif(50, 29.6, 30.0)),
#'   coords = c("x","y"), crs = 4326
#' )
#'
#' # Simple 5x5 square grid over the points' extent
#' bb  <- st_bbox(pts)
#' grd <- st_make_grid(st_as_sfc(bb), n = c(5, 5), what = "polygons")
#' grd <- st_as_sf(grd)
#'
#' # Assign points to grid cells (ID created as 'polygon_id' if absent)
#' res <- assign_features_to_polygons(pts, grd, polygon_id_col = "polygon_id")
#' table(res$polygon_id)
#'
#' # Keep unassigned points and consider boundary-touching points as inside
#' res2 <- assign_features_to_polygons(
#'   pts, grd, keep_unassigned = TRUE, predicate = sf::st_intersects
#' )
#' }
#'
#' @seealso [harmonize_crs()], [sf::st_join()], [sf::st_within()], [sf::st_intersects()]
#' @family tessellation utilities
#' @export
#' @importFrom sf st_join st_within st_intersects st_crs st_transform st_as_sf
assign_features_to_polygons <- function(
    features_sf,
    polygons_sf,
    polygon_id_col = "polygon_id",
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
  id_candidates <- c(polygon_id_col, "polygon_id", "id", "cell_id", "grid_id", "poly_id")
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
  if (id_col != polygon_id_col) {
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

# ---- internal elbow helper (numeric WSS -> elbow) ----
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

#' Determine optimal tessellation levels from data (UAT-compatible)
#'
#' Computes a k-means WSS curve on point coordinates and applies an elbow
#' heuristic to select candidate numbers of cells/seeds (levels).
#'
#' @param data_sf sf; input features. Non-POINT geometries are pointized; data are
#'   projected to a local metric CRS for stable distances.
#' @param max_levels integer; maximum K to evaluate (default 12).
#' @param top_n integer; return up to this many candidate K around the elbow (default 3).
#' @param sample_n integer; subsample size for speed if there are many points (default 1500).
#' @param set_seed integer; seed for k-means reproducibility (default 123).
#'
#' @return Integer vector of candidate levels (e.g., c(k-1, k, k+1) clipped to [1, max_levels]).
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

#' Fit a GWR model with adaptive bandwidth (and AICc)
#'
#' Fits geographically weighted regression using an adaptive kernel. Handles
#' non-point geometries (via \code{coerce_to_points()}), drops rows with
#' non-finite values in modeling columns, reprojects lon/lat to a local
#' projected CRS for stable distances, selects an adaptive bandwidth with
#' \code{spgwr::gwr.sel}, and returns AICc when available.
#'
#' @param data_sf        `sf` with response/predictors and POINT geometry
#'                       (non-points are coerced via `coerce_to_points()`).
#' @param response_var   string; response column name.
#' @param predictor_vars character vector; predictor column names.
#'
#' @return list with elements:
#'   \itemize{
#'     \item \code{model}     : the `spgwr::gwr` fit object
#'     \item \code{bandwidth} : selected adaptive bandwidth in (0,1]
#'     \item \code{adaptive}  : always `TRUE` (this routine uses adaptive kernels)
#'     \item \code{AICc}      : extracted AICc if available, else `NA_real_`
#'   }
#' @export
fit_gwr_model <- function(data_sf, response_var, predictor_vars) {
  if (!inherits(data_sf, "sf")) stop("fit_gwr_model(): `data_sf` must be an sf object.")
  # Ensure point geometry
  if (!all(sf::st_geometry_type(data_sf, by_geometry = TRUE) %in% c("POINT","MULTIPOINT"))) {
    data_sf <- coerce_to_points(data_sf, "auto")
  }
  
  # Keep only complete cases on modeling columns
  keep_cols <- unique(c(response_var, predictor_vars))
  miss <- setdiff(keep_cols, names(data_sf))
  if (length(miss)) stop(sprintf("fit_gwr_model(): missing columns: %s", paste(miss, collapse = ", ")))
  
  cc <- stats::complete.cases(sf::st_drop_geometry(data_sf)[, keep_cols, drop = FALSE])
  if (!all(cc)) {
    if (requireNamespace("logger", quietly = TRUE)) {
      logger::log_warn("fit_gwr_model(): dropping %d rows with NA in modeling columns.", sum(!cc))
    } else {
      warning(sprintf("fit_gwr_model(): dropping %d rows with NA in modeling columns.", sum(!cc)))
    }
  }
  dat <- data_sf[cc, , drop = FALSE]
  
  # ---- CRS guard: project lon/lat to a metric CRS for stable distances ----
  if (.is_longlat(dat)) {
    dat <- ensure_projected(dat)
  }
  
  # Build formula
  fml <- stats::as.formula(paste(response_var, "~", paste(predictor_vars, collapse = " + ")))
  
  # Convert to sp for spgwr
  sp_dat <- methods::as(dat, "Spatial")
  
  # Helper: fixed-distance fallback (kept for potential future use)
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
    eps <- 1e-8 * (mean(diag(V)) + .Machine$double.eps)
    diag(V) <- diag(V) + eps
    out <- mvtnorm::dmvnorm(y, mean = mu, sigma = V, log = TRUE)
  }
  out
}

#' Bayesian spatial regression via spBayes (robust recovery + safe Matérn)
#'
#' Fits a Gaussian spatial regression using \pkg{spBayes} and returns a compact
#' summary plus a finite DIC-like criterion for relative comparisons.
#'
#' @param data_sf        sf POINT layer with response/predictor columns.
#' @param response_var   Character; numeric response column name.
#' @param predictor_vars Character vector; numeric predictor column names.
#' @param n.samples      Integer MCMC iterations for spLM.
#' @param cov_model      One of "exponential","spherical","matern".
#' @param seed           Integer RNG seed.
#' @param verbose        Logical; progress every ~5% if TRUE.
#' @return list(model=list(fit, recover), beta_mean, theta_mean, DIC)
fit_bayesian_spatial_model <- function(
    data_sf,
    response_var,
    predictor_vars,
    n.samples = 1000,
    cov_model = c("exponential","spherical","matern"),
    seed = 123,
    verbose = FALSE
) {
  # --- explicit validation so UAT negative-path pattern matches ---
  allowed <- c("exponential","spherical","matern")
  cm <- tolower(if (length(cov_model)) cov_model[1] else "")
  if (!cm %in% allowed) {
    stop(sprintf(
      "Invalid cov_model '%s'. Supported cov_model values are: %s.",
      as.character(cov_model)[1], paste(allowed, collapse = ", ")
    ))
  }
  cov_model <- cm
  use_nu <- identical(cov_model, "matern")
  
  if (!inherits(data_sf, "sf")) {
    stop("fit_bayesian_spatial_model(): 'data_sf' must be an sf object with POINT geometry.")
  }
  if (!requireNamespace("spBayes", quietly = TRUE)) {
    stop("fit_bayesian_spatial_model(): package 'spBayes' is required. Please install it.")
  }
  
  # helper to safely coerce spRecover outputs to matrices
  as_matrix_safe <- function(x) {
    if (is.null(x)) return(NULL)
    mx <- tryCatch(as.matrix(x), error = function(e) NULL)
    if (!is.null(mx)) return(mx)
    if (is.vector(x) && is.numeric(x)) return(matrix(x, nrow = 1L))
    NULL
  }
  
  gaussian_dic <- function(y, mu, var_eff, draws_mu = NULL, draws_var = NULL) {
    var_eff <- pmax(var_eff, .Machine$double.eps)
    ll_hat <- sum(stats::dnorm(y, mean = mu, sd = sqrt(var_eff), log = TRUE))
    if (!is.null(draws_mu) && !is.null(draws_var) && length(draws_mu) > 0) {
      draws_var <- pmax(draws_var, .Machine$double.eps)
      ll_draw <- numeric(length(draws_mu))
      for (i in seq_along(draws_mu)) {
        ll_draw[i] <- sum(stats::dnorm(y, mean = draws_mu[[i]], sd = sqrt(draws_var[i]), log = TRUE))
      }
      Dbar <- -2 * mean(ll_draw)
      Dhat <- -2 * ll_hat
      dic  <- 2 * Dbar - Dhat
      if (is.finite(dic)) return(dic)
    }
    -2 * ll_hat
  }
  
  coords <- tryCatch(sf::st_coordinates(sf::st_geometry(data_sf)), error = function(e) NULL)
  if (is.null(coords) || ncol(coords) < 2) {
    stop("fit_bayesian_spatial_model(): geometry must be POINT-like; couldn't extract XY coordinates.")
  }
  df <- sf::st_drop_geometry(data_sf)
  if (!response_var %in% names(df)) stop("Response variable not found: ", response_var)
  miss_pred <- setdiff(predictor_vars, names(df))
  if (length(miss_pred) > 0) stop("Predictor(s) not found: ", paste(miss_pred, collapse = ", "))
  
  y   <- as.numeric(df[[response_var]])
  Xdf <- df[, predictor_vars, drop = FALSE]
  
  keep <- is.finite(y)
  if (ncol(Xdf) > 0) keep <- keep & apply(as.matrix(Xdf), 1, function(r) all(is.finite(r)))
  keep <- keep & apply(coords[, 1:2, drop = FALSE], 1, function(r) all(is.finite(r)))
  y      <- y[keep]
  coords <- coords[keep, 1:2, drop = FALSE]
  Xdf    <- Xdf[keep, , drop = FALSE]
  if (length(y) < 5L) stop("Too few valid (finite) observations after cleaning.")
  
  dat <- data.frame(y = y, Xdf, check.names = FALSE)
  rhs <- if (length(predictor_vars)) paste(predictor_vars, collapse = " + ") else "1"
  frm <- stats::as.formula(paste("y ~", rhs))
  X   <- stats::model.matrix(frm, data = dat)
  
  d_nonzero <- as.numeric(dist(coords)); d_nonzero <- d_nonzero[d_nonzero > 0]
  dmax <- if (length(d_nonzero)) max(d_nonzero) else 1
  dmin <- if (length(d_nonzero)) min(d_nonzero) else dmax / 10
  phi_unif <- c(3 / dmax, 3 / dmin)
  
  yvar <- stats::var(y); if (!is.finite(yvar) || yvar <= 0) yvar <- 1
  
  starting <- list(phi = mean(phi_unif), "sigma.sq" = 0.5 * yvar, "tau.sq" = 0.5 * yvar)
  tuning   <- list(phi = 0.1 * starting$phi, "sigma.sq" = 0.1 * starting$`sigma.sq`, "tau.sq" = 0.1 * starting$`tau.sq`)
  priors   <- list("phi.Unif" = phi_unif, "sigma.sq.IG" = c(2, 1), "tau.sq.IG" = c(2, 1))
  if (use_nu) {
    starting$nu <- 1.0
    tuning$nu   <- 0.1
    priors$"nu.Unif" <- c(0.5, 2.5)
  }
  
  set.seed(seed)
  n.report <- if (isTRUE(verbose)) max(50L, floor(n.samples / 20)) else n.samples + 1L
  fit <- spBayes::spLM(
    formula   = frm,
    data      = dat,
    coords    = coords,
    starting  = starting,
    tuning    = tuning,
    priors    = priors,
    cov.model = cov_model,
    n.samples = n.samples,
    n.report  = n.report,
    verbose   = FALSE
  )
  
  burn <- max(100L, floor(n.samples / 2)); burn <- min(burn, n.samples - 2L); if (burn < 0L) burn <- 0L
  rec <- tryCatch(spBayes::spRecover(fit, start = burn, thin = 1, verbose = FALSE), error = function(e) NULL)
  
  beta_samp  <- if (!is.null(rec)) rec$p.beta.samples  else NULL
  if (is.null(beta_samp) && !is.null(rec) && !is.null(rec$beta.samples))  beta_samp  <- rec$beta.samples
  theta_samp <- if (!is.null(rec)) rec$p.theta.samples else NULL
  if (is.null(theta_samp) && !is.null(rec) && !is.null(rec$theta.samples)) theta_samp <- rec$theta.samples
  
  beta_mat  <- as_matrix_safe(beta_samp)
  theta_mat <- as_matrix_safe(theta_samp)
  
  if (is.null(beta_mat) || is.null(theta_mat) || nrow(beta_mat) < 1 || ncol(beta_mat) < 1) {
    ols <- stats::lm(frm, data = dat)
    beta_mean <- stats::coef(ols); beta_mean[!is.finite(beta_mean)] <- 0
    s2 <- stats::sigma(ols)^2
    theta_mean <- list(phi = mean(phi_unif), "sigma.sq" = 0.5 * s2, "tau.sq" = 0.5 * s2)
    if (use_nu) theta_mean$nu <- 1.0
    mu_hat <- as.numeric(X %*% beta_mean)
    var_eff <- theta_mean$tau.sq + 0.25 * theta_mean$`sigma.sq`
    DIC <- gaussian_dic(y, mu_hat, var_eff)
    return(list(model = list(fit = fit, recover = rec),
                beta_mean = beta_mean, theta_mean = theta_mean, DIC = as.numeric(DIC)))
  }
  
  beta_mean_vec  <- colMeans(beta_mat,  na.rm = TRUE)
  theta_mean_vec <- colMeans(theta_mat, na.rm = TRUE)
  theta_mean <- list(
    phi        = unname(theta_mean_vec[["phi"]]),
    "sigma.sq" = unname(theta_mean_vec[["sigma.sq"]]),
    "tau.sq"   = unname(theta_mean_vec[["tau.sq"]])
  )
  if (use_nu && "nu" %in% colnames(theta_mat)) theta_mean$nu <- unname(theta_mean_vec[["nu"]])
  
  mu_hat <- as.numeric(X %*% beta_mean_vec)
  m_all <- nrow(beta_mat); use_m <- min(50L, m_all); take <- unique(round(seq(1, m_all, length.out = use_m)))
  beta_draws  <- beta_mat[take, , drop = FALSE]
  theta_draws <- theta_mat[take, , drop = FALSE]
  eff_var_draws <- as.numeric(theta_draws[, "tau.sq"]) + 0.25 * as.numeric(theta_draws[, "sigma.sq"])
  draws_mu <- lapply(seq_len(nrow(beta_draws)), function(i) as.numeric(X %*% beta_draws[i, ]))
  DIC <- gaussian_dic(y, mu_hat, mean(eff_var_draws), draws_mu, eff_var_draws)
  
  list(
    model      = list(fit = fit, recover = rec),
    beta_mean  = beta_mean_vec,
    theta_mean = theta_mean,
    DIC        = as.numeric(DIC)
  )
}

# -----------------------------------------------------------------------------
# Tessellation Pipeline (Voronoi / Hex / Square)
# -----------------------------------------------------------------------------

#' Build tessellations and assign features to cells
#'
#' Constructs a tessellation over a study area using one of three methods
#' (Voronoi, hexagonal grid, or square grid) at one or more target levels, and
#' assigns input features to the resulting polygons. Robustly harmonizes CRS,
#' pointizes non-point geometries, and (optionally) clips to a supplied
#' boundary. Ensures the returned polygons carry a \code{poly_id} column
#' (1..nrow(polygons)).
#'
#' @param features_sf An \code{sf} object of features to be assigned to
#'   tessellation cells. May contain POINT/LINE/POLYGON geometries; non-point
#'   geometries are converted to representative points via
#'   \code{\link{coerce_to_points}}.
#' @param levels Integer or numeric vector; target numbers of cells (for grids)
#'   or seeds (for Voronoi). Each value is processed independently; the return
#'   is a named list keyed by \code{as.character(level)}.
#' @param method One of \code{"voronoi"}, \code{"hex"}, or \code{"square"}.
#'   \itemize{
#'     \item \strong{voronoi}: Seeds are generated and Voronoi tiles are clipped
#'       to \code{boundary} (if provided) or the features' bounding box.
#'     \item \strong{hex}/\strong{square}: A grid is generated over \code{boundary}
#'       (recommended; if omitted, a bounding-box polygon is used).
#'   }
#' @param boundary Optional \code{sf} or \code{sfc} POLYGON/MULTIPOLYGON used to
#'   constrain/clip the tessellation. If supplied, CRSs are harmonized and
#'   \code{features_sf} is projected to match. If omitted:
#'   \itemize{
#'     \item Voronoi: tiles are clipped to the features' bounding box.
#'     \item Grids: a bounding-box polygon is derived from the features.
#'   }
#' @param seeds Seeding strategy for Voronoi; one of \code{"kmeans"},
#'   \code{"random"}, or \code{"provided"}. Ignored for hex/square grids.
#'   \itemize{
#'     \item \code{"kmeans"}: cluster the (pointized) features into
#'       \code{k = level} groups and use cluster centroids as seeds.
#'     \item \code{"random"}: sample \code{k = level} points uniformly within
#'       \code{boundary} (or the features' bbox if \code{boundary} is missing).
#'     \item \code{"provided"}: use \code{provided_seed_points} as-is (see below).
#'   }
#' @param provided_seed_points Optional \code{sf} POINT/MULTIPOINT object used
#'   when \code{seeds = "provided"}. If the count differs from \code{level},
#'   the function proceeds with the provided count and logs a warning.
#' @param pointize Strategy for converting non-point geometries to points prior
#'   to seeding/assignment; one of \code{"auto"}, \code{"surface"},
#'   \code{"centroid"}, or \code{"line_midpoint"}. See
#'   \code{\link{coerce_to_points}} for details. Default \code{"auto"} chooses a
#'   sensible method by geometry type.
#' @param seed_kmeans_from Deprecated compatibility parameter; currently treated
#'   as \code{"features"} and ignored. Kept to avoid breaking existing code.
#'
#' @details
#' \strong{CRS handling:} All inputs are passed through
#' \code{\link{ensure_projected}} to avoid long/lat distance artifacts. When
#' \code{boundary} is supplied, \code{features_sf} is projected to match;
#' otherwise a local projected CRS is chosen heuristically.
#'
#' \strong{Voronoi path:} Seeds are prepared according to \code{seeds}.
#' Tessellation is built via \code{\link{create_voronoi_polygons}} and clipped
#' to \code{boundary} (if provided) or the features' bounding box. Exact
#' duplicate seeds are dropped upstream and a tiny deterministic jitter is
#' applied internally for numerical stability. If fewer than two seeds remain,
#' a single polygon equal to the clip target is returned. The resulting polygons
#' are ensured to include \code{poly_id}.
#'
#' \strong{Grid path:} Hex/square grids are generated with
#' \code{\link{create_grid_polygons}}, adjusting cell size iteratively to
#' approach \code{target_cells = level}. The output is clipped to
#' \code{boundary} (or a bbox polygon) and restricted to POLYGON geometries.
#' A \code{poly_id} column is added if missing.
#'
#' \strong{Assignment:} Features (pointized if needed) are assigned to polygons
#' using \code{\link{assign_features_to_polygons}} with \code{mode = "intersects"}
#' (rows with no assignment are dropped). The assigned features carry
#' \code{polygon_id} referencing polygons' \code{poly_id}.
#'
#' \strong{Logs \& warnings:} Informative messages are emitted via
#' \pkg{logger} about CRS choices, seed counts, and grid sizing. If
#' \code{seeds = "provided"} and \code{nrow(provided_seed_points) != level},
#' a warning is logged and the provided seeds are used.
#'
#' @return A named \code{list} where each element corresponds to a value in
#'   \code{levels}. Each element is a \code{list} with components:
#'   \itemize{
#'     \item \code{polygons}: \code{sf} POLYGON layer with column \code{poly_id}.
#'     \item \code{data}: \code{sf} POINT layer of assigned features with
#'       \code{polygon_id} referencing \code{poly_id}.
#'     \item \code{method}: the tessellation method used.
#'     \item \code{level}: the processed level (target K / target cells).
#'   }
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
#' vt <- build_tessellation(
#'   features_sf = pts,
#'   levels = 8,
#'   method = "voronoi",
#'   boundary = st_sf(geometry = aoi),
#'   seeds = "kmeans"
#' )
#' plot(st_geometry(vt[["8"]]$polygons))
#' points(st_coordinates(vt[["8"]]$data), pch = 16, cex = 0.6)
#'
#' # Hex grid with ~20 cells
#' ht <- build_tessellation(
#'   features_sf = pts,
#'   levels = 20,
#'   method = "hex",
#'   boundary = st_sf(geometry = aoi)
#' )
#' }
#'
#' @seealso
#' \code{\link{coerce_to_points}},
#' \code{\link{ensure_projected}},
#' \code{\link{harmonize_crs}},
#' \code{\link{create_voronoi_polygons}},
#' \code{\link{create_grid_polygons}},
#' \code{\link{assign_features_to_polygons}},
#' \code{\link{voronoi_seeds_kmeans}},
#' \code{\link{voronoi_seeds_random}}
#'
#' @export


# -----------------------------------------------------------------------------
# Modeling Pipeline Wrapper
# -----------------------------------------------------------------------------

#' Evaluate spatial models across tessellation designs
#'
#' @param data_sf sf; input features with geometry (POINTs preferred). Non-point geometries are pointized via \code{coerce_to_points()}.
#' @param response_var character; name of response column.
#' @param predictor_vars character vector; names of predictor columns.
#' @param levels integer vector or NULL; numbers of cells/seeds by level. If NULL, uses \code{determine_optimal_levels()} if available.
#' @param tessellation character vector; any of \code{c("voronoi","hex","square")}.
#' @param boundary sf POLYGON/MULTIPOLYGON or NULL; if NULL, a bbox fallback is used inside \code{build_tessellation()}.
#' @param seeds character; for Voronoi: \code{"kmeans"}, \code{"random"}, or \code{"provided"} (ignored for grids).
#' @param provided_seed_points sf POINT; required if \code{seeds="provided"}.
#' @param models character vector; any of \code{c("GWR","Bayesian")}.
#' @param n.samples integer; samples for the Bayesian model (passed to \code{fit_bayesian_spatial_model()}).
#' @param cov_model character; covariance model for Bayesian: \code{"exponential"}, \code{"spherical"}, or \code{"matern"}.
#' @param pointize character; strategy for \code{coerce_to_points()} when \code{data_sf} is not POINT, default \code{"auto"}.
#'
#' @return A list with:
#' \describe{
#'   \item{results}{data.frame summarizing metrics by tessellation, level, and model.}
#'   \item{tessellations}{(invisible) the raw tessellation objects returned by \code{build_tessellation()} for reuse.}
#' }
#'
#' @examples
#' # res <- evaluate_models(pts, "resp", c("x1","x2","x3"),
#' #                        levels = c(5,10),
#' #                        tessellation = c("voronoi","hex","square"),
#' #                        boundary = bnd, seeds = "kmeans",
#' #                        models = c("GWR","Bayesian"),
#' #                        n.samples = 400, cov_model = "exponential")
evaluate_models <- function(
    data_sf,
    response_var,
    predictor_vars,
    levels = NULL,
    tessellation = c("voronoi","hex","square"),
    boundary = NULL,
    seeds = c("kmeans","random","provided"),
    provided_seed_points = NULL,
    models = c("GWR","Bayesian"),
    n.samples = 1000L,
    cov_model = c("exponential","spherical","matern"),
    pointize = "auto"
){
  # ---- Arguments & groundwork -------------------------------------------------
  if (!inherits(data_sf, "sf")) stop("data_sf must be an sf object.")
  tessellation <- match.arg(tessellation, several.ok = TRUE)
  seeds        <- match.arg(seeds)
  cov_model    <- match.arg(cov_model)
  
  if (!is.character(response_var) || length(response_var) != 1L)
    stop("response_var must be a single column name.")
  if (!is.character(predictor_vars) || length(predictor_vars) < 1L)
    stop("predictor_vars must be a non-empty character vector.")
  
  # Ensure columns exist
  req_cols <- c(response_var, predictor_vars)
  miss <- setdiff(req_cols, names(data_sf))
  if (length(miss)) stop("Missing required columns in data_sf: ", paste(miss, collapse = ", "))
  
  # Coerce to points if needed (keeps attributes)
  if (!all(sf::st_geometry_type(data_sf, by_geometry = TRUE) %in% c("POINT","MULTIPOINT"))) {
    data_sf <- coerce_to_points(data_sf, pointize)
  }
  
  # Project to a stable planar CRS (align to boundary if given)
  if (!is.null(boundary)) {
    boundary <- ensure_projected(boundary)
    data_sf  <- ensure_projected(data_sf, boundary)
  } else {
    data_sf  <- ensure_projected(data_sf)
  }
  
  # Determine levels when NULL
  if (is.null(levels)) {
    if (exists("determine_optimal_levels", mode = "function")) {
      levels <- determine_optimal_levels(data_sf, max_levels = 12, top_n = 3L)
    } else {
      stop("levels is NULL and determine_optimal_levels() is not available.")
    }
  }
  if (!is.numeric(levels) || any(!is.finite(levels)) || any(levels <= 0))
    stop("levels must be positive finite integers.")
  levels <- as.integer(levels)
  
  # ---- Helpers ----------------------------------------------------------------
  # Drop rows with non-finite predictors/response
  .finite_rows <- function(df, y, x) {
    num <- sf::st_drop_geometry(df)[, c(y, x), drop = FALSE]
    ok1 <- stats::complete.cases(num)
    ok2 <- apply(num, 1L, function(r) all(is.finite(r)))
    ok1 & ok2
  }
  
  .safe_gwr <- function(dat, y, x) {
    out <- try(fit_gwr_model(dat, y, x), silent = TRUE)
    if (inherits(out, "try-error")) return(list(error = attr(out, "condition")$message))
    out
  }
  
  .safe_bayes <- function(dat, y, x, ns, covm) {
    out <- try(fit_bayesian_spatial_model(dat, y, x, n.samples = ns, cov_model = covm), silent = TRUE)
    if (inherits(out, "try-error")) return(list(error = attr(out, "condition")$message))
    out
  }
  
  # ---- Main loop --------------------------------------------------------------
  res_rows <- list()
  tess_store <- list()
  
  for (meth in tessellation) {
    bt <- build_tessellation(
      data_sf = data_sf,
      levels = levels,
      method = meth,
      boundary = boundary,
      seeds = seeds,
      provided_seed_points = provided_seed_points,
      pointize = "auto"
    )
    tess_store[[meth]] <- bt
    
    for (lvl_chr in names(bt)) {
      lvl_obj <- bt[[lvl_chr]]
      assigned <- lvl_obj$data
      
      # Prepare modeling data (keep geometry)
      keep <- .finite_rows(assigned, response_var, predictor_vars)
      dat_model <- assigned[keep, c(predictor_vars, response_var, attr(assigned, "sf_column")), drop = FALSE]
      
      n_points <- nrow(dat_model)
      n_cells  <- nrow(lvl_obj$polygons)
      
      # GWR ---------------------------------------------------------------------
      if ("GWR" %in% models) {
        gwr <- .safe_gwr(dat_model, response_var, predictor_vars)
        if (!is.null(gwr$error)) {
          res_rows[[length(res_rows)+1L]] <- data.frame(
            tessellation = meth,
            level = as.integer(lvl_chr),
            model = "GWR",
            metric = "AICc",
            score = NA_real_,
            bandwidth = NA_real_,
            cov_model = NA_character_,
            n_points = n_points,
            n_cells = n_cells,
            status = "error",
            message = as.character(gwr$error),
            stringsAsFactors = FALSE
          )
        } else {
          res_rows[[length(res_rows)+1L]] <- data.frame(
            tessellation = meth,
            level = as.integer(lvl_chr),
            model = "GWR",
            metric = "AICc",
            score = as.numeric(gwr$AICc),
            bandwidth = if (!is.null(gwr$bandwidth)) as.numeric(gwr$bandwidth) else NA_real_,
            cov_model = NA_character_,
            n_points = n_points,
            n_cells = n_cells,
            status = "ok",
            message = NA_character_,
            stringsAsFactors = FALSE
          )
        }
      }
      
      # Bayesian ---------------------------------------------------------------
      if ("Bayesian" %in% models) {
        bay <- .safe_bayes(dat_model, response_var, predictor_vars, n.samples, cov_model)
        if (!is.null(bay$error)) {
          res_rows[[length(res_rows)+1L]] <- data.frame(
            tessellation = meth,
            level = as.integer(lvl_chr),
            model = "Bayesian",
            metric = "DIC",
            score = NA_real_,
            bandwidth = NA_real_,
            cov_model = cov_model,
            n_points = n_points,
            n_cells = n_cells,
            status = "error",
            message = as.character(bay$error),
            stringsAsFactors = FALSE
          )
        } else {
          res_rows[[length(res_rows)+1L]] <- data.frame(
            tessellation = meth,
            level = as.integer(lvl_chr),
            model = "Bayesian",
            metric = "DIC",
            score = as.numeric(bay$DIC),
            bandwidth = NA_real_,
            cov_model = cov_model,
            n_points = n_points,
            n_cells = n_cells,
            status = "ok",
            message = NA_character_,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  results <- if (length(res_rows)) do.call(rbind, res_rows) else
    data.frame(tessellation = character(),
               level = integer(),
               model = character(),
               metric = character(),
               score = numeric(),
               bandwidth = numeric(),
               cov_model = character(),
               n_points = integer(),
               n_cells = integer(),
               status = character(),
               message = character(),
               stringsAsFactors = FALSE)
  
  # Return shape the UAT expects
  out <- list(results = results)
  attr(out, "tessellations") <- tess_store
  out
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

#' Plot a tessellation colored by count or an outcome summary
#'
#' @param polygons   sf POLYGON/MULTIPOLYGON; one row per cell. If missing
#'                   \code{polygon_id}, a sequential id will be created.
#' @param assignments Optional sf of features assigned to cells (points, lines-as-points, etc.)
#'                   that includes a \code{polygon_id} column. Required for
#'                   \code{fill_mode = "stat"} or when you want count-based fill/labels.
#' @param title      Optional plot title.
#' @param boundary   Optional sf polygon(s) to outline on top of the fill (e.g., study area).
#' @param basemap    Optional sf layer to draw beneath (e.g., counties).
#' @param fill_mode  Either \code{"count"} (default) or \code{"stat"}.
#'                   - "count": fills by number of assigned features per cell.
#'                   - "stat" : fills by a per-cell summary of \code{fill_var}.
#' @param fill_var   Name of numeric column in \code{assignments} to summarize
#'                   when \code{fill_mode = "stat"}.
#' @param fill_stat  Summary to use when \code{fill_mode = "stat"}.
#'                   One of \code{"mean"}, \code{"median"}, \code{"sum"}. Default "mean".
#' @param fill_title Legend title. If \code{NULL}, a sensible default is used.
#' @param show_counts Logical; if \code{TRUE}, adds text labels with per-cell counts
#'                   (from \code{assignments}). Ignored when \code{label} is given.
#' @param label      Optional column name (present in the plotted \code{polygons} after joins)
#'                   to draw as text labels (e.g., "polygon_id", "count", "mean_resp").
#' @param label_round Integer digits for rounding numeric labels. Default 1.
#' @param na_fill_color Color for cells with missing fill values. Default "grey90".
#'
#' @return A \code{ggplot} object.
#' @examples
#' # p <- plot_tessellation_map(polys, assigns, fill_mode="stat", fill_var="resp")
#' # print(p)
plot_tessellation_map <- function(
    polygons,
    data_sf = NULL,                # points/observations assigned to polygons (ideally with polygon_id)
    title = NULL,
    boundary = NULL,               # overlay outline (sf/sfc ok)
    basemap  = NULL,               # overlay features (sf/sfc ok), drawn under polygons
    value_var = NULL,              # e.g. "resp" — if NULL, color by count
    summary_fun = c("mean","median","sum","min","max","sd"),  # how to aggregate value_var
    show_counts = FALSE,           # draw text labels with counts if label_var is NULL
    label_var = NULL,              # which column to label polygons with (e.g. "polygon_id","count","value")
    label_method = c("surface","centroid","bbox_center"),
    fill_alpha = 0.85,
    line_color = "white",
    line_size  = 0.2,
    palette = "viridis",           # "viridis","magma","plasma","inferno","cividis" OR a vector of colors
    legend_title = NULL
) {
  # ---- helpers ---------------------------------------------------------------
  to_sf <- function(x) {
    if (is.null(x)) return(NULL)
    if (inherits(x, "sfc") || inherits(x, "sfg")) x <- sf::st_as_sf(x)
    if (!inherits(x, "sf")) stop("Overlay inputs must be 'sf' or 'sfc'.")
    x
  }
  # choose label geometry method
  label_method <- match.arg(label_method)
  # summary function
  summary_fun <- match.arg(summary_fun)
  
  # ---- inputs to sf + polygons cleanup --------------------------------------
  polygons <- to_sf(polygons)
  # keep only polygonal geometry
  if (!all(sf::st_geometry_type(polygons) %in% c("POLYGON","MULTIPOLYGON"))) {
    polygons <- sf::st_collection_extract(polygons, "POLYGON")
  }
  # project polygons to a sensible local CRS if needed
  polygons <- ensure_projected(polygons)
  
  # target CRS for overlays
  target_crs <- sf::st_crs(polygons)
  tfm <- function(x) {
    x <- to_sf(x)
    if (is.null(x)) return(NULL)
    if (!is.na(sf::st_crs(x))) x <- sf::st_transform(x, target_crs)
    x
  }
  boundary <- tfm(boundary)
  basemap  <- tfm(basemap)
  
  # ensure polygon_id exists
  if (!("polygon_id" %in% names(polygons))) {
    polygons$polygon_id <- seq_len(nrow(polygons))
    warning("polygons lacked 'polygon_id'; created sequential IDs.")
  }
  
  # ---- prepare data and join summaries --------------------------------------
  fill_col <- NULL  # the column we will map to fill
  if (!is.null(data_sf)) {
    data_sf <- to_sf(data_sf)
    # align crs
    if (!is.na(sf::st_crs(data_sf))) data_sf <- sf::st_transform(data_sf, target_crs)
    
    # If data has no polygon_id, assign by spatial join
    if (!("polygon_id" %in% names(data_sf))) {
      data_sf <- suppressWarnings(
        sf::st_join(data_sf, polygons[, c("polygon_id")], left = TRUE)
      )
    }
    
    # aggregate
    if (!is.null(value_var)) {
      if (!(value_var %in% names(data_sf)))
        stop("value_var '", value_var, "' not found in data_sf.")
      # choose aggregation
      fun <- switch(summary_fun,
                    mean = function(z) mean(z, na.rm = TRUE),
                    median = function(z) stats::median(z, na.rm = TRUE),
                    sum = function(z) sum(z, na.rm = TRUE),
                    min = function(z) suppressWarnings(min(z, na.rm = TRUE)),
                    max = function(z) suppressWarnings(max(z, na.rm = TRUE)),
                    sd  = function(z) stats::sd(z, na.rm = TRUE))
      agg <- data_sf |>
        sf::st_drop_geometry() |>
        dplyr::group_by(.data$polygon_id) |>
        dplyr::summarise(
          value = fun(.data[[value_var]]),
          count = dplyr::n(),
          .groups = "drop"
        )
      polygons <- dplyr::left_join(polygons, agg, by = "polygon_id")
      fill_col <- "value"
      if (is.null(legend_title)) legend_title <- paste0(summary_fun, "(", value_var, ")")
    } else {
      # color by count
      agg <- data_sf |>
        sf::st_drop_geometry() |>
        dplyr::count(.data$polygon_id, name = "count")
      polygons <- dplyr::left_join(polygons, agg, by = "polygon_id")
      fill_col <- "count"
      if (is.null(legend_title)) legend_title <- "count"
    }
  } else {
    # no data: just outline polygons (no fill)
    fill_col <- NULL
  }
  
  # ---- label points & label column ------------------------------------------
  lab_sf <- NULL
  lab_col <- NULL
  need_labels <- !is.null(label_var) || isTRUE(show_counts)
  
  if (need_labels) {
    # derive a label column
    if (!is.null(label_var) && (label_var %in% names(polygons))) {
      lab_col <- label_var
    } else if (isTRUE(show_counts) && ("count" %in% names(polygons))) {
      lab_col <- "count"
    } else if ("polygon_id" %in% names(polygons)) {
      lab_col <- "polygon_id"
    }
    
    if (!is.null(lab_col)) {
      # compute label positions
      lab_geom <- switch(
        label_method,
        surface = sf::st_point_on_surface(polygons),
        centroid = suppressWarnings(sf::st_centroid(polygons)),
        bbox_center = {
          bb <- sf::st_bbox(polygons)
          # fallback: per-polygon bbox centers
          sf::st_as_sf(
            data.frame(polygon_id = polygons$polygon_id),
            coords = c("X","Y"), crs = target_crs
          )
        }
      )
      # lab_geom inherits columns of polygons; keep only the label col
      lab_sf <- lab_geom[, c("polygon_id", lab_col)]
      names(lab_sf)[names(lab_sf) == lab_col] <- ".__label__"
    }
  }
  
  # ---- ggplot ---------------------------------------------------------------
  p <- ggplot2::ggplot()
  
  if (!is.null(basemap)) {
    p <- p + ggplot2::geom_sf(data = basemap, fill = NA, color = "grey80", linewidth = 0.25)
  }
  if (!is.null(boundary)) {
    p <- p + ggplot2::geom_sf(data = boundary, fill = NA, color = "grey20", linewidth = 0.4)
  }
  
  if (!is.null(fill_col) && (fill_col %in% names(polygons))) {
    p <- p + ggplot2::geom_sf(
      data = polygons,
      ggplot2::aes(fill = .data[[fill_col]]),
      color = line_color, linewidth = line_size, alpha = fill_alpha
    )
    # fill scale
    if (is.character(palette) && length(palette) == 1L &&
        palette %in% c("viridis","magma","plasma","inferno","cividis")) {
      p <- p + ggplot2::scale_fill_viridis_c(option = palette, na.value = "grey90")
    } else if (is.character(palette) && length(palette) >= 2L) {
      p <- p + ggplot2::scale_fill_gradientn(colours = palette, na.value = "grey90")
    } else {
      # simple fallback
      p <- p + ggplot2::scale_fill_gradient(low = "grey90", high = "black", na.value = "grey95")
    }
    p <- p + ggplot2::labs(fill = legend_title %||% fill_col)
  } else {
    # outline only
    p <- p + ggplot2::geom_sf(
      data = polygons, fill = NA, color = line_color, linewidth = line_size
    )
  }
  
  if (!is.null(lab_sf)) {
    p <- p + ggplot2::geom_sf_text(
      data = lab_sf, ggplot2::aes(label = .data$.__label__), size = 3
    )
  }
  
  p +
    ggplot2::coord_sf(crs = target_crs) +
    ggplot2::labs(title = title) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}


#' Create cross-validation folds for point-referenced data
#'
#' Builds train/test splits for three CV schemes:
#' \itemize{
#'   \item \code{random_kfold}: standard random k-fold on points.
#'   \item \code{block_kfold}: spatially blocked k-fold using a regular grid; whole blocks
#'         are held out to reduce spatial leakage.
#'   \item \code{buffered_loo}: leave-one-out with an exclusion buffer; for each index i,
#'         the test set is \{i\} and the training set excludes all points within \code{buffer}.
#' }
#'
#' @param points_sf \code{sf} object. Prefer POINT geometry; non-points will be
#'   coerced to representative points via \code{coerce_to_points()}.
#' @param k Integer number of folds (ignored for \code{buffered_loo}, where
#'   \eqn{k = nrow(points_sf)}).
#' @param method One of \code{c("random_kfold","block_kfold","buffered_loo")}.
#' @param seed Optional integer for reproducibility (used to shuffle indices and
#'   to break ties when balancing blocks across folds).
#' @param block_nx,block_ny Optional integers for \code{block_kfold}: grid
#'   dimensions in x/y. If omitted, a near-square grid is chosen automatically.
#' @param block_multiplier Numeric \eqn{>= 1}. When \code{block_nx/block_ny} are
#'   not supplied, the target number of blocks is \code{round(block_multiplier * k)}
#'   to provide multiple blocks per fold. Default \code{3}.
#' @param boundary Optional \code{sf}/\code{sfc} polygon to define block extent
#'   for \code{block_kfold}. If \code{NULL}, the points' bbox is used.
#' @param buffer Numeric distance for \code{buffered_loo} (in the linear units of
#'   the data CRS). Required for \code{buffered_loo}.
#' @param drop_empty_blocks Logical, default \code{TRUE}. Drop grid cells (blocks)
#'   with no points before assigning blocks to folds.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{method}}{The requested method.}
#'   \item{\code{k}}{Number of folds actually produced.}
#'   \item{\code{folds}}{List of length \code{k}; each element is a list with integer
#'         vectors \code{train} and \code{test} (row indices into \code{points_sf}).}
#'   \item{\code{assignment}}{A tibble with one row per point and columns:
#'         \code{row_id}, \code{fold}. For \code{buffered_loo}, \code{fold = row_id}.}
#'   \item{\code{params}}{Echo of key parameters used (grid dims, buffer, seed, etc.).}
#' }
#'
#' @details
#' \strong{CRS and distances.} For methods that rely on distances or buffers
#' (\code{block_kfold} during spatial assignment and \code{buffered_loo}), the
#' data are projected via \code{ensure_projected()} if in lon/lat so that units
#' are linear and numerically stable.
#'
#' \strong{block\_kfold algorithm.} A regular grid is laid over the chosen extent,
#' points are assigned to blocks, and blocks are then greedily allocated to folds
#' to balance (approximately) the number of points per fold (largest blocks placed
#' first onto the currently smallest fold). This version robustly checks whether
#' points lie inside a multi-part \code{boundary} by using an \emph{any-overlap}
#' test across all polygons, not just a single-union-row check.
#'
#' @seealso \code{ensure_projected()}, \code{coerce_to_points()}
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
  
  if (!inherits(points_sf, "sf")) stop("make_folds(): `points_sf` must be an sf object.")
  # Ensure POINTs (preserve attributes) and a projected CRS where needed
  if (!all(sf::st_geometry_type(points_sf, by_geometry = TRUE) %in% c("POINT","MULTIPOINT"))) {
    points_sf <- coerce_to_points(points_sf, "auto")
  }
  # Keep an index to refer back to original row order
  points_sf$..row_id <- seq_len(nrow(points_sf))
  
  # RNG
  if (!is.null(seed)) set.seed(seed)
  
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
      logger::log_warn("make_folds(random_kfold): k (%d) > n (%d); reducing k to n.", k, n)
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
    
    logger::log_info("make_folds(block_kfold): making grid %d x %d (≈ %d blocks) over region.", nx, ny, nx*ny)
    
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
      logger::log_warn("make_folds(block_kfold): %d points not assigned to any block; assigning to nearest grid centroid.", n_na)
      cent <- sf::st_centroid(sf::st_geometry(grid_sf))
      na_idx <- which(is.na(pts$..block_id))
      dmat <- sf::st_distance(sf::st_geometry(pts[na_idx, ]), cent)
      pts$..block_id[na_idx] <- apply(as.matrix(dmat), 1, which.min)
    }
    
    B <- max(pts$..block_id, na.rm = TRUE)
    if (!is.finite(B) || B < 1) stop("make_folds(block_kfold): no usable blocks found.")
    if (B < k) {
      logger::log_warn("make_folds(block_kfold): number of non-empty blocks (%d) < k (%d); reducing k to %d.", B, k, B)
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
    
    logger::log_info("make_folds(buffered_loo): building neighbor lists within buffer = %g (CRS units).", buffer)
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

#' Cross-validated GWR with fold-wise training/prediction
#'
#' Runs geographically weighted regression (GWR) across supplied CV folds:
#' fits on training-only data for each fold, predicts at the fold's test
#' locations, and aggregates RMSE/MAE/R^2 metrics. Returns per-point predictions
#' aligned to original row IDs and per-fold diagnostics (including bandwidth used).
#'
#' @param data_sf An \code{sf} object containing the response/predictors and
#'   point geometries (non-points will be coerced with \code{coerce_to_points()}).
#'   Prefer a projected CRS; otherwise it will be projected via \code{ensure_projected()}.
#' @param response_var String; name of the response column.
#' @param predictor_vars Character vector of predictor column names.
#' @param folds A fold object from \code{make_folds()} (must contain a
#'   \code{folds} list with \code{train}/\code{test} integer indices that refer
#'   to the original \code{data_sf} row order).
#' @param bandwidth Either \code{"auto"} (default) to select per-fold bandwidth
#'   via \code{spgwr::gwr.sel}, or a numeric value. If \code{kernel_adaptive=TRUE},
#'   numeric is interpreted as an adaptive proportion in \eqn{(0,1]} (values
#'   \eqn{\ge} 1 are clamped to 0.99); if \code{kernel_adaptive=FALSE}, numeric
#'   is a fixed distance in CRS units.
#' @param kernel_adaptive Logical; \code{TRUE} for adaptive kernel (recommended),
#'   \code{FALSE} for fixed-bandwidth kernel.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{metrics_overall}}{Tibble with RMSE, MAE, R2 across all test points.}
#'   \item{\code{metrics_by_fold}}{Tibble with fold-level RMSE, MAE, R2, n, and bandwidth used.}
#'   \item{\code{predictions}}{Tibble with \code{row_id} (original row index),
#'         \code{fold}, \code{observed}, \code{predicted}.}
#'   \item{\code{params}}{Echo of key settings (adaptive vs fixed, bandwidth mode).}
#' }
#'
#' @details
#' \strong{CRS \& distances:} GWR relies on Euclidean distances. If \code{data_sf}
#' is in lon/lat, it is projected using \code{ensure_projected()} before modeling.
#'
#' \strong{Bandwidth selection:} With \code{bandwidth="auto"}, each fold calls
#' \code{spgwr::gwr.sel} on the training subset with \code{adapt = kernel_adaptive}.
#' If selection fails (non-finite), fallbacks are used: 0.75 (adaptive) or a robust
#' fixed-distance heuristic based on median pairwise spacing.
#'
#' \strong{Predictions:} Local fits are computed at the test coordinates by
#' passing them as \code{fit.points} to \code{spgwr::gwr}. Predicted values are read
#' from the model's SDF \code{pred} column when present; otherwise, they are
#' reconstructed from the local coefficients and the test design matrix.
#'
#' Rows with missing values in any modeling column are removed before CV, and fold
#' indices are re-mapped to the filtered data via original \code{row_id}.
#'
#' @seealso \code{make_folds()}, \code{fit_gwr_model()}, \code{ensure_projected()}, \code{coerce_to_points()}
#' @export
cv_gwr <- function(data_sf,
                   response_var,
                   predictor_vars,
                   folds,
                   bandwidth = "auto",
                   kernel_adaptive = TRUE) {
  
  # ---- Prep: coerce to points, ensure projected CRS, drop NAs ----
  if (!inherits(data_sf, "sf")) stop("cv_gwr(): `data_sf` must be an sf object.")
  if (!all(sf::st_geometry_type(data_sf, by_geometry = TRUE) %in% c("POINT","MULTIPOINT"))) {
    data_sf <- coerce_to_points(data_sf, "auto")
  }
  data_sf$..row_id <- seq_len(nrow(data_sf))
  
  keep_cols <- unique(c(response_var, predictor_vars))
  miss <- setdiff(keep_cols, names(data_sf))
  if (length(miss)) stop(sprintf("cv_gwr(): missing columns: %s", paste(miss, collapse = ", ")))
  
  cc <- stats::complete.cases(sf::st_drop_geometry(data_sf)[, keep_cols, drop = FALSE])
  if (!all(cc)) {
    logger::log_warn("cv_gwr(): dropping %d rows with NA in modeling columns.", sum(!cc))
  }
  keep_idx <- which(cc)
  dat <- data_sf[cc, , drop = FALSE]
  dat <- ensure_projected(dat)
  
  # ---- Helpers ----
  fml <- stats::as.formula(paste(response_var, "~", paste(predictor_vars, collapse = " + ")))
  
  # robust fixed-bandwidth fallback (median interpoint distance on a sample)
  .fixed_bw_fallback <- function(spobj) {
    coords <- sp::coordinates(spobj)
    n <- nrow(coords)
    s <- if (n > 1000) sample.int(n, 1000) else seq_len(n)
    d <- as.matrix(stats::dist(coords[s, , drop = FALSE]))
    stats::median(d[upper.tri(d)], na.rm = TRUE)
  }
  
  # Try to extract predictions directly; else rebuild from local betas
  .extract_predictions <- function(gwr_obj, new_sf) {
    sdf <- gwr_obj$SDF
    df <- as.data.frame(sdf)
    if ("pred" %in% names(df)) {
      return(df$pred)
    }
    # Reconstruct yhat = X_new %*% beta_local(new)
    X_new <- model.matrix(fml, sf::st_drop_geometry(new_sf))
    mm_cols <- colnames(X_new)
    
    # Normalize names for matching (lowercase, alnum only)
    .norm <- function(x) gsub("[^[:alnum:]]", "", tolower(x))
    # map normalized name -> column position in df
    norm_map <- setNames(seq_along(df), .norm(names(df)))
    
    yhat <- numeric(nrow(new_sf))
    beta_mat <- matrix(NA_real_, nrow(new_sf), ncol(X_new))
    colnames(beta_mat) <- mm_cols
    
    for (j in seq_along(mm_cols)) {
      target <- .norm(mm_cols[j])
      cand <- c(
        target,
        paste0("x", target),
        paste0("x.", target),
        paste0(target, "coef"),
        paste0("coef", target),
        # common intercept variants:
        if (target %in% c("intercept","xintercept","x.intercept.")) c("intercept","xintercept","x.intercept.") else character()
      )
      hit <- cand[cand %in% names(norm_map)]
      if (length(hit)) {
        pos <- as.integer(norm_map[hit[1]])
        colname <- names(df)[pos]
        beta_mat[, j] <- df[[colname]]
      } else {
        # Last resort: exact (unnormalized) name match as seen in some spgwr builds
        alt <- mm_cols[j]
        if (alt %in% names(df)) beta_mat[, j] <- df[[alt]]
      }
    }
    
    if (any(!is.finite(beta_mat))) {
      stop("cv_gwr(): could not reconstruct predictions from local coefficients; 'pred' column not found.")
    }
    as.vector(rowSums(beta_mat * X_new))
  }
  
  # ---- Validate folds ----
  if (is.null(folds) || is.null(folds$folds)) {
    stop("cv_gwr(): `folds` must be the object returned by make_folds() (with $folds[[i]]$train/test).")
  }
  
  k <- length(folds$folds)
  all_pred <- list()
  fold_metrics <- vector("list", k)
  
  for (i in seq_len(k)) {
    tr_orig <- folds$folds[[i]]$train
    te_orig <- folds$folds[[i]]$test
    if (is.null(tr_orig) || is.null(te_orig)) next
    
    # Map original indices -> filtered data indices
    train_idx <- which(keep_idx %in% tr_orig)
    test_idx  <- which(keep_idx %in% te_orig)
    
    if (!length(test_idx)) {
      logger::log_warn("cv_gwr(): fold %d has empty test set after filtering; skipping.", i)
      next
    }
    if (!length(train_idx)) {
      logger::log_warn("cv_gwr(): fold %d has empty train set after filtering; skipping.", i)
      next
    }
    
    train_sf <- dat[train_idx, , drop = FALSE]
    test_sf  <- dat[test_idx,  , drop = FALSE]
    
    # sp objects for spgwr
    train_sp <- methods::as(train_sf, "Spatial")
    test_pts <- methods::as(test_sf,  "Spatial")
    
    # --- Bandwidth (per fold) ---
    bw_used <- NA_real_
    if (identical(bandwidth, "auto")) {
      sel <- try(suppressWarnings(spgwr::gwr.sel(fml, data = train_sp, adapt = kernel_adaptive)),
                 silent = TRUE)
      if (inherits(sel, "try-error") || !is.finite(sel) || is.na(sel)) {
        if (isTRUE(kernel_adaptive)) {
          bw_used <- 0.75
          logger::log_warn("cv_gwr(): fold %d bandwidth selection failed; using adaptive=0.75.", i)
        } else {
          bw_used <- .fixed_bw_fallback(train_sp)
          logger::log_warn("cv_gwr(): fold %d bandwidth selection failed; using fixed=%g.", i, bw_used)
        }
      } else {
        bw_used <- as.numeric(sel)
        if (isTRUE(kernel_adaptive) && (bw_used >= 1 || bw_used <= 0 || !is.finite(bw_used))) {
          bw_used <- 0.99
          logger::log_warn("cv_gwr(): fold %d adaptive bw out of bounds; clamping to 0.99.", i)
        }
      }
    } else if (is.numeric(bandwidth)) {
      if (isTRUE(kernel_adaptive)) {
        bw_used <- as.numeric(bandwidth)
        if (!is.finite(bw_used) || bw_used <= 0) bw_used <- 0.75
        if (bw_used >= 1) {
          logger::log_warn("cv_gwr(): fold %d numeric adaptive bw >= 1; clamping to 0.99.", i)
          bw_used <- 0.99
        }
      } else {
        bw_used <- as.numeric(bandwidth)
        if (!is.finite(bw_used) || bw_used <= 0) {
          bw_used <- .fixed_bw_fallback(train_sp)
          logger::log_warn("cv_gwr(): fold %d invalid fixed bw; using fallback=%g.", i, bw_used)
        }
      }
    } else {
      stop("cv_gwr(): `bandwidth` must be \"auto\" or numeric.")
    }
    
    # --- Fit at TEST locations (fit.points) and predict ---
    fit <- try({
      if (isTRUE(kernel_adaptive)) {
        spgwr::gwr(fml, data = train_sp, adapt = bw_used, hatmatrix = FALSE,
                   fit.points = sp::coordinates(test_pts))
      } else {
        spgwr::gwr(fml, data = train_sp, bandwidth = bw_used, hatmatrix = FALSE,
                   fit.points = sp::coordinates(test_pts))
      }
    }, silent = TRUE)
    
    if (inherits(fit, "try-error")) {
      stop(sprintf("cv_gwr(): fold %d GWR fit failed: %s", i, as.character(fit)))
    }
    
    # Extract predictions at test locations
    y_obs <- sf::st_drop_geometry(test_sf)[[response_var]]
    y_hat <- try(.extract_predictions(fit, test_sf), silent = TRUE)
    if (inherits(y_hat, "try-error")) {
      stop(sprintf("cv_gwr(): fold %d could not obtain predictions: %s",
                   i, as.character(y_hat)))
    }
    
    # Metrics for this fold
    resid <- y_obs - y_hat
    rmse <- sqrt(mean(resid^2, na.rm = TRUE))
    mae  <- mean(abs(resid), na.rm = TRUE)
    r2   <- {
      ss_res <- sum(resid^2, na.rm = TRUE)
      ss_tot <- sum((y_obs - mean(y_obs, na.rm = TRUE))^2, na.rm = TRUE)
      if (ss_tot > 0) 1 - ss_res/ss_tot else NA_real_
    }
    
    fold_metrics[[i]] <- tibble::tibble(
      fold = i, n = length(y_obs),
      RMSE = rmse, MAE = mae, R2 = r2,
      bandwidth = bw_used,
      adaptive = kernel_adaptive
    )
    
    all_pred[[i]] <- tibble::tibble(
      row_id    = test_sf$..row_id,
      fold      = i,
      observed  = as.numeric(y_obs),
      predicted = as.numeric(y_hat)
    )
  }
  
  preds <- dplyr::bind_rows(all_pred)
  if (!nrow(preds)) stop("cv_gwr(): no predictions were produced (empty folds after filtering?).")
  
  # Overall metrics
  resid_all <- preds$observed - preds$predicted
  overall <- tibble::tibble(
    RMSE = sqrt(mean(resid_all^2, na.rm = TRUE)),
    MAE  = mean(abs(resid_all), na.rm = TRUE),
    R2   = {
      ss_res <- sum(resid_all^2, na.rm = TRUE)
      ss_tot <- sum((preds$observed - mean(preds$observed, na.rm = TRUE))^2, na.rm = TRUE)
      if (ss_tot > 0) 1 - ss_res/ss_tot else NA_real_
    },
    n = nrow(preds)
  )
  
  list(
    metrics_overall = overall,
    metrics_by_fold = dplyr::bind_rows(fold_metrics),
    predictions     = preds,
    params = list(
      kernel_adaptive = kernel_adaptive,
      bandwidth       = bandwidth
    )
  )
}

#' Cross-validated Bayesian spatial modeling (spBayes)
#'
#' Fits a Bayesian spatial regression (via \pkg{spBayes}) on each CV fold's
#' training set, predicts at the fold's test coordinates with \code{spPredict()},
#' and aggregates RMSE/MAE/R^2. Returns per-point predictions aligned to the
#' original row order. If \pkg{spBayes} isn't available or a fold fit/predict
#' fails, a regression-only fallback is used for that fold.
#'
#' @param data_sf An \code{sf} with response, predictors, and geometry. Non-POINT
#'   geometries are pointized with \code{coerce_to_points()} and data are projected
#'   with \code{ensure_projected()} for sensible distances.
#' @param response_var String; name of response column.
#' @param predictor_vars Character vector; predictor column names.
#' @param folds A folds object from \code{make_folds()} with
#'   \code{$folds[[i]]$train} / \code{$folds[[i]]$test} integer indices that refer
#'   to the original \code{data_sf} row order.
#' @param cov_model One of \code{"exponential"}, \code{"spherical"}, \code{"matern"}.
#' @param n.samples Integer; MCMC samples per fold (used if \pkg{spBayes} is present).
#' @param verbose Logical; forwarded to \code{spBayes::spPredict()}.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{metrics_overall}}{Tibble with overall RMSE, MAE, R2, n.}
#'   \item{\code{metrics_by_fold}}{Tibble with per-fold RMSE, MAE, R2, n and DIC (if computed).}
#'   \item{\code{predictions}}{Tibble with \code{row_id}, \code{fold}, \code{observed}, \code{predicted}, \code{pred_sd} (if available).}
#'   \item{\code{params}}{List echoing \code{cov_model}, \code{n.samples}, \code{verbose}.}
#' }
#' @export
cv_bayes <- function(data_sf,
                     response_var,
                     predictor_vars,
                     folds,
                     cov_model = c("exponential","spherical","matern"),
                     n.samples = 1000L,
                     verbose = FALSE) {
  
  cov_model <- match.arg(cov_model, c("exponential","spherical","matern"))
  
  # ---- Prep: ensure sf/points/CRS, drop NAs, keep row_id for back-mapping ----
  if (!inherits(data_sf, "sf")) stop("cv_bayes(): `data_sf` must be an sf object.")
  if (!all(sf::st_geometry_type(data_sf, by_geometry = TRUE) %in% c("POINT","MULTIPOINT"))) {
    data_sf <- coerce_to_points(data_sf, "auto")
  }
  data_sf$..row_id <- seq_len(nrow(data_sf))
  
  keep_cols <- unique(c(response_var, predictor_vars))
  miss <- setdiff(keep_cols, names(data_sf))
  if (length(miss)) stop(sprintf("cv_bayes(): missing columns: %s", paste(miss, collapse = ", ")))
  
  cc <- stats::complete.cases(sf::st_drop_geometry(data_sf)[, keep_cols, drop = FALSE])
  if (!all(cc)) {
    logger::log_warn("cv_bayes(): dropping %d rows with NA in modeling columns.", sum(!cc))
  }
  keep_idx <- which(cc)
  dat <- data_sf[cc, , drop = FALSE]
  dat <- ensure_projected(dat)
  
  # ---- Validate folds ----
  if (is.null(folds) || is.null(folds$folds)) {
    stop("cv_bayes(): `folds` must be the object returned by make_folds() (with $folds[[i]]$train/test).")
  }
  
  # ---- Helpers ---------------------------------------------------------------
  fml_sp  <- stats::as.formula(paste("y ~ -1 +", paste(predictor_vars, collapse = " + "))) # no intercept for spLM
  fml_lm  <- stats::as.formula(paste(response_var, "~", paste(predictor_vars, collapse = " + "))) # with intercept for fallback
  
  .corr_mat <- function(D, phi, model = cov_model, nu = 1.5) {
    if (model == "exponential") {
      R <- exp(-phi * D)
    } else if (model == "spherical") {
      r <- pmin(phi * D, 1.0)
      R <- (1 - 1.5*r + 0.5*r^3)
      R[D == 0] <- 1
    } else { # matern
      x <- pmax(phi * D, .Machine$double.eps)
      cst <- (2^(1 - nu))/gamma(nu)
      R <- cst * (x^nu) * besselK(x, nu)
      diag(R) <- 1
      R[!is.finite(R)] <- 1
    }
    R
  }
  
  .loglik_gauss <- function(y, mu, Sigma) {
    cholS <- tryCatch(chol(Sigma), error = function(e) NULL)
    if (is.null(cholS)) return(-Inf)
    z <- backsolve(cholS, y - mu, transpose = TRUE)
    -0.5 * (length(y) * log(2*pi) + 2*sum(log(diag(cholS))) + sum(z^2))
  }
  
  # ---- Fold loop -------------------------------------------------------------
  k <- length(folds$folds)
  all_pred  <- vector("list", k)
  fold_tbls <- vector("list", k)
  
  have_spBayes <- requireNamespace("spBayes", quietly = TRUE)
  
  for (i in seq_len(k)) {
    tr_orig <- folds$folds[[i]]$train
    te_orig <- folds$folds[[i]]$test
    if (is.null(tr_orig) || is.null(te_orig)) next
    
    # Map original indices -> filtered data indices
    train_idx <- which(keep_idx %in% tr_orig)
    test_idx  <- which(keep_idx %in% te_orig)
    
    if (!length(test_idx)) {
      logger::log_warn("cv_bayes(): fold %d has empty test set after filtering; skipping.", i)
      next
    }
    if (!length(train_idx)) {
      logger::log_warn("cv_bayes(): fold %d has empty train set after filtering; skipping.", i)
      next
    }
    
    train_sf <- dat[train_idx, , drop = FALSE]
    test_sf  <- dat[test_idx,  , drop = FALSE]
    
    y_tr <- sf::st_drop_geometry(train_sf)[[response_var]]
    X_tr <- model.matrix(fml_sp, sf::st_drop_geometry(train_sf))
    X_te <- model.matrix(fml_sp, sf::st_drop_geometry(test_sf))
    coords_tr <- sf::st_coordinates(train_sf)
    coords_te <- sf::st_coordinates(test_sf)
    
    n_tr <- length(y_tr)
    
    # Default outputs (in case we fall back)
    y_hat  <- rep(NA_real_, nrow(test_sf))
    y_sd   <- rep(NA_real_, nrow(test_sf))
    dic_f  <- NA_real_
    
    if (have_spBayes) {
      # ---- spBayes fit per fold ---------------------------------------------
      # Rough prior bounds for phi from interpoint distances (train only)
      Dtr   <- as.matrix(dist(coords_tr))
      dpos  <- Dtr[Dtr > 0]
      if (length(dpos)) {
        dqs   <- stats::quantile(dpos, probs = c(0.05, 0.95), names = FALSE)
        phi_lo <- 3 / max(dqs, na.rm = TRUE)
        phi_hi <- 3 / max(1e-6, min(dqs, na.rm = TRUE))
      } else {
        phi_lo <- 0.001; phi_hi <- 1
      }
      
      ysd <- stats::sd(y_tr)
      starting <- list("sigma.sq" = max(ysd^2, 1e-4),
                       "tau.sq"   = max(0.1*ysd^2, 1e-6),
                       "phi"      = sqrt(phi_lo * phi_hi))
      tuning   <- list("sigma.sq" = 0.1,
                       "tau.sq"   = 0.1,
                       "phi"      = 0.1)
      priors <- list(
        "beta.Norm"   = list(rep(0, ncol(X_tr)), diag(1e6, ncol(X_tr))),
        "sigma.sq.IG" = c(2, 1),
        "tau.sq.IG"   = c(2, 1),
        "phi.Unif"    = c(phi_lo, phi_hi)
      )
      if (cov_model == "matern") {
        starting$nu <- 1.5
        tuning$nu   <- 0.05
        priors$nu.Unif <- c(0.5, 2.5)
      }
      
      df_fit <- data.frame(y = y_tr, X_tr)
      colnames(df_fit) <- c("y", colnames(X_tr))
      
      fit <- try(
        spBayes::spLM(
          formula   = stats::as.formula(paste("y ~ -1 +", paste(colnames(X_tr), collapse = "+"))),
          data      = df_fit,
          coords    = coords_tr,
          starting  = starting,
          tuning    = tuning,
          priors    = priors,
          cov.model = cov_model,
          n.samples = as.integer(n.samples),
          verbose   = FALSE
        ),
        silent = TRUE
      )
      
      if (!inherits(fit, "try-error")) {
        # ---- Predict at test coords -----------------------------------------
        burn <- max(1L, floor(0.5 * n.samples))
        pred_obj <- try(
          spBayes::spPredict(
            fit,
            pred.covars = X_te,
            pred.coords = coords_te,
            start = burn,
            thin  = 1,
            verbose = verbose
          ),
          silent = TRUE
        )
        
        if (!inherits(pred_obj, "try-error")) {
          Ydraw <- as.matrix(pred_obj$p.y.predictive.samples)
          if (ncol(Ydraw) == nrow(test_sf)) {
            y_hat <- colMeans(Ydraw)
            y_sd  <- apply(Ydraw, 2, stats::sd)
          } else if (nrow(Ydraw) == nrow(test_sf)) {
            y_hat <- rowMeans(Ydraw)
            y_sd  <- apply(Ydraw, 1, stats::sd)
          } else {
            y_hat <- as.numeric(Ydraw)
            y_sd  <- rep(NA_real_, length(y_hat))
          }
        } else {
          logger::log_warn("cv_bayes(): fold %d spPredict() failed; using regression-only fallback.", i)
          # Fallback: regression-only mean using OLS on training data
          fit_lm <- stats::lm(fml_lm, data = sf::st_drop_geometry(train_sf))
          X_te_lm <- model.matrix(fml_lm, sf::st_drop_geometry(test_sf))
          b <- stats::coef(fit_lm)
          b_full <- setNames(numeric(ncol(X_te_lm)), colnames(X_te_lm))
          b_full[names(b)] <- b
          y_hat <- as.vector(X_te_lm %*% b_full)
          y_sd  <- rep(NA_real_, length(y_hat))
        }
        
        # ---- DIC (optional, if samples are available) -----------------------
        # Use parameter samples if present in fit object
        theta_samps <- try(fit$p.theta.samples, silent = TRUE)
        beta_samps  <- try(fit$p.beta.samples,  silent = TRUE)
        if (!inherits(theta_samps, "try-error") && !inherits(beta_samps, "try-error")) {
          theta_samps <- as.matrix(theta_samps)
          beta_samps  <- as.matrix(beta_samps)
          S <- min(nrow(theta_samps), nrow(beta_samps))
          take <- if (S > 200L) unique(round(seq(1, S, length.out = 200L))) else seq_len(S)
          
          dev_sample <- numeric(length(take))
          for (tix in seq_along(take)) {
            s  <- take[tix]
            th <- theta_samps[s, , drop = TRUE]
            sigma.sq <- as.numeric(th[["sigma.sq"]])
            tau.sq   <- as.numeric(th[["tau.sq"]])
            phi      <- as.numeric(th[["phi"]])
            nu       <- if ("nu" %in% colnames(theta_samps)) as.numeric(th[["nu"]]) else 1.5
            R   <- .corr_mat(as.matrix(dist(coords_tr)), phi, cov_model, nu)
            Sig <- sigma.sq * R + diag(tau.sq, n_tr)
            mu  <- as.numeric(X_tr %*% beta_samps[s, ])
            dev_sample[tix] <- -2 * .loglik_gauss(y_tr, mu, Sig)
          }
          
          # Posterior means
          th_bar <- colMeans(theta_samps, na.rm = TRUE)
          beta_bar <- colMeans(beta_samps, na.rm = TRUE)
          sigma.sq.h <- as.numeric(th_bar[["sigma.sq"]])
          tau.sq.h   <- as.numeric(th_bar[["tau.sq"]])
          phi.h      <- as.numeric(th_bar[["phi"]])
          nu.h       <- if ("nu" %in% names(th_bar)) as.numeric(th_bar[["nu"]]) else 1.5
          R.h   <- .corr_mat(as.matrix(dist(coords_tr)), phi.h, cov_model, nu.h)
          Sig.h <- sigma.sq.h * R.h + diag(tau.sq.h, n_tr)
          mu.h  <- as.numeric(X_tr %*% beta_bar)
          dev.h <- -2 * .loglik_gauss(y_tr, mu.h, Sig.h)
          
          dic_f <- if (length(dev_sample)) 2 * mean(dev_sample) - dev.h else NA_real_
        }
      } else {
        logger::log_warn("cv_bayes(): fold %d spLM() failed; using regression-only fallback.", i)
        fit_lm <- stats::lm(fml_lm, data = sf::st_drop_geometry(train_sf))
        X_te_lm <- model.matrix(fml_lm, sf::st_drop_geometry(test_sf))
        b <- stats::coef(fit_lm)
        b_full <- setNames(numeric(ncol(X_te_lm)), colnames(X_te_lm))
        b_full[names(b)] <- b
        y_hat <- as.vector(X_te_lm %*% b_full)
        y_sd  <- rep(NA_real_, length(y_hat))
      }
    } else {
      logger::log_warn("cv_bayes(): spBayes not available; using regression-only fallback for fold %d.", i)
      fit_lm <- stats::lm(fml_lm, data = sf::st_drop_geometry(train_sf))
      X_te_lm <- model.matrix(fml_lm, sf::st_drop_geometry(test_sf))
      b <- stats::coef(fit_lm)
      b_full <- setNames(numeric(ncol(X_te_lm)), colnames(X_te_lm))
      b_full[names(b)] <- b
      y_hat <- as.vector(X_te_lm %*% b_full)
      y_sd  <- rep(NA_real_, length(y_hat))
    }
    
    # ---- Metrics for this fold ----------------------------------------------
    y_obs <- sf::st_drop_geometry(test_sf)[[response_var]]
    resid <- y_obs - y_hat
    rmse <- sqrt(mean(resid^2, na.rm = TRUE))
    mae  <- mean(abs(resid), na.rm = TRUE)
    r2   <- {
      ss_res <- sum(resid^2, na.rm = TRUE)
      ss_tot <- sum((y_obs - mean(y_obs, na.rm = TRUE))^2, na.rm = TRUE)
      if (ss_tot > 0) 1 - ss_res/ss_tot else NA_real_
    }
    
    fold_tbls[[i]] <- tibble::tibble(
      fold = i, n = length(y_obs),
      RMSE = rmse, MAE = mae, R2 = r2,
      cov_model = cov_model,
      DIC = as.numeric(dic_f)
    )
    
    all_pred[[i]] <- tibble::tibble(
      row_id    = test_sf$..row_id,
      fold      = i,
      observed  = as.numeric(y_obs),
      predicted = as.numeric(y_hat),
      pred_sd   = as.numeric(y_sd)
    )
  }
  
  predictions <- dplyr::bind_rows(all_pred)
  if (!nrow(predictions)) stop("cv_bayes(): no predictions were produced (empty folds after filtering?).")
  
  # ---- Overall metrics ----
  resid_all <- predictions$observed - predictions$predicted
  overall <- tibble::tibble(
    RMSE = sqrt(mean(resid_all^2, na.rm = TRUE)),
    MAE  = mean(abs(resid_all), na.rm = TRUE),
    R2   = {
      ss_res <- sum(resid_all^2, na.rm = TRUE)
      ss_tot <- sum((predictions$observed - mean(predictions$observed, na.rm = TRUE))^2, na.rm = TRUE)
      if (ss_tot > 0) 1 - ss_res/ss_tot else NA_real_
    },
    n = nrow(predictions)
  )
  
  list(
    metrics_overall = overall,
    metrics_by_fold = dplyr::bind_rows(fold_tbls),
    predictions     = predictions,
    params = list(
      cov_model = cov_model,
      n.samples = n.samples,
      verbose   = verbose
    )
  )
}

#' Build a tessellation and assign features to cells
#'
#' @param data_sf  sf object; will be pointized via `coerce_to_points(pointize)`.
#' @param levels   Integer vector. For voronoi = number of seeds; for grids = target cell count.
#' @param method   One of "voronoi","hex","square".
#' @param boundary sf **or** sfc polygon (POLYGON/MULTIPOLYGON) or NULL.
#' @param seeds    One of "kmeans","random","provided" (voronoi only).
#' @param provided_seed_points Optional sf POINTS if `seeds='provided'`.
#' @param pointize Mode passed to `coerce_to_points()`.
#' @return A named list; names are levels, each element is a list with `$polygons` and `$data`.
build_tessellation <- function(
    data_sf,
    levels,
    method = c("voronoi","hex","square"),
    boundary = NULL,
    seeds = c("kmeans","random","provided"),
    provided_seed_points = NULL,
    pointize = c("auto","centroid","point_on_surface","bbox_center","line_midpoint")
) {
  method <- match.arg(tolower(method), c("voronoi","hex","square"))
  seeds  <- match.arg(tolower(seeds),  c("kmeans","random","provided"))
  pointize <- match.arg(tolower(pointize),
                        c("auto","centroid","point_on_surface","bbox_center","line_midpoint"))
  
  # --- normalize boundary: accept sfc or sf; require polygonal if provided ----
  normalize_boundary <- function(bnd) {
    if (is.null(bnd)) return(NULL)
    if (inherits(bnd, "sfc")) {
      bnd <- sf::st_as_sf(bnd)  # wrap sfc into an sf data frame
    }
    if (!inherits(bnd, "sf")) {
      stop("boundary must be an sf polygon layer or NULL.")
    }
    gtypes <- as.character(sf::st_geometry_type(bnd))
    if (!all(gtypes %in% c("POLYGON","MULTIPOLYGON"))) {
      stop("boundary must have POLYGON or MULTIPOLYGON geometry.")
    }
    bnd
  }
  
  boundary <- normalize_boundary(boundary)
  
  # --- coerce data to points & harmonize CRS ---------------------------------
  pts <- coerce_to_points(data_sf, mode = pointize)
  if (!is.null(boundary)) {
    hh <- harmonize_crs(pts, boundary)
    pts      <- hh$a
    boundary <- hh$b
  } else {
    # Build a gentle bbox polygon in the pts' CRS for grid/voronoi clipping
    bb  <- sf::st_as_sfc(sf::st_bbox(pts))
    bb  <- sf::st_buffer(bb, dist = 0)  # zero-buffer to close tiny issues
    boundary <- sf::st_as_sf(bb)
  }
  
  out <- list()
  levs <- as.integer(levels)
  if (any(!is.finite(levs) | levs < 1)) stop("All 'levels' must be positive integers.")
  
  for (k in levs) {
    if (method == "voronoi") {
      # choose/prepare seeds
      seed_pts <- switch(
        seeds,
        "kmeans"   = voronoi_seeds_kmeans(pts, k = k),
        "random"   = voronoi_seeds_random(boundary, k = k),
        "provided" = {
          if (is.null(provided_seed_points) || !inherits(provided_seed_points, "sf")) {
            stop("seeds='provided' requires 'provided_seed_points' as sf POINTS.")
          }
          # harmonize CRSs and take first k (or sample) deterministically
          hh <- harmonize_crs(provided_seed_points, boundary)
          sp <- hh$a
          if (nrow(sp) >= k) sp[seq_len(k), , drop = FALSE] else sp
        }
      )
      polys <- create_voronoi_polygons(seed_pts, clip_with = boundary)
    } else {
      # grid methods use 'target_cells' = k
      polys <- create_grid_polygons(boundary, target_cells = k, type = method)
    }
    
    assigned <- assign_features_to_polygons(pts, polys)
    out[[as.character(k)]] <- list(polygons = polys, data = assigned)
  }
  
  out
}
