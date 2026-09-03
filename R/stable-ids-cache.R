# -----------------------------------------------------------------------------
# Cache environment
# -----------------------------------------------------------------------------

#' Internal package cache environment
#' @format An environment.
#' @keywords internal
#' @noRd
.gmt_cache <- new.env(parent = emptyenv())

# -----------------------------------------------------------------------------
# Stable polygon IDs
# -----------------------------------------------------------------------------

#' Create deterministic, stable polygon IDs based on spatial sort keys
#'
#' Ensures that a polygon layer has a reproducible, deterministic identifier
#' column by sorting features using representative point coordinates (and
#' secondary tie-breakers) and then assigning sequential IDs.
#'
#' @param polygons_sf An sf or sfc object containing polygonal features.
#' @param id_col Character scalar; name of the identifier column.
#' @param method One of "centroid", "surface_point", "bbox_center".
#' @param make_valid Logical; apply st_make_valid() first. Default TRUE.
#' @param transform_for_sort CRS used only for computing sort-key coordinates.
#'   Default 4326. Set to NULL to disable.
#' @return An sf polygon layer re-ordered with sequential IDs in id_col.
#'   Non-polygonal rows are **dropped** (with a warning), so the result can
#'   have fewer rows than the input; if no polygonal rows remain, an error is
#'   raised.
#' @export
ensure_stable_poly_id <- function(polygons_sf,
                                  id_col = "poly_id",
                                  method = c("centroid", "surface_point",
                                             "bbox_center"),
                                  make_valid = TRUE,
                                  transform_for_sort = 4326) {
  # Normalize to sf
  if (inherits(polygons_sf, "sfc")) polygons_sf <- sf::st_as_sf(polygons_sf)
  if (!inherits(polygons_sf, "sf"))
    stop("ensure_stable_poly_id(): `polygons_sf` must be an sf/sfc object.")

  # Keep only polygon rows
  gtypes <- as.character(sf::st_geometry_type(polygons_sf, by_geometry = TRUE))
  keep   <- gtypes %in% c("POLYGON", "MULTIPOLYGON")
  if (any(!keep)) {
    dropped <- table(gtypes[!keep])
    .warn_and_log(
      "ensure_stable_poly_id(): dropping %d non-polygon row(s) (%s); only POLYGON/MULTIPOLYGON features are given IDs.",
      sum(!keep),
      paste(sprintf("%s: %d", names(dropped), as.integer(dropped)), collapse = ", ")
    )
  }
  polygons_sf <- polygons_sf[keep, , drop = FALSE]
  if (nrow(polygons_sf) == 0L)
    stop("ensure_stable_poly_id(): no polygon rows found.")

  if (isTRUE(make_valid))
    polygons_sf <- .safe_make_valid(polygons_sf)

  method <- match.arg(method)

  # Work on a transformed copy for the sort key.
  #
  # The transform is the whole mechanism by which the IDs are stable: sorting
  # on a common CRS is what makes the same layer get the same IDs whichever
  # projection it arrives in.  Falling back to the untransformed geometry
  # silently therefore does not degrade the result, it defeats the function's
  # purpose -- the IDs stop being comparable with any other run -- so say so.
  sort_sf <- polygons_sf
  if (!is.null(transform_for_sort) && !is.na(sf::st_crs(sort_sf)))
    sort_sf <- tryCatch(
      sf::st_transform(sort_sf, transform_for_sort),
      error = function(e) {
        .log_warn(paste0("ensure_stable_poly_id(): could not transform to the ",
                         "sort CRS (%s) -- %s. Sorting in the layer's own CRS ",
                         "instead, so the IDs assigned here are NOT comparable ",
                         "with IDs assigned to the same features in another ",
                         "projection."),
                  paste(format(transform_for_sort), collapse = " "),
                  conditionMessage(e))
        sort_sf
      })

  # Representative points — all paths produce an sfc_POINT vector
  rep_sfc <- switch(method,
    centroid      = suppressWarnings(sf::st_geometry(sf::st_centroid(sort_sf))),
    surface_point = sf::st_geometry(sf::st_point_on_surface(sort_sf)),
    bbox_center   = {
      geoms <- sf::st_geometry(sort_sf)
      sf::st_sfc(
        lapply(geoms, function(g) {
          bb <- sf::st_bbox(g)
          sf::st_point(c((bb["xmin"] + bb["xmax"]) / 2,
                         (bb["ymin"] + bb["ymax"]) / 2))
        }),
        crs = sf::st_crs(sort_sf)
      )
    }
  )

  # Sort key: x, y, area, original index
  xy <- suppressWarnings(sf::st_coordinates(rep_sfc))
  if (!is.matrix(xy) || nrow(xy) != nrow(polygons_sf))
    stop("ensure_stable_poly_id(): failed to compute representative coordinates.")

  area <- suppressWarnings(as.numeric(sf::st_area(sort_sf)))
  area[!is.finite(area)] <- 0
  idx0 <- seq_len(nrow(polygons_sf))

  ord <- do.call(order, list(xy[, 1], xy[, 2], area, idx0))

  out <- polygons_sf[ord, , drop = FALSE]
  out[[id_col]] <- seq_len(nrow(out))
  out
}

# -----------------------------------------------------------------------------
# Cache key builder
# -----------------------------------------------------------------------------

#' Build a deterministic cache key from geometry, CRS, and parameters
#'
#' Uses digest on the binary WKB representation of geometry rather than WKT
#' text for much better performance on complex geometries.
#'
#' @param boundary An sf or sfc object.
#' @param type Character grid/tessellation type.
#' @param target_cells Approximate desired cells.
#' @param ... Additional parameters affecting the grid.
#' @return A length-1 character vector.
#' @keywords internal
#' @noRd
.cache_key <- function(boundary, type, target_cells, ...) {
  crs_obj <- sf::st_crs(boundary)
  crs_token <- if (!is.null(crs_obj) && !is.na(crs_obj)) {
    inp <- crs_obj$input
    eps <- crs_obj$epsg
    if (!is.null(inp) && !is.na(inp) && nzchar(as.character(inp))) as.character(inp)
    else if (!is.null(eps) && !is.na(eps)) as.character(eps)
    else "NA_CRS"
  } else "NA_CRS"

  # Use binary (WKB) digest for geometry — much faster than WKT for complex shapes
  geom_hash <- tryCatch(
    digest::digest(sf::st_as_binary(sf::st_union(sf::st_geometry(boundary)))),
    error = function(e) {
      # Fallback to WKT if binary fails
      digest::digest(sf::st_as_text(sf::st_union(boundary)))
    }
  )

  dots <- list(...)
  if (length(dots)) {
    nms <- names(dots)
    if (!is.null(nms)) {
      nms[is.na(nms) | !nzchar(nms)] <- ""
      dots <- dots[order(nms)]
    }
  }

  # `target_cells` goes into the hashed payload rather than into the key text:
  # as.integer() truncated it, so 25.2 and 25.7 collided on one key (the second
  # call silently got the first one's grid), and a NULL target_cells collapsed
  # paste0() to character(0), which crashes the exists() lookup downstream.
  paste0(type, "::",
         digest::digest(list(geom_hash = geom_hash, crs = crs_token,
                             target_cells = target_cells, args = dots)))
}

# -----------------------------------------------------------------------------
# Cached grid builder
# -----------------------------------------------------------------------------

#' Create and cache grid polygons over a boundary
#'
#' Builds a grid via \code{create_grid_polygons()} and memoizes the result
#' so repeated calls with the same inputs return instantly.
#'
#' @param boundary An sf or sfc polygonal object.
#' @param target_cells Approximate desired number of cells.
#' @param type Grid type: `"square"` (the default) or `"hex"`, matching
#'   [create_grid_polygons()].
#' @param ... Additional arguments forwarded to create_grid_polygons().
#' @param cache_env Environment for memoized grids. Default .gmt_cache.
#' @return An sf data frame with a stable poly_id column.
#' @export
create_grid_polygons_cached <- function(boundary,
                                        target_cells,
                                        type = c("square", "hex"),
                                        ...,
                                        cache_env = .gmt_cache) {
  type <- match.arg(type)

  bnd <- if (inherits(boundary, "sfc")) sf::st_as_sf(boundary) else boundary
  if (!inherits(bnd, "sf"))
    stop("create_grid_polygons_cached(): 'boundary' must be sf/sfc POLYGON/MULTIPOLYGON.")
  bnd <- ensure_projected(bnd)

  key <- .cache_key(bnd, type, target_cells, ...)

  if (exists(key, envir = cache_env, inherits = FALSE))
    return(get(key, envir = cache_env, inherits = FALSE))

  out <- create_grid_polygons(bnd, target_cells = target_cells, type = type, ...)
  out <- ensure_stable_poly_id(out)

  assign(key, out, envir = cache_env)
  out
}


#' Clear the in-session grid cache
#'
#' Removes all memoized grid results from the internal cache environment.
#'
#' @param cache_env Environment to clear. Default .gmt_cache.
#' @return Invisibly, the number of entries removed.
#' @export
clear_grid_cache <- function(cache_env = .gmt_cache) {
  keys <- ls(envir = cache_env, all.names = TRUE)
  rm(list = keys, envir = cache_env)
  invisible(length(keys))
}
