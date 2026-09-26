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
#'   Default 4326. This is the whole mechanism by which the IDs are stable
#'   (sorting in one common CRS is what makes the same layer get the same IDs
#'   whichever projection it arrives in), so if the transform fails the
#'   function says so rather than quietly sorting in the input's own CRS. The
#'   sort key is rounded to 7 decimal degrees (about 1 cm) before ordering, so
#'   the floating-point noise of a round trip through a different projection
#'   does not usually reverse two neighbouring cells. It can where two cells'
#'   centres lie within about that step of the same longitude, as fine cells
#'   stacked north-south near a projection's central meridian do: 36 of 2,500
#'   100 m cells straddling a UTM central meridian changed ID after a
#'   transform to EPSG:3035. No rounding step removes that, so to match cells
#'   computed in different projections, join them on geometry rather than on
#'   the ID. The key is computed on the sphere (s2) whether or not
#'   \code{sf::sf_use_s2()} is on, so the session setting does not change the
#'   IDs. Set to NULL to sort in the input CRS, which gives IDs that are
#'   reproducible but not comparable across projections.
#' @return An sf polygon layer re-ordered with sequential IDs in id_col.
#'   Non-polygonal rows are **dropped** (with a warning), so the result can
#'   have fewer rows than the input; if no polygonal rows remain, an error is
#'   raised.
#' @family tessellation
#' @examples
#' library(sf)
#' bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
#'   c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
#' ))), crs = 32632))
#' g <- create_grid_polygons(bnd, target_cells = 9)
#' # Reverse the rows and re-derive the IDs: the SAME cell gets the same ID,
#' # which is the property a row-position ID does not have.  The check joins
#' # the two layers on geometry, because comparing sorted ID vectors would
#' # pass for any two permutations of 1:9.
#' fwd <- ensure_stable_poly_id(g)
#' rev <- ensure_stable_poly_id(g[nrow(g):1, ])
#' same_cell <- match(st_as_text(st_geometry(rev)), st_as_text(st_geometry(fwd)))
#' all(rev$poly_id == fwd$poly_id[same_cell])
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
  #
  # The key is measured in lon/lat, where sf routes st_centroid() and
  # st_area() to s2, or with sf_use_s2(FALSE) to lwgeom (not a dependency:
  # every Voronoi tessellation, projected ones included, died with "package
  # lwgeom required") and to planar arithmetic on degrees, which differs from
  # the spherical centroid by far more than the rounding step below, so a few
  # near-tied cells took different IDs in s2-on and s2-off sessions (4 of
  # 2,000 Voronoi cells).  s2 is switched on for this function only and
  # restored on exit, so the key is the same whatever sf_use_s2() says.
  if (!isTRUE(sf::sf_use_s2())) {
    suppressMessages(sf::sf_use_s2(TRUE))
    on.exit(suppressMessages(sf::sf_use_s2(FALSE)), add = TRUE)
  }
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

  # Validity is a property of the geometry in the CRS it is being measured
  # in, so the make_valid above (in the layer's own CRS) does not carry over.
  # Two vertices a centimetre apart in a projected CRS can land on the same
  # longitude and latitude, and s2 rejects the ring as degenerate: on a
  # clipped hex tessellation of North Carolina, 2 of 18 cells that are valid
  # projected are invalid once transformed, and st_centroid() on one of them
  # aborted this function rather than returning IDs.  Only the sort copy is
  # repaired; the geometry that comes back is the caller's own.
  if (isTRUE(make_valid))
    sort_sf <- .safe_make_valid(sort_sf)

  # Representative points — all paths produce an sfc_POINT vector
  rep_sfc <- switch(method,
    centroid      = suppressWarnings(sf::st_geometry(sf::st_centroid(sort_sf))),
    surface_point = sf::st_geometry(sf::st_point_on_surface(.drop_empty_parts(sort_sf))),
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

  # Round the sort key before ordering.  The whole point of this function is
  # that "the same layer gets the same IDs whichever projection it arrives in",
  # and a raw double comparison cannot deliver that: cells in one column of a
  # grid share an exact x only in the CRS the grid was built in, and after a
  # round trip through another projection they differ by ~1e-11 degrees.  The
  # within-column order was then decided by that noise -- 14 of 16 cells got a
  # different ID depending on which CRS the layer arrived in, which is exactly
  # the failure the function exists to prevent.  7 decimals is about a
  # centimetre of longitude; the transform_for_sort default puts the key in
  # degrees, and the tie-break on area then index keeps the result total.
  # Rounding moves the problem to the step boundaries rather than removing
  # it: two cells whose longitudes differ by less than a step (fine cells in
  # one column near a central meridian) can still round apart in one CRS and
  # together in another -- 36 of 2,500 100 m cells did via EPSG:3035, and 6
  # decimals is no better -- which is why the documentation says "usually".
  kx <- round(xy[, 1], 7L)
  ky <- round(xy[, 2], 7L)
  ord <- do.call(order, list(kx, ky, signif(area, 9L), idx0))

  out <- polygons_sf[ord, , drop = FALSE]
  out[[id_col]] <- seq_len(nrow(out))
  out
}

# -----------------------------------------------------------------------------
# Cache key builder
# -----------------------------------------------------------------------------

#' Build a deterministic cache key from geometry, CRS, and parameters
#'
#' Uses digest on the binary WKB representation of geometry, which is much
#' faster on complex geometries than hashing WKT text.
#'
#' @param boundary An sf or sfc object.
#' @param type Character grid/tessellation type.
#' @param target_cells Approximate desired cells.
#' @param ... Additional parameters affecting the grid.
#' @param version The package version, hashed into the key so that a cache
#'   that outlives a package upgrade (one persisted by a user, say) cannot
#'   hand a grid built by an older \code{create_grid_polygons()} to a newer
#'   one.  An argument only so that tests can vary it.
#' @return A length-1 character vector.
#' @keywords internal
#' @noRd
.cache_key <- function(boundary, type, target_cells, ...,
                       version = .spatialkit_version()) {
  # The CRS's full WKT, not its `input` name.  A layer read from a file with
  # a custom CRS reports a generic name such as "unknown", so two different
  # site-centred CRSs with the same local boundary coordinates shared a key,
  # and the second site was handed the first one's grid, 11,000 km away.
  # The cost is a rebuild when one CRS arrives written two ways (an EPSG
  # code and the equivalent proj string).
  crs_obj <- sf::st_crs(boundary)
  crs_token <- if (!is.null(crs_obj) && !is.na(crs_obj) &&
                   !is.null(crs_obj$wkt) && nzchar(crs_obj$wkt)) crs_obj$wkt
               else "NA_CRS"

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
  # The "spatialkit_grid::" prefix marks the entry as ours, so
  # clear_grid_cache() can tell its own bindings from anything else living in
  # the environment it is handed.
  paste0("spatialkit_grid::", type, "::",
         digest::digest(list(geom_hash = geom_hash, crs = crs_token,
                             target_cells = target_cells, args = dots,
                             version = as.character(version))))
}

# The installed version string; a separate function so that the key can be
# computed with the package loaded any way (installed or via pkgload).
.spatialkit_version <- function() {
  v <- tryCatch(as.character(utils::packageVersion("spatialkit")),
                error = function(e) NA_character_)
  if (is.na(v)) "unknown" else v
}

# Insertion order of each cache environment's entries, so that the oldest can
# be evicted when the cache is full.  Kept here, keyed by the environment's
# identity, rather than as a binding inside the cache environment: a user
# who hands over their own environment gets nothing written into it but the
# grids.  A recycled address is harmless because the order is re-derived
# from what is still bound before it is used, and .set_cache_order() registers
# a finalizer so an entry does not outlive the environment it describes.
.gmt_cache_meta <- new.env(parent = emptyenv())

.cache_env_id <- function(cache_env) format(cache_env)

.cache_order <- function(cache_env) {
  o <- get0(.cache_env_id(cache_env), envir = .gmt_cache_meta, inherits = FALSE)
  if (is.character(o)) o else character(0)
}

.set_cache_order <- function(cache_env, order) {
  id <- .cache_env_id(cache_env)
  # First time we record anything for this environment, arrange for its entry
  # to die with it.  Without this the registry grew one permanent character
  # vector per cache environment ever used -- the environments themselves are
  # collected, so no caller could name them to clear them, and
  # clear_grid_cache() only removes the entry for the environment it is
  # handed.  The finalizer takes the environment as its ARGUMENT rather than
  # closing over it, so registering it does not keep the environment alive.
  if (!exists(id, envir = .gmt_cache_meta, inherits = FALSE))
    reg.finalizer(cache_env, function(e) {
      eid <- format(e)
      if (exists(eid, envir = .gmt_cache_meta, inherits = FALSE))
        rm(list = eid, envir = .gmt_cache_meta)
    }, onexit = FALSE)
  assign(id, order, envir = .gmt_cache_meta)
  invisible(order)
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
#' @param target_cells Approximate desired number of cells. Default `NULL`,
#'   as in [create_grid_polygons()], so the grid can be sized by `cellsize`
#'   or `n` passed through `...` instead.
#' @param type Grid type: `"square"` (the default) or `"hex"`, matching
#'   [create_grid_polygons()].
#' @param ... Additional arguments forwarded to create_grid_polygons().
#' @param cache_env Environment for memoized grids. Default .gmt_cache.
#' @param max_entries Maximum number of grids the cache holds.  Default 50.
#'   Once full, adding a grid evicts the one added earliest, so a loop over
#'   many boundaries holds at most this many grids (about 2 MB per 2,500-cell
#'   grid) rather than every grid it ever built for the life of the session.
#'   \code{\link{clear_grid_cache}} empties it outright.
#' @return An sf data frame with a stable poly_id column.  The rows are
#'   re-ordered and re-numbered by \code{\link{ensure_stable_poly_id}},
#'   which \code{\link{create_grid_polygons}} does not do: the same cell
#'   therefore carries a different \code{poly_id} depending on which of the two
#'   builders produced it.  Use one builder throughout an analysis; joining a
#'   summary keyed on IDs from one onto geometries from the other draws the
#'   values on the wrong polygons.
#' @family tessellation
#' @family package options and caches
#' @examples
#' library(sf)
#' bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
#'   c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
#' ))), crs = 32632))
#' g <- create_grid_polygons_cached(bnd, target_cells = 16, type = "hex")
#' nrow(g)
#' # The IDs come from ensure_stable_poly_id(), so they follow the geometry:
#' # the same request with the boundary's vertices in another order gives the
#' # same ID to the same cell.
#' bnd2 <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
#'   c(100, 100), c(0, 100), c(0, 0), c(100, 0), c(100, 100)
#' ))), crs = 32632))
#' g2 <- create_grid_polygons_cached(bnd2, target_cells = 16, type = "hex")
#' same_cell <- match(st_as_text(st_geometry(g2)), st_as_text(st_geometry(g)))
#' all(g2$poly_id == g$poly_id[same_cell])
#' @export
create_grid_polygons_cached <- function(boundary,
                                        target_cells = NULL,
                                        type = c("square", "hex"),
                                        ...,
                                        cache_env = .gmt_cache,
                                        max_entries = 50L) {
  type <- match.arg(type)
  if (!is.numeric(max_entries) || length(max_entries) != 1L ||
      !is.finite(max_entries) || max_entries < 1)
    stop("create_grid_polygons_cached(): `max_entries` must be a single number >= 1.",
         call. = FALSE)
  max_entries <- as.integer(max_entries)

  bnd <- if (inherits(boundary, "sfc")) sf::st_as_sf(boundary) else boundary
  if (!inherits(bnd, "sf"))
    stop("create_grid_polygons_cached(): 'boundary' must be sf/sfc POLYGON/MULTIPOLYGON.")
  # The same projection create_grid_polygons() makes, so a cached grid is laid
  # in the CRS an uncached one would be (see .project_for_grid()).
  bnd <- .project_for_grid(bnd, "create_grid_polygons_cached")

  key <- .cache_key(bnd, type, target_cells, ...)

  if (exists(key, envir = cache_env, inherits = FALSE))
    return(get(key, envir = cache_env, inherits = FALSE))

  # ensure_stable_poly_id() is applied here and NOT in create_grid_polygons(),
  # so the same cell carried a different poly_id depending on which of the two
  # functions built it; mixing them (assign with one, join geometry with the
  # other) drew cell statistics on the wrong polygons for 8 of 9 cells.  The
  # cached builder is documented as memoizing create_grid_polygons(), so the
  # renumbering has to be visible in its documentation -- see @return -- and
  # callers must not mix the two builders for one analysis.
  out <- create_grid_polygons(bnd, target_cells = target_cells, type = type, ...)
  out <- ensure_stable_poly_id(out)

  # Evict the oldest entries first when the cache is full.  The order vector
  # is re-derived from what is actually bound, so an entry removed behind
  # our back (rm() by the user) does not count against the cap.
  order <- .cache_order(cache_env)
  order <- order[vapply(order, exists, logical(1), envir = cache_env,
                        inherits = FALSE)]
  while (length(order) >= max_entries) {
    rm(list = order[1L], envir = cache_env)
    order <- order[-1L]
  }
  assign(key, out, envir = cache_env)
  .set_cache_order(cache_env, c(order, key))
  out
}


#' Clear the in-session grid cache
#'
#' Removes all memoized grid results from the internal cache environment.
#'
#' @param cache_env Environment to clear. Default .gmt_cache.
#' @return Invisibly, the number of entries removed.
#' @family package options and caches
#' @examples
#' library(sf)
#' bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
#'   c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
#' ))), crs = 32632))
#' g1 <- create_grid_polygons_cached(bnd, target_cells = 9)
#' g2 <- create_grid_polygons_cached(bnd, target_cells = 9)   # cache hit
#' clear_grid_cache()   # returns the number of entries removed, invisibly
#' print(clear_grid_cache())                                   # 0: already empty
#' @export
clear_grid_cache <- function(cache_env = .gmt_cache) {
  # Only OUR entries.  rm(ls()) wiped every binding in the environment it was
  # given -- a user who passed their own workspace lost unrelated objects and
  # was told they were "entries removed".
  keys <- ls(envir = cache_env, all.names = TRUE)
  keys <- keys[startsWith(keys, "spatialkit_grid::")]
  if (length(keys)) rm(list = keys, envir = cache_env)
  if (exists(.cache_env_id(cache_env), envir = .gmt_cache_meta, inherits = FALSE))
    rm(list = .cache_env_id(cache_env), envir = .gmt_cache_meta)
  invisible(length(keys))
}
