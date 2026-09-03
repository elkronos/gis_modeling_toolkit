#' Generate seed points for Voronoi tessellation
#'
#' Creates an sf POINT layer of "seed" locations. Multiple strategies are
#' supported: user-provided points, uniform random sampling within a boundary,
#' or k-means clustering of a sampling cloud.
#'
#' @param boundary Optional polygonal sf object defining the sampling area.
#' @param method One of "kmeans", "random", "provided".
#' @param n Integer; number of seeds to return. Required for
#'   `method = "kmeans"` and `method = "random"`. **Ignored** for
#'   `method = "provided"`, where every row of `seeds` is returned; a mismatch
#'   between `n` and `nrow(seeds)` is reported as a warning.
#'
#'   For `method = "kmeans"` it is an upper bound rather than a guarantee:
#'   k-means cannot produce more centres than there are distinct positions in
#'   the sampling cloud, nor as many centres as there are rows. When `n`
#'   exceeds either ceiling it is clamped, with a warning naming the count
#'   actually used — `n = nrow(sample_points)` is the common case, and yields
#'   `nrow(sample_points) - 1` seeds. Check `nrow()` on the result rather than
#'   assuming `n`.
#' @param seeds sf POINT object of user-provided seeds (method = "provided").
#' @param sample_points Optional sf POINT cloud for k-means clustering.
#' @param kmeans_nstart Integer; nstart for kmeans(). Default 10.
#' @param kmeans_iter Integer; iter.max for kmeans(). Default 100.
#' @param set_seed Optional integer RNG seed.
#' @return An sf POINT object with seed_id and method columns.
#' @family tessellation
#' @examples
#' library(sf)
#' bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
#'   c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
#' ))), crs = 32632))
#' get_voronoi_seeds(bnd, method = "random", n = 5, set_seed = 1)
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

  # --- validation ---
  if (!is.null(boundary)) .assert_sf(boundary, c("POLYGON", "MULTIPOLYGON"), "boundary")
  if (method %in% c("random", "kmeans") && is.null(n))
    stop("get_voronoi_seeds(): argument 'n' is required for method = '", method,
         "'.", call. = FALSE)
  if (method == "provided") {
    if (is.null(seeds))
      stop("get_voronoi_seeds(): provide 'seeds' (sf POINT) for method = 'provided'.",
           call. = FALSE)
    .assert_sf(seeds, "POINT", "seeds")
    if (!is.null(n) && !identical(as.integer(n), nrow(seeds))) {
      .warn_and_log(
        "get_voronoi_seeds(): method = 'provided' ignores 'n'; returning all %d row(s) of 'seeds' rather than the %d requested.",
        nrow(seeds), as.integer(n)
      )
    }
  }
  if (!is.null(sample_points)) .assert_sf(sample_points, "POINT", "sample_points")

  # --- RNG handling (always restored) ---
  cleanup <- .with_seed(set_seed)
  on.exit(cleanup(), add = TRUE)

  boundary_union <- function(b) sf::st_union(.safe_make_valid(b))

  out <- switch(
    method,

    "provided" = {
      s <- seeds
      s$seed_id <- seq_len(nrow(s))
      s$method  <- "provided"
      s
    },

    "random" = {
      if (is.null(boundary))
        stop("get_voronoi_seeds(): random seeds require 'boundary'.", call. = FALSE)
      b <- boundary_union(boundary)
      pts <- .robust_st_sample(b, n)
      pts_sfc <- sf::st_sfc(pts, crs = sf::st_crs(boundary))
      sf::st_sf(seed_id = seq_along(pts_sfc), method = "random",
                geometry = pts_sfc)
    },

    "kmeans" = {
      if (!is.null(sample_points)) {
        cloud <- if (!is.null(boundary)) .align_crs(sample_points, boundary) else sample_points
      } else {
        if (is.null(boundary))
          stop("get_voronoi_seeds(): k-means seeds require 'boundary' (or provide 'sample_points').",
               call. = FALSE)
        b <- boundary_union(boundary)
        cloud_n <- max(2000L, 50L * as.integer(n))
        cloud_geom <- .robust_st_sample(b, cloud_n)
        cloud_sfc <- sf::st_sfc(cloud_geom, crs = sf::st_crs(boundary))
        cloud <- sf::st_sf(geometry = cloud_sfc)
      }

      # If cloud is in lon/lat, project to a local CRS before k-means so
      # that clustering is distance-faithful (k-means in degrees is
      # distorted except in very small areas).
      cloud_for_km <- cloud
      cloud_is_ll <- .is_longlat(cloud)
      if (cloud_is_ll) {
        cloud_for_km <- ensure_projected(cloud)
        .log_info("get_voronoi_seeds(kmeans): projecting cloud from lon/lat before k-means clustering.")
      }
      xy <- sf::st_coordinates(cloud_for_km)
      
      # Two separate ceilings, both of which stats::kmeans() enforces with a
      # raw message the caller cannot act on.  `n_uniq` is "more cluster
      # centers than distinct data points"; nrow(xy) - 1L is "number of cluster
      # centres must lie between 1 and nrow(x)", which fires at k == nrow(x)
      # exactly -- so n = nrow(sample_points) used to die on a raw kmeans error.
      n_uniq <- nrow(unique(round(xy, 10)))
      k_max  <- min(n_uniq, nrow(xy) - 1L)
      k_use  <- max(1L, min(as.integer(n), k_max))
      if (k_use < n) {
        .warn_and_log("get_voronoi_seeds(kmeans): requested %d seeds but the sampling cloud supports at most %d (%d unique position(s) among %d point(s)); clamping.",
                  as.integer(n), k_use, n_uniq, nrow(xy))
      }

      km <- stats::kmeans(x = xy, centers = k_use, iter.max = kmeans_iter,
                          nstart = kmeans_nstart)
      centers <- km$centers
      centers_sfc <- sf::st_sfc(
        lapply(seq_len(nrow(centers)), function(i) sf::st_point(centers[i, ])),
        crs = sf::st_crs(cloud_for_km)
      )
      s <- sf::st_sf(seed_id = seq_len(k_use), method = "kmeans", geometry = centers_sfc)
      # Transform back to original cloud CRS if we projected for k-means
      if (cloud_is_ll && !identical(sf::st_crs(s), sf::st_crs(cloud))) {
        s <- sf::st_transform(s, sf::st_crs(cloud))
      }
      s
    }
  )

  # Single alignment point for every branch: the seeds are returned in the
  # boundary's CRS whenever one is available.  (The branches above no longer
  # align individually, which made this block unreachable.)
  .align_crs(out, boundary)
}


#' Robust wrapper around sf::st_sample that handles exact= failures
#'
#' Tries exact sampling first, then falls back to iterative padding to
#' reach the desired count.
#'
#' @param geom An sfc geometry to sample within.
#' @param n Integer; desired number of sample points.
#' @return An sfc_POINT of length n.
#' @keywords internal
#' @noRd
.robust_st_sample <- function(geom, n) {
  pts <- try(sf::st_sample(geom, size = n, type = "random", exact = TRUE),
             silent = TRUE)
  if (inherits(pts, "try-error")) {
    pts <- sf::st_sample(geom, size = n, type = "random")
  }
  # Pad or trim to exactly n
  max_attempts <- 10L
  attempt <- 0L
  while (length(pts) < n && attempt < max_attempts) {
    attempt <- attempt + 1L
    extra <- sf::st_sample(geom, size = n - length(pts), type = "random")
    pts <- c(pts, extra)
  }
  if (length(pts) > n) pts <- pts[seq_len(n)]
  if (length(pts) < n)
    .log_warn(".robust_st_sample(): could only generate %d of %d requested points.",
              length(pts), n)
  pts
}


#' K-means seed generation from point coordinates
#'
#' Places `k` seed points at k-means cluster centres of the observed
#' coordinates, so seeds — and the Voronoi cells built from them — follow the
#' sampling density: clusters of observations attract seeds, empty ground gets
#' none. Reach for this when you want cells that each carry a comparable number
#' of observations, which is what makes per-cell aggregates in
#' [summarize_by_cell()] similarly precise. Use [voronoi_seeds_random()]
#' instead when you want coverage of the study area rather than of the data,
#' and [get_voronoi_seeds()] to pick between them by name.
#'
#' Lon/lat input is projected first so the k-means distances are metric rather
#' than degrees. Rows with empty or non-finite coordinates are dropped with a
#' warning, and `k` is clamped to the number of distinct positions.
#'
#' @param points_sf An sf object with POINT geometries.
#' @param k Integer; requested number of clusters, and an upper bound rather
#'   than a guarantee. It is clamped, with a warning, to whichever is smaller
#'   of the number of distinct point positions and `nrow(points_sf) - 1` —
#'   k-means can produce neither more centres than there are distinct points
#'   nor as many centres as there are rows. Check `nrow()` on the result.
#' @param set_seed Optional integer RNG seed. Default 456.
#' @return An sf object of **at most** `k` cluster-centre POINTs (fewer when
#'   `k` exceeds the number of distinct positions), with `seed_id` and
#'   `method = "kmeans"` columns matching [get_voronoi_seeds()].
#' @export
voronoi_seeds_kmeans <- function(points_sf, k, set_seed = 456) {
  .assert_sf(points_sf, "POINT", "points_sf")

  # Project to metric CRS if lon/lat to make k-means distance-faithful
  pts_for_km <- points_sf
  pts_is_ll <- .is_longlat(points_sf)
  if (pts_is_ll) {
    pts_for_km <- ensure_projected(points_sf)
  }
  coords <- sf::st_coordinates(pts_for_km)
  if (ncol(coords) >= 2L) coords <- coords[, 1:2, drop = FALSE]
  # st_coordinates() yields one all-NA row per EMPTY POINT rather than zero
  # rows, so a row-count check alone would let empty geometries through to
  # kmeans(), which fails with "NA/NaN/Inf in foreign function call".
  finite_rows <- stats::complete.cases(coords) &
    apply(coords, 1L, function(r) all(is.finite(r)))
  n_drop <- sum(!finite_rows)
  if (n_drop > 0L) {
    .warn_and_log("voronoi_seeds_kmeans(): dropping %d point(s) with empty or non-finite coordinates.",
              n_drop)
    coords <- coords[finite_rows, , drop = FALSE]
  }
  n <- nrow(coords)
  if (n == 0L)
    stop("voronoi_seeds_kmeans(): `points_sf` has no usable coordinates; ",
         "nothing to cluster.", call. = FALSE)

  # As in get_voronoi_seeds(): stats::kmeans() refuses k > distinct rows AND
  # k >= nrow(x), the latter with "number of cluster centres must lie between 1
  # and nrow(x)".  k = nrow(points_sf) therefore used to error rather than
  # clamp, despite `@return` promising "at most k".
  n_uniq <- nrow(unique(round(coords, 10)))
  k_max  <- min(n_uniq, n - 1L)
  k_use  <- max(1L, min(as.integer(k), k_max))
  if (k_use < k) {
    .warn_and_log("voronoi_seeds_kmeans(): requested %d seeds but only %d unique positions among %d point(s); clamping to %d.",
              as.integer(k), n_uniq, n, k_use)
  }

  cleanup <- .with_seed(set_seed)
  on.exit(cleanup(), add = TRUE)
  km <- stats::kmeans(coords, centers = k_use, iter.max = 50, nstart = 10)
  cent <- as.data.frame(km$centers); names(cent) <- c("x", "y")
  result <- sf::st_as_sf(cent, coords = c("x", "y"), crs = sf::st_crs(pts_for_km))
  # Transform back to original CRS if we projected
  if (pts_is_ll) {
    result <- sf::st_transform(result, sf::st_crs(points_sf))
  }
  # Same output contract as get_voronoi_seeds(), so the three seeding
  # functions are drop-in interchangeable.
  result$seed_id <- seq_len(nrow(result))
  result$method  <- "kmeans"
  result
}


#' Random seed generation within a polygonal boundary
#'
#' Draws `k` seed points uniformly at random inside `boundary`, ignoring where
#' the observations are. Reach for this when the cells should cover the study
#' area evenly — so that sparsely sampled ground still gets its own cells and
#' is visibly under-sampled in the results — rather than concentrating
#' resolution where the data already are, which is what
#' [voronoi_seeds_kmeans()] does. It is also the honest choice for a null or
#' sensitivity comparison: re-running an analysis over several random seedings
#' shows how much of a result depends on one particular tessellation.
#'
#' Sampling is by rejection inside the polygon, so an awkward geometry can
#' return fewer than `k` seeds; that shortfall is warned about rather than
#' silently padded.
#'
#' @param boundary An sf or sfc polygonal object.
#' @param k Integer; number of random seeds.
#' @param set_seed Integer RNG seed. Default 456.
#' @return An sf object of **at most** `k` random POINTs (rejection sampling
#'   inside an awkward geometry can fall short of `k`, which is warned about),
#'   with `seed_id` and `method = "random"` columns matching
#'   [get_voronoi_seeds()].
#' @export
voronoi_seeds_random <- function(boundary, k, set_seed = 456) {
  # `@param boundary` documents sf *or* sfc, and .assert_sf() only accepts sf.
  if (inherits(boundary, "sfc")) boundary <- sf::st_as_sf(boundary)
  .assert_sf(boundary, c("POLYGON", "MULTIPOLYGON"), "boundary")

  cleanup <- .with_seed(set_seed)
  on.exit(cleanup(), add = TRUE)
  geom <- sf::st_union(boundary)
  pts <- .robust_st_sample(geom, k)
  out <- sf::st_sf(geometry = pts) |> sf::st_set_crs(sf::st_crs(boundary))
  # Same output contract as get_voronoi_seeds(), so the three seeding
  # functions are drop-in interchangeable.
  out$seed_id <- seq_len(nrow(out))
  out$method  <- "random"
  out
}
