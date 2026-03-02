#' Generate seed points for Voronoi tessellation
#'
#' Creates an sf POINT layer of "seed" locations. Multiple strategies are
#' supported: user-provided points, uniform random sampling within a boundary,
#' or k-means clustering of a sampling cloud.
#'
#' @param boundary Optional polygonal sf object defining the sampling area.
#' @param method One of "kmeans", "random", "provided".
#' @param n Integer; number of seeds to return.
#' @param seeds sf POINT object of user-provided seeds (method = "provided").
#' @param sample_points Optional sf POINT cloud for k-means clustering.
#' @param kmeans_nstart Integer; nstart for kmeans(). Default 10.
#' @param kmeans_iter Integer; iter.max for kmeans(). Default 100.
#' @param set_seed Optional integer RNG seed.
#' @return An sf POINT object with seed_id and method columns.
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
    stop("Argument 'n' is required for method = '", method, "'.", call. = FALSE)
  if (method == "provided") {
    if (is.null(seeds)) stop("Provide 'seeds' (sf POINT) for method = 'provided'.", call. = FALSE)
    .assert_sf(seeds, "POINT", "seeds")
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
      if (!is.null(boundary)) s <- .align_crs(s, boundary)
      s$seed_id <- seq_len(nrow(s))
      s$method  <- "provided"
      s
    },

    "random" = {
      if (is.null(boundary)) stop("Random seeds require 'boundary'.", call. = FALSE)
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
          stop("K-means seeds require 'boundary' (or provide 'sample_points').", call. = FALSE)
        b <- boundary_union(boundary)
        cloud_n <- max(2000L, 50L * as.integer(n))
        cloud_geom <- .robust_st_sample(b, cloud_n)
        cloud_sfc <- sf::st_sfc(cloud_geom, crs = sf::st_crs(boundary))
        cloud <- sf::st_sf(geometry = cloud_sfc)
      }

      xy <- sf::st_coordinates(cloud)
      
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
      
      n_uniq <- nrow(unique(round(xy, 10)))
      k_use <- max(1L, min(as.integer(n), n_uniq))
      if (k_use < n) {
        .log_warn("get_voronoi_seeds(kmeans): requested %d seeds but only %d unique positions; clamping.",
                  n, n_uniq)
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
      if (!is.null(boundary)) s <- .align_crs(s, boundary)
      s
    }
  )

  if (!is.null(boundary) &&
      !is.na(sf::st_crs(boundary)) &&
      !is.na(sf::st_crs(out)) &&
      !(sf::st_crs(out) == sf::st_crs(boundary))) {
    out <- sf::st_transform(out, sf::st_crs(boundary))
  }
  out
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
#' @param points_sf An sf object with POINT geometries.
#' @param k Integer; requested number of clusters.
#' @param set_seed Optional integer RNG seed. Default 456.
#' @return An sf object of k cluster center POINTs.
#' @export
voronoi_seeds_kmeans <- function(points_sf, k, set_seed = 456) {
  # Project to metric CRS if lon/lat to make k-means distance-faithful
  pts_for_km <- points_sf
  pts_is_ll <- .is_longlat(points_sf)
  if (pts_is_ll) {
    pts_for_km <- ensure_projected(points_sf)
  }
  coords <- sf::st_coordinates(pts_for_km)
  n <- nrow(coords)
  
  n_uniq <- nrow(unique(round(coords, 10)))
  k <- max(1L, min(k, n_uniq))

  cleanup <- .with_seed(set_seed)
  on.exit(cleanup(), add = TRUE)
  km <- stats::kmeans(coords, centers = k, iter.max = 50, nstart = 10)
  cent <- as.data.frame(km$centers); names(cent) <- c("x", "y")
  result <- sf::st_as_sf(cent, coords = c("x", "y"), crs = sf::st_crs(pts_for_km))
  # Transform back to original CRS if we projected
  if (pts_is_ll) {
    result <- sf::st_transform(result, sf::st_crs(points_sf))
  }
  result
}


#' Random seed generation within a polygonal boundary
#'
#' @param boundary An sf or sfc polygonal object.
#' @param k Integer; number of random seeds.
#' @param set_seed Integer RNG seed. Default 456.
#' @return An sf object of k random POINTs.
#' @export
voronoi_seeds_random <- function(boundary, k, set_seed = 456) {
  cleanup <- .with_seed(set_seed)
  on.exit(cleanup(), add = TRUE)
  geom <- sf::st_union(boundary)
  pts <- .robust_st_sample(geom, k)
  sf::st_sf(geometry = pts) |> sf::st_set_crs(sf::st_crs(boundary))
}
