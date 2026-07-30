#' Prepare and sanitize an sf dataset for spatial modeling
#'
#' Ensures point geometry, projected CRS, and removes rows with missing or
#' non-finite values in modeling columns.  All non-POINT geometries
#' (including MULTIPOINT) are coerced to representative points via
#' \code{coerce_to_points()}, so downstream coordinate extraction always
#' aligns one row per observation.
#'
#' @param data_sf An sf object.
#' @param response_var Response variable column name.
#' @param predictor_vars Predictor column names.
#' @param boundary Optional sf/sfc for CRS alignment.
#' @param pointize Strategy for non-point geometry coercion.
#' @param require_response Logical; if FALSE the response column is not
#'   required to be present (useful for out-of-sample prediction where the
#'   response is unknown).  Default TRUE.
#' @return An sf object (points) in a projected CRS, cleaned.
#' @family model fitting
#' @examples
#' library(sf)
#' dat <- st_as_sf(
#'   data.frame(x = 1:5, y = 5:1,
#'              resp = c(1, 2, NA, 4, 5),
#'              pred = c(1, 2, 3, 4, Inf)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' prep_model_data(dat, "resp", "pred")  # drops rows 3 (NA) and 5 (Inf)
#' @export
prep_model_data <- function(data_sf, response_var, predictor_vars,
                            boundary = NULL,
                            pointize = c("auto", "surface", "centroid",
                                         "line_midpoint", "bbox_center"),
                            require_response = TRUE) {
  if (!inherits(data_sf, "sf"))
    stop("prep_model_data(): 'data_sf' must be an sf object.")
  pointize <- match.arg(pointize)
  if (!is.character(response_var) || length(response_var) != 1L)
    stop("prep_model_data(): 'response_var' must be a single column name.")
  if (!is.character(predictor_vars) || length(predictor_vars) < 1L)
    stop("prep_model_data(): 'predictor_vars' must be a non-empty character vector.")

  req_cols <- if (require_response) c(response_var, predictor_vars) else predictor_vars
  miss <- setdiff(req_cols, names(data_sf))
  if (length(miss))
    stop("prep_model_data(): missing required column(s): ", paste(miss, collapse = ", "))

  # Coerce unless every geometry is already a plain POINT.  MULTIPOINT must
  # be coerced too: a multi-vertex MULTIPOINT survives st_coordinates() with
  # one row per *vertex*, misaligning coordinates with data rows downstream
  # (fit_bayesian_spatial_model()'s cbind, GWmodel's Spatial coercion).
  # coerce_to_points() maps MULTIPOINT to its centroid, so single-vertex
  # MULTIPOINTs (common in shapefiles) are unchanged in practice.
  if (!all(sf::st_geometry_type(data_sf, by_geometry = TRUE) == "POINT")) {
    data_sf <- coerce_to_points(data_sf, pointize)
  }

  if (!is.null(boundary)) {
    bnd <- if (inherits(boundary, "sfc")) sf::st_as_sf(boundary) else boundary
    if (!inherits(bnd, "sf"))
      stop("prep_model_data(): 'boundary' must be sf/sfc when supplied.")
    bnd <- ensure_projected(bnd)
    data_sf <- ensure_projected(data_sf, bnd)
  } else {
    data_sf <- ensure_projected(data_sf)
  }

  # For cleaning, only check columns that are actually present and relevant.
  # When require_response = FALSE (prediction mode), exclude the response

  # variable entirely so that NAs in the response column of reused datasets
  # do not silently drop valid prediction rows.
  check_vars <- if (require_response) c(response_var, predictor_vars) else predictor_vars
  clean_cols <- intersect(check_vars, names(data_sf))
  df <- sf::st_drop_geometry(data_sf)[, clean_cols, drop = FALSE]
  ok_cc <- stats::complete.cases(df)
  num_mask <- vapply(df, is.numeric, logical(1))
  ok_fin <- if (any(num_mask)) {
    apply(as.matrix(df[, num_mask, drop = FALSE]), 1L,
          function(r) all(is.finite(r)))
  } else rep(TRUE, nrow(df))

  keep <- ok_cc & ok_fin
  dropped <- sum(!keep)
  if (dropped > 0)
    .log_warn("prep_model_data(): dropping %d row(s) with non-finite or missing values.", dropped)

  data_sf[keep, , drop = FALSE]
}


#' Heuristic length-scale bounds for a squared-exponential GP
#'
#' Computes sensible prior bounds for the GP length-scale parameter \eqn{\ell}
#' of a squared-exponential (exponentiated-quadratic) kernel,
#' \eqn{k(h) = \exp(-h^2 / (2\ell^2))}{k(h) = exp(-h^2 / (2 l^2))}.
#' The "effective range" where correlation drops to ~5\% is
#' \eqn{\ell \sqrt{2 \ln 20} \approx 2.45\,\ell}{l * sqrt(2 log(20)) ≈ 2.45 l}.
#'
#' Subsamples large datasets to avoid O(n^2) memory and time cost.
#'
#' @param coords_xy Numeric matrix of (x, y) coordinates.
#' @param q_small Numeric quantile for the lower bound. Default 0.25.
#'   Previous versions used 0.1, but the 10th percentile can be
#'   dominated by within-cluster spacing in clustered data, producing
#'   a misleadingly small lower bound.
#' @param max_n Maximum number of points to use in distance computation.
#'   Default 1000. Set to \code{Inf} to use all points.
#' @return Named numeric vector \code{c(lower, upper)} on the length-scale.
#' @examples
#' set.seed(1)
#' xy <- cbind(runif(50), runif(50))
#' gp_lengthscale_bounds(xy)
#' @export
gp_lengthscale_bounds <- function(coords_xy, q_small = 0.25, max_n = 1000L) {
  d <- .safe_dist(coords_xy, max_n = max_n)
  d <- d[d > 0]
  if (!length(d)) return(c(lower = 0.001, upper = 1))
  dmax  <- max(d)
  dq    <- stats::quantile(d, probs = q_small, names = FALSE, type = 7)
  # For SE kernel: 5% correlation at distance d  =>  l = d / sqrt(2 ln 20)
  se_range_factor <- sqrt(2 * log(20))          # ≈ 2.448
  lower <- dq   / se_range_factor
  upper <- dmax / se_range_factor
  upper <- max(upper, lower * 1.2)               # ensure separation
  c(lower = as.numeric(lower), upper = as.numeric(upper))
}

#' @rdname gp_lengthscale_bounds
#' @usage phi_prior_bounds(coords_xy, q_small = 0.25, max_n = 1000L)
#' @details \code{phi_prior_bounds} is a deprecated alias retained for
#'   backward compatibility.  Previous versions computed bounds calibrated
#'   for an exponential covariance (\code{3 / d}); the current
#'   implementation delegates to \code{gp_lengthscale_bounds} which is
#'   calibrated for the squared-exponential kernel used by
#'   \code{brms::gp()}.
#' @export
phi_prior_bounds <- function(coords_xy, q_small = 0.25, max_n = 1000L) {
  .Deprecated("gp_lengthscale_bounds")
  gp_lengthscale_bounds(coords_xy, q_small = q_small, max_n = max_n)
}

# Internal alias for backward compatibility
.phi_prior_bounds <- phi_prior_bounds


# NOTE: the previously defined .matern_rho() and .safe_ll() helpers were
# removed as dead code (no callers anywhere in the package); mvtnorm was
# dropped from Suggests along with .safe_ll().
