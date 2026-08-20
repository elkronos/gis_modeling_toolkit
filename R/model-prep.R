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


#' Choose GP basis count and boundary factor from the length-scale/domain ratio
#'
#' Implements the practical recommendations of Riutort-Mayol et al. (2023,
#' Statistics and Computing 33:1) for the squared-exponential kernel used by
#' \code{brms::gp()}:
#' \preformatted{
#'   c >= 3.2 * (ell/S),  c >= 1.2
#'   m >= 1.75 * c / (ell/S)
#' }
#' where \code{S} is the half-range of the (scaled) coordinate domain.
#'
#' \code{m} is the count PER DIMENSION.  brms expands a full tensor grid over
#' the GP covariates, so the fitted model carries \code{m^D} basis functions
#' (\code{D = 2} for \code{gp(..x, ..y)}) -- see brms/R/data-predictor.R,
#' \code{data_gp()}:  \code{out[[paste0("NBgp", pi)]] <- k ^ D}.  Because
#' \code{brms::gp()} accepts a single \code{k} shared across dimensions, the
#' per-dimension recommendations are collapsed into one value.
#'
#' The boundary factor is set from the UPPER length-scale bound -- it must
#' contain the longest plausible correlation range -- while the basis count is
#' set from the LOWER bound, because it must resolve the shortest.  Deriving
#' both from the same length-scale estimate is the point: the previous rule
#' scaled the basis count with the number of observations, which made the
#' approximation cost grow as n while adding no resolution the data supported.
#'
#' @param coords_xy Numeric matrix of (already scaled) coordinates.
#' @param ls_bounds Named numeric \code{c(lower, upper)}, as returned by
#'   \code{gp_lengthscale_bounds()}.
#' @param k_min Integer floor on the per-dimension basis count.
#' @param max_basis Integer cap on the TOTAL basis count (\code{k^2}).  The
#'   per-dimension ceiling is derived from this as \code{floor(sqrt(max_basis))},
#'   so there is a single cap rather than two that can contradict each other.
#' @return A list with \code{k} (integer, per dimension), \code{c} (numeric),
#'   \code{S} (numeric half-range) and \code{capped} (logical).
#' @keywords internal
#' @noRd
.gp_basis_spec <- function(coords_xy, ls_bounds,
                           k_min = 10L, max_basis = 2500L) {
  S <- max(apply(coords_xy, 2, function(z) diff(range(z)) / 2))
  if (!is.finite(S) || S <= 0) S <- 1

  # Ratios of length-scale to domain half-range.
  r_lo <- max(ls_bounds[["lower"]] / S, .Machine$double.eps)  # must resolve
  r_hi <- max(ls_bounds[["upper"]] / S, r_lo)                 # must contain

  c_val <- max(3.2 * r_hi, 1.2)
  k_raw <- ceiling(1.75 * c_val / r_lo)

  k_max  <- as.integer(floor(sqrt(max_basis)))
  k_val  <- min(max(k_raw, k_min), k_max)
  capped <- k_raw > k_max

  list(k = as.integer(k_val), c = as.numeric(c_val),
       S = as.numeric(S), capped = capped)
}


# NOTE: the previously defined .matern_rho() and .safe_ll() helpers were
# removed as dead code (no callers anywhere in the package); mvtnorm was
# dropped from Suggests along with .safe_ll().
