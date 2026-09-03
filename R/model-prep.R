#' Prepare and sanitize an sf dataset for spatial modeling
#'
#' Ensures point geometry, projected CRS, and removes rows with missing or
#' non-finite values in modeling columns -- \emph{including} rows whose
#' geometry is empty or whose coordinates are not finite, which no model
#' backend can use.  All non-POINT geometries (including MULTIPOINT) are
#' coerced to representative points via \code{coerce_to_points()}, so
#' downstream coordinate extraction always aligns one row per observation.
#'
#' The response may not appear in \code{predictor_vars}.  Using it as its own
#' predictor is leakage no backend catches -- an out-of-bag R^2 near 1 in the
#' random forest, a silently reduced design matrix in GWR, duplicated rows and a
#' phantom \code{<none>} entry in the GWR selection table -- so it is refused
#' here.
#'
#' Column names must be syntactically valid R names (\code{make.names(x) == x}).
#' Every backend builds a model formula from these names, and a name R parses
#' as an expression -- \code{"B5-B4"} is \code{B5 - B4}, \code{"log(a)"} is a
#' function call -- would fit a different model from the one requested while
#' the fit object still recorded the name you asked for.  Rename the column
#' (for example with \code{make.names()}) before fitting.
#'
#' @param data_sf An sf object.
#' @param response_var Response variable column name.
#' @param predictor_vars Predictor column names.  May be \code{character(0)}
#'   for an intercept-only model (\code{fit_bayesian_spatial_model()} supports
#'   one; the GWR and random-forest backends do not and reject it themselves).
#' @param boundary Optional sf/sfc for CRS alignment.
#' @param pointize Strategy for non-point geometry coercion, passed to
#'   \code{\link{coerce_to_points}}.  One of \code{"auto"} (default),
#'   \code{"centroid"}, \code{"point_on_surface"}, \code{"surface"} (an alias
#'   for \code{"point_on_surface"}), \code{"line_midpoint"} or
#'   \code{"bbox_center"}.
#' @param require_response Logical; if FALSE the response column is not
#'   required to be present (useful for out-of-sample prediction where the
#'   response is unknown).  Default TRUE.
#' @return An sf object with POINT geometry, cleaned of rows carrying missing
#'   or non-finite values in the modelling columns or in the coordinates.  The CRS is projected
#'   whenever one can be established.  A CRS-less layer is decided by the
#'   lon/lat heuristic (see \code{\link{ensure_projected}}): if its bounding
#'   box fits the lon/lat envelope \emph{and} it spans more than one unit on
#'   some axis — or carries decimal-degree-like precision — it is read as
#'   EPSG:4326 and projected, with a warning.  A small planar survey inside that
#'   envelope is included in that, deliberately; only coordinates the heuristic
#'   declines are passed through as-is.  Set the CRS on \code{data_sf} if the
#'   data are planar.
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
                            pointize = c("auto", "surface", "point_on_surface",
                                         "centroid", "line_midpoint",
                                         "bbox_center"),
                            require_response = TRUE) {
  # call. = FALSE on every one of these, as elsewhere in the package: these are
  # argument-validation messages, and the call frame R would otherwise append
  # is an internal one the user did not write (prep_model_data() is reached
  # through fit_gwr_model(), the CV internals and every predict() method), so
  # printing it buries the message that names the actual problem.
  if (!inherits(data_sf, "sf"))
    stop("prep_model_data(): 'data_sf' must be an sf object.", call. = FALSE)
  pointize <- match.arg(pointize)
  if (!is.character(response_var) || length(response_var) != 1L)
    stop("prep_model_data(): 'response_var' must be a single column name.",
         call. = FALSE)
  # character(0) is allowed: an intercept-only spatial GP is a legitimate
  # model (see fit_bayesian_spatial_model()).  Backends that genuinely need a
  # predictor reject an empty set themselves, where the message can say why.
  if (!is.character(predictor_vars))
    stop("prep_model_data(): 'predictor_vars' must be a character vector.",
         call. = FALSE)

  # Every backend turns these names into a model formula -- reformulate() for
  # GWR, the model-selection sweep and the random forest's display formula,
  # sprintf() + as.formula() for the brms term.  R parses the result, so a
  # name that is not a syntactic name is not a name at all: "B5-B4" becomes
  # B5 - B4 (a DIFFERENT model, silently, with $predictor_vars still reporting
  # "B5-B4"), "log(a)" becomes a function call, and "band 4" dies inside
  # str2lang() with a parser error that names no column.  Backticking is not a
  # fix here: the sp coercion GWmodel needs runs the names through
  # make.names() anyway, so the formula and the data would disagree.  Refuse
  # early, name the columns, and say what to do about it.
  bad_nm <- setdiff(
    unique(c(if (require_response) response_var, predictor_vars)),
    make.names(unique(c(if (require_response) response_var, predictor_vars)))
  )
  if (length(bad_nm))
    stop("prep_model_data(): column name(s) ",
         paste(sprintf("'%s'", bad_nm), collapse = ", "),
         " are not syntactically valid R names, and every model backend builds",
         " a formula from them. Rename the column(s) -- e.g. names(x) <-",
         " make.names(names(x)) -- before fitting.", call. = FALSE)

  # The response among its own predictors is leakage in the random forest
  # (OOB R2 near 1, response top of the importance table), a silently reduced
  # model in GWR (model.matrix drops the duplicated column) and duplicated
  # rows plus a phantom "<none>" entry in the GWR selection table.  None of
  # the backends catches it, so catch it here.
  if (require_response && response_var %in% predictor_vars)
    stop("prep_model_data(): response '", response_var, "' is also listed in",
         " 'predictor_vars'. A model cannot use its own response as a",
         " predictor.", call. = FALSE)

  req_cols <- if (require_response) c(response_var, predictor_vars) else predictor_vars
  miss <- setdiff(req_cols, names(data_sf))
  if (length(miss))
    stop("prep_model_data(): missing required column(s): ",
         paste(miss, collapse = ", "), call. = FALSE)

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
      stop("prep_model_data(): 'boundary' must be sf/sfc when supplied.",
           call. = FALSE)
    bnd <- ensure_projected(bnd)
    data_sf <- ensure_projected(data_sf, .crs_or_null(bnd))
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

  # The attribute columns are not the whole story: st_geometry_type() reports
  # an EMPTY POINT as "POINT", so the POINT checks above (and the one in
  # fit_bayesian_spatial_model()) pass it straight through, and nothing here
  # ever looked at the coordinates themselves.  Such a row reaches GWmodel as
  # a raw sp coercion error naming no row, ranger(include_coords = TRUE) as
  # "Missing data in columns: ..x, ..y", and -- worst -- brms at predict time,
  # where one infinite coordinate makes max(Xgp) infinite, so the GP boundary
  # is infinite and EVERY basis function on the surface is NaN.  make_folds()
  # already drops these rows and says how many; do the same here.
  xy_mat  <- sf::st_coordinates(sf::st_geometry(data_sf))
  ok_geom <- rep(TRUE, nrow(data_sf))
  if (nrow(xy_mat) == nrow(data_sf) && ncol(xy_mat) >= 2L) {
    ok_geom <- is.finite(xy_mat[, 1L]) & is.finite(xy_mat[, 2L])
  } else if (nrow(data_sf) > 0L) {
    # Empty geometries drop out of st_coordinates() entirely, so the row
    # counts disagree; fall back to the per-feature test.
    ok_geom <- !sf::st_is_empty(sf::st_geometry(data_sf))
    xy_ok <- sf::st_coordinates(sf::st_geometry(data_sf)[ok_geom])
    if (nrow(xy_ok) == sum(ok_geom) && ncol(xy_ok) >= 2L)
      ok_geom[ok_geom] <- is.finite(xy_ok[, 1L]) & is.finite(xy_ok[, 2L])
  }

  keep <- ok_cc & ok_fin & ok_geom
  dropped     <- sum(!keep)
  dropped_geo <- sum(!ok_geom)
  if (dropped > 0)
    .log_warn(paste0("prep_model_data(): dropping %d row(s) with non-finite or ",
                     "missing values (%d of them with an empty or non-finite ",
                     "geometry)."), dropped, dropped_geo)

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
#' @param coords_xy Numeric matrix or data.frame of coordinates with at least
#'   two columns; the first two are used, and replicated rows are collapsed
#'   before the distance quantiles are taken — \code{brms::gp()} defaults to
#'   \code{gr = TRUE} and reduces its covariates to unique rows, so a
#'   heavily-sampled station would otherwise weight the quantile by how often it
#'   was measured rather than by where the sites are.
#' @param q_small Numeric quantile for the lower bound, a single number in
#'   \code{[0, 1]}. Default 0.25.
#'   Previous versions used 0.1, but the 10th percentile can be
#'   dominated by within-cluster spacing in clustered data, producing
#'   a misleadingly small lower bound.
#' @param max_n Maximum number of \emph{distinct} points to use in the distance
#'   computation.  Default 1000. Set to \code{Inf} to use all of them.
#' @return Named numeric vector \code{c(lower, upper)} on the length-scale;
#'   \code{c(lower = 0.001, upper = 1)} when fewer than two distinct locations
#'   or no positive distances remain.
#' @examples
#' set.seed(1)
#' xy <- cbind(runif(50), runif(50))
#' gp_lengthscale_bounds(xy)
#' @export
gp_lengthscale_bounds <- function(coords_xy, q_small = 0.25, max_n = 1000L) {
  # Without these, a vector `coords_xy` fails inside .safe_dist() with
  # "argument is of length zero" (nrow(NULL)) and an out-of-range `q_small`
  # fails inside quantile() -- neither message naming the argument at fault.
  if (!(is.matrix(coords_xy) || is.data.frame(coords_xy)))
    stop("gp_lengthscale_bounds(): 'coords_xy' must be a matrix or data.frame ",
         "of coordinates, not a ", class(coords_xy)[1L], ".", call. = FALSE)
  if (ncol(coords_xy) < 2L)
    stop("gp_lengthscale_bounds(): 'coords_xy' must have at least two columns ",
         "(x and y); it has ", ncol(coords_xy), ".", call. = FALSE)
  if (!is.numeric(q_small) || length(q_small) != 1L || !is.finite(q_small) ||
      q_small < 0 || q_small > 1)
    stop("gp_lengthscale_bounds(): 'q_small' must be a single number in [0, 1].",
         call. = FALSE)

  # "the first two are used" is a documented promise, and .safe_dist() ->
  # stats::dist() uses EVERY column.  sf::st_coordinates() on POINT Z geometry
  # returns X, Y, Z, so the documented entry point silently returned 3-D
  # bounds in mixed units (elevation metres mixed into map metres).
  coords_xy <- as.matrix(coords_xy)[, 1:2, drop = FALSE]

  # brms::gp() defaults to gr = TRUE and reduces its covariates to their
  # UNIQUE rows before doing anything else, so replicated locations (repeat
  # visits, a heavily-sampled reference station) would otherwise weight the
  # distance quantile by how often a site was measured rather than by where
  # the sites are.
  coords_xy <- coords_xy[!duplicated(coords_xy), , drop = FALSE]
  if (nrow(coords_xy) < 2L) return(c(lower = 0.001, upper = 1))

  d <- .safe_dist(coords_xy, max_n = max_n)
  d <- d[d > 0]
  if (!length(d)) return(c(lower = 0.001, upper = 1))
  dmax  <- max(d)
  dq    <- stats::quantile(d, probs = q_small, names = FALSE, type = 7)
  # For SE kernel: 5% correlation at distance d  =>  l = d / sqrt(2 ln 20)
  se_range_factor <- sqrt(2 * log(20))          # approx. 2.448
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
#'   c >= 3.2 * (ell/S),  c >= 1.25
#'   m >= 1.75 * c / (ell/S)
#' }
#' where \code{S} is \strong{the same domain measure \code{brms::gp(c = )}
#' multiplies}, namely the full pooled range of the column-centred coordinates.
#'
#' \strong{Why the range and not the half-range.}  Riutort-Mayol et al. state
#' their inequalities against the domain half-range, but both are really
#' constraints on the \emph{boundary} \eqn{L}: contain the longest plausible
#' range (\eqn{L \ge 3.2\,\ell_{upper}}) and resolve the shortest
#' (\eqn{m \ge 1.75\,L/\ell_{lower}}).  brms builds that boundary as
#' \preformatted{
#'   choose_L <- function(x, c) c * max(1, max(x) - min(x))
#' }
#' over the column-centred covariate matrix, pooled across dimensions
#' (\code{brms:::.data_gp()}) -- the FULL range, about twice the per-axis
#' half-range.  Deriving \code{c} against the half-range and handing the result
#' to \code{brms::gp()} therefore built a boundary twice as wide as intended:
#' \code{k} was sized for a boundary half the real one, so the GP was
#' systematically under-resolved, and the smallest resolvable length-scale
#' (\code{gp_ell_min}) was understated by the same factor -- making the post-fit
#' adequacy diagnostic, whose whole job is to catch under-resolution, twice too
#' lenient to fire.  Working in brms's own units removes both.
#'
#' The floor is \code{1.25} rather than Riutort-Mayol's \code{1.2} because the
#' floor is convention-dependent too: it is brms's own default (\code{c = 5/4}),
#' and on the range convention it is the more generous of the two.
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
#' @param coords_xy Numeric matrix of (already scaled) coordinates.  Only the
#'   first two columns are used, and replicated rows are collapsed, mirroring
#'   \code{brms::gp(gr = TRUE)}.
#' @param ls_bounds Named numeric \code{c(lower, upper)}, as returned by
#'   \code{gp_lengthscale_bounds()}.
#' @param k_min Integer floor on the per-dimension basis count.
#' @param max_basis Integer cap on the TOTAL basis count (\code{k^2}).  The
#'   per-dimension ceiling is derived from this as \code{floor(sqrt(max_basis))},
#'   so there is a single cap rather than two that can contradict each other.
#' @return A list with \code{k} (integer, per dimension), \code{c} (numeric),
#'   \code{S} (numeric; the pooled full range of the column-centred coordinates
#'   AFTER collapsing replicated rows, i.e. exactly what \code{brms::gp(c = )}
#'   multiplies under its default \code{gr = TRUE}) and \code{capped}
#'   (logical).
#' @keywords internal
#' @noRd
.gp_basis_spec <- function(coords_xy, ls_bounds,
                           k_min = 10L, max_basis = 2500L) {
  # Reproduce brms::choose_L()'s domain measure exactly: centre each column,
  # then take the range over the POOLED matrix.  na.rm mirrors brms.
  xy <- as.matrix(coords_xy)[, 1:2, drop = FALSE]
  # brms::gp() defaults to gr = TRUE, and brms:::.data_gp() reduces Xgp to its
  # UNIQUE rows BEFORE taking cmeans and the range.  Centring and ranging over
  # all rows therefore disagreed with the boundary brms actually builds
  # whenever two observations share a location: with one heavily-sampled
  # station the boundary came out at 0.68 of the intended one, and everything
  # derived from S -- gp_c, and the gp_ell_min threshold the post-fit adequacy
  # warning compares posterior length-scale draws against -- was wrong by the
  # same factor.  Mirror the reduction.
  xy <- xy[!duplicated(xy), , drop = FALSE]
  if (!nrow(xy)) xy <- matrix(0, nrow = 1L, ncol = 2L)
  Xc <- sweep(xy, 2L, colMeans(xy, na.rm = TRUE))
  S  <- suppressWarnings(
    max(1, max(Xc, na.rm = TRUE) - min(Xc, na.rm = TRUE)))
  if (!is.finite(S) || S <= 0) S <- 1

  # Ratios of length-scale to the domain measure brms will use.
  r_lo <- max(ls_bounds[["lower"]] / S, .Machine$double.eps)  # must resolve
  r_hi <- max(ls_bounds[["upper"]] / S, r_lo)                 # must contain

  # 1.25, not 1.2: the floor is stated on the half-range convention by
  # Riutort-Mayol et al., and 5/4 is brms's own default on this one.
  c_val <- max(3.2 * r_hi, 1.25)
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
