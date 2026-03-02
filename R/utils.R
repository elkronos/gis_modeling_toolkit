# ===========================================================================
# utils.R — Shared utilities, internal helpers, and core metric functions
#
# Migrated from 00_setup.R. All library() calls removed; dependencies
# declared in DESCRIPTION Imports.
# ===========================================================================

#' Null-coalescing operator
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

#' Check whether an sf/sfc object has a geographic (lon/lat) CRS
#'
#' @param x An sf or sfc object.
#' @return Logical; TRUE if x has a geographic CRS, FALSE otherwise.
#' @keywords internal
#' @noRd
.is_longlat <- function(x) {
  cr <- sf::st_crs(x)
  if (is.na(cr)) return(FALSE)
  isTRUE(tryCatch(sf::st_is_longlat(cr), error = function(e) FALSE))
}


#' Structured warning via logger
#' @keywords internal
#' @noRd
.log_warn <- function(fmt, ...) {
  logger::log_warn(sprintf(fmt, ...), namespace = "spatialkit")
}


#' Structured info message via logger
#' @keywords internal
#' @noRd
.log_info <- function(fmt, ...) {
  logger::log_info(sprintf(fmt, ...), namespace = "spatialkit")
}


#' Make geometries valid with best-available method
#'
#' Uses sf::st_make_valid() (sf >= 1.0) or falls back to lwgeom, then
#' st_buffer(, 0).
#'
#' @param g An sf, sfc, or sfg object.
#' @return The same object with repaired geometries.
#' @keywords internal
#' @noRd
.safe_make_valid <- function(g) {
  if ("st_make_valid" %in% getNamespaceExports("sf"))
    return(suppressWarnings(sf::st_make_valid(g)))
  if (requireNamespace("lwgeom", quietly = TRUE) &&
      "st_make_valid" %in% getNamespaceExports("lwgeom"))
    return(suppressWarnings(lwgeom::st_make_valid(g)))
  suppressWarnings(sf::st_buffer(g, 0))
}


#' Align the CRS of `a` to match `b` (transform if both defined & differ)
#'
#' @param a sf or sfc object to re-project.
#' @param b sf or sfc reference whose CRS takes precedence.
#' @return `a`, possibly transformed.
#' @keywords internal
#' @noRd
.align_crs <- function(a, b) {
  if (is.null(a) || is.null(b)) return(a)
  if (is.na(sf::st_crs(a)) || is.na(sf::st_crs(b))) return(a)
  if (sf::st_crs(a) == sf::st_crs(b)) return(a)
  sf::st_transform(a, sf::st_crs(b))
}


#' Remove duplicate points by rounded coordinate key
#'
#' Works on both sf objects with POINT geometry and bare sfc_POINT vectors.
#'
#' @param g An sf object or sfc vector with POINT geometry.
#' @param digits Rounding precision. Default 10.
#' @return `g` with duplicates removed (same class as input).
#' @keywords internal
#' @noRd
.dedup_points <- function(g, digits = 10L) {
  m <- sf::st_coordinates(g)
  if (nrow(m) == 0L) {
    if (inherits(g, "sfc")) return(g[integer(0)])
    return(g[FALSE, , drop = FALSE])
  }
  key <- paste0(round(m[, 1], digits), "_", round(m[, 2], digits))
  mask <- !duplicated(key)
  if (inherits(g, "sfc")) g[mask] else g[mask, , drop = FALSE]
}


#' Count unique points by rounded coordinate key
#'
#' @param g An sf object or sfc vector with POINT geometry.
#' @param digits Rounding precision. Default 10.
#' @return Integer count of unique coordinate positions.
#' @keywords internal
#' @noRd
.n_unique_points <- function(g, digits = 10L) {
  m <- sf::st_coordinates(g)
  if (nrow(m) == 0L) return(0L)
  key <- paste0(round(m[, 1], digits), "_", round(m[, 2], digits))
  length(unique(key))
}


#' Assert that an object is sf with one of the expected geometry types
#'
#' @param x An object.
#' @param what Character vector of acceptable geometry type names.
#' @param label Label used in error messages.
#' @keywords internal
#' @noRd
.assert_sf <- function(x, what = c("POINT", "POLYGON", "MULTIPOLYGON"),
                       label = deparse(substitute(x))) {
  if (!inherits(x, "sf"))
    stop(sprintf("Expected an sf object for `%s`.", label), call. = FALSE)
  gcls <- unique(as.character(sf::st_geometry_type(x, by_geometry = TRUE)))
  if (!any(gcls %in% what))
    stop(sprintf("`%s` geometry must be one of: %s (found: %s).",
                 label, paste(what, collapse = ", "), paste(gcls, collapse = ", ")),
         call. = FALSE)
}


#' Temporarily set the RNG seed and restore it on exit
#'
#' Returns an on.exit-compatible cleanup expression. Call inside a function:
#'   cleanup <- .with_seed(seed); on.exit(cleanup(), add = TRUE)
#'
#' @param seed Integer seed, or NULL to skip.
#' @return A zero-argument function that restores the previous RNG state.
#' @keywords internal
#' @noRd
.with_seed <- function(seed) {
  if (is.null(seed)) return(function() invisible(NULL))
  old_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (old_exists) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(seed)
  function() {
    if (old_exists)
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv))
      rm(".Random.seed", envir = .GlobalEnv)
  }
}


#' Compute basic regression metrics (RMSE, MAE, MAPE, SMAPE, R-squared, Adjusted R-squared)
#'
#' Shared across model-prep, cross-validation, and evaluation modules.
#'
#' @param y Numeric vector of observed values.
#' @param yhat Numeric vector of predicted values.
#' @param p Integer number of predictors (for Adjusted R-squared). Default NULL
#'   (Adjusted R-squared omitted).
#' @param y_train_mean Baseline mean for R-squared computation. When NULL (default),
#'   the mean of `y` is used. A scalar value (e.g., per-fold training mean) is
#'   used directly as the baseline. A per-observation vector (length matching `y`)
#'   is filtered in parallel with `y` and `yhat` to remove non-finite cases.
#' @return A data.frame with n, RMSE, MAE, MAPE, SMAPE, R2, and optionally Adj_R2.
#' @keywords internal
#' @noRd
.compute_reg_metrics <- function(y, yhat, p = NULL, y_train_mean = NULL) {
  ok <- is.finite(y) & is.finite(yhat)

  # Distinguish scalar baseline (per-fold) from per-observation vector (pooled).
  if (!is.null(y_train_mean) && length(y_train_mean) == 1L) {
    # Scalar training mean — use directly as baseline (no subsetting needed).
    y_train_mean <- y_train_mean
  } else if (!is.null(y_train_mean) && length(y_train_mean) == length(ok)) {
    # Per-observation vector — filter to match the finite-obs mask.
    y_train_mean <- y_train_mean[ok]
  }

  y <- y[ok]; yhat <- yhat[ok]
  n <- length(y)
  if (n == 0L) return(data.frame(n = 0L, RMSE = NA_real_, MAE = NA_real_,
                                 MAPE = NA_real_, SMAPE = NA_real_,
                                 R2 = NA_real_, Adj_R2 = NA_real_))
  rss  <- sum((y - yhat)^2)

  baseline <- if (!is.null(y_train_mean)) y_train_mean else mean(y)
  tss  <- sum((y - baseline)^2)
  rmse <- sqrt(rss / n)
  mae  <- mean(abs(y - yhat))

  nz <- abs(y) > .Machine$double.eps * 100
  mape <- if (any(nz)) mean(abs((y[nz] - yhat[nz]) / y[nz])) * 100 else NA_real_

  denom <- abs(y) + abs(yhat)
  smape_ok <- denom > .Machine$double.eps * 100
  smape <- if (any(smape_ok)) {
    mean(2 * abs(y[smape_ok] - yhat[smape_ok]) / denom[smape_ok]) * 100
  } else NA_real_

  r2 <- if (tss > .Machine$double.eps * n) 1 - rss / tss else NA_real_

  adj_r2 <- NA_real_
  if (!is.null(p) && is.finite(r2) && n > (p + 1L)) {
    adj_r2 <- 1 - (1 - r2) * (n - 1) / (n - p - 1)
  }

  data.frame(n = n, RMSE = rmse, MAE = mae, MAPE = mape, SMAPE = smape,
             R2 = r2, Adj_R2 = adj_r2)
}


#' Compute a distance matrix or vector with optional subsampling
#'
#' @param xy Numeric matrix of coordinates.
#' @param max_n Maximum number of rows to use. Default 1000.
#' @param seed RNG seed for reproducible subsampling. Default 42.
#' @return Numeric vector of pairwise distances.
#' @keywords internal
#' @noRd
.safe_dist <- function(xy, max_n = 1000L, seed = 42L) {
  n <- nrow(xy)
  if (n <= 1L) return(numeric(0))
  if (n > max_n) {
    cleanup <- .with_seed(seed)
    on.exit(cleanup(), add = TRUE)
    xy <- xy[sample.int(n, max_n), , drop = FALSE]
  }
  as.numeric(stats::dist(xy))
}
