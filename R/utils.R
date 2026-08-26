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

#' Resolve an object's CRS, or NULL when it has none
#'
#' `ensure_projected(target_crs = )` deliberately errors on an unusable
#' `target_crs`, so that a typo cannot silently turn it into a no-op.  Internal
#' callers, though, routinely derive the target from another object
#' (`sf::st_crs(training_data)`), and that object is legitimately allowed to
#' carry no CRS.  Passing `NA_crs_` straight through would turn a supported
#' CRS-less workflow into a hard error, so those call sites funnel the value
#' through here and get `NULL` — "no target, pick one automatically" — instead.
#'
#' @param x An sf/sfc object or anything `sf::st_crs()` accepts.
#' @return An `sf::crs` object, or NULL when the CRS is missing.
#' @keywords internal
#' @noRd
.crs_or_null <- function(x) {
  cr <- tryCatch(sf::st_crs(x), error = function(e) NA)
  if (length(cr) == 0L || all(is.na(cr))) NULL else cr
}

#' Reproject to a target CRS, or stamp it when the source has none
#'
#' `sf::st_transform()` refuses a CRS-less object ("cannot transform sfc object
#' with missing crs"), so any `if (!is.null(crs)) st_transform(x, crs)` turns a
#' `crs =` argument into a hard error for exactly the users most likely to pass
#' it — those whose data carries no CRS. Reprojection is impossible there, but
#' assumption is not: stamp the target and say so loudly, matching what
#' `ensure_projected()` already documents for the same situation.
#'
#' @param x An sf or sfc object.
#' @param crs Target CRS (anything `sf::st_crs()` accepts).
#' @param what Label naming `x` in the warning.
#' @param caller Calling function name, for the warning.
#' @return `x` reprojected to, or stamped with, `crs`.
#' @keywords internal
#' @noRd
.transform_or_stamp <- function(x, crs, what = "input", caller = "spatialkit") {
  if (is.null(crs)) return(x)
  if (is.na(sf::st_crs(x))) {
    .log_warn(
      "%s(): `%s` has no CRS, so it cannot be reprojected; stamping the supplied `crs` WITHOUT reprojection. Verify the coordinates are already expressed in that CRS, or set the input CRS with sf::st_crs().",
      caller, what)
    return(sf::st_set_crs(x, crs))
  }
  sf::st_transform(x, crs)
}

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


#' Make geometries valid
#'
#' sf >= 1.0 (a hard dependency of this package) provides st_make_valid()
#' natively, so no fallback path is needed.  (An earlier lwgeom fallback was
#' removed: modern lwgeom no longer exports st_make_valid, and the branch
#' was unreachable anyway with sf >= 1.0 installed.)
#'
#' @param g An sf, sfc, or sfg object.
#' @return The same object with repaired geometries.
#' @keywords internal
#' @noRd
.safe_make_valid <- function(g) {
  suppressWarnings(sf::st_make_valid(g))
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


#' Remove duplicate features by rounded coordinate key
#'
#' Works on sf objects and bare sfc vectors with POINT geometry, and also
#' handles MULTIPOINT (or mixed) input by keying each feature on all of its
#' vertices, keeping the deduplication mask aligned per feature.
#'
#' @param g An sf object or sfc vector with POINT/MULTIPOINT geometry.
#' @param digits Rounding precision. Default 10.
#' @return `g` with duplicates removed (same class as input).
#' @keywords internal
#' @noRd
.dedup_points <- function(g, digits = 10L) {
  geom <- if (inherits(g, "sfc")) g else sf::st_geometry(g)
  n <- length(geom)
  if (n == 0L) {
    if (inherits(g, "sfc")) return(g[integer(0)])
    return(g[FALSE, , drop = FALSE])
  }

  gtypes <- as.character(sf::st_geometry_type(geom, by_geometry = TRUE))
  if (all(gtypes == "POINT")) {
    # Fast path: one coordinate row per feature, so the key aligns 1:1.
    m <- sf::st_coordinates(geom)
    key <- paste0(round(m[, 1], digits), "_", round(m[, 2], digits))
  } else {
    # MULTIPOINT (or mixed) input: st_coordinates() returns one row per
    # *vertex*, which would misalign a per-feature mask.  Build one key
    # per feature from all of its vertices instead.
    key <- vapply(seq_len(n), function(i) {
      m <- sf::st_coordinates(geom[i])
      if (nrow(m) == 0L) return("<empty>")
      paste(round(m[, 1], digits), round(m[, 2], digits),
            sep = "_", collapse = ";")
    }, character(1))
  }

  mask <- !duplicated(key)
  if (inherits(g, "sfc")) g[mask] else g[mask, , drop = FALSE]
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
  # `all`, not `any`: a mixed-geometry layer with one acceptable type used to
  # pass, contradicting the error text below and letting e.g. a POINT/POLYGON
  # mix through a POINT-only check.  The explicit length check keeps a zero-row
  # layer failing as it did under `any()` -- all() is vacuously TRUE on the
  # empty set.
  if (length(gcls) == 0L || !all(gcls %in% what))
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
#'   Any other length is an error — recycling it would silently produce a wrong
#'   R-squared.
#' @return A data.frame with n, RMSE, MAE, MAPE, SMAPE, R2 and Adj_R2. `Adj_R2`
#'   is always present, and is `NA` when `p` is NULL or `n <= p + 1`.
#' @keywords internal
#' @noRd
.compute_reg_metrics <- function(y, yhat, p = NULL, y_train_mean = NULL) {
  ok <- is.finite(y) & is.finite(yhat)

  # A scalar baseline (per-fold training mean) is used as-is; a per-observation
  # vector (pooled) is filtered in parallel with `y`/`yhat`.  Anything else
  # would recycle against the filtered `y` and silently distort R-squared.
  if (!is.null(y_train_mean) && length(y_train_mean) != 1L) {
    if (length(y_train_mean) != length(ok))
      stop(sprintf(
        ".compute_reg_metrics(): `y_train_mean` must be length 1 (a scalar baseline) or length %d (one per observation); got %d.",
        length(ok), length(y_train_mean)
      ), call. = FALSE)
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


#' Sanitize a worker-core count to a safe positive integer
#'
#' \code{parallel::detectCores()} is documented to return \code{NA} on
#' platforms where the core count cannot be determined; arithmetic on that
#' NA (e.g. \code{max(1L, detectCores() - 1L)}) propagates it and later
#' crashes \code{if (cores > 1L)} checks or gets passed to backends.
#' This helper collapses NA/non-finite/invalid input to \code{fallback}.
#'
#' @param cores Candidate core count (any type).
#' @param fallback Integer used when `cores` is unusable. Default 1L.
#' @return A positive integer scalar.
#' @keywords internal
#' @noRd
.sanitize_core_count <- function(cores, fallback = 1L) {
  cores <- suppressWarnings(as.integer(cores[1L]))
  if (length(cores) != 1L || is.na(cores) || cores < 1L) fallback else cores
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
