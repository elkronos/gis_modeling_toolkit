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
#' assumption is not. The lon/lat heuristic runs first, exactly as
#' `ensure_projected()` does: coordinates that look like degrees are taken as
#' EPSG:4326 and REPROJECTED; only when they do not is the target stamped. Both
#' branches warn.
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
  # An NA target means "leave the coordinates in whatever space they are in":
  # build_tessellation() passes it for CRS-less points that the lon/lat
  # heuristic declined to interpret, so that the grid is built in the same
  # unnamed space as the points rather than being projected on its own.
  if (inherits(crs, "crs") && is.na(crs)) return(x)
  if (is.na(sf::st_crs(x))) {
    # Stamping is the LAST resort, not the first.  ensure_projected() -- which
    # this helper's own documentation says it matches -- runs the lon/lat
    # heuristic first and REPROJECTS when the coordinates look like degrees.
    # Stamping them instead put identical input 5,400 km apart depending on
    # which of the two paths a caller happened to take, and it is how
    # summarize_by_cell(deff = "variogram") came to compare degree separations
    # against a metre range (deff = n for every cell) and how
    # assign_features_to_polygons() came to return zero rows for points that
    # build_tessellation() had just tessellated correctly.
    ll <- .looks_like_lonlat(x)
    if (isTRUE(ll$lonlat)) {
      .warn_and_log(
        paste0("%s(): `%s` has no CRS; its coordinates look like lon/lat ",
               "(xmin=%.2f, xmax=%.2f, ymin=%.2f, ymax=%.2f), so they are ",
               "taken as EPSG:4326 and reprojected to the target CRS. Set the ",
               "CRS explicitly with sf::st_crs() to suppress this."),
        caller, what, ll$bb[["xmin"]], ll$bb[["xmax"]],
        ll$bb[["ymin"]], ll$bb[["ymax"]])
      return(sf::st_transform(sf::st_set_crs(x, 4326), crs))
    }
    .warn_and_log(
      "%s(): `%s` has no CRS and its coordinates do not look like lon/lat, so it cannot be reprojected; stamping the supplied `crs` WITHOUT reprojection. Verify the coordinates are already expressed in that CRS, or set the input CRS with sf::st_crs().",
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

#' Log a warning AND raise it as an R condition
#'
#' A logger line is invisible to \code{tryCatch(warning = )},
#' \code{withCallingHandlers()}, \code{testthat::expect_warning()} and
#' \code{options(warn = 2)}.  For a situation the caller should be able to
#' catch or escalate -- data dropped, a CRS assumed, a result degraded -- the
#' log line is not enough on its own.  Used wherever the reference manual says
#' the function \emph{warns}; purely methodological cautions stay
#' \code{.log_warn()} and their documentation says "logged".
#'
#' @keywords internal
#' @noRd
.warn_and_log <- function(fmt, ...) {
  msg <- sprintf(fmt, ...)
  logger::log_warn(msg, namespace = "spatialkit")
  warning(msg, call. = FALSE)
  invisible(msg)
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
#' @param caller Name of the user-facing function, for the message prefix.
#'   Defaults to the name of the function that called this one.
#' @keywords internal
#' @noRd
.assert_sf <- function(x, what = c("POINT", "POLYGON", "MULTIPOLYGON"),
                       label = deparse(substitute(x)),
                       caller = NULL) {
  if (is.null(caller)) {
    cl <- sys.call(-1L)
    caller <- if (is.null(cl)) "spatialkit" else
      sub("^.*:::?", "", deparse(cl[[1L]])[1L])   # drop a pkg:: prefix
  }
  if (!inherits(x, "sf"))
    stop(sprintf("%s(): `%s` must be an sf object%s.", caller, label,
                 if (is.list(x) && !is.null(x$cells))
                   " (this looks like a build_tessellation() result; pass its `$cells`)"
                 else ""),
         call. = FALSE)
  gcls <- unique(as.character(sf::st_geometry_type(x, by_geometry = TRUE)))
  # `all`, not `any`: a mixed-geometry layer with one acceptable type used to
  # pass, contradicting the error text below and letting e.g. a POINT/POLYGON
  # mix through a POINT-only check.  The explicit length check keeps a zero-row
  # layer failing as it did under `any()` -- all() is vacuously TRUE on the
  # empty set.
  if (length(gcls) == 0L || !all(gcls %in% what))
    stop(sprintf("%s(): `%s` geometry must be one of: %s (found: %s).",
                 caller, label, paste(what, collapse = ", "),
                 paste(gcls, collapse = ", ")),
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
#'
#'   `n` counts the finite (y, yhat) pairs. It is **not** the number of rows
#'   `MAPE` and `SMAPE` were averaged over: both have a denominator that can be
#'   zero, and each silently drops the rows where its own denominator vanishes
#'   (`MAPE` where `y == 0`, `SMAPE` where `|y| + |yhat| == 0`), returning `NA`
#'   only when no row qualifies. That subsetting is not reported anywhere in the
#'   return value, which is why the user-facing help
#'   (see `model_metrics()`'s "Percentage errors on responses with zeros")
#'   tells callers to prefer RMSE/MAE/R2 on a response that can be zero.
#'   Reporting the per-metric row count would change this frame's column set,
#'   so it is deferred rather than done here --- see `dev/BACKLOG.md`.
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


#' Refuse arguments in `...` that nothing downstream will read
#'
#' A function whose `...` is forwarded to another, or documented as
#' "Ignored", has no terminal check anywhere: a name that matches no formal
#' vanishes without a condition.  For the evaluation functions the one
#' argument they exist to take is `newdata`, and a one-character slip
#' (`newdta = hold`) silently turned an out-of-sample RMSE of 25.24 into the
#' in-sample 1.086 with the same return shape.  `create_grid_polygons()`, which
#' has no `...`, rejects the same typo.
#'
#' @param dots `list(...)` from the caller.
#' @param fn Name of the user-facing function, for the message.
#' @param allowed Character vector of names that ARE consumed downstream.
#' @return `dots`, invisibly, when every name is allowed.
#' @keywords internal
#' @noRd
.check_dots <- function(dots, fn, allowed = character(0)) {
  if (!length(dots)) return(invisible(dots))
  nm <- names(dots)
  if (is.null(nm)) nm <- rep("", length(dots))
  unnamed <- !nzchar(nm)
  bad     <- !unnamed & !(nm %in% allowed)
  if (any(unnamed))
    stop(fn, "(): ", sum(unnamed), " unnamed argument(s) in `...` would be ",
         "ignored. Name every argument.", call. = FALSE)
  if (any(bad))
    stop(fn, "(): unused argument(s) ",
         paste(sprintf("`%s`", nm[bad]), collapse = ", "),
         if (length(allowed)) paste0(" (this method accepts ",
                                     paste(sprintf("`%s`", allowed), collapse = ", "),
                                     " through `...`)")
         else " (this method takes nothing through `...`)",
         ". A misspelt `newdata` here would silently return in-sample values.",
         call. = FALSE)
  invisible(dots)
}


#' Catch a misspelt `newdata` on the evaluation functions
#'
#' `model_metrics()`, `evaluate_insample()` and `compare_models()` forward
#' `...` to `predict()`, which checks it -- but only on the out-of-sample
#' branch.  A misspelt `newdata` leaves `newdata` NULL and the typo in `...`,
#' the in-sample branch never calls `predict()`, and the in-sample metrics
#' come back with the out-of-sample return shape.  So: arguments in `...`
#' with no `newdata` is an error here.
#'
#' @param dots `list(...)`.
#' @param newdata The caller's `newdata`.
#' @param fn User-facing function name.
#' @keywords internal
#' @noRd
.check_dots_newdata <- function(dots, newdata, fn) {
  if (is.null(newdata) && length(dots)) {
    nm <- names(dots)
    if (is.null(nm)) nm <- rep("", length(dots))
    nm <- ifelse(nzchar(nm), sprintf("`%s`", nm), "<unnamed>")
    stop(fn, "(): argument(s) ", paste(nm, collapse = ", "),
         " were supplied but `newdata` was not; nothing reads them in the ",
         "in-sample case, and the in-sample metrics would be returned as if ",
         "they were out-of-sample. Did you mean `newdata = `?", call. = FALSE)
  }
  invisible(dots)
}


#' Belsley condition index of a design matrix
#'
#' The collinearity diagnostic the literature actually thresholds: the ratio
#' of the largest to the smallest singular value of the design matrix after
#' each column has been scaled to unit Euclidean length (Belsley, Kuh & Welsch
#' 1980; for GWR, Wheeler & Tiefelsdorf 2005).  Columns are scaled but NOT
#' centred, so an intercept column stays in and a column that is constant
#' inside a window shows up as collinear with it.  The conventional threshold
#' is 30.  \code{kappa()} on the raw matrix depends on the predictors' units
#' -- a design with condition index 1322, whose local coefficients ran from
#' -86 to +150 around a true value of 2, had a raw kappa under 1e6 and raised
#' nothing.
#'
#' @param X Numeric matrix.
#' @return A non-negative number; \code{Inf} for an exactly singular design or
#'   a zero column.
#' @keywords internal
#' @noRd
.condition_index <- function(X) {
  X <- as.matrix(X)
  if (!nrow(X) || !ncol(X)) return(NA_real_)
  if (nrow(X) < ncol(X)) return(Inf)
  nrm <- sqrt(colSums(X^2))
  if (any(!is.finite(nrm)) || any(nrm == 0)) return(Inf)
  sv <- tryCatch(svd(sweep(X, 2L, nrm, "/"), nu = 0, nv = 0)$d,
                 error = function(e) NULL)
  if (is.null(sv) || !length(sv)) return(Inf)
  if (min(sv) <= .Machine$double.eps * max(sv)) return(Inf)
  max(sv) / min(sv)
}


#' Ellipsoidal (WGS84) geodesic distance between lon/lat pairs
#'
#' Vincenty's inverse formula, vectorised over pairs.  \code{sf::st_distance()}
#' on lon/lat geometry uses s2's SPHERE (R = 6371 km), and the sphere-to-WGS84
#' gap of 0.24-0.56\% is the same size as the projection distortions
#' \code{.crs_distance_error()} compares -- so the "measured error X\% vs Y\%"
#' figures were off by up to half a percentage point and the least-distorting
#' candidate was mis-ranked in 16 of 40 random wide extents.  The ellipsoidal
#' distance needs no Suggests package (\pkg{lwgeom} is not a dependency).
#'
#' @param lon1,lat1,lon2,lat2 Numeric vectors in decimal degrees, recycled.
#' @return Distances in metres.  Nearly antipodal pairs, where the iteration
#'   does not converge, fall back to the spherical distance.
#' @keywords internal
#' @noRd
.geod_distance <- function(lon1, lat1, lon2, lat2) {
  a <- 6378137; f <- 1 / 298.257223563; b <- a * (1 - f)
  rad <- pi / 180
  U1 <- atan((1 - f) * tan(lat1 * rad)); U2 <- atan((1 - f) * tan(lat2 * rad))
  L  <- (lon2 - lon1) * rad
  sU1 <- sin(U1); cU1 <- cos(U1); sU2 <- sin(U2); cU2 <- cos(U2)
  n <- max(length(U1), length(U2), length(L))
  lambda <- rep_len(L, n); U1 <- rep_len(U1, n); U2 <- rep_len(U2, n)
  sU1 <- rep_len(sU1, n); cU1 <- rep_len(cU1, n)
  sU2 <- rep_len(sU2, n); cU2 <- rep_len(cU2, n)
  L   <- rep_len(L, n)
  sinSigma <- cosSigma <- sigma <- cosSqAlpha <- cos2SigmaM <- numeric(n)
  active <- rep(TRUE, n)
  for (iter in seq_len(200L)) {
    sl <- sin(lambda[active]); cl <- cos(lambda[active])
    ss <- sqrt((cU2[active] * sl)^2 +
               (cU1[active] * sU2[active] - sU1[active] * cU2[active] * cl)^2)
    cs <- sU1[active] * sU2[active] + cU1[active] * cU2[active] * cl
    sg <- atan2(ss, cs)
    sinAlpha <- ifelse(ss == 0, 0, cU1[active] * cU2[active] * sl / ss)
    csa <- 1 - sinAlpha^2
    c2sm <- ifelse(csa == 0, 0, cs - 2 * sU1[active] * sU2[active] / csa)
    C  <- f / 16 * csa * (4 + f * (4 - 3 * csa))
    lambda_new <- L[active] + (1 - C) * f * sinAlpha *
      (sg + C * ss * (c2sm + C * cs * (-1 + 2 * c2sm^2)))
    sinSigma[active] <- ss; cosSigma[active] <- cs; sigma[active] <- sg
    cosSqAlpha[active] <- csa; cos2SigmaM[active] <- c2sm
    done <- abs(lambda_new - lambda[active]) < 1e-12
    lambda[active] <- lambda_new
    active[active] <- !done
    if (!any(active)) break
  }
  uSq <- cosSqAlpha * (a^2 - b^2) / b^2
  A <- 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
  B <- uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))
  dSigma <- B * sinSigma * (cos2SigmaM + B / 4 *
    (cosSigma * (-1 + 2 * cos2SigmaM^2) -
     B / 6 * cos2SigmaM * (-3 + 4 * sinSigma^2) * (-3 + 4 * cos2SigmaM^2)))
  d <- b * A * (sigma - dSigma)
  # Non-converged (near-antipodal) pairs: spherical fallback.
  if (any(active)) {
    la1 <- rep_len(lat1, n)[active] * rad; la2 <- rep_len(lat2, n)[active] * rad
    dl  <- L[active]
    d[active] <- 6371008.8 * acos(pmin(1, pmax(-1,
      sin(la1) * sin(la2) + cos(la1) * cos(la2) * cos(dl))))
  }
  d[rep_len(lon1, n) == rep_len(lon2, n) & rep_len(lat1, n) == rep_len(lat2, n)] <- 0
  d
}
