# -----------------------------------------------------------------------------
# Internal: shared CV helpers
# -----------------------------------------------------------------------------

#' Safe tibble/data.frame constructor
#' @keywords internal
#' @noRd
.safe_tibble <- function(...) {
  if (requireNamespace("tibble", quietly = TRUE)) {
    tibble::tibble(...)
  } else {
    data.frame(..., stringsAsFactors = FALSE)
  }
}


#' Re-map pre-built fold indices to the subset of rows that survive prep
#'
#' @param folds List of list(train, test) using original row IDs.
#' @param keep_idx Integer vector of surviving row IDs.
#' @param k Fallback fold count if folds is NULL.
#' @param seed RNG seed for random fold creation.
#' @return List of list(train, test) with IDs restricted to keep_idx.
#' @keywords internal
#' @noRd
.remap_folds <- function(folds, keep_idx, k = 5L, seed = 123L) {
  if (is.null(folds)) {
    .log_warn(
      ".remap_folds(): no fold specification provided; falling back to random k-fold CV (k=%d). Random folds leak spatial autocorrelation and overstate out-of-sample performance.",
      k
    )
    warning(
      ".remap_folds(): falling back to random k-fold CV. For spatial data, use make_folds(method='block_kfold') to avoid optimistic performance estimates.",
      call. = FALSE
    )
    cleanup <- .with_seed(seed)
    on.exit(cleanup(), add = TRUE)
    fold_id <- sample(rep(seq_len(k), length.out = length(keep_idx)))
    out <- vector("list", k)
    for (i in seq_len(k))
      out[[i]] <- list(train = keep_idx[fold_id != i], test = keep_idx[fold_id == i])
    return(out)
  }

  if (!is.null(folds$folds) && is.list(folds$folds)) {
    folds <- folds$folds
  }

  remapped <- lapply(folds, function(f) {
    list(
      train = keep_idx[stats::na.omit(match(f$train, keep_idx))],
      test  = keep_idx[stats::na.omit(match(f$test, keep_idx))]
    )
  })

  empty_folds <- vapply(remapped, function(f) length(f$test) == 0L, logical(1))
  if (any(empty_folds)) {
    .log_warn(".remap_folds(): %d fold(s) have empty test sets after remapping.",
              sum(empty_folds))
  }

  remapped[!empty_folds]
}


#' Compute overall metrics from a predictions data.frame
#'
#' When the predictions data.frame contains a \code{y_train_mean} column
#' (one value per observation, equal to the training-fold mean for the fold
#' that produced that prediction), R² is computed against the correct
#' out-of-sample baseline: each observation's contribution to TSS uses the
#' training mean from its fold, not the pooled test-set mean.
#'
#' Adjusted R² is intentionally reported as \code{NA} for pooled CV metrics
#' because there are no well-defined degrees of freedom when predictions are
#' aggregated across independently fitted folds.  Per-fold Adj_R² is also
#' set to \code{NA} for GWR (whose effective parameter count is not the
#' global predictor count p) and for the Bayesian GP model (complex
#' effective degrees of freedom).  Adj_R² is only meaningful for models
#' with a fixed, global number of coefficients.
#'
#' @keywords internal
#' @noRd
.cv_overall_metrics <- function(preds) {
  ok <- is.finite(preds$y) & is.finite(preds$yhat)
  if (!any(ok))
    return(data.frame(RMSE = NA_real_, MAE = NA_real_, MAPE = NA_real_,
                      SMAPE = NA_real_, R2 = NA_real_, Adj_R2 = NA_real_,
                      n_pred = 0L))

  # Pass per-observation training-fold means when available so that
  # .compute_reg_metrics() uses the correct out-of-sample R² baseline.
  ytm <- if ("y_train_mean" %in% names(preds)) preds$y_train_mean else NULL

  # Intentionally omit `p` here (pass NULL) so that pooled Adj_R² is NA.

  # Adjusted R² requires well-defined degrees of freedom (n, p) from a
  # single model fit.  When predictions are pooled across CV folds, `n` is
  # the total number of held-out observations while each fold was fit
  # independently with its own training set — the resulting Adj_R² has no
  # clean statistical interpretation.  Per-fold Adj_R² (computed in
  # .cv_run_folds()) is valid and should be used instead.
  met <- .compute_reg_metrics(preds$y, preds$yhat, p = NULL, y_train_mean = ytm)
  data.frame(RMSE = met$RMSE, MAE = met$MAE, MAPE = met$MAPE, SMAPE = met$SMAPE,
             R2 = met$R2, Adj_R2 = met$Adj_R2, n_pred = met$n,
             stringsAsFactors = FALSE)
}


# -----------------------------------------------------------------------------
# Shared fold runner
# -----------------------------------------------------------------------------

#' Fit-predict a single CV fold
#'
#' Encapsulates the per-fold work so it can be called sequentially or in
#' parallel.  Returns \code{NULL} on failure so the caller can filter.
#'
#' @keywords internal
#' @noRd
.cv_fit_one_fold <- function(i, dat_sf, response_var, remapped_fold,
                             keep_idx, fit_one, fold_info_fn, predict_args,
                             p) {
  tr_pos <- stats::na.omit(match(remapped_fold$train, keep_idx))
  te_pos <- stats::na.omit(match(remapped_fold$test, keep_idx))
  if (length(tr_pos) < 2L || length(te_pos) < 1L) return(NULL)

  train_sf <- dat_sf[tr_pos, , drop = FALSE]
  test_sf  <- dat_sf[te_pos, , drop = FALSE]

  # Drop geometry once per fold to avoid redundant copies on wide data frames
  train_df <- sf::st_drop_geometry(train_sf)
  test_df  <- sf::st_drop_geometry(test_sf)

  # Fit model on training fold
  fit_obj <- try(fit_one(train_sf), silent = TRUE)
  if (inherits(fit_obj, "try-error")) {
    .log_warn(".cv_run_folds(): fold %d fit failed; skipping.", i)
    return(NULL)
  }
  if (!inherits(fit_obj, "spatial_fit")) {
    .log_warn(".cv_run_folds(): fold %d did not return a spatial_fit; skipping.", i)
    return(NULL)
  }

  # Predict on test fold via the S3 generic
  y_true <- test_df[[response_var]]
  y_hat  <- try(
    do.call(predict, c(list(object = fit_obj, newdata = test_sf), predict_args)),
    silent = TRUE
  )
  if (inherits(y_hat, "try-error") || !is.numeric(y_hat)) {
    .log_warn(".cv_run_folds(): fold %d predict failed; skipping.", i)
    return(NULL)
  }

  # Training-set mean: the correct null-model baseline for out-of-sample
  # R².  Using the test-set mean instead would give the null model credit
  # for knowing information that was not available at prediction time,
  # systematically inflating CV R².
  y_train <- train_df[[response_var]]
  y_train_mean <- mean(y_train[is.finite(y_train)], na.rm = TRUE)

  met <- .compute_reg_metrics(y_true, y_hat, p = p,
                              y_train_mean = y_train_mean)
  if (met$n == 0L) return(NULL)

  # Base fold stats
  fs <- data.frame(
    fold = i, n_train = length(tr_pos), n_test = length(y_true),
    n_pred = met$n,
    RMSE = met$RMSE, MAE = met$MAE, MAPE = met$MAPE, SMAPE = met$SMAPE,
    R2 = met$R2, Adj_R2 = met$Adj_R2, stringsAsFactors = FALSE
  )

  # Append model-specific per-fold info (bandwidth, gp_k, CRPS, coverage …)
  if (!is.null(fold_info_fn)) {
    extra <- try(fold_info_fn(fit_obj, test_sf, y_true, y_hat), silent = TRUE)
    if (!inherits(extra, "try-error") && is.list(extra)) {
      for (cn in names(extra)) fs[[cn]] <- extra[[cn]]
    }
  }

  # Prediction rows — include the training-fold mean so that
  # .cv_overall_metrics() can compute the pooled R² with the
  # correct per-observation baseline.
  pr <- data.frame(
    `..row_id` = test_sf$`..row_id`, fold = i,
    y = as.numeric(y_true), yhat = as.numeric(y_hat),
    y_train_mean = y_train_mean,
    stringsAsFactors = FALSE
  )

  list(pred_row = pr, fold_stat = fs)
}


#' Resolve parallel settings into a usable core count
#'
#' Returns 1L for sequential execution or an integer >= 2 for parallel.
#' On Windows, \code{parallel::mclapply()} falls back to serial, so we
#' warn and return 1L.
#'
#' @param parallel Logical or positive integer.  \code{TRUE} auto-detects,
#'   an integer sets the core count explicitly, \code{FALSE} is sequential.
#' @param n_cores Deprecated alias kept for backwards compatibility.
#'   Overrides \code{parallel} when not \code{NULL}.
#' @return Integer number of worker cores to use (1 = sequential).
#' @keywords internal
#' @noRd
.resolve_n_cores <- function(parallel = FALSE, n_cores = NULL) {
  if (!is.null(n_cores)) {
    n_cores <- as.integer(n_cores)
    if (n_cores < 1L) n_cores <- 1L
    eff <- n_cores
  } else if (isTRUE(parallel)) {
    eff <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
  } else if (is.numeric(parallel) && length(parallel) == 1L && parallel > 1) {
    eff <- as.integer(parallel)
  } else {
    return(1L)
  }
  if (.Platform$OS.type == "windows" && eff > 1L) {
    message("cv parallel: forked parallelism (mclapply) is not available on Windows. ",
            "Falling back to sequential execution. For Windows parallelism, consider ",
            "the 'future' / 'future.apply' packages.")
    return(1L)
  }
  eff
}


#' Run a fit-predict loop across CV folds
#'
#' Fits a model on each training fold and predicts onto the test fold
#' using the standardised \code{predict()} S3 method.  This is the
#' single implementation that both \code{cv_gwr()} and \code{cv_bayes()}
#' delegate to.
#'
#' When \code{parallel = TRUE} (or an integer > 1), folds are fitted in
#' parallel using \code{parallel::mclapply()}, which yields near-linear
#' speedup on macOS and Linux.  On Windows, forked parallelism is not
#' available and execution falls back to sequential with a message.
#'
#' @param dat_sf Prepared sf data (projected, clean).
#' @param response_var Character(1).
#' @param predictor_vars Character vector.
#' @param remapped_folds List of list(train, test).
#' @param keep_idx Integer vector of surviving row IDs.
#' @param fit_one A function(train_sf, ...) that returns a \code{spatial_fit}.
#' @param fold_info_fn Optional function(spatial_fit) -> named list of extra
#'   per-fold metadata columns (e.g. bandwidth, gp_k).
#' @param predict_args Named list of extra arguments for predict().
#' @param p Integer number of predictors for Adj R² (NULL to skip).
#'   Only meaningful when the model uses a fixed, global set of p
#'   coefficients (e.g. a linear model).  For models with spatially
#'   varying coefficients (GWR) or complex effective degrees of freedom
#'   (GP-based models), pass NULL so that per-fold Adj_R² is reported as
#'   NA rather than a misleadingly favourable value.
#' @param parallel Logical or positive integer.  If \code{TRUE},
#'   auto-detect the number of cores; if an integer > 1, use that many
#'   cores; if \code{FALSE} (default), run sequentially.
#' @param n_cores \emph{Deprecated.}
#'   Explicit core count; overrides \code{parallel} when set.
#' @return List with pred_rows and fold_stats.
#' @keywords internal
#' @noRd
.cv_run_folds <- function(dat_sf, response_var, predictor_vars,
                          remapped_folds, keep_idx, fit_one,
                          fold_info_fn = NULL, predict_args = list(),
                          p = NULL, parallel = FALSE, n_cores = NULL) {
  cores <- .resolve_n_cores(parallel, n_cores)

  fold_worker <- function(i) {
    .cv_fit_one_fold(
      i = i, dat_sf = dat_sf, response_var = response_var,
      remapped_fold = remapped_folds[[i]], keep_idx = keep_idx,
      fit_one = fit_one, fold_info_fn = fold_info_fn,
      predict_args = predict_args, p = p
    )
  }

  if (cores > 1L) {
    message(sprintf("cv: running %d folds in parallel on %d cores.",
                    length(remapped_folds), cores))
    results <- parallel::mclapply(
      seq_along(remapped_folds), fold_worker, mc.cores = cores
    )
  } else {
    results <- lapply(seq_along(remapped_folds), fold_worker)
  }

  # Filter NULLs (failed / skipped folds) and unpack
  results <- Filter(Negate(is.null), results)
  pred_rows  <- lapply(results, `[[`, "pred_row")
  fold_stats <- lapply(results, `[[`, "fold_stat")

  list(pred_rows = pred_rows, fold_stats = fold_stats)
}


# -----------------------------------------------------------------------------
# Spatial autocorrelation range estimation
# -----------------------------------------------------------------------------

#' Estimate the spatial autocorrelation range from data
#'
#' Fits exponential (or spherical) variogram models and returns the
#' \emph{effective range} — the distance at which the semivariance reaches
#' ~95 \% of the sill.
#'
#' To guard against anisotropy, the function first estimates directional
#' variograms at 0\u00b0 (N–S) and 90\u00b0 (E–W) azimuths (tolerance 22.5\u00b0,
#' which avoids double-counting point pairs near the 45\u00b0 diagonal but
#' requires denser point clouds for stable estimates).
#' When both fits succeed the \strong{maximum} of the two directional ranges
#' is returned, which is the conservative choice for spatial block CV —
#' blocks must be at least as large as the longest autocorrelation range to
#' avoid information leakage.
#'
#' If either directional fit fails (e.g., too few point pairs in a direction),
#' the function falls back to an omnidirectional (isotropic) variogram.
#'
#' A log warning is emitted when notable anisotropy is detected (ratio of
#' directional ranges > 1.5).
#'
#' The returned range is in the coordinate units of the (projected) data and
#' can be passed directly to \code{make_folds(block_size = ...)} to ensure
#' that CV blocks are at least as wide as the autocorrelation range.
#'
#' @param points_sf An sf object with point geometries (will be projected
#'   automatically if in geographic CRS).
#' @param response_var Character(1) name of the response column.
#' @param predictor_vars Optional character vector.  When supplied, an OLS
#'   residual variogram is fitted instead of a raw-response variogram, which
#'   better reflects the autocorrelation that the spatial model must handle.
#' @param n_max Maximum number of points to subsample before fitting.
#'   Variogram estimation is O(n²) so this keeps runtime bounded.
#' @param cutoff Fraction of the maximum inter-point distance to use as
#'   the variogram lag cutoff.  Default 0.5.
#' @param seed Optional RNG seed for subsampling reproducibility.
#' @return A single positive numeric value (the effective range in projected
#'   coordinate units), or \code{NA_real_} if estimation fails.
#' @export
estimate_sac_range <- function(points_sf, response_var,
                               predictor_vars = NULL,
                               n_max = 5000L, cutoff = 0.5,
                               seed = NULL) {
  if (!requireNamespace("gstat", quietly = TRUE)) {
    .log_warn("estimate_sac_range(): package 'gstat' is required for variogram estimation; returning NA.")
    return(NA_real_)
  }
  if (!inherits(points_sf, "sf"))
    stop("estimate_sac_range(): `points_sf` must be an sf object.", call. = FALSE)
  if (!all(sf::st_geometry_type(points_sf, by_geometry = TRUE) %in%
           c("POINT", "MULTIPOINT")))
    points_sf <- coerce_to_points(points_sf, "auto")

  pts <- ensure_projected(points_sf)

  # Subsample if large
  n <- nrow(pts)
  if (n > n_max) {
    cleanup <- .with_seed(seed)
    on.exit(cleanup(), add = TRUE)
    idx <- sample.int(n, n_max)
    pts <- pts[idx, , drop = FALSE]
    n <- n_max
  }

  if (n < 30L) {
    .log_warn("estimate_sac_range(): fewer than 30 points; variogram estimate unreliable. Returning NA.")
    return(NA_real_)
  }

  # Build the variable to model: raw response or OLS residuals
  y <- sf::st_drop_geometry(pts)[[response_var]]
  if (!is.null(predictor_vars) && length(predictor_vars) > 0L) {
    df <- sf::st_drop_geometry(pts)
    ok_preds <- intersect(predictor_vars, names(df))
    if (length(ok_preds) > 0L) {
      fml <- stats::reformulate(ok_preds, response_var)
      lm_fit <- try(stats::lm(fml, data = df, na.action = stats::na.exclude),
                     silent = TRUE)
      if (!inherits(lm_fit, "try-error")) {
        resid <- stats::residuals(lm_fit)
        if (length(resid) != nrow(pts)) {
          .log_warn("estimate_sac_range(): OLS residual length (%d) does not match data rows (%d); using raw response.",
                    length(resid), nrow(pts))
        } else {
          y <- resid
        }
      }
    }
  }

  pts$..sac_var <- as.numeric(y)
  pts <- pts[is.finite(pts$..sac_var), , drop = FALSE]
  if (nrow(pts) < 30L) {
    .log_warn("estimate_sac_range(): too few finite values after filtering; returning NA.")
    return(NA_real_)
  }

  # Empirical variogram
  max_dist <- try({
    bb <- sf::st_bbox(pts)
    sqrt((bb["xmax"] - bb["xmin"])^2 + (bb["ymax"] - bb["ymin"])^2)
  }, silent = TRUE)

  if (inherits(max_dist, "try-error") || !is.finite(max_dist) || max_dist <= 0)
    return(NA_real_)

  cutoff_dist <- as.numeric(cutoff * max_dist)

  # --- Helper: fit a variogram model and return the effective range --------
  .fit_vgm_range <- function(vg) {
    if (inherits(vg, "try-error") || !inherits(vg, "data.frame") || NROW(vg) < 3L) return(NA_real_)
    vgm_model <- try(
      gstat::fit.variogram(vg, gstat::vgm(model = "Exp")),
      silent = TRUE
    )
    if (inherits(vgm_model, "try-error")) {
      vgm_model <- try(
        gstat::fit.variogram(vg, gstat::vgm(model = "Sph")),
        silent = TRUE
      )
    }
    if (inherits(vgm_model, "try-error") || !is.data.frame(vgm_model))
      return(NA_real_)
    spatial_rows <- vgm_model[vgm_model$model != "Nug", , drop = FALSE]
    if (nrow(spatial_rows) == 0L || all(!is.finite(spatial_rows$range)))
      return(NA_real_)
    raw_range  <- max(spatial_rows$range, na.rm = TRUE)
    model_type <- spatial_rows$model[which.max(spatial_rows$range)]
    eff <- if (identical(as.character(model_type), "Exp")) 3 * raw_range else raw_range
    if (!is.finite(eff) || eff <= 0) NA_real_ else eff
  }

  # --- Directional variograms (0° and 90°, tolerance 22.5°) ----------------
  # gstat uses azimuth in degrees clockwise from north.  0° = N-S, 90° = E-W.
  dir_ranges <- vapply(c(0, 90), function(az) {
    vg_dir <- try(
      gstat::variogram(..sac_var ~ 1, data = pts,
                       cutoff = cutoff_dist,
                       alpha = az, tol.hor = 22.5),
      silent = TRUE
    )
    .fit_vgm_range(vg_dir)
  }, numeric(1))

  dir_success <- sum(is.finite(dir_ranges)) == 2L

  # --- Select the effective range ------------------------------------------
  if (dir_success) {
    effective_range <- max(dir_ranges)
    ratio <- max(dir_ranges) / min(dir_ranges)
    if (ratio > 1.5) {
      .log_warn(
        "estimate_sac_range(): notable anisotropy detected (range ratio %.1f). Directional ranges: 0\u00b0 = %.1f, 90\u00b0 = %.1f. Using the maximum.",
        ratio, dir_ranges[1], dir_ranges[2]
      )
    }
  } else {
    # --- Isotropic variogram (fallback when directional fits fail) ----------
    vg_iso <- try(
      gstat::variogram(..sac_var ~ 1, data = pts, cutoff = cutoff_dist),
      silent = TRUE
    )
    iso_range <- .fit_vgm_range(vg_iso)
    if (is.finite(iso_range)) {
      effective_range <- iso_range
    } else {
      # Neither directional nor isotropic succeeded
      .log_warn("estimate_sac_range(): variogram model fit failed; returning NA.")
      return(NA_real_)
    }
  }

  if (!is.finite(effective_range) || effective_range <= 0) {
    .log_warn("estimate_sac_range(): estimated range is non-positive or non-finite; returning NA.")
    return(NA_real_)
  }

  effective_range
}


#' Compute grid dimensions that respect a minimum block size
#'
#' Given a bounding box and a minimum block edge length, returns nx/ny
#' values such that each cell is at least \code{block_size} wide and tall.
#'
#' @param bb An sf bbox.
#' @param block_size Positive numeric minimum block edge length (CRS units).
#' @param k Integer fold count (floor for nx * ny).
#' @return Named list with \code{nx} and \code{ny}.
#' @keywords internal
#' @noRd
.block_dims_from_size <- function(bb, block_size) {
  w <- as.numeric(bb["xmax"] - bb["xmin"])
  h <- as.numeric(bb["ymax"] - bb["ymin"])
  nx <- max(1L, floor(w / block_size))
  ny <- max(1L, floor(h / block_size))
  list(nx = nx, ny = ny)
}


# -----------------------------------------------------------------------------
# Fold Construction
# -----------------------------------------------------------------------------

#' Create spatial cross-validation folds
#'
#' Builds train/test splits using random K-fold, spatial block K-fold, or
#' buffered leave-one-out strategies.
#'
#' For \code{block_kfold}, the default grid sizing is purely geometric and
#' unrelated to the autocorrelation range of the data.  When blocks are
#' smaller than the autocorrelation range, spatially correlated observations
#' leak across folds and CV metrics become optimistic.  Use
#' \code{block_size} to set a minimum block edge length (in CRS units), or
#' set \code{auto_range = TRUE} to estimate the range from an empirical
#' variogram and enforce it automatically.
#'
#' @param points_sf An sf object.
#' @param k Integer; number of folds.
#' @param method One of "random_kfold", "block_kfold", "buffered_loo".
#' @param seed Optional integer RNG seed.
#' @param block_nx,block_ny Optional grid dimensions for block_kfold.
#'   Ignored when \code{block_size} or \code{auto_range} override them.
#' @param block_multiplier Numeric; target blocks multiplier. Default 3.
#' @param block_size Optional positive numeric minimum block edge length
#'   (in projected CRS units).  When supplied, grid dimensions are clamped
#'   so that every block is at least this wide and tall.  Takes precedence
#'   over \code{block_nx}/\code{block_ny} and \code{block_multiplier}.
#' @param auto_range Logical.  If \code{TRUE}, the spatial autocorrelation
#'   range is estimated via \code{estimate_sac_range()} — which fits
#'   directional variograms to account for anisotropy — and used as the
#'   minimum \code{block_size}.  Requires \code{response_var}.  An explicit
#'   \code{block_size} takes precedence.  Default \code{FALSE}.
#' @param response_var Character(1) response column name.  Required when
#'   \code{auto_range = TRUE}.
#' @param predictor_vars Optional character vector of predictor column names.
#'   Passed to \code{estimate_sac_range()} for residual variogram estimation.
#' @param boundary Optional polygonal sf/sfc for block_kfold.
#' @param buffer Positive numeric distance for buffered_loo.
#' @param drop_empty_blocks Logical. Default TRUE.
#' @return A list with method, k, folds, assignment, params.
#' @export
make_folds <- function(points_sf, k,
                       method = c("random_kfold", "block_kfold", "buffered_loo"),
                       seed = NULL, block_nx = NULL, block_ny = NULL,
                       block_multiplier = 3, block_size = NULL,
                       auto_range = FALSE, response_var = NULL,
                       predictor_vars = NULL, boundary = NULL,
                       buffer = NULL,
                       drop_empty_blocks = TRUE) {
  method <- match.arg(method)

  cleanup <- .with_seed(seed)
  on.exit(cleanup(), add = TRUE)

  if (!inherits(points_sf, "sf")) stop("make_folds(): `points_sf` must be an sf object.")
  if (!all(sf::st_geometry_type(points_sf, by_geometry = TRUE) %in%
           c("POINT", "MULTIPOINT")))
    points_sf <- coerce_to_points(points_sf, "auto")
  if (!("..row_id" %in% names(points_sf)))
    points_sf$..row_id <- seq_len(nrow(points_sf))

  .ret <- function(method, k, folds, assignment, params)
    list(method = method, k = k, folds = folds, assignment = assignment, params = params)

  # ---- RANDOM K-FOLD ----
  if (method == "random_kfold") {
    n <- nrow(points_sf)
    if (k < 2) k <- 2L
    if (k > n) { .log_warn("make_folds(random_kfold): k > n; reducing."); k <- n }
    idx <- sample.int(n, n)
    sizes <- rep(floor(n / k), k)
    remainder <- n - sum(sizes)
    if (remainder > 0) sizes[seq_len(remainder)] <- sizes[seq_len(remainder)] + 1L
    splits <- vector("list", k); start <- 1L; assign_vec <- integer(n)
    for (j in seq_len(k)) {
      stop_i <- start + sizes[j] - 1L
      test_idx <- idx[start:stop_i]; train_idx <- setdiff(idx, test_idx)
      splits[[j]] <- list(train = train_idx, test = test_idx)
      assign_vec[test_idx] <- j; start <- stop_i + 1L
    }
    return(.ret(method, k, splits,
                .safe_tibble(row_id = points_sf$..row_id, fold = assign_vec),
                list(seed = seed)))
  }

  # ---- BLOCK K-FOLD ----
  if (method == "block_kfold") {
    if (k < 2) k <- 2L
    pts <- ensure_projected(points_sf)
    reg <- if (!is.null(boundary)) {
      b <- ensure_projected(boundary, sf::st_crs(pts))
      if (inherits(b, "sfc")) b <- sf::st_sf(geometry = b)
      b <- .safe_make_valid(sf::st_union(b))
      mat <- sf::st_intersects(pts, b, sparse = FALSE)
      inside_any <- apply(mat, 1L, any)
      if (!all(inside_any)) {
        bb_pts <- sf::st_as_sfc(sf::st_bbox(pts)) |> sf::st_set_crs(sf::st_crs(pts))
        b <- suppressWarnings(
          sf::st_union(.safe_make_valid(sf::st_sf(geometry = c(b, bb_pts)))))
      }
      b
    } else {
      sf::st_as_sfc(sf::st_bbox(pts)) |> sf::st_set_crs(sf::st_crs(pts)) |> sf::st_sf()
    }

    # --- Autocorrelation-aware block sizing ---
    sac_range <- NA_real_
    if (isTRUE(auto_range) && !is.null(response_var)) {
      sac_range <- estimate_sac_range(pts, response_var = response_var,
                                      predictor_vars = predictor_vars,
                                      seed = seed)
      if (is.finite(sac_range) && sac_range > 0) {
        message(sprintf(
          "make_folds(block_kfold): estimated spatial autocorrelation range = %.1f CRS units; using as minimum block size.",
          sac_range
        ))
        # auto_range sets block_size only if the caller didn't supply one
        if (is.null(block_size)) {
          block_size <- sac_range
        } else if (block_size < sac_range) {
          .log_warn(
            "make_folds(block_kfold): supplied block_size (%.1f) is smaller than the estimated autocorrelation range (%.1f). Spatial CV may still leak correlated information.",
            block_size, sac_range
          )
          warning(
            sprintf("make_folds(): block_size (%.1f) < estimated autocorrelation range (%.1f). Consider increasing block_size to reduce information leakage across folds.",
                    block_size, sac_range),
            call. = FALSE
          )
        }
      } else {
        .log_warn("make_folds(block_kfold): auto_range requested but estimation returned NA; falling back to geometric blocks.")
      }
    } else if (isTRUE(auto_range) && is.null(response_var)) {
      .log_warn("make_folds(block_kfold): auto_range = TRUE but response_var is NULL; cannot estimate range. Falling back to geometric blocks.")
      warning("make_folds(): auto_range requires response_var; ignoring.", call. = FALSE)
    }

    bb <- sf::st_bbox(reg)

    # Determine grid dimensions: block_size constrains nx/ny
    if (!is.null(block_size) && is.numeric(block_size) && block_size > 0) {
      size_dims <- .block_dims_from_size(bb, block_size)
      # If the caller also supplied explicit block_nx/block_ny, warn about override
      if (!is.null(block_nx) || !is.null(block_ny)) {
        .log_warn(
          "make_folds(block_kfold): block_size (%.1f) overrides explicit block_nx/block_ny.",
          block_size
        )
      }
      nx <- size_dims$nx
      ny <- size_dims$ny

      # Ensure at least k blocks so each fold can get one
      if (nx * ny < k) {
        .log_warn(
          "make_folds(block_kfold): block_size produces only %d blocks (< k = %d). Reducing k to match.",
          nx * ny, k
        )
        k <- max(2L, nx * ny)
      }
    } else if (is.null(block_nx) || is.null(block_ny)) {
      w  <- as.numeric(bb["xmax"] - bb["xmin"])
      h  <- as.numeric(bb["ymax"] - bb["ymin"])
      ratio <- if (h > 0) w / h else 1
      target_blocks <- max(1L, round(block_multiplier * k))
      nx <- max(1L, round(sqrt(target_blocks * ratio)))
      ny <- max(1L, round(max(1, target_blocks / nx)))

      # Diagnostic: warn if resulting block size is small relative to SAC range
      if (is.finite(sac_range) && sac_range > 0) {
        cell_w <- w / nx
        cell_h <- h / ny
        min_cell <- min(cell_w, cell_h)
        if (min_cell < sac_range) {
          .log_warn(
            "make_folds(block_kfold): geometric block size (%.1f) is smaller than the estimated autocorrelation range (%.1f). Consider setting block_size >= %.0f or auto_range = TRUE to avoid information leakage.",
            min_cell, sac_range, sac_range
          )
          warning(
            sprintf("make_folds(): block dimension (%.1f) < autocorrelation range (%.1f). Spatial CV may leak correlated information across folds. Pass block_size = %.0f or auto_range = TRUE.",
                    min_cell, sac_range, ceiling(sac_range)),
            call. = FALSE
          )
        }
      }
    } else {
      nx <- as.integer(block_nx); ny <- as.integer(block_ny)

      # Diagnostic: warn if user-supplied nx/ny yield blocks smaller than SAC
      if (is.finite(sac_range) && sac_range > 0) {
        w  <- as.numeric(bb["xmax"] - bb["xmin"])
        h  <- as.numeric(bb["ymax"] - bb["ymin"])
        cell_w <- w / nx; cell_h <- h / ny
        min_cell <- min(cell_w, cell_h)
        if (min_cell < sac_range) {
          .log_warn(
            "make_folds(block_kfold): user-supplied grid (%dx%d) yields blocks of ~%.1f units, smaller than estimated autocorrelation range (%.1f).",
            nx, ny, min_cell, sac_range
          )
          warning(
            sprintf("make_folds(): block_nx/block_ny yield blocks smaller than autocorrelation range (%.1f). Consider using block_size = %.0f.",
                    sac_range, ceiling(sac_range)),
            call. = FALSE
          )
        }
      }
    }

    grid <- sf::st_make_grid(reg, n = c(nx, ny), what = "polygons", square = TRUE)
    grid <- .safe_make_valid(grid)
    reg_union <- .safe_make_valid(sf::st_union(reg))
    grid <- suppressWarnings(sf::st_intersection(grid, reg_union))
    grid_sf <- sf::st_as_sf(grid)
    hits <- sf::st_intersects(pts, grid_sf)
    block_id <- vapply(hits, function(ix) if (length(ix)) ix[1] else NA_integer_, 1L)
    pts$..block_id <- block_id

    if (drop_empty_blocks) {
      used_blocks <- sort(unique(pts$..block_id[!is.na(pts$..block_id)]))
      grid_sf <- grid_sf[used_blocks, , drop = FALSE]
      pts$..block_id <- match(pts$..block_id, used_blocks)
    }
    if (anyNA(pts$..block_id)) {
      cent <- suppressWarnings(sf::st_centroid(sf::st_geometry(grid_sf)))
      na_idx <- which(is.na(pts$..block_id))
      dmat <- sf::st_distance(sf::st_geometry(pts[na_idx, ]), cent)
      pts$..block_id[na_idx] <- apply(as.matrix(dmat), 1, which.min)
    }
    B <- max(pts$..block_id, na.rm = TRUE)
    if (B < k) { .log_warn("make_folds(block_kfold): blocks < k; reducing k."); k <- B }

    blk_sizes <- as.integer(table(factor(pts$..block_id, levels = seq_len(B))))
    order_blk <- order(blk_sizes, decreasing = TRUE)
    fold_loads <- integer(k); fold_blocks <- vector("list", k)
    for (i in seq_along(order_blk)) {
      j <- which(fold_loads == min(fold_loads))
      if (length(j) > 1) j <- sample(j, 1L)
      fold_blocks[[j]] <- c(fold_blocks[[j]], seq_len(B)[order_blk[i]])
      fold_loads[j] <- fold_loads[j] + blk_sizes[order_blk[i]]
    }

    if (min(fold_loads) > 0L && max(fold_loads) > 3 * min(fold_loads)) {
      .log_warn("make_folds(block_kfold): fold size imbalance — largest fold has %d obs vs %d in smallest.",
                max(fold_loads), min(fold_loads))
    }

    splits <- vector("list", k); assign_vec <- integer(nrow(pts))
    for (j in seq_len(k)) {
      test_idx <- which(pts$..block_id %in% fold_blocks[[j]])
      splits[[j]] <- list(train = setdiff(seq_len(nrow(pts)), test_idx), test = test_idx)
      assign_vec[test_idx] <- j
    }
    return(.ret(method, k, splits,
                .safe_tibble(row_id = pts$..row_id, fold = assign_vec),
                list(seed = seed, grid_nx = nx, grid_ny = ny, blocks_used = B,
                     block_multiplier = block_multiplier,
                     block_size = block_size,
                     sac_range = sac_range,
                     auto_range = auto_range,
                     boundary_supplied = !is.null(boundary))))
  }

  # ---- BUFFERED LOO ----
  if (method == "buffered_loo") {
    if (is.null(buffer) || !is.numeric(buffer) || buffer <= 0)
      stop("make_folds(buffered_loo): `buffer` (positive numeric) is required.")
    pts <- ensure_projected(points_sf)
    n <- nrow(pts)
    if (n > 20000L) {
      stop(sprintf(
        "make_folds(buffered_loo): n = %d exceeds safety threshold of 20000. Although st_is_within_distance uses spatial indexing (GEOS STRtree), buffered LOO still produces n fold-splits and can be slow for very large datasets. Consider using 'block_kfold' instead, or subset your data.",
        n
      ), call. = FALSE)
    } else if (n > 5000L) {
      .log_warn(
        "make_folds(buffered_loo): n = %d is large; buffered LOO produces n fold-splits which may be slow.",
        n
      )
    }
    nb <- sf::st_is_within_distance(sf::st_geometry(pts), sf::st_geometry(pts),
                                     dist = buffer)
    splits <- vector("list", n)
    for (i in seq_len(n)) {
      excl <- sort(unique(nb[[i]]))
      splits[[i]] <- list(train = setdiff(seq_len(n), excl), test = i)
    }
    return(.ret(method, n, splits,
                .safe_tibble(row_id = pts$..row_id, fold = seq_len(n)),
                list(buffer = buffer)))
  }

  stop("make_folds(): unsupported method.")
}


# -----------------------------------------------------------------------------
# GWR Cross-Validation (now delegates to .cv_run_folds)
# -----------------------------------------------------------------------------

#' K-fold cross-validation for GWR
#'
#' @param data_sf An sf object.
#' @param response_var Response column name.
#' @param predictor_vars Predictor column names.
#' @param folds Optional list of fold definitions.
#' @param k Number of folds. Default 5.
#' @param seed RNG seed. Default 123.
#' @param adaptive Logical; use adaptive bandwidth. Default TRUE.
#' @param bandwidth Optional bandwidth value.
#' @param kernel Kernel function type.
#' @param boundary Optional polygonal sf/sfc for CRS alignment.
#' @param pointize Geometry coercion strategy.
#' @param block_size Optional minimum block edge length for spatial CV blocks
#'   (projected CRS units).  Ensures blocks are at least as large as the
#'   spatial autocorrelation range.
#' @param auto_range Logical.  If \code{TRUE} and \code{folds} is \code{NULL},
#'   estimate the autocorrelation range and use it as the minimum block size.
#'   Default \code{FALSE}.
#' @param parallel Logical or positive integer.  If \code{TRUE},
#'   auto-detect the number of cores and fit folds in parallel via
#'   \code{parallel::mclapply()} (macOS / Linux; falls back to sequential
#'   on Windows).  If an integer > 1, use that many cores.  Default
#'   \code{FALSE} (sequential).
#' @return A list with overall, fold_metrics, predictions, folds, formula, adaptive.
#' @export
cv_gwr <- function(data_sf, response_var, predictor_vars,
                   folds = NULL, k = 5, seed = 123,
                   adaptive = TRUE, bandwidth = NULL,
                   kernel = c("bisquare", "gaussian", "tricube",
                              "boxcar", "exponential"),
                   boundary = NULL, pointize = "auto",
                   block_size = NULL, auto_range = FALSE,
                   parallel = FALSE) {
  if (!inherits(data_sf, "sf")) stop("cv_gwr(): `data_sf` must be an sf object.")
  if (!requireNamespace("GWmodel", quietly = TRUE))
    stop("cv_gwr(): package 'GWmodel' is required.", call. = FALSE)
  if (!requireNamespace("sp", quietly = TRUE))
    stop("cv_gwr(): package 'sp' is required (for GWmodel interop).", call. = FALSE)

  kernel <- match.arg(kernel)
  kernel <- .validate_kernel(kernel)

  if (!("..row_id" %in% names(data_sf))) data_sf$`..row_id` <- seq_len(nrow(data_sf))
  dat_sf <- prep_model_data(data_sf, response_var, predictor_vars, boundary, pointize)
  if (!("..row_id" %in% names(dat_sf)))
    stop("cv_gwr(): `prep_model_data()` must preserve `..row_id`.")

  keep_idx <- dat_sf$`..row_id`

  if (is.null(folds)) {
    message("cv_gwr(): no folds supplied \u2014 using spatial block k-fold CV (k=", k, ").")
    folds <- make_folds(dat_sf, k = k, method = "block_kfold",
                        seed = seed, boundary = boundary,
                        block_size = block_size, auto_range = auto_range,
                        response_var = response_var,
                        predictor_vars = predictor_vars)
  }

  remapped_folds <- .remap_folds(folds, keep_idx, k, seed)

  cleanup_cv <- .with_seed(seed)
  on.exit(cleanup_cv(), add = TRUE)

  # Per-fold fitting function: fit_gwr_model returns a gwr_fit S3 object
  # Data is already prepped by the outer cv_gwr() call; skip the redundant
  # prep_model_data() pass inside fit_gwr_model() for each fold.
  fit_one <- function(train_sf) {
    fit_gwr_model(
      data_sf = train_sf, response_var = response_var,
      predictor_vars = predictor_vars,
      adaptive = adaptive, bandwidth = bandwidth, kernel = kernel,
      .already_prepped = TRUE
    )
  }

  # Extra per-fold info: bandwidth

  fold_info_fn <- function(fit_obj, test_sf, y_true, y_hat) {
    list(bandwidth = fit_obj$info$bandwidth %||% NA_real_)
  }

  # GWR fits (p+1) local coefficients at every training location, so the

  # effective number of parameters is >> the global predictor count p.
  # Passing p here would yield a per-fold Adj_R² that drastically overstates
  # parsimony.  We set p = NULL so that per-fold Adj_R² is reported as NA,
  # consistent with the pooled metric and with cv_bayes().
  res <- .cv_run_folds(
    dat_sf = dat_sf, response_var = response_var,
    predictor_vars = predictor_vars,
    remapped_folds = remapped_folds, keep_idx = keep_idx,
    fit_one = fit_one, fold_info_fn = fold_info_fn,
    p = NULL, parallel = parallel
  )

  preds <- if (length(res$pred_rows)) do.call(rbind, res$pred_rows) else
    data.frame(`..row_id` = integer(), fold = integer(),
               y = numeric(), yhat = numeric(), y_train_mean = numeric())
  folds_df <- if (length(res$fold_stats)) do.call(rbind, res$fold_stats) else
    data.frame(fold = integer(), n_train = integer(), n_test = integer(),
               n_pred = integer(),
               RMSE = numeric(), MAE = numeric(), MAPE = numeric(),
               SMAPE = numeric(), R2 = numeric(), Adj_R2 = numeric(),
               bandwidth = numeric())

  n_attempted <- length(remapped_folds)
  n_succeeded <- length(res$fold_stats)
  if (n_succeeded == 0L && n_attempted > 0L) {
    .log_warn("cv_gwr(): all %d folds failed to produce predictions; results are empty.", n_attempted)
    warning("cv_gwr(): all folds failed; cross-validation results contain no predictions.", call. = FALSE)
  } else if (n_succeeded < n_attempted) {
    .log_warn("cv_gwr(): %d of %d folds produced predictions.", n_succeeded, n_attempted)
  }

  list(overall = .cv_overall_metrics(preds), fold_metrics = folds_df,
       predictions = preds, folds = remapped_folds,
       formula = deparse(stats::reformulate(predictor_vars, response_var)),
       adaptive = adaptive)
}


# -----------------------------------------------------------------------------
# Bayesian Cross-Validation (now delegates to .cv_run_folds)
# -----------------------------------------------------------------------------

#' K-fold cross-validation for the Bayesian spatial model
#'
#' @param data_sf An sf object.
#' @param response_var Response column name.
#' @param predictor_vars Predictor column names.
#' @param folds Optional list of fold definitions.
#' @param k Number of folds. Default 5.
#' @param seed RNG seed. Default 123.
#' @param boundary Optional polygonal sf/sfc for CRS alignment.
#' @param pointize Geometry coercion strategy.
#' @param fit_args Named list of extra arguments for fit_bayesian_spatial_model().
#' @param summary "mean" or "median" for posterior predictions.
#' @param compute_pred_intervals Logical; compute predictive intervals.
#' @param coverage_levels Numeric vector of coverage levels.
#' @param block_size Optional minimum block edge length for spatial CV blocks
#'   (projected CRS units).
#' @param auto_range Logical.  If \code{TRUE} and \code{folds} is \code{NULL},
#'   estimate the autocorrelation range and use it as the minimum block size.
#'   Default \code{FALSE}.
#' @param parallel Logical or positive integer.  If \code{TRUE},
#'   auto-detect the number of cores and fit folds in parallel via
#'   \code{parallel::mclapply()} (macOS / Linux; falls back to sequential
#'   on Windows).  If an integer > 1, use that many cores.  Default
#'   \code{FALSE} (sequential).  Bayesian folds with full MCMC runs
#'   are the primary beneficiary of this option.
#' @return A list with overall, fold_metrics, predictions, folds, formula,
#'   and predictive_coverage.
#' @export
cv_bayes <- function(data_sf, response_var, predictor_vars,
                     folds = NULL, k = 5, seed = 123, boundary = NULL,
                     pointize = "auto", fit_args = list(),
                     summary = c("mean", "median"),
                     compute_pred_intervals = TRUE,
                     coverage_levels = c(0.50, 0.80, 0.95),
                     block_size = NULL, auto_range = FALSE,
                     parallel = FALSE) {
  summary <- match.arg(summary)
  if (!inherits(data_sf, "sf")) stop("cv_bayes(): `data_sf` must be an sf object.")

  if (!("..row_id" %in% names(data_sf))) data_sf$`..row_id` <- seq_len(nrow(data_sf))
  dat_sf <- prep_model_data(data_sf, response_var, predictor_vars, boundary, pointize)
  if (!("..row_id" %in% names(dat_sf)))
    stop("cv_bayes(): `prep_model_data()` must preserve `..row_id`.")

  keep_idx <- dat_sf$`..row_id`
  n_pred   <- length(predictor_vars)

  if (is.null(folds)) {
    message("cv_bayes(): no folds supplied \u2014 using spatial block k-fold CV (k=", k, ").")
    folds <- make_folds(dat_sf, k = k, method = "block_kfold",
                        seed = seed, boundary = boundary,
                        block_size = block_size, auto_range = auto_range,
                        response_var = response_var,
                        predictor_vars = predictor_vars)
  }

  remapped_folds <- .remap_folds(folds, keep_idx, k, seed)

  cleanup_cv <- .with_seed(seed)
  on.exit(cleanup_cv(), add = TRUE)

  # Per-fold fitting function: returns a bayesian_fit S3 object
  # Data is already prepped by the outer cv_bayes() call; skip the redundant
  # prep_model_data() pass inside fit_bayesian_spatial_model() for each fold.
  fit_one <- function(train_sf) {
    fold_fit_args <- modifyList(fit_args, list(
      compute_loo = FALSE, boundary = boundary,
      pointize = pointize, gp_k = NULL,
      .already_prepped = TRUE
    ))
    do.call(fit_bayesian_spatial_model,
            c(list(data_sf = train_sf, response_var = response_var,
                   predictor_vars = predictor_vars), fold_fit_args))
  }

  # Extra per-fold info: gp_k, CRPS, predictive coverage
  fold_info_fn <- function(fit_obj, test_sf, y_true, y_hat) {
    extras <- list(
      gp_k    = as.integer(fit_obj$info$gp_k %||% NA_integer_),
      n_draws = NA_integer_,
      CRPS    = NA_real_
    )

    # Full posterior predictive draws for intervals + CRPS
    if (isTRUE(compute_pred_intervals)) {
      ppred_draws <- try(
        predict(fit_obj, newdata = test_sf, type = "predict", draws = TRUE),
        silent = TRUE
      )
      if (!inherits(ppred_draws, "try-error") && is.matrix(ppred_draws)) {
        extras$n_draws <- nrow(ppred_draws)

        # Coverage at each level
        for (cl in coverage_levels) {
          alpha <- (1 - cl) / 2
          lwr <- apply(ppred_draws, 2L, stats::quantile, probs = alpha)
          upr <- apply(ppred_draws, 2L, stats::quantile, probs = 1 - alpha)
          in_interval <- y_true >= lwr & y_true <= upr
          extras[[sprintf("coverage_%.0f", cl * 100)]] <- mean(in_interval, na.rm = TRUE)
        }

        # Empirical CRPS via the NRG (energy) form, vectorised over
        # observations.  CRPS = E|X - y| - 0.5 * E|X - X'|, where the
        # spread term 0.5 * E|X - X'| = (1/m^2) * sum_i x_{(i)} * (2i - m - 1)
        # for m equally-weighted posterior draws (Gneiting & Raftery, 2007,
        # Eq. 21; see also the scoringRules package documentation).
        #
        # sort all columns in one call, then use matrix ops
        m <- nrow(ppred_draws)
        mae_terms   <- colMeans(abs(sweep(ppred_draws, 2L, y_true)))
        sorted_draws <- if (requireNamespace("matrixStats", quietly = TRUE)) {
          matrixStats::colSort(ppred_draws)            # m x n_test
        } else {
          apply(ppred_draws, 2L, sort)                 # m x n_test
        }
        weights      <- (2 * seq_len(m) - m - 1) / (m * m)
        spread_terms <- colSums(sorted_draws * weights)
        crps_per_obs <- mae_terms - spread_terms
        extras$CRPS  <- mean(crps_per_obs, na.rm = TRUE)
      }
    }

    extras
  }

  res <- .cv_run_folds(
    dat_sf = dat_sf, response_var = response_var,
    predictor_vars = predictor_vars,
    remapped_folds = remapped_folds, keep_idx = keep_idx,
    fit_one = fit_one, fold_info_fn = fold_info_fn,
    predict_args = list(summary = summary),
    p = NULL, parallel = parallel
  )

  # Assemble predictions — add yhat_sd column via a second pass
  pred_rows_with_sd <- lapply(res$pred_rows, function(pr) {
    pr$yhat_sd <- NA_real_  # placeholder; could be enriched
    pr
  })

  preds <- if (length(pred_rows_with_sd)) do.call(rbind, pred_rows_with_sd) else
    data.frame(`..row_id` = integer(), fold = integer(),
               y = numeric(), yhat = numeric(), y_train_mean = numeric(),
               yhat_sd = numeric())
  folds_df <- if (length(res$fold_stats)) do.call(rbind, res$fold_stats) else
    data.frame(fold = integer(), n_train = integer(), n_test = integer(),
               n_pred = integer(),
               RMSE = numeric(), MAE = numeric(), MAPE = numeric(),
               SMAPE = numeric(), R2 = numeric(), Adj_R2 = numeric(),
               CRPS = numeric(), gp_k = integer(), n_draws = integer())

  n_attempted <- length(remapped_folds)
  n_succeeded <- length(res$fold_stats)
  if (n_succeeded == 0L && n_attempted > 0L) {
    .log_warn("cv_bayes(): all %d folds failed to produce predictions; results are empty.", n_attempted)
    warning("cv_bayes(): all folds failed; cross-validation results contain no predictions.", call. = FALSE)
  } else if (n_succeeded < n_attempted) {
    .log_warn("cv_bayes(): %d of %d folds produced predictions.", n_succeeded, n_attempted)
  }

  list(overall = .cv_overall_metrics(preds), fold_metrics = folds_df,
       predictions = preds, folds = remapped_folds,
       formula = deparse(stats::reformulate(predictor_vars, response_var)),
       predictive_coverage = if (nrow(folds_df) > 0L && "CRPS" %in% names(folds_df)) {
         cov_cols <- grep("^coverage_", names(folds_df), value = TRUE)
         cov_means <- vapply(cov_cols, function(cn) mean(folds_df[[cn]], na.rm = TRUE), numeric(1))
         c(as.list(cov_means), mean_CRPS = mean(folds_df$CRPS, na.rm = TRUE))
       } else NULL)
}


# -----------------------------------------------------------------------------
# Unified CV (model-agnostic)
# -----------------------------------------------------------------------------

#' Model-agnostic spatial cross-validation
#'
#' Run K-fold CV for any model that returns a \code{spatial_fit} object.
#' This is the extensibility point: to plug in a new model type, supply
#' a \code{fit_fn(train_sf)} that returns a \code{spatial_fit}.
#'
#' @param data_sf An sf object.
#' @param response_var Response column name.
#' @param predictor_vars Predictor column names.
#' @param fit_fn A function(train_sf) that returns a \code{spatial_fit}.
#' @param folds Optional fold definitions. Built via block_kfold if NULL.
#' @param k Number of folds.
#' @param seed RNG seed.
#' @param boundary Optional boundary for fold construction.
#' @param pointize Geometry coercion strategy.
#' @param predict_args Extra arguments for predict().
#' @param fold_info_fn Optional function for per-fold extras.
#' @param p Number of predictors for Adj R² (NULL to skip).  Only meaningful
#'   for models with a fixed global parameter count; pass NULL for models
#'   with spatially varying coefficients (e.g. GWR).
#' @param block_size Optional minimum block edge length for spatial CV blocks
#'   (projected CRS units).
#' @param auto_range Logical.  If \code{TRUE} and \code{folds} is \code{NULL},
#'   estimate the autocorrelation range and use it as the minimum block size.
#'   Default \code{FALSE}.
#' @param parallel Logical or positive integer.  If \code{TRUE},
#'   auto-detect the number of cores and fit folds in parallel via
#'   \code{parallel::mclapply()} (macOS / Linux; falls back to sequential
#'   on Windows).  If an integer > 1, use that many cores.  Default
#'   \code{FALSE} (sequential).
#' @return A list with overall, fold_metrics, predictions, folds.
#' @export
cv_spatial <- function(data_sf, response_var, predictor_vars,
                       fit_fn, folds = NULL, k = 5, seed = 123,
                       boundary = NULL, pointize = "auto",
                       predict_args = list(), fold_info_fn = NULL,
                       p = NULL, block_size = NULL,
                       auto_range = FALSE, parallel = FALSE) {
  if (!inherits(data_sf, "sf")) stop("cv_spatial(): `data_sf` must be an sf object.")
  if (!is.function(fit_fn)) stop("cv_spatial(): `fit_fn` must be a function.")

  if (!("..row_id" %in% names(data_sf))) data_sf$`..row_id` <- seq_len(nrow(data_sf))
  dat_sf <- prep_model_data(data_sf, response_var, predictor_vars, boundary, pointize)
  keep_idx <- dat_sf$`..row_id`

  if (is.null(folds)) {
    message("cv_spatial(): no folds supplied \u2014 using spatial block k-fold CV (k=", k, ").")
    folds <- make_folds(dat_sf, k = k, method = "block_kfold",
                        seed = seed, boundary = boundary,
                        block_size = block_size, auto_range = auto_range,
                        response_var = response_var,
                        predictor_vars = predictor_vars)
  }

  remapped_folds <- .remap_folds(folds, keep_idx, k, seed)

  cleanup_cv <- .with_seed(seed)
  on.exit(cleanup_cv(), add = TRUE)

  res <- .cv_run_folds(
    dat_sf = dat_sf, response_var = response_var,
    predictor_vars = predictor_vars,
    remapped_folds = remapped_folds, keep_idx = keep_idx,
    fit_one = fit_fn, fold_info_fn = fold_info_fn,
    predict_args = predict_args, p = p,
    parallel = parallel
  )

  preds <- if (length(res$pred_rows)) do.call(rbind, res$pred_rows) else
    data.frame(`..row_id` = integer(), fold = integer(),
               y = numeric(), yhat = numeric(), y_train_mean = numeric())
  folds_df <- if (length(res$fold_stats)) do.call(rbind, res$fold_stats) else
    data.frame()

  list(overall = .cv_overall_metrics(preds),
       fold_metrics = folds_df, predictions = preds,
       folds = remapped_folds)
}
