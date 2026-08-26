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
#' @return List of \code{list(train, test, fold_id)} with IDs restricted to
#'   \code{keep_idx}.  \code{fold_id} is the fold's index in the ORIGINAL
#'   \code{folds} object, carried through so that dropping an unusable fold
#'   does not renumber the survivors: downstream \code{fold} columns then still
#'   agree with \code{make_folds()$assignment$fold}.
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
    assign_vec <- sample(rep(seq_len(k), length.out = length(keep_idx)))
    out <- vector("list", k)
    for (i in seq_len(k))
      out[[i]] <- list(train = keep_idx[assign_vec != i],
                       test  = keep_idx[assign_vec == i],
                       fold_id = i)
    return(out)
  }

  if (!is.null(folds$folds) && is.list(folds$folds)) {
    folds <- folds$folds
  }

  # Carry the ORIGINAL fold index alongside each split.  Dropping an unusable
  # fold below removes an element from the list, and without this every later
  # fold would be silently renumbered, so fold_metrics$fold and
  # predictions$fold would no longer line up with make_folds()$assignment$fold.
  remapped <- lapply(seq_along(folds), function(j) {
    f <- folds[[j]]
    list(
      train = keep_idx[stats::na.omit(match(f$train, keep_idx))],
      test  = keep_idx[stats::na.omit(match(f$test, keep_idx))],
      fold_id = j
    )
  })

  # An empty TRAINING set is just as fatal as an empty test set and was not
  # checked at all: buffered_loo with a buffer spanning the data, or
  # random_kfold on a single row, produces folds nothing can be fitted on.
  # They used to sail through here and get dropped one at a time inside
  # .cv_fit_one_fold(), so the only symptom was a generic "all folds failed"
  # warning at the very end.
  no_test  <- vapply(remapped, function(f) length(f$test) == 0L, logical(1))
  no_train <- vapply(remapped, function(f) length(f$train) < 2L, logical(1))
  if (any(no_test))
    .log_warn(".remap_folds(): %d fold(s) have empty test sets after remapping.",
              sum(no_test))
  if (any(no_train))
    .log_warn(".remap_folds(): %d fold(s) have fewer than 2 training rows after remapping and cannot be fitted.",
              sum(no_train))

  remapped[!(no_test | no_train)]
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
  # The label to report: the fold's index in the ORIGINAL fold object, so that
  # dropping an unusable fold in .remap_folds() does not renumber the rest.
  # Falls back to the loop position for hand-built folds carrying no fold_id.
  fold_lab <- remapped_fold$fold_id %||% i

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
    .log_warn(".cv_run_folds(): fold %d fit failed; skipping.", fold_lab)
    return(NULL)
  }
  if (!inherits(fit_obj, "spatial_fit")) {
    .log_warn(".cv_run_folds(): fold %d did not return a spatial_fit; skipping.", fold_lab)
    return(NULL)
  }

  # Predict on test fold via the S3 generic
  y_true <- test_df[[response_var]]
  y_hat  <- try(
    do.call(predict, c(list(object = fit_obj, newdata = test_sf), predict_args)),
    silent = TRUE
  )
  if (inherits(y_hat, "try-error") || !is.numeric(y_hat)) {
    .log_warn(".cv_run_folds(): fold %d predict failed; skipping.", fold_lab)
    return(NULL)
  }
  # cv_spatial() is documented as the extensibility point for arbitrary
  # learners, so a fit_fn whose predict() returns the wrong length is
  # user-reachable.  Nothing downstream notices: .compute_reg_metrics() and
  # the data.frame() below both RECYCLE silently, so a length-2 y_hat against
  # 4 test rows yields a 4-row frame with the predictions repeated and metrics
  # computed against fabricated pairs.
  if (length(y_hat) != length(y_true)) {
    .log_warn(".cv_run_folds(): fold %d predicted %d value(s) for %d test row(s); skipping.",
              fold_lab, length(y_hat), length(y_true))
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
    fold = fold_lab, n_train = length(tr_pos), n_test = length(y_true),
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
    `..row_id` = test_sf$`..row_id`, fold = fold_lab,
    y = as.numeric(y_true), yhat = as.numeric(y_hat),
    y_train_mean = y_train_mean,
    stringsAsFactors = FALSE
  )

  list(pred_row = pr, fold_stat = fs)
}


#' Empirical CRPS via the energy (NRG) form for equally weighted draws
#'
#' Computes the continuous ranked probability score per observation:
#' CRPS(F, y) = E|X - y| - 0.5 * E|X - X'|, where X, X' are independent
#' draws from the predictive distribution.  The spread term uses the
#' Gini-mean-difference identity on sorted draws:
#' 0.5 * E|X - X'| = (1/m^2) * sum_i x_(i) * (2i - m - 1)
#' for m equally weighted draws (Gneiting & Raftery 2007, Eq. 21; see also
#' the scoringRules package documentation).
#'
#' @param draws Numeric matrix (m draws x n observations) of posterior
#'   predictive draws.
#' @param y Numeric vector of length n of observed values.
#' @return Numeric vector of length n with per-observation CRPS.
#' @keywords internal
#' @noRd
.crps_energy <- function(draws, y) {
  m <- nrow(draws)
  mae_terms <- colMeans(abs(sweep(draws, 2L, y)))
  # Vectorised within-column sort in base R: order by column, then by value.
  sorted_draws <- matrix(draws[order(col(draws), draws)], nrow = m)
  weights <- (2 * seq_len(m) - m - 1) / (m * m)
  spread_terms <- colSums(sorted_draws * weights)
  mae_terms - spread_terms
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
    eff <- .sanitize_core_count(n_cores)
  } else if (isTRUE(parallel)) {
    # detectCores() can return NA on some platforms; sanitize before use.
    eff <- .sanitize_core_count(parallel::detectCores(logical = FALSE) - 1L)
  } else if (is.numeric(parallel) && length(parallel) == 1L &&
             !is.na(parallel) && parallel > 1) {
    eff <- .sanitize_core_count(parallel)
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
#' @param seed Integer RNG seed, or \code{NULL} to leave fold RNG unseeded.
#'   When supplied, one seed per fold is drawn in the parent process and
#'   applied inside the fold worker, so results depend on (seed, fold index)
#'   alone and \code{parallel = TRUE} reproduces \code{parallel = FALSE}
#'   exactly.
#' @return List with pred_rows and fold_stats.
#' @keywords internal
#' @noRd
.cv_run_folds <- function(dat_sf, response_var, predictor_vars,
                          remapped_folds, keep_idx, fit_one,
                          fold_info_fn = NULL, predict_args = list(),
                          p = NULL, parallel = FALSE, n_cores = NULL,
                          seed = NULL) {
  cores   <- .resolve_n_cores(parallel, n_cores)
  n_folds <- length(remapped_folds)

  # Draw one seed per fold in the parent so that each fold's RNG stream is a
  # function of (seed, fold index) only -- never of the execution path.  This
  # is what makes parallel and sequential runs bit-identical.  Forked workers
  # cannot guarantee that on their own: mclapply() seeds each child from the
  # current time and process ID unless the L'Ecuyer-CMRG generator is in use,
  # so without this the parallel path is irreproducible for any fit_fn that
  # consumes RNG (see cv_spatial(), which accepts an arbitrary learner).
  fold_seeds <- if (is.null(seed)) {
    rep(NA_integer_, n_folds)
  } else {
    cleanup_draw <- .with_seed(seed)
    on.exit(cleanup_draw(), add = TRUE)
    sample.int(.Machine$integer.max, n_folds)
  }

  fold_worker <- function(i) {
    # .with_seed() saves and restores .Random.seed, so seeding a fold never
    # leaks into the caller's RNG state.
    if (!is.na(fold_seeds[i])) {
      cleanup_i <- .with_seed(fold_seeds[i])
      on.exit(cleanup_i(), add = TRUE)
    }
    .cv_fit_one_fold(
      i = i, dat_sf = dat_sf, response_var = response_var,
      remapped_fold = remapped_folds[[i]], keep_idx = keep_idx,
      fit_one = fit_one, fold_info_fn = fold_info_fn,
      predict_args = predict_args, p = p
    )
  }

  if (cores > 1L) {
    message(sprintf("cv: running %d folds in parallel on %d cores.",
                    n_folds, cores))
    results <- parallel::mclapply(
      seq_along(remapped_folds), fold_worker, mc.cores = cores
    )
  } else {
    results <- lapply(seq_along(remapped_folds), fold_worker)
  }

  # mclapply() hands back a try-error OBJECT (not NULL) when a child errors or
  # is killed, and Negate(is.null) keeps it -- the subsequent [[ "pred_row" ]]
  # then fails with "subscript out of bounds", destroying the real diagnosis.
  failed <- vapply(results, function(z) inherits(z, "try-error"), logical(1))
  if (any(failed)) {
    # conditionMessage() has no method for a "try-error" object -- calling it
    # on one throws, which would destroy the diagnosis in precisely the branch
    # that only runs when a worker has already failed.  The condition object
    # hangs off the try-error as an attribute; the deparsed string is the
    # fallback for the rare try-error carrying none.
    msgs <- vapply(results[failed], function(z) {
      cond <- attr(z, "condition")
      if (is.null(cond)) trimws(paste(as.character(z), collapse = " "))
      else conditionMessage(cond)
    }, character(1))
    .log_warn("cv: %d fold(s) failed in a parallel worker: %s",
              sum(failed), paste(unique(msgs), collapse = "; "))
  }
  results <- results[!failed]
  results <- Filter(Negate(is.null), results)
  pred_rows  <- lapply(results, `[[`, "pred_row")
  fold_stats <- lapply(results, `[[`, "fold_stat")

  list(pred_rows = pred_rows, fold_stats = fold_stats)
}


# -----------------------------------------------------------------------------
# Spatial autocorrelation range estimation
# -----------------------------------------------------------------------------

#' Nearest-neighbour distance from each feature of \code{query} to \code{data}
#'
#' Uses \code{sf::st_nearest_feature()} so the search is indexed rather than a
#' dense cross-distance matrix, which would be prohibitive at the sizes NNDM is
#' guarded to.
#'
#' @param query,data sf layers in a common CRS.
#' @return Numeric vector, one distance per row of \code{query}.
#' @keywords internal
#' @noRd
.nn_dist_to <- function(query, data) {
  idx <- sf::st_nearest_feature(query, data)
  as.numeric(sf::st_distance(sf::st_geometry(query),
                             sf::st_geometry(data)[idx],
                             by_element = TRUE))
}


#' Estimate the spatial autocorrelation range from data
#'
#' Fits exponential (or spherical) variogram models and returns the
#' \emph{effective range} — the distance at which the semivariance reaches
#' ~95 \% of the sill.
#'
#' To guard against anisotropy, the function first estimates directional
#' variograms at 0° (N–S) and 90° (E–W) azimuths (tolerance 22.5°,
#' which avoids double-counting point pairs near the 45° diagonal but
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
#' @param range_frac Positive numeric.  A fitted range exceeding
#'   \code{range_frac * cutoff * max_dist} -- that is, beyond the longest lag
#'   the empirical variogram was actually fitted over -- is treated as
#'   unidentified and \code{NA_real_} is returned.
#'   \code{gstat::fit.variogram()} yields a finite number even when the
#'   variogram never reaches a sill, and such a value is extrapolation past the
#'   observed lags rather than a long autocorrelation range.  Passing it to
#'   \code{make_folds(auto_range = TRUE)} would collapse the block grid to a
#'   single block.  Default 1.0; raise it to accept ranges extrapolated beyond
#'   the fitted lags.
#' @param seed Optional RNG seed for subsampling reproducibility.
#' @return A single value of class \code{sac_range}, which behaves as an
#'   ordinary number.  There are three shapes, and they carry different
#'   attributes:
#'   \describe{
#'     \item{Success}{A positive effective range in projected coordinate units,
#'       with the fit attached as attributes \code{directional} (the 0° and 90°
#'       ranges), \code{anisotropy} (their ratio), \code{max_dist},
#'       \code{cutoff_dist}, \code{variogram} (the empirical variogram) and
#'       \code{variogram_model} (the fitted \code{gstat} model), so the fit can
#'       be inspected rather than trusted.}
#'     \item{Rejected range}{\code{NA_real_} when a range was fitted but exceeds
#'       \code{range_frac * cutoff * max_dist} and is therefore unidentified
#'       (see \code{range_frac}).  It is classed \code{sac_range} as well, so it
#'       prints as a bare \code{NA} rather than dumping its attributes, and it
#'       carries \code{max_dist}, \code{cutoff_dist}, \code{variogram} and
#'       \code{variogram_model} --- the evidence for the rejection --- plus
#'       \code{rejected_range} (the value that was refused) and
#'       \code{rejected_reason}.  It does \strong{not} carry \code{directional}
#'       or \code{anisotropy}.}
#'     \item{No fit}{A bare, attribute-less \code{NA_real_} when estimation
#'       could not be attempted or produced nothing at all (\pkg{gstat}
#'       missing, fewer than 30 finite values, a degenerate extent, or a
#'       singular variogram fit).}
#'   }
#'   Attributes and the class do not affect \code{is.na()} or
#'   \code{is.finite()}, so every downstream guard treats all three the same
#'   way it always did.
#' @family cross-validation
#' @examples
#' if (requireNamespace("gstat", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 80
#'   xy <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000))
#'   xy$z <- sin(xy$x / 200) + rnorm(n, sd = 0.2)
#'   pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
#'   estimate_sac_range(pts, response_var = "z")
#' }
#' @export
estimate_sac_range <- function(points_sf, response_var,
                               predictor_vars = NULL,
                               n_max = 5000L, cutoff = 0.5,
                               range_frac = 1.0, seed = NULL) {
  if (!requireNamespace("gstat", quietly = TRUE)) {
    .log_warn("estimate_sac_range(): package 'gstat' is required for variogram estimation; returning NA.")
    return(NA_real_)
  }
  if (!inherits(points_sf, "sf"))
    stop("estimate_sac_range(): `points_sf` must be an sf object.", call. = FALSE)
  # MULTIPOINT must be coerced too, not merely admitted: st_coordinates()
  # returns one row per VERTEX, so any multi-vertex feature makes xy[i, ] a
  # different feature than row i, and every fold below misaligns silently.
  # coerce_to_points() takes centroids, matching prep_model_data().
  if (!all(sf::st_geometry_type(points_sf, by_geometry = TRUE) == "POINT"))
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
  # gstat::fit.variogram() signals a failed fit by setting attr(., "singular")
  # to TRUE and emitting a warning -- it returns normally rather than throwing.
  # Testing only for a try-error therefore (a) made the "Sph" fallback below
  # unreachable, and (b) let a singular fit's `range` flow straight out as the
  # estimated autocorrelation range, which make_folds(auto_range = TRUE) then
  # sizes spatial blocks from.  The `max_supported` guard further down only
  # catches ranges beyond the fitted lags, not fitting artefacts inside them.
  .fit_one_vgm <- function(vg, model_type) {
    m <- try(gstat::fit.variogram(vg, gstat::vgm(model = model_type)),
             silent = TRUE)
    if (inherits(m, "try-error") || !is.data.frame(m)) return(NULL)
    if (isTRUE(attr(m, "singular"))) return(NULL)
    m
  }

  .fit_vgm_range <- function(vg) {
    if (inherits(vg, "try-error") || !inherits(vg, "data.frame") || NROW(vg) < 3L) return(NA_real_)
    vgm_model <- .fit_one_vgm(vg, "Exp")
    if (is.null(vgm_model)) vgm_model <- .fit_one_vgm(vg, "Sph")
    if (is.null(vgm_model)) {
      # Both the exponential and the spherical fit were singular or errored:
      # there is no identified range, so say so rather than returning one.
      .log_warn("estimate_sac_range(): both the exponential and the spherical fit to this variogram are singular; no range is identified from it.")
      return(NA_real_)
    }
    spatial_rows <- vgm_model[vgm_model$model != "Nug", , drop = FALSE]
    if (nrow(spatial_rows) == 0L || all(!is.finite(spatial_rows$range)))
      return(NA_real_)
    raw_range  <- max(spatial_rows$range, na.rm = TRUE)
    model_type <- spatial_rows$model[which.max(spatial_rows$range)]
    eff <- if (identical(as.character(model_type), "Exp")) 3 * raw_range else raw_range
    if (!is.finite(eff) || eff <= 0) return(NA_real_)
    # Carry the fitted model out so callers can inspect the fit rather than
    # trust a bare number.
    structure(eff, vgm_model = vgm_model)
  }

  # --- Directional variograms (0° and 90°, tolerance 22.5°) ----------------
  # gstat uses azimuth in degrees clockwise from north.  0° = N-S, 90° = E-W.
  dir_fits <- lapply(c(0, 90), function(az) {
    vg_dir <- try(
      gstat::variogram(..sac_var ~ 1, data = pts,
                       cutoff = cutoff_dist,
                       alpha = az, tol.hor = 22.5),
      silent = TRUE
    )
    list(vg = vg_dir, fit = .fit_vgm_range(vg_dir))
  })
  dir_ranges <- vapply(dir_fits, function(f) as.numeric(f$fit), numeric(1))

  dir_success <- sum(is.finite(dir_ranges)) == 2L
  anisotropy  <- NA_real_
  vg_used     <- NULL
  vgm_used    <- NULL

  # --- Select the effective range ------------------------------------------
  if (dir_success) {
    effective_range <- max(dir_ranges)
    anisotropy <- max(dir_ranges) / min(dir_ranges)
    winner   <- which.max(dir_ranges)
    vg_used  <- dir_fits[[winner]]$vg
    vgm_used <- attr(dir_fits[[winner]]$fit, "vgm_model")
    if (anisotropy > 1.5) {
      .log_warn(
        "estimate_sac_range(): notable anisotropy detected (range ratio %.1f). Directional ranges: 0\u00b0 = %.1f, 90\u00b0 = %.1f. Using the maximum.",
        anisotropy, dir_ranges[1], dir_ranges[2]
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
      # as.numeric() strips the vgm_model attribute .fit_vgm_range() attaches;
      # it is re-attached under its documented name below, and leaving both
      # would ship the same object under two attribute names.
      effective_range <- as.numeric(iso_range)
      vg_used  <- vg_iso
      vgm_used <- attr(iso_range, "vgm_model")
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

  # --- Reject ranges the data cannot actually support -----------------------
  # gstat::fit.variogram() returns a finite number even when the empirical
  # variogram never reaches a sill.  The range is then unidentified and the
  # value is a fitting artefact, not an autocorrelation range.
  #
  # The test is against the CUTOFF distance, not the extent: the empirical
  # variogram is only computed out to `cutoff * max_dist`, so a fitted range
  # beyond that is extrapolation past the longest lag ever observed.  No amount
  # of data at shorter lags identifies it.  (With the default cutoff = 0.5 this
  # is half the extent, but the rationale is the fitted lag limit rather than an
  # arbitrary fraction of the study area.)
  #
  # Returning such a value silently is worse than returning NA, because
  # make_folds(auto_range = TRUE) sizes spatial blocks from it: a range
  # spanning the data yields one block covering everything, which silently
  # defeats blocked cross-validation.
  max_supported <- range_frac * cutoff_dist
  if (is.finite(max_supported) && effective_range > max_supported) {
    .log_warn(
      paste0("estimate_sac_range(): fitted range (%.0f) exceeds the largest ",
             "lag the variogram was fitted over (%.4g = %s x cutoff %.0f); the ",
             "empirical variogram never reached a sill, so the range is ",
             "unidentified rather than long. Returning NA. Raise `cutoff` to ",
             "fit longer lags, supply `predictor_vars` to detrend, or set a ",
             "block size explicitly."),
      effective_range, max_supported, format(range_frac), cutoff_dist
    )
    # The VALUE is NA -- the range is genuinely unidentified and must not be
    # used to size blocks -- but the variogram that justified the rejection has
    # already been computed, and throwing it away leaves the user no way to see
    # WHY.  It is exactly the case where plot(type = "variogram") is most worth
    # looking at: a curve that never reaches a sill.  Attributes do not affect
    # is.finite()/is.na(), so every downstream guard behaves as before.
    #
    # The class is set here too, exactly as on the success return.  Without it
    # print.sac_range() never fired for a rejected range, so printing one
    # dumped the whole empirical variogram and the fitted gstat model to the
    # console as raw attributes instead of showing "NA".
    return(structure(
      NA_real_,
      class           = c("sac_range", "numeric"),
      max_dist        = as.numeric(max_dist),
      cutoff_dist     = as.numeric(cutoff_dist),
      variogram       = vg_used,
      variogram_model = vgm_used,
      rejected_range  = as.numeric(effective_range),
      rejected_reason = "fitted range exceeds the largest lag fitted"
    ))
  }

  structure(
    effective_range,
    class           = c("sac_range", "numeric"),
    directional     = stats::setNames(dir_ranges, c("0", "90")),
    anisotropy      = anisotropy,
    max_dist        = as.numeric(max_dist),
    cutoff_dist     = as.numeric(cutoff_dist),
    variogram       = vg_used,
    variogram_model = vgm_used
  )
}


#' Print a spatial autocorrelation range
#'
#' Prints the effective range as a plain number, with the directional fit
#' summarised beneath it when one is available.
#'
#' @param x An object of class \code{sac_range}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.sac_range <- function(x, ...) {
  cat(format(as.numeric(x)), "\n")
  d <- attr(x, "directional")
  a <- attr(x, "anisotropy")
  if (!is.null(d) && all(is.finite(d))) {
    cat(sprintf("  directional: 0 deg = %s, 90 deg = %s",
                format(d[[1]]), format(d[[2]])))
    if (is.finite(a)) cat(sprintf("  (ratio %.2f)", a))
    cat("\n")
  }
  invisible(x)
}


#' Compute grid dimensions that respect a minimum block size
#'
#' Given a bounding box and a minimum block edge length, returns nx/ny
#' values such that each cell is at least \code{block_size} wide and tall.
#'
#' @param bb An sf bbox.
#' @param block_size Positive numeric minimum block edge length (CRS units).
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
#' @param k Integer; number of folds.  Required for every method, but not every
#'   method honours it: \code{"buffered_loo"} and \code{"nndm"} are
#'   leave-one-out schemes and always return \code{k = n} regardless of what
#'   was asked for, while \code{"block_kfold"} lowers it when the grid yields
#'   fewer than \code{k} non-empty blocks and \code{"leave_location_out"}
#'   lowers it when there are fewer than \code{k} distinct groups.  Read the
#'   \code{k} element of the returned list rather than assuming the requested
#'   value; a reduction is logged.
#' @param method One of \code{"random_kfold"}, \code{"block_kfold"},
#'   \code{"buffered_loo"}, \code{"leave_location_out"} or \code{"nndm"}.  See
#'   \strong{Details} for what each one does and when it is appropriate.
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
#' @param range_frac Passed through to \code{estimate_sac_range()} when
#'   \code{auto_range = TRUE}.  A fitted range beyond the longest lag the
#'   empirical variogram was fitted over is rejected as unidentified, and block
#'   sizing falls back to geometry rather than collapsing to a single block.
#'   Default 1.0.
#' @param response_var Character(1) response column name.  Required when
#'   \code{auto_range = TRUE}.
#' @param predictor_vars Optional character vector of predictor column names.
#'   Passed to \code{estimate_sac_range()} for residual variogram estimation.
#' @param boundary Optional polygonal sf/sfc for block_kfold.
#' @param buffer Positive numeric distance for buffered_loo.
#' @param group_var Character(1) naming a column of \code{points_sf} that
#'   identifies the location each observation belongs to.  Required for
#'   \code{method = "leave_location_out"}, which keeps every observation from a
#'   location together in the same fold.  Repeated measurements at the same
#'   site otherwise get split across folds, and the model is scored partly on
#'   sites it has already seen -- which random k-fold reports as excellent
#'   performance.
#' @param prediction_points Optional \code{sf} layer of the locations you
#'   actually intend to predict onto.  Required for \code{method = "nndm"}.
#'   The grid from \code{predict_surface()} is the natural choice.
#' @details
#' \strong{Fold methods.}
#' \code{"random_kfold"} ignores geography entirely and will overstate
#' performance on autocorrelated data.  \code{"block_kfold"} separates folds
#' geographically.  \code{"buffered_loo"} holds out one point at a time and
#' excludes everything within a fixed \code{buffer}.
#'
#' \code{"leave_location_out"} groups by \code{group_var}, so all
#' observations from a location share a fold.
#'
#' \code{"nndm"} implements the distance-matching principle of Milà et al.
#' (2022): rather than choosing a buffer arbitrarily, it sizes the exclusion
#' around each held-out point so that the resulting training-to-test distance
#' distribution approaches the distribution of distances from your actual
#' prediction locations to the training data.
#'
#' Each held-out point draws a target distance from the empirical
#' prediction-distance CDF, then the nearest training points are excluded up to
#' the one whose distance is closest to that target.  Matching is therefore
#' \emph{as close as the training configuration permits}, not exact: the
#' achievable distances are the order statistics of that point's neighbour
#' distances, which are discrete and can leave gaps.  Clustered data is the
#' hard case -- excluding everything inside a hard radius would remove a whole
#' cluster and leave the nearest survivor in the next one, overshooting badly,
#' which is why the nearest achievable distance is used instead.
#'
#' This differs from the paper's iterative exclusion procedure.  Compare
#' \code{params$target_median} against \code{params$realised_median} to see
#' how well the matching worked on your data.
#' @references
#' Mila, C., Mateu, J., Pebesma, E. and Meyer, H. (2022). Nearest neighbour
#' distance matching Leave-One-Out Cross-Validation for map validation.
#' \emph{Methods in Ecology and Evolution} \strong{13}, 1304-1316.
#' \doi{10.1111/2041-210X.13851}
#' @param drop_empty_blocks Logical. Default TRUE.
#' @return A list with method, k, folds, assignment, params.  The
#'   \code{train}/\code{test} elements of each fold contain \code{..row_id}
#'   values (equal to row positions when the input has no pre-existing
#'   \code{..row_id} column), consistent with the \code{assignment} tibble.
#'   The returned \code{k} is the number of folds actually built, which is not
#'   always the \code{k} that was requested (see the \code{k} argument above),
#'   and \code{length(folds)} always matches it.
#' @family cross-validation
#' @examples
#' library(sf)
#' set.seed(1)
#' pts <- st_as_sf(
#'   data.frame(x = runif(30, 0, 1000), y = runif(30, 0, 1000)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' folds <- make_folds(pts, k = 3, method = "block_kfold", seed = 42)
#' folds$assignment          # fold membership per row
#' lengths(folds$folds[[1]]) # train/test row-ID splits
#'
#' # Buffered leave-one-out: neighbours within 100 units excluded from training
#' loo <- make_folds(pts, k = 1, method = "buffered_loo", buffer = 100)
#' @export
make_folds <- function(points_sf, k,
                       method = c("random_kfold", "block_kfold", "buffered_loo",
                                  "leave_location_out", "nndm"),
                       seed = NULL, block_nx = NULL, block_ny = NULL,
                       block_multiplier = 3, block_size = NULL,
                       auto_range = FALSE, range_frac = 1.0, response_var = NULL,
                       group_var = NULL, prediction_points = NULL,
                       predictor_vars = NULL, boundary = NULL,
                       buffer = NULL,
                       drop_empty_blocks = TRUE) {
  method <- match.arg(method)

  cleanup <- .with_seed(seed)
  on.exit(cleanup(), add = TRUE)

  if (!inherits(points_sf, "sf")) stop("make_folds(): `points_sf` must be an sf object.")
  # Guard zero rows here rather than letting each method fail in its own way:
  # block_kfold reaches sf::st_bbox(), which returns all-NA for an empty layer
  # and makes st_as_sfc() abort with an opaque "!anyNA(x) is not TRUE".
  if (nrow(points_sf) == 0L)
    stop("make_folds(): `points_sf` has no rows; there is nothing to split into folds.",
         call. = FALSE)
  # MULTIPOINT must be coerced too, not merely admitted: st_coordinates()
  # returns one row per VERTEX, so any multi-vertex feature makes xy[i, ] a
  # different feature than row i, and every fold below misaligns silently.
  # coerce_to_points() takes centroids, matching prep_model_data().
  if (!all(sf::st_geometry_type(points_sf, by_geometry = TRUE) == "POINT"))
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
    row_ids <- points_sf$..row_id
    idx <- sample.int(n, n)
    sizes <- rep(floor(n / k), k)
    remainder <- n - sum(sizes)
    if (remainder > 0) sizes[seq_len(remainder)] <- sizes[seq_len(remainder)] + 1L
    splits <- vector("list", k); start <- 1L; assign_vec <- integer(n)
    for (j in seq_len(k)) {
      stop_i <- start + sizes[j] - 1L
      test_idx <- idx[start:stop_i]; train_idx <- setdiff(idx, test_idx)
      # Express splits in ..row_id values (not positions) so that
      # .remap_folds() / .cv_fit_one_fold() interpret them correctly even
      # when the input rows carry non-sequential IDs (e.g. after
      # prep_model_data() has dropped rows inside cv_gwr()/cv_bayes()).
      splits[[j]] <- list(train = row_ids[train_idx], test = row_ids[test_idx])
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
      b <- ensure_projected(boundary, .crs_or_null(pts))
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
                                      range_frac = range_frac,
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
    # One block means one fold whose training set is empty -- blocked CV
    # silently degenerating into nothing at all.  It happens whenever the block
    # size exceeds half the extent, which an accepted autocorrelation range can
    # do (estimate_sac_range() only rejects ranges above half the DIAGONAL, a
    # factor of sqrt(2) looser than what .block_dims_from_size() needs).
    if (B < 2L) {
      # block_size is NULL when the caller set block_nx/block_ny directly, and
      # sprintf() on a zero-length argument yields character(0) -- an empty
      # error message.
      how <- if (is.null(block_size))
        "the requested block_nx/block_ny" else
        sprintf("the block size (%s)", format(block_size))
      stop(sprintf(paste0("make_folds(block_kfold): %s produces a single block ",
                          "covering the whole extent, so there is no spatial ",
                          "split to make and the one fold would have an empty ",
                          "training set. Pass a smaller `block_size` (or more ",
                          "blocks), or set auto_range = FALSE if the size came ",
                          "from the estimated autocorrelation range."), how),
           call. = FALSE)
    }
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
      .log_warn("make_folds(block_kfold): fold size imbalance \u2014 largest fold has %d obs vs %d in smallest.",
                max(fold_loads), min(fold_loads))
    }

    splits <- vector("list", k); assign_vec <- integer(nrow(pts))
    row_ids <- pts$..row_id
    for (j in seq_len(k)) {
      test_idx <- which(pts$..block_id %in% fold_blocks[[j]])
      # Row-ID (not positional) splits — see comment in the random_kfold branch.
      splits[[j]] <- list(train = row_ids[setdiff(seq_len(nrow(pts)), test_idx)],
                          test  = row_ids[test_idx])
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
    row_ids <- pts$..row_id
    for (i in seq_len(n)) {
      excl <- sort(unique(nb[[i]]))
      # Row-ID (not positional) splits — see comment in the random_kfold branch.
      splits[[i]] <- list(train = row_ids[setdiff(seq_len(n), excl)],
                          test  = row_ids[i])
    }

    # A buffer that spans the data leaves EVERY fold with an empty training
    # set.  Those folds used to pass straight through .remap_folds() and get
    # dropped one at a time inside .cv_fit_one_fold(), so the only symptom was
    # a generic "all folds failed" warning at the very end of the CV run --
    # long after the cause could be acted on.  Say it here instead.
    n_train_each <- vapply(splits, function(s) length(s$train), integer(1))
    if (!any(n_train_each >= 2L))
      stop(sprintf(paste0("make_folds(buffered_loo): a buffer of %s excludes ",
                          "so much of the data that no fold retains 2 training ",
                          "points (largest training set: %d of %d). The buffer ",
                          "spans the data; use a smaller one, or 'block_kfold'."),
                   format(buffer), max(n_train_each), n), call. = FALSE)

    return(.ret(method, n, splits,
                .safe_tibble(row_id = pts$..row_id, fold = seq_len(n)),
                list(buffer = buffer)))
  }

  # ---- LEAVE-LOCATION-OUT (grouped) ----
  if (method == "leave_location_out") {
    if (is.null(group_var) || !(group_var %in% names(points_sf)))
      stop("make_folds(leave_location_out): `group_var` must name a column of ",
           "`points_sf` identifying the location each observation belongs to.",
           call. = FALSE)

    grp <- as.character(sf::st_drop_geometry(points_sf)[[group_var]])
    if (anyNA(grp))
      stop("make_folds(leave_location_out): `group_var` contains NA; every ",
           "observation must belong to a location.", call. = FALSE)

    row_ids <- points_sf$..row_id
    ug <- unique(grp)
    n_g <- length(ug)
    if (n_g < 2L)
      stop("make_folds(leave_location_out): need at least 2 distinct groups.",
           call. = FALSE)

    # k groups per fold; k >= n_g degenerates to leave-one-group-out.
    if (k > n_g) {
      .log_warn("make_folds(leave_location_out): k = %d exceeds the %d distinct groups; using leave-one-group-out.", k, n_g)
      k <- n_g
    }
    if (k < 2L) k <- 2L

    grp_fold <- sample(rep(seq_len(k), length.out = n_g))
    names(grp_fold) <- ug
    fold_of_row <- grp_fold[grp]

    splits <- vector("list", k)
    for (j in seq_len(k)) {
      test_rows <- row_ids[fold_of_row == j]
      splits[[j]] <- list(train = row_ids[fold_of_row != j], test = test_rows)
    }

    return(.ret(method, k, splits,
                .safe_tibble(row_id = row_ids, fold = unname(fold_of_row)),
                list(seed = seed, group_var = group_var, n_groups = n_g)))
  }

  # ---- NNDM: nearest-neighbour distance matching ----
  if (method == "nndm") {
    pts <- ensure_projected(points_sf)
    n <- nrow(pts)
    if (n > 5000L)
      stop(sprintf("make_folds(nndm): n = %d exceeds the safety threshold of 5000. Fold construction sorts distances from every point to every other (O(n^2) time), and NNDM then produces n leave-one-out folds, so the model is refitted n times. Use 'block_kfold' instead, or subset your data.", n),
           call. = FALSE)
    if (n < 3L)
      stop("make_folds(nndm): need at least 3 points.", call. = FALSE)

    if (is.null(prediction_points))
      stop("make_folds(nndm): `prediction_points` is required -- NNDM matches ",
           "the training-to-test distance distribution to the distance ",
           "distribution from your actual prediction locations to the training ",
           "data. Use the grid you intend to predict onto (see ",
           "predict_surface()).", call. = FALSE)
    pred <- ensure_projected(prediction_points, target_crs = .crs_or_null(pts))

    # Target: distance from each prediction location to its nearest training
    # point.  This is the distance regime the model will actually face, and the
    # one an arbitrary fixed buffer has no reason to reproduce.
    g_target <- .nn_dist_to(pred, pts)
    g_target <- g_target[is.finite(g_target)]
    if (!length(g_target))
      stop("make_folds(nndm): could not compute prediction-to-training ",
           "distances.", call. = FALSE)

    # Give each held-out point a buffer drawn from the target distribution, so
    # the realised test-to-train nearest-neighbour distances reproduce it by
    # construction.  Exclusion uses st_is_within_distance (GEOS STRtree), not a
    # dense n x n matrix -- at the 20000-point guard that would be over 3 GB.
    # Re-seed, but do NOT register a second restore handler.  The .with_seed()
    # at the top of make_folds() already owns save/restore for the whole
    # function; a nested one is wrong twice over.  Reusing the name `cleanup`
    # made both on.exit() expressions resolve to the inner closure, so the
    # caller's state was never restored -- and merely renaming it is not enough
    # either, because on.exit(add = TRUE) APPENDS: the outer handler would
    # restore the caller's seed and the inner one would then immediately
    # overwrite it with the intermediate post-set.seed(seed) state.
    #
    # All of which holds only for a NON-NULL seed.  With the default seed =
    # NULL, .with_seed(NULL) returns a no-op that neither saves nor restores,
    # so there is no outer handler to lean on -- and set.seed(NULL) does not
    # mean "leave the RNG alone", it RE-INITIALISES the generator from the
    # clock and process ID.  The ordinary call make_folds(pts, k, "nndm",
    # prediction_points = grid) therefore used to destroy the caller's RNG
    # state, silently.  Unseeded means unseeded: leave the stream alone and let
    # runif() advance it the way any other unseeded function would.
    if (!is.null(seed)) set.seed(seed)
    radii <- stats::quantile(g_target,
                             probs = stats::runif(n), names = FALSE, type = 7)

    # For each held-out point, exclude the nearest training points up to the
    # one whose distance is CLOSEST to that point's target radius -- rather
    # than excluding everything inside a hard radius.
    #
    # The hard-radius version is what a naive reading of "match the
    # distribution" suggests, but it systematically overshoots on clustered
    # data: excluding everything within r of a point in a tight cluster removes
    # the whole cluster, so the nearest survivor is in the NEXT cluster, far
    # beyond r.  On a two-cluster test that produced a realised median of 596
    # against a target of 193.  Choosing the order statistic nearest the target
    # instead matches as closely as the training configuration allows, and
    # degrades gracefully when no achievable distance is near the target.
    xy <- sf::st_coordinates(pts)[, 1:2, drop = FALSE]
    row_ids <- pts$..row_id
    splits <- vector("list", n)
    n_excluded <- integer(n)
    realised  <- rep(NA_real_, n)

    for (i in seq_len(n)) {
      d <- sqrt((xy[, 1] - xy[i, 1])^2 + (xy[, 2] - xy[i, 2])^2)
      d[i] <- Inf
      o  <- order(d)
      ds <- d[o]
      keep_max <- n - 1L                      # all but the held-out point

      # index of the achievable distance closest to this point's target
      j <- which.min(abs(ds[seq_len(keep_max)] - radii[i]))
      # never strip the training set below two points
      if (j > keep_max - 1L) j <- max(1L, keep_max - 1L)

      drop_idx <- if (j > 1L) o[seq_len(j - 1L)] else integer(0)
      train_i  <- setdiff(seq_len(n), c(i, drop_idx))
      n_excluded[i] <- length(drop_idx)
      realised[i]   <- ds[j]
      splits[[i]] <- list(train = row_ids[train_i], test = row_ids[i])
    }

    # Exclusion is not always needed.  When prediction locations sit no further
    # from the training data than training points sit from each other, plain
    # LOO already reproduces the target distribution and nothing is removed --
    # that is the correct answer, not a failure.  Conversely NNDM cannot pull
    # training points closer, so it cannot match a target that is shorter than
    # the training nearest-neighbour distances.
    .log_info(paste0("make_folds(nndm): %d folds. Target prediction-to-training ",
                     "distance: median %.1f. Realised test-to-training distance: ",
                     "median %.1f. Median buffer %.1f, excluding a median of %.1f ",
                     "training point(s) per fold."),
              n, stats::median(g_target), stats::median(realised),
              stats::median(radii), stats::median(n_excluded))

    return(.ret(method, n, splits,
                .safe_tibble(row_id = row_ids, fold = seq_len(n)),
                list(seed = seed, n_prediction_points = nrow(pred),
                     median_buffer   = stats::median(radii),
                     median_excluded = stats::median(n_excluded),
                     target_median   = stats::median(g_target),
                     realised_median = stats::median(realised),
                     target_distances   = g_target,
                     realised_distances = realised)))
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
#' @return A list with \code{overall}, \code{fold_metrics},
#'   \code{predictions}, \code{folds}, \code{n_folds_attempted},
#'   \code{n_folds_succeeded}, \code{formula} and \code{adaptive}.  The two
#'   fold counts make a run where every fold failed visible in the return
#'   value rather than only in a warning, since \code{overall} is a
#'   well-formed all-\code{NA} row either way.
#' @family cross-validation
#' @examples
#' \donttest{
#' if (requireNamespace("GWmodel", quietly = TRUE) &&
#'     requireNamespace("sp", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 60
#'   dat <- st_as_sf(
#'     data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), elev = rnorm(n)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   dat$price <- 10 + 0.01 * st_coordinates(dat)[, 1] + 2 * dat$elev + rnorm(n)
#'   cv <- cv_gwr(dat, "price", "elev", k = 3, bandwidth = 30)
#'   cv$overall
#'   cv$fold_metrics
#' }
#' }
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
    p = NULL, parallel = parallel, seed = seed
  )

  preds <- if (length(res$pred_rows)) do.call(rbind, res$pred_rows) else
    data.frame(`..row_id` = integer(), fold = integer(),
               y = numeric(), yhat = numeric(), y_train_mean = numeric())
  folds_df <- if (length(res$fold_stats)) as.data.frame(dplyr::bind_rows(res$fold_stats)) else
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
       # Reported so that a run where every fold failed is visible in the
       # return value, not only in a warning: `overall` is a well-formed
       # all-NA row either way.  cv_spatial() and cv_rf() report the same two
       # fields, so all four CV entry points share one shape.
       n_folds_attempted = n_attempted, n_folds_succeeded = n_succeeded,
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
#'   A user-supplied \code{gp_k} is respected in every fold; when omitted, the
#'   GP rank is auto-selected per training fold.  \code{compute_loo},
#'   \code{boundary}, and \code{pointize} are always overridden by the CV
#'   internals.
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
#' @return A list with \code{overall}, \code{fold_metrics},
#'   \code{predictions}, \code{folds}, \code{n_folds_attempted},
#'   \code{n_folds_succeeded}, \code{formula} and \code{predictive_coverage}.
#'   The two fold counts make a run where every fold failed visible in the
#'   return value rather than only in a warning.
#'   The \code{predictive_coverage} entries (one per
#'   \code{coverage_levels} value, plus \code{mean_CRPS}) are averages across
#'   folds \strong{weighted by each fold's \code{n_pred}}, because the per-fold
#'   values in \code{fold_metrics} are themselves means over that fold's test
#'   rows; an unweighted average would not be the pooled quantity when fold
#'   sizes differ, which for spatially blocked folds they routinely do.
#' @family cross-validation
#' @examples
#' \dontrun{
#' if (requireNamespace("brms", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 60
#'   dat <- st_as_sf(
#'     data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), elev = rnorm(n)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   dat$price <- 10 + 0.01 * st_coordinates(dat)[, 1] + 2 * dat$elev + rnorm(n)
#'   cv <- cv_bayes(dat, "price", "elev", k = 2,
#'                  fit_args = list(chains = 2, iter = 500))
#'   cv$overall
#'   cv$predictive_coverage  # coverage at 50/80/95% plus mean CRPS
#' }
#' }
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
    # NOTE: a user-supplied fit_args$gp_k is respected as-is.  When absent,
    # fit_bayesian_spatial_model()'s default (gp_k = NULL) auto-selects a
    # rank per training fold.  (Previously `gp_k = NULL` was passed through
    # modifyList(), which *removes* NULL elements and therefore silently
    # deleted the user's gp_k.)
    fold_fit_args <- modifyList(fit_args, list(
      compute_loo = FALSE, boundary = boundary,
      pointize = pointize,
      .already_prepped = TRUE
    ))
    do.call(fit_bayesian_spatial_model,
            c(list(data_sf = train_sf, response_var = response_var,
                   predictor_vars = predictor_vars), fold_fit_args))
  }

  # Extra per-fold info: gp_k, CRPS, predictive coverage
  fold_info_fn <- function(fit_obj, test_sf, y_true, y_hat) {
    extras <- list(
      gp_k       = as.integer(fit_obj$info$gp_k %||% NA_integer_),
      gp_n_basis = as.integer(fit_obj$info$gp_n_basis %||% NA_integer_),
      n_draws    = NA_integer_,
      CRPS    = NA_real_
    )
    # Pre-initialise coverage columns so every fold emits the same schema
    # even when the posterior-draw step fails for some folds; heterogeneous
    # per-fold columns would otherwise break row-binding of fold_metrics.
    for (cl in coverage_levels) {
      extras[[sprintf("coverage_%.0f", cl * 100)]] <- NA_real_
    }

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
        # observations — see .crps_energy() for the formula and reference.
        crps_per_obs <- .crps_energy(ppred_draws, y_true)
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
    p = NULL, parallel = parallel, seed = seed
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
  # bind_rows() (rather than rbind) tolerates folds whose extras differ,
  # filling missing columns with NA instead of erroring.
  folds_df <- if (length(res$fold_stats)) as.data.frame(dplyr::bind_rows(res$fold_stats)) else
    data.frame(fold = integer(), n_train = integer(), n_test = integer(),
               n_pred = integer(),
               RMSE = numeric(), MAE = numeric(), MAPE = numeric(),
               SMAPE = numeric(), R2 = numeric(), Adj_R2 = numeric(),
               CRPS = numeric(), gp_k = integer(),
               gp_n_basis = integer(), n_draws = integer())

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
       # Reported so that a run where every fold failed is visible in the
       # return value, not only in a warning -- which matters most here,
       # where a fold dies on Stan sampling rather than on bad input.
       n_folds_attempted = n_attempted, n_folds_succeeded = n_succeeded,
       formula = deparse(stats::reformulate(predictor_vars, response_var)),
       # Weighted by each fold's n_pred.  Per-fold coverage and CRPS are
       # themselves means over that fold's test rows, so an unweighted mean
       # across folds is NOT the pooled quantity as soon as fold sizes differ
       # -- and block_kfold tolerates a 3:1 imbalance before it even warns, so
       # they routinely do.  (The `overall` metrics above are already pooled
       # correctly, from the stacked prediction rows; this was the odd one out.)
       predictive_coverage = if (nrow(folds_df) > 0L && "CRPS" %in% names(folds_df)) {
         w_all <- if (!is.null(folds_df$n_pred))
           suppressWarnings(as.numeric(folds_df$n_pred)) else rep(1, nrow(folds_df))
         wm <- function(v) {
           v  <- suppressWarnings(as.numeric(v))
           ok <- is.finite(v) & is.finite(w_all) & w_all > 0
           if (!any(ok)) return(NA_real_)
           sum(v[ok] * w_all[ok]) / sum(w_all[ok])
         }
         cov_cols  <- grep("^coverage_", names(folds_df), value = TRUE)
         cov_means <- vapply(cov_cols, function(cn) wm(folds_df[[cn]]), numeric(1))
         c(as.list(cov_means), mean_CRPS = wm(folds_df$CRPS))
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
#' @return A list with \code{overall}, \code{fold_metrics}, \code{predictions},
#'   \code{folds}, and the two fold counts \code{n_folds_attempted} and
#'   \code{n_folds_succeeded}.  The counts are reported deliberately: a
#'   \code{fit_fn} that fails on every fold otherwise looks like a successful
#'   run that happened to score \code{NA}, so compare them before trusting
#'   \code{overall}.  The \code{fold} column of \code{fold_metrics} and
#'   \code{predictions} carries the fold's index in the \code{folds} object
#'   that was supplied, so it lines up with
#'   \code{make_folds()$assignment$fold} even when some folds were unusable
#'   and dropped.
#' @family cross-validation
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
    parallel = parallel, seed = seed
  )

  preds <- if (length(res$pred_rows)) do.call(rbind, res$pred_rows) else
    data.frame(`..row_id` = integer(), fold = integer(),
               y = numeric(), yhat = numeric(), y_train_mean = numeric())
  folds_df <- if (length(res$fold_stats)) as.data.frame(dplyr::bind_rows(res$fold_stats)) else
    data.frame()

  # cv_gwr() and cv_bayes() both raise a real condition here; cv_spatial() used
  # to return an all-NA `overall` and an empty data.frame with nothing at R
  # condition level, so a fit_fn that failed on every fold looked like a
  # successful run that happened to score NA.
  n_attempted <- length(remapped_folds)
  n_succeeded <- length(res$fold_stats)
  if (n_succeeded == 0L && n_attempted > 0L) {
    .log_warn("cv_spatial(): all %d folds failed to produce predictions; results are empty.",
              n_attempted)
    warning("cv_spatial(): all folds failed; cross-validation results contain ",
            "no predictions.", call. = FALSE)
  } else if (n_succeeded < n_attempted) {
    .log_warn("cv_spatial(): %d of %d folds produced predictions.",
              n_succeeded, n_attempted)
  }

  list(overall = .cv_overall_metrics(preds),
       fold_metrics = folds_df, predictions = preds,
       folds = remapped_folds,
       n_folds_attempted = n_attempted, n_folds_succeeded = n_succeeded)
}
