# ---------------------------------------------------------------------------
# Constructor
# ---------------------------------------------------------------------------

#' Build a spatial_fit S3 object
#'
#' Low-level constructor used by \code{fit_gwr_model()} and
#' \code{fit_bayesian_spatial_model()}.  Users should not call this directly.
#'
#' @param subclass Character scalar: "gwr_fit" or "bayesian_fit".
#' @param engine   The raw model object.
#' @param formula  A formula.
#' @param response_var  Character(1).
#' @param predictor_vars Character vector.
#' @param data_sf  An sf object used for fitting.
#' @param info     Named list of model-specific extras.
#' @return An object of class \code{c(subclass, "spatial_fit")}.
#' @export
new_spatial_fit <- function(subclass, engine, formula, response_var,
                            predictor_vars, data_sf, info = list()) {
  stopifnot(is.character(subclass), length(subclass) == 1L)
  # Use an environment for the cache so it has reference semantics —

  # mutations persist across calls without triggering R copy-on-modify.
  if (is.null(info$.cache)) info$.cache <- new.env(parent = emptyenv())
  obj <- list(
    engine         = engine,
    formula        = formula,
    response_var   = response_var,
    predictor_vars = predictor_vars,
    data_sf        = data_sf,
    n              = nrow(data_sf),
    info           = info
  )
  class(obj) <- c(subclass, "spatial_fit")
  obj
}

# ---------------------------------------------------------------------------
# print / summary
# ---------------------------------------------------------------------------

#' @export
print.spatial_fit <- function(x, ...) {
  subclass <- class(x)[1L]
  label <- switch(subclass,
    gwr_fit      = "GWR (GWmodel)",
    bayesian_fit = "Bayesian Spatial GP (brms)",
    subclass
  )
  cat(sprintf("<%s> spatial model fit\n", label))
  cat(sprintf("  Formula : %s\n", deparse(x$formula)))
  cat(sprintf("  n       : %d\n", x$n))

  if (subclass == "gwr_fit") {
    cat(sprintf("  Bandwidth: %.4g (%s, %s kernel)\n",
                x$info$bandwidth,
                if (isTRUE(x$info$adaptive)) "adaptive" else "fixed",
                x$info$kernel %||% "bisquare"))
    if (is.finite(x$info$AICc %||% NA_real_))
      cat(sprintf("  AICc    : %.2f\n", x$info$AICc))
  }
  if (subclass == "bayesian_fit") {
    cat(sprintf("  GP rank : %d\n", x$info$gp_k %||% NA_integer_))
    if (is.finite(x$info$looic %||% NA_real_))
      cat(sprintf("  LOOIC   : %.2f\n", x$info$looic))
    if (!isTRUE(x$info$convergence_ok))
      cat("  ** Convergence warnings present — see $info$convergence_diagnostics\n")
  }
  invisible(x)
}


#' @export
summary.spatial_fit <- function(object, ...) {
  fit_vals <- fitted(object)
  y_obs    <- sf::st_drop_geometry(object$data_sf)[[object$response_var]]
  # GWR fits locally varying coefficients at every observation, so the
  # global predictor count is not a valid effective-parameter count for
  # Adj R².  Pass p = NULL to suppress it, consistent with cv_gwr().
  # Bayesian fits likewise use p = NULL.
  met      <- .compute_reg_metrics(y_obs, fit_vals, p = NULL)

  out <- list(
    class          = class(object)[1L],
    formula        = object$formula,
    n              = object$n,
    response_var   = object$response_var,
    predictor_vars = object$predictor_vars,
    info           = object$info,
    in_sample      = met
  )
  class(out) <- "summary.spatial_fit"
  out
}


#' @export
print.summary.spatial_fit <- function(x, ...) {
  cat(sprintf("Summary of <%s> fit (n = %d)\n\n", x$class, x$n))
  cat(sprintf("  Formula: %s\n", deparse(x$formula)))
  cat("\n  In-sample metrics:\n")
  m <- x$in_sample
  cat(sprintf("    RMSE  = %.4f\n", m$RMSE))
  cat(sprintf("    MAE   = %.4f\n", m$MAE))
  cat(sprintf("    R²    = %.4f\n", m$R2))
  if (is.finite(m$Adj_R2 %||% NA_real_))
    cat(sprintf("    Adj R²= %.4f\n", m$Adj_R2))
  if (is.finite(m$SMAPE %||% NA_real_))
    cat(sprintf("    SMAPE = %.2f%%\n", m$SMAPE))
  invisible(x)
}

# ---------------------------------------------------------------------------
# model_metrics generic
# ---------------------------------------------------------------------------

#' Compute goodness-of-fit metrics for a spatial model
#'
#' @param object A \code{spatial_fit} object.
#' @param newdata Optional sf object for out-of-sample evaluation. If NULL,
#'   in-sample (fitted) values are used.
#' @param ... Additional arguments passed to predict().
#' @return A data.frame with n, RMSE, MAE, MAPE, SMAPE, R2, Adj_R2.
#' @export
model_metrics <- function(object, ...) UseMethod("model_metrics")


#' @export
model_metrics.spatial_fit <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    y_hat <- fitted(object)
    y_obs <- sf::st_drop_geometry(object$data_sf)[[object$response_var]]
  } else {
    if (!(object$response_var %in% names(newdata)))
      stop("model_metrics(): 'newdata' must contain the response variable '",
           object$response_var,
           "' to compute evaluation metrics. predict() can be used without it.",
           call. = FALSE)
    y_hat <- predict(object, newdata = newdata, ...)
    y_obs <- sf::st_drop_geometry(newdata)[[object$response_var]]
  }
  # Adj R² is suppressed (p = NULL) because GWR's effective parameter count
  # far exceeds the global predictor count, and Bayesian GP models likewise
  # lack a simple p.  This is consistent with the CV evaluation path.
  .compute_reg_metrics(y_obs, y_hat, p = NULL)
}

# ---------------------------------------------------------------------------
# NA-safe prediction helpers
# ---------------------------------------------------------------------------

#' Identify rows that would survive prep_model_data() cleaning
#'
#' Mirrors the complete-cases + finite-check logic in prep_model_data() so
#' that predict methods can track which rows are dropped and later expand
#' the result back to the original length with NA fill.
#'
#' @param data_sf An sf object.
#' @param response_var Character(1) response column name.
#' @param predictor_vars Character vector of predictor column names.
#' @param require_response Logical; when FALSE (prediction mode) the response
#'   column is excluded from the completeness check so that NAs in a reused
#'   dataset do not silently drop valid rows.  Default TRUE.
#' @return Logical vector of length \code{nrow(data_sf)}; TRUE = row is clean.
#' @keywords internal
#' @noRd
.clean_row_mask <- function(data_sf, response_var, predictor_vars,
                            require_response = TRUE) {
  check_vars <- if (require_response) c(response_var, predictor_vars) else predictor_vars
  req_cols <- intersect(check_vars, names(data_sf))
  df <- sf::st_drop_geometry(data_sf)[, req_cols, drop = FALSE]
  ok_cc <- stats::complete.cases(df)
  num_mask <- vapply(df, is.numeric, logical(1))
  ok_fin <- if (any(num_mask)) {
    apply(as.matrix(df[, num_mask, drop = FALSE]), 1L,
          function(r) all(is.finite(r)))
  } else {
    rep(TRUE, nrow(df))
  }
  ok_cc & ok_fin
}


#' Expand a prediction vector back to the original newdata length
#'
#' Places predictions into the positions indicated by \code{clean_idx},
#' filling remaining positions with NA.
#'
#' @param preds Numeric vector of predictions on the clean subset.
#' @param clean_idx Logical vector from \code{.clean_row_mask()}.
#' @param n_orig Integer original number of rows.
#' @return Numeric vector of length \code{n_orig}.
#' @keywords internal
#' @noRd
.expand_predictions <- function(preds, clean_idx, n_orig) {
  if (sum(clean_idx) == n_orig) return(preds)
  out <- rep(NA_real_, n_orig)
  out[which(clean_idx)] <- preds
  out
}


# ---------------------------------------------------------------------------
# predict — GWR
# ---------------------------------------------------------------------------

#' Predict from a GWR spatial model
#'
#' When \code{newdata} is NULL, returns the in-sample fitted values.
#' Otherwise uses \code{GWmodel::gwr.predict()} on the new locations.
#'
#' @param object A \code{gwr_fit} object.
#' @param newdata An sf object with the same predictors.
#'   The response variable need not be present (true out-of-sample prediction
#'   is supported).  NULL = fitted values.
#' @param ... Ignored.
#' @return Numeric vector of predictions.
#' @export
predict.gwr_fit <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(fitted(object))

  if (!requireNamespace("GWmodel", quietly = TRUE))
    stop("predict.gwr_fit(): package 'GWmodel' is required.", call. = FALSE)
  if (!requireNamespace("sp", quietly = TRUE))
    stop("predict.gwr_fit(): package 'sp' is required.", call. = FALSE)

  # Record original length so we can return an aligned vector with NA fill.
  n_orig <- nrow(newdata)

  # Use a sentinel row ID column so that after prep_model_data()'s single-pass

  # cleaning we can derive exactly which original rows survived — no duplicated
  # filtering logic.
  newdata$..orig_row_id.. <- seq_len(n_orig)

  newdata <- prep_model_data(
    newdata, object$response_var, object$predictor_vars,
    pointize = "auto", require_response = FALSE
  )
  n_new <- nrow(newdata)

  clean_idx <- seq_len(n_orig) %in% newdata$..orig_row_id..
  newdata$..orig_row_id.. <- NULL

  # .to_sp() uses intersect() internally, so it gracefully handles newdata
  # that lacks the response column (true out-of-sample prediction).
  needed_cols <- unique(c(object$response_var, object$predictor_vars))
  sp_train <- .to_sp(object$data_sf, needed_cols)
  sp_new   <- .to_sp(newdata, needed_cols)

  bw <- object$info$bandwidth
  if (isTRUE(object$info$adaptive)) bw <- as.integer(round(bw))

  pred_obj <- tryCatch(
    suppressWarnings(
      GWmodel::gwr.predict(
        object$formula, data = sp_train, predictdata = sp_new,
        bw = bw, kernel = object$info$kernel %||% "bisquare",
        adaptive = object$info$adaptive %||% TRUE
      )
    ),
    error = function(e) {
      .log_warn("predict.gwr_fit(): gwr.predict() failed: %s", conditionMessage(e))
      NULL
    }
  )

  preds_clean <- if (is.null(pred_obj)) {
    rep(NA_real_, n_new)
  } else {
    .extract_gwr_values(pred_obj, newdata, object$formula, n_new,
                         object$response_var)
  }

  # Expand back to original length, filling dropped rows with NA.
  .expand_predictions(preds_clean, clean_idx, n_orig)
}


# ---------------------------------------------------------------------------
# predict — Bayesian
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Private helper: build the prediction data.frame expected by brms from an
# sf object, applying the coordinate scaling and predictor standardisation
# stored in a bayesian_fit's $info slot.
# ---------------------------------------------------------------------------
.prepare_brms_pred_df <- function(object, data_sf) {
  pred_df <- sf::st_drop_geometry(data_sf)
  crd     <- sf::st_coordinates(data_sf)

  cs <- object$info$coord_scaling
  if (!is.null(cs)) {
    pred_df$`..x` <- (crd[, 1] - cs$x_center) / cs$x_scale
    pred_df$`..y` <- (crd[, 2] - cs$y_center) / cs$y_scale
  } else {
    pred_df$`..x` <- crd[, 1]
    pred_df$`..y` <- crd[, 2]
  }

  ps <- object$info$predictor_scaling
  if (!is.null(ps) && length(ps) > 0L) {
    for (pv in names(ps)) {
      if (pv %in% names(pred_df) && is.numeric(pred_df[[pv]])) {
        pred_df[[pv]] <- (pred_df[[pv]] - ps[[pv]]$center) / ps[[pv]]$scale
      }
    }
  }

  pred_df
}


#' Predict from a Bayesian spatial GP model
#'
#' Applies the same newdata preparation pipeline as \code{predict.gwr_fit()}:
#' non-point geometries are coerced to points, the data is projected to the
#' CRS used during fitting (via \code{ensure_projected()}), and rows with
#' missing or non-finite values are dropped.  Coordinate scaling and predictor
#' standardisation stored at fit time are then applied before delegating to
#' \code{brms::posterior_epred()} or \code{brms::posterior_predict()}.
#'
#' @param object A \code{bayesian_fit} object.
#' @param newdata An sf object with the same predictors.  The response variable
#'   need not be present (true out-of-sample prediction is supported).
#'   NULL = fitted values.
#' @param summary "mean" (default) or "median" over posterior draws.
#' @param type "epred" (default) for expected predictions (no obs noise),
#'   or "predict" for full posterior predictive draws (includes obs noise).
#' @param draws If TRUE, return the full posterior draw matrix instead of a
#'   point summary.  Default FALSE.
#' @param ... Ignored.
#' @return Numeric vector (or matrix when draws=TRUE).
#' @export
predict.bayesian_fit <- function(object, newdata = NULL,
                                 summary = c("mean", "median"),
                                 type = c("epred", "predict"),
                                 draws = FALSE, ...) {
  summary <- match.arg(summary)
  type    <- match.arg(type)

  if (is.null(newdata)) {
    if (draws) {
      # Return fitted draw matrix from training data
      newdata <- object$data_sf
    } else {
      return(fitted(object))
    }
  }

  if (!requireNamespace("brms", quietly = TRUE))
    stop("predict.bayesian_fit(): package 'brms' is required.", call. = FALSE)

  model_obj <- object$engine
  if (!inherits(model_obj, "brmsfit"))
    stop("predict.bayesian_fit(): engine is not a brmsfit object.", call. = FALSE)

  # ---- Early validation: ensure all predictor columns are present ----
  missing_preds <- setdiff(object$predictor_vars, names(newdata))
  if (length(missing_preds) > 0L)
    stop(sprintf(
      "predict.bayesian_fit(): newdata is missing required predictor column(s): %s",
      paste(missing_preds, collapse = ", ")
    ), call. = FALSE)

  # ---- Preprocessing: match the pipeline used during fitting ----
  # Coerce non-point geometry to points (same as prep_model_data()).
  if (!all(sf::st_geometry_type(newdata, by_geometry = TRUE) %in%
           c("POINT", "MULTIPOINT"))) {
    newdata <- coerce_to_points(newdata, "auto")
  }

  # Ensure newdata is in the same projected CRS that was used for training.
  # The training data stored in `object$data_sf` is already projected by
  # prep_model_data(), so we use its CRS as the target.
  training_crs <- sf::st_crs(object$data_sf)
  newdata <- ensure_projected(newdata, target_crs = training_crs)

  # Track which rows are clean so we can return an aligned result.
  n_orig    <- nrow(newdata)
  clean_idx <- .clean_row_mask(newdata, object$response_var,
                               object$predictor_vars,
                               require_response = FALSE)
  if (sum(clean_idx) < n_orig) {
    .log_warn("predict.bayesian_fit(): dropping %d row(s) with non-finite or missing values.",
              n_orig - sum(clean_idx))
    newdata <- newdata[clean_idx, , drop = FALSE]
  }

  # Build prediction data.frame with scaled coordinates & standardised predictors
  pred_df <- .prepare_brms_pred_df(object, newdata)

  # Draw from posterior
  draw_fn <- if (type == "epred") brms::posterior_epred else brms::posterior_predict
  draw_mat <- try(draw_fn(model_obj, newdata = pred_df), silent = TRUE)

  if (inherits(draw_mat, "try-error") || !is.matrix(draw_mat)) {
    .log_warn("predict.bayesian_fit(): posterior draw failed.")
    return(rep(NA_real_, n_orig))
  }

  if (draws) {
    # Expand draw matrix: insert NA columns for dropped rows.
    if (sum(clean_idx) < n_orig) {
      full_mat <- matrix(NA_real_, nrow = nrow(draw_mat), ncol = n_orig)
      full_mat[, which(clean_idx)] <- draw_mat
      return(full_mat)
    }
    return(draw_mat)
  }

  preds_clean <- if (summary == "mean") colMeans(draw_mat) else apply(draw_mat, 2L, stats::median)

  # Expand back to original length, filling dropped rows with NA.
  .expand_predictions(preds_clean, clean_idx, n_orig)
}


# ---------------------------------------------------------------------------
# fitted
# ---------------------------------------------------------------------------

#' @export
fitted.gwr_fit <- function(object, ...) {
  .extract_gwr_values(
    object$engine, object$data_sf, object$formula,
    object$n, object$response_var
  )
}


#' @export
fitted.bayesian_fit <- function(object, ...) {
  # --- Lazy cache: posterior_epred() is O(draws × n) and expensive.
  #     summary(), residuals(), model_metrics(), and compare_models() all
  #     call fitted() independently.  Caching avoids redundant passes.
  #     The cache lives in an environment (reference semantics) so it

  #     persists even though R lists are copy-on-modify.
  cache <- object$info$.cache
  if (!is.null(cache) && exists(".fitted_values", envir = cache, inherits = FALSE)) {
    return(get(".fitted_values", envir = cache, inherits = FALSE))
  }

  if (!requireNamespace("brms", quietly = TRUE))
    stop("fitted.bayesian_fit(): package 'brms' is required.", call. = FALSE)

  model_obj <- object$engine
  pred_df   <- .prepare_brms_pred_df(object, object$data_sf)

  draws <- try(brms::posterior_epred(model_obj, newdata = pred_df), silent = TRUE)
  if (inherits(draws, "try-error") || !is.matrix(draws))
    return(rep(NA_real_, object$n))

  fitted_vals <- colMeans(draws)

  # Store in cache for subsequent calls
  if (!is.null(cache)) {
    assign(".fitted_values", fitted_vals, envir = cache)
  }

  fitted_vals
}


#' Clear cached fitted values for a Bayesian spatial model
#'
#' Removes the lazily-cached \code{fitted()} result so that the next call
#' recomputes from the posterior.  This is only necessary if the underlying
#' \code{brmsfit} engine or training data has been manually mutated after
#' fitting — normal usage never requires it.
#'
#' @param object A \code{bayesian_fit} object.
#' @return \code{object}, invisibly (called for side effect).
#' @export
clear_fitted_cache <- function(object) {
  cache <- object$info$.cache
  if (!is.null(cache) && exists(".fitted_values", envir = cache, inherits = FALSE)) {
    rm(".fitted_values", envir = cache)
  }
  invisible(object)
}


# ---------------------------------------------------------------------------
# residuals
# ---------------------------------------------------------------------------

#' @export
residuals.gwr_fit <- function(object, ...) {
  y_obs <- sf::st_drop_geometry(object$data_sf)[[object$response_var]]
  y_obs - fitted(object)
}


#' @export
residuals.bayesian_fit <- function(object, ...) {
  y_obs <- sf::st_drop_geometry(object$data_sf)[[object$response_var]]
  y_obs - fitted(object)
}


# ---------------------------------------------------------------------------
# coef
# ---------------------------------------------------------------------------

#' Extract GWR local coefficients
#'
#' @param object A \code{gwr_fit} object.
#' @param ... Ignored.
#' @return A data.frame of local coefficient estimates (one row per obs).
#' @export
coef.gwr_fit <- function(object, ...) {
  sdf <- tryCatch(object$engine$SDF, error = function(e) NULL)
  if (is.null(sdf)) return(NULL)
  if (inherits(sdf, "Spatial")) sdf@data else sf::st_drop_geometry(sdf)
}


#' Extract Bayesian model fixed-effect summaries
#'
#' @param object A \code{bayesian_fit} object.
#' @param ... Ignored.
#' @return A data.frame of fixed-effect posterior summaries.
#' @export
coef.bayesian_fit <- function(object, ...) {
  if (!requireNamespace("brms", quietly = TRUE)) return(NULL)
  tryCatch(brms::fixef(object$engine), error = function(e) NULL)
}
