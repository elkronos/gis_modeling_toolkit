# ---------------------------------------------------------------------------
# Constructor
# ---------------------------------------------------------------------------

#' Build a spatial_fit S3 object
#'
#' Low-level constructor used by \code{fit_gwr_model()},
#' \code{fit_bayesian_spatial_model()} and \code{fit_rf_model()}.  Users should
#' not call this directly.
#'
#' @section The coef() contract:
#' \code{coef()} on a \code{spatial_fit} either returns the coefficients or
#' signals an error -- it never returns \code{NULL}.  \code{coef.rf_fit()}
#' always errors, because a forest has no coefficients; \code{coef.gwr_fit()}
#' and \code{coef.bayesian_fit()} error when the backend cannot supply them
#' (a missing package, an engine without the expected component).  A
#' \code{NULL} return would be indistinguishable from "this model genuinely
#' has no fixed effects", so \code{lapply(fits, coef)} would quietly produce a
#' shorter answer than the caller expected.  Wrap in \code{try()} or
#' \code{tryCatch()} when sweeping over a heterogeneous list of fits.
#'
#' @param subclass Character scalar: "gwr_fit", "bayesian_fit" or "rf_fit".
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

#' Print a fitted spatial model
#'
#' @param x A \code{spatial_fit} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly (called for its side effect).
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
    cat(sprintf("  GP basis: %d per dim (%d total)\n",
                x$info$gp_k %||% NA_integer_,
                x$info$gp_n_basis %||% NA_integer_))
    if (is.finite(x$info$looic %||% NA_real_))
      cat(sprintf("  LOOIC   : %.2f\n", x$info$looic))
    if (!isTRUE(x$info$convergence_ok))
      cat("  ** Convergence warnings present \u2014 see $info$convergence_diagnostics\n")
  }
  invisible(x)
}


#' Summarise a fitted spatial model
#'
#' Computes goodness-of-fit metrics from \code{fitted(object)} against the
#' observed response.
#'
#' @section What the metrics are computed on:
#' For a \code{gwr_fit} or a \code{bayesian_fit} these are \strong{in-sample}
#' metrics: \code{fitted()} returns values computed at the training locations
#' from the model that saw them.  For an \code{rf_fit} they are
#' \strong{out-of-bag}, because \code{fitted.rf_fit()} returns out-of-bag
#' predictions rather than in-sample ones (see \code{\link{fit_rf_model}}).
#' The two are not comparable, and \code{print()} on the result labels which
#' one it is holding, driven by \code{$info$fitted_are_oob}.  Use
#' \code{\link{compare_models_cv}} to compare backends.
#'
#' Adjusted R-squared is suppressed: GWR's effective parameter count far
#' exceeds the global predictor count, and a GP model has no simple \code{p}.
#'
#' @param object A \code{spatial_fit} object.
#' @param ... Ignored.
#' @return An object of class \code{summary.spatial_fit}: a list with
#'   \code{class}, \code{formula}, \code{n}, \code{response_var},
#'   \code{predictor_vars}, \code{info} and \code{in_sample} (the metric
#'   data.frame, out-of-bag for an \code{rf_fit}).
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
  # An rf_fit's fitted() values are out-of-bag, not in-sample.  Labelling them
  # "In-sample" understates a random hold-out as a memorised fit, which is the
  # opposite of the truth and the opposite of the caveat ?fit_rf_model carries.
  if (isTRUE(x$info$fitted_are_oob))
    cat("\n  Out-of-bag metrics (NOT in-sample; see ?fit_rf_model):\n")
  else
    cat("\n  In-sample metrics:\n")
  m <- x$in_sample
  cat(sprintf("    RMSE  = %.4f\n", m$RMSE))
  cat(sprintf("    MAE   = %.4f\n", m$MAE))
  cat(sprintf("    R\u00b2    = %.4f\n", m$R2))
  if (is.finite(m$Adj_R2 %||% NA_real_))
    cat(sprintf("    Adj R\u00b2= %.4f\n", m$Adj_R2))
  if (is.finite(m$SMAPE %||% NA_real_))
    cat(sprintf("    SMAPE = %.2f%%\n", m$SMAPE))
  invisible(x)
}

# ---------------------------------------------------------------------------
# model_metrics generic
# ---------------------------------------------------------------------------

#' Compute goodness-of-fit metrics for a spatial model
#'
#' @section What the metrics are computed on:
#' With \code{newdata = NULL} the metrics come from \code{fitted(object)}.
#' That is \strong{in-sample} for a \code{gwr_fit} or a \code{bayesian_fit},
#' but \strong{out-of-bag} for an \code{rf_fit}, whose \code{fitted()} method
#' returns out-of-bag predictions (see \code{\link{fit_rf_model}}).  The
#' returned data.frame carries no label distinguishing the two, so check
#' \code{object$info$fitted_are_oob} before comparing numbers across backends
#' -- or use \code{\link{compare_models_cv}}, which scores every backend the
#' same way.
#'
#' @param object A \code{spatial_fit} object.
#' @param newdata Optional sf object for out-of-sample evaluation. If NULL,
#'   fitted values are used (in-sample, or out-of-bag for an \code{rf_fit};
#'   see above).
#' @param ... Additional arguments passed to predict().
#' @return A data.frame with n, RMSE, MAE, MAPE, SMAPE, R2, Adj_R2.
#' @export
model_metrics <- function(object, ...) UseMethod("model_metrics")


#' @rdname model_metrics
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

#' Expand a prediction vector back to the original newdata length
#'
#' Places predictions into the positions indicated by \code{clean_idx},
#' filling remaining positions with NA.
#'
#' \code{clean_idx} is derived from the \code{..orig_row_id..} sentinel column
#' every predict method attaches before calling \code{prep_model_data()}, so
#' there is exactly one implementation of "which rows survived cleaning".
#'
#' @param preds Numeric vector of predictions on the clean subset.  Its length
#'   must equal \code{sum(clean_idx)}; a mismatch is an error, never a recycle.
#' @param clean_idx Logical vector marking the rows that survived cleaning.
#' @param n_orig Integer original number of rows.
#' @return Numeric vector of length \code{n_orig}.
#' @keywords internal
#' @noRd
.expand_predictions <- function(preds, clean_idx, n_orig) {
  n_clean <- sum(clean_idx)
  # Without this check `out[which(clean_idx)] <- preds` RECYCLES silently: 4
  # clean slots and a length-2 `preds` fills 1 2 . 1 . 2 and returns it as a
  # prediction.  The equal-length fast path below would likewise hand back a
  # vector of the wrong length unchecked.
  if (length(preds) != n_clean)
    stop(sprintf(paste0(".expand_predictions(): %d prediction value(s) for %d ",
                        "row(s) that survived cleaning; the two cannot be ",
                        "aligned."),
                 length(preds), n_clean), call. = FALSE)
  if (n_clean == n_orig) return(preds)
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
#' \code{newdata} is first transformed to the CRS used during fitting
#' (via \code{ensure_projected()}), so predictions are computed in a
#' single coordinate system regardless of the CRS newdata arrives in.
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

  # Align newdata to the CRS used during fitting BEFORE prep_model_data().
  # prep's own ensure_projected() call has no target and would leave
  # already-projected newdata in *its* CRS (or auto-pick a UTM zone for
  # lon/lat input independent of training), after which gwr.predict()
  # would silently mix coordinates from two different systems.  This
  # mirrors predict.bayesian_fit().
  newdata <- ensure_projected(newdata, target_crs = .crs_or_null(object$data_sf))

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

  # Degrade to an all-NA vector when nothing survived cleaning, matching
  # predict.rf_fit() and predict.bayesian_fit().  Without this the .to_sp()
  # call below -- which sits outside the tryCatch -- would surface a raw
  # sf-to-Spatial coercion error instead.
  if (n_new == 0L) {
    .log_warn(paste0("predict.gwr_fit(): every row of newdata was dropped as ",
                     "missing or non-finite; returning %d NA(s)."), n_orig)
    return(rep(NA_real_, n_orig))
  }

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
      # Silently skipping a predictor here would hand brms an un-centred,
      # un-scaled column while the model was fitted on a standardised one --
      # wrong predictions with no diagnostic.  A numeric column that came back
      # as character from a CSV round-trip is the common way to get here.
      if (!(pv %in% names(pred_df)))
        stop(sprintf(paste0("predict.bayesian_fit(): predictor '%s' was ",
                            "standardised at fit time but is absent from ",
                            "newdata; the same standardisation cannot be ",
                            "applied."), pv), call. = FALSE)
      if (!is.numeric(pred_df[[pv]]))
        stop(sprintf(paste0("predict.bayesian_fit(): predictor '%s' was ",
                            "standardised at fit time but is %s in newdata, ",
                            "not numeric. Convert it back to numeric (a CSV ",
                            "round-trip is the usual cause)."),
                     pv, class(pred_df[[pv]])[1L]), call. = FALSE)
      pred_df[[pv]] <- (pred_df[[pv]] - ps[[pv]]$center) / ps[[pv]]$scale
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
#' @return Numeric vector of length \code{nrow(newdata)}, or a
#'   \code{n_draws x nrow(newdata)} matrix when \code{draws = TRUE} (a 1-row
#'   all-\code{NA} matrix if the posterior draw fails).  With
#'   \code{newdata = NULL} the cached \code{fitted()} values are returned only
#'   for the default \code{summary = "mean"}, \code{type = "epred"},
#'   \code{draws = FALSE} combination; any other combination is recomputed
#'   against the training data, because the cache holds epred column means and
#'   nothing else.
#' @export
predict.bayesian_fit <- function(object, newdata = NULL,
                                 summary = c("mean", "median"),
                                 type = c("epred", "predict"),
                                 draws = FALSE, ...) {
  summary <- match.arg(summary)
  type    <- match.arg(type)

  if (is.null(newdata)) {
    # fitted() caches colMeans(posterior_epred(...)) and nothing else, so it
    # answers exactly one argument combination.  Short-circuiting for any
    # other one returned means where medians were asked for, and expected
    # values where observation noise was asked for -- silently.
    if (!draws && summary == "mean" && type == "epred") return(fitted(object))
    newdata <- object$data_sf
  }

  # ---- Early validation: ensure all predictor columns are present ----
  # Done BEFORE the brms availability check: this is a structural check on
  # newdata that needs no backend, and it gives the most informative error.
  missing_preds <- setdiff(object$predictor_vars, names(newdata))
  if (length(missing_preds) > 0L)
    stop(sprintf(
      "predict.bayesian_fit(): newdata is missing required predictor column(s): %s",
      paste(missing_preds, collapse = ", ")
    ), call. = FALSE)

  if (!requireNamespace("brms", quietly = TRUE))
    stop("predict.bayesian_fit(): package 'brms' is required.", call. = FALSE)

  model_obj <- object$engine
  if (!inherits(model_obj, "brmsfit"))
    stop("predict.bayesian_fit(): engine is not a brmsfit object.", call. = FALSE)

  # ---- Preprocessing: match the pipeline used during fitting ----
  # Ensure newdata is in the same projected CRS that was used for training,
  # BEFORE prep_model_data(): prep's own ensure_projected() has no target and
  # would leave already-projected newdata in *its* CRS.  This mirrors
  # predict.gwr_fit() and predict.rf_fit().
  n_orig       <- nrow(newdata)
  training_crs <- sf::st_crs(object$data_sf)
  newdata      <- ensure_projected(newdata, target_crs = .crs_or_null(training_crs))

  # Sentinel row ID: after prep_model_data()'s single-pass cleaning (which also
  # coerces non-POINT geometry, MULTIPOINT included, so st_coordinates() in
  # .prepare_brms_pred_df() stays one row per observation) we can derive
  # exactly which original rows survived -- no duplicated filtering logic.
  newdata$..orig_row_id.. <- seq_len(n_orig)
  newdata <- prep_model_data(
    newdata, object$response_var, object$predictor_vars,
    pointize = "auto", require_response = FALSE
  )
  clean_idx <- seq_len(n_orig) %in% newdata$..orig_row_id..
  newdata$..orig_row_id.. <- NULL

  if (nrow(newdata) == 0L) {
    .log_warn(paste0("predict.bayesian_fit(): every row of newdata was dropped ",
                     "as missing or non-finite; returning %d NA(s)."), n_orig)
    return(if (draws) matrix(NA_real_, nrow = 1L, ncol = n_orig)
           else rep(NA_real_, n_orig))
  }

  # Build prediction data.frame with scaled coordinates & standardised predictors
  pred_df <- .prepare_brms_pred_df(object, newdata)

  # Draw from posterior
  draw_fn <- if (type == "epred") brms::posterior_epred else brms::posterior_predict
  draw_mat <- try(draw_fn(model_obj, newdata = pred_df), silent = TRUE)

  if (inherits(draw_mat, "try-error") || !is.matrix(draw_mat)) {
    .log_warn("predict.bayesian_fit(): posterior draw failed.")
    # Honour the documented return shape: a matrix when draws = TRUE, so a
    # caller that indexes columns is not handed a vector on the failure path.
    return(if (draws) matrix(NA_real_, nrow = 1L, ncol = n_orig)
           else rep(NA_real_, n_orig))
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

#' In-sample fitted values from a GWR fit
#'
#' Reads the fitted values out of the GWmodel result, which stores them under
#' one of several names depending on version and entry point; the extraction
#' falls back through the local coefficients and the residuals when no direct
#' column is present.  These are \strong{in-sample} values -- each observation
#' was inside its own bandwidth window -- so \code{summary()} on a
#' \code{gwr_fit} reports an optimistic fit.  Use \code{\link{cv_gwr}} for a
#' spatially blocked estimate.
#'
#' @param object A \code{gwr_fit}.
#' @param ... Ignored.
#' @return Numeric vector of length \code{object$n} (\code{NA} where extraction
#'   failed).
#' @export
fitted.gwr_fit <- function(object, ...) {
  .extract_gwr_values(
    object$engine, object$data_sf, object$formula,
    object$n, object$response_var
  )
}


#' In-sample fitted values from a Bayesian spatial GP fit
#'
#' Posterior expectation at the training locations: the column means of
#' \code{brms::posterior_epred()}.  These are \strong{in-sample} values.
#'
#' @section The result is cached:
#' \code{posterior_epred()} is O(draws x n), and \code{summary()},
#' \code{residuals()}, \code{model_metrics()} and \code{compare_models()} each
#' call \code{fitted()} independently, so the value is memoised in an
#' environment carried in \code{object$info$.cache} (reference semantics, so it
#' survives R's copy-on-modify).  The cache holds epred column means only,
#' which is why \code{predict(object, summary = "median")} and
#' \code{predict(object, type = "predict")} recompute rather than reuse it.
#' Call \code{\link{clear_fitted_cache}} if the engine or the training data has
#' been mutated by hand after fitting.
#'
#' @param object A \code{bayesian_fit}.
#' @param ... Ignored.
#' @return Numeric vector of length \code{object$n} (all \code{NA} if the
#'   posterior draw failed).
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

#' In-sample residuals from a GWR fit
#'
#' Observed response minus \code{\link{fitted.gwr_fit}}, so these are
#' \strong{in-sample} residuals.
#'
#' @param object A \code{gwr_fit}.
#' @param ... Ignored.
#' @return Numeric vector of length \code{object$n}.
#' @export
residuals.gwr_fit <- function(object, ...) {
  y_obs <- sf::st_drop_geometry(object$data_sf)[[object$response_var]]
  y_obs - fitted(object)
}


#' In-sample residuals from a Bayesian spatial GP fit
#'
#' Observed response minus \code{\link{fitted.bayesian_fit}} (the cached
#' posterior expectation), so these are \strong{in-sample} residuals.
#'
#' @param object A \code{bayesian_fit}.
#' @param ... Ignored.
#' @return Numeric vector of length \code{object$n}.
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
#'   Never \code{NULL}: when the engine carries no \code{SDF} component this
#'   errors, following the \code{coef()} contract described in
#'   \code{\link{new_spatial_fit}}.
#' @export
coef.gwr_fit <- function(object, ...) {
  sdf <- tryCatch(object$engine$SDF, error = function(e) NULL)
  # Returning NULL here would be indistinguishable from "this GWR has no
  # coefficients", which is never true.  See ?new_spatial_fit for the contract.
  if (is.null(sdf))
    stop("coef.gwr_fit(): the fit carries no GWmodel `SDF` component, so the ",
         "local coefficients cannot be extracted.", call. = FALSE)
  if (inherits(sdf, "Spatial")) sdf@data else sf::st_drop_geometry(sdf)
}


#' Extract Bayesian model fixed-effect summaries
#'
#' @param object A \code{bayesian_fit} object.
#' @param ... Ignored.
#' @return A matrix of fixed-effect posterior summaries, as returned by
#'   \code{brms::fixef()}.  Never \code{NULL}: a missing 'brms' or a failing
#'   \code{fixef()} call errors, following the \code{coef()} contract described
#'   in \code{\link{new_spatial_fit}}.
#' @export
coef.bayesian_fit <- function(object, ...) {
  # A NULL return would be indistinguishable from "this model has no fixed
  # effects", while predict()/fitted() on the same object stop() informatively.
  # See ?new_spatial_fit for the contract.
  if (!requireNamespace("brms", quietly = TRUE))
    stop("coef.bayesian_fit(): package 'brms' is required.", call. = FALSE)
  fx <- tryCatch(brms::fixef(object$engine), error = function(e) e)
  if (inherits(fx, "error"))
    stop(sprintf("coef.bayesian_fit(): brms::fixef() failed: %s",
                 conditionMessage(fx)), call. = FALSE)
  fx
}
