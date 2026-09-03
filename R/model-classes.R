# ---------------------------------------------------------------------------
# Constructor
# ---------------------------------------------------------------------------

#' Build a spatial_fit S3 object
#'
#' The constructor for the \code{spatial_fit} class, and the public entry point
#' for plugging your own model backend into this package.  The three built-in
#' fitters -- \code{\link{fit_gwr_model}()},
#' \code{\link{fit_bayesian_spatial_model}()} and \code{\link{fit_rf_model}()}
#' -- all end by calling it, and so should a custom \code{fit_fn} written for
#' \code{\link{cv_spatial}()}: wrapping your model in a \code{spatial_fit} is
#' what lets it use the package's folds, metrics, comparison and
#' area-of-applicability machinery unchanged.
#'
#' There are two obligations.  Return an object built here from your
#' \code{fit_fn}, and define a \code{predict()} method for the \code{subclass}
#' you chose -- \code{\link{cv_spatial}()} scores folds by calling the
#' \code{predict()} generic on the fit, so without a matching
#' \code{predict.<subclass>()} every fold fails.  Methods for
#' \code{\link{fitted}()}, \code{\link{residuals}()} and \code{\link{coef}()}
#' are optional; supply them if you want the corresponding helpers to work on
#' your fits too.
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
#' @param subclass Character scalar naming the class to stamp on the object:
#'   one of the built-ins \code{"gwr_fit"}, \code{"bayesian_fit"} or
#'   \code{"rf_fit"}, or any name of your own for a custom backend (say
#'   \code{"lm_fit"}).  It is the S3 dispatch key: \code{predict()},
#'   \code{fitted()}, \code{residuals()} and \code{coef()} on the result all
#'   dispatch on it, so a custom \code{subclass} \strong{requires} a matching
#'   \code{predict.<subclass>()} method to be usable with
#'   \code{\link{cv_spatial}()}.
#' @param engine   The raw model object your backend produced (an \code{lm},
#'   a \code{ranger} object, a \code{brmsfit}, ...).  Nothing inspects it
#'   except your own methods.
#' @param formula  A formula.
#' @param response_var  Character(1).
#' @param predictor_vars Character vector.
#' @param data_sf  An sf object used for fitting.
#' @param info     Named list of model-specific extras.  Set
#'   \code{fitted_are_oob = TRUE} if your \code{fitted()} values are held out
#'   rather than in-sample, so \code{summary()} labels them honestly.
#' @return An object of class \code{c(subclass, "spatial_fit")}.
#' @seealso \code{\link{cv_spatial}()}, which consumes a custom \code{fit_fn};
#'   \code{\link{fit_rf_model}()} for a worked built-in fitter.
#' @family model fitting
#' @examples
#' library(sf)
#' set.seed(1)
#' n <- 80
#' site <- st_as_sf(
#'   data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), elev = rnorm(n)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' site$price <- 10 + 0.01 * st_coordinates(site)[, 1] + 2 * site$elev + rnorm(n)
#'
#' # A custom backend: an ordinary linear model behind the spatial_fit interface.
#' lm_fit <- function(train_sf) {
#'   new_spatial_fit(
#'     subclass       = "lm_fit",
#'     engine         = lm(price ~ elev, st_drop_geometry(train_sf)),
#'     formula        = price ~ elev,
#'     response_var   = "price",
#'     predictor_vars = "elev",
#'     data_sf        = train_sf
#'   )
#' }
#'
#' # Required: cv_spatial() scores each fold through the predict() generic,
#' # which dispatches on the subclass named above.
#' predict.lm_fit <- function(object, newdata = NULL, ...) {
#'   if (is.null(newdata)) newdata <- object$data_sf
#'   as.numeric(stats::predict(object$engine, st_drop_geometry(newdata)))
#' }
#' registerS3method("predict", "lm_fit", predict.lm_fit)
#'
#' cv <- cv_spatial(site, "price", "elev", fit_fn = lm_fit, k = 3, seed = 1)
#' cv$overall
#' # Always check these two agree before trusting the metrics above.
#' c(attempted = cv$n_folds_attempted, succeeded = cv$n_folds_succeeded)
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

#' Fetch fitted() from a spatial_fit and check it against the fit's own n
#'
#' \code{new_spatial_fit()} is the documented extension point, so a subclass
#' whose \code{fitted()} method is missing or returns the wrong length is
#' user-reachable -- and nothing downstream notices.  With no method at all,
#' \code{stats::fitted()} finds \code{object$fitted} (absent) and returns
#' \code{NULL}, so \code{.compute_reg_metrics()} reports \code{n = 0} and
#' all-\code{NA}; with a method returning 60 values for 120 rows, the metric
#' code drops the pairs it cannot align and reports the fit's \code{n = 120}
#' against all-\code{NA} numbers -- a plausible row count over a silently
#' mis-indexed comparison.  \code{.cv_run_folds()} guards exactly this on the
#' prediction side; these two paths did not.
#'
#' @param object A \code{spatial_fit}.
#' @param .caller Name used in the error message.
#' @return \code{fitted(object)}, guaranteed to be length \code{object$n}.
#' @keywords internal
#' @noRd
.fitted_checked <- function(object, .caller = "model_metrics") {
  fit_vals <- fitted(object)
  n_exp    <- object$n
  if (is.null(fit_vals) || length(fit_vals) != n_exp)
    stop(sprintf(paste0("%s(): fitted() returned %s for a fit of %d ",
                        "observation(s). Define a fitted.%s() method that ",
                        "returns one value per row of the fit's `data_sf`, in ",
                        "the same order (see ?new_spatial_fit)."),
                 .caller,
                 if (is.null(fit_vals)) "NULL"
                 else sprintf("%d value(s)", length(fit_vals)),
                 n_exp, class(object)[1L]),
         call. = FALSE)
  fit_vals
}


# ---------------------------------------------------------------------------
# print / summary
# ---------------------------------------------------------------------------

#' Print a fitted spatial model
#'
#' Shows the one-screen summary of any \code{spatial_fit}: backend, formula,
#' number of observations, CRS, and the few backend-specific numbers worth
#' seeing immediately (GWR bandwidth, GP basis size, forest settings).  It is
#' what you get by typing the object's name, and the quickest way to confirm a
#' fit used the data, predictors and CRS you meant.  For fit quality use
#' \code{\link{model_metrics}()} or \code{\link{summary}()} instead --
#' nothing printed here is an out-of-sample score.
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
  fit_vals <- .fitted_checked(object, .caller = "summary")
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
#' Reports RMSE, MAE, MAPE, SMAPE, \eqn{R^2} and adjusted \eqn{R^2} for any
#' \code{spatial_fit}, in one row and on one scale, so that fits from different
#' backends can be read side by side.  Reach for it to score a model on data you
#' hold out yourself (pass it as \code{newdata}), or to get a quick in-sample
#' reading of how closely a fit tracks its training data.
#'
#' It is not a substitute for cross-validation.  With \code{newdata = NULL} the
#' numbers are in-sample for a \code{gwr_fit} or \code{bayesian_fit} -- and a
#' GWR can reach a near-perfect in-sample \eqn{R^2} at a small bandwidth
#' without predicting anything.  For a figure you can report, use
#' \code{\link{cv_gwr}()}, \code{\link{cv_bayes}()}, \code{\link{cv_rf}()}
#' or \code{\link{compare_models_cv}()}.
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
    y_hat <- .fitted_checked(object, .caller = "model_metrics")
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
  # A response that is not numeric cannot be scored, and is.finite() on a
  # character column is FALSE everywhere -- so 100 perfectly good rows used
  # to come back as n = 0 with every metric NA, silently.
  if (is.logical(y_obs)) y_obs <- as.numeric(y_obs)
  if (!is.numeric(y_obs))
    stop(sprintf(paste0("model_metrics(): response '%s' is %s, not numeric, so ",
                        "no metric can be computed from it. Convert the column ",
                        "first."),
                 object$response_var, class(y_obs)[1L]), call. = FALSE)
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
  # CRS-less newdata first gets the interpretation the TRAINING data got
  # (recorded by ensure_projected() when it assumed lon/lat), so the same rows
  # land where they did at fit time -- see .replay_crs_assumption().
  newdata <- .replay_crs_assumption(newdata, object$data_sf, "predict.gwr_fit")
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
      # A real warning(), not only a logger line.  A logger line is invisible
      # to tryCatch(warning = ), to withCallingHandlers(), to
      # testthat::expect_warning() and to R CMD check, so a predict() that
      # returns nothing but NA left no trace a caller could act on.  The
      # commonest cause is a factor or character predictor: gwr.basic() expands
      # contrasts via model.matrix() and fits, gwr.predict() does not and fails
      # here.  fit_gwr_model() now rejects those at fit time, so reaching this
      # generally means the fit was built by other means.
      .log_warn("predict.gwr_fit(): gwr.predict() failed: %s", conditionMessage(e))
      warning(sprintf(paste0("predict.gwr_fit(): GWmodel::gwr.predict() ",
                             "failed, so every prediction is NA. Cause: %s"),
                      conditionMessage(e)), call. = FALSE)
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
  newdata      <- .replay_crs_assumption(newdata, object$data_sf, "predict.bayesian_fit")
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
#' Call \code{\link{clear_fitted_cache}} if the engine has been mutated by hand
#' after fitting.
#'
#' @section The cache is shared by copies, and validated:
#' An environment has reference semantics, which is what makes the memo survive
#' R's copy-on-modify -- but it also means \code{fit2 <- fit} gives the two
#' objects \emph{the same} cache.  Assigning a different \code{data_sf} to the
#' copy would then have returned the original's cached values, at the original's
#' length, which \code{residuals()} silently recycled against the copy's shorter
#' response.  The entry therefore carries the \code{n} and a digest of the
#' training data it was computed from, and is recomputed whenever either fails
#' to match, so a copy with different data recomputes instead of reading the
#' original's answer.
#'
#' Two consequences of the shared environment remain and cannot be removed from
#' here: \code{\link{clear_fitted_cache}} on one copy empties the cache both
#' share (harmless -- the other simply recomputes), and \code{identical()}
#' cannot distinguish two fits by their caches.  The digest covers
#' \code{data_sf} only, not \code{$engine}: a hand-mutated \code{brmsfit} is
#' what \code{\link{clear_fitted_cache}} is for.
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
  #     persists even though R lists are copy-on-modify -- and, for the same
  #     reason, is SHARED by every copy of the fit.  See the @section above:
  #     the entry is stamped with the n and a digest of the data it was
  #     computed from, and anything that does not match is recomputed.
  cache <- object$info$.cache
  key   <- .fitted_cache_key(object)
  if (!is.null(cache) && exists(".fitted_values", envir = cache, inherits = FALSE)) {
    hit <- get(".fitted_values", envir = cache, inherits = FALSE)
    if (is.list(hit) && identical(hit$n, object$n) &&
        identical(hit$key, key) &&
        is.numeric(hit$values) && length(hit$values) == object$n)
      return(hit$values)
    # Stale (a copy carrying different data, or an entry written by an older
    # version of this package).  Drop it rather than returning it.
    rm(list = ".fitted_values", envir = cache)
  }

  if (!requireNamespace("brms", quietly = TRUE))
    stop("fitted.bayesian_fit(): package 'brms' is required.", call. = FALSE)

  model_obj <- object$engine
  pred_df   <- .prepare_brms_pred_df(object, object$data_sf)

  draws <- try(brms::posterior_epred(model_obj, newdata = pred_df), silent = TRUE)
  # An error here is an error: a posterior that cannot be drawn used to come
  # back as all-NA fitted values with nothing said, and summary() and
  # model_metrics() then reported n = 0 as though the data were missing.
  # predict.bayesian_fit() has always raised; this matches it.
  if (inherits(draws, "try-error"))
    stop(sprintf(paste0("fitted.bayesian_fit(): brms::posterior_epred() failed ",
                        "on the training data: %s"),
                 conditionMessage(attr(draws, "condition"))), call. = FALSE)
  if (!is.matrix(draws) || ncol(draws) != object$n)
    stop(sprintf(paste0("fitted.bayesian_fit(): brms::posterior_epred() returned ",
                        "%s where a draws x %d matrix was expected."),
                 if (is.matrix(draws)) paste(dim(draws), collapse = " x ")
                 else class(draws)[1L], object$n), call. = FALSE)

  fitted_vals <- colMeans(draws)

  # Store in cache for subsequent calls, stamped so a copy carrying different
  # data cannot read it back.  A wrong-length result is never cached.
  if (!is.null(cache) && length(fitted_vals) == object$n) {
    assign(".fitted_values",
           list(n = object$n, key = key, values = fitted_vals),
           envir = cache)
  }

  fitted_vals
}


#' Cheap fingerprint of the training data a cached fitted() was computed from
#'
#' Covers exactly what \code{.prepare_brms_pred_df()} reads: the attribute
#' columns and the coordinates.  \code{$engine} is deliberately excluded --
#' digesting a \code{brmsfit} with all its draws would cost more than the
#' posterior pass the cache exists to avoid, and a hand-mutated engine is what
#' \code{\link{clear_fitted_cache}} is for.  Returns \code{NA_character_} if a
#' digest cannot be taken, which still fingerprints consistently (a stored
#' \code{NA} matches a computed \code{NA}), leaving the \code{n} and length
#' checks as the guard.
#'
#' @param object A \code{spatial_fit}.
#' @return A length-one character.
#' @keywords internal
#' @noRd
.fitted_cache_key <- function(object) {
  tryCatch(
    digest::digest(list(
      n      = object$n,
      data   = sf::st_drop_geometry(object$data_sf),
      coords = unname(sf::st_coordinates(object$data_sf))
    )),
    error = function(e) NA_character_
  )
}


#' Clear cached fitted values for a Bayesian spatial model
#'
#' Removes the lazily-cached \code{fitted()} result so that the next call
#' recomputes from the posterior.  This is only necessary if the underlying
#' \code{brmsfit} engine has been manually mutated after fitting -- a change to
#' \code{data_sf} invalidates the entry on its own, because the cached value
#' carries a digest of the data it was computed from (see
#' \code{\link{fitted.bayesian_fit}}).  Normal usage never requires it.
#'
#' The cache environment is shared by every copy of a fit, so clearing it
#' through one copy clears it for all of them.  That is harmless: the others
#' recompute.
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
#' Returns the whole surface of coefficients -- one row per observation, one
#' column per term -- rather than the single global vector \code{coef()}
#' returns for an \code{lm}.  That table is the point of fitting a GWR at all:
#' inspect the spread of a predictor's column to see where, and by how much,
#' its relationship with the response changes across the study area, and join
#' it back to \code{object$data_sf} to map it.  Use
#' \code{\link{plot.spatial_fit}()} for a quick look at that map.
#'
#' @section What is and is not returned:
#' Only the model terms -- the intercept and one column per predictor.
#' GWmodel's \code{SDF} data slot carries a good deal more alongside them
#' (standard errors, t-values, the observed response, the fitted values, the
#' residuals, \code{Local_R2}): 15 columns for a two-predictor fit, of which 3
#' are coefficients.  Returning the whole slot would have made
#' \code{coef(fit)$Local_R2} and \code{coef(fit)$a_SE} read like coefficients
#' and \code{ncol(coef(fit))} a meaningless number.  Reach for
#' \code{object$engine$SDF} when you want the rest; it is the unmodified
#' GWmodel object.
#'
#' If the model terms cannot be located in the \code{SDF} -- a GWmodel that
#' names its coefficient columns differently -- the whole slot is returned with
#' a warning saying so, rather than an error or a silently short table.
#'
#' @param object A \code{gwr_fit} object.
#' @param ... Ignored.
#' @return A data.frame of local coefficient estimates: one row per
#'   observation, one column per model term.
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
  dat <- if (inherits(sdf, "Spatial")) sdf@data else sf::st_drop_geometry(sdf)

  # Keep only the model terms.  gwr.basic()'s SDF is
  #   c(colnames(betas), "y", "yhat", "residual", "CV_Score", "Stud_residual",
  #     paste0(colnames(betas), "_SE"), paste0(colnames(betas), "_TV"),
  #     "Local_R2")
  # so the coefficients are the leading block and everything after it is a
  # diagnostic.  Match on the model-matrix column names (normalised, so
  # "(Intercept)" and GWmodel's "Intercept" agree) and take the FIRST unclaimed
  # SDF column for each -- coefficients come first, so first-match is the
  # coefficient even when a predictor shares a name with a diagnostic column.
  want <- tryCatch(
    colnames(stats::model.matrix(object$formula,
                                 data = sf::st_drop_geometry(object$data_sf))),
    error = function(e) NULL
  )
  .norm  <- function(x) tolower(gsub("[^[:alnum:]_.]", "", x))
  if (!is.null(want) && length(want) > 0L) {
    sdf_n  <- .norm(names(dat))
    pos    <- integer(0)
    for (w in .norm(want)) {
      cand <- setdiff(which(sdf_n == w), pos)
      if (length(cand) > 0L) pos <- c(pos, cand[[1L]])
    }
    if (length(pos) == length(want)) return(dat[, pos, drop = FALSE])
    .log_warn(paste0("coef.gwr_fit(): matched only %d of %d model term(s) to a ",
                     "column of GWmodel's SDF, so the full SDF data slot is ",
                     "returned instead of the coefficients alone. Columns such ",
                     "as *_SE, *_TV, residual and Local_R2 are diagnostics, ",
                     "not coefficients."),
              length(pos), length(want))
  }
  dat
}


#' Extract Bayesian model fixed-effect summaries
#'
#' Returns the posterior summary of the global (non-spatial) regression terms:
#' estimate, error and credible interval per predictor, as
#' \code{brms::fixef()} reports them.  Reach for it to read the average effect
#' of a predictor with its uncertainty attached -- the Bayesian counterpart to
#' a coefficient table -- remembering that the Gaussian-process term has
#' already absorbed the spatially structured part of the signal, so these are
#' effects net of location.
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
