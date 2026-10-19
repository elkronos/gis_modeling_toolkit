# ---------------------------------------------------------------------------
# Constructor
# ---------------------------------------------------------------------------

#' Build a spatial_fit S3 object
#'
#' The constructor for the \code{spatial_fit} class, and the public entry point
#' for plugging your own model backend into this package.  The three built-in
#' fitters (\code{\link{fit_gwr_model}()},
#' \code{\link{fit_bayesian_spatial_model}()} and \code{\link{fit_rf_model}()})
#' all end by calling it, and so should a custom \code{fit_fn} written for
#' \code{\link{cv_spatial}()}: wrapping your model in a \code{spatial_fit} is
#' what lets it use the package's folds, metrics, comparison and
#' area-of-applicability machinery unchanged.
#'
#' There are two obligations.  Return an object built here from your
#' \code{fit_fn}, and define a \code{predict()} method for the \code{subclass}
#' you chose.  \code{\link{cv_spatial}()} scores folds by calling the
#' \code{predict()} generic on the fit, so without a matching
#' \code{predict.<subclass>()} every fold fails.  A
#' \code{fitted.<subclass>()} method returning one value per row of the fit's
#' \code{data_sf}, in the same order, is \strong{required} by
#' \code{\link{summary}()} and \code{\link{model_metrics}()}: both error,
#' naming the method to define, without one.  \code{\link{residuals}()} and
#' \code{\link{coef}()} methods are optional.
#'
#' @section The coef() contract:
#' \code{coef()} on one of the three built-in backends either returns the
#' coefficients or signals an error.  It never returns \code{NULL}.  A custom
#' subclass inherits \code{stats::coef.default()}, which returns \code{NULL},
#' so define a \code{coef.<subclass>()} that errors when your backend has no
#' coefficients; otherwise the hazard described below applies to your own fits.
#' \code{coef.rf_fit()}
#' always errors, because a forest has no coefficients; \code{coef.gwr_fit()}
#' and \code{coef.bayesian_fit()} error when the backend cannot supply them
#' (a missing package, an engine without the expected component).  A
#' \code{NULL} return would be indistinguishable from "this model has no
#' fixed effects", so \code{lapply(fits, coef)} would quietly produce a
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
#'   \code{fitted_are_oob = TRUE} when your \code{fitted()} values are held
#'   out instead of in-sample, so \code{summary()} labels them correctly.
#' @return An object of class \code{c(subclass, "spatial_fit")}.
#' @seealso \code{\link{cv_spatial}()}, which consumes a custom \code{fit_fn};
#'   \code{\link{fit_rf_model}()} for a worked built-in fitter.
#' @family model fitting
#' @examples
#' library(sf)
#' set.seed(1)
#' n <- 80
#' site <- st_as_sf(
#'   data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
#'              elev = rnorm(n)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' site$price <- 10 + 0.01 * (st_coordinates(site)[, 1] - 5e5) +
#'   2 * site$elev + rnorm(n)
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
  # ALWAYS a fresh one: an `info` list handed in by the caller can carry
  # another fit's `.cache` (summary() used to hand one out), and two fits
  # sharing a cache is exactly the condition under which one can read the
  # other's fitted values back.  Reusing it saved nothing -- the entry is
  # keyed anyway -- and cost correctness.
  info$.cache <- new.env(parent = emptyenv())
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
#' user-reachable, and nothing downstream notices.  With no method at all,
#' \code{stats::fitted()} finds \code{object$fitted} (absent) and returns
#' \code{NULL}, so \code{.compute_reg_metrics()} reports \code{n = 0} and
#' all-\code{NA}; with a method returning 60 values for 120 rows, the metric
#' code drops the pairs it cannot align and reports the fit's \code{n = 120}
#' against all-\code{NA} numbers, a plausible row count over a silently
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
#' Shows the one-screen summary of a \code{gwr_fit}, a \code{bayesian_fit} or a
#' custom subclass: backend, formula, number of observations, CRS, and the few
#' backend-specific numbers worth seeing immediately (GWR bandwidth; for a
#' \code{bayesian_fit} the GP basis size, the smallest effective range that
#' basis resolves in map units, the slope prior, and a flagged line when the
#' sampler's diagnostics or the basis-adequacy check failed).  An
#' \code{rf_fit} has its own method (see
#' \code{\link{print.rf_fit}}), which shows the same header plus the forest
#' settings.  It is what you get by typing the object's name, and the quickest
#' way to confirm a fit used the data, predictors and CRS you meant.  For fit
#' quality use \code{\link{model_metrics}()} or \code{\link{summary}()}
#' instead.  Nothing printed here is an out-of-sample score.
#'
#' For a \code{gwr_fit} the bandwidth and AICc are followed by the effective
#' number of parameters and, when any local window is of deficient rank, by
#' the number of such windows.
#'
#' @param x A \code{spatial_fit} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly (called for its side effect).
#' @family methods on a fitted model
#' @export
print.spatial_fit <- function(x, ...) {
  subclass <- class(x)[1L]
  # The parenthesis names the package that computes the fit.  A GWR is solved
  # by this package (see ?fit_gwr_model, "The local solve"), so it has none;
  # "GWR (GWmodel)" was the header while GWmodel::gwr.basic() made the fit.
  label <- switch(subclass,
    gwr_fit      = "GWR",
    bayesian_fit = "Bayesian Spatial GP (brms)",
    subclass
  )
  cat(sprintf("<%s> spatial model fit\n", label))
  # deparse(), not deparse1(): a long formula comes back as SEVERAL strings and
  # sprintf() vectorises over them, printing "Formula :" once per line with the
  # continuation indented as if it were a new field.  Every bayesian_fit hits
  # this, because the gp() term this package builds is always long.
  cat(sprintf("  Formula : %s\n", paste(deparse(x$formula), collapse = " ")))
  cat(sprintf("  n       : %d\n", x$n))
  # The CRS is what ?print.spatial_fit promises and never printed, and it is
  # the one thing a user cannot otherwise read off the object -- getting it
  # wrong silently changes every distance in the model.
  cat(sprintf("  CRS     : %s\n", .fold_crs_label(x$data_sf)))

  if (subclass == "gwr_fit") {
    # "%.4g" printed a fixed bandwidth of 122372 m as "1.224e+05", rounded and
    # without the unit the help page tells the reader to check.  A count of
    # neighbours for adaptive, a distance with the CRS's unit for fixed.
    bw <- suppressWarnings(as.numeric(x$info$bandwidth %||% NA_real_))[1L]
    bw_txt <- if (!is.finite(bw)) "unknown"
      else if (isTRUE(x$info$adaptive)) sprintf("%s neighbours", format(bw))
      else {
        u <- tryCatch(sf::st_crs(x$data_sf)$units_gdal, error = function(e) NULL)
        sprintf("%s %s", format(signif(bw, 6), big.mark = ",", scientific = FALSE),
                if (length(u) != 1L || is.na(u) || !nzchar(u)) "CRS units" else u)
      }
    cat(sprintf("  Bandwidth: %s (%s, %s kernel)\n", bw_txt,
                if (isTRUE(x$info$adaptive)) "adaptive" else "fixed",
                x$info$kernel %||% "bisquare"))
    if (is.finite(x$info$AICc %||% NA_real_))
      cat(sprintf("  AICc    : %.2f\n", x$info$AICc))
    # The effective number of parameters, 2 tr(S) - tr(S'S): what the local
    # regressions cost in degrees of freedom, against the handful of global
    # coefficients the formula shows.
    enp <- suppressWarnings(as.numeric(x$info$enp %||% NA_real_))[1L]
    if (is.finite(enp))
      cat(sprintf("  Eff. parameters: %.1f of %d observations\n", enp, x$n))
    # A count the coefficient table needs read beside it: windows that could
    # not estimate every coefficient return NA for those they could not.
    n_def <- suppressWarnings(as.integer(x$info$n_rank_deficient %||% 0L))[1L]
    if (isTRUE(n_def > 0L))
      cat(sprintf(paste0("  Rank-deficient windows: %d of %d (their ",
                         "unestimable coefficients are NA; see ",
                         "$info$nonfinite_coef)\n"), n_def, x$n))
  }
  if (subclass == "bayesian_fit") {
    cat(sprintf("  GP basis: %d per dim (%d total)\n",
                x$info$gp_k %||% NA_integer_,
                x$info$gp_n_basis %||% NA_integer_))
    for (ln in .bayes_gp_print_lines(x$info)) cat(ln, "\n", sep = "")
    if (is.finite(x$info$looic %||% NA_real_))
      cat(sprintf("  LOOIC   : %.2f\n", x$info$looic))
    # The slope prior the model used.  A standardised fit sets one, and it is
    # part of what the coefficients mean.
    bp <- x$info$b_prior
    if (is.character(bp) && length(bp) == 1L && !is.na(bp))
      cat(sprintf("  Slope prior: %s%s\n", bp,
                  if (is.character(x$info$b_prior_source) &&
                      !is.na(x$info$b_prior_source))
                    sprintf(" (%s)", x$info$b_prior_source) else ""))
    # The coefficients of a standardised fit are per SD of each predictor, and
    # nothing coef() returns says so.
    if (length(x$info$predictor_scaling) > 0L)
      cat(sprintf(paste0("  Predictors standardised: %s (coef() is per SD; ",
                         "see $info$predictor_scaling)\n"),
                  paste(names(x$info$predictor_scaling), collapse = ", ")))
    # NA is "not checked" (check_convergence = FALSE, or no diagnostic could be
    # read), which is neither a pass nor a failure; NULL (nothing recorded)
    # still reads as a failure.
    cv_ok <- x$info$convergence_ok
    if (length(cv_ok) == 1L && is.na(cv_ok))
      cat("  Convergence: NOT CHECKED (fitted with check_convergence = FALSE?)\n")
    else if (!isTRUE(cv_ok))
      cat("  ** Convergence warnings present -- see $info$convergence_diagnostics\n")
    else {
      # Passed, with a caution: R-hat between 1.01 and 1.05 on parameters
      # that have enough effective draws.  Not flagged, and not hidden.
      n_caution <- length(x$info$convergence_diagnostics$rhat_caution)
      if (n_caution > 0L)
        cat(sprintf(paste0("  Convergence: passed; %d parameter(s) with R-hat ",
                           "between 1.01 and 1.05 and enough effective draws ",
                           "(see $info$convergence_diagnostics$rhat_caution)\n"),
                    n_caution))
    }
  }
  invisible(x)
}


#' What print() says about a bayesian_fit's Gaussian process
#'
#' The smallest effective range the basis resolves, in map units; the
#' posterior effective range when the fit identified one inside that; and one
#' flagged line for every way the basis was found short: before fitting
#' (prior mass under the floor, the response's own range under it) and after
#' (the adequacy status).  The range estimated before fitting is a variogram
#' estimate, and what answers it is the residuals of the fitted model.  It is
#' printed unflagged only when the status is \code{"ok"} and that status
#' rests on a residual test that found no autocorrelation at short range;
#' the line then gives the test.  A status of \code{"ok"} on the posterior
#' alone does not unflag it: a Gaussian process that carries a share of the
#' variance at a long range says nothing about an 8 km range beside it, and
#' "the fitted model does not bear that out" was printed over such fits.  A
#' fit saved by an earlier version carries none of these fields and prints
#' none of these lines.
#'
#' @param info A \code{bayesian_fit}'s \code{$info}, or a summary's.
#' @return Character vector of lines, possibly empty.
#' @keywords internal
#' @noRd
.bayes_gp_print_lines <- function(info) {
  out  <- character(0)
  fmt  <- function(v) format(signif(v, 3), big.mark = ",", scientific = FALSE,
                             trim = TRUE)
  res  <- suppressWarnings(as.numeric(info$gp_resolution))
  ad   <- info$gp_adequacy
  unit <- ad$unit %||% info$gp_prefit$unit %||% "CRS units"
  st   <- ad$status
  if (length(res) == 2L && all(is.finite(res)))
    out <- c(out, sprintf(
      "  GP resolves: effective ranges from %s (x) and %s (y) %s; shorter structure is counted as noise",
      fmt(res[1L]), fmt(res[2L]), unit))
  pm <- suppressWarnings(as.numeric(ad$prior_mass_below %||% NA_real_))[1L]
  if (is.finite(pm) && pm > 0.10)
    out <- c(out, sprintf(
      "  ** %.0f%% of the length-scale prior lies below what this basis resolves (gp_k too small for the prior)",
      100 * pm))
  pf <- info$gp_prefit
  # The residual test the status was split on, when there is one.
  rm_test <- ad$residual_moran
  rm_ok   <- is.list(rm_test) && is.numeric(rm_test$I) && is.numeric(rm_test$p) &&
    length(rm_test$I) == 1L && length(rm_test$p) == 1L &&
    is.finite(rm_test$I) && is.finite(rm_test$p)
  rm_txt  <- if (rm_ok)
    sprintf(" (Moran's I = %.3f, one-sided p = %.2g)", rm_test$I, rm_test$p) else ""
  if (identical(pf$status, "below_floor")) {
    out <- c(out, if (identical(st, "ok") && rm_ok) sprintf(
      "  Range before fitting: a variogram put it at %s %s, below what the basis resolves; after fitting, the residuals show no autocorrelation at short range%s",
      fmt(pf$range), unit, rm_txt)
    else if (identical(st, "ok")) sprintf(
      "  ** Range before fitting: a variogram put it at %s %s, below what the basis resolves (an estimate; the residuals were not tested after fitting)",
      fmt(pf$range), unit)
    else sprintf(
      "  ** Range before fitting: a variogram put it at %s %s, below what the basis resolves (an estimate; see GP adequacy)",
      fmt(pf$range), unit))
  } else if (identical(pf$status, "below_floor_long_range")) {
    # A short range beside structure that reaches the floor.  The fitted
    # model says whether the short range was there as well: in residuals
    # that are still autocorrelated, or in length-scale draws below the floor.
    out <- c(out, if (isTRUE(st %in% c("unresolved_short_range", "below_resolution"))) sprintf(
      "  ** Range before fitting: a variogram put it at %s %s, below what the basis resolves; the response is also autocorrelated at the floor's distance (see GP adequacy)",
      fmt(pf$range), unit)
    else sprintf(
      "  Range before fitting: a variogram put it at %s %s, but the response is still autocorrelated at the floor's distance, so that is not read as its range",
      fmt(pf$range), unit))
  } else if (identical(pf$status, "not_estimated"))
    out <- c(out, sprintf("  Range before fitting: not estimated (%s)",
                          pf$reason %||% "no reason recorded"))
  eff <- ad$effective_range
  # The share of the variance the GP carries, when that is what the status
  # was decided on; what the length-scales said, when it was not.
  sh  <- ad$spatial_share
  sh_txt <- if (is.numeric(sh) && length(sh) == 3L && all(is.finite(sh)))
    sprintf("the GP's share of the variance is %.3f to %.3f (95%% interval)",
            sh[[1L]], sh[[3L]])
    else if (identical(ad$identified, FALSE))
      "the data did not move the GP's length-scale from its prior"
    else "the GP fitted a surface"
  if (identical(st, "ok") && is.matrix(eff) && all(dim(eff) == c(2L, 3L)) &&
      all(is.finite(eff)))
    out <- c(out, sprintf(
      "  GP range: %s (x, 95%% interval %s to %s) and %s (y, %s to %s) %s, posterior median effective range",
      fmt(eff[1L, 2L]), fmt(eff[1L, 1L]), fmt(eff[1L, 3L]),
      fmt(eff[2L, 2L]), fmt(eff[2L, 1L]), fmt(eff[2L, 3L]), unit))
  else if (identical(st, "below_resolution"))
    out <- c(out, sprintf(paste0(
      "  ** GP adequacy: below_resolution -- %.0f%% of length-scale draws fall ",
      "below what the basis resolves; the fitted range is too long"),
      100 * (ad$frac_below %||% NA_real_)))
  else if (identical(st, "unresolved_short_range"))
    out <- c(out, sprintf(paste0(
      "  ** GP adequacy: unresolved_short_range -- %s, and the residuals are ",
      "autocorrelated at short range%s"), sh_txt, rm_txt))
  else if (identical(st, "no_signal_above_floor"))
    out <- c(out, sprintf(paste0(
      "  GP adequacy: no_signal_above_floor -- %s, and the residuals %s"),
      sh_txt,
      if (is.null(ad$residual_moran)) "were not tested"
      else "show no autocorrelation at short range"))
  if (!is.null(st) && !(st %in% c("ok", "not_checked")))
    out <- c(out, "     (see $info$gp_adequacy)")
  out
}


#' Coerce and check a response column before scoring it
#'
#' \code{is.finite()} on a character column is \code{FALSE} everywhere, so an
#' unchecked non-numeric response silently produced \code{n = 0} and all-\code{NA}
#' metrics from perfectly good data; a factor died inside \code{abs()} instead.
#' Logical is a legitimate 0/1 response and is coerced.
#'
#' @param y The response values.
#' @param nm The response column's name, for the message.
#' @param caller Function name, for the message.
#' @return \code{y}, numeric.
#' @keywords internal
#' @noRd
.checked_response <- function(y, nm, caller) {
  if (is.logical(y)) return(as.numeric(y))
  if (!is.numeric(y))
    stop(sprintf(paste0("%s(): response '%s' is %s, not numeric, so no metric ",
                        "can be computed from it. Convert the column first."),
                 caller, nm, class(y)[1L]), call. = FALSE)
  y
}


#' Summarise a fitted spatial model
#'
#' Computes goodness-of-fit metrics from \code{fitted(object)} against the
#' observed response.  A non-numeric response is an error: a character or
#' factor response cannot be scored, and used to come back as \code{n = 0} with
#' every metric \code{NA}; a logical response is treated as 0/1.
#'
#' @section What the metrics are computed on:
#' For a \code{gwr_fit} or a \code{bayesian_fit} these are \strong{in-sample}
#' metrics: \code{fitted()} returns values computed at the training locations
#' from the model that saw them.  For an \code{rf_fit} they are
#' \strong{out-of-bag}, because \code{fitted.rf_fit()} returns out-of-bag
#' predictions in place of in-sample ones (see \code{\link{fit_rf_model}}).
#' The two are not comparable, and \code{print()} on the result labels which
#' one it is holding, driven by \code{$info$fitted_are_oob}.  Use
#' \code{\link{compare_models_cv}} to compare backends.
#'
#' Adjusted R-squared is suppressed: GWR's effective parameter count far
#' exceeds the global predictor count, and a GP model has no simple \code{p}.
#'
#' For a \code{gwr_fit}, \code{print()} on the result adds what the
#' in-sample metrics are to be read against: the effective number of
#' parameters and residual degrees of freedom, the number of local windows
#' of deficient rank, the number of locations when observations share one,
#' and the AICc and RMSE of one global least-squares
#' regression on the same rows beside the fit's own AICc
#' (\code{$info$global_ols}; see \code{\link{fit_gwr_model}}).
#'
#' @param object A \code{spatial_fit} object.
#' @param ... Ignored.
#' @inheritSection model_metrics Percentage errors on responses with zeros
#' @return An object of class \code{summary.spatial_fit}: a list with
#'   \code{class}, \code{formula}, \code{n}, \code{response_var},
#'   \code{predictor_vars}, \code{info} and \code{in_sample} (the metric
#'   data.frame, out-of-bag for an \code{rf_fit}).
#' @family methods on a fitted model
#' @export
summary.spatial_fit <- function(object, ...) {
  fit_vals <- .fitted_checked(object, .caller = "summary")
  y_obs    <- sf::st_drop_geometry(object$data_sf)[[object$response_var]]
  # The same guard model_metrics() applies, for the same reason and on the same
  # columns.  Without it summary() still did what model_metrics() was fixed to
  # stop doing: a character response gave n = 0 and every metric NA with
  # nothing said, and a factor died in `abs()`.  new_spatial_fit() is the
  # documented entry point for a user's own backend and validates neither.
  y_obs <- .checked_response(y_obs, object$response_var, "summary")
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
    # `.cache` is an ENVIRONMENT, so carrying it here would not copy it: the
    # summary would hold the fit's live cache.  That made a summary something
    # other than a value snapshot (its contents changed when anyone later
    # called fitted() on the fit), made clear_fitted_cache(summary) empty the
    # FIT's cache, and made saveRDS()/identical() on two summaries of
    # structurally identical fits disagree.  Drop it.
    info           = object$info[setdiff(names(object$info), ".cache")],
    in_sample      = met
  )
  class(out) <- "summary.spatial_fit"
  out
}


#' @export
print.summary.spatial_fit <- function(x, ...) {
  cat(sprintf("Summary of <%s> fit (n = %d)\n\n", x$class, x$n))
  cat(sprintf("  Formula: %s\n", paste(deparse(x$formula), collapse = " ")))
  # An rf_fit's fitted() values are out-of-bag, not in-sample.  Labelling them
  # "In-sample" understates a random hold-out as a memorised fit, which is the
  # opposite of the truth and the opposite of the caveat ?fit_rf_model carries.
  if (isTRUE(x$info$fitted_are_oob))
    cat("\n  Out-of-bag metrics (NOT in-sample; see ?fit_rf_model):\n")
  else
    cat("\n  In-sample metrics:\n")
  m <- x$in_sample
  # The metrics use the rows with a finite fitted value, which is not always
  # all of them: an rf_fit has no out-of-bag prediction for a row every tree
  # sampled (ranger returns NaN), so a 5-tree forest printed "n = 200" above
  # an R^2 computed on 180 rows.
  n_m <- m$n %||% NA_integer_
  if (length(n_m) == 1L && is.finite(n_m) && is.finite(x$n) && n_m < x$n)
    cat(sprintf("    (computed on %d of %d rows; the rest have no finite fitted value)\n",
                n_m, x$n))
  # ASCII on purpose: a superscript two rendered as R<U+00B2> on every
  # non-UTF-8 console, and the labels were not aligned.
  cat(sprintf("    RMSE    = %.4f\n", m$RMSE))
  cat(sprintf("    MAE     = %.4f\n", m$MAE))
  cat(sprintf("    R^2     = %.4f\n", m$R2))
  if (is.finite(m$Adj_R2 %||% NA_real_))
    cat(sprintf("    Adj R^2 = %.4f\n", m$Adj_R2))
  if (is.finite(m$SMAPE %||% NA_real_)) {
    # SMAPE skips rows where observation and prediction are both zero, so a
    # value averaged over a subset says so (see ?model_metrics).
    n_s <- m$n_SMAPE %||% NA_integer_
    sub <- if (is.finite(n_s) && is.finite(m$n) && n_s < m$n)
      sprintf("  (over %d of %d rows)", n_s, m$n) else ""
    cat(sprintf("    SMAPE   = %.2f%%%s\n", m$SMAPE, sub))
  }
  # What the local regressions cost and how many could not estimate every
  # coefficient: the two numbers the in-sample metrics of a GWR are read
  # against.  From the fit's own decomposition (see ?fit_gwr_model).
  if (identical(x$class, "gwr_fit")) {
    enp <- suppressWarnings(as.numeric(x$info$enp %||% NA_real_))[1L]
    edf <- suppressWarnings(as.numeric(x$info$edf %||% NA_real_))[1L]
    if (is.finite(enp) && is.finite(edf))
      cat(sprintf(paste0("\n  Effective parameters: %.1f (residual degrees of ",
                         "freedom %.1f)\n"), enp, edf))
    n_def <- suppressWarnings(as.integer(x$info$n_rank_deficient %||% 0L))[1L]
    if (isTRUE(n_def > 0L))
      cat(sprintf(paste0("  Rank-deficient windows: %d of %d (their ",
                         "unestimable coefficients are NA)\n"), n_def, x$n))
    # Observations that share a location: the AICc lines below and the
    # standard errors count them as independent.
    n_loc <- suppressWarnings(as.integer(x$info$n_locations %||% NA_integer_))[1L]
    if (isTRUE(n_loc < x$n))
      cat(sprintf(paste0("  Locations: %d for %d observations (AICc and the ",
                         "standard errors count\n    the observations at one ",
                         "location as independent)\n"), n_loc, x$n))
    # One global regression on the same rows, by the same AICc: what the GWR
    # has to beat to have shown that its coefficients vary.
    ols  <- x$info$global_ols
    aicc <- suppressWarnings(as.numeric(x$info$AICc %||% NA_real_))[1L]
    if (is.list(ols) && is.finite(ols$AICc %||% NA_real_)) {
      cat(sprintf(paste0("  Global least squares on the same rows: AICc %.2f, ",
                         "RMSE %.4f\n"), ols$AICc, ols$RMSE))
      if (is.finite(aicc))
        cat(sprintf("  This fit: AICc %.2f (%s the global fit's by %.2f)\n",
                    aicc, if (aicc < ols$AICc) "below" else "not below",
                    abs(aicc - ols$AICc)))
    }
  }
  # The same convergence verdict print() on the fit gives.  The summary carried
  # it in $info and never showed it, so metrics from a posterior that had not
  # converged printed exactly like metrics from one that had.
  if (identical(x$class, "bayesian_fit") && "convergence_ok" %in% names(x$info)) {
    cv_ok <- x$info$convergence_ok
    if (length(cv_ok) == 1L && is.na(cv_ok))
      cat("\n  Convergence: NOT CHECKED (fitted with check_convergence = FALSE?)\n")
    else if (!isTRUE(cv_ok))
      cat("\n  ** Convergence warnings present -- see $info$convergence_diagnostics\n")
    else {
      # The caution print() on the fit gives: passed, with R-hat between 1.01
      # and 1.05 on parameters that have enough effective draws.
      n_caution <- length(x$info$convergence_diagnostics$rhat_caution)
      if (n_caution > 0L)
        cat(sprintf(paste0("\n  Convergence: passed; %d parameter(s) with R-hat ",
                           "between 1.01 and 1.05 and enough effective draws ",
                           "(see $info$convergence_diagnostics$rhat_caution)\n"),
                    n_caution))
    }
  }
  # What the basis could not hold changes how the metrics above read: a GP
  # that counted short-range structure as noise scores like one that found
  # none.  The same lines print() on the fit gives, for any status but "ok".
  if (identical(x$class, "bayesian_fit")) {
    st <- x$info$gp_adequacy$status
    if (!is.null(st) && !(st %in% c("ok", "not_checked"))) {
      cat("\n")
      for (ln in .bayes_gp_print_lines(x$info)) cat(ln, "\n", sep = "")
    }
  }
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
#' numbers are in-sample for a \code{gwr_fit} or \code{bayesian_fit}, and a
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
#' data.frame \code{model_metrics()} returns carries no label distinguishing
#' the two, so check \code{object$info$fitted_are_oob} before comparing
#' numbers across backends; \code{\link{evaluate_insample}()} and
#' \code{\link{compare_models}()} record it per model in a
#' \code{metric_basis} column.  \code{\link{compare_models_cv}} scores every
#' backend the same way.
#'
#' \eqn{R^2} is \eqn{1 - RSS/TSS} with the total sum of squares taken about
#' the mean of the response the model was \emph{fitted} to.  In sample that
#' is the ordinary \eqn{R^2}.  With \code{newdata} it is out-of-sample
#' \eqn{R^2}, the convention every \code{cv_*()} function uses: the model is
#' measured against the prediction it had to beat, the training mean, not
#' against the new rows' own mean, which it could not have known.  It is
#' below 0 when the model predicts the new rows worse than the training mean
#' does, and it is \code{NA} when the response does not vary about that
#' baseline by more than rounding error (100 machine epsilons of its
#' magnitude, whatever its units).
#'
#' @section Percentage errors on responses with zeros:
#' \code{MAPE} divides by the observed value and \code{SMAPE} by
#' \eqn{|y| + |\hat{y}|}, so neither is defined where its denominator is zero.
#' Neither returns \code{Inf} or \code{NaN}.  Both are averaged over the rows
#' whose denominator is non-zero, and are \code{NA} when no row qualifies.
#' Non-zero is judged at the scale of the data: a denominator no larger
#' than 100 machine epsilons times the largest one counts as zero, so the
#' rule does not depend on the units of the response.
#' The \code{n_MAPE} and \code{n_SMAPE} columns record how many rows that was;
#' the \code{n} column counts finite observation/prediction pairs.  Read a
#' percentage error next to its count: when \code{n_MAPE < n}, \code{MAPE} is
#' an average over a subset of the data, whatever its value.
#'
#' This bites on any response taking exact zeros: counts, rainfall,
#' abundance, claim amounts.  On a zero-inflated response with 62 zeros out of
#' 120, \code{MAPE} is an average over the 58 non-zero rows, which
#' \code{n_MAPE = 58} now says.  \code{SMAPE} fails differently and more
#' subtly: it drops the rows where observation and prediction are both near
#' zero (which on a well-fitted zero-inflated model are the rows it got
#' \emph{right}), so it averages the harder rows only and reads worse than
#' the fit deserves; \code{n_SMAPE} shows how many rows it kept, and the
#' count is only a label, not a repair.
#'
#' \code{RMSE}, \code{MAE} and \eqn{R^2} use every finite row and are
#' unaffected; prefer them whenever the response can be zero.  For a Bayesian
#' fit, \code{\link{cv_bayes}()} additionally reports CRPS and interval
#' coverage, which are proper scoring rules and have no such failure mode.
#'
#' @param object A \code{spatial_fit} object.
#' @param newdata Optional sf object for out-of-sample evaluation. If NULL,
#'   fitted values are used (in-sample, or out-of-bag for an \code{rf_fit};
#'   see above).
#' @param ... Additional arguments passed to predict().  The metrics are
#'   computed from one predicted value per row, so \code{draws = TRUE} (a
#'   Bayesian fit) and \code{se.fit = TRUE} (a GWR), with which
#'   \code{predict()} returns a matrix or a list, are refused.
#' @section Which metrics survive a non-Gaussian response:
#' RMSE and MAE are defined for any numeric response and are what to read for
#' a count, a rate or a bounded outcome.  MAPE and SMAPE assume a response
#' that is rarely zero (see the previous section), and R-squared and adjusted
#' R-squared compare residual variance to total variance, which is the right
#' comparison for a Gaussian response and a loose one for anything whose
#' variance tracks its mean.  None of the four is wrong to compute; each is
#' Gaussian-shaped thinking, and on a Poisson or zero-inflated response should
#' be read as a rough summary rather than a score.
#'
#' For the Bayesian backend, \code{\link{cv_bayes}()} additionally reports
#' CRPS and interval coverage at 50, 80 and 95 percent.  Both are proper
#' scoring rules computed from posterior draws, so they are meaningful for any
#' \code{family} that predicts one number per row (a count, a rate, a binary
#' or bounded outcome), and they are the numbers to compare when the response
#' is not Gaussian.  A categorical or ordinal family predicts a probability
#' per response category instead, so \code{cv_bayes()} refuses one before
#' fitting anything.  When every fold fails, the
#' \code{fold_metrics} frame \code{cv_bayes()} returns carries the CRPS column
#' but not the \code{coverage_*} columns, so code that reads those columns
#' must tolerate their absence.
#'
#' @return A data.frame with n, RMSE, MAE, MAPE, SMAPE, R2, Adj_R2, n_MAPE and
#'   n_SMAPE (the last two are the rows each percentage error was averaged
#'   over; see "Percentage errors on responses with zeros").
#'   \code{Adj_R2} is always \code{NA}: GWR's effective parameter count far
#'   exceeds the global predictor count and a GP model has no simple \code{p},
#'   so it is deliberately suppressed.  A non-numeric response is an error: a
#'   character or factor response cannot be scored, and used to come back as
#'   \code{n = 0} with every metric \code{NA}; a logical response is treated
#'   as 0/1.
#' @family model evaluation
#' @examples
#' if (requireNamespace("ranger", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   pts <- st_as_sf(
#'     data.frame(x = 5e5 + runif(60, 0, 1000), y = 5e6 + runif(60, 0, 1000),
#'                a = rnorm(60)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   pts$z <- 2 * pts$a + rnorm(60, 0, 0.3)
#'   # Fit on 40 rows and keep 20 back, so the second call really is out of
#'   # sample; scoring the training rows again would only re-read the fit.
#'   fit <- fit_rf_model(pts[1:40, ], "z", "a", num_trees = 50, seed = 1)
#'   print(model_metrics(fit))                    # in-sample (out-of-bag for RF)
#'   model_metrics(fit, newdata = pts[41:60, ])   # on the 20 held-out rows
#' }
#' @export
model_metrics <- function(object, ...) UseMethod("model_metrics")


#' @family methods on a fitted model
#' @rdname model_metrics
#' @export
model_metrics.spatial_fit <- function(object, newdata = NULL, ...) {
  .check_dots_newdata(list(...), newdata, "model_metrics")
  # The metrics are computed from one predicted value per row.  Two predict()
  # arguments make a method return something else: `draws = TRUE` (a
  # bayesian_fit: a draws x rows matrix) and `se.fit = TRUE` (a gwr_fit: a
  # list).  Either reached .compute_reg_metrics() and stopped with "missing
  # value where TRUE/FALSE needed", or with an error about a list.
  dots <- list(...)
  shaped <- intersect(names(dots)[vapply(dots, isTRUE, logical(1))],
                      c("draws", "se.fit"))
  if (length(shaped))
    stop("model_metrics(): ", paste(sprintf("`%s = TRUE`", shaped), collapse = " and "),
         " cannot be passed through `...`: the metrics are computed from one ",
         "predicted value per row, and with it predict() returns ",
         if (identical(shaped, "draws")) "a matrix of posterior draws"
         else if (identical(shaped, "se.fit")) "a list with the standard errors"
         else "more than that",
         ". Call predict() yourself for those, or cv_bayes() for CRPS and ",
         "interval coverage.", call. = FALSE)
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
  y_obs <- .checked_response(y_obs, object$response_var, "model_metrics")
  # R² on newdata is measured against the TRAINING mean, as every cv_*()
  # measures it: the null prediction the model had to beat.  It used the
  # held-out rows' own mean, so the same predictions scored R² -0.89 here and
  # 0.35 from cv_spatial() on a trend split.  In sample the two means are the
  # same number.  A fit whose data_sf lacks a usable response keeps the
  # held-out mean.
  ytm <- NULL
  if (!is.null(newdata)) {
    y_tr <- tryCatch(suppressWarnings(as.numeric(
      sf::st_drop_geometry(object$data_sf)[[object$response_var]])),
      error = function(e) NULL)
    if (length(y_tr) && any(is.finite(y_tr))) ytm <- mean(y_tr[is.finite(y_tr)])
  }
  # Adj R² is suppressed (p = NULL) because GWR's effective parameter count
  # far exceeds the global predictor count, and Bayesian GP models likewise
  # lack a simple p.  This is consistent with the CV evaluation path.
  .compute_reg_metrics(y_obs, y_hat, p = NULL, y_train_mean = ytm)
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
#' Otherwise solves one local regression at each new location, with the
#' fit's kernel and bandwidth (an adaptive bandwidth counts neighbours among
#' the training points), and returns \eqn{x^\top\hat\beta(u)}{x'beta(u)}.
#' The local regressions are weighted least squares solved by a pivoted QR,
#' the solve \code{\link{fit_gwr_model}()} uses; on a window of full rank the
#' values equal those of \code{GWmodel::gwr.predict()} to rounding error, and
#' at a training location with that row's own predictors they equal
#' \code{fitted()}.  Each location stands alone: one that cannot be
#' estimated does not affect the others.
#' \code{newdata} is first transformed to the CRS used during fitting
#' (via \code{ensure_projected()}), so predictions are computed in a
#' single coordinate system regardless of the CRS newdata arrives in.
#'
#' @section When a prediction is NA:
#' A prediction is returned when it is estimable from the training points
#' that carry weight at the new location, and is \code{NA} otherwise:
#' \itemize{
#'   \item no training point carries weight there (a location beyond a fixed
#'     bisquare, tricube or boxcar bandwidth);
#'   \item the training points in the window cannot place the new row: a
#'     predictor is constant among them at a value other than the new row's
#'     (a regional indicator, at a location whose window lies wholly in
#'     another region), or fewer of them carry weight than the new row needs.
#' }
#' Where a predictor is constant in the window at the new row's own value,
#' its coefficient is not needed and the prediction is returned.  The
#' decision uses the rule and the tolerance of the fit (see
#' \code{\link{fit_gwr_model}}, \strong{The local solve}).  A warning of
#' class \code{"spatialkit_gwr_not_estimable"} counts the \code{NA}
#' predictions.
#'
#' @section Locations beyond the reach of the fit:
#' An adaptive bandwidth always finds its neighbours, and a Gaussian or an
#' exponential kernel weights every training point, so a location far from
#' the data still gets a number: the local regression of training points
#' that are all far away, extrapolated.  When a new location's nearest
#' training point is farther than the widest kernel used at any training
#' location reaches, the prediction is returned and a warning of class
#' \code{"spatialkit_gwr_extrapolation"} counts those locations and gives
#' the range of their distances.  The widest kernel is the largest distance
#' from a training point to its \code{bandwidth}-th nearest neighbour
#' (\code{object$info$max_kernel_width}), or the bandwidth itself when it is
#' fixed.  A bisquare, tricube or boxcar kernel reaches as far as its width;
#' a Gaussian or exponential kernel has no edge, and its reach is taken as
#' three widths.  The check is on distance alone and is lenient: a location
#' just outside the data is not flagged.
#'
#' The warning is raised once per call, with that call's counts.
#' \code{\link{predict_surface}()}, which predicts a grid in chunks, and
#' \code{\link{cv_gwr}()} and \code{\link{cv_spatial}()}, which predict
#' fold by fold, collect the warnings of their pieces and raise one under
#' their own name, with the totals and the number of calls they are over
#' (\code{calls}).  The folds of a cross-validation are separate fits, so
#' a cross-validation gives the smallest and the largest reach among them.
#' \code{\link{cv_block_size_sweep}()} and
#' \code{\link{select_features_forward}()} run many cross-validations of
#' one learner and raise one warning for the sweep.  A sweep predicts the
#' same rows from several fits, so its totals count a location once for each
#' fit that put it out of reach, and the message says so.
#' \code{model_metrics(fit, newdata)} is
#' one call and raises the warning once.  A loop of your own over
#' \code{predict()} gets one warning per call; the condition carries its
#' numbers as fields (\code{n_far},
#' \code{n}, \code{reach}, and \code{nearest}, the smallest and largest
#' distance to a nearest training point), so a handler can total them:
#' \code{withCallingHandlers(..., spatialkit_gwr_extrapolation = function(w)
#' \{ total <<- total + w$n_far; invokeRestart("muffleWarning") \})}.  To
#' drop it where extrapolation is intended, use
#' \code{suppressWarnings(predict(fit, newdata), classes =
#' "spatialkit_gwr_extrapolation")}.  The warning that counts \code{NA}
#' predictions (class \code{"spatialkit_gwr_not_estimable"}) carries
#' \code{n_na}, \code{n} and \code{n_empty} and is collected the same way.
#'
#' @param object A \code{gwr_fit} object.
#' @param newdata An sf object with the same predictors.
#'   The response variable need not be present (true out-of-sample prediction
#'   is supported).  NULL = fitted values.
#' @param se.fit Logical; also return standard errors.  Default \code{FALSE}.
#' @param ... Ignored.
#' @return With \code{se.fit = FALSE}, a numeric vector aligned to
#'   \code{nrow(newdata)}, with \code{NA} for rows dropped as missing or
#'   non-finite, and for locations whose local regression cannot be
#'   estimated (see \strong{When a prediction is NA}); a warning counts
#'   those.  If the design matrix for \code{newdata} cannot be built, every
#'   value is \code{NA} and a warning says why.  CRS-less
#'   \code{newdata} first receives the interpretation the training data got, so
#'   the same rows land where they did at fit time.
#'
#'   With \code{se.fit = TRUE}, a list: \code{fit} (that vector),
#'   \code{se.fit} (the standard error of the fitted value,
#'   \eqn{\hat\sigma\sqrt{s^\top s}} with \eqn{s} the weights the local
#'   regression puts on the training responses), \code{se.pred} (the standard
#'   error of a new observation there, \eqn{\hat\sigma\sqrt{1 + s^\top s}},
#'   the square root of \code{GWmodel::gwr.predict()}'s prediction variance),
#'   \code{residual.scale} (\eqn{\hat\sigma}, from the fit's residual sum of
#'   squares over its effective degrees of freedom) and \code{df} (those
#'   degrees of freedom).  The standard errors treat the bandwidth as known
#'   and the local model as true at each location; they do not measure the
#'   bias of a local fit and are \code{NA} where the prediction is.
#' @family methods on a fitted model
#' @export
predict.gwr_fit <- function(object, newdata = NULL, se.fit = FALSE, ...) {
  .check_dots(list(...), "predict.gwr_fit")
  if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit))
    stop("predict.gwr_fit(): `se.fit` must be TRUE or FALSE.", call. = FALSE)
  if (is.null(newdata)) {
    if (!se.fit) return(fitted(object))
    # The standard errors need each row's weights on the training responses,
    # which the fit does not keep; the training rows are predicted afresh.
    newdata <- object$data_sf
  }
  .check_predict_newdata(newdata, object, "predict.gwr_fit")

  if (!requireNamespace("sp", quietly = TRUE))
    stop("predict.gwr_fit(): package 'sp' is required.", call. = FALSE)

  # Record original length so we can return an aligned vector with NA fill.
  n_orig <- nrow(newdata)

  # Align newdata to the CRS used during fitting BEFORE prep_model_data().
  # prep's own ensure_projected() call has no target and would leave
  # already-projected newdata in *its* CRS (or auto-pick a UTM zone for
  # lon/lat input independent of training), after which the local regressions
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

  # What comes back when nothing can be predicted, in the shape asked for.
  na_result <- function() {
    v <- rep(NA_real_, n_orig)
    if (!se.fit) return(v)
    list(fit = v, se.fit = v, se.pred = v,
         residual.scale = NA_real_, df = NA_real_)
  }

  # Degrade to an all-NA vector when nothing survived cleaning, matching
  # predict.rf_fit() and predict.bayesian_fit().  Without this a zero-row
  # layer would reach the local regressions below and come back as a warning
  # about them, which are not what went wrong.
  if (n_new == 0L) {
    .log_warn(paste0("predict.gwr_fit(): every row of newdata was dropped as ",
                     "missing or non-finite; returning %d NA(s)."), n_orig)
    return(na_result())
  }

  needed_cols <- unique(c(object$response_var, object$predictor_vars))
  sp_train <- .to_sp(object$data_sf, needed_cols)

  bw <- object$info$bandwidth
  adaptive <- object$info$adaptive %||% TRUE
  if (isTRUE(adaptive)) bw <- as.integer(round(bw))
  kernel <- object$info$kernel %||% "bisquare"

  # One local regression per new location (see .gwr_predict_at()), each
  # solved for itself, so a location whose window is empty or cannot place
  # its row (a grid cell beyond a fixed bandwidth, a point inside a held-out
  # block) is NA and the others are what they would be without it.
  res <- tryCatch(
    .gwr_predict_at(object$formula, sp_train,
                    newdata_df = sf::st_drop_geometry(newdata),
                    rp = sf::st_coordinates(newdata)[, 1:2, drop = FALSE],
                    bw = bw, kernel = kernel, adaptive = adaptive,
                    se = se.fit),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    # A real warning(), not only a logger line.  A logger line is invisible
    # to tryCatch(warning = ), to withCallingHandlers(), to
    # testthat::expect_warning() and to R CMD check, so a predict() that
    # returns nothing but NA left no trace a caller could act on.
    .warn_and_log(paste0("predict.gwr_fit(): the local regressions could not ",
                         "be evaluated, so every prediction is NA. Cause: %s"),
                  conditionMessage(res))
    return(na_result())
  }
  preds_clean <- res$fitted

  # The two counted warnings carry their counts as fields of the condition
  # (see .gwr_tally_predict_warnings(), which a caller that predicts in
  # chunks or folds uses to raise each once with totals).
  if (anyNA(preds_clean)) {
    is_na   <- is.na(preds_clean)
    n_na    <- sum(is_na)
    n_empty <- sum(is_na & res$n_window == 0L)
    bw_text <- sprintf("%s bandwidth (%s)",
                       if (isTRUE(adaptive)) "adaptive" else "fixed", format(bw))
    .gwr_warn(
      "spatialkit_gwr_not_estimable", "%s",
      .gwr_not_estimable_msg("predict.gwr_fit", n_na, n_new, n_empty, bw_text),
      data = list(n_na = n_na, n = n_new, n_empty = n_empty,
                  bw_text = bw_text, calls = 1L),
      tallied = TRUE)
  }

  # A prediction at a location no window of the fit reached.  An adaptive
  # bandwidth always finds its neighbours and a Gaussian kernel weights every
  # point, so far from the data a number still comes back; it is counted
  # here, against the reach of the widest kernel at any training location.
  kr  <- .gwr_kernel_reach(object, sp::coordinates(sp_train), bw, adaptive,
                           kernel)
  far <- !is.na(preds_clean) & is.finite(res$nearest) &
    is.finite(kr$reach) & res$nearest > kr$reach
  if (any(far)) {
    rng <- range(res$nearest[far])
    reach_text <- sprintf(
      "%s in the units of the fit's CRS%s", format(signif(kr$reach, 4)),
      if (kernel %in% c("gaussian", "exponential"))
        sprintf(", three times the widest %s bandwidth of %s", kernel,
                format(signif(kr$width, 4))) else "")
    .gwr_warn(
      "spatialkit_gwr_extrapolation", "%s",
      .gwr_extrapolation_msg("predict.gwr_fit", sum(far), n_new, reach_text,
                             rng),
      data = list(n_far = sum(far), n = n_new, reach = kr$reach,
                  reach_text = reach_text, nearest = rng, calls = 1L),
      tallied = TRUE)
  }

  # Expand back to original length, filling dropped rows with NA.
  fit_vals <- .expand_predictions(preds_clean, clean_idx, n_orig)
  if (!se.fit) return(fit_vals)

  sigma2 <- suppressWarnings(as.numeric(object$info$sigma2 %||% NA_real_))[1L]
  if (!is.finite(sigma2))
    .warn_and_log(paste0("predict.gwr_fit(): the fit carries no residual ",
                         "variance (it has no effective residual degrees of ",
                         "freedom, a fitted value of its own is NA, or it was ",
                         "not made by fit_gwr_model()), so `se.fit` and ",
                         "`se.pred` are NA."))
  vf <- res$var_fit
  vf[is.na(preds_clean)] <- NA_real_
  list(fit            = fit_vals,
       se.fit         = .expand_predictions(sqrt(sigma2 * vf), clean_idx, n_orig),
       se.pred        = .expand_predictions(sqrt(sigma2 * (1 + vf)), clean_idx,
                                            n_orig),
       residual.scale = sqrt(sigma2),
       df             = suppressWarnings(
         as.numeric(object$info$edf %||% NA_real_))[1L])
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

  # Hand brms each categorical predictor with exactly the levels it was
  # fitted with.  brms compares the LEVEL SET of a newdata factor with the
  # levels it kept, not the values the rows carry, and it drops unused levels
  # when fitting.  So a factor that merely listed an extra level no row had,
  # and the training data itself whenever one of its levels had no rows (any
  # ordinary subset), were refused with "New factor levels are not allowed":
  # fitted() and summary() stopped, and predict() returned NA for every row.
  # A value outside the fitted levels becomes NA here, and brms would predict
  # the reference level for it; predict.bayesian_fit() removes those rows
  # before this is reached.
  lv <- .bayes_fit_levels(object)
  for (cn in intersect(names(lv), names(pred_df)))
    pred_df[[cn]] <- factor(as.character(pred_df[[cn]]),
                            levels = lv[[cn]]$levels,
                            ordered = isTRUE(lv[[cn]]$ordered))

  pred_df
}


#' The levels a bayesian_fit's categorical predictors were fitted with
#'
#' Stored at fit time as \code{$info$predictor_levels}.  A fit saved by an
#' earlier version has none, and gets them the way they were made: from the
#' values that occur in its training rows.
#'
#' @param object A \code{bayesian_fit}.
#' @return A named list, one element per factor or character predictor, each
#'   a list with \code{levels} and \code{ordered}.
#' @keywords internal
#' @noRd
.bayes_fit_levels <- function(object) {
  lv <- object$info$predictor_levels
  if (is.list(lv)) return(lv)
  df <- tryCatch(sf::st_drop_geometry(object$data_sf), error = function(e) NULL)
  if (is.null(df)) return(list())
  .bayes_predictor_levels(df, object$predictor_vars)
}


#' Rows of newdata that carry a level the model was not fitted with
#'
#' @param object A \code{bayesian_fit}.
#' @param newdata The prediction layer.
#' @return A list with \code{unseen} (logical, one per row), \code{bad} (a
#'   named list of per-column logical masks, only for columns with such
#'   rows) and \code{named} (one string per such column, naming the levels).
#'   A missing value is not an unseen level; the row cleaning handles it.
#' @keywords internal
#' @noRd
.bayes_unseen_levels <- function(object, newdata) {
  lv  <- .bayes_fit_levels(object)
  out <- list(unseen = rep(FALSE, nrow(newdata)), bad = list(),
              named = character(0))
  for (cn in intersect(names(lv), names(newdata))) {
    val <- as.character(newdata[[cn]])
    bad <- !is.na(val) & !(val %in% lv[[cn]]$levels)
    if (!any(bad)) next
    out$unseen    <- out$unseen | bad
    out$bad[[cn]] <- bad
    out$named     <- c(out$named, sprintf(
      "'%s': %s (fitted with %s)", cn,
      paste(sQuote(unique(val[bad]), q = FALSE), collapse = ", "),
      paste(sQuote(lv[[cn]]$levels, q = FALSE), collapse = ", ")))
  }
  out
}


#' Hold the GP boundary at the value the model was fitted with
#'
#' brms 2.17 to 2.22 store \code{Xgp}, \code{dmax} and \code{cmeans} in a fit's
#' GP basis but \strong{not} the Hilbert-space boundary \code{L}, so
#' \code{brms:::.data_gp()} recomputes
#' \code{L = c * max(1, diff(range(centred Xgp)))} from whatever rows
#' \code{predict()} is handed.  Every eigenfunction and eigenvalue of the
#' approximation therefore moves with the newdata bounding box while the fitted
#' basis coefficients stay put.  brms 2.23.0 stores \code{L} in the basis and
#' reuses it, so there the two repairs below are unnecessary and harmless.
#'
#' Measured on a fitted model: \code{L} was 5.57 at fit time, 4.02 for a
#' five-row \code{newdata} and 3.63 for one row; \code{predict_surface()} on the
#' same 10,000-cell grid differed by 1.78 between \code{chunk_size = 5000} (the
#' documented default) and one single call; and \code{cv_bayes()}, which
#' predicts each test fold separately, scored every fold against a basis the
#' model was never fitted with.
#'
#' The first repair makes the pooled centred range of the rows brms sees at
#' least the training one: append two synthetic rows at the training extrema,
#' predict, then drop their columns.  \code{cmeans} already comes from the
#' stored basis, so for newdata inside the training envelope that reproduces
#' the fitted \code{L} exactly.  Appending rows cannot narrow the range,
#' though: one row past the pooled extremes still widened \code{L}, and with it
#' the prediction of EVERY row in the call (on brms 2.20.4 one row 100 m past
#' the bbox moved interior predictions by up to 0.08, and a grid padded 15\%
#' past it moved in-bbox cells by 0.31 on average).  The second repair is for
#' that case: \code{c_scale = S_fit / S_new} (\code{.gp_c_scale()}), by which
#' \code{predict()} multiplies the \code{c} of the gp() term
#' (\code{.scale_gp_c()}), so brms rebuilds
#' \code{L = c * c_scale * S_new = c * S_fit}, the fitted value.
#'
#' Rows further than \code{L} from the training centre on either axis are
#' flagged in \code{beyond}.  The Dirichlet eigenfunctions vanish at
#' \code{+/- L} and continue past it as an odd reflection of the fitted
#' surface, so a prediction there means nothing on any brms version.
#'
#' @param object A \code{bayesian_fit}.
#' @param pred_df The prediction data.frame from
#'   \code{.prepare_brms_pred_df()}.
#' @return A list with \code{df} (possibly with rows appended), \code{n_pad}
#'   (how many were appended, to be dropped from the draw matrix),
#'   \code{c_scale} (the factor for the gp() term's \code{c}; 1 when the rows
#'   do not widen the range or the fit lacks what computing it needs) and
#'   \code{beyond} (logical, one per row of \code{pred_df}).
#' @keywords internal
#' @noRd
.pin_gp_boundary_rows <- function(object, pred_df) {
  none <- list(df = pred_df, n_pad = 0L, c_scale = 1,
               beyond = rep(FALSE, nrow(pred_df)))
  rng <- object$info$gp_xy_range
  if (is.null(rng) || !all(c("..x", "..y") %in% names(pred_df)) ||
      !nrow(pred_df))
    return(none)
  xr <- suppressWarnings(as.numeric(rng$x))
  yr <- suppressWarnings(as.numeric(rng$y))
  if (length(xr) != 2L || length(yr) != 2L || !all(is.finite(c(xr, yr))))
    return(none)

  # Say when newdata reaches outside the training envelope: with the boundary
  # held, those predictions no longer depend on the other rows, but they are
  # still extrapolation from a basis fitted without data there.
  ox <- range(pred_df[["..x"]], na.rm = TRUE)
  oy <- range(pred_df[["..y"]], na.rm = TRUE)
  if (any(is.finite(c(ox, oy))) &&
      (ox[1L] < xr[1L] || ox[2L] > xr[2L] ||
       oy[1L] < yr[1L] || oy[2L] > yr[2L]))
    .log_info(paste0("predict.bayesian_fit(): `newdata` reaches outside the ",
                     "training coordinate envelope; those predictions are ",
                     "extrapolation."))

  pad <- pred_df[c(1L, 1L), , drop = FALSE]
  pad[["..x"]] <- xr
  pad[["..y"]] <- yr
  rownames(pad) <- NULL
  out <- rbind(pred_df, pad)
  rownames(out) <- NULL
  res <- list(df = out, n_pad = 2L, c_scale = 1, beyond = none$beyond)

  # The centre brms uses.  A fit saved before gp_cmeans was stored gets it the
  # way brms got it: from the unique training rows it was fitted on.
  cm <- object$info$gp_cmeans
  ed <- object$engine$data
  if (is.null(cm) && is.data.frame(ed) && all(c("..x", "..y") %in% names(ed)))
    cm <- colMeans(unique(cbind(ed[["..x"]], ed[["..y"]])))
  cm    <- suppressWarnings(as.numeric(cm))
  S_fit <- suppressWarnings(as.numeric(object$info$gp_S))
  c_fit <- suppressWarnings(as.numeric(object$info$gp_c))
  if (length(cm) != 2L || !all(is.finite(cm)) ||
      length(S_fit) != 1L || !isTRUE(S_fit > 0) ||
      length(c_fit) != 1L || !isTRUE(c_fit > 0))
    return(res)

  L_fit <- c_fit * S_fit
  res$beyond  <- abs(pred_df[["..x"]] - cm[1L]) > L_fit |
                 abs(pred_df[["..y"]] - cm[2L]) > L_fit
  res$c_scale <- .gp_c_scale(cbind(out[["..x"]], out[["..y"]]), cm, S_fit)
  res
}


#' The factor that makes brms rebuild the fitted GP boundary from new rows
#'
#' brms 2.17 to 2.22 build \code{L = c * S_new} at predict time, where
#' \code{S_new = max(1, pooled range)} of the unique newdata rows centred on
#' the training \code{cmeans} (\code{brms:::.data_gp()},
#' \code{brms:::choose_L()}).  Handing brms \code{c * S_fit / S_new} in place
#' of \code{c} therefore gives back \code{L = c * S_fit}.  Pure arithmetic, so
#' it is tested without Stan.
#'
#' @param xy Numeric matrix; its first two columns are the scaled coordinates
#'   of every row brms will be handed.
#' @param cmeans Length-2 training column means.
#' @param S_fit The pooled centred training range (\code{$info$gp_S}).
#' @return \code{S_fit / S_new}; 1 when no coordinate is finite.
#' @keywords internal
#' @noRd
.gp_c_scale <- function(xy, cmeans, S_fit) {
  xy <- unique(as.matrix(xy)[, 1:2, drop = FALSE])
  Xc <- sweep(xy, 2L, cmeans)
  if (!any(is.finite(Xc))) return(1)
  S_new <- max(1, max(Xc[is.finite(Xc)]) - min(Xc[is.finite(Xc)]))
  S_fit / S_new
}


#' Scale the boundary factor of the gp() term in a brmsfit's formula
#'
#' brms re-reads \code{c} from the gp() term of \code{$formula} at every
#' predict call, so a copy of the fit with a rescaled \code{c} is how
#' \code{predict()} holds the boundary on brms < 2.23.0; see
#' \code{.pin_gp_boundary_rows()}.  Only the terms of the right-hand-side sum
#' are searched, which is where \code{fit_bayesian_spatial_model()} puts it.
#'
#' @param model_obj A \code{brmsfit}, or anything with \code{$formula$formula}.
#' @param c_scale Positive factor.
#' @return \code{model_obj} with the gp() term's \code{c} multiplied by
#'   \code{c_scale}, or \code{NULL} when there is no such term.
#' @keywords internal
#' @noRd
.scale_gp_c <- function(model_obj, c_scale) {
  found <- FALSE
  rw <- function(e) {
    if (is.call(e) && identical(e[[1L]], as.name("gp")) && !is.null(e[["c"]])) {
      e[["c"]] <- eval(e[["c"]], baseenv()) * c_scale
      found <<- TRUE
    } else if (is.call(e) && identical(e[[1L]], as.name("+"))) {
      for (i in seq_along(e)[-1L]) e[[i]] <- rw(e[[i]])
    }
    e
  }
  f <- model_obj$formula$formula
  if (!inherits(f, "formula") || length(f) != 3L) return(NULL)
  f[[3L]] <- rw(f[[3L]])
  if (!found) return(NULL)
  model_obj$formula$formula <- f
  model_obj
}


#' Refuse a per-category posterior_epred() with a message that says why
#'
#' For an ordinal or categorical family \code{brms::posterior_epred()} returns
#' a draws x rows x categories \emph{array}: a probability per category, not
#' one expected value per row.  \code{predict.bayesian_fit()} took anything
#' that was not a matrix for a failed draw, so a real \code{cumulative()} fit
#' returned all-\code{NA} predictions under "posterior draw failed", and
#' \code{fitted()} said only that it got an array where it wanted a matrix --
#' though the fit's own documentation listed ordinal families as supported.
#'
#' @param draws What \code{posterior_epred()} returned.
#' @param engine The \code{brmsfit}, to name its family.
#' @param caller Method name for the message.
#' @param hint What to do instead, appended to the message.
#' @return \code{NULL}, invisibly, when \code{draws} is not a 3-D array.
#' @keywords internal
#' @noRd
.stop_if_category_epred <- function(draws, engine, caller, hint) {
  if (!is.array(draws) || length(dim(draws)) != 3L) return(invisible(NULL))
  fam <- .brms_family_name(tryCatch(engine$family, error = function(e) NULL))
  stop(sprintf(paste0("%s(): %s gives a probability per response category, ",
                      "so brms::posterior_epred() returned a %s array (draws x ",
                      "rows x categories), not one expected value per row. %s"),
               caller,
               if (is.na(fam)) "this fit's family"
               else sprintf("the %s family", sQuote(fam)),
               paste(dim(draws), collapse = " x "), hint), call. = FALSE)
}

# What predict.bayesian_fit() offers instead for a category family: the
# posterior predictive draws, which go through the method's own newdata
# pipeline.  A draw of Y is category k with the posterior mean probability of
# k, so the share of draws in k estimates exactly what posterior_epred()
# would average to.
.category_draws_hint <- paste0(
  "Use type = \"predict\", draws = TRUE for posterior draws of the predicted ",
  "category (as category indices); the share of draws in each category ",
  "estimates its probability. brms::posterior_epred(<fit>$engine) gives the ",
  "probabilities for the training rows.")


#' Predict from a Bayesian spatial GP model
#'
#' @section The GP boundary is held at its fitted value:
#' brms 2.17 to 2.22 do not store the Hilbert-space boundary \eqn{L} in a
#' fitted GP basis, so \code{brms:::.data_gp()} recomputes it from whatever
#' rows \code{predict()} is handed, which moved every eigenfunction of the
#' approximation with the newdata bounding box while the fitted basis
#' coefficients stayed put.  Two synthetic rows at the training coordinate
#' extrema are therefore appended before the posterior draw and dropped from
#' the result, and when \code{newdata} reaches past the training range the
#' \code{c} of the \code{gp()} term is scaled down by as much as the range
#' grew, so brms rebuilds exactly the boundary the model was fitted with.  A
#' prediction therefore does not depend on which other rows share the call:
#' chunked, fold-wise and single-call predictions agree, and
#' \code{\link{predict_surface}()} does not depend on \code{chunk_size}.
#' brms 2.23.0 and later store \eqn{L} and reuse it, so there \code{c} is left
#' alone and the two extra rows change nothing.
#'
#' Predictions outside the training coordinate envelope are extrapolation (a
#' notice is written to the log).  A row further than \eqn{L} from the centre
#' of the training coordinates on either axis is past the edge of the basis,
#' where the approximate GP is an odd reflection of the fitted surface rather
#' than an estimate of anything, so it is returned as \code{NA} (a column of
#' \code{NA} with \code{draws = TRUE}) with a warning.  The default boundary
#' factor puts that edge well outside the training data, so only
#' \code{newdata} reaching far past it is affected.
#'
#' @section Levels of a categorical predictor:
#' A factor or character predictor is handed to brms with the levels the model
#' was fitted with: the values that occur in the training rows
#' (\code{object$info$predictor_levels}).  brms compares the level \emph{set}
#' of a \code{newdata} factor with those, so without this a factor that only
#' listed an extra level, or a model fitted on a subset whose factor kept an
#' unused level, was refused outright.  A row whose value is not among the
#' fitted levels has no prediction: it is returned as \code{NA}, the other
#' rows are predicted, and one warning of class \code{spatialkit_new_levels}
#' names the predictor, the levels and the number of rows.  A prediction
#' therefore does not depend on which other rows share the call, and
#' \code{\link{predict_surface}()} returns the same surface at any
#' \code{chunk_size}.
#'
#' @description
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
#'   \code{n_draws x nrow(newdata)} matrix when \code{draws = TRUE}.  Rows
#'   with a missing or non-finite predictor or coordinate, and rows carrying
#'   a level of a categorical predictor the model was not fitted with, are
#'   \code{NA} (an \code{NA} column with \code{draws = TRUE}); the other
#'   rows are unaffected.  If the
#'   posterior draw fails the result is all \code{NA} (a 1-row matrix for
#'   \code{draws = TRUE}) with a warning of class
#'   \code{spatialkit_predict_failed} that gives the cause.  An ordinal or categorical
#'   family is not a failed draw and is an error under
#'   \code{type = "epred"}: its expected value is a probability per response
#'   category, not one number per row.  Use \code{type = "predict",
#'   draws = TRUE} for posterior draws of the predicted category, as category
#'   indices; the share of draws in each category estimates its probability.
#'   Without \code{draws = TRUE}, \code{type = "predict"} returns the mean (or
#'   median) category index, an expected rank for an ordinal family and an
#'   error for \code{brms::categorical()}, whose categories have no order.
#'   With \code{newdata = NULL} the cached \code{fitted()} values are returned
#'   only for the default \code{summary = "mean"}, \code{type = "epred"},
#'   \code{draws = FALSE} combination; any other combination is recomputed
#'   against the training data, because the cache holds epred column means and
#'   nothing else.
#' @family methods on a fitted model
#' @export
predict.bayesian_fit <- function(object, newdata = NULL,
                                 summary = c("mean", "median"),
                                 type = c("epred", "predict"),
                                 draws = FALSE, ...) {
  .check_dots(list(...), "predict.bayesian_fit")
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

  # ---- Early validation: the shared newdata contract ----
  # Done BEFORE the brms availability check: this is a structural check on
  # newdata that needs no backend, and it gives the most informative error.
  .check_predict_newdata(newdata, object, "predict.bayesian_fit")

  if (!requireNamespace("brms", quietly = TRUE))
    stop("predict.bayesian_fit(): package 'brms' is required.", call. = FALSE)

  model_obj <- object$engine
  if (!inherits(model_obj, "brmsfit"))
    stop("predict.bayesian_fit(): engine is not a brmsfit object.", call. = FALSE)

  # type = "predict" draws category INDICES for a category family.  Their
  # mean (or median) is an expected rank for an ordinal family, but nothing
  # at all for brms::categorical(), whose categories have no order: it came
  # back as 1.46, 1.97, 1.48, ... with no word.
  if (type == "predict" && !isTRUE(draws) &&
      identical(.brms_family_name(tryCatch(model_obj$family,
                                           error = function(e) NULL)),
                "categorical"))
    stop(paste0("predict.bayesian_fit(): the 'categorical' family's categories ",
                "have no order, so the ", summary, " of the predicted category ",
                "indices is not a prediction. ", .category_draws_hint),
         call. = FALSE)

  # ---- Preprocessing: match the pipeline used during fitting ----
  # Ensure newdata is in the same projected CRS that was used for training,
  # BEFORE prep_model_data(): prep's own ensure_projected() has no target and
  # would leave already-projected newdata in *its* CRS.  This mirrors
  # predict.gwr_fit() and predict.rf_fit().
  n_orig       <- nrow(newdata)
  training_crs <- sf::st_crs(object$data_sf)
  newdata      <- .replay_crs_assumption(newdata, object$data_sf, "predict.bayesian_fit")
  newdata      <- ensure_projected(newdata, target_crs = .crs_or_null(training_crs))

  # ---- Levels the model was not fitted with ----
  # One such row used to turn every prediction of the call into NA, with a
  # log line and no R condition: brms refuses the whole newdata, and the
  # failed draw was returned as all NA.  The rows are set aside here, the
  # rest are predicted, and one warning says which levels and how many rows.
  # Setting the value to NA lets the row cleaning below drop them, so they
  # come back NA through the same path as any incomplete row.  They must not
  # reach brms as NA: it predicts the reference level for an NA factor value
  # (1.07 where the row's own level gave 3.11).
  lv_chk <- .bayes_unseen_levels(object, newdata)
  if (any(lv_chk$unseen)) {
    for (cn in names(lv_chk$bad)) newdata[[cn]][lv_chk$bad[[cn]]] <- NA
    .warn_classed(
      "spatialkit_new_levels",
      paste0("predict.bayesian_fit(): %d of %d row(s) of `newdata` carry a ",
             "level the model was not fitted with (%s). Their predictions ",
             "are NA; the other %d row(s) are predicted as usual."),
      sum(lv_chk$unseen), n_orig, paste(lv_chk$named, collapse = "; "),
      n_orig - sum(lv_chk$unseen))
  }

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
                     "as missing, non-finite or carrying an unfitted level; ",
                     "returning %d NA(s)."), n_orig)
    return(if (draws) matrix(NA_real_, nrow = 1L, ncol = n_orig)
           else rep(NA_real_, n_orig))
  }

  # Build prediction data.frame with scaled coordinates & standardised predictors
  pred_df <- .prepare_brms_pred_df(object, newdata)

  # Draw from posterior.  The two padding rows and, on brms < 2.23.0, the
  # rescaled gp(c = ) hold the GP boundary at the one the model was fitted
  # with -- see .pin_gp_boundary_rows() -- and the padding columns are dropped
  # again immediately.  brms >= 2.23.0 reuses the fitted boundary itself.
  pinned <- .pin_gp_boundary_rows(object, pred_df)
  if (pinned$c_scale < 1 && utils::packageVersion("brms") < "2.23.0") {
    held <- .scale_gp_c(model_obj, pinned$c_scale)
    if (is.null(held))
      .warn_and_log(paste0("predict.bayesian_fit(): `newdata` widens the GP ",
                           "boundary and the engine's formula has no gp(c = ) ",
                           "term to hold it with, so these predictions depend ",
                           "on which rows share the call."))
    else model_obj <- held
  }
  if (any(pinned$beyond))
    .warn_and_log(paste0("predict.bayesian_fit(): %d of %d row(s) of `newdata` ",
                         "lie beyond the GP boundary the model was fitted with ",
                         "(further than L = %.3g scaled units from the centre ",
                         "of the training coordinates), where the approximate ",
                         "GP means nothing; they are returned as NA."),
                  sum(pinned$beyond), length(pinned$beyond),
                  object$info$gp_c * object$info$gp_S)
  draw_fn <- if (type == "epred") brms::posterior_epred else brms::posterior_predict
  draw_mat <- try(draw_fn(model_obj, newdata = pinned$df), silent = TRUE)
  if (is.matrix(draw_mat) && pinned$n_pad > 0L &&
      ncol(draw_mat) == nrow(pinned$df))
    draw_mat <- draw_mat[, seq_len(ncol(draw_mat) - pinned$n_pad), drop = FALSE]
  # An ordinal or categorical epred is a draws x rows x categories array; the
  # padding rows are dropped from it too, so the message below counts the
  # caller's rows (it said "150 x 7 x 3" for five).
  if (is.array(draw_mat) && length(dim(draw_mat)) == 3L && pinned$n_pad > 0L &&
      dim(draw_mat)[2L] == nrow(pinned$df))
    draw_mat <- draw_mat[, seq_len(dim(draw_mat)[2L] - pinned$n_pad), ,
                         drop = FALSE]
  if (is.matrix(draw_mat) && ncol(draw_mat) == length(pinned$beyond))
    draw_mat[, pinned$beyond] <- NA_real_

  # Not a failed draw: an ordinal or categorical family's epred.  Raised, not
  # returned as NA, because no retry will produce one number per row.  The
  # probabilities for new rows are not offered through
  # brms::posterior_epred(<fit>$engine, newdata = ): the engine needs the
  # scaled ..x/..y (and standardised predictors) this method builds, and it
  # refused the user's newdata.  The share of predicted-category draws in
  # each category estimates the same posterior mean probability.
  if (type == "epred")
    .stop_if_category_epred(draw_mat, model_obj, "predict.bayesian_fit",
                            hint = .category_draws_hint)
  if (inherits(draw_mat, "try-error") || !is.matrix(draw_mat)) {
    # A real warning with a class, not only a log line: every prediction of
    # the call is NA, and a log line is invisible to tryCatch(),
    # options(warn = 2), knitr and spatialkit_quiet().
    .warn_classed(
      "spatialkit_predict_failed",
      paste0("predict.bayesian_fit(): the posterior draw failed, so all %d ",
             "prediction(s) are NA. Cause: %s"),
      n_orig,
      if (inherits(draw_mat, "try-error")) .try_error_message(draw_mat)
      else sprintf("brms returned a %s, not a draws x rows matrix",
                   class(draw_mat)[1L]))
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
#' column is present.  These are \strong{in-sample} values (each observation
#' was inside its own bandwidth window), so \code{summary()} on a
#' \code{gwr_fit} reports an optimistic fit.  Use \code{\link{cv_gwr}} for a
#' spatially blocked estimate.
#'
#' @param object A \code{gwr_fit}.
#' @param ... Ignored.
#' @return Numeric vector of length \code{object$n} (\code{NA} where extraction
#'   failed).
#' @family methods on a fitted model
#' @export
fitted.gwr_fit <- function(object, ...) {
  .extract_gwr_values(
    object$engine, object$data_sf, object$formula,
    object$n, object$response_var, mode = "basic"
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
#' \code{predict(object, type = "predict")} recompute it from scratch.
#' Call \code{\link{clear_fitted_cache}} if the engine has been mutated by hand
#' after fitting.
#'
#' @section The cache is shared by copies, and validated:
#' An environment has reference semantics, which is what makes the memo survive
#' R's copy-on-modify.  But it also means \code{fit2 <- fit} gives the two
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
#' share (harmless, since the other simply recomputes), and \code{identical()}
#' cannot distinguish two fits by their caches.  The digest covers
#' \code{data_sf} only.  The entry is also tied to the engine that computed
#' it -- a refit or \code{update()} of the \code{brmsfit} is a different
#' sampling run and recomputes -- but a \code{brmsfit} edited by hand in place
#' is what \code{\link{clear_fitted_cache}} is for.  The entry holds only the
#' values and a small identifier of the sampling run, so a fit saved with
#' \code{saveRDS()} after \code{fitted()} is no larger for it.
#'
#' @param object A \code{bayesian_fit}.
#' @param ... Ignored.
#' @return Numeric vector of length \code{object$n}.  A posterior that cannot
#'   be drawn is an error, as is a family with a probability per response
#'   category (ordinal, categorical), which has no single fitted value per row.
#' @family methods on a fitted model
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
    # `hit$engine` is the identity check the data stamp cannot make.  The key
    # digests n, the attribute table and the coordinates, so two fits over the
    # SAME data hash identically however different their engines are -- and
    # because the cache is shared by every copy of a fit, `refit <- fit;
    # refit$engine <- <re-estimated>` then read the original engine's fitted
    # values back out.  It holds .fitted_engine_token(), not necessarily the
    # engine itself: see there for why holding the engine doubled saveRDS().
    # identical() is cheap here: for the common case it is the same object,
    # which R settles by pointer.
    if (is.list(hit) && identical(hit$n, object$n) &&
        identical(hit$key, key) &&
        identical(hit$engine, .fitted_engine_token(object$engine)) &&
        is.numeric(hit$values) && length(hit$values) == object$n)
      return(hit$values)
    # Stale: a copy carrying different data or a different engine, or an entry
    # written by an older version of this package.  It is not returned -- the
    # test above already refused it -- and it is not deleted either: the cache
    # is shared by every copy of a fit, so deleting here threw away an entry
    # that is still valid for the fit that wrote it, and a recompute that then
    # failed left that fit with nothing.  The write below replaces it.
  }

  if (!requireNamespace("brms", quietly = TRUE))
    stop("fitted.bayesian_fit(): package 'brms' is required.", call. = FALSE)

  model_obj <- object$engine
  pred_df   <- .prepare_brms_pred_df(object, object$data_sf)

  draws <- try(brms::posterior_epred(model_obj, newdata = pred_df), silent = TRUE)
  # An error here is an error: a posterior that cannot be drawn used to come
  # back as all-NA fitted values with nothing said, and summary() and
  # model_metrics() then reported n = 0 as though the data were missing.
  # (predict.bayesian_fit() differs on purpose: it documents an all-NA result
  # for a failed draw on newdata, and logs the cause.)
  if (inherits(draws, "try-error"))
    stop(sprintf(paste0("fitted.bayesian_fit(): brms::posterior_epred() failed ",
                        "on the training data: %s"),
                 conditionMessage(attr(draws, "condition"))), call. = FALSE)
  .stop_if_category_epred(draws, model_obj, "fitted.bayesian_fit",
                          hint = paste0("fitted(), residuals(), summary() and ",
                                        "model_metrics() need one number per ",
                                        "row; use predict(<fit>, type = ",
                                        "\"predict\", draws = TRUE) for ",
                                        "predicted categories, or ",
                                        "brms::posterior_epred(<fit>$engine) ",
                                        "for the probabilities."))
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
           list(n = object$n, key = key,
                engine = .fitted_engine_token(object$engine),
                values = fitted_vals),
           envir = cache)
  }

  fitted_vals
}


#' What a fitted() cache entry keeps to tell one engine from another
#'
#' The entry must belong to the engine that produced it (see
#' \code{fitted.bayesian_fit}), and it used to hold the engine itself.  That
#' costs nothing in memory, but \code{serialize()} tracks environments by
#' reference and lists not at all, so \code{saveRDS()} on a fit whose cache was
#' warm -- after any \code{summary()}, \code{residuals()} or
#' \code{model_metrics()} -- wrote the whole brmsfit a second time: a 40 MB
#' engine saved as 80 MB before \code{fitted()} and 120 MB after.
#'
#' A brmsfit carries a stanfit, and every stanfit carries an environment
#' (\code{@.MISC}) that rstan's sampler, and brms's reader for CmdStan output,
#' create afresh for each run.  Copies of the engine share it, a refit or an
#' \code{update()} gets a new one, and it is serialised once however many
#' references point at it, so it tells engines apart as well as the engine
#' does and adds nothing to a saved fit.  Any other engine (a custom backend,
#' a test double) is its own token, as before.
#'
#' @param engine The fit's \code{$engine}.
#' @return An environment, or \code{engine}.
#' @keywords internal
#' @noRd
.fitted_engine_token <- function(engine) {
  misc <- tryCatch(engine$fit@.MISC, error = function(e) NULL)
  if (is.environment(misc)) misc else engine
}


#' Cheap fingerprint of the training data a cached fitted() was computed from
#'
#' Covers exactly what \code{.prepare_brms_pred_df()} reads: the attribute
#' columns and the coordinates.  \code{$engine} is deliberately excluded:
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
#' \code{brmsfit} engine has been manually mutated after fitting.  A change to
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
#' @family package options and caches
#' @examples
#' # Only a bayesian_fit carries the cache; on any other fit this is a no-op.
#' if (requireNamespace("ranger", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   pts <- st_as_sf(
#'     data.frame(x = 5e5 + runif(60, 0, 1000), y = 5e6 + runif(60, 0, 1000),
#'                a = rnorm(60)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   pts$z <- 2 * pts$a + rnorm(60, 0, 0.3)
#'   fit <- fit_rf_model(pts, "z", "a", num_trees = 50, seed = 1)
#'   clear_fitted_cache(fit)
#' }
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
#' @family methods on a fitted model
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
#' @family methods on a fitted model
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
#' Returns the whole surface of coefficients, one row per observation and one
#' column per term, against the single global vector \code{coef()} returns
#' for an \code{lm}.  That table is the point of fitting a GWR at all:
#' inspect the spread of a predictor's column to see where, and by how much,
#' its relationship with the response changes across the study area, and join
#' it back to \code{object$data_sf} to map it.  Use
#' \code{\link{plot.spatial_fit}()} for a quick look at that map.
#'
#' @section What is and is not returned:
#' Only the model terms: the intercept and one column per predictor.
#' The engine's \code{SDF} data slot carries a good deal more alongside them
#' (standard errors, t-values, the observed response, the fitted values, the
#' residuals, \code{Local_R2}): 15 columns for a two-predictor fit, of which 3
#' are coefficients.  Returning the whole slot would have made
#' \code{coef(fit)$Local_R2} and \code{coef(fit)$a_SE} read like coefficients
#' and \code{ncol(coef(fit))} a meaningless number.  Ask for the standard
#' errors and t values with \code{se = TRUE}; reach for
#' \code{object$engine$SDF} when you want the rest.
#'
#' A coefficient is \code{NA} at a location whose local window cannot
#' estimate it: a predictor that does not vary among the observations that
#' carry weight there has no slope, and the intercept is then not identified
#' either unless the predictor is 0 there.  \code{object$info$nonfinite_coef}
#' is the mask, and \code{\link{fit_gwr_model}} raises a warning that counts
#' them.
#'
#' If the model terms cannot be located in the \code{SDF} (an engine that
#' names its coefficient columns differently), the whole slot is returned with
#' a warning saying so.  The call neither errors nor returns a silently short
#' table.
#'
#' @section Standard errors:
#' With \code{se = TRUE} the local standard errors and t values come with the
#' estimates.  A standard error is
#' \eqn{\hat\sigma\sqrt{(C C^\top)_{jj}}}{sigma * sqrt(diag(C C'))}, with
#' \eqn{C = (X^\top W X)^{-1} X^\top W} the weights the local regression puts
#' on the responses and \eqn{\hat\sigma^2} the residual sum of squares over
#' the fit's effective degrees of freedom: the definition
#' \code{GWmodel::gwr.basic()} uses, computed from the QR decomposition that
#' gave the estimate.  They treat the bandwidth as known.  The t values are
#' one test per term and location, tens to thousands of them and strongly
#' dependent, so a map of \eqn{|t| > 2} overstates what is significant; no
#' p-value, adjusted or not, is returned.  What an adjustment for the number
#' of tests needs is on the fit: the t values here, the effective number of
#' parameters \code{object$info$enp} and the residual degrees of freedom
#' \code{object$info$edf}.  When observations share a location the standard
#' errors are too small as well (see \code{\link{fit_gwr_model}},
#' \strong{Observations that share a location}).
#'
#' @param object A \code{gwr_fit} object.
#' @param se Logical; also return the standard errors and t values.  Default
#'   \code{FALSE}.
#' @param ... Ignored.
#' @return With \code{se = FALSE}, a data.frame of local coefficient
#'   estimates: one row per observation, one column per model term.  With
#'   \code{se = TRUE}, a list of three such data.frames: \code{estimate},
#'   \code{se} and \code{t}.
#'   Never \code{NULL}: when the engine carries no \code{SDF} component this
#'   errors, following the \code{coef()} contract described in
#'   \code{\link{new_spatial_fit}}.
#' @family methods on a fitted model
#' @export
coef.gwr_fit <- function(object, se = FALSE, ...) {
  if (!is.logical(se) || length(se) != 1L || is.na(se))
    stop("coef.gwr_fit(): `se` must be TRUE or FALSE.", call. = FALSE)
  sdf <- tryCatch(object$engine$SDF, error = function(e) NULL)
  # Returning NULL here would be indistinguishable from "this GWR has no
  # coefficients", which is never true.  See ?new_spatial_fit for the contract.
  if (is.null(sdf))
    stop("coef.gwr_fit(): the fit carries no `SDF` component, so the ",
         "local coefficients cannot be extracted.", call. = FALSE)
  dat <- if (inherits(sdf, "Spatial")) sdf@data else sf::st_drop_geometry(sdf)

  # Keep only the model terms.  The SDF is
  #   c(colnames(betas), "y", "yhat", "residual", "CV_Score", "Stud_residual",
  #     paste0(colnames(betas), "_SE"), paste0(colnames(betas), "_TV"),
  #     "Local_R2")
  # so the coefficients are the leading block and everything after it is a
  # diagnostic.  Match on the model-matrix column names (normalised, so
  # "(Intercept)" and the SDF's "Intercept" agree) and take the FIRST unclaimed
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
    if (length(pos) == length(want)) {
      est <- dat[, pos, drop = FALSE]
      if (!se) return(est)
      # The standard errors and t values sit in two blocks of the same width
      # after the five diagnostic columns; they are taken by position, checked
      # against their names, so a predictor called `a_SE` cannot stand in.
      k    <- length(pos)
      i_se <- k + 5L + seq_len(k)
      i_tv <- 2L * k + 5L + seq_len(k)
      laid_out <- identical(pos, seq_len(k)) && ncol(dat) >= 3L * k + 5L &&
        identical(names(dat)[i_se], paste0(names(dat)[pos], "_SE")) &&
        identical(names(dat)[i_tv], paste0(names(dat)[pos], "_TV"))
      if (!laid_out)
        stop("coef.gwr_fit(): the fit's `SDF` does not carry standard errors ",
             "and t values in the layout fit_gwr_model() writes ",
             "(<term>_SE and <term>_TV after the coefficients), so ",
             "`se = TRUE` cannot be answered.", call. = FALSE)
      se_df <- dat[, i_se, drop = FALSE]
      tv_df <- dat[, i_tv, drop = FALSE]
      names(se_df) <- names(est)
      names(tv_df) <- names(est)
      return(list(estimate = est, se = se_df, t = tv_df))
    }
    .log_warn(paste0("coef.gwr_fit(): matched only %d of %d model term(s) to a ",
                     "column of the engine's SDF, so the full SDF data slot is ",
                     "returned instead of the coefficients alone. Columns such ",
                     "as *_SE, *_TV, residual and Local_R2 are diagnostics, ",
                     "not coefficients."),
              length(pos), length(want))
  }
  if (se)
    stop("coef.gwr_fit(): the model terms could not be located in the fit's ",
         "`SDF`, so `se = TRUE` cannot be answered.", call. = FALSE)
  dat
}


#' Extract Bayesian model fixed-effect summaries
#'
#' Returns the posterior summary of the global (non-spatial) regression terms:
#' estimate, error and credible interval per predictor, as
#' \code{brms::fixef()} reports them.  Reach for it to read the average effect
#' of a predictor with its uncertainty attached (the Bayesian counterpart to
#' a coefficient table), remembering that the Gaussian-process term has
#' already absorbed the spatially structured part of the signal, so these are
#' effects net of location.
#'
#' @section Standardised predictors:
#' The summaries are on the scale the model was fitted on.  A fit made with
#' \code{standardize_predictors = TRUE} was fitted on centred and scaled
#' numeric predictors, so each slope is the change in the linear predictor per
#' \emph{standard deviation} of its predictor and the intercept is its value
#' at the predictor \emph{means}, not the raw-unit numbers \code{stats::lm()}
#' reports on the same formula.  Nothing on the returned matrix says so;
#' \code{print()} on the fit does, and the centre and scale of each predictor
#' are in \code{object$info$predictor_scaling}.  To put a slope back in raw
#' units divide its \code{Estimate}, \code{Est.Error} and interval bounds by
#' that predictor's \code{scale}.  The intercept's \code{Estimate} follows by
#' linearity (subtract each raw-unit slope times its predictor's
#' \code{center}), but its \code{Est.Error} and interval depend on the
#' posterior covariance of the coefficients: transform the draws from
#' \code{brms::as_draws_df(object$engine)} for those, or refit without
#' standardising.
#'
#' @param object A \code{bayesian_fit} object.
#' @param ... Ignored.
#' @return A matrix of fixed-effect posterior summaries, as returned by
#'   \code{brms::fixef()}, on the fitted scale (per standard deviation of each
#'   predictor under \code{standardize_predictors = TRUE}; see above).  Never
#'   \code{NULL}: a missing 'brms' or a failing
#'   \code{fixef()} call errors, following the \code{coef()} contract described
#'   in \code{\link{new_spatial_fit}}.
#' @family methods on a fitted model
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


#' One `newdata` contract for every `predict()` method
#'
#' The three methods documented one contract and enforced three.  A predictor
#' that arrived as character was refused by name by \code{rf_fit} and
#' \code{bayesian_fit} and returned all-\code{NA} with a generic backend
#' warning from \code{gwr_fit}; a bare \code{sfc} died inside two of them with
#' R's "argument must be coercible to non-negative integer"; a missing column
#' was reported by two different functions in two wordings.  This runs first
#' in all three, so the message is the same whichever fit is behind it.
#'
#' @param newdata Whatever the caller passed.
#' @param object The \code{spatial_fit}.
#' @param caller The method name for the message.
#' @return \code{newdata}, invisibly.
#' @keywords internal
#' @noRd
.check_predict_newdata <- function(newdata, object, caller) {
  if (inherits(newdata, "sfc"))
    stop(caller, "(): `newdata` is a bare geometry column (sfc) and carries no ",
         "predictor columns; pass an sf object with the ",
         "predictors ", paste(sQuote(object$predictor_vars), collapse = ", "),
         ".", call. = FALSE)
  if (!inherits(newdata, "sf"))
    stop(caller, "(): `newdata` must be an sf object -- the prediction ",
         "locations are its geometry; got ",
         paste(class(newdata), collapse = "/"),
         if (is.data.frame(newdata))
           ". Build one with sf::st_as_sf(newdata, coords = c(\"x\", \"y\"), crs = )"
         else "", ".", call. = FALSE)
  pv <- object$predictor_vars
  miss <- setdiff(pv, names(newdata))
  if (length(miss))
    stop(caller, "(): `newdata` is missing predictor column(s) ",
         paste(sQuote(miss), collapse = ", "), ".", call. = FALSE)
  # Type agreement with the training data.  A CSV round-trip that turned a
  # numeric predictor into character must be named, not silently coerced or
  # passed on to a backend that fails for its own reasons.
  tr <- tryCatch(sf::st_drop_geometry(object$data_sf), error = function(e) NULL)
  if (!is.null(tr)) {
    for (cn in intersect(pv, names(tr))) {
      a <- tr[[cn]]; b <- newdata[[cn]]
      if ((is.numeric(a) || is.logical(a)) && !(is.numeric(b) || is.logical(b)))
        stop(sprintf(paste0("%s(): predictor '%s' was %s at fit time but is %s ",
                            "in `newdata`. Convert it (as.numeric() after ",
                            "checking for a stray thousands separator or ",
                            "\"NA\" string) rather than letting it be coerced."),
                     caller, cn, if (is.logical(a)) "logical" else "numeric",
                     class(b)[1L]), call. = FALSE)
    }
  }
  invisible(newdata)
}
