# =============================================================================
# Random forest backend (ranger)
# =============================================================================


#' Build the predictor frame a ranger forest expects
#'
#' Fitting and prediction must produce identical column names in identical
#' order, so both go through this.
#'
#' @param .caller Name used in the error message; the same helper is reached
#'   from \code{fit_rf_model()} and from \code{predict.rf_fit()}, and naming
#'   the wrong one sends the reader to the wrong argument.
#' @keywords internal
#' @noRd
.rf_frame <- function(data_sf, predictor_vars, include_coords = FALSE,
                      .caller = "fit_rf_model") {
  df <- sf::st_drop_geometry(data_sf)
  missing_v <- setdiff(predictor_vars, names(df))
  if (length(missing_v) > 0L)
    stop(.caller, "(): predictor(s) absent from the data: ",
         paste(sQuote(missing_v), collapse = ", "), call. = FALSE)
  out <- df[, predictor_vars, drop = FALSE]
  if (isTRUE(include_coords)) {
    crd <- sf::st_coordinates(data_sf)
    out[["..x"]] <- crd[, 1L]
    out[["..y"]] <- crd[, 2L]
  }
  out
}


#' Coerce a ranger scalar diagnostic to a single finite number or NA
#'
#' \code{ranger} can leave \code{prediction.error} / \code{r.squared} unset
#' (for instance when \code{oob.error = FALSE} is forwarded through
#' \code{...}).  \code{as.numeric(NULL)} is \code{numeric(0)}, which
#' \code{\%||\%} does not catch because it is not \code{NULL}: the downstream
#' \code{if (is.finite(x))} then raises "argument is of length zero", and
#' \code{sprintf()} handed a length-zero argument prints nothing at all.
#'
#' @param x Any value.
#' @return A length-one numeric: the value if it is a single finite number,
#'   \code{NA_real_} otherwise.
#' @keywords internal
#' @noRd
.num1 <- function(x) {
  v <- suppressWarnings(as.numeric(x))
  if (length(v) != 1L || !is.finite(v)) return(NA_real_)
  v
}


#' Check that categorical predictor levels in newdata were seen at fit time
#'
#' \code{ranger} matches factor predictors by level, so a level that was not
#' present when the forest was grown has no split to follow.  Depending on the
#' ranger version this either errors -- which \code{predict.rf_fit()}'s
#' \code{tryCatch} would turn into an all-\code{NA} vector plus one log line --
#' or, as in ranger 0.16, silently returns a plausible-looking number.  Both
#' are worse than an error naming the level.
#'
#' @param X Prediction frame from \code{.rf_frame()}.
#' @param train_X Training frame from \code{.rf_frame()}.
#' @return \code{X}, with categorical columns coerced to the training levels.
#' @keywords internal
#' @noRd
.rf_align_levels <- function(X, train_X) {
  for (cn in intersect(names(X), names(train_X))) {
    tr <- train_X[[cn]]
    if (!(is.factor(tr) || is.character(tr))) next
    lv  <- if (is.factor(tr)) levels(tr) else sort(unique(as.character(tr)))
    val <- as.character(X[[cn]])
    unseen <- setdiff(unique(val[!is.na(val)]), lv)
    if (length(unseen) > 0L)
      stop(sprintf(paste0("predict.rf_fit(): predictor '%s' in newdata has ",
                          "level(s) the forest was never grown with: %s. ",
                          "Trained levels: %s."),
                   cn, paste(sQuote(unseen), collapse = ", "),
                   paste(sQuote(lv), collapse = ", ")), call. = FALSE)
    # Coerce with the TRAINING levels so the codes ranger sees are the codes
    # it was built with, whatever order newdata happens to carry.
    X[[cn]] <- if (is.factor(tr)) factor(val, levels = lv) else val
  }
  X
}


# -----------------------------------------------------------------------------
# Exported
# -----------------------------------------------------------------------------

#' Fit a random forest via ranger
#'
#' Fits a regression random forest on an \code{sf} dataset and returns it as a
#' \code{spatial_fit}, so it works with \code{\link{cv_spatial}},
#' \code{\link{predict_surface}}, \code{\link{area_of_applicability}} and the
#' \code{plot()} method like any other backend.
#'
#' @section Coordinates are not predictors by default:
#' Handing a random forest the x and y coordinates lets it reproduce the
#' training surface almost exactly by memorising location, and then fail badly
#' anywhere it has not seen. Random cross-validation will not catch this --
#' nearby points leak between folds, so the memorised surface scores well --
#' which is how the practice became common. Meyer et al. (2019) show the
#' collapse directly. \code{include_coords} therefore defaults to
#' \code{FALSE}, and setting it to \code{TRUE} warns. If you do use it, score
#' the model with \code{\link{cv_spatial}} and blocked folds, never with the
#' out-of-bag error.
#'
#' @section The out-of-bag error is a random hold-out:
#' \code{ranger}'s OOB error holds each observation out of the trees that did
#' not sample it. That is a random hold-out, so under spatial autocorrelation
#' it is optimistic for exactly the reason random k-fold is: the trees that
#' "did not see" a point almost certainly saw its neighbours. It is reported
#' as \code{$info$oob_rmse} and \code{$info$oob_r_squared} and labelled as OOB
#' everywhere it appears. Use \code{\link{cv_rf}} for a spatial estimate.
#'
#' @section What fitted() returns:
#' \code{fitted()} on an \code{rf_fit} returns \strong{out-of-bag} predictions,
#' not in-sample ones, following the convention of the random forest packages
#' themselves. In-sample predictions from a forest are close to memorisation
#' and would make \code{summary()} report a fictitious R-squared. The
#' consequence is that \code{summary()} means something different here than for
#' a \code{gwr_fit} or \code{bayesian_fit}, whose fitted values are in-sample:
#' do not compare the two directly. \code{\link{compare_models_cv}} exists for
#' that.
#'
#' @param data_sf An sf object with response, predictors and geometry.
#' @param response_var Response column name.
#' @param predictor_vars Predictor column names.
#' @param num_trees Number of trees. Default 500.
#' @param mtry Predictors sampled per split. \code{NULL} uses ranger's default.
#' @param min_node_size Minimum node size. \code{NULL} uses ranger's default.
#' @param importance \code{"permutation"} (default), \code{"impurity"} or
#'   \code{"none"}. Impurity importance is biased toward continuous and
#'   high-cardinality predictors (Strobl et al. 2007), which is why permutation
#'   is the default despite costing more.
#' @param include_coords Add the coordinates as predictors. Default
#'   \code{FALSE}; see above.
#' @param seed Seed passed to ranger. Default 123.
#' @param num_threads Threads for ranger. \code{NULL} uses ranger's default.
#' @param .already_prepped Internal; skip \code{prep_model_data()} because the
#'   caller has already projected and cleaned the data.
#' @param ... Passed to \code{ranger::ranger()}.
#'
#' @return An \code{rf_fit} object (inherits from \code{spatial_fit}).
#'   \code{$info} carries \code{num_trees}, \code{mtry}, \code{min_node_size},
#'   \code{importance_type}, \code{importance} (a named numeric, or
#'   \code{NULL} when \code{importance = "none"}), \code{include_coords},
#'   \code{oob_rmse} and \code{oob_r_squared} (each \code{NA_real_} when
#'   ranger did not compute it -- forwarding \code{oob.error = FALSE} through
#'   \code{...} is one way to get there), \code{fitted_are_oob} (always
#'   \code{TRUE}; \code{summary()} reads it to label its metrics) and
#'   \code{seed}. The raw forest is in \code{$engine}.
#'
#' @references
#' Meyer, H., Reudenbach, C., Wöllauer, S. and Nauss, T. (2019). Importance of
#' spatial predictor variable selection in machine learning applications --
#' moving from data reproduction to spatial prediction. \emph{Ecological
#' Modelling} 411, 108815. \doi{10.1016/j.ecolmodel.2019.108815}
#'
#' Strobl, C., Boulesteix, A.-L., Zeileis, A. and Hothorn, T. (2007). Bias in
#' random forest variable importance measures: illustrations, sources and a
#' solution. \emph{BMC Bioinformatics} 8, 25. \doi{10.1186/1471-2105-8-25}
#'
#' @seealso \code{\link{cv_rf}} for a spatially blocked performance estimate,
#'   \code{\link{area_of_applicability}}, which can take
#'   \code{weights = pmax(fit$info$importance, 0)}.
#' @family model fitting
#' @examples
#' \donttest{
#' if (requireNamespace("ranger", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 150
#'   dat <- st_as_sf(
#'     data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
#'                a = rnorm(n), b = rnorm(n)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   dat$z <- 2 * dat$a - dat$b + rnorm(n, 0, 0.3)
#'   fit <- fit_rf_model(dat, "z", c("a", "b"))
#'   fit
#'   fit$info$importance
#' }
#' }
#' @export
fit_rf_model <- function(data_sf, response_var, predictor_vars,
                         num_trees = 500L, mtry = NULL, min_node_size = NULL,
                         importance = c("permutation", "impurity", "none"),
                         include_coords = FALSE, seed = 123L,
                         num_threads = NULL, .already_prepped = FALSE, ...) {
  if (!inherits(data_sf, "sf"))
    stop("fit_rf_model(): `data_sf` must be an sf object.", call. = FALSE)
  if (!requireNamespace("ranger", quietly = TRUE))
    stop("fit_rf_model(): package 'ranger' is required. Install with ",
         "install.packages('ranger').", call. = FALSE)
  importance <- match.arg(importance)
  # prep_model_data() accepts character(0) so an intercept-only spatial GP can
  # be fitted; a forest with nothing to split on is not a model, so reject it
  # here where the message can say why.
  if (!is.character(predictor_vars) || length(predictor_vars) < 1L)
    stop("fit_rf_model(): `predictor_vars` must name at least one predictor; ",
         "a forest has nothing to split on otherwise.", call. = FALSE)

  # Check the response type BEFORE prep_model_data(), which would otherwise
  # trip over a character column first and report something less useful.
  raw_y <- if (response_var %in% names(data_sf))
    sf::st_drop_geometry(data_sf)[[response_var]] else NULL
  if (!is.null(raw_y) && !is.numeric(raw_y))
    stop(sprintf(paste0("fit_rf_model(): response '%s' is not numeric. This ",
                        "fits a regression forest; use ranger::ranger() ",
                        "directly for classification."), response_var),
         call. = FALSE)

  dat <- if (isTRUE(.already_prepped)) data_sf else
    prep_model_data(data_sf = data_sf, response_var = response_var,
                    predictor_vars = predictor_vars, pointize = "auto")

  if (isTRUE(include_coords))
    .log_warn(paste0("fit_rf_model(): include_coords = TRUE. A forest given ",
                     "the coordinates can reproduce the training surface by ",
                     "memorising location and then fail wherever it has not ",
                     "been; random cross-validation will not detect this. ",
                     "Score the result with cv_rf() and blocked folds, not ",
                     "with the out-of-bag error (Meyer et al. 2019)."))

  y <- sf::st_drop_geometry(dat)[[response_var]]
  if (!is.numeric(y))
    stop(sprintf("fit_rf_model(): response '%s' is not numeric.", response_var),
         call. = FALSE)
  n_unique <- length(unique(y[is.finite(y)]))
  if (n_unique <= 2L)
    .log_warn(paste0("fit_rf_model(): response '%s' has %d distinct values. A ",
                     "regression forest estimates the conditional mean, which ",
                     "for a 0/1 response is a probability, but a probability ",
                     "forest (ranger::ranger(probability = TRUE)) is usually ",
                     "the better tool."), response_var, n_unique)

  X <- .rf_frame(dat, predictor_vars, include_coords)
  if (nrow(X) < 2L)
    stop("fit_rf_model(): fewer than two usable rows after cleaning.",
         call. = FALSE)

  fit <- tryCatch(
    ranger::ranger(x = X, y = y, num.trees = as.integer(num_trees),
                   mtry = mtry, min.node.size = min_node_size,
                   importance = importance, seed = seed,
                   num.threads = num_threads, ...),
    error = function(e)
      stop(sprintf("fit_rf_model(): ranger() failed: %s", conditionMessage(e)),
           call. = FALSE)
  )

  imp <- if (identical(importance, "none")) NULL else fit$variable.importance
  # .num1() collapses an unset ranger field to NA_real_ rather than the
  # numeric(0) that as.numeric(NULL) produces -- see ?.num1.
  oob_mse <- .num1(fit$prediction.error)

  new_spatial_fit(
    subclass       = "rf_fit",
    engine         = fit,
    # Show the coordinates in the formula when they are predictors, so
    # print()ing the fit does not hide them.  It is display-only: the forest is
    # built through ranger's x/y interface, never from this formula.
    formula        = stats::reformulate(
      termlabels = if (isTRUE(include_coords))
        c(predictor_vars, "..x", "..y") else predictor_vars,
      response = response_var),
    response_var   = response_var,
    predictor_vars = predictor_vars,
    data_sf        = dat,
    info = list(
      num_trees        = as.integer(num_trees),
      mtry             = fit$mtry,
      min_node_size    = fit$min.node.size,
      importance_type  = importance,
      importance       = imp,
      include_coords   = isTRUE(include_coords),
      oob_rmse         = if (is.finite(oob_mse)) sqrt(oob_mse) else NA_real_,
      oob_r_squared    = .num1(fit$r.squared),
      fitted_are_oob   = TRUE,
      seed             = seed
    )
  )
}


#' Cross-validate a random forest with spatial folds
#'
#' Thin wrapper over \code{\link{cv_spatial}} that refits a \code{ranger}
#' forest on each training fold. Unlike the out-of-bag error, this holds out
#' whole spatial blocks, so neighbours of a held-out point are not sitting in
#' the training set.
#'
#' @param data_sf An sf object.
#' @param response_var Response column name.
#' @param predictor_vars Predictor column names.
#' @param folds Optional fold definitions from \code{\link{make_folds}}.
#' @param k Number of folds when \code{folds} is \code{NULL}. Default 5.
#' @param seed RNG seed. Default 123.
#' @param parallel Passed to \code{\link{cv_spatial}}. Default \code{FALSE}.
#' @param block_size,auto_range,boundary Passed to \code{\link{cv_spatial}}.
#' @param ... Passed to \code{\link{fit_rf_model}} on every fold.
#' @return The \code{\link{cv_spatial}} result.
#' @family cross-validation
#' @examples
#' \donttest{
#' if (requireNamespace("ranger", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 200
#'   dat <- st_as_sf(
#'     data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), a = rnorm(n)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   dat$z <- 2 * dat$a + rnorm(n, 0, 0.3)
#'   cv_rf(dat, "z", "a", k = 4)$overall
#' }
#' }
#' @export
cv_rf <- function(data_sf, response_var, predictor_vars, folds = NULL, k = 5,
                  seed = 123, parallel = FALSE, block_size = NULL,
                  auto_range = FALSE, boundary = NULL, ...) {
  if (!requireNamespace("ranger", quietly = TRUE))
    stop("cv_rf(): package 'ranger' is required.", call. = FALSE)
  fit_fn <- function(train_sf)
    fit_rf_model(train_sf, response_var, predictor_vars,
                 .already_prepped = TRUE, ...)
  cv_spatial(data_sf, response_var, predictor_vars, fit_fn = fit_fn,
             folds = folds, k = k, seed = seed, boundary = boundary,
             block_size = block_size, auto_range = auto_range,
             parallel = parallel)
}


# -----------------------------------------------------------------------------
# S3 methods
# -----------------------------------------------------------------------------

#' Predict from a random forest fit
#'
#' With \code{newdata = NULL} this returns \strong{out-of-bag} predictions, not
#' in-sample ones -- see \code{\link{fit_rf_model}}.
#'
#' @param object An \code{rf_fit}.
#' @param newdata Optional sf object carrying the same predictors. It is
#'   transformed to the CRS used at fitting time first, so a forest that
#'   includes the coordinates is not fed a different coordinate system.
#'   Categorical predictors must not carry levels absent from the training
#'   data; an unseen level is an error, not a guess.
#' @param ... Passed to \code{ranger}'s predict method. Arguments that make
#'   \code{ranger} return a matrix rather than one value per row --
#'   \code{predict.all = TRUE}, \code{type = "quantiles"}, \code{type = "se"}
#'   with \code{predict.all} -- are rejected, because this method's contract is
#'   one number per row of \code{newdata}. Call
#'   \code{predict(fit$engine, data = ...)} directly for those.
#' @return Numeric vector, aligned to \code{nrow(newdata)} with \code{NA} for
#'   rows dropped as incomplete.
#' @export
predict.rf_fit <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(fitted(object))
  if (!requireNamespace("ranger", quietly = TRUE))
    stop("predict.rf_fit(): package 'ranger' is required.", call. = FALSE)

  n_orig  <- nrow(newdata)
  newdata <- ensure_projected(newdata, target_crs = .crs_or_null(object$data_sf))
  newdata$..orig_row_id.. <- seq_len(n_orig)
  newdata <- prep_model_data(newdata, object$response_var,
                             object$predictor_vars, pointize = "auto",
                             require_response = FALSE)
  clean_idx <- seq_len(n_orig) %in% newdata$..orig_row_id..
  newdata$..orig_row_id.. <- NULL
  n_new <- nrow(newdata)

  include_coords <- isTRUE(object$info$include_coords)
  X <- .rf_frame(newdata, object$predictor_vars, include_coords,
                 .caller = "predict.rf_fit")
  # A factor level the forest was never grown with has no split to follow;
  # ranger 0.16 returns a plausible-looking number for it rather than erroring.
  # The helper intersects on name, so the whole training frame can be handed
  # over -- the added coordinate columns are numeric and simply skipped.
  X <- .rf_align_levels(X, sf::st_drop_geometry(object$data_sf))

  p <- tryCatch(
    stats::predict(object$engine, data = X, ...)$predictions,
    error = function(e) {
      .log_warn("predict.rf_fit(): ranger predict failed: %s",
                conditionMessage(e))
      NULL
    }
  )
  # predict.all = TRUE and type = "quantiles" make ranger return a matrix.
  # as.numeric() would flatten it column-major into a vector of the wrong
  # length, which .expand_predictions() would then either reject or (at the
  # equal-length fast path) hand back as if it were a prediction per row.
  if (!is.null(dim(p)))
    stop(sprintf(paste0("predict.rf_fit(): ranger returned a %s matrix rather ",
                        "than one value per row. Arguments such as ",
                        "`predict.all = TRUE` or `type = \"quantiles\"` are ",
                        "not supported here, because this method must return ",
                        "one prediction per row of newdata. Call ",
                        "predict(<fit>$engine, data = ...) directly instead."),
                 paste(dim(p), collapse = " x ")), call. = FALSE)
  preds <- if (is.null(p)) rep(NA_real_, n_new) else as.numeric(p)
  .expand_predictions(preds, clean_idx, n_orig)
}


#' Out-of-bag predictions from a random forest fit
#'
#' Returns out-of-bag predictions rather than in-sample ones. See
#' \code{\link{fit_rf_model}} for why, and for what it means for
#' \code{summary()}.
#'
#' @param object An \code{rf_fit}.
#' @param ... Ignored.
#' @return Numeric vector of length \code{object$n}.
#' @export
fitted.rf_fit <- function(object, ...) {
  p <- object$engine$predictions
  if (is.null(p) || length(p) != object$n) {
    .log_warn("fitted.rf_fit(): no usable out-of-bag predictions; returning NA.")
    return(rep(NA_real_, object$n))
  }
  as.numeric(p)
}


#' Out-of-bag residuals from a random forest fit
#'
#' @param object An \code{rf_fit}.
#' @param ... Ignored.
#' @return Numeric vector of length \code{object$n}.
#' @export
residuals.rf_fit <- function(object, ...) {
  y <- sf::st_drop_geometry(object$data_sf)[[object$response_var]]
  as.numeric(y) - fitted(object)
}


#' Coefficients are undefined for a random forest
#'
#' Consistent with \code{coef.gwr_fit()} and \code{coef.bayesian_fit()}, which
#' also error rather than returning \code{NULL} when they cannot supply
#' coefficients -- see the \code{coef()} contract in
#' \code{\link{new_spatial_fit}}.
#'
#' @param object An \code{rf_fit}.
#' @param ... Ignored.
#' @return Never returns; always signals an error.
#' @export
coef.rf_fit <- function(object, ...) {
  stop("coef(): a random forest has no coefficients. For per-predictor ",
       "influence use `$info$importance` (", object$info$importance_type,
       " importance).", call. = FALSE)
}


#' Print a random forest fit
#'
#' @param x An \code{rf_fit}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.rf_fit <- function(x, ...) {
  cat("<Random Forest (ranger)> spatial model fit\n")
  cat(sprintf("  Formula : %s\n", deparse(x$formula)))
  cat(sprintf("  n       : %d\n", x$n))
  cat(sprintf("  Trees   : %d (mtry = %s, min node = %s)\n",
              x$info$num_trees %||% NA_integer_,
              format(x$info$mtry %||% NA), format(x$info$min_node_size %||% NA)))
  cat(sprintf("  Coords as predictors: %s\n",
              if (isTRUE(x$info$include_coords)) "YES - see ?fit_rf_model" else "no"))
  # .num1() rather than %||%: an unset ranger field arrives as numeric(0),
  # which is not NULL, so is.finite() would error and sprintf() would print
  # nothing at all.
  oob_rmse <- .num1(x$info$oob_rmse)
  if (is.finite(oob_rmse))
    cat(sprintf("  OOB RMSE: %.4f   OOB R\u00b2: %.4f\n",
                oob_rmse, .num1(x$info$oob_r_squared)))
  imp <- x$info$importance
  if (!is.null(imp) && length(imp) > 0L) {
    top <- utils::head(sort(imp, decreasing = TRUE), 5L)
    cat(sprintf("  Importance (%s): %s\n", x$info$importance_type,
                paste(sprintf("%s=%.4g", names(top), top), collapse = ", ")))
  }
  cat("\n  OOB is a random hold-out and is optimistic under spatial\n")
  cat("  autocorrelation; use cv_rf() for a spatial estimate.\n")
  invisible(x)
}
