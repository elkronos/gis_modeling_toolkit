# =============================================================================
# Random forest backend (ranger)
# =============================================================================



#' The include_coords caveat, raised once per session
#'
#' Kept here rather than inline because BOTH \code{fit_rf_model()} and
#' \code{cv_rf()} raise it, and they must not drift apart.  \code{cv_rf()}
#' raises it in the parent BEFORE dispatching folds: \code{.log_warn_once()}
#' records the key in a package-level environment, and under
#' \code{parallel = TRUE} that write happens in a forked child and dies with
#' it -- so the caveat printed once per worker, the parent's registry stayed
#' empty, and the next sequential fit printed it again.  A fork COPIES the
#' parent's memory, so setting the flag before the fork makes every child
#' inherit it and stay quiet, and leaves the parent correctly stamped.
#'
#' @return \code{TRUE} if this call logged, \code{FALSE} if it was already
#'   raised in this session.  Invisibly.
#' @keywords internal
#' @noRd
.warn_include_coords <- function() {
  .log_warn_once("rf_include_coords",
                 paste0("fit_rf_model(): include_coords = TRUE. The forest ",
                        "splits on location and can fit the training surface ",
                        "by memorising location, so its error grows quickly ",
                        "with the distance from the data. The out-of-bag ",
                        "error and random cross-validation measure it among ",
                        "the samples: the error of a map over an area the ",
                        "samples cover evenly, which reads low for ",
                        "predictions made farther from the data. How ",
                        "much the coordinates help depends on ",
                        "how far from the training data the predictions are: ",
                        "the gain is largest between sampled points, fades ",
                        "with distance, and can reverse in a region the data ",
                        "never visited (Meyer et al. 2019). Score the forest ",
                        "with cv_rf() on blocked or distance-matched folds ",
                        "at the distances the map will predict at ",
                        "(make_folds(method = \"block_kfold\") or \"nndm\"), ",
                        "and compare it on the same folds with the forest ",
                        "fitted without coordinates: build the folds once ",
                        "with make_folds() and pass them as `folds` to two ",
                        "cv_rf() calls, one with include_coords = TRUE and ",
                        "one without. ?cv_rf shows the two calls."))
}


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


# Bounds of a response the forest is grown on as it is given.  ranger scores
# a split by sum_left^2 / n_left + sum_right^2 / n_right, the square of a sum
# of up to n values: values up to 1e100 square below the largest double on
# any number of rows a forest can be grown on, and a spread down to 1e-100
# squares well above the smallest normal double.  Beyond them the squares
# overflow or lose their bits, and the response is divided by a power of two.
.RF_Y_MAX <- 1e100
.RF_Y_MIN <- 1e-100


#' The response a forest is grown on, taken about its centre and scaled
#'
#' ranger scores a split by a sum of squares it forms from sums of the
#' response and their squares.  At an offset of 1e8 (1e8 plus a spread of
#' 1) those squares are near 1e17 and their differences read the rounding,
#' so the trees split elsewhere than on the response less the offset; at
#' 1e160 the squares overflow and at 1e-160 they lose their bits.  The
#' forest's splits do not change under a shift or a scale of the response,
#' so the response is taken less the centre of \code{.spread_scaled()}
#' where that centre is more than \code{.PRED_CENTRE_RATIO} times its
#' spread (the predictors' gate in \code{.centre_columns()}), and divided
#' by the power of two \code{.spread_scaled()} gives where its spread is
#' below \code{.RF_Y_MIN} or its largest absolute value, after the centre,
#' above \code{.RF_Y_MAX}.  Every other response, a constant one among
#' them, is returned as it is, so its forest is the one it was to the last
#' bit.  Not under \code{splitrule = "beta"}, which needs the response as
#' given in (0, 1).
#'
#' @param y The response.
#' @param splitrule The \code{splitrule} passed to ranger, or \code{NULL}.
#' @return A list: \code{y}, the response the forest is grown on, which is
#'   \code{(y - centre) / scale}; \code{centre} (0 where none is taken) and
#'   \code{scale} (1 where none is applied).
#' @keywords internal
#' @noRd
.rf_response <- function(y, splitrule = NULL) {
  out <- list(y = y, centre = 0, scale = 1)
  if (identical(splitrule, "beta")) return(out)
  f <- y[is.finite(y)]
  if (length(f) < 2L || !(max(abs(f)) > 0)) return(out)
  sc <- .spread_scaled(f)
  if (!(sc$spread > 0)) return(out)
  centre <- if (abs(sc$centre) > .PRED_CENTRE_RATIO * sc$spread) sc$centre else 0
  yc <- if (centre != 0) y - centre else y
  a  <- max(abs(yc[is.finite(yc)]))
  scale <- if (sc$spread < .RF_Y_MIN || a > .RF_Y_MAX) sc$scale else 1
  if (scale != 1) yc <- yc / scale
  list(y = yc, centre = centre, scale = scale)
}


#' Values of the forest's response put back in the response's units
#'
#' The forest in \code{$engine} is grown on \code{(y - centre) / scale}
#' (\code{.rf_response()}); its predictions are \code{centre + scale * p}.
#' A fit with no centre or scale recorded (one saved before they were) is
#' returned as it is.
#'
#' @param p Predictions of the forest.
#' @param info The fit's \code{$info}.
#' @return \code{p} in the response's units.
#' @keywords internal
#' @noRd
.rf_unscale <- function(p, info) {
  centre <- info$response_centre %||% 0
  scale  <- info$response_scale %||% 1
  if (scale != 1) p <- p * scale
  if (centre != 0) p <- p + centre
  p
}


#' Check that categorical predictor levels in newdata were seen at fit time
#'
#' \code{ranger} matches factor predictors by level, so a level that was not
#' present when the forest was grown has no split to follow.  Depending on the
#' ranger version this either errors or, as in ranger 0.16, silently returns a
#' plausible-looking number.  An error from ranger gives its own reason, which
#' does not name the column.  Both are worse than an error naming the level.
#'
#' @param X Prediction frame from \code{.rf_frame()}.
#' @param train_X Training frame from \code{.rf_frame()}.
#' @return \code{X}, with categorical columns coerced to the training levels.
#' @keywords internal
#' @noRd
.rf_align_levels <- function(X, train_X) {
  for (cn in intersect(names(X), names(train_X))) {
    tr <- train_X[[cn]]
    # A column the forest was grown on as NUMERIC must arrive numeric.  ranger
    # does not check: it factor-codes a character column and applies the
    # numeric split thresholds to the codes, and predicts confidently from the
    # result.  Measured: a numeric-at-fit predictor supplied as text at predict
    # time gave cor(correct, returned) = 0.41 and R2 -1.91 against 0.98, with
    # no condition raised.  predict.bayesian_fit() refuses this; so does this.
    # is.logical(tr) too: a LOGICAL-at-fit predictor fell through all three
    # branches (is.numeric / is.factor / is.character are each FALSE) and was
    # skipped entirely, so a "TRUE"/"FALSE" character column -- exactly what a
    # CSV round trip produces -- was factor-coded 1/2 against splits built on
    # 0/1, sending every row to the TRUE side.  Measured: cor(correct,
    # returned) = 0.21 and FALSE rows predicted +1.5 where the truth was -1.75,
    # with nothing raised.  ranger sees a logical as 0/1, so treat it as
    # numeric here.
    if ((is.numeric(tr) || is.logical(tr)) &&
        !is.numeric(X[[cn]]) && !is.logical(X[[cn]]))
      stop(sprintf(paste0("predict.rf_fit(): predictor '%s' was %s when ",
                          "the forest was grown but is %s in `newdata`. ranger ",
                          "would silently code it as a factor and apply numeric ",
                          "split thresholds to the codes; convert the column ",
                          "first."),
                   cn, if (is.logical(tr)) "logical" else "numeric",
                   class(X[[cn]])[1L]), call. = FALSE)
    if (!(is.factor(tr) || is.character(tr))) next
    # Levels that OCCUR in the training rows, not the level SET.  An ordinary
    # subset (dat[dat$g != "c", ], or a spatial-CV training fold in which a
    # class sits entirely in the held-out block) keeps the unused level, so a
    # level with zero training rows passed the guard and ranger returned some
    # other level's prediction for it (49.9 where the truth was ~100).  The
    # documentation promises "an unseen level is an error, not a guess", and an
    # unseen level is one the forest never saw a row of.
    lv  <- sort(unique(as.character(tr[!is.na(tr)])))
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
    X[[cn]] <- if (is.factor(tr)) factor(val, levels = levels(tr)) else val
  }
  X
}


# -----------------------------------------------------------------------------
# Exported
# -----------------------------------------------------------------------------

#' Fit a random forest with ranger
#'
#' Fits a regression random forest on an \code{sf} layer and returns it as a
#' \code{spatial_fit}, so it works with \code{\link{cv_spatial}()},
#' \code{\link{predict_surface}()}, \code{\link{area_of_applicability}()}
#' and the \code{plot()} method like any other backend.  A forest suits a
#' response that is non-linear in its predictors when predictive accuracy
#' matters more than an interpretable model.
#'
#' @section Coordinates as predictors:
#' \code{include_coords = TRUE} adds the x and y coordinates to the
#' predictors.  They are the only spatial information a forest can be given:
#' without them it sees the covariates alone, and \code{print()} says "no
#' spatial information used".  The default is \code{FALSE}.
#'
#' The gain from the coordinates depends on where the predictions are
#' wanted.  A forest splits on location as on any column, so it follows the
#' training surface closely between sampled points, and has only the nearest
#' edge of the data to go on beyond them.  Between sampled points the
#' coordinates can lower the error, most where the covariates carry little
#' of the spatial structure themselves, and the gain fades with distance from
#' the training data.  Beyond the autocorrelation range, or in an area the
#' data never visited, the coordinates can make predictions worse (Meyer et
#' al. 2019).  With covariates that carry the spatial structure themselves,
#' expect less from the coordinates.
#'
#' A forest with coordinates fits location closely, so its error grows
#' quickly with the distance from the data, and a cross-validated score is
#' the error at the distance its folds put between the held-out points and
#' their training data.  Random folds, and the out-of-bag error, a random
#' hold-out too, measure it among the samples: the error of a map over an
#' area the samples cover evenly, which reads low for predictions made
#' farther from the data.  Score the forest with \code{\link{cv_rf}()} on
#' folds at the distances the map will be used at (folds matched to the
#' prediction distance, \code{make_folds(method = "nndm")} or
#' \code{"knndm"}, or blocks, \code{method = "block_kfold"}, for predictions
#' into unsampled gaps), and compare the forests with and without
#' coordinates on the same folds.  Build the folds once and pass them to
#' \code{cv_rf()} with and without the coordinates:
#' \preformatted{
#' folds   <- make_folds(dat, k = 5, method = "block_kfold")
#' plain   <- cv_rf(dat, "z", c("a", "b"), folds = folds)
#' with_xy <- cv_rf(dat, "z", c("a", "b"), folds = folds,
#'                  include_coords = TRUE)
#' c(plain = plain$overall$RMSE, with_xy = with_xy$overall$RMSE)
#' }
#' The example on \code{\link{cv_rf}()} runs this comparison and reads it
#' fold by fold.  Setting \code{include_coords = TRUE} logs this advice once
#' per session, as a log line rather than an R warning, and not once for
#' every fold of a cross-validation.
#'
#' @section Categorical predictors:
#' A predictor that is an unordered factor or a character column is split
#' with \code{respect.unordered.factors = "order"}: ranger orders its levels
#' by mean response once, before growing, and splits on that order, the
#' handling ranger's own help recommends.  ranger's default for a regression
#' forest, \code{"ignore"}, splits on the order the levels are stored in
#' (alphabetical, for most data), which for a nominal variable is arbitrary.
#' Pass \code{respect.unordered.factors} through \code{...} to choose
#' otherwise (\code{"ignore"}, or \code{"partition"} for every two-way split
#' of the levels, which is slow beyond a few levels).  Ordered factors are
#' left alone and split on their own order, and with
#' \code{splitrule = "extratrees"} ranger's own default for that rule,
#' \code{"partition"}, is kept.  The setting is recorded in
#' \code{$info$respect_unordered_factors} and printed with the fit.
#'
#' @section Out-of-bag error:
#' \pkg{ranger}'s out-of-bag (OOB) error holds each observation out of the
#' trees that did not sample it.  That is a random hold-out: the trees that
#' did not see a point almost always saw its neighbours, so under spatial
#' autocorrelation it is the error of predicting among the samples, as
#' random k-fold is, and reads low for predictions made farther from the
#' data than the samples sit from each other (see \code{\link{make_folds}()}).
#' It is reported as \code{$info$oob_rmse} and \code{$info$oob_r_squared}
#' and labelled as OOB wherever it appears.  Use \code{\link{cv_rf}()} for a
#' spatial estimate.
#'
#' @section Fitted values are out-of-bag:
#' \code{fitted()} on an \code{rf_fit} returns out-of-bag predictions, as the
#' random forest packages themselves do.  In-sample predictions from a forest
#' come close to reproducing the training data, and their \eqn{R^2} says
#' little about new data.  \code{summary()} therefore reports out-of-bag
#' metrics for a forest and in-sample metrics for a GWR or a Bayesian GP
#' fit, and the two
#' do not compare; \code{\link{compare_models_cv}()} scores every backend on
#' the same folds.
#'
#' A row that every tree sampled has no out-of-bag prediction, and ranger
#' reports \code{NaN} for it.  That is every row under \code{replace = FALSE}
#' with \code{sample_fraction = 1}, and a few under a small
#' \code{num_trees}.  The fit warns with the count; \code{fitted()} and
#' \code{residuals()} are \code{NaN} on those rows, and \code{summary()}
#' says how many rows its metrics were computed on.  With no row out of bag
#' at all, the OOB error (\code{NA} in \code{$info$oob_rmse} and
#' \code{$info$oob_r_squared}) and the permutation importance (\code{NaN})
#' are undefined too, and \code{print()} says so.  \code{\link{cv_rf}()}
#' scores its fold forests on the held-out rows, never out of bag, so it
#' warns once with the number of folds affected rather than once per fold.
#'
#' @param data_sf An sf object with response, predictors and geometry.
#' @param response_var Response column name.
#' @param predictor_vars Predictor column names.
#' @param num_trees Number of trees. Default 500.
#' @param mtry Predictors sampled per split. \code{NULL} uses ranger's default.
#' @param min_node_size Minimum node size. \code{NULL} uses ranger's default.
#' @param importance \code{"permutation"} (default), \code{"impurity"} or
#'   \code{"none"}.  Impurity importance is biased towards continuous and
#'   high-cardinality predictors (Strobl et al. 2007), so permutation
#'   importance is the default although it costs more.
#' @param include_coords Add the coordinates as predictors. Default
#'   \code{FALSE}, which leaves the forest with no spatial information; see
#'   \emph{Coordinates as predictors}.
#' @param replace Logical; grow each tree on a bootstrap sample drawn
#'   \emph{with} replacement (\code{TRUE}, ranger's default and this one) or
#'   on a subsample drawn without replacement (\code{FALSE}).  Strobl et al.
#'   (2007) show that bootstrap sampling with replacement is itself a source
#'   of bias in variable importance toward predictors with many distinct
#'   values or categories, and recommend subsampling without replacement.
#'   The forest's predictions change slightly with the choice.  The setting
#'   is recorded in \code{$info} and printed with the fit.
#' @param sample_fraction Fraction of rows drawn for each tree.  \code{NULL}
#'   (default) uses ranger's rule: all rows when \code{replace = TRUE},
#'   0.632 (the expected share of distinct rows in a bootstrap sample)
#'   when \code{replace = FALSE}.  A single number in (0, 1] overrides it.
#'   \code{replace = FALSE} with \code{sample_fraction = 1} grows every tree
#'   on every row, so nothing is out of bag: the fit warns (see "Fitted
#'   values are out-of-bag").
#' @param seed Seed passed to ranger. Default 123.
#' @param num_threads Threads for ranger.  Default \code{NULL} means
#'   \code{getOption("mc.cores", 1L)}: one thread unless the session has
#'   set more, the convention \pkg{brms} and \pkg{parallel} use.
#'   Predictions do not depend on the thread count; pass
#'   \code{parallel::detectCores()} to use every core.
#' @param .already_prepped Internal; leave at the default.
#' @param ... Passed to \code{ranger::ranger()}.  ranger's own spellings of
#'   the arguments this function already sets (\code{num.trees},
#'   \code{min.node.size}, \code{num.threads}, \code{mtry},
#'   \code{importance}, \code{seed}, \code{replace},
#'   \code{sample.fraction}, \code{x}, \code{y}) are rejected with a
#'   message naming the wrapper argument to use instead.
#'   \code{respect.unordered.factors} passed here overrides the
#'   \code{"order"} this function sets when a predictor is categorical; see
#'   \emph{Categorical predictors}.
#'
#' @return An \code{rf_fit} object (inherits from \code{spatial_fit}).
#'   \code{$info} carries \code{num_trees}, \code{mtry}, \code{min_node_size},
#'   \code{importance_type}, \code{importance} (a named numeric, or
#'   \code{NULL} when \code{importance = "none"}), \code{include_coords},
#'   \code{categorical_predictors} (the unordered factor and character
#'   predictors, possibly none) and \code{respect_unordered_factors} (how
#'   ranger was told to split them: \code{"order"} unless overridden through
#'   \code{...}, and \code{NA} when there is none and nothing was passed),
#'   \code{replace} and \code{sample_fraction} (the sampling each tree was
#'   grown on, with \code{sample_fraction} resolved to the number ranger
#'   used), \code{oob_rmse} and \code{oob_r_squared} (each \code{NA_real_} when
#'   ranger did not compute it; forwarding \code{oob.error = FALSE} through
#'   \code{...} is one way to get there), \code{fitted_are_oob} (always
#'   \code{TRUE}; \code{summary()} reads it to label its metrics),
#'   \code{seed} and \code{n_dropped} (the rows \code{prep_model_data()}
#'   removed for missing or non-finite values or a bad geometry, so \code{$n}
#'   can be read against \code{nrow(data_sf)}), and \code{response_centre}
#'   and \code{response_scale} (below). The raw forest is in
#'   \code{$engine}.
#'
#'   The forest is grown on the response taken about its centre where that
#'   centre is more than 1000 times the response's spread (as at 1e8 plus
#'   a spread of 1), and divided by a power of two where the response is so
#'   large or its spread so small that its squares overflow or lose their
#'   bits (beyond about 1e100 or 1e-100).  The splits of a forest do not
#'   change under a shift or a scale of the response, and ranger's split
#'   scores, formed from squares of the response, read the rounding of the
#'   offset otherwise.  \code{$info$response_centre} (0 where none is taken)
#'   and \code{$info$response_scale} (1 where none is applied) record them,
#'   and they are put back on \code{fitted()}, \code{predict()},
#'   \code{residuals()}, \code{oob_rmse} and \code{importance} (squared
#'   units, so multiplied by the scale twice; beyond the range of doubles it
#'   is \code{Inf} or loses its bits, and a WARN line says so).  The
#'   permutation importance under \code{scale.permutation.importance = TRUE}
#'   and the impurity importance under \code{splitrule = "maxstat"} have no
#'   units and are left as ranger reports them.  The forest
#'   in \code{$engine} predicts \code{(y - response_centre) /
#'   response_scale}, and its own \code{predictions} and
#'   \code{prediction.error} are in those units.  Every other response is
#'   passed to ranger as it is.
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
#' @seealso \code{\link{cv_rf}()} to score the forest on spatial folds;
#'   \code{predict()} (\link{spatial_fit-predict}) and
#'   \code{\link{predict_surface}()} to predict;
#'   \code{\link{area_of_applicability}()} to mark where the predictions
#'   apply (it weights the predictors by this forest's importance).
#' @family model fitting
#' @examplesIf requireNamespace("ranger", quietly = TRUE)
#' library(sf)
#' set.seed(1)
#' n <- 150
#' dat <- st_as_sf(
#'   data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
#'              a = rnorm(n), b = rnorm(n)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' dat$z <- 2 * dat$a - dat$b + rnorm(n, 0, 0.3)
#' fit <- fit_rf_model(dat, "z", c("a", "b"))
#' fit
#' fit$info$importance     # permutation importance of a and b
#' @export
fit_rf_model <- function(data_sf, response_var, predictor_vars,
                         num_trees = 500L, mtry = NULL, min_node_size = NULL,
                         importance = c("permutation", "impurity", "none"),
                         include_coords = FALSE, replace = TRUE,
                         sample_fraction = NULL, seed = 123L,
                         num_threads = NULL, .already_prepped = FALSE, ...) {
  if (!inherits(data_sf, "sf"))
    stop("fit_rf_model(): `data_sf` must be an sf object.", call. = FALSE)
  if (!is.logical(replace) || length(replace) != 1L || is.na(replace))
    stop("fit_rf_model(): `replace` must be TRUE or FALSE.", call. = FALSE)
  # ranger's own default for sample.fraction is ifelse(replace, 1, 0.632);
  # resolve it here so the number that reached the forest is recorded.
  if (is.null(sample_fraction)) {
    sample_fraction <- if (isTRUE(replace)) 1 else 0.632
  } else if (!is.numeric(sample_fraction) || length(sample_fraction) != 1L ||
             !is.finite(sample_fraction) || sample_fraction <= 0 ||
             sample_fraction > 1) {
    stop("fit_rf_model(): `sample_fraction` must be a single number in (0, 1], ",
         "or NULL for ranger's default.", call. = FALSE)
  }
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
  # Rows prep_model_data() removed on the way in, read from the record it
  # leaves -- which a layer prepared by the caller may still carry.  0 when
  # there is none, and when the record no longer describes this layer: a
  # fold's subset inside cv_*() has neither, and cv_*() reports its own.
  n_dropped <- as.integer(.get_row_record(dat, "dropped")$n %||% 0L)

  if (isTRUE(include_coords)) .warn_include_coords()

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

  # `...` is documented as "passed to ranger::ranger()", and ranger's own
  # spellings of the arguments this wrapper already fixes reach it TWICE:
  # num.trees, min.node.size, num.threads, x, y, importance, seed.  R then
  # refuses the call with "formal argument ... matched by multiple actual
  # arguments" -- and through cv_rf() or compare_models_cv(rf_args =
  # list(num.trees = 50)) that means every fold fails and the result is an
  # all-NA row.  Name the wrapper argument the caller should have used.
  dot_names <- names(list(...))
  ranger_dupes <- c(num.trees = "num_trees", min.node.size = "min_node_size",
                    num.threads = "num_threads", mtry = "mtry",
                    importance = "importance", seed = "seed",
                    replace = "replace", sample.fraction = "sample_fraction",
                    x = NA_character_, y = NA_character_)
  hit_dupes <- intersect(dot_names, names(ranger_dupes))
  if (length(hit_dupes)) {
    repl <- ranger_dupes[hit_dupes]
    stop(sprintf(paste0("fit_rf_model(): %s already set by this function and ",
                        "cannot be passed through `...` as well (ranger would ",
                        "receive it twice). Use %s instead."),
                 paste(sprintf("`%s` is", hit_dupes), collapse = ", "),
                 paste(ifelse(is.na(repl), "the x/y interface",
                              sprintf("`%s`", repl)), collapse = ", ")),
         call. = FALSE)
  }

  # NULL reaching ranger means "every core" -- ranger resolves it to
  # hardware_concurrency().  Resolve it here instead, to the session's
  # mc.cores opt-in or one thread.
  num_threads <- .sanitize_core_count(
    if (is.null(num_threads)) getOption("mc.cores", 1L) else num_threads)

  # Unordered factors and character columns.  ranger's default for a
  # regression forest, respect.unordered.factors = "ignore", splits them on
  # the order their levels are stored in, which for a nominal variable is
  # usually the alphabet: a 40-level factor with random level effects gave a
  # test RMSE of 1.05 that way against 0.57 with "order" (levels ordered by
  # mean response once, before growing), lower in 10 draws of 10.  ranger's
  # help recommends "order", so that is what a categorical predictor gets
  # unless the caller says otherwise through `...`.  Ordered factors carry an
  # order of their own and are untouched by the setting.  The one case left
  # to ranger is splitrule = "extratrees", where its own default is already
  # "partition" and not "ignore".
  cat_pred <- names(X)[vapply(X, function(v)
    (is.factor(v) && !is.ordered(v)) || is.character(v), logical(1))]
  ruf_set    <- "respect.unordered.factors" %in% dot_names
  extratrees <- "splitrule" %in% dot_names &&
    identical(list(...)[["splitrule"]], "extratrees")
  ruf_order  <- !ruf_set && !extratrees && length(cat_pred) > 0L
  ruf <- if (ruf_set) {
    v <- list(...)[["respect.unordered.factors"]]
    if (is.null(v)) NA_character_
    else if (is.logical(v)) (if (isTRUE(v)) "order" else "ignore")
    else as.character(v)[1L]
  } else if (length(cat_pred) == 0L) NA_character_
  else if (extratrees) "partition" else "order"

  # The forest is grown on the response taken about its centre where it
  # lies far from zero, and divided by a power of two where its squares
  # would overflow or lose their bits; see .rf_response().  The centre and
  # scale are put back on everything the fit returns.  `y` keeps its name,
  # which ranger records in its call.
  splitrule <- if ("splitrule" %in% dot_names) list(...)[["splitrule"]]
  ry <- .rf_response(y, splitrule)
  y <- ry$y

  # .call_capturing_stderr(): ranger diagnoses a bad argument on stderr from
  # C++ and then throws "User interrupt or internal error.", so the message
  # this reported named nothing a caller could act on while the real reason
  # printed itself to the console past every handler.
  # Two spelled-out calls rather than do.call(): ranger keeps its call, and
  # do.call() would put the whole predictor frame inside it.
  rr <- .call_capturing_stderr(function()
    if (ruf_order)
      ranger::ranger(x = X, y = y, num.trees = as.integer(num_trees),
                     mtry = mtry, min.node.size = min_node_size,
                     importance = importance, seed = seed,
                     replace = replace, sample.fraction = sample_fraction,
                     num.threads = num_threads,
                     respect.unordered.factors = "order", ...)
    else
      ranger::ranger(x = X, y = y, num.trees = as.integer(num_trees),
                     mtry = mtry, min.node.size = min_node_size,
                     importance = importance, seed = seed,
                     replace = replace, sample.fraction = sample_fraction,
                     num.threads = num_threads, ...))
  if (!is.null(rr$error))
    stop(sprintf("fit_rf_model(): ranger() failed: %s",
                 .stderr_reason(rr$error, rr$stderr)), call. = FALSE)
  fit <- rr$value

  imp <- if (identical(importance, "none")) NULL else fit$variable.importance
  # .num1() collapses an unset ranger field to NA_real_ rather than the
  # numeric(0) that as.numeric(NULL) produces -- see ?.num1.
  oob_mse <- .num1(fit$prediction.error)
  oob_rmse <- if (is.finite(oob_mse)) sqrt(oob_mse) else NA_real_
  # The importance is in squared units of the response, and so is the OOB
  # error: they are multiplied by the scale twice (scale^2 itself can pass
  # the largest double), and the OOB RMSE by the scale once.  A shift
  # changes neither.  Two importances ranger reports without units are left
  # as they are: permutation importance under scale.permutation.importance
  # = TRUE (divided by its standard error) and impurity importance under
  # splitrule = "maxstat" (a rank statistic).
  imp_unitless <-
    (identical(importance, "permutation") &&
       isTRUE(list(...)[["scale.permutation.importance"]])) ||
    (identical(importance, "impurity") &&
       identical(splitrule, "maxstat"))
  if (ry$scale != 1) {
    oob_rmse <- oob_rmse * ry$scale
    if (!is.null(imp) && !imp_unitless) {
      imp_fit <- imp
      imp <- imp * ry$scale * ry$scale
      lost <- is.finite(imp_fit) & imp_fit != 0 &
        (!is.finite(imp) | abs(imp) < .Machine$double.xmin)
      # Not for a cv_rf() fold forest, whose importance the run never reads
      # (as for the out-of-bag warning below).
      if (any(lost) && !isTRUE(.already_prepped))
        .log_warn(paste0("fit_rf_model(): the %s importance of %s is in ",
                         "squared units of the response and lies %s, so it ",
                         "is reported as %s. Rescale the response to compare ",
                         "the importances."),
                  importance, paste(sQuote(names(imp)[which(lost)]), collapse = ", "),
                  if (ry$scale > 1) "above the largest double"
                  else "below the smallest normal double",
                  if (ry$scale > 1) "Inf"
                  else "a number with few or none of its bits")
    }
  }

  # A row every tree sampled has no out-of-bag prediction, and ranger returns
  # NaN for it: every row under replace = FALSE with sample_fraction = 1
  # (each tree is grown on all of them), a few under a small num_trees.  Then
  # fitted() and residuals() are NaN there, summary() and model_metrics() score
  # the remaining rows while summary() heads its output "n = <all of them>",
  # and with no row out of bag at all the OOB error and the permutation
  # importance are NaN as well -- all with nothing said.  The forest itself is
  # sound and predict(newdata =) and cv_rf() are unaffected, so warn rather
  # than refuse.
  #
  # Not for a cv_rf() fold forest (.already_prepped = TRUE): cv_rf() scores
  # predict() on the held-out rows and never reads a fold forest's out-of-bag
  # predictions, so this fired once per fold about nothing the run used --
  # and told the user to score the forest with cv_rf() from inside cv_rf().
  # cv_rf() counts the folds itself and says it once (.rf_cv_oob_note()).
  oob_pred <- fit$predictions
  n_no_oob <- if (is.numeric(oob_pred) && length(oob_pred) == nrow(X))
    sum(!is.finite(oob_pred)) else 0L
  fold_fit <- isTRUE(.already_prepped)
  if (!fold_fit && n_no_oob == nrow(X)) {
    .warn_and_log(paste0("fit_rf_model(): no row is out of bag for any tree ",
                         "(%s), so fitted(), residuals(), the out-of-bag error%s ",
                         "are undefined (NaN) and summary() and model_metrics() ",
                         "have nothing to score.%s Use replace = TRUE or a ",
                         "sample_fraction below 1, or score the forest with ",
                         "cv_rf()."),
                  if (!isTRUE(replace) && sample_fraction >= 1)
                    "replace = FALSE with sample_fraction = 1 grows every tree on every row"
                  else sprintf("%d tree(s)", as.integer(num_trees)),
                  if (identical(importance, "permutation"))
                    " and the permutation importance" else "",
                  # pmax(NaN, 0) is NaN, so the weights ?fit_rf_model
                  # suggests for area_of_applicability() are refused too.
                  if (identical(importance, "permutation"))
                    paste0(" area_of_applicability() cannot be weighted by ",
                           "that importance either: pass weights = NULL, or ",
                           "leave its default weights = \"auto\", which uses ",
                           "equal weights for this fit.")
                  else "")
  } else if (!fold_fit && n_no_oob > 0L) {
    .warn_and_log(paste0("fit_rf_model(): %d of %d rows were sampled by every ",
                         "one of the %d tree(s) and so have no out-of-bag ",
                         "prediction: fitted() and residuals() are NaN there, ",
                         "and summary(), model_metrics() and the out-of-bag ",
                         "error use the other %d. Raise num_trees to cover ",
                         "every row."),
                  n_no_oob, nrow(X), as.integer(num_trees), nrow(X) - n_no_oob)
  }

  new_spatial_fit(
    subclass       = "rf_fit",
    engine         = fit,
    # Show the coordinates in the formula when they are predictors, so
    # print()ing the fit does not hide them.  It is display-only: the forest is
    # built through ranger's x/y interface, never from this formula.  Hence
    # env = globalenv(): reformulate()'s default is this frame, which holds the
    # forest (`fit`, `rr`), the data and the predictor frame, and a formula
    # serialises its environment -- so saveRDS() on an rf_fit wrote the forest
    # out a second time (1.62 MB for a 100-tree forest of 0.72 MB).
    formula        = stats::reformulate(
      termlabels = if (isTRUE(include_coords))
        c(predictor_vars, "..x", "..y") else predictor_vars,
      response = response_var, env = globalenv()),
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
      # The unordered factor and character predictors, and how ranger was
      # told to split them (NA: none, and nothing passed).
      categorical_predictors    = cat_pred,
      respect_unordered_factors = ruf,
      replace          = isTRUE(replace),
      sample_fraction  = sample_fraction,
      oob_rmse         = oob_rmse,
      oob_r_squared    = .num1(fit$r.squared),
      fitted_are_oob   = TRUE,
      # The forest in `engine` is grown on (y - response_centre) /
      # response_scale; see .rf_response().
      response_centre  = ry$centre,
      response_scale   = ry$scale,
      seed             = seed,
      n_dropped        = n_dropped
    )
  )
}


#' Out-of-bag coverage of the cv_rf() fold forests, reported once per run
#'
#' \code{fit_rf_model()} does not warn about rows no tree left out of bag
#' when it grows a fold forest (\code{.already_prepped = TRUE}): cv_rf()
#' scores \code{predict()} on the held-out rows and never reads a fold
#' forest's out-of-bag predictions.  It warned once per fold, and told the
#' user to score the forest with cv_rf() from inside cv_rf().  What the folds
#' share is still worth one warning, because a forest grown on all the data
#' with the same settings has the same gap.
#'
#' The counts travel back through \code{cv_spatial()}'s \code{fold_info_fn}
#' as two \code{fold_metrics} columns, which \code{.rf_cv_oob_note()} takes
#' out again.  Warnings cannot carry them: one raised in a forked worker
#' reaches the parent only once per distinct text, so under
#' \code{parallel = TRUE} the folds could not be counted.
#'
#' @param fit_obj A fold's \code{rf_fit}.
#' @param res The \code{cv_spatial()} result.
#' @param dots The arguments cv_rf() passes to \code{fit_rf_model()}.
#' @return \code{.rf_fold_oob_info()}: a list of two counts.
#'   \code{.rf_cv_oob_note()}: \code{res} without the two columns.
#' @keywords internal
#' @noRd
.rf_fold_oob_info <- function(fit_obj, ...) {
  p <- fit_obj$engine$predictions
  n <- if (is.numeric(p) && is.null(dim(p))) length(p) else 0L
  list(..rf_n_no_oob = if (n > 0L) as.integer(sum(!is.finite(p))) else 0L,
       ..rf_n_fit = as.integer(n))
}

# The run's half of .rf_fold_oob_info(): see there.
.rf_cv_oob_note <- function(res, dots) {
  fm <- res$fold_metrics
  if (!is.data.frame(fm) ||
      !all(c("..rf_n_no_oob", "..rf_n_fit") %in% names(fm)))
    return(res)
  none <- fm$..rf_n_no_oob
  n    <- fm$..rf_n_fit
  fm$..rf_n_no_oob <- NULL
  fm$..rf_n_fit    <- NULL
  res$fold_metrics <- fm
  ok     <- is.finite(none) & is.finite(n) & n > 0
  n_all  <- sum(ok & none == n)
  n_some <- sum(ok & none > 0 & none < n)
  if (n_all > 0L) {
    # The one cause fit_rf_model() can name; anything else (an `inbag`
    # through `...`, say) is the caller's own doing.
    sf    <- dots[["sample_fraction"]]
    whole <- isFALSE(dots[["replace"]]) && is.numeric(sf) && isTRUE(sf >= 1)
    imp   <- dots[["importance"]]
    .warn_and_log(paste0("cv_rf(): no training row was out of bag for any ",
                         "tree in %d of %d fold forest(s)%s. The ",
                         "cross-validation is unaffected, because it scores ",
                         "each forest on its held-out fold; but a forest grown ",
                         "on all the data with the same settings has no ",
                         "%s.%s"),
                  n_all, nrow(fm),
                  if (whole) paste0(" (replace = FALSE with sample_fraction = ",
                                    "1 grows every tree on every row)") else "",
                  if (is.null(imp) || identical(imp, "permutation"))
                    paste0("out-of-bag error, fitted() values or permutation ",
                           "importance")
                  else "out-of-bag error or fitted() values",
                  if (whole) paste0(" Use replace = TRUE or a sample_fraction ",
                                    "below 1 for those.") else "")
  } else if (n_some > 0L) {
    .warn_and_log(paste0("cv_rf(): in %d of %d fold forest(s) some training ",
                         "rows were sampled by every tree and have no ",
                         "out-of-bag prediction. The cross-validation is ",
                         "unaffected, because it scores each forest on its ",
                         "held-out fold; but a forest grown on all the data ",
                         "with the same num_trees may leave rows without one ",
                         "too. Raise num_trees if you will read its out-of-bag ",
                         "error or fitted() values."),
                  n_some, nrow(fm))
  }
  res
}


#' Cross-validate a random forest with spatial folds
#'
#' A wrapper over \code{\link{cv_spatial}()} that refits a \pkg{ranger}
#' forest on each training fold.  By default it holds out whole spatial
#' blocks, so most neighbours of a held-out point are not in the training
#' set, where the out-of-bag error holds out single rows; \code{folds} says
#' what the score estimates.
#'
#' @section Comparing the forest with and without coordinates:
#' Whether \code{include_coords = TRUE} lowers the error depends on how far
#' from the training data the predictions are (see \emph{Coordinates as
#' predictors} in \code{\link{fit_rf_model}}), so it is a comparison to make
#' on your own data, at the block size the map will be used at. Score the
#' two forests on the same folds: build the folds once with
#' \code{\link{make_folds}} and pass the result as \code{folds} to two
#' calls, one with \code{include_coords = TRUE} and one without. Each fold
#' then holds out the same rows in both, so the two \code{overall} rows
#' compare like with like and the two \code{fold_metrics} tables pair up
#' fold by fold.  A difference in \code{overall} that is small beside the
#' spread of the per-fold differences is not evidence for either forest.
#' The example below runs the comparison.  \code{\link{compare_models_cv}()}
#' takes one set of forest arguments (\code{rf_args}) per call, so it scores
#' only one of the two forests.
#'
#' @param data_sf An sf object.
#' @param response_var Response column name.
#' @param predictor_vars Predictor column names.
#' @param folds Optional fold definitions, in any of three shapes: a
#'   \code{\link{make_folds}()} return value; a bare list of
#'   \code{list(train =, test =)} pairs of \code{..row_id} values; or a vector
#'   of fold labels, one per row, which becomes leave-that-label-out splits.
#'   The label vector is how folds built by another package are used here,
#'   since \code{blockCV::cv_spatial()}'s \code{$folds_list} holds two
#'   \emph{unnamed} vectors per fold and is refused by name.  That function
#'   returns a label vector as \code{$folds_ids}.  Train and test must be
#'   disjoint: a fold that trains on its own test rows is not a
#'   cross-validation split and is refused with an error.  IDs naming no row in
#'   the prepared data are dropped with a logged count.  Built via
#'   \code{block_kfold} when \code{NULL}, as in \code{\link{cv_spatial}()}.
#'   With \code{folds = NULL} the score estimates the error at the distance
#'   the block grid built here puts between the held-out points and their
#'   training data, which overstates the error of a map where the samples
#'   cover the area it predicts over.  The excess is how much the model's
#'   error grows from the map's distances to the data to the block grid's,
#'   and is small for a model whose error changes little with the distance
#'   from the data.  \code{\link{fold_separation}()}
#'   with \code{prediction_points} compares the held-out distances with those
#'   of the locations the model will predict at, and
#'   \code{make_folds(method = "knndm", prediction_points = )} builds folds
#'   matched to them.
#' @param k Number of folds when \code{folds} is \code{NULL}. Default 5.
#' @param seed RNG seed.  Default 123.  It seeds fold construction and,
#'   through a per-fold draw, each fold's forest, so two different seeds give
#'   different results even on identical \code{folds}. To grow every fold's
#'   forest from one fixed ranger seed instead, call \code{\link{cv_spatial}}
#'   with your own \code{fit_fn} wrapping \code{fit_rf_model(seed = )}.
#' @param parallel Passed to \code{\link{cv_spatial}}. Default \code{FALSE}.
#'   Under forked workers each fold's forest runs single-threaded unless
#'   \code{num_threads} is passed explicitly, so \code{parallel = 4} means
#'   four threads in total, where it would otherwise mean four times the
#'   session's \code{mc.cores}.
#' @param block_size,auto_range,boundary Passed to \code{\link{cv_spatial}}.
#'   \code{block_size} and \code{auto_range = TRUE} size the blocks built
#'   when \code{folds} is \code{NULL}; beside supplied \code{folds} they
#'   have nothing to size and raise a warning of class
#'   \code{spatialkit_fold_args_ignored}.
#' @param pointize How non-POINT geometry is reduced to a point before
#'   fitting; passed to \code{\link{cv_spatial}}. Default \code{"auto"}.
#' @param metrics Optional scoring function of your own, passed to
#'   \code{\link{cv_spatial}}; see \strong{Your own metrics} there.
#' @param ... Passed to \code{\link{fit_rf_model}} on every fold:
#'   \code{num_trees}, \code{mtry}, \code{importance}, \code{include_coords}
#'   and so on.  \code{data_sf}, \code{response_var}, \code{predictor_vars}
#'   and \code{.already_prepped} are set by this function and must not be
#'   passed here (every fold would fail with "matched by multiple actual
#'   arguments").  \code{seed} is this function's own argument and never
#'   reaches \code{fit_rf_model()} through here; see \code{seed} above for
#'   growing every fold's forest from one fixed seed.
#' @return The \code{\link{cv_spatial}()} result.  \code{MAPE} and
#'   \code{SMAPE} skip rows whose denominator is zero; \code{n_MAPE} and
#'   \code{n_SMAPE} count the rows they used (see "Percentage errors on
#'   responses with zeros" in \code{\link{model_metrics}()}).
#' @seealso \code{\link{fit_rf_model}()} for the forest's arguments, and
#'   \code{\link{plot_cv_metrics}()} to draw the per-fold scores.
#' @family cross-validation
#' @examplesIf requireNamespace("ranger", quietly = TRUE)
#' library(sf)
#' set.seed(1)
#' n <- 200
#' dat <- st_as_sf(
#'   data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
#'              a = rnorm(n)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' # An offset keeps the response away from zero, where MAPE is undefined.
#' dat$z <- 10 + 2 * dat$a + rnorm(n, 0, 0.3)
#' cv_rf(dat, "z", "a", k = 4)$overall
#'
#' # The forest with and without the coordinates, on the same blocked
#' # folds.  The second call logs the include_coords caution once per
#' # session.  This response has no spatial structure, so here the
#' # coordinates can only cost accuracy.
#' folds   <- make_folds(dat, k = 4, method = "block_kfold", seed = 1)
#' plain   <- cv_rf(dat, "z", "a", folds = folds, num_trees = 200)
#' with_xy <- cv_rf(dat, "z", "a", folds = folds, num_trees = 200,
#'                  include_coords = TRUE)
#' c(plain = plain$overall$RMSE, with_xy = with_xy$overall$RMSE)
#' # Fold by fold: positive where the coordinates lowered the error.
#' plain$fold_metrics$RMSE - with_xy$fold_metrics$RMSE
#' @export
cv_rf <- function(data_sf, response_var, predictor_vars, folds = NULL, k = 5,
                  seed = 123, parallel = FALSE, block_size = NULL,
                  auto_range = FALSE, boundary = NULL, pointize = "auto",
                  metrics = NULL, ...) {
  if (!requireNamespace("ranger", quietly = TRUE))
    stop("cv_rf(): package 'ranger' is required.", call. = FALSE)
  # Raised HERE, in the parent, so that a forked worker inherits the
  # already-warned flag instead of each one raising its own copy.  See
  # .warn_include_coords().
  if (isTRUE(list(...)$include_coords)) .warn_include_coords()
  # `seed` has to reach the FOREST, not just the fold construction.
  # fit_rf_model() carries its own seed = 123L default and `...` did not
  # override it, so every fold's forest was grown with ranger seed 123
  # whatever the caller passed -- cv_rf(seed = 1) and cv_rf(seed = 2) returned
  # bit-identical predictions on fixed folds.  .cv_run_folds() runs each fold
  # under its own seeded stream, so drawing here gives a distinct,
  # reproducible seed per (seed, fold).  (`seed` is a formal of this function,
  # so it can never arrive in `...`; the branch that once tested for that was
  # dead code, and the doc that offered it as a way to fix one forest seed for
  # every fold described something that could not happen.)
  #
  # Thread policy under forked workers: each worker would otherwise read the
  # session's mc.cores opt-in for itself, so parallel = 4 with mc.cores = 8
  # meant 32 threads.  One thread per worker unless the caller says otherwise.
  dots <- list(...)
  # Quietly: cv_spatial() resolves `parallel` again below and says what it
  # does with it (a request above the machine's core count is capped, with a
  # message), so resolving it here aloud printed that message twice.
  n_workers <- suppressMessages(.resolve_n_cores(parallel))
  if (n_workers > 1L && !("num_threads" %in% names(dots)))
    dots$num_threads <- 1L
  fit_fn <- function(train_sf)
    do.call(fit_rf_model,
            c(list(train_sf, response_var, predictor_vars,
                   .already_prepped = TRUE,
                   seed = sample.int(.Machine$integer.max, 1L)),
              dots))
  # `pointize` is named here rather than left to `...`: it belongs to
  # cv_spatial(), and through `...` it reached ranger() as an unused argument.
  # `metrics` is named for the same reason as `pointize`: it belongs to
  # cv_spatial(), and through `...` it would reach ranger().
  # fold_info_fn carries each fold forest's out-of-bag count back, for one
  # warning per run in place of one per fold; see .rf_fold_oob_info().
  res <- cv_spatial(data_sf, response_var, predictor_vars, fit_fn = fit_fn,
                    .caller = "cv_rf",
                    folds = folds, k = k, seed = seed, boundary = boundary,
                    pointize = pointize,
                    block_size = block_size, auto_range = auto_range,
                    parallel = parallel, metrics = metrics,
                    fold_info_fn = .rf_fold_oob_info)
  .rf_cv_oob_note(res, dots)
}


# -----------------------------------------------------------------------------
# S3 methods
# -----------------------------------------------------------------------------

#' @rdname spatial_fit-predict
#' @export
predict.rf_fit <- function(object, newdata = NULL, ...) {
  if (!requireNamespace("ranger", quietly = TRUE))
    stop("predict.rf_fit(): package 'ranger' is required.", call. = FALSE)
  # Only ranger's own predict arguments pass; `data` and `object` are set
  # here, and anything else would vanish inside ranger's `...`.
  ranger_ok <- setdiff(names(formals(utils::getS3method("predict", "ranger",
                                                        envir = asNamespace("ranger")))),
                       c("object", "data", "..."))
  .check_dots(list(...), "predict.rf_fit", allowed = ranger_ok)
  if (is.null(newdata)) return(fitted(object))
  .check_predict_newdata(newdata, object, "predict.rf_fit")

  n_orig  <- nrow(newdata)
  # newdata is read in x and y when it mixes XY and XYZ (or XYM) geometries:
  # sf::st_transform() reads such a layer in the dimension of its first
  # geometry, so an XYZ first row stopped the move into the fit's CRS with
  # "OGR error".  Only predictions are returned, so nothing is said.
  newdata <- .xy_if_mixed(newdata)
  # A row with a coordinate that is not finite or lies beyond .POOL_XY_MAX
  # has no location, and is emptied here, as given: the move into the fit's
  # CRS below can bring it back finite (a northing of 1e200 in EPSG:32616
  # lands at a finite point in EPSG:32617), and it was then predicted there.
  # prep_model_data() drops it, so it comes back NA as the same row emptied.
  newdata <- .empty_unlocated(newdata)
  newdata <- .replay_crs_assumption(newdata, object$data_sf, "predict.rf_fit")
  newdata <- ensure_projected(newdata, target_crs = .crs_or_null(object$data_sf))
  newdata$..orig_row_id.. <- seq_len(n_orig)
  newdata <- prep_model_data(newdata, object$response_var,
                             object$predictor_vars, pointize = "auto",
                             require_response = FALSE)
  clean_idx <- seq_len(n_orig) %in% newdata$..orig_row_id..
  newdata$..orig_row_id.. <- NULL
  n_new <- nrow(newdata)

  # Nothing survived cleaning: every row is NA, as for any other dropped row
  # (and as predict.gwr_fit() and predict.bayesian_fit() do).  ranger would
  # otherwise be handed a zero-row frame and fail with a reason about
  # `sample_fraction`, which predict_surface() turns into an abort for a chunk
  # that merely lies outside the covariates' coverage.
  if (n_new == 0L) {
    .log_warn(paste0("predict.rf_fit(): every row of newdata was dropped as ",
                     "missing or non-finite; returning %d NA(s)."), n_orig)
    return(rep(NA_real_, n_orig))
  }

  include_coords <- isTRUE(object$info$include_coords)
  X <- .rf_frame(newdata, object$predictor_vars, include_coords,
                 .caller = "predict.rf_fit")
  # A factor level the forest was never grown with has no split to follow;
  # ranger 0.16 returns a plausible-looking number for it rather than erroring.
  # The helper intersects on name, so the whole training frame can be handed
  # over -- the added coordinate columns are numeric and simply skipped.
  X <- .rf_align_levels(X, sf::st_drop_geometry(object$data_sf))

  # Supply a seed unless the caller passed one.
  #
  # ranger:::predict.ranger.forest does
  #     if (is.null(seed)) seed <- runif(1, 0, .Machine$integer.max)
  # so every seedless predict() draws one uniform from the GLOBAL stream.
  # predict_surface() calls predict() once per chunk, which makes
  # `chunk_size` -- a pure performance knob that must not change any answer --
  # shift every subsequent random draw in the session: fold assignments
  # measurably differ between a script that called predict_surface() and one
  # that did not.  For a regression forest the seed changes nothing about the
  # prediction (it seeds ranger's own RNG, which only quantile/probability
  # prediction consumes), so pinning it costs nothing and makes this method
  # RNG-neutral.
  dots <- list(...)
  if (!("seed" %in% names(dots))) dots$seed <- 1L
  # Same thread policy as the fit: ranger's predict also defaults to every
  # core when num.threads is unset.
  if (!("num.threads" %in% names(dots)))
    dots$num.threads <- .sanitize_core_count(getOption("mc.cores", 1L))
  # A failure is an error naming ranger's reason.  It was caught, logged and
  # returned as an all-NA vector, so type = "quantiles" on a forest grown
  # without quantreg = TRUE, and type = "se" without keep.inbag = TRUE, both of
  # which the documentation promises are rejected, came back as NA with no R
  # condition, and model_metrics(newdata =) then reported n = 0.  Every caller
  # already handles an error: the cv_*() fold loop records it as the fold's
  # cause, and predict_surface() stops naming the rows.
  rr <- .call_capturing_stderr(function()
    do.call(stats::predict, c(list(object$engine, data = X), dots))$predictions)
  if (!is.null(rr$error))
    stop(sprintf("predict.rf_fit(): ranger's predict() failed: %s",
                 sub("^Error:\\s*", "", .stderr_reason(rr$error, rr$stderr))),
         call. = FALSE)
  p <- rr$value
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
  preds <- if (is.null(p)) rep(NA_real_, n_new)
           else .rf_unscale(as.numeric(p), object$info)
  .expand_predictions(preds, clean_idx, n_orig)
}


#' @rdname spatial_fit-fitted
#' @export
fitted.rf_fit <- function(object, ...) {
  p <- object$engine$predictions
  if (is.null(p) || length(p) != object$n) {
    .log_warn("fitted.rf_fit(): no usable out-of-bag predictions; returning NA.")
    return(rep(NA_real_, object$n))
  }
  .rf_unscale(as.numeric(p), object$info)
}


#' @rdname spatial_fit-residuals
#' @export
residuals.rf_fit <- function(object, ...) {
  y <- sf::st_drop_geometry(object$data_sf)[[object$response_var]]
  as.numeric(y) - fitted(object)
}


#' @rdname spatial_fit-coef
#' @export
coef.rf_fit <- function(object, ...) {
  stop("coef.rf_fit(): a random forest has no coefficients. For per-predictor ",
       "influence use `$info$importance` (", object$info$importance_type,
       " importance).", call. = FALSE)
}


#' @rdname print.spatial_fit
#' @export
print.rf_fit <- function(x, ...) {
  cat("<Random Forest (ranger)> spatial model fit\n")
  # paste(deparse(...)) and a CRS line, matching print.spatial_fit(): this
  # method takes precedence for an rf_fit, so without them the one backend
  # whose fits can carry the coordinates as predictors was the one that never
  # showed which space they were in, and a long formula printed as two mangled
  # fields.
  cat(sprintf("  Formula : %s\n", paste(deparse(x$formula), collapse = " ")))
  cat(sprintf("  n       : %d\n", x$n))
  cat(sprintf("  CRS     : %s\n", .fold_crs_label(x$data_sf)))
  cat(sprintf("  Trees   : %d (mtry = %s, min node = %s)\n",
              x$info$num_trees %||% NA_integer_,
              format(x$info$mtry %||% NA), format(x$info$min_node_size %||% NA)))
  cat(sprintf("  Coords as predictors: %s\n",
              if (isTRUE(x$info$include_coords))
                "YES - score with blocked folds; see ?fit_rf_model"
              else "no - no spatial information used; see ?fit_rf_model"))
  # Only when the fit recorded it and there is something to say.
  catp <- x$info$categorical_predictors
  ruf  <- x$info$respect_unordered_factors
  if (length(catp) > 0L && is.character(ruf) && length(ruf) == 1L && !is.na(ruf))
    cat(sprintf("  Categorical: %s split by %s (respect.unordered.factors = \"%s\")\n",
                paste(catp, collapse = ", "),
                switch(ruf,
                       order     = "levels ordered by mean response",
                       ignore    = "stored level order",
                       partition = "every two-way partition of the levels",
                       "ranger's rule"),
                ruf))
  # Only when the fit recorded it: an rf_fit built before these fields
  # existed prints as it did.
  if (!is.null(x$info$replace)) {
    frac <- .num1(x$info$sample_fraction)
    cat(sprintf("  Sampling: %s%s\n",
                if (isTRUE(x$info$replace)) "bootstrap, with replacement"
                else "subsample, without replacement",
                if (is.finite(frac)) sprintf(" (%.1f%% of rows per tree)", 100 * frac)
                else ""))
  }
  # .num1() rather than %||%: an unset ranger field arrives as numeric(0),
  # which is not NULL, so is.finite() would error and sprintf() would print
  # nothing at all.
  oob_rmse <- .num1(x$info$oob_rmse)
  # A forest with no row out of bag (replace = FALSE, sample_fraction = 1)
  # has no OOB error and NaN permutation importance.  The OOB line used to
  # vanish and the importance line to print empty, because sort() drops NaN.
  oob_p   <- x$engine$predictions
  no_oob  <- is.numeric(oob_p) && length(oob_p) > 0L && !any(is.finite(oob_p))
  if (is.finite(oob_rmse)) {
    cat(sprintf("  OOB RMSE: %.4f   OOB R^2: %.4f\n",
                oob_rmse, .num1(x$info$oob_r_squared)))
  } else if (no_oob) {
    cat("  OOB RMSE: undefined (no row is out of bag)\n")
  }
  imp <- x$info$importance
  if (!is.null(imp) && length(imp) > 0L) {
    # NaN is undefined; Inf is an importance whose squared units lie above
    # the largest double (fit_rf_model() says so in a WARN line) and is
    # printed as such.
    imp_ok <- imp[!is.na(imp)]
    cat(sprintf("  Importance (%s): %s\n", x$info$importance_type,
                if (!length(imp_ok)) {
                  if (no_oob) "undefined (no row is out of bag)" else "undefined"
                } else {
                  top <- utils::head(sort(imp_ok, decreasing = TRUE), 5L)
                  paste0(paste(sprintf("%s=%.4g", names(top), top), collapse = ", "),
                         if (any(is.infinite(top)))
                           " (Inf: above the largest double in squared units of the response)"
                         else "")
                }))
  }
  # What the out-of-bag error is the error of, and for which map it reads
  # low, as ?fit_rf_model says it: a direction with no target named reads as
  # true of every map.  A forest with no row out of bag has no such error, so
  # its print says only how to score it.
  if (no_oob) {
    cat("\n  With no out-of-bag error to read, score the forest with cv_rf() on\n")
    cat("  folds at the distances the map will predict at (?fold_separation).\n")
  } else {
    cat("\n  OOB is a random hold-out: under spatial autocorrelation it is the error\n")
    cat("  of predicting among the samples, and reads low for predictions made\n")
    cat("  farther from the data than the samples sit from each other. Score a map\n")
    cat("  that predicts farther with cv_rf() on folds at its distances\n")
    cat("  (?fold_separation).\n")
  }
  invisible(x)
}
