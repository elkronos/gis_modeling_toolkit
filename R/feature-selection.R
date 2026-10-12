# =============================================================================
# Forward feature selection with spatially blocked inner folds
# =============================================================================

#' Greedy forward feature selection with spatially blocked inner folds
#'
#' Selects predictors by repeatedly adding whichever candidate most improves a
#' cross-validated score, stopping when the best candidate's improvement is
#' no larger than the fold-to-fold noise in it (\code{rule = "se"}, the
#' default) or no larger than \code{tol} (\code{rule = "tol"}). The first
#' variable has a baseline to beat with every learner, so the selection can
#' be empty.
#'
#' \strong{The inner folds must be spatial, and that is the entire point.}
#' Nested selection is only worth doing if the inner loop is blocked the same
#' way the outer one is.  Random inner folds inside blocked outer folds select
#' variables that look predictive only because nearby points leak between
#' train and test.  The outer loop then reports honest-looking numbers for a
#' dishonestly chosen feature set, which is worse than not selecting at all,
#' because the dishonesty is now hidden behind a defensible-looking validation.
#' \code{method} therefore defaults to \code{"block_kfold"} and logs a loud
#' caution if set to \code{"random_kfold"} (a deliberate choice, so it is not
#' raised as an R warning).
#'
#' Call this \emph{inside} the \code{fit_fn} you pass to \code{cv_spatial()}.
#' \code{.cv_fit_one_fold()} calls \code{fit_fn(train_sf)} on the training
#' slice only, so anything done inside it is automatically nested and
#' leak-free; no extra plumbing is needed.  The cost grows fast: a sweep over
#' \code{p} candidates costs roughly \code{p^2 / 2 * k} model fits, and nesting
#' that inside \code{n} outer leave-one-out folds multiplies it by \code{n}.
#' \code{max_fits} guards against that.
#'
#' @section When a variable is accepted:
#' Each candidate set is scored by cross-validation on the same inner folds,
#' and the candidate with the best pooled score is the one considered. Under
#' \code{rule = "se"} it is accepted when two things hold: the pooled score
#' improves by more than \code{tol}, and the mean of the per-fold
#' improvements (the metric computed in each fold for the set with and
#' without the candidate, and differenced) exceeds \code{se_mult} times the
#' standard error of that mean across the folds. This is the
#' one-standard-error idea of Breiman et al. (1984) applied to a step: a gain
#' the folds do not agree on is not evidence that the variable predicts. The
#' fold is the unit because under spatially blocked folds the errors within a
#' block are correlated, so a standard error over rows would understate the
#' noise. With \code{tol = 0} and any gain accepted (the previous default,
#' still available as \code{rule = "tol"}), a GWR sweep over one real
#' predictor and six noise candidates added a noise variable in 11 runs of
#' 20; the default rule added one in 4 of 20 and kept the real predictor in
#' all 20. The rule costs power on weak predictors: a second real predictor
#' with a coefficient of 0.3 against 1 was kept in 14 runs of 20, against 18
#' of 20 when any gain was accepted. With few folds the standard error is
#' itself noisy, so raise \code{k} before lowering \code{se_mult}.
#' \code{history} records each candidate's pooled gain, mean per-fold gain
#' and its standard error, and \code{fold_scores} the per-fold scores
#' themselves, so every decision can be checked.
#'
#' @section The baseline for the first variable:
#' The first variable is compared with a model that has no predictors, which
#' is recorded as step 0 of \code{history}. Which model that is depends on
#' the learner, and \code{baseline$type} says which was used:
#' \describe{
#'   \item{\code{"null model"}}{\code{fit_fn(train_sf, character(0))}, when
#'     the learner can fit it. Step 1 is tested against it like any other
#'     step.}
#'   \item{\code{"null_fn"}}{the model \code{null_fn(train_sf)} returns, for
#'     a learner that cannot fit an empty set but has a natural
#'     no-predictor model. Step 1 is tested against it.}
#'   \item{\code{"fold mean"}}{the mean response of each fold's training
#'     rows, when neither is available; \code{\link{fit_rf_model}} refuses
#'     an empty predictor set. Here
#'     the first variable is \emph{provisional}: it is added without a test,
#'     the sweep continues by the rule, and the finished set is kept only if
#'     it beats the fold mean by the same rule; otherwise nothing is
#'     selected, and \code{baseline$path} holds the set that was reached.
#'     The test is at the end of the path, not at step 1, because a forest
#'     on a single variable can lose to the mean when two variables together
#'     beat it clearly (at a true \eqn{R^2} of 0.33 the best single-variable
#'     forest had RMSE 1.276 against 1.206 for the mean and 1.101 for the
#'     pair), so a test at step 1 would select nothing there. The first
#'     variable of a kept path then gets a test of its own; see \emph{The
#'     first variable under a fold-mean baseline}.}
#' }
#' Without a baseline, as before, at least one variable was always selected:
#' with a forest, in 10 draws of 10 on a response unrelated to every
#' candidate.
#'
#' @section The first variable under a fold-mean baseline:
#' The test at the end of the path is of the whole set, so it does not show
#' that the first variable belongs in it: two real predictors and an
#' unrelated candidate that happened to score best alone beat the mean as
#' the two real ones do. A kept path of two or more variables is therefore
#' scored once more without its first variable (one more set, \code{k} more
#' fits). That variable is dropped from \code{selected} when the set is no
#' better with it by the rule and the rest of the path still beats the fold
#' mean by the rule. \code{baseline$first} records the test and
#' \code{baseline$path} keeps the set as the sweep reached it, so
#' \code{selected} is \code{baseline$path[-1]} after a drop. A path of one
#' variable needs no second test: the test against the fold mean was the
#' test of that variable.
#'
#' This test is the one every later variable passed, an improvement of the
#' set by the rule, and it has the same limit. In 20 forest sweeps at a true
#' \eqn{R^2} of 0.33, with two real predictors and three unrelated
#' candidates, the path began with an unrelated candidate in 3: it was
#' dropped in 2 and kept in 1, where the forest scored better with it than
#' without (\emph{Selection with a random forest} says why it can). At a
#' true \eqn{R^2} of 0.15, 9 selections of 20 were kept, 3 of them began
#' with an unrelated candidate, and all 3 kept it for that reason. So in a
#' selection of two or more variables from a learner that falls back to the
#' fold mean, \code{selected[1]} has been shown to improve the set, not to
#' predict on its own. The test also costs power on a weak first predictor,
#' as the rule does at every step: with one real predictor of coefficient
#' 0.3 and five unrelated candidates, a linear learner that refuses the
#' empty set kept the real predictor in 77 sweeps of 100, against 80 before
#' this test existed and 77 for the same learner fitting its own null model
#' (in the 3 sweeps that differ it was the first variable, its gain was
#' within one standard error, and the learner with a null model selected
#' nothing). When the set is no better with the first variable but
#' the rest of the path does not beat the fold mean by the rule either, the
#' two cannot be told apart: the set is kept whole, a caution is logged, and
#' \code{baseline$first$note} says the variable is unproven. A learner that
#' fits its own null model gets no such test: its first variable passed one
#' at step 1.
#'
#' @section The baseline for a GWR:
#' For a GWR the fold mean is a weak baseline, because every GWR carries a
#' local intercept and beats a global mean on any smooth response whatever
#' its predictors are: on a smooth response unrelated to six candidates, a
#' GWR sweep against the fold mean selected one of them in 20 runs of 20.
#' Each of those paths had one variable, which the test of the first
#' variable above does not reach. The comparison that matters there is with
#' the kernel-weighted local mean at the same bandwidth, a GWR with a local
#' intercept and no predictors: against it the same sweeps selected a
#' variable in 4 runs of 20. \code{\link{fit_gwr_model}} fits that model
#' when its \code{predictor_vars} is \code{character(0)}, so a
#' \code{fit_fn} that hands its variables to \code{fit_gwr_model()} gets it
#' as the null model, and the first variable is tested against it at step 1
#' (\code{baseline$type} is \code{"null model"}). Give such a \code{fit_fn}
#' a fixed \code{bandwidth}, so that every set is fitted under the kernel
#' the baseline has; with \code{bandwidth = NULL} each set, the empty one
#' included, searches for its own. A \code{fit_fn} that refuses an empty
#' set falls back to the fold mean, and the sweep logs a caution: read a GWR
#' selection of one variable made that way as unproven, or supply the local
#' mean as \code{null_fn}, as the second half of the example at the end of
#' this page does for a learner written by hand.
#'
#' \code{predict()} on a GWR warns about held-out locations it cannot
#' estimate and about those beyond the reach of a fold's fit (classes
#' \code{"spatialkit_gwr_not_estimable"} and
#' \code{"spatialkit_gwr_extrapolation"}; see
#' \code{\link{predict.gwr_fit}}).  The sweep raises each once, under this
#' function's name, with the totals over the candidate sets it scored; every
#' set is predicted on the same rows, so a location can be counted once per
#' set.
#'
#' @section Selection with a random forest:
#' A forest confounds selection with tuning. Adding a predictor that carries
#' no signal still changes how the trees split, and with ranger's default
#' node size that regularises the forest: in one test the cross-validated
#' RMSE was 1.458 with the real predictor alone, 1.340 with a noise variable
#' added, and 1.341 with the real predictor alone and
#' \code{min_node_size = 20}. The noise variable bought what tuning the node
#' size buys. That gain is real, so no acceptance rule refuses it: with one
#' real predictor and six noise candidates a forest sweep kept the real one
#' in 20 runs of 20 and added 1.35 noise variables per run under the default
#' rule (2.4 per run, over 10 runs, when any gain was accepted). Setting the
#' node size in
#' \code{fit_fn} reduces this without removing it (0.8 per run at
#' \code{min_node_size = 20}). So a variable in a forest selection may be
#' there as a tuning effect wherever in the order it was added: later
#' additions with small gains are the usual case, and the first variable is
#' another, because a forest sweep starts from the fold mean and its first
#' variable is tested as an addition to the others (see \emph{The first
#' variable under a fold-mean baseline}). Compare the selected set with a
#' tuned forest on the variables you have other reasons to believe in:
#' \preformatted{
#' fit_fn <- function(train_sf, vars)
#'   fit_rf_model(train_sf, "z", vars, min_node_size = 20)
#' }
#' On a response unrelated to every candidate the same sweeps selected
#' nothing in 20 runs of 20 (19 of 20 with the larger node size).
#'
#' @section The score is not a performance estimate:
#' \code{$score} is the cross-validated \code{metric} of the winning set at
#' the final step: the best of every candidate set the sweep scored.  That
#' is the number the selection optimised, and a number optimised over many
#' candidates is optimistically biased by construction: Cawley and Talbot
#' (2010) show the bias can exceed the genuine differences between the models
#' being compared.  Quote it as the selection criterion, not as the
#' performance of the selected model.  An honest performance estimate needs
#' data the selection never saw.  Two ways to get one: run this function
#' inside the \code{fit_fn} of \code{\link{cv_spatial}()}, so the outer
#' folds score a model whose predictors were chosen on the inner ones alone;
#' or pass \code{select_on = "split"}, which selects on one spatial half of
#' \code{train_sf} and returns the selected set's score on the other as
#' \code{score_holdout}.
#'
#' @param train_sf Training data (\code{sf}).
#' @param response_var Character(1).
#' @param candidate_vars Character vector of predictors to choose among.
#' @param fit_fn A function \code{(train_sf, predictor_vars)} returning a
#'   \code{spatial_fit}.  The signature takes two arguments because selection
#'   has to refit with different predictor sets.
#' @param k Inner fold count. Default 5.
#' @param method Inner fold method. Default \code{"block_kfold"}.
#' @param block_size Passed to \code{make_folds()}; inherit the outer block
#'   size so inner and outer blocks are on the same spatial scale.
#' @param metric Score to optimise: \code{"RMSE"}, \code{"MAE"} (minimised) or
#'   \code{"R2"} (maximised). Default \code{"RMSE"}.
#' @param tol Minimum improvement of the pooled score required to accept a
#'   variable, in the units of \code{metric}.  Default 0.  It is a floor
#'   under both rules: with \code{rule = "tol"} it is the whole test, and
#'   with \code{rule = "se"} the per-fold test applies as well. The first
#'   variable is judged against a baseline with every learner; see \emph{The
#'   baseline for the first variable}.  Any single number is accepted:
#'   \code{Inf} leaves nothing selected, and \code{-Inf} puts no floor under
#'   a step.
#' @param max_vars Optional cap on how many predictors to select.
#' @param max_fits Abort if the sweep would exceed this many model fits.
#'   Default 5000.
#' @param seed RNG seed.  It governs both the inner fold construction and the
#'   cross-validation itself: it is forwarded to \code{cv_spatial(seed = )},
#'   which draws one RNG stream per fold from it, so it also seeds the
#'   \emph{learner} inside every fold.  A stochastic \code{fit_fn} is therefore
#'   reproducible from this one value.
#' @param quiet Logical; suppress this function's progress \code{message()}s.
#'   It does not silence R warnings, nor the package's console log echo
#'   (see \code{\link{spatialkit_quiet}} for that). Default \code{FALSE}.
#' @param auto_range Logical.  If \code{TRUE} and \code{method} is
#'   \code{"block_kfold"}, the autocorrelation range is estimated from the
#'   response (detrended on the candidates) and used as the minimum block
#'   size of the inner folds, exactly as in \code{\link{make_folds}()};
#'   \code{block_size} still applies as a floor.  Default \code{FALSE}, which
#'   keeps the geometric blocks.  Inner blocks smaller than the range let a
#'   candidate be selected for spatial proximity to the response rather than
#'   for predicting it, which is the same failure the \code{random_kfold}
#'   caution above exists to prevent, so the leakage warning
#'   \code{make_folds()} raises applies here with more force than usual.
#' @param select_on \code{"all"} (default) runs the sweep on every row of
#'   \code{train_sf}.  \code{"split"} runs it on one spatially blocked half,
#'   then fits the selected set on that half and scores it on the other:
#'   \code{score_holdout} is then the selected model's \code{metric} on rows
#'   whose response the sweep never read (the sweep's own \code{score} is
#'   not; see "The score is not a performance estimate").  It is one
#'   estimate from one region: the halves share a border with no buffer, so
#'   rows near it are still correlated with the selection half.
#'   Both halves come back in \code{$split}.  See the "Post-selection
#'   inference" section of \code{\link{determine_optimal_levels}} for the
#'   trade: coverage for half the sample.
#' @param rule How a step is judged. \code{"se"} (default) accepts the best
#'   candidate when its pooled gain exceeds \code{tol} and its mean per-fold
#'   gain exceeds \code{se_mult} standard errors of that mean;
#'   \code{"tol"} accepts any pooled gain above \code{tol}, the behaviour
#'   before this argument existed. See \emph{When a variable is accepted}.
#' @param se_mult Number of standard errors the mean per-fold gain must
#'   exceed under \code{rule = "se"}. Default 1; 0 asks only that the mean
#'   per-fold gain be positive.
#' @param null_fn Optional function \code{(train_sf)} returning a fitted
#'   model with no predictors, used as the baseline for the first variable
#'   when \code{fit_fn} cannot fit an empty predictor set. Default
#'   \code{NULL}: \code{fit_fn(train_sf, character(0))} is tried, and the
#'   training-fold mean is used when that fails. The example builds one.
#' @return A list of class \code{"feature_selection"} (so that
#'   \code{\link{plot.feature_selection}()} draws the selection path) with
#'   \code{selected} (the chosen predictors, in the order they were added;
#'   possibly none), \code{score}, \code{score_holdout}, \code{history},
#'   \code{params}, \code{split}, \code{baseline} and \code{fold_scores}.
#'   \code{score} is the winning set's cross-validated \code{metric} at the
#'   final step: the \strong{selection-internal} optimum, optimistically
#'   biased because it was chosen as the best of many (see the section above),
#'   and \code{NA} when nothing was selected.  \code{history} is a data.frame
#'   with \code{step}, \code{variable}, \code{score}, \code{n_pred},
#'   \code{gain}, \code{fold_gain} and \code{fold_gain_se},
#'   holding every candidate evaluated at every step, and a \code{step = 0}
#'   row named \code{"<none>"} giving the baseline, so the first variable's
#'   gain can be read off directly.  \code{gain} is the candidate set's
#'   improvement in the pooled score over the set it was added to (the
#'   baseline at step 1), \code{fold_gain} the mean of the per-fold
#'   improvements and \code{fold_gain_se} its standard error: the three
#'   numbers the rule reads.  \code{fold_scores} is a data.frame with
#'   \code{step}, \code{variable}, \code{fold} and \code{score}: the metric
#'   of every set in every fold.  \code{baseline} is a list with
#'   \code{type} (\code{"null model"}, \code{"null_fn"} or
#'   \code{"fold mean"}), \code{score}, \code{path} (the set the sweep
#'   reached), \code{gain}, \code{fold_gain} and \code{fold_gain_se} (that
#'   set against the baseline), \code{kept} (\code{FALSE} when the path
#'   did not beat a fold-mean baseline and \code{selected} is empty for that
#'   reason; \code{NA} when the path is empty) and \code{first}.
#'   \code{first} is \code{NULL} unless the baseline is the fold mean and a
#'   path of two or more variables was kept; it is then the test of the
#'   first variable, a list with \code{variable}; \code{gain},
#'   \code{fold_gain} and \code{fold_gain_se} (the path against the path
#'   without that variable); \code{rest}, a list with \code{vars},
#'   \code{score}, \code{gain}, \code{fold_gain} and \code{fold_gain_se}
#'   (the path without it, and that set against the baseline);
#'   \code{dropped} (\code{TRUE} when the variable was taken out:
#'   \code{selected} is then \code{path[-1]} and \code{score} is
#'   \code{rest$score}); and \code{note}, the verdict in words.  The set
#'   scored for that test is not a step of the sweep, so it is not in
#'   \code{history} or \code{fold_scores}.  Every set is scored on the
#'   same rows: those the null model's cross-validation predicted or, when
#'   the baseline is the fold mean, those any step-1 set predicted;
#'   \code{params$n_scored} counts them (a warning says so when that is
#'   fewer than all).  \code{n_pred} is how many
#'   rows the set's cross-validation predicted.  A set that left some of the
#'   scored rows unpredicted, because a fold failed for it, has \code{score}
#'   \code{NA}, with a warning naming it: scored on the rows it did predict
#'   it would be compared on fewer, usually easier, rows than its rivals.  A
#'   factor with a level found in one spatial block only is the usual case,
#'   and cannot be selected.
#'   \code{score_holdout} is \code{NA} unless \code{select_on = "split"}, and
#'   then the selected set's \code{metric} when fitted on the selection half
#'   and predicted on the estimation half (\eqn{R^2} against the selection
#'   half's mean, the out-of-sample convention); \code{NA} when nothing was
#'   selected or the prediction failed.  \code{split} is \code{NULL} or a
#'   list with \code{selection} and \code{estimation}, integer row positions
#'   in \code{train_sf} as passed; rows the completeness filter above dropped
#'   are in neither.
#'   \code{params} records \code{metric}, \code{method}, \code{k},
#'   \code{tol}, \code{seed}, \code{auto_range}, \code{select_on},
#'   \code{n_candidates}, \code{estimated_fits}, \code{n_scored},
#'   \code{rule}, \code{se_mult} and \code{baseline} (its type).
#' @references
#' Breiman, L., Friedman, J. H., Olshen, R. A. and Stone, C. J. (1984).
#' \emph{Classification and Regression Trees}. Wadsworth.
#'
#' Cawley, G. C. and Talbot, N. L. C. (2010). On over-fitting in model
#' selection and subsequent selection bias in performance evaluation.
#' \emph{Journal of Machine Learning Research}, 11, 2079--2107.
#' \url{https://jmlr.org/papers/v11/cawley10a.html}
#' @family cross-validation
#' @examples
#' if (requireNamespace("sp", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 120
#'   pts <- st_as_sf(
#'     data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
#'                a = rnorm(n), b = rnorm(n), noise = rnorm(n)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   pts$resp <- 3 * pts$a + 2 * pts$b + rnorm(n, 0, 0.3)
#'
#'   # The learner: a GWR with an adaptive bisquare kernel reaching the 30th
#'   # nearest training point.  With no variables it is the kernel-weighted
#'   # local mean, which is the baseline: a GWR's local intercept beats a
#'   # global mean on any smooth response, so the first variable has to beat
#'   # the local mean under the same kernel.
#'   fit_fn <- function(tr, vars) fit_gwr_model(tr, "resp", vars, bandwidth = 30)
#'   sel <- select_features_forward(pts, "resp", c("a", "b", "noise"), fit_fn,
#'                                  k = 3, quiet = TRUE)
#'   sel$selected
#'   sel$baseline$type        # "null model": step 1 was tested against it
#'   sel$history[sel$history$step <= 1, c("step", "variable", "score", "gain")]
#'
#'   # A learner that cannot fit an empty set takes its no-predictor model
#'   # through `null_fn`.  Here the same local mean, written as a spatial_fit
#'   # of its own with a predict() method, for a learner that insists on a
#'   # variable.  Without `null_fn` that sweep would start from the mean of
#'   # each training fold.
#'   needs_one <- function(tr, vars) {
#'     if (!length(vars)) stop("this learner needs a predictor")
#'     fit_gwr_model(tr, "resp", vars, bandwidth = 30)
#'   }
#'   local_mean <- function(tr, bandwidth = 30) {
#'     new_spatial_fit("local_mean_fit",
#'                     engine = list(xy = st_coordinates(tr), y = tr$resp,
#'                                   bandwidth = bandwidth),
#'                     formula = resp ~ 1, response_var = "resp",
#'                     predictor_vars = character(0), data_sf = tr)
#'   }
#'   predict.local_mean_fit <- function(object, newdata = NULL, ...) {
#'     e  <- object$engine
#'     at <- if (is.null(newdata)) e$xy else st_coordinates(newdata)
#'     vapply(seq_len(nrow(at)), function(i) {
#'       d <- sqrt((e$xy[, 1] - at[i, 1])^2 + (e$xy[, 2] - at[i, 2])^2)
#'       h <- sort(d)[min(e$bandwidth, length(d))]  # the kernel's reach here
#'       w <- ifelse(d > h, 0, (1 - d^2 / h^2)^2)   # bisquare weights
#'       sum(w * e$y) / sum(w)
#'     }, numeric(1))
#'   }
#'   registerS3method("predict", "local_mean_fit", predict.local_mean_fit)
#'   sel2 <- select_features_forward(pts, "resp", c("a", "b", "noise"), needs_one,
#'                                   k = 3, quiet = TRUE, null_fn = local_mean)
#'   sel2$baseline$type       # "null_fn"
#'   all.equal(sel2$baseline$score, sel$baseline$score)  # the same baseline
#' }
#' @export
select_features_forward <- function(train_sf, response_var, candidate_vars,
                                    fit_fn, k = 5,
                                    method = c("block_kfold", "random_kfold"),
                                    block_size = NULL,
                                    metric = c("RMSE", "MAE", "R2"),
                                    tol = 0, max_vars = NULL,
                                    max_fits = 5000L, seed = 123,
                                    quiet = FALSE, auto_range = FALSE,
                                    select_on = c("all", "split"),
                                    rule = c("se", "tol"), se_mult = 1,
                                    null_fn = NULL) {
  # Every candidate set is a cv_spatial() of the same rows, and each prepares
  # them for itself.  A candidate that takes one value in every row is
  # reported by each set that holds it, in the same words (class
  # spatialkit_constant_predictor): raised by the first, not repeated.
  end_run <- .once_per_run(.prep_warning_classes)
  on.exit(end_run(), add = TRUE)
  # A GWR learner counts, in every cross-validation of the sweep, the
  # held-out locations it cannot estimate and those beyond the reach of a
  # fold's fit, and cv_spatial() raises each count once per candidate set.
  # Collected over the sweep, each is raised once, with the totals.  Every
  # set is predicted on the same rows, so a location can be counted once per
  # set.  No other learner signals either class.
  .gwr_tally_predict_warnings(
    .select_features_forward(
      train_sf, response_var, candidate_vars, fit_fn, k = k, method = method,
      block_size = block_size, metric = metric, tol = tol,
      max_vars = max_vars, max_fits = max_fits, seed = seed, quiet = quiet,
      auto_range = auto_range, select_on = select_on, rule = rule,
      se_mult = se_mult, null_fn = null_fn),
    caller = "select_features_forward", refits = TRUE)
}

# The sweep itself; select_features_forward() is this with the counted
# warnings of a GWR learner collected over it.
.select_features_forward <- function(train_sf, response_var, candidate_vars,
                                     fit_fn, k = 5,
                                     method = c("block_kfold", "random_kfold"),
                                     block_size = NULL,
                                     metric = c("RMSE", "MAE", "R2"),
                                     tol = 0, max_vars = NULL,
                                     max_fits = 5000L, seed = 123,
                                     quiet = FALSE, auto_range = FALSE,
                                     select_on = c("all", "split"),
                                     rule = c("se", "tol"), se_mult = 1,
                                     null_fn = NULL) {
  method <- match.arg(method)
  metric <- match.arg(metric)
  select_on <- match.arg(select_on)
  rule <- match.arg(rule)
  .check_scalar(se_mult, "se_mult", "select_features_forward", min = 0,
                what = "a single non-negative number")
  # Any one number that is not missing.  `tol = Inf` asks for nothing to be
  # selected and `-Inf` puts no floor under a step; both were accepted before
  # `tol` was checked at all, so .check_scalar(), which wants a finite value,
  # is not the test.
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol))
    stop(sprintf(paste0("select_features_forward(): `tol` must be a single ",
                        "number (Inf and -Inf are allowed); got %s."),
                 if (length(tol) != 1L)
                   sprintf("%s of length %d", class(tol)[1L], length(tol))
                 else if (is.numeric(tol)) format(tol)
                 else class(tol)[1L]), call. = FALSE)
  if (!is.null(null_fn) && !is.function(null_fn))
    stop("select_features_forward(): `null_fn` must be NULL or a function of ",
         "(train_sf) returning a fitted model with no predictors.",
         call. = FALSE)
  .msg <- function(...) if (!quiet) message(...)

  if (!inherits(train_sf, "sf"))
    stop("select_features_forward(): `train_sf` must be an sf object.", call. = FALSE)
  if (!is.function(fit_fn))
    stop("select_features_forward(): `fit_fn` must be a function of ",
         "(train_sf, predictor_vars).", call. = FALSE)
  candidate_vars <- unique(candidate_vars)
  missing_v <- setdiff(candidate_vars, names(train_sf))
  if (length(missing_v) > 0L)
    stop("select_features_forward(): candidate(s) absent from `train_sf`: ",
         paste(sQuote(missing_v), collapse = ", "), call. = FALSE)
  if (length(candidate_vars) == 0L)
    stop("select_features_forward(): no candidate variables supplied.", call. = FALSE)

  # Every candidate set must be scored on the SAME observations.  cv_spatial()
  # re-runs prep_model_data() per candidate, which drops rows that are NA in
  # THAT candidate's columns -- so without this, a step compares an RMSE over
  # 150 rows against one over 120, and can prefer a variable for having an
  # easier surviving subset rather than for predicting better.
  #
  # The test must match prep_model_data()'s exactly.  complete.cases() alone
  # lets Inf through, while prep_model_data() also drops on all(is.finite(r))
  # over the numeric columns -- so a candidate carrying a single Inf would
  # still reach cv_spatial() with a smaller row set than its rivals, which is
  # precisely the comparison this filter exists to prevent.
  fs_df <- sf::st_drop_geometry(train_sf)[, c(response_var, candidate_vars),
                                          drop = FALSE]
  keep <- stats::complete.cases(fs_df)
  num_mask <- vapply(fs_df, is.numeric, logical(1))
  if (any(num_mask))
    keep <- keep & apply(as.matrix(fs_df[, num_mask, drop = FALSE]), 1L,
                         function(r) all(is.finite(r)))
  if (!all(keep)) {
    if (sum(keep) < 2L)
      stop("select_features_forward(): fewer than two rows are complete and ",
           "finite across the response and every candidate.", call. = FALSE)
    .log_warn(paste0("select_features_forward(): dropping %d row(s) incomplete ",
                     "or non-finite across the response and candidates, so ",
                     "every candidate set is scored on the same %d ",
                     "observations."),
              sum(!keep), sum(keep))
    train_sf <- train_sf[keep, , drop = FALSE]
  }

  # Sample splitting: the sweep sees the selection half only; the estimation
  # half scores the chosen set afterwards.  .spatial_half_split()'s positions
  # index train_sf as it stands here, after the completeness filter.
  split <- NULL
  holdout_sf <- NULL
  if (identical(select_on, "split")) {
    if (!all(sf::st_geometry_type(train_sf, by_geometry = TRUE) == "POINT"))
      train_sf <- coerce_to_points(train_sf, "auto")
    split <- .spatial_half_split(train_sf, seed = seed,
                                 caller = "select_features_forward")
    holdout_sf <- train_sf[split$estimation, , drop = FALSE]
    train_sf   <- train_sf[split$selection, , drop = FALSE]
    # Returned as positions in the layer the CALLER passed, as print() says
    # and as determine_optimal_levels() returns them.  They were positions
    # after the completeness filter, so with 7 incomplete rows
    # pts[fs$split$estimation, ] held 47 selection-half rows and 4 of the
    # dropped ones: the one use of the split, estimating on rows the
    # selection never saw, got rows it had seen.
    kept <- which(keep)
    split$selection  <- kept[split$selection]
    split$estimation <- kept[split$estimation]
    .msg(sprintf("select_features_forward(): selecting on %d points, scoring the result on the other %d.",
                 nrow(train_sf), nrow(holdout_sf)))
  }

  if (identical(method, "random_kfold"))
    .log_warn(paste0("select_features_forward(): random_kfold inner folds ",
                     "leak spatial autocorrelation between train and test. ",
                     "Variables can then be selected for being spatially ",
                     "close to the response rather than predictive of it, and ",
                     "an outer blocked loop will report honest-looking numbers ",
                     "for a dishonestly chosen feature set. Use ",
                     "'block_kfold' unless you know why you are not."))

  p <- length(candidate_vars)
  max_steps <- if (is.null(max_vars)) p else min(p, as.integer(max_vars))
  # sweep cost: p + (p-1) + ... over at most max_steps rounds, times k folds
  est_fits <- sum(p - seq_len(max_steps) + 1L) * k
  if (est_fits > max_fits)
    stop(sprintf(paste0("select_features_forward(): this sweep would fit about ",
                        "%d models (%d candidates, %d folds), above max_fits = %d. ",
                        "Reduce candidates, lower k, set max_vars, or raise ",
                        "max_fits deliberately. Nesting a sweep inside ",
                        "leave-one-out outer folds multiplies this by n."),
                 est_fits, p, k, max_fits), call. = FALSE)

  better <- if (metric == "R2") function(a, b) a > b else function(a, b) a < b
  worst  <- if (metric == "R2") -Inf else Inf

  # Build the inner folds ONCE, before the sweep.
  #
  # Fold construction does not depend on the candidate set: make_folds() reads
  # `predictor_vars` only to detrend the range estimate (under auto_range, and
  # for the leakage diagnostic), which uses the full candidate list whatever
  # the sweep has selected so far, and every other argument is fixed across
  # the sweep.  Rebuilding them inside score_set() therefore produced the same
  # splits p^2/2 times -- repeating every block-size warning as many times --
  # and, worse, put make_folds() OUTSIDE the try() below, so a single
  # fold-construction failure (block_kfold raising when the geometry collapses
  # to one block, say) killed the whole sweep instead of the candidate being
  # scored NA.  It cannot be a per-candidate NA in any case: if the folds
  # cannot be built, no candidate is scorable, so this is one informative
  # error instead of p^2/2 silent ones.
  folds <- try(make_folds(train_sf, k = k, method = method, seed = seed,
                          block_size = block_size, auto_range = auto_range,
                          response_var = response_var,
                          predictor_vars = candidate_vars), silent = TRUE)
  if (inherits(folds, "try-error"))
    stop("select_features_forward(): could not build the inner CV folds, so no ",
         "candidate can be scored: ",
         trimws(conditionMessage(attr(folds, "condition"))),
         " Adjust `method`, `k` or `block_size`.", call. = FALSE)

  # Every set is scored on ONE fixed row set, `ref_ids`: the rows the null
  # model predicted, or, when there is no null model (a forest refuses an
  # empty predictor set), the rows any step-1 candidate predicted.
  # cv_spatial()'s `overall` pools whichever folds survived, so a set whose
  # fit or predict failed on a fold -- a factor level found in one block
  # only, the ordinary case -- was scored on fewer rows, usually the easier
  # ones, and could win for that alone: a noise factor was chosen over the
  # true driver (RF RMSE 2.35 on 192 rows against 2.63 on 250).  A set that
  # leaves a reference row unpredicted is scored NA instead, and said so.
  # Rows the reference runs did not predict (a fold that fails for every set,
  # from geometry rather than predictors) drop out of every score alike.
  ref_ids <- NULL
  # The classes of the models the learner returns, noted as they are fitted:
  # the fold-mean baseline means something different for a GWR (see below).
  fit_classes <- character(0)
  run_set <- function(vars, inner_fit = function(tr) fit_fn(tr, vars)) {
    noting_fit <- function(tr) {
      f <- inner_fit(tr)
      fit_classes <<- union(fit_classes, class(f))
      f
    }
    # cv_spatial()'s own partial-failure warning is muffled: the consequence
    # for the sweep is reported below, once per step, in the sweep's terms.
    res <- try(withCallingHandlers(
      suppressMessages(
        cv_spatial(train_sf, response_var, vars, fit_fn = noting_fit,
                   folds = folds, seed = seed)),
      warning = function(w)
        if (.is_failed_folds_warning(w)) invokeRestart("muffleWarning")),
      silent = TRUE)
    if (inherits(res, "try-error") || !is.data.frame(res$predictions))
      return(NULL)
    pr <- res$predictions
    pr[is.finite(pr$y) & is.finite(pr$yhat), , drop = FALSE]
  }
  n_pred_of <- function(pr) if (is.null(pr)) 0L else nrow(pr)
  covers_ref <- function(pr) !is.null(pr) && all(ref_ids %in% pr$`..row_id`)
  score_on <- function(pr) {
    if (!length(ref_ids) || !covers_ref(pr)) return(NA_real_)
    val <- .cv_overall_metrics(pr[pr$`..row_id` %in% ref_ids, , drop = FALSE])[[metric]]
    if (is.null(val) || !is.finite(val)) NA_real_ else as.numeric(val)
  }
  # The same metric fold by fold, on the reference rows: what the "se" rule
  # reads.  Named by fold; NULL for a set that does not cover the rows.
  fold_scores_of <- function(pr) {
    if (!length(ref_ids) || !covers_ref(pr) || is.null(pr$fold)) return(NULL)
    sub <- pr[pr$`..row_id` %in% ref_ids, , drop = FALSE]
    vapply(split(seq_len(nrow(sub)), sub$fold), function(ix) {
      v <- .cv_overall_metrics(sub[ix, , drop = FALSE])[[metric]]
      if (is.null(v) || !is.finite(v)) NA_real_ else as.numeric(v)
    }, numeric(1))
  }
  # Paired per-fold improvement of set `a` over set `b`: its mean and the
  # standard error of that mean across the folds both sets scored.
  paired_gain <- function(pr_a, pr_b) {
    fa <- fold_scores_of(pr_a); fb <- fold_scores_of(pr_b)
    out <- list(mean = NA_real_, se = NA_real_, n = 0L)
    if (is.null(fa) || is.null(fb)) return(out)
    cm <- intersect(names(fa), names(fb))
    d  <- if (metric == "R2") fa[cm] - fb[cm] else fb[cm] - fa[cm]
    d  <- d[is.finite(d)]
    out$n <- length(d)
    if (length(d) >= 1L) out$mean <- mean(d)
    if (length(d) >= 2L) out$se <- stats::sd(d) / sqrt(length(d))
    out
  }
  # Does a gain count?  Always more than `tol` on the pooled score; under
  # rule = "se" also a mean per-fold gain above se_mult standard errors.
  # With fewer than two folds scored there is no standard error, and the
  # step is judged on `tol` alone, said once.
  se_missing_said <- FALSE
  accepts <- function(gain, pg) {
    if (!(is.finite(gain) && gain > tol)) return(FALSE)
    if (identical(rule, "tol")) return(TRUE)
    if (!is.finite(pg$se) || !is.finite(pg$mean)) {
      if (!se_missing_said) {
        se_missing_said <<- TRUE
        .warn_and_log(paste0(
          "select_features_forward(): rule = \"se\" needs the %s of at least ",
          "two folds to compare two sets, and %d could be scored for both; ",
          "such steps are judged on `tol` alone."), metric, pg$n)
      }
      return(TRUE)
    }
    pg$mean > se_mult * pg$se
  }
  set_ref <- function(ids, from) {
    ref_ids <<- ids
    if (length(ids) && length(ids) < nrow(train_sf))
      .warn_and_log(paste0(
        "select_features_forward(): %s predicted only %d of the %d rows, so ",
        "every candidate set is scored on those %d; the rest drop out of ",
        "every score alike."),
        from, length(ids), nrow(train_sf), length(ids))
  }

  selected  <- character(0)
  remaining <- candidate_vars
  history   <- list()
  fold_hist <- list()
  keep_folds <- function(step, variable, pr) {
    fs <- fold_scores_of(pr)
    if (is.null(fs)) return(invisible(NULL))
    fold_hist[[length(fold_hist) + 1L]] <<- data.frame(
      step = step, variable = variable, fold = names(fs), score = unname(fs),
      stringsAsFactors = FALSE)
  }

  # Step 0: something for the first variable to beat.  `best` used to start
  # at Inf/-Inf with the stopping test gated on is.finite(best), which made
  # the first step unconditional whenever the empty set could not be scored
  # -- and fit_rf_model() refuses a zero-length predictor_vars, as
  # fit_gwr_model() did.  With either as the learner a variable was therefore
  # always selected: in 10 draws of 10 on a response unrelated to every
  # candidate, 2.6 to 3.1 variables on average.  There are three baselines,
  # tried in this order:
  #   "null_fn"     the caller's own no-predictor model;
  #   "null model"  fit_fn(train_sf, character(0)), when the learner can;
  #   "fold mean"   the training-fold mean, which needs no model at all.
  # The first two gate step 1 like any other step.  The fold mean does not:
  # a one-variable forest loses to it even when two variables together beat
  # it clearly (true R2 0.33: 1.276 against 1.206 for the best single
  # variable, 1.101 for the pair), so gating step 1 on it turned "always
  # selects something" into "never selects anything" at moderate signal.
  # There the first variable is provisional, and the baseline is applied to
  # the finished path instead: the set is kept only if it beats the fold
  # mean by the same rule every step is held to, and its first variable
  # only if the set is better with it (the block after the sweep).
  #
  # A backend that cannot fit an intercept-only model makes cv_spatial()
  # report "all folds failed", which is expected here.  suppressWarnings()
  # silences the R condition, but the k "fold i fit failed" lines and the
  # "all k folds failed" line are LOGGER records, which no condition handler
  # touches: every successful RF/GWR run printed them, identical to a
  # genuinely failed run's, even with quiet = TRUE.  Raise the console
  # threshold for the probe alone; the file trace (index 1) keeps the lines,
  # where a diagnostic belongs.  The rows the null model predicts, when it
  # can be fitted, are the rows every set is scored on.
  baseline_type <- NA_character_
  null_score    <- NA_real_
  base_pr       <- NULL
  if (!is.null(null_fn)) {
    null_run <- run_set(character(0), inner_fit = function(tr) null_fn(tr))
    if (n_pred_of(null_run) == 0L)
      stop("select_features_forward(): `null_fn` gave no cross-validated ",
           "prediction on any fold. It must be a function of (train_sf) that ",
           "returns a fitted model with a predict() method, as `fit_fn` does.",
           call. = FALSE)
    baseline_type <- "null_fn"
  } else {
    null_run <- logger::with_log_threshold(
      suppressWarnings(run_set(character(0))),
      threshold = logger::FATAL, namespace = "spatialkit", index = 2)
    if (n_pred_of(null_run) > 0L) baseline_type <- "null model"
  }
  if (!is.na(baseline_type)) {
    set_ref(sort(unique(null_run$`..row_id`)),
            if (identical(baseline_type, "null_fn")) "the `null_fn` model"
            else "the null (intercept-only) model")
    null_score <- score_on(null_run)
    if (is.finite(null_score)) base_pr <- null_run
    else baseline_type <- NA_character_      # fall through to the fold mean
  }
  best    <- if (is.finite(null_score)) null_score else worst
  best_pr <- base_pr
  step0 <- function() {
    history <<- c(list(data.frame(
      step = 0L, variable = "<none>", score = unname(null_score),
      n_pred = n_pred_of(base_pr), gain = NA_real_, fold_gain = NA_real_,
      fold_gain_se = NA_real_, stringsAsFactors = FALSE)), history)
    fs <- fold_scores_of(base_pr)
    if (!is.null(fs))
      fold_hist <<- c(list(data.frame(
        step = 0L, variable = "<none>", fold = names(fs), score = unname(fs),
        stringsAsFactors = FALSE)), fold_hist)
  }
  if (!is.na(baseline_type)) step0()

  repeat {
    if (length(remaining) == 0L || length(selected) >= max_steps) break

    runs <- lapply(remaining, function(v) run_set(c(selected, v)))
    if (is.null(ref_ids))
      set_ref(sort(unique(unlist(lapply(runs, `[[`, "..row_id")))),
              "the step-1 candidate sets together")

    # No baseline from the learner: the training-fold mean, read off the
    # runs just made.  cv_spatial() records each held-out row's training-fold
    # mean of the response beside its prediction, and every candidate set
    # was scored on the same rows and folds, so this is that mean scored on
    # the sweep's own folds with no further fit.
    if (is.na(baseline_type)) {
      bp <- .fs_fold_mean_predictions(runs)
      sc <- score_on(bp)
      if (is.finite(sc)) {
        baseline_type <- "fold mean"
        base_pr <- bp; null_score <- sc
        best <- sc; best_pr <- bp
        step0()
        .msg(sprintf(paste0(
          "select_features_forward(): `fit_fn` cannot fit a model with no ",
          "predictors, so the baseline is the training-fold mean (%s = %.4f). ",
          "The first variable is provisional: the finished set is kept only ",
          "if it beats that baseline, and its first variable only if the set ",
          "is better with it."), metric, sc))
        # Every GWR carries a local intercept, so it beats a global mean on
        # any smooth response whatever its predictors are: against the fold
        # mean a GWR sweep kept one of six unrelated candidates in 20 runs
        # of 20, against 4 of 20 with the kernel-weighted local mean as the
        # baseline.  fit_gwr_model() fits that model for character(0), so a
        # GWR sweep arrives here only when its `fit_fn` refuses the empty
        # set itself.  The caller is told what the baseline does and does
        # not show, and the two ways to the right one.
        if ("gwr_fit" %in% fit_classes)
          .log_warn(paste0(
            "select_features_forward(): the learner is a GWR and the ",
            "baseline is the training-fold mean. A GWR beats a global mean ",
            "on any smooth response through its local intercept alone, so ",
            "passing this baseline does not show that the first variable ",
            "predicts. The comparison that does is with the kernel-weighted ",
            "local mean at the same bandwidth (an intercept-only GWR): let ",
            "`fit_fn` pass character(0) on to fit_gwr_model(), which fits ",
            "it, or pass a function fitting it as `null_fn`. The example in ",
            "?select_features_forward builds one."))
      }
    }

    step_scores <- vapply(runs, score_on, numeric(1))
    names(step_scores) <- remaining
    n_pred <- vapply(runs, n_pred_of, integer(1))
    step_gain <- if (metric == "R2") step_scores - best else best - step_scores
    step_gain[!is.finite(step_gain)] <- NA_real_
    step_pg <- lapply(runs, paired_gain, pr_b = best_pr)
    history[[length(history) + 1L]] <- data.frame(
      step = length(selected) + 1L, variable = remaining,
      score = unname(step_scores), n_pred = n_pred,
      gain = unname(step_gain),
      fold_gain = vapply(step_pg, `[[`, numeric(1), "mean"),
      fold_gain_se = vapply(step_pg, `[[`, numeric(1), "se"),
      stringsAsFactors = FALSE
    )
    for (j in seq_along(runs))
      keep_folds(length(selected) + 1L, remaining[j], runs[[j]])
    # Sets that predicted some reference rows but not all.  One predicting
    # none has already raised cv_spatial()'s "all folds failed" warning.
    short <- which(n_pred > 0L & !vapply(runs, covers_ref, logical(1)))
    if (length(short))
      .warn_and_log(paste0(
        "select_features_forward(): step %d: %s left some of the %d scored ",
        "rows unpredicted (a fold failed: a factor level found in one block ",
        "only, say), so %s scored NA rather than on fewer, easier rows."),
        length(selected) + 1L,
        paste(vapply(short, function(j) sprintf(
          "{%s} (%d of them predicted)",
          paste(c(selected, remaining[j]), collapse = ", "),
          sum(ref_ids %in% runs[[j]]$`..row_id`)), character(1)),
          collapse = ", "),
        length(ref_ids), if (length(short) == 1L) "it is" else "they are")

    if (all(is.na(step_scores))) {
      .msg("select_features_forward(): every candidate failed to score at step ",
           length(selected) + 1L, "; stopping.")
      break
    }

    # The commonest way to call this wrongly is with a learner written for
    # cv_spatial(), which takes ONE argument: `function(train_sf, ...)` swallows
    # `vars` into the dots and fits the same model whatever it is handed.  The
    # training layer still carries every column, so nothing errors; every
    # candidate simply scores bit-identically to the null model, nothing
    # improves on it, and the function returns an empty selection and an NA
    # score with -- until this -- no word about why.  Two different predictor
    # sets do not produce the same cross-validated metric to the last digit,
    # so the signature is unambiguous.
    if (length(selected) == 0L && length(step_scores) >= 2L &&
        identical(baseline_type, "null model") &&
        is.finite(null_score) && all(is.finite(step_scores)) &&
        isTRUE(all.equal(unname(step_scores), rep(null_score, length(step_scores)),
                         tolerance = 1e-12)))
      .warn_and_log(paste0(
        "select_features_forward(): every candidate scored exactly what the ",
        "intercept-only model scored (%s = %.6g), so `fit_fn` appears to be ",
        "ignoring its second argument. It must be a function of ",
        "(train_sf, predictor_vars) that fits only the variables it is given; ",
        "a learner written for cv_spatial() takes one argument and needs a ",
        "wrapper."), metric, null_score)

    idx  <- if (metric == "R2") which.max(step_scores) else which.min(step_scores)
    cand <- remaining[idx]
    cand_score <- step_scores[[idx]]

    gain <- if (metric == "R2") cand_score - best else best - cand_score
    # The first variable is not tested against the fold mean; see above.
    provisional <- identical(baseline_type, "fold mean") && length(selected) == 0L
    if (is.finite(best) && !provisional && !accepts(gain, step_pg[[idx]])) {
      .msg(if (identical(rule, "se") && is.finite(gain) && gain > tol)
        sprintf(paste0(
          "select_features_forward(): the best candidate, %s, improves %s by ",
          "%.4g, but its mean gain per fold (%.4g) is within %g standard ",
          "error(s) of zero (SE %.4g over %d folds); stopping."),
          cand, metric, gain, step_pg[[idx]]$mean, se_mult, step_pg[[idx]]$se,
          step_pg[[idx]]$n)
      else
        sprintf("select_features_forward(): no candidate improves %s by more than %g; stopping.",
                metric, tol))
      break
    }

    selected  <- c(selected, cand)
    remaining <- setdiff(remaining, cand)
    best      <- cand_score
    best_pr   <- runs[[idx]]
    .msg(sprintf("select_features_forward(): + %s  (%s = %.4f)%s",
                 cand, metric, cand_score,
                 if (provisional) "  [provisional]" else ""))
  }

  # The finished path against the baseline.  For a learner that scored its
  # own null model every step has already beaten the one before, so this
  # only reports the total.  For the fold-mean baseline it is the test the
  # first step was spared: a path that does not beat the mean by the rule is
  # not a selection.
  path <- selected
  path_gain <- NA_real_
  path_pg   <- list(mean = NA_real_, se = NA_real_, n = 0L)
  kept <- NA
  if (length(path) > 0L && is.finite(null_score)) {
    path_gain <- if (metric == "R2") best - null_score else null_score - best
    path_pg   <- paired_gain(best_pr, base_pr)
    kept <- TRUE
    if (identical(baseline_type, "fold mean") && !accepts(path_gain, path_pg)) {
      kept <- FALSE
      .msg(sprintf(paste0(
        "select_features_forward(): the set reached, {%s}, does not beat the ",
        "training-fold mean by the rule (%s %.4f against %.4f; mean gain per ",
        "fold %.4g, SE %.4g), so nothing is selected."),
        paste(path, collapse = ", "), metric, best, null_score,
        path_pg$mean, path_pg$se))
      selected <- character(0)
    }
  }

  # The first variable of a path that started from the fold mean.  It was
  # added without a test, and the test above is of the whole set: two real
  # predictors and a candidate that happened to score best alone beat the
  # mean as clearly as the two real ones do.  At a true R2 of 0.33 a forest
  # sweep's kept selection began with an unrelated candidate in 3 runs of 20.
  # So a kept path of two or more is scored once more without its first
  # variable (one more set, k more fits), and that variable stays only if
  # the set is better with it by the rule every other variable was held to.
  # It is dropped only when the rest of the path still beats the fold mean
  # by the rule on its own; otherwise the set is kept whole and the variable
  # is reported as unproven.  A path of one needs no such step: the test
  # against the fold mean above was the test of that variable.
  first <- NULL
  if (identical(baseline_type, "fold mean") && isTRUE(kept) &&
      length(path) >= 2L) {
    v1         <- path[1L]
    rest       <- path[-1L]
    rest_pr    <- run_set(rest)
    rest_score <- score_on(rest_pr)
    first <- list(variable = v1, gain = NA_real_, fold_gain = NA_real_,
                  fold_gain_se = NA_real_,
                  rest = list(vars = rest, score = rest_score, gain = NA_real_,
                              fold_gain = NA_real_, fold_gain_se = NA_real_),
                  dropped = FALSE, note = NA_character_)
    brace <- function(v) paste0("{", paste(v, collapse = ", "), "}")
    if (!is.finite(rest_score)) {
      first$note <- "not tested: the set without it could not be scored"
      .warn_and_log(paste0(
        "select_features_forward(): %s was added first, without a test, and ",
        "the set without it, %s, could not be scored on the %d scored rows ",
        "(%d of them predicted), so it could not be tested at the end ",
        "either. It stays in the selection, unproven."),
        v1, brace(rest), length(ref_ids),
        sum(ref_ids %in% rest_pr$`..row_id`))
    } else {
      g1  <- if (metric == "R2") best - rest_score else rest_score - best
      pg1 <- paired_gain(best_pr, rest_pr)
      g0  <- if (metric == "R2") rest_score - null_score else null_score - rest_score
      pg0 <- paired_gain(rest_pr, base_pr)
      first$gain <- g1
      first$fold_gain <- pg1$mean; first$fold_gain_se <- pg1$se
      first$rest$gain <- g0
      first$rest$fold_gain <- pg0$mean; first$rest$fold_gain_se <- pg0$se
      if (accepts(g1, pg1)) {
        first$note <- "kept: the set is better with it by the rule"
        .msg(sprintf(paste0(
          "select_features_forward(): %s, added first without a test, is ",
          "kept: without it %s is %.4f against %.4f (mean gain per fold ",
          "%.4g, SE %.4g)."),
          v1, metric, rest_score, best, pg1$mean, pg1$se))
      } else if (accepts(g0, pg0)) {
        first$dropped <- TRUE
        first$note <- "dropped: the set is no better with it by the rule"
        .msg(sprintf(paste0(
          "select_features_forward(): %s was added first, without a test. ",
          "Without it, %s scores %s %.4f against %.4f with it (mean gain per ",
          "fold from %s %.4g, SE %.4g), so it is dropped."),
          v1, brace(rest), metric, rest_score, best, v1, pg1$mean, pg1$se))
        selected <- rest
        best     <- rest_score
        best_pr  <- rest_pr
      } else {
        first$note <- paste0("kept, unproven: the set is no better with it by ",
                             "the rule, and without it the set does not beat ",
                             "the baseline by the rule")
        .log_warn(paste0(
          "select_features_forward(): %s was added first, without a test, ",
          "and this selection does not show that it predicts. Without it, ",
          "%s scores %s %.4f against %.4f with it: a mean gain per fold of ",
          "%.4g with SE %.4g, which the rule does not count. It is kept ",
          "because %s does not beat the training-fold mean (%.4f) by the ",
          "rule on its own (mean gain per fold %.4g, SE %.4g). Read %s as ",
          "unproven."),
          v1, brace(rest), metric, rest_score, best, pg1$mean, pg1$se,
          brace(rest), null_score, pg0$mean, pg0$se, v1)
      }
    }
  }

  # The honest score: the selected set fitted on the selection half and
  # predicted on the estimation half, which the sweep never saw.
  score_holdout <- NA_real_
  if (!is.null(holdout_sf) && length(selected) > 0L) {
    score_holdout <- tryCatch({
      fit <- fit_fn(train_sf, selected)
      yh  <- as.numeric(do.call(stats::predict, list(object = fit, newdata = holdout_sf)))
      y   <- as.numeric(sf::st_drop_geometry(holdout_sf)[[response_var]])
      if (length(yh) != length(y))
        stop(sprintf("predict() returned %d values for %d rows", length(yh), length(y)))
      y_tr <- as.numeric(sf::st_drop_geometry(train_sf)[[response_var]])
      met <- .compute_reg_metrics(y, yh, p = NULL,
                                  y_train_mean = mean(y_tr[is.finite(y_tr)]))
      v <- met[[metric]]
      if (is.null(v) || !is.finite(v)) NA_real_ else as.numeric(v)
    }, error = function(e) {
      .warn_and_log("select_features_forward(): the hold-out score could not be computed (%s); score_holdout is NA.",
                    conditionMessage(e))
      NA_real_
    })
  }

  # Classed so that plot() finds plot.feature_selection(); the object is still
  # the same list, and `$`, `[[`, names() and printing all behave as before.
  structure(list(
    selected = selected,
    # Never hand back the `worst` sentinel as if it were a score: when nothing
    # was selected there is no score, and Inf / -Inf reads as a real number to
    # any caller that compares it.
    score    = if (length(selected) == 0L) NA_real_ else best,
    score_holdout = score_holdout,
    history  = if (length(history)) do.call(rbind, history) else
      data.frame(step = integer(0), variable = character(0), score = numeric(0),
                 n_pred = integer(0), gain = numeric(0), fold_gain = numeric(0),
                 fold_gain_se = numeric(0)),
    params   = list(metric = metric, method = method, k = k, tol = tol,
                    seed = seed, auto_range = isTRUE(auto_range),
                    select_on = select_on,
                    n_candidates = p, estimated_fits = est_fits,
                    n_scored = length(ref_ids),
                    rule = rule, se_mult = se_mult,
                    baseline = baseline_type),
    split    = split,
    # What step 1 had to beat, and what the finished path made of it.
    baseline = list(type = baseline_type, score = null_score,
                    path = path, gain = path_gain,
                    fold_gain = path_pg$mean, fold_gain_se = path_pg$se,
                    kept = kept,
                    # The end-of-path test of a first variable added without
                    # one; NULL when there was no such test to make.
                    first = first),
    # The metric of every set in every fold, for auditing the rule.
    fold_scores = if (length(fold_hist)) do.call(rbind, fold_hist) else
      data.frame(step = integer(0), variable = character(0),
                 fold = character(0), score = numeric(0))
  ), class = c("feature_selection", "list"))
}


#' The training-fold mean as a set of cross-validated predictions
#'
#' Built from candidate runs that \code{cv_spatial()} already made: each of
#' their prediction rows carries \code{y_train_mean}, the mean response of
#' the rows its fold trained on.  Used as \code{yhat} it is the prediction
#' of a model with no predictors, on the same rows and the same folds, at no
#' further fit.
#'
#' @param runs A list of prediction frames (or \code{NULL}s) from the sweep.
#' @return A prediction frame with one row per row any run predicted, or
#'   \code{NULL} when no run predicted anything.
#' @keywords internal
#' @noRd
.fs_fold_mean_predictions <- function(runs) {
  need <- c("..row_id", "fold", "y", "y_train_mean")
  ok <- Filter(function(pr) is.data.frame(pr) && nrow(pr) > 0L &&
                 all(need %in% names(pr)), runs)
  if (!length(ok)) return(NULL)
  bp <- do.call(rbind, lapply(ok, function(pr) pr[, need, drop = FALSE]))
  bp <- bp[!duplicated(bp[["..row_id"]]), , drop = FALSE]
  bp <- bp[is.finite(bp$y) & is.finite(bp$y_train_mean), , drop = FALSE]
  if (!nrow(bp)) return(NULL)
  bp$yhat <- bp$y_train_mean
  rownames(bp) <- NULL
  bp
}
