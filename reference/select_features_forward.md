# Greedy forward feature selection with spatially blocked inner folds

Selects predictors by repeatedly adding whichever candidate most
improves a cross-validated score, stopping when the best candidate's
improvement is no larger than the fold-to-fold noise in it
(`rule = "se"`, the default) or no larger than `tol` (`rule = "tol"`).
The first variable has a baseline to beat with every learner, so the
selection can be empty.

## Usage

``` r
select_features_forward(
  train_sf,
  response_var,
  candidate_vars,
  fit_fn,
  k = 5,
  method = c("block_kfold", "random_kfold"),
  block_size = NULL,
  metric = c("RMSE", "MAE", "R2"),
  tol = 0,
  max_vars = NULL,
  max_fits = 5000L,
  seed = 123,
  quiet = FALSE,
  auto_range = FALSE,
  select_on = c("all", "split"),
  rule = c("se", "tol"),
  se_mult = 1,
  null_fn = NULL
)
```

## Arguments

- train_sf:

  Training data (`sf`).

- response_var:

  Character(1).

- candidate_vars:

  Character vector of predictors to choose among.

- fit_fn:

  A function `(train_sf, predictor_vars)` returning a `spatial_fit`. The
  signature takes two arguments because selection has to refit with
  different predictor sets.

- k:

  Inner fold count. Default 5.

- method:

  Inner fold method. Default `"block_kfold"`.

- block_size:

  Passed to
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md);
  inherit the outer block size so inner and outer blocks are on the same
  spatial scale.

- metric:

  Score to optimise: `"RMSE"`, `"MAE"` (minimised) or `"R2"`
  (maximised). Default `"RMSE"`.

- tol:

  Minimum improvement of the pooled score required to accept a variable,
  in the units of `metric`. Default 0. It is a floor under both rules:
  with `rule = "tol"` it is the whole test, and with `rule = "se"` the
  per-fold test applies as well. The first variable is judged against a
  baseline with every learner; see *The baseline for the first
  variable*. Any single number is accepted: `Inf` leaves nothing
  selected, and `-Inf` puts no floor under a step.

- max_vars:

  Optional cap on how many predictors to select.

- max_fits:

  Abort if the sweep would exceed this many model fits. Default 5000.

- seed:

  RNG seed. It governs both the inner fold construction and the
  cross-validation itself: it is forwarded to `cv_spatial(seed = )`,
  which draws one RNG stream per fold from it, so it also seeds the
  *learner* inside every fold. A stochastic `fit_fn` is therefore
  reproducible from this one value.

- quiet:

  Logical; suppress this function's progress
  [`message()`](https://rdrr.io/r/base/message.html)s. It does not
  silence R warnings, nor the package's console log echo (see
  [`spatialkit_quiet`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)
  for that). Default `FALSE`.

- auto_range:

  Logical. If `TRUE` and `method` is `"block_kfold"`, the
  autocorrelation range is estimated from the response (detrended on the
  candidates) and used as the minimum block size of the inner folds,
  exactly as in
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md);
  `block_size` still applies as a floor. Default `FALSE`, which keeps
  the geometric blocks. Inner blocks smaller than the range let a
  candidate be selected for spatial proximity to the response rather
  than for predicting it, which is the same failure the `random_kfold`
  caution above exists to prevent, so the leakage warning
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  raises applies here with more force than usual.

- select_on:

  `"all"` (default) runs the sweep on every row of `train_sf`. `"split"`
  runs it on one spatially blocked half, then fits the selected set on
  that half and scores it on the other: `score_holdout` is then the
  selected model's `metric` on rows whose response the sweep never read
  (the sweep's own `score` is not; see "The score is not a performance
  estimate"). It is one estimate from one region: the halves share a
  border with no buffer, so rows near it are still correlated with the
  selection half. Both halves come back in `$split`. See the
  "Post-selection inference" section of
  [`determine_optimal_levels`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  for the trade: coverage for half the sample.

- rule:

  How a step is judged. `"se"` (default) accepts the best candidate when
  its pooled gain exceeds `tol` and its mean per-fold gain exceeds
  `se_mult` standard errors of that mean; `"tol"` accepts any pooled
  gain above `tol`, the behaviour before this argument existed. See
  *When a variable is accepted*.

- se_mult:

  Number of standard errors the mean per-fold gain must exceed under
  `rule = "se"`. Default 1; 0 asks only that the mean per-fold gain be
  positive.

- null_fn:

  Optional function `(train_sf)` returning a fitted model with no
  predictors, used as the baseline for the first variable when `fit_fn`
  cannot fit an empty predictor set. Default `NULL`:
  `fit_fn(train_sf, character(0))` is tried, and the training-fold mean
  is used when that fails. The example builds one.

## Value

A list of class `"feature_selection"` (so that
[`plot.feature_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.feature_selection.md)
draws the selection path) with `selected` (the chosen predictors, in the
order they were added; possibly none), `score`, `score_holdout`,
`history`, `params`, `split`, `baseline` and `fold_scores`. `score` is
the winning set's cross-validated `metric` at the final step: the
**selection-internal** optimum, optimistically biased because it was
chosen as the best of many (see the section above), and `NA` when
nothing was selected. `history` is a data.frame with `step`, `variable`,
`score`, `n_pred`, `gain`, `fold_gain` and `fold_gain_se`, holding every
candidate evaluated at every step, and a `step = 0` row named `"<none>"`
giving the baseline, so the first variable's gain can be read off
directly. `gain` is the candidate set's improvement in the pooled score
over the set it was added to (the baseline at step 1), `fold_gain` the
mean of the per-fold improvements and `fold_gain_se` its standard error:
the three numbers the rule reads. `fold_scores` is a data.frame with
`step`, `variable`, `fold` and `score`: the metric of every set in every
fold. `baseline` is a list with `type` (`"null model"`, `"null_fn"` or
`"fold mean"`), `score`, `path` (the set the sweep reached), `gain`,
`fold_gain` and `fold_gain_se` (that set against the baseline), `kept`
(`FALSE` when the path did not beat a fold-mean baseline and `selected`
is empty for that reason; `NA` when the path is empty) and `first`.
`first` is `NULL` unless the baseline is the fold mean and a path of two
or more variables was kept; it is then the test of the first variable, a
list with `variable`; `gain`, `fold_gain` and `fold_gain_se` (the path
against the path without that variable); `rest`, a list with `vars`,
`score`, `gain`, `fold_gain` and `fold_gain_se` (the path without it,
and that set against the baseline); `dropped` (`TRUE` when the variable
was taken out: `selected` is then `path[-1]` and `score` is
`rest$score`); and `note`, the verdict in words. The set scored for that
test is not a step of the sweep, so it is not in `history` or
`fold_scores`. Every set is scored on the same rows: those the null
model's cross-validation predicted or, when the baseline is the fold
mean, those any step-1 set predicted; `params$n_scored` counts them (a
warning says so when that is fewer than all). `n_pred` is how many rows
the set's cross-validation predicted. A set that left some of the scored
rows unpredicted, because a fold failed for it, has `score` `NA`, with a
warning naming it: scored on the rows it did predict it would be
compared on fewer, usually easier, rows than its rivals. A factor with a
level found in one spatial block only is the usual case, and cannot be
selected. `score_holdout` is `NA` unless `select_on = "split"`, and then
the selected set's `metric` when fitted on the selection half and
predicted on the estimation half (\\R^2\\ against the selection half's
mean, the out-of-sample convention); `NA` when nothing was selected or
the prediction failed. `split` is `NULL` or a list with `selection` and
`estimation`, integer row positions in `train_sf` as passed; rows the
completeness filter above dropped are in neither. `params` records
`metric`, `method`, `k`, `tol`, `seed`, `auto_range`, `select_on`,
`n_candidates`, `estimated_fits`, `n_scored`, `rule`, `se_mult` and
`baseline` (its type).

## Details

**The inner folds must be spatial, and that is the entire point.**
Nested selection is only worth doing if the inner loop is blocked the
same way the outer one is. Random inner folds inside blocked outer folds
select variables that look predictive only because nearby points leak
between train and test. The outer loop then reports honest-looking
numbers for a dishonestly chosen feature set, which is worse than not
selecting at all, because the dishonesty is now hidden behind a
defensible-looking validation. `method` therefore defaults to
`"block_kfold"` and logs a loud caution if set to `"random_kfold"` (a
deliberate choice, so it is not raised as an R warning).

Call this *inside* the `fit_fn` you pass to
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md).
`.cv_fit_one_fold()` calls `fit_fn(train_sf)` on the training slice
only, so anything done inside it is automatically nested and leak-free;
no extra plumbing is needed. The cost grows fast: a sweep over `p`
candidates costs roughly `p^2 / 2 * k` model fits, and nesting that
inside `n` outer leave-one-out folds multiplies it by `n`. `max_fits`
guards against that.

## When a variable is accepted

Each candidate set is scored by cross-validation on the same inner
folds, and the candidate with the best pooled score is the one
considered. Under `rule = "se"` it is accepted when two things hold: the
pooled score improves by more than `tol`, and the mean of the per-fold
improvements (the metric computed in each fold for the set with and
without the candidate, and differenced) exceeds `se_mult` times the
standard error of that mean across the folds. This is the
one-standard-error idea of Breiman et al. (1984) applied to a step: a
gain the folds do not agree on is not evidence that the variable
predicts. The fold is the unit because under spatially blocked folds the
errors within a block are correlated, so a standard error over rows
would understate the noise. With `tol = 0` and any gain accepted (the
previous default, still available as `rule = "tol"`), a GWR sweep over
one real predictor and six noise candidates added a noise variable in 11
runs of 20; the default rule added one in 4 of 20 and kept the real
predictor in all 20. The rule costs power on weak predictors: a second
real predictor with a coefficient of 0.3 against 1 was kept in 14 runs
of 20, against 18 of 20 when any gain was accepted. With few folds the
standard error is itself noisy, so raise `k` before lowering `se_mult`.
`history` records each candidate's pooled gain, mean per-fold gain and
its standard error, and `fold_scores` the per-fold scores themselves, so
every decision can be checked.

## The baseline for the first variable

The first variable is compared with a model that has no predictors,
which is recorded as step 0 of `history`. Which model that is depends on
the learner, and `baseline$type` says which was used:

- `"null model"`:

  `fit_fn(train_sf, character(0))`, when the learner can fit it. Step 1
  is tested against it like any other step.

- `"null_fn"`:

  the model `null_fn(train_sf)` returns, for a learner that cannot fit
  an empty set but has a natural no-predictor model. Step 1 is tested
  against it.

- `"fold mean"`:

  the mean response of each fold's training rows, when neither is
  available;
  [`fit_rf_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
  refuses an empty predictor set. Here the first variable is
  *provisional*: it is added without a test, the sweep continues by the
  rule, and the finished set is kept only if it beats the fold mean by
  the same rule; otherwise nothing is selected, and `baseline$path`
  holds the set that was reached. The test is at the end of the path,
  not at step 1, because a forest on a single variable can lose to the
  mean when two variables together beat it clearly (at a true \\R^2\\ of
  0.33 the best single-variable forest had RMSE 1.276 against 1.206 for
  the mean and 1.101 for the pair), so a test at step 1 would select
  nothing there. The first variable of a kept path then gets a test of
  its own; see *The first variable under a fold-mean baseline*.

Without a baseline, as before, at least one variable was always
selected: with a forest, in 10 draws of 10 on a response unrelated to
every candidate.

## The first variable under a fold-mean baseline

The test at the end of the path is of the whole set, so it does not show
that the first variable belongs in it: two real predictors and an
unrelated candidate that happened to score best alone beat the mean as
the two real ones do. A kept path of two or more variables is therefore
scored once more without its first variable (one more set, `k` more
fits). That variable is dropped from `selected` when the set is no
better with it by the rule and the rest of the path still beats the fold
mean by the rule. `baseline$first` records the test and `baseline$path`
keeps the set as the sweep reached it, so `selected` is
`baseline$path[-1]` after a drop. A path of one variable needs no second
test: the test against the fold mean was the test of that variable.

This test is the one every later variable passed, an improvement of the
set by the rule, and it has the same limit. In 20 forest sweeps at a
true \\R^2\\ of 0.33, with two real predictors and three unrelated
candidates, the path began with an unrelated candidate in 3: it was
dropped in 2 and kept in 1, where the forest scored better with it than
without (*Selection with a random forest* says why it can). At a true
\\R^2\\ of 0.15, 9 selections of 20 were kept, 3 of them began with an
unrelated candidate, and all 3 kept it for that reason. So in a
selection of two or more variables from a learner that falls back to the
fold mean, `selected[1]` has been shown to improve the set, not to
predict on its own. The test also costs power on a weak first predictor,
as the rule does at every step: with one real predictor of coefficient
0.3 and five unrelated candidates, a linear learner that refuses the
empty set kept the real predictor in 77 sweeps of 100, against 80 before
this test existed and 77 for the same learner fitting its own null model
(in the 3 sweeps that differ it was the first variable, its gain was
within one standard error, and the learner with a null model selected
nothing). When the set is no better with the first variable but the rest
of the path does not beat the fold mean by the rule either, the two
cannot be told apart: the set is kept whole, a caution is logged, and
`baseline$first$note` says the variable is unproven. A learner that fits
its own null model gets no such test: its first variable passed one at
step 1.

## The baseline for a GWR

For a GWR the fold mean is a weak baseline, because every GWR carries a
local intercept and beats a global mean on any smooth response whatever
its predictors are: on a smooth response unrelated to six candidates, a
GWR sweep against the fold mean selected one of them in 20 runs of 20.
Each of those paths had one variable, which the test of the first
variable above does not reach. The comparison that matters there is with
the kernel-weighted local mean at the same bandwidth, a GWR with a local
intercept and no predictors: against it the same sweeps selected a
variable in 4 runs of 20.
[`fit_gwr_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
fits that model when its `predictor_vars` is `character(0)`, so a
`fit_fn` that hands its variables to
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
gets it as the null model, and the first variable is tested against it
at step 1 (`baseline$type` is `"null model"`). Give such a `fit_fn` a
fixed `bandwidth`, so that every set is fitted under the kernel the
baseline has; with `bandwidth = NULL` each set, the empty one included,
searches for its own. A `fit_fn` that refuses an empty set falls back to
the fold mean, and the sweep logs a caution: read a GWR selection of one
variable made that way as unproven, or supply the local mean as
`null_fn`, as the second half of the example at the end of this page
does for a learner written by hand.

[`predict()`](https://rdrr.io/r/stats/predict.html) on a GWR warns about
held-out locations it cannot estimate and about those beyond the reach
of a fold's fit (classes `"spatialkit_gwr_not_estimable"` and
`"spatialkit_gwr_extrapolation"`; see
[`predict.gwr_fit`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.gwr_fit.md)).
The sweep raises each once, under this function's name, with the totals
over the candidate sets it scored; every set is predicted on the same
rows, so a location can be counted once per set.

## Selection with a random forest

A forest confounds selection with tuning. Adding a predictor that
carries no signal still changes how the trees split, and with ranger's
default node size that regularises the forest: in one test the
cross-validated RMSE was 1.458 with the real predictor alone, 1.340 with
a noise variable added, and 1.341 with the real predictor alone and
`min_node_size = 20`. The noise variable bought what tuning the node
size buys. That gain is real, so no acceptance rule refuses it: with one
real predictor and six noise candidates a forest sweep kept the real one
in 20 runs of 20 and added 1.35 noise variables per run under the
default rule (2.4 per run, over 10 runs, when any gain was accepted).
Setting the node size in `fit_fn` reduces this without removing it (0.8
per run at `min_node_size = 20`). So a variable in a forest selection
may be there as a tuning effect wherever in the order it was added:
later additions with small gains are the usual case, and the first
variable is another, because a forest sweep starts from the fold mean
and its first variable is tested as an addition to the others (see *The
first variable under a fold-mean baseline*). Compare the selected set
with a tuned forest on the variables you have other reasons to believe
in:


    fit_fn <- function(train_sf, vars)
      fit_rf_model(train_sf, "z", vars, min_node_size = 20)

On a response unrelated to every candidate the same sweeps selected
nothing in 20 runs of 20 (19 of 20 with the larger node size).

## The score is not a performance estimate

`$score` is the cross-validated `metric` of the winning set at the final
step: the best of every candidate set the sweep scored. That is the
number the selection optimised, and a number optimised over many
candidates is optimistically biased by construction: Cawley and Talbot
(2010) show the bias can exceed the genuine differences between the
models being compared. Quote it as the selection criterion, not as the
performance of the selected model. An honest performance estimate needs
data the selection never saw. Two ways to get one: run this function
inside the `fit_fn` of
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
so the outer folds score a model whose predictors were chosen on the
inner ones alone; or pass `select_on = "split"`, which selects on one
spatial half of `train_sf` and returns the selected set's score on the
other as `score_holdout`.

## References

Breiman, L., Friedman, J. H., Olshen, R. A. and Stone, C. J. (1984).
*Classification and Regression Trees*. Wadsworth.

Cawley, G. C. and Talbot, N. L. C. (2010). On over-fitting in model
selection and subsequent selection bias in performance evaluation.
*Journal of Machine Learning Research*, 11, 2079–2107.
<https://jmlr.org/papers/v11/cawley10a.html>

## See also

Other cross-validation:
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md),
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
[`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md),
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md),
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
[`fold_separation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fold_separation.md),
[`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md),
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md),
[`sac_nugget()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md)

## Examples

``` r
if (requireNamespace("sp", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  n <- 120
  pts <- st_as_sf(
    data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
               a = rnorm(n), b = rnorm(n), noise = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
  pts$resp <- 3 * pts$a + 2 * pts$b + rnorm(n, 0, 0.3)

  # The learner: a GWR with an adaptive bisquare kernel reaching the 30th
  # nearest training point.  With no variables it is the kernel-weighted
  # local mean, which is the baseline: a GWR's local intercept beats a
  # global mean on any smooth response, so the first variable has to beat
  # the local mean under the same kernel.
  fit_fn <- function(tr, vars) fit_gwr_model(tr, "resp", vars, bandwidth = 30)
  sel <- select_features_forward(pts, "resp", c("a", "b", "noise"), fit_fn,
                                 k = 3, quiet = TRUE)
  sel$selected
  sel$baseline$type        # "null model": step 1 was tested against it
  sel$history[sel$history$step <= 1, c("step", "variable", "score", "gain")]

  # A learner that cannot fit an empty set takes its no-predictor model
  # through `null_fn`.  Here the same local mean, written as a spatial_fit
  # of its own with a predict() method, for a learner that insists on a
  # variable.  Without `null_fn` that sweep would start from the mean of
  # each training fold.
  needs_one <- function(tr, vars) {
    if (!length(vars)) stop("this learner needs a predictor")
    fit_gwr_model(tr, "resp", vars, bandwidth = 30)
  }
  local_mean <- function(tr, bandwidth = 30) {
    new_spatial_fit("local_mean_fit",
                    engine = list(xy = st_coordinates(tr), y = tr$resp,
                                  bandwidth = bandwidth),
                    formula = resp ~ 1, response_var = "resp",
                    predictor_vars = character(0), data_sf = tr)
  }
  predict.local_mean_fit <- function(object, newdata = NULL, ...) {
    e  <- object$engine
    at <- if (is.null(newdata)) e$xy else st_coordinates(newdata)
    vapply(seq_len(nrow(at)), function(i) {
      d <- sqrt((e$xy[, 1] - at[i, 1])^2 + (e$xy[, 2] - at[i, 2])^2)
      h <- sort(d)[min(e$bandwidth, length(d))]  # the kernel's reach here
      w <- ifelse(d > h, 0, (1 - d^2 / h^2)^2)   # bisquare weights
      sum(w * e$y) / sum(w)
    }, numeric(1))
  }
  registerS3method("predict", "local_mean_fit", predict.local_mean_fit)
  sel2 <- select_features_forward(pts, "resp", c("a", "b", "noise"), needs_one,
                                  k = 3, quiet = TRUE, null_fn = local_mean)
  sel2$baseline$type       # "null_fn"
  all.equal(sel2$baseline$score, sel$baseline$score)  # the same baseline
}
#> [1] TRUE
```
