# Greedy forward feature selection with spatially blocked inner folds

Selects predictors by repeatedly adding whichever candidate most
improves a cross-validated score, stopping when no candidate improves it
by more than `tol`.

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
  select_on = c("all", "split")
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

  Minimum improvement required to accept a variable. Default 0, meaning
  any improvement is accepted. The first variable is judged against the
  null (intercept-only) model, so `tol` bites from step 1, but only when
  that null model can be scored. Backends that refuse a zero-length
  `predictor_vars`
  ([`fit_rf_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
  and
  [`fit_gwr_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  both do) have no null score, and there the first variable is accepted
  unconditionally.

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
  that half and scores it on the other: `score_holdout` is then an
  honest estimate of the selected model's `metric` on data the selection
  never saw (the sweep's own `score` is not; see "The score is not a
  performance estimate"). Both halves come back in `$split`. See the
  "Post-selection inference" section of
  [`determine_optimal_levels`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  for the trade: coverage for half the sample.

## Value

A list of class `"feature_selection"` (so that
[`plot.feature_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.feature_selection.md)
draws the selection path) with `selected` (the chosen predictors, in the
order they were added), `score`, `score_holdout`, `history`, `params`
and `split`. `score` is the winning set's cross-validated `metric` at
the final step: the **selection-internal** optimum, optimistically
biased because it was chosen as the best of many (see the section
above), and `NA` when nothing was selected. `history` is a data.frame
with `step`, `variable` and `score`, holding every candidate evaluated
at every step; when the null model could be scored it also carries a
`step = 0` row named `"<none>"` giving that baseline, so the first
variable's gain can be read off directly. `score_holdout` is `NA` unless
`select_on = "split"`, and then the selected set's `metric` when fitted
on the selection half and predicted on the estimation half (\\R^2\\
against the selection half's mean, the out-of-sample convention); `NA`
when nothing was selected or the prediction failed. `split` is `NULL` or
a list with `selection` and `estimation`, integer row positions in
`train_sf` after the completeness filter above. `params` records
`metric`, `method`, `k`, `tol`, `seed`, `auto_range`, `select_on`,
`n_candidates` and `estimated_fits`.

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
[`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md),
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)

## Examples

``` r
if (requireNamespace("GWmodel", quietly = TRUE) &&
    requireNamespace("sp", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  n <- 120
  pts <- st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
               a = rnorm(n), b = rnorm(n), noise = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
  pts$resp <- 3 * pts$a + 2 * pts$b + rnorm(n, 0, 0.3)

  fit_fn <- function(tr, vars) fit_gwr_model(tr, "resp", vars, bandwidth = 30)
  sel <- select_features_forward(pts, "resp", c("a", "b", "noise"), fit_fn,
                                 k = 3, quiet = TRUE)
  sel$selected
}
#> [1] "a" "b"
```
