# Cross-validate any learner with spatial folds

You supply `fit_fn(train_sf)`, which fits on one fold's training rows
and returns a `spatial_fit`; `cv_spatial()` builds or takes the folds,
calls `fit_fn` once per fold, predicts the held-out rows and returns
pooled and per-fold scores.
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
is a wrapper over it.

## Usage

``` r
cv_spatial(
  data_sf,
  response_var,
  predictor_vars,
  fit_fn,
  folds = NULL,
  k = 5,
  seed = 123,
  boundary = NULL,
  pointize = "auto",
  predict_args = list(),
  fold_info_fn = NULL,
  p = NULL,
  block_size = NULL,
  auto_range = FALSE,
  parallel = FALSE,
  metrics = NULL,
  .caller = "cv_spatial"
)
```

## Arguments

- data_sf:

  An sf layer with the response and the predictors.

- response_var:

  Name of the response column.

- predictor_vars:

  Character vector of predictor column names.

- fit_fn:

  A function of one argument, the training slice of `data_sf`, returning
  a `spatial_fit` built with
  [`new_spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md).
  It is called once per fold on the training rows only, so anything done
  inside it (scaling, tuning, an inner variable sweep) is already nested
  and leak-free. The `subclass` it stamps must have a
  `predict.<subclass>()` method registered, since that is how each fold
  is scored. The fit models `response_var` (its `response_var` is the
  call's): a fold whose fit models another column fails, with the two
  names in `fold_status`, since its predictions are scored against
  `response_var`.

- folds:

  Optional fold definitions, in one of three shapes: a
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  result; a bare list of `list(train =, test =)` pairs of `..row_id`
  values; or a vector of fold labels, one per row, which becomes
  leave-that-label-out splits. The label vector is how folds built by
  another package are used here: `blockCV::cv_spatial()` returns one as
  `$folds_ids`, while its `$folds_list` holds two unnamed vectors per
  fold and is refused by name. Train and test must be disjoint; a fold
  that trains on its own test rows is not a cross-validation split and
  is refused with an error. IDs naming no row of the prepared data are
  dropped with a logged count (expected when rows were removed for
  missing values; a sign the folds came from other data when they were
  not). When `NULL`, the folds are `make_folds(method = "block_kfold")`,
  built with `k`, `seed`, `boundary`, `block_size` and `auto_range`.
  With `folds = NULL` the score estimates the error at the distance the
  block grid built here puts between the held-out points and their
  training data, which overstates the error of a map where the samples
  cover the area it predicts over. The excess is how much the model's
  error grows from the map's distances to the data to the block grid's,
  and is small for a model whose error changes little with the distance
  from the data.
  [`fold_separation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fold_separation.md)
  with `prediction_points` compares the held-out distances with those of
  the locations the model will predict at, and
  `make_folds(method = "knndm", prediction_points = )` builds folds
  matched to them.

- k:

  Number of folds to build when `folds` is `NULL`; not read when `folds`
  is supplied. Default 5.

- seed:

  Seed for building the folds when `folds` is `NULL`, and for the one
  seed drawn per fold before `fit_fn` runs, so a learner that draws
  random numbers gives the same scores on every run, sequential or
  parallel. With a seed, the caller's random stream is restored
  afterwards. Default 123; `NULL` draws the per-fold seeds from the
  session's RNG state, and advances it.

- boundary:

  Optional polygon layer (sf or sfc). Data without a CRS are aligned to
  it, and the blocks built when `folds` is `NULL` are clipped to it (see
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  and
  [spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)).

- pointize:

  How a feature that is not a point is reduced to one before fitting: a
  mode of
  [`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md).
  Default `"auto"`.

- predict_args:

  A list of extra arguments for
  [`predict()`](https://rdrr.io/r/stats/predict.html) on each fold's
  fit.

- fold_info_fn:

  Optional `function(fit, test_sf, y, yhat)` returning a named list (or
  a named vector) of per-fold extras (a bandwidth, a tuning value,
  anything read off the fitted object), added as columns of
  `fold_metrics`. It sees the fit and the held-out layer, which
  `metrics` does not; it is applied per fold only, and its values are
  not pooled. An element `..per_row` that is a data frame with one row
  per held-out observation is spliced into `predictions` instead (`NA`
  in the rows of a fold that returned none; one of the wrong length is
  dropped and logged); its columns must be named, once, and not reuse a
  column `predictions` already has (`..row_id`, `fold`, `y`, `yhat`,
  `y_train_mean`). Every other element must be named, hold one value,
  and not reuse a column `fold_metrics` already has; anything else is an
  error. A `fold_info_fn` that throws on a fold is logged, its columns
  are `NA` for that fold, and `fold_status$message` says so; the fold is
  kept.

- p:

  Number of predictors for the per-fold adjusted \\R^2\\, or `NULL` to
  skip it. It applies only to models with one global parameter count;
  pass `NULL` for models with spatially varying coefficients, such as
  GWR.

- block_size:

  Optional minimum block edge length, in the units of the CRS the folds
  are built in, for the blocks built when `folds` is `NULL`; see
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md).
  It is a size you choose and is compared with the autocorrelation range
  only when one can be estimated, in which case a smaller size is warned
  about. Beside supplied `folds` it has nothing to size and raises a
  warning of class `spatialkit_fold_args_ignored`.

- auto_range:

  Logical, default `FALSE`. `TRUE` with `folds = NULL` estimates the
  autocorrelation range and uses it as the minimum block size. `TRUE`
  beside supplied `folds` raises the same warning as `block_size`.

- parallel:

  Logical or positive integer. `TRUE` detects the number of cores and
  fits the folds in parallel with
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  (macOS and Linux; sequential on Windows); an integer above 1 uses that
  many cores. Default `FALSE` (sequential). A learner that runs OpenMP
  code (a function of GWmodel, or an xgboost built with GNU libgomp) can
  hang the forked workers once it has run in the session; keep such a
  `fit_fn` sequential.

- metrics:

  Optional scoring function of your own; see the section "Your own
  metrics". Default `NULL`: the built-in metrics only.

- .caller:

  For package wrappers; do not set it.

## Value

A list with:

- `overall`:

  One row of metrics pooled over the held-out rows of every fold that
  produced predictions.

- `fold_metrics`:

  One row per fold that produced predictions: `fold`, `n_train`,
  `n_test`, `n_pred`, the metrics, and the columns of `fold_info_fn` and
  `metrics`.

- `predictions`:

  One row per held-out row scored: `..row_id`, `fold`, `y`, `yhat` and
  `y_train_mean`, the mean response of that fold's training rows.

- `folds`:

  The splits as run, each with its `fold_id`.

- `n_folds_attempted`, `n_folds_succeeded`:

  The folds supplied (or built) and the folds scored. A `fit_fn` that
  fails on every fold otherwise looks like a run that scored `NA`, so
  compare the two before reading `overall`.

- `fold_status`:

  A data frame with one row per fold supplied: `fold`, `status` and
  `message`. `status` is one of

  - `"ok"`: `message` is empty unless `fold_info_fn` failed on the fold;

  - `"error"`: the fit or its
    [`predict()`](https://rdrr.io/r/stats/predict.html) threw, and
    `message` is the error text;

  - `"skipped"`: nothing to score (too few matched rows, a prediction of
    the wrong length, or no finite pair of observed and predicted
    values);

  - `"dropped"`: an empty test set, or fewer than two training rows once
    incomplete rows were removed, so the fold never reached the fitter;

  - `"worker_error"`: a parallel worker died, for example killed for
    lack of memory.

  Every fold missing from `fold_metrics` has its reason here, after the
  console output of a long run is gone.

- `orphan_rows`:

  The `..row_id`s of rows in the data that no fold names: they enter no
  training set and are never scored. It is empty unless the folds were
  built on a different or subsetted layer.

- `n_unknown_ids`:

  The number of distinct row IDs the folds name that the data do not
  have (expected when rows were removed for missing values). Each such
  row is named by every fold, once as a test row and once in each other
  fold's training set, and this counts the row, not the mentions.

- `n_dropped`:

  The number of rows
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
  removed for missing or non-finite values or a bad geometry before any
  fold was fitted.

A failure costs that fold only: with `parallel = TRUE` each fold runs in
its own worker, and a sequential run catches the error per fold. An
error that stops a sequential run (a `metrics` or `fold_info_fn` return
value of the wrong shape) stops a parallel one too, naming the fold.
When some folds, but not all, end as `"error"`, `"skipped"` or
`"worker_error"`, the function warns, naming them and how many rows
`overall` covers: it is pooled over the folds that produced predictions,
and the fold that fails is often the hardest to predict, so it may
flatter the model. A `"dropped"` fold has its own warning.

The `fold` column of `fold_metrics`, `predictions` and `fold_status`
carries the fold's index in the `folds` object that was supplied, so it
lines up with `make_folds()$assignment$fold` even when some folds were
unusable and dropped. Splits that already carry a `fold_id`, as the
`folds` of a `cv_*()` result do, keep it: handing one run's `folds` to
another labels each fold as the first run and
[`fold_separation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fold_separation.md)
do, a dropped fold's gap included.

`R2` is the out-of-sample \\R^2\\: the total sum of squares is taken
about the mean of the training rows (in `overall`, each held-out row
about its own fold's training mean), the null prediction available when
the fold is predicted, and not about the held-out rows' own mean, a null
model that would know the test data. On spatial blocks the two can
differ widely; `R2` is below 0 when the model predicts worse than the
training mean. `model_metrics(newdata = )` uses the same baseline.
`overall$Adj_R2` is always `NA`: the pooled out-of-sample predictions
come from `k` separately fitted models and have no single parameter
count to adjust for. The per-fold `fold_metrics$Adj_R2` carries the
adjusted value when `p` is supplied, and is `NA` otherwise. A `fit_fn`
that returns a `gwr_fit` of a response at a large offset gives
`predictions` the columns `y_centre`, `yhat_centred` and
`y_train_mean_centred` (see
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md)).

## Name collision with blockCV

The blockCV package exports a function of the same name that does the
opposite job: `blockCV::cv_spatial()` builds spatial folds, where this
function runs a cross-validation over folds it is given. With both
packages attached, whichever was attached last masks the other;
`spatialkit::cv_spatial()` always resolves to this one. The two work
together: `blockCV::cv_spatial()` returns its fold assignment as
`$folds_ids`, a vector of fold labels, and that vector is accepted as
the `folds` argument here and in every other `cv_*()` function.

## Your own metrics

The built-in columns (`RMSE`, `MAE`, `MAPE`, `SMAPE`, `R2`, `Adj_R2`)
are the Gaussian regression set; see "Which metrics survive a
non-Gaussian response" in
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
for which of them hold up on a count or a bounded response. `metrics`
scores what they cannot: a `function(y, yhat)` that returns a named
numeric vector (a named list of scalars, or a one-row data frame, also
serve), for example a Poisson deviance, a log score on a probability, or
a loss with your own weights. It is applied twice, as the built-in
metrics are: once per fold, to that fold's held-out rows, so each name
becomes a column of `fold_metrics`; and once to the pooled out-of-sample
predictions of every fold, so each name becomes a column of `overall`.
Only the pairs the built-in metrics use reach the function (both `y` and
`yhat` finite), so its columns describe the same rows as `RMSE`, and
`n_pred` counts them.

The contract: every element named, names unique and not one of the
built-in column names, one number per name. The built-in names include
the per-fold extras of the backend or of your `fold_info_fn`
(`bandwidth` for
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md);
`CRPS`, `coverage_*`, `gp_k`, `gp_n_basis` and `n_draws` for
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)),
and `mean_CRPS`, which
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
writes. Anything else is an error, since a scoring function that returns
the wrong shape is a mistake to report, not a fold to skip. A function
that throws on a fold is logged and its columns are `NA` for that fold
(and for `overall`, if it throws on the pooled predictions); a fold is
never dropped for it. When no fold produced a prediction the empty
`fold_metrics` frame still carries the function's columns, typed,
provided the function can be called on zero-length input.

`fold_info_fn` is the per-fold half of the same mechanism, with access
to the fitted object and the held-out layer; `metrics` sees only the two
vectors but is also pooled.
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
hands one `metrics` to every backend, so the columns are comparable
across the rows of its `overall`.

## Warnings a learner raises in every fold

A warning from a fold's fit or from its
[`predict()`](https://rdrr.io/r/stats/predict.html) reaches the caller
as it is raised, so a learner that warns in every fold warns once per
fold. Two groups of the package's own warnings are collected over the
folds instead and raised once per call, under this function's name and
with the class they had:

- the two in which [`predict()`](https://rdrr.io/r/stats/predict.html)
  on a GWR counts locations (`"spatialkit_gwr_not_estimable"` and
  `"spatialkit_gwr_extrapolation"`), with the totals over the folds that
  had any, as
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md)
  raises them (see
  [`predict.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-predict.md));

- the four
  [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
  raises about a fit, naming the folds, as
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
  raises them (see the section "Warnings from the folds" there).

With `parallel` above one core a fold's warnings are raised after the
folds have run, once per distinct text; one of the package's own classes
keeps its class and its fields.

## Percentage errors on responses with zeros

`MAPE` divides by the observed value and `SMAPE` by \\\|y\| +
\|\hat{y}\|\\, so neither is defined where its denominator is zero.
Neither returns `Inf` or `NaN`. Both are averaged over the rows whose
denominator is non-zero, and are `NA` when no row qualifies. Non-zero is
judged at the scale of the data: a denominator no larger than 100
machine epsilons times the largest one counts as zero, so the rule does
not depend on the units of the response. The `n_MAPE` and `n_SMAPE`
columns record how many rows that was; the `n` column (`n_pred` in a
cross-validation result) counts finite observation/prediction pairs.
Read a percentage error next to that count: when `n_MAPE` is smaller,
`MAPE` is an average over a subset of the data, whatever its value.

Any response taking exact zeros is affected: counts, rainfall,
abundance, claim amounts. On a zero-inflated response with 62 zeros out
of 120, `MAPE` averages the 58 non-zero rows and `n_MAPE` is 58. `SMAPE`
drops the rows where observation and prediction are both near zero,
which on a well-fitted zero-inflated model are the rows it got right, so
it averages the harder rows only and reads worse than the fit deserves.
`n_SMAPE` shows how many rows it kept; it does not correct the average.

`RMSE`, `MAE` and \\R^2\\ use every finite row; prefer them whenever the
response can be zero. For a Bayesian GP fit,
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
also reports CRPS and interval coverage, which are proper scoring rules
and are defined at zero.

## See also

[`new_spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md)
for the constructor a `fit_fn` must use;
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
and
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
for the built-in backends.

Other cross-validation:
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md),
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md),
[`plot_calibration()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_calibration.md),
[`plot_cv_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_cv_metrics.md)

## Examples

``` r
library(sf)
set.seed(1)
n <- 80
site <- st_as_sf(
  data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
             elev = rnorm(n)),
  coords = c("x", "y"), crs = 32632
)
site$price <- 10 + 0.01 * (st_coordinates(site)[, 1] - 5e5) +
  2 * site$elev + rnorm(n)

# 1. A fit_fn returning a spatial_fit of your own subclass.
lm_fit <- function(train_sf) {
  new_spatial_fit(
    subclass       = "lm_fit",
    engine         = lm(price ~ elev, st_drop_geometry(train_sf)),
    formula        = price ~ elev,
    response_var   = "price",
    predictor_vars = "elev",
    data_sf        = train_sf
  )
}

# 2. The predict() method cv_spatial() scores each fold with. Without it
#    every fold fails and `overall` comes back all-NA.
predict.lm_fit <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) newdata <- object$data_sf
  as.numeric(stats::predict(object$engine, st_drop_geometry(newdata)))
}
registerS3method("predict", "lm_fit", predict.lm_fit)

cv <- cv_spatial(site, "price", "elev", fit_fn = lm_fit, k = 3, seed = 1)
#> cv_spatial(): no folds supplied -- using spatial block k-fold CV (k=3).
cv$overall
#>       RMSE      MAE     MAPE    SMAPE        R2 Adj_R2 n_pred n_MAPE n_SMAPE
#> 1 3.472328 3.009688 21.29953 20.27029 0.2341257     NA     80     80      80
# Compare these before reading the metrics above; fold_status says why
# any fold is missing from fold_metrics.
c(attempted = cv$n_folds_attempted, succeeded = cv$n_folds_succeeded)
#> attempted succeeded 
#>         3         3 
cv$fold_status
#>   fold status message
#> 1    1     ok        
#> 2    2     ok        
#> 3    3     ok        

# 3. A metric of your own beside the built-in ones: per fold and pooled.
med_ae <- function(y, yhat) c(MedAE = stats::median(abs(y - yhat)))
cv2 <- cv_spatial(site, "price", "elev", fit_fn = lm_fit, k = 3, seed = 1,
                  metrics = med_ae)
#> cv_spatial(): no folds supplied -- using spatial block k-fold CV (k=3).
cv2$overall$MedAE
#> [1] 2.981173
cv2$fold_metrics[, c("fold", "RMSE", "MedAE")]
#>   fold     RMSE    MedAE
#> 1    1 3.793601 3.536960
#> 2    2 3.198481 2.590833
#> 3    3 3.359840 2.678648
```
