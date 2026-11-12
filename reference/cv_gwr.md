# Cross-validate a GWR with spatial folds

Refits a geographically weighted regression from scratch on each
training fold and scores it on the held-out fold, so the reported error
is the model's at locations it did not see. Use it for an accuracy
figure for a GWR: the in-sample \\R^2\\ that
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
reports on a `gwr_fit` says little, because a local regression with a
small bandwidth can track the training points almost exactly. The
bandwidth is selected again in each fold unless `bandwidth` fixes it, so
the selection itself is inside the cross-validation, with no tuning on
the full data first.

## Usage

``` r
cv_gwr(
  data_sf,
  response_var,
  predictor_vars,
  folds = NULL,
  k = 5,
  seed = 123,
  adaptive = TRUE,
  bandwidth = NULL,
  kernel = c("bisquare", "gaussian", "tricube", "boxcar", "exponential"),
  boundary = NULL,
  pointize = "auto",
  block_size = NULL,
  auto_range = FALSE,
  parallel = FALSE,
  metrics = NULL
)
```

## Arguments

- data_sf:

  An sf layer with the response and the predictors.

- response_var:

  Name of the response column.

- predictor_vars:

  Character vector of predictor column names.

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
  seed drawn per fold before it is fitted. With a seed, the caller's
  random stream is restored afterwards. Default 123.

- adaptive:

  Logical; use an adaptive bandwidth, a number of neighbours, rather
  than a fixed distance. Default `TRUE`.

- bandwidth:

  Optional bandwidth, applied to every fold. For `adaptive = TRUE` a
  whole number of neighbours; for `adaptive = FALSE` a distance in the
  units of the projected CRS the folds are fitted in (lon/lat input is
  projected first, so a value in degrees is read as metres). `NULL`
  (default) selects a bandwidth in each fold from that fold's training
  rows, by the search
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  runs (its section "Bandwidth selection": AICc, or the
  leave-one-location-out sum of squares when observations share a
  location), and reports it in `fold_metrics$bandwidth`.

- kernel:

  The kernel that weights each location's neighbours: `"bisquare"`
  (default), `"gaussian"`, `"tricube"`, `"boxcar"` or `"exponential"`,
  as in
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md).

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

  Accepted so that every `cv_*()` function takes the same arguments, but
  the GWR folds always run one after another in this R process, where
  the warnings [`predict()`](https://rdrr.io/r/stats/predict.html)
  raises for each fold are collected into one. A value asking for more
  than one core raises a warning saying so. Default `FALSE`.

- metrics:

  Optional scoring function of your own, a `function(y, yhat)` returning
  a named numeric vector; its names become columns of `fold_metrics`
  (per fold) and `overall` (pooled) beside the built-in ones. See the
  section "Your own metrics" in
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
  for the contract.

## Value

A list with the elements
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
returns (`overall`, `fold_metrics`, `predictions`, `folds`,
`n_folds_attempted`, `n_folds_succeeded`, `fold_status`, `orphan_rows`,
`n_unknown_ids` and `n_dropped`; see its Value for each, and for the
five fold statuses), plus:

- `formula`:

  The model formula, as text.

- `adaptive`:

  Whether the bandwidth was adaptive.

`fold_metrics` also carries each fold's `bandwidth`. `Adj_R2` is `NA` in
both `overall` and `fold_metrics`: the pooled predictions have no single
parameter count, and a GWR's effective parameter count is not its
predictor count.

A held-out location whose local regression cannot be estimated has
`yhat` `NA`, and one beyond the reach of its fold's fit is predicted and
scored; each case raises one warning for the whole run (classes
`"spatialkit_gwr_not_estimable"` and `"spatialkit_gwr_extrapolation"`),
with the totals over the folds that had any. See
[`predict.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-predict.md).
Where a fold's fit took the response about a centre (a response at a
large offset; see
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)),
`predictions` also carries `y_centre` (that centre, 0 in the rows of a
fold fitted on the raw response), `yhat_centred` and
`y_train_mean_centred` (`yhat` and `y_train_mean` less the centre, which
keep digits that rounding at the response's size would lose), and the
metrics are those of the centred values, as
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
computes them.

## Details

Folds default to spatial blocks (`make_folds(method = "block_kfold")`),
not random ones. A random split leaves a held-out point's neighbours in
the training set, so its score is the error of predicting among the
samples; blocks hold out whole areas, and `folds` says what their score
estimates. Use
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
for the same treatment of a Bayesian GP,
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
for a random forest, and
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
to score several backends on one set of folds.

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

Other cross-validation:
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md),
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md),
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
[`plot_calibration()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_calibration.md),
[`plot_cv_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_cv_metrics.md)

## Examples

``` r
library(sf)
set.seed(1)
n <- 60
dat <- st_as_sf(
  data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
             elev = rnorm(n)),
  coords = c("x", "y"), crs = 32632
)
dat$price <- 10 + 0.01 * (st_coordinates(dat)[, 1] - 5e5) +
  2 * dat$elev + rnorm(n)
cv <- cv_gwr(dat, "price", "elev", k = 3, bandwidth = 30)
#> cv_gwr(): no folds supplied -- using spatial block k-fold CV (k=3).
cv$overall          # pooled over the held-out rows of every fold
#>       RMSE      MAE     MAPE    SMAPE       R2 Adj_R2 n_pred n_MAPE n_SMAPE
#> 1 1.970881 1.640318 11.65992 11.14015 0.620296     NA     60     60      60
cv$fold_metrics     # and fold by fold
#>   fold n_train n_test n_pred    RMSE      MAE      MAPE     SMAPE        R2
#> 1    1      39     21     21 1.34471 1.167010  7.532286  7.606738 0.6918865
#> 2    2      41     19     19 2.31357 1.985458 14.905985 13.676105 0.5117606
#> 3    3      40     20     20 2.16090 1.809408 12.910182 12.441061 0.6691330
#>   Adj_R2 n_MAPE n_SMAPE bandwidth
#> 1     NA     21      21        30
#> 2     NA     19      19        30
#> 3     NA     20      20        30
```
