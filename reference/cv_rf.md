# Cross-validate a random forest with spatial folds

A wrapper over
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
that refits a ranger forest on each training fold. By default it holds
out whole spatial blocks, so most neighbours of a held-out point are not
in the training set; the forest's out-of-bag error holds out single
rows. The folds you pass decide what the score measures (see "Fold
methods" in
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)).

## Usage

``` r
cv_rf(
  data_sf,
  response_var,
  predictor_vars,
  folds = NULL,
  k = 5,
  seed = 123,
  parallel = FALSE,
  block_size = NULL,
  auto_range = FALSE,
  boundary = NULL,
  pointize = "auto",
  metrics = NULL,
  ...
)
```

## Arguments

- data_sf:

  An sf object.

- response_var:

  Response column name.

- predictor_vars:

  Predictor column names.

- folds:

  Optional fold definitions, in any of three shapes: a
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  return value; a bare list of `list(train =, test =)` pairs of
  `..row_id` values; or a vector of fold labels, one per row, which
  becomes leave-that-label-out splits. The label vector is how folds
  built by another package are used here, since
  `blockCV::cv_spatial()`'s `$folds_list` holds two *unnamed* vectors
  per fold and is refused by name. That function returns a label vector
  as `$folds_ids`. Train and test must be disjoint: a fold that trains
  on its own test rows is not a cross-validation split and is refused
  with an error. IDs naming no row in the prepared data are dropped with
  a logged count. Built via `block_kfold` when `NULL`, as in
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md).
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

  Number of folds when `folds` is `NULL`. Default 5.

- seed:

  RNG seed. Default 123. It seeds fold construction and, through a
  per-fold draw, each fold's forest, so two different seeds give
  different results even on identical `folds`. To grow every fold's
  forest from one fixed ranger seed instead, call
  [`cv_spatial`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
  with your own `fit_fn` wrapping `fit_rf_model(seed = )`.

- parallel:

  Passed to
  [`cv_spatial`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md).
  Default `FALSE`. Under forked workers each fold's forest runs
  single-threaded unless `num_threads` is passed explicitly, so
  `parallel = 4` means four threads in total, where it would otherwise
  mean four times the session's `mc.cores`.

- block_size, auto_range, boundary:

  Passed to
  [`cv_spatial`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md).
  `block_size` and `auto_range = TRUE` size the blocks built when
  `folds` is `NULL`; beside supplied `folds` they have nothing to size
  and raise a warning of class `spatialkit_fold_args_ignored`.

- pointize:

  How non-POINT geometry is reduced to a point before fitting; passed to
  [`cv_spatial`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md).
  Default `"auto"`.

- metrics:

  Optional scoring function of your own, passed to
  [`cv_spatial`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md);
  see **Your own metrics** there.

- ...:

  Passed to
  [`fit_rf_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
  on every fold: `num_trees`, `mtry`, `importance`, `include_coords` and
  so on. `data_sf`, `response_var`, `predictor_vars` and
  `.already_prepped` are set by this function and must not be passed
  here (every fold would fail with "matched by multiple actual
  arguments"). `seed` is this function's own argument and never reaches
  [`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
  through here; see `seed` above for growing every fold's forest from
  one fixed seed.

## Value

The
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
result. `MAPE` and `SMAPE` skip rows whose denominator is zero; `n_MAPE`
and `n_SMAPE` count the rows they used (see "Percentage errors on
responses with zeros" in
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)).

## Comparing the forest with and without coordinates

Whether `include_coords = TRUE` lowers the error depends on how far from
the training data the predictions are (see "Coordinates as predictors"
in
[`fit_rf_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)),
so it is a comparison to make on your own data, at the block size the
map will be used at. Score the two forests on the same folds: build the
folds once with
[`make_folds`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
and pass the result as `folds` to two calls, one with
`include_coords = TRUE` and one without. Each fold then holds out the
same rows in both, so the two `overall` rows compare like with like and
the two `fold_metrics` tables pair up fold by fold. A difference in
`overall` that is small beside the spread of the per-fold differences is
not evidence for either forest. The example below runs the comparison.
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
takes one set of forest arguments (`rf_args`) per call, so it scores
only one of the two forests.

## See also

[`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
for the forest's arguments, and
[`plot_cv_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_cv_metrics.md)
to draw the per-fold scores.

Other cross-validation:
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md),
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
[`plot_calibration()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_calibration.md),
[`plot_cv_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_cv_metrics.md)

## Examples

``` r
library(sf)
set.seed(1)
n <- 200
dat <- st_as_sf(
  data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
             a = rnorm(n)),
  coords = c("x", "y"), crs = 32632
)
# An offset keeps the response away from zero, where MAPE is undefined.
dat$z <- 10 + 2 * dat$a + rnorm(n, 0, 0.3)
cv_rf(dat, "z", "a", k = 4)$overall
#> cv_rf(): no folds supplied -- using spatial block k-fold CV (k=4).
#>       RMSE       MAE     MAPE    SMAPE        R2 Adj_R2 n_pred n_MAPE n_SMAPE
#> 1 0.417048 0.3264843 3.424117 3.405887 0.9587525     NA    200    200     200

# The forest with and without the coordinates, on the same blocked
# folds.  The second call logs the include_coords caution once per
# session.  This response has no spatial structure, so here the
# coordinates can only cost accuracy.
folds   <- make_folds(dat, k = 4, method = "block_kfold", seed = 1)
plain   <- cv_rf(dat, "z", "a", folds = folds, num_trees = 200)
with_xy <- cv_rf(dat, "z", "a", folds = folds, num_trees = 200,
                 include_coords = TRUE)
c(plain = plain$overall$RMSE, with_xy = with_xy$overall$RMSE)
#>     plain   with_xy 
#> 0.4208380 0.7656067 
# Fold by fold: positive where the coordinates lowered the error.
plain$fold_metrics$RMSE - with_xy$fold_metrics$RMSE
#> [1] -0.5017832 -0.2412238 -0.3906285 -0.1638759
```
