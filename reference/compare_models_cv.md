# Compare models on the same cross-validation folds

Cross-validates one or more backends on one set of folds and returns a
comparison table. Unlike
[`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md),
this fits the models, once per fold.

## Usage

``` r
compare_models_cv(
  data_sf,
  response_var,
  predictor_vars,
  models = c("GWR", "Bayesian"),
  k = 5,
  seed = 123,
  folds = NULL,
  boundary = NULL,
  pointize = "auto",
  gwr_args = list(),
  bayes_args = list(),
  rf_args = list(),
  summary = c("mean", "median"),
  quiet = FALSE,
  block_size = NULL,
  auto_range = FALSE,
  metrics = NULL
)
```

## Arguments

- data_sf:

  An sf object. When the folds are built here, a feature that is not a
  point is reduced to one, and the layer is read in x and y (see
  [spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)).

- response_var:

  Response column name.

- predictor_vars:

  Predictor column names.

- models:

  Character vector: any subset of `c("GWR", "Bayesian", "RF")`, in any
  order. Each is cross-validated on the same `folds`. Names outside that
  set raise a warning and are dropped; if nothing recognised remains,
  this is an error. A recognised model whose backend package is not
  installed is dropped with a message, so the call still returns the
  models that could run. If none of the requested backends is installed,
  the call stops with `"no viable models."`. Guard with
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) when the
  model set is not known in advance.

- k:

  Number of folds to build when `folds` is `NULL`; not read when `folds`
  is supplied. Default 5.

- seed:

  RNG seed. Default 123.

- folds:

  Optional fold definitions, in any of three shapes: a
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  return value; a bare list of `list(train =, test =)` pairs of
  `..row_id` values; or a vector of fold labels, one per row, which
  becomes leave-that-label-out splits (`blockCV::cv_spatial()` returns
  one as `$folds_ids`). Train and test must be disjoint (a fold that
  trains on its own test rows is not a cross-validation split and is
  refused with an error), and IDs naming no row in the prepared data are
  dropped with a logged count (expected when rows were removed for
  missing values; a sign the folds came from other data when they were
  not). When it is `NULL`, one set of `block_kfold` folds is built and
  every model is scored on it. With `folds = NULL` the score estimates
  the error at the distance the block grid built here puts between the
  held-out points and their training data, which overstates the error of
  a map where the samples cover the area it predicts over. The excess is
  how much the model's error grows from the map's distances to the data
  to the block grid's, and is small for a model whose error changes
  little with the distance from the data.
  [`fold_separation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fold_separation.md)
  with `prediction_points` compares the held-out distances with those of
  the locations the model will predict at, and
  `make_folds(method = "knndm", prediction_points = )` builds folds
  matched to them.

- boundary:

  Optional polygon sf/sfc. When the folds are built here, the block grid
  is clipped to it, as in
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md).
  It is also passed to each backend's `cv_*()` function, which projects
  the data to its CRS before fitting, with supplied `folds` too.

- pointize:

  How a feature that is not a point is reduced to one (see
  [`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md));
  default `"auto"`. It also decides where a polygon or line row falls in
  the shared blocks, so each row is placed by the point every model is
  fitted at.

- gwr_args:

  Extra arguments for
  [`cv_gwr`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md).
  Only names that are formal arguments of
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md)
  are forwarded (it has no `...`), so entries meant for
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  alone (e.g. `longlat`) cannot be passed this way. Anything dropped is
  named in a warning; call
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md)
  directly if you need it.

- bayes_args:

  Extra arguments for
  [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md).
  Forwarded whole as `cv_bayes(fit_args = )`, so an unrecognised name
  raises an error from
  [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
  and is not dropped silently. `compute_loo`, `boundary` and `pointize`
  are overridden by the CV internals.

- rf_args:

  Extra arguments for
  [`cv_rf`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md),
  which passes anything it does not recognise on to
  [`fit_rf_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
  and thence to
  [`ranger::ranger()`](http://imbs-hl.github.io/ranger/reference/ranger.md).

- summary:

  `"mean"` or `"median"` for the Bayesian GP's predictions.

- quiet:

  Logical; suppress this function's progress
  [`message()`](https://rdrr.io/r/base/message.html)s. It does not
  silence R warnings, nor the package's console log echo (see
  [`spatialkit_quiet`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)
  for that). Default `FALSE`.

- block_size:

  Optional minimum block edge length for the shared spatial CV blocks
  (projected CRS units), passed to
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  when `folds` is `NULL`. Default `NULL`. Beside supplied `folds` it has
  nothing to size and raises a warning (class
  `spatialkit_fold_args_ignored`); so does `auto_range = TRUE`.

- auto_range:

  Logical. If `TRUE` and `folds` is `NULL`, the autocorrelation range of
  the response (detrended on `predictor_vars`) is estimated and used as
  the minimum block size of the shared folds, as in
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md).
  Default `FALSE`: geometric blocks. Either way the fold set is built
  once and every backend is scored on it. When it cannot be built (a
  `block_size` or estimated range that leaves a single block, say) the
  call is an error, as it is for each backend on its own; no model is
  scored on a design other than the one asked for.

- metrics:

  Optional scoring function of your own, handed to every backend's
  `cv_*()`: a `function(y, yhat)` returning a named numeric vector,
  applied per fold and to each backend's pooled predictions, whose names
  become columns of `by_fold` and `overall` beside the built-in ones.
  See "Your own metrics" on
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
  for the contract. Because the backends are scored on the same folds
  and, in `overall`, on the same rows (see Value), the columns are
  comparable across rows of `overall`.

## Value

A list with these elements:

- `overall`:

  One row per model with the pooled metrics, `n_pred` (the rows they are
  computed on, see "Rows the overall table is computed on"), the
  coverage and CRPS columns described in "Coverage and CRPS in the
  overall table" when a Bayesian GP ran, and `model` as its last column.
  `MAPE` and `SMAPE` skip rows whose denominator is zero; `n_MAPE` and
  `n_SMAPE` count the rows they used, and for a count or other
  non-Gaussian response RMSE, MAE and the Bayesian GP's CRPS are the
  numbers to compare (see "Percentage errors on responses with zeros"
  and "Which metrics survive a non-Gaussian response" in
  [`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)).

- `by_fold`:

  The `fold_metrics` of every model that ran, stacked, with a `model`
  column.

- `gwr_cv`, `bayes_cv`, `rf_cv`:

  The full
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
  and
  [`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
  result of each model that ran.

## Coverage and CRPS in the overall table

A model can predict well on average and still be badly calibrated, so a
comparison read from RMSE alone can prefer the model whose uncertainty
is wrong. Heaton et al. (2019) found good point prediction routinely
alongside poor interval coverage. When `"Bayesian"` is among the models
that ran, `overall` therefore also carries the columns of
`cv_bayes()$predictive_coverage`: `coverage_50`, `coverage_80`,
`coverage_95` (the share of held-out rows inside the posterior
predictive interval at each level, averaged across folds weighted by
each fold's `n_pred`) and `mean_CRPS` (the continuous ranked probability
score, lower is better). They are `NA` on the GWR and RF rows, which
produce point predictions and no draws, and absent when no Bayesian GP
ran. Read coverage against its nominal level: 0.95 at `coverage_95` is
calibrated, well below it is overconfident, well above it is wider than
it needs to be.

## Rows the overall table is computed on

Shared folds do not guarantee shared rows: a model that fails on a fold
(GWR with a fixed bandwidth across a gap in the data, say) or predicts
`NA` for some rows pools fewer rows, usually without the hardest ones.
When the models predicted different rows, the function warns and
recomputes every model's pooled metrics, your own `metrics` included, on
the rows all of them predicted, so `n_pred` is the same on every row
that has predictions; a model that predicted nothing stays an `NA` row.
Each model's metrics over all the rows it predicted stay in its `*_cv`
element and in `attr(overall, "all_rows")`, a table of the same shape.
`by_fold` and the Bayesian coverage and CRPS columns are per fold and
are not recomputed.

Every model is cross-validated on the layer as passed, so the warnings
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
raises about that layer (class `spatialkit_rows_dropped`, counting the
incomplete rows dropped, and `spatialkit_constant_predictor`) are raised
once for the comparison, not once per model.

Only the models that ran appear, so check which names are present: a
backend whose package is missing is dropped with a message. When no
requested backend is available the function stops with
`"no viable models."` rather than returning an empty comparison.

## References

Heaton, M. J., Datta, A., Finley, A. O., Furrer, R., Guinness, J.,
Guhaniyogi, R., Gerber, F., Gramacy, R. B., Hammerling, D., Katzfuss,
M., Lindgren, F., Nychka, D. W., Sun, F. and Zammit-Mangion, A. (2019).
A case study competition among methods for analyzing large spatial data.
*Journal of Agricultural, Biological and Environmental Statistics*,
24(3), 398–425.
[doi:10.1007/s13253-018-00348-w](https://doi.org/10.1007/s13253-018-00348-w)

## See also

[`plot_cv_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_cv_metrics.md)
to draw the per-fold scores of every model, and
[`plot_calibration()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_calibration.md)
for the Bayesian GP's intervals.

Other cross-validation:
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md),
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
[`plot_calibration()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_calibration.md),
[`plot_cv_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_cv_metrics.md)

## Examples

``` r
library(sf)
set.seed(1)
n <- 120
dat <- st_as_sf(
  data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
             elev = rnorm(n)),
  coords = c("x", "y"), crs = 32632
)
dat$price <- 10 + 0.01 * (st_coordinates(dat)[, 1] - 5e5) +
  2 * dat$elev + rnorm(n)
# The random forest, and GWR beside it when sp is installed, on the same
# folds.  "Bayesian" is left out here because a Stan fit per fold takes
# minutes; add it to `models` for a run to report.
models <- c("RF", if (requireNamespace("sp", quietly = TRUE)) "GWR")
cmp <- compare_models_cv(dat, "price", "elev", models = models, k = 3,
                         rf_args = list(num_trees = 100),
                         gwr_args = list(bandwidth = 30))
#> compare_models_cv(): running CV for GWR ...
#> compare_models_cv(): running CV for RF ...
cmp$overall[, c("model", "RMSE", "R2", "n_pred")]
#>   model     RMSE        R2 n_pred
#> 1   GWR 2.281133 0.6738657    120
#> 2    RF 3.711068 0.1368361    120
```
