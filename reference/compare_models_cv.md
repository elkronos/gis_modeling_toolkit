# Cross-validated comparison of spatial models

Fits and cross-validates one or more model types, returning a unified
comparison table. Unlike
[`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md),
this function does perform fitting (inside CV folds), because CV
inherently requires repeated fitting.

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

  An sf object.

- response_var:

  Response column name.

- predictor_vars:

  Predictor column names.

- models:

  Character vector: any subset of `c("GWR", "Bayesian", "RF")`, in any
  order. Each is cross-validated on the same `folds`. Names outside that
  set raise a warning and are dropped; if nothing recognised remains,
  this is an error. There is no silent fallback. A recognised model
  whose backend package is not installed is dropped with a message so
  the call still returns the models that could run. But if *none* of the
  requested backends is installed, nothing is left to compare and the
  call errors with `"no viable models."`. Guard with
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) when the
  model set is not known in advance.

- k:

  Number of folds. Default 5.

- seed:

  RNG seed. Default 123.

- folds:

  Optional fold definitions: a
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  return value, or a bare list of `list(train =, test =)` pairs of
  `..row_id` values. Train and test must be disjoint (a fold that trains
  on its own test rows is not a cross-validation split and is refused
  with an error), and IDs naming no row in the prepared data are dropped
  with a logged count (expected when rows were removed for missing
  values; a sign the folds came from other data when they were not).

- boundary:

  Optional polygon sf/sfc.

- pointize:

  Geometry coercion strategy.

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

  "mean" or "median" for Bayesian predictions.

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
  when `folds` is `NULL`. Default `NULL`.

- auto_range:

  Logical. If `TRUE` and `folds` is `NULL`, the autocorrelation range of
  the response (detrended on `predictor_vars`) is estimated and used as
  the minimum block size of the shared folds, as in
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md).
  Default `FALSE`: geometric blocks, as before this argument existed.
  Either way the fold set is built once and every backend is scored on
  it.

- metrics:

  Optional scoring function of your own, handed to every backend's
  `cv_*()`: a `function(y, yhat)` returning a named numeric vector,
  applied per fold and to each backend's pooled predictions, whose names
  become columns of `by_fold` and `overall` beside the built-in ones.
  See **Your own metrics** on
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
  for the contract. Because the three backends are scored on the same
  folds, the columns are comparable across rows of `overall`.

## Value

A list with overall, by_fold, and per-model cv_results (`gwr_cv`,
`bayes_cv`, `rf_cv` for the models that ran). `overall` has one row per
model with the pooled metrics, the coverage and CRPS columns described
above when a Bayesian model ran, and `model` as its last column. Only
the models that actually ran appear, so check which names are present;
there is not always one entry per requested model, because a backend
whose package is missing is dropped with a message. When **no**
requested backend is available there is nothing to return and the
function errors with `"no viable models."` instead of returning an empty
comparison.

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
produce point predictions and no draws, and absent when no Bayesian
model ran. Read coverage against its nominal level: 0.95 at
`coverage_95` is calibrated, well below it is overconfident, well above
it is wider than it needs to be.

## Percentage errors on responses with zeros

`MAPE` divides by the observed value and `SMAPE` by \\\|y\| +
\|\hat{y}\|\\, so neither is defined where its denominator is zero.
Neither returns `Inf` or `NaN`. Both are averaged over the rows whose
denominator is non-zero, and are `NA` when no row qualifies. The
`n_MAPE` and `n_SMAPE` columns record how many rows that was; the `n`
column counts finite observation/prediction pairs. Read a percentage
error next to its count: when `n_MAPE < n`, `MAPE` is an average over a
subset of the data, whatever its value.

This bites on any response taking exact zeros: counts, rainfall,
abundance, claim amounts. On a zero-inflated response with 62 zeros out
of 120, `MAPE` is an average over the 58 non-zero rows, which
`n_MAPE = 58` now says. `SMAPE` fails differently and more subtly: it
drops the rows where observation and prediction are both near zero
(which on a well-fitted zero-inflated model are the rows it got
*right*), so it averages the harder rows only and reads worse than the
fit deserves; `n_SMAPE` shows how many rows it kept, and the count is
only a label, not a repair.

`RMSE`, `MAE` and \\R^2\\ use every finite row and are unaffected;
prefer them whenever the response can be zero. For a Bayesian fit,
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
additionally reports CRPS and interval coverage, which are proper
scoring rules and have no such failure mode.

## Which metrics survive a non-Gaussian response

RMSE and MAE are defined for any numeric response and are what to read
for a count, a rate or a bounded outcome. MAPE and SMAPE assume a
response that is rarely zero (see the previous section), and R-squared
and adjusted R-squared compare residual variance to total variance,
which is the right comparison for a Gaussian response and a loose one
for anything whose variance tracks its mean. None of the four is wrong
to compute; each is Gaussian-shaped thinking, and on a Poisson or
zero-inflated response should be read as a rough summary rather than a
score.

For the Bayesian backend,
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
additionally reports CRPS and interval coverage at 50, 80 and 95
percent. Both are proper scoring rules computed from posterior draws, so
they are meaningful for any `family` the backend accepts, and they are
the numbers to compare when the response is not Gaussian. When every
fold fails, the `fold_metrics` frame
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
returns carries the CRPS column but not the `coverage_*` columns, so
code that reads those columns must tolerate their absence.

## References

Heaton, M. J., Datta, A., Finley, A. O., Furrer, R., Guinness, J.,
Guhaniyogi, R., Gerber, F., Gramacy, R. B., Hammerling, D., Katzfuss,
M., Lindgren, F., Nychka, D. W., Sun, F. and Zammit-Mangion, A. (2019).
A case study competition among methods for analyzing large spatial data.
*Journal of Agricultural, Biological and Environmental Statistics*,
24(3), 398–425.
[doi:10.1007/s13253-018-00348-w](https://doi.org/10.1007/s13253-018-00348-w)

## See also

Other model evaluation:
[`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md),
[`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md),
[`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)

## Examples

``` r
if (requireNamespace("ranger", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  n <- 120
  dat <- st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), elev = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
  dat$price <- 10 + 0.01 * st_coordinates(dat)[, 1] + 2 * dat$elev + rnorm(n)
  cmp <- compare_models_cv(dat, "price", "elev", models = "RF", k = 3,
                           rf_args = list(num_trees = 100))
  cmp$overall
}
#> compare_models_cv(): running CV for RF ...
#>       RMSE      MAE     MAPE    SMAPE        R2 Adj_R2 n_pred n_MAPE n_SMAPE
#> 1 3.711068 3.128498 21.71902 20.99972 0.1368361     NA    120    120     120
#>   model
#> 1    RF
```
