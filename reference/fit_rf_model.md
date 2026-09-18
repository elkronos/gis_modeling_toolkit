# Fit a random forest via ranger

Fits a regression random forest on an `sf` dataset and returns it as a
`spatial_fit`, so it works with
[`cv_spatial`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
[`predict_surface`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md),
[`area_of_applicability`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
and the [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
like any other backend.

## Usage

``` r
fit_rf_model(
  data_sf,
  response_var,
  predictor_vars,
  num_trees = 500L,
  mtry = NULL,
  min_node_size = NULL,
  importance = c("permutation", "impurity", "none"),
  include_coords = FALSE,
  replace = TRUE,
  sample_fraction = NULL,
  seed = 123L,
  num_threads = NULL,
  .already_prepped = FALSE,
  ...
)
```

## Arguments

- data_sf:

  An sf object with response, predictors and geometry.

- response_var:

  Response column name.

- predictor_vars:

  Predictor column names.

- num_trees:

  Number of trees. Default 500.

- mtry:

  Predictors sampled per split. `NULL` uses ranger's default.

- min_node_size:

  Minimum node size. `NULL` uses ranger's default.

- importance:

  `"permutation"` (default), `"impurity"` or `"none"`. Impurity
  importance is biased toward continuous and high-cardinality predictors
  (Strobl et al. 2007), which is why permutation is the default despite
  costing more.

- include_coords:

  Add the coordinates as predictors. Default `FALSE`; see above.

- replace:

  Logical; grow each tree on a bootstrap sample drawn *with* replacement
  (`TRUE`, ranger's default and this one) or on a subsample drawn
  without replacement (`FALSE`). Strobl et al. (2007) show that
  bootstrap sampling with replacement is itself a source of bias in
  variable importance toward predictors with many distinct values or
  categories, and recommend subsampling without replacement. The
  forest's predictions change slightly with the choice; the default is
  kept at ranger's so that an existing script fits the same forest, and
  the setting is recorded in `$info` and printed with the fit.

- sample_fraction:

  Fraction of rows drawn for each tree. `NULL` (default) uses ranger's
  rule: all rows when `replace = TRUE`, 0.632 — the expected share of
  distinct rows in a bootstrap sample — when `replace = FALSE`. A single
  number in (0, 1\] overrides it.

- seed:

  Seed passed to ranger. Default 123.

- num_threads:

  Threads for ranger. Default `NULL` means `getOption("mc.cores", 1L)`:
  one thread unless the session has opted in to more, the convention
  brms and parallel use. ranger's own default is every core on the
  machine, which is the wrong default for a package function (a check
  farm limits jobs to two cores, and `cv_rf(parallel = )` would multiply
  it by the worker count). Predictions do not depend on the thread
  count, only speed does; pass
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)
  to use them all.

- .already_prepped:

  Internal; skip
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
  because the caller has already projected and cleaned the data.

- ...:

  Passed to
  [`ranger::ranger()`](http://imbs-hl.github.io/ranger/reference/ranger.md).
  ranger's own spellings of the arguments this function already sets
  (`num.trees`, `min.node.size`, `num.threads`, `mtry`, `importance`,
  `seed`, `replace`, `sample.fraction`, `x`, `y`) are rejected with a
  message naming the wrapper argument to use instead – passing them here
  would reach `ranger()` twice and fail the call.

## Value

An `rf_fit` object (inherits from `spatial_fit`). `$info` carries
`num_trees`, `mtry`, `min_node_size`, `importance_type`, `importance` (a
named numeric, or `NULL` when `importance = "none"`), `include_coords`,
`replace` and `sample_fraction` (the sampling each tree was grown on,
with `sample_fraction` resolved to the number ranger used), `oob_rmse`
and `oob_r_squared` (each `NA_real_` when ranger did not compute it –
forwarding `oob.error = FALSE` through `...` is one way to get there),
`fitted_are_oob` (always `TRUE`;
[`summary()`](https://rdrr.io/r/base/summary.html) reads it to label its
metrics), `seed` and `n_dropped` (the rows
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
removed for missing or non-finite values or a bad geometry, so `$n` can
be read against `nrow(data_sf)`). The raw forest is in `$engine`.

## Coordinates are not predictors by default

Handing a random forest the x and y coordinates lets it reproduce the
training surface almost exactly by memorising location, and then fail
badly anywhere it has not seen. Random cross-validation will not catch
this – nearby points leak between folds, so the memorised surface scores
well – which is how the practice became common. Meyer et al. (2019) show
the collapse directly. `include_coords` therefore defaults to `FALSE`,
and setting it to `TRUE` logs a caution (it is a deliberate choice, so
it is not raised as an R warning). If you do use it, score the model
with
[`cv_spatial`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
and blocked folds, never with the out-of-bag error.

## The out-of-bag error is a random hold-out

`ranger`'s OOB error holds each observation out of the trees that did
not sample it. That is a random hold-out, so under spatial
autocorrelation it is optimistic for exactly the reason random k-fold
is: the trees that "did not see" a point almost certainly saw its
neighbours. It is reported as `$info$oob_rmse` and `$info$oob_r_squared`
and labelled as OOB everywhere it appears. Use
[`cv_rf`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
for a spatial estimate.

## What fitted() returns

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) on an `rf_fit`
returns **out-of-bag** predictions, not in-sample ones, following the
convention of the random forest packages themselves. In-sample
predictions from a forest are close to memorisation and would make
[`summary()`](https://rdrr.io/r/base/summary.html) report a fictitious
R-squared. The consequence is that
[`summary()`](https://rdrr.io/r/base/summary.html) means something
different here than for a `gwr_fit` or `bayesian_fit`, whose fitted
values are in-sample: do not compare the two directly.
[`compare_models_cv`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
exists for that.

## References

Meyer, H., Reudenbach, C., Wöllauer, S. and Nauss, T. (2019). Importance
of spatial predictor variable selection in machine learning applications
– moving from data reproduction to spatial prediction. *Ecological
Modelling* 411, 108815.
[doi:10.1016/j.ecolmodel.2019.108815](https://doi.org/10.1016/j.ecolmodel.2019.108815)

Strobl, C., Boulesteix, A.-L., Zeileis, A. and Hothorn, T. (2007). Bias
in random forest variable importance measures: illustrations, sources
and a solution. *BMC Bioinformatics* 8, 25.
[doi:10.1186/1471-2105-8-25](https://doi.org/10.1186/1471-2105-8-25)

## See also

[`cv_rf`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
for a spatially blocked performance estimate,
[`area_of_applicability`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md),
which can take `weights = pmax(fit$info$importance, 0)`.

Other model fitting:
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md),
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md),
[`new_spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md),
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)

## Examples

``` r
if (requireNamespace("ranger", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  n <- 150
  dat <- st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
               a = rnorm(n), b = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
  dat$z <- 2 * dat$a - dat$b + rnorm(n, 0, 0.3)
  fit <- fit_rf_model(dat, "z", c("a", "b"))
  fit
  fit$info$importance
}
#>        a        b 
#> 6.886377 1.583079 
```
