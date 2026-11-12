# Fit a random forest with ranger

Fits a regression random forest on an `sf` layer and returns it as a
`spatial_fit`, so it works with
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
[`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md),
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
and the [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
like any other backend. A forest suits a response that is non-linear in
its predictors when predictive accuracy matters more than an
interpretable model.

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
  importance is biased towards continuous and high-cardinality
  predictors (Strobl et al. 2007), so permutation importance is the
  default although it costs more.

- include_coords:

  Add the coordinates as predictors. Default `FALSE`, which leaves the
  forest with no spatial information; see *Coordinates as predictors*.

- replace:

  Logical; grow each tree on a bootstrap sample drawn *with* replacement
  (`TRUE`, ranger's default and this one) or on a subsample drawn
  without replacement (`FALSE`). Strobl et al. (2007) show that
  bootstrap sampling with replacement is itself a source of bias in
  variable importance toward predictors with many distinct values or
  categories, and recommend subsampling without replacement. The
  forest's predictions change slightly with the choice. The setting is
  recorded in `$info` and printed with the fit.

- sample_fraction:

  Fraction of rows drawn for each tree. `NULL` (default) uses ranger's
  rule: all rows when `replace = TRUE`, 0.632 (the expected share of
  distinct rows in a bootstrap sample) when `replace = FALSE`. A single
  number in (0, 1\] overrides it. `replace = FALSE` with
  `sample_fraction = 1` grows every tree on every row, so nothing is out
  of bag: the fit warns (see "Fitted values are out-of-bag").

- seed:

  Seed passed to ranger. Default 123.

- num_threads:

  Threads for ranger. Default `NULL` means `getOption("mc.cores", 1L)`:
  one thread unless the session has set more, the convention brms and
  parallel use. Predictions do not depend on the thread count; pass
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)
  to use every core.

- .already_prepped:

  Internal; leave at the default.

- ...:

  Passed to
  [`ranger::ranger()`](http://imbs-hl.github.io/ranger/reference/ranger.md).
  ranger's own spellings of the arguments this function already sets
  (`num.trees`, `min.node.size`, `num.threads`, `mtry`, `importance`,
  `seed`, `replace`, `sample.fraction`, `x`, `y`) are rejected with a
  message naming the wrapper argument to use instead.
  `respect.unordered.factors` passed here overrides the `"order"` this
  function sets when a predictor is categorical; see *Categorical
  predictors*.

## Value

An `rf_fit` object (inherits from `spatial_fit`). `$info` carries
`num_trees`, `mtry`, `min_node_size`, `importance_type`, `importance` (a
named numeric, or `NULL` when `importance = "none"`), `include_coords`,
`categorical_predictors` (the unordered factor and character predictors,
possibly none) and `respect_unordered_factors` (how ranger was told to
split them: `"order"` unless overridden through `...`, and `NA` when
there is none and nothing was passed), `replace` and `sample_fraction`
(the sampling each tree was grown on, with `sample_fraction` resolved to
the number ranger used), `oob_rmse` and `oob_r_squared` (each `NA_real_`
when ranger did not compute it; forwarding `oob.error = FALSE` through
`...` is one way to get there), `fitted_are_oob` (always `TRUE`;
[`summary()`](https://rdrr.io/r/base/summary.html) reads it to label its
metrics), `seed` and `n_dropped` (the rows
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
removed for missing or non-finite values or a bad geometry, so `$n` can
be read against `nrow(data_sf)`), and `response_centre` and
`response_scale` (below). The raw forest is in `$engine`.

The forest is grown on the response taken about its centre where that
centre is more than 1000 times the response's spread (as at 1e8 plus a
spread of 1), and divided by a power of two where the response is so
large or its spread so small that its squares overflow or lose their
bits (beyond about 1e100 or 1e-100). The splits of a forest do not
change under a shift or a scale of the response, and ranger's split
scores, formed from squares of the response, read the rounding of the
offset otherwise. `$info$response_centre` (0 where none is taken) and
`$info$response_scale` (1 where none is applied) record them, and they
are put back on
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html), `oob_rmse` and
`importance` (squared units, so multiplied by the scale twice; beyond
the range of doubles it is `Inf` or loses its bits, and a WARN line says
so). The permutation importance under
`scale.permutation.importance = TRUE` and the impurity importance under
`splitrule = "maxstat"` have no units and are left as ranger reports
them. The forest in `$engine` predicts
`(y - response_centre) / response_scale`, and its own `predictions` and
`prediction.error` are in those units. Every other response is passed to
ranger as it is.

## Coordinates as predictors

`include_coords = TRUE` adds the x and y coordinates to the predictors.
They are the only spatial information a forest can be given: without
them it sees the covariates alone, and
[`print()`](https://rdrr.io/r/base/print.html) says "no spatial
information used". The default is `FALSE`.

The gain from the coordinates depends on where the predictions are
wanted. A forest splits on location as on any column, so it follows the
training surface closely between sampled points, and has only the
nearest edge of the data to go on beyond them. Between sampled points
the coordinates can lower the error, most where the covariates carry
little of the spatial structure themselves, and the gain fades with
distance from the training data. Beyond the autocorrelation range, or in
an area the data never visited, the coordinates can make predictions
worse (Meyer et al. 2019). With covariates that carry the spatial
structure themselves, expect less from the coordinates.

Score the forests with and without coordinates on the same folds, at the
distances the map will be used at (see "Fold methods" in
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md);
"Out-of-bag error" below says why the out-of-bag error does not serve).
Build the folds once and pass them to
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
with and without the coordinates:


    folds   <- make_folds(dat, k = 5, method = "block_kfold")
    plain   <- cv_rf(dat, "z", c("a", "b"), folds = folds)
    with_xy <- cv_rf(dat, "z", c("a", "b"), folds = folds,
                     include_coords = TRUE)
    c(plain = plain$overall$RMSE, with_xy = with_xy$overall$RMSE)

The example on
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
runs this comparison and reads it fold by fold. Setting
`include_coords = TRUE` logs this advice once per session, as a log line
rather than an R warning, and not once for every fold of a
cross-validation.

## Categorical predictors

A predictor that is an unordered factor or a character column is split
with `respect.unordered.factors = "order"`: ranger orders its levels by
mean response once, before growing, and splits on that order, the
handling ranger's own help recommends. ranger's default for a regression
forest, `"ignore"`, splits on the order the levels are stored in
(alphabetical, for most data), which for a nominal variable is
arbitrary. Pass `respect.unordered.factors` through `...` to choose
otherwise (`"ignore"`, or `"partition"` for every two-way split of the
levels, which is slow beyond a few levels). Ordered factors are left
alone and split on their own order, and with `splitrule = "extratrees"`
ranger's own default for that rule, `"partition"`, is kept. The setting
is recorded in `$info$respect_unordered_factors` and printed with the
fit.

## Out-of-bag error

ranger's out-of-bag (OOB) error holds each observation out of the trees
that did not sample it. That is a random hold-out: the trees that did
not see a point almost always saw its neighbours, so under spatial
autocorrelation it is the error of predicting among the samples, as
random k-fold is, and reads low for predictions made farther from the
data than the samples sit from each other (see
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)).
It is reported as `$info$oob_rmse` and `$info$oob_r_squared` and
labelled as OOB wherever it appears. Use
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
for a spatial estimate.

## Fitted values are out-of-bag

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) on an `rf_fit`
returns out-of-bag predictions, as the random forest packages themselves
do. In-sample predictions from a forest come close to reproducing the
training data, and their \\R^2\\ says little about new data.
[`summary()`](https://rdrr.io/r/base/summary.html) therefore reports
out-of-bag metrics for a forest and in-sample metrics for a GWR or a
Bayesian GP fit, and the two do not compare;
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
scores every backend on the same folds.

A row that every tree sampled has no out-of-bag prediction, and ranger
reports `NaN` for it. That is every row under `replace = FALSE` with
`sample_fraction = 1`, and a few under a small `num_trees`. The fit
warns with the count;
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) are `NaN` on
those rows, and [`summary()`](https://rdrr.io/r/base/summary.html) says
how many rows its metrics were computed on. With no row out of bag at
all, the OOB error (`NA` in `$info$oob_rmse` and `$info$oob_r_squared`)
and the permutation importance (`NaN`) are undefined too, and
[`print()`](https://rdrr.io/r/base/print.html) says so.
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
scores its fold forests on the held-out rows, never out of bag, so it
warns once with the number of folds affected rather than once per fold.

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

[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
to score the forest on spatial folds;
[`predict()`](https://rdrr.io/r/stats/predict.html)
([spatial_fit-predict](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-predict.md))
and
[`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md)
to predict;
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
to mark where the predictions apply (it weights the predictors by this
forest's importance).

Other model fitting:
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md),
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md),
[`gp_lengthscale_bounds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gp_lengthscale_bounds.md),
[`new_spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md)

## Examples

``` r
library(sf)
set.seed(1)
n <- 150
dat <- st_as_sf(
  data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
             a = rnorm(n), b = rnorm(n)),
  coords = c("x", "y"), crs = 32632
)
dat$z <- 2 * dat$a - dat$b + rnorm(n, 0, 0.3)
fit <- fit_rf_model(dat, "z", c("a", "b"))
fit
#> <Random Forest (ranger)> spatial model fit
#>   Formula : z ~ a + b
#>   n       : 150
#>   CRS     : EPSG:32632
#>   Trees   : 500 (mtry = 1, min node = 5)
#>   Coords as predictors: no - no spatial information used; see ?fit_rf_model
#>   Sampling: bootstrap, with replacement (100.0% of rows per tree)
#>   OOB RMSE: 0.6995   OOB R^2: 0.9093
#>   Importance (permutation): a=6.886, b=1.583
#> 
#>   OOB is a random hold-out: under spatial autocorrelation it is the error
#>   of predicting among the samples, and reads low for predictions made
#>   farther from the data than the samples sit from each other. Score a map
#>   that predicts farther with cv_rf() on folds at its distances
#>   (?fold_separation).
fit$info$importance     # permutation importance of a and b
#>        a        b 
#> 6.886377 1.583079 
```
