# Print a fitted model

Prints the backend, the formula, the number of observations and the CRS,
then lines for the backend:

- a random forest (`rf_fit`): the number of trees, `mtry`, the node
  size, whether the coordinates were used as predictors, how each tree's
  rows were sampled, the out-of-bag error and the largest importances. A
  forest fitted without them prints "no spatial information used": it
  sees the covariates alone. One fitted with `include_coords = TRUE`
  prints a note on how to score it (see "Coordinates as predictors" and
  "Out-of-bag error" in
  [`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)).
  When a predictor is an unordered factor or a character column, a line
  says how ranger was told to split it;

- a GWR (`gwr_fit`): the bandwidth (a count of neighbours, or a distance
  with the CRS's unit), the kernel, the AICc, the effective number of
  parameters and, when any local window is of deficient rank, the number
  of such windows;

- a Bayesian GP (`bayesian_fit`): the GP basis size, the smallest
  effective range the basis resolves in CRS units, the LOOIC when the
  fit computed it (`compute_loo = TRUE`, the default of
  [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)),
  the slope prior, the standardised predictors, and a flagged line when
  the sampler's diagnostics or the basis check failed.

A custom subclass prints the common lines only. Typing the object's name
calls this method; use it to confirm that a fit used the data,
predictors and CRS you meant. For fit quality use
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
or
[`summary.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.spatial_fit.md).
Apart from a forest's out-of-bag error and a Bayesian GP's LOOIC, an
approximate leave-one-out score of single rows, nothing printed here is
an out-of-sample score, and neither of those is a spatial one.

## Usage

``` r
# S3 method for class 'spatial_fit'
print(x, ...)

# S3 method for class 'rf_fit'
print(x, ...)
```

## Arguments

- x:

  A `spatial_fit` object.

- ...:

  Not used.

## Value

`x`, invisibly.

## See also

[`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md),
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
and
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md),
which make the objects.

Other fit methods:
[`plot.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.spatial_fit.md),
[`spatial_fit-coef`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-coef.md),
[`spatial_fit-fitted`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-fitted.md),
[`spatial_fit-predict`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-predict.md),
[`spatial_fit-residuals`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-residuals.md),
[`summary.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.spatial_fit.md)

## Examples

``` r
library(sf)
set.seed(1)
pts <- st_as_sf(
  data.frame(x = 5e5 + runif(60, 0, 1000), y = 5e6 + runif(60, 0, 1000),
             a = rnorm(60)),
  coords = c("x", "y"), crs = 32632
)
pts$z <- 10 + 2 * pts$a + rnorm(60, 0, 0.3)
fit <- fit_rf_model(pts, "z", "a", num_trees = 50, seed = 1)
fit                     # print.rf_fit()
#> <Random Forest (ranger)> spatial model fit
#>   Formula : z ~ a
#>   n       : 60
#>   CRS     : EPSG:32632
#>   Trees   : 50 (mtry = 1, min node = 5)
#>   Coords as predictors: no - no spatial information used; see ?fit_rf_model
#>   Sampling: bootstrap, with replacement (100.0% of rows per tree)
#>   OOB RMSE: 0.4094   OOB R^2: 0.9502
#>   Importance (permutation): a=6.034
#> 
#>   OOB is a random hold-out: under spatial autocorrelation it is the error
#>   of predicting among the samples, and reads low for predictions made
#>   farther from the data than the samples sit from each other. Score a map
#>   that predicts farther with cv_rf() on folds at its distances
#>   (?fold_separation).
if (requireNamespace("sp", quietly = TRUE))
  print(fit_gwr_model(pts, "z", "a", bandwidth = 30))   # print.spatial_fit()
#> <GWR> spatial model fit
#>   Formula : z ~ a
#>   n       : 60
#>   CRS     : EPSG:32632
#>   Bandwidth: 30 neighbours (adaptive, bisquare kernel)
#>   AICc    : 47.44
#>   Eff. parameters: 13.2 of 60 observations
```
