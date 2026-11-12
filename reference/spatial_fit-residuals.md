# Extract the residuals of a fitted model

[`residuals()`](https://rdrr.io/r/stats/residuals.html) on a random
forest, a GWR or a Bayesian GP fit returns the observed response minus
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html), one value per
training row:

- for a random forest these are out-of-bag residuals, since its fitted
  values are out-of-bag predictions. Out-of-bag is a random hold-out:
  use
  [`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
  on folds at the distances the map will be used at for the forest's
  accuracy;

- for a GWR they are in-sample. Where
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  took the response about its centre (a response at a large offset), the
  response and the fitted values are both taken about that centre, so
  the residuals keep digits that rounding at the response's size would
  lose (see "The local solve" in
  [`?fit_gwr_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md));

- for a Bayesian GP they are in-sample, from the cached posterior
  expectation
  ([spatial_fit-fitted](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-fitted.md)).

[`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)
tests them for spatial structure the model did not capture.

## Usage

``` r
# S3 method for class 'gwr_fit'
residuals(object, ...)

# S3 method for class 'bayesian_fit'
residuals(object, ...)

# S3 method for class 'rf_fit'
residuals(object, ...)
```

## Arguments

- object:

  A fitted `rf_fit`, `gwr_fit` or `bayesian_fit`.

- ...:

  Not used.

## Value

A numeric vector of length `object$n`.

## See also

[`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md).

Other fit methods:
[`plot.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.spatial_fit.md),
[`print.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.spatial_fit.md),
[`spatial_fit-coef`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-coef.md),
[`spatial_fit-fitted`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-fitted.md),
[`spatial_fit-predict`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-predict.md),
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
rf <- fit_rf_model(pts, "z", "a", num_trees = 50, seed = 1)
summary(residuals(rf))   # out-of-bag residuals
#>       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#> -0.9261726 -0.2510723  0.0065013  0.0008667  0.3091709  0.9810524 
```
