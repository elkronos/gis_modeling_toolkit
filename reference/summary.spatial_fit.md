# Summarise a fitted model

Computes goodness-of-fit metrics from `fitted(object)` against the
observed response, and prints them with what they are to be read
against. A character or factor response is an error; a logical response
is scored as 0/1.

## Usage

``` r
# S3 method for class 'spatial_fit'
summary(object, ...)
```

## Arguments

- object:

  A `spatial_fit` object.

- ...:

  Not used.

## Value

An object of class `summary.spatial_fit`: a list with `class`,
`formula`, `n`, `response_var`, `predictor_vars`, `info` and `in_sample`
(the metric data.frame, out-of-bag for an `rf_fit`).

## What the metrics are computed on

For a GWR or a Bayesian GP fit these are in-sample metrics:
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns values
computed at the training locations from the model that saw them. For a
random forest they are out-of-bag, because
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) on an `rf_fit`
returns out-of-bag predictions (see
[`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)).
The two are not comparable, and
[`print()`](https://rdrr.io/r/base/print.html) on the result says which
one it holds, from `$info$fitted_are_oob`. Use
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
to compare backends.

Adjusted \\R^2\\ is not computed: a GWR's effective number of parameters
far exceeds the global predictor count, and a GP model has no simple
`p`.

For a GWR, [`print()`](https://rdrr.io/r/base/print.html) on the result
adds what the in-sample metrics are to be read against: the effective
number of parameters and residual degrees of freedom, the number of
local windows of deficient rank, the number of locations when
observations share one, and the AICc and RMSE of one global
least-squares regression on the same rows beside the fit's own AICc
(`$info$global_ols`; see
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)).

`MAPE` and `SMAPE` skip rows whose denominator is zero; `n_MAPE` and
`n_SMAPE` count the rows they used (see "Percentage errors on responses
with zeros" in
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)).

## See also

[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
for the metrics alone, and on new data.

Other fit methods:
[`plot.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.spatial_fit.md),
[`print.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.spatial_fit.md),
[`spatial_fit-coef`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-coef.md),
[`spatial_fit-fitted`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-fitted.md),
[`spatial_fit-predict`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-predict.md),
[`spatial_fit-residuals`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-residuals.md)

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
summary(fit)            # out-of-bag metrics for a forest
#> Summary of <rf_fit> fit (n = 60)
#> 
#>   Formula: z ~ a
#> 
#>   Out-of-bag metrics (NOT in-sample; see ?fit_rf_model):
#>     RMSE    = 0.4094
#>     MAE     = 0.3270
#>     R^2     = 0.9493
#>     SMAPE   = 3.33%
```
