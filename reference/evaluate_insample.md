# Compute the metrics of several fitted models

Runs
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
on each fit in a list and binds the rows, with a `model` column and a
`metric_basis` column saying whether each row is in-sample, out-of-bag
or on `newdata`. It takes a single `spatial_fit` too. No model is
refitted: the fitted values are read with
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html), and `newdata`
is scored with [`predict()`](https://rdrr.io/r/stats/predict.html).
[`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md)
adds information criteria and residual Moran's I to the same table.

## Usage

``` r
evaluate_insample(fits, newdata = NULL, ...)
```

## Arguments

- fits:

  A `spatial_fit` object, or a named list of them (e.g.
  `list(GWR = gwr_obj, Bayesian = bayes_obj)`). The names are used as
  the model labels and every element must have one; an unnamed list is
  an error, and so are duplicated names. `model` is the key the
  comparison table is assembled on, so two fits sharing a name cannot be
  told apart in the output.

- newdata:

  Optional sf object for out-of-sample evaluation. Must contain the
  response variable and all predictors. If NULL, in-sample metrics are
  computed.

- ...:

  Passed to [`predict()`](https://rdrr.io/r/stats/predict.html) when
  `newdata` is given. Without `newdata` nothing reads them, so any
  argument here is an error, which catches a misspelt `newdata`.

## Value

A data.frame with one row per model and columns for model name, all
regression metrics, and `metric_basis`: what the row's metrics were
computed on, `"in-sample"` (fitted values), `"out-of-bag"` (an
`rf_fit`'s fitted values, see "What the metrics are computed on") or
`"newdata"`. Rows with different bases do not compare like for like. An
element that is not a `spatial_fit` is skipped, with a logged warning,
and has no row; a list in which no element is a `spatial_fit` is an
error.

## Details

The metrics are those of
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md):
see its sections "What the metrics are computed on" and "Percentage
errors on responses with zeros".

## See also

Other in-sample checks:
[`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md),
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md),
[`print.morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.morans_i.md),
[`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)

## Examples

``` r
library(sf)
set.seed(1)
pts <- st_as_sf(
  data.frame(x = 5e5 + runif(60, 0, 1000), y = 5e6 + runif(60, 0, 1000),
             a = rnorm(60)),
  coords = c("x", "y"), crs = 32632
)
# An offset keeps the response away from zero, where MAPE is undefined.
pts$z <- 10 + 2 * pts$a + rnorm(60, 0, 0.3)
# Fit on 40 rows and keep 20 back, so the second call is out of sample.
fit <- fit_rf_model(pts[1:40, ], "z", "a", num_trees = 50, seed = 1)
evaluate_insample(fit)                           # out-of-bag for a forest
#>    model  n     RMSE       MAE     MAPE    SMAPE        R2 Adj_R2 n_MAPE
#> 1 rf_fit 40 0.409009 0.3236812 3.425689 3.390414 0.9541711     NA     40
#>   n_SMAPE metric_basis
#> 1      40   out-of-bag
evaluate_insample(fit, newdata = pts[41:60, ])   # on the 20 held-out rows
#>    model  n      RMSE       MAE    MAPE    SMAPE        R2 Adj_R2 n_MAPE
#> 1 rf_fit 20 0.5080838 0.3836504 3.56621 3.649174 0.9022065     NA     20
#>   n_SMAPE metric_basis
#> 1      20      newdata
```
