# Print a random forest fit

Shows the forest's shape (formula, n, number of trees, `mtry`, node
size) along with the out-of-bag error and, prominently, whether the
coordinates were used as predictors. That line is the one to check. A
forest fitted without them prints "no spatial information used": it sees
the covariates alone. One fitted with `include_coords = TRUE` fits
location closely, so its out-of-bag error overstates its accuracy and it
has to be scored with blocked or distance-matched folds; see
*Coordinates as predictors* in
[`fit_rf_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md).
When a predictor is an unordered factor or a character column, a line
says how ranger was told to split it.

## Usage

``` r
# S3 method for class 'rf_fit'
print(x, ...)
```

## Arguments

- x:

  An `rf_fit`.

- ...:

  Ignored.

## Value

`x`, invisibly.

## See also

Other methods on a fitted model:
[`coef.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.bayesian_fit.md),
[`coef.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.gwr_fit.md),
[`coef.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.rf_fit.md),
[`fitted.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.bayesian_fit.md),
[`fitted.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.gwr_fit.md),
[`fitted.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.rf_fit.md),
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md),
[`predict.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.bayesian_fit.md),
[`predict.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.gwr_fit.md),
[`predict.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.rf_fit.md),
[`print.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.spatial_fit.md),
[`residuals.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.bayesian_fit.md),
[`residuals.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.gwr_fit.md),
[`residuals.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.rf_fit.md),
[`summary.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.spatial_fit.md)
