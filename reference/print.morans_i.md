# Print a residual Moran's I result

Prints the statistic, the null its moments come from and the evidence
behind them. The weight matrix the result carries is described in one
line instead of printed.

## Usage

``` r
# S3 method for class 'morans_i'
print(x, ...)
```

## Arguments

- x:

  An object of class `morans_i`, from
  [`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md).

- ...:

  Not used.

## Value

`x`, invisibly.

## See also

[`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md),
which makes the object.

Other in-sample checks:
[`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md),
[`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md),
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md),
[`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)

## Examples

``` r
library(sf)
set.seed(1)
dat <- st_as_sf(
  data.frame(x = 5e5 + runif(80, 0, 1000), y = 5e6 + runif(80, 0, 1000),
             elev = rnorm(80)),
  coords = c("x", "y"), crs = 32632
)
dat$price <- 10 + 2 * dat$elev + rnorm(80)
fit <- fit_rf_model(dat, "price", "elev", num_trees = 100, seed = 1)
print(residual_morans_i(fit))
#> Residual Moran's I = 0.0313   (E[I] = -0.0127, sd = 0.0494)
#>   z = 0.891, p = 0.373
#>   null: randomisation moments, approximate for these residuals; n = 80, df = 79
#>   residual kurtosis 2.584 (3 = Gaussian)
#>   weights: 80 x 80 dgCMatrix, 8 neighbour(s) per row; not retained, keep_weights = TRUE to keep it
```
