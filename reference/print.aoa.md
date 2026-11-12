# Print an area-of-applicability result

Prints how many prediction locations fall inside the area of
applicability and how many outside, the dissimilarity threshold that
separated them, and the predictors the index was computed over (naming
any dropped for having no usable variance). It also says where the
weights came from (the model's importance, a vector the caller supplied,
or equal weights and why), whether location was measured and how, how
categorical predictors were handled, and how many training rows repeat
another's predictor values. With a location index the rows outside are
split into those outside on the predictors and those outside on
location. The share of locations outside is the number to report: a map
that extrapolates over much of its extent rests on predictions the
training data cannot support.

## Usage

``` r
# S3 method for class 'aoa'
print(x, ...)
```

## Arguments

- x:

  An `aoa` object.

- ...:

  Not used.

## Value

`x`, invisibly.

## See also

[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md),
which makes the object, and
[`plot.aoa()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.aoa.md).

Other prediction:
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md),
[`plot.aoa()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.aoa.md),
[`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md)

## Examples

``` r
library(sf)
set.seed(1)
train <- st_as_sf(
  data.frame(x = 5e5 + runif(60, 0, 1000), y = 5e6 + runif(60, 0, 1000),
             a = rnorm(60)),
  coords = c("x", "y"), crs = 32632
)
# 15 new points take values of a seen in training; five lie far above
# its range, and print() counts those five outside.
newpts <- st_as_sf(
  data.frame(x = 5e5 + runif(20, 0, 1000), y = 5e6 + runif(20, 0, 1000),
             a = c(train$a[1:15], rnorm(5, 8))),
  coords = c("x", "y"), crs = 32632
)
res <- area_of_applicability(newpts, train_sf = train, predictor_vars = "a")
print(res)
#> Area of applicability (Meyer & Pebesma 2021)
#> 
#>   predictors  : 1 (a)
#>   weighted    : no (all predictors count equally: no model was given)
#>   training    : 60 points
#>   reference   : nearest other training point (no folds supplied)
#>   normaliser  : 1.1346 (mean pairwise distance)
#>   threshold   : 0.0487 (outlier-removed max of training DI)
#>   location    : not measured
#> 
#>   15 of 20 prediction points inside the AOA (75.0%)
#> 
#> Predictions outside the AOA are extrapolations; the cross-validated
#> performance estimate does not cover them.
```
