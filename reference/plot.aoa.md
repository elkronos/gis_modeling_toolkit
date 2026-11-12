# Plot the dissimilarity distribution behind an area of applicability

Draws the dissimilarity index of the prediction locations against that
of the training data, with the threshold marked, so you can see whether
the prediction set sits well inside, near the threshold or outside.
`n_outside` gives only how many prediction locations fall outside the
area of applicability, not how close to the threshold the rest sit or
how far outside the outsiders are. The training DI is cross-validated
over the `folds` passed to
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md),
or without them is each training point's distance to its nearest other
training point; the legend and the caption say which. The training curve
is the reference the threshold was derived from: the threshold is the
largest training DI inside an outlier fence, so the curve reaches it
exactly when no training value was fenced off and runs past it, by the
tail the fence removed, when some were. A prediction location outside on
a predictor dropped for having no training variance, or under
`categorical = "indicator"` on a level no training row has, has
`DI = Inf`: it counts in the prediction curve, which then tops out below
1, and the caption says how many are off the axis. A location with a
missing predictor (`DI = NA`) is neither inside nor outside; it is left
out of the curve and the caption counts it.

## Usage

``` r
# S3 method for class 'aoa'
plot(x, type = c("ecdf", "histogram"), ...)
```

## Arguments

- x:

  An `aoa` object from
  [`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md).

- type:

  `"ecdf"` (default), the two empirical distribution functions on one
  axis, or `"histogram"`, the prediction DI as bars with the training DI
  as an outline.

- ...:

  Not used.

## Value

A `ggplot` object.

## Details

When location has an index of its own (a `gwr_fit`, a `bayesian_fit`, a
forest fitted with `include_coords = TRUE`, or `coords = "yes"`; see
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)),
the curves and the dashed line are the predictor index alone (for the
forest, its `..x` and `..y` columns included). The subtitle then counts
three things: the locations beyond the predictor threshold, those beyond
the location threshold, and those outside the area of applicability,
which is either.

## See also

[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md),
which makes the object, and
[`print.aoa()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.aoa.md).

Other prediction:
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md),
[`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md),
[`print.aoa()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.aoa.md)

## Examples

``` r
library(sf)
set.seed(2)
n <- 200
train <- st_as_sf(
  data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
             a = rnorm(n), b = rnorm(n)),
  coords = c("x", "y"), crs = 32632)
train$z <- train$a - train$b + rnorm(n, 0, 0.3)
# Prediction locations whose predictor `a` drifts beyond the training range.
new <- st_as_sf(
  data.frame(x = 5e5 + runif(100, 0, 1000), y = 5e6 + runif(100, 0, 1000),
             a = rnorm(100, mean = 2), b = rnorm(100)),
  coords = c("x", "y"), crs = 32632)
aoa <- area_of_applicability(new, train_sf = train, predictor_vars = c("a", "b"))
plot(aoa)                     # the two ECDFs

plot(aoa, type = "histogram")
```
