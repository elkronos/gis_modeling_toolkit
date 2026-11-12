# Plot a GWR model selection

[`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md)
ranks every model it evaluated on AICc. This draws each model's
criterion against its number of predictors, the winner in red, so the
gap between the best model and the runners-up (which the ranked table
shows only as numbers) is read as a shape. A winner well below the rest
is preferred; one within 2 AICc units of others is not separated from
them (see "Printing" in
[`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md)).
The criterion is in-sample, and the caption says so.

## Usage

``` r
# S3 method for class 'gwr_model_selection'
plot(x, ...)
```

## Arguments

- x:

  A `gwr_model_selection` object.

- ...:

  Not used.

## Value

A `ggplot` object.

## See also

[`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md),
which makes the object.

Other predictor selection:
[`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md),
[`plot.feature_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.feature_selection.md),
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)

## Examples

``` r
library(sf)
set.seed(1)
n <- 80
dat <- st_as_sf(
  data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
             a = rnorm(n), b = rnorm(n), noise = rnorm(n)),
  coords = c("x", "y"), crs = 32632)
dat$z <- 2 * dat$a - dat$b + rnorm(n, 0, 0.5)
sel <- gwr_model_selection(dat, "z", c("a", "b", "noise"), bandwidth = 30)
# Every model tried, AICc against its size; the winner (a and b) marked.
plot(sel)
```
