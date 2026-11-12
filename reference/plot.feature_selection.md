# Plot the path of a forward feature selection

[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)
scores every candidate at every step and keeps the best; its `history`
holds all of them. This draws the accepted variable's score at each step
as the path, every other candidate's score at that step as a faint
point, and the step at which the selection stopped in red. When it
stopped because no candidate passed the acceptance rule (`rule` in
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)),
the path runs one step further to the best of the rejected candidates,
drawn hollow and labelled "not added", so the stop reads as a flattening
rather than a cut. The picture then shows whether the last variable was
a clear gain or a marginal one, and whether the runner-up would have
done as well. The scores are the selection's own cross-validated
criterion, optimistically biased by the selection (see "Selection score"
in
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md));
when a hold-out score was computed (`select_on = "split"`) it is drawn
as a separate mark at the final step and named in the caption.

## Usage

``` r
# S3 method for class 'feature_selection'
plot(x, ...)
```

## Arguments

- x:

  The list returned by
  [`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md).

- ...:

  Not used.

## Value

A `ggplot` object.

## Details

With a learner that cannot fit an empty predictor set the first variable
is added without a test and can be dropped at the end
(`baseline$first$dropped`). The path drawn is then the one the sweep
walked, with that variable hollow and labelled "dropped at the end", and
the red mark is the selected set, scored without it.

## See also

[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md),
which makes the object.

Other predictor selection:
[`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md),
[`plot.gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.gwr_model_selection.md),
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)

## Examples

``` r
library(sf)
set.seed(4)
n <- 150
dat <- st_as_sf(
  data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
             a = rnorm(n), b = rnorm(n), c = rnorm(n)),
  coords = c("x", "y"), crs = 32632)
dat$z <- 2 * dat$a - dat$b + rnorm(n, 0, 0.5)
fit_fn <- function(train_sf, vars)
  fit_rf_model(train_sf, "z", vars, num_trees = 80, seed = 1)
sel <- select_features_forward(dat, "z", c("a", "b", "c"), fit_fn = fit_fn,
                               k = 3, quiet = TRUE)
plot(sel)
```
