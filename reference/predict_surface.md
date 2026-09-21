# Predict a fitted spatial model onto a regular grid

Builds a prediction surface over the extent of the training data (or
over a grid you supply), predicts in chunks, and returns an `sf` layer.

## Usage

``` r
predict_surface(
  object,
  grid = NULL,
  cell_size = NULL,
  n_cells = 10000L,
  boundary = NULL,
  covariates = NULL,
  chunk_size = 5000L,
  se = FALSE,
  ...
)
```

## Arguments

- object:

  A `spatial_fit` (e.g. from
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  or
  [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)).

- grid:

  Optional `sf` POINT layer to predict onto. When `NULL`, a regular grid
  is built over the training extent. Must have at least one row. It is
  brought into the fit's CRS first: a CRS-less grid is given the
  interpretation the training data got (the assumption recorded on the
  fit), with a warning, and then reprojected. Otherwise a CRS-less grid
  can land thousands of kilometres from the covariates and every cell
  takes the same nearest feature.

- cell_size:

  Grid resolution in CRS units. Ignored when `grid` is supplied; when
  `NULL`, derived from `n_cells`. A value that would produce more than
  5,000,000 cells is refused, naming the implied count and the CRS
  units. The usual cause is a value in the wrong unit. A `cell_size`
  wider than the extent yields a single centred cell.

- n_cells:

  Approximate cell count used to derive `cell_size`. Default 10000. Must
  be a single positive finite number and at most 5,000,000; anything
  else is an error. Also ignored when `grid` is supplied. The grid you
  pass is used verbatim.

- boundary:

  Optional polygonal `sf`/`sfc`; grid points outside it are dropped. Put
  through the same CRS replay and reprojection as `grid`.

- covariates:

  Optional `sf` layer carrying the model's predictors. Required when the
  model has predictors and `grid` does not already contain them. Values
  are taken from the nearest feature.

- chunk_size:

  Rows per prediction call. Default 5000. A pure performance knob for
  the GWR and random-forest backends, whose rows do not interact. For a
  `bayesian_fit` it is also that, *provided* the grid stays inside the
  training extent. Beyond it the GP boundary has to grow and predictions
  depend on which rows share the call; see
  [`predict.bayesian_fit`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.bayesian_fit.md).

- se:

  Logical; also return a standard-error/posterior-SD column where the
  backend supports it. Default FALSE.

- ...:

  Passed to [`predict()`](https://rdrr.io/r/stats/predict.html).

## Value

An `sf` POINT layer with a `.pred` column (and `.pred_se` when
`se = TRUE` and available). For an auto-generated grid the resolution is
attached as attribute `"cell_size"`. For a user-supplied `grid` it is
only whatever `"cell_size"` attribute that object already carried. That
is usually `NULL`, and `NULL` for certain if the grid had to be
re-projected, since
[`st_transform()`](https://r-spatial.github.io/sf/reference/st_transform.html)
does not preserve custom attributes. The resolution of a grid you built
is not this function's to infer.

## Details

[`predict()`](https://rdrr.io/r/stats/predict.html) on a `spatial_fit`
requires `newdata` to be constructed by hand, which makes the most
common downstream task (produce a map) more work than it should be. This
wraps the grid construction, covariate join, chunking and CRS handling.

Prediction over a grid is embarrassingly parallel in the sense that rows
do not interact, so it is chunked: for `bayesian_fit` the posterior draw
matrix is `n_draws x n_newdata`, which will exhaust memory on a fine
grid long before the fit itself would.

## See also

Other prediction:
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)

## Examples

``` r
# Any spatial_fit works here; a forest keeps the example free of the
# optional GWR/Stan backends.
if (requireNamespace("ranger", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  n <- 120
  pts <- st_as_sf(
    data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
               elev = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
  pts$price <- 10 + 0.01 * (st_coordinates(pts)[, 1] - 5e5) +
    2 * pts$elev + rnorm(n)
  fit  <- fit_rf_model(pts, "price", "elev", num_trees = 100, seed = 1)
  surf <- predict_surface(fit, n_cells = 500, covariates = pts)
  print(surf[".pred"])        # one prediction per grid cell, as an sf layer
  # Check where that surface is extrapolating before mapping it.  The grid
  # took its covariates from the nearest observation, so here nothing is
  # outside; a grid with its own covariate raster is where this bites.
  area_of_applicability(surf, model = fit)
}
#> Simple feature collection with 484 features and 1 field
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 500034.7 ymin: 5000057 xmax: 500943.7 ymax: 5000966
#> Projected CRS: WGS 84 / UTM zone 32N
#> First 10 features:
#>       .pred                 geometry
#> 1  11.53496 POINT (500034.7 5000057)
#> 2  11.71447   POINT (500078 5000057)
#> 3  11.71447 POINT (500121.3 5000057)
#> 4  12.06628 POINT (500164.6 5000057)
#> 5  12.06628 POINT (500207.9 5000057)
#> 6  12.06628 POINT (500251.1 5000057)
#> 7  13.07745 POINT (500294.4 5000057)
#> 8  13.07745 POINT (500337.7 5000057)
#> 9  13.07745   POINT (500381 5000057)
#> 10 11.63669 POINT (500424.3 5000057)
#> Area of applicability (Meyer & Pebesma 2021)
#> 
#>   predictors  : 1 (elev)
#>   weighted    : no (all predictors count equally)
#>   training    : 120 points
#>   reference   : nearest other training point (no folds supplied)
#>   normaliser  : 1.1283 (mean pairwise distance)
#>   threshold   : 0.0480 (outlier-removed max of training DI)
#> 
#>   484 of 484 prediction points inside the AOA (100.0%)
```
