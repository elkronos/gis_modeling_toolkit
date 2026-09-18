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
  fit), with a warning, and then reprojected – otherwise a CRS-less grid
  can land thousands of kilometres from the covariates and every cell
  takes the same nearest feature.

- cell_size:

  Grid resolution in CRS units. Ignored when `grid` is supplied; when
  `NULL`, derived from `n_cells`. A value that would produce more than
  5,000,000 cells is refused, naming the implied count and the CRS units
  – the usual cause is a value in the wrong unit. A `cell_size` wider
  than the extent yields a single centred cell.

- n_cells:

  Approximate cell count used to derive `cell_size`. Default 10000. Must
  be a single positive finite number and at most 5,000,000; anything
  else is an error. Also ignored when `grid` is supplied – the grid you
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
  training extent – beyond it the GP boundary has to grow and
  predictions depend on which rows share the call; see
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
only whatever `"cell_size"` attribute that object already carried –
usually `NULL`, and `NULL` for certain if the grid had to be
re-projected, since
[`st_transform()`](https://r-spatial.github.io/sf/reference/st_transform.html)
does not preserve custom attributes. The resolution of a grid you built
is not this function's to infer.

## Details

[`predict()`](https://rdrr.io/r/stats/predict.html) on a `spatial_fit`
requires `newdata` to be constructed by hand, which makes the most
common downstream task – produce a map – more work than it should be.
This wraps the grid construction, covariate join, chunking and CRS
handling.

Prediction over a grid is embarrassingly parallel in the sense that rows
do not interact, so it is chunked: for `bayesian_fit` the posterior draw
matrix is `n_draws x n_newdata`, which will exhaust memory on a fine
grid long before the fit itself would.

## Examples

``` r
# Any spatial_fit works here; a forest keeps the example free of the
# optional GWR/Stan backends.
if (requireNamespace("ranger", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  n <- 120
  pts <- st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), elev = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
  pts$price <- 10 + 0.01 * st_coordinates(pts)[, 1] + 2 * pts$elev + rnorm(n)
  fit  <- fit_rf_model(pts, "price", "elev", num_trees = 100, seed = 1)
  surf <- predict_surface(fit, n_cells = 500, covariates = pts)
  surf[".pred"]
  # Check where that surface is extrapolating before mapping it.
  area_of_applicability(surf, model = fit)
}
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
