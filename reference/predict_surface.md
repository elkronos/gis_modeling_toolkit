# Predict a fitted model onto a regular grid

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

  A `spatial_fit` (from
  [`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md),
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  or
  [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)).

- grid:

  Optional `sf` POINT layer to predict onto, with at least one row. When
  `NULL`, a regular grid is built over the training extent. A grid of
  polygons is reduced to points; see "Coordinate reference systems and
  geometry".

- cell_size:

  Grid resolution in CRS units. Ignored when `grid` is supplied; when
  `NULL`, derived from `n_cells`. A value that would produce more than
  5,000,000 cells is refused, naming the implied count and the CRS
  units. The usual cause is a value in the wrong unit. The grid is
  centred on the training bounding box and covers it: when the extent is
  not a whole number of cells, it overhangs the box by less than one
  cell, split evenly between the two sides. A `cell_size` wider than the
  extent yields a single centred cell.

- n_cells:

  Approximate cell count used to derive `cell_size`. Default 10000. Must
  be a single positive finite number and at most 5,000,000; anything
  else is an error. Also ignored when `grid` is supplied. The grid you
  pass is used verbatim.

- boundary:

  Optional polygonal `sf`/`sfc`; grid points outside it are dropped.
  Aligned to the fit's CRS as `grid` is.

- covariates:

  Optional `sf` layer carrying the model's predictors. Required when the
  model has predictors and `grid` does not already contain them. Values
  are taken from the nearest feature (see "Ties between nearest
  covariate features"). Aligned to the fit's CRS as `grid` is, with the
  same warning when it has no CRS.

- chunk_size:

  Rows per prediction call. Default 5000. A performance setting only:
  rows do not interact, and for a `bayesian_fit` the GP boundary is held
  at its fitted value whatever the chunk holds (see
  [spatial_fit-predict](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-predict.md)).
  The two warnings in which
  [`predict()`](https://rdrr.io/r/stats/predict.html) on a GWR fit
  counts locations (classes `"spatialkit_gwr_extrapolation"` and
  `"spatialkit_gwr_not_estimable"`) are raised once for the grid, with
  its totals, whatever the number of chunks.

- se:

  Logical; also return a standard error (or posterior SD) column where
  the backend supports it. Default FALSE. For a `bayesian_fit` this is
  the SD of the posterior draws
  [`predict()`](https://rdrr.io/r/stats/predict.html) returns, and those
  are of the expected value by default (`type = "epred"`): the
  uncertainty of the mean surface, not of a new observation, which also
  carries the observation noise. For the predictive SD, the one that
  goes with prediction intervals and
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)'s
  calibration, pass `type = "predict"` as well.

- ...:

  Passed to [`predict()`](https://rdrr.io/r/stats/predict.html), e.g.
  `type = "predict"` for a `bayesian_fit`. Not `draws`, which this
  function sets itself and refuses here.

## Value

An `sf` POINT layer with a `.pred` column, and `.pred_se` when
`se = TRUE` and the backend gives one. A `.pred_se` column that a
supplied `grid` carries from an earlier surface is removed when this
call adds none. For an auto-generated grid the resolution is attached as
attribute `"cell_size"`. For a user-supplied `grid` it is only whatever
`"cell_size"` attribute that object already carried. That is usually
`NULL`, and `NULL` for certain if the grid had to be re-projected, since
[`st_transform()`](https://r-spatial.github.io/sf/reference/st_transform.html)
does not preserve custom attributes.

## Details

`predict_surface()` builds the grid, joins the covariates, aligns the
CRSs and calls [`predict()`](https://rdrr.io/r/stats/predict.html) in
chunks of `chunk_size` rows. Rows do not interact, so chunking changes
memory use and not the predictions. For a `bayesian_fit` the draw matrix
is `n_draws` by the number of rows, which on a fine grid outgrows memory
long before the fit does.

## Coordinate reference systems and geometry

`grid`, `boundary` and `covariates` are brought into the fit's CRS
first. A layer without a CRS is read the way the training data were (the
assumption recorded on the fit), with a warning, and then reprojected;
when that reading is longitude and latitude, a layer with an x outside
-180 to 360 or a y beyond 90 north or south cannot be in it, and the
call stops. Without this a CRS-less grid could land thousands of
kilometres from the covariates, and every cell would take the same
nearest feature. A layer still without a CRS after that is taken as
EPSG:4326 and reprojected when its coordinates look like longitude and
latitude, and otherwise stamped with the fit's CRS, with a warning
either way; a fit's CRS in longitude and latitude is not stamped on
coordinates outside longitude -180 to 360 and latitude -90 to 90, which
is refused. See
[spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)
for these rules in general.

A grid of polygons
([`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md)
output, say) is reduced to one representative point per cell, as
[`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md)
does, so covariates are taken at the location predicted for; `boundary`
then keeps the cells whose point falls inside it. All three layers are
read in x and y: one whose geometries mix XY with XYZ or XYM, or carry M
(XYM or XYZM), is read without its Z and M, with no warning for a
boundary or covariates. A POINT grid whose coordinates are not moved
(one in the fit's CRS, or one without a CRS that is stamped with it) is
returned with its rows as given, Z and M included; one that mixes
dimensions and is reprojected into the fit's CRS comes back in XY with a
warning of class `"spatialkit_zm_dropped"`. The points a grid of
polygons is reduced to are in XY when its geometries mix dimensions or
carry M.

A grid row without a location (an empty geometry, or a coordinate that
is not finite or lies beyond 1e150) gets an `NA` prediction; in a grid
that is reprojected, or reduced from polygons, it comes back as an empty
point. A covariate row without a location is never the nearest feature.

## Ties between nearest covariate features

A grid point the same distance from several features of `covariates`
takes one of them by a fixed rule, so a grid laid between the points of
a lattice takes the same covariates in another unit or with another
origin, and with or without the FNN package.

- Points. Distances count as equal when the squared distance is no more
  than the square of the least distance lengthened by a margin, plus one
  part in ten million of the least squared distance. Among them the
  point lowest in x is taken, then the one lowest in y, then the first
  in row order. x values count as one when they lie within one part in
  ten million of the diagonal of the points' bounding box, or within the
  margin where that is more; y values are read the same way. A grid
  point near two points that lie within the margin of each other (two
  records of one station that differ only in the last bits of their
  coordinates) takes the same one however near it lies.

- Lines and polygons, or a layer that mixes them with points (on the
  edge two polygons share, say). Distances within one part in ten
  million of the least, plus the margin, count as the same. Among them
  the feature whose bounding-box centre is lowest in x is taken, then
  lowest in y, x values and y values within one part in ten million of
  the diagonal of the bounding box of `covariates` counting as one, then
  the first in row order.

- A layer in degrees. The nearest feature is
  [`sf::st_nearest_feature()`](https://r-spatial.github.io/sf/reference/st_nearest_feature.html)'s,
  which settles a tie in the order its search meets the features (GEOS,
  or s2 while
  [`sf::sf_use_s2()`](https://r-spatial.github.io/sf/reference/s2.html)
  is on).

- The margin. It is one part in 1e10 of the largest coordinate of the
  features that have a location, at most a millionth of the diagonal of
  their bounding box and at least four times the spacing of doubles at
  that coordinate. It grows with the largest coordinate: about 0.46 mm
  at a UTM northing near 4,550 km, about half a micrometre once the
  origin lies next to data 5 km across. So a grid point whose distances
  to two features differ by less than the margin with one origin and by
  more with another takes the one the tie rule picks with the first and
  the nearer with the second.

## See also

[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
to mark where the surface applies, and
[spatial_fit-predict](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-predict.md)
for the predictions of each backend.

Other prediction:
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md),
[`plot.aoa()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.aoa.md),
[`print.aoa()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.aoa.md)

## Examples

``` r
# Any spatial_fit works here; a forest keeps the example free of the
# optional GWR and Stan backends.
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
surf[".pred"]               # one prediction per grid cell, as an sf layer
#> Simple feature collection with 529 features and 1 field
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 500026.7 ymin: 5000038 xmax: 500979 ymax: 5000990
#> Projected CRS: WGS 84 / UTM zone 32N
#> First 10 features:
#>       .pred                 geometry
#> 1  11.53496 POINT (500026.7 5000038)
#> 2  11.71447   POINT (500070 5000038)
#> 3  11.71447 POINT (500113.3 5000038)
#> 4  12.06628 POINT (500156.6 5000038)
#> 5  12.06628 POINT (500199.9 5000038)
#> 6  12.06628 POINT (500243.2 5000038)
#> 7  13.07745 POINT (500286.5 5000038)
#> 8  13.07745 POINT (500329.7 5000038)
#> 9  13.07745   POINT (500373 5000038)
#> 10 11.63669 POINT (500416.3 5000038)
# Check where the surface extrapolates before mapping it.  The grid took
# its covariates from the nearest observation, so here nothing is
# outside; covariates from another source, such as a raster converted to
# sf points, are where this bites.
area_of_applicability(surf, model = fit)
#> Area of applicability (Meyer & Pebesma 2021)
#> 
#>   predictors  : 1 (elev)
#>   weighted    : yes, by the model's permutation importance, floored at 0
#>   training    : 120 points
#>   reference   : nearest other training point (no folds supplied)
#>   normaliser  : 1.1283 (mean pairwise distance)
#>   threshold   : 0.0480 (outlier-removed max of training DI)
#>   location    : not measured
#> 
#>   529 of 529 prediction points inside the AOA (100.0%)
```
