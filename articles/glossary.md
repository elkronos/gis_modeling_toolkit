# Terms used on these pages

The articles and help pages use the terms below in these senses.

The terms, A to Z: [area of applicability
(AOA)](#area-of-applicability), [autocorrelation
range](#autocorrelation-range), [backend](#backend), [block](#block),
[cell](#cell), [cell’s true mean](#cells-true-mean), [CRS
units](#crs-units), [design effect (`deff`)](#design-effect),
[dissimilarity index (DI)](#dissimilarity-index), [flat
region](#flat-region), [fold](#fold), [level](#level),
[lon/lat](#lon-lat), [nugget, sill, range](#nugget-sill-range), [pooled
value](#pooled-value), [range floor](#range-floor), [row ID](#row-id),
[size](#size), [support ceiling](#support-ceiling),
[variogram](#variogram).

- area of applicability (AOA):

  The prediction locations where a model’s cross-validated score can be
  expected to hold: those whose dissimilarity index is at or below a
  threshold, the largest dissimilarity index of the training rows that
  is not an outlier.
  [`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
  computes it.

- autocorrelation range:

  The distance beyond which two observations are close to uncorrelated,
  in the units of the CRS the variogram is fitted in.
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  returns it as the effective range of an exponential variogram: three
  times the fitted range parameter, where the semivariance has covered
  about 95 percent of the rise from the nugget to the sill. A block or a
  buffer has to exceed it to keep held-out rows apart from their
  training data.

- backend:

  One of the three model fitters behind the `spatial_fit` class: the
  random forest
  ([`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)),
  GWR (geographically weighted regression,
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md))
  and the Bayesian GP (a Bayesian model with a Gaussian-process term in
  the coordinates,
  [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)).
  [`new_spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md)
  wraps a model of your own in the same class.

- block:

  A piece of the study area that a fold holds out whole: a cell of a
  rectangular grid, at least `block_size` on each side when that is
  given, or one of the polygons passed in `blocks`.
  `make_folds(method = "block_kfold")` builds the grid.

- cell:

  A polygon of a tessellation: a Voronoi cell, a hexagon, a square or a
  Delaunay triangle.
  [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
  puts each observation in one cell, and the `poly_id` column names it.

- cell’s true mean:

  What a cell value estimates: the expected value of one further
  observation located in the cell the way its observations were. It is
  not the average of a surface over the cell’s area, which would give
  weight to ground that holds no observation. “Pooled values” in
  [`?summarize_by_cell`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  defines it.

- CRS units:

  The unit of the coordinate reference system a function works in, such
  as the metre for a UTM zone. Distances, block sizes, bandwidths and
  ranges are given in it. A function that measures distances moves a
  lon/lat layer into a projected CRS that
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  chooses, whose unit is the metre
  ([`?"spatialkit-geometry"`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)).

- design effect (`deff`):

  The factor by which correlation among a cell’s observations inflates
  the variance of the cell mean as an estimate of the population mean;
  `n / deff` is then the cell’s effective sample size.
  `summarize_by_cell(deff = )` applies it to the standard errors. At the
  default, `deff = 1`, the standard error is that of the cell’s own
  mean, which needs no such correction when the cell’s observations are
  spread through it.

- dissimilarity index (DI):

  The distance from a prediction location to the nearest training row in
  the space of the predictors, each centred and scaled by the training
  data’s own mean and standard deviation and weighted by importance when
  the model carries it, divided by the mean distance between training
  rows. Meyer and Pebesma (2021) define it, and
  [`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
  computes it.

- flat region:

  The levels whose criterion value is within `tol` of the optimum, which
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  returns with the count it chooses. The criteria do not rise or fall
  steadily with the count, so the flat region is a set of levels that
  can leave out some in between. Report it with the count.

- fold:

  One split of a cross-validation: the rows held out and the rows the
  model is trained on to predict them.
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  returns the folds, each a `train` and a `test` set of row IDs.

- level:

  A candidate cell count on the ladder that
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  scores.
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  chooses one of the levels.

- lon/lat:

  Coordinates in degrees of longitude and latitude, in a geographic CRS.
  See [CRS units](#crs-units) for how distances are measured on such a
  layer.

- nugget, sill, range:

  The variogram’s value at zero distance (the noise), the value it
  levels off at (the total variance), and the range, the distance over
  which it rises toward the sill. The exponential model
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  fits approaches the sill without reaching it; see [autocorrelation
  range](#autocorrelation-range).
  [`sac_nugget()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md)
  reads the nugget of an
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  result.

- pooled value:

  A cell’s mean pulled toward the mean of the observations in the 8
  nearest cells (nearest by the centroid of their observations). The own
  weight, the share of the pooled value that is the cell’s own mean,
  grows with the cell’s count.
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  returns pooled values beside the plain means by default
  (`estimate = "pooled"`), and “Pooled values” in its help gives the
  formula.

- range floor:

  The fewest cells the autocorrelation range allows: the area of the
  points’ convex hull divided by the square of the range, rounded up,
  and at least 2. A cell wider than the range averages over more than
  one patch of the field.
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  reports the floor and by default does not apply it;
  `select_resolution(criterion = "reliability")` reads the levels from
  the floor up.

- row ID:

  `..row_id`, the column by which folds, cross-validated predictions and
  the record of dropped rows name rows. The `cv_*()` functions number
  the rows in it when a layer has none.
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  names rows by this column when the layer has it, and by row position
  otherwise.

- size:

  The mean number of observations a cell is to hold, `obs_per_cell` in
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md).
  A size `s` gives `max(2, round(n / s))` cells, with `n` the rows that
  carry a response, or every point without one. The default profile
  scores one size, 4.5, and that sets the default count of cells for
  pooled values.

- support ceiling:

  The most cells the data support on a ladder of counts: the scored rows
  (those with a response and every predictor, or every point without a
  response) divided by `min_cell_n`, 9 by default, rounded down. On a
  ladder of sizes the finest level takes its place. When a criterion’s
  optimum sits on the ceiling or on the range floor, that bound sets the
  count, and
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  names it in `edge`.

- variogram:

  The semivariance of a variable against the distance between pairs of
  points, where the semivariance is half the mean squared difference of
  the pairs at that distance.
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  fits an exponential model with a nugget to it, and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on the result
  ([`?plot.sac_range`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.sac_range.md))
  draws the variogram, with the fitted model where a range was
  identified.

## Next

This is the last article in reading order. The help pages that give
these terms in full are
[`?area_of_applicability`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md),
[`?estimate_sac_range`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
[`?make_folds`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md),
[`?resolution_profile`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
[`?select_resolution`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md),
[`?summarize_by_cell`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
and
[`?"spatialkit-geometry"`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md).
[`?spatialkit`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-package.md)
lists every step of the pipeline with its functions.

## References

Meyer, H. and Pebesma, E. (2021). Predicting into unknown space?
Estimating the area of applicability of spatial prediction models.
*Methods in Ecology and Evolution* 12(9), 1620–1633.
<https://doi.org/10.1111/2041-210X.13650>
