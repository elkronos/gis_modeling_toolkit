# Block-kriging adequacy diagnostics for a set of cells

[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
aggregates by plain means inside cell boundaries. A block-kriging
aggregator would instead weight observations by their spatial
correlation with the cell, and would return an estimate and a variance
for every cell, thin or empty. Whether that is worth having on a given
layer is a question with a measurable answer, and this function measures
it, changing no cell value: for every cell it reports the block-kriging
estimate and variance implied by a fitted variogram, that variance as a
share of the variance the cell's mean would have with no data at all,
and, where the cell has points, whether it exceeds the design-based
variance of the plain mean, \\s^2/n\\; and it scores the variogram
itself by blocked cross-validation.

## Usage

``` r
kriging_adequacy(
  assigned_points_sf,
  response_var,
  cells_sf,
  id_col = "poly_id",
  sac = NULL,
  folds = NULL,
  k = 5L,
  seed = 123L,
  nmax = 50L,
  max_neighbours = 2000L,
  max_box_ratio = 1000,
  quiet = TRUE
)
```

## Arguments

- assigned_points_sf:

  Points with a cell identifier column, as
  [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
  returns.

- response_var:

  The response column.

- cells_sf:

  The cell polygons, with the matching ID column. A layer with no CRS is
  taken to be in the points' CRS (and points with none in the cells'),
  with a warning.

- id_col:

  Preferred name of the ID column, found as
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  finds it: on the points the first of `id_col`, `"poly_id"`,
  `"polygon_id"` and `"cell_id"`; on the cells the first of that column,
  `"poly_id"`, `"polygon_id"`, `"id"`, `"cell_id"` and `"grid_id"`. IDs
  are matched as text, whole numbers written out in full, so a double
  `1e5` matches an integer `100000`; a point whose ID matches no cell is
  counted in no cell, with a warning. Default `"poly_id"`.

- sac:

  Optional `sac_range` carrying a variogram model.

- folds:

  Optional
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  result on `assigned_points_sf` for the cross-validation statistic;
  built here with `block_kfold` when `NULL`. As in `cv_*()`, folds whose
  recorded rows sit at other locations in `assigned_points_sf` (built on
  another layer, such as the points before
  [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
  dropped some) are refused.

- k, seed:

  Folds and seed for that construction.

- nmax:

  The number of neighbours each kriging system uses (`gstat`'s `nmax`):
  the locations nearest the cell's centre, or, for a cell those leave
  some of its own locations out of, all of its own plus this many
  outside it (see "The kriging neighbourhood"). The cross-validation
  kriges each held-out point from its `nmax` nearest training locations.
  Default 50.

- max_neighbours:

  The largest kriging system a cell is given when its neighbourhood has
  to grow to hold all its own locations; a cell that would need more is
  left out (`kr_` columns `NA`) with a warning. Never below `nmax`.
  Default 2000.

- max_box_ratio:

  A cell whose bounding box is more than this many times its area is
  left out (`kr_` columns `NA`) with a warning, because gstat's
  discretisation of it costs memory in proportion. Default 1000.

- quiet:

  Suppress progress messages. Default `TRUE`.

## Value

An `sf` object of class `"kriging_adequacy"`, one row per cell with the
cell geometry and: the ID column, `n` (points in the cell), `mean` (the
plain mean), `se` (its naive standard error), `kr_pred`, `kr_var`,
`kr_ratio`, `kr_exceeds_design`, `kr_shift` and `kr_n_used` (how many
locations the cell was kriged from; `NA` for a cell left out).
Attributes: `variogram` (the model frame), `sill`, `nugget`, `range`,
`range_identified`, `rejected_reason` (why the range was refused, from
`sac`; `NA` when it was identified), `cv` (a list: `zscore_var`,
`zscore_mean`, `rmse`, `n_pred`, `k`, `method`), `nmax`, `n_points` (the
points used), `n_locations` (the distinct locations among them, which
the kriging used) and `cells_left_out` (counts of cells left out, named
`shape` for `max_box_ratio` and `size` for `max_neighbours`).

## Reading the columns

- `kr_ratio`:

  The block-kriging variance over the cell's prior variance, in \\\[0,
  1\]\\. The prior variance is the variance the cell's mean would have
  with no data at all, \\\bar C(B,B)\\: the covariance averaged over
  pairs of points in the cell, on the discretisation gstat block-kriges
  with, and without the nugget, which averages out over a block (gstat
  leaves it out of the block variance too). Each cell has its own: a
  cell's mean varies less than a single point does, and far less once
  the cell is wider than the range, so the point sill is not the scale.
  It is the coverage score, and it needs no hand-set threshold in metres
  or point counts: as it approaches 1 the estimate carries almost no
  information from the data about the cell and is reverting to the
  estimated mean. Ordinary kriging adds the variance of that estimated
  mean, so a cell the data do not reach comes out at or above its prior
  variance and reads 1. A cell at 0.05 is well determined; a cell at 0.8
  is mostly prior. `NA` when the model is a pure nugget, where a cell
  mean has no prior variance to be a share of.

- `kr_exceeds_design`:

  `TRUE` where the kriging variance is larger than \\s^2/n\\ from the
  cell's own points: kriging is not earning its keep there, and that is
  said per cell instead of globally. `NA` for cells with fewer than two
  points, where \\s^2\\ does not exist.

- `kr_shift`:

  The kriged estimate minus the plain mean, in units of the plain mean's
  standard error (`NA` where that is not defined). How much the
  aggregator would move the value, against the precision the value has.

## The cross-validation statistic

The kriging variance is only as good as the variogram. With blocked
folds, each held-out point is kriged from the training folds and the
standardised error \\(z - \hat z) / \sqrt{kv}\\ is recorded; the
variance of those over all held-out points should be about 1. Its
departure from 1 measures how badly the kriging variance is understated
(or overstated): 1.4 means the variances are roughly 40 percent too
small. Measured on simulated exponential fields (n = 300 on a 1000-unit
extent, sill 1, nugget 0.2, five blocked folds, eight draws per
configuration): 0.93–1.07 with the true variogram, 0.85–1.01 with the
variogram estimated from the same points. With the nugget understated
tenfold it moved only to 0.95–1.24. That is a property of the folds and
not a weakness of the statistic: under blocked folds every held-out
point is far from the training data, where the kriging variance is close
to the sill whatever the nugget, so the blocked statistic checks the
sill and range. To check the nugget, pass random folds
(`make_folds(method = "random_kfold")`) as `folds`: the held-out points
are then close to their neighbours, where the nugget decides the
variance. The statistic is computed fold by fold with
[`gstat::krige()`](https://r-spatial.github.io/gstat/reference/krige.html)
on the splits
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
built, each held-out point kriged from that split's own training set, so
the folds carry the same separation the package uses everywhere else:
`"buffered_loo"` and `"nndm"` keep the points they exclude around each
held-out one out of its kriging, and
[`print()`](https://rdrr.io/r/base/print.html) names the scheme that
ran. A vector of fold labels is run as k-fold, each fold kriged from all
the others.

## What it said about block kriging as an aggregator

On the same simulated fields, with 16, 36 and 64 square cells: under
uniform sampling the kriged and plain cell means differed by more than
one standard error in 11–24 percent of cells and by more than two in 0–7
percent, the kriging variance exceeded \\s^2/n\\ in 10–26 percent of
populated cells, and at most one cell was empty. Under clustered
sampling (eight clusters of 60-unit spread) the two aggregators parted:
shifts above one standard error in 34–63 percent of cells and above two
in 9–25 percent, the kriging variance below \\s^2/n\\ in only 13–52
percent of populated cells, and 3–27 of the cells empty, each with a
kriged estimate and variance where the plain mean has nothing. So a
block-kriging aggregator earns its place on clustered layers and rarely
on uniform ones, and this function says which kind a layer is.

## What this needs

A variogram model. `sac` is a
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
result carrying one (`attr(, "variogram_model")`); when `NULL` it is
estimated here from the response. The model families are the ones the
package interprets elsewhere: exponential, spherical and Gaussian
components with a nugget. Anything else is refused by name. A model
whose range was not identified (an `NA` estimate with the model
attached) is used with a warning that says why it was refused, and the
reason is kept as `attr(, "rejected_reason")`: a range past the fitted
lags means the sill was never reached and the ratios rest on an
extrapolation; a fit that did not converge stopped wherever the
optimiser halted; a variogram that falls with distance, or a range below
the shortest lag, describes the data poorly at some lags.

The model has to be of the response itself. A variogram of residuals
(`estimate_sac_range(predictor_vars = ...)`, `attr(, "detrended")`
`TRUE`) leaves out the variance the predictors explain, while the
response is kriged here without them, so `kr_var`, `kr_ratio` and the
cross-validation statistic come out too small (the statistic at 4.3–5.2
against 0.67–1.53 with a spatially structured covariate); it is used
with a warning. The points are put in the CRS the variogram was fitted
in (`attr(sac, "crs")`), because its range is a length in that CRS's
units. Requires gstat.

## The kriging neighbourhood

gstat kriges a cell from the `nmax` locations nearest its centre. A cell
holding more locations than that would be estimated from its middle
alone, which describes the middle rather than the cell (in one simulated
case a cell of 1,500 points came out 0.43 off at `nmax = 50`, against
0.07 from all of them and 0.02 for its plain mean). So wherever the
`nmax` locations nearest a cell's centre leave out any of the cell's own
locations, the cell is kriged from all of its own locations plus the
`nmax` nearest outside it; `kr_n_used` says how many locations each cell
was kriged from. The cost of a kriging system grows with the cube of its
size (about 1 s at 2,000 locations and 17 s at 5,000), so a cell that
would need more than `max_neighbours` is left out with a warning, its
`kr_` columns `NA`.

gstat also discretises each cell into 500 points on a regular grid laid
over the cell's whole bounding box, keeping those inside, so the memory
a cell takes grows with the ratio of that box to its area: about 54 MB
more for a thin diagonal strip at a ratio of 708, and 592 MB at 7,072. A
cell whose bounding box exceeds its area more than `max_box_ratio` times
(a sliver, parts far apart, a cell that is mostly hole) is left out the
same way. `attr(, "cells_left_out")` counts both kinds, and
[`print()`](https://rdrr.io/r/base/print.html) says how many cells have
no estimate and why.

## Repeat measurements at one location

Two observations at the same coordinates (visits to a station, records
geocoded to one address) make a kriging system singular, because gstat
gives them the full sill, nugget included, as their covariance, as if
they were one observation. The kriging and its cross-validation
therefore use one observation per location, the mean of its replicates,
with a warning; `n` and `mean` still count every point. How much of the
nugget \\c_0\\ a mean of \\m\\ replicates keeps is read off the
replicates: their pooled within-location variance \\s_w^2\\, capped at
\\c_0\\, is the part that differs from visit to visit and averages down,
and the rest is micro-scale variation the visits share, so the mean
carries error variance \\c_0 - s_w^2 + s_w^2/m\\. That goes to gstat as
a known measurement error (its `weights`) on the model with the nugget
set to zero. For replicates that differ only by measurement error this
is exactly the kriging of every observation, and for identical
replicates it is the kriging of one; a location seen once is kriged as
before. In the cross-validation a location is held out whole, under the
fold of its first row, and its error variance is part of its
standardised error. A cell or held-out location gstat still cannot krige
is reported `NA` with a warning, and
[`print()`](https://rdrr.io/r/base/print.html) says how many.

## References

Cressie, N. (1993). *Statistics for Spatial Data*, revised edition.
Wiley. (Block kriging, chapter 3; cross-validation of the kriging
variance, section 2.6.4.)

## See also

Other aggregation:
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md),
[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md),
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md),
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md),
[`summary.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.resolution_profile.md)

## Examples

``` r
if (requireNamespace("gstat", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  n <- 200
  x <- 5e5 + runif(n, 0, 1000); y <- 5e6 + runif(n, 0, 1000)
  d <- as.matrix(dist(cbind(x, y)))
  z <- as.numeric(t(chol(0.8 * exp(-d / 100) + diag(0.2, n))) %*% rnorm(n))
  pts <- st_as_sf(data.frame(x = x, y = y, z = z), coords = c("x", "y"), crs = 32632)
  bnd <- st_sf(geometry = st_as_sfc(st_bbox(pts)))
  cells <- create_grid_polygons(bnd, target_cells = 16, type = "square")
  asg <- assign_features_to_polygons(pts, cells)
  ka <- kriging_adequacy(asg, "z", cells, k = 4)
  print(ka)                   # the report; print() because only a block's
                              # last value shows on its own
  attr(ka, "cv")$zscore_var   # about 1 when the variogram is right
}
#> Registered S3 method overwritten by 'stars':
#>   method                  from
#>   st_interpolate_aw.stars sf  
#> Block-kriging adequacy over 20 cells (200 points, nmax 50)
#>   variogram: Nug(0.168, 0) + Exp(1.23, 138); sill 1.39, nugget 0.168 (12%), range 414.6
#>   kriging variance / no-data variance of the cell mean (1 = nothing from the data): median 0.066, range 0.034-0.309; 0 cell(s) above 0.5
#>   kriging variance exceeds s^2/n in 0 of the 16 cell(s) with two or more points
#>   kriged minus plain mean: |shift| > 1 SE in 2 of 16 cell(s), > 2 SE in 0
#>   empty cells: 1 (kriged estimate and variance available for each)
#>   blocked CV (block_kfold, 4 folds, 200 points): var of standardised error 0.95 (1 = kriging variance correct; above 1 = understated), mean -0.12, RMSE 0.972
#> [1] 0.947184
```
