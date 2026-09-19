# Block-kriging adequacy diagnostics for a set of cells

[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
aggregates by plain means inside cell boundaries. A block-kriging
aggregator would instead weight observations by their spatial
correlation with the cell, and would return an estimate and a variance
for every cell, thin or empty. Whether that is worth having on a given
layer is a question with a measurable answer, and this function measures
it, changing no cell value: for every cell it reports the block-kriging
estimate and variance implied by a fitted variogram, that variance as a
share of the total sill, and, where the cell has points, whether it
exceeds the design-based variance of the plain mean, \\s^2/n\\; and it
scores the variogram itself by blocked cross-validation.

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

  The cell polygons, with the matching ID column.

- id_col:

  Preferred name of the ID column. Default `"poly_id"`.

- sac:

  Optional `sac_range` carrying a variogram model.

- folds:

  Optional
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  result on `assigned_points_sf` for the cross-validation statistic;
  built here with `block_kfold` when `NULL`.

- k, seed:

  Folds and seed for that construction.

- nmax:

  The largest number of neighbours each kriging system uses (`gstat`'s
  `nmax`). Default 50.

- quiet:

  Suppress progress messages. Default `TRUE`.

## Value

An `sf` object of class `"kriging_adequacy"`, one row per cell with the
cell geometry and: the ID column, `n` (points in the cell), `mean` (the
plain mean), `se` (its naive standard error), `kr_pred`, `kr_var`,
`kr_ratio`, `kr_exceeds_design` and `kr_shift`. Attributes: `variogram`
(the model frame), `sill`, `nugget`, `range`, `range_identified`, `cv`
(a list: `zscore_var`, `zscore_mean`, `rmse`, `n_pred`, `k`, `method`),
`nmax` and `n_points`.

## Reading the columns

- `kr_ratio`:

  The block-kriging variance over the total sill, in \\\[0, 1\]\\. It is
  the coverage score, and it needs no hand-set threshold in metres or
  point counts: as it approaches 1 the estimate carries almost no
  information from the data and is reverting to the global mean. A cell
  at 0.05 is well determined; a cell at 0.8 is mostly prior.

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
variance.
[`gstat::krige.cv()`](https://r-spatial.github.io/gstat/reference/krige.cv.html)
computes the statistic on the fold labels
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
built, so the folds carry the same separation the package uses
everywhere else.

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
whose range was not identified (a bare `NA` estimate with the model
attached) is used with a warning: its sill was never reached by the
data, so the ratios rest on an extrapolation. Requires gstat.

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
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)

## Examples

``` r
if (requireNamespace("gstat", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  n <- 200
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- as.matrix(dist(cbind(x, y)))
  z <- as.numeric(t(chol(0.8 * exp(-d / 100) + diag(0.2, n))) %*% rnorm(n))
  pts <- st_as_sf(data.frame(x = x, y = y, z = z), coords = c("x", "y"), crs = 32632)
  bnd <- st_sf(geometry = st_as_sfc(st_bbox(pts)))
  cells <- create_grid_polygons(bnd, target_cells = 16, type = "square")
  asg <- assign_features_to_polygons(pts, cells)
  ka <- kriging_adequacy(asg, "z", cells, k = 4)
  ka
  attr(ka, "cv")$zscore_var   # about 1 when the variogram is right
}
#> Registered S3 method overwritten by 'stars':
#>   method                  from
#>   st_interpolate_aw.stars sf  
#> [1] 0.947184
```
