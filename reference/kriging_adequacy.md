# Check whether block kriging would improve cell means

[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
aggregates by plain means inside cell boundaries. A block-kriging
aggregator would instead weight observations by their spatial
correlation with the cell, and would return an estimate and a variance
for every cell, thin or empty. This function measures what block kriging
would add on a given layer, without changing any cell value. For every
cell it reports the block-kriging estimate and variance implied by a
fitted variogram (`kr_pred`, `kr_var`), that variance as a share of the
variance the cell's mean would have with no data at all (`kr_ratio`),
and, where the cell has points, the error variance of its plain mean
under the same variogram (`plain_var`) and the share of that variance
kriging removes (`kr_gain`). It also scores the variogram itself by
blocked cross-validation.

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
  returns. Any Z or M dimension is dropped, on every layer, because
  [`gstat::krige()`](https://r-spatial.github.io/gstat/reference/krige.html)
  uses every coordinate dimension of its locations and the variogram is
  fitted in x and y
  ([`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  drops Z and M too): the points are kriged in x and y. Multipoints
  whose geometries carry M (XYM or XYZM) are reduced to points in x and
  y, since GEOS does not read M.

- response_var:

  The response column: numeric or logical. A logical column is kriged as
  0/1, so its cell mean is a share. Any other class (a factor, an
  ordered factor, character, a date) is an error naming the column and
  its class, and is never converted: a factor's level codes are an
  arbitrary relabelling, and kriging them gave a variogram, kriged
  estimates and plain means that all changed with the order of the
  levels. For a category, pass an indicator (`x == "level"`).

- cells_sf:

  The cell polygons, with the matching ID column. A layer with no CRS is
  taken to be in the points' CRS (and points with none in the cells'),
  with a warning. Cells that do not all have the same dimensions are
  read in x and y, without their Z and M. So are cells whose geometries
  carry M (XYM or XYZM), which GEOS does not read.

- id_col:

  Name of the ID column, found as
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  finds it: on the points the first of `id_col`, `"poly_id"`,
  `"polygon_id"` and `"cell_id"`; on the cells the first of that column,
  `"poly_id"`, `"polygon_id"`, `"id"`, `"cell_id"` and `"grid_id"`. IDs
  are matched as text, whole numbers written out in full, so a double
  `1e5` matches an integer `100000`; a point whose ID matches no cell is
  counted in no cell, with a warning. `NA` and `NaN` are no ID: a cell
  without one holds no point (`n` 0, `mean` and `se` `NA`) and is kriged
  as any other, and a `cells_sf` in which no cell has an ID is an error.
  Two cells with one ID are an error too, since the points that carry it
  cannot be placed in one of them. Default `"poly_id"`.

- sac:

  Optional `sac_range` carrying a variogram model.

- folds:

  Optional
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  result on `assigned_points_sf` for the cross-validation statistic, or
  a vector of fold labels, one per row of `assigned_points_sf` (the
  labels of rows dropped for a missing response or location are not
  used; one label per row kept is read as well); built here with
  `block_kfold` when `NULL`. As in `cv_*()`, folds whose recorded rows
  sit at other locations in `assigned_points_sf` (built on another
  layer, such as the points before
  [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
  dropped some) are refused. So are folds whose test sets share a
  location: each location is kriged once, under one fold.

- k, seed:

  The number of blocked folds built when `folds` is `NULL`, and the
  random seed for that construction and for the variogram estimated when
  `sac` is `NULL`. Defaults 5 and 123.

- nmax:

  The number of neighbours each kriging system uses (gstat's `nmax`):
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

  Logical; suppress this function's progress
  [`message()`](https://rdrr.io/r/base/message.html)s. It does not
  silence R warnings, nor the package's console log echo (see
  [`spatialkit_quiet`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)
  for that). Default `TRUE`.

## Value

An `sf` object of class `"kriging_adequacy"`, one row per cell with the
cell geometry and: the ID column, `n` (points in the cell), `mean` (the
plain mean), `se` (its naive standard error), `plain_var` (the plain
mean's error variance under the variogram), `kr_pred`, `kr_var`,
`kr_ratio`, `kr_gain`, `kr_shift` and `kr_n_used` (how many locations
the cell was kriged from; `NA` for a cell left out). Attributes:
`variogram` (the model frame), `sill`, `nugget`, `range`,
`range_identified`, `rejected_reason` (why the range was refused, from
`sac`; `NA` when it was identified), `cv` (a list: `zscore_var`,
`zscore_mean`, `rmse`, `n_pred`, `k`, `method`), `nmax`, `n_points` (the
points used), `n_locations` (the distinct locations among them, which
the kriging used), `cells_left_out` (counts of cells left out, named
`shape` for `max_box_ratio` and `size` for `max_neighbours`) and
`replicates` (`NULL` when no two points share a location; otherwise a
list with `n_points`, the points that share one, `within_var`, their
pooled within-location variance, `df`, its degrees of freedom, `nugget`,
the variogram's, `nugget_used`, the larger of the two, which is the
error variance a location seen once is kriged with, and `above_nugget`,
`TRUE` when the replicate variance exceeds the variogram's nugget).

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

- `plain_var`:

  The error variance of the cell's plain mean as an estimate of the
  cell's average, under the same variogram and the same discretisation
  of the cell as `kr_var`: \\a' C a - 2 a' \bar c(B) + \bar C(B,B)\\,
  with \\a\\ the plain mean's weights on the cell's observations, \\C\\
  their covariances (the nugget on the diagonal), \\\bar c(B)\\ their
  covariances with the cell and \\\bar C(B,B)\\ the covariance within
  it. It is the model's counterpart of `se^2`. \\s^2/n\\ is the variance
  the plain mean would have if the cell's points were independent draws;
  under spatial correlation it understates the error, most where a
  cell's points sit close together, so `se^2 / plain_var` below 1 says
  the naive standard error is too small if the variogram is right.
  Defined from one point up. `NA` for a cell with no point and for a
  cell left out of the kriging.

- `kr_gain`:

  `1 - kr_var / plain_var`: the share of the plain mean's error variance
  that kriging removes in the cell. Both variances come from the one
  variogram, so the gain does not move when sill and nugget are
  multiplied by a constant (unless locations are repeated, whose
  replicate variance is read off the data and does not scale with the
  model). It is at least 0: the kriging neighbourhood always holds the
  cell's own points, the plain mean is one unbiased combination of them,
  and kriging is the best such combination. The exception is a cell with
  rows at a location whose rows name more than one cell: that location
  is no cell's own (see "The kriging neighbourhood"), and the
  neighbourhood need not hold it. It is what the variogram expects, not
  what happened: read it beside `attr(, "cv")$zscore_var`, which says
  whether the variogram's variances are the right size. `NA` where
  `plain_var` is. `kr_var` is not compared with \\s^2/n\\. That
  comparison sets a variance from the model against one from the cell's
  own spread, so it moves with the fitted sill, and it does not mark the
  cells where kriging fails to pay.

- `kr_shift`:

  The kriged estimate minus the plain mean, over the standard deviation
  the variogram gives that difference, `sqrt(plain_var - kr_var)` (the
  kriging error is uncorrelated with the difference, so the two
  variances subtract). Under the right variogram it has variance 1:
  about 32 percent of cells beyond 1 and 5 percent beyond 2. More than
  that says the plain means and the variogram disagree, so the counts
  [`print()`](https://rdrr.io/r/base/print.html) gives are a check of
  the variogram on the cells themselves, beside the cross-validation
  statistic. In units of the naive standard error \\s/\sqrt{n}\\ the
  same shift does not have variance 1 under the right variogram, so
  counts beyond 1 and 2 in that unit are not this check. `NA` for a cell
  with no point, and where the two variances agree to within a
  ten-thousandth, which is where kriging returns the plain mean.

## The cross-validation statistic

The kriging variance is only as good as the variogram. With blocked
folds, each held-out point is kriged from the training folds and the
standardised error \\(z - \hat z) / \sqrt{kv}\\ is recorded; the
variance of those over all held-out points should be about 1. Its
departure from 1 measures how badly the kriging variance is understated
(or overstated): 1.4 means the variances are roughly 40 percent too
small. Under blocked folds every held-out point is far from the training
data, where the kriging variance is close to the sill whatever the
nugget, so the blocked statistic checks the sill and range and says
little about the nugget. To check the nugget, pass random folds
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
extrapolation; a variogram that falls with distance, or a range below
the shortest lag, describes the data poorly at some lags; and a response
in which
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
found no spatial structure comes with a pure-nugget model, which gives
kriging no correlation to use.

The model has to be of the response itself. A variogram of residuals
(`estimate_sac_range(predictor_vars = ...)`, `attr(, "detrended")`
`TRUE`) leaves out the variance the predictors explain, while the
response is kriged here without them, so `kr_var`, `kr_ratio` and the
cross-validation statistic come out too small; it is used with a
warning. The points are put in the CRS the variogram was fitted in
(`attr(sac, "crs")`), because its range is a length in that CRS's units.
Requires gstat.

## The kriging neighbourhood

gstat kriges a cell from the `nmax` locations nearest its centre. On
gridded locations several of them are often the same distance away, and
which of those the last bits of the coordinates make nearer moves with
the unit and the origin. So the `nmax` are chosen by the package's tie
rule: squared distances within one part in ten million of each other
count as equal, and a tie goes to the location of lower first
coordinate, then lower second, coordinates within that part of the
distance counting as equal. Where the `nmax`-th nearest location ties
with the next, the cell is kriged from the locations the rule picks;
elsewhere the `nmax` nearest are one set and gstat's own search finds
it. Each held-out location of the cross-validation is kriged from its
`nmax` nearest training locations the same way. The kriged figures, and
the cross-validation statistic on the same folds, are then the same in
another unit or with another origin, to rounding, as long as no
coordinate moves by more than about a two-hundred-millionth of the
distance between neighbouring locations. The rows are put in order by
their coordinates, then their cell ID, then their response before any of
this, so the result does not depend on the order of the rows either.

A cell holding more locations than `nmax` would be estimated from its
middle alone, which describes the middle rather than the cell. So
wherever the `nmax` locations nearest a cell's centre leave out any of
the cell's own locations (more than `nmax` locations are as near as its
furthest own one, under the same tolerance), the cell is kriged from all
of its own locations plus the `nmax` nearest outside it, chosen by the
same rule; `kr_n_used` says how many locations each cell was kriged
from. The cost of a kriging system grows with the cube of its size, so a
cell that would need more than `max_neighbours` is left out with a
warning, its `kr_` columns `NA`. A cell's own locations are those whose
rows name it. A location whose rows name more than one cell (rows at one
place carrying different IDs) is no cell's own, with a warning that
gives the count of such locations: it is kriged from as any other
location, and each of its rows still counts in `n`, `mean`, `se` and
`plain_var` of the cell it names.

gstat also discretises each cell into 500 points on a regular grid laid
over the cell's whole bounding box, keeping those inside, so the memory
a cell takes grows with the ratio of that box to its area. A cell whose
bounding box exceeds its area more than `max_box_ratio` times (a sliver,
parts far apart, a cell that is mostly hole) is left out the same way.
`attr(, "cells_left_out")` counts both kinds, and
[`print()`](https://rdrr.io/r/base/print.html) says how many cells have
no estimate and why. An empty cell, or one with a coordinate that is not
finite or lies beyond 1e150 (read as given), has no location and is left
out as well, with a warning of its own, and comes back with an empty
geometry; it is counted with the cells of the second kind (`shape`).

## Repeat measurements at one location

Two observations at the same coordinates (visits to a station, records
geocoded to one address) make a kriging system singular, because gstat
gives them the full sill, nugget included, as their covariance, as if
they were one observation. The kriging and its cross-validation
therefore use one observation per location, the mean of its replicates,
with a warning; `n` and `mean` still count every point. How much of the
nugget \\c_0\\ a mean of \\m\\ replicates keeps is read off the
replicates: their pooled within-location variance \\s_w^2\\ is the part
that differs from visit to visit and averages down, and what the nugget
holds beyond it is micro-scale variation the visits share, so the mean
carries error variance \\\max(c_0 - s_w^2, 0) + s_w^2/m\\. That goes to
gstat as a known measurement error (its `weights`) on the model with the
nugget set to zero. For replicates that differ only by measurement error
this is exactly the kriging of every observation, and for identical
replicates it is the kriging of one. `plain_var` uses the same error
variances: the observations a cell has at one location share the
micro-scale part and average the rest.

When \\s_w^2\\ is above \\c_0\\ the replicates show more variation
between visits than the model has room for (with a zero nugget, any at
all). The replicate variance is then used as it is: a mean of \\m\\ is
kriged with error variance \\s_w^2/m\\, and a location seen once with
\\s_w^2\\ where the model gives it \\c_0\\. Capping it at \\c_0\\ would
krige the means as more exact than their replicates show them to be. The
warning then carries the class `"spatialkit_replicate_var_above_nugget"`
and both numbers, which `attr(, "replicates")` records and
[`print()`](https://rdrr.io/r/base/print.html) repeats. The replicate
variance is an estimate (its degrees of freedom are in the warning and
the record), and the variogram still says \\c_0\\: refit it with a
nugget of at least \\s_w^2\\.

In the cross-validation a location is held out whole, under the fold of
its first row, and its error variance is part of its standardised error.
A cell or held-out location gstat still cannot krige is reported `NA`
with a warning, and [`print()`](https://rdrr.io/r/base/print.html) says
how many.

## References

Cressie, N. (1993). *Statistics for Spatial Data*, revised edition.
Wiley. (Block kriging, chapter 3; cross-validation of the kriging
variance, section 2.6.4.)

## See also

Other assign and summarise:
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md),
[`cell_values_at()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cell_values_at.md),
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)

## Examples

``` r
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
#> Registered S3 method overwritten by 'stars':
#>   method                  from
#>   st_interpolate_aw.stars sf  
print(ka)                   # the report
#> Block-kriging adequacy over 20 cells (200 points, nmax 50)
#>   variogram: Nug(0.184, 0) + Exp(1.23, 149); sill 1.42, nugget 0.184 (13%),
#>     range 445.9
#>   kriging variance / no-data variance of the cell mean (1 = nothing from the
#>     data): median 0.062, range 0.032-0.293; 0 cell(s) above 0.5
#>   kr_gain, the share of the plain mean's error variance that kriging removes
#>     under this variogram (1 - kr_var / plain_var): median 0.51, range 0.21-0.81
#>     over the 19 cell(s) with points and a kriged estimate
#>   kriged minus plain mean, over the standard deviation the variogram gives that
#>     difference (sqrt(plain_var - kr_var)): |shift| > 1 in 6 of 19 cell(s), > 2
#>     in 0 (about 6 and 1 when the variogram is right)
#>   naive s^2/n below plain_var (the naive standard error is too small if the
#>     variogram is right) in 8 of the 16 cell(s) with two or more points; median
#>     s^2/n over plain_var 1.00
#>   empty cells: 1 (kriged estimate and variance available for each)
#>   blocked CV (block_kfold, 4 folds, 200 points): var of standardised error 0.95
#>     (1 = kriging variance correct; above 1 = understated), mean -0.11, RMSE
#>     0.968
# Variance of the standardised errors: about 1 when the variogram is right
attr(ka, "cv")$zscore_var
#> [1] 0.950956
```
