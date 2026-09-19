# Score every candidate number of cells on several criteria at once

A tessellation's resolution is the number of cells the points are cut
into, and no single number settles it: the within-cluster sum of squares
of the coordinates has an elbow, a piecewise-constant approximation of
the response has a Mallows \\C_p\\, the residual autocorrelation of the
cell means has a \\z\\, and the cell means have a reliability. This
function computes all of them at every level of a ladder the data bound,
and returns the table, so that the level chosen, by
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
or by eye, can be defended with the whole profile in place of one
criterion's argmin. Mallows (1973) presents \\C_p\\ itself as a display
of the bias-variance trade-off across candidates; he does not offer it
as a rule that picks one. This is that display, with the other criteria
alongside.

## Usage

``` r
resolution_profile(
  data_sf,
  response_var = NULL,
  predictor_vars = NULL,
  levels = NULL,
  n_levels = 20L,
  min_cell_n = 9L,
  sample_n = 1500L,
  nstart = 25L,
  seed = 123L,
  sac = NULL,
  quiet = TRUE,
  select_on = c("all", "split")
)
```

## Arguments

- data_sf:

  An sf object of points (other geometries are reduced to representative
  points).

- response_var:

  Optional response column name (numeric or logical). Enables `cp` and
  `moran_z`. A variogram estimated from it also sets the floor of the
  ladder and `reliability`.

- predictor_vars:

  Optional predictor column names (numeric or logical). With them, `cp`
  scores the OLS residuals of the response on the predictors, the
  variogram is estimated from those residuals, and `moran_z` regresses
  the cell means on the cell-mean predictors.

- levels:

  Optional integer vector of level counts to score, replacing the
  ladder; values outside `[2, n - 1]` are dropped.

- n_levels:

  Number of levels on the ladder. Default 20.

- min_cell_n:

  Minimum average number of points per cell that a level must keep; sets
  the ceiling. Default 9.

- sample_n:

  Points are subsampled to this many before anything is fitted, as in
  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md).
  Default 1500. The support columns describe the subsample.

- nstart:

  k-means++ restarts per level. Default 25.

- seed:

  RNG seed for the subsample and the restarts; restored afterwards.
  Default 123.

- sac:

  Optional `sac_range` object from
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  to take the range, nugget and correlation function from. Pass one
  fitted with `detrend = "reml"`, say, or on a residual field of your
  choosing. When `NULL` and a response is given, one is estimated on the
  subsample with the same `predictor_vars`.

- quiet:

  Logical; suppress this function's progress
  [`message()`](https://rdrr.io/r/base/message.html)s. Default `TRUE`.

- select_on:

  `"all"` (default) profiles every point; `"split"` profiles one
  spatially blocked half and returns the other half as the set to
  estimate on, in the `"split"` attribute. See the "Post-selection
  inference" section of
  [`determine_optimal_levels`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md);
  the profile reads the response whenever `response_var` is given.

## Value

A data.frame of class `resolution_profile` with one row per level and
columns `levels`, `wss`, `wss_spread` (relative spread of WSS across the
restarts), `elbow`, `cell_n_min`, `cell_n_median`, `cell_diam_median`
(twice the median RMS radius of the cells, in coordinate units), `rss`,
`cp`, `moran_i`, `moran_z` and `reliability`; columns a missing input
leaves undefined are `NA`. Attributes: `bounds` (a list with `floor`,
`ceiling`, `supported`, `area`, `range`, `n`, `min_cell_n`), `variogram`
(a list with `nugget`, `psill`, `range`, `model`; `NULL` when none was
usable), `variable` (`"response"`, `"residuals"` or `NA`), `wss_bumps`,
`nstart`, `sac` (the range object used) and, with `select_on = "split"`,
`split` (a list with `selection` and `estimation`, integer row positions
in `data_sf`).

## The ladder and its bounds

Cells are k-means clusters of the (projected) coordinates, fitted at
each level as the best of `nstart` k-means++ restarts (see
[`determine_optimal_levels`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
for why). Levels are spaced logarithmically, because cell diameter
scales as \\L^{-1/2}\\: a unit step wastes fits at large \\L\\ and
starves resolution at small. The ladder runs from a floor to a ceiling
the data impose. The ceiling is `floor(n / min_cell_n)`: cells with
fewer than `min_cell_n` points on average have too little support, and
the model-aware criteria are not computable below nine cells in any
case. The floor is `ceiling(area / range^2)` when an autocorrelation
range is available: cells wider than the range average over more than
one patch of the field. When the floor exceeds the ceiling the data
cannot support a tessellation that respects their own correlation
structure; that is reported as a finding (a logged warning, and
`attr(x, "bounds")$supported` is `FALSE`) and the ladder runs from 2 to
the ceiling anyway, so the profile still shows what each level costs.

## The criteria, and how each behaved when measured

- `elbow`:

  The signed distance of the WSS curve below the chord from its first to
  its last level, the classical elbow statistic (larger is better).
  Geometry only; it knows nothing of the response.

- `cp`:

  Mallows' \\C_p\\ of the piecewise-constant approximation of the
  response (or of its OLS residuals on `predictor_vars`) by cell means:
  \\RSS(L)/n + 2 \tau^2 L / n\\, with \\\tau^2\\ the nugget of the
  fitted variogram (lower is better). **Measured on simulated
  exponential fields (600 points on a 1000-unit extent, sill 1, 20
  replicates): with a nugget of 0.3 its minimum sat at the support
  ceiling in every replicate at effective ranges of 90, 300 and 900;
  with a nugget of 2 it was interior (median 44, range 8–66).** On a
  smooth field with little noise the approximation keeps improving as
  cells shrink and the penalty is too small to stop it, so \\C_p\\ says
  "as fine as the support allows" and `min_cell_n` is what is choosing;
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  says so when that happens. It becomes a genuine interior criterion
  only when the nugget is a large share of the sill.

- `moran_z`:

  The standardised deviate of Moran's I on the residuals of the cell
  means regressed on the cell-mean predictors (an intercept alone when
  there are none): how much spatial structure the tessellation has left
  unexplained (\\\|z\|\\ smaller is better). Calibrated and flat in
  \\L\\ on a response with no structure (see
  [`determine_optimal_levels`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)),
  so it separates levels only where structure remains. `NA` at nine
  cells or fewer.

- `reliability`:

  The between-cell signal's share of the spread in the cell means, from
  the fitted variogram alone via Krige's additivity relation (Cressie
  1996), for square cells of the level's average area with the level's
  average point count (larger is better). This is the shrinkage factor
  of Fay and Herriot (1979). It has an interior optimum, and a broad
  one: validated against the empirical reliability of true block means
  on simulated fields, the analytic and empirical optima agreed to
  within a level or two where the empirical estimate was stable, and the
  band within 2 percent of the maximum spanned a factor of 3–6 in \\L\\.
  Read the flat region, not the argmax. `NA` without a usable variogram.

`cp` and `reliability` answer different questions: how well the cells
represent the field, and whether the cell values are distinguishable
from noise. They can disagree, and both are shown so the choice between
them is made knowingly.

## References

Cressie, N. (1996). Change of support and the modifiable areal unit
problem. *Geographical Systems*, 3(2–3), 159–180.

Fay, R. E. and Herriot, R. A. (1979). Estimates of income for small
places: an application of James-Stein procedures to census data.
*Journal of the American Statistical Association*, 74(366), 269–277.
[doi:10.1080/01621459.1979.10482505](https://doi.org/10.1080/01621459.1979.10482505)

Mallows, C. L. (1973). Some comments on \\C_p\\. *Technometrics*, 15(4),
661–675.
[doi:10.1080/00401706.1973.10489103](https://doi.org/10.1080/00401706.1973.10489103)

## See also

[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
to read a level and its flat region off the profile;
[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
for the integer-vector interface;
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
and
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
which accept the profile or a
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
result directly as `approx_n_cells` and `n`.

Other aggregation:
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md),
[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md),
[`kriging_adequacy()`](https://elkronos.github.io/gis_modeling_toolkit/reference/kriging_adequacy.md),
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md),
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)

## Examples

``` r
if (requireNamespace("gstat", quietly = TRUE)) {
  library(sf)
  set.seed(2)
  n <- 400
  xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
  D  <- as.matrix(dist(xy))
  xy$z <- as.numeric(t(chol(exp(-D / 100) + diag(0.3, n))) %*% rnorm(n))
  pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
  prof <- resolution_profile(pts, response_var = "z", n_levels = 12)
  prof
  select_resolution(prof, criterion = "reliability")
}
#> Resolution by reliability: 23 cells
#>   flat region : 23 to 31 (6 of 12 levels)
#>   note        : the optimum is the first level of the ladder; the floor
#>                is choosing, not the criterion.
```
