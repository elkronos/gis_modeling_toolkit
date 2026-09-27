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
  select_on = c("all", "split"),
  range_floor = TRUE
)
```

## Arguments

- data_sf:

  An sf object of points (other geometries are reduced to representative
  points). Features with empty or non-finite coordinates are dropped
  with a warning.

- response_var:

  Optional response column name (numeric or logical). Enables `cp` and
  `moran_z`. A variogram estimated from it also sets the floor of the
  ladder and `reliability`. Rows where it, or a predictor, is missing or
  non-finite stay in the geometry and are left out of the OLS fit, the
  RSS, `cp` and `moran_z`; a logged warning gives their number.

- predictor_vars:

  Optional predictor column names (numeric or logical). With them, `cp`
  scores the OLS residuals of the response on the predictors, the
  variogram is estimated from those residuals, and `moran_z` regresses
  the cell means on the cell-mean predictors.

- levels:

  Optional integer vector of level counts to score, replacing the
  ladder; values below 2, above the number of distinct locations, or at
  or above the number of points, are dropped (k-means cannot place more
  centres than there are distinct points, and
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) refuses as
  many centres as points).

- n_levels:

  Number of levels on the ladder. Default 20.

- min_cell_n:

  Minimum average number of points per cell that a level must keep; sets
  the ceiling. Default 9.

- sample_n:

  Points are subsampled to this many before the k-means fits, as in
  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md).
  Default 1500. The columns read off the fitted cells (`wss`, the
  `cell_` columns, `rss`, `moran_i`, `moran_z`) describe the subsample;
  the bounds, `supported`, the variance term of `cp` and `reliability`
  describe every point of the layer, so the answer does not change with
  `sample_n` except through the fits.

- nstart:

  k-means++ restarts per level. Default 25.

- seed:

  RNG seed for the subsample and the restarts; restored afterwards.
  Default 123. The rows are put in coordinate order before either, so
  the profile does not depend on the order they come in.

- sac:

  Optional `sac_range` object from
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  to take the range, nugget and correlation function from. Pass one
  fitted with `detrend = "reml"`, say. It must describe the variable the
  profile scores: the raw response without `predictor_vars`, the
  residuals on them with; a sac whose `detrended` attribute says
  otherwise is used with a warning. Its range is read in its own CRS
  (`attr(sac, "crs")`), to which the points are transformed first, so
  the area, the floor and `cell_diam_median` are then in that CRS's
  units. A sac whose range was rejected (`NA` with a `rejected_reason`)
  gives `cp` its nugget, with a warning, and leaves `reliability` `NA`,
  since that needs the range; one whose model did not converge, or whose
  range is below the shortest lag fitted (a structure that cannot be
  told from a nugget, so the nugget is not identified either), gives
  neither. Under `select_on = "split"` the sac must come from the
  selection half alone: run the profile once without it, fit the sac on
  `data_sf[attr(p, "split")$selection, ]` and pass it to a second call
  with the same `seed`, which makes the same split. When `NULL` and a
  response is given, one is estimated on the subsample (its selection
  half under `"split"`) with the same `predictor_vars`. A plain number
  is taken as the range alone, in the units of the CRS the profile is
  computed in (metres for lon/lat input): it sets the floor, and `cp`
  and `reliability`, which need a fitted model, are `NA`. A `units`
  object is refused rather than read as a number in whatever unit it was
  written in.

- select_on:

  `"all"` (default) profiles every point; `"split"` reads the response
  on one spatially blocked half only (the OLS fit, the variogram, `rss`,
  `cp` and `moran_z`) and returns the other half as the set to estimate
  on, in the `"split"` attribute. The cells, `wss`, `elbow`, and the
  extent and point counts behind the bounds, `cp`'s variance term and
  `reliability` still come from every point, because the tessellation
  the count is for is built on every point: the levels are cell counts
  for the whole layer, and the estimation half's response never touches
  them. See the "Post-selection inference" section of
  [`determine_optimal_levels`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md);
  the profile reads the response whenever `response_var` is given.

- range_floor:

  `TRUE` (default) starts the ladder at the range floor when the data
  support it; `FALSE` starts it at 2 whatever the range, and the floor
  is only reported in the bounds and the print. Use `FALSE` to compare
  profiles whose range estimates differ (see "The ladder and its
  bounds").

## Value

A data.frame of class `resolution_profile` with one row per level and
columns `levels`, `wss`, `wss_spread` (relative spread of WSS across the
restarts), `elbow`, `cell_n_min`, `cell_n_median`, `cell_diam_median`
(twice the median RMS radius of the cells, in coordinate units: about
0.8 of the side of a square cell of the same area, and less for finer
cells, so a size to compare levels by rather than a width), `rss`, `cp`,
`moran_i`, `moran_z` and `reliability`; columns a missing input leaves
undefined are `NA`. Attributes: `bounds` (a list with `floor`,
`ceiling`, `ceiling_from` (`"min_cell_n"`, `"distinct locations"` or
`"sample_n"`, whichever bound it), `supported`, `area`, `range`, `n`
(the points in the layer), `n_sample` (the points the k-means fits ran
on), `n_distinct`, `min_cell_n`, `range_floor`), `variogram` (a list
with `nugget`, `psill`, `range` (`NA` when rejected), `model` and
`detrended` (whether the sac says it is a variogram of residuals; `NA`
when it does not say); `NULL` when none was usable), `variable`
(`"response"`, `"residuals"` or `NA`), `wss_bumps`, `nstart`, `sac` (the
range object used) and, with `select_on = "split"`, `split` (a
`spatialkit_split`: `selection` and `estimation`, integer row positions
in `data_sf`, with the `method` and `seed` that made them).

## The ladder and its bounds

Cells are k-means clusters of the (projected) coordinates, fitted at
each level as the best of `nstart` k-means++ restarts (see
[`determine_optimal_levels`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
for why). Levels are spaced logarithmically, because cell diameter
scales as \\L^{-1/2}\\: a unit step wastes fits at large \\L\\ and
starves resolution at small. The ladder runs from a floor to a ceiling
the data impose. The ceiling is `floor(n / min_cell_n)`, with \\n\\
every point of the layer (not the subsample): cells with fewer than
`min_cell_n` points on average have too little support, and Moran's z is
not computable at nine cells or fewer in any case. It is also held to
the number of distinct locations, one cell on each, and to one short of
the number of points, which
[`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) needs. The
floor is `ceiling(area / range^2)` when an autocorrelation range is
available: cells wider than the range average over more than one patch
of the field. When the floor exceeds the ceiling the data cannot support
a tessellation that respects their own correlation structure; that is
reported as a finding (a logged warning, and
`attr(x, "bounds")$supported` is `FALSE`) and the ladder runs from 2 to
the ceiling anyway, so the profile still shows what each level costs. On
a layer larger than `sample_n` the ceiling is also held to half the
subsample (two subsample points per cell, the least a fitted cell can be
scored on); `ceiling_from` is then `"sample_n"`, a floor above that is
logged with a request to raise `sample_n`, and it does not make
`supported` `FALSE`.

The floor moves with the range estimate, and the ceiling with \\n\\, so
two profiles of similar data (the folds of a cross-validation, say) can
sit on either side of the point where the floor applies: one ladder then
starts at the floor and the other at 2, and the criteria that sit near
the bottom of the ladder (`reliability`, which routinely peaks at the
floor, and `elbow`) can differ between them by a factor of 10 or more.
To compare profiles, pass the same `levels` to each, or
`range_floor = FALSE` to start every ladder at 2 while the floor is
still reported.

## The criteria, and how each behaved when measured

- `elbow`:

  How far the WSS curve sags below a power law: on log-log axes,
  \\\log\\ WSS below the straight line from \\k = 1\\ (the total sum of
  squares) to the last level, in natural-log units (larger is better).
  Points with no cluster structure have a WSS close to \\c/k\\, which is
  straight on those axes, so the column is `NA` at every level unless
  the largest sag reaches \\\log 1.25\\, and the print says there is no
  elbow (the rule and its calibration are in
  [`determine_optimal_levels`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)).
  The classical chord on linear axes found a "knee" on such a layer
  anyway, at about \\\sqrt{L\_{first} L\_{last}}\\, where the ladder's
  ends put it. A level whose WSS is 0 (to within \\10^{-12}\\ of the
  total), one cell on every distinct location, which the ladder reaches
  when locations repeat, is left out of the line; when the other levels
  have no elbow, that fall to zero is the elbow, and its sag is measured
  with the WSS floored at \\10^{-12}\\ of the total, which puts it far
  above any other level's. Geometry only; it knows nothing of the
  response.

- `cp`:

  Mallows' \\C_p\\ of the piecewise-constant approximation of the
  response (or of its OLS residuals on `predictor_vars`) by cell means:
  \\RSS(L)/n + 2 \tau^2 L / n\\, with \\\tau^2\\ the nugget of the
  fitted variogram (lower is better). It estimates the error of
  predicting a new observation by the mean of its cell. When the cells
  are fitted to a subsample of \\m\\ of the \\N\\ points with a
  response, the penalty is split between the two: \\RSS(L)/m + \tau^2
  L_m / m + \tau^2 L / N\\, where the first two terms estimate the
  approximation error from the subsample (adding back the optimism of
  its own cell means, over the \\L_m\\ cells its scored points fall in)
  and the last is the variance of cell means built from all \\N\\, which
  is what the tessellation will carry. With no subsample it is the
  formula above. **Measured on simulated exponential fields (600 points
  on a 1000-unit extent, sill 1, 20 replicates): with a nugget of 0.3
  its minimum sat at the support ceiling in every replicate at effective
  ranges of 90, 300 and 900; with a nugget of 2 it was interior (median
  44, range 8–66).** On a smooth field with little noise the
  approximation keeps improving as cells shrink and the penalty is too
  small to stop it, so \\C_p\\ says "as fine as the support allows" and
  `min_cell_n` is what is choosing;
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  says so when that happens. It becomes a genuine interior criterion
  only when the nugget is a large share of the sill. With a nugget of 0
  (under \\10^{-4}\\ of the sill; usually a fit clipped at its lower
  bound) the penalty is 0 and \\C_p\\ descends to the ceiling whatever
  the field; the profile warns.

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
  1996), for square cells of the level's average area holding the
  level's average share of the layer's points with a response (larger is
  better). This is the shrinkage factor of Fay and Herriot (1979). It
  has an interior optimum, and a broad one: validated against the
  empirical reliability of true block means on simulated fields, the
  analytic and empirical optima agreed to within a level or two where
  the empirical estimate was stable, and the band within 2 percent of
  the maximum spanned a factor of 3–6 in \\L\\. Read the flat region,
  not the argmax. The domain term is taken over the convex hull the area
  is measured on, so rotating the layer does not move it. `NA` without a
  usable variogram, and when the variogram's range was rejected (see
  `sac`).

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
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md),
[`summary.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.resolution_profile.md)

## Examples

``` r
if (requireNamespace("gstat", quietly = TRUE)) {
  library(sf)
  # An exponential field with range parameter 200 (true effective range
  # 600 m) on a 1 km square, with a nugget of 0.6 on a unit sill: enough
  # noise for Mallows' Cp to have an interior optimum rather than descend
  # to the ceiling.
  set.seed(4)
  n <- 400
  xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
  D  <- as.matrix(dist(xy))
  xy$z <- as.numeric(t(chol(exp(-D / 200) + diag(0.6, n))) %*% rnorm(n))
  pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
  prof <- resolution_profile(pts, response_var = "z", n_levels = 12)
  print(prof)               # one row per level; print() because only the
                            # last value of a braced block is shown
  select_resolution(prof, criterion = "cp")
}
#> Resolution profile: 12 levels on 400 points
#>   ladder      : 4 to 44 cells (floor 4 from range 517, ceiling 44 from min_cell_n = 9)
#>   variogram   : nugget 0.582, partial sill 0.819, range 517
#>   scored on   : response; 25 k-means++ restarts per level; WSS rises at 0 step(s)
#>   elbow       : none; the WSS curve falls as it does with no cluster structure
#> 
#>  levels      wss wss_spread elbow cell_n_min cell_n_median cell_diam_median rss
#>       4 16700000   0.000836    NA         92          97.5              409 480
#>       5 13900000   0.032200    NA         56          72.0              365 445
#>       6 11500000   0.055800    NA         58          65.5              348 428
#>       8  8090000   0.177000    NA         33          50.5              284 433
#>      10  6240000   0.130000    NA         34          37.5              246 386
#>      12  5130000   0.121000    NA         19          33.5              223 373
#>      15  4010000   0.119000    NA         20          27.0              196 364
#>      18  3300000   0.106000    NA         14          23.5              180 343
#>      23  2510000   0.099100    NA         10          18.0              156 358
#>      28  2010000   0.127000    NA          8          14.0              136 329
#>      35  1490000   0.172000    NA          7          11.0              120 328
#>      44  1140000   0.145000    NA          5           9.0              104 327
#>     cp moran_i moran_z reliability
#>  1.210      NA      NA       0.863
#>  1.130      NA      NA       0.860
#>  1.090      NA      NA       0.855
#>  1.110      NA      NA       0.844
#>  0.993 -0.1140 -0.0764       0.833
#>  0.968 -0.0241  1.0700       0.821
#>  0.954 -0.0424  0.4020       0.805
#>  0.909 -0.0256  0.4340       0.789
#>  0.963  0.0635  1.4600       0.765
#>  0.904  0.0643  1.3800       0.742
#>  0.921  0.1930  3.2500       0.714
#>  0.945  0.2630  4.4500       0.681
#> Resolution by cp: 28 cells
#>   flat region : 18, 28 to 35 (3 of 12 levels)
```
