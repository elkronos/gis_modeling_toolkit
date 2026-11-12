# Score candidate cell counts on several criteria

Scores candidate numbers of cells (levels) for a layer of points, so
that the count chosen, by
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
or by eye, can be defended with the whole profile in place of one
criterion's optimum.

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
  nstart = NULL,
  seed = 123L,
  sac = "auto",
  select_on = c("all", "split"),
  range_floor = FALSE,
  partition = c("density", "kmeans"),
  init = c("kmeans++", "random"),
  density_power = 1,
  density_k = 8L,
  iter_max = NULL,
  estimate = "auto",
  obs_per_cell = "auto"
)

# S3 method for class 'resolution_profile'
print(x, digits = 3L, ...)
```

## Arguments

- data_sf:

  An sf object of points (other geometries are reduced to representative
  points with
  [`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md)).
  Features without a location, as given, once projected or in the CRS of
  `sac`, are dropped with a warning, and a layer that mixes XY and XYZ
  geometries or carries M is read in x and y (see
  [spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)).

- response_var:

  Name of the response column (numeric or logical), or `NULL`. Enables
  `cp` and `moran_z`, and, when it is numeric, the pooled estimate the
  default profile scores by `loo`. A variogram estimated from it also
  gives the range floor the bounds report and `reliability`; it is
  estimated on a ladder of counts or when `sac` asks for it, and not on
  the default ladder of one size, so `cp` and `reliability` are `NA`
  there. Rows where it, or a predictor, is missing or non-finite stay in
  the geometry and are left out of the OLS fit, the RSS, `cp`, `moran_z`
  and the count the support ceiling divides by `min_cell_n`; a logged
  warning gives their number, and an R warning says so when they lower
  the ceiling.

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
  many centres as points). When none is left the call stops with an
  error, before the variogram is estimated or read, so it raises no
  warning about the variogram.

- n_levels:

  Number of levels on the ladder. Default 20. The ladder has fewer when
  rounding puts two of them on one count, and one more when the range
  floor is added as a level (see "The ladder and its bounds"). Supplying
  it, or `min_cell_n`, `sample_n`, `levels` or `range_floor`, gives that
  ladder in place of the ladder of sizes `obs_per_cell = "auto"` would
  set.

- min_cell_n:

  Minimum average number of rows per cell that a level must keep (scored
  rows when `response_var` is given: those with a response and every
  predictor; points otherwise); sets the ceiling. Default 9. A layer
  that cannot give two cells that many is profiled at two cells, or at
  the `levels` asked for, with `supported = FALSE` and a warning. It is
  floored to a whole number, and it is not read on a ladder
  `obs_per_cell` sets (see "The ladder and its bounds").

- sample_n:

  Points are subsampled to this many before the k-means fits, as in
  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md),
  on the ladder of "The ladder and its bounds"; on a ladder of sizes,
  the default, the centres are fitted to every row. Default 1500. The
  columns read off the fitted cells (`wss`, the `cell_` and `cells_`
  columns, `rss`, `moran_i`, `moran_z`) describe the subsample; the
  bounds, `supported`, the variance term of `cp` and `reliability`
  describe the whole layer, so the answer does not change with
  `sample_n` except through the fits. The fits are those of the
  subsample, which `seed` draws, so on a layer of more than `sample_n`
  points the count read off the profile can move with `seed` and with
  `sample_n`. The level
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  chooses off such a profile without an estimate is refitted to every
  point before its seeds are handed on (see `refit_rows` under "Value");
  a profile with an estimate scores it on every row in the cells of the
  centres fitted, and hands those centres on.

- nstart:

  Starts per level; the best is kept. They are k-means++ starts unless
  `init` says otherwise. `NULL` (default) is the kind's own number: 3
  for the kind `"density"` and 25 for the kind `"kmeans"`. A number is
  run under either `partition`. `attr(, "nstart")` holds the number that
  was run.

- seed:

  RNG seed for the subsample and the restarts; restored afterwards.
  Default 123. The rows are put in coordinate order before either, so
  the profile does not depend on the order they come in. Each level's
  restarts run on a stream of their own, seeded from `seed` and the
  level, so for one `seed` a level gets the same cells on any ladder
  that holds it and is fitted to the same rows, and with them the same
  `wss`, `rss`, `moran_z` and the other columns read off its cells. On a
  layer of more than `sample_n` rows a ladder of counts is fitted to a
  subsample and a ladder of sizes to every row, so there the two can
  give a level other cells. `cp` also reads its noise variance, which
  without a usable variogram comes from the ladder's finest level
  (`attr(, "cp_noise")`).

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
  gives `cp` its nugget (held to the semivariance of nearest-neighbour
  pairs, and not at all when it is at the fit's lower bound), with a
  warning, and leaves `reliability` `NA`, since that needs the range;
  one whose model did not converge (its fit has no minimum inside the
  ranges searched, so its nugget and sill are not fitted values), or
  whose range is below the shortest lag fitted (a structure that cannot
  be told from a nugget, so the nugget is not identified either), gives
  neither, and `cp` takes its noise variance from the finest level (see
  `cp` above). Whether the model converged is read from the model
  (`attr(attr(sac, "variogram_model"), "converged")`, which
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  sets and describes), not from the reason: the reason names one ground
  for the refusal, and such a fit is refused as past the largest lag.
  The warning then names both. A sac that carries no such flag (one made
  by hand, or built from REML parameters) is judged by its reason. Under
  `select_on = "split"` the sac must come from the selection half alone:
  run the profile once without it, fit the sac on
  `data_sf[attr(p, "split")$selection, ]` and pass it to a second call
  with the same `seed`, which makes the same split. When `NULL` and a
  response is given, one is estimated on the subsample (its selection
  half under `"split"`) with the same `predictor_vars`. A plain number
  is taken as the range alone, in the units of the CRS the profile is
  computed in (metres for lon/lat input): it gives the range floor,
  `reliability`, which needs a fitted model, is `NA`, and `cp` takes its
  noise variance from the finest level whose cells hold at least two
  scored rows on average (see `cp` above). A sac refused because the
  response shows no spatial structure (either of the two
  `rejected_reason` values of
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  that say so; the model it attaches is a pure nugget) gives `cp` no
  nugget either, and
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  then returns `structured = FALSE`. A `units` object is refused rather
  than read as a number in whatever unit it was written in. The string
  `"auto"`, the default, is `NULL` on every ladder but one that
  `obs_per_cell` sets. On that ladder no variogram is estimated or read:
  `cp` and `reliability` are `NA` at every level, `cp` takes no noise
  variance from the finest level in its place, `attr(, "cp_noise")` is a
  list whose `source` is `"none"` and whose `reason` says why, and the
  `variogram` and `sac` attributes are `NULL`. The other columns do not
  read the variogram and are the numbers `sac = NULL` gives.
  `sac = NULL` estimates the variogram on a ladder of sizes too, and
  gives `cp` and `reliability` there.

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

  `FALSE` (default) starts the ladder at 2 whatever the range, and the
  range floor is reported in the bounds and the print without being
  applied. `TRUE` starts the ladder at the floor when the data support
  it and the floor is below the ceiling; an estimated range then sets
  the fewest cells any criterion can choose (see "The ladder and its
  bounds").

- partition:

  What the centres of the cells minimise, at every level of the ladder.
  `"density"` (default): the sum of squared distances from the points to
  their nearest centre, each point weighted by its local density to the
  power `density_power`; compact cells that, at the default power, hold
  about the same number of points each once there are more cells than
  dense places, and counts about as uneven as `"kmeans"` gives below
  that. `"kmeans"`: the same sum unweighted; compact cells whose counts
  are several times higher in dense places than in sparse ones, from 25
  starts under `nstart = NULL`. See "The partition".

- init:

  How each start is drawn. `"kmeans++"` (default): the centres one at a
  time, each in proportion to its squared distance from those already
  drawn, times its weight under `partition = "density"`. `"random"`: `k`
  distinct locations at random, as
  [`kmeans()`](https://rdrr.io/r/stats/kmeans.html) draws them.

- density_power:

  For `partition = "density"`: the power the local density is raised to
  for a point's weight, a number from 0 to 2. Default 1, which evens out
  the cell counts once the cells outnumber the dense places (see "The
  partition"). `0` weights every point alike, the objective of
  `"kmeans"` (reached by Lloyd's iteration, so the centres are close to
  that kind's and not equal to them); values between give counts between
  the two; values above 1 give the cells of dense places fewer points
  than those of sparse places. Setting it under `partition = "kmeans"`,
  which does not read it, is a warning.

- density_k:

  For `partition = "density"`: the number of nearest neighbours the
  local density is read from. Default 8. It is held to one fewer than
  the distinct locations. A small number follows the density closely and
  noisily; a large one smooths it. Setting it under
  `partition = "kmeans"` is a warning.

- iter_max:

  The most iterations a start runs for: a single positive whole number,
  or `NULL` (default) for the kind's own limit, 50 for
  `partition = "kmeans"` from k-means++ starts and 100 for random starts
  and for `partition = "density"`. A kept start that was cut off at the
  limit is a warning of class `"spatialkit_not_converged"`, except for
  the kind `"kmeans"` from k-means++ starts under `iter_max = NULL`,
  which reports nothing.

- estimate:

  The cell value the profile scores, which is the one
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  is to report: `"auto"` (default; see below), `NULL` for none, `"mean"`
  for the plain mean of a cell's rows, `"pooled"` for the mean pulled
  toward the mean of the 8 nearest cells. With an estimate the profile
  gains five columns and `attr(, "estimate")` (see "Value"), and every
  column and attribute it has without one is unchanged. `loo` is the
  leave-one-out mean squared error of that estimate: each row predicted
  by the value of its cell computed without it. Under `"mean"` a row in
  a cell of two or more is predicted by the mean of the others, and a
  row alone in its cell by the mean of the cell of its nearest other
  row. Under `"pooled"` the cell's mean and its own weight are computed
  again without the row, and the spread inside the cells, the spread
  between a cell and its target, and the targets are held. It is
  computed on every row of the layer with a finite response (of the
  selection half under `select_on = "split"`), not on the subsample the
  centres were fitted to, and each row is scored in the cell of its
  nearest centre, the cell
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
  hands on. It needs no variogram. Rows at one location are left out one
  at a time. `"auto"` is `"pooled"` when there is a numeric
  `response_var` with a finite value and no `predictor_vars`, and `NULL`
  otherwise. `NULL` scores none and gives the profile without the five
  columns. `"mean"` and `"pooled"` are errors without a `response_var`,
  with `predictor_vars` (the profile would then score residuals, and the
  summary reports cell values of the response), with a logical response,
  and with a response that has no finite value. Under
  `select_on = "split"` a selection half with no finite response is not
  an error: the five columns are `NA` at every level, with a warning of
  class `spatialkit_no_scored_rows`.

- obs_per_cell:

  Rows per cell: `"auto"` (default; see below); `NULL` for the ladder
  described under "The ladder and its bounds", whatever the other
  arguments are; or one number or several, each at least 1 (fractions
  are allowed), for a ladder of sizes. The level of a size \\s\\ is
  `max(2, round(n / s))` cells, with \\n\\ the rows that have a response
  (and, with `predictor_vars`, every predictor), or every point when
  there is no `response_var` or no row has a response. Under
  `select_on = "split"` \\n\\ counts the rows of both halves, since the
  cells are built on every point. On a ladder of sizes the centres are
  fitted to every row, no subsample is drawn, and `n_levels`,
  `min_cell_n`, `sample_n` and `range_floor` are not read: supplying
  `n_levels`, `min_cell_n`, `sample_n` or `range_floor = TRUE` beside a
  numeric `obs_per_cell` is a warning of class
  `spatialkit_ladder_unused`, and `levels` beside it is an error. A
  level is held to the most centres k-means can place (the distinct
  locations, and one short of the points). `"auto"` is `NULL` when
  `levels` is given or the caller supplies any of those four arguments;
  otherwise 4.5 when the estimate is `"pooled"`, the seven sizes 18, 12,
  9, 6, 4.5, 3 and 2 when it is `"mean"`, and `NULL` without an
  estimate. The cost of a ladder of sizes is bounded before any fit (see
  "Conditions").

- x:

  A `resolution_profile`, for
  [`print()`](https://rdrr.io/r/base/print.html).

- digits:

  For [`print()`](https://rdrr.io/r/base/print.html): significant digits
  for the table. Default 3.

- ...:

  Ignored.

## Value

A data.frame of class `resolution_profile`, one row per level. Columns a
missing input leaves undefined are `NA`.

- `levels`:

  The number of cells.

- `wss`:

  The within-cell sum of squares of the coordinates, every point counted
  once, under either `partition`.

- `wss_spread`:

  The relative spread of WSS across the restarts; under
  `partition = "density"`, of the weighted sum the starts minimise.

- `elbow`:

  The sag of the WSS curve below a power law (see "The criteria").

- `cell_n_min`, `cell_n_median`:

  The fewest and the median points per cell.

- `cell_diam_median`:

  Twice the median RMS radius of the cells, in coordinate units: about
  0.8 of the side of a square cell of the same area, and less for finer
  cells, so a size to compare levels by rather than a width.

- `rss`, `cp`:

  The residual sum of squares of the scored rows about their cell means,
  and Mallows' \\C_p\\.

- `moran_i`, `moran_z`:

  Moran's I of the cell means and its standardised deviate.

- `reliability`:

  The between-cell signal's share of the spread in the cell means.

- `cells_resp`:

  The cells that hold at least one scored row.

- `cell_resp_mean`:

  The scored rows per such cell; compare it with `min_cell_n`.

Like the other columns read off the fitted cells, `cells_resp` and
`cell_resp_mean` count rows of the subsample when the layer is larger
than `sample_n`: the cells of the layer then hold about
`n_response / n_scored` times as many scored rows as `cell_resp_mean`
shows (under `select_on = "split"`, `split$n_response / n_scored`: both
count the rows of the selection half).

With an `estimate` the frame has five more columns, after
`cell_resp_mean`:

- `scored_per_cell`:

  The rows `loo` read over the cells that hold one; the layer's
  counterpart of `cell_resp_mean`, which counts the subsample's rows.

- `loo`:

  The leave-one-out mean squared error of the estimate; lower is better.

- `loo_se`:

  The standard error over rows of that level's mean squared error. It is
  not the standard error of a difference between two levels, whose
  errors come from the same rows.

- `share_alone`:

  The share of the scored rows that are the only scored row of their
  cell.

- `own_weight`:

  The median over cells of the weight a pooled value gives the cell's
  own mean; `NA` under `"mean"`.

`loo`, `loo_se` and `own_weight` are `NA` at a level where no cell holds
two scored rows, at every level when the response is constant over the
rows `loo` reads, and at every level when those rows' values lie beyond
what the arithmetic of the estimate holds
(`spatialkit_pool_unavailable`). All five are `NA` at every level when
no row of the selection half has a response (`select_on = "split"`;
`spatialkit_no_scored_rows`).

Attributes:

- `bounds`:

  A list: `floor`, `ceiling`, `ceiling_from` (`"min_cell_n"`,
  `"min_cell_n (rows with a response)"`, `"distinct locations"`,
  `"sample_n"`, `"fewer than two supported cells"` or
  `"fewer than two distinct locations"`, whichever bound it),
  `supported`, `area`, `range`, `n` (the points in the layer),
  `n_sample` (the points the k-means fits ran on), `n_distinct`,
  `min_cell_n`, `range_floor`, `n_response` (the rows of the layer that
  carry a response and, with `predictor_vars`, every predictor, which
  the support ceiling counts; `NA` without a `response_var`), `n_scored`
  (the rows of the subsample the response criteria read; `NA` without a
  `response_var`), `levels_given` (whether `levels` set the ladder),
  `predictors` (whether `predictor_vars` were given) and, where rows
  have a response and none of them every predictor, `n_response_only`
  (the rows with a finite response). On a ladder that `obs_per_cell`
  sets, `bounds` describes that ladder: `ceiling` is the finest level;
  `ceiling_from` is `"obs_per_cell"`, `"distinct locations"`,
  `"fewer than two cells of the size"` or
  `"fewer than two distinct locations"`; `supported` is `FALSE` when the
  finest size gives fewer than two cells or the range floor is above the
  finest level; `levels_given` is `TRUE`; `n_sample` is every point;
  `min_cell_n` holds the argument, which is not read; and two fields are
  added, `obs_per_cell` (the sizes, largest first) and `size_levels`
  (the level of each, named by size; `NA` for a size of `"auto"` that
  gives fewer than two cells and was left out).

- `variogram`:

  A list with `nugget`, `psill`, `range` (`NA` when rejected), `model`
  and `detrended` (whether the sac says it is a variogram of residuals;
  `NA` when it does not say); `NULL` when none was usable.

- `variable`:

  `"response"`, `"residuals"` or `NA`.

- `wss_bumps`, `nstart`:

  The number of rising steps on the WSS curve, and the starts run per
  level.

- `sac`:

  The range object used.

- `split`:

  With `select_on = "split"`: a `spatialkit_split` with `selection` and
  `estimation` (integer row positions in `data_sf`), the `method` and
  `seed` that made them, and with a response `n_response`, the rows of
  the selection half with a finite response and every predictor, counted
  on the layer before any subsample.

- `cp_noise`:

  A list with `source` (`"variogram nugget"`,
  `"finest-level residual mean square"` or `"short-lag semivariance"`),
  `value` (the \\\tau^2\\ used) and `level` (the level it was read at,
  or `NA`). Under the short-lag semivariance it also has `cap`, that
  semivariance, and `uncapped`, a list of the `source`, `value` and
  `level` it held. When the finest level whose cells hold two scored
  rows on average has a residual mean square of 0 to within rounding,
  `source` and `value` are `NA`, and the list has `level` (that level),
  `rows_per_cell` (its scored rows a cell), `rms` (its residual mean
  square: 0, or the rounding it was read as) and `reason`, which
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  gives. On a ladder of sizes under `sac = "auto"` it is a list whose
  `source` is `"none"` and whose `reason` says why. `NULL` when `cp` has
  no noise variance for another reason; see `cp` under "The criteria".

- `one_cell`:

  What one cell for the whole layer scores, the reference a level has to
  beat: a list with `n` (the scored rows), `rss` (their total sum of
  squares), `cp` (the same formula at one cell, or `NA` when `cp` has no
  noise variance), `constant` (whether the scored variable is constant
  to within rounding), `cp_levels` (the `cp` column it belongs to) and
  `f_p` (named by level, the p-value of the F test of that level's cell
  means against one mean, \\(TSS - RSS) / (L - 1)\\ over \\RSS / (m -
  L)\\ on the \\L\\ cells that hold a scored row and the \\m\\ scored
  rows). `NULL` without a response.

- `centres`, `centres_crs`, `hull`:

  A list named by level of the two-column matrices of the scored
  partitions' centres, the CRS they are in, and the convex hull of the
  points profiled (an `sfc` in that CRS). A k-means partition is the
  Voronoi partition of its centres, so these seed the cells each level
  was scored on:
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  returns the chosen level's as `$seeds`, and
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
  and
  [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)
  use them when given the selection or the profile. The centres carry
  the hull as their `"cover"` attribute, which those two functions hand
  on so that Voronoi cells built on the seeds reach the edge of the
  points.

- `refit_rows`:

  When the centres were fitted to a subsample and the profile has no
  estimate: a list with `xy`, the coordinates of every point in the CRS
  of the centres, and `iter_max`, the iteration limit of the fits, from
  which
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  refits the level it chooses to every point. The seeds are then the
  chosen level's centres refitted to every point, whose cells are not
  the cells scored (see
  [`select_resolution`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md),
  "Value").

- `partition`:

  On a profile built with `partition = "density"` or `init = "random"`:
  a list with `kind`, `init`, `density_power` and `density_k` (both `NA`
  for `kind = "kmeans"`); for `kind = "density"`, `objective`, the
  weighted sum the centres minimise, named by level; `iter_max`, the
  iteration limit; and `iter` and `converged`, the iterations each
  level's kept start took and whether it stopped before the limit, named
  by level. The selection carries the first four as `$partition`, its
  seeds and the seeds the seeding functions return as
  `attr(, "partition")`, and a Voronoi tessellation built on them as
  `params$partition`. A profile of the kind `"kmeans"` from k-means++
  starts has no such attribute, and neither has anything made from it.

- `estimate`:

  With an `estimate`: a list with `value` (`"mean"` or `"pooled"`),
  `neighbours` (8, or `NA` under `"mean"`), `response` (the column), `n`
  (the rows `loo` read), `rows` (`"all"` or `"selection half"`),
  `features` (`"points"`, or `"other"` when the layer held features that
  were reduced to points), `seed`, `constant`, `loo_one_cell` and
  `loo_one_cell_se` (the same error with one mean for all those rows,
  the reference a level has to beat; `NA` for a constant response, for
  fewer than two rows and for values beyond the arithmetic of the
  estimate), `loo_levels` (the `loo` column those references belong to,
  named by level), and, named by level, `sigma2` (the spread of the rows
  about their cell's mean), `tau2` (the mean squared distance between a
  cell's level and its target; `NA` under `"mean"`), `cells` (the cells
  that hold a scored row) and `f_p` (the p-value of the F test of those
  cells' means against one mean). `constant` says whether the response
  is constant over those rows by the estimate's test, which reads the
  raw values and the values centred on their median and scaled, so that
  1e15 plus a spread of 1 is not constant while values that differ only
  by their rounding, 0.1 beside 0.3 / 3, are; it is `FALSE` for fewer
  than two rows. `one_cell$constant`, the test for `cp`, reads the
  scored rows of the subsample the same way, and residuals against the
  rounding of the regression and of the response.

## Details

Called with a numeric response, no predictors and none of `levels`,
`n_levels`, `min_cell_n`, `sample_n` or `range_floor`, the profile is
one level: the count that gives 4.5 rows per cell
(`obs_per_cell = "auto"`), scored for pooled values
(`estimate = "auto"`) by the leave-one-out error `loo`, with no
variogram fitted (`sac = "auto"`), in cells of the density kind
(`partition = "density"`) kept from 3 k-means++ starts
(`nstart = NULL`). `cp` and `reliability` need the variogram and `elbow`
needs a curve, so the three are `NA` there; `moran_z` is computed.
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
reads the count off it by the size.

Any of those five ladder arguments, `obs_per_cell = NULL`, or a call
with no estimate (no response, `predictor_vars`, or a logical or wholly
missing response) gives a ladder of counts (see "The ladder and its
bounds"). At every level the profile then scores the within-cluster sum
of squares of the coordinates and its elbow, Mallows' \\C_p\\ of the
piecewise-constant approximation of the response by cell means, Moran's
\\z\\ of the cell means and their reliability, and, with an estimate,
the leave-one-out error of the cell value. Mallows (1973) presents
\\C_p\\ itself as a display of the bias-variance trade-off across
candidates, and the profile is that display, with the other criteria
alongside.

[`print()`](https://rdrr.io/r/base/print.html) gives a header that
states the ladder and its bounds, the variogram (or why none was
fitted), what was scored and on how many starts, what one cell for the
whole layer scores, the estimate, whether the WSS curve has an elbow
and, under `select_on = "split"`, the rows of each half; then the table,
one row per level, at `digits` significant digits.

## The partition

The cells at every level are the Voronoi cells of a set of centres: a
point belongs to the centre nearest to it. `partition` says what the
centres minimise, and every level of the ladder is built with the same
kind.

- `"kmeans"`:

  The sum, over the points, of the squared distance to the nearest
  centre. The centres follow the square root of the point density. The
  cells are compact, smaller where the points are dense and larger where
  they are sparse, and a cell in a dense place holds several times the
  points of a cell in a sparse place.

- `"density"` (default):

  The same sum with each point weighted by its local density to the
  power `density_power`. The local density at a point is the number of
  points in the disc that reaches its `density_k`-th nearest neighbour,
  over the area of that disc. With `density_power = 1` the centres
  follow the point density itself. The cells are compact: many small
  ones where the points are dense, a few large ones where they are
  sparse. They hold about the same number of points each once there are
  more cells than dense places (clusters of points). With fewer cells
  than dense places a cell takes in whole clusters, large and small, and
  the counts are about as uneven as under `"kmeans"`.
  `density_power = 0` weights every point alike, which is the objective
  of `"kmeans"`; a power above 1 goes past equal counts and gives the
  cells of dense places fewer points than those of sparse places.

A ladder starts at 2 cells, so its low levels are usually below the
number of dense places, and the columns `cell_n_min` and `cell_n_median`
show where the two kinds part.

The density kind is the default. The default size, 4.5 rows per cell,
puts the count above the number of dense places on a layer whose
clusters hold more than a few observations each, so the default kind and
the default count go together. A coarse ladder the caller sets (with
`levels`, large sizes, or the low levels of the automatic ladder) may
not, and at those levels the density kind does not even the counts.
`partition = "kmeans"` gives plain k-means cells, from 25 k-means++
starts unless `nstart` says otherwise.

`init` says how each of the `nstart` starts is drawn before it is
refined; the best start is kept. `"kmeans++"` draws the centres one at a
time, each a point taken in proportion to its squared distance from the
centres already drawn (times its weight, under `"density"`), which
spreads the start over the points. `"random"` draws `k` distinct
locations at random, as
[`kmeans()`](https://rdrr.io/r/stats/kmeans.html) does.

A start of the kind `"kmeans"` is refined by
[`kmeans()`](https://rdrr.io/r/stats/kmeans.html) (Hartigan and Wong),
for at most 50 iterations from a k-means++ start and 100 from a random
one. A start of the kind `"density"` is refined by Lloyd's iteration,
for at most 100 steps: each point goes to its nearest centre, each
centre moves to the weighted mean of its points, and a centre left
without a point stays where it is. `iter_max` sets another limit. The
start with the lowest sum is kept (under the density kind a later start
replaces an earlier one only when its sum is lower by more than one part
in ten million, so of starts whose sums agree to rounding the first is
kept), and a kept start that was cut off at the limit is reported (see
"Conditions"): its centres are those of the last step. A large layer can
need more steps than the default limit.

The density kind depends on the points alone. It works on the distinct
locations in coordinate order (see below for first coordinates equal but
for their last bits), so the order of the rows does not matter. Its
tolerances do not depend on the unit of the coordinates: the neighbour
distance is not taken below a millionth of the extent of the points (the
diagonal of their bounding box), the iteration stops when no centre
moves by \\10^{-12}\\ of it, and the tie rules below are fractions of
the extent or of the quantities they compare. On gridded locations many
points are equally near two neighbours or two centres, and the last bits
of the coordinates, which move with the unit and the origin, would
settle those ties. So squared distances within one part in ten million
of each other count as equal, a tie between neighbours goes to the lower
first coordinate and then the lower second, a tie between centres goes
to the lower-numbered one (between centres a squared distance ties with
the least when it is no more than the square of the least distance
lengthened by a margin, plus one part in ten million of the least
squared distance; the margin is one part in 1e10 of the largest
coordinate of the centres, at most a millionth of the diagonal of their
bounding box and at least four times the spacing of doubles at that
coordinate), first coordinates within that part of the extent count as
equal when the rows and the locations are put in order, and so do
weights within that part of each other when the first k-means++ centre
is drawn. The cells are then the same in another unit or with another
origin as long as no coordinate moves by more than about a
two-hundred-millionth of the distance between neighbouring locations,
but for one case: the margin of the centres grows with their largest
coordinate (about 0.46 mm at a UTM northing near 4,550 km, about half a
micrometre once the origin lies next to data 5 km across), and the rows
are given to the centres in the layer's own coordinates, so a row whose
distances to two centres differ by less than the margin with one origin
and by more with another goes to the lower-numbered centre with the
first and to the nearer with the second. The rounding PROJ leaves when
it changes the unit or the false origin, about \\10^{-9}\\ m, is below
that for locations more than about a quarter of a metre apart.
Coordinates rounded more coarsely can give other cells: a lattice at 250
m with its coordinates stored in feet to four decimals does. Points at
one location count as that many points at it and share a cell; the
neighbour distance is measured to the `density_k`-th nearest *other*
location, so it is not zero where locations repeat. The neighbour
searches use FNN when it is installed. The result is the same without
it; the search then takes time that grows with the square of the number
of points.

The kind `"kmeans"` meets the same ties inside
[`kmeans()`](https://rdrr.io/r/stats/kmeans.html), which settles them by
the last bits of the coordinates and takes no tolerance. So it is handed
the coordinates centred on the middle of their bounding box, divided by
its diagonal and rounded to \\2^{-30}\\ of it (about a billionth); the
starts are drawn from those coordinates and ranked by their sums, and
the cells it returns are read back in the layer's own coordinates, so
the centres, `wss` and every column are in its unit. The rows are put in
order as for the density kind. The cells are then the same in another
unit or with another origin unless the change moves a coordinate across
the middle of a rounding step, which the rounding PROJ leaves (about
\\10^{-9}\\ m) does more often the smaller the extent. Distinct
locations that would round to one point keep their coordinates
unrounded, so they stay distinct.
[`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) also settles
such ties by the last bits of its own compiled arithmetic, which differ
between platforms, so on a layer where ties decide cells the kind
`"kmeans"` can return other cells on another platform from the same
coordinates and `seed`, and with them every column read off the cells.

The columns keep their meaning under `"density"`, because each is read
off the cells and not off the way they were drawn. `wss` is the
within-cell sum of squares of the coordinates with every point counted
once, taken about each cell's mean, so `wss` and `elbow` compare between
kinds; the weighted sum the centres minimise is kept in
`attr(, "partition")$objective`, and `wss_spread` is its spread over the
starts. `cp`, `moran_z` and the test against one cell need cells drawn
without the response, which both kinds give. `reliability` is computed
for square cells of the level's average area and count and never reads
the cells, so it is the same number under either kind. The `elbow` is
found more readily under `"density"` where the points cluster: once each
cluster has its centres, further centres go into the clusters, the
sparse ground between them keeps its sum of squares, and the curve
flattens.

A profile of the density kind, or of either kind from random starts,
carries `attr(, "partition")`, and
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md),
the seeding functions and
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
hand that record on with the cells (see "Value"). A result with no such
record was built by plain k-means: from k-means++ starts when it comes
from a profile, and from random starts when a seeding function fitted
it.

## The ladder and its bounds

This is the ladder of a call that sets `levels`, `n_levels`,
`min_cell_n`, `sample_n` or `range_floor`, or passes
`obs_per_cell = NULL`, and of a call with no estimate: one without a
response, with `predictor_vars`, or with a logical or wholly missing
response. Cells are clusters of the (projected) coordinates of the kind
`partition` names, fitted at each level as the best of `nstart`
k-means++ starts (see
[`determine_optimal_levels`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
for why) unless `init` says otherwise (see "The partition"). Levels are
spaced logarithmically, because cell diameter scales as \\L^{-1/2}\\: a
unit step wastes fits at large \\L\\ and starves resolution at small.
The ladder runs from 2 to a ceiling the data impose. The ceiling is
`floor(n / min_cell_n)`: cells with fewer than `min_cell_n` rows on
average have too little support. With a `response_var`, \\n\\ is the
number of rows of the layer (not of the subsample) that are scored:
those with a response and, with `predictor_vars`, every predictor, since
the cell means the criteria read are built from those rows alone.
`ceiling_from` is then `"min_cell_n (rows with a response)"` when some
rows are not scored (the label is the same when it is a predictor they
lack), and a warning of class `spatialkit_response_ceiling` gives both
counts and names the bound in force. Without a `response_var`, or when
no row has a response, \\n\\ is every point. The ceiling is also held to
the number of distinct locations, one cell on each, and to one short of
the number of points, which
[`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) needs. The
count \\n / L\\ is the rows per cell only when the scored rows are
spread over the whole extent. Where they sit in one part of it, the
cells elsewhere hold none and the occupied cells hold more than \\n /
L\\, so the ceiling is lower than those cells need: `cells_resp` and
`cell_resp_mean` give the cells that hold a scored row and the scored
rows per such cell at each level (of the subsample, when the cells are
fitted to one), and profiling the scored rows alone bounds the ladder by
the area they cover.

When fewer than two cells meet `min_cell_n`
(`floor(n / min_cell_n) < 2`), or the layer has one distinct location,
there is no ladder: the profile scores two cells so that level is still
visible (or the `levels` asked for, as given), `supported` is `FALSE`,
`ceiling_from` is `"fewer than two supported cells"` (or
`"fewer than two distinct locations"`) and a warning of class
`spatialkit_not_supported` names both numbers.

`min_cell_n` sets the ceiling of this ladder. It is a mean number of
rows per cell and not a minimum: a cell of a level at the ceiling can
hold fewer (`cell_n_min`). It is floored to a whole number, so 4.5 is
read as 4. With `levels` it does not choose the levels scored, and it
still sets `bounds$ceiling`, `supported` and the
`spatialkit_not_supported` warning.

The range floor is `ceiling(area / range^2)` when an autocorrelation
range is available: cells wider than the range average over more than
one patch of the field. It is reported in the bounds and the print, and
by default it is not applied: the ladder starts at 2. When the layer
supports the floor and it lies between the ladder's ends, the floor is
one of the ladder's levels all the same (so the ladder can have one more
level than `n_levels`): `select_resolution(criterion = "reliability")`
reads from the floor up and returns the floor itself, not the next level
the spacing happened to put above it. The range is an estimate, and a
response with no spatial structure can still be given one, so a ladder
started at the floor can leave out the coarse counts such a response
supports. `range_floor = TRUE` starts the ladder at the floor when the
layer supports it and it lies below the ceiling. When the floor exceeds
the ceiling the data cannot support a tessellation that respects their
own correlation structure; that is reported as a finding (a warning of
class `spatialkit_not_supported`, and `attr(x, "bounds")$supported` is
`FALSE`) and the ladder runs from 2 to the ceiling (or the `levels`
asked for are scored, as given), so the profile still shows what each
level costs. A floor equal to the ceiling would leave one level and
nothing to compare it with, so that ladder starts at 2 as well. On a
layer larger than `sample_n` the ceiling is also held to half the
subsample's scored rows (two per cell, the least a fitted cell can be
scored on; half the subsample without a response); `ceiling_from` is
then `"sample_n"`, a floor above that is logged with a request to raise
`sample_n`, and it does not make `supported` `FALSE`.

The floor moves with the range estimate, and the ceiling with \\n\\. To
compare profiles of similar data (the folds of a cross-validation, say),
pass the same `levels` to each. `reliability` is the one criterion that
still reads the floor:
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
holds it to the levels at or above the floor by default, so its pick can
differ widely between two such profiles when the floor is supported on
one and not on the other; pass `range_floor = FALSE` there to read the
whole ladder.

A ladder that `obs_per_cell` sets is a ladder of sizes (see
"Observations per cell"). It has one level per size \\s\\,
`max(2, round(n / s))` cells, held to the most centres k-means can
place: the distinct locations, and one short of the points. Its centres
are fitted to every row. No subsample is drawn, and `n_levels`,
`min_cell_n`, `sample_n` and `range_floor` are not read. The finest size
takes the part `min_cell_n` has on the ladder above: `bounds$ceiling` is
the finest level, and `supported` is `FALSE` when the finest size gives
fewer than two cells or the range floor is above the finest level. The
range floor is reported and not applied.

A ladder of sizes under `sac = "auto"` fits no variogram, so `cp` and
`reliability` are `NA` at every level there, and `attr(, "cp_noise")`
says why. `sac = NULL`, or a `sac` of the caller's own, gives them.
`sac = "auto"` and `obs_per_cell = "auto"` are the defaults, so the
default profile has neither column.

## Observations per cell

`obs_per_cell` sets the ladder by a size: the mean number of rows a cell
is to hold. The level of a size \\s\\ is `max(2, round(n / s))` cells,
with \\n\\ the rows the size divides (see the argument). A size says how
many observations a cell value is to stand on, and that is the caller's
choice: no criterion makes it. Under a size neither the cells, which are
fitted to the coordinates, nor their number, which is a count of rows,
depend on the values of the response. Under
`select_resolution(criterion = "loo")` the count does.

For pooled values the leave-one-out error tends to keep falling as the
cells shrink, so its minimum tends to be the finest size offered, and
the end of the ladder would then choose the count (see "How this was
checked"). `select_resolution(criterion = "size")` therefore sets the
count of pooled values from a size, 4.5 unless another is given, and
`criterion = "loo"` reads the lowest error for plain means.

Two columns say what a size leaves a cell value to stand on.
`share_alone` is the share of the scored rows that are the only scored
row of their cell: the plain mean of such a cell is one observation, and
its pooled value has the smallest own weight a cell can have.
`own_weight` is the median weight a pooled value gives the cell's own
mean: at 1 the pooled values are the plain means, and at 0 each is the
mean of the cell's neighbours.

`obs_per_cell = "auto"`, the default, is resolved against the other
ladder arguments, in this order. Whether the caller supplied an argument
is what counts, so a value a wrapper passes on counts as supplied.

1.  With `levels`, `"auto"` is `NULL`: the levels are scored as given. A
    numeric `obs_per_cell` beside `levels` is an error, since each sets
    the ladder.

2.  With `n_levels`, `min_cell_n`, `sample_n` or `range_floor` supplied,
    `"auto"` is `NULL`: those arguments describe the ladder of "The
    ladder and its bounds", and the caller gets it, with the estimate's
    columns when there is an estimate.

3.  Otherwise `"auto"` is 4.5 when the estimate is `"pooled"`, the seven
    sizes 18, 12, 9, 6, 4.5, 3 and 2 when it is `"mean"`, and `NULL`
    without an estimate.

4.  A numeric `obs_per_cell` beside `n_levels`, `min_cell_n`, `sample_n`
    or `range_floor = TRUE` sets the ladder, with a warning of class
    `spatialkit_ladder_unused` that names the arguments that are not
    read.

On a ladder of sizes the centres are fitted to every row and the count
of cells grows with the rows, so a step that compares every point with
every centre grows with the square of the rows. `partition = "kmeans"`
has such a step in every pass and in a k-means++ start.
`partition = "density"` has it in a k-means++ start and, without FNN, in
every pass; with FNN and `init = "random"` it has none. The profile
stops before any fit, with an error of class `spatialkit_work_limit`,
when the points times the cells of the finest level times `nstart`
(times `iter_max` for the density kind without FNN) is above
`getOption("spatialkit.fit_work_max", 1e10)`. The message names the ways
on: fewer starts, which keep the best of fewer; a coarser size;
`partition = "density"` with `init = "random"`, which leaves a cell
without a point more often (see `spatialkit_empty_cells` under
"Conditions") and, where FNN could not be loaded, needs it reinstalled
first; or a higher limit. The work and the limit are written in full, so
that a work just above the limit does not print as the limit. The limit
counts the finest level, and a ladder of several sizes fits every one of
its levels. At the default size the work grows with the square of the
rows, and a call just below the limit runs for a long time and prints
nothing while it runs.
[spatialkit-options](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-options.md)
describes both limits.

Without FNN the density kind also finds its weights by comparing every
pair of distinct locations. A profile on a ladder of sizes stops, with
the same class, when the square of their number is above
`getOption("spatialkit.exact_pairs_max", 2e9)`, which is about 44,700
locations. The result of every search is the same with and without the
package. `moran_i` and `moran_z` need FNN and Matrix at a level where
more than 5,000 cells hold a scored row. Without them a ladder of sizes
has `NA` in both columns there, with a warning of class
`spatialkit_moran_skipped`, and so does a ladder set any other way (see
"Conditions"). At a level of up to 5,000 such cells, when either package
is missing, the neighbour weights of the cell means are held as a dense
matrix of cells by cells.

## The criteria

[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
reads `cp` by default on a profile without an estimate; on one with an
estimate it reads the size (pooled values) or `loo` (plain means).
`reliability` needs a variogram range and is `NA` without one; `elbow`
is `NA` unless the points cluster; and `moran_z` is a diagnostic of the
cell means, which no function reads as a count rule by default. With an
`estimate` the profile has a fifth column, `loo`, which is the last
entry. Each entry below says what the column measures and what it
returns.

- `elbow`:

  How far the WSS curve sags below a power law: on log-log axes,
  \\\log\\ WSS below the straight line from \\k = 1\\ (the total sum of
  squares) to the last level, in natural-log units (larger is better).
  Points with no cluster structure have a WSS close to \\c/k\\, which is
  straight on those axes, so the column is `NA` at every level unless
  the largest sag reaches \\\log 1.25\\, and the print says there is no
  elbow (the rule is given in
  [`determine_optimal_levels`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)).
  The classical chord on linear axes found a "knee" on such a layer
  anyway, at about \\\sqrt{L\_{first} L\_{last}}\\, where the ladder's
  ends put it. A level whose WSS is 0 (to within \\10^{-12}\\ of the
  total), one cell on every distinct location, which the ladder reaches
  when locations repeat, is left out of the line; when the other levels
  have no elbow, that fall to zero is the elbow, and its sag is measured
  with the WSS floored at \\10^{-12}\\ of the total, which puts it far
  above any other level's. Geometry only; it knows nothing of the
  response. It is `NA` unless the points cluster, and where they do it
  names the number of clusters, which is usually the first level or two
  of the ladder.

- `cp`:

  Mallows' \\C_p\\ of the piecewise-constant approximation of the
  response (or of its OLS residuals on `predictor_vars`) by cell means:
  \\RSS(L)/n + 2 \tau^2 L / n\\, with \\\tau^2\\ the nugget of the
  fitted variogram (lower is better). It estimates the error of
  re-measuring the sampled sites: each scored row observed again at its
  own location, with fresh noise of variance \\\tau^2\\, and predicted
  by the mean of its cell. For an observation at a new location the
  penalty needs the spread of the field inside a cell as well, which the
  nugget leaves out. This \\\tau^2\\ is the nugget, the noise of one
  measurement; it is not the `tau2` of a pooled value, which is the mean
  squared distance between a cell's level and its target. When the cells
  are fitted to a subsample of \\m\\ of the \\N\\ points with a
  response, the penalty is split between the two: \\RSS(L)/m + \tau^2
  L_m / m + \tau^2 L / N\\, where the first two terms estimate the
  approximation error from the subsample (adding back the optimism of
  its own cell means, over the \\L_m\\ cells its scored points fall in)
  and the last is the variance of cell means built from all \\N\\, which
  is what the tessellation will carry. With no subsample it is the
  formula above. On a smooth field with little noise the approximation
  keeps improving as cells shrink and the penalty is too small to stop
  it, so the minimum of \\C_p\\ is the support ceiling and `min_cell_n`
  is what is choosing;
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  says so when that happens. \\C_p\\ has an interior minimum when the
  nugget is a large share of the sill. A nugget of 0, or of at most
  \\10^{-4}\\ of the sill, which
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  marks as `nugget_at_bound`, is the fit's lower bound and not an
  estimate, and as \\\tau^2\\ it would leave no penalty, so that \\C_p\\
  would descend to the ceiling whatever the field: \\C_p\\ does not take
  it. When the variogram gives no nugget (none was fitted, its nugget is
  at that bound, its fit has no minimum inside the ranges searched, its
  range is below the shortest lag, it found no spatial structure, or
  `sac` is a range alone), \\\tau^2\\ is Mallows' own choice instead:
  the residual mean square \\RSS / (m - L_m)\\ of the finest level whose
  cells hold at least two scored rows on average. It counts the
  structure within those cells as noise too, so on average it is no
  smaller than the nugget and errs towards fewer cells. When no level's
  cells hold that many, or the finest that does has a residual mean
  square of 0 to within the rounding of the values scored (as when each
  of its cells holds one location whose rows carry one value), `cp` is
  `NA` at every level, as it is for a constant response, which the
  profile warns about. Whatever its source, \\\tau^2\\ is held to the
  semivariance of nearest-neighbour pairs the variogram carries
  (`attr(sac, "short_lag_semivariance")`, see
  [`sac_nugget`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md)):
  a semivariance does not fall as the lag shrinks, so it bounds the
  nugget from above, and a larger \\\tau^2\\ is contradicted by the
  closest pairs and stops \\C_p\\ early. A semivariance of 0 holds
  nothing. When \\\tau^2\\ is not the nugget as it stands, a warning of
  class `spatialkit_cp_noise` says where it came from, and
  `attr(, "cp_noise")` records the source and, when the semivariance
  held it, the number it replaced. Its minimum is often a bound, the
  support ceiling when the nugget is a small share of the sill, and
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  then names the bound in `edge`. The same formula at one cell, with the
  total sum of squares as its RSS, is kept in `attr(, "one_cell")`, with
  the p-value at every level of the F test of that level's cell means
  against one mean. A level has to beat one cell:
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  returns `structured = FALSE` when \\C_p\\ at one cell is within `tol`
  of the optimum, or when the cell means at the optimum differ from one
  mean by no more than chance allows, because one mean for the whole
  layer then predicts as well as the cells and the count is not
  identified (see
  [`select_resolution`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  for the rule).

- `moran_z`:

  A diagnostic of the cell means: the standardised deviate of Moran's I
  on the cell means, or with `predictor_vars` on the residuals of the
  cell means regressed on the cell-mean predictors, over the 8 nearest
  neighbours of each cell. It tests whether neighbouring cell means are
  independent. It does not measure what the cells leave unexplained, and
  it does not say how well they represent the response. On any
  autocorrelated response it grows as the cells get finer than the
  range, because neighbouring cells then resemble each other, and it has
  almost no power below about 20 cells, where the 8-neighbour graph
  between cells is close to complete and Moran's I sits near \\-1/(L -
  1)\\ whatever the data. The level with the smallest \\\|z\|\\ is
  therefore among the coarsest it is computed at on a structured
  response, whatever the layer and its ceiling, and arbitrary on a
  response with no structure, where \\z\\ is calibrated and flat in
  \\L\\.
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  refuses to read it as a criterion unless asked with
  `diagnostic = TRUE`, and neither
  [`summary()`](https://rdrr.io/r/base/summary.html) nor the functions
  that take a profile as a cell count read it by default. `NA` at nine
  cells or fewer.

- `reliability`:

  The between-cell signal's share of the spread in the cell means, from
  the fitted variogram alone via Krige's additivity relation (Cressie
  1996), for square cells of the level's average area holding the
  level's average share of the layer's points with a response (larger is
  better). This is the shrinkage factor of Fay and Herriot (1979). The
  curve has an interior optimum, a broad one, which can lie below the
  range floor.
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  reads it at the levels at or above the floor (its `range_floor`
  argument), and what it returns is then often the floor, which is a
  level of every ladder the profile builds: the coarsest count the range
  allows, with the range choosing rather than the curve. Report the flat
  region with the count. The domain term is taken over the convex hull
  the area is measured on, so rotating the layer does not move it. `NA`
  without a usable variogram, and when the variogram's range was
  rejected (see `sac`).

- `loo`:

  With an `estimate`: the leave-one-out mean squared error of that cell
  value (lower is better). It estimates the squared error of predicting
  one further observation, at a location drawn the way the data's
  locations were, by the value of the cell it falls in. That error
  depends on the cell value, so the column is computed for the estimate
  the profile is given and for no other. Each scored row of the layer is
  predicted by the value of its cell computed without it, in the cell of
  its nearest centre, which is the cell the seeds will build. Under
  `"mean"` a row in a cell of two or more is predicted by the mean of
  the others, and a row alone in its cell by the mean of the cell of its
  nearest other scored row. Under `"pooled"` the cell's mean and its own
  weight are computed again without the row; the spread inside the
  cells, the spread between a cell and its target and the targets are
  held; and a cell left empty takes its target. It needs no variogram.
  It describes locations among the observations and says nothing of a
  location away from them: it is a criterion for comparing levels where
  data are near, not an estimate of the error of a map. `loo_se` is the
  standard error, over rows, of that level's mean squared error. It is
  not the standard error of a difference between two levels, whose
  errors come from the same rows. `scored_per_cell` counts the rows
  `loo` read, every scored row of the layer, over the cells that hold
  one; `cell_resp_mean` counts the rows the other response criteria
  read, which are the subsample's when the cells are fitted to one. For
  pooled values its minimum tends to be the finest size offered (see
  "How this was checked"), so a size sets the count, and
  `select_resolution(criterion = "size")` reads it.

`cp` and `reliability` answer different questions: how well the cells
represent the field, and whether the cell values are distinguishable
from noise. They can disagree, and the profile shows both.

## Conditions

What changes the meaning of a result is raised as an R warning with a
class, so [`tryCatch()`](https://rdrr.io/r/base/conditions.html),
[`withCallingHandlers()`](https://rdrr.io/r/base/conditions.html) and
`options(warn = 2)` see it without matching message text; the rest is
logged (see
[`spatialkit_quiet`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)).
These are the classes this function raises;
[spatialkit-conditions](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-conditions.md)
lists every class the package raises.

- `spatialkit_not_supported`:

  `supported` is `FALSE`, because the range floor exceeds the ceiling or
  because fewer than two cells meet `min_cell_n`.

- `spatialkit_response_ceiling`:

  Rows that are not scored (no response, or with `predictor_vars` a
  predictor missing) lowered the ceiling. The message names the bound in
  force, `min_cell_n` or, on a subsample, `sample_n`; the condition
  carries `n`, `n_response`, `ceiling`, `ceiling_from`, `n_sample` and
  `n_scored`.

- `spatialkit_constant_response`:

  The scored variable is constant to within rounding (a response the
  predictors fit exactly leaves such residuals), so `cp`, `moran_i` and
  `moran_z` are `NA` at every level. On a profile with an `estimate`
  whose rows are constant too, so are `loo` and `loo_se`, and for pooled
  values `own_weight`; the message names them.
  [`summary()`](https://rdrr.io/r/base/summary.html) of such a profile
  stops with an error that says the response is constant. Both tests
  read the raw values and, where those show no spread, the values
  centred on their median and scaled, so 1e15 plus a spread of 1 is not
  constant while values that differ only by their rounding (0.1 beside
  0.3 / 3) are; residuals are read against the rounding of the
  regression and of the response, and where the residuals of the raw
  values show no spread, they are those of a fit to the response less
  its median. The test for `cp` reads the scored rows of the subsample
  and the estimate's every row of the layer (of its selection half under
  `select_on = "split"`) with a value: where the first finds the
  response constant and the second does not, the message says that `loo`
  is scored, and the print does not call the response a constant.

- `spatialkit_empty_cells`:

  Under `partition = "density"`: a centre holds no point, so there are
  fewer occupied cells than the count. A centre that loses its points
  stays where it is, and the start that is kept is the one with the
  lowest sum, whether or not all its centres hold points. It happens at
  high counts, mostly from `init = "random"`. The condition carries
  `level`, the count, and `empty`.

- `spatialkit_not_converged`:

  The start that was kept was cut off at the iteration limit, so its
  centres are those of the last step and a longer run can still move
  points between cells. Not raised for the kind `"kmeans"` from
  k-means++ starts under `iter_max = NULL`. The condition carries
  `level` (every level it happened at), `iter` and `iter_max`. Raise
  `iter_max`.

- `spatialkit_ladder_unused`:

  A numeric `obs_per_cell` sets the ladder, and the caller also supplied
  `n_levels`, `min_cell_n`, `sample_n` or `range_floor = TRUE`, which
  are then not read. The condition carries `unused` (the argument names)
  and `obs_per_cell`.

- `spatialkit_size_not_reached`:

  On a ladder of sizes: a size does not give the count it asks for.
  Either a size given as a number gives fewer than two cells while a
  finer one does not, and its level is 2 (`reason` is
  `"fewer than two cells"`); or a size asks for more cells than k-means
  can place, and its level is the most that can (`reason` is
  `"more cells than can be placed"`). It is raised once for each reason,
  naming every size held for it. The condition carries `obs_per_cell`
  and `wanted` (one element for each such size, in the order of the
  ladder), `used` and `reason`. When the finest size gives fewer than
  two cells the warning is `spatialkit_not_supported`, with
  `obs_per_cell` among its fields and `min_cell_n` `NA`; and
  `spatialkit_response_ceiling` is not raised on such a ladder, since a
  size already counts the rows with a response.

- `spatialkit_repeated_locations`:

  With an `estimate`: scored rows share a location. Each row is left out
  alone and the others at its location stay, so `loo` then describes a
  further observation at a place already observed. Raised once, after
  the levels are scored, and not when `loo` is `NA` at every level. The
  condition carries `n` (the rows that share a location), `total` (the
  rows `loo` read) and `column`.

- `spatialkit_no_scored_rows`:

  Under `select_on = "split"`: no row of the selection half has a
  response, and the rows of the layer with one are all in the estimation
  half. `rss`, `cp`, `moran_i` and `moran_z`, which read the selection
  half, are `NA` at every level. With an `estimate`, whose check that it
  has rows to read is made on the layer, `loo` and the estimate's other
  columns read no row and are `NA` at every level too. With
  `predictor_vars` a row counts when it has a response and every
  predictor. Not raised when no row of the layer has a response. Raised
  once, after the levels are scored. The condition carries `reason`
  (`"no response in the selection half"`), `n` (0, the rows of the
  selection half scored), `selection` (the rows of the selection half),
  `n_response` (the rows of the layer with a response) and `column`.
  Raised too, with or without the split, when rows of the layer have a
  response and none of them every predictor: then `reason` is
  `"no row with a response and every predictor"`, `n_response` the rows
  with a finite response, and `selection` `NA` without the split.

- `spatialkit_pool_unavailable`:

  With an `estimate`: the values of the response lie beyond what the
  arithmetic of the estimate holds, values so large, or so close to
  their median, that the squares of many such values can overflow or
  lose their precision, so `loo`, `loo_se` and `own_weight` are `NA` at
  every level. Where the test for `cp` finds such a response constant,
  `cp`, `moran_i` and `moran_z` are `NA` for the same reason, and this
  message names them in place of `spatialkit_constant_response`. The
  message gives the bounds; the response rescaled by a power of ten is
  scored. Raised once, after the levels are scored. The condition
  carries `reason` (`"values out of range"`) and `columns`, the
  response, as
  [`summarize_by_cell`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)`(estimate = "pooled")`
  raises the class for its pooled columns.

- `spatialkit_moran_skipped`:

  With a level of more than 5,000 cells, on any ladder, when FNN cannot
  be loaded or Matrix is not installed: `moran_i` and `moran_z` are `NA`
  at a level where more than 5,000 cells hold a scored row, and every
  other column is computed. Raised once, before the first fit, when such
  a level can occur: the cells that hold a scored row are counted after
  a level's fit, so where rows without a response leave fewer of them
  the two columns are computed after all. The condition carries `limit`,
  `levels` (the levels above the limit) and `missing`.

- `spatialkit_variogram_unusable`:

  The variogram, the one estimated or `sac`, gives no usable range
  (`rejected_reason` of
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)),
  so `reliability` is `NA`, and the message says whether `cp` may still
  take its nugget. The condition carries `reason`, `converged` (whether
  the variogram model converged; `NA` when it does not say), `nugget`
  and `nugget_used`.

- `spatialkit_cp_noise`:

  `cp` does not take the variogram's nugget as it stands. Its noise
  variance is the residual mean square of the finest level (the
  variogram gives no nugget, or its nugget is at the fit's lower bound),
  or the semivariance of nearest-neighbour pairs, which holds the nugget
  or that residual mean square where they are larger; or there is none
  and `cp` is `NA` at every level: the finest level whose cells hold two
  scored rows on average has a residual mean square of 0 to within
  rounding, or, with the nugget at its bound, no level's cells hold that
  many. The condition carries `source`, `value` and `level` (as
  `attr(, "cp_noise")` holds them; with no noise variance `source` and
  `value` are `NA`, and `level` is the level whose residual mean square
  is 0 to within rounding, or `NA`), `cap` (the semivariance, `NA` when
  the variogram carries none), `uncapped_source` and `uncapped_value`
  (what `cp` would take without the hold, and with no noise variance
  what it refused: the nugget at its bound, or that residual mean
  square).

The other warnings of this function (points dropped, a variogram of
another variable) are R warnings without a class; rows left out for a
missing response, a rising WSS curve and a subsample too small to reach
the floor are logged only.

One refusal is an error with a class, `spatialkit_work_limit`:
`resolution_profile()` stops before a step whose count of comparisons is
above a limit. The condition carries `work` (the count), `limit` and
`what` (the step). On a ladder that `obs_per_cell` sets, the count of
cells grows with the rows, so two costs are bounded before the variogram
is estimated and before any fit. The fits: the points times the cells of
the finest level times `nstart` may not pass
`getOption("spatialkit.fit_work_max", 1e10)` (`what = "partition fit"`;
the condition also carries `kind`, `init`, `nstart` and `fnn`). That
product counts the step in which every point is compared with every
centre, which `partition = "kmeans"` has in every pass and in a
k-means++ start, and `partition = "density"` has in a k-means++ start.
Without FNN the density kind compares every point with every centre in
every pass, and the product is multiplied by `iter_max`; with FNN and
`init = "random"` it has no such step and nothing is counted. The
density kind's weights: without FNN they come from a search that
compares every pair of distinct locations, and the square of their
number may not pass `getOption("spatialkit.exact_pairs_max", 2e9)`
(`what = "density weights"`); `density_power = 0` has no weights and no
such search. With an `estimate` and without FNN, on any ladder, the same
limit bounds the three exact searches the estimate's columns make: each
scored row's nearest centre, each cell's nearest cells, and, under
`"mean"`, the nearest other location of a row alone in its cell. The
result of every search is the same with and without FNN.

The warnings of the pooled values and of the count record, which
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
raises, are listed on its page.

## How this was checked

The sizes for pooled values were compared in a simulation of 200 layers,
each 1,500 clustered observations of a field with covariance
`exp(-d / 100)` plus noise on a square of side 1,000
([`dev/sim-pooled-values.R`](https://github.com/elkronos/gis_modeling_toolkit/blob/main/dev/sim-pooled-values.R)
in the package's GitHub repository;
[`summarize_by_cell`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
describes it), profiled at the sizes 18, 12, 9, 6, 4.5, 3 and 2 with
`partition = "density"`. The lowest `loo` of the pooled value was at the
size 2 on 176 layers, at 3 on 22 and at 4.5 on 2. That of the plain mean
was at 2 on none: it was at 6 on 86 layers, at 4.5 on 51, at 9 on 49, at
12 on 12, and at 18 and at 3 on one each. These figures are of those
seven sizes and say nothing of another ladder.

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
[`summary.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.resolution_profile.md)
for every criterion's pick side by side;
[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
for a count from the coordinates alone;
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
and
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
which accept the profile or a
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
result directly as `approx_n_cells` and `n`.

Other resolution:
[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md),
[`plot.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.resolution_profile.md),
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md),
[`summary.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.resolution_profile.md)

## Examples

``` r
library(sf)
# An exponential field with range parameter 200 (true effective range
# 600 m) on a 1 km square, with a nugget of 0.6 on a unit sill.
set.seed(4)
n <- 400
xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
D  <- as.matrix(dist(xy))
xy$z <- as.numeric(t(chol(exp(-D / 200) + diag(0.6, n))) %*% rnorm(n))
pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)

# The default profile: one level at 4.5 rows per cell, scored for
# pooled values, in density-weighted cells, with no variogram.
prof1 <- resolution_profile(pts, response_var = "z")
prof1
#> Resolution profile: 1 level on 400 points
#>   ladder      : 89 cells (obs_per_cell = 4.5 on the 400 rows with a response; centres fitted to every row; floor 2 (no range))
#>   variogram   : none fitted: a size ladder fits one only when `sac` is given (cp and reliability are NA)
#>   scored on   : response; 3 density-weighted k-means++ starts (power 1, 8 neighbours) per level
#>   one cell    : rss 490 (one mean for the layer; a level has to beat it)
#>   estimate    : pooled toward the 8 nearest cells; loo on 400 rows of the layer; one cell: loo 1.23
#> 
#>  levels    wss wss_spread elbow cell_n_min cell_n_median cell_diam_median rss
#>      89 502000     0.0698    NA          1             4               64 255
#>  cp moran_i moran_z reliability cells_resp cell_resp_mean scored_per_cell   loo
#>  NA   0.389    8.28          NA         89           4.49            4.49 0.903
#>  loo_se share_alone own_weight
#>  0.0645        0.01      0.397
select_resolution(prof1)   # the count the size gives
#> Resolution by size: 89 cells
#>   size        : 4.5 rows per cell, the default (obs_per_cell); 400 rows / 89 cells = 4.49
#>   cells       : density-weighted k-means++ starts (power 1, 8 neighbours)
#>   estimate    : pooled toward the 8 nearest cells
#>   at this size: 1.0% of the rows are alone in their cell; median own weight 0.40
#>   loo         : 0.903 (one cell for the layer: 1.23)
#>   note        : the count is set by obs_per_cell = 4.5, the default; no
#>                 criterion chose it.

# A ladder of counts for plain means: n_levels asks for the ladder and
# estimate = NULL scores no cell value, so cp is the rule that goes
# with it.
prof <- resolution_profile(pts, response_var = "z", n_levels = 12,
                           estimate = NULL)
prof
#> Resolution profile: 12 levels on 400 points
#>   ladder      : 2 to 44 cells (floor 4 from range 516.8, not applied, ceiling 44 from min_cell_n = 9)
#>   variogram   : nugget 0.582, partial sill 0.819, range 516.8
#>   scored on   : response; 3 density-weighted k-means++ starts (power 1, 8 neighbours) per level
#>   one cell    : rss 490, cp 1.23 (one mean for the layer; a level has to beat it)
#>   elbow       : none; the WSS curve falls as it does with no cluster structure
#> 
#>  levels      wss wss_spread elbow cell_n_min cell_n_median cell_diam_median rss
#>       2 41800000   0.000000    NA        182         200.0              646 482
#>       3 25400000   0.268000    NA        122         134.0              473 480
#>       4 16900000   0.000112    NA         77         106.0              410 469
#>       5 14100000   0.048400    NA         65          74.0              373 445
#>       6 11600000   0.125000    NA         53          68.5              347 432
#>       8  8550000   0.057700    NA         39          51.0              288 413
#>      11  5950000   0.084200    NA         28          36.0              239 394
#>      14  4430000   0.084100    NA         19          28.0              207 367
#>      19  3360000   0.028100    NA         11          22.0              168 343
#>      25  2530000   0.178000    NA          7          16.0              155 338
#>      33  1790000   0.110000    NA          7          12.0              130 348
#>      44  1250000   0.063900    NA          3           9.0              106 314
#>     cp moran_i moran_z reliability cells_resp cell_resp_mean
#>  1.210      NA      NA       0.839          2         200.00
#>  1.210      NA      NA       0.861          3         133.00
#>  1.180      NA      NA       0.863          4         100.00
#>  1.130      NA      NA       0.860          5          80.00
#>  1.100      NA      NA       0.855          6          66.70
#>  1.060      NA      NA       0.844          8          50.00
#>  1.020 -0.0248   1.400       0.827         11          36.40
#>  0.957 -0.0668   0.137       0.810         14          28.60
#>  0.914  0.0123   0.915       0.784         19          21.10
#>  0.918  0.0142   0.743       0.756         25          16.00
#>  0.965  0.1870   3.110       0.722         33          12.10
#>  0.914  0.2490   4.240       0.681         44           9.09
# Three of the criteria, one row per level.
as.data.frame(prof)[, c("levels", "cp", "reliability", "moran_z")]
#>    levels        cp reliability   moran_z
#> 1       2 1.2097988   0.8393165        NA
#> 2       3 1.2097026   0.8607721        NA
#> 3       4 1.1829553   0.8626250        NA
#> 4       5 1.1273510   0.8597290        NA
#> 5       6 1.0976257   0.8551213        NA
#> 6       8 1.0550201   0.8441615        NA
#> 7      11 1.0158319   0.8269075 1.4018098
#> 8      14 0.9569884   0.8101464 0.1367780
#> 9      19 0.9136533   0.7840438 0.9149028
#> 10     25 0.9183700   0.7556101 0.7431945
#> 11     33 0.9648884   0.7217720 3.1141763
#> 12     44 0.9139798   0.6811883 4.2441007
attr(prof, "bounds")$floor   # the range floor: reported, not applied
#> [1] 4
select_resolution(prof, criterion = "cp")
#> Resolution by cp: 19 cells
#>   flat region : 19 to 25, 44 (3 of 12 levels)
#>   cells       : density-weighted k-means++ starts (power 1, 8 neighbours)
```
