# Score every candidate number of cells on several criteria at once

A tessellation's resolution is the number of cells the points are cut
into, and no single number settles it: the within-cluster sum of squares
of the coordinates has an elbow, a piecewise-constant approximation of
the response has a Mallows \\C_p\\, the residual autocorrelation of the
cell means has a \\z\\, the cell means have a reliability, and a cell
value has a leave-one-out error. On a ladder of counts this function
computes all of them at every level the data bound, and returns the
table, so that the level chosen, by
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
```

## Arguments

- data_sf:

  An sf object of points (other geometries are reduced to representative
  points). Features with empty or non-finite coordinates are dropped
  with a warning.

- response_var:

  Optional response column name (numeric or logical). Enables `cp` and
  `moran_z`, and, when it is numeric, the pooled estimate the default
  profile scores by `loo`. A variogram estimated from it also gives the
  range floor the bounds report and `reliability`; it is estimated on a
  ladder of counts or when `sac` asks for it, and not on the default
  ladder of one size, so `cp` and `reliability` are `NA` there. Rows
  where it, or a predictor, is missing or non-finite stay in the
  geometry and are left out of the OLS fit, the RSS, `cp`, `moran_z` and
  the count the support ceiling divides by `min_cell_n`; a logged
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
  many centres as points).

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
  `sample_n` except through the fits.

- nstart:

  Starts per level; the best is kept. They are k-means++ starts unless
  `init` says otherwise. `NULL` (default) is the kind's own number: 3
  for the kind `"density"` and 25 for the kind `"kmeans"`, the number
  earlier versions ran. A number is run under either `partition`:
  `nstart = 25` gives the density kind the starts it had when 25 was the
  default. `attr(, "nstart")` holds the number that was run.

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
  since that needs the range; one whose model did not converge (its fit
  has no minimum inside the ranges searched, so its nugget and sill are
  not fitted values), or whose range is below the shortest lag fitted (a
  structure that cannot be told from a nugget, so the nugget is not
  identified either), gives neither, and `cp` takes its noise variance
  from the finest level (see `cp` above). Whether the model converged is
  read from the model
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
  noise variance from the finest level. A sac refused because the
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
  `sac = NULL`, the default of earlier versions, estimates the variogram
  on a ladder of sizes too, and gives `cp` and `reliability` there.

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
  it and the floor is below the ceiling, as the default did up to
  version 2.0.0; an estimated range then sets the fewest cells any
  criterion can choose (see "The ladder and its bounds").

- partition:

  What the centres of the regions minimise, at every level of the
  ladder. `"density"` (default): the sum of squared distances from the
  points to their nearest centre, each point weighted by its local
  density to the power `density_power`; compact cells that, at the
  default power, hold about the same number of points each once there
  are more cells than dense places, and counts about as uneven as
  `"kmeans"` gives below that. `"kmeans"`: the same sum unweighted;
  compact cells whose counts are several times higher in dense places
  than in sparse ones. `"kmeans"` was the default of earlier versions,
  and with `nstart = NULL` it runs their 25 starts. See "The partition".

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
  which reports nothing, as it always has.

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
  otherwise. `NULL`, the default of earlier versions, scores none and
  gives the profile without the five columns. `"mean"` and `"pooled"`
  are errors without a `response_var`, with `predictor_vars` (the
  profile would then score residuals, and the summary reports cell
  values of the response), with a logical response, and with a response
  that has no finite value.

- obs_per_cell:

  Rows per cell: `"auto"` (default; see below); `NULL` for the ladder
  described under "The ladder and its bounds", the default of earlier
  versions, whatever the other arguments are; or one number or several,
  each at least 1 (fractions are allowed), for a ladder of sizes. The
  level of a size \\s\\ is `max(2, round(n / s))` cells, with \\n\\ the
  rows that have a response (and, with `predictor_vars`, every
  predictor), or every point when there is no `response_var` or no row
  has a response. Under `select_on = "split"` \\n\\ counts the rows of
  both halves, since the cells are built on every point. On a ladder of
  sizes the centres are fitted to every row, no subsample is drawn, and
  `n_levels`, `min_cell_n`, `sample_n` and `range_floor` are not read:
  supplying `n_levels`, `min_cell_n`, `sample_n` or `range_floor = TRUE`
  beside a numeric `obs_per_cell` is a warning of class
  `spatialkit_ladder_unused`, and `levels` beside it is an error. A
  level is held to the most centres k-means can place (the distinct
  locations, and one short of the points). `"auto"` is `NULL` when
  `levels` is given or the caller supplies any of those four arguments;
  otherwise 4.5 when the estimate is `"pooled"`, the seven sizes 18, 12,
  9, 6, 4.5, 3 and 2 when it is `"mean"`, and `NULL` without an
  estimate. The cost of a ladder of sizes is bounded before any fit (see
  "Conditions").

## Value

A data.frame of class `resolution_profile` with one row per level and
columns `levels`, `wss` (the within-cell sum of squares of the
coordinates, every point counted once, under either `partition`),
`wss_spread` (relative spread of WSS across the restarts; under
`partition = "density"`, of the weighted sum the starts minimise),
`elbow`, `cell_n_min`, `cell_n_median`, `cell_diam_median` (twice the
median RMS radius of the cells, in coordinate units: about 0.8 of the
side of a square cell of the same area, and less for finer cells, so a
size to compare levels by rather than a width), `rss`, `cp`, `moran_i`,
`moran_z`, `reliability`, `cells_resp` (the cells that hold at least one
scored row) and `cell_resp_mean` (the scored rows per such cell; compare
it with `min_cell_n`). Like the other columns read off the fitted cells,
the last two count rows of the subsample when the layer is larger than
`sample_n`: the cells of the layer then hold about
`n_response / n_scored` times as many scored rows as `cell_resp_mean`
shows. Columns a missing input leaves undefined are `NA`. Attributes:
`bounds` (a list with `floor`, `ceiling`, `ceiling_from`
(`"min_cell_n"`, `"min_cell_n (rows with a response)"`,
`"distinct locations"`, `"sample_n"`, `"fewer than two supported cells"`
or `"fewer than two distinct locations"`, whichever bound it),
`supported`, `area`, `range`, `n` (the points in the layer), `n_sample`
(the points the k-means fits ran on), `n_distinct`, `min_cell_n`,
`range_floor`, `n_response` (the rows of the layer that carry a response
and, with `predictor_vars`, every predictor, which the support ceiling
counts; `NA` without a `response_var`), `n_scored` (the rows of the
subsample the response criteria read; `NA` without a `response_var`),
`levels_given` (whether `levels` set the ladder) and `predictors`
(whether `predictor_vars` were given)), `variogram` (a list with
`nugget`, `psill`, `range` (`NA` when rejected), `model` and `detrended`
(whether the sac says it is a variogram of residuals; `NA` when it does
not say); `NULL` when none was usable), `variable` (`"response"`,
`"residuals"` or `NA`), `wss_bumps`, `nstart`, `sac` (the range object
used) and, with `select_on = "split"`, `split` (a `spatialkit_split`:
`selection` and `estimation`, integer row positions in `data_sf`, with
the `method` and `seed` that made them), `cp_noise` (a list with
`source`, `"variogram nugget"` or `"finest-level residual mean square"`,
`value`, the \\\tau^2\\ used, and `level`, the level it was read at or
`NA`; `NULL` when `cp` has neither), `one_cell` (what one cell for the
whole layer scores, the reference a level has to beat: a list with `n`,
the scored rows, `rss`, their total sum of squares, `cp`, the same
formula at one cell or `NA` when `cp` has no noise variance, `constant`,
whether the scored variable is constant to within rounding, `cp_levels`,
the `cp` column it belongs to, and `f_p`, named by level, the p-value of
the F test of that level's cell means against one mean, \\(TSS - RSS) /
(L - 1)\\ over \\RSS / (m - L)\\ on the \\L\\ cells that hold a scored
row and the \\m\\ scored rows; `NULL` without a response), `centres` (a
list named by level of the two-column matrices of the scored partitions'
centres), `centres_crs` (the CRS they are in) and `hull` (the convex
hull of the points profiled, an `sfc` in that CRS). A k-means partition
is the Voronoi partition of its centres, so these seed the very cells
each level was scored on:
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
returns the chosen level's as `$seeds`, and
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
and
[`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)
use them when given the selection or the profile. The centres carry the
hull as their `"cover"` attribute, which those two functions hand on so
that Voronoi cells built on the seeds reach the edge of the points. A
profile built with `partition = "density"` or `init = "random"` also
carries `partition`: a list with `kind`, `init`, `density_power` and
`density_k` (both `NA` for `kind = "kmeans"`); for `kind = "density"`,
`objective`, the weighted sum the centres minimise, named by level;
`iter_max`, the iteration limit; and `iter` and `converged`, the
iterations each level's kept start took and whether it stopped before
the limit, named by level. The selection carries the first four as
`$partition`, its seeds and the seeds the seeding functions return as
`attr(, "partition")`, and a Voronoi tessellation built on them as
`params$partition`. A profile of the kind `"kmeans"` from k-means++
starts has no such attribute, and neither has anything made from it.

With an `estimate` the frame has five more columns, after
`cell_resp_mean`: `scored_per_cell` (the rows `loo` read over the cells
that hold one; the layer's counterpart of `cell_resp_mean`, which counts
the subsample's rows), `loo` (the leave-one-out mean squared error of
the estimate; lower is better), `loo_se` (the standard error over rows
of that level's mean squared error; it is not the standard error of a
difference between two levels, whose errors come from the same rows),
`share_alone` (the share of the scored rows that are the only scored row
of their cell) and `own_weight` (the median over cells of the weight a
pooled value gives the cell's own mean; `NA` under `"mean"`). `loo`,
`loo_se` and `own_weight` are `NA` at a level where no cell holds two
scored rows, and at every level when the response is constant over the
rows `loo` reads. It also carries `estimate`: a list with `value`
(`"mean"` or `"pooled"`), `neighbours` (8, or `NA` under `"mean"`),
`response` (the column), `n` (the rows `loo` read), `rows` (`"all"` or
`"selection half"`), `features` (`"points"`, or `"other"` when the layer
held features that were reduced to points), `seed`, `constant` (whether
the response is constant over those rows; `FALSE` for fewer than two),
`loo_one_cell` and `loo_one_cell_se` (the same error with one mean for
all those rows, the reference a level has to beat), `loo_levels` (the
`loo` column those references belong to, named by level), and, named by
level, `sigma2` (the spread of the rows about their cell's mean), `tau2`
(the mean squared distance between a cell's level and its target; `NA`
under `"mean"`), `cells` (the cells that hold a scored row) and `f_p`
(the p-value of the F test of those cells' means against one mean).

On a ladder that `obs_per_cell` sets, `bounds` describes that ladder:
`ceiling` is the finest level; `ceiling_from` is `"obs_per_cell"`,
`"distinct locations"`, `"fewer than two cells of the size"` or
`"fewer than two distinct locations"`; `supported` is `FALSE` when the
finest size gives fewer than two cells or the range floor is above the
finest level; `levels_given` is `TRUE`; `n_sample` is every point;
`min_cell_n` holds the argument, which is not read; and two fields are
added, `obs_per_cell` (the sizes, largest first) and `size_levels` (the
level of each, named by size; `NA` for a size of `"auto"` that gives
fewer than two cells and was left out).

## Details

Called with a numeric response and no predictors and none of `levels`,
`n_levels`, `min_cell_n`, `sample_n` or `range_floor`, the profile is
one level: 4.5 rows per cell (`obs_per_cell = "auto"`), scored for
pooled values (`estimate = "auto"`) by the leave-one-out error `loo`,
with no variogram fitted (`sac = "auto"`), in regions of the density
kind (`partition = "density"`) kept from 3 k-means++ starts
(`nstart = NULL`). `cp` and `reliability` need the variogram and `elbow`
needs a curve, so the three are `NA` there; `moran_z` is computed.
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
reads the count off it by the size. Any of those five ladder arguments
gives the ladder of counts described below, scored for pooled values as
well.
`estimate = NULL, obs_per_cell = NULL, sac = NULL, partition = "kmeans"`
gives the profile of earlier versions: the ladder, its columns and its
k-means regions from 25 starts.

## The partition

The regions at every level are the Voronoi cells of a set of centres: a
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

Two simulations show both ends. The measure is the coefficient of
variation of the cell counts (their standard deviation over their mean,
0 when every cell holds the same number), as medians over 20 layers. On
1,200 points, 900 of them in four clusters over a sparse background, cut
into 60 cells, it was 1.2 under `"kmeans"`, where the fullest cell held
42 times the points of the smallest, and 0.29 under `"density"`, where
it held 4 times. On 2,000 points with 28 dense places (half of the
points in 3 large clusters, 30 percent in 25 small ones, the rest spread
out) it was 0.90 under `"kmeans"` and 0.97 under `"density"` at 10 cells
and 1.08 and 1.01 at 20 cells; at 60 cells it was 0.95 and 0.31, and at
130 cells 1.05 and 0.43. A ladder starts at 2 cells, so its low levels
are usually below the number of dense places, and the counts in
`cell_n_min` and `cell_n_median` show where the two kinds part.

The density kind is the default. The default size, 4.5 rows per cell,
puts the count above the number of dense places on a layer whose
clusters hold more than a few observations each, so the default kind and
the default count go together. A coarse ladder the caller sets (with
`levels`, large sizes, or the low levels of the automatic ladder) may
not, and at those levels the density kind does not even the counts.
`partition = "kmeans"` gives the k-means regions, with 25 k-means++
starts, as earlier versions built them.

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
start with the lowest sum is kept, and a kept start that was cut off at
the limit is reported (see "Conditions"): its centres are those of the
last step. Large layers can need more than 100 steps: on a simulated
layer of 30,000 points no start of six at 100 and 300 cells took more
than 68, and on one of 60,000 points three of six took 102 to 152.

The density kind depends on the points alone. It works on the distinct
locations in coordinate order, so the order of the rows does not matter.
Its two tolerances are fractions of the extent of the points (the
diagonal of their bounding box), so the unit of the coordinates does not
matter: the neighbour distance is not taken below a millionth of the
extent, and the iteration stops when no centre moves by \\10^{-12}\\ of
it. Points at one location count as that many points at it and share a
cell; the neighbour distance is measured to the `density_k`-th nearest
*other* location, so it is not zero where locations repeat. The
neighbour searches use FNN when it is installed. The result is the same
without it; the search then takes time that grows with the square of the
number of points, about thirty times as long for 10,000 points and 2,000
centres.

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
flattens. On simulated layers of 800 points in 3 to 8 clusters over a
sparse background an elbow was read on 24 of 30 under `"density"` and on
10 under `"kmeans"`; on 100 uniform layers, on none under either.

A profile of the density kind, or of either kind from random starts,
carries `attr(, "partition")`, and
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md),
the seeding functions and
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
hand that record on with the regions (see "Value"). A result with no
such record was built by plain k-means: from k-means++ starts when it
comes from a profile, and from random starts when a seeding function
fitted it.

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
response with no spatial structure can still be given one: over 100
simulated white-noise layers of 900 points the variogram returned a
range on 29, and started at the floor, 11 of those ladders began at 11
to 98 cells of a ceiling of 100 and \\C_p\\ took the floor or a level
beside it, where the best count was 2 to 8. `range_floor = TRUE` starts
the ladder at the floor, as the default did up to version 2.0.0. When
the floor exceeds the ceiling the data cannot support a tessellation
that respects their own correlation structure; that is reported as a
finding (a warning of class `spatialkit_not_supported`, and
`attr(x, "bounds")$supported` is `FALSE`) and the ladder runs from 2 to
the ceiling, so the profile still shows what each level costs. A floor
equal to the ceiling would leave one level and nothing to compare it
with, so that ladder starts at 2 as well. On a layer larger than
`sample_n` the ceiling is also held to half the subsample's scored rows
(two per cell, the least a fitted cell can be scored on; half the
subsample without a response); `ceiling_from` is then `"sample_n"`, a
floor above that is logged with a request to raise `sample_n`, and it
does not make `supported` `FALSE`.

The floor moves with the range estimate, and the ceiling with \\n\\. To
compare profiles of similar data (the folds of a cross-validation, say),
pass the same `levels` to each. `reliability` is the one criterion that
still reads the floor:
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
holds it to the levels at or above the floor by default, so its pick can
differ between two such profiles by a factor of 10 or more when the
floor is supported on one and not on the other; pass
`range_floor = FALSE` there to read the whole ladder.

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
the end of the ladder would then choose the count. In a simulation of
200 layers, each 1,500 clustered observations of a field with covariance
`exp(-d / 100)` plus noise on a square of side 1,000
(`dev/sim-pooled-values.R` in the package's repository;
[`summarize_by_cell`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
describes it), profiled at the sizes 18, 12, 9, 6, 4.5, 3 and 2 with
`partition = "density"`, the lowest `loo` of the pooled value was at the
size 2 on 176 layers, at 3 on 22 and at 4.5 on 2. That of the plain mean
was at 2 on none: it was at 6 on 86 layers, at 4.5 on 51, at 9 on 49, at
12 on 12, and at 18 and at 3 on one each.
`select_resolution(criterion = "size")` therefore sets the count of
pooled values from a size, 4.5 unless another is given, and
`criterion = "loo"` reads the lowest error for plain means. These
figures are of those seven sizes and say nothing of another ladder.

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
`partition = "density"` with `init = "random"` and FNN installed, which
leaves a cell without a point more often (see `spatialkit_empty_cells`
under "Conditions"); or a higher limit. The limit counts the finest
level, and a ladder of several sizes fits every one of its levels. At
the default size the work grows with the square of the rows, and the
limit is set where the default call's fit takes about a minute: a call
just below it runs that long, and prints nothing while it runs.

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

## The criteria, and what each returns in practice

The four columns are not four interchangeable count rules.
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
reads `cp` by default on a profile without an estimate; on one with an
estimate it reads the size (pooled values) or `loo` (plain means).
`reliability` needs a variogram range and is `NA` without one; `elbow`
is `NA` unless the points cluster; and `moran_z` is a diagnostic and not
a count rule. With an `estimate` the profile has a fifth column, `loo`,
which is the last entry. Each entry below says what the column measures
and what it returns in practice.

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
  response. **In practice it is `NA` unless the points cluster**, and
  where they do it names the number of clusters, which is usually the
  first level or two of the ladder.

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
  the field; the profile warns. When the variogram gives no nugget (none
  was fitted, its fit has no minimum inside the ranges searched, its
  range is below the shortest lag, it found no spatial structure, or
  `sac` is a range alone), \\\tau^2\\ is Mallows' own choice instead:
  the residual mean square \\RSS / (m - L_m)\\ of the finest level whose
  cells hold at least two scored rows on average. It counts the
  structure within those cells as noise too, so on average it is no
  smaller than the nugget and errs towards fewer cells; the profile
  warns, and `attr(, "cp_noise")` records which was used. When no
  level's cells hold that many, `cp` is `NA` at every level, as it is
  for a constant response, which the profile warns about. **In practice
  its minimum is often a bound**, the support ceiling when the nugget is
  a small share of the sill, and
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
  for the rule and its measured error rates).

- `moran_z`:

  **A diagnostic of the cell means, not a count rule.** The standardised
  deviate of Moran's I on the cell means, or with `predictor_vars` on
  the residuals of the cell means regressed on the cell-mean predictors,
  over the 8 nearest neighbours of each cell. It tests whether
  neighbouring cell means are independent. It does not measure what the
  cells leave unexplained, and it does not say how well they represent
  the response. On any autocorrelated response it grows as the cells get
  finer than the range, because neighbouring cells then resemble each
  other, and it has almost no power below about 20 cells, where the
  8-neighbour graph between cells is close to complete and Moran's I
  sits near \\-1/(L - 1)\\ whatever the data. The level with the
  smallest \\\|z\|\\ is therefore among the coarsest it is computed at
  on a structured response, whatever the layer and its ceiling, and
  arbitrary on a response with no structure, where \\z\\ is calibrated
  and flat in \\L\\.
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
  curve has an interior optimum, and a broad one: validated against the
  empirical reliability of true block means on simulated fields, the
  analytic and empirical optima agreed to within a level or two where
  the empirical estimate was stable, and the band within 2 percent of
  the maximum spanned a factor of 3–6 in \\L\\. **That optimum usually
  lies below the range floor**: on 20 simulated fields of 900 points the
  curve read from the floor up was highest at the floor on 19.
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  therefore reads it at the levels at or above the floor (its
  `range_floor` argument), and what it returns is then the floor, which
  is a level of every ladder the profile builds: the coarsest count the
  range allows, with the range choosing and not the curve. Read the flat
  region, not the argmax. The domain term is taken over the convex hull
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
  read, which are the subsample's when the cells are fitted to one. **In
  the simulation described under "Observations per cell", on a ladder of
  the seven sizes 18 to 2, its minimum was at the finest size on 176 of
  200 layers for pooled values and on none for plain means**: for pooled
  values a size sets the count, and
  `select_resolution(criterion = "size")` reads it.

`cp` and `reliability` answer different questions: how well the cells
represent the field, and whether the cell values are distinguishable
from noise. They can disagree, and both are shown so the choice between
them is made knowingly.

## Conditions

What changes the meaning of a result is raised as an R warning with a
class, so [`tryCatch()`](https://rdrr.io/r/base/conditions.html),
[`withCallingHandlers()`](https://rdrr.io/r/base/conditions.html) and
`options(warn = 2)` see it without matching message text; the rest is
logged (see
[`spatialkit_quiet`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)).

- `spatialkit_not_supported`:

  `resolution_profile()`: `supported` is `FALSE`, because the range
  floor exceeds the ceiling or because fewer than two cells meet
  `min_cell_n`.

- `spatialkit_response_ceiling`:

  `resolution_profile()`: rows that are not scored (no response, or with
  `predictor_vars` a predictor missing) lowered the ceiling. The message
  names the bound in force, `min_cell_n` or, on a subsample, `sample_n`;
  the condition carries `n`, `n_response`, `ceiling`, `ceiling_from`,
  `n_sample` and `n_scored`.

- `spatialkit_constant_response`:

  `resolution_profile()`: the scored variable is constant to within
  rounding (a response the predictors fit exactly leaves such
  residuals), so `cp`, `moran_i` and `moran_z` are `NA` at every level.

- `spatialkit_empty_cells`:

  `resolution_profile()`,
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
  and
  [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)
  under `partition = "density"`: a centre holds no point, so there are
  fewer occupied cells than the count. A centre that loses its points
  stays where it is, and the start that is kept is the one with the
  lowest sum, whether or not all its centres hold points. It happens at
  high counts, mostly from `init = "random"`. On simulated layers of 800
  points (in 4 clusters, in 28, and uniform; 20 of each) no partition
  into 80 cells or fewer had such a centre. At 400 cells, from random
  starts with the best of 10 kept, 2, 5 and 1 of the 20 layers had one;
  from k-means++ starts none of 840 partitions at 40 to 600 cells had
  one, though it can happen. The condition carries `level`, the count,
  and `empty`.

- `spatialkit_not_converged`:

  `resolution_profile()`,
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
  and
  [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md),
  for every call except the plain k-means each has always run: the start
  that was kept was cut off at the iteration limit, so its centres are
  those of the last step and a longer run can still move points between
  cells. The condition carries `level` (the count; for a profile every
  level it happened at), `iter` and `iter_max`. Raise `iter_max`
  (`kmeans_iter` in
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)).

- `spatialkit_partition_unused`:

  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
  and
  [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md),
  given a profile or a
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  result together with a `partition`, an `init` or density settings
  other than those the profile was built with: the scored centres are
  returned, and what the call asked for is not used. The condition
  carries `kind` and `init`, those of the centres returned.

- `spatialkit_bound_chosen`:

  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
  [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)
  and
  [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md),
  given a profile or a
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  result whose optimum sits on a bound of the ladder: the bound chose
  the count. The condition carries `edge`, `criterion` and `n`.

- `spatialkit_no_structure`:

  The same functions, given a count whose `structured` is `FALSE`: one
  cell does as well, so the count is not identified. A count set by a
  size stands, since no criterion chose it, and the warning then says
  that the cell values will say little that one mean for the layer does
  not.
  [`summarize_by_cell`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)`(estimate = "pooled")`
  raises the class for its own cells.

- `spatialkit_no_elbow`:

  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md):
  the WSS curve has no elbow, or the ladder is too short to show one, so
  the count returned is the ladder's and not a finding. The condition
  carries `k`, `sag` and `structured`.

- `spatialkit_diagnostic_ranking`:

  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  under `criterion = "morans_i"` or `"combined"`: Moran's z of the cell
  means, a diagnostic and not a count rule, ordered the candidates
  returned. The condition carries `k`, `criterion`, `first_computable`,
  `window` and `upgraded`.

- `spatialkit_ladder_unused`:

  `resolution_profile()`: a numeric `obs_per_cell` sets the ladder, and
  the caller also supplied `n_levels`, `min_cell_n`, `sample_n` or
  `range_floor = TRUE`, which are then not read. The condition carries
  `unused` (the argument names) and `obs_per_cell`.

- `spatialkit_size_not_reached`:

  `resolution_profile()` on a ladder of sizes: a size does not give the
  count it asks for. Either a size given as a number gives fewer than
  two cells while a finer one does not, and its level is 2 (`reason` is
  `"fewer than two cells"`); or a size asks for more cells than k-means
  can place, and its level is the most that can (`reason` is
  `"more cells than can be placed"`). The condition carries
  `obs_per_cell`, `wanted`, `used` and `reason`. When the finest size
  gives fewer than two cells the warning is `spatialkit_not_supported`,
  with `obs_per_cell` among its fields and `min_cell_n` `NA`; and
  `spatialkit_response_ceiling` is not raised on such a ladder, since a
  size already counts the rows with a response.

- `spatialkit_repeated_locations`:

  `resolution_profile()` with an `estimate`: scored rows share a
  location. Each row is left out alone and the others at its location
  stay, so `loo` then describes a further observation at a place already
  observed. Raised once, after the levels are scored, and not when `loo`
  is `NA` at every level. The condition carries `n` (the rows that share
  a location), `total` (the rows `loo` read) and `column`.

- `spatialkit_moran_skipped`:

  `resolution_profile()` with a level of more than 5,000 cells, on any
  ladder, when FNN or Matrix is not installed: `moran_i` and `moran_z`
  are `NA` at a level where more than 5,000 cells hold a scored row, and
  every other column is computed. Raised once, before the first fit,
  when such a level can occur: the cells that hold a scored row are
  counted after a level's fit, so where rows without a response leave
  fewer of them the two columns are computed after all. The condition
  carries `limit`, `levels` (the levels above the limit) and `missing`.

Two refusals are errors with a class: `spatialkit_diagnostic_criterion`
(`select_resolution(criterion = "moran_z")` without `diagnostic = TRUE`)
and `spatialkit_no_default_criterion` (a profile passed as a cell count
whose `cp` is `NA` at every level while another criterion is finite).
The remaining warnings of these functions (points dropped, a variogram
that gives `cp` no nugget, a zero nugget, a variogram of another
variable) are R warnings without a class; rows left out for a missing
response, a rising WSS curve and a subsample too small to reach the
floor are logged only.

A third refusal is an error of class `spatialkit_work_limit`:
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
raises, are listed under "Conditions" on its page.

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
[`cell_values_at()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cell_values_at.md),
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
  # 600 m) on a 1 km square, with a nugget of 0.6 on a unit sill.
  set.seed(4)
  n <- 400
  xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
  D  <- as.matrix(dist(xy))
  xy$z <- as.numeric(t(chol(exp(-D / 200) + diag(0.6, n))) %*% rnorm(n))
  pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
  # n_levels asks for the ladder of counts.
  prof <- resolution_profile(pts, response_var = "z", n_levels = 12)
  print(prof)               # one row per level; print() because only the
                            # last value of a braced block is shown
  attr(prof, "bounds")$floor   # the range floor: reported, not applied
  print(select_resolution(prof, criterion = "cp"))

  # The default profile: one level at 4.5 rows per cell, scored for
  # pooled values, in density-weighted cells, with no variogram.
  prof1 <- resolution_profile(pts, response_var = "z")
  print(prof1)
  select_resolution(prof1)   # the count the size gives
}
#> Resolution profile: 12 levels on 400 points
#>   ladder      : 2 to 44 cells (floor 4 from range 517, not applied, ceiling 44 from min_cell_n = 9)
#>   variogram   : nugget 0.582, partial sill 0.819, range 517
#>   scored on   : response; 3 density-weighted k-means++ starts (power 1, 8 neighbours) per level; WSS rises at 0 step(s)
#>   one cell    : rss 490, cp 1.23 (one mean for the layer; a level has to beat it)
#>   estimate    : pooled toward the 8 nearest cells; loo on 400 rows of the layer; one cell: loo 1.23
#>   elbow       : none; the WSS curve falls as it does with no cluster structure
#> 
#>  levels      wss wss_spread elbow cell_n_min cell_n_median cell_diam_median rss
#>       2 41800000   0.00e+00    NA        182         200.0              646 482
#>       3 25400000   6.21e-05    NA        122         134.0              473 480
#>       4 16900000   2.80e-01    NA         78         105.0              410 469
#>       5 14400000   4.31e-02    NA         64          78.0              359 457
#>       6 11600000   8.99e-02    NA         56          68.0              348 436
#>       8  8320000   1.33e-01    NA         40          49.5              297 437
#>      11  5670000   1.07e-01    NA         30          34.0              229 381
#>      14  4480000   1.14e-01    NA         19          29.0              211 370
#>      19  3350000   3.71e-03    NA         11          20.0              175 368
#>      25  2300000   1.82e-01    NA          9          16.0              153 345
#>      33  1820000   1.08e-01    NA          6          12.0              133 334
#>      44  1260000   2.74e-02    NA          5           9.0              110 297
#>     cp moran_i moran_z reliability cells_resp cell_resp_mean scored_per_cell
#>  1.210      NA      NA       0.839          2         200.00          200.00
#>  1.210      NA      NA       0.861          3         133.00          133.00
#>  1.180      NA      NA       0.863          4         100.00          100.00
#>  1.160      NA      NA       0.860          5          80.00           80.00
#>  1.110      NA      NA       0.855          6          66.70           66.70
#>  1.120      NA      NA       0.844          8          50.00           50.00
#>  0.984 -0.0566   0.812       0.827         11          36.40           36.40
#>  0.965 -0.0429   0.496       0.810         14          28.60           28.60
#>  0.974  0.0195   1.020       0.784         19          21.10           21.10
#>  0.935  0.1240   2.240       0.756         25          16.00           16.00
#>  0.931  0.1610   2.710       0.722         33          12.10           12.10
#>  0.871  0.2500   4.280       0.681         44           9.09            9.09
#>    loo loo_se share_alone own_weight
#>  1.220 0.0987           0      0.930
#>  1.220 0.0984           0      0.843
#>  1.190 0.0972           0      0.884
#>  1.170 0.0891           0      0.890
#>  1.120 0.0848           0      0.917
#>  1.130 0.0911           0      0.871
#>  1.000 0.0716           0      0.904
#>  0.988 0.0695           0      0.905
#>  1.010 0.0706           0      0.858
#>  0.973 0.0697           0      0.842
#>  0.970 0.0706           0      0.790
#>  0.905 0.0631           0      0.747
#> Resolution by cp: 44 cells
#>   flat region : 44 (1 of 12 levels)
#>   cells       : density-weighted k-means++ starts (power 1, 8 neighbours)
#>   note        : the optimum is the support ceiling (n / min_cell_n); the
#>                 bound is choosing, not the criterion. Lower min_cell_n to see
#>                 whether the criterion keeps going.
#>   note        : cp scores plain means, and this profile's estimate is the
#>                 pooled value; criterion = "size" is the rule for it.
#> Resolution profile: 1 levels on 400 points
#>   ladder      : 89 cells (obs_per_cell = 4.5 on the 400 rows with a response; centres fitted to every row; floor 2 (no range))
#>   variogram   : none fitted: a size ladder fits one only when `sac` is given (cp and reliability are NA)
#>   scored on   : response; 3 density-weighted k-means++ starts (power 1, 8 neighbours) per level; WSS rises at 0 step(s)
#>   one cell    : rss 490 (one mean for the layer; a level has to beat it)
#>   estimate    : pooled toward the 8 nearest cells; loo on 400 rows of the layer; one cell: loo 1.23
#> 
#>  levels    wss wss_spread elbow cell_n_min cell_n_median cell_diam_median rss
#>      89 519000     0.0338    NA          1             4             62.3 253
#>  cp moran_i moran_z reliability cells_resp cell_resp_mean scored_per_cell   loo
#>  NA   0.402    8.53          NA         89           4.49            4.49 0.896
#>  loo_se share_alone own_weight
#>  0.0636       0.005      0.447
#> Resolution by size: 89 cells
#>   size        : 4.5 rows per cell asked (obs_per_cell); 400 rows / 89 cells = 4.49
#>   cells       : density-weighted k-means++ starts (power 1, 8 neighbours)
#>   estimate    : pooled toward the 8 nearest cells
#>   at this size: 0.5% of the rows are alone in their cell; median own weight 0.45
#>   loo         : 0.896 (one cell for the layer: 1.23)
#>   note        : the count is set by obs_per_cell = 4.5, the caller's choice;
#>                 no criterion chose it.
```
