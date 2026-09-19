# Create spatial cross-validation folds

Builds train/test splits using random K-fold, spatial block K-fold, or
buffered leave-one-out strategies.

## Usage

``` r
make_folds(
  points_sf,
  k,
  method = c("random_kfold", "block_kfold", "buffered_loo", "leave_location_out", "nndm"),
  seed = NULL,
  block_nx = NULL,
  block_ny = NULL,
  block_multiplier = 3,
  block_size = NULL,
  auto_range = FALSE,
  range_frac = 1,
  response_var = NULL,
  group_var = NULL,
  prediction_points = NULL,
  predictor_vars = NULL,
  boundary = NULL,
  buffer = NULL,
  min_train = 0.5,
  phi = NULL,
  drop_empty_blocks = TRUE,
  blocks = NULL,
  balance_tol = 3
)
```

## Arguments

- points_sf:

  An sf object. Any Z or M dimension is dropped before folding:
  [`sf::st_distance()`](https://r-spatial.github.io/sf/reference/geos_measures.html)
  uses every coordinate dimension, so an XYZ layer would otherwise have
  elevation folded into every buffer, block and neighbour distance.
  CRS-less points are aligned to a `boundary` or `prediction_points`
  that carries a CRS. They are reprojected when the coordinates look
  like lon/lat and otherwise stamped without reprojection, with a
  warning either way.

- k:

  Integer; number of folds. Must be a single whole number \>= 1. A
  fraction, `NA` or a vector is an error, because a non-integer used to
  truncate silently and leave the last rows in no test set at all. Not
  every method honours it. `"buffered_loo"` and `"nndm"` are
  leave-one-out schemes and always return `k = n` regardless of what was
  asked for; `"block_kfold"` lowers it when the grid yields fewer than
  `k` non-empty blocks, and `"leave_location_out"` lowers it when there
  are fewer than `k` distinct groups. Read the `k` element of the
  returned list, and do not assume the requested value. A reduction is
  written to the package log and raises no R warning, so
  `tryCatch(warning = )` will not see it and
  [`suppressWarnings()`](https://rdrr.io/r/base/warning.html) will not
  hide it.

- method:

  One of `"random_kfold"`, `"block_kfold"`, `"buffered_loo"`,
  `"leave_location_out"` or `"nndm"`. See **Details** for what each one
  does and when it is appropriate.

- seed:

  Optional integer RNG seed.

- block_nx, block_ny:

  Optional grid dimensions for block_kfold. Ignored when `block_size` or
  `auto_range` override them.

- block_multiplier:

  Numeric, default 3. When neither `block_size` nor
  `block_nx`/`block_ny` is given, the automatic grid aims for
  `block_multiplier * k` blocks over the extent (aspect-preserving), so
  each fold holds out about `block_multiplier` blocks. With 1, every
  fold is one contiguous region and the score depends heavily on which
  region each fold happened to get; with many, the blocks shrink towards
  single points and the scheme drifts back towards random k-fold. 3 is a
  compromise between those two, not a published constant. The block size
  that matters for leakage is the autocorrelation range, which is what
  `block_size` and `auto_range` control.

- block_size:

  Optional positive numeric minimum block edge length, **in the units of
  the CRS the folds are built in**. When supplied, grid dimensions are
  clamped so that every block is at least this wide and tall. Takes
  precedence over `block_nx`/`block_ny` and `block_multiplier`.

  Which CRS that is depends on the input. Projected input is used as it
  stands, so `block_size` is in your own CRS's units. Geographic
  (lon/lat) input is projected first by
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md),
  which picks a local UTM zone or, at wide extents, an equal-area
  projection. That is a CRS you did not choose, whose units are metres
  but whose identity varies with the data. `block_size` is then
  interpreted in *that* CRS. The CRS actually used is recorded in
  `params$crs` of the returned list; project the data yourself before
  calling if you want to fix the units in advance.

  A `block_size` in the wrong unit asks for an enormous grid, so a
  request above 1,000,000 blocks is refused with an error naming the
  grid dimensions, the extent and the CRS's units.

- auto_range:

  Logical. If `TRUE`, the spatial autocorrelation range is estimated via
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  (which fits directional variograms to account for anisotropy) and used
  as the minimum `block_size`. Requires `response_var`. An explicit
  `block_size` takes precedence. Default `FALSE`. Sizing blocks from the
  autocorrelation range is the recommendation of Roberts et al. (2017)
  and what blockCV (Valavi et al. 2019) automates. blockCV takes the
  fitted variogram's range *parameter* as the block size, whereas this
  uses the *effective* range
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  returns (three times that parameter for an exponential fit), so its
  blocks are larger than blockCV's from the same variogram.

- range_frac:

  Passed through to
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  when `auto_range = TRUE`. A fitted range beyond the longest lag the
  empirical variogram was fitted over is rejected as unidentified, and
  block sizing falls back to geometry, so the grid does not collapse to
  a single block. Default 1.0.

- response_var:

  Character(1) response column name. Required when `auto_range = TRUE`.

- group_var:

  Character(1) naming a column of `points_sf` that identifies the
  location each observation belongs to. Required for
  `method = "leave_location_out"`, which keeps every observation from a
  location together in the same fold. Repeated measurements at the same
  site otherwise get split across folds, and the model is scored partly
  on sites it has already seen, which random k-fold reports as excellent
  performance.

- prediction_points:

  Optional `sf` layer of the locations you actually intend to predict
  onto. Required for `method = "nndm"`. The grid from
  [`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md)
  is the natural choice; a non-POINT layer (grid cells, polygons) is
  reduced to representative points first, so the target distances are
  point-to-point. A point-to-polygon distance is zero for every cell
  that contains a training point, which pulls the target distribution
  towards zero and degenerates the CV towards plain leave-one-out.

- predictor_vars:

  Optional character vector of predictor column names. Passed to
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  for residual variogram estimation.

- boundary:

  Optional polygonal sf/sfc for block_kfold.

- buffer:

  Positive numeric distance for buffered_loo.

- min_train:

  For `method = "nndm"`: the smallest fraction of the data any fold's
  training set may be reduced to by neighbour exclusion. Default `0.5`,
  as in `CAST::nndm()`.

- phi:

  For `method = "nndm"`: the distance up to which the two
  nearest-neighbour distance distributions are matched, in the CRS the
  folds are built in; the exclusion never pushes a held-out point's
  nearest neighbour beyond it. In Mila et al. (2022), and in
  `CAST::nndm()`, \\\phi\\ is the autocorrelation range of the outcome:
  beyond it observations are effectively independent, so matching is
  unnecessary.
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  gives such a value. Default `NULL` = the largest
  prediction-to-training distance, which matches everywhere (`CAST`'s
  `phi = "max"`).

- drop_empty_blocks:

  Logical. Default TRUE.

- blocks:

  Optional polygon layer (`sf` or `sfc`, POLYGON or MULTIPOLYGON, at
  least two features) to use as the blocks of `"block_kfold"` in place
  of the grid this function would otherwise build: the `$cells` of a
  [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
  result, hexagons, watersheds, administrative units, the `$blocks` of
  `blockCV::cv_spatial()`. Each point takes the block that contains it,
  and the blocks are then assigned to folds exactly as grid cells are.
  `block_size`, `block_nx`/`block_ny`, `block_multiplier` and `boundary`
  have nothing to act on and are ignored (logged); `auto_range` only
  compares the estimated range against the blocks. See **Supplied
  blocks** below for the CRS, overlap and coverage rules. An error for
  any other `method`.

- balance_tol:

  A number of at least 1, default 3. For `"block_kfold"`: the ratio of
  the largest fold's point count to the smallest's above which the folds
  are reported as imbalanced, with a warning (an R condition, also
  logged) that names both counts. `Inf` disables the check. The value is
  the tolerance of a check, not a target the packing aims for: see
  **Fold balance** below for what the packing can and cannot do.
  `params$balance_ratio` carries the ratio achieved.

## Value

A list with method, k, folds, assignment, params. The `train`/`test`
elements of each fold contain `..row_id` values (equal to row positions
when the input has no pre-existing `..row_id` column), consistent with
the `assignment` tibble. The returned `k` is the number of folds
actually built, which is not always the `k` that was requested (see the
`k` argument above), and `length(folds)` always matches it.

For `"block_kfold"` the block design is returned with the folds.
`assignment` has a third column, `block_id`: the block each point fell
in, numbered as the rows of `params$blocks`, an `sf` layer of the block
polygons in the CRS the folds were built in. Those rows are **not** the
rows of the layer the blocks came from: `drop_empty_blocks = TRUE` (the
default) removes the blocks that hold no point and renumbers the rest,
so with supplied `blocks` a nine-polygon layer of which three are empty
comes back as six rows numbered 1 to 6. `params$blocks$source_row` is
the row each one came from in that original layer (a grid cell's index
in the full `grid_nx` by `grid_ny` grid, or the row of the `blocks`
argument), so `blocks[params$blocks$source_row, ]` recovers them with
their own columns and in their own order. It runs from 1 to
`params$n_blocks` and is the identity when nothing was dropped.
`params$block_sizes` is the number of points in each block, indexed by
`block_id` (zeros are empty blocks that `drop_empty_blocks = FALSE`
kept), and `params$fold_blocks` is a list with one integer vector per
fold naming the blocks packed into it. Between them the folds account
for every block exactly once, empty ones included, so a fold's territory
on the map is all of its blocks and not merely the ones that happen to
hold points. So `table(assignment$fold)` can be traced back to the
blocks it is made of, a fold can be seen to be one contiguous region or
several, and the blocks can be drawn over the data
([`plot_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_folds.md)
does so).

For the methods that work in projected space (`"block_kfold"`,
`"buffered_loo"` and `"nndm"`), `params` carries a
`params$blocks_supplied` that says whether the blocks came from `blocks`
or from a grid built here, and `params$boundary_supplied` whether a
`boundary` was given; `params$row_probe` is a small sample of row IDs
and coordinates that every `cv_*()` compares against the data it is
handed, so folds built from a different layer of the same size are
refused, never applied silently.

For the methods that work in projected space, `params` also carries a
`crs` element naming the CRS the folds were built in (an `"EPSG:code"`
string where there is one, otherwise the CRS's input definition). Every
length in `params` (`block_size`, `sac_range`, `buffer`,
`median_buffer`) is in that CRS's units, which for geographic input is a
CRS
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
chose for you, and not one you passed.

Rows whose geometry is empty or has non-finite coordinates are dropped
before folding, with a logged warning naming the count; they appear in
no fold and in no `assignment` row.

## Details

For `block_kfold`, the default grid sizing is purely geometric and
unrelated to the autocorrelation range of the data. When blocks are
smaller than the autocorrelation range, spatially correlated
observations leak across folds and CV metrics become optimistic. Use
`block_size` to set a minimum block edge length (in CRS units), or set
`auto_range = TRUE` to estimate the range from an empirical variogram
and enforce it automatically.

**Fold methods.** `"random_kfold"` ignores geography entirely and will
overstate performance on autocorrelated data. `"block_kfold"` separates
folds geographically. `"buffered_loo"` holds out one point at a time and
excludes everything within a fixed `buffer`.

`"leave_location_out"` groups by `group_var`, so all observations from a
location share a fold.

`"nndm"` implements the distance-matching principle of Milà et al.
(2022): it sizes the exclusion around each held-out point, with no
arbitrary buffer, so that the resulting training-to-test distance
distribution approaches the distribution of distances from your actual
prediction locations to the training data.

The procedure is the paper's own (as in `CAST::nndm()`), and it is
deterministic. Let \\G\_{ij}\\ be the empirical distribution of
prediction-to-nearest-training distances and \\G_j^\*\\ the distribution
of each held-out point's nearest remaining training point. Starting from
plain leave-one-out, the point with the smallest \\G_j^\*\\ at which the
realised distribution exceeds the target (\\G_j^\*(r) \> G\_{ij}(r)\\)
has its nearest training neighbour removed, and this repeats until no
such point remains, subject to two limits: a point's nearest-neighbour
distance is never pushed beyond `phi` (default: the largest prediction
distance, since a training point already further than every prediction
distance has nothing to match), and no fold's training set is stripped
below `min_train` of the data.

The realised distribution is then never *more optimistic* than the
target: \\G_j^\*(r) \le G\_{ij}(r)\\ up to the granularity of the
neighbour distances, which is the property the method exists to deliver.
An earlier version of this package drew one random radius per point from
\\G\_{ij}\\ and excluded up to the order statistic *closest* to it,
which rounds down half the time: on a two-cluster layout the realised
distribution exceeded the target by up to 0.17 (13\\ nearest training
point within 50 m against a target of 9\\ *optimistic* cross-validation.
`params$max_ecdf_excess` reports the largest remaining excess; compare
`params$target_median` with `params$realised_median` as well.

## Supplied blocks

A polygon layer passed as `blocks` is aligned to the points the way
`boundary` is: CRS-less points are aligned to blocks that carry a CRS
(reprojected if they look like lon/lat, otherwise stamped, warning
either way), and the blocks are then brought into the CRS the folds are
built in. A point inside more than one block is given the first (lowest
row) that contains it, as for a point on the shared edge of two grid
cells. When the blocks that caught such a point share area instead of an
edge, the layer overlaps and is not a partition, and this is warned
about. A point inside no block is assigned to the nearest one, by
distance to the polygon itself, and the count of such points is warned
about, unless they sit within a millionth of the extent of a block,
which is an edge that reprojection or clipping moved by a rounding
error. Blocks that hold no point are dropped when
`drop_empty_blocks = TRUE`, and `k` is lowered to the number of blocks
that hold points when that is smaller. `params$n_blocks` is the number
of blocks before empties were dropped, `params$blocks_used` the number
after (so with `drop_empty_blocks = FALSE` the two are equal, and
`sum(params$block_sizes > 0)` is how many of them hold points),
`params$grid_nx` and `params$grid_ny` are `NA`, and `params$block_scale`
is the median over blocks that hold points of the side of the square
with the block's area. That is the length compared against the
autocorrelation range for the leakage warning, since a polygon has no
single edge length.

The connection to the rest of the package is
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md):
every shape it builds can be a block design here, including Voronoi
cells around
[`get_voronoi_seeds`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)`(method = "kmeans")`
seeds, which adapt to the density of the points.

## Fold balance

Only `"block_kfold"` balances the number of points per fold, and it does
so by packing: blocks are taken largest first and each goes to the fold
with the fewest points so far, ties broken at random. This is the
longest-processing-time rule for multiway partitioning. Measured against
the optimum by enumeration (two folds, up to ten blocks, heavy-tailed
block sizes) it is optimal in 72 percent of cases and within two points
of optimal on average; a local search with 30 random restarts moved the
largest-to-smallest ratio by 0.004 on average and never brought a
packing above the 3:1 tolerance below it. An imbalance past the
tolerance is therefore in the points per block, which no assignment of
blocks to folds can even out, and this function offers no search over
packings. The remedy is the block design: smaller blocks, or blocks that
adapt to the density of the points passed through `blocks`. On clustered
layouts where the geometric grid exceeded 3:1 in 22 percent of cases
(median ratio 1.8, worst 6.3), Voronoi cells around 15 k-means seeds
never exceeded 1.3 (median 1.09). Density-adaptive blocks are smaller
where points are dense, so check `params$block_scale` against the
autocorrelation range as you would a grid.

The other methods do not balance point counts. `"random_kfold"` is
balanced by construction (fold sizes differ by at most one);
`"buffered_loo"` and `"nndm"` hold out one point per fold;
`"leave_location_out"` gives each fold the same number of *locations*
(to within one), so folds differ by as much as the locations' sizes do.

## References

Mila, C., Mateu, J., Pebesma, E. and Meyer, H. (2022). Nearest neighbour
distance matching Leave-One-Out Cross-Validation for map validation.
*Methods in Ecology and Evolution* **13**, 1304-1316.
[doi:10.1111/2041-210X.13851](https://doi.org/10.1111/2041-210X.13851)

Roberts, D. R., Bahn, V., Ciuti, S., Boyce, M. S., Elith, J.,
Guillera-Arroita, G., Hauenstein, S., Lahoz-Monfort, J. J., Schroder,
B., Thuiller, W., Warton, D. I., Wintle, B. A., Hartig, F. and Dormann,
C. F. (2017). Cross-validation strategies for data with temporal,
spatial, hierarchical, or phylogenetic structure. *Ecography* **40**,
913-929. [doi:10.1111/ecog.02881](https://doi.org/10.1111/ecog.02881)

Valavi, R., Elith, J., Lahoz-Monfort, J. J. and Guillera-Arroita, G.
(2019). blockCV: An R package for generating spatially or
environmentally separated folds for k-fold cross-validation of species
distribution models. *Methods in Ecology and Evolution* **10**, 225-232.
[doi:10.1111/2041-210X.13107](https://doi.org/10.1111/2041-210X.13107)

## See also

Other cross-validation:
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md),
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
[`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md),
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md),
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
[`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md),
[`sac_nugget()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md),
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)

## Examples

``` r
library(sf)
set.seed(1)
pts <- st_as_sf(
  data.frame(x = 5e5 + runif(30, 0, 1000), y = 5e6 + runif(30, 0, 1000)),
  coords = c("x", "y"), crs = 32632
)
folds <- make_folds(pts, k = 3, method = "block_kfold", seed = 42)
folds$assignment          # fold and block membership per row
#> # A tibble: 30 × 3
#>    row_id  fold block_id
#>     <int> <int>    <int>
#>  1      1     1        4
#>  2      2     3        8
#>  3      3     1        5
#>  4      4     2        3
#>  5      5     3        7
#>  6      6     2        9
#>  7      7     2        9
#>  8      8     3        2
#>  9      9     3        8
#> 10     10     1        4
#> # ℹ 20 more rows
lengths(folds$folds[[1]]) # train/test row-ID splits
#> train  test 
#>    20    10 
folds$params$block_sizes  # points per block; params$fold_blocks packs them
#> [1] 2 2 2 3 5 3 4 4 5

# Buffered leave-one-out: neighbours within 100 units excluded from training
loo <- make_folds(pts, k = 1, method = "buffered_loo", buffer = 100)
```
