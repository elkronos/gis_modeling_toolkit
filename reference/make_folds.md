# Create spatial cross-validation folds

Builds train/test splits by random k-fold, spatial block k-fold,
leave-location-out, buffered leave-one-out, nearest-neighbour distance
matching (NNDM) leave-one-out or its k-fold version (kNNDM).

## Usage

``` r
make_folds(
  points_sf,
  k,
  method = c("random_kfold", "block_kfold", "buffered_loo", "leave_location_out", "nndm",
    "knndm"),
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
  balance_tol = 3,
  max_test = NULL
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

  Integer; number of folds. Must be a single whole number \>= 1, and is
  required except for the two leave-one-out methods. A fraction, `NA` or
  a vector is an error, because a non-integer used to truncate silently
  and leave the last rows in no test set at all. Not every method
  honours it. `"buffered_loo"` and `"nndm"` are leave-one-out schemes
  and always return one fold per held-out point (`k = n`, or `max_test`
  when that is given) regardless of what was asked for; `"block_kfold"`
  lowers it when the grid yields fewer than `k` non-empty blocks, and
  `"leave_location_out"` lowers it when there are fewer than `k`
  distinct groups. `"random_kfold"` and `"knndm"` lower it only when the
  layer has fewer than `k` rows. `k = 1` is raised to 2 by
  `"random_kfold"`, `"block_kfold"`, `"leave_location_out"` and
  `"knndm"`, since one fold has no training set. Read the `k` element of
  the returned list, and do not assume the requested value;
  `params$k_requested` keeps the value asked for. Returning fewer folds
  than asked raises a warning of class `"spatialkit_k_reduced"`, whose
  message gives both counts, the reason and a setting that works; the
  condition carries them as `k_requested` and `k`, so
  `tryCatch(spatialkit_k_reduced = )` or
  [`withCallingHandlers()`](https://rdrr.io/r/base/conditions.html) can
  act on it without reading the text. (Before this version a reduction
  was a log line only.)

- method:

  One of `"random_kfold"`, `"block_kfold"`, `"buffered_loo"`,
  `"leave_location_out"`, `"nndm"` or `"knndm"`. See **Details** for
  what each one does and when it is appropriate.

- seed:

  Optional integer RNG seed. `"nndm"` is deterministic and reads it only
  to draw the `max_test` subset; `"buffered_loo"` likewise. `"knndm"`
  uses it for the random partition among its candidates, for ties in the
  packing, and above 4,000 points for the k-means starts, so its folds
  are reproducible only with a seed.

- block_nx, block_ny:

  Optional grid dimensions for block_kfold, each a single whole number
  \>= 1. Give both, or give one and the other is derived from the
  extent's aspect ratio so that the blocks are roughly square. Ignored
  when `block_size` or `auto_range` override them.

- block_multiplier:

  A single positive number, default 3. When neither `block_size` nor
  `block_nx`/`block_ny` is given, the automatic grid aims for
  `block_multiplier * k` blocks over the extent (aspect-preserving), so
  each fold holds out about `block_multiplier` blocks. An extent more
  than about `block_multiplier * k` times as wide as it is tall (or as
  tall as it is wide) gets a single row (or column) of that many blocks.
  With 1, every fold is one contiguous region and the score depends
  heavily on which region each fold happened to get; with many, the
  blocks shrink towards single points and the scheme drifts back towards
  random k-fold. 3 is a compromise between those two, not a published
  constant. The block size that matters for leakage is the
  autocorrelation range, which is what `block_size` and `auto_range`
  control.

- block_size:

  Optional positive numeric minimum block edge length, **in the units of
  the CRS the folds are built in**. Takes precedence over
  `block_nx`/`block_ny` and `block_multiplier`. When supplied, a whole
  number of blocks is fitted across each axis of the extent: as many as
  leave every block at least this long on that axis. On an axis at least
  as long as `block_size` the blocks built are therefore from one to
  just under two times this size. An axis shorter than `block_size`
  cannot hold a block that long: a single row or column spans it, and
  its blocks are as long as the extent on that axis, which is less than
  the size asked (a 39.8 km by 5 km extent with `block_size = 10000`
  gets three blocks 13.3 km wide and 5 km tall). `params$block_width`
  and `params$block_height` record what was built. A size that leaves
  fewer than `k` blocks holding points returns that many folds, with a
  warning (see `k`); one that leaves a single block is an error. Both
  name the largest size that still gives `k` folds.

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
  (the omnidirectional, all-pairs range; the directional ranges it also
  fits are a diagnostic only, so on a field known to be anisotropic pass
  `block_size = max(attr(r, "directional_fitted"), na.rm = TRUE)`
  yourself after checking `attr(r, "directional_status")`, see
  [`estimate_sac_range`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md))
  and used as the minimum `block_size`. Requires `response_var`. An
  explicit `block_size` takes precedence. Default `FALSE`. Sizing blocks
  from the autocorrelation range is the recommendation of Roberts et
  al. (2017) and what blockCV (Valavi et al. 2019) automates. blockCV
  takes the fitted variogram's range *parameter* as the block size,
  whereas this uses the *effective* range
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  returns (three times that parameter for an exponential fit), so its
  blocks are larger than blockCV's from the same variogram. When no
  range is identified, the geometric grid is used instead, with a
  warning that gives the reason. Range-sized blocks do not fit every
  layout. On a square extent at `k = 5` a range up to a third of the
  side gives five folds; up to half, a 2 x 2 grid and four folds,
  returned with the `"spatialkit_k_reduced"` warning; and above half the
  extent in both directions there is room for a single block, and
  `make_folds()` stops with an error naming the size, rather than return
  one fold with an empty training set. The warning and the error name
  the largest block size that still gives `k` folds on this layout and
  how far below the range it is; `method = "nndm"` is the design that
  does not depend on a block size. The message printed when the range is
  used states the grid and the block size built.

- range_frac:

  Passed through to
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  when `auto_range = TRUE`. A fitted range beyond the longest lag the
  empirical variogram was fitted over is rejected as unidentified, and
  block sizing falls back to geometry (with a warning). A range within
  that bound can still be too wide for two blocks; see `auto_range`.
  Default 1.0.

- response_var:

  Character(1) response column name. Required when `auto_range = TRUE`.
  For `"block_kfold"` a name that is not a column of `points_sf` is an
  error, whether or not `auto_range` is set: with it off the response
  still feeds the leakage warning, which a misspelt name would silently
  disable.

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
  onto. Required for `method = "nndm"` and `method = "knndm"`. The grid
  from
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

  Optional polygonal sf/sfc for block_kfold. The grid is clipped to it;
  a cell that the boundary only touches at a corner or along an edge is
  not a block. A boundary without a CRS is brought into the points' CRS:
  reprojected from EPSG:4326 when its coordinates look like lon/lat,
  otherwise stamped, with a warning either way. The same goes for
  `prediction_points` and `blocks`.

- buffer:

  For `"buffered_loo"`: a single positive number, the distance within
  which the held-out point's neighbours are excluded from its training
  set. Like `block_size` it is in the units of the CRS the folds are
  built in (`params$crs`), which for geographic (lon/lat) input is the
  metre CRS
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  chooses: 0.1 means 0.1 m, not 0.1 degrees. A buffer that excludes no
  neighbour from any fold makes the scheme plain leave-one-out, and is
  warned about.

- min_train:

  For `method = "nndm"`: the smallest fraction of the data any fold's
  training set may be reduced to by neighbour exclusion. Default `0.5`,
  as in `CAST::nndm()`. Where it binds, the distance matching stops
  short and the cross-validation stays optimistic; see **Details**. For
  `method = "knndm"`: the smallest fraction of the data a fold may leave
  for training, so a candidate whose largest fold holds more than
  `1 - min_train` of the rows is set aside (the `maxp = 0.5` of
  `CAST::knndm()` is `min_train = 0.5`). Equal folds leave `1 - 1/k` of
  the data, so a value above that is an error.

- phi:

  For `method = "nndm"`: the distance up to which the two
  nearest-neighbour distance distributions are matched, in the CRS the
  folds are built in. Matching is attempted only while a held-out
  point's nearest-neighbour distance is at most `phi`; the exclusion
  that takes it past `phi` is the last one, so a realised distance can
  exceed `phi` by up to one neighbour step, as in `CAST::nndm()`. In
  Mila et al. (2022), and in `CAST::nndm()`, \\\phi\\ is the
  autocorrelation range of the outcome: beyond it observations are
  effectively independent, so matching is unnecessary.
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  gives such a value. Default `NULL` = the largest
  prediction-to-training distance, which matches everywhere (`CAST`'s
  `phi = "max"`).

- drop_empty_blocks:

  Logical. Default TRUE. With `FALSE` the blocks that hold no point are
  kept and packed into folds too, but `k` is still lowered to the number
  of blocks that hold points, so no fold is left without test points.

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
  `params$balance_ratio` carries the ratio achieved. For `"knndm"` it is
  a limit and not a check: a candidate whose largest fold holds more
  than `balance_tol` times the points of its smallest is set aside, so
  the folds returned are always within it.

- max_test:

  For `"buffered_loo"` and `"nndm"`: optional whole number of at
  least 2. Folds are built for a random subset of that many held-out
  points (drawn under `seed`), and every row stays available for
  training, so each fold returned is the fold the full design has for
  that point. The model is then refitted `max_test` times and not `n`
  times, and the folds hold `4 * n * max_test` bytes of row indices in
  place of `4 * n^2`, which lets `"buffered_loo"` run on layers its
  storage limit otherwise refuses. It does not lift the 5,000-row limit
  of `"nndm"`, whose cost is the matching itself; for that there is
  `"knndm"`. Rows that are not held out have `NA` in `assignment$fold`.
  For `"nndm"` the matching still runs over all `n` points and `params`
  describes all of them.
  [`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
  needs every training row in a test fold and refuses a subset. Default
  `NULL`: every point is held out once. An error for the k-fold methods,
  whose number of folds is `k`.

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
`params$n_blocks` and is the identity when nothing was dropped. For a
grid `params$n_blocks` is `grid_nx * grid_ny`, the cells a `boundary`
clips away included, and `params$blocks_used` is the number of blocks
returned. `params$block_sizes` is the number of points in each block,
indexed by `block_id` (zeros are empty blocks that
`drop_empty_blocks = FALSE` kept), and `params$fold_blocks` is a list
with one integer vector per fold naming the blocks packed into it.
Between them the folds account for every block exactly once, empty ones
included, so a fold's territory on the map is all of its blocks and not
merely the ones that happen to hold points. So `table(assignment$fold)`
can be traced back to the blocks it is made of, a fold can be seen to be
one contiguous region or several, and the blocks can be drawn over the
data
([`plot_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_folds.md)
does so).

For `"block_kfold"`, `params` carries a `params$blocks_supplied` that
says whether the blocks came from `blocks` or from a grid built here,
and `params$boundary_supplied` whether a `boundary` was given.
`params$block_width` and `params$block_height` are the width and height
of the grid cells built, before any clipping to a `boundary` (`NA` for
supplied `blocks`); they are what `block_size` turned into: at least
that size on an axis long enough to hold it, and the whole extent on an
axis that is not (see `block_size`). Every method's `params` carries
`params$k_requested`, the `k` that was asked for (`NA` when none was
given), beside the `k` returned, and `params$row_probe`, a small sample
of row IDs and coordinates that every `cv_*()` compares against the data
it is handed, so folds built from a different layer of the same size are
refused, never applied silently.

For `"knndm"`, `params` carries `chosen` (`"random"` or `"clusters"`),
`n_clusters` (the number of clusters packed into the folds; `NA` for the
random partition), `clustering` (`"ward"` or `"kmeans+ward"`; `NA` when
no cluster count was tried), `candidates` (a data frame with one row per
candidate, the random partition first: `candidate`, `n_clusters`,
`largest_fold`, `balance_ratio`, `admissible` and `W`, which is `NA` for
a candidate set aside), `wasserstein`, `wasserstein_random` and
`wasserstein_loo` (the Wasserstein distance to the target of the folds
returned, of the random partition among the candidates and of plain
leave-one-out, in distance units), `fold_sizes`, `balance_ratio`,
`target_distances` and `realised_distances` (one per prediction point,
and one per row: its distance to the nearest point of another fold) with
their medians and means, and the diagnostics `"nndm"` reports under the
same names: `max_ecdf_excess`, `max_ecdf_deficit`, `loo_ecdf_deficit`,
`ecdf_deficit_tol`, `target_separated` and `target_out_of_reach`, with
`short_of_target` for the opposite case (see **Details**).

With `max_test`, `"buffered_loo"` and `"nndm"` return that many folds,
`assignment$fold` is `NA` for the rows not held out, and
`params$max_test` records the value.

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
and enforce it automatically. When a response is given and the blocks
come out smaller than the estimated range, the warning says what
range-sized blocks would yield on this layout before recommending them:
the setting and the grid it builds when that gives `k` folds, and
otherwise how many folds it gives (or that it leaves a single block),
the largest block size that still gives `k` folds, and `method = "nndm"`
as the design that needs no block size.

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

The procedure is the paper's, and it is deterministic. Let \\G\_{ij}\\
be the empirical distribution of prediction-to-nearest-training
distances and \\G_j^\*\\ the distribution of each held-out point's
nearest remaining training point. Starting from plain leave-one-out, the
point with the smallest \\G_j^\*\\ at which the realised distribution
exceeds the target (\\G_j^\*(r) \> G\_{ij}(r)\\) has its nearest
training neighbour removed, and this repeats until no such point
remains, subject to two limits: a point is matched only while its
nearest-neighbour distance is at most `phi` (default: the largest
prediction distance, since a training point already further than every
prediction distance has nothing to match), so the exclusion that takes
it past `phi` is its last and a realised distance can exceed `phi` by up
to one neighbour step; and no fold's training set is stripped below
`min_train` of the data.

It differs from `CAST::nndm()` in two details, so the folds agree
closely with CAST's but not exactly. The removal rule is strict: a
neighbour is removed whenever the realised distribution exceeds the
target, whereas CAST removes one only while the realised distribution,
less the point about to move, is still at or above the target. The rule
here therefore removes up to one point more per distance value (a few
percent more removals in all on clustered layouts), erring on the
pessimistic side. And ties in \\G_j^\*\\ are broken by the points'
coordinates, not by their row index as in CAST, so the folds do not
depend on the order of the rows.

Where neither limit binds, the realised distribution of
*nearest-neighbour distances* is then never shorter than the target:
\\G_j^\*(r) \le G\_{ij}(r)\\ up to the granularity of the neighbour
distances, which is the property the method exists to deliver. It is a
statement about the distance to the one nearest training point. A
learner that draws on several neighbours (a forest, a
k-nearest-neighbour mean, a kriging neighbourhood) also depends on how
far the second and later neighbours are, which NNDM does not match, so
on clustered samples its cross-validated error can still read low with
the nearest-neighbour distances matched exactly. Beyond `phi` no
matching is attempted, by design. `min_train` is different: when the
prediction locations lie further from the samples than a fold can be
made to hold out (clustered samples and a prediction domain well beyond
them, the layout NNDM is meant for), it stops the matching early and the
realised distances stay too short. `params$n_at_min_train` counts the
folds that were held at the floor while still closer than the target
allows, and `make_folds()` warns when that leaves more than one point's
worth of excess at or below `phi`. Lower `min_train` to match further,
or read the cross-validated score as an upper bound on performance at
the prediction locations.

The match has a second side. `params$max_ecdf_deficit` is the largest
amount by which the target distribution runs above the realised one, the
mirror of `params$max_ecdf_excess`: the share of prediction points
within some distance of the data that the folds do not reproduce. Part
of it is the granularity of the matching. Excluding a neighbour moves a
held-out point out to its next remaining neighbour, and on clustered
samples that can be the next cluster, well past the prediction
distances.

The other part no setting reaches. Excluding neighbours can only
lengthen a held-out distance, so when the prediction points sit *nearer*
to the data than the samples sit to each other (prediction points among
the samples, or a regular sample predicted between its nodes) the folds
returned are plain leave-one-out and are still farther from their
training data than the prediction points are. The cross-validated error
is then pessimistic for those locations: an upper bound on the error
there. `params$loo_ecdf_deficit` is the largest amount by which the
target distribution runs above that of plain leave-one-out. When it
exceeds `params$ecdf_deficit_tol`, three times \\\sqrt{1/n + 1/m}\\ for
\\n\\ folds and \\m\\ prediction points (about three times the chance
difference between two samples of one distribution), `make_folds()`
raises a warning of class `"spatialkit_nndm_out_of_reach"` and sets
`params$target_out_of_reach` to `TRUE`.

On small samples no deficit can exceed that tolerance. A deficit is a
difference of two shares, so it is at most 1, and the tolerance is 1 or
more whenever \\nm/(n+m) \le 9\\: for nine folds or fewer whatever the
number of prediction points, and for up to eighteen folds with as many
prediction points. There the same warning is raised in one case only,
complete separation: every prediction point nearer to the data than any
sample is to its nearest other sample (a deficit of exactly 1), provided
the two sets are large enough for that ordering to mean something, which
is taken as `choose(n + m, n) >= 1000` ways of interleaving them (17
prediction points for 3 folds, 7 for 6 folds, 4 for 12 folds).
`params$target_separated` records the ordering whether or not it was
enough to warn, and the warning's condition carries it as `separated`.
With very few samples it can hold for a grid over the samples' own
extent, because a few points hug the box they span; the folds are then
farther than the grid all the same. A partial shortfall on such a sample
(most prediction points nearer than leave-one-out reaches, not all) is
not warned about: the tolerance is 1 or more at these sizes, which no
difference of shares exceeds, so only complete separation can be
reported. Compare `params$target_median` with `params$realised_median`.

`params$wasserstein` is the Wasserstein distance between the target and
realised distributions, in distance units, and `params$wasserstein_loo`
the same for plain leave-one-out: equal when nothing was excluded,
smaller when the matching helped.

**kNNDM.** `"knndm"` is the k-fold version of that idea (Linnenbrink et
al. 2024): `k` folds, and so `k` model fits, whose held-out distances
are as close to the same target as a set of candidates allows. It has no
row limit, and is the method to use above the 5,000 rows `"nndm"` takes
or when a model is too costly to refit once per point. It is not an
equal of NNDM below that: a fold design has far fewer degrees of freedom
than one exclusion set per point, so the match is coarser.

The candidates are a random k-fold partition (the one `"random_kfold"`
builds with the same `seed`) and, for each of up to 50 cluster counts
spaced evenly on the log scale from `k` to one fewer than the number of
distinct locations, the samples clustered on their coordinates and the
clusters merged into `k` folds. Few large clusters give folds that are
whole regions, far from their training data; many small ones approach
the random partition. Up to 4,000 points the clustering is Ward's
(`stats::hclust(method = "ward.D2")`), cut at each count. Above that,
where the n by n distance matrix it needs is no longer affordable,
k-means first reduces the samples to 4,000 groups (Lloyd's algorithm
from random sample locations, at most ten passes) and Ward's method is
continued from those groups; a count above 4,000 is a k-means run of its
own. `params$clustering` says which (`"ward"` or `"kmeans+ward"`).
Clusters are merged by the packing `"block_kfold"` uses for its blocks,
largest first and each to the fold with the fewest points so far (see
**Fold balance**). `CAST::knndm()` merges clusters by a rule of its own,
which this was *not* checked against, so the folds are not expected to
equal CAST's.

Each candidate is scored by the Wasserstein distance \\W\\ between its
held-out distances (each point's distance to the nearest point of
another fold) and the target distances, and the smallest \\W\\ is
returned. A candidate is set aside without being scored when its largest
fold leaves less than `min_train` of the rows for training or holds more
than `balance_tol` times the points of its smallest; the random
partition is always kept. So `params$wasserstein` is never above
`params$wasserstein_random`, the \\W\\ of that random partition, and
`params$candidates` lists every candidate with its cluster count, fold
sizes and \\W\\. `params$chosen` is `"random"` when no clustering came
closer, which is the right answer for samples spread evenly over the
prediction area.

Read \\W\\ against the mean target distance. On 400 points in 12
clusters predicted onto a grid over their extent (12 draws, 5 and 10
folds) the folds returned had \\W\\ at 0.32 of the mean target distance,
against 0.70 to 0.72 for a random partition and 0.01 for `"nndm"`. A
5-nearest-neighbour learner scored 9 to 11 percent high on them, where
random folds scored 23 to 25 percent low and `"nndm"` 6 percent low. On
2,000 points in 40 clusters of unequal size the ratio was 0.07 to 0.10,
against 0.56 to 0.59 for a random partition and 0.003 to 0.006 for
`"nndm"`; on 30,000 points of that layout it was 0.11 to 0.15, and on
100,000 points 0.27. So the method removes most of the mismatch of a
random partition and leaves more of it than NNDM does, in either
direction: \\W\\ has no sign. On 400 evenly spread points the random
partition or one close to it came back and the two scored alike (3.3 and
3.5 percent high). With few folds a random partition is itself farther
than the target, since each fold trains on `1 - 1/k` of the points: at
`k = 2` on evenly spread points the mean held-out distance was 9 to 53
percent above the target mean (60 draws of 100 to 1,000 points), and no
candidate is nearer.

Construction took 2 s at 8,000 points and 10 s at 30,000 with ten folds
(37 s at 100,000) with FNN installed, and three to four times that
without it; it grows with the number of folds.

Two cases are warned about, each with a condition class. When the
prediction points sit nearer to the data than a sample is to its nearest
other sample, no fold design reaches them, since a held-out point is
never nearer to its training data than that; the test is the one
`"nndm"` applies, the warning has the same class,
`"spatialkit_nndm_out_of_reach"`, and what comes back is the random
partition or a candidate as fine as one. In the other direction, the
closest candidate can still be short of the target: prediction points
well beyond the samples lie farther from the data than `k` folds of
those samples can be held apart. When the mean held-out distance of the
folds returned is below the mean target distance by more than a quarter
of it and by at least half of \\W\\ (the test
[`fold_separation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fold_separation.md)
words as `"closer"`), `make_folds()` raises a warning of class
`"spatialkit_knndm_short_of_target"` and sets `params$short_of_target`:
the cross-validated error then tends to read low for those locations.
Fewer folds hold out larger regions. Short of either warning the folds
can still differ from the target at some distances;
`params$max_ecdf_excess` and `params$max_ecdf_deficit` give the largest
difference each way, and
[`fold_separation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fold_separation.md)
reads the folds against the same prediction points fold by fold.

An earlier version of this package drew one random radius per point from
\\G\_{ij}\\ and excluded up to the order statistic *closest* to it,
which rounds down half the time: on a two-cluster layout the realised
distribution exceeded the target by up to 0.17 (13% of folds had a
nearest training point within 50 m against a target of 9%), an
*optimistic* cross-validation. `params$max_ecdf_excess` reports the
largest remaining excess; compare `params$target_median` with
`params$realised_median` as well.

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

`"block_kfold"` and `"knndm"` balance the number of points per fold, and
they do so by packing: blocks (for `"knndm"`, clusters of the samples)
are taken largest first and each goes to the fold with the fewest points
so far, ties broken at random. This is the longest-processing-time rule
for multiway partitioning. Measured against the optimum by enumeration
(two folds, up to ten blocks, heavy-tailed block sizes) it is optimal in
72 percent of cases and within two points of optimal on average; a local
search with 30 random restarts moved the largest-to-smallest ratio by
0.004 on average and never brought a packing above the 3:1 tolerance
below it. An imbalance past the tolerance is therefore in the points per
block, which no assignment of blocks to folds can even out, and this
function offers no search over packings. The remedy is the block design:
smaller blocks, or blocks that adapt to the density of the points passed
through `blocks`. On clustered layouts where the geometric grid exceeded
3:1 in 22 percent of cases (median ratio 1.8, worst 6.3), Voronoi cells
around 15 k-means seeds never exceeded 1.3 (median 1.09).
Density-adaptive blocks are smaller where points are dense, so check
`params$block_scale` against the autocorrelation range as you would a
grid.

`"knndm"` chooses its clusters, so for it `balance_tol` is a limit: a
clustering whose packing exceeds it is not a candidate, and the folds
returned are within it (`params$fold_sizes`, `params$balance_ratio`).

The other methods do not balance point counts. `"random_kfold"` is
balanced by construction (fold sizes differ by at most one);
`"buffered_loo"` and `"nndm"` hold out one point per fold;
`"leave_location_out"` gives each fold the same number of *locations*
(to within one), so folds differ by as much as the locations' sizes do.

## References

Linnenbrink, J., Mila, C., Ludwig, M. and Meyer, H. (2024). kNNDM CV:
k-fold nearest-neighbour distance matching cross-validation for map
accuracy estimation. *Geoscientific Model Development* **17**,
5897-5912.
[doi:10.5194/gmd-17-5897-2024](https://doi.org/10.5194/gmd-17-5897-2024)

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
[`fold_separation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fold_separation.md),
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

# Three folds matched to the distances from a prediction grid to the data
grid <- st_make_grid(pts, n = c(10, 10), what = "centers")
km <- make_folds(pts, k = 3, method = "knndm", prediction_points = grid,
                 seed = 42)
km$params[c("chosen", "n_clusters", "wasserstein", "wasserstein_random")]
#> $chosen
#> [1] "random"
#> 
#> $n_clusters
#> [1] NA
#> 
#> $wasserstein
#> [1] 20.95776
#> 
#> $wasserstein_random
#> [1] 20.95776
#> 
```
