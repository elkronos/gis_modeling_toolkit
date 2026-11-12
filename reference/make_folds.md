# Create spatial cross-validation folds

Builds train and test splits of a layer for cross-validation, by one of
six methods: random k-fold, spatial block k-fold, buffered
leave-one-out, leave-location-out, nearest-neighbour distance matching
(NNDM) leave-one-out, or its k-fold version (kNNDM). It returns the
splits (`folds`), one row per observation naming the fold that holds it
out (`assignment`), and the settings and diagnostics of the build
(`params`). Pass the result as `folds` to
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
or another `cv_*()` function, map it with
[`plot_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_folds.md),
and measure how far the held-out points sit from their training data
with
[`fold_separation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fold_separation.md).

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

  An sf layer of the observations. Z and M are dropped before folding,
  since distances are measured in x and y. A layer without a CRS is
  aligned to a `boundary`, `blocks` or `prediction_points` that has one,
  by the rule in
  [spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md).
  Rows without a location are dropped (see Value).

- k:

  The number of folds, a single whole number of at least 1; required
  except for `"buffered_loo"` and `"nndm"`, which build one fold per
  held-out point (`k = n`, or `max_test` of them when it is given)
  whatever `k` is. A fraction, `NA` or a vector is an error.
  `"block_kfold"` lowers `k` when fewer than `k` blocks hold points,
  `"leave_location_out"` when there are fewer than `k` groups, and
  `"random_kfold"` and `"knndm"` when the layer has fewer than `k` rows;
  those four raise `k = 1` to 2, since one fold has no training set.
  Read the `k` element of the result rather than assume the value asked
  for, which `params$k_requested` keeps. Fewer folds than asked raise a
  warning of class `"spatialkit_k_reduced"`, whose message gives both
  counts, the reason and, where there is one, a setting that works; the
  condition carries them as `k_requested` and `k`, so
  `tryCatch(spatialkit_k_reduced = )` or
  [`withCallingHandlers()`](https://rdrr.io/r/base/conditions.html) can
  act on it without reading the text.

- method:

  One of the six methods in the section "Fold methods": `"random_kfold"`
  (default), `"block_kfold"`, `"buffered_loo"`, `"leave_location_out"`,
  `"nndm"` or `"knndm"`.

- seed:

  Optional integer RNG seed; `NULL` (default) leaves the draws to the
  session's RNG state. It draws the partition of `"random_kfold"`, the
  assignment of groups to folds of `"leave_location_out"`, the packing
  ties of `"block_kfold"` and `"knndm"`, the random candidate of
  `"knndm"` and, above 4,000 points, its k-means starts, and the
  `max_test` subset of `"buffered_loo"` and `"nndm"`, which are
  otherwise deterministic. With `auto_range = TRUE` it is also the seed
  of
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
  which keeps its own default of 123 when `seed` is `NULL`.

- block_nx, block_ny:

  Optional grid dimensions for `"block_kfold"`, each a single whole
  number of at least 1. Give both, or give one and the other is derived
  from the extent's aspect ratio so that the blocks are about square.
  Ignored when `block_size` or `auto_range` set the grid.

- block_multiplier:

  A single positive number, default 3. When neither `block_size` nor
  `block_nx`/`block_ny` is given, the grid aims for
  `block_multiplier * k` blocks over the extent, keeping its aspect, so
  each fold holds out about `block_multiplier` blocks; an extent more
  than about `block_multiplier * k` times as wide as it is tall (or as
  tall as it is wide) gets a single row (or column) of that many blocks.
  With 1, every fold is one contiguous area and the score depends on
  which area each fold gets; with many, the blocks shrink towards single
  points and the folds approach random k-fold. 3 is a compromise between
  the two, not a published constant.

- block_size:

  Optional positive number, the minimum block edge length, in the units
  of the CRS the folds are built in (`params$crs`): your own CRS's units
  for projected input, and for lon/lat input those of the metre CRS
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  chooses, which varies with the data (project the data yourself to fix
  the units in advance). It takes precedence over `block_nx`/`block_ny`
  and `block_multiplier`. The blocks built are at least this long on an
  axis that can hold one, and up to just under twice as long;
  `params$block_width` and `params$block_height` record them. See the
  section "Block size and the autocorrelation range" for how the grid is
  laid and for a size that leaves too few blocks.

- auto_range:

  Logical, default `FALSE`. `TRUE` estimates the autocorrelation range
  with
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  and uses it as the minimum `block_size`; it needs `response_var`, and
  an explicit `block_size` takes precedence. With supplied `blocks` the
  range is only compared with them. See the section "Block size and the
  autocorrelation range".

- range_frac:

  Passed to
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  when `auto_range = TRUE`: a fitted range beyond the longest lag the
  empirical variogram was fitted over is refused as unidentified, and
  the grid falls back to geometry, with a warning. A range within that
  bound can still be too wide for two blocks. Default 1.

- response_var:

  Name of the response column. Required when `auto_range = TRUE`. For
  `"block_kfold"` a name that is not a column of `points_sf` is an error
  whether or not `auto_range` is set, since the response also feeds the
  warning on blocks smaller than the range, which a misspelt name would
  switch off.

- group_var:

  Name of a column of `points_sf` that identifies the location each
  observation belongs to. Required for `method = "leave_location_out"`,
  which keeps every observation from a location in one fold. Repeated
  measurements at a site split across folds score the model partly on
  sites it has already seen.

- prediction_points:

  Optional sf layer of the locations the model will predict at, such as
  the grid from
  [`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md);
  required for `"nndm"` and `"knndm"`. A layer of cells or polygons is
  reduced to representative points first, so the target distances are
  point to point: a point-to-polygon distance is zero for every cell
  that holds a training point, which would pull the target towards zero
  and the folds towards plain leave-one-out. A layer that mixes
  dimensions or carries M is read in x and y, before and after it is
  reduced to points
  ([spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)).
  A prediction point without a location, as given or once brought into
  the CRS of the folds, takes no part in the target and is not counted
  in `params$n_prediction_points`.

- predictor_vars:

  Optional character vector of predictor column names, passed to
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  so that the range is that of the residuals on them.

- boundary:

  Optional polygon layer (sf or sfc) that the block grid of
  `"block_kfold"` is clipped to. A cell the boundary only touches at a
  corner or along an edge is not a block, and a point on the boundary,
  to the margin in the section "Numerical tolerances", is inside it. See
  [spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)
  for how a boundary without a CRS is read, and for the rows of a
  boundary that are left out or refused.

- buffer:

  For `"buffered_loo"`: a single positive number, the distance within
  which the held-out point's neighbours are excluded from its training
  set. Like `block_size` it is in the units of the CRS the folds are
  built in (`params$crs`), which for lon/lat input is the metre CRS
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  chooses: 0.1 means 0.1 m, not 0.1 degrees. A neighbour at the buffer
  distance is excluded. A buffer that excludes no neighbour from any
  fold makes the scheme plain leave-one-out, with a warning.

- min_train:

  For `"nndm"`: the smallest fraction of the data any fold's training
  set may be reduced to by neighbour exclusion; default 0.5, as in
  `CAST::nndm()`. Where it binds, the distance matching stops short and
  the held-out points stay closer to their training data than the
  prediction points are (see the section "NNDM"). For `"knndm"`: the
  smallest fraction of the data a fold may leave for training, so a
  candidate whose largest fold holds more than `1 - min_train` of the
  rows is set aside (the `maxp = 0.5` of `CAST::knndm()` is
  `min_train = 0.5`). Equal folds leave `1 - 1/k` of the data, so a
  value above that is an error.

- phi:

  For `"nndm"`: the distance up to which the two nearest-neighbour
  distance distributions are matched, in the CRS the folds are built in.
  Matching is attempted only while a held-out point's nearest-neighbour
  distance is at most `phi`; the exclusion that takes it past `phi` is
  the last one, so a realised distance can exceed `phi` by up to one
  neighbour step, as in `CAST::nndm()`. In Milà et al. (2022) and in
  `CAST::nndm()`, \\\phi\\ is the autocorrelation range of the outcome:
  beyond it observations are close to independent, so matching adds
  nothing.
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  gives such a value. Default `NULL`: the largest prediction-to-training
  distance, which matches everywhere (`phi = "max"` in CAST).

- drop_empty_blocks:

  Logical, default `TRUE`. With `FALSE` the blocks that hold no point
  are kept and packed into folds too, but `k` is still lowered to the
  number of blocks that hold points, so no fold is left without test
  points.

- blocks:

  Optional polygon layer (sf or sfc, POLYGON or MULTIPOLYGON, at least
  two features) to use as the blocks of `"block_kfold"` in place of the
  grid: the `$cells` of a
  [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
  result, hexagons, watersheds, administrative units, or the `$blocks`
  of `blockCV::cv_spatial()`. Each point takes the block that contains
  it, and the blocks are then assigned to folds as grid cells are.
  `block_size`, `block_nx`/`block_ny`, `block_multiplier` and `boundary`
  have nothing to act on and are ignored (logged); `auto_range` only
  compares the estimated range with the blocks. See the section
  "Supplied blocks". Blocks that mix dimensions or carry M are read in x
  and y, and `params$blocks` holds them in XY. An error for any other
  `method`.

- balance_tol:

  A number of at least 1, default 3. For `"block_kfold"`: the ratio of
  the largest fold's point count to the smallest's above which the folds
  are reported as imbalanced, with a warning (an R condition, also
  logged) that names both counts; `Inf` turns the check off. It is the
  tolerance of a check, not a target the packing aims for (see the
  section "Fold balance"), and `params$balance_ratio` carries the ratio
  achieved. For `"knndm"` it is a limit: a candidate whose largest fold
  holds more than `balance_tol` times the points of its smallest is set
  aside, so the folds returned are within it.

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
  `"knndm"`. For `"nndm"` the matching still runs over all `n` points
  and `params` describes all of them.
  [`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
  needs every training row in a test fold and refuses a subset. Default
  `NULL`: every point is held out once. An error for the k-fold methods,
  whose number of folds is `k`.

## Value

A list with:

- `method`:

  The method used.

- `k`:

  The number of folds built, which can differ from the `k` asked for
  (see `k`); `length(folds)` always equals it.

- `folds`:

  A list with one element per fold, each a list of `train` and `test`
  row IDs: the values of the `..row_id` column when the input has one,
  otherwise row positions.

- `assignment`:

  A tibble with one row per observation kept: `row_id` (the same IDs)
  and `fold`, and for `"block_kfold"` a third column, `block_id`, the
  row of `params$blocks` the point fell in. With `max_test`, `fold` is
  `NA` for the rows not held out.

- `params`:

  The settings and diagnostics of the build, below.

Every method's `params` carries `seed`, `k_requested` (the `k` asked
for, `NA` when none was given) and `row_probe`, a small sample of row
IDs and coordinates that every `cv_*()` function compares with the data
it is handed, so folds built from a different layer of the same size are
refused rather than applied. The methods that work in projected space
add `crs`, the CRS the folds were built in (an `"EPSG:code"` string
where there is one, otherwise the CRS's input definition). Every length
in `params` (`block_size`, `sac_range`, `buffer`, `median_buffer`) is in
that CRS's units, which for lon/lat input is a CRS
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
chose, not one you passed.

For `"block_kfold"` the block design comes with the folds:

- `blocks`: an sf layer of the block polygons in the CRS the folds were
  built in. Its rows are not the rows of the layer the blocks came from:
  `drop_empty_blocks = TRUE` (the default) removes the blocks that hold
  no point and renumbers the rest, so with supplied `blocks` a
  nine-polygon layer of which three are empty comes back as six rows
  numbered 1 to 6. `blocks$source_row` is the row each one came from (a
  grid cell's index in the full `grid_nx` by `grid_ny` grid, or the row
  of the `blocks` argument), so `blocks[params$blocks$source_row, ]`
  recovers them with their own columns and in their own order. It runs
  from 1 to `n_blocks` and is the identity when nothing was dropped.

- `n_blocks` and `blocks_used`: for a grid, `n_blocks` is
  `grid_nx * grid_ny`, the cells a `boundary` clips away included, and
  `blocks_used` the number of blocks returned.

- `block_sizes`: the number of points in each block, indexed by
  `block_id` (zeros are empty blocks that `drop_empty_blocks = FALSE`
  kept); `fold_blocks`: a list with one integer vector per fold naming
  the blocks packed into it. Between them the folds account for every
  block once, empty ones included, so a fold's territory on the map is
  all of its blocks, and `table(assignment$fold)` can be traced back to
  the blocks it is made of.
  [`plot_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_folds.md)
  draws the blocks over the data.

- `block_width` and `block_height`: the width and height of the grid
  cells built, before any clipping to a `boundary` (`NA` for supplied
  `blocks`). They are what `block_size` turned into: at least that size
  on an axis long enough to hold it, and the whole extent on an axis
  that is not.

- `blocks_supplied` and `boundary_supplied`: whether the blocks came
  from `blocks`, and whether a `boundary` was given.

- `block_scale`, `balance_ratio`, `balance_tol`, `block_size`,
  `block_multiplier`, `sac_range` and `auto_range`: the length compared
  with the range (see the section "Supplied blocks"), the balance
  achieved and allowed, and the sizing settings.

For `"nndm"`, `params` carries the diagnostics named in the section
"NNDM". For `"knndm"` it carries `chosen` (`"random"` or `"clusters"`),
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
`short_of_target` for the opposite case (see the section "kNNDM"). With
`max_test`, `params$max_test` records the value.

Rows without a location (an empty point, or a coordinate that is not
finite or lies beyond 1e150) are dropped before folding, with a logged
warning naming the count; they appear in no fold and in no `assignment`
row. The bound is read on the coordinates as given, and again after the
projection for the methods that measure distance.

## Fold methods

|  |  |  |  |
|----|----|----|----|
| **Method** | **Held out in each fold** | **Needs** | **Suits a map that** |
| `"random_kfold"` | a random 1/k of the rows | `k` | predicts among the samples, over an area they cover evenly |
| `"block_kfold"` | whole blocks of a grid, or of `blocks` | `k`; `block_size` or `auto_range` to size the grid | predicts into unsampled gaps about as wide as the blocks |
| `"buffered_loo"` | one point, with its neighbours within `buffer` left out of training | `buffer` | predicts at a known distance from the data |
| `"leave_location_out"` | every row of one or more locations | `group_var` | predicts at sites with no measurements, from repeated measurements at others |
| `"nndm"` | one point, with neighbours left out to match the prediction distances | `prediction_points`; at most 5,000 rows | predicts at known locations, with a model cheap enough to refit once per point |
| `"knndm"` | `k` folds matched to the prediction distances | `prediction_points` | predicts at known locations, from a large layer or with a costly model |

A cross-validated error is the error at the distances the folds put
between the held-out points and their training data. It reads low for
predictions made farther from the data than the held-out points sit, and
high for predictions made closer, so which folds to use depends on where
the model will predict.
[`fold_separation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fold_separation.md)
with `prediction_points` compares the two distances, and
`method = "nndm"` and `"knndm"` build folds matched to them.

`"random_kfold"` ignores geography: a held-out point keeps most of its
neighbours in the training data, so the score is the error of predicting
among the samples. That is the error of a map over an area the samples
cover evenly, and it reads low for predictions made farther from the
data, such as a map of a whole region from clustered samples.
`"block_kfold"` holds out whole blocks, so the held-out points sit
farther from their training data the wider the blocks are, and the score
is the error of predicting into unsampled gaps of that size; for a map
of an area the samples cover, it reads high. `"buffered_loo"` holds out
one point at a time and excludes everything within a fixed `buffer` from
its training set. `"leave_location_out"` groups by `group_var`, so all
observations from a location share a fold.

## Block size and the autocorrelation range

With `block_size`, each axis of the extent gets as many blocks as leave
every block at least that long, so on an axis at least as long as
`block_size` the blocks are from one to just under two times this size.
An axis shorter than `block_size` is spanned by one row or column of
blocks as long as the extent there: a 39.8 km by 5 km extent with
`block_size = 10000` gets three blocks 13.3 km wide and 5 km tall. A
size that leaves fewer than `k` blocks holding points returns that many
folds, with a warning (see `k`); one that leaves a single block is an
error. Both name the largest size that still gives `k` folds. A size in
the wrong unit asks for an enormous grid, so a request above 1,000,000
blocks is refused with an error naming the grid dimensions, the extent
and the CRS's units.

Without `block_size` or `auto_range`, the grid is sized from the extent
and `block_multiplier` (or `block_nx` and `block_ny`) alone, with no
regard to the autocorrelation range of the data. The wider the blocks,
the farther the held-out points sit from their training data. When
blocks are smaller than the autocorrelation range, each held-out point
still has training data within the range. The score is then the error of
predicting that close to the data. It reads low for a map whose
locations lie farther from the samples, and is the error of a map whose
locations lie as close.

`block_size` sets a minimum block edge length in CRS units, and
`auto_range = TRUE` estimates the range with
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
and uses it as that minimum. When a response is given and the blocks are
smaller than the estimated range, `make_folds()` warns and gives the
block size of the range. When blocks of that size leave fewer than `k`
folds on this layout, the warning also names the largest size that still
gives `k` folds and a method that needs no block size (`"nndm"`, or
`"knndm"` above 5,000 rows).

Sizing blocks from the autocorrelation range is the recommendation of
Roberts et al. (2017) and what blockCV (Valavi et al. 2019) automates.
blockCV takes the fitted variogram's range parameter as the block size;
`auto_range` takes the effective range that
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
returns, three times that parameter for an exponential fit, so its
blocks are larger than blockCV's from the same variogram. The estimate
is the omnidirectional range. On a field known to be anisotropic, pass
`block_size = max(attr(r, "directional_fitted"), na.rm = TRUE)` yourself
after checking `attr(r, "directional_status")` (see the section
"Directional ranges and anisotropy" in
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)).

Range-sized blocks do not fit every layout. On a square extent at
`k = 5`, a range up to a third of the side gives five folds; up to half,
a 2 by 2 grid and four folds, returned with the `"spatialkit_k_reduced"`
warning. Above half the extent in both directions there is room for a
single block, and `make_folds()` stops with an error naming the size
rather than return one fold with an empty training set. The warning and
the error name the largest block size that still gives `k` folds on this
layout and how far below the range it is. When no range is identified,
the geometric grid is used instead, with a warning that gives the
reason. When the range is used, a message states the grid and the block
size built.

## NNDM

`"nndm"` implements the distance-matching principle of Milà et al.
(2022): it sizes the exclusion around each held-out point, with no fixed
buffer, so that the distribution of distances from the held-out points
to their training data approaches the distribution of distances from the
prediction locations to the training data.

The procedure is the paper's, and it is deterministic. Let \\G\_{ij}\\
be the empirical distribution of prediction-to-nearest-training
distances and \\G_j^\*\\ the distribution of each held-out point's
distance to its nearest remaining training point. Starting from plain
leave-one-out, the point with the smallest \\G_j^\*\\ at which the
realised distribution exceeds the target (\\G_j^\*(r) \> G\_{ij}(r)\\)
has its nearest training neighbour removed, and this repeats until no
such point remains, within two limits. A point is matched only while its
nearest-neighbour distance is at most `phi` (default: the largest
prediction distance, since a training point already farther than every
prediction distance has nothing to match), so the exclusion that takes
it past `phi` is its last and a realised distance can exceed `phi` by up
to one neighbour step. And no fold's training set is stripped below
`min_train` of the data.

It differs from `CAST::nndm()` in two details, so the folds agree
closely with CAST's but not exactly. The removal rule is strict: a
neighbour is removed whenever the realised distribution exceeds the
target, where CAST removes one only while the realised distribution,
less the point about to move, is still at or above the target. The rule
here therefore removes up to one point more per distance value, erring
towards held-out distances longer than the target. And ties in
\\G_j^\*\\ are broken by the points' coordinates, not by their row index
as in CAST, so the folds do not depend on the order of the rows (see the
section "Numerical tolerances").

Where neither limit binds, the realised distribution of
nearest-neighbour distances is never shorter than the target:
\\G_j^\*(r) \le G\_{ij}(r)\\ up to the granularity of the neighbour
distances, which is the property the method exists to deliver. It is a
statement about the distance to the one nearest training point. A
learner that draws on several neighbours (a forest, a
k-nearest-neighbour mean, a kriging neighbourhood) also depends on how
far the second and later neighbours are, which NNDM does not match, so
on clustered samples its cross-validated error can still read low for
the prediction points with the nearest-neighbour distances matched
exactly. Beyond `phi` no matching is attempted. `min_train` is
different: when the prediction locations lie farther from the samples
than a fold can be made to hold out (clustered samples and a prediction
domain well beyond them, the layout NNDM is meant for), it stops the
matching early and the realised distances stay too short.
`params$n_at_min_train` counts the folds held at the floor while still
closer than the target allows, and `make_folds()` warns when that leaves
more than one point's worth of excess at or below `phi`. Lower
`min_train` to match further, or read the cross-validated score as an
upper bound on performance at the prediction locations.
`params$max_ecdf_excess` reports the largest remaining excess; compare
`params$target_median` with `params$realised_median` as well.

The match has a second side. `params$max_ecdf_deficit` is the largest
amount by which the target distribution runs above the realised one, the
mirror of `params$max_ecdf_excess`: the share of prediction points
within some distance of the data that the folds do not reproduce. Part
of it is the granularity of the matching. Excluding a neighbour moves a
held-out point out to its next remaining neighbour, and on clustered
samples that can be the next cluster, well past the prediction
distances.

The other part no setting reaches. Excluding neighbours can only
lengthen a held-out distance, so when the prediction points sit nearer
to the data than the samples sit to each other (prediction points among
the samples, or a regular sample predicted between its nodes) the folds
returned are plain leave-one-out and are still farther from their
training data than the prediction points are. The cross-validated error
is then pessimistic for those prediction points: an upper bound on the
error there. `params$loo_ecdf_deficit` is the largest amount by which
the target distribution runs above that of plain leave-one-out. When it
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
sample is to its nearest other sample, by more than one part in ten
million (a deficit of 1), provided the two sets are large enough for
that ordering to mean something, taken as `choose(n + m, n) >= 1000`
ways of interleaving them (17 prediction points for 3 folds, 7 for 6
folds, 4 for 12 folds). `params$target_separated` records the ordering
whether or not it was enough to warn, and the warning's condition
carries it as `separated`. With very few samples it can hold for a grid
over the samples' own extent, because a few points hug the box they
span; the folds are then farther than the grid all the same. A partial
shortfall on such a sample (most prediction points nearer than
leave-one-out reaches, not all) is not warned about, since no difference
of shares exceeds a tolerance of 1 or more. Compare
`params$target_median` with `params$realised_median`.

`params$wasserstein` is the Wasserstein distance between the target and
realised distributions, in distance units, and `params$wasserstein_loo`
the same for plain leave-one-out: equal when nothing was excluded,
smaller when the matching helped.

## kNNDM

`"knndm"` is the k-fold version of NNDM (Linnenbrink et al. 2024): `k`
folds, and so `k` model fits, whose held-out distances are as close to
the same target as a set of candidates allows. It has no row limit, and
is the method to use above the 5,000 rows `"nndm"` takes or when a model
is too costly to refit once per point. Below that it matches less
closely than NNDM, since a fold design has far fewer degrees of freedom
than one exclusion set per point.

The candidates are a random k-fold partition (the one `"random_kfold"`
builds with the same `seed`) and, for each of up to 50 cluster counts
spaced evenly on the log scale from `k` to one fewer than the number of
distinct locations, the samples clustered on their coordinates and the
clusters merged into `k` folds. Few large clusters give folds that are
whole areas, far from their training data; many small ones approach the
random partition. Up to 4,000 points the clustering is Ward's
(`stats::hclust(method = "ward.D2")`), cut at each count. Above that,
where the n by n distance matrix it needs is too large, k-means first
reduces the samples to 4,000 groups (Lloyd's algorithm from random
sample locations, at most ten passes) and Ward's method is continued
from those groups; a count above 4,000 is a k-means run of its own.
`params$clustering` says which (`"ward"` or `"kmeans+ward"`). Clusters
are merged by the packing `"block_kfold"` uses for its blocks, largest
first and each to the fold with the fewest points so far (see the
section "Fold balance"). `CAST::knndm()` merges clusters by a rule of
its own, which this was not checked against, so the folds are not
expected to equal CAST's.

Each candidate is scored by the Wasserstein distance \\W\\ between its
held-out distances (each point's distance to the nearest point of
another fold) and the target distances, in the unit of the layer, and
the smallest \\W\\ is returned; of the candidates within one part in ten
million of the smallest, the earliest in `params$candidates`. A
candidate is set aside without being scored when its largest fold leaves
less than `min_train` of the rows for training or holds more than
`balance_tol` times the points of its smallest; the random partition is
always kept. So `params$wasserstein` is never above
`params$wasserstein_random`, the \\W\\ of that random partition, and
`params$candidates` lists every candidate with its cluster count, fold
sizes and \\W\\. `params$chosen` is `"random"` when no clustering came
closer, as for samples spread evenly over the prediction area. Read
\\W\\ against the mean target distance; it has no sign, so it measures
the mismatch in either direction. With few folds a random partition can
itself be farther than the target, since each fold trains on `1 - 1/k`
of the points. Construction time grows with the number of points and of
folds.

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
the cross-validated error then tends to read low for the prediction
points. Fewer folds hold out larger areas. Short of either warning the
folds can still differ from the target at some distances;
`params$max_ecdf_excess` and `params$max_ecdf_deficit` give the largest
difference each way, and
[`fold_separation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fold_separation.md)
reads the folds against the same prediction points fold by fold.

## Supplied blocks

A polygon layer passed as `blocks` is aligned to the points the way
`boundary` is, and then brought into the CRS the folds are built in (see
[spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)).
A point inside more than one block is given the first (lowest row) that
contains it, as for a point on the shared edge of two grid cells, and a
point within the margin of a block's edge (see the section "Numerical
tolerances") counts as inside it. When the blocks that caught such a
point share area instead of an edge, the layer overlaps and is not a
partition, and `make_folds()` warns. A point inside no block is assigned
to the nearest one, by distance to the polygon itself; a point the same
distance from two blocks, to one part in ten million or that margin,
goes to the lower row. The count of such points is warned about, unless
they sit within a millionth of the extent of a block, an edge that
reprojection or clipping moved by a rounding error. Blocks that hold no
point are dropped when `drop_empty_blocks = TRUE`, and `k` is lowered to
the number of blocks that hold points when that is smaller.

`params$n_blocks` is the number of blocks before empties were dropped
and `params$blocks_used` the number after (so with
`drop_empty_blocks = FALSE` the two are equal, and
`sum(params$block_sizes > 0)` is how many hold points); `params$grid_nx`
and `params$grid_ny` are `NA`. `params$block_scale` is the median, over
the blocks that hold points, of the side of the square with the block's
area. A polygon has no single edge length, so that is the length
compared with the autocorrelation range for the warning on blocks
smaller than the range.

Every set of cells
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
builds can be a block design here, including Voronoi cells around
`get_voronoi_seeds(method = "kmeans")` seeds, which adapt to the density
of the points.

## Fold balance

`"block_kfold"` and `"knndm"` balance the number of points per fold by
packing: blocks (for `"knndm"`, clusters of the samples) are taken
largest first and each goes to the fold with the fewest points so far,
ties broken at random. This is the longest-processing-time rule for
multiway partitioning, and `make_folds()` offers no search over other
packings. The remedy for an imbalance past `balance_tol` is the block
design: smaller blocks, or blocks that adapt to the density of the
points, passed through `blocks`. Density-adaptive blocks are smaller
where points are dense, so check `params$block_scale` against the
autocorrelation range as you would a grid.

`"knndm"` chooses its clusters, so for it `balance_tol` is a limit: a
clustering whose packing exceeds it is not a candidate, and the folds
returned are within it (`params$fold_sizes`, `params$balance_ratio`).

The other methods do not balance point counts. `"random_kfold"` is
balanced by construction (fold sizes differ by at most one);
`"buffered_loo"` and `"nndm"` hold out one point per fold;
`"leave_location_out"` gives each fold the same number of locations (to
within one), so folds differ by as much as the locations' sizes do.

## Numerical tolerances

A point on the edge two grid cells share goes to the lower one, the cell
with the lower index in the full grid (`params$blocks$source_row`; the
grid is numbered along each row from the bottom left, so the cell to the
left or below). A point is on an edge while it lies within one part in
1e10 of the largest coordinate of the points and the blocks from it
(0.44 mm at a UTM northing of 4,400 km), at most a millionth of the
diagonal of their bounding box and at least four times the spacing of
doubles at that coordinate, so points on the grid lines of a lattice are
blocked alike in any unit or origin. A point in no cell of a clipped
grid goes to the cell with the nearest centroid, and to the lower of two
whose centroids are as near to one part in ten million or that margin.
The same margin decides whether a point lies on a `boundary` or inside a
supplied block.

A block shorter than `block_size` by no more than one part in ten
million counts as that long, so a size of exactly the axis over a whole
number gives that number of blocks in any unit. A neighbour of a
`"buffered_loo"` point at the buffer distance is excluded, and so is one
over it by rounding: a squared distance within one part in ten million
of the buffer's square (a distance about one part in twenty million over
the buffer) counts as at it, so a lattice whose spacing divides the
buffer gives the same folds in any unit or origin.

NNDM counts distances within one part in ten million of each other as
equal, when the held-out distances are counted, compared with the target
and `phi`, and put in order. Equal ones are taken by the column of the
point (first coordinates within that part of the extent of the points
count as one), then its second coordinate, counted the same way, then
its row. Each point's neighbours are put in order by the same rule, so
of neighbours that are equally far away the same one is excluded first.
The folds are then the same in another unit or with another origin as
long as no coordinate moves by more than about a two-hundred-millionth
of the distance between neighbouring locations.

The clusterings of kNNDM read the coordinates centred on the middle of
their bounding box, divided by its diagonal and rounded to 2^-30 of it
(about a billionth), as the `"kmeans"` partition hands them to
[`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html); on a lattice
many merges and many assignments tie exactly, and the last bits of the
coordinates, which move with the unit and the origin, would settle them.
A point as near two k-means centres to one part in ten million goes to
the lower-numbered one: a squared distance ties with the least when it
is no more than the square of the least distance lengthened by a margin,
plus one part in ten million of the least squared distance, the margin
one part in 1e10 of the largest coordinate of the centres in that frame,
at most a millionth of the diagonal of their bounding box and at least
four times the spacing of doubles at that coordinate. The clusters are
then the same in another unit or with another origin, but for a location
within its own rounding noise of the middle of a step.

Points so far from the origin of their coordinates that doubles there
are spaced more than a millionth of the points' extent apart are held to
that spacing, and so are the edges of the blocks laid over them (see the
section "Coordinates far from the origin" in
[spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)).
With `method = "block_kfold"` or `"buffered_loo"` the folds are made,
and a warning of class `"spatialkit_xy_imprecise"` gives the spacing and
the extent.

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

[`plot_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_folds.md)
to map the folds,
[`fold_separation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fold_separation.md)
to measure how far the held-out points sit from their training data,
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
to score a model on them, and
[`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md)
to cross-validate at a ladder of block sizes.

Other folds:
[`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md),
[`fold_separation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fold_separation.md),
[`plot_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_folds.md)

## Examples

``` r
library(sf)
set.seed(1)
pts <- st_as_sf(
  data.frame(x = 5e5 + runif(30, 0, 1000), y = 5e6 + runif(30, 0, 1000)),
  coords = c("x", "y"), crs = 32632
)
folds <- make_folds(pts, k = 3, method = "block_kfold", seed = 42)
folds$assignment          # fold and block of each row
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
lengths(folds$folds[[1]]) # the train and test row IDs of fold 1
#> train  test 
#>    20    10 
folds$params$block_sizes  # points per block
#> [1] 2 2 2 3 5 3 4 4 5
folds$params$fold_blocks  # the blocks each fold holds
#> [[1]]
#> [1] 5 4 1
#> 
#> [[2]]
#> [1] 9 6 3
#> 
#> [[3]]
#> [1] 7 8 2
#> 

# Buffered leave-one-out: neighbours within 100 m left out of training.
loo <- make_folds(pts, method = "buffered_loo", buffer = 100)
loo$k                     # one fold per point
#> [1] 30

# Three folds matched to the distances from a prediction grid to
# clustered samples.
set.seed(2)
cx <- runif(6, 0, 1e4); cy <- runif(6, 0, 1e4)
xy <- do.call(rbind, lapply(1:6, function(i)
  cbind(cx[i] + rnorm(15, 0, 200), cy[i] + rnorm(15, 0, 200))))
cl <- st_as_sf(data.frame(x = 5e5 + xy[, 1], y = 5e6 + xy[, 2]),
               coords = c("x", "y"), crs = 32632)
grid <- st_make_grid(cl, n = c(20, 20), what = "centers")
km <- make_folds(cl, k = 3, method = "knndm", prediction_points = grid,
                 seed = 42)
km$params[c("chosen", "n_clusters", "wasserstein", "wasserstein_random")]
#> $chosen
#> [1] "clusters"
#> 
#> $n_clusters
#> [1] 18
#> 
#> $wasserstein
#> [1] 741.4627
#> 
#> $wasserstein_random
#> [1] 1108.152
#> 
```
