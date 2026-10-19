# Cell values at other locations

Gives each location of `newdata` the value of the cell it lies in, read
from a
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
result: the pooled value when the summary holds one, the plain mean
otherwise, with the interval for one observation when the summary holds
it. Beside each value it returns the location's distance from the
nearest observation, since a cell's value describes the part of the cell
where its observations lie.

## Usage

``` r
cell_values_at(
  newdata,
  values,
  points,
  response_var,
  id_col = "poly_id",
  seeds = NULL
)
```

## Arguments

- newdata:

  The locations: an sf layer. A feature that is not a point is reduced
  to one with
  [`coerce_to_points(, "auto")`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md).
  No column of it is read. With no CRS while `points` has one, the call
  stops.

- values:

  The
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  result for `points`. It must hold the ID column and a cell value of
  `response_var`: the pooled value (`estimate = "pooled"`) or the mean
  (the default `agg_funs`). With `seeds = NULL` it must also carry the
  cell polygons, which it does when it was made with `cells_sf`.

- points:

  The assigned layer `values` was computed from: an sf layer with the ID
  column and `response_var`, and at least one scored row.

- response_var:

  Name of the summarised column wanted: the response, or a predictor.

- id_col:

  Name of the ID column, found as
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  finds it: when `points` has no column of that name, the first of
  `"poly_id"`, `"polygon_id"` and `"cell_id"` it has. It cannot have the
  name of another column of the result (`value`, `estimate`, `obs_lo`,
  `obs_hi`, `filled`, `filled_from`, `dist_nearest`, `dist_ratio`).

- seeds:

  `NULL`, or the seeds the cells were built on: an sf POINT layer with a
  column named as the ID column that gives each seed's cell. For cells
  from `tess <- build_tessellation(seeds)` that is
  `seeds$cell_id <- tess$index` when the ID column is `"cell_id"`, after
  which a seed outside a boundary that is not convex needs its cell's ID
  (see "Seeds outside the boundary"). Two seeds at different locations
  with one ID are an error; seeds that share a location share its cell,
  as `index` gives them. `NA` marks a seed without a cell.

## Value

A data frame with one row per row of `newdata`, in its order:

- the ID column:

  The cell the location lies in, under the name of the ID column; `NA`
  for a location in no cell.

- `value`:

  The cell's value, or for a filled row the value of the cell in
  `filled_from`.

- `estimate`:

  `"pooled"` or `"mean"`: which value of the summary `value` is.

- `obs_lo`, `obs_hi`:

  The interval for one observation, when `values` holds it.

- `filled`:

  `TRUE` when the location has no cell, or its cell neither a scored row
  nor a value, and the value is that of the cell of the nearest scored
  row.

- `filled_from`:

  For a filled row, the ID of the cell the value came from; `NA`
  otherwise.

- `dist_nearest`:

  Distance to the nearest location of a scored row, in the units of the
  working CRS.

- `dist_ratio`:

  `dist_nearest` over the spacing; `NA` when the scored rows sit on one
  location.

A row of `newdata` with an empty geometry, or with a coordinate that is
not finite or lies beyond 1e150, has `NA` in every column. The two ID
columns have the class of the ID column of `values`; when a cell that
`values` has no row for has an ID that class cannot hold (a string
beside numbers), both are character.

What produced the frame is in `attr(, "values_at")`, a list of class
`"spatialkit_values_at"` with a print method: `assignment` (`"seeds"` or
`"polygons"`), `crs` (the working CRS), `column` (`response_var`),
`estimate`, `obs_level` (the level of the interval columns; `NULL`
without them, and `NA` when `values` holds them and no longer says at
what level), `n` (the scored rows), `locations` (their distinct
locations), `spacing` (what `dist_ratio` divides by; `NA` with one
location), `filled` (how many rows of `newdata` were filled) and
`count_rule` (the summary's record of how the cell count was set, or
`NULL`). It describes the call: a subset of the frame's rows keeps it as
it is.

The function raises no classed warning of its own. It passes on what
[`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md),
a CRS transformation and
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
raise, but for one warning: the one that function raises when no
location at all falls in a polygon is neither raised nor logged, since
every such row has `filled = TRUE`, which says the same row by row.

## Details

The function works in six steps.

1.  *The working CRS.* The CRS the summary measured in: the one in
    `attr(values, "count_rule")` when that record names a projected CRS,
    else the one in `attr(values, "pooled")`, else the CRS of `points`
    when it is projected and the one
    [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
    picks for it when it is not. `points`, `newdata` and `seeds` are
    transformed to it. When neither `points` nor `newdata` has a CRS the
    coordinates are used as they are, unless those of `points` look like
    longitude and latitude: both are then read as EPSG:4326, as
    [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
    reads `points`, with
    [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)'s
    warning.

2.  *The scored rows.* The rows of `points` with a finite
    `response_var`, a cell ID and finite coordinates. Their distinct
    locations are the observations every distance below is taken from.

3.  *The cell of a location.* With `seeds`, the cell of the nearest
    seed. That is the Voronoi cell before its clip to the cover, so
    every location has one, unless its nearest seed has an `NA` ID. A
    location nearest such a seed has no cell, even where the polygon
    that holds it has scored rows.
    [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
    gives `NA` in `index` to a seed that lies in no cell: one whose cell
    the clip removed, and one outside a boundary that is not convex,
    whose clipped cell exists (see "Seeds outside the boundary").
    Without `seeds`, the polygon of `values` that holds the location, by
    [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
    with its default `tie_break` and the summary's record of the count,
    so a location on the edge between two cells the profile scored is
    placed as that function places such a row (see "Ties"); a location
    in no polygon has no cell. Polygons without a CRS are read as that
    function read them for the rows of `points`: the locations are
    handed to it in the CRS of `points`, whatever CRS `newdata` arrived
    in.

4.  *The value.* The cell's pooled value of `response_var`,
    `..pooled_<prefix><response_var>` with the prefix `resp_` or
    `pred_`, when `values` has that column; otherwise its mean,
    `<prefix>mean_<response_var>`. The column `estimate` says which. For
    a column summarised as response and as predictor the response's
    columns are read. The interval is the cell's
    `..obs_lo_<prefix><response_var>` and
    `..obs_hi_<prefix><response_var>` when `values` has them.

5.  *The fill.* A location with no cell, or whose cell holds no scored
    row and has no value in `values`, takes the value and interval of
    the cell of its nearest scored row. `filled` is then `TRUE` and
    `filled_from` is that cell's ID. The ID column still holds the
    location's own cell, as step 3 finds it, so the two columns together
    say where a location lies and where its value came from. A cell that
    holds a scored row is not filled. Where `values` has no value for
    such a cell (the summary could not compute the pooled columns, or
    has no row for the cell), its locations have `value = NA` and
    `filled = FALSE`; a filled value is `NA` in the same case, when the
    cell it is taken from has none.

6.  *The distance.* `dist_nearest` is the distance from the location to
    the nearest location of a scored row, in the units of the working
    CRS. `dist_ratio` is that distance over the layer's spacing: the
    median, over the distinct locations of the scored rows, of the
    distance to the nearest other one.

Cell IDs are matched as text. The ID column of `values` is character
when the summary joined a `cells_sf` whose IDs were of another class
than the assigned layer's, while `points` and `seeds` keep the class
they had; the IDs of all three are written out as text before any match,
and both ID columns returned have the class of the ID column of `values`
wherever that class can hold every ID (see Value).

## Seeds outside the boundary

A k-means seed can lie outside a boundary that is not convex, in a bay
or a bend of a coast (see
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)).
Its cell, clipped to the boundary, holds ground and observations, but
the seed lies in no cell, so `seeds$cell_id <- tess$index` gives it
`NA`. The locations nearest it then have `NA` in the ID column and are
filled from the cell of their nearest scored row, which can be another
cell than the one that holds them. Every point of a Voronoi cell is
nearer its own seed than any other, so a point inside each cell that no
seed lies in finds that cell's seed. In a projected CRS:

    seeds$cell_id <- tess$index
    own  <- tess$cells[!tess$cells$cell_id %in% tess$index, ]
    near <- sf::st_nearest_feature(sf::st_point_on_surface(sf::st_geometry(own)),
                                   seeds)
    seeds$cell_id[near] <- own$cell_id

Each location inside the boundary then has the cell that holds it, and
none is filled for want of a seed. A seed whose cell the clip removed
keeps `NA`.

## Ties

A location the same distance from two seeds takes the one that comes
first in `seeds`. For the seeds of a
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
result that is the lower-numbered centre, the rule
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
scores with. When the summary's record of the count holds the centres
the profile scored and `seeds` are those centres (the seeds a selection
carries, or those a seeding function returns for it), in any order and
any CRS, the distances are taken to the centres as the record holds
them, a tie goes to the lower-numbered centre, and each centre gives the
ID of the seed that stands at it. Seeds count as the centres when they
are as many and each lies at a different centre, within a thousandth of
its distance from the next nearest one, so a seed that came back from
another CRS with its coordinates rounded, or seeds sorted by their cell,
do not settle a tie another way. A location on the edge between two such
cells then takes the cell
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
gives a row there, with `seeds` and without. Other seeds, among them a
subset of the centres, are searched as given. A location the same
distance from two locations of scored rows takes the one with the lower
x coordinate, then the lower y; among scored rows at one location that
lie in different cells, the cell that comes first in the order of the
IDs (ascending for numbers, the C-locale order for strings, the order of
the levels for a factor). Neither rule depends on the order of the rows
of `points`, or on whether the FNN package is installed. A tie is a tie
between distances computed from the coordinates in the working CRS, so
on gridded locations the same layer in another CRS or unit can put a
location on a boundary in the other cell.

## Away from the observations

The value and the interval describe locations among the observations.
The interval's two ends are quantiles of leave-one-out errors, each the
error of predicting an observation from the others around it, so it says
nothing of a location far from every observation: a gap in the layer, or
ground beyond its edge. Away from the observations no interval here is
calibrated. `dist_nearest` and `dist_ratio` say how far each location is
from the data. No threshold is applied to either and no location is
flagged or left out.
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
and
[`fold_separation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fold_separation.md)
are the tools for deciding where a prediction applies.

A simulation shows how the interval's coverage falls with `dist_ratio`
(`dev/sim-pooled-values.R` in the package's repository; "Pooled values"
in
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
describes it). Each of 200 layers has 1,500 clustered observations in
333 cells, 300 further observations held out at random among them and
300 held out as a block, which leaves a gap in the layer. Over the
held-out observations of all layers, the interval at `obs_level = 0.90`
held 91.5 percent of those with a `dist_ratio` below 1, 89.2 percent
from 1 to 2, 83.9 percent from 2 to 4, 77.9 percent from 4 to 8 and 69.5
percent above 8. It held 89.6 percent of the observations held out at
random and 75.1 percent of those in the block. A second simulation, on
fields whose level jumps at boundaries, gave this function's interval to
held-out observations; "Intervals: which one answers which question" in
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
gives its coverage, which in cells of one training observation fell
short of `obs_level` on one kind of field.

## Without the FNN package

The function makes three nearest-neighbour searches: each location's
nearest scored location, each location's nearest seed, and each scored
location's nearest other one (for the spacing). With FNN each uses a
kd-tree. Without it each compares every pair, and a search that would
compare more than `getOption("spatialkit.exact_pairs_max", 2e9)` pairs
stops the call first, with an error of class `spatialkit_work_limit`
that names the search (the fields `work`, `limit` and `what`). The
search for the spacing does not depend on the number of locations asked
for: above about 44,700 distinct scored locations the call stops without
FNN however few they are. The result is the same with and without the
package.

## See also

[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md),
which computes `values`;
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md),
whose `index` gives each seed's cell;
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
and
[`fold_separation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fold_separation.md)
for where a prediction applies.

Other aggregation:
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md),
[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md),
[`kriging_adequacy()`](https://elkronos.github.io/gis_modeling_toolkit/reference/kriging_adequacy.md),
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md),
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md),
[`summary.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.resolution_profile.md)

## Examples

``` r
library(sf)
set.seed(1)
n <- 400
east  <- runif(n, 0, 1e4)
north <- runif(n, 0, 1e4)
all_pts <- st_as_sf(
  data.frame(x = 5e5 + east, y = 5e6 + north,
             income = 30 + 5 * sin(east / 2000) + 5 * cos(north / 1500) +
               rnorm(n)),
  coords = c("x", "y"), crs = 32632
)
pts    <- all_pts[1:350, ]
newpts <- all_pts[351:400, ]

# The count from a size, the cells, the observations in them, the values,
# all at the defaults: pooled values at 4.5 observations per cell.
prof  <- resolution_profile(pts, "income")
sel   <- select_resolution(prof)
seeds <- get_voronoi_seeds(method = "kmeans", n = sel, sample_points = pts)
tess  <- build_tessellation(seeds, method = "voronoi")
#> create_voronoi_polygons(): deriving boundary from the cover the seeds carry (its hull or polygons, grown by 2% of the diagonal).
asg   <- assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id")
out   <- summarize_by_cell(asg, "income", id_col = "cell_id",
                           cells_sf = tess$cells, obs_level = 0.90)

# The value at other locations, each placed in the cell of its nearest
# seed.  `tess$index` gives the cell of each seed.
seeds$cell_id <- tess$index
at <- cell_values_at(newpts, out, asg, "income", id_col = "cell_id",
                     seeds = seeds)
head(at)
#>   cell_id    value estimate   obs_lo   obs_hi filled filled_from dist_nearest
#> 1      41 28.48055   pooled 26.02536 31.09205  FALSE          NA     238.4435
#> 2       5 37.30667   pooled 34.85149 39.91818  FALSE          NA     115.6546
#> 3      64 22.93542   pooled 20.48023 25.54692  FALSE          NA     271.2353
#> 4      28 38.91665   pooled 36.46147 41.52815  FALSE          NA     367.1291
#> 5      75 20.96892   pooled 18.51374 23.58042  FALSE          NA     227.2415
#> 6      10 34.21084   pooled 31.75566 36.82234  FALSE          NA     351.8853
#>   dist_ratio
#> 1  0.9403784
#> 2  0.4561211
#> 3  1.0697034
#> 4  1.4478917
#> 5  0.8961997
#> 6  1.3877727
attr(at, "values_at")
#> Cell values at other locations: income, the pooled value
#>   cell        : of the nearest seed (the Voronoi cell before its clip to the cover)
#>   interval    : 90% for ONE OBSERVATION (obs_lo, obs_hi), centred on the value
#>   scored rows : 350 rows at 350 locations, measured in EPSG:32632
#>   spacing     : 253.6, the median distance from a scored location to the
#>                 nearest other; dist_ratio is dist_nearest over it
#>   filled      : no location
#>   count       : 78 cells, set by size (the summary's record, in `count_rule`)
#> Away from the observations no interval here is calibrated.

# Without the seeds a location is placed in the polygon that holds it.  One
# beyond the layer's edge lies in no polygon: it is filled from the cell of
# the nearest observation, and `dist_ratio` says how far that is.
beyond <- st_sf(geometry = st_sfc(st_point(c(5e5 + 2e4, 5e6 + 5e3)), crs = 32632))
cell_values_at(rbind(st_sf(geometry = st_geometry(newpts[1:2, ])), beyond),
               out, asg, "income", id_col = "cell_id")
#>   cell_id    value estimate   obs_lo   obs_hi filled filled_from dist_nearest
#> 1      41 28.48055   pooled 26.02536 31.09205  FALSE          NA     238.4435
#> 2       5 37.30667   pooled 34.85149 39.91818  FALSE          NA     115.6546
#> 3      NA 25.14089   pooled 22.68570 27.75239   TRUE          78   10204.7354
#>   dist_ratio
#> 1  0.9403784
#> 2  0.4561211
#> 3 40.2456524
```
