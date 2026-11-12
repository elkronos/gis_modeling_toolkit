# Build Voronoi cells around points

Gives every location in the study area to its nearest input point, so
each point gets one cell. Use it when the observations themselves define
the cells (sampling sites, monitoring stations, service points): cell
size then follows sampling density, with small cells where the points
are dense and large ones where they are sparse. Use
[`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md)
instead for equal-area cells, or for a resolution that does not depend
on where the data are.

## Usage

``` r
create_voronoi_polygons(
  points_sf,
  boundary = NULL,
  expand = 0,
  clip = TRUE,
  keep_duplicates = FALSE,
  crs = NULL,
  quiet = FALSE,
  cover = NULL
)
```

## Arguments

- points_sf:

  An sf object with POINT/MULTIPOINT geometries. Points with no CRS
  whose coordinates look like lon/lat (the heuristic
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  applies, with its warning) are taken as EPSG:4326 and projected, as
  lon/lat points are. Points that mix XY and XYZ (or XYM), or that carry
  M, are read in x and y, since sf reads a layer in the dimension of its
  first geometry and GEOS refuses M; the attributes the seeding
  functions set on them are kept.

- boundary:

  Optional polygonal sf object. When exactly one of `points_sf` and
  `boundary` has a CRS, the other is interpreted in it, with a warning.
  CRS-less points, and a CRS-less boundary given with projected points,
  are read as
  [`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md)
  does (lon/lat-looking coordinates are reprojected from EPSG:4326,
  others are stamped); CRS-less points that do not look like lon/lat
  cannot take a geographic boundary's CRS, and are refused with an
  error. A CRS-less boundary given with lon/lat points (or with CRS-less
  points taken as lon/lat) is read as lon/lat when its coordinates fit
  the lon/lat envelope, and refused with an error otherwise. Rows
  without a location, rows the repair leaves without a polygon, vertices
  at a latitude beyond 90 and mixed dimensions are handled as
  "Boundaries" in
  [spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)
  says.

- expand:

  Numeric; absolute distance, in the working CRS's units, by which the
  boundary (or the hull derived from the points) is grown before the
  diagram is built. With `clip = TRUE` the cells are clipped to the
  grown boundary, so they reach `expand` beyond the study area, and a
  point up to `expand` outside it gets a cell and an `index` value. The
  grown boundary is the one returned as `boundary`.

- clip:

  Logical; clip the cells to the boundary. Default `TRUE`.

- keep_duplicates:

  Logical. With `FALSE` (the default), points are merged before the
  diagram is built when they round to one point in coordinates centred
  on the middle of their bounding box, divided by its diagonal and
  rounded to 2^-30 of it, so points within about a billionth of the
  extent of each other count as coincident in any unit or origin and
  share one cell. With `TRUE` only points at exactly one location are
  merged. Every point is indexed to the cell that holds it.

- crs:

  Optional target CRS: anything
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html)
  accepts, including an sf or sfc layer, whose CRS is used. A projected
  CRS is the working CRS. A geographic one (EPSG:4326, say) is the CRS
  the result is returned in: the cells are built in the local projected
  CRS
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  picks for the points, so they are nearest-point cells on the ground,
  and are then transformed, with long edges densified. Seeds that are
  the centres a profile scored are built in the CRS the profile worked
  in instead (see "Seeds a profile scored"). Where four seeds stand the
  same distance from one point, as on gridded locations, two cells can
  meet there at two corners that differ only by rounding, which the
  transform to degrees can put at one point, and sf's spherical engine
  (s2) refuses such a cell.
  [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
  joins the cells of the centres a profile scored in the CRS they were
  built in, while they carry the record that says so
  (`scored_cells = TRUE`).
  [`sf::st_area()`](https://r-spatial.github.io/sf/reference/geos_measures.html)
  and
  [`sf::st_union()`](https://r-spatial.github.io/sf/reference/geos_combine.html)
  stop on such a cell in degrees with s2's error. The package's other
  joins and measures of it stop with an error of class
  `"spatialkit_cells_degenerate"` ("Polygons in degrees" in
  [spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)
  gives the exceptions). Transform the cells to a projected CRS first.

- quiet:

  Logical; suppress this function's progress
  [`message()`](https://rdrr.io/r/base/message.html)s. It does not
  silence R warnings, nor the package's console log echo (see
  [`spatialkit_quiet`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)
  for that). Default `FALSE`.

- cover:

  Optional sf or sfc layer the cells have to hold, read when there is no
  `boundary`: the observations (points), or a polygon. The boundary is
  then derived from it instead of from `points_sf`: the convex hull of
  the cover's points, or its polygons, grown by 2 percent of the
  diagonal. This is what seeds need, because seeds lie inside the data
  they stand for and their own hull stops short of it (see "Details").
  When `NULL`, the cover the seeding functions attach is used
  (`attr(points_sf, "cover")`, set by
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
  [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)
  and
  [`voronoi_seeds_random()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_random.md)).
  A cover with no CRS is read in the coordinates of `points_sf`, and one
  in another CRS is transformed. A generator outside the cover keeps its
  ground (the hull of the generators is added); a cover that holds none
  of the generators is set aside with a warning. With a `boundary`, the
  boundary sets the cells and `cover` is only checked against them: the
  number of its points that fall outside every cell is returned as
  `params$cover_outside`, with a warning of class
  `"spatialkit_unassigned"` when it is not 0 (a polygon cover is not
  checked). A cover, given or attached, that mixes XY and XYZ (or XYM)
  geometries, or that carries M, is read in x and y. A row of the cover
  with a coordinate that is not finite or lies beyond 1e150, read as
  given before any move, is read as an empty row.

## Value

A list with these elements:

- `cells`:

  The cells, one polygon per unique point, with a stable `cell_id`.

- `index`:

  One `cell_id` per row of `points_sf`, and `NA` for a point that falls
  outside every cell: outside the study area (grown by `expand` when it
  is positive), and farther than a thousandth of the median cell width
  from the nearest cell, so a summary built from it counts only the
  points the tessellation covers. Such points are counted in
  `params$n_outside`, and a warning of class `"spatialkit_unassigned"`
  gives their number. For seeds (see "Details") `index` holds the cell
  each seed generates, whether the seed lies inside the boundary or not,
  and `NA` only for a seed whose cell the clip removed whole or without
  a location (an empty geometry, or a coordinate that is not finite or
  lies beyond 1e150); `params$n_outside` counts those, the warning names
  them as seeds, and `params$snapped` lists only the seeds whose entry
  is the cell they were snapped to.

- `boundary`:

  The boundary the cells were built in: the one supplied or derived,
  grown by `expand`.

- `method`:

  `"voronoi"`.

- `params`:

  The settings and diagnostics of the build, among them:

  - `boundary_source`: where the boundary came from, `"boundary"` (the
    argument), `"cover"` (the `cover` argument), `"seeds' cover"` (the
    cover attached to `points_sf`) or `"points"` (the hull of
    `points_sf` itself);

  - `cover_outside`: the number of `cover` points outside every cell
    when a point `cover` was given with a `boundary`, and `NULL`
    otherwise;

  - `partition`: present when `points_sf` is a layer of seeds that
    records the kind of partition that placed it
    (`attr(points_sf, "partition")`, see
    [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)),
    the same list, so the cells say what built them. Seeds of plain
    k-means record nothing, and `params` then has no such element;
    neither has it for points that are not seeds;

  - `count_rule`: present when `points_sf` carries the record of how its
    count was set (`attr(points_sf, "count_rule")`: the seeds of a
    [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
    result, or seeds a seeding function placed from one), the same list,
    with `scored_cells = FALSE` when the seeds are fewer or more than
    its count. The cells carry it too, as `attr(cells, "count_rule")`;
    [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
    copies it to the rows it assigns, and
    [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
    checks its result against it.

## Details

[`sf::st_voronoi()`](https://r-spatial.github.io/sf/reference/geos_unary.html)
builds the diagram. This function projects lon/lat input, builds and
buffers an envelope so edge cells are bounded, clips to `boundary`,
restores the point-to-cell correspondence that
[`st_voronoi()`](https://r-spatial.github.io/sf/reference/geos_unary.html)
scrambles, and numbers the cells with stable `cell_id` values.

The generators are the points' vertices, not the features. A MULTIPOINT
feature with several vertices therefore gets one cell per vertex, and
its `index` entry is the smallest `cell_id` among the cells it touches;
the others are referenced by no feature. Other functions in the package
([`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md),
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md))
reduce such a feature to its centroid instead, so cast to POINT, or take
centroids, first if one cell per feature is what you want.

With no `boundary` the cells are clipped to a boundary derived from the
data. For observations that is their own convex hull grown by 2 percent
of its diagonal, so every observation has a cell. Seeds are different:
k-means centres, random points in a study area and the centres a
resolution profile scored all lie inside the ground they stand for, and
cells clipped to the hull of the seeds leave the edge of that ground
uncovered. The boundary is therefore derived from the `cover` when there
is one: the argument, or the cover the seeding functions attach to their
result. A layer that looks like seeds and has neither a boundary nor a
cover (seeds bound with [`rbind()`](https://rdrr.io/r/base/cbind.html),
which drops the attribute, or built by hand) gets the seeds' own hull
and a warning of class `"spatialkit_no_cover"`.

A seed stands for its cell, so `index` gives each seed the cell it
generates. A k-means seed is a mean of observations and can lie outside
a `boundary` that is not convex, in a bay or a bend of a coast (see
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)).
Its cell, the part of the study area nearer it than any other seed, is
built, and its `index` entry is that cell's `cell_id`, although the seed
lies in no cell. A seed is `NA` only when every part of the study area
is nearer another seed, so the clip removed its cell whole, or when it
has no location (its geometry is empty, or a coordinate is not finite or
lies beyond 1e150). The points are seeds when a seeding function
returned them
([`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
[`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md),
[`voronoi_seeds_random()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_random.md),
[`rbind()`](https://rdrr.io/r/base/cbind.html) of such layers included),
when they are the seeds of a
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
result, and when a `cover` is passed with them. A selection's seeds are
seeds by the attributes they carry, which
[`rbind()`](https://rdrr.io/r/base/cbind.html) and a write to a file
drop. Any other points are indexed by the cell that holds them, so a
point outside the study area is `NA`, unless it lies within a thousandth
of the median cell width of a cell, whose ID it then takes (see `index`
in
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md));
that includes observations that carry a selection's record of the count,
as
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
returns them. Seeds placed by hand, and a selection's seeds that lost
their attributes, are seeds once they pass through
`get_voronoi_seeds(method = "provided", seeds = )`.

## Seeds a profile scored

The seeds of a
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
result, and the seeds
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
or
[`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)
return for one, are the centres
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
scored, and their record of the count (`attr(, "count_rule")`) says so
(`scored_cells = TRUE`) and holds the centres in the CRS the profile
worked in. (On a profile without an estimate fitted to a subsample the
seeds are those centres refitted to every point, which carry
`attr(, "refit")`; their record says `scored_cells = FALSE`, and their
cells are built as other seeds' are. So are seeds fewer or more than the
record's count, such as a selection's seeds cut down to some of its
rows: they are not the centres the profile scored, and their record is
handed on with `scored_cells = FALSE`, as
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
(`method = "provided"`) marks such seeds.) Unless a projected `crs` is
given, their cells are built in that CRS, whatever CRS the seeds arrive
in (degrees, or the CRS of a boundary they were placed in), and returned
in it, or in a geographic `crs` when one is given. When the seeds are
those centres (as many, each at a different centre, within a thousandth
of its distance from the next nearest), each is first put at its centre
as the record holds it, so rounding the seeds took on in another CRS
does not move an edge. The cells are then the cells each row was scored
in, and
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
gives a point on an edge between two of them the cell the profile scored
it in.

## Points far from the origin of their CRS

Points so far from the origin of their CRS that doubles there are spaced
more than a millionth of their extent apart (0.125 at 1e15, beside a
layer a few kilometres across) have the cells' edges and their own cells
read at that rounding, so two points that close can fall in one cell.
The cells are built, and a warning of class `"spatialkit_xy_imprecise"`
gives the spacing and the extent.

## See also

Other cell building:
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md),
[`clip_target_for()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clip_target_for.md),
[`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md),
[`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md),
[`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md),
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
[`plot_tessellation_map()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_tessellation_map.md),
[`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md),
[`voronoi_seeds_random()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_random.md)

## Examples

``` r
library(sf)
set.seed(1)
pts <- st_as_sf(
  data.frame(x = 5e5 + runif(15, 0, 100), y = 5e6 + runif(15, 0, 100)),
  coords = c("x", "y"), crs = 32632
)
res <- create_voronoi_polygons(pts, quiet = TRUE)
res$cells   # one polygon per unique point, with stable cell_id
#> Simple feature collection with 15 features and 1 field
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 500003.5 ymin: 4999999 xmax: 500097.1 ymax: 5000102
#> Projected CRS: WGS 84 / UTM zone 32N
#> First 10 features:
#>                          geometry cell_id
#> 1  POLYGON ((500024.2 5000020,...       1
#> 2  POLYGON ((500030.5 5000080,...       2
#> 3  POLYGON ((500007.1 5000040,...       3
#> 4  POLYGON ((500024.2 5000020,...       4
#> 5  POLYGON ((500048.1 5000053,...       5
#> 6  POLYGON ((500025 5000064, 5...       6
#> 7  POLYGON ((500055.9 5000081,...       7
#> 8  POLYGON ((500037.4 5000019,...       8
#> 9  POLYGON ((500071.5 5000085,...       9
#> 10 POLYGON ((500044.6 5000033,...      10
res$index   # cell_id assignment for each input point
#>  [1]  5  6  9 15  2 13 14 11  8  1  3  4 10  7 12
```
