# Build a tessellation (Voronoi, Delaunay triangles, hex grid, or square grid)

The single entry point for turning a point pattern into analysis
regions, and the first step of the package's pipeline. It wraps the four
tessellation methods behind one interface that handles CRS projection,
clipping and stable cell identifiers consistently, and returns the cell
layer together with the point-to-cell index that
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
and
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
consume. Prefer it to the individual constructors whenever you might
want to compare methods: the return shape does not change with `method`,
so swapping `"voronoi"` for `"hex"` costs one argument.

## Usage

``` r
build_tessellation(
  points_sf,
  boundary = NULL,
  method = c("voronoi", "triangles", "hex", "square"),
  approx_n_cells = NULL,
  cellsize = NULL,
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

  An sf object with POINT/MULTIPOINT geometry.

- boundary:

  Polygonal sf/sfc study area. **Required** for `method = "hex"` and
  `method = "square"`, which have no extent of their own to lay a grid
  over and error without it; supply the study-area polygon, or build one
  from the points with
  [`clip_target_for()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clip_target_for.md).
  **Optional** for `method = "voronoi"` and `method = "triangles"`,
  which derive their extent from the points themselves (for `"voronoi"`,
  from `cover` when there is one) and use `boundary` only to clip the
  result when `clip = TRUE`. When `points_sf` holds seeds and not the
  observations, pass the study area here or the observations as `cover`:
  seeds lie inside the data, and cells clipped to the seeds' own hull
  leave its edge uncovered. When exactly one of `points_sf` and
  `boundary` has a CRS, the other is interpreted in it, with a warning.
  CRS-less points, and a CRS-less boundary given with projected points,
  are read as
  [`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md)
  does; CRS-less points that do not look like lon/lat cannot take a
  geographic boundary's CRS and are refused with an error. A CRS-less
  boundary given with lon/lat points is read as lon/lat when its
  coordinates fit the lon/lat envelope, and refused with an error
  otherwise. When neither has one, both are read by the lon/lat
  heuristic of
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md):
  taken as EPSG:4326 and projected when the points look like degrees (a
  boundary whose coordinates do not fit the lon/lat envelope is then
  refused with an error), otherwise left in the same unnamed planar
  space.

- method:

  One of "voronoi", "triangles", "hex", "square".

- approx_n_cells:

  Approximate number of cells. Read by `method = "hex"` and `"square"`
  only: `"voronoi"` grows one cell per input point and `"triangles"` one
  triangle per neighbouring triple, so neither has a count to set, and
  both warn that the argument was ignored. For a Voronoi cell count,
  place the seeds with
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
  and tessellate those, passing the study area as `boundary` or the
  observations as `cover` (seeds made from a boundary or from
  `sample_points` carry their cover with them). For hex grids the target
  is adjusted for packing density; the actual count after clipping to an
  irregular boundary may differ noticeably. Besides a number, this
  accepts what the level-selection step returned: the integer vector of
  ranked candidates from
  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  (its first element is used), a
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  result (its `$best`), or a
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md).
  A profile built without an `estimate` is read with
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  by `cp`; one built with `estimate` is read by `criterion = "auto"`,
  the rule that goes with its estimate: a size for pooled values, the
  lowest leave-one-out error for plain means. The count used is returned
  as `params$approx_n_cells` and where it came from as
  `params$approx_n_cells_from` (`NULL` for a plain number). A count read
  off a selection or a profile raises a warning of class
  `"spatialkit_bound_chosen"` when its optimum sits on a bound of the
  ladder and one of class `"spatialkit_no_structure"` when one cell for
  the whole layer scores as well (for a count set by a size: when cells
  of that size predict no better than one mean; the count stands); a
  profile without an estimate and with no `cp` to read is an error
  (class `"spatialkit_no_default_criterion"` when another criterion is
  finite). See
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md).
  A count read off a profile or selection is a number of k-means cells:
  every one occupied, and small where the points are dense. A lattice
  lays that many equal cells over the whole boundary, so on clustered
  points many of them hold no point (about half, on six clusters in a
  square); on evenly spread points it matches. `params$cells_occupied`
  and `params$cells_empty` report how the points filled the grid, and a
  count that came from a profile or selection warns when fewer than
  three quarters of it are occupied. For cells that follow the points,
  seed a Voronoi tessellation with
  `get_voronoi_seeds(method = "kmeans", n = <the selection>, sample_points = <the points>)`:
  given the
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  result or the profile itself, it returns the centres of the partition
  the profile scored, whose Voronoi cells are the cells the criteria
  judged. Seeds made that way carry the hull of `sample_points` as their
  cover, so the cells built on them hold every point; made with neither
  `sample_points` nor a boundary, they carry the hull of the points the
  profile was computed on.

- cellsize:

  Numeric cell size, in the units of the working CRS. Read by
  `method = "hex"` and `"square"` only; the other two methods warn that
  it was ignored. When both `cellsize` and `approx_n_cells` are given,
  `cellsize` wins and `approx_n_cells` is ignored with a logged warning;
  supply one or the other. With a geographic `crs` it is in that CRS's
  degrees, and the grid is laid in degrees.

- expand:

  Buffer distance, in the working CRS's units, by which the Voronoi
  boundary is grown before the diagram is built. Applied by
  `method = "voronoi"` only; the `"hex"`, `"square"` and `"triangles"`
  methods ignore it (the value you passed is still echoed back in
  `params$expand`). With `clip = TRUE` the cells are clipped to the
  grown boundary, which is the one returned as `boundary`: cells reach
  `expand` beyond the study area, and a point up to `expand` outside it
  is indexed.

- clip:

  Logical; clip to boundary.

- keep_duplicates:

  Logical. Has no effect on the cells or the index: coincident points
  are merged before a Voronoi diagram or a Delaunay triangulation is
  built either way, and every one of them is indexed to the cell they
  share.

- crs:

  Optional target CRS: anything
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html)
  accepts, including an sf or sfc layer, whose CRS is used. A projected
  CRS is the working CRS. A geographic one (EPSG:4326, say) is the CRS
  the result is returned in: the cells are built in the local projected
  CRS
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  picks for the points, indexed there, and then transformed with long
  edges densified, so Voronoi cells are nearest-point cells on the
  ground and grid cells are laid in metres rather than degrees. The
  exception is a hex or square grid sized by `cellsize`, which is in
  degrees and so is laid in degrees. Whenever that local CRS is picked
  for lon/lat points, or CRS-less ones taken as lon/lat (no `crs`, or a
  geographic one), a hex or square grid with a boundary is laid in it
  unless it distorts areas across the boundary by more than 1 percent
  (Web Mercator over a near-global extent, say); the grid is then laid,
  and the points indexed, in the equal-area CRS
  `ensure_projected(purpose = "area")` picks for the boundary, with a
  logged warning, as
  [`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md)
  does, so the cells stay equal-area. Seeds that carry the record of how
  their count was set (`attr(, "count_rule")`) were scored in the CRS
  the profile worked in, each row in the cell of its nearest centre.
  When the record says they are the centres it scored and holds them
  (the seeds of a
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  result, or those a seeding function returns for one),
  `method = "voronoi"` builds their cells in that CRS, with no `crs` or
  a geographic one, whatever CRS the seeds arrive in, and returns them
  there or in the geographic `crs`; see
  [`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md),
  section "Seeds a profile scored". Voronoi cells built in another
  projected CRS, named here, are the nearest-seed cells of that CRS, so
  a row near a boundary can fall in the neighbouring cell, and
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  then reports that the cells are not the ones scored.

- quiet:

  Logical; suppress this function's progress
  [`message()`](https://rdrr.io/r/base/message.html)s. It does not
  silence R warnings, nor the package's console log echo (see
  [`spatialkit_quiet`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)
  for that). Default `FALSE`.

- cover:

  Optional sf or sfc layer the cells have to hold: the observations
  (points), or a polygon. Read by `method = "voronoi"` when there is no
  `boundary`: the boundary is then derived from the cover (the convex
  hull of its points, or its polygons, grown by 2 percent of the
  diagonal) instead of from the hull of `points_sf`. Pass it when
  `points_sf` holds seeds. When `NULL`, the cover the seeding functions
  attach to their result is used (`attr(points_sf, "cover")`: the
  boundary the seeds were drawn in, or the hull of the points they were
  clustered from), so `build_tessellation(voronoi_seeds_kmeans(obs, k))`
  holds every observation without it;
  [`rbind()`](https://rdrr.io/r/base/cbind.html) of seed layers drops
  that attribute. Seeds with neither a boundary nor a cover are
  tessellated over their own hull, with a warning of class
  `"spatialkit_no_cover"`. A cover with no CRS is read in the
  coordinates of `points_sf`; one in another CRS is transformed. With a
  `boundary` (and for the other methods) the cover sizes nothing: when
  it is a point layer, the number of its points that fall outside every
  cell is returned as `params$cover_outside`, with a warning of class
  `"spatialkit_unassigned"` when it is not 0. See
  [`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md).

## Value

A list with components:

- `cells`:

  An sf polygon layer, one row per cell. It always carries a `cell_id`
  column; the `"hex"` and `"square"` methods additionally carry
  `poly_id`, which holds the same values.

- `index`:

  Integer vector of `cell_id` values, one per row of `points_sf`, and
  `NA` for a point that falls inside no cell, that is, one outside the
  study area. Only a point within a thousandth of the median cell width
  of a cell is snapped to it. That covers points sitting exactly on a
  shared edge, and leaves points outside the study area as `NA`. A
  summary built from `index` therefore counts only the points the
  tessellation actually covers. When any point is `NA`, a warning of
  class `"spatialkit_unassigned"` gives the number
  (`tryCatch(spatialkit_unassigned = )` catches it; the condition
  carries `n` and `total`).

- `boundary`:

  The boundary used (possibly derived and/or reprojected, and for
  `"voronoi"` grown by `expand`).

- `method`:

  The method actually used.

- `params`:

  The parameters the tessellation was built with, plus `snapped`, the
  record of that nearest-cell repair: a list with `n`, `which` (row
  positions in `points_sf`) and `distance` (how far outside every cell
  each sat, in CRS units); `n_outside`, the number of points of
  `points_sf` that fall in no cell (the `NA`s of `index`); and
  `cover_outside`, the number of `cover` points in no cell when a point
  `cover` was checked against the cells (`NULL` otherwise). For
  `"voronoi"` also `boundary_source`, where the boundary came from:
  `"boundary"`, `"cover"`, `"seeds' cover"` or `"points"`, and, when
  `points_sf` is a layer of seeds that records the kind of partition
  that placed it (see
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)),
  `partition`: that record, a list with `kind`, `init`, `density_power`
  and `density_k` and, for seeds a seeding function placed itself,
  `iter`, `iter_max` and `converged`. The cells are the Voronoi cells of
  the seeds whatever kind placed them. For `"hex"` and `"square"` also
  `cells_occupied` and `cells_empty`, the number of cells that hold at
  least one point and that hold none. `count_rule`, the record of how
  the count was set, is present on two paths, and the cells then carry
  it too as `attr(cells, "count_rule")`. For `"voronoi"`, when
  `points_sf` carries one (`attr(points_sf, "count_rule")`: the seeds of
  a
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  result, or seeds
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
  or
  [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)
  placed from one): the record as the seeds carry it, whose
  `scored_cells` says whether these are the cells the profile scored.
  For `"hex"` and `"square"`, when `approx_n_cells` is a selection that
  carries a record, or a profile built with an `estimate`, and no
  `cellsize` is given: that record with `scored_cells` `FALSE`, since a
  lattice at that count is not the partition the profile scored.
  [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
  copies the record to the rows it assigns, and
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  checks its result against it.

## Details

Which method to reach for. `"voronoi"` gives one cell per point, so
resolution follows sampling density. That is the choice when the
observations themselves define the regions. `"hex"` and `"square"` give
equal-area cells on a fixed grid, so cell size is a decision you make,
and the data does not make it for you; hexagons avoid the axis-aligned
artefacts of squares and have uniform neighbour distances. `"triangles"`
returns the Delaunay triangulation, useful for interpolation and
adjacency work; it is not meant as an aggregation unit.
[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
will suggest a cell count from the spatial structure of the data.

`"voronoi"` and `"triangles"` are built on the points' vertices: a
MULTIPOINT feature with several vertices gets one cell (or triangle
corner) per vertex, and its `index` entry is the smallest `cell_id`
among the cells it touches. See
[`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md).

## See also

Other tessellation:
[`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md),
[`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md),
[`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md),
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
  data.frame(x = 5e5 + runif(20, 0, 100), y = 5e6 + runif(20, 0, 100)),
  coords = c("x", "y"), crs = 32632
)
tess <- build_tessellation(pts, method = "voronoi", quiet = TRUE)
tess$cells
#> Simple feature collection with 20 features and 1 field
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 500003.6 ymin: 4999999 xmax: 500101.8 ymax: 5000096
#> Projected CRS: WGS 84 / UTM zone 32N
#> First 10 features:
#>                          geometry cell_id
#> 1  POLYGON ((500016.9 5000038,...       1
#> 2  POLYGON ((500021.3 5000077,...       2
#> 3  POLYGON ((500016.9 5000038,...       3
#> 4  POLYGON ((500007.8 5000047,...       4
#> 5  POLYGON ((500044.6 5000090,...       5
#> 6  POLYGON ((500021.3 5000077,...       6
#> 7  POLYGON ((500043.4 5000044,...       7
#> 8  POLYGON ((500043.6 5000044,...       8
#> 9  POLYGON ((500053.2 5000027,...       9
#> 10 POLYGON ((500070.6 5000087,...      10
```
