# Random seed generation within a polygonal boundary

Draws `k` seed points uniformly at random inside `boundary`, ignoring
where the observations are. Reach for this when the cells should cover
the study area evenly (so that sparsely sampled ground still gets its
own cells and is visibly under-sampled in the results) instead of
concentrating resolution where the data already are, which is what
[`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)
does. It is also the honest choice for a null or sensitivity comparison:
re-running an analysis over several random seedings shows how much of a
result depends on one particular tessellation.

## Usage

``` r
voronoi_seeds_random(boundary, k, set_seed = NULL)
```

## Arguments

- boundary:

  An sf or sfc polygonal object. One that mixes XY and XYZ (or XYM)
  geometries is read without Z and M, since sf reads a layer in the
  dimension of its first geometry. So is one whose geometries carry M
  (XYM or XYZM), which GEOS refuses. Rows that GEOS (in a projected CRS)
  or s2 (in longitude and latitude) refuses are repaired before the rows
  are unioned, as in
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md);
  valid rows are used as they are. A row that is empty, or has a
  coordinate that is not finite or lies beyond 1e150 (read as given, and
  again after any move into another CRS), covers no ground and is left
  out; a boundary with no such row left is refused with an error that
  says so. A boundary in a geographic CRS with a vertex at a latitude
  beyond 90 degrees, which is no place on the Earth, is refused with an
  error of class `"spatialkit_cells_degenerate"` that names its rows and
  says to check its CRS. A row the repair leaves without a polygon (a
  ring with no height: every vertex on one line, or in degrees on one
  parallel) is left out as well, and a boundary with no polygon left is
  refused. In longitude and latitude, a boundary whose bounding box is
  over 10,000 times its area (small pieces far apart), which sf does not
  sample whole, is sampled piece by piece, the points shared among the
  pieces in proportion to their areas; a piece whose own box is over
  10,000 times its area stops the call with an error that says to
  transform `boundary` to a projected CRS. A boundary whose bounding box
  has no height or no width (three vertices on one parallel) is refused.
  In longitude and latitude s2 reads the edges as great circles, which
  leave the bounding box read from the coordinates, and sf (without the
  lwgeom package) draws in that box: a boundary whose box has less than
  half its area on the sphere (three vertices on or near one parallel)
  is refused too. That error, and the error for a box with no height
  whose vertices lie on one parallel other than the equator, say to
  transform `boundary` to a projected CRS in which a parallel is a
  curve, such as the one
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  picks for it (a polar one, such as EPSG:3413 or 3031, for a ring
  around a pole), where the region has an area and lies in its box. sf
  draws the points in the bounding box and keeps those inside, so it
  draws `k` times the ratio of the box to the area. In any CRS, a
  boundary that would take more than 10^6 points drawn in its box, or 10
  `k` if more (small pieces far apart, a thin strip at a slant), is
  sampled piece by piece in the same way (or whole, in draws of at most
  10^6 points in its box, where that takes less time), and a piece that
  would still take more is drawn in such draws. Where the points drawn
  in all would exceed 2e7 (or 10 `k`, if more), the call stops before
  any draw with an error that gives the ratio and says to ask for fewer
  seeds.

- k:

  Integer; number of random seeds.

- set_seed:

  Optional integer RNG seed. Default `NULL`: the seeds are drawn from
  the session's random-number stream, so consecutive calls give
  different seedings and
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before a call makes
  it reproducible. Pass a number to get the same seeds whatever that
  stream holds; the caller's stream is then left as it was. The default
  used to be 456, which made every call return the same "random"
  seeding, even inside a loop over
  [`set.seed()`](https://rdrr.io/r/base/Random.html).

## Value

An sf object of `k` random POINTs (fewer only in the degenerate case
above), with `seed_id` and `method = "random"` columns matching
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
and `attr(, "cover")`, the polygons of `boundary`:
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
clips the cells to it (grown by 2 percent of its diagonal) when no
`boundary` is passed there, so the cells fill the area the seeds were
drawn in and not only the seeds' own hull.

## Details

Seeds are drawn uniformly inside the polygon. When a draw falls short of
`k` it is topped up with further uniform draws from the same polygon, so
the result has exactly `k` seeds; only a geometry that still yields too
few after ten top-ups returns fewer, and that shortfall is logged.

## See also

Other tessellation:
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md),
[`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md),
[`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md),
[`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md),
[`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md),
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
[`plot_tessellation_map()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_tessellation_map.md),
[`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)

## Examples

``` r
library(sf)
bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
  c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
))), crs = 32632))
set.seed(1)
seeds <- voronoi_seeds_random(bnd, k = 10)
nrow(seeds)   # 10
#> [1] 10
seeds
#> Simple feature collection with 10 features and 2 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 6.178627 ymin: 17.65568 xmax: 94.46753 ymax: 99.19061
#> Projected CRS: WGS 84 / UTM zone 32N
#>                     geometry seed_id method
#> 1  POINT (26.55087 20.59746)       1 random
#> 2  POINT (37.21239 17.65568)       2 random
#> 3  POINT (57.28534 68.70228)       3 random
#> 4  POINT (90.82078 38.41037)       4 random
#> 5  POINT (20.16819 76.98414)       5 random
#> 6  POINT (89.83897 49.76992)       6 random
#> 7  POINT (94.46753 71.76185)       7 random
#> 8  POINT (66.07978 99.19061)       8 random
#> 9   POINT (62.9114 38.00352)       9 random
#> 10 POINT (6.178627 77.74452)      10 random
# Another call is another seeding; set_seed pins one.
identical(st_coordinates(voronoi_seeds_random(bnd, k = 10, set_seed = 7)),
          st_coordinates(voronoi_seeds_random(bnd, k = 10, set_seed = 7)))
#> [1] TRUE
```
