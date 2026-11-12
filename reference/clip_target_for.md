# Build the area to tessellate within, from points or a boundary

Builds the area that the hex and square cells of
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
are clipped to. With a `boundary` it is that boundary, row for row;
without one it is the axis-aligned bounding box of `points_sf` in the
working CRS. Either can be buffered by `expand`.

## Usage

``` r
clip_target_for(points_sf, boundary = NULL, expand = 0, quiet = FALSE)
```

## Arguments

- points_sf:

  An sf object with POINT/MULTIPOINT geometry. Points that mix XY and
  XYZ (or XYM), or that carry M, are read in x and y, since sf reads a
  layer in the dimension of its first geometry and GEOS refuses M. A
  point that is empty, or has a coordinate that is not finite or lies
  beyond 1e150, has no location and takes no part in the bounding box.
  With no `boundary`, points none of which has a location (as given, or
  after the move of lon/lat points into a projected CRS) are refused
  with an error that says so. Points that take the buffer, so far from
  the origin of their CRS that doubles there are spaced more than a
  millionth of the target's size apart (0.125 at 1e15, beside a disk of
  radius 1), give a buffer drawn at that rounding, with a warning of
  class `"spatialkit_xy_imprecise"` that gives the spacing and the
  radius. Where the buffer comes back empty or with no area, as a radius
  below that spacing can make it, the call stops with an error of class
  `"spatialkit_xy_imprecise_error"`. Shift such points nearer the
  origin, or pass `boundary`.

- boundary:

  Optional polygonal sf object. One with no CRS, given with points that
  have one, is interpreted in the points' own CRS, with a warning. With
  points in a projected CRS it is read as
  [`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md)
  does: coordinates that look like lon/lat are taken as EPSG:4326 and
  reprojected, anything else is stamped with the points' CRS. With
  lon/lat points it is read as lon/lat when its coordinates fit the
  lon/lat envelope, and refused with an error otherwise. Rows without a
  location, rows the repair leaves without a polygon, vertices at a
  latitude beyond 90 and mixed dimensions are handled as "Boundaries" in
  [spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)
  says.

- expand:

  A single non-negative number. Below 1 it is a fraction of the longer
  side of the bounding box (of `boundary` when one is given, otherwise
  of `points_sf`); 1 or more is a distance in the units of the CRS the
  target is built in. Default 0 (no buffer). Lon/lat input is first
  projected to a local projected CRS (see Value), so both forms are in
  that CRS's units.

- quiet:

  Logical; suppress this function's progress
  [`message()`](https://rdrr.io/r/base/message.html)s. It does not
  silence R warnings, nor the package's console log echo (see
  [`spatialkit_quiet`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)
  for that). Default `FALSE`.

## Value

An sf polygon layer representing the clip target. For lon/lat input the
layer is returned in the automatically selected local projected CRS, not
the input CRS; a message reports this unless `quiet = TRUE`.

## Details

Use it to build the `boundary` that `method = "hex"` and `"square"`
require, to check that a study-area polygon contains the observations
before tessellating, or to pass the same envelope to
[`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md)
and
[`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md)
so that two tessellations of one dataset cover the same ground.

It is not the target
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
derives on its own: `method = "voronoi"` without a boundary clips to the
convex hull of the points (or of their `cover`, for seeds) buffered by 2
percent of its diagonal, and there `expand` is always a distance. A
bounding box over a non-rectangular point cloud includes corners with no
data, so a hex or square grid laid over it has cells that hold no
points; pass the study-area polygon when there is one.

## Points with no extent

Points that all share one x or one y, or nearly so (the short side of
their bounding box below a millionth of the long side, or below about
1.5e-154, a length whose square doubles cannot hold), get a small buffer
around them instead of a box. When every point sits at one location, or
the share of their spread the buffer would take (5 percent, or a
fractional `expand`) is below that length, the box has no extent to take
a share of, and the target is a disk of radius 1 unit of the CRS around
it (1 m in a metre CRS, 1 km in a kilometre CRS), unless `expand` is 1
or more, which is then the radius; the message says so. Pass `boundary`
for a target of real size.

## See also

Other cell building:
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md),
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
  data.frame(x = 5e5 + runif(30, 0, 100), y = 5e6 + runif(30, 0, 100)),
  coords = c("x", "y"), crs = 32632
)
# No boundary: the bounding box, expanded by 10% of the extent
box <- clip_target_for(pts, expand = 0.1, quiet = TRUE)
st_area(box)
#> 12054.18 [m^2]
```
