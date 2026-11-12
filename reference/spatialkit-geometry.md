# How input geometry is read

The rules the package applies to the spatial layers it is given: how two
layers are brought into one CRS, how a layer without a CRS is read, how
mixed dimensions and rows without a location are handled, what a
boundary needs, and what s2 refuses in longitude and latitude. A
function's own page says how it applies them: which of its arguments
they cover, and whether a row it leaves out is reported by a warning or
in the log.

## Coordinate reference systems

A function that takes two layers, such as points and a boundary or
features and cells, brings them into one CRS before it reads either.
[`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md)
is the rule in its plain form: when both layers have a CRS, one is
transformed into the other's; a layer without one is read by the test in
"Layers without a CRS" and then reprojected or stamped, with a warning.

A function that measures distances moves a lon/lat layer into the
projected CRS
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
chooses from the layer's extent: the UTM zone for a local extent,
otherwise the candidate that distorts distances least (its page gives
the cases). A layer that already has a projected CRS is used as it is.
Distances, and arguments given as distances, are in the units of the CRS
the function works in, which its page names. For a lon/lat layer the
function projects, that is the CRS
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
chose, whose unit is the metre.

## Layers without a CRS

A layer without a CRS is read by a test of its coordinates. When their
bounding box lies within longitude -180 to 180 and latitude -90 to 90,
and they either span more than one unit on some axis or carry the
decimal places of degrees, they are taken as EPSG:4326, with a warning.
Otherwise they are left as they are. A small planar survey inside that
box is read as degrees, so set the CRS of such a layer with
`sf::st_crs(x) <- ...` before passing it.

Given a CRS to read the layer in (`target_crs`, or the CRS of a second
layer), coordinates the test takes for degrees are reprojected from
EPSG:4326, and the others are stamped with that CRS and not moved. A
stamp is a relabel: check that the coordinates are in that CRS. A CRS in
longitude and latitude is never stamped on coordinates that cannot be in
it. An x outside -180 to 360, or a y beyond 90 north or south, stops the
call with an error that gives the bounding box of the coordinates.

## Mixed dimensions, Z and M

sf reads a layer in the dimension of its first geometry, and GEOS
refuses M. A function that moves a layer into another CRS and returns
it, or reduces it to points, reads a layer that mixes XY and XYZ (or
XYM) geometries, as [`rbind()`](https://rdrr.io/r/base/cbind.html) of an
XY and an XYZ layer gives, in x and y. It returns that layer in XY, with
a warning of class `spatialkit_zm_dropped` (its page names any
exception). A layer in one dimension, XYM and XYZM included, is moved
with its Z and M, and so is a layer whose only geometries of another
dimension than the first are empty and not points.

A layer a function only reads for its x and y, such as a boundary or the
points and cover of a builder, is read in x and y when it mixes
dimensions or its geometries carry M (XYM or XYZM).
[`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md)
reduces a geometry that carries M in x and y as well, with the same
warning; its page gives the exceptions.

## Rows without a location

A row has no location when its geometry is empty, or one of its
coordinates is not finite or lies beyond 1e150. Such a row takes no part
in the choice of a projection, and its coordinates none in the lon/lat
test. Moved into another CRS, it comes back empty, in its own geometry
type.
[`sf::st_transform()`](https://r-spatial.github.io/sf/reference/st_transform.html)
alone would drop the vertex it cannot move, or put a coordinate beyond
1e150 at a finite place, and return the rest of the feature.

In a geographic CRS a vertex at a latitude beyond 90 is no place on the
Earth. A row with one takes no part in the choice of a projection
either.
[`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md),
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
and the fits, the points of the k-means seeders, the features of
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md),
the pooled columns of
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
and
[`cell_values_at()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cell_values_at.md),
and
[`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md)
read such a row as a row without a location.
[`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md)
and
[`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md)
say so with a warning of class `spatialkit_lat_beyond_90` that names the
rows.

A row that has a location and that the target CRS cannot place in full
comes back empty, or without the vertices that CRS cannot place, with a
warning of class `spatialkit_rows_unplaced` that names the rows.

What a function does with such a row (drops it, keeps it with `NA`, or
leaves it out of one step) is on its page.

## Boundaries

A boundary is a polygonal sf or sfc layer: the study area that cells are
clipped to, that seeds are drawn in, or that fold blocks are cut from.

- When exactly one of the data layer and the boundary has a CRS, the
  other is read in it, with a warning. A boundary without a CRS given
  with lon/lat points is read as lon/lat when its coordinates fit
  longitude and latitude, and refused with an error otherwise.

- A row that is empty, or has a coordinate that is not finite or lies
  beyond 1e150, read as given and again after any move into another CRS,
  covers no ground and is left out. A boundary with no such row left is
  refused with an error that says so.

- A boundary in a geographic CRS with a vertex at a latitude beyond 90
  is refused with an error of class `spatialkit_cells_degenerate` that
  names its rows and says to check its CRS.

- A row the repair of the boundary leaves without a polygon (a ring with
  no height: every vertex on one line, or in degrees on one parallel) is
  left out as well, and a boundary with no polygon left is refused.

- A boundary that mixes dimensions or carries M is read in x and y.

[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
adds the rules for drawing points inside a boundary, and
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
the rule for a block the boundary only touches.

## Polygons in degrees

With s2 on
([`sf::sf_use_s2()`](https://r-spatial.github.io/sf/reference/s2.html),
sf's default), sf reads polygons in longitude and latitude on the
sphere, and s2 refuses a polygon that is not valid there: one whose
edges cross, or with two corners it reads as one point. Cells returned
in lon/lat can hold two such corners. This happens where four Voronoi
seeds stand the same distance from one point, and the transform to
degrees puts two corners that differ only by rounding at one point (see
`crs` in
[`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md)).

A function that joins, measures or unions such a polygon in degrees
stops with an error that names its rows, quotes s2's reason and says to
transform the layer to a projected CRS, such as the one
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
picks for it. The class is `spatialkit_cells_degenerate` for cells and
boundaries, and `spatialkit_features_degenerate` for features being
assigned or reduced to points. Some steps repair the polygons first,
such as the largest-overlap join of
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
and the union of a boundary in
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
and stop only when s2 refuses a repaired polygon.

A hex or square grid built with a geographic `crs` is laid in degrees
([`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md)).
Voronoi cells asked for in a geographic `crs` are built in a projected
CRS and then transformed.

## Coordinates far from the origin

Far from the origin of a CRS, doubles are spaced more widely: 0.125
apart at 1e15. Where they are spaced more than a millionth of a layer's
extent apart, edges and points are held to that spacing, and a point
within it of an edge can fall on the other side. The functions that
build cells, blocks and joins still run, with a warning of class
`spatialkit_xy_imprecise` that gives the spacing and the extent;
[`clip_target_for()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clip_target_for.md)
stops with an error of class `spatialkit_xy_imprecise_error` where its
buffer comes back empty or with no area. Those that search for
neighbours
([`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md),
[`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md))
give the same warning where a distance lies past a tie by less than the
rounding. Subtract a constant from the coordinates to bring them nearer
the origin.

## See also

[spatialkit-conditions](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-conditions.md)
for every condition class named here,
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
for the choice of a projected CRS, and
[`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md)
for aligning two layers.

Other data preparation:
[`[.spatialkit_rows()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sub-.spatialkit_rows.md),
[`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md),
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md),
[`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md),
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
