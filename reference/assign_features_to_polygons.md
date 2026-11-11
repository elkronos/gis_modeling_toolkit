# Assign features to polygons and attach a polygon ID

Joins an sf layer of input features to a polygon layer via spatial join.

## Usage

``` r
assign_features_to_polygons(
  features_sf,
  polygons_sf,
  polygon_id_col = "poly_id",
  keep_unassigned = FALSE,
  predicate = sf::st_intersects,
  largest = TRUE,
  tie_break = c("smallest_area", "first"),
  quiet = FALSE
)
```

## Arguments

- features_sf:

  An sf object containing features to assign. A layer whose geometries
  mix dimensions or carry M (XYM or XYZM) is joined in x and y, and the
  rows returned keep their Z and M (see Details).

- polygons_sf:

  An sf or sfc polygonal layer. Polygons whose geometries mix dimensions
  or carry M (XYM or XYZM) are read in x and y.

- polygon_id_col:

  Name of the polygon identifier column: the column of `polygons_sf` the
  IDs are read from, and the name they are returned under. Default
  "poly_id". When `polygons_sf` has no column of that name, the IDs are
  read from the first of `"poly_id"`, `"polygon_id"`, `"id"`,
  `"cell_id"` and `"grid_id"` it carries, and a message names the column
  that was used; when it has none of them, the polygons are numbered
  `1, 2, ...` in row order, and a message says so. The column the IDs
  are read from must hold one distinct, non-missing value per polygon: a
  duplicated or missing ID is an error, because features in different
  polygons would come back under one ID. That includes a unit stored as
  several rows that share its ID (a region split into `POLYGON` parts):
  merge the parts into one row per ID first, for example with
  `dplyr::summarise(dplyr::group_by(polygons_sf, <ID column>))`, which
  keeps the unit one cell. Numbering the rows would make each part a
  cell of its own.

- keep_unassigned:

  Logical; retain features that fall inside no polygon, carrying `NA` in
  the ID column. Default FALSE, which drops them with a warning of class
  `"spatialkit_unassigned"` that gives their number. Features that
  `tie_break` gives no polygon (see there) are kept or dropped with
  them.

- predicate:

  Binary spatial predicate function. Default sf::st_intersects. Not used
  when `largest` applies: sf then assigns polygon features by overlap
  area and never calls the predicate. With
  [`sf::st_intersects`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html)
  a POINT or MULTIPOINT feature with a point within a margin of a
  polygon is in it: one part in 1e10 of the largest coordinate of the
  polygons, 0.44 mm at a UTM northing of 4,400 km, at most a millionth
  of the diagonal of their bounding box and at least four times the
  spacing of doubles at that coordinate (in longitude and latitude, that
  many degrees of arc on the sphere while
  [`sf::sf_use_s2()`](https://r-spatial.github.io/sf/reference/s2.html)
  is on, and no margin while it is off). A point on the edge two
  polygons share is then in both in any unit and origin of the
  coordinates, and `tie_break` decides it, where the last bits of its
  coordinates put it inside one polygon or the other.

- largest:

  Logical; when `features_sf` is itself polygonal, keep the polygon with
  the largest overlap. Default TRUE. Ignored for point and line
  features. A feature that only touches the polygon layer (shares an
  edge or a corner with it, with no overlap area) has no largest overlap
  and is unassigned, in any CRS; with `largest = FALSE` the default
  `st_intersects` counts touching, so such a feature is assigned. A
  feature that overlaps two or more polygons by exactly the same area
  (to 9 significant digits; a square split evenly across a cell edge) is
  counted in `"ties"` and given to the one `tie_break` names, so the
  choice does not depend on the order of the polygon rows; where
  `"smallest_area"` names none, the feature is given no polygon (see
  `tie_break`). Invalid geometries (usually a self-intersecting ring),
  whose overlap is undefined, are repaired with
  [`sf::st_make_valid()`](https://r-spatial.github.io/sf/reference/valid.html)
  for the join, with a warning, and returned as they arrived. If the
  overlap still cannot be computed the function stops: falling back to
  `predicate` and `tie_break` would change the rule for every feature in
  the layer, so pass `largest = FALSE` to ask for that. In longitude and
  latitude with s2 on (sf's default), where s2 refuses a polygon of
  either layer as it arrived, the join by `predicate` meets that refusal
  too, and the error names the polygons and the transform instead (see
  Details). In longitude and latitude the validity check, the repair and
  the overlaps run on s2 whatever
  [`sf::sf_use_s2()`](https://r-spatial.github.io/sf/reference/s2.html)
  says; with s2 off the join by `predicate` runs with GEOS, and the
  error that names the polygons also offers `largest = FALSE`.

- tie_break:

  Strategy for resolving features that match multiple polygons (with
  `largest`, that overlap several polygons equally): `"smallest_area"`
  (default) keeps the polygon with the smallest area and, among polygons
  of equal area (a point on the shared edge of two grid cells), the one
  whose bounding-box centre is lowest, then leftmost, then the one whose
  vertices' mean is lowest, then leftmost (the two triangles of a split
  square share a bounding box), so the choice does not depend on the
  order of the rows. Areas are compared to 9 significant digits, and
  centres and means within one part in 1e10 of the largest coordinate of
  the polygons, at most a millionth of the diagonal of their bounding
  box and at least four times the spacing of doubles at that coordinate
  (with any `predicate`, and in longitude and latitude whatever
  [`sf::sf_use_s2()`](https://r-spatial.github.io/sf/reference/s2.html)
  says), so the choice does not depend on the unit or the origin of the
  coordinates either. Where polygons equal on every key are left (one
  polygon under two IDs, or a square and a rhombus of its area about its
  centre), the rule names no polygon: the feature is given none, `NA` in
  the ID column, counted in `"ties"` and listed in its `unresolved`, and
  a warning of class `"spatialkit_ties_unresolved"` gives how many such
  features there are (the condition carries `n` and `total`); it is kept
  or dropped as `keep_unassigned` says. Areas in longitude and latitude
  are measured on the sphere with s2, whatever
  [`sf::sf_use_s2()`](https://r-spatial.github.io/sf/reference/s2.html)
  says. A polygon among those compared whose geometry is invalid (a ring
  whose edges cross, whose lobes GEOS measures as cancelling) is
  measured on a copy repaired with
  [`sf::st_make_valid()`](https://r-spatial.github.io/sf/reference/valid.html),
  with a warning, and returned as it arrived; in degrees, one s2 refuses
  even repaired (a bowtie) stops the call with the error of class
  `"spatialkit_cells_degenerate"` (see Details). `"first"` keeps the
  first match, the polygon that comes first in the rows of `polygons_sf`
  (original order-dependent behavior). On cells built from the centres a
  profile scored, a point on the edge between two of them is placed by
  the profile's rule whichever is given (see Value).

- quiet:

  Logical; suppress this function's
  [`message()`](https://rdrr.io/r/base/message.html)s (the one naming a
  fallback ID column). It does not silence R warnings. Default `FALSE`.

## Value

An sf object with `polygon_id_col` attached, one row per input feature
(fewer if `keep_unassigned = FALSE` dropped unmatched ones, and features
`tie_break` gives no polygon), in the CRS `features_sf` arrived in. A
column of `features_sf` already called `polygon_id_col` is dropped
before the spatial join (with a warning), so re-assigning an
already-assigned layer replaces the old IDs and does not fail. Every
other column is kept, including one named like the polygons' own ID
column when that is read from a fallback such as `"id"` (a site `id`
joined to cells keyed by `id`). If *no* feature falls inside any polygon
the result is empty (or all-`NA` with `keep_unassigned = TRUE`) and a
warning is raised, since the usual cause is two layers in different
places (a CRS that could only be stamped, not reprojected). A layer
without a CRS is brought into the other layer's CRS as
[`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md)
does; beside a layer in longitude and latitude, one with an x outside
-180 to 360 or a y beyond 90 north or south cannot be in degrees, and
the call stops with an error that names `features_sf` or `polygons_sf`
and gives its box.

Features that fall inside no polygon are counted on the attribute
`"unassigned"`: a list with `n`, `which` (their row positions in
`features_sf`) and `n_rows`. With `keep_unassigned = FALSE` they are not
in the result, and a warning of class `"spatialkit_unassigned"` gives
their number (`tryCatch(spatialkit_unassigned = )` catches it without
matching the message; the condition carries `n` and `total`). The usual
causes are a polygon layer that stops short of the data (a coastline
drawn at another scale, a tessellation of seeds built without a
`boundary` or `cover`; see
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md))
and points outside the study area. A feature with empty geometry, or
with a coordinate that is not finite or lies beyond 1e150, or in a
geographic CRS a vertex at a latitude beyond 90 degrees (read as given,
before any move into the polygons' CRS), has no location to assign: it
matches no polygon, and the warning says how many of them there are.
With `keep_unassigned = TRUE` they are kept with an `NA` ID and no
warning is raised, unless every feature is unassigned. A polygon of
`polygons_sf` with a coordinate that is not finite or lies beyond 1e150
(read as given) has no location either: it takes no feature, and the
features inside the rest of its ring are unassigned.

The attribute `"ties"` records how many features matched more than one
polygon (with `largest`, overlapped several by exactly the same area)
and had the `tie_break` rule decide for them (on cells built from the
centres a profile scored, the profile's rule then places the points
among them, as the last paragraph says): a list with `n`, `which` (their
row positions in `features_sf`), `rule`, `n_rows` (the number of rows
returned, which the record was made for) and, when the rule gave some of
them no polygon (see `tie_break`), `unresolved`, the row positions of
those among `which`. A large `n` means the polygon layer overlaps, and
per-cell counts built from the result depend on the rule. The record
describes the rows this call returned and does not survive subsetting:
`joined[i, ]`, like
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html),
`slice()` or `arrange()` of it, is a plain layer with no `"ties"`
attribute, so nothing reports the parent's count for a different set of
rows.
[`sf::st_drop_geometry()`](https://r-spatial.github.io/sf/reference/st_geometry.html)
keeps the record, since the rows are the same; see
[`[.spatialkit_rows`](https://elkronos.github.io/gis_modeling_toolkit/reference/sub-.spatialkit_rows.md)
for what binding such data frames does. The `"unassigned"` record is
stamped and removed the same way.

When `polygons_sf` carries `attr(, "count_rule")`, the record of how its
cell count was set (cells from
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
on seeds placed from a
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
result or from a profile built with an `estimate`), the result carries
the same record.
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
reads it and checks that the rows and cells it is given are the ones the
count was scored on. It is a plain attribute: a row subset keeps it and
[`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html),
[`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html)
and
[`dplyr::arrange()`](https://dplyr.tidyverse.org/reference/arrange.html)
drop it, which is why the summary also looks for it on its `cells_sf`.
The profile scores a feature that is not a point at one point of it
(`coerce_to_points(, "auto")`), in the cell of the nearest centre; this
function gives a polygon feature to the cell it overlaps most. A feature
that straddles a cell boundary can therefore be assigned to another cell
than the one it was scored in, and the summary does not compare the two
for such a layer. To assign the features as they were scored, assign
their points:
`assign_features_to_polygons(coerce_to_points(features_sf, "auto"), polygons_sf)`.

When that record says the cells are the ones the profile scored
(`scored_cells` is `TRUE`: the cells
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
builds on a selection's seeds, or on the seeds a seeding function
returns for one) and holds the centres it scored
(`count_rule$scored$centres`), a point on the edge between two or more
of those cells is given the cell
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
scored it in: that of its nearest centre, measured in the CRS the record
names, and of centres the same distance away the lower-numbered one (a
distance counts as the same as the least when its square is no more than
the square of the least distance lengthened by a margin, plus one part
in ten million of the least squared distance; the margin is one part in
1e10 of the largest coordinate of the centres, at most a millionth of
the diagonal of their bounding box and at least four times the spacing
of doubles at that coordinate). Such points occur on gridded locations
(raster cell centres, a regular sampling grid), where a row can lie the
same distance from two centres and the polygons hold it in either cell
or in both. The rows each cell holds are then the rows its count was
scored on, and
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)'s
check of them reads `"match"`. A point counts as on the edge when it
lies within a margin of both cells: one part in 1e10 of the largest
coordinate of the cells, 0.44 mm at a UTM northing of 4,400 km, at most
a millionth of the diagonal of their bounding box and at least four
times the spacing of doubles at that coordinate. For a point that two
centres tie for while the cell holding it is that of the strictly nearer
one, the margin also takes in the point's distance from the boundary
between the two, measured in the CRS the record names: midway between
centres far apart the tie reaches about 1.25e-8 of their distance from
the boundary, which can exceed that margin.
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
draws such cells in the CRS the record names, around the centres as the
record holds them, whatever CRS the seeds came back in (unless its `crs`
names another projected CRS), so their corners carry only the rounding
of that computation, far less than the margin. Cells it returns in a
geographic `crs` are moved back into the record's CRS for the join (see
Details). In a UTM zone the round trip through degrees moves a corner by
about 2e-9 m; in the equal-area CRS the profile works in for a layer
wider than a zone it can move one by more than the margin (how far
depends on the PROJ build), so a point on such an edge can keep the
polygon that holds it, and the summary's check then reads `"differ"`. A
point farther from its nearest centre's cell keeps the polygon that
holds it, and so do features that are not points, points without a CRS
beside a record with one or with a CRS beside a record without one, and
a join in longitude and latitude (polygons in degrees beside a record in
degrees). A centre's cell is the one cell that holds it and no other
centre; a centre outside a boundary that clipped its cell away from it
has the cell, among those that hold no centre, whose inside is nearest
it. Cells built from other seeds leave most centres without a cell, and
their points to the polygons. Polygons without the record, such as
polygons the caller brings, and cells whose record says `scored_cells`
is `FALSE` or `NA` (other cells built from the count, such as a lattice
sized by it) are placed by the polygons alone.

## Details

This is the second step of the package's pipeline: it labels every
observation with the cell it falls in, which is what
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
then aggregates over. Prefer it to
[`sf::st_join()`](https://r-spatial.github.io/sf/reference/st_join.html)
when the join has to be *unambiguous*. It resolves features matching
several polygons by an explicit `tie_break` rule instead of silently
duplicating rows, so the assigned layer keeps one row per input feature
and cell-level counts mean what they say.

The join runs in the CRS of `polygons_sf` whenever that CRS is
projected, so cell edges are the straight lines the cells were drawn
with and overlap areas are planar. Cells in lon/lat whose record of
their count says they are the cells a profile scored, and holds its
centres in a projected CRS (the Voronoi cells
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
builds on the centres a profile scored and returns in a geographic
`crs`; see Value), were drawn in the record's CRS, and the join runs
there. A copy of `features_sf` is transformed for it, and the features
come back with the coordinates they arrived with. Otherwise (other
polygons in lon/lat, such as a hex or square lattice sized by a
selection or the cells of centres refitted to every point, or polygons
without a CRS) the join runs in the CRS of `features_sf`. In longitude
and latitude, with s2 on
([`sf::sf_use_s2()`](https://r-spatial.github.io/sf/reference/s2.html),
sf's default), sf joins on the sphere with s2, which refuses a polygon
that is not valid there, such as one with two corners it reads as one
point or with edges that cross; cells returned in lon/lat can hold two
corners that differ only by rounding and fall on one point in degrees
(see `crs` in
[`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md)).
A join by `predicate` then stops with an error of class
`"spatialkit_cells_degenerate"` that gives the rows of the polygons s2
refuses (the condition carries them as `rows`), names the cause and says
to transform `polygons_sf` to a projected CRS, where the join runs. The
largest-overlap join repairs such polygons first (see `largest`), and
stops with that error when s2 still refuses a polygon it repaired, such
as one whose hole crosses its outer ring. A polygon of `features_sf`
that s2 refuses stops the join by `predicate` with an error of class
`"spatialkit_features_degenerate"`, with the same fields, that gives its
rows, quotes s2's reason and says to transform `features_sf` to a
projected CRS, where the join runs; so does one s2 refuses even after
the repair of the largest-overlap join. When that repair makes s2 accept
the polygons of both layers and s2 then refuses a piece of an overlap of
the two layers (a figure eight, one vertex listed twice, is repaired
into two parts that touch at a point, and the piece can hold that point
twice), the join stops with the error for the layer whose polygons s2
refused as they arrived; with s2 on (sf's default), the join by
`predicate` meets those polygons too. When s2 refused no polygon of
either layer as it arrived, the error says to pass `largest = FALSE` or
to project both layers. Whether s2 refuses such a piece follows the last
bits of the vertices s2 computes for it, which differ between platforms
and s2 versions: a figure eight pinched on the edge of a polygon of the
other layer can give a piece that s2 refuses on macOS on arm64 and
accepts on Linux on x86-64. Where s2 accepts the piece, the join runs
and gives each feature the polygon of its largest overlap. In longitude
and latitude the largest-overlap join runs on s2, and the areas
`tie_break` compares are measured on the sphere with s2, whatever
[`sf::sf_use_s2()`](https://r-spatial.github.io/sf/reference/s2.html)
says, as the package's other areas in degrees are; the join by
`predicate` runs with the engine the session names. With s2 off that
engine is GEOS, which reads these polygons in degrees, so the
largest-overlap join's error for a polygon s2 refuses also offers
`largest = FALSE`.

sf reads a layer in the dimension of its first geometry, so s2, GEOS and
[`sf::st_transform()`](https://r-spatial.github.io/sf/reference/st_transform.html)
read a layer that mixes XY and XYZ (or XYM) geometries wrongly;
[`rbind()`](https://rdrr.io/r/base/cbind.html) of an XY and an XYZ layer
gives one. When either layer mixes dimensions, the join runs on a copy
of it without Z and M, and the rows come back with the geometry they
arrived with, Z and M included. The same holds, with no warning, for a
layer whose geometries carry M (XYM or XYZM), which GEOS does not read:
the join is defined in x and y. Features without a CRS that mix
dimensions and are moved from lon/lat into the polygons' CRS come back
moved and in XY, with a warning of class `"spatialkit_zm_dropped"`,
unless their only geometries of another dimension than the first are
empty and not points. Features without a CRS in one dimension that carry
M come back with their own geometry: with the CRS they were given, or
moved from lon/lat with
[`sf::st_transform()`](https://r-spatial.github.io/sf/reference/st_transform.html),
which keeps their M.

Layers so far from the origin of their coordinates that doubles there
are spaced more than a millionth of the extent of the two layers apart
(0.125 at 1e15, beside a layer a few kilometres across) are joined with
their edges and points held to that spacing, so a point that lies within
it of an edge near the origin can fall in the polygon beside it. The
join runs, and a warning of class `"spatialkit_xy_imprecise"` gives the
spacing and the extent (the condition carries them as `spacing` and
`extent`).

## See also

[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
to build the polygon layer;
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
for the aggregation step that consumes the result.

Other aggregation:
[`cell_values_at()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cell_values_at.md),
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
pts <- st_as_sf(
  data.frame(x = 5e5 + runif(20, 0, 100), y = 5e6 + runif(20, 0, 100),
             val = rnorm(20)),
  coords = c("x", "y"), crs = 32632
)
bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
  c(5e5, 5e6), c(5e5 + 100, 5e6), c(5e5 + 100, 5e6 + 100),
  c(5e5, 5e6 + 100), c(5e5, 5e6)
))), crs = 32632))
grid <- create_grid_polygons(bnd, target_cells = 9, type = "square")
assigned <- assign_features_to_polygons(pts, grid)
table(assigned$poly_id)
#> 
#> 1 2 3 4 5 6 7 8 9 
#> 1 2 3 3 2 3 1 3 2 
```
