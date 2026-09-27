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
  tie_break = c("smallest_area", "first")
)
```

## Arguments

- features_sf:

  An sf object containing features to assign.

- polygons_sf:

  An sf or sfc polygonal layer.

- polygon_id_col:

  Name of the polygon identifier column. Default "poly_id".

- keep_unassigned:

  Logical; retain features that fall inside no polygon, carrying `NA` in
  the ID column. Default FALSE, which drops them.

- predicate:

  Binary spatial predicate function. Default sf::st_intersects. Not used
  when `largest` applies: sf then assigns polygon features by overlap
  area and never calls the predicate.

- largest:

  Logical; when `features_sf` is itself polygonal, keep the polygon with
  the largest overlap. Default TRUE. Ignored for point and line
  features. A feature that only touches the polygon layer (shares an
  edge or a corner with it, with no overlap area) has no largest overlap
  and is unassigned, in any CRS; with `largest = FALSE` the default
  `st_intersects` counts touching, so such a feature is assigned. A
  feature that overlaps two or more polygons by exactly the same area
  (to 9 significant digits; a square split evenly across a cell edge) is
  given to one of them by `tie_break` and counted in `"ties"`, so the
  choice does not depend on the order of the polygon rows. Invalid
  geometries (usually a self-intersecting ring), whose overlap is
  undefined, are repaired with
  [`sf::st_make_valid()`](https://r-spatial.github.io/sf/reference/valid.html)
  for the join, with a warning, and returned as they arrived. If the
  overlap still cannot be computed the function stops: falling back to
  `predicate` and `tie_break` would change the rule for every feature in
  the layer, so pass `largest = FALSE` to ask for that.

- tie_break:

  Strategy for resolving features that match multiple polygons (with
  `largest`, that overlap several polygons equally): `"smallest_area"`
  (default) keeps the polygon with the smallest area and, among polygons
  of equal area (a point on the shared edge of two grid cells), the one
  whose bounding-box centre is lowest, then leftmost, so the choice does
  not depend on the order of the rows; `"first"` keeps the first match
  (original order-dependent behavior).

## Value

An sf object with `polygon_id_col` attached, one row per input feature
(fewer if `keep_unassigned = FALSE` dropped unmatched ones), in the CRS
`features_sf` arrived in. A column of `features_sf` already called
`polygon_id_col` is dropped before the spatial join (with a warning), so
re-assigning an already-assigned layer replaces the old IDs and does not
fail. Every other column is kept, including one named like the polygons'
own ID column when that is read from a fallback such as `"id"` (a site
`id` joined to cells keyed by `id`). If *no* feature falls inside any
polygon the result is empty (or all-`NA` with `keep_unassigned = TRUE`)
and a warning is raised, since the usual cause is two layers in
different places (a CRS that could only be stamped, not reprojected).
The attribute `"ties"` records how many features matched more than one
polygon (with `largest`, overlapped several by exactly the same area)
and had the `tie_break` rule decide for them: a list with `n`, `which`
(their row positions in `features_sf`), `rule` and `n_rows` (the number
of rows returned, which the record was made for). A large `n` means the
polygon layer overlaps, and per-cell counts built from the result depend
on the rule. The record describes the rows this call returned and does
not survive subsetting: `joined[i, ]`, like
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html),
`slice()` or `arrange()` of it, is a plain layer with no `"ties"`
attribute, so nothing reports the parent's count for a different set of
rows.
[`sf::st_drop_geometry()`](https://r-spatial.github.io/sf/reference/st_geometry.html)
keeps the record, since the rows are the same; see
[`[.spatialkit_rows`](https://elkronos.github.io/gis_modeling_toolkit/reference/sub-.spatialkit_rows.md)
for what binding such data frames does.

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
with and overlap areas are planar. A copy of `features_sf` is
transformed for it, and the features come back with the coordinates they
arrived with. Otherwise (the polygons are in lon/lat, or carry no CRS)
the join runs in the CRS of `features_sf`.

## See also

[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
to build the polygon layer;
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
for the aggregation step that consumes the result.

Other aggregation:
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
