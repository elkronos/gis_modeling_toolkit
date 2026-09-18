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

  Binary spatial predicate function. Default sf::st_intersects.

- largest:

  Logical; when `features_sf` is itself polygonal, keep the polygon with
  the largest overlap. Default TRUE. Ignored for point and line
  features, and silently dropped if the `predicate` does not support it
  ([`sf::st_intersects`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html)
  does).

- tie_break:

  Strategy for resolving features that match multiple polygons:
  `"smallest_area"` (default) keeps the polygon with the smallest area,
  `"first"` keeps the first match (original order-dependent behavior).

## Value

An sf object with `polygon_id_col` attached, one row per input feature
(fewer if `keep_unassigned = FALSE` dropped unmatched ones), in the CRS
`features_sf` arrived in. Any column of `features_sf` whose name would
collide with the polygon ID column is dropped before the spatial join
(with a warning), so re-assigning an already-assigned layer replaces the
old IDs rather than failing. If *no* feature falls inside any polygon
the result is empty (or all-`NA` with `keep_unassigned = TRUE`) and a
warning is raised, since the usual cause is two layers in different
places — a CRS that could only be stamped, not reprojected. The
attribute `"ties"` records how many features matched more than one
polygon and had the `tie_break` rule decide for them: a list with `n`,
`which` (their row positions in `features_sf`) and `rule`. A large `n`
means the polygon layer overlaps, and per-cell counts built from the
result depend on the rule. The record describes the rows this call
returned and does not survive subsetting: `joined[i, ]` is a plain layer
with no `"ties"` attribute, rather than one reporting the parent's count
against row positions that no longer resolve.

## Details

This is the second step of the package's pipeline: it labels every
observation with the cell it falls in, which is what
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
then aggregates over. Reach for it directly (rather than for
[`sf::st_join()`](https://r-spatial.github.io/sf/reference/st_join.html))
when the join has to be *unambiguous* — it resolves features matching
several polygons by an explicit `tie_break` rule instead of silently
duplicating rows, so the assigned layer keeps one row per input feature
and cell-level counts mean what they say.

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
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)

## Examples

``` r
library(sf)
set.seed(1)
pts <- st_as_sf(
  data.frame(x = runif(20, 0, 100), y = runif(20, 0, 100), val = rnorm(20)),
  coords = c("x", "y"), crs = 32632
)
bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
  c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
))), crs = 32632))
grid <- create_grid_polygons(bnd, target_cells = 9, type = "square")
assigned <- assign_features_to_polygons(pts, grid)
table(assigned$poly_id)
#> 
#> 1 2 3 4 5 6 7 8 9 
#> 1 2 3 3 2 3 1 3 2 
```
