# Create deterministic, stable polygon IDs based on spatial sort keys

Ensures that a polygon layer has a reproducible, deterministic
identifier column by sorting features using representative point
coordinates (and secondary tie-breakers) and then assigning sequential
IDs.

## Usage

``` r
ensure_stable_poly_id(
  polygons_sf,
  id_col = "poly_id",
  method = c("centroid", "surface_point", "bbox_center"),
  make_valid = TRUE,
  transform_for_sort = 4326
)
```

## Arguments

- polygons_sf:

  An sf or sfc object containing polygonal features.

- id_col:

  Character scalar; name of the identifier column.

- method:

  One of "centroid", "surface_point", "bbox_center".

- make_valid:

  Logical; apply st_make_valid() first. Default TRUE.

- transform_for_sort:

  CRS used only for computing sort-key coordinates. Default 4326. This
  is the whole mechanism by which the IDs are stable (sorting in one
  common CRS is what makes the same layer get the same IDs whichever
  projection it arrives in), so if the transform fails the function says
  so rather than quietly sorting in the input's own CRS. The sort key is
  rounded to 7 decimal degrees (about 1 cm) before ordering, so the
  floating-point noise of a round trip through a different projection
  cannot reverse two neighbouring cells. Set to NULL to sort in the
  input CRS, which gives IDs that are reproducible but not comparable
  across projections.

## Value

An sf polygon layer re-ordered with sequential IDs in id_col.
Non-polygonal rows are **dropped** (with a warning), so the result can
have fewer rows than the input; if no polygonal rows remain, an error is
raised.

## See also

Other tessellation:
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md),
[`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md),
[`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md),
[`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md),
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
[`plot_tessellation_map()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_tessellation_map.md),
[`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md),
[`voronoi_seeds_random()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_random.md)

## Examples

``` r
library(sf)
bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
  c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
))), crs = 32632))
g <- create_grid_polygons(bnd, target_cells = 9)
# Reverse the rows: the IDs come back in the same spatial order regardless
ids_fwd <- ensure_stable_poly_id(g)$poly_id
ids_rev <- ensure_stable_poly_id(g[nrow(g):1, ])$poly_id
identical(sort(ids_fwd), sort(ids_rev))
#> [1] TRUE
```
