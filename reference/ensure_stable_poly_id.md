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

  Logical; apply st_make_valid() first. Default TRUE. The copy the sort
  key is measured on is repaired as well, after it has been transformed
  (see Details); with `FALSE` neither is.

- transform_for_sort:

  CRS used only for computing sort-key coordinates. Default 4326. This
  is the whole mechanism by which the IDs are stable (sorting in one
  common CRS is what makes the same layer get the same IDs whichever
  projection it arrives in), so if the transform fails the function says
  so rather than quietly sorting in the input's own CRS. The sort key is
  rounded to 7 decimal degrees (about 1 cm) before ordering, so the
  floating-point noise of a round trip through a different projection
  does not usually reverse two neighbouring cells. It can where two
  cells' centres lie within about that step of the same longitude, as
  fine cells stacked north-south near a projection's central meridian
  do: 36 of 2,500 100 m cells straddling a UTM central meridian changed
  ID after a transform to EPSG:3035. No rounding step removes that, so
  to match cells computed in different projections, join them on
  geometry rather than on the ID. The key is computed on the sphere (s2)
  whether or not
  [`sf::sf_use_s2()`](https://r-spatial.github.io/sf/reference/s2.html)
  is on, so the session setting does not change the IDs. Set to NULL to
  sort in the input CRS, which gives IDs that are reproducible but not
  comparable across projections.

## Value

An sf polygon layer re-ordered with sequential IDs in id_col.
Non-polygonal rows are **dropped** (with a warning), so the result can
have fewer rows than the input; if no polygonal rows remain, an error is
raised.

## Details

A feature that is valid in its own CRS is not always one s2 can measure.
Only the vertices are transformed to the sort CRS, and on the sphere
they are joined by great-circle arcs, so a long edge that is straight in
the layer's CRS and passes about a metre from another vertex of the same
ring can end up on the other side of that vertex. The ring then crosses
itself, and s2 refuses it (1 of 291 Voronoi cells of Texas clipped to
the state outline in EPSG:5070; "Loop 1 is not valid: Edge 36 crosses
edge 52"). The repair of the sort copy does not split crossing edges, so
a feature s2 still refuses is handled on its own. Its sort copy is
repaired a second time with the crossing edges split. If s2 refuses that
as well, or with `make_valid = FALSE`, which asks for no repair, the
feature's centroid and area are measured as plane geometry in the
layer's own CRS and the centroid is transformed to the sort CRS. With
the methods `"surface_point"` and `"bbox_center"` only the area needs
this, because those points are not taken with s2.

The function warns once, saying how many features took each route. Their
keys are close to the spherical ones and not equal to them. The second
repair moved the Texas cell's area by about 40 square metres in 2,949
square kilometres. Its centroid taken in the plane lay 20 m from the one
the second repair gave on the sphere, and the two kinds of centroid lay
5 m apart at the median and 61 m at most for the other 290 cells. A
feature measured in the plane can therefore sort on the other side of a
neighbour whose centre is that close to its own in longitude. Every
other feature's key is unaffected, and the geometry returned is never
the sort copy.

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
# Reverse the rows and re-derive the IDs: the SAME cell gets the same ID,
# which is the property a row-position ID does not have.  The check joins
# the two layers on geometry, because comparing sorted ID vectors would
# pass for any two permutations of 1:9.
fwd <- ensure_stable_poly_id(g)
rev <- ensure_stable_poly_id(g[nrow(g):1, ])
same_cell <- match(st_as_text(st_geometry(rev)), st_as_text(st_geometry(fwd)))
all(rev$poly_id == fwd$poly_id[same_cell])
#> [1] TRUE
```
