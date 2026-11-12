# Plot a tessellation map with optional boundary, seeds, and features

Builds a layered ggplot2 map of the cells, with optional overlays of a
study boundary, seed points and other features. Requires the ggplot2
package.

## Usage

``` r
plot_tessellation_map(
  tessellation_sf,
  boundary = NULL,
  seeds_sf = NULL,
  features_sf = NULL,
  fill_col = NULL,
  palette = "viridis",
  na_fill = "grey90",
  tile_alpha = 0.9,
  outline_col = "white",
  outline_size = 0.2,
  features_col = "#333333",
  features_size = 0.5,
  seeds_col = "#1f77b4",
  seeds_size = 1.5,
  boundary_col = "#111111",
  boundary_size = 0.6,
  labels = FALSE,
  label_col = NULL,
  label_size = 2.7,
  legend = TRUE,
  legend_title = NULL,
  theme = NULL,
  target_crs = NULL,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  xlim = NULL,
  ylim = NULL,
  expand = TRUE
)
```

## Arguments

- tessellation_sf:

  An sf POLYGON/MULTIPOLYGON layer. Required.

- boundary:

  Optional sf/sfc polygon outline layer.

- seeds_sf:

  Optional sf/sfc point layer of seed locations.

- features_sf:

  Optional sf/sfc layer of additional features.

- fill_col:

  Name of the COLUMN in `tessellation_sf` to map to fill; `NULL` for no
  fill. A numeric column gets a continuous scale, as does a Date or
  POSIXct column (on a date or time axis) and a `units` or `difftime`
  column (drawn as numbers, with the unit in the legend title); anything
  else a discrete one. `fill_col` and `label_col` name columns, while
  `outline_col`, `features_col`, `seeds_col` and `boundary_col` are
  colours.

- palette:

  Viridis palette, passed as `option` to
  [`ggplot2::scale_fill_viridis_c()`](https://ggplot2.tidyverse.org/reference/scale_viridis.html)
  or
  [`ggplot2::scale_fill_viridis_d()`](https://ggplot2.tidyverse.org/reference/scale_viridis.html):
  `"viridis"` (default), `"magma"`, `"inferno"`, `"plasma"`,
  `"cividis"`, `"rocket"`, `"mako"` or `"turbo"`.

- na_fill:

  Fill for NA values. Default "grey90".

- tile_alpha:

  Alpha for filled polygons. Default 0.9.

- outline_col, outline_size:

  Colour and line width of the cell outlines. Defaults `"white"` and
  0.2.

- features_col, features_size:

  Colour and point size of the feature overlay. Defaults `"#333333"` and
  0.5.

- seeds_col, seeds_size:

  Colour and point size of the seed overlay. Defaults `"#1f77b4"` and
  1.5.

- boundary_col, boundary_size:

  Colour and line width of the boundary outline. Defaults `"#111111"`
  and 0.6.

- labels:

  Logical; draw per-cell labels. Default FALSE. A cell is labelled when
  its label is not NA and it has a location; when no cell has both, the
  labels are skipped with a log line.

- label_col:

  Name of the COLUMN holding the label text. Default `NULL`: the first
  of `"poly_id"`, `"polygon_id"`, `"id"`, `"cell_id"` and `"grid_id"`
  the layer has, the order in which
  [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
  reads a polygon layer's IDs, so the cells of
  [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
  and the output of
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  are labelled without naming one, by the ID the summary is keyed on. A
  `units` or `difftime` column is drawn formatted, with its unit.

- label_size:

  Label text size. Default 2.7.

- legend:

  Logical; show fill legend. Default TRUE.

- legend_title:

  Optional legend title.

- theme:

  A ggplot2 theme. Default `NULL`, which uses
  [`ggplot2::theme_void()`](https://ggplot2.tidyverse.org/reference/ggtheme.html).

- target_crs:

  Optional CRS for plotting. When `NULL` (the default) the
  tessellation's own CRS is used, or (if the tessellation has none) a
  CRS borrowed from the first overlay that carries one, so a CRS-less
  grid drawn with located features still lines up. Any layer that
  arrives without a CRS is brought into the plot's CRS instead of being
  dropped: reprojected when its coordinates look like
  longitude/latitude, otherwise stamped with a warning, since a stamp
  assumes the coordinates were already in that CRS. A plot CRS in
  longitude and latitude is not stamped on coordinates that cannot be in
  it (an x outside -180 to 360, or a y beyond 90 north or south): the
  call stops and asks for the layer's CRS or a projected `target_crs`.

- title, subtitle, caption:

  Plot annotations.

- xlim, ylim:

  Optional numeric vectors of length 2 for coordinate limits (in the
  plot CRS). Default NULL (auto).

- expand:

  Logical; expand plot area slightly beyond data limits. Default TRUE.

## Value

A ggplot2 object.

## Details

Every layer is drawn in x and y: a layer whose geometries mix XY with
XYZ or XYM (as [`rbind()`](https://rdrr.io/r/base/cbind.html) of an XY
and an XYZ layer gives) is drawn without its Z and M. The labels are
placed in x and y, also on cells whose geometries carry M (XYM or XYZM).

A row of any layer with a coordinate that is not finite or lies beyond
1e150 has no location: it is left out of the map, with one warning that
says how many rows of which argument were not drawn. A row the plot's
CRS cannot place in full is drawn as the move leaves it (empty, or
without the vertices it cannot place), with a warning of class
`spatialkit_rows_unplaced` that names it.

## See also

Other cell building:
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md),
[`clip_target_for()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clip_target_for.md),
[`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md),
[`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md),
[`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md),
[`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md),
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
[`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md),
[`voronoi_seeds_random()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_random.md)

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  n <- 60
  pts <- st_as_sf(
    data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
               z = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
  # Cells from the points, the response aggregated onto them, and the map
  # shaded by the cell mean rather than by an arbitrary ID.  Each cell
  # holds one point, so there is nothing to pool: estimate = "mean".
  tess  <- build_tessellation(pts, method = "voronoi", quiet = TRUE)
  asg   <- assign_features_to_polygons(pts, tess$cells,
                                       polygon_id_col = "cell_id")
  cells <- summarize_by_cell(asg, response_var = "z", id_col = "cell_id",
                             cells_sf = tess$cells, estimate = "mean")
  plot_tessellation_map(cells, features_sf = pts, fill_col = "resp_mean_z",
                        legend_title = "mean z")
}
```
