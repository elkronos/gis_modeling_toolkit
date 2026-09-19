# K-means seed generation from point coordinates

Places `k` seed points at k-means cluster centres of the observed
coordinates, so seeds (and the Voronoi cells built from them) follow the
sampling density: clusters of observations attract seeds, empty ground
gets none. Reach for this when you want cells that each carry a
comparable number of observations, which is what makes per-cell
aggregates in
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
similarly precise. Use
[`voronoi_seeds_random()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_random.md)
instead when you want coverage of the study area rather than of the
data, and
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
to pick between them by name.

## Usage

``` r
voronoi_seeds_kmeans(points_sf, k, set_seed = 456)
```

## Arguments

- points_sf:

  An sf object with POINT geometries.

- k:

  Integer; requested number of clusters, treated as an upper bound only.
  It is clamped, with a warning, to whichever is smaller of the number
  of distinct point positions and `nrow(points_sf) - 1`, because k-means
  can produce neither more centres than there are distinct points nor as
  many centres as there are rows. Check
  [`nrow()`](https://rdrr.io/r/base/nrow.html) on the result.

- set_seed:

  Optional integer RNG seed. Default 456.

## Value

An sf object of **at most** `k` cluster-centre POINTs (fewer when `k`
exceeds the number of distinct positions), with `seed_id` and
`method = "kmeans"` columns matching
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md).

## Details

Lon/lat input is projected first so the k-means distances are metric and
not degrees. Rows with empty or non-finite coordinates are dropped with
a warning, and `k` is clamped to the number of distinct positions.

## Examples

``` r
library(sf)
set.seed(1)
pts <- st_as_sf(
  data.frame(x = 5e5 + runif(100, 0, 1000), y = 5e6 + runif(100, 0, 1000)),
  coords = c("x", "y"), crs = 32632
)
seeds <- voronoi_seeds_kmeans(pts, k = 8)
nrow(seeds)   # at most 8
#> [1] 8
```
