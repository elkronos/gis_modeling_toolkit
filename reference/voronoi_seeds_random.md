# Random seed generation within a polygonal boundary

Draws `k` seed points uniformly at random inside `boundary`, ignoring
where the observations are. Reach for this when the cells should cover
the study area evenly — so that sparsely sampled ground still gets its
own cells and is visibly under-sampled in the results — rather than
concentrating resolution where the data already are, which is what
[`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)
does. It is also the honest choice for a null or sensitivity comparison:
re-running an analysis over several random seedings shows how much of a
result depends on one particular tessellation.

## Usage

``` r
voronoi_seeds_random(boundary, k, set_seed = 456)
```

## Arguments

- boundary:

  An sf or sfc polygonal object.

- k:

  Integer; number of random seeds.

- set_seed:

  Integer RNG seed. Default 456.

## Value

An sf object of **at most** `k` random POINTs (rejection sampling inside
an awkward geometry can fall short of `k`, which is warned about), with
`seed_id` and `method = "random"` columns matching
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md).

## Details

Sampling is by rejection inside the polygon, so an awkward geometry can
return fewer than `k` seeds; that shortfall is warned about rather than
silently padded.

## Examples

``` r
library(sf)
bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
  c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
))), crs = 32632))
seeds <- voronoi_seeds_random(bnd, k = 10)
nrow(seeds)   # at most 10
#> [1] 10
```
