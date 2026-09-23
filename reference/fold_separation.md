# How far the held-out points actually sit from the training data

Blocked cross-validation exists to put distance between a test point and
the training points that could tell you its value. Whether it succeeded
is a measurement, and until now the package only offered a proxy for it:
the block size compared against the estimated autocorrelation range,
reported by
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
as a warning. That comparison is about the *design*. This is about the
*result*: for every held-out point, the distance to its nearest training
point, summarised per fold.

## Usage

``` r
fold_separation(folds, data_sf, sac = NULL)
```

## Arguments

- folds:

  A
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  result, or its `$folds` element (a list of `train`/`test` splits).

- data_sf:

  The layer the folds were built on. Row identifiers are matched through
  `..row_id` when the layer carries one, and by row position otherwise,
  which is what
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  and every `cv_*()` do.

- sac:

  Optional: an
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  result or a single number, in the CRS units of `data_sf`. Defaults to
  the range the folds carry, if any. Supplying one adds the
  `within_range` column and the closing verdict.

## Value

A data.frame of class `fold_separation`, one row per fold: `fold`,
`n_train`, `n_test`, `n_blocks` (`NA` for a scheme with no blocks),
`min_dist` and `median_dist` (distance from a held-out point to its
nearest training point, in CRS units), and `within_range` (the share of
held-out points closer to training data than `sac`; `NA` without one).
Attributes: `method`, `sac_range`, `crs` and `n_unknown_ids`.

## Details

The two can disagree, and the direction is not obvious. Blocks wider
than the range still leak wherever a test point sits near a block edge
with training data just across it, which is most of the points in a fine
block grid; conversely a fold whose blocks are narrower than the range
can still separate well if the points inside them are clustered. The
share of held-out points closer to training data than the correlation
range is the number that settles it, and it is the last column here.

Nothing is estimated: the distances come from the geometry, and the
range, when one is shown, is the one the folds already carry
(`make_folds( auto_range = TRUE)` records it) or the one you pass as
`sac`.

## See also

[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
for the fold schemes and the block sizing this measures the outcome of;
[`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md)
for choosing a block size by cross-validated error instead.

Other cross-validation:
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md),
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
[`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md),
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md),
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
[`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md),
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md),
[`sac_nugget()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md),
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)

## Examples

``` r
library(sf)
set.seed(1)
n <- 200
pts <- st_as_sf(
  data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000)),
  coords = c("x", "y"), crs = 32632
)

# Random folds put a training point almost on top of every held-out one.
random  <- make_folds(pts, k = 4, method = "random_kfold", seed = 1)
print(fold_separation(random, pts, sac = 200))
#> Fold separation: random_kfold, 4 fold(s), 200 held-out point(s) (EPSG:32632)
#>   autocorrelation range: 200
#> 
#>  fold n_train n_test min_dist median_dist within_range
#>     1     150     50    4.696       36.95         100%
#>     2     150     50    3.835       42.14         100%
#>     3     150     50    3.835       38.52         100%
#>     4     150     50    4.696       45.79         100%
#> 
#>   100% of held-out points sit closer to a training point than the
#>   correlation range (200), and the closest is 3.83 away. Most of the
#>   hold-out is inside the range of its own training data, so this score is
#>   optimistic: widen the blocks.

# Blocked folds hold out whole neighbourhoods, so the distances grow.
blocked <- make_folds(pts, k = 4, method = "block_kfold",
                      block_size = 250, seed = 1)
print(fold_separation(blocked, pts, sac = 200))
#> Fold separation: block_kfold, 4 fold(s), 200 held-out point(s) (EPSG:32632)
#>   autocorrelation range: 200
#> 
#>  fold n_train n_test n_blocks min_dist median_dist within_range
#>     1     152     48        2    23.02       111.9         100%
#>     2     143     57        3    29.04       169.5          67%
#>     3     152     48        2    21.03       118.3          96%
#>     4     153     47        2    21.03       115.9          85%
#> 
#>   86% of held-out points sit closer to a training point than the
#>   correlation range (200), and the closest is 21 away. Most of the
#>   hold-out is inside the range of its own training data, so this score is
#>   optimistic: widen the blocks.
```
