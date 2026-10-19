# How far the held-out points sit from the training data

Blocked cross-validation exists to put distance between a test point and
the training points that could tell you its value. Whether it succeeded
is a measurement, and
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
only offers a proxy for it: the block size compared against the
estimated autocorrelation range, reported as a warning. That comparison
is about the *design*. This is about the *result*: for every held-out
point, the distance to its nearest training point, summarised per fold.

## Usage

``` r
fold_separation(folds, data_sf, sac = NULL, prediction_points = NULL)
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
  and every `cv_*()` do. As in `cv_*()`, a
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  result whose recorded rows sit at other locations here (folds built on
  another layer, such as the points before
  [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
  dropped some) is refused; the location check is skipped when one of
  the two layers is POINT and the other is not.

- sac:

  Optional: an
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  result or a single number. Defaults to the range the folds carry, if
  any. Supplying one adds the `mean_cor` column and, with a prediction
  target, the gap on the correlation scale. A range that records its CRS
  (the folds' own, or an
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  result) is compared with distances measured in that CRS, whatever CRS
  `data_sf` is in. A bare number is taken to be in the units the
  distances are otherwise measured in: those of `data_sf` if it is
  projected, and for geographic (lon/lat) input metres, in the CRS
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  chooses (as for
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)'s
  `block_size`), not degrees. A `units` object is refused.

- prediction_points:

  Optional `sf` or `sfc` layer of the locations the model will predict
  at: the grid from
  [`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md),
  or a grid over the region of interest
  (`sf::st_make_grid(region, what = "centers")`). A non-POINT layer is
  reduced to one representative point per feature, as in
  `make_folds(method = "nndm")`, so pass grid cells or points, not one
  polygon of the region. It is brought into the CRS the distances are
  measured in. Default `NULL`: no target and no reading, except for NNDM
  and kNNDM folds (`method = "nndm"` and `"knndm"`), which carry the
  distances from the prediction points they were built for and are
  compared against those, measured to all sites. The carried distances
  are lengths in the CRS the folds were built in (`folds$params$crs`).
  When the distances here are measured in another CRS (`data_sf` is in
  another projected CRS, or `sac` records another), the two are not
  compared, [`print()`](https://rdrr.io/r/base/print.html) says so, and
  passing `prediction_points` has the target measured here.

## Value

A data.frame of class `fold_separation`, one row per fold: `fold` (the
fold's number: for the `$folds` of a `cv_*()` result, the `fold_id` its
`fold_metrics` use, which differs from the list position once a fold has
been dropped), `n_train`, `n_test`, `n_blocks` (`NA` for a scheme with
no blocks), `min_dist` and `median_dist` (distance from a held-out point
to its nearest training point, in the units of the CRS the `crs`
attribute names), `mean_cor`, `target_median_dist` (the median distance
from the prediction points to the fold's training data) and `dist_gap`
(the fold's mean held-out distance minus its mean prediction distance);
the last two are `NA` without a target, and
[`print()`](https://rdrr.io/r/base/print.html) shows them under the
headings `target` and `gap`.

`mean_cor` is a mean correlation, `NA` without a range: the mean over
the fold's held-out points of \\\rho(d)\\, the correlation the range
leaves between a point and its nearest training point at distance \\d\\.
It is 1 for a point on top of training data and falls towards 0 beyond
the range, so a point at a twentieth of the range counts for more than
one at nine tenths of it. \\\rho(d)\\ is \\\exp(-3 d / range)\\, the
exponential model whose effective range is the range given (0.05 at the
range), or the spherical model when the range comes from an
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
fit that was spherical; the `range_model` attribute says which. Earlier
development versions had a column `within_range` in its place: the share
of held-out points closer than the range, which counted every such point
the same.

Attributes: `method`, `sac_range`, `range_model`, `crs` (the CRS the
distances were measured in), `n_unknown_ids` (the number of distinct
rows the folds name that `data_sf` does not have, each counted once
however many folds mention it, as in the `cv_*()` field of the same
name), `pooled` (a list with the number, minimum, median and mean of the
held-out distances over all folds, and the pooled `mean_cor`) and
`target`: `NULL` without a prediction target, otherwise a list with
`n_prediction_points`, `source`, `basis` (`"each fold's training data"`
or `"all sites"`), `held_out_mean`, `held_out_median`, `target_mean`,
`target_median`, `target_all_mean`, `target_all_median` (the prediction
distances to all sites), `mean_diff`, `w1`, `cor_gap`, `cor_w1` (`NA`
without a range), `scale` (`"correlation"` or `"distance"`),
`tolerance`, `reading` (see **Reading the gap**) and `nndm_out_of_reach`
(for NNDM and kNNDM folds, the `params$target_out_of_reach` that
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
recorded; `NA` for other schemes). `not_compared` is `NULL`, or, for
NNDM or kNNDM folds whose carried target was not compared because the
distances here are measured in another CRS, one sentence naming the two
CRSs; `target` is then `NULL` and
[`print()`](https://rdrr.io/r/base/print.html) gives the sentence where
it would otherwise say that no target is given.

## Details

The distances do not say by themselves which way a score errs. A
cross-validated error reads low when the held-out points sit closer to
their training data than the locations the model will predict at, and
high when they sit farther, so the same held-out distances can go with
either, depending on the target. Random folds on evenly sited data,
scored for a map of the sampled area, are close to unbiased although
nearly every held-out point has a training neighbour inside the range;
range-sized blocks on the same data read high; random folds on clustered
data scored for a map of the whole region read low. So without a
prediction target the function reports the distances and draws no
conclusion. Give it `prediction_points` and it sets the held-out
distances beside the distances from those points to the same training
data, and words its reading by the signed gap between the two.

Nothing is estimated: the distances come from the geometry, and the
range, when one is shown, is the one the folds already carry
(`make_folds( auto_range = TRUE)` records it) or the one you pass as
`sac`.

## Reading the gap

With a target, the held-out distances and the prediction distances are
compared as two distributions. Each fold's prediction distances are
measured to that fold's own training rows, so that a five-fold design is
not read as too far only because each fold trains on four fifths of the
sites; the distances to all sites are reported beside them. (With more
than 50 folds, the leave-one-out schemes, every training set is nearly
the whole sample and the prediction points are measured against all
sites once.) Folds count by their share of the held-out points.

Two numbers summarise the comparison: the signed difference of the
means, held-out minus prediction (the signed area between the two
cumulative distributions), and the Wasserstein distance \\W_1\\ (the
unsigned area). When a range is available both are also given on the
correlation scale, after each distance \\d\\ is replaced by \\1 -
\rho(d)\\, with \\\rho\\ the correlation the range implies (see
`mean_cor` under Value). That scale counts a difference between two
distances well inside the range and ignores one between two distances
beyond it.

The reading is taken on the correlation scale when there is a range,
with a tolerance of 0.1, and otherwise on the mean distance, with a
tolerance of a quarter of the prediction mean. (When the prediction
points are measured to all sites and every one lies on a site, as when
the sites themselves are passed as `prediction_points` to a
leave-one-out design, the prediction mean and with it that tolerance are
0: held-out distances above 0 then read `"farther"`, and held-out
distances that are all 0 as well read `"matched"`. With 50 folds or
fewer the points are measured to each fold's own training rows, which
leave out that fold's held-out sites, so the mean is above 0 and the
quarter rule gives the reading.)

- `"closer"`:

  The gap is below minus the tolerance and at least half of \\W_1\\ in
  size: the held-out points sit closer to their training data than the
  prediction points do, and the cross-validated error tends to read low
  for that target.

- `"farther"`:

  The gap is above the tolerance and at least half of \\W_1\\: the
  held-out points sit farther, and the error tends to read high.

- `"matched"`:

  The gap is within the tolerance and \\W_1\\ within twice the
  tolerance.

- `"mixed"`:

  The two distributions cross, closer at some distances and farther at
  others (the gap is under half of \\W_1\\, or it is within the
  tolerance while \\W_1\\ is not), so the distances give no single
  direction.

The reading is a statement about nearest-neighbour distances, not a
measured bias. A learner that draws on several neighbours can still
score low on clustered samples whose nearest-neighbour distances match.

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
  data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000)),
  coords = c("x", "y"), crs = 32632
)

# Random folds put a training point close to every held-out one.
random  <- make_folds(pts, k = 4, method = "random_kfold", seed = 1)
print(fold_separation(random, pts, sac = 200))
#> Fold separation: random_kfold, 4 fold(s), 200 held-out point(s) (EPSG:32632)
#>   autocorrelation range: 200 (in metres, like the distances)
#> 
#>  fold n_train n_test min_dist median_dist mean_cor
#>     1     150     50    4.696       36.95     0.57
#>     2     150     50    3.835       42.14     0.55
#>     3     150     50    3.835       38.52     0.58
#>     4     150     50    4.696       45.79     0.53
#> 
#>   Held-out points sit a median of 41.61 from their nearest training
#>   point; the closest is 3.835 away. Under the range of 200 that leaves a
#>   mean correlation of 0.56 between a held-out point and its nearest
#>   training point (exp(-3 d / range) averaged over the hold-out: 1 at
#>   distance 0, 0.05 at the range).
#>   No prediction target is given, so these are distances and no conclusion
#>   about the score: a cross-validated error reads low when the hold-out
#>   sits closer to its training data than the locations to be predicted do,
#>   and high when it sits farther. Pass `prediction_points` to compare the
#>   two.

# Blocked folds hold out whole neighbourhoods, so the distances grow.
blocked <- make_folds(pts, k = 4, method = "block_kfold",
                      block_size = 250, seed = 1)
print(fold_separation(blocked, pts, sac = 200))
#> Fold separation: block_kfold, 4 fold(s), 200 held-out point(s) (EPSG:32632)
#>   autocorrelation range: 200 (in metres, like the distances)
#> 
#>  fold n_train n_test n_blocks min_dist median_dist mean_cor
#>     1     152     48        2    23.02       111.9     0.23
#>     2     143     57        3    29.04       169.5     0.14
#>     3     152     48        2    21.03       118.3     0.22
#>     4     153     47        2    21.03       115.9     0.23
#> 
#>   Held-out points sit a median of 123.5 from their nearest training
#>   point; the closest is 21.03 away. Under the range of 200 that leaves a
#>   mean correlation of 0.20 between a held-out point and its nearest
#>   training point (exp(-3 d / range) averaged over the hold-out: 1 at
#>   distance 0, 0.05 at the range).
#>   No prediction target is given, so these are distances and no conclusion
#>   about the score: a cross-validated error reads low when the hold-out
#>   sits closer to its training data than the locations to be predicted do,
#>   and high when it sits farther. Pass `prediction_points` to compare the
#>   two.

# Which of the two suits a map of the sampled square?  Compare each with
# the distances from the cells of that map to the training data.
cells <- st_make_grid(pts, n = c(20, 20), what = "centers")
print(fold_separation(random, pts, sac = 200, prediction_points = cells))
#> Fold separation: random_kfold, 4 fold(s), 200 held-out point(s) (EPSG:32632)
#>   autocorrelation range: 200 (in metres, like the distances)
#>   prediction target: 400 point(s), measured to each fold's training data
#> 
#>  fold n_train n_test min_dist median_dist target     gap mean_cor
#>     1     150     50    4.696       36.95  37.82 -0.9056     0.57
#>     2     150     50    3.835       42.14  38.02  0.2229     0.55
#>     3     150     50    3.835       38.52  36.00 -1.8310     0.58
#>     4     150     50    4.696       45.79  37.28  3.9360     0.53
#> 
#>   Held-out points sit a median of 41.61 from their nearest training
#>   point; the closest is 3.835 away. Under the range of 200 that leaves a
#>   mean correlation of 0.56 between a held-out point and its nearest
#>   training point (exp(-3 d / range) averaged over the hold-out: 1 at
#>   distance 0, 0.05 at the range).
#>   The prediction points sit a median of 37.48 from the same training data
#>   (32.16 from the nearest of all sites). Held-out minus prediction:
#>   +0.3554 in mean distance (W1 3.371); +0.01 on the correlation scale (W1
#>   0.02). The held-out points sit about as far from their training data as
#>   the prediction points do. Matched nearest-neighbour distances do not
#>   make a score unbiased by themselves: a learner that draws on several
#>   neighbours can still read low on clustered samples.
print(fold_separation(blocked, pts, sac = 200, prediction_points = cells))
#> Fold separation: block_kfold, 4 fold(s), 200 held-out point(s) (EPSG:32632)
#>   autocorrelation range: 200 (in metres, like the distances)
#>   prediction target: 400 point(s), measured to each fold's training data
#> 
#>  fold n_train n_test n_blocks min_dist median_dist target   gap mean_cor
#>     1     152     48        2    23.02       111.9  39.37 58.97     0.23
#>     2     143     57        3    29.04       169.5  44.55 89.14     0.14
#>     3     152     48        2    21.03       118.3  40.95 65.42     0.22
#>     4     153     47        2    21.03       115.9  38.61 67.19     0.23
#> 
#>   Held-out points sit a median of 123.5 from their nearest training
#>   point; the closest is 21.03 away. Under the range of 200 that leaves a
#>   mean correlation of 0.20 between a held-out point and its nearest
#>   training point (exp(-3 d / range) averaged over the hold-out: 1 at
#>   distance 0, 0.05 at the range).
#>   The prediction points sit a median of 40.55 from the same training data
#>   (32.16 from the nearest of all sites). Held-out minus prediction:
#>   +71.05 in mean distance (W1 71.15); +0.30 on the correlation scale (W1
#>   0.30). The held-out points sit farther from their training data than
#>   the prediction points do, so the cross-validated error tends to read
#>   high for these prediction points.
```
