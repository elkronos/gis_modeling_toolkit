# Area of applicability of a spatial prediction model

Computes the dissimilarity index (DI) of Meyer & Pebesma (2021) for a
set of prediction locations and flags those that fall inside the model's
area of applicability (AOA), the region of predictor space where the
model's cross-validated performance estimate can be expected to hold.

## Usage

``` r
area_of_applicability(
  newdata,
  model = NULL,
  train_sf = NULL,
  predictor_vars = NULL,
  weights = "auto",
  folds = NULL,
  threshold = NULL,
  normalizer_max_n = 5000L,
  seed = 123L,
  chunk_size = NULL,
  use_fnn = requireNamespace("FNN", quietly = TRUE),
  coords = c("auto", "yes", "no"),
  categorical = c("refuse", "indicator", "drop")
)
```

## Arguments

- newdata:

  Prediction locations: an `sf` object (typically from
  [`predict_surface`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md))
  or a data.frame, carrying the predictor columns. When location is
  measured (a `gwr_fit`, a `bayesian_fit`, a forest fitted with
  `include_coords = TRUE`, or `coords = "yes"`) an `sf` object is
  required, non-`POINT` geometry is reduced to representative points,
  and a CRS mismatch with the training data is reconciled; see *Models
  whose predictions depend on location*.

- model:

  A fitted `spatial_fit`, supplying the training data and predictor
  names, the importance `weights = "auto"` reads, and whether its
  predictions depend on location. Optional if `train_sf` and
  `predictor_vars` are given directly, which lets this be used with any
  model.

- train_sf:

  Training data, if not taken from `model`.

- predictor_vars:

  Predictor names, if not taken from `model`.

- weights:

  `"auto"` (default), `NULL`, or a named numeric vector of predictor
  importances. `"auto"` uses `pmax(model$info$importance, 0)` when
  `model` carries a finite importance with at least one positive value,
  naming every predictor in use, from a model with out-of-bag skill;
  otherwise every predictor counts equally and a logged note says why
  (see *Weights*). `NULL` is the explicit request for equal weights. For
  a vector of your own, any positive scale works. When `model` was
  fitted with `include_coords = TRUE`, the coordinates `"..x"` and
  `"..y"` are measured as predictors too (see below), and you are not
  expected to supply importances for them: any you leave out default to
  the mean of the weights you did supply, so location counts about as
  much as a typical predictor. Naming them explicitly overrides that. An
  unnamed vector may have one value per predictor either with or without
  the two coordinate columns. Weights must be finite and non-negative,
  so pass permutation importance as `pmax(importance, 0)`. (A forest
  with no out-of-bag rows, `replace = FALSE` with `sample_fraction = 1`,
  has `NaN` importance, which
  [`pmax()`](https://rdrr.io/r/base/Extremes.html) keeps and which is
  refused when passed by hand; `"auto"` uses equal weights for it.)
  `pmax(importance, 0)` is all zero when the model found no predictor
  useful, and then the weights cannot say anything: with a single
  predictor any weight gives the same index and zero is accepted; with
  several, all of them are weighted equally, as with `weights = NULL`,
  and a warning says so. The coordinate default above is the mean of the
  supplied weights, so a zero weight on the only covariate of a
  coordinate-using model zeroes the coordinates too, and that equal
  weighting applies.

- folds:

  Cross-validation folds: a
  [`make_folds`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  result, a list of `train`/`test` splits, or a vector of fold labels
  with one entry per training row. Default `NULL` (plain nearest
  neighbour). The folds you passed to `cv_*()`, built on the layer
  `model` was fitted from, may name rows that
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
  removed (a missing or non-finite value, an empty geometry): as in
  `cv_*()`, they are dropped from the folds, and a label vector with one
  entry per row of that layer loses theirs. Also as in `cv_*()`, a
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  result built on other data – another layer, or these rows in another
  order – is refused. That check is skipped when the folds were built on
  polygons and the training data are the points a fit reduced them to.

- threshold:

  Optional numeric override for the DI threshold of the predictor index.
  The location index, when there is one, always takes its threshold from
  the rule.

- normalizer_max_n:

  Subsample the training data to this many points when computing the
  mean pairwise distance, which is quadratic. Default 5000.

- seed:

  Seed for that subsample. Default 123.

- chunk_size:

  Query rows per distance block on the dense path. Default `NULL`
  (chosen from the training size). Otherwise a single number of at least
  1; a fractional value is truncated to a whole number of rows.

- use_fnn:

  Use FNN for nearest-neighbour search when available. Exposed so the
  dense fallback can be tested. The two paths return the same distances:
  the dense one searches with a fast formula and then recomputes the
  nearest distance exactly, so an exact copy of a training row is at
  distance 0 on both.

- coords:

  Whether location is measured: `"auto"` (default) measures it for a
  model whose predictions depend on it (a `gwr_fit`, a `bayesian_fit`, a
  forest fitted with `include_coords = TRUE`) and not otherwise; `"yes"`
  measures it whatever the model, and with no `model`; `"no"` never
  does, and logs a caution for a model that depends on it. See *Models
  whose predictions depend on location*.

- categorical:

  What to do with factor and character predictors: `"refuse"` (default)
  stops with a message naming them and the two alternatives;
  `"indicator"` codes each level as a 0/1 column and marks a prediction
  row on a level no training row has as outside (`DI = Inf`); `"drop"`
  computes the index on the other predictors and logs which were left
  out. See *Limitations* for the scaling `"indicator"` uses.

## Value

An object of class `aoa`: a list with

- `aoa`: `newdata` with a numeric `DI` column and a logical `AOA` column
  added, and between them a numeric `DI_location` column when location
  has an index of its own. For a forest fitted with
  `include_coords = TRUE` it is the object the computation ran on:
  `newdata` after pointizing, CRS reconciliation and the addition of the
  `"..x"` and `"..y"` columns. With a location index it is `newdata` as
  passed. A row whose predictors are not all finite gets `NA` in `DI`,
  and a row outside the training range of a predictor in `dropped_vars`,
  or on a level no training row has, gets `DI = Inf` and `AOA = FALSE`
  (see *Limitations*). `AOA` is `DI <= threshold`, and with a location
  index also `DI_location <= location$threshold`: a row outside on
  either index is `FALSE`, and a row with one index missing and the
  other inside is `NA`.

- `threshold`: the DI cut-off used.

- `train_DI`: the training points' own DI values.

- `normalizer`: the mean pairwise training distance.

- `weights`: the weight vector actually applied, named by
  `predictor_vars`.

- `predictor_vars`: the predictors used, including `"..x"`/`"..y"` when
  a forest uses the coordinates and excluding `dropped_vars` and, under
  `categorical = "drop"`, `categorical_vars`.

- `dropped_vars`: predictors dropped for negligible variance.

- `scaling`: a list with `center` and `scale`, each named by
  `predictor_vars`: the training means and standard deviations the index
  is computed in, so a location's DI can be traced to the predictor that
  put it outside. A predictor coded as indicators has one entry per
  level, named `"<predictor>=<level>"`.

- `n_outliers`: the number of training DI values above the
  `Q3 + 1.5 * IQR` fence, which the default threshold rule sets aside
  (the "outlier-removed" in its name); computed whether or not
  `threshold` was supplied.

- `n_duplicates`: with `folds = NULL`, the number of training rows that
  repeat an earlier row's predictor values exactly and so share its
  reference distance; `NA` when folds were supplied (see *Training rows
  with identical predictor values*).

- `location`: `NULL`, or for a location index a list with `threshold`,
  `train_DI`, `normalizer` (the mean distance between training points,
  in map units, so `threshold * normalizer` is the cut-off as a
  distance), `n_duplicates` (training rows at an earlier row's location;
  `NA` when folds were supplied) and `n_outside` (prediction rows beyond
  the location threshold).

- `categorical_vars`: the factor and character predictors found, coded
  or left out as `params$categorical` says.

- `n_train`, `n_new`, `n_inside`, `n_outside`, `n_na`: row counts.
  `n_train` counts the training rows that survived the finite-value
  filter; `n_new` is every row of `newdata`, so
  `n_new = n_inside + n_outside + n_na`, with `n_na` the rows whose
  `AOA` is `NA`. `n_outside_predictors` counts the rows beyond the
  predictor threshold, which is `n_outside` unless there is a location
  index.

- `params` records the call: `folds_supplied`, `n_folds`,
  `folds_method`, `threshold_supplied`, `normalizer_max_n`,
  `normalizer_n_used`, `normalizer_subsampled`, `weights_supplied` (the
  caller passed a vector), `weights_source` (`"model importance"`,
  `"supplied"` or `"equal"`), `weights_note` (the kind of importance, or
  why the weights are equal), `coords` (`"location index"`,
  `"predictors"` or `"none"`), `categorical` and `seed`. `folds_method`
  is the `method` of a
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  result (`"block_kfold"`, `"random_kfold"`, ...), `"labels"` for a
  vector of fold labels, `"splits"` for a bare list of train/test
  splits, and `NA` when no folds were supplied. It is printed with the
  object, because the threshold's meaning depends on it: the
  cross-validated DI that sets it comes from the same kind of hold-out
  as the CV error it should be quoted beside, so an AOA built on
  `random_kfold` folds is logged as a caution and does not belong next
  to a blocked `cv_*()` result.

## Why a map alone is not enough

A fitted model will return a number for any location you hand it,
including locations whose predictor values look nothing like anything it
was trained on. Those predictions are extrapolations dressed as
interpolations, and a cross-validation score says nothing about them,
because the held-out folds were drawn from the same predictor
distribution as the training data. The AOA marks where the score
applies.

## How it is computed

Predictors are centred and scaled using the training data's own means
and standard deviations, then weighted by variable importance when the
model carries one (see *Weights*). For a prediction point \\p\\, the DI
is the distance to its nearest training point in that space, divided by
the mean pairwise distance among training points. The same quantity is
computed for the training data itself, using each point's nearest
neighbour *among the training rows of the fold that holds it out*. That
means everything outside its own fold for random and block folds, and
the smaller training set that buffered and NNDM folds actually leave
(see *The fold scheme changes the answer*); without folds it is the
nearest training row with different predictor values. The threshold is
then the largest training DI that is not an upper outlier, i.e. not
above the fence `Q3 + 1.5 * IQR` of the training DI, with the quartiles
of [`stats::quantile()`](https://rdrr.io/r/stats/quantile.html)'s
default type 7. Prediction points at or below that threshold are inside
the AOA.

That is the paper's "outlier-removed maximum". CAST, the reference
implementation, computes the same fence with the same quartiles but uses
the fence itself as the threshold, capped at the largest training DI.
The two agree whenever no training DI lies above the fence (`n_outliers`
is 0); otherwise CAST's threshold is the larger, and so is its AOA.
(Earlier CAST releases used
[`grDevices::boxplot.stats()`](https://rdrr.io/r/grDevices/boxplot.stats.html),
which gives the rule used here but with Tukey's hinges as the quartiles,
so they can also differ when the number of training points is even.) To
apply the current CAST rule to the same training DI, pass
`threshold = min(quantile(res$train_DI, 0.75) + 1.5 * IQR(res$train_DI), max(res$train_DI))`
for an earlier result `res`.

The DI is invariant to the overall scale of `weights`: the numerator and
the normaliser carry the same factor. Importance values can be passed
as-is.

## Training rows with identical predictor values

Repeat visits to a site with static covariates, or covariates read off a
raster coarser than the sampling, give training rows that are exact
copies of each other in predictor space. A copy is the same point of
that space, not a neighbour. Counted as a neighbour it gives a training
DI of 0, and once about three quarters of the rows have one the
threshold is 0 and only exact copies of a training row count as inside
(60 sites visited 4 times put 30 of 200 new points inside; the same
sites counted once put 188).

With `folds = NULL` each training row's reference is therefore its
nearest row with *different* values: the search runs on the distinct
rows and every row takes its point's distance. `n_duplicates` counts the
rows that repeat an earlier one, and a layer with none is computed
exactly as before. The normaliser stays the mean pairwise distance over
all the training rows, so a prediction point's DI is the same with and
without folds.

With `folds` supplied the definition is applied as it stands: a held-out
row whose copy is on the training side of its fold has a DI of 0,
because that hold-out did have an identical row to learn from. If the
threshold comes out 0 that way, a warning of class
`spatialkit_aoa_zero_threshold` gives the count, so a script or a
cross-validation loop can catch it with
`tryCatch(..., spatialkit_aoa_zero_threshold = )`. Folds that keep the
copies together
(`make_folds(method = "leave_location_out", group_var = ...)`) give the
threshold its meaning back. The same applies to the location index when
training rows share a location.

## The fold scheme changes the answer, and should

With `folds = NULL` the training reference is each point's nearest
neighbour anywhere in the training data, which for clustered data is
very close, giving a small threshold and a conservative AOA. Passing the
folds you actually validated with makes the reference distances larger
and the AOA correspondingly wider. That is not a loophole. The AOA is
defined relative to a performance estimate, and a spatially blocked
estimate is a claim about predicting further away. Pass the same
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
result you passed to
[`cv_spatial`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md).
Buffered and NNDM folds use the training set they actually left
available, not merely "everything outside the fold".

## Weights

Meyer and Pebesma weight the predictors by the model's variable
importance, and with the default `weights = "auto"` so does this
function whenever `model` carries an importance that means something:
`pmax(model$info$importance, 0)`, which a forest from
[`fit_rf_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
has. Equal weights are used, with a logged note, when it does not: no
`model`; a model with no importance (`importance = "none"`, a GWR, a
Bayesian fit); an importance that is not finite (`NaN` for a forest with
no out-of-bag rows); one that names fewer predictors than are in use; no
positive value; or a forest with no out-of-bag skill (`oob_r_squared`
missing or at most 0), whose positive importances are noise.
[`print()`](https://rdrr.io/r/base/print.html) says which it was.
`weights = NULL` asks for equal weights whatever the model carries.

What the weights buy: with 2 predictors the forest uses and 8 it does
not, equal weights put about 70\\ (4 times the in-sample RMSE) inside
the AOA and importance weights about 26\\ 43\\ so its square enters the
distance, and a weak but real predictor is down-weighted by the square
of its relative importance. A shift along one can then go unseen (94\\
test, against 19\\ `weights = sqrt(pmax(importance, 0))` is the gentler
alternative: it kept most of the gain on used predictors and put 46\\
shift inside, at the cost of 77\\ shifts.

## Limitations

Categorical predictors (factor or character columns) have no distance of
their own, and by default they are refused and never silently
dummy-coded. `categorical = "indicator"` codes each level that occurs in
the training rows as a 0/1 column. The columns of one predictor share
one scale, `sqrt(1 - sum(p^2))` with `p` the levels' shares of the
training rows, so two rows on different levels are the same distance
apart whichever levels they are, two on the same level are 0 apart, and
the predictor adds as much to the mean squared distance between training
rows as one standardised numeric predictor. That scaling is a choice,
not a result, which is why it is not the default. A prediction row on a
level no training row has gets `DI = Inf`, is outside the AOA, and a
warning gives the count. `categorical = "drop"` leaves the categorical
predictors out and computes the index on the others; a level the
training data never has then goes unseen.

Predictors whose variance is negligible *relative to their own
magnitude* (the test is `sd < sqrt(.Machine$double.eps) * max(abs(x))`,
so the same variable in metres and in gigametres is treated identically)
are dropped from the distance. A prediction point taking a value there
outside the training range (with the same relative tolerance) is
extrapolation along a direction the training data never varied in: its
scaled distance along it is infinite, so it gets `DI = Inf`, is outside
the AOA, and a warning gives the count. A point missing that value is
judged on the other predictors. With equal weights every predictor
counts the same, which overstates dissimilarity along directions the
model barely uses; see *Weights*.

## Models whose predictions depend on location

For three kinds of model the prediction is a function of where the point
is, not only of its covariates, so the index has to measure location
too. Without that a prediction point far outside the training extent but
with ordinary covariate values reads as *inside* the area of
applicability, which is exactly the extrapolation this index exists to
catch. `coords = "auto"` (the default) measures location for all three;
`coords = "no"` turns it off, with a logged caution; `coords = "yes"`
measures it for any model, or for none.

**A geographically weighted regression** (`gwr_fit`) has coefficients
that vary over the map, and **a Bayesian spatial model**
(`bayesian_fit`) has a Gaussian-process term in the coordinates. For
these, and under `coords = "yes"`, location gets an index of its own,
built like the predictor index but in map space: `DI_location` is the
distance to the nearest training point over the mean distance between
training points; the training points' own values use the same folds; and
the threshold is the same outlier-removed maximum. A prediction point is
inside the AOA only if it is inside on both indices, and both DI columns
are returned so the cause can be read off. Location is kept apart from
the covariates because, as two more equally weighted columns, it is
diluted as covariates are added: with 12 of them, 89\\ beyond a
1000-unit training extent still read as inside, against none with the
separate index. In a test with a GWR whose slope rises eastward, points
one half to one and a half extents beyond the data had 13 times the
in-extent RMSE; 84\\ with location measured.

Inside the extent the index marks an unsampled gap only where the gap is
wider than the distance its threshold allows, and that distance depends
on `folds`. With `folds = NULL` it is the spacing of the training points
(the outlier-removed largest nearest-neighbour distance). With blocked
folds it is the largest distance from a held-out point to the nearest
row its fold trained on, outliers removed, which grows with the block
size; a gap narrower than about twice that distance has no cell beyond
it and reads as inside. In a test with 400 training points in two strips
400 units apart, the cut-off without folds was 46 to 54 units and no
grid cell in the middle half of the gap, 100 units or more from either
strip, was inside. That held although the prediction grid took its
covariates from the nearest observation, so that every cell's covariate
DI was 0. With five block folds the cut-off was 199 to 251 units and 76
to 100\\ three draws). This follows from the definition and is not a
fault of the index: the cross-validated error it is paired with was
earned on hold-outs that far from their training rows. It does mean that
with blocked folds the location index says nothing about a gap narrower
than about twice the cut-off; `location$threshold * location$normalizer`
is the cut-off in map units, and
[`print()`](https://rdrr.io/r/base/print.html) shows it.

**A forest fitted with `include_coords = TRUE`** splits on the
coordinates like on any other column, so there they join the predictor
index as `"..x"` and `"..y"` and are centred, scaled and weighted with
the rest.

Either way geometry is needed on both sides, so `train_sf` and `newdata`
must both be `sf` objects; a data.frame is refused rather than quietly
measured without location. Non-`POINT` `newdata` (grid polygons, say) is
reduced to representative points first, as
[`coerce_to_points`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md)
would. If exactly one side carries a CRS the other is brought into it:
reprojected when its coordinates look like longitude/latitude, stamped
otherwise, with a warning either way. This is done because degrees fed
into a metre-space index silently understate the distances. The location
index is a map distance and refuses training data in a geographic CRS.

## References

Meyer, H. and Pebesma, E. (2021). Predicting into unknown space?
Estimating the area of applicability of spatial prediction models.
*Methods in Ecology and Evolution* 12(9), 1620–1633.
[doi:10.1111/2041-210X.13650](https://doi.org/10.1111/2041-210X.13650)

## See also

[`predict_surface`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md)
to build the grid,
[`make_folds`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
for the fold scheme.

Other cross-validation:
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
[`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md),
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md),
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
[`fold_separation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fold_separation.md),
[`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md),
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md),
[`sac_nugget()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md),
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)

Other prediction:
[`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md)

## Examples

``` r
library(sf)
#> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE
set.seed(1)
n <- 120
train <- st_as_sf(
  data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
             a = rnorm(n), b = rnorm(n)),
  coords = c("x", "y"), crs = 32632
)
train$z <- 2 * train$a - train$b + rnorm(n, 0, 0.3)

# Prediction points, some of them well outside the training predictor range
newpts <- st_as_sf(
  data.frame(x = 5e5 + runif(50, 0, 1000), y = 5e6 + runif(50, 0, 1000),
             a = c(rnorm(40), rnorm(10, 8)), b = rnorm(50)),
  coords = c("x", "y"), crs = 32632
)

res <- area_of_applicability(newpts, train_sf = train,
                             predictor_vars = c("a", "b"))
res
#> Area of applicability (Meyer & Pebesma 2021)
#> 
#>   predictors  : 2 (a, b)
#>   weighted    : no (all predictors count equally: no model was given)
#>   training    : 120 points
#>   reference   : nearest other training point (no folds supplied)
#>   normaliser  : 1.7777 (mean pairwise distance)
#>   threshold   : 0.2899 (outlier-removed max of training DI)
#>   location    : not measured
#> 
#>   39 of 50 prediction points inside the AOA (78.0%)
#> 
#> Predictions outside the AOA are extrapolations; the cross-validated
#> performance estimate does not cover them.
table(res$aoa$AOA)
#> 
#> FALSE  TRUE 
#>    11    39 
```
