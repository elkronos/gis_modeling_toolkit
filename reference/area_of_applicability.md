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
  weights = NULL,
  folds = NULL,
  threshold = NULL,
  normalizer_max_n = 5000L,
  seed = 123L,
  chunk_size = NULL,
  use_fnn = requireNamespace("FNN", quietly = TRUE)
)
```

## Arguments

- newdata:

  Prediction locations: an `sf` object (typically from
  [`predict_surface`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md))
  or a data.frame, carrying the predictor columns. For a model fitted
  with `include_coords = TRUE` an `sf` object is required, non-`POINT`
  geometry is reduced to representative points, and a CRS mismatch with
  the training data is reconciled; see *Models fitted with the
  coordinates as predictors*.

- model:

  A fitted `spatial_fit`, supplying the training data and predictor
  names. Optional if `train_sf` and `predictor_vars` are given directly,
  which lets this be used with any model.

- train_sf:

  Training data, if not taken from `model`.

- predictor_vars:

  Predictor names, if not taken from `model`.

- weights:

  Optional named numeric vector of predictor importances. Any positive
  scale works. When `model` was fitted with `include_coords = TRUE`, the
  coordinates `"..x"` and `"..y"` are measured as predictors too (see
  below), and you are not expected to supply importances for them: any
  you leave out default to the mean of the weights you did supply, so
  location counts about as much as a typical predictor. Naming them
  explicitly overrides that. An unnamed vector may have one value per
  predictor either with or without the two coordinate columns.

- folds:

  Cross-validation folds: a
  [`make_folds`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  result, a list of `train`/`test` splits, or a vector of fold labels
  with one entry per training row. Default `NULL` (plain nearest
  neighbour).

- threshold:

  Optional numeric override for the DI threshold.

- normalizer_max_n:

  Subsample the training data to this many points when computing the
  mean pairwise distance, which is quadratic. Default 5000.

- seed:

  Seed for that subsample. Default 123.

- chunk_size:

  Query rows per distance block on the dense path. Default `NULL`
  (chosen from the training size).

- use_fnn:

  Use FNN for nearest-neighbour search when available. Exposed so the
  dense fallback can be tested.

## Value

An object of class `aoa`: a list with

- `aoa`: `newdata` with a numeric `DI` column and a logical `AOA` column
  added. This is the object the computation ran on, which for a
  coordinate-using model is `newdata` after pointizing, CRS
  reconciliation and the addition of the `"..x"` and `"..y"` columns. A
  row whose predictors are not all finite gets `NA` in both columns.

- `threshold`: the DI cut-off used.

- `train_DI`: the training points' own DI values.

- `normalizer`: the mean pairwise training distance.

- `weights`: the weight vector actually applied, named by
  `predictor_vars`.

- `predictor_vars`: the predictors used, including `"..x"`/`"..y"` when
  the model uses coordinates and excluding `dropped_vars`.

- `dropped_vars`: predictors dropped for negligible variance.

- `scaling`: a list with `center` and `scale`, each named by
  `predictor_vars`: the training means and standard deviations the index
  is computed in, so a location's DI can be traced to the predictor that
  put it outside.

- `n_outliers`: the number of training DI values above the
  `Q3 + 1.5 * IQR` fence, which the default threshold rule sets aside
  (the "outlier-removed" in its name); computed whether or not
  `threshold` was supplied.

- `n_train`, `n_new`, `n_inside`, `n_outside`, `n_na`: row counts;
  `n_train` and `n_new` count the rows that survived the finite-value
  filter.

- `params` records the call: `folds_supplied`, `n_folds`,
  `folds_method`, `threshold_supplied`, `normalizer_max_n`,
  `normalizer_n_used`, `normalizer_subsampled`, `weights_supplied` and
  `seed`. `folds_method` is the `method` of a
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
and standard deviations, then optionally weighted by variable
importance. For a prediction point \\p\\, the DI is the distance to its
nearest training point in that space, divided by the mean pairwise
distance among training points. The same quantity is computed for the
training data itself, using each point's nearest neighbour *among the
training rows of the fold that holds it out*. That means everything
outside its own fold for random and block folds, and the smaller
training set that buffered and NNDM folds actually leave (see the next
section). The threshold is then the largest training DI that is not an
upper outlier. Prediction points at or below that threshold are inside
the AOA.

The DI is invariant to the overall scale of `weights`: the numerator and
the normaliser carry the same factor. Importance values can be passed
as-is.

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

## Limitations

Predictors must be numeric; categorical variables are refused and never
silently dummy-coded. Predictors whose variance is negligible *relative
to their own magnitude* (the test is
`sd < sqrt(.Machine$double.eps) * max(abs(x))`, so the same variable in
metres and in gigametres is treated identically) are dropped, and a
prediction point taking a different value there is a form of
extrapolation this index cannot express. Without `weights` every
predictor counts equally, which overstates dissimilarity along
directions the model barely uses.

## Models fitted with the coordinates as predictors

When `model` was fitted with `include_coords = TRUE` the model splits on
location, so the dissimilarity index has to measure location too: the
coordinates are added to both sides as the predictors `"..x"` and
`"..y"` and are then centred, scaled and weighted like any other column.
Without this a prediction point far outside the training extent but with
ordinary covariate values reads as *inside* the area of applicability.
That is exactly the extrapolation this index exists to catch.

This path needs geometry on both sides, so `train_sf` and `newdata` must
both be `sf` objects; a data.frame is refused rather than quietly
measured without location. Non-`POINT` `newdata` (grid polygons, say) is
reduced to representative points first, as
[`coerce_to_points`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md)
would. If exactly one side carries a CRS the other is brought into it:
reprojected when its coordinates look like longitude/latitude, stamped
otherwise, with a warning either way. This is done because degrees fed
into a metre-space index silently understate the distances.

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
[`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md),
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md),
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)

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
#>   weighted    : no (all predictors count equally)
#>   training    : 120 points
#>   reference   : nearest other training point (no folds supplied)
#>   normaliser  : 1.7777 (mean pairwise distance)
#>   threshold   : 0.2899 (outlier-removed max of training DI)
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
