# The nugget of an estimated autocorrelation range

The nugget variance of the variogram model behind a
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
result: the model's semivariance as the separation shrinks to zero, i.e.
measurement error plus variation at scales shorter than the first lag of
the empirical variogram. It is extrapolated to zero from that lag, not
observed, so it is only as good as the first lag is short:
`attr(x, "first_lag")` says where the extrapolation starts. With
`lags = "adaptive"` (the default) that is the separation of the closest
pairs of points; with `lags = "equal"` it is the middle of a bin
`cutoff * max_dist / n_lags` wide, which on a dense layer can be many
times the spacing of close pairs, and structure shorter than that bin is
then counted as nugget.

## Usage

``` r
sac_nugget(x)
```

## Arguments

- x:

  A `sac_range` object, or anything else.

## Value

A single number: the nugget in the units of the response's variance;
`NA_real_` when `x` carries no fitted model (a bare `NA` from a run that
could not fit anything, a rejected result with too few lag bins to fit,
or an object that is not a `sac_range`).

## Details

Three more attributes say how far to trust it.
`attr(x, "nugget_at_bound")` is `TRUE` when the fit put the nugget at
its lower bound (0, or within 1e-4 of the sill): that is the bound and
not an estimate, and a criterion that uses it as a noise variance has no
penalty left. `attr(x, "short_lag_semivariance")` is the semivariance of
nearest-neighbour pairs, computed with no model: each location is paired
with the nearest other location (`short_lag_pairs` and `short_lag_dist`
give the number of such pairs and their median separation). A
semivariance does not fall as the lag shrinks, so this value bounds the
nugget from above, up to its sampling error; a fitted nugget well above
it is contradicted by the closest pairs.
`attr(x, "zero_lag_semivariance")` exists where locations were visited
more than once (`zero_lag_pairs` is not 0): it is the semivariance of
every pair of rows that share a location, that is, the variance between
visits to one place. That is the measurement part of the nugget,
observed and not extrapolated, without the variation at scales below the
spacing of the locations, so it bounds the nugget from below. Such pairs
are at no lag and are in no bin of the variogram: the model, and so the
nugget, describe positive separations. The variogram fit uses the bound,
though: where it rests on 30 degrees of freedom or more (rows less
locations) the nugget is not let under it, and `attr(x, "nugget_floor")`
holds the value (0 when there is no floor). With `detrend = "reml"` the
model is built from the REML parameters and takes no floor: its nugget
can come out below the variance between visits, and `nugget_floor` is 0.
A nugget equal to its floor is the variance between visits, a
measurement, and `nugget_at_bound` stays `FALSE` for it.

When `attr(x, "lag_fallback")` is `TRUE` the result is the equal-lag one
although adaptive lags were asked for (see "The lags" in
[`estimate_sac_range`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)),
and its nugget is extrapolated from the first equal bin.

The nugget is carried as the `nugget` attribute of every classed result,
identified or rejected, because it is the number a resolution criterion
for a tessellation needs (the short-lag variance that no cell can
average away). On a result refused for want of spatial structure it is
the variance of the modelled variable, which is the nugget of the model
that refusal supports.

## See also

[`estimate_sac_range`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
which produces the object.

Other cross-validation:
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md),
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
[`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md),
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md),
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
[`fold_separation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fold_separation.md),
[`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md),
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md),
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)

## Examples

``` r
if (requireNamespace("gstat", quietly = TRUE)) {
  library(sf)
  # A field with a real nugget: half a unit of white noise on a unit sill.
  set.seed(3)
  n <- 250
  xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
  D  <- as.matrix(dist(xy))
  xy$z <- as.numeric(t(chol(exp(-D / 150) + diag(0.5, n))) %*% rnorm(n))
  r <- estimate_sac_range(st_as_sf(xy, coords = c("x", "y"), crs = 32632), "z")
  print(sac_nugget(r))  # the fitted nugget variance, on the sill's scale
  attr(r, "nugget_at_bound")         # FALSE: an estimate, not the bound
  attr(r, "short_lag_semivariance")  # what the closest pairs say, no model
  attr(r, "zero_lag_pairs")          # 0: no location was visited twice
  sac_nugget(NA)        # nothing fitted: NA
}
#> [1] 0.4253438
#> [1] NA
```
