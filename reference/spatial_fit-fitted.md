# Extract the fitted values of a fitted model

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) on a random
forest, a GWR or a Bayesian GP fit returns one value per training row,
in the order of `object$data_sf`. What the value is depends on the
backend, and so does what
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
without `newdata` report:

- a random forest returns out-of-bag predictions: each row predicted by
  the trees that did not sample it (see "Fitted values are out-of-bag"
  in
  [`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)).
  A row that every tree sampled has no out-of-bag prediction and is
  `NaN`;

- a GWR returns the values
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  stored on the fit: at each training location, the value of the local
  regression centred there. These are in-sample values (each observation
  sits inside its own window), so
  [`summary()`](https://rdrr.io/r/base/summary.html), and
  [`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
  without `newdata`, report an optimistic fit; use
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md)
  for a spatially blocked estimate. An engine made elsewhere (a GWmodel
  result) is read from its `SDF` table instead;

- a Bayesian GP returns the posterior expectation at the training
  locations, the column means of
  [`brms::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html),
  in the units of the response (see "The response's centre and scale" in
  [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)).
  These are in-sample values.

## Usage

``` r
# S3 method for class 'gwr_fit'
fitted(object, ...)

# S3 method for class 'bayesian_fit'
fitted(object, ...)

# S3 method for class 'rf_fit'
fitted(object, ...)
```

## Arguments

- object:

  A fitted `rf_fit`, `gwr_fit` or `bayesian_fit`.

- ...:

  Not used.

## Value

A numeric vector of length `object$n`. For a GWR, `NA` where the values
cannot be read from the engine. For a Bayesian GP, a posterior that
cannot be drawn is an error, as is a family with a probability per
response category (ordinal, categorical), which has no single fitted
value per row.

## The Bayesian GP cache

[`brms::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
takes time in proportion to the number of draws times the number of
rows, and [`summary()`](https://rdrr.io/r/base/summary.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
and
[`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md)
each call [`fitted()`](https://rdrr.io/r/stats/fitted.values.html), so
the values are kept in an environment carried in `object$info$.cache`.
The cache holds the expected values only, so
`predict(object, summary = "median")` and
`predict(object, type = "predict")` recompute from the posterior.

The entry records the `n` and a digest of the training data it was
computed from, and is recomputed when either differs, so a copy given
other data does not read the original's values. It is also tied to the
sampling run: a refit or
[`update()`](https://rdrr.io/r/stats/update.html) of the `brmsfit`
recomputes. A `brmsfit` edited in place is what
[`clear_fitted_cache()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_fitted_cache.md)
is for. A copy made with `fit2 <- fit` shares the environment, so
[`clear_fitted_cache()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_fitted_cache.md)
on one empties the cache of both (the other recomputes), and
[`identical()`](https://rdrr.io/r/base/identical.html) cannot tell two
fits apart by their caches. The entry holds the values and a small
identifier of the sampling run, so a fit saved with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) after
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) is no larger
for it.

## See also

[`residuals()`](https://rdrr.io/r/stats/residuals.html) on the same fits
([spatial_fit-residuals](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-residuals.md)),
and
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md).

Other fit methods:
[`plot.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.spatial_fit.md),
[`print.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.spatial_fit.md),
[`spatial_fit-coef`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-coef.md),
[`spatial_fit-predict`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-predict.md),
[`spatial_fit-residuals`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-residuals.md),
[`summary.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.spatial_fit.md)

## Examples

``` r
library(sf)
set.seed(1)
pts <- st_as_sf(
  data.frame(x = 5e5 + runif(60, 0, 1000), y = 5e6 + runif(60, 0, 1000),
             a = rnorm(60)),
  coords = c("x", "y"), crs = 32632
)
pts$z <- 10 + 2 * pts$a + rnorm(60, 0, 0.3)
rf <- fit_rf_model(pts, "z", "a", num_trees = 50, seed = 1)
gw <- fit_gwr_model(pts, "z", "a", bandwidth = 30)
# Out-of-bag for the forest, in-sample for the GWR.
head(cbind(z = pts$z, rf = fitted(rf), gwr = fitted(gw)))
#>              z        rf       gwr
#> [1,] 14.651448 14.029639 14.619169
#> [2,] 10.324432  9.960532  9.882031
#> [3,] 11.315105 11.204729 11.347377
#> [4,] 10.002137  9.993318 10.070932
#> [5,]  8.483396  8.534352  8.563605
#> [6,] 10.591384 10.239645 10.411187
```
