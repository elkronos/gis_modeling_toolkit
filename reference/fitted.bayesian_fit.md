# In-sample fitted values from a Bayesian spatial GP fit

Posterior expectation at the training locations: the column means of
[`brms::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html).
These are **in-sample** values.

## Usage

``` r
# S3 method for class 'bayesian_fit'
fitted(object, ...)
```

## Arguments

- object:

  A `bayesian_fit`.

- ...:

  Ignored.

## Value

Numeric vector of length `object$n`. A posterior that cannot be drawn is
an error, as is a family with a probability per response category
(ordinal, categorical), which has no single fitted value per row.

## The result is cached

`posterior_epred()` is O(draws x n), and
[`summary()`](https://rdrr.io/r/base/summary.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
and
[`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md)
each call [`fitted()`](https://rdrr.io/r/stats/fitted.values.html)
independently, so the value is memoised in an environment carried in
`object$info$.cache` (reference semantics, so it survives R's
copy-on-modify). The cache holds epred column means only, which is why
`predict(object, summary = "median")` and
`predict(object, type = "predict")` recompute it from scratch. Call
[`clear_fitted_cache`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_fitted_cache.md)
if the engine has been mutated by hand after fitting.

## The cache is shared by copies, and validated

An environment has reference semantics, which is what makes the memo
survive R's copy-on-modify. But it also means `fit2 <- fit` gives the
two objects *the same* cache. Assigning a different `data_sf` to the
copy would then have returned the original's cached values, at the
original's length, which
[`residuals()`](https://rdrr.io/r/stats/residuals.html) silently
recycled against the copy's shorter response. The entry therefore
carries the `n` and a digest of the training data it was computed from,
and is recomputed whenever either fails to match, so a copy with
different data recomputes instead of reading the original's answer.

Two consequences of the shared environment remain and cannot be removed
from here:
[`clear_fitted_cache`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_fitted_cache.md)
on one copy empties the cache both share (harmless, since the other
simply recomputes), and
[`identical()`](https://rdrr.io/r/base/identical.html) cannot
distinguish two fits by their caches. The digest covers `data_sf` only.
The entry is also tied to the engine that computed it – a refit or
[`update()`](https://rdrr.io/r/stats/update.html) of the `brmsfit` is a
different sampling run and recomputes – but a `brmsfit` edited by hand
in place is what
[`clear_fitted_cache`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_fitted_cache.md)
is for. The entry holds only the values and a small identifier of the
sampling run, so a fit saved with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) after
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) is no larger
for it.

## See also

Other methods on a fitted model:
[`coef.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.bayesian_fit.md),
[`coef.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.gwr_fit.md),
[`coef.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.rf_fit.md),
[`fitted.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.gwr_fit.md),
[`fitted.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.rf_fit.md),
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md),
[`predict.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.bayesian_fit.md),
[`predict.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.gwr_fit.md),
[`predict.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.rf_fit.md),
[`print.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.rf_fit.md),
[`print.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.spatial_fit.md),
[`residuals.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.bayesian_fit.md),
[`residuals.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.gwr_fit.md),
[`residuals.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.rf_fit.md),
[`summary.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.spatial_fit.md)
