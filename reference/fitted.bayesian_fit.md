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

Numeric vector of length `object$n` (all `NA` if the posterior draw
failed).

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
`predict(object, type = "predict")` recompute rather than reuse it. Call
[`clear_fitted_cache`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_fitted_cache.md)
if the engine has been mutated by hand after fitting.

## The cache is shared by copies, and validated

An environment has reference semantics, which is what makes the memo
survive R's copy-on-modify – but it also means `fit2 <- fit` gives the
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
on one copy empties the cache both share (harmless – the other simply
recomputes), and [`identical()`](https://rdrr.io/r/base/identical.html)
cannot distinguish two fits by their caches. The digest covers `data_sf`
only, not `$engine`: a hand-mutated `brmsfit` is what
[`clear_fitted_cache`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_fitted_cache.md)
is for.
