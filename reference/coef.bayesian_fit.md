# Extract Bayesian model fixed-effect summaries

Returns the posterior summary of the global (non-spatial) regression
terms: estimate, error and credible interval per predictor, as
[`brms::fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html)
reports them. Reach for it to read the average effect of a predictor
with its uncertainty attached (the Bayesian counterpart to a coefficient
table), remembering that the Gaussian-process term has already absorbed
the spatially structured part of the signal, so these are effects net of
location.

## Usage

``` r
# S3 method for class 'bayesian_fit'
coef(object, ...)
```

## Arguments

- object:

  A `bayesian_fit` object.

- ...:

  Ignored.

## Value

A matrix of fixed-effect posterior summaries, as returned by
[`brms::fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html).
Never `NULL`: a missing 'brms' or a failing `fixef()` call errors,
following the [`coef()`](https://rdrr.io/r/stats/coef.html) contract
described in
[`new_spatial_fit`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md).
