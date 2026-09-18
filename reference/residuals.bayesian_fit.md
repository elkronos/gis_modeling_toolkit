# In-sample residuals from a Bayesian spatial GP fit

Observed response minus
[`fitted.bayesian_fit`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.bayesian_fit.md)
(the cached posterior expectation), so these are **in-sample**
residuals.

## Usage

``` r
# S3 method for class 'bayesian_fit'
residuals(object, ...)
```

## Arguments

- object:

  A `bayesian_fit`.

- ...:

  Ignored.

## Value

Numeric vector of length `object$n`.
