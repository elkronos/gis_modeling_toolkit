# Out-of-bag residuals from a random forest fit

Observed response minus
[`fitted.rf_fit`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.rf_fit.md),
which for a forest is the **out-of-bag** prediction – each observation
predicted only by the trees that did not see it. These are therefore
already held-out residuals, unlike
[`residuals.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.gwr_fit.md)
and
[`residuals.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.bayesian_fit.md),
which are in-sample. Feed them to
[`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)
to test whether spatial structure the forest failed to capture is still
sitting in the residuals. Out-of-bag is not a substitute for spatial CV:
use
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
for an honest map-accuracy figure.

## Usage

``` r
# S3 method for class 'rf_fit'
residuals(object, ...)
```

## Arguments

- object:

  An `rf_fit`.

- ...:

  Ignored.

## Value

Numeric vector of length `object$n`.
