# Predict from a GWR spatial model

When `newdata` is NULL, returns the in-sample fitted values. Otherwise
uses
[`GWmodel::gwr.predict()`](https://rdrr.io/pkg/GWmodel/man/gwr.predict.html)
on the new locations. `newdata` is first transformed to the CRS used
during fitting (via
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)),
so predictions are computed in a single coordinate system regardless of
the CRS newdata arrives in.

## Usage

``` r
# S3 method for class 'gwr_fit'
predict(object, newdata = NULL, ...)
```

## Arguments

- object:

  A `gwr_fit` object.

- newdata:

  An sf object with the same predictors. The response variable need not
  be present (true out-of-sample prediction is supported). NULL = fitted
  values.

- ...:

  Ignored.

## Value

Numeric vector aligned to `nrow(newdata)`, with `NA` for rows dropped as
missing or non-finite. If
[`GWmodel::gwr.predict()`](https://rdrr.io/pkg/GWmodel/man/gwr.predict.html)
fails, every value is `NA` and a warning says why. CRS-less `newdata`
first receives the interpretation the training data got, so the same
rows land where they did at fit time.
