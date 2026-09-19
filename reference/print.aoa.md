# Print an area-of-applicability result

Summarises where the model may be trusted: how many prediction locations
fall inside the area of applicability and how many outside, the
dissimilarity threshold that separated them, and the predictors the
index was computed over (naming any dropped for having no usable
variance). The proportion outside is the headline number. A map that
extrapolates over much of its extent is reporting predictions its
training data cannot support, whatever the cross-validation score said.

## Usage

``` r
# S3 method for class 'aoa'
print(x, ...)
```

## Arguments

- x:

  An `aoa` object.

- ...:

  Ignored.

## Value

`x`, invisibly.

## See also

Other print methods:
[`print.gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.gwr_model_selection.md),
[`print.morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.morans_i.md),
[`print.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.resolution_profile.md),
[`print.sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.sac_range.md)
