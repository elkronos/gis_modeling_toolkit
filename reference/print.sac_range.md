# Print a spatial autocorrelation range

Prints the effective range as a plain number, with the directional fit
summarised beneath it when one is available. A direction whose fit was
unusable is labelled with why (`directional_status`) and the range its
fit reported (`directional_fitted`) when the object carries them, and
`unidentified` otherwise. A last line names the unit and the CRS the
range is a length in (`attr(x, "crs")`, which for lon/lat input is the
projected CRS the estimate chose; for a layer with no CRS, it says the
range is in that layer's own coordinate units), and whether the
variogram is of the response or of its residuals on `predictor_vars`
(`detrended`, `detrend_method`).

## Usage

``` r
# S3 method for class 'sac_range'
print(x, ...)
```

## Arguments

- x:

  An object of class `sac_range`.

- ...:

  Ignored.

## Value

`x`, invisibly.

## See also

Other print methods:
[`print.aoa()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.aoa.md),
[`print.gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.gwr_model_selection.md),
[`print.morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.morans_i.md),
[`print.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.resolution_profile.md)
