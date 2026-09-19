# Print a GWR model selection result

Shows the forward-selection trail: the response, the candidate
predictors, and the top-ranked models with their criterion values, so
you can see both which model won and by how much. A shallow gap between
the first few rows means the ranking is not well identified and the
choice of predictors should not be treated as settled. Check the gap
before reporting one model as the selected one.

## Usage

``` r
# S3 method for class 'gwr_model_selection'
print(x, n = 10L, ...)
```

## Arguments

- x:

  A `gwr_model_selection` object.

- n:

  Number of top-ranked models to show. Default 10.

- ...:

  Ignored.

## Value

`x`, invisibly.

## See also

Other print methods:
[`print.aoa()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.aoa.md),
[`print.morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.morans_i.md),
[`print.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.resolution_profile.md),
[`print.sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.sac_range.md)
