# Print a GWR model selection result

Shows the forward-selection trail: the response, the candidate
predictors, and the top-ranked models with their criterion values and,
in the `delta` column, how far each is from the best. Below the table it
lists the models within 2 of the best, which a difference in AICc does
not separate (Burnham and Anderson 2002), and names the smallest of
them. When several models are that close the ranking is not well
identified and the choice of predictors should not be treated as
settled; when the model without predictors (`<none>`) is among them, no
candidate has been shown to help.

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

## References

Burnham, K. P. and Anderson, D. R. (2002). *Model Selection and
Multimodel Inference: A Practical Information-Theoretic Approach*, 2nd
edition. Springer.

## See also

Other print methods:
[`print.aoa()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.aoa.md),
[`print.morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.morans_i.md),
[`print.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.resolution_profile.md),
[`print.sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.sac_range.md)
