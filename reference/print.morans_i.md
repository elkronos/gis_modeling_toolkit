# Print a residual Moran's I result

Shows the statistic, the null its moments come from and the evidence
behind them. The weight matrix the result carries gets a one-line
description in place of the matrix: it is \\n \times n\\, and
autoprinting it buried the statistic under a thousand lines of matrix.

## Usage

``` r
# S3 method for class 'morans_i'
print(x, ...)
```

## Arguments

- x:

  An object of class `morans_i`, from
  [`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md).

- ...:

  Ignored.

## Value

`x`, invisibly.
