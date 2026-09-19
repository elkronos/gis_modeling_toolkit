# Print a random forest fit

Shows the forest's shape (formula, n, number of trees, `mtry`, node
size) along with the out-of-bag error and, prominently, whether the
coordinates were used as predictors. That last line is the one to check:
a forest fitted with `include_coords = TRUE` can memorise location and
score well out-of-bag while failing everywhere it has not been.

## Usage

``` r
# S3 method for class 'rf_fit'
print(x, ...)
```

## Arguments

- x:

  An `rf_fit`.

- ...:

  Ignored.

## Value

`x`, invisibly.
