# Print an area-of-applicability result

Summarises where the model may be trusted: how many prediction locations
fall inside the area of applicability and how many outside, the
dissimilarity threshold that separated them, and the predictors the
index was computed over (naming any dropped for having no usable
variance). The proportion outside is the headline number – a map that
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
