# Print a spatial autocorrelation range

Prints the effective range as a plain number, with the directional fit
summarised beneath it when one is available. A direction whose fit was
unusable is labelled with why (`directional_status`) and the range its
fit reported (`directional_fitted`) when the object carries them, and
`unidentified` otherwise.

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
