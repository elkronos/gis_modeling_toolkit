# The nugget of an estimated autocorrelation range

The nugget variance of the variogram model behind a
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
result: the semivariance at zero separation, i.e. measurement error plus
variation at scales shorter than the closest pair. It is carried as the
`nugget` attribute of every classed result, identified or rejected,
because it is the number a resolution criterion for a tessellation needs
(the short-lag variance that no cell can average away).

## Usage

``` r
sac_nugget(x)
```

## Arguments

- x:

  A `sac_range` object, or anything else.

## Value

A single number: the nugget in the units of the response's variance;
`NA_real_` when `x` carries no fitted model (a bare `NA` from a run that
could not fit anything, a rejected result whose fits were all singular,
or an object that is not a `sac_range`).

## See also

[`estimate_sac_range`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
which produces the object.

## Examples

``` r
if (requireNamespace("gstat", quietly = TRUE)) {
  library(sf)
  set.seed(9)
  n <- 150
  xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
  D  <- as.matrix(dist(xy))
  xy$z <- as.numeric(t(chol(exp(-D / 100) + diag(0.1, n))) %*% rnorm(n))
  r <- estimate_sac_range(st_as_sf(xy, coords = c("x", "y"), crs = 32632), "z")
  sac_nugget(r)
  sac_nugget(NA)        # nothing fitted: NA
}
#> [1] NA
```
