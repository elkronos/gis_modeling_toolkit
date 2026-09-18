# Clear cached fitted values for a Bayesian spatial model

Removes the lazily-cached
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) result so that
the next call recomputes from the posterior. This is only necessary if
the underlying `brmsfit` engine has been manually mutated after fitting
– a change to `data_sf` invalidates the entry on its own, because the
cached value carries a digest of the data it was computed from (see
[`fitted.bayesian_fit`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.bayesian_fit.md)).
Normal usage never requires it.

## Usage

``` r
clear_fitted_cache(object)
```

## Arguments

- object:

  A `bayesian_fit` object.

## Value

`object`, invisibly (called for side effect).

## Details

The cache environment is shared by every copy of a fit, so clearing it
through one copy clears it for all of them. That is harmless: the others
recompute.

## Examples

``` r
# Only a bayesian_fit carries the cache; on any other fit this is a no-op.
if (requireNamespace("ranger", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  pts <- st_as_sf(
    data.frame(x = runif(60, 0, 1000), y = runif(60, 0, 1000), a = rnorm(60)),
    coords = c("x", "y"), crs = 32632
  )
  pts$z <- 2 * pts$a + rnorm(60, 0, 0.3)
  fit <- fit_rf_model(pts, "z", "a", num_trees = 50, seed = 1)
  clear_fitted_cache(fit)
}
```
