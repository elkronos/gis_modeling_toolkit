# Clear the cached fitted values of a Bayesian GP fit

Removes the cached
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) values of a
`bayesian_fit`, so the next call recomputes them from the posterior. You
need it only after editing the `brmsfit` engine in place. A change to
`data_sf`, or a refit, invalidates the entry on its own (see "The
Bayesian GP cache" in
[spatial_fit-fitted](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-fitted.md)).
On any other fit it does nothing.

## Usage

``` r
clear_fitted_cache(object)
```

## Arguments

- object:

  A `bayesian_fit` object.

## Value

`object`, invisibly.

## Details

Every copy of a fit shares the cache, so clearing it through one copy
clears it for all of them; the others recompute.

## See also

Other options and caches:
[`clear_grid_cache()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_grid_cache.md),
[`spatialkit-conditions`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-conditions.md),
[`spatialkit-options`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-options.md),
[`spatialkit_quiet()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Needs brms and a Stan toolchain, and takes minutes.
library(sf)
set.seed(1)
n <- 60
dat <- st_as_sf(
  data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
             elev = rnorm(n)),
  coords = c("x", "y"), crs = 32632
)
dat$price <- 10 + 2 * dat$elev + rnorm(n)
fit <- fit_bayesian_spatial_model(dat, "price", "elev", chains = 2,
                                  iter = 1000, compute_loo = FALSE)
f1 <- fitted(fit)          # computed from the posterior and cached
clear_fitted_cache(fit)
f2 <- fitted(fit)          # computed again
all.equal(f1, f2)
} # }
```
