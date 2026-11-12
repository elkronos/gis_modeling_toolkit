# Extract the nugget of an estimated autocorrelation range

The nugget variance of the variogram model behind an
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
result: the model's semivariance as the separation shrinks to zero,
which is measurement error plus variation at scales shorter than the
first lag of the empirical variogram.

## Usage

``` r
sac_nugget(x)
```

## Arguments

- x:

  A `sac_range` object, or any other value.

## Value

A single number, the nugget in the units of the response's variance, or
`NA_real_` when `x` carries no fitted model: a bare `NA` from a run that
could not fit anything, a rejected result with too few lag bins to fit,
or an object that is not a `sac_range`.

## Details

The nugget is extrapolated to zero from the first lag, not observed, so
it is only as good as that lag is short; `attr(x, "first_lag")` says
where the extrapolation starts. With `lags = "adaptive"` (the default)
the first lag is the separation of the closest pairs of points. With
`lags = "equal"` it is the middle of a bin `cutoff * max_dist / n_lags`
wide, which on a dense layer can be many times the spacing of close
pairs, and structure shorter than that bin is then counted as nugget.
When `attr(x, "lag_fallback")` is `TRUE` the result is the equal-lag one
although adaptive lags were asked for (see the section "The lags" in
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)),
and its nugget is extrapolated from the first equal bin.

Every classed result, identified or rejected, carries the nugget as its
`nugget` attribute, since a resolution criterion for a tessellation
needs it: it is the short-lag variance that no cell can average away. On
a result refused for want of spatial structure it is the variance of the
modelled variable, the nugget of the pure-nugget model that the refusal
supports.

## How far to trust it

Four more attributes of the result bear on the nugget.

- `nugget_at_bound` is `TRUE` when the fit put the nugget at its lower
  bound (0, or at most 1e-4 of the sill). That value is the bound and
  not an estimate, and a criterion that uses it as a noise variance has
  no penalty left. A nugget equal to its floor (see `nugget_floor`) is
  the variance between visits, a measurement, and `nugget_at_bound`
  stays `FALSE` for it unless the floor itself is at most 1e-4 of the
  sill.

- `short_lag_semivariance` is the semivariance of nearest-neighbour
  pairs of locations, computed with no model; `short_lag_pairs` and
  `short_lag_dist` give the number of pairs and their median separation.
  Each location is paired with the nearest other location, and of
  several whose squared distances agree to one part in ten million, with
  the one of lower x coordinate, then lower y. On a square lattice that
  pairs almost every location with its neighbour along x, so the value
  measures the field along x, and with the axes of the CRS swapped it
  can differ. A semivariance does not fall as the lag shrinks, so this
  value bounds the nugget from above, up to its sampling error: a fitted
  nugget well above it is contradicted by the closest pairs.

- `zero_lag_semivariance` exists where locations were visited more than
  once (`zero_lag_pairs` is not 0). It is the semivariance of every pair
  of rows that share a location, the variance between visits to one
  place. That is the measurement part of the nugget, observed and not
  extrapolated, without the variation at scales below the spacing of the
  locations, so it bounds the nugget from below. Such pairs are at no
  lag and in no bin of the variogram, since the model, and so the
  nugget, describe positive separations.

- `nugget_floor` is the value the fit did not let the nugget under:
  `zero_lag_semivariance` where it rests on 30 degrees of freedom or
  more (rows less locations), and 0 otherwise. With `detrend = "reml"`
  the model is built from the REML parameters and takes no floor: its
  nugget can come out below the variance between visits, and
  `nugget_floor` is 0.

## See also

[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
which produces the object.

Other spatial autocorrelation:
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
[`plot.sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.sac_range.md)

## Examples

``` r
library(sf)
# A field with a nugget: half a unit of white noise on a unit sill.
set.seed(3)
n <- 250
xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
D  <- as.matrix(dist(xy))
xy$z <- as.numeric(t(chol(exp(-D / 150) + diag(0.5, n))) %*% rnorm(n))
r <- estimate_sac_range(st_as_sf(xy, coords = c("x", "y"), crs = 32632), "z")
sac_nugget(r)                      # the fitted nugget, on the sill's scale
#> [1] 0.4253438
attr(r, "nugget_at_bound")         # FALSE: an estimate, not the bound
#> [1] FALSE
attr(r, "short_lag_semivariance")  # what the closest pairs say, no model
#> [1] 0.7541593
attr(r, "zero_lag_pairs")          # 0: no location was visited twice
#> [1] 0
sac_nugget(NA)                     # nothing fitted: NA
#> [1] NA
```
