# Plot the interval calibration of a Bayesian cross-validation

[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
scores every fold's posterior predictive intervals at each of
`coverage_levels`: the share of held-out observations the interval of
that nominal level contained. This draws observed coverage against
nominal with the diagonal, one point per level for the fold-weighted
pooled value (from `predictive_coverage`) and one faint point per fold,
so you can see systematic over-confidence (points below the line) or
intervals wider than they need to be (above it). Three levels is a thin
curve; pass `coverage_levels = seq(0.1, 0.9, by = 0.1)` to
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
for a full one. Each level is drawn at the nominal value
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
records in `coverage_levels`, so whatever was computed is drawn where it
belongs (0.975 at 0.975, not rounded).

## Usage

``` r
plot_calibration(cv, ...)
```

## Arguments

- cv:

  The list returned by
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
  or a
  [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
  result that ran the Bayesian backend (its `$bayes_cv` is used).

- ...:

  Not used.

## Value

A `ggplot` object.

## See also

[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
which makes the result.

Other cross-validation:
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md),
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md),
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
[`plot_cv_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_cv_metrics.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Needs brms and a Stan toolchain, and takes minutes: run it, do not check it.
library(sf)
set.seed(1)
n <- 80
dat <- st_as_sf(
  data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
             a = rnorm(n)),
  coords = c("x", "y"), crs = 32632)
dat$z <- 2 * dat$a + rnorm(n)
# Nine coverage levels give a curve rather than three points; chains and
# iterations are kept small to keep this to a few minutes, so the intervals
# will be rough.
cv <- cv_bayes(dat, "z", "a", k = 3, coverage_levels = seq(0.1, 0.9, 0.1),
               fit_args = list(chains = 2, iter = 1000))
plot_calibration(cv)
} # }
```
