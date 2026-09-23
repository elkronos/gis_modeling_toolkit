# Estimate the spatial autocorrelation range from data

Fits exponential (or spherical) variogram models and returns the
*effective range*: for the exponential model, three times the fitted
range parameter, which is where the semivariance reaches ~95 % of the
sill; for the spherical model (fitted only when the exponential fit is
singular) the fitted range itself, which is where the spherical
semivariance reaches its sill exactly. Both are the distance beyond
which two observations are (near) uncorrelated, which is what a block or
a buffer has to exceed.

## Usage

``` r
estimate_sac_range(
  points_sf,
  response_var,
  predictor_vars = NULL,
  n_max = 5000L,
  cutoff = 0.5,
  range_frac = 1,
  seed = 123L,
  detrend = c("ols", "reml"),
  reml_max_n = 400L,
  keep_directional_fits = FALSE
)
```

## Arguments

- points_sf:

  An sf object with point geometries (will be projected automatically if
  in geographic CRS). Non-POINT geometry is reduced to representative
  points; any Z or M dimension is dropped, because
  [`gstat::variogram()`](https://r-spatial.github.io/gstat/reference/variogram.html)
  uses every coordinate dimension and an XYZ layer would otherwise
  return a range in 3-D while every consumer of it works in 2-D map
  distance; and rows with empty or non-finite coordinates are dropped
  with a logged count.

- response_var:

  Character(1) name of the response column. For a count or other
  response whose variance tracks its mean, see the section on
  non-Gaussian responses: the range is still estimated, but it is a less
  reliable number than for a Gaussian response.

- predictor_vars:

  Optional character vector. When supplied, the trend on these
  predictors is removed first and the variogram describes the residual
  autocorrelation, the part a spatial model has to handle once the
  covariates have done their work. How the trend is removed is set by
  `detrend`, and it matters: see "Detrending and the residual-variogram
  bias".

- n_max:

  Maximum number of points to subsample before fitting. Variogram
  estimation is O(n²) so this keeps runtime bounded.

- cutoff:

  Fraction of the maximum inter-point distance to use as the variogram
  lag cutoff. That distance is the farthest pair, found on the convex
  hull, and not the bounding-box diagonal, which depends on how the axes
  are oriented. Default 0.5.

- range_frac:

  Positive numeric. A fitted range exceeding
  `range_frac * cutoff * max_dist` (that is, beyond the longest lag the
  empirical variogram was actually fitted over) is treated as
  unidentified and `NA_real_` is returned.
  [`gstat::fit.variogram()`](https://r-spatial.github.io/gstat/reference/fit.variogram.html)
  yields a finite number even when the variogram never reaches a sill,
  and such a value extrapolates past the observed lags instead of
  measuring a long autocorrelation range. Passing it to
  `make_folds(auto_range = TRUE)` would collapse the block grid to a
  single block. Default 1.0; raise it to accept ranges extrapolated
  beyond the fitted lags.

- seed:

  RNG seed for the `n_max` subsample, restored afterwards so the
  caller's random stream is untouched. Default `123L`: the subsample is
  an internal approximation and no part of the answer, and leaving it
  unseeded made the returned range differ between runs on identical
  input (19531, 19589, 19605 on three calls) and silently advanced the
  caller's RNG. Pass `NULL` for the old unseeded behaviour, or a
  different number to check how sensitive the estimate is to the
  subsample. Ignored when `nrow(points_sf) <= n_max`, where nothing is
  sampled. The `reml_max_n` subsample uses the same seed.

- detrend:

  How the trend on `predictor_vars` is removed; ignored when there are
  none. `"ols"` (default) fits it by ordinary least squares and fits the
  variogram to the residuals. This is the long-standing behaviour, which
  underestimates the range (see the section below). `"reml"` fits the
  trend and an exponential-plus-nugget covariance together by residual
  maximum likelihood with
  [`nlme::gls()`](https://rdrr.io/pkg/nlme/man/gls.html), returns the
  REML range, and attaches the empirical variogram of the REML residuals
  for inspection. It needs nlme, costs \\O(n^3)\\ (about 2 s at 300
  points, 12 s at 500, 45 s at 800), and so runs on at most `reml_max_n`
  points; when it does not converge the `"ols"` path runs instead with
  an R warning saying so.

- reml_max_n:

  Positive integer, at least 30. With `detrend = "reml"`, the trend and
  covariance are fitted on a seeded random subsample of this many points
  when the layer has more (exact duplicate locations are dropped first);
  the fitted trend is then removed from every point. Default 400. Raise
  it for a better-determined fit at the cost above.

- keep_directional_fits:

  Logical. Attach each direction's empirical variogram and fitted model
  as `directional_fits`? Defaults to `FALSE`: the four variograms are
  most of the object's size (42.1 KB of 59.3 KB at \\n = 400\\, and the
  difference between a 52.4 KB and a 94.6 KB
  `make_folds(auto_range = TRUE)` result), while the numbers read from
  them (`directional`, `directional_fitted`, `directional_status`,
  `anisotropy`) are attached either way, and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
  effective variogram from its own attribute. Set `TRUE` to inspect the
  directional curves.

## Value

A single number, of class `sac_range` in the first two of the three
shapes below and a bare `NA` in the third; all three behave as an
ordinary number. The shapes carry different attributes:

- Success:

  A positive effective range in projected coordinate units, with the fit
  attached as attributes `directional` (the 0°, 45°, 90° and 135°
  ranges, named by azimuth; `NA` where that direction's fit was
  unusable), `anisotropy` (largest over smallest), `anisotropy_used`
  (logical: `TRUE` only when the all-pairs fit was unusable and the
  directional maximum stands in for it), `directional_status` (per
  azimuth, why a direction is `NA` in `directional`: `"ok"`,
  `"over_cutoff"` (its range ran past the largest lag fitted),
  `"not_converged"` or `"no_fit"`), `directional_fitted` (the range each
  direction's fit reported whether or not it was usable, so a refused
  directional range stays recoverable) and `directional_fits` (a list by
  azimuth of each direction's empirical `variogram` and fitted `model`,
  `NULL` where there is none, and `NULL` altogether unless
  `keep_directional_fits = TRUE`), `detrended` (logical: whether the
  variogram is of the residuals on `predictor_vars` or of the raw
  response. A missing predictor is an error, and a failed detrending fit
  warns and falls back to the raw response with this set to `FALSE`),
  `detrend_method` (`"ols"` or `"reml"` when detrended, `NA` otherwise),
  `reml` (with `detrend = "reml"`: a list with `n_used`, `subsampled`,
  `nugget_prop` and `sigma2` from the REML fit; `NULL` otherwise), `crs`
  (the projected CRS the variogram was fitted in: the unit of the
  range), `max_dist`, `cutoff_dist`, `variogram` (the empirical
  variogram), `variogram_model` (the fitted `gstat` model, or with
  `detrend = "reml"` a `gstat` model built from the REML parameters) and
  `nugget` (that model's nugget variance; see
  [`sac_nugget`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md)),
  so the fit can be inspected and need not be taken on trust.

- Rejected range:

  `NA_real_` when a range was fitted but is not identified: it exceeds
  `range_frac * cutoff * max_dist` (see `range_frac`); or the model did
  not converge; or the empirical variogram *decreases* with distance
  over its shorter lags (a net fall of more than 15 percent of the mean
  semivariance there, weighted by pairs), which is the shape of a
  periodic, hole-effect structure or of a variance that differs between
  a dense cluster and the rest of the layer (an unremoved trend instead
  makes the variogram rise without a sill, and the first test catches
  that); or the fitted range is non-positive. It is classed `sac_range`
  as well, so it prints as a bare `NA` without dumping its attributes,
  and it carries `max_dist`, `cutoff_dist`, `variogram`,
  `variogram_model` and `nugget` (the evidence for the rejection), plus
  `rejected_range` (the value that was refused), `rejected_reason` (one
  of `"fitted range exceeds the largest lag fitted"`,
  `"variogram model did not converge"`,
  `"empirical variogram decreases with distance"`,
  `"fitted range is non-positive or non-finite"`,
  `"no variogram model could be fitted (singular fits)"`), `crs` (so the
  units the rejected number was in stay recoverable, which is what
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) labels its
  axis from) and `detrend_method`. It carries `directional`,
  `anisotropy`, `anisotropy_used`, `directional_status`,
  `directional_fitted` and, with `keep_directional_fits = TRUE`,
  `directional_fits` as well: the directional sweep runs whatever
  becomes of the all-pairs fit, and its per-azimuth outcome is what says
  whether any direction reached a sill the pooled variogram did not, or
  whether every direction ran past the fitted lags alike. The same
  shape, with `rejected_range = NA`, `variogram_model = NULL` and
  `nugget = NA`, is returned when no variogram model could be fitted at
  all (both the exponential and the spherical fit singular, which is
  what a flat, nugget-only variogram produces); `rejected_reason` says
  so and the empirical variogram is still attached.

- No fit:

  A bare, attribute-less `NA_real_` when estimation could not be
  attempted at all: gstat missing, fewer than 30 finite values, a
  variable with no variance, or a degenerate extent. Without gstat
  nothing is fitted, so none of the attributes above exist either.

Attributes and the class do not affect
[`is.na()`](https://rdrr.io/r/base/NA.html) or
[`is.finite()`](https://rdrr.io/r/base/is.finite.html), so every
downstream guard treats all three the same way it always did.

## Details

The estimate is the **omnidirectional** (all-pairs) fit. Directional
variograms are fitted as well, at 0° (N–S), 45°, 90° (E–W) and 135°
azimuths with a ±22.5° tolerance. Those four windows tile all 180
distinct azimuths exactly once, and their ranges are returned in the
`directional` attribute, with their largest-over-smallest ratio in
`anisotropy`. They are a diagnostic, not the answer, for two reasons.
Each direction sees about a quarter of the point pairs, and the maximum
of four quarter-sample fits is biased upward: on simulated *isotropic*
fields it came in about 40% above the truth, and no hurdle placed in
front of it (all four directions fitted, ratio above 1.5, maximum above
1.5× the all-pairs fit) kept it out. One isotropic field rotated in 10°
steps "established" anisotropy in 14 of 18 orientations. And the windows
are fixed to the coordinate axes, so any answer built from them changes
when the layer is rotated, which a property of the field must not do.
The all-pairs fit is the best-powered estimate available and is
invariant to rotation.

Where a field is *known* to be anisotropic, blocks must be at least as
large as the longest autocorrelation range to avoid leakage, and the
conservative choice is to size them from
`max(attr(range, "directional"))` explicitly. A ratio above 1.5 is
logged so the case is not missed, with that advice. Only when the
omnidirectional fit is itself unusable is the directional maximum
returned in its place, and `anisotropy_used` is `TRUE` in that case
alone.

A direction whose fit fails, does not converge, or reports a range
beyond the longest fitted lag is excluded and recorded as `NA` in the
`directional` attribute.

Every variogram model is fitted **with a nugget**. A nugget-free model
forces the curve through the origin, and on any real measurement (which
has one) gstat's default N/h² weights buy that constraint by collapsing
the range: with a 50% nugget the fitted range came back at about 0.45 of
the truth, so `make_folds(auto_range = TRUE)` built blocks less than
half the correlation length it reported.

A log warning is emitted when the directional maximum is used; where the
all-pairs estimate is available it names both the ratio and that
estimate. A log note is emitted instead when the directional ranges vary
but the spread is consistent with sampling noise.

The returned range is in the coordinate units of the (projected) data
and can be passed directly to `make_folds(block_size = ...)` so that CV
blocks are at least as wide as the autocorrelation range.

## Detrending and the residual-variogram bias

Fitting a variogram to the residuals of a least-squares trend
underestimates both the sill and the range, because the trend fit
absorbs part of the long-wavelength spatial variation (Lark, Cullis and
Welham 2006). Blocks sized from that range are then too small and a
blocked validation is less conservative than it claims. How large the
effect is depends on how smooth the trend terms are in space, and on how
many there are. Measured for this estimator on simulated exponential
fields (n = 300, true effective range 300, nugget 0.2, 40–60 draws), as
the median ratio of the estimate from the trend-removed data to the
estimate from the true field:

- a white-noise covariate: OLS 1.00, REML 0.99 (no bias to speak of);

- a spatially smooth covariate (a random field with range 300 or 1000):
  OLS 0.97, REML 0.95–0.98;

- a linear trend in the coordinates: OLS 0.92, REML 1.04;

- a quadratic trend in the coordinates (five terms): OLS 0.75, REML
  1.06.

So for ordinary covariates the OLS bias is a few percent, and for trend
surfaces in the coordinates it is large. Iterating between a GLS trend
fit and a variogram refit (Neuman and Jacobson 1984) recovers only part
of it (0.80 in the quadratic case), because the variogram of GLS
residuals is biased too; REML does not fit a variogram to residuals at
all, which is why `detrend = "reml"` returns its own range estimate. Its
price is a single family (exponential with nugget), a cubic cost in `n`,
and a somewhat wider sampling spread. The default stays `"ols"` so that
existing scripts return what they did; a script that detrends on smooth
or coordinate-based terms should pass `detrend = "reml"`, or size its
blocks with a margin.

## Count and other non-Gaussian responses

The empirical variogram assumes second-order stationarity: a variance
that is the same everywhere, so that semivariance depends on separation
alone. A count response breaks that assumption by construction, because
its variance tracks its mean, and so does any response with a
mean-variance relationship (rates, proportions, skewed amounts). On such
data the variogram mixes distance-dependence with mean-dependence, and
the fitted range is not the thing it claims to be: where the mean is
high the semivariance is inflated whatever the distance, which flattens
the curve and can lengthen or shorten the apparent range depending on
where the high-mean regions sit.

Supplying `predictor_vars` to detrend helps, because it removes the part
of the mean the covariates explain, but it does not fix the variance
structure: the residuals of an OLS fit to a count still have a variance
that tracks the fitted mean. The principled remedy is a variogram of
Pearson residuals from a model in the right family, which this function
does not compute. Until it does, treat the range from a count response
as no more than an order of magnitude, size blocks conservatively from
it, and prefer
[`make_folds`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)`(method = "nndm")`,
which does not depend on a fitted range at all.

## References

Lark, R. M., Cullis, B. R. and Welham, S. J. (2006). On spatial
prediction of soil properties in the presence of a spatial trend: the
empirical best linear unbiased predictor (E-BLUP) with REML. *European
Journal of Soil Science*, 57(6), 787–799.
[doi:10.1111/j.1365-2389.2005.00768.x](https://doi.org/10.1111/j.1365-2389.2005.00768.x)

Neuman, S. P. and Jacobson, E. A. (1984). Analysis of nonintrinsic
spatial variability by residual kriging with application to regional
groundwater levels. *Mathematical Geology*, 16(5), 499–521.
[doi:10.1007/BF01886329](https://doi.org/10.1007/BF01886329)

## See also

[`sac_nugget`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md)
for the nugget behind the estimate,
[`plot.sac_range`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.sac_range.md)
to see the variogram the estimate rests on.

Other cross-validation:
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md),
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
[`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md),
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md),
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
[`fold_separation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fold_separation.md),
[`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md),
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md),
[`sac_nugget()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md),
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)

## Examples

``` r
if (requireNamespace("gstat", quietly = TRUE)) {
  library(sf)
  # A Gaussian random field with an exponential covariance: range
  # parameter 100, so the true effective range is 3 x 100 = 300, plus a
  # small nugget.
  set.seed(9)
  n <- 150
  xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
  D  <- as.matrix(dist(xy))
  xy$z <- as.numeric(t(chol(exp(-D / 100) + diag(0.1, n))) %*% rnorm(n))
  pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
  r <- estimate_sac_range(pts, response_var = "z")
  # print() throughout, because only the last value of a braced block is
  # shown on its own, and every one of these is worth reading.
  print(r)                              # the effective range, in metres
  print(attr(r, "directional"))         # the four directional ranges
  print(attr(r, "variogram_model"))     # the fitted gstat model behind it

  # A field whose range the data cannot pin down: the variogram never
  # reaches a sill within the lags fitted, so the answer is NA with the
  # refused value attached rather than a long range asserted.
  xy$trend <- sin(xy$x / 400) + rnorm(n, sd = 0.2)
  r2 <- estimate_sac_range(st_as_sf(xy, coords = c("x", "y"), crs = 32632),
                           response_var = "trend")
  print(r2)
  attr(r2, "rejected_range")
}
#> 304.4623 
#>   directional: 0 deg = 267.5557, 45 deg = 439.1852, 90 deg = 370.9347, 135 deg = 328.2169  (ratio 1.64)
#>        0       45       90      135 
#> 267.5557 439.1852 370.9347 328.2169 
#>   model    psill    range
#> 1   Nug 0.000000   0.0000
#> 2   Exp 1.286699 101.4874
#> NA 
#>   directional: 0 deg = 13190.44 (not converged), 45 deg = 42407.17 (not converged), 90 deg = 114115.2 (not converged), 135 deg = 45807.9 (not converged)
#> [1] 47261.27
```
