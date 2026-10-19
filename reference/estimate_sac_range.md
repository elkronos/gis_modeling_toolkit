# Estimate the spatial autocorrelation range from data

Fits an exponential variogram model with a nugget and returns the
*effective range*: three times the fitted range parameter, which is
where the semivariance reaches about 95 percent of the sill. It is the
distance beyond which two observations are (near) uncorrelated, which is
what a block or a buffer has to exceed.

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
  keep_directional_fits = FALSE,
  lags = c("adaptive", "equal"),
  n_lags = 15L,
  structure_alpha = 0.05
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

  Character(1) name of the response column. Rows where it is missing or
  infinite are left out before anything else is done, so they neither
  take part in the choice of projection for lon/lat input nor use up the
  `n_max` subsample. For a count or other response whose variance tracks
  its mean, see the section on non-Gaussian responses: the range is
  still estimated, but it is a less reliable number than for a Gaussian
  response.

- predictor_vars:

  Optional character vector. When supplied, the trend on these
  predictors is removed first and the variogram describes the residual
  autocorrelation, the part a spatial model has to handle once the
  covariates have done their work. How the trend is removed is set by
  `detrend`, and it matters: see "Detrending and the residual-variogram
  bias". Rows whose response or predictor is missing or infinite are
  left out of that fit, and so of the variogram, with a logged count.
  When that leaves fewer than 30 rows the result is the unclassed `NA`
  that says how many were left (see Value). One case differs, as it
  always has: when no row at all has every predictor (one is missing on
  every row) no trend can be fitted, and the variogram is of the raw
  response, with an R warning and `detrended` set to `FALSE`.

- n_max:

  Maximum number of points the variogram is computed on, at least 30;
  `Inf` for no limit. Variogram estimation is O(n²) so this keeps
  runtime bounded. A layer with more rows that carry a value is cut to a
  seeded random subsample of this size, drawn after rows with nothing to
  model have been left out and with the rows in a fixed order, so the
  subsample does not depend on how the layer is sorted. The counts are
  returned as the attributes `n_used`, `n_missing` and `n_subsampled`.

- cutoff:

  Fraction of the maximum inter-point distance to use as the variogram
  lag cutoff, above 0 and at most 1. That distance is the farthest pair,
  found on the convex hull, and not the bounding-box diagonal, which
  depends on how the axes are oriented. Default 0.5.

- range_frac:

  Positive numeric. A fitted range exceeding
  `range_frac * cutoff * max_dist` (that is, beyond the longest lag the
  empirical variogram was actually fitted over) is treated as
  unidentified and `NA_real_` is returned. A model fit yields a finite
  number even when the variogram never reaches a sill, and such a value
  extrapolates past the observed lags instead of measuring a long
  autocorrelation range. Passing it to `make_folds(auto_range = TRUE)`
  would collapse the block grid to a single block. Default 1.0; raise it
  to accept ranges extrapolated beyond the fitted lags. The bound does
  not guarantee room for two blocks: at the defaults it is half the
  farthest-pair distance, about 0.71 of the side of a square layer, and
  a block grid needs a range below half the width of the bounding box in
  one direction or the other. An accepted range between the two leaves
  `make_folds(auto_range = TRUE)` room for a single block of that size
  (see its `auto_range` argument for what it does then). Lowering
  `range_frac` to fit the grid would turn those estimates into `NA` and
  the blocks into geometric ones smaller than the range.

- seed:

  RNG seed for the `n_max` subsample, restored afterwards so the
  caller's random stream is untouched. Default `123L`: the subsample is
  an internal approximation and no part of the answer, and leaving it
  unseeded would make the returned range differ between runs on
  identical input and silently advance the caller's RNG. Pass `NULL` for
  an unseeded subsample, or a different number to check how sensitive
  the estimate is to the subsample. Ignored when no more than `n_max`
  rows carry a value, where nothing is sampled. The `reml_max_n`
  subsample uses the same seed, and so does the draw of one row a
  location that `detrend = "reml"` makes on a layer with repeat visits.

- detrend:

  How the trend on `predictor_vars` is removed; ignored when there are
  none. `"ols"` (default) fits it by ordinary least squares and fits the
  variogram to the residuals, which underestimates the range (see the
  section below). `"reml"` fits the trend and an exponential-plus-nugget
  covariance together by residual maximum likelihood with
  [`nlme::gls()`](https://rdrr.io/pkg/nlme/man/gls.html), returns the
  REML range, and attaches the empirical variogram of the REML residuals
  for inspection. It needs nlme, costs \\O(n^3)\\, and so runs on at
  most `reml_max_n` points; when it does not converge the `"ols"` path
  runs instead with an R warning saying so. The start of the fit is
  chosen from the profile REML likelihood on a grid of ranges and nugget
  proportions, and when more than one start is fitted the fit with the
  largest log-likelihood is kept.

- reml_max_n:

  Positive integer, at least 30. With `detrend = "reml"`, the trend and
  covariance are fitted on a seeded random subsample of this many points
  when the layer has more; the fitted trend is then removed from every
  point. The fit takes one row a location: where a location has several
  rows (repeat visits), one of them is drawn at random under `seed`, so
  that the nugget fitted is the variance between visits and not that of
  their minimum or of their mean. Default 400. Raise it for a
  better-determined fit at a cubic cost.

- keep_directional_fits:

  `TRUE` or `FALSE`. Attach each direction's empirical variogram and
  fitted model as `directional_fits`? Defaults to `FALSE`: the four
  variograms are most of the object's size, while the numbers read from
  them (`directional`, `directional_fitted`, `directional_status`,
  `anisotropy`) are attached either way, and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
  effective variogram from its own attribute. Set `TRUE` to inspect the
  directional curves.

- lags:

  How the lags are binned: `"adaptive"` (default) or `"equal"`; see "The
  lags" above.

- n_lags:

  A whole number from 5 to 100, default 15. With `lags = "equal"`, the
  number of bins. With `lags = "adaptive"`, the most bins spaced
  geometrically past the first; fewer are used when that many would hold
  under 30 pairs each.

- structure_alpha:

  The level of the test for spatial structure, above 0 and at most 1;
  default 0.05. A layer whose adjusted p-value is not below it gets no
  range. `1` switches the refusal off; the test is still computed and
  attached.

## Value

A single number, of class `sac_range` in the first two of the three
shapes below and an unclassed `NA` in the third; all three behave as an
ordinary number. The shapes carry different attributes:

- Success:

  A positive effective range in projected coordinate units, with the fit
  attached as attributes `directional` (the 0°, 45°, 90° and 135°
  ranges, named by azimuth; `NA` where that direction's fit was
  unusable), `anisotropy` (largest over smallest), `anisotropy_used`
  (always `FALSE`: the all-pairs fit is the estimate; kept for code that
  reads it), `directional_status` (per azimuth, why a direction is `NA`
  in `directional`: `"ok"`, `"over_cutoff"` (its range ran past the
  largest lag fitted), `"not_converged"` (its fit has no minimum inside
  the ranges searched) or `"no_fit"` (a flat variogram, or too few
  bins)), `directional_fitted` (the range each direction's fit reported
  whether or not it was usable, so a refused directional range stays
  recoverable) and `directional_fits` (a list by azimuth of each
  direction's empirical `variogram` and fitted `model`, `NULL` where
  there is none, and `NULL` altogether unless
  `keep_directional_fits = TRUE`), `detrended` (logical: whether the
  variogram is of the residuals on `predictor_vars` or of the raw
  response. A missing predictor is an error, and a failed detrending fit
  warns and falls back to the raw response with this set to `FALSE`),
  `detrend_method` (`"ols"` or `"reml"` when detrended, `NA` otherwise),
  `reml` (with `detrend = "reml"`: a list with `n_used`, `subsampled`,
  `nugget_prop`, `sigma2`, `logLik` (the REML log-likelihood of the fit
  kept), `logLik_gain` (its gain over the same trend with independent
  errors) and `n_fits` (the starts fitted); `NULL` otherwise), `crs`
  (the projected CRS the variogram was fitted in: the unit of the
  range), `max_dist`, `cutoff_dist`, `variogram` (the empirical
  variogram, in the units of the response's variance), `variogram_model`
  (the fitted model, in gstat's layout, or with `detrend = "reml"` one
  built from the REML parameters), `nugget` (that model's nugget
  variance; see
  [`sac_nugget`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md)),
  `nugget_at_bound` (`TRUE` when the nugget is at its lower bound,
  within 1e-4 of the sill, and so is not an estimate), `first_lag` (the
  mean separation in the first bin: the shortest lag resolved, from
  which the nugget is extrapolated), `bins_below_range` (how many bins
  lie inside the fitted range; with fewer than three the split between
  nugget and structure rests on very little), `short_lag_semivariance`,
  `short_lag_pairs` and `short_lag_dist` (the semivariance of
  nearest-neighbour pairs of locations with no model, their number and
  their median separation: an upper bound on the nugget),
  `zero_lag_semivariance` and `zero_lag_pairs` (the semivariance of the
  pairs of rows that share a location, which is the variance between
  visits to one place and a lower bound on the nugget, and their number;
  `NA` and 0 when no location has two rows), `nugget_floor` (the value
  the fitted nugget was not let under: `zero_lag_semivariance` when the
  rows outnumber the locations by 30 or more, else 0, and always 0 with
  `detrend = "reml"`), `structure_p`, `morans_i` and `structure_test`
  (the adjusted p-value of the test for spatial structure, Moran's I on
  the 8 nearest neighbours, and the statistic, z and p for each
  neighbourhood size), `structures` and `short_range` (how many
  structures the model has besides the nugget, and the effective range
  of the short one when there are two, `NA` otherwise), `lags`,
  `lag_boundaries` and `lag_fallback` (the binning behind the result,
  its upper bin edges, and `TRUE` when adaptive lags were asked for and
  the equal-lag result is returned; see "The lags"), and `n_used`,
  `n_missing` and `n_subsampled` (the points the variogram was computed
  on, the rows left out for a missing response or predictor, and the
  rows left out by the `n_max` subsample), so the fit can be inspected
  and need not be taken on trust.

- Rejected range:

  `NA_real_` when no range is identified. It is classed `sac_range` as
  well, so it prints as `NA` without dumping its attributes, and it
  carries every attribute of a success (the evidence for the rejection),
  plus `rejected_range` (the value that was refused, `NA` when none was
  fitted) and `rejected_reason`, one of:

  - `"no spatial structure detected at the nearest-neighbour scale"`:
    the test above. The model attached is a pure nugget at the variance
    of the modelled variable, and nothing fitted is handed on.

  - `"the empirical variogram is flat: no spatial structure at the lags resolved"`:
    no model with a range fits the variogram better than a constant. The
    model attached is that constant, a pure nugget.

  - `"empirical variogram decreases with distance"`: the semivariance
    falls over the shorter lags (a net fall of more than 15 percent of
    the mean semivariance there, weighted by pairs), which is the shape
    of a periodic, hole-effect structure or of a variance that differs
    between a dense cluster and the rest of the layer. Sampling noise in
    the short-lag bins of a small sample can make that fall too.

  - `"fitted range exceeds the largest lag fitted"`: the range is past
    `range_frac * cutoff * max_dist` (see `range_frac`). An unremoved
    trend ends here when the variogram rises past the fitted lags; a
    milder trend only lengthens the fitted range and passes, which is
    what `predictor_vars` is for.

  - `"fitted range is below the shortest lag fitted"`: too few pairs of
    points lie inside the range to identify it, because it is shorter
    than the shortest lag the empirical variogram resolves (`first_lag`)
    or, with `detrend = "reml"`, than the distance within which 30 pairs
    of the points the REML fit used lie, when that is shorter (the REML
    range is fitted to the point pairs, not to the bins). A structure
    that short cannot be told from a nugget. `range_floor` holds the
    distance the refused range fell short of.

  - `"fitted range is non-positive or non-finite"`.

  - `"no variogram model could be fitted (singular fits)"`: fewer than
    three lag bins hold point pairs. `variogram_model` is `NULL` and
    `nugget` is `NA`.

  The reason names one ground for the refusal, in the order above.
  `attr(attr(range, "variogram_model"), "converged")` says one thing
  more about the model: it is `FALSE` when the least-squares criterion
  is still falling at an edge of the ranges searched (a variogram with
  no sill inside four cutoffs), and the nugget and sill of such a model
  are that edge and not fitted values. For a model with two structures
  the flag is that of the equal-lag fit the long range came from: when
  that fit has no minimum, the long structure's range and partial sill
  are the edge of its search and the flag is `FALSE`, although the
  nugget and the short structure were fitted at the adaptive lags. A
  model built from REML parameters carries no such attribute.
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  reads it, and takes no nugget from a model marked `FALSE`. Results
  stored by earlier versions can carry
  `"variogram model did not converge"` as their reason and
  `anisotropy_used = TRUE`; neither is produced any more.

- No fit:

  An unclassed `NA_real_` when estimation could not be attempted at all,
  whose one attribute, `rejected_reason`, says why:
  `"package 'gstat', which the variogram needs, is not installed"`,
  `"<n> points, fewer than the 30 a variogram range is estimated from"`,
  `"<n> point(s) with a finite value to model, fewer than the 30 a variogram range is estimated from"`
  (the count is of the whole layer's rows with a value),
  `"the response is constant"`,
  `"the residuals on predictor_vars are constant: the predictors explain the response exactly"`
  or
  `"the points have no extent (the largest distance between them is zero or could not be computed)"`.
  Nothing is fitted, so none of the other attributes above exist: no
  `variogram`, `variogram_model` or `rejected_range`, which is what
  tells it from a range that was fitted and refused.

Attributes and the class do not affect
[`is.na()`](https://rdrr.io/r/base/NA.html) or
[`is.finite()`](https://rdrr.io/r/base/is.finite.html), so every
downstream guard treats all three the same way it always did.

## Details

The result is a function of the data alone. The same layer gives the
same answer whatever the unit of the response or of the coordinates and
whatever the order of its rows, and a change in the last bits of the
values changes the answer by about as much. Three things make it so. The
rows are sorted by their coordinates as supplied, then by the response
and the predictors, before anything else, so the projection chosen for
lon/lat input, the `n_max` subsample and the pair sums do not depend on
their order. (Rows that share a location are therefore in the order of
their values, and no step treats the first of them as special.) The
variogram is computed on the standardised variable, and a response is
called constant only when its standard deviation is under 1e-12 of its
largest absolute value, which no change of unit alters. And the model is
fitted without iteration: the criterion is gstat's weighted least
squares (weights \\N / h^2\\) on lags divided by the cutoff; for each
trial range the nugget and the partial sill are solved exactly, with
both kept at or above zero; and the range is found by a search over a
logarithmic grid from a quarter of the first lag to four cutoffs,
refined in one dimension. There is no start value, no iteration limit
and no fit that fails to converge.

Only the exponential family is fitted. On fields smoother than
exponential that makes the range long, which is the safe side (blocks
too large, cross-validation pessimistic). Choosing the family by the
smaller weighted sum of squares would correct a spherical field but send
exponential fields low, the direction that leaks.

**Spatial structure is tested first.** A variogram model can be fitted
to anything, and on a layer with no structure the fit returns a finite
range for a large share of draws. So Moran's I of the modelled variable
is computed on its 2, 4, 8 and 16 nearest neighbours, each standardised
with the moments that fit the variable (those of regression residuals
when `predictor_vars` are given), and the smallest one-sided p-value is
multiplied by four. When that adjusted p-value is not below
`structure_alpha`, the result is `NA` with `rejected_reason`
`"no spatial structure detected at the nearest-neighbour scale"`, a
pure-nugget model and no fitted range. The test uses no random numbers.
`structure_p`, `morans_i` and `structure_test` are attached to every
classed result. Measured on a 1000 m square with 300 points, the test
refuses 388 of 400 white-noise layers and 8 get a finite range, against
196 with the test off. Its power is limited where few neighbours lie
inside the range. With a nugget of a fifth of the sill it finds
structure in 88 percent of fields of effective range 60 on 300 points,
and in 21 percent at range 30, which is about the spacing of those
points. A short range on a sparse layer can therefore be refused;
`structure_alpha = 1` switches the refusal off. Rows that share a
location (repeat visits) are one place and not each other's neighbours,
so on such a layer the test is made on the means of the locations: what
it asks is whether locations near each other are alike, which a
location's own effect, however large, does not answer. With very unequal
numbers of visits the p-value is approximate.

**The lags.** With `lags = "adaptive"` (the default) the first bin runs
from zero to the distance within which the closest `max(50, n / 4)`
pairs of distinct locations lie, and up to `n_lags` further bins are
spaced geometrically from there to the cutoff; a bin holding fewer than
30 pairs is merged with the next. The short lags are then resolved at
the spacing of the points, which is what the nugget and a short range
are read from. With `lags = "equal"` there are `n_lags` bins of equal
width out to the cutoff, as gstat bins by default; each is
`cutoff * max_dist / n_lags` wide, a structure shorter than the first of
them is invisible and is counted as nugget, and the range fitted is that
of whatever longer structure the layer has. Below about a hundred points
neither binning resolves much and the adaptive bins bring no gain.

Two rows at one location are a pair at zero separation, and are in no
bin, whichever the binning: the variogram and the model describe
positive separations. Their semivariance, the variance between visits to
a location, is returned on its own as `zero_lag_semivariance`, and the
nugget fitted to the variogram is not let under it when the layer has 30
or more rows beyond one a location to estimate it from (a model from
`detrend = "reml"` takes no such floor); see
[`sac_nugget`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md).

The adaptive lags put little weight on the long ones, and a weak, long
structure under much noise can be clear at equal lags while the adaptive
fit runs past the largest lag. So when the adaptive range is refused as
past the fitted lags and the equal-lag fit has a range it would itself
return, the result is the equal-lag one, whole, exactly as
`lags = "equal"` gives it, with `lags` saying `"equal"` and
`lag_fallback` set to `TRUE`. Its nugget is then the equal-lag nugget.
No other refusal falls back.

**Structure at two scales.** One structure fitted to lags that resolve
the spacing of the points follows the shortest structure of the layer. A
layer can have a weaker structure several times longer as well, and
blocks sized from the short one would leave the long one across folds.
So with `lags = "adaptive"` the equal-lag fit is made too, and when its
range is at least twice the adaptive one, a model with a nugget and both
structures, at those two ranges, is fitted to the adaptive variogram.
That model is returned (`structures` is 2 and `short_range` holds the
range of the short structure) and the range is the long one, unless the
long structure carries under 3 percent of the model's semivariance at
the largest lag: then the lags that resolve the layer give the equal-lag
range nothing, and the single structure stands. The bar is low because
the two mistakes are not alike. A long structure that is missed sizes
blocks too small, which leaks. One that is believed in error returns the
equal-lag range, which is too long and makes a validation pessimistic,
and `short_range` still holds the short one. When the equal-lag range is
one that `lags = "equal"` would itself refuse (past the fitted lags, as
under an unremoved trend), the model with it is refused as well, the
answer is taken away and not lengthened, and the long structure has to
carry 10 percent.

Measured on a clustered layer of 1200 points over a 1000 m square, 30
draws each. With structures of range 20 and 200 both are kept in 20
draws, and the range returned is under half the equal-lag one in 1; with
40 and 400, in 18 and 1; with 30 and 300, in 21 and 0; with 15 and 450,
in 28 and 1. With one structure of range 40 a second is kept in 4, and
the range is within a factor 2 of the truth in 27 draws (19 with equal
lags); with one of range 30, in 7 and 25 (14). The nugget comes from the
adaptive lags either way: where it is 0.3, the median is 0.32 to 0.35,
and 0.49 to 0.65 with equal lags. A long structure whose equal-lag range
is under twice the adaptive one is not separated: the range is then that
of the single fit, which lies between the two.

The estimate is the **omnidirectional** (all-pairs) fit. Directional
variograms are fitted as well, at 0° (N–S), 45°, 90° (E–W) and 135°
azimuths with a ±22.5° tolerance. Those four windows tile all 180
distinct azimuths exactly once, and their ranges are returned in the
`directional` attribute, with their largest-over-smallest ratio in
`anisotropy`. They are a diagnostic and never the answer, for two
reasons. Each direction sees about a quarter of the point pairs, and the
maximum of four quarter-sample fits is biased upward. And the windows
are fixed to the coordinate axes, so any answer built from them changes
when the layer is rotated, which a property of the field must not do.
The all-pairs fit is the best-powered estimate available and is
invariant to rotation. When it identifies no range, no range is
returned: the largest directional range does not stand in for it.

Where a field is *known* to be anisotropic, blocks must be at least as
large as the longest autocorrelation range to avoid leakage, and the
conservative choice is to size them from the longest directional range
explicitly. Read it from `directional_fitted`, not `directional`: on a
strongly anisotropic field the major axis is the direction most likely
to run past the fitted lags, which leaves it `NA` in `directional`, so
[`max()`](https://rdrr.io/r/base/Extremes.html) of that is `NA`, or with
`na.rm = TRUE` the second-longest range. Check `directional_status`
first: a major axis marked `"over_cutoff"` has no identified range at
all, and a longer `cutoff` or
[`make_folds`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)`(method = "nndm")`
is the way on. A ratio above 1.5 is written to the package log at INFO
level with that advice, which reaches the session log file but not the
console (the ratio passes 1.5 on most isotropic fields too);
[`print()`](https://rdrr.io/r/base/print.html) shows the directional
ranges and the ratio, and `attr(range, "anisotropy")` holds it.

A direction whose variogram is flat or has too few bins, whose range
runs beyond the longest fitted lag, or whose fit has no minimum inside
the ranges searched is excluded and recorded as `NA` in the
`directional` attribute.

The model is fitted **with a nugget**. A nugget-free model forces the
curve through the origin, and on any real measurement (which has one)
the \\N / h^2\\ weights buy that constraint by collapsing the range. The
fitted nugget can sit at its lower bound of zero; the result then
carries `nugget_at_bound = TRUE`, and
[`sac_nugget`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md)
says what else is attached to judge a nugget by.

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
many there are: for ordinary covariates it is a few percent, and for
trend surfaces in the coordinates it is large. Iterating between a GLS
trend fit and a variogram refit (Neuman and Jacobson 1984) recovers only
part of it, because the variogram of GLS residuals is biased too; REML
does not fit a variogram to residuals at all, which is why
`detrend = "reml"` returns its own range estimate. Its price is a single
family (exponential with nugget), a cubic cost in `n`, and a somewhat
wider sampling spread. The default stays `"ols"` so that existing
scripts keep their method; a script that detrends on smooth or
coordinate-based terms should pass `detrend = "reml"`, or size its
blocks with a margin. Detrending removes a trend the predictors explain.
It does not make an unidentified range identified.

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
does not compute. Treat the range from a count response as no more than
an order of magnitude, size blocks conservatively from it, and prefer
[`make_folds`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)`(method = "nndm")`,
which does not depend on a fitted range at all.

## References

Cliff, A. D. and Ord, J. K. (1981). *Spatial Processes: Models and
Applications*. Pion, London.

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
  print(attr(r, "variogram_model"))     # the fitted model behind it
  print(attr(r, "structure_p"))         # the test for spatial structure

  # A field whose range the data cannot pin down: the variogram never
  # reaches a sill within the lags fitted, so the answer is NA with the
  # refused value attached rather than a long range asserted.
  xy$trend <- sin(xy$x / 400) + rnorm(n, sd = 0.2)
  r2 <- estimate_sac_range(st_as_sf(xy, coords = c("x", "y"), crs = 32632),
                           response_var = "trend")
  print(r2)
  attr(r2, "rejected_reason")
  attr(r2, "rejected_range")

  # White noise: no structure, so no range.
  xy$noise <- rnorm(n)
  r3 <- estimate_sac_range(st_as_sf(xy, coords = c("x", "y"), crs = 32632),
                           response_var = "noise")
  attr(r3, "rejected_reason")
}
#> 292.1814 
#>   directional: 0 deg = 248.1608, 45 deg = 413.7318, 90 deg = 617.9676, 135 deg = 426.5880  (ratio 2.49)
#>   in metres of EPSG:32632; variogram of the response itself
#>        0       45       90      135 
#> 248.1608 413.7318 617.9676 426.5880 
#>   model    psill    range
#> 1   Nug 0.000000  0.00000
#> 2   Exp 1.273642 97.39379
#> [1] 8.510444e-22
#> NA 
#>   directional: 0 deg = 2672.119 (past the fitted lags), 45 deg = 2672.119 (past the fitted lags), 90 deg = 2672.119 (past the fitted lags), 135 deg = 2672.119 (past the fitted lags)
#>   in metres of EPSG:32632; variogram of the response itself
#> [1] "no spatial structure detected at the nearest-neighbour scale"
```
