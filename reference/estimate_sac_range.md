# Estimate the spatial autocorrelation range from data

Fits an exponential variogram model with a nugget and returns the
effective range: three times the fitted range parameter, where the
semivariance reaches about 95 percent of the sill. Beyond it two
observations are close to uncorrelated, so it is the distance a block or
a buffer has to exceed. The range is in the units of the projected CRS
the variogram is fitted in. Pass it as the `block_size` of
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md),
or let `make_folds(auto_range = TRUE)` estimate it, so that the blocks
are at least as wide as the autocorrelation range.

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

# S3 method for class 'sac_range'
print(x, ...)
```

## Arguments

- points_sf:

  An sf layer of points. A lon/lat layer is projected first
  ([spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md));
  other geometry is reduced to representative points, and Z and M are
  dropped, since
  [`gstat::variogram()`](https://r-spatial.github.io/gstat/reference/variogram.html)
  would otherwise measure distances in three dimensions. Rows without a
  location are dropped with a logged count.

- response_var:

  Name of the response column. Rows where it is missing or infinite are
  left out before anything else is done, so they take no part in the
  choice of projection for lon/lat input and do not use up the `n_max`
  subsample. For a count or other response whose variance tracks its
  mean, see the section "Count and other non-Gaussian responses": the
  range is still estimated, but it is less reliable than for a Gaussian
  response.

- predictor_vars:

  Optional character vector of predictor column names. When given, the
  trend on these predictors is removed first and the variogram describes
  the residual autocorrelation, the part a spatial model has to handle
  once the covariates have done their work. How the trend is removed is
  set by `detrend` (see the section "Detrending and the
  residual-variogram bias"). Rows whose response or predictor is missing
  or infinite are left out of that fit, and so of the variogram, with a
  logged count; when that leaves fewer than 30 rows, or fewer than 30
  distinct locations, the result is the unclassed `NA` that says how
  many were left (see Value). When no row has every predictor (one is
  missing on every row) no trend can be fitted, and the variogram is of
  the raw response, with an R warning and `detrended` set to `FALSE`.
  Where nothing is fitted then either (fewer than 30 distinct locations,
  a constant response, or points with no extent), the result is the `NA`
  with its reason, and there is no such warning.

- n_max:

  The most points the variogram is computed on, at least 30; `Inf` for
  no limit. Variogram estimation is `O(n^2)`, so this bounds the run
  time. A layer with more rows that carry a value is cut to a seeded
  random subsample of this size, drawn after the rows with nothing to
  model are left out and with the rows in a fixed order, so the
  subsample does not depend on how the layer is sorted. The counts are
  returned as the attributes `n_used`, `n_missing` and `n_subsampled`.
  The floor of 30 distinct locations applies to the rows the variogram
  is computed on, and a subsample of a layer with many rows a location
  can hold fewer locations than the layer. When the layer itself has 30
  locations or more and the subsample holds fewer, the reason of the
  `NA` says `"of the n_max subsample"`, and a larger `n_max` keeps more
  of them; when the layer has fewer than 30, the reason gives the
  layer's counts.

- cutoff:

  Fraction of the largest distance between two points to use as the
  variogram lag cutoff, above 0 and at most 1. That distance is the
  farthest pair, found on the convex hull, not the diagonal of the
  bounding box, which depends on how the axes are oriented. Default 0.5.

- range_frac:

  A positive number, default 1. A fitted range beyond
  `range_frac * cutoff * max_dist` (at the default, beyond the longest
  lag the empirical variogram was fitted over) is treated as
  unidentified and the result is `NA`. A model fit gives a finite number
  even when the variogram never reaches a sill, and such a value
  extrapolates past the observed lags; passed to
  `make_folds(auto_range = TRUE)` it would collapse the block grid to a
  single block. Raise `range_frac` to accept ranges extrapolated beyond
  the fitted lags. The bound does not guarantee room for two blocks: at
  the defaults it is half the farthest-pair distance, about 0.71 of the
  side of a square layer, and a block grid needs a range below half the
  width of the bounding box in one direction or the other. An accepted
  range between the two leaves `make_folds(auto_range = TRUE)` room for
  a single block of that size (see its `auto_range` argument). Lowering
  `range_frac` to fit the grid would turn those estimates into `NA` and
  the blocks into geometric ones smaller than the range.

- seed:

  RNG seed for the `n_max` subsample, the `reml_max_n` subsample, and
  the draw of one row a location that `detrend = "reml"` makes on a
  layer with repeat visits. With a seed, the caller's random stream is
  restored afterwards. Default 123, so the same layer gives the same
  range on every run; `NULL` leaves the draws unseeded, and another
  number shows how sensitive the estimate is to the subsample. It has no
  effect where nothing is drawn. On a layer whose rows crowd onto few
  locations the subsample can also decide whether 30 distinct locations
  remain, and so whether a range is estimated at all (see `n_max`).

- detrend:

  How the trend on `predictor_vars` is removed; ignored when there are
  none. `"ols"` (default) fits it by ordinary least squares and fits the
  variogram to the residuals, which underestimates the range (see the
  section "Detrending and the residual-variogram bias"). `"reml"` fits
  the trend and an exponential-plus-nugget covariance together by
  residual maximum likelihood with
  [`nlme::gls()`](https://rdrr.io/pkg/nlme/man/gls.html), returns the
  REML range, and attaches the empirical variogram of the REML residuals
  for inspection. It needs nlme, costs `O(n^3)` and so runs on at most
  `reml_max_n` points, and takes one row a location. When it is not
  fitted (no nlme), cannot be fitted (the design of the trend cannot be
  built, as on a factor with one level) or does not converge, the
  `"ols"` path runs instead with an R warning that says which. The start
  of the fit is chosen from the profile REML likelihood on a grid of
  ranges and nugget proportions, and when more than one start is fitted
  the fit with the largest log-likelihood is kept.

- reml_max_n:

  A whole number of at least 30, default 400. With `detrend = "reml"`,
  the trend and covariance are fitted on a seeded random subsample of
  this many points when the layer has more, and the fitted trend is then
  removed from every point. The fit takes one row a location: where a
  location has several rows (repeat visits), one of them is drawn at
  random under `seed`, so the nugget fitted is the variance between
  visits and not that of their minimum or of their mean. Raise it for a
  better-determined fit at a cubic cost.

- keep_directional_fits:

  `TRUE` or `FALSE`, default `FALSE`. Attach each direction's empirical
  variogram and fitted model as `directional_fits`. The four variograms
  are most of the object's size; the numbers read from them
  (`directional`, `directional_fitted`, `directional_status`,
  `anisotropy`) are attached either way, and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
  all-pairs variogram from its own attribute.

- lags:

  How the lags are binned: `"adaptive"` (default) or `"equal"`; see the
  section "The lags".

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

- x:

  A `sac_range` object, for
  [`print()`](https://rdrr.io/r/base/print.html).

- ...:

  Ignored.

## Value

A single number, in one of three shapes, each of which behaves as an
ordinary number. A range that was fitted, and a range fitted and
refused, are of class `sac_range`; a call that could fit nothing returns
an unclassed `NA`. The attributes and the class do not affect
[`is.na()`](https://rdrr.io/r/base/NA.html) or
[`is.finite()`](https://rdrr.io/r/base/is.finite.html), so code that
tests for a missing range treats all three alike.

- Success:

  A positive effective range in the units of the projected CRS, with the
  fit attached as attributes:

  - `crs`: the projected CRS the variogram was fitted in, which gives
    the unit of the range; `max_dist` and `cutoff_dist`: the farthest
    pair and the lag cutoff.

  - `variogram`: the empirical variogram, in the units of the response's
    variance; `variogram_model`: the fitted model, in gstat's layout, or
    with `detrend = "reml"` one built from the REML parameters.

  - `lags`, `lag_boundaries` and `lag_fallback`: the binning behind the
    result, its upper bin edges, and `TRUE` when adaptive lags were
    asked for and the equal-lag result is returned (see the section "The
    lags").

  - `first_lag`: the mean separation in the first bin, the shortest lag
    resolved, from which the nugget is extrapolated; `bins_below_range`:
    how many bins lie inside the fitted range (with fewer than three,
    the split between nugget and structure rests on very little).

  - `structures` and `short_range`: how many structures the model has
    besides the nugget, and the effective range of the short one when
    there are two (`NA` otherwise).

  - `nugget` (the model's nugget variance), `nugget_at_bound`,
    `nugget_floor`, `short_lag_semivariance`, `short_lag_pairs`,
    `short_lag_dist`, `zero_lag_semivariance` and `zero_lag_pairs`: the
    nugget and the evidence to judge it by, described in
    [`sac_nugget()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md).
    `zero_lag_semivariance` is `NA` and `zero_lag_pairs` 0 when no
    location has two rows.

  - `structure_p`, `morans_i` and `structure_test`: the adjusted p-value
    of the test for spatial structure, Moran's I on the 8 nearest
    neighbours, and the statistic, z and p for each neighbourhood size.

  - `directional`: the ranges at azimuths 0, 45, 90 and 135 degrees,
    named by azimuth, `NA` where that direction's fit was unusable;
    `directional_status`: per azimuth, why a direction is `NA` in
    `directional`: `"ok"`, `"over_cutoff"` (its range ran past the
    largest lag fitted), `"not_converged"` (its fit has no minimum
    inside the ranges searched) or `"no_fit"` (a flat variogram, or too
    few bins); `directional_fitted`: the range each direction's fit
    reported, usable or not, so a refused directional range stays
    recoverable; `directional_fits`: with
    `keep_directional_fits = TRUE`, a list by azimuth of each
    direction's empirical `variogram` and fitted `model` (`NULL` where
    there is none), and `NULL` otherwise; `anisotropy`: the largest
    directional range over the smallest; `anisotropy_used`: always
    `FALSE`, since the all-pairs fit is the estimate.

  - `detrended`: whether the variogram is of the residuals on
    `predictor_vars` (`TRUE`) or of the raw response. A predictor that
    is not a column is an error, and a detrending fit that fails warns
    and falls back to the raw response with this set to `FALSE`.
    `detrend_method`: `"ols"` or `"reml"` when detrended, `NA`
    otherwise. `reml`: with `detrend = "reml"`, a list with `n_used`,
    `subsampled`, `nugget_prop`, `sigma2`, `logLik` (the REML
    log-likelihood of the fit kept), `logLik_gain` (its gain over the
    same trend with independent errors) and `n_fits` (the starts
    fitted); `NULL` otherwise.

  - `n_used`, `n_missing` and `n_subsampled`: the points the variogram
    was computed on, the rows left out for a missing response or
    predictor, and the rows the `n_max` subsample left out.

- Rejected range:

  `NA_real_` when no range is identified. It is classed `sac_range` as
  well, so it prints as `NA` without dumping its attributes, and it
  carries every attribute of a success (the evidence for the rejection),
  plus `rejected_range` (the value that was refused, `NA` when none was
  fitted) and `rejected_reason`, one of:

  - `"no spatial structure detected at the nearest-neighbour scale"`:
    the test in the section "The test for spatial structure". The model
    attached is a pure nugget at the variance of the modelled variable,
    and nothing fitted is handed on.

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

  - `"no variogram model could be fitted: fewer than three lag bins hold point pairs"`.
    `variogram_model` is then `NULL` and `nugget` is `NA`.

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
  reads it, and takes no nugget from a model marked `FALSE`.

- No fit:

  An unclassed `NA_real_` when estimation could not be attempted at all,
  whose one attribute, `rejected_reason`, says why:
  `"package 'gstat', which the variogram needs, is not installed"`,
  `"<n> points, fewer than the 30 a variogram range is estimated from"`,
  `"<n> point(s) with a finite value to model, fewer than the 30 a variogram range is estimated from"`
  (the count is of the whole layer's rows with a value),
  `"<k> distinct location(s) among the <n> point(s) with a finite value to model, fewer than the 30 a variogram range is estimated from"`
  (rows that share a location count once; when the layer has 30 or more
  locations and the `n_max` subsample kept fewer, the counts are the
  subsample's and `"of the n_max subsample"` stands in place of
  `"with a finite value to model"`), `"the response is constant"`,
  `"the residuals on predictor_vars are constant: the predictors explain the response exactly"`
  or
  `"the points have no extent (the largest distance between them is zero or could not be computed)"`.
  Nothing is fitted, so none of the other attributes exist: no
  `variogram`, `variogram_model` or `rejected_range`, which is what
  tells it from a range that was fitted and refused.

## Details

The layer is first tested for spatial structure, and a range the data
cannot identify is refused: the answer is then `NA`, with the reason and
the evidence attached (see Value).
[`plot.sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.sac_range.md)
draws the variogram behind the estimate, and
[`sac_nugget()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md)
reads its nugget.

## What is estimated

The estimate is the omnidirectional (all-pairs) fit of an exponential
model with a nugget. Only the exponential family is fitted. On fields
smoother than exponential that makes the range long, and blocks sized
from it come out larger: they still hold the test points at least the
field's range from their training data, which blocks sized from a short
range do not. Choosing the family by the smaller weighted sum of squares
would correct a spherical field but send exponential fields short, and
blocks sized from them smaller than the range.

The model is fitted without iteration. The criterion is gstat's weighted
least squares (weights `N / h^2`) on lags divided by the cutoff. For
each trial range the nugget and the partial sill are solved exactly,
both kept at or above zero, and the range is found by a search over a
logarithmic grid from a quarter of the first lag to four cutoffs,
refined in one dimension. There is no start value, no iteration limit
and no fit that fails to converge.

## The test for spatial structure

A variogram model can be fitted to anything, and on a layer with no
structure the fit returns a finite range for a large share of draws. So
Moran's I of the modelled variable is computed first, on its 2, 4, 8 and
16 nearest neighbours, each standardised with the moments that fit the
variable (those of regression residuals when `predictor_vars` are
given), and the smallest one-sided p-value is multiplied by four. When
that adjusted p-value is not below `structure_alpha`, the result is `NA`
with `rejected_reason`
`"no spatial structure detected at the nearest-neighbour scale"`, a
pure-nugget model and no fitted range. The test uses no random numbers,
and `structure_p`, `morans_i` and `structure_test` are attached to every
classed result. Its power is limited where few neighbours lie inside the
range, so a short range on a sparse layer can be refused;
`structure_alpha = 1` switches the refusal off.

Neighbours tied at the k-th distance, those whose squared distances are
within one part in ten million of the k-th's, share its slot's weight
equally, as in
[`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md).
The graphs then do not depend on the unit or the origin of the
coordinates wherever the doubles hold them to better than about one part
in 1e8 of the distances between neighbours. Where they do not, and a
distance lies past a tie by less than that rounding, the call warns with
class `"spatialkit_xy_imprecise"` (see
[spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)),
as it does for the nearest other location of the short-lag
semivariances. The results are those the coordinates as stored give.

Rows that share a location (repeat visits) are one place and not each
other's neighbours, so on such a layer the test is made on the means of
the locations: it asks whether locations near each other are alike,
which a location's own effect, however large, does not answer. With very
unequal numbers of visits the p-value is approximate. For the same
reason a layer needs 30 distinct locations, not only 30 rows; below that
the result is the unclassed `NA` whose reason gives the count (see
Value), and nothing is fitted. Five locations visited 20 times give ten
pairs of locations and five means, from which neither a range nor a
verdict on structure can be read. Two rows are at one location when both
their coordinates are equal to the last bit after the projection, so
repeat visits recorded with slightly different coordinates count as
separate locations; round the coordinates first to pool them.

## The lags

With `lags = "adaptive"` (the default) the first bin runs from zero to
the distance within which the closest `max(50, n / 4)` pairs of distinct
locations lie, and up to `n_lags` further bins are spaced geometrically
from there to the cutoff; a bin holding fewer than 30 pairs is merged
with the next. The short lags, which the nugget and a short range are
read from, are then resolved at the spacing of the points. With
`lags = "equal"` there are `n_lags` bins of equal width out to the
cutoff, as gstat bins by default. Each is `cutoff * max_dist / n_lags`
wide; a structure shorter than the first of them is counted as nugget,
and the range fitted is that of whatever longer structure the layer has.
Below about a hundred points neither binning resolves much, and the
adaptive bins bring no gain.

The adaptive lags put little weight on the long ones, and a weak, long
structure under much noise can be clear at equal lags while the adaptive
fit runs past the largest lag. So when the adaptive range is refused as
past the fitted lags and the equal-lag fit has a range it would itself
return, the result is the equal-lag one, whole, as `lags = "equal"`
gives it, with `lags` saying `"equal"` and `lag_fallback` set to `TRUE`.
Its nugget is then the equal-lag nugget. No other refusal falls back.

## Structure at two scales

One structure fitted to lags that resolve the spacing of the points
follows the shortest structure of the layer. A layer can have a weaker
structure several times longer as well, and blocks sized from the short
one would leave the long one across folds. So with `lags = "adaptive"`
the equal-lag fit is made too, and when its range is at least twice the
adaptive one, a model with a nugget and both structures, at those two
ranges, is fitted to the adaptive variogram. That model is returned
(`structures` is 2 and `short_range` holds the range of the short
structure) and the range is the long one, unless the long structure
carries under 3 percent of the model's semivariance at the largest lag;
then the single structure stands.

The bar is low because the two mistakes are not alike. A long structure
that is missed is lost: the range returned is the short one alone. One
that is believed in error returns the equal-lag range, which is too
long, and `short_range` still holds the short one. When the equal-lag
range is one that `lags = "equal"` would itself refuse (past the fitted
lags, as under an unremoved trend), the model with it is refused as
well, the answer is taken away and not lengthened, and the long
structure has to carry 10 percent. The nugget comes from the adaptive
lags either way. A long structure whose equal-lag range is under twice
the adaptive one is not separated: the range is then that of the single
fit, which lies between the two.

## Directional ranges and anisotropy

Directional variograms are fitted as well, at azimuths of 0
(north-south), 45, 90 (east-west) and 135 degrees, each with a tolerance
of 22.5 degrees either way. Those four windows tile all 180 distinct
azimuths once, and their ranges are returned in the `directional`
attribute, with their largest-over-smallest ratio in `anisotropy`. They
are a diagnostic and never the answer, for two reasons. Each direction
sees about a quarter of the point pairs, and the maximum of four
quarter-sample fits is biased upward. And the windows are fixed to the
coordinate axes, so any answer built from them changes when the layer is
rotated, which a property of the field must not do. The all-pairs fit is
the best-powered estimate available and does not change under rotation.
When it identifies no range, no range is returned: the largest
directional range does not stand in for it.

A direction whose variogram is flat or has too few bins, whose range
runs beyond the longest fitted lag, or whose fit has no minimum inside
the ranges searched is left out and recorded as `NA` in `directional`.

Where a field is known to be anisotropic and blocks are to be at least
as large as its autocorrelation range in every direction, size them from
the longest directional range yourself. Read it from
`directional_fitted`, not `directional`: on a strongly anisotropic field
the major axis is the direction most likely to run past the fitted lags,
which leaves it `NA` in `directional`, so
[`max()`](https://rdrr.io/r/base/Extremes.html) of that is `NA`, or with
`na.rm = TRUE` the second-longest range. Check `directional_status`
first: a major axis marked `"over_cutoff"` has no identified range at
all, and a longer `cutoff` or `make_folds(method = "nndm")` is the way
on. A ratio above 1.5 is written to the package log at INFO level with
that advice, which reaches the session log file but not the console (the
ratio passes 1.5 on most isotropic fields too).

## The nugget

The model is fitted with a nugget. A nugget-free model forces the curve
through the origin, and on a measurement that has a nugget the `N / h^2`
weights buy that constraint by collapsing the range. The fitted nugget
can sit at its lower bound of zero; the result then carries
`nugget_at_bound = TRUE`, and
[`sac_nugget()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md)
lists what else is attached to judge a nugget by.

Two rows at one location are a pair at zero separation and are in no
bin, whichever the binning, since the variogram and the model describe
positive separations. Their semivariance, the variance between visits to
a location, is returned on its own as `zero_lag_semivariance`. The
nugget fitted to the variogram is not let under it when the layer has 30
or more rows beyond one a location to estimate it from; a model from
`detrend = "reml"` takes no such floor.

## Detrending and the residual-variogram bias

Fitting a variogram to the residuals of a least-squares trend
underestimates both the sill and the range, because the trend fit
absorbs part of the long-wavelength spatial variation (Lark, Cullis and
Welham 2006). Blocks sized from that range are then narrower than the
field's range, and hold the test points less far from their training
data than range-sized blocks are meant to. How large the effect is
depends on how smooth the trend terms are in space and on how many there
are: it is small for ordinary covariates and large for trend surfaces in
the coordinates. Iterating between a GLS trend fit and a variogram refit
(Neuman and Jacobson 1984) recovers only part of it, because the
variogram of GLS residuals is biased too. REML fits no variogram to
residuals, which is why `detrend = "reml"` returns its own range
estimate. Its price is a single family (exponential with nugget), a
cubic cost in `n`, and a somewhat wider sampling spread. The default is
`"ols"`; a script that detrends on smooth or coordinate-based terms
should pass `detrend = "reml"`, or size its blocks with a margin above
the range returned, to make up the shortfall. Detrending removes a trend
the predictors explain; it does not make an unidentified range
identified.

## Count and other non-Gaussian responses

The empirical variogram assumes second-order stationarity: a variance
that is the same everywhere, so that semivariance depends on separation
alone. A count response breaks that assumption by construction, because
its variance tracks its mean, and so does any response with a
mean-variance relationship (rates, proportions, skewed amounts). On such
data the variogram mixes distance-dependence with mean-dependence: where
the mean is high the semivariance is inflated whatever the distance,
which flattens the curve and can lengthen or shorten the apparent range
depending on where the high-mean areas sit.

Detrending on `predictor_vars` removes the part of the mean the
covariates explain, but it does not fix the variance structure: the
residuals of an OLS fit to a count still have a variance that tracks the
fitted mean. The remedy is a variogram of Pearson residuals from a model
in the right family, which this function does not compute. Treat the
range from a count response as no more than an order of magnitude, since
it can come out long or short, and prefer `make_folds(method = "nndm")`,
which needs no fitted range.

## Numerical behaviour

The result is a function of the data alone. The same layer gives the
same answer whatever the unit of the response or of the coordinates and
whatever the order of its rows, and a change in the last bits of the
values changes the answer by about as much.

The rows are sorted by their coordinates as supplied, then by the
response and the predictors, before anything else, so the projection
chosen for lon/lat input, the `n_max` subsample and the pair sums do not
depend on their order. First coordinates within 1e-7 of the diagonal of
the points' bounding box count as one column, sorted by the second
coordinate, so a lattice whose first coordinates in one column differ in
their last bits after a change of unit or origin gets the same order and
the same subsample in each. Rows that share a location are therefore in
the order of their values, and no step treats the first of them as
special.

The variogram is computed on the standardised variable. A response is
called constant only when its standard deviation is under 1e-12 of its
largest absolute value, which no change of unit alters, and the values
centred on their median and scaled have a standard deviation of at most
three units in the last place of that value. Residuals of the predictors
are called constant when they are no more than the rounding of the
regression, taken on the size of the fit's terms, or three such units of
the response, whether or not the floor finds spread in them, as the test
for `cp` of
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
reads them. So 1e15 plus a spread of 1, stored to an eighth, is not
constant, and the residuals of a response the predictors fit exactly
are, on a predictor far from 0 (a northing) as well. A response that
loses its spread on the raw values is centred on its median before the
predictors are fitted to it and before the test for spatial structure,
neither of which changes under a shift; so is a response whose OLS
residuals on the raw values lose their spread (1e15 plus 1e5 times a
predictor plus a spread of 1, on that predictor).

The test for spatial structure, the semivariances and the REML fit are
taken on the response divided by a power of two near its spread. A
response near 1e-160 or 1e160, whose squares lie past what doubles hold,
then gives the range, the test and the folds it gives in ordinary units.
The variogram, the nugget, the short-lag semivariances and the REML
variance are returned in the response's units squared, which near 1e160
are `Inf` and near 1e-160 are subnormal doubles that keep few bits, or
0. A response so far from zero that doubles there are spaced more than a
millionth of its standard deviation apart (a standard deviation of 1 at
about 1e10, 2^33, or beyond) is read as stored, rounded to that spacing,
with a warning of class `spatialkit_response_imprecise`.

With `detrend = "reml"` the REML fit reads the coordinates of the rows
it uses centred on the middle of their bounding box, divided by its
diagonal and rounded to 2^-30 of it, and the response and each numeric
predictor centred on the middle of their range on those rows, divided by
it and rounded the same way. Its range, coefficients and variance are
taken back to the layer's units. Where the REML range is not identified
the likelihood is flat, and whether the fit converges or falls back to
OLS is set by the last bits of the coordinates, the response and the
predictors; in that frame it does not depend on their units or origins.
A value within its own rounding noise of the middle of a step, about a
billionth of the diagonal or the range, can still fall on either side of
it.

## Printing

[`print()`](https://rdrr.io/r/base/print.html) shows the effective range
as a plain number, and beneath it the directional ranges and their ratio
when the directional fit is available. A direction whose fit was
unusable is labelled with why (`directional_status`) and the range its
fit reported (`directional_fitted`) when the object carries them, and
`unidentified` otherwise. A last line names the unit and the CRS the
range is a length in (`attr(x, "crs")`, which for lon/lat input is the
projected CRS the estimate chose; for a layer with no CRS, it says the
range is in that layer's own coordinate units), and whether the
variogram is of the response or of its residuals on `predictor_vars`
(`detrended`, `detrend_method`).
[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.

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

[`sac_nugget()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md)
for the nugget behind the estimate, and
[`plot.sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.sac_range.md)
to see the variogram the estimate rests on.

Other spatial autocorrelation:
[`plot.sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.sac_range.md),
[`sac_nugget()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md)

## Examples

``` r
library(sf)
# A Gaussian random field with an exponential covariance: range parameter
# 100, so the true effective range is 3 x 100 = 300, plus a small nugget.
set.seed(9)
n <- 150
xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
D  <- as.matrix(dist(xy))
xy$z <- as.numeric(t(chol(exp(-D / 100) + diag(0.1, n))) %*% rnorm(n))
pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
r <- estimate_sac_range(pts, response_var = "z")
r                            # the effective range, in metres
#> 292.1814 
#>   directional: 0 deg = 248.1608, 45 deg = 413.7318, 90 deg = 617.9676,
#>                135 deg = 426.588  (ratio 2.49)
#>   in metres of EPSG:32632; variogram of the response itself
attr(r, "directional")       # the four directional ranges
#>        0       45       90      135 
#> 248.1608 413.7318 617.9676 426.5880 
attr(r, "variogram_model")   # the fitted model behind it
#>   model    psill    range
#> 1   Nug 0.000000  0.00000
#> 2   Exp 1.273642 97.39379
attr(r, "nugget_at_bound")   # TRUE: the fit put the nugget at its bound;
#> [1] TRUE
                             # ?sac_nugget says how to read that
attr(r, "structure_p")       # the test for spatial structure
#> [1] 8.510444e-22

# A field whose range the data cannot pin down: the variogram never
# reaches a sill within the lags fitted, so the answer is NA, with the
# refused value and the reason attached.
xy$trend <- sin(xy$x / 400) + rnorm(n, sd = 0.2)
r2 <- estimate_sac_range(st_as_sf(xy, coords = c("x", "y"), crs = 32632),
                         response_var = "trend")
r2
#> NA 
#>   directional: 0 deg = 2672.119 (past the fitted lags),
#>                45 deg = 2672.119 (past the fitted lags),
#>                90 deg = 2672.119 (past the fitted lags),
#>                135 deg = 2672.119 (past the fitted lags)
#>   in metres of EPSG:32632; variogram of the response itself
attr(r2, "rejected_reason")
#> [1] "fitted range exceeds the largest lag fitted"
attr(r2, "rejected_range")
#> [1] 2672.119

# White noise: no structure, so no range.
xy$noise <- rnorm(n)
r3 <- estimate_sac_range(st_as_sf(xy, coords = c("x", "y"), crs = 32632),
                         response_var = "noise")
attr(r3, "rejected_reason")
#> [1] "no spatial structure detected at the nearest-neighbour scale"
```
