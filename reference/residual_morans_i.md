# Test the residuals of a fitted model for spatial autocorrelation

Takes the residuals and the observation coordinates of a `spatial_fit`,
builds a spatial weight matrix, and computes Moran's I with its
expectation and variance under a stated null. A z-score and a p-value
say whether spatial autocorrelation remains after fitting. The default
weights are the row-standardised k-nearest-neighbour matrix described
under `weights`.

## Usage

``` r
residual_morans_i(
  fit,
  alternative = c("two.sided", "greater", "less"),
  weights = NULL,
  k = 8L,
  null = c("auto", "randomisation", "residual"),
  keep_weights = FALSE
)
```

## Arguments

- fit:

  A `spatial_fit` object (from `fit_gwr_model`,
  `fit_bayesian_spatial_model` or `fit_rf_model`). The locations are
  read from its `data_sf`, a feature that is not a point reduced to one
  as
  [`coerce_to_points`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md)`(, "auto")`
  does. A layer whose geometries mix XY with XYZ or XYM, or carry M (XYM
  or XYZM), is read in x and y, with no warning (see
  [spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)).

- alternative:

  Character: `"two.sided"` (default), `"greater"` (positive
  autocorrelation), or `"less"`.

- weights:

  Optional user-supplied n x n weight matrix: a base matrix or a
  Matrix-package matrix (e.g. a sparse dgCMatrix). When `NULL` (the
  default), a row-standardised k-nearest-neighbour weight matrix (k = 8)
  is built from the observation coordinates. Ties at the k-th distance
  (every regular grid, every site with repeat visits) share that slot's
  weight equally, with no tie-break by row order or by which backend
  found them, so the matrix is a function of the geometry alone. A
  distance is tied with the k-th when their squares are within one part
  in ten million of each other, so the matrix is the same in any unit
  and at any origin where the doubles hold the coordinates to better
  than about one part in 1e8 of the distances between neighbours.
  Farther from the origin the doubles round the coordinates more
  coarsely, and a lattice whose coordinates are not on that rounding can
  lose its ties. Where a distance lies past a tie by less than that
  rounding, the call warns, with class `"spatialkit_xy_imprecise"`,
  naming the rounding and the distance; the matrix is the one the
  coordinates as stored give. On distinct, untied coordinates it equals
  spdep's `knearneigh()` + `nb2listw(style = "W")` exactly. If a
  non-row-standardised matrix is supplied (i.e. rows do not all sum to
  1), the Cliff & Ord variance formula is still valid for general W and
  the computation proceeds; a note is logged (not raised as an R
  warning) because the magnitude of I is not directly comparable to
  results obtained with row-standardised weights. A `weights` argument
  that is not an n x n matrix is an error.

- k:

  Integer number of nearest neighbours used when building the default
  weight matrix (ignored when `weights` is supplied). Default 8.
  `k >= n - 1` is a complete graph (\\n(n-1)\\ weights however they are
  stored) and is refused above 20 million of them.

- null:

  Which null distribution the expectation, variance and p-value are
  computed against. One of:

  `"auto"` (default)

  :   `"residual"` when the design matrix can be rebuilt *and* the fit's
      residuals are the OLS residuals on it, `"randomisation"`
      otherwise.

  `"randomisation"`

  :   Always the exchangeable moments.

  `"residual"`

  :   Always the Cliff & Ord residual moments. Falls back to
      `"randomisation"` with a logged warning if the design cannot be
      rebuilt, and logs a warning (but proceeds) if the residuals are
      not the OLS residuals on it, in which case the moments are
      approximate.

  Both `"auto"` and `"residual"` also fall back to `"randomisation"`
  when the residual degrees of freedom \\n - p\\ are below 4 (the
  residual variance divides by \\(n-p)(n-p+2)\\ and the normal
  approximation means nothing there); `"residual"` logs a warning when
  it does, `"auto"` does not, and `df` in the result is then \\n - 1\\.
  Read `null` in the result to see which it was. See the section "Which
  null, and when it is approximate".

- keep_weights:

  Logical. Return the \\n \times n\\ weight matrix in the result?
  Defaults to `FALSE`: the matrix makes up most of the object's size,
  while most uses read only the statistic and its moments, and
  `weights_summary` says what it was. Set `TRUE` when you need the
  matrix itself.

## Value

A list with components:

- observed:

  Numeric scalar, Moran's I statistic.

- expected:

  Expected I under the null named by `null`: \\-1/(n-1)\\ for
  `"randomisation"`, \\(n/S_0)\mathrm{tr}(MW)/(n-p)\\ for `"residual"`.

- sd:

  Standard deviation of I under that same null.

- z:

  Standardised z-score, \\(I - E\[I\]) / sd(I)\\.

- p_value:

  Two-sided (or one-sided) p-value from the normal approximation.

- n:

  Number of observations used.

- null:

  The null actually used, `"randomisation"` or `"residual"`. Check it,
  since `"auto"` chooses per fit and `"residual"` can fall back.

- df:

  Residual degrees of freedom behind the moments: \\n - p\\ for
  `"residual"` (where \\p\\ is the rank of the design matrix), \\n - 1\\
  for `"randomisation"`.

- weights:

  The \\n \times n\\ weight matrix the statistic was computed with,
  either the row-standardised k-nearest-neighbour matrix built here
  (sparse when Matrix is installed) or the supplied `weights` after the
  diagonal was zeroed. It is `NULL` unless `keep_weights = TRUE`.

- kurtosis:

  The residual kurtosis \\m_4 / m_2^2\\, which the randomisation
  variance conditions on and which says how far the residuals sit from
  the Gaussian case the residual moments assume (3 for a normal sample).

- p:

  The rank of the rebuilt design matrix behind the `"residual"` moments;
  `NA` under `"randomisation"`.

- exact:

  Logical: `TRUE` when the moments are exact for these residuals (the
  residual null on OLS residuals of the response on the rebuilt design),
  `FALSE` when they are an approximation (the residual null forced onto
  a non-OLS backend, or the randomisation null, whose exchangeable
  moments model residuals do not satisfy).

- weights_summary:

  What the weight matrix was, present whether or not the matrix itself
  was kept: `n`, `storage` (its class), `neighbours` (the smallest and
  largest number of neighbours any row has, and `NA` for a dense matrix,
  where counting them would allocate a second one), `kept` (whether
  `weights` holds the matrix) and `desc`, the one line
  [`print()`](https://rdrr.io/r/base/print.html) shows.

The list is classed `"morans_i"` and has a
[`print()`](https://rdrr.io/r/base/print.html) method, so the console
shows the statistic and its null without printing the \\n \times n\\
`weights` matrix; `[` drops the class, and `$`, `[[` and
[`unlist()`](https://rdrr.io/r/base/unlist.html) are unaffected. A
custom fit whose class has no
[`residuals()`](https://rdrr.io/r/stats/residuals.html) method (which
[`new_spatial_fit`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md)
calls optional) is scored on the observed response minus
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html), as
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) does for it.
Returns `NULL` with a warning if computation fails, saying why:
[`residuals()`](https://rdrr.io/r/stats/residuals.html) raised an error
(its message is quoted), or returned `NULL` and the response minus
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) could not be
formed either (the reason is quoted), fewer than 4 valid residuals, a
residual vector whose length does not match the fit's `data_sf`, weights
that sum to less than `.Machine$double.eps`, or residuals that leave I a
ratio of rounding: they do not vary, or vary by no more than the
rounding of the response (a standard deviation of 3 units in the last
place of its largest value, or a root mean square of \\10\epsilon\\
times that of the response about its centre), or are the OLS residuals
of an exact fit of the response on the predictors (a root mean square of
\\10\epsilon\\ times that of the fit's terms). The test reads the
residuals against the response, so residuals in any unit are scored or
refused alike, and I, z and the p-value do not depend on the unit.

## Which null, and when it is approximate

Two nulls are available, and the one used is reported in the `null`
element of the return value.

`"randomisation"` is the classical exchangeable null: \\E\[I\] =
-1/(n-1)\\ with the Cliff and Ord randomisation variance, conditioning
on the observed kurtosis. These are the moments of I for a vector whose
elements are equally likely in any order.

Model residuals are not exchangeable. They are orthogonal to the design
matrix, which pushes \\E\[I\]\\ below \\-1/(n-1)\\ whenever the
covariates are spatially smooth, and further the more covariates there
are. Under the exchangeable null the z-score of OLS residuals is then
shifted below 0, and the test loses power.

`"residual"` therefore uses the Cliff and Ord (1981, section 8.3)
moments for regression residuals, with \\M = I - X(X'X)^{-1}X'\\ rebuilt
from `predictor_vars` and `data_sf`: \$\$E\[I\] =
(n/S_0)\\\mathrm{tr}(MW)/(n-p)\$\$ \$\$Var\[I\] =
(n/S_0)^2\[\mathrm{tr}(MWMW') + \mathrm{tr}((MW)^2) +
(\mathrm{tr}MW)^2\]/\[(n-p)(n-p+2)\] - E\[I\]^2\$\$ These assume normal
errors and do not condition on the observed kurtosis. They agree with
[`spdep::lm.morantest()`](https://r-spatial.github.io/spdep/reference/lm.morantest.html)
to machine precision.

The residual moments are exact for \\e = My\\ and for nothing else, so
`null = "auto"` does not guess from the fit's class: it rebuilds `X`,
regresses the response on it, and uses the residual moments only when
the supplied residuals are those OLS residuals to numerical tolerance. A
GWR wide enough to have collapsed to global OLS passes that test; the
same GWR at a working bandwidth does not. The response, and each
predictor whose centre is more than 1000 times its spread, are taken
about their centres before the regression, which leaves its residuals as
they are, so a timestamp or a response at 1e10 gives the null, moments
and rank \\p\\ of the same columns at 0. The supplied residuals pass
where none differs from the OLS residuals by more than 1e-7 of the
response's standard deviation, or where the differences are no more than
the rounding a regression of the raw columns leaves: a root mean square
of \\10\epsilon\\ (\\\epsilon\\ = `.Machine$double.eps`) times that of
the fit's terms, among them the offset the fit regressed the response at
(for a GWR that took the response about its centre, the offset left
after that centre). They also pass where none differs by more than 1e-7
of the response's standard deviation from the residuals
[`lm()`](https://rdrr.io/r/stats/lm.html) gives on the raw columns,
where the raw design keeps the rank of the centred one: that rounding
grows with the number of rows.

For the flexible backends neither null is exact, and `"auto"` leaves
them on `"randomisation"`: forcing the OLS moments on them does not
correct them, and for a forest it raises the false-positive rate. A
random forest's residuals here are the out-of-bag ones
([spatial_fit-residuals](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-residuals.md)).
The random forest is anticonservative under both: its residuals are not
a linear projection of the response and they are spatially
heteroscedastic, so neither set of moments describes their null
distribution and the variance is understated whichever is used. GWR is
conservative under both, because it removes far more structure than a
rank-\\p\\ projection does. Treat the p-value from those backends as a
rough indicator, and prefer spatially blocked cross-validated residuals
or an explicit spatial covariance model when the answer has to carry
weight.

A permutation null does not help: permuting the residual vector destroys
the orthogonality that causes the bias, so its mean is the exchangeable
\\-1/(n-1)\\ by construction and it reproduces the randomisation null.

## Reading a large I

A large positive I with a small p-value means the residuals keep spatial
structure the model did not capture (for a forest, read the p-value as
rough; see above): a covariate it lacks, or a scale it cannot represent.
Add the missing covariate if there is one, compare a backend that models
space
([`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
or
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)),
and score the model with blocked folds
([`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)).
A z near 0 with a large p-value means the test found no structure left.

## References

Cliff, A. D. and Ord, J. K. (1981) *Spatial Processes: Models and
Applications*. Pion, London. Section 8.3.

## See also

[`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md),
which adds the statistic to a comparison of fits.

Other in-sample checks:
[`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md),
[`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md),
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md),
[`print.morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.morans_i.md)

## Examples

``` r
# Works on any spatial_fit; a forest keeps the example free of the
# optional GWR and Stan backends.
library(sf)
set.seed(1)
n <- 120
dat <- st_as_sf(
  data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
             elev = rnorm(n)),
  coords = c("x", "y"), crs = 32632
)
# An east-west trend the predictor cannot explain.
dat$price <- 10 + 0.02 * (st_coordinates(dat)[, 1] - 5e5) +
  2 * dat$elev + rnorm(n)
fit <- fit_rf_model(dat, "price", "elev", num_trees = 100, seed = 1)
# A large positive I and z: the residuals carry the east-west trend the
# forest could not see.  A z near 0 with a large p-value would mean no
# structure was left.
residual_morans_i(fit)
#> Residual Moran's I = 0.6437   (E[I] = -0.0084, sd = 0.0419)
#>   z = 15.548, p = 1.653e-54
#>   null: randomisation moments, approximate for these residuals; n = 120, df = 119
#>   residual kurtosis 2.371 (3 = Gaussian)
#>   weights: 120 x 120 dgCMatrix, 8 neighbour(s) per row; not retained, keep_weights = TRUE to keep it
```
