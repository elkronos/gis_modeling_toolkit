# Compute Moran's I on the residuals of a fitted spatial model

Given a `spatial_fit` object, extracts the residuals and the observation
coordinates, builds a spatial weight matrix, and computes Moran's I
together with its analytical expectation and variance under a stated
null. A z-score and two-sided p-value are provided so the caller can
assess whether statistically significant spatial autocorrelation remains
after fitting.

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
  `fit_bayesian_spatial_model` or `fit_rf_model`).

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
  found them, so the matrix is a function of the geometry alone; on
  distinct, untied coordinates it equals `spdep`'s `knearneigh()` +
  `nb2listw(style = "W")` exactly. If a non-row-standardised matrix is
  supplied (i.e. rows do not all sum to 1), the Cliff & Ord variance
  formula is still valid for general W and the computation proceeds; a
  note is logged (not raised as an R warning) because the magnitude of I
  is not directly comparable to results obtained with row-standardised
  weights. A `weights` argument that is not an n x n matrix is an error.

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
      rebuilt, and warns (but proceeds) if the residuals are not the OLS
      residuals on it, in which case the moments are approximate.

  Both `"auto"` and `"residual"` also fall back to `"randomisation"`
  when the residual degrees of freedom \\n - p\\ are below 4 (the
  residual variance divides by \\(n-p)(n-p+2)\\ and the normal
  approximation means nothing there); `"residual"` logs a warning when
  it does, `"auto"` does not, and `df` in the result is then \\n - 1\\.
  Read `null` in the result to see which it was. See **Which null, and
  when it is approximate** above.

- keep_weights:

  Logical. Return the \\n \times n\\ weight matrix in the result?
  Defaults to `FALSE`: the matrix dominates the object's size (50.3 KB
  of a 53.5 KB result at \\n = 500\\ in its sparse form, and 191 MB at
  the \\n = 5000\\ the dense fallback is capped at) while most uses read
  only the statistic and its moments, and `weights_summary` says what it
  was. Set `TRUE` when you need the matrix itself.

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
[`unlist()`](https://rdrr.io/r/base/unlist.html) are unaffected. Returns
`NULL` with a warning if computation fails (e.g. fewer than 4 valid
residuals).

## Details

By default, weights are constructed as a k-nearest-neighbour (k = 8)
binary matrix, row-standardised. Users may supply their own weight
matrix via the `weights` argument.

## Which null, and when it is approximate

Two nulls are available, and the one actually used is reported back in
the `null` element of the return value.

`"randomisation"` is the classical exchangeable null: \\E\[I\] =
-1/(n-1)\\ with the Cliff & Ord randomisation variance, conditioning on
the observed kurtosis. These are the moments of I for a vector whose
elements are equally likely in any order.

**Model residuals are not exchangeable.** They are orthogonal to the
design matrix, which pushes \\E\[I\]\\ materially below \\-1/(n-1)\\
whenever the covariates are spatially smooth, and pushes it further the
more covariates there are. In a simulation with \\n = 120\\, six smooth
covariates and *independent* errors (so the truth is "no residual
autocorrelation"), OLS residuals had mean \\I = -0.031\\ against the
exchangeable \\E\[I\] = -0.008\\; the z-score averaged \\-0.54\\ with
\\sd = 0.90\\ instead of 0 and 1. The cost is power, which is the point
of the test: at a moderate residual autocorrelation the exchangeable
null rejected 13\\

`"residual"` therefore uses the Cliff & Ord (1981) sec. 8.3 moments for
regression residuals, with \\M = I - X(X'X)^{-1}X'\\ rebuilt from
`predictor_vars` and `data_sf`: \$\$E\[I\] =
(n/S_0)\\\mathrm{tr}(MW)/(n-p)\$\$ \$\$Var\[I\] =
(n/S_0)^2\[\mathrm{tr}(MWMW') + \mathrm{tr}((MW)^2) +
(\mathrm{tr}MW)^2\]/\[(n-p)(n-p+2)\] - E\[I\]^2\$\$ These assume normal
errors and do not condition on the observed kurtosis. On the simulation
above they restored the z-score to mean \\-0.09\\, \\sd = 1.03\\, and
the rejection rate to 4.3\\ nominal 5\\ precision.

**These moments are exact for \\e = My\\ and for nothing else**, so
`null = "auto"` does not guess from the fit's class: it rebuilds `X`,
regresses the response on it, and uses the residual moments only when
the supplied residuals *are* those OLS residuals to numerical tolerance.
A GWR wide enough to have collapsed to global OLS passes that test; the
same GWR at a working bandwidth does not.

**For the flexible backends neither null is exact**, and `"auto"` leaves
them on `"randomisation"` because forcing the OLS moments on them
measurably makes matters worse, not better. Measured on null data (\\n =
120\\, three smooth covariates, independent errors; nominal 5\\
one-sided):

|               |                   |              |
|---------------|-------------------|--------------|
| **backend**   | **randomisation** | **residual** |
| OLS           | 0.035             | 0.060        |
| random forest | 0.128             | 0.200        |
| GWR           | 0.000             | 0.000        |

The random forest is anticonservative under both. Its residuals here are
the **out-of-bag** ones
([`residuals.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.rf_fit.md)),
not in-sample fits. They are inflated instead of shrunk (measured sd
1.08 against a true 1.00), but they are not a linear projection of the
response and they are spatially heteroscedastic, so neither set of
moments describes their null distribution and the variance is
understated whichever is used (\\sd(z) \approx 1.3\\). GWR is
conservative under both, because it removes far more structure than a
rank-p projection does. Treat the p-value from those backends as a rough
indicator, and prefer spatially-blocked cross-validated residuals or an
explicit spatial covariance model when the answer has to carry weight.

A permutation null was considered and rejected: permuting the residual
vector destroys exactly the orthogonality that causes the bias, so its
mean is the exchangeable \\-1/(n-1)\\ by construction (measured:
\\-0.00840\\ against \\-1/(n-1) = -0.00840\\) and it reproduces the
randomisation null without correcting it.

## References

Cliff, A. D. and Ord, J. K. (1981) *Spatial Processes: Models and
Applications*. Pion, London. Section 8.3.

## See also

Other model evaluation:
[`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md),
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md),
[`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md)

## Examples

``` r
# Works on any spatial_fit; a forest keeps the example free of the optional
# GWR/Stan backends.
if (requireNamespace("ranger", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  n <- 120
  dat <- st_as_sf(
    data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
               elev = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
  # A strong east-west trend the predictor cannot explain: the residuals
  # should still carry spatial structure, and this is what detects it.
  dat$price <- 10 + 0.02 * (st_coordinates(dat)[, 1] - 5e5) +
    2 * dat$elev + rnorm(n)
  fit <- fit_rf_model(dat, "price", "elev", num_trees = 100, seed = 1)
  # I = 0.64, z = 15.5: strong positive residual autocorrelation, exactly
  # as constructed.  A z near 0 with a large p-value would be the opposite
  # verdict -- no structure the model failed to capture.
  residual_morans_i(fit)
}
#> Residual Moran's I = 0.6437   (E[I] = -0.0084, sd = 0.0419)
#>   z = 15.548, p = 1.653e-54
#>   null: randomisation moments, approximate for these residuals; n = 120, df = 119
#>   residual kurtosis 2.371 (3 = Gaussian)
#>   weights: 120 x 120 dgCMatrix, 8 neighbour(s) per row; not retained, keep_weights = TRUE to keep it
```
