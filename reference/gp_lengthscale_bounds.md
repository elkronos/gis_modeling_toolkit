# Heuristic length-scale bounds for a squared-exponential GP

Computes sensible prior bounds for the GP length-scale parameter
\\\ell\\ of a squared-exponential (exponentiated-quadratic) kernel,
\\k(h) = \exp(-h^2 / (2\ell^2))\\. The "effective range" where
correlation drops to ~5% is \\\ell \sqrt{2 \ln 20} \approx 2.45\\\ell\\.

## Usage

``` r
gp_lengthscale_bounds(coords_xy, q_small = 0.25, max_n = 1000L)
```

## Arguments

- coords_xy:

  Numeric matrix or data.frame of coordinates with at least two columns;
  the first two are used, and replicated rows are collapsed before the
  distance quantiles are taken.
  [`brms::gp()`](https://paulbuerkner.com/brms/reference/gp.html)
  defaults to `gr = TRUE` and reduces its covariates to unique rows, so
  a heavily-sampled station would otherwise weight the quantile by how
  often it was measured rather than by where the sites are.

- q_small:

  Numeric quantile for the lower bound, a single number in `[0, 1]`.
  Default 0.25. Previous versions used 0.1, but the 10th percentile can
  be dominated by within-cluster spacing in clustered data, producing a
  misleadingly small lower bound.

- max_n:

  Maximum number of *distinct* points to use in the distance
  computation. Default 1000. Set to `Inf` to use all of them.

## Value

Named numeric vector `c(lower, upper)` on the length-scale;
`c(lower = 0.001, upper = 1)` when fewer than two distinct locations or
no positive distances remain.

## Details

Subsamples large datasets to avoid O(n^2) memory and time cost.

These are the bounds a length-scale *prior* is calibrated over, not the
scales a fitted model can resolve: that depends on the basis size
(`gp_k` in
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md),
which reports it as `$info$gp_ell_min`). Both bounds are fixed fractions
of the spread of pairwise distances, so they do not shrink as points are
added to the same area. A surface whose range sits below what the basis
resolves needs a larger `gp_k`: more points help the data identify a
short range, but they make neither these bounds nor the derived basis
finer.

## See also

Other model fitting:
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md),
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md),
[`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md),
[`new_spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md),
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)

## Examples

``` r
set.seed(1)
xy <- cbind(runif(50), runif(50))
gp_lengthscale_bounds(xy)
#>     lower     upper 
#> 0.1253801 0.4778780 
```
