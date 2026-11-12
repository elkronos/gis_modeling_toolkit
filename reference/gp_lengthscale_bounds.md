# Compute length-scale bounds for a squared-exponential GP

Computes lower and upper bounds for the length-scale \\\ell\\ of a
squared-exponential (exponentiated-quadratic) kernel, \\k(h) = \exp(-h^2
/ (2\ell^2))\\, from the spread of the pairwise distances between the
points.
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
calibrates its length-scale prior over them. The effective range, where
the correlation drops to 5 percent, is \\\ell \sqrt{2 \ln 20} \approx
2.45\\\ell\\.

## Usage

``` r
gp_lengthscale_bounds(coords_xy, q_small = 0.25, max_n = 1000L)
```

## Arguments

- coords_xy:

  Numeric matrix or data.frame of coordinates with at least two columns;
  the first two are used, and replicated rows are collapsed before the
  distance quantiles are taken, because
  [`brms::gp()`](https://paulbuerkner.com/brms/reference/gp.html)
  (`gr = TRUE`, its default) reduces its covariates to unique rows: a
  station sampled many times would otherwise weight the quantile by how
  often it was measured rather than by where the sites are. A row with a
  coordinate that is missing, not finite or beyond 1e150 has no
  location: it is left out, with a warning of class
  `spatialkit_rows_dropped` that gives the count.

- q_small:

  Quantile of the pairwise distances for the lower bound, a single
  number in `[0, 1]`. Default 0.25. On clustered data a low quantile
  reflects the spacing within clusters.

- max_n:

  Maximum number of *distinct* points to use in the distance
  computation. Default 1000. Set to `Inf` to use all of them. Above it,
  a subsample is drawn with a fixed seed from the distinct points in
  coordinate order, so the bounds do not depend on the row order of
  `coords_xy`, and the caller's random-number stream is left where it
  was.

## Value

Named numeric vector `c(lower, upper)`: length-scales in the units of
`coords_xy`. `lower` is the `q_small` quantile of the distinct pairwise
distances divided by \\\sqrt{2 \ln 20}\\ (about 2.45); `upper` is the
largest pairwise distance divided by the same factor, and at least 1.2
times `lower`. `c(lower = 0.001, upper = 1)` when fewer than two
distinct located rows or no positive distances remain.
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
standardises the coordinates before it calls this function, so its
bounds are in those standardised units.

## Details

The lower bound is the `q_small` quantile of the pairwise distances (by
default the 25th percentile) taken as an effective range, and the upper
bound is the largest pairwise distance taken as an effective range. The
inverse-gamma prior that
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
calibrates on them puts 1 percent of its mass below the lower bound and
1 percent above the upper. Both are fixed fractions of the spread of
pairwise distances, so they do not shrink as points are added to the
same area: more points help the data identify a short range, but they
make neither these bounds nor the basis derived from them finer.

These are the bounds a length-scale prior is calibrated over, not the
scales a fitted model can resolve, which depend on the basis size
(`gp_k` in
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md),
which reports it as `$info$gp_ell_min`, and in CRS units as
`$info$gp_resolution`). A surface whose range is shorter than the lower
bound needs a basis fine enough to represent it and a prior that allows
it. A larger `gp_k` under the prior calibrated here changes almost
nothing, because the prior still holds the length-scale above this
bound.
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
therefore moves the prior's lower tail down to the length-scale the
basis resolves whenever `gp_resolution` or an explicit `gp_k` asks for a
finer basis than the one derived from these bounds.

Large layers are subsampled (`max_n`) to bound the cost of the pairwise
distances, which grows with the square of the points.

## See also

Other model fitting:
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md),
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md),
[`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md),
[`new_spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md)

## Examples

``` r
# 500 points spread evenly over a square of side 100 km, in metres.
set.seed(1)
xy <- cbind(runif(500, 0, 1e5), runif(500, 0, 1e5))
b <- gp_lengthscale_bounds(xy)
b
#>    lower    upper 
#> 13445.22 55051.07 
# The bounds as effective ranges: the 25th percentile of the pairwise
# distances, about a third of the side, and the largest distance.
b * sqrt(2 * log(20))
#>     lower     upper 
#>  32910.48 134751.08 
```
