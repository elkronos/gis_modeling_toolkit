# Summarize features by polygon/cell ID

Aggregates an sf point dataset into one row per cell. By default
computes counts and means, but the aggregation function is configurable.

## Usage

``` r
summarize_by_cell(
  assigned_points_sf,
  response_var = NULL,
  predictor_vars = NULL,
  id_col = "poly_id",
  agg_funs = list(mean = function(x) mean(x, na.rm = TRUE)),
  cells_sf = NULL,
  deff = 1,
  sac = NULL,
  deff_max_n = 500L,
  quiet = TRUE,
  conf_level = NULL,
  area = FALSE
)
```

## Arguments

- assigned_points_sf:

  An sf object with a cell identifier column.

- response_var:

  Optional response column name for per-cell aggregation.

- predictor_vars:

  Optional predictor column names for per-cell aggregation.

- id_col:

  Preferred name of the polygon/cell ID column.

- agg_funs:

  Named list of aggregation functions. Default
  `list(mean = \(x) mean(x, na.rm = TRUE))`. Additional common options:
  `median`, `sum`, `sd`.

- cells_sf:

  Optional polygon sf layer to join cell geometries onto the output.
  When supplied, the return value is an sf object with the polygon
  geometry from cells_sf, with one row per cell in `cells_sf`. Cells
  that no feature fell in are kept, with `NA` summaries. Duplicate ID
  values in `cells_sf` would multiply those rows, so they are reported
  with a warning. When NULL (default), a plain data.frame/tibble is
  returned (previous behaviour).

- deff:

  Design-effect adjustment for standard errors. One of:

  `1` (default)

  :   No adjustment; the classic IID standard error `sd / sqrt(n)`,
      which assumes the observations within a cell are independent. No
      `"deff_applied"` attribute is attached, so a result carrying none
      was computed this way.

  `"variogram"`

  :   Compute a per-cell design effect from a fitted variogram: for `n`
      points in a cell with correlation matrix `R`, the effective sample
      size of the mean is `n^2 / sum(R)`, so `deff = sum(R) / n`. This
      generalises Kish (substituting a constant off-diagonal correlation
      recovers `1 + (n - 1) * rho` exactly) but lets correlation decay
      with distance, which matters increasingly as cells get larger and
      Kish's single-`rho` assumption degrades. Supply the fit via `sac`,
      or it is estimated when `response_var` is given and 'gstat' is
      available. Exponential, spherical and Gaussian models are
      supported, with a nugget and with several structured components
      (each weighted by its partial sill); a model of any other family
      falls back to `deff = 1` with a warning naming it.

  `"kish"`

  :   Estimate per-variable-type intra-class correlations (ICCs) from
      the grouped data using a one-way random-effects ANOVA
      decomposition (one ICC for the response variable and a separate
      ICC for the predictor variables), then apply Kish's formula per
      cell: `deff_i = 1 + (n_i - 1) * rho`. The ICC is the ANOVA
      (method-of-moments) estimator with Donner's `n0` for unequal cell
      sizes, not the REML estimate a mixed model returns: on a single
      unbalanced draw the two can differ by 0.1–0.2 (one check with cell
      sizes 3 to 77 and a true ICC of 0.5 gave 0.33 against REML's
      0.51), while balanced designs agree to about 0.01. Neither is
      wrong, so do not read the difference as a defect. When multiple
      columns are pooled for a single ICC estimate (e.g. several
      predictor variables), each column is z-scored before pooling so
      that variables with different scales contribute equally to the
      variance decomposition. The response-specific ICC is used for
      response SEs and the predictor-specific ICC for predictor SEs. The
      response's ICC is never applied to predictor columns or vice
      versa, and it is the response's ICC (not the predictors') that
      sets `cell_weight` whenever a response was given. Requires at
      least 2 cells with 2+ observations and at least 2 residual degrees
      of freedom (`N - k >= 2`); the ICC is taken as 0 (no correction,
      no `"deff_applied"` attribute) otherwise, and likewise when the
      estimate itself comes out at or below 0.

  A positive number

  :   Applied as a uniform design effect to every cell, as
      `sd * sqrt(deff / n)`, exactly `sqrt(deff)` times the naive SE.
      Use when you have an external estimate of the design effect.
      Anything that is not a single number `>= 1` (including a value
      below 1, which would *shrink* the standard errors) is refused with
      a warning and replaced by 1.

- sac:

  Optional `sac_range` object from
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
  used when `deff = "variogram"`. Supplying one avoids re-fitting the
  variogram and lets you inspect the fit the design effect is based on.
  A `sac_range` whose fit was *rejected* (its `status` is not `"ok"`)
  carries no usable correlation function, so `deff` falls back to 1 with
  a warning and does not correct by a shape that was not trusted enough
  to report a range.

- deff_max_n:

  Cells with more than this many points are subsampled before forming
  the `n x n` correlation matrix used by `deff = "variogram"`. Default
  500.

- quiet:

  Logical; suppress this function's progress
  [`message()`](https://rdrr.io/r/base/message.html)s. It does not
  silence R warnings, nor the package's console log echo (see
  [`spatialkit_quiet`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)
  for that). Default `TRUE`, unlike the tessellation functions, whose
  default is `FALSE`.

- conf_level:

  Optional confidence level in (0, 1), such as `0.95`. When given, every
  numeric response and predictor column also gets `..neff_*`, `..df_*`,
  `..ci_lo_*` and `..ci_hi_*` (see "Confidence intervals"). Default
  `NULL`: no interval columns, and the frame is exactly what it was
  before this argument existed.

- area:

  Logical, default `FALSE`. With `TRUE`, and `cells_sf` supplied, the
  result gains `cell_area` (each cell's planar area in the squared units
  of `cells_sf`'s CRS) and `n_per_area` (the count of rows in the cell
  over that area, a point density; a rate of anything else is that
  thing's `agg_funs` sum over `cell_area`). The request is **refused**
  with an error, not answered with a number, when the cells' CRS
  distorts areas across them by more than 1 percent, measured as the
  spread of planar-to-geodesic area ratios over the cells: a density is
  a comparison between cells and is meaningless where the map scale
  differs from one cell to the next. Inside a UTM zone the spread is
  under 0.3 percent and the request goes through; the conterminous
  United States forced into one zone (14 percent), or a few degrees of
  latitude in Web Mercator (4 percent at 48N), does not.
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  with `purpose = "area"` chooses an equal-area CRS for lon/lat input;
  build the cells in it. The measured spread is attached as
  `attr(, "area_error")`.

## Value

A tibble/data.frame (or sf if cells_sf given) with per-cell summaries:
the ID column, `n` (rows in the cell), one column per `agg_funs` entry
per variable, `..sd_*` / `..se_*` for every numeric response and
predictor, `..neff_*` / `..df_*` / `..ci_lo_*` / `..ci_hi_*` for the
same columns when `conf_level` is given, `cell_weight`, and `cell_area`
/ `n_per_area` when `area = TRUE`. An input column also called `n` is
not allowed to shadow the count.

When a correction was actually applied, an attribute `"deff_applied"` is
attached recording it: `method` plus `icc_resp`/`icc_pred` for `"kish"`,
`deff`/`deff_rows`/`rbar`/`crs`/`max_n` for `"variogram"` (`deff` is the
design effect at the primary variable's non-missing count per cell,
`deff_rows` at the cell's row count, which is the vector the log line
summarises as a median and a max), and `deff` alone for a fixed number.
When `cells_sf` is supplied, *every* per-cell vector in that attribute
(`deff`, `deff_rows` and `rbar` alike) is realigned to the joined row
order, so `deff[i]` and `rbar[i]` still describe row `i`; cells with no
observations carry `NA`. No attribute is attached when no correction was
applied: `deff = 1`, a `deff = "kish"` ICC of 0, or a `"variogram"`
request that could not be fitted. A `deff = "kish"` request always
records the ICCs it estimated on an attribute `"icc"` (`resp` and
`pred`, `NA` for a variable type with no numeric column), whether or not
they were positive enough to apply, so a result with no `"deff_applied"`
still says what the ICC came out as.

The ID column keeps its input type when `cells_sf`'s ID column and the
summarised IDs already have the same class. When the classes differ,
both are coerced to character in order to join (logged as a warning),
and the returned ID column is therefore character.

## Details

This is the third step of the package's pipeline, taking the labelled
layer from
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
down to cell level. What distinguishes it from a plain
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html) +
`summarise()` is that it carries the *uncertainty* of each aggregate
with it: alongside every mean it returns a within-cell standard
deviation, a standard error and an observation count, and it can correct
that standard error for within-cell spatial autocorrelation via `deff`.
Reach for it whenever the cell-level values will be modelled or mapped,
because a cell mean over 2 observations and one over 200 are not the
same measurement and nothing downstream can tell them apart otherwise.

In addition to user-specified aggregation functions, this function
always computes within-cell standard deviation (`..sd_<var>`) and
standard error (`..se_<var>`) for every numeric response/predictor
column, plus an `n` column (rows falling in the cell) and a
`cell_weight` column. These columns let downstream models account for
the fact that a cell with 2 observations carries more aggregation
uncertainty than one with 200.

`cell_weight` is the *effective* sample size of the primary variable:
the response when one was supplied, otherwise the first predictor. It
counts that variable's non-missing rows, not all rows (a cell of 10 rows
with 3 finite responses carries 3 observations' worth of information
about the response, not 10), and it is divided by that cell's design
effect when `deff` applied one. With `deff = 1` and no missing values it
equals `n`. Pass it as the `weights` argument of a downstream
regression.

## Spatial autocorrelation and standard-error bias

By default (`deff = 1`), the `..se_*` columns are computed as
`sd / sqrt(n)`, which assumes observations within each cell are
independent. When data are spatially autocorrelated (the common case for
the spatial workflows this package supports), within-cell observations
are typically positively correlated, so the effective sample size is
smaller than `n`. The naive SE is therefore **anticonservative** (too
small), and downstream weighted regressions using `cell_weight` or
`..se_*` columns will produce overconfident standard errors for cells
with strong intra-cell correlation.

Setting `deff = "kish"` applies an approximate correction using Kish's
design effect. Separate intra-class correlations (ICCs) are estimated
for response and predictor variables via a one-way random-effects
decomposition across all cells. Each variable type's ICC is used for its
own SE adjustment, and each cell's effective sample size is reduced to
`n_i / (1 + (n_i - 1) * rho)`. This is a first-order correction that
does not require a full spatial covariance model but does require enough
cells and observations for a stable ICC estimate. You may also pass a
fixed numeric design effect (e.g. `deff = 2`) to uniformly inflate
standard errors: an externally supplied constant is applied as
`sd * sqrt(deff / n)`, exactly `sqrt(deff)` times the naive SE in every
cell. (The `E[s^2]` correction that the estimated design effects also
apply is derived from within-cell correlation and would not be justified
for a number the caller chose.)

Even with the Kish correction, the adjusted SE is an approximation. For
rigorous inference under spatial dependence, consider fitting an
explicit spatial covariance model (e.g. via
[`fit_bayesian_spatial_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)).

## What the standard error estimates

The `..se_*` columns are the standard error of the cell mean **as an
estimate of the population (grand) mean**: the unconditional quantity,
in which the cell's own realised deviation is part of the error. That is
the right quantity when cells are treated as samples from a common
population, and the design-effect correction is calibrated for it:
measured 95% interval coverage of the grand mean is 0.95 with
`deff = "kish"` (and 0.29 with the naive SE) on exchangeable within-cell
correlation, and 0.93 with `deff = "variogram"` on a simulated Gaussian
field.

It is **not** the standard error of the cell's own mean (the block
average over that cell), which is what a cell-level map or a regression
on cell values usually wants. For that quantity the naive `sd / sqrt(n)`
is the better of the two on offer: measured coverage 0.95 under
exchangeable within-cell correlation, against very nearly 1.00 for the
design-effect-corrected SE, which is about five times too wide. That
0.95 is exact under the exchangeable model and holds under a spatial
covariance model only when the cell's points are spread through the
cell; with *clustered* sampling inside a cell it is anticonservative for
the block average too (measured 0.58), and the right quantity there is a
block-kriging variance, which this function does not compute. Use `deff`
when the cell means feed a population-level inference; leave it at 1
when they are measurements of the cells themselves and the sampling
within cells is reasonably uniform.

## Design effects and variable types

`deff = "kish"` estimates a separate ICC for response and predictor
variables and applies each to its own columns. `deff = "variogram"` fits
or accepts **one** correlation function and applies it to every numeric
column, because a variogram is a property of the field being modelled
rather than of a variable type; a predictor whose spatial structure
differs markedly from the response's will have its SE corrected by the
response's correlation. The internally estimated variogram is fitted to
the **response itself**, never to OLS residuals, whatever
`predictor_vars` holds: the `..se_resp_*` columns estimate the grand
mean of the response, so the correlation to correct for is the
response's own. (A residual variogram, whose correlation is that of the
part the predictors do not explain, is weaker; using it here dropped
grand-mean coverage from 0.93 to 0.51 the moment a predictor was
listed.) Pass `sac` explicitly when you want a different variogram, such
as a residual one from `estimate_sac_range(..., predictor_vars = )`, and
check `attr(sac, "detrended")` to know which you have.

## Confidence intervals

With `conf_level` set, every numeric response and predictor column gains
four more columns: `..neff_*`, that column's effective sample size in
the cell (its non-missing count over its design effect, the per-column
version of `cell_weight`); `..df_*`, the degrees of freedom the interval
uses; and `..ci_lo_*` / `..ci_hi_*`, a t interval for the **cell mean as
an estimate of the grand mean**. The estimand is the same as `..se_*`'s,
so everything in "What the standard error estimates" applies to it,
including that it is not an interval for the cell's own block average.
The interval is `mean +/- qt((1 + conf_level) / 2, df) * se`, centred on
the plain mean of the column's non-missing values whatever `agg_funs`
computes, and is `NA` wherever the standard error is (a single
observation; complete redundancy under `deff`).

The degrees of freedom are **not** `neff - 1`. The interval's spread
comes from the within-cell sample variance, and under exchangeable
correlation (`deff = "kish"`) that variance keeps its `n - 1` degrees of
freedom whatever the design effect. The design effect inflates the
mean's variance and biases `s^2`, both of which the standard error
already corrects, and the resulting pivot is exactly t on `n - 1` df.
Measured 95% coverage of the grand mean on the Kish path, 20 cells of 20
at an ICC of 0.2 / 0.6 / 0.9: 0.954 / 0.953 / 0.952 with `n - 1`,
against 0.992 / 1.000 / 1.000 with `neff - 1`, which is not an interval
so much as a refusal to say anything. The effective-sample-size degrees
of freedom of Faes et al. (2009) belong to a mean estimated across many
correlated units whose variance is estimated from the spread between
them; they do not transfer to a single cell's mean with a variance
estimated from inside it. `..df_*` is therefore `n - 1` at `deff = 1`,
for a numeric `deff` and for `"kish"`. Under `deff = "variogram"`
correlation decays with distance and the within-cell variance loses
degrees of freedom to it; `..df_*` is then the Satterthwaite (1946)
moment-matched df of `s^2` under the fitted correlation matrix, a
fraction of `n - 1` that shrinks as the range grows. Measured on a
Gaussian field with an exponential range of 150 and 400 on a 1000-unit
domain, 16 cells of about 25 points: 0.960 and 0.958 with that df,
against 0.931 and 0.918 with `n - 1`, and 0.34 and 0.19 for the naive
`deff = 1` interval. The interval is only as good as the design effect
under it: a mis-specified variogram, or an ICC estimated from too few
cells, moves the coverage with it.

## References

Faes, C., Molenberghs, G., Aerts, M., Verbeke, G. and Kenward, M. G.
(2009). The effective sample size and an alternative small-sample
degrees-of-freedom method. *The American Statistician*, 63(4), 389–399.
[doi:10.1198/tast.2009.08196](https://doi.org/10.1198/tast.2009.08196)

Satterthwaite, F. E. (1946). An approximate distribution of estimates of
variance components. *Biometrics Bulletin*, 2(6), 110–114.
[doi:10.2307/3002019](https://doi.org/10.2307/3002019)

## See also

[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md),
which produces the input layer;
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
for the cells themselves.

Other aggregation:
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md),
[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md),
[`kriging_adequacy()`](https://elkronos.github.io/gis_modeling_toolkit/reference/kriging_adequacy.md),
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)

## Examples

``` r
library(sf)
set.seed(1)
n <- 200
east  <- runif(n, 0, 100)
north <- runif(n, 0, 100)
# A response with spatial structure, so the within-cell ICC is not zero.
pts <- st_as_sf(
  data.frame(x = 5e5 + east, y = 5e6 + north,
             val = 0.05 * east + 0.05 * north + rnorm(n, sd = 0.5)),
  coords = c("x", "y"), crs = 32632
)
bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
  c(5e5, 5e6), c(5e5 + 100, 5e6), c(5e5 + 100, 5e6 + 100),
  c(5e5, 5e6 + 100), c(5e5, 5e6)
))), crs = 32632))
grid <- create_grid_polygons(bnd, target_cells = 9, type = "square")
assigned <- assign_features_to_polygons(pts, grid)

# IID standard errors (default) vs Kish design-effect adjustment
naive <- summarize_by_cell(assigned, response_var = "val")
kish  <- summarize_by_cell(assigned, response_var = "val", deff = "kish")
data.frame(n = naive$n,
           se_naive = naive$..se_resp_val,
           se_kish  = kish$..se_resp_val)
#>    n  se_naive  se_kish
#> 1 25 0.1282574 1.576773
#> 2 28 0.1687903 2.195277
#> 3 25 0.1586372 1.950257
#> 4 21 0.1797135 2.026192
#> 5 23 0.2010768 2.371742
#> 6 21 0.1974045 2.225651
#> 7 10 0.2315705 1.809441
#> 8 27 0.1499232 1.914966
#> 9 20 0.2162204 2.379509
attr(kish, "deff_applied")   # method, icc_resp, icc_pred, per-cell deff
#> $method
#> [1] "kish"
#> 
#> $icc_resp
#> [1] 0.8572554
#> 
#> $icc_pred
#> [1] NA
#> 
#> $deff
#> [1] 21.574129 24.145895 21.574129 18.145108 19.859618 18.145108  8.715298
#> [8] 23.288640 17.287852
#> 

# A 95% interval for each cell mean as an estimate of the grand mean, on
# the Kish-corrected standard error and n - 1 degrees of freedom.
ci <- summarize_by_cell(assigned, response_var = "val", deff = "kish",
                        conf_level = 0.95)
ci[, c("poly_id", "n", "resp_mean_val", "..neff_resp_val", "..df_resp_val",
       "..ci_lo_resp_val", "..ci_hi_resp_val")]
#> # A tibble: 9 × 7
#>   poly_id     n resp_mean_val ..neff_resp_val ..df_resp_val ..ci_lo_resp_val
#>     <int> <int>         <dbl>           <dbl>         <dbl>            <dbl>
#> 1       1    25          1.74            1.16            24          -1.52  
#> 2       2    28          3.25            1.16            27          -1.26  
#> 3       3    25          5.02            1.16            24           1.000 
#> 4       4    21          3.37            1.16            20          -0.857 
#> 5       5    23          5.01            1.16            22           0.0947
#> 6       6    21          7.00            1.16            20           2.35  
#> 7       7    10          5.43            1.15             9           1.34  
#> 8       8    27          6.84            1.16            26           2.90  
#> 9       9    20          7.94            1.16            19           2.96  
#> # ℹ 1 more variable: ..ci_hi_resp_val <dbl>
```
