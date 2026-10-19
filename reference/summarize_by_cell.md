# Summarize features by polygon/cell ID

Aggregates an sf point dataset into one row per cell. By default
computes counts and means, but the aggregation function is configurable.
It also returns a pooled value for each cell, its mean pulled toward the
mean of the cells around it (`estimate = "pooled"`, the default;
`estimate = "mean"` leaves it out), and can return a range for one
further observation.

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
  area = FALSE,
  estimate = c("pooled", "mean"),
  obs_level = NULL
)
```

## Arguments

- assigned_points_sf:

  An sf object with a cell identifier column.

- response_var:

  Optional response column name for per-cell aggregation: a single
  character string. A name that is not a column, or a column that is not
  numeric, is skipped with a warning.

- predictor_vars:

  Optional predictor column names for per-cell aggregation. Names that
  are not columns, and columns that are not numeric, are skipped with a
  warning.

- id_col:

  Preferred name of the polygon/cell ID column.

- agg_funs:

  Named list of aggregation functions. Default
  `list(mean = \(x) mean(x, na.rm = TRUE))`. Additional common options:
  `median`, `sum`, `sd`. A single function (`agg_funs = median`) or a
  character vector of function names (`c("median", "sum")`) is also
  accepted. A single function is named after the expression passed: a
  name gives that name (`median`, or `f` for a variable `f` holding a
  function), [`stats::median`](https://rdrr.io/r/stats/median.html)
  gives `median`, and any other expression gives `agg1`; pass a named
  list to choose the name. Anything else falls back to the default mean
  with a warning. The functions receive each column with its infinite
  values already set to `NA` (see "Details"), so they need
  `na.rm = TRUE` as they do for missing values.

- cells_sf:

  Optional polygon sf layer to join cell geometries onto the output.
  When supplied, the return value is an sf object with the polygon
  geometry from cells_sf, with one row per cell in `cells_sf`. Cells
  that no feature fell in are kept, with `NA` summaries. Duplicate ID
  values in `cells_sf` would multiply those rows, so they are reported
  with a warning. Its ID column is the first of `id_col`, `"poly_id"`,
  `"polygon_id"`, `"id"`, `"cell_id"` and `"grid_id"` it carries, the
  list and order
  [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
  reads the polygons' IDs from. A summarised ID that matches no cell is
  reported with a warning, since the join drops it with its points; a
  `cells_sf` with none of those columns, or one that is not an sf
  object, gives a warning and a plain data frame (an error with
  `area = TRUE`). When NULL (default), a plain data.frame/tibble is
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
      (each weighted by its partial sill) and with gstat's 2-D geometric
      anisotropy (`vgm(..., anis = c(angle, ratio))`), applied as gstat
      applies it. A request that cannot be met falls back to `"kish"`,
      with a warning of class `"spatialkit_deff_fallback"` that gives
      the reason and the ICCs: a model of any other family (named in the
      warning), or no usable model (none supplied and none could be
      estimated, or a rejected fit that could not be replaced by an
      estimate). The standard errors, intervals and `cell_weight` are
      then Kish's, and the uncorrected ones only where Kish's own ICC is
      0; the `deff_method` column says which on every row (see "Value").
      The variogram route needs an identified range, which a layer does
      not always have, so read `deff_method` for the fallback. A model
      with no structured component (a pure nugget) implies that distinct
      observations are uncorrelated, so it is applied as a design effect
      of 1 in every cell, not treated as a fallback. That includes the
      model
      [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
      returns when it finds no spatial structure (an `NA` range whose
      `rejected_reason` says so): the answer is then a design effect of
      1, `deff_method` reads `"variogram"`, and no fallback warning is
      raised. Points with empty geometry count towards their cells'
      values but not towards the correlation, with a warning.

  `"kish"`

  :   Estimate an intra-class correlation (ICC) for each summarised
      column from the grouped data using a one-way random-effects ANOVA
      decomposition (one ICC for the response and one for each
      predictor), then apply Kish's formula per cell:
      `deff_i = 1 + (n_i - 1) * rho`. The ICC is the ANOVA
      (method-of-moments) estimator with Donner's `n0` for unequal cell
      sizes, not the REML estimate a mixed model returns: on a single
      unbalanced draw the two can differ by 0.1–0.2 (one check with cell
      sizes 3 to 77 and a true ICC of 0.5 gave 0.33 against REML's
      0.51), while balanced designs agree to about 0.01. Neither is
      wrong, so do not read the difference as a defect. Each column's
      ICC is used for that column's SEs and for no other: the response's
      ICC is never applied to a predictor, nor one predictor's to
      another. It is the response's ICC that sets `cell_weight` whenever
      a response was given, and the first predictor's own ICC otherwise.
      Requires at least 2 cells with 2+ observations and at least 2
      residual degrees of freedom (`N - k >= 2`); the ICC is taken as 0
      (no correction for that column) otherwise, and likewise when the
      estimate itself comes out at or below 0. The `"deff_applied"`
      attribute is attached when any ICC is positive, so when all are 0
      there is none.

  A positive number

  :   Applied as a uniform design effect to every cell, as
      `sd * sqrt(deff / n)`, exactly `sqrt(deff)` times the naive SE.
      Use when you have an external estimate of the design effect.
      Anything that is not a single finite number `>= 1` (including a
      value below 1, which would *shrink* the standard errors, and `NA`
      or `Inf`) is refused with a warning and replaced by 1.

- sac:

  Optional `sac_range` object from
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
  used when `deff = "variogram"`. Supplying one avoids re-fitting the
  variogram and lets you inspect the fit the design effect is based on.
  A `sac_range` whose fit was *rejected* (its value is `NA` and it
  carries a `rejected_reason` attribute) carries no usable correlation
  function, so it is set aside rather than correcting by a shape that
  was not trusted enough to report a range: the variogram is then
  estimated as if no `sac` had been given, with a plain warning saying
  so, or, where that is not possible, `deff` falls back to `"kish"` with
  the fallback warning, which names the rejection. The exception is a
  `sac` refused because no spatial structure was found, whose model is a
  pure nugget: it is used, and gives a design effect of 1. A `sac` with
  no `variogram_model` attribute – a plain number or a `units` object,
  say – is a range without a correlation function, and is set aside the
  same way. A `sac` fitted to residuals (`attr(sac, "detrended")`
  `TRUE`) is used as given, with a warning when it corrects response
  columns (see "Design effects and variable types").

- deff_max_n:

  Cells with more than this many points are subsampled before forming
  the `n x n` correlation matrix used by `deff = "variogram"`.
  Default 500. It must be a single number of at least 2 when
  `deff = "variogram"`; anything else is an error.

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
  `NULL`: no such columns.

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
  build the cells in it. Lon/lat cells are measured geodesically
  instead: their `cell_area` is
  [`sf::st_area()`](https://r-spatial.github.io/sf/reference/geos_measures.html)'s
  area in square metres (on the sphere with s2, sf's default), not a
  planar area in squared degrees, and the distortion check passes by
  construction; with s2 switched off, sf needs the lwgeom package for
  that area and the request is refused without it. The measured spread
  is attached as `attr(, "area_error")`.

- estimate:

  Which value of a cell the frame reports: `"pooled"` (default) or
  `"mean"`. With `"mean"` the frame is the one described above, the
  frame of earlier versions, whose default this was. With `"pooled"`
  every numeric response and predictor column also gets a pooled value,
  `w * m + (1 - w) * t`: the cell's own mean `m` and a target `t`, the
  mean of the rows of the 8 cells whose observations' centroids are
  nearest (all the other cells when there are fewer than 10; the
  polygons are not read). The own weight is
  `w = tau2 / (tau2 + sigma2 / n)`, where `n` is the cell's count,
  `sigma2` the spread of the observations about their cell's mean and
  `tau2` the mean squared distance between a cell's level and its
  target. Each column is pooled on its own rows: those with a value, a
  cell and a location. `agg_funs`, `deff`, `sac` and `conf_level` do not
  enter the pooled columns, and every other column is computed as with
  `"mean"`. Lon/lat coordinates, and coordinates without a CRS that look
  like degrees, are projected for the distances between centroids, as
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  projects them. The projection is chosen from every row of the layer,
  those without a value or a cell included, so on such input the rows a
  column does not use can still decide the CRS its distances are
  measured in; a layer projected beforehand is used in its own CRS. When
  the input carries the record of its cells' count (see "Value") and
  that record names a projected CRS, the one the profile worked in, the
  distances are measured in it, whatever CRS the layer is in. The value
  is matched as [`match.arg()`](https://rdrr.io/r/base/match.arg.html)
  matches a choice, so an abbreviation (`"m"`) is accepted and `NULL` is
  the default. See "Pooled values" for what the pooled value estimates
  and what its numbers rest on, "Value" for the columns and the
  attributes, and "Conditions" for the warnings.

- obs_level:

  Optional level in (0, 1), such as `0.90`, of an interval for one
  further observation in a cell, located like the data. It needs
  `estimate = "pooled"`: the interval is the pooled value plus the
  `(1 - obs_level) / 2` and `(1 + obs_level) / 2` quantiles of the
  leave-one-out residuals of the pooled values, so its two offsets are
  the same in every cell. It is a statement about observations on
  average, not about each cell. It is a separate argument from
  `conf_level`, whose interval is for a mean and is not a range for an
  observation. Any level strictly between 0 and 1 is computed the same
  way; only 0.90 has been checked (see "Intervals: which one answers
  which question"). Default `NULL`: no such interval.

## Value

A tibble/data.frame (or sf if cells_sf given) with per-cell summaries:
the ID column, `n` (rows in the cell), one column per `agg_funs` entry
per variable, `..sd_*` / `..se_*` for every numeric response and
predictor, `..neff_*` / `..df_*` / `..ci_lo_*` / `..ci_hi_*` for the
same columns when `conf_level` is given, `cell_weight`, `deff_applied`
and `deff_method` when a design effect was requested (any `deff` other
than 1), and `cell_area` / `n_per_area` when `area = TRUE`. An input
column also called `n` is not allowed to shadow the count.

`deff_applied` is `TRUE` on every row when a design-effect correction
was applied (exactly when the `"deff_applied"` attribute below is
attached; under `"kish"`, when any standard-error column was corrected)
and `FALSE` when the standard errors are the uncorrected ones: a refused
`deff`, or `"kish"` ICCs of 0 for every column (and `NA` on a `cells_sf`
row no point fell in). `deff_method` names the correction the standard
errors carry, since their column names do not change with it:
`"variogram"`, `"kish"` or `"fixed"` when the one requested was applied,
`"none"` when none was, and `"kish (variogram fallback)"` or
`"none (variogram fallback)"` when a `"variogram"` request could not be
met and Kish stood in for it (the second when Kish's ICCs were 0 too).
So `deff_applied` is `TRUE` on a Kish fallback row: test
`deff_method == "variogram"` to keep only rows the variogram corrected.
Unlike the attribute, both columns survive
[`rbind()`](https://rdrr.io/r/base/cbind.html) and
[`dplyr::bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html)
of many results. A fallback is also signalled by a warning of class
`"spatialkit_deff_fallback"` (a Kish ICC of 0 under a `"kish"` request
is reported on `attr(, "icc")` instead), which
`tryCatch(spatialkit_deff_fallback = )` catches without matching the
message. It is raised when the correction requested was not the one
applied: a refused number (the standard errors are then the uncorrected
ones) or a `"variogram"` request that fell back to Kish; the condition
of the second carries `requested` and `used`.

When a correction was actually applied, an attribute `"deff_applied"` is
attached recording it: `method` plus `icc_resp`/`icc_pred` and `deff`
for `"kish"` (`icc_pred` is a vector named by predictor column, one ICC
each; `deff` is the primary variable's per-cell design effect, all 1
when only a predictor's ICC was positive; `requested` and
`fallback_reason` are added when a `"variogram"` request fell back to
it), `deff`/`deff_rows`/`rbar`/`crs`/`max_n` for `"variogram"` (`deff`
is the design effect at the primary variable's non-missing count per
cell, `deff_rows` at the cell's row count, which is the vector the log
line summarises as a median and a max), and `deff` alone for a fixed
number. When `cells_sf` is supplied, *every* per-cell vector in that
attribute (`deff`, `deff_rows` and `rbar` alike) is realigned to the
joined row order, so `deff[i]` and `rbar[i]` still describe row `i`;
cells with no observations carry `NA`. No attribute is attached when no
correction was applied: `deff = 1`, or a `deff = "kish"` request (or
Kish fallback) whose ICCs are all 0. A `deff = "kish"` request, and a
`"variogram"` request that fell back to it, always records the ICCs it
estimated on an attribute `"icc"` (`resp`, and `pred`, a vector named by
predictor column; `NA` for a variable type with no numeric column),
whether or not they were positive enough to apply, so a result with no
`"deff_applied"` still says what the ICCs came out as.

The ID column keeps its input type when `cells_sf`'s ID column and the
summarised IDs already have the same class. When the classes differ,
both are coerced to character in order to join (logged as a warning),
and the returned ID column is therefore character. Whole numbers are
written out in full for that (`"100000"`, never `"1e+05"`), so an
integer and a double ID of the same cell still match.

With `estimate = "pooled"`, each numeric column `v` with its prefix
(`resp_` or `pred_`) gains `..pooled_<prefix><v>` (the pooled value),
`..pooled_w_<prefix><v>` (the own weight, 0 to 1; it is not a regression
weight) and `..pooled_target_<prefix><v>` (the target), and with
`obs_level` also `..obs_lo_<prefix><v>` and `..obs_hi_<prefix><v>`. They
follow every column listed above and precede `cell_area` and
`n_per_area`, the response's first, then each predictor's. A cell with
no row that has a value, a cell and a location gets `NA`, never a value
made from its neighbours. Two attributes are attached.
`attr(, "pooled")`, of class `"spatialkit_pooled"`, records what was
fitted: `toward`, `neighbours`, `crs` (the CRS the centroid distances
were measured in), `obs_level`, `conf_level`, `intervals` and `columns`,
a list with one entry per pooled column holding `column`, `role`
(`"response"` or `"predictor"`), `n` (rows used), `cells`, `df` (values
beyond one per cell, which `sigma2` rests on), `neighbours` (the cells
in each target), `sigma2`, `tau2`, `own_weight` (the median weight),
`share_alone` (the share of rows that are the only one in their cell),
`loo` and `loo_se` (the leave-one-out mean squared error of the pooled
value, and its standard error over rows), `loo_one_cell` (the same error
with one mean for all rows), `structured` and `structure_p` (whether the
cells predict a left-out observation better than one mean, and the
p-value of the F test of the cell means against one mean, which has to
be below 0.05 as well), `offsets` (the two offsets of the interval),
`repeated`, `skipped`, `constant` and `status`. `attr(, "intervals")` is
the same data frame as `attr(, "pooled")$intervals`: one row per family
of interval columns in the frame, saying what it covers, what it is
centred on, its level and its method. Away from the observations no
interval here is calibrated.

Where the pooled columns of a column cannot be computed they are `NA`,
with a warning of class `"spatialkit_pool_unavailable"` and not an
error, and the plain columns are returned. "Conditions" lists that class
and the other warnings of the pooled columns, with what each condition
carries.

When `assigned_points_sf` carries `attr(, "count_rule")`, or `cells_sf`
does (the assigned layer's is read first), the result carries that
record with one more element, `check`, whatever `estimate` is. The
record is the one
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
wrote: `rule`, `estimate` (the cell value the count was chosen for),
`cells`, `obs_per_cell`, `response`, `scored` (what the profile scored
at the count: the rows, the cells that hold one, the leave-one-out error
among them, and the centres the rows were given to), `partition`, `crs`
and `scored_cells`. It has class `"spatialkit_count_rule"` and a print
method. `check$estimate` is `"match"` or `"mismatch"`: whether the
summary reports the estimate the count was chosen for (`NA` for a count
chosen for neither). A mismatch is a warning of class
`"spatialkit_estimate_mismatch"` (fields `chosen_for`, `reported`,
`rule`, `cells`), and the summary is returned as asked. `check$cells` is
`"match"`, `"differ"` or `"not checked: <reason>"`: the response is
fitted again on the rows and cells given, with the estimate the profile
scored, and its scored rows, its occupied cells and its leave-one-out
error are compared with the record's, the error to a relative tolerance
of 1e-8. The first that differs is a warning of class
`"spatialkit_cells_differ"` (fields `what`: `"n"`, `"cells"` or `"loo"`;
`scored`; `found`): rows were dropped or moved between cells, or the
layer is another one. The same class with `what = "record"` (fields
`field`, `scored`, `found`) is raised when the assigned layer and
`cells_sf` carry different records, so the rows were assigned to other
cells than the `cells_sf` given: the cells of another selection, or the
scored cells of a selection beside other cells built from its count,
such as a lattice sized by it. `check$cells` is then `"differ"`,
whatever the comparison of the rows finds. The comparison is not made,
and `check$cells` says why, when the response is another column or was
not summarised, the profile read the selection half
(`select_on = "split"`), the cells are not the ones the profile scored
(`scored_cells` is `FALSE` or `NA`), the profile scored no estimate, the
input has no coordinates, the features are not points (see
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)),
or the fit was stopped by the limit on an exact search without the 'FNN'
package. Cells and rows without a record, such as polygons the caller
brings, are not checked and raise neither warning.

A row that is the same distance from two centres lies on the boundary
between their cells, which happens on gridded locations.
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
scores such a row in the cell of the lower-numbered centre, and
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md),
on cells built from the centres the profile scored, gives a point on
that boundary the same cell, measured in the CRS the record names, so
the comparison reads `"match"` on gridded locations too.
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
draws those cells in the CRS the record names, whatever CRS the seeds
arrive in, unless its `crs` names another projected CRS. Cells built in
another projected CRS are the nearest-centre cells of that CRS, so a row
near a boundary, on it or not, can sit in another cell than the profile
put it in. So can a row joined in longitude and latitude (rows and cells
both in degrees), and a row that reached the record's CRS by another
transformation than the profile's, which can settle a tie the other way.
The comparison then reports `"differ"`, with `what = "loo"` when it
finds the same rows in the same number of cells with another
leave-one-out error. The summary's own numbers describe the cells as
assigned.

The comparison uses the fit of the pooled columns when that is a fit of
the estimate the profile scored, and fits the response once more
otherwise. Without the 'FNN' package that fit compares every pair of
cell centroids (pooled values) or, for each row alone in its cell, every
other location (plain means), and above
`options(spatialkit.exact_pairs_max = )` pairs it is not made
(`"not checked: work limit"`). A summary is not compared when neither
layer carries the record: `attr(x, "count_rule") <- NULL` on the
assigned layer and on `cells_sf` removes it.

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

An infinite value (`Inf`, `-Inf`) in a response or predictor column is
treated as missing, as `NA` and `NaN` are, in every statistic of that
column: the `agg_funs` results, `..sd_*`, `..se_*`, the interval
columns, `..neff_*`, `cell_weight` and the design effect. `n` still
counts the row. A warning of class `"spatialkit_nonfinite"` names the
columns and gives the number of rows and cells affected (the condition
carries `n` and `columns`).
[`kriging_adequacy()`](https://elkronos.github.io/gis_modeling_toolkit/reference/kriging_adequacy.md),
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
and
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
leave such rows out too, so the steps of one pipeline agree on a cell's
mean. A function in `agg_funs` therefore never sees an infinite value:
`max` returns the largest finite value, and infinities cannot be counted
through `agg_funs`.

## Spatial autocorrelation and standard-error bias

By default (`deff = 1`), the `..se_*` columns are computed as
`sd / sqrt(n)`, which treats the observations within each cell as
independent. When data are spatially autocorrelated (the common case for
the spatial workflows this package supports), within-cell observations
are typically positively correlated: they share the cell's departure
from the population mean, so as an estimate of the **population (grand)
mean** a cell mean has an effective sample size smaller than `n`. For
that estimand the naive SE is **anticonservative** (too small), and a
downstream weighted regression that uses `cell_weight` or the `..se_*`
columns for population-level inference will produce overconfident
standard errors for cells with strong intra-cell correlation. For the
cell's **own** mean the naive SE is the right one when the cell's points
are spread through it, and the corrected SE is too wide; see "What the
standard error estimates" before setting `deff`.

Setting `deff = "kish"` applies an approximate correction using Kish's
design effect. An intra-class correlation (ICC) is estimated for every
summarised column on its own, the response and each predictor, via a
one-way random-effects decomposition across all cells. Each column's ICC
is used for that column's SE adjustment, and each cell's effective
sample size is reduced to `n_i / (1 + (n_i - 1) * rho)`. This is a
first-order correction that does not require a full spatial covariance
model but does require enough cells and observations for a stable ICC
estimate.

A design effect estimated from the data (`"kish"` or `"variogram"`)
comes with a second, small-sample correction. The within-cell
correlation that inflates the variance of the mean to
`sigma^2 * deff / n` also biases the within-cell sample variance
downward: under exchangeable correlation `rho` (Kish's own assumption),
with `deff = 1 + (n - 1) * rho`,
`E[s^2] = sigma^2 * (n - deff) / (n - 1)`, so `s^2` understates
`sigma^2` by very nearly the factor by which `deff` inflates the mean's
variance, and the two errors compound rather than cancel. The standard
error is therefore `s * sqrt(deff / n) * sqrt((n - 1) / (n - deff))`,
and `NA` where `deff >= n` (the cell then holds one observation's worth
of information and `s` carries none about `sigma`). This is the
package's own derivation, not taken from a reference. Measured 95%
interval coverage at `n = 30` over 20,000 replicates: 0.921, 0.844 and
0.628 at `rho` = 0.2, 0.5 and 0.8 with `s * sqrt(deff / n)` alone,
against 0.948, 0.950 and 0.949 with the rescaling.

You may also pass a fixed numeric design effect (e.g. `deff = 2`) to
uniformly inflate standard errors: an externally supplied constant is
applied as `sd * sqrt(deff / n)`, exactly `sqrt(deff)` times the naive
SE in every cell. (The `E[s^2]` correction that the estimated design
effects also apply is derived from within-cell correlation and would not
be justified for a number the caller chose.)

Even with the Kish correction, the adjusted SE is an approximation. For
rigorous inference under spatial dependence, consider fitting an
explicit spatial covariance model (e.g. via
[`fit_bayesian_spatial_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)).

## What the standard error estimates

What the `..se_*` columns estimate follows `deff`. At the default,
`deff = 1`, `sd / sqrt(n)` is the standard error of **the cell's own
mean**, the mean of the units the cell holds, and it is the right one
when the cell's observations are spread through the cell. Under
`deff = "kish"`, `"variogram"` or a number other than 1 the columns are
the standard error of the cell mean **as an estimate of the population
(grand) mean**: the unconditional quantity, in which the cell's own
realised deviation is part of the error, widened for the correlation
inside the cell. That is the right quantity when cells are treated as
samples from a common population: measured 95% interval coverage of the
grand mean is 0.95 with `deff = "kish"` (and 0.29 with the naive SE) on
exchangeable within-cell correlation, and 0.93 with `deff = "variogram"`
on a simulated Gaussian field.

The corrected standard error is **not** that of the cell's own mean (the
block average over that cell), which is what a cell-level map or a
regression on cell values usually wants. For that quantity the naive
`sd / sqrt(n)` is the better of the two on offer: measured coverage 0.95
under exchangeable within-cell correlation, against very nearly 1.00 for
the design-effect-corrected SE, which is too wide by the factor
`sqrt(deff / (1 - rho))` (4.6 at 20 points a cell and `rho = 0.5`). That
0.95 is exact under the exchangeable model and holds under a spatial
covariance model only when the cell's points are spread through the
cell; with *clustered* sampling inside a cell it is anticonservative for
the block average too (measured 0.58), and the right quantity there is a
block-kriging variance, which this function does not compute. Use `deff`
when the cell means feed a population-level inference; leave it at 1
when they are measurements of the cells themselves and the sampling
within cells is reasonably uniform.

## Design effects and variable types

`deff = "kish"` estimates a separate ICC for each summarised column, the
response and every predictor, and applies each to that column alone:
listing or removing a predictor leaves every other column's standard
error, effective n and interval as they were. (One ICC used to be pooled
over all the predictors, which was right for a column only when they all
shared it: beside a predictor with an ICC of 0.87, one with an ICC of 0
had its standard error multiplied by 3.1.) `deff = "variogram"` fits or
accepts **one** correlation function and applies it to every numeric
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
as a residual one from `estimate_sac_range(..., predictor_vars = )`. A
`sac` whose `attr(sac, "detrended")` is `TRUE` is used as given, but
when it corrects response columns a warning says that their standard
errors are understated, as
[`kriging_adequacy()`](https://elkronos.github.io/gis_modeling_toolkit/reference/kriging_adequacy.md)
warns about the same mismatch.

## Confidence intervals

With `conf_level` set, every numeric response and predictor column gains
four more columns: `..neff_*`, that column's effective sample size in
the cell (its non-missing count over its design effect, the per-column
version of `cell_weight`); `..df_*`, the degrees of freedom the interval
uses; and `..ci_lo_*` / `..ci_hi_*`, a t interval for a mean. Which mean
follows `deff`, as for `..se_*`: at the default, `deff = 1`, **the
cell's own mean**, when its observations are spread through the cell;
under `deff = "kish"`, `"variogram"` or a number other than 1, the
**cell mean as an estimate of the population (grand) mean**. Everything
in "What the standard error estimates" applies to it, including that the
corrected interval is not one for the cell's own block average. The
interval is `mean +/- qt((1 + conf_level) / 2, df) * se`, centred on the
plain mean of the column's non-missing values whatever `agg_funs`
computes, and is `NA` wherever the standard error is (a single
observation; complete redundancy under `deff`). `..neff_*` and `..df_*`
are `NA` where the column has one non-missing value or none, like the
standard error; `cell_weight` still counts such a cell (1, or `1 / deff`
for a numeric `deff`), so the two differ there.

At no setting is `..ci_*` a range for one observation: it is an interval
for a mean, and it narrows as a cell's count grows where the spread of
single observations does not. The interval for one further observation
in a cell has its own argument and columns, `obs_level` with
`estimate = "pooled"` and `..obs_lo_*` / `..obs_hi_*` (see "Intervals:
which one answers which question"). In the simulation described there,
`..ci_*` at `conf_level = 0.90` held 58.9 percent of the observations
held out among the others, and the interval for one observation at
`obs_level = 0.90` held 89.6 percent. Away from the observations no
interval here is calibrated.

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

## Pooled values

With `estimate = "pooled"` every numeric column also gets a pooled
value. That is the default; `estimate = "mean"` leaves the pooled
columns out and returns the frame of earlier versions. A cell's value,
plain or pooled, estimates the cell's *level*: the expected value of one
further observation located in the cell the way its observations were.
It is not the average of a surface over the cell's area, which would
give weight to ground that holds no observation.

For one column, write `n`, `S` and `m = S / n` for the count, the sum
and the mean of a cell's rows, `N` for the rows in all and `K` for the
cells that hold one. Then

- `sigma2 = sum((y - m)^2) / (N - K)` is the spread of the observed
  values about their cell's mean;

- the target `t` of a cell is the mean of the rows of its neighbours,
  `sum(S) / sum(n)` over them. The neighbours are the 8 cells whose
  centroids are nearest, the centroid of a cell being the mean of the
  coordinates of its rows. The polygons are not read, so two cells are
  neighbours whether or not they touch, and a cell is never its own
  neighbour;

- `tau2 = max(0, (sum(n * (m - t)^2) - K * sigma2) / N)` is the mean
  squared distance between a cell's level and its target, the target's
  own error included;

- the own weight is `w = tau2 / (tau2 + sigma2 / n)`;

- the pooled value is `w * m + (1 - w) * t`.

The own weight is the share of the pooled value that is the cell's own
mean. It grows with the cell's count and with `tau2`. At 1 the pooled
value is the plain mean; at 0 it is the target and holds none of the
cell's own values. Pooled values are for maps and for predicting
observations. `..pooled_w_*` is not a regression weight. The plain
columns are the same under either `estimate`, and they are the ones to
use, with `cell_weight`, in a regression on cell values: each column is
pulled toward its neighbours by its own amount, so a regression of
pooled responses on pooled predictors is not a regression of the levels.

The cells can be any polygons the rows were assigned to, the package's
own or ones the caller brings, with IDs of any type: only each row's
cell and location are read. Each column is pooled on its own rows, those
with a value, a cell and a location, with its own centroids, neighbours,
`sigma2` and `tau2`; a missing value in one column changes nothing in
another. A cell with one such row gets a value. A cell with none gets
`NA`, never a value made from its neighbours. The pooled columns do not
read `agg_funs` or `deff`: the mean that is pooled is the arithmetic
mean whatever `agg_funs` computes, and no design effect enters `sigma2`,
`tau2`, the weights or the interval of `obs_level`.

Three cases change what the numbers rest on, and each is a warning (see
"Conditions"). With 2 to 9 cells the neighbours of a cell are all the
others, so its target is the mean of every row outside it, and the
overall mean lies between that target and the cell's own mean. A cell
whose own weight is below its share of the rows, `n / N`, then gets a
pooled value on the far side of the overall mean from its plain mean,
and two cells can change order. When few cells hold two rows, `sigma2`
is estimated from few values and every weight rests on it: `df` in
`attr(, "pooled")` is the number of values beyond one per cell. And when
`tau2` comes out at 0, every own weight is 0 and each pooled value is
its target; with two cells that is the other cell's mean, so the two
values change places.

A tie between two candidate neighbours goes to the centroid with the
lower x, then the lower y, and the centroids are summed with the rows in
coordinate order, so a cell's neighbours do not depend on the order of
the rows. The values then agree to rounding: a cell's sums are taken in
the order of its rows. The tie is between distances computed from the
coordinates in the CRS the distances are measured in (see `estimate`),
so on gridded locations the same layer in another CRS or unit can give a
cell another neighbour.

The figures on this page that describe pooled values and the two
intervals come from one simulation, `dev/sim-pooled-values.R` in the
package's repository, whose header holds its settings and seeds. It
draws 200 layers on a square of side 1,000. Each has 2,100 clustered
locations (four fifths of them round 25 parent points with a standard
deviation of 40, the rest uniform) and a response that is a field with
covariance `exp(-d / 100)` plus noise of variance 0.3. The 300 locations
nearest a random centre are held out as a block, 300 of the others are
held out at random, and the 1,500 left are the layer: 333 cells at 4.5
observations a cell from
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
with `partition = "density"`, `select_resolution(criterion = "size")`,
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
and
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md),
summarised with `estimate = "pooled"`. A held-out observation is given
the value of its cell by
[`cell_values_at()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cell_values_at.md).
At the observations held out at random the root mean squared error of
the pooled value was 0.96 of the plain mean's: 0.84 in cells that hold
one of the layer's observations, 0.95 in cells of 2 to 4, 0.98 in cells
of 5 to 9 and 0.99 in cells of 10 or more. In the held-out block it was
0.95.

`structured` in `attr(, "pooled")` says whether these cells predict a
left-out observation better than one mean for all rows. It is `FALSE`
when the leave-one-out error of one mean (`loo_one_cell`) is within 2
percent of that of the pooled values (`loo`) or below it, or when the F
test of the cell means against one mean has `p >= 0.05` (`structure_p`).
`FALSE` raises a warning of class `"spatialkit_no_structure"`, and the
values are returned. On a response without spatial structure they still
make a smooth map, because each is mostly the mean of its neighbours.
With white noise as the response on the simulation's 200 layers the
median own weight was 0.11, `tau2` was 0 on 20 layers, the warning was
raised on all 200, and a cell's pooled value had a correlation of 0.79
with the mean of the pooled values of its 8 nearest cells, where the
plain means had one of -0.02.

The pooling model implies an error variance for a pooled value as an
estimate of the cell's level, `(1 - w) * tau2`. It is not returned, and
neither is an interval for the level, because no check has confirmed an
interval built from it; `sigma2`, `tau2` and `..pooled_w_*` are in the
result for a reader who wants to form it. In the simulation, with a
cell's level taken as the mean of the noise-free field over further
locations drawn in it like the data's, the squared error of the pooled
value was 0.98 of `(1 - w) * tau2` among the observations and 6.9 times
it in the held-out block: away from the observations it is too small.

Beside a sharp boundary the pooled value can be further from a cell's
level than the plain mean, because its target takes in cells across the
boundary. A second simulation measured this on fields whose level jumps
at county boundaries: one level for each county, alone and with a smooth
field added, on two layouts of real locations (892 and 2,634 points,
each in a county), with noise of 5, 20 or 50 percent of an observation's
variance, 200 layers of each, five random folds, and on each fold the
cells the package's defaults build from the training observations. A
held-out observation is a boundary row when every training observation
of its cell lies in its own county and more than half of those of the
cell's neighbours lie in other counties. Boundary rows were 21.75 to
28.04 percent of the held-out observations. Over them the root mean
squared error of the pooled value against the noise-free field, as a
geometric mean over the layers of its ratio to the plain mean's, was
from 9.16 percent lower to 56.36 percent higher. On county levels alone
it was 43.21 and 56.36 percent higher with 5 percent noise, 11.23 and
20.58 percent higher with 20 percent, and from 3.79 percent lower to
1.47 percent higher with 50 percent; with the smooth field added it was
from 9.16 percent lower to 1.59 percent higher. Over all held-out
observations the pooled value's error was the lower of the two in every
setting, by 1.75 to 18.73 percent.

## Intervals: which one answers which question

Two families of interval columns can be in the frame. Each has its own
argument, and neither changes the other.

|  |  |  |  |  |
|----|----|----|----|----|
| Family | Columns | Target | Centre | Asked for by |
| observation | `..obs_lo_*`, `..obs_hi_*` | one further observation in the cell | pooled value | `obs_level` |
| mean | `..ci_lo_*`, `..ci_hi_*` | a mean, never one observation | plain mean | `conf_level` |

The observation of the first is one located in the cell the way the
data's were. `obs_level` needs `estimate = "pooled"`: without it the
call is an error, not an interval round the plain mean. The mean the
second covers follows `deff`: at `deff = 1` the cell's own mean, when
its observations are spread through the cell; under `deff = "kish"`,
`"variogram"` or a number other than 1 the cell mean as an estimate of
the population mean (see "Confidence intervals", which also describes
`..neff_*` and `..df_*`). No interval and no standard error is returned
for a cell's level as the pooled value estimates it.
`attr(, "intervals")` says the same of the frame in hand, one row per
family present.

The interval for one observation is the pooled value plus two quantiles
of the leave-one-out residuals of the pooled values. The residual `e` of
a row is its value minus the pooled value of its cell computed without
it: the cell's mean and own weight are computed again without the row,
and `sigma2`, `tau2` and the targets are held. With
`p = (1 - obs_level) / 2`, `..obs_lo = pooled + quantile(e, p)` and
`..obs_hi = pooled + quantile(e, 1 - p)`.

- It is a statement about observations on average: of further
  observations located like the data, about `obs_level` of them fall
  inside the interval of their cell. It does not hold cell by cell (the
  second simulation below gives the coverage in cells of one
  observation).

- The two offsets are the same in every cell. The interval's centre
  varies over the map and its width does not. It does not narrow where a
  cell holds more rows: the offsets are quantiles of the errors of
  predicting single observations, and the spread of observations inside
  a cell is in every one of them.

- It is wide. One observation is uncertain.

- Only the level 0.90 has been checked. Another level is computed the
  same way.

- Each residual is the error of predicting an observation from the other
  rows of its cell and from the cells around it, so the interval
  describes locations among the observations. A cell's value and
  interval describe the part of the cell where its observations lie, and
  a Voronoi cell can reach far past them toward the edge of the study
  area. The interval is too narrow for a cell with no observation, a gap
  in the layer or ground beyond the layer's edge.

In the simulation described under "Pooled values", the interval at
`obs_level = 0.90` held 89.6 percent of the observations held out at
random among the others (83.7 to 94.0 percent in single layers) and 75.1
percent of those in the held-out block. `..ci_*` at `conf_level = 0.90`,
read as a range for one observation, held 58.9 and 49.2 percent of them,
over the held-out observations whose cell has that interval.
[`cell_values_at()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cell_values_at.md)
returns each location's distance from the nearest observation, and its
help gives the coverage by that distance.

The second simulation described under "Pooled values" read the interval
at `obs_level = 0.90` on smooth fields of two ranges and on county
levels with a smooth field, with 20 or 50 percent noise, on both
layouts: twelve settings. Setting by setting it held 89.38 to 89.79
percent of the held-out observations. Counted by the training
observations in an observation's cell (1, 2 to 4, 5 to 9, 10 or more),
it held at least 85.65 percent in every class but one: on county levels
with a smooth field, on the layout of 892 points with 20 percent noise,
it held 83.74 percent of the 4,878 observations whose cell held one
training observation. On fields of that kind it is not a 90 percent
interval in such cells: there it held 83.74 to 87.83 percent with 20 or
50 percent noise, and 82.64 and 82.94 percent with 5 percent.

Away from the observations no interval here is calibrated.

## Conditions

What the pooled columns and the count record have to say is raised as R
warnings with a class, so `tryCatch(<class> = )` and
[`withCallingHandlers()`](https://rdrr.io/r/base/conditions.html) catch
each without matching its text. A warning of the pooled columns is
raised once per call and reason, for all the columns it applies to: its
text names them and the condition carries them in `columns`, in the
order of the pooled columns. A field said to be by column is a vector
named by those columns, and the numbers of each column are in
`attr(, "pooled")`.

- `"spatialkit_pool_unavailable"`: the pooled columns of a column are
  `NA`. The plain columns are returned, and this is never an error.
  `reason` is `"no coordinates"` (a layer without geometry),
  `"no located value"` (no row with a value, a cell and a location),
  `"no cell holds two values"` (the condition also carries `n` and
  `cells` by column) or `"work limit"`. The last is the case of many
  cells without the 'FNN' package: the neighbour search then compares
  every pair of cells and is not made above
  `getOption("spatialkit.exact_pairs_max", 2e9)` pairs, which is about
  44,700 cells (the condition also carries `work`, `limit`, and `cells`
  by column). With 'FNN' there is no such limit, and the result is the
  same.

- `"spatialkit_pool_degenerate"`: the pooled values were computed and
  carry no pooling. `reason` is `"one cell"` (the pooled value is the
  plain mean; also `n` by column), `"no spread inside any cell"` (every
  own weight is 1; also `constant` by column) or
  `"no signal between cells"` (`tau2` is 0 and every own weight is 0;
  also `sigma2` by column).

- `"spatialkit_pool_few_cells"`: 2 to 9 cells, so each target is the
  mean of all the other cells (`cells`, and `own_weight`, the smallest,
  by column).

- `"spatialkit_pool_few_df"`: `sigma2` rests on fewer than 30 values
  beyond one per cell (`needed`, and `df`, `n` and `cells` by column).

- `"spatialkit_no_structure"`: `structured` is `FALSE` (see "Pooled
  values"). The text says which of the two tests failed, and `beats` is
  `TRUE` for a column whose pooled values beat one mean while the F test
  does not tell the cell means from one mean (also `criterion = "loo"`,
  `estimate = "pooled"`, and by column `n`, the cells, `one_cell`, `loo`
  and `structure_p`).

- `"spatialkit_pool_rows_skipped"`: rows with a value and a cell have
  empty or non-finite coordinates. They count in `n` and the plain
  columns and take no part in the pooled ones (`n`, `cells`).

- `"spatialkit_repeated_locations"`, with `obs_level`: rows with a value
  share a location. A residual leaves one row out and keeps the others
  at its location, so the interval then describes a further observation
  at a place already observed (`n` and `total` by column).

- `"spatialkit_few_residuals"`, with `obs_level`: fewer than five
  residuals are expected beyond each end of the interval, which at 0.90
  is fewer than 100 residuals (`needed`, `level`, and `n` by column).

- `"spatialkit_estimate_mismatch"` and `"spatialkit_cells_differ"`: the
  two checks of the count record, described under "Value".

The two about few cells and few values are raised for columns whose
weights were estimated, and the two of `obs_level` for columns whose
interval is read from residuals, so a constant column raises
`"spatialkit_pool_degenerate"` and no other.

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
[`cell_values_at()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cell_values_at.md),
[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md),
[`kriging_adequacy()`](https://elkronos.github.io/gis_modeling_toolkit/reference/kriging_adequacy.md),
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md),
[`summary.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.resolution_profile.md)

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
             val = 5 + 0.02 * east + 0.02 * north + rnorm(n, sd = 0.5)),
  coords = c("x", "y"), crs = 32632
)
bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
  c(5e5, 5e6), c(5e5 + 100, 5e6), c(5e5 + 100, 5e6 + 100),
  c(5e5, 5e6 + 100), c(5e5, 5e6)
))), crs = 32632))
grid <- create_grid_polygons(bnd, target_cells = 9, type = "square")
assigned <- assign_features_to_polygons(pts, grid)

# Nine cells are too few for pooled values (each would be pulled toward
# all the other cells), so these summaries report plain means alone:
# estimate = "mean".
# IID standard errors (default) vs Kish design-effect adjustment.  The
# correction is large here, and meant to be: about 22 points per cell that
# all share the cell's part of the trend carry far fewer than 22
# independent pieces of information about the cell mean.
naive <- summarize_by_cell(assigned, response_var = "val", estimate = "mean")
kish  <- summarize_by_cell(assigned, response_var = "val", deff = "kish",
                           estimate = "mean")
data.frame(n = naive$n,
           se_naive = naive$..se_resp_val,
           se_kish  = kish$..se_resp_val)
#>    n   se_naive   se_kish
#> 1 25 0.09706775 0.7210227
#> 2 28 0.11815220 0.9279034
#> 3 25 0.11573069 0.8596517
#> 4 21 0.11556510 0.7881134
#> 5 23 0.13192931 0.9407000
#> 6 21 0.12707224 0.8665880
#> 7 10 0.18215988 0.8673238
#> 8 27 0.09534639 0.7355266
#> 9 20 0.14402822 0.9590654
attr(kish, "deff_applied")   # method, icc_resp, icc_pred, per-cell deff
#> $method
#> [1] "kish"
#> 
#> $icc_resp
#> [1] 0.6842466
#> 
#> $icc_pred
#> [1] NA
#> 
#> $deff
#> [1] 17.42192 19.47466 17.42192 14.68493 16.05343 14.68493  7.15822 18.79041
#> [9] 14.00069
#> 

# A 95% interval for each cell mean as an estimate of the grand mean, on
# the Kish-corrected standard error and n - 1 degrees of freedom.
ci <- summarize_by_cell(assigned, response_var = "val", deff = "kish",
                        conf_level = 0.95, estimate = "mean")
ci[, c("poly_id", "n", "resp_mean_val", "..neff_resp_val", "..df_resp_val",
       "..ci_lo_resp_val", "..ci_hi_resp_val")]
#> # A tibble: 9 × 7
#>   poly_id     n resp_mean_val ..neff_resp_val ..df_resp_val ..ci_lo_resp_val
#>     <int> <int>         <dbl>           <dbl>         <dbl>            <dbl>
#> 1       1    25          5.63            1.43            24             4.14
#> 2       2    28          6.30            1.44            27             4.40
#> 3       3    25          6.93            1.43            24             5.16
#> 4       4    21          6.42            1.43            20             4.77
#> 5       5    23          7.11            1.43            22             5.15
#> 6       6    21          7.89            1.43            20             6.09
#> 7       7    10          7.16            1.40             9             5.19
#> 8       8    27          7.77            1.44            26             6.25
#> 9       9    20          8.16            1.43            19             6.15
#> # ℹ 1 more variable: ..ci_hi_resp_val <dbl>

# Pooled values, and a range for one observation: a smooth field plus
# noise at 400 points, on an 8 x 8 grid, so that a cell holds about six.
set.seed(2)
east  <- runif(400, 0, 8000)
north <- runif(400, 0, 8000)
pts2 <- st_as_sf(
  data.frame(x = 5e5 + east, y = 5e6 + north,
             val = 10 + 3 * sin(east / 1500) + 3 * cos(north / 1200) +
               rnorm(400)),
  coords = c("x", "y"), crs = 32632
)
bnd2 <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
  c(5e5, 5e6), c(5e5 + 8000, 5e6), c(5e5 + 8000, 5e6 + 8000),
  c(5e5, 5e6 + 8000), c(5e5, 5e6)
))), crs = 32632))
grid2 <- create_grid_polygons(bnd2, target_cells = 64, type = "square")
assigned2 <- assign_features_to_polygons(pts2, grid2)
# estimate = "pooled" is the default.
pooled <- summarize_by_cell(assigned2, response_var = "val", obs_level = 0.90)
# The plain mean, the pooled value and the weight on the cell's own mean.
head(pooled[, c("poly_id", "n", "resp_mean_val", "..pooled_resp_val",
                "..pooled_w_resp_val")])
#> # A tibble: 6 × 5
#>   poly_id     n resp_mean_val ..pooled_resp_val ..pooled_w_resp_val
#>     <int> <int>         <dbl>             <dbl>               <dbl>
#> 1       1     3          13.6              13.3               0.652
#> 2       2    15          15.0              14.8               0.904
#> 3       3     2          16.6              15.4               0.555
#> 4       4     2          14.0              13.8               0.555
#> 5       5     8          13.4              13.3               0.833
#> 6       6     9          10.9              10.9               0.849
# The range for one further observation in a cell.
head(pooled[, c("poly_id", "..obs_lo_resp_val", "..obs_hi_resp_val")])
#> # A tibble: 6 × 3
#>   poly_id ..obs_lo_resp_val ..obs_hi_resp_val
#>     <int>             <dbl>             <dbl>
#> 1       1             11.3               15.4
#> 2       2             12.8               16.9
#> 3       3             13.4               17.5
#> 4       4             11.8               15.9
#> 5       5             11.3               15.4
#> 6       6              8.86              12.9
attr(pooled, "pooled")   # what was fitted, and what the interval covers
#> Cell values: val, 400 rows in 64 cells
#>   value       : ..pooled_resp_val, pooled toward the 8 nearest cells (centroids of the rows)
#>   variances   : inside cells 1.41e+00 (from 336 values beyond one per cell);
#>                 between a cell and its target 8.80e-01
#>   own weight  : median 0.79 (..pooled_w_resp_val)
#>   alone       : 0.5% of the rows are the only one in their cell
#>   loo         : 1.61e+00 (leave-one-out mean squared error of the pooled value);
#>                 one mean for all rows: 1.08e+01
#>   interval    : 90% for ONE OBSERVATION (..obs_lo_resp_val, ..obs_hi_resp_val):
#>                 value - 2 to value + 2.09, quantiles of 400 leave-one-out residuals
#> Away from the observations no interval here is calibrated.
```
