# Options that bound the work of a search or a fit

The package reads two options. Each bounds a step whose cost grows with
the square of the rows or of the cells, so that such a step stops, or is
skipped, before it starts. Set them with
[`options()`](https://rdrr.io/r/base/options.html).

## The options

- `spatialkit.exact_pairs_max`:

  Default `2e9` pairs. The most pairs a search made without the FNN
  package may compare. Without FNN, a nearest-neighbour search (each
  cell's neighbours, each row's nearest centre, each location's nearest
  other location) compares every pair, and a search over the pairs of
  about 44,700 cells or locations reaches the default. With FNN, which
  installs with spatialkit, these searches use a kd-tree and the limit
  is not read. The result of every search is the same with and without
  FNN. Above the limit:

  - [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
    stops with an error of class `spatialkit_work_limit` before the
    search: on a ladder of sizes, for the weights of
    `partition = "density"`; with an `estimate`, on any ladder, for the
    estimate's searches.

  - [`cell_values_at()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cell_values_at.md)
    stops with the same error before the search.

  - [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
    returns `NA` pooled columns with a warning of class
    `spatialkit_pool_unavailable` (`reason = "work limit"`), and does
    not run the check of the count record (`"not checked: work limit"`).

  - [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
    hands on the centres fitted to the subsample without refitting them
    to every point, with a logged warning.

- `spatialkit.fit_work_max`:

  Default `1e10` point-to-centre comparisons. The most work the
  partition fits of a ladder of sizes may do: the points times the cells
  of the finest level times `nstart`, times `iter_max` for
  `partition = "density"` without FNN. With FNN and `init = "random"`
  the density kind has no step that compares every point with every
  centre, and nothing is counted. Above the limit,
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  stops with an error of class `spatialkit_work_limit` before any fit.
  At the default size the work grows with the square of the rows.
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  refits the centres of a profile fitted to a subsample to every point,
  and counts the points times the cells of the count it selects, times
  `iter_max` for the density kind without FNN, whatever the start. Above
  the limit it hands on the centres fitted to the subsample, with a
  logged warning.

## Setting a limit

Each option takes a single number, 0 or more; `Inf` removes the limit.
Any other value is an error when the option is read. The error of class
`spatialkit_work_limit` carries `work` (the count), `limit` and `what`
(the step), and its message says what to change: for a fit, fewer
starts, a coarser `obs_per_cell`, random starts of the density kind or a
higher limit; for a search, reinstalling FNN. A profile just below
`spatialkit.fit_work_max` prints nothing while its fits run.

The console log is set with
[`spatialkit_quiet()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md),
not with an option.

## See also

[spatialkit-conditions](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-conditions.md)
for the classes named here, and
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
for the work its fits do.

Other options and caches:
[`clear_fitted_cache()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_fitted_cache.md),
[`clear_grid_cache()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_grid_cache.md),
[`spatialkit-conditions`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-conditions.md),
[`spatialkit_quiet()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)

## Examples

``` r
getOption("spatialkit.exact_pairs_max", 2e9)
#> [1] 2e+09
getOption("spatialkit.fit_work_max", 1e10)
#> [1] 1e+10

# Raise one limit for a session, then put it back.
old <- options(spatialkit.fit_work_max = 2e10)
getOption("spatialkit.fit_work_max")
#> [1] 2e+10
options(old)
```
