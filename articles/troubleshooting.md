# Troubleshooting and performance

This article lists the messages you are most likely to meet, what each
means and what to do about it. It then covers large layers and long
runs: the work limits, parallel folds, the two caches and the log.
[`?"spatialkit-conditions"`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-conditions.md)
lists every warning and error class the package raises, with the fields
a handler can read.

Each message below is quoted as R prints it, up to the point where it
starts to vary. Angle brackets stand for the names and counts that
differ from call to call.

## A layer for the examples

Two examples below run on a simulated layer: 300 points over a 20 km
square in UTM zone 32N, whose unit is the metre.

``` r

suppressPackageStartupMessages(library(sf))
library(spatialkit)

set.seed(1)
n   <- 300
pts <- st_as_sf(data.frame(x = 5e5 + runif(n, 0, 20000),
                           y = 5e6 + runif(n, 0, 20000)),
                coords = c("x", "y"), crs = 32632)
pts$w <- rnorm(n)
pts$z <- 2 * pts$w + rnorm(n)
st_crs(pts)$units_gdal
#> [1] "metre"
```

## Messages you are likely to meet

### Preparing the data

- `` build_tessellation(): `boundary` is required for hex/square grids. ``:

  A hex or square grid has no extent of its own. Pass `boundary =` a
  polygon of the study area:
  [`st_union()`](https://r-spatial.github.io/sf/reference/geos_combine.html)
  of your source polygons is one way to get it, and
  [`clip_target_for()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clip_target_for.md)
  builds the bounding box of the points. Voronoi cells and `"triangles"`
  do not need a boundary, since they take their extent from the points,
  or from `cover` when the points are seeds.

- `prep_model_data(): missing required column(s): <names>`:

  A name in `response_var` or `predictor_vars` is not a column of the
  data. The usual causes are a typo, a difference in case, or a column
  that [`read.csv()`](https://rdrr.io/r/utils/read.table.html) renamed
  under `check.names = TRUE` (`pop density` becomes `pop.density`).
  Check `names(data_sf)`. The geometry column is not a predictor.

- `prep_model_data(): column name(s) '<names>' are not syntactically valid R names`:

  Every backend builds a formula from the column names, and a name such
  as `B5-B4` or `band 4` does not survive that. Rename the columns
  before fitting, for example with `names(x) <- make.names(names(x))`.

- `prep_model_data(): response '<name>' is also listed in 'predictor_vars'.`:

  A model cannot use its own response as a predictor, and none of the
  backends would stop it, so
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
  refuses the call before any fit. Drop the response from
  `predictor_vars`.

### Fitting a model

- `fit_rf_model(): response '<name>' is not numeric.`:

  The random forest and GWR backends fit a regression, so they refuse a
  logical, factor or character response. The GWR message starts
  `fit_gwr_model(): response '<name>' is not numeric (it is <class>).` A
  number read from a file as text needs
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html) first; check for
  a thousands separator or an `"NA"` string if that gives `NA`s.
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  also refuses a response with two whole-number values, and names
  [`GWmodel::ggwr.basic()`](https://rdrr.io/pkg/GWmodel/man/ggwr.basic.html)
  for it.
  [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
  takes a factor or character response under
  [`brms::categorical()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
  or an ordinal family, and refuses it under any other family with a
  message that starts
  `fit_bayesian_spatial_model(): response '<name>' is <class>, not numeric`.

- `` fit_rf_model(): `num.trees` is already set by this function and cannot be passed through `...` as well ``:

  ranger’s own spelling of an argument that
  [`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
  sets would reach
  [`ranger::ranger()`](http://imbs-hl.github.io/ranger/reference/ranger.md)
  twice. Use the wrapper’s spelling, which the message names, and it
  passes the value on, for example `num_trees` for `num.trees` or
  `sample_fraction` for `sample.fraction`.
  [`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
  and `compare_models_cv(rf_args = )` pass their arguments to
  [`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md),
  so the same rule applies there.

### Units and distances

Distances, bandwidths and block sizes are in the units of the CRS a
function works in. When one looks absurd, check that unit first with
`st_crs(x)$units_gdal`. A block size of 1,000 meant in metres is about
305 m on a layer in US survey feet, and a layer in lon/lat is moved into
a projected CRS that
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
chooses, whose unit is the metre.
[`?"spatialkit-geometry"`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)
gives the rules, and
[`vignette("getting-started")`](https://elkronos.github.io/gis_modeling_toolkit/articles/getting-started.md)
shows them on a real layer.

- `make_folds(block_kfold): the requested grid is <nx> x <ny> = <n> cells, above the 1,000,000 this function will build.`:

  A `block_size` in the wrong unit asks for a grid far finer than meant.
  The message gives the extent and the CRS unit to compare it with:

  ``` r

  # A block size of 10, meant in kilometres, on a layer in metres.
  # The message is wrapped to fit the page.
  tryCatch(
    make_folds(pts, k = 5, method = "block_kfold", block_size = 10),
    error = function(e) writeLines(strwrap(conditionMessage(e), 76))
  )
  #> make_folds(block_kfold): the requested grid is 1959 x 1988 = 3,894,492
  #> cells, above the 1,000,000 this function will build. Check that
  #> `block_size` (10) is expressed in the data's CRS units (metre) over an
  #> extent of 19590 x 19880; a value in the wrong unit is the usual cause.
  ```

  [`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md)
  refuses a grid of more than 5,000,000 cells for the same reason, and
  “Large layers” below lists the other limits.

### Cross-validation

- `cv_*(): all folds failed (all <n> folds failed to produce predictions); cross-validation results contain no predictions.`:

  The call returns with this warning, `$overall` all `NA` and
  `n_pred = 0`, and the log line of each fold names its cause. A missing
  backend is one cause:
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
  without brms fails every fold this way, while
  [`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
  and
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md)
  stop at once when ranger or sp is missing. A fit that stops on one
  fold’s training rows fails that fold the same way. Here the learner
  stops on every fold:

  ``` r

  folds <- make_folds(pts, k = 5, method = "block_kfold", seed = 1)
  cv <- cv_spatial(pts, "z", "w", folds = folds,
                   fit_fn = function(train) stop("the backend is not installed"))
  #> WARN [2026-11-12 12:26:23] cross-validation: fold 1 fit failed; skipping. Cause: the backend is not installed
  #> WARN [2026-11-12 12:26:23] cross-validation: fold 2 fit failed; skipping. Cause: the backend is not installed
  #> WARN [2026-11-12 12:26:23] cross-validation: fold 3 fit failed; skipping. Cause: the backend is not installed
  #> WARN [2026-11-12 12:26:23] cross-validation: fold 4 fit failed; skipping. Cause: the backend is not installed
  #> WARN [2026-11-12 12:26:23] cross-validation: fold 5 fit failed; skipping. Cause: the backend is not installed
  #> Warning: cv_spatial(): all folds failed (all 5 folds failed to produce
  #> predictions); cross-validation results contain no predictions. First error: the
  #> backend is not installed
  c(attempted = cv$n_folds_attempted, succeeded = cv$n_folds_succeeded)
  #> attempted succeeded 
  #>         5         0
  cv$overall
  #>   RMSE MAE MAPE SMAPE R2 Adj_R2 n_pred n_MAPE n_SMAPE
  #> 1   NA  NA   NA    NA NA     NA      0      0       0
  ```

- `` cv_*(): <k> of <n> fold(s) failed (<causes>), so `overall` pools the other folds only ``:

  Some folds failed and the rest were scored. `$overall` then pools the
  folds that ran, and covers fewer rows than the data hold;
  `fold_status` gives the cause of each failure. Every `cv_*()` result
  records `n_folds_attempted` and `n_folds_succeeded`, so compare the
  two before you read `$overall`, on runs that look fine as well.

- `compare_models_cv(): no recognised model requested.`:

  None of the names in `models` is `"GWR"`, `"Bayesian"` or `"RF"`,
  which are matched exactly, case included. The warning
  `ignoring unknown model(s)` printed with the error names the ones it
  ignored.

- `compare_models_cv(): no viable models.`:

  At least one name was recognised, but no recognised backend is
  installed. Each prints
  `dropping <name> (package/function unavailable)` before the error, and
  any name not recognised is named in a warning. Install the backend, or
  ask for one you have.

### Autocorrelation and resolution

- [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  returned `NA` :

  The range was not identified, so none is reported.
  `attr(x, "rejected_reason")` says why, and
  [`?estimate_sac_range`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  lists the reasons. When the variogram was computed,
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on the result
  ([`?plot.sac_range`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.sac_range.md))
  draws it, so you can see why the range was refused. An `NA` returned
  before the variogram was computed (gstat missing, or fewer than 30
  points or locations, for example) carries the reason alone.
  [`vignette("spatial-cross-validation")`](https://elkronos.github.io/gis_modeling_toolkit/articles/spatial-cross-validation.md)
  covers choosing a block size.

- `determine_optimal_levels(): the model-aware criteria score only the elbow's neighbourhood, k = <a> to <b>, and carry no information at nine cells or fewer`:

  A logged warning. Every candidate around the elbow has nine cells or
  fewer, where Moran’s I carries no information, so the call falls back
  to the geometric ranking. On points with no cluster structure, under
  `on_no_elbow = "warn"`, the elbow sits near `sqrt(max_levels)`, and
  the neighbourhood reaches ten cells from about `max_levels = 40`. At
  the default, `on_no_elbow = "na"`, such points return `NA` with a
  warning of class `spatialkit_no_elbow` before any neighbourhood is
  scored. A line logged before the k-means sweep says the same when the
  ladder itself ends at nine cells or fewer.
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  scores Moran’s z at every level of its ladder; see
  [`vignette("resolution")`](https://elkronos.github.io/gis_modeling_toolkit/articles/resolution.md).

- `determine_optimal_levels(): Moran's I could not be computed at any k from <a> to <b> in the elbow's neighbourhood`:

  A logged warning. The neighbourhood reaches past nine cells, but too
  few cells hold a row with a response and every predictor, the
  regression of the cell means on the predictors is singular (a
  predictor constant or collinear across cells), or the cell means
  differ by rounding alone. Check for missing values first. Under
  `criterion = "morans_i"` the call then falls back to the geometric
  ranking; under `"combined"` it ranks the counts as usual, with
  `moran_i` and `moran_z` `NA` in the diagnostics.

### Checking a fit

- `n = <n> requires FNN for k-NN weights, and Matrix to hold them sparsely`:

  [`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)
  on more than 5,000 observations needs FNN to find the neighbours and
  Matrix to store the weights; without them it would allocate a dense n
  by n matrix, so it stops instead. FNN installs with spatialkit, so the
  missing one is usually Matrix: `install.packages("Matrix")`. The
  message names the one it could not load.

## Large layers

Most steps take time and memory that grow with the rows or the cells.
The table lists where a function stops, subsamples or leaves a column
out on a large layer. A limit given as an argument is that argument’s
default.

| function | limit | at the limit | to change it |
|----|----|----|----|
| `make_folds(method = "block_kfold")`, and the `cv_*()` functions that build blocks | a grid of 1,000,000 blocks | error naming the extent and the CRS unit | a larger `block_size` in the CRS unit, fewer `block_nx` and `block_ny`, or a smaller `block_multiplier` |
| [`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md) | a grid of 1,000,000 blocks at a size | error for the sizes you pass; a size on the default ladder is dropped, with a log line | larger `block_sizes` |
| `make_folds(method = "nndm")` | 5,000 rows | error | `method = "knndm"`, which has no row limit, or `"block_kfold"` |
| [`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md) | 5,000,000 grid cells | error | a larger `cell_size` or a smaller `n_cells`; `chunk_size`, 5,000 rows, sets the rows per prediction call |
| [`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md), and [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md) for hex and square cells | `max_cells`, 1,000,000 cells, estimated from the bounding box | error | a larger `cellsize`; `create_grid_polygons(max_cells = )` raises the limit, which [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md) keeps at its default |
| [`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md) | 5,000 observations without FNN and Matrix; 20,000,000 weights (`k` times the observations) | error | install Matrix; a smaller `k` |
| [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md) | a level with more than 5,000 cells holding a scored row, without FNN and Matrix | `moran_i` and `moran_z` `NA` at that level, with a warning of class `spatialkit_moran_skipped` | install Matrix |
| [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md) on a ladder of sizes, the default | `spatialkit.fit_work_max`, 1e10 point-to-centre comparisons | error of class `spatialkit_work_limit` before any fit | a smaller `nstart`, a coarser `obs_per_cell`, `init = "random"`, or a higher `spatialkit.fit_work_max` |
| searches made without FNN in [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md), [`cell_values_at()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cell_values_at.md), [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md) and [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md) | `spatialkit.exact_pairs_max`, 2e9 pairs | by function: an error, `NA` pooled columns with a warning, or the subsample’s centres with a log line | reinstall FNN, or a higher `spatialkit.exact_pairs_max` |
| [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md) | `n_max`, 5,000 rows | a seeded random subsample, counted in the attributes `n_used` and `n_subsampled` | `n_max`, or `Inf` for no limit |
| [`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md) | `normalizer_max_n`, 5,000 training rows | the mean pairwise distance is taken on a subsample | `normalizer_max_n` |
| `summarize_by_cell(deff = "variogram")` | `deff_max_n`, 500 points in a cell | the cell’s correlation matrix is built on a seeded subsample of its points | `deff_max_n` |
| [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md), and [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md) on a ladder of counts | `sample_n`, 1,500 points | the k-means fits run on a subsample | `sample_n` |
| `fit_bayesian_spatial_model(gp_resolution = )` | 50 basis functions per dimension, 2,500 in all | error giving the basis the range needs | a longer `gp_resolution`; the cost grows with the square of `gp_k` |

[`?"spatialkit-options"`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-options.md)
describes the two options and what each step does above its limit.

## Parallel folds

[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md),
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
and
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
take `parallel`. `TRUE` uses the physical cores less one, or
`getOption("mc.cores")` when that is lower, and a whole number above 1
uses that many workers, up to the cores the machine has. The folds then
run in forked workers through
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) on
macOS and Linux. On Windows the call runs them one after another, with a
message.
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md)
takes the argument and runs its folds in turn, with a warning when more
than one worker is asked for. Parallel folds help most in
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
where every fold compiles and samples a Stan model.

``` r

cv <- cv_rf(pts, "z", "w", k = 5, parallel = 4L)       # four workers
cv <- cv_bayes(pts, "z", "w", k = 5, parallel = TRUE)  # detects the cores
```

The result depends on `seed` and the fold index alone. The call draws
one seed per fold in the parent process, so `parallel = TRUE` returns
what `parallel = FALSE` returns. With `seed = NULL` the per-fold seeds
come from the session’s random stream, and
[`set.seed()`](https://rdrr.io/r/base/Random.html) before the call makes
both runs reproducible.

## Caches

The package keeps two results in the session so that it does not compute
them twice.

[`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md)
keeps each grid it builds, keyed on the boundary’s geometry and CRS, the
cell count and the other arguments, and returns the kept grid when the
same call comes again. It holds 50 grids (`max_entries`) and drops the
oldest first;
[`clear_grid_cache()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_grid_cache.md)
empties it. Its grids pass through
[`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md),
which reorders and renumbers the cells, so the same cell can carry a
different `poly_id` from
[`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md).
Use one of the two builders throughout an analysis: a summary keyed on
IDs from one and joined to the geometry of the other can draw the values
on the wrong cells.

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) on a Bayesian
GP fit
([`?fitted.bayesian_fit`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-fitted.md))
keeps the posterior expectations it computes, in an environment the fit
carries, because [`summary()`](https://rdrr.io/r/base/summary.html)
([`?summary.spatial_fit`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.spatial_fit.md)),
[`residuals()`](https://rdrr.io/r/stats/residuals.html)
([`?residuals.bayesian_fit`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-residuals.md)),
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
and
[`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md)
each call it. `clear_fitted_cache(fit)` drops them. You need it only
after editing the `brmsfit` engine in place; a change to `data_sf`, or a
refit, invalidates the entry on its own.

## The log

The package keeps a log through `logger`, apart from R’s warnings. Lines
at the INFO level and above go to
`file.path(tempdir(), "spatialkit_model_log.log")`, a file that lasts as
long as the R session, and lines at WARN and above are also printed to
the console. The log has its own `"spatialkit"` namespace: it leaves
your global `logger` configuration alone, and a global configuration set
before the package loads does not carry over into it.

The last line of the file is the warning the failed cross-validation
above ended on:

``` r

log_file <- file.path(tempdir(), "spatialkit_model_log.log")
# The line is wrapped to fit the page.
writeLines(strwrap(tail(readLines(log_file), 1), 76))
#> WARN [2026-11-12 12:26:23] cv_spatial(): all folds failed (all 5 folds
#> failed to produce predictions); cross-validation results contain no
#> predictions. First error: the backend is not installed
```

[`spatialkit_quiet()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)
turns the console echo off and returns the threshold that was in force,
which puts it back:

``` r

old <- spatialkit_quiet()
spatialkit_quiet(old)
```

The file and the console are two `logger` appenders, index 1 and index
2.
[`logger::log_appender()`](https://daroczig.github.io/logger/reference/log_appender.html)
and
[`logger::log_threshold()`](https://daroczig.github.io/logger/reference/log_threshold.html)
change index 1, the file, unless you name index 2. Each line below is a
separate choice:

``` r

# Write the log file to my_analysis.log (index 1)
logger::log_appender(logger::appender_file("my_analysis.log"),
                     namespace = "spatialkit")
# Send the console lines to a file instead of the console (index 2)
logger::log_appender(logger::appender_file("my_console.log"),
                     namespace = "spatialkit", index = 2)
# Or print only FATAL lines on the console (index 2)
logger::log_threshold(logger::FATAL, namespace = "spatialkit", index = 2)
```

A log line is not an R condition: `tryCatch(warning = )` does not catch
it and [`suppressWarnings()`](https://rdrr.io/r/base/warning.html) does
not silence it. “Warnings and log lines” in
[`?"spatialkit-conditions"`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-conditions.md)
says how the help pages tell the two apart.

knitr does not capture what is written to the console’s error stream, so
while a document is knitted (R Markdown, Quarto, pkgdown) each console
line is also sent as an R message, as the examples above show. A line
the package also raises as an R warning is not repeated, and the chunk
option `message = FALSE` keeps the rest out of the document.

## Next

[`vignette("glossary")`](https://elkronos.github.io/gis_modeling_toolkit/articles/glossary.md)
comes next: short definitions of the terms the articles and help pages
use. This article relied on
[`?"spatialkit-conditions"`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-conditions.md),
[`?"spatialkit-options"`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-options.md),
[`?"spatialkit-geometry"`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md),
[`?spatialkit_quiet`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md),
[`?make_folds`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md),
[`?cv_spatial`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
[`?compare_models_cv`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md),
[`?estimate_sac_range`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
[`?determine_optimal_levels`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md),
[`?residual_morans_i`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md),
[`?predict_surface`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md),
[`?create_grid_polygons_cached`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md)
and
[`?clear_fitted_cache`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_fitted_cache.md).
