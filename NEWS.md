# spatialkit (development version)

## Breaking changes

Defaults and results that differ from 2.0.0.  Each entry says how to get the
old behaviour where there is a way.

* `resolution_profile()` scores pooled values by default (`estimate =
  "auto"`) at 4.5 observations per cell (`obs_per_cell = "auto"`), with the
  centres fitted to every row.  A call with a numeric response, no
  predictors and none of `levels`, `n_levels`, `min_cell_n`, `sample_n` and
  `range_floor` returns one level where it returned a ladder of about 20,
  and its count is larger: on 3,000 rows, 667 cells where the ladder stopped
  at 333.  `elbow` is `NA` on one level.  The profile gains
  `scored_per_cell`, `loo`, `loo_se`, `share_alone` and `own_weight`,
  `attr(, "estimate")` and a print line, and a layer whose scored rows share
  locations raises `spatialkit_repeated_locations`.  A call that sets one of
  those five arguments keeps the ladder they describe and gains the same
  columns, and since it now carries the pooled estimate,
  `select_resolution()` with no criterion reads it by the size: it returns
  the level nearest `round(n / 4.5)`, on the automatic ladder its top level,
  with `spatialkit_size_not_reached`, where 2.0.0 returned Cp's pick.
  `criterion = "cp"`, or `estimate = NULL` in the profile, gives Cp's count,
  with the noise variance of the entry on Cp below.
  On the default ladder no variogram is fitted (`sac = "auto"`): `cp` and
  `reliability` are `NA`, `attr(, "cp_noise")` says why, and nothing is said
  about a variogram; on any other ladder it is estimated as before.  The
  default ladder's fits are bounded: from the 3 starts of the density kind
  the call stops above about 122,000 rows with an error of class
  `spatialkit_work_limit` that names the ways on (the limit is set where the
  fit takes about a minute, so a call just below it runs that long and
  prints nothing while it runs), and without the Matrix package `moran_i`
  and `moran_z` are `NA` above 22,500 rows, with `spatialkit_moran_skipped`.
  A profile passed straight as a count (`get_voronoi_seeds(n = prof)`,
  `voronoi_seeds_kmeans(k = prof)`, `build_tessellation(approx_n_cells =
  prof)`) is read by the size and no longer by `cp`.  `estimate = NULL,
  obs_per_cell = NULL` gives the previous ladder and columns, except `cp`
  wherever the entry on Cp's noise variance below applies; `sac = NULL`
  estimates the variogram on a ladder of sizes too.

* `resolution_profile()`, `get_voronoi_seeds()` and `voronoi_seeds_kmeans()`
  build density-weighted regions by default (`partition = "density"`), kept
  from the best of 3 k-means++ starts (`nstart = NULL`, `kmeans_nstart =
  NULL`).  Their cells are small where observations are dense, and hold
  about the same number of observations each once there are more cells than
  dense places.  Every default call of the three functions returns other
  centres and other cells.  The profile, its selections and its seeds carry
  the partition record (`attr(, "partition")`, `sel$partition`), and the
  prints name the kind.  In the profile `wss` stays the unweighted sum,
  `wss_spread` is the spread of the weighted one, and a start runs for at
  most 100 steps of Lloyd's iteration.  The seeding functions put the rows
  in coordinate order first, so their seeds no longer depend on the order
  of the rows, and `get_voronoi_seeds()`'s `attr(, "kmeans")` gains
  `converged` and `objective`.  `density_power` and `density_k` are read,
  and a default call can raise `spatialkit_empty_cells` and
  `spatialkit_not_converged`.  A cloud `get_voronoi_seeds()` draws inside a
  `boundary` is uniform, so the kind has no density to follow there.
  `method = "random"` and `"provided"` are unchanged, and a selection or a
  profile passed as the count returns the centres it scored, as before,
  except off a profile without an estimate fitted to a subsample, whose
  chosen level is refitted to every point (see the entry below on the
  refit of the chosen level).  `partition = "kmeans"` gives the previous
  regions with their previous starts: the profile's best of 25 k-means++
  starts, the seeding functions' `stats::kmeans()` from 10 random starts on
  the rows as they come.  A call
  that named `partition = "density"` and no number of starts runs 3 starts
  where it ran 25 (the profile) or 10 (the seeding functions); `nstart =
  25`, `kmeans_nstart = 10` or `nstart = 10` restores that.

* `get_voronoi_seeds()` given a number of starts that is not whole, such as
  `kmeans_nstart = 2.5`, and no `partition` now stops: the density kind
  refuses it, with "`kmeans_nstart` must be NULL or a single positive whole
  number of starts; got 2.5."  With `partition = "kmeans"` the call returns
  as before.

* `select_resolution()` reads a profile by `criterion = "auto"`: the size
  for pooled values, the leave-one-out minimum for plain means, `cp` for a
  profile built with `estimate = NULL`.  On a profile with an estimate the
  count is the level nearest `round(n / 4.5)`, not Cp's minimum: on the
  default profile that count itself, and on a ladder of counts the nearest
  level, with `spatialkit_size_not_reached` when it is not the count.  A
  response with no structure is still reported, by another test: `structured
  = FALSE`, and `spatialkit_no_structure` where the count becomes cells,
  from the leave-one-out error against one mean.  The count is then the
  size's where Cp's was 2 cells or the like, since a size does not shrink
  for want of structure, and `spatialkit_bound_chosen` is not raised for a
  size.  On a profile without an estimate the criteria are read as before,
  and `criterion = "cp"` reads `cp` on a profile that has it; what is new
  there is the `cp` column itself where the entry on Cp's noise variance
  below applies, the record a cp selection carries, and the seeds of a
  profile fitted to a subsample (see the entry below on the refit of the
  chosen level).  The default profile, one level with no variogram, has
  none, and `criterion = "cp"` or `"reliability"` on it is an error that
  says what to pass.  For Cp's count as 2.0.0 read it, build the profile
  with `estimate = NULL, obs_per_cell = NULL` (and `partition = "kmeans"`
  for the same regions): Cp's noise variance is then the variogram's
  nugget, as in 2.0.0, wherever that nugget is above the fit's lower bound
  and no more than the semivariance of nearest-neighbour pairs.

* `summarize_by_cell()` adds the pooled columns by default (`estimate =
  "pooled"`): `..pooled_*`, `..pooled_w_*` and `..pooled_target_*` for
  every numeric response and predictor, and the attributes `"pooled"` and
  `"intervals"`.  Every existing column is unchanged.  The default call can
  now raise the warnings of the pooled values, each once per call:
  `spatialkit_pool_unavailable` (a layer with no geometry, which the call
  accepts; no row with a value, a cell and a location; no two rows in any
  cell; values beyond what the arithmetic of pooled values holds),
  `spatialkit_pool_degenerate` (one cell, no spread inside any cell,
  or no signal between cells), `spatialkit_pool_few_cells` (2 to 9 cells),
  `spatialkit_pool_few_df`, `spatialkit_no_structure` (a response the cells
  do not help with), `spatialkit_pool_rows_skipped`, and, on cells that
  carry a count record, `spatialkit_estimate_mismatch` and
  `spatialkit_cells_differ`.  Small layers and layers without spatial
  structure raise several.  `estimate = "mean"` returns the previous frame.

* `resolution_profile()` no longer stops when a ladder set with `levels`
  has a level of more than 5,000 cells and the Matrix package is not
  installed: `moran_i` and `moran_z` are `NA` at such a level, with a
  warning of class `spatialkit_moran_skipped`, and every other column is
  computed.  Installing Matrix gives the two columns.

* The FNN package moves from Suggests to Imports and installs with
  spatialkit.  The nearest-neighbour searches use it, among them the one
  behind the weights of the density kind, now the default partition.  The
  package's own searches, which compare every pair of locations and give
  the same results in a time that grows with the square of their number,
  remain for a session in which FNN cannot be loaded.

* `estimate_sac_range()` fits its variogram model without iteration, so the
  same layer gives the same verdict, range and nugget whatever the unit of
  the response or of the coordinates and whatever the order of its rows.
  The criterion is gstat's (squared error weighted by N/h^2); for each trial
  range the nugget and partial sill are solved exactly, and the range is
  searched on a grid and refined.  Every range and nugget moves a little.
  Only the exponential family is fitted: the spherical fallback is gone.
  `rejected_reason` no longer takes `"variogram model did not converge"`;
  `attr(attr(x, "variogram_model"), "converged")` is `FALSE` when the
  criterion has no minimum inside the ranges searched.  Nothing is printed
  to the console.

* `estimate_sac_range()` tests for spatial structure before it returns a
  range: Moran's I of the modelled variable on its 2, 4, 8 and 16 nearest
  neighbours, with the smallest one-sided p-value multiplied by four.  When
  that is not below the new argument `structure_alpha` (default 0.05) the
  result is `NA` with `rejected_reason` `"no spatial structure detected at
  the nearest-neighbour scale"` and a pure-nugget model.  On white noise
  (300 points) a finite range came back for 24 of 60 layers and now comes
  back for 1.  The test has little power where few neighbours lie inside the
  range.  `structure_alpha = 1` switches the refusal off.  A variogram that
  no model with a range fits better than a constant is refused as `"the
  empirical variogram is flat: no spatial structure at the lags resolved"`.

* `estimate_sac_range()` gains `lags = c("adaptive", "equal")` and `n_lags`,
  and `"adaptive"` is the default: the first bin ends where the closest
  `max(50, n / 4)` pairs lie, and the rest widen geometrically to the
  cutoff.  Fifteen equal bins hid structure shorter than the first bin: on a
  clustered layer with a true range of 30 the median estimate was 77 and is
  now 29.  A longer second structure seen at equal lags is kept in the model
  and sets the range (`structures`, `short_range`), and a range the adaptive
  lags leave past the largest lag falls back to the equal-lag result
  (`lag_fallback`).  Ranges, nuggets and the block sizes of
  `make_folds(auto_range = TRUE)` change.  `lags = "equal"` keeps the old
  bins.

* `estimate_sac_range()` never returns the largest directional range in
  place of an all-pairs fit that identifies no range: the result is then
  `NA`.  The directional ranges are still attached, and `anisotropy_used`
  is always `FALSE`.

* `resolution_profile(range_floor = FALSE)` is the default.  The ladder
  starts at 2 cells; the floor the autocorrelation range implies is reported
  in the bounds and the print, and is a level of the ladder when the layer
  supports it.  A range is an estimate and can be spurious, and a ladder
  that started at it could not show the counts below.
  `select_resolution(criterion = "reliability")` reads the levels at or
  above the floor, so it returns the floor, as it did when the ladder
  started there; its new argument
  `range_floor = FALSE` reads every level.  `resolution_profile(range_floor
  = TRUE)` starts the ladder at the floor.

* `resolution_profile()` takes its support ceiling from the rows that carry
  a response (and every predictor, when predictors are given), not from
  every point.  A layer with a response on a quarter of its rows gets a
  quarter of the ceiling, `ceiling_from` reads `"min_cell_n (rows with a
  response)"`, and a warning of class `spatialkit_response_ceiling` gives
  both counts.  New columns `cells_resp` and `cell_resp_mean` count the
  cells that hold a scored row and the scored rows per such cell.

* `moran_z` is a diagnostic of the cell means and is no longer read as a
  count rule.  It grows as cells get finer than the range and is near zero
  at the coarsest counts it is computed at, whatever the data.
  `select_resolution(criterion = "moran_z")` is an error of class
  `spatialkit_diagnostic_criterion` unless `diagnostic = TRUE`; `summary()`
  leaves it out unless it is named in `criteria`; `plot()` draws its curve
  with no selected level.  `determine_optimal_levels()` still ranks by it
  under `criterion = "morans_i"`, and raises a warning of class
  `spatialkit_diagnostic_ranking` when it does; under `"combined"` it is
  reported and no longer ranks (below).

* `get_voronoi_seeds()`, `voronoi_seeds_kmeans()` and `build_tessellation()`
  read a `resolution_profile()` built without an estimate by `cp` or not at
  all (one built with an estimate, the default, is read by
  `criterion = "auto"`).  A profile whose `cp` is `NA` was read by the next
  finite criterion with a log line; it is now an error (class
  `spatialkit_no_default_criterion` when another criterion is finite).  A
  geometry-only profile whose only finite criterion is `elbow` is still read
  by it.  Pass `select_resolution(profile, criterion = )` to choose another
  criterion.

* `select_features_forward()` accepts a candidate by a new default rule,
  `rule = "se"`: its pooled gain has to exceed `tol` and the mean of its
  per-fold gains has to exceed `se_mult` (default 1) standard errors of that
  mean.  With `tol = 0` any gain was accepted.  `rule = "tol"` is the old
  test.  Every learner now has a baseline, so the selection can be empty
  (see Bug fixes).

* `area_of_applicability()` weights the predictors by the model's variable
  importance by default (`weights = "auto"`): `pmax(importance, 0)` when the
  model carries a finite importance for every predictor and has out-of-bag
  skill, equal weights otherwise, with the reason logged.
  `params$weights_source` says which.  `weights = NULL` gives equal weights
  as before.

* `area_of_applicability()` measures location for a `gwr_fit`, a
  `bayesian_fit` and a forest fitted with `include_coords = TRUE`, which
  predict as a function of location.  Location has an index of its own
  (`DI_location`, element `location`): the distance to the nearest training
  point over the mean distance between training points, with the same folds
  and threshold rule.  A point is inside only when it is inside on both
  indices; `n_outside_predictors` counts the predictor index alone, and a
  `threshold` the caller passes sets that index alone.  The forest keeps its
  coordinates in the predictor index as `..x` and `..y`, so for it
  `DI <= threshold` is the old answer under the same weights; its
  `params$coords` reads `"location index and predictors"`, and `print()`
  names both.  Three more things change for the forest.  A row outside on
  location is `FALSE` also when a predictor is missing, where it was `NA`,
  so `n_na` can fall and `n_outside` rise.  Folds that split repeat visits
  to a site can give it a location threshold of 0, with the
  `spatialkit_aoa_zero_threshold` warning, and then only points at a
  training location are inside.  A `train_sf` passed in a geographic CRS
  (longitude/latitude), which the predictor index measured in degrees, is
  refused: project it, or pass `coords = "no"`.  The default `train_sf` is
  the forest's own data, which `fit_rf_model()` projects.  For a GWR or a
  Bayesian fit `coords = "no"` gives the old result; for the forest it
  leaves the coordinates out of both indices.  `coords = "yes"` measures
  location for any model.

* `fit_rf_model()` splits an unordered factor or character predictor with
  `respect.unordered.factors = "order"`, which ranger's help recommends;
  ranger's default, `"ignore"`, splits on the order the levels are stored
  in.  Forests with such a predictor change.  Pass
  `respect.unordered.factors` through `...` to set it yourself.
  `$info$categorical_predictors` and `$info$respect_unordered_factors` record
  what was done, and `print()` shows it.  `print()` also reads "no spatial
  information used" for a forest fitted without the coordinates.

* The local regressions of a GWR are solved by the package, by pivoted QR on
  predictors centred in each window, and GWmodel is no longer needed to fit,
  predict, cross-validate or select (`sp` still is).  A window of deficient
  rank neither stops the fit nor returns coefficients made of rounding
  error: a coefficient it cannot estimate is `NA`, and a warning of class
  `spatialkit_gwr_rank_deficient` counts them.  On windows of full rank the
  results equal `GWmodel::gwr.basic()`'s to rounding error.  `fit$engine` is
  a list in GWmodel's layout (`SDF`, `GW.diagnostic`, `GW.arguments`, `lm`)
  and no longer an object of class `gwrm`.  `print()` heads the fit `<GWR>`,
  not `<GWR (GWmodel)>`.

* With `bandwidth = NULL`, `fit_gwr_model()`, `gwr_model_selection()` and
  `cv_gwr()` select the bandwidth by a search of their own: 16 log-spaced
  bandwidths are scored and the best is refined.  `GWmodel::bw.gwr()`'s
  golden-section search assumes one minimum and returned a near-global
  bandwidth when AICc had two.  The scored curve is kept in
  `info$bandwidth_search`.  Selected bandwidths can differ from 2.0.0's;
  supply `bandwidth` to fix one.  Below 20 observations the bandwidth is
  searched like any other.

* `fit_gwr_model()`, `gwr_model_selection()` and `cv_gwr()` select the
  default bandwidth by the leave-one-location-out sum of squares, not AICc,
  when rows share a location (repeat visits, replicate samples).  AICc
  counts such rows as independent and is lowest where the windows reproduce
  each location's mean.  A warning of class
  `spatialkit_gwr_repeated_locations` gives the counts;
  `info$n_locations` and `info$bandwidth_search$approach` record it.  A
  layer with no repeated location is searched by AICc.

* `fit_bayesian_spatial_model(standardize_predictors = TRUE)` scales its
  automatic slope prior to the response: `normal(0, 2.5 * sd(y))` under an
  identity link and `normal(0, 2.5 * sd(log(y)))` for the lognormal
  families.  The fixed `normal(0, 5)` shrank a slope of about 14,000 per SD
  on a response in dollars to 0.06.  Other links keep `normal(0, 5)`.  The
  prior used is stored as `$info$b_prior`; pass `prior` to set your own.

* `$info$convergence_ok` of a Bayesian fit follows Stan's own diagnostics
  more closely.  It is `FALSE` on any divergent transition, an R-hat above
  1.05, an R-hat above 1.01 with an effective sample size under 100 per
  chain, a bulk or tail effective sample size under 100 per chain, or an
  E-BFMI under 0.3 in a chain.  Short runs that passed can now fail.  A
  failed verdict raises one warning of class `spatialkit_convergence` that
  names each failed check and its remedy.

* `fit_bayesian_spatial_model()` samples at `adapt_delta = 0.99` by default,
  where it sampled at 0.9, and so do `cv_bayes()`, `compare_models_cv()`
  and every other call that leaves `control` unset.  At 0.9 ordinary data
  left a few divergent transitions, and one is enough to set
  `convergence_ok` to `FALSE` and raise `spatialkit_convergence`: four
  layers of 250 points fitted at the default length from three seeds each
  had 1 to 24 in every fit, and none at 0.99.  Sampling takes about two to
  three times as long, and a fit at a given seed draws other values, so its
  estimates, intervals and predictions change.
  `control = list(adapt_delta = 0.9)` restores the previous setting.  The
  warning's advice prints `adapt_delta` as it was set; it printed 0.975 as
  0.97, 0.995 as 0.99 and 0.999 as 1.00.

* `summarize_by_cell(deff = "variogram")` falls back to `"kish"`, not to
  the uncorrected standard error, when no usable variogram exists.  The
  warning of class `spatialkit_deff_fallback` gives the reason.  A new
  column `deff_method`, present whenever a design effect is requested, names
  the correction each row carries.  A variogram refused because the
  response shows no spatial structure is not a fallback: its pure-nugget
  model gives a design effect of 1.  It gives that design effect to every
  column, and a warning of class `spatialkit_deff_no_structure` names the
  predictor columns summarised under it, since the verdict is about the
  variable the variogram was fitted to; `deff = "kish"` corrects each
  predictor by its own ICC.

* `summarize_by_cell(deff = "kish")` estimates the ICC of each predictor on
  its own.  One ICC was pooled over every predictor listed, so listing a
  second predictor changed the first one's standard error.
  `attr(, "icc")$pred` and `attr(, "deff_applied")$icc_pred` are vectors
  named by predictor column.

* `summarize_by_cell(deff = "kish")`, and a `"variogram"` request that falls
  back to it, return `NA` standard errors, in every cell of two or more
  values, for a column whose ICC cannot be estimated: its values in one
  cell, fewer than two cells holding two or more of them, or values that
  fall in no cell or differ only there (rows with no cell id give no cell to
  the ICC).  They returned the uncorrected ones, as for an ICC estimated
  at or below 0, under the same label.  `..neff_*`, `..ci_*` and, for the
  column that sets it, `cell_weight` are `NA` in those cells too, so a
  regression weighted by `cell_weight` leaves them out; a cell of one
  value keeps its `cell_weight` of 1.  A warning of class
  `spatialkit_icc_unavailable` names the columns and the reason,
  `attr(, "icc")` holds `NA` for each, with `reason`, and under
  `estimate = "pooled"` the `method` of `attr(, "intervals")` names them.
  An estimate at or below 0 still gives an ICC of 0 and the uncorrected
  standard errors, and a column whose values are all equal, those of rows
  with no cell id included, keeps standard errors of 0 (its ICC is `NA`,
  with the reason "the column does not vary").  `deff = 1` gives the
  uncorrected standard errors in every case.

* Voronoi cells built on seeds reach the edge of the data the seeds stand
  for.  `build_tessellation(seeds, method = "voronoi")` without a boundary
  clipped the cells to the seeds' own hull, and the observations outside it
  were dropped by the assignment: a quarter of 3,000 uniform points at 25
  k-means seeds.  The seeding functions now attach the boundary they drew
  in, or the hull of the points they clustered, as `attr(, "cover")`, and
  the cells are built over it.  Pass `boundary` to set the extent yourself.

* `assign_features_to_polygons()` stops on a duplicated or missing value in
  the polygon ID column.  Nine polygons with a zone attribute called `id`
  came back as three.  When the rows that share an ID are parts of one unit,
  merge them first with `dplyr::summarise(dplyr::group_by(polygons_sf, <ID
  column>))`.

* Conditions that were log lines are R warnings with a class, so
  `tryCatch()`, `withCallingHandlers()` and `options(warn = 2)` see them:
  rows dropped before a fit (`spatialkit_rows_dropped`), a constant
  predictor (`spatialkit_constant_predictor`), features in no polygon
  (`spatialkit_unassigned`), fewer folds than asked for
  (`spatialkit_k_reduced`), a cell count chosen by a bound of the ladder
  (`spatialkit_bound_chosen`) and a profile the layer does not support
  (`spatialkit_not_supported`).  `?resolution_profile` lists the classes of
  the resolution functions.

* `determine_optimal_levels()` returns `NA` where its WSS curve has no
  elbow (`on_no_elbow = "na"`, the default).  The count it returned there
  is where the chord rule lands on a curve with no bend, about the square
  root of `max_levels` whatever the points, so the ladder chose it and not
  the data.  A ladder too short to show an elbow gives `NA` too, and so
  does a layer of one or two points or of one location, where it gave 1.
  The `NA` comes with a warning of class `spatialkit_no_elbow` that says
  why and carries the count as `k`, and with `attr(, "structured")`,
  `FALSE`, and `attr(, "sag")`, the largest log-log sag of the curve.
  Every result at the default carries the two attributes, `structured`
  being `TRUE` when the WSS curve has an elbow, and the counts themselves
  are unchanged.  Under the model-aware criteria a curve with no elbow
  leaves no window to score: the call returns `NA` without ranking counts
  by Moran's z, its `NA` carries the `"diagnostics"` of the sweep (every
  Moran's I `NA`, `eval_ks` empty), and its log lines no longer say the
  window is scored or the call falls back to the geometric ranking.
  `build_tessellation()` and `get_voronoi_seeds()` stop on the `NA` and say
  where it came from.  `on_no_elbow = "warn"` returns the count, with the
  same warning and without the `structured` and `sag` attributes.

* `determine_optimal_levels()` keeps the criterion it is given.  Supplying
  both `response_var` and `predictor_vars` switched `criterion =
  "geometric"`, the default, to `"combined"`, also when `"geometric"` was
  named, so the count read the response unasked.  Such a call now returns
  what the call without the two variables returns, with no `"diagnostics"`
  attribute and no `spatialkit_diagnostic_ranking` warning, and a line
  logged at the INFO level says that `criterion = "combined"` reports
  Moran's I of the cell means.  The columns must still exist; a factor or
  character response is refused only under a criterion that reads it.
  Name `criterion = "combined"` for the path the switch took, whose ranking
  the next entry changes.

* `determine_optimal_levels(criterion = "combined")` ranks the counts
  around the elbow by the WSS curve's log-log sag alone.  Moran's z of the
  cell means is dropped from the ranking, not replaced by another
  statistic: averaged with the sag, it put first ten, the lowest count it
  can be computed at, or the count after it.  On a curve with no elbow the
  sag is flat and the counts below ten, which Moran's z cannot score, took
  the last place on its axis; with an elbow its |z| is near 0 at the
  coarsest counts it scores whatever the response.  On 600 uniform points
  with a smooth field, under `on_no_elbow = "warn"`, the call returned
  `10 6 5` where `"geometric"` returns `6 5 7`, and it now returns `6 5 7`;
  on eleven clusters with a smooth field it returned `10 11 12` and now
  returns `11 12 10`, the elbow, 11, first.  The sag is that of the sweep's
  WSS curve, the one the elbow is read from, so with an elbow the elbow
  comes first and the counts after it follow by sag; they can differ from
  the two neighbours `"geometric"` returns.  With no elbow the order is
  `"geometric"`'s.  The k-means re-run behind Moran's I no longer ranks:
  where it landed in a worse optimum at the elbow than the sweep did,
  another count went first, and on eighteen loose clusters the elbow, 15,
  was not among the counts returned.  Its WSS stays in `attr(,
  "diagnostics")` as `wss_eval`, beside `moran_i` and `moran_z`, and
  `combined_rank` is the rank of the sweep's sag.  Where Moran's I can be
  computed at no count of a window that reaches ten cells (a constant
  response, too few cells holding a row with one, or too few in the half
  `select_on = "split"` reads), `"combined"` returned the geometric ranking
  without diagnostics; it now ranks as it does otherwise, with `moran_i`
  and `moran_z` `NA` and a logged warning, so neither the response nor the
  split changes the counts it returns.  `"morans_i"` still falls back
  there.  `"combined"` no longer raises `spatialkit_diagnostic_ranking`.
  `criterion = "morans_i"`, which still ranks by Moran's z, does; the
  condition no longer carries `upgraded`, and its message names
  `resolution_profile()` and `select_resolution()`, which run on the
  default profile, where it named `select_resolution(criterion = "cp")`,
  an error there.  No argument gives the earlier ranking.

* Cp's noise variance in `resolution_profile()` is no longer a nugget at
  the variogram fit's lower bound (0, or at most 1e-4 of the sill), and it
  is held to the semivariance of nearest-neighbour pairs the variogram
  carries (`attr(sac, "short_lag_semivariance")`); Bug fixes has the
  details.  Wherever either applies, the `cp` column, Cp at one cell
  (`attr(, "one_cell")`) and Cp's count differ from 2.0.0's on every ladder
  that fits a variogram, `estimate = NULL, obs_per_cell = NULL` included,
  and no argument gives the earlier noise variance back.
  `attr(, "cp_noise")` names the source used and, under the hold, the one
  it replaced.  A selection by `cp` off a profile without an estimate gains
  `count_rule`, which its seeds, the cells built on them and the assigned
  layer carry on.

* `select_resolution()` refits the level it chooses off a profile without
  an estimate fitted to a subsample (a ladder of counts on more than
  `sample_n` points) to every point before it hands the seeds on (Bug fixes
  has the details).  So those seeds, and the seeds `get_voronoi_seeds()`
  and `voronoi_seeds_kmeans()` return for such a selection or profile, are
  not the centres 2.0.0 returned, `partition = "kmeans"` and `estimate =
  NULL, obs_per_cell = NULL` included.  The centres the profile scored,
  which 2.0.0 returned, are
  `attr(profile, "centres")[[as.character(sel$best)]]`, in the CRS
  `attr(profile, "centres_crs")`.

* `resolution_profile()` raises one warning of class
  `spatialkit_size_not_reached` for all the sizes of a ladder held for the
  same reason, where it raised one for each size.  Its `obs_per_cell` and
  `wanted` hold one element for each such size, in the order of the ladder,
  so a script that read them as single numbers gets vectors where several
  sizes are held; `used` and `reason` are single values as before.  One
  size held gives the condition it gave.

* `estimate_sac_range()` gives the `rejected_reason` "no variogram model
  could be fitted: fewer than three lag bins hold point pairs" where it gave
  "no variogram model could be fitted (singular fits)": nothing was fitted
  in that case.  A script that compared the reason with the old text can
  compare its start, "no variogram model could be fitted", which both share.

* `estimate_sac_range()` needs 30 distinct locations as well as 30 rows
  (Bug fixes has the details).  A layer of 30 or more rows at 2 to 29
  distinct locations, to which 2.0.0 fitted a variogram, gets the unclassed
  `NA` whose `rejected_reason` begins with the count of its distinct
  locations.  So on such layers a range that came back is `NA`,
  `summarize_by_cell(deff = "variogram")` returns Kish standard errors with
  `spatialkit_deff_fallback`, and `make_folds(auto_range = TRUE)` builds
  geometric blocks.  No argument gives the earlier result back: a range
  read off fewer than 30 locations is not identified.

* `new_spatial_fit()` stops on a fit its methods cannot use, with an error
  that names the argument and what it got: `data_sf` must be an sf object,
  `formula` a formula, `response_var` the name of one attribute column of
  `data_sf` (a column other than its geometry), `predictor_vars` a
  character vector of its attribute columns (`character(0)` for none),
  `subclass` one class name and `info` a list.  It also stops when
  `predictor_vars` names the response, when the left side of `formula`
  does not hold `response_var` (a transform such as `log(price) ~ elev`
  holds it), and when an object classed sf has no geometry column (its
  `"sf_column"` attribute missing, not one name, or naming a column that is
  not an sfc).
  It built a fit from any of them: the response among the predictors was
  measured in `area_of_applicability()`'s distances and had to be present
  in the new rows, a formula such as `elev ~ price` with
  `response_var = "price"` printed one model while `summary()` and
  `model_metrics()` scored the fitted values against `price`, the object
  classed sf without geometry made `summary()`, `model_metrics()` and
  `area_of_applicability()` stop with "missing value where TRUE/FALSE
  needed", a column number or a factor as `response_var` made `summary()`
  and `model_metrics()` score the fitted values against the column at that
  position, a misspelt name or the geometry column's name failed later in
  another function with a message about something else, and a `data_sf`
  without geometry made `residual_morans_i()` blame a missing `residuals()`
  method.  A `fit_fn` that passed `st_drop_geometry(train_sf)` has to pass
  `train_sf`.  One whose fit names a column the training rows lack (a
  predictor derived inside `fit_fn`, say) fails every fold of
  `cv_spatial()`, with the message in `fold_status`: add the column to
  `train_sf` before building the fit, or leave it out of `predictor_vars`,
  as the README's `trend_fit` does with the coordinates it derives.  A
  `brms::bf()` formula is passed as its `$formula`.  In
  `select_features_forward()`, a `fit_fn` whose fit for no predictors is
  refused (`predictor_vars = NULL` where `character(0)` is meant, say) can
  no longer give the null-model baseline: the baseline is the training-fold
  mean, as for any learner that cannot fit the empty set, and
  `baseline$type` says `"fold mean"`.

* `new_spatial_fit()` stops when the right side of `formula` and
  `predictor_vars` disagree: when `predictor_vars` names a column the right
  side does not use, or the right side uses an attribute column of
  `data_sf`, other than the response, that `predictor_vars` does not name.
  Both are read by `all.vars()`, so transforms (`log(elev)`), interactions,
  `.`, names that are not columns of `data_sf` and the reserved coordinate
  names `..x` and `..y` build; the columns of an `offset()` or a grouping
  term count as used.  A `fit_fn` that derives coordinates `x` and `y` for
  a trend surface, on a layer that keeps its coordinates as columns `x` and
  `y`, is refused: the error says to name such columns in `predictor_vars`
  or to derive them under names that are not columns of `data_sf`.  It
  built such a fit, whose `print()` showed
  one model while `area_of_applicability()` measured distances in other
  predictors and `residual_morans_i()` rebuilt the OLS design from them: a
  fit of `price ~ elev` with `predictor_vars = "slope"` put rows 10 standard
  deviations outside the training range of `elev` inside the area of
  applicability, and its Moran test switched to the randomisation null in
  silence.  A `fit_fn` whose formula and `predictor_vars` disagree fails
  every fold of `cv_spatial()`, with the message in `fold_status`.

## New features

* `resolution_profile()`, `get_voronoi_seeds(method = "kmeans")` and
  `voronoi_seeds_kmeans()` gain `partition = c("density", "kmeans")`.
  `"kmeans"` is what each function did: centres that minimise the sum of
  squared distances from the points.  `"density"`, the default (see
  Breaking changes), weights each point by its local density (the points in
  the disc that reaches its `density_k`-th nearest neighbour, default 8,
  over the disc's area) to the power `density_power` (0 to 2, default 1).
  The cells are compact, and at power 1 they hold about the same number of
  points each once there are more
  cells than dense places; with fewer cells the counts are about as uneven
  as under `"kmeans"`.  On 2,000 simulated points with 28 dense places the
  coefficient of variation of the cell counts was 0.90 under `"kmeans"` and
  0.97 under `"density"` at 10 cells, and 0.95 and 0.31 at 60 cells.  The
  result does not depend on the order of the rows, the unit of the
  coordinates or whether FNN is installed, and points at one location share
  a cell.

* The same three functions gain `init`: `"kmeans++"` draws a start's centres
  one at a time, each in proportion to its squared distance from those
  already drawn, and `"random"` draws distinct locations at random, as
  `stats::kmeans()` does.  `resolution_profile()` defaults to `"kmeans++"`,
  as before.  The two seeding functions default to `NULL`: random starts for
  `partition = "kmeans"`, the run they have always made, and k-means++ for
  `partition = "density"`.

* `resolution_profile()` and `voronoi_seeds_kmeans()` gain `iter_max`, the
  most iterations a start runs for.  `NULL`, the default, is 50 for the
  k-means each has always run and 100 for random starts in the profile and
  for `partition = "density"`.  A kept start that was cut off at the limit,
  or at `kmeans_iter` in `get_voronoi_seeds()`, is a warning of class
  `spatialkit_not_converged` that carries `level`, `iter` and `iter_max`.
  The k-means kind from its own starts reports as before.

* A result built with `partition = "density"`, or with an `init` other
  than the k-means kind's own, says so: `attr(, "partition")` on a profile
  and on seeds, `$partition` on a `select_resolution()` result,
  `params$partition` on a tessellation from
  `build_tessellation()` or `create_voronoi_polygons()`.  The record holds
  `kind`, `init`, `density_power` and `density_k`, and for a profile and for
  seeds a seeding function placed, the iterations the kept start took and
  whether it converged.  A selection or profile passed as the count gives
  the centres it scored, so the regions built are the regions scored.
  Results of the k-means kind from its own starts carry no record and are
  unchanged.

* Under `partition = "density"` the profile's `wss` is the within-cell sum
  of squares with every point counted once, so `wss` and `elbow` compare
  between kinds; the weighted sum the centres minimise is
  `attr(, "partition")$objective`.  `print()` of the profile and of a
  selection, and `summary()` of the profile, name the kind.
  `determine_optimal_levels()` keeps plain k-means.

* Two more warnings have a class.  `spatialkit_partition_unused`: a
  selection or profile passed as the count together with a `partition`,
  `init` or density setting it was not built with; the scored centres are
  returned.  `spatialkit_empty_cells`: a density partition in which a centre
  holds no point, which happens at counts from about a fifth of the points
  up, mostly from random starts.

* `summarize_by_cell()` gains `estimate = "pooled"`.  Beside each mean it
  returns a pooled value, `..pooled_*`: the cell's mean pulled toward the
  count-weighted mean of the 8 cells whose observations' centroids are
  nearest, `..pooled_target_*`, with the weight `..pooled_w_*` on its own
  mean.  The weight grows with the cell's count and with how much the cells
  differ from their neighbours.  A cell with one observation gets a value; a
  cell with none gets `NA`.  It works on any polygons.  Every present column
  is unchanged, and `estimate = "mean"` returns the frame it always did;
  `"pooled"` is the default (see Breaking changes).  `attr(, "pooled")`
  records the fit, and whether the cells predict a left-out observation
  better than one mean for all rows.  On 200
  simulated layers of 1,500 clustered observations in cells of 4.5 on
  average, the root mean squared error of the pooled value at held-out
  observations among the others was 0.96 of the plain mean's
  (`dev/sim-pooled-values.R`; `?summarize_by_cell` describes the
  simulation).

* `summarize_by_cell()` gains `obs_level`: `..obs_lo_*` and `..obs_hi_*`, a
  range for one further observation in the cell, as the pooled value plus
  two quantiles of its leave-one-out residuals.  It is a statement about
  observations on average and describes observations among those
  summarised; it is too narrow for a place with none nearby.  `conf_level`
  and the `..ci_*` columns are unchanged; that interval is not a range for
  an observation.  On the same layers the range at 0.90 held 89.6 percent of
  the observations held out among the others and 75.1 percent of those in a
  held-out block, and `..ci_*` at 0.90 read that way held 58.9 percent of
  the first.  `attr(, "intervals")` says what each interval in the frame
  covers.

* `resolution_profile()` gains `estimate` and `obs_per_cell`.
  `obs_per_cell` sets the ladder in observations per cell and accepts
  fractions; the centres are then fitted to every row.  With `estimate` the
  profile adds `loo`, the leave-one-out error of that cell value on every
  row of the layer with a response, in the cells the seeds will build, with
  `loo_se`, `scored_per_cell`, `share_alone` and `own_weight`.  It needs no
  variogram.  `print()`, `summary()` and `plot()` of such a profile show
  it.

* `resolution_profile()`'s `sac` accepts `"auto"`: the variogram is
  estimated as before unless `obs_per_cell` sets the ladder, where none is
  fitted, `cp` and `reliability` are `NA`, and `attr(, "cp_noise")` says
  why.  `"auto"` is the default (see Breaking changes).

* `resolution_profile()`'s `nstart`, `get_voronoi_seeds()`'s `kmeans_nstart`
  and `voronoi_seeds_kmeans()`'s `nstart` accept `NULL`, the kind's own
  number of starts: the earlier default for `partition = "kmeans"` and 3
  for `partition = "density"`.  `NULL` is the default (see Breaking
  changes).  On a ladder set by `obs_per_cell` a fit whose cost grows with
  the square of the rows is stopped above a limit, with an error of class
  `spatialkit_work_limit` that names the ways on; so is the search behind
  the weights of `partition = "density"` when the FNN package is not
  installed.

* `select_resolution()` gains `criterion = "size"` (the level at
  `obs_per_cell` observations per cell, 4.5 unless given), `"loo"` (the
  lowest leave-one-out error) and `"auto"` (size for pooled values,
  leave-one-out for plain means).  A size does not test the response, so
  `structured` is computed for it from the leave-one-out error against one
  mean for the layer.  A profile built with `estimate`, passed straight to a
  seeding function as the count, is read by `"auto"`, which is the default
  (see Breaking changes).

* A selection made by these rules carries `count_rule`, and the seeds, the
  tessellation, the assigned layer and the summary carry it on.
  `summarize_by_cell()` warns (`spatialkit_estimate_mismatch`) when the
  count was chosen for one cell value and the summary reports the other, and
  (`spatialkit_cells_differ`) when the cells or rows it is given are not the
  ones the count was scored on.

* The record holds the centres the profile scored (`scored$centres`).
  `build_tessellation(method = "voronoi")` and `create_voronoi_polygons()`
  build the cells of those centres in the CRS the profile worked in, with no
  `crs` or a geographic one, whatever CRS the seeds arrive in (degrees, or
  the CRS of a boundary they were placed in), each seed put at its centre as
  the record holds it.  On such cells `assign_features_to_polygons()` gives
  a point on the edge between two cells the cell `resolution_profile()`
  scored it in: its nearest centre in the record's CRS, and of two the same
  distance away (squared distances within one part in ten million of each
  other) the lower-numbered one.  On gridded locations (raster cell
  centres, a regular sampling grid) the rows each cell holds are then the
  rows its count was scored on, and the summary's check reads `match`.
  `cell_values_at()` places a location on such an edge in the same cell,
  without `seeds` and with the seeds in any order.  Polygons without the
  record, cells whose record says they are other cells built from the
  count, and features that are not points are placed by the polygons alone.

* New `cell_values_at()`: the value and interval of the cell at other
  locations, with each location's distance from the nearest observation.

* Each warning of the pooled path of `summarize_by_cell()` is raised once
  per call and reason and names the columns it applies to, in its text and in
  a `columns` field.  Where the pooled columns cannot be computed they are
  `NA` with a warning and the plain columns are returned; the function
  raises no error for them.

* New condition classes: `spatialkit_pool_unavailable`,
  `spatialkit_pool_degenerate`, `spatialkit_pool_few_cells`,
  `spatialkit_pool_few_df`, `spatialkit_pool_rows_skipped`,
  `spatialkit_repeated_locations`, `spatialkit_few_residuals`,
  `spatialkit_estimate_mismatch`, `spatialkit_cells_differ`,
  `spatialkit_size_not_reached`, `spatialkit_ladder_unused`,
  `spatialkit_moran_skipped`, and the error class `spatialkit_work_limit`.
  `spatialkit_no_structure` is also raised by `summarize_by_cell(estimate =
  "pooled")`, and by the seeding functions for a count set by a size.

* `make_folds()` gains `method = "knndm"`, the k-fold version of
  nearest-neighbour distance matching (Linnenbrink et al. 2024).  It builds
  `k` folds whose held-out distances are as close as a set of candidates
  allows to the distances from `prediction_points` to the data, and has no
  row limit, where `"nndm"` takes at most 5,000 rows and one fit per point.
  The candidates are the random partition and, for up to 50 cluster counts,
  the samples clustered on their coordinates and packed into `k` folds; the
  one with the smallest Wasserstein distance to the target is returned, and
  `params$candidates` lists them all.  The match is coarser than NNDM's.
  Two cases raise a classed warning: prediction points nearer to the data
  than any fold can hold a point out (`spatialkit_nndm_out_of_reach`), and
  folds still closer to their training data than the prediction points are
  to the samples (`spatialkit_knndm_short_of_target`).  The rule that packs
  clusters into folds is this package's and was not checked against
  `CAST::knndm()`.

* `make_folds()` gains `max_test` for `"buffered_loo"` and `"nndm"`: folds
  are built for a seeded random subset of that many held-out points, with
  every row still available for training.  `"buffered_loo"` then runs on
  layers its storage limit refuses.  It does not lift the 5,000-row limit of
  `"nndm"`; `"knndm"` is the design above it.

* `fold_separation()` gains `prediction_points`: the locations to be
  predicted.  The held-out distances are then set beside the distances from
  those points to the training data, and the result's `target` attribute
  carries a `reading`: `"closer"` (the cross-validated error tends to read
  low for that target), `"farther"` (it tends to read high), `"matched"` or
  `"mixed"`.  NNDM and kNNDM folds are read against the prediction distances
  they carry.  The table gains `target_median_dist` and `dist_gap`.

* Every classed `estimate_sac_range()` result says how far its nugget can be
  trusted.  `nugget_at_bound` is `TRUE` when the fit put the nugget at zero,
  which is the bound of the search and not an estimate; `first_lag` is the
  shortest lag resolved; `short_lag_semivariance` is the semivariance of
  nearest-neighbour pairs with no model, an upper bound on the nugget.
  `?sac_nugget` describes them.  The result also records `n_used`,
  `n_missing`, `n_subsampled`, `structure_p`, `morans_i`, `lags` and
  `lag_boundaries`.

* `estimate_sac_range()` reads rows that share a location as one place.  The
  test for spatial structure is made on the location means, pairs at zero
  separation are in no lag bin, and their semivariance, the variance between
  visits, is returned as `zero_lag_semivariance`.  When the rows outnumber
  the locations by 30 or more, the nugget fitted to the variogram is not let
  under it (`nugget_floor`).  With `detrend = "reml"` the fit takes one row
  a location, drawn under `seed`, and takes no such floor.  A layer needs
  30 distinct locations, not only 30 rows.

* `select_resolution(criterion = "cp")` says when the count is not
  identified.  `structured` is `FALSE` when one cell for the whole layer
  scores no more than `tol` above the optimum (below it included), or when
  the cell means at the optimum differ from one mean by no more than chance
  (an F test at 5 percent, adjusted for the levels searched;
  `structure_p`).  The pick is unchanged.
  On white noise of 300 points it was `FALSE` on 145 of 150 layers.  The
  functions that take the count raise a warning of class
  `spatialkit_no_structure`.  `resolution_profile()` keeps what one cell
  scores in `attr(, "one_cell")`.

* `determine_optimal_levels()` gains `on_no_elbow = c("na", "warn",
  "error")`; its default is in Breaking changes.  `"na"` returns
  `NA_integer_` in place of the chord rule's answer and `"warn"` returns
  that answer, each with a warning of class `spatialkit_no_elbow` that
  carries the count and `structured = FALSE`; `"error"` stops with an error
  of that class.
  `criterion = "morans_i"` and `"combined"` run with a response and no
  predictors, where they fell back to the geometric ranking.  A layer of one
  or two points, or with every point on one location, is read by the same
  rule (`NA` at the default, 1 under `"warn"`); it returned 1 with nothing
  said.

* A GWR fit gives its standard errors.  `coef(fit, se = TRUE)` returns the
  local estimates with their standard errors and t values, and
  `predict(fit, newdata, se.fit = TRUE)` the standard error of the fitted
  value (`se.fit`) and of a new observation (`se.pred`).  They follow
  GWmodel's definitions, treat the bandwidth as known, and are not adjusted
  for the number of local tests.  `print(summary(fit))` shows the effective
  number of parameters and the AICc and RMSE of a global least-squares fit
  on the same rows (`info$global_ols`): a GWR whose AICc is no lower has not
  shown that its coefficients vary.

* `fit_gwr_model(predictor_vars = character(0))` fits the model with a local
  intercept alone, the kernel-weighted mean of the response.  It is the
  baseline a GWR with predictors has to beat, and
  `select_features_forward()` fits it as the null model of a GWR sweep.

* `predict()` on a GWR says when a location is beyond the reach of the fit.
  When a new location's nearest training point is farther than the widest
  kernel at any training location reaches, the prediction is returned and a
  warning of class `spatialkit_gwr_extrapolation` gives the count.  The
  warning that counts `NA` predictions has class
  `spatialkit_gwr_not_estimable`.  Both carry their counts as fields.
  `predict_surface()`, `cv_gwr()`, `cv_spatial()`, `cv_block_size_sweep()`
  and `select_features_forward()` raise each once per call, with totals.

* `gwr_model_selection()` scores the model without predictors as `<none>`,
  and `best` is `character(0)` when it wins.  The table gains `delta`, each
  model's criterion minus the best one's, and `print()` lists the models
  within 2 of the best and names the smallest (`within_2`,
  `smallest_within_2`).  The sweep runs on the package's own solve, so a
  candidate that is constant inside some windows no longer stops it.
  `dmat_max_n` is accepted and ignored.

* `area_of_applicability()` gains `categorical = c("refuse", "indicator",
  "drop")` for factor and character predictors.  `"refuse"` is the default,
  as before.  `"indicator"` codes each level as a 0/1 column with one shared
  scale per predictor, and a row on a level no training row has gets
  `DI = Inf`.  `"drop"` computes the index on the numeric predictors.

* `select_features_forward()` gains `null_fn`, a function that fits the
  no-predictor model for a learner that cannot.  `history` always starts
  with the step-0 row and gains `gain`, `fold_gain` and `fold_gain_se`; the
  result gains `baseline` (its type, score and the path reached) and
  `fold_scores`.

* `fit_bayesian_spatial_model()` reports what its GP basis resolves.  The
  reduced-rank basis cannot represent structure below a floor, 14 to 33
  percent of the extent along each axis at the defaults, stored in map units
  as `$info$gp_resolution` and shown by `print()`.  The new argument
  `gp_resolution` asks for a resolution in map units and derives `gp_k` and
  the length-scale prior from it, up to 50 functions per dimension.  Before
  compiling, a variogram range below the floor can raise a warning of class
  `spatialkit_gp_resolution`.  After sampling, `$info$gp_adequacy$status`
  is `"ok"`, `"below_resolution"`, `"unresolved_short_range"`,
  `"no_signal_above_floor"` or `"not_checked"`, from the length-scale
  posterior, the process's share of the variance and a short-lag Moran's I
  of the residuals; the two that say the basis is too coarse raise
  `spatialkit_gp_adequacy`.  `convergence_ok` stays about the sampler.

* `cv_bayes()`'s `fold_metrics` and `compare_models()` gain a `gp_adequacy`
  column beside `convergence_ok`.  `cv_bayes()` and `cv_spatial()` collect
  the classed warnings a Bayesian fit raises in each fold and raise each
  class once, naming the folds (fields `folds`, `fold_messages`).

* `build_tessellation()` and `create_voronoi_polygons()` gain `cover`: the
  observations, or a polygon, that the cells have to hold when the
  generators are seeds made elsewhere.  Seeds with neither a boundary nor a
  cover raise a warning of class `spatialkit_no_cover`.  `params` gains
  `boundary_source`, `n_outside` and `cover_outside`.

* `assign_features_to_polygons()` records the features no polygon holds as
  `attr(, "unassigned")` (`n`, `which`), beside `"ties"`, and names the
  column the polygon IDs were read from when `polygon_id_col` is absent.
  The new argument `quiet` silences that message.

* `fold_separation()` measures what blocked cross-validation is for: the
  distance from every held-out point to its nearest training point,
  summarised per fold.  Until now the only evidence a fold scheme worked was
  a proxy --- the block size compared against the estimated autocorrelation
  range, which `make_folds()` warns about.  That is a statement about the
  design; this is the result.  The returned table carries `n_train`,
  `n_test`, `n_blocks`, `min_dist`, `median_dist` and `mean_cor`, the
  mean correlation: the mean of `exp(-3 d / range)` over the
  held-out points, so a point at a twentieth of the range weighs more than
  one at nine tenths of it.  Nothing is estimated:
  the distances come from the geometry and the range is the one the folds
  already carry or the one you pass.  A range that records its CRS (the
  folds' own `params$sac_range` and `params$crs`, or an
  `estimate_sac_range()` result) is compared with distances measured in that
  CRS, so a copy of the layer in other units gives the same answer; `sac`
  and the distances are in the units of the CRS the `crs` attribute names,
  metres for lon/lat input.
  `fold` is the `fold_id` a `cv_*()` result's `$folds` carries, so a join on
  it pairs each fold's error with its own distances.  A `units` object for
  `sac` is refused.  As in `cv_*()`, a `make_folds()` result whose recorded
  rows sit at other locations in `data_sf` is refused; the location check is
  skipped when one layer is POINT and the other is not, so folds built on
  polygons still measure their pointized copy.  `print()` reports the
  distances and says which way the score tends to err only when it is given
  the prediction locations (`prediction_points`, below); it shows the first
  ten folds of a design with more than twenty.  The range line gives the
  CRS's unit ("in metres").

* `summary()` on a `resolution_profile()` puts every criterion's pick in one
  table: the level each prefers, the flat region around it, whether a ladder
  bound is doing the choosing, and the levels that lie in every band.  Reading
  the criteria off a profile took one `select_resolution()` call per criterion,
  and the comparison had to be assembled by hand.  The intersection of the
  bands is often empty, which is a result rather than a failure: it says the
  field has no single resolution that satisfies every way of asking.  Nothing
  in the table chooses, and the help page says so.  Each region comes back in
  full on the `"bands"` attribute.

* Every `cv_*()` result now says what became of each fold.  `fold_status`
  is a data.frame with one row per fold supplied --- `fold`, `status`,
  `message` --- where `status` is `"ok"`, `"error"` (the fit or its
  `predict()` threw; `message` is the error text), `"skipped"` (nothing
  scorable: too few matched rows, a prediction of the wrong length, no
  finite observed/predicted pair), `"dropped"` (an empty test set or fewer
  than two training rows once incomplete rows were removed, so the fold
  never reached the fitter) or `"worker_error"` (a parallel worker died).
  The per-fold error text was already collected and thrown away except when
  *every* fold failed, so a partial failure --- "3 of 5 folds produced
  predictions" --- left its causes only in console scrollback, which a
  script, a `callr` job or a knitted document does not keep.  This is worth
  most where a run is expensive: a `cv_bayes()` fold whose sampler failed
  now names the reason in the returned object.  Beside it, `orphan_rows`
  holds the row IDs no fold names (they enter no training set and are never
  scored --- non-empty only when the folds were built on a different or
  subsetted layer), `n_unknown_ids` counts the distinct rows the folds name
  that the data does not have, and `n_dropped` the rows `prep_model_data()` removed before
  any fold was fitted.  The four together account for every row and every
  fold, so `n_folds_attempted - n_folds_succeeded` never has to be explained
  from the log.  A parallel worker that was killed outright (for lack of
  memory, say) is a `"worker_error"` too and enters the first-error text; it
  used to show as `"skipped"` ("no result returned") and was left out of
  that text.  An `"ok"` row carries a `message` when the fold's
  `fold_info_fn` failed.

* `prep_model_data()` records what it removed.  `attr(x, "dropped")` is a
  list with `n`, `n_geometry`, `which` (positions in the input), `row_id`
  (when the layer carries `..row_id`) and `reason`, one of `"geometry"`,
  `"missing"` or `"non_finite"` per dropped row.  The three masks behind
  that decision were already computed and collapsed into a log line; the row
  identities reached nothing, and through the eight-plus internal call sites
  even the count was invisible, so a fit's `$n` was the post-cleaning row
  count with nothing saying how many rows were lost or why.  Every fit now
  carries the count as `$info$n_dropped`, and every `cv_*()` result as
  `n_dropped`.  Silently losing a third of the rows is a classic cause of a
  suspiciously good score.  The layer carries the class `"spatialkit_rows"`
  after `"sf"` (`c("sf", "spatialkit_rows", "data.frame")`): `[` returns a
  plain layer without the record, and `dplyr::bind_rows()`,
  `vctrs::vec_rbind()` and `dplyr::union_all()` return a plain `sf` (with
  the class ahead of `"sf"` they failed on two layers with the same record,
  'attr(obj, "sf_column") does not point to a geometry column').
  `dplyr::filter()`, `slice()`, `arrange()` and `distinct()` now drop the
  record as `[` does (they kept it, so `filter()` down to 3 of 5 rows still
  reported the parent's counts).  Each record carries `n_rows`, the number
  of rows it was made for.  `sf::st_drop_geometry()` keeps the rows and the
  record; binding such data frames keeps the first one's record, whose
  `n_rows` then no longer matches, and the package's readers ignore it.

* `make_folds(method = "block_kfold")` returns the block design it built the
  folds from.  `assignment` gains a third column, `block_id`; `params` gains
  `blocks` (an `sf` layer of the block polygons in the CRS the folds were
  built in, numbered to match), `block_sizes` (points per block, indexed by
  `block_id`, so a block kept empty by `drop_empty_blocks = FALSE` shows as
  a zero) and `fold_blocks` (which blocks were packed into each fold).
  `blocks$source_row` is the row each block came from in the layer it
  originated in --- a cell's index in the full `grid_nx` by `grid_ny` grid,
  or the row of the `blocks` argument --- because dropping the empty blocks
  renumbers the rest: nine supplied blocks of which three are empty come
  back as six rows numbered 1 to 6, and a join by row position would
  mis-attribute every block after the first gap.
  `blocks[params$blocks$source_row, ]` recovers them with their own columns
  and in their own order.  The folds account for every block exactly once,
  empty ones included, so a fold's territory on the map is all of its blocks
  rather than only those holding points; and `blocks_used` is the number of
  blocks the design has, which equals `nrow(params$blocks)` and
  `length(params$block_sizes)`, with `sum(params$block_sizes > 0)` giving how
  many of them hold points.  The
  grid **is** the design of a blocked cross-validation: without it a user
  could not draw the blocks over their data, see that 40 of 64 blocks were
  empty, or tell whether a fold is one contiguous region or several.
  `plot_folds()` now draws those outlines under the points when the folds
  carry them, and takes `blocks = FALSE` to suppress them.  Its subtitle
  states the parameter that sets how far the held-out points sit from
  their training data (the block size, the buffer, the number of location
  groups or the median NNDM exclusion), with the units of the folds' CRS;
  folds built from points without a CRS get no units rather than "(NA
  units)".

* `estimate_sac_range()`'s returns are uniform.  The four directional
  variograms and their fits run unconditionally on every call, and the
  rejected-range paths used to discard them --- precisely where a user most
  needs to know whether the field is anisotropic.  Every classed return now
  carries `directional`, `anisotropy` and `anisotropy_used`, and three new
  attributes report the sweep rather than collapsing it: `directional_status`
  (per azimuth, why that direction is `NA` in `directional` ---
  `"ok"`, `"over_cutoff"`, `"not_converged"` or `"no_fit"`, which were
  indistinguishable before), `directional_fitted` (the range each direction's
  fit reported whether or not it was usable, so a refused directional range
  --- the most informative number in an anisotropic failure --- stays
  recoverable) and, under `keep_directional_fits = TRUE`, `directional_fits`
  (each direction's empirical variogram and fitted model).  Those four
  objects are off by default because they dominate the result when present
  --- 42.2 KB of a 58.8 KB object at n = 400, against 16.6 KB without them
  --- and `make_folds(auto_range = TRUE)` calls this on every build.
  `print()` names the reason and the refused value for a direction it
  cannot use.

* A roster of quantities the package already computed and dropped are now
  returned.  `summarize_by_cell()` attaches the Kish ICCs it estimated as
  `attr(, "icc")` whether or not either was large enough to apply --- the
  case with no `"deff_applied"` is exactly the one where a user wants to
  know what the ICC came out as --- and the `"variogram"` path adds
  `deff_rows`, the per-cell design effect at the cell's row count that the
  log line reduced to a median and a max.  A `summarize_by_cell()` call that
  requests a design effect (any `deff` other than 1) also returns a logical
  `deff_applied` column, `TRUE` on every row when a correction was applied
  (a Kish fallback included; `deff_method` names the correction) and `FALSE`
  when none was (the standard errors are then the uncorrected ones, or `NA`
  where a Kish ICC cannot be estimated); unlike
  the `"deff_applied"` attribute it survives `rbind()` and
  `dplyr::bind_rows()` of many results.  At the default `deff = 1` the
  column is not added.
  `assign_features_to_polygons()` reports the features that matched more
  than one polygon and had the `tie_break` rule decide for them, as
  `attr(, "ties")` and a log line; a tie-break firing on a third of the
  features means the polygon layer overlaps and every cell count built from
  it is suspect.  The layer carries the `"spatialkit_rows"` class after
  `"sf"`, as `prep_model_data()`'s does, so binding such layers returns a
  plain `sf`.  The `"ties"` record is dropped by the same `dplyr` row verbs
  as `prep_model_data()`'s, and with `largest = TRUE` it also counts
  features that overlap two or more polygons by exactly the same area.
  `ensure_projected()` attaches `crs_choice`, the projections
  it considered with each one's measured worst-case distance error.  Every
  path that picks a local projection reports what it picked, the two that
  compare nothing included: a UTM zone on a local extent, and the
  equal-area projection chosen for a layer straddling the antimeridian.
  `residual_morans_i()` returns the residual `kurtosis` the randomisation
  variance conditions on, the design rank `p` behind the residual moments,
  `exact`, whether those moments are exact for these residuals, and
  `weights_summary`, a description of the weight matrix it used: `n`,
  `storage` (the matrix class), `neighbours` (the smallest and largest
  number of neighbours any row has, `NA` for a dense matrix, where counting
  them would allocate a second one), `kept` and `desc`, the line `print()`
  shows.  The matrix itself comes back as `weights` under
  `keep_weights = TRUE` and is `NULL` otherwise, because it is n by n: at
  n = 500 it is 50.3 KB of a 53.5 KB object in its sparse form, 112.7 KB at
  n = 120 when the dense fallback is taken (going sparse needs both **FNN**
  and **Matrix**, so a no-Suggests install always falls back), and 191 MB
  at the n = 5000 that fallback is capped at --- against the 3.2 KB
  everything else occupies --- and scoring a list of fits would hold one
  matrix per fit.  The result is classed `"morans_i"` and prints through a
  `print()` method that shows the statistic, its null and the weights line;
  `[` drops the class, and `$`, `[[` and `unlist()` read the result exactly
  as for a plain list.  `fit_gwr_model()` keeps `info$nonfinite_coef`, the
  per-row, per-term mask behind `n_local_singular`, and
  `fit_bayesian_spatial_model()` keeps `rhat_failed` / `neff_failed`, the
  parameters that failed each convergence check by name --- "max R-hat 1.09"
  is not actionable where "`sdgp_gp..x..y` has R-hat 1.09" is.
  `determine_optimal_levels()`'s model-aware diagnostics gain `knee_k` and
  `failed_k` (whose interpolated WSS entries are not measurements).
  `build_tessellation()` records the points whose cell assignment was
  repaired by nearest-cell snapping, and how far outside each sat, as
  `params$snapped` --- a comment had long said it should.
  `area_of_applicability()` returns the `scaling` (per-predictor training
  centre and SD) the dissimilarity index is computed in, without which a
  location's DI cannot be traced to the predictor that put it outside, and
  `n_outliers`, the training DI values the threshold's fence set aside.
  `get_voronoi_seeds(method = "kmeans")` returns the clustering as
  `attr(, "kmeans")` --- which cloud points fed which seed, the cluster
  sizes and the within-cluster sums of squares.  `gwr_model_selection()`
  reports `criterion_by_name`, `criterion_column` and `criterion_verified`,
  so a script can gate on the case its log calls "unverified" instead of
  reading the label.

* `estimate_sac_range()`'s result gains a `plot()` method.
  `plot(estimate_sac_range(pts, "z"))` draws the empirical variogram, with
  the fitted model and the effective range overlaid where a range was
  identified, and a subtitle saying why not where it was not: the variogram
  never reached a sill, both model fits were singular, the optimiser
  halted, or a `range_frac` below 1 refused a range inside the fitted lags.
  Nothing is recomputed --- the plot reads the attributes the
  estimate already carries --- and the same drawing routine now serves
  `plot.spatial_fit(type = "variogram")`, so the two pictures agree.  The
  "attached for inspection" messages `estimate_sac_range()` logs when it
  returns `NA` used to point at `plot(type = "variogram")`, which is the
  method for a fitted model and could not take the estimate; they now point
  at `plot()` on the returned value.

* `summarize_by_cell()` gains `conf_level`.  With `conf_level = 0.95`, every
  numeric response and predictor column gets four more columns beside its
  `..sd_*` and `..se_*`: `..neff_*`, the column's effective sample size in
  the cell (its non-missing count over its design effect --- the per-column
  version of `cell_weight`); `..df_*`, the degrees of freedom the interval
  uses; and `..ci_lo_*` / `..ci_hi_*`, a t interval for the cell mean as an
  estimate of the grand mean, built on the design-effect-corrected standard
  error.  The default `conf_level = NULL` adds none of these columns.  The
  degrees of freedom are `n - 1` at `deff = 1`, for a numeric
  `deff` and for `deff = "kish"`, because the interval's spread is estimated
  from the within-cell variance, whose `n - 1` df survive exchangeable
  correlation whatever the design effect (measured 95% coverage on the Kish
  path at an ICC of 0.2 / 0.6 / 0.9: 0.954 / 0.953 / 0.952; an
  effective-sample-size df of `neff - 1` gives 0.992 / 1.000 / 1.000 and is
  not used).  Under `deff = "variogram"` the df are the Satterthwaite
  (1946) moment-matched df of the within-cell variance under the fitted
  correlation, a fraction of `n - 1` that shrinks with the range: 0.960 and
  0.958 coverage at exponential ranges of 150 and 400 on a 1000-unit
  domain, against 0.931 and 0.918 with `n - 1`.  The help page's
  "Confidence intervals" section has the reasoning and the numbers.
  `..neff_*`, like `..df_*` and the interval, is `NA` for a column with a
  single non-missing value in the cell, where `cell_weight` still counts it.

* `compare_models_cv()` gains `block_size` and `auto_range`, and now hands
  `response_var` and `predictor_vars` to `make_folds()` when it builds the
  shared fold set.  With the defaults the folds are the same geometric
  blocks as before, so an existing comparison does not move; what changes is
  that the fold-leakage diagnostic (above) can now fire for the one function
  that compares models, which until now was the one whose folds could never
  be checked against the range.  `select_features_forward()` gains
  `auto_range` for its inner folds for the same reason, and records it in
  `$params`.  A shared fold set that cannot be built is an error.  It used
  to be a log line and a fall-back to each backend's default blocks, which
  were given neither argument: `block_size = 1e6` (one block, which
  `cv_rf()` refuses) returned a five-fold comparison on geometric blocks,
  and `auto_range = TRUE` the small blocks it exists to prevent, with no R
  condition.
* `compare_models_cv()$overall` carries the Bayesian backend's calibration
  when a Bayesian model ran: `coverage_50`, `coverage_80`, `coverage_95` and
  `mean_CRPS`, the same fold-weighted summary `cv_bayes()` returns as
  `predictive_coverage`, with `NA` on the GWR and RF rows.  A model that
  predicts well on average and covers badly (Heaton et al. 2019) is now
  visible in the table a user picks from, not only in `$bayes_cv`.  `model`
  stays the last column.
* `fit_rf_model()` gains `replace` and `sample_fraction`, which reach
  `ranger::ranger()` as `replace` and `sample.fraction`; both were already
  accepted through `...`, but are now recorded in `$info` and printed with
  the fit ("Sampling: bootstrap, with replacement (100.0% of rows per
  tree)").  The defaults are ranger's, so no forest changes.  The help page
  carries Strobl et al.'s (2007) case for `replace = FALSE`.  Passing the
  ranger spellings through `...` is now refused like the other arguments the
  wrapper sets.  `replace = FALSE` with `sample_fraction = 1` grows every
  tree on every row, so no row is out of bag.  The fit now warns that
  `fitted()`, the OOB error and the permutation importance are all `NaN`
  (they were `NaN` in silence), and `print()` on such a fit says the OOB
  error and the permutation importance are undefined, where it dropped the
  OOB line and printed the importance line empty.
  `area_of_applicability(weights = pmax(imp, 0))` failed on that importance
  with a hint to use the very `pmax()` it had been given, because `pmax()`
  keeps `NaN`; it now names the predictors whose weight is `NaN`, says why
  (no row is out of bag), and says what to do instead: refit with out-of-bag
  rows, or pass `weights = NULL`.  With permutation importance, the fit's
  own warning now says as much: `area_of_applicability()` cannot be
  weighted by that importance, so pass `weights = NULL`.
* `area_of_applicability()` records the method of the folds its threshold
  came from as `$params$folds_method` (`"block_kfold"`, `"random_kfold"`,
  ... from a `make_folds()` result; `"labels"` or `"splits"` when the input
  cannot say) and prints it, because the threshold is a statistic of a
  cross-validated hold-out and pairs with the CV error from the same kind
  of hold-out.  An AOA built on `random_kfold` folds logs a caution saying
  it does not belong beside a blocked `cv_*()` result.

* `estimate_sac_range()` gains `detrend = c("ols", "reml")` and
  `reml_max_n`.  A variogram fitted to least-squares residuals
  underestimates the range, because the trend fit absorbs part of the
  long-wavelength variation (Lark, Cullis and Welham 2006).  Measured for
  this estimator on simulated fields (n = 300, true effective range 300),
  as the median ratio to the estimate from the true field: a white-noise
  covariate 1.00; a spatially smooth covariate 0.97; a linear trend in the
  coordinates 0.92; a quadratic one 0.75.  `detrend = "reml"` fits the trend
  and an exponential-plus-nugget covariance together by REML with
  `nlme::gls()` and returns the REML range --- 0.95 to 1.06 on the same
  designs --- with the empirical variogram of the REML residuals attached
  for inspection.  It is cubic in `n`, so it runs on at most `reml_max_n`
  (400) points; a fit that does not converge falls back to OLS with a
  warning.  The default stays `"ols"`, so nothing changes unless asked; the
  help page's new section carries the numbers.  Iterating GLS trend fits
  against variogram refits (Neuman and Jacobson 1984) was measured too and
  recovers only part of the bias (0.80 in the quadratic case), so it was not
  added.  `nlme` joins Suggests.  A fitted range too short for enough pairs
  of points to lie inside it is refused, with
  `rejected_reason = "fitted range is below the shortest lag fitted"` and
  the bound it fell short of as the `range_floor` attribute.  For the
  least-squares fits the bound is the shortest lag the empirical variogram
  resolves (the mean separation in its first bin); it did not fire on white
  noise or on fields with a 60 m range.  The REML range is fitted to the
  point pairs, not to those bins, so its bound is the distance within which
  30 pairs of the points the REML fit used lie, or that first lag when it is
  shorter.  On white noise (n = 300, 30 draws) the REML fit returned ranges
  of 0.18--23.6 m in 19 draws (one of 0.27 m sized a 3642 x 3676 block grid,
  refused as a unit mistake); the 16 of them up to 12.5 m are refused and
  the estimate is finite in 13 of 30.  REML estimates of a true 30 m range
  that the first lag alone would have refused (14.5--28.8 m, in 10 of 20
  draws at n = 400 and 6 of 15 at n = 300) are returned, and at n = 400 the
  same holds at `cutoff = 0.5` and `0.1`.  The bound is about
  identification, not a test for spatial structure.  A refused REML result
  keeps its `reml` list.
* `sac_nugget()` returns the nugget variance behind an
  `estimate_sac_range()` result, and every classed result --- identified or
  rejected --- now carries it as a `nugget` attribute (`NA` when no model
  could be fitted).  It was reachable before only by reading the `Nug` row
  of the attached `gstat` model.  Results also record `detrend_method`
  (`"ols"`, `"reml"` or `NA`), and with `detrend = "reml"` a `reml` list
  (`n_used`, `subsampled`, `nugget_prop`, `sigma2`).
* `estimate_sac_range()` refuses a range fitted through an empirical
  variogram that *decreases* with distance over its shorter lags: a net fall
  of more than 15% of the mean semivariance there, weighted by pairs.  A
  model that rises to a sill has nothing to identify on such a curve, and
  the result is `NA` with `rejected_reason = "empirical variogram decreases
  with distance"`, the refused value in `rejected_range`, and the variogram
  attached; `plot()` captions it.  The shape is what a periodic
  (hole-effect) structure or a variance that differs between a dense cluster
  and the rest of the layer produces --- measured on 60 draws each: 98% of
  fields with a periodic component and 100% of the clustered case are
  flagged, against 0% of an exponential field with a range of 100 or more
  on a 1000-unit extent, 2--3% at very short ranges, 2% of white noise ---
  and *not* what a trend produces, which is a variogram that rises without a
  sill and is refused as before.  The one other path that returned a bare
  `NA` from inside a completed fit (a non-positive fitted range) now returns
  the classed, inspectable shape too.

* `resolution_profile()` and `select_resolution()`: the number of cells,
  scored on every criterion at once.  Given a ladder argument (`n_levels`,
  for one; at the defaults it scores one size, see Breaking changes),
  `resolution_profile()` runs a log-spaced ladder of level counts from 2 to
  a ceiling the support implies
  (`floor(n / min_cell_n)`, where `n` is the rows of the layer that carry a
  response, not the `sample_n` subsample the k-means runs on), reports the
  floor the autocorrelation range implies (`ceiling(area / range^2)`),
  fits each level as the best of its partition's k-means++ starts (25 for
  `partition = "kmeans"`), and returns a data.frame with the WSS elbow
  statistic (read on log-log axes, and `NA` when the points have no cluster
  structure to bend the curve), Mallows'
  C_p of the piecewise-constant approximation of the response (or of
  its residuals on the predictors, fitted on the rows with a complete
  response and predictors, with the nugget from the fitted variogram as the
  noise variance and the penalty set for the whole layer), the
  standardised residual Moran's
  z of the cell means (a diagnostic, not a count rule), and an analytic
  reliability of the cell means
  --- the share of their spread that is between-cell signal rather than
  sampling noise, from the variogram alone via Krige's additivity relation
  (Cressie 1996), the shrinkage factor of Fay and Herriot (1979) --- plus
  the cell-support and cell-diameter columns and the between-restart
  spread.  A floor above the ceiling is reported as a finding rather than
  resolved silently.  `select_resolution()` reads a level off one criterion
  together with its flat region and says when a bound, not the criterion,
  is choosing.  The flat region is a set, not an interval: the criterion
  curves are not monotone, so a region can skip a rung of the ladder, and it
  is printed as the runs the criterion accepts ("19 to 21, 26, 31 to 33")
  rather than a range that would quietly include the levels it rejected.
  Two things measured before this shipped, both on the help page: on smooth
  fields with a small nugget C_p descends to the support ceiling (every
  replicate at effective ranges 90--900 with nugget 0.3 on a unit sill;
  interior only at nugget 2), and the reliability optimum agrees with the
  empirical one from true block means on simulated fields but is broad ---
  flat to within 2 percent over a factor of 3--6 in the number of cells.
  Read the flat region.
  `determine_optimal_levels()` is unchanged in shape and keeps its
  integer-vector interface.  A supplied `sac` whose range was refused (an
  `NA` with a `rejected_reason`) is not read as accepted: `cp` keeps that
  fit's nugget and warns, naming the reason, `reliability` is `NA` (it
  needs the range), and a fit that did not converge, or whose range is below
  the shortest lag fitted, gives neither its nugget nor its range (the whole refused fit used to be
  used, and a correlation function whose range could be many times the
  extent pinned the reliability optimum to the first level, with no R
  warning).  A range below the shortest lag means the structure cannot be
  told from a nugget, so the nugget is not identified either: on white noise
  detrended by REML it was 6e-7 on a sill of 0.99, and C_p, with no penalty,
  ran to the support ceiling (33 cells).  A nugget of 0 (under 1e-4 of the
  sill, since a REML fit stops short of its bound: 6e-7 passed a test for
  exactly 0), on which C_p has no penalty and falls to the ceiling, is
  warned about, and so is a `sac` whose `detrended` flag does not match the
  variable scored (a residual variogram on the raw response moved the C_p
  pick from 2--4 cells to the ceiling of 44 in five of five simulated
  fields); `attr(x, "variogram")` records `detrended`.  When no `sac` is
  passed, these warnings name the variogram the profile estimated (kept in
  `attr(x, "sac")`), not a `sac` argument the caller never gave.  A supplied
  `sac` is read in its own CRS: the points are transformed to
  `attr(sac, "crs")` first, as `summarize_by_cell()` does (a range in US
  feet put the floor of a metre layer at 2 where the same range in metres
  put it at 8).  A `sac` given as a `units` object (`set_units(1.5, "km")`)
  or a character string is refused by name; a `units` object used to be
  read as a number in the CRS units, so 1.5 km became a range of 1.5 m and
  a floor of 41 million cells.  A plain number is taken as the range alone:
  it sets the floor, and `reliability` is `NA`.  Wherever the variogram
  gives no nugget (none was fitted, its fit did not converge, its range is
  below the shortest lag, or `sac` is a range alone), `cp` takes Mallows'
  own noise variance, the residual mean square RSS / (m - L) of the finest
  level whose cells hold two scored rows each on average, and the profile
  warns; `attr(x, "cp_noise")` records which variance was used, and the
  print says so.  It counts the structure within those cells as noise too,
  so on average it is no smaller than the nugget and errs towards fewer
  cells.  `cp` used to be `NA` at every level in all four cases, so a
  workflow that reads C_p by default (`build_tessellation(approx_n_cells =
  <profile>)`, `summary()`) fell through to another criterion with a log
  line at most.
  Reliability's domain term is taken over the convex hull the
  area is measured on, not the bounding box: on a 3000 x 120 strip the
  reliability pick is 6 cells whether the strip lies axis-aligned or rotated
  by 45 degrees (it was 8 and 2), and 10 either way on a square (it was 10
  and 8), and reliability values shift slightly on every profile.  Row order
  does not change the profile (see the `determine_optimal_levels()` item
  under Bug fixes); over permutations of one 2000-point layer, WSS used to
  move by up to 2.6 percent and the C_p pick across its flat region (222,
  173 and 135 cells).  The ceiling reaches the number of distinct locations
  when locations repeat (and explicit `levels` up to it are kept) instead of
  stopping one short; with no location repeated, the print says the ceiling
  is one short of the points rather than crediting the distinct locations.
  Points with empty or non-finite coordinates are dropped with an R warning,
  as in `determine_optimal_levels()`; it was a log line alone.
  `range_floor = FALSE`, the default, starts the ladder at 2 whatever
  the range and only reports the floor, so profiles whose range estimates
  differ (cross-validation folds) have comparable ladders; with the floor
  applying on one fold and not the next, reliability and the elbow moved by
  a factor of 10--15 between folds.  `range_floor = TRUE` starts the ladder
  at the floor.  Under
  `select_on = "split"` a supplied `sac` is flagged in the log, since it
  must be fitted on the selection half, and the help page shows the
  two-call workflow.

* `select_on = c("all", "split")` on `determine_optimal_levels()`,
  `resolution_profile()` and `select_features_forward()`.  Whenever a
  selection reads the response --- a level count chosen under
  `criterion = "morans_i"`, or a predictor set chosen by a
  forward sweep --- what is estimated afterwards on the result is
  post-selection, and its standard errors are descriptive rather than at
  nominal coverage (Gao, Bien and Witten 2022).  `select_on = "split"` is
  sample splitting: the layer is cut into two spatially blocked halves
  (`make_folds(k = 2, method = "block_kfold")`), the selection reads the
  response on the first only, and the row positions of both come back (as a
  `"split"` attribute on the first two functions, as `$split` on the third)
  so the estimation can be done on the half the selection never saw.  The
  positions index the layer as passed, for all three functions
  (`select_features_forward()`'s were positions after its completeness
  filter, so with seven incomplete rows `pts[fs$split$estimation, ]` held 47
  rows of the selection half and 4 of the dropped ones).  When the default
  six-block split leaves fewer than 10 points in a half (a small layer, or a
  small group far from the rest), it is retried on grids of about 16, 36
  and 100 blocks, with a warning, instead of stopping, and the split records
  `grid`, `n_blocks`, `balance` (a 2:1 split on clustered layers used to
  pass silently) and each half's `extent`.  The halves are not independent:
  blocking reduces the dependence across their shared border but does not
  remove it (85 percent of the estimation points of the package's test
  layer lie within the fitted range of a selection point), and the seed
  decides only which side selects, not where the cut falls.  The level-count
  functions still draw their cells on every point, so the count they return
  is a count for the whole layer it will be applied to.
  `select_features_forward()` also returns `score_holdout`: the selected set
  fitted on the selection half and scored on the other, the honest number
  its selection-internal `score` is not.  The cost is precision --- half the
  points estimate, and a contiguous spatial half is less efficient than an
  exchangeable one (García Rasines and Young 2023).  Data thinning (Neufeld
  et al. 2024) and data fission (Leiner et al. 2023) keep the whole sample
  and are noted on the help page, not implemented.  Under
  `criterion = "combined"`, which no longer ranks by the response (see
  Breaking changes), the split changes only the rows its diagnostics read.

* `make_folds()` gains `blocks`: a polygon layer (`sf` or `sfc`) to use as
  the blocks of `method = "block_kfold"` in place of the grid it would
  otherwise build --- the `$cells` of a `build_tessellation()` result,
  hexagons, watersheds, administrative units, the `$blocks` of
  `blockCV::cv_spatial()`.  Each point takes the block that contains it and
  the blocks are assigned to folds exactly as grid cells are, so the fold
  builder can now consume every shape the tessellation half of the package
  produces.  The grid-sizing arguments and `boundary` are ignored with a
  log line, `auto_range` compares the estimated range against the blocks
  instead of resizing them, and the leakage warning uses the median over
  blocks of the side of the square with the block's area
  (`params$block_scale`).  Points inside no block are assigned to the
  nearest one with a warning that counts them, except points within a
  millionth of the extent of a block --- an edge that reprojection moved by
  a rounding error; points inside more than one block take the first, with
  a warning when the blocks concerned overlap in area rather than share an
  edge.  `params` gains `n_blocks` (before empties were dropped; for a grid,
  `grid_nx * grid_ny`, the cells a `boundary` clips away included),
  `blocks_supplied` and `block_scale` on every `block_kfold` result;
  `grid_nx`/`grid_ny` are `NA` for supplied blocks.  A `blocks` layer
  without a CRS is brought into the points' CRS with an R warning naming
  it, as `boundary` is.

* `make_folds()` gains `balance_tol`, the largest-to-smallest fold size
  ratio above which `block_kfold` reports its folds as imbalanced.  The
  check was always there at a hard-coded 3:1, and it was a log line only;
  it is now an R warning (`Inf` disables it), the ratio achieved is
  returned as `params$balance_ratio`, and the help page says which
  methods balance what: only `block_kfold` balances point counts, by
  packing blocks largest first into the fold with the fewest points so far.
  No search over packings was added, because the packing is not where the
  imbalance comes from.  Measured against the optimum by enumeration (two
  folds, up to ten blocks, heavy-tailed sizes), the greedy packing is
  optimal in 72 percent of cases and within two points of optimal on
  average; a local search with 30 random restarts moved the ratio by 0.004
  on average over 480 block-size vectors and never brought one of the 91
  above 3:1 below it.  The remedy for an imbalance past the tolerance is
  the block design --- the warning now says so --- and `blocks` is how to
  supply one: on clustered layouts where the geometric grid exceeded 3:1
  in 22 percent of draws (median 1.8, worst 6.3), Voronoi cells around 15
  `get_voronoi_seeds(method = "kmeans")` seeds never exceeded 1.3 (median
  1.09).  With the default tolerance the folds of every existing call are
  unchanged.

* `metrics` on `cv_spatial()`, `cv_gwr()`, `cv_bayes()`, `cv_rf()` and
  `compare_models_cv()`: a scoring function of your own, `function(y,
  yhat)` returning a named numeric vector (a named list of scalars or a
  one-row data frame also serve), applied the way the built-in metrics are
  --- once per fold, so each name becomes a column of `fold_metrics`, and
  once to the pooled out-of-sample predictions, so each name becomes a
  column of `overall` --- on the same finite pairs `RMSE` uses.  This is
  the way to score what the Gaussian set cannot: a Poisson deviance, a log
  score on a probability, a weighted loss.  The contract is strict where it
  should be (every element named, names unique and not a built-in column,
  one number per name --- anything else is an error, because a scoring
  function of the wrong shape is a mistake to surface) and forgiving where
  it should be (a function that throws on a fold is logged and its columns
  are `NA` there; a fold is never dropped for it).  Nor may a name reuse a
  per-fold extra (`CRPS`, `coverage_*`, `gp_k`, `gp_n_basis`, `n_draws`,
  `bandwidth`, or a name the `fold_info_fn` returns) or `mean_CRPS`: each
  used to replace the package's value silently, and
  `cv_bayes()$predictive_coverage` then reported the user's number.  The
  empty frames of a run where every fold failed carry the columns, typed,
  when the function can be called on zero-length input.
  `compare_models_cv()` hands one function to every backend and protects it
  like the fold arguments, so the columns of its `overall` are comparable
  across rows.  `fold_info_fn` is
  documented as the per-fold half of the same mechanism, with access to the
  fitted object and the held-out layer.

* `ensure_projected()` gains `purpose = c("distance", "area")`.  The
  default is what it always did: for lon/lat input, the candidate that
  distorts distances least.  With `"area"` --- densities or rates per cell
  are going to be computed --- the choice is made among equal-area
  projections only (a Lambert azimuthal centred on the data, or an Albers
  conic where its parallels do not degenerate, whichever distorts distances
  less), which a UTM zone never enters, and global coverage gets Equal
  Earth rather than Web Mercator.  Already-projected input is still
  returned untouched, but its area distortion over the extent is now
  measured --- the spread of planar-to-geodesic area ratios over probe
  polygons --- and logged as a warning above 1 percent.  Measured: a UTM
  zone edge to edge 0.25 percent, a 2.5-degree extent inside one 0.04
  percent, the conterminous United States forced into one zone 14 percent,
  Web Mercator over 2.5 degrees of latitude at 48N 4 percent, an
  equal-area projection a few tenths of a percent (the sphere the geodesic
  areas are computed on against the ellipsoid).  The measurement works
  whether `sf_use_s2()` is on or off (with it off and no lwgeom it returned
  `NA`, so `summarize_by_cell(area = TRUE)` refused every grid and
  `purpose = "area"` skipped its warning), and its probe polygons are
  densified before their geodesic area is taken, so an equal-area grid over
  a near-global extent no longer measures 10 percent.  With
  `purpose = "area"`, a single polygon's candidates are scored rather than
  falling back to the Lambert azimuthal unmeasured.

* `summarize_by_cell()` gains `area = TRUE`: with `cells_sf`, the result
  carries `cell_area` (planar, in the squared units of the cells' CRS; for
  lon/lat cells the geodesic area in square metres, not a planar area in
  squared degrees) and `n_per_area`, a point density; a rate of anything
  else is its `agg_funs` sum over `cell_area`.  The request is refused with
  an error --- not answered with a number --- when the cells' CRS distorts
  areas across them by more than 1 percent by the measurement above, because
  a density is a comparison between cells and means nothing where the map
  scale differs from one cell to the next; the message names the CRS, the
  figure and the remedy.  A cell with no observations gets `NA`, not zero.
  Lon/lat cells pass the distortion check by construction; with s2 switched
  off their area needs lwgeom, and the request is refused without it.  The
  measured spread is attached as `attr(, "area_error")`.  A `cells_sf`
  with no ID column the summaries can be joined on is an error under
  `area = TRUE`, since the area columns cannot be produced.

* `build_tessellation(approx_n_cells = )` and `get_voronoi_seeds(n = )`
  accept what the level-selection step returned: the integer vector of
  ranked candidates from `determine_optimal_levels()` (its first element is
  used), a `select_resolution()` result (its `$best`), or a
  `resolution_profile()` (read by `cp`; a profile with no `cp` is an
  error).  Both were hard errors before, so nothing that worked
  changes; the count used and where it came from are recorded as
  `params$approx_n_cells` / `params$approx_n_cells_from` and as
  `attr(seeds, "n_from")`.  The two functions the pipeline documents as a
  pair are now connected: `get_voronoi_seeds(n = determine_optimal_levels(pts))`
  needs no number carried between the calls by hand.  A hex or square
  lattice records `params$cells_occupied` and `params$cells_empty`, and a
  count taken from a profile or a selection warns when fewer than three
  quarters of it end up occupied: the count is of k-means cells, all of
  them occupied, and on clustered points a lattice leaves about half its
  cells empty.

* Five diagnostic plots that show the curve behind a chosen point, the
  folds behind a pooled number, or the distribution behind a count.  None
  recomputes anything; each draws what the result already carries.

  - `plot_cv_metrics(cv, metric)`: one point per fold, sized by the
    held-out rows it contributed, with the pooled value from `overall` as
    a dashed line; a `compare_models_cv()` result gets one panel per model
    on a shared scale.  Any column of `fold_metrics` can be drawn,
    including backend extras and columns a `metrics` function added; a
    column that is `NA` in every fold is refused with the reason (`Adj_R2`
    without `p`, coverage without draws) rather than drawn empty, and a
    per-fold extra with no pooled counterpart draws without the line and
    says so.  So does a count (`n_pred`, `n_MAPE`, `n_SMAPE`), whose
    `overall` value is the total over the folds; it was drawn as the pooled
    line, at 150 against folds of 30.  A model with no finite per-fold value
    gets no panel and no pooled line, and the caption names it, whether or
    not `overall` has a value for it (RF's `bandwidth` was dropped without a
    mention).  A `compare_models_cv()` result keeps its model strip when
    only one model is left to draw.
  - `plot.aoa()`: the dissimilarity index of the prediction locations
    against the training DI --- cross-validated over the `folds` passed, or
    else each point's distance to its nearest other training point; the
    legend says which --- as ECDFs, or a histogram with the training curve,
    threshold marked, with the share outside and how close the inside ones
    run to the edge in the subtitle, and whether the threshold came from
    cross-validated folds in the caption.  Prediction locations outside on a
    predictor dropped for having no training variance (`DI = Inf`) count in
    the prediction curve, which then tops out below 1, and the caption says
    how many are off the axis.  It counts the `DI = NA` rows (a missing
    predictor) as well.  The curve used to leave the `Inf` rows out: it read
    0.97 inside at the threshold while the subtitle counted 11 of 40
    outside.  A result with every row at `DI = Inf` was refused as "every
    row had a missing or non-finite predictor", and the error now gives the
    true reason.
  - `plot.spatial_fit(type = "variogram")` overlays the response's own
    variogram (hollow points, dashed fit) on the residual variogram, on the
    same points and lags, so the structure the model absorbed is the gap
    between the two curves.  The response curve is whichever variogram
    `estimate_sac_range()` returns for the response, which is its widest
    single direction when the all-pairs fit is unusable (a response with a
    trend, whose residuals are fine).  A single-direction response curve is
    labelled with its azimuth, and the caption compares the sills only when
    both ranges were identified and the two curves cover the same point
    pairs (such a curve used to be labelled plainly as the response, its
    sill set against the all-pairs residual curve's: "Residual sill is 45%
    of the response sill" from a quarter of the pairs).  `response = FALSE`
    restores the residual curve alone.
  - `plot_calibration(cv)`: observed against nominal coverage of
    `cv_bayes()`'s posterior predictive intervals, pooled (blue) and per
    fold (grey), with the diagonal and a one-line verdict.  Each level is
    drawn at the nominal value `cv_bayes()` records in `coverage_levels`
    (0.995 was drawn at 1.00), so `coverage_levels =
    seq(0.1, 0.9, by = 0.1)` gives a full curve; the default three levels
    are unchanged.  Its error for all-`NA` coverage names
    `compute_pred_intervals = FALSE` as well as failed draws.
  - One sweep drawer behind three methods: `plot.resolution_profile()` (a
    panel per criterion, the level each selects marked, its flat region
    shaded, a note when a bound rather than the criterion is choosing);
    `plot.feature_selection()` (the accepted variable's score at each step
    as the path, every other candidate faint, the stop in red, the best
    candidate at the step after the stop, which was scored and not added,
    hollow and labelled "not added" rather than drawn like the accepted
    variables, the hold-out score as a separate mark when
    `select_on = "split"` computed one --- and a caption saying whether the
    intercept-only model was scored, since for the RF and GWR backends it
    usually is not, so the path starts at the first variable);
    `plot.gwr_model_selection()` (every model's AICc against its size, the
    best of each size joined, the winner marked, its lead over the
    runner-up in the subtitle).
    `select_features_forward()`'s result now carries class
    `"feature_selection"` so `plot()` finds the method; it is the same list
    otherwise.

* `cv_block_size_sweep()`: the same cross-validation at a ladder of block
  sizes, with random folds as the reference, returned as a table with the
  fold-to-fold spread at each size and the estimated autocorrelation range
  alongside; `plot()` draws the curve with the range marked.  The curve
  rises as the blocks widen because the held-out points sit farther from
  their training data: its height at a size is the error at that distance
  from the data, which is the error of a map only at locations that far
  from it.  On simulated fields with a 100--130-unit range and a random
  forest with coordinates, the curve levels off from one to two times the
  estimated range.  Each size is a full `cv_spatial()`, so a fit budget
  (`max_fits`, default 60) refuses to start rather than run past it.  The
  default ladder runs up
  to the largest size, at most half the shorter side, whose grid still holds
  `k` cells.  At the default `k = 5` half the side is a 2 x 2 grid of four
  cells on any extent less than 1.5 times as long as it is wide, so that
  rung used to be dropped every time.  `n_sizes = 6` then ran five
  cross-validations, the last at 0.30 of the side, and missed ranges up to a
  third of it that a 3 x 3 grid reaches.  The top is now about a third of
  the side on such an extent.  Sizes the caller passes in `block_sizes`
  whose grid holds fewer than `k` cells are not run, and a warning names
  them and the largest size that gives `k` blocks; they were dropped with
  only a log line.  A `units` object for `block_sizes` is refused by name.
  On clustered data `make_folds()` can still lower `k` at the top sizes,
  which the `k` column shows; the sweep warns once, naming those sizes.  The
  default ladder also skips sizes whose
  grid would exceed the 1,000,000 blocks `make_folds()` builds (a 6 km x 2 m
  transect used to abort with an error about `block_size`), runs along the
  line for points on one axis-parallel line (they were refused as having "no
  extent"), and warns when every rung is below the estimated range while
  longer blocks would still fit `k` times along the longer side (a 10 km x
  100 m corridor: rungs of 4 to 50 m against a 1.7 km range, a flat curve
  that read as no leakage).  On a roughly square extent at `k = 5` that
  warning also fired, as in the block-size tour script on a 996 m square,
  with two wrong numbers.  It called the top rung it had run (300) "half the
  shorter side" (498).  It said blocks only up to the longer side over `k`
  (199, smaller than rungs already run) still gave `k` blocks, when blocks
  up to 332 did.  It now fires only where a longer block would still give
  `k` blocks, and it names the top rung run and the largest size that gives
  `k` blocks.  User-supplied `block_sizes` over the grid cap
  are refused before any fit, and a `sac` estimated in another CRS is
  converted to the sweep's units with a warning (a metre range on a US-foot
  axis was drawn 3.3 times too short).  The plot's caption reads every
  `coverage_*` column as "closer to the nominal level is better" (only
  `coverage_50`, `coverage_80` and `coverage_95` were known, as
  higher-is-better, so `coverage_97.5` from `cv_bayes()`'s full-precision
  names was captioned "lower is better").

* `fit_gwr_model()` keeps its local collinearity survey.  Every fitting
  window --- not a sample of 30 --- has scaled condition indices of its
  kernel-weighted local design computed, the way Wheeler and Tiefelsdorf
  diagnose GWR collinearity, and the fit carries them as
  `info$local_collinearity` (one row per observation: coordinates, window
  size, and two condition indices: `cn`, Belsley's uncentred index with the
  intercept, and `cn_slopes`, the predictors centred in the window and
  scaled by their study-area standard deviation), with
  `info$n_local_collinear`, `info$n_local_singular` and the global
  `info$condition_index` (on the centred predictors) beside `AICc`.
  `n_local_collinear` and the warning count the windows whose slopes are
  collinear (`cn_slopes` above 30, or singular), so a predictor's
  origin (degrees C or kelvin) does not change the verdict.  The warning is
  now the exact fraction of locations rather than a sampled one, and its
  thresholds (a quarter of the locations, or any) are unchanged; it says
  that a coefficient an exactly singular window cannot estimate is `NA`.
  The weighted survey sees what the
  unweighted spot-check could not: a bisquare window's edge points
  contribute almost nothing to the fit, so they contribute almost nothing to
  its conditioning.  The survey also runs with a single numeric predictor,
  and its kernel weights equal `GWmodel::gw.weight()` exactly: a boxcar
  keeps a point at the kernel's edge, and a zero-width adaptive kernel gives
  `NaN`, counted as singular, instead of weight 1 at the co-located points.
  (A grid with a boxcar bandwidth equal to its spacing used to be reported
  collinear at every location while GWmodel fitted windows of 3 to 5
  points.)

* `plot(fit, type = "coefficients")` for a GWR fit maps one local
  coefficient (`term`) at the training locations, which is the reason to
  fit GWR at all --- and masks the locations where it is not to be
  believed: a collinear local design (for a slope, the slope condition index
  above 30 or singular; for the intercept, the condition index with the
  intercept above 30) or a non-finite coefficient is drawn hollow and grey,
  counted in the subtitle, because the smooth surface a naive map draws over
  them is the picture of an unstable estimate.  `mask = FALSE` draws them
  anyway and says how many it is drawing.  A fit with a duplicated
  predictor name is drawn rather than refused as "every location is masked".
  A slope map in kelvin is drawn like the one in degrees C; before, every
  location was masked and the map refused.

* `kriging_adequacy()`: what a block-kriging aggregator would deliver on a
  set of cells, computed beside the plain means and changing none of them.
  Per cell, from a fitted variogram (`estimate_sac_range()`'s, or estimated
  here): the block-kriging estimate and variance, that variance as a share
  of the variance the cell's mean would have with no data at all
  (`kr_ratio`, the coverage score --- near 1 the data tell the cell nothing;
  a share of the point sill never came near 1 for cells larger than the
  range), the error variance of the cell's plain mean under the same
  variogram (`plain_var`), the share of it kriging removes (`kr_gain`), and
  the kriged-minus-plain shift over the standard deviation the variogram
  gives that difference (`kr_shift`); plus the variance of the standardised
  errors from blocked cross-validation (`attr(, "cv")$zscore_var`), which
  is about 1 when the kriging variance is right.  Measured on simulated
  exponential fields: 0.93--1.07 with the true variogram, 0.85--1.01 with
  the estimated one; under blocked folds it checks the sill and range
  rather than the nugget (0.95--1.24 with the nugget understated tenfold),
  and random folds are the instrument for the nugget.  The comparison the
  function exists for, against known cell averages on those fields: kriging
  cut the plain mean's squared error by 54--56 percent under uniform and
  clustered sampling alike, where `kr_var` over `plain_var` predicted
  55--58.  A naive `s^2/n` below `plain_var` says the naive standard error
  is too small under the variogram; `kr_var` is not compared with `s^2/n`,
  which moves with the fitted sill.  Repeat visits to one location are
  kriged from their mean, with their pooled within-location variance
  divided by their count; kriging them as separate rows made every kriging system
  singular and returned `NA` everywhere.  Any cell or held-out location
  gstat still cannot solve is counted in a warning.  This is the first
  kriging path in the package (`gstat::krige()`); its model families are
  the ones the package interprets elsewhere, and any other is refused by
  name.  The cross-validation runs each `make_folds()` split on its own
  training set, so `"buffered_loo"` and `"nndm"` keep their exclusion zones
  (reduced to fold labels they ran as plain leave-one-out: with a 250 m
  buffer, a standardised-error variance of 1.02 and RMSE 0.757, against
  0.873 and 0.944 with the buffer kept), and `print()` names the scheme
  instead of calling every one "blocked CV".  A cell holding locations that
  the `nmax` nearest its centre leave out is kriged from all of its own
  locations plus the `nmax` nearest outside it (a cell of 1,500 points was
  otherwise estimated from its middle 50: 0.43 off against 0.07), and the
  column `kr_n_used` says how many locations each cell was kriged from.
  `max_neighbours` (default 2000) leaves out, with a warning, a cell that
  would need a larger kriging system, and `max_box_ratio` (default 1000) a
  cell whose bounding box exceeds its area that many times, because gstat
  discretises over the whole box (+592 MB for one strip at 7,072); both are
  counted in `attr(, "cells_left_out")` and by `print()`.  The points are
  put in the CRS the variogram was fitted in (`attr(sac, "crs")`), as
  `summarize_by_cell()` does: a `sac` fitted in metres used on points in km
  gave a CV statistic of 3.05 against 0.85.  A `sac` whose variogram is of
  residuals (`detrended = TRUE`) is used with a warning that the response is
  kriged without its predictors and the variances come out too small
  (4.3--5.2 against 0.67--1.53).  `attr(, "rejected_reason")` records why a
  `sac`'s range was refused, and the warning and `print()` say it instead of
  "sill never reached" for every refusal.  The cell ID is found and matched
  as `summarize_by_cell()` finds it (cells keyed by `id` or `grid_id` too,
  and a double ID of 1e5 matches an integer 100000, where it used to leave
  that cell with n = 0), a point whose ID matches no cell is counted in a
  warning, and a layer with no CRS is taken to be in the other's.  Folds
  whose recorded rows sit at other locations in `assigned_points_sf` are
  refused, as in `cv_*()`.  Folds built before
  `assign_features_to_polygons()` dropped points used to be applied by
  position, holding out the wrong points (cross-validation RMSE 1.82
  against 2.11 with folds built on the assigned layer), with nothing said.

* `MAPE` and `SMAPE` now say how many rows they were averaged over.  Every
  metrics frame --- `model_metrics()`, `summary()`, `evaluate_insample()`,
  `compare_models()`, and the `overall` and `fold_metrics` of `cv_gwr()`,
  `cv_bayes()`, `cv_spatial()`, `cv_rf()` and `compare_models_cv()` --- gains
  two trailing integer columns, `n_MAPE` and `n_SMAPE`: the rows each
  percentage error actually used once those where its denominator is zero
  were dropped (`y` for MAPE, `|y| + |yhat|` for SMAPE, zero meaning no
  larger than 100 machine epsilons times the data's own magnitude).  They
  equal `n` (`n_pred` in the CV frames) when nothing was dropped, and are
  `0` in an empty frame.  The values themselves are unchanged, apart from
  what now counts as zero (see Bug fixes): a MAPE over 58 of 120 rows is the
  same number 2.0.0 reported, but it now arrives labelled, where before
  nothing in the frame recorded that it was a subset average.
  `print(summary(fit))` appends "(over k of n rows)" to its SMAPE line when
  the two differ.  The columns sit after `Adj_R2` so code addressing the seven
  metric columns by position is unaffected; code pinning the exact column set
  needs the two names added.

* `voronoi_seeds_kmeans()` gains `nstart`, the number of starts (10 for
  `partition = "kmeans"`, as before).

* `cv_bayes()` and `compare_models()` say when a Bayesian fit did not
  converge.  Both scored such a fit like any other, and only the fit's WARN
  log lines (which `tryCatch()` and knitr never see, and
  `spatialkit_quiet()` hides) said that its posterior was not to be
  trusted.  `cv_bayes()`'s `fold_metrics` gains `convergence_ok`: `TRUE` or
  `FALSE` as `fit_bayesian_spatial_model()` judged that fold's sampler
  (divergent transitions, R-hat, effective sample size, E-BFMI), and `NA`
  when `fit_args` sets `check_convergence = FALSE`.  A run with any `FALSE`
  raises one warning of class `spatialkit_convergence` naming those folds
  and the checks each failed.  `compare_models()` gains the same column and
  warns once for each model that did not converge.  Code that pins the
  exact column set of either table needs the new name added.

* `select_resolution()` returns `$seeds`, the centres of the partition the
  profile scored at the chosen level, and `get_voronoi_seeds(method =
  "kmeans")` and `voronoi_seeds_kmeans()` return those centres when given
  the selection or the profile, instead of running a k-means of their own.
  A k-means partition is the Voronoi partition of its centres, so these
  seeds rebuild the very cells the criteria judged.  The fresh k-means the
  seeding functions ran (Hartigan-Wong, 10 starts, on whatever cloud they
  were given) was a different partition: on 300 points it put 5--16 percent
  of them in a different cell from the one they were scored in, and on
  2000 points (profiled on a 1500-point subsample) 26--29 percent, so the
  count was defended on cells nobody built.  The profile keeps every
  level's centres in `attr(x, "centres")` (with `attr(x, "centres_crs")`),
  the seeds come back in the CRS of the boundary or points they are for,
  and `attr(seeds, "kmeans")` assigns `sample_points` to them, each point
  to its nearest seed.  A count passed as a number (`n = sel$best`) still
  runs a fresh k-means, as does a selection made before seeds were kept.
  `voronoi_seeds_kmeans()` did not accept a selection or profile as `k`
  before, and the advice in `build_tessellation()`'s occupancy warning now
  points at the selection rather than the bare count.

## Bug fixes

* **A point layer in an `sfc_GEOMETRY` column stopped the fits, the
  cross-validations, the folds, the cell summaries, the tessellations, the
  seeders and the surface.**  A geometry column of the general
  class `sfc_GEOMETRY` can hold points only: `sf::st_cast(x, "GEOMETRY")`
  builds one, and so do `c()` of such a column with points and a point set
  with `[[<-` in place of the last other geometry of a mixed column.
  `sf::st_coordinates()` refuses that class, so `fit_gwr_model()`,
  `fit_rf_model()`, `fit_bayesian_spatial_model()`,
  `gwr_model_selection()`, `prep_model_data()`, `predict()` of a fit,
  `model_metrics()`, `evaluate_insample()`, `compare_models()`,
  `cv_spatial()`, `cv_gwr()`, `cv_rf()`, `compare_models_cv()`,
  `cv_block_size_sweep()`, `area_of_applicability()`,
  `summarize_by_cell()`, `cell_values_at()` (`newdata`, `points` and
  `seeds`), `determine_optimal_levels()`, `resolution_profile()` (and so
  `select_resolution()`), `create_voronoi_polygons()`,
  `build_tessellation(method = "triangles")`,
  `get_voronoi_seeds(method = "kmeans")` with `sample_points`,
  `voronoi_seeds_kmeans()` and `predict_surface()` (`grid` and
  `covariates`) stopped with sf's "not implemented for objects of class
  sfc_GEOMETRY".  On the column `sf::st_cast()` returns,
  `sf::st_dimension()`, `sf::st_union()` and `sf::st_transform()` stop with
  "Not compatible with STRSXP": `make_folds()` stopped with it, as
  `points_sf` and as `prediction_points`, and so did `fold_separation()`,
  `ensure_projected()`, `harmonize_crs()`, `assign_features_to_polygons()`
  and `create_voronoi_polygons()`, and every call that moves such a layer
  to another CRS, such as `model_metrics()` with it as `newdata` in another
  CRS than the fit.  `plot_tessellation_map()` (`features_sf`, `seeds_sf`
  or `boundary`), `plot_folds()` and `plot()` of a fit whose layer is such
  a column returned a plot that stopped with that message when it was
  built or printed, and `create_grid_polygons_cached()` stopped with it on
  such points as `boundary`, where `create_grid_polygons()` refuses a
  boundary that is not polygonal.
  `select_features_forward()` reported that stop as a failure to build its
  inner folds, and on the other form of the column it scored no candidate
  and returned an empty selection.  When every geometry of such a column
  is a POINT, these functions now read it as the same points in
  `sfc_POINT` and give the answer those points give, and when every
  geometry is a MULTIPOINT, as the same features in `sfc_MULTIPOINT`.  A
  layer moved to another CRS comes back in that class, and so do the rows
  `assign_features_to_polygons()` returns; `area_of_applicability()` still
  returns the rows of `newdata`, and `predict_surface()` a `grid` it does
  not move, in the class they came in.  A layer that
  mixes points with other geometry types is read as before.
  `build_tessellation(method = "triangles")` also stopped with sf's message
  on a layer that mixes POINT and MULTIPOINT features, which sf holds in an
  `sfc_GEOMETRY` column; it now reads their vertices as it reads a
  MULTIPOINT layer's.

* **`resolution_profile()` with a variogram stopped on a study area much
  longer than wide.**  The mean correlation over the study area, which
  `reliability` reads, came from a grid over the hull's bounding box at a
  spacing that follows the hull's area.  On a strip along an axis narrower
  than half that spacing (300 points over 20000 x 10 m, in metres or in
  kilometres), the call stopped with "wrong sign in 'by' argument" whenever
  it held a fitted variogram (`sac` from `estimate_sac_range()`, or a
  ladder set with `levels` and `sac = "auto"`).  On a thin strip at a
  slant, the grid's rows crossed part of its length, or fewer than 10 nodes
  fell inside and the bounding box's value was used: on the strip at 30
  degrees the mean correlation was 0.058 where the strip's own is about
  0.096.  A hull narrower than two steps of the grid is now laid out in
  the frame of its rectangle of least area, in bands along it with the
  nodes spread over each band's width, so a hull that tapers (a wedge, a
  lens) is read over its area, and its value does not move with its angle.
  Every other hull keeps its grid and its results.

* **Areas in a CRS whose unit sf has no name for were read in square
  metres.**  For a projected CRS whose `+to_meter` unit sf has no name for
  (a unit of 1000 km, of 3.28 m or of 0.5 mm), `sf::st_area()` converts to
  square metres.  `summarize_by_cell(area = TRUE)` returned that as
  `cell_area`, which its help gives in the CRS's unit squared, so a square
  of 1000 m in units of 1000 km got 1e6 where its area is 1e-6, and
  `n_per_area` was a count per square metre.  `resolution_profile()`
  divided the hull's area in square metres by the square of a range in CRS
  units, so a range floor of 3 came out as 31 in units of 3.28 m, and in
  units of 1000 km the floor was not finite or the profile stopped with
  "wrong sign in 'by' argument".  Three checks read such an area beside
  lengths in the CRS's unit as well.  `kriging_adequacy()` left out every
  square cell as a sliver in units of 0.5 mm (a box 4e6 times its area)
  and let a real sliver through in units of 1000 km.  `make_folds()` with
  supplied `blocks` read their scale in metres against a range in CRS
  units, so in units of 1000 km it gave no warning for blocks smaller
  than the range, and in units of 0.5 mm it quoted the wrong scale.
  `create_grid_polygons()` chose the edge strips to drop after its clip by
  comparing such an area with a bound in the CRS's unit squared.  All of
  them now take the planar area in the CRS's unit squared.  Layers in
  metres, in a unit sf names (such as US survey feet or kilometres) and in
  lon/lat, whose cells keep their geodesic square metres, give the same
  results as before.

* **Line midpoints in a CRS whose unit sf has no name for were off the
  middle.**  `coerce_to_points()` modes `"auto"` and `"line_midpoint"`, and
  so `prep_model_data()` and the fits that pointize lines, took a line's
  midpoint with `sf::st_line_sample()`, which measures the line in metres
  for such a CRS and reads half that length as a distance along the line
  in the CRS's unit.  A straight line from x = 0 to x = 20 got its point at
  x = 20, its end, in a unit of 2 m, and at x = 5 in a unit of 0.5 m; a
  MULTILINESTRING's longest part the same.  The midpoint is now taken at
  half the planar length in the CRS's own unit, at x = 10 in every unit.
  Lines in metres, in a unit sf names (such as US survey feet or
  kilometres) and in lon/lat with `tmp_project = TRUE` get the same points
  as before.

* **Voronoi cells with no boundary changed with the unit of the
  coordinates.**  `create_voronoi_polygons()` and
  `build_tessellation(method = "voronoi")` clip the cells to the points'
  convex hull, or to the cover, grown by 2 percent of its diagonal, as
  their help says; the margin was floored at 1 unit of the CRS.  A layer
  whose diagonal is under 1 unit, such as a 500 m layer in a kilometre CRS,
  grew by 0.02 units (20 m) where in metres it grows by 15.9 m, and its
  edge cells covered 2.9 percent more ground.  The margin is now 2 percent
  of the diagonal in every unit.  A layer whose diagonal is 1 unit or more
  gets the same cells as before, and one location still grows by 0.02
  units.

* **`clip_target_for()` returned disks around points whose box it should
  return, depending on where they lie.**  With no `boundary`, the points'
  bounding box is the target, and a buffer around the points only when
  the short side of the box is under a millionth of the long side.  A
  side whose two ends `all.equal()` called equal also took the buffer, and
  that test compares the ends relative to their size, so it read the
  origin: a 3 km x 2 km box at x = 1e12 got three disks of 150 m (an area
  of 211,961 for 6,000,000), and a 5 cm x 3 cm box at y = 4,000,000 got
  disks too.  The sides are now read against the span only, so the box
  comes back at any offset.

* **`clip_target_for()` stopped with R's "missing value where TRUE/FALSE
  needed" on points none of which has a location.**  With no `boundary`,
  points that are all empty, or all have a coordinate that is not finite or
  lies beyond 1e150, or that the move of lon/lat points into a projected
  CRS empties (a latitude beyond 90), have no bounding box.  The call now
  refuses them with an error that says no point has a location.  With a
  `boundary`, such points still get the boundary as their target.

* **`clip_target_for()` said it used a "small buffer" around points at one
  location, and gave neither its size nor its unit.**  With no `boundary`
  and every point at one location the bounding box has no extent to take a
  share of, and the target is a disk of radius 1 unit of the CRS: 1 m in a
  metre CRS and 1 km in a kilometre CRS.  The target is unchanged.  The
  message now gives the radius in the layer's unit (in the units of the
  coordinates when the points have no CRS) and says to pass `boundary` for
  a target of real size, or an `expand` of 1 or more, which is taken as the
  radius; `?clip_target_for` says the same.

* **Block folds, the join to polygons, the Voronoi cells and the fill of
  `cell_values_at()` changed with the origin of the coordinates.**  A point
  counts as on an edge while it lies within a margin of it, one part in
  1e10 of the largest coordinate: 0.44 mm at a UTM northing of 4,400 km,
  but 100 m on a 4 km layer shifted to 1e12 and 100 km at 1e15.  At 1e12
  points within 100 m of a block edge went to the lower block, so
  `make_folds()` and the functions that make folds with it
  (`cv_spatial()` and `compare_models_cv()` without `folds`,
  `fold_separation()`, `cv_block_size_sweep()`, `kriging_adequacy()`,
  `area_of_applicability()` and `select_features_forward()`) gave other
  folds than at the origin; `assign_features_to_polygons()` gave points
  tens of metres inside a polygon to its neighbour; the 60 points of a
  layer fell in 39 Voronoi cells of `create_voronoi_polygons()` and
  `build_tessellation(method = "voronoi")`; and `cell_values_at()` gave
  other IDs.  At 1e15 `make_folds()` refused with a false "single block"
  and the join called every point a tie.  The margin is now at most a
  millionth of the diagonal of the coordinates' bounding box, and at
  least four times the spacing of doubles at the largest coordinate, so
  at any offset where doubles hold the layer these results are the ones
  near the origin for every point farther than that from an edge.  Where
  the diagonal is at least a ten-thousandth of the largest coordinate
  (400 m at a UTM northing of 4,000 km) the margin is unchanged; on a
  smaller layer it is the millionth of its diagonal (28 micrometres on a
  layer 20 m across).  The help pages that state the margin say so.

* **Coordinates so far from their origin that doubles cannot hold the
  layer were read without a word.**  At 1e15 doubles are 0.125 apart, more
  than a millionth of a layer a few kilometres across, and the edges
  `make_folds()` lays for blocks and `create_grid_polygons()` lays for its
  grid are rounded to that spacing, so a point near an edge could fall in
  the next block with nothing said, and two points 0.25 m apart fell in one
  Voronoi cell.  `make_folds()` (`"block_kfold"` and `"buffered_loo"`),
  `assign_features_to_polygons()`, `create_grid_polygons()`,
  `create_voronoi_polygons()` and `build_tessellation()` (`"voronoi"` and
  `"triangles"`; `"hex"` and `"square"` through the grid) now warn once, with class
  `spatialkit_xy_imprecise`, where the spacing of doubles at the
  coordinates is more than a millionth of the layer's extent; the warning
  gives both and says to shift the layer nearer the origin or project it.
  A function that makes folds several times on one layer says it once.

* **`clip_target_for()` drew its buffer around points on one line or at
  one location at the rounding of their coordinates, without a word.**
  With no `boundary`, GEOS draws the buffer at the size of the
  coordinates: at an offset of 1e15, where doubles are 0.125 apart, the
  1-metre disk around one location had an area of 2.6875 for pi, and two
  points 0.5 m apart on a line got an empty LINESTRING for a target.  The
  function now warns, with class `spatialkit_xy_imprecise`, where the
  spacing of doubles at the points is more than a millionth of the
  target's extent, and refuses, with an error of class
  `spatialkit_xy_imprecise_error`, where the buffer comes back empty or
  with no area; both give the spacing and the radius and say to shift the
  layer nearer the origin or pass `boundary`.  `quiet = TRUE` does not
  silence the warning.

* **`clip_target_for()` returned a point, a MULTIPOINT, an empty line or a
  disk of subnormal area for points less than about 3e-153 apart.**  With
  no `boundary`, points on one line were buffered by 5 percent of their
  spread, and a radius whose square is below the smallest normal double
  (about 1.5e-154) gave GEOS a disk whose area underflowed: points 1e-320,
  1e-310 and 1e-160 apart got a MULTIPOINT, an empty LINESTRING and a
  disk of area 1.6e-322, with the message "using small buffer around
  points"; a box 1e-320 on a side came back as a POINT with no message.
  A side of the box, or a share of the spread, below that length now
  counts as none, so such points get the disk of radius 1 unit of the CRS
  that points at one location get, and the message gives their spread.

* **`summarize_by_cell(area = TRUE)` stopped on cells that carry M.**  Cells
  that are XYM or XYZM throughout reached GEOS, which does not read M, and
  the call stopped with GEOS's "GEOS does not support XYM or XYZM
  geometries".  The cells are now checked and measured in x and y, as cells
  that mix dimensions were, and are returned with their own geometry.

* **`estimate_sac_range()` drew its `n_max` subsample over rows with nothing
  to model.**  300 valued rows among 20,300 were estimated from the 74 that
  fell in the draw, and 100 among 20,100 came back `NA`.  Rows with a missing
  or infinite response or predictor are now dropped first, and the subsample
  no longer depends on the order of the rows.

* **`estimate_sac_range()` counted rows, not locations, against its floor of
  30.**  Rows that share a location are one place to the variogram and to
  its test for spatial structure, yet five locations visited 20 times passed
  the floor, got a variogram of ten pairs of locations, and were judged by a
  test of five location means, which called most such layers unstructured
  even where the range was three times the extent.  The floor now counts
  distinct locations as well, among the rows the variogram is computed on.
  Below 30 the result is the unclassed `NA` whose `rejected_reason` gives
  the count ("5 distinct location(s) among the 100 point(s) with a finite
  value to model, fewer than the 30 a variogram range is estimated from"),
  and nothing is fitted.  When the layer has 30 or more locations and the
  `n_max` subsample keeps fewer, the reason gives the subsample's counts
  and says "of the n_max subsample".  The functions that estimate a
  variogram read it as they read the other early refusals:
  `make_folds(auto_range = TRUE)` builds geometric blocks and names it in
  its warning, `summarize_by_cell(deff = "variogram")` falls back to Kish
  with `spatialkit_deff_fallback` where the verdict of no structure had
  given a design effect of 1 in every cell, `kriging_adequacy()` stops with
  it, and the Bayesian fit's check before fitting reports the range as not
  estimated.  With `detrend = "reml"` such a layer got a warning that the
  REML fit did not converge when nothing had been fitted; it gets the same
  `NA` and no warning.  Where a predictor is missing on every row of such a
  layer, there is no warning that the variogram is fitted to the raw
  response, since none is fitted.  A layer whose rows all lie at one
  location gets this reason in place of "the points have no extent".

* **`estimate_sac_range(detrend = "reml")` kept the first fit to converge,
  not the best.**  When the true range is short against the extent the
  starts of `nlme::gls()` lay in a valley of the likelihood: a layer with a
  10 km range came back at 502 km and was refused.  The start is now the
  best point of the profile likelihood on a grid, and the fit with the
  largest log-likelihood is kept (`attr(, "reml")` gains `logLik`,
  `logLik_gain` and `n_fits`).

* **`estimate_sac_range()` returned `NA` or an R error for a bad argument
  without naming it.**  `cutoff`, `range_frac`, `n_max`, `seed`,
  `keep_directional_fits`, `n_lags` and `structure_alpha` are checked, and a
  bad value is an error that names the argument.  `cutoff` above 1 and
  `n_max` below 30 are errors; `n_max = Inf` means no subsample.

* **`estimate_sac_range()` called a response constant by an absolute
  floor.**  A response whose standard deviation was under 1.5e-8 in its own
  unit was refused, and so was one of standard deviation 1 on a mean of
  1e8.  The test is now relative: under 1e-12 of the largest absolute value.

* **The projection chosen for lon/lat input depended on the order of the
  rows.**  `ensure_projected()` scores candidate projections on a sample of
  the points, taken at even spacing along the rows: 400 points over the
  conterminous United States got an Albers CRS in some row orders and a
  Lambert azimuthal one in others.  The sample is now taken along the points
  sorted by coordinate.

* **The projection chosen for a few detailed outlines, and the area check
  on more than 200 polygons, depended on the order of the features.**  A
  layer of fewer than 40 polygons or lines is scored on its outlines'
  vertices too, and more than 160 vertices were first thinned at even
  spacing along the vertices as the layer lists them: five outlines of 200
  vertices over the conterminous United States got an Albers CRS in three
  of eight feature orders and a Lambert azimuthal one in the other five.
  The area check of `summarize_by_cell(area = TRUE)`, and of
  `ensure_projected(purpose = "area")` on projected input, measures at most
  200 of a layer's polygons, and they were taken at even spacing along the
  rows: 300 cells in a UTM zone, one of them 700 km east of its central
  meridian, were refused in five of eight row orders and given a density in
  the other three, and `attr(, "area_error")` changed with the order.  The
  vertices are now thinned along the outlines' distinct vertices sorted by
  coordinate, as the points of a point layer are sampled: a vertex listed
  twice (a ring's closing vertex, a border two features share) counts once,
  so a ring listed from another vertex is scored the same, and a thinned
  outline is not densified between the vertices kept.  The polygons are
  taken along the polygons sorted by their bounding boxes, and polygons
  that share a box, as the Delaunay triangles of gridded points do, by the
  means of their vertices; the polygons that reach furthest along x and y,
  either way, are always measured.  For a given layer the CRS, the area
  error and the refusal can differ from before.

* **A small lon/lat layer whose features share edges, or meet another
  feature's edge at a corner, got its projection chosen on rounding
  noise.**  When a polygon or line layer has fewer than 40 distinct
  feature points and fewer than 40 distinct vertices, `ensure_projected()`
  scores the candidates on points added along each edge of its outlines as
  well.
  An edge two features share was filled once from each end, some copies of
  an added point differed in their last bits, and the pair, about 1e-9 m
  apart, set the worst-case distance error of every candidate.  Two boxes of
  12 by 15 degrees, one above the other, scored 76.9 percent in UTM zone 13
  and 86.8 percent in EPSG:5070, where their outline as one box scores 0.5
  and 10.9 percent, and the zone was kept with "76.92% vs 76.92%".  A 2 by 2
  grid of boxes moved to a Lambert azimuthal CRS with "11.22% against the
  zone's 13.50%", where its outline gives 0.06% and 0.26%.  The same held
  where a point added along one feature's edge falls on a corner of another
  (a T-junction, or a corner that touches the edge) or on a point added
  along an edge that overlaps it in part.  Each edge is now filled from the
  same end whichever feature lists it, and an added point within about 1e-9
  degrees (0.1 mm) of a vertex or of another added point is left out; the
  vertices are kept as they are.  The two layers above score within a tenth
  of a percentage point of their outlines, and the log quotes the outlines'
  figures.  A layer on which no added point lies within 1e-9 degrees of
  another point scores as before, but for the last bits.  For the others
  the CRS, `attr(, "crs_choice")` and the logged figures can differ from
  before, and the points added along an edge no longer depend on the
  direction in which its ring is listed.

* **The projection chosen for a lon/lat polygon or line layer depended on
  how many rows repeat one feature.**  `ensure_projected()` scores the
  candidates on one point per row, plus the outline's vertices when there
  are fewer than 40 such points, and it thinned an outline of more than 160
  vertex rows without filling its edges.  A layer that lists one polygon
  once per year was scored by the number of years: a box of 52 by 22
  degrees over the conterminous United States got a Lambert azimuthal CRS
  with "1.76% against the zone's 10.87%" on one row, "1.28% against the
  zone's 8.33%" on 33 rows, and from 40 rows every candidate scored `NA`
  and UTM zone 15 was kept "because the candidate projections could not be
  scored on this layer".  Copies of one point now count once, before the
  count of 40 and in the sample of 40 points, and an outline is thinned
  only when it has more than 160 distinct vertices.  A layer whose features
  lie on distinct points, and whose outline lists more than 160 vertex rows
  only when it has 40 distinct vertices or more, scores as before.  For
  the others the CRS, `attr(, "crs_choice")` and the logged figures can
  differ from before.

* **The projection chosen for a lon/lat point layer depended on how many
  rows list one location.**  `ensure_projected()` scores the candidates on
  40 points taken at even spacing along the rows sorted by coordinate, and
  repeat visits to one station took places in that sample from the other
  stations.  39 stations scored 5.07 percent in UTM zone 14 and 0.61
  percent in the Albers CRS chosen; with station 1 listed 100 times, 4.59
  and 0.45 percent, and listed 1,001 times, 2.32 and 0.35 percent, the
  sample holding 4 of the 39 stations.  The sample is now taken over the
  distinct locations.  A layer with no location listed twice scores as
  before; for the others the CRS, `attr(, "crs_choice")` and the logged
  figures can differ from before.

* **The area check of `summarize_by_cell(area = TRUE)` and of
  `ensure_projected(purpose = "area")` on projected input depended on how
  many rows repeat one polygon.**  The check takes the spread of the
  ratios of planar to geodesic area over the layer's polygons, and it
  counted rows: copies of a polygon weighted the median by their number,
  took places in the sample of 200 polygons, and two copies of one polygon
  passed as two probes.  Three 50 km squares in UTM zone 14 measured 0.93
  percent, and 1.04 percent with the left one listed three times, so
  `ensure_projected()` warned and `summarize_by_cell()` refused the cells
  only then; 150 squares of 1 degree in UTM zone 15 measured 24.98 percent,
  and 30.79 percent with one of them listed 300 more times; a box in Web
  Mercator listed twice measured 0, where alone it is unmeasured.  The
  check now measures each distinct polygon once, the first of its exact
  copies kept, and a layer whose polygons are all copies of one is
  measured on the grid over its bounding box, as one polygon is.  A layer with no polygon listed
  twice gets the same figure as before.

* **Two listings of one place a few bits apart, or at longitudes 180 and
  -180, set the distance error of every candidate CRS.**
  `ensure_projected()` scores each candidate on the worst ratio of planar
  to geodesic length over pairs of sampled points, and it took every pair
  longer than 0.  A corner two features share, written `-107.2` in one and
  `-110.6 + 3.4` in the other, a point listed twice a few bits apart, and
  a layer whose features meet at longitude 180, which lists (180, y) and
  (-180, y), each gave a pair about 1e-9 m long whose planar length is
  rounding noise, or the width of the world in Web Mercator.  Two boxes
  side by side, typed that way, scored 13.7 percent in UTM zone 13 and
  24.5 percent in EPSG:5070, where their outline as one box scores 0.26
  and 0.97 percent, and got an Albers CRS.  Two boxes either side of 180
  reported a distance error of 1 for their one candidate.  Pairs under 1
  mm are now left out: the two boxes score 0.25 and 0.97 percent and get
  their outline's Lambert azimuthal CRS, and the boxes either side of 180
  report 0.35 percent.  A layer with no pair 1 mm or longer is scored `NA`,
  as a layer of one location is.  A layer whose sampled pairs are all 1 mm
  or longer scores as before; for the others the CRS, `attr(,
  "crs_choice")` and the logged figures can differ from before.

* **The geodesic distances that score candidate projections came out
  short on east-west pairs of a few millimetres.**  `ensure_projected()`
  compares planar distances with WGS84 ellipsoidal ones from Vincenty's
  formula, whose iteration stopped once a pass moved its longitude on the
  auxiliary sphere by less than 1e-12 radians.  On an east-west pair under
  about 2 mm at the equator, 2.5 mm at latitude 40 and 11 mm at latitude
  80, the first pass already moved it by less, and the distance came out
  up to 0.34 percent short; pairs up to about a metre long stopped after
  the second pass, up to 1.1e-5 short.  The shortest pair can set the
  worst-case distance error of every candidate: 30 points between
  latitudes 33 and 41 score 0.04 percent in UTM zone 13, and scored 0.20
  percent with a 31st point 1.8 mm east of one of them.  The iteration now
  stops once a pass moves the longitude by no more than 1e-12 of the
  pair's longitude difference, taken the short way round (of 1 radian for
  pairs 1 radian or more apart), and such pairs agree with the ellipsoid
  to 1e-9, across longitude 180 too.  Nearly antipodal pairs still fall
  back to the spherical distance.  Pairs a degree or more long move by
  less than 1e-10 of their length, and other pairs by up to 1.1e-5 of
  theirs, the error the old test left on them; the distance errors in
  `attr(, "crs_choice")` and the log move by about as much, and the CRS
  chosen changes only where two candidates score that close.

* **One row of a lon/lat layer at a coordinate that is not finite made
  every candidate projection score `NA`.**  The geodesic distances stopped
  with "NAs are not allowed in subscripted assignments" on the pairs of
  that row, the scorer turned the stop into `NA`, and 39 stations over the
  conterminous United States with a 40th row at longitude `Inf` were put
  in UTM zone 14 with an `NA` in `attr(, "crs_choice")`, where the 39 get
  an Albers CRS.  Such a row now takes no part in the scores, as an empty
  row does, and the zone's figure is the stations'.

* **`ensure_projected()` kept a UTM zone that cannot place every location
  of a lon/lat layer, and `coerce_to_points()` then ended the R session.**
  PROJ gives no coordinates for a point a UTM zone cannot place, about 90
  degrees or more from its central meridian on the equator, and the
  distance score left such points out, so the zone was scored on the rest.
  `LINESTRING (0 0, 180 0)` kept UTM zone 46, which holds neither end of
  the line (278.82% against 1169.31% for the Lambert azimuthal CRS);
  `ensure_projected()` returned the line empty, and five points on the
  equator at longitudes 0 to 180 with two of them empty, without a word;
  and `coerce_to_points()` in mode `"auto"` (the default, also of
  `prep_model_data()` and of the functions that reduce a layer to points)
  or `"line_midpoint"` handed the empty line to `sf::st_line_sample()`,
  which segfaulted.  A candidate is now used only if it places every
  located vertex of the layer (a vertex at a latitude beyond 90, which no
  CRS places, takes no part in this test): one that does not gets
  `distance_error` `NA` in `attr(, "crs_choice")`, is not chosen and is
  named in the log line, and when no candidate places the layer the global
  fallback (Web Mercator, Equal Earth for `purpose = "area"`) is used.  A
  sampled point that a candidate cannot place makes its score `NA`.  The
  line and the five points now get the Lambert azimuthal CRS, and the
  line's midpoint is (90, 0).  When a move out of lon/lat still leaves a located row empty,
  or without some of its vertices (`sf::st_transform()` drops a vertex it
  cannot move), `ensure_projected()` warns with class
  `spatialkit_rows_unplaced`, naming the rows.  `coerce_to_points()` never
  samples a line that its temporary projection did not place in full: it
  stops with class `spatialkit_features_degenerate` and
  `site = "projection"`, naming the rows and offering `mode = "centroid"`
  or `tmp_project = FALSE`, and the functions that reduce a layer to
  points raise it in their own words.  A lon/lat line with a coordinate
  that is not finite or lies beyond 1e150, which crashed the same way,
  gets an empty point, as an empty line does.  Layers on which every
  candidate places each vertex and each point the score samples get the
  CRS, the scores and the points they got.  A candidate that places every
  vertex but not a point the score samples along an edge now scores `NA`:
  on `LINESTRING (0 0, 180 0)` with `purpose = "area"`, the Lambert
  azimuthal CRS, which is still used.

* **`harmonize_crs()` and the functions that move a layer into the CRS
  they work in returned a row the target CRS cannot place in full, emptied
  or short of a vertex, without a word.**  GDAL returns a point it cannot
  move as an empty point and drops such a vertex from a line or a ring.
  `harmonize_crs(a, b, prefer = "b")` with `b` in UTM zone 46 returned the
  point (10, 0.1) of `a` empty and `LINESTRING (10 0, 90 0, 170 0)`
  without its vertex at longitude 10, and a ring with a vertex at latitude
  95 lost that vertex, through `harmonize_crs()` and through the
  `boundary` of `predict_surface()`; `ensure_projected(target_crs =)` did
  the same in silence from a projected layer, such as one in EPSG:3857.
  Every such move now warns with class `spatialkit_rows_unplaced`, naming
  the rows, the function called and the argument: `harmonize_crs()` (`a`
  or `b`, not after the `on_transform_error = "set_crs"` stamp),
  `ensure_projected(target_crs =)` from any CRS, the seeders and the
  tessellation builders (`boundary`, `seeds`, `sample_points`, `cover`,
  and the centres of a profile passed as `n` or `k`),
  `assign_features_to_polygons()` (`features_sf`, `polygons_sf`),
  `kriging_adequacy()`, `cell_values_at()` (`values`, and `newdata` where
  it is placed in the polygons of `values`), and the functions that move a
  layer they are given into the CRS of their data: `make_folds()`,
  `cv_spatial()`, `cv_gwr()`, `cv_bayes()`, `area_of_applicability()`,
  `predict_surface()`, `fold_separation()`, `resolution_profile()`,
  `summarize_by_cell()`, `create_grid_polygons()`, `plot_folds()`, the
  `predict()` methods of the GWR, random forest and Bayesian fits, and
  `plot_tessellation_map()` for every layer it moves into the plot's CRS,
  with a CRS of its own or without one.  So does a row of
  a projected layer with a coordinate beyond 1e10 that a move to a CRS in
  degrees empties.  The warning names the CRS the layer was in ("a
  location in EPSG:3857"), and longitude and latitude for a lon/lat layer
  or for one without a CRS that is read as lon/lat.
  `cell_values_at()` and the pooled readers of `summarize_by_cell()` still
  report the locations and seeds they place in their own words.  The rows
  come back as they did; layers the target CRS places in full raise no new
  warning and give the same results.

* **`plot_tessellation_map()` and `plot_folds()` drew rows without a
  location.**  A cell, boundary, seed or feature with a coordinate at `Inf`
  gave the map an infinite x range in silence, and one at 1e200 (beyond
  the package's 1e150 bound) a map 1e200 wide where every other row shrank
  to a point, with ggplot2's bare "Position guide is perpendicular to the
  intended axis".  `plot_folds()` drew the point `make_folds()` had just
  dropped as without a location, at that coordinate.  Such a row is now
  left out of its layer, read as given before any move into the plot's
  CRS, with one warning naming the argument and the count ("not drawing 1
  row(s) of `features_sf` without a location").  `plot_folds()` also
  leaves out a `boundary` row with such a coordinate, and a point with no
  fold that `make_folds()` left out as empty or because the CRS the folds
  were built in cannot place it (a latitude beyond 90), which it drew at
  that latitude.  A move into another CRS emptied the same rows before,
  so the map no longer depends on whether a move happens.  Layers without
  such rows are drawn as before.

* **`plot_tessellation_map(labels = TRUE)` handed ggplot2 a label for a
  cell with an NA ID or an empty geometry**, and ggplot2 warned at print
  "Removed 1 row containing missing values or values outside the scale
  range (`geom_text()`)" (12 rows when every ID was NA).  A cell is now
  labelled when it has a label and a location; when no cell has both, the
  labels are skipped with a log line, as they are when there is no ID
  column.

* **A layer without a CRS whose coordinates cannot be longitude and
  latitude was stamped with a CRS in longitude and latitude.**
  `plot_tessellation_map(target_crs = 4326)` on CRS-less cells in metres
  found that they do not look like longitude and latitude, then stamped
  EPSG:4326 on them all the same: the map was drawn at "latitude"
  4,000,000, and sf added its bare "bounding box has potentially an
  invalid value range for longlat data".  A layer without a CRS is now
  refused a CRS in longitude and latitude when it has an x outside -180
  to 360 or a y beyond 90 north or south.  The refusal names the
  argument, its coordinates' box and the target CRS, and asks for the
  layer's CRS (or, in `plot_tessellation_map()`, `harmonize_crs(target_crs
  =)` and `ensure_projected(target_crs =)`, a projected `target_crs`).  It
  applies where a layer without a CRS is brought into the CRS of another
  layer or argument: `harmonize_crs()`,
  `ensure_projected(target_crs =)`, `assign_features_to_polygons()` (which
  names `features_sf` or `polygons_sf`; features in metres beside cells in
  degrees fell in no polygon, and polygons in metres stopped on an s2
  refusal), `plot_tessellation_map()`, `plot_folds()`, `make_folds()` and
  the `cv_*()` functions, `area_of_applicability()`, `predict_surface()`,
  `fold_separation()`, `kriging_adequacy()`, `resolution_profile()`,
  `summarize_by_cell()`, `cell_values_at()`,
  `create_voronoi_polygons()`, `create_grid_polygons()` and
  `build_tessellation()`.  It applies as well where a fit's reading of its
  CRS-less training data as longitude and latitude is given to CRS-less
  new data: `predict()` on a random forest, GWR or Bayesian fit, and the
  `grid`, `boundary` and `covariates` of `predict_surface()`.  Every
  prediction came back `NA` with a warning that the rows have a location
  in longitude and latitude, and such a `boundary` held no grid point.
  So does `cell_values_at()`'s reading of
  `newdata` or `seeds` without a CRS beside `points` without one in
  degrees, where every row came back `NA` with no word about `newdata`.
  In `create_grid_polygons(crs =)` and `build_tessellation(crs =)` a
  `boundary` with a y between 90 and 180 is still stamped, and refused by
  the row's latitude as before.  A target CRS that is projected is still
  stamped, with a warning (with a log line in
  `ensure_projected(target_crs =)`).

* **A row of a lon/lat layer with a vertex at a latitude beyond 90 shaped
  the projection `ensure_projected()` chose for the other rows.**  No CRS
  places such a row, but the choice read it in the extent, the longitudes
  of the antimeridian test and the centre.  39 stations over the
  conterminous United States get an Albers CRS (standard parallels 32.9
  and 43.0); with a 40th row at (-96, 95) they got a Lambert azimuthal CRS
  centred on (-96.7, 40.7), and the Albers candidate's parallels moved to
  41.2 and 84.2.  A row at (-100, -95), or a swapped row at (40, -96), did
  the same.  Such a row now takes no part in the choice, the whole row, as
  a row with a coordinate beyond 1e150: the stations get the CRS and the
  `crs_choice` they get alone.  The row is still moved, comes back empty or
  without the vertices the CRS cannot place, and the
  `spatialkit_rows_unplaced` warning names it.  A layer each row of which
  has such a vertex gets the global fallback, and the log line says so.
  Layers with no latitude beyond 90 get the CRS they got.

* **The readers of a lon/lat data layer read a row with a latitude beyond
  90 as given.**  sf warns "bounding box has potentially an invalid value
  range for longlat data" each time it reads the box of such a layer, and
  that warning reached the caller in sf's words from `coerce_to_points()`,
  `prep_model_data()` and the fits on polygons, the k-means seeders
  (`get_voronoi_seeds(method = "kmeans", sample_points =)`,
  `voronoi_seeds_kmeans()`) and `assign_features_to_polygons()` on its
  features.  Where the row lasted to the end it was read where it was
  written: a square at 35 N with one vertex at 95 got a box centre at
  65 N under `pointize = "bbox_center"`, about 3,300 km from its other
  vertices, and `ensure_stable_poly_id()` keyed it and returned it after
  s2's repair, without a word.  Such a row of a layer in a geographic CRS
  now has no location in these readers, as a row with a coordinate beyond
  1e150: `coerce_to_points()` gives it an EMPTY POINT in every mode,
  `prep_model_data()` and the fits drop it and count it, the seeders drop
  it from the cloud, `assign_features_to_polygons()` leaves it unassigned,
  and `summarize_by_cell()` and `cell_values_at()` leave it out of the
  pooled and scored rows, each in the words they use for a row without a
  location, which now name "a latitude beyond 90 degrees" for a lon/lat
  layer.  `coerce_to_points()` and `ensure_stable_poly_id()` say so with a
  warning of class `spatialkit_lat_beyond_90` that names the rows;
  `ensure_stable_poly_id()` keys such a row as an empty polygon, after the
  polygons with a location, and returns it as it came.  A point layer with
  such a row no longer gets the projection's warning that it cannot place
  the row beside the drop line.  `ensure_stable_poly_id(method =
  "centroid")` and `method = "surface_point"` on a layer whose every
  polygon is empty stopped with base R's "missing value where TRUE/FALSE
  needed"; they now give each row an ID.  Layers with no latitude beyond 90 give the same results.

* **`ensure_projected()` stopped with "invalid crs" on a lon/lat layer
  whose box has no height and whose centre lies within 5 degrees of a
  pole.**  On such a box the Albers candidate's standard parallels are 5
  degrees either side of the centre on the sphere, and a line along one
  parallel that spans 180 degrees of longitude, such as `LINESTRING (0 45,
  180 45)`, runs over the pole, so its centre is the pole and the
  parallels were 85 and 95, which PROJ refuses.  With `purpose = "area"`
  every such layer stopped, a single point included; with `purpose =
  "distance"` such a layer stopped where it reaches more than 5 degrees
  from its UTM zone's central meridian.  Every function that moves such a
  layer through `ensure_projected()` stopped with it.  The Albers
  candidate is now built only when both standard parallels lie between
  -90 and 90; otherwise it is left out of the comparison, as a degenerate
  cone is.  The line gets the Lambert azimuthal CRS centred on the pole,
  and a point at (10, 89) with `purpose = "area"` gets the one centred on
  it.  Layers whose parallels lie in that range get the CRS and the
  `crs_choice` they got.

* **`ensure_projected()`'s warning on a lon/lat layer across a pole, or
  one whose longitudes are given outside -180 to 180, gave an offset from
  the UTM zone's central meridian beyond 180 degrees.**  The zone comes
  from the centre on the sphere, whose longitude can lie far outside the
  box when the layer runs across a pole or its longitudes are given
  outside -180 to 180, and the offset was the difference of the box's ends
  from the meridian, not wrapped.  For two lines from (-10, 70) to
  (170, 70) and from (-10, 75) to (170, 75) the warning said the extent
  "reaches 257.0 deg from the central meridian of UTM zone 16", and for
  points at longitudes 170, 185 and 200 it said "reaches 377.0 deg from
  the central meridian of UTM zone 1".  The offset is now the largest
  distance on the circle from the meridian to a longitude of the box,
  which is 180.0 and 17.0 there, and the 5-degree test that sends a layer
  to the comparison of candidates reads the same figure.  That test
  changes its answer only for a layer whose longitudes are given outside
  -180 to 180 and lie within 5 degrees of the zone's central meridian:
  points at longitudes 181 and 183, 2 degrees from the central meridian
  of UTM zone 1, were said to reach 360.0 degrees from it and got an
  Albers CRS, and now keep UTM zone 1, as a layer at -179 and -177 does;
  points at -181 and -183 now keep UTM zone 60 in the same way, as a
  layer at 179 and 177 does.  Other layers get the CRS and the
  `crs_choice` they got.

* **A row of a lon/lat layer with a coordinate that is not finite or lies
  beyond 1e150 shaped the projection chosen for the other rows.**
  `ensure_projected()` read the extent, the longitudes of the antimeridian
  test and the centre on every row.  Among 300 points near (-80.4, 45.6), a
  row at longitude 1e200 was taken for a layer across the antimeridian, and
  the call stopped with sf's "invalid crs" for a `lon_0` of -2.95e190, as
  did every function that projects lon/lat input through it
  (`resolution_profile()` among them); a row at latitude 1e200 under
  `purpose = "area"` gave a Lambert azimuthal projection in place of the
  Albers conic the other 299 points get, and moved the second point by
  21 km; the 39 stations above with a 40th row at longitude `Inf` were put
  in UTM zone 14.  Such rows also raised sf's "bounding box has
  potentially an invalid value range for longlat data".  The projection is
  now chosen with those rows empty, so the layer gets the CRS, and the other
  rows the coordinates, it gets with them empty, without that warning.  A
  layer in which no row has a location gets Web Mercator (Equal Earth for
  `purpose = "area"`), and the log line says so, where it
  quoted a centre at "lon = NaN, lat = NaN".  Moved to another CRS, the
  rows come back empty in their own geometry type: by `ensure_projected()`
  from lon/lat, with or without `target_crs`, and from a projected CRS with
  `target_crs`; by `harmonize_crs()`; and wherever a function moves a
  `boundary`, `cells_sf`, `polygons_sf` or `cover`, or the polygons of
  `cell_values_at()`'s `values`, into the CRS of another layer.  A
  polygon or line with one such vertex lost that vertex in
  `sf::st_transform()` and came back as a smaller finite shape, a
  half-degree square near (-93, 45) as a triangle of about half its area:
  `create_grid_polygons()` laid 55 cells over that triangle with nothing
  said, where it refuses the same row empty, and
  `assign_features_to_polygons()` joined a point inside the triangle to the
  square.  A vertex at a northing of 1e200 in EPSG:32616 came back at a
  finite northing in EPSG:32615.

* **One coordinate that is not finite or lies beyond 1e150 made a layer
  without a CRS read as planar.**  The lon/lat test of a layer without a
  CRS read its bounding box, which one such row made infinite or put off
  the lon/lat envelope.  Five points in degrees near (-93, 45) and a sixth
  at x = `Inf` were left in degrees by `ensure_projected()` with no
  warning, and `prep_model_data()` returned them in degrees for the fit;
  with a target CRS, or wherever a layer without a CRS takes the CRS of
  another (`harmonize_crs()`, a `boundary` beside projected points), the
  target was stamped on the degrees with a warning that they "do not look
  like lon/lat"; `build_tessellation()` on them stopped with GEOS's
  "encountered NaN/Inf numbers".  A row at latitude 1e200 did the same.  A
  boundary in metres without a CRS, given with points in EPSG:4326, was
  given EPSG:4326 when one vertex was infinite, and the call then stopped
  with "`boundary` must be polygonal".  The test now reads the box of the
  coordinates that have a location: the points are taken as EPSG:4326 and
  projected, with the warning, and the boundary is refused with the error
  that names both layers.  A boundary without a CRS that has no such
  coordinate (no rows, or every coordinate infinite), given with points in
  a geographic CRS, is now refused with an error that says so; it was given
  the points' CRS and refused later as "not polygonal".

* **A polygon or line with a vertex at `NA`, and a row edited in place to a
  coordinate beyond 1e150, were read as located.**  The test of a row's
  location read the layer's bounding box first and read no vertex when the
  box was finite. `sf::st_bbox()` skips `NA`, so a ring with one vertex at
  (`NA`, y), which a file can carry as `nan`, passed: `prep_model_data()`
  kept it, the fits and `predict_surface()` gave it a value, `make_folds()`
  and `fold_separation()` read it as a point, `harmonize_crs()` moved it to
  a finite shape, and the readers of `boundary` and `cover` used the rest
  of its ring, where the same row with a vertex at `Inf` is left out. sf
  also keeps the box it stored on a geometry column when a row is replaced
  with `[[<-`, so a point edited to (1e200, y) or (`Inf`, y) passed:
  `make_folds()` put it in a fold with nothing said, where the same layer
  built whole leaves it out. The test now reads the coordinates
  themselves. `kriging_adequacy()` now reads `cells_sf` by the same rule: a
  cell with a vertex at `NA`, `Inf` or beyond 1e150 is left out (`kr_`
  columns `NA`) with a warning that says why, as an empty cell is. It was
  repaired into the rest of its ring and kriged, and at 1e200 the call
  stopped with sp's "invalid Polygon object". `print()` of the result
  counts empty cells apart from thin ones.

* **`coerce_to_points()` gave a polygon with a vertex that is not finite or
  lies beyond 1e150 a point, and `assign_features_to_polygons()` joined
  such a polygon of `polygons_sf`.**  A line with such a vertex got an
  EMPTY POINT, but a polygon got a ring vertex under `"auto"` and
  `"point_on_surface"`, (`Inf`, y) or (5e199, y) under `"bbox_center"`, and
  (`Inf`, y) under `"centroid"` at 1e200. A geometry of any type with such
  a vertex now gets an EMPTY POINT in every mode. In
  `assign_features_to_polygons()` a polygon with a vertex at `NA` or `Inf`
  caught the points inside the rest of its ring, and with a vertex at `Inf`
  some of its neighbours' too (15 points in a square that holds 5), and one
  at 1e200
  stopped the call with GEOS's "orientationIndex encountered NaN/Inf
  numbers". Such a polygon now takes no point, and the points in it are
  unassigned, with the warning about unassigned features.

* **`coerce_to_points()` passed on sf's warning "st_point_on_surface may
  not give correct results for longitude/latitude data", and the
  summaries muffled it also where it was true.**  Every function that
  reduces a polygon layer in longitude and latitude to points
  (`prep_model_data()` and the fits, the `cv_*()` functions,
  `make_folds()`, `fold_separation()`, `estimate_sac_range()`,
  `resolution_profile()`, `determine_optimal_levels()` and others) raised
  it in sf's words, while `summarize_by_cell()` and `cell_values_at()`
  muffled it. The point lies on the polygon read on the plane in degrees,
  the reading `?coerce_to_points` states, and the warning is no longer
  raised; the points do not move. A polygon written across longitude 180
  (from 179.5 to -179.5) GEOS reads the long way round the globe, and its
  point lay about 12,000 km from it: the summaries read it there with
  nothing said. Such a polygon now gets a point on it, taken with its
  longitudes read from 0 to 360; one whose point lies on it neither way, or
  that s2 cannot read, keeps GEOS's point with a warning of class
  `spatialkit_surface_point_off`. A layer without a CRS whose coordinates
  look like longitude and latitude, which the functions that project it
  read as EPSG:4326, is read the same way, and `coerce_to_points()` says
  which rows that reading moves. The labels of `plot_tessellation_map()`
  on such cells are placed on them too.

* **`coerce_to_points(mode = "bbox_center")`, and `mode = "centroid"` with
  s2 off, put a lon/lat geometry written across longitude 180 near
  longitude 0.**  The box of a square from 179.5 to -179.5 spans the globe
  on the plane, and its centre, and GEOS's centroid in degrees, lay at
  longitude 0, about 17,600 km from it at latitude 10, with nothing said.
  `prep_model_data(pointize = "bbox_center")` then read a layer of such
  squares as one that circles the North Pole and projected it in a Lambert
  azimuthal centred on the pole, where `"auto"` centres the projection on
  longitude 180.  Such a geometry now gets the centre of its box, or GEOS's
  centroid, with the longitudes below 0 read from 180 to 360, brought back
  into -180 to 180, as the points on the surface are read; a geometry
  wider than 180 degrees in both readings keeps the planar point.  A layer
  without a CRS whose coordinates look like longitude and latitude is read
  the same way, with the warning that names the rows that reading moves.
  The same holds for the centroids `"line_midpoint"` takes of the rows
  that are not a LINESTRING, and `"auto"` of a MULTIPOINT and of the rows
  that are not points, lines or polygons, such as geometry collections: with
  s2 off, `"auto"` put a MULTIPOINT of (179.5, 10) and (-179.5, 10.5) at
  (0, 10.25).  The centroids s2 takes and every other row's point are
  unchanged.

* **`prep_model_data()` kept a row with a coordinate beyond 1e150 as a
  located row.**  It dropped a row whose geometry was empty or whose
  coordinates were not finite, but a point at (5e5, 1e200) in EPSG:32617
  was kept, and so was the same point in EPSG:32616 beside a boundary in
  EPSG:32617, which the move put at a finite northing of 5,766,270.  The
  fits read it as a location (`fit_gwr_model()` fitted 200 locations where
  the same layer with the row empty gives 199), and the grid of
  `predict_surface()` and the extent check of `cv_block_size_sweep()` read
  the extent of the prepared layer with the row in it.  The bound is now
  read on the coordinates as given, before the projection, and again after
  it, so such a row is dropped with the rows whose geometry is empty, as
  the same row emptied is, and every fitter and `cv_*()` function that
  prepares its data here reads the layer without it.  A polygon or line
  with such a vertex is dropped the same way.  The `predict()` methods of
  `rf_fit`, `gwr_fit` and `bayesian_fit` moved `newdata` into the fit's CRS
  before this reading, so a row of `newdata` at (5e5, 1e200) in EPSG:32616
  beside a fit in EPSG:32617 was predicted at the finite point the move gave
  it (0.391 from a forest on the coordinates, where the same row emptied
  gives `NA`).  They now empty such rows as given, before the move, and
  return `NA` for them; so does `predict_surface()` for such a row of a
  `grid` that is moved into the fit's CRS or reduced from polygons (the row
  comes back empty), and a row of `covariates` without a location is no
  longer taken as the nearest feature of a grid point.  The drop warning, and
  the log line of `require_response = FALSE`, count these rows as rows
  "without a location (an empty point, or a coordinate that is not finite
  or lies beyond 1e150)", where they said "an empty or non-finite
  geometry".

* **`estimate_sac_range()`, `make_folds()`, `resolution_profile()`,
  `determine_optimal_levels()` and `fold_separation()` kept a row at `Inf`
  or beyond 1e150, or read it after a move.**  They dropped a row whose
  geometry was empty or whose coordinates were missing, or not finite after
  the projection.  A point at x = `Inf` or at a northing of 1e200 in
  EPSG:32617 made `estimate_sac_range()` return `NA` with the reason "the
  points have no extent", where the same row empty gives a range; at
  latitude 1e200 in EPSG:4326 it stopped with FNN's "Data include NAs".
  `make_folds(method = "block_kfold")` stopped with base R's "'to' must be
  a finite number" on the `Inf` row and gave the 1e200 row a block of its
  own (2 folds, of 299 points and 1, where the same row empty gives 5),
  `random_kfold` put either row in a test fold, and `cv_block_size_sweep()`
  stopped in `make_folds()` on either.  `resolution_profile()` returned a
  WSS of `NA` at every level with the 1e200 row, and with `sac` in
  EPSG:32617 the same row in EPSG:32616 or EPSG:3857 was moved to a finite
  point and the profile built around it; `determine_optimal_levels()`
  stopped with base R's "missing value where TRUE/FALSE needed".  A
  prediction point of `fold_separation()` or of `make_folds(method =
  "nndm")` and `"knndm"` at 1e200 in EPSG:32616 was moved into the data's
  EPSG:32617 and measured there, and `fold_separation()` counted such a
  point in its print ("prediction target: 60 point(s)" where 59 have a
  location).  A latitude beyond 90, finite as given, came back `NaN` from
  the projection and stopped `estimate_sac_range()` and `make_folds(method
  = "nndm")` and `"knndm"` with "Data include NAs".  The bound is now read
  on the coordinates as given, before the reduction to points and before
  any move, and again after the move: such a row is dropped with the rows
  whose geometry is empty, as the same row emptied is, and so is a lon/lat
  polygon with such a vertex.  The drop warnings and log lines count these
  rows as points "without a location (an empty point, or a coordinate that
  is not finite or lies beyond 1e150)", where they said "with empty or
  non-finite coordinates".  A row that loses its location in the
  projection is dropped too: `estimate_sac_range()` and `make_folds()` log
  a line that says it has none once projected, `resolution_profile()` and
  `determine_optimal_levels()` count it in their drop warning, and
  `fold_separation()` leaves it out of the distances as it leaves out an
  empty row.  `resolution_profile()` drops a row that loses its location in
  the move into the CRS of `sac` with a warning that says so, where the
  move stopped the call with "some points have no coordinates in the CRS of
  `sac`"; with `select_on = "split"` the halves lose that row.
  The record of the rows' locations that `make_folds()` keeps, and that the
  `cv_*()` functions and `fold_separation()` check folds against, leaves a
  row at `Inf` or beyond 1e150 out as it leaves out an empty one: in
  EPSG:3857 a row at x = 1e200 made `make_folds()` wait without end on
  PROJ's move of it to EPSG:4326, and in a layer without a CRS a row at
  1e200 widened the check's tolerance to 1e191, so folds built on other
  data of the same size were accepted.  `cv_block_size_sweep()` now runs
  on such a layer.  An empty prediction point is no longer counted in
  `make_folds()`'s `params$n_prediction_points`.  Layers whose rows all
  have a location give the same results as before.

* **`make_folds()`, the checks of supplied folds and `compare_models()`
  stalled on a projected row far beyond the Earth.**  The record of the
  rows' locations that `make_folds()` keeps, and the record of a fit's rows
  that `compare_models()` checks its AICc and LOOIC with, move the rows to
  EPSG:4326, and PROJ moves a point out of Web Mercator in time that grows
  with x: 1.6 s at x = 1e17, more than 20 s at 1e19, and no return within a
  minute at 1e20.  Such a row is within 1e150 and has a location, so
  `make_folds()` on every method and the functions that build folds with it,
  `compare_models_cv()`, every `cv_*()` function, `fold_separation()`,
  `kriging_adequacy()` and `area_of_applicability()` given folds, and
  `compare_models()` with two or more AICc or LOOIC values did not return on
  an EPSG:3857 layer with a row at x = 1e20.  A row of a projected layer
  with a coordinate beyond 1e10 in the CRS's own unit, further out than any
  projected CRS places a point of the Earth, is now left out of both
  records, as an empty row is: the folds, the fits and the rows each
  function keeps are unchanged.  The same move stalled
  `ensure_projected(purpose = "area")` and `summarize_by_cell(area = TRUE)`
  on polygons with such a vertex, which now find the area distortion
  unmeasured, and `get_voronoi_seeds()` and `voronoi_seeds_kmeans()` given a
  resolution profile and points in degrees, with or without a CRS, which now
  return such a centre of the profile empty.  `harmonize_crs()`, and the
  shared step by which the package's functions move one layer into the CRS
  of another or of a `crs` argument, now give such a row back empty when the
  target is in degrees: `get_voronoi_seeds()` with `method = "kmeans"` and a
  boundary in degrees stalled on such a sample point, and drops it now with
  its warning about points without a location, as it drops an empty one.
  Layers in degrees, without a CRS, or whose rows all lie within 1e10 of
  the origin, and moves into a projected CRS, give the same results as
  before.

* **The seeding functions, the Voronoi cells, `clip_target_for()`,
  `assign_features_to_polygons()`, the variogram design effect of
  `summarize_by_cell()` and `kriging_adequacy()` read a point at `Inf` or
  beyond 1e150 as a location, or read it after a move.**  With one of 300
  sample points at a northing of 1e200, `voronoi_seeds_kmeans()` and
  `get_voronoi_seeds(method = "kmeans")` stopped with "no start gave a
  partition of the 300 point(s) into 15 cells"; the same point in
  EPSG:32616 beside a boundary in EPSG:32617 was moved to a finite point and
  clustered with nothing said.  `create_voronoi_polygons()` and
  `build_tessellation(method = "voronoi")` with one of 30 seeds at x =
  `Inf` or a northing of 1e200 stopped with GEOS's "Could not locate
  vertex", "encountered NaN/Inf numbers" or "distance must be a finite
  value", and with the seed in another CRS than the boundary's with "Points
  of LinearRing do not form a closed linestring"; `clip_target_for()`
  stopped the same way at `Inf` and returned a target of infinite area at
  1e200.  One point of a `cover` at `Inf` or 1e200 set the whole cover
  aside, with the warning that it gives no polygon, and the cells covered
  5,334 km^2 where the same point empty gives 10,516 km^2; a polygon
  `cover` in EPSG:32616 with a vertex at a northing of 1e200 was moved
  into EPSG:32617 as a boundary of 68,562 km^2, where the same row empty
  gives 8,983 km^2.  `assign_features_to_polygons()` joined a point at 1e200 in
  EPSG:32616 to the polygon in EPSG:32617 where the move put it, and its
  warning said that a row at `Inf` or 1e200 falls inside no polygon and
  advised checking that the polygons cover the features.
  `summarize_by_cell(deff = "variogram")` correlated a point at 1e200, as
  given or moved into the variogram's CRS, with the others of its cell (a
  standard error of 0.543 where the same row empty gives 0.664), and
  `kriging_adequacy()` stopped in `gstat::krige()` on a point at `Inf` and
  kriged from one at 1e200, as given or moved to a finite point (moved, its
  folds held 399 points and 1).  The
  bound is now read on the coordinates as given, before any move, and
  again after it, and such a point is read as an empty point: dropped by
  the seeders with their warning, `NA` in `index` and counted with the
  empty seeds and points, left out of the clip target's box and of the
  cover, matched to no
  polygon, and left out of the correlation and of the kriging.  The
  warnings and log lines count these points as points "without a location"
  or with "no location (an empty point, or a coordinate that is not finite
  or lies beyond 1e150)", where they said "with empty or non-finite
  coordinates" or "empty geometry", and the refusal of `kriging_adequacy()`
  counts the points with a finite response and a location.  Layers whose
  rows all have a location give the same results as before.

* **`area_of_applicability()` gave a location index to a row beyond
  1e150, and refused training data with one for a false reason.**  With
  location measured (`coords = "yes"`, or a fit whose predictions depend
  on it), a row of `newdata` at (5e5, 1e200) in the training CRS got a
  `DI_location` of 2.6e149 and `AOA = FALSE`, and the same row in
  EPSG:32616 beside training data in EPSG:32617, moved to a finite point,
  got 15.1; the same row emptied gets `NA` in both.  Such a row in
  `train_sf` stopped the call with "every training row is at the same
  location", which was false.  The bound is now read on the coordinates as
  given and again after the move, so such a row of `newdata` gets `NA` in
  `DI_location` and `AOA`, and such a training row is left out of the
  reference distances, as the same rows emptied are.  For a forest fitted
  with `include_coords = TRUE` the row's `"..x"` and `"..y"` are `NA`
  (they held 1e200, or `Inf` for an infinite coordinate), and a row of
  `newdata` reprojected into the training CRS comes back with an empty
  geometry.  The refusals for too few usable training rows, and for
  `folds` beside training rows that are left out, now name a missing
  location as a cause where location is measured; the second said
  "missing predictor values" for a row with an infinite predictor or no
  location.

* **`area_of_applicability()` dropped a predictor at a large offset as
  having no variance.**  Its test kept a predictor only where the standard
  deviation was more than 1.5e-8 of the largest absolute value, which reads
  the column's distance from 0, not its spread.  A predictor with standard
  deviation 1 at 1e8, or a timestamp in seconds, was dropped with the log
  line "dropping predictor(s) with no variance in the training data", the
  index was computed on the other predictors, and new rows were read on
  that predictor through the slack of a dropped column instead: at 1e8,
  rows more than 1.5 beyond the training range got `DI = Inf`, the others
  were judged on the other predictors alone, and the threshold fell from
  0.306 to 0.048; at 1.7e9, 19 of 20 rows 4 standard deviations out were
  inside the AOA, and at 1e15, 20 of 20 rows 6 out.  A predictor the test drops is now
  kept where its values differ by more than their rounding (3 units in the
  last place of the largest absolute value), with its training mean and
  standard deviation in `scaling`, so `dropped_vars`, `DI`, the threshold
  and `AOA` are those of the same predictor at 0.  The same reading keeps
  a predictor in units whose squares vanish or overflow (standard deviation
  1e-170 or 1e200), which was dropped too.  A predictor still
  dropped reads a new row as a different value when it lies outside the
  training range by more than 3 units in the last place of the training
  values, where the slack was 1.5e-8 of their size: 25 for a constant
  timestamp at 1.7e9, so a row 25 seconds away was not extrapolation.  For
  a predictor that is 0 throughout, the slack is 3 units in the last place
  of 1 (6.7e-16) where it was 1.5e-8.  Predictors the first test keeps
  give what they gave.

* **`gp_lengthscale_bounds()` read every row of `coords_xy`.**  A row of
  `NA` stopped it with base R's "missing values and NaN's not allowed if
  'na.rm' is FALSE", and a row at `Inf` or 1e200 made the upper bound
  `Inf` with no warning.  A row with a coordinate that is missing, not
  finite or beyond 1e150 is now left out, with a warning of class
  `spatialkit_rows_dropped` that gives the count, and the bounds are those
  of the other rows.

* **The projection chosen for a small lon/lat MULTIPOINT layer depended on
  the order of the points in each feature.**  A layer of features that
  are not points, with fewer than 40 distinct locations, is scored on its
  features' vertices too, with points added along the edges between them.  A MULTIPOINT has no
  edges, but points were added on the segments between consecutive members
  of each feature, at places the layer does not have, and they moved with
  the order of the members: two MULTIPOINTs of five points each scored
  0.75 percent in UTM zone 14 in one order and 0.78 percent in another.  A
  MULTIPOINT layer is now scored on its member points with nothing added
  between them, which is the score of the same points as a POINT layer.
  For such a layer the CRS, `attr(, "crs_choice")` and the logged figures
  can differ from before; MULTIPOINT layers with 40 distinct locations or
  more score as before.

* **A lon/lat layer that mixes geometry types or holds
  GEOMETRYCOLLECTIONs was scored without its outline, and across 180
  degrees `ensure_projected()` stopped on it.**  `sf::st_coordinates()`
  refuses such a layer (class `sfc_GEOMETRY` or `sfc_GEOMETRYCOLLECTION`,
  as `rbind()`, `c()`, `sf::st_make_valid()` and `sf::st_intersection()`
  can give), so the vertices that score a layer of fewer than 40 features
  were not read.  Two boxes, one listed as a MULTIPOLYGON, were scored on
  the one pair between their interior points and got a Lambert azimuthal
  CRS with "0.12% against the zone's 0.47%", where the same boxes as
  MULTIPOLYGONs get an Albers CRS with "1.13% against the zone's 6.57%";
  a box over the conterminous United States listed as one
  GEOMETRYCOLLECTION scored `NA` for every candidate and kept UTM zone 15.
  The test for a layer that straddles longitude 180 read the coordinates
  the same way, outside any handler, so `ensure_projected()`, and the
  builders and fits that project a lon/lat layer through it, stopped with
  sf's "not implemented for objects of class sfc_GEOMETRY" on such a layer
  across 180.  Such a layer is now read feature by feature, each collection
  unpacked into its members, so it scores, chooses and logs as the same
  features in one type.  The points of a POINT or MULTIPOINT among them
  are not filled between either.  `sf::st_is_empty()` calls a
  GEOMETRYCOLLECTION, MULTIPOLYGON or MULTILINESTRING empty when its first
  part is empty, and such a feature was left out of the score: a box
  listed after an empty line in one GEOMETRYCOLLECTION scored `NA` for
  every candidate.  Empty parts are now dropped before empty features, so
  such a feature is scored on its other parts.  Layers of one geometry
  type with no empty part score as before.

* **`ensure_projected()` read a small CRS-less layer of mixed geometry
  types as planar, and gave a layer of mixed types that s2 refuses
  EPSG:3857.**  With no CRS, a layer under one degree across is taken as
  lon/lat when its coordinates have fractional degrees, and that test
  read the coordinates with `sf::st_coordinates()`, which refuses a layer
  of mixed types or of collections: two boxes of 0.3 degrees in Germany,
  one listed as a MULTIPOLYGON, were returned with no CRS and
  `crs_assumed` `"none"`, where the same boxes as MULTIPOLYGONs are taken
  as EPSG:4326 and projected to UTM zone 32.  A layer whose only feature
  has an empty first part, which `sf::st_is_empty()` calls empty, was read
  as planar too.  When s2
  refuses a lon/lat layer even after the repair (a bowtie polygon), its
  centre is the mean of its vertices, which were read the same way: a
  bowtie beside a MULTIPOLYGON got no centre and EPSG:3857, where the same
  features as MULTIPOLYGONs get a CRS centred on them.  Both now read such
  a layer feature by feature, as the same features in one type.

* **A lon/lat layer with Z or M coordinates was not scored as its XY
  copy.**  The outline of a layer of fewer than 40 features is read ring by
  ring, and the key of a ring read Z and M as well, so an edge between two
  vertices with other Z or M was not filled: a box over the conterminous
  United States whose corners differ in Z got a Lambert azimuthal CRS with
  "0.97% against the zone's 5.42%", where its XY copy gets an Albers CRS
  with "1.13% against the zone's 6.71%".  A line with Z kept its Z on the
  point that represents it beside XY vertices, and GEOS refuses M, so on a
  LINESTRING Z, or on a layer of lines, polygons or MULTIPOINTs with M,
  every candidate scored `NA` and the UTM zone was kept.  The layer is now
  scored in x and y, and the ring key reads the ring and part columns
  alone, so such a layer scores, chooses and logs as its `sf::st_zm()`
  copy.  The layer returned keeps its Z and
  M.  Layers in XY score as before.

* **`summarize_by_cell(area = TRUE)` refused projected cells when s2
  refused one of them in degrees.**  The area check compares the cells'
  planar areas with their geodesic areas, which s2 measures on the cells'
  copy in EPSG:4326, and it took the geodesic areas in one call.
  Two corners a fraction of a nanometre apart in the cells' CRS can fall on
  one point in degrees, where s2 refuses the ring ("duplicate vertex"), and
  then the call failed for every cell.  A vertex listed twice, which GEOS
  accepts, does the same.  The k-means cells of a 30 by 20 lattice at 30 m,
  built and kept in EPSG:32615 from rows in degrees, are all valid in
  EPSG:32615 and s2 refuses 4 of the 133; the request was refused because
  "the area distortion of `cells_sf`'s CRS could not be measured".
  `ensure_projected(purpose = "area")` on projected input skipped its check
  in silence the same way: 20 cells of 1 degree from 30N to 69N in Web
  Mercator, which distort areas by 214 percent, logged nothing when one
  cell had such a pair of corners.  When the one call fails, the geodesic
  areas are now taken one cell at a time, and a cell s2 refuses is measured
  again without the vertices that repeat the one before them in degrees,
  which leaves its points and its area as they were.  The lattice's cells
  get a density, with `attr(, "area_error")` the spread over all 133
  cells, and the warning is logged for the Web Mercator cells.  A cell s2
  refuses even then, such as one whose edges cross, still leaves the
  spread over the cells unmeasured: `summarize_by_cell(area = TRUE)`
  refuses the request (since GEOS calls such a cell invalid, with the
  error on invalid cells below), and `ensure_projected()` measures the
  grid over the extent instead.  A layer whose cells s2 all accepts gets
  the same figure as before.

* **`ensure_projected(purpose = "area")` skipped its check in silence on a
  projected layer of one polygon, and `summarize_by_cell(area = TRUE)`
  refused one cell for a cause it did not have.**  One polygon gives one
  ratio of planar to geodesic area, which has no spread, so the distortion
  came out unmeasured and nothing was logged: a box from 100W to 90W and
  from 30N to 60N in Web Mercator logged nothing, where the same box split
  in two logs 26.2 percent.  A layer with a polygon s2 refuses in degrees,
  such as one whose edges cross, was passed over the same way: 20 cells of
  1 degree from 30N to 69N in Web Mercator, one of them a bowtie, logged
  nothing.  `summarize_by_cell(area = TRUE)` refused a single cell, a 10 km
  square in UTM zone 15 among them, because "the area distortion of
  `cells_sf`'s CRS could not be measured (no polygon areas, or geodesic
  areas unavailable)", though the cell has both; a single lon/lat cell was
  refused the same way.  Both functions now measure such a layer on a 6 by
  6 grid over its bounding box, as they measure a layer without polygons.
  `ensure_projected()` logs the grid's figure when it exceeds 1 percent
  (65.2 percent for the box, 124.4 percent for the cells with the bowtie),
  and logs that the distortion could not be measured when the grid gives
  no figure either (a single point, or a line along one axis).
  `summarize_by_cell()` gives a single cell its density, with
  `attr(, "area_error")` the grid's figure, or refuses it for the
  distortion.  A layer whose polygons give a figure gets the same figure
  and the same log as before.

* **`summarize_by_cell(area = TRUE)` gave a projected cell whose edges
  cross a density over GEOS's area of its ring.**  The area check measures
  up to 200 of the cells, together with the cells that reach furthest; a
  cell s2 refused among them made the request refused, but one outside them
  got `sf::st_area()`'s area, which for a ring whose edges cross is the
  signed sum of its lobes.  Of 400 squares of 10 km in EPSG:32615 with a
  bowtie of 4.25e7 square metres, the bowtie got a `cell_area` of 2e7
  square metres at cell 148 and the request was refused with the bowtie at
  cell 1, so the outcome turned on which cells the check took.  Every
  non-empty projected cell is now checked with `sf::st_is_valid()` in its
  CRS before any summary, and a cell GEOS calls invalid makes the request
  refused, with the first such cells, GEOS's reason and the remedy: repair
  the cells with `sf::st_make_valid()` and pass them again.  A cell GEOS
  calls invalid that s2 measures in degrees, such as one with a hole
  outside its shell, is refused by the same check.  Valid cells get the
  densities and `attr(, "area_error")` they got before.
  `ensure_projected(purpose = "area")` on projected input read the same
  planar area as a figure about the CRS: three 10 km squares of UTM zone
  15 beside a square with a hole outside its shell logged a distortion of
  40.0 percent, the difference between the shell minus the hole and the
  shell plus the hole.  A layer with a polygon GEOS calls invalid is now
  measured on the 6 by 6 grid over its extent, which logs nothing for that
  layer.

* **`resolution_profile()` reported a layer too small for a ladder as
  supported.**  With fewer than two cells at `min_cell_n` it raised the
  ceiling to 2 with `supported = TRUE`.  It still scores two cells, with
  `supported = FALSE`, `ceiling_from = "fewer than two supported cells"` and
  a warning of class `spatialkit_not_supported` that says what would work.

* **A constant response got a finite `cp`.**  A response of 1e6 plus
  rounding noise was scored on its rounding error and `select_resolution()`
  returned 2 cells for it.  `cp`, `moran_i` and `moran_z` are now `NA` at
  every level, with a warning of class `spatialkit_constant_response`.

* **`determine_optimal_levels()` returned 1 for a layer with no rows.**  It
  is now an error.

* **The leakage warning of `make_folds(method = "block_kfold")` recommended
  settings that stopped with an error or returned fewer folds.**  It ended
  "Pass block_size = <range> or auto_range = TRUE" whatever the layout.  The
  warning now builds the grid the advice asks for and says what comes back:
  the setting when that gives `k` folds, and otherwise the largest block
  size that does and `method = "nndm"` (`"knndm"` above 5,000 rows) as the
  design that needs no block size.  The single-block error names that size too.

* **`make_folds()` returned fewer folds than asked with only a log line.**
  `"random_kfold"` with more folds than rows, `"block_kfold"` with fewer
  than `k` occupied cells and `"leave_location_out"` with more folds than
  groups now raise a warning of class `spatialkit_k_reduced`, which carries
  `k_requested` and `k` and, where there is one, gives a setting that works.
  `params$k_requested` records the `k` asked for, and `params$block_width`
  and `params$block_height` the cells built.  `cv_block_size_sweep()` gathers
  the reductions of its sizes into one warning, and its table gains the
  columns `block_width` and `block_height`.

* **The default ladder of `cv_block_size_sweep()` gave other blocks in
  another unit.**  Its bottom rung is a twenty-fifth of the shorter side,
  which that side holds exactly 25 times.  The rung was computed as
  `exp(log(side / 25))`, a hair above or below it, and the blocks per side
  as `floor(side / size)`, so the last bits of the unit decided between 24
  and 25 blocks.  On an 8 x 8 lattice at 250 m the bottom rung's blocks
  were 72.9 m wide in metres and 70 m in units of 1e-6 m, and
  `blocks_used` changed with the unit on other layers.  The ladder now
  starts at `side / 25` and ends at its top exactly.
  `make_folds(block_size = )`, which builds each rung's grid, now reads a
  block size within one part in ten million above an extent over a whole
  number as that number of blocks, called from the sweep or directly.  The
  bottom rung gives 25 blocks along the shorter side in every unit.

* **`cv_block_size_sweep()` stated its fits wrongly in four texts.**  With
  `block_sizes` supplied, the refusal over `max_fits` advised lowering
  `n_sizes`, which only the default ladder reads; it now advises fewer
  `block_sizes` or a smaller `k`.  On either ladder, when the count
  includes the random reference, it also names `include_random = FALSE`,
  which at `n_sizes = 1` and `k = 2` is the only way to cut the count.
  When every size of the default ladder needs a grid of more than
  1,000,000 blocks (points on a line a hair wide), the sweep refused with
  "no block size leaves at least k = 5 blocks ...; pass smaller
  `block_sizes`", where the sizes held too many blocks and smaller ones
  are refused; it now says so and names the largest size whose grid holds
  `k` blocks, to pass as `block_sizes`.  The header of `print()` and the
  caption of `plot()` gave the fits planned before any fold was built, and
  the caption said `k` folds per size, where a size whose grid has fewer
  than `k` cells holding points builds fewer.  Both now count the fits made
  from the `k` column, and the caption says "up to `k` folds per size" when
  a size built fewer.  `attr(, "n_fits")` stays the count planned, which
  the `max_fits` check reads.

* **`fold_separation()` printed "this score is optimistic" from a share that
  cannot give a direction.**  Which way a score errs depends on how the
  held-out distances compare with the distances from the prediction
  locations to the data.  Without `prediction_points` the function now
  reports the distances and draws no conclusion.  The column `within_range`
  (the share of held-out points closer than the range) is replaced by
  `mean_cor`, the mean correlation between a held-out point and its nearest
  training point; code that reads the old column has to be changed.
  `n_unknown_ids` counts each row once; it counted fold entries, so five
  missing rows read as 25 at `k = 5`.

* **`make_folds(method = "nndm")` returned plain leave-one-out without a
  warning when the prediction points were nearer to the data than
  leave-one-out can reach.**  Excluding neighbours only lengthens a held-out
  distance.  A warning of class `spatialkit_nndm_out_of_reach` now says that
  the cross-validated error is an upper bound for those locations, and
  `params` gains `max_ecdf_deficit`, `loo_ecdf_deficit`,
  `target_out_of_reach` and the Wasserstein distances.

* **`cv_spatial()`, `cv_rf()`, `cv_gwr()`, `cv_bayes()` and
  `compare_models_cv()` dropped `block_size` and `auto_range = TRUE` without
  a word when `folds` was supplied.**  They raise a warning of class
  `spatialkit_fold_args_ignored`.

* **`cv_spatial()` scored a fold's predictions against `response_var`
  without checking that the fold's fit models it.**  A `fit_fn` that fits
  `price`, called with `response_var = "slope"`, had its predictions of
  `price` scored against `slope` on every fold, and the result was reported
  as the cross-validated score of a model of `slope`, with every fold
  `ok`.  A fold whose fit's `response_var` is not the call's now fails,
  with both names in `fold_status`, as a fold whose `fit_fn` returns no
  `spatial_fit` does; when every fold fails so, the all-folds-failed
  warning gives that cause.  The same holds in `cv_block_size_sweep()` and
  `select_features_forward()`, which run their `fit_fn` through the same
  folds; the hold-out score of `select_features_forward(select_on =
  "split")` is `NA` for such a fit, with a warning that gives both names.
  A named `response_var`, such as `c(r = "price")`, is compared by the
  column it names.
  The package's own fitters fit the call's response and are unchanged.

* **`compare_models(fits, newdata = )` labelled training-residual Moran's I
  as a `newdata` metric.**  A new column `resid_morans_basis` says
  `"in-sample"`, or `"out-of-bag"` for a forest.

* **`compare_models_cv()`, `cv_block_size_sweep()` and
  `select_features_forward()` repeated the same data warning for every
  model, size or candidate set.**  The warnings `prep_model_data()` raises
  about a layer are raised once per call.  With `parallel` above one core, a
  fold's classed warning keeps its class and fields.

* **`model_metrics(fit, newdata, draws = TRUE)` stopped with "missing value
  where TRUE/FALSE needed".**  `draws = TRUE` and `se.fit = TRUE` are
  refused by name: the metrics take one predicted value per row.

* **A predictor far from zero made `fit_gwr_model()` report every window as
  collinear.**  Windows whose condition index with the intercept (`cn`) was
  above 1e6 were counted, a rule for a solve of the uncentred normal
  equations.  A window's slopes now count as collinear on the slope index
  alone (`cn_slopes` above 30, or singular).

* **The GWR bandwidth search could select a bandwidth at which a term is
  estimable at no location.**  Such a bandwidth is not selected; when every
  bandwidth is like that the widest is used, with a warning and
  `info$bandwidth_is_fallback = TRUE`.

* **`area_of_applicability()` returned a threshold of 0 for repeat visits.**
  With `folds = NULL` a training row that repeats another's predictor values
  was its own nearest neighbour at distance 0: 60 sites visited 4 times put
  30 of 200 new points inside.  Each row's reference is now its nearest row
  with different values, and the same layout puts 188 inside
  (`n_duplicates` counts the repeats).  With `folds`, a threshold of 0 is a
  warning of class `spatialkit_aoa_zero_threshold`.

* **Without FNN, `area_of_applicability()` gave exact copies of training
  rows a small non-zero distance.**  About one copy in ten came out near
  1e-8 and could fall outside a threshold of 0.  The distance returned is
  now recomputed exactly.

* **`area_of_applicability()` counted points beyond the training extent as
  inside for a forest fitted with coordinates.**  Location was measured only
  through the forest's `..x` and `..y` columns, which under
  `weights = "auto"` took the forest's own importance, so a forest that
  barely split on its coordinates gave location almost no weight: in a test
  where half of 50 prediction points lay one half to one and a half extents
  beyond the data, 41 read as inside, 16 of them beyond the extent.  Such a
  forest now gets the location index a GWR and a Bayesian fit get, and none
  of the 25 is inside.  The predictor index is unchanged.

* **`select_features_forward()` always selected at least one variable with a
  forest or a GWR.**  A learner that could not fit an empty predictor set
  had no baseline, so the first variable was accepted whatever it scored.
  The baseline is now the learner's own null model (a GWR fits one),
  `null_fn`, or the mean of each training fold.  Under a fold-mean baseline
  the first variable is tested at the end and can be dropped
  (`baseline$first`); `plot()` draws the path as walked.  A forest sweep can
  still add noise variables after the real ones, because an extra variable
  regularises the forest; `?select_features_forward` says how to read that.
  `tol` takes any single number, `Inf` included.

* **`plot()` on an area of applicability with a location index counted the
  wrong rows beside the predictor threshold.**  It printed the rows outside
  on either index.  The subtitle now counts the predictor index, location
  and the total separately.

* **One row of `newdata` with a factor level the model was not fitted with
  made `predict()` on a `bayesian_fit` return `NA` for every row.**  Rows
  with such a value are `NA`, the rest are predicted, and a warning of class
  `spatialkit_new_levels` names the predictor and the levels.  Any other
  failed posterior draw is a warning of class `spatialkit_predict_failed`,
  where it was a log line.

* **`prep_model_data()` dropped incomplete rows with a log line only.**
  When preparing data for a fit it raises a warning of class
  `spatialkit_rows_dropped` with the number and share of rows, the count by
  reason and by column, and how far the response of the dropped rows sits
  from that of the kept rows.  A predictor with one value in every kept row
  raises `spatialkit_constant_predictor`.  Preparing `newdata` for
  `predict()` raises neither.  `attr(, "dropped")` gains `by_column`,
  `response_std_diff`, `constant_predictors` and `n_duplicate_locations`.

* **`prep_model_data()` measured the shift of the dropped rows on the
  rounding of a response at a large offset.**  `response_std_diff` and the
  `spatialkit_rows_dropped` warning's "a difference of X SD" were a
  difference of the raw means of the dropped and the kept rows, and a
  double holds each mean of a response at 1e15 only to an eighth.  On 40
  rows with 6 dropped, the same stored values gave 0.87 SD at offset 0 and
  0.93 SD at 1e15, and only the second added "more than two standard
  errors, so the rows are not missing at random".  Where the raw values of
  the response lose their spread, the means and the standard deviation are
  now taken less the response's centre, as `estimate_sac_range()` takes
  it, and the centre is added back to the means the message quotes.  The
  shift and the verdict are those of the same values at 0.  A response of
  one value up to rounding, whose standard deviation is its last bit, now
  gets `NA` and no shift sentence, as a response of one value does: with
  the dropped rows one unit in the last place above most kept rows, it
  read "a difference of 1.61 SD, more than two standard errors, so the
  rows are not missing at random" at 1.23, 4.5e6 and 1e15 alike.  Every
  other response is read as it was, to the last bit.

* `fit_bayesian_spatial_model()` stops before compiling when the response
  takes a single value, and raises a warning of class
  `spatialkit_gp_basis_size` when there are fewer rows than a tenth of the
  basis functions.

* **`fit_bayesian_spatial_model()` handed brms the response as given, so
  its fit depended on the response's origin and units.**  Stan starts every
  parameter between -2 and 2, and brms floors its default priors for the
  intercept, the residual SD and the GP's SD at a scale of 2.5 in the
  response's units.  At an offset of 1e8 the posterior mean fitted values
  missed the response by the offset; in units of 1e-20 they were 1e15
  response SDs off, after divergent transitions; near 1e160 the fit had no
  draws and said nothing; and `cv_bayes()` folds failed or scored the same
  way.  Under the gaussian, student, skew_normal, exgaussian and
  asym_laplace families with an identity link, a response whose centre is
  more than 1,000 times its spread from zero is now handed to brms less
  that centre, and one whose SD lies outside 2^-20 to 2^20 is divided by
  the power of two nearest it.  Under exgaussian and asym_laplace the rule
  is narrower, a centre more than 4 times the spread and an SD outside 2^-2
  to 2^2, because brms's default prior on exgaussian's `beta` and Stan's
  starting values do not follow the response: an exgaussian fit at an SD
  near 1,000, or 300 SDs from zero, had no draws, and an asym_laplace fit
  300 SDs from zero or at an SD near 1e-4 came back with `sigma` at twice
  its value or more.
  `fitted()`, `predict()` and its draws, `residuals()`, `coef()`,
  `$info$loo`, `$info$looic`, `$info$gp_adequacy$sdgp` and `$info$b_prior`
  are given in the response's units, and `$info$response_scaling` records
  the centre and scale.  Every other response is handed over as it is, and
  its draws do not change.  A `prior` on a class on the response's scale
  keeps the response as given, with a warning of class
  `spatialkit_response_scale`.  A response whose values differ only by
  rounding (1e15 and the next double) is refused as a constant one is; it
  was sampled.  `cv_bayes()` takes `yhat_sd` at any size, where `sd()` of
  draws near 1e160 overflowed.  The spatial share in `$info$gp_adequacy`
  is taken at any size of the response: its squares overflowed near 1e160
  and underflowed below about 1e-162, so the share was `NA` and the
  length-scale rule decided whether the GP is identified in its place.

* **`fit_bayesian_spatial_model()` returned a fit with no posterior draws
  and said nothing.**  When Stan rejected every starting value, rstan
  returned an empty fit, `$info$convergence_ok` was `NA`, `$info$looic`
  was `NA`, no warning was raised, and `fitted()` stopped later.  Such a
  fit is now an error that gives the response's median and SD as brms was
  given it and the family.

* **`summarize_by_cell()` let one infinite value through.**  An `Inf` in a
  response or predictor column made its cell's mean `Inf` and, under
  `deff = "kish"`, set the ICC to 0 for every cell.  Infinite values are now
  left out of every statistic of their column, as `NA` is, with a warning of
  class `spatialkit_nonfinite`.

* **`kriging_adequacy()` kriged a factor on its level codes.**  A response
  that is not numeric or logical is an error that names the column and its
  class.  A logical response is kriged as 0/1.

* **`kriging_adequacy()` capped the replicates' variance at the nugget.**
  Under a model with too small a nugget the means of repeat observations
  were kriged as more exact than their replicates show.  Their pooled
  within-location variance `sw2` is now used as it is: a mean of `m`
  replicates has error variance `max(c0 - sw2, 0) + sw2 / m` for a nugget
  `c0`.  A warning of class `spatialkit_replicate_var_above_nugget` gives
  both numbers, and `attr(, "replicates")` records them.

* **A subset of an assigned layer kept its parent's `"unassigned"` record.**
  `[` and the dplyr row verbs drop it, as they drop `"dropped"` and
  `"ties"`.

* `plot.resolution_profile()` says "fewer than two cells meet `min_cell_n`"
  of a profile that is unsupported for that reason; it said the floor
  exceeds the ceiling.  `plot_folds()` and `print()` on a
  `kriging_adequacy()` result label kNNDM folds.

* **`ensure_stable_poly_id()` stopped on a cell that is valid in its own CRS
  and crosses itself in lon/lat, and took `build_tessellation(method =
  "voronoi")` with it.**  The sort key is measured on a lon/lat copy with s2.
  Only the vertices are transformed, and s2 joins them with great-circle
  arcs, so a long straight edge that passes about a metre from another vertex
  of the same ring can come out on the other side of that vertex.  s2 refuses
  such a ring, the repair of the lon/lat copy does not split crossing edges,
  and `st_centroid()` and `st_area()` stopped with "Loop 1 is not valid: Edge
  36 crosses edge 52".  Met on 1 of 291 Voronoi cells of Texas clipped to the
  state outline in EPSG:5070, where a 36.9 km cell edge passes 1.47 m from a
  vertex of the coastline and the great circle through its ends passes 0.37 m
  beyond that vertex: no cell of the state got an ID.  A feature s2 still
  refuses after the repair now gets a second repair of its own sort copy with
  the crossing edges split, and if s2 refuses that too its centroid and area
  are measured as plane geometry in the layer's own CRS and the centroid is
  transformed to the sort CRS.  The function warns once, saying how many
  features took each route, and `create_voronoi_polygons()` and
  `build_tessellation(method = "voronoi")` pass the warning on.  Those
  features' keys are close to the spherical ones and not equal to them: the
  second repair moved the Texas cell's area by about 40 square metres in
  2,949 square kilometres, and its centroid taken in the plane of EPSG:5070
  lay 20 m from the one the second repair gave on the sphere.  Every other
  feature keeps exactly the key it had, so the other 290 Texas cells are
  numbered in the order they always would have been, the geometry returned
  is still the caller's own, and the same 291 cells get the same IDs in
  EPSG:5070, EPSG:3083, UTM zone 14N, Web Mercator and lon/lat.  A feature
  keyed in the plane is measured in the CRS the layer arrived in (the Texas
  cell's plane centroid lies 9 to 73 m from the spherical one depending on
  that CRS), so it can take another place in the order when the same layer
  arrives in another projection, and the IDs between the two places move
  with it; the warning says so whenever that route is taken.  With
  `make_valid = FALSE` no repair is tried, as the caller asked, and a
  feature s2 refuses goes straight to the plane: such a call used to stop
  with s2's error and now returns IDs with the warning.  All three `method`s
  are covered: for `"surface_point"` and `"bbox_center"` the point never
  went to s2, but the area, the third key of each, did, and both stopped
  there.  An error that is not s2 refusing a feature is raised as it was.

* **`resolution_profile()` took the noise variance of C_p from the nugget of
  a variogram fit that had not converged, whenever the range was refused for
  another reason.**  Its documentation says a model that did not converge
  gives `cp` no nugget, the nugget and sill being where the optimiser
  stopped, and the profile decided that by looking for "converge" in
  `rejected_reason`.  `estimate_sac_range()` gives a refused range one
  reason, and a falling variogram and a range past the largest lag come
  before the optimiser in it, so a fit that stopped at gstat's iteration
  limit with a range past the fitted lags is recorded as `"fitted range
  exceeds the largest lag fitted"`.  That is most fits that do not converge:
  over 200 simulated layers (trends, ranges longer than the extent, 40--60
  points, white noise) 96 of the 102 all-pairs fits that did not converge
  carried that reason and 6 `"variogram model did not converge"`.  With an
  east-west trend on an exponential field (300 points, nugget 0.1, 25 draws)
  every all-pairs fit stopped at the limit, every one was refused as past
  the lags, and C_p was given nuggets of 1.18--1.62; in the first draw the
  range reported was 29,180 against a largest lag of 647.  Whether the model
  converged is now read from the model
  (`attr(attr(sac, "variogram_model"), "converged")`, which
  `estimate_sac_range()` already set); a sac that carries no such flag (made
  by hand, or built from REML parameters) is judged by its reason, as
  before, and a reason that says the fit did not converge is believed
  whatever the flag says.  For a refused range whose model did not converge,
  `cp` takes its noise variance from the finest level whose cells hold at
  least two scored rows on average, `attr(x, "cp_noise")$source` reads
  `"finest-level residual mean square"`, `attr(x, "variogram")` is `NULL`,
  and the warning names both facts: why the range was refused, and that the
  model did not converge.  The existing notice that names the level and the
  variance `cp` was given follows it, so such a profile raises two warnings
  where it raised one.  The finest-level value can be larger or smaller than
  the nugget it replaces (larger in 41 of 52 simulated layers, smaller in
  11), so the `cp` pick can move either way: it stayed in 49, moved to fewer
  cells in 2 and to more in 1.  When no level holds two scored rows a cell,
  `cp` is `NA` at every level for such a sac, where it was finite; the
  warning's last clause says so, and the notice is not raised.  A fit that
  converged to a range past the lags, or on a falling variogram, still gives
  `cp` its nugget.  (`estimate_sac_range()` no longer fits with gstat's
  iterative fit: the flag is now `FALSE` when the criterion has no minimum
  inside the ranges searched, and the reason `"variogram model did not
  converge"` is carried only by results stored by earlier versions.)

* **A clipped cell that also touched the boundary from outside was dropped,
  and its points left with no cell.**  `st_intersection()` returns such a
  cell as a GEOMETRYCOLLECTION of its area and a line or point (on an
  L-shaped boundary, the cell over the inner corner), and the Voronoi, grid
  and Delaunay builders kept only POLYGON and MULTIPOLYGON rows, so the
  whole cell went, area and all.  On a 4 x 4 L with square cells of side 2,
  ten of thirty points had no cell and the cells covered 8 of the
  boundary's 10 square units.  Each collection is now reduced to its
  polygonal part, one row per cell still; a piece with no area is dropped
  as before.

* **`build_tessellation()` gave `NA` in `index` to a seed outside a boundary
  that is not convex, although its cell was built.**  A k-means centre is a
  mean of observations and can lie over a bay or a bend of a coast; its
  cell, the part of the study area nearer it than any other seed, holds
  ground and observations.  `index` gave each point the cell it lies in, and
  such a seed lies in none, so it was `NA` and the warning told the user to
  pass a boundary that holds it.  With `seeds$cell_id <- tess$index`,
  `cell_values_at(seeds = )` then put every location nearest that seed in
  no cell and filled it from the cell of its nearest scored row, which is
  not always the cell that holds it: on a notched square with one seed in
  the notch, 89 of 591 locations were filled, 6 of them from another cell.
  `build_tessellation(method = "voronoi")` and `create_voronoi_polygons()`
  now carry each cell's generator through the clip, and for seeds `index`
  holds the cell each one generates, so on that square no location is
  filled for want of a seed.  Seeds are a layer a seeding function returned
  (`rbind()` of such layers included), the seeds of a `select_resolution()`
  result, points that carry a record of the count and are its seeds (a
  `seed_id` column, or the record's centres), and points given with a
  `cover`.  Other points are indexed by the cell that holds them, as
  before, so an observation outside the boundary is still `NA`, also one
  that carries a selection's record as `assign_features_to_polygons()`
  returns it.  A seed is `NA` only when the clip removed its cell whole or
  its geometry is empty, and the warning then says that the points are
  seeds and why they have no cell.  The cells are unchanged.  The cell that
  holds each seed, which `index` gave before, is
  `assign_features_to_polygons(seeds, tess$cells, polygon_id_col =
  "cell_id", keep_unassigned = TRUE)$cell_id`.

* **`cell_values_at(seeds = )` filled the locations nearest a point placed
  by hand outside a boundary that is not convex.**  Points that are not
  seeds are indexed by the cell that holds them, so such a point used as a
  seed had `NA` in `index`, and the locations nearest it had no cell and
  were filled: the same 89 of 591 locations for the six points of the
  notched square.  A selection's seeds bound with `rbind()`, which drops the
  attributes that make them seeds, did the same.  When `values` carries the
  cell polygons, a seed with an `NA` ID that lies in none of them now takes
  the ID of the cell it generates, the polygon no seed names that a point
  inside it finds as its nearest seed, so with `index` alone every location
  of the notched square takes the cell that holds it and none is filled
  for want of a seed.  A seed with an `NA` ID inside a polygon keeps it, as
  does one whose cell the clip removed.  Without the polygons in `values`
  the locations nearest such a seed still have no cell;
  `get_voronoi_seeds(method = "provided", seeds = )` before the
  tessellation makes the points seeds, and `index` then gives each its
  cell.

* **For points placed by hand that share a cell in `index`, the error of
  `cell_values_at(seeds = )` told the caller to do what had been done, and
  `?cell_values_at` said such a point is `NA`.**  `build_tessellation()`
  gives a point that is not a seed and lies in no cell, but within a
  thousandth of the median cell width of one, the ID of the nearest cell.
  Beside a notch in the boundary that can be the cell of another point, so
  two points had one ID in `index`, and with `seeds$cell_id <- tess$index`
  the call stopped with "2 seeds have the cell of another", which advised
  setting the column from `index`.  The error now says why `index` repeats
  the ID, and that points passed through `get_voronoi_seeds(method =
  "provided", seeds = )` before `build_tessellation()` each take the cell
  they generate.  `?cell_values_at` and `?create_voronoi_polygons` now say
  such a point takes the ID of the nearest cell.  The call still stops: two
  seeds at different locations cannot share a cell.

* **`create_voronoi_polygons()` counted the points and seeds in no cell in
  the name of `build_tessellation()`, and `build_tessellation()` counted
  the `cover` points in no cell in the name of `create_voronoi_polygons()`.**
  Called directly, the `spatialkit_unassigned` warnings of
  `create_voronoi_polygons()` on points and seeds began
  "build_tessellation():", a function the caller had not called, as did the
  log line on points and `cover` points assigned to the nearest cell.
  Through `build_tessellation(method = "voronoi")` the warning on `cover`
  points outside every cell began "create_voronoi_polygons():", where the
  grid and triangle methods name `build_tessellation()`.  Each now names
  the function called.

* **`voronoi_seeds_kmeans()` gave the distinct positions as the reason for
  a clamp set by the number of points.**  Asked for as many seeds as there
  are points, each at its own position, it returns one fewer, since
  `stats::kmeans()` refuses as many centres as rows, and warned "requested 10
  seeds but only 10 unique positions among 10 point(s); clamping to 9".  The
  warning now takes the words of `get_voronoi_seeds()` for the same clamp:
  "the points support at most 9 (10 unique position(s) among 10
  point(s))".  The help says `k` is clamped to the distinct positions or
  one fewer than the points with a location, whichever is smaller; `k`
  said `nrow(points_sf) - 1`, which counts the rows without a location
  that are dropped first.  A count that is not whole, such as 9.5, is now
  given as requested in this warning and in the same warning of
  `get_voronoi_seeds()`, which said "requested 9 seeds" for a clamp to 9.

* **`cell_values_at(seeds = )` filled every location, and raised nothing,
  when no seed had the ID of a cell of `values`.**  Seeds whose ID column is
  not the one the cells were built with (a column of that name from
  elsewhere, or the `index` of another tessellation) gave every location a
  cell that `values` has no row for, and each location took the value of
  the cell of its nearest scored row: on a layer of 107 cells, seeds
  numbered 1001 to 1107 filled all 120 locations asked for, and 18 of the
  values differed from those of the locations' own cells.  Seeds whose IDs
  are all `NA` did the same unless a seed outside the cell polygons of
  `values` took the cell it generates, and a `values` with no cell ID (no
  row, or every ID `NA`) gave every location the value `NA`.  Each is now
  an error that names the argument to set: the seeds' ID column, with how
  to set it from the tessellation, or `values`.  The cells that seeds
  outside the polygons take count as their IDs, so points placed by hand
  outside the boundary, all `NA` in `index`, are still the seeds of their
  cells when `values` carries the polygons.  One seed with a location and
  the ID of a cell of `values` is enough for the call to run, and a seed
  whose ID `values` has no row for still stands for a cell without a value.

* **`cell_values_at(seeds = )` filled every location, with no error, when
  the seeds with the ID of a cell of `values` had no location.**  The check
  that some seed has the ID of a cell of `values` counted seeds with an
  empty point, or a coordinate that is not finite or lies beyond 1e150,
  though such a seed is the nearest seed of no location, and it read the
  coordinates only as given, so seeds that leave the globe once placed in
  the working CRS (metres given EPSG:4326, a latitude beyond 90) counted
  too.  Seeds that kept their IDs and lost their points left every location
  without a cell, and so did one such seed with a cell's ID beside seeds
  with `NA` IDs; beside seeds with other IDs it gave every location a cell
  that `values` has no row for.  Each location then took the value of the
  cell of its nearest scored row: the seeds of a tessellation of 107 cells,
  every point emptied, filled all 120 locations asked for.  Only a seed with
  a location, as given and in the working CRS, now counts, and when every
  seed with the ID of a cell has none the call stops with an error that
  says so; when they had a location as given, it names the working CRS and
  the CRS of `seeds`.  A seed with a coordinate beyond 1e150 is no longer
  searched where a move between CRSs puts it (a northing of 1e200 in
  EPSG:32618 comes back finite in EPSG:32617).  One seed with a location
  and the ID of a cell is still enough for the call to run.  A cell whose
  ID only seeds without a location carry is now one no seed names, which
  step 3 offers to a seed with an `NA` ID outside the polygons: points
  placed by hand outside the boundary take their cells beside such a seed,
  where the locations in the cell it named were filled.

* **A row of `newdata` or `points` with a coordinate beyond 1e150 in
  another CRS than the working one was read where the move to that CRS put
  it.**  Such a coordinate stands for no location, but the bound was read
  only after the move, and a move can bring it back finite: a northing of
  1e200 in EPSG:32616 lands at (87304, 5766268) in EPSG:32617, and one in
  EPSG:3857 at (5e5, 9997965).  In `cell_values_at()` such a row of
  `newdata` took a cell, a value and a distance, where the help says it has
  `NA` in every column; a row of `points` was a scored location, which
  moved the spacing and so `dist_ratio` at every location.  In
  `summarize_by_cell(estimate = "pooled")` such a row of `points` took part
  in the pooled columns and their record with no
  `spatialkit_pool_rows_skipped` warning.  Such a row of `points` also took
  part in choosing how the other rows are read: in EPSG:4326 a longitude of
  1e200 made `summarize_by_cell()` without a count record stop with an sf
  "invalid crs" error, and a latitude of 1e200 raised two sf warnings about
  the bounding box; in a layer without a CRS whose coordinates are lon/lat,
  one such row made every other row read as planar degrees, so that
  `summarize_by_cell()` recorded no CRS and `cell_values_at()` took the
  degrees for metres in the summary's CRS, with no warning.  The bound is
  now read on the coordinates as given, and again after the move, and the
  rows of `points` beyond it as given are emptied before the CRS is chosen:
  each such row gives what the same row with an empty point gives.  Rows in
  the working CRS itself are read as before.  The text of `spatialkit_pool_rows_skipped`, and its
  entry in `?summarize_by_cell`, said such rows "have empty or non-finite
  coordinates", which a row beyond 1e150 does not; they now read "no
  location (an empty point, or a coordinate that is not finite or lies
  beyond 1e150)".

* **A layer of `newdata`, `seeds` or `points` that mixes XY and XYZ points
  stopped `cell_values_at()` and `summarize_by_cell(estimate = "pooled")`
  with a bare message from base R or sf, or was read at the wrong
  coordinates.**  sf reads a layer in the dimension of its first point,
  and `rbind()` of an XY and an XYZ layer builds such a layer without a
  word.  With an XYZ point first, `sf::st_coordinates()` laid the values
  out of step, so that every later row was read at another row's place,
  with only base R's "data length is not a sub-multiple" warning: locations
  took another cell and distance, seeds another cell, and the pooled
  columns of every cell moved; in another CRS than the working one the call
  stopped with sf's "OGR error".  With an XY point first the call stopped
  with base R's "length of 'dimnames' [2] not equal to array extent".
  The same held for a layer XYZ throughout whose features are not all
  points: the point taken on the surface of a POLYGON Z, and the empty
  point of an empty LINESTRING Z, are XY beside the POINT Z rows kept.
  These layers are now read in x and y: when a layer's points, or the
  points taken from its other features, do not all have the same
  dimensions, Z and M are dropped from the copy read, and the call gives
  what the same rows in XY give.  An XYM point is read as an XYZ point is.
  A layer of points in one dimension, XY or XYZ throughout, is read as
  before.  Neither function returns these layers' geometry, so the
  caller's points keep their Z.

* **A layer that mixed XY and XYZ (or XYM) geometries was joined, measured
  and projected wrongly, depending on which row came first.**  sf reads a
  layer in the dimension of its first geometry: s2 drops Z and M only when
  the first geometry carries them, and GEOS and `sf::st_transform()` read
  an XYZ ring among XY rows as XY pairs.  `rbind()` of an XY and an XYZ
  layer gives such a layer without a word.  With a clean XYZ square in row
  3 of six lon/lat squares, `assign_features_to_polygons()` left the point
  in it unassigned under the warning that it fell inside no polygon and
  counted three others as ties; an XYZ square whose ring repeats a vertex
  in x and y passed the validity test and misplaced its points, where the
  same square alone gets the `spatialkit_cells_degenerate` error; polygon
  features that mixed dimensions went to the wrong cell under the
  largest-overlap join; in a projected CRS the join stopped with GEOS's
  "Points of LinearRing do not form a closed linestring".
  `summarize_by_cell(area = TRUE)` gave the XYZ square 1.6e14 square metres
  for its 95 square km, and refused mixed projected cells as not
  measurable.  `ensure_projected()` chose UTM zone 25 South for a layer in
  Missouri and returned the XYZ row as a ring of other points, and stopped
  with "OGR error" on points with an XYZ point first; `harmonize_crs()`
  did the same; both took such a layer in degrees without a CRS, a degree
  or less across, for planar coordinates, as their lon/lat test read its
  coordinates with `sf::st_coordinates()`, which stops on it;
  `ensure_projected(purpose = "area")` skipped its area distortion warning
  for such a layer in a projected CRS, as the distortion measured on it
  was `NA`; `create_grid_polygons()`, `create_voronoi_polygons()`,
  `build_tessellation()` and `clip_target_for()` stopped with GEOS's
  message on such a `boundary`.  Each geometry is now tested, and a layer
  that mixes dimensions is handed to sf without Z and M: both layers of
  the join, the cells whose areas and distortion `summarize_by_cell()`
  measures, the layer `ensure_projected()` chooses a projection from, and
  the builders' `boundary`.  The rows the join returns and the cells the
  summary returns keep the geometry they arrived with, Z and M included.
  `ensure_projected()` and `harmonize_crs()` return a transformed copy of a
  mixed layer in XY, with a warning of class `spatialkit_zm_dropped`, and
  so does the join for features without a CRS that are moved from lon/lat.
  A layer in one dimension, XY or XYZ throughout, is read and transformed
  as before, and so is a layer whose only geometries of another dimension
  than the first are empty and not points, which sf transforms with its Z.

* **A polygon layer without a CRS, in degrees and a degree or less
  across, was taken as planar when it held an empty polygon.**  The lon/lat test of
  `ensure_projected()` and `harmonize_crs()` reads the coordinates of such
  a layer to see whether they have the precision of decimal degrees;
  `sf::st_coordinates()` stops on a polygon layer with an empty polygon,
  so the test found no precision and the layer was not moved from lon/lat.
  The test now leaves empty geometries out.

* **`coerce_to_points()` stopped with a bare message from GEOS or sf, or
  returned points at the wrong place, on a polygon or line layer that
  mixed XY and XYZ (or XYM) geometries, and so did every function that
  reduces the caller's layer to points first.**  sf hands GEOS and s2 every
  geometry of a layer in the dimension of its first, and `rbind()` of an
  XY and an XYZ layer gives such a layer without a word.  On projected
  polygons the modes `"auto"`, `"centroid"` and `"point_on_surface"`
  stopped with GEOS's "ParseException: Unexpected EOF parsing WKB" or
  "Points of LinearRing do not form a closed linestring"; on lines
  `"auto"` and `"line_midpoint"` stopped with sf's "OGR error" when the
  first line was XYZ, and gave a later XYZ line another midpoint in
  silence; on squares in degrees `"centroid"` put an XYZ square's point
  8,666 km from the square.  `prep_model_data()`, and through it the fits,
  the cross-validations and `predict()` on `newdata` in the fit's CRS,
  `make_folds()`,
  `estimate_sac_range()`, `compare_models_cv()` and
  `select_features_forward(select_on = "split")` stopped the same way.
  `make_folds()`, `compare_models_cv()` and
  `select_features_forward(select_on = "split")` also stopped with "OGR
  error" on a layer XYZ throughout that mixes POINT Z rows with POLYGON Z
  or LINESTRING Z rows, as the points taken from the polygons and lines
  are XY, and `make_folds()` stopped the same way on a point layer with an
  XYZ point first among XY points.
  `coerce_to_points()` now reads such a layer in x and y and returns XY
  points, with a warning of class `spatialkit_zm_dropped` that carries the
  rows with Z or M as `rows`.  This includes a layer whose geometry types
  each have one dimension, such as XY polygons beside a POINT Z row: it
  was read right, and its points came back with Z on that row only.
  `prep_model_data()`, `make_folds()`, `estimate_sac_range()`,
  `compare_models_cv()` and `select_features_forward(select_on = "split")`
  read the layer in x and y with no warning, and give what the same rows
  in XY give.  A layer in one dimension, XY or XYZ throughout, is
  read as before, and so is a layer whose only geometries of another
  dimension than the first are empty and not points.  The mode
  `"bbox_center"` read each geometry's box and gave the right points; it
  still does, with no warning.

* **`coerce_to_points(mode = "auto")` gave a MULTILINESTRING whose longest
  parts are of equal length the midpoint of the first of them in the
  feature's part order**, so the same geometry with its parts listed the
  other way round got a point 1 km away in a two-part example; congruent
  parts were also told apart by the last bits of their lengths, so the
  point moved with a shift of the origin.  Among parts of equal greatest
  length, to one part in ten million, the midpoint lowest in x, then in y,
  is now taken, x values within one part in ten million of the diagonal of
  the midpoints' bounding box counting as one: two congruent parts whose
  midpoints have one x, reached from other coordinates, differ in their
  last bits once the origin moves, and y now decides between them, so the
  point is the same at every origin.  A feature with one longest part gets
  its midpoint as before.
  `?coerce_to_points` states the rule.

* **`summarize_by_cell(deff = "variogram")`, `kriging_adequacy()` and
  `cell_values_at()` stopped with a bare message from sf or GEOS, or read
  other coordinates, on layers that mixed XY and XYZ (or XYM)
  geometries.**  sf reads a layer in the dimension of its first geometry.
  The variogram design effect of `summarize_by_cell()` stopped with "OGR
  error" on points with an XYZ point first when `sac` records its CRS;
  when it records none, it read other coordinates with base R's "data
  length is not a sub-multiple" warning as the only sign, or stopped with
  base R's "length of 'dimnames' [2] not equal to array extent".
  `kriging_adequacy()` stopped or read other coordinates the same ways on
  such points; with `sac = NULL` the other coordinates came with a warning
  that points share a location, though none did.  It stopped on such
  `cells_sf` with GEOS's "ParseException: Unexpected EOF parsing WKB",
  "Points of LinearRing do not form a closed linestring" or "GEOS does not
  support XYM", or with "OGR error" for cells in degrees.
  `cell_values_at()` with `seeds` whose ID is `NA` stopped with "OGR error"
  or the LinearRing message on such polygons in `values`.  These layers,
  and the points `coerce_to_points()` takes from such polygon features or
  multipoints, are now read in x and y, with no warning, as none of the
  three functions returns them, and give what the same rows in XY give.
  A layer XYZ throughout that puts POINT Z rows beside MULTIPOINT Z or
  POLYGON Z rows stopped with "OGR error" or read other coordinates, as the
  points taken from its multipoints and polygons are XY; it is now read in
  x and y in the same way.  Any other layer in one dimension, XY or XYZ
  throughout, is read as before.

* **The seeding functions, `ensure_stable_poly_id()` and
  `create_grid_polygons_cached()` stopped with a bare message from base R,
  sf, GEOS or s2, or gave other seeds in silence, on layers that mixed XY
  and XYZ (or XYM) geometries.**  sf reads a layer in the dimension of its
  first geometry.  `get_voronoi_seeds(method = "kmeans")` and
  `voronoi_seeds_kmeans()` clustered `sample_points` or `points_sf` with
  an XYZ point first on other coordinates, with base R's "data length is
  not a sub-multiple" warning as the only sign, and stopped with base R's
  "length of 'dimnames' [2] not equal to array extent" when the XYZ point
  came later; in degrees they stopped with "OGR error" or gave other seeds.
  With a profile or a selection as `n`, `get_voronoi_seeds()` stopped with
  GEOS's "ParseException: Unexpected EOF parsing WKB" or the `dimnames`
  message.  `voronoi_seeds_random()` and `get_voronoi_seeds()` stopped on
  such a `boundary` with GEOS's message ("ParseException", "Points of
  LinearRing do not form a closed linestring" or "GEOS does not support
  XYM") or s2's "Loop 0 is not valid"; in degrees an XYZ polygon after the
  first gave random seeds that fell outside the boundary, with no warning.
  `ensure_stable_poly_id()` stopped on such polygons with GEOS's message,
  sf's "OGR error" or wk's parse error, or keyed them under a false warning
  that s2 refused most of them.  `create_grid_polygons_cached()` stopped
  with GEOS's message on such a `boundary` while it built the cache key.
  These layers are now read in x and y, before any transform, union or
  repair, and give what the same rows in XY give.  The seeders return new
  points and the cached grid is new geometry, so they say nothing; a mixed
  boundary shares its cache key with the same boundary in XY.
  `ensure_stable_poly_id(make_valid = FALSE)` returns the caller's rows
  with their Z, keyed in XY; with `make_valid = TRUE` (the default) it
  returns its repaired copy, which for such polygons is made in XY, with a
  warning of class `spatialkit_zm_dropped`.  A layer in one dimension, XY
  or XYZ throughout, is read as before, save that
  `ensure_stable_poly_id(make_valid = TRUE)` now repairs polygons that all
  carry Z or M in XY with a warning (see below).

* **`get_voronoi_seeds(method = "provided")` stopped with sf's "OGR
  error" on `seeds` in another CRS than `boundary` whose first point had a
  Z or M coordinate among XY points.**  The seeds are returned in the
  boundary's CRS, and `sf::st_transform()` reads a layer in the dimension
  of its first geometry.  Such seeds are now returned in the boundary's
  CRS in XY, with a warning of class `spatialkit_zm_dropped` that carries
  the rows with Z or M as `rows`.  Seeds in the boundary's CRS, or given
  without a boundary, come back as they were given, Z included.

* **`create_voronoi_polygons()`, `build_tessellation()` and
  `clip_target_for()` stopped with a bare message from base R, sf or GEOS
  on `points_sf` that mixed XY and XYZ (or XYM) points, and set a mixed
  `cover` aside under a false warning.**  sf reads a layer in the dimension
  of its first geometry.  With an XYZ point first among XY points,
  `create_voronoi_polygons()` stopped with GEOS's "ParseException:
  Unexpected EOF parsing WKB" (or "GEOS does not support XYM"), and
  `build_tessellation()` with sf's "OGR error" for every method; with the
  XYZ point later, `create_voronoi_polygons()` and `method = "triangles"`
  stopped with base R's "length of 'dimnames' [2] not equal to array
  extent".  Points in degrees ran, but raised `ensure_projected()`'s
  `spatialkit_zm_dropped` warning about an XY copy none of the three
  functions returns.  A polygon `cover` that mixed dimensions, given as the
  argument or carried by the seeds, raised `spatialkit_no_cover` ("gives no
  polygon in the working CRS", or "could not be put in the working CRS")
  and the cells were built on the points' own hull; a point `cover` with an
  XYZ point first, checked against cells built in a `boundary`, stopped
  with "OGR error".  The points and the cover are now read in x and y,
  before they are moved or read, with no warning, as the cells and the
  clip target are new geometry, and give what the same rows in XY give.
  The seeds keep their attributes, their cover and count record among
  them.  A layer XY or XYZ throughout is read as before.

* **The tessellation builders, the grid builders, the seeders and
  `clip_target_for()` stopped with GEOS's "GEOS does not support XYM or
  XYZM geometries" on a projected `boundary` whose polygons all carry M,
  where the same boundary with one XY polygon among them ran.**  In
  degrees the grids, the builders and `clip_target_for()` stopped the same
  way, and the seeders ran but attached the cover with M.
  `get_voronoi_seeds(method = "kmeans")` given `sample_points` ran in
  either CRS, and a builder given its seeds set their cover aside under a
  false `spatialkit_no_cover` warning or stopped.  The cells and seeds
  are defined in x and y and M is read nowhere.  Such a boundary is now
  read in x and y by `create_grid_polygons()`,
  `create_grid_polygons_cached()`, `create_voronoi_polygons()`,
  `build_tessellation()`, `clip_target_for()`, `get_voronoi_seeds()` and
  `voronoi_seeds_random()`, with no warning, and gives what the same
  boundary in XY gives; `create_grid_polygons_cached()` keeps one cache
  entry for the two, and the seeders' cover is the boundary in XY.  So are
  `points_sf` and `cover` of `create_voronoi_polygons()` and
  `build_tessellation()` (and `points_sf` of `clip_target_for()`) that
  carry M: points XYM throughout stopped the Voronoi builder and the grids
  with the same message, and a cover XYM throughout was set aside with the
  false warning above.  A layer XY or XYZ throughout is read as before;
  the cover the seeders attach to a boundary XYZ throughout keeps its Z.

* **`ensure_stable_poly_id(make_valid = TRUE)`, the default, returned
  polygons that all carry Z in XY without a word, and stopped with GEOS's
  "GEOS does not support XYM or XYZM geometries" on projected polygons
  that all carry M (in degrees they came back in XY without a word).**  The rows returned are the repaired polygons, and
  `sf::st_make_valid()` on GEOS 3.13 drops Z from every polygon it repairs
  and refuses M.  Only a layer that mixed XY and XYZ (or XYM) polygons
  raised a warning, which let a reader take a layer in one dimension to
  keep its Z.  Z and M are now dropped before the repair whenever any
  polygon carries them, so the rows returned are XY on every GEOS version,
  and a warning of class `spatialkit_zm_dropped` names `polygons_sf`, says
  how many of its geometries carry Z or M, and carries their rows as
  `rows`.  The same holds for an XYZ layer with an empty XY polygon after
  the first.  An empty polygon in Z or M has no coordinates to drop, so it
  is not counted, and XY polygons with such an empty one among them raise
  no warning.  Mixed layers keep their warning.  With `make_valid = FALSE`
  the rows come back with their Z and M, with no warning, as before.  The
  cells `create_voronoi_polygons()`, `build_tessellation()` and
  `create_grid_polygons_cached()` key are new geometry in XY and raise no
  such warning.

* **`resolution_profile()`, `determine_optimal_levels()`, the fold
  provenance check, `make_folds()`'s `boundary`, `blocks` and
  `prediction_points`, and `fold_separation()` stopped with a bare message
  from base R, sf or GEOS, read other coordinates, or let foreign folds
  through, on layers that mixed XY and XYZ (or XYM) geometries.**  sf
  reads a layer in the dimension of its first geometry.
  `resolution_profile()` stopped on points with an XYZ point first with
  GEOS's "ParseException: Unexpected EOF parsing WKB", after base R's
  "data length is not a sub-multiple" warning, and with base R's "length
  of 'dimnames' [2] not equal to array extent" when the XYZ point came
  later; `determine_optimal_levels()` clustered other coordinates with
  that base R warning as the only sign, or stopped with the `dimnames`
  message.  In degrees both stopped with sf's "OGR error" when the first
  point had Z.  The
  provenance check that `cv_rf()`, `cv_gwr()`, `cv_bayes()`,
  `cv_spatial()`, `compare_models_cv()`, `fold_separation()`,
  `kriging_adequacy()` and `area_of_applicability(folds =)` run on supplied
  folds could not take its sample of the layer when the first row had a Z
  or M coordinate, and passed folds built on other data with no word; on
  a layer without a CRS it read the coordinates out of step and refused
  folds built on the same rows as built from different data.
  `make_folds(method = "block_kfold")` stopped on such a `boundary` or
  `blocks` with sf's "OGR error" when the first polygon had Z or M, and
  with GEOS's "Points of LinearRing do not form a closed linestring" when
  a later one had, also through the `boundary` of the `cv_*()` functions.
  `make_folds(method = "nndm")` and `"knndm"` stopped with "OGR error" on
  `prediction_points` with a Z or M point first, and `fold_separation()`
  did the same on its `prediction_points`, and on `data_sf` when it moves
  it into the CRS of the range or projects it from degrees.  Polygon
  layers that mixed dimensions stopped these functions with GEOS's
  message.  These layers are now read in x and y, before they are reduced
  to points, moved or repaired, and again after they are reduced to
  points, with no warning, as none of these functions returns them; they
  give what the same rows in XY give.  The blocks a `make_folds()` result
  records as `params$blocks` are its working copy and come back in XY for
  such `blocks`.  The points a layer is reduced to are read in x and y
  too: those of a layer XYZ throughout that holds points beside polygons
  mix XY and XYZ, since the point on the surface of a POLYGON Z is XY.
  With a point first, `resolution_profile()` stopped on such a layer with
  "ParseException: Unexpected EOF parsing WKB", `determine_optimal_levels()`
  clustered other coordinates with base R's "data length" warning as the
  only sign, `make_folds(method = "nndm")`, `"knndm"` and
  `fold_separation()` stopped with "OGR error" on such `prediction_points`,
  and `fold_separation()` did the same on such `data_sf` when it moves it
  into the CRS of the range.  With a polygon first, `resolution_profile()`
  and `determine_optimal_levels()` stopped with the `dimnames` message.
  Either way `fold_separation()` raised `ensure_projected()`'s
  `spatialkit_zm_dropped` warning when it projects such `data_sf` from
  degrees.  Any other layer XY or XYZ throughout is read as before.

* **`resolution_profile()`, `determine_optimal_levels()`,
  `estimate_sac_range()`, `prep_model_data()`, `make_folds()` and
  `fold_separation()` stopped with GEOS's "GEOS does not support XYM or
  XYZM geometries" on layers whose geometries all carry M.**  The profile
  stopped on points or polygons XYM; `determine_optimal_levels()`,
  `estimate_sac_range()` and `prep_model_data()`, and through it the
  `cv_*()` functions, on polygons; `make_folds()` on such a `boundary`
  (also through the `boundary` of the `cv_*()` functions), such `blocks`,
  and polygon `points_sf` or `prediction_points`; `fold_separation()` on
  polygon `data_sf` or `prediction_points`.  These functions read these
  layers only for their x and y and do not return them, so the layers are
  now read in x and y, as the tessellation builders read a `boundary` that
  carries M, with no warning, and give what the same rows in XY give.  A
  layer XY or XYZ throughout is read as before.

* **`prep_model_data()` raised `ensure_projected()`'s
  `spatialkit_zm_dropped` warning, which says the layer is returned in XY,
  for a `boundary` in degrees that mixed XY and XYZ (or XYM) geometries.**
  Only the boundary's CRS is taken and the boundary is not returned.  Such
  a boundary is now read in x and y before it is projected, with no
  warning, and the data are projected as for the same boundary in XY.

* **`area_of_applicability()`, `predict_surface()`, `predict()`,
  `plot_tessellation_map()`, `plot_folds()`, `plot()` of a fit and
  `residual_morans_i()` stopped with a bare message from base R, sf or
  GEOS, or read other coordinates, on layers that mixed XY and XYZ (or
  XYM) geometries.**  sf reads a layer in the dimension of its first
  geometry.
  `area_of_applicability()`, wherever it measures location (a forest
  fitted with `include_coords = TRUE`, a GWR or Bayesian fit, or
  `coords = "yes"`), read the coordinates of such points in `newdata` or
  `train_sf` out of step, with base R's "data length is not a
  sub-multiple" warning as the only sign, when the first point had Z or
  M, and stopped with base R's "length of 'dimnames' [2] not equal to
  array extent" when a later one had; `newdata` in another CRS with a Z
  or M point first stopped it with sf's "OGR error".
  `predict_surface()` stopped with "OGR error" or "Points of LinearRing
  do not form a closed linestring" on such a `boundary`, with
  "ParseException: Unexpected EOF parsing WKB" (or "OGR error", given a
  `boundary`) on such `covariates` or a POINT `grid` with an XYZ point
  first, and with "GEOS does not support XYM or XYZM geometries" with an
  XYM point first, and raised
  `ensure_projected()`'s warning about a copy it does not return for
  such `covariates` in another CRS.  `predict()` on a forest, a GWR or a
  Bayesian fit, and through it `evaluate_insample()`, `compare_models()`
  and `model_metrics()`, stopped with "OGR error" on `newdata` in another
  CRS with an XYZ or XYM row first; once `ensure_projected()` read such a
  layer in x and y (above), they raised its warning about a copy they do
  not return.
  `plot_tessellation_map()` stopped with "OGR error" on any such layer
  with a Z or M row first, and drew a later XYZ cell or boundary far from
  its place, or stopped with GEOS's message when it labels the cells;
  `plot_folds()` stopped with "OGR error" when the plot was printed, or
  drew the boundary elsewhere, and so did `plot()` of a custom fit (built
  with `new_spatial_fit()`) on its `data_sf`.  A polygon or line layer
  that mixed dimensions, given to `area_of_applicability()` or as the
  `grid` of `predict_surface()`, and the layer of a custom fit given to
  `residual_morans_i()`, raised the warning `coerce_to_points()` gives
  about its `x`, an argument of a function the caller did not call.
  These layers are now read in x and y, and give what the same rows in XY
  give, with no warning, save where a copy of the caller's rows comes
  back: `newdata` that `area_of_applicability()` reprojects into the
  training CRS for a forest fitted with `include_coords = TRUE`, and a
  POINT `grid` that `predict_surface()` reprojects into the fit's CRS, are
  returned in XY with a warning of class `spatialkit_zm_dropped`, as
  `ensure_projected()` returns such a layer.  Rows returned in their own
  CRS, or without a CRS and stamped with the fit's, keep their geometry,
  Z and M included.  A layer XY or XYZ
  throughout is read as before.

* **`coerce_to_points()` stopped with GEOS's "GEOS does not support XYM or
  XYZM geometries" on geometries that carry M, and so did the functions
  that reduce the caller's polygons to points.**  In a projected CRS
  `"auto"`, `"centroid"` and `"point_on_surface"` stopped on polygons
  whose geometries all carry M (XYM or XYZM), `"auto"` on such
  multipoints, `"centroid"` and `"point_on_surface"` on such points and
  lines, and `"line_midpoint"` on such geometries that are not lines; in
  degrees `"auto"` stopped on such polygons and `"point_on_surface"` on
  all of them.  `"auto"` and `"line_midpoint"` gave lines XY midpoints
  without a word, and so did `"centroid"` in degrees for every geometry.
  `plot_tessellation_map(labels = TRUE)` stopped the same way on cells
  that carry M.
  `compare_models_cv()`, `select_features_forward(select_on = "split")`,
  `cell_values_at()` (polygon `newdata` or `points`),
  `summarize_by_cell()` with `deff = "variogram"` or `estimate =
  "pooled"`, `area_of_applicability()` where it measures location, and
  `predict_surface()` with a polygon `grid` stopped the same way on such
  polygons, `predict_surface()` also on such a `boundary`, `covariates`
  or POINT `grid`, `kriging_adequacy()` on such multipoints, and
  `residual_morans_i()` on a custom fit with such polygons returned
  `NULL` with "coordinate extraction failed".  The same layers
  XYZ throughout ran, as GEOS's points of a POLYGON Z are XY.  When the
  mode reduces a geometry that carries M, `coerce_to_points()` now drops
  Z and M from the layer and returns XY points with a warning of class
  `spatialkit_zm_dropped` that carries the rows with M as `rows`; empty
  geometries are not counted, and when they are the only ones the mode
  reduces that carry M, Z and M are dropped from them alone, without a
  warning.  POINT rows under `"auto"` are returned as they are, with
  their M, and `"bbox_center"` is unchanged.  The
  functions above read such a layer in x and y before they reduce it, with
  no warning, and give what the same rows in XY give; a POINT `grid` of
  `predict_surface()` in the fit's CRS keeps its M.  The labels of
  `plot_tessellation_map()` are taken in x and y.  A layer XY or XYZ
  throughout is read as before.

* **`assign_features_to_polygons()`, `cell_values_at()` and
  `kriging_adequacy()` stopped with GEOS's "GEOS does not support XYM or
  XYZM geometries" on layers whose geometries all carry M.**
  `assign_features_to_polygons()` stopped on such `features_sf`, points or
  polygons, and on such `polygons_sf`, with `largest = TRUE` or `FALSE`,
  wherever the join runs with GEOS (in a projected CRS, without a CRS, or
  in degrees with s2 off); in degrees on the sphere s2 drops M and the join
  ran.  `cell_values_at()` stopped on such
  `values` without `seeds` or beside a seed with an `NA` ID, and on point
  `newdata` without `seeds`; `kriging_adequacy()` on such `cells_sf`.  The
  same layers XYZ throughout ran, and a layer that mixed XY and XYM rows
  was already joined in x and y.  The join and the
  cells are defined in x and y and M is read nowhere, so these layers are
  now read in x and y, with no warning, and give what the same rows in XY
  give.  The rows `assign_features_to_polygons()` returns keep the
  geometry they arrived with, M included; features without a CRS in one
  dimension that are moved from lon/lat into the polygons' CRS come back
  moved with `sf::st_transform()`, which keeps their M.  A layer XY or XYZ
  throughout is read as before.

* **The smallest-area tie-break of `assign_features_to_polygons()` gave a
  feature to a polygon whose ring crosses itself as the smallest.**  GEOS
  measures such a ring as the signed sum of its lobes, which cancel: a
  bowtie whose two lobes cover 180,000 square metres measured 0, so a point
  inside a valid 100 m square and inside the bowtie went to the bowtie,
  without a warning, where the same bowtie repaired by the caller, or a
  valid polygon in its place, lost the tie.  The largest-overlap join
  already repaired such polygons for its overlaps.  The tie-break now
  checks the polygons it compares and measures the invalid ones on copies
  repaired with `sf::st_make_valid()`, with a warning in the words of the
  largest-overlap join's; the polygons and features are returned as they
  arrived.  Layers without ties, or whose tied polygons are valid, are
  assigned as before.

* **`assign_features_to_polygons()` on layers in longitude and latitude
  stopped with sf's "package lwgeom required, please install it first"
  when `sf::sf_use_s2()` was `FALSE`.**  sf hands areas in degrees to
  lwgeom when s2 is off, and the package does not depend on lwgeom: the
  smallest-area tie-break stopped on a point in two polygons, and the
  largest-overlap join of polygon features stopped with "the largest-overlap
  join failed (package lwgeom required ...). Pass largest = FALSE ...",
  whose advice stopped the same way.  In degrees the largest-overlap join
  (its validity check, repair and overlaps) and the areas the tie-break
  compares now run on s2 whatever the session says, as the package's other
  areas in degrees do, and give what they give with s2 on; the join by
  `predicate` keeps the session's engine.  The largest-overlap join no
  longer advises `largest = FALSE` for an error that the join by
  `predicate` meets too, such as a missing lwgeom.  With s2 off, a polygon
  in a tie that s2 refuses even repaired (a bowtie in degrees) stops the
  call with the `spatialkit_cells_degenerate` error, which names it and the
  transform.  With s2 off, the error of the largest-overlap join for a
  polygon s2 refuses offers `largest = FALSE`, whose join by `predicate`
  runs with GEOS, and these errors name s2 as the engine the package uses
  there, not as the one sf uses in the session.

* **`cell_values_at()` without `seeds` gave every location the value `NA`,
  marked filled, when `values` had no cell ID.**  A `values` whose ID
  column is all `NA` placed no location in a cell, so every location was
  filled from the cell of its nearest scored row, which `values` has no row
  for, and took `NA`; nothing was raised.  The same `values` beside `seeds`
  was refused.  It is now refused on both paths with the error the seeds'
  path gave ("`values` has no cell ID (it has no row, or its ID column ...
  is all NA)"), whatever the locations asked for.  IDs that are `NaN`
  count as `NA`.  A `values` of no rows, or one without polygons whose IDs
  are all `NA`, without `seeds`, refused before because it carries no
  polygons, gets that error too.  A `values` with some IDs `NA` runs as
  before.

* **`cell_values_at(seeds = )` gave a seed with an `NA` ID outside the
  polygons of `values` the first of several cells by the row order of
  `values`.**  Such a seed takes the cell no seed names whose inner point
  finds it as its nearest seed.  When the inner points of several such
  cells found it, it took the first of them in the order of the rows of
  `values`, and every location nearest it took that cell's value, not
  filled: one point placed by hand outside the boundary, passed alone as
  `seeds`, put all 120 locations of a layer of 107 cells in cell 1, and in
  cell 107 with the rows of `values` reversed.  Such a seed generates at
  most one of those cells, and nothing says which, so it now keeps its
  `NA` ID: the locations nearest it have no cell and are filled, and a
  warning gives the count of such seeds.  When no seed is then left with
  the ID of a cell of `values`, as for that lone point, the call stops with
  the error for seeds whose IDs are all `NA`; its text, with the polygons
  in `values` and without, now says when a seed outside the polygons is
  given its cell.  A seed that the inner point
  of one such cell alone finds, such as a point placed by hand in the bend
  of a boundary that is not convex, takes that cell as before.

* **`cell_values_at()` refused `seeds` without a CRS beside `points` with
  one, with an error that did not say which CRS to set.**  For rows without
  a CRS whose coordinates look like longitude and latitude,
  `get_voronoi_seeds()` returns the seeds without a CRS, in degrees, and
  `assign_features_to_polygons()` returns the rows in the cells' projected
  CRS, so on such rows `cell_values_at()` stopped on the seeds even once
  `newdata` was given its CRS.  Such seeds are now read as
  `build_tessellation()` reads a boundary without a CRS, with a warning:
  beside `points` in a projected CRS, coordinates that look like longitude
  and latitude are taken as EPSG:4326 and transformed, and others are taken
  to be in the CRS of `points`, which the warning names; beside `points` in
  a geographic CRS, coordinates within the ranges of longitude and latitude
  are taken to be in it, and others are an error that names the CRS the
  summary measured in.  `newdata` without a CRS beside `points` with one is
  still an error.

* **`assign_features_to_polygons()` stopped on cells returned in lon/lat
  with "duplicate vertex".**  `build_tessellation(crs = 4326)` builds the
  cells in a projected CRS and transforms them.  On gridded locations four
  centres can stand the same distance from one point, and two cells meet
  there at two corners that differ only by rounding, less than a nanometre
  apart; the transform puts them at one point, or at two that sf's
  spherical engine (s2) cannot tell apart, and s2 refused the cell.  Rows
  in lon/lat were joined on the sphere, so the join stopped: on 16 lattices
  at spacings of 1 to 1,000 m with rows in EPSG:4326, 27 of 32 paths
  stopped under the two kinds, and 3 of the 5 others put rows on an edge in
  another cell than the profile's.  Cells in lon/lat whose record says a
  profile scored them are now joined in the record's projected CRS, where
  they were drawn: every row is in the cell the profile scored it in, a
  point on an edge goes to the cell of its nearest centre, and the
  summary's check reads `match` on all 32.  This holds for the cells s2
  accepts too, so features that are not points get their cell from
  overlaps and areas measured in that CRS, and rows without a CRS, which
  also stopped, are joined there.  `cell_values_at()` without `seeds`
  places locations on such cells there as well.  Other cells in lon/lat are
  joined as before: cells without the record (from a count given as a
  number, or after a verb that drops attributes) and cells whose record
  says they are not the cells the profile scored (`scored_cells = FALSE`,
  such as the cells of centres refitted to every point).  s2 can still
  refuse them, as it refuses `sf::st_area()` and `sf::st_union()` on such
  cells: transform them to a projected CRS first.

* **`create_voronoi_polygons()` and `build_tessellation(method = "voronoi")`
  handed on a seeds layer's record of its count unchanged when the seeds
  were fewer or more than that count.**  The record survives row
  subsetting, so a selection's 32 seeds cut down to 16 gave 16 cells whose
  record said they were the 32 cells the profile scored, and
  `summarize_by_cell()` warned that rows were dropped or moved between
  cells, by points outside every cell, a clipping boundary, a subset of the
  rows or another layer, when every row was assigned.  Such seeds now carry
  the record with `scored_cells = FALSE`, as `get_voronoi_seeds(method =
  "provided")` marks them, so their cells are built as other seeds' are and
  the summary reads "not checked: other cells were built from this count".
  A selection's seeds as they come are unchanged.

* **`assign_features_to_polygons()` stopped with s2's bare message on
  polygons in lon/lat that s2 refuses.**  A join in longitude and latitude
  runs on the sphere, where s2 refuses a polygon with two corners it reads
  as one point, and the join stopped with "Loop 0 is not valid: Edge 20 is
  degenerate (duplicate vertex)", which names neither the polygon nor a way
  on.  Cells returned in lon/lat hold such corners on gridded locations, and
  the joins that do not run in the record's CRS met them: s2 refuses 4 of
  the 133 k-means cells of a 30 by 20 lattice at 30 m returned in
  EPSG:4326, and the join of rows in degrees stopped once the cells had lost
  their record or their CRS; the cells of centres refitted to every point
  (`scored_cells = FALSE`) on a 50 by 40 lattice at 250 m stopped it the
  same way.  The join by `predicate` now stops with an error of class
  `spatialkit_cells_degenerate` that gives the number and rows of the
  polygons s2 refuses (`rows` on the condition), quotes s2, names the cause
  and says to transform `polygons_sf` to a projected CRS, such as the one
  the cells were built in, where the join runs.  Every refusal by s2 of a
  polygon of the layer gets this error, among them a hole that crosses its
  outer ring and two holes that share an edge, which also stopped with s2's
  bare message.  The largest-overlap join of polygon features still repairs
  the polygons s2 refuses, with its warning, and runs on most of them; on
  a polygon s2 refuses even after the repair (a hole that crosses its outer
  ring) it stops with the same error, where it advised `largest = FALSE`,
  which meets the same refusal.  When the repair makes s2 accept a polygon
  and s2 then refuses a piece of its overlap with the other layer (a figure eight,
  repaired into two triangles that touch at a point), the join gave the
  same advice, and with s2 on the join by `predicate` then stopped on the
  polygon as it arrived; it now stops with the error that names that
  polygon and the transform.  Whether s2 refuses such a piece follows the
  last bits of the vertices s2 computes for it, which differ between
  platforms and s2 versions: a figure eight pinched on the edge of two
  cells can give a piece that s2 refuses on macOS on arm64 and accepts on
  Linux on x86-64, and where s2 accepts the piece the join runs and gives
  each feature the polygon of its largest overlap.  Where s2 refuses a
  polygon of either layer even after the repair, the error names that
  polygon, and not a polygon of the other layer that the repair made s2
  accept.  A polygon of `features_sf` that s2 refuses stopped the join by
  `predicate` with s2's bare message, such as "Loop 0 is not valid: Edge 0
  crosses edge 2" for a bowtie; it now stops with an error of class
  `spatialkit_features_degenerate`, with the fields of the cells' class,
  that gives its rows, quotes s2 and says to transform `features_sf` to a
  projected CRS, and, where `sf::st_make_valid()` makes s2 accept the rows,
  that repairing them works too.  Where s2 refuses polygons of both layers,
  the error for `polygons_sf` quotes each row's own reason, where it quoted
  s2's message for a feature as a cell's, and says that s2 also refuses rows
  of `features_sf` (`features_rows` on the condition).
  `cell_values_at()` without `seeds`, which places its locations with this
  join, gives the error in its own words: it names
  `values` and the rows of `values`.  Polygons s2 accepts are joined as
  before, and an error of another kind keeps its own message.

* **`summarize_by_cell(area = TRUE)` stopped without the cause on lon/lat
  cells that s2 refuses.**  The areas of lon/lat cells are s2's, and s2
  refuses a cell with two corners it reads as one point, which cells
  returned in lon/lat hold on gridded locations.  When the area check
  measured such a cell, the request was refused because "the area
  distortion of `cells_sf`'s CRS could not be measured", which names
  neither the cell nor a way on; above 200 cells, when the check passed
  over it, the request ran through the summary and stopped at the cells'
  areas with s2's bare "Loop 0 is not valid: Edge 1 is degenerate
  (duplicate vertex)".  `?summarize_by_cell` said the check passes by
  construction on lon/lat cells.  Every cell's area is now taken before
  any summary, and a cell s2 refuses stops the request with the error of
  class `spatialkit_cells_degenerate` that `assign_features_to_polygons()`
  gives: it gives the rows of `cells_sf` s2 refuses, names the cause and
  says to transform the cells to a projected CRS, such as the one they
  were built in or the equal-area one `ensure_projected(purpose = "area")`
  picks.  Lon/lat cells s2 accepts get the density and
  `attr(, "area_error")` they got before, and with s2 off the request
  needs lwgeom as before.

* **`coerce_to_points(mode = "centroid")` and `prep_model_data(pointize =
  "centroid")` stopped with s2's bare message on a polygon in longitude and
  latitude that s2 refuses.**  With s2 on, sf takes the centroids of a layer
  in degrees with s2, which refuses a ring with a vertex listed twice in a
  row, although GEOS accepts it, and the call stopped with "Loop 0 is not
  valid: Edge 1 is degenerate (duplicate vertex)", naming neither the row
  nor a way on; a bowtie and a figure eight stopped it the same way.  The
  same held for the rows other than LINESTRINGs under `"line_midpoint"`,
  for collections under `"auto"`, for the `cv_*()` and `fit_*()`
  functions that pass `pointize = "centroid"`, and for the functions that
  reduce a layer they are given to points with `"auto"`: `make_folds()`,
  `fold_separation()`, `estimate_sac_range()`, `resolution_profile()`,
  `determine_optimal_levels()`, `cell_values_at()`, `summarize_by_cell()`
  and `select_features_forward()` (the last wrapped the message in a plain
  error of its own).  A polygon s2 refuses for a
  repeated vertex, on its own or inside a collection, now gets the centroid
  of the same ring without the repeat, whose points are unchanged.  A
  polygon s2 refuses even then stops
  the call with an error of class `spatialkit_features_degenerate` that
  gives its rows (`rows` on the condition), quotes s2, and says to
  transform the layer to a projected CRS, or to pass another `mode` or
  `pointize`, or, in a function that takes neither, to reduce the layer to
  points first.  The `mode` or `pointize` it offers runs on the rows it
  names: `"auto"` when each is a polygon or multipolygon, which `"auto"`
  reduces with a point on the surface, and `"point_on_surface"` when one
  is a collection, which `"auto"` takes the centroid of and refuses
  again.  Each function raises it in its own name, naming its
  argument and that argument's rows: `prep_model_data()`,
  `compare_models_cv()`, `resolution_profile()` and
  `determine_optimal_levels()` name `data_sf`, `make_folds()` `points_sf`
  or `prediction_points`, `fold_separation()` `data_sf` or
  `prediction_points`, `estimate_sac_range()` `points_sf`,
  `cell_values_at()` `newdata` or `points`, `summarize_by_cell()`
  `assigned_points_sf` and `select_features_forward()` `train_sf`, with
  the rows as passed.  Layers s2 accepts get the centroids they got, and
  with s2 off the centroids are GEOS's, as before, but for those of the
  next entry.

* **`coerce_to_points()` returned, without a word, a centroid far outside a
  polygon whose edges cross, in a projected CRS.**  GEOS, which takes the
  centroids in a projected CRS, without a CRS, and in degrees with s2 off,
  weights each part of a ring by its signed area.  Where the edges of a
  ring cross (a bowtie), or a hole lies outside its shell, parts of the
  feature count with opposite signs: a bowtie 10 km across, one corner 1 m
  off symmetric, got a point over 33,000 km away, and the functions that
  reduce a layer to points ran on it.  The transform that the refusal of
  the entry above names for such a bowtie in degrees led to that point.
  This held under `"centroid"`, for the rows other than LINESTRINGs under
  `"line_midpoint"`, and for collections under `"auto"`.  A centroid
  outside the bounding box of a feature GEOS reports invalid now stops the
  call with the error of class `spatialkit_features_degenerate`, which
  gives the rows (`rows`), quotes GEOS's reason for each (`reasons`),
  carries `site = "geos_centroid"`, and says to take a point on the surface
  or to repair the features with `sf::st_make_valid()`.  The functions of
  the entry above, `area_of_applicability()` (`train_sf` or `newdata`) and
  `predict_surface()` (`grid`) raise it in their own name, naming their
  argument and its rows.  A centroid inside the box is GEOS's, as before,
  and so is every centroid of a valid feature.

* **The projection `ensure_projected()` chose for a layer in longitude and
  latitude moved, in silence, when s2 refused one of its polygons even
  after the repair.**  The centre that places the projection is the
  centroid of the layer's union on the sphere.  When s2 refused the layer
  and its repair (a bowtie, which s2's repair leaves as it is), the centre
  fell back to the mean of every vertex of the layer, and nothing said so:
  on 100 one-degree squares, one of them drawn with 4,000 vertices, one
  bowtie moved the centre 4 degrees to that square, and the distance and
  area CRSs became an Albers CRS centred there, where the layer without the
  bowtie gets UTM zone 14 and a Lambert azimuthal CRS.  The mean counted
  each ring's closing vertex twice, so the same bowtie listed from each of
  its four corners gave four centres and four CRSs.  The features s2
  refuses are now repaired one by one (without their repeated vertices,
  then with crossing edges split), and the centre is taken over the
  features s2 then accepts; a log line names any feature s2 still refuses,
  which is left out.  The vertex mean is used only when s2 accepts no
  feature, counts each closing vertex once, and is logged.  Every function
  that projects a lon/lat layer through `ensure_projected()` reads this
  centre.  Layers s2 accepts, or whose repair it accepts, get the centre
  and the CRS they got.

* **The distance score of a candidate projection put every feature that
  crosses 180 degrees on the far side of the globe when s2 refused one of
  them.**  Such a feature is scored at its centroid on the sphere, because
  GEOS reads it the long way round.  The centroids of these features were
  taken in one call, and when s2 refused one of them (a vertex listed twice
  in a row) every one kept GEOS's point: a box around Fiji with one vertex
  listed twice scored 25% in the Lambert azimuthal CRS chosen for it, where
  the same box without the repeat scores 0.02%, and `crs_choice` and the
  log reported that figure.  They are now taken feature by feature: a
  feature s2 refuses for a repeated vertex gets the centroid of the same
  ring, one whose edges cross that of its rebuild with the crossing edges
  split, and only a feature s2 refuses after both keeps GEOS's point.
  Features s2 accepts score as before.

* **Folds saved by spatialkit 2.0.0 were applied to other data in silence
  when that data held a lon/lat polygon s2 refuses.**  The `cv_*()`
  functions, `area_of_applicability()` and `fold_separation()` check that
  `folds` were built on the data they are given by retaking the location
  of the rows the folds recorded.  A record saved by 2.0.0 (without
  `kind`) is retaken with the centroid sf takes in the session, and with
  s2 on a polygon s2 refuses stopped it; the check then passed without a
  word, so folds built on other data of the same size were applied.  Such
  a record on such data can only have been taken with s2 off, and it is now
  retaken that way, so folds from other data are refused as they are with
  s2 off.  When no location can be taken for the rows at all, the check is
  skipped with a log line that gives the cause, as it is when one side has
  a CRS and the other has none.

* **A grid laid in degrees (`cellsize` with a geographic `crs`) returned
  cells s2 refuses, and `build_tessellation(method = "square")` or `"hex"`
  stopped on them with s2's bare message.**  With s2 on (the default) sf
  clips such a grid with s2, which writes some clipped cells with a corner
  listed twice in a row, or with two corners that differ by rounding (7e-15
  degrees), which s2 reads as one point: on a 0.5 degree square at
  `cellsize = 0.1`, one of 25 square cells and three of 39 hexagons, and at
  `cellsize = 0.05` seven of 137 hexagons; `build_tessellation()` stopped
  with "Loop 0 is not valid: Edge 5 is degenerate (duplicate vertex)".
  `sf::st_area()` of the grid stopped the same way, and
  `assign_features_to_polygons()` gave a cause that speaks of seeds and of
  cells built in a projected CRS.  `create_grid_polygons()` now returns
  such a cell without the repeat, and without a corner that lies within
  1e-12 degrees of the one kept before it, which moves no corner and changes
  its area by rounding, and `build_tessellation()` runs on the grid; cells
  s2 accepts are returned as they were.  A cell s2
  refuses even then, or a cell of the grid before the clip that s2 refuses
  (a grid laid past latitude 90, on which the clip stopped with s2's
  message), stops the call with an error of class
  `"spatialkit_cells_degenerate"` that names the cells and s2's reason and
  says to pass a projected `crs` or to size the grid by `target_cells` or
  `n` (`approx_n_cells` in `build_tessellation()`).  The index of the points
  in a grid laid in degrees gets the same error where s2 refuses a cell.

* **`voronoi_seeds_random()`, `create_grid_polygons_cached()` and
  `make_folds(method = "block_kfold", boundary = )` unioned the rows of
  `boundary` before repairing them.**  A boundary that `get_voronoi_seeds()`,
  `create_grid_polygons()` and the tessellation builders repair and use
  stopped these three with GEOS's or s2's bare message: in degrees, a
  square with a vertex listed twice or a figure eight stopped
  `voronoi_seeds_random()` ("Loop 0 is not valid: Edge 1 is degenerate
  (duplicate vertex)"), with s2 on or off; in a projected CRS, a bowtie row
  beside a square stopped all three with "TopologyException: side location
  conflict".  The rows GEOS (in a projected CRS) or s2 (in degrees) refuses
  are now repaired before the union; valid rows are unioned as they
  arrived, so a boundary of valid rows gives the seeds, the cache key and
  the folds it gave.  The key of the cache no longer stops on a boundary it
  cannot union: it then hashes the rows, and `create_grid_polygons()` gives
  the outcome it gives uncached.

* **A lon/lat boundary with crossing edges (a bowtie) stopped
  `get_voronoi_seeds(method = "random")` and `"kmeans"`, and
  `create_grid_polygons(cellsize = , crs = <geographic>)`, with s2's bare
  "Loop 0 is not valid: Edge 0 crosses edge 2".**  In degrees, with s2 on
  (the default), `sf::st_make_valid()` is s2's rebuild, which leaves
  crossing edges as they are, where GEOS splits them in a projected CRS, so
  the same boundary ran
  through `create_grid_polygons(target_cells = )` and
  `create_voronoi_polygons()`.  The rows s2 still refuses after the repair
  are now rebuilt with crossing edges split, so the seeds and the grid in
  degrees are placed in the two lobes, as the projected builders place
  them; `voronoi_seeds_random()` repairs them the same way.  A row s2
  refuses even then stops the call with an error of class
  `"spatialkit_cells_degenerate"` that names `boundary`, its rows and s2's
  reason, and says to transform `boundary` to a projected CRS (for the
  grid, to pass a projected `crs` or to size the grid by `target_cells` or
  `n`).

* **`voronoi_seeds_random()` and `get_voronoi_seeds(method = "random")` and
  `"kmeans"` stopped with sf's bare "sampling box is N times larger than
  sampling region" on a lon/lat boundary of small pieces far apart.**  In
  longitude and latitude sf does not sample a region whose bounding box is
  over 10,000 times its area, so two small squares on two continents were
  refused, while the same boundary in a projected CRS was sampled.  Such a
  boundary is now sampled piece by piece: the seeds (or the k-means cloud)
  are shared among the pieces in proportion to their areas in one draw from
  the random stream, so they stay uniform over the whole, and each piece is
  sampled as before.  A piece whose own bounding box is over 10,000 times
  its area stops the call with an error that names `boundary` and says to
  transform it to a projected CRS.  A boundary sf samples whole gives the
  seeds it gave.

* **`voronoi_seeds_random()` and `get_voronoi_seeds(method = "random")` and
  `"kmeans"` exhausted memory on a boundary whose bounding box is far larger
  than its area.**  sf draws the points asked for times the ratio of the box
  to the area in the box and keeps those inside, and the package passed that
  draw on without a bound.  The k-means cloud of two small squares on two
  continents in a projected CRS, or of a strip 50 km long and 1 m wide at a
  slant, stopped with R's "vector memory limit ... reached", and a few
  random seeds there took most of a minute.  A boundary that would take more
  than 10^6 points drawn in its box (or 10 times the points asked for, if
  more) is now sampled piece by piece in any CRS, the points shared among
  its polygons by their areas (or whole, in draws of at most 10^6 points in
  its box, where that takes less time), and a polygon that would still take
  more is drawn in such draws.  Where the points drawn in all would exceed
  2e7 (or 10 times the points asked for, if more), the call stops before
  any draw with an error that names `boundary`, the ratio of its box to its
  area and the points it would take, and says to pass `sample_points` or
  ask for fewer seeds.  The cloud of the 1 m strip (5e7 points in its box)
  is refused in this way; the squares on two continents are sampled.  A
  boundary within the first of these limits gives the seeds it gave.

* **A lon/lat `boundary` with a vertex at a latitude beyond 90 was read as a
  shape folded over the pole.**  `voronoi_seeds_random()` and
  `get_voronoi_seeds()` placed seeds 180 degrees of longitude away, a grid
  laid in degrees put its cells there, `clip_target_for()` returned a target
  a few metres high, `create_voronoi_polygons()` returned 0 cells or
  stopped with "`boundary` must be polygonal" on a polygon, and
  `make_folds(method = "block_kfold")` said none of the points fell inside
  `boundary`, all without saying why.  Such a vertex is no place on the
  Earth, usually metres stamped with a geographic CRS or x and y swapped.
  These functions, `create_grid_polygons()`, `create_grid_polygons_cached()`
  and `build_tessellation()` now refuse such a boundary, read as given,
  with an error of class `"spatialkit_cells_degenerate"` (`site =
  "boundary"`) that names its rows and their latitudes and says to check
  its CRS.  A boundary without a CRS that `create_grid_polygons(crs =)` or
  `build_tessellation(crs =)` stamps with a geographic CRS is tested in the
  same way.  A boundary that reaches latitude 90 and not beyond is read as
  before.

* **A `boundary` with no polygon that has a location was read in silence or
  stopped with a bare or false message.**  A row that is empty, or has a
  vertex with a coordinate that is not finite or lies beyond 1e150, covers
  no ground, but `clip_target_for()` returned a target of area 0 (or 5e204
  with a vertex at 1e200), `create_voronoi_polygons()` returned 0 cells
  with a warning that blamed the points, the random and k-means seeders
  stopped with R's "no applicable method for 'st_geometry'" or put every
  seed near y = 1e199, a vertex at `Inf` moved the cells of the grids, and
  `build_tessellation(method = "hex")` and `"square"` refused a square
  beside an empty row as "not polygonal".  The readers of `boundary`
  (`clip_target_for()`, `create_voronoi_polygons()`, `build_tessellation()`,
  `create_grid_polygons()`, `create_grid_polygons_cached()`, the seeders and
  `make_folds(method = "block_kfold")`) now leave such rows out, read as
  given and again after a move into another CRS, which empties a row it
  cannot place, so a square beside such a row gives what the square alone
  gives.  A boundary with no row left stops the call with an error that
  says `boundary` has no polygon with a location, and why.

* **The random and k-means seeders stopped with R's "no applicable method
  for 'st_geometry'" on a lon/lat boundary whose repair leaves no polygon,
  and never returned on one whose bounding box has no height.**  s2 refuses
  a ring with no height (every vertex on one parallel), and its repair
  leaves nothing; such a row is now left out like a row without a location,
  and a boundary with no polygon left stops the call with the error above,
  whose reason says the repair left no polygon with an area.  A ring of
  three vertices on one parallel is valid under s2, which reads its edges
  as great circles, but its bounding box has no height, and sf drew in that
  box without end: `voronoi_seeds_random()` and `get_voronoi_seeds()` now
  refuse a boundary whose bounding box has no height or no width at once,
  as `create_grid_polygons()` does.

* **The seeders' refusal of a lon/lat boundary whose bounding box has no
  height gave no way forward, and on a boundary sampled piece by piece it
  said the box of `boundary` had no height when only one of its polygons
  had none.**  In a projected CRS such as the one `ensure_projected()`
  picks, a parallel other than the equator is a curve, so a ring with every
  vertex on it has an area there.  The error of `voronoi_seeds_random()`
  and `get_voronoi_seeds()` on such a boundary now says to transform
  `boundary` to a projected CRS in which a parallel is a curve, such as the
  one `ensure_projected()` picks for it (a polar one, such as EPSG:3413 or
  3031, for a ring around a pole), and where only one of the polygons of
  `boundary` has a box with no height, it says so.  In Web Mercator, World
  Mercator or plate carree a parallel stays straight and the ring has no
  area.  A ring on the equator or on a meridian holds no ground, and a
  projected boundary with no height has no area; their errors offer no
  remedy, as before.

* **The seeders still never returned on most lon/lat rings on a parallel,
  and on slivers along one.**  The s2 union of `boundary` moves the
  latitudes of such a ring by a few units in the last place, so its box
  has a height of about 1e-14 degrees and passed the test above, and a
  sliver 1e-7 degrees tall did the same; which latitudes passed differed
  between platforms.  s2 reads the edges as great circles, which leave the
  bounding box, and sf (without the lwgeom package) draws in the box and
  keeps the points in the region, so it drew without end.  A lon/lat
  boundary whose box has less than half its area on the sphere is now
  refused at once, with the remedy above: in such a box sf drew in less
  than half the region, so a triangle 120 degrees wide at 80 N, whose box
  has 0.42 times its area, took seeds only in the part inside the box.
  With `sf_use_s2(FALSE)`, `get_voronoi_seeds()` repaired every lon/lat
  row in the plane, which turned a ring on a parallel into a line and
  refused it as a row s2 refuses; it now repairs them on the sphere, as
  with s2 on.  `voronoi_seeds_random()` on a ring along one meridian, whose
  union on the sphere is a line, returned MULTIPOINT seeds, some empty; it
  now refuses it, as `get_voronoi_seeds()` does, with the reason that the
  union of its rows leaves no polygon with an area, and the reason
  `get_voronoi_seeds()` gives no longer says that s2 refuses the row.

* **A `boundary` row whose repair leaves no polygon made the builders
  refuse a valid boundary or return a line.**  A ring with no height (every
  vertex on one line, or in degrees on one parallel) comes out of the
  repair of `boundary` as a line or as nothing.  Beside a valid square,
  `build_tessellation(method = "hex")` and `"square"` refused the boundary
  as "not polygonal", `clip_target_for()` returned a target with a line in
  it, and the grids and Voronoi cells changed; alone, `clip_target_for()`
  returned a line and `create_voronoi_polygons()` returned cells that held
  none of the points, with a warning that blamed them.  Every reader of
  `boundary` now leaves such a row out before any move, and again after its
  own repair, so a square beside it gives what the square alone gives, and
  a boundary with no polygon left stops the call with the error above.

* **`make_folds(method = "block_kfold")` stopped with R's "no applicable
  method for 'st_geometry'" on a `boundary` that is not an sf or sfc
  layer** (a bounding box, a data frame).  It now refuses such a boundary
  by name, as `build_tessellation()` does.

* **`create_grid_polygons_cached()` handed a grid returned in degrees to a
  session with the other `sf::sf_use_s2()` setting.**  A grid with a
  geographic `crs` and a `cellsize` is laid in degrees and clipped by s2
  with s2 on and by GEOS with s2 off, and `ensure_stable_poly_id()`
  rebuilds the cells of a grid returned in degrees that s2 or GEOS calls
  invalid with s2's rebuild with s2 on and repairs every cell with GEOS's
  with s2 off, so the two settings give two grids.
  The cache key did not hold the setting, so a call got the grid of
  whichever setting filled the cache first.  The key of a grid returned in
  degrees now holds the setting, and each setting gets its own entry; every
  other key is unchanged, so a grid laid and returned in a projected CRS is
  still found in a cache filled before this change.

* **`create_grid_polygons_cached()` returned a different grid from
  `create_grid_polygons()` for a boundary in degrees with `crs` given.**  It
  projected the boundary to the CRS `create_grid_polygons()` picks without
  `crs`, and the builder then transformed that copy to `crs`, where the
  uncached builder transforms the boundary once.  The cell counts and IDs
  agreed, but every corner of the grid differed in the last bits, so no
  cached cell matched an uncached one on geometry.  With `crs` given the
  builder is now handed the boundary as the caller gave it, and the cached
  grid is the one `ensure_stable_poly_id(create_grid_polygons())` returns.
  The same holds for a boundary with no CRS and `crs` given.  The cache
  keys are unchanged.

* **`ensure_stable_poly_id()` moved the corners of valid polygons in
  degrees with s2 on.**  With `make_valid = TRUE` (the default), on a layer
  in longitude and latitude with `sf::sf_use_s2()` on, every row went
  through s2's rebuild, valid or not: valid rows came back with their
  corners moved by up to 7e-8 degrees and a POLYGON as a MULTIPOLYGON,
  where with s2 off or in a projected CRS no corner moved.  The rows
  reached the user through `create_grid_polygons_cached()` with a
  geographic `crs`.  Only the rows s2 or GEOS calls invalid are rebuilt
  now, and the other rows are returned as they arrived.  The IDs are
  unchanged.

* **`ensure_stable_poly_id()` returned rows with crossing edges unrepaired
  in degrees with s2 on.**  With `make_valid = TRUE` (the default), s2's
  rebuild does not split crossing edges, so a bowtie, a figure eight or a
  MULTIPOLYGON whose parts overlap came back invalid under s2 and GEOS,
  where with s2 off GEOS returns it valid.  `sf::st_area()` then stopped
  on the row, and its copy in a projected CRS stopped GEOS's overlay.  A
  row s2's rebuild leaves invalid is now repaired by GEOS on its
  coordinates, as with s2 off.  Such a row is keyed on that repair, so its
  ID is the one s2 off gives.  With s2 on it was keyed from a second repair
  that drops the overlap of overlapping parts, and in a layer that holds
  one, the ID of that row and of its neighbours in the sort order can
  differ from before.  A row GEOS's repair does not make valid under s2
  either, such as a polygon whose edges cross only as great circles, is
  returned as GEOS repaired it, with a warning that names its rows.

* **`ensure_stable_poly_id()` said the geometry returned was the caller's
  own when it was the repair.**  When s2 refuses a feature of the lon/lat
  copy the sort key is measured on, the warning ended "the geometry
  returned is the caller's own".  With `make_valid = TRUE`, the default,
  the rows returned are the repaired polygons, made in the layer's own CRS,
  so a row the repair rewrote came back with another shape than the one
  given.  The warning now says so with `make_valid = TRUE`, and keeps its
  words with `make_valid = FALSE`.  No key, ID or row changes.

* **`ensure_stable_poly_id(method = "surface_point")` passed on sf's
  warning about longitude and latitude.**  The points on the surface are
  taken on the copy in EPSG:4326 the sort key is measured on, so every
  call warned "st_point_on_surface may not give correct results for
  longitude/latitude data" in sf's words, on any polygon layer moved to
  that CRS.  The key is a sort order in degrees, as
  `?ensure_stable_poly_id` says, and the cell summaries already muffled the
  same text.  The warning is no longer
  passed on.  The points, the IDs and the rows are unchanged.

* **`ensure_stable_poly_id()` numbered a layer far from its CRS's origin
  by the row order.**  PROJ returns a row of a projected layer it cannot
  move to the sort CRS (EPSG:4326 by default) as an empty row, or a ring
  without the vertices it cannot place, instead of stopping.  A UTM layer
  shifted by 1e8, 1e12 or 1e15 then had a NaN key for every row, and its
  IDs followed the row order, so the same polygon took another ID when the
  rows came in another order; with `method = "surface_point"` and a ring
  left with two vertices the call stopped with GEOS's "Invalid number of
  points in LinearRing".  `create_grid_polygons_cached()` numbered its
  cells the same way at such offsets.  A move that cannot place every
  polygon with a location now takes the documented fallback a failed move
  takes: the whole layer is sorted in its own CRS, and the log says so.  A
  move that stops no longer passes on GDAL's warning beside that log line.
  A key in a CRS that is not geographic (this fallback, a layer without a
  CRS, `transform_for_sort = NULL` or a projected sort CRS) is now rounded
  to a step that grows with the magnitude of the coordinates: the power of
  two at or above 2^-44 of the layer's largest, never under 1e-7 of its
  unit, with ties ordered at the step of each polygon's own largest
  coordinate.  Cells 1e12 m out were ordered by the last bits of their
  centres, so a hex grid took other IDs there than at the origin, with or
  without a CRS.  A layer within about 1.7e6 units of its origin is still
  rounded to 7 decimals.
  A feature whose point or area the sort copy cannot give is keyed by the
  centre of its box with area 0, as an empty polygon is, where the call
  stopped with GEOS's message.  A ring of fewer than four points encloses
  no area and is dropped from the key and from the repair: GEOS stopped
  the repair on a ring of two ("Invalid number of points in LinearRing"),
  and in degrees sf read such a ring as the whole globe and every call
  stopped with "To set the crs to NA, first remove the full polygons".  A
  polygon left with no ring is keyed by the centre of the box of its points.
  Layers keyed in degrees get the IDs they got; a layer keyed in its own
  unit beyond about 1.7e6 units from its origin gets other IDs only where
  two centres lie within the new step of each other.

* **`ensure_stable_poly_id(make_valid = TRUE)` returned rows its repair
  turned into lines.**  GEOS's repair gives an empty polygon back as
  `LINESTRING EMPTY`, and a ring with no area (every vertex at one
  latitude, or on one line) as a `GEOMETRYCOLLECTION` of lines or an
  empty one, and those rows came back in a layer given as polygons,
  without a word.  An empty polygon is no longer repaired and comes back
  as it arrived.  A row the repair leaves with no polygon, a row of rings
  of fewer than four points among them, comes back as it arrived, with one
  warning of class `spatialkit_cells_degenerate` that gives its rows of
  `polygons_sf`; a row the repair turns into a collection keeps its
  polygonal parts.  Valid layers come back as they
  did.

* **Cells built on seeds placed from lon/lat points could leave observations
  near the edge of a wide layer in no cell.**  The seeding functions attach
  the convex hull of the points they cluster as the seeds' cover, and
  `build_tessellation()` and `create_voronoi_polygons()` move the cover into
  the CRS the cells are built in with its edges read as straight lines in
  degrees, then grow it by 2 percent of its diagonal.  For points in a
  geographic CRS the hull was the one sf builds with s2, whose edges are great
  circles, and away from the equator a great circle bows toward the pole.  Of
  800 points uniform over the box of the conterminous United States, 24 lay
  outside the hull read in degrees; one of them, near the northern edge, lay
  130 km outside it in the cells' CRS against a margin of 127 km, and fell in
  no cell at 44 and at 178 seeds.  `assign_features_to_polygons()` said so
  (`spatialkit_unassigned`), and the summary of the cells built on a
  selection's seeds then raised `spatialkit_cells_differ`.  Round a pole the
  hull, read in degrees, held few of the points: of 600 points north of 70
  degrees, 230 fell in no cell at 10 seeds, where it held none of the seeds
  and was set aside with `spatialkit_no_cover`, and 17 at 60 seeds, where it
  was used.  For points spread round the globe it was the whole sphere, and
  was set aside too.  For points on both sides of 180 degrees, read in
  degrees, it ran the long way round the globe: between 170 E and 170 W and 45
  S and 45 N the cells on 10 seeds were 11 times the area of the points' hull,
  and of 600 points between 105 E and 105 W and 20 N and 70 N, 239 fell in no
  cell at 10 seeds.  The hull is now taken on the coordinates in degrees, as
  it already was with `sf::sf_use_s2(FALSE)` and for points with no CRS, and
  for points on both sides of 180 degrees with their longitudes taken from 0
  to 360 (-170 as 190), the short way round.  Where the cells are built in a
  CRS that splits the map at 180 degrees, as Web Mercator, the CRS for points
  spread round the globe, does, a cover whose longitudes run past 180 is cut
  there into two pieces, one at each side of the map.  On all these layers
  every observation now has a cell, and between 170 E and 170 W the cells are
  1.3 times the points' hull.  `attr(, "cover")` of seeds placed from lon/lat
  points is now that hull, and its longitudes run past 180 for points on both
  sides of 180; s2 reads its edges as great circles, so tested against it with
  s2, points near a long edge can lie outside it.  Cells on such seeds can
  differ at their outer edge, and differ throughout where the cover was set
  aside or ran the long way round.  Layers in a projected CRS, seeds drawn in
  a boundary and the hull a profile keeps are unchanged.

* **A selection by size called the default size the caller's choice.**  On
  a default call, where nothing names a size, `select_resolution()`'s note
  read "the count is set by obs_per_cell = 4.5, the caller's choice; no
  criterion chose it."  It now reads "obs_per_cell = 4.5, the default"
  whenever the size is 4.5, and "the caller's choice" for any other size,
  which only a caller gives, to the selection or to the profile.
  `?select_resolution` and `vignette("resolution")` said that the count of
  a selection by size stands "since the caller set it"; they now say that it
  stands since a size set it and no criterion.

* **The print of `attr(, "pooled")` hid the fit of a column whose `tau2` is
  0.**  `summarize_by_cell()` fits such a column, and every pooled value is
  then its neighbours' mean, but the print showed only "status : no signal
  between cells", also for a response whose cells predict a left-out row
  better than one mean, and listed the column's interval under "Also in the
  frame" as if it belonged to nothing above.  The block is now printed in
  full (value, variances, own weight, the leave-one-out errors and the
  interval), with a line that `tau2` is 0, every own weight is 0 and each
  value is its neighbours' mean.  Whether the cells do better than one mean
  is said on the `loo` line, as for any other column.  The status in the
  record, and the warning, are unchanged.

* **The warning `spatialkit_pool_few_cells` did not say when the pooled
  values of two cells change places.**  With 2 to 9 cells a cell's target
  is the mean of every row outside it, so the overall mean lies between
  the target and the cell's own mean, and a cell whose own weight is below
  its share of the column's rows gets a pooled value on the far side of
  the overall mean.  With two cells and `tau2` at 0 every own weight is 0,
  and each cell's pooled value is the other cell's mean.  The warning said
  only that a cell "with a small own weight" crosses the overall mean and
  that two cells can change order.  It now gives the share of the column's
  rows as the bound, and names each column with two cells and `tau2` at 0,
  whose two values change places.  The "Conditions" entry of
  `?summarize_by_cell` says the same; its "Pooled values" section already
  did.  The values and the condition's fields are unchanged: such a column
  has `cells` 2 and `own_weight` 0.

* **On gridded locations the pooling's neighbours depended on the unit and
  the origin of the coordinates.**  A cell's 8 neighbours are the cells
  whose centroids are nearest, and a tie between two goes to the lower x,
  then the lower y.  Only exact ties counted: centroids of cells on a
  lattice tie exactly in one unit and differ in their last bits in another,
  so the same layer in metres and in feet, or with the origin moved, gave
  some cells other neighbours, and since `tau2` is fitted over the layer
  every pooled value of `summarize_by_cell()` moved, and so did the
  leave-one-out error of `resolution_profile()` on the same cells.  The
  nearest other location of a row alone in its cell, which the error of a
  plain mean reads, had the same fault.  Two squared distances within one
  part in ten million of each other are now a tie, coordinates within that
  part of the distance count as equal, and of cells still tied the one of
  the lower number is taken (in `summarize_by_cell()`, the one whose ID
  sorts first), so the same cells in another unit or with another origin get
  the same neighbours and the same values, cells whose centroids coincide
  among them.  That holds while the spacing of the lattice is above about a
  ten-millionth of its coordinates.  Where three to eight cells' centroids
  coincide the values agree to rounding.  Layers without such ties are
  unchanged.  Another projection still gives other distances.

* **Where more cells' centroids shared a point than a cell has neighbours,
  which of them each took depended on the unit and the origin of the
  coordinates.**  A cell whose centroid shares its point with 8 other cells
  or more takes the first 8 of them by ID, without a search.  The point was
  found by exact equality.  In another unit, or with the origin moved,
  rounding parts such centroids by a few units in the last place, and the
  cells were then searched for among themselves by those distances.  Twelve
  cells of 1 to 12 rows on one point beside a 10 by 10 lattice of one-row
  cells at 250 m gave 3 cells other neighbours in feet and 12 with the
  origin moved by 0.1 m.  Since `tau2` is fitted over the layer, every
  pooled value of `summarize_by_cell()` could move, by up to 0.21 on that
  layer, and so could the leave-one-out error of `resolution_profile()` on
  the same cells.  Centroids that lie within a ten-millionth of the largest
  coordinate of one another, in x and in y, now count as one point; that is
  the reach of the tie rule above.  Where other centroids that close chain on
  further than that, the chain is taken again at half the reach, and so on,
  until its parts span no more than the reach they were taken at.  Such
  cells take the same neighbours, and get the same values, in another unit
  and with the origin moved by little against the coordinates.  Cells whose
  centroids lie that close without coinciding take each other by ID too, not
  by distance.  The reach follows the largest coordinate, so with the origin
  moved far, to the corner of the data say, such centroids can lie beyond it
  and their cells then take the nearest.  A layer
  whose centroids are either equal or further apart than that reach in x or
  in y gets the same neighbours and values as before.

* **On gridded locations the density partition depended on the unit and the
  origin of the coordinates.**  The regions `resolution_profile()`,
  `get_voronoi_seeds()` and `voronoi_seeds_kmeans()` draw under `partition =
  "density"`, and with them `rss`, `loo` and the count chosen, differed
  between metres and feet, or with the origin moved by 0.1 m: on a 40 by 30
  lattice at 250 m the default profile gave `loo` 0.2978 in EPSG:26915,
  0.3068 in UTM 15N in US feet and 0.3006 with the false origin moved.  On a
  lattice many distances tie exactly, and the last bits of the coordinates,
  which move with the unit and the origin, settled them in four places: the
  neighbour searches of the density weights and the nearest-centre search of
  Lloyd's iteration counted exact ties alone; the first k-means++ centre was
  drawn from weights equal but for their last bits, which `sample.int()`
  reads as different; the locations, which the random draws index, were put
  in order by first coordinates equal but for their last bits; and of two
  starts whose sums agree to rounding, the rounding chose.  Now squared
  distances within one part in ten million of each other are a tie in those
  searches, as in the pooling's: a tie between neighbours goes to the lower
  first coordinate, then the lower second, and a tie between centres to the
  lower-numbered one.  The first centre is drawn on the weights with those
  within one part in ten million of each other made equal, first coordinates
  within that part of the extent count as equal when the rows and the
  locations are put in order (so a profile that draws a subsample of the
  rows draws the same locations in every frame), and a later start is kept
  only when its sum is lower by more than that part.  The regions are the
  same in another unit or with another origin as long as no coordinate moves
  by more than about a two-hundred-millionth of the distance between
  neighbouring locations; coordinates rounded more coarsely can still give
  other regions.
  The lattice above gives `loo` 0.3066 and the same regions in all three.
  The scoring of a profile's levels, the refit of `select_resolution()`, the
  edge rule of `assign_features_to_polygons()` and `cell_values_at()` take
  the same tie between centres, so a row the profile scores in the cell of
  the lower-numbered of two tied centres is assigned there too: for such a
  row the assignment's margin also takes in its distance from the boundary
  of the two cells, which midway between centres far apart can exceed one
  part in 1e10 of the coordinates.  Layers on which no two distances, no two
  weights, no two first coordinates and no two starts' sums come within one
  part in ten million of each other are partitioned as before; on others the
  regions, and every figure drawn from them, can differ from before for the
  same seed.  `resolution_profile()` puts the rows in this order under
  `partition = "kmeans"` too, so there a layer with two first coordinates
  that close can also get other regions for the same seed.

* **On gridded locations the k-means partition depended on the unit and the
  origin of the coordinates.**  The regions `resolution_profile()`,
  `get_voronoi_seeds()` and `voronoi_seeds_kmeans()` draw under `partition =
  "kmeans"`, the refit of such a profile by `select_resolution()` and the
  sweep of `determine_optimal_levels()` differed between metres and feet, or
  with the origin moved by 0.1 m: on a 40 by 30 lattice at 250 m the
  default profile under `partition = "kmeans"` gave `loo` 0.2963 in
  EPSG:26915, 0.2998 in UTM 15N in US feet and 0.3000 with the false origin
  moved, and `voronoi_seeds_kmeans(k = 40, partition = "kmeans")` put 391
  and 363 of the 1,200 points in other cells.  `stats::kmeans()` (Hartigan
  and Wong) settles a tie between two centres by the last bits of the
  coordinates, which move with the unit and the origin, and it takes no
  tolerance.  It is now handed the coordinates centred on the middle of
  their bounding box, divided by its diagonal and rounded to 2^-30 of it
  (about a billionth).  The starts are drawn from those coordinates and the
  best start is chosen by its sum on them.  Distinct locations that would
  round to one point keep their coordinates unrounded, so distinct
  locations stay distinct and `stats::kmeans()` refuses no count the
  callers allow.  The cells come back in the order of the rows given, and the
  centres, the within-cell sums of squares and every figure read off the
  cells are taken in the layer's own coordinates.  The seeding functions with `init = "kmeans++"` and
  `determine_optimal_levels()` put the rows in order as `resolution_profile()`
  does, first coordinates equal but for their last bits counting as one
  column.  The lattice above gives the same regions in all three.
  `stats::kmeans()` also settles ties by the last bits of its own compiled
  arithmetic, which differ between platforms, so on a layer where ties
  decide cells the regions, and every figure read off them, can differ
  between platforms for the same coordinates and seed: this lattice gets
  other regions on Linux on x86-64 than on macOS on arm64.  A change of
  unit or origin can still move a coordinate across the middle of a
  rounding step, and does so more often the smaller the extent:
  over 60 lattices, each in five other units or origins, the regions were
  the same on every lattice that spans a kilometre or more and differed in
  205 of 1,560 partitions of lattices that span less ("The partition" in
  `?resolution_profile`).  On layers where no tie decides a cell the regions
  are those of earlier versions, and the centres and sums of squares differ
  from before in their last bits; where two starts reach the same regions
  with sums equal to rounding, the start kept can be the other one, so the
  centres, the seeds and the cells numbered by them can come in another
  order.  On gridded locations the regions, and every figure drawn from
  them, can differ from before for the same seed.

* **The count of distinct locations that caps the ladders and the seeds
  depended on the unit of the coordinates.**  `resolution_profile()` (its
  support ceiling, under both kinds), `determine_optimal_levels()` (the top
  of its ladder) and `get_voronoi_seeds()` and `voronoi_seeds_kmeans()` (the
  most seeds they return) counted the distinct locations after rounding the
  coordinates to 8 or 10 decimals of the layer's unit.  Forty locations,
  each listed twice 2e-8 m apart, counted 80 in metres and 40 in UTM 15N in
  kilometres: the profile's ladder stopped at 40 cells in kilometres and at
  79 in metres, `voronoi_seeds_kmeans(k = 79)` returned 58 seeds in
  kilometres and 79 in metres, and `determine_optimal_levels()` returned
  other counts.  The locations are now counted on the coordinates
  `stats::kmeans()` is handed (centred, divided by the diagonal of their
  bounding box and rounded to 2^-30 of it, distinct locations kept
  distinct), so the count is the same in every unit: 80 in metres, in
  kilometres and in US feet.  Locations closer than 1e-8 of the layer's
  unit (1e-10 for the seeds) could count as one and now count as distinct,
  so such a layer can get a higher ceiling and more seeds than before.

* **`predict_surface()` took a grid point's covariates from whichever of
  several equally near features of `covariates` GEOS met first, so the
  surface moved with the unit and the origin of the coordinates.**  A grid
  laid between the points of a lattice is the same distance from four of
  them, and `sf::st_nearest_feature()` settled the tie by the last bits of
  the distances.  On a 20 by 15 lattice at 250 m, with the grid at the
  midpoints, 147 of 266 grid points took another covariate row in UTM 15N
  in US feet than in EPSG:26915, and 189 with the false origin moved by
  0.1 m; `.pred` moved by up to 0.03.  For point covariates, squared
  distances within one part in ten million of the least, plus a margin, are
  now a tie, and
  the tie goes to the point lowest in x (x values within one part in ten
  million of the diagonal of the points' bounding box, or within one part
  in 1e10 of their largest coordinate where that is more, counting as
  one), then lowest in y (y values read the same way), then first in row
  order.  The grid above takes the same covariates in all three frames,
  with or without the FNN package.  Two points that differ only in the
  last bits of their coordinates (one station recorded twice, one
  billionth of a metre apart) now take the tie the same way with another
  origin, where a grid point the same distance from both, 570 m away, took
  one row in metres and US feet and the other with the false origin moved
  by 0.1 m; the first in row order takes it now, and so does a grid point
  near the pair, however near (see the entry on a point a hair from two
  records of one station).  A grid point that ties
  with no other point takes the row it took before.  Lines and polygons
  had the same fault: a grid point on the edge two covariate
  polygons share, or at the shared end of two lines, is at distance 0 from
  both, and GEOS gave it either.  On 10 by 8 square cells of 500 m with
  the grid on their shared edges, 34 of 72 grid points took another cell
  with the false origin moved by 0.1 m, and 29 with the rows of
  `covariates` shuffled.  Distances within one part in ten million of the
  least, plus one part in 1e10 of the largest coordinate of `covariates`,
  are now a tie, and the tie goes to the feature whose bounding-box centre
  is lowest in x, then lowest in y (x values and y values within one part
  in ten million of the diagonal of the bounding box of `covariates`
  counting as one), then first in row order: each of the 72 takes the cell
  to its left in every frame and row order.  Features whose box centres
  all coincide (a polygon with a hole and the island that fills it, alone
  in the layer, the polygon first) now take the tie the same way with
  another origin, where a grid point on the hole's ring took the polygon in
  metres and US feet and the island with the false origin moved by 0.1 m.
  A grid point with one nearest feature takes the row it took before.
  Layers in degrees are searched as before.
  `?predict_surface` (`covariates`) states the rule.  A point of a grid
  kept as given with a coordinate that is not finite or lies beyond 1e150
  now gets `NA` covariates, whatever the layer, where it took those of a
  row GEOS gave it.  With point covariates its `.pred` was `NA` as well;
  with lines or polygons, or with a layer in degrees, it got a finite
  `.pred` from that row where `?predict_surface` (`grid`) states `NA`, and
  now gets `NA`.

* **A point a hair from two records of one station took the nearer record
  by the last bits of the coordinates, which move with the unit and the
  origin.**  The nearest-centre search of `predict_surface()` (point
  covariates), `cell_values_at()` (the nearest scored location, seed or
  centre), `resolution_profile()` (the cell of each scored row and the
  steps of the density partition), `make_folds()` (`"knndm"`),
  `assign_features_to_polygons()` (a point on the edge between two cells
  a profile scored) and `get_voronoi_seeds()` (the seeds of a selection)
  counted two centres as a tie when their squared distances agreed to
  within one part in ten million.  Two records one billionth of a metre
  apart have squared distances from a point r away that differ by about
  2r times that billionth, within one part in ten million of r^2 only
  while r is more than about 2 cm.  Closer, the nearer by the last bits
  was taken: of 320 grid points 0.1 mm to 2 cm from such a pair, 48 to 60
  took another covariate row in UTM 15N in US feet or with the false
  origin moved by 0.1 m than in EPSG:26915, and `cell_values_at()` filled
  20 of 160 locations 1 mm from the station from another cell.  A squared
  distance now ties with the least when it is no more than the square of
  the least distance lengthened by a margin, plus one part in ten million
  of the least squared distance.  The margin is one part in 1e10 of the
  largest coordinate of the centres, at most a millionth of the diagonal
  of their bounding box and at least four times the spacing of doubles at
  that coordinate, as for the lines and polygons of `covariates`.  No grid
  point and no location of those layouts moves between the frames.  A
  point whose two nearest centres differ in distance by less than the
  margin, and that took the nearer, now takes the lower-numbered one (for
  `predict_surface()` and `cell_values_at()`, the one their column order
  puts first).  The margin grows with the largest coordinate, as that of
  lines and polygons does: about 0.46 mm at a UTM northing near 4,550 km,
  about half a micrometre once the origin lies next to data 5 km across.
  A point whose distances to two centres differ by less than the margin
  with one origin and by more with another takes the lower-numbered centre
  with the first and the nearer with the second, in the searches above
  that work in the layer's own coordinates (the k-means steps of the
  density partition and of `"knndm"` run in a frame centred on the data and
  do not).  `?predict_surface`, `?cell_values_at`,
  `?resolution_profile`, `?make_folds`, `?summarize_by_cell` and
  `?assign_features_to_polygons` state the rule.

* **The `attr(, "kmeans")` record of `get_voronoi_seeds()` for the scored
  centres of a profile or selection gave a point the same distance from two
  seeds to whichever `sf::st_nearest_feature()` returned, so `cluster`,
  `size` and `withinss` moved with the unit and the origin, and could
  differ from the cells the profile scored.**  On a 40 by 30 lattice at
  250 m with the 267 centres of `resolution_profile(seed = 1)`, `cluster`
  differed from the one in EPSG:26915 in 5 rows in UTM 15N in US feet and
  in 2 with the false origin moved by 0.1 m, and in each frame 1 to 4 rows
  were not in the cell the profile scored them in.  Each point now goes to
  its nearest seed by the profile's rule: squared distances within one part
  in ten million of the least are a tie, and a tie goes to the lower
  `seed_id`.  The record then holds the profile's cells in every frame, and
  `tot_withinss` scales with the square of the unit.  `?get_voronoi_seeds`
  states the rule.

* **`build_tessellation()` let the unit and the origin of the coordinates
  choose the cell of a lattice point on a cell edge, and of a point snapped
  to two cells the same distance away.**  A point on the edge two cells
  share went to the smaller `cell_id` only while the edge passed through
  the point exactly; the last bits of the coordinates put it a hair inside
  one cell or the other.  On a 40 x 30 lattice at 250 m with square cells
  of 1000 m, 80 of 1200 rows were in another cell in UTM 15N in US feet than
  in EPSG:26915 and 226 with the false origin moved by 0.1 m; with
  hexagons, 156 with the origin moved.  A point within a margin of one part
  in 1e10 of the largest coordinate of the cells from a cell's edge is now
  in that cell, and a point in two or more cells takes the smallest
  `cell_id`.  A point just outside every cell, snapped to the nearest cell,
  went to whichever cell `sf::st_nearest_feature()` returned when two were
  the same distance away: of 180 points beside the corners of boundary
  cells, 4 to 27 went to another cell in US feet, international feet,
  kilometres and with moved origins.  Among the cells as near to one part
  in ten million it now goes to the smaller `cell_id`.  The index is the
  same in every frame on both layers.  `params$snapped` no longer lists a
  point within the margin of an edge, which is in that cell.  The same index
  serves `create_voronoi_polygons()` and the check of a `cover`.
  `?build_tessellation` (`index`) states the rule.

* **`build_tessellation(method = "triangles")` split the squares of a
  lattice by the last bits of the coordinates.**  The four corners of each
  square lie on one circle, and qhull split the square by the rounding of
  the coordinates it was handed, which moves with the unit and the origin.
  A 37 x 25 lattice at 250 m gave 1728 triangles in EPSG:26915 and in UTM
  15N in US feet, 860 of them other triangles in US feet, and 1827 with the
  false origin moved by 0.1 m, 99 of them flat.  qhull, and GEOS where the
  geometry package is not installed, are now handed the points centred on
  the middle of their bounding box, divided by its diagonal and rounded to
  2^-30 of it, and each triangle is drawn on the points' own coordinates.
  The lattice gives the same 1728 triangles in every frame, none flat, but
  where a point lies within its own rounding noise of the middle of such a
  step.  A jittered lattice gives the triangles it gave.
  `?build_tessellation` says so.

* **`create_voronoi_polygons()` and `build_tessellation(method = "voronoi"
  | "triangles")` merged near-coincident points by ten decimals of the
  layer's unit, so the number of cells moved with the unit and the
  origin.**  60 points with copies of 20 of them moved by 3e-8 m gave 80
  cells in metres and 65 in kilometres; moved by 1e-9 m, 80, 62, 79 and 60
  cells in metres, in US feet, with the false origin moved by 0.1 m and in
  kilometres.  With `keep_duplicates = FALSE` points are now merged when
  they round to one point in coordinates centred on the middle of their
  bounding box, divided by its diagonal and rounded to 2^-30 of it, so
  points within about a billionth of the extent of each other count as
  coincident in any unit or origin: both sets of copies give 60 cells in
  every frame.  Points that close but more than ten decimals of the unit
  apart, which were kept, can now be merged: copies 1e-6 m apart on that
  5 km layer gave 80 cells and give 64.  `?create_voronoi_polygons` and
  `?build_tessellation` (`keep_duplicates`) said the argument had no effect,
  which was false for points that close; they now say what each value
  merges.

* **`create_grid_polygons(cellsize = )` laid an extra row or column of
  sliver cells over a boundary whose side is a whole number of cells, in
  some units and origins.**  sf counts the cells per side as
  `ceiling(width / cellsize)`, and rounding put a 9000 m side at
  9.0000000000001741 cells of 1000 m with the false origin moved by 0.1 m.
  The square grid had 60 cells there against 54 in EPSG:26915 and in UTM
  15N in US feet, the smallest of 5e-8 m^2, and the hexagonal grid 79
  against 76.  A side that exceeds a whole number of square cells by no
  more than a margin of one part in 1e10 of the largest coordinate of the
  boundary now gets that number, the side `target_cells` derives as well,
  and with `clip = TRUE` a cell that meets the boundary only within that
  margin of its edge is dropped.  The grids above have 54 and 76 cells in
  every frame.  `build_tessellation(method = "square")` on a 37 x 25
  lattice at 250 m in its own box gave 70 cells with the origin moved
  against 54, and `approx_n_cells = 54` gave 63; both give 54 in every
  frame, and so does `create_grid_polygons_cached()`.  A side more than
  that margin over a whole number keeps its row or column of narrow cells,
  however many cells it holds, as before.  `?create_grid_polygons`
  (`cellsize`, `clip`) states the rule.

* **`cell_values_at()` without `seeds` settled a tie between two scored
  locations of one lattice column by their exact x, so the fill from the
  nearest scored location moved with the origin.**  Of two locations of one
  column the same distance from a location asked for, the one whose x was
  lower in its last bits won, and a lattice moved to another false origin
  no longer has one x per column.  On a 20 x 15 lattice at 250 m, 2 of 14
  locations midway between two points of the first column took the upper
  point with the false origin moved by 0.1 m, where EPSG:26915 and US feet
  take the lower, and their `filled_from` and `value` moved with it.  The
  scored locations are now in column order (x values within one part in
  ten million of the diagonal of their bounding box, or within one part in
  1e10 of their largest coordinate where that is more, counting as one),
  then in row order (y values read the same way), then by the first cell in
  ID order among their rows, so the tie goes to the lower point in every
  frame, as "Ties" in `?cell_values_at` said of another unit or origin.
  Two scored locations that differ only in the last bits of their
  coordinates (one station recorded twice, in two cells) now take the tie
  the same way with another origin, from the cell first in ID order as
  rows at one location do, where the fill came from one cell in metres and
  US feet and from the other with the false origin moved by 0.1 m.  That
  section now states the column rule.

* **`summarize_by_cell()` read columns at extreme scales as having no
  spread.**  The test compared the spread inside the cells with the size of
  the values, so 1e15 plus a spread of 1 was called constant, and so were
  values near 1e200 and 1e-170, whose squares overflow and underflow.  The
  warning said their values "do not vary inside any cell", and the pooled
  values were the plain means.  The test now also reads the column centred
  on its median and scaled by its largest distance from it, and a column
  has no spread only where neither reading finds any.  On that column a
  standard deviation, over all the rows or inside the cells, of at most
  three units in the last place of the column's largest absolute value is
  rounding: a column whose values differ only by their rounding (0.1 beside
  0.3 / 3, 2 beside sqrt(2)^2) has no spread by either reading, as before,
  while 1e15 plus a spread of 1, at eight such units, has spread.  Where the
  raw values lose their spread to rounding, the fit is made on the scaled
  column and put back in the column's units: for 1e15 plus a spread of 1 the
  weights, `sigma2`, `tau2` and the leave-one-out errors are those of the
  spread alone, and the pooled values are those of the spread plus 1e15, to
  the rounding of a value stored at 1e15.  The warning for no spread inside
  any cell now says the values do not vary there "to the precision of the
  values".  A column that is not constant, with an absolute value above
  1e150 or a largest distance from its median below 1e-150, gets `NA`
  pooled columns, with a `spatialkit_pool_unavailable` warning whose
  `reason` is `"values out of range"`, and its plain columns are returned:
  beyond these bounds the squares of many values can overflow or lose their
  precision.  A few columns just beyond them were still fitted and are now
  refused: at 800 rows, 1e150 to 1e152 times a spread of 1, and 1e-151 to
  1e-160 times a spread of 1.  Inside the bounds, every column whose raw
  values show their spread gives the numbers it gave, except the standard
  errors of the leave-one-out errors (`loo_se` in the record and in
  `resolution_profile()`), which square the squares: for residuals above
  about 1e77 they were `Inf`, and below about 1e-80 they lost their
  precision and then were 0.  They are now taken on the residuals divided
  by a power of two.  `resolution_profile()` reads its estimate's columns
  the same way: `attr(, "estimate")$constant` and `loo_one_cell` follow the
  engine's test, so `loo` is computed for 1e15 plus a spread of 1, and
  values beyond those bounds leave `loo`, `loo_se` and `own_weight` `NA` at
  every level with the same warning; `print()`, `summary()` and
  `select_resolution(criterion = "loo")` of such a profile say why.  At
  values beyond the bounds the test of a constant response for `cp`
  (`attr(, "one_cell")$constant`) finds the response constant, one
  `spatialkit_pool_unavailable` warning names `cp`, `moran_i` and `moran_z`
  as well, and `summary()` says the values are out of range rather than
  constant.  The test for `cp` now reads the values the same way as the
  estimate's (see the entry on a response at a large offset).

* **`summarize_by_cell(deff = "kish")` called a column in small units
  constant.**  The ICC skipped a column as one that does not vary when its
  standard deviation was below `.Machine$double.eps` (about 2.2e-16),
  whatever the size of its values.  A column that varies, in units that make
  its values of order 1e-16 or smaller, got an ICC of `NA` with the reason
  "the column does not vary", no design effect and the uncorrected standard
  errors, without a warning, while the same column in larger units got its
  correction: on 800 points in 178 cells the ICC was 0.478 at a scale of 1
  and of 1e-14 and `NA` at 1e-17, where the standard errors came out 2.2
  times smaller.  Where only the rows with a cell id were that small, the
  ICC was `NA` with "its values differ only on rows that fall in no cell",
  and its standard errors `NA`.  Where the standard deviation is below
  that floor or not finite, the column is now read again centred on its
  median and scaled, as for the pooled values, and it does not vary only
  where neither reading finds spread; a column with spread gets its ICC
  from the scaled column, which is the ICC of the column.  A column whose
  values are all equal still has an ICC of `NA` and standard errors of 0.
  A column whose values differ only by their rounding (0.1 beside
  0.3 / 3) still has an ICC of `NA` with "the column does not vary"; its
  standard errors are 0 in a cell whose values are equal, and of the size
  of that rounding in a cell whose values differ.  Such a column was read
  that way only while its standard deviation was below the floor: the same
  values multiplied by 16 or more, whose rounding is above it, got an ICC
  of 0.395 from that rounding and a design effect on their standard errors,
  while the pooled summary of the same call said the column does not vary.
  Now the column is also read again where its sum of squares shows no
  spread on the scale of its values, so a column that differs only by its
  rounding does not vary at any size.  A column whose standard deviation
  overflows, with values near 1e200, had an ICC of `NA` with "the variance
  decomposition is undefined"; its ICC is now estimated.  Columns with
  spread whose standard deviation is finite and at or above the floor give
  the ICC they gave.

* **`summarize_by_cell()` gave standard errors of `NA` or 0 in silence for a
  column in very large or very small units.**  The standard deviation inside
  a cell squared the deviations of the raw values from the cell mean.  For
  values above about 1e154 the squares overflowed, and in each cell where
  they did `..sd_*` was `Inf` and the standard errors and intervals were
  `NA`, under any `deff`, with no warning.  For a spread below about 1e-154 the squares lost
  their bits, and at a spread of 1e-170 `..sd_*` and the standard errors
  were 0, which says that the mean has no error.  Where the raw standard
  deviation is not finite or is below 1e-150, it is now taken on the cell's
  values centred on their median and divided by powers of two, so
  `..sd_*`, `..se_*` and the intervals of a column follow its unit.  A cell whose standard
  deviation is beyond the largest double (a spread of order 1e308) still has
  `..sd_*` of `Inf` and standard errors of `NA`, now with a warning of class
  `"spatialkit_sd_overflow"` that names the columns and gives the number of
  such cells of each (fields `columns` and `cells`).  Every other standard
  deviation is the one it was.

* **`estimate_sac_range()` and the test for `cp` in `resolution_profile()`
  called a response at a large offset constant.**  Both read the raw values
  alone.  `estimate_sac_range()` refused a response whose standard
  deviation is under 1e-12 of its largest absolute value, so 1e15 plus a
  spread of 1, eight units in the last place of its values, was "the
  response is constant" and, with `predictor_vars`, "the residuals on
  predictor_vars are constant: the predictors explain the response
  exactly".  The profile left `cp`, `moran_i` and `moran_z` `NA` at every
  level while it scored `loo`, and `select_resolution(criterion = "cp")`
  refused.  Where the raw values show no spread, both now read them again,
  as `summarize_by_cell()` does: a response has spread when, centred on its
  median and scaled, its standard deviation is above three units in the
  last place of its largest absolute value; residuals are constant when
  they are within the rounding of the regression or within three such
  units of the response.  Such a response is centred on its median before
  the predictors are fitted to it and before the test for spatial structure
  and, in `resolution_profile()`, Moran's I of the cell means, none of which
  change under a shift.  A response whose raw values show their spread
  through a trend, such as 1e15 plus 1e5 times a predictor plus a spread of
  1, is fitted on its raw values; where the OLS residuals of that fit show
  no spread, the predictors are fitted again to the response less its
  median, and those residuals are read and scored.  On the raw values the
  rounding of these regressions at 1e15 was larger than the spread:
  `resolution_profile()` with `predictor_vars` scored `cp` on residuals
  with a sum of squares about 150 times the spread's, and gave a `cp` to a
  response the predictors fit exactly.  A response of one
  value, values that differ only by their rounding (0.1 beside 0.3 / 3) and
  the residuals of an exact fit stay constant.  Every response whose raw
  values show their spread, and with `predictor_vars` whose OLS residuals
  on the raw values show theirs, gives the numbers it gave, except the
  residuals of an exact fit that the bound on the raw values took for
  spread (see the next entry).  Values whose
  spread is between three units in the last place and 1e-12 of their size,
  such as 3 times (1 + 4 eps times noise), were constant to
  `estimate_sac_range()` and now reach the variogram fit, and so do
  residuals whose spread is between the rounding of the regression and
  1e-12 of the response's size, such as those of 1e3 times a predictor
  plus 1e-11 times a field, on that predictor.  Where the test for `cp`,
  which reads the scored rows of the subsample, finds the response constant
  and the estimate, which reads every row of the layer (of its selection
  half under `select_on = "split"`) with a value, does not, the warning,
  `print()` and `select_resolution(criterion = "cp")` now say that.  They
  said that the test for `cp` reads the raw values alone and the estimate's
  also reads them centred and scaled, which was untrue of that case.

* **The test for `cp` in `resolution_profile()` scored the residuals of some
  exact fits, where `estimate_sac_range()` refused the same layer.** With
  `predictor_vars`, the test read the OLS residuals against 100
  `.Machine$double.eps^2` times the response's sum of squares first, and
  only where that bound found no spread against the rounding of the
  regression.  `estimate_sac_range()` read the bound on the regression's
  rounding only where its floor, a standard deviation of 1e-12 of the
  response's largest value, found no spread.  `lm.fit()` and `lm()` leave
  more rounding than either bound on some exact fits, and more with more
  rows.  On 800 points, 7 w - 2 w2 on w and w2 left a sum of squares more
  than twice the first bound, so `cp` was scored on rounding with no
  `spatialkit_constant_response` warning and `select_resolution(criterion =
  "cp")` returned 2, while `estimate_sac_range()` refused the layer as an
  exact fit by a margin of 4%.  Of 100 exact fits on one to three predictors
  of that layer, at offsets 0 and 1e6, both scored 8 of the 200: the profile
  gave them a `cp`, and the estimate ran its test for spatial structure on
  the rounding and returned "no spatial structure detected".  On a predictor
  far from 0, such as a northing of 4.5e6 plus the y coordinate, `lm()`
  leaves the residuals of an exact fit at about 2 `.Machine$double.eps`
  times the size of the terms, which is thousands of times the response, so
  the floor found spread in them and the estimate ran its test for spatial
  structure on the rounding.  Both functions now read one bound on the
  residuals alone, whether or not the earlier bounds find spread, taken on
  the size of the fit's terms (each row's absolute intercept plus the
  absolute products of its predictors and their coefficients) and read on
  the residuals taken again from the coefficients, which carry the rounding
  of the response and not that of the fit's arithmetic.  Residuals with a
  root mean square of at most 10 `.Machine$double.eps` times that of the
  terms, or of three units in the last place of the response, are constant,
  whatever the number of predictors and rows.  So the residuals of an exact
  fit are constant in both: `estimate_sac_range()` says the predictors
  explain the response exactly, `cp`, `moran_i` and `moran_z` are `NA` with
  the constant-response warning, and `select_resolution(criterion = "cp")`
  refuses as it does for any constant response.  The same bound refuses
  Moran's I of the cell means in `determine_optimal_levels(criterion =
  "morans_i")` where it regresses them again about the response's centre, so
  the cell means of an exact fit on an integer predictor over 300 cells are
  no longer scored.  Residuals above the bound are scored as before: those
  of 1e3 w + 1e-11 f on w, 45 `.Machine$double.eps` times the size of the
  terms, reach the variogram fit and `cp`.

* **At a large offset, Moran's I of the cell means scored the residuals of
  an exact fit and refused cell means with real spread, so the counts of
  `determine_optimal_levels(criterion = "morans_i")` and the `moran_z` of
  `"morans_i"` and `"combined"` depended on the response's origin.**  Where
  the raw values lost their spread and the response was taken about its
  centre (1e12 or 1e15 plus a spread of a few units), the residuals of the
  cell means were tested only against the centred cell means, and the
  residuals of an exact fit, which are the rounding of the offset, were
  scored: on 800 points, 30 of 30 random fits a + b w + c w2 at 1e15.
  Where the cell means were regressed again about the response's centre,
  the residuals were refused when within three units in the last place of
  the raw response, the rounding of a single row, while the rounding and
  the spread of a mean of n rows are both those of a row over sqrt(n):
  1e15 + 1e5 w + f on w gave `NA` at 20 to 80 cells, and the call fell back
  to the geometric ranking with a log line that said the cell means differ
  by rounding alone.  The rows are now fitted first, and where their
  residuals are within the rounding of a row, as `resolution_profile()`
  reads its own, I and z are `NA`.  Both passes then test the residuals of
  the cell means against the rounding of cell means: three units in the
  last place of the raw response over the square root of each cell's rows,
  plus the rounding of the sums each mean is formed from, summed over the
  cells on the residual degrees of freedom, or the rounding of the
  regression on the size of the terms of the cell design.  The exact fits
  give `NA` at offsets 0, 1e6, 1e12 and 1e15, and so do exact fits on
  cells of 1,200 to 2,200 rows and those on a predictor along which the
  response at 1e15 changes by a unit in the last place over the layer.
  The spread above gives the same I and z at 0 and at 1e15, so the counts
  and `moran_z` are the same at both offsets.  Residuals that the rounding
  of a row holds (1e15 + 1e5 w + 0.3 f) give `NA`, and residuals above the
  bounds are scored as before.

* **A predictor at a large offset was left out of the OLS fits in silence,
  so `rss`, `cp`, Moran's I and z, the nugget and the range depended on the
  predictor's origin.**  The QR decomposition behind `lm.fit()` and `lm()`
  drops a column whose part orthogonal to the earlier columns is under 1e-7
  of its norm, so a predictor whose spread is under about 1e-7 of its size,
  such as a timestamp in seconds or 1e8 plus a spread of 1, was taken as
  collinear with the intercept, and the residuals scored were those of the
  intercept alone.  A little below that size the column was kept, but the
  fit worked on the rounding of the offset, and with a predictor that
  explains most of the response the numbers moved with the origin by up to
  1e-4 of their size.  This held in the OLS fit of `resolution_profile()`,
  in the regression of the cell means behind `moran_i` and `moran_z` there
  and in `determine_optimal_levels(criterion = "morans_i")` and
  `"combined"`, and in the OLS detrending of `estimate_sac_range()`, which
  `resolution_profile(sac = NULL)` reads.  On 600 points with a predictor
  of unit spread, the profile's `rss` at its coarsest level was 548 at
  offset 0 and 2768 at 1e8, and on a clustered layer the estimate's nugget
  was 0.010 and 4.01; the warning of `determine_optimal_levels()` said the
  cell means were regressed on the cell-mean predictors.  A predictor whose
  centre is more than 1000 times its spread is now taken less its centre
  before these fits (for the cell means, before they are formed), and a
  column of times, dates or time differences that `estimate_sac_range()`
  is given is read as the numbers `lm()` reads.  The fit of a model with an
  intercept does not change under a shift of a predictor, so these numbers
  are now those of the predictor at 0, to rounding.  Every other predictor
  is fitted as it was, to the last bit, and a predictor that is constant,
  or a linear combination of the others, is still left out.  The test for
  an exact fit still reads the size of the terms of the predictors as
  given, so the residuals of a response computed from a northing are still
  constant.  The REML detrending of `estimate_sac_range()` has its own
  entry below.

* **`estimate_sac_range(detrend = "reml")` fell back to OLS on a predictor
  at a large offset, with the reason that the REML fit did not converge.**
  `nlme::gls()` fails on a predictor such as 1e8 plus a spread of 1, and the
  estimate warned that "the REML detrending on p did not converge",
  detrended by OLS and returned `detrend_method` `"ols"` and the range and
  nugget of the OLS residuals, where the same predictor at 0 gives the REML
  ones: on 600 points with a predictor of unit spread, a range of 1429 and
  a nugget of 0 at 1e8 and 1e15, against 1562 and 0.0214 at 0.  At 1e5 the
  REML fit ran on the rounding of the offset and moved in the seventh
  digit.  A predictor whose centre is more than 1000 times its spread is
  now taken less its centre before the REML fit, as in the OLS detrending,
  and the test for spatial structure and the test for an exact fit read
  the same design, the latter with the centres added back, as there.  The
  REML fit of a model with an intercept does not change under a shift of a
  predictor, so the range and the nugget are now those of the predictor at
  0, to rounding.  Every other predictor is fitted as it was, to the last
  bit, and a fit that fails still falls back to OLS with the warning.

* **`estimate_sac_range(detrend = "reml")` drew its subsample of a lattice
  by the last bits of the coordinates.**  The REML fit runs on at most
  `reml_max_n` locations and one row a location, drawn by position from
  the rows in the exact order of their first coordinate.  After another
  unit or origin the first coordinates of one column of a lattice differ in
  their last bits, so the draw took other rows: on a 40 by 30 lattice at
  250 m with `reml_max_n = 300` the nugget was 0.206 in EPSG:26915 and in
  UTM 15N in US feet and 0.241 with the false origin moved by 0.1 m.  The
  rows are now put in order as `resolution_profile()` puts them, first
  coordinates within one part in ten million of the extent counting as one
  column, so the draws take the same points in every unit and origin.
  Layers on which no two first coordinates come that close get the order,
  the subsample and the result they got before.

* **`estimate_sac_range(detrend = "reml")` fell back to OLS on a response
  at a large offset, with the reason that the REML fit did not converge.**
  `nlme::gls()` fails on a response at 1e8 plus a spread of 2 as it does on
  such a predictor: on 300 points the estimate warned, detrended by OLS and
  returned a range of 1400 and a nugget of 0, where the same response at 0
  gives the REML range of 1361 and a nugget of 0.0142, and at 1e7 the REML
  range moved in the fifth digit.  A response whose centre is more than
  1000 times its spread is now taken less its centre before the REML fit,
  as a predictor is, and the centre is added back to the intercept that the
  test for an exact fit reads.  Every other response is fitted as it was,
  to the last bit.

* **`estimate_sac_range()` drew its `n_max` subsample of a lattice by the
  last bits of the coordinates.**  The subsample is drawn by position from
  the rows in the order of their coordinates, and after another unit or
  origin the first coordinates of one column of a lattice differ in their
  last bits: on a 40 by 30 lattice at 250 m with `n_max = 300` the range
  was 647 m in EPSG:26915 and in UTM 15N in US feet and 795 m with the
  false origin moved by 0.1 m.  The rows are now put in order by columns
  of the first coordinate, as in the REML fit above, so the help's
  statement that the answer does not depend on the unit or the origin of
  the coordinates holds for the subsample.  Layers on which no two first
  coordinates come that close get the order, the subsample and the result
  they got before.

* **`estimate_sac_range(detrend = "reml")` gave REML in one unit or origin
  and fell back to OLS in another where the REML range is not
  identified.**  On 100 to 150 points of a lattice the REML range can run
  to 1e8 m and more with a nugget proportion near 1e-6.  The likelihood is
  flat there, and whether `nlme::gls()` returned a fit or failed on the
  same rows was set by the last bits of the coordinates it was handed as
  supplied.  On a 20 by 15 lattice at 250 m with `reml_max_n = 100`,
  EPSG:26915 gave `detrend_method` `"reml"` and a nugget of 0.217, and UTM
  15N in US feet and the false origin moved by 0.1 m gave `"ols"`, a nugget
  of 0.154 and the warning that the REML detrending did not converge.  On a
  24 by 18 lattice with `reml_max_n = 120` it was the other way round.  The
  REML fit now reads the coordinates of the rows it uses centred on the
  middle of their bounding box, divided by its diagonal and rounded to
  2^-30 of it (about a billionth), as the partition hands them to
  `stats::kmeans()`, and its range is multiplied back into the layer's
  unit.  The response and each numeric predictor, whose unit and last bits
  set the outcome the same way (the response times 1000, or changed by two
  units in the last place, gave `"reml"` where the response as given gave
  `"ols"`), are read in the
  same frame of their range, and the coefficients, the variance and the
  log-likelihoods are taken back to their units.  The extent that bounds
  the starting ranges no longer has a floor of 1 in the layer's unit, which
  on a layer under 1 km across in km set other starts than in metres.
  `gls()` then gets the same input in every unit and origin, and every
  configuration tried gives one `detrend_method`, one warning and one
  nugget in each.  The range is `NA` on both lattices, before and after.
  On a 40 by 30 lattice whose REML range is identified, and on the same
  lattice jittered, the range changed in the sixth digit or less.

* **`kriging_adequacy()` kriged a response at a large offset on the
  rounding of the offset.**  At 1e15 plus a spread of 1 a double holds the
  values to an eighth, and the replicate means, the kriged and plain cell
  means and the cross-validation folds were computed on the raw values, so
  their rounding was of the order of the difference between a kriged and a
  plain mean.  On one layer of 800 points and 36 cells, the same stored
  values gave shifts beyond 1 in 10 cells at 0 and in 14 at 1e15, and
  `kr_shift` moved by up to 1.  Where the raw values lose their spread, the
  response is now centred as `estimate_sac_range()` centres it before any
  of these are computed, and the centre is added back to `mean` and
  `kr_pred` alone.  `kr_shift`, the counts `print()` gives, the replicates'
  within-location variance and the cross-validation statistic are those of
  the same values at 0; `mean` and `kr_pred` are in the response's units,
  to the precision a double holds them to there.  Every response whose raw
  values show their spread gives the numbers it gave.

* **`kriging_adequacy()` let the unit, the origin and the row order of the
  points choose a lattice's kriging neighbourhoods.**  On gridded locations
  several locations are often the same distance from a cell's centre or a
  held-out point, and `gstat::krige()` took its `nmax` among them by the
  last bits of the coordinates and by the order of the rows.  On a 32 x 24
  lattice at 250 m, the same points and cells in US feet or with the false
  origin moved by 0.1 m moved `kr_pred` in about half of 48 cells, by up to
  0.03, and `kr_shift` by up to 1; on a 30 x 30 lattice the rows reversed
  moved `kr_pred` in 5 of 9 cells and the cross-validation statistic in its
  fourth digit.  The `nmax` are now chosen by the package's tie rule, as its
  other neighbour searches choose them: squared distances within one part in
  ten million count as equal, and a tie goes to the lower first coordinate,
  then the lower second.  A cell or held-out location whose `nmax`-th
  nearest ties with the next is kriged on its own from the locations the
  rule picks; the test for a cell whose own locations reach beyond its
  `nmax` nearest counts with the same tolerance, and the locations added
  outside it are chosen the same way.  The rows are put in order by
  coordinates, cell ID and response before anything is computed, and a
  location is still held out under the fold of its first row as given.  A
  location whose rows name more than one cell was given the cell of its
  first row; it is now no cell's own, with a warning that gives the count of
  such locations, and each of its rows still counts in the cell it names.
  Layers with no tied distances give the same figures, to rounding.  A
  lattice takes longer: one call to gstat per tied cell and held-out
  location.

* **`make_folds()` let the unit and the origin of the coordinates choose the
  block of a lattice point on a block edge, and the neighbours of a
  buffered split at the buffer distance.**  A point on the edge two blocks
  share went to the first block `sf::st_intersects()` found it in, which was
  the lower block only while the edge passed through the point exactly; the
  last bits of the coordinates put it a hair inside one block or the other.
  On a 37 x 25 lattice at 250 m, 3 x 2 blocks put 469 of 925 rows in
  another fold with the false origin moved by 0.1 m, and 4 x 3 blocks 448 in
  UTM 15N in US feet; the selection half of `select_on = "split"` on a 40 x
  30 lattice moved by 786 rows, and with it the `rss` of
  `resolution_profile(select_on = "split")`.  A point within a margin of
  one part in 1e10 of the largest coordinate of the points and the blocks
  from a block's edge is now in that block, and a point in two blocks goes
  to the lower one (for the grid, the cell to its left or below; for
  supplied `blocks`, the lower row, as `?make_folds` said), in every unit
  and origin.  A point within the same margin of `boundary` is inside it:
  lattice points on a `boundary` drawn through them were read as outside
  in other units and origins (8 to 50 of 469 points of a lattice clipped to
  a triangle), so the call warned about them and extended the block region
  to the points' bounding box.  The block counts behind the leakage
  advice read an edge point the same way, and `cv_block_size_sweep()` and
  the cross-validation functions given these folds score on them.  A
  point outside every supplied block, or outside the cells of a clipped
  grid, as near two blocks to one part in ten million went to whichever
  `which.min()` met first in the last bits of the distances (3 to 10 of 72
  such points moved on probe layers in feet, kilometres and moved origins);
  it now goes to the lower row, and the count of points outside the design
  and its warning are unchanged.  `method = "buffered_loo"` excluded a
  neighbour exactly at the buffer distance in one frame and kept it in
  another: on a 20 x 15 lattice with `buffer` 500 m, 295 of 300 training
  sets differed in US feet and 253 with the origin moved.  A squared
  distance within one part in ten million of the buffer's square now counts
  as at the buffer, and such neighbours are excluded in every frame, as they
  were in metres.  Layers whose points sit off the block edges and the
  buffer distance get the same folds as before.

* **`make_folds(method = "nndm")` and `method = "knndm"` let the unit and the
  origin of the coordinates choose the folds of a lattice.**  On gridded
  samples many neighbour distances tie.  NNDM took each point's neighbours
  in the order FNN's kd-tree or `order()` gave equidistant ones, broke ties
  between held-out distances by the exact order of the coordinates, and
  counted, compared with the target and inserted its distances exactly, so
  the last bits of the coordinates settled each tie.  On a 20 x 15 lattice
  at 250 m predicted onto a 30 x 25 lattice at 400 m, 195 of the 300
  training sets in UTM 15N in US feet differed from those in EPSG:26915,
  and 234 with the false origin moved by 0.1 m.  kNNDM handed Ward's method
  and its k-means the coordinates as given and assigned each point to the
  centre FNN or `sf::st_nearest_feature()` returned, so on the same layers
  with `k = 5` 170 and 290 of the 300 rows were in another fold, and the
  `wasserstein` reported was 1191.21, 1223.84 and 1288.05 m.  NNDM now
  orders each point's neighbours by the package's tie rule (squared
  distances within one part in ten million count as equal, a tie goes to
  the lower column, then the lower second coordinate, then the lower row),
  and treats distances within that part of each other as equal in its
  sweep and in its comparisons with the target and with `phi`.  kNNDM's
  Ward clustering and k-means read the coordinates centred on the middle of
  their bounding box, divided by its diagonal and rounded to 2^-30 of it,
  as the `"kmeans"` partition does; a point as near two centres to one part
  in ten million goes to the lower-numbered one, with or without FNN; and
  of the candidates whose W is within that part of the least, the earliest
  is returned.  The held-out distances and W are in the layer's unit as
  before.  `params$max_ecdf_excess`, `params$max_ecdf_deficit` and
  `params$loo_ecdf_deficit` count distances within the same part as equal
  too: with prediction points beside a lattice, prediction and held-out
  distances share lengths, and the step between them changed these figures
  from frame to frame (with the same NNDM folds in all three,
  `max_ecdf_excess` 0.42, 0.42 and 0.413 and `max_ecdf_deficit` 0, 0 and
  0.003).  `params$target_separated`, and the out-of-reach warning it can
  raise on small samples, take a prediction distance within that part of
  a leave-one-out distance as not shorter.  The folds and these figures are
  now the same in every frame.
  They need not be the ones any frame gave before: in metres 62 of the 300
  NNDM training sets exclude other neighbours at the same distance, with
  60670 training rows in all as before, and kNNDM can return another
  candidate.  Layers on which no two distances lie within one part in ten
  million of each other, such as that lattice jittered by up to 60 m, get
  the same folds as before.  On scattered points such pairs occur by
  chance: on 800 uniform points 2 of the 800 NNDM training sets changed.

* **`assign_features_to_polygons()` let the unit and the origin of the
  coordinates choose the cell of a point on a shared edge, and the row
  order of `polygons_sf` choose between polygons its tie-break could not
  tell apart.**  A point on the edge two polygons share was in both only
  while the edge passed through it exactly; the last bits of its
  coordinates put it a hair inside one or the other, and the bounding-box
  centres the smallest-area rule compares were compared exactly, so their
  last bits settled the rest.  On a 37 x 25 lattice at 250 m and a 10 x 7
  grid of square cells over it, 25 of 925 rows went to another cell in UTM
  15N in US feet and 11 with the false origin moved by 0.1 m.  With the
  default `predicate`, a POINT or MULTIPOINT feature is now in every
  polygon one of its points lies within a margin of (one part in 1e10 of
  the largest coordinate of the polygons; in longitude and latitude only
  while `sf::sf_use_s2()` is on), the centres are compared within the same
  margin whatever the predicate and `sf::sf_use_s2()`, and every row goes
  to the same cell in every frame, under `tie_break = "smallest_area"` and
  `"first"`.
  Polygons of one area and one bounding-box centre went to the first in
  the rows of `polygons_sf`, where `?assign_features_to_polygons` said the
  choice does not depend on their order: the two triangles of a split
  square, a square and a rhombus of its area about its centre, one polygon
  under two IDs, in the join by `predicate` and in the largest-overlap
  join.  The rule now compares the mean of the polygons' vertices next,
  lowest, then leftmost, which tells the triangles apart, and a feature
  whose polygons are equal on every key is given none: its ID is `NA`, it
  is counted in `attr(, "ties")` and listed in its new `unresolved`, it is
  kept or dropped as `keep_unassigned` says, and a warning of class
  `spatialkit_ties_unresolved` gives how many such features there are.
  `tie_break = "first"` keeps the first polygon in row order.  Points off
  the edges, and features whose candidates differ in area or centre, are
  assigned as before.

* **`kriging_adequacy()` kriged points XYZ throughout in three
  dimensions.**  `gstat::krige()` uses every coordinate dimension of its
  locations, while the variogram is fitted in x and y, so the heights moved
  every column: on 400 points with heights of 0 to 1000 m, `kr_var` by 13
  percent on average.  The points are now read in x and y on every layer,
  as `estimate_sac_range()` reads them, and a layer XYZ or XYM throughout
  gives the result of its XY copy.

* **Four subsamples were drawn by row position, so their results moved with
  the row order of the input.**  Each draws, or thins, from the rows in a
  fixed order now, as `estimate_sac_range()` and
  `determine_optimal_levels()` order their rows before their draws: by
  coordinates, then the other values that tell rows at one location apart,
  and for the predictor index of `area_of_applicability()` by the
  predictors in turn:
  - `summarize_by_cell(deff = "variogram")` on a cell of more than
    `deff_max_n` points.  One cell of 1500 points had a standard error of
    0.261 with its rows as given and 0.255 with them reversed, and a
    `cell_weight` of 12.8 against 13.4; `..neff_*`, `..df_*` and `..ci_*`
    moved with them.
  - `area_of_applicability()` with more training rows than
    `normalizer_max_n`, for the normaliser of the predictor index and of
    the location index.  1000 training rows at `normalizer_max_n = 300` gave
    a normaliser of 1.018 as given and 1.155 reversed, and DI and the
    threshold with it; with a `threshold` supplied, the inside or outside
    call moved too.
  - `gp_lengthscale_bounds()` above `max_n` distinct locations, and the GP
    prior and basis that `fit_bayesian_spatial_model()` and `cv_bayes()`
    build from it.  1500 locations gave a lower bound of 675.7 as given and
    684.0 reversed.
  - The GP pre-fit check of `fit_bayesian_spatial_model()` and `cv_bayes()`
    on more than 1500 rows, whose Moran's I at the basis floor thins the
    rows at even steps.  On 2400 rows, 10 row orders gave I from 0.031 to
    0.061 at one floor, so the status, and the pre-fit message, was
    `"below_floor_long_range"` in 4 and `"below_floor"` in 6.  The thinning
    still draws no random number.

  Each result is now the same for the input in any row order; under the cap
  nothing is drawn and the results are those of before.  Above it the
  subsample is another one, so the figures differ from before for the rows
  as given.

* **`kriging_adequacy()` kept the last fold's prediction for a location held
  out by several folds.**  Each location is kriged once, so a hand-built
  `folds` whose test sets shared rows wrote over a location's prediction
  with each fold that held it out, and the cross-validation statistic moved
  with the order of the folds: 200 points with rows 1 to 50 in every test
  set gave a variance of the z-scores of 0.978 and 0.989 with the folds in
  the two orders.  Such folds are now refused before anything is fitted,
  with an error that names the location's row and the two folds, as
  `area_of_applicability()` refuses them.  `make_folds()` results and
  vectors of fold labels never share a test row.

* **`kriging_adequacy(folds =)` refused a vector of fold labels with one per
  row of `assigned_points_sf` once a row was dropped for a missing response
  or no location.**  The labels were read against the rows kept, and the
  call stopped, after the variogram fit, with a message saying the labels
  must have one entry per row, which they had.  A vector of labels is now
  read as every `cv_*()` reads one: one label per row of
  `assigned_points_sf`, the i-th for the i-th row, the labels of dropped
  rows unused.  One label per row kept is read as before.  A vector of any
  other length stops before anything is fitted, with an error that gives
  its length, the number of rows and, when rows were dropped, the number
  kept.

* **`residual_morans_i()` read its response, design and residuals at their
  origin and unit.**  It rebuilt the design and regressed the response on
  the raw columns to test whether the residuals are OLS residuals, and gave
  the same design to the residual moments.  A predictor at 1e8 plus a spread
  of 1 was dropped by the QR, so the moments had rank 1, and a response at
  1e8 left the rounding of the offset in the reference residuals.  On 300
  points with OLS residuals that are the same at every offset, `null =
  "auto"` used the randomisation null with nothing said (`exact` `FALSE`, z
  28.24 where offset 0 gives 28.25), and `null = "residual"` logged that the
  residuals were not OLS residuals.  The tolerance of that test was 1e-7
  times `max(sd(y), .Machine$double.eps)`, so for a response whose standard
  deviation is below about 1e-22 it passed residuals that are not OLS
  residuals.  The test for residuals with no spread was an absolute floor on
  their sum of squares: residuals of real spread whose sum of squares was
  below `.Machine$double.eps` (a standard deviation of 1e-9 on 200 rows)
  were refused with the warning "degenerate weights or zero residual
  variance", and `compare_models()` gave `NA` Moran's I columns for them,
  while the residuals of an exact fit multiplied by 2^20 were scored.  Now
  the response is taken about its centre and the predictors at an offset
  about theirs, as in the other OLS fits, and the moments read that
  design.  Residuals pass as OLS residuals within 1e-7 of the response's
  standard deviation, read at any size, of the OLS residuals about the
  centres or of those `lm()` gives on the raw columns where the raw design
  keeps the rank, or within the rounding a regression of the raw columns
  leaves, so a backend that regressed the raw response at an offset passes
  too, on any number of rows.  Residuals are refused where they do not vary beyond
  the rounding of the response, or where they are the OLS residuals of an
  exact fit, the bound `resolution_profile()` and `estimate_sac_range()`
  read; the warning names that cause, and degenerate weights have their
  own.  The sums are taken on the residuals divided by a power of two, so I,
  z and the p-value are the same in every unit.  `compare_models()` and the
  residual check of `fit_bayesian_spatial_model()` read the same function.  A
  fit whose columns show their spread, with residuals above the bound, gives
  the numbers it gave, except that the kurtosis, and through it the
  randomisation variance, can differ in the last bit: R takes the fourth
  power with `pow()`, which does not commute with the division.  The
  rounding bound reads the offset the fit regressed the response at, so a
  GWR that `fit_gwr_model()` took about the response's centre, at a
  working bandwidth with the response at 1e15, is no longer read as OLS:
  it took the residual null with `exact` `TRUE` where offset 0 gives the
  randomisation null.

* **The k-nearest-neighbour weights tied neighbour distances within 1e-9
  map units on any layer.**  The weights behind `residual_morans_i()`,
  `estimate_sac_range()`, the Moran tests of `determine_optimal_levels()`
  and `resolution_profile()` and the range check of
  `fit_bayesian_spatial_model()` took two distances as one when they differed
  by less than 1e-9 of the coordinates' size, and that size was floored at
  1 map unit.  A layer whose coordinates span less than about 1e-5 map
  units therefore got other neighbours than the same layer in metres, more
  of them as the unit shrank, up to every other point with equal weight:
  on 120 points in a 5 km square scaled by 2^-50, Moran's I of the
  residuals was -0.008 with p 1 where the layer in metres gives 0.277, and
  `estimate_sac_range()` returned `NA` with "no spatial structure detected
  at the nearest-neighbour scale" where it finds a range.  Ties now read
  the distances alone (see the next entry), so the weights, I, z, p and the
  range scale with the coordinates.  `compare_models()`
  had the same floor where it checks that fits with an information
  criterion share their rows, on the coordinates of a layer with no CRS and
  on the response: rows whose coordinates differ by less than 1e-9 units,
  or responses that differ by less than 1e-10 (a response at 1e-12 and 1.5
  times it), passed as the same, and the AICc or LOOIC column was kept.
  The check that `folds` passed to a `cv_*()` function,
  `area_of_applicability()`, `fold_separation()` or `kriging_adequacy()`
  were built on the data had the same floor on a layer with no CRS: folds
  built on other data whose rows sit less than 1e-9 units from this data's
  rows were applied in silence.  Both checks dropped the floor (their
  margin is now a part of the spread; see the entry on these checks at a
  large offset).  The margin of one part in 1e10 of the
  coordinates, by which `build_tessellation()` and `cell_values_at()` read
  seeds as the centres a count record holds and
  `assign_features_to_polygons()` reads a point as on the edge between
  scored cells, had a floor of 1e-10 units: on a layer with no CRS whose
  centres are 1000 units apart and whose seeds lie 10 units from them, all
  scaled by 2^-40, the seeds were read as the centres, and the cells were
  drawn around the centres in place of the seeds.  The margin is now one
  part in 1e10 of the coordinates in any unit, as the help of
  `assign_features_to_polygons()` says.

* **The k-nearest-neighbour weights moved with the origin of the
  coordinates.**  The weights behind `residual_morans_i()`, the residual
  Moran's I of `compare_models()`, the structure test of
  `estimate_sac_range()` and `sac_nugget()`, the Moran tests of
  `determine_optimal_levels()`, `resolution_profile()` and
  `select_resolution(criterion = "moran_z")`, and the range check of
  `fit_bayesian_spatial_model()` took two neighbour distances as
  one when they differed by less than 1e-9 of the largest coordinate: 4 mm
  at a UTM northing, 1000 m at an offset of 1e12 and 1e6 m at 1e15.  The
  same points moved by a whole number of map units therefore got other
  neighbours.  On 80 points 2 km across, each point had 46 to 79 neighbours
  in place of 8 at an offset of 1e12; at 1e15 every point was every other
  point's neighbour, Moran's I was exactly -1 / (n - 1), and
  `estimate_sac_range()` returned `NA` with "no spatial structure detected
  at the nearest-neighbour scale".  Two distances are now tied when their
  squares are within one part in ten million of each other, the rule of the
  package's other neighbour searches, with no term in the size of the
  coordinates.  The weights, and every I, z and p read from them, are the
  same at any origin and in any unit wherever the doubles hold the
  coordinates to better than about one part in 1e8 of the distances
  between neighbours.  A lattice ties as before.  On a layer in UTM metres,
  distances less than 4 mm apart but more than that part apart no longer
  share a slot, and a point a few units in the last place from two
  co-located points no longer shares their slot at k = 1.

* **The checks that folds and fits describe the same data passed other data
  at a large offset.**  On a layer with no CRS, the check that `folds`
  passed to a `cv_*()` function, `area_of_applicability()`,
  `fold_separation()` or `kriging_adequacy()` were built on the data, and
  the check by which `compare_models()` keeps AICc or LOOIC only for fits
  of the same rows, took two locations as one when they differed by less
  than 1e-9 of the largest coordinate, and `compare_models()` took two
  responses as one when they differed by less than 1e-10 of the largest
  response.  At an offset of 1e12 that is 1000 units: folds built on a
  layer were applied to the same layer moved by 500 units, which offset 0
  refuses, and a response near 1e12 and the same plus 50 were read as one
  response.  The margins are now a part of the spread of the values plus
  eight units in the last place of their size, so they cover the rounding
  of the values and do not grow with the origin.

* **`fit_gwr_model()` read a predictor or a response at a large offset on
  the rounding of the offset.**  The global least-squares fit the GWR is
  read against was `lm()` on the raw columns, which drops a predictor whose
  spread is under about 1e-7 of its size (1e8 plus a spread of 1, a
  timestamp in seconds) as collinear with the intercept, in silence.  On
  300 points with a predictor of unit spread, `info$global_ols` and
  `summary()`'s line "Global least squares on the same rows" gave an AICc
  of 881.32 at offset 0 and 1364.41, that of the intercept alone, from
  1e8.  The local regressions divided a predictor by its standard
  deviation before they centred it, so from about 1e10 the slopes, fitted
  values and AICc moved with the predictor's origin: at 1e15 the slopes by
  2 percent and the AICc from 489.16 to 496.03.  A response at 1e15 plus a
  spread of 1 was regressed on its raw values, which a double holds to an
  eighth: the residual sum of squares was 111.4 where the same values at 0
  give 73.6, the bandwidth search chose 71 neighbours where it chose 41,
  and `gwr_model_selection()`'s table and the residuals `residuals()` and
  `residual_morans_i()` read were those of the same rounding.  The
  collinearity survey centred a predictor at 1e15 about its rounded mean,
  so `info$condition_index` and the slope index `cn_slopes` of
  `info$local_collinearity` moved with its origin (a maximum `cn_slopes` of
  2.05 where the same values at 0 give 2.15).  A predictor whose centre
  is more than 1000 times its spread is now taken less its centre before
  the bandwidth search, the local regressions, the global fit and the
  collinearity survey, and so is a response whose raw values, or whose
  least-squares residuals on the predictors, lose their spread to
  rounding, as `estimate_sac_range()` takes it.  The centres are added
  back only to the fitted values, the local intercepts, the intercept and
  fitted values of `fit$engine$lm` and `predict()`'s values, and
  `predict()` takes a new row's predictors about the training centres.
  Residuals, slopes and their standard errors, AICc, the bandwidth chosen,
  the global fit, the condition indices of the slopes and the selection
  table are now those of the same values at 0, to rounding; fitted values
  and predictions of a response at 1e15 are in its units, to the eighth a
  double holds there, and `residuals()` is the centred response less the
  centred fitted values.  `model_metrics()`, `summary()`,
  `evaluate_insample()`, `compare_models()`, `cv_gwr()`, `cv_spatial()`
  with a GWR learner and `select_features_forward()`'s `score_holdout`
  scored such a fit on its fitted values and predictions in the response's
  units, which carry that eighth: on 300 points at 1e15 RMSE was 0.49555
  where `residuals()` give 0.49543, `cv_gwr()`'s RMSE 0.65316 where the
  values at 0 give 0.65533, and R-squared was `NA`.  They now score it on
  the response, the fitted values and the predictions less the centre, so
  RMSE, MAE and R-squared are those of the values at 0, and MAPE and SMAPE
  divide the same errors by the values in the response's units.  The
  `predictions` of a cross-validation whose fold fits took the response
  about a centre gain the columns `y_centre`, `yhat_centred` and
  `y_train_mean_centred`, which the pooled scores read; `yhat` stays in
  the response's units.  `fit$engine$lm` is solved on the centred columns
  and returned as the fit at the data's origin, with the QR factor of the
  raw columns, so `summary()`, `vcov()`, `confint()` and `predict()` on it
  give the intercept at 0 and its standard error; with a predictor more
  than about 1e7 times its spread from 0, `predict(se.fit = TRUE)` on it
  stops with "singular matrix".  `fit$engine$centres` gives the centres.
  Every other predictor is fitted as it was, to the last bit (the next
  entry gives the response), and a predictor that is constant, or a linear
  combination of the others, is still left out of the global fit.

* **A GWR at a global window on a response at 1e8 or 1e10 failed
  `residual_morans_i()`'s check for OLS residuals.**  `fit_gwr_model()`
  took a response about its centre only where its raw values, or its
  least-squares residuals on the predictors, lose their spread to rounding,
  so 1e8 plus a spread of 1 was regressed on its raw values.  The local
  solve's residuals then carry the rounding of the offset, which grows with
  the rows.  On 6000 rows with `bandwidth = 6000`, `adaptive = TRUE` and
  `kernel = "boxcar"`, which is global OLS, the residuals missed the check
  at 1e8 and 1e10, and `residual_morans_i()` with `null = "auto"` used
  randomisation, with `exact = FALSE` and a z of 126.387931 where the same
  values at 0 give the residual null and 126.391564, in silence.  A
  response at -1e8, and a model of the intercept alone, failed the same
  way.  `fit_gwr_model()`, `predict()` and `gwr_model_selection()` now also take
  the response less its centre where that centre is more than 1000 times
  its spread, the gate a predictor already has.  The global GWR then gives
  the OLS residuals, the residual null and the z of the values at 0 at any
  offset, and `residuals()` is the centred response less the centred
  fitted values.  A response whose centre is within 1000 times its spread
  is fitted as it was, to the last bit.

* **`select_features_forward()`'s training-fold-mean baseline read the
  rounding of a response at a large offset.**  Where the learner cannot
  fit a model with no predictors (a forest, or a GWR `fit_fn` that refuses
  `character(0)`) and no `null_fn` is given, the baseline is the mean
  response of each fold's training rows.  It was scored on the response
  and those means in the response's units, which a double holds at 1e15
  to an eighth.  On 300 points with a GWR learner the baseline's RMSE was
  2.4737792 at 1e15 and 2.4701674 at 1e12, where the same values at 0
  give 2.4701727, so `baseline$score`, the step-0 rows of `history` and
  of `fold_scores`, which the "se" rule reads, and the figure in the
  message "the baseline is the training-fold mean" moved with the
  response's origin (2.4521 at 1e15 against 2.4480 at 0), while every
  candidate set scored the same at both.  The baseline is now scored on
  the response less a centre, with each fold's training mean taken about
  the same centre: the centre of the GWR fit where a fold's fit took the
  response about one, and otherwise, where a fold's raw training response
  loses its spread to rounding, the centre `estimate_sac_range()` takes.
  Its RMSE and MAE are then those of the values at 0.  Where no fold has
  a centre the baseline is scored as before, to the last bit, and
  `cv_spatial()` returns what it returned for every call.

* **`fit_gwr_model()` logged "predictor has near-zero variance; GWR may be
  unstable" for a predictor in small units.**  The test was a standard
  deviation under 100 times `.Machine$double.eps`, so a predictor of unit
  spread multiplied by 2^-46 had the note and the same predictor at 2^-40
  did not, though the fit is the same at every scale.  Where the standard
  deviation is under that floor, the column is now read again on its own
  scale, and the note is logged only for a column without spread beyond
  the rounding of its values.

* **The range check of `fit_bayesian_spatial_model()` before sampling read
  its Moran's I tests at the origin of a column, and scored the rounding of
  an exact fit.**  The short-lag and floor-distance tests regressed the raw
  response on the raw predictors.  A predictor at 1e8 plus a spread of 1
  was dropped by the QR, and a response at 1e15 left the rounding of the
  offset in the residuals: on 600 points, the short-lag Moran's I went from
  0.770 to 0.177 with the predictor at 1e8 and to 0.023 (p 0.18) with the
  response at 1e15, and a layer whose range is below the basis floor went
  from "below_floor" to "below_floor_no_autocorrelation".  The only test
  for residuals without spread was a sum of squares of 0, so where
  `estimate_sac_range()` refused the layer as an exact fit of the response
  on `predictor_vars`, the check still computed Moran's I on the rounding,
  and could log in the same line that the response "is autocorrelated at
  short range".  Now the predictors are fitted about their centres where
  the centre is more than 1,000 times the spread, as in the other OLS fits,
  and the response by the same rule; the moments read that design.  A
  response whose residuals lose their spread through a trend is fitted
  again about its centre, as `estimate_sac_range()` fits it.  Residuals within the rounding
  of an exact fit, the bound `estimate_sac_range()` and
  `residual_morans_i()` read, give no test (`$info$gp_prefit$moran` and
  `moran_floor` are `NULL`), and the log says nothing about
  autocorrelation.  The sums are taken on the residuals divided by a power
  of two, so I, z and the p-value are the same in every unit.  Every other
  layer gives the numbers it gave, to the last bit.

* **`fit_bayesian_spatial_model(standardize_predictors = TRUE)` did not
  standardise a predictor in small units.**  A predictor's scale was its
  standard deviation floored at `.Machine$double.eps`, so a predictor with a
  standard deviation below 2.2e-16 was handed to `brms::brm()` with a
  standard deviation of 0.063 (at 1.4e-17) down to 3.8e-6 (at 8.5e-22), while
  the slope prior is set for a predictor with a standard deviation of 1; the
  slope and its posterior depended on the predictor's unit.  The scales of
  the coordinates had the same floor.  Now a predictor or a coordinate whose
  standard deviation is below `.Machine$double.eps` is read again: one with
  spread is scaled by its standard deviation, taken at any size.  One
  without spread got a scale of `.Machine$double.eps`, which turned the
  rounding of a constant column into values of order 1; now a predictor
  without spread gets a scale of 1, and a coordinate axis without spread
  (a transect along the other axis) takes the other axis's scale, so the
  basis's resolution and the floor of the range check are those of the
  axis with extent.  The scale used is stored in
  `$info$predictor_scaling` and `$info$coord_scaling`, and `predict()`
  reads it there.  A predictor or a coordinate axis whose values differ
  only by their rounding was read as one without spread only while its
  standard deviation was below the floor: at 4.5e6 such a predictor was
  divided by a scale of 7.5e-10 and handed to `brms::brm()` with a standard
  deviation of 1, and a transect at x = 4.5e6 with jitter in the last bit
  of x got an x scale of 7.7e-10.  Such a column is now read as without
  spread at any size.  A predictor or a coordinate with spread whose
  standard deviation is at least `.Machine$double.eps` is scaled as
  before.

* **The print of a profile of one level read "1 levels" and "WSS rises at 0
  step(s)".**  Every default profile has one level.  The print now reads "1
  level", and leaves out the count of steps the WSS curve rises at when
  fewer than two levels have a WSS, since there is no step between levels
  to count.  A subset of a profile's rows printed the count of the whole
  ladder; it now counts the steps between its own levels.  `plot()` of a
  profile of one level has the subtitle "1 level of 133 cells" where it had
  "1 levels from 133 to 133 cells".

* **Messages asked for FNN, which installs with spatialkit, and the
  work-limit errors printed two equal numbers.**  The texts that appear when
  FNN cannot be loaded (the work limits of the exact searches of
  `resolution_profile()`, `summarize_by_cell()` and `cell_values_at()`, the
  density kind's log line, `area_of_applicability(use_fnn = TRUE)` and the
  weights of Moran's z) said to install it, and `spatialkit_moran_skipped`
  said that it is not installed; they now say that it could not be loaded
  (and is to be reinstalled).  The weights' error asks for Matrix alone when
  FNN is loaded.  The work and the limit are printed in full, with as many
  digits as tell them apart: a work 0.04 percent over the limit read "1e+10
  ... above the limit of 1e+10".  The fit guard of `resolution_profile()`
  named "the FNN package installed" in a way on; it now says that the random
  starts of the density kind need FNN only where FNN could not be loaded.

* **Texts of `resolution_profile()`, `select_resolution()`, `summary()` and
  the seeding functions that misstated the case.**  The note on a count that
  one mean matches said "within 2% of the optimum" also where one mean was
  lower; it says "is not above the optimum by more than 2%", which is the
  rule.  After Cp's note, the note of a variogram that found no structure
  was pasted in lower case with its close said twice; it is a sentence of
  its own.  At the support ceiling the hint said to lower `min_cell_n`,
  which moves every level of the ladder; it names `levels` that run past the
  count.  A selection by size called the default 4.5 "asked", a selection on
  one level read "(1 of 1 levels)", the summary of one level "(133 to 133
  cells)" and a profile of one start "1 k-means++ restarts".  `summary()`
  said that one cell does as well where the cells beat one cell and only the
  F test (or, for `cp`, the profile's variogram) did not show structure; it
  now says which, and the result then carries `attr(, "not_identified")`.
  On a constant response the warning left out `loo`, `loo_se` and
  `own_weight`, which are `NA` too, and `summary()` stopped with an error
  that named neither the constant response nor `loo`; both now name them.
  The print of a profile whose nugget is at the fit's lower bound, with no
  level to take a noise variance from, marks the nugget and says that `cp`
  is `NA` at every level, as the warning does.  The sizes of a ladder held
  for the same reason raise one `spatialkit_size_not_reached` (Breaking
  changes has the fields).  The numbers of one printed line or note share
  one notation, where "90200000 (one cell for the layer: 1.09e+08)" mixed
  two, and the summary's closing line reads "every pick is 60 cells" where
  it read "picks span 60 to 60 cells (1.0x)".  `voronoi_seeds_kmeans()` and
  `voronoi_seeds_random()` refuse the `NA` of `determine_optimal_levels()`
  with the note `get_voronoi_seeds()` gives, which names `on_no_elbow`.

* **Texts of `summarize_by_cell()`, `cell_values_at()` and
  `assign_features_to_polygons()` that misstated the case.**  The
  no-structure warning said of predictors that a map of their pooled values
  looks smooth "not because the response varies smoothly"; it says "the
  column" or "the columns" unless the response alone is named.
  `cell_values_at()` on a logical, factor or character column asked for a
  summary that `summarize_by_cell()`, which summarises numeric columns only,
  cannot make; it says how to convert the column, and for a factor that
  `as.numeric()` gives the level codes, not the labels.  Features with empty
  geometry were said to fall outside the polygons, with advice on their
  cover or the CRS; the warning counts them and says that they have no
  location to assign.

* **Texts of the folds, the variogram, the seeding, the fits and kriging
  that misstated the case.**  `make_folds(auto_range = TRUE)` with a
  `block_size` said that the range was used as the minimum block size, where
  the `block_size` given is; it says so, with the grid of blocks that gives,
  and `?select_features_forward` no longer calls `block_size` a floor under
  the range.  On one row `make_folds(method = "random_kfold")` wrote "1
  folds" and offered fewer folds.  `get_voronoi_seeds(method = "random")` or
  `"provided"`, given a density setting under the k-means kind, said that
  the density kind would read it; it says that the method does not.
  `model_metrics(se.fit = TRUE)` on a fit that is not Bayesian pointed to
  `cv_bayes()`.  The `rejected_reason` of a variogram with fewer than three
  lag bins holding pairs said "(singular fits)" (Breaking changes has the
  new text).  `detrend = "reml"` without nlme said that the fit did not
  converge; it says that it was not fitted, and why (on fewer than 30
  locations no fit is tried: see the entry on the floor of 30 locations).
  `kriging_adequacy()`'s print gave a pure nugget a gain with no word on it;
  it says that there the gain is that of averaging the locations a cell is
  kriged from (`kr_n_used`), its own among them, over averaging its own
  alone.

* **`make_folds()` advised fold counts that give the same folds or none,
  and kNNDM warned of folds it then refused to build.**  Where it cut `k` to
  the number of rows (`"random_kfold"`, `"knndm"`) or of groups
  (`"leave_location_out"`), the `spatialkit_k_reduced` warning advised `k <=`
  that number, which gives the same one-per-fold folds: 5 rows and `k = 10`
  said "pass k <= 5".  It now names one fewer ("pass k <= 4"), and with 2
  rows or 2 groups, where no smaller `k` exists, it gives no advice.  Under
  `"knndm"` it names the range from the smallest `k` that `min_train`
  allows: 8 rows, `k = 10` and `min_train = 0.85` say "pass k = 7 for
  fewer, larger folds (min_train = 0.85 needs k >= 7)", and where no `k`
  below the number of rows meets `min_train` it gives no advice.
  `make_folds(method = "knndm")` with more folds than rows and a `min_train`
  above `1 - 1/n` warned "returning 8 folds" and then stopped with an error
  that ended "Use min_train <= 0.875, or more folds.", though no more folds
  were possible.  It now tests `min_train` before it warns, so such a call
  stops with the error alone; the error names the `k` asked for and the
  rows it was cut to, and offers more folds only where one fold per row
  would meet `min_train` ("Use min_train <= 0.875." otherwise).  The folds
  returned, `k` and the condition's fields are unchanged.

* **Cp took a nugget at the fit's lower bound as its noise variance, and one
  that the closest pairs of points contradict.**  `resolution_profile()`
  read Cp's noise variance off the variogram's nugget, and a nugget of 0 (or
  of at most 1e-4 of the sill), which `estimate_sac_range()` marks as
  `nugget_at_bound` since it is the bound and not an estimate, left Cp no
  penalty: it fell with every finer level to the ceiling whatever the field,
  with an unclassed warning.  Such a nugget is no longer taken: Cp's noise
  variance comes from the finest level whose cells hold two scored rows on
  average, as when the variogram gives no nugget, and the range still serves
  `reliability`.  Whatever its source, the noise variance is now held to the
  semivariance of nearest-neighbour pairs the variogram carries
  (`attr(sac, "short_lag_semivariance")`), which bounds a nugget from above:
  a nugget or a residual mean square above it is replaced by it.
  `attr(, "cp_noise")$source` names the source used, `"short-lag
  semivariance"` among them, and under the hold the record adds `cap`, the
  semivariance, and `uncapped`, the source, value and level it replaced; the
  print's `variogram` line says which.  The `cp` column, its value at one cell and
  Cp's count change wherever either applies.  A profile on a ladder of
  sizes under `sac = "auto"` fits no variogram and is unchanged.

* **On gridded locations the nearest-neighbour pairs of
  `estimate_sac_range()` depended on whether FNN could be loaded.**  Each
  location is paired with its nearest other location, and on a lattice two
  to four are equally near.  FNN's kd-tree and the package's own search took
  different ones, so `short_lag_pairs` and `short_lag_semivariance`, and
  with them Cp's noise variance where that semivariance holds it, the `cp`
  column and Cp at one cell, differed between a session with FNN and one
  without (on a 12 x 12 lattice, 80 pairs against 143), and moved with the
  unit of the coordinates.  Of several locations whose squared distances
  agree to within one part in ten million, the pair now takes the one with
  the lower x coordinate, then the lower y, with or without FNN, in any
  order of the rows and in any unit or origin.  On a square lattice that
  pairs almost every location with its neighbour along x, so the
  semivariance measures the field along x.  Where the ties are exact, as on
  a lattice in whole metres, that is the pick the search without FNN made,
  so on gridded layers a session with FNN gets other pairs and other `cp`
  values than before, and a session without FNN does where rounding had
  broken the ties.  Layers on which each location's nearest other location
  is nearer than the next by more than that part in squared distance are
  unchanged.

* **The two variogram warnings of `resolution_profile()` had no class.**
  The warning that the variogram gives no usable range now has class
  `spatialkit_variogram_unusable` (fields `reason`, `converged`, `nugget`,
  `nugget_used`), and the warning on where Cp's noise variance came from
  has class `spatialkit_cp_noise` (fields `source`, `value`, `level`, `cap`,
  `uncapped_source`, `uncapped_value`), so a script can catch each by its
  class.  The warning that the nugget is 0 is no longer raised, since Cp
  does not take such a nugget: `spatialkit_cp_noise` says where its noise
  variance came from instead.  Points dropped and a variogram of another
  variable are still warned about without a class.

* **A count chosen by `cp` carried no record of how it was set.**  On a
  profile without an estimate (one with `predictor_vars` or a logical
  response, or one built with `estimate = NULL`), `select_resolution(criterion
  = "cp")` returned no `count_rule`, so `summarize_by_cell(estimate =
  "pooled")`, the default, pooled on cells whose count was chosen for plain
  means and said nothing.  A cp selection now carries the record
  (`estimate = "mean"`, no `scored` list), and so do its seeds, the cells
  built on them and the assigned layer, and a summary of pooled values on
  them raises `spatialkit_estimate_mismatch`.  The selection gains
  `count_rule` and nothing else.  A profile without an estimate passed
  straight as the count is read by `cp` and carries it too.  Selections by
  `reliability`, `elbow` and `moran_z` off such a profile still carry none.

* **The seeds of a profile without an estimate fitted to a subsample were
  the subsample's centres.**  On a ladder of counts the centres are fitted
  to `sample_n` (1,500) points, and `select_resolution()` handed the chosen
  level's centres on as they were, so on a larger layer the cells built on
  them were fitted to a few rows each, far from the cells a fit to every
  point gives (their within-cell sum of squares over every point is
  higher).  Such a profile now keeps every point's coordinates
  (`attr(, "refit_rows")`), and `select_resolution()` refits the level it
  chooses to every point before its seeds are handed on: the kind's own
  iteration (Hartigan and Wong's for `"kmeans"`, Lloyd's on every point's
  density weights for `"density"`) from the centres scored, with the
  profile's iteration limit and no random draw; where Hartigan and Wong's
  algorithm refuses the centres, Lloyd's iteration with every weight one
  runs in its place, with a log line.  The seeds carry
  `attr(, "refit")` (`rows`, `subsample`, `kind`, `iter`, `iter_max`,
  `converged`, `wss_scored`, `wss`), the seeding functions hand them on, and
  the selection's print says so.  The count is still the one the subsample
  gave, and the cells are not the ones scored: the record on the seeds says
  `scored_cells = FALSE`.  The refit costs one fit of that level to every
  point; `summary()` and `plot()` of a profile, and `build_tessellation()`
  given a profile as the count of a lattice, refit nothing, and a refit that
  would pass the work limits of a profile's own fits is not run (a logged
  warning says so), and a profile saved by an earlier version, which keeps
  no coordinates, is not refitted (a log line says so).  The centres the
  profile scored are `attr(profile, "centres")[[as.character(sel$best)]]`.
  A profile with an estimate is not refitted: it scores the estimate on every
  row of the layer in the cells of the centres fitted to the subsample, so
  those centres are its seeds, and `summarize_by_cell()` checks the cells
  built on them against its record.

* **The print of a profile fitted to a subsample called density cells
  "k-means fitted to a subsample".**  The header now names the kind:
  "density-weighted k-means fitted to a subsample of 1500" for the density
  kind, the default, and "k-means fitted to a subsample of 1500" for
  `partition = "kmeans"`.

* **`resolution_profile(select_on = "split")` returned a `loo` read off no
  row, with nothing said, when no row of the selection half had a
  response.**  The estimate is scored on the selection half, and the check
  that `estimate` has rows to read is made on the layer, so a layer whose
  rows with a response all lay in the estimation half passed it, under a
  named estimate and under `"auto"`, which became `"pooled"`.  `loo` and
  the estimate's other columns are still `NA` at every level, and a warning
  of class `spatialkit_no_scored_rows` now gives the reason; the condition
  carries `reason`, `n`, `selection`, `n_response` and `column`.  The
  profile's print says "loo NA: no row of the selection half has a
  response" where it said "fewer than two rows", and so does the print of
  the count record, which read "in NA cells".  `select_resolution(criterion
  = "loo")`, the criterion a profile of plain means is read by, stops with
  the same reason where it said that `loo` needs a level at which some cell
  holds two rows with a response, and `summary()` of the profile, whose
  error said that "cp needs a response", gives it too.  No number changes.

* **`resolution_profile(select_on = "split")` without an estimate returned
  `rss`, `cp`, `moran_i` and `moran_z` `NA` at every level, with nothing
  said, when no row of the selection half had a response.**  The warning of
  class `spatialkit_no_scored_rows`, which gave the reason only with an
  estimate, is now raised without one too, and with `predictor_vars` it
  counts the rows with a response and every predictor.  It is not raised
  when no row of the layer has a response.  The split the profile carries
  records the rows of the selection half with a response (`n_response`).
  `summary()` of the profile, which said that "cp needs a response", and
  `select_resolution(criterion = "cp")` and `"moran_z"`, and `"auto"` where
  it reads `cp`, which said what each criterion needs, now give the reason,
  with or without an estimate.  On a size ladder that fitted no variogram,
  `"cp"`, and `"auto"` where it reads `cp`, give the size ladder's reason,
  as before.  No number changes.

* **The warnings of `make_folds(method = "block_kfold")` on blocks smaller
  than the autocorrelation range said "Spatial CV may leak correlated
  information across folds".**  That read blocks below the range as a
  defect whatever the map.  A cross-validated error is the error at the
  distance the folds put between the held-out points and their training
  data, and it reads low only for predictions made farther from the data
  than that.  The five warnings now say so, name
  `fold_separation(prediction_points = )` as the comparison, and say that
  wider blocks give the error farther from the data before their advice on
  the block size, which is unchanged.  The subtitles of `plot_folds()` and
  of `plot()` on a `cv_block_size_sweep()` result no longer say that blocks
  below the range leak, and the sweep's caption and progress message no
  longer call random folds "the leaky reference".  When no range can be
  estimated (no gstat, or fewer than 30 points), the line `make_folds()`
  logs says that the check of the blocks against the range is skipped,
  not a "leakage check".

* **Three run-time texts said which way an error errs without saying where
  the model will predict.**  `print()` on an `rf_fit` said the out-of-bag
  error "is optimistic under spatial autocorrelation"; the caution
  `fit_rf_model(include_coords = TRUE)` logs said the out-of-bag error and
  random cross-validation "overstate its accuracy"; and the caution
  `select_features_forward(method = "random_kfold")` logs said random inner
  folds "leak spatial autocorrelation" and that a blocked outer loop reports
  "honest-looking numbers for a dishonestly chosen feature set".  They now
  read as `?fit_rf_model` and `?select_features_forward` do: the out-of-bag
  error and random folds are the error of predicting among the samples,
  which is the error of a map over an area the samples cover evenly and
  reads low for predictions made farther from the data, so a map that
  predicts farther is scored with `cv_rf()` on folds at its distances;
  random inner folds choose the variables that predict best among the
  samples, and a blocked outer loop scores the chosen set at the blocks'
  distance, for which it was not chosen.  A forest with no row out of bag
  (`replace = FALSE`, `sample_fraction = 1`) has no out-of-bag error, and
  its print no longer describes one: it says to score the forest with
  `cv_rf()` on folds at the distances the map will predict at.  No result
  changes.

* **The `spatialkit_cells_differ` warning named a cause the tie rule
  removed.**  With `what = "loo"` (the same rows in the same number of
  cells, another leave-one-out error) it said that rows exactly on a
  boundary between two cells (gridded locations) can be assigned to the
  other cell.  `assign_features_to_polygons()` gives such a row the cell
  the profile scored it in on the cells the profile scored, kept in the
  CRS it measured the rows in and carrying the record of their count, for
  rows that have a CRS exactly when the record names one.  The warning now
  names the cases that remain: rows the polygons alone placed (rows without
  a CRS beside a count record that names one or the reverse, rows assigned
  to cells without their count record, such as cells read back from a
  file, and rows given their cells by another join),
  cells built in, or transformed to, another CRS than the one the profile
  measured the rows in, cells that are not the Voronoi cells of the
  centres the profile scored (seeds moved, cells edited), and coordinates
  that reached that CRS by another transformation.  `?summarize_by_cell`
  names the same cases.  No result changes.

* **`coerce_to_points()` crashed R on an empty line feature.**  sf's
  `st_cast()` turns an empty MULTILINESTRING into one empty LINESTRING, not
  zero parts, and `st_line_sample()` on it segfaulted and took the session
  with it.  A null geometry in a line layer loads exactly like this from a
  GeoPackage or a shapefile, and `prep_model_data()`, every `cv_*()` function
  and `fold_separation()` go through this path.  GEOS's
  `st_point_on_surface()` crashed the same way on a line feature holding an
  empty part beside real ones, which `ensure_projected()` reached on ordinary
  lon/lat input.  Empty parts are now removed before either call, and an
  empty feature becomes an empty POINT in its own row, which
  `prep_model_data()` and `make_folds()` drop like any empty geometry.  An
  empty LINESTRING used to raise an error; it now gives an empty POINT like
  every other geometry type.

* **`cv_gwr(parallel = n)` hung forever once a GWR had been fitted in the
  session.**  GWmodel is built with OpenMP, and GNU libgomp is not fork-safe:
  after `fit_gwr_model()`, a sequential `cv_gwr()` or a bare
  `GWmodel::bw.gwr()`, the forked `mclapply()` workers blocked on a futex and
  never returned, with no timeout.  That is the ordinary fit-then-validate
  order on Linux.  `cv_gwr()` now runs its folds one after another whenever
  `parallel` would fork, and says so in a warning.  This gives up the
  speed-up parallel folds had in a fresh session; the other `cv_*()`
  functions still fork, since ranger does not use libgomp.  A `cv_spatial()`
  `fit_fn` that calls GWmodel can still hang with `parallel`, which its help
  page now says.

* **When some folds failed, `overall` quietly left them out.**  A partial
  failure was only logged, so `overall` pooled the surviving folds with no R
  condition, and the folds that fail are usually the hardest to predict (a
  region or a factor level no training fold covers).  `cv_gwr()`,
  `cv_bayes()`, `cv_spatial()` and `cv_rf()` now warn, naming each failed
  fold and how many rows `overall` covers.  Folds dropped before fitting are
  already warned about and are not counted twice.

* **`select_features_forward()` could pick a variable for making folds
  fail.**  Each candidate was scored on whatever rows its CV run predicted,
  so a set whose fit failed on a fold (a factor level found in one block, an
  ordinary case under block folds) was scored on fewer, easier rows.  A
  pure-noise factor beat the true driver: RMSE 2.35 on 192 rows against 2.63
  on 250, although the driver scored 1.73 on those same 192 rows.  Every
  candidate is now scored on one fixed row set (the rows the null model
  predicts, or, with no null model, the rows the step-1 sets predict), a set
  that leaves any of them unpredicted scores `NA` with a warning, and
  `history` gains `n_pred`.

* **`compare_models_cv()` ranked models scored on different rows.**  Shared
  folds guarantee the same splits, not the same scored rows: when one
  backend lost a fold or returned `NA` predictions its `overall` row pooled a
  subset, and the help page called the columns comparable.  A fixed-bandwidth
  GWR scored on 158 of 200 rows ranked above a random forest scored on all
  200, although the forest was 45 percent better on the 158 rows both
  predicted.  When the predicted row sets differ, every model is now rescored
  on the rows all of them predicted, with a warning, `overall` carries
  `n_pred`, and the all-rows numbers stay in `attr(overall, "all_rows")`.
  The per-fold table and the Bayesian coverage columns are not rescored.

* **`determine_optimal_levels()` read an elbow into points that have none.**
  The elbow was the level furthest below the chord of the WSS curve on linear
  axes.  For points with no cluster structure WSS falls like `c / k`, and the
  furthest point below that chord is exactly `sqrt(a * b)` for a ladder from
  `a` to `b`, so the answer was set by the ladder's ends: 1,500 uniform points
  gave 4, 4, 7, 9 and 13 for `max_levels` of 12, 20, 40, 80 and 160.  At the
  default `max_levels = 12` it also missed well-separated clusters (four
  clusters came back as 3).  The elbow is now read on log-log axes, where
  `c / k` is a straight line, and counts only when the curve sags clearly
  below it.  Two, three, four, five, eight and ten well-separated clusters are
  now recovered at the default on every seed tried (six came back as five on
  one seed in five).  With
  no elbow the function returns `NA` (see Breaking changes), and under
  `on_no_elbow = "warn"` its linear-axis answer, each with a warning that
  the ladder chose that answer; `build_tessellation()` and
  `get_voronoi_seeds()` refuse a geometry-only `resolution_profile()` with
  no elbow instead of drawing a count from it.  A ladder of two levels (`max_levels` of 1 or 2, or three
  points) has no line to test.  The warning now says the ladder is too short
  to read an elbow from and names what ended it.  It used to say the curve
  fell in a straight line "as it does for points with no cluster structure",
  which two clusters 90 m apart at `max_levels = 2` were told.

* **One invalid polygon changed the assignment rule for the whole layer.**
  `assign_features_to_polygons()` wrapped `st_join(largest = TRUE)` in a
  retry without `largest` on any error.  The comment blamed a predicate that
  cannot take `largest`, but sf never calls the predicate on that path; the
  retry fired when GEOS threw on an invalid ring, and every straddling
  feature was then assigned by `tie_break` instead: 95 of 200 buffered
  parcels changed cell, and a polygon 91 percent inside one cell went to its
  neighbour, with no warning.  Invalid features and cells are now repaired
  with `sf::st_make_valid()` for the join only, with a warning counting them,
  and a join that still fails stops and names `largest = FALSE`.

* **Lon/lat polygons were assigned to projected cells on bent cell edges.**
  `assign_features_to_polygons()` moved the cells into the features' CRS, so
  lon/lat features pulled the package's own projected cells into lon/lat,
  and with s2 the largest-overlap join failed and fell into the silent
  retry above: 55 of North Carolina's 100 counties went to a cell other than
  their largest overlap on a 36-cell grid.  The join now runs in the cells'
  CRS whenever it is projected, and the features come back with the
  geometry they arrived with.  Lon/lat points near a cell edge can change
  cell as a result; they now agree with a join done in the projected CRS.

* **`build_tessellation(method = "triangles")` dropped most points at UTM
  coordinates.**  qhull lifts each point onto x^2 + y^2, and at projected
  magnitudes (a northing near 5e6) the lift had no precision left to
  separate points a few metres apart, so they never became vertices: 200
  points over 100 m gave 26 triangles instead of 386, and the help page's own
  example lost 8 of its 20 points.  The points are now centred before
  triangulation.  Triangles that were already right are the same triangles,
  but qhull returns them in a different order, so triangle `cell_id` values
  change.

* **Data around a pole were projected to Web Mercator.**  A layer spanning
  more than 180 degrees of longitude with no gap was treated as global
  coverage.  Antarctic stations came out with worst-case distance errors near
  20,000 percent (the South Pole at y = -2.4e8 m), and Voronoi cells put 9.5
  percent of Arctic locations in a station's cell that was not their
  nearest.  A layer that lies wholly on one side of the equator now gets a
  Lambert azimuthal equal-area projection centred on its pole whenever that
  measures a smaller distance error than the global fallback: about 2
  percent on the same stations.

* **GWR mixed elevation into its distances.**  `prep_model_data()` kept the
  Z (and M) coordinate of POINT Z input, which GPS layers and
  `st_as_sf(coords = c("x", "y", "z"))` produce.  `predict.gwr_fit()` then
  handed GWmodel three coordinate columns, which it reshaped into two,
  scrambling the prediction locations with no warning;
  `gwr_model_selection()` ranked models on 3-D distances; and
  `fit_gwr_model()` and `cv_gwr()` failed.  Z and M are now dropped in
  `prep_model_data()` and again where a GWR takes its coordinates.

* **`predict.gwr_fit()` lost every prediction to one bad location.**  It
  went through `GWmodel::gwr.predict()`, which returns nothing for any row if
  one location's window is empty or singular, so fixed-bandwidth block CV
  failed every fold.  It also never assigns its distance matrix once
  training and new rows together exceed 10,000 (a 100 by 100 prediction grid
  came back all `NA`), and it built the full training hat matrix for a
  variance it then discarded, which is cubic in the training size.
  Predictions now come from the package's own local regressions, one per
  new location, and only a location with no estimable prediction is `NA`,
  with a warning counting them.

* **On brms 2.17 to 2.22, one far-off row moved every Bayesian prediction in
  the call.**  brms rebuilt the Hilbert-space GP's boundary from the rows
  being predicted, and the package's padding rows could only widen it, so a
  single row outside the training envelope changed the basis for every row,
  interior ones included: a 40 by 40 grid padded 15 percent moved the
  in-bbox cells by a mean of 22 percent of the surface's standard deviation,
  and `predict_surface()` depended on `chunk_size`.  On those versions the
  GP term's boundary factor is now rescaled per call so the boundary stays
  at its fitted value, and a row beyond the boundary, where the basis means
  nothing, is `NA` with a warning.  brms 2.23 stores the boundary itself;
  there only the `NA` rule applies.  The help page no longer says the
  boundary "has to grow".

* **A `bayesian_fit`'s cached fitted values could come from another model.**
  The cache lives in an environment, so it is shared by every copy of a fit,
  and the entry was stamped with the row count and a digest of the training
  data only.  Two fits over the same data hash identically however different
  their engines are, so `refit <- fit; refit$engine <- <re-estimated>` made
  `fitted()`, `residuals()`, `model_metrics()` and `summary()` on *either*
  object return the other's numbers --- and `clear_fitted_cache()`, which the
  help page offers for exactly this case, could not fix it, because clearing
  through one copy cleared the one shared entry and the next call re-wrote it.
  The entry now carries the environment rstan and brms create for each
  sampling run (`@.MISC` on the stanfit), which tells engines apart as well
  as the engine itself does and is written once, and is used only for the
  engine it came from (`identical()`, which settles the common case by
  pointer, so nothing is slower and no memory is held that the fit did not
  already hold, in a session or on disk).
  A stale entry is no longer deleted on a miss either: the fit that wrote it
  still wants it.  `new_spatial_fit()` now always builds a fresh cache rather
  than adopting one that arrived in `info`, and `summary()` no longer carries
  `.cache` at all --- a summary was not a value snapshot (its contents changed
  when anyone later called `fitted()` on the fit), `clear_fitted_cache()` on a
  summary emptied the *fit's* cache, and `saveRDS()` on one serialised an
  environment holding the full n-vector.

* **Every CRPS was `NA` above 46,340 posterior draws.**  `.crps_energy()`
  formed its weights with `m * m`, where `m` is `nrow(draws)` and therefore an
  integer, so the product overflowed and took `CRPS`, `mean_CRPS` and the
  whole `predictive_coverage` summary to `NA` behind one
  `"NAs produced by integer overflow"` warning.  48,000 draws is an ordinary
  `cv_bayes(fit_args = list(chains = 4, iter = 13000))` run.  The arithmetic
  is now done in double precision; the values it produces were, and remain,
  exact against the closed-form Gaussian CRPS.

* **Fold numbers depended on the machine's collation.**  Character fold labels
  are turned into a factor, whose levels `as.factor()` sorts under
  `LC_COLLATE`, and a fold's number is its level's position.  Labels differing
  only in case or punctuation --- `"north"` and `"North"` --- therefore landed
  in a different order under `C` than under `en_US`, so `fold_metrics$fold`,
  `predictions$fold` and `fold_status$fold` named different groups on
  different machines from the same data and the same seed.  The partition was
  never affected, so pooled scores were right.  Character levels are now
  sorted with `method = "radix"`, which is always C collation, making the
  numbering a property of the labels alone.  The C order applies to
  character labels only: numeric labels are numbered in numeric order, as in
  2.0.0, and a factor by its own level order (an earlier development build
  sorted numbers as strings too, so with ten or more numeric labels fold 2
  was the user's label 10).  `area_of_applicability()` given a vector of
  fold labels now numbers them by the same rule; it sorted character labels
  with `as.factor()` under the session's collation, so a message could name
  a different fold from the one `cv_*()` named for the same label.  Its
  partition and threshold were never affected.

* **`residual_morans_i(k = )` silently answered a different question.**  `k`
  reached the weight builder unvalidated, where `min()` collapses a vector:
  `k = c(4, 8)` built the `k = 4` matrix and returned a statistic for
  neighbours the caller never asked for, with no condition raised, and
  `k = NA` aborted on `"missing value where TRUE/FALSE needed"`.  `k` is now
  checked.  Two neighbouring gaps are closed with it: a user-supplied
  `weights` matrix holding any `NA`, `NaN` or `Inf` is refused by name
  instead of aborting inside a guard after a log line reading
  `"row sums range from NA to NA"`, and a row whose geometry is empty is now
  dropped with a count, as `make_folds()` and `estimate_sac_range()` already
  do --- its residual is perfectly finite, so it used to survive into
  `FNN::get.knn()` and abort with `"Data include NAs"`.

* **A scalar argument that reaches `as.integer()` is validated.**  Seven
  exported functions took a count or a distance and passed it straight into
  `as.integer()` or an `if ()` test, where `NA`, `Inf`, a length-2 vector or
  any value above `.Machine$integer.max` aborts with
  `"missing value where TRUE/FALSE needed"` or
  `"'length = 2' in coercion to 'logical(1)'"` --- errors that name nothing
  the caller passed.  `fit_gwr_model(bandwidth = )` under `adaptive = TRUE`,
  `predict_surface(chunk_size = )` (`Inf` being the natural way to ask for one
  chunk), `voronoi_seeds_kmeans(k = )`, `voronoi_seeds_random(k = )` and the
  count resolved by `build_tessellation()` / `get_voronoi_seeds()` now report
  the argument and the bound.  `voronoi_seeds_kmeans(k = 0)` and a negative
  `k` used to return *one* seed in silence, which "at most `k`" does not
  describe.

* **`expand` was ignored rather than refused.**  `clip_target_for(expand =
  c(0.05, 0.05))` returned a clip target byte-identical to `expand = 0` with
  no condition raised, and `create_voronoi_polygons()` did the same for a
  non-numeric `expand`, because both tested `is.numeric(expand) && expand > 0`
  and short-circuited to FALSE.  A malformed `expand` is now an error naming
  the argument; the same values on a degenerate bounding box used to abort
  inside the expansion helper instead.

* **`summarize_by_cell(deff = "variogram")` aborted when given no value
  column.**  With neither `response_var` nor `predictor_vars` the internal
  primary column is `NULL` by design, and `df[[NULL]]` raised
  `"attempt to select less than one element in get1index"`.  The variogram
  design effect is a function of the cell's coordinates and the fitted
  correlation, not of any column's values, so every point in the cell now
  counts and the call returns its counts and `cell_weight` as documented.

* **`compare_models()` aborted on a list holding no `spatial_fit`.**
  `evaluate_insample()` warns and skips a non-fit, and returned `NULL` when
  every element was skipped (it is now an error, below); the `NULL` then
  became a bare list and `seq_len(nrow(NULL))` raised
  `"argument must be coercible to non-negative integer"`.  It now says which
  argument is wrong and what belongs there.

* **`fit_rf_model(include_coords = TRUE)`'s caveat said "once per session" and
  was not.**  `.log_warn_once()` records the key in a package-level
  environment, and under `cv_rf(parallel = )` that write happens inside a
  forked worker and dies with it: the paragraph printed once per worker, the
  parent's registry stayed empty, and the next sequential fit printed it
  again.  `cv_rf()` now raises it in the parent before dispatching any fold,
  so each worker inherits the already-warned flag and stays quiet.  Measured
  on two workers: three occurrences before, one after.

* **The grid cache's order registry outlived the environments it described.**
  Insertion order is kept in a package-level environment keyed by the cache
  environment's printed address, and `clear_grid_cache()` removes only the
  entry for the environment it is handed.  Passing a fresh `cache_env` per
  call therefore left one permanent character vector per call, for
  environments that had since been garbage-collected and could no longer be
  named.  A finalizer now removes an entry with its environment.

* `plot(fit, type = "variogram")` no longer runs its subtitle off the edge of
  the figure.  ggplot2 clips a label that is wider than the plot instead of
  wrapping it, and the sentence saying why no range was identified is up to
  150 characters, so on a six-inch figure it was cut mid-word.  Labels built
  from a fit's own numbers are now wrapped at draw time.

* `ensure_stable_poly_id()` could not give IDs to a tessellation this package
  had just built.  It repaired the geometry in the layer's own CRS and then
  transformed it to the sort CRS, but validity is a property of the geometry
  in the CRS it is measured in: two vertices a centimetre apart in a projected
  CRS can land on one longitude and latitude, and s2 calls the ring
  degenerate.  On a clipped hex tessellation of North Carolina, 2 of 18 cells
  that are valid projected are invalid once transformed, and `st_centroid()`
  on one of them aborted the call with "Loop 0 is not valid: Edge 1 is
  degenerate".  The sort copy is now repaired after the transform as well.  The
  geometry returned is still the caller's own, and the same cell gets the same
  ID whether the layer arrives projected, in lon/lat or in Web Mercator,
  except where two fine cells' centres lie within the sort key's rounding
  step of the same longitude (36 of 2,500 100 m cells straddling a UTM
  central meridian changed ID via EPSG:3035); join such layers on geometry.
* `select_features_forward()` now says when `fit_fn` is ignoring the
  variables it is handed.  The learner it takes is a function of
  `(train_sf, predictor_vars)`, but the one `cv_spatial()` takes is a function
  of `train_sf` alone, and a learner written for that --- `function(train_sf,
  ...)` --- swallows the second argument and fits the same model every time.
  Nothing errored, because the training layer still carries every column: each
  candidate scored exactly what the intercept-only model scored, no candidate
  improved on it, and the result was an empty `selected` and an `NA` `score`
  with no word about why.  Two different predictor sets do not produce the
  same cross-validated metric to the last digit, so that pattern is now
  recognised at the first step and reported as a warning naming the fix.  The
  result is still returned.
* `create_grid_polygons()` warns when both `cellsize` and `target_cells` are
  supplied.  It already did for `cellsize` and `n`, and the documentation says
  to supply exactly one of the three, but `target_cells` was dropped in
  silence when `cellsize` was present.  Through `build_tessellation()` that
  meant `method = "hex", approx_n_cells = 25, cellsize = 10` returned
  however many cells a 10-unit lattice holds and said nothing about the 25.
  `cellsize` still wins; the override is now logged like its sibling.
* `build_tessellation()` warns when `approx_n_cells` or `cellsize` is supplied
  with `method = "voronoi"` or `"triangles"`.  Both arguments size the hex and
  square lattices and nothing else --- Voronoi grows one cell per input point
  and Delaunay one triangle per neighbouring triple --- and both used to be
  dropped in silence.  So `build_tessellation(pts, method = "voronoi",
  approx_n_cells = 25)` returned one cell per observation: the degenerate
  nearest-neighbour case, where every cell holds a single point, there is no
  within-cell variation and every standard error is `NA`.  The warning names
  the argument and points at `get_voronoi_seeds()`, which is where a Voronoi
  cell count is actually set.  It adds that `params` does not record the
  request either, so a saved result carries no sign of it: the Voronoi
  branch returns `create_voronoi_polygons()`'s own list, which has no slot
  for the argument, and the triangles branch no longer echoes
  `approx_n_cells` back (below).  The warning fires under `quiet = TRUE`,
  which gates this function's `message()`s and is documented not to silence
  R warnings.

* `determine_optimal_levels()` fits each k as the best of 25 k-means++
  restarts (Arthur and Vassilvitskii 2007; Fränti and Sieranoja 2019;
  Steinley 2003) instead of `stats::kmeans(nstart = 5)`.  The WSS curve is
  read for its shape, and with a handful of random restarts it carried
  optimisation noise: on eight-cluster layouts a sweep over k = 1..30
  rose at one or two steps in three of five draws, and an earlier form of
  the elbow rule once selected such a bump.  With the new budget the same
  sweeps rose at no step.  A curve that still rises is now logged as a
  warning naming the number of rising steps, and the model-aware
  diagnostics carry it as `wss_bumps` beside `wss_spread` (the relative
  spread of WSS across restarts at each k) and `nstart`.  A selection
  made on a curve that had a bump can differ from before --- those were the
  cases that were wrong; a clean curve gives the same answer.  Under a
  model-aware criterion the function also now warns *before* the sweep when
  `max_levels` leaves no k above the nine-cell floor, rather than
  fitting every k first and falling back afterwards.  The seeding draws each
  centre by inverting the cumulative squared distance rather than with
  `sample.int(prob = )`: the same law, and a `resolution_profile()` of 3000
  points takes 10 s instead of 35 s (`determine_optimal_levels()` on 3000
  points, 4 s instead of 13 s), but it gives different centres for the same
  seed than earlier development builds did.

* `make_folds(method = "block_kfold")` can now raise its "block dimension <
  autocorrelation range" warning.  The comparison was always there, but the
  range it compared against was estimated only under `auto_range = TRUE`,
  the one setting in which the blocks had already been sized from that range
  and the warning could never fire; on every default call the diagnostic was
  dead code.  With `auto_range = FALSE` (the default) and a `response_var`
  to hand --- always the case when a `cv_*()` function builds the folds ---
  a range is now estimated for the diagnostic alone.  It sizes nothing: the
  blocks are the same geometric blocks as before and the folds do not
  change.  A hand-set `block_size` below the range raises the same warning.
  The estimate's own log lines stay off the console, and the check is
  skipped (with an INFO log line saying so) when `gstat` is not installed or
  there are fewer than 30 points.  Only a grid dimension that is split is
  compared with the range, because a single row (or column) of blocks
  borders no other block across its width.  On points along a line the
  comparison was with 0, so the warning fired on every call that had a
  response: a 10 km line had 15 blocks 667 m long against a 499 m range.
  Its advice, `block_size = 499`, made the blocks shorter (19 of 523 m).  A
  10 km x 100 m corridor was compared with its 100 m width in the same way.
  A `response_var` that names no column of `points_sf` is now an error for
  `block_kfold`, whether or not `auto_range` is set.  With `auto_range` off,
  a misspelt name used to switch the check off silently.

* `fit_rf_model(include_coords = TRUE)` logs its caution once per session
  rather than once per fit.  Inside a five-fold `cv_rf()` or a twenty-fit
  `cv_block_size_sweep()` the same paragraph printed on every fit, which
  reads as twenty problems rather than one decision; the message now says
  it will not repeat.

* `fit_rf_model()` reports what ranger actually objected to.  ranger diagnoses
  a bad argument in its C++ layer, writes the diagnosis straight to stderr and
  then throws "User interrupt or internal error." --- so `mtry = 99` on a
  two-predictor forest printed "mtry can not be larger than number of
  variables in data. Ranger will EXIT now." to the console and raised an error
  naming neither the argument nor the problem.  That line is not an R
  condition, so `suppressMessages()`, `withCallingHandlers()` and `tryCatch()`
  all missed it: it escaped every handler to the console (and into CI logs,
  where it reads as an error from a test that is passing), while `cv_*()`
  recorded the placeholder as the fold's cause.  The message stream is now
  diverted for the duration of the call, so the diagnosis becomes the reported
  reason --- `fold_status$message` included --- and nothing is printed behind
  the caller's back.  Output from a call that succeeds is passed through
  unchanged, and when the stream is already diverted (under testthat, knitr or
  `capture.output(type = "message")`, where only one sink is permitted) the
  call runs exactly as before.  So does a call made after the session temp
  directory has been deleted, when there is nowhere to divert the stream to;
  it used to fail with "cannot open the connection".

* **One empty point made `determine_optimal_levels()` return 1.**  The
  function reduced every feature to a point and projected it but never
  dropped an empty or non-finite one, so a single `POINT EMPTY` made the
  `k = 1` WSS `NA`, k-means failed at every `k`, and the failure handler
  shrank the sweep to nothing: two well-separated clusters that gave
  `2 1 3` gave `1` once one empty row was added, with only a log line about
  interpolating the WSS.  Such rows are now dropped with a warning that
  gives their number, and under `select_on = "split"` the returned
  positions still index the layer as passed.

* **A few missing predictor values took Moran's I away from
  `determine_optimal_levels()`.**  A cell mean over a row with one missing
  predictor was `NA` and dropped the whole cell, and how many cells that
  removed depended on the points per cell, so `z` was computed on a
  different subset of cells at each `k`.  Three missing values in 400 rows
  made every `z` in the evaluated window `NA`, and the call fell back to the
  geometric ranking with a warning that did not mention missing values.
  Rows with a missing or non-finite response or predictor now stay in the
  WSS sweep and the cells and are left out of Moran's I, with a logged
  count, so every cell mean uses the same rows.

* **A misspelt `response_var` or `predictor_vars` in
  `determine_optimal_levels()` silently changed the criterion.**  A column
  that was not there counted as "no model variables": with the default
  criterion a typo kept `"geometric"` where supplying both variables then
  upgraded to `"combined"` (on one 400-point layer `7 6 8` instead of
  `11 10 7`, with no warning and no diagnostics), and under `"morans_i"` or
  `"combined"` the one log line said the variables were required although
  both had been supplied.  A named column that is not in the layer is now
  an error, as it is in `resolution_profile()`, and a model-aware criterion
  given no `predictor_vars` (or no `response_var`) falls back with an R
  warning that says which is missing, rather than a log line.

* **`determine_optimal_levels()` stopped one level short when locations
  repeat.**  The sweep was capped one short of the number of distinct
  locations, the bound `stats::kmeans()` needs only when no location
  repeats: five stations visited thirty times each could not reach `k = 5`.
  The cap is now the number of distinct locations, and still one short of
  the number of points.  A `k` whose WSS is 0 (to within 1e-12 of the total:
  a cell on every location) is left out of the log-log elbow line, and when
  the rest of the curve has no elbow, the fall to zero is the elbow.
  Leaving the level out and reading the rest made the call warn that the
  five stations had no cluster structure and return `3 2 4` (one metre of
  jitter gave `5 4 6`); it now returns `5 4`.  Two stations visited thirty
  times each give `2 1`, not `1 2`.  Zero is relative because k-means leaves
  floating-point residue: two groups of ten stations visited ten times each
  have a WSS of 7.8e-17 at `k = 20`, which dragged the whole line down and,
  at `max_levels = 30`, reported no cluster structure (still answering 2).
  `resolution_profile()` reads its `elbow` column the same way.  The
  no-elbow warning names the bound that ended the ladder (the distinct
  locations, the points or `max_levels`); it named `max_levels` whichever
  bound it was.

* **The same layer with its rows in another order gave
  `determine_optimal_levels()` another answer.**  The subsample and every
  k-means start index rows, so a permutation of the input moved the WSS
  curve and could move the count.  The rows are now put in coordinate order
  (response and predictors breaking ties) before either, so any
  permutation of a layer gives the same result.  Results for a given seed
  differ from 2.0.0's for that reason too.

* **Under `criterion = "morans_i"`, `top_n` could change the count
  `determine_optimal_levels()` ranks first.**  The model-aware criteria
  re-run k-means at each count of the window around the elbow, the elbow
  plus or minus `max(4, top_n)`, and the re-runs drew from one random
  stream in order.  A `top_n` above 4 starts the window lower, so more
  re-runs came before each count, and the cells of every count in the
  window changed, and with them `wss_eval`, `moran_i` and `moran_z`.  On
  480 points in twelve clusters with `max_levels = 30`, `"morans_i"` put 15
  first at `top_n` 1 to 4 and 16 at `top_n = 5`, both counts in both
  windows, and k = 17 had a |z| of 0.0013 at `top_n = 6` and of 0.30 at
  `top_n = 7`.  Under `"combined"` the ranking stayed and the diagnostics
  moved.  Each re-run is now seeded with one integer drawn before the sweep
  plus its count, and the stream is put back after it, so for one
  `set_seed` a count gets the same cells, `wss_eval` and Moran's I in any
  window that holds it.  A
  wider window still scores more counts, and under `"morans_i"` one of
  them can rank first: on the same layer the first count is 16 at `top_n`
  1 to 4, 17 at 5 and 6 and 19 at 7, each the count of smallest |z| in its
  window, and 17 and 19 are counts the wider windows add.  For a given
  seed, `moran_i`, `moran_z`, `wss_eval` and the `"morans_i"` ranking can
  differ from before, and a `"morans_i"` call can fall back to the
  geometric ranking where it ranked, or rank where it fell back; the
  sweep, the elbow and the `"combined"` ranking are as before.  With
  `set_seed = NULL` the re-runs no longer draw from the caller's stream,
  so the stream after the call no longer depends on `top_n`.

* **Moran's I of the cell means was `NA` for a response in small units, and
  in `determine_optimal_levels()` moved with the response's origin.**
  Moran's I and its z do not change under a change of unit or origin, but
  the statistic was refused whenever the residuals of the cell means had a
  sum of squares under `.Machine$double.eps`, whatever their unit.  A
  response of standard deviation 0.5 in units of 1e-9 got `moran_i` and
  `moran_z` `NA` at every level of `resolution_profile()`, with no warning
  that names them, while `cp` was the same as in units.  Under
  `determine_optimal_levels(criterion = "morans_i")` every count was `NA`,
  and the call fell back to the geometric ranking with no diagnostics and a
  log line that blamed too few cells with a response or a singular
  regression.  `determine_optimal_levels()` also formed the cell means from
  the raw response and regressed them, so at an offset of 1e15 they carried
  the rounding of the offset: on 2000 points in 40 clusters, `"morans_i"`
  ranked 43 35 38 at that offset and 36 37 35 at 0, from the same
  spread.  The residuals are now refused when they are within the rounding of
  the cell means they come from, a bound relative to those means, and a
  response whose raw values lose their spread is taken about its median
  first, as `resolution_profile()` already takes it.  A response whose
  largest absolute value lies outside 2^-255 to 2^255 is divided by a power
  of two first, which leaves its I and z as they are; in units of 1e200
  Moran's I came out `NaN`, and in units of 1e-200 `NA`, so
  `determine_optimal_levels()` fell back there too.  The layer above gets the
  same counts and the same `moran_z` in units, in units of 1e-9 and at 1e15,
  and the profile the same `moran_i` and `moran_z` in units and in units of
  1e-9 and 1e-12.  A constant response, one the predictors fit exactly and
  cell means that differ by rounding alone still give `NA`, now in any unit:
  under `determine_optimal_levels(criterion = "morans_i")` a response of
  1e15 in every row, and one of 1e10 times (3 + 2 w) with `w` the predictor,
  were ranked by Moran's z of their rounding, and now fall back to the
  geometric ranking.  Both fallback log lines name cell means that differ by
  rounding alone among their reasons.  A response whose raw values show their
  spread, through a predictor (1e15 plus 1e5 times the predictor plus a
  spread of 1) or at an offset such as 1e11 plus a spread of 1, can still
  leave residuals small beside the cell means, and the regression of the raw
  cell means then carried the rounding of the offset: with a predictor the
  call fell back at 1e15, and Moran's I moved by a relative 4e-4 to 2e-3
  at 1e13.  Where
  the residuals' sum of squares is under 1e-12 of the cell means', or they
  show no spread, the cell means are now formed and regressed again about
  the response's centre, as `resolution_profile()` takes its fit, and are
  refused only if those residuals are within the rounding of the regression
  or of the raw values.  Where a response's residuals are larger than that
  and its values lie between 2^-255 and 2^255, it gets the same Moran's I as
  before, unless its cell means differ by rounding alone.

* **`determine_optimal_levels()` blamed the wrong cause when its
  model-aware criteria had nothing to score, and its help page understated
  the fix.**  The model-aware pass scores only the elbow's neighbourhood,
  and on points with no cluster structure the elbow sits near
  `sqrt(max_levels)`, so every candidate stayed at or below the nine-cell
  floor until `max_levels` was about 40 (measured on 1000 uniform points:
  12, 20 and 30 all fell back, and 40 scored `k` = 10 and 11 alone).  The
  help page said "above roughly 10", and the fallback was logged as "Moran's
  I could not be computed".  The logged warning now names the window and
  the floor and points to `resolution_profile()`, and the help page gives
  the measured numbers.  Which candidates are scored is unchanged.

* **The lines `determine_optimal_levels()` logs where the window around
  the elbow ends at nine cells or fewer told a `"combined"` call that the
  response would take part with a larger `max_levels`, and named
  `max_levels` whatever bound ended the ladder.**  The line logged before
  the sweep (`k_max` of 9 or fewer) and the one logged after it said that
  the model-aware criteria carry no information there and that the call
  falls back to the geometric ranking for that reason, and the first said
  to raise `max_levels` "for the criteria to contribute".  `"combined"`
  ranks by the sag of the WSS curve and reads no response at any count: it
  returns the geometric ranking there by a fixed rule, and a larger
  `max_levels` gives the sag's order.  For `"combined"` the two lines now
  say that it returns the geometric ranking when the window ends at nine
  cells or fewer and that Moran's I, which it reports beside its ranking,
  is computed from ten cells.  The line before the sweep also began
  "max_levels leaves k_max = ..." and advised raising `max_levels` when
  the distinct locations or the points ended the ladder: on 200 rows at 8
  locations it said so at `max_levels` 12 and 30 alike, and no
  `max_levels` takes `k_max` above 8 there.  It now names the bound that
  ended the ladder ("the 8 distinct locations leave k_max = 8", "one short
  of the 10 points leaves k_max = 9") and advises raising `max_levels` only
  where a larger `max_levels` takes `k_max` above 9; otherwise it says that
  no `max_levels` takes `k_max` higher on these points.  Where
  `max_levels` equals the other bound (8 locations at `max_levels = 8`)
  the line names `max_levels` and no longer advises raising it.  On a
  subsample (more rows than `sample_n`) the points and the distinct
  locations that end the ladder are the subsample's, and the line, the
  warning for a curve with no elbow and the warning for a ladder too short
  to read named them as the layer's ("one short of the 9 points" on 400
  rows at `sample_n = 9`): they now say "one short of the sample_n = 9
  points drawn from 400" or "the 8 distinct locations in the sample_n = 30
  rows drawn from 361", and the advice is to raise `sample_n` where that
  takes `k_max` higher.  The warning for a ladder too short to read
  advised raising `max_levels` to 3 wherever `max_levels` named the bound,
  as on 3 distinct points at `max_levels = 2`, where no `max_levels`
  lengthens the ladder; it now says to choose the count on other grounds
  there.  The `"morans_i"` lines bound by `max_levels` alone, the counts
  and the diagnostics are as before.

* **`build_tessellation()` could not tessellate CRS-less planar points
  inside a CRS-less boundary.**  `ensure_projected()` marks such points
  `crs_assumed = "none"` ("planar, leave alone"), and `build_tessellation()`
  read that mark as a CRS name: `st_crs("none")` failed with "invalid crs:
  none" for every method, including the documented
  `boundary = clip_target_for(pts)`, so CRS-less planar data could not be
  gridded at all.  Only a real assumption (EPSG:4326) is now given to the
  boundary, which is refused if its coordinates cannot be degrees (below);
  with no assumption both stay in the same unnamed space.

* **When only one of the points and the boundary had a CRS, the
  tessellation builders stopped on sf's bare "st_crs(x) == st_crs(y) is not
  TRUE".**  UTM points read from a CSV with a UTM boundary failed in all four
  methods of `build_tessellation()` and in `create_voronoi_polygons()`;
  projected points with a boundary that had lost its `.prj` failed in
  `create_voronoi_polygons()` and `method = "triangles"`, and
  `clip_target_for()` returned a target with no CRS --- for lon/lat points
  with a CRS-less lon/lat boundary, in degrees, with `expand = 20` buffering
  by 20 degrees.  The side without a CRS is now interpreted in the other's,
  as `harmonize_crs()` does, with a warning: lon/lat-looking coordinates are
  reprojected from EPSG:4326, others are stamped.  CRS-less points that do
  not look like lon/lat are refused, with a message saying what to do, when
  the boundary is geographic, because stamping degrees on them would be
  wrong; a CRS-less boundary beside geographic points is read as lon/lat or
  refused (below).

* **A `boundary` without a CRS got a log line in `make_folds()`, `cv_*()`
  and `predict_surface()`, where every other function raises an R warning,
  and `cv_rf()` warned twice about one lon/lat boundary.**  With projected
  points, `build_tessellation()`, `clip_target_for()`, `plot_folds()` and
  the rest said "`boundary` has no CRS ... stamping" as an R warning, while
  `make_folds()`, `cv_*()` and `predict_surface()` stamped it with a "WARN
  ensure_projected(): input has no CRS" log line that `tryCatch()` and
  knitr never see and `spatialkit_quiet()` hides, naming neither the
  function nor the argument.  A boundary whose coordinates looked like
  lon/lat got two R warnings from one `cv_rf()` call, one from
  `prep_model_data()` and one from `make_folds()`, both naming
  `ensure_projected()`.  All three now warn as the others do, naming
  themselves and `boundary` (`make_folds()` also `prediction_points`,
  `predict_surface()` also `grid` and `covariates`), and a `cv_*()` call
  warns once, naming the `cv_*()` function.  Which CRS the layer ends up in
  is unchanged.  For `predict_surface()` this matters most on `grid` and
  `covariates`: a layer stamped with the wrong CRS puts every covariate
  lookup in the wrong place.  The stamping warning of every function
  now names the CRS it stamps ("stamping the target CRS ('EPSG:32632')
  WITHOUT reprojection") instead of "the supplied `crs`", an argument most
  of them do not have.

* **A geographic `crs` made every tessellation method work in degrees.**
  `build_tessellation()`, `create_voronoi_polygons()` and
  `create_grid_polygons()` took `crs = 4326` as the CRS to compute in, so
  Voronoi cells stopped being a nearest-point partition (150 points at
  53-57N: 20 percent of sampled locations lay in another point's cell), grid
  cells were neither square nor equal-area, and clipping them under s2
  stopped with "Edge 0 is degenerate" (hex) or left a point inside the
  boundary with an `NA` index (square).  A geographic `crs` is now the CRS
  the result is returned in: the cells are built and indexed in the local
  projected CRS `ensure_projected()` picks, then transformed with long edges
  densified.  A hex or square grid sized by an explicit `cellsize` is still
  laid in degrees, since that is the unit `cellsize` is in.

* **`ensure_projected()` stopped on lon/lat polygons that s2 rejects, and
  chose its UTM zone differently when `sf_use_s2()` was off.**  The centre
  that places the zone was `st_centroid(st_union())`: with s2 on, a polygon
  with a repeated vertex (valid for GEOS, common in shapefiles) stopped
  `ensure_projected()`, `create_grid_polygons()` and
  `prep_model_data(boundary =)` with "Edge 1 is degenerate (duplicate
  vertex)"; with s2 off it was planar in degrees, so two clusters at 0N and
  60N near -78 got UTM zone 17 in one session and zone 18 in another, and
  sf's warning and message about it got past `quiet = TRUE`.  The centre is
  now always taken on the sphere, a geometry s2 rejects is repaired first,
  and a mean of unit vectors is the last resort.

* **A single study-area polygon was never scored, so `ensure_projected()`
  kept a UTM zone at any extent.**  A polygon layer is scored on one point
  per feature, and one point makes no pair: every candidate scored `NA` and
  the selector fell back to the zone.  A CONUS outline stayed in UTM zone 15
  (12.8 percent worst-case distance error) where a Lambert azimuthal scores
  2.1 percent, and `prep_model_data(boundary =)` moved the whole analysis
  into the zone with it.  Layers with fewer than 40 features are now scored
  on their outline's vertices too.

* **`create_grid_polygons()` laid a grid over a near-global lon/lat boundary
  in Web Mercator.**  The cells were equal on the map and not on the ground:
  the true areas of whole cells differed nearly five-fold.  When the CRS
  picked for distances distorts areas across the boundary by more than 1
  percent, the grid is now laid in the equal-area CRS
  `ensure_projected(purpose = "area")` picks (Equal Earth here), with a
  logged warning; a local extent keeps its UTM zone.
  `create_grid_polygons_cached()` makes the same choice.

* **With `sf_use_s2(FALSE)` the package needed lwgeom, which it does not
  depend on.**  The stable-ID sort key is measured in lon/lat, so every
  Voronoi tessellation --- projected ones included, and the examples of
  `create_voronoi_polygons()` and `build_tessellation()` --- failed with
  "package lwgeom required", as did `ensure_stable_poly_id()`,
  `create_grid_polygons_cached()` on a cache miss, and random and k-means
  seeding on a lon/lat boundary.  These measurements now run on the sphere
  (s2) whatever the session's setting, which is restored afterwards; the
  IDs are the ones s2-on sessions always got.

* **Hex cell size depended on which way the boundary lay.**
  `st_make_grid()` builds hexagons from `cellsize[1]` alone, and that was the
  box's width over a rounded column count: a 1 x 1000 strip at
  `target_cells = 9` got 1734 hexagons where the same strip lying flat got
  89.  The size is now counted along the longer side, which leaves every
  boundary at least as wide as it is tall with exactly the grid it had.

* **Voronoi with `expand > 0` returned the boundary before it was grown.**
  The cells are clipped to the grown boundary, so they covered 1.93 km^2
  against a returned `boundary` of 1 km^2, and points up to `expand` outside
  it were indexed.  `create_voronoi_polygons()` and
  `build_tessellation(method = "voronoi")` now return the grown boundary,
  and the documentation says `expand` grows the study area too.

* **Collinear points gave an empty Delaunay triangulation.**
  `build_tessellation(method = "triangles")` on a transect returned no cells
  and an index of `NA`s, after logging that `delaunayn()` had failed, which
  it had not.  It now stops with a message that says the points are
  collinear and points to `method = "voronoi"`; the fallback's log line
  names the reason that applies.

* **k-means seeding clustered CRS-less lon/lat points as if degrees were
  metres.**  `voronoi_seeds_kmeans()` and `get_voronoi_seeds(method =
  "kmeans")` projected only when a CRS said lon/lat, so 800 CRS-less points
  89 km wide and 111 km tall were split east-west where the same points
  tagged EPSG:4326 were split north-south.  They now apply the lon/lat
  heuristic `ensure_projected()` applies (with its warning) and return the
  seeds in the input's own coordinates.

* **`voronoi_seeds_random()` returned the same seeding on every call.**  Its
  default `set_seed = 456` reset the random-number stream inside the call,
  so five calls under `set.seed(1)` to `set.seed(5)` gave one draw, and the
  sensitivity comparison its help page recommends compared a seeding with
  itself.  The default is now `NULL`, as in `get_voronoi_seeds()`: the draw
  comes from the session's stream.  Pass `set_seed` for a fixed seeding.
  This changes the default.

* **`harmonize_crs()` refused a layer as `target_crs`.**  A multi-row sf
  failed with "the condition has length > 1" and a one-row one with "cannot
  create a crs from an object of class sf", where `ensure_projected()`
  accepts both.  An sf or sfc target now means its CRS.

* **`summarize_by_cell()` failed on `deff = NA` and mis-recorded
  `deff = Inf`.**  The check on a numeric `deff` came out `NA` for
  `NA_real_` and `NaN`, so the call stopped with "missing value where
  TRUE/FALSE needed" instead of the documented warning and fallback to 1.
  `Inf` passed the check: the standard errors were the uncorrected ones,
  yet `cell_weight` was 0 in every cell and `deff_applied` recorded
  `deff = Inf`.  Any `deff` that is not a single finite number of at least 1
  now falls back to 1 with the warning.

* **`summarize_by_cell(deff = "variogram")` fell back to uncorrected
  standard errors with no R warning when no model could be fitted.**  With
  predictors only and no `sac`, without gstat, or when
  `estimate_sac_range()` returned no fit (fewer than 30 points), the only
  signal was a log line: no `tryCatch()`, `withCallingHandlers()` or
  `warnings()` saw it, and `spatialkit_quiet()` hid it, while the standard
  errors came out about 5 times smaller than the corrected ones in one
  check.  It is now a warning that names the reason.  Every design-effect
  fallback (a refused `deff`, a rejected or unsupported variogram, no
  model) is raised with class `"spatialkit_deff_fallback"`, so a loop over
  many summaries can catch exactly that case.  A rejected `sac` that a
  variogram estimated from `response_var` replaces is not a fallback: it
  gets a plain warning, and the classed warning is raised only when nothing
  replaces it, once per call, naming every reason.  A pure-nugget model (no
  structured component) implies that distinct observations are uncorrelated,
  so it is applied as a design effect of 1 in every cell, with
  `deff_applied = TRUE`; it was reported as "the supplied model could not be
  read", with the fallback warning.  (The fallback is now to `"kish"`; see
  Breaking changes.)

* **An empty point switched off `summarize_by_cell(deff = "variogram")` for
  its cell.**  Its missing coordinates made the cell's mean correlation
  `NA`, and the cell silently got a design effect of 1: the standard error
  of that cell was a quarter of the corrected one in one check.  (In the
  development version the same input stopped the call with "missing value
  where TRUE/FALSE needed".)  Points with empty or non-finite coordinates
  now count towards their cell's values but not towards its correlation,
  with a warning.

* **`summarize_by_cell(deff = "variogram")` read an anisotropic variogram
  as isotropic.**  Only `model`, `psill` and `range` were read, so
  `vgm(0.8, "Exp", 300, 0.2, anis = c(0, 0.2))` gave a correlation of 0.677
  at 50 m east-west where gstat's is 0.348, and a median design effect of
  14.1 against 7.7: standard errors too wide by a factor of 1.7.  A 2-D
  geometric anisotropy (`ang1`, `anis1`) is now applied as gstat applies it.
  (`resolution_profile()`, which uses the same correlation function on
  distances alone, still reads the major range.)

* **A misspelt or non-numeric `response_var` was dropped in silence.**
  `summarize_by_cell()` reported it through a progress message, which the
  default `quiet = TRUE` suppresses, and returned a frame with no `resp_*`
  column and `cell_weight` equal to `n`.  A missing or non-numeric
  `response_var` or predictor is now a warning, and a `response_var` of
  more than one name is an error instead of "the condition has length > 1".

* **A double cell ID of 100000 lost its cell in `summarize_by_cell(cells_sf
  = )`.**  When the points' and the cells' ID columns had different classes
  (an integer `poly_id` against a double from a GeoPackage Integer64 field
  or a CSV), both were converted with `as.character()`, which writes
  `1e+05` for the double and `100000` for the integer.  That cell came back
  with `NA` summaries and its points left the result: 11 of 30 points in one
  check, with no R warning.  Whole numbers are now written out in full, and
  a summarised ID that matches no cell is reported with a warning.

* **`summarize_by_cell(cells_sf = )` did not read the ID columns
  `assign_features_to_polygons()` writes from.**  Cells keyed by `id` or
  `grid_id` (common in shapefiles) were assigned cleanly and then summarised
  to a plain table with no geometry, `cell_area` or `n_per_area`, and a
  layer with `id` and a differently numbered `cell_id` was joined on
  `cell_id`, putting 89 percent of the summaries on the wrong polygons in
  one check.  The cells are now searched in the order the assignment used,
  and a `cells_sf` that cannot be joined is a warning (an error with
  `area = TRUE`) rather than a log line.  `agg_funs = median` or `"median"`
  is honoured instead of being replaced by the mean.  A single function is
  named after the expression passed, so `stats::median` gives
  `resp_median_*` as `median` does (it gave `resp_agg1_*`).

* **`assign_features_to_polygons(largest = TRUE)` assigned polygon features
  that only touch the cells.**  sf keeps the largest intersection piece
  without checking its area, so under GEOS a feature sharing only an edge or
  a corner with the cell layer went to that cell with zero overlap, while
  under s2 (lon/lat) the same feature was unassigned: sf's nc counties
  against 50 of them as cells gave 70 rows projected and 50 in lon/lat.  A
  feature with no overlap area is now unassigned in every CRS.

* **`tie_break = "smallest_area"` depended on the row order after all.**
  Candidates of equal area --- the cells of any regular grid, for a point on
  a shared edge --- fell through to the first row, so reversing the cells'
  rows moved every edge point to the neighbouring cell (points at x = 100
  went to cells 1, 4 and 7, or to 2, 5 and 8), and on a cached grid, whose
  IDs do not run row by row, one cell could take both of its edges.  Equal
  areas (to 9 digits) are now decided by the lowest, then leftmost,
  bounding-box centre.  On a `create_grid_polygons()` square grid that is the
  cell row order already picked (1,681 of 1,681 lattice points unchanged); on
  a hex grid 3 of 56 shared vertices move.

* **`create_grid_polygons_cached()` could return another site's grid.**  The
  cache key used the CRS's `input` name, which is `"unknown"` for any custom
  CRS read back from a GeoPackage or shapefile, so two site-centred CRSs
  with the same local boundary coordinates shared an entry, and the second
  site received the first one's grid 11,000 km away.  The key now hashes the
  CRS's WKT; the only cost is a rebuild when one CRS arrives written two
  ways.  `target_cells` now defaults to `NULL`, as in
  `create_grid_polygons()`, so `cellsize =` or `n =` work without it.

* **`ensure_stable_poly_id()` could number cells differently with s2 off.**
  The sort key's centroid was taken planar in degrees when
  `sf::sf_use_s2(FALSE)`, which moves it by far more than the key's rounding
  step, so near-tied cells swapped IDs between s2-on and s2-off sessions (4
  of 2,000 Voronoi cells), and without lwgeom the s2-off call stopped at the
  area.  The key is now taken on the sphere for the sort copy only, whatever
  the session setting.

* **`summarize_by_cell(deff = "variogram")` built each cell's correlation
  matrix up to five times per column.**  The same call with two partly
  missing columns and `conf_level` now builds 15 where it built 45.

* **`estimate_sac_range()` returned the widest direction that reached a sill
  when the all-pairs variogram had run past the fitted lags.**  Under an
  unremoved trend the pooled variogram rises without a sill, which the help
  page said the fitted-lag bound catches; but when two of the four
  directional variograms (those across the slope) did reach one, their
  maximum came back as the range with `anisotropy_used = TRUE` and nothing
  on the console.  On an exponential field of range 150 with an east-west
  trend, 12 of 30 draws returned 131--596 this way and 16 others `NA`, so
  `make_folds(auto_range = TRUE)` switched between a 2 x 2 grid of 430 m and
  geometric blocks from one draw to the next.  The directions that reach a
  sill are the shorter ones, so their maximum is a lower bound, not an
  estimate.  A converged all-pairs fit past the fitted lags is now refused
  whatever the directions found (`rejected_reason = "fitted range exceeds
  the largest lag fitted"`, the directional ranges still attached), and
  the directional maximum does not stand in for any other all-pairs fit
  either (see Breaking changes).  On stationary fields whose range is close to the
  cutoff this also turns a few draws from a directional maximum into `NA`
  (3 of 30 at an effective range of 570 on a 1000 m square).

* **`print()` on an `estimate_sac_range()` result showed a bare number.**
  For lon/lat input the number is in metres of a CRS the estimate picked,
  and with `predictor_vars` it is the range of the residuals rather than of
  the response; neither showed, so a mismatch with the layer it was about
  to be used on could not be seen.  A last line now names the unit and the
  CRS (`in metres of EPSG:32617`) and whether the variogram is of the
  response or of its residuals (and by which `detrend` method).  For a layer
  with no CRS the line says the range is in that layer's own coordinate
  units (`in the coordinate units of a layer with no CRS`), still with what
  was modelled; it used to be left out.

* **`estimate_sac_range(predictor_vars = )` fitted its variogram to the raw
  response, trend and all, when a single predictor or response value was
  infinite.**  `lm()`'s `na.exclude` drops `NA` but not `Inf`, so one `Inf`
  among 200 rows stopped the detrending ("NA/NaN/Inf in 'x'"), and the
  function fell back to the raw response with a warning.  The range came
  out at 2950 instead of 2168 (the answer with that row removed), and
  `make_folds(auto_range = TRUE, predictor_vars = )` built blocks 37
  percent larger.  `resolution_profile()`, which hides that warning, then
  warned that the variogram it had estimated itself was of the raw
  response.  With `detrend = "reml"`, a `-Inf` predictor (a `log(0)`
  covariate) first gave a false warning that the REML fit "did not
  converge", and then the OLS fallback failed the same way.  Rows with a
  missing or non-finite response or predictor are now left out of the
  detrending fit and the variogram, with a logged count, as
  `prep_model_data()` and `resolution_profile()` already do.  One `Inf` now
  gives the same range as that row set to `NA` (2167.5 detrended on the
  example; 1598.4 under REML).

* **`make_folds(drop_empty_blocks = FALSE)` could return folds with no test
  points.**  `k` was lowered only when the highest block id holding a point
  was below it, and with empty blocks kept that id says nothing about how
  many blocks hold points: two clusters on a 4 x 4 grid gave id 16, so
  `k = 5` was kept for 2 occupied blocks and three folds came back empty,
  with no warning (the imbalance check skipped an empty fold).  The `cv_*()`
  functions then ran on the two folds that had points and
  `area_of_applicability(folds =)` failed.  `k` is now lowered to the number
  of blocks that hold points, as `?make_folds` already said, the log line
  names the numbers, and every point in a single one of several blocks is
  the single-block error it always was with the default.

* **The automatic block grid gave an east-west corridor more than twice as
  many blocks as the same corridor running north-south.**  Only the row count
  was capped, so a layer more than about `block_multiplier * k` times as
  wide as it is tall got `round(sqrt(15 * w/h))` columns at `k = 5`: 39 x 1
  on a 10 km x 100 m corridor against 1 x 15 turned on its side, blocks less
  than half as long, and a scheme drifting towards random k-fold (1-NN CV
  RMSE 0.77 against 0.94 on the same values; lower in 16 of 20 seeds).
  Points on one horizontal line were treated as a square and got a 4 x 4
  grid that collapsed onto the line, lowering `k` from 5 to 4.  The column
  count is now capped at `block_multiplier * k` as the row count was, so
  both orientations and both lines get 15 blocks at `k = 5`.  Folds change
  only for extents more than about `block_multiplier * k + 1` times as wide
  as they are tall.

* **`make_folds()` ignored `block_nx` or `block_ny` given alone, and
  accepted invalid ones.**  Giving one dimension sent the call to the
  automatic grid without a word (`block_nx = 10` alone gave a 3 x 4 grid);
  0, a negative, `NA` or a vector failed inside sf or base R, and 2.7 was
  truncated to 2.  The dimension given is now used and the other derived
  from the extent's aspect ratio (roughly square blocks), and each must be
  a single whole number >= 1.

* **With a non-rectangular `boundary`, `make_folds(block_kfold)` kept
  zero-area slivers as blocks.**  Where the boundary only touches a grid
  cell at a corner or along an edge, the clipped cell is a POINT or a
  LINESTRING, and it was kept as a block: packed into a fold under
  `drop_empty_blocks = FALSE` (3 of 13 blocks under a triangular boundary),
  and able to catch a data point lying exactly on the boundary as a
  one-point block of its own.  Only the areal part of the grid is kept now.

* **`make_folds()` failed with R's own errors on a missing `k` or an
  invalid `buffer`.**  `k = NULL` (or no `k`) reached `if (k < 2)` and failed
  with "argument is of length zero"; a `buffer` of `NA`, `numeric(0)` or
  length 2 failed with "missing value where TRUE/FALSE needed" or a length
  error, and an `NA` is exactly what `estimate_sac_range()` returns when no
  range is identified.  Both are now refused by name (the leave-one-out
  methods still need no `k`), and so is a `units` object passed as
  `buffer`, `block_size` or `block_nx`/`block_ny`, which used to fail inside
  the units package without naming the argument.

* **`make_folds(method = "buffered_loo")` said nothing when the buffer
  excluded no neighbour.**  The buffer is in the units of the CRS the folds
  are built in, which for lon/lat input is metres, so a buffer in degrees
  (0.1) excluded nothing and the scheme was plain leave-one-out: 79 of 79
  training points in every fold of an 80-point layer, no condition raised.
  It now warns when no fold excludes any neighbour, and `?make_folds` says
  what unit `buffer` is in.

* **NNDM folds could be far more optimistic than their target with only a
  log-file line to say so.**  When `min_train` stops the matching --- samples
  clustered well inside the prediction domain, the layout NNDM is meant
  for --- the realised distances stay short: one cluster predicted onto a
  20 km grid kept a median of 1171 m against a target of 6704 m, 96 of 100
  folds held at the floor, while `?make_folds` said the result is never more
  optimistic than the target.  `make_folds()` now warns when the floor
  leaves more than one point's worth of excess below `phi`, records
  `params$n_at_min_train`, and the documentation states the guarantee only
  where neither `phi` nor `min_train` binds.  The warning names the
  statistic that fires it, the largest excess of the realised
  nearest-neighbour ECDF over the target at distances up to `phi`, and
  gives the two medians only as context: they summarise all the distances,
  and the realised median can sit above the target's while the short
  distances are over-represented.

* **`make_folds(auto_range = TRUE)` fell back to geometric blocks with only
  a log line.**  When no range was identified (an unremoved trend, a range
  past the fitted lags, fewer than 30 points, gstat missing), the blocks the
  caller asked to be sized from the data were not, and under knitr,
  `spatialkit_quiet` or `tryCatch()` nothing showed it.  This is now a
  warning that gives the rejection reason.  `estimate_sac_range()` can give
  up before it fits anything: fewer than 30 points, fewer than 30 finite
  values, a constant response (or residuals, when the predictors explain
  the response exactly), points with no extent, or gstat missing.  It used
  to return a bare `NA` then, with the reason only in a log line, and the
  warning said only "estimate_sac_range() returned NA".  That `NA` now
  carries a `rejected_reason` attribute saying which (it is still unclassed,
  with no other attribute).  The warning quotes it, as do
  `kriging_adequacy()`'s no-model error and `summarize_by_cell()`'s
  fallback warning.

* NNDM fold construction releases FNN's copy of the neighbour tables as soon
  as it has them, so a second `n` x `n/2` pair is no longer held through the
  sweep and the construction of the folds.  The peak inside `get.knn()` is
  unchanged.

* **`cv_bayes()` rounded coverage levels to a whole percent, so close levels
  overwrote each other.**  Coverage columns were named
  `sprintf("coverage_%.0f", 100 * level)`: `coverage_levels = c(0.5, 0.975,
  0.985, 0.995)` gave three columns for four levels, `coverage_98` holding
  the 0.985 value and the 0.975 value lost, with no condition raised.  Levels
  given as percentages, `c(50, 80, 95)`, made every fold throw away its
  CRPS, `n_draws` and coverage in silence.  Columns are now named at full
  precision (`coverage_97.5`; the default 50/80/95 names are unchanged), a
  level outside (0, 1) or given twice is an error that suggests dividing by
  100 where that fits, and the result carries `coverage_levels`, the nominal
  level of each column.

* **An error in a `fold_info_fn` threw away all of that fold's extras
  without a word.**  `cv_spatial()` caught it and dropped the whole list, so
  the columns were `NA` with `fold_status` `"ok"` and nothing logged; in
  `cv_bayes()` one failing quantile cost `gp_k`, `n_draws`, CRPS and every
  coverage column.  The failure is now logged and named in
  `fold_status$message`, and `cv_bayes()` computes coverage and CRPS in a
  step of their own, so `gp_k`, `n_draws` and `yhat_sd` survive it.  A
  `fold_info_fn` that returns the wrong shape (a vector where a value
  belongs, an unnamed or duplicated element, or a name `fold_metrics`
  already has, such as `RMSE`) is an error that says so; `RMSE = -5` used to
  overwrite the fold's real RMSE.

* **`cv_spatial()` failed after fitting every fold when `..per_row` came
  back from some folds only.**  The prediction rows were stacked with
  `rbind()`, which died on "numbers of columns of arguments do not match"
  once every fold had been fitted, so the work was lost and the message did
  not point at the cause.  A fold without `..per_row` (returned
  conditionally, of the wrong length, or from a `fold_info_fn` that threw) now
  gets `NA` in those columns, and one of the wrong length is logged.

* **Parallel cross-validation discarded folds that shared a core with a
  failure.**  `mclapply()` ran prescheduled, handing each core a chunk of
  folds: an error that escaped one fold was copied to every fold of its
  chunk, and a worker killed for lack of memory took all of its folds with
  it.  With four folds on two cores, a failure on fold 2 lost fold 4 as
  well, and `overall` pooled 40 of 80 rows.  Each fold now runs in its own
  worker, so a failure costs that fold only, and an error that stops a
  sequential run (a `metrics` or `fold_info_fn` return value of the wrong
  shape) stops a parallel one too, naming the fold.  A run in which nothing
  fails gives the same numbers as before, since the per-fold seeds are drawn
  before forking.

* **`model_metrics(newdata = )` measured R-squared against a different
  baseline from every `cv_*()` function.**  It took the total sum of squares
  about the new rows' own mean, where cross-validation takes it about the
  training mean, so the same predictions on a split across a trend scored
  R-squared -0.89 here and 0.35 from `cv_spatial()`.  `model_metrics()`,
  `evaluate_insample()` and `compare_models()` with `newdata` now use the
  training mean too, the out-of-sample convention, and the help pages say
  which baseline R-squared uses.  In-sample numbers do not change.

* **R-squared and MAPE depended on the units of the response.**  The
  thresholds below which a total sum of squares or a percentage-error
  denominator counted as zero were absolute, so a response with standard
  deviation below about 1.5e-8 got R-squared `NA` while RMSE and MAE were
  fine (and `select_features_forward(metric = "R2")` selected nothing), and
  one on a 1e-15 scale lost MAPE and SMAPE as well.  "Zero" is now 100
  machine epsilons of the data's own magnitude for every metric: a rescaled
  response gets the same R-squared and MAPE, a constant one still gets `NA`,
  and on a response spanning many orders of magnitude a row whose
  denominator is no larger than 100 epsilons times the largest (such as 1e-9
  against 1e6) no longer enters MAPE.

* **R-squared was `NA` for a response at a large offset.**
  `model_metrics()`, `summary()`, `evaluate_insample()`, `compare_models()`,
  the `cv_*()` functions and `select_features_forward(metric = "R2")` took
  R-squared only where the total sum of squares showed spread on the scale
  of the values, a test that reads their size: a response of 1e15 plus a
  spread of 2.3 failed it, and R-squared and adjusted R-squared were `NA`
  with RMSE and MAE computed.  On 200 rows R-squared was 0.804 at offset 0
  and 1e14 and `NA` at 1e15.  Where that test finds no spread, the
  response is now read again centred and scaled, as the pooled summaries
  read a constant column, and R-squared is computed where it varies; a
  response of one value up to rounding still gets `NA`.  Every response
  the first test passes gets the R-squared it got.

* **RMSE and R-squared were wrong for a response whose squares underflow
  or overflow.**  `model_metrics()`, `summary()`, `evaluate_insample()`,
  `compare_models()`, every fold and the pooled scores of the `cv_*()`
  functions and `compare_models_cv()`, and `select_features_forward()` took
  the sums of squared errors and deviations on the values as they stand.
  For 200 rows with R-squared 0.768, scaled by 2^-540 every squared error
  fell below the smallest double: RMSE was 0 and R-squared 1, and at 2^-560
  R-squared was `NA`.  Scaled by 2^510 or more, the sums passed the largest
  double: RMSE was `Inf` and R-squared `NA`.  One prediction of 1e200
  gave RMSE `Inf`.  Where a sum is not finite, or is below `n * 1e-300`
  with errors or deviations that are not all 0, that sum is now taken
  again on its errors or deviations divided by a power of two, which is
  exact, and RMSE is multiplied back: RMSE over the scale, R-squared and
  adjusted R-squared are those of the unscaled response, and RMSE is `Inf`
  only where it passes the largest double.  MAE, MAPE and SMAPE do not
  change, and where both sums lie between `n * 1e-300` and the largest
  double every score is as it was.

* **Out-of-sample R-squared was `NA` for training responses near the
  largest double on platforms whose long double is a double.**  Every
  out-of-sample R-squared is scored against the mean of the training
  responses: each fold of the `cv_*()` functions and their pooled score,
  `model_metrics()` with `newdata`, and `select_features_forward()`'s
  `score_holdout` and training-fold-mean baseline.  `mean()` sums in a
  long double, and where that is no wider than a double (as on macOS on
  Apple silicon) the sum
  overflows for values near the largest double: `mean(c(x, x, x / 2, -x))`
  is `-Inf` for `x = .Machine$double.xmax`, whose mean is `0.375 * x`.  On
  120 points scaled by 2^1020 (largest value 6.1e307) one fold's training
  mean was `-Inf`, that fold's R-squared and the pooled one `NA`, where
  every scale up to 2^1019 gives 0.7655072; `model_metrics(newdata =)` of a
  fit trained on `c(x, x, x / 2, -x)` gave `NA` for 0.9894806, and
  `select_features_forward()` dropped that fold's rows from its
  training-fold-mean baseline.  Where `mean()` is not finite the training
  mean is now taken on the values divided by a power of two, which is
  exact.  Where `mean()` is finite it is used as before, so every other
  score is the same to the last bit.

* **`select_features_forward()` warned that a working `fit_fn` "appears to
  be ignoring its second argument" on a response in small units.**  The
  check compared every candidate's score with the intercept-only score by
  `all.equal(tolerance = 1e-12)`, which compares absolute differences when
  the scores are below 1e-12.  On 60 points with the response in units of
  1e-20 (RMSE 1.2e-20) a linear learner selected the same predictor as in
  ordinary units and also drew the warning.  The scores are now compared
  relative to the intercept-only score, so the warning comes at every
  scale for a `fit_fn` that ignores its second argument and at none for
  one that does not.

* **`cv_block_size_sweep()`'s `fold_sd` and the standard error of
  `select_features_forward(rule = "se")` depended on the response's
  units.**  Both were `sd()` of fold scores or of their per-fold
  differences, whose squared deviations fall into the subnormal range near
  1e-160 and pass the largest double near 1e160.  For a response scaled by 2^-530, `fold_sd` over the scale
  was off by up to 5.7e-4 relative; scaled by 2^530 it was `Inf` for every
  block size.  At 2^530 the "se" rule read an infinite standard error,
  warned that fewer than two folds could be scored while saying 5 were,
  judged the step on `tol` alone and selected a predictor the rule rejects
  on the same data in ordinary units.  The standard deviation is now taken
  on the scores divided by a power of two wherever `sd()` is not finite or
  is below 1e-150, and where `sd()` gives 0 for values that differ (their
  deviations all below about 1e-162, as in `cv_bayes()`'s `yhat_sd` for
  draws that small).  Elsewhere it is `sd()` as before, to the last bit.

* **`prep_model_data()` left out the shift of the dropped rows' response
  for a response whose values are near 1e160 or near 1e-160.**  The shift is read in
  standard deviations of the response, and `sd()` was `Inf` for a response
  scaled by 2^530 and 0 for one scaled by 2^-560.  On 60 points with 6 rows
  missing a predictor, the warning at scale 1 says the dropped rows'
  response sits 1.18 SD from the kept rows', more than two standard errors;
  at both scales that sentence was missing and
  `attr(, "dropped")$response_std_diff` was `NA`.  The standard deviation is
  now read as in `cv_block_size_sweep()`'s `fold_sd`, so both scales give
  1.18 SD and the same sentence.  Where `sd()` is finite and at least
  1e-150 the reading is as before, to the last bit.

* **`compare_models()` set out-of-bag random-forest metrics beside in-sample
  ones without saying so.**  Without `newdata` an `rf_fit`'s fitted values
  are out-of-bag and every other backend's are in-sample, and the table could
  rank the models the wrong way round: GWR RMSE 0.77 in-sample against RF
  0.82 out-of-bag, where the forest's in-sample RMSE was 0.40.
  `evaluate_insample()` and `compare_models()` now carry a `metric_basis`
  column (`"in-sample"`, `"out-of-bag"` or `"newdata"`), `compare_models()`
  logs a note when a table mixes them, and both help pages say what the
  metrics are computed on.

* **`compare_models()` put LOOIC and AICc side by side for fits on different
  rows.**  Both are sums over the rows a model was fitted to, so a model
  that lost 20 rows to a predictor's missing values showed LOOIC 32.2
  against 54.8 for the model on all 70, and looked 22.6 better while it was
  worse on the rows they share.  A column whose models were fitted to
  different rows is now set to `NA`, with a warning naming each model's `n`.
  The same happens to fits of the same rows with different responses (a
  response and its log, say), since an information criterion compares models
  of one response only, and the warning now says so.  It used to say they
  were "fitted to different rows (raw: n = 80, logged: n = 80)" and to
  "Refit them on the same rows".

* **`compare_models()` read significantly negative residual autocorrelation
  as missed spatial structure.**  The caution fired on a two-sided p-value
  whatever the sign, so the alternating in-sample residuals of a GP or a
  small-bandwidth GWR (Moran's I -0.13, p = 0.02) were logged as "may not
  fully capture the spatial structure", the opposite diagnosis.  A negative
  z is now logged as what it usually means, a model tracking its data
  closely.

* **`compare_models_cv()` placed polygon rows in its shared blocks by a
  different point than every model is fitted at.**  The shared folds reduced
  polygons and lines to their point-on-surface whatever `pointize` said,
  while each backend fitted them at the `pointize` point: with
  `pointize = "centroid"`, 119 of 150 L-shaped parcels were in a different
  fold from a standalone `cv_gwr()` run.  The shared blocks now use
  `pointize`; with the default `"auto"` nothing changes.

* **Saved folds on lon/lat polygons were refused after `sf_use_s2()` was
  toggled.**  The provenance check located each probed row by its centroid,
  which on a geographic CRS is spherical with s2 on and planar with it off;
  the two differ by up to 5e-4 degrees on county polygons, 500 times the
  tolerance.  Folds built before `sf_use_s2(FALSE)` (a common workaround for
  invalid polygons), or saved and read in a session set the other way, were
  rejected by every `cv_*()` as "built from different data".  The probe now
  always takes the planar centroid; folds saved by an older version are
  checked the way they were made.

* **`residual_morans_i()` gave the wrong reason when it could not use a
  fit's residuals.**  An error from `residuals()` was thrown away, and it and
  a fit with no `residuals()` method (the `?new_spatial_fit` example has
  none) were both reported as "could not extract enough residuals (n < 4)",
  on a 100-row fit; a residual vector of the wrong length was reported as
  "coordinate extraction failed".  Each now has its own warning, quoting the
  error where there is one, except that a fit with no `residuals()` method
  is now scored on the response minus `fitted()` instead (below).

* **A GWR whose local regressions interpolate the data won on AICc.**
  GWmodel's AICc is defined only while the effective number of parameters,
  tr(S), is below n - 2; past that its penalty turns negative.  A small
  adaptive bandwidth reached it, and so did any bandwidth raised to the old
  floor, which for the bisquare and tricube kernels (they give the farthest
  neighbour in a window weight 0) fitted every window exactly.  On 100
  points `compare_models()` listed R2 = 1 and an AICc of -15033 against 292
  for the automatic bandwidth, and at n = 60 `gwr_model_selection()`
  selected the real predictor plus three noise variables with -66689.
  `fit_gwr_model()` now reports such an AICc as `NA` and
  `gwr_model_selection()` ranks such models last, each with a warning giving
  tr(S).  The adaptive floor is one neighbour higher for bisquare and
  tricube, and raising a supplied bandwidth to it is now a warning, not a log
  line.  `bandwidth = NULL` was not affected above 20 points.  The raised
  floor is enough unless several neighbours tie at the kernel's edge (a
  regular grid); the warning now says so.  Where an adaptive bandwidth is
  already every observation, the undefined-AICc warning suggests fewer
  predictors, more observations or a gaussian or exponential kernel instead
  of a larger bandwidth.

* **`fit_gwr_model()` never checked a one-predictor model for local
  collinearity.**  The check ran only with two or more numeric predictors,
  but every local design includes the intercept, and a predictor nearly
  constant inside a window is collinear with it.  A regional covariate
  nearly constant within each of four clusters gave local slopes from -97 to
  221 around a true 3 with no warning, while adding a noise predictor to the
  same data warned at every location.  One numeric predictor is now enough.

* **GWR said a singular window came back as `NaN` coefficients; it stopped
  the fit.**  GWmodel's matrix inverse throws on an exactly singular window,
  so a 0/1 indicator constant within clusters failed the whole fit with a
  bare "inv(): matrix is singular", while the help page and the collinearity
  warning promised masked `NaN` coefficients.  The package now solves the
  local regressions itself: a coefficient a window cannot estimate is `NA`
  and the fit goes on (see Breaking changes).  A zero-width adaptive kernel
  at co-located points still leaves a window with no weighted row, and the
  warning names the co-located points when they are the cause.

* **An adaptive GWR bandwidth above the number of points was capped in
  silence.**  `bandwidth = 1500`, meant as metres with `adaptive` left at
  `TRUE`, became a 200-neighbour, near-global fit on 200 points without a
  word; its local slopes varied less than half as much as the intended
  fixed-distance fit's.  `fit_gwr_model()` and `gwr_model_selection()` now
  warn, naming n and pointing to `adaptive = FALSE`.  `cv_gwr()` repeats the
  warning in each fold whose training set is smaller than the bandwidth.

* **Below 20 points, `bandwidth = NULL` did not fit at the bandwidth
  `bw.gwr()` chose.**  GWmodel searches adaptive bandwidths from 20
  neighbours up to n, so with fewer points its choice exceeded n and was
  capped at n, a different kernel with a worse AICc, without a word.  The
  package's own search starts at the smallest bandwidth the model allows,
  so a layer this small is searched like any other.

* **A GWR predictor named twice raised false collinearity warnings.**
  `fit_gwr_model(predictor_vars = c("a", "b", "a"))` fitted correctly, but
  its collinearity checks ran on the duplicated column and warned "exactly
  singular" and "100% of locations collinear", once in every `cv_gwr()`
  fold, and the doubled count raised the bandwidth floor.  Names are now
  collapsed on entry, as `gwr_model_selection()` already did.

* **A `gp_c` you set did not change the `gp_k` derived for it.**  The basis
  count was always sized for the boundary factor the package would have
  chosen, so a wider boundary got the same number of basis functions and
  could no longer resolve the lower length-scale bound it was sized for.  The
  advice under `gp_c` is to raise it for a long-range surface, which is
  exactly the case that coarsened the basis.  On 200 uniform points,
  `gp_c = 3` fitted with `gp_k = 23` where the rule gives 43, and `gp_c = 5`
  with 23 where the rule gives 70 (capped at 50); the cap warning could never
  fire on this path.  With `gp_k = NULL` the derived `gp_k` is now sized for
  the `gp_c` actually used, and a capped value is logged.  An explicit `gp_k`
  still passes through untouched.

* **`fit_bayesian_spatial_model()` could not fit `brms::categorical()` or a
  `mixture()` family.**  Those families give each distributional parameter
  its own GP under the same coefficient names, and the automatic length-scale
  prior kept only the names, so every coefficient got two identical rows and
  brms stopped with "Duplicated prior specifications are not allowed" before
  sampling.  A user's `lscale` prior restricted to one `dpar` failed the same
  way, because it was copied onto every category's coefficients.  The prior
  now carries each coefficient's `dpar`, `nlpar` and `resp`, and a global or
  `dpar`-level `lscale` prior is expanded only onto the coefficients it
  addresses and never over a coefficient-level one the user already gave.
  With `standardize_predictors = TRUE` they still failed, on the automatic
  `normal(0, 5)` slope prior, which carried no `dpar` and so matched no
  slope of either family (brms: "The following priors do not correspond to
  any model parameter: b ~ normal(0, 5)", a prior the user never wrote).
  That prior is now set on each distributional parameter's slopes, as the
  length-scale prior is; a family with one `mu` gets the same single row as
  before.

* **A two-level factor response under `brms::bernoulli()` fitted, and then
  nothing could score it.**  The response check refused a non-numeric
  response only under gaussian, and the gaussian refusal itself pointed at
  `bernoulli()`.  brms fits the factor, but `residuals()` came back all `NA`,
  `summary()`, `model_metrics()` and `compare_models()` stopped on "response
  is factor", and `cv_bayes()` ran a full fold of MCMC before aborting in the
  fold scoring (with `parallel = 2`, every fold ended as `worker_error`).  A
  factor or character response is now refused, before anything is compiled,
  under every family except `categorical()` and the ordinal ones (cumulative,
  sratio, cratio, acat), with a message saying to convert it to 0/1.  Numeric
  and logical 0/1 responses are unaffected.

* **`predict()` on an ordinal or categorical `bayesian_fit` returned all
  `NA` as a "posterior draw failed".**  brms returns `posterior_epred()` for
  those families as a draws x rows x categories array, which the method took
  for a failed draw: a real `cumulative()` fit returned `NA` for all five
  new rows, with only a log line, while `fitted()`, `summary()` and
  `model_metrics()` said merely that they got an array.  `predict()` under
  its default `type = "epred"` and `fitted()` now stop, saying the family has
  a probability per category and pointing at
  `type = "predict", draws = TRUE`, whose share of draws in each category
  estimates its probability for any rows, and at
  `brms::posterior_epred(<fit>$engine)` for the training rows
  (`posterior_epred(<fit>$engine, newdata = )`, which the message used to
  suggest, refuses new rows without the scaled coordinates the method
  builds).  The message now counts the caller's rows, where it counted the
  two GP-boundary rows as well ("150 x 7 x 3" for five rows).
  `type = "predict"` without `draws = TRUE` on a `brms::categorical()` fit,
  which returned the mean of unordered category indices (1.46, 1.97, ...),
  is now an error; for an ordinal family it is the expected category index,
  as documented.  A genuinely failed draw still returns `NA` as
  documented, and the log line now carries the cause.

* **`predict()` on an `rf_fit` turned every ranger error into an all-`NA`
  vector.**  `type = "quantiles"` on a forest grown without
  `quantreg = TRUE`, and `type = "se"` without `keep.inbag = TRUE`, which the
  help page says are rejected, returned ten `NA`s for ten rows with no R
  condition, and `model_metrics(newdata =, type = "se")` then reported
  `n = 0`.  A failure in ranger's predict method is now an error naming
  ranger's reason; the `cv_*()` fold loop records it as the fold's cause and
  `predict_surface()` stops naming the rows, as they already did for other
  backends.  A `newdata` with no complete row still returns all `NA` with a
  log line, as for the other backends, rather than reaching ranger as a
  zero-row frame; that had made a `predict_surface()` chunk outside the
  covariates' coverage abort the whole surface.

* **`check_convergence = FALSE` returned `convergence_ok = TRUE`.**  The flag
  started out `TRUE`, so a fit whose checks never ran (its max R-hat was 1.28)
  claimed to have passed them over an empty diagnostics list, and `print()`
  had nothing to caveat.  It is now `NA` when nothing was checked, and
  `print()` on the fit and on its `summary()` says "Convergence: NOT
  CHECKED"; `summary()`'s printout also repeats the "Convergence warnings
  present" flag, which it carried and never showed.  A failed PSIS-LOO is now
  logged with its cause instead of "LOO computation failed." alone, which had
  left `compare_models()` showing `LOOIC` `NA` with nothing saying why.

* **The convergence check raised dozens of "The ESS has been capped"
  warnings.**  `brms::neff_ratio()` runs posterior's ESS over every GP basis
  weight, and posterior warns once per well-mixed one: an `n = 80` fit raised
  44 R warnings, 41 of them this one.  R keeps only the first 50 warnings, so
  a warning that mattered and came later, loo's Pareto-k among them, could
  be dropped.  That one message is now muffled around the R-hat and ESS
  accessors; every other warning passes through, and the ratios are
  unchanged.

* **A saved `rf_fit` or `bayesian_fit` carried its engine twice.**  The
  model formula was built in the fitting function's frame and so captured
  it, and that frame holds the forest or the `brmsfit` itself; a formula
  serialises its environment, so `saveRDS()` wrote the engine a second time
  (1.62 MB for a 100-tree forest of 0.72 MB; about 80 MB for a 40 MB
  `brmsfit`, which brms's own copy of the formula doubled even in
  `saveRDS(fit$engine)`).  The formulas now carry the global environment, as
  a formula typed at the console does.

* **A forest with rows out of every tree's bag said nothing.**  ranger returns
  `NaN` as the out-of-bag prediction of a row every tree sampled, so with
  `num_trees = 5` 20 of 200 rows had `NaN` fitted values and `summary()`
  printed "n = 200" over an R-squared computed on 180.  `fit_rf_model()` now
  warns with the count, and `summary()` prints "(computed on 180 of 200
  rows ...)" when its metrics use fewer rows than the fit has.  `cv_rf()`
  does not use its fold forests' out-of-bag predictions, so it warns once
  per run with the number of fold forests affected, instead of once per fold
  (each of which told the user to score the forest with `cv_rf()`).

* **`area_of_applicability()` counted rows that differ on a dropped
  zero-variance predictor as inside the AOA.**  A predictor constant in the
  training data --- a land-cover dummy absent from the training region --- is
  dropped from the distance, so new rows taking another value there were
  judged on the other predictors alone: 37 of 40 urban rows came out inside,
  where a single urban training row would have kept the predictor and put 1
  inside.  Such rows now get `DI = Inf` (their scaled distance along that
  predictor is infinite), are counted outside, and a warning gives the count;
  `print()` says how many.

* **`area_of_applicability()` returned a threshold of 0 from duplicated
  training rows without saying why.**  Each training row's reference is its
  nearest other row, so exact duplicates in predictor space --- repeat visits
  to a site with static covariates, covariates from a raster coarser than the
  sampling --- have a training DI of 0; past about three quarters of the rows
  the threshold is 0 and only exact copies count as inside (30 sites visited
  four times: 0 of 200 new points inside, against 197 after deduplication).
  Without folds a row's reference is now its nearest row with different
  predictor values, so repeats no longer set the threshold; with folds a
  threshold of 0 is a classed warning that gives the remedy
  (leave-location-out folds, or deduplication).  See the entry on repeat
  visits above.

* **`area_of_applicability()` refused all-zero weights, which the advice
  `pmax(importance, 0)` produces whenever the model found no useful
  predictor.**  A one-predictor forest gave all-zero weights in 9 fits of 20
  when its predictor carried no signal, so the AOA was lost in the folds where
  it mattered most.  With one predictor the weight cannot change the index
  (it is scale-invariant) and zero is accepted silently; with several, all
  are weighted equally, as `weights = NULL` would, with a warning.

* **A fractional `chunk_size` in `area_of_applicability()` marked
  extrapolation as inside the AOA.**  On the dense path (`use_fnn = FALSE`,
  or FNN not installed) block starts became fractional and the rows between
  blocks kept an initial DI of 0: with `chunk_size = 2.5`, 5 of 25 far-out
  points were reported inside and 16 training DI of 0 moved the threshold.
  `chunk_size` is now validated by name and truncated to whole rows.

* **`area_of_applicability(model = fit, folds = folds)` stopped with "fold 1
  refers to rows outside 1:n" whenever `prep_model_data()` had dropped a
  row, for the same folds `cv_*()` accepted.**  One missing response in 200
  rows was enough: `make_folds()` numbers the rows of the layer it is given,
  the fit keeps only the 199 rows `prep_model_data()` returned, and the fold
  IDs were read as positions in those.  The documented workflow ("pass the
  same `make_folds()` result you passed to `cv_spatial()`") therefore failed
  on any layer with a missing or non-finite modelling value or an empty
  geometry.  The rows the fit's `"dropped"` record names are now taken out
  of the folds, as `cv_*()` take them out, with a log line giving the count;
  a label vector with one label per row of the layer fitted from loses those
  labels.  The threshold is the one you get by removing the rows from the
  folds by hand.  Folds built on the model's own training data (the
  `prep_model_data()` output) are still read as positions in it, and a fold
  ID naming a row the data never had is still an error.

* **`area_of_applicability()` applied folds built on other rows without a
  word.**  Fold splits are row IDs, and `cv_*()` compare the sample of row
  locations `make_folds()` records against the data, refusing folds built
  on another layer.  `area_of_applicability()` did not: a model fitted on
  the same rows in another order took the folds anyway and moved the
  threshold (0.2432 against 0.2404).  It now makes the same check on the
  training data and refuses such folds with the `cv_*()` message.  The
  check is skipped (logged) when the folds were built on polygons and the
  training data are the points a fit reduced them to, so a model fitted on
  polygons with folds built on those polygons keeps working.

* **`predict_surface()` filled a polygon grid with covariates from an
  arbitrary point inside each cell.**  `st_nearest_feature()` returns the
  first zero-distance match the spatial index yields, so a
  `create_grid_polygons()` grid got covariates that did not match the
  location predicted at (predictions off by up to 2.6 on a 0-30 response) and
  that changed with the row order of `covariates` (by up to 4.7), and the
  result was a polygon layer where the manual promises points.  The grid is
  now reduced to one representative point per cell first.

* **`predict_surface(..., draws = TRUE)` flattened the draw matrix into
  `.pred`.**  The argument was forwarded to `predict()`, and a backend that
  honours it returned an `n_draws x n` matrix that became `.pred` column by
  column (correlation with the right values: 0.03); with `se = TRUE` the
  duplicated argument was reported as "backend does not expose posterior
  draws".  `draws` is now refused by name, and a `predict()` returning the
  wrong number of values is an error.

* **`predict_surface()`'s automatic grid could lose a whole column or row.**
  When the extent is an exact multiple of the cell size, `floor()` of the
  ratio landed one short through rounding (0.3 / 0.1 gives 2 cells), leaving
  a cell-wide strip uncovered; at the default `n_cells` this hit 1197 of
  10000 random squares.  The cell count now has a relative tolerance, so an
  exact multiple gets exactly that many cells.

* **`predict_surface()` kept a reused grid's old `.pred_se`.**  Passing an
  earlier surface as `grid` left that model's `.pred_se` beside the new
  `.pred` unless this call replaced it, even after logging "returning
  predictions only".  It is now removed unless this call computes it.

* **A `logger` configuration made before loading the package could still
  abort its functions, and received its log lines.**  `logger` seeds a new
  namespace by copying *every* index of the user's global configuration, and
  2.0.0 pinned the formatter on index 1 only.  A user with two global indices
  set up before `library(spatialkit)` got `formatter_sprintf` or
  `formatter_glue` on the console echo, so a `%` or a `{` in a message
  (`"fold 2 skipped: object 'cov_{x' not found"`) aborted the function that
  logged it and the R warning the manual promises never arrived; a third
  global index kept the user's own appender and received spatialkit's WARN and
  INFO lines in the user's log file.  Both indices now have formatter, layout,
  appender and threshold pinned, every message is marked
  `logger::skip_formatter()`, and copied indices beyond the second are deleted
  (on `logger` 0.2.2, which cannot delete one, switched off).  The global
  configuration is still never touched.

* **Deleting the session temp directory made every function that logs fail
  until the package was reloaded.**  The trace file's path was fixed in
  `tempdir()` at load time, so after an OS cleaner or `unlink(tempdir())`
  every log call failed with "cannot open the connection", and a documented R
  warning (`ensure_projected()`'s CRS assumption, say) became that error;
  `tempdir(check = TRUE)`, R's own recovery, did not help.  The trace now
  resolves its path when a line is written, recreates the directory if it has
  gone, and drops a line it cannot write; no logging failure aborts the caller
  any more, so the warning always arrives.

* **Logged cautions were missing from knitted documents.**  The console echo
  wrote to stderr, which knitr does not capture, so an R Markdown, Quarto or
  pkgdown document showed the package's R warnings but none of its logged
  cautions (`compare_models()`'s significant residual autocorrelation, for
  one).  While knitr is running the line is now also sent as an R message, so
  it appears in the output and `message = FALSE` hides it; a line that is
  raised as a warning too is not repeated.  stderr gets exactly what it got
  before, and nothing changes outside knitr.

* **`plot_tessellation_map(labels = TRUE)` drew no labels on any layer this
  package builds.**  `label_col` defaulted to `"grid_id"`, a column no
  function produces: Voronoi and Delaunay cells carry `cell_id`, grids
  `poly_id` and `cell_id`, and `summarize_by_cell()` output `poly_id`, so the
  map came back unlabelled with only a log line to say why.  `label_col` now
  defaults to `NULL`, which takes the first of `grid_id`, `cell_id`,
  `poly_id`, `polygon_id` and `id` the layer has; `grid_id` stays first, so a
  layer that has one is labelled as before, and naming a column still works.

* **`plot_tessellation_map()` failed at print on a units, Date, POSIXct or
  difftime fill column.**  The fill scale was chosen with `is.numeric()`, so
  Date, POSIXct and difftime columns got a discrete scale ("Continuous value
  supplied to a discrete scale"), and an `st_area()` column (class units)
  passed the test and then broke the viridis scale's arithmetic.  The
  function returned normally and the error came only when the plot was drawn.
  Date and POSIXct now get the continuous scale on a date or time axis, and
  units and difftime columns are drawn as numbers with the unit in the legend
  title (`"area [m^2]"`) unless `legend_title` is given.

* **`plot_folds()` failed at print when its layers disagreed on having a
  CRS.**  A CRS-less boundary beside projected points, or the reverse,
  aborted inside `coord_sf()` with sf's "cannot transform sfc object with
  missing crs".  Since `plot_folds()` began drawing the block outlines it
  also failed on the very layer the folds were built from: `make_folds()`
  projects CRS-less lon/lat points to a UTM zone, and stamps CRS-less points
  with a boundary's CRS, so the stored blocks carry a CRS the points do not.
  A CRS-less layer is now brought into the points' CRS, or failing that the
  folds' own, or the first layer's that has one: reprojected when its
  coordinates look like lon/lat, and stamped with a warning otherwise, as
  `make_folds()` does.

* **`plot()` on a custom `spatial_fit` without a `residuals()` method
  stopped with "could not extract residuals".**  `?new_spatial_fit` calls
  that method optional, but `residuals.default()` returns `NULL` for a
  `spatial_fit`, so the residual map, the observed-against-predicted plot and
  the residual variogram all refused a backend that had the required
  `fitted()` method.  The residuals are now the response minus `fitted()`,
  which is what the built-in backends return; without a `fitted()` method the
  error names the method to define.

* **`citation("spatialkit")` gave the year it was called in, not the year of
  the release.**  DESCRIPTION has no `Date` field, so `inst/CITATION` fell
  back to `Sys.Date()`.  CRAN installs got the right year only by accident:
  `meta$Date` partially matched `Date/Publication`.  The year is now looked
  up the way `utils::citation()` does it, by exact field name: the CRAN
  publication date, then `Date`, then the date `R CMD build` packaged the
  source (which a GitHub install via remotes or pak has).  Only an install
  straight from a source directory records none of these, and only then does
  the current year appear.

* **The plots' size arguments did nothing on ggplot2 older than 3.4.0.**
  The line layers pass `linewidth =`, which ggplot2 3.4.0 introduced; older
  versions warn "Ignoring unknown parameters" and draw at the default width,
  so `outline_size`, `boundary_size` and the other size arguments were
  silently ignored.  Suggests now asks for `ggplot2 (>= 3.4.0)`.

* The test suite calls `local_mocked_bindings()`, which testthat added in
  3.1.7, but Suggests allowed 3.1.5.  On 3.1.5 or 3.1.6 every test that mocks
  a function failed with "could not find function".  Suggests now asks for
  `testthat (>= 3.1.7)`.

* **`determine_optimal_levels(criterion = "combined")` could put first a
  cell count that Moran's I never scored, chosen by the rule the elbow had
  stopped using.**  On eight separated clusters (800 points,
  `max_levels = 40`, four seeds) it returned `10 7 6`, `10 7 6`, `6 5 10`
  and `6 10 7`.  The geometric axis was still the chord on linear axes
  across the elbow's window, which ranked 6 or 7 above the elbow of 8.  The
  candidates below the nine-cell floor, which Moran's I cannot score, shared
  an average rank that shrank as more of them went unscored (6 of 9 for
  seven of them), although the help page said they ranked last.  The
  geometric axis is now the log-log sag the elbow is read from, and exact
  ties go to the `k` nearest the elbow.  And an elbow below ten cells, a
  count Moran's I cannot score, is no longer ranked against the counts it
  can: ranking it put ten, the smallest count Moran's I scores, first
  whatever the response did (a response of noise and one varying by
  cluster gave the same answer).  There `"combined"` returns the geometric
  ranking, logs why, and records it in the diagnostics (`criterion =
  "geometric"`, `fallback`).  The same layers now give `8 7 9`, `8 7 9`,
  `7 6 8` and `8 7 9`, the geometric answer.  Moran's z no longer takes
  part in this ranking at all, and supplying both `response_var` and
  `predictor_vars` no longer selects this criterion (see Breaking changes).

* **`build_tessellation(method = "hex")` or `"square"` laid its lattice over
  a near-global lon/lat boundary in Web Mercator.**  The points' CRS is
  chosen for distances, and handed on as the grid's CRS it skipped the area
  check `create_grid_polygons()` makes: on a boundary from 170W to 170E and
  60S to 70N, the full hexagons differed 5.75-fold in true area, where
  `create_grid_polygons()` on the same boundary used Equal Earth (0.7
  percent).  The lattice is now laid where `create_grid_polygons()` lays it:
  in the CRS picked for the points unless that CRS distorts areas across the
  boundary by more than 1 percent, and otherwise in the equal-area CRS
  `ensure_projected(purpose = "area")` picks for the boundary, with a logged
  warning, the points indexed in the same CRS.  This applies with no `crs`
  and with a geographic one; a local extent keeps its UTM zone.

* **A CRS-less study area given with lon/lat points could be read as a
  one-metre square.**  `build_tessellation()` resolved a boundary without a
  CRS against the UTM zone picked for the points, after projecting them, so
  a one-degree tile with integer corners (which the lon/lat heuristic
  declines) was stamped with that zone: one Voronoi cell, or 27 hexagons,
  and all 50 points indexed `NA`.  A boundary in British National Grid
  metres was stamped with the UTM zone too.  Such a boundary is now read in
  the points' own CRS when its coordinates fit the lon/lat envelope, with a
  warning, and refused with an error naming both layers when they do not;
  `create_voronoi_polygons()` and `clip_target_for()` read it the same way.

* **CRS-less lon/lat points with a CRS-less boundary in metres failed with
  "`boundary` must be polygonal".**  `build_tessellation()` stamped
  EPSG:4326 on the boundary without looking at its coordinates, so a UTM
  polygon was transformed to nothing and the error named its geometry type.
  It now stops with an error saying the two layers cannot be placed in one
  space.

* **`create_voronoi_polygons()` tessellated CRS-less lon/lat points in
  degrees, silently.**  It projected only when a CRS said lon/lat, so 60
  CRS-less points at 55N got cells in which 17.7 percent of sampled
  locations were not nearest to their cell's point, while
  `build_tessellation()` on the same points warned, took them as EPSG:4326
  and projected them.  It now applies the same lon/lat heuristic, with its
  warning, and returns what `build_tessellation()` returns (0.1 percent, at
  the cell edges).  `?build_tessellation` no longer says a CRS-less pair
  "stays in the same unnamed planar space" whatever its coordinates.

* **Delaunay triangles returned in a geographic `crs` did not contain their
  own points.**  `build_tessellation(method = "triangles", crs = 4326)` on
  150 points left 9 of them touching no returned triangle, so a spatial join
  on the result did not reproduce `index`.  The corners of the returned
  triangles are now put back on the input points after the round trip
  through the working projection, and every point lies in its indexed
  triangle.

* **The tessellation builders failed on an sf layer as `crs`.**
  `build_tessellation()`, `create_voronoi_polygons()` and
  `create_grid_polygons()` stopped with an error from sf ("the condition has
  length > 1", or "cannot create a crs from an object of class sf").  They
  now take the layer's CRS, as `ensure_projected(target_crs =)` and
  `harmonize_crs()` do.

* **Random and k-means seeding on a lon/lat boundary warned "install package
  lwgeom" on every call.**  sf raises "coordinate ranges not computed along
  great circles" for each lon/lat draw when lwgeom, which this package does
  not depend on, is absent: one R warning per `voronoi_seeds_random()` or
  `get_voronoi_seeds(method = "random")` call and two per k-means call.
  That warning is muffled, and the draw is unchanged.  The boundary's union
  is taken on the sphere as well, so with `sf_use_s2(FALSE)` sf no longer
  prints its planar `st_union()` message and the seeds are the ones an s2-on
  session gets.

* **A transect with sub-millimetre scatter got a sliver study area and a
  166,536-cell grid for `approx_n_cells = 25`.**  `clip_target_for()` called
  a bounding box degenerate only when its two ends were equal to rounding,
  so 30 points along 1000 m with a y scatter of 1e-6 got a 1000 x 9e-7
  rectangle, over which the square grid had 166,536 cells (16 seconds) and
  the hexagonal one 154,980 (28 seconds); a little thinner, and
  `create_grid_polygons()` stopped at `max_cells` telling the user to check
  the units of a `cellsize` they had not passed.  A box whose short side is
  below a millionth of the long side is now degenerate too (a buffer around
  the points: 34 squares or 45 hexagons for 25), and the `max_cells` error
  names the argument the size came from, `target_cells` (`approx_n_cells`),
  `n` or `cellsize`.

* **A whole `build_tessellation()` result passed as `boundary` failed with
  sf's `no applicable method for 'st_geometry' applied to an object of class
  "list"`.**  `build_tessellation()`, `create_voronoi_polygons()` and
  `clip_target_for()` now say that the object looks like a
  `build_tessellation()` result and to pass its `$boundary`, as
  `make_folds(blocks = )` already did for `$cells`;
  `create_grid_polygons()`, `create_grid_polygons_cached()` and
  `ensure_stable_poly_id()` add the same hint to their type errors.

* **`build_tessellation(method = "triangles")` recorded an `approx_n_cells`
  it had ignored.**  `params$approx_n_cells` held the ignored count where
  the help page says the count used is kept; it is now `NULL`, as for
  Voronoi, and the warning says so for both methods.

* **`summarize_by_cell(deff = "variogram")` said it was "Falling back to
  deff = 1" for a rejected `sac`, and then applied a design effect.**  A
  `sac` whose fit `estimate_sac_range()` had rejected was set aside with
  that warning, after which a variogram estimated from `response_var` was
  fitted and applied: in one check every row came back corrected, with a
  median design effect of 5.2, and in the development version the warning
  carried the fallback class, so `tryCatch(spatialkit_deff_fallback = )`
  threw the corrected result away.  When the estimate was rejected too, the
  one fallback raised two R warnings.  A `sac` that an estimate replaces now
  gets a plain warning saying what replaced it, and the fallback warning is
  raised once per call, only when the correction applied is not the one
  asked for (Kish, or none), naming every reason, the rejected `sac`
  included.

* **`summarize_by_cell(deff = "variogram")` ignored a `sac` that carried no
  variogram model without saying so.**  A plain number, or
  `units::set_units(1.5, "km")`, was passed over and the design effect came
  from a variogram estimated from `response_var`, with no condition; the
  caller could not tell that the value given had not been used.  Such a
  `sac` is now set aside with a warning, as a rejected `sac` is: a plain
  warning when a variogram is estimated instead, and the classed
  `spatialkit_deff_fallback` warning, naming it, when `deff` falls back
  to `"kish"`.

* **`summarize_by_cell(deff = "kish")` recorded no correction when only the
  predictor standard errors were corrected.**  The `"deff_applied"`
  attribute followed the primary variable's ICC alone, so with an
  unclustered response (ICC 0) and a clustered predictor (ICC 0.82) no
  attribute was attached, and in the development version every row said
  `deff_applied = FALSE`, while the predictor standard errors had been
  inflated elevenfold.  The attribute is now attached whenever either ICC is
  positive; its `deff` stays the primary variable's (all 1 in that case).

* **`summarize_by_cell()` corrected the response's standard errors with a
  residual variogram without a word.**  A `sac` from
  `estimate_sac_range(..., predictor_vars = )` describes what the predictors
  leave unexplained, a weaker correlation than the response's own: the
  response standard errors came out at 0.60 of those from the response
  variogram in one check, while `kriging_adequacy()` warned about the same
  object.  Such a `sac` is still used as given, as documented, but a warning
  now says the response standard errors are understated.

* **`assign_features_to_polygons()` dropped the features' own `id` column
  when the cells were keyed by `id`.**  The polygons' ID went through the
  spatial join under its own name, so a site `id` collided with the cells'
  `id` and was dropped, with a warning about a collision the result never
  had: it only gains `polygon_id_col`.  Only a column named `polygon_id_col`
  is replaced now.

* **`assign_features_to_polygons(largest = TRUE)` let the polygon row order
  decide an exact tie in overlap.**  sf keeps the first of equal largest
  overlaps, so a 40 m square split evenly across the edge of cells 1 and 2
  went to cell 1, or to cell 2 with the polygon rows reversed.  An exact tie
  (areas equal to 9 significant digits) is now decided by `tie_break` among
  the equally largest polygons, and counted in `attr(, "ties")`.

* **Every largest-overlap assignment of polygon features raised sf's
  "attribute variables are assumed to be spatially constant" warning.**
  `sf::st_join(largest = TRUE)` adds grouping columns of its own before
  intersecting, so the warning came with every ordinary call (the reporting
  vignette hid it with `warning = FALSE`), said nothing about the data, and
  could not be avoided with `st_agr()`.  That warning alone is now muffled.

* **`summarize_by_cell(deff = "variogram")` accepted a `deff_max_n` of
  less than 2.**  A value of 1 or 0 subsampled every cell to one point or
  none, so the mean correlation came out `NaN`, the design effect 1 and the
  standard errors uncorrected (4.1 times smaller than with the default in
  one check) with nothing to say so; `NA` stopped the call with "missing
  value where TRUE/FALSE needed".  It must now be a single number of at
  least 2, and anything else is an error that names it.

* **The rows with no cell ID got a `cell_weight` of 0.**  A layer assigned
  with `keep_unassigned = TRUE` is summarised with its unassigned rows as a
  group whose ID is `NA`, and `summarize_by_cell()` lost that group's count
  of non-missing values: `cell_weight` was 0 beside `n = 5` and a finite
  standard error.  The group is now counted like any other.

* **`attr(, "deff_applied")$deff` turned into a vector of `NA`s for a fixed
  `deff` and one populated cell.**  With `cells_sf`, the realignment to the
  joined rows took a fixed `deff = 2` for a per-cell vector whenever exactly
  one cell was summarised, and recorded `c(2, NA, NA, ...)`.  A fixed design
  effect is now recorded as the number.

* **`make_folds(method = "block_kfold")`'s refusal of a grid above 1,000,000
  blocks told the caller to check a `block_size` they never passed, and a
  `block_size` hundreds of orders of magnitude too small slipped past it.**
  With `block_nx = 2000, block_ny = 1000` (or `block_multiplier = 1e6`) the
  error read "Check that `block_size` (unset) is expressed in the data's CRS
  units".  With `block_size = 1e-200` the cell count overflowed to `Inf`,
  which the guard let through, and `st_make_grid()` failed with "result
  would be too long a vector".  The message now names what produced the
  grid: `block_size`, the range `auto_range` estimated,
  `block_nx`/`block_ny`, or `block_multiplier` x `k`.  A count too large to
  represent is refused like any other.

* **`make_folds()` accepted any `block_multiplier`, and a `units` object for
  `phi` or `min_train` failed with an error that named no argument.**
  `block_multiplier = NA` died inside `st_make_grid()` with "'length.out'
  must be a non-negative number".  `c(1, 3)` silently used 3, and
  `units::set_units(100, m)` was read as 100, giving a 32 x 16 grid aimed at
  500 blocks.  `phi = units::set_units(100, m)` failed inside the units
  package with 'both operands of the expression should be "units" objects'.
  `block_multiplier` must now be a single positive number, and `phi` and
  `min_train` refuse a `units` object by name, as `block_size`, `buffer` and
  `block_nx`/`block_ny` already did.

* **`residual_morans_i()` and `compare_models()` had no residual Moran's I
  for a custom fit without a `residuals()` method.**  `?new_spatial_fit`
  says `residuals()` is optional, and such a fit falls through to
  `residuals.default()` and gets `NULL`, so `residual_morans_i()` returned
  `NULL` with a warning and `compare_models()` reported all-`NA`
  `resid_morans_*` columns.  On 80 points with an east-west trend the
  predictor could not explain, that hid a residual Moran's I of 0.754 (z =
  15.5, p = 4e-54).  Such a fit is now scored on the observed response minus
  `fitted()`, which is what the built-in backends' residuals are and what
  `plot()` already used for it.  `NULL` is returned only when that cannot be
  formed either, and the warning quotes the reason.

* **Re-using a `cv_*()` result's `$folds` renumbered every fold after a
  dropped one.**  The result's `$folds` holds the splits that survived, each
  carrying the `fold_id` it was reported under, so five folds with fold 3
  dropped are labelled 1, 2, 4 and 5.  Handed to a second `cv_*()` call, to
  score another learner on the same splits, they were numbered by position
  as 1, 2, 3 and 4, so the second run's fold 3 was the first run's fold 4
  (RMSE 2.90 in both), and a join on `fold` with the first run or with
  `fold_separation()` paired different folds.  Splits that carry a distinct
  whole-number `fold_id` now keep it.  Splits without one are still numbered
  by position, as are splits whose ids repeat.

* **A `fold_info_fn` whose `..per_row` reused a `predictions` column name
  corrupted `predictions`.**  A `..per_row` column called `yhat` (or `y`,
  `fold`, `..row_id` or `y_train_mean`) was bound on beside the original,
  and `predictions` came back with two `yhat` columns.  In the development
  version, `dplyr::bind_rows()` then renamed both (`yhat...4`, `yhat...6`)
  with only a message, so `overall` found no `yhat` and was all `NA` with
  `n_pred = 0`, beside per-fold RMSEs of 2.4 to 3.4 and no warning.  Such a
  name is now an error naming the column, as a clashing scalar extra already
  was.  Duplicated or empty `..per_row` names are errors too, in a parallel
  run as in a sequential one.

* **A `fold_info_fn` that returned a named vector lost its values without a
  word.**  `c(slope = 1.9)` in place of `list(slope = 1.9)` (the shape
  `metrics` accepts) added no column, and `fold_status` read `"ok"` with an
  empty message on every fold.  A named vector is now taken as the list it
  stands for.  Any other return value that is not a list or `NULL` (a
  function, an environment) is now an error.

* **Folds built on a pointized copy of a polygon layer were refused as
  "built from different data".**  `make_folds(coerce_to_points(parcels))`
  applied to `parcels` itself, with the same rows and IDs, was refused on 80
  L-shaped parcels because "64 of 64 checked row IDs sit at a different
  location here".  The error blamed "folds from another dataset".  The
  provenance check compares each polygon's centroid with the point
  `coerce_to_points()` gave it, and under `"auto"` that is
  `st_point_on_surface()`, which differs for every non-convex feature.  The
  folds' `params$row_probe` now records whether they were built on POINT
  geometry.  When that differs from the data, the error says so and tells
  you to build the folds with `make_folds()` on the layer passed, which
  reduces polygons to points itself.  Folds from data that really differs
  get the old message.

* **`evaluate_insample()` returned `NULL` when no element of `fits` was a
  `spatial_fit`.**  Its help promises a data frame with one row per model,
  and the only sign that every element had been skipped was a log line per
  element, which `spatialkit_quiet()` hides and `tryCatch(warning = )` never
  sees.  It is now an error that names `fits`.

* **`cv_rf(parallel = )` printed its core-count message twice when asked for
  more workers than the machine has.**  `cv_rf(parallel = 16)` on a 4-core
  machine printed "cv parallel: 16 workers requested on a machine with 4
  cores; using 4." twice.  It resolved the count once for its thread policy
  and once more when the folds ran.  It now prints the message once.

* **`fit_gwr_model()` called local regressions unstable because of where a
  predictor's units start.**  A temperature field in kelvin beside a second
  predictor drew "global design (intercept + predictors) has scaled
  condition index 230 ... (collinearity risk)" and "100% of 30 sampled
  locations have a collinear local design ... Local regressions there are
  unstable", while the same field in degrees C drew no global warning and 6
  of 30; the local slopes of the two fits agree to 1e-10.  Both indices were
  computed on the uncentred design, where a predictor far from 0 against its
  spread (kelvin, a year, elevation in feet) is collinear with the
  intercept.  That makes the local intercept an extrapolation to 0, but a
  slope's precision does not depend on where the predictor's origin is.  The
  global index is now computed on the centred predictors, so it is 1 for a
  single predictor and a change of origin does not move it.  The local
  survey keeps Belsley's uncentred index with the intercept as
  `info$local_collinearity$cn` and adds `cn_slopes`: the predictors centred
  at their weighted mean in the window, each divided by its standard
  deviation over the study area.  A window's slopes count as collinear when
  `cn_slopes` is above 30 or singular.
  `n_local_collinear` and the warning count those windows.  The kelvin and
  degrees C fits now get the same verdict (0 of 200 locations), and a
  regional covariate nearly constant within clusters is still flagged
  everywhere (200 of 200).  A window where only `cn` is above 30 is logged,
  not warned about.

* **`gwr_model_selection()` checked an adaptive bandwidth less strictly than
  `fit_gwr_model()`.**  `bandwidth = 3e9` with `adaptive = TRUE`, a distance
  passed as a neighbour count, failed with "NAs introduced by coercion to
  integer range" and then a bare "missing value where TRUE/FALSE needed".
  `bandwidth = 0.5` was rounded to 0 neighbours and quietly raised to the
  floor, where `fit_gwr_model()` refuses it.  `gwr_model_selection()` now
  runs `fit_gwr_model()`'s check before it prepares anything: an adaptive
  count must be at least 1 and at most R's largest integer.  `cv_gwr()` runs
  the same check once, up front, instead of failing it in every fold and
  returning "all folds failed".

* **`gwr_model_selection()` with a fixed bandwidth in the wrong units failed
  with a bare "inv(): matrix is singular".**  On lon/lat data projected to
  EPSG:32617, `bandwidth = 0.2, adaptive = FALSE` is 0.2 metres against an
  extent of 44720 metres, so every window is empty.  `fit_gwr_model()`
  warned about exactly this and explained the singular window, but the sweep
  said nothing.  It now raises the same "less than a ten-thousandth of the
  data's extent" warning, and its error explains a singular window as
  `fit_gwr_model()`'s does.

* **`print()` on a GWR fit showed a fixed bandwidth in scientific notation
  and without its unit.**  A fixed bandwidth of 122372 m printed as
  "Bandwidth: 1.224e+05 (fixed, bisquare kernel)", so the value was rounded
  to four digits and its unit was missing, although the help page tells the
  reader to check it.  It now prints "122,372 metre (fixed, bisquare
  kernel)", and an adaptive one as "42 neighbours (adaptive, ...)".

* **`fit_gwr_model()` did not refuse a character or factor response, as the
  README says every model function does.**  A response read from a CSV as
  text went into `GWmodel::bw.gwr()`, which failed twice with "Not
  compatible with requested type".  That drew the arbitrary-fallback
  bandwidth warning, and the fit then stopped with "'x' must contain finite
  values only", naming neither the column nor the cause.  It now stops first
  with "response 'zc' is not numeric (it is character)", as `fit_rf_model()`
  and `fit_bayesian_spatial_model()` do.

* **`cv_bayes()` under an ordinal or categorical family sampled every fold
  and then scored none of them.**  Such a family predicts a probability per
  response category, and cross-validation scores one number per row, so each
  fold compiled and sampled a full model and was then discarded: under
  `brms::cumulative()`, `k = 2` on 50 rows took 4.35 minutes to return "all
  2 folds failed to produce predictions".  `cv_bayes()` now refuses
  `brms::categorical()` and the ordinal families (cumulative, sratio,
  cratio, acat) before fitting anything, with a message naming the family,
  and the "Which metrics survive a non-Gaussian response" section (in
  `?model_metrics`, `?cv_bayes` and `?compare_models_cv`) no longer calls
  CRPS and interval coverage meaningful for "any family the backend
  accepts".

* **`predict_surface()`'s automatic grid left up to a cell of the training
  extent uncovered on the east and north.**  It took
  `floor(extent / cell_size)` cells from the lower-left corner, so whenever
  the extent was not a whole number of cells the rest of it got no
  prediction.  `cell_size = 100` on a 980 x 956 extent covered 900 x 900:
  13.5 percent of the box was uncovered, and 14 of 120 training points lay
  in no cell.  `cell_size = 334` gave 2 x 2 cells and left 65 of the 120
  points out.  The grid now has enough cells to cover the box and is centred
  on it.  It overhangs the box by less than one cell, split evenly between
  the two sides, so every cell centre still lies inside the box.  An extent
  that is an exact multiple of the cell size gets the same grid as before.
  At the default `n_cells`, a grid usually gains one column or row (102 x 99
  cells instead of 101 x 98 on the extent above), and its centres move by
  less than half a cell.

* **`plot_tessellation_map(labels = TRUE)` warned at print on every lon/lat
  layer, and failed at print on a units label column.**  The label points
  were computed with `st_point_on_surface()`, and `geom_sf_text()` ran it
  again on those points when the plot was drawn.  On longitude/latitude
  cells that raised "st_point_on_surface may not give correct results for
  longitude/latitude data", although nothing was wrong.  A `label_col`
  holding `st_area()` values (class units) failed with "units package is not
  attached", as a units fill column did.  The label points are now used as
  computed, and a units or difftime label is drawn formatted, with its unit
  (`"15073.393 [m^2]"`).  A `fill_col` or `label_col` naming more than one
  column, or `NA`, is now refused by name.  It used to fail with R's
  "'length = 2' in coercion to 'logical(1)'".

* **`resolution_profile()` with `predictor_vars` on which no row is
  complete returned `rss`, `cp` and the Moran columns `NA` with base R's
  bare warning, and `summary()` and `select_resolution()` blamed a missing
  response.**  On 400 points that all carry a response and a predictor
  that is `NA` on every row, the OLS fit was taken on zero rows, which
  raised "number of rows of result is not a multiple of vector length",
  and the log said the raw response was scored, which no row was.
  `summary()` said "cp needs a response" and `select_resolution(criterion =
  "cp")` that `cp` "needs a response". The fit is now skipped when no row
  is complete, and the warning of class `spatialkit_no_scored_rows`, which
  `select_on = "split"` raised when only the selection half had no complete
  row, is raised here too, with or without the split: no row of the layer
  has a finite value of the response and of every predictor, so `rss`,
  `cp`, `moran_i` and `moran_z` are `NA` at every level. The bounds record
  then carries `n_response_only`, the rows with a finite response, and
  `summary()` and
  `select_resolution()` (`cp`, `moran_z`) name the predictors as the cause.
  The split's own case keeps its text, now without the base R warning.

* **`cell_values_at()` without `seeds` filled every location in silence
  when every polygon of `values` was empty.**  The check that `values`
  carries cell polygons read the geometry type alone, so a `values` of
  `POLYGON EMPTY` rows passed it: no location was in a cell, and each was
  filled from the cell of its nearest scored row with nothing said, where
  seeds without a location and a `values` without polygons are refused.
  A `values` none of whose polygons has a location now meets the same
  refusal ("`values` carries no cell polygons"). A `values` with some empty
  polygons runs as before.

* **`create_voronoi_polygons()` and `build_tessellation(method =
  "voronoi")` on points none of which has a location returned zero cells
  with sf's bare warning, or stopped with GEOS's message.**  With a
  boundary the call gave zero cells, `index` `NA` for every point and only
  "x contains no geometries of specified type"; without one, and for seeds
  from `get_voronoi_seeds(method = "provided")`, it stopped with
  "BufferOp::getResultGeometry distance must be a finite value". Both now
  stop with "no point of `points_sf` has a location ..., so no cell can be
  built", as the triangles and the k-means seeding refuse such points.
  Points with a location that the move to the CRS the cells are built in
  leaves without one (latitudes beyond 90 moved to EPSG:3857) gave the
  same bare messages; they now stop with "the move of `points_sf` to the
  CRS the cells are built in placed none of its points", named by the
  function that made the move. The
  index of points in cells skipped its `spatialkit_unassigned` warning when
  there was no cell; it now gives the count there too, as
  `?build_tessellation` (`index`) says it does whenever a point is `NA`.

* **`get_voronoi_seeds()` ignored `seeds` under `method = "kmeans"` and
  `"random"`, and `sample_points` under `method = "random"`, with nothing
  said.**  It reports every other argument a method does not read
  (`partition`, `init`, the density settings, `n`). The seed of the random
  stream is `set_seed`, so `seed = 1` was taken as `seeds = 1` and the call
  drew from the session's stream in silence: two such calls gave two
  different seed layouts. Each now warns that the method ignores the
  argument, and returns what the call without it returns. A call the
  method refuses (`"random"` without `boundary`, `"kmeans"` with neither
  `boundary` nor `sample_points`) is told the refusal only.
  `?get_voronoi_seeds` (`seeds`, `sample_points`) says so.

* **`create_grid_polygons()` returned an empty grid without a word when the
  clip left no cell.**  The refusal of zero cells ran on the lattice before
  the clip only. In degrees with s2 on, a boundary whose edges span more
  than 180 degrees of longitude is read as the shorter great-circle arcs,
  across the antimeridian, and the strip s2 clips with can miss the lattice
  whole: the box from (-179, -80) to (179, 80) at `cellsize = c(358 / 3,
  10)` gave 0 cells, and `build_tessellation(method = "square")` with it
  gave 0 cells and `index` `NA` for every point. The call now stops with
  "the clip to `boundary` left zero cells", and in degrees names the cause
  and the remedy (a projected `crs`, or `target_cells` or `n`).
  `?create_grid_polygons` (`clip`) says so.

* **`build_tessellation()` on cells in degrees held its snap tolerance in
  degrees against distances in metres, and stopped with "there is no
  package called 'lwgeom'" with s2 off.**  A point outside every cell is
  snapped to the nearest cell within a thousandth of the median cell width,
  but for cells in degrees that width stayed in degrees while the distance
  was in metres on the sphere, so a 0.1 degree grid snapped within 1e-4 m:
  a lattice point 1.2 m outside its cell (s2 reads the cell's edge as an
  arc) was `NA` in `index`. The width is now taken as an arc of the sphere,
  as the margin of an edge already was, about 11 m on that grid, and the
  point is snapped. With `sf::sf_use_s2(FALSE)` sf handed those distances
  to lwgeom, which spatialkit does not depend on, and the call stopped once
  a point lay in no cell. The nearest cell and the distances of the snap
  are now taken on the sphere whatever `sf::sf_use_s2()` says; the join of
  points in cells keeps the session's engine. Projected and CRS-less cells
  are indexed as before. `?build_tessellation` (`index`, `params$snapped`)
  says so.

* **`determine_optimal_levels()` scored a count on other cells when
  `max_levels` changed and the window around the elbow did not.**  Under
  `criterion = "morans_i"` and `"combined"` the integer that seeds the
  re-runs in the window was drawn after the sweep, which takes more draws
  the more counts it runs, so `max_levels` moved it. On 480 points in
  twelve clusters, `max_levels` 30 and 31 both put the elbow at 12 and the
  window at 8 to 16, and `moran_z` at 14, 15 and 16 cells was -0.45, -0.45
  and -0.59 against -1.27, -0.41 and 0.53; `"morans_i"` returned 15, 14, 16
  against 15, 16, 10. The integer is now read before the sweep and the
  stream is put back after the read, so for one `set_seed` a count gets the
  same cells, `wss_eval`, `moran_i` and `moran_z` in any window that holds
  it, whatever `max_levels` and `top_n` are. The sweep, the elbow and the
  geometric ranking are as before; for a given seed the re-runs' scores
  and the `"morans_i"` ranking can differ from before.
  `?determine_optimal_levels` (Details) says so.

* **`resolution_profile()` gave a level other cells when the ladder held
  other levels.**  Every level's partition drew from one random stream in
  ladder order, so a level's cells, and `wss`, `rss`, `cp`, `moran_z`,
  `loo` and the other columns read off them, depended on the levels before
  it: on 480 points in twelve clusters at `seed = 123`, the 20-cell level
  had a `wss` of 106,786 on the ladder 10, 20, 30, 40 and of 106,030 on 20,
  30, 40, and a `moran_z` of 0.65 against -1.02. Dropping low levels that
  did not win moved the count `select_resolution()` chose: on the ladders
  of even counts to 36 from 6, 8, 10 and 12, `criterion = "cp"` chose 34,
  36, 30 and 28. A ladder of sizes moved the same way with the sizes it
  held. Each level now runs on a stream of its own, seeded with one integer
  drawn after the subsample plus the level, and the stream is put back
  after it, so for one `seed` a level gets the same cells on any ladder
  that holds it and is fitted to the same rows (on a layer of more than
  `sample_n` rows a ladder of counts is fitted to a subsample and a ladder
  of sizes to every row); on those four ladders `"cp"` and `"loo"` each choose one
  count. Profiles for a given seed differ from before, and so do the
  selections, seeds and cells built from them. With `seed = NULL` the
  levels no longer draw from the caller's stream; the integer is one draw
  from it. `?resolution_profile` (`seed`) says so.

* **`fit_rf_model()` grew another forest on a response at a large offset,
  and on one whose squares overflow or underflow.**  ranger scores a split
  by sums of squares it forms from the response, so at 1e8 plus a spread of
  1 those squares are near 1e17 and their differences read the rounding:
  the same forest (seed, rows, trees) less the offset gave fitted values up
  to 1.3 from those at the origin, and `cv_rf()` and `compare_models_cv()`
  scored those forests in every fold. At 1e160 the squares overflow, the
  forest was grown on infinite scores and its fitted values were off by
  1.5e160 with no OOB error at all; at 1e-160 the OOB error lost its bits.
  The forest is now grown on the response taken about its centre where
  that centre is more than 1000 times its spread, and divided by a power
  of two where it lies beyond 1e100 or its spread below 1e-100; the centre
  and scale are recorded in `$info$response_centre` and
  `$info$response_scale` and put back on `fitted()`, `predict()`,
  `residuals()`, the OOB RMSE and the importance. An importance whose
  squared units lie beyond the range of doubles is `Inf`, or keeps few of
  its bits, with a WARN line that says so, and `print()` shows it as `Inf`
  with that reason. The permutation importance under
  `scale.permutation.importance = TRUE` and the impurity importance under
  `splitrule = "maxstat"` have no units and are left as ranger reports
  them. Any other response gives the forest it gave. `fit$engine`
  predicts the response less the centre over the scale; `?fit_rf_model`
  (Value) and `?predict.rf_fit` say so.

* **`fit_gwr_model()`, `gwr_model_selection()` and `cv_gwr()` squared a
  response whose squares overflow or underflow.**  The local regressions
  are linear in the response, so the fitted values of a response times
  1e160 or 1e-160 were those of the response times the factor, but the
  residual sum of squares, AICc and the bandwidth search squared the raw
  residuals. On 50 rows times 2^-530, AICc was -36665.629 where -36665.615
  is right, and `gwr_model_selection()`'s second gap was 6.509 where 6.484
  is right. Times 2^530, AICc was `NA`, the search took 50 neighbours where
  the unscaled response gets 23, the fitted values moved by up to 1.01
  times the factor, `cv_gwr()` searched each fold the same way, and
  `gwr_model_selection()` stopped with "no model produced a finite
  criterion". The response, less any centre, is now divided by a power of
  two where its values reach beyond 1e100 or its spread is below 1e-100,
  and the factor is put back on what is returned: the bandwidth is that of
  the unscaled response, the fitted values, residuals, coefficients,
  standard errors, `$engine$lm` and the predictions are those of the
  unscaled response times the factor, AICc (that of `$bandwidth_search`
  too) gains 2 n log of it, and a cross-validation score of the search is
  multiplied by its square. The residual sum of squares, `info$sigma2` and
  such a score, in the response's units squared, are `NA` with a warning
  where a double cannot hold them, and `predict(se.fit = TRUE)` then takes
  its standard errors from the residual standard deviation.
  `fit$engine$centres$scale` gives the factor. Any other response is
  fitted as it was, to the last bit; `?fit_gwr_model` (The local solve)
  says so.

* **`gwr_model_selection()` on a layer with no rows stopped with sp's
  "error in evaluating the argument 'x' in selecting a method for function
  'addAttrToGeom': subscript out of bounds".**  Its own check of the
  number of rows ran after the rows were converted, and only with
  `adaptive = TRUE`: with a fixed bandwidth one or three rows ended in "no
  model produced a finite criterion". The check now runs first, at any
  bandwidth, and gives "0 observations cannot support a local fit of 2
  candidate predictors plus an intercept." With a fixed bandwidth it now
  refuses, as with an adaptive one, any layer with no more rows than the
  model with every candidate has parameters. A constant response, which
  ended in the same "no model produced a finite criterion; nothing can be
  selected.", is refused in the words of `fit_gwr_model()` ("response 'y'
  is constant"), and where no model has a criterion the refusal now names
  the bandwidth and the reason.

* **The boxcar kernel's window moved with the unit of the coordinates on a
  regular grid.**  With an adaptive bandwidth a boxcar window holds every
  observation up to the distance of the k-th nearest, and on a lattice
  many tie with it. Computed from coordinates in US survey feet or in
  millionths of a metre, their distances and the k-th round apart by a
  unit in the last place, and the test `d <= h` left some of them out: on
  a 7 x 7 lattice with 12 neighbours, 16 of 49 fitted values moved, by up
  to 0.026, between metres and US survey feet. A squared distance within
  one part in ten million of the edge now counts as at the edge, the
  package's tie rule, in the local regressions, the bandwidth search and
  the local collinearity survey of `fit_gwr_model()`,
  `gwr_model_selection()`, `cv_gwr()` and `predict()`. The other kernels
  are unchanged.

* **`resolution_profile()` warned what `cp` and `reliability` take from the
  variogram, then stopped with "no usable value in `levels`".**  A
  `levels` with no value the layer can hold (`levels = 500` on 120
  points), with a `sac` that draws a warning or with the variogram
  estimated, raised that warning about columns the call never returned.
  `levels` is now checked before the variogram is estimated or read, and
  the work limit on the search of each scored row's nearest centre (a
  ladder of counts with an `estimate` and without FNN) is applied there
  too: a call either stops with its error alone, or returns the profile
  the warnings describe. A call that returns a profile is unchanged.

* **`estimate_sac_range()` warned that "the variogram is fitted to the RAW
  response instead" on calls that fitted none.**  Where the detrending
  could not run (a predictor missing on every row, or an `lm()` that
  fails, as on a factor with one level) and the response was then
  constant or the points had no extent, the result was the `NA` for that
  reason together with the warning of a raw-response fit. The warning now
  comes after those refusals, on the calls that go on to the raw
  response.

* **`estimate_sac_range(detrend = "reml")` said the REML fit "did not
  converge" where `nlme::gls()` could not fit the trend at all.**  On a
  factor predictor with one level `gls()` fails on its contrasts before
  any iteration, and the warning of the fallback to OLS gave
  non-convergence as the reason. Where no start gives a fit and the
  design of the trend cannot be built on the rows of the fit, the warning
  now says the REML detrending "could not be fitted" and gives the
  design's error.

* **`estimate_sac_range()` and `kriging_adequacy()` squared a response near
  1e160 or 1e-160 as given.**  The squares of such a response lie past
  what doubles hold. Near 1e160 the structure test of
  `estimate_sac_range()` was not made and the semivariances overflowed,
  so the range was refused with "fewer than three lag bins hold point
  pairs", which was false, and `make_folds(auto_range = TRUE)` built
  other blocks; `kriging_adequacy()` refused with the same reason. Near
  1e-160 Moran's I of the structure test and the range moved in their
  fourth or fifth digit, `kriging_adequacy()` kriged no cell, and its
  `se` moved by up to 0.3 percent. With `detrend = "reml"` the REML fit
  failed near 1e160 and the range came from the OLS fallback, with a
  warning that REML "did not converge". The structure test now reads the
  response divided by a power of two near its spread; the semivariances,
  the model and the nugget of `estimate_sac_range()`, its REML fit, and
  the variogram, the kriging and the per-cell variances of
  `kriging_adequacy()`, are computed in units of a power of two near the
  response's spread wherever that lies outside 1e-150 to 1e150, and put
  back in the response's units. A response times 2^-530 or 2^530 gives
  the range, the test, the folds and the kriged values of the response
  itself, times the factor where they carry its unit; values in the
  response's units squared (the variogram, the nugget, the REML variance,
  `kr_var`, `plain_var`) are `Inf` near 1e160, and near 1e-160 are
  subnormal doubles that keep few bits, or 0. A variogram model whose
  sill is `Inf`, passed to `kriging_adequacy()` as `sac`, is refused with
  that reason, where it was refused as having "no positive sill" or for
  "fewer than three lag bins". A response of ordinary size gives the same
  results as before.

* **A response far from zero was read without a word where doubles there
  do not hold a millionth of its spread.**  At 1e15 doubles are 0.125
  apart, so a response of standard deviation 1 stored there is rounded to
  eighths, and so is every mean and prediction computed at that size:
  `make_folds(auto_range = TRUE)` built other blocks from another range,
  and `select_features_forward()` selected a predictor where the same
  response near zero selects none, both in silence. `prep_model_data()`,
  and so every fit and `cv_*()` function, `estimate_sac_range()` and
  `kriging_adequacy()` (with or without `sac`) now warn, with class
  `spatialkit_response_imprecise`, where one unit in the last place of
  the largest absolute value of the response is more than a millionth of
  its standard deviation (a standard deviation of 1 at about 1e10, 2^33,
  or beyond), and give that spacing and the standard deviation. The
  values are read as given. A response whose spread is a few units in the
  last place, which the package reads as constant, is left to the tests
  for a constant. `compare_models_cv()`, `cv_block_size_sweep()` and
  `select_features_forward()` raise it once for the run, and
  `cv_spatial()` once for its folds when its `fit_fn` prepares each
  fold's rows again.

* **`resolution_profile()` and `select_resolution()` gave a false reason
  for `cp` being `NA` where the level its noise variance comes from has a
  residual mean square of 0.**  Without a nugget `cp` takes as its noise
  variance the residual mean square of the finest level whose cells hold at
  least two scored rows on average, and is `NA` at every level when that
  mean square is 0, as it is when each cell of that level holds one
  location whose rows carry one value.  With the nugget at the fit's lower
  bound, the `spatialkit_cp_noise` warning and the print said that no
  level's cells hold two scored rows on average, and
  `select_resolution(criterion = "cp")` named needs the layer met; on 40
  locations of three rows each, at 5 and 40 cells, the 40-cell level holds
  three a cell.  With no usable variogram nothing said that `cp` is `NA`,
  and the warning or log line on the variogram said it takes its noise
  variance from that level.  The warning, now raised in both cases, the print and
  `select_resolution()` name that level, its scored rows a cell and its
  mean square of 0; the warning's `level` is that level, and
  `attr(, "cp_noise")`, `NULL` there before, records it with `source` and
  `value` `NA`, `level`, `rows_per_cell`, `rms` and `reason`.  With
  `predictor_vars` the residuals of such rows differ in their last bits,
  and `cp` took that rounding (a mean square of 2.4e-31 on values of 1) as
  its noise variance, so it had no penalty and `select_resolution()` picked
  the finest level.  A mean square within the rounding of the values scored
  now counts as 0 there, and `cp` is `NA` at every level, as it is without
  predictors.

* **The warnings and log lines of `resolution_profile()` on a variogram it
  cannot use said `cp` "takes its noise variance from the finest level"
  where no level's cells hold two scored rows on average and `cp` is `NA`
  at every level.**  A range below the shortest lag, a variogram with no
  spatial structure, a refused range with the nugget at its bound, a pure
  nugget and a `sac` with no model now give the condition, as the warning
  for a fit that did not converge did: `cp` takes its noise variance "from
  the finest level whose cells hold at least two scored rows on average,
  and is NA at every level when none does".

* **`kriging_adequacy()` stopped with base R's "replacement has 3 rows,
  data has 4" on a `cells_sf` with a cell whose ID is `NA`.**  The
  per-cell counts were tabulated over the cells' IDs, which hold no level
  for `NA`, so they came out one short for each such cell, and none when
  every ID was `NA`. Where the cells with an ID divided the rows evenly
  (two `NA` IDs among four cells), the shorter counts were recycled over
  the rows instead, and the call returned with no error: cells took other
  cells' `n`, `mean`, `se` and `kr_shift`. The counts are now tabulated
  over the cells that have an ID: a cell whose ID is `NA` comes back with
  `n` 0 and `mean`, `se`, `plain_var`, `kr_gain` and `kr_shift` `NA`, and
  is kriged as any other cell. A `cells_sf` in which no cell has an ID
  (`NA`, `NaN` or a factor's `NA` level) is refused in words that name
  `cells_sf`, its ID column and `id_col`. A double `NaN` ID is no ID
  either, on the cells and on the points: such a cell holds no point,
  where a cell and points that all had `NaN` were matched as one cell,
  and every cell `NaN` stopped with "factor level [2] is duplicated". A
  point whose ID is `NaN` is counted in no cell without a warning, as a
  point whose ID is `NA` is; it was warned about as an ID that matches no
  cell. Two cells that share an ID, which stopped with base R's "factor
  level [5] is duplicated", are refused in words that name the column and
  the IDs, as `assign_features_to_polygons()` refuses them. Layers whose
  cells all have their own ID give the same results as before.

* **Under `select_on = "split"` on a subsample, the warning of class
  `spatialkit_response_ceiling` called the selection half's rows of the
  subsample that carry a response all the subsample's rows that carry one,
  and with none of them it printed "about NA times as many".**  The
  criteria read the selection half, so the ceiling counts that half's rows:
  "0 of the 500 points of the subsample carry a response" was said of a
  subsample in which about half did, all in the estimation half, and with
  a response on half the rows at random it gave 141 where 271 carry one.
  The message now says that the count is of the selection half and gives
  the layer's rows with a response in either half. The factor by which the
  layer's cells hold more scored rows than `cell_resp_mean` shows is taken
  from the layer's selection half (`split$n_response`), where it took both
  halves and came out about twice as large, and the sentence that gives it
  is left out when no row of the subsample is scored, with or without the
  split, and `?resolution_profile` gives that factor under the split. The
  ceiling and the condition's fields do not change.

* **`resolution_profile()`'s warning of class `spatialkit_not_supported`
  for a range floor above the ceiling said "the profile runs from 2 to"
  the ceiling over the levels given in `levels`, and named `min_cell_n` or
  the size as the bound where the distinct locations set it.**  With
  `levels = c(4, 5)` it now says the levels asked for are scored as given;
  on 60 rows at 5 stations it says the 5 distinct locations allow at most
  5 cells, where it said "60 points at min_cell_n = 9 support at most 5"
  of a `min_cell_n` that allows 6; and on a ladder of sizes whose finest
  size asks for more cells than k-means can place it names the size's
  count and the locations, where it called the 5 cells the locations
  allow "the 5 cells that size gives" of a size that asks for 20. The log
  line for a subsample too small to reach the floor says the same of given
  `levels`, and so does `?resolution_profile`. The rule, the levels and the
  condition's fields do not change.

* **Where `min_cell_n` and the distinct locations set the same ceiling,
  an optimum at the ceiling was labelled "the support ceiling (n /
  min_cell_n)", and the note advised `levels` past it and a lower
  `min_cell_n`, neither of which moves it.**  On 45 rows at 5 stations
  the ceiling is 5 cells by both bounds: `levels = c(2, 5, 7, 9)` scores 2
  and 5 only, and `min_cell_n = 3` leaves the ceiling at 5. The `edge` of
  `select_resolution()`, its print, `summary()` and the
  `spatialkit_bound_chosen` warning of the functions that take the count
  now name "the support ceiling (as many cells as the distinct locations
  allow)", which offers no remedy. `attr(, "bounds")$ceiling_from` is
  unchanged.

* **Three texts called a default something the caller asked for.**  On a
  call that names no size, `print(sel$count_rule)` (and the print of the
  count record on the seeds and on a summary) read "set by size (4.5 rows
  per cell asked; ...)" while `print(sel)` read "4.5 rows per cell, the
  default"; it now reads "4.5 rows per cell, the default". The subtitle of
  `plot()` of a profile whose finest size gives fewer than two cells read
  "fewer than two cells of the size asked for"; it now names the default
  size as "(obs_per_cell = 4.5, the default)" and the finest of the plain
  mean's default sizes as "(obs_per_cell = 2, the finest of the default
  sizes)". `select_resolution(profile, obs_per_cell = ...)` with
  `criterion` left at "auto", on a profile it reads by "loo" or "cp", said
  that criterion = "auto" "was asked for"; it now says that "auto", the
  default, reads the profile by that criterion. As in the selection's
  print, a size or a criterion is named by its value: 4.5, the ladder of
  the seven default sizes and "auto" are the default whoever names them,
  and any other is the one asked for. The refusal, the counts and the records do
  not change.

* **`estimate_sac_range(cutoff = 1)` advised raising `cutoff` for a range
  past the fitted lags, and a `cutoff` above 1 is refused.**  At
  `cutoff = 1` the lags already run to the farthest pair of points (of
  the `n_max` subsample, when there is one), and the warning now says so
  and keeps the rest of its advice: detrend, or set a block size. Below 1
  it still advises a larger `cutoff`. The value (NA) and `rejected_reason`
  do not change.

* **With `range_frac` below 1, `estimate_sac_range()` said a range inside
  the fitted lags exceeded the largest lag and that the variogram never
  reached a sill.**  The bound is `range_frac` times the largest lag, so a
  range between the two is refused although it lies within the lags. The
  warning now says the range exceeds the largest range `range_frac`
  accepts but lies within the largest lag, and advises a larger
  `range_frac` or a block size set explicitly, as `plot()`'s caption
  already did. A range past the largest lag keeps the old warning. The
  value (NA), `rejected_range` and `rejected_reason` do not change.

## Documentation

* The README, the vignettes, the help's examples, the scripts in
  `inst/scripts/` and the README's figures follow the new defaults, and each
  document says, where it brings in the pipeline, how to get the earlier
  behaviour.  `?summarize_by_cell` gives two results of a second
  simulation, on fields whose level jumps at boundaries: the error of the
  pooled value at observations beside such a boundary, and the coverage of
  the interval for one observation in cells of one observation.

* `?resolution_profile` describes `cp` as the error of re-measuring the
  sampled sites.  It no longer says `cp` estimates the error of predicting a
  new observation by its cell's mean.  The numbers are unchanged.

* `?resolution_profile` says under `sample_n` that the fits of a ladder of
  counts on more than `sample_n` points are those of the subsample, which
  `seed` draws, so the count read off such a profile can move with `seed`
  and with `sample_n`.

* The help, the vignettes and the README no longer say which way a
  cross-validated score errs without saying where the model will predict.
  "Random folds flatter autocorrelated data" (`?spatialkit`), "the blocked
  estimate is the one to report" (`vignette("spatialkit_nc_demo")`), "the
  optimistic bound" (`vignette("spatial-cross-validation")`), "the score is
  still optimistic" (`vignette("reporting")`) and the like in `?make_folds`,
  `?cv_gwr`, `?fit_rf_model`, `?print.rf_fit`, `?residuals.rf_fit`,
  `?select_features_forward`, `?estimate_sac_range`, `?plot_folds`,
  `vignette("diagnostics")` and `vignette("getting-started")` now read as
  `?fold_separation` does: a cross-validated error is the error at the
  distance the folds put between the held-out points and their training
  data, and it reads low for predictions made farther from the data and
  high for predictions made closer.  Random folds give the error of
  predicting among the samples, which is that of a map over an area the
  samples cover evenly; blocks give the error of predicting into unsampled
  gaps of their size.  The vignettes on cross-validation, on reporting and
  on North Carolina read their folds against the locations of a map with
  `fold_separation()`, and the README lists it.  The README no longer says
  that cross-validation holds out regions sized to the autocorrelation
  range: the blocks of the `cv_*()` functions are a geometric grid unless
  `block_size` or `auto_range` sizes them.  No page presents
  `cv_block_size_sweep()` as the way to choose a block size, and
  `?spatialkit` no longer calls the autocorrelation range the number that
  should set it: the sweep gives the error at each block size, and
  `fold_separation()` says which size matches a map.

* The tour scripts `03-folds.R`, `04-block-size.R`, `08-feature-select.R`
  and `example_nc_demo.R` say the same.  `example_nc_demo.R`, the runnable
  version of `vignette("spatialkit_nc_demo")`, no longer marks the blocked
  score as the one to report, and it and `03-folds.R` read random folds and
  blocks against a map of their area with `fold_separation()`, as the
  vignettes do; `04-block-size.R` reads the folds of each size of its sweep
  the same way.

* `?cv_spatial`, `?cv_gwr`, `?cv_bayes`, `?cv_rf` and `?compare_models_cv`
  say what `folds = NULL` estimates: the error at the distance the block
  grid built there puts between the held-out points and their training
  data, which overstates the error of a map where the samples cover the area
  it predicts over by as much as the model's error grows between the two
  distances.  No default changes.

* `?cv_block_size_sweep` and `?plot.block_size_sweep` say that the rise of
  the curve with block size measures the error at that distance from the
  data, which is the error of a map only at locations that far from it, so
  the curve does not by itself say which size to report.  They no longer
  call random folds the leaky reference or say that the height of the rise
  is what the random-fold number overstated.

* `?summarize_by_cell` says what `..se_*` and `..ci_*` cover at each `deff`:
  at the default `deff = 1`, the cell's own mean, when its observations are
  spread through the cell; under `"kish"`, `"variogram"` or a number, the
  cell mean as an estimate of the population mean.  It no longer names the
  grand mean at every `deff`.  The numbers are unchanged.

* The README, the vignettes and the tour scripts describe the behaviour of
  this version: GWR needs `sp` alone, `moran_z` is a diagnostic, the ladder
  starts at 2 with the range floor reported, `kriging_adequacy()` returns
  `plain_var` and `kr_gain`, `make_folds()` has six methods, and
  `estimate_sac_range()` has seven refusals and no directional stand-in.

* `?fit_rf_model` and the caution logged for `include_coords = TRUE` no
  longer say that a forest given the coordinates fails wherever it has not
  been.  The gain from coordinates is largest between sampled points, fades
  with distance and can reverse in a region the data never visited.  The
  advice that stays is to score such a forest with blocked or
  distance-matched folds, never out of bag, and to compare it with the
  forest without coordinates on the same folds; `?cv_rf` runs that
  comparison in its example.

* `?make_folds` and `?cv_block_size_sweep` no longer say that blocks built
  from `block_size` are one to two times that size on each axis.  That holds
  on an axis at least as long as `block_size`; a shorter axis is spanned by
  a single row or column of blocks as long as the extent there.
  `?make_folds` says that the NNDM guarantee concerns nearest-neighbour
  distances, so a learner that draws on several neighbours can still score
  low on clustered samples.

* `?area_of_applicability` has a Weights section with the cost of
  importance weighting (a weak but real predictor is down-weighted by the
  square of its relative importance), and says that the location index marks
  a gap only when it is wider than the distance its threshold allows.

* `?resolution_profile`, `?select_resolution` and
  `?determine_optimal_levels` say what each criterion returns in practice,
  and `?resolution_profile` lists which conditions are classed warnings,
  classed errors, plain warnings and log lines.

* Every figure in the vignettes carries alt text, which is what a screen
  reader announces and the only thing a reader gets when an image fails to
  load.  Each one states what the picture shows and what it is there to
  demonstrate, rather than naming the axes.

* Every help page now ends with a "See also" that leads somewhere; more than
  a third had none, every S3 method among them.  Four families are new:
  **spatial data preparation** (`ensure_projected()`, `harmonize_crs()`,
  `coerce_to_points()`, `clip_target_for()`, `prep_model_data()`), **package
  options and caches** (`spatialkit_quiet()`, the two cache clearers,
  `create_grid_polygons_cached()`), **methods on a fitted model** (the
  `predict()`, `fitted()`, `residuals()`, `coef()`, `print()` and `summary()`
  methods of the three backends, which document one contract between them) and
  **print methods** for the other result objects.  The seed generators, the
  cached grid constructor and `ensure_stable_poly_id()` join **tessellation**,
  `model_metrics()` joins **model evaluation**, `sac_nugget()` joins
  **cross-validation**, `gp_lengthscale_bounds()` joins **model fitting**, and
  `area_of_applicability()` is now in **prediction** as well as
  cross-validation.  The website's reference index moves the two cache
  clearers into the same group, so it and the help pages agree.

* Ten numbered scripts are installed with the package, in
  `system.file("scripts", package = "spatialkit")`.  `00-run-all.R` runs them
  in order; each of `01-` to `10-` is self-contained and covers one topic:
  tessellations, resolution, fold schemes, block sizing, fitting and
  diagnosing, model comparison, prediction surfaces and the area of
  applicability, feature selection, GWR and the Bayesian GP.  They run on a
  simulated field with known structure, print what they are doing, and state
  what to look for in a figure before drawing it.  No script asserts a number
  it did not compute, so a run that finishes is one whose claims held on the
  machine that ran it.  `SPATIALKIT_TOUR_OUTPUT` writes the figures to a folder
  instead of the device; `SPATIALKIT_TOUR_PAUSE = "no"` skips the per-figure
  pause.  A script skips the part that needs an absent package with a message
  naming it, and scripts 02, 09 and 10 skip themselves when gstat, sp or
  brms is absent.
* Four vignettes join `spatialkit_nc_demo`: `getting-started` (installing, what
  the coordinates are in, the pipeline from points to a scored model and a map,
  and a glossary of the terms that recur), `resolution` (the ladder, the
  criteria and why they disagree), `spatial-cross-validation` (the six
  fold schemes, block sizing, and reading a CV result down to the last row) and
  `diagnostics` (residual autocorrelation, aggregation standard errors, kriging
  adequacy, area of applicability, and two ways to leak).  Each is executed at
  build time and gates itself on the optional packages it needs.
* The README is generated from `README.Rmd`, so its figures and every number in
  it are computed when it is built rather than pasted in.  It is about half its
  previous length: the material that had accumulated in it moved to the
  vignettes above, and what remains is what the package is for, a quick start
  that shows a real cross-validation gap, the troubleshooting list, and pointers
  to the rest.
* New vignette `reporting`: what leaves the session at the end of a run.  The
  first half is the regions as a file someone else can use.  Grouping a layer
  you already have with `assign_features_to_polygons(largest = TRUE)` (100
  North Carolina counties into a hex grid of 18 regions, 12 of them
  populated, one row per county), answering
  whether a location falls in one with `keep_unassigned = TRUE`, IDs that
  survive a reprojection, and why the aggregates go into a GeoPackage instead
  of a shapefile: of the 10 columns `summarize_by_cell()` produces, 2 come
  back from a shapefile with their names intact.  The second half is a worked
  report of six numbers, each read out of an object the run already produced,
  with what each one is there to stop a reader believing (Roberts et al. 2017;
  Meyer and Pebesma 2021; Heaton et al. 2019).  Its example reports a run that
  fails its own checks, which is the case the section exists for.
* The documentation is published as a website at
  <https://elkronos.github.io/gis_modeling_toolkit/>: the README, every help
  page grouped by pipeline step, the six vignettes as articles, and this
  changelog.  It is rebuilt from `main` on every push, so it describes the
  development version; the "development version" heading at the top of this
  file lists what the CRAN release does not have yet.

* `select_features_forward()` now says what `$score` is: the cross-validated
  metric of the winning set at the final step, which is the selection
  criterion and is optimistically biased by the selection itself (Cawley
  and Talbot 2010) --- not a performance estimate of the selected model.
  The honest estimate comes from running the selection inside
  `cv_spatial()`'s `fit_fn`, which the page now spells out.
* `fit_bayesian_spatial_model()` documents which `brms` families `family =`
  takes --- zero-inflated and hurdle counts, negative binomial, Bernoulli,
  beta, ordinal, categorical and mixture families all reach `brms::brm()`
  with the spatial GP term intact, and a factor response is taken only
  under `categorical()` and the ordinal families --- and gains a worked
  zero-inflated Poisson example.
  The same page now records a trap: a family object `brms` cannot name skips
  the response-type check entirely rather than falling back to the gaussian
  rule, so a malformed `family` buys less validation, not more.
* `fit_bayesian_spatial_model()` gains a "Spatial confounding" section: a
  coefficient estimated beside a spatial random effect is a different estimand
  from the non-spatial one (Zimmerman and Ver Hoef 2022), can shrink toward
  zero when the response is smoother than the covariate (Bolin and Wallin
  2025), and the honest diagnostic is to report both side by side.  The
  section names both sides of the restricted-spatial-regression dispute
  (Hughes and Haran 2013; Hanks et al. 2015; Khan and Calder 2022) and the
  remedies expressible in a `brms` formula (Marques, Kneib and Klein 2022;
  Guan et al. 2023).
* `model_metrics()` gains a "Which metrics survive a non-Gaussian response"
  section, inherited by `cv_bayes()` and `compare_models_cv()`: RMSE and MAE
  are defined for any numeric response; MAPE, SMAPE and R-squared are
  Gaussian-shaped; CRPS and interval coverage from `cv_bayes()` are the proper
  scores for a count or bounded response.  It also records that the
  all-folds-failed `fold_metrics` frame omits the `coverage_*` columns.
* `estimate_sac_range()` documents what a count or other mean-variance-linked
  response does to the variogram, why detrending with `predictor_vars` helps
  but does not fix it, and what to do instead.
* `?spatialkit` now reads its nine-step pipeline as an argument --- steps 1
  to 4 are the claim, steps 5, 7 and 9 the evidence that makes it checkable
  --- and gains a "Defaults and their sources" section listing which defaults
  rest on a citation and which were chosen, so the two are not mistaken for
  each other.
* `cv_spatial()` documents its name collision with `blockCV::cv_spatial()`,
  which builds folds where this one runs them, and that `blockCV`'s
  `$folds_ids` is accepted directly as `folds` everywhere.

* `?determine_optimal_levels` no longer calls the Cliff and Ord moments
  behind Moran's z exact: they assume cell means of equal variance, and with
  single-point cells beside cells of 70 or more points `z` ran slightly high
  (mean 0.2--0.34, 7--8% rejection at 5%); `resolution_profile()` says
  `cell_diam_median` is about 0.8 of the side of an equal-area square cell,
  not a width (and `vignette("resolution")` no longer calls it one); the
  Post-selection inference section and the vignette say that the estimation
  rows of a split leave cells in the selection half empty and cells across
  the border estimated from part of their points (coverage 0.88 against
  0.97 in simulations), and how to find the cells to trust; `sac_nugget()`
  says the nugget is extrapolated from the first lag bin (about
  `max_dist / 30` at the default cutoff), not observed below the closest
  pair, and can be exactly 0 at the fit's bound.

* The tessellation help pages are corrected: `clip_target_for()` returns
  the points' bounding box, not their convex hull, and is not the target
  `build_tessellation(method = "voronoi")` derives;
  `create_voronoi_polygons()` and `build_tessellation()` say that a
  multi-vertex MULTIPOINT feature gets one cell per vertex and that
  `keep_duplicates` has no effect; `voronoi_seeds_random()` tops a short draw
  up to exactly `k`; `voronoi_seeds_kmeans()` and `get_voronoi_seeds()` say
  that their `stats::kmeans()` partition is not the k-means++ one
  `resolution_profile()` scored, and no longer promise equal counts per
  cell.

* Documentation corrected: `summarize_by_cell()`'s `sac` no longer
  identifies a rejected fit by a `status` attribute no `sac_range` carries
  (it is an `NA` value with a `rejected_reason`), and says a rejected `sac`
  is set aside and the variogram estimated as if none had been given;
  `ensure_stable_poly_id()` and the reporting vignette say IDs are the same
  across projections except for fine cells within the rounding step of one
  longitude (36 of 2,500 100 m cells via EPSG:3035), not always;
  `summarize_by_cell()`'s "Confidence intervals" says `..neff_*` is `NA` for
  a single observation; the resolution vignette says the blocked split
  reduces the leak between halves rather than removing it.

* `?estimate_sac_range` now says that only the exponential family is fitted
  and so overestimates the range on smoother fields, and why the family is
  not chosen by fit; that the "decreases with distance" refusal also fires
  on small samples of ordinary fields, which its warning says below 100
  points; that the longest directional range is read from
  `directional_fitted` after checking `directional_status`, not
  `max(attr(, "directional"))`, which is `NA` exactly when the major axis
  ran past the fitted lags; and that an accepted range can still exceed half
  the width of the layer, leaving `make_folds(auto_range = TRUE)` room for
  one block (`range_frac`).  It has sections on the test for spatial
  structure, the lags and structure at two scales, and lists the seven
  `rejected_reason` values; the README counts seven refusals.

* `?make_folds`: the NNDM details now say how the procedure differs from
  `CAST::nndm()` (a strict removal rule, one point more conservative per
  distance value, and ties broken by coordinates rather than row index)
  instead of calling it the same; the n > 5000 refusal gives the worst-case
  cost, O(n^3) time where `min_train` binds (about nine minutes at
  n = 3000), not O(n^2); `drop_empty_blocks`, `boundary` and
  `block_multiplier` say what they do on the cases above.

* `?cv_rf` no longer says a `seed` passed through `...` overrides the
  per-fold forest seed (`seed` is the function's own argument and never
  reaches `...`); `?cv_bayes` says `parallel = n` compiles `n` Stan models at
  once, at several GB each, and what a fold killed for memory reports;
  `?residual_morans_i` and `?compare_models` say their methodological
  cautions are logged, not raised as warnings.

* `?fit_gwr_model`'s "Collinearity diagnostics" section described the
  30-location unweighted spot check the survey replaced, a global index on
  the predictors alone, and a caveat that only a subset of locations is
  examined; it now describes the code (every location, kernel-weighted, a
  global index on the centred predictors, and at each location one index
  with the intercept and one for the slopes alone).  `?fit_gwr_model` and
  `?gwr_model_selection` now state the adaptive bandwidth's floor and cap.

* `?fit_bayesian_spatial_model`: `check_convergence` lists the checks that
  set `$info$convergence_ok` and the warning a failed verdict raises; a
  section "What the basis resolves" describes the floor of the basis, the
  check before compiling and the adequacy status after sampling, which
  changes neither `convergence_ok`; a Warnings section lists the condition
  classes; the `family` argument and
  the non-Gaussian section no longer promise that any response type brms can
  fit works here; the spatial confounding section warns that coefficients
  under `standardize_predictors = TRUE` are per standard deviation before
  comparing with `lm()`.  `?coef.bayesian_fit` gains a section on
  standardised predictors, and `print()` on such a fit names them.
  `?fitted.bayesian_fit` no longer says a failed posterior draw returns `NA`
  (it has been an error since before this release).
  `?gp_lengthscale_bounds` says its bounds are the prior's calibration range
  and do not shrink with `n`.

* `?area_of_applicability` states the outlier rule (type-7 quartiles) and how
  the threshold differs from CAST's (the fence itself, capped at the largest
  training DI, which is larger whenever `n_outliers > 0`), with the
  `threshold =` value that reproduces it; the internal note that CAST uses
  `boxplot.stats()` was out of date, and the package page no longer implies
  the threshold matches CAST.  `n_new` is documented as every row of
  `newdata` (`n_inside + n_outside + n_na`), not the rows that passed the
  finite-value filter.  `?predict_surface` and the README say that
  `se = TRUE` on a `bayesian_fit` gives the SD of the mean surface and that
  `type = "predict"` gives the predictive SD.

* The diagnostics vignette's area-of-applicability figure alt text no longer
  calls the training curve cross-validated (no folds are passed there);
  `?plot.aoa` and `?plot.feature_selection` describe the training DI and the
  rejected last step as drawn.

* The examples that stamp EPSG:32632 (UTM zone 32N) on simulated points now
  put the points inside that zone.  The README quick start, the
  `fold_separation()` example, the resolution, diagnostics and spatial
  cross-validation vignettes, and the fixture the scripts in `inst/scripts/`
  share placed them at x and y between 0 and 1000, which is on the equator
  about 4.5 degrees east, outside the zone.  They are now offset by 500000 m
  east and 5000000 m north, as the other examples already were.  Every
  printed result is unchanged, since the package works in planar units; only
  the coordinates themselves, and the graticule on the maps, differ.

* The examples of `resolution_profile()`, `select_resolution()`, `summary()`
  and `plot()` on a profile use `set.seed(4)`, and their comments describe
  either outcome for C_p (an interior level or the support ceiling) and
  reliability on the range floor.  The `plot()` example no longer promises
  four panels on points that have no elbow.

* Tour script 02 runs to the end again: step 02.6 drew the elbow's pick,
  which its evenly spread points no longer have, and now draws the C_p pick.
  It prints the bound each criterion's optimum sits on, points to
  `min_cell_n` and `range_floor` rather than `n_levels` for moving one, and
  computes what the `select_on = "split"` comparison shows (both picks on
  the range floor) instead of calling the difference tuning.  Script 08 no
  longer calls its held-out half untouched or its score gap the selection
  effect, and says its block size of 300 is under the range of about 357.

* The README's resolution figure labels its middle count as the one
  `resolution_profile()` rates most reliable, with the width of that
  criterion's flat region; it was `determine_optimal_levels()`'s count,
  which on those evenly spread sites the function now says the ladder chose.
  The troubleshooting entry quotes the fallback messages
  `determine_optimal_levels()` now gives, with their causes.

* The resolution vignette's criteria table describes the elbow as the
  log-log sag it is, `NA` on points with no cluster structure, and its
  introduction says `determine_optimal_levels()` returns `NA` there.

* `?coerce_to_points` (`tmp_project`) states the rule a CRS-less layer is
  read by: the lon/lat heuristic of `ensure_projected()`, not merely lying
  inside the lon/lat envelope.

* `?summarize_by_cell` now gives the derivation and the measured coverage of
  the small-sample rescaling that goes with every data-derived design
  effect, which the package page said were there; the package page lists it
  among the defaults that were chosen, not cited.

* `?estimate_sac_range` now says that the largest directional range never
  stands in for the all-pairs fit, and how the adaptive lags resolve a range
  shorter than one equal lag bin, which equal bins returned several times
  too long.

* `?make_folds` no longer says `auto_range` "fits directional variograms to
  account for anisotropy".  It sizes blocks from the omnidirectional range,
  and for a field known to be anisotropic the page now points to
  `directional_fitted`.  An accepted range too wide for two blocks makes
  `make_folds()` stop with an error, which the page now says; it had claimed
  the grid "does not collapse to a single block".  The log line announcing a
  lowered `k` just before that error is gone.  The page no longer says NNDM
  never pushes a point's nearest-neighbour distance past `phi`: the last
  exclusion can take it past by one neighbour step, as in `CAST::nndm()`.
  The description lists all six methods.  The page now says that `k = 1` is
  raised to 2 by the k-fold methods, and that only `block_kfold`
  returns `params$blocks_supplied` and `params$boundary_supplied`.

* The README's entry for "response 'y' is not numeric" notes that
  `fit_bayesian_spatial_model()` takes a factor or character response under
  `brms::categorical()` and the ordinal families.

* The README, `?summarize_by_cell` and the getting-started, diagnostics,
  reporting and North Carolina vignettes now say which estimand the
  design-effect-corrected standard error is for.  It is the standard error
  of a cell mean as an estimate of the population mean.  For a cell's own
  mean, which is what a map reports, the default `deff = 1` standard error
  is the right one when the points are spread through the cell.  Several of
  these pages had presented the correction as the right standard error for
  the cells themselves, and there it is too wide by `sqrt(deff / (1 - rho))`
  (4.6 at 20 points a cell and `rho = 0.5`).  The getting-started pipeline,
  which maps the cells, now aggregates at `deff = 1`.

* Smaller corrections: the README's troubleshooting list adds
  `"compare_models_cv(): no recognised model requested."`, which is the
  error when no requested name is recognised, and says that
  `"no viable models."` means every recognised backend is uninstalled.  The
  getting-started install table no longer says `patchwork` is needed for the
  `plot_*()` functions.  The North Carolina vignette explains why some
  local designs of its GWR fit have a high condition index with the
  intercept (an uncentred `elevation`, nearly collinear with the intercept
  inside each window), and why the fit raises no collinearity warning: the
  slope index stays below 30, so the slopes are well determined and only
  the local intercepts are not.  Its fold-map alt text now describes each
  blocked fold as whole blocks in separate parts of the state, not as one
  contiguous area.

* `vignette("diagnostics")`: the "Two ways to leak" example of selection
  inside the folds leaked itself (blocks about 250 m across against an
  autocorrelation range of about 330 m) and its learner could not fit the
  intercept-only model, so `tol` did not apply to the first variable.  It
  now passes `block_size = 400`, fits `z ~ 1` for an empty predictor set,
  and shows `sel$history`.

* `?fit_rf_model` recommended
  `area_of_applicability(weights = pmax(fit$info$importance, 0))` without
  condition.  It now says this works only when that importance is finite,
  which it is not when no row is out of bag.

* `vignette("diagnostics")` said that `assign_features_to_polygons()`
  gives a point on the boundary between two scored cells the cell the
  profile scored it in whenever the cells carry the record of the centres
  scored.  It does so for rows that have a CRS exactly when the record
  names one.  The vignette now says so and, as `?summarize_by_cell` does,
  names the rows the polygons alone place: rows without a CRS beside a
  record that names one, rows with a CRS beside a record that names none,
  rows assigned to the cells after they lost the record, and rows given
  their cells by another join.  `?summarize_by_cell` and the
  `spatialkit_cells_differ` warning said the polygons alone place every row
  without a CRS; rows without a CRS on cells whose record names none take
  the profile's cell, and the help and the warning now give the condition.

* `vignette("diagnostics")` said the inner folds of
  `select_features_forward()` "default to the scheme you pass".  They are
  `block_kfold` unless the function's `method` names another scheme,
  whatever folds the outer loop uses, and the vignette now says so.

* `?print.spatial_fit` said that nothing the print shows is an
  out-of-sample score, and did not list the LOOIC a `bayesian_fit` prints
  when the fit computed it (`compute_loo = TRUE`, the default).  It now
  lists the LOOIC and says it is an approximate leave-one-out score of
  single rows, not a spatial one.

* `vignette("getting-started")` said the UTM zone distorts the distances of
  `nc.shp` by 0.20 percent at worst, where `ensure_projected()` measures
  and logs 0.21 percent.  The vignette now writes both figures from
  `attr(counties, "crs_choice")`, rounded as the log rounds them.

* `?estimate_sac_range` and `?sac_nugget` said `nugget_at_bound` marks a
  nugget "within 1e-4 of the sill".  It marks a nugget of 0, or of at most
  10^-4 of the sill, and both pages now say so, as `?resolution_profile`
  does.

* `?sac_nugget` said a nugget equal to its floor, the variance between
  visits, keeps `nugget_at_bound` `FALSE`.  The mark does not read the
  floor: a nugget at a floor of at most 10^-4 of the sill is marked
  `TRUE`, and the page now says so.  `?estimate_sac_range` described a
  range past `range_frac * cutoff * max_dist` as one beyond the longest
  lag fitted, which holds only at `range_frac = 1`; it now says so.

# spatialkit 2.0.0

Everything below is relative to **1.0.0** (published on CRAN 2026-08-07).
The major bump is warranted: three exported functions are removed, and several
defaults change the result of a fit or a comparison, so the same script can get
a different answer. Both are under "Breaking changes" — read that section
before upgrading a running analysis.

These notes describe what changed **for a user of 1.0.0**. A good deal of this
release was written after 1.0.0 and then revised before shipping; defects that
existed only between those points are not listed, since no released version
behaved that way. The commit history has that record in full.

Throughout, *raises a warning* means a genuine R `warning()` — one
`tryCatch(warning = )` catches, `suppressWarnings()` suppresses and
`options(warn = 2)` escalates. *Logs a warning* means a `logger` message in the
`"spatialkit"` namespace, which none of those touch.

## Breaking changes

### Corrections that change results

Each item below was measured and independently reproduced before it was
touched; the figures quoted are from those reproductions, so you can judge
whether an item affects an analysis you have already run.

* **`residual_morans_i()` no longer puts weight on a point's own residual.**
  `FNN::get.knn()` reports a point's OWN index among its neighbours whenever
  exact duplicate coordinates are present, which put `1/k` on the diagonal of a
  matrix Moran's I is only defined for with a zero diagonal. On 40 sites x 4
  repeats with a response carrying **no** spatial structure, 120 of 160 rows
  gained a self-weight, mean I came out at +0.086 against E[I] = -0.0063, and
  **77% of samples were "significant" at p < 0.05** against a nominal 5%.
  Repeat observations at one site are exactly what
  `make_folds(method = "leave_location_out")` is for, so this was a mainstream
  input. The dense fallback never had the fault, so the statistic also depended
  silently on whether **FNN** happened to be installed; the two paths now share
  one neighbour lookup and agree exactly.

  Requesting `k + 1` neighbours and dropping self is **not** sufficient on its
  own — the slot self occupied displaced a genuine co-located neighbour and left
  a farther point standing in for it (75 of 400 retained pairs sat at distance
  121 where a neighbour at distance 0 existed). Duplicate coordinates are now
  grouped and answered exactly.

* **`residual_morans_i()` gains a `null` argument, defaulting to `"auto"`.**
  Model residuals are not exchangeable — they are orthogonal to the design
  matrix — so the classical randomisation moments are wrong for them. At
  n = 120 with six smooth covariates and independent errors, OLS residuals had
  mean I = -0.031 against the exchangeable E[I] = -0.008, and the z-score
  averaged -0.54 with sd 0.90 instead of 0 and 1. The Cliff & Ord (1981,
  sec. 8.3) regression-residual moments restore mean z = -0.09, sd 1.03 and a
  4.3% rejection rate against a nominal 5%, and **agree with
  `spdep::lm.morantest()` to machine precision** (verified at 1e-16 through the
  public function). `"auto"` applies them only when the fit's residuals really
  are the OLS residuals on the rebuilt design, which a forest's and a working
  GWR's are not; the null actually used is reported in the return value.

* **`summarize_by_cell()` standard errors under a design effect were too
  small.** `s / sqrt(n / deff)` corrects the mean's variance for clustering but
  leaves `s^2` biased low by the same clustering: for exchangeable correlation
  rho, `E[s^2] = sigma^2 (n - deff)/(n - 1)`. The two errors compound. Measured
  95% CI coverage at n = 20: **0.905 at rho = 0.3, 0.796 at rho = 0.6, 0.632 at
  rho = 0.8**; after rescaling by `sqrt((n-1)/(n-deff))`, 0.952 / 0.952 / 0.953.
  Applies to `deff = "kish"`, `deff = "variogram"` and a fixed numeric `deff`.
  **The default `deff = 1` path is bit-identical to before.**

* **`determine_optimal_levels()` ranks model-aware candidates on the
  standardised deviate, not on |Moran's I|.** E[I] and Var(I) both depend on the
  cell count, so |I| shrinks as k grows whether or not the finer tessellation
  captures anything. Over 300 replicates of a response with **no** spatial
  structure, mean |I| fell monotonically from 0.114 at k = 10 to 0.050 at
  k = 60 — an |I| ranking prefers the largest candidate for arithmetic reasons
  alone. Candidates are now ordered by |z| using the Cliff & Ord residual
  moments (exact here, since the cell-level residuals are OLS residuals by
  construction); over the same runs z had mean ~0, sd ~1 and a 5% rejection rate
  of 0.040-0.057 at every k. The `"diagnostics"` attribute now carries
  `moran_z` alongside `moran_i`.

* **`estimate_sac_range()` sweeps four azimuths, not two.** A +/-22.5 degree
  window around 0 and 90 covers exactly **90 of the 180 distinct azimuths** —
  every direction between 23 and 67 degrees, and between 113 and 157, fell into
  neither. On simulated fields with 3:1 anisotropy and a true major-axis range
  of 300, the estimate came back at 255 and 249 for major axes at 0 and 90
  degrees but **151 and 147 at 45 and 135**. Since
  `make_folds(auto_range = TRUE)` sizes blocks from this number, a diagonally
  oriented field silently got blocks half as wide as the correlation they were
  meant to separate. `c(0, 45, 90, 135)` tiles all 180 azimuths; the same fields
  now return 255 / 245 / 249 / 228. A direction whose variogram never reaches a
  sill is excluded rather than taken as a long range, and the `directional`
  attribute now has four named entries.

* **`fit_bayesian_spatial_model()`'s calibrated length-scale prior never
  reached Stan.** `brms::set_prior(spec, class = "lscale")` with no `coef` is a
  *global* prior, and brms applies a global prior only to coefficients with no
  individual prior of their own — every `lscale` coefficient always has one.
  brms dropped it with a note and Stan received brms's defaults, which made
  `gp_lengthscale_bounds()`, the tail calibration and `$info$gp_lscale_prior`
  dead weight. Confirmed with `brms::make_stancode()`: the requested prior is
  absent under the global form and present under the coefficient-level form,
  which is now used. `$info$gp_lscale_prior` is read back from
  `brms::validate_prior()`, so it records what brms will actually use.

* **The GP basis was sized against the wrong domain measure.** brms builds the
  boundary as `choose_L(x, c) = c * max(1, max(x) - min(x))` over the pooled,
  column-centred covariates — the **full range**, not the per-axis half-range in
  which Riutort-Mayol et al. state their inequalities. Recovering the boundary
  from `make_standata()`'s eigenvalues confirms `L = c * full range` exactly at
  every `c`, so the old convention built a boundary **twice as wide** as `gp_k`
  was sized for: the GP was under-resolved, and `$info$gp_ell_min` — the
  diagnostic meant to catch exactly that — was twice too lenient to fire. The
  `c` floor is now brms's own default 1.25 rather than 1.2.

* **`fitted()` on a `gwr_fit` could return a coefficient surface.** The search
  for GWmodel's fitted-value column matched the whole candidate name vector with
  `%in%` and took the first hit in the *SDF's* column order — and the local
  coefficients come first. A predictor named `fit`, `pred`, `prediction`,
  `fitted` or `yhat` therefore returned its own coefficient column, silently:
  executed in-sample R^2 was **-1.18** against a true 0.981, and `residuals()`,
  `summary()`, `model_metrics()`, `compare_models()` and every `cv_gwr()` fold
  consumed it without a warning. The search now runs in preference order and
  excludes any candidate that is also a model term; all five colliding names now
  give R^2 = 0.981, identical to the renamed control.

* **`coef.gwr_fit()` returned GWmodel's whole SDF data slot** — 15 columns for a
  two-predictor fit, of which 3 are coefficients and the rest are standard
  errors, t-values, the response, the fitted values, residuals and `Local_R2`.
  It now returns the model terms only; reach for `object$engine$SDF` for the
  rest.

* **`estimate_sac_range()` is reproducible, and no longer disturbs the caller's
  RNG.** `seed` now defaults to `123L` rather than `NULL`. The `n_max`
  subsample is an internal approximation, not part of the answer, and leaving it
  unseeded made the returned range differ between runs on identical input
  (19531 / 19589 / 19605 on three calls) while silently advancing the caller's
  stream — and `make_folds(auto_range = TRUE)` sizes its blocks from that
  number. Pass `seed = NULL` for the old behaviour.

* **`estimate_sac_range()` rejects a non-numeric response.** `as.numeric()` on a
  factor returns its level codes, so a factor response produced a variogram of
  an arbitrary integer relabelling of the categories and the estimated range
  changed when the levels were reordered (3700 against 2497 on the same data).
  Factors and character columns are now an error naming the column; logicals are
  read as 0/1.

* **Fold sets built from a different dataset are refused.** Fold splits are
  lists of `..row_id` values, and row IDs are `seq_len(nrow())` unless supplied,
  so passing `cv_gwr()` a `folds` object built from another dataset of the same
  size applied cleanly — every ID matched, every fold was populated, and the
  model was scored on splits describing other observations. `make_folds()` now
  records a small projection-invariant row fingerprint in
  `params$row_probe`, and `cv_gwr()`, `cv_bayes()`, `cv_spatial()` and `cv_rf()`
  error rather than proceed. Fold objects from earlier versions carry no
  fingerprint and are passed through unchecked.

* **`evaluate_insample()` rejects duplicated names in `fits`.** `model` is the
  key `compare_models()` joins its metric and Moran's I tables on, so two fits
  called `"GWR"` produced a 2x2 cross-join: four rows, every one carrying the
  first fit's numbers, with the second fit never scored at all.

* **`fit_gwr_model()` rejects a non-numeric predictor.** `gwr.basic()` expands
  contrasts via `model.matrix()` and fits, but `gwr.predict()` does not and
  fails, so the model appeared to fit and then silently predicted all `NA`.

* **`fit_gwr_model()` no longer rejects a two-valued continuous response.** The
  "binary" error is now gated on the response being integer-like. A
  left-censored or saturated measurement (every observation at a detection limit
  or a ceiling) has two distinct values and is perfectly continuous; it now
  warns instead. The guard also runs once per fold inside `cv_gwr()`, where a
  small training fold can legitimately hold only two distinct values.

* `fitted()` returning the wrong length, or nothing, is now an error in
  `summary()` and `model_metrics()` rather than a plausible row count over an
  all-`NA` comparison. `new_spatial_fit()` is the documented extension point, so
  a subclass with a missing or mis-sized `fitted()` method is user-reachable.

* The cached `fitted()` on a `bayesian_fit` is stamped with the `n` and a digest
  of the data it was computed from. The cache environment has reference
  semantics — which is what makes it survive copy-on-modify — so `fit2 <- fit`
  gave both objects the *same* cache, and assigning different data to the copy
  returned the original's values at the original's length.

* `make_folds()` drops rows with empty or non-finite coordinates, with a logged
  warning naming the count, rather than letting an `EMPTY` POINT reach
  `block_kfold`'s nearest-block rescue and die with "replacement has length
  zero".

* When every fold fails, the warning now names the first underlying error.
  Previously "all 5 folds failed" was the whole diagnosis even when the cause
  was simply that **brms** or **GWmodel** was not installed.

* `make_folds()` records the CRS the folds were built in as `params$crs`.
  Geographic input is projected by `ensure_projected()` to a CRS the caller
  never chose, and `block_size`, `sac_range` and `buffer` are lengths in *that*
  CRS.

* `.morans_i_for_k()` returns `NA` at or below nine cells, where every cell
  neighbours every other and Moran's I collapses to exactly `-1/(k-1)` for any
  residual vector — a function of the cell count alone.

* `residual_morans_i(fit, k = 1)` works on machines without **FNN**. `apply()`
  simplified the length-1 result to a vector, making the neighbour index a
  1 x n matrix and every row after the first out of bounds.

* **`summarize_by_cell(deff = "kish")` under-estimated the predictor ICC by
  about a factor of `m`.** The pooled one-way ANOVA grouped the `m` z-scored
  predictor columns under the same cell label, so independent per-column cell
  effects averaged away in the shared cell mean and the between-cell sum of
  squares shrank by ~1/m. Measured at true rho = 0.5 with m = 4: pooled ICC
  **0.12** against 0.495 per column, so every predictor SE was ~44% too small.
  The pooled group is now (variable, cell), which recovers 0.49.

* **Design effects are built from what a column actually observes.** A cell
  of 10 rows with 2 finite responses had its response SE formed at the
  10-row design effect, then applied to a 2-observation mean: adding 8
  NA-response rows moved the SE from 8.46 to **26.38**. Each column's design
  effect now uses its own non-missing count, with the mean pairwise
  correlation recomputed over the observed locations when a cell has NAs;
  `cell_weight` is the effective count of the primary variable, not of rows
  (`n` still counts rows).

* **CRS-less coordinates get ONE interpretation, wherever they enter.**
  `prep_model_data()` assumed EPSG:4326 for CRS-less data that looked like
  lon/lat and projected it, while every `predict()` method passed the fit's
  CRS as a target — a branch that *stamped* it onto the raw numbers. The same
  rows sat in two different places, and `predict(fit, newdata = training
  rows)` disagreed with `fitted(fit)` by up to one response SD (R² 0.98
  in-sample, 0.64 via newdata). The heuristic is now a single function used
  by both branches; a fit records the assumption it was built under and
  replays it on CRS-less `newdata`. Two further symptoms of the same split —
  CRS-less LINESTRINGs aborting in `coerce_to_points()` ("crs not found"),
  and hex/square `build_tessellation()` refusing input voronoi accepted — are
  fixed with it. The assumption is now announced with a real R warning.

* **`residual_morans_i()` refuses a malformed `weights` matrix** instead of
  silently substituting the default k-NN(8) matrix (I = 0.874 returned for
  four malformed shapes against 0.805 for the weights actually supplied).

* **The fold-provenance fingerprint no longer refuses the caller's own data.**
  Three defects in the version introduced last pass: a character `..row_id`
  was coerced to all-`NA` (and matched row 1 everywhere); coordinates were
  compared as `"%.7g"` strings and flipped on the ~1 in 5000 that a
  reprojection moved by 5e-9°; and polygon input was probed *after*
  pointization, so a different `pointize` in the cv call read as different
  data. Both sides now probe the geometry as supplied, keep IDs in their own
  type, and compare numerically within 1e-6°.

* **`n_folds_attempted` counts the folds supplied.** A fold whose test rows
  were all removed as incomplete vanished before fitting and was absent from
  both counts, so five supplied folds reported `4/4`. It is announced with a
  real warning.

* **`determine_optimal_levels()`'s elbow uses the signed deviation below the
  chord.** `abs()` let a concave bump *above* the chord — a k where k-means
  fell into a worse local optimum — win with the same magnitude.

* **`estimate_sac_range()` returns `NA` for a constant variable** (an exactly
  explained response, or a constant one) instead of a fitted "range" of 168
  or 673 from a variogram that is identically zero. `make_folds(auto_range =
  TRUE)` no longer re-opens the unseeded subsample by forwarding its own
  `seed = NULL`.

* **`fit_bayesian_spatial_model()` attaches a user-supplied *global*
  `lscale` prior at coefficient level**, the same way it does its own, so
  `set_prior(..., class = "lscale")` reaches Stan instead of being discarded.

* **`fitted()` on a `gwr_fit` returns the prediction when a predictor is
  named `prediction`.** The model-term exclusion added last pass was applied
  to `gwr.predict()`'s SDF too, where coefficients are suffixed `_coef` and
  the column literally named `prediction` *is* the prediction; `predict()`
  returned all `NA`.

* Smaller: a logical response meets the same binary-response guard as `0/1`;
  `model_metrics()` errors on a non-numeric response instead of returning
  `n = 0`; `fitted.bayesian_fit()` errors when the posterior cannot be drawn
  instead of returning silent `NA`; `create_grid_polygons()` refuses a grid
  above `max_cells` (default 1e6) up front;
  `make_folds(method = "buffered_loo")` states its guard in bytes (splits are
  ~4n² bytes; the old n = 20000 cap admitted 1.6 GB); `-0` and `0` are the same
  coordinate in the duplicate-aware k-NN; `.morans_i_for_k()` returns the `NA`
  pair whenever the moments are unavailable.

* **Documented warnings are now R warnings.** Eight paths the manual
  described as warning only wrote a logger line, invisible to
  `tryCatch(warning = )`, `expect_warning()` and `options(warn = 2)`:
  `residual_morans_i()` returning `NULL`, the CRS assumption in
  `ensure_projected()` and stamping in `harmonize_crs()`, GWR collinearity,
  seed clamping in `voronoi_seeds_kmeans()` / `get_voronoi_seeds()`, dropped
  rows in `ensure_stable_poly_id()`, and the dropped column in
  `assign_features_to_polygons()`. Three deliberate methodological cautions
  (`include_coords = TRUE`, `random_kfold` feature selection, non-standardised
  Moran weights) stay logged and their documentation now says so.

* **`predict()` on a Bayesian GP fit depended on which other rows shared the
  call.** brms 2.x stores `Xgp`, `dmax` and `cmeans` in a fit's GP basis but
  **not** the Hilbert-space boundary `L`, so `brms:::.data_gp()` recomputes it
  from whatever rows `predict()` is handed. Every eigenfunction of the
  approximation therefore moved with the newdata bounding box while the fitted
  basis coefficients stayed put. Measured: `L` was 5.57 at fit time, 4.02 for a
  five-row `newdata` and 3.63 for one row; `predict_surface()` on the same
  10,000-cell grid differed by **1.78** between `chunk_size = 5000` (the
  documented default) and a single call; and `cv_bayes()`, which predicts each
  test fold separately, scored **every** fold against a basis the model was
  never fitted with. `predict.bayesian_fit()` now pins the boundary by
  appending the training coordinate extrema and dropping them again, so the
  chunked, fold-wise and single-call answers are identical.

* **A predictor whose name is not a syntactic R name fitted a different
  model.** Every backend builds its formula from these names, so a column
  called `"B5-B4"` was fitted as `B5 - B4` — a different model, silently, with
  `$predictor_vars` still reporting `"B5-B4"` — and `"band 4"` died inside
  `str2lang()` with a parser error naming no column. `prep_model_data()` now
  refuses such names up front and says which column and what to do about it.
  (Backticking is not a fix: the `sp` coercion GWmodel needs runs the names
  through `make.names()` anyway, after which the formula and the data
  disagree.)

* **`summarize_by_cell(deff = "kish")` weighted the response with a
  predictor's ICC.** The fallback to the predictor ICC is for the case where
  no response was supplied; applying it whenever the response's own ICC came
  out non-positive meant that merely adding a predictor to the call changed the
  response's regression weight by a factor of **20**, while the response SE was
  (correctly) left at `deff = 1`.

* **A column literally named `n` was summarised as the row count.**
  `dplyr::summarise()` makes each new column visible to the expressions after
  it, so the row count shadowed a response or predictor named `n`:
  `resp_mean_n` came back equal to `n`, with `sd` and `se` `NA`, and no
  message.

* **Variogram models are fitted with a nugget.** A nugget-free model forces the
  curve through the origin, and gstat's default `N/h^2` weights buy that
  constraint by collapsing the range: with a 50% nugget the fitted range came
  back at about **0.45** of the truth, so `make_folds(auto_range = TRUE)` built
  blocks less than half the correlation length it reported.

* **The directional maximum is used only when anisotropy is established.**
  Splitting 180 degrees four ways leaves each directional variogram about a
  quarter of the point pairs, and the maximum of four noisy estimates is biased
  upward: on isotropic simulated fields the returned range ran about **40%**
  above the truth and the "notable anisotropy" warning fired on the majority of
  them. An omnidirectional variogram is now fitted alongside and is the default
  answer; the directional maximum is used when all four directions fit, their
  ratio exceeds 1.5, **and** the widest stands more than 1.5x above the
  all-pairs estimate. On the package's own test field (true range 80) the old
  rule returned 248 and the new one returns 84; genuine 3:1 anisotropy is still
  recovered.

* **`predict()` replays the CRS decision the fit was made under, including a
  negative one.** When CRS-less training data were passed through as planar,
  nothing recorded that,
  so every `predict()` re-ran the lon/lat heuristic on `newdata` alone — and a
  **subset** of those same training rows, whose own bounding box sits inside
  the lon/lat envelope, was judged differently from the whole: taken for
  degrees, reprojected, and predicted about 1e6 m from where it was fitted.
  `predict(fit, training_subset)` disagreed with `fitted(fit)[subset]` by more
  than the response's standard deviation while `predict(fit, full_data)` agreed
  exactly.

* **Antimeridian data are no longer flattened onto EPSG:3857.** A bounding box
  cannot distinguish global coverage from a layer straddling +/-180 degrees, but
  the coordinates can (one very large gap in the sorted longitudes). Web
  Mercator *splits* such a layer: two stations 41 km apart came out **40,068
  km** apart, destroying every distance downstream. Only genuinely global
  coverage now falls back to EPSG:3857.

* **The projection for a wide extent is chosen by measurement, not by rule of
  thumb.** Albers standard parallels came from the bounding box while the
  conic/azimuthal branch came from the centroid latitude, so a
  trans-equatorial extent put the parallels either side of the equator: at
  `lat_1 = -lat_2` PROJ refused the string outright and `ensure_projected()`,
  `make_folds()` and `build_tessellation()` all aborted with an
  internal-looking "invalid crs", and just short of it the projection distorted
  distances by **15.6%** against the single UTM zone's 1.65%. The candidates
  are now scored by projecting a sample of the data's own points and comparing
  planar with geodesic distances — WGS84 ellipsoidal distances (Vincenty), not
  `sf::st_distance()`'s s2 sphere, whose 0.24–0.56% gap from the ellipsoid is
  the size of the errors being ranked and mis-ordered the candidates on 16 of
  40 random wide extents; the message reports both figures.

* **`fitted.gwr_fit()` could return a coefficient surface.** Which GWmodel call
  produced an SDF was inferred from its column names, so a *predictor* ending
  in `_coef` made a `gwr.basic()` SDF look like a `gwr.predict()` one and
  switched off the model-term exclusion; with a second predictor named `yhat`,
  `fitted()` then returned that predictor's local coefficient surface (R2
  **-2.28** against 0.986) and `residuals()`, `summary()` and `model_metrics()`
  followed it. The mode is now passed explicitly.

* **A logical predictor supplied as text was silently mis-coded.**
  `is.numeric`, `is.factor` and `is.character` are all `FALSE` for a logical, so
  the type guard skipped it entirely: a `"TRUE"`/`"FALSE"` character column —
  what a CSV round trip produces — was factor-coded 1/2 against splits built on
  0/1, sending every row to the `TRUE` side (correlation **0.21** with the
  correct predictions).

* **`compare_models_cv()` scores every model on identical folds.** The
  documented guarantee only held when `folds` was supplied: with `folds = NULL`
  each backend built its own from `k`, `seed`, `block_size` and the rest, so a
  per-model `block_size` put **104 of 150 rows** in different folds for GWR and
  RF, and `seed = NULL` did the same with no overrides at all. The folds are now
  built once, and the arguments that decide the split are protected.

* **User weights with a non-zero diagonal.** Every moment of Moran's I assumes
  no observation is its own neighbour. A self-inclusive row-standardised kNN
  matrix — an easy thing to build by hand — rejected the null on **75%** of
  white-noise residuals at a nominal 5%, with no condition raised. The diagonal
  is now zeroed with a warning.

* **`build_tessellation()$index` no longer snaps outside points to the nearest
  cell.** Points outside every cell were assigned to whichever was closest, so a
  summary built from the index counted all 40 points of a layer whose study
  area held 10, while `assign_features_to_polygons()` on the same cells
  correctly reported 30 misses. Only a point within a thousandth of a cell
  width is snapped now; the rest are `NA`, and the count is logged.

* **`ensure_stable_poly_id()` is stable across projections again.** The sort
  key was the raw double centroid, so cells sharing an exact `x` in one CRS
  differed by ~1e-11 degrees after a round trip through another: **14 of 16**
  cells got a different ID depending on which CRS the layer arrived in — the
  exact failure the function exists to prevent. The key is now rounded to 7
  decimals (about a centimetre).

* **A numeric `deff` is applied as the uniform inflation it is documented to
  be.** The `E[s^2]` correction the estimated design effects apply is derived
  from within-cell correlation and is unjustified for a constant the caller
  chose: `deff = 2` doubled a 3-point cell's SE instead of multiplying it by
  `sqrt(2)`, and returned `NA` for every cell with `n <= deff`.

* **The Kish ICC guard matches its documentation.** The docs promise "at least
  2 cells with 2+ observations; falls back to `deff = 1` otherwise", and the
  check was only `k >= 2, N >= 4`: all-singleton cells give a within-group sum
  of squares of 0 and therefore an ICC of exactly **1** — a design effect of
  `n` from data carrying no within-cell information at all.

* **`gp_lengthscale_bounds()` and the GP basis use two dimensions and unique
  locations.** `stats::dist()` uses every column, so POINT Z geometry gave
  3-D length-scale bounds in mixed units; and `brms::gp()` defaults to
  `gr = TRUE` and reduces its covariates to unique rows before taking the
  boundary, so replicated locations made the package's `S` — and with it
  `gp_c` and the basis-adequacy threshold `gp_ell_min` — wrong by that factor
  (0.68x with one heavily-sampled station).

* **`estimate_sac_range()` and `make_folds()` drop the Z dimension.**
  `gstat::variogram()` and `sf::st_distance()` use every coordinate dimension,
  so an XYZ layer had its elevation folded into each lag: the returned "range"
  was a length in 3-D (413.6 against 136.8 for the same stations) while the
  block grid, the buffered-LOO buffer, NNDM's neighbour distances and
  `summarize_by_cell()` all work in 2-D map distance.

* **`determine_optimal_levels()` returns the elbow first.** Under the
  geometric criterion — the default, and the fallback every model-aware call
  takes below the nine-cell floor — the candidates came back sorted ascending,
  so the knee sat in the middle: `k[1]` and `top_n = 1`, both documented as
  "the top-ranked candidate", returned **knee − 1** on every such call, and the
  help example answered `1` for two clearly separated clusters. The vector is
  now knee first, then its lower and upper neighbours, so position 1 means the
  same thing on both paths. The quick-start data's answer moves from `3 4 5`
  to `4 3 5`. Same code in 1.0.0.

* **`estimate_sac_range()` is invariant to rotating the layer.** Two things
  were not. The lag cutoff was a fraction of the bounding-box *diagonal*, a
  property of the axes rather than of the points, which grows by up to √2 when
  the same layer is rotated 45°; it is now a fraction of the farthest pair
  (found on the convex hull), as the documentation always said. And the
  directional maximum is **never** preferred over a usable all-pairs fit. The
  1.0.0 code returned the widest of two axis-aligned directions; the three
  hurdles later put in front of it (all four directions fitted, ratio above
  1.5, maximum above 1.5× the all-pairs fit) still let the noise through — on
  one isotropic field rotated in 10° steps they "established" anisotropy in
  **14 of 18** orientations and returned ranges from 225 to 529, a 2.35×
  spread produced by nothing but the direction the axes pointed. The four
  directional ranges are still reported, as a diagnostic, in the
  `directional` attribute; `anisotropy_used` is `TRUE` only when the all-pairs
  fit itself failed. A field *known* to be anisotropic should size its blocks
  from `max(attr(range, "directional"))` explicitly, and the log line says so.

  The variogram fit itself no longer depends on gstat's single starting
  value. `fit.variogram()` starts the optimiser at a third of the longest lag;
  for a field whose range is a small fraction of the extent that start is ten
  times too long, and whether the iteration landed or collapsed to a singular
  model depended on floating-point details — the same 250-point field fitted
  on Linux and came back singular on an arm64 Mac, where the estimate then
  fell through to the directional maximum (315 against a true range of 80).
  Shorter and longer starting ranges are tried as well, the converged,
  non-singular fit with the smallest weighted sum of squares (gstat's own
  criterion) is kept, and a converged spherical fit is preferred over an
  exponential one that did not converge. When no model fits at all — a flat,
  nugget-only variogram — the `NA` now carries the empirical variogram, so
  `plot(fit, type = "variogram")` draws it and says why there is no range
  instead of refusing.

* **`residual_morans_i()` splits tied neighbour distances.** The k-nearest
  neighbour weights broke ties at the k-th distance by whichever point came
  first — in row order on the dense path, in kd-tree order under **FNN** — so
  the same data in a different row order gave a different I, z and p whenever
  observations shared a location (repeat visits), and on gridded data the
  answer depended on whether FNN was installed. Points tied at the k-th
  distance now share the remaining weight equally, on both paths, which is the
  only rule that is a function of the geometry alone; the matrix is
  row-standardised as before. Where no distances tie the weights equal
  `spdep`'s k-NN weights exactly; on repeat-visit data every co-located twin
  is retained with its share, and the statistic no longer changes when the
  rows are shuffled or when **FNN** is installed. Same code in 1.0.0.

* **`fit_gwr_model()`'s collinearity diagnostic is the scaled condition
  index, thresholded at 30.** It was `kappa()` on the raw, unscaled predictor
  matrix at 1e6 — a number that depends on the predictors' units, so
  rescaling a column changed it and the threshold was a threshold on nothing:
  a design with a scaled condition index of 1322, whose local coefficients ran
  from −86 to +150 around a true value of 2, raised nothing. Belsley's index
  (every column scaled to unit length, intercept included, the ratio of the
  largest to the smallest singular value) is now computed for the global
  design and for each sampled local window, and the conventional 30 is the
  threshold in both. Expect the warning on designs that were silent before.

### API and default changes

* `estimate_sac_range()` returns `NA` instead of a number when the fitted range
  runs past the longest observed lag, or the optimiser stopped at its iteration
  limit. Such a range is *unidentified*, not long — the empirical variogram
  never reached a sill — and sizing blocks or design effects with it is worse
  than declining to. The refusal carries `rejected_range` and
  `rejected_reason` attributes and keeps the fitted variogram, so
  `plot(type = "variogram")` can draw exactly the case worth looking at; the
  result is classed `sac_range`, so it prints as a bare `NA` rather than
  dumping the fit to the console. Callers that fed the old number straight into
  `make_folds(block_size = )` now need to handle `NA` — that is the point.

* Removed the legacy wrappers `evaluate_models()`, `evaluate_models_cv()` and
  `phi_prior_bounds()`. Use `compare_models()`, `compare_models_cv()` and
  `gp_lengthscale_bounds()`.

* `compare_models_cv()` gains an `"RF"` branch and an `rf_args` argument, so a
  `ranger` forest can be compared against GWR and the Bayesian GP on identical
  folds. Unrecognised model names now raise a warning and are dropped, and a
  request with nothing recognised left is an error. Previously a bare
  `intersect()` discarded anything outside `c("GWR", "Bayesian")` and fell back
  to GWR, so `models = "RF"` silently ran GWR and reported it as the answer.

* `coef()` on a `spatial_fit` now either returns coefficients or errors; it
  never returns `NULL`. `coef.gwr_fit()` and `coef.bayesian_fit()` used to
  return `NULL` on failure, which is indistinguishable from "this model has no
  fixed effects", so `lapply(fits, coef)` quietly produced a short answer.
  `coef.rf_fit()` errors as before — a forest has no coefficients; use
  `fit$info$importance`. See `?new_spatial_fit`.

* Every `predict()` method errors when the number of predictions does not match
  the number of rows that survived cleaning. It used to recycle silently: two
  predictions for four clean rows produced a four-row answer.

* `create_grid_polygons_cached()`'s default `type` is now `"square"`, matching
  `create_grid_polygons()`. The two disagreed, so cached and uncached calls
  built different grids from the same arguments.

* `build_tessellation()` keeps `cell_id` on `"hex"` and `"square"` grids instead
  of deleting it after indexing, so all four methods return the same ID column
  and `plot_tessellation_map(fill_col = "cell_id")` works on a grid. `poly_id`
  is retained alongside it. `params$expand` now echoes the value passed rather
  than a hard-coded `0`; the grid and triangle methods still ignore `expand`,
  which is now documented.

* `clip_target_for()` projects lon/lat input before applying `expand`. The
  fraction-of-extent form was computed from a bounding box in degrees and handed
  to `sf::st_buffer()`, which reads `dist` as metres. The returned clip target
  is therefore in the projected CRS, not the input CRS, and says so.

* Voronoi, grid and triangle clipping union the boundary first. Against a
  multi-feature boundary, `st_intersection()` split every straddling cell into
  one row per boundary feature and grafted the boundary's attribute columns onto
  the result.

* `predict.bayesian_fit(newdata = NULL)` honours `summary`, `type` and `draws`.
  It short-circuited to `fitted()`, which caches epred column means and nothing
  else, so `summary = "median"` returned means and `type = "predict"` returned
  expected values — silently.

* `prep_model_data()` accepts `predictor_vars = character(0)`, making
  intercept-only spatial GP models reachable from
  `fit_bayesian_spatial_model()`. `fit_gwr_model()` and `fit_rf_model()` reject
  an empty set explicitly.

* `fit_bayesian_spatial_model(control = )` is *merged* over the package defaults
  (`adapt_delta = 0.9`, `max_treedepth = 12`) rather than replacing them.
  Passing `list(max_treedepth = 15)` silently dropped `adapt_delta` — the
  setting the divergence warning tells you to raise.

* `cv_bayes()`'s `predictive_coverage` is averaged across folds **weighted by
  each fold's `n_pred`**. The per-fold values are means over that fold's test
  rows, so an unweighted average is not the pooled quantity once fold sizes
  differ — and `block_kfold` tolerates a 3:1 imbalance before it even logs a
  warning.

* The three seeding functions (`get_voronoi_seeds()`, `voronoi_seeds_kmeans()`,
  `voronoi_seeds_random()`) all emit `seed_id` and `method` columns, so they are
  drop-in interchangeable. `get_voronoi_seeds()` returns seeds in the boundary's
  CRS for every method; per-branch alignment had made the final alignment block
  unreachable.

* Geometry-type checks require **every** geometry to be an accepted type, not
  merely one of them. A mixed POINT/POLYGON layer passed a POINT-only check.

* `ensure_projected()` errors on a `target_crs` that does not resolve to a
  usable CRS, rather than returning the input unchanged and letting unprojected
  coordinates flow into distance and area computations.

* `fit_bayesian_spatial_model()` derives the GP basis count (`gp_k`) and
  boundary factor (`gp_c`) from the ratio of the estimated length-scale to the
  domain size, rather than from the number of observations.

  `brms::gp()` builds a full tensor grid, so `gp(..x, ..y, k = gp_k)` carries
  `gp_k^2` basis functions — `gp_k` is the count *per dimension*. The previous
  rule reduced to `max(15, floor(sqrt(n)))` for any `n` above 45, making
  `gp_k^2` identically `n`: at n = 10,000 the model carried 10,000 basis
  functions and an n × n design matrix, at which point the approximation was no
  longer approximating anything.

  Across the scenarios in `dev/baseline-structural.rds` the derived value is
  22–24 per dimension and largely independent of `n`: at n = 2,000 the basis
  count falls from 1,936 to 576, at n = 10,000 from 10,000 to 529, and at
  n = 200 it *rises*, 15 to 23 — a correction, not an optimisation, so results
  move in both directions. `gp_c` was hard-coded at 1.5, too small whenever the
  length-scale exceeds roughly half the domain half-range; the derived value
  ranges 2.85–3.59 over the same scenarios. Pass `gp_k` and `gp_c` explicitly
  to restore the old behaviour. For scale, a 4-fold cross-validated fit at
  n = 2,000 with 2 chains × 1,000 iterations takes 1,186 s on the reference
  machine at cross-validated R² 0.927 (`dev/baseline-accuracy.rds`, 2026-08-20);
  no comparable timing was captured for 1.0.0, so none is quoted.

* The GP term is built with `scale = FALSE`, changing the default result of
  every Bayesian fit. `brms::gp()` otherwise rescales its covariates so the
  maximum pairwise distance is 1 and reports `lscale` in that space, while this
  package standardises the coordinates itself and expresses the length-scale
  prior, `gp_c` and the basis adequacy threshold in those units. The two
  normalisations differed by roughly the maximum pairwise distance (~4.9 for
  standardised 2D coordinates), leaving the automatic prior about five times too
  diffuse — a likely contributor to divergent transitions and rejected initial
  values. There is now exactly one coordinate scaling.

* The GP fits one length-scale per coordinate axis (`gp_iso = FALSE`), a second
  change to default Bayesian results. Coordinates are standardised per axis, so
  a single shared length-scale made the kernel anisotropic in the original CRS
  by the ratio `sd(X)/sd(Y)` — a property of the sampling layout, not of the
  process. Pass `gp_iso = TRUE` for the previous behaviour; cost is unchanged,
  since the tensor grid is `gp_k^2` either way.

* The automatic GP length-scale prior is a calibrated inverse-gamma rather than
  `normal(0, sd)`, a third change to default Bayesian results. A half-normal on
  a positive parameter puts its mode at zero, so most of its mass sat at
  length-scales shorter than the basis can resolve — where the Hilbert-space
  approximation develops a funnel and the sampler diverges. The replacement pins
  1% of its mass below the estimated lower bound and 1% above the upper. The
  two tail conditions have one exact solution — a one-dimensional root in the
  shape — and that is how it is found, so the calibration succeeds for any
  bounds with `upper > lower`; degenerate bounds fall back to the half-normal
  with a logged note. The prior applied is recorded in
  `$info$gp_lscale_prior`.

* `ensure_projected()` no longer forces continental-extent data into a single
  UTM zone. Transverse Mercator scale error grows quadratically with distance
  from the central meridian, so data spanning the contiguous United States
  carried distance errors of roughly 7.5% near the extent edge, propagating
  silently into `estimate_sac_range()`, `make_folds(block_kfold)` block sizing,
  GWR bandwidth selection and the GP length-scale. Extents reaching more than 5°
  from the candidate zone's central meridian now receive an equal-area
  projection centred on the data, with a logged explanation. Only longitude
  offset triggers the switch, since `cos(lat)` shrinks the distance from the
  central meridian and a tall narrow north-south extent is UTM's design case.
  (Which equal-area projection is no longer decided by latitude band — see
  *The projection for a wide extent is chosen by measurement* below, which
  also replaced the EPSG:3857 fallback for wide bounding boxes with
  antimeridian detection.) Pass `target_crs` to
  override.

* Core counts follow the session's `mc.cores` opt-in, and are capped.
  `fit_bayesian_spatial_model()`'s documented default was
  `cores = max(1L, parallel::detectCores() - 1L)`, and `cv_*(parallel = TRUE)`
  auto-detected the same way with no cap — 63 workers on a 64-core host, and a
  hard error wherever `_R_CHECK_LIMIT_CORES_` is set. The Bayesian default is
  now `getOption("mc.cores", 1L)`; the auto-detect path is capped by that
  option when it is set; and every worker count, explicit or not, is capped at
  the machine's core count (with a message) and at two under `R CMD check`.

### Guards, messages and stricter input handling

* `prep_model_data()` now drops rows whose geometry is empty or whose
  coordinates are not finite, and counts them in its existing log line.
  `st_geometry_type()` calls an EMPTY POINT a "POINT", so nothing ever looked
  at the coordinates: such a row reached GWmodel as a raw `sp` coercion error
  naming no row, `ranger(include_coords = TRUE)` as "Missing data in columns",
  and brms at predict time as an infinite GP boundary that made **every** basis
  function on the surface `NaN`. It also refuses a response listed among its
  own predictors, which was leakage in the forest, a silently reduced model in
  GWR, and duplicated rows plus a phantom `<none>` entry in the GWR selection
  table.

* `make_folds()` validates `k`. A non-integer `k` truncated inside
  `rep(floor(n/k), k)`: only `floor(k)` folds were built, the last rows of the
  permutation landed in **no** test set, and `k` was echoed back unchanged, so
  `length(folds) != k`.

* `make_folds(block_kfold)` refuses a grid it cannot build. There was no cap on
  `nx * ny`, so `block_size` in the wrong unit asked for 1e8-1e11 cells and
  exhausted memory; the message now names the implied cell count, the extent
  and the CRS units, as `create_grid_polygons()` already did.
  `predict_surface()` gained the same guard on `cell_size` and `n_cells`.

* Hand-built `folds` are checked. Train and test overlapping is not
  cross-validation — the model is fitted and scored on the same rows and the
  result is reported as a CV score (RMSE 0.50 against 0.97 for the same data
  properly split) — and is now refused; fold IDs that name no row were dropped
  silently by `na.omit(match())` and are now counted and logged.

* The fold provenance probe tolerates a row it could not measure. It was taken
  before the empty-geometry filter, so it carried `NaN` coordinates and every
  later `cv_*()` call on the same data died with R's internal "missing value
  where TRUE/FALSE needed".

* `estimate_sac_range()` drops rows with unusable coordinates instead of
  returning `NA` for the whole layer under the message "variogram model fit
  failed", which blamed the fit rather than the row and disagreed with
  `make_folds()` on the same data.

* GWR's local collinearity spot-check now runs **after** the bandwidth is
  chosen, includes the intercept column, and treats a non-finite condition
  number as extreme. On the default path (`bandwidth = NULL`) it used a
  stand-in window, so the documented warning never fired; without the intercept
  it could not see an indicator constant inside a window; and
  `is.finite(cn) && cn > 1e6` discarded exactly singular designs. A new post-fit
  warning counts local regressions that returned non-finite coefficients —
  previously `fitted()`, `summary()` and `model_metrics()` silently reported
  metrics computed from the survivors (n = 18 of 200, R2 = 0.96).

* `fit_gwr_model()` warns when a fixed bandwidth is implausibly small for the
  data's extent. The bandwidth is a distance in the CRS the fit runs in, which
  `prep_model_data()` may have chosen: 0.2 supplied for lon/lat data is 0.2
  metres, and every local window came back empty with nothing raised. The
  argument's documentation now says so.

* `summary.spatial_fit()` applies the same response-type guard as
  `model_metrics()`. A character response still produced `n = 0` and all-`NA`
  metrics there, and a factor died inside `abs()`.

* `print.spatial_fit()` prints the CRS (which its documentation has always
  promised) and one `Formula` line. `sprintf()` vectorised over a multi-element
  `deparse()`, so every `bayesian_fit` printed its formula as two mangled
  fields — the `gp()` term this package builds is always long.

* `summarize_by_cell()` warns when `cells_sf` carries duplicated IDs (each
  summary row is then repeated per matching cell, so `sum(n)` exceeds the
  number of points), and realigns **every** per-cell vector in
  `deff_applied` after the join, not only `deff`: `$rbar` was left in pre-join
  order, so `deff[i]` and `rbar[i]` described different cells.

* `assign_features_to_polygons()` warns when no feature falls in any polygon
  instead of returning an empty layer silently.

* `create_grid_polygons(target_cells = )` builds square cells for
  `type = "square"`. `cellsize = c(w/nx, h/ny)` forced an exact bbox tiling, so
  the cells were rectangles (aspect 10.5 on a 1000:1 strip). For `type = "hex"`
  a differing `cellsize[2]` is now collapsed with a warning before the
  `max_cells` estimate, which used both components and was therefore off by
  their ratio.

* `create_grid_polygons_cached()`'s renumbering is documented: it applies
  `ensure_stable_poly_id()` and `create_grid_polygons()` does not, so the same
  cell carries a different `poly_id` from the two builders.

* `get_voronoi_seeds(method = "kmeans")` clusters in two dimensions and drops
  rows with unusable coordinates, matching `voronoi_seeds_kmeans()`. It used
  every column of `st_coordinates()`, so POINT Z geometry was clustered in 3-D
  with elevation dominating, and an EMPTY POINT crashed inside `kmeans()`.

* `determine_optimal_levels()` refuses a factor or character response, which
  `as.numeric()` silently turned into level codes: re-ordering the levels of
  the same factor changed the chosen number of levels and every `moran_z`.

* `plot_tessellation_map()` brings CRS-less layers into the plot's CRS instead
  of passing them through to fail inside `ggplot_build()` at print time with
  sf's message, naming no layer; and it borrows a CRS from an overlay when the
  tessellation itself has none.

* Documentation corrected where it did not match behaviour: `ensure_projected()`
  states how it chooses a projection for a wide extent (by centroid latitude
  when that entry was written, by measured distortion in this release) and that
  the choice is announced rather than silent;
  `.looks_like_lonlat()`'s two tests are a disjunction and the extent test
  decides first, so a small planar survey inside the lon/lat envelope IS taken
  for degrees — the trade and its reasoning are now stated;
  `harmonize_crs()` no longer claims to match `ensure_projected()` while doing
  something else; `summarize_by_cell()` states which estimand its standard
  errors are for (the grand mean, where measured coverage is 0.95, not the
  cell's own mean, where the naive SE is the better estimate) and that the
  variogram path applies one correlation function to every column.

* **A misspelt `newdata` is an error, not an in-sample answer.**
  `model_metrics()`, `evaluate_insample()` and `compare_models()` forward `...`
  to `predict()`, which checks it only on the out-of-sample branch, so
  `model_metrics(fit, newdta = hold)` silently took the in-sample branch and
  returned an RMSE of **1.086** where the held-out answer was **25.24**, with
  the same return shape. Arguments in `...` with no `newdata` are now refused
  by name. `predict()` on a `gwr_fit` or a `bayesian_fit` likewise refuses
  unknown arguments instead of swallowing them; `predict.rf_fit()` accepts
  only ranger's own predict arguments through `...`.

* **`predict()` enforces one `newdata` contract on all three fit classes.** A
  bare `sfc` died inside two of them with R's "argument must be coercible to
  non-negative integer"; a numeric-at-fit predictor that arrived as character
  (a CSV round-trip) was refused by name by `rf_fit` and `bayesian_fit` and
  returned all-`NA` with a generic backend warning from `gwr_fit`; a missing
  column was reported by two different functions in two wordings. All three
  now run the same check first, so the message is the same whichever fit is
  behind it.

* `make_folds()` validates `block_size` the way it validates `k`: a single
  finite positive number, else an error naming the argument. `NA` and a
  length-2 vector used to die as internal R errors, and a negative, zero or
  character value was silently ignored — yet echoed back in
  `params$block_size` as if it had been used. The grid-size guard also formats
  its own message: a `block_size` in the wrong unit could ask for a grid past
  2³¹ cells on a side, which `%d` refused with "invalid format" instead of the
  documented refusal.

* `estimate_sac_range()` says which variable it modelled. A predictor name
  absent from the data was dropped silently by `intersect()`, and with every
  name unknown the **raw** response was modelled — so
  `make_folds(auto_range = TRUE)` sized blocks from a range of 88.5 instead of
  the residual range of 362.2 with nothing said. An unknown predictor is now
  an error, as it is everywhere else; a detrending fit that fails (an all-`NA`
  column) raises a warning and falls back to the raw response; and a new
  `detrended` attribute records which was used.

* `make_folds(block_kfold, boundary = )` refuses a boundary containing none of
  the points — almost always two layers in different places, a CRS that could
  only be stamped — and raises a warning counting the points that fall outside
  a boundary that contains some, since the region is silently extended to
  cover them. The single-block error names what produced the grid (the block
  size, `block_nx`/`block_ny`, or the automatic grid) rather than always
  blaming `block_nx`/`block_ny`, which the caller may never have passed.

* Rows that no fold names are reported. The `folds` ↔ data guard was
  one-directional: fold IDs naming no row were counted and dropped, but rows in
  the data that appear in no fold's train or test set passed with no condition
  at any level — a `folds` object built on `site[1:45, ]` and applied to all
  90 rows scored 45 of them and reported `n_folds_attempted =
  n_folds_succeeded = 3`. Every `cv_*()` now raises a warning with the count
  and an example row ID.

* User-facing conditions name the function the user called. The cross-
  validation path leaked two internal names into ordinary console output —
  `.remap_folds():` and `.cv_run_folds():` — and `cv_rf()`, a wrapper around
  `cv_spatial()`, reported every message, warning and error in
  `cv_spatial()`'s name. The `sf`-input assertion shared by the tessellation
  and seeding functions said "Expected an sf object" with no function named;
  it now names the caller and, when handed a whole `build_tessellation()`
  result, says to pass its `$cells`. `coerce_to_points(mode =
  "line_midpoint")`'s MULTILINESTRING refusal is prefixed like every other.

* `clear_grid_cache(cache_env = )` removes only its own entries. It removed
  every binding in the environment it was handed and counted them all as
  "entries removed", so a user who passed a project environment lost
  unrelated objects. Cache keys now carry a `spatialkit_grid::` prefix and
  nothing else is touched.

* `summarize_by_cell(deff = "variogram")` now uses **every** structured
  component of a nested variogram model, each weighted by its partial sill,
  which is the correlation the model implies (`1 - gamma(h) / sill`). It read
  the single largest component, so a user-built `Nug + Exp + Sph` model gave
  a correlation of 0.108 at 200 m where `gstat::variogramLine()` implies
  0.197, and the design effects and standard errors with it. Models from
  `estimate_sac_range()` are single-component and are unaffected. A model of
  a family the function does not implement (Matern, power, circular, ...)
  was silently read as exponential; it now falls back to `deff = 1` with a
  warning that names the family. Exponential, spherical and Gaussian are
  supported.

* `cv_bayes()$predictions$yhat_sd` is the posterior predictive standard
  deviation of each held-out row, from the same draws that give the coverage
  columns. It was an unconditional `NA` placeholder. It stays `NA` when
  `compute_pred_intervals = FALSE` or the draws failed for a fold, and it is
  documented.

* `fit_gwr_model()` refuses `n <= p + 1` observations with one error, before
  touching the backend. Such a fit cannot have a residual degree of freedom
  in any window; `n = 2` used to warn three times ("only 2 observations",
  "fallback bandwidth", "2 of 2 local regressions singular") and return a fit
  whose fitted values were all `NA`.

* `create_grid_polygons_cached()` gains `max_entries` (default 50): once the
  cache holds that many grids, adding one evicts the earliest-added. It never
  evicted, so a loop over a thousand boundaries held every grid (about 2 MB
  per 2,500 cells) for the life of the session. Nothing but the grids is
  written into a caller-supplied `cache_env`; the insertion order lives
  inside the package. The cache key also hashes the package version, so a
  cache that outlives an upgrade cannot serve a grid built by an older
  `create_grid_polygons()`.

* `spatialkit_quiet()` accepts a `logger` threshold as well as `TRUE`/`FALSE`,
  and the value it returns can be passed back: `old <- spatialkit_quiet();
  spatialkit_quiet(old)` restores exactly the level that was in force.
  `spatialkit_quiet(FALSE)` put back the package default (WARN) whatever had
  been set, and the returned value was refused as "must be TRUE or FALSE".

* `.onUnload()` disarms both `logger` appenders, so the `spatialkit` logger
  namespace no longer keeps pointing at the session's temp-file path after
  `unloadNamespace("spatialkit")`. `.onLoad()` re-registers them.

* `residual_morans_i()` documents the second reason `"residual"` (and
  `"auto"`) falls back to the randomisation null: fewer than four residual
  degrees of freedom, where the residual variance formula divides by
  `(n - p)(n - p + 2)`. `"residual"` logs a warning when it does, `"auto"`
  does not, and `df` is then `n - 1`.

* `estimate_sac_range()` states what "effective range" is for each model --
  three times the range parameter for the exponential fit (95% of the sill),
  the range parameter itself for the spherical fallback (100%) -- and its
  return-value documentation matches the code: the no-model case (both fits
  singular) returns the classed `NA` with `rejected_reason` set, and only the
  cannot-even-start cases (no `gstat`, too few values, no variance, a
  degenerate extent) return a bare `NA`. Its example now demonstrates a
  fitted range on a simulated field with a known one (3 x 100 = 300), and a
  refusal on a field whose range the data cannot pin down.

* `make_folds()` documents what `block_multiplier` does (the automatic grid
  aims for `block_multiplier * k` blocks, so each fold holds out about that
  many; 3 is a compromise, not a published constant), cites Roberts et al.
  (2017) and `blockCV` (Valavi et al. 2019) for sizing blocks from the
  autocorrelation range, and notes that `blockCV` takes the fitted variogram's
  range *parameter* where `auto_range` takes the *effective* range -- three
  times that parameter for an exponential fit, so larger blocks. `phi` for
  `method = "nndm"` is explained as Mila et al. (2022) define it: the
  autocorrelation range beyond which matching is unnecessary, which
  `estimate_sac_range()` supplies.

* `summarize_by_cell(deff = "kish")` says which ICC estimator it is (ANOVA
  with Donner's `n0`, not REML) and how far the two can differ on an
  unbalanced draw, so the difference is not read as a defect;
  `area_of_applicability()`'s training-DI sentence now says what the code
  does (each fold's actual training rows, not "everything outside the fold");
  `?spatialkit` no longer lists `coef()` among the methods all three backends
  share (a forest has none); the internal elbow helper no longer claims to
  match Kneedle, which it does not on shouldered curves; the GP basis
  diagnostic's code comment attributes its 10% posterior-mass trigger to this
  package rather than to Riutort-Mayol et al. (2023).

* `summary()` on a fit prints `R^2` and `Adj R^2` in ASCII with aligned
  labels (the superscript two rendered as `R<U+00B2>` on non-UTF-8 consoles,
  and `Adj R²=` had no space); `print()` on a random forest likewise; the six
  console messages that carried an em dash use `--`. `coef()` on an `rf_fit`
  prefixes its error `coef.rf_fit():` like its siblings.

* `get_voronoi_seeds(method = "kmeans")` sizes its candidate cloud in double
  precision; `50L * as.integer(n)` overflowed to `NA` above 42,949,672 seeds.

* README: the opening leakage example says it needs `ranger`; the
  no-viable-models example is assigned so that it does not print every fold
  table on a machine that has the backends; the `auto_range` and `cv_bayes()`
  failure examples show every line the console actually prints; the
  installation table no longer suggests installing `loo` separately (`brms`
  installs and calls it); the "logged note" from `determine_optimal_levels()`
  is identified as an INFO-level line in the session log file, not console
  output; the roxygen2 sentence no longer names a version or a `RoxygenNote`
  field.

* The memory-guard test for `make_folds(block_kfold)` stubs out
  `sf::st_make_grid()` for the two refused calls, so a regression of the guard
  fails fast by name instead of attempting an 8000 x 8000 grid.

* Every `cv_*()` refuses a fold list whose splits carry no `train`/`test`
  element, with an error that names the problem. It read `f$train` / `f$test`
  straight, got `NULL` for both, and built empty folds -- so the row-coverage
  warning fired and blamed folds "built on a different or subsetted layer",
  which was not the cause, and the run returned an all-`NA` `overall` with
  `n_folds_succeeded = 0`. `area_of_applicability()` had always refused the
  same input by name; the two now agree. When the splits look positional (two
  unnamed vectors each) the error says so and points at the fold label vector
  instead.

* The `folds` argument of every `cv_*()` documents all three shapes it has
  accepted since the label vector was added earlier in this pass -- a
  `make_folds()` result, a list of `list(train =, test =)` splits, or a vector
  of fold labels -- where the help listed only the first two. The label vector
  is what makes folds from another package usable directly:
  `blockCV::cv_spatial()` returns one as `$folds_ids` (its `$folds_list` holds
  two unnamed vectors per fold, which is the shape now refused by name).

* `MAPE` and `SMAPE` are documented as what they are: averages over the rows
  whose denominator is non-zero. Both have a denominator that can vanish --
  `MAPE` where the observation is zero, `SMAPE` where observation and
  prediction are both zero -- and each drops those rows rather than returning
  `Inf`, which is the right arithmetic but was reported nowhere. On a response
  taking exact zeros (counts, rainfall, abundance) the consequence is
  material: with 62 zeros out of 120, `MAPE` is an average over 58 rows
  presented as though it covered 120, and `SMAPE` drops precisely the rows a
  well-fitted model got right, so it reads worse than the fit deserves. The
  new "Percentage errors on responses with zeros" section on
  `model_metrics()` -- inherited by `evaluate_insample()`, `compare_models()`,
  `compare_models_cv()`, `summary()` and all four `cv_*()` -- says so, notes
  that the `n` column is the finite-pair count and not the row count either
  percentage error used, and points at RMSE/MAE/R-squared (and, for a Bayesian
  fit, CRPS and interval coverage) as the metrics unaffected by it. No
  computed value changes; returning the per-metric row count would alter the
  metric frame's column set and is deferred.

## Bug fixes

* Data carrying **no CRS** works again throughout. `ensure_projected()` now
  rejects a `target_crs` that does not resolve to a usable CRS (previously a
  typo silently made the call a no-op), but internal callers derive that target
  from another object — `st_crs(training_data)` — and that object is allowed to
  have no CRS. Passing `NA_crs_` through turned every CRS-less workflow into a
  hard error: `predict()` on all three backends, `make_folds(method = "nndm")`,
  `make_folds(boundary = )`, `prep_model_data(boundary = )` and
  `predict_surface()`. `cv_rf()` was worse than an error — the per-fold
  `predict()` threw, so every fold "failed" and `$overall` came back with
  `n_pred = 0` and `RMSE = NA` behind a generic warning. Internal call sites now
  pass `NULL` ("choose one automatically") when the source has no CRS; the
  user-facing validation is unchanged.

* `build_tessellation(crs = )`, `create_voronoi_polygons(crs = )` and
  `create_grid_polygons(crs = )` errored with sf's "cannot transform sfc object
  with missing crs" whenever the input had no CRS — exactly the users most
  likely to pass `crs =`. Reprojection is impossible there, but assumption is
  not: the target CRS is now stamped on with a loud warning, matching what
  `ensure_projected()` already documents. Input that *does* carry a CRS is still
  reprojected, not relabelled.

* `compare_models_cv()` built its argument list with `c(list(...), rf_args)`, so
  any `gwr_args`/`rf_args` entry whose name collided with one the function sets
  itself produced two entries of that name and `do.call()` died with "formal
  argument 'seed' matched by multiple actual arguments". Since `cv_rf()` has
  both `k` and `seed` as formals, `rf_args = list(seed = 3)` — straight from the
  documented usage — was enough to trigger it. Extras now replace base entries
  by name. `data_sf`, `response_var`, `predictor_vars` and `folds` are protected
  and dropped with a warning, because a per-model override of those would
  silently make the models incomparable.

* `compare_models()` given a single bare `spatial_fit` reported all-`NA`
  Moran's I columns and logged "'fit' is not a spatial_fit object" once per
  component. A `spatial_fit` is itself a list, so it passed the `is.list()`
  check and the loops then iterated the fit's own components as though they
  were models. It is now wrapped into a one-element named list, exactly as
  `evaluate_insample()` already did.

* `residual_morans_i()` failed on its own documented fast path. With `FNN` and
  `Matrix` installed the weights are a sparse `Matrix`, and `base::crossprod()`
  does not S4-dispatch on the `dgeMatrix` that `W %*% resid` produces, so the
  call died with "requires numeric/complex matrix/vector arguments" — taking
  `compare_models()`, which calls it automatically, down with it. Rewritten as
  `sum(resid_c * (W %*% resid_c))`, which is numerically identical and uses only
  dispatching primitives. `determine_optimal_levels()` carried the same bug.

* `residual_morans_i()` no longer errors on constant non-zero residuals: the
  degeneracy guard tested the raw sum of squares where Moran's I is a function
  of the *centred* residuals, so `VI` came out `NaN` and `if (VI > 0)` raised
  "missing value where TRUE/FALSE needed". A non-finite `VI` is handled too.

* `.build_knn_weights()`'s n > 5,000 guard tests for both `FNN` **and**
  `Matrix`. Keyed on `FNN` alone, an unbounded dense n × n allocation went
  through whenever `FNN` was present but `Matrix` was not.

* `assign_features_to_polygons()` drops columns of `features_sf` that would
  collide with the polygon ID column, with a logged warning. `sf::st_join()`
  suffixed them (`poly_id.x` / `poly_id.y`), which defeated the rename
  afterwards and left the result with **no rows** — reachable simply by
  re-assigning already-assigned points. A join that still fails to produce the
  ID column now errors and names the columns it did produce.

* `summarize_by_cell()` keeps the `"deff_applied"` attribute when `cells_sf` is
  supplied; `dplyr::left_join()` rebuilds attributes from its `x` template and
  dropped it. The per-cell vector is remapped onto the joined row order, `NA`
  for cells holding no observations.

* `summarize_by_cell()` joins on the native ID type when both sides agree.
  Coercing unconditionally made the returned ID type depend on an unrelated
  argument and turned integer IDs into `"1"`, `"10"`, `"2"`, …. A genuine class
  mismatch still coerces both to character and logs why.

* `summarize_by_cell()` coerces non-POINT geometry before computing a
  variogram-based design effect. `sf::st_coordinates()` returns one row per
  *vertex*, so a POLYGON or multi-vertex MULTIPOINT feature misaligned the
  coordinate matrix with the data and fed the wrong points into every cell.

* `make_folds(method = "buffered_loo")` errors when the buffer excludes so much
  of the data that no fold retains two training points. Those folds used to sail
  through and be dropped one at a time inside the CV loop, so the only symptom
  was a generic "all folds failed" warning at the very end.

* `make_folds(method = "block_kfold")` refuses a block size yielding a single
  block covering the whole extent — one fold with an empty training set,
  reported as a run that merely happened to score `NA`. An accepted
  autocorrelation range could trigger it: `estimate_sac_range()` rejects ranges
  above half the bounding-box *diagonal* while block construction needs half the
  *width*.

* `make_folds()` coerces MULTIPOINT geometry rather than merely accepting it,
  for the `st_coordinates()` reason above; every fold was misaligned silently.

* Cross-validation no longer renumbers folds. `.remap_folds()` dropped unusable
  folds from a list, shifting every later fold's index, so `fold_metrics$fold`
  and `predictions$fold` stopped lining up with `make_folds()$assignment$fold`.
  The original index is carried through. Folds left with fewer than two training
  rows are detected there and logged, instead of failing one at a time deeper
  in.

* `cv_spatial()` rejects a `fit_fn` whose `predict()` returns the wrong number
  of values. Both the metric computation and the prediction frame recycled
  silently, so two predictions against four test rows yielded a four-row frame
  with metrics computed against fabricated pairs.

* Cross-validation under `parallel = TRUE` reports a fold that died in a worker.
  `parallel::mclapply()` returns a `try-error` rather than `NULL`, which the
  `NULL` filter kept, and the failure surfaced as "subscript out of bounds".
  `conditionMessage()` has no method for a `try-error`, so the diagnostic branch
  itself threw; the condition is now taken from the object's attribute.

* Cross-validation under `parallel = TRUE` is reproducible from `seed` and gives
  results identical to `parallel = FALSE`. `.cv_run_folds()` called
  `parallel::mclapply()` without seeding the fork streams, so each worker seeded
  itself from the clock and process ID. One seed per fold is now drawn in the
  parent, making each fold's stream a function of `(seed, fold index)` alone.

* `estimate_sac_range()` rejects a singular variogram fit.
  `gstat::fit.variogram()` signals failure by setting `attr(., "singular")` and
  returning normally, so testing only for a `try-error` made the spherical
  fallback unreachable *and* let a singular fit's `range` flow out as the
  estimated autocorrelation range — which `make_folds(auto_range = TRUE)` then
  sizes spatial blocks from.

* `.extract_gwr_values()` requires **every** model-matrix column to match a
  column of GWmodel's `SDF` before multiplying the local coefficients through. A
  partial match reconstructed a linear predictor missing one or more terms and
  returned it as the fitted value — plausible numbers that were simply wrong,
  feeding `fitted()`, `residuals()`, `summary()` and every metric with no
  warning. A non-numeric coefficient column is refused rather than coerced.

* `fit_gwr_model()` separates the three degenerate response cases. Folded
  together, an all-dropped dataset was reported as "binary (0 unique values)"
  and a constant response as "binary (1 unique value)", while a genuinely binary
  non-integer response (1.5 / 2.5) failed the integer-like gate and passed
  unremarked.

* `fit_gwr_model()` and `gwr_model_selection()` validate `bandwidth`.
  Unvalidated, `NA` gave "missing value where TRUE/FALSE needed", a length-2
  vector gave "the condition has length > 1", and with `adaptive = FALSE` a zero
  or negative distance reached GWmodel untouched.

* `fit_gwr_model()`'s local-collinearity spot-check no longer touches the RNG
  at all. It sampled its 30 locations from the global stream and fires only
  when `n > 30` with at least two numeric predictors, so the same script
  produced different fold assignments depending on how many predictors a model
  happened to carry; `cv_gwr()` calls it once per fold. The 30 locations are
  now evenly spaced ranks of the observations ordered by x, then y --
  reproducible, independent of the row order, and drawing no random numbers.

* `predict.gwr_fit()` returns an all-`NA` vector when every row of `newdata` is
  dropped as incomplete, matching the other two backends, rather than surfacing
  a raw sf-to-`Spatial` coercion error.

* `predict.bayesian_fit()` transforms `newdata` to the training CRS *before*
  cleaning it, and derives the surviving rows from one sentinel column instead
  of a second, separately-maintained copy of the cleaning rules. It errors when
  a predictor standardised at fit time is absent from `newdata` or has arrived
  as character — silently skipping it handed brms an unscaled column against a
  model fitted on a scaled one. Its failure path returns a matrix when
  `draws = TRUE`, honouring the documented return shape.

* `plot()` on a `spatial_fit` errors when there are no finite residuals, instead
  of producing a uniformly grey map from `limits = c(Inf, -Inf)`. A *perfect*
  fit is handled too: all-zero residuals gave `limits = c(0, 0)`, a degenerate
  diverging scale whose breaks collapse onto one value.

* `plot_tessellation_map()` logs a warning for a `fill_col` that is not present,
  instead of drawing an unfilled outline map with nothing to say anything had
  gone wrong — a mistyped `label_col` already warned. `xlim`/`ylim` are
  validated, and the `theme` default moved out of the formals so a Suggests
  package never appears in an exported function's default arguments.

* `harmonize_crs()` announces when it *stamps* a CRS rather than reprojecting.
  `sf::st_set_crs()` only relabels; the coordinates do not move.
  `ensure_projected()` already made that assumption loudly.

* `coerce_to_points()` rejects an EMPTY LINESTRING rather than misaligning the
  result. `st_line_sample()` yields no midpoint for one (and segfaults in sf
  1.0.x), so the sampled midpoints stopped corresponding 1:1 with the rows they
  are scattered back into. A count check backstops any other divergence.

* `evaluate_insample()` errors on an unnamed list. The loop is over
  `names(fits)`, so an unnamed list iterated zero times and returned `NULL`
  silently; `compare_models()` then died in `seq_len(nrow(...))` nowhere near
  the cause.

* `determine_optimal_levels()` coerces MULTIPOINT geometry rather than admitting
  it, and errors on a factor or character predictor by name instead of dying
  inside `colMeans()` with "'x' must be numeric".

* `create_grid_polygons()` passes both `cellsize` and `n` to
  `sf::st_make_grid()` when both are known; `st_make_grid()` does not ignore `n`
  in the presence of `cellsize` for square grids, and omitting it made sf
  recompute `nx = ceiling(w / cellsize)`, which floating-point division pushes
  one past the intended count. `n` is parsed and validated once, up front,
  instead of being silently coerced to `NULL` in one branch and erroring in the
  other.

* The grid cache key no longer truncates `target_cells`. `as.integer()` made
  25.2 and 25.7 collide on one key, so the second call silently received the
  first one's grid, and a `NULL` `target_cells` collapsed `paste0()` to
  `character(0)`, crashing the lookup.

* `build_tessellation()` normalises a CRS-less `points_sf` to `NULL` rather than
  `NA_crs_`, which is a list and so was not treated as "no CRS supplied"
  downstream. Hex and square grids are built in the points' CRS, so the grid and
  the points no longer end up in different CRSs and break the point-to-cell
  index.

* `build_tessellation(method = "triangles")` triangulates the **point set** when
  `geometry` is unavailable, via `sf::st_triangulate()` on the unioned points.
  The fallback previously triangulated the convex hull *polygon*, discarding
  every interior point. The result is still the Delaunay triangulation of the
  input; only the resolution of degenerate configurations can differ from
  qhull's, and the logged warning now says so.

* `ensure_stable_poly_id()` logs a warning naming the geometry types when it
  drops non-polygonal rows, which it silently did before.

* `voronoi_seeds_kmeans()` and `voronoi_seeds_random()` validate their inputs
  (`voronoi_seeds_random()` also accepts the `sfc` its documentation always
  promised), and clamping `k` to the number of distinct positions is logged.
  `get_voronoi_seeds(method = "provided")` logs a warning when `n` disagrees
  with `nrow(seeds)`, which it ignores.

* `gp_lengthscale_bounds()` validates `coords_xy` and `q_small`. A vector
  `coords_xy` failed inside `.safe_dist()` with "argument is of length zero" and
  an out-of-range `q_small` inside `quantile()`, neither naming the argument.

* `fit_bayesian_spatial_model()` validates the response before handing it to
  Stan, where nothing points back at the column, and validates `gp_k`, `gp_c`
  and `control`. The inverse-gamma prior is written with `%.10g` rather than
  `%.6f`: a small scale rounded to the literal `"0.000000"` and Stan rejected
  `inv_gamma(a, 0)` from deep inside the model block. Tightly clustered
  coordinates get there. The half-normal fallback's scale is guarded the same
  way.

* `compare_models_cv()` names, in a warning, any `gwr_args` entry it drops.
  `cv_gwr()` has no `...`, so entries meant for `fit_gwr_model()` alone (e.g.
  `longlat`) were discarded silently and simply had no effect.

* `.compute_reg_metrics()` errors on a `y_train_mean` that is neither a scalar
  baseline nor one value per observation, instead of recycling it against the
  filtered response and silently distorting R².

* **`create_grid_polygons()` no longer truncates the grid when `cellsize` and
  `n` are both supplied. This changes results.** `sf::st_make_grid()` does not
  ignore `n` when `cellsize` is given: for square grids it takes the cell
  dimensions from `cellsize` *and* the counts from `nx = n[1]`, `ny = n[2]`,
  anchored at the bounding-box corner. `cellsize = 25` with `n = 2` on a
  100 × 100 boundary therefore produced 4 cells covering 2,500 of 10,000 square
  units and silently left three quarters of the study area with no cells at
  all — and because `clip = TRUE` had nothing outside the boundary to discard,
  the result looked like an ordinary, complete grid. `cellsize` now wins, `n`
  is dropped with a logged warning naming what it would have done, and the same
  call returns 16 cells covering the whole boundary. `n` is still forwarded
  when the *package* derived `cellsize` from it or from `target_cells`, which
  is what the original code was written for: omitting it there lets sf
  recompute `ceiling(w / cellsize)` and floating-point division pushes the
  count one past the intended value.

* **`fit_gwr_model()` no longer refuses a continuous response that happens to
  take two values. This changes results: fits that used to error now run.** The
  guard rejected any response with exactly two distinct finite values as
  "binary" and pointed at `GWmodel::ggwr.basic(family = "binomial")`. Two
  distinct values is not the same thing as binary: a measurement censored at a
  detection limit or saturated at a ceiling (0.0031 / 12.7401) is perfectly
  continuous, Gaussian GWR on it is a well-defined least-squares problem, and
  the advice to switch to a binomial family is nonsense for such values. The
  hard stop is now gated on the response also being integer-like, which is what
  the surrounding code already used to separate coded categories from
  measurements. A two-valued non-integer response raises a `warning()` naming
  the two values and asking you to confirm it is genuinely continuous, then
  fits. This also mattered inside `cv_gwr()`, where the guard runs once per
  fold and a small training fold can legitimately hold only two distinct
  values.

* **`determine_optimal_levels()` no longer reports a Moran's I that is
  arithmetically fixed. This changes which cell counts it returns.**
  `.morans_i_for_k()` builds a `min(8, n_cells - 1)`-nearest-neighbour weight
  matrix, so at nine cells or fewer every cell neighbours every other one. The
  row-standardised matrix is then complete, `W %*% e = -e/(n - 1)` for *any*
  mean-zero residual vector, and Moran's I collapses to exactly
  `-1/(n_cells - 1)` whatever the data are. That is not merely uninformative:
  `|I| = 1/(n_cells - 1)` falls monotonically in the number of cells, so
  `criterion = "morans_i"` ranked the largest evaluated candidate first every
  time, and `"combined"` carried the same tilt at half weight. Candidates below
  the floor now return `NA_real_` and are excluded from the model-aware
  ranking; when none clears it — the usual outcome at the default
  `max_levels = 12`, since the search evaluates a window around the elbow — the
  call falls back to the geometric ranking and logs a warning. Raise
  `max_levels` above roughly 10 for the model-aware criteria to contribute at
  all. `predictor_vars` also accepts logical columns now, read as 0/1, matching
  `fit_rf_model()`/`cv_rf()`/`predict()`; factor and character predictors are
  still refused by name.

* `residual_morans_i(fit, k = 1)` no longer errors with "subscript out of
  bounds" on a machine without `FNN`. In the dense fallback the inner function
  returns a scalar at `k = 1`, so `apply()` simplified the neighbour table to a
  length-n vector and `t()` made it a 1 × n matrix; indexing `nn_idx[i, ]`
  then failed for every `i > 1`. The result is now forced to `n × k`.

* `make_folds()` no longer dies on an empty or non-finite geometry.
  `st_coordinates()` yields one all-`NA` row per EMPTY POINT rather than zero
  rows, so a row-count check let them through: `block_kfold`'s
  `st_intersects()` returned `integer(0)`, `..block_id` went `NA`, and the
  nearest-block rescue aborted with "replacement has length zero". Unusable
  rows are now dropped with a warning naming the count, after `..row_id` is
  stamped so the survivors keep their original row identities, and for every
  method rather than just `block_kfold` — `random_kfold` would otherwise put an
  unplottable point in a fold, and `nndm` and `buffered_loo` both feed the
  coordinates to distance code. The rescue itself uses `vapply()` rather than
  `apply()`, so a point whose distances are all `NA` keeps its `NA` instead of
  collapsing the assignment. `points_sf` with no usable coordinates at all is
  an error naming that, not a downstream one.

* Every cross-validation wrapper names the cause when folds fail.
  `.cv_run_folds()` returns each fold's error text rather than a bare `NULL`,
  and `cv_gwr()`, `cv_bayes()` and `cv_spatial()` append `First error: ...` to
  both the logged and the R-level "all N folds failed" message. Running
  `cv_bayes()` without `brms` installed previously produced five `fold N fit
  failed` warnings and an all-`NA` `$overall` with `n_pred = 0` in which the
  word "brms" never appeared.


* **The package's log lines no longer depend on the user's global `logger`
  configuration.** `logger` seeds a new namespace from the global one, so the
  `"spatialkit"` namespace inherited whatever formatter the user had set before
  loading — and every logging helper hands `logger` an *already formatted*
  string. Under a user's `formatter_sprintf`, every package message containing
  a literal `%` — the CRS distortion figures in `ensure_projected()`, the local
  collinearity percentage in `fit_gwr_model()` — hard-errored with "too few
  arguments", and because the helper logs *before* it raises the R warning,
  the warning the manual promises died with it. Under the default
  `formatter_glue` a `{...}` inside a fold error was re-evaluated. The
  namespace's formatter is now pinned to `formatter_paste`, so the message
  logged is the message written.

* `fit_bayesian_spatial_model(backend = "auto")` chooses **cmdstanr** only
  when a CmdStan build is actually available. It chose it whenever the
  **cmdstanr** *package* could be loaded — a thin interface that is often
  installed without the toolchain it drives — so on such a machine every fit
  died inside the sampler with "CmdStan path has not been set yet. See
  ?set_cmdstan_path". The package's own weekly `check-brms` job was one such
  machine: it installs **cmdstanr** to satisfy Suggests and never builds
  CmdStan, and every scheduled run since the Stan smoke tests landed failed
  there. "auto" now falls back to **rstan**, which **brms** always brings,
  and logs the choice; an explicit `backend = "cmdstanr"` with no usable
  build is an error that says to run `cmdstanr::install_cmdstan()`.

* `cv_bayes(seed = )` reaches the sampler. `fit_bayesian_spatial_model()`
  carries `seed = 123` and the per-fold `fit_args` never set it, so every fold
  of every run sampled from Stan seed 123 and changing `seed` changed nothing
  on fixed folds. Each fold now draws its own sampler seed from the fold's
  seeded stream, as `cv_rf()` does for the forest; a `seed` in `fit_args`
  still overrides it for every fold.

* `cv_*(seed = NULL)` is reproducible from `set.seed()` under `parallel > 1`,
  as the README promised without qualification. With `seed = NULL` no per-fold
  seeds were drawn, so each forked worker was seeded by `mclapply()` from the
  clock and the process ID — three runs after the same `set.seed(777)` gave
  0.5687, 0.5649 and 0.5623 while the sequential call was reproducible. The
  per-fold seeds are now drawn from the caller's current stream (advancing it,
  as any RNG-consuming call would), so the sequential and parallel paths are
  the same function of the state `set.seed()` left.

* Warnings raised inside a fold reach the caller from the parallel path. R
  conditions do not cross a fork, so under `parallel > 1` every warning raised
  by the model — including `fit_gwr_model()`'s documented integer-response
  warning, raised in every fold — reached nobody, while the numbers came back
  identical and the run looked like a clean version of the same analysis. The
  worker now collects them and the parent re-raises each distinct message once.

* `cv_spatial()` (and therefore `cv_rf()`) returns the same typed, zero-row
  `fold_metrics` frame as `cv_gwr()` and `cv_bayes()` when every fold failed,
  so `subset(fold_metrics, RMSE < 5)` works instead of erroring on a missing
  column.

* `residual_morans_i()` refuses `k` large enough to make the neighbour matrix
  dense. The only size guard (n > 5000) applied to the dense fallback; with
  **FNN** and **Matrix** present, `k >= n - 1` allocated `n (n − 1)` pairs
  unguarded. Requests above 2e7 pairs are now an error naming `k` and `n`.

## New features

* New `fit_rf_model()` and `cv_rf()`: a `ranger` random forest as a first-class
  backend, returning an `rf_fit` that works with `cv_spatial()`,
  `predict_surface()`, `area_of_applicability()` and `plot()` like any other
  model. Three defaults are opinionated: `include_coords = FALSE` (a forest
  given the coordinates memorises location and fails wherever it has not been —
  Meyer et al. 2019, <https://doi.org/10.1016/j.ecolmodel.2019.108815> — and random CV does
  not catch it); `fitted()` returns **out-of-bag** predictions, so `summary()`
  on an `rf_fit` is not comparable with the other backends and says so
  (`$info$fitted_are_oob`); and importance defaults to permutation rather than
  impurity, which is biased toward continuous and high-cardinality predictors
  (Strobl et al. 2007, <https://doi.org/10.1186/1471-2105-8-25>). Compare backends with
  `compare_models_cv()`, which now has an RF branch.

  `predict()` on an `rf_fit` refuses the type confusions `ranger` would
  otherwise absorb silently: a numeric-at-fit predictor supplied as text (which
  ranger factor-codes, then applies numeric split thresholds to the codes), a
  logical-at-fit predictor supplied as text, and a categorical level the forest
  was never *grown* with — the level set is not enough, since a spatial fold
  holding out a whole class leaves a level with no training rows. Arguments
  that make ranger return a matrix (`predict.all = TRUE`, `type = "quantiles"`)
  are rejected rather than flattened column-major. A constant seed is supplied
  to ranger's predict unless the caller passes one, so prediction does not
  consume the global RNG and `predict_surface(chunk_size = )` — a performance
  knob — cannot shift later random draws. `cv_rf(seed = )` reaches the forest
  in every fold, and gains `pointize`. Passing ranger's own spelling of an
  argument the wrapper already sets (`num.trees`, `min.node.size`,
  `num.threads`, `mtry`, `importance`, `seed`) through `...` is an error naming
  the wrapper argument to use, rather than reaching `ranger()` twice.
  `num_threads` defaults to `getOption("mc.cores", 1L)` — one thread unless
  the session has opted in — for both the fit and `predict()`, rather than
  ranger's own default of every core on the machine, and `cv_rf(parallel = )`
  runs each forked fold's forest on one thread unless told otherwise, so the
  worker count is never multiplied by a thread count. See `?fit_rf_model`.

* Every `cv_*()` and `compare_models_cv()` accept `folds` as a vector of fold
  labels, one per row — `make_folds()$assignment$fold`, the object most
  naturally to hand — in addition to a `make_folds()` result and a list of
  `train`/`test` splits, the three shapes `area_of_applicability()` already
  took. The label vector used to fail with R's "$ operator is invalid for
  atomic vectors".

* New `area_of_applicability()`, implementing the dissimilarity index of Meyer &
  Pebesma (2021, <https://doi.org/10.1111/2041-210X.13650>). Predictors are centred and
  scaled on the training data's own statistics, optionally weighted by variable
  importance — by the importance itself, not its square root, matching `CAST`.
  A prediction point's DI is its distance to the nearest training point in that
  space over the mean pairwise training distance, and the threshold is the
  outlier-removed maximum of the training data's own DI. Pass the
  `make_folds()` result you actually validated with — the area is defined
  relative to a performance estimate, and a blocked estimate is a claim about
  predicting further away.

  A model fitted with `include_coords = TRUE` is measured in coordinate space,
  since an index that ignores location would report a point far outside the
  training extent as *inside* on ordinary covariate values alone; weights for
  the two coordinate columns default to the mean of those supplied, as the
  caller has never seen them. Non-`POINT` `newdata` is reduced to points, and a
  CRS present on one side is applied to the other. The zero-variance test is
  relative to each column's magnitude rather than an absolute tolerance, so a
  predictor is not dropped for the **unit** it was recorded in. Categorical
  predictors are refused rather than dummy-coded; logicals are read as 0/1.
  A `make_folds()` result is resolved by its `..row_id` **values**, which
  coincide with row positions only when the input carried no prior IDs.
  See `?area_of_applicability`.

* New `select_features_forward()`: greedy forward feature selection with
  **spatially blocked inner folds**, which is the whole point of having it.
  Random inner folds inside blocked outer folds select variables that look
  predictive only because nearby points leak between train and test, and the
  outer loop then reports honest-looking numbers for a dishonestly chosen
  feature set. `method` defaults to `"block_kfold"` and logs a warning if set
  to `"random_kfold"`. The empty set is scored first where the backend can fit
  it, so the first variable is judged against a null-model baseline rather than
  accepted unconditionally, and `history` carries that baseline as a `step = 0`
  row. Every candidate set is scored on the same observations — the
  completeness filter matches `prep_model_data()` exactly, finiteness test
  included, so a candidate carrying a single `Inf` cannot be preferred for
  having an easier subset — and the inner folds are built once, before the
  sweep, rather than rebuilt per candidate. A `max_fits` budget guards against
  nesting a sweep inside leave-one-out outer folds. Where the backend cannot
  fit the empty set at all — `fit_rf_model()` and `fit_gwr_model()` both refuse
  a zero-length `predictor_vars` — the probe is silent on the console: its
  per-fold failures go to the file trace only, rather than printing the same
  lines a genuinely failed run prints.

* New `gwr_model_selection()`: wraps `GWmodel::gwr.model.selection()` (Lu et al.
  2014, <https://doi.org/10.1080/10095020.2014.917453>) and returns a ranked table instead
  of two loosely-coupled lists. It is the fast, in-sample counterpart to
  `select_features_forward()` — the same forward search scored by **AICc**,
  read from the documented `c(bandwidth, AIC, AICc, RSS)` layout of GWmodel's
  `GWR.df`, which carries no column names; the result records whether the table
  arrived in that shape. Candidates must be numeric, and `dmat_max_n = Inf`
  means *always precompute* the distance matrix. Both limitations are
  documented rather than papered over: one bandwidth is shared by every
  candidate (which is what makes the criteria comparable), and the null model
  is never evaluated, so the result always names at least one predictor. When
  it disagrees with the blocked estimate, believe the blocked estimate.
  See `?gwr_model_selection`.

* New `predict_surface()`: builds a regular grid over the training extent (or a
  grid you supply), joins covariates, predicts in chunks and returns `sf`.
  Supports `boundary` clipping, `cell_size` or approximate `n_cells`, and
  `se = TRUE` for a posterior-SD surface where the backend exposes draws.

* New `plot()` method for `spatial_fit`, with `type = "residuals"`,
  `"observed_predicted"` and `"variogram"` (the empirical residual variogram
  with the fitted model and effective range overlaid, so the fit can be judged
  rather than trusted). The variogram's distance axis is labelled in the units
  of the CRS it was actually fitted in — metres of an auto-chosen zone for a
  lon/lat fit, not the caller's degrees — it names the azimuth when a single
  direction is drawn, and a fit that did not converge says so in the caption.
  New `plot_folds()` maps a fold scheme, which is the fastest way to see
  whether spatial blocks separate the data or are smaller than the
  autocorrelation range and therefore leaking.

* `make_folds()` gains `method = "leave_location_out"`, which keeps every
  observation from a location (named by the new `group_var`) in the same fold.
  Repeated measurements at one site were previously unrepresentable: random
  k-fold splits them across folds, so the model is scored partly on sites it
  trained on.

* `make_folds()` gains `method = "nndm"`, implementing the distance-matching
  principle of Milà et al. (2022, <https://doi.org/10.1111/2041-210X.13851>), as in
  `CAST::nndm()`. Rather than choosing a `buffer` with nothing to justify it,
  the exclusion around each held-out point is sized so the training-to-test
  distance distribution reproduces the distances from your actual prediction
  locations (the new `prediction_points`) to the training data. The procedure
  follows the paper's iterative exclusion removal for removal and is
  deterministic: no random numbers are drawn, so the caller's RNG is untouched,
  and ties in the nearest-neighbour distance — every mutual-nearest-neighbour
  pair, all of a regular grid — are broken by the point's position rather than
  by its row index, so identical data give identical folds whatever order the
  rows arrive in (`CAST` breaks them by row).
  `params$target_median`, `params$realised_median` and
  `params$max_ecdf_excess` record how close the match came, and `min_train`
  (default 0.5) and `phi` control it. Matching is as close as the training
  configuration permits — the achievable distances are discrete order
  statistics. When prediction locations sit no further from the training data
  than training points sit from each other, plain leave-one-out already
  reproduces the target and nothing is excluded; that is the correct outcome.
  A non-POINT `prediction_points` layer is reduced to points first, since
  point-to-polygon distances are zero for any cell containing a training point
  and would collapse the scheme towards plain LOO.

* `summarize_by_cell()` gains `deff = "variogram"`, computing a per-cell design
  effect from a fitted variogram rather than one pooled intra-class correlation.
  For `n` points in a cell with correlation matrix `R` the effective sample size
  of the mean is `n^2 / sum(R)`, so `deff = sum(R) / n`. This generalises the
  Kish option — a constant off-diagonal correlation recovers
  `1 + (n - 1) * rho` exactly — but lets correlation decay with distance, which
  is what having fitted a variogram is for. Pass the fit via the new `sac`
  argument, or it is estimated when `response_var` is supplied — on the
  **response**, not on OLS residuals, even when `predictor_vars` are listed: the
  `..se_resp_*` columns estimate the SE of the cell mean as an estimate of the
  response's grand mean, so the correlation to correct for is the response's
  own (measured grand-mean coverage 0.93 with the response variogram against
  0.51 with the residual one on a field with a smooth predictor). Pass a
  residual variogram through `sac` if that is the field you want. Large cells
  are subsampled at `deff_max_n` (default 500), with the correlation scaled
  back to the cell's own size. A `sac_range` whose fit was *rejected* carries no
  usable correlation function, so both the supplied and the internally estimated
  path fall back to `deff = 1` and say so rather than saturating the correlation
  at every within-cell distance. One correlation function is fitted and applied
  to every numeric column, response and predictors alike, because a variogram is
  a property of the field rather than of a variable type.

* `fit_bayesian_spatial_model()` supports intercept-only models
  (`predictor_vars = character(0)`): the response is explained by the intercept
  and the spatial GP alone, the natural null for asking how much of a surface is
  spatial structure rather than covariate effect.

* `fit_bayesian_spatial_model()` checks the posterior length-scale against the
  smallest scale the chosen basis can resolve and logs a warning when more than
  10% of the posterior mass falls below it — the adequacy diagnostic recommended
  by Riutort-Mayol et al. (2023, <https://doi.org/10.1007/s11222-022-10167-2>), and what
  makes the smaller default `gp_k` safe rather than merely cheaper. `$info`
  gains `gp_c`, `gp_n_basis`, `gp_ell_min` and `gp_lengthscale_bounds`, and
  `print()` on a `bayesian_fit` and `cv_bayes()`'s `fold_metrics` report the
  total basis count alongside the per-dimension rank.

* `cv_spatial()` raises a condition when folds fail, matching `cv_gwr()` and
  `cv_bayes()`; an all-failing `fit_fn` previously returned an all-`NA`
  `overall` and an empty `fold_metrics` with nothing at R condition level. The
  result records `n_folds_attempted` and `n_folds_succeeded` — compare them
  before trusting `overall`.

* `make_folds()` records the CRS the folds were built in as `params$crs`
  (`"EPSG:32632"`, an input string, or a WKT). `block_size` and `sac_range` are
  lengths in *that* CRS, which is not necessarily the one the caller passed:
  geographic input is projected by `ensure_projected()` to a CRS chosen for the
  extent. Without the label the units of a recorded block size were not
  recoverable from the result.

* `spatialkit_quiet()` is a new exported helper. Both
  `logger::log_appender()` and `logger::log_threshold()` default to `index = 1`,
  which is the temp-file trace, so the two-line recipe in the README could not
  redirect or quieten the **console** echo (index 2) — there was no documented
  way to silence the package. The README now says so too.

## Documentation

* `estimate_sac_range()` documents its three return shapes (a range, a rejected
  range, and no fit at all) and which attributes each carries.

* `make_folds()` documents that `k` is not always honoured: `buffered_loo` and
  `nndm` always return `k = n`, and `block_kfold` and `leave_location_out` lower
  it when the geometry or the grouping cannot support the request. Read
  `folds$k`.

* `new_spatial_fit()` documents the `coef()` contract; `summary.spatial_fit()`
  and `model_metrics()` document that their metrics are in-sample for a
  `gwr_fit` and a `bayesian_fit` but out-of-bag for an `rf_fit`;
  `prep_model_data()` documents that the projected CRS is not an unconditional
  guarantee, since `ensure_projected()` passes a CRS-less dataset through
  unchanged when its coordinates do not look like lon/lat.

* The vignette and `inst/scripts/example_nc_demo.R` read fit quality from
  `fit$metrics$r_squared` and CV results from `cv$summary$rmse`. Neither field
  has ever existed. Because `sprintf()` returns `character(0)` when any argument
  has length zero, the reporting lines printed *nothing* rather than erroring,
  so the shipped vignette silently omitted every number it claimed to show. Both
  now use `model_metrics()` and `$overall`.

* The demo's Voronoi tessellation was built from all 300 observations rather
  than from the 40 k-means seeds it computed one line earlier — one cell per
  observation, a nearest-neighbour interpolation rather than an aggregation,
  compared side by side against two ~50-cell grids. The seeds are now used.

* The vignette builds as `rmarkdown::html_vignette` rather than
  `html_document`, guards its `ggplot2` and `geometry` use, demonstrates
  `summarize_by_cell()` instead of reimplementing it with
  `group_by()`/`summarise()`, and adds a spatial cross-validation section
  contrasting `block_kfold` against `random_kfold` on the same data.

* The package-level help page (`?spatialkit`) gains "The pipeline, in order"
  and "Where to start" sections, so `help(package = "spatialkit")` leads
  somewhere rather than presenting 40 exports in alphabetical order.

* Every exported function's description now says *when to reach for it* rather
  than only what it does, and `@family` / `@seealso` links connect each step of
  the pipeline to the one before and after it — `assign_features_to_polygons()`
  to `summarize_by_cell()`, `determine_optimal_levels()` to
  `build_tessellation()`, `new_spatial_fit()` to `cv_spatial()`, and the two
  seeding functions to each other. `create_voronoi_polygons()` versus
  `create_grid_polygons()`, and `voronoi_seeds_kmeans()` versus
  `voronoi_seeds_random()`, each say which to pick and why.

* `build_tessellation()` documents that `boundary` is **required** for
  `method = "hex"` and `method = "square"` — the grid methods have no extent of
  their own — and optional for `"voronoi"` and `"triangles"`, which derive one
  from the points. The error existed; the requirement was not written down
  anywhere.

* `create_grid_polygons()` documents that `target_cells`, `cellsize` and `n`
  are three ways of sizing one grid and that exactly one should be supplied,
  that `cellsize` is in the units of the working CRS, and that `cellsize` takes
  precedence over `n`.

* `determine_optimal_levels()` documents the nine-cell resolution floor on the
  model-aware criteria, why it exists, and that the whole call falls back to
  the geometric ranking when no candidate clears it.

* `compare_models_cv()` documents that dropping every requested backend is an
  error (`"no viable models."`) rather than an empty comparison, and that the
  returned frame carries only the models that actually ran, so callers should
  check which names are present rather than assuming one row per request.

* `new_spatial_fit()` documents the two obligations on a custom backend:
  return an object built by the constructor, and define a
  `predict.<subclass>()` method — `cv_spatial()` scores folds through the
  `predict()` generic, so without one every fold fails.

* **README.** A new "Your own data" section shows both entry points —
  `st_read()` for a spatial file and `read.csv()` + `st_as_sf()` for a table of
  coordinates — using the `nc.shp` demo shapefile shipped with `sf` so it runs
  anywhere. The README previously manufactured every example inline with a
  hard-coded `crs = 32632` and never showed data entering the package at all.
  A companion "CRS: what the numbers are in" subsection states that block
  sizes, buffers, bandwidths, variogram ranges and `expand` distances are in
  the units of the working CRS; that geographic input is projected
  automatically to a CRS chosen for the extent; and how to pin one.

* **README.** New guidance where none existed: how to choose among the four
  tessellation methods, how `k` and `block_size` trade off against the
  autocorrelation range, what to do when `estimate_sac_range()` returns `NA`,
  how to read a design effect, which model backend to reach for (with the
  recorded cost of each), and a "Troubleshooting" section covering the errors a
  new user actually hits first. A worked hex-grid example replaces the previous
  picture-only coverage of the grid methods.

* **README.** Three corrections. The `estimate_sac_range()` example showed a
  rejected range printing its attributes, which `print.sac_range()` has not
  done since the attribute dump was removed; it now shows the bare `NA` and
  reads the attributes explicitly. The `determine_optimal_levels()` passage
  claimed the residual-autocorrelation criterion was doing work at cell counts
  where it is arithmetically degenerate. The test-suite paragraph said "exactly
  one" test guards on `brms`; six do, five of them additionally gated behind
  `SPATIALKIT_TEST_BRMS` so they never run in the matrix.

* **`inst/scripts/example_nc_demo.R`** said EPSG:2264 was projected "so
  distances are metric". Its unit is the US survey foot, which is what the
  script's own "Autocorrelation range: %.0f ft" line reports. The comment now
  says planar, and names the unit every distance, bandwidth and block size in
  the script is in.

* **Vignette.** `print(rf_fit)` and `summary(rf_fit)` report the same OOB RMSE
  but different R² (0.4733 against 0.4715). The vignette now explains why:
  `print.rf_fit()` echoes `ranger`'s `r.squared` (`1 - MSE/var(y)`, unbiased
  n − 1 variance) while `summary()` recomputes `1 - SS_res/SS_tot` from the
  same out-of-bag predictions with an n denominator, so the unexplained
  fractions differ by exactly n/(n − 1).

* Every exported function has runnable examples: the eleven that shipped
  without any — `clear_fitted_cache()`, `clear_grid_cache()`,
  `clip_target_for()`, `compare_models()`, `create_grid_polygons_cached()`,
  `ensure_stable_poly_id()`, `evaluate_insample()`, `harmonize_crs()`,
  `model_metrics()`, `voronoi_seeds_kmeans()` and `voronoi_seeds_random()` —
  gained one, and the two `\dontrun{}` blocks say why they cannot be run (a
  Stan toolchain and minutes of MCMC).

* `residual_morans_i()`'s default null is described correctly: `null =
  "auto"` uses the Cliff & Ord *regression-residual* moments whenever the
  residuals are OLS residuals on the rebuilt design, and the randomisation null
  otherwise; the README said "the randomisation variance" without
  qualification. The type-I error of a random forest's residual test is
  attributed to what the package actually feeds it — out-of-bag residuals,
  which are honest out-of-sample errors with their own spatial structure — not
  to "shrunk in-sample residuals". `determine_optimal_levels()` gives the real
  reason for its nine-cell floor (the standardised deviate is 0/0 there, so the
  criterion carries no information) rather than an argument from `|I|` that its
  own details section had just called wrong. `summarize_by_cell()` notes that
  the "use the naive SE for the cell's own mean" advice is calibrated under
  uniform within-cell sampling. `cv_gwr(bandwidth = )` states its units and
  semantics like its siblings. `quiet` is documented as "suppress this
  function's progress messages" everywhere, with a pointer to
  `spatialkit_quiet()` for the console log echo it does not touch.

* README: the square-grid call returns 36 cells, not 32; the installation
  section no longer promises a specific version from CRAN; the resolution
  figure and the quick-start output are regenerated for the elbow-first
  ordering (`4 3 5`, `k = 4`).

* The DESCRIPTION now cites the methods it implements -- Lu et al. (2014) for
  `GWmodel`, Riutort-Mayol et al. (2023) for the Hilbert space Gaussian
  process, Strobl et al. (2007) for the permutation importance, Mila et al.
  (2022) for NNDM folds and Meyer and Pebesma (2021) for the area of
  applicability -- each with its DOI, and quotes only software names. The
  Riutort-Mayol reference on `fit_bayesian_spatial_model()`'s help page gives
  the article number (33, 17) rather than "33, 1".

* No example is wrapped in `\donttest{}` any more: the fifteen that were --
  every fit, cross-validation, comparison, plotting and surface example that
  needs a Suggests package -- run unconditionally behind their
  `requireNamespace()` guards, the slowest in under 2 s. The two Stan examples
  (`fit_bayesian_spatial_model()`, `cv_bayes()`) keep `\dontrun{}` because they
  need a C++ toolchain and minutes of MCMC, and say so in a leading comment.

# spatialkit 1.0.0

First CRAN release, published 2026-08-07.

* CRS management: `ensure_projected()`, `harmonize_crs()`,
  `coerce_to_points()`, `prep_model_data()`.
* Voronoi, hexagonal, square and Delaunay tessellation
  (`build_tessellation()` and the `create_*_polygons()` functions), with
  boundary clipping, stable reproducible cell IDs (`ensure_stable_poly_id()`)
  and a memoised grid builder (`create_grid_polygons_cached()`).
* Seeding (`get_voronoi_seeds()`, `voronoi_seeds_kmeans()`,
  `voronoi_seeds_random()`) and resolution selection
  (`determine_optimal_levels()`).
* Feature-to-polygon assignment (`assign_features_to_polygons()`) and
  cell-level aggregation with design-effect-corrected standard errors
  (`summarize_by_cell()`).
* GWR (`fit_gwr_model()`) and Bayesian spatial Gaussian process
  (`fit_bayesian_spatial_model()`) backends behind a common `spatial_fit` S3
  class.
* Spatial cross-validation: `make_folds()` with `random_kfold`, `block_kfold`
  and `buffered_loo`; `estimate_sac_range()`; `cv_gwr()`, `cv_bayes()` and
  `cv_spatial()`.
* Model comparison and diagnostics: `compare_models()`, `compare_models_cv()`,
  `evaluate_insample()`, `residual_morans_i()`.
* Tessellation mapping (`plot_tessellation_map()`) and scoped logging.
