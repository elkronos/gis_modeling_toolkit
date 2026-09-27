# spatialkit (development version)

## New features

* `fold_separation()` measures what blocked cross-validation is for: the
  distance from every held-out point to its nearest training point,
  summarised per fold.  Until now the only evidence a fold scheme worked was
  a proxy --- the block size compared against the estimated autocorrelation
  range, which `make_folds()` warns about.  That is a statement about the
  design; this is the result, and the two can disagree in the direction that
  matters.  Measured on 300 points over a 1000-unit square with a fitted
  range of 292.5: blocks of 343 units are wider than the range, so
  `make_folds()` raises no warning at all, and yet 76 percent of the
  held-out points still sit closer to a training point than the range, the
  nearest of them 23.3 units away.  Blocks wider than the range leak wherever
  a test point sits near a block edge with training data just across it,
  which in a fine grid is most of them.  The returned table carries
  `n_train`, `n_test`, `n_blocks`, `min_dist`, `median_dist` and
  `within_range` (the share inside the range), and its `print()` method
  closes with what that share means for the score.  Nothing is estimated:
  the distances come from the geometry and the range is the one the folds
  already carry or the one you pass.  A range that records its CRS (the
  folds' own `params$sac_range` and `params$crs`, or an
  `estimate_sac_range()` result) is compared with distances measured in that
  CRS, so a copy of the layer in other units gives the same answer (a
  US-foot copy used to report 19 percent of the hold-out inside a metre
  range where the metre layer reported 89); `sac` and the distances are in
  the units of the CRS the `crs` attribute names, metres for lon/lat input.
  `fold` is the `fold_id` a `cv_*()` result's `$folds` carries, so a join on
  it pairs each fold's error with its own distances (it was the list
  position, which after a dropped fold paired one fold's error with
  another's distances).  A `units` object for `sac` is refused.  As in
  `cv_*()`, a `make_folds()` result whose recorded rows sit at other
  locations in `data_sf` is refused.  Folds built on `prep_model_data()`'s
  200 points and measured on the 195 that `assign_features_to_polygons()`
  kept used to be matched by position, and reported 100 percent of the
  hold-out inside the range against 23--48 percent on the right layer.  The
  location check is skipped when one layer is POINT and the other is not,
  so folds built on polygons still measure their pointized copy.  The
  `print()` verdict now depends on the fold scheme: "widen the blocks" is
  said only of blocked folds, buffered leave-one-out is told to widen the
  buffer, and random, leave-location-out and hand-made splits are told to
  use blocked or buffered folds.  NNDM folds are no longer called
  optimistic: they are built to reproduce the prediction-to-data distances,
  so the share describes the prediction task.  The range line gives the
  CRS's unit ("in metres"), where it used to print the CRS code as if it
  were a unit ("in EPSG:32617 units").

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
  subsetted layer), `n_unknown_ids` counts fold entries naming rows the data
  does not have, and `n_dropped` the rows `prep_model_data()` removed before
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
  states the parameter that decides whether the scheme leaks (the block
  size, the buffer, the number of location groups or the median NNDM
  exclusion), with the units of the folds' CRS; folds built from points
  without a CRS get no units rather than "(NA units)".

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
  `deff_applied` column, `TRUE` on every row when the correction was applied
  and `FALSE` when it fell back to the uncorrected standard errors; unlike
  the `"deff_applied"` attribute it survives `rbind()` and
  `dplyr::bind_rows()` of many results.  The default frame is unchanged.
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
  error.  The default `conf_level = NULL` returns exactly the frame it always
  did.  The degrees of freedom are `n - 1` at `deff = 1`, for a numeric
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
  scored on every criterion at once.  `resolution_profile()` runs a
  log-spaced ladder of level counts from a floor the autocorrelation range
  implies (`ceiling(area / range^2)`) to a ceiling the support implies
  (`floor(n / min_cell_n)`, where `n` is every point in the layer, not the
  `sample_n` subsample the k-means runs on), fits each level as the best of
  25 k-means++ restarts, and returns a data.frame with the WSS elbow
  statistic (read on log-log axes, and `NA` when the points have no cluster
  structure to bend the curve), Mallows'
  C_p of the piecewise-constant approximation of the response (or of
  its residuals on the predictors, fitted on the rows with a complete
  response and predictors, with the nugget from the fitted variogram as the
  noise variance and the penalty set for the whole layer), the
  standardised residual Moran's
  z of the cell means, and an analytic reliability of the cell means
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
  `range_floor = FALSE` starts the ladder at 2 whatever
  the range and only reports the floor, so profiles whose range estimates
  differ (cross-validation folds) have comparable ladders; with the floor
  applying on one fold and not the next, reliability and the elbow moved by
  a factor of 10--15 between folds.  The default, `range_floor = TRUE`,
  keeps the floor, and the help page describes that regime switch.  Under
  `select_on = "split"` a supplied `sac` is flagged in the log, since it
  must be fitted on the selection half, and the help page shows the
  two-call workflow.

* `select_on = c("all", "split")` on `determine_optimal_levels()`,
  `resolution_profile()` and `select_features_forward()`.  Whenever a
  selection reads the response --- a level count chosen with
  `response_var` and `predictor_vars`, or a predictor set chosen by a
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
  and are noted on the help page, not implemented.  The help page also now
  says plainly that supplying both `response_var` and `predictor_vars`
  upgrades the level-selection criterion to `"combined"`, so the selection
  depends on the response without that having been asked for.

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
  `resolution_profile()` (read with `select_resolution()` at its default
  criterion).  Both were hard errors before, so nothing that worked
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
  sizes, with random folds as the leaky reference, returned as a table
  with the fold-to-fold spread at each size and the estimated
  autocorrelation range alongside; `plot()` draws the curve with the range
  marked.  Blocks smaller than the range leak, so the curve rises from the
  random-fold value towards the range and plateaus beyond it, and the
  height of the rise is what the random-fold number overstated --- measured
  on simulated fields with a 100--130-unit range and a random forest with
  coordinates, the plateau begins at one to two times the estimated range.
  Each size is a full `cv_spatial()`, so a fit budget (`max_fits`, default
  60) refuses to start rather than run past it.  The default ladder runs up
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
  which the `k` column shows.  The default ladder also skips sizes whose
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
  collinear (`cn_slopes` above 30, or `cn` above 1e6), so a predictor's
  origin (degrees C or kelvin) does not change the verdict.  The warning is
  now the exact fraction of locations rather than a sampled one, and its
  thresholds (a quarter of the locations, or any) are unchanged; above a
  quarter it now says that an exactly singular window stops the fit, instead
  of promising non-finite coefficients.  The weighted survey sees what the
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
  range), whether it exceeds the design-based `s^2/n` of the
  plain mean (`kr_exceeds_design`), and the kriged-minus-plain shift in
  standard errors (`kr_shift`); plus the variance of the standardised
  errors from blocked cross-validation (`attr(, "cv")$zscore_var`), which
  is about 1 when the kriging variance is right.  Measured on simulated
  exponential fields: 0.93--1.07 with the true variogram, 0.85--1.01 with
  the estimated one; under blocked folds it checks the sill and range
  rather than the nugget (0.95--1.24 with the nugget understated tenfold),
  and random folds are the instrument for the nugget.  The comparison the
  function exists for: under uniform sampling kriged and plain means
  differed by more than one standard error in 11--24 percent of cells; under
  clustered sampling in 34--63 percent, with 3--27 of 16--64 cells empty
  and kriged anyway.  Repeat visits to one location are kriged from their
  mean, with the part of the nugget that varies between visits divided by
  their count; kriging them as separate rows made every kriging system
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

* `voronoi_seeds_kmeans()` gains `nstart`, the number of `stats::kmeans()`
  starts (default 10, as before).

* `cv_bayes()` and `compare_models()` say when a Bayesian fit did not
  converge.  Both scored such a fit like any other, and only the fit's WARN
  log lines (which `tryCatch()` and knitr never see, and
  `spatialkit_quiet()` hides) said that its posterior was not to be
  trusted.  `cv_bayes()`'s `fold_metrics` gains `convergence_ok`: `TRUE` or
  `FALSE` as `fit_bayesian_spatial_model()` judged that fold's sampler
  (R-hat, effective sample size, divergences), and `NA` when `fit_args`
  sets `check_convergence = FALSE`.  A run with any `FALSE` raises one
  warning naming those folds.  `compare_models()` gains the same column and
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
  no elbow the function still returns its linear-axis answer, but warns that
  the ladder chose it; `build_tessellation()` and `get_voronoi_seeds()` refuse
  a geometry-only `resolution_profile()` with no elbow instead of drawing a
  count from it.  A ladder of two levels (`max_levels` of 1 or 2, or three
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
  `prep_model_data()` and again where the data are handed to GWmodel.

* **`predict.gwr_fit()` lost every prediction to one bad location.**  It
  went through `GWmodel::gwr.predict()`, which returns nothing for any row if
  one location's window is empty or singular, so fixed-bandwidth block CV
  failed every fold.  It also never assigns its distance matrix once
  training and new rows together exceed 10,000 (a 100 by 100 prediction grid
  came back all `NA`), and it built the full training hat matrix for a
  variance it then discarded, which is cubic in the training size.
  Predictions now come from `gwr.basic(regression.points = )` in chunks, a
  failing chunk is redone location by location, and only a truly singular
  location is `NA`, with a warning counting them.  Where the old path
  worked, the values are identical.

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
  criterion a typo kept `"geometric"` where supplying both variables
  upgrades to `"combined"` (on one 400-point layer `7 6 8` instead of
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
  read", with the fallback warning.

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
  the largest lag fitted"`, the directional ranges still attached); the
  directional maximum stands in only for an all-pairs fit that is singular
  or did not converge.  On stationary fields whose range is close to the
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

* **GWR said a singular window came back as `NaN` coefficients; it stops
  the fit.**  GWmodel's matrix inverse throws on an exactly singular window,
  so a 0/1 indicator constant within clusters failed the whole fit with a
  bare "inv(): matrix is singular", while the help page and the collinearity
  warning promised masked `NaN` coefficients.  The non-finite coefficients
  GWmodel does return come from co-located points, where an adaptive
  bandwidth no larger than the number of observations at a site gives the
  kernel zero width (160 of 160 at 40 sites of 4 observations, 4
  neighbours), and the warning blamed singular windows for those.  The fit
  error now says a window is singular and how many the collinearity check
  found, and the non-finite warning names co-located points when they are
  the cause.

* **An adaptive GWR bandwidth above the number of points was capped in
  silence.**  `bandwidth = 1500`, meant as metres with `adaptive` left at
  `TRUE`, became a 200-neighbour, near-global fit on 200 points without a
  word; its local slopes varied less than half as much as the intended
  fixed-distance fit's.  `fit_gwr_model()` and `gwr_model_selection()` now
  warn, naming n and pointing to `adaptive = FALSE`.  `cv_gwr()` repeats the
  warning in each fold whose training set is smaller than the bandwidth.

* **Below 20 points, `bandwidth = NULL` did not fit at the bandwidth
  `bw.gwr()` chose.**  GWmodel searches adaptive bandwidths from 20
  neighbours up to n, so with fewer points its choice exceeds n (18 for 12
  points) and was capped at n, a different kernel with a worse AICc (18.3
  against 9.7), without a word.  The cap stays and now raises a warning in
  `fit_gwr_model()` and `gwr_model_selection()`; supply `bandwidth` for data
  this small.

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
  The rule is unchanged; the result is now logged with the remedy
  (leave-location-out folds, or deduplication) and `print()` shows the count
  of zero training DI.

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
  geometric axis is now the log-log sag the elbow is read from, and it is
  flat when the curve has no elbow, so Moran's z alone orders the window.
  Every unscored candidate takes the last place on the Moran axis, and exact
  ties go to the `k` nearest the elbow.  And an elbow below ten cells, a
  count Moran's I cannot score, is no longer ranked against the counts it
  can: ranking it put ten, the smallest count Moran's I scores, first
  whatever the response did (a response of noise and one varying by
  cluster gave the same answer).  There `"combined"` returns the geometric
  ranking, logs why, and records it in the diagnostics (`criterion =
  "geometric"`, `fallback`).  The same layers now give `8 7 9`, `8 7 9`,
  `7 6 8` and `8 7 9`, the geometric answer, and 800 uniform points, which
  have no elbow and are ordered by Moran's z, `10 11 7` where they gave
  `5 4 10`.  Supplying both `response_var` and `predictor_vars` selects
  this criterion by default.

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
  raised once per call, only when the standard errors really are the
  uncorrected ones, naming every reason, the rejected `sac` included.

* **`summarize_by_cell(deff = "variogram")` ignored a `sac` that carried no
  variogram model without saying so.**  A plain number, or
  `units::set_units(1.5, "km")`, was passed over and the design effect came
  from a variogram estimated from `response_var`, with no condition; the
  caller could not tell that the value given had not been used.  Such a
  `sac` is now set aside with a warning, as a rejected `sac` is: a plain
  warning when a variogram is estimated instead, and the classed
  `spatialkit_deff_fallback` warning, naming it, when `deff` falls back
  to 1.

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
  `cn_slopes` is above 30 or singular, or when `cn` is above 1e6, where
  GWmodel's uncentred solve starts losing precision in the slopes.
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

## Documentation

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
  naming it, and scripts 02, 09 and 10 skip themselves when gstat, GWmodel or
  brms is absent.
* Four vignettes join `spatialkit_nc_demo`: `getting-started` (installing, what
  the coordinates are in, the pipeline from points to a scored model and a map,
  and a glossary of the terms that recur), `resolution` (the ladder, the
  four criteria and why they disagree), `spatial-cross-validation` (the five
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

* `?estimate_sac_range` now says that the exponential model is kept
  whenever it converges and so overestimates the range on smoother fields
  (about 1.8--2.1 times a Gaussian practical range, 1.3--1.4 times a
  spherical one), and why the family is not chosen by fit (that biases
  exponential fields low, to 0.82); that 15 lag bins make a range spanning
  one or two of them come out long (60 m returned 89--102 m); that nothing
  tests for spatial structure (white noise gave a finite range in 8 of 30
  draws); that the "decreases with distance" refusal also fires on small
  samples of ordinary fields (7--9 of 60 at n = 30, 0--1 at n = 100), which
  its warning now says below 100 points; that the anisotropy note goes to
  the log file at INFO, not the console, and that the longest directional
  range is read from `directional_fitted` after checking
  `directional_status`, not `max(attr(, "directional"))`, which is `NA`
  exactly when the major axis ran past the fitted lags (the note itself, the
  nc_demo vignette and the help now say so); that an accepted range can
  still exceed half the width of the layer, leaving
  `make_folds(auto_range = TRUE)` room for one block (`range_frac`); and a
  flat variogram's refusal message no longer claims the data "never reached
  a sill" without mentioning that a structureless variogram ends there too.
  The README counts six refusals, not five.

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

* `?fit_bayesian_spatial_model`: `check_convergence` says the checks write
  WARN log lines and set `$info$convergence_ok` rather than "issue warnings",
  and that under cmdstanr nothing is raised as an R warning; the basis
  adequacy check is described as logged, and as changing neither
  `convergence_ok` nor `print()` (the argument's text listed it among the
  checks that set `convergence_ok` to `FALSE`); the `family` argument and
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
  and `plot()` on a profile use `set.seed(4)`, on which their comments hold
  (C_p interior at 28 cells, reliability on the range floor); with
  `set.seed(2)` C_p had moved to the support ceiling.  The `plot()` example
  no longer promises four panels on points that have no elbow.

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
  introduction says `determine_optimal_levels()` still returns a count
  there, with a warning.

* `?coerce_to_points` (`tmp_project`) states the rule a CRS-less layer is
  read by: the lon/lat heuristic of `ensure_projected()`, not merely lying
  inside the lon/lat envelope.

* `?summarize_by_cell` now gives the derivation and the measured coverage of
  the small-sample rescaling that goes with every data-derived design
  effect, which the package page said were there; the package page lists it
  among the defaults that were chosen, not cited.

* `?estimate_sac_range` now says when the directional maximum is logged
  (only when it stands in for a singular or non-converged all-pairs fit, at
  WARN only above a ratio of 1.5, naming the directional ranges), and that a
  range shorter than the first lag bin can come back several times too long
  (a true 24 m returned 93--479 m at n = 1500, 19--32 m at `cutoff = 0.1`),
  so a variogram at its sill in the first one or two bins calls for a
  smaller `cutoff` whatever range was fitted.

* `?make_folds` no longer says `auto_range` "fits directional variograms to
  account for anisotropy".  It sizes blocks from the omnidirectional range,
  and for a field known to be anisotropic the page now points to
  `directional_fitted`.  An accepted range too wide for two blocks makes
  `make_folds()` stop with an error, which the page now says; it had claimed
  the grid "does not collapse to a single block".  The log line announcing a
  lowered `k` just before that error is gone.  The page no longer says NNDM
  never pushes a point's nearest-neighbour distance past `phi`: the last
  exclusion can take it past by one neighbour step, as in `CAST::nndm()`.
  The description lists all five methods.  The page now says that `k = 1` is
  raised to 2 by the three k-fold methods, and that only `block_kfold`
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
