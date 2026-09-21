# Changelog

## spatialkit (development version)

### New features

- [`summary()`](https://rdrr.io/r/base/summary.html) on a
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  puts every criterion’s pick in one table: the level each prefers, the
  flat region around it, whether a ladder bound is doing the choosing,
  and the levels that lie in every band. Reading the criteria off a
  profile took one
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  call per criterion, and the comparison had to be assembled by hand.
  The intersection of the bands is often empty, which is a result rather
  than a failure: it says the field has no single resolution that
  satisfies every way of asking. Nothing in the table chooses, and the
  help page says so. Each region comes back in full on the `"bands"`
  attribute.

- Every `cv_*()` result now says what became of each fold. `fold_status`
  is a data.frame with one row per fold supplied — `fold`, `status`,
  `message` — where `status` is `"ok"`, `"error"` (the fit or its
  [`predict()`](https://rdrr.io/r/stats/predict.html) threw; `message`
  is the error text), `"skipped"` (nothing scorable: too few matched
  rows, a prediction of the wrong length, no finite observed/predicted
  pair), `"dropped"` (an empty test set or fewer than two training rows
  once incomplete rows were removed, so the fold never reached the
  fitter) or `"worker_error"` (a parallel worker died). The per-fold
  error text was already collected and thrown away except when *every*
  fold failed, so a partial failure — “3 of 5 folds produced
  predictions” — left its causes only in console scrollback, which a
  script, a `callr` job or a knitted document does not keep. This is
  worth most where a run is expensive: a
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
  fold whose sampler failed now names the reason in the returned object.
  Beside it, `orphan_rows` holds the row IDs no fold names (they enter
  no training set and are never scored — non-empty only when the folds
  were built on a different or subsetted layer), `n_unknown_ids` counts
  fold entries naming rows the data does not have, and `n_dropped` the
  rows
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
  removed before any fold was fitted. The four together account for
  every row and every fold, so `n_folds_attempted - n_folds_succeeded`
  never has to be explained from the log.

- [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
  records what it removed. `attr(x, "dropped")` is a list with `n`,
  `n_geometry`, `which` (positions in the input), `row_id` (when the
  layer carries `..row_id`) and `reason`, one of `"geometry"`,
  `"missing"` or `"non_finite"` per dropped row. The three masks behind
  that decision were already computed and collapsed into a log line; the
  row identities reached nothing, and through the eight-plus internal
  call sites even the count was invisible, so a fit’s `$n` was the
  post-cleaning row count with nothing saying how many rows were lost or
  why. Every fit now carries the count as `$info$n_dropped`, and every
  `cv_*()` result as `n_dropped`. Silently losing a third of the rows is
  a classic cause of a suspiciously good score.

- `make_folds(method = "block_kfold")` returns the block design it built
  the folds from. `assignment` gains a third column, `block_id`;
  `params` gains `blocks` (an `sf` layer of the block polygons in the
  CRS the folds were built in, numbered to match), `block_sizes` (points
  per block, indexed by `block_id`, so a block kept empty by
  `drop_empty_blocks = FALSE` shows as a zero) and `fold_blocks` (which
  blocks were packed into each fold). `blocks$source_row` is the row
  each block came from in the layer it originated in — a cell’s index in
  the full `grid_nx` by `grid_ny` grid, or the row of the `blocks`
  argument — because dropping the empty blocks renumbers the rest: nine
  supplied blocks of which three are empty come back as six rows
  numbered 1 to 6, and a join by row position would mis-attribute every
  block after the first gap. `blocks[params$blocks$source_row, ]`
  recovers them with their own columns and in their own order. The folds
  account for every block exactly once, empty ones included, so a fold’s
  territory on the map is all of its blocks rather than only those
  holding points; and `blocks_used` is the number of blocks the design
  has, which equals `nrow(params$blocks)` and
  `length(params$block_sizes)`, with `sum(params$block_sizes > 0)`
  giving how many of them hold points. The grid **is** the design of a
  blocked cross-validation: without it a user could not draw the blocks
  over their data, see that 40 of 64 blocks were empty, or tell whether
  a fold is one contiguous region or several.
  [`plot_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_folds.md)
  now draws those outlines under the points when the folds carry them,
  and takes `blocks = FALSE` to suppress them.

- [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)’s
  returns are uniform. The four directional variograms and their fits
  run unconditionally on every call, and the rejected-range paths used
  to discard them — precisely where a user most needs to know whether
  the field is anisotropic. Every classed return now carries
  `directional`, `anisotropy` and `anisotropy_used`, and three new
  attributes report the sweep rather than collapsing it:
  `directional_status` (per azimuth, why that direction is `NA` in
  `directional` — `"ok"`, `"over_cutoff"`, `"not_converged"` or
  `"no_fit"`, which were indistinguishable before), `directional_fitted`
  (the range each direction’s fit reported whether or not it was usable,
  so a refused directional range — the most informative number in an
  anisotropic failure — stays recoverable) and, under
  `keep_directional_fits = TRUE`, `directional_fits` (each direction’s
  empirical variogram and fitted model). Those four objects are off by
  default because they dominate the result when present — 42.2 KB of a
  58.8 KB object at n = 400, against 16.6 KB without them — and
  `make_folds(auto_range = TRUE)` calls this on every build.
  [`print()`](https://rdrr.io/r/base/print.html) names the reason and
  the refused value for a direction it cannot use.

- A roster of quantities the package already computed and dropped are
  now returned.
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  attaches the Kish ICCs it estimated as `attr(, "icc")` whether or not
  either was large enough to apply — the case with no `"deff_applied"`
  is exactly the one where a user wants to know what the ICC came out as
  — and the `"variogram"` path adds `deff_rows`, the per-cell design
  effect at the cell’s row count that the log line reduced to a median
  and a max.
  [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
  reports the features that matched more than one polygon and had the
  `tie_break` rule decide for them, as `attr(, "ties")` and a log line;
  a tie-break firing on a third of the features means the polygon layer
  overlaps and every cell count built from it is suspect.
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  attaches `crs_choice`, the projections it considered with each one’s
  measured worst-case distance error. Every path that picks a local
  projection reports what it picked, the two that compare nothing
  included: a UTM zone on a local extent, and the equal-area projection
  chosen for a layer straddling the antimeridian.
  [`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)
  returns the residual `kurtosis` the randomisation variance conditions
  on, the design rank `p` behind the residual moments, `exact`, whether
  those moments are exact for these residuals, and `weights_summary`, a
  description of the weight matrix it used: `n`, `storage` (the matrix
  class), `neighbours` (the smallest and largest number of neighbours
  any row has, `NA` for a dense matrix, where counting them would
  allocate a second one), `kept` and `desc`, the line
  [`print()`](https://rdrr.io/r/base/print.html) shows. The matrix
  itself comes back as `weights` under `keep_weights = TRUE` and is
  `NULL` otherwise, because it is n by n: at n = 500 it is 50.3 KB of a
  53.5 KB object in its sparse form, 112.7 KB at n = 120 when the dense
  fallback is taken (going sparse needs both **FNN** and **Matrix**, so
  a no-Suggests install always falls back), and 191 MB at the n = 5000
  that fallback is capped at — against the 3.2 KB everything else
  occupies — and scoring a list of fits would hold one matrix per fit.
  The result is classed `"morans_i"` and prints through a
  [`print()`](https://rdrr.io/r/base/print.html) method that shows the
  statistic, its null and the weights line; `[` drops the class, and
  `$`, `[[` and [`unlist()`](https://rdrr.io/r/base/unlist.html) read
  the result exactly as for a plain list.
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  keeps `info$nonfinite_coef`, the per-row, per-term mask behind
  `n_local_singular`, and
  [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
  keeps `rhat_failed` / `neff_failed`, the parameters that failed each
  convergence check by name — “max R-hat 1.09” is not actionable where
  “`sdgp_gp..x..y` has R-hat 1.09” is.
  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)’s
  model-aware diagnostics gain `knee_k` and `failed_k` (whose
  interpolated WSS entries are not measurements).
  [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
  records the points whose cell assignment was repaired by nearest-cell
  snapping, and how far outside each sat, as `params$snapped` — a
  comment had long said it should.
  [`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
  returns the `scaling` (per-predictor training centre and SD) the
  dissimilarity index is computed in, without which a location’s DI
  cannot be traced to the predictor that put it outside, and
  `n_outliers`, the training DI values the threshold’s fence set aside.
  `get_voronoi_seeds(method = "kmeans")` returns the clustering as
  `attr(, "kmeans")` — which cloud points fed which seed, the cluster
  sizes and the within-cluster sums of squares.
  [`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md)
  reports `criterion_by_name`, `criterion_column` and
  `criterion_verified`, so a script can gate on the case its log calls
  “unverified” instead of reading the label.

- [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)’s
  result gains a
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method.
  `plot(estimate_sac_range(pts, "z"))` draws the empirical variogram,
  with the fitted model and the effective range overlaid where a range
  was identified, and a subtitle saying why not where it was not: the
  variogram never reached a sill, both model fits were singular, or the
  optimiser halted. Nothing is recomputed — the plot reads the
  attributes the estimate already carries — and the same drawing routine
  now serves `plot.spatial_fit(type = "variogram")`, so the two pictures
  agree. The “attached for inspection” messages
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  logs when it returns `NA` used to point at `plot(type = "variogram")`,
  which is the method for a fitted model and could not take the
  estimate; they now point at
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on the
  returned value.

- [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  gains `conf_level`. With `conf_level = 0.95`, every numeric response
  and predictor column gets four more columns beside its `..sd_*` and
  `..se_*`: `..neff_*`, the column’s effective sample size in the cell
  (its non-missing count over its design effect — the per-column version
  of `cell_weight`); `..df_*`, the degrees of freedom the interval uses;
  and `..ci_lo_*` / `..ci_hi_*`, a t interval for the cell mean as an
  estimate of the grand mean, built on the design-effect-corrected
  standard error. The default `conf_level = NULL` returns exactly the
  frame it always did. The degrees of freedom are `n - 1` at `deff = 1`,
  for a numeric `deff` and for `deff = "kish"`, because the interval’s
  spread is estimated from the within-cell variance, whose `n - 1` df
  survive exchangeable correlation whatever the design effect (measured
  95% coverage on the Kish path at an ICC of 0.2 / 0.6 / 0.9: 0.954 /
  0.953 / 0.952; an effective-sample-size df of `neff - 1` gives 0.992 /
  1.000 / 1.000 and is not used). Under `deff = "variogram"` the df are
  the Satterthwaite

  1946. moment-matched df of the within-cell variance under the fitted
        correlation, a fraction of `n - 1` that shrinks with the range:
        0.960 and 0.958 coverage at exponential ranges of 150 and 400 on
        a 1000-unit domain, against 0.931 and 0.918 with `n - 1`. The
        help page’s “Confidence intervals” section has the reasoning and
        the numbers.

- [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
  gains `block_size` and `auto_range`, and now hands `response_var` and
  `predictor_vars` to
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  when it builds the shared fold set. With the defaults the folds are
  the same geometric blocks as before, so an existing comparison does
  not move; what changes is that the fold-leakage diagnostic (above) can
  now fire for the one function that compares models, which until now
  was the one whose folds could never be checked against the range.
  [`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)
  gains `auto_range` for its inner folds for the same reason, and
  records it in `$params`.

- `compare_models_cv()$overall` carries the Bayesian backend’s
  calibration when a Bayesian model ran: `coverage_50`, `coverage_80`,
  `coverage_95` and `mean_CRPS`, the same fold-weighted summary
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
  returns as `predictive_coverage`, with `NA` on the GWR and RF rows. A
  model that predicts well on average and covers badly (Heaton et
  al. 2019) is now visible in the table a user picks from, not only in
  `$bayes_cv`. `model` stays the last column.

- [`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
  gains `replace` and `sample_fraction`, which reach
  [`ranger::ranger()`](http://imbs-hl.github.io/ranger/reference/ranger.md)
  as `replace` and `sample.fraction`; both were already accepted through
  `...`, but are now recorded in `$info` and printed with the fit
  (“Sampling: bootstrap, with replacement (100.0% of rows per tree)”).
  The defaults are ranger’s, so no forest changes. The help page carries
  Strobl et al.’s (2007) case for `replace = FALSE`. Passing the ranger
  spellings through `...` is now refused like the other arguments the
  wrapper sets.

- [`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
  records the method of the folds its threshold came from as
  `$params$folds_method` (`"block_kfold"`, `"random_kfold"`, … from a
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  result; `"labels"` or `"splits"` when the input cannot say) and prints
  it, because the threshold is a statistic of a cross-validated hold-out
  and pairs with the CV error from the same kind of hold-out. An AOA
  built on `random_kfold` folds logs a caution saying it does not belong
  beside a blocked `cv_*()` result.

- [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  gains `detrend = c("ols", "reml")` and `reml_max_n`. A variogram
  fitted to least-squares residuals underestimates the range, because
  the trend fit absorbs part of the long-wavelength variation (Lark,
  Cullis and Welham 2006). Measured for this estimator on simulated
  fields (n = 300, true effective range 300), as the median ratio to the
  estimate from the true field: a white-noise covariate 1.00; a
  spatially smooth covariate 0.97; a linear trend in the coordinates
  0.92; a quadratic one 0.75. `detrend = "reml"` fits the trend and an
  exponential-plus-nugget covariance together by REML with
  [`nlme::gls()`](https://rdrr.io/pkg/nlme/man/gls.html) and returns the
  REML range — 0.95 to 1.06 on the same designs — with the empirical
  variogram of the REML residuals attached for inspection. It is cubic
  in `n`, so it runs on at most `reml_max_n`

  400. points; a fit that does not converge falls back to OLS with a
       warning. The default stays `"ols"`, so nothing changes unless
       asked; the help page’s new section carries the numbers. Iterating
       GLS trend fits against variogram refits (Neuman and
       Jacobson 1984) was measured too and recovers only part of the
       bias (0.80 in the quadratic case), so it was not added. `nlme`
       joins Suggests.

- [`sac_nugget()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md)
  returns the nugget variance behind an
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  result, and every classed result — identified or rejected — now
  carries it as a `nugget` attribute (`NA` when no model could be
  fitted). It was reachable before only by reading the `Nug` row of the
  attached `gstat` model. Results also record `detrend_method` (`"ols"`,
  `"reml"` or `NA`), and with `detrend = "reml"` a `reml` list
  (`n_used`, `subsampled`, `nugget_prop`, `sigma2`).

- [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  refuses a range fitted through an empirical variogram that *decreases*
  with distance over its shorter lags: a net fall of more than 15% of
  the mean semivariance there, weighted by pairs. A model that rises to
  a sill has nothing to identify on such a curve, and the result is `NA`
  with
  `rejected_reason = "empirical variogram decreases with distance"`, the
  refused value in `rejected_range`, and the variogram attached;
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) captions it.
  The shape is what a periodic (hole-effect) structure or a variance
  that differs between a dense cluster and the rest of the layer
  produces — measured on 60 draws each: 98% of fields with a periodic
  component and 100% of the clustered case are flagged, against 0% of an
  exponential field with a range of 100 or more on a 1000-unit extent,
  2–3% at very short ranges, 2% of white noise — and *not* what a trend
  produces, which is a variogram that rises without a sill and is
  refused as before. The one other path that returned a bare `NA` from
  inside a completed fit (a non-positive fitted range) now returns the
  classed, inspectable shape too.

- [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  and
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md):
  the number of cells, scored on every criterion at once.
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  runs a log-spaced ladder of level counts from a floor the
  autocorrelation range implies (`ceiling(area / range^2)`) to a ceiling
  the support implies (`floor(n / min_cell_n)`), fits each level as the
  best of 25 k-means++ restarts, and returns a data.frame with the WSS
  elbow statistic, Mallows’ C_p of the piecewise-constant approximation
  of the response (or of its residuals on the predictors, with the
  nugget from the fitted variogram as the noise variance), the
  standardised residual Moran’s z of the cell means, and an analytic
  reliability of the cell means — the share of their spread that is
  between-cell signal rather than sampling noise, from the variogram
  alone via Krige’s additivity relation (Cressie 1996), the shrinkage
  factor of Fay and Herriot (1979) — plus the cell-support and
  cell-diameter columns and the between-restart spread. A floor above
  the ceiling is reported as a finding rather than resolved silently.
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  reads a level off one criterion together with its flat region and says
  when a bound, not the criterion, is choosing. The flat region is a
  set, not an interval: the criterion curves are not monotone, so a
  region can skip a rung of the ladder, and it is printed as the runs
  the criterion accepts (“19 to 21, 26, 31 to 33”) rather than a range
  that would quietly include the levels it rejected. Two things measured
  before this shipped, both on the help page: on smooth fields with a
  small nugget C_p descends to the support ceiling (every replicate at
  effective ranges 90–900 with nugget 0.3 on a unit sill; interior only
  at nugget 2), and the reliability optimum agrees with the empirical
  one from true block means on simulated fields but is broad — flat to
  within 2 percent over a factor of 3–6 in the number of cells. Read the
  flat region.
  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  is unchanged in shape and keeps its integer-vector interface.

- `select_on = c("all", "split")` on
  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md),
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  and
  [`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md).
  Whenever a selection reads the response — a level count chosen with
  `response_var` and `predictor_vars`, or a predictor set chosen by a
  forward sweep — what is estimated afterwards on the result is
  post-selection, and its standard errors are descriptive rather than at
  nominal coverage (Gao, Bien and Witten 2022). `select_on = "split"` is
  sample splitting: the layer is cut into two spatially blocked halves
  (`make_folds(k = 2, method = "block_kfold")`), the selection runs on
  the first, and the row positions of both come back (as a `"split"`
  attribute on the first two functions, as `$split` on the third) so the
  estimation can be done on the half the selection never saw.
  [`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)
  also returns `score_holdout`: the selected set fitted on the selection
  half and scored on the other, the honest number its selection-internal
  `score` is not. The cost is precision — half the points estimate, and
  a contiguous spatial half is less efficient than an exchangeable one
  (García Rasines and Young 2023). Data thinning (Neufeld et al. 2024)
  and data fission (Leiner et al. 2023) keep the whole sample and are
  noted on the help page, not implemented. The help page also now says
  plainly that supplying both `response_var` and `predictor_vars`
  upgrades the level-selection criterion to `"combined"`, so the
  selection depends on the response without that having been asked for.

- [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  gains `blocks`: a polygon layer (`sf` or `sfc`) to use as the blocks
  of `method = "block_kfold"` in place of the grid it would otherwise
  build — the `$cells` of a
  [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
  result, hexagons, watersheds, administrative units, the `$blocks` of
  `blockCV::cv_spatial()`. Each point takes the block that contains it
  and the blocks are assigned to folds exactly as grid cells are, so the
  fold builder can now consume every shape the tessellation half of the
  package produces. The grid-sizing arguments and `boundary` are ignored
  with a log line, `auto_range` compares the estimated range against the
  blocks instead of resizing them, and the leakage warning uses the
  median over blocks of the side of the square with the block’s area
  (`params$block_scale`). Points inside no block are assigned to the
  nearest one with a warning that counts them, except points within a
  millionth of the extent of a block — an edge that reprojection moved
  by a rounding error; points inside more than one block take the first,
  with a warning when the blocks concerned overlap in area rather than
  share an edge. `params` gains `n_blocks` (before empties were
  dropped), `blocks_supplied` and `block_scale` on every `block_kfold`
  result; `grid_nx`/`grid_ny` are `NA` for supplied blocks.

- [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  gains `balance_tol`, the largest-to-smallest fold size ratio above
  which `block_kfold` reports its folds as imbalanced. The check was
  always there at a hard-coded 3:1, and it was a log line only; it is
  now an R warning (`Inf` disables it), the ratio achieved is returned
  as `params$balance_ratio`, and the help page says which methods
  balance what: only `block_kfold` balances point counts, by packing
  blocks largest first into the fold with the fewest points so far. No
  search over packings was added, because the packing is not where the
  imbalance comes from. Measured against the optimum by enumeration (two
  folds, up to ten blocks, heavy-tailed sizes), the greedy packing is
  optimal in 72 percent of cases and within two points of optimal on
  average; a local search with 30 random restarts moved the ratio by
  0.004 on average over 480 block-size vectors and never brought one of
  the 91 above 3:1 below it. The remedy for an imbalance past the
  tolerance is the block design — the warning now says so — and `blocks`
  is how to supply one: on clustered layouts where the geometric grid
  exceeded 3:1 in 22 percent of draws (median 1.8, worst 6.3), Voronoi
  cells around 15 `get_voronoi_seeds(method = "kmeans")` seeds never
  exceeded 1.3 (median 1.09). With the default tolerance the folds of
  every existing call are unchanged.

- `metrics` on
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
  [`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
  and
  [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md):
  a scoring function of your own, `function(y, yhat)` returning a named
  numeric vector (a named list of scalars or a one-row data frame also
  serve), applied the way the built-in metrics are — once per fold, so
  each name becomes a column of `fold_metrics`, and once to the pooled
  out-of-sample predictions, so each name becomes a column of `overall`
  — on the same finite pairs `RMSE` uses. This is the way to score what
  the Gaussian set cannot: a Poisson deviance, a log score on a
  probability, a weighted loss. The contract is strict where it should
  be (every element named, names unique and not a built-in column, one
  number per name — anything else is an error, because a scoring
  function of the wrong shape is a mistake to surface) and forgiving
  where it should be (a function that throws on a fold is logged and its
  columns are `NA` there; a fold is never dropped for it). The empty
  frames of a run where every fold failed carry the columns, typed, when
  the function can be called on zero-length input.
  [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
  hands one function to every backend and protects it like the fold
  arguments, so the columns of its `overall` are comparable across rows.
  `fold_info_fn` is documented as the per-fold half of the same
  mechanism, with access to the fitted object and the held-out layer.

- [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  gains `purpose = c("distance", "area")`. The default is what it always
  did: for lon/lat input, the candidate that distorts distances least.
  With `"area"` — densities or rates per cell are going to be computed —
  the choice is made among equal-area projections only (a Lambert
  azimuthal centred on the data, or an Albers conic where its parallels
  do not degenerate, whichever distorts distances less), which a UTM
  zone never enters, and global coverage gets Equal Earth rather than
  Web Mercator. Already-projected input is still returned untouched, but
  its area distortion over the extent is now measured — the spread of
  planar-to-geodesic area ratios over probe polygons — and logged as a
  warning above 1 percent. Measured: a UTM zone edge to edge 0.25
  percent, a 2.5-degree extent inside one 0.04 percent, the conterminous
  United States forced into one zone 14 percent, Web Mercator over 2.5
  degrees of latitude at 48N 4 percent, an equal-area projection a few
  tenths of a percent (the sphere the geodesic areas are computed on
  against the ellipsoid).

- [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  gains `area = TRUE`: with `cells_sf`, the result carries `cell_area`
  (planar, in the squared units of the cells’ CRS) and `n_per_area`, a
  point density; a rate of anything else is its `agg_funs` sum over
  `cell_area`. The request is refused with an error — not answered with
  a number — when the cells’ CRS distorts areas across them by more than
  1 percent by the measurement above, because a density is a comparison
  between cells and means nothing where the map scale differs from one
  cell to the next; the message names the CRS, the figure and the
  remedy. A cell with no observations gets `NA`, not zero. The measured
  spread is attached as `attr(, "area_error")`.

- `build_tessellation(approx_n_cells = )` and `get_voronoi_seeds(n = )`
  accept what the level-selection step returned: the integer vector of
  ranked candidates from
  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  (its first element is used), a
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  result (its `$best`), or a
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  (read with
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  at its default criterion). Both were hard errors before, so nothing
  that worked changes; the count used and where it came from are
  recorded as `params$approx_n_cells` / `params$approx_n_cells_from` and
  as `attr(seeds, "n_from")`. The two functions the pipeline documents
  as a pair are now connected:
  `get_voronoi_seeds(n = determine_optimal_levels(pts))` needs no number
  carried between the calls by hand.

- Five diagnostic plots that show the curve behind a chosen point, the
  folds behind a pooled number, or the distribution behind a count. None
  recomputes anything; each draws what the result already carries.

  - `plot_cv_metrics(cv, metric)`: one point per fold, sized by the
    held-out rows it contributed, with the pooled value from `overall`
    as a dashed line; a
    [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
    result gets one panel per model on a shared scale. Any column of
    `fold_metrics` can be drawn, including backend extras and columns a
    `metrics` function added; a column that is `NA` in every fold is
    refused with the reason (`Adj_R2` without `p`, coverage without
    draws) rather than drawn empty, and a per-fold extra with no pooled
    counterpart draws without the line and says so.
  - [`plot.aoa()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.aoa.md):
    the dissimilarity index of the prediction locations against the
    cross-validated training DI (ECDFs, or a histogram with the training
    curve), threshold marked, with the share outside and how close the
    inside ones run to the edge in the subtitle, and whether the
    threshold came from cross-validated folds in the caption.
  - `plot.spatial_fit(type = "variogram")` overlays the response’s own
    variogram (hollow points, dashed fit) on the residual variogram, on
    the same points and lags, so the structure the model absorbed is the
    gap between the two curves. The caption compares the sills only when
    both ranges were identified; `response = FALSE` restores the
    residual curve alone.
  - `plot_calibration(cv)`: observed against nominal coverage of
    [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)’s
    posterior predictive intervals, pooled (blue) and per fold (grey),
    with the diagonal and a one-line verdict. The levels are read off
    the `coverage_*` column names, so
    `coverage_levels = seq(0.1, 0.9, by = 0.1)` gives a full curve; the
    default three levels are unchanged.
  - One sweep drawer behind three methods:
    [`plot.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.resolution_profile.md)
    (a panel per criterion, the level each selects marked, its flat
    region shaded, a note when a bound rather than the criterion is
    choosing);
    [`plot.feature_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.feature_selection.md)
    (the accepted variable’s score at each step as the path, every other
    candidate faint, the stop in red, the hold-out score as a separate
    mark when `select_on = "split"` computed one — and a caption saying
    whether the intercept-only model was scored, since for the RF and
    GWR backends it usually is not, so the path starts at the first
    variable);
    [`plot.gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.gwr_model_selection.md)
    (every model’s AICc against its size, the best of each size joined,
    the winner marked, its lead over the runner-up in the subtitle).
    [`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)’s
    result now carries class `"feature_selection"` so
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html) finds the
    method; it is the same list otherwise.

- [`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md):
  the same cross-validation at a ladder of block sizes, with random
  folds as the leaky reference, returned as a table with the
  fold-to-fold spread at each size and the estimated autocorrelation
  range alongside;
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
  curve with the range marked. Blocks smaller than the range leak, so
  the curve rises from the random-fold value towards the range and
  plateaus beyond it, and the height of the rise is what the random-fold
  number overstated — measured on simulated fields with a 100–130-unit
  range and a random forest with coordinates, the plateau begins at one
  to two times the estimated range. Each size is a full
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
  so a fit budget (`max_fits`, default

  60. refuses to start rather than run past it, and the ladder drops
      sizes at which the grid holds fewer than `k` blocks so every point
      on the curve is a `k`-fold cross-validation of the same shape.

- [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  keeps its local collinearity survey. Every fitting window — not a
  sample of 30 — has its kernel-weighted local design’s scaled condition
  index computed, the way Wheeler and Tiefelsdorf diagnose GWR
  collinearity, and the fit carries it as `info$local_collinearity` (one
  row per observation: coordinates, window size, condition index), with
  `info$n_local_collinear`, `info$n_local_singular` and the global
  `info$condition_index` beside `AICc`. The warning is now the exact
  fraction of locations rather than a sampled one; its wording and its
  thresholds (a quarter of the locations, or any) are unchanged. The
  weighted survey sees what the unweighted spot-check could not: a
  bisquare window’s edge points contribute almost nothing to the fit, so
  they contribute almost nothing to its conditioning.

- `plot(fit, type = "coefficients")` for a GWR fit maps one local
  coefficient (`term`) at the training locations, which is the reason to
  fit GWR at all — and masks the locations where it is not to be
  believed: a collinear local design (condition index above 30, or
  singular) or a non-finite coefficient is drawn hollow and grey,
  counted in the subtitle, because the smooth surface a naive map draws
  over them is the picture of an unstable estimate. `mask = FALSE` draws
  them anyway and says how many it is drawing.

- [`kriging_adequacy()`](https://elkronos.github.io/gis_modeling_toolkit/reference/kriging_adequacy.md):
  what a block-kriging aggregator would deliver on a set of cells,
  computed beside the plain means and changing none of them. Per cell,
  from a fitted variogram
  ([`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)’s,
  or estimated here): the block-kriging estimate and variance, that
  variance as a share of the sill (`kr_ratio`, the coverage score — near
  1 the estimate is the global mean), whether it exceeds the
  design-based `s^2/n` of the plain mean (`kr_exceeds_design`), and the
  kriged-minus-plain shift in standard errors (`kr_shift`); plus the
  variance of the standardised errors from blocked cross-validation
  (`attr(, "cv")$zscore_var`), which is about 1 when the kriging
  variance is right. Measured on simulated exponential fields: 0.93–1.07
  with the true variogram, 0.85–1.01 with the estimated one; under
  blocked folds it checks the sill and range rather than the nugget
  (0.95–1.24 with the nugget understated tenfold), and random folds are
  the instrument for the nugget. The comparison the function exists for:
  under uniform sampling kriged and plain means differed by more than
  one standard error in 11–24 percent of cells; under clustered sampling
  in 34–63 percent, with 3–27 of 16–64 cells empty and kriged anyway.
  This is the first kriging path in the package
  ([`gstat::krige()`](https://r-spatial.github.io/gstat/reference/krige.html)
  and
  [`gstat::krige.cv()`](https://r-spatial.github.io/gstat/reference/krige.cv.html));
  its model families are the ones the package interprets elsewhere, and
  any other is refused by name.

- `MAPE` and `SMAPE` now say how many rows they were averaged over.
  Every metrics frame —
  [`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md),
  [`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md),
  and the `overall` and `fold_metrics` of
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
  [`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
  and
  [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
  — gains two trailing integer columns, `n_MAPE` and `n_SMAPE`: the rows
  each percentage error actually used once those where its denominator
  is zero were dropped (`y == 0` for MAPE; `|y| + |yhat| == 0` for
  SMAPE). They equal `n` (`n_pred` in the CV frames) when nothing was
  dropped, and are `0` in an empty frame. The values themselves are
  unchanged: a MAPE over 58 of 120 rows is the same number 2.0.0
  reported, but it now arrives labelled, where before nothing in the
  frame recorded that it was a subset average. `print(summary(fit))`
  appends “(over k of n rows)” to its SMAPE line when the two differ.
  The columns sit after `Adj_R2` so code addressing the seven metric
  columns by position is unaffected; code pinning the exact column set
  needs the two names added.

### Bug fixes

- `plot(fit, type = "variogram")` no longer runs its subtitle off the
  edge of the figure. ggplot2 clips a label that is wider than the plot
  instead of wrapping it, and the sentence saying why no range was
  identified is up to 150 characters, so on a six-inch figure it was cut
  mid-word. Labels built from a fit’s own numbers are now wrapped at
  draw time.

- [`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md)
  could not give IDs to a tessellation this package had just built. It
  repaired the geometry in the layer’s own CRS and then transformed it
  to the sort CRS, but validity is a property of the geometry in the CRS
  it is measured in: two vertices a centimetre apart in a projected CRS
  can land on one longitude and latitude, and s2 calls the ring
  degenerate. On a clipped hex tessellation of North Carolina, 2 of 18
  cells that are valid projected are invalid once transformed, and
  [`st_centroid()`](https://r-spatial.github.io/sf/reference/geos_unary.html)
  on one of them aborted the call with “Loop 0 is not valid: Edge 1 is
  degenerate”. The sort copy is now repaired after the transform as
  well. The geometry returned is still the caller’s own, and the same
  cell gets the same ID whether the layer arrives projected, in lon/lat
  or in Web Mercator.

- [`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)
  now says when `fit_fn` is ignoring the variables it is handed. The
  learner it takes is a function of `(train_sf, predictor_vars)`, but
  the one
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
  takes is a function of `train_sf` alone, and a learner written for
  that — `function(train_sf, ...)` — swallows the second argument and
  fits the same model every time. Nothing errored, because the training
  layer still carries every column: each candidate scored exactly what
  the intercept-only model scored, no candidate improved on it, and the
  result was an empty `selected` and an `NA` `score` with no word about
  why. Two different predictor sets do not produce the same
  cross-validated metric to the last digit, so that pattern is now
  recognised at the first step and reported as a warning naming the fix.
  The result is still returned.

- [`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md)
  warns when both `cellsize` and `target_cells` are supplied. It already
  did for `cellsize` and `n`, and the documentation says to supply
  exactly one of the three, but `target_cells` was dropped in silence
  when `cellsize` was present. Through
  [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
  that meant `method = "hex", approx_n_cells = 25, cellsize = 10`
  returned however many cells a 10-unit lattice holds and said nothing
  about the 25. `cellsize` still wins; the override is now logged like
  its sibling.

- [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
  warns when `approx_n_cells` or `cellsize` is supplied with
  `method = "voronoi"` or `"triangles"`. Both arguments size the hex and
  square lattices and nothing else — Voronoi grows one cell per input
  point and Delaunay one triangle per neighbouring triple — and both
  used to be dropped in silence. So
  `build_tessellation(pts, method = "voronoi", approx_n_cells = 25)`
  returned one cell per observation: the degenerate nearest-neighbour
  case, where every cell holds a single point, there is no within-cell
  variation and every standard error is `NA`. The warning names the
  argument and points at
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
  which is where a Voronoi cell count is actually set. Under
  `method = "voronoi"` it adds that `params` does not record the request
  either, so a saved result carries no sign of it; that branch returns
  [`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md)’s
  own list, which has no slot for the argument, whereas the triangles
  branch does echo `approx_n_cells` back. The warning fires under
  `quiet = TRUE`, which gates this function’s
  [`message()`](https://rdrr.io/r/base/message.html)s and is documented
  not to silence R warnings.

- [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  fits each k as the best of 25 k-means++ restarts (Arthur and
  Vassilvitskii 2007; Fränti and Sieranoja 2019; Steinley 2003) instead
  of `stats::kmeans(nstart = 5)`. The WSS curve is read for its shape,
  and with a handful of random restarts it carried optimisation noise:
  on eight-cluster layouts a sweep over k = 1..30 rose at one or two
  steps in three of five draws, and an earlier form of the elbow rule
  once selected such a bump. With the new budget the same sweeps rose at
  no step. A curve that still rises is now logged as a warning naming
  the number of rising steps, and the model-aware diagnostics carry it
  as `wss_bumps` beside `wss_spread` (the relative spread of WSS across
  restarts at each k) and `nstart`. A selection made on a curve that had
  a bump can differ from before — those were the cases that were wrong;
  a clean curve gives the same answer. Under a model-aware criterion the
  function also now warns *before* the sweep when `max_levels` leaves no
  k above the nine-cell floor, rather than fitting every k first and
  falling back afterwards.

- `make_folds(method = "block_kfold")` can now raise its “block
  dimension \< autocorrelation range” warning. The comparison was always
  there, but the range it compared against was estimated only under
  `auto_range = TRUE`, the one setting in which the blocks had already
  been sized from that range and the warning could never fire; on every
  default call the diagnostic was dead code. With `auto_range = FALSE`
  (the default) and a `response_var` to hand — always the case when a
  `cv_*()` function builds the folds — a range is now estimated for the
  diagnostic alone. It sizes nothing: the blocks are the same geometric
  blocks as before and the folds do not change. A hand-set `block_size`
  below the range raises the same warning. The estimate’s own log lines
  stay off the console, and the check is skipped (with an INFO log line
  saying so) when `gstat` is not installed or there are fewer than 30
  points.

- `fit_rf_model(include_coords = TRUE)` logs its caution once per
  session rather than once per fit. Inside a five-fold
  [`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
  or a twenty-fit
  [`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md)
  the same paragraph printed on every fit, which reads as twenty
  problems rather than one decision; the message now says it will not
  repeat.

- [`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
  reports what ranger actually objected to. ranger diagnoses a bad
  argument in its C++ layer, writes the diagnosis straight to stderr and
  then throws “User interrupt or internal error.” — so `mtry = 99` on a
  two-predictor forest printed “mtry can not be larger than number of
  variables in data. Ranger will EXIT now.” to the console and raised an
  error naming neither the argument nor the problem. That line is not an
  R condition, so
  [`suppressMessages()`](https://rdrr.io/r/base/message.html),
  [`withCallingHandlers()`](https://rdrr.io/r/base/conditions.html) and
  [`tryCatch()`](https://rdrr.io/r/base/conditions.html) all missed it:
  it escaped every handler to the console (and into CI logs, where it
  reads as an error from a test that is passing), while `cv_*()`
  recorded the placeholder as the fold’s cause. The message stream is
  now diverted for the duration of the call, so the diagnosis becomes
  the reported reason — `fold_status$message` included — and nothing is
  printed behind the caller’s back. Output from a call that succeeds is
  passed through unchanged, and when the stream is already diverted
  (under testthat, knitr or `capture.output(type = "message")`, where
  only one sink is permitted) the call runs exactly as before.

### Documentation

- Every figure in the vignettes carries alt text, which is what a screen
  reader announces and the only thing a reader gets when an image fails
  to load. Each one states what the picture shows and what it is there
  to demonstrate, rather than naming the axes.

- Every help page now ends with a “See also” that leads somewhere; more
  than a third had none, every S3 method among them. Four families are
  new: **spatial data preparation**
  ([`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md),
  [`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md),
  [`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md),
  [`clip_target_for()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clip_target_for.md),
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)),
  **package options and caches**
  ([`spatialkit_quiet()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md),
  the two cache clearers,
  [`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md)),
  **methods on a fitted model** (the
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods of the
  three backends, which document one contract between them) and **print
  methods** for the other result objects. The seed generators, the
  cached grid constructor and
  [`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md)
  join **tessellation**,
  [`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
  joins **model evaluation**,
  [`sac_nugget()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md)
  joins **cross-validation**,
  [`gp_lengthscale_bounds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gp_lengthscale_bounds.md)
  joins **model fitting**, and
  [`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
  is now in **prediction** as well as cross-validation. The website’s
  reference index moves the two cache clearers into the same group, so
  it and the help pages agree.

- Ten numbered scripts are installed with the package, in
  `system.file("scripts", package = "spatialkit")`. `00-run-all.R` runs
  them in order; each of `01-` to `10-` is self-contained and covers one
  topic: tessellations, resolution, fold schemes, block sizing, fitting
  and diagnosing, model comparison, prediction surfaces and the area of
  applicability, feature selection, GWR and the Bayesian GP. They run on
  a simulated field with known structure, print what they are doing, and
  state what to look for in a figure before drawing it. No script
  asserts a number it did not compute, so a run that finishes is one
  whose claims held on the machine that ran it. `SPATIALKIT_TOUR_OUTPUT`
  writes the figures to a folder instead of the device;
  `SPATIALKIT_TOUR_PAUSE = "no"` skips the per-figure pause. A script
  skips the part that needs an absent package with a message naming it,
  and scripts 02, 09 and 10 skip themselves when gstat, GWmodel or brms
  is absent.

- Four vignettes join `spatialkit_nc_demo`: `getting-started`
  (installing, what the coordinates are in, the pipeline from points to
  a scored model and a map, and a glossary of the terms that recur),
  `resolution` (the ladder, the four criteria and why they disagree),
  `spatial-cross-validation` (the five fold schemes, block sizing, and
  reading a CV result down to the last row) and `diagnostics` (residual
  autocorrelation, aggregation standard errors, kriging adequacy, area
  of applicability, and two ways to leak). Each is executed at build
  time and gates itself on the optional packages it needs.

- The README is generated from `README.Rmd`, so its figures and every
  number in it are computed when it is built rather than pasted in. It
  is about half its previous length: the material that had accumulated
  in it moved to the vignettes above, and what remains is what the
  package is for, a quick start that shows a real cross-validation gap,
  the troubleshooting list, and pointers to the rest.

- New vignette `reporting`: what leaves the session at the end of a run.
  The first half is the regions as a file someone else can use. Grouping
  a layer you already have with
  `assign_features_to_polygons(largest = TRUE)` (100 North Carolina
  counties into a hex grid of 18 regions, 12 of them populated, one row
  per county), answering whether a location falls in one with
  `keep_unassigned = TRUE`, IDs that survive a reprojection, and why the
  aggregates go into a GeoPackage instead of a shapefile: of the 10
  columns
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  produces, 2 come back from a shapefile with their names intact. The
  second half is a worked report of six numbers, each read out of an
  object the run already produced, with what each one is there to stop a
  reader believing (Roberts et al. 2017; Meyer and Pebesma 2021; Heaton
  et al. 2019). Its example reports a run that fails its own checks,
  which is the case the section exists for.

- The documentation is published as a website at
  <https://elkronos.github.io/gis_modeling_toolkit/>: the README, every
  help page grouped by pipeline step, the six vignettes as articles, and
  this changelog. It is rebuilt from `main` on every push, so it
  describes the development version; the “development version” heading
  at the top of this file lists what the CRAN release does not have yet.

- [`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)
  now says what `$score` is: the cross-validated metric of the winning
  set at the final step, which is the selection criterion and is
  optimistically biased by the selection itself (Cawley and Talbot 2010)
  — not a performance estimate of the selected model. The honest
  estimate comes from running the selection inside
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)’s
  `fit_fn`, which the page now spells out.

- [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
  documents that `family =` accepts any `brms` family — zero-inflated
  and hurdle counts, negative binomial, Bernoulli, beta and ordinal
  responses all reach
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) with
  the spatial GP term intact — and gains a worked zero-inflated Poisson
  example. The same page now records a trap: a family object `brms`
  cannot name skips the response-type check entirely rather than falling
  back to the gaussian rule, so a malformed `family` buys less
  validation, not more.

- [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
  gains a “Spatial confounding” section: a coefficient estimated beside
  a spatial random effect is a different estimand from the non-spatial
  one (Zimmerman and Ver Hoef 2022), can shrink toward zero when the
  response is smoother than the covariate (Bolin and Wallin 2025), and
  the honest diagnostic is to report both side by side. The section
  names both sides of the restricted-spatial-regression dispute (Hughes
  and Haran 2013; Hanks et al. 2015; Khan and Calder 2022) and the
  remedies expressible in a `brms` formula (Marques, Kneib and Klein
  2022; Guan et al. 2023).

- [`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
  gains a “Which metrics survive a non-Gaussian response” section,
  inherited by
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
  and
  [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md):
  RMSE and MAE are defined for any numeric response; MAPE, SMAPE and
  R-squared are Gaussian-shaped; CRPS and interval coverage from
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
  are the proper scores for a count or bounded response. It also records
  that the all-folds-failed `fold_metrics` frame omits the `coverage_*`
  columns.

- [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  documents what a count or other mean-variance-linked response does to
  the variogram, why detrending with `predictor_vars` helps but does not
  fix it, and what to do instead.

- [`?spatialkit`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-package.md)
  now reads its nine-step pipeline as an argument — steps 1 to 4 are the
  claim, steps 5, 7 and 9 the evidence that makes it checkable — and
  gains a “Defaults and their sources” section listing which defaults
  rest on a citation and which were chosen, so the two are not mistaken
  for each other.

- [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
  documents its name collision with `blockCV::cv_spatial()`, which
  builds folds where this one runs them, and that `blockCV`’s
  `$folds_ids` is accepted directly as `folds` everywhere.

## spatialkit 2.0.0

CRAN release: 2026-09-11

Everything below is relative to **1.0.0** (published on CRAN
2026-08-07). The major bump is warranted: three exported functions are
removed, and several defaults change the result of a fit or a
comparison, so the same script can get a different answer. Both are
under “Breaking changes” — read that section before upgrading a running
analysis.

These notes describe what changed **for a user of 1.0.0**. A good deal
of this release was written after 1.0.0 and then revised before
shipping; defects that existed only between those points are not listed,
since no released version behaved that way. The commit history has that
record in full.

Throughout, *raises a warning* means a genuine R
[`warning()`](https://rdrr.io/r/base/warning.html) — one
`tryCatch(warning = )` catches,
[`suppressWarnings()`](https://rdrr.io/r/base/warning.html) suppresses
and `options(warn = 2)` escalates. *Logs a warning* means a `logger`
message in the `"spatialkit"` namespace, which none of those touch.

### Breaking changes

#### Corrections that change results

Each item below was measured and independently reproduced before it was
touched; the figures quoted are from those reproductions, so you can
judge whether an item affects an analysis you have already run.

- **[`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)
  no longer puts weight on a point’s own residual.**
  [`FNN::get.knn()`](https://rdrr.io/pkg/FNN/man/get.knn.html) reports a
  point’s OWN index among its neighbours whenever exact duplicate
  coordinates are present, which put `1/k` on the diagonal of a matrix
  Moran’s I is only defined for with a zero diagonal. On 40 sites x 4
  repeats with a response carrying **no** spatial structure, 120 of 160
  rows gained a self-weight, mean I came out at +0.086 against E\[I\] =
  -0.0063, and **77% of samples were “significant” at p \< 0.05**
  against a nominal 5%. Repeat observations at one site are exactly what
  `make_folds(method = "leave_location_out")` is for, so this was a
  mainstream input. The dense fallback never had the fault, so the
  statistic also depended silently on whether **FNN** happened to be
  installed; the two paths now share one neighbour lookup and agree
  exactly.

  Requesting `k + 1` neighbours and dropping self is **not** sufficient
  on its own — the slot self occupied displaced a genuine co-located
  neighbour and left a farther point standing in for it (75 of 400
  retained pairs sat at distance 121 where a neighbour at distance 0
  existed). Duplicate coordinates are now grouped and answered exactly.

- **[`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)
  gains a `null` argument, defaulting to `"auto"`.** Model residuals are
  not exchangeable — they are orthogonal to the design matrix — so the
  classical randomisation moments are wrong for them. At n = 120 with
  six smooth covariates and independent errors, OLS residuals had mean I
  = -0.031 against the exchangeable E\[I\] = -0.008, and the z-score
  averaged -0.54 with sd 0.90 instead of 0 and 1. The Cliff & Ord (1981,
  sec. 8.3) regression-residual moments restore mean z = -0.09, sd 1.03
  and a 4.3% rejection rate against a nominal 5%, and **agree with
  [`spdep::lm.morantest()`](https://r-spatial.github.io/spdep/reference/lm.morantest.html)
  to machine precision** (verified at 1e-16 through the public
  function). `"auto"` applies them only when the fit’s residuals really
  are the OLS residuals on the rebuilt design, which a forest’s and a
  working GWR’s are not; the null actually used is reported in the
  return value.

- **[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  standard errors under a design effect were too small.**
  `s / sqrt(n / deff)` corrects the mean’s variance for clustering but
  leaves `s^2` biased low by the same clustering: for exchangeable
  correlation rho, `E[s^2] = sigma^2 (n - deff)/(n - 1)`. The two errors
  compound. Measured 95% CI coverage at n = 20: **0.905 at rho = 0.3,
  0.796 at rho = 0.6, 0.632 at rho = 0.8**; after rescaling by
  `sqrt((n-1)/(n-deff))`, 0.952 / 0.952 / 0.953. Applies to
  `deff = "kish"`, `deff = "variogram"` and a fixed numeric `deff`.
  **The default `deff = 1` path is bit-identical to before.**

- **[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  ranks model-aware candidates on the standardised deviate, not on
  \|Moran’s I\|.** E\[I\] and Var(I) both depend on the cell count, so
  \|I\| shrinks as k grows whether or not the finer tessellation
  captures anything. Over 300 replicates of a response with **no**
  spatial structure, mean \|I\| fell monotonically from 0.114 at k = 10
  to 0.050 at k = 60 — an \|I\| ranking prefers the largest candidate
  for arithmetic reasons alone. Candidates are now ordered by \|z\|
  using the Cliff & Ord residual moments (exact here, since the
  cell-level residuals are OLS residuals by construction); over the same
  runs z had mean ~0, sd ~1 and a 5% rejection rate of 0.040-0.057 at
  every k. The `"diagnostics"` attribute now carries `moran_z` alongside
  `moran_i`.

- **[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  sweeps four azimuths, not two.** A +/-22.5 degree window around 0 and
  90 covers exactly **90 of the 180 distinct azimuths** — every
  direction between 23 and 67 degrees, and between 113 and 157, fell
  into neither. On simulated fields with 3:1 anisotropy and a true
  major-axis range of 300, the estimate came back at 255 and 249 for
  major axes at 0 and 90 degrees but **151 and 147 at 45 and 135**.
  Since `make_folds(auto_range = TRUE)` sizes blocks from this number, a
  diagonally oriented field silently got blocks half as wide as the
  correlation they were meant to separate. `c(0, 45, 90, 135)` tiles all
  180 azimuths; the same fields now return 255 / 245 / 249 / 228. A
  direction whose variogram never reaches a sill is excluded rather than
  taken as a long range, and the `directional` attribute now has four
  named entries.

- **[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)’s
  calibrated length-scale prior never reached Stan.**
  `brms::set_prior(spec, class = "lscale")` with no `coef` is a *global*
  prior, and brms applies a global prior only to coefficients with no
  individual prior of their own — every `lscale` coefficient always has
  one. brms dropped it with a note and Stan received brms’s defaults,
  which made
  [`gp_lengthscale_bounds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gp_lengthscale_bounds.md),
  the tail calibration and `$info$gp_lscale_prior` dead weight.
  Confirmed with
  [`brms::make_stancode()`](https://paulbuerkner.com/brms/reference/stancode.html):
  the requested prior is absent under the global form and present under
  the coefficient-level form, which is now used. `$info$gp_lscale_prior`
  is read back from
  [`brms::validate_prior()`](https://paulbuerkner.com/brms/reference/validate_prior.html),
  so it records what brms will actually use.

- **The GP basis was sized against the wrong domain measure.** brms
  builds the boundary as `choose_L(x, c) = c * max(1, max(x) - min(x))`
  over the pooled, column-centred covariates — the **full range**, not
  the per-axis half-range in which Riutort-Mayol et al. state their
  inequalities. Recovering the boundary from `make_standata()`’s
  eigenvalues confirms `L = c * full range` exactly at every `c`, so the
  old convention built a boundary **twice as wide** as `gp_k` was sized
  for: the GP was under-resolved, and `$info$gp_ell_min` — the
  diagnostic meant to catch exactly that — was twice too lenient to
  fire. The `c` floor is now brms’s own default 1.25 rather than 1.2.

- **[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) on a
  `gwr_fit` could return a coefficient surface.** The search for
  GWmodel’s fitted-value column matched the whole candidate name vector
  with `%in%` and took the first hit in the *SDF’s* column order — and
  the local coefficients come first. A predictor named `fit`, `pred`,
  `prediction`, `fitted` or `yhat` therefore returned its own
  coefficient column, silently: executed in-sample R^2 was **-1.18**
  against a true 0.981, and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md),
  [`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md)
  and every
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md)
  fold consumed it without a warning. The search now runs in preference
  order and excludes any candidate that is also a model term; all five
  colliding names now give R^2 = 0.981, identical to the renamed
  control.

- **[`coef.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.gwr_fit.md)
  returned GWmodel’s whole SDF data slot** — 15 columns for a
  two-predictor fit, of which 3 are coefficients and the rest are
  standard errors, t-values, the response, the fitted values, residuals
  and `Local_R2`. It now returns the model terms only; reach for
  `object$engine$SDF` for the rest.

- **[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  is reproducible, and no longer disturbs the caller’s RNG.** `seed` now
  defaults to `123L` rather than `NULL`. The `n_max` subsample is an
  internal approximation, not part of the answer, and leaving it
  unseeded made the returned range differ between runs on identical
  input (19531 / 19589 / 19605 on three calls) while silently advancing
  the caller’s stream — and `make_folds(auto_range = TRUE)` sizes its
  blocks from that number. Pass `seed = NULL` for the old behaviour.

- **[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  rejects a non-numeric response.**
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html) on a factor
  returns its level codes, so a factor response produced a variogram of
  an arbitrary integer relabelling of the categories and the estimated
  range changed when the levels were reordered (3700 against 2497 on the
  same data). Factors and character columns are now an error naming the
  column; logicals are read as 0/1.

- **Fold sets built from a different dataset are refused.** Fold splits
  are lists of `..row_id` values, and row IDs are `seq_len(nrow())`
  unless supplied, so passing
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md)
  a `folds` object built from another dataset of the same size applied
  cleanly — every ID matched, every fold was populated, and the model
  was scored on splits describing other observations.
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  now records a small projection-invariant row fingerprint in
  `params$row_probe`, and
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
  and
  [`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
  error rather than proceed. Fold objects from earlier versions carry no
  fingerprint and are passed through unchecked.

- **[`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md)
  rejects duplicated names in `fits`.** `model` is the key
  [`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md)
  joins its metric and Moran’s I tables on, so two fits called `"GWR"`
  produced a 2x2 cross-join: four rows, every one carrying the first
  fit’s numbers, with the second fit never scored at all.

- **[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  rejects a non-numeric predictor.** `gwr.basic()` expands contrasts via
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) and
  fits, but `gwr.predict()` does not and fails, so the model appeared to
  fit and then silently predicted all `NA`.

- **[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  no longer rejects a two-valued continuous response.** The “binary”
  error is now gated on the response being integer-like. A left-censored
  or saturated measurement (every observation at a detection limit or a
  ceiling) has two distinct values and is perfectly continuous; it now
  warns instead. The guard also runs once per fold inside
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
  where a small training fold can legitimately hold only two distinct
  values.

- [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returning the
  wrong length, or nothing, is now an error in
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
  rather than a plausible row count over an all-`NA` comparison.
  [`new_spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md)
  is the documented extension point, so a subclass with a missing or
  mis-sized [`fitted()`](https://rdrr.io/r/stats/fitted.values.html)
  method is user-reachable.

- The cached [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) on
  a `bayesian_fit` is stamped with the `n` and a digest of the data it
  was computed from. The cache environment has reference semantics —
  which is what makes it survive copy-on-modify — so `fit2 <- fit` gave
  both objects the *same* cache, and assigning different data to the
  copy returned the original’s values at the original’s length.

- [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  drops rows with empty or non-finite coordinates, with a logged warning
  naming the count, rather than letting an `EMPTY` POINT reach
  `block_kfold`’s nearest-block rescue and die with “replacement has
  length zero”.

- When every fold fails, the warning now names the first underlying
  error. Previously “all 5 folds failed” was the whole diagnosis even
  when the cause was simply that **brms** or **GWmodel** was not
  installed.

- [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  records the CRS the folds were built in as `params$crs`. Geographic
  input is projected by
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  to a CRS the caller never chose, and `block_size`, `sac_range` and
  `buffer` are lengths in *that* CRS.

- `.morans_i_for_k()` returns `NA` at or below nine cells, where every
  cell neighbours every other and Moran’s I collapses to exactly
  `-1/(k-1)` for any residual vector — a function of the cell count
  alone.

- `residual_morans_i(fit, k = 1)` works on machines without **FNN**.
  [`apply()`](https://rdrr.io/r/base/apply.html) simplified the length-1
  result to a vector, making the neighbour index a 1 x n matrix and
  every row after the first out of bounds.

- **`summarize_by_cell(deff = "kish")` under-estimated the predictor ICC
  by about a factor of `m`.** The pooled one-way ANOVA grouped the `m`
  z-scored predictor columns under the same cell label, so independent
  per-column cell effects averaged away in the shared cell mean and the
  between-cell sum of squares shrank by ~1/m. Measured at true rho = 0.5
  with m = 4: pooled ICC **0.12** against 0.495 per column, so every
  predictor SE was ~44% too small. The pooled group is now (variable,
  cell), which recovers 0.49.

- **Design effects are built from what a column actually observes.** A
  cell of 10 rows with 2 finite responses had its response SE formed at
  the 10-row design effect, then applied to a 2-observation mean: adding
  8 NA-response rows moved the SE from 8.46 to **26.38**. Each column’s
  design effect now uses its own non-missing count, with the mean
  pairwise correlation recomputed over the observed locations when a
  cell has NAs; `cell_weight` is the effective count of the primary
  variable, not of rows (`n` still counts rows).

- **CRS-less coordinates get ONE interpretation, wherever they enter.**
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
  assumed EPSG:4326 for CRS-less data that looked like lon/lat and
  projected it, while every
  [`predict()`](https://rdrr.io/r/stats/predict.html) method passed the
  fit’s CRS as a target — a branch that *stamped* it onto the raw
  numbers. The same rows sat in two different places, and
  `predict(fit, newdata = training rows)` disagreed with `fitted(fit)`
  by up to one response SD (R² 0.98 in-sample, 0.64 via newdata). The
  heuristic is now a single function used by both branches; a fit
  records the assumption it was built under and replays it on CRS-less
  `newdata`. Two further symptoms of the same split — CRS-less
  LINESTRINGs aborting in
  [`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md)
  (“crs not found”), and hex/square
  [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
  refusing input voronoi accepted — are fixed with it. The assumption is
  now announced with a real R warning.

- **[`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)
  refuses a malformed `weights` matrix** instead of silently
  substituting the default k-NN(8) matrix (I = 0.874 returned for four
  malformed shapes against 0.805 for the weights actually supplied).

- **The fold-provenance fingerprint no longer refuses the caller’s own
  data.** Three defects in the version introduced last pass: a character
  `..row_id` was coerced to all-`NA` (and matched row 1 everywhere);
  coordinates were compared as `"%.7g"` strings and flipped on the ~1 in
  5000 that a reprojection moved by 5e-9°; and polygon input was probed
  *after* pointization, so a different `pointize` in the cv call read as
  different data. Both sides now probe the geometry as supplied, keep
  IDs in their own type, and compare numerically within 1e-6°.

- **`n_folds_attempted` counts the folds supplied.** A fold whose test
  rows were all removed as incomplete vanished before fitting and was
  absent from both counts, so five supplied folds reported `4/4`. It is
  announced with a real warning.

- **[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)’s
  elbow uses the signed deviation below the chord.**
  [`abs()`](https://rdrr.io/r/base/MathFun.html) let a concave bump
  *above* the chord — a k where k-means fell into a worse local optimum
  — win with the same magnitude.

- **[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  returns `NA` for a constant variable** (an exactly explained response,
  or a constant one) instead of a fitted “range” of 168 or 673 from a
  variogram that is identically zero. `make_folds(auto_range = TRUE)` no
  longer re-opens the unseeded subsample by forwarding its own
  `seed = NULL`.

- **[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
  attaches a user-supplied *global* `lscale` prior at coefficient
  level**, the same way it does its own, so
  `set_prior(..., class = "lscale")` reaches Stan instead of being
  discarded.

- **[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) on a
  `gwr_fit` returns the prediction when a predictor is named
  `prediction`.** The model-term exclusion added last pass was applied
  to `gwr.predict()`’s SDF too, where coefficients are suffixed `_coef`
  and the column literally named `prediction` *is* the prediction;
  [`predict()`](https://rdrr.io/r/stats/predict.html) returned all `NA`.

- Smaller: a logical response meets the same binary-response guard as
  `0/1`;
  [`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
  errors on a non-numeric response instead of returning `n = 0`;
  [`fitted.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.bayesian_fit.md)
  errors when the posterior cannot be drawn instead of returning silent
  `NA`;
  [`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md)
  refuses a grid above `max_cells` (default 1e6) up front;
  `make_folds(method = "buffered_loo")` states its guard in bytes
  (splits are ~4n² bytes; the old n = 20000 cap admitted 1.6 GB); `-0`
  and `0` are the same coordinate in the duplicate-aware k-NN;
  `.morans_i_for_k()` returns the `NA` pair whenever the moments are
  unavailable.

- **Documented warnings are now R warnings.** Eight paths the manual
  described as warning only wrote a logger line, invisible to
  `tryCatch(warning = )`, `expect_warning()` and `options(warn = 2)`:
  [`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)
  returning `NULL`, the CRS assumption in
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  and stamping in
  [`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md),
  GWR collinearity, seed clamping in
  [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)
  /
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
  dropped rows in
  [`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md),
  and the dropped column in
  [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md).
  Three deliberate methodological cautions (`include_coords = TRUE`,
  `random_kfold` feature selection, non-standardised Moran weights) stay
  logged and their documentation now says so.

- **[`predict()`](https://rdrr.io/r/stats/predict.html) on a Bayesian GP
  fit depended on which other rows shared the call.** brms 2.x stores
  `Xgp`, `dmax` and `cmeans` in a fit’s GP basis but **not** the
  Hilbert-space boundary `L`, so `brms:::.data_gp()` recomputes it from
  whatever rows [`predict()`](https://rdrr.io/r/stats/predict.html) is
  handed. Every eigenfunction of the approximation therefore moved with
  the newdata bounding box while the fitted basis coefficients stayed
  put. Measured: `L` was 5.57 at fit time, 4.02 for a five-row `newdata`
  and 3.63 for one row;
  [`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md)
  on the same 10,000-cell grid differed by **1.78** between
  `chunk_size = 5000` (the documented default) and a single call; and
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
  which predicts each test fold separately, scored **every** fold
  against a basis the model was never fitted with.
  [`predict.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.bayesian_fit.md)
  now pins the boundary by appending the training coordinate extrema and
  dropping them again, so the chunked, fold-wise and single-call answers
  are identical.

- **A predictor whose name is not a syntactic R name fitted a different
  model.** Every backend builds its formula from these names, so a
  column called `"B5-B4"` was fitted as `B5 - B4` — a different model,
  silently, with `$predictor_vars` still reporting `"B5-B4"` — and
  `"band 4"` died inside
  [`str2lang()`](https://rdrr.io/r/base/parse.html) with a parser error
  naming no column.
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
  now refuses such names up front and says which column and what to do
  about it. (Backticking is not a fix: the `sp` coercion GWmodel needs
  runs the names through
  [`make.names()`](https://rdrr.io/r/base/make.names.html) anyway, after
  which the formula and the data disagree.)

- **`summarize_by_cell(deff = "kish")` weighted the response with a
  predictor’s ICC.** The fallback to the predictor ICC is for the case
  where no response was supplied; applying it whenever the response’s
  own ICC came out non-positive meant that merely adding a predictor to
  the call changed the response’s regression weight by a factor of
  **20**, while the response SE was (correctly) left at `deff = 1`.

- **A column literally named `n` was summarised as the row count.**
  [`dplyr::summarise()`](https://dplyr.tidyverse.org/reference/summarise.html)
  makes each new column visible to the expressions after it, so the row
  count shadowed a response or predictor named `n`: `resp_mean_n` came
  back equal to `n`, with `sd` and `se` `NA`, and no message.

- **Variogram models are fitted with a nugget.** A nugget-free model
  forces the curve through the origin, and gstat’s default `N/h^2`
  weights buy that constraint by collapsing the range: with a 50% nugget
  the fitted range came back at about **0.45** of the truth, so
  `make_folds(auto_range = TRUE)` built blocks less than half the
  correlation length it reported.

- **The directional maximum is used only when anisotropy is
  established.** Splitting 180 degrees four ways leaves each directional
  variogram about a quarter of the point pairs, and the maximum of four
  noisy estimates is biased upward: on isotropic simulated fields the
  returned range ran about **40%** above the truth and the “notable
  anisotropy” warning fired on the majority of them. An omnidirectional
  variogram is now fitted alongside and is the default answer; the
  directional maximum is used when all four directions fit, their ratio
  exceeds 1.5, **and** the widest stands more than 1.5x above the
  all-pairs estimate. On the package’s own test field (true range 80)
  the old rule returned 248 and the new one returns 84; genuine 3:1
  anisotropy is still recovered.

- **[`predict()`](https://rdrr.io/r/stats/predict.html) replays the CRS
  decision the fit was made under, including a negative one.** When
  CRS-less training data were passed through as planar, nothing recorded
  that, so every [`predict()`](https://rdrr.io/r/stats/predict.html)
  re-ran the lon/lat heuristic on `newdata` alone — and a **subset** of
  those same training rows, whose own bounding box sits inside the
  lon/lat envelope, was judged differently from the whole: taken for
  degrees, reprojected, and predicted about 1e6 m from where it was
  fitted. `predict(fit, training_subset)` disagreed with
  `fitted(fit)[subset]` by more than the response’s standard deviation
  while `predict(fit, full_data)` agreed exactly.

- **Antimeridian data are no longer flattened onto EPSG:3857.** A
  bounding box cannot distinguish global coverage from a layer
  straddling +/-180 degrees, but the coordinates can (one very large gap
  in the sorted longitudes). Web Mercator *splits* such a layer: two
  stations 41 km apart came out **40,068 km** apart, destroying every
  distance downstream. Only genuinely global coverage now falls back to
  EPSG:3857.

- **The projection for a wide extent is chosen by measurement, not by
  rule of thumb.** Albers standard parallels came from the bounding box
  while the conic/azimuthal branch came from the centroid latitude, so a
  trans-equatorial extent put the parallels either side of the equator:
  at `lat_1 = -lat_2` PROJ refused the string outright and
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md),
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  and
  [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
  all aborted with an internal-looking “invalid crs”, and just short of
  it the projection distorted distances by **15.6%** against the single
  UTM zone’s 1.65%. The candidates are now scored by projecting a sample
  of the data’s own points and comparing planar with geodesic distances
  — WGS84 ellipsoidal distances (Vincenty), not
  [`sf::st_distance()`](https://r-spatial.github.io/sf/reference/geos_measures.html)’s
  s2 sphere, whose 0.24–0.56% gap from the ellipsoid is the size of the
  errors being ranked and mis-ordered the candidates on 16 of 40 random
  wide extents; the message reports both figures.

- **[`fitted.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.gwr_fit.md)
  could return a coefficient surface.** Which GWmodel call produced an
  SDF was inferred from its column names, so a *predictor* ending in
  `_coef` made a `gwr.basic()` SDF look like a `gwr.predict()` one and
  switched off the model-term exclusion; with a second predictor named
  `yhat`, [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) then
  returned that predictor’s local coefficient surface (R2 **-2.28**
  against 0.986) and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
  followed it. The mode is now passed explicitly.

- **A logical predictor supplied as text was silently mis-coded.**
  `is.numeric`, `is.factor` and `is.character` are all `FALSE` for a
  logical, so the type guard skipped it entirely: a `"TRUE"`/`"FALSE"`
  character column — what a CSV round trip produces — was factor-coded
  1/2 against splits built on 0/1, sending every row to the `TRUE` side
  (correlation **0.21** with the correct predictions).

- **[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
  scores every model on identical folds.** The documented guarantee only
  held when `folds` was supplied: with `folds = NULL` each backend built
  its own from `k`, `seed`, `block_size` and the rest, so a per-model
  `block_size` put **104 of 150 rows** in different folds for GWR and
  RF, and `seed = NULL` did the same with no overrides at all. The folds
  are now built once, and the arguments that decide the split are
  protected.

- **User weights with a non-zero diagonal.** Every moment of Moran’s I
  assumes no observation is its own neighbour. A self-inclusive
  row-standardised kNN matrix — an easy thing to build by hand —
  rejected the null on **75%** of white-noise residuals at a nominal 5%,
  with no condition raised. The diagonal is now zeroed with a warning.

- **`build_tessellation()$index` no longer snaps outside points to the
  nearest cell.** Points outside every cell were assigned to whichever
  was closest, so a summary built from the index counted all 40 points
  of a layer whose study area held 10, while
  [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
  on the same cells correctly reported 30 misses. Only a point within a
  thousandth of a cell width is snapped now; the rest are `NA`, and the
  count is logged.

- **[`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md)
  is stable across projections again.** The sort key was the raw double
  centroid, so cells sharing an exact `x` in one CRS differed by ~1e-11
  degrees after a round trip through another: **14 of 16** cells got a
  different ID depending on which CRS the layer arrived in — the exact
  failure the function exists to prevent. The key is now rounded to 7
  decimals (about a centimetre).

- **A numeric `deff` is applied as the uniform inflation it is
  documented to be.** The `E[s^2]` correction the estimated design
  effects apply is derived from within-cell correlation and is
  unjustified for a constant the caller chose: `deff = 2` doubled a
  3-point cell’s SE instead of multiplying it by `sqrt(2)`, and returned
  `NA` for every cell with `n <= deff`.

- **The Kish ICC guard matches its documentation.** The docs promise “at
  least 2 cells with 2+ observations; falls back to `deff = 1`
  otherwise”, and the check was only `k >= 2, N >= 4`: all-singleton
  cells give a within-group sum of squares of 0 and therefore an ICC of
  exactly **1** — a design effect of `n` from data carrying no
  within-cell information at all.

- **[`gp_lengthscale_bounds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gp_lengthscale_bounds.md)
  and the GP basis use two dimensions and unique locations.**
  [`stats::dist()`](https://rdrr.io/r/stats/dist.html) uses every
  column, so POINT Z geometry gave 3-D length-scale bounds in mixed
  units; and
  [`brms::gp()`](https://paulbuerkner.com/brms/reference/gp.html)
  defaults to `gr = TRUE` and reduces its covariates to unique rows
  before taking the boundary, so replicated locations made the package’s
  `S` — and with it `gp_c` and the basis-adequacy threshold `gp_ell_min`
  — wrong by that factor (0.68x with one heavily-sampled station).

- **[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  and
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  drop the Z dimension.**
  [`gstat::variogram()`](https://r-spatial.github.io/gstat/reference/variogram.html)
  and
  [`sf::st_distance()`](https://r-spatial.github.io/sf/reference/geos_measures.html)
  use every coordinate dimension, so an XYZ layer had its elevation
  folded into each lag: the returned “range” was a length in 3-D (413.6
  against 136.8 for the same stations) while the block grid, the
  buffered-LOO buffer, NNDM’s neighbour distances and
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  all work in 2-D map distance.

- **[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  returns the elbow first.** Under the geometric criterion — the
  default, and the fallback every model-aware call takes below the
  nine-cell floor — the candidates came back sorted ascending, so the
  knee sat in the middle: `k[1]` and `top_n = 1`, both documented as
  “the top-ranked candidate”, returned **knee − 1** on every such call,
  and the help example answered `1` for two clearly separated clusters.
  The vector is now knee first, then its lower and upper neighbours, so
  position 1 means the same thing on both paths. The quick-start data’s
  answer moves from `3 4 5` to `4 3 5`. Same code in 1.0.0.

- **[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  is invariant to rotating the layer.** Two things were not. The lag
  cutoff was a fraction of the bounding-box *diagonal*, a property of
  the axes rather than of the points, which grows by up to √2 when the
  same layer is rotated 45°; it is now a fraction of the farthest pair
  (found on the convex hull), as the documentation always said. And the
  directional maximum is **never** preferred over a usable all-pairs
  fit. The 1.0.0 code returned the widest of two axis-aligned
  directions; the three hurdles later put in front of it (all four
  directions fitted, ratio above 1.5, maximum above 1.5× the all-pairs
  fit) still let the noise through — on one isotropic field rotated in
  10° steps they “established” anisotropy in **14 of 18** orientations
  and returned ranges from 225 to 529, a 2.35× spread produced by
  nothing but the direction the axes pointed. The four directional
  ranges are still reported, as a diagnostic, in the `directional`
  attribute; `anisotropy_used` is `TRUE` only when the all-pairs fit
  itself failed. A field *known* to be anisotropic should size its
  blocks from `max(attr(range, "directional"))` explicitly, and the log
  line says so.

  The variogram fit itself no longer depends on gstat’s single starting
  value. `fit.variogram()` starts the optimiser at a third of the
  longest lag; for a field whose range is a small fraction of the extent
  that start is ten times too long, and whether the iteration landed or
  collapsed to a singular model depended on floating-point details — the
  same 250-point field fitted on Linux and came back singular on an
  arm64 Mac, where the estimate then fell through to the directional
  maximum (315 against a true range of 80). Shorter and longer starting
  ranges are tried as well, the converged, non-singular fit with the
  smallest weighted sum of squares (gstat’s own criterion) is kept, and
  a converged spherical fit is preferred over an exponential one that
  did not converge. When no model fits at all — a flat, nugget-only
  variogram — the `NA` now carries the empirical variogram, so
  `plot(fit, type = "variogram")` draws it and says why there is no
  range instead of refusing.

- **[`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)
  splits tied neighbour distances.** The k-nearest neighbour weights
  broke ties at the k-th distance by whichever point came first — in row
  order on the dense path, in kd-tree order under **FNN** — so the same
  data in a different row order gave a different I, z and p whenever
  observations shared a location (repeat visits), and on gridded data
  the answer depended on whether FNN was installed. Points tied at the
  k-th distance now share the remaining weight equally, on both paths,
  which is the only rule that is a function of the geometry alone; the
  matrix is row-standardised as before. Where no distances tie the
  weights equal `spdep`’s k-NN weights exactly; on repeat-visit data
  every co-located twin is retained with its share, and the statistic no
  longer changes when the rows are shuffled or when **FNN** is
  installed. Same code in 1.0.0.

- \*\*[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)‘s
  collinearity diagnostic is the scaled condition index, thresholded at
  30.\*\* It was [`kappa()`](https://rdrr.io/r/base/kappa.html) on the
  raw, unscaled predictor matrix at 1e6 — a number that depends on the
  predictors’ units, so rescaling a column changed it and the threshold
  was a threshold on nothing: a design with a scaled condition index of
  1322, whose local coefficients ran from −86 to +150 around a true
  value of 2, raised nothing. Belsley’s index (every column scaled to
  unit length, intercept included, the ratio of the largest to the
  smallest singular value) is now computed for the global design and for
  each sampled local window, and the conventional 30 is the threshold in
  both. Expect the warning on designs that were silent before.

#### API and default changes

- [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  returns `NA` instead of a number when the fitted range runs past the
  longest observed lag, or the optimiser stopped at its iteration limit.
  Such a range is *unidentified*, not long — the empirical variogram
  never reached a sill — and sizing blocks or design effects with it is
  worse than declining to. The refusal carries `rejected_range` and
  `rejected_reason` attributes and keeps the fitted variogram, so
  `plot(type = "variogram")` can draw exactly the case worth looking at;
  the result is classed `sac_range`, so it prints as a bare `NA` rather
  than dumping the fit to the console. Callers that fed the old number
  straight into `make_folds(block_size = )` now need to handle `NA` —
  that is the point.

- Removed the legacy wrappers `evaluate_models()`,
  `evaluate_models_cv()` and `phi_prior_bounds()`. Use
  [`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md),
  [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
  and
  [`gp_lengthscale_bounds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gp_lengthscale_bounds.md).

- [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
  gains an `"RF"` branch and an `rf_args` argument, so a `ranger` forest
  can be compared against GWR and the Bayesian GP on identical folds.
  Unrecognised model names now raise a warning and are dropped, and a
  request with nothing recognised left is an error. Previously a bare
  [`intersect()`](https://generics.r-lib.org/reference/setops.html)
  discarded anything outside `c("GWR", "Bayesian")` and fell back to
  GWR, so `models = "RF"` silently ran GWR and reported it as the
  answer.

- [`coef()`](https://rdrr.io/r/stats/coef.html) on a `spatial_fit` now
  either returns coefficients or errors; it never returns `NULL`.
  [`coef.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.gwr_fit.md)
  and
  [`coef.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.bayesian_fit.md)
  used to return `NULL` on failure, which is indistinguishable from
  “this model has no fixed effects”, so `lapply(fits, coef)` quietly
  produced a short answer.
  [`coef.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.rf_fit.md)
  errors as before — a forest has no coefficients; use
  `fit$info$importance`. See
  [`?new_spatial_fit`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md).

- Every [`predict()`](https://rdrr.io/r/stats/predict.html) method
  errors when the number of predictions does not match the number of
  rows that survived cleaning. It used to recycle silently: two
  predictions for four clean rows produced a four-row answer.

- [`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md)’s
  default `type` is now `"square"`, matching
  [`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md).
  The two disagreed, so cached and uncached calls built different grids
  from the same arguments.

- [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
  keeps `cell_id` on `"hex"` and `"square"` grids instead of deleting it
  after indexing, so all four methods return the same ID column and
  `plot_tessellation_map(fill_col = "cell_id")` works on a grid.
  `poly_id` is retained alongside it. `params$expand` now echoes the
  value passed rather than a hard-coded `0`; the grid and triangle
  methods still ignore `expand`, which is now documented.

- [`clip_target_for()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clip_target_for.md)
  projects lon/lat input before applying `expand`. The
  fraction-of-extent form was computed from a bounding box in degrees
  and handed to
  [`sf::st_buffer()`](https://r-spatial.github.io/sf/reference/geos_unary.html),
  which reads `dist` as metres. The returned clip target is therefore in
  the projected CRS, not the input CRS, and says so.

- Voronoi, grid and triangle clipping union the boundary first. Against
  a multi-feature boundary,
  [`st_intersection()`](https://r-spatial.github.io/sf/reference/geos_binary_ops.html)
  split every straddling cell into one row per boundary feature and
  grafted the boundary’s attribute columns onto the result.

- `predict.bayesian_fit(newdata = NULL)` honours `summary`, `type` and
  `draws`. It short-circuited to
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html), which caches
  epred column means and nothing else, so `summary = "median"` returned
  means and `type = "predict"` returned expected values — silently.

- [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
  accepts `predictor_vars = character(0)`, making intercept-only spatial
  GP models reachable from
  [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md).
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  and
  [`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
  reject an empty set explicitly.

- `fit_bayesian_spatial_model(control = )` is *merged* over the package
  defaults (`adapt_delta = 0.9`, `max_treedepth = 12`) rather than
  replacing them. Passing `list(max_treedepth = 15)` silently dropped
  `adapt_delta` — the setting the divergence warning tells you to raise.

- [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)’s
  `predictive_coverage` is averaged across folds **weighted by each
  fold’s `n_pred`**. The per-fold values are means over that fold’s test
  rows, so an unweighted average is not the pooled quantity once fold
  sizes differ — and `block_kfold` tolerates a 3:1 imbalance before it
  even logs a warning.

- The three seeding functions
  ([`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
  [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md),
  [`voronoi_seeds_random()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_random.md))
  all emit `seed_id` and `method` columns, so they are drop-in
  interchangeable.
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
  returns seeds in the boundary’s CRS for every method; per-branch
  alignment had made the final alignment block unreachable.

- Geometry-type checks require **every** geometry to be an accepted
  type, not merely one of them. A mixed POINT/POLYGON layer passed a
  POINT-only check.

- [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  errors on a `target_crs` that does not resolve to a usable CRS, rather
  than returning the input unchanged and letting unprojected coordinates
  flow into distance and area computations.

- [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
  derives the GP basis count (`gp_k`) and boundary factor (`gp_c`) from
  the ratio of the estimated length-scale to the domain size, rather
  than from the number of observations.

  [`brms::gp()`](https://paulbuerkner.com/brms/reference/gp.html) builds
  a full tensor grid, so `gp(..x, ..y, k = gp_k)` carries `gp_k^2` basis
  functions — `gp_k` is the count *per dimension*. The previous rule
  reduced to `max(15, floor(sqrt(n)))` for any `n` above 45, making
  `gp_k^2` identically `n`: at n = 10,000 the model carried 10,000 basis
  functions and an n × n design matrix, at which point the approximation
  was no longer approximating anything.

  Across the scenarios in `dev/baseline-structural.rds` the derived
  value is 22–24 per dimension and largely independent of `n`: at n =
  2,000 the basis count falls from 1,936 to 576, at n = 10,000 from
  10,000 to 529, and at n = 200 it *rises*, 15 to 23 — a correction, not
  an optimisation, so results move in both directions. `gp_c` was
  hard-coded at 1.5, too small whenever the length-scale exceeds roughly
  half the domain half-range; the derived value ranges 2.85–3.59 over
  the same scenarios. Pass `gp_k` and `gp_c` explicitly to restore the
  old behaviour. For scale, a 4-fold cross-validated fit at n = 2,000
  with 2 chains × 1,000 iterations takes 1,186 s on the reference
  machine at cross-validated R² 0.927 (`dev/baseline-accuracy.rds`,
  2026-08-20); no comparable timing was captured for 1.0.0, so none is
  quoted.

- The GP term is built with `scale = FALSE`, changing the default result
  of every Bayesian fit.
  [`brms::gp()`](https://paulbuerkner.com/brms/reference/gp.html)
  otherwise rescales its covariates so the maximum pairwise distance is
  1 and reports `lscale` in that space, while this package standardises
  the coordinates itself and expresses the length-scale prior, `gp_c`
  and the basis adequacy threshold in those units. The two
  normalisations differed by roughly the maximum pairwise distance (~4.9
  for standardised 2D coordinates), leaving the automatic prior about
  five times too diffuse — a likely contributor to divergent transitions
  and rejected initial values. There is now exactly one coordinate
  scaling.

- The GP fits one length-scale per coordinate axis (`gp_iso = FALSE`), a
  second change to default Bayesian results. Coordinates are
  standardised per axis, so a single shared length-scale made the kernel
  anisotropic in the original CRS by the ratio `sd(X)/sd(Y)` — a
  property of the sampling layout, not of the process. Pass
  `gp_iso = TRUE` for the previous behaviour; cost is unchanged, since
  the tensor grid is `gp_k^2` either way.

- The automatic GP length-scale prior is a calibrated inverse-gamma
  rather than `normal(0, sd)`, a third change to default Bayesian
  results. A half-normal on a positive parameter puts its mode at zero,
  so most of its mass sat at length-scales shorter than the basis can
  resolve — where the Hilbert-space approximation develops a funnel and
  the sampler diverges. The replacement pins 1% of its mass below the
  estimated lower bound and 1% above the upper. The two tail conditions
  have one exact solution — a one-dimensional root in the shape — and
  that is how it is found, so the calibration succeeds for any bounds
  with `upper > lower`; degenerate bounds fall back to the half-normal
  with a logged note. The prior applied is recorded in
  `$info$gp_lscale_prior`.

- [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  no longer forces continental-extent data into a single UTM zone.
  Transverse Mercator scale error grows quadratically with distance from
  the central meridian, so data spanning the contiguous United States
  carried distance errors of roughly 7.5% near the extent edge,
  propagating silently into
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
  `make_folds(block_kfold)` block sizing, GWR bandwidth selection and
  the GP length-scale. Extents reaching more than 5° from the candidate
  zone’s central meridian now receive an equal-area projection centred
  on the data, with a logged explanation. Only longitude offset triggers
  the switch, since `cos(lat)` shrinks the distance from the central
  meridian and a tall narrow north-south extent is UTM’s design case.
  (Which equal-area projection is no longer decided by latitude band —
  see *The projection for a wide extent is chosen by measurement* below,
  which also replaced the EPSG:3857 fallback for wide bounding boxes
  with antimeridian detection.) Pass `target_crs` to override.

- Core counts follow the session’s `mc.cores` opt-in, and are capped.
  [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)’s
  documented default was
  `cores = max(1L, parallel::detectCores() - 1L)`, and
  `cv_*(parallel = TRUE)` auto-detected the same way with no cap — 63
  workers on a 64-core host, and a hard error wherever
  `_R_CHECK_LIMIT_CORES_` is set. The Bayesian default is now
  `getOption("mc.cores", 1L)`; the auto-detect path is capped by that
  option when it is set; and every worker count, explicit or not, is
  capped at the machine’s core count (with a message) and at two under
  `R CMD check`.

#### Guards, messages and stricter input handling

- [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
  now drops rows whose geometry is empty or whose coordinates are not
  finite, and counts them in its existing log line.
  [`st_geometry_type()`](https://r-spatial.github.io/sf/reference/st_geometry_type.html)
  calls an EMPTY POINT a “POINT”, so nothing ever looked at the
  coordinates: such a row reached GWmodel as a raw `sp` coercion error
  naming no row, `ranger(include_coords = TRUE)` as “Missing data in
  columns”, and brms at predict time as an infinite GP boundary that
  made **every** basis function on the surface `NaN`. It also refuses a
  response listed among its own predictors, which was leakage in the
  forest, a silently reduced model in GWR, and duplicated rows plus a
  phantom `<none>` entry in the GWR selection table.

- [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  validates `k`. A non-integer `k` truncated inside
  `rep(floor(n/k), k)`: only `floor(k)` folds were built, the last rows
  of the permutation landed in **no** test set, and `k` was echoed back
  unchanged, so `length(folds) != k`.

- `make_folds(block_kfold)` refuses a grid it cannot build. There was no
  cap on `nx * ny`, so `block_size` in the wrong unit asked for 1e8-1e11
  cells and exhausted memory; the message now names the implied cell
  count, the extent and the CRS units, as
  [`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md)
  already did.
  [`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md)
  gained the same guard on `cell_size` and `n_cells`.

- Hand-built `folds` are checked. Train and test overlapping is not
  cross-validation — the model is fitted and scored on the same rows and
  the result is reported as a CV score (RMSE 0.50 against 0.97 for the
  same data properly split) — and is now refused; fold IDs that name no
  row were dropped silently by `na.omit(match())` and are now counted
  and logged.

- The fold provenance probe tolerates a row it could not measure. It was
  taken before the empty-geometry filter, so it carried `NaN`
  coordinates and every later `cv_*()` call on the same data died with
  R’s internal “missing value where TRUE/FALSE needed”.

- [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  drops rows with unusable coordinates instead of returning `NA` for the
  whole layer under the message “variogram model fit failed”, which
  blamed the fit rather than the row and disagreed with
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  on the same data.

- GWR’s local collinearity spot-check now runs **after** the bandwidth
  is chosen, includes the intercept column, and treats a non-finite
  condition number as extreme. On the default path (`bandwidth = NULL`)
  it used a stand-in window, so the documented warning never fired;
  without the intercept it could not see an indicator constant inside a
  window; and `is.finite(cn) && cn > 1e6` discarded exactly singular
  designs. A new post-fit warning counts local regressions that returned
  non-finite coefficients — previously
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
  silently reported metrics computed from the survivors (n = 18 of 200,
  R2 = 0.96).

- [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  warns when a fixed bandwidth is implausibly small for the data’s
  extent. The bandwidth is a distance in the CRS the fit runs in, which
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
  may have chosen: 0.2 supplied for lon/lat data is 0.2 metres, and
  every local window came back empty with nothing raised. The argument’s
  documentation now says so.

- [`summary.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.spatial_fit.md)
  applies the same response-type guard as
  [`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md).
  A character response still produced `n = 0` and all-`NA` metrics
  there, and a factor died inside
  [`abs()`](https://rdrr.io/r/base/MathFun.html).

- [`print.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.spatial_fit.md)
  prints the CRS (which its documentation has always promised) and one
  `Formula` line. [`sprintf()`](https://rdrr.io/r/base/sprintf.html)
  vectorised over a multi-element
  [`deparse()`](https://rdrr.io/r/base/deparse.html), so every
  `bayesian_fit` printed its formula as two mangled fields — the `gp()`
  term this package builds is always long.

- [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  warns when `cells_sf` carries duplicated IDs (each summary row is then
  repeated per matching cell, so `sum(n)` exceeds the number of points),
  and realigns **every** per-cell vector in `deff_applied` after the
  join, not only `deff`: `$rbar` was left in pre-join order, so
  `deff[i]` and `rbar[i]` described different cells.

- [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
  warns when no feature falls in any polygon instead of returning an
  empty layer silently.

- `create_grid_polygons(target_cells = )` builds square cells for
  `type = "square"`. `cellsize = c(w/nx, h/ny)` forced an exact bbox
  tiling, so the cells were rectangles (aspect 10.5 on a 1000:1 strip).
  For `type = "hex"` a differing `cellsize[2]` is now collapsed with a
  warning before the `max_cells` estimate, which used both components
  and was therefore off by their ratio.

- [`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md)’s
  renumbering is documented: it applies
  [`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md)
  and
  [`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md)
  does not, so the same cell carries a different `poly_id` from the two
  builders.

- `get_voronoi_seeds(method = "kmeans")` clusters in two dimensions and
  drops rows with unusable coordinates, matching
  [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md).
  It used every column of
  [`st_coordinates()`](https://r-spatial.github.io/sf/reference/st_coordinates.html),
  so POINT Z geometry was clustered in 3-D with elevation dominating,
  and an EMPTY POINT crashed inside
  [`kmeans()`](https://rdrr.io/r/stats/kmeans.html).

- [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  refuses a factor or character response, which
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html) silently turned
  into level codes: re-ordering the levels of the same factor changed
  the chosen number of levels and every `moran_z`.

- [`plot_tessellation_map()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_tessellation_map.md)
  brings CRS-less layers into the plot’s CRS instead of passing them
  through to fail inside
  [`ggplot_build()`](https://ggplot2.tidyverse.org/reference/ggplot_build.html)
  at print time with sf’s message, naming no layer; and it borrows a CRS
  from an overlay when the tessellation itself has none.

- Documentation corrected where it did not match behaviour:
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  states how it chooses a projection for a wide extent (by centroid
  latitude when that entry was written, by measured distortion in this
  release) and that the choice is announced rather than silent;
  `.looks_like_lonlat()`’s two tests are a disjunction and the extent
  test decides first, so a small planar survey inside the lon/lat
  envelope IS taken for degrees — the trade and its reasoning are now
  stated;
  [`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md)
  no longer claims to match
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  while doing something else;
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  states which estimand its standard errors are for (the grand mean,
  where measured coverage is 0.95, not the cell’s own mean, where the
  naive SE is the better estimate) and that the variogram path applies
  one correlation function to every column.

- **A misspelt `newdata` is an error, not an in-sample answer.**
  [`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md),
  [`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md)
  and
  [`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md)
  forward `...` to [`predict()`](https://rdrr.io/r/stats/predict.html),
  which checks it only on the out-of-sample branch, so
  `model_metrics(fit, newdta = hold)` silently took the in-sample branch
  and returned an RMSE of **1.086** where the held-out answer was
  **25.24**, with the same return shape. Arguments in `...` with no
  `newdata` are now refused by name.
  [`predict()`](https://rdrr.io/r/stats/predict.html) on a `gwr_fit` or
  a `bayesian_fit` likewise refuses unknown arguments instead of
  swallowing them;
  [`predict.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.rf_fit.md)
  accepts only ranger’s own predict arguments through `...`.

- **[`predict()`](https://rdrr.io/r/stats/predict.html) enforces one
  `newdata` contract on all three fit classes.** A bare `sfc` died
  inside two of them with R’s “argument must be coercible to
  non-negative integer”; a numeric-at-fit predictor that arrived as
  character (a CSV round-trip) was refused by name by `rf_fit` and
  `bayesian_fit` and returned all-`NA` with a generic backend warning
  from `gwr_fit`; a missing column was reported by two different
  functions in two wordings. All three now run the same check first, so
  the message is the same whichever fit is behind it.

- [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  validates `block_size` the way it validates `k`: a single finite
  positive number, else an error naming the argument. `NA` and a
  length-2 vector used to die as internal R errors, and a negative, zero
  or character value was silently ignored — yet echoed back in
  `params$block_size` as if it had been used. The grid-size guard also
  formats its own message: a `block_size` in the wrong unit could ask
  for a grid past 2³¹ cells on a side, which `%d` refused with “invalid
  format” instead of the documented refusal.

- [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  says which variable it modelled. A predictor name absent from the data
  was dropped silently by
  [`intersect()`](https://generics.r-lib.org/reference/setops.html), and
  with every name unknown the **raw** response was modelled — so
  `make_folds(auto_range = TRUE)` sized blocks from a range of 88.5
  instead of the residual range of 362.2 with nothing said. An unknown
  predictor is now an error, as it is everywhere else; a detrending fit
  that fails (an all-`NA` column) raises a warning and falls back to the
  raw response; and a new `detrended` attribute records which was used.

- `make_folds(block_kfold, boundary = )` refuses a boundary containing
  none of the points — almost always two layers in different places, a
  CRS that could only be stamped — and raises a warning counting the
  points that fall outside a boundary that contains some, since the
  region is silently extended to cover them. The single-block error
  names what produced the grid (the block size, `block_nx`/`block_ny`,
  or the automatic grid) rather than always blaming
  `block_nx`/`block_ny`, which the caller may never have passed.

- Rows that no fold names are reported. The `folds` ↔︎ data guard was
  one-directional: fold IDs naming no row were counted and dropped, but
  rows in the data that appear in no fold’s train or test set passed
  with no condition at any level — a `folds` object built on
  `site[1:45, ]` and applied to all 90 rows scored 45 of them and
  reported `n_folds_attempted = n_folds_succeeded = 3`. Every `cv_*()`
  now raises a warning with the count and an example row ID.

- User-facing conditions name the function the user called. The cross-
  validation path leaked two internal names into ordinary console output
  — `.remap_folds():` and `.cv_run_folds():` — and
  [`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md),
  a wrapper around
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
  reported every message, warning and error in
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)’s
  name. The `sf`-input assertion shared by the tessellation and seeding
  functions said “Expected an sf object” with no function named; it now
  names the caller and, when handed a whole
  [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
  result, says to pass its `$cells`.
  `coerce_to_points(mode = "line_midpoint")`’s MULTILINESTRING refusal
  is prefixed like every other.

- `clear_grid_cache(cache_env = )` removes only its own entries. It
  removed every binding in the environment it was handed and counted
  them all as “entries removed”, so a user who passed a project
  environment lost unrelated objects. Cache keys now carry a
  `spatialkit_grid::` prefix and nothing else is touched.

- `summarize_by_cell(deff = "variogram")` now uses **every** structured
  component of a nested variogram model, each weighted by its partial
  sill, which is the correlation the model implies
  (`1 - gamma(h) / sill`). It read the single largest component, so a
  user-built `Nug + Exp + Sph` model gave a correlation of 0.108 at 200
  m where
  [`gstat::variogramLine()`](https://r-spatial.github.io/gstat/reference/variogramLine.html)
  implies 0.197, and the design effects and standard errors with it.
  Models from
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  are single-component and are unaffected. A model of a family the
  function does not implement (Matern, power, circular, …) was silently
  read as exponential; it now falls back to `deff = 1` with a warning
  that names the family. Exponential, spherical and Gaussian are
  supported.

- `cv_bayes()$predictions$yhat_sd` is the posterior predictive standard
  deviation of each held-out row, from the same draws that give the
  coverage columns. It was an unconditional `NA` placeholder. It stays
  `NA` when `compute_pred_intervals = FALSE` or the draws failed for a
  fold, and it is documented.

- [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  refuses `n <= p + 1` observations with one error, before touching the
  backend. Such a fit cannot have a residual degree of freedom in any
  window; `n = 2` used to warn three times (“only 2 observations”,
  “fallback bandwidth”, “2 of 2 local regressions singular”) and return
  a fit whose fitted values were all `NA`.

- [`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md)
  gains `max_entries` (default 50): once the cache holds that many
  grids, adding one evicts the earliest-added. It never evicted, so a
  loop over a thousand boundaries held every grid (about 2 MB per 2,500
  cells) for the life of the session. Nothing but the grids is written
  into a caller-supplied `cache_env`; the insertion order lives inside
  the package. The cache key also hashes the package version, so a cache
  that outlives an upgrade cannot serve a grid built by an older
  [`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md).

- [`spatialkit_quiet()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)
  accepts a `logger` threshold as well as `TRUE`/`FALSE`, and the value
  it returns can be passed back:
  `old <- spatialkit_quiet(); spatialkit_quiet(old)` restores exactly
  the level that was in force. `spatialkit_quiet(FALSE)` put back the
  package default (WARN) whatever had been set, and the returned value
  was refused as “must be TRUE or FALSE”.

- `.onUnload()` disarms both `logger` appenders, so the `spatialkit`
  logger namespace no longer keeps pointing at the session’s temp-file
  path after `unloadNamespace("spatialkit")`. `.onLoad()` re-registers
  them.

- [`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)
  documents the second reason `"residual"` (and `"auto"`) falls back to
  the randomisation null: fewer than four residual degrees of freedom,
  where the residual variance formula divides by `(n - p)(n - p + 2)`.
  `"residual"` logs a warning when it does, `"auto"` does not, and `df`
  is then `n - 1`.

- [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  states what “effective range” is for each model – three times the
  range parameter for the exponential fit (95% of the sill), the range
  parameter itself for the spherical fallback (100%) – and its
  return-value documentation matches the code: the no-model case (both
  fits singular) returns the classed `NA` with `rejected_reason` set,
  and only the cannot-even-start cases (no `gstat`, too few values, no
  variance, a degenerate extent) return a bare `NA`. Its example now
  demonstrates a fitted range on a simulated field with a known one (3 x
  100 = 300), and a refusal on a field whose range the data cannot pin
  down.

- [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  documents what `block_multiplier` does (the automatic grid aims for
  `block_multiplier * k` blocks, so each fold holds out about that many;
  3 is a compromise, not a published constant), cites Roberts et al.

  2017. and `blockCV` (Valavi et al. 2019) for sizing blocks from the
        autocorrelation range, and notes that `blockCV` takes the fitted
        variogram’s range *parameter* where `auto_range` takes the
        *effective* range – three times that parameter for an
        exponential fit, so larger blocks. `phi` for `method = "nndm"`
        is explained as Mila et al. (2022) define it: the
        autocorrelation range beyond which matching is unnecessary,
        which
        [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
        supplies.

- `summarize_by_cell(deff = "kish")` says which ICC estimator it is
  (ANOVA with Donner’s `n0`, not REML) and how far the two can differ on
  an unbalanced draw, so the difference is not read as a defect;
  [`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)’s
  training-DI sentence now says what the code does (each fold’s actual
  training rows, not “everything outside the fold”);
  [`?spatialkit`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-package.md)
  no longer lists [`coef()`](https://rdrr.io/r/stats/coef.html) among
  the methods all three backends share (a forest has none); the internal
  elbow helper no longer claims to match Kneedle, which it does not on
  shouldered curves; the GP basis diagnostic’s code comment attributes
  its 10% posterior-mass trigger to this package rather than to
  Riutort-Mayol et al. (2023).

- [`summary()`](https://rdrr.io/r/base/summary.html) on a fit prints
  `R^2` and `Adj R^2` in ASCII with aligned labels (the superscript two
  rendered as `R<U+00B2>` on non-UTF-8 consoles, and `Adj R²=` had no
  space); [`print()`](https://rdrr.io/r/base/print.html) on a random
  forest likewise; the six console messages that carried an em dash use
  `--`. [`coef()`](https://rdrr.io/r/stats/coef.html) on an `rf_fit`
  prefixes its error `coef.rf_fit():` like its siblings.

- `get_voronoi_seeds(method = "kmeans")` sizes its candidate cloud in
  double precision; `50L * as.integer(n)` overflowed to `NA` above
  42,949,672 seeds.

- README: the opening leakage example says it needs `ranger`; the
  no-viable-models example is assigned so that it does not print every
  fold table on a machine that has the backends; the `auto_range` and
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
  failure examples show every line the console actually prints; the
  installation table no longer suggests installing `loo` separately
  (`brms` installs and calls it); the “logged note” from
  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  is identified as an INFO-level line in the session log file, not
  console output; the roxygen2 sentence no longer names a version or a
  `RoxygenNote` field.

- The memory-guard test for `make_folds(block_kfold)` stubs out
  [`sf::st_make_grid()`](https://r-spatial.github.io/sf/reference/st_make_grid.html)
  for the two refused calls, so a regression of the guard fails fast by
  name instead of attempting an 8000 x 8000 grid.

- Every `cv_*()` refuses a fold list whose splits carry no
  `train`/`test` element, with an error that names the problem. It read
  `f$train` / `f$test` straight, got `NULL` for both, and built empty
  folds – so the row-coverage warning fired and blamed folds “built on a
  different or subsetted layer”, which was not the cause, and the run
  returned an all-`NA` `overall` with `n_folds_succeeded = 0`.
  [`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
  had always refused the same input by name; the two now agree. When the
  splits look positional (two unnamed vectors each) the error says so
  and points at the fold label vector instead.

- The `folds` argument of every `cv_*()` documents all three shapes it
  has accepted since the label vector was added earlier in this pass – a
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  result, a list of `list(train =, test =)` splits, or a vector of fold
  labels – where the help listed only the first two. The label vector is
  what makes folds from another package usable directly:
  `blockCV::cv_spatial()` returns one as `$folds_ids` (its `$folds_list`
  holds two unnamed vectors per fold, which is the shape now refused by
  name).

- `MAPE` and `SMAPE` are documented as what they are: averages over the
  rows whose denominator is non-zero. Both have a denominator that can
  vanish – `MAPE` where the observation is zero, `SMAPE` where
  observation and prediction are both zero – and each drops those rows
  rather than returning `Inf`, which is the right arithmetic but was
  reported nowhere. On a response taking exact zeros (counts, rainfall,
  abundance) the consequence is material: with 62 zeros out of 120,
  `MAPE` is an average over 58 rows presented as though it covered 120,
  and `SMAPE` drops precisely the rows a well-fitted model got right, so
  it reads worse than the fit deserves. The new “Percentage errors on
  responses with zeros” section on
  [`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
  – inherited by
  [`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md),
  [`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md),
  [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md),
  [`summary()`](https://rdrr.io/r/base/summary.html) and all four
  `cv_*()` – says so, notes that the `n` column is the finite-pair count
  and not the row count either percentage error used, and points at
  RMSE/MAE/R-squared (and, for a Bayesian fit, CRPS and interval
  coverage) as the metrics unaffected by it. No computed value changes;
  returning the per-metric row count would alter the metric frame’s
  column set and is deferred.

### Bug fixes

- Data carrying **no CRS** works again throughout.
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  now rejects a `target_crs` that does not resolve to a usable CRS
  (previously a typo silently made the call a no-op), but internal
  callers derive that target from another object —
  `st_crs(training_data)` — and that object is allowed to have no CRS.
  Passing `NA_crs_` through turned every CRS-less workflow into a hard
  error: [`predict()`](https://rdrr.io/r/stats/predict.html) on all
  three backends, `make_folds(method = "nndm")`,
  `make_folds(boundary = )`, `prep_model_data(boundary = )` and
  [`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md).
  [`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
  was worse than an error — the per-fold
  [`predict()`](https://rdrr.io/r/stats/predict.html) threw, so every
  fold “failed” and `$overall` came back with `n_pred = 0` and
  `RMSE = NA` behind a generic warning. Internal call sites now pass
  `NULL` (“choose one automatically”) when the source has no CRS; the
  user-facing validation is unchanged.

- `build_tessellation(crs = )`, `create_voronoi_polygons(crs = )` and
  `create_grid_polygons(crs = )` errored with sf’s “cannot transform sfc
  object with missing crs” whenever the input had no CRS — exactly the
  users most likely to pass `crs =`. Reprojection is impossible there,
  but assumption is not: the target CRS is now stamped on with a loud
  warning, matching what
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  already documents. Input that *does* carry a CRS is still reprojected,
  not relabelled.

- [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
  built its argument list with `c(list(...), rf_args)`, so any
  `gwr_args`/`rf_args` entry whose name collided with one the function
  sets itself produced two entries of that name and
  [`do.call()`](https://rdrr.io/r/base/do.call.html) died with “formal
  argument ‘seed’ matched by multiple actual arguments”. Since
  [`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
  has both `k` and `seed` as formals, `rf_args = list(seed = 3)` —
  straight from the documented usage — was enough to trigger it. Extras
  now replace base entries by name. `data_sf`, `response_var`,
  `predictor_vars` and `folds` are protected and dropped with a warning,
  because a per-model override of those would silently make the models
  incomparable.

- [`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md)
  given a single bare `spatial_fit` reported all-`NA` Moran’s I columns
  and logged “‘fit’ is not a spatial_fit object” once per component. A
  `spatial_fit` is itself a list, so it passed the
  [`is.list()`](https://rdrr.io/r/base/list.html) check and the loops
  then iterated the fit’s own components as though they were models. It
  is now wrapped into a one-element named list, exactly as
  [`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md)
  already did.

- [`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)
  failed on its own documented fast path. With `FNN` and `Matrix`
  installed the weights are a sparse `Matrix`, and
  [`base::crossprod()`](https://rdrr.io/r/base/crossprod.html) does not
  S4-dispatch on the `dgeMatrix` that `W %*% resid` produces, so the
  call died with “requires numeric/complex matrix/vector arguments” —
  taking
  [`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md),
  which calls it automatically, down with it. Rewritten as
  `sum(resid_c * (W %*% resid_c))`, which is numerically identical and
  uses only dispatching primitives.
  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  carried the same bug.

- [`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)
  no longer errors on constant non-zero residuals: the degeneracy guard
  tested the raw sum of squares where Moran’s I is a function of the
  *centred* residuals, so `VI` came out `NaN` and `if (VI > 0)` raised
  “missing value where TRUE/FALSE needed”. A non-finite `VI` is handled
  too.

- `.build_knn_weights()`’s n \> 5,000 guard tests for both `FNN` **and**
  `Matrix`. Keyed on `FNN` alone, an unbounded dense n × n allocation
  went through whenever `FNN` was present but `Matrix` was not.

- [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
  drops columns of `features_sf` that would collide with the polygon ID
  column, with a logged warning.
  [`sf::st_join()`](https://r-spatial.github.io/sf/reference/st_join.html)
  suffixed them (`poly_id.x` / `poly_id.y`), which defeated the rename
  afterwards and left the result with **no rows** — reachable simply by
  re-assigning already-assigned points. A join that still fails to
  produce the ID column now errors and names the columns it did produce.

- [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  keeps the `"deff_applied"` attribute when `cells_sf` is supplied;
  [`dplyr::left_join()`](https://dplyr.tidyverse.org/reference/mutate-joins.html)
  rebuilds attributes from its `x` template and dropped it. The per-cell
  vector is remapped onto the joined row order, `NA` for cells holding
  no observations.

- [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  joins on the native ID type when both sides agree. Coercing
  unconditionally made the returned ID type depend on an unrelated
  argument and turned integer IDs into `"1"`, `"10"`, `"2"`, …. A
  genuine class mismatch still coerces both to character and logs why.

- [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  coerces non-POINT geometry before computing a variogram-based design
  effect.
  [`sf::st_coordinates()`](https://r-spatial.github.io/sf/reference/st_coordinates.html)
  returns one row per *vertex*, so a POLYGON or multi-vertex MULTIPOINT
  feature misaligned the coordinate matrix with the data and fed the
  wrong points into every cell.

- `make_folds(method = "buffered_loo")` errors when the buffer excludes
  so much of the data that no fold retains two training points. Those
  folds used to sail through and be dropped one at a time inside the CV
  loop, so the only symptom was a generic “all folds failed” warning at
  the very end.

- `make_folds(method = "block_kfold")` refuses a block size yielding a
  single block covering the whole extent — one fold with an empty
  training set, reported as a run that merely happened to score `NA`. An
  accepted autocorrelation range could trigger it:
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  rejects ranges above half the bounding-box *diagonal* while block
  construction needs half the *width*.

- [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  coerces MULTIPOINT geometry rather than merely accepting it, for the
  [`st_coordinates()`](https://r-spatial.github.io/sf/reference/st_coordinates.html)
  reason above; every fold was misaligned silently.

- Cross-validation no longer renumbers folds. `.remap_folds()` dropped
  unusable folds from a list, shifting every later fold’s index, so
  `fold_metrics$fold` and `predictions$fold` stopped lining up with
  `make_folds()$assignment$fold`. The original index is carried through.
  Folds left with fewer than two training rows are detected there and
  logged, instead of failing one at a time deeper in.

- [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
  rejects a `fit_fn` whose
  [`predict()`](https://rdrr.io/r/stats/predict.html) returns the wrong
  number of values. Both the metric computation and the prediction frame
  recycled silently, so two predictions against four test rows yielded a
  four-row frame with metrics computed against fabricated pairs.

- Cross-validation under `parallel = TRUE` reports a fold that died in a
  worker.
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  returns a `try-error` rather than `NULL`, which the `NULL` filter
  kept, and the failure surfaced as “subscript out of bounds”.
  [`conditionMessage()`](https://rdrr.io/r/base/conditions.html) has no
  method for a `try-error`, so the diagnostic branch itself threw; the
  condition is now taken from the object’s attribute.

- Cross-validation under `parallel = TRUE` is reproducible from `seed`
  and gives results identical to `parallel = FALSE`. `.cv_run_folds()`
  called
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  without seeding the fork streams, so each worker seeded itself from
  the clock and process ID. One seed per fold is now drawn in the
  parent, making each fold’s stream a function of `(seed, fold index)`
  alone.

- [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  rejects a singular variogram fit.
  [`gstat::fit.variogram()`](https://r-spatial.github.io/gstat/reference/fit.variogram.html)
  signals failure by setting `attr(., "singular")` and returning
  normally, so testing only for a `try-error` made the spherical
  fallback unreachable *and* let a singular fit’s `range` flow out as
  the estimated autocorrelation range — which
  `make_folds(auto_range = TRUE)` then sizes spatial blocks from.

- `.extract_gwr_values()` requires **every** model-matrix column to
  match a column of GWmodel’s `SDF` before multiplying the local
  coefficients through. A partial match reconstructed a linear predictor
  missing one or more terms and returned it as the fitted value —
  plausible numbers that were simply wrong, feeding
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and every metric
  with no warning. A non-numeric coefficient column is refused rather
  than coerced.

- [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  separates the three degenerate response cases. Folded together, an
  all-dropped dataset was reported as “binary (0 unique values)” and a
  constant response as “binary (1 unique value)”, while a genuinely
  binary non-integer response (1.5 / 2.5) failed the integer-like gate
  and passed unremarked.

- [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  and
  [`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md)
  validate `bandwidth`. Unvalidated, `NA` gave “missing value where
  TRUE/FALSE needed”, a length-2 vector gave “the condition has length
  \> 1”, and with `adaptive = FALSE` a zero or negative distance reached
  GWmodel untouched.

- [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)’s
  local-collinearity spot-check no longer touches the RNG at all. It
  sampled its 30 locations from the global stream and fires only when
  `n > 30` with at least two numeric predictors, so the same script
  produced different fold assignments depending on how many predictors a
  model happened to carry;
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md)
  calls it once per fold. The 30 locations are now evenly spaced ranks
  of the observations ordered by x, then y – reproducible, independent
  of the row order, and drawing no random numbers.

- [`predict.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.gwr_fit.md)
  returns an all-`NA` vector when every row of `newdata` is dropped as
  incomplete, matching the other two backends, rather than surfacing a
  raw sf-to-`Spatial` coercion error.

- [`predict.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.bayesian_fit.md)
  transforms `newdata` to the training CRS *before* cleaning it, and
  derives the surviving rows from one sentinel column instead of a
  second, separately-maintained copy of the cleaning rules. It errors
  when a predictor standardised at fit time is absent from `newdata` or
  has arrived as character — silently skipping it handed brms an
  unscaled column against a model fitted on a scaled one. Its failure
  path returns a matrix when `draws = TRUE`, honouring the documented
  return shape.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
  `spatial_fit` errors when there are no finite residuals, instead of
  producing a uniformly grey map from `limits = c(Inf, -Inf)`. A
  *perfect* fit is handled too: all-zero residuals gave
  `limits = c(0, 0)`, a degenerate diverging scale whose breaks collapse
  onto one value.

- [`plot_tessellation_map()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_tessellation_map.md)
  logs a warning for a `fill_col` that is not present, instead of
  drawing an unfilled outline map with nothing to say anything had gone
  wrong — a mistyped `label_col` already warned. `xlim`/`ylim` are
  validated, and the `theme` default moved out of the formals so a
  Suggests package never appears in an exported function’s default
  arguments.

- [`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md)
  announces when it *stamps* a CRS rather than reprojecting.
  [`sf::st_set_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html)
  only relabels; the coordinates do not move.
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  already made that assumption loudly.

- [`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md)
  rejects an EMPTY LINESTRING rather than misaligning the result.
  [`st_line_sample()`](https://r-spatial.github.io/sf/reference/st_line_sample.html)
  yields no midpoint for one (and segfaults in sf 1.0.x), so the sampled
  midpoints stopped corresponding 1:1 with the rows they are scattered
  back into. A count check backstops any other divergence.

- [`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md)
  errors on an unnamed list. The loop is over `names(fits)`, so an
  unnamed list iterated zero times and returned `NULL` silently;
  [`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md)
  then died in `seq_len(nrow(...))` nowhere near the cause.

- [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  coerces MULTIPOINT geometry rather than admitting it, and errors on a
  factor or character predictor by name instead of dying inside
  [`colMeans()`](https://rdrr.io/r/base/colSums.html) with “‘x’ must be
  numeric”.

- [`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md)
  passes both `cellsize` and `n` to
  [`sf::st_make_grid()`](https://r-spatial.github.io/sf/reference/st_make_grid.html)
  when both are known;
  [`st_make_grid()`](https://r-spatial.github.io/sf/reference/st_make_grid.html)
  does not ignore `n` in the presence of `cellsize` for square grids,
  and omitting it made sf recompute `nx = ceiling(w / cellsize)`, which
  floating-point division pushes one past the intended count. `n` is
  parsed and validated once, up front, instead of being silently coerced
  to `NULL` in one branch and erroring in the other.

- The grid cache key no longer truncates `target_cells`.
  [`as.integer()`](https://rdrr.io/r/base/integer.html) made 25.2 and
  25.7 collide on one key, so the second call silently received the
  first one’s grid, and a `NULL` `target_cells` collapsed
  [`paste0()`](https://rdrr.io/r/base/paste.html) to `character(0)`,
  crashing the lookup.

- [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
  normalises a CRS-less `points_sf` to `NULL` rather than `NA_crs_`,
  which is a list and so was not treated as “no CRS supplied”
  downstream. Hex and square grids are built in the points’ CRS, so the
  grid and the points no longer end up in different CRSs and break the
  point-to-cell index.

- `build_tessellation(method = "triangles")` triangulates the **point
  set** when `geometry` is unavailable, via
  [`sf::st_triangulate()`](https://r-spatial.github.io/sf/reference/geos_unary.html)
  on the unioned points. The fallback previously triangulated the convex
  hull *polygon*, discarding every interior point. The result is still
  the Delaunay triangulation of the input; only the resolution of
  degenerate configurations can differ from qhull’s, and the logged
  warning now says so.

- [`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md)
  logs a warning naming the geometry types when it drops non-polygonal
  rows, which it silently did before.

- [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)
  and
  [`voronoi_seeds_random()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_random.md)
  validate their inputs
  ([`voronoi_seeds_random()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_random.md)
  also accepts the `sfc` its documentation always promised), and
  clamping `k` to the number of distinct positions is logged.
  `get_voronoi_seeds(method = "provided")` logs a warning when `n`
  disagrees with `nrow(seeds)`, which it ignores.

- [`gp_lengthscale_bounds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gp_lengthscale_bounds.md)
  validates `coords_xy` and `q_small`. A vector `coords_xy` failed
  inside `.safe_dist()` with “argument is of length zero” and an
  out-of-range `q_small` inside
  [`quantile()`](https://rdrr.io/r/stats/quantile.html), neither naming
  the argument.

- [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
  validates the response before handing it to Stan, where nothing points
  back at the column, and validates `gp_k`, `gp_c` and `control`. The
  inverse-gamma prior is written with `%.10g` rather than `%.6f`: a
  small scale rounded to the literal `"0.000000"` and Stan rejected
  `inv_gamma(a, 0)` from deep inside the model block. Tightly clustered
  coordinates get there. The half-normal fallback’s scale is guarded the
  same way.

- [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
  names, in a warning, any `gwr_args` entry it drops.
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md)
  has no `...`, so entries meant for
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  alone (e.g. `longlat`) were discarded silently and simply had no
  effect.

- `.compute_reg_metrics()` errors on a `y_train_mean` that is neither a
  scalar baseline nor one value per observation, instead of recycling it
  against the filtered response and silently distorting R².

- **[`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md)
  no longer truncates the grid when `cellsize` and `n` are both
  supplied. This changes results.**
  [`sf::st_make_grid()`](https://r-spatial.github.io/sf/reference/st_make_grid.html)
  does not ignore `n` when `cellsize` is given: for square grids it
  takes the cell dimensions from `cellsize` *and* the counts from
  `nx = n[1]`, `ny = n[2]`, anchored at the bounding-box corner.
  `cellsize = 25` with `n = 2` on a 100 × 100 boundary therefore
  produced 4 cells covering 2,500 of 10,000 square units and silently
  left three quarters of the study area with no cells at all — and
  because `clip = TRUE` had nothing outside the boundary to discard, the
  result looked like an ordinary, complete grid. `cellsize` now wins,
  `n` is dropped with a logged warning naming what it would have done,
  and the same call returns 16 cells covering the whole boundary. `n` is
  still forwarded when the *package* derived `cellsize` from it or from
  `target_cells`, which is what the original code was written for:
  omitting it there lets sf recompute `ceiling(w / cellsize)` and
  floating-point division pushes the count one past the intended value.

- **[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  no longer refuses a continuous response that happens to take two
  values. This changes results: fits that used to error now run.** The
  guard rejected any response with exactly two distinct finite values as
  “binary” and pointed at `GWmodel::ggwr.basic(family = "binomial")`.
  Two distinct values is not the same thing as binary: a measurement
  censored at a detection limit or saturated at a ceiling (0.0031 /
  12.7401) is perfectly continuous, Gaussian GWR on it is a well-defined
  least-squares problem, and the advice to switch to a binomial family
  is nonsense for such values. The hard stop is now gated on the
  response also being integer-like, which is what the surrounding code
  already used to separate coded categories from measurements. A
  two-valued non-integer response raises a
  [`warning()`](https://rdrr.io/r/base/warning.html) naming the two
  values and asking you to confirm it is genuinely continuous, then
  fits. This also mattered inside
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
  where the guard runs once per fold and a small training fold can
  legitimately hold only two distinct values.

- **[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  no longer reports a Moran’s I that is arithmetically fixed. This
  changes which cell counts it returns.** `.morans_i_for_k()` builds a
  `min(8, n_cells - 1)`-nearest-neighbour weight matrix, so at nine
  cells or fewer every cell neighbours every other one. The
  row-standardised matrix is then complete, `W %*% e = -e/(n - 1)` for
  *any* mean-zero residual vector, and Moran’s I collapses to exactly
  `-1/(n_cells - 1)` whatever the data are. That is not merely
  uninformative: `|I| = 1/(n_cells - 1)` falls monotonically in the
  number of cells, so `criterion = "morans_i"` ranked the largest
  evaluated candidate first every time, and `"combined"` carried the
  same tilt at half weight. Candidates below the floor now return
  `NA_real_` and are excluded from the model-aware ranking; when none
  clears it — the usual outcome at the default `max_levels = 12`, since
  the search evaluates a window around the elbow — the call falls back
  to the geometric ranking and logs a warning. Raise `max_levels` above
  roughly 10 for the model-aware criteria to contribute at all.
  `predictor_vars` also accepts logical columns now, read as 0/1,
  matching
  [`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)/[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)/[`predict()`](https://rdrr.io/r/stats/predict.html);
  factor and character predictors are still refused by name.

- `residual_morans_i(fit, k = 1)` no longer errors with “subscript out
  of bounds” on a machine without `FNN`. In the dense fallback the inner
  function returns a scalar at `k = 1`, so
  [`apply()`](https://rdrr.io/r/base/apply.html) simplified the
  neighbour table to a length-n vector and
  [`t()`](https://rdrr.io/r/base/t.html) made it a 1 × n matrix;
  indexing `nn_idx[i, ]` then failed for every `i > 1`. The result is
  now forced to `n × k`.

- [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  no longer dies on an empty or non-finite geometry.
  [`st_coordinates()`](https://r-spatial.github.io/sf/reference/st_coordinates.html)
  yields one all-`NA` row per EMPTY POINT rather than zero rows, so a
  row-count check let them through: `block_kfold`’s
  [`st_intersects()`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html)
  returned `integer(0)`, `..block_id` went `NA`, and the nearest-block
  rescue aborted with “replacement has length zero”. Unusable rows are
  now dropped with a warning naming the count, after `..row_id` is
  stamped so the survivors keep their original row identities, and for
  every method rather than just `block_kfold` — `random_kfold` would
  otherwise put an unplottable point in a fold, and `nndm` and
  `buffered_loo` both feed the coordinates to distance code. The rescue
  itself uses [`vapply()`](https://rdrr.io/r/base/lapply.html) rather
  than [`apply()`](https://rdrr.io/r/base/apply.html), so a point whose
  distances are all `NA` keeps its `NA` instead of collapsing the
  assignment. `points_sf` with no usable coordinates at all is an error
  naming that, not a downstream one.

- Every cross-validation wrapper names the cause when folds fail.
  `.cv_run_folds()` returns each fold’s error text rather than a bare
  `NULL`, and
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
  and
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
  append `First error: ...` to both the logged and the R-level “all N
  folds failed” message. Running
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
  without `brms` installed previously produced five `fold N fit failed`
  warnings and an all-`NA` `$overall` with `n_pred = 0` in which the
  word “brms” never appeared.

- **The package’s log lines no longer depend on the user’s global
  `logger` configuration.** `logger` seeds a new namespace from the
  global one, so the `"spatialkit"` namespace inherited whatever
  formatter the user had set before loading — and every logging helper
  hands `logger` an *already formatted* string. Under a user’s
  `formatter_sprintf`, every package message containing a literal `%` —
  the CRS distortion figures in
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md),
  the local collinearity percentage in
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  — hard-errored with “too few arguments”, and because the helper logs
  *before* it raises the R warning, the warning the manual promises died
  with it. Under the default `formatter_glue` a `{...}` inside a fold
  error was re-evaluated. The namespace’s formatter is now pinned to
  `formatter_paste`, so the message logged is the message written.

- `fit_bayesian_spatial_model(backend = "auto")` chooses **cmdstanr**
  only when a CmdStan build is actually available. It chose it whenever
  the **cmdstanr** *package* could be loaded — a thin interface that is
  often installed without the toolchain it drives — so on such a machine
  every fit died inside the sampler with “CmdStan path has not been set
  yet. See ?set_cmdstan_path”. The package’s own weekly `check-brms` job
  was one such machine: it installs **cmdstanr** to satisfy Suggests and
  never builds CmdStan, and every scheduled run since the Stan smoke
  tests landed failed there. “auto” now falls back to **rstan**, which
  **brms** always brings, and logs the choice; an explicit
  `backend = "cmdstanr"` with no usable build is an error that says to
  run
  [`cmdstanr::install_cmdstan()`](https://mc-stan.org/cmdstanr/reference/install_cmdstan.html).

- `cv_bayes(seed = )` reaches the sampler.
  [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
  carries `seed = 123` and the per-fold `fit_args` never set it, so
  every fold of every run sampled from Stan seed 123 and changing `seed`
  changed nothing on fixed folds. Each fold now draws its own sampler
  seed from the fold’s seeded stream, as
  [`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
  does for the forest; a `seed` in `fit_args` still overrides it for
  every fold.

- `cv_*(seed = NULL)` is reproducible from
  [`set.seed()`](https://rdrr.io/r/base/Random.html) under
  `parallel > 1`, as the README promised without qualification. With
  `seed = NULL` no per-fold seeds were drawn, so each forked worker was
  seeded by `mclapply()` from the clock and the process ID — three runs
  after the same `set.seed(777)` gave 0.5687, 0.5649 and 0.5623 while
  the sequential call was reproducible. The per-fold seeds are now drawn
  from the caller’s current stream (advancing it, as any RNG-consuming
  call would), so the sequential and parallel paths are the same
  function of the state
  [`set.seed()`](https://rdrr.io/r/base/Random.html) left.

- Warnings raised inside a fold reach the caller from the parallel path.
  R conditions do not cross a fork, so under `parallel > 1` every
  warning raised by the model — including
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)’s
  documented integer-response warning, raised in every fold — reached
  nobody, while the numbers came back identical and the run looked like
  a clean version of the same analysis. The worker now collects them and
  the parent re-raises each distinct message once.

- [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
  (and therefore
  [`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md))
  returns the same typed, zero-row `fold_metrics` frame as
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md)
  and
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
  when every fold failed, so `subset(fold_metrics, RMSE < 5)` works
  instead of erroring on a missing column.

- [`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)
  refuses `k` large enough to make the neighbour matrix dense. The only
  size guard (n \> 5000) applied to the dense fallback; with **FNN** and
  **Matrix** present, `k >= n - 1` allocated `n (n − 1)` pairs
  unguarded. Requests above 2e7 pairs are now an error naming `k` and
  `n`.

### New features

- New
  [`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
  and
  [`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md):
  a `ranger` random forest as a first-class backend, returning an
  `rf_fit` that works with
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
  [`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md),
  [`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
  and [`plot()`](https://rdrr.io/r/graphics/plot.default.html) like any
  other model. Three defaults are opinionated: `include_coords = FALSE`
  (a forest given the coordinates memorises location and fails wherever
  it has not been — Meyer et al. 2019,
  <https://doi.org/10.1016/j.ecolmodel.2019.108815> — and random CV does
  not catch it);
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns
  **out-of-bag** predictions, so
  [`summary()`](https://rdrr.io/r/base/summary.html) on an `rf_fit` is
  not comparable with the other backends and says so
  (`$info$fitted_are_oob`); and importance defaults to permutation
  rather than impurity, which is biased toward continuous and
  high-cardinality predictors (Strobl et al. 2007,
  <https://doi.org/10.1186/1471-2105-8-25>). Compare backends with
  [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md),
  which now has an RF branch.

  [`predict()`](https://rdrr.io/r/stats/predict.html) on an `rf_fit`
  refuses the type confusions `ranger` would otherwise absorb silently:
  a numeric-at-fit predictor supplied as text (which ranger
  factor-codes, then applies numeric split thresholds to the codes), a
  logical-at-fit predictor supplied as text, and a categorical level the
  forest was never *grown* with — the level set is not enough, since a
  spatial fold holding out a whole class leaves a level with no training
  rows. Arguments that make ranger return a matrix
  (`predict.all = TRUE`, `type = "quantiles"`) are rejected rather than
  flattened column-major. A constant seed is supplied to ranger’s
  predict unless the caller passes one, so prediction does not consume
  the global RNG and `predict_surface(chunk_size = )` — a performance
  knob — cannot shift later random draws. `cv_rf(seed = )` reaches the
  forest in every fold, and gains `pointize`. Passing ranger’s own
  spelling of an argument the wrapper already sets (`num.trees`,
  `min.node.size`, `num.threads`, `mtry`, `importance`, `seed`) through
  `...` is an error naming the wrapper argument to use, rather than
  reaching `ranger()` twice. `num_threads` defaults to
  `getOption("mc.cores", 1L)` — one thread unless the session has opted
  in — for both the fit and
  [`predict()`](https://rdrr.io/r/stats/predict.html), rather than
  ranger’s own default of every core on the machine, and
  `cv_rf(parallel = )` runs each forked fold’s forest on one thread
  unless told otherwise, so the worker count is never multiplied by a
  thread count. See
  [`?fit_rf_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md).

- Every `cv_*()` and
  [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
  accept `folds` as a vector of fold labels, one per row —
  `make_folds()$assignment$fold`, the object most naturally to hand — in
  addition to a
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  result and a list of `train`/`test` splits, the three shapes
  [`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
  already took. The label vector used to fail with R’s “\$ operator is
  invalid for atomic vectors”.

- New
  [`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md),
  implementing the dissimilarity index of Meyer & Pebesma (2021,
  <https://doi.org/10.1111/2041-210X.13650>). Predictors are centred and
  scaled on the training data’s own statistics, optionally weighted by
  variable importance — by the importance itself, not its square root,
  matching `CAST`. A prediction point’s DI is its distance to the
  nearest training point in that space over the mean pairwise training
  distance, and the threshold is the outlier-removed maximum of the
  training data’s own DI. Pass the
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  result you actually validated with — the area is defined relative to a
  performance estimate, and a blocked estimate is a claim about
  predicting further away.

  A model fitted with `include_coords = TRUE` is measured in coordinate
  space, since an index that ignores location would report a point far
  outside the training extent as *inside* on ordinary covariate values
  alone; weights for the two coordinate columns default to the mean of
  those supplied, as the caller has never seen them. Non-`POINT`
  `newdata` is reduced to points, and a CRS present on one side is
  applied to the other. The zero-variance test is relative to each
  column’s magnitude rather than an absolute tolerance, so a predictor
  is not dropped for the **unit** it was recorded in. Categorical
  predictors are refused rather than dummy-coded; logicals are read as
  0/1. A
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  result is resolved by its `..row_id` **values**, which coincide with
  row positions only when the input carried no prior IDs. See
  [`?area_of_applicability`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md).

- New
  [`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md):
  greedy forward feature selection with **spatially blocked inner
  folds**, which is the whole point of having it. Random inner folds
  inside blocked outer folds select variables that look predictive only
  because nearby points leak between train and test, and the outer loop
  then reports honest-looking numbers for a dishonestly chosen feature
  set. `method` defaults to `"block_kfold"` and logs a warning if set to
  `"random_kfold"`. The empty set is scored first where the backend can
  fit it, so the first variable is judged against a null-model baseline
  rather than accepted unconditionally, and `history` carries that
  baseline as a `step = 0` row. Every candidate set is scored on the
  same observations — the completeness filter matches
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
  exactly, finiteness test included, so a candidate carrying a single
  `Inf` cannot be preferred for having an easier subset — and the inner
  folds are built once, before the sweep, rather than rebuilt per
  candidate. A `max_fits` budget guards against nesting a sweep inside
  leave-one-out outer folds. Where the backend cannot fit the empty set
  at all —
  [`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
  and
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  both refuse a zero-length `predictor_vars` — the probe is silent on
  the console: its per-fold failures go to the file trace only, rather
  than printing the same lines a genuinely failed run prints.

- New
  [`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md):
  wraps
  [`GWmodel::gwr.model.selection()`](https://rdrr.io/pkg/GWmodel/man/gwr.model.selection.html)
  (Lu et al. 2014, <https://doi.org/10.1080/10095020.2014.917453>) and
  returns a ranked table instead of two loosely-coupled lists. It is the
  fast, in-sample counterpart to
  [`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)
  — the same forward search scored by **AICc**, read from the documented
  `c(bandwidth, AIC, AICc, RSS)` layout of GWmodel’s `GWR.df`, which
  carries no column names; the result records whether the table arrived
  in that shape. Candidates must be numeric, and `dmat_max_n = Inf`
  means *always precompute* the distance matrix. Both limitations are
  documented rather than papered over: one bandwidth is shared by every
  candidate (which is what makes the criteria comparable), and the null
  model is never evaluated, so the result always names at least one
  predictor. When it disagrees with the blocked estimate, believe the
  blocked estimate. See
  [`?gwr_model_selection`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md).

- New
  [`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md):
  builds a regular grid over the training extent (or a grid you supply),
  joins covariates, predicts in chunks and returns `sf`. Supports
  `boundary` clipping, `cell_size` or approximate `n_cells`, and
  `se = TRUE` for a posterior-SD surface where the backend exposes
  draws.

- New [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
  for `spatial_fit`, with `type = "residuals"`, `"observed_predicted"`
  and `"variogram"` (the empirical residual variogram with the fitted
  model and effective range overlaid, so the fit can be judged rather
  than trusted). The variogram’s distance axis is labelled in the units
  of the CRS it was actually fitted in — metres of an auto-chosen zone
  for a lon/lat fit, not the caller’s degrees — it names the azimuth
  when a single direction is drawn, and a fit that did not converge says
  so in the caption. New
  [`plot_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_folds.md)
  maps a fold scheme, which is the fastest way to see whether spatial
  blocks separate the data or are smaller than the autocorrelation range
  and therefore leaking.

- [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  gains `method = "leave_location_out"`, which keeps every observation
  from a location (named by the new `group_var`) in the same fold.
  Repeated measurements at one site were previously unrepresentable:
  random k-fold splits them across folds, so the model is scored partly
  on sites it trained on.

- [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  gains `method = "nndm"`, implementing the distance-matching principle
  of Milà et al. (2022, <https://doi.org/10.1111/2041-210X.13851>), as
  in `CAST::nndm()`. Rather than choosing a `buffer` with nothing to
  justify it, the exclusion around each held-out point is sized so the
  training-to-test distance distribution reproduces the distances from
  your actual prediction locations (the new `prediction_points`) to the
  training data. The procedure follows the paper’s iterative exclusion
  removal for removal and is deterministic: no random numbers are drawn,
  so the caller’s RNG is untouched, and ties in the nearest-neighbour
  distance — every mutual-nearest-neighbour pair, all of a regular grid
  — are broken by the point’s position rather than by its row index, so
  identical data give identical folds whatever order the rows arrive in
  (`CAST` breaks them by row). `params$target_median`,
  `params$realised_median` and `params$max_ecdf_excess` record how close
  the match came, and `min_train` (default 0.5) and `phi` control it.
  Matching is as close as the training configuration permits — the
  achievable distances are discrete order statistics. When prediction
  locations sit no further from the training data than training points
  sit from each other, plain leave-one-out already reproduces the target
  and nothing is excluded; that is the correct outcome. A non-POINT
  `prediction_points` layer is reduced to points first, since
  point-to-polygon distances are zero for any cell containing a training
  point and would collapse the scheme towards plain LOO.

- [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  gains `deff = "variogram"`, computing a per-cell design effect from a
  fitted variogram rather than one pooled intra-class correlation. For
  `n` points in a cell with correlation matrix `R` the effective sample
  size of the mean is `n^2 / sum(R)`, so `deff = sum(R) / n`. This
  generalises the Kish option — a constant off-diagonal correlation
  recovers `1 + (n - 1) * rho` exactly — but lets correlation decay with
  distance, which is what having fitted a variogram is for. Pass the fit
  via the new `sac` argument, or it is estimated when `response_var` is
  supplied — on the **response**, not on OLS residuals, even when
  `predictor_vars` are listed: the `..se_resp_*` columns estimate the SE
  of the cell mean as an estimate of the response’s grand mean, so the
  correlation to correct for is the response’s own (measured grand-mean
  coverage 0.93 with the response variogram against 0.51 with the
  residual one on a field with a smooth predictor). Pass a residual
  variogram through `sac` if that is the field you want. Large cells are
  subsampled at `deff_max_n` (default 500), with the correlation scaled
  back to the cell’s own size. A `sac_range` whose fit was *rejected*
  carries no usable correlation function, so both the supplied and the
  internally estimated path fall back to `deff = 1` and say so rather
  than saturating the correlation at every within-cell distance. One
  correlation function is fitted and applied to every numeric column,
  response and predictors alike, because a variogram is a property of
  the field rather than of a variable type.

- [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
  supports intercept-only models (`predictor_vars = character(0)`): the
  response is explained by the intercept and the spatial GP alone, the
  natural null for asking how much of a surface is spatial structure
  rather than covariate effect.

- [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
  checks the posterior length-scale against the smallest scale the
  chosen basis can resolve and logs a warning when more than 10% of the
  posterior mass falls below it — the adequacy diagnostic recommended by
  Riutort-Mayol et al. (2023,
  <https://doi.org/10.1007/s11222-022-10167-2>), and what makes the
  smaller default `gp_k` safe rather than merely cheaper. `$info` gains
  `gp_c`, `gp_n_basis`, `gp_ell_min` and `gp_lengthscale_bounds`, and
  [`print()`](https://rdrr.io/r/base/print.html) on a `bayesian_fit` and
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)’s
  `fold_metrics` report the total basis count alongside the
  per-dimension rank.

- [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
  raises a condition when folds fail, matching
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md)
  and
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md);
  an all-failing `fit_fn` previously returned an all-`NA` `overall` and
  an empty `fold_metrics` with nothing at R condition level. The result
  records `n_folds_attempted` and `n_folds_succeeded` — compare them
  before trusting `overall`.

- [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  records the CRS the folds were built in as `params$crs`
  (`"EPSG:32632"`, an input string, or a WKT). `block_size` and
  `sac_range` are lengths in *that* CRS, which is not necessarily the
  one the caller passed: geographic input is projected by
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  to a CRS chosen for the extent. Without the label the units of a
  recorded block size were not recoverable from the result.

- [`spatialkit_quiet()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)
  is a new exported helper. Both
  [`logger::log_appender()`](https://daroczig.github.io/logger/reference/log_appender.html)
  and
  [`logger::log_threshold()`](https://daroczig.github.io/logger/reference/log_threshold.html)
  default to `index = 1`, which is the temp-file trace, so the two-line
  recipe in the README could not redirect or quieten the **console**
  echo (index 2) — there was no documented way to silence the package.
  The README now says so too.

### Documentation

- [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  documents its three return shapes (a range, a rejected range, and no
  fit at all) and which attributes each carries.

- [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  documents that `k` is not always honoured: `buffered_loo` and `nndm`
  always return `k = n`, and `block_kfold` and `leave_location_out`
  lower it when the geometry or the grouping cannot support the request.
  Read `folds$k`.

- [`new_spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md)
  documents the [`coef()`](https://rdrr.io/r/stats/coef.html) contract;
  [`summary.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.spatial_fit.md)
  and
  [`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
  document that their metrics are in-sample for a `gwr_fit` and a
  `bayesian_fit` but out-of-bag for an `rf_fit`;
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
  documents that the projected CRS is not an unconditional guarantee,
  since
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  passes a CRS-less dataset through unchanged when its coordinates do
  not look like lon/lat.

- The vignette and `inst/scripts/example_nc_demo.R` read fit quality
  from `fit$metrics$r_squared` and CV results from `cv$summary$rmse`.
  Neither field has ever existed. Because
  [`sprintf()`](https://rdrr.io/r/base/sprintf.html) returns
  `character(0)` when any argument has length zero, the reporting lines
  printed *nothing* rather than erroring, so the shipped vignette
  silently omitted every number it claimed to show. Both now use
  [`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
  and `$overall`.

- The demo’s Voronoi tessellation was built from all 300 observations
  rather than from the 40 k-means seeds it computed one line earlier —
  one cell per observation, a nearest-neighbour interpolation rather
  than an aggregation, compared side by side against two ~50-cell grids.
  The seeds are now used.

- The vignette builds as
  [`rmarkdown::html_vignette`](https://pkgs.rstudio.com/rmarkdown/reference/html_vignette.html)
  rather than `html_document`, guards its `ggplot2` and `geometry` use,
  demonstrates
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  instead of reimplementing it with
  [`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html)/[`summarise()`](https://dplyr.tidyverse.org/reference/summarise.html),
  and adds a spatial cross-validation section contrasting `block_kfold`
  against `random_kfold` on the same data.

- The package-level help page
  ([`?spatialkit`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-package.md))
  gains “The pipeline, in order” and “Where to start” sections, so
  [`help(package = "spatialkit")`](https://elkronos.github.io/gis_modeling_toolkit/reference)
  leads somewhere rather than presenting 40 exports in alphabetical
  order.

- Every exported function’s description now says *when to reach for it*
  rather than only what it does, and `@family` / `@seealso` links
  connect each step of the pipeline to the one before and after it —
  [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
  to
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md),
  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  to
  [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md),
  [`new_spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md)
  to
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
  and the two seeding functions to each other.
  [`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md)
  versus
  [`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md),
  and
  [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)
  versus
  [`voronoi_seeds_random()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_random.md),
  each say which to pick and why.

- [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
  documents that `boundary` is **required** for `method = "hex"` and
  `method = "square"` — the grid methods have no extent of their own —
  and optional for `"voronoi"` and `"triangles"`, which derive one from
  the points. The error existed; the requirement was not written down
  anywhere.

- [`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md)
  documents that `target_cells`, `cellsize` and `n` are three ways of
  sizing one grid and that exactly one should be supplied, that
  `cellsize` is in the units of the working CRS, and that `cellsize`
  takes precedence over `n`.

- [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  documents the nine-cell resolution floor on the model-aware criteria,
  why it exists, and that the whole call falls back to the geometric
  ranking when no candidate clears it.

- [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
  documents that dropping every requested backend is an error
  (`"no viable models."`) rather than an empty comparison, and that the
  returned frame carries only the models that actually ran, so callers
  should check which names are present rather than assuming one row per
  request.

- [`new_spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md)
  documents the two obligations on a custom backend: return an object
  built by the constructor, and define a `predict.<subclass>()` method —
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
  scores folds through the
  [`predict()`](https://rdrr.io/r/stats/predict.html) generic, so
  without one every fold fails.

- **README.** A new “Your own data” section shows both entry points —
  [`st_read()`](https://r-spatial.github.io/sf/reference/st_read.html)
  for a spatial file and
  [`read.csv()`](https://rdrr.io/r/utils/read.table.html) +
  [`st_as_sf()`](https://r-spatial.github.io/sf/reference/st_as_sf.html)
  for a table of coordinates — using the `nc.shp` demo shapefile shipped
  with `sf` so it runs anywhere. The README previously manufactured
  every example inline with a hard-coded `crs = 32632` and never showed
  data entering the package at all. A companion “CRS: what the numbers
  are in” subsection states that block sizes, buffers, bandwidths,
  variogram ranges and `expand` distances are in the units of the
  working CRS; that geographic input is projected automatically to a CRS
  chosen for the extent; and how to pin one.

- **README.** New guidance where none existed: how to choose among the
  four tessellation methods, how `k` and `block_size` trade off against
  the autocorrelation range, what to do when
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  returns `NA`, how to read a design effect, which model backend to
  reach for (with the recorded cost of each), and a “Troubleshooting”
  section covering the errors a new user actually hits first. A worked
  hex-grid example replaces the previous picture-only coverage of the
  grid methods.

- **README.** Three corrections. The
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  example showed a rejected range printing its attributes, which
  [`print.sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.sac_range.md)
  has not done since the attribute dump was removed; it now shows the
  bare `NA` and reads the attributes explicitly. The
  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  passage claimed the residual-autocorrelation criterion was doing work
  at cell counts where it is arithmetically degenerate. The test-suite
  paragraph said “exactly one” test guards on `brms`; six do, five of
  them additionally gated behind `SPATIALKIT_TEST_BRMS` so they never
  run in the matrix.

- **`inst/scripts/example_nc_demo.R`** said EPSG:2264 was projected “so
  distances are metric”. Its unit is the US survey foot, which is what
  the script’s own “Autocorrelation range: %.0f ft” line reports. The
  comment now says planar, and names the unit every distance, bandwidth
  and block size in the script is in.

- **Vignette.** `print(rf_fit)` and `summary(rf_fit)` report the same
  OOB RMSE but different R² (0.4733 against 0.4715). The vignette now
  explains why:
  [`print.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.rf_fit.md)
  echoes `ranger`’s `r.squared` (`1 - MSE/var(y)`, unbiased n − 1
  variance) while [`summary()`](https://rdrr.io/r/base/summary.html)
  recomputes `1 - SS_res/SS_tot` from the same out-of-bag predictions
  with an n denominator, so the unexplained fractions differ by exactly
  n/(n − 1).

- Every exported function has runnable examples: the eleven that shipped
  without any —
  [`clear_fitted_cache()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_fitted_cache.md),
  [`clear_grid_cache()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_grid_cache.md),
  [`clip_target_for()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clip_target_for.md),
  [`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md),
  [`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md),
  [`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md),
  [`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md),
  [`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md),
  [`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md),
  [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)
  and
  [`voronoi_seeds_random()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_random.md)
  — gained one, and the two `\dontrun{}` blocks say why they cannot be
  run (a Stan toolchain and minutes of MCMC).

- [`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)’s
  default null is described correctly: `null = "auto"` uses the Cliff &
  Ord *regression-residual* moments whenever the residuals are OLS
  residuals on the rebuilt design, and the randomisation null otherwise;
  the README said “the randomisation variance” without qualification.
  The type-I error of a random forest’s residual test is attributed to
  what the package actually feeds it — out-of-bag residuals, which are
  honest out-of-sample errors with their own spatial structure — not to
  “shrunk in-sample residuals”.
  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  gives the real reason for its nine-cell floor (the standardised
  deviate is 0/0 there, so the criterion carries no information) rather
  than an argument from `|I|` that its own details section had just
  called wrong.
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  notes that the “use the naive SE for the cell’s own mean” advice is
  calibrated under uniform within-cell sampling. `cv_gwr(bandwidth = )`
  states its units and semantics like its siblings. `quiet` is
  documented as “suppress this function’s progress messages” everywhere,
  with a pointer to
  [`spatialkit_quiet()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)
  for the console log echo it does not touch.

- README: the square-grid call returns 36 cells, not 32; the
  installation section no longer promises a specific version from CRAN;
  the resolution figure and the quick-start output are regenerated for
  the elbow-first ordering (`4 3 5`, `k = 4`).

- The DESCRIPTION now cites the methods it implements – Lu et al. (2014)
  for `GWmodel`, Riutort-Mayol et al. (2023) for the Hilbert space
  Gaussian process, Strobl et al. (2007) for the permutation importance,
  Mila et al.

  2022. for NNDM folds and Meyer and Pebesma (2021) for the area of
        applicability – each with its DOI, and quotes only software
        names. The Riutort-Mayol reference on
        [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)’s
        help page gives the article number (33, 17) rather than “33, 1”.

- No example is wrapped in `\donttest{}` any more: the fifteen that were
  – every fit, cross-validation, comparison, plotting and surface
  example that needs a Suggests package – run unconditionally behind
  their [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html)
  guards, the slowest in under 2 s. The two Stan examples
  ([`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md),
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md))
  keep `\dontrun{}` because they need a C++ toolchain and minutes of
  MCMC, and say so in a leading comment.

## spatialkit 1.0.0

CRAN release: 2026-08-07

First CRAN release, published 2026-08-07.

- CRS management:
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md),
  [`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md),
  [`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md),
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md).
- Voronoi, hexagonal, square and Delaunay tessellation
  ([`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
  and the `create_*_polygons()` functions), with boundary clipping,
  stable reproducible cell IDs
  ([`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md))
  and a memoised grid builder
  ([`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md)).
- Seeding
  ([`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
  [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md),
  [`voronoi_seeds_random()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_random.md))
  and resolution selection
  ([`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)).
- Feature-to-polygon assignment
  ([`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md))
  and cell-level aggregation with design-effect-corrected standard
  errors
  ([`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)).
- GWR
  ([`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md))
  and Bayesian spatial Gaussian process
  ([`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md))
  backends behind a common `spatial_fit` S3 class.
- Spatial cross-validation:
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  with `random_kfold`, `block_kfold` and `buffered_loo`;
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md);
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
  and
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md).
- Model comparison and diagnostics:
  [`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md),
  [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md),
  [`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md),
  [`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md).
- Tessellation mapping
  ([`plot_tessellation_map()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_tessellation_map.md))
  and scoped logging.
