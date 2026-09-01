# spatialkit 1.0.0.9000

This is the **development version**; nothing below is on CRAN yet.
Everything is relative to **1.0.0**, the version currently on CRAN
(published 2026-08-07).

Release will need a **major** version bump: three exported functions have been
removed, and several defaults change the result of a fit or a comparison. Both
are under "Breaking changes".

Throughout, *raises a warning* means a genuine R `warning()` — one
`tryCatch(warning = )` catches, `suppressWarnings()` suppresses and
`options(warn = 2)` escalates. *Logs a warning* means a `logger` message in the
`"spatialkit"` namespace, which none of those touch.

## Breaking changes

### Statistical corrections that change results (third audit pass)

Each item below was measured, and the measurement is quoted so you can judge
whether it affects an analysis you have already run.

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

* **`summarize_by_cell(deff = "variogram")` reported the design effect of the
  subsample, not of the cell.** A cell larger than `deff_max_n` was subsampled
  and then answered with `sum(R)/n_used` — the design effect of a cell of
  `deff_max_n` points. Measured on 4000 points with an exponential correlation
  of range 60: true deff **1821.8**, reported **228.6**. The mean off-diagonal
  correlation is now estimated from the subsample and scaled to the cell's own
  size, giving 1825.2. Un-subsampled cells are algebraically unchanged.

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

* **`gwr_model_selection()` ranked on AIC while labelling the result AICc.**
  GWmodel's `GWR.df` is documented and built as `c(bandwidth, AIC, AICc, RSS)`
  and is `rbind()`ed over unnamed vectors, so it never carries column names —
  making the positional read the path every real call takes, not a rare
  fallback. Column 2 is the uncorrected AIC. Executed against a faithful
  GWmodel stub on 200 points with two real predictors and one pure-noise
  variable: **column 2 selects the model containing the noise variable; column 3
  selects `a + b`.** The result now also reports whether the table had the
  documented four columns.

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

* **`predict.rf_fit()` no longer consumes the global RNG.** `ranger`'s predict
  draws a seed from the global stream when none is supplied, so the number of
  `predict()` calls a script happened to make — via `predict_surface()`'s
  `chunk_size`, a pure performance knob — shifted every later random draw,
  including fold assignments. A constant seed is now supplied unless the caller
  passes one; it does not affect a regression forest's predictions.

* **`fit_gwr_model()` rejects a non-numeric predictor.** `gwr.basic()` expands
  contrasts via `model.matrix()` and fits, but `gwr.predict()` does not and
  fails, so the model appeared to fit and then silently predicted all `NA`.

* **`fit_gwr_model()` no longer rejects a two-valued continuous response.** The
  "binary" error is now gated on the response being integer-like. A
  left-censored or saturated measurement (every observation at a detection limit
  or a ceiling) has two distinct values and is perfectly continuous; it now
  warns instead. The guard also runs once per fold inside `cv_gwr()`, where a
  small training fold can legitimately hold only two distinct values.

* **`gwr_model_selection()` works above `dmat_max_n`.** `bw.gwr()` branches on
  `missing(dMat)`, not `is.null(dMat)`, so passing `dMat = NULL` explicitly took
  the *supplied* branch and died — meaning **every dataset above 2000 points
  silently fell back to the arbitrary fallback bandwidth**, with the message
  blaming GWmodel. The argument is now omitted rather than passed as `NULL`. A
  fixed bandwidth with no distance matrix, which GWmodel cannot start, is now
  refused with a message naming both remedies.

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

* `n_models` in `gwr_model_selection()` is computed in double precision; the
  integer form overflowed to `NA` above 46341 candidates, so the guard meant to
  refuse an impossible sweep errored on exactly the inputs it existed to catch.

* `.morans_i_for_k()` returns `NA` at or below nine cells, where every cell
  neighbours every other and Moran's I collapses to exactly `-1/(k-1)` for any
  residual vector — a function of the cell count alone.

* `residual_morans_i(fit, k = 1)` works on machines without **FNN**. `apply()`
  simplified the length-1 result to a vector, making the neighbour index a
  1 x n matrix and every row after the first out of bounds.

### API and default changes

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
  1% of its mass below the estimated lower bound and 1% above the upper. Where
  the bounds are too wide or degenerate to calibrate against, the half-normal is
  used and the fallback is logged. The prior applied is recorded in
  `$info$gp_lscale_prior`.

* `ensure_projected()` no longer forces continental-extent data into a single
  UTM zone. Transverse Mercator scale error grows quadratically with distance
  from the central meridian, so data spanning the contiguous United States
  carried distance errors of roughly 7.5% near the extent edge, propagating
  silently into `estimate_sac_range()`, `make_folds(block_kfold)` block sizing,
  GWR bandwidth selection and the GP length-scale. Extents reaching more than 5°
  from the candidate zone's central meridian now receive an equal-area
  projection centred on the data — Albers for mid-latitudes, Lambert Azimuthal
  Equal Area otherwise — with a logged explanation. Only longitude offset
  triggers the switch, since `cos(lat)` shrinks the distance from the central
  meridian and a tall narrow north-south extent is UTM's design case. Bounding
  boxes spanning more than 180° of longitude cannot be distinguished from data
  straddling the antimeridian and fall back to EPSG:3857 with a logged warning.
  Pass
  `target_crs` to override.

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
  suffixed them (`poly_id.x` / `poly_id.y`), which defeated the rename afterwards
  and left the result with **no rows** — reachable simply by re-assigning
  already-assigned points. A join that still fails to produce the ID column now
  errors and names the columns it did produce.

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

* `make_folds(method = "nndm")` no longer destroys the caller's random number
  stream. It called `set.seed(seed)` unconditionally, and `set.seed(NULL)` does
  not mean "leave the RNG alone" — it re-initialises the generator from the clock
  and process ID. Since `seed` defaults to `NULL`, the ordinary call
  `make_folds(pts, k, "nndm", prediction_points = grid)` reseeded the session
  every time. It also no longer rebinds `cleanup`, the name the enclosing
  `.with_seed()` handler uses: `on.exit()` stores the *expression*, so both
  handlers resolved to the inner closure and the outer one — the only holder of
  the caller's pre-call seed — never ran.

* `make_folds(method = "nndm")` no longer errors when the median number of
  excluded training points is not a whole number; a `%d` format was applied to
  `median()`, which returns a double for an even-length vector.

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
  rows are detected there and logged, instead of failing one at a time deeper in.

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

* `area_of_applicability()` resolves a `make_folds()` result correctly. Every
  `make_folds()` branch emits `..row_id` **values** in its `train`/`test` slots,
  and those coincide with row positions only when the input carried no
  pre-existing `..row_id` — which `cv_gwr()`, `cv_bayes()` and `cv_spatial()` all
  stamp on before prepping. Indexing the predictor matrix with them computed the
  training reference distances for the wrong rows whenever the IDs landed inside
  `1:n`, silently shifting the threshold. IDs are now resolved with `match()`,
  and an ID absent from the training data is an error.

* `area_of_applicability()` measures a model fitted with `include_coords = TRUE`
  in coordinate space. `rf_fit` stores `predictor_vars` without `..x`/`..y`, so
  the index ignored location entirely and a point far outside the training
  extent with ordinary covariate values was reported *inside* the area of
  applicability. Weights for those two columns are filled in from the mean of
  the weights supplied, since the caller has never seen them.

* `area_of_applicability()` weights the scaled columns by the importance itself
  rather than its square root, matching CAST; drops a training *row* carrying
  `Inf` rather than the whole predictor; applies the same finiteness test to
  `newdata` as to the training data; rejects a fold with a row in both its train
  and test slots; does not treat unused factor levels as empty folds; and errors
  on an explicit `use_fnn = TRUE` when `FNN` is absent.

* `estimate_sac_range()` rejects a singular variogram fit.
  `gstat::fit.variogram()` signals failure by setting `attr(., "singular")` and
  returning normally, so testing only for a `try-error` made the spherical
  fallback unreachable *and* let a singular fit's `range` flow out as the
  estimated autocorrelation range — which `make_folds(auto_range = TRUE)` then
  sizes spatial blocks from.

* `estimate_sac_range()`'s rejected-range return is classed `sac_range`, so
  `print()` shows a bare `NA` instead of dumping the empirical variogram and the
  fitted `gstat` model to the console as raw attributes. It keeps the fitted
  variogram either way, so `plot(type = "variogram")` can draw exactly the case
  worth looking at: a curve that never reaches a sill.

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

* `fit_gwr_model()`'s local-collinearity spot-check no longer disturbs the
  caller's RNG. It sampled from the global stream and fires only when `n > 30`
  with at least two numeric predictors, so the same script produced different
  fold assignments depending on how many predictors a model happened to carry.
  `cv_gwr()` calls it once per fold.

* `select_features_forward()` builds the inner folds once, before the sweep,
  instead of once per candidate. Fold construction does not depend on the
  candidate set, so the old code produced the same splits `p^2/2` times, repeated
  every block-size warning as many times, and put `make_folds()` outside the
  `try()`, so one fold-construction failure killed the whole sweep.

* `select_features_forward()` scores every candidate set on the same rows. Its
  completeness filter now matches `prep_model_data()` exactly, including the
  finiteness test: a candidate carrying a single `Inf` still reached
  `cv_spatial()` with a smaller row set than its rivals — precisely the
  comparison the filter exists to prevent.

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

* `predict.rf_fit()` errors on a categorical predictor level the forest was
  never grown with; ranger 0.16 returns a plausible-looking number instead.
  Levels are coerced to the training levels so the codes ranger sees are the
  codes it was built with. Arguments that make ranger return a matrix
  (`predict.all = TRUE`, `type = "quantiles"`) are rejected, because
  `as.numeric()` would flatten them column-major into a vector of the wrong
  length.

* `fit_rf_model()` collapses an unset ranger diagnostic to `NA_real_`.
  `as.numeric(NULL)` is `numeric(0)`, which `%||%` does not catch, so
  `is.finite()` raised "argument is of length zero" and `sprintf()` printed
  nothing at all — reachable by forwarding `oob.error = FALSE` through `...`.

* `predict_surface()` handles a `cell_size` wider than the data extent, which
  aborted inside `seq()` with "wrong sign in 'by' argument"; the length-zero
  guard it was supposed to hit was unreachable. An empty `grid` and an empty
  `covariates` layer are rejected by name. With `se = TRUE` it draws the
  posterior once per chunk instead of twice.

* `predict_surface()` detects a partially-matched `summary` argument. R matches
  `summ = "median"` onto the `summary` formal of `predict()`, but the guard
  tested `list(...)$summary`, which does not partial-match in that direction —
  so `.pred` silently became `colMeans(draws)` for a caller who asked for
  medians.

* `plot()` on a `spatial_fit` errors when there are no finite residuals, instead
  of producing a uniformly grey map from `limits = c(Inf, -Inf)`. A *perfect*
  fit is handled too: all-zero residuals gave `limits = c(0, 0)`, a degenerate
  diverging scale whose breaks collapse onto one value.

* `plot_folds()` requires `method` and `k` on the folds object. Both are read
  into the title, and a list missing either collapsed `sprintf()` to
  `character(0)`, which ggplot2 renders as a plot with no title rather than
  erroring.

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

* `gwr_model_selection()` reports which column it read when the diagnostic table
  is labelled but not with `AICc`, rather than calling it unlabelled.

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

* `area_of_applicability()` accepts logical predictors, read as 0/1, instead of
  refusing them alongside factors and characters. `fit_rf_model()` fits a
  logical predictor, `cv_rf()` cross-validates it and `predict()` predicts with
  it, so refusing to compute the dissimilarity index for the same model was
  self-inconsistent. A logical's 0/1 standard deviation is meaningful, which is
  the property the scaling needs; factor and character predictors, which would
  need an arbitrary one-hot scaling, are still refused.

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

## New features

* New `fit_rf_model()` and `cv_rf()`: a `ranger` random forest as a first-class
  backend, returning an `rf_fit` that works with `cv_spatial()`,
  `predict_surface()`, `area_of_applicability()` and `plot()` like any other
  model. Three defaults are opinionated: `include_coords = FALSE` (a forest
  given the coordinates memorises location and fails wherever it has not been —
  Meyer et al. 2019, <doi:10.1016/j.ecolmodel.2019.108815> — and random CV does
  not catch it); `fitted()` returns **out-of-bag** predictions, so `summary()`
  on an `rf_fit` is not comparable with the other backends and says so
  (`$info$fitted_are_oob`); and importance defaults to permutation rather than
  impurity, which is biased toward continuous and high-cardinality predictors
  (Strobl et al. 2007, <doi:10.1186/1471-2105-8-25>). Compare backends with
  `compare_models_cv()`, which now has an RF branch. See `?fit_rf_model`.

* New `area_of_applicability()`, implementing the dissimilarity index of Meyer &
  Pebesma (2021, <doi:10.1111/2041-210X.13650>). Predictors are centred and
  scaled on the training data's own statistics, optionally weighted by variable
  importance; a prediction point's DI is its distance to the nearest training
  point in that space over the mean pairwise training distance, and the
  threshold is the outlier-removed maximum of the training data's own DI. Pass
  the `make_folds()` result you actually validated with — the area is defined
  relative to a performance estimate, and a blocked estimate is a claim about
  predicting further away. Categorical predictors are refused rather than
  dummy-coded. See `?area_of_applicability`.

* New `select_features_forward()`: greedy forward feature selection with
  **spatially blocked inner folds**, which is the whole point of having it.
  Random inner folds inside blocked outer folds select variables that look
  predictive only because nearby points leak between train and test, and the
  outer loop then reports honest-looking numbers for a dishonestly chosen
  feature set. `method` defaults to `"block_kfold"` and logs a warning if set
  to `"random_kfold"`. A `max_fits` budget guards against nesting a sweep inside
  leave-one-out outer folds.

* New `gwr_model_selection()`: wraps `GWmodel::gwr.model.selection()` (Lu et al.
  2014, <doi:10.1080/10095020.2014.917453>) and returns a ranked table instead
  of two loosely-coupled lists. It is the fast, in-sample counterpart to
  `select_features_forward()` — the same forward search scored by AICc. Both
  limitations are documented rather than papered over: one bandwidth is shared
  by every candidate (which is what makes the criteria comparable), and the null
  model is never evaluated, so the result always names at least one predictor.
  When it disagrees with the blocked estimate, believe the blocked estimate.
  See `?gwr_model_selection`.

* New `predict_surface()`: builds a regular grid over the training extent (or a
  grid you supply), joins covariates, predicts in chunks and returns `sf`.
  Supports `boundary` clipping, `cell_size` or approximate `n_cells`, and
  `se = TRUE` for a posterior-SD surface where the backend exposes draws.

* New `plot()` method for `spatial_fit`, with `type = "residuals"`,
  `"observed_predicted"` and `"variogram"` (the empirical residual variogram
  with the fitted model and effective range overlaid, so the fit can be judged
  rather than trusted). New `plot_folds()` maps a fold scheme, which is the
  fastest way to see whether spatial blocks separate the data or are smaller
  than the autocorrelation range and therefore leaking.

* `make_folds()` gains `method = "leave_location_out"`, which keeps every
  observation from a location (named by the new `group_var`) in the same fold.
  Repeated measurements at one site were previously unrepresentable: random
  k-fold splits them across folds, so the model is scored partly on sites it
  trained on.

* `make_folds()` gains `method = "nndm"`, implementing the distance-matching
  principle of Milà et al. (2022, <doi:10.1111/2041-210X.13851>). Rather than
  choosing a `buffer` with nothing to justify it, the exclusion around each
  held-out point is sized so the training-to-test distance distribution
  reproduces the distances from your actual prediction locations (the new
  `prediction_points`) to the training data. `params$target_median` and
  `params$realised_median` record how close the match came. Matching is as close
  as the training configuration permits — the achievable distances are discrete
  order statistics — and the procedure differs from the paper's iterative
  exclusion. When prediction locations sit no further from the training data
  than training points sit from each other, plain leave-one-out already
  reproduces the target and nothing is excluded; that is the correct outcome.

* `summarize_by_cell()` gains `deff = "variogram"`, computing a per-cell design
  effect from a fitted variogram rather than one pooled intra-class correlation.
  For `n` points in a cell with correlation matrix `R` the effective sample size
  of the mean is `n^2 / sum(R)`, so `deff = sum(R) / n`. This generalises the
  Kish option — a constant off-diagonal correlation recovers
  `1 + (n - 1) * rho` exactly — but lets correlation decay with distance, which
  is what having fitted a variogram is for. Pass the fit via the new `sac`
  argument, or it is estimated when `response_var` is supplied. Large cells are
  subsampled at `deff_max_n` (default 500).

* `fit_bayesian_spatial_model()` supports intercept-only models
  (`predictor_vars = character(0)`): the response is explained by the intercept
  and the spatial GP alone, the natural null for asking how much of a surface is
  spatial structure rather than covariate effect.

* `fit_bayesian_spatial_model()` checks the posterior length-scale against the
  smallest scale the chosen basis can resolve and logs a warning when more than
  10% of the posterior mass falls below it — the adequacy diagnostic recommended by
  Riutort-Mayol et al. (2023, <doi:10.1007/s11222-022-10167-2>), and what makes
  the smaller default `gp_k` safe rather than merely cheaper. `$info` gains
  `gp_c`, `gp_n_basis`, `gp_ell_min` and `gp_lengthscale_bounds`, and `print()`
  on a `bayesian_fit` and `cv_bayes()`'s `fold_metrics` report the total basis
  count alongside the per-dimension rank.

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

## Documentation

* Every runnable example is now runnable. Four exported functions
  (`select_features_forward()`, `predict_surface()`, `plot.spatial_fit()`,
  `plot_folds()`) carried `\dontrun{}` stubs referring to an undefined `pts`;
  they are now self-contained `\donttest{}` examples guarded on the backend they
  need. The only remaining `\dontrun{}` blocks are the two that run full MCMC
  via brms (`fit_bayesian_spatial_model()`, `cv_bayes()`).
  `compare_models_cv()` gained a runnable RF example.

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
