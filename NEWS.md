# spatialkit 1.0.0.9000

The heading has to carry a parseable version: R's `NEWS.md` reader keys entries
off it, and a heading such as "(development version)" makes `R CMD check`
report "No news entries found".

This is the **development version**; nothing below has been released.
Everything is relative to **1.0.0**, the released version (2026-08-07).

When this is released it will need a **major** version bump rather than a minor
one: three exported functions have been removed, and the default results of
every Bayesian fit change. Both are documented under "Breaking changes".

## Breaking changes

* `fit_bayesian_spatial_model()` now derives the number of GP basis functions
  (`gp_k`) and the boundary factor (`gp_c`) from the ratio of the estimated
  length-scale to the domain size, rather than from the number of
  observations.

  `brms::gp()` builds a full tensor grid over its covariates, so a
  `gp(..x, ..y, k = gp_k)` term carries `gp_k^2` basis functions — `gp_k` is
  the count *per dimension*, not the total rank. The previous rule reduced to
  `max(15, floor(sqrt(n)))` for any `n` above 45, which made `gp_k^2`
  identically `n`: a model at n = 10,000 carried 10,000 basis functions and an
  n × n design matrix, at which point the approximation was no longer
  approximating anything.

  The derived value is typically 21–25 per dimension and is largely
  independent of `n`. At n = 2,000 this cuts the basis count from 1,936 to
  576 and a two-chain fit from 3,445 s to 795 s — a 4.3x speedup,
  superlinear in the 3.4x basis reduction because a larger basis costs more
  per gradient evaluation *and* worsens the posterior geometry. Note that at
  small `n` the derived `gp_k` is slightly *larger* than before (22 versus 15
  at n = 200) — the change is a correction, not an optimisation, and results
  will differ in both directions. Pass `gp_k` and `gp_c` explicitly to restore
  the previous behaviour.

  `gp_c` was previously hard-coded at 1.5, which is too small whenever the
  length-scale exceeds roughly half the domain half-range; the derived value
  is typically 3.2–3.6. Too small a boundary truncates the domain and degrades
  the approximation for exactly the smooth, long-range surfaces spatial models
  are usually fitted to.

* Removed the legacy wrappers `evaluate_models()`, `evaluate_models_cv()`, and
  `phi_prior_bounds()`. Use `compare_models()`, `compare_models_cv()`, and
  `gp_lengthscale_bounds()` respectively.

## Bug fixes

* `make_folds(method = "nndm")` no longer overwrites the caller's random number
  stream. The branch rebound `cleanup`, the same name the enclosing
  `.with_seed()` handler used; `on.exit()` stores the *expression*, so both
  registered handlers resolved to the inner closure and the outer one — the only
  holder of the caller's pre-call `.Random.seed` — never ran.

* `make_folds(method = "nndm")` no longer errors when the median number of
  excluded training points is not a whole number. A `%d` format was applied to
  `median()`, which returns a double for an even-length vector, so the call died
  before returning.

* `make_folds()` now coerces MULTIPOINT geometry rather than merely accepting
  it. `sf::st_coordinates()` returns one row per *vertex*, so any multi-vertex
  feature made row `i` of the coordinate matrix a different feature than row `i`
  of the data, misaligning every fold silently.

* `make_folds(method = "block_kfold")` now refuses a block size that yields a
  single block covering the whole extent. That produced one fold with an empty
  training set — blocked cross-validation degenerating into nothing, reported as
  a run that merely happened to score `NA`. An accepted autocorrelation range
  could trigger it, because `estimate_sac_range()` rejects ranges above half the
  bounding-box *diagonal* while block construction needs half the *width*.

* `cv_spatial()` now raises a condition when folds fail, matching `cv_gwr()` and
  `cv_bayes()`. Previously an all-failing `fit_fn` returned an all-`NA` `overall`
  and an empty `fold_metrics` with nothing at R condition level. The result also
  records `n_folds_attempted` and `n_folds_succeeded`.

* Cross-validation under `parallel = TRUE` now reports a fold that died in a
  worker. `parallel::mclapply()` returns a `try-error` object rather than `NULL`,
  which the `NULL` filter kept, and the failure surfaced as
  "subscript out of bounds" with the real diagnosis lost.

* `estimate_sac_range()` keeps the fitted variogram on the `NA` it returns when
  the range exceeds the largest lag fitted. The value stays `NA` — the range is
  genuinely unidentified and must not size blocks — but the diagnostic was being
  computed and discarded, which left `plot(type = "variogram")` unable to draw
  exactly the case worth looking at: a curve that never reaches a sill. That
  plot now draws it and labels why no range is marked.

* `plot()` on a `spatial_fit` errors when there are no finite residuals instead
  of producing a uniformly grey map. All-`NA` residuals gave
  `limits = c(Inf, -Inf)`, which ggplot renders without complaint by mapping
  every point to `na.value`.

* `area_of_applicability()` now measures a model fitted with
  `include_coords = TRUE` in coordinate space. `rf_fit` stores `predictor_vars`
  without `..x`/`..y`, so the dissimilarity index ignored location entirely and
  a point far outside the training extent with ordinary covariate values was
  reported *inside* the area of applicability — the exact extrapolation the
  index exists to catch.

* `area_of_applicability()` weights the scaled columns by the importance itself
  rather than its square root, matching CAST, the reference implementation.

* `area_of_applicability()` drops a training *row* carrying `Inf` rather than
  the whole predictor. `complete.cases()` passes `Inf` through, and it then
  poisoned that predictor's mean and standard deviation, so the predictor was
  dropped with a warning falsely reporting it as having no variance.

* `area_of_applicability()` rejects a fold with a row in both its train and test
  slots (each was its own nearest neighbour at distance zero, collapsing the
  threshold), and no longer treats unused factor levels as empty folds.

* `predict_surface()` handles a `cell_size` wider than the data extent, which
  previously aborted inside `seq()` with "wrong sign in 'by' argument"; the
  length-zero guard it was supposed to hit was unreachable. With `se = TRUE` it
  now draws the posterior once per chunk instead of twice, where the caller has
  left the summary at its default.

* `select_features_forward()` scores every candidate set on the same
  observations. `cv_spatial()` re-runs `prep_model_data()` per candidate, which
  drops rows missing in *that* candidate's columns, so a step could compare an
  RMSE over 150 rows against one over 120 and prefer a variable for having an
  easier surviving subset. It also no longer returns `Inf` as a score when
  nothing was selected.

* `gwr_model_selection()` reports which column it read when the diagnostic table
  is labelled but not with `AICc`, rather than calling it unlabelled.


* Cross-validation with `parallel = TRUE` is now reproducible from `seed`, and
  produces results identical to `parallel = FALSE`. Previously `.cv_run_folds()`
  called `parallel::mclapply()` without seeding the fork streams, so each
  worker seeded itself from the current time and process ID. The package's own
  backends were unaffected (`brms::brm()` receives a fixed seed, GWR bandwidth
  selection is deterministic), but `cv_spatial()` accepts an arbitrary
  user-supplied `fit_fn`, so any stochastic learner gave irreproducible
  results. One seed per fold is now drawn in the parent process, making each
  fold's RNG stream a function of `(seed, fold index)` alone rather than of the
  execution path.

* `ensure_projected()` no longer forces continental-extent data into a single
  UTM zone. Transverse Mercator scale error grows quadratically with distance
  from the central meridian, so a dataset spanning the contiguous United
  States projected into one zone carried distance errors of roughly 7.5% near
  the extent edge. Those errors were silent and propagated into
  `estimate_sac_range()`, `make_folds(block_kfold)` block sizing, GWR
  bandwidth selection, and the GP length-scale. Extents reaching more than 5°
  from the central meridian of the candidate zone now receive an equal-area
  projection centred on the data — Albers for mid-latitudes, Lambert Azimuthal
  Equal Area for equatorial and high-latitude extents — with a logged
  explanation. Pass `target_crs` to override.

  Only longitude offset triggers the switch. Transverse Mercator error scales
  with distance from the central meridian, and `cos(lat)` shrinks that
  distance, so a tall narrow north-south extent is UTM's design case and keeps
  its zone: a 4°-wide, 60°-tall strip peaks around +0.02%, whereas an
  equal-area projection would distort it by about −1.1% at 30° from centre.

  Bounding boxes spanning more than 180° of longitude cannot be distinguished
  from data straddling the antimeridian, where the centroid lands on the
  opposite side of the planet. These now fall back to EPSG:3857 with a warning
  rather than centring a projection 180° from the data.

* The GP term is now built with `scale = FALSE`. `brms::gp()` otherwise
  rescales its covariates so the maximum Euclidean distance between two points
  is 1 and reports `lscale` in that space, while this package standardises the
  coordinates itself and expresses the length-scale prior, `gp_c` and the basis
  adequacy threshold in those standardised units. The two normalisations
  differed by roughly the maximum pairwise distance (~4.9 for standardised 2D
  coordinates), leaving the automatic `normal(0, sd)` prior on `lscale` about
  five times too diffuse — a likely contributor to divergent transitions and
  rejected initial values. There is now exactly one coordinate scaling.

* The Gaussian process now fits one length-scale per coordinate axis
  (`gp_iso = FALSE`). Coordinates are standardised per axis before fitting, so
  a single shared length-scale made the kernel anisotropic in the original CRS
  by the ratio `sd(X)/sd(Y)` — a property of how the sampling locations happen
  to be laid out, not of the process being modelled. The model now estimates
  directional structure from the data instead of inheriting it from the
  standardisation. Pass `gp_iso = TRUE` for the previous behaviour. This does
  not change cost: `brms::gp()` builds a tensor grid either way, so the basis
  count is `gp_k^2` regardless.

* The automatic GP length-scale prior is now a calibrated inverse-gamma rather
  than `normal(0, sd)`. A half-normal on a positive parameter puts its mode at
  zero, so most of its mass sat at length-scales shorter than the basis can
  resolve — precisely where the Hilbert-space approximation develops a funnel
  and the sampler divergences. The replacement pins 1% of its mass below the
  estimated lower bound and 1% above the upper, concentrating the prior on
  length-scales the data can actually identify. Where the bounds are too wide
  or degenerate to calibrate against, the previous half-normal is used and the
  fallback is logged. The prior actually applied is recorded in
  `$info$gp_lscale_prior`.

## New features

* New `select_features_forward()`: greedy forward feature selection with
  **spatially blocked inner folds**, which is the whole point of having it.
  Random inner folds inside blocked outer folds select variables that look
  predictive only because nearby points leak between train and test, and the
  outer loop then reports honest-looking numbers for a dishonestly chosen
  feature set — worse than not selecting at all, because the problem is now
  hidden behind a defensible-looking validation. `method` therefore defaults to
  `"block_kfold"` and warns if set to `"random_kfold"`.

  No new plumbing was required: `.cv_fit_one_fold()` calls `fit_fn(train_sf)`
  on the training slice only, so anything done inside `fit_fn` is already
  nested and leak-free. A `max_fits` budget guards against the cost of nesting
  a sweep inside leave-one-out outer folds, which multiplies the roughly
  `p^2 / 2 * k` fits by `n`.

* New `fit_rf_model()` and `cv_rf()`: a `ranger` random forest as a first-class
  backend, returning an `rf_fit` that works with `cv_spatial()`,
  `predict_surface()`, `area_of_applicability()` and `plot()` like any other
  model. Three choices are opinionated and documented rather than left to the
  caller to get wrong:

  `include_coords` defaults to `FALSE`. Handing a forest the x and y
  coordinates lets it reproduce the training surface almost exactly by
  memorising location, then fail badly anywhere it has not seen. Random
  cross-validation does not catch this — neighbouring points leak between
  folds, so the memorised surface scores well — which is how the practice
  became common. Meyer et al. (2019) show the collapse directly. Setting it to
  `TRUE` warns.

  `fitted()` returns **out-of-bag** predictions, not in-sample ones, following
  the convention of the random forest packages themselves. In-sample
  predictions from a forest are close to memorisation and would make
  `summary()` report a fictitious R-squared. The consequence is that
  `summary()` means something different for an `rf_fit` than for a `gwr_fit` or
  `bayesian_fit`, whose fitted values are genuinely in-sample; do not compare
  the two directly, use `compare_models_cv()`. This is recorded as
  `$info$fitted_are_oob`.

  The out-of-bag error itself is reported but labelled everywhere as what it
  is: a *random* hold-out, and therefore optimistic under spatial
  autocorrelation for exactly the reason random k-fold is — the trees that did
  not sample a point almost certainly sampled its neighbours. `cv_rf()` gives
  the blocked estimate.

  Importance defaults to permutation rather than impurity, which is biased
  toward continuous and high-cardinality predictors (Strobl et al. 2007). It is
  exposed as `$info$importance` and feeds
  `area_of_applicability(weights = pmax(fit$info$importance, 0))` — the `pmax`
  because permutation importance goes slightly negative for predictors that do
  not help, and the weight validator now says so when it rejects them.

  `coef()` on an `rf_fit` errors rather than returning something invented; a
  forest has no coefficients.

* New `area_of_applicability()`, implementing the dissimilarity index and area
  of applicability of Meyer & Pebesma (2021). A fitted model returns a number
  for any location you hand it, including locations whose predictor values look
  nothing like anything it was trained on. Those predictions are extrapolations
  dressed as interpolations, and a cross-validation score says nothing about
  them — the held-out folds were drawn from the same predictor distribution as
  the training data. The AOA marks where the score applies.

  Predictors are centred and scaled using the training data's own statistics
  (never the prediction data's, which would re-centre a far-away block and make
  it look familiar), optionally weighted by variable importance. The DI of a
  prediction point is its distance to the nearest training point in that space
  divided by the mean pairwise training distance; the threshold is the
  outlier-removed maximum of the training data's own DI. The index is invariant
  to the overall scale of `weights`, so importance values can be passed as-is.
  The weighting convention is the one CAST uses -- the scaled column is
  multiplied by the weight, so a weight `w` contributes `w^2` to the squared
  distance -- so numbers cross-check against the reference implementation.

  Pass the `make_folds()` result you actually validated with. Without it the
  training reference is each point's nearest neighbour anywhere in the data,
  which for clustered data is very close, giving a conservative AOA; with it
  the reference distances are larger and the AOA is correspondingly wider. That
  is not a loophole — the AOA is defined relative to a performance estimate,
  and a spatially blocked estimate is a claim about predicting further away.
  Buffered and NNDM folds contribute the training set they actually left
  available rather than "everything outside the fold", so the exclusion they
  were built to enforce is not silently undone.

  Categorical predictors are refused rather than dummy-coded, since a one-hot
  column has no meaningful standard deviation to scale by. Zero-variance
  training predictors are dropped with a warning that a prediction point taking
  a different value there is extrapolation the index cannot represent. Works
  with any model via `train_sf` and `predictor_vars`, not just a `spatial_fit`.

* New `gwr_model_selection()`: wraps `GWmodel::gwr.model.selection()` and
  returns a ranked table instead of the two loosely-coupled lists GWmodel
  produces, with the model list and the diagnostic matrix already joined. It
  is the fast, in-sample counterpart to `select_features_forward()` — the same
  forward search, scored by AICc rather than by a spatially blocked
  cross-validated estimate. Both limitations of the method are documented
  rather than papered over: one bandwidth is shared by every candidate model
  (which is what makes the criteria comparable, and also means the bandwidth is
  not optimal for any model but the largest), and the null model is never
  evaluated, so the result always names at least one predictor and cannot tell
  you that none of the candidates help. When it disagrees with
  `select_features_forward()`, believe the blocked estimate.

  The bandwidth is clamped against the *largest* model in the sweep rather than
  the smallest, since a bandwidth that supports a one-predictor local fit can
  leave the full model underdetermined partway through the sweep. An `n × n`
  distance matrix is computed once and reused across every candidate fit and
  the bandwidth search, capped by `dmat_max_n` so the memory cost stays
  bounded, and `max_models` refuses a sweep that would fit more GWR models than
  asked for — the count grows as `p * (p + 1) / 2`.

  The criterion column is read by name when GWmodel labels it and positionally
  otherwise (which is what GWmodel's own documentation does); the result
  records which, so a silent change in GWmodel's return shape shows up in
  `$criterion` rather than in the ranking. In practice GWmodel does not label
  it, so the positional path is the live one.

  GWmodel's progress output is discarded by default (`quiet = TRUE`). It is
  written with bare `cat()`, so neither `suppressMessages()` nor
  `suppressWarnings()` touches it, and it scales with the square of the
  candidate count — 19 candidates would print 190 "Now calibrating the model"
  blocks plus the full bandwidth-search trace. Pass `quiet = FALSE` to watch a
  long sweep.

* `make_folds()` gains `method = "leave_location_out"`, which keeps every
  observation from a location (named by the new `group_var`) in the same fold.
  Repeated measurements at the same site were previously unrepresentable:
  random k-fold splits them across folds, so the model is scored partly on
  sites it trained on and reports excellent performance for having memorised
  them.

* `make_folds()` gains `method = "nndm"`, implementing the distance-matching
  principle of Milà et al. (2022). Rather than choosing a `buffer` with nothing
  to justify it, the exclusion around each held-out point is sized so the
  resulting training-to-test distance distribution reproduces the distances
  from your actual prediction locations (the new `prediction_points`, naturally
  the grid from `predict_surface()`) to the training data. Each held-out point
  draws a target distance from the empirical prediction-distance CDF, and the
  nearest training points are excluded up to the one whose distance is closest
  to it. Matching is as close as the training configuration permits rather than
  exact — the achievable distances are discrete order statistics — and the
  procedure differs from the paper's iterative exclusion.

  The returned `params` record `target_median` and `realised_median` (plus the
  full `target_distances` and `realised_distances`) so the matching can be
  checked rather than assumed. Note that when prediction locations sit no
  further from the training data than training points sit from each other,
  plain leave-one-out already reproduces the target and nothing is excluded —
  that is the correct outcome, not a failure. NNDM also cannot pull training
  points closer, so it cannot match a target shorter than the training
  nearest-neighbour distances.

* New `plot()` method for `spatial_fit` objects, with `type = "residuals"`
  (mapped at the training locations, where visible structure means unmodelled
  autocorrelation), `"observed_predicted"`, and `"variogram"` (the empirical
  residual variogram with the fitted model and effective range overlaid, so the
  fit can be judged rather than trusted). The package previously had `print()`
  and `summary()` methods but no `plot()`.

* New `plot_folds()` maps a cross-validation fold scheme, which is the fastest
  way to see whether spatial blocks actually separate the data or are smaller
  than the autocorrelation range and therefore leaking.

* New `predict_surface()`: builds a regular grid over the training extent (or
  a grid you supply), joins covariates, predicts in chunks and returns `sf`.
  `predict()` previously required `newdata` to be constructed by hand, which
  made producing a map — the most common thing anyone wants from a fitted
  spatial model — more work than it should be. Supports `boundary` clipping,
  `cell_size` or approximate `n_cells` resolution, and `se = TRUE` for a
  posterior-SD surface where the backend exposes draws. Chunking matters for
  `bayesian_fit`, where the draw matrix is `n_draws x n_newdata` and a fine
  grid would exhaust memory long before the fit itself would.

* `summarize_by_cell()` gains `deff = "variogram"`, computing a per-cell design
  effect from a fitted variogram rather than a single pooled intra-class
  correlation. For `n` points in a cell with correlation matrix `R`, the
  effective sample size of the mean is `n^2 / sum(R)`, so `deff = sum(R) / n`.
  This generalises the existing Kish option — substituting a constant
  off-diagonal correlation recovers `1 + (n - 1) * rho` exactly — but lets
  correlation decay with distance, which is what having fitted a variogram is
  for. Kish assumes every pair in a cell is equally correlated regardless of
  separation, and that assumption degrades as cells grow. Pass the fit via the
  new `sac` argument, or it is estimated when `response_var` is supplied.
  Large cells are subsampled at `deff_max_n` (default 500) before forming the
  correlation matrix.

* `fit_bayesian_spatial_model()` checks the posterior length-scale against the
  smallest scale the chosen basis can resolve and warns when more than 10% of
  the posterior mass falls below it. This is the adequacy diagnostic
  recommended by Riutort-Mayol et al. (2023) and is what makes the smaller
  default `gp_k` safe rather than merely cheaper. The threshold is stored as
  `$info$gp_ell_min` and the result as
  `$info$convergence_diagnostics$gp_lscale_below_resolution`.

* `print()` on a `bayesian_fit` and the `fold_metrics` frame from `cv_bayes()`
  now report the total GP basis count (`gp_n_basis`) alongside the
  per-dimension rank, so the real cost of a fit is visible.

* `$info` on a `bayesian_fit` now includes `gp_c`, `gp_n_basis`, and
  `gp_ell_min`.

* `fit_bayesian_spatial_model()` now validates `gp_k` and `gp_c` instead of
  letting an invalid value surface later as an opaque brms formula error.

## Infrastructure

* Added GitHub Actions workflows: `R-CMD-check` across macOS, Windows, and
  Linux on R release, devel, and oldrel-1, with a dedicated job that installs
  the optional backends so their code paths are actually exercised rather than
  skipped; and a weekly `check-brms` workflow for the Stan-dependent paths.

* Added regression tests for parallel/sequential CV equivalence, GP basis
  selection, and CRS selection at wide extents — the three areas the previous
  suite did not cover. The parallel tests use a deliberately stochastic
  learner, since a deterministic one passes whether or not fold seeding works.

## References

Lu, B., Harris, P., Charlton, M. and Brunsdon, C. (2014). The GWmodel R
package: further topics for exploring spatial heterogeneity using
geographically weighted models. *Geo-spatial Information Science* **17**(2),
85–101. <doi:10.1080/10095020.2014.917453>

Meyer, H. and Pebesma, E. (2021). Predicting into unknown space? Estimating the
area of applicability of spatial prediction models. *Methods in Ecology and
Evolution* **12**(9), 1620–1633. <doi:10.1111/2041-210X.13650>

Milà, C., Mateu, J., Pebesma, E. and Meyer, H. (2022). Nearest neighbour
distance matching leave-one-out cross-validation for map validation. *Methods
in Ecology and Evolution* **13**(6), 1304–1316.
<doi:10.1111/2041-210X.13851>

Riutort-Mayol, G., Bürkner, P.-C., Andersen, M. R., Solin, A. and Vehtari, A.
(2023). Practical Hilbert space approximate Bayesian Gaussian processes for
probabilistic programming. *Statistics and Computing* **33**, 1.
<doi:10.1007/s11222-022-10167-2>
