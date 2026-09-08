# cran-comments

**Before submitting**, replace every "<fill in>" below with the result of the
run it names: the local macOS check, the two win-builder runs, and the three
GitHub Actions workflows. Everything else in this file describes the tree as
it is.

## What this submission is

**This is an update to `spatialkit` 1.0.0, published on CRAN 2026-08-07.**

The version being submitted is **2.0.0**. A major bump is the honest choice:
three exported functions are removed, and several defaults change the result of
a fit or a comparison — the same script gets a different block size, a
different projection and different standard errors than it did under 1.0.0.
See "Breaking changes" in `NEWS.md`.

### Reverse dependencies

The CRAN page for 1.0.0 lists no reverse depends, imports or suggests, so there
is nothing to break. Confirm with `revdepcheck::revdep_check()` before
submitting rather than relying on this note. (This container has no network
access, so nothing here was checked against CRAN.)

## Summary

`spatialkit` builds tessellations over point data, aggregates observations to
cells, and cross-validates spatial models with folds that respect the data's
autocorrelation structure. Three model backends (GWR via `GWmodel`, a Bayesian
Gaussian process via `brms`, random forests via `ranger`) sit behind one S3
class so they can be scored on identical folds.

## Test environments

Ubuntu 24.04.4 LTS, R 4.3.3, x86_64-pc-linux-gnu — the environment the results
below were produced on. Installed: sf 1.0.15, dplyr 1.1.4, logger 0.2.2, digest
0.6.34, testthat 3.2.1, sp 2.1.2, ranger 0.16.0, tibble 3.2.1, geometry 0.4.7,
gstat 2.1.1, ggplot2 3.4.4, patchwork 1.2.0, FNN 1.1.4, Matrix 1.6.5,
roxygen2 7.3.1, knitr 1.45, rmarkdown 2.25, pkgload 1.3.4, spdep 1.3.1,
brms 2.20.4.

**`GWmodel` was present only as a faithful API stub built from CRAN sources** --
argument names and order, return shapes, SDF column naming and the `GWR.df`
column order, with real weighted least squares -- because the real package would
not install in this container. `brms` is installed, and `rstan` was made to
compile by pointing it at the system Boost headers, which is how the
`backend = "auto"` fix was verified: the five Stan smoke tests in
`test-bayes-smoke.R` sampled and passed here (34 expectations). That is the only
sampling done in this container; every other brms finding was established with
`make_stancode()`, `make_standata()`, `get_prior()` and `validate_prior()`,
which need no compilation. **Every GWR and brms result recorded below must be
re-confirmed on a machine with the real packages before submission** -- the
full suite has passed on macOS with both installed (6,162 expectations before
the Low-list pass, below). Not installed: `cmdstanr`.

Still to run before any submission:

* local: macOS (aarch64-apple-darwin), R release — <fill in>
* win-builder: R-devel — <fill in>
* win-builder: R-release — <fill in>
* GitHub Actions `R-CMD-check.yaml`, five jobs (not a cross-product):
  macOS-latest R release; windows-latest R release; ubuntu-latest R devel;
  ubuntu-latest R release; ubuntu-latest R oldrel-1. R-devel and oldrel-1 are
  tested on Linux only. — <fill in>
* GitHub Actions `backends` job (ubuntu-latest, R release) with `sp`, `GWmodel`,
  `gstat`, `FNN`, `Matrix`, `geometry`, `ranger`, `tibble` and `spdep`
  installed, so the optional code paths actually execute rather than skip.
  — <fill in>
* GitHub Actions `check-brms.yaml` (weekly, ubuntu-latest, R release) with
  `brms` and the Stan toolchain, which is the only job that sets
  `SPATIALKIT_TEST_BRMS` and therefore the only one that runs the five Stan
  smoke tests. Its scheduled runs of 2026-08-31 and 2026-09-07 failed because
  `backend = "auto"` chose `cmdstanr` on the package's presence alone while
  the runner has no CmdStan build (fixed in this tree; see `NEWS.md`), so it
  must be re-run by hand (`workflow_dispatch`) on the submitted commit.
  — <fill in>

## R CMD check results

<fill in: `R CMD check --as-cran` on a release-R machine with LaTeX and
network access, so the manual and the network-dependent incoming checks (URL
validity, the `Additional_repositories` lookup for `cmdstanr`) are exercised for
real.>

On the environment above, `R CMD check --as-cran --run-donttest --timings
--no-manual` on the built tarball -- with the vignette built, the `GWmodel`
stub and `brms` installed, `_R_CHECK_FORCE_SUGGESTS_=false` because `cmdstanr`
is not installed, and the incoming-feasibility check pointed at a local copy
of CRAN's package index because this container cannot reach CRAN -- reports
**3 NOTEs and nothing else**:

* **NOTE -- "checking CRAN incoming feasibility".** Two parts. The "Possibly
  misspelled words in DESCRIPTION" list, discussed under "Comments for
  reviewers"; and six URLs reported as 403 or "CONNECT tunnel failed" -- the
  package's own GitHub and issues pages, the README badge and licence links,
  r-spatial.github.io and www.r-project.org -- because every outbound
  connection from this container is refused by its proxy. Expect the spelling
  part on CRAN and nothing else.

* **NOTE -- "Package suggested but not available for checking: 'cmdstanr'".**
  Environmental, and expected: `cmdstanr` is not on CRAN. It is reached through
  `Additional_repositories`, is used strictly conditionally via
  `requireNamespace()`, and the `check-brms` CI job installs it.

* **NOTE -- "checking for future file timestamps ... unable to verify current
  time".** The container cannot reach the clock service the check consults.
  It does not occur on a networked machine.

`checking examples` runs all 42 examples in 7.0 s elapsed on the latest run
(6.5--8.7 s across runs); the slowest, `select_features_forward()`, takes
1.6--2.3 s. `checking tests` takes 72--93 s and `checking re-building of
vignette outputs` 18--26 s.

Everything previously recorded here as a blocker is resolved:

* The **ERROR in "checking tests" is gone.** Those 15 failures were tests
  asserting the pre-fix behaviour of `.morans_i_for_k()`. They have been
  rewritten against the corrected contract, along with the tests that encoded
  the pre-fix GP domain measure, GWR criterion column, Kish standard error and
  subsampled design effect. `checking tests` passes.

* The **non-ASCII WARNING** is gone: `checking R files for non-ASCII characters`
  passes. The literal U+2014 EM DASH that sat inside a string literal in
  `ensure_projected()`'s `target_crs` error message has been replaced with
  `--`. Em dashes elsewhere in `R/` sit in comments and roxygen blocks, which
  the check tolerates.

* The **syntax-error WARNING** is gone: `checking R files for syntax errors`
  passes with no output. It previously carried only a failed
  `Sys.setlocale("LC_CTYPE", "en_US.UTF-8")` from the check's own locale switch.

Two artefacts of this container are worth recognising if they appear again.
`checking package dependencies` emits
`Warning: unable to access index for repository ...` lines, because there is no
outbound network access here; it is not a check condition and does not affect
the status line. And building with `--no-build-vignettes` adds two WARNINGs
("Files in the 'vignettes' directory but no files in 'inst/doc'" and
"Directory 'inst/doc' does not exist") that a normal `R CMD build` does not
produce.

`R CMD build` produces a 1.2 MB tarball, of which the built vignette HTML is the
bulk. `checking running R code from vignettes` passes.

`testthat` reports **6,225 passing, 0 failures, 0 errors, 0 warnings, 9 skips**
in this container with `NOT_CRAN=true` and both backends present (the `GWmodel`
stub, real `brms`), and no test runs with zero assertions. On macOS with the
real `GWmodel` the count before the Low-list pass was 6,162 with the same
nine skips; the macOS `<fill in>` above carries the current one. The nine
skips are:

* 5 Stan smoke tests in `test-bayes-smoke.R`, skipped because
  `SPATIALKIT_TEST_BRMS` is unset -- `skip_if_not_installed("brms")` alone was
  not enough, since these compile Stan models and would otherwise run in any
  matrix job that happened to have `brms`
* 1 Windows-only fallback path that cannot run on Linux
* 3 that skip *because* an optional backend is installed: they assert the
  behaviour seen when `GWmodel` or `brms` is absent

Under `R CMD check`, where `NOT_CRAN` is unset, nine further tests skip on
purpose -- the `parallel::mclapply()` fork tests in `test-cv-parallel.R` and
`test-audit-pass6.R`, and two slow simulation checks -- all `skip_on_cran()`.
The check's own count under `--as-cran` is **6,200 passing, 0 failures, 18
skips**; `--as-cran` sets `_R_CHECK_LIMIT_CORES_`, which switches off two
expectations in `test-core-count.R` that read the machine's core count, so a
plain `R CMD check` counts two more.

With the optional backends absent (`sp`, `GWmodel`, `ranger`, `brms`, `gstat`,
`FNN`, `geometry`, `loo`, `patchwork`, `spdep` hidden, as in the CI matrix
jobs), `R CMD check` reports the Suggests-not-available NOTE and nothing else:
all 42 examples run, since each guards its optional packages with
`requireNamespace()`, the vignette builds, and the suite reports **2,234
passing, 0 failures, 139 skips** -- every test that needs a backend skips rather
than fails.

`README.md` does not restate these counts, precisely so the two cannot drift
apart.

## Breaking changes since the 1.0.0 tag

Three exported functions have been removed: `evaluate_models()`,
`evaluate_models_cv()` and `phi_prior_bounds()`. All three were thin wrappers
the package's own documentation described as legacy, and each has a documented
replacement (`compare_models()`, `compare_models_cv()` and
`gp_lengthscale_bounds()`). `NAMESPACE` now exports 41 objects.

Default results from `fit_bayesian_spatial_model()` also change, as a
consequence of three corrections to the Gaussian process path (the basis count,
the coordinate scaling, and the length-scale prior). The previous behaviour
remains reachable by passing `gp_k`, `gp_c` and `gp_iso` explicitly.

Three further corrections change results without being deliberate default
changes, so they are recorded under "Bug fixes" rather than "Breaking changes",
but a user upgrading should know about them:

* `create_grid_polygons()` supplied with both `cellsize` and `n` used to return
  a grid truncated to `n[1]` x `n[2]` cells anchored at the bounding-box corner,
  silently covering only part of the boundary. `cellsize` now wins and `n` is
  dropped with a logged warning.
* `fit_gwr_model()` used to refuse any response with exactly two distinct finite
  values as "binary". A censored or saturated continuous measurement has two
  distinct values and is a well-defined Gaussian GWR problem; the hard stop is
  now gated on the response also being integer-like, and the non-integer case
  warns and fits. Fits that used to error now run.
* `determine_optimal_levels()` no longer reports a Moran's I that is fixed by
  arithmetic. Below a nine-cell floor the criterion returns `NA` and the call
  falls back to the geometric ranking, which changes which cell counts it
  returns at the default `max_levels`.

`NEWS.md` is the full record relative to 1.0.0: 57 corrections that change
results (across four audit passes), 22 API and default changes, 26 new guards
and message changes plus the 14 entries of the sixth pass's Low list (two of
which change a number: nested variogram models in the design effect, and
`cv_bayes()`'s `yhat_sd` column), 54 bug fixes, 15 new features and 21
documentation entries.

## What was wrong in 1.0.0

In descending order of user impact:

1. `residual_morans_i()` failed on its own documented fast path. With `FNN` and
   `Matrix` installed the weights are a sparse `Matrix`, and `base::crossprod()`
   does not S4-dispatch on the `dgeMatrix` that `W %*% resid` produces, so the
   call errored — taking `compare_models()`, which invokes it automatically,
   down with it. This is a live failure on the recommended configuration, not a
   corner case.

2. The automatic length-scale prior for the Gaussian process was expressed in
   the wrong units. `brms::gp()` defaults to `scale = TRUE`, which rescales its
   covariates so the maximum pairwise distance is 1 and reports `lscale` in that
   space. This package standardises the coordinates itself and derived the prior
   in *those* units, so the two normalisations differed by roughly the maximum
   pairwise distance and the prior was about five times too diffuse. In practice
   that produced rejected initial values ("Gradient evaluated at the initial
   value is not finite") on most chains. The GP term now sets `scale = FALSE`.

3. The number of GP basis functions was chosen from the number of observations
   rather than from the spatial structure of the data. `brms::gp()` builds a
   tensor grid, so `gp(x, y, k = k)` carries `k^2` basis functions, and the
   previous rule reduced to `max(15, floor(sqrt(n)))` for any n above 45 —
   making the basis count identically n. A model at n = 10,000 carried 10,000
   basis functions. The count is now derived from the length-scale-to-domain
   ratio following Riutort-Mayol et al. (2023, Statistics and Computing 33:17).
   Measured on `dev/baseline-structural.rds`: at n = 2,000, `gp_k` 44 to 24 and
   the basis count 1,936 to 576 on the elongated layout (44 to 22 and 1,936 to
   484 on the clustered one); at n = 10,000, 100 to 23 and 10,000 to 529; at
   n = 200 the basis is *larger* than before (225 to 529), which is the expected
   consequence of a correction rather than an optimisation.

4. `ensure_projected()` forced data of any extent into a single UTM zone.
   Transverse Mercator scale error grows quadratically with distance from the
   central meridian, so continental-extent data carried distance errors of
   several percent, propagating silently into variogram ranges, spatial block
   sizes, GWR bandwidth selection and GP length-scales. Wide extents now receive
   an equal-area projection centred on the data.

5. `make_folds(method = "nndm")` called `set.seed(seed)` unconditionally, and
   `seed` defaults to `NULL`. `set.seed(NULL)` re-initialises the RNG from the
   clock and process ID, so an ordinary call destroyed the caller's random
   number stream. Every seeded path is now guarded by an internal
   `.with_seed()` helper that is a no-op when `seed` is `NULL`.

6. Cross-validation under `parallel = TRUE` was not reproducible; forked workers
   seeded themselves from the current time and process ID.

A third audit pass found a further group of defects that returned a **plausible
wrong number** rather than failing, which is why they are called out here as
well as in `NEWS.md`. In each case the measurement is recorded there.

7. `gwr_model_selection()` ranked models on AIC while labelling the answer
   AICc: GWmodel's diagnostic table is built by `rbind()` over unnamed vectors
   so it never carries column names, making the positional read the normal path
   rather than a fallback, and column 2 is the uncorrected AIC. Executed, the
   old column selected a model containing a pure-noise predictor that AICc drops.

8. The calibrated GP length-scale prior never reached Stan. A `class = "lscale"`
   prior with no `coef` is a *global* prior, which brms applies only to
   coefficients that lack an individual prior — and every `lscale` coefficient
   has one. Confirmed with `make_stancode()`.

9. The GP basis was sized against the per-axis half-range while `brms::gp(c = )`
   multiplies the full pooled range, so the boundary was twice as wide as
   `gp_k` was sized for and the diagnostic meant to catch under-resolution was
   twice too lenient.

10. `fitted()` on a `gwr_fit` returned a local coefficient surface when a
    predictor was named `fit`, `pred`, `prediction`, `fitted` or `yhat`;
    executed in-sample R^2 was −1.18 against a true 0.981.

11. `residual_morans_i()` put weight on a point's own residual whenever
    coordinates were duplicated, and applied an exchangeable null to model
    residuals. Both inflate significance; the corrected moments agree with
    `spdep::lm.morantest()` to machine precision.

12. `summarize_by_cell()`'s design-effect standard errors were too small
    (95% coverage 0.63 at rho = 0.8), and a subsampled variogram design effect
    answered for a cell of `deff_max_n` points rather than the cell's own size.

13. `estimate_sac_range()` swept only two azimuths at ±22.5°, leaving half of
    all directions covered by neither, and its `n_max` subsample was unseeded —
    so the range was irreproducible above `n_max` and the caller's RNG was
    advanced.

A fourth, adversarial pass then set three reviewers the task of breaking the
package on valid input. The corrections it produced are listed under "Fourth
audit pass" in `NEWS.md`; the ones that changed results were the predictor ICC
pooled under shared cell labels (under-estimated by ~1/m), the variogram design
effect evaluated at degree distances on lon/lat input (SEs out by a factor of
354), CRS-less coordinates interpreted one way at fit time and another at
predict time (predictions at the training rows off by one response SD), and
NNDM folds built by an approximation that left the realised distance
distribution up to 0.17 above the target — now the paper's own deterministic
procedure, verified removal for removal against a transcription of it.

A fifth pass set seven reviewers on the package one area each, and a sixth set
eight on it with lenses the earlier passes had not used — differential testing
against reference implementations (`spdep`, `CAST`, `brms`'s own basis
construction, `GWmodel`), invariance under rotation, translation and row order,
mutation testing of the test suite, and a CRAN-policy read. Both are recorded
under their own headings in `NEWS.md`. The sixth pass's corrections that change
results are: `determine_optimal_levels()` returned the elbow's *neighbour* as
its top-ranked candidate (the candidates were sorted ascending, so `k[1]` and
`top_n = 1` were knee − 1 — in 1.0.0 too); `estimate_sac_range()` was not
rotation invariant (its lag cutoff came from the bounding-box diagonal, and the
axis-anchored directional sweep declared anisotropy on an isotropic field in 14
of 18 orientations — the all-pairs fit is now the estimate); `residual_morans_i()`
broke k-nearest-neighbour distance ties by row order, so repeat-visit and
gridded data gave a different statistic in a different row order or with a
different optional package installed (ties are now split); and the GWR
collinearity diagnostic was `kappa()` on the raw matrix at 1e6, a threshold
that depended on the predictors' units (now Belsley's scaled condition index
at 30). The mutation-testing reviewer measured the suite's mutation score at
79% (66 of 84 non-equivalent mutants killed) and named every survivor; each
now has a test that kills it, most of them against an independent reference
rather than the package's own formula.

## Reverse dependencies

None listed on the CRAN page for 1.0.0. Re-confirm with
`revdepcheck::revdep_check()` before submitting.

## Comments for reviewers

* The incoming check's "Possibly misspelled words in DESCRIPTION" NOTE lists
  Delaunay, Voronoi, the surnames Riutort-Mayol, Mila, Pebesma and Strobl, and
  "et al" (reproduced here with CRAN's own aspell configuration: `en_US` plus
  `en_GB` and the `en_stats` dictionary, quoted names and `<doi:...>` targets
  ignored). All are the names of the two tessellations and the authors of the
  five references cited in the Description; every DOI was resolved against
  Crossref and matches its citation.

* `cmdstanr` (Suggests) is not on CRAN; it is available from the repository
  declared in `Additional_repositories` (https://stan-dev.r-universe.dev). It is
  an optional backend for `brms`, used only conditionally via
  `requireNamespace()`; the `rstan` backend that ships with `brms` is used
  otherwise.

* Optional model backends (`GWmodel`, `brms`, `ranger`) and other heavy
  dependencies live in Suggests and are used strictly conditionally. All package
  code, examples, tests **and the vignette** guard their use with
  `requireNamespace()` and skip or degrade gracefully when the package is
  absent. The vignette resolves every optional backend in its setup chunk and
  gates the relevant chunks on the result; the `ggplot2` gate is global, since
  every chunk in it either draws something or feeds something that does, so on a
  machine without `ggplot2` the vignette builds as code without output rather
  than failing `R CMD build`. The vignette built and ran here with `ranger`,
  `gstat`, `geometry` and `patchwork` present and `GWmodel` absent.

* Exactly two examples are wrapped in `\dontrun{}`: `fit_bayesian_spatial_model()`
  and `cv_bayes()`. Both cannot be executed as examples: they compile a Stan
  model, which needs a C++ toolchain (or a CmdStan build) that neither the
  package nor `brms` can supply, and then run minutes of MCMC. The block opens
  with a comment saying so. `brms` itself wraps its fitting examples the same
  way. No example uses `\donttest{}`: every other example runs unconditionally,
  behind a `requireNamespace()` guard where it needs a Suggests package, and
  the slowest of the 42 takes 1.8 s (7.2 s for all of them together, see the
  timings below), so none is close to the 5 s threshold. `checking examples`
  passes.

* Logging writes INFO+ to a session `tempdir()` file and WARN+ to the console
  (see `.onLoad` in `R/zzz.R`), all within a package-specific `logger` namespace
  so the user's global `logger` configuration is never modified. Nothing is
  written outside `tempdir()`. The demo script in `inst/scripts/` also writes
  its PNGs to `tempdir()` unless the user sets an environment variable naming
  somewhere else.

* RNG state is saved and restored around seeded operations via an internal
  `.with_seed()` helper — the per-fold streams used by the parallel
  cross-validation path, the k-means seeding functions, the NNDM and
  subsampled paths. No function body calls `set.seed()` with a literal: every
  seed is an argument the caller can change or set to `NULL`, and
  `fit_gwr_model()`'s local-collinearity spot-check, which used to draw its
  sample under a constant seed, now draws no random numbers at all. Where
  `seed` is `NULL`, nothing is seeded and
  nothing is restored: unseeded functions advance the caller's stream the way
  any other unseeded R function does, rather than re-initialising it. That
  distinction is the subject of fix 5 above.

* Core use is bounded. No function defaults to `parallel::detectCores()`:
  `fit_bayesian_spatial_model(cores = )` and `fit_rf_model(num_threads = )`
  default to `getOption("mc.cores", 1L)`, `predict()` on a random forest passes
  the same, `cv_*(parallel = TRUE)` is capped by that option and by the machine,
  and every worker count is capped at two when `_R_CHECK_LIMIT_CORES_` is set.
  `cv_rf(parallel = )` runs each forked fold's forest on one thread, so a worker
  count is never multiplied by a thread count.

* `NEWS.md` is long. The package changed substantially since the 1.0.0 tag, and
  several of the changes alter results that users may have already reported, so
  each is recorded with enough detail to tell whether it affects a given
  analysis.
