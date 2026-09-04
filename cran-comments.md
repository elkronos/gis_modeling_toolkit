# cran-comments

**Status: draft. Nothing is being submitted.** This file describes the working
tree so it does not drift out of date while development continues. Every
"<fill in>" below must be replaced with a real result before this is used for an
actual submission.

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
not install in this container. And although `brms` itself is installed, **Stan
could not compile here**, so no posterior was ever sampled: the brms findings
were established with `make_stancode()`, `make_standata()`, `get_prior()` and
`validate_prior()`, which need no compilation. **Every GWR and brms result
recorded below must be re-confirmed on a machine with the real packages before
submission.** Not installed: `cmdstanr`, `loo`.

Still to run before any submission:

* local: macOS (aarch64-apple-darwin), R release — <fill in>
* win-builder: R-devel — <fill in>
* win-builder: R-release — <fill in>
* GitHub Actions `R-CMD-check.yaml`, five jobs (not a cross-product):
  macOS-latest R release; windows-latest R release; ubuntu-latest R devel;
  ubuntu-latest R release; ubuntu-latest R oldrel-1. R-devel and oldrel-1 are
  tested on Linux only. — <fill in>
* GitHub Actions `backends` job (ubuntu-latest, R release) with `sp`, `GWmodel`,
  `gstat`, `FNN`, `Matrix`, `geometry`, `ranger` and `tibble` installed, so the
  optional code paths actually execute rather than skip. — <fill in>
* GitHub Actions `check-brms.yaml` (weekly, ubuntu-latest, R release) with
  `brms` and the Stan toolchain, which is the only job that sets
  `SPATIALKIT_TEST_BRMS` and therefore the only one that runs the five Stan
  smoke tests. — <fill in>

## R CMD check results

<fill in: `R CMD check --as-cran` has not been run against the current tree on a
release-R machine, and neither has the manual (no LaTeX here).>

On the environment above, `R CMD check --no-manual` on the built tarball --
with the vignette built, the `GWmodel` stub and `brms` installed, and
`_R_CHECK_CRAN_INCOMING_=false` because `--as-cran` needs network access to
CRAN that this container does not have -- reports **1 NOTE and nothing else**.

* **NOTE -- "Package suggested but not available for checking: 'cmdstanr'".**
  Environmental, and expected: `cmdstanr` is not on CRAN. It is reached through
  `Additional_repositories`, is used strictly conditionally via
  `requireNamespace()`, and the `check-brms` CI job installs it.

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

`R CMD build` produces a 994 KB tarball, of which the built vignette HTML is the
bulk. `checking running R code from vignettes` passes.

`testthat` reports **3479 passing, 0 failures, 0 errors, 0 warnings, 9 skips**
with `NOT_CRAN=true` and both backends present, and no test runs with zero
assertions. The nine skips are:

* 5 Stan smoke tests in `test-bayes-smoke.R`, skipped because
  `SPATIALKIT_TEST_BRMS` is unset -- `skip_if_not_installed("brms")` alone was
  not enough, since these compile Stan models and would otherwise run in any
  matrix job that happened to have `brms`
* 1 Windows-only fallback path that cannot run on Linux
* 3 that skip *because* an optional backend is installed: they assert the
  behaviour seen when `GWmodel` or `brms` is absent

Under `R CMD check`, where `NOT_CRAN` is unset, four further tests skip on
purpose -- the `parallel::mclapply()` fork tests in `test-cv-parallel.R`, which
are `skip_on_cran()`.

`README.md` does not restate these counts, precisely so the two cannot drift
apart.

## Breaking changes since the 1.0.0 tag

Three exported functions have been removed: `evaluate_models()`,
`evaluate_models_cv()` and `phi_prior_bounds()`. All three were thin wrappers
the package's own documentation described as legacy, and each has a documented
replacement (`compare_models()`, `compare_models_cv()` and
`gp_lengthscale_bounds()`). `NAMESPACE` now exports 40 objects.

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

`NEWS.md` is the full record: 20 breaking changes, 62 bug fixes, 13 new
features and 19 documentation entries relative to 1.0.0.

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
   ratio following Riutort-Mayol et al. (2023, Statistics and Computing 33:1).
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

## Reverse dependencies

None listed on the CRAN page for 1.0.0. Re-confirm with
`revdepcheck::revdep_check()` before submitting.

## Comments for reviewers

* Words flagged by the spell checker in DESCRIPTION (Voronoi, Delaunay, GWR,
  backends) are correctly spelled domain terms or standard abbreviations.

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
  and `cv_bayes()`, both of which run full MCMC via `brms` and require Stan
  compilation. Every other example runs, using `\donttest{}` plus a
  `requireNamespace()` guard where it needs an optional backend. `checking
  examples` passes.

* Logging writes INFO+ to a session `tempdir()` file and WARN+ to the console
  (see `.onLoad` in `R/zzz.R`), all within a package-specific `logger` namespace
  so the user's global `logger` configuration is never modified. Nothing is
  written outside `tempdir()`. The demo script in `inst/scripts/` also writes
  its PNGs to `tempdir()` unless the user sets an environment variable naming
  somewhere else.

* RNG state is saved and restored around seeded operations via an internal
  `.with_seed()` helper — the per-fold streams used by the parallel
  cross-validation path, the k-means seeding functions, and `fit_gwr_model()`'s
  local-collinearity spot-check. Where `seed` is `NULL`, nothing is seeded and
  nothing is restored: unseeded functions advance the caller's stream the way
  any other unseeded R function does, rather than re-initialising it. That
  distinction is the subject of fix 5 above.

* `NEWS.md` is long. The package changed substantially since the 1.0.0 tag, and
  several of the changes alter results that users may have already reported, so
  each is recorded with enough detail to tell whether it affects a given
  analysis.
