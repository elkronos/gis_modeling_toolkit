# cran-comments

**Status: draft. Nothing is being submitted.** This file describes the working
tree so it does not drift out of date while development continues. Every
"<fill in>" below must be replaced with a real result before this is used for an
actual submission.

## What this submission is

**This is an update to `spatialkit` 1.0.0, published on CRAN 2026-08-07.**

The version to submit is **not** `1.0.0.9000` — that is the development string.
Pick a release number before submitting and set `Version:` in `DESCRIPTION`
accordingly. A **major** bump to `2.0.0` is the honest choice: three exported
functions are removed, and several defaults change the result of a fit or a
comparison (see "Breaking changes" in `NEWS.md`). Update this file's version
references at the same time.

### Reverse dependencies

The CRAN page for 1.0.0 lists no reverse depends, imports or suggests, so there
is nothing to break. Confirm with `revdepcheck::revdep_check()` before
submitting rather than relying on this note.

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
roxygen2 7.3.1, knitr 1.45, rmarkdown 2.25. **Not installed: `GWmodel`, `brms`,
`cmdstanr`, `loo`.**

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
  `brms` and the Stan toolchain. — <fill in>

## R CMD check results

<fill in: `R CMD check --as-cran` has not been run against the current tree on a
release-R machine, and neither has the manual (no LaTeX here).>

On the environment above, `R CMD check --no-manual --no-build-vignettes` on the
built tarball reports **2 WARNINGs, 1 NOTE**. All three are accounted for:

* **NOTE — "Packages suggested but not available for checking: 'GWmodel',
  'brms', 'cmdstanr', 'loo'".** Environmental. These are the optional backends
  this machine does not have; the `backends` and `check-brms` CI jobs exist to
  cover them.

* **WARNING — "checking R files for non-ASCII characters": `crs-geometry.R`.**
  Real, and must be fixed before submission. `R/crs-geometry.R` line 188 carries
  a literal U+2014 EM DASH inside a string literal in the `ensure_projected()`
  error message:

  ```
  "object that carries a CRS — or omit `target_crs` to let a ",
  ```

  Replace it with `--` or the escape `—`. Em dashes elsewhere in `R/` sit
  in comments and roxygen blocks, which the check tolerates; this is the only
  one in code.

* **WARNING — "checking R files for syntax errors":** the only content is
  `Warning in Sys.setlocale("LC_CTYPE", "en_US.UTF-8") : OS reports request to
  set locale to "en_US.UTF-8" cannot be honored`. Environmental: this container
  ships only the C locale, so the check's own locale switch fails. No syntax
  error was reported. Expect this WARNING to disappear on any machine with a
  UTF-8 locale.

`R CMD build` produces a 994 KB tarball, of which the built vignette is 857 KB.
Both `checking tests` and `checking running R code from vignettes` pass.

`testthat` reports **1752 passing, 0 failures, 0 errors, 7 skips** on that
environment (`NOT_CRAN=true`). The seven skips are:

* 5 skipped for `GWmodel` (`skip_if_not_installed`), which is not installed here
* 1 skipped for `brms`, likewise
* 1 Windows-only fallback path that cannot run on Linux

With `GWmodel` and `brms` installed there would be none. Without `NOT_CRAN` set,
four further tests skip on purpose: they exercise `parallel::mclapply()` fork
behaviour, which is not appropriate to run on CRAN.

This figure and the one in `README.md` describe the same run; the README does
not restate the count, precisely so the two cannot drift apart.

## Breaking changes since the 1.0.0 tag

Three exported functions have been removed: `evaluate_models()`,
`evaluate_models_cv()` and `phi_prior_bounds()`. All three were thin wrappers
the package's own documentation described as legacy, and each has a documented
replacement (`compare_models()`, `compare_models_cv()` and
`gp_lengthscale_bounds()`).

Default results from `fit_bayesian_spatial_model()` also change, as a
consequence of three corrections to the Gaussian process path (the basis count,
the coordinate scaling, and the length-scale prior). The previous behaviour
remains reachable by passing `gp_k`, `gp_c` and `gp_iso` explicitly. `NEWS.md`
is the full record.

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
   Measured on the recorded baselines: at n = 2,000, `gp_k` 44 to 24 and the
   basis count 1,936 to 576; at n = 10,000, 100 to 23 and 10,000 to 529; at
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
   number stream.

6. Cross-validation under `parallel = TRUE` was not reproducible; forked workers
   seeded themselves from the current time and process ID.

`NEWS.md` records roughly 100 further user-facing fixes and 14 new features.

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
  than failing `R CMD build`.

* Exactly two examples are wrapped in `\dontrun{}`: `fit_bayesian_spatial_model()`
  and `cv_bayes()`, both of which run full MCMC via `brms` and require Stan
  compilation. Every other example runs, using `\donttest{}` plus a
  `requireNamespace()` guard where it needs an optional backend.

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
