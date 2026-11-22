# cran-comments

## What this submission is

This is an update from 2.0.0 to 3.0.0. It adds a choice of the number of cells
from the data (`resolution_profile()`, `select_resolution()`), pooled cell
values with intervals as the default of `summarize_by_cell()`, density-weighted
seeds as the default partition, and a GWR the package solves itself. NEWS.md
lists the defaults that change results under "Breaking changes".

## Test environments

* macOS (aarch64-apple-darwin23), R 4.6.1, `R CMD check --as-cran --run-donttest
  --timings`: 0 errors, 0 warnings, 1 note. The check takes 2.3 minutes, of
  which the tests take 17 seconds.
* win-builder, R-devel and R-release: <!-- result after the run -->
* macOS builder: unavailable on 21 November 2026 (HTTP 502 on upload, twice).
  macOS is covered by the local check above and the GitHub Actions macOS job.
* GitHub Actions `R-CMD-check.yaml`, five jobs: macOS-latest R release;
  windows-latest R release; ubuntu-latest R devel; ubuntu-latest R release;
  ubuntu-latest R oldrel-1: <!-- result after the run -->
* GitHub Actions `backends` job (ubuntu-latest, R release) with `ggplot2`, `sp`,
  `GWmodel`, `gstat`, `FNN`, `Matrix`, `geometry`, `ranger`, `tibble` and
  `spdep` installed, so those code paths execute rather than skip. It leaves
  out `brms`, `cmdstanr`, `loo` and `patchwork`: <!-- result after the run -->
* GitHub Actions `check-brms.yaml`, dispatched by hand on the submitted commit:
  <!-- result after the run -->

The last job runs on a schedule and on `workflow_dispatch` only, so "Run
workflow" is the step that tests the submitted commit. The Stan tests in 11 test
files, among them the five smoke tests in `test-bayes-smoke.R`, are gated on
`SPATIALKIT_TEST_BRMS`, which nothing else sets, so it is the only CI job that
confirms the Bayesian backend fits a posterior end to end.

## R CMD check results

0 errors | 0 warnings | 1 note on macOS. The note is that this machine lacks a
recent HTML Tidy and `V8` for "checking HTML version of manual"; win-builder
checks both manuals.

**Two check messages recur and are not the package.**

* *"checking CRAN incoming feasibility": possibly misspelled words in
  DESCRIPTION.* Linnenbrink, Mayol, Milà, Pebesma, Riutort and Strobl are
  surnames of the authors of the references cited in the Description, with
  Riutort-Mayol split on its hyphen; "et" and "al" come from "et al.", and
  "lon" from "lon/lat". Delaunay and Voronoi appear under some aspell
  configurations. All six DOIs were resolved against Crossref for this release
  and match their citations. `inst/WORDLIST` affects the `spelling` package and
  leaves the aspell run inside `R CMD check` as it is.
* *"Package suggested but not available for checking: 'cmdstanr'"*, on any
  machine that does not have `cmdstanr` installed, whether or not it can reach
  <https://stan-dev.r-universe.dev>. R 4.6.1 reports this as INFO rather than a
  NOTE. Where the repository is reachable the incoming check also reports
  `cmdstanr   yes`. See the standing note below.

## Check results of the current release

2.0.0 shows an ERROR on r-devel-linux-x86_64-fedora-gcc and fedora-clang and is
OK on the other 11 flavors. In the fedora-gcc log the tests stop after about 51
minutes (91 of CPU) with no failure reported, near the tests of
`estimate_sac_range()`, which calls `gstat`. `gstat` 2.1-6 also shows an ERROR
on both flavors, with one of its own tests stopping after 49 minutes (90 of
CPU). I could not reproduce this on macOS, Windows or Ubuntu.

## Reverse dependencies

spatialkit has no reverse dependencies on CRAN
(`tools::package_dependencies("spatialkit", reverse = TRUE, which = "all")`, 21
November 2026).

## Standing notes

* Test time. Under CRAN settings the tests run the 56 of 306 test files listed
  in `tests/cran-tests.txt`, which together call every exported function. Every
  CI job but ubuntu release sets `NOT_CRAN=true` and runs all 306.

* `cmdstanr` (Suggests) is not on CRAN; it is available from the repository
  declared in `Additional_repositories` (https://stan-dev.r-universe.dev). It is
  an optional backend for `brms`, used only conditionally via
  `requireNamespace()`; the `rstan` backend that ships with `brms` is used
  otherwise.

* The optional model backends (`brms`, `ranger`), the `sp` classes that
  `fit_gwr_model()` needs and other heavy dependencies live in Suggests and are
  used strictly conditionally. `GWmodel` is also in Suggests, but no package
  code calls it: the GWR is implemented here, the tests compare against
  `GWmodel` when it is installed, and the help and messages name its functions.
  All package code, examples, tests **and the nine vignettes** guard their use
  with `requireNamespace()` and skip or degrade gracefully when the package is
  absent. Each vignette resolves the optional packages it uses in its setup
  chunk and gates the relevant chunks on the result. Where every chunk depends
  on a package (`ggplot2` in the North Carolina demo; `gstat` and `ggplot2` in
  the resolution and spatial cross-validation vignettes; `ranger`, `sp`, `gstat`
  and `ggplot2` in the models vignette), the gate is global, so on a machine
  without it the vignette builds as code without output rather than failing
  `R CMD build`. The `R-CMD-check.yaml` matrix builds the vignettes with the
  hard dependencies plus `ggplot2` and `gstat` (and what `gstat` pulls in), so
  the chunks gated on `ranger`, `geometry` or `patchwork` are skipped there, and
  so is the whole models vignette, which needs `ranger`; the pkgdown workflow
  builds them with every package they gate on.

* Exactly four examples are wrapped in `\dontrun{}`:
  `fit_bayesian_spatial_model()`, `cv_bayes()`, `plot_calibration()`, whose
  example has to run `cv_bayes()` to have something to plot, and
  `clear_fitted_cache()`, whose example has to fit a model to have a cache to
  clear. These are the "missing additional software" case the CRAN cookbook
  gives for `\dontrun{}`: all four compile a Stan model, which needs a working
  C++ toolchain (or a CmdStan build) that neither this package nor `brms` can
  supply, and then run minutes of MCMC. Each block opens with a comment saying
  so, and `brms` itself wraps its own fitting examples the same way.

  `\donttest{}` would put a Stan compilation on CRAN's machines, because
  `--run-donttest` is exercised on several CRAN platforms: minutes of C++ per
  example, on shared infrastructure. If the preference is nonetheless for
  `\donttest{}`, say so and it will be changed.

  No other example uses either tag. Every other example runs unconditionally,
  behind a `requireNamespace()` guard where it needs a Suggests package. The
  longest takes about 1 s on macOS (`--timings`).

* Logging writes INFO+ to a session `tempdir()` file and WARN+ to the console
  (see `.onLoad` in `R/zzz.R`), all within a package-specific `logger` namespace
  so the user's global `logger` configuration is never modified. Logging writes
  nothing outside `tempdir()`. The demo script `inst/scripts/example_nc_demo.R`
  writes its PNGs to `tempdir()` unless `SPATIALKIT_DEMO_OUTPUT` names somewhere
  else. The numbered tour scripts draw to the active graphics device, which in a
  non-interactive session is `Rplots.pdf` in the working directory, unless
  `SPATIALKIT_TOUR_OUTPUT` names a folder for their PNGs. Script 07 writes its
  surface files to that folder or to `tempdir()`.

* RNG state is saved and restored around seeded operations via an internal
  `.with_seed()` helper: the per-fold streams used by the parallel
  cross-validation path, the k-means seeding functions, the NNDM and subsampled
  paths. No function body calls `set.seed()` with a literal, and most seeds are
  arguments the caller can change or set to `NULL`. There are two exceptions.
  `fit_bayesian_spatial_model(seed = )` must be a number, because `brms::brm()`
  refuses `NULL`. `gp_lengthscale_bounds()`, which
  `fit_bayesian_spatial_model()` calls, and
  `summarize_by_cell(deff = "variogram")` subsample large inputs under an
  internal seed of 42 that has no argument, and `.with_seed()` restores the
  caller's stream afterwards. `fit_gwr_model()`'s local-collinearity spot-check,
  which used to draw its sample under a constant seed, now draws no random
  numbers at all. Where `seed` is `NULL`, nothing is seeded and nothing is
  restored: unseeded functions advance the caller's stream the way any other
  unseeded R function does, rather than re-initialising it.

* Core use is bounded. No function defaults to `parallel::detectCores()`:
  `fit_bayesian_spatial_model(cores = )` and `fit_rf_model(num_threads = )`
  default to `getOption("mc.cores", 1L)`, `predict()` on a random forest passes
  the same, `cv_*(parallel = TRUE)` is capped by that option and by the machine,
  and every forked worker count is capped at two when `_R_CHECK_LIMIT_CORES_` is
  set. That cap does not apply to `cores` or `num_threads`, which take the
  caller's value or `mc.cores`. `cv_rf(parallel = )` runs each forked fold's
  forest on one thread unless the caller passes `num_threads`.
  `cv_bayes(parallel = )` does not reduce the sampler's cores, so with
  `options(mc.cores = 4)` and `parallel = 2` each of the two workers uses four
  cores. With `mc.cores` unset, each uses one.
