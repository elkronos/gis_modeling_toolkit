# cran-comments

A template plus the notes that recur every submission. Sections marked
`<!-- per release -->` are rewritten each time; everything under **Standing
notes** is true of the package rather than of one release.

The full record for a shipped release is archived under `dev/` — see
`dev/cran-comments-2.0.0.md`, which carries the 2.0.0 submission including what
changed since 1.0.0 and why. `dev/` is in `.Rbuildignore`, so none of it ships.

## What this submission is

<!-- per release -->

## Test environments

<!-- per release: name each machine and its result -->

The set to run, and the order to run it in:

* macOS (aarch64-apple-darwin), R release — `R CMD build .` then
  `_R_CHECK_FORCE_SUGGESTS_=false R CMD check --as-cran --run-donttest
  --timings` on the resulting tarball. Delete any `*.Rcheck` directory and old
  tarball first; `.Rbuildignore` excludes neither.
* win-builder, R-devel and R-release (`devtools::check_win_devel()`,
  `check_win_release()`). Both upload the same filename, so a second submission
  while one is still pending fails with an FTP 550 — that is a duplicate, not a
  problem with the package.
* GitHub Actions `R-CMD-check.yaml`, five jobs: macOS-latest R release;
  windows-latest R release; ubuntu-latest R devel; ubuntu-latest R release;
  ubuntu-latest R oldrel-1.
* GitHub Actions `backends` job (ubuntu-latest, R release) with the optional
  packages installed, so those code paths execute rather than skip.
* GitHub Actions `check-brms.yaml`, **dispatched by hand on the submitted
  commit** — it runs on a schedule and on `workflow_dispatch` only, never on
  push, so "Re-run job" replays the old commit and "Run workflow" is the one
  that matters.

That last job is load-bearing out of proportion to its size: the five Stan
smoke tests in `test-bayes-smoke.R` are gated on `SPATIALKIT_TEST_BRMS`, which
nothing else sets, so it is the only confirmation anywhere that the Bayesian
backend fits a posterior end to end.

## R CMD check results

<!-- per release -->

**Two NOTEs recur and are not the package.**

* *"checking CRAN incoming feasibility" — possibly misspelled words in
  DESCRIPTION.* Mayol, Mila, Pebesma, Riutort and Strobl are surnames of the
  authors of the references cited in the Description, with Riutort-Mayol split
  on its hyphen; "et" and "al" come from "et al." Delaunay and Voronoi appear
  under some aspell configurations. Every DOI was resolved against Crossref and
  matches its citation. Nothing can be done about this from the package side:
  `inst/WORDLIST` affects the `spelling` package, not the aspell run inside
  `R CMD check`.
* *"Package suggested but not available for checking: 'cmdstanr'"*, on any
  machine that cannot reach <https://stan-dev.r-universe.dev>. Where the
  repository is reachable the incoming check resolves it and reports
  `cmdstanr   yes`. See the standing note below.

A checking machine lacking a recent HTML Tidy or `V8` also NOTEs at "checking
HTML version of manual". That is the machine; win-builder checks both manuals
cleanly.

## Reverse dependencies

<!-- per release: re-confirm with revdepcheck::revdep_check(), do not assume -->

## Standing notes

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
  than failing `R CMD build`. It has been built both with every optional
  backend present and with `GWmodel` absent, and re-builds cleanly either way.

* Exactly two examples are wrapped in `\dontrun{}`: `fit_bayesian_spatial_model()`
  and `cv_bayes()`. These are the "missing additional software" case the CRAN
  cookbook gives for `\dontrun{}`: both compile a Stan model, which needs a
  working C++ toolchain (or a CmdStan build) that neither this package nor
  `brms` can supply, and then run minutes of MCMC. Each block opens with a
  comment saying so, and `brms` itself wraps its own fitting examples the same
  way.

  `\donttest{}` would be the wrong tag here rather than a more conservative
  one: `--run-donttest` is exercised on several CRAN platforms, so tagging
  these `\donttest{}` would have CRAN's own machines attempt a Stan
  compilation -- minutes of C++ per example, on shared infrastructure -- which
  is precisely the cost `\dontrun{}` exists to avoid. If the preference is
  nonetheless for `\donttest{}`, say so and it will be changed.

  No other example uses either tag. Every other example runs unconditionally,
  behind a `requireNamespace()` guard where it needs a Suggests package, and
  none comes near the 5 s threshold at which wrapping would be worth
  considering. Quote the `--timings` figures for the release being submitted
  under "R CMD check results" above.

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
