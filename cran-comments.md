# cran-comments

**Status: draft. Nothing is being submitted.** This file describes the working
tree relative to the released 1.0.0 (2026-08-07) so it does not drift out of
date while development continues. Every "<fill in>" below must be replaced with
a real result before this is used for an actual submission.

## Submission

This is an update to spatialkit 1.0.0.

It is a substantial release. The published 1.0.0 carried three defects in the
Bayesian model path, one of which affected every fit; correcting them changes
default results, and correcting them properly required work that grew into a
broader set of changes. `NEWS.md` is the full record; the summary below is what
a reviewer needs.

## Test environments

* local: macOS (aarch64-apple-darwin), R 4.6.x — <fill in: R CMD check result>
* win-builder: R-devel — <fill in>
* win-builder: R-release — <fill in>
* GitHub Actions: macOS / Windows / Ubuntu, R release, devel and oldrel-1 —
  <fill in>

## R CMD check results

<fill in. `devtools::test()` currently reports 856 passing, 0 failures, 1 skip
(a Windows-only fallback path that cannot run on macOS). A full R CMD check has
not been run against the current tree.>

## Breaking changes

Three exported functions have been removed: `evaluate_models()`,
`evaluate_models_cv()` and `phi_prior_bounds()`. All three were thin wrappers
that the package's own documentation described as legacy, and each has a
documented replacement (`compare_models()`, `compare_models_cv()` and
`gp_lengthscale_bounds()` respectively).

Default results from `fit_bayesian_spatial_model()` also change, as a
consequence of the corrections below. This is recorded in `NEWS.md`, and the
previous behaviour remains available by passing `gp_k` and `gp_c` explicitly.

Because of both, the version bump is major rather than minor.

## What was wrong in 1.0.0

Three corrections, in descending order of user impact:

1. The automatic length-scale prior for the Gaussian process was expressed in
   the wrong units. `brms::gp()` defaults to `scale = TRUE`, which rescales its
   covariates so the maximum Euclidean distance between two points is 1 and
   reports the `lscale` parameter in that space. This package standardises the
   coordinates itself and derived the prior in *those* units, so the two
   normalisations differed by roughly the maximum pairwise distance and the
   prior was about five times too diffuse. In practice this produced rejected
   initial values ("Gradient evaluated at the initial value is not finite") on
   most chains. The GP term now sets `scale = FALSE` so there is exactly one
   coordinate scaling.

2. The number of Gaussian process basis functions was chosen from the number of
   observations rather than from the spatial structure of the data. Because
   `brms::gp()` builds a tensor grid, a term `gp(x, y, k = k)` carries `k^2`
   basis functions, and the previous rule reduced to `max(15, floor(sqrt(n)))`
   for any n above 45 — making the basis count identically n. A model at
   n = 10,000 therefore carried 10,000 basis functions, at which point the
   approximation is no longer an approximation. The count is now derived from
   the length-scale-to-domain ratio following Riutort-Mayol et al. (2023,
   Statistics and Computing 33:1).

   Measured on the same simulated field and seeds, 1.0.0 versus the current
   tree: at n = 2,000, `gp_k` 44 to 24 and the basis count 1,936 to 576, with
   cross-validated R-squared 0.9269 to 0.9272 — unchanged within noise — and
   elapsed time 10,908 s to 1,186 s. At n = 300 the derived basis is *larger*
   than before (289 to 529) and the fit correspondingly slower, which is the
   expected consequence of a correction rather than an optimisation.

3. `ensure_projected()` forced data of any extent into a single UTM zone.
   Transverse Mercator scale error grows quadratically with distance from the
   central meridian, so continental-extent data carried distance errors of
   several percent, which propagated silently into variogram ranges, spatial
   block sizes, GWR bandwidth selection and GP length-scales. Wide extents now
   receive an equal-area projection centred on the data.

A fourth fix makes cross-validation reproducible under `parallel = TRUE`;
forked workers previously seeded themselves from the current time and process
ID.

## Scope of the release

Beyond the four items above, `NEWS.md` records roughly 20 further bug fixes and
14 new features. The new user-facing functions are
`area_of_applicability()`, `fit_rf_model()`, `cv_rf()`,
`gwr_model_selection()`, `select_features_forward()`, `predict_surface()`,
`plot_folds()`, and a `plot()` method for `spatial_fit`; `make_folds()` gains
`leave_location_out` and `nndm` methods, and `summarize_by_cell()` gains a
variogram-based design effect.

## Reverse dependencies

<fill in: run revdepcheck::revdep_check(). At 1.0.0 there were none.>

## Comments for reviewers

* Words flagged by the spell checker in DESCRIPTION (Voronoi, Delaunay,
  GWR, backends) are correctly spelled domain terms / standard
  abbreviations.
* 'cmdstanr' (Suggests) is not on CRAN; it is available from the
  repository declared in Additional_repositories
  (https://stan-dev.r-universe.dev). It is an optional backend for 'brms',
  used only conditionally via requireNamespace(); the 'rstan' backend that
  ships with 'brms' is used otherwise.
* Optional model backends ('GWmodel', 'brms', 'ranger') and other heavy
  dependencies live in Suggests and are used strictly conditionally: all
  package code, examples, tests, and the vignette guard their use with
  requireNamespace() and skip or degrade gracefully when the packages are
  unavailable.
* Examples that run full MCMC via 'brms' are wrapped in \dontrun{} because
  they require Stan compilation and long sampling times; fast conditional
  examples for the 'GWmodel' and 'ranger' backends use \donttest{}.
* Logging writes INFO+ to a session tempdir() file and WARN+ to the console
  (see .onLoad in R/zzz.R), all within a package-specific 'logger' namespace so
  the user's global logger configuration is never modified. Nothing is written
  outside tempdir(). RNG state is saved and restored around all seeded
  operations, including the per-fold streams used by the parallel
  cross-validation path and every method of `make_folds()`.
