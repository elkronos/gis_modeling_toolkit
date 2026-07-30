# cran-comments

## Submission

This is a new submission (first release of spatialkit).

## Test environments

* local: macOS (aarch64-apple-darwin23), R 4.6.1 — 0 errors, 0 warnings, 0 notes
* win-builder: R-devel (2026-07-29 r90317 ucrt), Windows Server 2022 x64 —
  1 NOTE (CRAN incoming feasibility: new submission; see comments below)
* win-builder: R-release — 1 NOTE (as above)

## R CMD check results

0 errors | 0 warnings | 1 note

The single NOTE is the standard "CRAN incoming feasibility" NOTE for a new
submission; its contents are addressed in the comments below.

## Comments for reviewers

* Words flagged by the spell checker in DESCRIPTION (Voronoi, Delaunay,
  GWR, backends) are correctly spelled domain terms / standard
  abbreviations.
* 'cmdstanr' (Suggests) is not on CRAN; it is available from the
  repository declared in Additional_repositories
  (https://stan-dev.r-universe.dev). It is an optional backend for 'brms',
  used only conditionally via requireNamespace(); the 'rstan' backend that
  ships with 'brms' is used otherwise.

* Optional model backends ('GWmodel', 'brms') and other heavy dependencies
  live in Suggests and are used strictly conditionally: all package code,
  examples, tests, and the vignette guard their use with requireNamespace()
  and skip or degrade gracefully when the packages are unavailable.
* Examples that run full MCMC via 'brms' are wrapped in \dontrun{} because
  they require Stan compilation and long sampling times; fast conditional
  examples for the 'GWmodel' backend use \donttest{}.
* Logging writes only to tempdir(); RNG state is saved and restored around
  all seeded operations.
