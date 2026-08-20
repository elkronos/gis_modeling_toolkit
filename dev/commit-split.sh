#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Split the post-1.0.0 working tree into reviewable commits.
#
# Run from the package root:   bash dev/commit-split.sh
#
# Notes before you run it:
#   * It creates a branch first; main is left where it is.
#   * Each group is skipped if nothing is staged for it, so a path that does
#     not exist (or was never modified) is harmless.
#   * Only the FINAL commit is a verified state. The intermediate ones are
#     grouped for reviewability, not bisectability -- a test file may land one
#     commit before or after the code it exercises. If you need every commit
#     green, squash instead.
#   * Nothing here contacts CRAN or pushes anywhere.
# ---------------------------------------------------------------------------
set -euo pipefail

BRANCH="${1:-dev-2.0.0}"

if [ ! -f DESCRIPTION ]; then
  echo "Run this from the package root." >&2; exit 1
fi

echo "=== working tree before anything is staged ==="
git status --porcelain
echo
read -r -p "Proceed? [y/N] " reply
case "$reply" in y|Y) ;; *) echo "aborted."; exit 0 ;; esac
echo

current=$(git rev-parse --abbrev-ref HEAD)
echo "current branch: $current"
if [ "$current" != "$BRANCH" ]; then
  git checkout -b "$BRANCH"
fi

TRAILER=$'\n\nCo-Authored-By: Claude Opus 5 <noreply@anthropic.com>\nClaude-Session: https://claude.ai/code/session_014DkVoswUXWntSsYAJsxgPw'

commit_group () {
  local subject="$1"; shift
  local body="$1"; shift
  git add -A -- "$@" 2>/dev/null || true
  if git diff --cached --quiet; then
    echo "  (nothing staged) skip: $subject"
    return 0
  fi
  printf '%s\n\n%s%s\n' "$subject" "$body" "$TRAILER" | git commit -q -F -
  echo "  committed: $subject"
}

# --- 1 ---------------------------------------------------------------------
commit_group \
"Correct the GP length-scale prior and derive the basis count" \
"brms::gp() defaults to scale = TRUE, which rescales its covariates so the
maximum pairwise distance is 1 and reports lscale in that space. The package
standardised coordinates itself and derived the prior in those units, so the
two normalisations differed by roughly the maximum pairwise distance and the
prior was about five times too diffuse -- producing rejected initial values on
most chains. The GP term now sets scale = FALSE so there is one scaling.

Separately, gp_k was chosen from n rather than from spatial structure. Because
gp() builds a tensor grid, gp(x, y, k = k) carries k^2 basis functions, and the
old rule reduced to max(15, floor(sqrt(n))) -- making the basis count
identically n. It is now derived from the length-scale-to-domain ratio
following Riutort-Mayol et al. (2023).

Measured on the same field and seeds, n = 2000: gp_k 44 -> 24, basis
1936 -> 576, cross-validated R-squared 0.9269 -> 0.9272 (unchanged within
noise), elapsed 10908s -> 1186s. At n = 300 the basis is larger than before
(289 -> 529) and the fit slower: a correction, not an optimisation." \
  R/model-bayesian.R R/model-prep.R \
  tests/testthat/test-gp-basis.R tests/testthat/test-lscale-prior.R

# --- 2 ---------------------------------------------------------------------
commit_group \
"Stop forcing continental extents into a single UTM zone" \
"Transverse Mercator scale error grows quadratically with distance from the
central meridian, so data spanning the contiguous US in one zone carried
distance errors of several percent. Those errors were silent and propagated
into variogram ranges, block sizes, GWR bandwidth and GP length-scales.

Extents reaching more than 5 degrees from the candidate zone's central meridian
now get an equal-area projection centred on the data. The trigger measures
longitude offset only: latitude span does not drive TM error, and cos(lat)
shrinks it, so a tall narrow strip is fine under UTM and worse under LAEA.

Measured over CONUS: median distance error 1.525% -> 0.249%, max 11.68% ->
2.43%." \
  R/crs-geometry.R \
  tests/testthat/test-crs-selection.R tests/testthat/test-crs-projection.R \
  dev/verify-crs-distortion.R

# --- 3 ---------------------------------------------------------------------
commit_group \
"Make parallel CV reproducible; add fold methods and SAC guards" \
".cv_run_folds() called mclapply() without seeding the fork streams, so each
worker seeded itself from the time and PID. One seed per fold is now drawn in
the parent, making each fold's stream a function of (seed, fold index) alone.

make_folds() gains leave_location_out (repeated measurements at a site stay in
one fold) and nndm (Milà et al. 2022 distance matching, sizing each exclusion
so training-to-test distances reproduce the distances from actual prediction
locations).

estimate_sac_range() now rejects a fitted range that exceeds the longest lag
the variogram was fitted over -- such a range is unidentified rather than long,
and block sizing from it yields one block covering everything.

Also fixed here, from the full-tree audit: an on.exit() name collision that
silently replaced the caller's RNG stream; a %d format applied to a median that
is not always integral; MULTIPOINT accepted without coercion, so per-vertex
st_coordinates() rows misaligned every fold; a single block producing an empty
training set; mclapply() try-error objects surviving the NULL filter; and
cv_spatial() staying silent when every fold failed." \
  R/cross-validation.R \
  tests/testthat/test-cv-parallel.R tests/testthat/test-fold-methods.R \
  tests/testthat/test-sac-range.R tests/testthat/test-fold-construction.R \
  tests/testthat/test-make-folds-row-ids.R tests/testthat/test-core-count.R \
  tests/testthat/helper-bootlm.R tests/testthat/helper-logging.R

# --- 4 ---------------------------------------------------------------------
commit_group \
"Add a variogram design effect; make the kNN fallback reachable" \
"summarize_by_cell() gains deff = 'variogram', computing a per-cell design
effect from a fitted variogram rather than one pooled intra-class correlation.
For n points with correlation matrix R the effective sample size of the mean is
n^2 / sum(R), so deff = sum(R)/n. This generalises the Kish option -- a constant
off-diagonal correlation recovers 1 + (n-1)rho exactly -- but lets correlation
decay with distance, which is the point of having fitted a variogram.

.build_knn_weights() takes its backend availability as parameters rather than
calling requireNamespace() directly, so the dense fallback and its size guard
can be exercised on a machine that has FNN installed. The test for that
fallback had been skipping silently." \
  R/evaluation.R R/assignment.R \
  tests/testthat/test-knn-weights.R tests/testthat/test-deff-variogram.R \
  tests/testthat/test-summarize-by-cell.R tests/testthat/test-morans-weights.R \
  tests/testthat/test-morans-variance.R

# --- 5 ---------------------------------------------------------------------
commit_group \
"Add predict_surface(), plot() for spatial_fit, and plot_folds()" \
"predict() required newdata to be built by hand, which made producing a map --
the most common thing wanted from a fitted spatial model -- more work than it
should be. predict_surface() builds the grid, joins covariates by nearest
feature, clips to a boundary, predicts in chunks and returns sf. Chunking
matters for bayesian_fit, where the draw matrix is n_draws x n_newdata.

plot() on a spatial_fit gives residuals mapped at the training locations,
observed-vs-predicted, and the residual variogram with the fitted model and
effective range overlaid -- so the variogram fit can be judged rather than
trusted. plot_folds() maps a fold scheme, which is the fastest way to see
whether blocks actually separate the data or are smaller than the
autocorrelation range and therefore leaking." \
  R/predict-surface.R R/plotting-fits.R \
  tests/testthat/test-predict-surface.R tests/testthat/test-plotting-fits.R \
  tests/testthat/helper-lmfit.R

# --- 6 ---------------------------------------------------------------------
commit_group \
"Add forward feature selection and GWR model selection" \
"select_features_forward() selects predictors against a spatially blocked
inner-fold score. The blocking is the whole justification: random inner folds
inside blocked outer folds select variables that look predictive only because
nearby points leak between train and test, and the outer loop then reports
honest-looking numbers for a dishonestly chosen feature set. method therefore
defaults to block_kfold and warns otherwise.

gwr_model_selection() wraps GWmodel::gwr.model.selection() and returns a ranked
table instead of two loosely-coupled lists. It is the fast in-sample
counterpart -- same forward search, scored by AICc rather than by a blocked
estimate. Both limitations are documented rather than hidden: one bandwidth is
shared across candidate models, and the null model is never evaluated, so the
result always names at least one predictor.

GWmodel does not label its diagnostic table, so the AICc column is read
positionally (as GWmodel's own documentation does) and the result records which
way it was found. Its progress output is discarded by default: it is written
with bare cat(), so no suppressMessages() touches it, and it scales with the
square of the candidate count." \
  R/feature-selection.R R/model-selection-gwr.R \
  tests/testthat/test-feature-selection.R \
  tests/testthat/test-gwr-model-selection.R

# --- 7 ---------------------------------------------------------------------
commit_group \
"Add area_of_applicability()" \
"Implements the dissimilarity index and area of applicability of Meyer &
Pebesma (2021). A fitted model returns a number for any location, including
locations whose predictors look nothing like anything it trained on. Those
predictions are extrapolations dressed as interpolations, and a
cross-validation score says nothing about them -- the held-out folds came from
the same predictor distribution as the training data.

Predictors are centred and scaled with the TRAINING data's statistics, never
the prediction data's, which would re-centre a far-away block and make it look
familiar. Weighting follows CAST, the reference implementation: the scaled
column is multiplied by the weight, contributing w^2 to the squared distance.

Pass the folds you actually validated with. Without them the reference is each
point's nearest neighbour anywhere, which for clustered data is very close,
giving a conservative AOA. Buffered and NNDM folds contribute the training set
they actually left available rather than 'everything outside the fold', so the
exclusion they exist to enforce is not silently undone." \
  R/area-of-applicability.R \
  tests/testthat/test-area-of-applicability.R

# --- 8 ---------------------------------------------------------------------
commit_group \
"Add a ranger random forest backend" \
"fit_rf_model() and cv_rf() return an rf_fit that works with cv_spatial(),
predict_surface(), area_of_applicability() and plot() like any other backend.
Three choices are opinionated because each is where an RF spatial model usually
goes wrong.

include_coords defaults to FALSE. Giving a forest x and y lets it reproduce the
training surface by memorising location, then fail badly off it; random CV does
not catch that, which is how the practice became common (Meyer et al. 2019).

fitted() returns out-of-bag predictions, following the RF packages' own
convention. It matters because summary() calls fitted() to compute R-squared,
and in-sample forest predictions are near-memorisation, so the alternative is a
summary reporting a fictitious number. The cost is that summary() means
something different here than for a gwr_fit; compare_models_cv() is the
like-for-like path.

The OOB error is reported but labelled everywhere as a random hold-out, and so
optimistic under spatial autocorrelation for the same reason random k-fold is.
Importance defaults to permutation rather than impurity, which is biased toward
continuous and high-cardinality predictors (Strobl et al. 2007)." \
  R/model-rf.R tests/testthat/test-model-rf.R

# --- 9 ---------------------------------------------------------------------
commit_group \
"Regroup the test suite by subject" \
"test-fixes-misc.R and test-fixes-round2.R were named after the review round
that produced them; test-core.R held CRS handling, regression metrics, grid
geometry, fold construction, cell summaries, Moran's I weights and GWR
bandwidth in one file. Renaming alone would have produced a different lie, so
the 36 tests are regrouped by subject into files whose names describe their
contents, each opening with a comment on why the group exists.

No test was added or removed in the move: 36 in, 36 out." \
  tests/testthat/

# --- 10 --------------------------------------------------------------------
commit_group \
"Harden the dev scripts against measuring the installed package" \
"library(spatialkit) loads the INSTALLED package. devtools::test() loads from
source, so a working tree can be many changes ahead and a dev script using
library() silently measures the wrong code -- which is what happened, costing a
five-hour baseline run that reported gp_k = 44 at n = 2000 (exactly the rule
this branch replaces).

All dev scripts now load via pkgload::load_all(), print the namespace path they
loaded from, and say whether it is the working tree. The namespace path is the
ground truth; packageVersion() is not, since after load_all() it may still
report the installed DESCRIPTION.

dev/baseline-accuracy.R now checks its acceptance criteria at runtime instead
of leaving them in a trailing comment -- it detects gp_k == floor(sqrt(n))
specifically, because those numbers look entirely plausible otherwise.
dev/check-gp-live.R is new: one small fit, about two minutes, answering whether
the GP fix is live in the fit path before committing hours to the full
baseline." \
  dev/ .github/

# --- 11 --------------------------------------------------------------------
commit_group \
"Update generated docs, DESCRIPTION and release notes" \
"Regenerate NAMESPACE and man/ for the new exports and S3 methods. Add ranger
to Suggests and raise the testthat floor to 3.1.5 for expect_no_error().

NEWS.md documents everything relative to 1.0.0, the released version. Its
heading carries a parseable version because R's NEWS reader keys entries off
it: a heading like '(development version)' makes R CMD check report 'No news
entries found'.

cran-comments.md is marked draft and rewritten against 1.0.0; the previous text
described a 1.1.0 submission fixing three defects, and the tree now holds around
thirty changes. Everything unverified is an explicit placeholder so it cannot
quietly become a false claim.

R CMD check: 0 errors, 0 warnings, 0 notes on R 4.6.1 / macOS.
devtools::test(): 856 passing, 1 skip." \
  NAMESPACE man/ DESCRIPTION NEWS.md cran-comments.md README.md

# --- anything left ---------------------------------------------------------
# Deliberately NOT auto-committed: this script was written without sight of the
# working tree, so a catch-all `git add -A .` could sweep in .Rhistory, .rds
# baselines, .DS_Store or anything else not covered by .gitignore.
leftover=$(git status --porcelain)
if [ -n "$leftover" ]; then
  echo
  echo "=== NOT committed -- review and handle these yourself ==="
  printf '%s\n' "$leftover"
fi

echo
echo "done. review with:  git log --oneline ${current}..HEAD"
echo "undo everything with: git checkout ${current} && git branch -D ${BRANCH}"
