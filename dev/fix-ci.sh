#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Rewrite the two CI workflows.  The file bridge cannot write .github/workflows
# (workflow files execute code, so remote writes to them are blocked), hence a
# script.
#
#   bash dev/fix-ci.sh
#
# Changes, both files:
#   * setup-r gains extra-repositories for the stan-dev universe. cmdstanr is
#     in Suggests and is NOT on CRAN; R CMD check honours DESCRIPTION's
#     Additional_repositories but pak does not, so dependency resolution failed
#     with "Can't find package called cmdstanr" before any R code ran.
#   * _R_CHECK_FORCE_SUGGESTS_: false, so an optional backend that will not
#     install on one platform skips its tests instead of failing the run.
#   * actions/checkout v4 -> v5 (Node 20 deprecation warning).
# ---------------------------------------------------------------------------
set -euo pipefail
[ -f DESCRIPTION ] || { echo "Run from the package root." >&2; exit 1; }
mkdir -p .github/workflows

cat > .github/workflows/R-CMD-check.yaml <<'SKEOF'
# Workflow derived from https://github.com/r-lib/actions/tree/v2/examples
on:
  push:
    branches: [main, master]
  pull_request:
  workflow_dispatch:

name: R-CMD-check

permissions: read-all

jobs:
  R-CMD-check:
    runs-on: ${{ matrix.config.os }}
    name: ${{ matrix.config.os }} (${{ matrix.config.r }})

    strategy:
      fail-fast: false
      matrix:
        config:
          - {os: macos-latest,   r: 'release'}
          - {os: windows-latest, r: 'release'}
          - {os: ubuntu-latest,  r: 'devel', http-user-agent: 'release'}
          - {os: ubuntu-latest,  r: 'release'}
          - {os: ubuntu-latest,  r: 'oldrel-1'}

    env:
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}
      R_KEEP_PKG_SOURCE: yes
      # Every Suggests use in this package is requireNamespace()-guarded,
      # so an optional backend that will not install on one platform should
      # skip its tests, not fail the run.
      _R_CHECK_FORCE_SUGGESTS_: false

    steps:
      - uses: actions/checkout@v5

      - uses: r-lib/actions/setup-pandoc@v2

      - uses: r-lib/actions/setup-r@v2
        with:
          r-version: ${{ matrix.config.r }}
          http-user-agent: ${{ matrix.config.http-user-agent }}
          use-public-rspm: true
          # Load-bearing: cmdstanr is in Suggests and is NOT on CRAN. R CMD
          # check honours DESCRIPTION's Additional_repositories, but pak --
          # which setup-r-dependencies uses to resolve Suggests -- does not,
          # so without this the dependency graph fails to solve with
          # "Can't find package called cmdstanr" before any R code runs.
          extra-repositories: 'https://stan-dev.r-universe.dev'

      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          extra-packages: any::rcmdcheck
          needs: check

      - uses: r-lib/actions/check-r-package@v2
        with:
          upload-snapshots: true
          build_args: 'c("--no-manual","--compact-vignettes=gs+qpdf")'

  # ---------------------------------------------------------------------------
  # Optional backends live in Suggests and are guarded by requireNamespace().
  # Without this job the matrix above skips most model tests, so a green
  # matrix would prove very little.  GWmodel / gstat / FNN / Matrix are cheap
  # enough to install on every push; brms + Stan is not (see nightly workflow).
  # ---------------------------------------------------------------------------
  backends:
    runs-on: ubuntu-latest
    name: ubuntu-latest (release, with optional backends)

    env:
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}
      R_KEEP_PKG_SOURCE: yes
      # Every Suggests use in this package is requireNamespace()-guarded,
      # so an optional backend that will not install on one platform should
      # skip its tests, not fail the run.
      _R_CHECK_FORCE_SUGGESTS_: false

    steps:
      - uses: actions/checkout@v5

      - uses: r-lib/actions/setup-pandoc@v2

      - uses: r-lib/actions/setup-r@v2
        with:
          r-version: 'release'
          use-public-rspm: true
          # Load-bearing: cmdstanr is in Suggests and is NOT on CRAN. R CMD
          # check honours DESCRIPTION's Additional_repositories, but pak --
          # which setup-r-dependencies uses to resolve Suggests -- does not,
          # so without this the dependency graph fails to solve with
          # "Can't find package called cmdstanr" before any R code runs.
          extra-repositories: 'https://stan-dev.r-universe.dev'

      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          extra-packages: |
            any::rcmdcheck
            any::sp
            any::GWmodel
            any::gstat
            any::FNN
            any::Matrix
            any::geometry
          needs: check

      - name: Confirm optional backends are installed
        run: |
          for (p in c("sp", "GWmodel", "gstat", "FNN", "Matrix", "geometry")) {
            if (!requireNamespace(p, quietly = TRUE))
              stop("Backend package not installed: ", p)
          }
          cat("All optional backends present.\n")
        shell: Rscript {0}

      - uses: r-lib/actions/check-r-package@v2
        with:
          upload-snapshots: true
          build_args: 'c("--no-manual","--compact-vignettes=gs+qpdf")'
SKEOF

cat > .github/workflows/check-brms.yaml <<'SKEOF'
# Weekly check with the brms / Stan backend installed.
#
# brms pulls in a Stan toolchain and compiles models, which takes minutes per
# run -- too slow for every push.  The main R-CMD-check workflow therefore
# skips every brms-guarded code path.  This workflow closes that gap on a
# schedule, and can be triggered by hand before a CRAN submission.
on:
  schedule:
    - cron: '0 6 * * 1'   # Mondays 06:00 UTC
  workflow_dispatch:

name: check-brms

permissions: read-all

jobs:
  check-brms:
    runs-on: ubuntu-latest
    name: ubuntu-latest (release, with brms)

    env:
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}
      R_KEEP_PKG_SOURCE: yes
      # Every Suggests use in this package is requireNamespace()-guarded,
      # so an optional backend that will not install on one platform should
      # skip its tests, not fail the run.
      _R_CHECK_FORCE_SUGGESTS_: false

    steps:
      - uses: actions/checkout@v5

      - uses: r-lib/actions/setup-pandoc@v2

      - uses: r-lib/actions/setup-r@v2
        with:
          r-version: 'release'
          use-public-rspm: true
          # Load-bearing: cmdstanr is in Suggests and is NOT on CRAN. R CMD
          # check honours DESCRIPTION's Additional_repositories, but pak --
          # which setup-r-dependencies uses to resolve Suggests -- does not,
          # so without this the dependency graph fails to solve with
          # "Can't find package called cmdstanr" before any R code runs.
          extra-repositories: 'https://stan-dev.r-universe.dev'

      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          extra-packages: |
            any::rcmdcheck
            any::brms
            any::loo
            any::sp
            any::GWmodel
            any::gstat
            any::FNN
            any::Matrix
          needs: check

      - name: Confirm brms is installed
        run: |
          if (!requireNamespace("brms", quietly = TRUE))
            stop("brms not installed; this workflow has no purpose without it.")
          cat("brms ", as.character(utils::packageVersion("brms")), "\n")
        shell: Rscript {0}

      - uses: r-lib/actions/check-r-package@v2
        with:
          upload-snapshots: true
          build_args: 'c("--no-manual","--compact-vignettes=gs+qpdf")'
SKEOF

echo "wrote both workflows."
git --no-pager diff --stat -- .github/workflows/ || true
