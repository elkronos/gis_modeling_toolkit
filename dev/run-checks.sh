#!/usr/bin/env bash
# dev/run-checks.sh -- run the post-edit verification sequence.
#
#   cd /Users/appleair/Downloads/gis_modeling_toolkit-main
#   bash dev/run-checks.sh
#
# Stops at the first failure.  Each stage is cheap-to-expensive, so a problem
# surfaces as early as possible.
set -euo pipefail

PKG_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PKG_DIR"
echo "package: $PKG_DIR"

command -v Rscript >/dev/null 2>&1 || {
  echo "ERROR: Rscript not on PATH. Install R from https://cran.r-project.org"
  exit 1
}
Rscript -e 'cat("R", as.character(getRversion()), "\n")'

# --- 0. remove the staging folder --------------------------------------------
# The remote session could not delete files, so cleaned-up junk was moved here.
# R CMD check flags a non-standard top-level directory, so it must go first.
if [ -d "_to_delete" ]; then
  echo "==> removing _to_delete/ ($(ls _to_delete | wc -l | tr -d ' ') files)"
  rm -rf _to_delete
else
  echo "==> _to_delete/ already gone"
fi

# --- 1. dependencies ---------------------------------------------------------
echo "==> checking required packages"
Rscript -e '
need <- c("devtools", "roxygen2", "testthat", "sf", "dplyr", "logger", "digest")
opt  <- c("sp", "GWmodel", "gstat", "FNN", "Matrix", "geometry",
          "ggplot2", "patchwork", "tibble", "knitr", "rmarkdown", "loo")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  cat("MISSING (required):", paste(miss, collapse = ", "), "\n")
  cat("install with: install.packages(c(",
      paste(sprintf("\"%s\"", miss), collapse = ", "), "))\n", sep = "")
  quit(status = 1)
}
skip <- opt[!vapply(opt, requireNamespace, logical(1), quietly = TRUE)]
if (length(skip))
  cat("NOTE: optional packages absent, related tests will skip:",
      paste(skip, collapse = ", "), "\n")
cat("required packages present\n")
'

# --- 2. regenerate docs ------------------------------------------------------
# NAMESPACE and man/*.Rd are generated.  Phases 1-4 changed roxygen in
# model-bayesian.R, model-prep.R and crs-geometry.R, and removed three exported
# functions, so this must run before anything else reads them.
echo "==> devtools::document()"
Rscript -e 'devtools::document(roclets = c("rd", "collate", "namespace"))'

echo "==> confirming the removed exports are gone"
Rscript -e '
ns <- readLines("NAMESPACE")
bad <- grep("evaluate_models|phi_prior_bounds", ns, value = TRUE)
if (length(bad)) { cat("STILL EXPORTED:\n"); cat(bad, sep = "\n"); quit(status = 1) }
rd <- list.files("man", pattern = "[.]Rd$", full.names = TRUE)
stale <- rd[vapply(rd, function(f)
  any(grepl("phi_prior_bounds", readLines(f, warn = FALSE))), logical(1))]
if (length(stale)) {
  cat("STALE .Rd still references phi_prior_bounds:\n"); cat(stale, sep = "\n")
  quit(status = 1)
}
cat("NAMESPACE and man/ are clean\n")
'

# --- 3. structural baseline (fast, no brms) ----------------------------------
# The Phase 3 acceptance evidence: new_gp_k should be roughly constant across
# every scenario instead of tracking sqrt(n).
echo "==> dev/baseline-structural.R"
Rscript -e 'devtools::load_all(".", quiet = TRUE); source("dev/baseline-structural.R")'

# --- 4. tests ----------------------------------------------------------------
echo "==> devtools::test()"
Rscript -e 'devtools::test(stop_on_failure = TRUE)'

# --- 5. full check -----------------------------------------------------------
echo "==> devtools::check()"
Rscript -e 'devtools::check(document = FALSE, cran = TRUE, error_on = "warning")'

echo
echo "All stages passed."
echo "Remaining (slow, needs brms + a Stan toolchain):"
echo "  Rscript -e 'devtools::load_all(\".\"); source(\"dev/baseline-accuracy.R\")'"
