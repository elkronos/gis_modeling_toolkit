# dev/baseline-structural.R
# ---------------------------------------------------------------------------
# Structural baseline for the Phase 3 GP change.
#
# Records WHICH gp_k / gp_c get chosen and how many basis functions that
# implies, WITHOUT fitting anything.  brms builds a full tensor grid over the
# two GP covariates, so the model carries gp_k^2 basis functions -- see
# brms/R/data-predictor.R, data_gp():  out[[paste0("NBgp", pi)]] <- k ^ D
#
# Run before Phase 3, and again at Step 3.9.  `dev/` is in .Rbuildignore, so
# nothing here ships.
#
# Usage:  Rscript dev/baseline-structural.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(sf))

# ---------------------------------------------------------------------------
# Load the WORKING TREE, never whatever happens to be installed.
#
# library(spatialkit) loads the INSTALLED package.  devtools::test() loads from
# source, so a working tree can be many changes ahead of the installed copy and
# a dev script using library() will silently measure the wrong code -- which is
# exactly what happened once, costing a five-hour baseline run.
# ---------------------------------------------------------------------------
.sk_from_source <- FALSE
if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)   # internals needed (.gp_basis_spec)
  .sk_from_source <- TRUE
} else {
  library(spatialkit)
}

# The namespace PATH is the ground truth for which copy is loaded.
# packageVersion() is not: after load_all() it may still report the installed
# DESCRIPTION, which is precisely the kind of near-miss diagnostic that let the
# original mistake through.
.sk_path <- tryCatch(getNamespaceInfo(asNamespace("spatialkit"), "path"),
                     error = function(e) NA_character_)
message("spatialkit loaded from : ", .sk_path)
message("working tree (not installed): ", .sk_from_source)
if (!.sk_from_source)
  message("*** MEASURING THE INSTALLED PACKAGE -- your edits are NOT in this run. ***\n",
          "*** Run from the package root with pkgload installed.                 ***")


set.seed(42)

scenarios <- list(
  uniform_200   = function() matrix(runif(2 * 200),   ncol = 2),
  uniform_1000  = function() matrix(runif(2 * 1000),  ncol = 2),
  uniform_10000 = function() matrix(runif(2 * 10000), ncol = 2),
  clustered_20  = function() {
    ctr <- matrix(runif(2 * 20), ncol = 2)
    lab <- sample(20, 2000, replace = TRUE)
    ctr[lab, ] + matrix(rnorm(2 * 2000, 0, 0.02), ncol = 2)
  },
  elongated     = function() cbind(runif(2000, 0, 10), runif(2000, 0, 0.5))
)

# The heuristic in force before Phase 3 (R/model-bayesian.R:151-155).
current_gp_k <- function(n) {
  max(5L, as.integer(min(n / 3, max(15L, floor(sqrt(n))))))
}

# Both rules are computed here, side by side, so this script is a direct
# before/after comparison that can be run at any point -- the old rule is
# reimplemented inline above rather than read out of the package, so it stays
# available after the new derivation has landed.
baseline <- lapply(names(scenarios), function(nm) {
  xy <- scale(scenarios[[nm]]())        # mirrors the package's per-axis scaling
  n  <- nrow(xy)
  b  <- gp_lengthscale_bounds(xy)
  S  <- max(apply(xy, 2, function(z) diff(range(z)) / 2))

  k_old <- current_gp_k(n)
  spec  <- spatialkit:::.gp_basis_spec(xy, b)

  list(
    scenario     = nm,
    n            = n,
    # before: chosen from n, with gp_c hard-coded
    old_gp_k     = k_old,
    old_gp_c     = 1.5,
    old_n_basis  = k_old^2,
    # after: derived from the length-scale/domain ratio
    new_gp_k     = spec$k,
    new_gp_c     = round(spec$c, 3),
    new_n_basis  = spec$k^2,
    basis_ratio  = round(k_old^2 / spec$k^2, 2),
    ls_lower     = round(unname(b[["lower"]]), 4),
    ls_upper     = round(unname(b[["upper"]]), 4),
    S            = round(S, 4),
    # smallest length-scale each configuration can resolve
    old_ell_min  = round(1.75 * 1.5 * S / k_old, 4),
    new_ell_min  = round(1.75 * spec$c * S / spec$k, 4)
  )
})
names(baseline) <- names(scenarios)

df <- do.call(rbind, lapply(baseline, function(r) as.data.frame(r,
                            stringsAsFactors = FALSE)))
cols <- c("scenario", "n", "old_gp_k", "old_n_basis",
          "new_gp_k", "new_gp_c", "new_n_basis", "basis_ratio",
          "ls_lower", "new_ell_min")
print(df[, cols], row.names = FALSE)

attr(baseline, "captured_at") <- Sys.time()
attr(baseline, "package_version") <- as.character(utils::packageVersion("spatialkit"))
saveRDS(baseline, "dev/baseline-structural.rds")
cat("\nWrote dev/baseline-structural.rds\n")

# What to look for:
#
#   old_n_basis   225 / 961 / 10000 for the three uniform cases -- i.e. equal
#                 to n, because the old rule reduced to max(15, sqrt(n)).
#   new_gp_k      roughly constant (about 21-25) across EVERY scenario,
#                 including the clustered and elongated ones.  Constancy is the
#                 property being bought; the exact value is not.
#   new_gp_c      about 3.2-3.6, versus the hard-coded 1.5.
#   basis_ratio   ~0.5 at n=200 (the new basis count is LARGER -- expected, a
#                 correction rather than an optimisation), ~2 at n=1000, ~19 at
#                 n=10000.
#   new_ell_min   should sit at or below ls_lower.  If it is above, the basis
#                 cannot resolve the shortest length-scale the data supports and
#                 k_min needs raising.
