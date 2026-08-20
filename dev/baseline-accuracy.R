# dev/baseline-accuracy.R
# ---------------------------------------------------------------------------
# Accuracy baseline for the Phase 3 GP change.
#
# The structural baseline proves the COST changed.  Only fitted models prove
# ACCURACY did not degrade, which is the actual acceptance criterion at
# Step 3.9.
#
# Deliberately cheap so it can be re-run:
#   * n = 300 and n = 2000 only.  Full MCMC at n = 10000 is hours, and the
#     effect is already unambiguous at 2000 (current gp_k = 44 -> 1936 basis
#     functions, versus ~23 -> ~529 after the fix).
#   * chains = 2, iter = 1000.  This is a RELATIVE comparison between two
#     model specifications, not a publication fit.  The settings are stored in
#     the .rds so Step 3.9 re-runs identically.
#   * parallel = FALSE throughout.  Phase 2 must land before the parallel path
#     is trustworthy as a measuring device.
#
# cv_gwr() is the CONTROL: Phase 3 must not change it at all.
#
# Usage:  Rscript dev/baseline-accuracy.R        (expect this to take a while)
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
  pkgload::load_all(".", quiet = TRUE)
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


if (!requireNamespace("brms", quietly = TRUE))
  stop("brms is required for the accuracy baseline.")

SETTINGS <- list(chains = 2L, iter = 1000L, k_folds = 4L, seed = 123L,
                 sizes = c(300L, 2000L))

# --- Simulated field with genuine spatial structure ------------------------
make_field <- function(n, seed = 1) {
  set.seed(seed)
  x <- runif(n, 0, 10000); y <- runif(n, 0, 10000)
  w <- rnorm(n)
  # smooth spatial signal + covariate effect + noise
  z <- 2 * sin(x / 2500) * cos(y / 2500) + 0.5 * w + rnorm(n, 0, 0.3)
  st_as_sf(data.frame(x = x, y = y, z = z, w = w),
           coords = c("x", "y"), crs = 3857)
}

results <- list()

for (n in SETTINGS$sizes) {
  message("=== n = ", n, " ===")
  pts <- make_field(n)

  t_bayes <- system.time({
    cvb <- cv_bayes(pts, "z", "w", k = SETTINGS$k_folds, seed = SETTINGS$seed,
                    fit_args = list(chains = SETTINGS$chains,
                                    iter   = SETTINGS$iter),
                    parallel = FALSE)
  })

  fit <- fit_bayesian_spatial_model(
    pts, "z", "w", chains = SETTINGS$chains, iter = SETTINGS$iter
  )

  t_gwr <- system.time({
    cvg <- if (requireNamespace("GWmodel", quietly = TRUE)) {
      cv_gwr(pts, "z", "w", k = SETTINGS$k_folds, seed = SETTINGS$seed,
             parallel = FALSE)
    } else NULL
  })

  results[[as.character(n)]] <- list(
    n              = n,
    bayes_overall  = cvb$overall,
    bayes_folds    = cvb$fold_metrics,
    bayes_elapsed  = unname(t_bayes[["elapsed"]]),
    gp_k           = fit$info$gp_k,
    gp_c           = if (is.null(fit$info$gp_c)) NA_real_ else fit$info$gp_c,
    n_basis        = if (is.null(fit$info$gp_n_basis)) fit$info$gp_k^2
                     else fit$info$gp_n_basis,
    ell_min        = if (is.null(fit$info$gp_ell_min)) NA_real_
                     else fit$info$gp_ell_min,
    looic          = fit$info$looic,
    # CONTROL -- must be identical after Phase 3
    gwr_overall    = if (!is.null(cvg)) cvg$overall else NULL,
    gwr_elapsed    = unname(t_gwr[["elapsed"]])
  )
}

# --- Phase 4 control: narrow-extent SAC range and block sizing --------------
set.seed(7)
narrow <- st_as_sf(
  data.frame(lon = runif(400, -78.9, -78.3), lat = runif(400, 35.5, 36.1),
             z = rnorm(400)),
  coords = c("lon", "lat"), crs = 4326
)
results$phase4_control <- list(
  sac_range  = estimate_sac_range(narrow, "z", seed = 1),
  crs_epsg   = st_crs(ensure_projected(narrow))$epsg,
  fold_params = make_folds(narrow, k = 4, method = "block_kfold",
                           seed = 1)$params
)

attr(results, "settings")        <- SETTINGS
attr(results, "captured_at")     <- Sys.time()
attr(results, "package_version") <- as.character(utils::packageVersion("spatialkit"))
attr(results, "session")         <- utils::sessionInfo()
saveRDS(results, "dev/baseline-accuracy.rds")

cat("\nWrote dev/baseline-accuracy.rds\n")
for (n in as.character(SETTINGS$sizes)) {
  r <- results[[n]]
  cat(sprintf("n=%-5s gp_k=%-4d basis=%-6d RMSE=%.4f  R2=%.4f  %.1fs\n",
              n, r$gp_k, r$n_basis, r$bayes_overall$RMSE,
              r$bayes_overall$R2, r$bayes_elapsed))
}

# ---------------------------------------------------------------------------
# Check the acceptance criteria rather than leaving them in a comment.
#
# These numbers look perfectly plausible when the WRONG code is measured -- a
# five-hour run against the installed package once produced gp_k = 44 at
# n = 2000 and nobody noticed until the values were read closely.  The
# distinguishing signature is gp_k tracking floor(sqrt(n)), which is precisely
# the rule this release replaced.
# ---------------------------------------------------------------------------
cat("\n--- acceptance criteria ---\n")
ok <- TRUE
for (n in SETTINGS$sizes) {
  r     <- results[[as.character(n)]]
  old_k <- max(15L, as.integer(floor(sqrt(n))))
  if (identical(as.integer(r$gp_k), old_k)) {
    ok <- FALSE
    cat(sprintf("FAIL n=%d: gp_k = %d is exactly floor(sqrt(n)) -- the OLD rule.\n",
                n, r$gp_k))
  } else if (r$gp_k < 15L || r$gp_k > 40L) {
    ok <- FALSE
    cat(sprintf("FAIL n=%d: gp_k = %d is outside the expected 15-40 band.\n",
                n, r$gp_k))
  } else {
    cat(sprintf("ok   n=%d: gp_k = %d (old rule would give %d)\n",
                n, r$gp_k, old_k))
  }
  if (!is.na(r$gp_c) && (r$gp_c < 1.2 || r$gp_c > 8)) {
    ok <- FALSE
    cat(sprintf("FAIL n=%d: gp_c = %.2f is outside the expected range.\n",
                n, r$gp_c))
  }
}
if (!ok) {
  cat("\nThe GP basis rule in effect is NOT the one this release introduces.\n",
      "Almost always this means the script measured the INSTALLED package\n",
      "instead of the working tree. Version loaded: ",
      as.character(utils::packageVersion("spatialkit")), "\n", sep = "")
} else {
  cat("\nBasis rule looks right. Still check by eye:\n",
      "  RMSE / R2   not materially worse   <- the actual acceptance criterion\n",
      "  elapsed     down at n = 2000\n",
      "  gwr_overall identical              <- if this moved, something leaked\n",
      sep = "")
}
