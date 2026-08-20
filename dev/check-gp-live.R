# ---------------------------------------------------------------------------
# Is the GP basis fix actually live in the FIT path?  (~2 minutes)
#
# dev/baseline-structural.R calls .gp_basis_spec() directly, so it proves the
# helper computes the new rule -- but not that fit_bayesian_spatial_model()
# consumes it.  That distinction is not academic: a five-hour accuracy run once
# reported gp_k = 44 at n = 2000 (exactly floor(sqrt(n)), the rule this release
# replaced) because the script had loaded the INSTALLED package via library().
# Reading the source and seeing `if (gp_k_auto) gp_k <- gp_spec$k` is not the
# same as observing it.
#
# Run this BEFORE committing hours to dev/baseline-accuracy.R.
#
# Usage:  Rscript dev/check-gp-live.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(sf))

.sk_from_source <- FALSE
if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
  .sk_from_source <- TRUE
} else {
  library(spatialkit)
}
.sk_path <- tryCatch(getNamespaceInfo(asNamespace("spatialkit"), "path"),
                     error = function(e) NA_character_)

cat("spatialkit loaded from      : ", .sk_path, "\n", sep = "")
cat("working tree (not installed): ", .sk_from_source, "\n\n", sep = "")
if (!.sk_from_source)
  stop("This script is pointless against the installed package. Run it from ",
       "the package root with pkgload available.", call. = FALSE)

if (!requireNamespace("brms", quietly = TRUE))
  stop("brms is required.", call. = FALSE)

n <- 300L
set.seed(1)
x <- runif(n, 0, 10000); y <- runif(n, 0, 10000); w <- rnorm(n)
z <- 2 * sin(x / 2500) * cos(y / 2500) + 0.5 * w + rnorm(n, 0, 0.3)
pts <- st_as_sf(data.frame(x = x, y = y, z = z, w = w),
                coords = c("x", "y"), crs = 3857)

cat("Fitting one small model (n = 300, 1 chain, 400 iterations)...\n\n")
fit <- fit_bayesian_spatial_model(pts, "z", "w", chains = 1, iter = 400)

old_k <- max(15L, as.integer(floor(sqrt(n))))
new_k <- as.integer(fit$info$gp_k)
new_c <- as.numeric(fit$info$gp_c)

cat("\n--- what the FIT returned ---\n")
cat(sprintf("  gp_k        = %d\n", new_k))
cat(sprintf("  gp_c        = %.2f\n", new_c))
cat(sprintf("  gp_n_basis  = %s\n", format(fit$info$gp_n_basis)))
cat(sprintf("  old rule    = %d   (max(15, floor(sqrt(%d))))\n", old_k, n))
cat(sprintf("  old gp_c    = 1.50 (hard-coded)\n"))

ok <- TRUE
if (identical(new_k, old_k)) {
  ok <- FALSE
  cat("\nFAIL: gp_k equals the OLD rule exactly. The fix is not in this code path.\n")
} else if (new_k < 15L || new_k > 40L) {
  ok <- FALSE
  cat(sprintf("\nFAIL: gp_k = %d is outside the expected 15-40 band.\n", new_k))
}
if (!is.finite(new_c) || new_c < 1.2 || new_c > 8) {
  ok <- FALSE
  cat(sprintf("\nFAIL: gp_c = %.2f is outside the expected range.\n", new_c))
}

# The scale = FALSE / inverse-gamma prior fix exists to stop Stan rejecting
# initial values.  A run littered with "Gradient evaluated at the initial value
# is not finite" is the old prior, whatever gp_k says.
cat("\nAlso check the sampler output above: repeated\n",
    "  'Gradient evaluated at the initial value is not finite'\n",
    "means the OLD length-scale prior is in force. A clean fit shows none.\n",
    sep = "")

if (ok) {
  cat("\nPASS: the fit path is using the derived basis rule.\n",
      "dev/baseline-accuracy.R is now worth its runtime.\n", sep = "")
} else {
  cat("\nDo NOT run dev/baseline-accuracy.R until this passes.\n")
}
