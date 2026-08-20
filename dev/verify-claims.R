# dev/verify-claims.R
# ---------------------------------------------------------------------------
# Verifies two claims that the test suite and the accuracy baseline do not
# cover.  Neither needs a long MCMC run.
#
#   1. NEWS.md asserts fits are "substantially faster at n >= 1,000".  The
#      basis count dropped from 1936 to 576 at n = 2000, but the two were
#      never timed against each other.  Fits one model at each basis count
#      and reports the ratio.
#
#   2. The wide-extent CRS change is covered only by unit tests on synthetic
#      bounding boxes.  This runs the actual downstream consumers -- variogram
#      range and spatial block sizing -- on continental-extent data under a
#      forced single UTM zone versus the automatic equal-area choice, which is
#      the behaviour the fix exists to produce.
#
# Usage:  Rscript -e 'devtools::load_all("."); source("dev/verify-claims.R")'
# Expect roughly 12-20 minutes, almost all of it in part 1.
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


# =========================================================================
# 1. Does the smaller basis actually cost less?
# =========================================================================
if (!requireNamespace("brms", quietly = TRUE)) {
  message("brms not installed; skipping part 1.")
} else {
  set.seed(1); n <- 2000
  x <- runif(n, 0, 10000); y <- runif(n, 0, 10000); w <- rnorm(n)
  z <- 2 * sin(x / 2500) * cos(y / 2500) + 0.5 * w + rnorm(n, 0, 0.3)
  pts <- st_as_sf(data.frame(x = x, y = y, z = z, w = w),
                  coords = c("x", "y"), crs = 3857)

  # The old rule at n = 2000: max(15, floor(sqrt(2000))) = 44 -> 1936 basis.
  # gp_c = 1.5 was the old hard-coded default.
  t_old <- system.time(
    fit_old <- fit_bayesian_spatial_model(pts, "z", "w", gp_k = 44, gp_c = 1.5,
                                          chains = 2, iter = 1000,
                                          compute_loo = FALSE)
  )[["elapsed"]]

  t_new <- system.time(
    fit_new <- fit_bayesian_spatial_model(pts, "z", "w",
                                          chains = 2, iter = 1000,
                                          compute_loo = FALSE)
  )[["elapsed"]]

  cat("\n=========== 1. basis count vs wall clock (n = 2000) ===========\n")
  cat(sprintf("  old rule : gp_k = %-3d  basis = %-5d  %7.1f s\n",
              fit_old$info$gp_k, fit_old$info$gp_n_basis, t_old))
  cat(sprintf("  derived  : gp_k = %-3d  basis = %-5d  %7.1f s\n",
              fit_new$info$gp_k, fit_new$info$gp_n_basis, t_new))
  cat(sprintf("  speedup  : %.2fx  (basis ratio %.2fx)\n",
              t_old / t_new, fit_old$info$gp_n_basis / fit_new$info$gp_n_basis))
  cat("\n  -> If the speedup is well under 1.5x, soften the NEWS.md wording.\n")
}

# =========================================================================
# 2. Does the equal-area switch actually change the downstream numbers?
# =========================================================================
set.seed(2); m <- 600
lon <- runif(m, -124, -67); lat <- runif(m, 25, 49)          # CONUS extent
val <- sin(lon / 8) * cos(lat / 6) + rnorm(m, 0, 0.2)
conus <- st_as_sf(data.frame(lon = lon, lat = lat, val = val),
                  coords = c("lon", "lat"), crs = 4326)

auto <- ensure_projected(conus)                               # equal-area
utm  <- ensure_projected(conus, target_crs = st_crs(32615))    # forced 1 zone

bbw <- function(x) { b <- st_bbox(x); as.numeric(b["xmax"] - b["xmin"]) }

cat("\n=========== 2. wide extent: forced UTM vs equal-area ===========\n")
cat(sprintf("  auto CRS          : %s\n",
            substr(st_crs(auto)$proj4string, 1, 60)))
cat(sprintf("  bbox width  auto  : %10.0f m\n", bbw(auto)))
cat(sprintf("  bbox width  UTM15 : %10.0f m  (%+.1f%%)\n",
            bbw(utm), 100 * (bbw(utm) - bbw(auto)) / bbw(auto)))

if (requireNamespace("gstat", quietly = TRUE)) {
  r_auto <- estimate_sac_range(auto, "val", seed = 1)
  r_utm  <- estimate_sac_range(utm,  "val", seed = 1)
  cat(sprintf("  SAC range   auto  : %10.0f m\n", r_auto))
  cat(sprintf("  SAC range   UTM15 : %10.0f m  (%+.1f%%)\n",
              r_utm, 100 * (r_utm - r_auto) / r_auto))
}

f_auto <- make_folds(auto, k = 4, method = "block_kfold", seed = 1)
f_utm  <- make_folds(utm,  k = 4, method = "block_kfold", seed = 1)
cat(sprintf("  block grid  auto  : %d x %d\n",
            f_auto$params$grid_nx, f_auto$params$grid_ny))
cat(sprintf("  block grid  UTM15 : %d x %d\n",
            f_utm$params$grid_nx,  f_utm$params$grid_ny))
cat("\n  -> A multi-percent gap here is the error the fix removes.\n\n")
