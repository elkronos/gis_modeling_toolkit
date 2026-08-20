# dev/verify-crs-distortion.R
# ---------------------------------------------------------------------------
# Measures projection distance error against GEODESIC truth.
#
# An earlier version of this check compared bounding-box widths under Albers
# versus a forced UTM zone.  That was the wrong metric: Albers preserves area,
# not distance, so the comparison measured the gap between two distorted
# projections rather than the error of either one.  The honest reference is the
# true ellipsoidal distance, which sf computes via s2 on lon/lat geometry.
#
# Note on the error floor: st_distance() on lon/lat geometry uses s2, which
# models the Earth as a sphere, while UTM and Albers are built on the WGS84
# ellipsoid (flattening ~0.34%).  That leaves a systematic floor of roughly
# 0.15% in EVERY row -- visible in the narrow Wake County case, where true UTM
# scale distortion is only about 0.02%.  The floor affects both columns
# equally, so the auto-vs-forced comparison is valid; the absolute figures are
# not a pure measure of projection error.
#
# Runs in seconds.  No MCMC.
#
# Usage:  Rscript -e 'devtools::load_all("."); source("dev/verify-crs-distortion.R")'
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


report <- function(label, lon_rng, lat_rng, n = 200, seed = 42) {
  set.seed(seed)
  pts_ll <- st_as_sf(
    data.frame(lon = runif(n, lon_rng[1], lon_rng[2]),
               lat = runif(n, lat_rng[1], lat_rng[2])),
    coords = c("lon", "lat"), crs = 4326
  )

  # Truth: ellipsoidal distance on the sphere.
  d_true <- as.numeric(st_distance(pts_ll))

  # What the package picks automatically.
  auto   <- suppressWarnings(ensure_projected(pts_ll))
  d_auto <- as.numeric(st_distance(auto))

  # A single UTM zone chosen from the centroid -- the pre-1.1.0 behaviour.
  ctr  <- st_coordinates(st_centroid(st_union(st_geometry(pts_ll))))
  zone <- max(1, min(60, floor((ctr[1] + 180) / 6) + 1))
  epsg <- if (ctr[2] < 0) 32700 + zone else 32600 + zone
  utm   <- st_transform(pts_ll, epsg)
  d_utm <- as.numeric(st_distance(utm))

  keep <- d_true > 0
  err  <- function(d) 100 * abs(d[keep] - d_true[keep]) / d_true[keep]

  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("  auto CRS            : %s\n",
              if (!is.na(st_crs(auto)$epsg)) paste0("EPSG:", st_crs(auto)$epsg)
              else sub("^\\+proj=([a-z]+).*", "\\1", st_crs(auto)$proj4string)))
  cat(sprintf("  forced single zone  : EPSG:%d\n", epsg))
  cat(sprintf("  %-22s median %6.3f%%   90th pct %6.3f%%   max %6.3f%%\n",
              "distance error, auto  :", median(err(d_auto)),
              quantile(err(d_auto), 0.9), max(err(d_auto))))
  cat(sprintf("  %-22s median %6.3f%%   90th pct %6.3f%%   max %6.3f%%\n",
              "distance error, UTM   :", median(err(d_utm)),
              quantile(err(d_utm), 0.9), max(err(d_utm))))
  invisible(list(auto = median(err(d_auto)), utm = median(err(d_utm))))
}

cat("\n===== projection distance error vs geodesic truth =====\n")
cat("(the package switches to equal-area beyond 5 deg from the zone meridian)\n")

report("Wake County, NC  (narrow -- must stay UTM)", c(-78.9, -78.3), c(35.5, 36.1))
report("Texas            (~13 deg wide)",           c(-106, -93),    c(26, 36))
report("CONUS            (~57 deg wide)",           c(-124, -67),    c(25, 49))
report("Europe           (~50 deg wide)",           c(-10, 40),      c(36, 60))

cat("\n  -> For narrow extents both should be tiny and near-identical.\n")
cat("  -> For wide extents the forced single zone should be materially worse.\n\n")
