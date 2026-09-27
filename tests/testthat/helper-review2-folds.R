# tests/testthat/helper-review2-folds.R
# ---------------------------------------------------------------------------
# Shared by the second-review tests of make_folds(), fold_separation() and
# cv_block_size_sweep().
# ---------------------------------------------------------------------------

r2_pts <- function(x, y, crs = 32632, ...) {
  sf::st_as_sf(data.frame(x = x, y = y, ...), coords = c("x", "y"), crs = crs)
}
r2_fit <- function(train_sf) lm_spatial_fit(train_sf, "z", "a")
r2_quiet <- function(expr) {
  # Keep the console clear of the package's own WARN lines; R conditions
  # still reach the expectations.
  logger::with_log_threshold(expr, threshold = logger::FATAL,
                             namespace = "spatialkit", index = 2)
}
