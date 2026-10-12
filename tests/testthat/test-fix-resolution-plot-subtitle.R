# tests/testthat/test-fix-resolution-plot-subtitle.R
# ---------------------------------------------------------------------------
# plot.resolution_profile() says why a profile is not supported, and labels
# the Moran panel as the diagnostic it is.  Goes with the two-line change to
# R/plotting-diagnostics.R (lines 688 and 735).
# ---------------------------------------------------------------------------

test_that("the subtitle of an unsupported profile gives the reason its bounds record", {
  skip_if_not_installed("ggplot2")
  set.seed(1)
  n <- 200
  pts <- sf::st_as_sf(data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), z = rnorm(n)),
                      coords = c("x", "y"), crs = 32632)
  # Fewer than two cells at min_cell_n: the subtitle read "(floor 2, ceiling
  # 2); the floor exceeds the ceiling".
  few <- suppressWarnings(resolution_profile(pts, "z", min_cell_n = 150,
                                             levels = c(2, 4, 8, 12, 16), nstart = 2))
  expect_identical(attr(few, "bounds")$ceiling_from, "fewer than two supported cells")
  sub <- plot(few)$labels$subtitle
  expect_match(sub, "; fewer than two cells meet min_cell_n = 150$")
  expect_false(grepl("the floor exceeds the ceiling", sub))
  # A floor above the ceiling is still called that.
  vm <- data.frame(model = c("Nug", "Exp"), psill = c(0.3, 1), range = c(0, 100 / 3),
                   stringsAsFactors = FALSE)
  sac <- structure(100, class = c("sac_range", "numeric"), variogram_model = vm,
                   crs = sf::st_crs(pts))
  hi <- suppressWarnings(resolution_profile(pts, "z", sac = sac, n_levels = 5, nstart = 2))
  expect_false(attr(hi, "bounds")$supported)
  expect_match(plot(hi)$labels$subtitle, "; the floor exceeds the ceiling$")
  # A supported profile has no such clause.
  ok <- suppressWarnings(resolution_profile(pts, "z", sac = 300, n_levels = 5, nstart = 2))
  expect_match(plot(ok)$labels$subtitle, "\\)$")
})

test_that("the Moran panel is labelled as a diagnostic", {
  skip_if_not_installed("ggplot2")
  set.seed(2)
  n <- 300
  pts <- sf::st_as_sf(data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), z = rnorm(n)),
                      coords = c("x", "y"), crs = 32632)
  p <- suppressWarnings(resolution_profile(pts, "z", sac = 300, n_levels = 6, nstart = 2))
  g <- plot(p, criteria = c("cp", "moran_z"))
  panels <- as.character(ggplot2::ggplot_build(g)$layout$layout$panel)
  expect_true("|Moran's z| of cell means (a diagnostic, not a count rule)" %in% panels)
  expect_false(any(grepl("Moran.*lower is better", panels)))
})
