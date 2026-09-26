# tests/testthat/test-review2-plotting.R
# ---------------------------------------------------------------------------
# Second review pass over the plotting functions: plots that failed only when
# printed, labels that said something the data did not, and a pooled line
# drawn in the wrong panel.  Each test builds the plot with
# ggplot2::ggplot_build() and reads what is drawn, since a ggplot object is
# built lazily and constructing one proves nothing.
# ---------------------------------------------------------------------------

layer_geoms_r2 <- function(p) vapply(p$layers, function(l) class(l$geom)[1L], character(1))

r2_points <- function(n = 60, crs = 32632, seed = 1) {
  set.seed(seed)
  sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
                          z = rnorm(n)),
               coords = c("x", "y"), crs = crs)
}


# ---- plot_tessellation_map() ---------------------------------------------

test_that("labels = TRUE labels the package's own cells without naming a column", {
  skip_if_not_installed("ggplot2")
  # The default label_col was "grid_id", which no function produces, so
  # labels = TRUE drew nothing on Voronoi cells, grids or summarize_by_cell()
  # output.
  pts <- r2_points()
  tess <- build_tessellation(pts, method = "voronoi", quiet = TRUE)
  p <- plot_tessellation_map(tess$cells, labels = TRUE)
  expect_true("GeomText" %in% layer_geoms_r2(p))
  txt <- ggplot2::layer_data(p, which(layer_geoms_r2(p) == "GeomText"))
  expect_equal(sort(as.numeric(txt$label)), sort(tess$cells$cell_id))

  cells <- summarize_by_cell(assign_features_to_polygons(pts, tess$cells),
                             response_var = "z", cells_sf = tess$cells)
  expect_false("cell_id" %in% names(cells))
  p2 <- plot_tessellation_map(cells, labels = TRUE)
  txt2 <- ggplot2::layer_data(p2, which(layer_geoms_r2(p2) == "GeomText"))
  expect_equal(sort(as.numeric(txt2$label)), sort(cells$poly_id))

  # A grid_id column still wins, and an explicit label_col is honoured.
  cells$grid_id <- paste0("g", cells$poly_id)
  p3 <- plot_tessellation_map(cells, labels = TRUE)
  expect_true(all(grepl("^g", ggplot2::layer_data(p3, 2L)$label)))
  p4 <- plot_tessellation_map(cells, labels = TRUE, label_col = "n")
  expect_equal(sort(ggplot2::layer_data(p4, 2L)$label), sort(cells$n))
})

test_that("labels = TRUE on a layer with no ID column draws no labels and logs why", {
  skip_if_not_installed("ggplot2")
  g <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(
    c(xmin = 0, ymin = 0, xmax = 1000, ymax = 1000), crs = sf::st_crs(32632))), n = c(2, 2))
  cells <- sf::st_sf(v = 1:4, geometry = g)
  lines <- capture_spatialkit_log(p <- plot_tessellation_map(cells, labels = TRUE))
  expect_false("GeomText" %in% layer_geoms_r2(p))
  expect_true(log_has(lines, "no ID column"))
})

test_that("a units, Date, POSIXct or difftime fill column draws instead of failing at print", {
  skip_if_not_installed("ggplot2")
  # The scale was chosen with is.numeric(): Date, POSIXct and difftime got a
  # discrete scale ("Continuous value supplied to a discrete scale") and an
  # st_area() column a continuous one whose arithmetic failed on units --
  # both only once the plot was built.
  g <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(
    c(xmin = 0, ymin = 0, xmax = 1000, ymax = 1000), crs = sf::st_crs(32632))), n = c(3, 3))
  cells <- sf::st_sf(id = seq_along(g), geometry = g)
  cells$area  <- sf::st_area(cells) * seq(0.5, 1.5, length.out = 9)
  cells$when  <- as.Date("2024-01-01") + 0:8
  cells$stamp <- as.POSIXct("2024-01-01", tz = "UTC") + 3600 * (0:8)
  cells$lag   <- as.difftime(1:9, units = "days")
  for (col in c("area", "when", "stamp", "lag")) {
    p <- plot_tessellation_map(cells, fill_col = col)
    expect_no_error(b <- ggplot2::ggplot_build(p))
    fills <- b$data[[1]]$fill
    # A continuous scale: nine values, nine different colours, none NA.
    expect_equal(length(unique(fills)), 9L, label = col)
    expect_false(anyNA(fills), label = col)
  }
  # The unit goes in the legend title; an explicit title is left alone.
  expect_identical(plot_tessellation_map(cells, fill_col = "area")$scales$get_scales("fill")$name,
                   "area [m^2]")
  expect_identical(plot_tessellation_map(cells, fill_col = "lag")$scales$get_scales("fill")$name,
                   "lag [days]")
  expect_identical(plot_tessellation_map(cells, fill_col = "area", legend_title = "A")$scales$get_scales("fill")$name,
                   "A")
  # A Date legend is labelled with dates ("Jan 01" in English), not with the
  # day counts (19723) a numeric scale would print.
  b <- ggplot2::ggplot_build(plot_tessellation_map(cells, fill_col = "when"))
  labs <- b$plot$scales$get_scales("fill")$get_labels()
  expect_gt(length(labs), 1L)
  expect_false(any(grepl("^[0-9.e+]+$", labs)))
})
