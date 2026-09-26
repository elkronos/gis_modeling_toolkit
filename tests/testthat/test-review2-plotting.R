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


# ---- plot_folds() --------------------------------------------------------

test_that("plot_folds draws the layer the folds came from when only the blocks carry a CRS", {
  skip_if_not_installed("ggplot2")
  # make_folds() projects CRS-less lon/lat points to a UTM zone, so its blocks
  # carry a CRS the points do not; coord_sf() then aborted at print with
  # "cannot transform sfc object with missing crs".
  set.seed(1); n <- 80
  ll <- sf::st_as_sf(data.frame(x = 10 + runif(n, 0, 0.01), y = 45 + runif(n, 0, 0.01)),
                     coords = c("x", "y"))
  f <- suppressWarnings(make_folds(ll, k = 5, method = "block_kfold", block_size = 300))
  expect_false(is.na(sf::st_crs(f$params$blocks)))
  expect_warning(p <- plot_folds(f, ll), "plot_folds\\(\\): `points_sf` has no CRS; its coordinates look like lon/lat")
  expect_no_error(b <- ggplot2::ggplot_build(p))
  # The points are reprojected onto the blocks, not stamped as metres.
  blk <- sf::st_bbox(b$data[[1]]$geometry); pt <- sf::st_bbox(b$data[[2]]$geometry)
  expect_true(pt[["xmin"]] >= blk[["xmin"]] && pt[["xmax"]] <= blk[["xmax"]] &&
              pt[["ymin"]] >= blk[["ymin"]] && pt[["ymax"]] <= blk[["ymax"]])
  expect_equal(nrow(b$data[[2]]), n)
})

test_that("plot_folds aligns a CRS-less boundary or CRS-less points to the other layers", {
  skip_if_not_installed("ggplot2")
  pts <- r2_points(n = 80)
  f <- make_folds(pts, k = 5, method = "block_kfold", block_size = 300, seed = 1)
  bnd <- sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(pts)))
  bnd_na <- sf::st_set_crs(bnd, NA)
  expect_warning(p <- plot_folds(f, pts, boundary = bnd_na), "`boundary` has no CRS")
  expect_no_error(b <- ggplot2::ggplot_build(p))
  expect_length(b$data, 3L)

  # CRS-less points whose folds took a boundary's CRS inside make_folds().
  pts_na <- sf::st_set_crs(pts, NA)
  f2 <- suppressWarnings(make_folds(pts_na, k = 5, method = "block_kfold",
                                    block_size = 300, boundary = bnd, seed = 1))
  expect_warning(p2 <- plot_folds(f2, pts_na), "`points_sf` has no CRS")
  expect_no_error(ggplot2::ggplot_build(p2))
  # CRS-less points and folds, a boundary with a CRS.
  f3 <- suppressWarnings(make_folds(pts_na, k = 5, method = "block_kfold",
                                    block_size = 300, seed = 1))
  p3 <- suppressWarnings(plot_folds(f3, pts_na, boundary = bnd))
  expect_no_error(ggplot2::ggplot_build(p3))
  # No layer with a CRS: drawn as they are, with nothing to warn about.
  expect_no_warning(p4 <- plot_folds(f3, pts_na))
  expect_no_error(ggplot2::ggplot_build(p4))
})

test_that("plot_folds' subtitle gives no unit for folds built without a CRS", {
  skip_if_not_installed("ggplot2")
  # make_folds() records NA_character_ as the CRS and nzchar(NA) is TRUE, so
  # the subtitle read "Block size 300 (NA units)".
  pts_na <- sf::st_set_crs(r2_points(n = 80), NA)
  fb <- suppressWarnings(make_folds(pts_na, k = 4, method = "block_kfold",
                                    block_size = 300, seed = 1))
  expect_true(is.na(fb$params$crs))
  sub_b <- plot_folds(fb, pts_na)$labels$subtitle
  expect_match(sub_b, "Block size 300\n")
  expect_false(grepl("NA units", sub_b))
  fl <- suppressWarnings(make_folds(pts_na, k = 80, method = "buffered_loo", buffer = 150))
  sub_l <- plot_folds(fl, pts_na)$labels$subtitle
  expect_identical(sub_l, "Leave-one-out with a 150 buffer")
})
