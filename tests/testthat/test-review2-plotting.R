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


# ---- plot.spatial_fit() ---------------------------------------------------

test_that("plot.spatial_fit draws a custom fit that has no residuals() method", {
  skip_if_not_installed("ggplot2")
  # ?new_spatial_fit calls residuals() methods optional; without one,
  # residuals.default() returned NULL and every type but "coefficients"
  # stopped with "could not extract residuals".
  pts <- surf_test_points(n = 70)
  engine <- stats::lm(z ~ w, sf::st_drop_geometry(pts))
  # With neither method, the error names the one a custom subclass must have.
  bare <- new_spatial_fit("r2nomethods_fit", engine, z ~ w, "z", "w", pts)
  expect_error(plot(bare, type = "residuals"), "Define a fitted.r2nomethods_fit\\(\\) method")
  fit <- new_spatial_fit("r2noresid_fit", engine, z ~ w, "z", "w", pts)
  registerS3method("fitted", "r2noresid_fit",
                   function(object, ...) as.numeric(stats::fitted(object$engine)))
  expect_null(stats::residuals(fit))
  p <- plot(fit, type = "residuals")
  expect_no_error(b <- ggplot2::ggplot_build(p))
  expect_equal(nrow(b$data[[1]]), 70L)
  expect_equal(p$data$.resid, as.numeric(stats::residuals(engine)), tolerance = 1e-10)
  po <- plot(fit, type = "observed_predicted")
  bo <- ggplot2::ggplot_build(po)
  expect_equal(sort(bo$data[[2]]$x), sort(pts$z), tolerance = 1e-8)
})

test_that("the residual and response variograms are compared only over the same pairs", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gstat")
  # estimate_sac_range() returns the widest single direction when the
  # all-pairs fit is unusable -- as it often is for a response with a trend
  # whose residuals are fine -- and that overlay was labelled plainly as the
  # response, its sill set against the all-pairs residual curve's.
  set.seed(5); n <- 200
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- as.matrix(stats::dist(cbind(x, y)))
  z <- as.numeric(t(chol(exp(-d / 50) + diag(1e-4, n))) %*% rnorm(n))
  pts <- sf::st_as_sf(data.frame(x = x, y = y, z = z), coords = c("x", "y"), crs = 3857)
  sac <- estimate_sac_range(pts, "z")
  expect_true(is.finite(sac))
  expect_false(isTRUE(attr(sac, "anisotropy_used")))
  one_dir <- function(s, widest) {
    attr(s, "anisotropy_used") <- TRUE
    attr(s, "directional") <- c(`0` = 40, `45` = NA, `90` = 40, `135` = 40)
    attr(s, "directional")[widest] <- 60
    s
  }
  deg <- "\u00b0 \u00b1 22.5\u00b0"
  draw <- function(main, ov) .draw_sac_variogram(main, what = "Residual variogram",
                                                  overlay = ov, overlay_label = "Response (z)")
  cap <- function(p) gsub("\n", " ", p$labels$caption)

  # All-pairs residuals, one-direction response.
  p1 <- draw(sac, one_dir(sac, "90"))
  expect_no_error(ggplot2::ggplot_build(p1))
  expect_match(cap(p1), paste0("Response \\(z\\), the 90", deg, " direction only"))
  expect_match(cap(p1), "Sills not compared: the two variograms are not over the same point pairs")
  expect_false(grepl("Residual sill is", cap(p1)))
  # One-direction residuals, all-pairs response.
  p2 <- draw(one_dir(sac, "0"), sac)
  expect_match(p2$labels$title, paste0("Residual variogram, 0", deg))
  expect_match(cap(p2), "Sills not compared: the two variograms are not over the same point pairs")
  # Different directions are not the same pairs either; the same direction is.
  expect_match(cap(draw(one_dir(sac, "0"), one_dir(sac, "90"))), "Sills not compared")
  expect_match(cap(draw(one_dir(sac, "90"), one_dir(sac, "90"))), "Residual sill is [0-9]+% of the response sill")
  # Both all-pairs: compared, as before.
  expect_match(cap(draw(sac, sac)), "Residual sill is 100% of the response sill")
})

test_that("the no-sill subtitle does not claim a range range_frac refused ran past the lags", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gstat")
  # range_frac < 1 refuses a range below the largest lag fitted; the subtitle
  # printed that lag as the bound and said the variogram never reached a sill.
  set.seed(3); n <- 250
  xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
  D  <- as.matrix(stats::dist(xy))
  xy$z <- as.numeric(t(chol(exp(-D / 150) + diag(0.5, n))) %*% rnorm(n))
  pts <- sf::st_as_sf(xy, coords = c("x", "y"), crs = 32632)
  r <- suppressWarnings(estimate_sac_range(pts, "z", range_frac = 0.1))
  expect_identical(attr(r, "rejected_reason"), "fitted range exceeds the largest lag fitted")
  expect_true(attr(r, "rejected_range") <= attr(r, "cutoff_dist"))
  sub <- gsub("\n", " ", plot(r)$labels$subtitle)
  expect_false(grepl("never reached a sill", sub))
  expect_match(sub, "within the largest lag fitted .* `range_frac` accepts")
  # A range past the largest lag keeps the sill wording.
  r_over <- r
  attr(r_over, "rejected_range") <- 2 * attr(r, "cutoff_dist")
  expect_match(gsub("\n", " ", plot(r_over)$labels$subtitle),
               "exceeds the largest lag fitted .* never reached a sill")
})


# ---- plot(type = "coefficients") on a GWR fit ------------------------------

test_that("the GWR coefficient subtitle is right with mask = FALSE", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sp")
  # A clean survey is not a missing one.
  set.seed(3); n <- 200
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000); a <- rnorm(n); b <- rnorm(n)
  pts <- sf::st_as_sf(data.frame(x = x, y = y, a = a, b = b,
                                 z = (1 + 2 * x / 1000) * a - b + rnorm(n, 0, 0.3)),
                      coords = c("x", "y"), crs = 32632)
  fit <- suppressWarnings(fit_gwr_model(pts, "z", c("a", "b"), adaptive = TRUE, bandwidth = 60))
  expect_true(is.data.frame(fit$info$local_collinearity))
  expect_identical(fit$info$n_local_collinear, 0L)
  s0 <- plot(fit, type = "coefficients", mask = FALSE)$labels$subtitle
  expect_false(grepl("not surveyed", s0))
  expect_match(s0, "every local design is well conditioned")

  # Collinear windows and a non-finite coefficient: mask = FALSE still says
  # the collinear windows are drawn as if reliable.
  set.seed(5); cl <- rep(1:4, each = 50)
  d <- sf::st_as_sf(data.frame(
    x = c(runif(50, 0, 100), runif(50, 400, 500), runif(50, 0, 100), runif(50, 400, 500)),
    y = c(runif(50, 0, 100), runif(50, 0, 100), runif(50, 400, 500), runif(50, 400, 500)),
    a = rnorm(200), soil = c(0.5, 1.5, 0.5, 1.5)[cl] + rnorm(200, 0, 0.01)),
    coords = c("x", "y"), crs = 32632)
  d$z <- 2 * d$a + 3 * d$soil + rnorm(200, 0, 0.2)
  f70 <- suppressWarnings(fit_gwr_model(d, "z", c("a", "soil"), adaptive = TRUE, bandwidth = 70))
  f70$engine$SDF@data$a[1:3] <- NaN
  cn_bad <- with(f70$info$local_collinearity, !is.finite(cn) | cn > 30)
  n_cn <- sum(cn_bad[-(1:3)])
  expect_gt(n_cn, 0L)
  s1 <- plot(f70, type = "coefficients", term = "a", mask = FALSE)$labels$subtitle
  expect_match(s1, "^3 of 200 locations masked")
  expect_match(s1, "3 with a non-finite coefficient")
  expect_match(s1, sprintf("mask = FALSE: %d location\\(s\\) with a collinear local design are drawn as if reliable", n_cn))
  # mask = TRUE is unchanged.
  s2 <- plot(f70, type = "coefficients", term = "a")$labels$subtitle
  expect_match(s2, sprintf("^%d of 200 locations masked", n_cn + 3L))
  expect_false(grepl("mask = FALSE", s2))
})
