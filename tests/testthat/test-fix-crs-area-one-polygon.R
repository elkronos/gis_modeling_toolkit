# tests/testthat/test-fix-crs-area-one-polygon.R
# ---------------------------------------------------------------------------
# One polygon gives one ratio of planar to geodesic area, which has no
# spread, so .crs_area_error() gives NA for it, and a polygon s2 refuses in
# degrees (edges that cross, a spike) gives NA too.
# ensure_projected(purpose = "area") on projected input logged nothing on
# an NA figure, and summarize_by_cell(area = TRUE) refused a single cell as
# "could not be measured", though the cell has both areas.
# ensure_projected() now measures such a layer on the 6 by 6 grid over its
# bounding box and says so when that gives no figure either;
# summarize_by_cell() measures a layer of one distinct cell on that grid.
# ---------------------------------------------------------------------------

quiet <- function(expr) suppressWarnings(suppressMessages(expr))
square <- function(x0, y0, w, h = w)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + h), c(x0, y0 + h), c(x0, y0))))
# A box densified at half a degree, so that its edges in Web Mercator follow
# the parallels and meridians.
dense_box <- function(x0, y0, w, h) {
  xs <- seq(x0, x0 + w, by = 0.5); ys <- seq(y0, y0 + h, by = 0.5)
  sf::st_polygon(list(rbind(cbind(xs, y0), cbind(x0 + w, ys[-1]), cbind(rev(xs)[-1], y0 + h),
                            cbind(x0, rev(ys)[-1]))))
}
# A bowtie: a ring whose edges cross, with lobes of 1.125e7 and 3.125e7 m^2.
bowtie <- function(x0, y0)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + 1e4, y0 + 1e4), c(x0 + 1e4, y0),
                            c(x0, y0 + 6e3), c(x0, y0))))
with_repeat <- function(p) {
  m <- unclass(p)[[1]]
  sf::st_polygon(list(rbind(m[1:2, ], m[2, ], m[-(1:2), ])))
}
# Rows at the points on the cells' surfaces, one per non-empty cell.
rows_in <- function(cells) {
  g <- sf::st_geometry(cells)
  ne <- !sf::st_is_empty(g)
  sf::st_sf(cell_id = cells$cell_id[ne], z = seq_len(sum(ne)),
            geometry = quiet(sf::st_point_on_surface(g[ne])))
}
sbc <- function(cells)
  quiet(summarize_by_cell(rows_in(cells), "z", id_col = "cell_id", cells_sf = cells,
                          estimate = "mean", area = TRUE))
area_lines <- function(x) {
  lines <- capture_spatialkit_log(out <- quiet(ensure_projected(x, purpose = "area")))
  expect_identical(out, x)
  grep("ensure_projected(purpose = \"area\")", lines, fixed = TRUE, value = TRUE)
}
merc_one <- function()
  sf::st_sf(cell_id = 1, geometry = sf::st_transform(sf::st_sfc(dense_box(-100, 30, 10, 30),
                                                                crs = 4326), 3857))

test_that("one polygon in Web Mercator is measured on the grid over its box", {
  one <- merc_one()
  expect_true(is.na(.crs_area_error(one)))
  e <- .crs_area_error(one, grid = TRUE)
  expect_gt(e, 0.6)
  expect_identical(.crs_area_error(one, one_by_grid = TRUE), e)
  l <- area_lines(one)
  expect_length(l, 1L)
  expect_match(l, sprintf("distorts areas across this extent by up to %.1f%%", 100 * e),
               fixed = TRUE)
  expect_match(l, "65.2%", fixed = TRUE)
  expect_error(sbc(one), sprintf("distorts areas across the cells by up to %.1f%%", 100 * e),
               fixed = TRUE)
  # A bowtie over the same box: s2 refuses it, the grid measures the CRS.
  b <- sf::st_bbox(one)
  bt <- sf::st_sf(cell_id = 1, geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(b[[1]], b[[2]]), c(b[[3]], b[[4]]), c(b[[3]], b[[2]]), c(b[[1]], b[[4]]),
    c(b[[1]], b[[2]])))), crs = 3857))
  expect_false(sf::st_is_valid(bt))
  expect_true(is.na(.crs_area_error(bt)))
  l <- area_lines(bt)
  expect_length(l, 1L)
  expect_match(l, sprintf("by up to %.1f%%", 100 * .crs_area_error(bt, grid = TRUE)), fixed = TRUE)
  # summarize_by_cell() refuses it as a cell GEOS calls invalid.
  expect_error(sbc(bt), "the one cell of `cells_sf` is not a valid polygon in its CRS (cell_id 1;",
               fixed = TRUE)
})

test_that("a projected layer with a polygon s2 refuses is measured on the grid", {
  g <- sf::st_transform(sf::st_sfc(lapply(0:19, function(k) square(-100, 30 + 2 * k, 1)),
                                   crs = 4326), 3857)
  b <- sf::st_bbox(g[[20]])
  g2 <- g
  g2[[20]] <- sf::st_polygon(list(rbind(c(b[1], b[2]), c(b[3], b[4]), c(b[3], b[2]),
                                        c(b[1], b[2] + 0.6 * (b[4] - b[2])), c(b[1], b[2]))))
  lay <- sf::st_sf(cell_id = 1:20, geometry = sf::st_sfc(g2, crs = 3857))
  expect_identical(sum(sf::st_is_valid(lay)), 19L)
  expect_true(is.na(.crs_area_error(lay)))
  e <- .crs_area_error(lay, grid = TRUE)
  expect_gt(e, 1)
  l <- area_lines(lay)
  expect_length(l, 1L)
  expect_match(l, sprintf("by up to %.1f%%", 100 * e), fixed = TRUE)
  expect_match(l, "124.4%", fixed = TRUE)
  # summarize_by_cell() refuses the bowtie by name.
  expect_error(sbc(lay), "1 of the 20 cells of `cells_sf` is not a valid polygon in its CRS (cell_id 20;",
               fixed = TRUE)
})

test_that("a layer on which the grid gives no figure either says so", {
  p1 <- sf::st_sf(id = 1, geometry = sf::st_sfc(sf::st_point(c(4e5, 4.4e6)), crs = 32615))
  l1 <- sf::st_sf(id = 1, geometry = sf::st_sfc(sf::st_linestring(rbind(c(4e5, 4.4e6),
                                                                        c(4.1e5, 4.4e6))), crs = 32615))
  z1 <- sf::st_sf(id = 1, geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(4e5, 4.4e6), c(4.1e5, 4.4e6), c(4e5, 4.4e6)))), crs = 32615))
  for (x in list(p1, l1, z1)) {
    expect_true(is.na(.crs_area_error(x, grid = TRUE)))
    l <- area_lines(x)
    expect_length(l, 1L)
    expect_match(l, paste0("ensure_projected(purpose = \"area\"): the area distortion of the ",
                           "input's CRS (EPSG:32615) across this extent could not be measured"),
                 fixed = TRUE)
    expect_no_match(l, "distorts areas", fixed = TRUE)
  }
})

test_that("layers whose polygons give a figure log as before", {
  two <- sf::st_sf(cell_id = 1:2, geometry = sf::st_transform(sf::st_sfc(
    list(dense_box(-100, 30, 10, 15), dense_box(-100, 45, 10, 15)), crs = 4326), 3857))
  set.seed(2)
  conus <- sf::st_transform(sf::st_as_sf(data.frame(lon = runif(200, -124, -67),
                                                    lat = runif(200, 25, 49)),
                                         coords = c("lon", "lat"), crs = 4326), 32614)
  four <- sf::st_sf(cell_id = 1:4, geometry = sf::st_sfc(lapply(0:3, function(k)
    square(4e5 + 1e4 * k, 4.4e6, 1e4)), crs = 32615))
  msg <- function(x, e) sprintf(paste0(
    "ensure_projected(purpose = \"area\"): the input's CRS (%s) distorts areas across this ",
    "extent by up to %.1f%% (planar against geodesic, spread over probe polygons)"),
    spatialkit:::.fold_crs_label(x), 100 * e)
  e2 <- .crs_area_error(two)
  expect_equal(e2, 0.2619, tolerance = 1e-3)
  expect_match(area_lines(two), msg(two, e2), fixed = TRUE)
  ec <- .crs_area_error(conus)
  expect_equal(ec, 0.1288, tolerance = 1e-3)
  expect_match(area_lines(conus), msg(conus, ec), fixed = TRUE)
  expect_length(area_lines(four), 0L)
  # one_by_grid changes nothing on a layer of two polygons or more.
  for (x in list(two, four, conus))
    expect_identical(.crs_area_error(x, one_by_grid = TRUE), .crs_area_error(x))
})

test_that("summarize_by_cell(area = TRUE) gives one projected cell its density", {
  u1 <- sf::st_sf(cell_id = 1, geometry = sf::st_sfc(square(4e5, 4.4e6, 1e4), crs = 32615))
  g <- .crs_area_error(u1, grid = TRUE)
  expect_true(is.finite(g))
  expect_lt(g, .area_error_tol)
  r <- sbc(u1)
  expect_equal(r$cell_area, 1e8, tolerance = 1e-12)
  expect_equal(r$n_per_area, 1e-8, tolerance = 1e-12)
  expect_identical(attr(r, "area_error"), g)
  # Beside an empty row.
  u1e <- sf::st_sf(cell_id = 1:2, geometry = sf::st_sfc(list(square(4e5, 4.4e6, 1e4),
                                                             sf::st_polygon()), crs = 32615))
  r <- sbc(u1e)
  expect_equal(r$cell_area[r$cell_id == 1], 1e8, tolerance = 1e-12)
  expect_true(is.na(r$cell_area[r$cell_id == 2]))
  expect_identical(attr(r, "area_error"), g)
  # The same square with a vertex listed twice: GEOS accepts it, s2 refuses
  # it in degrees until the repeat is dropped.
  ur <- sf::st_sf(cell_id = 1, geometry = sf::st_sfc(with_repeat(square(4e5, 4.4e6, 1e4)),
                                                     crs = 32615))
  expect_true(sf::st_is_valid(ur))
  r <- sbc(ur)
  expect_equal(r$cell_area, 1e8, tolerance = 1e-12)
  expect_identical(attr(r, "area_error"), g)
  # The cell listed twice under two IDs is measured as the cell alone.
  u2 <- sf::st_sf(cell_id = 1:2, geometry = rep(sf::st_geometry(u1), 2))
  r <- sbc(u2)
  expect_identical(attr(r, "area_error"), g)
})

test_that("one bowtie cell is refused, and ensure_projected() measures the grid", {
  bt <- sf::st_sf(cell_id = 1, geometry = sf::st_sfc(bowtie(4e5, 4.4e6), crs = 32615))
  expect_true(is.na(.crs_area_error(bt, one_by_grid = TRUE)))
  expect_error(sbc(bt), "Repair the cells with sf::st_make_valid() and pass them again.",
               fixed = TRUE)
  # The grid over it gives 2.8e-05, under the tolerance: no line.
  expect_lt(.crs_area_error(bt, grid = TRUE), .area_error_tol)
  expect_length(area_lines(bt), 0L)
})

test_that("summarize_by_cell(area = TRUE) gives one lon/lat cell its density", {
  ll1 <- sf::st_sf(cell_id = 1, geometry = sf::st_sfc(square(-100, 40, 0.1), crs = 4326))
  expect_true(is.na(.crs_area_error(ll1)))
  r <- sbc(ll1)
  expect_identical(attr(r, "area_error"), 0)
  expect_identical(r$cell_area, as.numeric(sf::st_area(ll1)))
})

test_that("the help says how one polygon is measured", {
  root <- testthat::test_path("..", "..")
  if (!(file.exists(file.path(root, "DESCRIPTION")) && dir.exists(file.path(root, "man"))))
    skip("the help's source is not here to read")
  read <- function(...) gsub("[[:space:]]+", " ",
                             paste(readLines(file.path(root, ...), warn = FALSE), collapse = " "))
  ep <- read("man", "ensure_projected.Rd")
  expect_match(ep, paste0("a layer of one distinct polygon, which gives one ratio and no ",
                          "spread, and a layer with a polygon s2 cannot measure in degrees"),
               fixed = TRUE)
  expect_match(ep, "a warning says that the distortion could not be measured", fixed = TRUE)
  sb <- read("man", "summarize_by_cell.Rd")
  expect_match(sb, paste0("A layer of one distinct cell gives one ratio and no spread, so its ",
                          "distortion is measured on a 6 by 6 grid over the cell's bounding box."),
               fixed = TRUE)
})
