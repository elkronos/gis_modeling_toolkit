# tests/testthat/test-fix-sbc-lonlat-cells-degenerate.R
# ---------------------------------------------------------------------------
# summarize_by_cell(area = TRUE) takes the areas of lon/lat cells on the
# sphere, where s2 refuses a cell with two corners it reads as one point.
# When the area check measured such a cell, the request was refused as
# "could not be measured", which names neither the cell nor a way on; above
# 200 cells, when the check passed over it, the request ran through the
# summary and stopped at the cells' areas with s2's bare message.  Every
# cell's area is now taken before any summary, and a cell s2 refuses stops
# the request with the error of class "spatialkit_cells_degenerate", which
# names the rows, the cause and the transform to a projected CRS.
# ---------------------------------------------------------------------------

.sl_s2 <- function(on = TRUE, env = parent.frame()) {
  old <- suppressMessages(sf::sf_use_s2(on))
  withr::defer(suppressMessages(sf::sf_use_s2(old)), envir = env)
}
.sl_quiet <- function(expr) suppressWarnings(suppressMessages(expr))
.sl_catch <- function(expr) tryCatch(.sl_quiet(expr), error = function(e) e)
.sl_square <- function(x0, y0, w = 0.1)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + w), c(x0, y0 + w), c(x0, y0))))
.sl_repeat <- function(p) {
  m <- unclass(p)[[1]]
  sf::st_polygon(list(rbind(m[1:2, ], m[2, ], m[-(1:2), ])))
}
# `n` squares of 0.1 degree, `across` to a row, the second vertex of the
# squares `repeat_in` listed twice (GEOS accepts that ring, s2 refuses it).
.sl_cells <- function(n = 30L, across = 6L, repeat_in = integer(0), bowtie_in = integer(0)) {
  k <- seq_len(n) - 1L
  g <- lapply(seq_len(n), function(i) {
    x0 <- -100 + 0.1 * (k[i] %% across); y0 <- 40 + 0.1 * (k[i] %/% across)
    p <- .sl_square(x0, y0)
    if (i %in% repeat_in) p <- .sl_repeat(p)
    if (i %in% bowtie_in)
      p <- sf::st_polygon(list(rbind(c(x0, y0), c(x0 + 0.1, y0 + 0.1), c(x0 + 0.1, y0),
                                     c(x0, y0 + 0.06), c(x0, y0))))
    p
  })
  sf::st_sf(cell_id = seq_len(n), geometry = sf::st_sfc(g, crs = 4326))
}
# A row a little off the lower left corner of each cell.
.sl_rows <- function(cells) {
  bb <- do.call(rbind, lapply(sf::st_geometry(cells), function(p) as.numeric(sf::st_bbox(p))))
  sf::st_as_sf(data.frame(cell_id = cells$cell_id, z = seq_len(nrow(cells)),
                          x = bb[, 1] + 0.03, y = bb[, 2] + 0.04),
               coords = c("x", "y"), crs = 4326)
}
.sl_sbc <- function(rows, cells, area = TRUE)
  summarize_by_cell(rows, "z", id_col = "cell_id", cells_sf = cells, estimate = "mean",
                    area = area)
# The probes of the area check on a layer of more than 200 polygons.
.sl_probes <- function(g, max_n = 200L) {
  k  <- .polygon_sort_keys(g)
  rk <- order(k[, 1], k[, 2], k[, 3], k[, 4], k[, 5], k[, 6], k[, 7], k[, 8], k[, 9])
  ks <- k[rk, , drop = FALSE]
  at <- c(round(seq(1, length(g), length.out = max_n)), which.min(ks[, "xmin"]),
          which.min(ks[, "ymin"]), which.max(ks[, "xmax"]), which.max(ks[, "ymax"]))
  rk[sort(unique(at))]
}

test_that("lon/lat cells s2 refuses stop the request with the cells' error before any summary", {
  .sl_s2(TRUE)
  cells <- .sl_cells(repeat_in = 17L)
  rows <- .sl_rows(cells)
  expect_identical(which(!sf::st_is_valid(cells)), 17L)
  err <- .sl_catch(.sl_sbc(rows, cells))
  expect_s3_class(err, "spatialkit_cells_degenerate")
  expect_identical(err$rows, 17L)
  expect_identical(err$n, 30L)
  expect_match(err$reasons, "duplicate vertex", fixed = TRUE)
  expect_match(conditionMessage(err$parent), "Loop 0 is not valid", fixed = TRUE)
  expect_null(conditionCall(err))
  msg <- conditionMessage(err)
  expect_match(msg, paste0("summarize_by_cell(): s2, the engine sf uses for geometry in longitude ",
                           "and latitude, refuses 1 of the 30 polygons of `cells_sf` (row 17; "),
               fixed = TRUE)
  expect_match(msg, "so the measurement of the cells' areas (`area = TRUE`) cannot run in degrees.",
               fixed = TRUE)
  expect_match(msg, "It has two corners that s2 reads as one point.", fixed = TRUE)
  expect_match(msg, paste0("Transform `cells_sf` to a projected CRS, such as the one its cells ",
                           "were built in (sf::st_transform()) or the equal-area one ",
                           "ensure_projected(purpose = \"area\") picks for it, and call ",
                           "summarize_by_cell() again; the cells' areas are then planar areas ",
                           "in that CRS."), fixed = TRUE)
  expect_no_match(msg, "could not be measured", fixed = TRUE)
  expect_no_match(msg, "polygons_sf", fixed = TRUE)
  expect_no_match(msg, "assign again", fixed = TRUE)
  # Before any summary: the area check is not reached.
  local_mocked_bindings(.crs_area_error = function(...) stop("area check reached"))
  expect_s3_class(.sl_catch(.sl_sbc(rows, cells)), "spatialkit_cells_degenerate")
})

test_that("the remedy runs, and area = FALSE is not affected", {
  .sl_s2(TRUE)
  cells <- .sl_cells(repeat_in = 17L)
  rows <- .sl_rows(cells)
  expect_s3_class(.sl_quiet(.sl_sbc(rows, cells, area = FALSE)), "sf")
  for (pc in list(sf::st_transform(cells, 32615),
                  .sl_quiet(ensure_projected(cells, purpose = "area")))) {
    r <- .sl_quiet(.sl_sbc(rows, pc))
    expect_identical(nrow(r), 30L)
    expect_true(all(is.finite(r$cell_area) & r$cell_area > 0))
    expect_equal(r$n_per_area, r$n / r$cell_area)
    expect_lt(attr(r, "area_error"), .area_error_tol)
  }
})

test_that("above 200 cells, a refused cell gives the same error whether or not the check measures it", {
  .sl_s2(TRUE)
  plain <- .sl_cells(300L, 20L)
  probes <- .sl_probes(sf::st_geometry(plain))
  expect_false(30L %in% probes)
  expect_true(281L %in% probes)
  for (bad in c(30L, 281L)) {
    cells <- .sl_cells(300L, 20L, repeat_in = bad)
    rows <- .sl_rows(cells)
    err <- .sl_catch(.sl_sbc(rows, cells))
    expect_s3_class(err, "spatialkit_cells_degenerate")
    expect_identical(err$rows, bad)
    expect_identical(err$n, 300L)
    expect_match(conditionMessage(err), "refuses 1 of the 300 polygons of `cells_sf`", fixed = TRUE)
    r <- .sl_quiet(.sl_sbc(rows, sf::st_transform(cells, 32615)))
    expect_true(all(is.finite(r$cell_area)))
  }
  # Not a probe: the figure is 0 and the request stopped after the summary
  # with s2's bare message.
  expect_identical(.crs_area_error(.sl_cells(300L, 20L, repeat_in = 30L)), 0)
})

test_that("the k-means cells of a lattice returned in lon/lat stop with the error", {
  .sl_s2(TRUE)
  withr::local_seed(9001, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
                    .rng_sample_kind = "Rejection")
  g <- expand.grid(x = 4e5 + (0:29) * 30, y = 4.4e6 + (0:19) * 30)
  z <- sin(g$x / 180) + cos(g$y / 150) + stats::rnorm(nrow(g), 0, 0.4)
  pts <- sf::st_transform(sf::st_as_sf(data.frame(x = g$x, y = g$y, z = z),
                                       coords = c("x", "y"), crs = 32615), 4326)
  prof <- .sl_quiet(resolution_profile(pts, "z", partition = "kmeans", nstart = 3, sac = "auto",
                                       seed = 123, estimate = "pooled", obs_per_cell = 4.5))
  sel <- select_resolution(prof, criterion = "size")
  seeds <- .sl_quiet(get_voronoi_seeds(method = "kmeans", n = sel, sample_points = pts))
  cells <- .sl_quiet(build_tessellation(seeds, method = "voronoi", crs = 4326, quiet = TRUE))$cells
  bad <- which(!sf::st_is_valid(cells))
  skip_if(!length(bad), "s2 refuses none of these cells on this platform")
  asg <- .sl_quiet(assign_features_to_polygons(pts, cells, polygon_id_col = "cell_id", quiet = TRUE))
  err <- .sl_catch(summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = cells,
                                     estimate = "mean", area = TRUE))
  expect_s3_class(err, "spatialkit_cells_degenerate")
  expect_identical(err$rows, bad)
  expect_identical(err$n, nrow(cells))
  expect_match(conditionMessage(err), "two corners that s2 reads as one point.", fixed = TRUE)
  rec <- attr(cells, "count_rule")
  r <- .sl_quiet(summarize_by_cell(asg, "z", id_col = "cell_id",
                                   cells_sf = sf::st_transform(cells, rec$crs),
                                   estimate = "mean", area = TRUE))
  expect_true(is.finite(attr(r, "area_error")))
  expect_lt(attr(r, "area_error"), .area_error_tol)
})

test_that("lon/lat cells s2 accepts get the density and the figure they got", {
  .sl_s2(TRUE)
  plain <- .sl_cells()
  r <- .sl_quiet(.sl_sbc(.sl_rows(plain), plain))
  expect_identical(attr(r, "area_error"), .crs_area_error(plain))
  expect_identical(attr(r, "area_error"), 0)
  expect_identical(r$cell_area, as.numeric(sf::st_area(plain))[match(r$cell_id, plain$cell_id)])
})

test_that("a lon/lat bowtie gets the error with s2's reason", {
  .sl_s2(TRUE)
  cells <- .sl_cells(bowtie_in = 11L)
  err <- .sl_catch(.sl_sbc(.sl_rows(cells), cells))
  expect_s3_class(err, "spatialkit_cells_degenerate")
  expect_identical(err$rows, 11L)
  expect_match(conditionMessage(err), "crosses", fixed = TRUE)
  expect_no_match(conditionMessage(err), "two corners that s2 reads as one point", fixed = TRUE)
})

test_that("with s2 off and no lwgeom, lon/lat cells are refused as before", {
  skip_if(nzchar(system.file(package = "lwgeom")), "lwgeom is installed")
  .sl_s2(FALSE)
  for (cells in list(.sl_cells(), .sl_cells(repeat_in = 17L))) {
    err <- .sl_catch(.sl_sbc(.sl_rows(cells), cells))
    expect_false(inherits(err, "spatialkit_cells_degenerate"))
    expect_match(conditionMessage(err), "could not be measured", fixed = TRUE)
  }
  expect_false(sf::sf_use_s2())
})

test_that("the help and NEWS say what lon/lat cells s2 refuses get", {
  root <- testthat::test_path("..", "..")
  if (!(file.exists(file.path(root, "DESCRIPTION")) && dir.exists(file.path(root, "man"))))
    skip("the help's source is not here to read")
  read <- function(...) gsub("[[:space:]]+", " ",
                             paste(readLines(file.path(root, ...), warn = FALSE), collapse = " "))
  rd <- read("man", "summarize_by_cell.Rd")
  expect_match(rd, "on the cells s2 measures the distortion check passes by construction",
               fixed = TRUE)
  expect_match(rd, paste0("Every cell's area is taken before any summary, and when s2 refuses a ",
                          "cell that is not valid on the sphere"), fixed = TRUE)
  expect_no_match(rd, "distortion check passes by construction; with s2 switched off", fixed = TRUE)
  cv <- read("man", "create_voronoi_polygons.Rd")
  expect_match(cv, "with \\code{area = TRUE}, which takes the areas of such cells, stops with the same error before any summary.",
               fixed = TRUE)
  news <- read("NEWS.md")
  expect_match(news, paste0("`summarize_by_cell(area = TRUE)` stopped without the cause on lon/lat ",
                            "cells that s2 refuses."), fixed = TRUE)
})
