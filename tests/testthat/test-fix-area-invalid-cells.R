# tests/testthat/test-fix-area-invalid-cells.R
# ---------------------------------------------------------------------------
# summarize_by_cell(area = TRUE) divides by each cell's planar area, which
# GEOS takes on the cell's rings as they are: a ring whose edges cross gets
# the signed sum of its lobes.  The area check measures up to 200 cells, so
# a bowtie among them made the request refused and one outside them got a
# density over that area.  Every non-empty projected cell is now checked
# with sf::st_is_valid() in its CRS before any summary, and a cell GEOS
# calls invalid makes the request refused, naming the cells, GEOS's reason
# and the remedy, sf::st_make_valid().
# ---------------------------------------------------------------------------

.ai_quiet <- function(expr) suppressWarnings(suppressMessages(expr))
.ai_square <- function(x0, y0, w = 1e4)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + w), c(x0, y0 + w), c(x0, y0))))
# A bowtie of 10 km: lobes of 1.125e7 and 3.125e7 m^2; GEOS's area 2e7.
.ai_bowtie <- function(x0, y0)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + 1e4, y0 + 1e4), c(x0 + 1e4, y0),
                            c(x0, y0 + 6e3), c(x0, y0))))
.ai_spike <- function(x0, y0)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + 1e4, y0), c(x0 + 1e4, y0 + 1e4),
                            c(x0 + 5e3, y0 + 1e4), c(x0 + 5e3, y0 + 2e4),
                            c(x0 + 5e3, y0 + 1e4), c(x0, y0 + 1e4), c(x0, y0))))
# 400 squares of 10 km in EPSG:32615, 20 across; cells `odd` drawn by `f`.
.ai_grid <- function(odd = integer(0), f = .ai_bowtie, n = 400L, id = "cell_id") {
  k <- seq_len(n) - 1L
  g <- lapply(seq_len(n), function(i) {
    x0 <- 4e5 + 1e4 * (k[i] %% 20); y0 <- 4.4e6 + 1e4 * (k[i] %/% 20)
    if (i %in% odd) f(x0, y0) else .ai_square(x0, y0)
  })
  out <- sf::st_sf(id = seq_len(n), geometry = sf::st_sfc(g, crs = 32615))
  names(out)[1] <- id
  out
}
.ai_rows <- function(cells, id = "cell_id") {
  bb <- do.call(rbind, lapply(sf::st_geometry(cells), function(p) as.numeric(sf::st_bbox(p))))
  out <- sf::st_as_sf(data.frame(id = cells[[id]], z = seq_len(nrow(cells)),
                                 x = bb[, 1] + 900, y = bb[, 2] + 1100),
                      coords = c("x", "y"), crs = sf::st_crs(cells))
  names(out)[1] <- id
  out
}
.ai_sbc <- function(cells, id = "cell_id", rows = .ai_rows(cells, id), ...)
  .ai_quiet(summarize_by_cell(rows, "z", id_col = id, cells_sf = cells, estimate = "mean",
                              area = TRUE, ...))

test_that("a bowtie is refused whether or not the area check measures it", {
  msgs <- character()
  for (at in c(148L, 1L)) {
    cells <- .ai_grid(at)
    # The premise: the bowtie at 148 is not a probe (the figure is finite and
    # the request ran with cell_area 2e7), at 1 it is (the figure is NA).
    if (at == 148L) expect_true(is.finite(.crs_area_error(cells)))
    else expect_true(is.na(.crs_area_error(cells)))
    err <- tryCatch(.ai_sbc(cells), error = function(e) e)
    expect_s3_class(err, "error")
    msg <- conditionMessage(err)
    expect_match(msg, sprintf(paste0("summarize_by_cell(): `area = TRUE` refused: 1 of the 400 ",
                                     "cells of `cells_sf` is not a valid polygon in its CRS ",
                                     "(cell_id %d; GEOS's reason: \""), at), fixed = TRUE)
    expect_match(msg, "Self-intersection", fixed = TRUE)
    expect_match(msg, paste0("so its planar area, and a density over it, cannot be vouched for. ",
                             "Repair the cells with sf::st_make_valid() and pass them again."),
                 fixed = TRUE)
    expect_null(conditionCall(err))
    msgs <- c(msgs, sub("cell_id [0-9]+; GEOS's reason: \"[^\"]*\"", "", msg))
  }
  expect_identical(msgs[1], msgs[2])
  # Rows reversed: the same cell named.
  cells <- .ai_grid(148L)[400:1, ]
  expect_error(.ai_sbc(cells), "(cell_id 148; GEOS", fixed = TRUE)
  # Before any summary: the area check is not reached.
  local_mocked_bindings(.crs_area_error = function(...) stop("area check reached"))
  expect_error(.ai_sbc(.ai_grid(148L)), "is not a valid polygon in its CRS", fixed = TRUE)
})

test_that("small layers, spikes, overlapping parts and holes outside the shell are refused the same way", {
  four <- function(odd) {
    sq <- lapply(0:2, function(k) .ai_square(4e5 + 1e4 * k, 4.4e6))
    sf::st_sf(cell_id = 1:4, geometry = sf::st_sfc(c(sq, list(odd)), crs = 32615))
  }
  overlap <- sf::st_multipolygon(list(list(unclass(.ai_square(4.3e5, 4.4e6))[[1]]),
                                      list(unclass(.ai_square(4.35e5, 4.4e6))[[1]])))
  hole_out <- sf::st_polygon(list(unclass(.ai_square(4.3e5, 4.4e6))[[1]],
                                  unclass(.ai_square(4.5e5, 4.4e6, 1e3))[[1]]))
  for (odd in list(.ai_bowtie(4.3e5, 4.4e6), .ai_spike(4.3e5, 4.4e6), overlap, hole_out)) {
    cells <- four(odd)
    expect_identical(sf::st_is_valid(cells), c(TRUE, TRUE, TRUE, FALSE))
    expect_error(.ai_sbc(cells), "1 of the 4 cells of `cells_sf` is not a valid polygon in its CRS (cell_id 4;",
                 fixed = TRUE)
  }
  # Two invalid cells are both named; the reason is the first one's.
  cells <- .ai_grid(c(7L, 300L))
  expect_error(.ai_sbc(cells), paste0("2 of the 400 cells of `cells_sf` are not valid polygons ",
                                      "in its CRS (cell_id 7 and 300; GEOS's reason for cell_id 7: "),
               fixed = TRUE)
  expect_error(.ai_sbc(cells), "so each one's planar area", fixed = TRUE)
  # Eight name six and two more.
  cells <- .ai_grid(c(2L, 4L, 6L, 8L, 10L, 12L, 14L, 16L))
  expect_error(.ai_sbc(cells), "(cell_id 2, 4, 6, 8, 10, 12 and 2 more; ", fixed = TRUE)
  # Six are all named, with no "more"; seven name six and one more.
  err6 <- tryCatch(.ai_sbc(.ai_grid(5:10)), error = function(e) conditionMessage(e))
  expect_match(err6, "6 of the 400 cells of `cells_sf` are not valid polygons in its CRS ",
               fixed = TRUE)
  expect_match(err6, "(cell_id 5, 6, 7, 8, 9 and 10; GEOS's reason for cell_id 5: ", fixed = TRUE)
  expect_no_match(err6, "more", fixed = TRUE)
  expect_error(.ai_sbc(.ai_grid(5:11)), "(cell_id 5, 6, 7, 8, 9, 10 and 1 more; ", fixed = TRUE)
})

test_that("ensure_projected(purpose = \"area\") measures a layer with a polygon GEOS calls invalid on the grid", {
  # Three 10 km squares of UTM zone 15 beside a square with a hole outside
  # its shell.  GEOS takes the shell minus the hole (7.5e7 m^2), s2 the shell
  # plus the hole (1.25e8 m^2), so the polygons' spread was 40 percent: a
  # statement about the polygon's invalidity, not about the CRS.
  sq <- lapply(0:2, function(k) .ai_square(4e5 + 1e4 * k, 4.4e6))
  hole_out <- sf::st_polygon(list(unclass(.ai_square(4.3e5, 4.4e6))[[1]],
                                  unclass(.ai_square(4.5e5, 4.4e6, 5e3))[[1]]))
  l4 <- sf::st_sf(id = 1:4, geometry = sf::st_sfc(c(sq, list(hole_out)), crs = 32615))
  expect_identical(sf::st_is_valid(l4), c(TRUE, TRUE, TRUE, FALSE))
  expect_true(is.na(.crs_area_error(l4)))
  grid_err <- .crs_area_error(l4, grid = TRUE)
  expect_true(is.finite(grid_err) && grid_err < 0.01)
  log <- capture_spatialkit_log(out <- .ai_quiet(ensure_projected(l4, purpose = "area")))
  expect_identical(out, l4)
  expect_false(any(grepl("distorts areas", log, fixed = TRUE)))
  expect_false(any(grepl("could not be measured", log, fixed = TRUE)))
  # In Web Mercator from 30N to 60N the grid's figure is logged.
  bx <- sf::st_transform(sf::st_as_sfc(sf::st_bbox(c(xmin = -100, ymin = 30, xmax = -90, ymax = 60),
                                                   crs = 4326)), 3857)
  bxp <- unclass(bx[[1]])[[1]]
  bb <- sf::st_bbox(bx)
  w <- as.numeric(bb["xmax"] - bb["xmin"]) / 10
  hole <- rbind(c(bb["xmax"] + w, bb["ymin"]), c(bb["xmax"] + 2 * w, bb["ymin"]),
                c(bb["xmax"] + 2 * w, bb["ymin"] + w), c(bb["xmax"] + w, bb["ymin"] + w),
                c(bb["xmax"] + w, bb["ymin"]))
  bad <- sf::st_sf(id = 1, geometry = sf::st_sfc(sf::st_polygon(list(bxp, unname(hole))), crs = 3857))
  expect_false(sf::st_is_valid(bad))
  expect_true(is.na(.crs_area_error(bad)))
  log <- capture_spatialkit_log(.ai_quiet(ensure_projected(bad, purpose = "area")))
  expect_length(log, 1L)
  expect_match(log, sprintf("distorts areas across this extent by up to %.1f%%",
                            100 * .crs_area_error(bad, grid = TRUE)), fixed = TRUE)
  # The three valid squares alone are measured on themselves, as before.
  expect_true(is.finite(.crs_area_error(l4[1:3, ])))
})

test_that("cells are named by the ID column of the join, or by row", {
  cells <- .ai_grid(5L, id = "id")
  cells$id <- cells$id + 1000L
  expect_error(.ai_sbc(cells, id = "id"), "(id 1005; GEOS", fixed = TRUE)
  # The rows keyed by poly_id (the default id_col), the cells by id.
  rows <- .ai_rows(cells, "id"); names(rows)[1] <- "poly_id"
  expect_error(.ai_quiet(summarize_by_cell(rows, "z", cells_sf = cells, estimate = "mean",
                                           area = TRUE)),
               "(id 1005; GEOS", fixed = TRUE)
  # Cells with no ID column the join reads are named by row.
  bare <- cells; names(bare)[1] <- "other"
  expect_error(.ai_quiet(summarize_by_cell(rows, "z", cells_sf = bare, estimate = "mean",
                                           area = TRUE)),
               "(row 5; GEOS", fixed = TRUE)
  two <- .ai_grid(c(5L, 9L)); names(two)[1] <- "other"
  expect_error(.ai_quiet(summarize_by_cell(.ai_rows(.ai_grid(c(5L, 9L))), "z", id_col = "cell_id",
                                           cells_sf = two, estimate = "mean", area = TRUE)),
               "(rows 5 and 9; GEOS's reason for row 5: ", fixed = TRUE)
  # Valid cells without an ID column keep the join's own refusal.
  bare_ok <- .ai_grid(); names(bare_ok)[1] <- "other"
  expect_error(.ai_quiet(summarize_by_cell(.ai_rows(.ai_grid()), "z", id_col = "cell_id",
                                           cells_sf = bare_ok, estimate = "mean", area = TRUE)),
               "`area = TRUE` needs the cells' areas", fixed = TRUE)
})

test_that("the remedy runs, wherever the bowtie sits", {
  for (at in c(148L, 1L)) {
    cells <- .ai_grid(at)
    fixed <- sf::st_make_valid(cells)
    expect_true(all(sf::st_is_valid(fixed)))
    r <- .ai_sbc(fixed)
    expect_equal(r$cell_area[r$cell_id == at], 4.25e7, tolerance = 1e-9)
    expect_true(is.finite(attr(r, "area_error")))
  }
})

test_that("valid projected cells get the densities and the figure they got", {
  cells <- .ai_grid()
  r <- .ai_sbc(cells)
  expect_identical(r$cell_area, as.numeric(sf::st_area(cells))[match(r$cell_id, cells$cell_id)])
  expect_identical(attr(r, "area_error"), .crs_area_error(cells))
  # An empty cell among valid ones is not checked.
  g <- sf::st_geometry(cells); g[[3]] <- sf::st_polygon()
  withe <- sf::st_set_geometry(cells, g)
  r <- .ai_sbc(withe, rows = .ai_rows(cells))
  expect_true(is.na(r$cell_area[r$cell_id == 3]))
  expect_equal(r$cell_area[r$cell_id == 4], 1e8, tolerance = 1e-12)
  # A ring with a vertex listed twice is valid in GEOS and runs.
  m <- unclass(sf::st_geometry(cells)[[4]])[[1]]
  g <- sf::st_geometry(cells); g[[4]] <- sf::st_polygon(list(rbind(m[1:2, ], m[2, ], m[-(1:2), ])))
  r <- .ai_sbc(sf::st_set_geometry(cells, g))
  expect_equal(r$cell_area[r$cell_id == 4], 1e8, tolerance = 1e-12)
})

test_that("the help and NEWS say how invalid cells are refused", {
  root <- testthat::test_path("..", "..")
  if (!(file.exists(file.path(root, "DESCRIPTION")) && dir.exists(file.path(root, "man"))))
    skip("the help's source is not here to read")
  read <- function(...) gsub("[[:space:]]+", " ",
                             paste(readLines(file.path(root, ...), warn = FALSE), collapse = " "))
  rd <- read("man", "summarize_by_cell.Rd")
  expect_match(rd, paste0("Projected cells are first checked with \\code{sf::st_is_valid()} in ",
                          "their CRS, every non-empty cell before any summary"), fixed = TRUE)
  expect_match(rd, "\\code{sf::st_make_valid()}", fixed = TRUE)
  expect_match(rd, "Exact copies of one cell count once.", fixed = TRUE)
  ep <- read("man", "ensure_projected.Rd")
  expect_match(ep, paste0("or that \\code{\\link[sf:st_is_valid]{sf::st_is_valid()}} rejects in its CRS ",
                          "(a hole outside its shell), whose planar area is not the area it ",
                          "encloses."), fixed = TRUE)
  expect_match(ep, "(exact copies of a polygon count once)", fixed = TRUE)
  news <- read("NEWS.md")
  expect_match(news, paste0("`summarize_by_cell(area = TRUE)` gave a projected cell whose edges ",
                            "cross a density over GEOS's area of its ring."), fixed = TRUE)
  expect_match(news, paste0("A layer with a polygon GEOS calls invalid is now measured on the 6 ",
                            "by 6 grid over its extent"), fixed = TRUE)
})
