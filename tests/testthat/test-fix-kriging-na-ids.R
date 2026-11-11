# tests/testthat/test-fix-kriging-na-ids.R
# ---------------------------------------------------------------------------
# kriging_adequacy() on a `cells_sf` with a cell whose ID is NA.  The
# per-cell counts were tabulated over the cells' IDs, which hold no level
# for NA, so they came out one short for each such cell and the call stopped
# with base R's "replacement has 3 rows, data has 4" ("0 rows" when every ID
# was NA).  A cell without an ID now comes back with n 0 and NA plain
# summaries and is kriged as any other; a layer in which no cell has an ID
# is refused in words.  A double NaN is no ID either, as in cell_values_at().
# ---------------------------------------------------------------------------

skip_if_not_installed("gstat")

.kn_catch <- function(expr) {
  got <- list()
  value <- logger::with_log_threshold(
    withCallingHandlers(expr, warning = function(w) {
      got[[length(got) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }), threshold = logger::FATAL, namespace = "spatialkit", index = 2)
  list(value = value, texts = vapply(got, conditionMessage, character(1)))
}

# Four 500 m cells A to D in EPSG:32632 and 80 points assigned to them, with
# a smooth response; the folds are given as labels so every call runs the
# same cross-validation.
.kn_layers <- function() {
  sq <- function(x0, y0) sf::st_polygon(list(rbind(c(x0, y0), c(x0 + 500, y0),
                                                   c(x0 + 500, y0 + 500), c(x0, y0 + 500),
                                                   c(x0, y0))))
  cells <- sf::st_sf(poly_id = c("A", "B", "C", "D"),
                     geometry = sf::st_sfc(sq(0, 0), sq(500, 0), sq(0, 500), sq(500, 500),
                                           crs = 32632))
  xy <- withr::with_seed(3, cbind(stats::runif(80, 0, 1000), stats::runif(80, 0, 1000)))
  z <- sin(xy[, 1] / 300) + cos(xy[, 2] / 300) +
    withr::with_seed(4, stats::rnorm(80, 0, 0.2))
  pts <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2], z = z), coords = c("x", "y"),
                      crs = 32632)
  pts$poly_id <- cells$poly_id[unlist(sf::st_intersects(pts, cells))]
  list(cells = cells, pts = pts, folds = rep_len(1:4, 80))
}

.kn_run <- function(pts, cells, folds)
  .kn_catch(kriging_adequacy(pts, "z", cells, folds = folds, nmax = 20L))

.kn_plain <- c("n", "mean", "se", "plain_var")
.kn_kr <- c("kr_pred", "kr_var", "kr_ratio", "kr_gain", "kr_shift", "kr_n_used")

test_that("kriging_adequacy() returns a cell whose ID is NA with n 0 and NA plain summaries", {
  L <- .kn_layers()
  full <- .kn_run(L$pts, L$cells, L$folds)$value
  c1 <- L$cells; c1$poly_id[1] <- NA
  got <- .kn_run(L$pts, c1, L$folds)
  ka <- got$value
  expect_s3_class(ka, "kriging_adequacy")
  expect_identical(nrow(ka), 4L)
  expect_identical(ka$poly_id, c(NA, "B", "C", "D"))
  # The points of cell A name an ID no cell has now: said as before.
  expect_true(any(grepl("22 of the 80 point\\(s\\) with a cell ID carry one that matches no `cells_sf\\$poly_id` \\(A\\)",
                        got$texts)))
  # The cell without an ID holds no point, and is kriged from its neighbours.
  expect_identical(ka$n[1], 0L)
  expect_true(all(is.na(c(ka$mean[1], ka$se[1], ka$plain_var[1], ka$kr_gain[1],
                          ka$kr_shift[1]))))
  expect_true(is.finite(ka$kr_pred[1]) && is.finite(ka$kr_var[1]) && ka$kr_var[1] > 0)
  # The other cells' plain summaries are those of the layer with every ID.
  plain <- function(x) sf::st_drop_geometry(x)[2:4, .kn_plain]
  expect_identical(plain(ka), plain(full))
  # Their kriging is the same set of systems, solved by gstat in another
  # batch of cells; compared to 1e-10, well above the rounding of a solve
  # of 20 locations.
  kr <- function(x) sf::st_drop_geometry(x)[2:4, .kn_kr]
  expect_equal(kr(ka), kr(full), tolerance = 1e-10)
  expect_equal(attr(ka, "cv")[c("zscore_var", "zscore_mean", "rmse", "n_pred")],
               attr(full, "cv")[c("zscore_var", "zscore_mean", "rmse", "n_pred")],
               tolerance = 1e-10)
})

test_that("kriging_adequacy() counts no point with an NA ID in a cell whose ID is NA", {
  # match() pairs NA with NA, so points without an ID were taken for the
  # NA cell's own when its plain_var was computed.
  L <- .kn_layers()
  c1 <- L$cells; c1$poly_id[1] <- NA
  p1 <- L$pts; p1$poly_id[p1$poly_id == "A"] <- NA
  ka <- .kn_run(p1, c1, L$folds)$value
  expect_identical(ka$n, c(0L, 22L, 18L, 18L))
  expect_true(all(is.na(c(ka$mean[1], ka$se[1], ka$plain_var[1]))))
  expect_true(all(is.finite(ka$plain_var[2:4])))
})

test_that("kriging_adequacy() refuses a cells_sf in which no cell has an ID", {
  L <- .kn_layers()
  c0 <- L$cells; c0$poly_id <- NA_character_
  err <- tryCatch(.kn_run(L$pts, c0, L$folds), error = function(e) conditionMessage(e))
  expect_type(err, "character")
  expect_match(err, "kriging_adequacy(): `cells_sf` has no cell ID", fixed = TRUE)
  expect_match(err, "'poly_id'", fixed = TRUE)
  expect_match(err, "`id_col`", fixed = TRUE)
  # NaN in a double ID column is no ID either.
  c0$poly_id <- NaN
  err <- tryCatch(.kn_run(L$pts, c0, L$folds), error = function(e) conditionMessage(e))
  expect_type(err, "character")
  expect_match(err, "kriging_adequacy(): `cells_sf` has no cell ID", fixed = TRUE)
})

test_that("kriging_adequacy() reads a double NaN ID as no ID on the cells and the points", {
  L <- .kn_layers()
  cn <- L$cells; cn$poly_id <- c(NaN, 2, 3, 4)
  pn <- L$pts; pn$poly_id <- match(L$pts$poly_id, L$cells$poly_id) + 0
  pn$poly_id[pn$poly_id == 1] <- NaN
  got <- .kn_run(pn, cn, L$folds)
  ka <- got$value
  # Written as "NaN", the cell and its points were matched as one cell.
  expect_identical(ka$n, c(0L, 22L, 18L, 18L))
  expect_true(all(is.na(c(ka$mean[1], ka$se[1], ka$plain_var[1]))))
  expect_false(any(grepl("matches no `cells_sf", got$texts)))
  # The cells with an ID are those of the layer keyed by letters.
  full <- .kn_run(L$pts, L$cells, L$folds)$value
  expect_identical(sf::st_drop_geometry(ka)[2:4, .kn_plain],
                   sf::st_drop_geometry(full)[2:4, .kn_plain])
})

test_that("kriging_adequacy() keeps each cell's own counts when two of four cell IDs are NA", {
  # The counts of the two cells with an ID were recycled over the four rows,
  # with no error: cell B took 18 points and cell D's mean.
  L <- .kn_layers()
  full <- .kn_run(L$pts, L$cells, L$folds)$value
  c2 <- L$cells; c2$poly_id[c(1, 3)] <- NA
  ka <- .kn_run(L$pts, c2, L$folds)$value
  expect_identical(nrow(ka), 4L)
  expect_identical(ka$n, c(0L, 22L, 0L, 18L))
  expect_true(all(is.na(c(ka$mean[c(1, 3)], ka$se[c(1, 3)], ka$plain_var[c(1, 3)],
                          ka$kr_shift[c(1, 3)]))))
  plain <- function(x) sf::st_drop_geometry(x)[c(2, 4), .kn_plain]
  expect_identical(plain(ka), plain(full))
  # Kriging compared to 1e-10, as above: the same systems solved in another
  # batch of cells.
  kr <- function(x) sf::st_drop_geometry(x)[c(2, 4), .kn_kr]
  expect_equal(kr(ka), kr(full), tolerance = 1e-10)
})

test_that("kriging_adequacy() refuses a cells_sf in which only a factor's NA level is an ID", {
  # is.na() is FALSE on a factor's NA level, though the level is no ID.
  L <- .kn_layers()
  cf <- L$cells; cf$poly_id <- factor(rep(NA, 4), exclude = NULL)
  pf <- L$pts; pf$poly_id <- factor(pf$poly_id)
  err <- tryCatch(.kn_run(pf, cf, L$folds), error = function(e) conditionMessage(e))
  expect_type(err, "character")
  expect_match(err, "kriging_adequacy(): `cells_sf` has no cell ID", fixed = TRUE)
})

test_that("kriging_adequacy() refuses two cells that share an ID, in words", {
  L <- .kn_layers()
  cd <- rbind(L$cells, L$cells[1, ])
  err <- tryCatch(.kn_run(L$pts, cd, L$folds), error = function(e) conditionMessage(e))
  expect_type(err, "character")
  expect_match(err, "kriging_adequacy(): `cells_sf$poly_id` holds 1 duplicated value(s) (A)",
               fixed = TRUE)
  expect_match(err, "ensure_stable_poly_id()", fixed = TRUE)
  # Read as text, as the IDs are matched: a double 1 and 1 are one ID, and
  # two NA IDs are no duplicate.
  cn <- L$cells; cn$poly_id <- c(1, 1, 3, 4)
  pn <- L$pts; pn$poly_id <- match(L$pts$poly_id, L$cells$poly_id) + 0
  err <- tryCatch(.kn_run(pn, cn, L$folds), error = function(e) conditionMessage(e))
  expect_match(err, "holds 1 duplicated value(s) (1)", fixed = TRUE)
  c2 <- L$cells; c2$poly_id[c(1, 3)] <- NA
  expect_s3_class(.kn_run(L$pts, c2, L$folds)$value, "kriging_adequacy")
})
