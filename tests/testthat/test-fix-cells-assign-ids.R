# tests/testthat/test-fix-cells-assign-ids.R
# ---------------------------------------------------------------------------
# assign_features_to_polygons(): where the polygon IDs come from, and what
# happens to features no polygon holds.
#
# The polygons' ID was read from a fallback column without a word (nine
# polygons with a zone attribute called `id` came back as three), polygons
# with no ID column at all were numbered in silence (their `GEOID`s unused),
# a duplicated ID merged two polygons' features, and features outside the
# layer were removed with one INFO line in the temp log.
# ---------------------------------------------------------------------------

.fa_sq <- function(x0, y0, s) sf::st_polygon(list(rbind(
  c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s), c(x0, y0 + s), c(x0, y0))))

# A 3 x 3 grid of 100 m cells with no ID column, and 300 points inside it.
.fa_grid <- function() sf::st_sf(geometry = sf::st_sfc(lapply(0:8, function(i)
  .fa_sq(5e5 + 100 * (i %% 3), 5e6 + 100 * (i %/% 3), 100)), crs = 32632))
.fa_pts <- function(n = 300, seed = 9) {
  set.seed(seed)
  sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, 300), y = 5e6 + runif(n, 0, 300)),
               coords = c("x", "y"), crs = 32632)
}

# Value, warnings (with classes and conditions) and messages of a call.
.fa_run <- function(expr) {
  w <- character(0); cls <- list(); conds <- list(); m <- character(0)
  val <- withCallingHandlers(
    expr,
    warning = function(c) {
      w <<- c(w, conditionMessage(c)); cls[[length(cls) + 1L]] <<- class(c)
      conds[[length(conds) + 1L]] <<- c
      invokeRestart("muffleWarning")
    },
    message = function(c) { m <<- c(m, conditionMessage(c)); invokeRestart("muffleMessage") })
  list(value = val, warnings = w, classes = cls, conditions = conds, messages = m)
}


# --- which column the IDs came from ------------------------------------------

test_that("a fallback ID column is named in a message", {
  pts <- .fa_pts()
  g <- .fa_grid(); g$cell_id <- 9:1
  r <- .fa_run(assign_features_to_polygons(pts, g))
  expect_length(r$warnings, 0L)
  expect_length(r$messages, 1L)
  expect_match(r$messages, "no 'poly_id' column; the polygon IDs are read from its 'cell_id' column")
  # The IDs are the cells' own, under the name asked for.
  expect_setequal(r$value$poly_id, 9:1)
  expect_identical(nrow(r$value), 300L)

  # Naming the column is not a fallback, and says nothing.
  r2 <- .fa_run(assign_features_to_polygons(pts, g, polygon_id_col = "cell_id"))
  expect_length(r2$messages, 0L)
  expect_setequal(r2$value$cell_id, 9:1)
  # Nor is the default column when the polygons carry it.
  g$poly_id <- 1:9
  expect_length(.fa_run(assign_features_to_polygons(pts, g))$messages, 0L)

  # quiet = TRUE silences the message, not the fallback.
  g$poly_id <- NULL
  r3 <- .fa_run(assign_features_to_polygons(pts, g, quiet = TRUE))
  expect_length(r3$messages, 0L)
  expect_setequal(r3$value$poly_id, 9:1)
})

test_that("polygons with no ID column are numbered, and the message says so", {
  pts <- .fa_pts()
  g <- .fa_grid(); g$GEOID <- sprintf("37%03d", 1:9)
  r <- .fa_run(assign_features_to_polygons(pts, g))
  expect_length(r$messages, 1L)
  expect_match(r$messages, "no 'poly_id' column, and none of the fallback ID columns")
  expect_match(r$messages, "numbered 1 to 9 in row order")
  expect_match(r$messages, "`polygon_id_col`")
  expect_setequal(r$value$poly_id, 1:9)
  # The way the message points to: the GEOIDs, used.
  r2 <- .fa_run(assign_features_to_polygons(pts, g, polygon_id_col = "GEOID"))
  expect_length(r2$messages, 0L)
  expect_setequal(r2$value$GEOID, sprintf("37%03d", 1:9))
})

test_that("duplicated or missing polygon IDs are refused", {
  pts <- .fa_pts()
  # A zone attribute that happens to be called `id`: six cells used to be
  # merged into three "polygons", with no condition raised.
  g1 <- .fa_grid(); g1$id <- rep(1:3, each = 3)
  expect_error(suppressMessages(assign_features_to_polygons(pts, g1)),
               "polygon ID column 'id' \\(read because `polygons_sf` has no 'poly_id' column\\) holds 6 duplicated and 0 missing value\\(s\\) among 9 polygon\\(s\\)")
  # A `poly_id` repeated once gave 8 IDs for 9 polygons.
  g2 <- .fa_grid(); g2$poly_id <- c(1:8, 8L)
  expect_error(assign_features_to_polygons(pts, g2),
               "polygon ID column 'poly_id' \\(the column named by `polygon_id_col`\\) holds 1 duplicated and 0 missing")
  # A missing ID made the features inside that polygon look unassigned.
  g3 <- .fa_grid(); g3$poly_id <- c(1:8, NA)
  expect_error(assign_features_to_polygons(pts, g3),
               "holds 0 duplicated and 1 missing value\\(s\\) among 9 polygon\\(s\\)")
  # The repair the message names works.
  g1$poly_id <- seq_len(nrow(g1))
  expect_identical(nrow(assign_features_to_polygons(pts, g1)), 300L)
})

test_that("one unit stored in several rows: the error says how to keep it one cell", {
  pts <- .fa_pts()
  # Unit "A" is two of the nine squares (a region split into POLYGON rows).
  g <- .fa_grid(); g$poly_id <- c("A", "A", "B", "C", "D", "E", "F", "G", "H")
  err <- tryCatch(assign_features_to_polygons(pts, g), error = function(e) conditionMessage(e))
  expect_match(err, "holds 1 duplicated and 0 missing value\\(s\\) among 9 polygon\\(s\\)")
  # Numbering the rows, the first repair named, would split the unit in two;
  # the message says so and gives the repair that does not.
  expect_match(err, "parts of one unit")
  expect_match(err, "dplyr::summarise(dplyr::group_by(polygons_sf, poly_id))", fixed = TRUE)
  expect_match(err, "Numbering the rows would make each part a cell of its own", fixed = TRUE)
  # That repair works: eight polygons, the unit whole, every point assigned.
  merged <- dplyr::summarise(dplyr::group_by(g, poly_id))
  expect_identical(nrow(merged), 8L)
  out <- assign_features_to_polygons(pts, merged)
  expect_identical(nrow(out), 300L)
  in_a <- lengths(sf::st_intersects(pts, sf::st_union(g[1:2, ]))) > 0L
  expect_identical(out$poly_id == "A", in_a)
  # The column is named as it has to be typed, and only duplicates get the hint.
  names(g)[names(g) == "poly_id"] <- "unit id"
  expect_error(assign_features_to_polygons(pts, g, polygon_id_col = "unit id"),
               "dplyr::group_by(polygons_sf, `unit id`)", fixed = TRUE)
  g3 <- .fa_grid(); g3$poly_id <- c(1:8, NA)
  err_na <- tryCatch(assign_features_to_polygons(pts, g3), error = function(e) conditionMessage(e))
  expect_false(grepl("parts of one unit", err_na))
})


# --- features no polygon holds -----------------------------------------------

test_that("features outside the layer are dropped with a classed warning and a record", {
  set.seed(9)
  inside <- .fa_pts()
  outside <- sf::st_as_sf(data.frame(x = 5e5 + runif(30, 300.5, 305),
                                     y = 5e6 + runif(30, 0, 300)),
                          coords = c("x", "y"), crs = 32632)
  p2 <- rbind(inside, outside)
  g <- .fa_grid(); g$poly_id <- 1:9

  r <- .fa_run(assign_features_to_polygons(p2, g))
  expect_identical(nrow(r$value), 300L)
  expect_length(r$warnings, 1L)
  expect_true("spatialkit_unassigned" %in% r$classes[[1L]])
  expect_match(r$warnings, "30 of 330 feature\\(s\\) \\(9\\.1%\\) fall inside no polygon and are left out")
  expect_identical(r$conditions[[1L]]$n, 30L)
  expect_identical(r$conditions[[1L]]$total, 330L)
  rec <- attr(r$value, "unassigned")
  expect_identical(rec$n, 30L)
  expect_identical(rec$which, 301:330)
  expect_identical(rec$n_rows, 300L)
  # Caught by class, without matching the message.
  expect_identical(tryCatch(assign_features_to_polygons(p2, g),
                            spatialkit_unassigned = function(w) "caught"), "caught")

  # Kept on request: no warning, the same record, NA IDs.
  k <- .fa_run(assign_features_to_polygons(p2, g, keep_unassigned = TRUE))
  expect_length(k$warnings, 0L)
  expect_identical(nrow(k$value), 330L)
  expect_identical(which(is.na(k$value$poly_id)), 301:330)
  expect_identical(attr(k$value, "unassigned")$n, 30L)
  expect_identical(attr(k$value, "unassigned")$which, 301:330)

  # Nothing outside: a record of zero, and silence.
  ok <- .fa_run(assign_features_to_polygons(inside, g))
  expect_length(ok$warnings, 0L)
  expect_identical(attr(ok$value, "unassigned")$n, 0L)
  expect_identical(attr(ok$value, "unassigned")$which, integer(0))
  # The ties record is still there beside it.
  expect_identical(attr(ok$value, "ties")$n_rows, 300L)
})

test_that("a subset of the result is never read as carrying its parent's unassigned record", {
  set.seed(9)
  inside <- .fa_pts(); inside$v <- rnorm(300)
  outside <- sf::st_as_sf(data.frame(x = 5e5 + runif(30, 300.5, 305),
                                     y = 5e6 + runif(30, 0, 300), v = 0),
                          coords = c("x", "y"), crs = 32632)
  g <- .fa_grid(); g$poly_id <- 1:9
  a <- suppressWarnings(assign_features_to_polygons(rbind(inside, outside), g))
  expect_identical(spatialkit:::.get_row_record(a, "unassigned")$n, 30L)
  # The record is stamped with the rows it was made for, so the package's
  # reader returns nothing for a subset, whether or not `[` removed it.
  expect_null(spatialkit:::.get_row_record(a[1:10, ], "unassigned"))
  expect_null(spatialkit:::.get_row_record(dplyr::filter(a, v > 0), "unassigned"))
  # The same rows without their geometry are the same rows.
  expect_identical(attr(sf::st_drop_geometry(a), "unassigned")$n, 30L)
})

test_that("the unassigned record is removed by a subset, as the ties record beside it is", {
  # `[` removes the records named in .row_record_attrs (R/utils.R).  Before
  # "unassigned" was listed there a 10-row subset still held the parent's
  # n = 30 and n_rows = 300.
  expect_true("unassigned" %in% spatialkit:::.row_record_attrs)
  set.seed(9)
  inside <- .fa_pts(); inside$v <- rnorm(300)
  outside <- sf::st_as_sf(data.frame(x = 5e5 + runif(30, 300.5, 305),
                                     y = 5e6 + runif(30, 0, 300), v = 0),
                          coords = c("x", "y"), crs = 32632)
  g <- .fa_grid(); g$poly_id <- 1:9
  a <- suppressWarnings(assign_features_to_polygons(rbind(inside, outside), g))
  expect_null(attr(a[1:10, ], "ties"))
  expect_null(attr(a[1:10, ], "unassigned"))
  expect_null(attr(dplyr::filter(a, v > 0), "unassigned"))
  expect_null(attr(dplyr::slice(a, 1:5), "unassigned"))
})

test_that("no feature inside any polygon keeps its warning, now with the class", {
  g <- .fa_grid(); g$poly_id <- 1:9
  far <- sf::st_as_sf(data.frame(x = 6e5 + 1:5, y = 5e6 + 1:5), coords = c("x", "y"),
                      crs = 32632)
  for (keep in c(FALSE, TRUE)) {
    r <- .fa_run(assign_features_to_polygons(far, g, keep_unassigned = keep))
    expect_length(r$warnings, 1L)
    expect_match(r$warnings, "none of the 5 feature\\(s\\) fall inside any polygon")
    expect_true("spatialkit_unassigned" %in% r$classes[[1L]])
    expect_identical(attr(r$value, "unassigned")$n, 5L)
    expect_identical(nrow(r$value), if (keep) 5L else 0L)
  }
})

test_that("polygon features that only touch the layer are counted as unassigned", {
  g <- .fa_grid(); g$poly_id <- 1:9
  # One square inside cell 1, one sharing only an edge with the grid.
  f <- sf::st_sf(v = 1:2, geometry = sf::st_sfc(.fa_sq(5e5 + 10, 5e6 + 10, 20),
                                                .fa_sq(5e5 + 300, 5e6 + 10, 20),
                                                crs = 32632))
  r <- .fa_run(assign_features_to_polygons(f, g))
  expect_identical(nrow(r$value), 1L)
  expect_true("spatialkit_unassigned" %in% r$classes[[1L]])
  expect_identical(attr(r$value, "unassigned")$which, 2L)
})
