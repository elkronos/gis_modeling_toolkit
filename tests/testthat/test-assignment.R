# ===========================================================================
# assign_features_to_polygons(): the sole producer of the `poly_id` column
# that summarize_by_cell() -- and everything downstream of it -- consumes.
#
# It had no tests at all, including none for the collision case that silently
# returned ZERO ROWS.
# ===========================================================================

.sq <- function(x0, y0, x1, y1) {
  sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1),
                            c(x0, y0))))
}

# Two NESTED polygons: id 1 covers 0-10 (area 100), id 2 covers 5-8 (area 9).
# A point at (6, 6) falls in both, which is what makes tie_break observable.
.nested <- function(ids = c(1L, 2L)) {
  sf::st_sf(poly_id = ids,
            geometry = sf::st_sfc(.sq(0, 0, 10, 10), .sq(5, 5, 8, 8),
                                  crs = 32632))
}

# Two disjoint 10 x 10 cells side by side.
.two_cells <- function() {
  sf::st_sf(poly_id = c(1L, 2L),
            geometry = sf::st_sfc(.sq(0, 0, 10, 10), .sq(10, 0, 20, 10),
                                  crs = 32632))
}

.three_points <- function() {
  sf::st_sf(v = 1:3,
            geometry = sf::st_sfc(sf::st_point(c(6, 6)),      # in both nested
                                  sf::st_point(c(1, 1)),      # in the big one
                                  sf::st_point(c(20, 20)),    # in neither
                                  crs = 32632))
}


test_that("assign_features_to_polygons attaches the polygon id and drops misses", {
  out <- assign_features_to_polygons(.three_points(), .nested())

  expect_s3_class(out, "sf")
  # The unassigned point is dropped by default.
  expect_equal(nrow(out), 2L)
  expect_equal(out$v, c(1L, 2L))
  expect_true("poly_id" %in% names(out))
  expect_false(anyNA(out$poly_id))
  # Feature attributes survive, and no bookkeeping column leaks out.
  expect_false(any(grepl("^\\.\\.", names(out))))
  expect_equal(sf::st_crs(out), sf::st_crs(32632))

  # keep_unassigned retains the miss, with an NA id.
  kept <- assign_features_to_polygons(.three_points(), .nested(),
                                      keep_unassigned = TRUE)
  expect_equal(nrow(kept), 3L)
  expect_true(is.na(kept$poly_id[3]))
})


test_that("tie_break = 'smallest_area' is independent of polygon order", {
  # The point at (6, 6) is inside BOTH nested polygons.  The smallest-area
  # rule must pick the 9-unit polygon whichever row it occupies -- the previous
  # keep-the-first-duplicate behaviour depended on the join order.
  fwd <- assign_features_to_polygons(.three_points(), .nested(c(1L, 2L)))
  rev <- assign_features_to_polygons(.three_points(),
                                     sf::st_sf(poly_id = c(2L, 1L),
                                               geometry = sf::st_sfc(
                                                 .sq(5, 5, 8, 8),
                                                 .sq(0, 0, 10, 10),
                                                 crs = 32632)))

  # In both layers the SMALL polygon carries id 2 and the big one id 1.
  expect_equal(fwd$poly_id[fwd$v == 1L], 2L)
  expect_equal(rev$poly_id[rev$v == 1L], 2L)
  # The unambiguous point is unaffected either way.
  expect_equal(fwd$poly_id[fwd$v == 2L], 1L)
  expect_equal(rev$poly_id[rev$v == 2L], 1L)
})


test_that("tie_break = 'first' keeps the first match, which is order-dependent", {
  # Documented behaviour, and the reason "smallest_area" is the default: with
  # `first` the answer for the ambiguous point follows the row order.
  big_first   <- .nested(c(1L, 2L))                      # rows: big, small
  small_first <- sf::st_sf(poly_id = c(1L, 2L),
                           geometry = sf::st_sfc(.sq(5, 5, 8, 8),
                                                 .sq(0, 0, 10, 10),
                                                 crs = 32632))

  a <- assign_features_to_polygons(.three_points(), big_first,
                                   tie_break = "first")
  b <- assign_features_to_polygons(.three_points(), small_first,
                                   tie_break = "first")

  expect_equal(a$poly_id[a$v == 1L], 1L)      # big listed first -> big wins
  expect_equal(b$poly_id[b$v == 1L], 1L)      # small listed first -> small wins

  # The two strategies genuinely differ on this input.
  smallest <- assign_features_to_polygons(.three_points(), big_first)
  expect_false(identical(a$poly_id[a$v == 1L], smallest$poly_id[smallest$v == 1L]))

  expect_error(assign_features_to_polygons(.three_points(), big_first,
                                           tie_break = "nonsense"),
               "'arg' should be one of")
})


test_that("largest = TRUE resolves a straddling polygon by overlap area", {
  cells <- .two_cells()

  # x from 7 to 14: 3 units of width in cell 1, 4 in cell 2.
  straddler <- sf::st_sf(v = 1L,
                         geometry = sf::st_sfc(.sq(7, 2, 14, 5), crs = 32632))

  with_largest <- suppressWarnings(
    assign_features_to_polygons(straddler, cells, largest = TRUE))
  without      <- suppressWarnings(
    assign_features_to_polygons(straddler, cells, largest = FALSE))

  expect_equal(nrow(with_largest), 1L)
  expect_equal(with_largest$poly_id, 2L)   # bigger overlap
  # Without it the tie-break falls back to area (the cells are equal), so the
  # first match wins -- a different answer, which is what makes the assertion
  # above discriminating.
  expect_equal(without$poly_id, 1L)

  # Mirror image: 6 to 13 leaves more in cell 1.
  other <- sf::st_sf(v = 1L,
                     geometry = sf::st_sfc(.sq(6, 2, 13, 5), crs = 32632))
  expect_equal(suppressWarnings(
    assign_features_to_polygons(other, cells, largest = TRUE))$poly_id, 1L)

  # `largest` is only applicable to polygon features; point input must not
  # error on it (st_join rejects largest= for some predicates).
  expect_equal(nrow(assign_features_to_polygons(.three_points(), cells,
                                                largest = TRUE)), 2L)
})


test_that("re-assigning an already-assigned layer drops the colliding column", {
  # The just-fixed bug.  st_join() suffixes columns present on both sides
  # (poly_id.x / poly_id.y), so the rename below it never found `poly_id` and
  # the function returned ZERO ROWS -- silently, for the completely ordinary
  # case of re-assigning points that came out of a previous call.
  first  <- assign_features_to_polygons(.three_points(), .nested())
  expect_true("poly_id" %in% names(first))

  lines <- capture_spatialkit_log(
    second <- assign_features_to_polygons(first, .two_cells())
  )
  expect_true(log_has(lines, "already carries column\\(s\\) 'poly_id'"))
  expect_true(log_has(lines, "dropping them before the spatial join"))

  # The row count is the point: 2 in, 2 out, not 0.
  expect_equal(nrow(second), nrow(first))
  expect_equal(second$v, first$v)
  # Exactly one poly_id column, holding the NEW layer's ids.
  expect_equal(sum(names(second) == "poly_id"), 1L)
  expect_false(any(grepl("poly_id\\.[xy]", names(second))))
  expect_equal(second$poly_id, c(1L, 1L))      # both points are in cell 1

  # A custom id column name collides the same way.
  named <- first
  names(named)[names(named) == "poly_id"] <- "cell_id"
  lines2 <- capture_spatialkit_log(
    third <- assign_features_to_polygons(named, .two_cells(),
                                         polygon_id_col = "cell_id")
  )
  expect_true(log_has(lines2, "already carries column\\(s\\) 'cell_id'"))
  expect_equal(nrow(third), 2L)
  expect_equal(third$cell_id, c(1L, 1L))
})


test_that("assign_features_to_polygons invents ids when the layer has none", {
  bare <- sf::st_sf(geometry = sf::st_sfc(.sq(0, 0, 10, 10),
                                          .sq(10, 0, 20, 10), crs = 32632))
  out <- assign_features_to_polygons(.three_points(), bare)

  expect_true("poly_id" %in% names(out))
  expect_equal(out$poly_id, c(1L, 1L))         # both surviving points in cell 1
  expect_true(all(out$poly_id %in% seq_len(nrow(bare))))
})


test_that("assign_features_to_polygons harmonises CRS and restores the input's", {
  cells  <- .two_cells()                       # EPSG:32632
  pts_ll <- sf::st_transform(.three_points(), 4326)

  out <- assign_features_to_polygons(pts_ll, cells)

  # Result comes back in the CRS the features arrived in, not the polygons'.
  expect_equal(sf::st_crs(out), sf::st_crs(4326))
  # The assignment itself is the same one the projected input gets.
  ref <- assign_features_to_polygons(.three_points(), cells)
  expect_equal(out$poly_id, ref$poly_id)
  expect_equal(out$v, ref$v)
})


test_that("assign_features_to_polygons rejects unusable input", {
  expect_error(assign_features_to_polygons(sf::st_drop_geometry(.three_points()),
                                           .two_cells()),
               "`features_sf` must be an sf object")
  expect_error(assign_features_to_polygons(.three_points(), list()),
               "`polygons_sf` must be sf/sfc")
  expect_error(assign_features_to_polygons(.three_points(),
                                           .two_cells()[0, , drop = FALSE]),
               "has zero rows")

  # An sfc polygon layer is accepted and given sequential ids.
  sfc_cells <- sf::st_geometry(.two_cells())
  out <- assign_features_to_polygons(.three_points(), sfc_cells)
  expect_equal(nrow(out), 2L)
  expect_true("poly_id" %in% names(out))
})


test_that("summarize_by_cell consumes the ids assign_features_to_polygons emits", {
  # The contract that makes this function load-bearing: its output is the
  # input to the aggregation step, keyed on poly_id.
  set.seed(31)
  n <- 60
  pts <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 20), y = runif(n, 0, 10), val = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
  assigned <- assign_features_to_polygons(pts, .two_cells())
  cells    <- summarize_by_cell(assigned, response_var = "val")

  expect_equal(sum(cells$n), nrow(assigned))
  expect_equal(sort(cells$poly_id), sort(unique(assigned$poly_id)))
  # Per-cell means match a direct computation on the assigned points.
  for (i in seq_len(nrow(cells))) {
    id <- cells$poly_id[i]
    expect_equal(cells[["resp_mean_val"]][i],
                 mean(assigned$val[assigned$poly_id == id]))
  }
})
