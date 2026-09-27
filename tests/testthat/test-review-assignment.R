# ===========================================================================
# assign_features_to_polygons(largest = TRUE): regressions from the review.
#
# The largest-overlap join used to sit in a catch-all tryCatch that retried
# WITHOUT `largest` on any error.  No predicate ever rejects `largest` (sf
# ignores `join` on that path), so the retry only fired on geometry failures,
# and then every straddling feature in the layer was assigned by `tie_break`
# instead of by overlap, with no R warning.
# ===========================================================================

.rv_sq <- function(x0, y0, x1, y1) {
  sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1),
                            c(x0, y0))))
}

# Two 10 x 10 cells side by side, equal areas, so a tie-break on area alone
# falls to the first cell.
.rv_two_cells <- function(second = .rv_sq(10, 0, 20, 10)) {
  sf::st_sf(poly_id = c(1L, 2L),
            geometry = sf::st_sfc(.rv_sq(0, 0, 10, 10), second, crs = 32632))
}

# A self-crossing ring: GEOS cannot intersect it (TopologyException).
.rv_bowtie <- sf::st_polygon(list(rbind(c(1, 1), c(3, 3), c(3, 1), c(1, 3),
                                        c(1, 1))))

# Collect every R warning `expr` raises while returning its value.
.rv_warnings <- function(expr) {
  ws <- character(0)
  val <- withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = val, warnings = ws)
}


test_that("one invalid feature does not switch the layer to tie-break assignment", {
  # x 7..14: 3 units of width in cell 1, 4 in cell 2.
  straddler <- sf::st_sfc(.rv_sq(7, 2, 14, 5), crs = 32632)
  x <- sf::st_sf(v = 1:2,
                 geometry = c(straddler, sf::st_sfc(.rv_bowtie, crs = 32632)))
  expect_false(all(sf::st_is_valid(x)))

  lines <- capture_spatialkit_log(
    res <- .rv_warnings(assign_features_to_polygons(x, .rv_two_cells())))
  out <- res$value

  # The valid straddler keeps its largest-overlap cell.  It used to go to
  # cell 1: the bow-tie made GEOS throw, the whole join was rerun without
  # `largest`, and the equal-area tie-break took the first match.
  expect_equal(out$poly_id[out$v == 1L], 2L)
  # The bow-tie lies wholly in cell 1 either way.
  expect_equal(out$poly_id[out$v == 2L], 1L)
  expect_equal(attr(out, "ties")$n, 0L)
  # The repair is announced as a real R warning, with the count.
  expect_true(any(grepl("1 of 2 feature\\(s\\) and 0 of 2 polygon\\(s\\) have invalid geometry",
                        res$warnings)))
  expect_true(log_has(lines, "repaired with sf::st_make_valid\\(\\)"))
  # The repair is for the join only: the caller's geometry comes back as it
  # arrived.
  expect_identical(sf::st_geometry(out)[[2]], .rv_bowtie)
  expect_false(sf::st_is_valid(sf::st_geometry(out)[2]))
})


test_that("an invalid cell is repaired for the join instead of dropping `largest`", {
  # Cell 2 is a bow-tie with lobes left and right of (15, 5); its shoelace
  # area cancels to 0, so the old fallback's smallest-area tie-break always
  # picked it.  Repaired, the straddler x 6..13 has 12 units in cell 1 and
  # 8.5 in cell 2's left lobe.
  bow_cell <- sf::st_polygon(list(rbind(c(10, 0), c(20, 10), c(20, 0),
                                        c(10, 10), c(10, 0))))
  cells <- .rv_two_cells(second = bow_cell)
  x <- sf::st_sf(v = 1L, geometry = sf::st_sfc(.rv_sq(6, 2, 13, 5), crs = 32632))

  lines <- capture_spatialkit_log(
    res <- .rv_warnings(assign_features_to_polygons(x, cells)))

  expect_equal(res$value$poly_id, 1L)
  expect_true(any(grepl("0 of 1 feature\\(s\\) and 1 of 2 polygon\\(s\\) have invalid geometry",
                        res$warnings)))
})


test_that("a largest-overlap join that still fails stops instead of falling back", {
  # After repair nothing in sf's largest path should throw, but when it does
  # (s2 on a degenerate intersection piece, say) the answer must not be
  # recomputed under a different rule for the whole layer.
  real_st_join <- sf::st_join
  local_mocked_bindings(
    st_join = function(x, y, ..., largest = FALSE) {
      if (isTRUE(largest)) stop("TopologyException: simulated")
      real_st_join(x, y, ..., largest = largest)
    },
    .package = "sf"
  )
  x <- sf::st_sf(v = 1L, geometry = sf::st_sfc(.rv_sq(7, 2, 14, 5), crs = 32632))

  expect_error(suppressWarnings(assign_features_to_polygons(x, .rv_two_cells())),
               "largest-overlap join failed \\(TopologyException: simulated\\).*largest = FALSE")
  # Asking for the tie-break rule explicitly still works.
  expect_equal(suppressWarnings(
    assign_features_to_polygons(x, .rv_two_cells(), largest = FALSE))$poly_id, 1L)
})


# Evaluate `code` with sf_use_s2() set to `on`, restoring the session's value.
.rv_with_s2 <- function(on, code) {
  was <- suppressMessages(sf::sf_use_s2(on))
  on.exit(suppressMessages(sf::sf_use_s2(was)), add = TRUE)
  force(code)
}

test_that("lon/lat features are assigned in the cells' projected CRS", {
  # Two 200 km cells in Web Mercator sharing a horizontal edge near 45N.  The
  # features used to stay in lon/lat and the cells were moved to them, where
  # s2 reads that edge as a great-circle arc bulging about 0.8 km north
  # mid-way; with s2 off the largest-overlap join needed lwgeom and failed.
  y0 <- 5621521                          # ~45N in EPSG:3857
  w  <- 2e5 / cos(pi / 4)                # 200 km of ground
  mid <- w / 2
  cells <- sf::st_sf(poly_id = c(1L, 2L),
                     geometry = sf::st_sfc(.rv_sq(0, y0 - w, w, y0),
                                           .rv_sq(0, y0, w, y0 + w),
                                           crs = 3857))
  # 1200 of its 2000 units of height lie in cell 2, drawn in the cells' CRS.
  feat_ll <- sf::st_transform(
    sf::st_sf(v = 1L, geometry = sf::st_sfc(
      .rv_sq(mid - 5000, y0 - 800, mid + 5000, y0 + 1200), crs = 3857)),
    4326)
  # Points 500 units either side of the edge.
  pts_ll <- sf::st_transform(
    sf::st_sf(v = 1:2, geometry = sf::st_sfc(sf::st_point(c(mid, y0 + 500)),
                                             sf::st_point(c(mid, y0 - 500)),
                                             crs = 3857)),
    4326)

  for (s2 in c(TRUE, FALSE)) {
    .rv_with_s2(s2, {
      out <- suppressWarnings(assign_features_to_polygons(feat_ll, cells))
      pts <- assign_features_to_polygons(pts_ll, cells)
    })
    # Under s2 the old code put the polygon, and the point inside cell 2, in
    # cell 1.
    expect_equal(out$poly_id, 2L, info = paste("s2 =", s2))
    expect_equal(pts$poly_id, c(2L, 1L), info = paste("s2 =", s2))
    # Only a copy was moved into the cells' CRS: the caller's coordinates come
    # back untouched, not as a transform round trip.
    expect_identical(sf::st_geometry(out), sf::st_geometry(feat_ll))
    expect_identical(sf::st_geometry(pts), sf::st_geometry(pts_ll))
  }
})
