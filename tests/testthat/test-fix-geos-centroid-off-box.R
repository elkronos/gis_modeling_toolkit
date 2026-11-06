# tests/testthat/test-fix-geos-centroid-off-box.R
# ---------------------------------------------------------------------------
# A centroid GEOS puts outside its feature is refused.
#
# GEOS, which takes the centroids of coerce_to_points() in a projected CRS,
# without a CRS, and in degrees with s2 off, weights each part of a ring by
# its signed area.  Where the edges of a ring cross (a bowtie whose lobes
# differ a little in area), or a hole lies outside its shell, the point lay
# far outside the feature, and was returned without a word: a bowtie 10 km
# across, one corner 1 m off symmetric, at x = 33,850,000.  The transform the
# s2 refusal of the same bowtie in degrees names led to that point.  Such a
# centroid, outside the bounding box of a feature GEOS reports invalid, now
# stops the call with the error of class `spatialkit_features_degenerate`,
# naming its rows (review finding R-fb180-1 of batch found-bug-180).
# ---------------------------------------------------------------------------

.go_quiet <- function(expr) {
  r <- NULL
  capture_spatialkit_log(r <- suppressWarnings(suppressMessages(expr)))
  r
}
# A bowtie `d` across at (x0, y0); `e` moves one corner up, off symmetric.
.go_bow <- function(x0, y0, d, e = 0)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + d, y0 + d + e), c(x0 + d, y0), c(x0, y0 + d),
                            c(x0, y0))))
.go_sq <- function(x0, y0, d)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + d, y0), c(x0 + d, y0 + d), c(x0, y0 + d),
                            c(x0, y0))))
# A ring that touches itself at one vertex, two lobes of the same sense.
.go_fig8 <- function(x0, y0, d) {
  m <- d / 2
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + m, y0 + m), c(x0 + d, y0), c(x0 + d, y0 + d),
                            c(x0 + m, y0 + m), c(x0, y0 + d), c(x0, y0))))
}
.go_hole_out <- function(x0, y0, d)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + d, y0), c(x0 + d, y0 + d), c(x0, y0 + d),
                            c(x0, y0)),
                      rbind(c(x0 + 2 * d, y0), c(x0 + 2 * d, y0 + d / 2),
                            c(x0 + 3 * d, y0 + d / 2), c(x0 + 3 * d, y0), c(x0 + 2 * d, y0))))
# Eight rows in EPSG:32615: squares, with row 2 a bowtie 1 m off symmetric
# and row 4 a collection holding one.
.go_layer <- function(crs = 32615) {
  x0 <- 5.1e5; y0 <- 4.3e6; d <- 1e4
  g <- lapply(0:7, function(i) .go_sq(x0 + i * d, y0, d))
  g[[2]] <- .go_bow(x0 + d, y0, d, 1)
  g[[4]] <- sf::st_geometrycollection(list(.go_bow(x0 + 3 * d, y0, d, 1)))
  sf::st_sf(v = c(1, 3, 2, 5, 4, 6, 8, 7), x1 = c(2, 1, 4, 3, 6, 5, 8, 7),
            geometry = sf::st_sfc(g, crs = crs))
}

test_that("a bowtie off symmetric got a point far outside it, and is now refused", {
  x <- .go_layer()
  # GEOS's own centroid of row 2 lies far east of the layer.
  far <- sf::st_coordinates(suppressWarnings(sf::st_centroid(sf::st_geometry(x)[2])))
  expect_gt(far[1L, "X"], 3e7)
  e <- tryCatch(coerce_to_points(x, "centroid"), error = function(e) e)
  expect_s3_class(e, "spatialkit_features_degenerate")
  expect_identical(e$rows, c(2L, 4L))
  expect_identical(e$n, 8L)
  expect_identical(e$site, "geos_centroid")
  expect_null(e$parent)
  expect_match(e$reasons, "^Self-intersection\\[")
  msg <- conditionMessage(e)
  expect_match(msg, paste0("^coerce_to_points\\(\\): the centroids GEOS takes of 2 of the 8 ",
                           "features of `x` lie outside each feature's bounding box \\(row 2: ",
                           "\"Self-intersection\\[.*\"; row 4: \"Self-intersection\\["))
  expect_match(msg, "Pass mode = \"point_on_surface\"", fixed = TRUE)
  expect_match(msg, "repair the features with sf::st_make_valid() and call coerce_to_points() again.",
               fixed = TRUE)
  # Under "auto" a POLYGON gets a point on its surface; the collection takes
  # a centroid, and is refused alone.
  e <- tryCatch(coerce_to_points(x, "auto"), error = function(e) e)
  expect_s3_class(e, "spatialkit_features_degenerate")
  expect_identical(e$rows, 4L)
  expect_match(conditionMessage(e), paste0("^coerce_to_points\\(\\): the centroid GEOS takes ",
                                           "of 1 of the 8 features of `x` lies outside the ",
                                           "feature's bounding box \\(row 4: "))
  # Under "line_midpoint" the rows that are not lines take centroids.
  xl <- x[1:2, ]
  sf::st_geometry(xl)[1] <- sf::st_linestring(rbind(c(5.1e5, 4.3e6), c(5.2e5, 4.3e6)))
  e <- tryCatch(coerce_to_points(xl, "line_midpoint"), error = function(e) e)
  expect_s3_class(e, "spatialkit_features_degenerate")
  expect_identical(e$rows, 2L)
  # The remedies run, and give points inside the features.
  for (y in list(coerce_to_points(x, "point_on_surface"),
                 coerce_to_points(sf::st_make_valid(x), "centroid"))) {
    xy <- sf::st_coordinates(y)
    expect_true(all(xy[, "X"] >= 5.1e5 & xy[, "X"] <= 5.9e5))
  }
})

test_that("a hole outside its shell is refused, and the same holds without a CRS and with s2 off", {
  x <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(.go_hole_out(0, 0, 10), .go_sq(0, 0, 1),
                                                 crs = 32615))
  e <- tryCatch(coerce_to_points(x, "centroid"), error = function(e) e)
  expect_s3_class(e, "spatialkit_features_degenerate")
  expect_identical(e$rows, 1L)
  expect_match(e$reasons, "^Hole lies outside shell")
  nocrs <- sf::st_set_crs(.go_layer(), NA)
  e <- tryCatch(coerce_to_points(nocrs, "centroid"), error = function(e) e)
  expect_identical(e$rows, c(2L, 4L))
  # In degrees with s2 off, GEOS takes the centroids, planar in degrees.
  old <- suppressMessages(sf::sf_use_s2(FALSE))
  withr::defer(suppressMessages(sf::sf_use_s2(old)))
  ll <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(.go_bow(-94, 39, 0.1, 1e-4),
                                                  .go_sq(-93.9, 39, 0.1), crs = 4326))
  e <- tryCatch(.go_quiet(coerce_to_points(ll, "centroid")), error = function(e) e)
  expect_s3_class(e, "spatialkit_features_degenerate")
  expect_identical(e$rows, 1L)
  expect_identical(e$site, "geos_centroid")
})

test_that("centroids inside the box are GEOS's, to the bit", {
  # Valid squares, a C shape whose centroid lies outside it but inside its
  # box, a symmetric bowtie (its centre), a figure eight that touches
  # itself, and two squares that overlap: GEOS reports the last three
  # invalid, and their points lie in their boxes.
  cshape <- sf::st_polygon(list(rbind(c(0, 0), c(3, 0), c(3, 1), c(1, 1), c(1, 2), c(3, 2),
                                      c(3, 3), c(0, 3), c(0, 0))))
  over <- sf::st_multipolygon(list(unclass(.go_sq(0, 0, 10)), unclass(.go_sq(5, 5, 10))))
  g <- sf::st_sfc(.go_sq(0, 0, 1), cshape, .go_bow(5.1e5, 4.3e6, 1e4), .go_fig8(0, 0, 10),
                  over, sf::st_geometrycollection(list(.go_fig8(20, 0, 10))), crs = 32615)
  x <- sf::st_sf(id = seq_along(g), geometry = g)
  ref <- suppressWarnings(sf::st_centroid(g))
  for (m in c("centroid", "auto")) {
    got <- coerce_to_points(x, m)
    expect_identical(sf::st_geometry(got)[6], ref[6], label = m)
    if (m == "centroid") expect_identical(sf::st_geometry(got), ref)
  }
  expect_equal(unname(sf::st_coordinates(ref)[3, ]), c(515000, 4305000), tolerance = 1e-12)
})

test_that("prep_model_data() and compare_models_cv() refuse it in their words", {
  skip_if_not_installed("ranger")
  x <- .go_layer()
  for (pz in c("centroid", "auto")) {
    e <- tryCatch(.go_quiet(prep_model_data(x, "v", "x1", pointize = pz)),
                  error = function(e) e)
    expect_s3_class(e, "spatialkit_features_degenerate")
    expect_identical(e$site, "geos_centroid")
    expect_identical(e$rows, if (pz == "centroid") c(2L, 4L) else 4L)
    msg <- conditionMessage(e)
    expect_match(msg, "^prep_model_data\\(\\): the centroids? GEOS takes of [12] of the 8 features of `data_sf`")
    expect_match(msg, paste0("Pass pointize = \"point_on_surface\", which takes a point on each ",
                             "geometry with GEOS, or repair the features with sf::st_make_valid() ",
                             "and run the call again."), fixed = TRUE)
    for (no in c("`x`", "mode =", "coerce_to_points() again", "Transform"))
      expect_false(grepl(no, msg, fixed = TRUE), label = no)
  }
  e <- tryCatch(.go_quiet(compare_models_cv(x, "v", "x1", models = "RF", k = 3,
                                            pointize = "centroid", quiet = TRUE)),
                error = function(e) e)
  expect_s3_class(e, "spatialkit_features_degenerate")
  expect_identical(e$rows, c(2L, 4L))
  expect_match(conditionMessage(e), "^compare_models_cv\\(\\): the centroids GEOS takes of 2 of the 8 features of `data_sf`")
  # A point on the surface runs.
  expect_no_error(.go_quiet(prep_model_data(x, "v", "x1", pointize = "point_on_surface")))
})
