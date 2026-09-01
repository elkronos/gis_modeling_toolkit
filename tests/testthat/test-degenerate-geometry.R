# ===========================================================================
# Degenerate geometry: zero rows, a single point, and all-coincident points.
#
# Roughly a dozen explicit guards in R/ exist for exactly these inputs and
# fewer than half were exercised.  Degenerate geometry is the most common way
# real sf pipelines break, and the failure mode that matters is not "it
# errored" but "it returned something plausible and wrong" -- an empty fold, a
# zero-area clip target, a bounding box with no extent.
#
# Each block therefore pins WHICH answer comes back, not merely that a call
# survives.
# ===========================================================================

.dg_empty <- function() {
  sf::st_sf(z = numeric(0), w = numeric(0), geometry = sf::st_sfc(crs = 32632))
}

.dg_one <- function() {
  sf::st_sf(z = 1, w = 2,
            geometry = sf::st_sfc(sf::st_point(c(5, 5)), crs = 32632))
}

.dg_coincident <- function(n = 10, seed = 1) {
  set.seed(seed)
  sf::st_as_sf(data.frame(x = rep(5, n), y = rep(5, n),
                          z = rnorm(n), w = rnorm(n)),
               coords = c("x", "y"), crs = 32632)
}


# ---------------------------------------------------------------------------
# .assert_sf(): the shared gate in front of most of the others
# ---------------------------------------------------------------------------

test_that(".assert_sf refuses a zero-row layer", {
  # all() is vacuously TRUE on the empty set, so without the explicit length
  # check a zero-row layer would sail through every geometry-type gate in the
  # package.
  expect_error(spatialkit:::.assert_sf(.dg_empty(), c("POINT", "MULTIPOINT")),
               "geometry must be one of: POINT, MULTIPOINT")
  expect_silent(spatialkit:::.assert_sf(.dg_one(), "POINT"))
})


# ---------------------------------------------------------------------------
# make_folds()
# ---------------------------------------------------------------------------

test_that("make_folds(random_kfold) degrades rather than inventing folds", {
  # Zero rows are rejected outright, for every method.  Returning k = 0 with no
  # folds would sail into CV and surface only as a generic "all folds failed"
  # at the very end; block_kfold did not even get that far, dying inside
  # sf::st_as_sfc() on an all-NA bbox with "!anyNA(x) is not TRUE".
  expect_error(
    make_folds(.dg_empty(), k = 3, method = "random_kfold", seed = 1),
    "has no rows")
  expect_error(
    make_folds(.dg_empty(), k = 3, method = "block_kfold", seed = 1),
    "has no rows")

  # One row: exactly one leave-one-out fold, whose training set is (correctly)
  # empty -- there is nothing else to train on.
  f1 <- suppressWarnings(make_folds(.dg_one(), k = 3, method = "random_kfold",
                                    seed = 1))
  expect_equal(f1$k, 1)
  expect_length(f1$folds, 1L)
  expect_equal(f1$folds[[1]]$test, 1L)
  expect_length(f1$folds[[1]]$train, 0L)

  # Coincident points are still n distinct OBSERVATIONS; random k-fold does not
  # look at geometry, so all 10 must be partitioned across the 3 folds.
  fd <- make_folds(.dg_coincident(), k = 3, method = "random_kfold", seed = 1)
  expect_equal(fd$k, 3)
  sizes <- vapply(fd$folds, function(z) length(z$test), integer(1))
  expect_equal(sum(sizes), 10L)
  expect_equal(max(sizes) - min(sizes), 1L)
  expect_equal(sort(unlist(lapply(fd$folds, `[[`, "test"))), 1:10)
})


test_that("make_folds(block_kfold) refuses a collapsed extent", {
  # Coincident points have a zero-area bounding box, so every block grid
  # collapses to a single block -- one fold with an empty training set.  That
  # has to be an error: a CV run on it would report a score for a model that
  # was never fitted.
  expect_error(
    suppressMessages(make_folds(.dg_coincident(), k = 3,
                                method = "block_kfold", seed = 1)),
    "single block covering the whole extent")
  expect_error(
    suppressMessages(make_folds(.dg_one(), k = 3, method = "block_kfold",
                                seed = 1)),
    "single block covering the whole extent")

  # Zero rows must also fail rather than return empty folds.  Deliberately no
  # regexp here: the current failure comes from sf's bbox machinery
  # ("!anyNA(x) is not TRUE" out of st_as_sfc() on an all-NA bbox) rather than
  # from a named guard in make_folds(), so pinning the wording would pin a
  # message that ought to be improved.  What is asserted is the contract --
  # the call must not return -- which survives that improvement.
  expect_error(suppressMessages(make_folds(.dg_empty(), k = 3,
                                           method = "block_kfold", seed = 1)))
})


# ---------------------------------------------------------------------------
# coerce_to_points()
# ---------------------------------------------------------------------------

test_that("coerce_to_points passes degenerate input through unchanged", {
  empty <- coerce_to_points(.dg_empty(), "auto")
  expect_s3_class(empty, "sf")
  expect_equal(nrow(empty), 0L)
  expect_equal(sf::st_crs(empty), sf::st_crs(32632))

  one <- coerce_to_points(.dg_one(), "auto")
  expect_equal(nrow(one), 1L)
  expect_equal(unname(sf::st_coordinates(one)[1, 1:2]), c(5, 5))

  # Coincident points stay 10 rows -- coercion is per feature, never a dedup.
  dup <- coerce_to_points(.dg_coincident(), "auto")
  expect_equal(nrow(dup), 10L)
  expect_equal(unname(sf::st_coordinates(dup)[, 1]), rep(5, 10))
  # Attributes ride along untouched.
  expect_equal(dup$z, .dg_coincident()$z)

  # Every mode agrees on a POINT layer, which is what "representative point"
  # has to mean for a point.
  for (m in c("centroid", "point_on_surface", "bbox_center")) {
    got <- suppressWarnings(coerce_to_points(.dg_coincident(), m))
    expect_equal(nrow(got), 10L, info = m)
    expect_equal(unname(sf::st_coordinates(got)[, 1:2]),
                 unname(sf::st_coordinates(.dg_coincident())[, 1:2]), info = m)
  }
})


# ---------------------------------------------------------------------------
# create_voronoi_polygons()
# ---------------------------------------------------------------------------

test_that("create_voronoi_polygons handles one point and coincident points", {
  expect_error(create_voronoi_polygons(.dg_empty(), quiet = TRUE),
               "geometry must be one of: POINT, MULTIPOINT")

  # A single point tessellates the whole (buffered) envelope as one cell.
  one <- suppressMessages(create_voronoi_polygons(.dg_one(), quiet = TRUE))
  expect_equal(nrow(one$cells), 1L)
  expect_equal(one$index, 1L)
  expect_gt(as.numeric(sf::st_area(one$cells)), 0)

  # Ten coincident points are ONE Voronoi generator, so one cell -- but ten
  # index entries, one per input feature.  A cells/index length mismatch here
  # is what silently misaligns an attribute join downstream.
  dup <- suppressMessages(create_voronoi_polygons(.dg_coincident(),
                                                  quiet = TRUE))
  expect_equal(nrow(dup$cells), 1L)
  expect_length(dup$index, 10L)
  expect_equal(dup$index, rep(1L, 10L))

  # keep_duplicates does not manufacture cells out of coincident generators.
  keep <- suppressMessages(create_voronoi_polygons(.dg_coincident(),
                                                   keep_duplicates = TRUE,
                                                   quiet = TRUE))
  expect_equal(nrow(keep$cells), 1L)
  expect_length(keep$index, 10L)
})


# ---------------------------------------------------------------------------
# clip_target_for()
# ---------------------------------------------------------------------------

test_that("clip_target_for buffers a zero-extent point set instead of returning nothing", {
  expect_error(clip_target_for(.dg_empty(), quiet = TRUE),
               "geometry must be one of: POINT, MULTIPOINT")

  for (input in list(.dg_one(), .dg_coincident())) {
    tgt <- suppressMessages(clip_target_for(input))
    expect_s3_class(tgt, "sf")
    expect_equal(nrow(tgt), 1L)
    # A degenerate bbox is not a polygon; the fallback must produce real area
    # and must still cover the points.
    expect_gt(as.numeric(sf::st_area(tgt)), 0)
    expect_true(all(lengths(sf::st_covered_by(input, tgt)) > 0L))
  }
})


# ---------------------------------------------------------------------------
# ensure_stable_poly_id()
# ---------------------------------------------------------------------------

test_that("ensure_stable_poly_id refuses an empty layer and numbers a single polygon", {
  expect_error(
    ensure_stable_poly_id(sf::st_sf(geometry = sf::st_sfc(crs = 32632))),
    "no polygon rows found")

  one_poly <- sf::st_sf(v = 1L, geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)))),
    crs = 32632))
  out <- ensure_stable_poly_id(one_poly, transform_for_sort = NULL)
  expect_equal(nrow(out), 1L)
  expect_equal(out$poly_id, 1L)

  # Two IDENTICAL polygons are two features and get two distinct IDs; the
  # tie-breakers (area, then original index) keep the order deterministic.
  two_same <- rbind(one_poly, one_poly)
  two_same$v <- 1:2
  out2 <- ensure_stable_poly_id(two_same, transform_for_sort = NULL)
  expect_equal(out2$poly_id, 1:2)
  expect_equal(out2$v, 1:2)
})


# ---------------------------------------------------------------------------
# summarize_by_cell()
# ---------------------------------------------------------------------------

test_that("summarize_by_cell returns the full column set for degenerate input", {
  empty <- sf::st_sf(poly_id = integer(0), y = numeric(0),
                     geometry = sf::st_sfc(crs = 32632))
  out <- summarize_by_cell(empty, response_var = "y")
  expect_equal(nrow(out), 0L)
  # The columns must still be there, or rbind()ing per-region results ragged.
  expect_true(all(c("poly_id", "n", "resp_mean_y", "..sd_resp_y",
                    "..se_resp_y", "cell_weight") %in% names(out)))

  # A single observation in a cell: sd and se are undefined (NA), not 0 --
  # reporting 0 would claim perfect precision from one measurement.
  one <- sf::st_sf(poly_id = 1L, y = 3,
                   geometry = sf::st_sfc(sf::st_point(c(1, 1)), crs = 32632))
  o1 <- summarize_by_cell(one, response_var = "y")
  expect_equal(nrow(o1), 1L)
  expect_equal(o1$n, 1L)
  expect_equal(o1[["resp_mean_y"]], 3)
  expect_true(is.na(o1[["..sd_resp_y"]]))
  expect_true(is.na(o1[["..se_resp_y"]]))
  expect_equal(o1$cell_weight, 1L)

  # All observations in ONE cell at ONE location: still a single summarised
  # row over 10 values, with a real sd.
  dup <- .dg_coincident()
  dup$poly_id <- 1L
  od <- summarize_by_cell(dup, response_var = "z")
  expect_equal(nrow(od), 1L)
  expect_equal(od$n, 10L)
  expect_equal(od[["resp_mean_z"]], mean(dup$z))
  expect_equal(od[["..se_resp_z"]], stats::sd(dup$z) / sqrt(10))

  expect_error(summarize_by_cell(.dg_one(), response_var = "z"),
               "could not find an ID column")
})


# ---------------------------------------------------------------------------
# predict_surface()
# ---------------------------------------------------------------------------

test_that("predict_surface refuses to build a grid over no extent", {
  # A single training point (or ten coincident ones) has a zero-area bounding
  # box, so there is no grid to derive.  Silently returning a one-cell surface
  # would look like a prediction map.
  for (input in list(.dg_one(), .dg_coincident())) {
    fit <- lm_spatial_fit(input, predictor_vars = character(0))
    expect_error(predict_surface(fit, n_cells = 25),
                 "no finite extent")
  }

  fit <- lm_spatial_fit(surf_test_points(50), predictor_vars = character(0))

  # An empty user grid: named guard rather than seq(1L, 0L)'s "wrong sign in
  # 'by' argument".
  expect_error(
    predict_surface(fit, grid = sf::st_sf(geometry = sf::st_sfc(crs = 3857))),
    "`grid` has no rows")

  # A boundary that excludes the whole grid.
  far <- sf::st_sfc(sf::st_polygon(list(rbind(
    c(1e6, 1e6), c(1e6 + 1, 1e6), c(1e6 + 1, 1e6 + 1), c(1e6, 1e6 + 1),
    c(1e6, 1e6)))), crs = 3857)
  expect_error(predict_surface(fit, boundary = far),
               "no grid points fall inside `boundary`")

  # A fit whose training data lost its geometry.
  no_geom <- fit
  no_geom$data_sf <- sf::st_drop_geometry(no_geom$data_sf)
  expect_error(predict_surface(no_geom), "no training geometry")

  # The non-degenerate baseline, so the errors above are not masking a
  # function that never works.
  ok <- predict_surface(fit, n_cells = 25)
  expect_s3_class(ok, "sf")
  expect_gt(nrow(ok), 0L)
  expect_true(".pred" %in% names(ok))
  expect_false(anyNA(ok$.pred))
})


# ---------------------------------------------------------------------------
# An EMPTY geometry INSIDE a populated layer
#
# The blocks above cover zero-row layers, which every guard already handled.
# An EMPTY POINT sitting among ordinary ones is a different animal: nrow() is
# positive, st_geometry_type() says POINT, and st_coordinates() returns one
# ALL-NA row for it rather than none -- so it slips past every row-count and
# type check and reaches the distance code.  In block_kfold that meant
# st_intersects() -> integer(0) -> ..block_id NA -> an all-NA distance row ->
# which.min() returning integer(0) -> "replacement has length zero".
# ---------------------------------------------------------------------------

.dg_with_empty <- function(n = 60, seed = 9) {
  set.seed(seed)
  geom <- sf::st_sfc(
    c(lapply(seq_len(n),
             function(i) sf::st_point(c(runif(1, 0, 1000), runif(1, 0, 1000)))),
      list(sf::st_point())),
    crs = 32632)
  sf::st_sf(z = c(rnorm(n), 0), w = c(rnorm(n), 0), geometry = geom)
}

test_that("make_folds drops an EMPTY geometry from a populated layer", {
  pts <- .dg_with_empty()
  n   <- nrow(pts) - 1L
  expect_equal(sum(sf::st_is_empty(pts)), 1L)
  expect_true(all(as.character(sf::st_geometry_type(pts)) == "POINT"))
  # st_coordinates() returns a FULL-HEIGHT matrix with one all-NA row, not a
  # shorter one -- which is exactly why a row-count check cannot spot the
  # empty geometry and a complete.cases() check has to.
  xy <- suppressWarnings(sf::st_coordinates(pts))
  expect_equal(nrow(xy), n + 1L)
  expect_equal(sum(!stats::complete.cases(xy[, 1:2, drop = FALSE])), 1L)

  grid <- sf::st_as_sf(
    data.frame(x = stats::runif(15, 0, 1000), y = stats::runif(15, 0, 1000)),
    coords = c("x", "y"), crs = 32632)

  args <- list(
    random_kfold = list(k = 3, method = "random_kfold", seed = 1),
    block_kfold  = list(k = 3, method = "block_kfold",  seed = 1),
    buffered_loo = list(k = 3, method = "buffered_loo", buffer = 100, seed = 1),
    nndm         = list(method = "nndm", prediction_points = grid, seed = 1)
  )

  for (m in names(args)) {
    lines <- capture_spatialkit_log(
      f <- do.call(make_folds, c(list(points_sf = pts), args[[m]])))
    # Dropped, and the count is named rather than silently absorbed.
    expect_true(log_has(lines, "dropping 1 point\\(s\\) with empty or non-finite"),
                info = m)
    # The empty row appears in no fold and in no assignment row ...
    expect_equal(nrow(f$assignment), n, info = m)
    expect_false((n + 1L) %in% f$assignment$row_id, info = m)
    expect_setequal(f$assignment$row_id, seq_len(n))
    # ... and the survivors keep their ORIGINAL row identities, so the folds
    # still index the layer the caller passed in.
    all_rows <- sort(unique(unlist(lapply(f$folds, `[[`, "test"))))
    expect_true(all(all_rows <= n), info = m)
    expect_gt(length(f$folds), 0L)
    expect_true(all(vapply(f$folds, function(s) length(s$train) > 0L,
                           logical(1))), info = m)
  }
})

test_that("make_folds refuses a layer whose geometries are ALL empty", {
  allempty <- sf::st_sf(
    z = c(1, 2, 3), w = c(4, 5, 6),
    geometry = sf::st_sfc(sf::st_point(), sf::st_point(), sf::st_point(),
                          crs = 32632))
  expect_equal(nrow(allempty), 3L)
  expect_error(suppressWarnings(make_folds(allempty, k = 2,
                                           method = "random_kfold", seed = 1)),
               "no usable coordinates; there is nothing to split into folds")
})
