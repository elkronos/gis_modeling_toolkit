# ===========================================================================
# Regressions from the second review of assignment, aggregation, stable IDs
# and the spatial half-split.  Every test here failed on the code before the
# fix it names.
# ===========================================================================

# Collect every R warning `expr` raises (message and class) while returning
# its value.
.r2_warnings <- function(expr) {
  msgs <- character(0); classes <- list()
  val <- withCallingHandlers(expr, warning = function(w) {
    msgs    <<- c(msgs, conditionMessage(w))
    classes <<- c(classes, list(class(w)))
    invokeRestart("muffleWarning")
  })
  list(value = val, warnings = msgs, classes = classes)
}

.r2_sq <- function(x0, y0, s = 10) {
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s),
                            c(x0, y0 + s), c(x0, y0))))
}

# Four cells of five points, with a response.
.r2_points <- function() {
  set.seed(1)
  sf::st_as_sf(data.frame(x = runif(20, 0, 200), y = runif(20, 0, 200),
                          v = rnorm(20), poly_id = rep(1:4, each = 5)),
               coords = c("x", "y"), crs = 32632)
}

# Six cells of 25 points on a correlated Gaussian field (as in
# test-summarize-by-cell.R).
.r2_field <- function() {
  set.seed(77)
  n  <- 150
  x  <- rep(seq(0, 500, length.out = 6), each = 25) + runif(n, 0, 60)
  y  <- runif(n, 0, 60)
  d  <- as.matrix(stats::dist(cbind(x, y)))
  z  <- as.numeric(t(chol(exp(-d / 40) + diag(1e-6, n))) %*% rnorm(n))
  sf::st_as_sf(data.frame(x = x, y = y, z = z, poly_id = rep(1:6, each = 25)),
               coords = c("x", "y"), crs = 32632)
}

.r2_sac <- function(model = data.frame(model = c("Nug", "Exp"),
                                       psill = c(0.2, 0.8),
                                       range = c(0, 40))) {
  structure(120, class = c("sac_range", "numeric"), variogram_model = model,
            crs = sf::st_crs(32632))
}


# --- summarize_by_cell(): design effects ------------------------------------

test_that("a non-finite numeric deff falls back to 1 with the classed warning", {
  pts <- .r2_points()
  for (bad in list(NA_real_, NaN, Inf)) {
    # NA and NaN used to stop with "missing value where TRUE/FALSE needed";
    # Inf passed the check and gave uncorrected SEs beside cell_weight 0.
    expect_warning(
      out <- summarize_by_cell(pts, "v", deff = bad),
      "must be a single number >= 1", class = "spatialkit_deff_fallback")
    expect_equal(out$cell_weight, out$n)
    expect_null(attr(out, "deff_applied"))
    expect_identical(out$deff_applied, rep(FALSE, 4))
  }
})

test_that("deff_applied marks each row, and survives binding results together", {
  pts <- .r2_points()
  # The default frame is unchanged: no column.
  expect_false("deff_applied" %in% names(summarize_by_cell(pts, "v")))
  fixed <- summarize_by_cell(pts, "v", deff = 2)
  expect_identical(fixed$deff_applied, rep(TRUE, 4))
  expect_false(is.null(attr(fixed, "deff_applied")))
  fell <- suppressWarnings(summarize_by_cell(pts, "v", deff = 0.5))
  expect_identical(fell$deff_applied, rep(FALSE, 4))
  # Combined, the attribute is at best the first result's (bind_rows()) or
  # gone (rbind()); the column describes every row.
  expect_identical(dplyr::bind_rows(fixed, fell)$deff_applied,
                   rep(c(TRUE, FALSE), each = 4))
  expect_identical(rbind(fell, fixed)$deff_applied,
                   rep(c(FALSE, TRUE), each = 4))
})

test_that("a rejected sac falls back with the classed warning and FALSE rows", {
  pts <- .r2_field()
  rejected <- structure(NA_real_, class = c("sac_range", "numeric"),
                        rejected_reason = "fitted range exceeds the largest lag fitted",
                        variogram_model = data.frame(model = "Exp", psill = 1,
                                                     range = 1e6))
  expect_warning(
    out <- summarize_by_cell(pts, predictor_vars = "z", deff = "variogram",
                             sac = rejected),
    "fitted range exceeds", class = "spatialkit_deff_fallback")
  expect_identical(out$deff_applied, rep(FALSE, 6))
  # Caught by class, without matching the message.
  caught <- tryCatch(
    summarize_by_cell(pts, predictor_vars = "z", deff = "variogram",
                      sac = rejected),
    spatialkit_deff_fallback = function(w) "caught")
  expect_identical(caught, "caught")
  # Joined to cells, an empty cell carries NA rather than a verdict.
  cells <- sf::st_sf(poly_id = 1:7,
                     geometry = sf::st_sfc(lapply(0:6, function(i) .r2_sq(100 * i, 0, 60)),
                                           crs = 32632))
  ok <- summarize_by_cell(pts, "z", deff = "variogram", sac = .r2_sac(),
                          cells_sf = cells)
  expect_identical(ok$deff_applied, c(rep(TRUE, 6), NA))
})

test_that("a variogram fallback with no fit reports why, as an R warning", {
  pts <- .r2_field()
  skip_if_not_installed("gstat")
  # 25 points: estimate_sac_range() returns no model.  This fell back to
  # deff = 1 with only a log line.
  res <- .r2_warnings(summarize_by_cell(pts[1:25, ], "z", deff = "variogram"))
  hit <- vapply(res$classes, function(k) "spatialkit_deff_fallback" %in% k, logical(1))
  expect_true(any(hit))
  expect_match(res$warnings[hit][1], "requires a fitted variogram model")
  expect_false(any(res$value$deff_applied))
})

test_that("each cell's correlation matrix is built once, not once per column statistic", {
  pts <- .r2_field()
  sac <- .r2_sac()
  n_builds <- 0L
  real <- spatialkit:::.cor_stats_from_coords
  local_mocked_bindings(.cor_stats_from_coords = function(...) {
    n_builds <<- n_builds + 1L
    real(...)
  })
  summarize_by_cell(pts, "z", deff = "variogram", sac = sac)
  # One per cell; the design effect at the row count reuses them (it used to
  # rebuild all six).
  expect_identical(n_builds, 6L)

  n_builds <- 0L
  p2 <- pts
  p2$z[c(1, 30, 60)] <- NA         # one missing value in cells 1, 2 and 3
  p2$w <- p2$z + 1
  summarize_by_cell(p2, "z", predictor_vars = "w", deff = "variogram",
                    sac = sac, conf_level = 0.95)
  # 6 for the cells, then per column the three incomplete cells once each
  # (z and w), and cell_weight's three: 15.  The se/neff/df/ci closures used
  # to rebuild each of those five times: 45.
  expect_identical(n_builds, 15L)
})

test_that("an empty point no longer stops deff = 'variogram'", {
  pts <- .r2_field()
  g <- sf::st_geometry(pts)
  g[5] <- sf::st_point()
  empty <- pts
  sf::st_geometry(empty) <- g
  res <- .r2_warnings(summarize_by_cell(empty, "z", deff = "variogram",
                                        sac = .r2_sac()))
  expect_true(any(grepl("1 point\\(s\\) have empty or non-finite coordinates",
                        res$warnings)))
  out <- res$value
  expect_true(all(out$deff_applied))
  # Its value still counts in its cell; the other cells are untouched.
  expect_identical(out$n[1], 25L)
  expect_equal(out$resp_mean_z[1], mean(pts$z[1:25]))
  ref <- summarize_by_cell(pts, "z", deff = "variogram", sac = .r2_sac())
  expect_equal(out[["..se_resp_z"]][-1], ref[["..se_resp_z"]][-1])
  # The design effect of cell 1 is 1 + (25 - 1) * rbar over its 24 located
  # points.
  cor_fn <- spatialkit:::.vgm_correlation_fn(attr(.r2_sac(), "variogram_model"))
  xy <- sf::st_coordinates(pts)[setdiff(1:25, 5), 1:2]
  R  <- matrix(cor_fn(as.numeric(as.matrix(stats::dist(xy)))), 24, 24)
  diag(R) <- 1
  rbar <- (sum(R) - 24) / (24 * 23)
  expect_equal(attr(out, "deff_applied")$deff[1], 1 + 24 * rbar)
})

test_that("an anisotropic variogram is evaluated in its own geometry", {
  # vgm(0.8, "Exp", 300, 0.2, anis = c(0, 0.2)): the major axis runs north,
  # the east-west range is 60.  gstat's variogramLine() gives correlations of
  # 0.348 / 0.151 / 0.029 east-west at 50 / 100 / 200 m, and 0.677 / 0.573 /
  # 0.411 north-south.  Read as isotropic, every direction got the latter.
  m <- data.frame(model = c("Nug", "Exp"), psill = c(0.2, 0.8),
                  range = c(0, 300), ang1 = 0, anis1 = c(1, 0.2))
  fn <- spatialkit:::.vgm_correlation_fn(m)
  cor_xy <- attr(fn, "cor_xy")
  expect_true(is.function(cor_xy))
  xy <- rbind(c(0, 0), c(50, 0), c(100, 0), c(200, 0), c(0, 50), c(0, 200))
  R <- cor_xy(xy)
  expect_equal(unname(R[1, 2:4]), 0.8 * exp(-c(50, 100, 200) / 60))
  expect_equal(unname(R[1, 5:6]), 0.8 * exp(-c(50, 200) / 300))
  # Rotated: the major axis 30 degrees east of north.
  m2 <- transform(m, ang1 = 30, anis1 = c(1, 0.5))
  cx <- attr(spatialkit:::.vgm_correlation_fn(m2), "cor_xy")
  along  <- rbind(c(0, 0), 100 * c(sin(pi / 6), cos(pi / 6)))
  across <- rbind(c(0, 0), 100 * c(cos(pi / 6), -sin(pi / 6)))
  expect_equal(cx(along)[1, 2],  0.8 * exp(-100 / 300))
  expect_equal(cx(across)[1, 2], 0.8 * exp(-100 / 150))
  # An isotropic model carries no such attribute, so nothing else changes.
  expect_null(attr(spatialkit:::.vgm_correlation_fn(m[, 1:3]), "cor_xy"))

  # summarize_by_cell() uses it: one cell of points strung out east-west.
  pts <- sf::st_as_sf(data.frame(x = seq(0, 450, by = 50), y = 0, v = 1:10,
                                 poly_id = 1L),
                      coords = c("x", "y"), crs = 32632)
  sac <- structure(300, class = c("sac_range", "numeric"),
                   variogram_model = m, crs = sf::st_crs(32632))
  out <- summarize_by_cell(pts, "v", deff = "variogram", sac = sac)
  d <- as.matrix(stats::dist(cbind(seq(0, 450, by = 50), 0)))
  R <- 0.8 * exp(-d / 60); diag(R) <- 1
  expect_equal(attr(out, "deff_applied")$deff, sum(R) / 10)
})


# --- summarize_by_cell(): inputs and the join onto cells_sf ------------------

test_that("a response_var that cannot be summarised is a warning, not a silent skip", {
  pts <- .r2_points()
  expect_warning(out <- summarize_by_cell(pts, "vall"),
                 "response_var 'vall' is not a column")
  expect_false(any(grepl("^resp_", names(out))))
  pts$flag <- pts$v > 0
  expect_warning(summarize_by_cell(pts, "flag"), "'flag' are not numeric")
  expect_warning(summarize_by_cell(pts, "v", predictor_vars = c("v", "nope")),
                 "predictor_vars 'nope' are not columns")
  # Two names used to stop with "the condition has length > 1".
  expect_error(summarize_by_cell(pts, c("v", "flag")),
               "`response_var` must be a single column name")
})

test_that("a double ID of 100000 still joins an integer cell ID", {
  cells <- sf::st_sf(poly_id = c(99999L, 100000L, 100001L),
                     geometry = sf::st_sfc(.r2_sq(0, 0), .r2_sq(10, 0),
                                           .r2_sq(20, 0), crs = 32632))
  set.seed(2)
  pts <- sf::st_as_sf(data.frame(x = runif(30, 0, 30), y = runif(30, 0, 10),
                                 v = rnorm(30)),
                      coords = c("x", "y"), crs = 32632)
  a <- assign_features_to_polygons(pts, cells)
  a$poly_id <- as.double(a$poly_id)      # as from a GeoPackage Integer64 field
  res <- .r2_warnings(summarize_by_cell(a, "v", cells_sf = cells))
  out <- res$value
  # as.character(1e5) is "1e+05": that cell came back NA and its 11 points
  # were lost, with no R warning.
  expect_identical(out$poly_id, c("99999", "100000", "100001"))
  expect_identical(sum(out$n), 30L)
  expect_length(res$warnings, 0L)
  expect_identical(spatialkit:::.id_as_character(c(1e5, 2.5, NA, -0)),
                   c("100000", "2.5", NA, "0"))
})

test_that("summarised IDs that match no cell are reported", {
  pts <- .r2_points()
  cells <- sf::st_sf(poly_id = 1:3,
                     geometry = sf::st_sfc(lapply(0:2, function(i) .r2_sq(10 * i, 0)),
                                           crs = 32632))
  expect_warning(out <- summarize_by_cell(pts, "v", cells_sf = cells),
                 "1 of the 4 summarised cell ID\\(s\\) \\(4\\) match no")
  expect_identical(nrow(out), 3L)
})

test_that("cells keyed by 'id' are joined like assign_features_to_polygons() reads them", {
  cells <- sf::st_sf(id = 1:3, cell_id = 3:1,
                     geometry = sf::st_sfc(.r2_sq(0, 0), .r2_sq(10, 0),
                                           .r2_sq(20, 0), crs = 32632))
  set.seed(3)
  pts <- sf::st_as_sf(data.frame(x = c(runif(4, 0, 10), runif(6, 10, 20),
                                       runif(8, 20, 30)),
                                 y = runif(18, 0, 10), v = rnorm(18)),
                      coords = c("x", "y"), crs = 32632)
  a <- assign_features_to_polygons(pts, cells)      # takes the 'id' column
  out <- summarize_by_cell(a, "v", cells_sf = cells, area = TRUE)
  expect_s3_class(out, "sf")
  # Joined on 'id', not on the differently numbered 'cell_id': each count
  # sits on its own polygon.
  expect_identical(out$n, c(4L, 6L, 8L))
  expect_equal(sf::st_coordinates(sf::st_centroid(sf::st_geometry(out)))[, 1],
               c(5, 15, 25))
  expect_true(all(c("cell_area", "n_per_area") %in% names(out)))

  # No ID column at all: a warning and a plain table, or an error when an
  # area was asked for.  It used to be one log line and a plain table.
  nocol <- cells[, "geometry"]
  expect_warning(plain <- summarize_by_cell(a, "v", cells_sf = nocol),
                 "has none of the ID columns")
  expect_false(inherits(plain, "sf"))
  expect_error(summarize_by_cell(a, "v", cells_sf = nocol, area = TRUE),
               "`area = TRUE` needs the cells' areas")
  expect_warning(summarize_by_cell(a, "v", cells_sf = sf::st_drop_geometry(cells)),
                 "`cells_sf` is not an sf object")
})

test_that("agg_funs takes a bare function or function names", {
  pts <- .r2_points()
  med <- tapply(pts$v, pts$poly_id, median)
  out <- summarize_by_cell(pts, "v", agg_funs = median)
  # It used to become the mean, with only a log line.
  expect_equal(out$resp_median_v, as.numeric(med))
  out2 <- summarize_by_cell(pts, "v", agg_funs = c("median", "sum"))
  expect_equal(out2$resp_median_v, as.numeric(med))
  expect_equal(out2$resp_sum_v, as.numeric(tapply(pts$v, pts$poly_id, sum)))
  expect_warning(out3 <- summarize_by_cell(pts, "v", agg_funs = 3),
                 "Falling back to the mean")
  expect_true("resp_mean_v" %in% names(out3))
})
