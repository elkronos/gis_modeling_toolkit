# ===========================================================================
# Regressions from the third review of assignment and aggregation
# (assign_features_to_polygons(), summarize_by_cell(), the row records and
# the knitr echo of the log).  Every test here failed on the code before the
# fix it names.
# ===========================================================================

# Every R warning (message and class) and R message `expr` raises, with its
# value.  The console echo written to stderr is swallowed.
.r3_conditions <- function(expr) {
  w_msg <- character(0); w_cls <- list(); n_msg <- 0L
  utils::capture.output(
    val <- withCallingHandlers(expr,
      warning = function(w) {
        w_msg <<- c(w_msg, conditionMessage(w))
        w_cls <<- c(w_cls, list(class(w)))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        n_msg <<- n_msg + 1L
        invokeRestart("muffleMessage")
      }),
    type = "message")
  list(value = val, warnings = w_msg, classes = w_cls, messages = n_msg)
}

.r3_is_fallback <- function(res)
  vapply(res$classes, function(k) "spatialkit_deff_fallback" %in% k, logical(1))

.r3_sq <- function(x0, y0, s) {
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s),
                            c(x0, y0 + s), c(x0, y0))))
}

# A 3 x 3 grid of 100 m cells, numbered from the lower left a row at a time.
.r3_grid <- function() {
  sf::st_sf(poly_id = 1:9,
            geometry = sf::st_sfc(lapply(0:8, function(k)
              .r3_sq(5e5 + (k %% 3) * 100, 5e6 + (k %/% 3) * 100, 100)),
              crs = 32632))
}

# Six cells of 25 points on a correlated Gaussian field.
.r3_field <- function() {
  set.seed(77)
  n <- 150
  x <- rep(seq(0, 500, length.out = 6), each = 25) + runif(n, 0, 60)
  y <- runif(n, 0, 60)
  z <- as.numeric(t(chol(exp(-as.matrix(stats::dist(cbind(x, y))) / 40) +
                           diag(1e-6, n))) %*% rnorm(n))
  sf::st_as_sf(data.frame(x = x, y = y, z = z, w = z + rnorm(n),
                          poly_id = rep(1:6, each = 25)),
               coords = c("x", "y"), crs = 32632)
}

.r3_sac <- function(model = data.frame(model = c("Nug", "Exp"),
                                       psill = c(0.2, 0.8), range = c(0, 40)),
                    ...) {
  structure(120, class = c("sac_range", "numeric"), variogram_model = model,
            crs = sf::st_crs(32632), ...)
}

.r3_rejected <- function(reason = "fitted range exceeds the largest lag fitted") {
  structure(NA_real_, class = c("sac_range", "numeric"),
            rejected_reason = reason,
            variogram_model = data.frame(model = "Exp", psill = 1, range = 1e6),
            crs = sf::st_crs(32632))
}


# --- summarize_by_cell(): a rejected `sac` -----------------------------------

test_that("a rejected sac replaced by an estimate is not reported as a fallback", {
  skip_if_not_installed("gstat")
  pts <- .r3_field()
  # The internal estimate succeeds.  The rejected sac used to raise the
  # classed fallback warning ("Falling back to deff = 1") although the
  # estimate was then applied.
  local_mocked_bindings(estimate_sac_range = function(...) .r3_sac())
  res <- .r3_conditions(summarize_by_cell(pts, "z", deff = "variogram",
                                          sac = .r3_rejected()))
  expect_false(any(.r3_is_fallback(res)))
  expect_length(res$warnings, 1L)
  expect_match(res$warnings, "no usable range \\(fitted range exceeds")
  expect_match(res$warnings, "estimated from `response_var` instead")
  expect_true(all(res$value$deff_applied))
  expect_identical(attr(res$value, "deff_applied")$method, "variogram")
  # tryCatch() on the class no longer throws the corrected result away.
  kept <- suppressWarnings(tryCatch(
    summarize_by_cell(pts, "z", deff = "variogram", sac = .r3_rejected()),
    spatialkit_deff_fallback = function(w) "caught"))
  expect_s3_class(kept, "data.frame")
})

test_that("a rejected sac that cannot be replaced gives one fallback warning naming both reasons", {
  skip_if_not_installed("gstat")
  pts <- .r3_field()
  local_mocked_bindings(estimate_sac_range = function(...)
    .r3_rejected("variogram model did not converge"))
  res <- .r3_conditions(summarize_by_cell(pts, "z", deff = "variogram",
                                          sac = .r3_rejected()))
  expect_length(res$warnings, 1L)
  expect_true(.r3_is_fallback(res))
  expect_match(res$warnings, "supplied `sac` reports no usable range \\(fitted range exceeds")
  expect_match(res$warnings, "estimated from `response_var` reports no usable range \\(variogram model did not converge")
  expect_false(any(res$value$deff_applied))

  # No response to estimate from: the rejection is still named.
  res2 <- .r3_conditions(summarize_by_cell(pts, predictor_vars = "z",
                                           deff = "variogram",
                                           sac = .r3_rejected()))
  expect_length(res2$warnings, 1L)
  expect_true(.r3_is_fallback(res2))
  expect_match(res2$warnings, "fitted range exceeds.*no `response_var`")
})


# --- summarize_by_cell(): a detrended `sac` ----------------------------------

test_that("a residual (detrended) sac correcting response SEs is warned about", {
  pts <- .r3_field()
  res <- .r3_sac(detrended = TRUE, detrend_method = "ols")
  # Used as given, as documented, but no longer in silence.
  expect_warning(out <- summarize_by_cell(pts, "z", "w", deff = "variogram", sac = res),
                 "residuals on predictors \\(detrend = \"ols\"\\).*understated")
  expect_true(all(out$deff_applied))
  quiet <- summarize_by_cell(pts, "z", "w", deff = "variogram", sac = .r3_sac())
  expect_equal(out[["..se_resp_z"]], quiet[["..se_resp_z"]])
  # Nothing to warn about without response columns, or with a response variogram.
  expect_no_warning(summarize_by_cell(pts, predictor_vars = "w", deff = "variogram",
                                      sac = res))
  expect_no_warning(summarize_by_cell(pts, "z", deff = "variogram",
                                      sac = .r3_sac(detrended = FALSE)))
})


# --- summarize_by_cell(): the design-effect record ---------------------------

test_that("deff = 'kish' records a correction made to the predictor SEs only", {
  set.seed(5)
  n <- 400
  pts <- sf::st_as_sf(data.frame(x = runif(n, 0, 200), y = runif(n, 0, 200)),
                      coords = c("x", "y"), crs = 32632)
  xy <- sf::st_coordinates(pts)
  pts$poly_id <- 1L + (xy[, 1] %/% 50) + 4L * (xy[, 2] %/% 50)
  pts$v <- rnorm(n)                                    # response: unclustered
  pts$p <- rnorm(16, sd = 2)[pts$poly_id] + rnorm(n)   # predictor: clustered
  naive <- summarize_by_cell(pts, "v", "p")
  k <- summarize_by_cell(pts, "v", "p", deff = "kish")
  expect_identical(attr(k, "icc")$resp, 0)
  expect_gt(attr(k, "icc")$pred, 0.5)
  # The predictor SEs were inflated about elevenfold; the rows used to say
  # deff_applied = FALSE, with no attribute.
  expect_gt(stats::median(k$..se_pred_p / naive$..se_pred_p), 5)
  expect_true(all(k$deff_applied))
  rec <- attr(k, "deff_applied")
  expect_identical(rec$method, "kish")
  expect_gt(rec$icc_pred, 0.5)
  # `deff` and cell_weight belong to the primary variable, the response.
  expect_equal(rec$deff, rep(1, nrow(k)))
  expect_equal(k$cell_weight, k$n)
})

test_that("a pure-nugget variogram is applied as deff = 1, not reported as a fallback", {
  pts <- .r3_field()
  nug <- .r3_sac(model = data.frame(model = "Nug", psill = 1, range = 0))
  expect_no_warning(out <- summarize_by_cell(pts, "z", deff = "variogram", sac = nug))
  naive <- summarize_by_cell(pts, "z")
  expect_true(all(out$deff_applied))
  expect_equal(out[["..se_resp_z"]], naive[["..se_resp_z"]])
  expect_equal(attr(out, "deff_applied")$rbar, rep(0, nrow(out)))
  # So is a structured component with no sill.
  nug2 <- .r3_sac(model = data.frame(model = c("Nug", "Exp"), psill = c(1, 0),
                                     range = c(0, 40)))
  expect_no_warning(out2 <- summarize_by_cell(pts, "z", deff = "variogram", sac = nug2))
  expect_true(all(out2$deff_applied))
})

test_that("deff_max_n must be a number of at least 2 for deff = 'variogram'", {
  pts <- .r3_field()
  for (bad in list(1L, 0L, NA, c(10, 20))) {
    # 1 and 0 returned uncorrected SEs marked deff_applied = TRUE; NA stopped
    # with "missing value where TRUE/FALSE needed".
    expect_error(summarize_by_cell(pts, "z", deff = "variogram", sac = .r3_sac(),
                                   deff_max_n = bad),
                 "`deff_max_n` must be")
  }
  # Unused, so unchecked, by every other deff.
  expect_no_error(summarize_by_cell(pts, "z", deff = "kish", deff_max_n = 0))
})

test_that("the rows with no cell ID get a cell_weight like any other group", {
  set.seed(4)
  pts <- sf::st_as_sf(data.frame(x = 5e5 + runif(40, -100, 300),
                                 y = 5e6 + runif(40, 0, 300), v = rnorm(40)),
                      coords = c("x", "y"), crs = 32632)
  a <- assign_features_to_polygons(pts, .r3_grid(), keep_unassigned = TRUE)
  s <- summarize_by_cell(a, "v")
  na_row <- which(is.na(s$poly_id))
  expect_length(na_row, 1L)
  expect_gt(s$n[na_row], 0)
  # It was 0 beside n = 5 and a finite SE.
  expect_equal(s$cell_weight, s$n)
})

test_that("a fixed deff stays a scalar on the record when one cell is summarised", {
  pts <- sf::st_as_sf(data.frame(x = 5e5 + c(10, 20, 30), y = 5e6 + c(10, 20, 30),
                                 v = c(1, 2, 4)),
                      coords = c("x", "y"), crs = 32632)
  grid <- .r3_grid()
  out <- summarize_by_cell(assign_features_to_polygons(pts, grid), "v",
                           cells_sf = grid, deff = 2)
  # It came back as c(2, NA, NA, ...).
  expect_identical(attr(out, "deff_applied"), list(method = "fixed", deff = 2))
})

test_that("agg_funs = pkg::fn is named after the function", {
  pts <- .r3_field()
  out <- summarize_by_cell(pts, "z", agg_funs = stats::median)
  expect_true("resp_median_z" %in% names(out))
  expect_false("resp_agg1_z" %in% names(out))
})


# --- assign_features_to_polygons() ------------------------------------------

test_that("the features' own column named like the polygons' fallback ID column is kept", {
  cells <- .r3_grid()
  names(cells)[names(cells) == "poly_id"] <- "id"
  pts <- sf::st_as_sf(data.frame(id = c("site-A", "site-B", "site-C"),
                                 x = 5e5 + c(10, 150, 250), y = 5e6 + c(10, 150, 250),
                                 v = 1:3),
                      coords = c("x", "y"), crs = 32632)
  # The site IDs were dropped, with a warning about a collision that the
  # result never had.
  expect_no_warning(a <- assign_features_to_polygons(pts, cells))
  expect_identical(a$id, c("site-A", "site-B", "site-C"))
  expect_identical(a$poly_id, c(1L, 5L, 9L))
  # A column named `polygon_id_col` is still replaced, with the warning.
  expect_warning(b <- assign_features_to_polygons(a, cells), "'poly_id'")
  expect_identical(b$id, a$id)
  expect_identical(b$poly_id, a$poly_id)
})

test_that("an exact tie in overlap area does not depend on the polygon row order", {
  grid <- .r3_grid()
  # A 40 m square split evenly across the edge of cells 1 and 2, and a 20 m
  # square split evenly among cells 1, 2, 4 and 5.
  f <- sf::st_sf(fid = 1:3, geometry = sf::st_sfc(
    .r3_sq(5e5 + 80, 5e6 + 10, 40), .r3_sq(5e5 + 130, 5e6 + 130, 60),
    .r3_sq(5e5 + 90, 5e6 + 90, 20), crs = 32632))
  fwd <- assign_features_to_polygons(f, grid)
  rev <- assign_features_to_polygons(f, grid[9:1, ])
  # Reversed rows used to give cells 2 and 5, and ties$n was 0.
  expect_identical(fwd$poly_id, c(1L, 5L, 1L))
  expect_identical(rev$poly_id, fwd$poly_id)
  expect_identical(attr(fwd, "ties")$n, 2L)
  expect_identical(attr(fwd, "ties")$which, c(1L, 3L))
  expect_identical(attr(rev, "ties")$which, c(1L, 3L))
  # "first" keeps the polygon row order, as documented, and still counts.
  first <- assign_features_to_polygons(f, grid[9:1, ], tie_break = "first")
  expect_identical(first$poly_id, c(2L, 5L, 5L))
  expect_identical(attr(first, "ties")$n, 2L)
})

test_that("a largest-overlap assignment does not leak sf's attribute warning", {
  grid <- .r3_grid()
  f <- sf::st_sf(fid = 1:2, geometry = sf::st_sfc(
    .r3_sq(5e5 + 10, 5e6 + 10, 30), .r3_sq(5e5 + 130, 5e6 + 130, 60), crs = 32632))
  # "attribute variables are assumed to be spatially constant throughout all
  # geometries" came with every call.
  expect_no_warning(a <- assign_features_to_polygons(f, grid))
  expect_identical(a$poly_id, c(1L, 5L))
})


# --- the row records ---------------------------------------------------------

test_that("dplyr's row verbs drop a row record as `[` does", {
  grid <- .r3_grid()
  e <- sf::st_as_sf(data.frame(x = 5e5 + c(100, 20, 200, 60, 150),
                               y = 5e6 + c(50, 20, 150, 60, 200), v = 1:5),
                    coords = c("x", "y"), crs = 32632)
  a <- assign_features_to_polygons(e, grid)
  expect_identical(attr(a, "ties")$n_rows, 5L)
  # filter() on 5 rows returned 3 rows still reporting the parent's ties.
  for (r in list(dplyr::filter(a, v >= 3), dplyr::slice(a, 1:3),
                 dplyr::arrange(a, dplyr::desc(v)))) {
    expect_null(attr(r, "ties"))
    expect_s3_class(r, "sf")
    expect_false(inherits(r, "spatialkit_rows"))
  }
  # The same for prep_model_data()'s record.
  dat <- sf::st_as_sf(data.frame(x = 1:5, y = 5:1, resp = c(1, 2, NA, 4, 5),
                                 pred = c(1, 2, 3, 4, Inf)),
                      coords = c("x", "y"), crs = 32632)
  clean <- prep_model_data(dat, "resp", "pred")
  expect_null(attr(dplyr::filter(clean, resp > 1), "dropped"))
  # And for a geometry-free frame carrying a record (what st_drop_geometry()
  # of such a layer is).
  d <- spatialkit:::.set_row_record(data.frame(v = 1:5), "ties",
                                    list(n = 1L, which = 2L, rule = "first"))
  expect_null(attr(dplyr::filter(d, v >= 3), "ties"))
  expect_false(inherits(dplyr::filter(d, v >= 3), "spatialkit_rows"))
  # Binding such frames keeps the first one's record, as documented; its
  # stamp no longer matches, so the package's readers ignore it.
  expect_identical(attr(rbind(d, d), "ties")$n_rows, 5L)
  expect_null(spatialkit:::.get_row_record(rbind(d, d), "ties"))
})


# --- the knitr echo of the log ----------------------------------------------

test_that("under knitr a caution raised as a warning is not repeated as a message", {
  old <- spatialkit_quiet(FALSE)
  withr::defer(spatialkit_quiet(old))
  withr::local_options(knitr.in.progress = TRUE)
  set.seed(1)
  pts <- sf::st_as_sf(data.frame(x = runif(20, 0, 200), y = runif(20, 0, 200),
                                 v = rnorm(20), poly_id = rep(1:4, each = 5)),
                      coords = c("x", "y"), crs = 32632)
  # A design-effect fallback: it was a message and a warning.
  res <- .r3_conditions(summarize_by_cell(pts, "v", deff = 0.5))
  expect_identical(res$messages, 0L)
  expect_length(res$warnings, 1L)
  # The .log_warn() + warning() pairs of the cross-validation code: the
  # random-folds fallback.
  res <- .r3_conditions(spatialkit:::.remap_folds(NULL, 1:20))
  expect_identical(res$messages, 0L)
  expect_length(res$warnings, 1L)
  # A caution that is only logged still reaches the document.
  log_only <- function() { spatialkit:::.log_warn("only logged %d", 1L); invisible(NULL) }
  res <- .r3_conditions(log_only())
  expect_identical(res$messages, 1L)
  expect_length(res$warnings, 0L)
  # A warning() that does not follow the log line directly does not count.
  apart <- function() {
    spatialkit:::.log_warn("logged %d", 2L)
    x <- 1
    warning("raised")
  }
  res <- .r3_conditions(apart())
  expect_identical(res$messages, 1L)
  expect_length(res$warnings, 1L)
})
