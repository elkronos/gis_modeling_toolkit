# prep_model_data() says what it removes and marks what it lets through.
#
# It dropped incomplete rows with a log line only: 432 of 1,050 rows (41%)
# went with no R warning while the mean response moved from 0.03 to -0.62.  A
# predictor with one value in every row and repeated locations passed through
# unmarked.

.fbp_pts <- function(n = 1000, seed = 2) {
  set.seed(seed)
  d <- sf::st_as_sf(
    data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
               resp = stats::rnorm(n), p1 = stats::rnorm(n), const = 3,
               grp = sample(letters[1:3], n, TRUE)),
    coords = c("x", "y"), crs = 32632)
  d
}

.fbp_warnings <- function(expr) {
  w <- list()
  v <- withCallingHandlers(expr, warning = function(c) {
    w[[length(w) + 1L]] <<- c
    invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w,
       classes = vapply(w, function(c) class(c)[1L], character(1)))
}


test_that("dropped rows raise a classed warning with the counts", {
  d <- .fbp_pts()
  d$p1[d$resp > 0.25] <- NA                    # missing on the high responses
  n_na <- sum(is.na(d$p1))
  r <- .fbp_warnings(prep_model_data(d, "resp", c("p1", "grp")))
  expect_identical(nrow(r$value), 1000L - n_na)
  # There was no R warning at all.
  expect_identical(r$classes, "spatialkit_rows_dropped")
  msg <- conditionMessage(r$warnings[[1]])
  expect_match(msg, sprintf("dropped %d of 1000 row(s) (%.1f%%)", n_na, n_na / 10),
               fixed = TRUE)
  expect_match(msg, sprintf("%d with a missing value, 0 with a non-finite value, 0 with an empty or non-finite geometry",
                            n_na), fixed = TRUE)
  expect_match(msg, sprintf("Rows affected by column: p1 %d.", n_na), fixed = TRUE)
  expect_match(msg, sprintf("The fit uses the %d row(s) that remain", 1000L - n_na),
               fixed = TRUE)
  # The missingness follows the response, and the message says by how much.
  expect_match(msg, "not missing at random", fixed = TRUE)

  dr <- attr(r$value, "dropped")
  expect_identical(dr$n, n_na)
  expect_identical(dr$by_column, c(resp = 0L, p1 = n_na, grp = 0L))
  y <- d$resp
  expect_equal(dr$response_std_diff,
               (mean(y[is.na(d$p1)]) - mean(y[!is.na(d$p1)])) / stats::sd(y))
  expect_gt(dr$response_std_diff, 1)
  # It can be caught, and muffled, by class.
  expect_identical(
    tryCatch(prep_model_data(d, "resp", "p1"),
             spatialkit_rows_dropped = function(c) "caught"), "caught")
  expect_no_warning(withCallingHandlers(
    prep_model_data(d, "resp", "p1"),
    spatialkit_rows_dropped = function(c) invokeRestart("muffleWarning")))
})

test_that("the drop warning counts each reason and each column", {
  d <- sf::st_as_sf(data.frame(x = 1:10, y = 10:1,
                               resp = c(1, 2, NA, 4, 5, 6, 7, NaN, 9, 10),
                               pred = c(1, 2, 3, 4, Inf, 6, 7, 8, 9.5, 11)),
                    coords = c("x", "y"), crs = 32632)
  sf::st_geometry(d)[6] <- sf::st_sfc(sf::st_point(), crs = 32632)
  r <- .fbp_warnings(prep_model_data(d, "resp", "pred"))
  expect_identical(r$classes, "spatialkit_rows_dropped")
  msg <- conditionMessage(r$warnings[[1]])
  expect_match(msg, "dropped 4 of 10 row(s) (40.0%): 2 with a missing value, 1 with a non-finite value, 1 with an empty or non-finite geometry",
               fixed = TRUE)
  expect_match(msg, "Rows affected by column: resp 2, pred 1, (geometry) 1.", fixed = TRUE)
  dr <- attr(r$value, "dropped")
  expect_identical(dr$by_column, c(resp = 2L, pred = 1L, `(geometry)` = 1L))
  # The counts the record always had are unchanged.
  expect_identical(dr$n, 4L)
  expect_identical(dr$n_geometry, 1L)
  expect_identical(dr$which, c(3L, 5L, 6L, 8L))
  expect_identical(dr$reason, c("missing", "non_finite", "geometry", "missing"))
  # Two dropped rows with a response, too few to read anything into.
  expect_false(grepl("not missing at random", msg, fixed = TRUE))
})

test_that("nothing dropped, nothing raised", {
  d <- .fbp_pts(n = 60)
  r <- .fbp_warnings(prep_model_data(d, "resp", c("p1", "grp")))
  expect_length(r$warnings, 0L)
  dr <- attr(r$value, "dropped")
  expect_identical(dr$n, 0L)
  expect_true(is.na(dr$response_std_diff))
  expect_identical(dr$constant_predictors, character(0))
  expect_identical(dr$n_duplicate_locations, 0L)
  expect_true(all(dr$by_column == 0L))
})

test_that("a constant predictor is warned about by name and recorded", {
  d <- .fbp_pts(n = 80)
  d$one <- "a"
  r <- .fbp_warnings(prep_model_data(d, "resp", c("p1", "const", "one", "grp")))
  expect_identical(r$classes, "spatialkit_constant_predictor")
  msg <- conditionMessage(r$warnings[[1]])
  expect_match(msg, "const", fixed = TRUE)
  expect_match(msg, "one", fixed = TRUE)
  expect_match(msg, "take a single value in all 80 row(s) kept", fixed = TRUE)
  expect_false(grepl("p1", msg, fixed = TRUE))
  expect_identical(attr(r$value, "dropped")$constant_predictors, c("const", "one"))
  # The predictor is kept: the caller decides.
  expect_true(all(c("const", "one") %in% names(r$value)))

  # Constant only among the rows that survive cleaning.
  d2 <- .fbp_pts(n = 80)
  d2$v <- 1; d2$v[1:5] <- 2; d2$p1[1:5] <- NA
  r2 <- .fbp_warnings(prep_model_data(d2, "resp", c("p1", "v")))
  expect_setequal(r2$classes, c("spatialkit_rows_dropped",
                                "spatialkit_constant_predictor"))
  expect_identical(attr(r2$value, "dropped")$constant_predictors, "v")
})

test_that("repeated locations are counted in the record, and not warned about", {
  d <- .fbp_pts(n = 100)
  d <- rbind(d, d[1:27, ])
  logged <- capture_spatialkit_log(
    r <- .fbp_warnings(prep_model_data(d, "resp", "p1")))
  expect_length(r$warnings, 0L)
  expect_identical(attr(r$value, "dropped")$n_duplicate_locations, 27L)
  expect_true(log_has(logged, "27 of 127 row\\(s\\) kept share their coordinates"))
  # Counted after cleaning: a repeat that is dropped is not a repeat.
  d$p1[101:110] <- NA
  r2 <- suppressWarnings(prep_model_data(d, "resp", "p1"))
  expect_identical(attr(r2, "dropped")$n_duplicate_locations, 17L)
})

test_that("prediction mode logs and records, and raises neither warning", {
  d <- .fbp_pts(n = 60)
  d$p1[c(4, 9)] <- NA
  d$const2 <- 1
  logged <- capture_spatialkit_log(
    r <- .fbp_warnings(prep_model_data(d, "resp", c("p1", "const2"),
                                       require_response = FALSE)),
    level = logger::WARN)
  # predict() returns NA at exactly these rows, so nothing is lost silently,
  # and one row of newdata is "constant" in every column.
  expect_length(r$warnings, 0L)
  expect_identical(nrow(r$value), 58L)
  expect_identical(attr(r$value, "dropped")$n, 2L)
  expect_true(log_has(logged, "dropping 2 row\\(s\\)"))
})

test_that("the marks ride in the record and do not outlive the rows", {
  d <- .fbp_pts(n = 60)
  d <- rbind(d, d[1:3, ])
  out <- suppressWarnings(prep_model_data(d, "resp", c("p1", "const")))
  expect_identical(attr(out, "dropped")$n_duplicate_locations, 3L)
  expect_null(attr(out[1:10, ], "dropped"))
  expect_null(attr(out, "n_duplicate_locations"))
})

test_that("a fit still reports the rows prep removed as n_dropped", {
  skip_if_not_installed("brms")
  d <- .fbp_pts(n = 150)
  d$p1[1:7] <- NA
  local_mocked_bindings(brm = function(...) structure(list(), class = "fbp_stub"),
                        .package = "brms")
  r <- .fbp_warnings(suppressMessages(fit_bayesian_spatial_model(
    d, "resp", "p1", gp_k = 6, compute_loo = FALSE, check_convergence = FALSE)))
  expect_identical(r$value$info$n_dropped, 7L)
  expect_identical(r$value$n, 143L)
  # The fitter passes the warning on, once.
  expect_identical(sum(r$classes == "spatialkit_rows_dropped"), 1L)
})
