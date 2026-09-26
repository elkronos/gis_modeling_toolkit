# tests/testthat/test-review2-aoa.R
# ---------------------------------------------------------------------------
# area_of_applicability(): a fractional chunk_size, duplicated training rows,
# a predictor dropped for zero variance, and all-zero importance weights.
# ---------------------------------------------------------------------------

r2_aoa_pts <- function(n, a = stats::rnorm(n), b = stats::rnorm(n), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sf::st_as_sf(
    data.frame(x = seq_len(n) * 10, y = seq_len(n) * 5, a = a, b = b),
    coords = c("x", "y"), crs = 32632
  )
}


test_that("a fractional chunk_size cannot leave dense-path DI values at zero", {
  set.seed(11)
  tr <- r2_aoa_pts(80)
  # 25 prediction points far outside the training predictor range.
  nd <- r2_aoa_pts(25, a = stats::runif(25, 8, 12), b = stats::runif(25, 8, 12))
  ref <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"),
                               use_fnn = FALSE)
  expect_identical(ref$n_inside, 0L)

  # 2.5 used to give fractional block starts, leaving every fifth row at its
  # initial DI of 0 -- inside the AOA -- and 16 training DI of 0 that moved
  # the threshold.
  for (cs in c(2.5, 1.5, 12.5)) {
    res <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"),
                                 use_fnn = FALSE, chunk_size = cs)
    expect_identical(res$n_inside, 0L)
    expect_true(all(res$aoa$DI > 0))
    expect_equal(res$aoa$DI, ref$aoa$DI, tolerance = 1e-10)
    expect_equal(res$train_DI, ref$train_DI, tolerance = 1e-10)
    expect_equal(res$threshold, ref$threshold, tolerance = 1e-10)
  }

  # Values that cannot be a row count are refused by name.
  for (bad in list(0, NA_real_, c(10, 20), "10", -3, Inf))
    expect_error(area_of_applicability(nd, train_sf = tr,
                                       predictor_vars = c("a", "b"),
                                       use_fnn = FALSE, chunk_size = bad),
                 "area_of_applicability\\(\\): `chunk_size` must be")
})


test_that("a zero threshold from duplicated training rows is explained", {
  set.seed(5)
  # 30 sites visited four times each, with covariates that do not change
  # between visits: every training row has an exact twin.
  site <- data.frame(a = stats::rnorm(30), b = stats::rnorm(30))
  rep_rows <- site[rep(seq_len(30), each = 4), ]
  tr <- r2_aoa_pts(120, a = rep_rows$a, b = rep_rows$b)
  tr$site <- rep(seq_len(30), each = 4)
  nd <- r2_aoa_pts(200, a = stats::rnorm(200), b = stats::rnorm(200))

  lines <- capture_spatialkit_log(
    res <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b")))
  expect_identical(res$threshold, 0)
  expect_identical(res$n_inside, 0L)
  expect_true(log_has(lines, "DI threshold is 0 because 120 of 120 training rows have an exact duplicate"))
  expect_true(log_has(lines, "leave_location_out"))
  out <- utils::capture.output(print(res))
  expect_true(any(grepl("120 of 120 training DI are 0: exact duplicates", out,
                        fixed = TRUE)))

  # Folds that keep a site's visits together give the threshold its meaning
  # back, and say nothing.
  fo <- make_folds(tr, k = 5, method = "leave_location_out", group_var = "site",
                   seed = 1)
  lines2 <- capture_spatialkit_log(
    res2 <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"),
                                  folds = fo))
  expect_gt(res2$threshold, 0)
  expect_gt(res2$n_inside, 150L)
  expect_false(log_has(lines2, "DI threshold is 0"))
  expect_false(any(grepl("exact duplicates",
                         utils::capture.output(print(res2)), fixed = TRUE)))

  # A threshold the caller supplied is theirs; nothing to explain.
  lines3 <- capture_spatialkit_log(
    area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"),
                          threshold = 0.5))
  expect_false(log_has(lines3, "DI threshold is 0"))
})


test_that("a prediction row that differs on a zero-variance predictor is outside", {
  set.seed(7)
  # A land-cover dummy that is 0 everywhere in the training region.
  tr <- r2_aoa_pts(150)
  tr$urban <- 0
  nd <- r2_aoa_pts(80)
  nd$urban <- rep(c(0, 1), each = 40)
  nd$urban[80] <- NA          # cannot be compared: left to the other predictors

  expect_warning(
    res <- area_of_applicability(nd, train_sf = tr,
                                 predictor_vars = c("a", "b", "urban")),
    "39 of 80 prediction row\\(s\\) take a value the training data never has on .*urban")
  expect_identical(res$dropped_vars, "urban")
  urb <- which(nd$urban == 1)
  expect_true(all(is.infinite(res$aoa$DI[urb])))
  expect_false(any(res$aoa$AOA[urb]))
  expect_identical(res$n_inside + res$n_outside + res$n_na, 80L)
  expect_gte(res$n_outside, 39L)

  # Rows that agree with the training constant, or lack the value, are judged
  # on the other predictors exactly as before.
  ref <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"))
  keep <- c(1:40, 80)
  expect_equal(res$aoa$DI[keep], ref$aoa$DI[keep], tolerance = 1e-12)
  expect_true(is.finite(res$aoa$DI[80]))
  expect_true(any(grepl("39 outside on a dropped predictor (DI = Inf)",
                        utils::capture.output(print(res)), fixed = TRUE)))

  # Nothing differs, nothing is said.
  nd0 <- nd[1:40, ]
  expect_no_warning(
    res0 <- area_of_applicability(nd0, train_sf = tr,
                                  predictor_vars = c("a", "b", "urban")))
  expect_false(any(is.infinite(res0$aoa$DI)))
})


test_that("all-zero importance weights give an AOA instead of an error", {
  set.seed(9)
  tr <- r2_aoa_pts(60)
  nd <- r2_aoa_pts(30, a = c(stats::rnorm(20), stats::rnorm(10, 6)))

  # One predictor whose permutation importance was not positive: the index is
  # invariant to the weight's scale, so zero is accepted and changes nothing.
  ref1 <- area_of_applicability(nd, train_sf = tr, predictor_vars = "a")
  expect_no_warning(
    one <- area_of_applicability(nd, train_sf = tr, predictor_vars = "a",
                                 weights = pmax(c(a = -0.002), 0)))
  expect_equal(one$aoa$DI, ref1$aoa$DI, tolerance = 1e-12)
  expect_identical(one$aoa$AOA, ref1$aoa$AOA)

  # Several predictors, all zero: weighted equally, and said so.
  ref2 <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"))
  expect_warning(
    two <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"),
                                 weights = c(a = 0, b = 0)),
    "every weight is zero \\(a, b\\)")
  expect_equal(two$aoa$DI, ref2$aoa$DI, tolerance = 1e-12)
  expect_equal(unname(two$weights), c(1, 1))

  # The only non-zero weight sat on a predictor dropped for zero variance:
  # what is left is one predictor weighted zero, which cannot matter either.
  tr$const <- 3; nd$const <- 3
  expect_no_warning(
    dropped <- area_of_applicability(nd, train_sf = tr,
                                     predictor_vars = c("a", "const"),
                                     weights = c(a = 0, const = 1)))
  expect_equal(dropped$aoa$DI, ref1$aoa$DI, tolerance = 1e-12)

  # Negative weights are still refused, with the advice that now works.
  expect_error(area_of_applicability(nd, train_sf = tr, predictor_vars = "a",
                                     weights = c(a = -1)),
               "pmax\\(importance, 0\\)")
})


test_that("fold labels are numbered by the rule cv_*() uses", {
  splits <- spatialkit:::.aoa_fold_splits
  tests_of <- function(sp) lapply(sp, `[[`, "test")

  # Numbers in numeric order: fold 3 is label 10.
  expect_identical(tests_of(splits(rep(c(10, 2, 1), each = 3), 9)),
                   list(7:9, 4:6, 1:3))
  # A factor by its own levels, even when they are not sorted.
  expect_identical(tests_of(splits(factor(c("s", "s", "n", "n"),
                                          levels = c("s", "n")), 4)),
                   list(1:2, 3:4))
  # Anything else in C (radix) order, whatever the session's collation.
  lab <- c("north", "North", "south", "South")
  expect_identical(tests_of(splits(lab, 4)), list(2L, 4L, 1L, 3L))
  # ... which is how cv_*() numbers the same labels.
  cv <- spatialkit:::.folds_from_labels(lab, data.frame(..row_id = 1:4),
                                        "cv_spatial")
  expect_identical(tests_of(cv), tests_of(splits(lab, 4)))
})


test_that("character fold labels do not follow the session's collation", {
  # as.factor() sorted under LC_COLLATE, which puts "north" before "North"
  # in en_US and after it in C.  Only testable where en_US is installed.
  ok <- suppressWarnings(tryCatch({
    withr::local_collate("en_US.UTF-8")
    grepl("en_US", Sys.getlocale("LC_COLLATE"))
  }, error = function(e) FALSE))
  skip_if_not(ok, "the en_US.UTF-8 collation is not available")
  lab <- c("north", "North", "south", "South")
  expect_identical(lapply(spatialkit:::.aoa_fold_splits(lab, 4), `[[`, "test"),
                   list(2L, 4L, 1L, 3L))
})
