# make_folds(): the advice of the k-reduction warnings (BUG-43) and the
# order of kNNDM's warning and min_train refusal (BUG-159).
#
# k = n gives the same one-per-fold folds as a k cut to n, so the advice
# names k <= n - 1, and none where no smaller k exists.  kNNDM tests
# min_train on the cut k before it warns that it returns n folds, and its
# refusal offers more folds only where n folds would meet min_train.

fk_pts <- function(n, seed) {
  set.seed(seed)
  sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, 1e4), y = 5e6 + runif(n, 0, 1e4)),
               coords = c("x", "y"), crs = 32632)
}
fk_grid <- function(m = 8)
  sf::st_as_sf(expand.grid(x = 5e5 + seq(0, 1.5e4, length.out = m),
                           y = 5e6 + seq(0, 1.5e4, length.out = m)),
               coords = c("x", "y"), crs = 32632)
fk_collect <- function(expr) {
  w <- list()
  val <- withCallingHandlers(expr, warning = function(c) {
    w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning")
  })
  list(value = val,
       reduced = Filter(function(c) inherits(c, "spatialkit_k_reduced"), w))
}
fk_tests <- function(f) sort(vapply(f$folds, function(s) paste(sort(s$test), collapse = ","), ""))

test_that("random_kfold's advice names k <= n - 1, which gives fewer, larger folds", {
  p5 <- fk_pts(5, 3)
  res <- r2_quiet(fk_collect(make_folds(p5, k = 10, method = "random_kfold", seed = 1)))
  expect_length(res$reduced, 1L)
  msg <- conditionMessage(res$reduced[[1L]])
  expect_match(msg, "This is leave-one-out; pass k <= 4 for fewer, larger folds\\.$")
  expect_false(grepl("k <= 5", msg, fixed = TRUE))
  # The counts the condition carries are unchanged.
  expect_identical(res$reduced[[1L]]$k_requested, 10L)
  expect_identical(res$reduced[[1L]]$k, 5L)
  expect_identical(res$value$k, 5L)
  # k = 5, the bound the old text gave, returns the same folds as the cut;
  # k = 4, the bound now named, returns fewer, larger folds and no warning.
  five <- r2_quiet(fk_collect(make_folds(p5, k = 5, method = "random_kfold", seed = 1)))
  expect_identical(fk_tests(five$value), fk_tests(res$value))
  four <- r2_quiet(fk_collect(make_folds(p5, k = 4, method = "random_kfold", seed = 1)))
  expect_length(four$reduced, 0L)
  expect_identical(four$value$k, 4L)
  expect_identical(max(lengths(lapply(four$value$folds, `[[`, "test"))), 2L)
})

test_that("on 2 rows or 2 groups the k-reduction warning gives no advice", {
  p2 <- fk_pts(2, 4)
  res <- r2_quiet(fk_collect(make_folds(p2, k = 5, method = "random_kfold", seed = 1)))
  expect_length(res$reduced, 1L)
  expect_identical(conditionMessage(res$reduced[[1L]]), paste0(
    "make_folds(random_kfold): returning 2 folds, not the k = 5 asked for: ",
    "the layer has 2 rows with usable coordinates and every fold holds out ",
    "at least one. This is leave-one-out."))
  # No smaller k exists: k = 1 is raised to 2 and gives the same folds.
  one <- r2_quiet(fk_collect(make_folds(p2, k = 1, method = "random_kfold", seed = 1)))
  expect_identical(one$value$k, 2L)
  expect_identical(fk_tests(one$value), fk_tests(res$value))

  p6 <- fk_pts(6, 6); p6$g <- rep(c("a", "b"), 3)
  res <- r2_quiet(fk_collect(make_folds(p6, k = 5, method = "leave_location_out",
                                        group_var = "g", seed = 1)))
  expect_length(res$reduced, 1L)
  expect_identical(conditionMessage(res$reduced[[1L]]), paste0(
    "make_folds(leave_location_out): returning 2 folds, not the k = 5 asked ",
    "for: `g` has 2 distinct groups and a group is never split across folds. ",
    "This is leave-one-group-out."))
})

test_that("leave_location_out's advice names k <= n_g - 1, which holds out several groups", {
  p12 <- fk_pts(12, 5); p12$g <- rep(c("a", "b", "c", "d"), 3)
  res <- r2_quiet(fk_collect(make_folds(p12, k = 6, method = "leave_location_out",
                                        group_var = "g", seed = 1)))
  expect_length(res$reduced, 1L)
  msg <- conditionMessage(res$reduced[[1L]])
  expect_match(msg, "This is leave-one-group-out; pass k <= 3 to hold out several groups per fold\\.$")
  expect_false(grepl("k <= 4", msg, fixed = TRUE))
  three <- r2_quiet(fk_collect(make_folds(p12, k = 3, method = "leave_location_out",
                                          group_var = "g", seed = 1)))
  expect_length(three$reduced, 0L)
  expect_identical(three$value$k, 3L)
  groups_per_fold <- vapply(three$value$folds, function(s)
    length(unique(p12$g[s$test])), 1L)
  expect_identical(max(groups_per_fold), 2L)
})

test_that("knndm's advice names k <= n - 1, from the k that min_train allows", {
  p6 <- fk_pts(6, 7)
  res <- r2_quiet(fk_collect(make_folds(p6, k = 10, method = "knndm",
                                        prediction_points = fk_grid(), seed = 1)))
  expect_length(res$reduced, 1L)
  expect_match(conditionMessage(res$reduced[[1L]]),
               "This is leave-one-out; pass k <= 5 for fewer, larger folds\\.$")
  five <- r2_quiet(fk_collect(make_folds(p6, k = 5, method = "knndm",
                                         prediction_points = fk_grid(), seed = 1)))
  expect_length(five$reduced, 0L)
  expect_identical(five$value$k, 5L)

  # min_train = 0.85 needs k >= 7 (1 - 1/7 = 0.857); on 8 rows k = 7 is the
  # one k that gives fewer folds, and k = 6 is refused.
  p8 <- fk_pts(8, 1)
  res <- r2_quiet(fk_collect(make_folds(p8, k = 10, method = "knndm",
                                        prediction_points = fk_grid(), seed = 1,
                                        min_train = 0.85)))
  expect_identical(res$value$k, 8L)
  expect_length(res$reduced, 1L)
  expect_match(conditionMessage(res$reduced[[1L]]),
               "This is leave-one-out; pass k = 7 for fewer, larger folds \\(min_train = 0.85 needs k >= 7\\)\\.$")
  seven <- r2_quiet(fk_collect(make_folds(p8, k = 7, method = "knndm",
                                          prediction_points = fk_grid(), seed = 1,
                                          min_train = 0.85)))
  expect_identical(seven$value$k, 7L)
  expect_error(r2_quiet(make_folds(p8, k = 6, method = "knndm",
                                   prediction_points = fk_grid(), seed = 1,
                                   min_train = 0.85)),
               "min_train = 0.85 cannot be met with k = 6 folds")
  # On 12 rows the same min_train leaves k from 7 to 11.
  p12 <- fk_pts(12, 2)
  res <- r2_quiet(fk_collect(make_folds(p12, k = 20, method = "knndm",
                                        prediction_points = fk_grid(), seed = 1,
                                        min_train = 0.85)))
  expect_match(conditionMessage(res$reduced[[1L]]),
               "pass k from 7 to 11 for fewer, larger folds \\(min_train = 0.85 needs k >= 7\\)\\.$")
  # At min_train = 0.875 on 8 rows only k = 8 is allowed: no advice.
  res <- r2_quiet(fk_collect(make_folds(p8, k = 10, method = "knndm",
                                        prediction_points = fk_grid(), seed = 1,
                                        min_train = 0.875)))
  expect_identical(res$value$k, 8L)
  expect_match(conditionMessage(res$reduced[[1L]]), "This is leave-one-out\\.$")
})

test_that("knndm refuses min_train on the cut k without first warning that it returns n folds", {
  p8 <- fk_pts(8, 1)
  for (k in c(10, 20)) {
    w <- list()
    err <- tryCatch(
      withCallingHandlers(
        r2_quiet(make_folds(p8, k = k, method = "knndm", prediction_points = fk_grid(),
                            seed = 1, min_train = 0.9)),
        warning = function(c) { w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning") }),
      error = function(e) e)
    expect_s3_class(err, "error")
    expect_length(Filter(function(c) inherits(c, "spatialkit_k_reduced"), w), 0L)
    msg <- conditionMessage(err)
    expect_match(msg, sprintf(paste0(
      "min_train = 0.9 cannot be met with k = 8 folds \\(the k = %d asked for, ",
      "cut to the 8 rows with usable coordinates\\)"), k))
    # More folds than rows cannot be had: the advice is the lower min_train.
    expect_match(msg, "Use min_train <= 0.875\\.$")
    expect_false(grepl("more folds", msg, fixed = TRUE))
  }
  # k = n, nothing cut: the same advice, and no "(the k ... asked for" clause.
  err <- tryCatch(r2_quiet(make_folds(p8, k = 8, method = "knndm",
                                      prediction_points = fk_grid(), min_train = 0.9)),
                  error = function(e) e)
  expect_match(conditionMessage(err), "cannot be met with k = 8 folds: .*Use min_train <= 0.875\\.$")
  # Below n, more folds help where n folds meet min_train ...
  expect_error(r2_quiet(make_folds(p8, k = 4, method = "knndm", prediction_points = fk_grid(),
                                   min_train = 0.85)),
               "Use min_train <= 0.75, or more folds\\.$")
  ok <- r2_quiet(suppressWarnings(make_folds(p8, k = 8, method = "knndm",
                                             prediction_points = fk_grid(), seed = 1,
                                             min_train = 0.85)))
  expect_identical(ok$k, 8L)
  # ... and are not offered where even n folds fall short.
  err <- tryCatch(r2_quiet(make_folds(p8, k = 4, method = "knndm",
                                      prediction_points = fk_grid(), min_train = 0.9)),
                  error = function(e) e)
  expect_match(conditionMessage(err), "Use min_train <= 0.75\\.$")
  # At k_min == n more folds still help: min_train = 0.875 needs k = 8, and
  # 8 rows allow it, so the refusal at k = 4 offers more folds and k = 8 runs.
  expect_error(r2_quiet(make_folds(p8, k = 4, method = "knndm", prediction_points = fk_grid(),
                                   min_train = 0.875)),
               "Use min_train <= 0.75, or more folds\\.$")
  at_n <- r2_quiet(fk_collect(make_folds(p8, k = 8, method = "knndm",
                                         prediction_points = fk_grid(), seed = 1,
                                         min_train = 0.875)))
  expect_identical(at_n$value$k, 8L)
  expect_length(at_n$value$folds, 8L)
  expect_length(at_n$reduced, 0L)
})

test_that("the smallest k that kNNDM's min_train allows matches the refusal's test", {
  k_min <- spatialkit:::.knndm_k_min
  meets <- spatialkit:::.knndm_k_meets
  for (m in c(0.01, 0.3, 0.5, 0.6, 2 / 3, 0.75, 0.8, 0.85, 0.875, 0.9, 0.95,
              1 - 1 / 7, 1 - 1 / 13, 0.999)) {
    km <- k_min(m)
    expect_true(meets(km, m), info = format(m, digits = 17))
    if (km > 2) expect_false(meets(km - 1, m), info = format(m, digits = 17))
  }
  # Values at the test's 1e-9 tolerance edge, where the closed form can land
  # one below the smallest k that meets it: 0.8 + 1e-9 gives 5 there, which
  # the refusal stops (1 - 0.800000001 < 1/5 - 1e-9), so k_min moves up to 6.
  # The assertions are properties of the two functions, not last bits.
  for (m in c(0.8 + 1e-9, 1 - 1 / c(3, 4, 5, 7, 8, 10, 20, 100) + 1e-9,
              1 - 1 / c(3, 4, 5, 7, 8, 10, 20, 100) - 1e-9)) {
    km <- k_min(m)
    expect_true(meets(km, m), info = format(m, digits = 17))
    if (km > 2) expect_false(meets(km - 1, m), info = format(m, digits = 17))
  }
  expect_identical(k_min(0.5), 2)
  expect_identical(k_min(0.875), 8)
  expect_identical(k_min(1 - 1 / 7), 7)
})

test_that("?make_folds says the warning gives a setting only where there is one", {
  # The help page from man/ where the sources are at hand (the Rd_db
  # fallback otherwise), as test-core-docs.R's .docs_rd() reads it.
  root <- testthat::test_path("..", "..")
  flat <- function(x)
    paste(unlist(lapply(x, function(el) if (is.character(el)) el else flat(el))), collapse = "")
  file <- file.path(root, "man", "make_folds.Rd")
  db <- if (file.exists(file)) tools::parse_Rd(file)
        else tools::Rd_db("spatialkit")[["make_folds.Rd"]]
  expect_true(grepl("reason and, where there is one, a setting that works",
                    gsub("[[:space:]]+", " ", flat(db)), fixed = TRUE))
})
