# predict.bayesian_fit() and the levels of a categorical predictor.
#
# brms refuses a whole `newdata` whenever the LEVEL SET of one of its factors
# is not a subset of the levels brms kept at fit time, and predict() returned
# that failure as NA for every row, with a log line and no R condition.  Three
# things triggered it:
#   1. one row carrying an unseen value;
#   2. a factor that only lists an extra level no row carries;
#   3. a model fitted on a subset whose factor kept a level with no rows, where
#      even the training rows were refused (fitted() and summary() stopped).
#
# No Stan here except the last test (opt-in).  The engine is an empty brmsfit
# (brms::brm(empty = TRUE): no compilation, no draws), and posterior_epred()
# is replaced by a stand-in that runs brms's own newdata validation, so the
# level rule being tested is brms's and not a copy of it.

.fbl_pts <- function(n = 90, seed = 3) {
  set.seed(seed)
  d <- sf::st_as_sf(
    data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
               cover = factor(sample(c("forest", "crop", "urban"), n, TRUE))),
    coords = c("x", "y"), crs = 32632)
  d$z <- c(forest = 0, crop = 1, urban = 3)[as.character(d$cover)] +
    stats::rnorm(n, 0, 0.5)
  d
}

# A bayesian_fit over an empty brmsfit, built the way the fitter builds one.
.fbl_fit <- function(d, store_levels = TRUE) {
  xy <- sf::st_coordinates(d)
  cs <- list(x_center = mean(xy[, 1]), x_scale = stats::sd(xy[, 1]),
             y_center = mean(xy[, 2]), y_scale = stats::sd(xy[, 2]))
  df <- sf::st_drop_geometry(d)
  df$`..x` <- (xy[, 1] - cs$x_center) / cs$x_scale
  df$`..y` <- (xy[, 2] - cs$y_center) / cs$y_scale
  eng <- suppressMessages(brms::brm(
    z ~ cover + gp(..x, ..y, k = 5, c = 1.5, scale = FALSE), data = df,
    empty = TRUE))
  info <- list(coord_scaling = cs)
  # What the fitter stores: the levels that occur in the training rows.
  if (store_levels)
    info$predictor_levels <- list(cover = list(
      levels = levels(droplevels(as.factor(df$cover))), ordered = FALSE))
  new_spatial_fit("bayesian_fit", engine = eng, formula = z ~ cover,
                  response_var = "z", predictor_vars = "cover", data_sf = d,
                  info = info)
}

# posterior_epred() stand-in: brms's own validation of `newdata` against the
# fit, then one number per row that depends only on the row's level (10, 20,
# 30 by level, in the fitted order), so a reference-level prediction for a row
# that should have been NA would show.
.fbl_mock_epred <- function(seen) {
  function(object, newdata, ...) {
    seen$na_levels <- c(seen$na_levels, anyNA(newdata$cover))
    seen$n <- c(seen$n, nrow(newdata))
    v <- brms::validate_newdata(newdata, object, check_response = FALSE)
    matrix(rep(10 * as.integer(v$cover), each = 4L), nrow = 4L)
  }
}

.fbl_warnings <- function(expr) {
  w <- list()
  v <- withCallingHandlers(expr, warning = function(c) {
    w[[length(w) + 1L]] <<- c
    invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w)
}


test_that("the level helpers flag a value, not a listed level, and follow the training rows", {
  d  <- .fbl_pts()
  df <- sf::st_drop_geometry(d)
  lv <- spatialkit:::.bayes_predictor_levels(df, c("cover", "z"))
  # Only categorical predictors, with the levels that occur.
  expect_named(lv, "cover")
  expect_identical(lv$cover$levels, c("crop", "forest", "urban"))
  expect_false(lv$cover$ordered)
  # A training level with no rows is not a fitted level: brms drops it.
  lv2 <- spatialkit:::.bayes_predictor_levels(df[df$cover != "urban", ], "cover")
  expect_identical(lv2$cover$levels, c("crop", "forest"))
  # A character column gets its sorted distinct values, as brms makes them.
  dfc <- df; dfc$cover <- as.character(dfc$cover)
  expect_identical(spatialkit:::.bayes_predictor_levels(dfc, "cover")$cover$levels,
                   c("crop", "forest", "urban"))

  stub <- list(info = list(predictor_levels = lv), predictor_vars = "cover")
  nd <- df[1:6, ]
  # An unseen VALUE is flagged, by row.
  ndA <- nd; ndA$cover <- as.character(ndA$cover); ndA$cover[6] <- "water"
  a <- spatialkit:::.bayes_unseen_levels(stub, ndA)
  expect_identical(a$unseen, c(rep(FALSE, 5), TRUE))
  expect_match(a$named, "'cover': 'water'", fixed = TRUE)
  # An extra level no row carries is not.
  ndC <- nd; ndC$cover <- factor(as.character(ndC$cover),
                                 levels = c("crop", "forest", "urban", "water"))
  expect_false(any(spatialkit:::.bayes_unseen_levels(stub, ndC)$unseen))
  # A missing value is not an unseen level; the row cleaning handles it.
  ndN <- nd; ndN$cover[2] <- NA
  expect_false(any(spatialkit:::.bayes_unseen_levels(stub, ndN)$unseen))
  # A level the training factor lists but has no rows for is unseen.
  stub2 <- list(info = list(predictor_levels = lv2), predictor_vars = "cover")
  u <- spatialkit:::.bayes_unseen_levels(stub2, nd)
  expect_identical(u$unseen, as.character(nd$cover) == "urban")
})

test_that("one row with an unseen level is NA; the others are predicted, with one classed warning", {
  skip_if_not_installed("brms")
  d   <- .fbl_pts()
  fit <- .fbl_fit(d)
  seen <- new.env()
  local_mocked_bindings(posterior_epred = .fbl_mock_epred(seen), .package = "brms")

  ref <- suppressMessages(predict(fit, newdata = d[1:6, ]))
  expect_equal(ref, 10 * as.integer(d$cover[1:6]))

  ndA <- d[1:6, ]; ndA$cover <- as.character(ndA$cover); ndA$cover[6] <- "water"
  r <- .fbl_warnings(suppressMessages(predict(fit, newdata = ndA)))
  # It was NA for all six, with no R warning.
  expect_equal(r$value[1:5], ref[1:5])
  expect_true(is.na(r$value[6]))
  expect_length(r$warnings, 1L)
  expect_s3_class(r$warnings[[1]], "spatialkit_new_levels")
  msg <- conditionMessage(r$warnings[[1]])
  expect_match(msg, "1 of 6 row(s)", fixed = TRUE)
  expect_match(msg, "'cover': 'water'", fixed = TRUE)
  expect_match(msg, "the other 5 row(s) are predicted", fixed = TRUE)
  # The unseen row never reached brms, as a level or as an NA (brms predicts
  # the reference level for an NA factor value).
  expect_false(any(seen$na_levels))

  # The same with a factor column, and with draws = TRUE: an NA column.
  ndF <- d[1:6, ]
  ndF$cover <- factor(c(as.character(d$cover[1:5]), "water"))
  dr <- suppressWarnings(suppressMessages(predict(fit, newdata = ndF, draws = TRUE)))
  expect_identical(dim(dr), c(4L, 6L))
  expect_identical(which(colSums(is.na(dr)) > 0), 6L)
  expect_false(anyNA(dr[, 1:5]))

  # Every row unseen: all NA, and still the classed warning.
  ndW <- d[1:3, ]; ndW$cover <- "water"
  rw <- .fbl_warnings(suppressMessages(predict(fit, newdata = ndW)))
  expect_true(all(is.na(rw$value)))
  expect_s3_class(rw$warnings[[1]], "spatialkit_new_levels")
})

test_that("a factor that only lists an extra level predicts as before", {
  skip_if_not_installed("brms")
  d   <- .fbl_pts()
  fit <- .fbl_fit(d)
  seen <- new.env()
  local_mocked_bindings(posterior_epred = .fbl_mock_epred(seen), .package = "brms")
  ref <- suppressMessages(predict(fit, newdata = d[1:6, ]))
  ndC <- d[1:6, ]
  ndC$cover <- factor(as.character(ndC$cover),
                      levels = c("crop", "forest", "urban", "water"))
  # brms refused this newdata although no row carries "water": six NAs.
  expect_no_warning(got <- suppressMessages(predict(fit, newdata = ndC)))
  expect_equal(got, ref)
})

test_that("a model fitted on a subset with an unused level predicts its own rows", {
  skip_if_not_installed("brms")
  d  <- .fbl_pts()
  d2 <- d[d$cover != "urban", ]          # levels() still lists "urban"
  expect_true("urban" %in% levels(d2$cover))
  seen <- new.env()
  local_mocked_bindings(posterior_epred = .fbl_mock_epred(seen), .package = "brms")
  # With and without the levels stored at fit time (a fit saved by an earlier
  # version has none and derives them from its training rows).
  for (stored in c(TRUE, FALSE)) {
    fit2 <- .fbl_fit(d2, store_levels = stored)
    # fitted() stopped with brms's "New factor levels are not allowed".
    fv <- fitted(fit2)
    expect_length(fv, nrow(d2))
    expect_true(all(is.finite(fv)))
    # predict() on the training rows returned NA for every one of them.
    expect_no_warning(p <- suppressMessages(predict(fit2, newdata = d2[1:5, ])))
    expect_true(all(is.finite(p)))
    expect_no_warning(pm <- suppressMessages(predict(fit2, summary = "median")))
    expect_true(all(is.finite(pm)))
    # Rows of the parent layer that carry the level it was not fitted with.
    r <- .fbl_warnings(suppressMessages(predict(fit2, newdata = d[1:6, ])))
    urban <- as.character(d$cover[1:6]) == "urban"
    expect_identical(is.na(r$value), urban)
    expect_length(r$warnings, 1L)
    expect_match(conditionMessage(r$warnings[[1]]),
                 sprintf("%d of 6 row(s)", sum(urban)), fixed = TRUE)
  }
})

test_that("predict_surface() does not depend on chunk_size when one cell has a new level", {
  skip_if_not_installed("brms")
  d   <- .fbl_pts()
  fit <- .fbl_fit(d)
  seen <- new.env()
  local_mocked_bindings(posterior_epred = .fbl_mock_epred(seen), .package = "brms")
  g <- d[1:12, "cover"]; g$cover <- as.character(g$cover); g$cover[11] <- "water"
  s4 <- suppressWarnings(suppressMessages(predict_surface(fit, grid = g, chunk_size = 4)))
  sA <- suppressWarnings(suppressMessages(predict_surface(fit, grid = g, chunk_size = 5000)))
  # Rows 9-12 were NA at chunk_size = 4 and all twelve at 5000.
  expect_identical(which(is.na(s4$.pred)), 11L)
  expect_equal(s4$.pred, sA$.pred)
})

test_that("any other failed draw is a classed warning, not only a log line", {
  skip_if_not_installed("brms")
  d   <- .fbl_pts()
  fit <- .fbl_fit(d)
  local_mocked_bindings(
    posterior_epred = function(object, newdata, ...) stop("injected draw failure 9041"),
    .package = "brms")
  r <- .fbl_warnings(suppressMessages(predict(fit, newdata = d[1:4, ])))
  expect_true(all(is.na(r$value)))
  expect_length(r$value, 4L)
  expect_length(r$warnings, 1L)
  expect_s3_class(r$warnings[[1]], "spatialkit_predict_failed")
  expect_match(conditionMessage(r$warnings[[1]]), "injected draw failure 9041")
  expect_match(conditionMessage(r$warnings[[1]]), "all 4 prediction(s) are NA",
               fixed = TRUE)
})


# ---------------------------------------------------------------------------
# Real sampler (opt-in)
# ---------------------------------------------------------------------------

test_that("a real fit predicts the known-level rows and returns NA for the new level", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan tests")
  skip_if_not_installed("brms")
  d <- .fbl_pts()
  fit_one <- function(dat) {
    fit <- NULL
    utils::capture.output(suppressWarnings(suppressMessages(
      fit <- fit_bayesian_spatial_model(
        dat, "z", "cover", gp_k = 5, chains = 1, iter = 400, cores = 1,
        seed = 1, compute_loo = FALSE, check_convergence = FALSE,
        backend = "rstan"))), type = "output")
    fit
  }
  fit <- fit_one(d)
  ref <- suppressMessages(predict(fit, newdata = d[1:6, ]))
  expect_true(all(is.finite(ref)))

  # (a) five finite values, one NA, one warning
  ndA <- d[1:6, ]; ndA$cover <- as.character(ndA$cover); ndA$cover[6] <- "water"
  r <- .fbl_warnings(suppressMessages(predict(fit, newdata = ndA)))
  expect_equal(r$value[1:5], ref[1:5])
  expect_true(is.na(r$value[6]))
  expect_length(r$warnings, 1L)
  expect_s3_class(r$warnings[[1]], "spatialkit_new_levels")

  # (b) an unused extra level predicts identically to the original factor
  ndC <- d[1:6, ]
  ndC$cover <- factor(as.character(ndC$cover),
                      levels = c("crop", "forest", "urban", "water"))
  expect_equal(suppressMessages(predict(fit, newdata = ndC)), ref)

  # (d) predict_surface() gives the same surface at any chunk_size
  g <- d[1:12, "cover"]; g$cover <- as.character(g$cover); g$cover[11] <- "water"
  s4 <- suppressWarnings(suppressMessages(predict_surface(fit, grid = g, chunk_size = 4)))
  sA <- suppressWarnings(suppressMessages(predict_surface(fit, grid = g, chunk_size = 5000)))
  expect_identical(which(is.na(s4$.pred)), 11L)
  expect_equal(s4$.pred, sA$.pred)

  # (c) a fit on a subset with an unused level: finite fitted() and predict()
  d2   <- d[d$cover != "urban", ]
  fit2 <- fit_one(d2)
  expect_true(all(is.finite(fitted(fit2))))
  expect_true(all(is.finite(suppressMessages(predict(fit2, newdata = d2[1:5, ])))))
})
