# BUG-154, BUG-41: a call that stops, or returns an NA for a reason that
# fits nothing, raises no warning about what a later step would have done.
#
# resolution_profile() read the variogram, and warned what `cp` and
# `reliability` take from it, before it checked `levels` and bounded the
# search of each scored row's nearest centre; a `levels` with no value the
# layer can hold then stopped the call after the warning.  The checks now
# run before the variogram is estimated or read.
#
# estimate_sac_range() warned that "the variogram is fitted to the RAW
# response instead" before the refusals of a constant response and of
# points with no extent, which fit no variogram.  The warning now follows
# them.  And where nlme::gls() cannot fit the trend at all (a factor with
# one level), the warning of the fallback to OLS said the REML fit "did not
# converge"; it now says it could not be fitted, with the design's error.

# The R warnings `expr` raises, with its value or its error.
fwr_run <- function(expr) {
  ws <- list()
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(w) {
      ws[[length(ws) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  list(value = v, warnings = ws)
}
fwr_quiet <- function(expr)
  logger::with_log_threshold(expr, threshold = logger::FATAL,
                             namespace = "spatialkit", index = 2)
fwr_texts <- function(r) vapply(r$warnings, conditionMessage, "")

# 300 points with a trend on `w`, and a variogram of the residuals on `w`:
# passed to a profile of the raw response, it draws the warning that `cp`
# and `reliability` take their nugget, sill and range from another variable.
fwr_layer <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(4)
      n <- 300
      xy <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000))
      D <- as.matrix(stats::dist(xy))
      xy$w <- xy$x / 100 + rnorm(n)
      xy$z <- as.numeric(t(chol(exp(-D / 150) + diag(0.3, n))) %*% rnorm(n)) + xy$w
      p <- sf::st_as_sf(xy, coords = c("x", "y"), crs = 32632)
      sac <- fwr_quiet(suppressWarnings(estimate_sac_range(p, "z", predictor_vars = "w")))
      cache <<- list(p = p, sac = sac)
    }
    cache
  }
})

test_that("a `levels` the layer cannot hold stops with no warning about the variogram", {
  skip_if_not_installed("gstat")
  L <- fwr_layer()
  # The warning the control raises: the variogram is of residuals and the
  # profile scores the raw response.
  ctrl <- fwr_quiet(fwr_run(resolution_profile(L$p, "z", sac = L$sac, levels = c(5, 10),
                                               estimate = NULL)))
  expect_s3_class(ctrl$value, "resolution_profile")
  expect_identical(ctrl$value$levels, c(5L, 10L))
  expect_true(any(grepl("so `cp` and `reliability` take their nugget, sill and range from a different variable",
                        fwr_texts(ctrl), fixed = TRUE)))
  # No value between 2 and 299, with that `sac`: the error alone.
  for (lv in list(1, 5000, c(0, 1, 300))) {
    r <- fwr_quiet(fwr_run(resolution_profile(L$p, "z", sac = L$sac, levels = lv,
                                              estimate = NULL)))
    expect_s3_class(r$value, "error")
    expect_identical(conditionMessage(r$value),
                     "resolution_profile(): no usable value in `levels`.")
    expect_length(r$warnings, 0L)
  }
  # The variogram estimated by the profile (one that reports no usable
  # range on this layer, with spatialkit_variogram_unusable): not estimated
  # at all now, and no warning.
  r <- fwr_quiet(fwr_run(resolution_profile(L$p, "z", levels = c(400L, 500L),
                                            estimate = NULL, seed = 1)))
  expect_s3_class(r$value, "error")
  expect_identical(conditionMessage(r$value),
                   "resolution_profile(): no usable value in `levels`.")
  expect_length(r$warnings, 0L)
})

test_that("the nearest-centre work limit stops a ladder of counts before any warning", {
  skip_if_not_installed("gstat")
  L <- fwr_layer()
  # The ladder of counts as built, with FNN: its finest level is the one the
  # guard reads.
  full <- fwr_quiet(suppressWarnings(resolution_profile(L$p, "z", sac = L$sac, n_levels = 4)))
  finest <- max(full$levels)
  no_fnn <- function(pkg, ...) if (identical(pkg, "FNN")) FALSE else requireNamespace(pkg, quietly = TRUE)
  withr::local_options(spatialkit.exact_pairs_max = 1000)
  local_mocked_bindings(.has_backend = no_fnn)
  # `levels` given, the ladder of counts, and a layer whose min_cell_n
  # leaves fewer than two cells (which warns "spatialkit_not_supported" on a
  # profile that runs): each stops at the guard with no warning.
  calls <- list(
    given  = list(fn = function() resolution_profile(L$p, "z", sac = L$sac, levels = c(5, 10)),
                  work = 300 * 10),
    ladder = list(fn = function() resolution_profile(L$p, "z", sac = L$sac, n_levels = 4),
                  work = 300 * finest),
    few    = list(fn = function() resolution_profile(L$p, "z", sac = L$sac, levels = c(5, 10),
                                                     min_cell_n = 200),
                  work = 300 * 10))
  for (nm in names(calls)) {
    r <- fwr_quiet(fwr_run(calls[[nm]]$fn()))
    expect_s3_class(r$value, "spatialkit_work_limit")
    expect_identical(r$value$work, calls[[nm]]$work, info = nm)
    expect_length(r$warnings, 0L)
  }
  # Under sample_n the ladder ends at what the subsample can fit, two scored
  # rows a cell: 20 cells on 40 points, below the 33 the 300 points allow at
  # min_cell_n = 9.  The guard reads that ceiling, so its work is 300 * 20
  # (an integer product, the same on every platform): a limit of 6000 lets
  # the profile run to 20 levels, and one below it stops with that work.
  sub <- function() resolution_profile(L$p, "z", sac = L$sac, sample_n = 40, n_levels = 4)
  r <- withr::with_options(list(spatialkit.exact_pairs_max = 6000),
                           fwr_quiet(fwr_run(sub())))
  expect_s3_class(r$value, "resolution_profile")
  expect_identical(max(r$value$levels), 20L)
  for (lim in c(5999, 5000)) {
    r <- withr::with_options(list(spatialkit.exact_pairs_max = lim),
                             fwr_quiet(fwr_run(sub())))
    expect_s3_class(r$value, "spatialkit_work_limit")
    expect_identical(r$value$work, 300 * 20, info = paste("limit", lim))
    expect_length(r$warnings, 0L)
  }
})

test_that("a `levels` of one short of the points is kept, and of the points dropped", {
  skip_if_not_installed("gstat")
  L <- fwr_layer()
  r <- fwr_quiet(fwr_run(resolution_profile(L$p, "z", sac = 300, levels = 299, estimate = NULL)))
  expect_s3_class(r$value, "resolution_profile")
  expect_identical(r$value$levels, 299L)
  r <- fwr_quiet(fwr_run(resolution_profile(L$p, "z", sac = 300, levels = 300, estimate = NULL)))
  expect_s3_class(r$value, "error")
  expect_identical(conditionMessage(r$value),
                   "resolution_profile(): no usable value in `levels`.")
  expect_length(r$warnings, 0L)
})

# 60 points, the response constant, the predictor missing on every row.
fwr_constant <- function() {
  set.seed(21)
  sf::st_as_sf(data.frame(x = runif(60, 0, 1000), y = runif(60, 0, 1000), z = 7,
                          w = NA_real_),
               coords = c("x", "y"), crs = 32615)
}
# 40 locations on a lattice, a structured response, the predictor missing on
# every row: the raw response is fitted, so the warning is true.
fwr_structured <- function() {
  g <- expand.grid(x = seq(0, 900, by = 100), y = seq(0, 300, by = 100))
  set.seed(1)
  sf::st_as_sf(data.frame(g, z = sin(g$x / 300) + cos(g$y / 200) + rnorm(40, sd = 0.1),
                          w = NA_real_),
               coords = c("x", "y"), crs = 32632)
}

test_that("a constant response or no extent after a failed detrending gives the NA and no warning", {
  skip_if_not_installed("gstat")
  p <- fwr_constant()
  for (dt in c("ols", "reml")) {
    r <- fwr_quiet(fwr_run(estimate_sac_range(p, "z", predictor_vars = "w", detrend = dt)))
    expect_true(is.na(r$value))
    expect_identical(attr(r$value, "rejected_reason"), "the response is constant")
    expect_length(r$warnings, 0L)
  }
  # lm() fails on a factor with one level: the same.
  pf <- p
  pf$w <- factor("a")
  r <- fwr_quiet(fwr_run(estimate_sac_range(pf, "z", predictor_vars = "w")))
  expect_true(is.na(r$value))
  expect_identical(attr(r$value, "rejected_reason"), "the response is constant")
  expect_length(r$warnings, 0L)
  # 40 distinct locations 1e-300 apart: no extent.
  set.seed(2)
  ext <- sf::st_as_sf(data.frame(x = (0:39) * 1e-300, y = 0, z = rnorm(40), w = NA_real_),
                      coords = c("x", "y"))
  r <- fwr_quiet(fwr_run(estimate_sac_range(ext, "z", predictor_vars = "w")))
  expect_true(is.na(r$value))
  expect_match(attr(r$value, "rejected_reason"), "^the points have no extent")
  expect_length(r$warnings, 0L)
})

test_that("a raw response that is fitted after the failed detrending keeps the warning", {
  skip_if_not_installed("gstat")
  p <- fwr_structured()
  r <- fwr_quiet(fwr_run(estimate_sac_range(p, "z", predictor_vars = "w")))
  expect_false(is.null(attr(r$value, "variogram_model")))
  expect_false(attr(r$value, "detrended"))
  expect_identical(fwr_texts(r), paste0(
    "estimate_sac_range(): the OLS detrending on w failed (none of the 40 row(s) with a ",
    "response has every predictor present and finite); the variogram is fitted to the RAW ",
    "response instead, which includes the trend."))
  pf <- p
  pf$w <- factor("a")
  r <- fwr_quiet(fwr_run(estimate_sac_range(pf, "z", predictor_vars = "w")))
  expect_false(is.null(attr(r$value, "variogram_model")))
  expect_false(attr(r$value, "detrended"))
  expect_length(r$warnings, 1L)
  expect_match(fwr_texts(r), paste0("^estimate_sac_range\\(\\): the OLS detrending on w failed ",
                                    "\\(.*\\); the variogram is fitted to the RAW response ",
                                    "instead, which includes the trend\\.$"))
})

test_that("a REML trend gls() cannot fit is said to be not fitted, not to not converge", {
  skip_if_not_installed("gstat")
  skip_if_not_installed("nlme")
  # A factor with one level: the trend's design cannot be built, so gls()
  # failed before any iteration.  The warning gives the design's error, and
  # the OLS fallback it names still runs (and fails on the same factor).
  reml_said <- function(r) grep("REML detrending", fwr_texts(r), value = TRUE, fixed = TRUE)
  pattern <- paste0("^estimate_sac_range\\(\\): the REML detrending on w could not be ",
                    "fitted \\(.+\\); falling back to OLS detrending, whose residual ",
                    "variogram is biased toward a shorter range \\(see \\?estimate_sac_range\\)\\.$")
  pf <- fwr_constant()
  pf$w <- factor("a")
  r <- fwr_quiet(fwr_run(estimate_sac_range(pf, "z", predictor_vars = "w", detrend = "reml")))
  expect_identical(attr(r$value, "rejected_reason"), "the response is constant")
  expect_length(r$warnings, 1L)
  expect_match(reml_said(r), pattern)
  expect_false(any(grepl("did not converge", fwr_texts(r), fixed = TRUE)))
  ps <- fwr_structured()
  ps$w <- factor("a")
  r <- fwr_quiet(fwr_run(estimate_sac_range(ps, "z", predictor_vars = "w", detrend = "reml")))
  expect_length(reml_said(r), 1L)
  expect_match(reml_said(r), pattern)
  expect_false(any(grepl("did not converge", fwr_texts(r), fixed = TRUE)))
  # .reml_trend() hands the error back; a trend it can fit is fitted.
  mf <- data.frame(.sac_x = rep(0:5, 5) * 100, .sac_y = rep(0:4, each = 6) * 100,
                   z = withr::with_seed(3, stats::rnorm(30)), w = factor("a"))
  u <- .reml_trend(mf, z ~ w, extent = 500)
  expect_named(u, "unfitted")
  expect_true(is.character(u$unfitted) && nzchar(u$unfitted))
})

test_that("the help says when there is no warning", {
  root <- testthat::test_path("..", "..")
  flat <- function(x)
    paste(unlist(lapply(x, function(el) if (is.character(el)) el else flat(el))), collapse = "")
  rd_of <- function(name) {
    file <- file.path(root, "man", paste0(name, ".Rd"))
    db <- if (file.exists(file)) tools::parse_Rd(file)
          else tools::Rd_db("spatialkit")[[paste0(name, ".Rd")]]
    gsub("[[:space:]]+", " ", flat(db))
  }
  says <- function(text, x) expect_true(grepl(x, text, fixed = TRUE), info = x)
  says(rd_of("resolution_profile"),
       paste0("When none is left the call stops with an error, before the variogram is ",
              "estimated or read, so it raises no warning about the variogram."))
  says(rd_of("estimate_sac_range"),
       paste0("When it is not fitted (no nlme), cannot be fitted (the design of the trend ",
              "cannot be built, as on a factor with one level) or does not converge, the ",
              "\"ols\" path runs instead with an R warning that says which."))
  says(rd_of("estimate_sac_range"),
       paste0("Where nothing is fitted then either (fewer than 30 distinct locations, a ",
              "constant response, or points with no extent), the result is the NA with its ",
              "reason, and there is no such warning."))
})
