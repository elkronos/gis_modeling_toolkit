# tests/testthat/test-audit-pass6-low.R
# ---------------------------------------------------------------------------
# The sixth pass's Low list: the items that were code rather than prose.
# Each test names the note it closes.
# ---------------------------------------------------------------------------

.p6l_pts <- function(n = 60, seed = 1) {
  set.seed(seed)
  sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), w = rnorm(n)),
    coords = c("x", "y"), crs = 32632)
}


# --------------------------------------------------------------------------
# C (Low): cv_bayes()$predictions$yhat_sd was an unconditional NA placeholder
# --------------------------------------------------------------------------

test_that("cv_bayes() fills yhat_sd from the posterior predictive draws", {
  d <- .p6l_pts(n = 60, seed = 2)
  d$z <- 2 * d$w + rnorm(60, 0, 0.2)
  f <- make_folds(d, k = 3, method = "random_kfold", seed = 1)
  # A stand-in backend whose draws have a known spread: three draws at
  # mu - 1, mu, mu + 1, so the per-row SD is exactly 1.
  registerS3method("predict", "p6l_sd_probe",
                   function(object, newdata = NULL, draws = FALSE, ...) {
                     mu <- predict.lmsurf_fit(object, newdata = newdata)
                     if (!isTRUE(draws)) return(mu)
                     rbind(mu - 1, mu, mu + 1)
                   })
  local_mocked_bindings(
    fit_bayesian_spatial_model = function(data_sf, response_var, predictor_vars,
                                          ..., seed = 123) {
      fit <- lm_spatial_fit(data_sf, response_var, predictor_vars)
      class(fit) <- c("p6l_sd_probe", class(fit))
      fit
    },
    .package = "spatialkit")
  cv <- suppressWarnings(cv_bayes(d, "z", "w", folds = f, seed = 1))
  expect_true("yhat_sd" %in% names(cv$predictions))
  expect_equal(nrow(cv$predictions), 60L)
  expect_equal(cv$predictions$yhat_sd, rep(1, 60))
  # The per-fold extras are unaffected: coverage at 50/80/95% and CRPS are
  # still per-fold scalars, and ..per_row never leaks into fold_metrics.
  expect_false("..per_row" %in% names(cv$fold_metrics))
  expect_true(all(c("coverage_50", "coverage_80", "coverage_95", "CRPS") %in%
                    names(cv$fold_metrics)))
  # Without the draws the column is present and NA, as documented.
  cv0 <- suppressWarnings(cv_bayes(d, "z", "w", folds = f, seed = 1,
                                   compute_pred_intervals = FALSE))
  expect_true(all(is.na(cv0$predictions$yhat_sd)))
  expect_equal(nrow(cv0$predictions), 60L)
})


# --------------------------------------------------------------------------
# E (Low): fit_gwr_model() on n <= p + 1 warned three times and returned an
# all-NA fit; now one error, before any backend work
# --------------------------------------------------------------------------

test_that("fit_gwr_model() refuses n <= p + 1 observations with a single error", {
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sp")
  d <- .p6l_pts(n = 3, seed = 3)
  d$z <- c(1, 2, 4)
  # n = 2, p = 1: two parameters, two rows.
  expect_error(
    expect_no_warning(fit_gwr_model(d[1:2, ], "z", "w", bandwidth = 30)),
    "^fit_gwr_model\\(\\): 2 observation\\(s\\) for 2 parameters")
  # n = 3 clears the guard (p + 2) and reaches the backend, however it fares.
  res <- tryCatch(suppressWarnings(fit_gwr_model(d, "z", "w", bandwidth = 30)),
                  error = function(e) conditionMessage(e))
  if (is.character(res))
    expect_no_match(res, "observation\\(s\\) for 2 parameters")
  else
    expect_s3_class(res, "gwr_fit")
})


# --------------------------------------------------------------------------
# F (Low): the grid cache never evicted and its key ignored the package version
# --------------------------------------------------------------------------

test_that("create_grid_polygons_cached() evicts the oldest grid past max_entries", {
  bnd <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)))), crs = 32632))
  e <- new.env()
  n_grids <- function() sum(startsWith(ls(e, all.names = TRUE), "spatialkit_grid::"))
  g4  <- create_grid_polygons_cached(bnd, target_cells = 4,  cache_env = e, max_entries = 2)
  g9  <- create_grid_polygons_cached(bnd, target_cells = 9,  cache_env = e, max_entries = 2)
  expect_equal(n_grids(), 2L)
  g16 <- create_grid_polygons_cached(bnd, target_cells = 16, cache_env = e, max_entries = 2)
  expect_equal(n_grids(), 2L)                       # capped
  k4 <- spatialkit:::.cache_key(bnd, "square", 4)
  k9 <- spatialkit:::.cache_key(bnd, "square", 9)
  expect_false(exists(k4, envir = e, inherits = FALSE))   # the oldest went
  expect_true(exists(k9, envir = e, inherits = FALSE))
  # Eviction is transparent: the evicted grid is rebuilt identically.
  g4b <- create_grid_polygons_cached(bnd, target_cells = 4, cache_env = e, max_entries = 2)
  expect_equal(g4b, g4)
  expect_equal(n_grids(), 2L)
  expect_false(exists(k9, envir = e, inherits = FALSE))   # 9 was the oldest by then
  # A hit does not count as an insertion and does not evict anything.
  create_grid_polygons_cached(bnd, target_cells = 4, cache_env = e, max_entries = 2)
  expect_equal(n_grids(), 2L)
  # Nothing but the grids is written into the caller's environment, and
  # clear_grid_cache() forgets the insertion order along with them.
  expect_length(ls(e, all.names = TRUE), 2L)
  expect_equal(clear_grid_cache(e), 2L)
  expect_length(ls(e, all.names = TRUE), 0L)
  expect_length(spatialkit:::.cache_order(e), 0L)
  expect_error(create_grid_polygons_cached(bnd, target_cells = 4, cache_env = e,
                                           max_entries = 0),
               "`max_entries` must be a single number >= 1")
})

test_that("the grid cache key carries the package version", {
  bnd <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)))), crs = 32632))
  k1 <- spatialkit:::.cache_key(bnd, "square", 9, version = "1.0.0")
  k2 <- spatialkit:::.cache_key(bnd, "square", 9, version = "2.0.0")
  expect_false(identical(k1, k2))
  expect_identical(k1, spatialkit:::.cache_key(bnd, "square", 9, version = "1.0.0"))
  expect_match(spatialkit:::.spatialkit_version(), "^[0-9]+\\.[0-9]+")
})


# --------------------------------------------------------------------------
# F (Low): spatialkit_quiet(FALSE) restored the default, not the prior level,
# and the value it returned could not be passed back
# --------------------------------------------------------------------------

test_that("spatialkit_quiet() can put back exactly the level that was in force", {
  before <- logger::log_threshold(namespace = "spatialkit", index = 2)
  on.exit(logger::log_threshold(before, namespace = "spatialkit", index = 2),
          add = TRUE)
  logger::log_threshold(logger::ERROR, namespace = "spatialkit", index = 2)
  old <- spatialkit_quiet()
  expect_identical(old, logger::ERROR)
  expect_identical(logger::log_threshold(namespace = "spatialkit", index = 2),
                   logger::FATAL)
  # The returned value restores ERROR, where FALSE would have given WARN.
  spatialkit_quiet(old)
  expect_identical(logger::log_threshold(namespace = "spatialkit", index = 2),
                   logger::ERROR)
  spatialkit_quiet(FALSE)
  expect_identical(logger::log_threshold(namespace = "spatialkit", index = 2),
                   logger::WARN)
  # A level object is accepted directly; anything else is refused by name.
  spatialkit_quiet(logger::INFO)
  expect_identical(logger::log_threshold(namespace = "spatialkit", index = 2),
                   logger::INFO)
  expect_error(spatialkit_quiet("quiet"), "must be TRUE, FALSE or a logger threshold")
  expect_error(spatialkit_quiet(NA), "must be TRUE, FALSE or a logger threshold")
})


# --------------------------------------------------------------------------
# H (Low): no .onUnload, so the logger namespace kept pointing at the session
# temp file after unloadNamespace()
# --------------------------------------------------------------------------

test_that(".onUnload() disarms both logger appenders and .onLoad() re-arms them", {
  # Restore by re-running .onLoad(), never by putting back what the getter
  # returns: logger's getter hands back the appender's GENERATOR expression
  # (`appender_file(log_path)`), and storing that as an appender breaks every
  # later log call with "object 'log_path' not found".
  on.exit(spatialkit:::.onLoad(tempdir(), "spatialkit"), add = TRUE)
  spatialkit:::.onUnload(tempdir())
  a1 <- logger::log_appender(namespace = "spatialkit", index = 1)
  a2 <- logger::log_appender(namespace = "spatialkit", index = 2)
  expect_true(is.function(a1))
  expect_true(is.function(a2))
  expect_null(a1(c("a line")))                  # a no-op, not a file write
  expect_null(a2(c("a line")))
  expect_silent(logger::log_warn("nothing should print", namespace = "spatialkit"))
  # Re-loading re-registers the file trace (a generator-backed appender, which
  # the getter reports as its call) and the console echo at WARN.
  spatialkit:::.onLoad(tempdir(), "spatialkit")
  expect_true(is.call(logger::log_appender(namespace = "spatialkit", index = 1)))
  expect_identical(logger::log_threshold(namespace = "spatialkit", index = 2),
                   logger::WARN)
})


# --------------------------------------------------------------------------
# Fold shapes: a split list without `train`/`test` names is refused by name.
# area_of_applicability() always did; the cv_*() path read f$train / f$test,
# got NULL, built empty folds, and blamed "folds built on a different or
# subsetted layer".  blockCV::cv_spatial()$folds_list is the realistic way in.
# --------------------------------------------------------------------------

test_that("cv_*() refuse a positional fold list and name the real problem", {
  d <- .p6l_pts(n = 120, seed = 11)
  d$z <- 2 * d$w + rnorm(120, 0, 0.3)
  ids <- as.integer(cut(sf::st_coordinates(d)[, 1], breaks = 4, labels = FALSE))
  # The shape blockCV's $folds_list has: two UNNAMED vectors, train first.
  positional <- lapply(1:4, function(k)
    list(which(ids != k), which(ids == k)))

  msg <- paste0("^cross-validation: 4 of 4 fold\\(s\\) carry no `train`/`test` ",
                "element \\(e\\.g\\. fold 1, 2, 3\\)")
  expect_error(cv_spatial(d, "z", "w", fit_fn = lm_spatial_fit,
                          folds = positional), msg)
  expect_error(cv_spatial(d, "z", "w", fit_fn = lm_spatial_fit,
                          folds = positional), "folds_ids", fixed = TRUE)
  # The hint is specific to the positional shape, not printed for every refusal.
  half_named <- positional
  half_named[[2]] <- list(train = which(ids != 2), test = which(ids == 2))
  e <- tryCatch(cv_spatial(d, "z", "w", fit_fn = lm_spatial_fit,
                           folds = half_named), error = conditionMessage)
  expect_match(e, "^cross-validation: 3 of 4 fold\\(s\\) carry no")
  expect_false(grepl("folds_ids", e, fixed = TRUE))
  # Not a list of splits at all.
  expect_error(cv_spatial(d, "z", "w", fit_fn = lm_spatial_fit, folds = list()),
               "is not a recognised fold object")
  # The wording matches area_of_applicability()'s three-shape sentence.
  expect_error(area_of_applicability(d, train_sf = d, predictor_vars = "w",
                                     folds = positional),
               "must be a make_folds\\(\\) result, a list of train/test splits")
})

test_that("a fold label vector -- blockCV's $folds_ids shape -- runs end to end", {
  skip_if_not_installed("ranger")
  d <- .p6l_pts(n = 150, seed = 12)
  d$z <- 2 * d$w + rnorm(150, 0, 0.3)
  # One integer label per row, exactly what blockCV::cv_spatial() returns as
  # $folds_ids (and CAST's fold columns look the same).
  folds_ids <- as.integer(cut(sf::st_coordinates(d)[, 1], breaks = 5,
                              labels = FALSE))
  cv <- cv_rf(d, "z", "w", folds = folds_ids, num_trees = 50, seed = 1)
  expect_equal(cv$n_folds_attempted, 5L)
  expect_equal(cv$n_folds_succeeded, 5L)
  expect_equal(cv$overall$n_pred, 150L)          # every row scored exactly once
  expect_setequal(cv$fold_metrics$fold, 1:5)
  expect_setequal(cv$predictions$..row_id, seq_len(150))
  # Identical to spelling the same splits out by name.
  named <- lapply(1:5, function(k)
    list(train = which(folds_ids != k), test = which(folds_ids == k)))
  cv2 <- cv_rf(d, "z", "w", folds = named, num_trees = 50, seed = 1)
  expect_equal(cv2$overall, cv$overall)
  # And the same vector drives area_of_applicability()'s fold-aware threshold.
  aoa <- area_of_applicability(d, train_sf = d, predictor_vars = "w",
                               folds = folds_ids)
  expect_true(is.finite(aoa$threshold))
  # A make_folds() result still works, unchanged.
  f <- make_folds(d, k = 4, method = "random_kfold", seed = 1)
  expect_equal(cv_rf(d, "z", "w", folds = f, num_trees = 50,
                     seed = 1)$n_folds_succeeded, 4L)
})
