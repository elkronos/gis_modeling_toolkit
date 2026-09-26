# tests/testthat/test-review3-S6-cv-eval.R
# ---------------------------------------------------------------------------
# Regressions from the third review of the CV runners and the model
# evaluation functions.  Each test names the finding it closes.  The models
# are helper-lmfit.R's lm fit or a small custom subclass, so no optional
# backend is needed except where skipped.
# ---------------------------------------------------------------------------

.r3_pts <- function(n = 80, seed = 1, extent = 1000) {
  set.seed(seed)
  d <- sf::st_as_sf(
    data.frame(x = runif(n, 0, extent), y = runif(n, 0, extent), w = rnorm(n)),
    coords = c("x", "y"), crs = 32632)
  # A trend w cannot explain, so the residuals carry spatial structure.
  d$z <- 5 + 0.004 * sf::st_coordinates(d)[, 1] + 2 * d$w + rnorm(n, 0, 0.5)
  d
}

.r3_lm <- function(tr) lm_spatial_fit(tr, "z", "w")

.r3_warnings <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w)
}


# ---------------------------------------------------------------------------
# S6-CV-EVAL-1: residual Moran's I for a custom fit without residuals()
# ---------------------------------------------------------------------------

test_that("residual_morans_i() scores a fit with fitted() but no residuals() on y - fitted", {
  d <- .r3_pts(80, seed = 1)
  # The documented minimum: predict() and fitted(), no residuals() method.
  registerS3method("fitted", "r3_fitonly",
                   function(object, ...) as.numeric(stats::fitted(object$engine)))
  registerS3method("predict", "r3_fitonly", function(object, newdata = NULL, ...)
    as.numeric(stats::predict(object$engine, sf::st_drop_geometry(newdata))))
  bare <- new_spatial_fit(subclass = "r3_fitonly",
                          engine = stats::lm(z ~ w, sf::st_drop_geometry(d)),
                          formula = z ~ w, response_var = "z",
                          predictor_vars = "w", data_sf = d)
  full <- lm_spatial_fit(d, "z", "w")            # has a residuals() method
  expect_no_warning(mi <- residual_morans_i(bare))
  ref <- residual_morans_i(full)
  expect_s3_class(mi, "morans_i")
  expect_equal(mi$observed, ref$observed)
  expect_equal(mi$z, ref$z)
  expect_gt(mi$z, 3)                             # the trend left behind

  # compare_models() now reports it instead of all-NA columns.
  cmp <- .r3_warnings(compare_models(list(bare = bare)))
  expect_false(any(grepl("residuals\\(\\) returned NULL", cmp$warnings)))
  expect_equal(cmp$value$resid_morans_I, ref$observed)
  expect_false(is.na(cmp$value$resid_morans_null))
})

test_that("residual_morans_i() still says why when neither residuals() nor fitted() exists", {
  d <- .r3_pts(60, seed = 2)
  fit <- new_spatial_fit(subclass = "r3_nothing", engine = NULL,
                         formula = z ~ w, response_var = "z",
                         predictor_vars = "w", data_sf = d)
  expect_warning(out <- residual_morans_i(fit),
                 "has no residuals\\(\\) method, and the observed response minus fitted\\(\\) could not be formed either: fitted\\(\\) returned NULL")
  expect_null(out)
})


# ---------------------------------------------------------------------------
# S6-CV-EVAL-2: a cv_*() result's folds keep their fold_id when re-used
# ---------------------------------------------------------------------------

test_that("re-feeding a cv result's folds keeps the labels after a dropped fold", {
  d <- .r3_pts(80, seed = 3)
  f <- make_folds(d, k = 5, method = "random_kfold", seed = 1)
  holed <- d
  holed$z[f$folds[[3]]$test] <- NA         # fold 3 loses every test row
  cv1 <- suppressWarnings(cv_spatial(holed, "z", "w", fit_fn = .r3_lm, folds = f))
  expect_identical(cv1$fold_metrics$fold, c(1L, 2L, 4L, 5L))
  cv2 <- suppressWarnings(cv_spatial(holed, "z", "w", fit_fn = .r3_lm,
                                     folds = cv1$folds))
  expect_identical(cv2$fold_metrics$fold, c(1L, 2L, 4L, 5L))
  expect_equal(cv2$fold_metrics$RMSE, cv1$fold_metrics$RMSE)
  expect_identical(sort(unique(cv2$predictions$fold)), c(1L, 2L, 4L, 5L))
  expect_identical(cv2$fold_status$fold, c(1L, 2L, 4L, 5L))
  expect_identical(vapply(cv2$folds, `[[`, integer(1), "fold_id"), c(1L, 2L, 4L, 5L))
  # ... which is what fold_separation() calls them.
  fs <- fold_separation(cv1$folds, holed)
  expect_identical(as.integer(fs$fold), cv2$fold_metrics$fold)
})

test_that("splits without a usable, distinct fold_id are numbered by position", {
  d <- .r3_pts(60, seed = 4)
  lab <- rep(1:3, each = 20)
  sp <- lapply(1:3, function(j) list(train = which(lab != j), test = which(lab == j),
                                     fold_id = 7L))            # all the same
  cv <- cv_spatial(d, "z", "w", fit_fn = .r3_lm, folds = sp)
  expect_identical(cv$fold_metrics$fold, 1:3)
  sp[[1]]$fold_id <- NULL; sp[[2]]$fold_id <- 8L; sp[[3]]$fold_id <- 9L
  cv <- cv_spatial(d, "z", "w", fit_fn = .r3_lm, folds = sp)
  expect_identical(cv$fold_metrics$fold, 1:3)
  # A carried id also names the fold in the train/test overlap error.
  bad <- list(list(train = 1:40, test = 35:60, fold_id = 4L),
              list(train = 21:60, test = 1:20, fold_id = 6L))
  expect_error(cv_spatial(d, "z", "w", fit_fn = .r3_lm, folds = bad),
               "fold 4 has 6 row ID\\(s\\) in BOTH")
})


# ---------------------------------------------------------------------------
# S6-CV-EVAL-3: a ..per_row column may not reuse a predictions column
# ---------------------------------------------------------------------------

test_that("fold_info_fn's ..per_row may not reuse a column predictions has", {
  d <- .r3_pts(60, seed = 5)
  lab <- rep(1:3, each = 20)
  for (col in c("yhat", "fold", "y", "..row_id", "y_train_mean")) {
    info <- function(fit, test_sf, y, yhat) {
      pr <- data.frame(v = yhat + 0); names(pr) <- col; list(..per_row = pr)
    }
    expect_error(cv_spatial(d, "z", "w", fit_fn = .r3_lm, folds = lab,
                            fold_info_fn = info),
                 paste0("`..per_row` has columns that `predictions` already has: ",
                        col, "."), fixed = TRUE)
  }
  dup <- function(fit, test_sf, y, yhat)
    list(..per_row = data.frame(a = y, a = yhat, check.names = FALSE))
  expect_error(cv_spatial(d, "z", "w", fit_fn = .r3_lm, folds = lab,
                          fold_info_fn = dup), "duplicated column names: a")
  # A new name is spliced in and `overall` is pooled as usual.
  ok <- function(fit, test_sf, y, yhat) list(..per_row = data.frame(ae = abs(y - yhat)))
  cv <- cv_spatial(d, "z", "w", fit_fn = .r3_lm, folds = lab, fold_info_fn = ok)
  expect_equal(cv$predictions$ae, abs(cv$predictions$y - cv$predictions$yhat))
  expect_identical(cv$overall$n_pred, 60L)
})

test_that("a ..per_row name clash stops a parallel run too", {
  skip_on_os("windows")
  skip_on_cran()
  skip_if(isTRUE(parallel::detectCores(logical = TRUE) < 2L))
  d <- .r3_pts(60, seed = 6)
  info <- function(fit, test_sf, y, yhat) list(..per_row = data.frame(yhat = yhat))
  expect_error(suppressMessages(
    cv_spatial(d, "z", "w", fit_fn = .r3_lm, folds = rep(1:2, each = 30),
               fold_info_fn = info, parallel = 2)),
    "already has: yhat.*raised in a parallel worker")
})


# ---------------------------------------------------------------------------
# S6-CV-EVAL-4: fold_info_fn returning a named vector, or something else
# ---------------------------------------------------------------------------

test_that("fold_info_fn may return a named vector; anything else not a list is an error", {
  d <- .r3_pts(60, seed = 7)
  lab <- rep(1:3, each = 20)
  vec <- function(fit, test_sf, y, yhat)
    c(slope = unname(stats::coef(fit$engine)[2]), n_te = length(y))
  cv <- cv_spatial(d, "z", "w", fit_fn = .r3_lm, folds = lab, fold_info_fn = vec)
  expect_true(all(c("slope", "n_te") %in% names(cv$fold_metrics)))
  expect_equal(cv$fold_metrics$n_te, c(20, 20, 20))
  expect_true(all(is.finite(cv$fold_metrics$slope)))

  # NULL is still "no extras".
  cv0 <- cv_spatial(d, "z", "w", fit_fn = .r3_lm, folds = lab,
                    fold_info_fn = function(...) NULL)
  expect_identical(cv0$fold_status$status, rep("ok", 3))

  expect_error(cv_spatial(d, "z", "w", fit_fn = .r3_lm, folds = lab,
                          fold_info_fn = function(...) new.env()),
               "must return a named list \\(or a named vector\\).*class environment")
  expect_error(cv_spatial(d, "z", "w", fit_fn = .r3_lm, folds = lab,
                          fold_info_fn = function(...) 3),
               "every element is named")
})


# ---------------------------------------------------------------------------
# S6-CV-EVAL-5: folds built on a pointized copy of a polygon layer
# ---------------------------------------------------------------------------

test_that("folds built on a pointized copy are refused for the geometry, not the data", {
  d <- .r3_pts(60, seed = 8)
  xy <- sf::st_coordinates(d)
  ell <- function(x, y) sf::st_polygon(list(rbind(
    c(x, y), c(x + 60, y), c(x + 60, y + 15), c(x + 15, y + 15),
    c(x + 15, y + 60), c(x, y + 60), c(x, y))))
  poly <- sf::st_sf(sf::st_drop_geometry(d),
                    geometry = sf::st_sfc(lapply(seq_len(nrow(xy)), function(i)
                      ell(xy[i, 1], xy[i, 2])), crs = 32632))
  pp <- coerce_to_points(poly, "auto")           # st_point_on_surface()
  f_pts <- make_folds(pp, k = 3, method = "random_kfold", seed = 1)
  expect_true(f_pts$params$row_probe$points)
  err <- tryCatch(cv_spatial(poly, "z", "w", fit_fn = .r3_lm, folds = f_pts),
                  error = conditionMessage)
  expect_match(err, "built on POINT geometry and this data has non-POINT geometry")
  expect_match(err, "make_folds\\(\\) on the layer passed here")
  expect_no_match(err, "different data")

  # The remedy the message names works.
  f_poly <- make_folds(poly, k = 3, method = "random_kfold", seed = 1)
  expect_false(f_poly$params$row_probe$points)
  expect_no_error(cv_spatial(poly, "z", "w", fit_fn = .r3_lm, folds = f_poly))

  # A probe from before the field keeps the old refusal.
  legacy <- f_pts
  legacy$params$row_probe$points <- NULL
  expect_error(cv_spatial(poly, "z", "w", fit_fn = .r3_lm, folds = legacy),
               "built from different data")
})


# ---------------------------------------------------------------------------
# S6-CV-EVAL-7: an information criterion across different responses
# ---------------------------------------------------------------------------

test_that("compare_models() names a different response, not different rows", {
  d <- .r3_pts(70, seed = 9)
  d$lz <- log(d$z - min(d$z) + 1)
  stub <- function(data, resp, looic) {
    fit <- lm_spatial_fit(data, resp, "w")
    fit$info$looic <- looic
    class(fit) <- c("lmsurf_fit", "bayesian_fit", "spatial_fit")
    fit
  }
  a <- stub(d, "z", 54.8); b <- stub(d, "lz", 12.1)
  res <- .r3_warnings(compare_models(list(raw = a, logged = b)))
  expect_true(any(grepl(paste0("LOOIC .* fitted to the same rows but to different ",
                               "responses \\(raw: z, logged: lz\\)"), res$warnings)))
  expect_false(any(grepl("Refit them on the same rows", res$warnings)))
  expect_true(all(is.na(res$value$LOOIC)))

  # The same column name holding other values says that instead.
  d2 <- d; d2$z <- 2 * d2$z
  res2 <- .r3_warnings(compare_models(list(raw = a, doubled = stub(d2, "z", 60))))
  expect_true(any(grepl("the values of 'z' differ between raw, doubled",
                        res2$warnings)))
  # Different rows keep their own message.
  res3 <- .r3_warnings(compare_models(list(raw = a, sub = stub(d[1:50, ], "z", 30))))
  expect_true(any(grepl("different rows \\(raw: n = 70, sub: n = 50\\)",
                        res3$warnings)))
})


# ---------------------------------------------------------------------------
# S6-CV-EVAL-8: nothing to evaluate is an error, naming `fits`
# ---------------------------------------------------------------------------

test_that("evaluate_insample() errors when no element is a spatial_fit", {
  expect_error(suppressWarnings(
    evaluate_insample(list(a = 1, b = stats::lm(dist ~ speed, datasets::cars)))),
    "evaluate_insample\\(\\): no element of `fits` is a spatial_fit")
  expect_error(compare_models(list(a = 1, b = "x")),
               "compare_models\\(\\): no element of `fits` is a spatial_fit")
  # A mixed list still skips the stranger and scores the fit.
  d <- .r3_pts(40, seed = 10)
  out <- evaluate_insample(list(fit = lm_spatial_fit(d, "z", "w"), junk = 1))
  expect_s3_class(out, "data.frame")
  expect_identical(out$model, "fit")
})


# ---------------------------------------------------------------------------
# S6-CV-EVAL-9: cv_rf(parallel = more than the machine has) says so once
# ---------------------------------------------------------------------------

test_that("cv_rf() prints the worker cap message once", {
  skip_if_not_installed("ranger")
  skip_on_os("windows")
  skip_on_cran()                                 # forks the folds
  d <- .r3_pts(45, seed = 11)
  n_machine <- parallel::detectCores(logical = TRUE)
  skip_if(is.na(n_machine))
  msgs <- character(0)
  withCallingHandlers(
    cv_rf(d, "z", "w", folds = rep(1:3, each = 15), num_trees = 60, seed = 1,
          parallel = n_machine + 5L),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
    })
  expect_identical(sum(grepl("workers requested on a machine with", msgs)), 1L)
})
