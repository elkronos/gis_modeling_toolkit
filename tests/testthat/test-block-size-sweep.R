# tests/testthat/test-block-size-sweep.R
# ---------------------------------------------------------------------------
# cv_block_size_sweep(): the same cross-validation at a ladder of block
# sizes, with the random-fold reference and the autocorrelation range, so
# the leakage of small blocks is a curve rather than an argument.  The lm
# stand-in keeps the sweep fast and free of optional backends.
# ---------------------------------------------------------------------------

sweep_points <- function(n = 160, seed = 21) {
  set.seed(seed)
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- as.matrix(stats::dist(cbind(x, y)))
  field <- as.numeric(t(chol(exp(-d / 60) + diag(1e-6, n))) %*% rnorm(n))
  pts <- sf::st_as_sf(data.frame(x = x, y = y, a = rnorm(n)), coords = c("x", "y"), crs = 32632)
  pts$z <- field + 0.5 * pts$a + rnorm(n, 0, 0.2)
  pts
}
sweep_fit <- function(train_sf) lm_spatial_fit(train_sf, "z", "a")

test_that("the sweep runs one k-fold CV per size plus the random reference", {
  pts <- sweep_points()
  sw <- cv_block_size_sweep(pts, "z", "a", fit_fn = sweep_fit, k = 3, n_sizes = 3,
                            sac = NA, quiet = TRUE)
  expect_s3_class(sw, "block_size_sweep")
  expect_s3_class(sw, "data.frame")
  expect_equal(nrow(sw), 4L)
  expect_equal(sw$method, c("random_kfold", rep("block_kfold", 3L)))
  expect_true(is.na(sw$block_size[1L]))
  expect_equal(sw$block_size[-1L], sort(sw$block_size[-1L]))
  expect_true(all(sw$k == 3L))
  expect_true(all(sw$n_folds_succeeded == 3L))
  expect_true(all(is.finite(sw$value)))
  expect_true(all(sw$fold_min <= sw$value & sw$value <= sw$fold_max))
  expect_true(all(sw$blocks_used[-1L] >= 3L))
  expect_equal(attr(sw, "metric"), "RMSE")
  expect_equal(attr(sw, "n_fits"), 12L)
  expect_equal(attr(sw, "crs"), "EPSG:32632")
  expect_true(is.na(attr(sw, "sac_range")))
  expect_length(attr(sw, "results"), 4L)
  # Each row is exactly what cv_spatial() reports on those folds.
  res <- attr(sw, "results")
  expect_equal(sw$value[1L], res$random_kfold$overall$RMSE)
  expect_equal(sw$value[2L], res[[2L]]$overall$RMSE)
  expect_output(print(sw), "against block size")
})

test_that("the ladder respects the fit budget and the k-block floor", {
  pts <- sweep_points()
  expect_error(cv_block_size_sweep(pts, "z", "a", fit_fn = sweep_fit, k = 4, n_sizes = 20,
                                   quiet = TRUE),
               "above the budget max_fits = 60")
  expect_error(cv_block_size_sweep(pts, "z", "a", fit_fn = sweep_fit, k = 4, n_sizes = 20,
                                   quiet = TRUE),
               "20 block sizes x 4 folds \\+ 4 for the random reference = 84")
  # With k = 5 the two-by-two grid at the top of the ladder holds too few
  # blocks and is dropped before the budget is counted.
  expect_error(cv_block_size_sweep(pts, "z", "a", fit_fn = sweep_fit, k = 5, n_sizes = 20,
                                   quiet = TRUE),
               "17 block sizes x 5 folds")
  # Sizes whose grid holds fewer than k blocks are dropped (logged), and a
  # ladder with nothing left is an error naming the extent.
  lines <- capture_spatialkit_log(
    sw <- cv_block_size_sweep(pts, "z", "a", fit_fn = sweep_fit, k = 4,
                              block_sizes = c(150, 600), include_random = FALSE,
                              sac = NA, quiet = TRUE),
    level = logger::INFO)
  expect_equal(nrow(sw), 1L)
  expect_equal(sw$block_size, 150)
  expect_true(log_has(lines, "dropping 1 block size"))
  expect_error(cv_block_size_sweep(pts, "z", "a", fit_fn = sweep_fit, k = 4,
                                   block_sizes = 600, quiet = TRUE),
               "no block size leaves at least k = 4 blocks")
  # Arguments the sweep sets cannot arrive through `...`.
  expect_error(cv_block_size_sweep(pts, "z", "a", fit_fn = sweep_fit, k = 3, n_sizes = 2,
                                   folds = list(), quiet = TRUE),
               "cannot be passed through `...`")
  expect_error(cv_block_size_sweep(pts, "z", "a", fit_fn = sweep_fit, k = 3, n_sizes = 2,
                                   metric = "nope", sac = NA, quiet = TRUE),
               "'nope' is not a column of cv_spatial\\(\\)\\$overall")
  expect_error(cv_block_size_sweep(pts, "z", "a", fit_fn = "not a function"),
               "`fit_fn` must be a function")
  expect_error(cv_block_size_sweep(pts, "z", "a", fit_fn = sweep_fit, k = 1), "`k` must be")
})

test_that("a supplied or estimated range is carried, and the plot draws the curve", {
  skip_if_not_installed("ggplot2")
  pts <- sweep_points()
  sac_fake <- structure(120, class = "sac_range")
  sw <- cv_block_size_sweep(pts, "z", "a", fit_fn = sweep_fit, k = 3, n_sizes = 3,
                            sac = sac_fake, quiet = TRUE)
  expect_equal(attr(sw, "sac_range"), 120)
  p <- plot(sw)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  geoms <- vapply(p$layers, function(l) class(l$geom)[1L], character(1))
  expect_true(all(c("GeomRibbon", "GeomLine", "GeomPoint", "GeomHline", "GeomVline") %in% geoms))
  expect_equal(b$data[[which(geoms == "GeomVline")]]$xintercept, log10(120))
  expect_equal(b$data[[which(geoms == "GeomHline")]]$yintercept, sw$value[1L])
  expect_equal(nrow(b$data[[which(geoms == "GeomPoint")]]), 3L)
  expect_match(p$labels$subtitle, "range 120")
  expect_match(p$labels$caption, "random folds, RMSE = ")
  # No reference and no range: neither marker, and the subtitle says so.
  sw2 <- cv_block_size_sweep(pts, "z", "a", fit_fn = sweep_fit, k = 3, n_sizes = 2,
                             include_random = FALSE, sac = NA, quiet = TRUE)
  p2 <- plot(sw2)
  geoms2 <- vapply(p2$layers, function(l) class(l$geom)[1L], character(1))
  expect_false(any(c("GeomHline", "GeomVline") %in% geoms2))
  expect_match(p2$labels$subtitle, "not identified")
})

test_that("the estimated range comes from the response detrended on the predictors", {
  skip_if_not_installed("gstat")
  pts <- sweep_points(n = 200, seed = 22)
  sw <- cv_block_size_sweep(pts, "z", "a", fit_fn = sweep_fit, k = 3, n_sizes = 2,
                            include_random = FALSE, quiet = TRUE)
  r <- attr(sw, "sac_range")
  ref <- suppressWarnings(estimate_sac_range(prep_model_data(pts, "z", "a"), "z",
                                             predictor_vars = "a", seed = 123L))
  expect_equal(r, as.numeric(ref))
})
