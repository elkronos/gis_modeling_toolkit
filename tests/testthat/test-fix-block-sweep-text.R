# tests/testthat/test-fix-block-sweep-text.R
# ---------------------------------------------------------------------------
# The texts cv_block_size_sweep() writes about its fits.
#   * BUG-157: the max_fits refusal advised lowering `n_sizes`, which only
#     the default ladder reads; with `block_sizes` supplied the advice now
#     names fewer `block_sizes`, a smaller `k` and, when the count includes
#     it, include_random = FALSE.  On the default ladder it names
#     include_random = FALSE beside `n_sizes` and `k`, and when every rung
#     is too fine for make_folds() the refusal says so and advises larger
#     `block_sizes` instead of smaller ones.
#   * BUG-164: the print header and the plot caption gave the fits planned
#     before any fold was built, and "k folds per size", when a size built
#     fewer folds.  Both now count the fits from the table's `k` column,
#     the folds built, not the fits that succeeded.
# Each test here fails on the code before the fix.
# ---------------------------------------------------------------------------

bst_uniform <- function(n = 300, side = 1000, seed = 5) {
  set.seed(seed)
  d <- data.frame(x = runif(n, 0, side), y = runif(n, 0, side))
  d$a <- rnorm(n); d$z <- d$a + rnorm(n)
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
}

# Two clusters at opposite corners of a 5 km square: at block size 1600 the
# 3 x 3 grid has only two cells that hold points, so make_folds() builds 2
# folds there instead of 5.
bst_clusters <- function(seed = 5) {
  set.seed(seed)
  d <- data.frame(x = c(runif(100, 0, 1000), runif(100, 4000, 5000)),
                  y = c(runif(100, 0, 1000), runif(100, 4000, 5000)))
  d$a <- rnorm(200); d$z <- d$a + rnorm(200)
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32615)
}

bst_fit <- function(train_sf) lm_spatial_fit(train_sf, "z", "a")

bst_refusal <- function(...) {
  tryCatch({
    cv_block_size_sweep(..., quiet = TRUE)
    NA_character_
  }, error = function(e) conditionMessage(e))
}

test_that("the max_fits refusal names `block_sizes` for a supplied ladder (BUG-157)", {
  p <- bst_uniform()
  sizes <- seq(50, 230, by = 20)
  msg <- bst_refusal(p, "z", "a", fit_fn = bst_fit, block_sizes = sizes, k = 5,
                     max_fits = 20, sac = 100)
  expect_match(msg, "10 block sizes x 5 folds + 5 for the random reference = 55 model fits",
               fixed = TRUE)
  expect_match(msg, paste0("Raise `max_fits` if that cost is intended, or pass fewer ",
                           "`block_sizes`, a smaller `k` or include_random = FALSE."),
               fixed = TRUE)
  expect_false(grepl("n_sizes", msg, fixed = TRUE))
  # n_sizes is not read with a supplied ladder: the refusal is the same.
  expect_identical(bst_refusal(p, "z", "a", fit_fn = bst_fit, block_sizes = sizes,
                               n_sizes = 1, k = 5, max_fits = 20, sac = 100), msg)

  # Without the random reference the count does not include it, and the
  # advice does not offer to drop it.
  msg2 <- bst_refusal(p, "z", "a", fit_fn = bst_fit, block_sizes = sizes, k = 5,
                      include_random = FALSE, max_fits = 20, sac = 100)
  expect_match(msg2, "10 block sizes x 5 folds = 50 model fits", fixed = TRUE)
  expect_match(msg2, paste0("Raise `max_fits` if that cost is intended, or pass fewer ",
                            "`block_sizes` or a smaller `k`."), fixed = TRUE)
  expect_false(grepl("include_random", msg2, fixed = TRUE))
  expect_false(grepl("n_sizes", msg2, fixed = TRUE))

  # The default ladder keeps `n_sizes`, and lowering it is the remedy; the
  # count includes the random reference there too, which can be dropped.
  msg3 <- bst_refusal(p, "z", "a", fit_fn = bst_fit, n_sizes = 10, k = 5,
                      max_fits = 20, sac = 100)
  expect_match(msg3, "10 block sizes x 5 folds", fixed = TRUE)
  expect_match(msg3, paste0("Raise `max_fits` if that cost is intended, or lower ",
                            "`n_sizes` or `k`, or pass include_random = FALSE."),
               fixed = TRUE)
  expect_false(grepl("block_sizes`", msg3, fixed = TRUE))
  msg4 <- bst_refusal(p, "z", "a", fit_fn = bst_fit, n_sizes = 10, k = 5,
                      include_random = FALSE, max_fits = 20, sac = 100)
  expect_match(msg4, "10 block sizes x 5 folds = 50 model fits", fixed = TRUE)
  expect_match(msg4, "Raise `max_fits` if that cost is intended, or lower `n_sizes` or `k`.",
               fixed = TRUE)
  expect_false(grepl("include_random", msg4, fixed = TRUE))
  sw <- cv_block_size_sweep(p, "z", "a", fit_fn = bst_fit, n_sizes = 3, k = 5,
                            max_fits = 20, sac = 100, quiet = TRUE)
  expect_equal(nrow(sw), 4L)
})

test_that("the print header and plot caption count the fits made (BUG-164)", {
  q <- bst_clusters()
  expect_warning(
    sw <- cv_block_size_sweep(q, "z", "a", fit_fn = bst_fit, k = 5,
                              block_sizes = c(300, 1600), sac = NA, quiet = TRUE),
    "block size 1600: 2 folds", fixed = TRUE)
  expect_identical(sw$k, c(5L, 5L, 2L))
  # The budget the max_fits check read stays as it was.
  expect_identical(attr(sw, "n_fits"), 15L)
  out <- utils::capture.output(print(sw))
  expect_identical(out[1L], paste0("Cross-validation RMSE against block size ",
                                   "(12 fits, k = 5, sizes in EPSG:32615)"))
  skip_if_not_installed("ggplot2")
  cap <- plot(sw)$labels$caption
  expect_match(cap, "Band: fold-to-fold range; up to 5 folds per size, 12 fits in total; ",
               fixed = TRUE)
  expect_false(grepl("15 fits", cap, fixed = TRUE))
})

test_that("a sweep whose sizes all build k folds keeps its texts (BUG-164 control)", {
  set.seed(6)
  u <- sf::st_as_sf(data.frame(x = runif(200, 0, 5000), y = runif(200, 0, 5000),
                               a = rnorm(200)), coords = c("x", "y"), crs = 32615)
  u$z <- u$a + rnorm(200)
  sw <- cv_block_size_sweep(u, "z", "a", fit_fn = bst_fit, k = 5,
                            block_sizes = c(300, 1000), sac = NA, quiet = TRUE)
  expect_identical(sw$k, c(5L, 5L, 5L))
  out <- utils::capture.output(print(sw))
  expect_identical(out[1L], paste0("Cross-validation RMSE against block size ",
                                   "(15 fits, k = 5, sizes in EPSG:32615)"))
  # A row subset counts the fits of the rows it keeps.
  out2 <- utils::capture.output(print(sw[2:3, ]))
  expect_identical(out2[1L], paste0("Cross-validation RMSE against block size ",
                                    "(10 fits, k = 5, sizes in EPSG:32615)"))
  skip_if_not_installed("ggplot2")
  cap <- plot(sw)$labels$caption
  expect_match(cap, "Band: fold-to-fold range; 5 folds per size, 15 fits in total; ",
               fixed = TRUE)
  expect_false(grepl("up to", cap, fixed = TRUE))
})

test_that("the default-ladder refusal at the smallest n_sizes and k names include_random (BUG-157)", {
  p <- bst_uniform()
  # n_sizes = 1 and k = 2 cannot be lowered; the reference can be dropped.
  msg <- bst_refusal(p, "z", "a", fit_fn = bst_fit, n_sizes = 1, k = 2,
                     max_fits = 3, sac = 100)
  expect_match(msg, "1 block sizes x 2 folds + 2 for the random reference = 4 model fits",
               fixed = TRUE)
  expect_match(msg, "or pass include_random = FALSE.", fixed = TRUE)
  sw <- cv_block_size_sweep(p, "z", "a", fit_fn = bst_fit, n_sizes = 1, k = 2,
                            include_random = FALSE, max_fits = 3, sac = NA,
                            quiet = TRUE)
  expect_equal(nrow(sw), 1L)
})

test_that("a default ladder too fine for make_folds() advises larger block_sizes (BUG-157)", {
  # An E-W transect 6 km long and 1e-6 m wide: every rung of the default
  # ladder, built from the shorter side, needs a grid of trillions of blocks.
  set.seed(1)
  t <- data.frame(x = runif(200, 0, 6000), y = 5e6 + runif(200, 0, 1e-6))
  t$a <- rnorm(200); t$z <- t$a + rnorm(200)
  tp <- sf::st_as_sf(t, coords = c("x", "y"), crs = 32632)
  msg <- bst_refusal(tp, "z", "a", fit_fn = bst_fit, k = 5, sac = 100)
  expect_match(msg, "every size of the default ladder (", fixed = TRUE)
  expect_match(msg, "would need a grid of more than 1,000,000 blocks over the", fixed = TRUE)
  expect_match(msg, "which make_folds() refuses. Pass larger `block_sizes`; the largest ",
               fixed = TRUE)
  expect_false(grepl("smaller `block_sizes`", msg, fixed = TRUE))
  expect_false(grepl("fewer than|at least k", msg))
  # The size it names runs.
  b <- as.numeric(sub(".*holds k = 5 blocks is ([0-9.e+]+)\\.$", "\\1", msg))
  expect_true(is.finite(b) && b > 1000 && b < 6000 / 5)
  sw <- cv_block_size_sweep(tp, "z", "a", fit_fn = bst_fit, k = 5,
                            block_sizes = c(200, b), sac = 100, quiet = TRUE)
  expect_identical(sw$k, c(5L, 5L, 5L))
  # Smaller sizes, the old advice, are refused.
  expect_match(bst_refusal(tp, "z", "a", fit_fn = bst_fit, k = 5, block_sizes = 1e-8,
                           sac = 100),
               "would each need a grid of more than 1,000,000 blocks", fixed = TRUE)

  # A supplied ladder whose sizes all hold fewer than k blocks keeps its text.
  p <- bst_uniform()
  msg2 <- bst_refusal(p, "z", "a", fit_fn = bst_fit, k = 5, block_sizes = 900,
                      sac = 100)
  expect_match(msg2, "no block size leaves at least k = 5 blocks over the extent (",
               fixed = TRUE)
  expect_match(msg2, "; pass smaller `block_sizes` or a smaller `k`.", fixed = TRUE)
})

test_that("the header counts the folds built, not the fits that succeeded (BUG-164)", {
  p <- bst_uniform()
  p$flag <- seq_len(nrow(p)) == 1L
  # The fit fails on the one fold that holds out row 1, at every size.
  fail_fit <- function(train_sf) {
    if (!any(train_sf$flag)) stop("row 1 held out")
    lm_spatial_fit(train_sf, "z", "a")
  }
  sw <- suppressWarnings(
    cv_block_size_sweep(p, "z", "a", fit_fn = fail_fit, k = 5,
                        block_sizes = c(150, 300), sac = NA, quiet = TRUE))
  expect_identical(sw$k, c(5L, 5L, 5L))
  expect_true(all(sw$n_folds_succeeded == 4L))
  out <- utils::capture.output(print(sw))
  expect_identical(out[1L], paste0("Cross-validation RMSE against block size ",
                                   "(15 fits, k = 5, sizes in EPSG:32632)"))
  skip_if_not_installed("ggplot2")
  cap <- plot(sw)$labels$caption
  expect_match(cap, "Band: fold-to-fold range; 5 folds per size, 15 fits in total; ",
               fixed = TRUE)
})

test_that("the help says what n_fits is and what print() and plot() count", {
  # The page from man/ in the source tree, or the installed one under
  # R CMD check, which runs the tests where there is no man/.
  file <- test_path("..", "..", "man", "cv_block_size_sweep.Rd")
  rd <- if (file.exists(file)) tools::parse_Rd(file) else
    tools::Rd_db("spatialkit")[["cv_block_size_sweep.Rd"]]
  expect_s3_class(rd, "Rd")
  txt <- gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  expect_match(txt, paste0("n_fits (the fits planned, which the max_fits check reads; ",
                           "print() and plot() count the fits made from the k column)"),
               fixed = TRUE)
  expect_match(txt, paste0("lower n_sizes or k to cut it (with block_sizes supplied: pass ",
                           "fewer block_sizes or a smaller k), or pass include_random = ",
                           "FALSE when the count includes the reference."),
               fixed = TRUE)
  expect_match(txt, paste0("when that drops every size, the sweep refuses and names the ",
                           "largest size whose grid holds k blocks, to pass as block_sizes."),
               fixed = TRUE)
})
