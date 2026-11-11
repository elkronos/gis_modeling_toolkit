# The checks that folds were built on the data (.check_fold_probe()) and that
# compare_models()'s fits share their rows and response (.same_rowset_xy(),
# .same_rowset()) read the same values on a layer with no CRS at any origin.
#
# Their margin was a part of the values' size: 1e-9 of the largest
# coordinate, 1000 units at an offset of 1e12, and 1e-10 of the largest
# response, 100 at a response near 1e12.  So folds built on a layer were
# accepted on the same layer moved by 500 units at an offset of 1e12 and
# refused at offset 0, and a response and the same plus 50 were read as one
# response.  The margin is now a part of the values' spread plus their
# rounding, which does not grow with the origin.
#
# Every shift below is exact: whole units, and an offset of 1e15 holds them.

.cm_layer <- function(n = 60L) {
  withr::with_seed(1L, {
    xy <- cbind(round(stats::runif(n, 0, 4000)), round(stats::runif(n, 0, 3000)))
    y <- round(stats::rnorm(n), 3)
  })
  list(xy = xy, y = y)
}
.cm_sf <- function(xy, y) {
  sf::st_as_sf(data.frame(y = y, X = xy[, 1], Y = xy[, 2]), coords = c("X", "Y"))
}
.cm_print <- function(d) {
  xy <- sf::st_coordinates(d)
  list(n = nrow(d), y = sort(d$y), x1 = sort(xy[, 1]), x2 = sort(xy[, 2]),
       lonlat = FALSE)
}

test_that("folds built on a layer are refused on it moved by 500 units, at any offset", {
  L <- .cm_layer()
  n <- nrow(L$xy)
  for (off in c(0, 1e12, 1e15)) {
    a <- .cm_sf(L$xy + off, L$y)
    stopifnot(all(sf::st_coordinates(a)[, 1] - off == L$xy[, 1]))
    f <- suppressWarnings(suppressMessages(
      make_folds(a, k = 5, method = "random_kfold", seed = 1)))
    a$..row_id <- seq_len(n)
    expect_silent(.check_fold_probe(f, a, "cv_probe"))
    b <- .cm_sf(L$xy + off + 500, L$y)
    b$..row_id <- seq_len(n)
    expect_error(.check_fold_probe(f, b, "cv_probe"), "built from different data",
                 info = sprintf("offset %g", off))
  }
})

test_that("compare_models()'s fingerprint tells rows and responses apart at any offset", {
  L <- .cm_layer()
  for (off in c(0, 1e12, 1e15)) {
    a <- .cm_print(.cm_sf(L$xy + off, L$y))
    b <- .cm_print(.cm_sf(L$xy + off + 500, L$y))
    lbl <- sprintf("offset %g", off)
    expect_true(.same_rowset_xy(a, a), info = lbl)
    expect_true(.same_rowset(a, a), info = lbl)
    expect_false(.same_rowset_xy(a, b), info = lbl)
    # The response at the offset, and the same plus 50.
    ya <- .cm_print(.cm_sf(L$xy, L$y + off))
    yb <- .cm_print(.cm_sf(L$xy, L$y + off + 50))
    expect_true(.same_rowset(ya, ya), info = lbl)
    expect_false(.same_rowset(ya, yb), info = lbl)
  }
})

test_that("the margin is a part of the spread plus the rounding of the size", {
  expect_identical(.same_values_tol(c(0, 0), rel = 1e-9), 0)
  expect_identical(.same_values_tol(numeric(0), rel = 1e-9), 0)
  expect_identical(.same_values_tol(c(NA, Inf), rel = 1e-9), 0)
  expect_equal(.same_values_tol(c(1e12, 1e12 + 4000), c(1e12 + 10, 1e12 + 3000),
                                rel = 1e-9),
               4e-6 + 8 * .Machine$double.eps * (1e12 + 4000))
  # A spread near the largest double does not overflow.
  expect_true(is.finite(.same_values_tol(c(-1e308, 1e308), rel = 1e-9)))
})
