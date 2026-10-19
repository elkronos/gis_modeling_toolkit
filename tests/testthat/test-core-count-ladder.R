# tests/testthat/test-core-count-ladder.R
# ---------------------------------------------------------------------------
# resolution_profile(obs_per_cell = ): the ladder a size sets, its bounds and
# its warnings; how "auto" gives way to a ladder the caller describes;
# nstart = NULL and sac = "auto"; Moran's z on such a ladder; and the two
# bounds on its cost.
#
# Every call names `partition`, `nstart`, `estimate`, `obs_per_cell` and
# `sac`.  Unless a test says otherwise the regions are of the density kind
# from 3 k-means++ starts.  A layer is drawn under the seed its test states:
# the locations (x, then y), then the response.  A backend is switched off by
# replacing .has_backend(), the one function resolution_profile() asks.
# ---------------------------------------------------------------------------

# ---- Layers and helpers ------------------------------------------------------

.cl_seed <- function(seed, env = parent.frame())
  withr::local_seed(seed, .local_envir = env, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

# `n` locations uniform on a square of side `side`, with the field of the
# engine's tests: sin(x / 2e4) + cos(y / 1.5e4) plus normal noise of standard
# deviation 1.5.
.cl_layer <- function(n, seed, side = 1e5) {
  d <- local({
    .cl_seed(seed)
    x <- stats::runif(n, 0, side); y <- stats::runif(n, 0, side)
    data.frame(x = x, y = y, z = sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(n, 0, 1.5))
  })
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32617)
}

# A profile of the density kind from 3 k-means++ starts, every setting that
# has a default the switch moves written out.  What a test varies it names.
.cl_profile <- function(pts, ..., response = "z", partition = "density", nstart = 3,
                        sac = "auto", estimate = "pooled", obs_per_cell = 4.5)
  resolution_profile(pts, response, ..., partition = partition, nstart = nstart, sac = sac,
                     estimate = estimate, obs_per_cell = obs_per_cell)

# A call with its warnings kept and muffled.
.cl_catch <- function(expr) {
  got <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    got[[length(got) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = got)
}
.cl_classes <- function(r) vapply(r$warnings, function(w) class(w)[1L], character(1))
.cl_of <- function(r, class) Filter(function(w) inherits(w, class), r$warnings)
.cl_fields <- function(w) unclass(w)[setdiff(names(unclass(w)), c("message", "call"))]

# The error a call stops with.
.cl_error <- function(expr) tryCatch({ suppressWarnings(expr); NULL }, error = function(e) e)

# The columns of a profile, without its attributes.
.cl_cols <- function(p) stats::setNames(lapply(seq_along(p), function(j) p[[j]]), names(p))

.cl_without <- function(pkg, env = parent.frame())
  testthat::local_mocked_bindings(
    .has_backend = function(p)
      if (p %in% pkg) FALSE else requireNamespace(p, quietly = TRUE),
    .env = env)

.cl_store <- new.env(parent = emptyenv())
.cl_keep <- function(key, expr) {
  if (is.null(.cl_store[[key]])) .cl_store[[key]] <- expr
  .cl_store[[key]]
}
.cl_3000 <- function() .cl_keep("l3000", .cl_layer(3000, 301))
.cl_600  <- function() .cl_keep("l600", .cl_layer(600, 302))


# ---- The levels of a size ----------------------------------------------------

test_that("obs_per_cell gives the levels rung_k() gives", {
  pts   <- .cl_3000()
  sizes <- c(18, 12, 9, 6, 4.5, 3, 2)
  want  <- vapply(sizes, function(s) rung_k(3000L, s), integer(1))
  expect_identical(want, c(167L, 250L, 333L, 500L, 667L, 1000L, 1500L))

  r <- .cl_catch(.cl_profile(pts, estimate = "mean", obs_per_cell = sizes))
  expect_length(r$warnings, 0L)
  expect_identical(r$value$levels, want)
  expect_identical(attr(r$value, "bounds")$obs_per_cell, sizes)
  expect_identical(attr(r$value, "bounds")$size_levels,
                   stats::setNames(want, c("18", "12", "9", "6", "4.5", "3", "2")))
  # Every partition was built, on every row.
  expect_true(all(is.finite(r$value$wss)))
  expect_identical(vapply(attr(r$value, "centres"), nrow, integer(1), USE.NAMES = FALSE), want)

  # The order the sizes are given in, and a size given twice, change nothing.
  expect_identical(.cl_profile(pts, estimate = "mean", obs_per_cell = c(4.5, 18, 4.5)),
                   .cl_profile(pts, estimate = "mean", obs_per_cell = c(18, 4.5)))

  # Two sizes that give one count are one level.
  two <- .cl_profile(pts, obs_per_cell = c(4.5, 4.501))
  expect_identical(two$levels, 667L)
  expect_identical(unname(attr(two, "bounds")$size_levels), c(667L, 667L))
})

test_that("a size divides the rows with a response, and every point without one", {
  pts <- .cl_3000()
  # A response on 2,000 of the rows.
  part <- pts
  part$z[seq(3, 3000, by = 3)] <- NA
  expect_identical(sum(is.finite(part$z)), 2000L)
  p <- suppressWarnings(.cl_profile(part, obs_per_cell = 9))
  expect_identical(p$levels, rung_k(2000L, 9))
  expect_identical(p$levels, 222L)
  expect_identical(attr(p, "estimate")$n, 2000L)
  expect_identical(attr(p, "bounds")$n_response, 2000L)
  expect_identical(attr(p, "bounds")$n_scored, 2000L)
  # The centres are fitted to every point, with a response or without.
  expect_identical(attr(p, "bounds")$n_sample, 3000L)
  expect_identical(attr(p, "bounds")$n, 3000L)

  # No response: the points.
  geo <- .cl_profile(pts, response = NULL, estimate = NULL, obs_per_cell = 9)
  expect_identical(geo$levels, rung_k(3000L, 9))
  expect_identical(attr(geo, "bounds")$n_response, NA_integer_)
  expect_null(attr(geo, "estimate"))
  expect_false("loo" %in% names(geo))

  # A response with no finite value: the points, as the support ceiling counts.
  none <- pts
  none$z <- NA_real_
  empty <- suppressWarnings(.cl_profile(none, estimate = NULL, obs_per_cell = 9))
  expect_identical(empty$levels, rung_k(3000L, 9))
  expect_identical(attr(empty, "bounds")$n_response, 0L)
})

test_that("sample_n changes no number under a numeric obs_per_cell", {
  pts  <- .cl_layer(1200, 303)
  base <- .cl_catch(.cl_profile(pts, obs_per_cell = c(9, 4.5)))
  expect_length(base$warnings, 0L)
  expect_identical(attr(base$value, "bounds")$n_sample, 1200L)
  for (sn in c(3, 100, 400, 1199, 1e6)) {
    r <- .cl_catch(.cl_profile(pts, obs_per_cell = c(9, 4.5), sample_n = sn))
    expect_identical(r$value, base$value, info = paste("sample_n", sn))
    expect_identical(.cl_classes(r), "spatialkit_ladder_unused", info = paste("sample_n", sn))
  }
  # Nor without an estimate: no subsample is drawn on a ladder of sizes.
  plain <- .cl_profile(pts, estimate = NULL, obs_per_cell = c(9, 4.5))
  expect_identical(attr(plain, "bounds")$n_sample, 1200L)
  for (sn in c(100, 400)) {
    r <- .cl_catch(.cl_profile(pts, estimate = NULL, obs_per_cell = c(9, 4.5), sample_n = sn))
    expect_identical(r$value, plain, info = paste("no estimate, sample_n", sn))
    expect_identical(attr(r$value, "bounds")$n_sample, 1200L)
    expect_identical(attr(r$value, "bounds")$n_scored, 1200L)
  }
  # The same number draws a subsample on the ladder the function builds.
  sub <- suppressWarnings(.cl_profile(pts, obs_per_cell = NULL, sample_n = 400, n_levels = 4))
  expect_identical(attr(sub, "bounds")$n_sample, 400L)
})


# ---- How "auto" is resolved --------------------------------------------------

test_that("obs_per_cell = \"auto\" is one size for pooled values and seven for plain means", {
  pts <- .cl_600()
  pooled <- .cl_catch(.cl_profile(pts, estimate = "pooled", obs_per_cell = "auto"))
  expect_length(pooled$warnings, 0L)
  expect_identical(pooled$value$levels, 133L)
  expect_identical(attr(pooled$value, "bounds")$obs_per_cell, 4.5)
  expect_identical(pooled$value, .cl_profile(pts, estimate = "pooled", obs_per_cell = 4.5))

  plain <- .cl_catch(.cl_profile(pts, estimate = "mean", obs_per_cell = "auto"))
  expect_length(plain$warnings, 0L)
  expect_identical(nrow(plain$value), 7L)
  expect_identical(plain$value$levels, c(33L, 50L, 67L, 100L, 133L, 200L, 300L))
  expect_identical(attr(plain$value, "bounds")$obs_per_cell, c(18, 12, 9, 6, 4.5, 3, 2))
  expect_identical(plain$value,
                   .cl_profile(pts, estimate = "mean", obs_per_cell = c(18, 12, 9, 6, 4.5, 3, 2)))

  # estimate = "auto" is the pooled value on this layer, so its ladder is one size.
  expect_identical(.cl_profile(pts, estimate = "auto", obs_per_cell = "auto"), pooled$value)
})

test_that("obs_per_cell = \"auto\" gives way to `levels` and to the automatic ladder's arguments", {
  skip_if_not_installed("gstat")
  pts <- .cl_600()
  for (est in c("pooled", "mean")) {
    # `levels`: scored as given, on the present path.
    lv <- c(10L, 40L)
    a <- suppressWarnings(.cl_profile(pts, estimate = est, obs_per_cell = "auto", levels = lv))
    b <- suppressWarnings(.cl_profile(pts, estimate = est, obs_per_cell = NULL, levels = lv))
    expect_identical(a$levels, lv, info = est)
    expect_identical(attr(a, "bounds"), attr(b, "bounds"), info = est)
    expect_identical(a, b, info = est)
    expect_null(attr(a, "bounds")$obs_per_cell)
    expect_true(all(is.finite(a$loo)), info = est)

    # Each of the four arguments that describe the automatic ladder.
    for (arg in list(list(n_levels = 5L), list(min_cell_n = 12L), list(sample_n = 300L),
                     list(range_floor = TRUE), list(range_floor = FALSE))) {
      what <- paste(est, names(arg), format(arg[[1L]]))
      a <- suppressWarnings(do.call(.cl_profile,
                                    c(list(pts, estimate = est, obs_per_cell = "auto"), arg)))
      b <- suppressWarnings(do.call(.cl_profile,
                                    c(list(pts, estimate = est, obs_per_cell = NULL), arg)))
      expect_identical(a, b, info = what)
      expect_null(attr(a, "bounds")$obs_per_cell, info = what)
      expect_false(attr(a, "bounds")$levels_given, info = what)
      expect_gt(nrow(a), 2L)
      expect_true("loo" %in% names(a), info = what)
    }
  }
})

test_that("a wrapper that passes n_levels on counts as supplying it", {
  skip_if_not_installed("gstat")
  pts <- .cl_600()
  wrapper <- function(layer, n_levels = 6L)
    resolution_profile(layer, "z", n_levels = n_levels, partition = "density", nstart = 3,
                       sac = "auto", estimate = "pooled", obs_per_cell = "auto")
  w <- suppressWarnings(wrapper(pts))
  expect_null(attr(w, "bounds")$obs_per_cell)
  expect_gt(nrow(w), 2L)
  expect_identical(w, suppressWarnings(.cl_profile(pts, obs_per_cell = NULL, n_levels = 6L)))
  # Beside a size it is the argument that is not read.
  sized <- function(layer, n_levels = 6L)
    resolution_profile(layer, "z", n_levels = n_levels, partition = "density", nstart = 3,
                       sac = "auto", estimate = "pooled", obs_per_cell = 4.5)
  r <- .cl_catch(sized(pts))
  expect_identical(.cl_classes(r), "spatialkit_ladder_unused")
  expect_identical(r$value$levels, 133L)
})

test_that("a numeric obs_per_cell refuses `levels` and names the ladder arguments it does not read", {
  pts  <- .cl_600()
  base <- .cl_profile(pts, obs_per_cell = 4.5)
  e <- .cl_error(.cl_profile(pts, obs_per_cell = 4.5, levels = 20L))
  expect_identical(conditionMessage(e),
                   "resolution_profile(): pass one of `obs_per_cell` and `levels`: each sets the ladder.")
  expect_null(conditionCall(e))

  args <- list(n_levels = 5L, min_cell_n = 12L, sample_n = 300L, range_floor = TRUE)
  for (nm in names(args)) {
    r <- .cl_catch(do.call(.cl_profile, c(list(pts, obs_per_cell = 4.5), args[nm])))
    expect_identical(.cl_classes(r), "spatialkit_ladder_unused", info = nm)
    w <- r$warnings[[1L]]
    expect_identical(conditionMessage(w),
                     sprintf(paste0("resolution_profile(): `%s` is not read when `obs_per_cell` ",
                                    "sets the ladder: the levels are the sizes asked for, and the ",
                                    "centres are fitted to every row. Leave `obs_per_cell` NULL for ",
                                    "the ladder it describes."), nm))
    expect_identical(.cl_fields(w), list(unused = nm, obs_per_cell = 4.5), info = nm)
    expect_null(conditionCall(w))
    # The ladder is the size's, whatever the argument holds.
    expect_identical(.cl_cols(r$value), .cl_cols(base), info = nm)
    expect_identical(attr(r$value, "bounds")$size_levels, c("4.5" = 133L), info = nm)
  }
  # The four together are one warning.
  r <- .cl_catch(do.call(.cl_profile, c(list(pts, obs_per_cell = c(9, 4.5)), args)))
  expect_identical(.cl_classes(r), "spatialkit_ladder_unused")
  expect_identical(conditionMessage(r$warnings[[1L]]),
                   paste0("resolution_profile(): `n_levels`, `min_cell_n`, `sample_n` and ",
                          "`range_floor` are not read when `obs_per_cell` sets the ladder: the ",
                          "levels are the sizes asked for, and the centres are fitted to every ",
                          "row. Leave `obs_per_cell` NULL for the ladder they describe."))
  expect_identical(.cl_fields(r$warnings[[1L]]),
                   list(unused = names(args), obs_per_cell = c(9, 4.5)))
  # The condition carries the sizes of the ladder, largest first and each
  # once, whatever order they were given in.
  r <- .cl_catch(.cl_profile(pts, obs_per_cell = c(4.5, 9, 4.5), n_levels = 5L))
  expect_identical(.cl_fields(r$warnings[[1L]]),
                   list(unused = "n_levels", obs_per_cell = c(9, 4.5)))
  r <- .cl_catch(do.call(.cl_profile, c(list(pts, obs_per_cell = c(9, 4.5)), args)))
  # `min_cell_n` is kept in the bounds, and is not read.
  expect_identical(attr(r$value, "bounds")$min_cell_n, 12L)
  expect_identical(attr(r$value, "bounds")$range_floor, TRUE)
  # range_floor = FALSE asks for what a ladder of sizes does anyway.
  r <- .cl_catch(.cl_profile(pts, obs_per_cell = 4.5, range_floor = FALSE))
  expect_length(r$warnings, 0L)
  expect_identical(r$value, base)
})

test_that("without an estimate obs_per_cell = \"auto\" is NULL", {
  skip_if_not_installed("gstat")
  pts <- .cl_600()
  a <- suppressWarnings(.cl_profile(pts, estimate = NULL, obs_per_cell = "auto"))
  b <- suppressWarnings(.cl_profile(pts, estimate = NULL, obs_per_cell = NULL))
  expect_identical(a, b)
  expect_null(attr(a, "bounds")$obs_per_cell)
  expect_null(attr(a, "estimate"))
  expect_identical(names(attr(a, "bounds")),
                   c("floor", "ceiling", "ceiling_from", "supported", "area", "range", "n",
                     "n_sample", "n_distinct", "min_cell_n", "range_floor", "n_response",
                     "n_scored", "levels_given", "predictors"))
  # And without a response.
  expect_identical(.cl_profile(pts, response = NULL, estimate = NULL, obs_per_cell = "auto"),
                   .cl_profile(pts, response = NULL, estimate = NULL, obs_per_cell = NULL))
  expect_identical(.cl_profile(pts, response = NULL, estimate = "auto", obs_per_cell = "auto"),
                   .cl_profile(pts, response = NULL, estimate = NULL, obs_per_cell = NULL))
})


# ---- The bounds of a ladder of sizes -----------------------------------------

test_that("the bounds of a ladder of sizes follow the finest size", {
  skip_if_not_installed("gstat")
  pts <- .cl_3000()
  r <- .cl_catch(.cl_profile(pts, obs_per_cell = c(18, 4.5), min_cell_n = 7))
  expect_identical(.cl_classes(r), "spatialkit_ladder_unused")
  b <- attr(r$value, "bounds")

  # The formulas, written again.
  n <- 3000; s_min <- 4.5
  cap  <- min(b$n_distinct, b$n - 1L)
  want <- round(n / s_min)
  too_few <- min(want, cap) < 2
  expect_identical(cap, 2999L)
  expect_identical(want, 667)
  expect_false(too_few)
  expect_identical(b$ceiling, as.integer(max(2, min(want, cap))))
  expect_identical(b$ceiling, 667L)
  expect_identical(b$ceiling_from, "obs_per_cell")
  expect_identical(b$supported, !too_few && 2 <= max(2, min(want, cap)))
  expect_true(b$supported)
  expect_identical(b$obs_per_cell, c(18, 4.5))
  expect_identical(b$size_levels, c("18" = 167L, "4.5" = 667L))
  expect_identical(r$value$levels, c(167L, 667L))

  # The other fields.
  expect_identical(b$n_sample, 3000L)
  expect_identical(b$n_scored, 3000L)
  expect_true(b$levels_given)
  expect_identical(b$min_cell_n, 7L)
  expect_identical(b$floor, 2L)
  expect_identical(b$range, NA_real_)
  expect_identical(b$range_floor, FALSE)
  expect_identical(names(b),
                   c("floor", "ceiling", "ceiling_from", "supported", "area", "range", "n",
                     "n_sample", "n_distinct", "min_cell_n", "range_floor", "n_response",
                     "n_scored", "levels_given", "predictors", "obs_per_cell", "size_levels"))
  present <- attr(suppressWarnings(.cl_profile(pts, estimate = NULL, obs_per_cell = NULL,
                                               sac = NULL, n_levels = 3)), "bounds")
  for (f in c("area", "n", "n_distinct", "n_response", "predictors"))
    expect_identical(b[[f]], present[[f]], info = f)
  # The present ladder of the same layer fits a subsample and has no size.
  expect_identical(present$n_sample, 1500L)
  expect_false(present$levels_given)
})

test_that("rows without a response lower a size's count and raise no warning about a ceiling", {
  skip_if_not_installed("gstat")
  pts <- .cl_layer(900, 304)
  pts$z[1:450] <- NA
  # The present ladder says that the rows without a response moved its ceiling.
  present <- .cl_catch(.cl_profile(pts, estimate = NULL, obs_per_cell = NULL, sac = NULL,
                                   n_levels = 3))
  expect_true("spatialkit_response_ceiling" %in% .cl_classes(present))
  # A ladder of sizes counts those rows already.
  for (sac in list("auto", NULL)) {
    r <- .cl_catch(.cl_profile(pts, obs_per_cell = c(9, 4.5), sac = sac))
    expect_false("spatialkit_response_ceiling" %in% .cl_classes(r))
    expect_false(any(grepl("spatialkit_", .cl_classes(r))))
    expect_identical(r$value$levels, c(50L, 100L))
    expect_identical(attr(r$value, "bounds")$n_response, 450L)
    expect_identical(attr(r$value, "bounds")$n_sample, 900L)
    expect_identical(attr(r$value, "bounds")$ceiling_from, "obs_per_cell")
  }
})

test_that("fewer than two cells of the finest size is not supported, and says so once", {
  pts <- .cl_layer(6, 305)
  r <- .cl_catch(.cl_profile(pts, obs_per_cell = 4.5))
  expect_identical(.cl_classes(r), "spatialkit_not_supported")
  w <- r$warnings[[1L]]
  expect_identical(conditionMessage(w),
                   paste0("resolution_profile(): 6 points at obs_per_cell = 4.5 make 1 cell, and it ",
                          "takes two, so there is no cell count of that size. The profile scores 2 ",
                          "cells (3 rows per cell) so that level is still visible, and `supported` ",
                          "is FALSE. An `obs_per_cell` of at most 4 gives two cells here."))
  expect_identical(.cl_fields(w),
                   list(reason = "fewer than two cells of the size", n = 6L,
                        min_cell_n = NA_integer_, obs_per_cell = 4.5))
  b <- attr(r$value, "bounds")
  expect_false(b$supported)
  expect_identical(r$value$levels, 2L)
  expect_identical(b$ceiling, 2L)
  expect_identical(b$ceiling_from, "fewer than two cells of the size")
  expect_identical(b$size_levels, c("4.5" = 2L))
  expect_true(is.finite(r$value$wss))
  # The size the message offers does give two cells, and is supported.
  ok <- .cl_catch(.cl_profile(pts, obs_per_cell = 4))
  expect_length(ok$warnings, 0L)
  expect_true(attr(ok$value, "bounds")$supported)
  expect_identical(ok$value$levels, 2L)
  # Every size of several too coarse: still one warning, about the finest.
  r <- .cl_catch(.cl_profile(pts, obs_per_cell = c(18, 4.5)))
  expect_identical(.cl_classes(r), "spatialkit_not_supported")
  expect_identical(r$warnings[[1L]]$obs_per_cell, 4.5)
  expect_identical(attr(r$value, "bounds")$size_levels, c("18" = 2L, "4.5" = 2L))
  # The seven sizes of "auto" on rows too few for any of them: none is
  # dropped.  Two rows with a response among five points.
  two <- .cl_layer(5, 305)
  two$z[3:5] <- NA
  r <- .cl_catch(.cl_profile(two, estimate = "mean", obs_per_cell = "auto"))
  expect_identical(.cl_classes(r), "spatialkit_not_supported")
  expect_identical(conditionMessage(r$warnings[[1L]]),
                   paste0("resolution_profile(): 2 rows with a response at obs_per_cell = 2 make 1 ",
                          "cell, and it takes two, so there is no cell count of that size. The ",
                          "profile scores 2 cells (1 row per cell) so that level is still visible, ",
                          "and `supported` is FALSE. An `obs_per_cell` of at most 1.33 gives two ",
                          "cells here."))
  expect_identical(unname(attr(r$value, "bounds")$size_levels), rep(2L, 7L))
  expect_identical(r$value$levels, 2L)
  # One row with a response: no size gives two cells, and the words are for one row.
  one <- .cl_layer(5, 305)
  one$z[2:5] <- NA
  r <- .cl_catch(.cl_profile(one, obs_per_cell = 4.5))
  expect_identical(.cl_classes(r), "spatialkit_not_supported")
  expect_identical(conditionMessage(r$warnings[[1L]]),
                   paste0("resolution_profile(): 1 row with a response at obs_per_cell = 4.5 makes 0 ",
                          "cells, and it takes two, so there is no cell count of that size. The ",
                          "profile scores 2 cells (0.5 rows per cell) so that level is still ",
                          "visible, and `supported` is FALSE. Two cells take a response on at least ",
                          "two rows, whatever `obs_per_cell` is."))
  # One row is too few rows, not a constant response: the profile says so.
  e <- attr(r$value, "estimate")
  expect_identical(e$n, 1L)
  expect_false(e$constant)
  expect_identical(e$loo_one_cell, NA_real_)
  expect_identical(r$value$loo, NA_real_)
  out <- utils::capture.output(print(r$value))
  expect_identical(out[2L],
                   paste0("  ladder      : 2 cells (obs_per_cell = 4.5 on the 1 row with a response; ",
                          "centres fitted to every row; floor 2 (no range))  -- fewer than two cells ",
                          "of the size: not supported"))
  expect_true(paste0("  estimate    : pooled toward the 8 nearest cells; loo on 1 row of the layer; ",
                     "one cell: loo NA (fewer than two rows)") %in% out)
})

test_that("the words for rows per cell are singular for one", {
  expect_identical(.per_cell_text(3), "3 rows per cell")
  expect_identical(.per_cell_text(1), "1 row per cell")
  expect_identical(.per_cell_text(1.0004), "1 row per cell")
  expect_identical(.per_cell_text(1.03), "1.03 rows per cell")
  expect_identical(.per_cell_text(0.5), "0.5 rows per cell")
  expect_identical(.per_cell_text(1234.5), "1.23e+03 rows per cell")
})

test_that("the finest level the guards read is the ceiling the ladder reports", {
  # The guards run before the ladder is built, on the count the finest size
  # gives; .size_ladder() must build that count as its finest level.
  for (n in c(1L, 2L, 3L, 6L, 20L, 26L, 27L, 100L, 3000L))
    for (cap in c(1L, 2L, 5L, 99L, 2999L))
      for (sizes in list(4.5, c(18, 4.5), c(18, 12, 9, 6, 4.5, 3, 2), 1, c(400, 50), 2.5))
        for (auto in c(TRUE, FALSE)) {
          what <- paste(n, cap, paste(sizes, collapse = "/"), auto)
          sl <- suppressWarnings(.size_ladder(sizes, auto, n, max(n, 3L), cap, cap, 2, NA_real_,
                                              "points", "a response"))
          expect_identical(sl$ceiling, .finest_level(n, sizes, cap), info = what)
          expect_identical(max(sl$levels), sl$ceiling, info = what)
          expect_gte(sl$ceiling, 2L)
        }
  expect_identical(.finest_level(3000L, c(18, 4.5), 2999L), 667L)
  expect_identical(.finest_level(3000L, 4.5, 100L), 100L)
  expect_identical(.finest_level(6L, 4.5, 5L), 2L)
  expect_identical(.finest_level(40L, 1, 39L), 39L)
})

test_that("a coarser size of fewer than two cells is held to two when asked for, and dropped from \"auto\"", {
  pts <- .cl_layer(20, 306)
  r <- .cl_catch(.cl_profile(pts, obs_per_cell = c(18, 4.5)))
  expect_identical(.cl_classes(r), "spatialkit_size_not_reached")
  w <- r$warnings[[1L]]
  expect_identical(conditionMessage(w),
                   paste0("resolution_profile(): obs_per_cell = 18 on 20 rows is fewer than two ",
                          "cells; two cells are used (10 rows per cell)."))
  expect_identical(.cl_fields(w), list(obs_per_cell = 18, wanted = 1L, used = 2L,
                                       reason = "fewer than two cells"))
  expect_identical(r$value$levels, c(2L, 4L))
  expect_identical(attr(r$value, "bounds")$size_levels, c("18" = 2L, "4.5" = 4L))
  expect_true(attr(r$value, "bounds")$supported)

  # From "auto" the size is left out, and nothing is said.
  r <- .cl_catch(.cl_profile(pts, estimate = "mean", obs_per_cell = "auto"))
  expect_length(r$warnings, 0L)
  expect_identical(r$value$levels, c(2L, 3L, 4L, 7L, 10L))
  b <- attr(r$value, "bounds")
  expect_identical(b$obs_per_cell, c(18, 12, 9, 6, 4.5, 3, 2))
  expect_identical(b$size_levels,
                   c("18" = NA, "12" = 2L, "9" = 2L, "6" = 3L, "4.5" = 4L, "3" = 7L, "2" = 10L))
  expect_true(b$supported)
  expect_identical(b$ceiling, 10L)
  # The same sizes given as numbers keep the coarsest, at two cells, with the warning.
  r <- .cl_catch(.cl_profile(pts, estimate = "mean", obs_per_cell = c(18, 12, 9, 6, 4.5, 3, 2)))
  expect_identical(.cl_classes(r), "spatialkit_size_not_reached")
  expect_identical(attr(r$value, "bounds")$size_levels[["18"]], 2L)
  expect_identical(r$value$levels, c(2L, 3L, 4L, 7L, 10L))
})

test_that("a size that asks for more cells than can be placed is held to that many", {
  # 3,000 rows at 100 stations.
  pts <- .cl_3000()
  at  <- local({ .cl_seed(307); sample.int(100L, 3000L, replace = TRUE) })
  expect_identical(length(unique(at)), 100L)
  sf::st_geometry(pts) <- sf::st_geometry(pts)[at]
  r <- .cl_catch(.cl_profile(pts, obs_per_cell = 4.5))
  expect_identical(.cl_classes(r),
                   c("spatialkit_size_not_reached", "spatialkit_repeated_locations"))
  w <- r$warnings[[1L]]
  expect_identical(conditionMessage(w),
                   paste0("resolution_profile(): obs_per_cell = 4.5 asks for 667 cells, more than ",
                          "k-means can place on the layer's 100 distinct locations; 100 cells are ",
                          "used (30 rows per cell)."))
  expect_identical(.cl_fields(w), list(obs_per_cell = 4.5, wanted = 667L, used = 100L,
                                       reason = "more cells than can be placed"))
  expect_identical(r$value$levels, 100L)
  b <- attr(r$value, "bounds")
  expect_identical(b$ceiling, 100L)
  expect_identical(b$ceiling_from, "distinct locations")
  expect_true(b$supported)
  expect_identical(b$size_levels, c("4.5" = 100L))
  expect_identical(r$value$cells_resp, 100L)
  # The sizes of "auto" are held the same way, each with its warning.
  r <- .cl_catch(.cl_profile(pts, estimate = "mean", obs_per_cell = "auto"))
  held <- .cl_of(r, "spatialkit_size_not_reached")
  expect_identical(vapply(held, function(w) w$obs_per_cell, numeric(1)), c(18, 12, 9, 6, 4.5, 3, 2))
  expect_true(all(vapply(held, function(w) w$reason, "") == "more cells than can be placed"))
  expect_identical(r$value$levels, 100L)

  # A size that asks for exactly as many cells as can be placed is reached,
  # and nothing is said of it: 50 stations visited six times, at six rows a
  # cell.  One row a cell more and it is held.
  st <- .cl_layer(300, 325)
  sf::st_geometry(st) <- sf::st_geometry(st)[rep(1:50, 6)]
  r <- .cl_catch(.cl_profile(st, obs_per_cell = 6))
  expect_length(.cl_of(r, "spatialkit_size_not_reached"), 0L)
  expect_identical(r$value$levels, 50L)
  expect_identical(attr(r$value, "bounds")$ceiling_from, "obs_per_cell")
  expect_identical(attr(r$value, "bounds")$size_levels, c("6" = 50L))
  r <- .cl_catch(.cl_profile(st, obs_per_cell = 5))
  expect_length(.cl_of(r, "spatialkit_size_not_reached"), 1L)
  expect_identical(r$value$levels, 50L)
  expect_identical(attr(r$value, "bounds")$ceiling_from, "distinct locations")

  # No location repeated: one short of the points, which stats::kmeans() needs.
  few <- .cl_layer(40, 308)
  r <- .cl_catch(.cl_profile(few, estimate = NULL, response = NULL, obs_per_cell = 1,
                             partition = "kmeans", nstart = 2))
  expect_identical(.cl_classes(r), "spatialkit_size_not_reached")
  expect_identical(conditionMessage(r$warnings[[1L]]),
                   paste0("resolution_profile(): obs_per_cell = 1 asks for 40 cells, more than ",
                          "k-means can place on 40 points; 39 cells are used (1.03 rows per cell)."))
  expect_identical(r$value$levels, 39L)
})

test_that("a range floor above the finest level of a ladder of sizes is not supported", {
  pts <- .cl_3000()
  # A range of 2,000 on a 1e5 square: about 2,500 cells no wider than it,
  # against the 667 of the size.
  r <- .cl_catch(.cl_profile(pts, obs_per_cell = 4.5, sac = 2000))
  w <- .cl_of(r, "spatialkit_not_supported")
  expect_length(w, 1L)
  expect_false("spatialkit_size_not_reached" %in% .cl_classes(r))
  b <- attr(r$value, "bounds")
  expect_false(b$supported)
  expect_identical(b$ceiling_from, "obs_per_cell")
  expect_identical(b$floor, as.integer(ceiling(b$area / 2000^2)))
  expect_gt(b$floor, 667L)
  expect_identical(b$range, 2000)
  expect_identical(conditionMessage(w[[1L]]),
                   sprintf(paste0("resolution_profile(): cells no wider than the autocorrelation ",
                                  "range (2000) would need at least %d of them, but 3000 points at ",
                                  "obs_per_cell = 4.5 support at most 667. The data cannot support a ",
                                  "tessellation that respects their own correlation structure; the ",
                                  "profile scores the 667 cells that size gives so the cost of each ",
                                  "level is still visible."), b$floor))
  expect_false(grepl("min_cell_n", conditionMessage(w[[1L]]), fixed = TRUE))
  expect_identical(.cl_fields(w[[1L]]),
                   list(reason = "floor above ceiling", n = 3000L, min_cell_n = NA_integer_,
                        obs_per_cell = 4.5))
  # The floor is reported and not applied: the level is the size's.
  expect_identical(r$value$levels, 667L)
  # A range the size's cells are narrower than is supported.
  ok <- .cl_catch(.cl_profile(pts, obs_per_cell = 4.5, sac = 20000))
  expect_length(.cl_of(ok, "spatialkit_not_supported"), 0L)
  expect_true(attr(ok$value, "bounds")$supported)
  expect_identical(attr(ok$value, "bounds")$floor, as.integer(ceiling(b$area / 20000^2)))
})

test_that("points on one location are not supported on a ladder of sizes", {
  pts <- .cl_layer(12, 309)
  sf::st_geometry(pts) <- sf::st_geometry(pts)[rep(1L, 12L)]
  r <- .cl_catch(.cl_profile(pts, estimate = NULL, obs_per_cell = 4.5, partition = "kmeans",
                             nstart = 2))
  w <- .cl_of(r, "spatialkit_not_supported")
  expect_length(w, 1L)
  expect_identical(conditionMessage(w[[1L]]),
                   paste0("resolution_profile(): the 12 points sit on 1 distinct location, so there ",
                          "is no ladder of cell counts to choose from. The profile scores 2 cells so ",
                          "that level is still visible, and `supported` is FALSE. A tessellation ",
                          "needs points at two locations or more."))
  expect_identical(.cl_fields(w[[1L]]),
                   list(reason = "fewer than two distinct locations", n = 12L,
                        min_cell_n = NA_integer_, obs_per_cell = 4.5))
  expect_false("spatialkit_size_not_reached" %in% .cl_classes(r))
  expect_false(attr(r$value, "bounds")$supported)
  expect_identical(r$value$levels, 2L)
  # With an estimate: no level can be fitted, so `loo` is NA, and nothing is
  # said of what a leave-one-out error describes where rows share a location.
  r <- .cl_catch(.cl_profile(pts, estimate = "pooled", obs_per_cell = 4.5, partition = "kmeans",
                             nstart = 2))
  expect_identical(.cl_classes(r), "spatialkit_not_supported")
  expect_identical(r$value$loo, NA_real_)
  expect_false(attr(r$value, "estimate")$constant)
})


# ---- Strings are read by their value -----------------------------------------

test_that("\"mean\", \"pooled\" and \"auto\" are read by value, whatever name or class the string carries", {
  pts <- .cl_layer(300, 322)
  # One configuration vector, as a caller keeps one: each element carries its name.
  cfg  <- c(estimate = "pooled", size = "auto", sac = "auto")
  bare <- .cl_profile(pts, estimate = "pooled", obs_per_cell = "auto", sac = "auto")
  expect_identical(bare$levels, 67L)
  expect_true(is.finite(bare$loo))
  expect_true(is.finite(bare$own_weight))
  got <- .cl_catch(.cl_profile(pts, estimate = cfg["estimate"], obs_per_cell = cfg["size"],
                               sac = cfg["sac"]))
  expect_length(got$warnings, 0L)
  expect_identical(got$value, bare)
  expect_identical(attr(got$value, "estimate")$value, "pooled")
  expect_identical(attr(got$value, "estimate")$neighbours, 8L)
  expect_identical(attr(got$value, "cp_noise")$source, "none")
  expect_true(any(grepl("pooled toward the 8 nearest cells", utils::capture.output(print(got$value)),
                        fixed = TRUE)))

  # Each argument alone, a name and a class.
  glue <- function(s) structure(s, class = c("glue", "character"))
  for (est in c("pooled", "mean", "auto")) {
    want <- .cl_profile(pts, estimate = est, obs_per_cell = "auto")
    for (given in list(stats::setNames(est, "a"), glue(est), unlist(list(e = est)))) {
      p <- .cl_profile(pts, estimate = given, obs_per_cell = "auto")
      expect_identical(p, want, info = est)
      expect_identical(attr(p, "estimate")$value, if (est == "auto") "pooled" else est, info = est)
      expect_null(attributes(attr(p, "estimate")$value), info = est)
    }
    # The automatic ladder of the estimate: one size for pooled values, seven for plain means.
    expect_identical(nrow(want), if (est == "mean") 7L else 1L, info = est)
    expect_true(all(is.finite(want$loo)), info = est)
  }
  for (auto in list(c(a = "auto"), glue("auto"))) {
    expect_identical(.cl_profile(pts, estimate = "pooled", obs_per_cell = auto), bare)
    expect_identical(.cl_profile(pts, estimate = "pooled", obs_per_cell = 4.5, sac = auto), bare)
  }
  # A size that carries a name is the size.
  expect_identical(.cl_profile(pts, estimate = "pooled", obs_per_cell = c(fine = 4.5)), bare)
  expect_null(names(attr(bare, "bounds")$obs_per_cell))

  # A value that is not one of the strings is refused as before, name or none.
  expect_identical(conditionMessage(.cl_error(.cl_profile(pts, estimate = c(a = "x")))),
                   "resolution_profile(): `estimate` must be NULL, \"mean\", \"pooled\" or \"auto\"; got \"x\".")
  expect_s3_class(.cl_error(.cl_profile(pts, obs_per_cell = c(a = "Auto"))), "error")
  expect_s3_class(.cl_error(.cl_profile(pts, sac = c(a = "none"))), "error")
  expect_s3_class(.cl_error(.cl_profile(pts, estimate = factor("pooled"))), "error")
  expect_true(.is_auto("auto"))
  expect_true(.is_auto(c(a = "auto")))
  expect_false(.is_auto(c("auto", "auto")))
  expect_false(.is_auto(NA_character_))
  expect_false(.is_auto("Auto"))
  expect_false(.is_auto(NULL))
  expect_false(.is_auto(1))
})


# ---- Values the arguments refuse ---------------------------------------------

test_that("obs_per_cell and estimate refuse what they cannot read, in the words of the design", {
  pts <- .cl_600()
  size_text <- function(got)
    sprintf(paste0("resolution_profile(): `obs_per_cell` must be NULL, \"auto\", or numbers of at ",
                   "least 1 (observations per cell; fractions are allowed); got %s."), got)
  bad <- list(list(NA, "NA"), list(NA_real_, "NA"), list(0.5, "0.5"), list(numeric(0), "numeric(0)"),
              list("big", "\"big\""), list(c(4.5, 0.9), "c(4.5, 0.9)"), list(Inf, "Inf"),
              list(c(9, NA), "c(9, NA)"), list(TRUE, "TRUE"), list(list(4.5), "an object of class list"))
  # A `units` object is numeric too; it is refused in the argument's own words
  # and not in those of a comparison inside that package.
  if (requireNamespace("units", quietly = TRUE))
    bad <- c(bad, list(list(units::set_units(4.5, "1"), "an object of class units")))
  for (b in bad) {
    e <- .cl_error(.cl_profile(pts, obs_per_cell = b[[1L]]))
    expect_identical(conditionMessage(e), size_text(b[[2L]]), info = b[[2L]])
    expect_null(conditionCall(e))
  }

  e <- .cl_error(.cl_profile(pts, estimate = "x"))
  expect_identical(conditionMessage(e),
                   "resolution_profile(): `estimate` must be NULL, \"mean\", \"pooled\" or \"auto\"; got \"x\".")
  # No partial matching, and one value.
  for (bad_est in list("pool", "Pooled", c("mean", "pooled"), NA_character_, 1, TRUE))
    expect_match(conditionMessage(.cl_error(.cl_profile(pts, estimate = bad_est))),
                 "`estimate` must be NULL, \"mean\", \"pooled\" or \"auto\"; got ", fixed = TRUE)

  pts$w    <- pts$z + 1
  pts$flag <- pts$z > 0
  pts$none <- NA_real_
  for (est in c("pooled", "mean")) {
    expect_identical(
      conditionMessage(.cl_error(.cl_profile(pts, predictor_vars = "w", estimate = est,
                                             obs_per_cell = NULL))),
      paste0("resolution_profile(): `estimate` cannot be combined with `predictor_vars`. With ",
             "predictors the profile scores the residuals of the response on them, and the ",
             "summary reports region values of the response itself, so the criterion would be ",
             "computed for a number nobody reports. Drop `predictor_vars`, or leave `estimate` ",
             "NULL."))
    expect_identical(
      conditionMessage(.cl_error(.cl_profile(pts, response = NULL, estimate = est,
                                             obs_per_cell = NULL))),
      "resolution_profile(): `estimate` names the value of a response, so it needs a `response_var`.")
    expect_identical(
      conditionMessage(.cl_error(.cl_profile(pts, response = "flag", estimate = est,
                                             obs_per_cell = NULL))),
      paste0("resolution_profile(): `estimate` scores the cell value summarize_by_cell() reports, ",
             "and that function summarises numeric columns only; 'flag' is logical. Convert it ",
             "with as.numeric() first, or leave `estimate` NULL."))
    expect_identical(
      conditionMessage(.cl_error(.cl_profile(pts, response = "none", estimate = est,
                                             obs_per_cell = NULL))),
      "resolution_profile(): `estimate` needs rows with a value of 'none', and none has a finite one.")
  }
  # The checks come in the design's order: a logical column with no value at
  # all is refused as logical, and a numeric one whose only values are
  # infinite as having no finite value.
  pts$lna <- NA
  pts$inf <- rep(c(Inf, NA, -Inf), 200)
  expect_true(is.logical(pts$lna))
  expect_identical(
    conditionMessage(.cl_error(.cl_profile(pts, response = "lna", estimate = "pooled",
                                           obs_per_cell = NULL))),
    paste0("resolution_profile(): `estimate` scores the cell value summarize_by_cell() reports, ",
           "and that function summarises numeric columns only; 'lna' is logical. Convert it ",
           "with as.numeric() first, or leave `estimate` NULL."))
  expect_identical(
    conditionMessage(.cl_error(.cl_profile(pts, response = "inf", estimate = "pooled",
                                           obs_per_cell = NULL))),
    "resolution_profile(): `estimate` needs rows with a value of 'inf', and none has a finite one.")
})

test_that("estimate = \"auto\" is NULL wherever naming an estimate is an error", {
  skip_if_not_installed("gstat")
  pts <- .cl_600()
  pts$w    <- pts$z + sin(seq_len(600))
  pts$flag <- pts$z > 0
  pts$none <- NA_real_
  pts$inf  <- rep(c(Inf, NA, -Inf), 200)
  same <- function(...) {
    a <- suppressWarnings(.cl_profile(pts, ..., estimate = "auto", obs_per_cell = NULL,
                                      sac = NULL, n_levels = 4))
    b <- suppressWarnings(.cl_profile(pts, ..., estimate = NULL, obs_per_cell = NULL,
                                      sac = NULL, n_levels = 4))
    expect_identical(a, b)
    expect_null(attr(a, "estimate"))
    expect_false("loo" %in% names(a))
  }
  same(predictor_vars = "w")
  same(response = NULL)
  same(response = "flag")
  same(response = "none")
  # A column whose only values are infinite has no finite value either.
  same(response = "inf")
  # And where it is not an error, "auto" is the pooled value.
  a <- suppressWarnings(.cl_profile(pts, estimate = "auto", obs_per_cell = NULL, sac = NULL,
                                    n_levels = 4))
  expect_identical(attr(a, "estimate")$value, "pooled")
  expect_identical(a, suppressWarnings(.cl_profile(pts, estimate = "pooled", obs_per_cell = NULL,
                                                   sac = NULL, n_levels = 4)))
})

test_that("the present checks of the ladder's arguments run beside a size", {
  pts <- .cl_600()
  expect_identical(conditionMessage(.cl_error(.cl_profile(pts, obs_per_cell = 4.5, n_levels = 1))),
                   "resolution_profile(): `n_levels` must be a single number >= 2.")
  expect_identical(conditionMessage(.cl_error(.cl_profile(pts, obs_per_cell = 4.5, sample_n = 2))),
                   "resolution_profile(): `sample_n` must be a single number >= 3.")
  expect_identical(conditionMessage(.cl_error(.cl_profile(pts, obs_per_cell = 4.5, min_cell_n = 0))),
                   "resolution_profile(): `min_cell_n` must be a single number >= 1.")
  # They come before the checks of the new arguments.
  expect_identical(conditionMessage(.cl_error(.cl_profile(pts, obs_per_cell = "big", estimate = "x",
                                                          n_levels = 1))),
                   "resolution_profile(): `n_levels` must be a single number >= 2.")
})


# ---- nstart = NULL -----------------------------------------------------------

test_that("nstart = NULL is the kind's own number of starts", {
  skip_if_not_installed("gstat")
  pts <- .cl_layer(300, 310)
  for (kind in c("kmeans", "density")) {
    own <- if (kind == "kmeans") 25L else 3L
    a <- suppressWarnings(.cl_profile(pts, partition = kind, nstart = NULL, estimate = NULL,
                                      obs_per_cell = NULL, sac = NULL, levels = c(6L, 30L)))
    b <- suppressWarnings(.cl_profile(pts, partition = kind, nstart = own, estimate = NULL,
                                      obs_per_cell = NULL, sac = NULL, levels = c(6L, 30L)))
    expect_identical(a, b, info = kind)
    expect_identical(attr(a, "nstart"), own, info = kind)
    # Another number is another profile, so the comparison can fail.
    other <- suppressWarnings(.cl_profile(pts, partition = kind, nstart = 1, estimate = NULL,
                                          obs_per_cell = NULL, sac = NULL, levels = c(6L, 30L)))
    expect_false(identical(a$wss_spread, other$wss_spread), info = kind)
    # And on a ladder of sizes.
    s <- .cl_profile(pts, partition = kind, nstart = NULL, obs_per_cell = 9)
    expect_identical(attr(s, "nstart"), own, info = kind)
    expect_identical(s, .cl_profile(pts, partition = kind, nstart = own, obs_per_cell = 9), info = kind)
  }
  expect_identical(.NSTART_DENSITY, 3L)
  for (bad in list(0, "many", NA_real_, c(2, 3), Inf))
    expect_identical(conditionMessage(.cl_error(.cl_profile(pts, nstart = bad))),
                     "resolution_profile(): `nstart` must be NULL or a single number >= 1.")
})

test_that("the two present errors whose argument gains a value name that value, and say nothing else new", {
  # `nstart` gains NULL and `sac` gains "auto", and the error of each for a
  # value it cannot read lists what it can: the present text with the new
  # value added.  These are the only two texts of a call that names no new
  # argument that are not the baseline's, and the frozen fixture holds no
  # failing call, so they are held here in full, each beside the sentence it
  # was before.
  pts <- .cl_600()
  nstart_text <- "resolution_profile(): `nstart` must be NULL or a single number >= 1."
  for (bad in list(0, "many", c(1, 2), NA, NA_real_, Inf, -3, TRUE, list(5))) {
    e <- .cl_error(.cl_profile(pts, nstart = bad))
    expect_identical(conditionMessage(e), nstart_text)
    expect_identical(class(e), c("simpleError", "error", "condition"))
    expect_null(conditionCall(e))
  }
  expect_identical(sub("NULL or ", "", nstart_text, fixed = TRUE),
                   "resolution_profile(): `nstart` must be a single number >= 1.")

  sac_text <- function(got)
    paste0("resolution_profile(): `sac` must be NULL, \"auto\", an estimate_sac_range() result, ",
           "or a plain number in the units of the CRS the profile is computed in (metres for ",
           "lon/lat input); got ", got, ".")
  bad <- list(list("none", "an object of class character"),
              list(TRUE, "an object of class logical"),
              list(list(1), "an object of class list"),
              list(c("auto", "auto"), "an object of class character"),
              list(factor("auto"), "an object of class factor"))
  if (requireNamespace("units", quietly = TRUE))
    bad <- c(bad, list(list(units::set_units(1.5, "km"), "1.5 [km]")))
  for (b in bad) {
    e <- .cl_error(.cl_profile(pts, sac = b[[1L]]))
    expect_identical(conditionMessage(e), sac_text(b[[2L]]), info = b[[2L]])
    expect_identical(class(e), c("simpleError", "error", "condition"))
    expect_null(conditionCall(e))
  }
  expect_identical(sub("NULL, \"auto\", ", "", sac_text("an object of class list"), fixed = TRUE),
                   paste0("resolution_profile(): `sac` must be an estimate_sac_range() result, or a ",
                          "plain number in the units of the CRS the profile is computed in (metres ",
                          "for lon/lat input); got an object of class list."))
})


# ---- sac = "auto" ------------------------------------------------------------

test_that("sac = \"auto\" is NULL on the automatic ladder and on a `levels` ladder", {
  skip_if_not_installed("gstat")
  pts <- .cl_layer(400, 311)
  for (est in list(NULL, "pooled")) {
    expect_identical(suppressWarnings(.cl_profile(pts, estimate = est, obs_per_cell = NULL,
                                                  sac = "auto", n_levels = 5)),
                     suppressWarnings(.cl_profile(pts, estimate = est, obs_per_cell = NULL,
                                                  sac = NULL, n_levels = 5)))
    expect_identical(suppressWarnings(.cl_profile(pts, estimate = est, obs_per_cell = NULL,
                                                  sac = "auto", levels = c(5L, 25L))),
                     suppressWarnings(.cl_profile(pts, estimate = est, obs_per_cell = NULL,
                                                  sac = NULL, levels = c(5L, 25L))))
  }
  # The variogram was estimated there.
  a <- suppressWarnings(.cl_profile(pts, estimate = NULL, obs_per_cell = NULL, sac = "auto",
                                    n_levels = 5))
  expect_s3_class(attr(a, "sac"), "sac_range")
  expect_true(all(is.finite(a$cp)))
})

test_that("sac = \"auto\" fits no variogram on a ladder of sizes", {
  skip_if_not_installed("gstat")
  pts  <- .cl_layer(900, 312)
  with <- suppressWarnings(.cl_profile(pts, obs_per_cell = c(18, 4.5), sac = NULL))
  own  <- suppressWarnings(.cl_profile(pts, obs_per_cell = c(18, 4.5), sac = 20000))

  called <- 0L
  logged <- NULL
  r <- NULL
  local({
    testthat::local_mocked_bindings(estimate_sac_range = function(...) {
      called <<- called + 1L
      stop("estimate_sac_range() was reached")
    })
    logged <<- capture_spatialkit_log(
      r <<- .cl_catch(.cl_profile(pts, obs_per_cell = c(18, 4.5), sac = "auto")))
  })
  expect_identical(called, 0L)
  expect_length(r$warnings, 0L)
  expect_false(any(grepl("variogram|estimate_sac_range|`sac`", logged)))
  none <- r$value

  expect_identical(none$levels, c(50L, 200L))
  expect_true(all(is.na(none$cp)))
  expect_true(all(is.na(none$reliability)))
  expect_identical(attr(none, "cp_noise"),
                   list(source = "none", value = NA_real_, level = NA_integer_,
                        reason = "no variogram: a size ladder fits one only when `sac` is given"))
  expect_null(attr(none, "variogram"))
  expect_null(attr(none, "sac"))
  expect_identical(attr(none, "bounds")$range, NA_real_)
  expect_identical(attr(none, "bounds")$floor, 2L)
  one <- attr(none, "one_cell")
  expect_identical(one$cp, NA_real_)
  expect_identical(unname(one$cp_levels), c(NA_real_, NA_real_))
  expect_identical(names(one$cp_levels), c("50", "200"))
  expect_true(all(is.finite(one$f_p)))
  expect_identical(one[c("n", "rss", "constant", "f_p")],
                   attr(with, "one_cell")[c("n", "rss", "constant", "f_p")])

  # No other column reads the variogram, and the variogram step leaves the
  # random-number stream where it found it, so the fits draw the same starts.
  for (col in c("levels", "wss", "wss_spread", "elbow", "cell_n_min", "cell_n_median",
                "cell_diam_median", "rss", "moran_i", "moran_z", "cells_resp", "cell_resp_mean",
                "scored_per_cell", "loo", "loo_se", "share_alone", "own_weight"))
    expect_identical(none[[col]], with[[col]], info = col)
  expect_identical(attr(none, "centres"), attr(with, "centres"))
  expect_identical(attr(none, "estimate"), attr(with, "estimate"))

  # With sac = NULL, and with a `sac` of the caller's, `cp` is there.
  expect_true(all(is.finite(with$cp)))
  expect_true(all(is.finite(own$cp)))
  expect_false(identical(attr(with, "cp_noise")$source, "none"))
  expect_false(identical(attr(own, "cp_noise")$source, "none"))
  expect_identical(attr(own, "bounds")$range, 20000)

  out <- utils::capture.output(print(none))
  expect_true(paste0("  variogram   : none fitted: a size ladder fits one only when `sac` is given ",
                     "(cp and reliability are NA)") %in% out)
  expect_false(any(grepl("none fitted", utils::capture.output(print(with)), fixed = TRUE)))
})

test_that("sac = \"auto\" is not a variogram the caller supplied", {
  pts <- .cl_layer(400, 313)
  line <- "select_on = \"split\" with a supplied `sac`"
  auto <- capture_spatialkit_log(
    p <- suppressWarnings(.cl_profile(pts, obs_per_cell = 4.5, sac = "auto", select_on = "split")))
  expect_false(any(grepl(line, auto, fixed = TRUE)))
  expect_identical(attr(p, "cp_noise")$source, "none")
  expect_null(attr(p, "sac"))
  # A number is one, and the line is written for it.
  given <- capture_spatialkit_log(
    suppressWarnings(.cl_profile(pts, obs_per_cell = 4.5, sac = 20000, select_on = "split")))
  expect_true(any(grepl(line, given, fixed = TRUE)))
})

test_that("a `sac` that is another string is refused, and the message names \"auto\"", {
  pts <- .cl_600()
  e <- .cl_error(.cl_profile(pts, sac = "none"))
  expect_identical(conditionMessage(e),
                   paste0("resolution_profile(): `sac` must be NULL, \"auto\", an ",
                          "estimate_sac_range() result, or a plain number in the units of the CRS ",
                          "the profile is computed in (metres for lon/lat input); got an object of ",
                          "class character."))
  expect_null(conditionCall(e))
  expect_s3_class(.cl_error(.cl_profile(pts, sac = c("auto", "auto"))), "error")
  expect_s3_class(.cl_error(.cl_profile(pts, sac = "Auto")), "error")
})


# ---- Moran's z ---------------------------------------------------------------

test_that("grouped sums give the cell means of the loop", {
  d <- local({
    .cl_seed(314)
    n  <- 2000L
    xy <- cbind(stats::runif(n, 0, 1e5), stats::runif(n, 0, 1e5))
    list(xy = xy, y = sin(xy[, 1] / 2e4) + stats::rnorm(n),
         p = cbind(stats::rnorm(n), xy[, 2] / 1e5 + stats::rnorm(n)),
         cl = sample(c(3L, 7L, 11L:160L), n, replace = TRUE))
  })
  none <- matrix(numeric(0), nrow(d$xy), 0L)
  for (pred in list(none, d$p)) {
    loop    <- .morans_i_for_k(d$xy, d$y, pred, d$cl)
    grouped <- .morans_i_for_k(d$xy, d$y, pred, d$cl, grouped = TRUE)
    expect_true(all(is.finite(loop)))
    expect_equal(grouped, loop, tolerance = 1e-10)
    expect_identical(names(grouped), c("I", "z"))
  }
  # A missing response is left out of its cell's mean by both.
  y <- d$y
  y[c(5L, 90L, 91L, 1500L)] <- NA
  expect_equal(.morans_i_for_k(d$xy, y, none, d$cl, grouped = TRUE),
               .morans_i_for_k(d$xy, y, none, d$cl), tolerance = 1e-10)
  # Fewer than four cells: not computable, either way.
  expect_identical(.morans_i_for_k(d$xy, d$y, none, rep(1:3, length.out = 2000), grouped = TRUE),
                   c(I = NA_real_, z = NA_real_))
  # The default is the loop, whose numbers the present ladder keeps.
  expect_identical(formals(.morans_i_for_k)$grouped, FALSE)
})

test_that("on a ladder of sizes Moran's z is the present loop's on the same cells", {
  skip_if_not_installed("gstat")
  pts <- .cl_layer(1200, 315)
  for (kind in c("density", "kmeans")) {
    size <- suppressWarnings(.cl_profile(pts, partition = kind, obs_per_cell = c(18, 4.5),
                                         sac = NULL))
    # The same counts as `levels`: the same fits, from the same stream, scored
    # by the loop.
    present <- suppressWarnings(.cl_profile(pts, partition = kind, estimate = NULL,
                                            obs_per_cell = NULL, sac = NULL, sample_n = 1500,
                                            levels = size$levels))
    expect_identical(size$wss, present$wss, info = kind)
    expect_identical(attr(size, "centres"), attr(present, "centres"), info = kind)
    expect_true(all(is.finite(size$moran_z)), info = kind)
    expect_equal(size$moran_i, present$moran_i, tolerance = 1e-10, info = kind)
    expect_equal(size$moran_z, present$moran_z, tolerance = 1e-10, info = kind)
  }
})

test_that("the profile hands its two backends and the grouped form to the Moran step", {
  skip_if_not_installed("gstat")
  pts  <- .cl_layer(300, 316)
  real <- .morans_i_for_k
  force(real)
  seen <- list()
  spy  <- function(xy, response, predictors, cluster_ids, grouped = FALSE, use_fnn = NA,
                   use_matrix = NA) {
    seen[[length(seen) + 1L]] <<- list(grouped = grouped, use_fnn = use_fnn,
                                       use_matrix = use_matrix)
    real(xy, response, predictors, cluster_ids, grouped = grouped)
  }
  local({
    testthat::local_mocked_bindings(.morans_i_for_k = spy)
    .cl_profile(pts, obs_per_cell = c(18, 4.5))
  })
  expect_length(seen, 2L)
  expect_true(all(vapply(seen, function(s) isTRUE(s$grouped), NA)))
  expect_identical(seen[[1L]]$use_fnn, requireNamespace("FNN", quietly = TRUE))
  expect_identical(seen[[1L]]$use_matrix, requireNamespace("Matrix", quietly = TRUE))
  # The form goes with the ladder and not with the estimate: a ladder of
  # sizes without one takes the grouped sums too.
  seen <- list()
  local({
    testthat::local_mocked_bindings(.morans_i_for_k = spy)
    .cl_profile(pts, estimate = NULL, obs_per_cell = c(18, 4.5))
  })
  expect_length(seen, 2L)
  expect_true(all(vapply(seen, function(s) isTRUE(s$grouped), NA)))

  # A ladder the function builds, and one given as `levels`: the loop.
  seen <- list()
  local({
    testthat::local_mocked_bindings(.morans_i_for_k = spy)
    suppressWarnings(.cl_profile(pts, obs_per_cell = NULL, sac = NULL, levels = c(12L, 30L)))
  })
  expect_length(seen, 2L)
  expect_true(all(vapply(seen, function(s) isFALSE(s$grouped), NA)))

  # With a backend switched off the step is told so, on either ladder.
  for (pkg in c("FNN", "Matrix")) {
    seen <- list()
    local({
      testthat::local_mocked_bindings(.morans_i_for_k = spy)
      .cl_without(pkg)
      suppressWarnings(.cl_profile(pts, partition = "kmeans", obs_per_cell = NULL, sac = NULL,
                                   levels = c(12L, 30L)))
    })
    expect_identical(seen[[1L]]$use_fnn, pkg != "FNN" && requireNamespace("FNN", quietly = TRUE),
                     info = pkg)
    expect_identical(seen[[1L]]$use_matrix,
                     pkg != "Matrix" && requireNamespace("Matrix", quietly = TRUE), info = pkg)
  }
  expect_identical(.KNN_DENSE_MAX, 5000L)
})

# The limit of 5,000 cells lowered to a number a small layer reaches, so that
# the rule runs wherever the tests do.  The 23,000-row test below holds it at
# the limit itself.
test_that("Moran's z is left out above the limit and not at it, and the warning says where", {
  skip_if_not_installed("FNN")
  skip_if_not_installed("Matrix")
  # 720 rows at 4.5 a cell: 160 cells, from one random k-means start.
  pts  <- .cl_layer(720, 323)
  size <- function(limit, ..., layer = pts, off = "FNN", obs_per_cell = 4.5) {
    testthat::local_mocked_bindings(.KNN_DENSE_MAX = limit)
    if (length(off)) .cl_without(off)
    .cl_catch(.cl_profile(layer, ..., obs_per_cell = obs_per_cell, partition = "kmeans",
                          init = "random", nstart = 1))
  }
  text <- function(limit, reach, lacking)
    sprintf(paste0("resolution_profile(): `moran_i` and `moran_z` will be NA at a level where more ",
                   "than %d cells hold a scored row (this ladder reaches %d cells), because %s not ",
                   "installed: Moran's z of the cell means needs FNN to find each cell's neighbours ",
                   "and Matrix to hold them. Every other column is computed."), limit, reach, lacking)
  both <- size(159L, off = character(0))
  expect_length(both$warnings, 0L)
  expect_identical(both$value$levels, 160L)
  expect_identical(both$value$cells_resp, 160L)
  expect_true(is.finite(both$value$moran_z))

  # At the limit the step runs and nothing is said; one cell above it the two
  # columns are NA, with one warning that names the package that is missing.
  at <- size(160L)
  expect_length(at$warnings, 0L)
  expect_equal(at$value$moran_z, both$value$moran_z, tolerance = 1e-10)
  for (pkg in c("FNN", "Matrix")) {
    above <- size(159L, off = pkg)
    expect_identical(.cl_classes(above), "spatialkit_moran_skipped", info = pkg)
    w <- above$warnings[[1L]]
    expect_identical(conditionMessage(w), text(159L, 160L, paste(pkg, "is")), info = pkg)
    expect_identical(.cl_fields(w), list(limit = 159L, levels = 160L, missing = pkg), info = pkg)
    expect_null(conditionCall(w))
    expect_identical(above$value$moran_i, NA_real_, info = pkg)
    expect_identical(above$value$moran_z, NA_real_, info = pkg)
    for (col in setdiff(names(both$value), c("moran_i", "moran_z")))
      expect_identical(above$value[[col]], both$value[[col]], info = paste(pkg, col))
  }
  above <- size(159L, off = c("FNN", "Matrix"))
  expect_identical(conditionMessage(above$warnings[[1L]]), text(159L, 160L, "FNN and Matrix are"))
  expect_identical(above$warnings[[1L]]$missing, c("FNN", "Matrix"))

  # Three levels, two of them above the limit: the condition lists those two,
  # and the level under it keeps its Moran's z.
  r <- size(100L, obs_per_cell = c(18, 6, 4.5))
  expect_identical(r$value$levels, c(40L, 120L, 160L))
  expect_identical(.cl_classes(r), "spatialkit_moran_skipped")
  expect_identical(r$warnings[[1L]]$levels, c(120L, 160L))
  expect_identical(is.finite(r$value$moran_z), c(TRUE, FALSE, FALSE))
  expect_identical(is.finite(r$value$moran_i), c(TRUE, FALSE, FALSE))
  expect_true(all(is.finite(r$value$loo)))
  # A level at the limit beside one above it keeps its Moran's z as well.
  r <- size(160L, obs_per_cell = c(4.5, 3))
  expect_identical(r$value$levels, c(160L, 240L))
  expect_identical(r$value$cells_resp[1L], 160L)
  expect_identical(.cl_classes(r), "spatialkit_moran_skipped")
  expect_identical(r$warnings[[1L]]$levels, 240L)
  expect_identical(is.finite(r$value$moran_z), c(TRUE, FALSE))

  # The rule goes with the ladder and not with the estimate.
  r <- size(159L, estimate = NULL)
  expect_identical(.cl_classes(r), "spatialkit_moran_skipped")
  expect_identical(r$value$moran_z, NA_real_)

  # Where the step never runs nothing is said: no response, and a constant one.
  r <- size(100L, response = NULL, estimate = NULL)
  expect_length(r$warnings, 0L)
  flat <- pts
  flat$z <- 0.1
  r <- size(100L, layer = flat)
  expect_identical(.cl_classes(r), "spatialkit_constant_response")

  # The step is left out by the cells that hold a scored row, which are known
  # after a level's fit.  A response on the western half alone, at 2.25 rows a
  # cell: 160 centres over every point, about half of them with a scored row.
  # The warning is raised, since such a level can occur, and the two columns
  # are computed, since it did not.
  west <- pts
  west$z[sf::st_coordinates(pts)[, 1] > 5e4] <- NA
  n_west <- sum(is.finite(west$z))
  r <- suppressWarnings(size(1e6L, layer = west, obs_per_cell = n_west / 160))
  expect_identical(r$value$levels, 160L)
  expect_lt(r$value$cells_resp, 120L)
  r <- size(120L, layer = west, obs_per_cell = n_west / 160)
  w <- .cl_of(r, "spatialkit_moran_skipped")
  expect_length(w, 1L)
  expect_identical(conditionMessage(w[[1L]]), text(120L, 160L, "FNN is"))
  expect_true(is.finite(r$value$moran_z))
  expect_true(is.finite(r$value$moran_i))

  # Fewer scored rows than the limit: no level can hold that many cells with
  # one, so nothing is said.  Under the split the selection half alone is
  # scored, and at 1.5 rows a cell the ladder has more cells than that half
  # has rows.
  half <- suppressWarnings(size(1e6L, select_on = "split", obs_per_cell = 1.5))
  rows <- attr(half$value, "bounds")$n_scored
  expect_identical(half$value$levels, 480L)
  expect_lt(rows, 480L)
  expect_gt(rows, 200L)
  r <- size(rows, select_on = "split", obs_per_cell = 1.5)
  expect_length(.cl_of(r, "spatialkit_moran_skipped"), 0L)
  expect_true(is.finite(r$value$moran_z))
  r <- size(rows - 1L, select_on = "split", obs_per_cell = 1.5)
  expect_length(.cl_of(r, "spatialkit_moran_skipped"), 1L)

  # A call a guard stops says nothing of a column it never computes.
  r <- local({
    withr::local_options(spatialkit.fit_work_max = 1)
    testthat::local_mocked_bindings(.KNN_DENSE_MAX = 100L)
    .cl_without("FNN")
    got <- list()
    e <- tryCatch(withCallingHandlers(
      .cl_profile(pts, obs_per_cell = 4.5, partition = "kmeans", init = "random", nstart = 1),
      warning = function(w) { got[[length(got) + 1L]] <<- w; invokeRestart("muffleWarning") }),
      error = function(e) e)
    list(error = e, warnings = got)
  })
  expect_s3_class(r$error, "spatialkit_work_limit")
  expect_identical(r$error$what, "partition fit")
  expect_length(r$warnings, 0L)
})

test_that("a present ladder above the limit stops as it always has, until the rule holds on every ladder", {
  skip_if_not_installed("FNN")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("gstat")
  pts <- .cl_layer(720, 323)
  present <- function(any_ladder, pkg) {
    testthat::local_mocked_bindings(.KNN_DENSE_MAX = 159L,
                                    .moran_skip_any_ladder = function() any_ladder)
    .cl_without(pkg)
    got <- list()
    value <- tryCatch(withCallingHandlers(
      .cl_profile(pts, levels = 160, sample_n = 720, estimate = NULL, obs_per_cell = NULL,
                  sac = NULL, partition = "kmeans", init = "random", nstart = 1),
      warning = function(w) { got[[length(got) + 1L]] <<- w; invokeRestart("muffleWarning") }),
      error = function(e) e)
    list(value = value, warnings = got)
  }
  for (pkg in c("FNN", "Matrix")) {
    # The rule not applied: the unclassed error of .build_knn_weights().
    r <- present(FALSE, pkg)
    expect_identical(class(r$value), c("simpleError", "error", "condition"), info = pkg)
    expect_match(conditionMessage(r$value), "n = 160 requires FNN for k-NN weights", fixed = TRUE)
    expect_length(.cl_of(r, "spatialkit_moran_skipped"), 0L)
    # The rule applied: the profile returns, with the warning and the two columns NA.
    r <- present(TRUE, pkg)
    expect_s3_class(r$value, "resolution_profile")
    w <- .cl_of(r, "spatialkit_moran_skipped")
    expect_length(w, 1L)
    expect_identical(.cl_fields(w[[1L]]), list(limit = 159L, levels = 160L, missing = pkg), info = pkg)
    expect_identical(r$value$levels, 160L)
    expect_identical(r$value$moran_i, NA_real_, info = pkg)
    expect_identical(r$value$moran_z, NA_real_, info = pkg)
    expect_true(is.finite(r$value$wss))
    expect_true(is.finite(r$value$rss))
  }
})

test_that("above 5,000 cells Moran's z is left out without FNN or Matrix, and the profile returns", {
  skip_on_cran()
  skip_if_not_installed("FNN")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("gstat")
  # 23,000 rows at 4.5 a cell: 5,111 cells.  The kind is k-means from one
  # random start, which the fit guard lets through without FNN.
  big  <- .cl_layer(23000, 31)
  size <- function() .cl_profile(big, estimate = "pooled", obs_per_cell = 4.5, partition = "kmeans",
                                 init = "random", nstart = 1, sac = "auto", seed = 123)
  both <- .cl_catch(size())
  expect_length(.cl_of(both, "spatialkit_moran_skipped"), 0L)
  expect_identical(both$value$levels, 5111L)
  expect_gt(both$value$cells_resp, 5000L)
  expect_true(is.finite(both$value$moran_i))
  expect_true(is.finite(both$value$moran_z))
  expect_true(is.finite(both$value$loo))

  real <- .partition_centres
  force(real)
  for (pkg in c("FNN", "Matrix")) {
    # The warnings and the fits, in the order they come.
    events <- character(0)
    r <- NULL
    local({
      .cl_without(pkg)
      testthat::local_mocked_bindings(.partition_centres = function(...) {
        events <<- c(events, "fit")
        real(...)
      })
      got <- list()
      value <- withCallingHandlers(size(), warning = function(w) {
        events <<- c(events, class(w)[1L])
        got[[length(got) + 1L]] <<- w
        invokeRestart("muffleWarning")
      })
      r <<- list(value = value, warnings = got)
    })
    w <- .cl_of(r, "spatialkit_moran_skipped")
    expect_length(w, 1L)
    expect_identical(.cl_classes(r), "spatialkit_moran_skipped", info = pkg)
    # Said once, before the first fit.
    expect_identical(events, c("spatialkit_moran_skipped", "fit"), info = pkg)
    expect_identical(conditionMessage(w[[1L]]),
                     sprintf(paste0("resolution_profile(): `moran_i` and `moran_z` will be NA at a ",
                                    "level where more than 5000 cells hold a scored row (this ladder ",
                                    "reaches 5111 cells), because %s is not installed: Moran's z of ",
                                    "the cell means needs FNN to find each cell's neighbours and ",
                                    "Matrix to hold them. Every other column is computed."), pkg))
    expect_identical(.cl_fields(w[[1L]]), list(limit = 5000L, levels = 5111L, missing = pkg))
    expect_identical(r$value$moran_i, NA_real_, info = pkg)
    expect_identical(r$value$moran_z, NA_real_, info = pkg)
    expect_true(is.finite(r$value$loo), info = pkg)
    expect_identical(r$value$loo, both$value$loo, info = pkg)
    # Every other column is the run's with both.
    for (col in setdiff(names(both$value), c("moran_i", "moran_z")))
      expect_identical(r$value[[col]], both$value[[col]], info = paste(pkg, col))
  }
  # Both missing: one warning that names both.
  r <- local({ .cl_without(c("FNN", "Matrix")); .cl_catch(size()) })
  w <- .cl_of(r, "spatialkit_moran_skipped")
  expect_length(w, 1L)
  expect_match(conditionMessage(w[[1L]]), "because FNN and Matrix are not installed", fixed = TRUE)
  expect_identical(w[[1L]]$missing, c("FNN", "Matrix"))
  expect_identical(r$value$loo, both$value$loo)

  # The step is left out by the cells that hold a scored row, not by the
  # centres.  A response on the western half alone, at 2.25 rows a cell: more
  # than 5,000 centres over every point, about half of them with a scored row.
  # The warning says where the columns will be NA, which is all that is known
  # before the fits; here no level is such a level, and they are computed.
  west <- big
  west$z[sf::st_coordinates(big)[, 1] > 5e4] <- NA
  r <- local({
    .cl_without("FNN")
    .cl_catch(.cl_profile(west, estimate = "pooled", obs_per_cell = 2.25, partition = "kmeans",
                          init = "random", nstart = 1, sac = "auto", seed = 123))
  })
  expect_gt(r$value$levels, 5000L)
  expect_lt(r$value$cells_resp, 4000L)
  expect_length(.cl_of(r, "spatialkit_moran_skipped"), 1L)
  expect_true(is.finite(r$value$moran_i))
  expect_true(is.finite(r$value$moran_z))

  # A present ladder at such a level stops with the error it has always had
  # while .moran_skip_any_ladder() is FALSE ...  Both values are set here, so
  # the test does not read the one the function holds.
  present <- function() .cl_profile(big, levels = 5200, sample_n = 23000, estimate = NULL,
                                    obs_per_cell = NULL, sac = NULL, partition = "kmeans",
                                    init = "random", nstart = 1, seed = 123)
  for (pkg in c("FNN", "Matrix")) {
    e <- local({
      .cl_without(pkg)
      testthat::local_mocked_bindings(.moran_skip_any_ladder = function() FALSE)
      .cl_error(present())
    })
    expect_identical(class(e), c("simpleError", "error", "condition"), info = pkg)
    expect_match(conditionMessage(e), "n = 5200 requires FNN for k-NN weights", fixed = TRUE)
  }
  # ... and, once it says the rule holds on every ladder, returns.
  r <- local({
    .cl_without("FNN")
    testthat::local_mocked_bindings(.moran_skip_any_ladder = function() TRUE)
    .cl_catch(present())
  })
  expect_length(.cl_of(r, "spatialkit_moran_skipped"), 1L)
  expect_identical(.cl_of(r, "spatialkit_moran_skipped")[[1L]]$levels, 5200L)
  expect_identical(r$value$levels, 5200L)
  expect_identical(r$value$moran_i, NA_real_)
  expect_identical(r$value$moran_z, NA_real_)
  expect_true(is.finite(r$value$wss))
  expect_true(is.finite(r$value$rss))
})

test_that("Moran's z is computed below 5,001 cells whatever is installed, and nothing is said", {
  pts <- .cl_layer(600, 317)
  with <- .cl_catch(.cl_profile(pts, obs_per_cell = c(18, 4.5)))
  for (pkg in list("FNN", "Matrix", c("FNN", "Matrix"))) {
    r <- local({ .cl_without(pkg); .cl_catch(.cl_profile(pts, obs_per_cell = c(18, 4.5))) })
    expect_length(r$warnings, 0L)
    expect_true(all(is.finite(r$value$moran_z)))
    expect_equal(r$value$moran_i, with$value$moran_i, tolerance = 1e-10)
    expect_equal(r$value$moran_z, with$value$moran_z, tolerance = 1e-10)
  }
})


# ---- The fit guard -----------------------------------------------------------

# 1,000 points with a response on 500: at 4.5 rows a cell, 111 cells.
.cl_half <- function() .cl_keep("half", {
  pts <- .cl_layer(1000, 318)
  pts$z[seq(2, 1000, by = 2)] <- NA
  pts
})

test_that("the fit guard stops a ladder of sizes whose fits compare too many points with centres", {
  skip_if_not_installed("FNN")
  pts  <- .cl_half()
  stop_fit <- function(...) stop("a fit was reached")
  guard <- function(limit, ..., obs_per_cell = 4.5) {
    withr::local_options(spatialkit.fit_work_max = limit)
    testthat::local_mocked_bindings(.partition_centres = stop_fit)
    .cl_error(.cl_profile(pts, obs_per_cell = obs_per_cell, ...))
  }
  cases <- list(list(kind = "kmeans",  init = "kmeans++", nstart = 25L, with = "25 k-means++ starts"),
                list(kind = "kmeans",  init = "random",   nstart = 25L, with = "25 random starts"),
                list(kind = "density", init = "kmeans++", nstart = 3L,
                     with = "3 density-weighted k-means++ starts"))
  for (cs in cases) {
    what <- paste(cs$kind, cs$init)
    # The work is the points, with a response or without, times the cells of
    # the finest level, times the starts: `nstart = NULL` is the kind's own.
    work <- 1000 * 111 * cs$nstart
    e <- guard(work - 1, partition = cs$kind, init = cs$init, nstart = NULL)
    expect_s3_class(e, "spatialkit_work_limit")
    expect_identical(class(e), c("spatialkit_work_limit", "error", "condition"), info = what)
    expect_null(conditionCall(e))
    expect_identical(.cl_fields(e),
                     list(work = work, limit = work - 1, what = "partition fit", kind = cs$kind,
                          init = cs$init, nstart = cs$nstart, fnn = TRUE), info = what)
    expect_identical(conditionMessage(e),
                     sprintf(paste0("resolution_profile(): fitting 111 cells to 1000 rows with %s is ",
                                    "%.2g point-to-centre comparisons, above the limit of %.2g ",
                                    "(options(spatialkit.fit_work_max = )). Ways on: fewer starts ",
                                    "(`nstart`), which keep the best of fewer; a coarser ",
                                    "`obs_per_cell`; partition = \"density\" with init = \"random\" and ",
                                    "the FNN package installed, which has no step that compares every ",
                                    "point with every centre and leaves a cell without a point more ",
                                    "often; or a higher limit."), cs$with, work, work - 1),
                     info = what)
    # At the limit the fits run: the stand-in for the fit is what stops the call.
    e <- guard(work, partition = cs$kind, init = cs$init, nstart = NULL)
    expect_identical(conditionMessage(e), "a fit was reached", info = what)
    # A number of starts the caller gives is the number counted.
    e <- guard(1000 * 111 * 2 - 1, partition = cs$kind, init = cs$init, nstart = 2)
    expect_identical(e$work, 1000 * 111 * 2, info = what)
    expect_identical(e$nstart, 2L, info = what)
  }
  # One start is written in the singular.
  e <- guard(1, partition = "kmeans", init = "random", nstart = 1)
  expect_match(conditionMessage(e), "with 1 random start is", fixed = TRUE)

  # The finest level of several is the one counted.
  e <- guard(1, partition = "density", init = "kmeans++", nstart = 3, obs_per_cell = c(18, 9, 2))
  expect_identical(e$work, 1000 * 250 * 3)

  # The density kind from random starts has no such step with FNN: nothing
  # is counted at any limit.
  e <- guard(0, partition = "density", init = "random", nstart = 3)
  expect_identical(conditionMessage(e), "a fit was reached")

  # The guard bounds the fits, whatever is scored on them: without an
  # estimate, and without a response, where the size divides every point.
  e <- guard(1000 * 111 * 2 - 1, partition = "kmeans", init = "kmeans++", nstart = 2,
             estimate = NULL)
  expect_s3_class(e, "spatialkit_work_limit")
  expect_identical(e$work, 1000 * 111 * 2)
  e <- guard(1000 * 222 * 3 - 1, partition = "density", init = "kmeans++", nstart = 3,
             response = NULL, estimate = NULL)
  expect_s3_class(e, "spatialkit_work_limit")
  expect_identical(e$work, 1000 * 222 * 3)
  expect_identical(conditionMessage(guard(1000 * 222 * 3, partition = "density", init = "kmeans++",
                                          nstart = 3, response = NULL, estimate = NULL)),
                   "a fit was reached")
})

test_that("a call a guard stops has estimated no variogram and raised no warning", {
  skip_if_not_installed("gstat")
  pts <- .cl_layer(300, 326)
  # With sac = NULL a ladder of sizes estimates the variogram, whose warnings
  # are about `cp` and `reliability`: columns of a profile that a stopped
  # call never returns.  The guards come first.
  stopped <- function(..., fit_max = 1, pairs_max = NULL, off = character(0)) {
    reached <- 0L
    withr::local_options(spatialkit.fit_work_max = fit_max, spatialkit.exact_pairs_max = pairs_max)
    testthat::local_mocked_bindings(
      estimate_sac_range = function(...) { reached <<- reached + 1L; stop("the variogram was reached") },
      .partition_centres = function(...) stop("a fit was reached"))
    if (length(off)) .cl_without(off)
    got <- list()
    e <- tryCatch(withCallingHandlers(
      .cl_profile(pts, ..., sac = NULL),
      warning = function(w) { got[[length(got) + 1L]] <<- w; invokeRestart("muffleWarning") }),
      error = function(e) e)
    list(error = e, warnings = got, variogram = reached)
  }
  # The fit guard, with and without an estimate.
  for (est in list("pooled", NULL)) {
    r <- stopped(estimate = est, obs_per_cell = 4.5, partition = "kmeans")
    expect_s3_class(r$error, "spatialkit_work_limit")
    expect_identical(r$error$what, "partition fit")
    expect_identical(r$variogram, 0L)
    expect_length(r$warnings, 0L)
  }
  # A size that would be held to two cells is not reported either.
  r <- stopped(obs_per_cell = c(400, 4.5), partition = "kmeans")
  expect_s3_class(r$error, "spatialkit_work_limit")
  expect_length(r$warnings, 0L)
  # The guard on the density weights' search, and the one on each row's nearest centre.
  r <- stopped(obs_per_cell = 4.5, partition = "density", fit_max = 1e10, pairs_max = 10, off = "FNN")
  expect_identical(r$error$what, "density weights")
  expect_identical(r$variogram, 0L)
  expect_length(r$warnings, 0L)
  r <- stopped(obs_per_cell = 4.5, partition = "kmeans", fit_max = 1e10, pairs_max = 10, off = "FNN")
  expect_identical(r$error$what, "nearest centre")
  expect_identical(r$error$work, 300 * 67)
  expect_identical(r$variogram, 0L)
  expect_length(r$warnings, 0L)
  # The order of the three: the weights' search, the fits, each row's nearest centre.
  r <- stopped(obs_per_cell = 4.5, partition = "density", fit_max = 1, pairs_max = 10, off = "FNN")
  expect_identical(r$error$what, "density weights")
  r <- stopped(obs_per_cell = 4.5, partition = "kmeans", fit_max = 1, pairs_max = 10, off = "FNN")
  expect_identical(r$error$what, "partition fit")
  # Under the limits the variogram is estimated, as on every ladder under sac = NULL.
  r <- stopped(obs_per_cell = 4.5, partition = "kmeans", fit_max = 1e10)
  expect_identical(conditionMessage(r$error), "a fit was reached")
  expect_identical(r$variogram, 1L)
})

test_that("without FNN the fit guard counts every pass of the density kind", {
  pts <- .cl_half()
  guard <- function(limit, ...) {
    withr::local_options(spatialkit.fit_work_max = limit)
    testthat::local_mocked_bindings(.partition_centres = function(...) stop("a fit was reached"))
    .cl_without("FNN")
    .cl_error(.cl_profile(pts, obs_per_cell = 4.5, ...))
  }
  for (init in c("kmeans++", "random")) {
    # Up to 100 passes a start, each of which compares every point with every centre.
    work <- 1000 * 111 * 3 * 100
    e <- guard(work - 1, partition = "density", init = init, nstart = NULL)
    expect_s3_class(e, "spatialkit_work_limit")
    expect_identical(.cl_fields(e),
                     list(work = work, limit = work - 1, what = "partition fit", kind = "density",
                          init = init, nstart = 3L, fnn = FALSE), info = init)
    expect_match(conditionMessage(e),
                 "fitting 111 cells to 1000 rows with 3 starts of up to 100 passes each is 3.3e+07",
                 fixed = TRUE)
    expect_identical(conditionMessage(guard(work, partition = "density", init = init,
                                            nstart = NULL)),
                     "a fit was reached", info = init)
    # `iter_max` is the number of passes counted.
    e <- guard(1, partition = "density", init = init, nstart = 2, iter_max = 7)
    expect_identical(e$work, 1000 * 111 * 2 * 7, info = init)
    expect_match(conditionMessage(e), "with 2 starts of up to 7 passes each", fixed = TRUE)
  }
  # The k-means kind has no kd-tree to lose: its count is the same.
  e <- guard(1, partition = "kmeans", init = "kmeans++", nstart = NULL)
  expect_identical(e$work, 1000 * 111 * 25)
  expect_false(e$fnn)
})

test_that("the fit guard is not applied to a ladder the subsample bounds", {
  skip_if_not_installed("gstat")
  pts <- .cl_layer(300, 319)
  withr::local_options(spatialkit.fit_work_max = 0)
  for (kind in c("kmeans", "density")) {
    expect_s3_class(suppressWarnings(.cl_profile(pts, partition = kind, estimate = NULL,
                                                 obs_per_cell = NULL, sac = NULL, n_levels = 3)),
                    "resolution_profile")
    expect_s3_class(suppressWarnings(.cl_profile(pts, partition = kind, obs_per_cell = NULL,
                                                 sac = NULL, levels = c(5L, 20L))),
                    "resolution_profile")
  }
  # On a ladder of sizes the same limit stops the call.
  expect_s3_class(.cl_error(.cl_profile(pts, partition = "kmeans", obs_per_cell = 4.5)),
                  "spatialkit_work_limit")
})

test_that("the limit on a fit is 1e10 comparisons unless the option says otherwise", {
  expect_identical(.FIT_WORK_MAX, 1e10)
  withr::local_options(spatialkit.fit_work_max = NULL)
  expect_identical(.fit_work_max(), 1e10)
  part <- list(kind = "kmeans", init = "kmeans++")
  # 100,000 points into 22,222 cells from 5 starts is past it; from 4 it is not.
  expect_s3_class(tryCatch(.fit_guard(100000L, 22222L, 5L, part, 50L, TRUE), error = identity),
                  "spatialkit_work_limit")
  expect_null(.fit_guard(100000L, 22222L, 4L, part, 50L, TRUE))
  # Counts past the largest integer are counted as doubles.
  e <- tryCatch(.fit_guard(3e6L, 7e5L, 25L, part, 50L, TRUE), error = identity)
  expect_identical(e$work, 3e6 * 7e5 * 25)
  withr::local_options(spatialkit.fit_work_max = 5e11)
  expect_identical(.fit_work_max(), 5e11)
  expect_null(.fit_guard(100000L, 22222L, 25L, part, 50L, TRUE))
  for (bad in list("many", c(1, 2), NA_real_, -1)) {
    withr::local_options(spatialkit.fit_work_max = bad)
    expect_error(.fit_work_max(), "options(spatialkit.fit_work_max = ) must be a single number",
                 fixed = TRUE)
  }
})

test_that("one density start from k-means++ takes more than twice as long on twice the rows", {
  skip_on_cran()
  skip_if_not_installed("FNN")
  one <- function(n, seed) {
    xy <- local({
      .cl_seed(seed)
      cbind(stats::runif(n, 0, 1e5), stats::runif(n, 0, 1e5))
    })
    .cl_seed(123)
    system.time(.partition_centres(xy, as.integer(round(n / 4.5)), partition = "density",
                                   init = "kmeans++", nstart = 1L, iter.max = 100L,
                                   use_fnn = TRUE))[["elapsed"]]
  }
  t20 <- one(20000L, 201)
  t40 <- one(40000L, 301)
  expect_gt(t40, 2 * t20)
})


# ---- The guard on the density weights' search --------------------------------

test_that("without FNN a ladder of sizes bounds the search behind the density kind's weights", {
  # 60 points at 50 distinct locations.
  pts <- local({
    d <- sf::st_drop_geometry(.cl_layer(60, 320))
    g <- sf::st_geometry(.cl_layer(60, 320))
    sf::st_sf(d, geometry = g[c(1:50, 1:10)])
  })
  expect_identical(nrow(unique(sf::st_coordinates(pts))), 50L)
  call <- function(...) {
    testthat::local_mocked_bindings(.partition_centres = function(...) stop("a fit was reached"))
    .cl_error(suppressWarnings(.cl_profile(pts, obs_per_cell = 4.5, ...)))
  }
  withr::local_options(spatialkit.exact_pairs_max = 100)
  e <- local({ .cl_without("FNN"); call() })
  expect_identical(class(e), c("spatialkit_work_limit", "error", "condition"))
  # The distinct locations squared, not the points squared.
  expect_identical(.cl_fields(e), list(work = 2500, limit = 100, what = "density weights"))
  expect_identical(conditionMessage(e),
                   paste0("resolution_profile(): the neighbour search behind the weights of ",
                          "partition = \"density\", over 50 distinct locations, would compare ",
                          "2.5e+03 pairs without the FNN package, above the limit of 1e+02 ",
                          "(options(spatialkit.exact_pairs_max = )). Install FNN: the result is the ",
                          "same."))
  expect_null(conditionCall(e))
  # At the limit the search runs.
  local({
    withr::local_options(spatialkit.exact_pairs_max = 2500)
    .cl_without("FNN")
    expect_identical(conditionMessage(call()), "a fit was reached")
  })
  # The search belongs to the fit and not to the estimate: it is bounded
  # without one, and without a response.
  for (extra in list(list(estimate = NULL), list(response = NULL, estimate = NULL))) {
    e <- local({ .cl_without("FNN"); do.call(call, extra) })
    expect_identical(.cl_fields(e), list(work = 2500, limit = 100, what = "density weights"))
  }
  # It is bounded before the fits are: a call both guards would stop names the search.
  e <- local({
    .cl_without("FNN")
    withr::local_options(spatialkit.fit_work_max = 1)
    call()
  })
  expect_identical(e$what, "density weights")
  # With no weights there is no such search: another search stops the call.
  e <- local({ .cl_without("FNN"); call(density_power = 0) })
  expect_s3_class(e, "spatialkit_work_limit")
  expect_false(identical(e$what, "density weights"))
  # The k-means kind has no weights either.
  e <- local({ .cl_without("FNN"); call(partition = "kmeans") })
  expect_false(identical(e$what, "density weights"))
  # With FNN the weights come from a kd-tree, at any limit.
  skip_if_not_installed("FNN")
  withr::local_options(spatialkit.exact_pairs_max = 0)
  expect_identical(conditionMessage(call()), "a fit was reached")
})

test_that("the search behind the density weights is not bounded on a ladder the subsample bounds", {
  skip_if_not_installed("gstat")
  pts <- .cl_layer(120, 321)
  withr::local_options(spatialkit.exact_pairs_max = 0)
  .cl_without("FNN")
  # No estimate, so no search of the estimate's is made either.  The ladder
  # the function builds, and one given as `levels`.
  for (ladder in list(list(n_levels = 3), list(levels = c(4L, 10L)))) {
    p <- suppressWarnings(do.call(.cl_profile, c(list(pts, estimate = NULL, obs_per_cell = NULL,
                                                      sac = NULL), ladder)))
    expect_s3_class(p, "resolution_profile")
    expect_true(all(is.finite(p$wss)))
    expect_null(attr(p, "bounds")$obs_per_cell)
  }
})


# ---- The print ---------------------------------------------------------------

test_that("a ladder of sizes prints its sizes, the rows they divide and no ceiling", {
  pts <- .cl_600()
  one <- .cl_profile(pts, obs_per_cell = 4.5)
  out <- utils::capture.output(print(one))
  expect_identical(out[1L], "Resolution profile: 1 levels on 600 points")
  expect_identical(out[2L],
                   paste0("  ladder      : 133 cells (obs_per_cell = 4.5 on the 600 rows with a ",
                          "response; centres fitted to every row; floor 2 (no range))"))
  expect_identical(out[3L],
                   paste0("  variogram   : none fitted: a size ladder fits one only when `sac` is ",
                          "given (cp and reliability are NA)"))
  expect_identical(out[4L],
                   paste0("  scored on   : response; 3 density-weighted k-means++ starts (power 1, 8 ",
                          "neighbours) per level; WSS rises at 0 step(s)"))
  expect_match(out[5L], "^  one cell    : rss ")
  expect_identical(out[6L],
                   sprintf(paste0("  estimate    : pooled toward the 8 nearest cells; loo on 600 rows ",
                                  "of the layer; one cell: loo %s"),
                           format(signif(attr(one, "estimate")$loo_one_cell, 3))))
  # The table shows the new columns like the others.
  for (col in c("scored_per_cell", "loo", "loo_se", "share_alone", "own_weight"))
    expect_true(any(grepl(paste0("\\b", col, "\\b"), out)), info = col)

  several <- .cl_profile(pts, estimate = "mean", obs_per_cell = "auto")
  out <- utils::capture.output(print(several))
  expect_identical(out[2L],
                   paste0("  ladder      : 33 to 300 cells (obs_per_cell = 18, 12, 9, 6, 4.5, 3, 2 ",
                          "on the 600 rows with a response; centres fitted to every row; floor 2 ",
                          "(no range))"))
  expect_match(out[6L], "^  estimate    : plain means; loo on 600 rows of the layer; one cell: loo ")

  # A size of "auto" that gave fewer than two cells was left out of the
  # ladder, and is left out of its line: 18 rows a cell on 20 rows.
  few <- .cl_profile(.cl_layer(20, 306), estimate = "mean", obs_per_cell = "auto")
  expect_identical(attr(few, "bounds")$obs_per_cell, c(18, 12, 9, 6, 4.5, 3, 2))
  out <- utils::capture.output(print(few))
  expect_identical(out[2L],
                   paste0("  ladder      : 2 to 10 cells (obs_per_cell = 12, 9, 6, 4.5, 3, 2 on the ",
                          "20 rows with a response; centres fitted to every row; floor 2 (no ",
                          "range))"))

  # Rows without a response, a range, and no response at all.
  part <- pts
  part$z[1:150] <- NA
  out <- utils::capture.output(print(suppressWarnings(.cl_profile(part, obs_per_cell = 9, sac = 20000))))
  expect_identical(out[1L], "Resolution profile: 1 levels on 600 points (450 with a response)")
  expect_match(out[2L], paste0("^  ladder      : 50 cells \\(obs_per_cell = 9 on the 450 rows with a ",
                               "response; centres fitted to every row; floor [0-9]+ from range 20000, ",
                               "not applied\\)$"))
  out <- utils::capture.output(print(.cl_profile(pts, response = NULL, estimate = NULL,
                                                 obs_per_cell = 9)))
  expect_identical(out[2L],
                   paste0("  ladder      : 67 cells (obs_per_cell = 9 on the 600 points; centres ",
                          "fitted to every row; floor 2 (no range))"))
  expect_false(any(grepl("^  estimate", out)))

  # Not supported.
  out <- utils::capture.output(print(suppressWarnings(.cl_profile(.cl_layer(6, 305),
                                                                   obs_per_cell = 4.5))))
  expect_match(out[2L], "  -- fewer than two cells of the size: not supported$")
  out <- utils::capture.output(print(suppressWarnings(.cl_profile(.cl_3000(), obs_per_cell = 4.5,
                                                                   sac = 2000))))
  expect_match(out[2L], "not applied\\)  -- floor above the finest level: not supported$")
})

test_that("a profile without an estimate and without a size prints no new line", {
  skip_if_not_installed("gstat")
  pts <- .cl_600()
  out <- utils::capture.output(print(suppressWarnings(.cl_profile(pts, estimate = NULL,
                                                                   obs_per_cell = NULL, sac = NULL,
                                                                   n_levels = 4))))
  expect_false(any(grepl("^  estimate|obs_per_cell|none fitted|scored_per_cell", out)))
  expect_match(out[2L], "^  ladder      : 2 to 66 cells \\(floor ")
})
