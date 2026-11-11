# tests/testthat/test-fix-split-ceiling-text.R
# ---------------------------------------------------------------------------
# The ceiling texts of resolution_profile() say the bound the code applies.
#
# * Under select_on = "split" the response criteria read the selection half,
#   and on a subsample the support ceiling counts that half's rows of it.  The
#   warning of class spatialkit_response_ceiling called them the subsample's
#   rows that carry a response ("0 of the 500 points of the subsample carry a
#   response" of a subsample in which about half did, all in the estimation
#   half), set the layer's rows of both halves over them, and with none of
#   them printed "the cells of the layer hold about NA times as many", with or
#   without the split.
# * The warning for a range floor above the ceiling said "the profile runs
#   from 2 to" the ceiling over the levels given in `levels`, and named
#   min_cell_n, or the size, where the distinct locations set the ceiling.
# * Where min_cell_n and the distinct locations set the same ceiling, the
#   edge was "the support ceiling (n / min_cell_n)" and its hint offered
#   `levels` past it and a lower min_cell_n, neither of which moves it.
#
# The rules, the ceilings, the levels and the conditions' fields do not
# change; only the texts do.
# ---------------------------------------------------------------------------

fsc_layer <- function(n = 3000, seed = 21) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(n, 0, 5000), y = stats::runif(n, 0, 5000))
  d$z <- sin(d$x / 700) + cos(d$y / 900) + stats::rnorm(n, sd = 0.3)
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32615)
}
# `n_rep` rows at each of `n_st` stations on a 900 m square, a response of
# the station's mean plus noise.
fsc_stations <- function(n_st, n_rep, seed = 153) {
  set.seed(seed)
  st <- data.frame(x = stats::runif(n_st, 0, 900), y = stats::runif(n_st, 0, 900),
                   mu = 4 * seq_len(n_st))
  g <- st[rep(seq_len(n_st), each = n_rep), ]
  g$z <- g$mu + stats::rnorm(nrow(g), sd = 0.3)
  sf::st_as_sf(g, coords = c("x", "y"), crs = 32632)
}
fsc_distinct <- function(n, seed = 154) {
  set.seed(seed)
  sf::st_as_sf(data.frame(x = stats::runif(n, 0, 900), y = stats::runif(n, 0, 900),
                          z = stats::rnorm(n)), coords = c("x", "y"), crs = 32632)
}
# Every warning an expression raises, none of them left to reach the reporter.
fsc_catch <- function(expr) {
  ws <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    ws[[length(ws) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = ws)
}
fsc_of <- function(r, class) Filter(function(w) inherits(w, class), r$warnings)
fsc_ceiling <- function(r) fsc_of(r, "spatialkit_response_ceiling")
# A variogram range given as a number, so no variogram is fitted.
fsc_profile <- function(pts, sample_n = 500, n_levels = 4, ...)
  resolution_profile(pts, "z", estimate = NULL, obs_per_cell = NULL, sac = 2000,
                     n_levels = n_levels, sample_n = sample_n, nstart = 2, seed = 123, ...)
fsc_split <- function(pts)
  attr(suppressWarnings(resolution_profile(pts, estimate = NULL, obs_per_cell = NULL,
                                           levels = 2, nstart = 1, sample_n = 500,
                                           seed = 123, select_on = "split")), "split")

# ---- BUG-39: the response ceiling under the split -----------------------------

test_that("with no response in the selection half the warning names the half and has no NA", {
  pts <- fsc_layer()
  sp <- fsc_split(pts)
  pts$z[sp$selection] <- NA
  n_resp <- sum(is.finite(pts$z))
  r <- fsc_catch(fsc_profile(pts, select_on = "split"))
  w <- fsc_ceiling(r)
  expect_length(w, 1L)
  w <- w[[1L]]
  msg <- conditionMessage(w)
  expect_true(startsWith(msg, sprintf(paste0(
    "resolution_profile(): 0 of the 500 points of the subsample carry a response in the ",
    "selection half, the half the criteria read under select_on = \"split\" (%d of the ",
    "layer's 3000 carry a response in either half), so the ceiling is 2 cells, two scored ",
    "subsample rows per cell, not the 250 that half the subsample would give."), n_resp)))
  expect_false(grepl("\\bNA\\b", msg))
  expect_false(grepl("cell_resp_mean", msg, fixed = TRUE))
  expect_true(endsWith(msg, "and a larger `sample_n` raises the ceiling towards that."))
  # The ceiling and the condition's fields are the rule's, as they were.
  b <- attr(r$value, "bounds")
  expect_identical(b$ceiling, 2L)
  expect_identical(b$ceiling_from, "sample_n")
  expect_identical(w$n_scored, 0L)
  expect_identical(w$n_response, n_resp)
  expect_identical(w$n_sample, 500L)
  expect_identical(w$ceiling, 2L)
  expect_identical(r$value$cells_resp, rep(0L, nrow(r$value)))
  expect_true(all(is.na(r$value$cell_resp_mean)))
})

test_that("with a response in both halves the count is the selection half's, and so is the factor", {
  pts <- fsc_layer()
  set.seed(4); pts$z[sample(3000, 1800)] <- NA
  sp <- fsc_split(pts)
  n_sel <- sum(is.finite(pts$z[sp$selection]))
  r <- fsc_catch(fsc_profile(pts, select_on = "split"))
  w <- fsc_ceiling(r)
  expect_length(w, 1L)
  msg <- conditionMessage(w[[1L]])
  b <- attr(r$value, "bounds")
  expect_identical(b$ceiling_from, "sample_n")
  expect_identical(attr(r$value, "split")$n_response, n_sel)
  # The same subsample without the split: the rows of both halves that
  # carry a response, which the split's count is a part of.
  b_all <- attr(suppressWarnings(fsc_profile(pts)), "bounds")
  expect_lt(b$n_scored, b_all$n_scored)
  expect_true(startsWith(msg, sprintf(paste0(
    "resolution_profile(): %d of the 500 points of the subsample carry a response in the ",
    "selection half, the half the criteria read under select_on = \"split\" (1200 of the ",
    "layer's 3000 carry a response in either half), so the ceiling is %d cells"),
    b$n_scored, b$ceiling)))
  # `cell_resp_mean` counts the selection half's rows of the subsample, so
  # the layer's cells hold more of the selection half's rows, not of both.
  expect_true(endsWith(msg, sprintf(paste0(
    "`cell_resp_mean` counts the scored rows of the subsample per occupied cell; the cells ",
    "of the layer hold about %.3g times as many."), n_sel / b$n_scored)))
  expect_false(grepl(sprintf("about %.3g times", 1200 / b$n_scored), msg, fixed = TRUE))
  expect_equal(r$value$cell_resp_mean, b$n_scored / r$value$cells_resp)
})

test_that("without the split, and without a subsample, the texts are as they were", {
  pts <- fsc_layer()
  set.seed(4); pts$z[sample(3000, 1800)] <- NA
  r <- fsc_catch(fsc_profile(pts))
  msg <- conditionMessage(fsc_ceiling(r)[[1L]])
  b <- attr(r$value, "bounds")
  expect_identical(msg, sprintf(paste0(
    "resolution_profile(): %d of the 500 points of the subsample carry a response (1200 of ",
    "the layer's 3000), so the ceiling is %d cells, two scored subsample rows per cell, not ",
    "the 250 that half the subsample would give. `sample_n` sets this ceiling, not ",
    "`min_cell_n`: the layer's 1200 rows with a response support 133 cells at min_cell_n = ",
    "9, and a larger `sample_n` raises the ceiling towards that. `cell_resp_mean` counts the ",
    "scored rows of the subsample per occupied cell; the cells of the layer hold about %.3g ",
    "times as many."), b$n_scored, b$ceiling, 1200 / b$n_scored))

  # The split on a layer below `sample_n`: the layer's counts, no subsample.
  small <- fsc_layer(400)
  sp <- fsc_split(small)
  small$z[sp$selection] <- NA
  r <- fsc_catch(fsc_profile(small, select_on = "split", n_levels = 5))
  msg <- conditionMessage(fsc_ceiling(r)[[1L]])
  n_resp <- sum(is.finite(small$z))
  expect_true(startsWith(msg, sprintf(paste0(
    "resolution_profile(): %d of the 400 points carry a response, so the ceiling is %d ",
    "cells, not the 44 that every point at min_cell_n = 9 would give"),
    n_resp, attr(r$value, "bounds")$ceiling)))
  expect_false(grepl("selection half", msg, fixed = TRUE))
})

test_that("without the split, a subsample with no scored row leaves the factor out", {
  pts <- fsc_layer()
  pts$z[-c(11, 900, 1700, 2600)] <- NA
  # The first seed whose subsample of 500 draws none of the 4 rows with a
  # response, found by the profile itself.
  n_scored_at <- function(s)
    attr(suppressWarnings(resolution_profile(pts, "z", estimate = NULL, obs_per_cell = NULL,
                                             sac = 2000, levels = 2, min_cell_n = 2,
                                             nstart = 1, sample_n = 500, seed = s)),
         "bounds")$n_scored
  s0 <- Find(function(s) identical(n_scored_at(s), 0L), 1:40)
  expect_false(is.null(s0))
  r <- fsc_catch(resolution_profile(pts, "z", estimate = NULL, obs_per_cell = NULL, sac = 2000,
                                    n_levels = 3, min_cell_n = 2, sample_n = 500, nstart = 2,
                                    seed = s0))
  b <- attr(r$value, "bounds")
  # The case the test is about: no row of the subsample is scored.
  expect_identical(b$n_scored, 0L)
  msg <- conditionMessage(fsc_ceiling(r)[[1L]])
  expect_true(startsWith(msg, paste0(
    "resolution_profile(): 0 of the 500 points of the subsample carry a response (4 of the ",
    "layer's 3000), so the ceiling is 2 cells")))
  expect_true(endsWith(msg, "and a larger `sample_n` raises the ceiling towards that."))
  expect_false(grepl("\\bNA\\b", msg))
  expect_false(grepl("cell_resp_mean", msg, fixed = TRUE))
})

test_that("under the split, the factor where min_cell_n sets the ceiling is the selection half's", {
  # 4000 points, a response on 400 of them, a subsample of 3500: the rows
  # with a response set the ceiling at min_cell_n = 30, not `sample_n`.
  pts <- fsc_layer(4000, seed = 5)
  set.seed(5); pts$z[sample(4000, 3600)] <- NA
  r <- fsc_catch(resolution_profile(pts, "z", estimate = NULL, obs_per_cell = NULL, sac = 2000,
                                    n_levels = 4, sample_n = 3500, min_cell_n = 30, nstart = 2,
                                    seed = 1, select_on = "split"))
  w <- fsc_ceiling(r)
  expect_length(w, 1L)
  msg <- conditionMessage(w[[1L]])
  b <- attr(r$value, "bounds")
  n_sel <- attr(r$value, "split")$n_response
  expect_identical(b$ceiling_from, "min_cell_n (rows with a response)")
  expect_lt(b$n_scored, 400L)
  expect_lt(n_sel, 400L)
  # `n_scored` and `split$n_response` both count the selection half's rows,
  # of the subsample and of the layer; the layer's 400 count both halves.
  expect_true(grepl(sprintf(paste0("(rows of the subsample; the cells of the layer hold about ",
                                   "%.3g times as many)"), n_sel / b$n_scored),
                    msg, fixed = TRUE))
  expect_false(grepl(sprintf("about %.3g times", 400 / b$n_scored), msg, fixed = TRUE))
  # The ceiling counts the rows with a response in either half here, and
  # the text says so.
  expect_true(startsWith(msg, "resolution_profile(): 400 of the 4000 points carry a response, so the ceiling is"))
})

# ---- BUG-153: the range floor above the ceiling --------------------------------

fsc_floor_msg <- function(r) {
  w <- fsc_of(r, "spatialkit_not_supported")
  expect_length(w, 1L)
  expect_identical(w[[1L]]$reason, "floor above ceiling")
  conditionMessage(w[[1L]])
}

test_that("levels given are said to be scored as given, not a ladder from 2", {
  u <- fsc_distinct(60)
  r <- fsc_catch(resolution_profile(u, "z", sac = 10, levels = c(4, 5), estimate = NULL,
                                    obs_per_cell = NULL, nstart = 2))
  msg <- fsc_floor_msg(r)
  expect_identical(r$value$levels, c(4L, 5L))
  expect_true(endsWith(msg, paste0(
    "but 60 points at min_cell_n = 9 support at most 6. The data cannot support a ",
    "tessellation that respects their own correlation structure; the 2 levels asked for in ",
    "`levels` (4 to 5 cells) are scored as given so the cost of each level is still visible.")))
  expect_false(grepl("runs from 2", msg, fixed = TRUE))
  # Without `levels` the ladder the profile builds runs from 2, as it says.
  r <- fsc_catch(resolution_profile(u, "z", sac = 10, n_levels = 5, estimate = NULL,
                                    obs_per_cell = NULL, nstart = 2))
  msg <- fsc_floor_msg(r)
  expect_identical(attr(r$value, "bounds")$ceiling_from, "min_cell_n")
  expect_identical(min(r$value$levels), 2L)
  expect_identical(max(r$value$levels), 6L)
  expect_true(endsWith(msg, paste0(
    "but 60 points at min_cell_n = 9 support at most 6. The data cannot support a ",
    "tessellation that respects their own correlation structure; the profile runs from 2 ",
    "to 6 so the cost of each level is still visible.")))
})

test_that("the bound named is the distinct locations where they set the ceiling", {
  # 60 rows on 5 stations: min_cell_n = 9 allows 6 cells, the locations 5.
  p5 <- fsc_stations(5, 12)
  r <- fsc_catch(resolution_profile(p5, "z", sac = 10, n_levels = 5, estimate = NULL,
                                    obs_per_cell = NULL, nstart = 2))
  msg <- fsc_floor_msg(r)
  b <- attr(r$value, "bounds")
  expect_identical(b$ceiling, 5L)
  expect_identical(b$ceiling_from, "distinct locations")
  expect_true(grepl(paste0(
    "but the 60 points sit on 5 distinct locations, which allow at most 5 cells. The data ",
    "cannot support a tessellation that respects their own correlation structure; the ",
    "profile runs from 2 to 5 so"), msg, fixed = TRUE))
  expect_false(grepl("min_cell_n", msg, fixed = TRUE))
  # No location repeated and min_cell_n = 1: one short of the points binds.
  u <- fsc_distinct(30)
  r <- fsc_catch(resolution_profile(u, "z", sac = 10, n_levels = 4, min_cell_n = 1,
                                    estimate = NULL, obs_per_cell = NULL, nstart = 2))
  msg <- fsc_floor_msg(r)
  b <- attr(r$value, "bounds")
  expect_identical(b$ceiling, 29L)
  expect_identical(b$ceiling_from, "distinct locations")
  expect_true(grepl("but k-means can fit at most 29 cells to the 30 points, one short of the points.",
                    msg, fixed = TRUE))
})

test_that("a size held to the distinct locations is not called the cells that size gives", {
  p5 <- fsc_stations(5, 12)
  r <- fsc_catch(resolution_profile(p5, "z", sac = 10, obs_per_cell = 3, estimate = "mean",
                                    nstart = 2))
  msg <- fsc_floor_msg(r)
  expect_identical(r$value$levels, 5L)
  expect_identical(attr(r$value, "bounds")$ceiling_from, "distinct locations")
  expect_true(grepl(paste0(
    "but obs_per_cell = 3 asks for 20 cells of the 60 points, and k-means can place at most ",
    "5 on the layer's 5 distinct locations. The data cannot support a tessellation that ",
    "respects their own correlation structure; the profile scores the 5 cells the distinct ",
    "locations allow so the cost of each level is still visible."), msg, fixed = TRUE))
  expect_false(grepl("that size gives", msg, fixed = TRUE))
  # Two sizes, the finest held: the coarser one's count is its own.
  r <- fsc_catch(resolution_profile(p5, "z", sac = 10, obs_per_cell = c(20, 3),
                                    estimate = "mean", nstart = 2))
  msg <- fsc_floor_msg(r)
  expect_identical(r$value$levels, c(3L, 5L))
  expect_true(grepl(paste0(
    "the profile scores the 3 to 5 cells the sizes give, the finest held to the 5 the ",
    "distinct locations allow so the cost"), msg, fixed = TRUE))
})

test_that("where min_cell_n and the locations set the same ceiling, the floor warning names min_cell_n", {
  # 45 rows on 5 stations: min_cell_n = 9 allows 5 cells, and so do the
  # locations.  `ceiling_from` names min_cell_n at the tie, and the warning
  # takes the same test.
  p45 <- fsc_stations(5, 9)
  r <- fsc_catch(resolution_profile(p45, "z", sac = 10, n_levels = 5, estimate = NULL,
                                    obs_per_cell = NULL, nstart = 2))
  msg <- fsc_floor_msg(r)
  b <- attr(r$value, "bounds")
  expect_identical(b$ceiling, 5L)
  expect_identical(b$ceiling_from, "min_cell_n")
  expect_true(grepl(paste0(
    "but 45 points at min_cell_n = 9 support at most 5. The data cannot support a ",
    "tessellation that respects their own correlation structure; the profile runs from 2 ",
    "to 5 so"), msg, fixed = TRUE))
  expect_false(grepl("distinct locations", msg, fixed = TRUE))
})

test_that("a size held on points with no location repeated names the most k-means can fit", {
  # 40 distinct points at obs_per_cell = 1: the size asks for 40 cells, and
  # k-means fits at most 39, one short of the points.  The 40 locations
  # would allow 40, so the bound is not theirs.
  u <- fsc_distinct(40)
  r <- fsc_catch(resolution_profile(u, "z", sac = 10, obs_per_cell = 1, estimate = "mean",
                                    nstart = 2))
  msg <- fsc_floor_msg(r)
  b <- attr(r$value, "bounds")
  expect_identical(r$value$levels, 39L)
  expect_identical(b$ceiling_from, "distinct locations")
  expect_identical(b$n_distinct, 40L)
  expect_true(grepl(paste0(
    "but obs_per_cell = 1 asks for 40 cells of the 40 points, and k-means can place at most ",
    "39 on 40 points. The data cannot support a tessellation that respects their own ",
    "correlation structure; the profile scores the 39 cells k-means can fit to the points ",
    "so the cost of each level is still visible."), msg, fixed = TRUE))
  expect_false(grepl("distinct locations allow", msg, fixed = TRUE))
  expect_false(grepl("that size gives", msg, fixed = TRUE))
})

test_that("a size that asks for exactly the cells k-means can place is the size's count", {
  # 60 rows on 5 stations at obs_per_cell = 12: the size asks for 5 cells,
  # as many as the locations allow, so the size sets the ceiling and the
  # warning says the size's count, as `ceiling_from` does.
  p5 <- fsc_stations(5, 12)
  r <- fsc_catch(resolution_profile(p5, "z", sac = 10, obs_per_cell = 12, estimate = "mean",
                                    nstart = 2))
  msg <- fsc_floor_msg(r)
  expect_identical(r$value$levels, 5L)
  expect_identical(attr(r$value, "bounds")$ceiling_from, "obs_per_cell")
  expect_true(grepl(paste0(
    "but 60 points at obs_per_cell = 12 support at most 5. The data cannot support a ",
    "tessellation that respects their own correlation structure; the profile scores the 5 ",
    "cells that size gives so the cost of each level is still visible."), msg, fixed = TRUE))
  expect_false(grepl("asks for", msg, fixed = TRUE))
})

test_that("a subsample below the floor says levels given are scored as given", {
  pts <- fsc_layer(3000)
  # The layer supports 333 cells at min_cell_n = 9, the subsample of 100
  # fits 50, and a range of 400 asks for about 150.
  lines <- capture_spatialkit_log(
    pr <- suppressWarnings(resolution_profile(pts, "z", sac = 400, levels = c(10, 20),
                                              estimate = NULL, obs_per_cell = NULL,
                                              sample_n = 100, nstart = 1)),
    level = logger::WARN)
  b <- attr(pr, "bounds")
  expect_true(b$supported)
  expect_gt(b$floor, b$ceiling)
  line <- grep("raise `sample_n`", lines, value = TRUE, fixed = TRUE)
  expect_length(line, 1L)
  expect_match(line, paste0("raise `sample_n`. The 2 levels asked for in `levels` ",
                            "(10 to 20 cells) are scored as given."), fixed = TRUE)
  expect_false(grepl("runs from 2", line, fixed = TRUE))
  expect_identical(pr$levels, c(10L, 20L))
})

# ---- BUG-161: the edge where min_cell_n and the locations tie ------------------

test_that("an optimum at a ceiling the locations set as well names the locations, with no remedy", {
  # 45 rows on 5 stations: 45 / 9 = 5 = the distinct locations.
  set.seed(161)
  st <- data.frame(x = c(0, 2000, 0, 2000, 1000), y = c(0, 0, 2000, 2000, 1000),
                   mu = c(0, 4, 8, 12, 16))
  h <- st[rep(1:5, each = 9), ]
  h$z <- h$mu + stats::rnorm(45, sd = 0.3)
  p45 <- sf::st_as_sf(h, coords = c("x", "y"), crs = 32632)
  pr <- suppressWarnings(resolution_profile(p45, "z", estimate = NULL, obs_per_cell = NULL,
                                            n_levels = 5, sac = 1500, nstart = 2))
  b <- attr(pr, "bounds")
  # The record is as it was: min_cell_n names the tie.
  expect_identical(b$ceiling, 5L)
  expect_identical(b$ceiling_from, "min_cell_n")
  expect_identical(b$n_distinct, 5L)
  sel <- suppressWarnings(select_resolution(pr, criterion = "cp"))
  expect_identical(sel$best, 5L)
  label <- "the support ceiling (as many cells as the distinct locations allow)"
  expect_identical(sel$edge, label)
  out <- paste(utils::capture.output(print(sel)), collapse = " ")
  out <- gsub("\\s+", " ", out)
  expect_match(out, "the optimum is the support ceiling (as many cells as the distinct locations allow); the bound is choosing, not the criterion.",
               fixed = TRUE)
  expect_false(grepl("Pass `levels`", out, fixed = TRUE))
  expect_false(grepl("lower min_cell_n", out, fixed = TRUE))
  r <- fsc_catch(get_voronoi_seeds(method = "kmeans", n = sel, sample_points = p45))
  w <- fsc_of(r, "spatialkit_bound_chosen")
  expect_length(w, 1L)
  expect_identical(w[[1L]]$edge, label)
  expect_false(grepl("min_cell_n", conditionMessage(w[[1L]]), fixed = TRUE))
  # Neither remedy moves the ceiling: levels past it are dropped, and a
  # lower min_cell_n leaves it where the locations put it.
  pr1 <- suppressWarnings(resolution_profile(p45, "z", estimate = NULL, obs_per_cell = NULL,
                                             levels = c(2, 5, 7, 9), sac = 1500, nstart = 2))
  expect_identical(pr1$levels, c(2L, 5L))
  pr2 <- suppressWarnings(resolution_profile(p45, "z", estimate = NULL, obs_per_cell = NULL,
                                             n_levels = 5, min_cell_n = 3, sac = 1500,
                                             nstart = 2))
  expect_identical(attr(pr2, "bounds")$ceiling, 5L)
})

test_that("the edge label reads the tie from the bounds, and min_cell_n alone keeps its own", {
  lv <- c(2L, 5L)
  tie <- list(ceiling = 5L, ceiling_from = "min_cell_n", n_distinct = 5L, n = 45L,
              floor = 2L, supported = TRUE, range = NA_real_, range_floor = FALSE)
  loc <- "the support ceiling (as many cells as the distinct locations allow)"
  expect_identical(.ladder_edge(5L, lv, lv, tie), loc)
  expect_identical(.edge_hint(loc), "")
  # The rows with a response set the same ceiling as the locations.
  rsp <- utils::modifyList(tie, list(ceiling_from = "min_cell_n (rows with a response)"))
  expect_identical(.ladder_edge(5L, lv, lv, rsp), loc)
  # No location repeated, min_cell_n = 1: one short of the points.
  short <- utils::modifyList(tie, list(n_distinct = 6L, n = 6L))
  expect_identical(.ladder_edge(5L, lv, lv, short),
                   "the support ceiling (one short of the points, the most k-means can fit)")
  # min_cell_n alone binds: more locations than the ceiling.
  own <- utils::modifyList(tie, list(n_distinct = 8L))
  expect_identical(.ladder_edge(5L, lv, lv, own), "the support ceiling (n / min_cell_n)")
  expect_match(.edge_hint(.ladder_edge(5L, lv, lv, own)), "a lower min_cell_n raises the ceiling",
               fixed = TRUE)
  # A record without the locations keeps the min_cell_n label.
  old <- tie[setdiff(names(tie), "n_distinct")]
  expect_identical(.ladder_edge(5L, lv, lv, old), "the support ceiling (n / min_cell_n)")
})
