# tests/testthat/test-fix-aoa-rf-plot.R
# ---------------------------------------------------------------------------
# plot.aoa() and plot.feature_selection() against the results this group of
# fixes changed.
#
# (a) With a location index (a gwr_fit, a bayesian_fit, coords = "yes")
#     `n_outside` counts the rows outside on EITHER index.  plot.aoa() printed
#     it beside the predictor threshold: a GWR result with 10 of 100 rows
#     beyond that threshold and all 100 beyond the location one read "100 of
#     100 prediction locations outside (DI > 0.017)", then "90% of those
#     inside sit below 51% of the threshold".  The count beside "(DI > ...)"
#     has to be the predictor index's, `n_outside_predictors`.  And under
#     categorical = "indicator" a row on a level no training row has was
#     captioned "outside on a dropped predictor".
#
# (b) Under a fold-mean baseline select_features_forward() can drop the
#     first variable at the end, so `selected` is the walked path without
#     it.  plot.feature_selection() matched selected[s] to step s and so
#     drew sets that were never on the path.
#
# Both plot methods live in R/plotting-diagnostics.R and read the fields the
# two functions return: `n_outside_predictors`, `location`,
# `categorical_vars` and `params$categorical` of an aoa result, and
# `baseline$first` and `baseline$path` of a selection.  These tests fail
# when a plot method does not read them; they are not skipped for it,
# because a result drawn wrongly is the defect they are here to catch.
# What a plot says is compared as one value in each test, so such a plot
# method fails each test once and the difference shows every line it got
# wrong.
# ---------------------------------------------------------------------------

plt_pts <- function(m, x0, x1, seed) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(m, x0, x1), y = stats::runif(m, 0, 1000),
                  a = stats::rnorm(m))
  d$z <- d$a + stats::rnorm(m, 0, 0.3)
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
}
# Carries the class and the fields area_of_applicability() reads.
plt_fake_gwr <- function(train)
  structure(list(data_sf = train, predictor_vars = "a", response_var = "z",
                 info = list(), n = nrow(train)),
            class = c("gwr_fit", "spatial_fit"))
plt_lines <- function(s) strsplit(s, "\n", fixed = TRUE)[[1]]
# The subtitle, and the caption's line on location when it has one.  The
# share of the threshold in the subtitle's last line varies with the draw
# and is not what these tests are about, so it is written as N%.
plt_aoa_said <- function(p)
  gsub("below [0-9]+%", "below N%",
       c(plt_lines(p$labels$subtitle),
         grep("second index", plt_lines(p$labels$caption), value = TRUE)))
# The caption's sentence on the rows off the axis, on one line.
plt_off_axis <- function(p) {
  cap <- gsub("\n", " ", p$labels$caption, fixed = TRUE)
  regmatches(cap, regexpr("[0-9]+ prediction locations? outside on .* off the axis", cap))
}


test_that("plot.aoa() counts the predictor index beside the predictor threshold", {
  skip_if_not_installed("ggplot2")
  trn <- plt_pts(300, 0, 1000, seed = 3)
  far <- plt_pts(100, 1500, 2500, seed = 4)
  a <- area_of_applicability(far, model = plt_fake_gwr(trn))
  # All 100 are outside on location; far fewer on the predictor.
  n_pred <- sum(a$aoa$DI > a$threshold)
  expect_identical(a$n_outside, 100L)
  expect_identical(a$location$n_outside, 100L)
  expect_identical(a$n_outside_predictors, n_pred)
  expect_lt(n_pred, 50L)

  # Points across the edge of the training extent: some are outside on the
  # predictor, some on location, a few on both, so the three counts differ.
  edge <- plt_pts(100, 0, 1300, seed = 6)
  e <- area_of_applicability(edge, model = plt_fake_gwr(trn))
  e_pred <- sum(e$aoa$DI > e$threshold)
  e_loc  <- sum(e$aoa$DI_location > e$location$threshold)
  expect_identical(c(e$n_outside_predictors, e$location$n_outside), c(e_pred, e_loc))
  expect_gt(e_pred, 0L)
  expect_gt(e_loc, 0L)
  expect_gt(e$n_outside, max(e_pred, e_loc))
  expect_lt(e$n_outside, e_pred + e_loc)

  p  <- plot(a)
  pe <- plot(e)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_no_error(ggplot2::ggplot_build(pe))
  # Line 1 is the count that goes with "(DI > threshold)": the predictor
  # index's, beside the predictor threshold, not the rows outside on either
  # index.  Line 2 is location, and the total the AOA column gives.  Line 3
  # speaks of the predictor index, which is what it is computed from.  The
  # caption says what the curves are and gives the location threshold.
  should_say <- function(r, n_pred) c(
    sprintf("%d of 100 prediction locations outside on the predictors (DI > %.3g)",
            n_pred, r$threshold),
    sprintf("%d outside on location; %d outside the AOA (either)",
            r$location$n_outside, r$n_outside),
    "90% of those inside on the predictors sit below N% of that threshold",
    sprintf(paste0("The curves are the predictor index; location is a second ",
                   "index (DI_location, threshold %.3g)"), r$location$threshold))
  expect_identical(list(far = plt_aoa_said(p), edge = plt_aoa_said(pe)),
                   list(far = should_say(a, n_pred), edge = should_say(e, e_pred)))
  expect_no_error(ggplot2::ggplot_build(plot(a, type = "histogram")))

  # Without a location index nothing changes: one count, and it is n_outside.
  b <- area_of_applicability(far, train_sf = trn, predictor_vars = "a")
  pb <- plot(b)
  expect_match(pb$labels$subtitle,
               sprintf("^%d of 100 prediction locations outside \\(DI > ", b$n_outside))
  expect_no_match(pb$labels$subtitle, "outside on location")
  expect_no_match(pb$labels$caption, "second index")
})

test_that("plot.aoa() words DI = Inf as print.aoa() does", {
  skip_if_not_installed("ggplot2")
  set.seed(5)
  tr <- data.frame(a = stats::rnorm(80), g = rep(c("u", "f"), 40), k = 1)
  nd <- data.frame(a = stats::rnorm(30), g = c("zzz", "zzz", rep(c("u", "f"), 14)),
                   k = 1)
  a <- suppressWarnings(area_of_applicability(nd, train_sf = tr,
                                              predictor_vars = c("a", "g"),
                                              categorical = "indicator"))
  expect_identical(sum(is.infinite(a$aoa$DI)), 2L)
  # Every row on an unseen level: nothing to draw, and plot() stops.
  nd3 <- nd; nd3$g <- "zzz"
  d <- suppressWarnings(area_of_applicability(nd3, train_sf = tr,
                                              predictor_vars = c("a", "g"),
                                              categorical = "indicator"))
  expect_identical(sum(is.infinite(d$aoa$DI)), 30L)
  # The caption for the two rows, and the error for the thirty, each say
  # which kind of row DI = Inf is under this coding.
  expect_identical(
    c(caption = plt_off_axis(plot(a)),
      error = tryCatch({ plot(d); "plot() did not stop" }, error = conditionMessage)),
    c(caption = paste0("2 prediction locations outside on a level no training ",
                       "row has, or on a dropped predictor (DI = Inf) are off ",
                       "the axis"),
      error = paste0("plot.aoa(): no finite dissimilarity index to plot: 30 of ",
                     "30 prediction locations are outside on a level no ",
                     "training row has, or on a predictor dropped for having ",
                     "no training variance (DI = Inf).")))
  # A dropped numeric predictor, and nothing categorical: as before.
  nd2 <- nd; nd2$k[1] <- 5
  b <- suppressWarnings(area_of_applicability(nd2, train_sf = tr,
                                              predictor_vars = c("a", "k")))
  expect_identical(plt_off_axis(plot(b)),
                   "1 prediction location outside on a dropped predictor (DI = Inf) is off the axis")
})

test_that("plot.feature_selection() draws the walked path when the first variable was dropped", {
  skip_if_not_installed("ggplot2")
  # `proxy` is the best single predictor and of no use once a1 and a2 are in.
  set.seed(4001)
  n <- 240
  d <- data.frame(x = stats::runif(n, 0, 2000), y = stats::runif(n, 0, 2000),
                  a1 = stats::rnorm(n), a2 = stats::rnorm(n), n1 = stats::rnorm(n))
  d$proxy <- d$a1 + d$a2 + stats::rnorm(n)
  d$z     <- d$a1 + d$a2 + stats::rnorm(n, 0, 0.5)
  d <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  refusing <- function(tr, v) {
    if (!length(v)) stop("this backend needs at least one predictor", call. = FALSE)
    lm_spatial_fit(tr, "z", v)
  }
  sel <- suppressWarnings(suppressMessages(
    select_features_forward(d, "z", c("proxy", "a1", "a2", "n1"), refusing,
                            k = 4, seed = 1, quiet = TRUE)))
  expect_true(sel$baseline$first$dropped)
  path <- sel$baseline$path
  expect_identical(path[1], "proxy")
  expect_identical(sel$selected, path[-1])
  # The selected set is not a step of the walk, so its score is in no row of
  # the history: a plot that looks it up there marks another set.
  h <- sel$history
  expect_false(any(abs(h$score - sel$score) < 1e-12, na.rm = TRUE))

  p <- plot(sel)
  b <- ggplot2::ggplot_build(p)
  geoms <- vapply(p$layers, function(l) class(l$geom)[1L], character(1))
  lab <- b$data[[which(geoms == "GeomText")]]
  pts <- which(geoms == "GeomPoint")
  chosen <- b$data[[pts[length(pts)]]]
  # The labelled points are the sets the sweep walked: step s is path[s], at
  # the score the history gives that set, and the variable taken out says
  # so.  The red mark is the selected set: its size, and its own score.
  at <- lab[match(seq_along(path), lab$x), ]
  expect_equal(
    list(score = at$y, label = at$label, red_mark = c(chosen$x, chosen$y),
         caption_names_the_drop = grepl(
           "proxy was added first without a test and dropped at the end",
           p$labels$caption, fixed = TRUE)),
    list(score = vapply(seq_along(path), function(s)
           h$score[h$step == s & h$variable == path[s]], numeric(1)),
         label = c("proxy (dropped at the end)", path[-1]),
         red_mark = c(length(sel$selected), sel$score),
         caption_names_the_drop = TRUE))
  expect_match(p$labels$caption, sprintf("Selected: %s", paste(sel$selected, collapse = ", ")))

  # Nothing dropped: the path is the selection, as before.
  own <- suppressWarnings(suppressMessages(
    select_features_forward(d, "z", c("proxy", "a1", "a2", "n1"),
                            function(tr, v) lm_spatial_fit(tr, "z", v),
                            k = 4, seed = 1, quiet = TRUE)))
  po <- plot(own)
  expect_no_match(po$labels$caption, "dropped at the end")
  bo <- ggplot2::ggplot_build(po)
  lab_o <- bo$data[[which(vapply(po$layers, function(l) class(l$geom)[1L],
                                 character(1)) == "GeomText")]]
  expect_identical(lab_o$label[lab_o$x == 1], own$selected[1])
})
