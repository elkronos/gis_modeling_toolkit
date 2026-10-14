# tests/testthat/test-core-count-size.R
# ---------------------------------------------------------------------------
# select_resolution(criterion = "loo" | "size" | "auto"): the three rules,
# the record of how the count was set, what the present criteria say on a
# ladder of sizes, a profile passed straight as a count, and what summary()
# and plot() show for a ladder of one size.
#
# Every profile names `partition`, `nstart`, `estimate`, `obs_per_cell` and
# `sac`, and every selection its `criterion`.  Unless a test says otherwise
# the regions are of the density kind from 3 k-means++ starts under seed 123.
# A layer is drawn under the seed its test states: the locations (x, then
# y), then the response.
# ---------------------------------------------------------------------------

# ---- Layers and helpers ------------------------------------------------------

.cs_seed <- function(seed, env = parent.frame())
  withr::local_seed(seed, .local_envir = env, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

.cs_points <- function(x, y, z, crs = 32617)
  sf::st_as_sf(data.frame(x = x, y = y, z = z), coords = c("x", "y"), crs = crs)

# `n` locations uniform on a square of side `side`, with the field of the
# engine's tests, sin(x / 2e4) + cos(y / 1.5e4) plus normal noise of standard
# deviation 1.5, or with a standard normal response and no field.
.cs_layer <- function(n, seed, side = 1e5, field = TRUE) {
  d <- local({
    .cs_seed(seed)
    x <- stats::runif(n, 0, side); y <- stats::runif(n, 0, side)
    z <- if (field) sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(n, 0, 1.5)
         else stats::rnorm(n)
    data.frame(x = x, y = y, z = z)
  })
  .cs_points(d$x, d$y, d$z)
}

# A profile with every setting that has a default the switch moves written
# out.  What a test varies it names.
.cs_profile <- function(pts, ..., response = "z", partition = "density", nstart = 3,
                        sac = "auto", estimate = "pooled", obs_per_cell = 4.5, seed = 123)
  resolution_profile(pts, response, ..., partition = partition, nstart = nstart, sac = sac,
                     estimate = estimate, obs_per_cell = obs_per_cell, seed = seed)

# A call with its warnings kept and muffled.
.cs_catch <- function(expr) {
  got <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    got[[length(got) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = got)
}
.cs_classes <- function(r) vapply(r$warnings, function(w) class(w)[1L], character(1))
.cs_of <- function(r, class) Filter(function(w) inherits(w, class), r$warnings)
.cs_fields <- function(w) unclass(w)[setdiff(names(unclass(w)), c("message", "call"))]
.cs_error <- function(expr) tryCatch({ suppressWarnings(expr); NULL }, error = function(e) e)
.cs_print <- function(x) utils::capture.output(print(x))

.cs_store <- new.env(parent = emptyenv())
.cs_keep <- function(key, expr) {
  if (is.null(.cs_store[[key]])) .cs_store[[key]] <- expr
  .cs_store[[key]]
}
.cs_sizes <- c(18, 12, 9, 6, 4.5, 3, 2)
# 600 rows of the field (seed 4401), and the profiles most tests read.
.cs_600    <- function() .cs_keep("l600", .cs_layer(600, 4401))
.cs_pooled <- function() .cs_keep("pooled", .cs_profile(.cs_600()))
.cs_pooled7 <- function() .cs_keep("pooled7", .cs_profile(.cs_600(), obs_per_cell = .cs_sizes))
.cs_mean7  <- function()
  .cs_keep("mean7", .cs_profile(.cs_600(), estimate = "mean", obs_per_cell = .cs_sizes))
# A profile without an estimate on the automatic ladder, as a profile is
# built today: k-means, a `sac` the caller supplies so that no variogram is
# fitted.
.cs_plain <- function()
  .cs_keep("plain", suppressWarnings(
    .cs_profile(.cs_600(), n_levels = 6L, partition = "kmeans", nstart = 5, sac = 3e4,
                estimate = NULL, obs_per_cell = NULL)))
# Three levels the caller set, with an estimate.
.cs_levels <- function()
  .cs_keep("levels", suppressWarnings(
    .cs_profile(.cs_600(), levels = c(50, 133, 300), sac = 3e4, obs_per_cell = NULL)))

# The stations layer of the ladder's tests: 3,000 rows (seed 301) at 100
# stations (seed 307).
.cs_stations <- function()
  .cs_keep("stations", {
    pts <- .cs_layer(3000, 301)
    at  <- local({ .cs_seed(307); sample.int(100L, 3000L, replace = TRUE) })
    sf::st_geometry(pts) <- sf::st_geometry(pts)[at]
    pts
  })
# White noise on 300 rows (seed 4001), and its profile as a profile is built
# today: k-means, the automatic ladder, a `sac` the caller supplies, no
# estimate.
.cs_noise <- function() .cs_keep("noise", .cs_layer(300, 4001, field = FALSE))
.cs_noise_plain <- function()
  .cs_keep("noise_plain", suppressWarnings(
    .cs_profile(.cs_noise(), n_levels = 6L, partition = "kmeans", nstart = 5, sac = 3e4,
                estimate = NULL, obs_per_cell = NULL)))
.cs_geoms <- function(g) unname(vapply(g$layers, function(l) class(l$geom)[1L], character(1)))
.cs_panels <- function(g)
  as.character(unique(g$layers[[which(.cs_geoms(g) == "GeomPoint")[1L]]]$data$panel))
# What ggplot2 says while a plot is drawn.
.cs_drawn <- function(g) {
  said <- character(0)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  withCallingHandlers(print(g), message = function(m) {
    said <<- c(said, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  said
}

.cs_present <- c("best", "flat", "criterion", "value", "at_ceiling", "at_floor", "edge",
                 "n_levels", "values", "structured", "one_cell", "structure_p", "range_floor",
                 "note", "seeds")
.cs_added <- c("estimate", "obs_per_cell", "obs_per_cell_asked", "scored_per_cell", "loo",
               "loo_one_cell", "share_alone", "own_weight", "count_rule")
.cs_record <- c("version", "rule", "estimate", "cells", "obs_per_cell", "obs_per_cell_asked",
                "response", "n", "n_points", "n_fit", "seed", "features", "scored",
                "partition", "crs", "scored_cells")


# ---- What exists today is what it was ----------------------------------------

test_that("a selection by a present criterion from a profile without an estimate has today's names", {
  p <- .cs_plain()
  expect_null(attr(p, "estimate"))
  for (cr in c("cp", "moran_z")) {
    s <- select_resolution(p, criterion = cr, diagnostic = TRUE)
    expect_identical(names(s), .cs_present, info = cr)
    expect_null(s$count_rule, info = cr)
    expect_null(attr(s$seeds, "count_rule"), info = cr)
    expect_identical(s$criterion, cr)
  }
  # A profile of the density kind names its kind last, and nothing else.
  d <- suppressWarnings(.cs_profile(.cs_600(), n_levels = 6L, sac = 3e4, estimate = NULL,
                                    obs_per_cell = NULL))
  s <- select_resolution(d, criterion = "cp")
  expect_identical(names(s), c(.cs_present, "partition"))
  expect_null(attr(s$seeds, "count_rule"))
  # The names the baseline's selections have.
  base <- readRDS(test_path("fixtures", "core-frozen.rds"))
  expect_identical(names(base$records$select_default_cp$value), .cs_present)
  expect_identical(names(base$records$select_density_cp$value), c(.cs_present, "partition"))
})

test_that("a profile of today, read by cp, by \"auto\" or as a count, gives the baseline's results", {
  skip_on_cran()
  skip_if_not_installed("gstat")
  skip_if_not_installed("dplyr")
  base <- readRDS(test_path("fixtures", "core-frozen.rds"))
  # The baseline's calls with the defaults of today written out, and the new
  # readings of the same objects.  Each is held to the baseline's record of
  # the call named beside it.
  prof <- function(response)
    bquote(resolution_profile(pts, .(response), estimate = NULL, obs_per_cell = NULL,
                              sac = NULL, partition = "kmeans", nstart = 25L))
  calls <- list(
    profile_default      = prof("z"),
    select_default_cp    = quote(select_resolution(profile_default, criterion = "cp")),
    auto_default         = quote(select_resolution(profile_default, criterion = "auto")),
    subset_cp            = quote(select_resolution(profile_default[, c("levels", "cp")],
                                                   criterion = "cp")),
    subset_auto          = quote(select_resolution(profile_default[, c("levels", "cp")],
                                                   criterion = "auto")),
    seeds_profile        = quote(get_voronoi_seeds(method = "kmeans", n = profile_default,
                                                   sample_points = pts)),
    kmeans_profile       = quote(voronoi_seeds_kmeans(pts, profile_default)),
    cells_hex_profile    = quote(build_tessellation(pts, boundary = bnd, method = "hex",
                                                    approx_n_cells = profile_default)),
    profile_no_value     = prof("none"),
    count_seeds_no_value = quote(get_voronoi_seeds(method = "kmeans", n = profile_no_value,
                                                   sample_points = pts)),
    profile_noise        = prof("noise"),
    count_seeds_noise    = quote(get_voronoi_seeds(method = "kmeans", n = profile_noise,
                                                   sample_points = pts)))
  against <- c(profile_default = "profile_default", select_default_cp = "select_default_cp",
               auto_default = "select_default_cp", subset_cp = "select_column_subset_cp",
               subset_auto = "select_column_subset_cp", seeds_profile = "seeds_profile",
               kmeans_profile = "kmeans_profile", cells_hex_profile = "cells_hex_profile",
               profile_no_value = "profile_no_value",
               count_seeds_no_value = "count_seeds_no_value", profile_noise = "profile_noise",
               count_seeds_noise = "count_seeds_noise")
  expect_true(all(against %in% names(base$records)))
  got <- core_frozen_run(base$data, calls)
  same_platform <- core_frozen_same_platform(core_frozen_platform(), base$info)
  for (nm in names(against)) {
    want <- base$records[[against[[nm]]]]
    d <- core_frozen_differences(got[[nm]], want)
    expect(length(d) == 0L,
           paste(c(sprintf("`%s` is not the baseline's `%s`.", nm, against[[nm]]), d),
                 collapse = "\n"))
    if (same_platform)
      for (part in c("value", "conditions", "printed", "log"))
        expect_identical(got[[nm]][[part]], want[[part]],
                         label = sprintf("the %s of `%s`", part, nm))
  }
  # What those records hold: the error of a profile with no finite
  # criterion, and the two warnings of a count on a layer without structure.
  expect_identical(got$count_seeds_no_value$conditions[[1L]]$kind, "error")
  expect_match(got$count_seeds_no_value$conditions[[1L]]$message,
               "is a resolution profile with no criterion finite at any level", fixed = TRUE)
  expect_true("spatialkit_no_structure" %in%
                unlist(lapply(got$count_seeds_noise$conditions, function(c) c$class)))
  expect_false("count_rule" %in% names(got$select_default_cp$value))
  expect_null(attr(got$seeds_profile$value, "count_rule"))
})


test_that("a column called loo on a profile without an estimate is read only when it is named", {
  skip_if_not_installed("ggplot2")
  label <- "Leave-one-out error (lower is better)"
  # A column of the caller's own beside the criteria of a profile built as
  # profiles are today: the table and the plot are the ones without it.
  p <- .cs_plain()
  q <- p; q$loo <- (q$levels - 30)^2
  expect_null(attr(q, "estimate"))
  expect_false(.loo_by_default(q))
  expect_false(.loo_by_default(p))
  expect_identical(summary(q), summary(p))
  expect_identical(.cs_print(summary(q)), .cs_print(summary(p)))
  gq <- plot(q); gp <- plot(p)
  expect_identical(.cs_panels(gq), .cs_panels(gp))
  expect_false(label %in% .cs_panels(gq))
  expect_identical(.cs_geoms(gq), .cs_geoms(gp))
  expect_identical(gq$labels$caption, gp$labels$caption)
  # Named, it is read as the column of that name is.
  expect_identical(summary(q, criteria = "loo")$criterion, "loo")
  expect_identical(.cs_panels(plot(q, criteria = "loo")), label)
  # The column resolution_profile() wrote is read without being named: on the
  # profile and its row subsets, which keep the estimate, and on the column
  # subset that kept `loo` and no other count criterion.
  l <- .cs_levels()
  expect_true(.loo_by_default(l))
  expect_true(.loo_by_default(l[1:2, ]))
  expect_true("loo" %in% summary(l)$criterion)
  only <- l[, c("levels", "loo")]
  expect_null(attr(only, "estimate"))
  expect_true(.loo_by_default(only))
  expect_identical(summary(only)$criterion, "loo")
  expect_identical(.cs_panels(plot(only)), label)
  # A column subset that kept `cp` has lost the estimate, and is read as a
  # profile without one is.
  both <- l[, c("levels", "cp", "loo")]
  expect_false(.loo_by_default(both))
  expect_identical(summary(both)$criterion, "cp")
  expect_false(label %in% .cs_panels(plot(both)))
  # A subset whose `loo` is NA at every level is told what a subset without
  # a criterion column is told today.
  none <- only; none$loo <- NA_real_
  expect_identical(conditionMessage(.cs_error(summary(none))),
                   paste0("summary.resolution_profile(): `object` carries none of the criterion ",
                          "columns (cp, reliability, elbow, moran_z); it is a column subset of a ",
                          "resolution profile. Pass the profile resolution_profile() returned."))
})


# ---- "loo" ------------------------------------------------------------------

test_that("\"loo\" returns the level pick_rung() picks, with the flat region around it", {
  for (p in list(.cs_mean7(), .cs_pooled7())) {
    what <- attr(p, "estimate")$value
    expect_identical(p$levels, vapply(.cs_sizes, function(s) rung_k(600L, s), integer(1)),
                     info = what)
    s <- select_resolution(p, criterion = "loo", tol = 0.02)
    expect_identical(s$best, p$levels[pick_rung(p$loo)], info = what)
    expect_identical(s$criterion, "loo")
    expect_identical(s$value, min(p$loo))
    expect_identical(s$flat, p$levels[p$loo <= min(p$loo) * 1.02], info = what)
    expect_identical(s$values, stats::setNames(p$loo, p$levels))
    expect_identical(s$n_levels, 7L)
    expect_identical(s$range_floor, NA_integer_)
    expect_identical(s$one_cell, attr(p, "estimate")$loo_one_cell)
    expect_identical(s$estimate, what)
    expect_identical(names(s), c(.cs_present, .cs_added, "partition"))
  }
})

test_that("of levels with the same leave-one-out error \"loo\" takes the fewest cells", {
  p <- .cs_mean7()
  lv <- p$levels
  q <- p
  q$loo <- c(5, 3, 4, 3, 6, 7, 8)
  s <- select_resolution(q, criterion = "loo")
  expect_identical(s$best, lv[2L])
  expect_identical(s$best, lv[pick_rung(q$loo)])
  expect_identical(s$flat, lv[c(2L, 4L)])
  # Whatever order the rows of the profile are in: pick_rung() reads its
  # rungs coarse first.
  r <- q[order(-q$levels), ]
  expect_identical(r$levels, rev(lv))
  expect_identical(select_resolution(r, criterion = "loo")$best, lv[2L])
  # The column is no longer the one the references were recorded with, so
  # nothing is said about structure.
  expect_identical(s$structured, NA)
  expect_identical(s$one_cell, NA_real_)
  expect_identical(s$structure_p, NA_real_)
  expect_null(s$note)
})

test_that("the flat region and the edge of \"loo\" are the hand computation's", {
  p <- .cs_mean7()
  lv <- p$levels
  q <- p
  q$loo <- c(5, 4, 3, 2.99, 3.5, 4, 5)
  s <- select_resolution(q, criterion = "loo", tol = 0.02)
  expect_identical(s$best, lv[4L])
  expect_identical(s$flat, lv[3:4])              # 3 <= 2.99 * 1.02 = 3.0498
  expect_identical(c(s$at_ceiling, s$at_floor), c(FALSE, FALSE))
  expect_identical(s$edge, NA_character_)
  expect_identical(select_resolution(q, criterion = "loo", tol = 0.2)$flat, lv[3:5])
  expect_identical(select_resolution(q, criterion = "loo", tol = 0)$flat, lv[4L])
  # The ends of a ladder of sizes are the sizes offered.
  q$loo <- 7:1
  s <- select_resolution(q, criterion = "loo")
  expect_identical(s$best, lv[7L])
  expect_true(s$at_ceiling); expect_false(s$at_floor)
  expect_identical(s$edge, "the finest size offered (obs_per_cell)")
  q$loo <- 1:7
  s <- select_resolution(q, criterion = "loo")
  expect_identical(s$best, lv[1L])
  expect_true(s$at_floor); expect_false(s$at_ceiling)
  expect_identical(s$edge, "the coarsest size offered (obs_per_cell)")
  # Where the column stops short of the ladder's ends, that is the edge.
  q$loo <- c(NA, 1, 2, 3, 4, 5, NA)
  expect_identical(select_resolution(q, criterion = "loo")$edge,
                   "the first level the criterion is computable at")
  q$loo <- c(NA, 5, 4, 3, 2, 1, NA)
  expect_identical(select_resolution(q, criterion = "loo")$edge,
                   "the last level the criterion is computable at")
  # On a ladder the caller set with `levels` the ends are that ladder's.
  l <- .cs_levels()
  l$loo <- c(1, 2, 3)
  expect_identical(select_resolution(l, criterion = "loo")$edge, "the first level of the ladder")
})

test_that("the edge labels and their hints: a ladder of sizes has its own, every other keeps its own", {
  sized <- list(obs_per_cell = c(9, 4.5), ceiling = 133L, floor = 2L, supported = TRUE,
                range = NA_real_, range_floor = FALSE, ceiling_from = "obs_per_cell")
  lv <- c(67L, 133L)
  expect_identical(.ladder_edge(133L, lv, lv, sized), "the finest size offered (obs_per_cell)")
  expect_identical(.ladder_edge(67L, lv, lv, sized), "the coarsest size offered (obs_per_cell)")
  expect_identical(.ladder_edge(133L, 133L, 133L, sized), "the only level scored")
  expect_identical(.ladder_edge(100L, c(67L, 100L, 133L), c(67L, 100L, 133L), sized),
                   NA_character_)
  # A finest size the profile held to the cells k-means can place is the
  # layer's bound, not the caller's: the label of the support ceiling, which
  # has no hint about sizes.  The coarsest size is still the caller's.
  held <- utils::modifyList(sized, list(ceiling_from = "distinct locations", n_distinct = 133L))
  expect_identical(.ladder_edge(133L, lv, lv, held),
                   "the support ceiling (as many cells as the distinct locations allow)")
  expect_identical(.ladder_edge(67L, lv, lv, held), "the coarsest size offered (obs_per_cell)")
  short <- utils::modifyList(held, list(n_distinct = 134L))
  expect_identical(.ladder_edge(133L, lv, lv, short),
                   "the support ceiling (one short of the points, the most k-means can fit)")
  expect_identical(.edge_hint("the support ceiling (as many cells as the distinct locations allow)"),
                   "")
  # The same bounds without sizes: the labels of today.
  auto <- utils::modifyList(sized, list(obs_per_cell = NULL, ceiling_from = "min_cell_n"))
  expect_identical(.ladder_edge(133L, lv, lv, auto), "the support ceiling (n / min_cell_n)")
  expect_identical(.ladder_edge(67L, lv, lv, auto), "the first level of the ladder")
  expect_identical(.ladder_edge(133L, lv, lv, NULL), "the last level of the ladder")
  hint <- " `obs_per_cell` sets the sizes offered; how small a cell may be is the caller's choice."
  expect_identical(.edge_hint("the finest size offered (obs_per_cell)"), hint)
  expect_identical(.edge_hint("the coarsest size offered (obs_per_cell)"), hint)
  expect_identical(.edge_hint("the support ceiling (n / min_cell_n)"),
                   " Lower min_cell_n to see whether the criterion keeps going.")
  expect_identical(.edge_hint("the only level scored"),
                   " One level leaves nothing to compare it with.")
  expect_identical(.edge_hint(NA_character_), "")
})

test_that("the note for pooled values appears only at the finest level of a pooled profile", {
  text <- paste0("the lowest leave-one-out error of the pooled value is at the finest level ",
                 "scored, so the end of the ladder chose the count. criterion = \"size\" sets ",
                 "it by a size.")
  pp <- .cs_pooled7(); pm <- .cs_mean7()
  finest <- pp; finest$loo <- 7:1
  s <- select_resolution(finest, criterion = "loo")
  expect_identical(s$note[["finest"]], text)
  # The finest level scored, where the column stops short of the ladder's end.
  short <- pp; short$loo <- c(6, 5, 4, 3, 2, 1, NA)
  expect_identical(select_resolution(short, criterion = "loo")$note[["finest"]], text)
  # Not at another level, and not for plain means.
  inner <- pp; inner$loo <- c(5, 4, 3, 2.99, 3.5, 4, 5)
  expect_false("finest" %in% names(select_resolution(inner, criterion = "loo")$note))
  coarse <- pp; coarse$loo <- 1:7
  expect_false("finest" %in% names(select_resolution(coarse, criterion = "loo")$note))
  plain <- pm; plain$loo <- 7:1
  expect_false("finest" %in% names(select_resolution(plain, criterion = "loo")$note))
  # On the profile as it was built the note follows the level picked.
  s <- select_resolution(pp, criterion = "loo")
  expect_identical("finest" %in% names(s$note), s$best == max(pp$levels))
  # A ladder of one size: the level is the finest scored.
  expect_identical(select_resolution(.cs_pooled(), criterion = "loo")$note[["finest"]], text)
  # The finest level, whatever order the rows of the profile are in: with
  # the rows fine first it is the first row, and the coarsest is the last.
  turned <- finest[order(-finest$levels), ]
  expect_identical(turned$levels[1L], max(pp$levels))
  expect_identical(select_resolution(turned, criterion = "loo")$note[["finest"]], text)
  turned <- coarse[order(-coarse$levels), ]
  s <- select_resolution(turned, criterion = "loo")
  expect_identical(s$best, min(pp$levels))
  expect_false("finest" %in% names(s$note))
})

test_that("\"loo\" runs on a column subset that keeps `levels` and `loo`", {
  p <- .cs_pooled7()
  sub <- p[, c("levels", "loo")]
  expect_null(attr(sub, "estimate"))
  s <- select_resolution(sub, criterion = "loo")
  full <- select_resolution(p, criterion = "loo")
  expect_identical(s$best, full$best)
  expect_identical(s$flat, full$flat)
  expect_identical(s$structured, NA)
  expect_identical(s$estimate, NA_character_)
  expect_identical(s$one_cell, NA_real_)
  expect_identical(s$edge, .ladder_edge(s$best, sub$levels, sub$levels, NULL))
  expect_null(s$seeds)
  expect_identical(s$loo, full$loo)
  expect_identical(s$loo_one_cell, NA_real_)
  expect_identical(c(s$obs_per_cell, s$obs_per_cell_asked, s$scored_per_cell, s$share_alone,
                     s$own_weight), rep(NA_real_, 5L))
  rec <- s$count_rule
  expect_identical(names(rec), .cs_record)
  expect_identical(rec$rule, "loo")
  expect_identical(rec$estimate, NA_character_)
  expect_identical(list(rec$n, rec$n_points, rec$n_fit), list(NA_integer_, NA_integer_, NA_integer_))
  expect_null(rec$scored)
  # The ends of the column are the ladder's, with no bounds to name them.
  sub$loo <- 7:1
  expect_identical(select_resolution(sub, criterion = "loo")$edge, "the last level of the ladder")
  expect_null(select_resolution(sub, criterion = "loo")$note)
})

test_that("\"loo\" says what it lacks: the column, a response that varies, a scored level", {
  e <- .cs_error(select_resolution(.cs_plain(), criterion = "loo"))
  expect_identical(conditionMessage(e),
                   paste0("select_resolution(): `profile` has no `loo` column; build it with ",
                          "resolution_profile(estimate = \"mean\") or estimate = \"pooled\"."))
  expect_null(conditionCall(e))
  # A subset that kept the column and lost the levels is told so in the
  # present words.
  e <- .cs_error(select_resolution(.cs_mean7()[, "loo", drop = FALSE], criterion = "loo"))
  expect_identical(conditionMessage(e),
                   paste0("select_resolution(): `profile` has no `levels` column; it is a ",
                          "column subset of a resolution profile. Pass the profile ",
                          "resolution_profile() returned."))
  # A constant response: the present error, with `loo` named.
  pts <- .cs_600(); pts$k <- 0.1
  k <- suppressWarnings(.cs_profile(pts, response = "k", estimate = "mean",
                                    obs_per_cell = c(9, 4.5)))
  expect_true(all(is.na(k$loo)))
  e <- .cs_error(select_resolution(k, criterion = "loo"))
  expect_identical(conditionMessage(e),
                   paste0("select_resolution(): `loo` is NA at every level because the response ",
                          "is constant: its total sum of squares is rounding error at most, so ",
                          "no cell count fits it better than another."))
  # NA at every level for another reason.
  q <- .cs_mean7(); q$loo <- NA_real_
  e <- .cs_error(select_resolution(q, criterion = "loo"))
  expect_identical(conditionMessage(e),
                   paste0("select_resolution(): `loo` is NA at every level (it needs a level at ",
                          "which some cell holds two rows with a response)."))
})

test_that("an optimum of \"loo\" on an end of the sizes raises spatialkit_bound_chosen downstream", {
  pts <- .cs_600()
  q <- .cs_mean7(); q$loo <- 1:7
  s <- select_resolution(q, criterion = "loo")
  r <- .cs_catch(get_voronoi_seeds(method = "kmeans", n = s, sample_points = pts))
  w <- .cs_of(r, "spatialkit_bound_chosen")
  expect_length(w, 1L)
  expect_identical(.cs_fields(w[[1L]]),
                   list(edge = "the coarsest size offered (obs_per_cell)", criterion = "loo",
                        n = 33L))
  expect_identical(conditionMessage(w[[1L]]),
                   paste0("get_voronoi_seeds(): `n` is 33 cells, the optimum of `loo`, and that ",
                          "optimum is the coarsest size offered (obs_per_cell): the bound chose ",
                          "the count, not the criterion. `obs_per_cell` sets the sizes offered; ",
                          "how small a cell may be is the caller's choice. Pass a number as `n` ",
                          "to set the count yourself."))
  expect_identical(nrow(r$value), 33L)
  expect_identical(attr(r$value, "n_from"),
                   paste0("select_resolution(criterion = \"loo\"), an optimum at the coarsest ",
                          "size offered (obs_per_cell)"))
})


# ---- "size" -----------------------------------------------------------------

test_that("\"size\" returns the count of the size when the profile has that level", {
  p <- .cs_levels()
  expect_null(attr(p, "bounds")$obs_per_cell)
  r <- .cs_catch(select_resolution(p, criterion = "size", obs_per_cell = 4.5))
  expect_length(r$warnings, 0L)
  s <- r$value
  expect_identical(s$best, .count_at(600L, 4.5))
  expect_identical(s$best, rung_k(600L, 4.5))
  expect_identical(s$best, 133L)
  expect_identical(s$criterion, "size")
  expect_identical(s$flat, 133L)
  expect_identical(s$value, 600 / 133)
  expect_identical(s$values, stats::setNames(600 / p$levels, p$levels))
  expect_identical(s$n_levels, 3L)
  expect_identical(list(s$at_ceiling, s$at_floor, s$edge), list(NA, NA, NA_character_))
  expect_identical(s$range_floor, NA_integer_)
  expect_identical(s$obs_per_cell, 600 / 133)
  expect_identical(s$obs_per_cell_asked, 4.5)
  expect_identical(s$estimate, "pooled")
  expect_identical(s$scored_per_cell, p$scored_per_cell[2L])
  expect_identical(s$loo, p$loo[2L])
  expect_identical(s$loo_one_cell, attr(p, "estimate")$loo_one_cell)
  expect_identical(s$one_cell, attr(p, "estimate")$loo_one_cell)
  expect_identical(s$share_alone, p$share_alone[2L])
  expect_identical(s$own_weight, p$own_weight[2L])
  expect_identical(s$note, c(size = paste0("the count is set by obs_per_cell = 4.5, the ",
                                           "caller's choice; no criterion chose it.")))
  expect_identical(nrow(s$seeds), 133L)
  expect_identical(names(s), c(.cs_present, .cs_added, "partition"))
  # Without a size of its own the profile is read at the default size.
  expect_identical(.OBS_PER_CELL_DEFAULT, 4.5)
  expect_identical(select_resolution(p, criterion = "size"), s)
  # Each of the three levels is the count of a size.
  for (size in c(12, 4.5, 2)) {
    got <- .cs_catch(select_resolution(p, criterion = "size", obs_per_cell = size))
    expect_length(got$warnings, 0L)
    expect_identical(got$value$best, rung_k(600L, size), info = format(size))
  }
})

test_that("a size the profile was not built at gets the nearest level on the log scale, and a warning", {
  p <- .cs_levels()
  r <- .cs_catch(select_resolution(p, criterion = "size", obs_per_cell = 6))
  expect_identical(.cs_classes(r), "spatialkit_size_not_reached")
  w <- r$warnings[[1L]]
  expect_identical(conditionMessage(w),
                   paste0("select_resolution(): obs_per_cell = 6 asks for 100 cells, which is not ",
                          "a level of this profile; the nearest level, 133 cells (4.51 rows per ",
                          "cell), is used. Build the profile with obs_per_cell = 6 to score that ",
                          "count."))
  expect_identical(.cs_fields(w), list(obs_per_cell = 6, wanted = 100L, used = 133L,
                                       reason = "not a level"))
  expect_null(conditionCall(w))
  s <- r$value
  expect_identical(s$best, 133L)
  expect_identical(s$obs_per_cell_asked, 6)
  expect_identical(s$obs_per_cell, 600 / 133)
  expect_identical(s$note[["size"]],
                   "the count is set by obs_per_cell = 6, the caller's choice; no criterion chose it.")
  # Nearest on the log scale, not in cells: 40 cells lie 30 from 10 and 60
  # from 100, and a factor of 4 from 10 against 2.5 from 100.  Of two levels
  # as near, the one with fewer cells: 4 lies a factor of 2 from 2 and from 8.
  q <- suppressWarnings(.cs_profile(.cs_600(), levels = c(2, 8, 10, 100), sac = 3e4,
                                    obs_per_cell = NULL))
  expect_identical(.count_at(600L, 15), 40L)
  expect_identical(suppressWarnings(select_resolution(q, criterion = "size", obs_per_cell = 15))$best,
                   100L)
  expect_identical(.count_at(600L, 150), 4L)
  r <- .cs_catch(select_resolution(q, criterion = "size", obs_per_cell = 150))
  expect_identical(r$value$best, 2L)
  expect_identical(.cs_fields(r$warnings[[1L]])[c("wanted", "used")], list(wanted = 4L, used = 2L))
  # A level whose partition was not built does not stand for a size: the
  # count of 4.5 is 133, and with no partition there the nearest built level
  # is 300 (a factor of 2.26, against 2.66 for 50).
  failed <- p; failed$wss[2L] <- NA_real_
  r <- .cs_catch(select_resolution(failed, criterion = "size", obs_per_cell = 4.5))
  expect_identical(r$value$best, 300L)
  expect_identical(.cs_fields(r$warnings[[1L]]),
                   list(obs_per_cell = 4.5, wanted = 133L, used = 300L, reason = "not a level"))
  none <- p; none$wss <- NA_real_
  expect_error(suppressWarnings(select_resolution(none, criterion = "size")),
               "no level of `profile` has a partition", fixed = TRUE)
})

test_that("the level a profile built for a size stands for it only while the profile has it", {
  p <- .cs_pooled7()
  expect_identical(attr(p, "bounds")$size_levels[["4.5"]], 133L)
  words <- paste0("select_resolution(): obs_per_cell = 4.5 asks for 133 cells, which is not a ",
                  "level of this profile; the nearest level, 100 cells (6 rows per cell), is ",
                  "used. Build the profile with obs_per_cell = 4.5 to score that count.")
  # The fit at that level failed.  The nearest built level is 100 cells (a
  # factor of 1.33, against 1.50 for 200), with the warning of a size the
  # profile was not built at.
  failed <- p; failed$wss[5L] <- NA_real_
  # A row subset that has lost the level.
  cut <- p[-5L, ]
  expect_false(133L %in% cut$levels)
  expect_identical(attr(cut, "bounds")$size_levels[["4.5"]], 133L)
  for (q in list(failed, cut)) {
    r <- .cs_catch(select_resolution(q, criterion = "size", obs_per_cell = 4.5))
    expect_identical(.cs_classes(r), "spatialkit_size_not_reached")
    expect_identical(conditionMessage(r$warnings[[1L]]), words)
    expect_identical(.cs_fields(r$warnings[[1L]]),
                     list(obs_per_cell = 4.5, wanted = 133L, used = 100L, reason = "not a level"))
    expect_identical(r$value$best, 100L)
    expect_identical(names(r$value$note), "size")
    expect_identical(nrow(r$value$seeds), 100L)
  }
  # The other sizes still stand for their levels, and raise nothing.
  r <- .cs_catch(select_resolution(cut, criterion = "size", obs_per_cell = 6))
  expect_length(r$warnings, 0L)
  expect_identical(r$value$best, 100L)
})

test_that("with obs_per_cell = NULL \"size\" takes the profile's own size when it has one", {
  pts <- .cs_600()
  six <- .cs_profile(pts, obs_per_cell = 6)
  s <- select_resolution(six, criterion = "size")
  expect_identical(s$obs_per_cell_asked, 6)
  expect_identical(s$best, 100L)
  expect_identical(s$obs_per_cell, 6)
  expect_identical(select_resolution(six, criterion = "size", obs_per_cell = 6), s)
  # Several sizes are not one: the default size, which is one of them here.
  s <- select_resolution(.cs_pooled7(), criterion = "size")
  expect_identical(s$obs_per_cell_asked, 4.5)
  expect_identical(s$best, 133L)
  # Any of the profile's sizes is the level built for it, with no warning.
  for (size in .cs_sizes) {
    r <- .cs_catch(select_resolution(.cs_pooled7(), criterion = "size", obs_per_cell = size))
    expect_length(r$warnings, 0L)
    expect_identical(r$value$best, rung_k(600L, size), info = format(size))
    expect_length(r$value$note, 1L)
  }
})

test_that("\"size\" divides the rows with a response, and every point without one", {
  pts <- .cs_600()
  some <- pts; some$z[1:200] <- NA
  p <- .cs_profile(some, obs_per_cell = 4)
  s <- select_resolution(p, criterion = "size")
  expect_identical(s$best, rung_k(400L, 4))
  expect_identical(s$value, 4)
  expect_identical(s$count_rule$n, 400L)
  expect_identical(s$count_rule$n_points, 600L)
  # A profile of the geometry alone: the size counts the points, and nothing
  # is said of an estimate or of structure.
  geo <- resolution_profile(pts, partition = "density", nstart = 3, sac = "auto",
                            estimate = NULL, obs_per_cell = 6, seed = 123)
  s <- select_resolution(geo, criterion = "size")
  expect_identical(s$best, 100L)
  expect_identical(s$obs_per_cell_asked, 6)
  expect_identical(s$structured, NA)
  expect_identical(c(s$one_cell, s$structure_p, s$loo, s$loo_one_cell, s$scored_per_cell),
                   rep(NA_real_, 5L))
  expect_identical(s$estimate, NA_character_)
  expect_identical(names(s), c(.cs_present, .cs_added, "partition"))
  expect_identical(s$count_rule$estimate, NA_character_)
  expect_identical(s$count_rule$response, NA_character_)
  expect_null(s$count_rule$scored)
  expect_identical(s$count_rule$n, 600L)
  expect_identical(unname(s$note), "the count is set by obs_per_cell = 6, the caller's choice; no criterion chose it.")
  # A profile built as profiles are today.
  plain <- .cs_plain()
  r <- .cs_catch(select_resolution(plain, criterion = "size", obs_per_cell = 6))
  s <- r$value
  expect_identical(s$structured, NA)
  expect_identical(s$estimate, NA_character_)
  expect_true(s$best %in% plain$levels)
  expect_identical(names(s), c(.cs_present, .cs_added))
  expect_identical(s$count_rule$partition$kind, "kmeans")
  if (!(100L %in% plain$levels))
    expect_identical(.cs_classes(r), "spatialkit_size_not_reached")
})

test_that("\"size\" on a column subset is the error that names the row count", {
  p <- .cs_pooled()
  e <- .cs_error(select_resolution(p[, c("levels", "loo")], criterion = "size"))
  expect_identical(conditionMessage(e),
                   paste0("select_resolution(): criterion = \"size\" needs the number of rows the ",
                          "size divides, which a profile keeps in attr(, \"bounds\"); this object ",
                          "has none (it is a column subset of a resolution profile). Pass the ",
                          "profile resolution_profile() returned."))
  expect_null(conditionCall(e))
  expect_error(select_resolution(p[, c("levels", "loo")], criterion = "auto"), "has no `cp` column")
  # A present criterion still runs on the column subset that keeps its column.
  l <- .cs_levels()
  s <- select_resolution(l[, c("levels", "cp")], criterion = "cp")
  expect_identical(s$best, select_resolution(l, criterion = "cp")$best)
  expect_identical(names(s), .cs_present)
})

test_that("no bound chooses a size: no spatialkit_bound_chosen downstream", {
  pts <- .cs_600()
  p <- .cs_levels()
  # The last level of the ladder, and the only level of a ladder of one size.
  top <- select_resolution(p, criterion = "size", obs_per_cell = 2)
  expect_identical(top$best, max(p$levels))
  one <- select_resolution(.cs_pooled(), criterion = "size")
  for (s in list(top, one)) {
    r <- .cs_catch(get_voronoi_seeds(method = "kmeans", n = s, sample_points = pts))
    expect_false("spatialkit_bound_chosen" %in% .cs_classes(r))
    expect_identical(nrow(r$value), s$best)
    r <- .cs_catch(voronoi_seeds_kmeans(pts, s))
    expect_false("spatialkit_bound_chosen" %in% .cs_classes(r))
  }
  # The same level read by a criterion that searched is an optimum on a bound.
  l <- p; l$loo <- c(3, 2, 1)
  r <- .cs_catch(get_voronoi_seeds(method = "kmeans", n = select_resolution(l, criterion = "loo"),
                                   sample_points = pts))
  expect_true("spatialkit_bound_chosen" %in% .cs_classes(r))
})


# ---- A size the profile held ---------------------------------------------------

test_that("a size the profile held to the cells k-means can place is returned without a second warning", {
  # The stations layer of the ladder's tests: 3,000 rows (seed 301) at 100
  # stations (seed 307).
  pts <- .cs_layer(3000, 301)
  at  <- local({ .cs_seed(307); sample.int(100L, 3000L, replace = TRUE) })
  sf::st_geometry(pts) <- sf::st_geometry(pts)[at]
  built <- .cs_catch(.cs_profile(pts))
  expect_identical(.cs_classes(built),
                   c("spatialkit_size_not_reached", "spatialkit_repeated_locations"))
  p <- built$value
  expect_identical(p$levels, 100L)
  r <- .cs_catch(select_resolution(p, criterion = "size"))
  expect_length(r$warnings, 0L)
  s <- r$value
  expect_identical(s$best, 100L)
  expect_identical(s$obs_per_cell_asked, 4.5)
  expect_identical(s$obs_per_cell, 30)
  expect_identical(s$note[["held"]],
                   paste0("the profile holds this size to 100 cells, the most k-means can place ",
                          "on this layer; obs_per_cell = 4.5 asks for 667."))
  expect_identical(names(s$note)[1:2], c("size", "held"))
  expect_identical(select_resolution(p, criterion = "size", obs_per_cell = 4.5), s)
  # Neither the selection nor the profile, passed as the count, says it again.
  r <- .cs_catch(get_voronoi_seeds(method = "kmeans", n = s, sample_points = pts))
  expect_false("spatialkit_size_not_reached" %in% .cs_classes(r))
  expect_identical(nrow(r$value), 100L)
  r <- .cs_catch(get_voronoi_seeds(method = "kmeans", n = p, sample_points = pts))
  expect_false("spatialkit_size_not_reached" %in% .cs_classes(r))
  expect_identical(nrow(r$value), 100L)
  # A size it was not built at is the nearest level, with the warning.
  r <- .cs_catch(select_resolution(p, criterion = "size", obs_per_cell = 9))
  expect_identical(.cs_classes(r), "spatialkit_size_not_reached")
  expect_identical(r$value$best, 100L)
  expect_identical(.cs_fields(r$warnings[[1L]]),
                   list(obs_per_cell = 9, wanted = 333L, used = 100L, reason = "not a level"))
  expect_false("held" %in% names(r$value$note))
})


test_that("the seeds a size hands on are the centres of the level returned, held or nearest", {
  pts <- .cs_stations()
  p <- suppressWarnings(.cs_profile(pts))
  expect_identical(p$levels, 100L)
  centres <- unname(attr(p, "centres")[["100"]])
  # The level the profile held the size to (rule i), and the nearest level
  # for a size it was not built at (rule iii): 100 cells where 667 and 333
  # were asked for.
  held <- select_resolution(p, criterion = "size")
  near <- suppressWarnings(select_resolution(p, criterion = "size", obs_per_cell = 9))
  for (s in list(held, near)) {
    expect_identical(s$best, 100L)
    expect_identical(nrow(s$seeds), 100L)
    expect_identical(unname(sf::st_coordinates(s$seeds)[, 1:2]), centres)
    expect_identical(attr(s$seeds, "count_rule")$scored_cells, TRUE)
    expect_identical(attr(s$seeds, "count_rule")$cells, 100L)
    got <- suppressWarnings(get_voronoi_seeds(method = "kmeans", n = s, sample_points = pts))
    expect_identical(unname(sf::st_coordinates(got)[, 1:2]), centres)
  }
  # The stations' response has no structure cells of 30 rows can use.  The
  # warning at the seeding says why, and leaves out the note about the size
  # the profile held.
  expect_identical(names(held$note), c("size", "held", "structure"))
  r <- .cs_catch(get_voronoi_seeds(method = "kmeans", n = held, sample_points = pts))
  expect_identical(.cs_classes(r), "spatialkit_no_structure")
  expect_identical(conditionMessage(r$warnings[[1L]]),
                   paste0("get_voronoi_seeds(): `n` is 100 cells, set by obs_per_cell = 4.5, and ",
                          held$note[["structure"]],
                          " The count stands, since a size set it and no criterion; the cell ",
                          "values will say little that one mean for the layer does not."))
})

test_that("a count the layer cannot hold is said to be one, and no profile is advised for it", {
  nearest <- function(size, wanted, used, per, ending)
    sprintf(paste0("select_resolution(): obs_per_cell = %s asks for %d cells, which is not a ",
                   "level of this profile; the nearest level, %d cells (%s rows per cell), is ",
                   "used. %s"), size, wanted, used, per, ending)
  build <- function(size)
    sprintf("Build the profile with obs_per_cell = %s to score that count.", size)
  said <- function(p, size)
    conditionMessage(.cs_catch(select_resolution(p, criterion = "size",
                                                 obs_per_cell = size))$warnings[[1L]])
  # 3,000 rows at 100 stations: 9 rows a cell are 333 cells, and k-means
  # places 100.
  p <- suppressWarnings(.cs_profile(.cs_stations()))
  expect_identical(
    said(p, 9),
    nearest("9", 333L, 100L, "30",
            paste0("That count is more than k-means can place on the layer's 100 distinct ",
                   "locations (100 cells), so no profile scores it.")))
  # A count the stations can hold keeps the advice: 60 rows a cell are 50.
  expect_identical(said(p, 60), nearest("60", 50L, 100L, "30", build("60")))
  # No location repeated: one short of the points is the most.  600 cells
  # cannot be placed on 600 points, and 599 can.
  q <- .cs_pooled()
  expect_identical(
    said(q, 1),
    nearest("1", 600L, 133L, "4.51",
            paste0("That count is more than k-means can place on 600 points (599 cells), so ",
                   "no profile scores it.")))
  expect_identical(.count_at(600L, 1.002), 599L)
  expect_identical(said(q, 1.002), nearest("1.002", 599L, 133L, "4.51", build("1.002")))
  # A profile fitted to a subsample knows the points and not the layer's
  # distinct locations: one short of the points is still a bound.
  sub <- p; attr(sub, "bounds")$n_sample <- 1500L
  expect_identical(said(sub, 9), nearest("9", 333L, 100L, "30", build("9")))
  expect_identical(
    said(sub, 1),
    nearest("1", 3000L, 100L, "30",
            paste0("That count is more than k-means can place on 3000 points (2999 cells), ",
                   "so no profile scores it.")))
})

test_that("an optimum at a finest size the profile held is on the layer's bound, and is told so", {
  pts <- .cs_stations()
  # Four sizes, of which three ask for more cells than the 100 stations hold.
  p <- suppressWarnings(.cs_profile(pts, obs_per_cell = c(60, 30, 9, 4.5)))
  expect_identical(p$levels, c(50L, 100L))
  expect_identical(attr(p, "bounds")$ceiling_from, "distinct locations")
  q <- p; q$loo <- c(2, 1)
  s <- select_resolution(q, criterion = "loo")
  expect_identical(s$best, 100L)
  expect_identical(s$edge, "the support ceiling (as many cells as the distinct locations allow)")
  w <- .cs_of(.cs_catch(get_voronoi_seeds(method = "kmeans", n = s, sample_points = pts)),
              "spatialkit_bound_chosen")
  expect_length(w, 1L)
  expect_identical(conditionMessage(w[[1L]]),
                   paste0("get_voronoi_seeds(): `n` is 100 cells, the optimum of `loo`, and that ",
                          "optimum is the support ceiling (as many cells as the distinct ",
                          "locations allow): the bound chose the count, not the criterion. Pass ",
                          "a number as `n` to set the count yourself."))
  # The other end is the coarsest size the caller offered.
  q$loo <- c(1, 2)
  expect_identical(select_resolution(q, criterion = "loo")$edge,
                   "the coarsest size offered (obs_per_cell)")
  # A size held one short of the points, where no location repeats.
  one <- suppressWarnings(.cs_profile(.cs_600(), obs_per_cell = c(4.5, 1)))
  expect_identical(one$levels, c(133L, 599L))
  one$loo <- c(2, 1)
  expect_identical(select_resolution(one, criterion = "loo")$edge,
                   "the support ceiling (one short of the points, the most k-means can fit)")
})

# ---- The record ----------------------------------------------------------------

test_that("the record has the fields of the design and no other", {
  pts <- .cs_600()
  p <- .cs_pooled()
  s <- select_resolution(p, criterion = "size")
  rec <- s$count_rule
  est <- attr(p, "estimate")
  expect_identical(names(rec), .cs_record)
  expect_identical(class(rec), "spatialkit_count_rule")
  expect_identical(
    unclass(rec),
    list(version = 1L, rule = "size", estimate = "pooled", cells = 133L,
         obs_per_cell = 600 / 133, obs_per_cell_asked = 4.5, response = "z",
         n = 600L, n_points = 600L, n_fit = 600L, seed = 123, features = "points",
         scored = list(estimate = "pooled", n = 600L, cells = unname(est$cells[["133"]]),
                       loo = p$loo, loo_se = p$loo_se, loo_one_cell = est$loo_one_cell,
                       share_alone = p$share_alone, own_weight = p$own_weight,
                       structured = s$structured, rows = "all"),
         partition = list(kind = "density", init = "kmeans++", density_power = 1,
                          density_k = 8L, nstart = attr(p, "nstart")),
         crs = sf::st_crs(pts), scored_cells = NA))
  expect_identical(rec$partition$nstart, 3L)
  expect_identical(rec$scored$cells, 133L)
  expect_true(s$structured)
  # Beside the scored centres the same list says whose centres they are.
  on <- attr(s$seeds, "count_rule")
  expect_identical(on$scored_cells, TRUE)
  on$scored_cells <- NA
  expect_identical(on, rec)
  expect_identical(unname(sf::st_coordinates(s$seeds)[, 1:2]), unname(attr(p, "centres")[["133"]]))
  # The other attributes of the seeds are those of today.
  expect_setequal(setdiff(names(attributes(s$seeds)), "count_rule"),
                  names(attributes(.profile_seeds(p, 133L))))
})

test_that("the record's row counts tell a subsample fit from a fit to every row", {
  # A ladder of sizes is fitted to every row.
  rec <- select_resolution(.cs_pooled7(), criterion = "loo")$count_rule
  expect_identical(c(rec$n, rec$n_points, rec$n_fit), c(600L, 600L, 600L))
  # The automatic ladder of 3,000 points (seed 4402) is fitted to 1,500.
  big <- .cs_layer(3000, 4402)
  p <- suppressWarnings(.cs_profile(big, n_levels = 5L, sample_n = 1500L, sac = 3e4,
                                    obs_per_cell = NULL))
  s <- select_resolution(p, criterion = "loo")
  rec <- s$count_rule
  expect_identical(c(rec$n, rec$n_points, rec$n_fit), c(3000L, 3000L, 1500L))
  expect_identical(rec$rule, "loo")
  expect_identical(rec$estimate, "pooled")
  expect_identical(rec$obs_per_cell, 3000 / s$best)
  expect_identical(rec$obs_per_cell_asked, NA_real_)
  expect_identical(rec$scored$n, 3000L)
  expect_identical(rec$scored$structured, s$structured)
  expect_identical(s$obs_per_cell, 3000 / s$best)
  expect_identical(s$scored_per_cell, p$scored_per_cell[match(s$best, p$levels)])
  # This ladder has a range floor above its first level, which the
  # reliability is read from.  The leave-one-out error is read at every
  # level, and no floor is held for it.
  b <- attr(p, "bounds")
  expect_true(b$supported && is.finite(b$range) && !b$levels_given)
  expect_gt(b$floor, min(p$levels))
  expect_identical(s$range_floor, NA_integer_)
  low <- p; low$loo <- as.numeric(seq_along(p$levels))
  s <- select_resolution(low, criterion = "loo", range_floor = TRUE)
  expect_identical(s$best, min(p$levels))
  expect_identical(s$range_floor, NA_integer_)
  expect_identical(s$edge, "the first level of the ladder")
  # Half of the rows scored: the record says which.
  half <- .cs_profile(.cs_600(), select_on = "split")
  s <- select_resolution(half, criterion = "size")
  expect_identical(s$count_rule$scored$rows, "selection half")
  expect_identical(s$count_rule$scored$n, attr(half, "estimate")$n)
  # The count chosen is the centres; the cells that hold a scored row are
  # fewer, and are the record's of what was scored.
  expect_identical(s$count_rule$cells, 133L)
  expect_identical(s$count_rule$scored$cells, unname(attr(half, "estimate")$cells[["133"]]))
  expect_lt(s$count_rule$scored$cells, 133L)
  expect_lt(s$count_rule$scored$n, 600L)
  expect_identical(s$count_rule$n, 600L)
  expect_identical(s$obs_per_cell, 600 / 133)
  expect_identical(s$scored_per_cell, half$scored_per_cell)
  expect_false(isTRUE(all.equal(s$scored_per_cell, s$obs_per_cell)))
})

test_that("the record names the estimate its rule was computed for", {
  p <- .cs_levels()
  expect_identical(attr(p, "estimate")$value, "pooled")
  # Cp scores plain means, whatever the profile's estimate; the leave-one-out
  # error and a size are the profile's.
  cp <- select_resolution(p, criterion = "cp")
  expect_identical(cp$count_rule$rule, "cp")
  expect_identical(cp$count_rule$estimate, "mean")
  expect_identical(cp$count_rule$scored$estimate, "pooled")
  expect_identical(cp$estimate, "pooled")
  expect_identical(cp$obs_per_cell_asked, NA_real_)
  expect_identical(names(cp), c(.cs_present, .cs_added, "partition"))
  expect_identical(attr(cp$seeds, "count_rule")$scored_cells, TRUE)
  # What the profile scored at Cp's level, for its own estimate.
  expect_identical(cp$count_rule$scored$loo, p$loo[match(cp$best, p$levels)])
  expect_identical(cp$count_rule$scored$structured,
                   .estimate_structure(p, cp$best, 0.02, n_searched = 1L)$structured)
  # It is the estimate's verdict, not Cp's: with Cp's own F test made to
  # fail, `structured` of the selection is FALSE and the record's is as it
  # was.
  expect_true(cp$structured)
  expect_true(cp$count_rule$scored$structured)
  chance <- p; attr(chance, "one_cell")$f_p[] <- 0.5
  s <- select_resolution(chance, criterion = "cp")
  expect_false(s$structured)
  expect_true(s$count_rule$scored$structured)
  # And it follows the estimate's own F test, at one level searched.
  weak <- p; attr(weak, "estimate")$f_p[] <- 0.04
  expect_true(select_resolution(weak, criterion = "cp")$count_rule$scored$structured)
  weak <- p; attr(weak, "estimate")$f_p[] <- 0.06
  expect_false(select_resolution(weak, criterion = "cp")$count_rule$scored$structured)
  # Cp's own answer is the one a profile without the estimate gives.
  bare <- suppressWarnings(.cs_profile(.cs_600(), levels = c(50, 133, 300), sac = 3e4,
                                       estimate = NULL, obs_per_cell = NULL))
  expect_identical(cp[setdiff(.cs_present, c("note", "seeds"))],
                   select_resolution(bare, criterion = "cp")[setdiff(.cs_present, c("note", "seeds"))])
  expect_identical(.count_rule_new(p, "loo", 133L)$estimate, "pooled")
  expect_identical(.count_rule_new(p, "size", 133L, 4.5)$estimate, "pooled")
  for (rule in c("reliability", "elbow", "moran_z"))
    expect_identical(.count_rule_new(p, rule, 133L)$estimate, NA_character_, info = rule)
  mz <- select_resolution(p, criterion = "moran_z", diagnostic = TRUE)
  expect_identical(mz$count_rule$rule, "moran_z")
  expect_identical(mz$count_rule$estimate, NA_character_)
  expect_identical(.count_rule_new(.cs_mean7(), "loo", 133L)$estimate, "mean")
  # Without an estimate a size is computed for none, and nothing is known of
  # the features the profile read.
  expect_identical(.count_rule_new(.cs_plain(), "size", 50L, 4.5)$estimate, NA_character_)
  expect_identical(.count_rule_new(.cs_plain(), "size", 50L, 4.5)$features, NA_character_)
  # Features that are not points are recorded as such: 200 discs.
  discs <- sf::st_buffer(.cs_600()[1:200, ], 50)
  d <- suppressWarnings(.cs_profile(discs))
  expect_identical(attr(d, "estimate")$features, "other")
  expect_identical(select_resolution(d, criterion = "size")$count_rule$features, "other")
  # A profile that keeps no seed.
  expect_identical(.count_rule_new(.cs_plain(), "size", 50L, 4.5)$seed, NA_integer_)
  unseeded <- .cs_profile(.cs_600(), seed = NULL)
  expect_identical(select_resolution(unseeded, criterion = "size")$count_rule$seed, NA_integer_)
})

test_that("the record of plain means on a layer where some rows lack a response", {
  some <- .cs_600(); some$z[1:200] <- NA
  m <- .cs_profile(some, estimate = "mean", obs_per_cell = 4)
  s <- select_resolution(m, criterion = "size")
  rec <- s$count_rule
  expect_identical(s$best, 100L)
  # The size divides the rows with a response, not the points.
  expect_identical(c(rec$n, rec$n_points, rec$n_fit), c(400L, 600L, 600L))
  expect_identical(rec$obs_per_cell, 4)
  expect_identical(s$obs_per_cell, 4)
  expect_identical(rec$obs_per_cell_asked, 4)
  # The estimate is the profile's, in the rule and in what was scored.
  expect_identical(rec$estimate, "mean")
  expect_identical(rec$scored$estimate, "mean")
  expect_identical(rec$scored$n, 400L)
  expect_identical(rec$scored$loo, m$loo)
  expect_identical(rec$scored$own_weight, m$own_weight)
  expect_identical(rec$response, "z")
})

test_that("the pieces the record is read from", {
  b <- list(n = 600L, n_response = 400L)
  expect_identical(.count_rule_rows(b), 400L)
  expect_identical(.count_rule_rows(list(n = 600L, n_response = NA_integer_)), 600L)
  expect_identical(.count_rule_rows(list(n = 600L, n_response = 0L)), 600L)
  expect_identical(.count_rule_rows(NULL), NA_integer_)
  expect_identical(.count_rule_rows(list()), NA_integer_)
  p <- .cs_pooled7()
  expect_identical(.profile_at(p, 133L, "loo"), p$loo[5L])
  expect_identical(.profile_at(p, 134L, "loo"), NA_real_)
  expect_identical(.profile_at(p, 133L, "nothing"), NA_real_)
  # The structure test reads the profile's own numbers, and nothing once the
  # column is another.
  est <- attr(p, "estimate")
  st <- .estimate_structure(p, 133L, 0.02, n_searched = 7L)
  want <- .no_structure(p$loo[5L], est$loo_one_cell, est$f_p[["133"]], 0.02, n_searched = 7L)
  expect_identical(st[c("structured", "structure_p", "beats")],
                   want[c("structured", "structure_p", "beats")])
  expect_identical(st$loo, p$loo[5L])
  expect_identical(st$one_cell, est$loo_one_cell)
  one <- .estimate_structure(p, 133L, 0.02, n_searched = 1L)
  expect_identical(one$structure_p, min(1, est$f_p[["133"]]))
  expect_identical(st$structure_p, min(1, 7 * est$f_p[["133"]]))
  q <- p; q$loo[1L] <- q$loo[1L] + 1
  expect_identical(.estimate_structure(q, 133L, 0.02)$structured, NA)
  # Without the error of one mean there is nothing to compare with.
  none <- list(structured = NA, structure_p = NA_real_, beats = NA, loo = NA_real_,
               one_cell = NA_real_)
  q <- p; attr(q, "estimate")$loo_one_cell <- NA_real_
  expect_identical(.estimate_structure(q, 133L, 0.02), none)
  expect_identical(.estimate_structure(.cs_plain(), 50L, 0.02), none)
  expect_identical(.estimate_structure(p[, c("levels", "loo")], 133L, 0.02), none)
  expect_identical(.estimate_structure(.cs_plain(), 50L, 0.02)$structured, NA)
})


# ---- The arguments ---------------------------------------------------------------

test_that("obs_per_cell is NULL or one number of at least 1, and is read by \"size\" alone", {
  p <- .cs_pooled()
  refused <- list(list(NA, "NA"), list(0.5, "0.5"), list(c(3, 4.5), "c(3, 4.5)"),
                  list("auto", "\"auto\""), list(numeric(0), "numeric(0)"), list(Inf, "Inf"),
                  list(TRUE, "TRUE"))
  for (v in refused) {
    e <- .cs_error(select_resolution(p, criterion = "size", obs_per_cell = v[[1L]]))
    expect_identical(conditionMessage(e),
                     sprintf(paste0("select_resolution(): `obs_per_cell` must be NULL or a single ",
                                    "number of at least 1; got %s."), v[[2L]]))
    expect_null(conditionCall(e))
  }
  # Refused whatever the criterion.
  expect_error(select_resolution(p, criterion = "loo", obs_per_cell = 0.5),
               "must be NULL or a single number of at least 1; got 0.5.", fixed = TRUE)
  # A size beside a criterion that is not "size".
  l <- .cs_levels()
  for (cr in c("cp", "reliability", "elbow", "moran_z", "loo")) {
    e <- .cs_error(select_resolution(l, criterion = cr, obs_per_cell = 4.5))
    expect_identical(conditionMessage(e),
                     sprintf(paste0("select_resolution(): `obs_per_cell` is read by criterion = ",
                                    "\"size\" only, and criterion = \"%s\" was asked for."), cr))
  }
  # "auto" that does not resolve to "size" says what it resolved to.
  e <- .cs_error(select_resolution(.cs_mean7(), criterion = "auto", obs_per_cell = 4.5))
  expect_identical(conditionMessage(e),
                   paste0("select_resolution(): `obs_per_cell` is read by criterion = \"size\" ",
                          "only, and criterion = \"auto\" was asked for, which reads this ",
                          "profile by \"loo\"."))
  e <- .cs_error(select_resolution(.cs_plain(), criterion = "auto", obs_per_cell = 4.5))
  expect_match(conditionMessage(e), "which reads this profile by \"cp\".", fixed = TRUE)
  # And one that does is read.
  expect_identical(select_resolution(.cs_pooled7(), criterion = "auto", obs_per_cell = 9)$best,
                   rung_k(600L, 9))
  # 1 is a size; a whole number and a named one are read as numbers.
  expect_identical(suppressWarnings(
    select_resolution(p, criterion = "size", obs_per_cell = 1))$obs_per_cell_asked, 1)
  expect_identical(suppressWarnings(
    select_resolution(p, criterion = "size", obs_per_cell = c(a = 9L)))$obs_per_cell_asked, 9)
  expect_error(select_resolution(p, criterion = "sizes"), "should be one of")
})


test_that("what the selection reads before it applies a rule: the size, the subset, the estimate", {
  p <- .cs_pooled()
  # A number that carries a class is not a size.
  e <- .cs_error(select_resolution(p, criterion = "size",
                                   obs_per_cell = structure(4.5, class = "units")))
  expect_identical(conditionMessage(e),
                   paste0("select_resolution(): `obs_per_cell` must be NULL or a single number of ",
                          "at least 1; got an object of class units."))
  # An empty subset is told so under every criterion, a size included.
  for (cr in c("size", "loo", "cp", "auto"))
    expect_identical(conditionMessage(.cs_error(select_resolution(p[0, ], criterion = cr))),
                     paste0("select_resolution(): `profile` has no levels; it is an empty subset ",
                            "of a resolution profile."), info = cr)
  # The estimate is read by value.
  expect_identical(.profile_estimate(p), "pooled")
  expect_identical(.profile_estimate(.cs_mean7()), "mean")
  expect_identical(.profile_estimate(.cs_plain()), NA_character_)
  expect_identical(.profile_estimate(p[, c("levels", "loo")]), NA_character_)
  # An "estimate" attribute that names neither estimate is no estimate:
  # "auto" reads the profile by cp, and the selection is the one cp gives.
  q <- .cs_plain(); attr(q, "estimate") <- list(value = "median")
  expect_identical(.profile_estimate(q), NA_character_)
  s <- select_resolution(q, criterion = "auto")
  expect_identical(s$criterion, "cp")
  expect_identical(s, select_resolution(q, criterion = "cp"))
  expect_identical(names(s), .cs_present)
  # A profile built as today whose column is NA at every level has the
  # present error, not the one of a ladder of sizes.
  na <- .cs_plain(); na$cp <- NA_real_; na$reliability <- NA_real_
  expect_identical(conditionMessage(.cs_error(select_resolution(na, criterion = "cp"))),
                   paste0("select_resolution(): `cp` is NA at every level (it needs a response, ",
                          "and a level whose cells hold two scored rows each on average for its ",
                          "noise variance when the variogram gives no nugget)."))
  expect_identical(conditionMessage(.cs_error(select_resolution(na, criterion = "reliability"))),
                   paste0("select_resolution(): `reliability` is NA at every level (it needs a ",
                          "usable variogram with an identified range)."))
})


# ---- Structure under a size -----------------------------------------------------

test_that("a size on white noise is called structureless, and the seeding says so in the size's words", {
  skip_on_cran()
  # The layers of the engine's test 15: 300 locations uniform on a 1e5
  # square, x then y, then the response; white noise under seeds 4001 to
  # 4040, the field under 5001 to 5040.
  one <- function(seed, field) {
    pts <- .cs_layer(300, seed, field = field)
    p <- .cs_profile(pts)
    s <- select_resolution(p, criterion = "size")
    r <- .cs_catch(get_voronoi_seeds(method = "kmeans", n = s, sample_points = pts))
    tess <- build_tessellation(r$value, method = "voronoi", quiet = TRUE)
    asg  <- assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id")
    out  <- suppressWarnings(summarize_by_cell(asg, "z", id_col = "cell_id", estimate = "pooled"))
    list(sel = s, r = r, rows = nrow(asg), summary = attr(out, "pooled")$columns$resp_z)
  }
  noise <- lapply(1:40, function(r) one(4000 + r, FALSE))
  field <- lapply(1:40, function(r) one(5000 + r, TRUE))
  all80 <- c(noise, field)
  expect_true(all(vapply(all80, function(x) x$sel$best, 1L) == 67L))
  expect_false(anyNA(vapply(all80, function(x) x$sel$structured, NA)))
  # One run of these seeds: FALSE on 39 of the noise layers, TRUE on 40 of the
  # field layers.
  unstructured <- vapply(noise, function(x) isFALSE(x$sel$structured), NA)
  expect_gte(sum(unstructured), 36L)
  expect_gte(sum(vapply(field, function(x) isTRUE(x$sel$structured), NA)), 36L)
  # The warning is raised exactly where `structured` is FALSE, once.
  warned <- function(x) length(.cs_of(x$r, "spatialkit_no_structure"))
  expect_identical(vapply(all80, warned, 1L),
                   as.integer(vapply(all80, function(x) isFALSE(x$sel$structured), NA)))
  expect_false(any(vapply(all80, function(x) "spatialkit_bound_chosen" %in% .cs_classes(x$r), NA)))
  # In the words of a size: the count stands.
  for (x in noise[unstructured]) {
    w <- .cs_of(x$r, "spatialkit_no_structure")[[1L]]
    expect_identical(
      conditionMessage(w),
      paste0("get_voronoi_seeds(): `n` is 67 cells, set by obs_per_cell = 4.5, and ",
             x$sel$note[["structure"]],
             " The count stands, since a size set it and no criterion; the cell values will ",
             "say little that one mean for the layer does not."))
    expect_identical(.cs_fields(w), list(criterion = "size", n = 67L, one_cell = x$sel$one_cell))
    expect_identical(nrow(x$r$value), 67L)
    expect_identical(names(x$sel$note), c("size", "structure"))
  }
  # The two forms of the note: one mean does as well, or the F test does not
  # tell the cells from one mean.
  ends <- ": the response shows no spatial structure that cells of this size can use."
  for (x in noise[unstructured]) {
    s <- x$sel
    nt <- s$note[["structure"]]
    if (s$one_cell <= s$loo * 1.02)
      expect_identical(nt, sprintf(paste0("at 67 cells the leave-one-out error of the pooled ",
                                          "value (%s) is not below that of one mean for the ",
                                          "whole layer (%s) by more than 2%%%s"),
                                   format(signif(s$loo, 3)), format(signif(s$one_cell, 3)), ends))
    else
      expect_identical(nt, sprintf(paste0("the cell means at 67 cells differ from one mean for ",
                                          "the whole layer by no more than chance allows (F ",
                                          "test, p = %s)%s"),
                                   format(signif(s$structure_p, 2)), ends))
  }
  # What the summary of the same cells and rows records.
  expect_true(all(vapply(all80, function(x) x$rows, 1L) == 300L))
  for (x in all80) {
    expect_identical(x$summary$structured, x$sel$structured)
    expect_equal(x$summary$loo, x$sel$loo, tolerance = 1e-10)
    expect_equal(x$summary$loo_one_cell, x$sel$loo_one_cell, tolerance = 1e-10)
  }
})

test_that("a size without structure is said at the seeding, in every function that takes the count", {
  # White noise on 300 rows (seed 4001).
  pts <- .cs_layer(300, 4001, field = FALSE)
  bnd <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(pts)))
  p <- .cs_profile(pts)
  s <- select_resolution(p, criterion = "size")
  expect_false(s$structured)
  expect_identical(names(s$note), c("size", "structure"))
  text <- function(fn, arg)
    paste0(fn, "(): `", arg, "` is 67 cells, set by obs_per_cell = 4.5, and ",
           s$note[["structure"]],
           " The count stands, since a size set it and no criterion; the cell values will ",
           "say little that one mean for the layer does not.")
  expect_match(s$note[["structure"]],
               "^at 67 cells the leave-one-out error of the pooled value \\(")
  calls <- list(
    list("get_voronoi_seeds", "n",
         function(n) get_voronoi_seeds(method = "kmeans", n = n, sample_points = pts)),
    list("voronoi_seeds_kmeans", "k", function(n) voronoi_seeds_kmeans(pts, n)),
    list("build_tessellation", "approx_n_cells",
         function(n) build_tessellation(pts, boundary = bnd, method = "hex", approx_n_cells = n,
                                        quiet = TRUE)))
  for (cl in calls) for (count in list(s, p)) {
    r <- .cs_catch(cl[[3L]](count))
    w <- .cs_of(r, "spatialkit_no_structure")
    expect_length(w, 1L)
    expect_identical(conditionMessage(w[[1L]]), text(cl[[1L]], cl[[2L]]))
    expect_identical(.cs_fields(w[[1L]]), list(criterion = "size", n = 67L, one_cell = s$one_cell))
    expect_null(conditionCall(w[[1L]]))
    expect_false("spatialkit_bound_chosen" %in% .cs_classes(r))
  }
})

test_that("under a size `tol` is the margin of the structure test and one level is searched", {
  p <- .cs_pooled()
  est <- attr(p, "estimate")
  s <- select_resolution(p, criterion = "size")
  expect_true(s$structured)
  expect_identical(s$structure_p, min(1, est$f_p[["133"]]))
  expect_identical(s$one_cell, est$loo_one_cell)
  # One level, on a ladder of seven: the F test is not adjusted for the
  # levels a size did not search.
  p7 <- .cs_pooled7()
  s7 <- select_resolution(p7, criterion = "size")
  expect_identical(s7$structure_p, min(1, attr(p7, "estimate")$f_p[["133"]]))
  weak <- p7; attr(weak, "estimate")$f_p[] <- 0.02
  s7 <- select_resolution(weak, criterion = "size")
  expect_identical(s7$structure_p, 0.02)
  expect_true(s7$structured)
  # Where "loo" searched the seven, it is.
  expect_identical(select_resolution(weak, criterion = "loo")$structure_p, 0.14)
  expect_false(select_resolution(weak, criterion = "loo")$structured)
  # The F test alone can say no: the second form of the note.
  weak <- p7; attr(weak, "estimate")$f_p[] <- 0.2
  s7 <- select_resolution(weak, criterion = "size")
  expect_false(s7$structured)
  expect_identical(s7$note[["structure"]],
                   paste0("the cell means at 133 cells differ from one mean for the whole ",
                          "layer by no more than chance allows (F test, p = 0.2): the response ",
                          "shows no spatial structure that cells of this size can use."))
  # A margin the cells do not clear.
  wide <- est$loo_one_cell / p$loo - 1 + 0.01
  t <- select_resolution(p, criterion = "size", tol = wide)
  expect_false(t$structured)
  expect_identical(t$best, s$best)
  expect_identical(t$flat, s$flat)
  expect_match(t$note[["structure"]], sprintf("by more than %s%%", format(100 * wide)), fixed = TRUE)
  expect_true(select_resolution(p, criterion = "size", tol = wide - 0.02)$structured)
  # One mean's error equal to the cells' is not below it: the plain mean is
  # named for a profile of plain means.
  m <- .cs_profile(.cs_600(), estimate = "mean", obs_per_cell = 4.5)
  t <- select_resolution(m, criterion = "size", tol = 10)
  expect_false(t$structured)
  expect_match(t$note[["structure"]], "the leave-one-out error of the plain mean (", fixed = TRUE)
  # A constant response has a count and no verdict.
  pts <- .cs_600(); pts$k <- 0.1
  k <- suppressWarnings(.cs_profile(pts, response = "k"))
  s <- select_resolution(k, criterion = "size")
  expect_identical(s$best, 133L)
  expect_identical(s$structured, NA)
  expect_identical(names(s$note), "size")
})


# ---- "auto" ---------------------------------------------------------------------

test_that("\"auto\" is the rule that goes with the profile's estimate", {
  pooled <- .cs_pooled7(); plain <- .cs_mean7(); none <- .cs_plain()
  expect_identical(select_resolution(pooled, criterion = "auto"),
                   select_resolution(pooled, criterion = "size"))
  expect_identical(select_resolution(plain, criterion = "auto"),
                   select_resolution(plain, criterion = "loo"))
  expect_identical(select_resolution(none, criterion = "auto"),
                   select_resolution(none, criterion = "cp"))
  # A column subset has lost its estimate and is read by cp.
  l <- .cs_levels()
  expect_identical(select_resolution(l[, c("levels", "cp", "loo")], criterion = "auto"),
                   select_resolution(l[, c("levels", "cp", "loo")], criterion = "cp"))
  # `tol` goes with it.
  expect_identical(select_resolution(plain, criterion = "auto", tol = 0.2),
                   select_resolution(plain, criterion = "loo", tol = 0.2))
  # The criterion in the result is the one applied, and the record's rule.
  for (s in list(select_resolution(pooled, criterion = "auto"),
                 select_resolution(plain, criterion = "auto"))) {
    expect_true(s$criterion %in% c("size", "loo"))
    expect_identical(s$criterion, s$count_rule$rule)
    expect_identical(s$range_floor, NA_integer_)
    # On a density profile the kind is the last element, after the record.
    expect_identical(utils::tail(names(s), 2L), c("count_rule", "partition"))
    expect_identical(s$partition, list(kind = "density", init = "kmeans++", density_power = 1,
                                       density_k = 8L))
  }
  expect_identical(select_resolution(none, criterion = "auto")$criterion, "cp")
  expect_identical(select_resolution(pooled, criterion = "auto")$criterion, "size")
  expect_identical(select_resolution(plain, criterion = "auto")$criterion, "loo")
})


# ---- Structure under "loo" --------------------------------------------------------

test_that("\"loo\" calls white noise structureless and a field with a nugget structured", {
  skip_on_cran()
  # 300 locations uniform on a 1000 x 1000 square, x then y.  White noise: a
  # standard normal response, seeds 2001 to 2040.  A field: seeds 3001 to
  # 3020, an exponential covariance with range parameter 100, sill 1 and
  # nugget 0.3, drawn as the profile's test of new sites draws its own.
  one <- function(seed, field, estimate) {
    d <- local({
      .cs_seed(seed)
      x <- stats::runif(300, 0, 1000); y <- stats::runif(300, 0, 1000)
      z <- if (field) {
        D <- as.matrix(stats::dist(cbind(x, y)))
        as.numeric(t(chol(exp(-D / 100) + diag(1e-10, 300))) %*% stats::rnorm(300) +
                     stats::rnorm(300, 0, sqrt(0.3)))
      } else stats::rnorm(300)
      list(x = x, y = y, z = z)
    })
    p <- suppressWarnings(.cs_profile(.cs_points(d$x, d$y, d$z), estimate = estimate,
                                      obs_per_cell = .cs_sizes))
    select_resolution(p, criterion = "loo", tol = 0.02)$structured
  }
  for (estimate in c("mean", "pooled")) {
    noise <- vapply(1:40, function(r) one(2000 + r, FALSE, estimate), NA)
    field <- vapply(1:20, function(r) one(3000 + r, TRUE, estimate), NA)
    expect_false(anyNA(noise), info = estimate)
    expect_false(anyNA(field), info = estimate)
    # One run of these seeds: TRUE on 0 (plain means) and 1 (pooled values)
    # of the noise layers, and on 20 of the field layers for each.
    expect_lte(sum(noise), 6L)
    expect_gte(sum(field), 16L)
  }
})

test_that("the notes of \"loo\" on a response without structure are Cp's, for the leave-one-out error", {
  # White noise on 300 rows (seed 2003), the seven sizes, plain means.
  pts <- .cs_layer(300, 2003, side = 1000, field = FALSE)
  p <- suppressWarnings(.cs_profile(pts, estimate = "mean", obs_per_cell = .cs_sizes))
  est <- attr(p, "estimate")
  s <- select_resolution(p, criterion = "loo")
  expect_false(s$structured)
  expect_identical(s$one_cell, est$loo_one_cell)
  finite <- sum(is.finite(p$loo))
  expect_identical(s$structure_p, min(1, finite * est$f_p[[as.character(s$best)]]))
  # One mean does as well as any level.
  expect_lte(s$one_cell, s$value * 1.02)
  expect_identical(
    s$note,
    c(structure = sprintf(
      paste0("the leave-one-out error at one cell (%s) is within 2%% of the optimum (%s): one ",
             "mean for the whole layer predicts as well as any level scored, so the response ",
             "shows no spatial structure the cells can use and the count is not identified."),
      format(signif(s$one_cell, 3)), format(signif(s$value, 3)))))
  # The other ground: the cells beat one mean by the margin, and by no more
  # than chance.  A margin of 0 on a column made to beat one mean barely.
  q <- p
  q$loo <- est$loo_one_cell * c(0.999, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1)
  attr(q, "estimate")$loo_levels <- stats::setNames(q$loo, q$levels)
  attr(q, "estimate")$f_p[] <- 0.5
  t <- select_resolution(q, criterion = "loo", tol = 0)
  expect_false(t$structured)
  expect_identical(t$structure_p, 1)
  expect_identical(
    t$note[["structure"]],
    sprintf(paste0("the cell means at %d cells differ from one mean for the whole layer by no ",
                   "more than chance allows (F test of the cells against one mean, p = 1 after ",
                   "adjusting for the 7 levels searched), so the response shows no spatial ",
                   "structure the cells can use and the count is not identified. The ",
                   "leave-one-out error at one cell is %s against %s at the optimum, a gap that ",
                   "chance produces on a layer of this size."),
            q$levels[1L], format(signif(est$loo_one_cell, 3)), format(signif(q$loo[1L], 3))))
  # The adjustment is for the levels with a finite error.
  attr(q, "estimate")$f_p[] <- 0.004
  expect_true(select_resolution(q, criterion = "loo", tol = 0)$structured)
  attr(q, "estimate")$f_p[] <- 0.008
  expect_false(select_resolution(q, criterion = "loo", tol = 0)$structured)
  # No F test at the level: not known.
  attr(q, "estimate")$f_p[] <- NA_real_
  t <- select_resolution(q, criterion = "loo", tol = 0)
  expect_identical(t$structured, NA)
  expect_identical(
    t$note[["structure"]],
    paste0("the leave-one-out error at one cell is above the optimum by more than `tol`, but ",
           "the profile carries no F test of the cells against one mean at this level, so ",
           "whether that is more than chance is not known."))
  # Passed on as the count: the present warning, with `loo` as the criterion.
  r <- .cs_catch(get_voronoi_seeds(method = "kmeans", n = s, sample_points = pts))
  w <- .cs_of(r, "spatialkit_no_structure")
  expect_length(w, 1L)
  expect_identical(conditionMessage(w[[1L]]),
                   sprintf(paste0("get_voronoi_seeds(): `n` is %d cells, the optimum of `loo`, ",
                                  "but %s Pass a number as `n` if a count is wanted on other ",
                                  "grounds."), s$best, s$note[["structure"]]))
  expect_identical(.cs_fields(w[[1L]]),
                   list(criterion = "loo", n = s$best, one_cell = s$one_cell))
})

test_that("\"loo\" adjusts its F test for the levels that have an error, and the record holds its verdict", {
  p <- .cs_pooled7()
  full <- select_resolution(p, criterion = "loo")
  expect_true(full$structured)
  # Seven levels searched: 0.0075 at the level is 0.0525 after the
  # adjustment, and the count is not identified.
  weak <- p; attr(weak, "estimate")$f_p[] <- 0.0075
  s <- select_resolution(weak, criterion = "loo")
  expect_identical(s$structure_p, 0.0075 * 7)
  expect_false(s$structured)
  # A level without an error was not searched.  Six levels: 0.045.
  drop <- which(p$levels != full$best)[1L]
  six <- weak; six$loo[drop] <- NA_real_
  attr(six, "estimate")$loo_levels[drop] <- NA_real_
  s <- select_resolution(six, criterion = "loo")
  expect_identical(s$best, full$best)
  expect_identical(s$structure_p, 0.0075 * 6)
  expect_true(s$structured)
  # The record's verdict is the selection's, adjusted as it is: not the one
  # the level would get had it been the only level searched.
  two <- p; attr(two, "estimate")$f_p[] <- 0.02
  s <- select_resolution(two, criterion = "loo")
  expect_identical(s$structure_p, 0.14)
  expect_false(s$structured)
  expect_false(s$count_rule$scored$structured)
  expect_false(attr(s$seeds, "count_rule")$scored$structured)
  expect_true(.estimate_structure(two, s$best, 0.02, n_searched = 1L)$structured)
})

test_that("a selection whose notes the caller wrote is told in the caller's words at the seeding", {
  pts <- .cs_noise()
  s <- select_resolution(.cs_noise_plain(), criterion = "cp")
  expect_false(s$structured)
  expect_null(names(s$note))
  warned <- function(sel) {
    r <- .cs_catch(get_voronoi_seeds(method = "kmeans", n = sel, sample_points = pts))
    conditionMessage(.cs_of(r, "spatialkit_no_structure")[[1L]])
  }
  text <- function(why)
    sprintf(paste0("get_voronoi_seeds(): `n` is %d cells, the optimum of `cp`, but %s Pass a ",
                   "number as `n` if a count is wanted on other grounds."), s$best, why)
  stock <- "one cell for the whole layer does as well, so the count is not identified."
  # As select_resolution() returned it: every note.
  expect_identical(warned(s), text(paste(s$note, collapse = " ")))
  expect_match(s$note, "^Cp at one cell \\(")
  # What the caller put there is pasted as it always was: a note with a name
  # of the caller's, a note beside the selection's own, an emptied vector
  # (pasted as nothing), and no note at all (the stock sentence).
  named <- s; named$note <- c(remark = "the analyst's own remark.")
  expect_identical(warned(named), text("the analyst's own remark."))
  beside <- s; beside$note <- c(s$note, remark = "the analyst's own remark.")
  expect_identical(warned(beside), text(paste(s$note, "the analyst's own remark.")))
  empty <- s; empty$note <- character(0)
  expect_identical(warned(empty), text(""))
  none <- s; none$note <- NULL
  expect_identical(warned(none), text(stock))
  # The provenance of a selection is the line it always was, whatever the
  # selection holds as its criterion.
  seeds <- suppressWarnings(get_voronoi_seeds(method = "kmeans", n = s, sample_points = pts))
  expect_identical(attr(seeds, "n_from"),
                   sprintf("select_resolution(criterion = \"cp\"), an optimum at %s", s$edge))
  bare <- s; bare$criterion <- NULL
  expect_identical(suppressWarnings(.resolve_cell_count(bare, "n", "f"))$from, character(0))
  # Cp without structure on a profile of pooled values: Cp's own sentence has
  # no name and is said; the note about the estimate has one and is not.
  pooled <- suppressWarnings(.cs_profile(pts, levels = c(10, 30, 67), sac = 3e4,
                                         obs_per_cell = NULL))
  sp <- select_resolution(pooled, criterion = "cp")
  expect_false(sp$structured)
  expect_identical(names(sp$note), c("", "estimate"))
  r <- .cs_catch(get_voronoi_seeds(method = "kmeans", n = sp, sample_points = pts))
  expect_identical(conditionMessage(.cs_of(r, "spatialkit_no_structure")[[1L]]),
                   sprintf(paste0("get_voronoi_seeds(): `n` is %d cells, the optimum of `cp`, but ",
                                  "%s Pass a number as `n` if a count is wanted on other grounds."),
                           sp$best, sp$note[[1L]]))
  # The same on a ladder of one size, whose note has a name too.
  one <- suppressWarnings(.cs_profile(pts, sac = 3e4))
  so <- select_resolution(one, criterion = "cp")
  expect_false(so$structured)
  expect_identical(names(so$note), c("", "estimate", "ladder"))
  r <- .cs_catch(get_voronoi_seeds(method = "kmeans", n = so, sample_points = pts))
  expect_identical(conditionMessage(.cs_of(r, "spatialkit_no_structure")[[1L]]),
                   sprintf(paste0("get_voronoi_seeds(): `n` is %d cells, the optimum of `cp`, but ",
                                  "%s Pass a number as `n` if a count is wanted on other grounds."),
                           so$best, so$note[[1L]]))
  # A selection the caller relabelled is not one a size set: it carries no
  # size asked for, and is told in the present words.
  relabelled <- s; relabelled$criterion <- "size"
  r <- .cs_catch(get_voronoi_seeds(method = "kmeans", n = relabelled, sample_points = pts))
  expect_identical(conditionMessage(.cs_of(r, "spatialkit_no_structure")[[1L]]),
                   sprintf(paste0("get_voronoi_seeds(): `n` is %d cells, the optimum of `size`, ",
                                  "but %s Pass a number as `n` if a count is wanted on other ",
                                  "grounds."), s$best, paste(s$note, collapse = " ")))
  expect_identical(attr(r$value, "n_from"),
                   sprintf("select_resolution(criterion = \"size\"), an optimum at %s", s$edge))
  # The five names the warning leaves out, and no other.
  expect_identical(.NOTES_NOT_STRUCTURE, c("size", "held", "finest", "estimate", "ladder"))
  # A size without structure whose note of it was taken away: the size's
  # stock sentence, and nothing of the note that names the size.
  sized <- select_resolution(.cs_profile(pts), criterion = "size")
  expect_identical(names(sized$note), c("size", "structure"))
  sized$note <- sized$note["size"]
  r <- .cs_catch(get_voronoi_seeds(method = "kmeans", n = sized, sample_points = pts))
  expect_identical(conditionMessage(.cs_of(r, "spatialkit_no_structure")[[1L]]),
                   paste0("get_voronoi_seeds(): `n` is 67 cells, set by obs_per_cell = 4.5, and ",
                          "one mean for the whole layer predicts a left-out observation as well ",
                          "as cells of this size do. The count stands, since a size set it and ",
                          "no criterion; the cell values will say little that one mean for the ",
                          "layer does not."))
})

test_that("the warning about structure leaves out the notes that are about something else", {
  # White noise on 300 rows (seed 2003), pooled values, the seven sizes with
  # the column made to fall to the finest level.
  pts <- .cs_layer(300, 2003, side = 1000, field = FALSE)
  p <- suppressWarnings(.cs_profile(pts, obs_per_cell = .cs_sizes))
  est <- attr(p, "estimate")
  q <- p
  q$loo <- est$loo_one_cell * (1 + (7:1) / 1000)
  attr(q, "estimate")$loo_levels <- stats::setNames(q$loo, q$levels)
  s <- select_resolution(q, criterion = "loo")
  expect_false(s$structured)
  expect_identical(names(s$note), c("structure", "finest"))
  r <- .cs_catch(get_voronoi_seeds(method = "kmeans", n = s, sample_points = pts))
  w <- .cs_of(r, "spatialkit_no_structure")[[1L]]
  expect_identical(conditionMessage(w),
                   sprintf(paste0("get_voronoi_seeds(): `n` is %d cells, the optimum of `loo`, ",
                                  "but %s Pass a number as `n` if a count is wanted on other ",
                                  "grounds."), s$best, s$note[["structure"]]))
  expect_false(grepl("finest level scored", conditionMessage(w), fixed = TRUE))
})


# ---- cp and reliability on a ladder of sizes ---------------------------------------

test_that("cp and reliability on a ladder of sizes that fitted no variogram are an error that names the ways back", {
  p <- .cs_pooled()
  expect_identical(attr(p, "cp_noise")$source, "none")
  text <- function(criterion, from, sizes = "4.5")
    sprintf(paste0("select_resolution(): `%s` is NA at every level: this profile is on a size ",
                   "ladder (obs_per_cell = %s), which fits no variogram unless `sac` is given, ",
                   "and `%s` takes its %s from one. Build the profile with `sac = NULL` to have ",
                   "the variogram estimated, with a `sac` of your own, or with obs_per_cell = ",
                   "NULL for the automatic ladder; criterion = \"size\" reads this profile."),
            criterion, sizes, criterion, from)
  e <- .cs_error(select_resolution(p, criterion = "cp"))
  expect_identical(conditionMessage(e), text("cp", "noise variance"))
  expect_null(conditionCall(e))
  e <- .cs_error(select_resolution(p, criterion = "reliability"))
  expect_identical(conditionMessage(e), text("reliability", "range"))
  # Several sizes are named together.
  e <- .cs_error(select_resolution(.cs_mean7(), criterion = "cp"))
  expect_identical(conditionMessage(e),
                   text("cp", "noise variance", "c(18, 12, 9, 6, 4.5, 3, ...)"))
  # The two rules that need no variogram run.
  expect_identical(select_resolution(p, criterion = "size")$best, 133L)
  expect_identical(select_resolution(p, criterion = "loo")$best, 133L)
  # The other present criteria keep their own errors there.
  expect_error(select_resolution(p, criterion = "elbow"),
               "`elbow` is NA at every level (the WSS curve has no elbow", fixed = TRUE)
  # A constant response is said first: no variogram would change it.
  pts <- .cs_600(); pts$k <- 0.1
  k <- suppressWarnings(.cs_profile(pts, response = "k"))
  expect_error(select_resolution(k, criterion = "cp"),
               "`cp` is NA at every level because the response is constant", fixed = TRUE)
})

.cs_one_size_note <- function(criterion, sizes = "4.5")
  sprintf(paste0("the ladder is one size (obs_per_cell = %s), so %s had one level to choose ",
                 "from. Build the profile with obs_per_cell = NULL, or with several sizes, for ",
                 "a ladder %s can choose on."), sizes, criterion, criterion)
.cs_pooled_note <- paste0("cp scores plain means, and this profile's estimate is the pooled ",
                          "value; criterion = \"size\" is the rule for it.")

test_that("with sac = NULL a ladder of one size has a cp, which returns the level and says it had one", {
  skip_if_not_installed("gstat")
  # The variogram is estimated, as on any profile of today.
  p <- suppressWarnings(.cs_profile(.cs_600(), sac = NULL))
  expect_false(identical(attr(p, "cp_noise")$source, "none"))
  expect_true(is.finite(p$cp))
  s <- select_resolution(p, criterion = "cp")
  expect_identical(s$best, 133L)
  expect_identical(s$edge, "the only level scored")
  expect_identical(s$note[["ladder"]], .cs_one_size_note("cp"))
  expect_identical(s$note[["estimate"]], .cs_pooled_note)
  expect_identical(s$count_rule$rule, "cp")
  expect_identical(s$count_rule$estimate, "mean")
})

test_that("with a variogram cp on a ladder of one size returns the level, and says it had one", {
  pts <- .cs_600()
  note <- .cs_one_size_note
  # A `sac` of the caller's own.
  p <- suppressWarnings(.cs_profile(pts, sac = 3e4))
  expect_true(is.finite(p$cp))
  s <- select_resolution(p, criterion = "cp")
  expect_identical(s$best, 133L)
  expect_identical(s$edge, "the only level scored")
  expect_identical(s$note[c("estimate", "ladder")],
                   c(estimate = .cs_pooled_note, ladder = note("cp")))
  # Without an estimate the list is today's, with the note in it.
  bare <- suppressWarnings(.cs_profile(pts, sac = 3e4, estimate = NULL))
  s <- select_resolution(bare, criterion = "cp")
  expect_identical(names(s), c(.cs_present, "partition"))
  expect_identical(s$note[["ladder"]], note("cp"))
  expect_false("estimate" %in% names(s$note))
  expect_identical(s$edge, "the only level scored")
  # The reliability, where the profile has one, is told the same.
  r <- bare; r$reliability <- 0.5
  expect_identical(select_resolution(r, criterion = "reliability")$note,
                   c(ladder = note("reliability")))
  # The other criteria are not.
  expect_false("ladder" %in%
                 names(select_resolution(bare, criterion = "moran_z", diagnostic = TRUE)$note))
  # Several sizes are a ladder to choose on: no such note, and the edges are
  # the sizes offered.
  many <- suppressWarnings(.cs_profile(pts, sac = 3e4, obs_per_cell = c(9, 4.5, 3)))
  s <- select_resolution(many, criterion = "cp")
  expect_true("estimate" %in% names(s$note))
  expect_false("ladder" %in% names(s$note))
  expect_true(is.na(s$edge) || grepl("(obs_per_cell)", s$edge, fixed = TRUE))
  # Plain means are what cp scores: nothing is said of the estimate.
  m <- suppressWarnings(.cs_profile(pts, sac = 3e4, estimate = "mean", obs_per_cell = c(9, 4.5, 3)))
  expect_false("estimate" %in% names(select_resolution(m, criterion = "cp")$note))
})


test_that("the notes about the estimate and the numbers of the level are written where they apply", {
  # The reliability on a profile of pooled values: the note about a ladder
  # of one size, and nothing about the estimate cp scores.
  v <- suppressWarnings(.cs_profile(.cs_600(), sac = 3e4))
  v$reliability <- 0.5
  s <- select_resolution(v, criterion = "reliability")
  expect_identical(s$note, c(ladder = .cs_one_size_note("reliability")))
  l <- .cs_levels()
  expect_false("estimate" %in%
                 names(select_resolution(l, criterion = "moran_z", diagnostic = TRUE)$note))
  # A profile whose `loo` column is gone has no error to report, at the level
  # or for one cell.
  q <- l; q$loo <- NULL
  expect_identical(attr(q, "estimate")$value, "pooled")
  s <- select_resolution(q, criterion = "cp")
  expect_identical(s$estimate, "pooled")
  expect_identical(c(s$loo, s$loo_one_cell), c(NA_real_, NA_real_))
  expect_identical(s$count_rule$scored$structured, NA)
})


# ---- A profile as a count -----------------------------------------------------------

test_that("a profile with an estimate passed as a count is read by \"auto\"", {
  pts <- .cs_600()
  bnd <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(pts)))
  p <- .cs_pooled()
  # No criterion a profile is read by today is finite on it.
  expect_true(all(is.na(p$cp)) && all(is.na(p$reliability)) && all(is.na(p$elbow)))
  expect_identical(as.integer(round(600 / 4.5)), 133L)
  from <- "resolution_profile() read with select_resolution(criterion = \"size\", obs_per_cell = 4.5)"
  r <- .cs_catch(get_voronoi_seeds(method = "kmeans", n = p, sample_points = pts))
  expect_length(r$warnings, 0L)
  expect_identical(nrow(r$value), 133L)
  expect_identical(attr(r$value, "n_from"), from)
  expect_identical(unname(sf::st_coordinates(r$value)[, 1:2]), unname(attr(p, "centres")[["133"]]))
  r <- .cs_catch(voronoi_seeds_kmeans(pts, p))
  expect_length(r$warnings, 0L)
  expect_identical(nrow(r$value), 133L)
  expect_identical(attr(r$value, "n_from"), from)
  r <- .cs_catch(build_tessellation(pts, boundary = bnd, method = "hex", approx_n_cells = p,
                                    quiet = TRUE))
  expect_false(any(c("spatialkit_bound_chosen", "spatialkit_no_structure") %in% .cs_classes(r)))
  expect_identical(r$value$params$approx_n_cells, 133)
  expect_identical(r$value$params$approx_n_cells_from, from)
  # A profile of one size that has a cp is read by its size all the same.
  # Read by cp the level would be an optimum on a bound, and the count one
  # chosen for plain means.
  v <- suppressWarnings(.cs_profile(pts, sac = 3e4))
  expect_true(is.finite(v$cp))
  expect_identical(select_resolution(v, criterion = "cp")$edge, "the only level scored")
  r <- .cs_catch(get_voronoi_seeds(method = "kmeans", n = v, sample_points = pts))
  expect_length(r$warnings, 0L)
  expect_identical(nrow(r$value), 133L)
  expect_identical(attr(r$value, "n_from"), from)
  res <- .resolve_cell_count(v, "n", "get_voronoi_seeds")
  expect_identical(res$count_rule[c("rule", "estimate", "cells")],
                   list(rule = "size", estimate = "pooled", cells = 133L))
  expect_identical(res$edge, NA_character_)
  # The selection gives the same seeds, and names its size.
  s <- select_resolution(p, criterion = "size")
  viasel <- get_voronoi_seeds(method = "kmeans", n = s, sample_points = pts)
  expect_identical(sf::st_coordinates(viasel),
                   sf::st_coordinates(get_voronoi_seeds(method = "kmeans", n = p,
                                                        sample_points = pts)))
  expect_identical(attr(viasel, "n_from"),
                   "select_resolution(criterion = \"size\", obs_per_cell = 4.5)")
  # A profile of plain means is read by its leave-one-out error.
  m <- .cs_mean7()
  r <- .cs_catch(get_voronoi_seeds(method = "kmeans", n = m, sample_points = pts))
  loo <- select_resolution(m, criterion = "loo")
  expect_identical(nrow(r$value), loo$best)
  expect_match(attr(r$value, "n_from"),
               "resolution_profile() read with select_resolution(criterion = \"loo\")", fixed = TRUE)
  expect_identical("spatialkit_bound_chosen" %in% .cs_classes(r), !is.na(loo$edge))
  # The provenance names the bound an optimum sits on: one size is one level.
  m1 <- .cs_profile(pts, estimate = "mean", obs_per_cell = 4.5)
  r <- .cs_catch(get_voronoi_seeds(method = "kmeans", n = m1, sample_points = pts))
  expect_identical(attr(r$value, "n_from"),
                   paste0("resolution_profile() read with select_resolution(criterion = ",
                          "\"loo\"), an optimum at the only level scored"))
  expect_identical(.cs_classes(r), "spatialkit_bound_chosen")
})

test_that("a profile on a ladder of sizes without an estimate, passed as a count, is told what reads it", {
  pts <- .cs_600()
  text <- function(sizes, why = "", tail = "")
    sprintf(paste0("get_voronoi_seeds(): `n` is a resolution profile with no criterion finite ",
                   "at any level, so no cell count can be read off it. It is on a ladder of ",
                   "sizes (obs_per_cell = %s) and was built without an `estimate`, so it is read ",
                   "by `cp`, or by `elbow` when it has no `response_var`, and neither is finite ",
                   "on it%s. A size reads it: pass select_resolution(<profile>, criterion = ",
                   "\"size\") as `n`, or pass a number.%s"), sizes, why, tail)
  # The geometry alone at one size.  A size reads it; as a count it is read
  # as a profile without an estimate is, and has nothing such a profile is
  # read by.
  geo <- resolution_profile(pts, partition = "kmeans", nstart = 2, sac = "auto", estimate = NULL,
                            obs_per_cell = 6, seed = 123)
  expect_identical(select_resolution(geo, criterion = "size")$best, 100L)
  e <- .cs_error(get_voronoi_seeds(method = "kmeans", n = geo, sample_points = pts))
  expect_identical(conditionMessage(e), text("6"))
  expect_null(conditionCall(e))
  expect_false(inherits(e, "spatialkit_no_default_criterion"))
  # A response at two sizes, with no variogram fitted.
  two <- resolution_profile(pts, "z", partition = "kmeans", nstart = 2, sac = "auto",
                            estimate = NULL, obs_per_cell = c(9, 4.5), seed = 123)
  expect_true(all(is.na(two$cp)) && all(is.na(two$reliability)) && all(is.na(two$elbow)))
  expect_true(all(is.finite(two$moran_z)))
  e <- .cs_error(get_voronoi_seeds(method = "kmeans", n = two, sample_points = pts))
  expect_identical(
    conditionMessage(e),
    text("c(9, 4.5)",
         paste0(" (a ladder of sizes fits no variogram unless `sac` is given, and `cp` takes ",
                "its noise variance from one)"),
         paste0(" (`moran_z` is finite on it, but it is a diagnostic of the cell means and is ",
                "not read as a count.)")))
  # The selection the message names is a count.
  s <- select_resolution(two, criterion = "size")
  expect_identical(nrow(get_voronoi_seeds(method = "kmeans", n = s, sample_points = pts)), 133L)
  # summary() of that profile says the same of the ladder.
  present <- paste0("summary.resolution_profile(): no count criterion is finite at any level. ",
                    "cp needs a response, reliability a usable variogram; a geometry-only ",
                    "profile carries elbow alone, and elbow is NA when the WSS curve has no elbow.")
  expect_identical(
    conditionMessage(.cs_error(summary(two))),
    paste0(present, " moran_z is finite, and is a diagnostic: name it in `criteria` to see it.",
           " The profile is on a ladder of sizes (obs_per_cell = c(9, 4.5)), which fits no ",
           "variogram unless `sac` is given, so it has no cp and no reliability; ",
           "select_resolution(criterion = \"size\") reads it."))
  # With another criterion finite on it the refusal is the classed one, and
  # it names the size as one more reader.
  bent <- two; bent$elbow <- c(0.5, 0.3)
  e <- .cs_error(get_voronoi_seeds(method = "kmeans", n = bent, sample_points = pts))
  expect_s3_class(e, "spatialkit_no_default_criterion")
  expect_identical(
    conditionMessage(e),
    paste0("get_voronoi_seeds(): `n` is a resolution profile whose `cp`, the criterion a ",
           "profile is read by, is NA at every level, so no cell count is read off it. Finite ",
           "on it: elbow. To use one, pass select_resolution(<profile>, criterion = \"elbow\") ",
           "as `n`; or pass a number. The profile is on a ladder of sizes (obs_per_cell = c(9, ",
           "4.5)): criterion = \"size\" reads it by a size."))
  # On another ladder the classed refusal ends where it does today.
  today <- .cs_plain(); today$cp <- NA_real_; today$elbow[1L] <- 0.5
  e <- .cs_error(get_voronoi_seeds(method = "kmeans", n = today, sample_points = pts))
  expect_s3_class(e, "spatialkit_no_default_criterion")
  expect_match(conditionMessage(e), "as `n`; or pass a number\\.$")
  # A profile of the geometry on the automatic ladder keeps today's words.
  auto <- suppressWarnings(resolution_profile(pts, n_levels = 6L, partition = "kmeans", nstart = 2,
                                              sac = NULL, estimate = NULL, obs_per_cell = NULL,
                                              seed = 123))
  expect_true(all(is.na(auto$elbow)))
  expect_identical(
    conditionMessage(.cs_error(get_voronoi_seeds(method = "kmeans", n = auto, sample_points = pts))),
    paste0("get_voronoi_seeds(): `n` is a resolution profile with no criterion finite at any ",
           "level, so no cell count can be read off it. It is geometry-only, and its WSS curve ",
           "has no elbow: the points have no cluster structure to choose a count. Pass a ",
           "number, or profile with a `response_var`."))
  expect_identical(conditionMessage(.cs_error(summary(auto))), present)
})

test_that(".resolve_cell_count() returns the record with the count", {
  p <- .cs_pooled()
  res <- .resolve_cell_count(p, "n", "get_voronoi_seeds")
  expect_identical(names(res), c("n", "from", "seeds", "criterion", "edge", "at_floor",
                                 "at_ceiling", "structured", "count_rule"))
  expect_identical(res$n, 133)
  expect_identical(res$criterion, "size")
  expect_identical(class(res$count_rule), "spatialkit_count_rule")
  expect_identical(res$count_rule$rule, "size")
  expect_identical(res$count_rule$estimate, "pooled")
  expect_identical(res$count_rule$scored_cells, NA)
  expect_identical(res$count_rule, select_resolution(p, criterion = "size")$count_rule)
  expect_identical(attr(res$seeds, "count_rule")$scored_cells, TRUE)
  # A selection passes through with its own.
  s <- select_resolution(.cs_mean7(), criterion = "loo")
  res <- suppressWarnings(.resolve_cell_count(s, "k", "voronoi_seeds_kmeans"))
  expect_identical(res$count_rule, s$count_rule)
  expect_identical(res$from, paste0("select_resolution(criterion = \"loo\")",
                                    if (is.na(s$edge)) "" else sprintf(", an optimum at %s", s$edge)))
  # A number, and a selection of today, carry none.
  expect_null(.resolve_cell_count(12, "n", "f")$count_rule)
  expect_identical(.resolve_cell_count(12, "n", "f")$n, 12)
  plain <- select_resolution(.cs_plain(), criterion = "cp")
  expect_null(suppressWarnings(.resolve_cell_count(plain, "n", "f"))$count_rule)
  expect_null(suppressWarnings(.resolve_cell_count(.cs_plain(), "n", "f"))$count_rule)
  expect_identical(suppressWarnings(.resolve_cell_count(.cs_plain(), "n", "f"))$criterion, "cp")
})

test_that("a profile of a constant response passed as a count: a size has one, a minimum has none", {
  pts <- .cs_600(); pts$k <- 0.1
  pooled <- suppressWarnings(.cs_profile(pts, response = "k"))
  expect_true(is.na(pooled$loo))
  r <- .cs_catch(get_voronoi_seeds(method = "kmeans", n = pooled, sample_points = pts))
  expect_identical(nrow(r$value), 133L)
  expect_false(any(c("spatialkit_no_structure", "spatialkit_bound_chosen") %in% .cs_classes(r)))
  res <- .resolve_cell_count(pooled, "n", "get_voronoi_seeds")
  expect_identical(res$structured, NA)
  expect_identical(res$count_rule$scored$structured, NA)
  expect_identical(res$count_rule$scored$loo, NA_real_)
  # Plain means: select_resolution()'s own error, passed on as it is.
  plain <- suppressWarnings(.cs_profile(pts, response = "k", estimate = "mean",
                                        obs_per_cell = c(9, 4.5)))
  e <- .cs_error(get_voronoi_seeds(method = "kmeans", n = plain, sample_points = pts))
  expect_identical(conditionMessage(e),
                   paste0("select_resolution(): `loo` is NA at every level because the response ",
                          "is constant: its total sum of squares is rounding error at most, so ",
                          "no cell count fits it better than another."))
  expect_false(inherits(e, "spatialkit_no_default_criterion"))
  e <- .cs_error(voronoi_seeds_kmeans(pts, plain))
  expect_match(conditionMessage(e), "^select_resolution\\(\\): `loo` is NA at every level")
})


# ---- print() ----------------------------------------------------------------------

test_that("a selection by size prints its size, the estimate, the level's numbers and no flat region", {
  p <- .cs_pooled()
  s <- select_resolution(p, criterion = "size")
  out <- .cs_print(s)
  expect_identical(out[1:6], c(
    "Resolution by size: 133 cells",
    "  size        : 4.5 rows per cell asked (obs_per_cell); 600 rows / 133 cells = 4.51",
    "  cells       : density-weighted k-means++ starts (power 1, 8 neighbours)",
    "  estimate    : pooled toward the 8 nearest cells",
    sprintf("  at this size: %.1f%% of the rows are alone in their cell; median own weight %.2f",
            100 * p$share_alone, p$own_weight),
    sprintf("  loo         : %s (one cell for the layer: %s)", format(signif(p$loo, 3)),
            format(signif(attr(p, "estimate")$loo_one_cell, 3)))))
  expect_identical(out[7:8], c(
    "  note        : the count is set by obs_per_cell = 4.5, the caller's choice;",
    "                no criterion chose it."))
  expect_length(out, 8L)
  expect_false(any(grepl("flat region", out, fixed = TRUE)))
  expect_false(any(grepl("the bound is choosing", out, fixed = TRUE)))
  utils::capture.output(shown <- withVisible(print(s)))
  expect_false(shown$visible)
  expect_identical(shown$value, s)
  # Plain means have no own weight.
  m <- .cs_profile(.cs_600(), estimate = "mean", obs_per_cell = 4.5)
  out <- .cs_print(select_resolution(m, criterion = "size"))
  expect_identical(out[4:5], c(
    "  estimate    : plain means",
    sprintf("  at this size: %.1f%% of the rows are alone in their cell", 100 * m$share_alone)))
  # Without an estimate: the size and the note.
  geo <- resolution_profile(.cs_600(), partition = "kmeans", nstart = 2, sac = "auto",
                            estimate = NULL, obs_per_cell = 6, seed = 123)
  expect_identical(.cs_print(select_resolution(geo, criterion = "size")), c(
    "Resolution by size: 100 cells",
    "  size        : 6 rows per cell asked (obs_per_cell); 600 rows / 100 cells = 6.00",
    "  note        : the count is set by obs_per_cell = 6, the caller's choice; no",
    "                criterion chose it."))
})

test_that("under the split the print says which rows the level's numbers are of", {
  half <- .cs_profile(.cs_600(), select_on = "split")
  s <- select_resolution(half, criterion = "size")
  out <- .cs_print(s)
  n_half <- attr(half, "estimate")$n
  # The size line counts the rows of both halves.
  expect_identical(out[2L],
                   "  size        : 4.5 rows per cell asked (obs_per_cell); 600 rows / 133 cells = 4.51")
  # The rows the error read, the cells that hold one of them, and the rows
  # per such cell: 299 rows in 75 of the 133 cells.
  occupied <- unname(attr(half, "estimate")$cells[["133"]])
  expect_lt(occupied, 133L)
  expect_equal(n_half / occupied, half$scored_per_cell)
  line <- paste(sub("^ +", "", sub("^  at this size: ", "", out[5:6])), collapse = " ")
  expect_identical(line,
                   sprintf(paste0("in the selection half (%d rows in %d of the 133 cells, %.2f ",
                                  "per cell): %.1f%% of the rows are alone in their cell; median ",
                                  "own weight %.2f"),
                           n_half, occupied, half$scored_per_cell, 100 * half$share_alone,
                           half$own_weight))
  expect_match(out[5L], "^  at this size: in the selection half \\(")
  expect_match(out[6L], "^                [0-9a-z]")
  # Both numbers of the next line are of that half too.
  expect_identical(out[7L],
                   sprintf("  loo         : %s (one cell for that half: %s)",
                           format(signif(half$loo, 3)),
                           format(signif(attr(half, "estimate")$loo_one_cell, 3))))
})

test_that("the print says where centres hold no scored row, and leaves out what the profile lacks", {
  # One row in eight has a response: 75 rows at 1.2 a cell are 62 centres,
  # fitted to all 600 points, and not every one has a scored row nearest it.
  some <- .cs_600(); some$z[seq_len(600) %% 8 != 0] <- NA
  p <- suppressWarnings(.cs_profile(some, obs_per_cell = 1.2))
  s <- select_resolution(p, criterion = "size")
  expect_identical(s$best, 62L)
  occupied <- s$count_rule$scored$cells
  expect_lt(occupied, 62L)
  out <- .cs_print(s)
  expect_identical(out[2L],
                   "  size        : 1.2 rows per cell asked (obs_per_cell); 75 rows / 62 cells = 1.21")
  line <- paste(sub("^ +", "", sub("^  at this size: ", "", out[5:6])), collapse = " ")
  expect_identical(line,
                   sprintf(paste0("in the %d of the 62 cells that hold a scored row (%.2f per ",
                                  "cell): %.1f%% of the rows are alone in their cell; median own ",
                                  "weight %.2f"),
                           occupied, p$scored_per_cell, 100 * p$share_alone, p$own_weight))
  expect_match(out[7L], "^  loo         : .* \\(one cell for the layer: ")
  # Every centre holds a scored row: nothing is said of cells.
  expect_false(any(grepl("that hold a scored row", .cs_print(
    select_resolution(.cs_pooled(), criterion = "size")), fixed = TRUE)))
  # A constant response has no leave-one-out error at the level and none for
  # one mean: the line is left out, and the share of rows alone is kept.
  pts <- .cs_600(); pts$k <- 0.1
  k <- suppressWarnings(.cs_profile(pts, response = "k"))
  expect_true(is.na(k$loo) && is.na(attr(k, "estimate")$loo_one_cell))
  out <- .cs_print(select_resolution(k, criterion = "size"))
  expect_identical(out[4:5], c(
    "  estimate    : pooled toward the 8 nearest cells",
    sprintf("  at this size: %.1f%% of the rows are alone in their cell", 100 * k$share_alone)))
  expect_false(any(grepl("^  loo  ", out)))
  # Without the share of rows alone the line that carries it is left out.
  q <- .cs_pooled(); q$share_alone <- NA_real_
  out <- .cs_print(select_resolution(q, criterion = "size"))
  expect_false(any(grepl("at this size", out, fixed = TRUE)))
  expect_true(any(grepl("^  loo         : ", out)))
  # With one of the two errors the line gives it, and NA for the other.
  q <- .cs_pooled(); q$loo <- NA_real_
  out <- .cs_print(select_resolution(q, criterion = "size"))
  expect_true(sprintf("  loo         : NA (one cell for the layer: %s)",
                      format(signif(attr(q, "estimate")$loo_one_cell, 3))) %in% out)
  # A size of 1 is one row per cell.  The profile holds it to 599 cells.
  one <- suppressWarnings(.cs_profile(.cs_600(), obs_per_cell = c(4.5, 1)))
  out <- .cs_print(select_resolution(one, criterion = "size", obs_per_cell = 1))
  expect_identical(out[1:2], c(
    "Resolution by size: 599 cells",
    "  size        : 1 row per cell asked (obs_per_cell); 600 rows / 599 cells = 1.00"))
})

test_that("a selection by loo prints the present layout with the estimate's lines added", {
  p <- .cs_mean7()
  s <- select_resolution(p, criterion = "loo")
  out <- .cs_print(s)
  at <- match(s$best, p$levels)
  expect_identical(out[1:6], c(
    sprintf("Resolution by loo: %d cells", s$best),
    sprintf("  flat region : %s (%d of 7 levels)", .band_label(s$flat, p$levels), length(s$flat)),
    "  cells       : density-weighted k-means++ starts (power 1, 8 neighbours)",
    "  estimate    : plain means",
    sprintf("  at this size: %.1f%% of the rows are alone in their cell", 100 * p$share_alone[at]),
    sprintf("  loo         : %s (one cell for the layer: %s)", format(signif(p$loo[at], 3)),
            format(signif(attr(p, "estimate")$loo_one_cell, 3)))))
  # A column subset has no estimate to print.
  out <- .cs_print(select_resolution(p[, c("levels", "loo")], criterion = "loo"))
  expect_false(any(grepl("estimate    :|at this size:|loo         :", out)))
  # A selection by cp on a profile with an estimate prints as it does today,
  # with its note.
  l <- .cs_levels()
  out <- .cs_print(select_resolution(l, criterion = "cp"))
  expect_false(any(grepl("estimate    :|at this size:|loo         :", out)))
  expect_true(any(grepl("note        : cp scores plain means", out, fixed = TRUE)))
})


# ---- summary() and plot() ----------------------------------------------------------

test_that("summary() reads loo where the profile has it, and refuses a size", {
  l <- .cs_levels()
  sm <- summary(l)
  expect_true(all(c("cp", "loo") %in% sm$criterion))
  expect_identical(sm$best[sm$criterion == "loo"], select_resolution(l, criterion = "loo")$best)
  expect_identical(attr(sm, "bands")$loo, select_resolution(l, criterion = "loo")$flat)
  expect_null(attr(sm, "obs_per_cell"))
  # Named alone.
  expect_identical(summary(l, criteria = "loo")$criterion, "loo")
  # A profile without the column has the table it has today.
  expect_false("loo" %in% summary(.cs_plain())$criterion)
  expect_error(summary(.cs_plain(), criteria = "loo"), "has no `loo` column", fixed = TRUE)
  unknown <- function(cr)
    sprintf(paste0("summary.resolution_profile(): unknown criteria: %s. Choose from cp, ",
                   "reliability, elbow, moran_z, loo."), paste(cr, collapse = ", "))
  size_says <- paste0(" A size sets the count and has no flat region to compare; ",
                      "select_resolution(criterion = \"size\") reads it.")
  auto_says <- function(is)
    sprintf(paste0(" \"auto\" names no criterion of the table; on this profile ",
                   "select_resolution() takes it as \"%s\"."), is)
  expect_identical(conditionMessage(.cs_error(summary(l, criteria = "size"))),
                   paste0(unknown("size"), size_says))
  # "auto" is told what it is on the profile in hand: a size for pooled
  # values, the leave-one-out error for plain means, cp without an estimate.
  expect_identical(conditionMessage(.cs_error(summary(l, criteria = "auto"))),
                   paste0(unknown("auto"), auto_says("size"), size_says))
  expect_identical(conditionMessage(.cs_error(summary(.cs_mean7(), criteria = "auto"))),
                   paste0(unknown("auto"), auto_says("loo")))
  expect_identical(conditionMessage(.cs_error(summary(.cs_plain(), criteria = "auto"))),
                   paste0(unknown("auto"), auto_says("cp")))
  expect_identical(conditionMessage(.cs_error(summary(.cs_mean7(), criteria = c("size", "auto")))),
                   paste0(unknown(c("size", "auto")), auto_says("loo"), size_says))
  expect_identical(conditionMessage(.cs_error(summary(l, criteria = "nope"))),
                   paste0("summary.resolution_profile(): unknown criteria: nope. Choose from cp, ",
                          "reliability, elbow, moran_z, loo."))
  # A ladder of several sizes names the end an optimum sits on, as any ladder.
  m <- .cs_mean7(); m$loo <- 1:7
  out <- .cs_print(summary(m))
  expect_true(any(grepl("loo: the optimum is the coarsest size offered (obs_per_cell).", out,
                        fixed = TRUE)))
  expect_true(any(grepl("There the bound is choosing, not the criterion.", out, fixed = TRUE)))
})

test_that("for a ladder of one size summary() says the size set the count", {
  pts <- .cs_600()
  line <- c("  the ladder is one size (obs_per_cell = 4.5); the count is set by the size,",
            "    and no criterion chooses among levels.")
  # No variogram: the leave-one-out error is the one criterion finite there.
  p <- .cs_pooled()
  sm <- summary(p)
  expect_identical(sm$criterion, "loo")
  expect_identical(attr(sm, "obs_per_cell"), 4.5)
  expect_identical(sm$edge, "the only level scored")
  out <- .cs_print(sm)
  expect_identical(utils::tail(out, 2L), line)
  expect_false(any(grepl("bound is choosing|only level scored", out)))
  # With a variogram: one row for each criterion finite at the level.
  v <- suppressWarnings(.cs_profile(pts, sac = 3e4))
  finite <- c("cp", "reliability", "elbow", "loo")[
    vapply(c("cp", "reliability", "elbow", "loo"), function(cn) is.finite(v[[cn]]), NA)]
  expect_true(all(c("cp", "loo") %in% finite))
  sm <- summary(v)
  expect_identical(sm$criterion, finite)
  out <- .cs_print(sm)
  expect_identical(out[match(line[1L], out) + 0:1], line)
  expect_false(any(grepl("bound is choosing|only level scored", out)))
  # A count one cell matches is still said: white noise on 300 rows (seed
  # 4001).
  noise <- .cs_profile(.cs_layer(300, 4001, field = FALSE))
  sm <- summary(noise)
  expect_identical(sm$structured, FALSE)
  out <- .cs_print(sm)
  expect_identical(utils::tail(out, 4L), c(
    line,
    "  loo: one cell for the whole layer does as well, so the response shows no",
    "    spatial structure the cells can use and this count is not identified."))
  # A constant response: no criterion is finite, and the error ends with the
  # line and names the reader of the profile.
  pts$k <- 0.1
  k <- suppressWarnings(.cs_profile(pts, response = "k"))
  e <- .cs_error(summary(k))
  expect_match(conditionMessage(e),
               "^summary.resolution_profile\\(\\): no count criterion is finite at any level\\. ")
  ending <- paste0(" The ladder is one size (obs_per_cell = 4.5); the count is set by the size, ",
                   "and select_resolution(criterion = \"size\") reads it.")
  expect_identical(substring(conditionMessage(e), nchar(conditionMessage(e)) - nchar(ending) + 1L),
                   ending)
  # Several sizes, and every other ladder: no such line, and the present error.
  expect_null(attr(summary(.cs_mean7()), "obs_per_cell"))
  expect_false(any(grepl("the ladder is one size", .cs_print(summary(.cs_mean7())), fixed = TRUE)))
  k7 <- suppressWarnings(.cs_profile(pts, response = "k", obs_per_cell = c(9, 4.5)))
  e <- .cs_error(summary(k7))
  expect_false(grepl("The ladder is one size", conditionMessage(e), fixed = TRUE))
  expect_match(conditionMessage(e), "no count criterion is finite at any level", fixed = TRUE)
  # Several sizes are named as a ladder of sizes, with the missing variogram
  # where the profile fitted none and without it where the caller gave one.
  present <- paste0("summary.resolution_profile(): no count criterion is finite at any level. ",
                    "cp needs a response, reliability a usable variogram; a geometry-only ",
                    "profile carries elbow alone, and elbow is NA when the WSS curve has no elbow.")
  expect_identical(
    conditionMessage(e),
    paste0(present, " The profile is on a ladder of sizes (obs_per_cell = c(9, 4.5)), which ",
           "fits no variogram unless `sac` is given, so it has no cp and no reliability; ",
           "select_resolution(criterion = \"size\") reads it."))
  kv <- suppressWarnings(.cs_profile(pts, response = "k", obs_per_cell = c(9, 4.5), sac = 3e4))
  expect_false(identical(attr(kv, "cp_noise")$source, "none"))
  expect_identical(
    conditionMessage(.cs_error(summary(kv))),
    paste0(present, " The profile is on a ladder of sizes (obs_per_cell = c(9, 4.5)); ",
           "select_resolution(criterion = \"size\") reads it."))
})

test_that("plot() draws loo like a count criterion, and a ladder of one size as points", {
  skip_if_not_installed("ggplot2")
  geoms <- function(g) vapply(g$layers, function(l) class(l$geom)[1L], character(1))
  points_of <- function(g) g$layers[[which(geoms(g) == "GeomPoint")[1L]]]$data
  label <- "Leave-one-out error (lower is better)"
  # Several levels: the selected level and the band, as for cp.
  l <- .cs_levels()
  g <- plot(l, criteria = "loo")
  expect_s3_class(g, "ggplot")
  expect_true(all(c("GeomRect", "GeomVline") %in% geoms(g)))
  d <- points_of(g)
  expect_identical(as.character(unique(d$panel)), label)
  expect_identical(d$x, as.numeric(l$levels))
  expect_identical(d$y, l$loo)
  marked <- g$layers[[which(geoms(g) == "GeomVline")]]$data
  expect_identical(marked$x, select_resolution(l, criterion = "loo")$best)
  expect_match(g$labels$caption, "^Red: the level each criterion selects")
  # Drawn by default where the column is finite, after the other count criteria.
  panels <- as.character(unique(points_of(plot(l))$panel))
  expect_true(label %in% panels)
  expect_lt(match("Mallows' Cp (lower is better)", panels), match(label, panels))
  expect_false(label %in% as.character(unique(points_of(plot(.cs_plain()))$panel)))
  expect_error(plot(l, criteria = "size"),
               "unknown criteria: size. Choose from cp, reliability, elbow, loo, moran_z, wss.",
               fixed = TRUE)
  # One size: one point for each finite column, no band, no selected level.
  caption <- paste0("the ladder is one size (obs_per_cell = 4.5); the count is set by the size,",
                    "\nand no criterion chooses among levels.")
  v <- suppressWarnings(.cs_profile(.cs_600(), sac = 3e4))
  g <- plot(v)
  expect_false(any(c("GeomRect", "GeomVline") %in% geoms(g)))
  expect_identical(sum(geoms(g) == "GeomPoint"), 1L)
  d <- points_of(g)
  shown <- c("cp", "reliability", "elbow", "loo", "moran_z")
  shown <- shown[vapply(shown, function(cn) is.finite(v[[cn]]), NA)]
  expect_true(all(c("cp", "loo", "moran_z") %in% shown))
  expect_identical(nrow(d), length(shown))
  expect_identical(d$x, rep(133, length(shown)))
  expect_identical(d$y, vapply(shown, function(cn) abs(v[[cn]]), numeric(1), USE.NAMES = FALSE))
  expect_identical(g$labels$caption, caption)
  p <- .cs_pooled()
  g <- plot(p)
  expect_false(any(c("GeomRect", "GeomVline") %in% geoms(g)))
  expect_identical(g$labels$caption, caption)
  expect_identical(nrow(points_of(g)), sum(is.finite(c(p$loo, p$moran_z))))
  g <- plot(p, criteria = "loo")
  expect_identical(points_of(g)$y, p$loo)
  expect_false("GeomVline" %in% geoms(g))
  expect_identical(g$labels$subtitle, "1 levels from 133 to 133 cells (floor 2, ceiling 133)")
  # A size too coarse for the layer: the subtitle says why the profile is not
  # supported.  Six rows at 4.5 a cell.
  few <- suppressWarnings(.cs_profile(.cs_600()[1:6, ]))
  expect_false(attr(few, "bounds")$supported)
  expect_identical(plot(few, criteria = "wss")$labels$subtitle,
                   paste0("1 levels from 2 to 2 cells (floor 2, ceiling 2); fewer than two ",
                          "cells of the size asked for"))
})

test_that("a profile of one level that is not a ladder of one size is summarised and drawn as today", {
  skip_if_not_installed("ggplot2")
  for (estimate in list(NULL, "pooled")) {
    what <- if (is.null(estimate)) "no estimate" else estimate
    o <- suppressWarnings(.cs_profile(.cs_600(), levels = 40, partition = "kmeans", nstart = 5,
                                      sac = 3e4, estimate = estimate, obs_per_cell = NULL))
    expect_identical(o$levels, 40L, info = what)
    expect_null(attr(o, "bounds")$obs_per_cell, info = what)
    sm <- summary(o)
    expect_null(attr(sm, "obs_per_cell"), info = what)
    expect_true("cp" %in% sm$criterion, info = what)
    expect_identical("loo" %in% sm$criterion, !is.null(estimate), info = what)
    expect_true(all(sm$edge == "the only level scored"), info = what)
    out <- .cs_print(sm)
    expect_true(sprintf("  %s: the optimum is the only level scored.",
                        paste(sm$criterion, collapse = ", ")) %in% out, info = what)
    expect_true("  There the bound is choosing, not the criterion." %in% out, info = what)
    expect_false(any(grepl("the ladder is one size", out, fixed = TRUE)), info = what)
    g <- plot(o)
    expect_true(all(c("GeomRect", "GeomLine", "GeomVline") %in% .cs_geoms(g)), info = what)
    expect_match(g$labels$caption, "^Red: the level each criterion selects", info = what)
    expect_match(g$labels$caption, "cp: optimum at the only level scored (the bound is choosing)",
                 fixed = TRUE, info = what)
    expect_false(grepl("the ladder is one size", g$labels$caption, fixed = TRUE), info = what)
  }
})

test_that("plot() marks a level on a ladder of several sizes, and draws no line through one point", {
  skip_if_not_installed("ggplot2")
  # Seven sizes are a ladder to choose on, drawn as any other.
  m <- .cs_mean7()
  s <- select_resolution(m, criterion = "loo")
  g <- plot(m, criteria = "loo")
  expect_true(all(c("GeomRect", "GeomLine", "GeomVline") %in% .cs_geoms(g)))
  expect_identical(g$layers[[which(.cs_geoms(g) == "GeomVline")]]$data$x, s$best)
  expect_match(g$labels$caption, "^Red: the level each criterion selects")
  expect_false(grepl("the ladder is one size", g$labels$caption, fixed = TRUE))
  # One size: points alone.  A line through one observation draws nothing,
  # and ggplot2 says so for each panel when the plot is printed.
  one <- "Each group consists of only one observation"
  for (p in list(.cs_pooled(), suppressWarnings(.cs_profile(.cs_600(), sac = 3e4)))) {
    g <- plot(p)
    expect_identical(.cs_geoms(g), "GeomPoint")
    expect_false(any(grepl(one, .cs_drawn(g), fixed = TRUE)))
  }
})
