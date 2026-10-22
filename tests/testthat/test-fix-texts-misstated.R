# BUG-07: messages and prints that misstated what happened in particular
# cases.  One block per item of the backlog's list; each test fails on the
# code before the fix.  Item 3 (the pointer of determine_optimal_levels()'s
# ranking warning) was done with BUG-05 and is tested in
# test-fix-levels-combined.R.

# The layers of the frozen fixture: 400 points with a structured response
# `z`, a white-noise column `noise`, a logical `flag` and an empty `none`.
.tm_store <- new.env(parent = emptyenv())
.tm_layers <- function() {
  if (is.null(.tm_store$layers)) {
    fx <- readRDS(testthat::test_path("fixtures", "core-frozen.rds"))
    .tm_store$layers <- core_frozen_layers(fx$data)
  }
  .tm_store$layers
}
.tm_pts <- function() .tm_layers()$pts

# A call with its warnings and messages kept and muffled.
.tm_catch <- function(expr) {
  ws <- list(); ms <- character(0)
  value <- withCallingHandlers(expr,
    warning = function(w) { ws[[length(ws) + 1L]] <<- w; invokeRestart("muffleWarning") },
    message = function(m) { ms <<- c(ms, conditionMessage(m)); invokeRestart("muffleMessage") })
  list(value = value, warnings = ws, messages = ms,
       classes = vapply(ws, function(w) class(w)[1L], ""),
       text = vapply(ws, conditionMessage, ""))
}
.tm_of <- function(r, cls) Filter(function(w) inherits(w, cls), r$warnings)
.tm_error <- function(expr)
  tryCatch({ suppressWarnings(suppressMessages(expr)); NULL }, error = function(e) e)
.tm_quiet <- function(expr) suppressWarnings(suppressMessages(expr))
.tm_print <- function(x) utils::capture.output(print(x))
.tm_without_fnn <- function(env = parent.frame())
  testthat::local_mocked_bindings(
    .has_backend = function(p) if (identical(p, "FNN")) FALSE else requireNamespace(p, quietly = TRUE),
    .env = env)
.tm_fnn_remedy <- "FNN, which installs with spatialkit, could not be loaded: reinstall it; the result is the same."


# ---- 1. FNN installs with spatialkit ------------------------------------------------

test_that("no message asks the user to install FNN, which installs with spatialkit", {
  ns <- asNamespace("spatialkit")
  fns <- Filter(function(nm) is.function(get(nm, envir = ns)), ls(ns, all.names = TRUE))
  body_of <- function(nm) paste(deparse(get(nm, envir = ns)), collapse = "\n")
  # Nor says that it is not installed: it is, and could not be loaded.
  asks <- c("Install FNN", "install FNN", "install.packages(\\\"FNN\\\"",
            "install.packages(c(\\\"FNN\\\"", "the FNN package installed",
            "FNN is not installed", "FNN and Matrix are not installed")
  hits <- Filter(function(nm) any(vapply(asks, grepl, NA, x = body_of(nm), fixed = TRUE)), fns)
  expect_identical(hits, character(0))
  # The constant the work limit's errors end with.
  expect_identical(spatialkit:::.FNN_REMEDY, .tm_fnn_remedy)
})

test_that("the weights of Moran's z name what is missing, and ask to install Matrix alone", {
  xy <- local({ set.seed(71); matrix(stats::runif(2 * 5001), ncol = 2) })
  e <- .tm_error(spatialkit:::.build_knn_weights(xy, k = 8, use_fnn = TRUE, use_matrix = FALSE))
  msg <- conditionMessage(e)
  expect_match(msg, "^n = 5001 requires FNN for k-NN weights, and Matrix to hold them sparsely")
  expect_match(msg, "Matrix could not be loaded: install it with install.packages(\"Matrix\").",
               fixed = TRUE)
  expect_false(grepl("FNN\"", msg, fixed = TRUE))
  expect_false(grepl("reinstall", msg, fixed = TRUE))
  e <- .tm_error(spatialkit:::.build_knn_weights(xy, k = 8, use_fnn = FALSE, use_matrix = TRUE))
  msg <- conditionMessage(e)
  expect_match(msg, "FNN, which installs with spatialkit, could not be loaded: reinstall it.",
               fixed = TRUE)
  expect_false(grepl("install.packages", msg, fixed = TRUE))
  e <- .tm_error(spatialkit:::.build_knn_weights(xy, k = 8, use_fnn = FALSE, use_matrix = FALSE))
  msg <- conditionMessage(e)
  expect_match(msg, "could not be loaded: reinstall it. Matrix could not be loaded: install it",
               fixed = TRUE)
})


# ---- 2. The work limits print the work and the limit in full ------------------------

test_that("a work 0.04 percent over the limit prints two different numbers", {
  e <- .tm_error(spatialkit:::.work_guard(1.0004e10, 1e10, "cell neighbours", "f",
                                          search = "a search"))
  expect_s3_class(e, "spatialkit_work_limit")
  expect_identical(conditionMessage(e), paste0(
    "f(): a search would compare 10,004,000,000 pairs without the FNN package, above the ",
    "limit of 10,000,000,000 (options(spatialkit.exact_pairs_max = )). ", .tm_fnn_remedy))
  expect_identical(spatialkit:::.count_text(c(79800, 79799, 665733.7)),
                   c("79,800", "79,799", "665,733.7"))
})

test_that("a limit less than one below the work is written with the digits that tell them apart", {
  # Seven significant digits wrote both as "10,000,000,000".
  expect_identical(spatialkit:::.count_text(c(1e10, 1e10 - 0.4)),
                   c("10,000,000,000", "9,999,999,999.6"))
  e <- .tm_error(spatialkit:::.work_guard(1e10, 1e10 - 0.4, "cell neighbours", "f",
                                          search = "a search"))
  expect_identical(conditionMessage(e), paste0(
    "f(): a search would compare 10,000,000,000 pairs without the FNN package, above the ",
    "limit of 9,999,999,999.6 (options(spatialkit.exact_pairs_max = )). ", .tm_fnn_remedy))
  # Numbers that seven digits tell apart, or that are equal, are written as before.
  expect_identical(spatialkit:::.count_text(c(1.0004e10, 1e10)),
                   c("10,004,000,000", "10,000,000,000"))
  expect_identical(spatialkit:::.count_text(c(665733.7, 665733.7)), c("665,733.7", "665,733.7"))
})

test_that("the profile's fit guard prints the work and the limit in full, and names no installed FNN", {
  p1000 <- withr::with_seed(7, sf::st_as_sf(
    data.frame(x = stats::runif(1000, 0, 1e4), y = stats::runif(1000, 0, 1e4),
               z = stats::rnorm(1000)), coords = c("x", "y"), crs = 32632))
  # 222 cells of the default size, 3 density-weighted starts: 666,000.
  work <- 1000 * 222 * 3
  e <- withr::with_options(list(spatialkit.fit_work_max = work / 1.0004),
                           .tm_error(resolution_profile(p1000, "z")))
  expect_s3_class(e, "spatialkit_work_limit")
  expect_identical(e$work, work)
  msg <- conditionMessage(e)
  expect_match(msg, "is 666,000 point-to-centre comparisons, above the limit of 665,733.7 ",
               fixed = TRUE)
  expect_false(grepl("6.7e+05", msg, fixed = TRUE))
  expect_false(grepl("FNN package installed", msg, fixed = TRUE))
  expect_match(msg, paste0("partition = \"density\" with init = \"random\", which has no step ",
                           "that compares every point with every centre and leaves a cell ",
                           "without a point more often; or a higher limit."), fixed = TRUE)
  # Where FNN could not be loaded, the way on says that it needs FNN.
  local({
    .tm_without_fnn()
    work <- 1000 * 222 * 3 * 100
    e <- withr::with_options(list(spatialkit.fit_work_max = work - 1),
                             .tm_error(resolution_profile(p1000, "z")))
    expect_false(e$fnn)
    expect_match(conditionMessage(e),
                 paste0("is 66,600,000 point-to-centre comparisons, above the limit of ",
                        "66,599,999 "), fixed = TRUE)
    expect_match(conditionMessage(e),
                 paste0("more often (it needs FNN, which installs with spatialkit and could not ",
                        "be loaded: reinstall it); or a higher limit."), fixed = TRUE)
  })
})


# ---- 4. A constant response -----------------------------------------------------------

test_that("on a constant response the warning names loo, and summary() says the response is constant", {
  const <- .tm_pts(); const$z <- 5
  r <- .tm_catch(resolution_profile(const, "z"))
  w <- .tm_of(r, "spatialkit_constant_response")
  expect_length(w, 1L)
  expect_match(conditionMessage(w[[1L]]),
               "`cp`, `moran_i`, `moran_z`, `loo`, `loo_se` and `own_weight` are NA at every level.",
               fixed = TRUE)
  expect_true(all(is.na(r$value$loo)) && all(is.na(r$value$own_weight)))
  e <- .tm_error(summary(r$value))
  expect_identical(conditionMessage(e), paste0(
    "summary.resolution_profile(): no count criterion is finite at any level. The response is ",
    "constant over the 400 scored rows, so no criterion that reads it is finite: cp, loo and ",
    "moran_z are NA at every level, and reliability needs a usable variogram of it; elbow, read ",
    "off the points alone, is NA when the WSS curve has no elbow. The ladder is one size ",
    "(obs_per_cell = 4.5); the count is set by the size, and select_resolution(criterion = ",
    "\"size\") reads it."))
  # Plain means: `loo` and `loo_se`, and no `own_weight`, which plain means never have.
  r <- .tm_catch(resolution_profile(const, "z", estimate = "mean"))
  expect_match(conditionMessage(.tm_of(r, "spatialkit_constant_response")[[1L]]),
               "`cp`, `moran_i`, `moran_z`, `loo` and `loo_se` are NA at every level.",
               fixed = TRUE)
  # A profile with no response keeps the general sentence.
  none <- .tm_quiet(resolution_profile(.tm_pts(), "none", estimate = NULL, obs_per_cell = NULL,
                                       sac = NULL, partition = "kmeans"))
  e <- .tm_error(summary(none))
  expect_match(conditionMessage(e),
               "cp needs a response, reliability a usable variogram; a geometry-only profile",
               fixed = TRUE)
  expect_false(grepl("constant", conditionMessage(e), fixed = TRUE))
})


# ---- 5. The no-structure warnings ---------------------------------------------------

test_that("the no-structure warning calls a predictor a column, not the response", {
  pts <- .tm_pts()
  asg <- .tm_quiet(assign_features_to_polygons(pts, .tm_layers()$units, polygon_id_col = "GEOID"))
  # The response alone: "the response".
  r <- .tm_catch(summarize_by_cell(asg, "noise", id_col = "GEOID", estimate = "pooled"))
  w <- .tm_of(r, "spatialkit_no_structure")
  expect_length(w, 1L)
  expect_identical(w[[1L]]$columns, "noise")
  expect_match(conditionMessage(w[[1L]]), "not because the response varies smoothly.",
               fixed = TRUE)
  # White noise as a predictor beside a structured response: the predictor
  # alone is named, and the map is said not to show a smooth column.
  r <- .tm_catch(summarize_by_cell(asg, "z", "noise", id_col = "GEOID", estimate = "pooled"))
  w <- .tm_of(r, "spatialkit_no_structure")
  expect_length(w, 1L)
  expect_identical(w[[1L]]$columns, "noise")
  expect_identical(unname(w[[1L]]$beats), FALSE)
  expect_match(conditionMessage(w[[1L]]), "not because the column varies smoothly.", fixed = TRUE)
  expect_false(grepl("the response", conditionMessage(w[[1L]]), fixed = TRUE))
  # The response and a predictor of white noise, neither of which beats one
  # mean: "the columns".
  pts$p1 <- withr::with_seed(3, stats::rnorm(nrow(pts)))
  asg <- .tm_quiet(assign_features_to_polygons(pts, .tm_layers()$units, polygon_id_col = "GEOID"))
  r <- .tm_catch(summarize_by_cell(asg, "noise", "p1", id_col = "GEOID", estimate = "pooled"))
  w <- .tm_of(r, "spatialkit_no_structure")
  expect_length(w, 1L)
  expect_identical(w[[1L]]$beats, c(noise = FALSE, p1 = FALSE))
  expect_match(conditionMessage(w[[1L]]), "not because the columns vary smoothly.", fixed = TRUE)
  expect_false(grepl("the response varies", conditionMessage(w[[1L]]), fixed = TRUE))
})

test_that("the variogram's note after Cp's is a sentence of its own that does not repeat its close", {
  pts <- .tm_pts()
  pn <- .tm_quiet(resolution_profile(pts, "noise", estimate = NULL, obs_per_cell = NULL,
                                     sac = NULL, partition = "kmeans"))
  s <- select_resolution(pn, "cp")
  expect_false(s$structured)
  expect_length(s$note, 2L)
  expect_match(s$note[[2L]],
               "^The profile's variogram found no spatial structure in the response either \\(")
  r <- .tm_catch(get_voronoi_seeds(method = "kmeans", n = pn, sample_points = pts,
                                   partition = "kmeans"))
  msg <- conditionMessage(.tm_of(r, "spatialkit_no_structure")[[1L]])
  expect_false(grepl(". the ", msg, fixed = TRUE))
  expect_identical(lengths(regmatches(msg, gregexpr("not identified", msg, fixed = TRUE))), 1L)
  # Alone, the note keeps its words.
  q <- pn
  attr(q, "one_cell") <- NULL
  s <- select_resolution(q, "cp")
  expect_identical(length(s$note), 1L)
  expect_match(s$note, paste0("^the profile's variogram found no spatial structure in the ",
                              "response \\(.*\\), so the count is not identified\\.$"))
})

test_that("where Cp beats one cell and the profile's variogram found no structure, the notes and summary() say so once", {
  pts <- .tm_pts()
  # A variogram refused for showing no spatial structure, given to a profile
  # of the structured `z`: Cp beats one cell by far more than `tol`, and the
  # F test tells the cells from one mean.
  vm <- data.frame(model = c("Nug", "Exp"), psill = c(0.05, 0.95), range = c(0, 10),
                   stringsAsFactors = FALSE)
  why <- "no spatial structure detected at the nearest-neighbour scale"
  nosac <- structure(NA_real_, class = c("sac_range", "numeric"), variogram_model = vm,
                     crs = sf::st_crs(pts), rejected_reason = why)
  p <- .tm_quiet(resolution_profile(pts, "z", sac = nosac, estimate = NULL, obs_per_cell = NULL,
                                    partition = "kmeans", nstart = 3L))
  s <- select_resolution(p, "cp")
  expect_gt(s$one_cell, s$value * 1.02)
  expect_lt(s$structure_p, 0.05)
  expect_false(s$structured)
  # Alone, the variogram's note keeps its words.
  expect_identical(s$note, sprintf(paste0("the profile's variogram found no spatial structure in ",
                                          "the response (%s), so the count is not identified."),
                                   why))
  # summary() says why the count is not identified, in place of "one cell
  # does as well", which is untrue here.
  sm <- summary(p, criteria = "cp")
  expect_identical(sm$structured, FALSE)
  expect_identical(attr(sm, "not_identified"),
                   c(cp = "the profile's variogram found no spatial structure in the response"))
  out <- .tm_print(sm)
  said <- strwrap(paste0("cp: the cells beat one cell for the whole layer, but the profile's ",
                         "variogram found no spatial structure in the response, so no spatial ",
                         "structure the cells can use is shown and this count is not identified."),
                  width = 78, prefix = "  ", exdent = 2)
  expect_identical(utils::tail(out, length(said)), said)
  expect_false(any(grepl("does as well", out, fixed = TRUE)))

  # After a note that does not say "not identified" (no F test at the
  # optimum), the variogram's note is a sentence of its own with its close.
  q <- p
  one <- attr(q, "one_cell"); one$f_p <- NULL; attr(q, "one_cell") <- one
  s <- select_resolution(q, "cp")
  expect_length(s$note, 2L)
  expect_match(s$note[[1L]], paste0("^Cp at one cell is above the optimum by more than `tol`, but ",
                                    "the profile carries no F test"))
  expect_false(grepl("not identified", s$note[[1L]], fixed = TRUE))
  expect_identical(s$note[[2L]], sprintf(paste0("The profile's variogram found no spatial ",
                                                "structure in the response (%s), so the count ",
                                                "is not identified."), why))

  # After the F test's note, which says "not identified" already, the
  # variogram's note says "either" and does not say it again.  The F test's
  # p-value is set to 0.5 at every level: the gap is then one chance gives.
  q <- p
  one <- attr(q, "one_cell"); one$f_p[] <- 0.5; attr(q, "one_cell") <- one
  s <- select_resolution(q, "cp")
  expect_length(s$note, 2L)
  expect_match(s$note[[1L]], paste0("^the cell means at [0-9]+ cells differ from one mean for the ",
                                    "whole layer by no more than chance allows"))
  expect_identical(s$note[[2L]], sprintf(paste0("The profile's variogram found no spatial ",
                                                "structure in the response either (%s)."), why))
  all_notes <- paste(s$note, collapse = " ")
  expect_identical(lengths(regmatches(all_notes, gregexpr("not identified", all_notes,
                                                          fixed = TRUE))), 1L)
  # summary() names both reasons.
  sm <- summary(q, criteria = "cp")
  expect_identical(attr(sm, "not_identified"),
                   c(cp = paste0("the F test does not tell the cell means from one mean ",
                                 "(p = 1), and the profile's variogram found no spatial ",
                                 "structure in the response")))
})


# ---- 6. cell_values_at() on a logical column; empty geometry ---------------------------

test_that("cell_values_at() on a logical column says to convert it, not to summarise it", {
  lg <- .tm_pts(); lg$flag <- lg$z > 0
  sel  <- .tm_quiet(select_resolution(resolution_profile(lg, "z")))
  seeds <- .tm_quiet(get_voronoi_seeds(method = "kmeans", n = sel, sample_points = lg))
  tess <- .tm_quiet(build_tessellation(seeds, method = "voronoi"))
  asg  <- .tm_quiet(assign_features_to_polygons(lg, tess$cells, polygon_id_col = "cell_id"))
  vals <- .tm_quiet(summarize_by_cell(asg, "flag", id_col = "cell_id", cells_sf = tess$cells))
  expect_false(any(grepl("flag", names(vals), fixed = TRUE)))
  e <- .tm_error(cell_values_at(lg[1:3, ], vals, asg, "flag", id_col = "cell_id"))
  expect_identical(conditionMessage(e), paste0(
    "cell_values_at(): 'flag' is a logical column of `points`, and summarize_by_cell() ",
    "summarises numeric columns only, so `values` has no cell value of it. Convert it to ",
    "numbers first (as.numeric()) and summarise it again."))
})

test_that("cell_values_at() on a factor, character or date column says how to make numbers of it", {
  units <- .tm_layers()$units
  pts <- .tm_pts()
  labels <- c("10", "20", "30")[(seq_len(nrow(pts)) %% 3L) + 1L]
  pts$fac <- factor(labels)
  pts$ord <- factor(labels, ordered = TRUE)
  pts$lab <- labels
  pts$day <- as.Date("2026-01-01") + seq_len(nrow(pts))
  asg <- .tm_quiet(assign_features_to_polygons(pts, units, polygon_id_col = "GEOID"))
  vals <- .tm_quiet(summarize_by_cell(asg, "z", id_col = "GEOID", cells_sf = units))
  said <- function(col) conditionMessage(.tm_error(cell_values_at(pts[1:3, ], vals, asg, col,
                                                                  id_col = "GEOID")))
  head_of <- function(col, cls) sprintf(paste0(
    "cell_values_at(): '%s' is a %s column of `points`, and summarize_by_cell() summarises ",
    "numeric columns only, so `values` has no cell value of it. "), col, cls)
  # A factor: as.numeric() gives the level codes, which the plain advice
  # would have had summarised.
  expect_identical(as.numeric(pts$fac[1:3]), c(2, 3, 1))
  expect_identical(as.numeric(as.character(pts$fac[1:3])), c(20, 30, 10))
  factor_way <- paste0("On a factor as.numeric() gives the level codes, not the labels: where ",
                       "the labels are numbers, convert it with as.numeric(as.character()) and ",
                       "summarise it again; a factor of categories has no cell mean.")
  expect_identical(said("fac"), paste0(head_of("fac", "factor"), factor_way))
  expect_identical(said("ord"), paste0(head_of("ord", "factor"), factor_way))
  # Text: numbers written as text convert; labels have no mean.
  expect_identical(said("lab"), paste0(head_of("lab", "character"),
                                       "Where it holds numbers written as text, convert it with ",
                                       "as.numeric() and summarise it again; a column of labels ",
                                       "has no cell mean."))
  # Any other class: the general advice, with no function named.
  expect_identical(said("day"), paste0(head_of("day", "Date"),
                                       "Convert it to numbers first and summarise it again."))
})

test_that("features with empty geometry are said to have no location, not to fall outside the polygons", {
  units <- .tm_layers()$units
  e <- .tm_pts()[1:40, ]
  sf::st_geometry(e)[1:2] <- sf::st_point()
  r <- .tm_catch(assign_features_to_polygons(e, units, polygon_id_col = "GEOID"))
  w <- .tm_of(r, "spatialkit_unassigned")
  expect_length(w, 1L)
  expect_identical(conditionMessage(w[[1L]]), paste0(
    "assign_features_to_polygons(): 2 of 40 feature(s) (5.0%) fall inside no polygon and are ",
    "left out of the result; attr(, \"unassigned\")$which gives their rows in `features_sf`. ",
    "They have empty geometry, so there is no location to assign; pass keep_unassigned = TRUE ",
    "to keep them with an NA ID."))
  # One empty point beside one outside the polygons: both counted, the advice for the rest.
  f <- .tm_pts()[1:40, ]
  sf::st_geometry(f)[1] <- sf::st_point()
  sf::st_geometry(f)[2] <- sf::st_point(c(0, 0))
  r <- .tm_catch(assign_features_to_polygons(f, units, polygon_id_col = "GEOID"))
  expect_match(conditionMessage(.tm_of(r, "spatialkit_unassigned")[[1L]]),
               "rows in `features_sf`. 1 of them has empty geometry; for the rest, check that the polygons cover the features",
               fixed = TRUE)
  # Kept with an NA ID, the log line counts the empty ones too.
  lines <- capture_spatialkit_log(.tm_quiet(
    assign_features_to_polygons(f, units, polygon_id_col = "GEOID", keep_unassigned = TRUE)))
  expect_true(log_has(lines, paste0("2 of 40 feature\\(s\\) fall inside no polygon and are kept ",
                                    "with an NA ID \\(1 of them has empty geometry\\)\\.")))
  # Every feature empty: no word on the CRS.
  g <- .tm_pts()[1:3, ]
  sf::st_geometry(g)[1:3] <- sf::st_point()
  r <- .tm_catch(assign_features_to_polygons(g, units, polygon_id_col = "GEOID"))
  expect_identical(conditionMessage(.tm_of(r, "spatialkit_unassigned")[[1L]]), paste0(
    "assign_features_to_polygons(): none of the 3 feature(s) fall inside any polygon. They have ",
    "empty geometry, so there is no location to assign."))
})


# ---- 7. The note on one mean; one warning for the sizes held --------------------------

test_that("where one mean is below the optimum the notes do not say \"within 2%\"", {
  pts <- .tm_pts()
  s <- select_resolution(.tm_quiet(resolution_profile(pts, "noise", estimate = "mean")), "loo")
  expect_lt(s$one_cell, s$value)
  expect_false(any(grepl("within", s$note, fixed = TRUE)))
  expect_match(s$note[["structure"]], sprintf(
    "^the leave-one-out error at one cell \\(%s\\) is not above the optimum \\(%s\\) by more than 2%%: ",
    format(signif(s$one_cell, 3)), format(signif(s$value, 3))))
  # Cp's note, which the loo note copies.
  s <- select_resolution(.tm_quiet(resolution_profile(pts[1:5, ], "z", estimate = NULL,
                                                      obs_per_cell = NULL, sac = NULL,
                                                      partition = "kmeans")), "cp")
  expect_lt(s$one_cell, s$value)
  expect_false(any(grepl("within", s$note, fixed = TRUE)))
  expect_true(any(startsWith(s$note, sprintf(
    "Cp at one cell (%s) is not above the optimum (%s) by more than 2%%: ",
    format(signif(s$one_cell, 3)), format(signif(s$value, 3))))))
})

test_that("the sizes held to the same count raise one spatialkit_size_not_reached that names them", {
  # 50 stations visited six times.
  st <- withr::with_seed(11, {
    s <- data.frame(x = stats::runif(50, 0, 1e4), y = stats::runif(50, 0, 1e4))
    d <- s[rep(1:50, each = 6), ]
    d$z <- stats::rnorm(300)
    sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  })
  r <- .tm_catch(resolution_profile(st, "z", estimate = "mean", obs_per_cell = "auto"))
  w <- .tm_of(r, "spatialkit_size_not_reached")
  expect_length(w, 1L)
  expect_identical(conditionMessage(w[[1L]]), paste0(
    "resolution_profile(): obs_per_cell = 4.5, 3 and 2 ask for 67 to 150 cells, more than ",
    "k-means can place on the layer's 50 distinct locations; 50 cells are used for each ",
    "(6 rows per cell)."))
  expect_identical(w[[1L]]$obs_per_cell, c(4.5, 3, 2))
  expect_identical(w[[1L]]$wanted, c(67L, 100L, 150L))
  expect_identical(w[[1L]]$used, 50L)
  # Sizes given as numbers that each give fewer than two cells: one warning too.
  r <- .tm_catch(resolution_profile(.tm_pts()[1:20, ], "z", obs_per_cell = c(40, 18, 4.5),
                                    estimate = "mean"))
  w <- .tm_of(r, "spatialkit_size_not_reached")
  expect_length(w, 1L)
  expect_identical(conditionMessage(w[[1L]]), paste0(
    "resolution_profile(): obs_per_cell = 40 and 18 on 20 rows are each fewer than two cells; ",
    "two cells are used for each (10 rows per cell)."))
  expect_identical(w[[1L]]$obs_per_cell, c(40, 18))
  expect_identical(w[[1L]]$wanted, c(0L, 1L))
  expect_identical(w[[1L]]$reason, "fewer than two cells")
  # Two sizes that ask for the same count: that count, "each".
  r <- .tm_catch(resolution_profile(st, "z", estimate = "mean", obs_per_cell = c(3, 2.99)))
  w <- .tm_of(r, "spatialkit_size_not_reached")
  expect_length(w, 1L)
  expect_identical(conditionMessage(w[[1L]]), paste0(
    "resolution_profile(): obs_per_cell = 3 and 2.99 ask for 100 cells each, more than k-means ",
    "can place on the layer's 50 distinct locations; 50 cells are used for each (6 rows per ",
    "cell)."))
  expect_identical(w[[1L]]$wanted, c(100L, 100L))
})


# ---- 8. Older texts ---------------------------------------------------------------------

test_that("at the support ceiling the hint names `levels` past it", {
  s <- .tm_quiet(select_resolution(resolution_profile(.tm_pts(), "z", partition = "density",
                                                      estimate = NULL, obs_per_cell = NULL,
                                                      sac = NULL, nstart = 25L), "cp"))
  expect_identical(s$edge, "the support ceiling (n / min_cell_n)")
  out <- paste(sub("^ +(note +: )?", "", .tm_print(s)), collapse = " ")
  expect_match(out, paste0("Pass `levels` that run past this count to see whether the ",
                           "criterion keeps going; a lower min_cell_n raises the ceiling too, ",
                           "but moves every level below it."), fixed = TRUE)
  expect_false(grepl("Lower min_cell_n", out, fixed = TRUE))
})

test_that("make_folds() says that a `block_size` given is used beside auto_range", {
  skip_if_not_installed("gstat")
  r <- .tm_catch(make_folds(.tm_pts(), k = 3, method = "block_kfold", response_var = "z",
                            auto_range = TRUE, block_size = 300, seed = 1))
  said <- grep("estimated spatial autocorrelation range", r$messages, value = TRUE)
  expect_length(said, 1L)
  # The size given is a minimum edge, as the range is when it sets the
  # blocks, so the message says what is built from it.
  expect_match(said, paste0("CRS units; the `block_size` given \\(300\\) is used as minimum ",
                            "block size, so the range is only compared against it\\. Over the ",
                            ".* extent that is a ([0-9]+) x ([0-9]+) grid of blocks ",
                            "([0-9.]+) wide and ([0-9.]+) tall\\.\n$"))
  dims <- as.numeric(regmatches(said, regexec(paste0("a ([0-9]+) x ([0-9]+) grid of blocks ",
                                                     "([0-9.]+) wide and ([0-9.]+) tall"),
                                              said))[[1L]][-1L])
  bb <- sf::st_bbox(.tm_pts())
  expect_identical(dims[1:2], c(floor((bb[["xmax"]] - bb[["xmin"]]) / 300),
                                floor((bb[["ymax"]] - bb[["ymin"]]) / 300)))
  expect_gte(dims[3], 300)
  expect_gte(dims[4], 300)
  expect_false(grepl("using as minimum block size", said, fixed = TRUE))
})

test_that("make_folds() on one row says 1 fold and no way to fewer folds", {
  r <- .tm_catch(make_folds(.tm_pts()[1, ], k = 5, method = "random_kfold", seed = 1))
  w <- .tm_of(r, "spatialkit_k_reduced")
  expect_length(w, 1L)
  expect_identical(conditionMessage(w[[1L]]), paste0(
    "make_folds(random_kfold): returning 1 fold, not the k = 5 asked for: the layer has 1 row ",
    "with usable coordinates, so the fold holds it out and leaves no row to fit on."))
  expect_identical(w[[1L]]$k, 1L)
})

test_that("a random seeding given a density setting is told so in the words of the method", {
  r <- .tm_catch(get_voronoi_seeds(.tm_layers()$bnd, method = "random", n = 5, set_seed = 3,
                                   density_power = 0.5, partition = "kmeans"))
  expect_identical(r$text, paste0("get_voronoi_seeds(): `density_power` is read by method = ",
                                  "\"kmeans\" only; method = \"random\" ignores it."))
  # A k-means seeding keeps the words of the kind.
  r <- .tm_catch(get_voronoi_seeds(method = "kmeans", n = 5, sample_points = .tm_pts(),
                                   set_seed = 3, density_power = 0.5, partition = "kmeans"))
  expect_true(any(grepl("is read by partition = \"density\" only, and partition = \"kmeans\" ignores it",
                        r$text, fixed = TRUE)))
})

test_that("model_metrics() names cv_bayes() only where draws were asked for", {
  d <- surf_test_points(80)
  fit <- lm_spatial_fit(d, "z", "w")
  e <- .tm_error(model_metrics(fit, newdata = d, se.fit = TRUE))
  expect_true(endsWith(conditionMessage(e), "Call predict() yourself for those."))
  expect_false(grepl("cv_bayes", conditionMessage(e), fixed = TRUE))
  e <- .tm_error(model_metrics(fit, newdata = d, draws = TRUE))
  expect_match(conditionMessage(e),
               "Call predict() yourself for those, or cv_bayes() for CRPS and interval coverage.",
               fixed = TRUE)
  # A Bayesian fit is pointed to cv_bayes() for its standard errors too.  The
  # refusal reads only the class and the arguments, so a bare object of the
  # class reaches it.
  bayes <- structure(list(), class = c("bayesian_fit", "spatial_fit"))
  e <- .tm_error(model_metrics(bayes, newdata = d, se.fit = TRUE))
  expect_match(conditionMessage(e),
               "Call predict() yourself for those, or cv_bayes() for CRPS and interval coverage.",
               fixed = TRUE)
  # The standard errors of a GWR are not.
  skip_if_not_installed("sp")
  g <- withr::with_seed(4, {
    gd <- sf::st_as_sf(data.frame(x = 5e5 + stats::runif(150, 0, 1000),
                                  y = 5e6 + stats::runif(150, 0, 1000), a = stats::rnorm(150)),
                       coords = c("x", "y"), crs = 32632)
    gd$z <- 1 + 2 * gd$a + stats::rnorm(150, 0, 0.3)
    gd
  })
  gw <- .tm_quiet(fit_gwr_model(g, "z", "a", bandwidth = 30))
  e <- .tm_error(model_metrics(gw, newdata = g, se.fit = TRUE))
  expect_true(endsWith(conditionMessage(e), "Call predict() yourself for those."))
  expect_false(grepl("cv_bayes", conditionMessage(e), fixed = TRUE))
})

test_that("the reason for no variogram model names the lag bins, not singular fits", {
  skip_if_not_installed("gstat")
  s <- testthat::with_mocked_bindings(
    .tm_quiet(estimate_sac_range(.tm_pts(), "z")),
    .sac_fit_model = function(...) NULL)
  expect_identical(attr(s, "rejected_reason"),
                   "no variogram model could be fitted: fewer than three lag bins hold point pairs")
})

test_that("the REML detrending on fewer than 30 locations is said not to be fitted", {
  skip_if_not_installed("gstat")
  skip_if_not_installed("nlme")
  d <- .tm_pts()[rep(1:20, 2), ]
  d$z <- d$z + withr::with_seed(2, stats::rnorm(40, 0, 0.01))
  r <- .tm_catch(estimate_sac_range(d, "z", predictor_vars = "noise", detrend = "reml"))
  said <- grep("REML detrending", r$text, value = TRUE, fixed = TRUE)
  expect_length(said, 1L)
  expect_match(said, paste0("the REML detrending on noise was not fitted: the 40 rows lie at 20 ",
                            "locations, fewer than the 30 it needs; falling back to OLS"),
               fixed = TRUE)
  # A fit that returns nothing on enough locations did not converge.
  r <- testthat::with_mocked_bindings(
    .tm_catch(estimate_sac_range(.tm_pts(), "z", predictor_vars = "noise", detrend = "reml")),
    .reml_trend = function(...) NULL)
  said <- grep("REML detrending", r$text, value = TRUE, fixed = TRUE)
  expect_length(said, 1L)
  expect_match(said, "the REML detrending on noise did not converge; falling back to OLS",
               fixed = TRUE)
  # 30 locations are enough: a fit that returns nothing there did not converge.
  d30 <- .tm_pts()[rep(1:30, 2), ]
  d30$z <- d30$z + withr::with_seed(2, stats::rnorm(60, 0, 0.01))
  r <- testthat::with_mocked_bindings(
    .tm_catch(estimate_sac_range(d30, "z", predictor_vars = "noise", detrend = "reml")),
    .reml_trend = function(...) NULL)
  said <- grep("REML detrending", r$text, value = TRUE, fixed = TRUE)
  expect_length(said, 1L)
  expect_match(said, "the REML detrending on noise did not converge; falling back to OLS",
               fixed = TRUE)
  # Without nlme nothing was fitted, whatever the locations.
  r <- testthat::with_mocked_bindings(
    .tm_catch(estimate_sac_range(.tm_pts(), "z", predictor_vars = "noise", detrend = "reml")),
    .reml_trend = function(...) NULL,
    .has_backend = function(p) if (identical(p, "nlme")) FALSE else requireNamespace(p, quietly = TRUE))
  said <- grep("REML detrending", r$text, value = TRUE, fixed = TRUE)
  expect_length(said, 1L)
  expect_match(said, paste0("the REML detrending on noise was not fitted (package 'nlme' is not ",
                            "installed); falling back to OLS"), fixed = TRUE)
})

test_that("kriging_adequacy()'s print says why a pure nugget still has a gain", {
  skip_if_not_installed("gstat")
  units <- .tm_layers()$units
  asg <- .tm_quiet(assign_features_to_polygons(.tm_pts(), units, polygon_id_col = "GEOID"))
  ka <- .tm_quiet(kriging_adequacy(asg, "noise", units, id_col = "GEOID"))
  expect_true(all(attr(ka, "variogram")$model == "Nug"))
  # What the line says: each cell is kriged from the 50 locations nearest its
  # centre, its own among them, and under one mean for every place the gain
  # is that of averaging those over averaging its own.
  df <- sf::st_drop_geometry(ka)
  expect_true(all(df$kr_n_used == 50L))
  expect_equal(df$kr_gain, 1 - df$n / df$kr_n_used, tolerance = 1e-3)
  line <- grep("^  kr_gain", .tm_print(ka), value = TRUE)
  expect_length(line, 1L)
  expect_match(line, paste0("; under a pure nugget the variogram gives every place one mean, so ",
                            "this is the gain of averaging the kr_n_used locations a cell is ",
                            "kriged from, its own among them, over averaging its own alone ",
                            "\\(the nmax = 50 nearest its centre\\)$"))
  # A variogram with structure gets no such clause.
  kz <- .tm_quiet(kriging_adequacy(asg, "z", units, id_col = "GEOID"))
  expect_false(all(attr(kz, "variogram")$model == "Nug"))
  out <- .tm_print(kz)
  expect_length(grep("^  kr_gain", out), 1L)
  expect_false(any(grepl("pure nugget", out, fixed = TRUE)))
})

test_that("under a pure nugget the kriging print says that a cell of more than nmax locations is kriged from all of them", {
  skip_if_not_installed("gstat")
  # 300 points; one cell of 600 by 1000 m beside ten squares of 200 m, with
  # nmax = 30.  The large cell's own locations reach beyond the 30 nearest its
  # centre, so it is kriged from all of its own plus the 30 nearest outside it.
  pts <- withr::with_seed(12, sf::st_as_sf(
    data.frame(x = stats::runif(300, 0, 1000), y = stats::runif(300, 0, 1000),
               z = stats::rnorm(300)), coords = c("x", "y"), crs = 32632))
  sq <- function(x0, x1, y0, y1)
    sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0))))
  small <- sf::st_make_grid(sf::st_sfc(sq(600, 1000, 0, 1000), crs = 32632), cellsize = 200)
  cells <- sf::st_sf(poly_id = seq_len(length(small) + 1L),
                     geometry = c(sf::st_sfc(sq(0, 600, 0, 1000), crs = 32632), small))
  asg <- .tm_quiet(assign_features_to_polygons(pts, cells))
  nug <- structure(7500, class = c("sac_range", "numeric"),
                   variogram_model = gstat::vgm(psill = 1, model = "Nug", range = 0),
                   crs = sf::st_crs(32632))
  ka <- .tm_quiet(kriging_adequacy(asg, "z", cells, sac = nug, k = 3, nmax = 30))
  df <- sf::st_drop_geometry(ka)
  expect_equal(df$n[1L], 180)
  expect_equal(df$kr_n_used[1L], 210)
  expect_true(all(df$kr_n_used[-1L] == 30L))
  expect_equal(df$kr_gain, 1 - df$n / df$kr_n_used, tolerance = 1e-3)
  line <- grep("^  kr_gain", .tm_print(ka), value = TRUE)
  expect_length(line, 1L)
  expect_match(line, paste0("over averaging its own alone \\(the nmax = 30 nearest its centre, ",
                            "or, for 1 cell\\(s\\) with own locations beyond those, all of its ",
                            "own plus the 30 nearest outside it\\)$"))
  expect_false(grepl("in place of", line, fixed = TRUE))
})

test_that("the numbers of one printed line share one notation", {
  expect_identical(spatialkit:::.one_notation(c(9.02e7, 1.09e8)), c("9.02e+07", "1.09e+08"))
  # Lines whose numbers share a notation print as they did.
  expect_identical(spatialkit:::.one_notation(c(440, 1.1)), c("440", "1.1"))
  expect_identical(spatialkit:::.one_notation(c(9820000, 10100000)), c("9820000", "10100000"))
  expect_identical(spatialkit:::.one_notation(c(1.44e8, 1.44e8)), c("1.44e+08", "1.44e+08"))
  expect_identical(spatialkit:::.one_notation(c(NA, 1.2e8)), c("NA", "1.2e+08"))
  # Rewritten numbers keep no trailing zeros, as format(signif()) writes
  # every other number of these prints.
  expect_identical(spatialkit:::.one_notation(c(90200000, 1.1e8)), c("9.02e+07", "1.1e+08"))
  expect_identical(spatialkit:::.one_notation(c(1e-5, 0.5)), c("1e-05", "5e-01"))
  big <- .tm_pts(); big$z <- big$z * 9950
  s <- .tm_quiet(select_resolution(resolution_profile(big, "z", estimate = "mean"), "loo"))
  line <- grep("^  loo ", .tm_print(s), value = TRUE)
  expect_length(line, 1L)
  expect_match(line, "^  loo         : [0-9.]+e\\+[0-9]+ \\(one cell for the layer: [0-9.]+e\\+[0-9]+\\)$")
})

test_that("the size criterion's note writes its two errors in the notation of the loo line", {
  # The default path: pooled values at the default size, on a white-noise
  # column scaled so that one error comes out in scientific notation and the
  # other not.
  big <- .tm_pts(); big$noise <- big$noise * 10300
  s <- .tm_quiet(select_resolution(resolution_profile(big, "noise")))
  expect_identical(s$criterion, "size")
  expect_identical(c(format(signif(s$loo, 3)), format(signif(s$loo_one_cell, 3))),
                   c("1.01e+08", "98500000"))
  nums <- spatialkit:::.one_notation(c(s$loo, s$loo_one_cell))
  expect_identical(nums, c("1.01e+08", "9.85e+07"))
  expect_identical(s$note[["structure"]], sprintf(paste0(
    "at %d cells the leave-one-out error of the pooled value (%s) is not below that of one ",
    "mean for the whole layer (%s) by more than 2%%: the response shows no spatial structure ",
    "that cells of this size can use."), s$best, nums[1L], nums[2L]))
  out <- .tm_print(s)
  expect_identical(grep("^  loo ", out, value = TRUE),
                   sprintf("  loo         : %s (one cell for the layer: %s)", nums[1L], nums[2L]))
  expect_false(any(grepl("98500000", out, fixed = TRUE)))
})


# ---- 9. Left by the high bugs ---------------------------------------------------------

test_that("the prints of the default path and of one level or one start say what happened", {
  pts <- .tm_pts()
  pd <- .tm_quiet(resolution_profile(pts, "z"))
  # The default size is called the default, as the selection's note calls it.
  out <- .tm_print(.tm_quiet(select_resolution(pd)))
  expect_identical(out[2L],
                   "  size        : 4.5 rows per cell, the default (obs_per_cell); 400 rows / 89 cells = 4.49")
  # One level: one count in the summary's header.
  out <- .tm_print(summary(pd))
  expect_identical(out[1L], "Resolution picks: 1 criterion over 1 level (89 cells)")
  # A selection on one level: "1 of 1 level".
  s <- select_resolution(.tm_quiet(resolution_profile(pts[1:5, ], "z", estimate = NULL,
                                                      obs_per_cell = NULL, sac = NULL,
                                                      partition = "kmeans")), "cp")
  expect_identical(.tm_print(s)[2L], "  flat region : 2 (1 of 1 level)")
  # One start: "1 k-means++ restart".
  p1 <- .tm_quiet(resolution_profile(pts, "z", partition = "kmeans", nstart = 1))
  expect_true(any(grepl("; 1 k-means++ restart per level", .tm_print(p1), fixed = TRUE)))
})

test_that("a summary whose picks agree names one count", {
  # 60 stations visited five times, with a trend: one size, and two criteria
  # that read it.
  v <- withr::with_seed(23, {
    s <- data.frame(x = stats::runif(60, 0, 8e3), y = stats::runif(60, 0, 8e3))
    d <- s[rep(1:60, each = 5), ]
    d$z <- d$x / 8e3 + stats::rnorm(300, 0, 0.5)
    sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  })
  sm <- .tm_quiet(summary(.tm_quiet(resolution_profile(v, "z"))))
  expect_gte(nrow(sm), 2L)
  expect_length(unique(sm$best), 1L)
  out <- .tm_print(sm)
  expect_true(sprintf("  every pick is %d cells", sm$best[1L]) %in% out)
  expect_false(any(grepl("picks span", out, fixed = TRUE)))
  # Picks that differ keep their span.
  pk <- .tm_quiet(resolution_profile(.tm_pts(), "z", estimate = NULL, obs_per_cell = NULL,
                                     partition = "kmeans", nstart = 3L))
  sk <- .tm_quiet(summary(pk))
  expect_gt(length(unique(sk$best)), 1L)
  rng <- range(sk$best)
  expect_true(sprintf("  picks span %d to %d cells (%.1fx)", rng[1L], rng[2L], rng[2L] / rng[1L])
              %in% .tm_print(sk))
})

test_that("summary() says that the F test failed where the cells beat one mean", {
  # 60 points with a weak trend (seed 2): the pooled loo beats one mean by
  # more than 2 percent, and the F test does not tell the cells from one mean.
  p <- withr::with_seed(2, {
    d <- data.frame(x = stats::runif(60, 0, 1e4), y = stats::runif(60, 0, 1e4))
    d$z <- d$x / 1e4 * 0.8 + stats::rnorm(60)
    sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  })
  pr <- .tm_quiet(resolution_profile(p, "z"))
  s <- .tm_quiet(select_resolution(pr, "loo"))
  # That it is the case: one mean is worse by more than the margin, and p >= 0.05.
  expect_false(s$structured)
  expect_gt(s$one_cell, s$value * 1.02)
  expect_gte(s$structure_p, 0.05)
  sm <- summary(pr, criteria = "loo")
  expect_identical(sm$structured, FALSE)
  expect_identical(attr(sm, "not_identified"),
                   c(loo = sprintf("the F test does not tell the cell means from one mean (p = %s)",
                                   format(signif(s$structure_p, 2)))))
  out <- .tm_print(sm)
  expect_false(any(grepl("does as well", out, fixed = TRUE)))
  said <- strwrap(sprintf(paste0("loo: the cells beat one cell for the whole layer, but the F ",
                                 "test does not tell the cell means from one mean (p = %s), so ",
                                 "no spatial structure the cells can use is shown and this count ",
                                 "is not identified."), format(signif(s$structure_p, 2))),
                  width = 78, prefix = "  ", exdent = 2)
  expect_identical(utils::tail(out, length(said)), said)
  # Where one mean does as well, the summary has no such attribute and keeps its words.
  noise <- .tm_quiet(resolution_profile(.tm_pts(), "noise", estimate = "mean"))
  sm <- summary(noise, criteria = "loo")
  expect_identical(sm$structured, FALSE)
  expect_null(attr(sm, "not_identified"))
  expect_true(any(grepl("loo: one cell for the whole layer does as well", .tm_print(sm),
                        fixed = TRUE)))
})

test_that("the seeding functions refuse the NA of determine_optimal_levels() with its note", {
  u <- withr::with_seed(9, sf::st_as_sf(data.frame(x = stats::runif(300), y = stats::runif(300)),
                                        coords = c("x", "y"), crs = 32632))
  k <- .tm_quiet(determine_optimal_levels(u))
  expect_true(is.na(k))
  expect_false(attr(k, "structured"))
  note <- paste0(" It is what determine_optimal_levels() returns at its default, on_no_elbow = ",
                 "\"na\", when the points name no count: the WSS curve has no elbow, or there ",
                 "is no ladder long enough to show one. Pass a number.")
  expect_identical(conditionMessage(.tm_error(voronoi_seeds_kmeans(u, k))),
                   paste0("voronoi_seeds_kmeans(): `k` must be a single positive number of seeds; ",
                          "got NA.", note))
  bnd <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(u)))
  expect_identical(conditionMessage(.tm_error(voronoi_seeds_random(bnd, k))),
                   paste0("voronoi_seeds_random(): `k` must be a single positive number of seeds; ",
                          "got NA.", note))
  # get_voronoi_seeds() gives the same note, as it did.
  expect_true(endsWith(conditionMessage(.tm_error(get_voronoi_seeds(method = "kmeans", n = k,
                                                                    sample_points = u))), note))
  # An NA of no such origin keeps the present words.
  expect_identical(conditionMessage(.tm_error(voronoi_seeds_kmeans(u, NA_real_))),
                   "voronoi_seeds_kmeans(): `k` must be a single positive number of seeds; got NA.")
})
