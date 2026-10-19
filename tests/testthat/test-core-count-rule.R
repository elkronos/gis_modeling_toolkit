# tests/testthat/test-core-count-rule.R
# ---------------------------------------------------------------------------
# The record of how a cell count was set, from the selection to the summary:
# where it sits, what `scored_cells` says beside each layer of seeds, and the
# two checks summarize_by_cell() runs against it (the estimate; the cells
# and the rows), with every case in which the second is skipped.
#
# Every profile names `partition`, `nstart`, `estimate`, `obs_per_cell` and
# `sac`, every selection its `criterion`, every summary its `estimate`, and
# every seeding that places seeds of its own its kind and its starts.
# Unless a test says otherwise the regions are of the density kind from 3
# k-means++ starts under seed 123.  A layer is drawn under the seed its test
# states: the locations (x, then y), then the response.
# ---------------------------------------------------------------------------

# ---- Layers and helpers ------------------------------------------------------

.cr_seed <- function(seed, env = parent.frame())
  withr::local_seed(seed, .local_envir = env, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

# `n` locations uniform on a square of side 1e5 with its corner at (5e5,
# 5e6) in EPSG:32617, with the field of the engine's tests, sin(x / 2e4) +
# cos(y / 1.5e4) plus normal noise of standard deviation 1.5.
.cr_layer <- function(n, seed) {
  d <- local({
    .cr_seed(seed)
    x <- stats::runif(n, 0, 1e5); y <- stats::runif(n, 0, 1e5)
    data.frame(x = 5e5 + x, y = 5e6 + y,
               z = sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(n, 0, 1.5))
  })
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32617)
}

.cr_profile <- function(pts, ..., response = "z", partition = "density", nstart = 3,
                        sac = "auto", estimate = "pooled", obs_per_cell = 4.5, seed = 123)
  resolution_profile(pts, response, ..., partition = partition, nstart = nstart, sac = sac,
                     estimate = estimate, obs_per_cell = obs_per_cell, seed = seed)

# A call with its warnings kept and muffled.
.cr_catch <- function(expr) {
  got <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    got[[length(got) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = got)
}
.cr_classes <- function(r) vapply(r$warnings, function(w) class(w)[1L], character(1))
.cr_of <- function(r, class) Filter(function(w) inherits(w, class), r$warnings)
.cr_fields <- function(w) unclass(w)[setdiff(names(unclass(w)), c("message", "call"))]
.cr_print <- function(x) utils::capture.output(print(x))

.cr_store <- new.env(parent = emptyenv())
.cr_keep <- function(key, expr) {
  if (is.null(.cr_store[[key]])) .cr_store[[key]] <- expr
  .cr_store[[key]]
}
.cr_sizes <- c(18, 12, 9, 6, 4.5, 3, 2)
# The fields of the record, in order.
.cr_names <- c("version", "rule", "estimate", "cells", "obs_per_cell", "obs_per_cell_asked",
               "response", "n", "n_points", "n_fit", "seed", "features", "scored",
               "partition", "crs", "scored_cells")

# 600 rows (seed 5101).  Its pooled profile at 4.5 rows a cell has 133
# cells; its profile of plain means on the seven sizes has its lowest
# leave-one-out error at 50 cells, a level inside the ladder.  Neither path
# raises a warning from the profile to the summary.
.cr_600    <- function() .cr_keep("l600", .cr_layer(600, 5101))
.cr_pooled <- function() .cr_keep("pooled", .cr_profile(.cr_600()))
.cr_mean7  <- function()
  .cr_keep("mean7", .cr_profile(.cr_600(), estimate = "mean", obs_per_cell = .cr_sizes))
.cr_size <- function() .cr_keep("size", select_resolution(.cr_pooled(), criterion = "size"))
.cr_loo  <- function() .cr_keep("loo", select_resolution(.cr_mean7(), criterion = "loo"))

# The calls from a count to the assigned layer.
.cr_cells <- function(pts, n, ...) {
  seeds <- get_voronoi_seeds(method = "kmeans", n = n, sample_points = pts, ...)
  tess  <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)
  asg   <- assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id",
                                       quiet = TRUE)
  list(seeds = seeds, tess = tess, asg = asg)
}
.cr_A <- function() .cr_keep("A", .cr_cells(.cr_600(), .cr_size()))   # sized for pooled values
.cr_B <- function() .cr_keep("B", .cr_cells(.cr_600(), .cr_loo()))    # chosen for plain means
# The same layer with no response east of x = 550,000: the centres are fitted
# to all 600 points and the size divides the 299 rows with a response, so
# the count (66) is not the number of cells that hold a scored row (37 of
# them on the machine this was written on; about half on any).
.cr_gap <- function() .cr_keep("gap", {
  pts <- .cr_600()
  pts$z[sf::st_coordinates(pts)[, 1L] > 5.5e5] <- NA
  sel <- select_resolution(suppressWarnings(.cr_profile(pts)), criterion = "size")
  c(list(pts = pts, sel = sel), .cr_cells(pts, sel))
})
.cr_sum <- function(asg, cells, estimate, response = "z", ...)
  .cr_catch(summarize_by_cell(asg, response, id_col = "cell_id", cells_sf = cells,
                              estimate = estimate, ...))
.cr_check <- function(r) attr(r$value, "count_rule")$check
.cr_bnd <- function(pts) sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(pts)))


# ---- Where no record is written ------------------------------------------------

test_that("a count given as a number, or read by a present criterion off a profile without an estimate, carries no record", {
  pts <- .cr_600()
  bare <- function(x) is.null(attr(x, "count_rule"))
  none <- function(P) {
    expect_true(bare(P$seeds))
    expect_false("count_rule" %in% names(P$tess$params))
    expect_true(bare(P$tess$cells))
    expect_true(bare(P$asg))
    for (est in c("mean", "pooled")) {
      r <- suppressWarnings(.cr_sum(P$asg, P$tess$cells, est))
      expect_true(bare(r$value))
      expect_length(.cr_of(r, "spatialkit_cells_differ"), 0L)
      expect_length(.cr_of(r, "spatialkit_estimate_mismatch"), 0L)
    }
  }
  # A number.
  none(.cr_cells(pts, 40, set_seed = 1, partition = "kmeans", kmeans_nstart = 10))
  expect_true(bare(voronoi_seeds_kmeans(pts, 40, partition = "kmeans", nstart = 10)))
  expect_true(bare(get_voronoi_seeds(.cr_bnd(pts), method = "random", n = 40, set_seed = 1)))
  hex <- build_tessellation(pts, boundary = .cr_bnd(pts), method = "hex", approx_n_cells = 40,
                            quiet = TRUE)
  expect_false("count_rule" %in% names(hex$params))
  expect_true(bare(hex$cells))
  # A profile as it is built without the new arguments, read by cp: as a
  # selection, and as the count itself.
  plain <- suppressWarnings(
    .cr_profile(pts, n_levels = 6L, partition = "kmeans", nstart = 5, sac = 3e4,
                estimate = NULL, obs_per_cell = NULL))
  sel <- select_resolution(plain, criterion = "cp")
  expect_null(sel$count_rule)
  none(suppressWarnings(.cr_cells(pts, sel)))
  none(suppressWarnings(.cr_cells(pts, plain)))
  expect_true(bare(suppressWarnings(voronoi_seeds_kmeans(pts, sel))))
  hex <- suppressWarnings(build_tessellation(pts, boundary = .cr_bnd(pts), method = "hex",
                                             approx_n_cells = sel, quiet = TRUE))
  expect_false("count_rule" %in% names(hex$params))
  expect_true(bare(hex$cells))
})

test_that("cells the caller brings have no record and are not checked", {
  pts <- .cr_600()
  units <- sf::st_sf(GEOID = sprintf("u%02d", 1:16),
                     geometry = sf::st_make_grid(.cr_bnd(pts), n = c(4, 4)))
  asg <- assign_features_to_polygons(pts, units, polygon_id_col = "GEOID", quiet = TRUE)
  expect_null(attr(asg, "count_rule"))
  for (est in c("mean", "pooled")) {
    r <- .cr_catch(summarize_by_cell(asg, "z", id_col = "GEOID", cells_sf = units,
                                     estimate = est))
    expect_null(attr(r$value, "count_rule"))
    expect_length(.cr_of(r, "spatialkit_cells_differ"), 0L)
    expect_length(.cr_of(r, "spatialkit_estimate_mismatch"), 0L)
  }
})


# ---- Where the record sits ------------------------------------------------------

test_that("the record travels from the selection to the summary, and says beside the seeds that they are the scored centres", {
  pts <- .cr_600(); sel <- .cr_size(); A <- .cr_A()
  rec <- sel$count_rule
  expect_identical(names(rec), .cr_names)
  expect_s3_class(rec, "spatialkit_count_rule")
  # Beside no centres in the selection; beside the scored centres on its seeds.
  expect_identical(rec$scored_cells, NA)
  expect_identical(attr(sel$seeds, "count_rule")$scored_cells, TRUE)
  on_seeds <- rec; on_seeds$scored_cells <- TRUE
  expect_identical(attr(sel$seeds, "count_rule"), on_seeds)

  expect_identical(attr(A$seeds, "count_rule"), on_seeds)
  expect_identical(A$tess$params$count_rule, on_seeds)
  expect_identical(attr(A$tess$cells, "count_rule"), on_seeds)
  expect_identical(attr(A$asg, "count_rule"), on_seeds)
  # The seeds, the cells and the assigned rows are otherwise those of the
  # same selection without a record.
  no_rec <- sel; no_rec$count_rule <- NULL
  attr(no_rec$seeds, "count_rule") <- NULL
  N <- .cr_cells(pts, no_rec)
  expect_null(attr(N$seeds, "count_rule"))
  expect_false("count_rule" %in% names(N$tess$params))
  without <- function(x) { attr(x, "count_rule") <- NULL; x }
  expect_identical(without(A$seeds), N$seeds)
  expect_identical(without(A$tess$cells), N$tess$cells)
  expect_identical(without(A$asg), N$asg)
  expect_identical(A$tess$params[names(N$tess$params)], N$tess$params)
  expect_identical(setdiff(names(A$tess$params), names(N$tess$params)), "count_rule")
  expect_identical(utils::tail(names(A$tess$params), 2L), c("partition", "count_rule"))

  r <- .cr_sum(A$asg, A$tess$cells, "pooled")
  expect_identical(.cr_classes(r), character(0))
  got <- attr(r$value, "count_rule")
  expect_identical(names(got), c(.cr_names, "check"))
  expect_s3_class(got, "spatialkit_count_rule")
  expect_identical(got$check, list(estimate = "match", cells = "match"))
  got$check <- NULL
  expect_identical(got, on_seeds)
  # A frame without `cells_sf` carries it too.
  r <- .cr_catch(summarize_by_cell(A$asg, "z", id_col = "cell_id", estimate = "pooled"))
  expect_identical(attr(r$value, "count_rule")$check, list(estimate = "match", cells = "match"))
  expect_false(inherits(r$value, "sf"))

  # voronoi_seeds_kmeans() hands on the same centres with the same record.
  vk <- voronoi_seeds_kmeans(pts, sel)
  expect_identical(attr(vk, "count_rule"), on_seeds)
  expect_identical(sf::st_geometry(vk), sf::st_geometry(A$seeds))
  # create_voronoi_polygons() writes it where build_tessellation() does.
  vor <- create_voronoi_polygons(A$seeds, quiet = TRUE)
  expect_identical(vor$params$count_rule, on_seeds)
  expect_identical(attr(vor$cells, "count_rule"), on_seeds)
  # A profile with an estimate passed as the count is read by "auto", and
  # its seeds carry the record of that reading.
  P <- .cr_cells(pts, .cr_pooled())
  expect_identical(attr(P$seeds, "count_rule"), on_seeds)
  expect_identical(.cr_check(.cr_sum(P$asg, P$tess$cells, "pooled")),
                   list(estimate = "match", cells = "match"))
})

test_that("the record stays on the seeds and the cells through a boundary and a change of CRS", {
  pts <- .cr_600(); sel <- .cr_size()
  on_seeds <- attr(sel$seeds, "count_rule")
  bnd <- .cr_bnd(pts)
  # A boundary in another CRS: the seeds come back in it.
  ll <- get_voronoi_seeds(sf::st_transform(bnd, 4326), method = "kmeans", n = sel)
  expect_true(sf::st_is_longlat(ll))
  expect_identical(attr(ll, "count_rule"), on_seeds)
  expect_identical(attr(voronoi_seeds_kmeans(sf::st_transform(pts, 4326), sel), "count_rule"),
                   on_seeds)
  # Cells returned in a geographic CRS.
  tess <- build_tessellation(.cr_A()$seeds, method = "voronoi", crs = 4326, quiet = TRUE)
  expect_true(sf::st_is_longlat(tess$cells))
  expect_identical(attr(tess$cells, "count_rule"), on_seeds)
  expect_identical(tess$params$count_rule, on_seeds)
  # The record names the CRS the profile worked in, whatever the seeds are in.
  expect_identical(attr(ll, "count_rule")$crs, sf::st_crs(pts))
})


# ---- scored_cells ---------------------------------------------------------------

test_that("seeds a function places itself at the count are marked as not the scored centres", {
  pts <- .cr_600(); sel <- .cr_size(); bnd <- .cr_bnd(pts)
  rec_false <- sel$count_rule; rec_false$scored_cells <- FALSE
  # A selection whose seeds were removed: a partition of the function's own.
  bare_sel <- sel; bare_sel$seeds <- NULL
  fresh <- get_voronoi_seeds(method = "kmeans", n = bare_sel, sample_points = pts, set_seed = 1,
                             partition = "kmeans", kmeans_nstart = 10)
  expect_identical(nrow(fresh), 133L)
  expect_identical(attr(fresh, "count_rule"), rec_false)
  expect_identical(attr(voronoi_seeds_kmeans(pts, bare_sel, partition = "kmeans", nstart = 10),
                        "count_rule"), rec_false)
  # Seeds that are not at the count are no more the scored centres.
  short <- sel; short$seeds <- sel$seeds[1:50, ]
  expect_identical(attr(get_voronoi_seeds(method = "kmeans", n = short, sample_points = pts,
                                          set_seed = 1, partition = "kmeans",
                                          kmeans_nstart = 10), "count_rule"), rec_false)
  # method = "random" takes the count and places its own seeds.
  rnd <- get_voronoi_seeds(bnd, method = "random", n = sel, set_seed = 3)
  expect_identical(nrow(rnd), 133L)
  expect_identical(attr(rnd, "count_rule"), rec_false)
  # The cells built on such seeds say the same.
  tess <- build_tessellation(fresh, method = "voronoi", quiet = TRUE)
  expect_identical(tess$params$count_rule, rec_false)
  expect_identical(attr(tess$cells, "count_rule"), rec_false)
})

test_that("method = \"provided\" writes no record, and hands on the one its seeds carry", {
  pts <- .cr_600(); sel <- .cr_size()
  r <- .cr_catch(get_voronoi_seeds(method = "provided", n = sel, seeds = pts[1:20, ]))
  expect_null(attr(r$value, "count_rule"))
  expect_match(conditionMessage(r$warnings[[1L]]), "method = 'provided' ignores 'n'", fixed = TRUE)
  expect_null(attr(get_voronoi_seeds(method = "provided", seeds = pts[1:20, ]), "count_rule"))
  # The selection's own seeds, whole: still the scored centres.
  whole <- get_voronoi_seeds(method = "provided", seeds = sel$seeds)
  expect_identical(attr(whole, "count_rule"), attr(sel$seeds, "count_rule"))
  # Twenty of them carry a record of 133 cells, and are not those centres.
  part <- get_voronoi_seeds(method = "provided", seeds = sel$seeds[1:20, ])
  expect_identical(attr(sel$seeds[1:20, ], "count_rule")$scored_cells, TRUE)
  expect_identical(attr(part, "count_rule")$scored_cells, FALSE)
  expect_identical(attr(part, "count_rule")$cells, 133L)
  # More of them than the count are no more those centres than fewer are.
  more <- rbind(sel$seeds, sel$seeds[1:5, ])
  attr(more, "count_rule") <- attr(sel$seeds, "count_rule")
  expect_identical(nrow(more), 138L)
  expect_identical(attr(get_voronoi_seeds(method = "provided", seeds = more),
                        "count_rule")$scored_cells, FALSE)
  # The seeds are counted against the count, not against the cells that
  # hold a scored row: all 66 centres of a layer on which about half hold one.
  G <- .cr_gap()
  expect_identical(nrow(G$sel$seeds), 66L)
  expect_identical(attr(get_voronoi_seeds(method = "provided", seeds = G$sel$seeds),
                        "count_rule")$scored_cells, TRUE)
})

test_that("a lattice sized by a selection carries the record as not the scored cells; a cellsize sizes it alone", {
  pts <- .cr_600(); sel <- .cr_size(); bnd <- .cr_bnd(pts)
  rec_false <- sel$count_rule; rec_false$scored_cells <- FALSE
  for (method in c("hex", "square")) {
    grid <- suppressWarnings(build_tessellation(pts, boundary = bnd, method = method,
                                                approx_n_cells = sel, quiet = TRUE))
    expect_identical(grid$params$count_rule, rec_false, info = method)
    expect_identical(attr(grid$cells, "count_rule"), rec_false, info = method)
    expect_identical(utils::tail(names(grid$params), 1L), "count_rule")
  }
  sized <- suppressWarnings(build_tessellation(pts, boundary = bnd, method = "hex",
                                               approx_n_cells = sel, cellsize = 9000,
                                               quiet = TRUE))
  expect_false("count_rule" %in% names(sized$params))
  expect_null(attr(sized$cells, "count_rule"))
  # A Voronoi tessellation does not read the count, and takes no record from it.
  vor <- suppressWarnings(build_tessellation(pts, method = "voronoi", approx_n_cells = sel,
                                             quiet = TRUE))
  expect_false("count_rule" %in% names(vor$params))
  expect_null(attr(vor$cells, "count_rule"))
})

test_that("the selection's seeds, tessellated as they are, lead to a summary that is checked", {
  pts <- .cr_600(); sel <- .cr_size()
  expect_identical(sel$count_rule$scored_cells, NA)
  expect_identical(attr(sel$seeds, "count_rule")$scored_cells, TRUE)
  tess <- build_tessellation(sel$seeds, method = "voronoi", quiet = TRUE)
  expect_identical(attr(tess$cells, "count_rule")$scored_cells, TRUE)
  asg <- assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id", quiet = TRUE)
  r <- .cr_sum(asg, tess$cells, "pooled")
  expect_identical(.cr_check(r), list(estimate = "match", cells = "match"))
  expect_identical(.cr_classes(r), character(0))
  # Twenty of them are other cells, and the check says so.
  tess <- build_tessellation(sel$seeds[1:20, ], method = "voronoi", quiet = TRUE)
  asg <- assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id", quiet = TRUE)
  r <- .cr_sum(asg, tess$cells, "pooled")
  expect_identical(.cr_check(r)$cells, "differ")
  w <- .cr_of(r, "spatialkit_cells_differ")
  expect_length(w, 1L)
  expect_identical(w[[1L]]$what, "cells")
  expect_identical(w[[1L]]$found$cells, 20L)
})


# ---- The estimate ---------------------------------------------------------------

test_that("a summary of the estimate the count was chosen for is checked without a word", {
  A <- .cr_A(); B <- .cr_B()
  for (case in list(list(A, "pooled"), list(B, "mean"))) {
    r <- .cr_sum(case[[1L]]$asg, case[[1L]]$tess$cells, case[[2L]])
    expect_identical(.cr_classes(r), character(0))
    expect_identical(.cr_check(r), list(estimate = "match", cells = "match"))
  }
  expect_identical(attr(B$asg, "count_rule")$rule, "loo")
  expect_identical(attr(B$asg, "count_rule")$cells, 50L)
})

test_that("plain means on cells sized for pooled values are a mismatch, and the frame is the one asked for", {
  A <- .cr_A()
  r <- .cr_sum(A$asg, A$tess$cells, "mean")
  expect_identical(.cr_classes(r), "spatialkit_estimate_mismatch")
  w <- r$warnings[[1L]]
  expect_identical(conditionMessage(w), paste0(
    "summarize_by_cell(): the number of cells (133, set by size) was chosen for the pooled ",
    "value, and this summary reports the plain mean. At this cell size plain means of small ",
    "cells are noisy: pass estimate = \"pooled\", or choose the count for plain means ",
    "(resolution_profile(estimate = \"mean\"))."))
  expect_identical(.cr_fields(w),
                   list(chosen_for = "pooled", reported = "mean", rule = "size", cells = 133L))
  expect_null(conditionCall(w))
  # The cells and rows are still the ones scored: checked with a fit of the
  # estimate the profile scored.
  expect_identical(.cr_check(r), list(estimate = "mismatch", cells = "match"))
  # The frame is the one the same call returns on the same cells and rows
  # without a record.
  asg <- A$asg; cells <- A$tess$cells
  attr(asg, "count_rule") <- NULL; attr(cells, "count_rule") <- NULL
  plain <- .cr_sum(asg, cells, "mean")
  expect_identical(.cr_classes(plain), character(0))
  expect_null(attr(plain$value, "count_rule"))
  got <- r$value
  expect_identical(setdiff(names(attributes(got)), names(attributes(plain$value))), "count_rule")
  attr(got, "count_rule") <- NULL
  expect_identical(got, plain$value)
  # With the pooled values asked for, the same holds of the pooled frame.
  pooled <- .cr_sum(A$asg, A$tess$cells, "pooled")
  bare   <- .cr_sum(asg, cells, "pooled")
  got <- pooled$value
  attr(got, "count_rule") <- NULL
  expect_identical(got, bare$value)
})

test_that("pooled values on cells chosen for plain means are a mismatch in the other words", {
  B <- .cr_B()
  r <- .cr_sum(B$asg, B$tess$cells, "pooled")
  w <- .cr_of(r, "spatialkit_estimate_mismatch")
  expect_length(w, 1L)
  expect_identical(conditionMessage(w[[1L]]), paste0(
    "summarize_by_cell(): the number of cells (50, set by loo) was chosen for the plain mean, ",
    "and this summary reports the pooled value. Pooled values are valid on these cells; the ",
    "count that suits them is set by a size (criterion = \"size\")."))
  expect_identical(.cr_fields(w[[1L]]),
                   list(chosen_for = "mean", reported = "pooled", rule = "loo", cells = 50L))
  # The record scored plain means, so the check fits those, beside the
  # pooled values the call asked for.
  expect_identical(.cr_check(r), list(estimate = "mismatch", cells = "match"))
  expect_length(.cr_of(r, "spatialkit_cells_differ"), 0L)
})

test_that("where centres hold no scored row the count is the centres, and the record keeps both numbers", {
  G <- .cr_gap()
  rec <- G$sel$count_rule
  expect_identical(rec$cells, 66L)
  occupied <- rec$scored$cells
  expect_gt(occupied, 25L)
  expect_lt(occupied, 50L)
  expect_identical(c(rec$n, rec$scored$n, rec$n_points), c(299L, 299L, 600L))
  expect_identical(nrow(G$tess$cells), 66L)
  # The path is the one scored, under either estimate.
  r <- .cr_sum(G$asg, G$tess$cells, "pooled")
  expect_identical(.cr_classes(r), character(0))
  expect_identical(.cr_check(r), list(estimate = "match", cells = "match"))
  r <- .cr_sum(G$asg, G$tess$cells, "mean")
  expect_identical(.cr_classes(r), "spatialkit_estimate_mismatch")
  expect_identical(.cr_check(r), list(estimate = "mismatch", cells = "match"))
  # The warning names the count, in its text and in its field.
  w <- r$warnings[[1L]]
  expect_match(conditionMessage(w),
               "the number of cells (66, set by size) was chosen for the pooled value",
               fixed = TRUE)
  expect_identical(.cr_fields(w),
                   list(chosen_for = "pooled", reported = "mean", rule = "size", cells = 66L))
  # So does the print, with the cells that hold a scored row on their line.
  out <- .cr_print(attr(r$value, "count_rule"))
  expect_identical(out[1L],
    "Count: 66 cells, set by size (4.5 rows per cell asked; 299 rows / 66 cells = 4.53)")
  expect_match(out[6L], sprintf("on 299 rows in %d cells (one cell: ", occupied), fixed = TRUE)
  # Rows removed are reported against the occupied cells and the rows scored.
  r <- .cr_sum(G$asg[!is.na(G$asg$z), ][-(1:20), ], G$tess$cells, "pooled")
  w <- .cr_of(r, "spatialkit_cells_differ")[[1L]]
  expect_identical(w$what, "n")
  expect_identical(w$scored[c("n", "cells")], list(n = 299L, cells = occupied))
  expect_identical(w$found$n, 279L)
})

test_that("a count cp chose on a pooled profile was chosen for plain means", {
  pts <- .cr_600()
  prof <- suppressWarnings(.cr_profile(pts, levels = c(50, 133, 300), sac = 3e4,
                                       obs_per_cell = NULL))
  sel <- select_resolution(prof, criterion = "cp")
  expect_identical(sel$count_rule$rule, "cp")
  expect_identical(sel$count_rule$estimate, "mean")
  expect_identical(sel$count_rule$scored$estimate, "pooled")
  C <- suppressWarnings(.cr_cells(pts, sel))
  r <- .cr_sum(C$asg, C$tess$cells, "pooled")
  w <- .cr_of(r, "spatialkit_estimate_mismatch")
  expect_length(w, 1L)
  expect_identical(.cr_fields(w[[1L]]),
                   list(chosen_for = "mean", reported = "pooled", rule = "cp",
                        cells = as.integer(sel$best)))
  expect_match(conditionMessage(w[[1L]]),
               sprintf("the number of cells (%d, set by cp) was chosen for the plain mean",
                       sel$best), fixed = TRUE)
  # The profile scored pooled values at that level, and those are compared.
  expect_identical(.cr_check(r), list(estimate = "mismatch", cells = "match"))
  r <- .cr_sum(C$asg, C$tess$cells, "mean")
  expect_identical(.cr_check(r), list(estimate = "match", cells = "match"))
  expect_length(.cr_of(r, "spatialkit_estimate_mismatch"), 0L)
})

test_that("a count chosen for neither estimate has no estimate to match", {
  pts <- .cr_600()
  prof <- suppressWarnings(.cr_profile(pts, levels = c(50, 133, 300), sac = 3e4,
                                       obs_per_cell = NULL))
  sel <- suppressWarnings(select_resolution(prof, criterion = "moran_z", diagnostic = TRUE))
  expect_identical(sel$count_rule$rule, "moran_z")
  expect_identical(sel$count_rule$estimate, NA_character_)
  C <- suppressWarnings(.cr_cells(pts, sel))
  for (est in c("mean", "pooled")) {
    r <- .cr_sum(C$asg, C$tess$cells, est)
    expect_length(.cr_of(r, "spatialkit_estimate_mismatch"), 0L)
    expect_identical(.cr_check(r), list(estimate = NA_character_, cells = "match"))
  }
})


# ---- The cells and the rows -----------------------------------------------------

test_that("rows removed from the assigned layer are reported with what the profile scored and what was found", {
  A <- .cr_A()
  drop <- local({ .cr_seed(11); sample.int(600L, 30L) })       # 5 percent
  r <- .cr_sum(A$asg[-drop, ], A$tess$cells, "pooled")
  w <- .cr_of(r, "spatialkit_cells_differ")
  expect_length(w, 1L)
  w <- w[[1L]]
  rec <- attr(A$asg, "count_rule")
  found <- attr(r$value, "pooled")$columns$resp_z
  expect_identical(.cr_fields(w),
                   list(what = "n",
                        scored = list(n = 600L, cells = 133L, loo = rec$scored$loo),
                        found  = list(n = 570L, cells = found$cells, loo = found$loo)))
  expect_identical(conditionMessage(w), sprintf(paste0(
    "summarize_by_cell(): these are not the cells and rows the count was scored on (the ",
    "number of rows differs). The profile scored 600 rows of 'z' in 133 cells with a ",
    "leave-one-out error of %s; here 570 rows in %d cells give %s. Rows were dropped or moved ",
    "between cells: points outside every cell, a boundary that clips the cells, a subset of ",
    "the rows, or another layer."),
    format(rec$scored$loo, digits = 6), found$cells, format(found$loo, digits = 6)))
  expect_null(conditionCall(w))
  expect_identical(.cr_check(r), list(estimate = "match", cells = "differ"))
  # Under plain means too: one more fit of the response, and the same report.
  r <- .cr_sum(A$asg[-drop, ], A$tess$cells, "mean")
  expect_identical(.cr_classes(r), c("spatialkit_estimate_mismatch", "spatialkit_cells_differ"))
  expect_identical(.cr_of(r, "spatialkit_cells_differ")[[1L]]$found, w$found)
})

test_that("a boundary that leaves points in no cell is reported, kept or dropped", {
  pts <- .cr_600(); A <- .cr_A()
  bb <- sf::st_bbox(pts)
  small <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(
    c(xmin = bb[["xmin"]] + 5000, ymin = bb[["ymin"]], xmax = bb[["xmax"]],
      ymax = bb[["ymax"]]), crs = sf::st_crs(pts))))
  tess <- suppressWarnings(build_tessellation(A$seeds, boundary = small, method = "voronoi",
                                              quiet = TRUE))
  expect_identical(attr(tess$cells, "count_rule")$scored_cells, TRUE)
  dropped <- .cr_catch(assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id",
                                                   quiet = TRUE))
  expect_identical(.cr_classes(dropped), "spatialkit_unassigned")
  expect_lt(nrow(dropped$value), 600L)
  kept <- assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id",
                                      quiet = TRUE, keep_unassigned = TRUE)
  expect_identical(nrow(kept), 600L)
  for (asg in list(dropped$value, kept)) {
    r <- .cr_sum(asg, tess$cells, "pooled")
    w <- .cr_of(r, "spatialkit_cells_differ")
    expect_length(w, 1L)
    expect_identical(w[[1L]]$what, "n")
    expect_identical(w[[1L]]$found$n, nrow(dropped$value))
    expect_identical(.cr_check(r)$cells, "differ")
  }
})

test_that("another layer assigned to the cells is reported", {
  A <- .cr_A()
  other <- .cr_layer(600, 5102)
  asg <- assign_features_to_polygons(other, A$tess$cells, polygon_id_col = "cell_id",
                                     quiet = TRUE)
  expect_identical(nrow(asg), 600L)
  # The assigned layer takes the record of the cells it was assigned to.
  expect_identical(attr(asg, "count_rule"), attr(A$tess$cells, "count_rule"))
  r <- .cr_sum(asg, A$tess$cells, "pooled")
  w <- .cr_of(r, "spatialkit_cells_differ")
  expect_length(w, 1L)
  # As many rows; they fill other cells, or the same cells with other values.
  expect_identical(w[[1L]]$found$n, 600L)
  expect_true(w[[1L]]$what %in% c("cells", "loo"))
  expect_identical(.cr_check(r)$cells, "differ")
})

test_that("two rows exchanged between neighbouring cells are told from rows dropped", {
  A <- .cr_A()
  asg <- A$asg
  cells <- A$tess$cells
  nb <- sf::st_touches(cells)
  i <- which(asg$cell_id == cells$cell_id[1L])[1L]
  j <- which(asg$cell_id == cells$cell_id[nb[[1L]][1L]])[1L]
  swap <- asg$cell_id[c(j, i)]
  asg$cell_id[c(i, j)] <- swap
  expect_identical(sort(table(asg$cell_id)), sort(table(A$asg$cell_id)))    # the counts stand
  r <- .cr_sum(asg, cells, "pooled")
  w <- .cr_of(r, "spatialkit_cells_differ")
  expect_length(w, 1L)
  w <- w[[1L]]
  expect_identical(w$what, "loo")
  expect_identical(w$found[c("n", "cells")], w$scored[c("n", "cells")])
  expect_false(isTRUE(all.equal(w$found$loo, w$scored$loo, tolerance = 1e-8)))
  text <- conditionMessage(w)
  expect_match(text, "(the leave-one-out error differs)", fixed = TRUE)
  expect_match(text, paste0(
    "The same number of rows and cells was found, so rows sit in other cells than the ",
    "profile put them in, or values changed: rows exactly on a boundary between two cells ",
    "(gridded locations) can be assigned to the other cell, and coordinates that reached ",
    "this CRS by another transformation can settle a tie between neighbouring cells the ",
    "other way."), fixed = TRUE)
  expect_false(grepl("Rows were dropped or moved between cells", text, fixed = TRUE))
  expect_identical(.cr_check(r)$cells, "differ")
  # Values changed on the same rows and cells are the same finding.
  asg <- A$asg; asg$z[1L] <- asg$z[1L] + 1
  expect_identical(.cr_of(.cr_sum(asg, cells, "pooled"), "spatialkit_cells_differ")[[1L]]$what,
                   "loo")
})

test_that("cells of another selection than the assigned layer's are named by the field that differs", {
  pts <- .cr_600(); A <- .cr_A()
  other <- select_resolution(.cr_profile(pts, obs_per_cell = 9), criterion = "size")
  C <- .cr_cells(pts, other)
  r <- .cr_sum(A$asg, C$tess$cells, "pooled")
  w <- .cr_of(r, "spatialkit_cells_differ")
  expect_length(w, 1L)
  expect_identical(conditionMessage(w[[1L]]), paste0(
    "summarize_by_cell(): `assigned_points_sf` and `cells_sf` carry different count records ",
    "(cells: 133 against 67), so the rows were assigned to other cells than the `cells_sf` ",
    "given. The record of the assigned layer is used."))
  expect_identical(.cr_fields(w[[1L]]),
                   list(what = "record", field = "cells", scored = "133", found = "67"))
  # The record attached is the assigned layer's, and the cells are not its.
  got <- attr(r$value, "count_rule")
  expect_identical(got$cells, 133L)
  expect_identical(got$check, list(estimate = "match", cells = "differ"))
  # The same cells and the same record on both: nothing to say.
  expect_length(.cr_of(.cr_sum(A$asg, A$tess$cells, "pooled"), "spatialkit_cells_differ"), 0L)
  # Cells of a selection made for the other estimate differ in the rule first.
  B <- .cr_B()
  w <- .cr_of(.cr_sum(A$asg, B$tess$cells, "pooled"), "spatialkit_cells_differ")[[1L]]
  expect_identical(.cr_fields(w),
                   list(what = "record", field = "rule", scored = "size", found = "loo"))
  # A hex grid sized by the same selection is other cells of that selection:
  # the two records agree in everything but whose cells they sit beside.
  hex <- suppressWarnings(build_tessellation(pts, boundary = .cr_bnd(pts), method = "hex",
                                             approx_n_cells = .cr_size(), quiet = TRUE))
  expect_true("cell_id" %in% names(hex$cells))
  expect_false(nrow(hex$cells) == 133L)
  rec_hex <- attr(hex$cells, "count_rule")
  rec_asg <- attr(A$asg, "count_rule")
  expect_identical(rec_hex$scored_cells, FALSE)
  rec_hex$scored_cells <- TRUE
  expect_identical(rec_hex, rec_asg)
  r <- .cr_sum(A$asg, hex$cells, "pooled")
  expect_identical(.cr_classes(r), "spatialkit_cells_differ")
  expect_identical(conditionMessage(r$warnings[[1L]]), paste0(
    "summarize_by_cell(): `assigned_points_sf` and `cells_sf` carry different count records ",
    "(scored_cells: TRUE against FALSE), so the rows were assigned to other cells than the ",
    "`cells_sf` given. The record of the assigned layer is used."))
  expect_identical(.cr_fields(r$warnings[[1L]]),
                   list(what = "record", field = "scored_cells", scored = "TRUE",
                        found = "FALSE"))
  expect_identical(.cr_check(r), list(estimate = "match", cells = "differ"))
  # The other way round: rows assigned to the lattice, the scored cells given.
  on_hex <- assign_features_to_polygons(pts, hex$cells, polygon_id_col = "cell_id", quiet = TRUE)
  r <- .cr_sum(on_hex, A$tess$cells, "pooled")
  expect_identical(.cr_fields(.cr_of(r, "spatialkit_cells_differ")[[1L]]),
                   list(what = "record", field = "scored_cells", scored = "FALSE",
                        found = "TRUE"))
  expect_identical(.cr_check(r)$cells, "differ")
})

test_that("the warnings of the pooled columns come first, then the record's: the two records, the estimate, the cells", {
  A <- .cr_A(); B <- .cr_B()
  # Rows of the cells chosen for plain means, thirty of them removed, with a
  # constant column beside the response; summarised as pooled values, with
  # the cells of the other selection as `cells_sf`.
  asg <- B$asg[-(1:30), ]
  asg$k <- 1
  r <- .cr_catch(summarize_by_cell(asg, "z", predictor_vars = "k", id_col = "cell_id",
                                   cells_sf = A$tess$cells, estimate = "pooled"))
  expect_identical(.cr_classes(r),
                   c("spatialkit_pool_degenerate", "spatialkit_cells_differ",
                     "spatialkit_estimate_mismatch", "spatialkit_cells_differ"))
  expect_identical(vapply(.cr_of(r, "spatialkit_cells_differ"), function(w) w$what, ""),
                   c("record", "n"))
  expect_identical(.cr_check(r), list(estimate = "mismatch", cells = "differ"))
})

test_that("the check reads the cells from the ID column the summary found", {
  B <- .cr_B()
  # `id_col` names a column the layer lacks, and the summary falls back to
  # `cell_id`; the check's own fit has to read that column too.
  r <- .cr_catch(summarize_by_cell(B$asg, "z", id_col = "no_such_column", estimate = "mean"))
  expect_identical(.cr_classes(r), character(0))
  expect_true("cell_id" %in% names(r$value))
  expect_identical(.cr_check(r), list(estimate = "match", cells = "match"))
})

test_that("a constant response has no leave-one-out error at either end, and matches", {
  pts <- .cr_600()
  pts$z <- 3
  prof <- suppressWarnings(.cr_profile(pts))
  sel <- select_resolution(prof, criterion = "size")
  expect_identical(sel$count_rule$scored$loo, NA_real_)
  C <- suppressWarnings(.cr_cells(pts, sel))
  r <- .cr_sum(C$asg, C$tess$cells, "pooled")
  expect_identical(attr(r$value, "pooled")$columns$resp_z$loo, NA_real_)
  expect_identical(.cr_check(r), list(estimate = "match", cells = "match"))
  expect_length(.cr_of(r, "spatialkit_cells_differ"), 0L)
  # With rows removed the counts differ, whatever the errors are.
  r <- .cr_sum(C$asg[-(1:30), ], C$tess$cells, "pooled")
  expect_identical(.cr_of(r, "spatialkit_cells_differ")[[1L]]$what, "n")
})


# ---- Where the check is skipped -------------------------------------------------

test_that("the check of the cells is skipped, with its reason, where the record cannot speak for them", {
  pts <- .cr_600(); sel <- .cr_size(); A <- .cr_A(); bnd <- .cr_bnd(pts)
  skipped <- function(r, reason, estimate = "match") {
    expect_identical(.cr_check(r), list(estimate = estimate,
                                        cells = paste("not checked:", reason)))
    expect_length(.cr_of(r, "spatialkit_cells_differ"), 0L)
  }

  # Another response column.
  two <- A$asg; two$z2 <- rev(two$z)
  expect_false(is.null(attr(two, "count_rule")))
  skipped(.cr_sum(two, A$tess$cells, "pooled", response = "z2"), "the response is another column")
  # A response that was not summarised: not numeric, not named, not there.
  text <- A$asg; text$z <- as.character(text$z)
  skipped(.cr_sum(text, A$tess$cells, "pooled"), "the response was not summarised")
  r <- .cr_catch(summarize_by_cell(A$asg, NULL, predictor_vars = "z", id_col = "cell_id",
                                   cells_sf = A$tess$cells, estimate = "pooled"))
  skipped(r, "the response was not summarised")
  gone <- A$asg; gone$z <- NULL
  skipped(.cr_sum(gone, A$tess$cells, "pooled"), "the response was not summarised")

  # The profile read the selection half.
  half <- suppressWarnings(.cr_profile(pts, select_on = "split"))
  sel_half <- select_resolution(half, criterion = "size")
  expect_identical(sel_half$count_rule$scored$rows, "selection half")
  H <- suppressWarnings(.cr_cells(pts, sel_half))
  expect_identical(attr(H$asg, "count_rule")$scored_cells, TRUE)
  skipped(.cr_sum(H$asg, H$tess$cells, "pooled"), "the profile read the selection half")

  # Other centres at the count: a lattice, a partition of the seeding
  # function's own, random seeds.
  other <- "other cells were built from this count"
  hex <- suppressWarnings(build_tessellation(pts, boundary = bnd, method = "hex",
                                             approx_n_cells = sel, quiet = TRUE))
  asg <- assign_features_to_polygons(pts, hex$cells, polygon_id_col = "cell_id", quiet = TRUE)
  expect_identical(attr(asg, "count_rule")$scored_cells, FALSE)
  skipped(.cr_sum(asg, hex$cells, "pooled"), other)
  bare_sel <- sel; bare_sel$seeds <- NULL
  own <- .cr_cells(pts, bare_sel, set_seed = 1, partition = "kmeans", kmeans_nstart = 10)
  skipped(.cr_sum(own$asg, own$tess$cells, "pooled"), other)
  rnd <- build_tessellation(get_voronoi_seeds(bnd, method = "random", n = sel, set_seed = 3),
                            method = "voronoi", quiet = TRUE)
  asg <- suppressWarnings(assign_features_to_polygons(pts, rnd$cells, polygon_id_col = "cell_id",
                                                      quiet = TRUE))
  skipped(.cr_sum(asg, rnd$cells, "pooled"), other)

  # A record that stands for the count alone, attached by hand.
  by_hand <- own$tess$cells
  attr(by_hand, "count_rule") <- sel$count_rule
  asg <- assign_features_to_polygons(pts, by_hand, polygon_id_col = "cell_id", quiet = TRUE)
  expect_identical(attr(asg, "count_rule")$scored_cells, NA)
  skipped(.cr_sum(asg, by_hand, "pooled"), "the record does not say whose centres these are")

  # An input with no coordinates.  The pooled columns say so themselves.
  frame <- sf::st_drop_geometry(A$asg)
  expect_false(is.null(attr(frame, "count_rule")))
  r <- .cr_catch(summarize_by_cell(frame, "z", id_col = "cell_id", estimate = "pooled"))
  expect_identical(.cr_classes(r), "spatialkit_pool_unavailable")
  skipped(r, "the input has no coordinates")
  r <- .cr_catch(summarize_by_cell(frame, "z", id_col = "cell_id", estimate = "mean"))
  expect_identical(.cr_classes(r), "spatialkit_estimate_mismatch")
  skipped(r, "the input has no coordinates", estimate = "mismatch")

  # A selection from a profile without an estimate, by a size: no scored level.
  plain <- suppressWarnings(
    .cr_profile(pts, n_levels = 6L, partition = "kmeans", nstart = 5, sac = 3e4,
                estimate = NULL, obs_per_cell = NULL))
  sel_plain <- suppressWarnings(select_resolution(plain, criterion = "size"))
  expect_null(sel_plain$count_rule$scored)
  S <- suppressWarnings(.cr_cells(pts, sel_plain))
  expect_identical(attr(S$asg, "count_rule")$scored_cells, TRUE)
  r <- .cr_sum(S$asg, S$tess$cells, "mean")
  expect_identical(.cr_check(r), list(estimate = NA_character_,
                                      cells = "not checked: the profile scored no estimate"))
  expect_identical(.cr_classes(r), character(0))
})

test_that("a check the limit on an exact search stops is recorded, and nothing is raised for it", {
  A <- .cr_A(); B <- .cr_B()
  testthat::local_mocked_bindings(.has_backend = function(pkg) FALSE)
  withr::local_options(spatialkit.exact_pairs_max = 100)
  # Plain means on cells chosen for plain means: the check's own fit.
  r <- .cr_sum(B$asg, B$tess$cells, "mean")
  expect_identical(.cr_classes(r), character(0))
  expect_identical(.cr_check(r), list(estimate = "match", cells = "not checked: work limit"))
  # Pooled values on cells sized for them: the fit of the pooled columns,
  # which says for itself that they are NA.
  r <- .cr_sum(A$asg, A$tess$cells, "pooled")
  expect_identical(.cr_classes(r), "spatialkit_pool_unavailable")
  expect_identical(.cr_check(r), list(estimate = "match", cells = "not checked: work limit"))
  # Under a limit the search passes, the same calls are checked without FNN.
  withr::local_options(spatialkit.exact_pairs_max = 1e6)
  expect_identical(.cr_check(.cr_sum(B$asg, B$tess$cells, "mean")),
                   list(estimate = "match", cells = "match"))
  expect_identical(.cr_check(.cr_sum(A$asg, A$tess$cells, "pooled")),
                   list(estimate = "match", cells = "match"))
})


# ---- A record that was lost -----------------------------------------------------

test_that("a record lost on the assigned layer is read from the cells, and one lost on both raises nothing", {
  skip_if_not_installed("dplyr")
  A <- .cr_A()
  lost <- list(mutate  = dplyr::mutate(A$asg, w = 1),
               arrange = dplyr::arrange(A$asg, z),
               select  = dplyr::select(A$asg, z, cell_id))
  for (nm in names(lost)) {
    expect_null(attr(lost[[nm]], "count_rule"), info = nm)
    r <- .cr_sum(lost[[nm]], A$tess$cells, "pooled")
    expect_identical(.cr_check(r), list(estimate = "match", cells = "match"), info = nm)
    expect_identical(.cr_classes(r), character(0), info = nm)
    # Read from the cells, a subset of the rows is still caught.
    r <- .cr_sum(lost[[nm]][-(1:30), ], A$tess$cells, "pooled")
    expect_identical(.cr_of(r, "spatialkit_cells_differ")[[1L]]$what, "n", info = nm)
  }
  # A record read from the cells names the working CRS as one read from the
  # rows does: rows that lost theirs and sit in another projected CRS are
  # measured in the record's.
  moved <- dplyr::select(sf::st_transform(A$asg, 4087), z, cell_id)
  expect_null(attr(moved, "count_rule"))
  r <- .cr_sum(moved, A$tess$cells, "pooled")
  expect_identical(attr(r$value, "pooled")$crs, attr(A$tess$cells, "count_rule")$crs)
  expect_identical(.cr_check(r), list(estimate = "match", cells = "match"))
  expect_identical(.cr_classes(r), character(0))
  # A row subset and a filter keep the record, which is then stale: the
  # check is what says so.
  expect_false(is.null(attr(A$asg[1:300, ], "count_rule")))
  expect_false(is.null(attr(dplyr::filter(A$asg, z > 0), "count_rule")))
  expect_identical(.cr_check(.cr_sum(A$asg[1:300, ], NULL, "pooled"))$cells, "differ")

  bare <- A$tess$cells
  attr(bare, "count_rule") <- NULL
  for (est in c("mean", "pooled")) {
    r <- .cr_sum(lost$select, bare, est)
    expect_identical(.cr_classes(r), character(0))
    expect_null(attr(r$value, "count_rule"))
  }
  # Without `cells_sf` there is no second place to look.
  r <- .cr_catch(summarize_by_cell(lost$select, "z", id_col = "cell_id", estimate = "mean"))
  expect_identical(.cr_classes(r), character(0))
  expect_null(attr(r$value, "count_rule"))
})


# ---- The working CRS ------------------------------------------------------------

test_that("with a record the pooled values are measured in the CRS the profile worked in", {
  A <- .cr_A(); B <- .cr_B()
  rec_crs <- attr(A$asg, "count_rule")$crs
  expect_identical(rec_crs, sf::st_crs(32617))
  base <- .cr_sum(A$asg, A$tess$cells, "pooled")$value
  for (epsg in c(26917, 3857, 4326)) {
    moved <- sf::st_transform(A$asg, epsg)
    expect_false(is.null(attr(moved, "count_rule")))
    r <- .cr_sum(moved, A$tess$cells, "pooled")
    expect_identical(attr(r$value, "pooled")$crs, rec_crs, info = epsg)
    expect_identical(.cr_check(r), list(estimate = "match", cells = "match"), info = epsg)
    expect_identical(.cr_classes(r), character(0), info = epsg)
    # The same cells are neighbours, so the values are the ones of the
    # layer in the profile's own CRS.
    expect_equal(r$value$..pooled_resp_z, base$..pooled_resp_z, tolerance = 1e-9, info = epsg)
    # The frame keeps the CRS of the cells it was joined to.
    expect_identical(sf::st_crs(r$value), sf::st_crs(A$tess$cells), info = epsg)
    # Without the record the layer is measured where it is, or where
    # ensure_projected() puts it.
    plain <- moved; attr(plain, "count_rule") <- NULL
    cells <- A$tess$cells; attr(cells, "count_rule") <- NULL
    own <- attr(suppressWarnings(.cr_sum(plain, cells, "pooled"))$value, "pooled")$crs
    if (epsg != 4326) expect_identical(own, sf::st_crs(epsg), info = epsg)
  }
  # Under plain means the check's own fit is made there too.
  r <- .cr_sum(sf::st_transform(B$asg, 3857), B$tess$cells, "mean")
  expect_identical(.cr_check(r), list(estimate = "match", cells = "match"))
  # Plain means read the coordinates for one thing: the nearest other row of
  # a row alone in its cell.  At 2 rows a cell a sixth of the rows are alone,
  # and in a CRS that stretches east against north (EPSG:4087) some of them
  # have another nearest row.  Measured there the error is another number;
  # measured in the record's CRS it is the profile's.
  fine <- select_resolution(.cr_profile(.cr_600(), estimate = "mean", obs_per_cell = 2),
                            criterion = "size")
  expect_gt(fine$count_rule$scored$share_alone, 0.1)
  C <- suppressWarnings(.cr_cells(.cr_600(), fine))
  stretched <- sf::st_transform(C$asg, 4087)
  there <- .pool_column(stretched$z, stretched$cell_id,
                        unname(sf::st_coordinates(stretched)[, 1:2]), NULL, want = "mean",
                        use_fnn = .has_backend("FNN"))$fit$loo_mean
  expect_false(.same_loo(there, fine$count_rule$scored$loo))
  r <- .cr_sum(stretched, C$tess$cells, "mean")
  expect_identical(.cr_check(r), list(estimate = "match", cells = "match"))
  expect_length(.cr_of(r, "spatialkit_cells_differ"), 0L)
  # A record that names no CRS, or a geographic one, leaves the layer to be
  # measured where it is.
  for (crs in list(sf::NA_crs_, sf::st_crs(4326))) {
    asg <- sf::st_transform(A$asg, 26917)
    rec <- attr(asg, "count_rule"); rec$crs <- crs
    attr(asg, "count_rule") <- rec
    cells <- A$tess$cells; attr(cells, "count_rule") <- NULL
    r <- suppressWarnings(.cr_sum(asg, cells, "pooled"))
    expect_identical(attr(r$value, "pooled")$crs, sf::st_crs(26917))
    expect_identical(attr(r$value, "count_rule")$crs, crs)
  }
  # .pooled_xy() puts the points in the CRS it is handed.
  at <- .pooled_xy(sf::st_transform(A$asg, 3857), rec_crs)
  expect_identical(at$crs, rec_crs)
  expect_equal(at$xy, unname(sf::st_coordinates(A$asg)), tolerance = 1e-9)
  expect_identical(.pooled_xy(sf::st_transform(A$asg, 3857))$crs, sf::st_crs(3857))
})


# ---- The pieces -----------------------------------------------------------------

test_that(".count_rule_get() reads a record and nothing else; .count_rule_put() writes one", {
  rec <- .cr_size()$count_rule
  x <- .cr_600()[1:5, ]
  expect_null(.count_rule_get(x))
  expect_identical(.count_rule_get(.count_rule_put(x, rec)), rec)
  expect_identical(.count_rule_put(x, NULL), x)
  expect_identical(.count_rule_get(.count_rule_put(x, rec, scored_cells = FALSE))$scored_cells,
                   FALSE)
  expect_identical(.count_rule_get(.count_rule_put(x, rec, scored_cells = TRUE))$scored_cells,
                   TRUE)
  # The other fields are not touched, `scored` of a record without one included.
  marked <- .count_rule_get(.count_rule_put(x, rec, scored_cells = TRUE))
  marked$scored_cells <- NA
  expect_identical(marked, rec)
  none <- rec; none["scored"] <- list(NULL)
  expect_identical(names(.count_rule_get(.count_rule_put(x, none, scored_cells = TRUE))),
                   .cr_names)
  # A list that is not a record, and a record that is not a list, are not read.
  attr(x, "count_rule") <- unclass(rec)
  expect_null(.count_rule_get(x))
  attr(x, "count_rule") <- structure("size", class = "spatialkit_count_rule")
  expect_null(.count_rule_get(x))
  expect_null(.count_rule_get(NULL))
  # Nor is a record of another version, or one without a rule, a count or
  # an estimate: every reader takes those for granted.
  unreadable <- function(edit) {
    r <- edit(rec)
    attr(x, "count_rule") <- r
    .count_rule_get(x)
  }
  expect_null(unreadable(function(r) { r$version <- 2L; r }))
  expect_null(unreadable(function(r) { r$version <- NULL; r }))
  expect_null(unreadable(function(r) { r$cells <- NULL; r }))
  expect_null(unreadable(function(r) { r["cells"] <- list(NULL); r }))
  expect_null(unreadable(function(r) { r$cells <- NA_integer_; r }))
  expect_null(unreadable(function(r) { r$rule <- NULL; r }))
  expect_null(unreadable(function(r) { r$rule <- c("size", "loo"); r }))
  expect_null(unreadable(function(r) { r$estimate <- NULL; r }))
  # An estimate that is NA is a record's own: a rule that scores no cell value.
  expect_false(is.null(unreadable(function(r) { r$estimate <- NA_character_; r })))
  expect_false(is.null(unreadable(function(r) { r["scored"] <- list(NULL); r })))
  # A summary of rows that carry such a record is a summary without one.
  A <- .cr_A()
  asg <- A$asg
  cut <- attr(asg, "count_rule"); cut$cells <- NULL; cut$version <- 2L
  attr(asg, "count_rule") <- cut
  cells <- A$tess$cells; attr(cells, "count_rule") <- NULL
  r <- .cr_sum(asg, cells, "mean")
  expect_identical(.cr_classes(r), character(0))
  expect_null(attr(r$value, "count_rule"))
  # It is not handed on either.
  seeds <- A$seeds; attr(seeds, "count_rule") <- cut
  tess <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)
  expect_false("count_rule" %in% names(tess$params))
  expect_null(attr(tess$cells, "count_rule"))
})

test_that("the print method of a record is registered", {
  expect_true(is.function(utils::getS3method("print", "spatialkit_count_rule",
                                             envir = asNamespace("spatialkit"))))
  registered <- getNamespaceInfo(asNamespace("spatialkit"), "S3methods")
  expect_true(any(registered[, 1L] == "print" & registered[, 2L] == "spatialkit_count_rule"))
})

test_that(".count_rule_differs() names the first of six fields in which two records differ", {
  a <- .cr_size()$count_rule
  expect_null(.count_rule_differs(a, a))
  expect_identical(.count_rule_same_fields,
                   c("rule", "estimate", "cells", "crs", "scored$loo", "scored_cells"))
  # Fields the comparison does not read.
  b <- a; b$seed <- 9L; b$n <- 1L; b$partition$nstart <- 7L; b$features <- "other"
  expect_null(.count_rule_differs(a, b))
  # Whose cells the record sits beside: the scored ones, others, or none.
  b <- a; b$scored_cells <- TRUE
  expect_identical(.count_rule_differs(a, b),
                   list(field = "scored_cells", a = "NA", b = "TRUE"))
  f <- b; f$scored_cells <- FALSE
  expect_identical(.count_rule_differs(b, f),
                   list(field = "scored_cells", a = "TRUE", b = "FALSE"))
  expect_null(.count_rule_differs(f, f))
  # It is the last of the six: an error that differs is named before it.
  f$scored$loo <- a$scored$loo * 2
  expect_identical(.count_rule_differs(b, f)$field, "scored$loo")
  b <- a; b$rule <- "loo"
  expect_identical(.count_rule_differs(a, b), list(field = "rule", a = "size", b = "loo"))
  b <- a; b$estimate <- "mean"
  expect_identical(.count_rule_differs(a, b), list(field = "estimate", a = "pooled", b = "mean"))
  b <- a; b$estimate <- NA_character_
  expect_identical(.count_rule_differs(a, b)$b, "NA")
  b <- a; b$cells <- 67L
  expect_identical(.count_rule_differs(a, b), list(field = "cells", a = "133", b = "67"))
  b <- a; b$crs <- sf::st_crs(26917)
  expect_identical(.count_rule_differs(a, b),
                   list(field = "crs", a = "EPSG:32617", b = "EPSG:26917"))
  b <- a; b$crs <- sf::NA_crs_
  expect_identical(.count_rule_differs(a, b)$b, "NA")
  expect_null(.count_rule_differs(b, b))
  # The same CRS written another way is the same CRS.
  b <- a; b$crs <- sf::st_crs("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
  expect_identical(isTRUE(a$crs == b$crs), is.null(.count_rule_differs(a, b)))
  # The error the profile scored: to the tolerance of the check.
  b <- a; b$scored$loo <- a$scored$loo * (1 + 5e-9)
  expect_null(.count_rule_differs(a, b))
  b$scored$loo <- a$scored$loo * (1 + 5e-8)
  expect_identical(.count_rule_differs(a, b)$field, "scored$loo")
  b <- a; b["scored"] <- list(NULL)
  expect_identical(.count_rule_differs(a, b)[c("field", "b")],
                   list(field = "scored$loo", b = "none"))
  expect_null(.count_rule_differs(b, b))
  # The first that differs is the one named.
  b <- a; b$cells <- 67L; b$rule <- "loo"; b$scored$loo <- 1
  expect_identical(.count_rule_differs(a, b)$field, "rule")
})

test_that(".same_loo() takes two missing errors as one, and numbers to a relative 1e-8", {
  expect_identical(.RECORD_TOL, 1e-8)
  expect_true(.same_loo(NA_real_, NA_real_))
  expect_false(.same_loo(NA_real_, 2))
  expect_false(.same_loo(2, NA_real_))
  expect_true(.same_loo(2, 2))
  expect_true(.same_loo(0, 0))
  expect_true(.same_loo(2e8, 2e8 * (1 + 9e-9)))
  expect_false(.same_loo(2e8, 2e8 * (1 + 2e-8)))
  expect_true(.same_loo(2e-8, 2e-8 * (1 - 9e-9)))
  expect_false(.same_loo(2e-8, 2e-8 * (1 - 2e-8)))
  expect_false(.same_loo(0, 1e-300))
})

test_that(".count_rule_check() compares the rows, then the cells, then the error, and asks for a fit only when it compares", {
  rec <- attr(.cr_A()$asg, "count_rule")
  layer <- .cr_600()
  sc <- rec$scored
  asked <- character(0)
  check <- function(found = list(), ..., record = rec, estimate = "pooled", response = "z",
                    summarised = "z", lay = layer, other = NULL) {
    asked <<- character(0)
    got <- utils::modifyList(list(n = sc$n, cells = sc$cells, loo = sc$loo, status = "ok"),
                             found)
    .cr_catch(.count_rule_check(record, other = other, estimate = estimate,
                                response_var = response, summarised = summarised, layer = lay,
                                found = function(est) { asked <<- c(asked, est); got }))
  }
  cells_of <- function(r) r$value$check$cells
  what_of  <- function(r) vapply(.cr_of(r, "spatialkit_cells_differ"), function(w) w$what, "")

  r <- check()
  expect_identical(r$value$check, list(estimate = "match", cells = "match"))
  expect_identical(asked, "pooled")                 # the estimate the profile scored
  expect_identical(.cr_classes(r), character(0))
  # The record handed back is the record handed in, with `check` appended.
  got <- r$value; got$check <- NULL
  expect_identical(got, rec)

  # The first that differs is the one named.
  expect_identical(what_of(check(list(n = 599L, cells = 1L, loo = 9))), "n")
  expect_identical(what_of(check(list(cells = 132L, loo = 9))), "cells")
  expect_identical(what_of(check(list(loo = sc$loo * (1 + 2e-8)))), "loo")
  expect_identical(cells_of(check(list(loo = sc$loo * (1 + 9e-9)))), "match")
  expect_identical(what_of(check(list(loo = NA_real_))), "loo")
  expect_identical(what_of(check(list(n = 0L, cells = 0L, loo = NA_real_,
                                      status = "no located value"))), "n")
  na <- rec; na$scored$loo <- NA_real_
  expect_identical(cells_of(check(list(loo = NA_real_), record = na)), "match")
  expect_identical(what_of(check(record = na)), "loo")
  expect_match(conditionMessage(check(record = na)$warnings[[1L]]),
               "with a leave-one-out error of NA; here 600 rows in 133 cells give 2.2", fixed = TRUE)

  # The estimate asked for is the one scored, whatever the summary reports.
  mean_rec <- rec; mean_rec$scored$estimate <- "mean"
  check(record = mean_rec)
  expect_identical(asked, "mean")
  r <- check(estimate = "mean")
  expect_identical(asked, "pooled")
  expect_identical(r$value$check, list(estimate = "mismatch", cells = "match"))

  # No fit is asked for where the check is skipped; the reasons in their order.
  no_fit <- function(r, reason) {
    expect_identical(cells_of(r), paste("not checked:", reason))
    expect_identical(asked, character(0))
  }
  none <- rec; none["scored"] <- list(NULL)
  no_fit(check(record = none), "the profile scored no estimate")
  odd <- rec; odd$scored$estimate <- "median"
  no_fit(check(record = odd), "the profile scored no estimate")
  half <- rec; half$scored$rows <- "selection half"; half$scored_cells <- FALSE
  no_fit(check(record = half), "the profile read the selection half")
  placed <- rec; placed$scored_cells <- FALSE
  no_fit(check(record = placed, response = "other"), "other cells were built from this count")
  alone <- rec; alone$scored_cells <- NA
  no_fit(check(record = alone, response = NULL), "the record does not say whose centres these are")
  no_fit(check(response = NULL, lay = sf::st_drop_geometry(layer)),
         "the response was not summarised")
  no_fit(check(response = "other", summarised = "other"), "the response is another column")
  no_fit(check(summarised = character(0), lay = sf::st_drop_geometry(layer)),
         "the response was not summarised")
  no_fit(check(lay = sf::st_drop_geometry(layer)), "the input has no coordinates")
  squares <- sf::st_buffer(layer[1:5, ], 10)
  no_fit(check(lay = squares), "the features are not points")
  # One feature that is not a point is enough.
  mixed <- rbind(layer[1:5, ], squares[1L, ])
  expect_identical(as.character(sf::st_geometry_type(mixed)),
                   c(rep("POINT", 5L), "POLYGON"))
  no_fit(check(lay = mixed), "the features are not points")
  # The record's own `features` does not stop the check of a layer of points.
  from_polygons <- rec; from_polygons$features <- "other"
  expect_identical(cells_of(check(record = from_polygons)), "match")
  # A fit the limit stopped.
  r <- check(list(n = 600L, cells = 133L, loo = NA_real_, status = "work limit"))
  expect_identical(cells_of(r), "not checked: work limit")
  expect_identical(asked, "pooled")
  expect_identical(.cr_classes(r), character(0))

  # Two records: the warning names the field, and the cells are not the ones
  # scored whatever the rows say.
  two <- rec; two$cells <- 67L
  r <- check(other = two)
  expect_identical(what_of(r), "record")
  expect_identical(r$value$check, list(estimate = "match", cells = "differ"))
  r <- check(list(n = 1L), other = two)
  expect_identical(what_of(r), c("record", "n"))
  r <- check(other = two, record = placed)
  expect_identical(what_of(r), "record")
  expect_identical(cells_of(r), "differ")
  expect_identical(what_of(check(other = rec)), character(0))
})


# ---- The print method -----------------------------------------------------------

test_that("a record prints the count, the estimate, the cells and what the profile scored", {
  sel <- .cr_size()
  rec <- sel$count_rule
  sc <- rec$scored
  scored_line <- sprintf(
    "  scored      : loo %s on 600 rows in 133 cells (one cell: %s); %.1f%% of rows alone",
    sprintf("%.2e", sc$loo), sprintf("%.2e", sc$loo_one_cell), 100 * sc$share_alone)
  # The errors are written as the record of the pooled columns writes them.
  expect_match(scored_line, "loo [0-9]\\.[0-9]{2}e[+-][0-9]{2} on 600 rows in 133 cells \\(one cell: [0-9]\\.[0-9]{2}e[+-][0-9]{2}\\)")
  head <- c(
    "Count: 133 cells, set by size (4.5 rows per cell asked; 600 rows / 133 cells = 4.51)",
    "  chosen for  : the pooled value",
    "  cells       : density, 3 density-weighted k-means++ starts (power 1, 8",
    "                neighbours; seed 123), fitted to every row;")
  expect_identical(.cr_print(rec), c(
    head, "                the record stands for the count and sits beside no centres",
    scored_line))
  expect_identical(.cr_print(attr(sel$seeds, "count_rule")), c(
    head, "                these are the cells the profile scored", scored_line))
  A <- .cr_A()
  out <- .cr_sum(A$asg, A$tess$cells, "pooled")$value
  expect_identical(.cr_print(attr(out, "count_rule")), c(
    head, "                these are the cells the profile scored", scored_line,
    "  checked     : the cells and rows summarised are the ones scored"))
  expect_true(any(grepl(sprintf("loo         : %s ", sprintf("%.2e", sc$loo)),
                        .cr_print(attr(out, "pooled")), fixed = TRUE)))
  # print() returns its argument, unseen.
  utils::capture.output(back <- withVisible(print(rec)))
  expect_identical(back, list(value = rec, visible = FALSE))
  # Large errors and a large seed are written with their digits, both errors
  # of a line one way.
  big <- rec; big$seed <- 100000; big$scored$loo <- 71800000; big$scored$loo_one_cell <- 1.76e8
  expect_identical(.cr_print(big)[c(4L, 6L)], c(
    "                neighbours; seed 100000), fitted to every row;",
    sprintf("  scored      : loo 7.18e+07 on 600 rows in 133 cells (one cell: 1.76e+08); %.1f%% of rows alone",
            100 * sc$share_alone)))
  # One row, one cell.
  single <- rec; single$cells <- 1L; single$n <- 1L; single$scored$n <- 1L
  single$scored$cells <- 1L
  out1 <- .cr_print(single)
  expect_identical(out1[1L],
    "Count: 1 cell, set by size (4.5 rows per cell asked; 1 row / 1 cell = 1.00)")
  expect_match(out1[6L], " on 1 row in 1 cell (one cell: ", fixed = TRUE)

  # A mismatch, rows that differ, a check that was skipped.
  out <- suppressWarnings(.cr_sum(A$asg[-(1:30), ], A$tess$cells, "mean"))$value
  expect_identical(utils::tail(.cr_print(attr(out, "count_rule")), 2L), c(
    "  checked     : the cells and rows summarised are not the ones scored",
    "                the summary reports the plain mean, not the pooled value"))
  placed <- rec; placed$scored_cells <- FALSE
  placed$check <- list(estimate = "match", cells = "not checked: other cells were built from this count")
  expect_identical(.cr_print(placed)[c(5L, 7L)], c(
    "                other cells were built from this count: not the cells the profile scored",
    "  not checked : other cells were built from this count"))
})

test_that("the print names the rule, the kind, the seed and the rows the centres were fitted to", {
  B <- .cr_B()
  rec <- attr(B$asg, "count_rule")
  out <- .cr_print(rec)
  expect_identical(out[1:2], c("Count: 50 cells, set by loo (600 rows / 50 cells = 12.00)",
                               "  chosen for  : the plain mean"))
  # k-means from k-means++ starts, as a profile without a partition record is.
  km <- suppressWarnings(.cr_profile(.cr_600(), partition = "kmeans", nstart = 25))
  rec_km <- select_resolution(km, criterion = "size")$count_rule
  expect_identical(.cr_print(rec_km)[3:4], c(
    "  cells       : k-means, 25 k-means++ starts (seed 123), fitted to every row;",
    "                the record stands for the count and sits beside no centres"))
  # No seed kept, one start, random starts.
  r <- rec_km; r$seed <- NA_integer_; r$partition$nstart <- 1L; r$partition$init <- "random"
  expect_identical(.cr_print(r)[3L],
                   "  cells       : k-means, 1 random start, fitted to every row;")
  r <- rec_km; r$partition$nstart <- NA_integer_
  expect_identical(.cr_print(r)[3L],
                   "  cells       : k-means, k-means++ starts (seed 123), fitted to every row;")
  # A subsample, and row counts that are not known.
  r <- rec_km; r$n_fit <- 1500L; r$n_points <- 8251L
  expect_identical(.cr_print(r)[3:4], c(
    "  cells       : k-means, 25 k-means++ starts (seed 123), fitted to a subsample",
    "                of 1500 of the 8251 points;"))
  r <- rec_km; r$n_fit <- NA_integer_
  expect_identical(.cr_print(r)[3L], "  cells       : k-means, 25 k-means++ starts (seed 123);")
  r <- rec_km; r$n <- NA_integer_; r$obs_per_cell <- NA_real_
  expect_identical(.cr_print(r)[1L], "Count: 133 cells, set by size (4.5 rows per cell asked)")
  r$rule <- "loo"
  expect_identical(.cr_print(r)[1L], "Count: 133 cells, set by loo")

  # The selection half.
  half <- suppressWarnings(.cr_profile(.cr_600(), select_on = "split"))
  s <- select_resolution(half, criterion = "size")
  on <- attr(s$seeds, "count_rule")
  out <- .cr_print(on)
  expect_identical(out[5L], "                these are the cells the profile scored, on half the rows")
  expect_identical(out[6L], sprintf(
    "  scored      : loo %s on the %d rows of the selection half in %d cells (one cell: %s); %.1f%% of rows alone",
    sprintf("%.2e", on$scored$loo), on$scored$n, on$scored$cells,
    sprintf("%.2e", on$scored$loo_one_cell), 100 * on$scored$share_alone))
  expect_lt(on$scored$n, 600L)

  # A count no cell value was scored for.
  plain <- suppressWarnings(
    .cr_profile(.cr_600(), n_levels = 6L, partition = "kmeans", nstart = 5, sac = 3e4,
                estimate = NULL, obs_per_cell = NULL))
  r <- suppressWarnings(select_resolution(plain, criterion = "size"))$count_rule
  out <- .cr_print(r)
  expect_identical(out[2L], "  chosen for  : no cell value: the profile names no estimate")
  expect_length(out, 4L)                              # no `scored` line
  r$rule <- "elbow"
  expect_identical(.cr_print(r)[2L], "  chosen for  : no cell value: elbow scores none")
  # "loo" scores the profile's estimate as a size does; on a column subset,
  # which keeps no estimate, the record names none.
  subset_rec <- select_resolution(.cr_mean7()[, c("levels", "loo")], criterion = "loo")$count_rule
  expect_identical(subset_rec$rule, "loo")
  expect_identical(subset_rec$estimate, NA_character_)
  expect_identical(.cr_print(subset_rec)[1:2], c(
    "Count: 50 cells, set by loo",
    "  chosen for  : no cell value: the profile names no estimate"))
  for (rule in c("reliability", "moran_z")) {
    r$rule <- rule
    expect_identical(.cr_print(r)[2L],
                     sprintf("  chosen for  : no cell value: %s scores none", rule))
  }
  # A constant response: errors that are not numbers are written as NA.
  r <- rec_km; r$scored$loo <- NA_real_; r$scored$loo_one_cell <- NA_real_
  r$scored$share_alone <- NA_real_
  expect_identical(.cr_print(r)[5L],
                   "  scored      : loo NA on 600 rows in 133 cells (one cell: NA)")
})
