# tests/testthat/test-fix-streams-per-step.R
# ---------------------------------------------------------------------------
# Two loops drew each step from one random stream whose position depended on
# how many draws came before it.
#
# determine_optimal_levels(): under the model-aware criteria the integer that
# seeds the re-runs in the window around the elbow was drawn after the sweep,
# which takes more draws the more counts it runs, so max_levels 30 and 31 gave
# one knee and one window and other cells, Moran's I and ranking at the same
# counts.  The integer is now read before the sweep, and the stream is put
# back after the read, so the sweep starts where it did.
#
# resolution_profile(): every level's partition drew from the stream in
# ladder order, so a level's cells, and every column read off them, depended
# on the other levels of the ladder, and dropping low levels that did not win
# moved the count select_resolution() chose.  Each level now runs under a
# stream of its own, one integer drawn after the subsample plus the level,
# and the stream is put back after it.
# ---------------------------------------------------------------------------

# A smooth field: Gaussian covariance, range 250 on a 1000 x 1000 square.
.ss_smooth_z <- function(xy, seed, a = 250, nugget = 0.05) {
  set.seed(seed)
  D <- as.matrix(stats::dist(xy))
  as.numeric(t(chol(exp(-(D / a)^2) + diag(1e-6, nrow(xy)))) %*% rnorm(nrow(xy))) +
    rnorm(nrow(xy), sd = sqrt(nugget))
}
# k separated clusters of `per` points: an elbow at k.  w is noise and z the
# smooth field.
.ss_clustered <- function(k, per = 40, seed = 1, sd = 12) {
  set.seed(seed)
  ctr <- cbind(runif(k, 0, 1000), runif(k, 0, 1000))
  g <- rep(seq_len(k), each = per)
  d <- data.frame(x = ctr[g, 1] + rnorm(k * per, 0, sd), y = ctr[g, 2] + rnorm(k * per, 0, sd))
  d$w <- rnorm(nrow(d))
  d$z <- .ss_smooth_z(d[, c("x", "y")], seed + 100)
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
}
.ss_q <- function(expr) suppressWarnings(suppressMessages(expr))
.ss_cache <- new.env(parent = emptyenv())
.ss_layer <- function() {
  if (is.null(.ss_cache$layer)) .ss_cache$layer <- .ss_clustered(12, seed = 1)
  .ss_cache$layer
}
.ss_levels <- function(max_levels, criterion, set_seed = 123L, top_n = 3L) {
  key <- paste("dol", max_levels, criterion, set_seed, top_n, sep = "|")
  if (is.null(.ss_cache[[key]]))
    .ss_cache[[key]] <- .ss_q(determine_optimal_levels(
      .ss_layer(), max_levels = max_levels, response_var = "z", predictor_vars = "w",
      criterion = criterion, top_n = top_n, set_seed = set_seed))
  .ss_cache[[key]]
}
.ss_profile <- function(..., seed = 123L) {
  args <- list(...)
  key <- paste("rp", paste(names(args), vapply(args, paste, "", collapse = ","),
                           collapse = ";"), seed, sep = "|")
  if (is.null(.ss_cache[[key]]))
    .ss_cache[[key]] <- .ss_q(resolution_profile(.ss_layer(), response_var = "z", ...,
                                                 seed = seed))
  .ss_cache[[key]]
}
# The columns read off a level's cells, at the levels `at`, as a plain list.
.ss_cols <- c("wss", "wss_spread", "cell_n_min", "cell_n_median", "cell_diam_median",
              "cells_resp", "cell_resp_mean", "rss", "cp", "moran_i", "moran_z", "loo",
              "loo_se")
.ss_at <- function(p, at) {
  d <- as.data.frame(p)
  cols <- intersect(.ss_cols, names(d))
  rows <- match(at, d$levels)
  c(lapply(stats::setNames(cols, cols), function(cl) unname(d[[cl]][rows])),
    list(centres = unname(attr(p, "centres")[as.character(at)])))
}

# ---- determine_optimal_levels(): the re-run integer --------------------------------

test_that("max_levels that give one knee and one window give one score at each count", {
  for (ss in c(123L, 5L)) for (crit in c("morans_i", "combined")) {
    info <- sprintf("set_seed = %d, %s", ss, crit)
    a <- .ss_levels(30, crit, set_seed = ss); b <- .ss_levels(31, crit, set_seed = ss)
    da <- attr(a, "diagnostics"); db <- attr(b, "diagnostics")
    # One knee and one window: the two calls score the same counts.
    expect_identical(da$knee_k, db$knee_k, info = info)
    expect_identical(da$eval_ks, db$eval_ks, info = info)
    expect_true(any(is.finite(da$moran_z[da$eval_ks])), info = info)
    ks <- da$eval_ks
    for (col in c("moran_i", "moran_z", "wss_eval"))
      expect_identical(da[[col]][ks], db[[col]][ks], info = paste(info, col))
    # The sweep up to 30 is the same, and so is the answer.
    expect_identical(da$wss[1:30], db$wss[1:30], info = info)
    expect_identical(as.vector(a), as.vector(b), info = info)
  }
  # max_levels 28 to 32 under the default set_seed: one window, one ranking.
  runs <- lapply(28:32, function(ml) .ss_levels(ml, "morans_i"))
  d <- lapply(runs, attr, "diagnostics")
  for (j in 2:5) {
    info <- sprintf("max_levels %d against 28", 27L + j)
    expect_identical(d[[j]]$eval_ks, d[[1]]$eval_ks, info = info)
    expect_identical(d[[j]]$moran_z[d[[1]]$eval_ks], d[[1]]$moran_z[d[[1]]$eval_ks],
                     info = info)
    expect_identical(as.vector(runs[[j]]), as.vector(runs[[1]]), info = info)
  }
})

test_that("the re-runs start from the same stream whatever max_levels, and the sweep from where the geometric call starts it", {
  pts <- .ss_layer()
  real <- spatialkit:::.partition_centres
  seen <- list()
  local_mocked_bindings(.partition_centres = function(xy, k, ...) {
    seen[[length(seen) + 1L]] <<- list(k = as.integer(k),
                                        state = get(".Random.seed", envir = .GlobalEnv))
    real(xy, k, ...)
  })
  states <- list()
  for (run in c("morans_i_30", "morans_i_31", "geometric_30")) {
    crit <- sub("_[0-9]+$", "", run); ml <- as.integer(sub("^.*_", "", run))
    seen <- list()
    r <- .ss_q(determine_optimal_levels(pts, max_levels = ml, response_var = "z",
                                        predictor_vars = "w", criterion = crit))
    k_max <- ml
    ks <- vapply(seen, `[[`, integer(1), "k")
    n_sweep <- k_max - 1L                          # the sweep fits k = 2 to k_max
    expect_identical(ks[seq_len(n_sweep)], 2:k_max, info = run)
    states[[run]] <- list(
      sweep = lapply(seen[seq_len(n_sweep)], `[[`, "state"),
      rerun = stats::setNames(lapply(seen[-seq_len(n_sweep)], `[[`, "state"),
                              ks[-seq_len(n_sweep)]))
  }
  # The re-run at each count starts from the same stream at max_levels 30
  # and 31.
  shared <- intersect(names(states$morans_i_30$rerun), names(states$morans_i_31$rerun))
  expect_gt(length(shared), 0L)
  expect_identical(states$morans_i_31$rerun[shared], states$morans_i_30$rerun[shared])
  # The read before the sweep puts the stream back: the sweep of a
  # model-aware call draws what the geometric call's sweep draws.
  expect_identical(states$morans_i_30$sweep, states$geometric_30$sweep)
  expect_length(states$geometric_30$rerun, 0L)
})

test_that("with set_seed = NULL the read leaves the caller's stream where the geometric call leaves it", {
  pts <- .ss_layer()
  old <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  call <- function(crit, ml) {
    set.seed(99)
    r <- .ss_q(determine_optimal_levels(pts, max_levels = ml, response_var = "z",
                                        predictor_vars = "w", criterion = crit,
                                        set_seed = NULL))
    list(r = r, after = get(".Random.seed", envir = .GlobalEnv))
  }
  g <- call("geometric", 30)
  for (crit in c("morans_i", "combined")) {
    a <- call(crit, 30); b <- call(crit, 31)
    expect_identical(a$after, g$after, info = crit)
    da <- attr(a$r, "diagnostics"); db <- attr(b$r, "diagnostics")
    expect_identical(da$eval_ks, db$eval_ks, info = crit)
    expect_identical(da$moran_z[da$eval_ks], db$moran_z[da$eval_ks], info = crit)
    expect_identical(as.vector(a$r), as.vector(b$r), info = crit)
  }
})

# ---- resolution_profile(): a stream per level --------------------------------------

test_that("a level scores the same on any ladder that holds it, under both kinds", {
  ladders <- list(c(10, 20, 30, 40), c(20, 30, 40), c(5, 10, 20, 30, 40))
  for (part in c("kmeans", "density")) {
    ps <- lapply(ladders, function(lv) .ss_profile(levels = lv, partition = part))
    ref <- .ss_at(ps[[1]], c(20, 30, 40))
    expect_true(all(is.finite(ref$wss)) && all(is.finite(ref$moran_z)), info = part)
    for (j in 2:3) {
      got <- .ss_at(ps[[j]], c(20, 30, 40))
      expect_identical(names(got), names(ref))
      for (nm in names(ref))
        expect_identical(got[[nm]], ref[[nm]],
                         info = sprintf("%s, ladder %s, %s", part,
                                        paste(ladders[[j]], collapse = ","), nm))
    }
    # The ladder given in descending order is the same ladder.
    rev_p <- .ss_profile(levels = rev(ladders[[1]]), partition = part)
    expect_identical(.ss_at(rev_p, c(10, 20, 30, 40)), .ss_at(ps[[1]], c(10, 20, 30, 40)),
                     info = part)
  }
})

test_that("a ladder of sizes with and without its largest size scores its shared levels the same", {
  a <- .ss_profile(obs_per_cell = c(40, 20, 10))
  b <- .ss_profile(obs_per_cell = c(20, 10))
  shared <- intersect(a$levels, b$levels)
  expect_length(shared, 2L)
  expect_identical(.ss_at(b, shared), .ss_at(a, shared))
})

test_that("levels a ladder adds or drops that do not win leave the count chosen", {
  full <- seq(6, 36, by = 2)
  for (crit in c("cp", "loo")) {
    best <- vapply(c(6, 8, 10, 12), function(lo)
      as.integer(.ss_q(select_resolution(.ss_profile(levels = full[full >= lo],
                                                     partition = "kmeans"),
                                         criterion = crit))$best), integer(1))
    expect_true(all(best >= 12L), info = crit)
    expect_identical(best, rep(best[1L], 4L), info = crit)
  }
})

test_that("each level's partition starts from the same stream on every ladder, one of its own", {
  pts <- .ss_layer()
  real <- spatialkit:::.partition_centres
  seen <- list()
  local_mocked_bindings(.partition_centres = function(xy, k, ...) {
    seen[[length(seen) + 1L]] <<- list(k = as.integer(k),
                                        state = get(".Random.seed", envir = .GlobalEnv))
    real(xy, k, ...)
  })
  starts <- list()
  for (lv in list(c(10, 20, 30, 40), c(20, 30, 40))) {
    seen <- list()
    .ss_q(resolution_profile(pts, response_var = "z", levels = lv, partition = "kmeans",
                             seed = 123))
    ks <- vapply(seen, `[[`, integer(1), "k")
    expect_identical(ks, as.integer(lv))
    starts[[paste(lv, collapse = ",")]] <- stats::setNames(lapply(seen, `[[`, "state"), ks)
  }
  expect_identical(starts[["20,30,40"]], starts[["10,20,30,40"]][c("20", "30", "40")])
  expect_false(any(duplicated(starts[["10,20,30,40"]])))
})

test_that("the level streams depend on seed, and with seed = NULL on the caller's stream alone", {
  # Seeded by the level alone, every seed would give one set of cells.
  a <- .ss_at(.ss_profile(levels = c(20, 30, 40), partition = "kmeans"), c(20, 30, 40))
  b <- .ss_at(.ss_profile(levels = c(20, 30, 40), partition = "kmeans", seed = 5L),
              c(20, 30, 40))
  expect_false(identical(a$wss, b$wss))
  # With seed = NULL the levels do not draw from the caller's stream: two
  # ladders leave it in one place and score their shared levels the same.
  pts <- .ss_layer()
  old <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  call <- function(lv, start = 99L) {
    set.seed(start)
    p <- .ss_q(resolution_profile(pts, response_var = "z", levels = lv, partition = "kmeans",
                                  seed = NULL))
    list(p = p, after = get(".Random.seed", envir = .GlobalEnv))
  }
  x <- call(c(10, 20, 30, 40)); y <- call(c(20, 30, 40))
  expect_identical(y$after, x$after)
  expect_identical(.ss_at(y$p, c(20, 30, 40)), .ss_at(x$p, c(20, 30, 40)))
  z <- call(c(20, 30, 40), start = 100L)
  expect_false(identical(.ss_at(z$p, c(20, 30, 40))$wss, .ss_at(y$p, c(20, 30, 40))$wss))
  # Ladders that end at other levels: the stream is put back after each
  # level, so the caller's stream moves by the one draw of the integer and
  # not by the last level's partition.  The 480 rows draw no subsample.
  u <- call(c(10, 20, 30)); v <- call(c(10, 20, 40))
  set.seed(99L)
  invisible(sample.int(.Machine$integer.max - nrow(pts), 1L))
  one_draw <- get(".Random.seed", envir = .GlobalEnv)
  expect_identical(u$after, one_draw)
  expect_identical(v$after, one_draw)
  expect_identical(x$after, one_draw)
})

test_that("the figures read off the level streams and the re-run streams hold", {
  # The first count under "morans_i" at top_n 1 to 7 on the twelve
  # clusters, max_levels = 30, the default set_seed.
  first <- vapply(1:7, function(tn) as.integer(.ss_q(determine_optimal_levels(
    .ss_layer(), max_levels = 30, response_var = "z", predictor_vars = "w",
    criterion = "morans_i", top_n = tn))[1L]), integer(1))
  expect_identical(first, c(16L, 16L, 16L, 16L, 17L, 17L, 19L))
  # A lattice of 40 by 30 locations at 250 m, in EPSG:26915, in UTM 15N in
  # US feet and with the false origin moved by 0.1 m, the default profile
  # under each kind.
  g <- expand.grid(i = 0:39, j = 0:29)
  x <- 4.5e5 + g$i * 250; y <- 4.55e6 + g$j * 250
  z <- withr::with_seed(1, sin((x - min(x)) / 1500) + cos((y - min(y)) / 1250) +
                          stats::rnorm(length(x), 0, 0.5),
                        .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
                        .rng_sample_kind = "Rejection")
  lat <- sf::st_as_sf(data.frame(x = x, y = y, z = z), coords = c("x", "y"), crs = 26915)
  frames <- list(lat, sf::st_transform(lat, "+proj=utm +zone=15 +datum=NAD83 +units=us-ft +no_defs"),
                 sf::st_transform(lat, paste("+proj=tmerc +lat_0=0 +lon_0=-93 +k=0.9996",
                                             "+x_0=500000.1 +y_0=0.1 +datum=NAD83 +units=m +no_defs")))
  # Each kind gives one `loo` in the three frames (the same regions); the
  # figure is pinned for the density kind only.  stats::kmeans() (Hartigan
  # and Wong) settles ties on this lattice by the last bits of its own
  # compiled arithmetic, which differ between platforms: from the same frame
  # and the same starts it gives other regions on Linux on x86-64 (`loo`
  # 0.2910) than on macOS on arm64 (0.2929), so the kind "kmeans" has no
  # figure that holds on both.
  loo <- lapply(c(density = "density", kmeans = "kmeans"), function(part)
    unique(vapply(frames, function(l)
      sprintf("%.4f", .ss_q(resolution_profile(l, "z", partition = part))$loo[1L]), "")))
  expect_length(loo$kmeans, 1L)
  expect_identical(loo$density, "0.3066")
})

test_that("?resolution_profile says the kind \"kmeans\" can give other cells on another platform", {
  root <- testthat::test_path("..", "..")
  # The page's source where it is here, else the installed page as Rd text.
  rp <- if (file.exists(file.path(root, "man", "resolution_profile.Rd"))) {
    paste(readLines(file.path(root, "man", "resolution_profile.Rd"), warn = FALSE),
          collapse = "\n")
  } else {
    db <- tools::Rd_db("spatialkit")
    if (is.null(db[["resolution_profile.Rd"]])) skip("the help is not here to read")
    paste(as.character(db[["resolution_profile.Rd"]]), collapse = "")
  }
  rp <- gsub("[[:space:]]+", " ", rp)
  expect_match(rp, paste0(
    "\\code{stats::kmeans()} also settles such ties by the last bits of its own compiled ",
    "arithmetic, which differ between platforms, so on a layer where ties decide cells the kind ",
    "\\code{\"kmeans\"} can return other cells on another platform from the same coordinates ",
    "and \\code{seed}, and with them every column read off the cells."), fixed = TRUE)
})

test_that("the help says a count or level scores the same in any window or ladder that holds it", {
  root <- testthat::test_path("..", "..")
  if (!(file.exists(file.path(root, "DESCRIPTION")) && dir.exists(file.path(root, "man"))))
    skip("the help's source is not here to read")
  read <- function(...) gsub("[[:space:]]+", " ",
                             paste(readLines(file.path(root, ...), warn = FALSE), collapse = " "))
  rd <- read("man", "determine_optimal_levels.Rd")
  expect_match(rd, paste0("The re-run at each count draws its starts from a stream of its own, ",
                          "seeded before the sweep, so for one \\code{set_seed} a count gets the ",
                          "same cells, \\code{wss_eval} and Moran's I in any neighbourhood that ",
                          "holds it, whichever \\code{top_n} and \\code{max_levels} set the ",
                          "neighbourhood"), fixed = TRUE)
  expect_false(grepl("with the other arguments fixed", rd, fixed = TRUE))
  rp <- read("man", "resolution_profile.Rd")
  expect_match(rp, paste0("Each level's restarts run on a stream of their own, seeded from ",
                          "\\code{seed} and the level, so for one \\code{seed} a level gets the ",
                          "same cells on any ladder that holds it and is fitted to the same ",
                          "rows"), fixed = TRUE)
  expect_match(rp, paste0("On a layer of more than \\code{sample_n} rows a ladder of counts is ",
                          "fitted to a subsample and a ladder of sizes to every row, so there the ",
                          "two can give a level other cells."), fixed = TRUE)
})
