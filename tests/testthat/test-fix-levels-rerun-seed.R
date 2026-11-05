# tests/testthat/test-fix-levels-rerun-seed.R
# ---------------------------------------------------------------------------
# Under the model-aware criteria determine_optimal_levels() re-runs k-means
# at each count of the window around the elbow, the elbow plus or minus
# max(4, top_n).  The re-runs drew from one random stream in order, so a
# window that started lower (top_n above 4) moved the cells of every count
# in it, and with them wss_eval, moran_i and moran_z; under "morans_i" that
# moved the first count (15 at top_n 1 to 4 and 16 at top_n 5 on twelve
# clusters, both counts in both windows).  Each re-run is now seeded with
# one integer drawn after the sweep plus its count, so a count scores the
# same whichever top_n sets the window.  A wider window still adds counts, and
# under "morans_i" one of those can rank first.
# ---------------------------------------------------------------------------

# A smooth field: Gaussian covariance, range 250 on a 1000 x 1000 square.
.rs_smooth_z <- function(xy, seed, a = 250, nugget = 0.05) {
  set.seed(seed)
  D <- as.matrix(stats::dist(xy))
  as.numeric(t(chol(exp(-(D / a)^2) + diag(1e-6, nrow(xy)))) %*% rnorm(nrow(xy))) +
    rnorm(nrow(xy), sd = sqrt(nugget))
}
# k separated clusters of `per` points: an elbow at k.  w is noise and z the
# smooth field.
.rs_clustered <- function(k, per = 40, seed = 1, sd = 12) {
  set.seed(seed)
  ctr <- cbind(runif(k, 0, 1000), runif(k, 0, 1000))
  g <- rep(seq_len(k), each = per)
  d <- data.frame(x = ctr[g, 1] + rnorm(k * per, 0, sd), y = ctr[g, 2] + rnorm(k * per, 0, sd))
  d$w <- rnorm(nrow(d))
  d$z <- .rs_smooth_z(d[, c("x", "y")], seed + 100)
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
}
# The calls, made once for the tests below.
.rs_cache <- new.env(parent = emptyenv())
.rs_call <- function(layer, criterion, top_n, set_seed = 123L, select_on = "all") {
  key <- paste(layer, criterion, top_n, set_seed, select_on, sep = "|")
  if (is.null(.rs_cache[[key]])) {
    pts <- switch(layer,
                  c12 = .rs_clustered(12, seed = 1),
                  c18 = .rs_clustered(18, seed = 2, sd = 20))
    .rs_cache[[key]] <- suppressWarnings(suppressMessages(
      determine_optimal_levels(pts, max_levels = 30, response_var = "z",
                               predictor_vars = "w", criterion = criterion,
                               top_n = top_n, set_seed = set_seed,
                               select_on = select_on)))
  }
  .rs_cache[[key]]
}
# For every pair of calls, the columns at the counts both windows hold.
.rs_expect_shared_same <- function(calls, cols, info) {
  for (a in seq_along(calls)) for (b in seq_along(calls)) if (a < b) {
    da <- attr(calls[[a]], "diagnostics"); db <- attr(calls[[b]], "diagnostics")
    shared <- intersect(da$eval_ks, db$eval_ks)
    expect_gt(length(shared), 0L)
    for (col in cols)
      expect_identical(da[[col]][shared], db[[col]][shared],
                       info = sprintf("%s, %s, calls %s and %s", info, col,
                                      names(calls)[a], names(calls)[b]))
  }
}

test_that("a count scores the same in every window that holds it, under both model-aware criteria", {
  for (ss in c(123L, 5L)) for (crit in c("morans_i", "combined")) {
    info <- sprintf("set_seed = %d, %s", ss, crit)
    calls <- stats::setNames(lapply(c(1L, 3L, 5L, 7L), function(tn)
      .rs_call("c12", crit, tn, set_seed = ss)), paste0("top_n=", c(1, 3, 5, 7)))
    ks <- lapply(calls, function(r) attr(r, "diagnostics")$eval_ks)
    # The windows: one for top_n 1 and 3, wider at 5 and wider again at 7.
    expect_identical(ks[[1]], ks[[2]], info = info)
    expect_true(all(ks[[2]] %in% ks[[3]]) && length(ks[[3]]) > length(ks[[2]]), info = info)
    expect_true(all(ks[[3]] %in% ks[[4]]) && length(ks[[4]]) > length(ks[[3]]), info = info)
    # At the counts the windows share, the re-run's cells give the same WSS
    # and the same Moran's I and z.
    .rs_expect_shared_same(calls, c("moran_i", "moran_z", if (crit == "combined") "wss_eval"),
                           info)
    # The sweep does not depend on the window either.
    for (j in 2:4) {
      expect_identical(attr(calls[[j]], "diagnostics")$wss, attr(calls[[1]], "diagnostics")$wss,
                       info = info)
      expect_identical(attr(calls[[j]], "diagnostics")$knee_k,
                       attr(calls[[1]], "diagnostics")$knee_k, info = info)
    }
  }
  # Eighteen looser clusters, under the split too: the half the response is
  # read on changes which cells score, not the cells.
  for (so in c("all", "split")) for (crit in c("morans_i", "combined")) {
    info <- sprintf("eighteen clusters, select_on = %s, %s", so, crit)
    calls <- stats::setNames(lapply(c(1L, 7L), function(tn)
      .rs_call("c18", crit, tn, select_on = so)), paste0("top_n=", c(1, 7)))
    .rs_expect_shared_same(calls, c("moran_i", "moran_z", if (crit == "combined") "wss_eval"),
                           info)
  }
})

test_that("under \"morans_i\" a wider window moves the first count only to a count it adds", {
  for (ss in c(123L, 5L)) {
    res <- lapply(1:7, function(tn) .rs_call("c12", "morans_i", tn, set_seed = ss))
    d <- lapply(res, attr, "diagnostics")
    for (tn in 1:7) {
      info <- sprintf("set_seed = %d, top_n = %d", ss, tn)
      # The counts returned are the finite counts of the window by |z|.
      fin <- d[[tn]]$eval_ks[is.finite(d[[tn]]$moran_z[d[[tn]]$eval_ks])]
      ranked <- fin[order(abs(d[[tn]]$moran_z[fin]))]
      expect_identical(as.vector(res[[tn]]), as.integer(head(ranked, tn)), info = info)
    }
    # top_n 1 to 4 score one window: the same first count, and each result
    # the head of the next.
    for (tn in 2:4) {
      info <- sprintf("set_seed = %d, top_n = %d", ss, tn)
      expect_identical(d[[tn]]$eval_ks, d[[1]]$eval_ks, info = info)
      expect_identical(as.vector(res[[tn]])[seq_len(tn - 1L)], as.vector(res[[tn - 1L]]),
                       info = info)
    }
    # A wider window puts first either the narrower window's first count or
    # a count the narrower window does not hold, with a smaller |z|.
    for (a in 1:6) for (b in (a + 1L):7) {
      info <- sprintf("set_seed = %d, top_n %d and %d", ss, a, b)
      fa <- res[[a]][1L]; fb <- res[[b]][1L]
      if (fa != fb) {
        expect_false(fb %in% d[[a]]$eval_ks, info = info)
        expect_lt(abs(d[[b]]$moran_z[fb]), abs(d[[b]]$moran_z[fa]), label = info)
      } else {
        expect_identical(d[[a]]$moran_z[fa], d[[b]]$moran_z[fb], info = info)
      }
    }
  }
})

test_that("the re-run at a count starts from the same stream whichever window holds it", {
  pts <- .rs_clustered(12, seed = 1)
  real <- spatialkit:::.partition_centres
  seen <- list()
  local_mocked_bindings(.partition_centres = function(xy, k, ...) {
    seen[[length(seen) + 1L]] <<- list(k = as.integer(k),
                                        state = get(".Random.seed", envir = .GlobalEnv))
    real(xy, k, ...)
  })
  starts <- list()
  for (tn in c(1L, 5L, 7L)) {
    seen <- list()
    r <- suppressWarnings(suppressMessages(
      determine_optimal_levels(pts, max_levels = 30, response_var = "z",
                               predictor_vars = "w", criterion = "morans_i", top_n = tn)))
    d <- attr(r, "diagnostics")
    n_sweep <- length(d$wss) - 1L                     # the sweep fits k = 2 to k_max
    ks <- vapply(seen, `[[`, integer(1), "k")
    expect_identical(ks[seq_len(n_sweep)], 2:length(d$wss))
    expect_identical(ks[-seq_len(n_sweep)], as.integer(d$eval_ks))
    starts[[as.character(tn)]] <- list(
      sweep = lapply(seen[seq_len(n_sweep)], `[[`, "state"),
      rerun = stats::setNames(lapply(seen[-seq_len(n_sweep)], `[[`, "state"), d$eval_ks))
  }
  # The sweep's stream is the same for every window.
  expect_identical(starts[["5"]]$sweep, starts[["1"]]$sweep)
  expect_identical(starts[["7"]]$sweep, starts[["1"]]$sweep)
  # The re-run at each count starts where it does in every other window.
  for (b in c("5", "7")) {
    shared <- intersect(names(starts[["1"]]$rerun), names(starts[[b]]$rerun))
    expect_gt(length(shared), 0L)
    expect_identical(starts[[b]]$rerun[shared], starts[["1"]]$rerun[shared], info = b)
  }
  shared <- intersect(names(starts[["5"]]$rerun), names(starts[["7"]]$rerun))
  expect_identical(starts[["7"]]$rerun[shared], starts[["5"]]$rerun[shared])
  # Each count has a stream of its own.
  r7 <- starts[["7"]]$rerun
  expect_false(any(duplicated(r7)))
})

test_that("the re-runs leave the caller's stream as it was and do not read it", {
  pts <- .rs_clustered(12, seed = 1)
  call <- function() suppressWarnings(suppressMessages(
    determine_optimal_levels(pts, max_levels = 30, response_var = "z", predictor_vars = "w",
                             criterion = "morans_i", top_n = 7)))
  set.seed(42)
  before <- get(".Random.seed", envir = .GlobalEnv)
  first <- call()
  expect_identical(get(".Random.seed", envir = .GlobalEnv), before)
  # The same call on another stream, or on none, gives the same result.
  set.seed(7)
  expect_identical(call(), first)
  old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  rm(".Random.seed", envir = .GlobalEnv)
  expect_identical(call(), first)
})

test_that("the re-runs depend on set_seed", {
  # Seeded by the count alone, the re-runs would give the same cells under
  # every set_seed.  At the counts both windows hold, the scores differ.
  for (crit in c("morans_i", "combined")) {
    da <- attr(.rs_call("c12", crit, 1L, set_seed = 123L), "diagnostics")
    db <- attr(.rs_call("c12", crit, 1L, set_seed = 5L), "diagnostics")
    shared <- intersect(da$eval_ks, db$eval_ks)
    expect_gt(length(shared), 0L)
    # wss_eval is stored under "combined" only.
    for (col in c("moran_z", if (crit == "combined") "wss_eval"))
      expect_false(identical(da[[col]][shared], db[[col]][shared]),
                   info = sprintf("%s, %s", crit, col))
  }
})

test_that("with set_seed = NULL a count scores the same in every window and the stream after the call does not depend on top_n", {
  pts <- .rs_clustered(12, seed = 1)
  old <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  call <- function(crit, tn, start) {
    set.seed(start)
    r <- suppressWarnings(suppressMessages(
      determine_optimal_levels(pts, max_levels = 30, response_var = "z",
                               predictor_vars = "w", criterion = crit, top_n = tn,
                               set_seed = NULL)))
    list(r = r, after = get(".Random.seed", envir = .GlobalEnv))
  }
  for (crit in c("morans_i", "combined")) {
    runs <- stats::setNames(lapply(c(1L, 5L, 7L), function(tn) call(crit, tn, 99L)),
                            paste0("top_n=", c(1, 5, 7)))
    .rs_expect_shared_same(lapply(runs, `[[`, "r"),
                           c("moran_i", "moran_z", if (crit == "combined") "wss_eval"),
                           sprintf("set_seed = NULL, %s", crit))
    for (j in 2:3)
      expect_identical(runs[[j]]$after, runs[[1]]$after,
                       info = sprintf("%s, %s", crit, names(runs)[j]))
    # The caller's stream seeds the re-runs: another start gives other scores.
    other <- attr(call(crit, 1L, 100L)$r, "diagnostics")
    d1 <- attr(runs[[1]]$r, "diagnostics")
    shared <- intersect(d1$eval_ks, other$eval_ks)
    expect_gt(length(shared), 0L)
    expect_false(identical(d1$moran_z[shared], other$moran_z[shared]), info = crit)
  }
})

test_that("?determine_optimal_levels and NEWS say a count scores the same whichever top_n sets the window", {
  root <- testthat::test_path("..", "..")
  if (!(file.exists(file.path(root, "DESCRIPTION")) && dir.exists(file.path(root, "man"))))
    skip("the help's source is not here to read")
  read <- function(...) gsub("[[:space:]]+", " ",
                             paste(readLines(file.path(root, ...), warn = FALSE), collapse = " "))
  rd <- read("man", "determine_optimal_levels.Rd")
  expect_match(rd, paste0("Under those it also sets the counts scored, the elbow plus or minus ",
                          "\\code{max(4, top_n)}: a \\code{top_n} above 4 scores more counts, and ",
                          "under \\code{\"morans_i\"} one of them can rank first. A count scores ",
                          "the same whichever \\code{top_n} puts it in the window."), fixed = TRUE)
  expect_match(rd, paste0("The re-run at each count draws its starts from a stream of its own, so ",
                          "with the other arguments fixed a count gets the same cells, ",
                          "\\code{wss_eval} and Moran's I whichever \\code{top_n} sets the ",
                          "neighbourhood"), fixed = TRUE)
  news <- read("NEWS.md")
  expect_match(news, paste0("Each re-run is now seeded with one integer drawn after the sweep plus ",
                            "its count, and the stream is put back after it, so with the other ",
                            "arguments fixed a count gets the same cells, `wss_eval` and Moran's I ",
                            "whichever `top_n` sets the window."),
               fixed = TRUE)
  expect_match(news, paste0("can fall back to the geometric ranking where it ranked, or rank where ",
                            "it fell back"), fixed = TRUE)
  expect_match(news, paste0("With `set_seed = NULL` the re-runs take one draw from the caller's ",
                            "stream in place of the draws of every re-run, so the stream after the ",
                            "call no longer depends on `top_n`."), fixed = TRUE)
})
