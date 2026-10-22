# tests/testthat/test-fix-levels-no-elbow.R
# ---------------------------------------------------------------------------
# determine_optimal_levels() returns NA, not the ladder's count, when its WSS
# curve has no elbow, with a warning of class spatialkit_no_elbow that says
# why.  on_no_elbow = "na" is the default; "warn" returns the count the chord
# rule lands on, with the same warning, as the default did before.
# ---------------------------------------------------------------------------

fne_catch <- function(expr) {
  ws <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    ws[[length(ws) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = ws, classes = unlist(lapply(ws, class)))
}
# 600 uniform points: no cluster structure, so no elbow.
fne_uniform <- function(n = 600L, seed = 1L, response = FALSE) {
  set.seed(seed)
  d <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000))
  if (response) { d$w <- rnorm(n); d$z <- d$w + rnorm(n) }
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
}
# Eight separated clusters of 50 points: an elbow at 8.
fne_clusters <- function() {
  set.seed(6)
  ctr <- cbind(runif(8, 0, 1000), runif(8, 0, 1000))
  g <- rep(1:8, each = 50)
  sf::st_as_sf(data.frame(x = ctr[g, 1] + rnorm(400, 0, 8), y = ctr[g, 2] + rnorm(400, 0, 8)),
               coords = c("x", "y"), crs = 32632)
}
# 24 separated clusters of 25 points with a response that varies by cluster:
# an elbow at 24, above the nine cells Moran's I needs, so the model-aware
# criteria rank the window around it.
fne_clusters24 <- function() {
  set.seed(12)
  ctr <- cbind(runif(24, 0, 1e4), runif(24, 0, 1e4))
  g <- rep(1:24, each = 25)
  d <- data.frame(x = ctr[g, 1] + rnorm(600, 0, 40), y = ctr[g, 2] + rnorm(600, 0, 40))
  d$w <- rnorm(600)
  d$z <- d$w + g / 10 + rnorm(600)
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
}
# Two separated clusters of 25 points: an elbow at 2.
fne_two <- function(response = FALSE) {
  set.seed(1)
  d <- data.frame(x = 5e5 + c(runif(25, 0, 10), runif(25, 90, 100)),
                  y = 5e6 + c(runif(25, 0, 10), runif(25, 90, 100)))
  if (response) d$z <- rnorm(50)
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
}
# The root of the source tree, or NULL when the tests run from an installed
# package (README.Rmd and man/ are then not there to read).
fne_source <- function() {
  root <- testthat::test_path("..", "..")
  if (file.exists(file.path(root, "DESCRIPTION")) && dir.exists(file.path(root, "man")))
    normalizePath(root)
  else NULL
}


test_that("on_no_elbow is \"na\" by default", {
  expect_identical(eval(formals(determine_optimal_levels)$on_no_elbow),
                   c("na", "warn", "error"))
  src <- fne_source()
  if (is.null(src)) skip("the help's source is not here to read")
  rd <- paste(readLines(file.path(src, "man", "determine_optimal_levels.Rd"), warn = FALSE),
              collapse = " ")
  rd <- gsub("[[:space:]]+", " ", rd)
  expect_match(rd, "on_no_elbow = c(\"na\", \"warn\", \"error\")", fixed = TRUE)
  expect_match(rd, "\\code{\"na\"} (default) returns \\code{NA_integer_}", fixed = TRUE)
  expect_false(grepl("\\code{\"warn\"} (default)", rd, fixed = TRUE))
})

test_that("a layer with no elbow returns NA at the default, with the flag, the sag and the classed warning", {
  uni <- fne_uniform()
  lines <- capture_spatialkit_log(out <- fne_catch(determine_optimal_levels(uni)))
  k <- out$value
  # No count: NA, with structured = FALSE and the largest sag, and one
  # warning of class spatialkit_no_elbow that says why and that NA is
  # returned; the log has the same line.
  expect_identical(length(k), 1L)
  expect_true(is.na(k))
  expect_type(k, "integer")
  expect_identical(attr(k, "structured"), FALSE)
  expect_lt(attr(k, "sag"), log(1.25))
  expect_named(attributes(k), c("structured", "sag"))
  expect_length(out$warnings, 1L)
  cond <- out$warnings[[1L]]
  expect_s3_class(cond, "spatialkit_no_elbow")
  expect_identical(cond$structured, FALSE)
  expect_equal(cond$sag, attr(k, "sag"))
  expect_match(conditionMessage(cond), "the WSS curve has no elbow.*set by the ladder")
  expect_match(conditionMessage(cond), "Returning NA (on_no_elbow = \"na\").", fixed = TRUE)
  expect_true(log_has(lines, "the WSS curve has no elbow"))
  expect_true(log_has(lines, "Returning NA \\(on_no_elbow = \"na\"\\)"))
  na <- fne_catch(determine_optimal_levels(uni, on_no_elbow = "na"))
  expect_identical(na$value, k)
  expect_identical(na$warnings, out$warnings)
  # A handler for the class finds it, and a caller who turns warnings into
  # errors is stopped.
  expect_identical(tryCatch(determine_optimal_levels(uni),
                            spatialkit_no_elbow = function(w) "caught"), "caught")
  op <- options(warn = 2)
  stopped <- tryCatch(determine_optimal_levels(uni), error = function(e) conditionMessage(e))
  options(op)
  expect_type(stopped, "character")
  expect_match(stopped, "the WSS curve has no elbow")

  # "warn": the count the chord rule lands on, which the ladder set, with the
  # same warning, which carries it and the same sag, and no "Returning NA".
  w <- fne_catch(determine_optimal_levels(uni, on_no_elbow = "warn"))
  expect_true(all(is.finite(w$value)))
  expect_null(attributes(w$value))
  expect_identical(w$classes[w$classes == "spatialkit_no_elbow"], "spatialkit_no_elbow")
  wc <- w$warnings[[which(vapply(w$warnings, inherits, logical(1), "spatialkit_no_elbow"))]]
  expect_identical(wc$k, w$value[1L])
  expect_identical(wc$structured, FALSE)
  expect_equal(wc$sag, attr(k, "sag"))
  expect_match(conditionMessage(wc), "set by the ladder")
  expect_false(grepl("Returning NA", conditionMessage(wc), fixed = TRUE))
  # The default's condition carries the same count as `k`, the one "warn"
  # returns, and the same reason.
  expect_identical(cond$k, wc$k)
  expect_identical(conditionMessage(cond),
                   paste(conditionMessage(wc), "Returning NA (on_no_elbow = \"na\")."))

  # "error": the same reason, as an error.
  e <- tryCatch(determine_optimal_levels(uni, on_no_elbow = "error"), error = function(e) e)
  expect_s3_class(e, "spatialkit_no_elbow")
  expect_identical(e$k, w$value[1L])
})

test_that("a clustered layer returns the same count at the default as before, with structured = TRUE", {
  cl <- fne_clusters()
  out <- fne_catch(determine_optimal_levels(cl, max_levels = 14))
  k <- out$value
  expect_length(out$warnings, 0L)
  expect_identical(k[1L], 8L)
  # The count "warn" returns, which is what the default returned before; the
  # default adds the flag and the sag.
  old <- determine_optimal_levels(cl, max_levels = 14, on_no_elbow = "warn")
  expect_null(attributes(old))
  expect_identical(as.integer(k), old)
  expect_named(attributes(k), c("structured", "sag"))
  expect_identical(attr(k, "structured"), TRUE)
  expect_gte(attr(k, "sag"), log(1.25))
  expect_identical(determine_optimal_levels(cl, max_levels = 14, on_no_elbow = "na"), k)
  expect_identical(determine_optimal_levels(cl, max_levels = 14, on_no_elbow = "error"), old)
  # top_n and the split keep the flag.
  k1 <- determine_optimal_levels(cl, max_levels = 14, top_n = 1)
  expect_identical(as.integer(k1), 8L)
  expect_identical(attr(k1, "structured"), TRUE)
  ks <- suppressWarnings(determine_optimal_levels(cl, max_levels = 14, select_on = "split"))
  expect_identical(as.integer(ks), old)
  expect_identical(attr(ks, "structured"), TRUE)
  expect_false(is.null(attr(ks, "split")))
})

test_that("with no elbow the model-aware criteria score nothing at the default, and the NA carries the sweep's diagnostics", {
  uni <- fne_uniform(400L, seed = 5L, response = TRUE)
  for (args in list(list(criterion = "morans_i", response_var = "z"),
                    list(criterion = "combined", response_var = "z", predictor_vars = "w"))) {
    call <- function(...) do.call(determine_optimal_levels,
                                  c(list(uni, max_levels = 40), args, list(...)))
    lines <- capture_spatialkit_log(out <- fne_catch(call()))
    expect_true(is.na(out$value))
    expect_identical(attr(out$value, "structured"), FALSE)
    expect_true("spatialkit_no_elbow" %in% out$classes)
    expect_false("spatialkit_diagnostic_ranking" %in% out$classes)
    # The warning and the log say that no window is scored, not that one is.
    msg <- conditionMessage(out$warnings[[1L]])
    expect_match(msg, "the model-aware criteria have no window to score", fixed = TRUE)
    expect_false(grepl("evaluated around it", msg, fixed = TRUE))
    expect_false(log_has(lines, "evaluated around it"))
    # Under "warn" the ladder places the window; Moran's z ranks it under
    # "morans_i", with the classed warning, and "combined" ranks it by the
    # sag, with none.
    w <- fne_catch(call(on_no_elbow = "warn"))
    expect_true(all(is.finite(w$value)))
    expect_true("spatialkit_no_elbow" %in% w$classes)
    expect_identical("spatialkit_diagnostic_ranking" %in% w$classes,
                     identical(args$criterion, "morans_i"))
    expect_match(conditionMessage(w$warnings[[1L]]),
                 "the model-aware criteria are evaluated around it", fixed = TRUE)
    dw <- attr(w$value, "diagnostics")
    expect_identical(dw$structured, FALSE)
    # The NA carries the diagnostics of the same sweep, with no count scored.
    d <- attr(out$value, "diagnostics")
    expect_named(d, c("moran_i", "moran_z", "wss", "wss_spread", "wss_bumps", "nstart",
                      "knee_k", "failed_k", "eval_ks", "structured", "sag"))
    same <- c("wss", "wss_spread", "wss_bumps", "nstart", "knee_k", "failed_k", "structured")
    expect_identical(d[same], dw[same])
    expect_equal(d$sag, attr(out$value, "sag"))
    expect_identical(d$knee_k, out$warnings[[1L]]$k)
    expect_identical(d$eval_ks, integer(0))
    expect_true(all(is.na(d$moran_i)) && all(is.na(d$moran_z)))
    expect_length(d$moran_z, length(dw$moran_z))
  }
})

test_that("the line before the sweep promises the fallback only where there is an elbow to fall back to", {
  # max_levels = 8 leaves no count Moran's I can score.  At the default the
  # line says the call returns NA when the curve has no elbow, and it does.
  uni <- fne_uniform(400L, seed = 5L, response = TRUE)
  line <- function(lines) grep("max_levels leaves k_max = 8", lines, value = TRUE)
  lines <- capture_spatialkit_log(out <- fne_catch(
    determine_optimal_levels(uni, max_levels = 8, response_var = "z", predictor_vars = "w",
                             criterion = "combined")))
  expect_true(is.na(out$value))
  expect_length(line(lines), 1L)
  expect_match(line(lines), paste0("criterion = 'combined' will fall back to the geometric ",
                                   "ranking if the WSS curve has an elbow (with none the call ",
                                   "returns NA). Raise max_levels"), fixed = TRUE)
  lines <- capture_spatialkit_log(e <- tryCatch(
    determine_optimal_levels(uni, max_levels = 8, response_var = "z", predictor_vars = "w",
                             criterion = "combined", on_no_elbow = "error"),
    error = function(e) e))
  expect_s3_class(e, "spatialkit_no_elbow")
  expect_match(line(lines), "if the WSS curve has an elbow (with none the call stops). Raise",
               fixed = TRUE)
  # Under "warn" the call does fall back, and the line says so as it did.
  lines <- capture_spatialkit_log(out <- fne_catch(
    determine_optimal_levels(uni, max_levels = 8, response_var = "z", predictor_vars = "w",
                             criterion = "combined", on_no_elbow = "warn")))
  expect_true(all(is.finite(out$value)))
  expect_match(line(lines), "will fall back to the geometric ranking. Raise max_levels",
               fixed = TRUE)
  # A curve with an elbow falls back at the default, as the line says.
  lines <- capture_spatialkit_log(
    k <- determine_optimal_levels(fne_two(response = TRUE), max_levels = 6, response_var = "z",
                                  criterion = "morans_i"))
  expect_match(grep("max_levels leaves k_max = 6", lines, value = TRUE),
               "if the WSS curve has an elbow (with none the call returns NA)", fixed = TRUE)
  expect_identical(as.vector(k), c(2L, 1L, 3L))
  expect_identical(attr(k, "structured"), TRUE)
})

test_that("a model-aware count at the default carries the flag and the sag beside the diagnostics", {
  c24 <- fne_clusters24()
  for (args in list(list(criterion = "morans_i", response_var = "z"),
                    list(criterion = "combined", response_var = "z", predictor_vars = "w"))) {
    call <- function(...) suppressWarnings(do.call(determine_optimal_levels,
                                                   c(list(c24, max_levels = 40), args, list(...))))
    k <- call()
    old <- call(on_no_elbow = "warn")
    # The same ranking and diagnostics as "warn" (the default before), with
    # the flag and the sag added at the default.
    expect_identical(as.vector(k), as.vector(old))
    expect_identical(attr(k, "diagnostics"), attr(old, "diagnostics"))
    expect_identical(attr(k, "diagnostics")$knee_k, 24L)
    expect_identical(attr(k, "structured"), TRUE)
    expect_gte(attr(k, "sag"), log(1.25))
    expect_equal(attr(k, "sag"), attr(k, "diagnostics")$sag)
    expect_named(attributes(k), c("diagnostics", "structured", "sag"))
    expect_named(attributes(old), "diagnostics")
  }
})

test_that("a ladder too short to show an elbow returns NA at the default", {
  two <- fne_two()
  for (ml in 1:2) {
    out <- fne_catch(determine_optimal_levels(two, max_levels = ml))
    expect_identical(as.vector(out$value), NA_integer_)
    expect_identical(attr(out$value, "structured"), FALSE)
    expect_named(attributes(out$value), c("structured", "sag"))
    expect_identical(out$classes[out$classes == "spatialkit_no_elbow"], "spatialkit_no_elbow")
    expect_match(conditionMessage(out$warnings[[1L]]),
                 "too short to read an elbow from.*Returning NA")
    # "warn" returns the ladder's count, with the same warning.
    w <- fne_catch(determine_optimal_levels(two, max_levels = ml, on_no_elbow = "warn"))
    expect_true(all(is.finite(w$value)))
    expect_true("spatialkit_no_elbow" %in% w$classes)
  }
  # Three levels show the bend at two.
  expect_identical(as.vector(determine_optimal_levels(two, max_levels = 3))[1L], 2L)
})

test_that("the default's NA stops the functions that take a count, and they say where it came from", {
  uni <- fne_uniform()
  expect_warning(k <- determine_optimal_levels(uni), class = "spatialkit_no_elbow")
  bnd <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(uni)))
  why <- paste0("It is what determine_optimal_levels() returns at its default, ",
                "on_no_elbow = \"na\", when the points name no count")
  expect_error(get_voronoi_seeds(method = "kmeans", n = k, sample_points = uni), why,
               fixed = TRUE)
  expect_error(build_tessellation(uni, boundary = bnd, method = "hex", approx_n_cells = k),
               why, fixed = TRUE)
  # The note is about that NA alone: a bare NA (k[1] drops the attributes),
  # an NA typed by hand and a count of 0 are refused without it.
  for (bad in list(k[1], NA_integer_, 0L)) {
    e <- tryCatch(get_voronoi_seeds(bnd, method = "random", n = bad), error = function(e) e)
    expect_match(conditionMessage(e), "must resolve to a positive number of cells")
    expect_false(grepl("returns at its default", conditionMessage(e), fixed = TRUE))
  }
  # A count the elbow chose is taken as before: the flag and the sag change
  # neither the cells nor the record of where the count came from.
  cl <- fne_clusters()
  kc <- determine_optimal_levels(cl, max_levels = 14)
  kw <- determine_optimal_levels(cl, max_levels = 14, on_no_elbow = "warn")
  bc <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(cl)))
  s_new <- get_voronoi_seeds(method = "kmeans", n = kc, sample_points = cl, set_seed = 3)
  s_old <- get_voronoi_seeds(method = "kmeans", n = kw, sample_points = cl, set_seed = 3)
  expect_identical(s_new, s_old)
  expect_identical(nrow(s_new), 8L)
  expect_identical(attr(s_new, "n_from"), "the first of 3 ranked candidates (8, 7, 9)")
  t_new <- build_tessellation(cl, boundary = bc, method = "hex", approx_n_cells = kc, quiet = TRUE)
  t_old <- build_tessellation(cl, boundary = bc, method = "hex", approx_n_cells = kw, quiet = TRUE)
  expect_identical(t_new, t_old)
})

test_that("the README offers determine_optimal_levels() for what it does", {
  src <- fne_source()
  if (is.null(src)) skip("the README is not here to read")
  for (file in c("README.Rmd", "README.md")) {
    lines <- readLines(file.path(src, file), warn = FALSE, encoding = "UTF-8")
    row <- grep("^\\|.*\\| *`determine_optimal_levels\\(\\)` *\\|$", lines, value = TRUE)
    expect_length(row, 1L)
    task <- trimws(strsplit(row, "|", fixed = TRUE)[[1L]][2L])
    # It reads a count off the elbow of the WSS curve and returns NA where
    # there is none; it does not choose a count for any layer, quickly or
    # not, and it does not count clusters (an elongated extent bends too).
    expect_false(grepl("choose how many cells", task, fixed = TRUE), info = file)
    expect_false(grepl("count the clusters", task, fixed = TRUE), info = file)
    expect_identical(task, paste0("read a cell count off the elbow of the clustering of my ",
                                  "points (coordinates only; `NA` when there is none)"),
                     info = file)
  }
})

test_that("the help, ?spatialkit and the resolution vignette say what the NA means", {
  src <- fne_source()
  if (is.null(src)) skip("the help's source is not here to read")
  read <- function(...) gsub("[[:space:]]+", " ",
                             paste(readLines(file.path(src, ...), warn = FALSE), collapse = " "))
  # ?spatialkit and the vignette tie the NA to the elbow, not to clusters.
  pkg <- read("man", "spatialkit-package.Rd")
  expect_match(pkg, "returns \\code{NA} when the WSS curve of their k-means clustering has no elbow",
               fixed = TRUE)
  expect_false(grepl("when the points have no clusters", pkg, fixed = TRUE))
  vig <- read("vignettes", "resolution.Rmd")
  expect_match(vig, "when the WSS curve of a k-means clustering of the points has an elbow, and `NA` when it has none",
               fixed = TRUE)
  expect_false(grepl("when the points cluster, and `NA` when they do not", vig, fixed = TRUE))
  # The vignette names `lv`, which keeps the attributes the note reads, as
  # what the two functions take.
  expect_match(vig, "`lv` is what `get_voronoi_seeds()` takes as `n`", fixed = TRUE)
  expect_false(grepl("`lv[1]` is what", vig, fixed = TRUE))
  # The help says the NA comes with the class, and no longer that it is
  # only logged or that the default raises nothing.
  rd <- read("man", "determine_optimal_levels.Rd")
  expect_match(rd, "a warning of class \\code{spatialkit_no_elbow} says why", fixed = TRUE)
  expect_false(grepl("a logged warning says why", rd, fixed = TRUE))
  rp <- read("man", "resolution_profile.Rd")
  expect_false(grepl("raises neither", rp, fixed = TRUE))
  expect_match(rp, "A warning beside the \\code{NA} the function returns at the default",
               fixed = TRUE)
})
