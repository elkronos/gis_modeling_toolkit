# tests/testthat/test-fix-levels-combined.R
# ---------------------------------------------------------------------------
# determine_optimal_levels() keeps the criterion it is given: both
# response_var and predictor_vars used to switch "geometric" to "combined",
# also when "geometric" was named.  "combined" ranks the counts around the
# elbow by the sweep's log-log sag alone, so the elbow comes first; Moran's
# z of the cell means, which put first the lowest count it can be computed
# at on a smooth field, stays in the diagnostics, and where it can be
# computed nowhere the ranking is the same.  The
# spatialkit_diagnostic_ranking warning, now raised under "morans_i" alone,
# names a way on that runs on the default profile.
# ---------------------------------------------------------------------------

flc_catch <- function(expr) {
  ws <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    ws[[length(ws) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = ws,
       classes = vapply(ws, function(w) class(w)[1L], character(1)))
}
# A smooth field: Gaussian covariance, range 250 on a 1000 x 1000 square.
flc_smooth_z <- function(xy, seed, a = 250, nugget = 0.05) {
  set.seed(seed)
  D <- as.matrix(stats::dist(xy))
  as.numeric(t(chol(exp(-(D / a)^2) + diag(1e-6, nrow(xy)))) %*% rnorm(nrow(xy))) +
    rnorm(nrow(xy), sd = sqrt(nugget))
}
# k separated clusters of `per` points: an elbow at k.  w is noise and z the
# smooth field.
flc_clustered <- function(k, per = 40, seed = 1, sd = 12) {
  set.seed(seed)
  ctr <- cbind(runif(k, 0, 1000), runif(k, 0, 1000))
  g <- rep(seq_len(k), each = per)
  d <- data.frame(x = ctr[g, 1] + rnorm(k * per, 0, sd), y = ctr[g, 2] + rnorm(k * per, 0, sd))
  d$w <- rnorm(nrow(d))
  d$z <- flc_smooth_z(d[, c("x", "y")], seed + 100)
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
}
# Uniform points with an exponential field and a nugget (the field of
# test-fix-resolution-moran-diagnostic.R): no elbow.
flc_field <- function(n = 500, seed = 2, a = 100, nugget = 0.3) {
  set.seed(seed)
  xy <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000))
  D  <- as.matrix(stats::dist(xy))
  xy$w <- rnorm(n)
  xy$z <- as.numeric(t(chol(exp(-D / a) + diag(1e-8, n))) %*% rnorm(n)) +
    rnorm(n, sd = sqrt(nugget))
  sf::st_as_sf(xy, coords = c("x", "y"), crs = 32632)
}
# The first count of the window at which Moran's z is finite.
flc_first_computable <- function(d) min(d$eval_ks[is.finite(d$moran_z[d$eval_ks])])
flc_source <- function() {
  root <- testthat::test_path("..", "..")
  if (file.exists(file.path(root, "DESCRIPTION")) && dir.exists(file.path(root, "man")))
    normalizePath(root)
  else NULL
}


test_that("a call that names both model variables keeps the criterion it names, or the default", {
  cl <- flc_clustered(11, seed = 1)
  geo <- determine_optimal_levels(cl, max_levels = 24)
  expect_identical(as.vector(geo), c(11L, 10L, 12L))
  expect_identical(attr(geo, "structured"), TRUE)
  # Both variables with no criterion, and with "geometric" named: the
  # coordinates-only result, to the attribute, with nothing about Moran's z.
  # Both used to return 10 11 12 under "combined", with its diagnostics and
  # the ranking warning.
  for (crit in list(NULL, "geometric")) {
    lines <- capture_spatialkit_log(out <- flc_catch(do.call(
      determine_optimal_levels,
      c(list(cl, max_levels = 24, response_var = "z", predictor_vars = "w"),
        if (!is.null(crit)) list(criterion = crit)))))
    expect_identical(out$value, geo)
    expect_null(attr(out$value, "diagnostics"))
    expect_length(out$warnings, 0L)
    # A line says the variables are not read and what "combined" adds.
    expect_true(log_has(lines, paste0("criterion = 'geometric' reads the coordinates alone, ",
                                      "so `response_var` and `predictor_vars` are not read")))
    expect_true(log_has(lines, "criterion = 'combined' adds Moran's I of the cell means"))
    expect_false(log_has(lines, "using combined criterion"))
    # Logged at INFO, as the help says: not a caution, so not at WARN.
    expect_true(log_has(lines, "^INFO .*are not read"))
    warn_lines <- capture_spatialkit_log(suppressWarnings(do.call(
      determine_optimal_levels,
      c(list(cl, max_levels = 24, response_var = "z", predictor_vars = "w"),
        if (!is.null(crit)) list(criterion = crit)))), level = logger::WARN)
    expect_false(log_has(warn_lines, "are not read"))
  }
  # A response alone kept "geometric" before and keeps it; nothing to say.
  lines <- capture_spatialkit_log(r1 <- determine_optimal_levels(cl, max_levels = 24,
                                                                 response_var = "z"))
  expect_identical(r1, geo)
  expect_false(log_has(lines, "are not read"))
  # Named, "combined" runs: its diagnostics say so, and no line says the
  # variables are not read.
  lines <- capture_spatialkit_log(rc <- flc_catch(determine_optimal_levels(
    cl, max_levels = 24, response_var = "z", predictor_vars = "w", criterion = "combined")))
  expect_identical(attr(rc$value, "diagnostics")$criterion, "combined")
  expect_false(log_has(lines, "are not read"))
  # A factor response is not read under "geometric", so it is not refused
  # there; the criteria that read it refuse it as before.
  cl$zf <- cut(cl$z, 3, labels = c("low", "mid", "high"))
  expect_identical(determine_optimal_levels(cl, max_levels = 24, response_var = "zf",
                                            predictor_vars = "w"), geo)
  expect_error(determine_optimal_levels(cl, max_levels = 24, response_var = "zf",
                                        predictor_vars = "w", criterion = "combined"),
               "must be numeric or logical")
  # A column that is not there is still an error under every criterion.
  expect_error(determine_optimal_levels(cl, max_levels = 24, response_var = "zz",
                                        predictor_vars = "w"),
               "column 'zz' not found", fixed = TRUE)
})


test_that("on a smooth field \"combined\" does not put first the lowest count Moran's z is computed at", {
  # 600 uniform points with an autocorrelated field: no elbow, so the window
  # is ranked under on_no_elbow = "warn" only.  The window is 2 to 10 and
  # Moran's z is computed at 10 alone, so the rank average put 10 first:
  # 10 6 5 where "geometric" returns 6 5 7.
  pts <- flc_field(n = 600, seed = 12)
  expect_warning(geo <- determine_optimal_levels(pts, max_levels = 40, on_no_elbow = "warn"),
                 class = "spatialkit_no_elbow")
  expect_identical(as.vector(geo), c(6L, 5L, 7L))
  out <- flc_catch(determine_optimal_levels(pts, max_levels = 40, response_var = "z",
                                            predictor_vars = "w", criterion = "combined",
                                            on_no_elbow = "warn"))
  d <- attr(out$value, "diagnostics")
  expect_identical(d$criterion, "combined")
  fc <- flc_first_computable(d)
  expect_identical(fc, 10L)
  expect_false(out$value[1L] == fc)
  # The sag is flat with no elbow, so the window comes out in the
  # geometric order, around the chord rule's count.
  expect_identical(as.vector(out$value), as.vector(geo))
  expect_identical(out$value[1L], d$knee_k)
  expect_length(unique(d$combined_rank), 1L)
  # The only warning is the one that says there is no elbow.
  expect_identical(out$classes, "spatialkit_no_elbow")

  # Clusters at the default: the elbow, 11, comes first, not 10, the first
  # count Moran's z is computed at (10 11 12 on both layers before).  The
  # counts after it follow by sag, 12 then 10, where "geometric" gives
  # 11 10 12.
  for (seed in 1:2) {
    cl <- flc_clustered(11, seed = seed)
    rc <- flc_catch(determine_optimal_levels(cl, max_levels = 24, response_var = "z",
                                             predictor_vars = "w", criterion = "combined"))
    dc <- attr(rc$value, "diagnostics")
    expect_identical(dc$knee_k, 11L, info = paste("seed", seed))
    expect_identical(flc_first_computable(dc), 10L, info = paste("seed", seed))
    expect_identical(rc$value[1L], 11L, info = paste("seed", seed))
    expect_identical(as.vector(rc$value), c(11L, 12L, 10L), info = paste("seed", seed))
    expect_identical(attr(rc$value, "structured"), TRUE, info = paste("seed", seed))
    expect_length(rc$warnings, 0L)
  }
})


test_that("the ranking of \"combined\" is the rank of the sag, and Moran's z does not move it", {
  # Twelve separated clusters (test-review3-S1-resolution.R): the elbow is
  # 12 and its window 8 to 16.
  set.seed(1)
  K <- 12; n <- 720
  cen <- cbind(rep(c(0, 6000, 12000, 18000), 3), rep(c(0, 7000, 14000), each = 4))
  cl  <- rep(seq_len(K), length.out = n)
  d <- data.frame(x = 5e5 + cen[cl, 1] + rnorm(n, 0, 250),
                  y = 5e6 + cen[cl, 2] + rnorm(n, 0, 250), p = rnorm(n))
  d$z <- d$p + rnorm(n)
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  call <- function() suppressWarnings(determine_optimal_levels(
    pts, max_levels = 30, response_var = "z", predictor_vars = "p", criterion = "combined"))
  k <- call()
  dg <- attr(k, "diagnostics")
  expect_identical(dg$knee_k, 12L)
  expect_identical(k[1L], 12L)
  # combined_rank is the rank of the sag of the sweep's WSS curve, the one
  # the elbow is read from.
  sag <- spatialkit:::.elbow_read(seq_along(dg$wss), dg$wss)$sag[dg$eval_ks]
  expect_identical(names(dg$combined_rank), as.character(dg$eval_ks))
  expect_equal(unname(dg$combined_rank), rank(-sag, na.last = TRUE, ties.method = "average"))
  # Moran's z set to favour 10 strongly (0 there, 50 + k elsewhere): the
  # statistic is computed and kept, and the ranking does not move.  With
  # Moran's z in the ranking the same call returned 12 13 11 against
  # 12 11 13.
  zfake <- function(xy, response, predictors, cluster_ids, ...) {
    kk <- length(unique(cluster_ids))
    if (kk <= 9L) return(c(I = NA_real_, z = NA_real_))
    c(I = 0, z = if (kk == 10L) 0 else 50 + kk)
  }
  local_mocked_bindings(.morans_i_for_k = zfake)
  km <- call()
  dm <- attr(km, "diagnostics")
  expect_identical(dm$moran_z[10L], 0)
  expect_identical(dm$moran_z[12L], 62)
  expect_identical(as.vector(km), as.vector(k))
  expect_identical(dm$combined_rank, dg$combined_rank)
})


test_that("\"combined\" puts the elbow first, whatever the re-run clustering finds there", {
  # Loose clusters, where the k-means re-run at the elbow lands in a worse
  # optimum than the sweep did.  Ranked by the re-run's sag, eighteen
  # clusters gave 16 17 18, the elbow, 15, left out, and fifteen gave
  # 15 13 16 where the elbow is 14.
  for (cs in list(c(18L, 15L), c(15L, 14L))) {
    info <- paste(cs[1], "clusters")
    cl <- flc_clustered(cs[1], seed = 1, sd = 40)
    geo <- determine_optimal_levels(cl, max_levels = 24)
    expect_identical(geo[1L], cs[2], info = info)
    rc <- flc_catch(determine_optimal_levels(cl, max_levels = 24, response_var = "z",
                                             predictor_vars = "w", criterion = "combined"))
    d <- attr(rc$value, "diagnostics")
    expect_identical(d$knee_k, cs[2], info = info)
    expect_identical(rc$value[1L], d$knee_k, info = info)
    # combined_rank is the rank of the sweep's sag, whatever wss_eval holds.
    sag <- spatialkit:::.elbow_read(seq_along(d$wss), d$wss)$sag[d$eval_ks]
    expect_equal(unname(d$combined_rank), rank(-sag, na.last = TRUE, ties.method = "average"),
                 info = info)
    expect_length(rc$warnings, 0L)
  }

  # The same on any platform: the re-run's WSS at the elbow made a hundred
  # times the sweep's.  The sweep reads `wss` of the result and the re-run
  # `km$tot.withinss`, so only wss_eval moves, and the ranking does not.
  cl <- flc_clustered(11, seed = 1)
  call <- function() flc_catch(determine_optimal_levels(
    cl, max_levels = 24, response_var = "z", predictor_vars = "w", criterion = "combined"))
  plain <- call()
  real <- spatialkit:::.partition_centres
  local_mocked_bindings(.partition_centres = function(xy, k, ...) {
    kb <- real(xy, k, ...)
    if (!is.null(kb) && k == 11L) kb$km$tot.withinss <- 100 * kb$km$tot.withinss
    kb
  })
  worse <- call()
  dp <- attr(plain$value, "diagnostics")
  dw <- attr(worse$value, "diagnostics")
  expect_gt(dw$wss_eval[11L], 50 * dw$wss[11L])
  expect_identical(dw$wss, dp$wss)
  expect_identical(as.vector(worse$value), as.vector(plain$value))
  expect_identical(worse$value[1L], 11L)
  expect_identical(dw$combined_rank, dp$combined_rank)
  expect_length(worse$warnings, 0L)
})


test_that("where Moran's I can be computed at no count, \"combined\" still ranks by the sag", {
  # Eleven clusters at max_levels = 30: the window is 7 to 15.  The
  # selection half of select_on = "split" leaves too few cells holding a
  # row for Moran's I at any count of it, and so does a constant response.
  # "combined" fell back to the geometric ranking there, without
  # diagnostics, 11 10 12 against 11 12 10 for the whole layer, so the
  # response or the split decided which order came back.
  cl <- flc_clustered(11, seed = 1)
  cl$const <- 1
  call <- function(...) {
    lines <- capture_spatialkit_log(out <- flc_catch(determine_optimal_levels(
      cl, max_levels = 30, predictor_vars = "w", ...)), level = logger::WARN)
    c(out, list(log = lines))
  }
  all <- call(response_var = "z", criterion = "combined")
  da <- attr(all$value, "diagnostics")
  expect_identical(as.vector(all$value), c(11L, 12L, 10L))
  expect_true(any(is.finite(da$moran_z)))
  for (case in list(list(response_var = "z", select_on = "split"),
                    list(response_var = "const"))) {
    info <- paste(names(case), unlist(case), sep = " = ", collapse = ", ")
    r <- do.call(call, c(case, criterion = "combined"))
    d <- attr(r$value, "diagnostics")
    expect_identical(as.vector(r$value), as.vector(all$value), info = info)
    expect_identical(d$criterion, "combined", info = info)
    expect_identical(d$eval_ks, da$eval_ks, info = info)
    expect_true(all(is.na(d$moran_i)), info = info)
    expect_true(all(is.na(d$moran_z)), info = info)
    expect_identical(d$combined_rank, da$combined_rank, info = info)
    expect_true(log_has(r$log, "Moran's I could not be computed at any k from 10 to 15"),
                info = info)
    expect_true(log_has(r$log, paste0("criterion = 'combined' ranks the counts by the sag ",
                                      "as usual, with moran_i and moran_z NA")), info = info)
    expect_false(log_has(r$log, "falling back to geometric"), info = info)
    expect_length(r$warnings, 0L)
  }
  # The split changes the rows the diagnostics read and adds its attribute.
  sp <- call(response_var = "z", select_on = "split", criterion = "combined")
  spl <- attr(sp$value, "split")
  expect_true(all(c("selection", "estimation") %in% names(spl)))
  expect_lt(length(spl$selection), nrow(cl))
  # "morans_i" has nothing to rank by there and falls back, as before.
  rm <- call(response_var = "const", criterion = "morans_i")
  expect_identical(as.vector(rm$value), c(11L, 10L, 12L))
  expect_null(attr(rm$value, "diagnostics"))
  expect_true(log_has(rm$log, paste0("Moran's I could not be computed at any k from 10 to 15 ",
                                     "in the elbow's neighbourhood .*; falling back to ",
                                     "geometric")))
  # It falls back there, once, and does not go on to the ranking.
  expect_false(log_has(rm$log, "no candidate produced a finite Moran's I"))
})


test_that("the ranking warning comes with \"morans_i\" alone and names a way on that runs on the default profile", {
  catch_rank <- function(...) {
    out <- flc_catch(determine_optimal_levels(..., on_no_elbow = "warn"))
    hit <- vapply(out$warnings, inherits, logical(1), "spatialkit_diagnostic_ranking")
    list(value = out$value, n = sum(hit), w = if (any(hit)) out$warnings[[which(hit)[1L]]])
  }
  pts <- flc_field(seed = 701)
  r <- catch_rank(pts, max_levels = 40, response_var = "z", criterion = "morans_i")
  expect_identical(r$n, 1L)
  expect_named(r$w, c("message", "call", "k", "criterion", "first_computable", "window"),
               ignore.order = TRUE)
  msg <- conditionMessage(r$w)
  expect_match(msg, paste0("Use criterion = \"geometric\" for the elbow's ranking, or ",
                           "resolution_profile() and select_resolution() for a count chosen ",
                           "for summaries of the response."), fixed = TRUE)
  expect_false(grepl("criterion = \"cp\"", msg, fixed = TRUE))
  # The way on runs on the default profile, where criterion = "cp" stops.
  prof <- suppressWarnings(resolution_profile(pts, "z"))
  expect_no_error(sel <- suppressWarnings(select_resolution(prof)))
  expect_true(is.finite(sel$best))
  expect_error(suppressWarnings(select_resolution(prof, criterion = "cp")), "`cp` is NA")
  # With both variables the advice is the same: "geometric" no longer reads
  # them, so nothing has to be left out.
  rb <- catch_rank(pts, max_levels = 40, response_var = "z", predictor_vars = "w",
                   criterion = "morans_i")
  expect_identical(rb$n, 1L)
  mb <- conditionMessage(rb$w)
  expect_match(mb, "Use criterion = \"geometric\" for the elbow's ranking", fixed = TRUE)
  expect_false(grepl("leave out", mb, ignore.case = TRUE))
  expect_false(grepl("is read as", mb, fixed = TRUE))
  # "combined" ranks by the sag and raises no such warning.
  expect_identical(catch_rank(pts, max_levels = 40, response_var = "z", predictor_vars = "w",
                              criterion = "combined")$n, 0L)
})


test_that("the help and the vignette say \"combined\" ranks by the sag and the criterion stands", {
  src <- flc_source()
  if (is.null(src)) skip("the help's source is not here to read")
  read <- function(...) gsub("[[:space:]]+", " ",
                             paste(readLines(file.path(src, ...), warn = FALSE), collapse = " "))
  rd <- read("man", "determine_optimal_levels.Rd")
  expect_false(grepl("upgrade", rd, fixed = TRUE))
  expect_false(grepl("rank-average", rd, fixed = TRUE))
  expect_false(grepl("rank average", rd, fixed = TRUE))
  expect_match(rd, paste0("\\code{\"combined\"} (the candidates around the elbow in the order ",
                          "of the WSS curve's log-log sag"), fixed = TRUE)
  expect_match(rd, paste0("Supplying both \\code{response_var} and \\code{predictor_vars} does ",
                          "not change the criterion"), fixed = TRUE)
  expect_match(rd, "Moran's I does not take part in the ranking", fixed = TRUE)
  # The sag is the sweep's, so the elbow comes first, and a Moran's I that
  # cannot be computed changes the diagnostics, not the ranking.
  expect_match(rd, paste0("The sag is read off the sweep's WSS curve, the one the elbow is ",
                          "read from, so with an elbow the first count is the elbow"),
               fixed = TRUE)
  expect_match(rd, "its WSS is reported as \\code{wss_eval} and does not rank", fixed = TRUE)
  expect_match(rd, paste0("\\code{moran_i} and \\code{moran_z} are \\code{NA}, the log ",
                          "gives the reason at the WARN level, and the counts are ranked as ",
                          "they would be otherwise"), fixed = TRUE)
  expect_false(grepl("the geometric ranking is returned without diagnostics", rd, fixed = TRUE))
  expect_false(grepl("whether the geometric ranking is returned in its place", rd, fixed = TRUE))
  # The old average's reason: with no elbow the counts it could not score
  # ranked last, which is not the |z| being near 0.
  expect_match(rd, "the counts below the floor took the last place on its axis", fixed = TRUE)
  # The class list in ?resolution_profile.
  rp <- read("man", "resolution_profile.Rd")
  start <- "\\item{\\code{spatialkit_diagnostic_ranking}}"
  at <- regexpr(start, rp, fixed = TRUE)
  expect_true(at > 0L)
  rest <- substring(rp, at + nchar(start))
  item <- substring(rest, 1L, regexpr("\\item{", rest, fixed = TRUE) - 1L)
  expect_match(item, "under \\code{criterion = \"morans_i\"}:", fixed = TRUE)
  expect_false(grepl("combined", item, fixed = TRUE))
  expect_false(grepl("upgraded", item, fixed = TRUE))
  # The vignette no longer says "combined" ranks by Moran's z.
  vig <- read("vignettes", "resolution.Rmd")
  expect_false(grepl("`criterion = \"combined\"` settings rank", vig, fixed = TRUE))
  expect_match(vig, "`criterion = \"combined\"` keeps the elbow's ranking of those counts",
               fixed = TRUE)
  # The README's entry on the log line where Moran's I is computable nowhere
  # says what each criterion does next.
  for (file in c("README.Rmd", "README.md")) {
    rm <- read(file)
    expect_match(rm, paste0("Under `criterion = \"morans_i\"` the call then falls back to the ",
                            "geometric ranking; under `\"combined\"` it ranks the counts as ",
                            "usual, with `moran_i` and `moran_z` `NA` in the diagnostics."),
                 fixed = TRUE, info = file)
  }
})
