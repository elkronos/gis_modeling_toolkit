# BUG-20: estimate_sac_range()'s floor of 30 counts distinct locations, not
# only rows.
#
# Rows that share a location are one place to the variogram (a pair of them
# is at no lag) and to its test for spatial structure (made on the location
# means).  Five locations visited 20 times passed a floor of 30 rows, got a
# variogram of ten pairs of locations, and were judged by a test of five
# means, which called most such layers unstructured even where the range was
# three times the extent; summarize_by_cell(deff = "variogram") then gave
# every cell a design effect of 1.  Below 30 locations the result is the
# unclassed NA whose reason gives the count, and nothing is fitted.

# `sites` locations on a 1000 m square, `visits` rows at each (a number, or
# one per site).  The location effect has variance 1 and effective range R;
# the variance between visits is 0.3; `w` is a predictor.
fvl_layer <- function(seed, sites, visits, R = 300, crs = 32632) {
  set.seed(seed)
  sx <- runif(sites, 0, 1000); sy <- runif(sites, 0, 1000)
  eff <- drop(crossprod(chol(exp(-3 * as.matrix(stats::dist(cbind(sx, sy))) / R) +
                               diag(1e-8, sites)), rnorm(sites)))
  m <- if (length(visits) == 1L) rep(visits, sites) else visits
  site <- rep(seq_len(sites), m)
  d <- data.frame(x = 5e5 + sx[site], y = 5e6 + sy[site], site = site,
                  w = rnorm(length(site)))
  d$z <- eff[site] + 0.5 * d$w + rnorm(nrow(d), sd = sqrt(0.3))
  sf::st_as_sf(d, coords = c("x", "y"), crs = crs)
}
fvl_quiet <- function(expr)
  logger::with_log_threshold(expr, threshold = logger::FATAL,
                             namespace = "spatialkit", index = 2)
fvl_reason <- function(k, n, where = "with a finite value to model")
  sprintf(paste0("%d distinct location(s) among the %d point(s) %s, fewer than ",
                 "the 30 a variogram range is estimated from"), k, n, where)
# The early NA: unclassed, with its reason and nothing else.
fvl_expect_floor <- function(r, k, n, where = "with a finite value to model") {
  expect_true(is.na(r))
  expect_false(inherits(r, "sac_range"))
  expect_identical(names(attributes(r)), "rejected_reason")
  expect_identical(attr(r, "rejected_reason"), fvl_reason(k, n, where))
}
# Warnings kept and muffled, messages muffled.
fvl_catch <- function(expr) {
  ws <- list()
  value <- withCallingHandlers(expr,
    warning = function(w) { ws[[length(ws) + 1L]] <<- w; invokeRestart("muffleWarning") },
    message = function(m) invokeRestart("muffleMessage"))
  list(value = value, warnings = ws,
       classes = vapply(ws, function(w) class(w)[1L], ""),
       text = vapply(ws, conditionMessage, ""))
}

test_that("five locations visited 20 times give the NA with the count, and nothing is fitted", {
  skip_if_not_installed("gstat")
  p <- fvl_layer(101, sites = 5, visits = 20)
  # No variogram is computed and no model fitted: gstat is not reached.
  local_mocked_bindings(variogram = function(...) stop("gstat::variogram() was called"),
                        .package = "gstat")
  lines <- capture_spatialkit_log(
    expect_no_warning(r <- estimate_sac_range(p, "z")), level = logger::WARN)
  fvl_expect_floor(r, 5L, 100L)
  expect_true(log_has(lines, paste0("the 100 point\\(s\\) with a value to model lie at 5 ",
                                    "distinct location\\(s\\), fewer than 30")))
  # It is not a verdict that the variable has no spatial structure.
  expect_false(spatialkit:::.sac_no_structure(r))
  # The structure test is not asked either: switching it off changes nothing.
  fvl_expect_floor(fvl_quiet(estimate_sac_range(p, "z", structure_alpha = 1)), 5L, 100L)
  # Nor do the binning, the predictors or their detrending.
  fvl_expect_floor(fvl_quiet(estimate_sac_range(p, "z", lags = "equal")), 5L, 100L)
  fvl_expect_floor(fvl_quiet(estimate_sac_range(p, "z", predictor_vars = "w")), 5L, 100L)
  # The order of the rows changes nothing.
  set.seed(2); perm <- sample(nrow(p))
  expect_identical(fvl_quiet(estimate_sac_range(p[perm, ], "z")), r)
})

test_that("the floor is 30 locations: 29 are refused, 30 are fitted", {
  skip_if_not_installed("gstat")
  fvl_expect_floor(fvl_quiet(estimate_sac_range(fvl_layer(7, 29, 2, R = 600), "z")), 29L, 58L)
  for (v in 1:2) {
    r <- fvl_quiet(estimate_sac_range(fvl_layer(8, 30, v, R = 600), "z", structure_alpha = 1))
    # A classed result, refused or not: the variogram was fitted.
    expect_s3_class(r, "sac_range")
    expect_identical(attr(r, "n_used"), 30L * v)
    expect_false(isTRUE(grepl("distinct location", attr(r, "rejected_reason"))))
  }
  # Unequal visits: 29 locations hold 1 to 5 rows each.
  set.seed(3); m <- sample(1:5, 29, replace = TRUE)
  fvl_expect_floor(fvl_quiet(estimate_sac_range(fvl_layer(9, 29, m), "z")), 29L, sum(m))
})

test_that("the floors on rows come first, and the locations counted are those with a value", {
  skip_if_not_installed("gstat")
  # 20 rows at 5 locations: the floor on rows, as before.
  p20 <- fvl_layer(11, 5, 4)
  expect_identical(attr(fvl_quiet(estimate_sac_range(p20, "z")), "rejected_reason"),
                   "20 points, fewer than the 30 a variogram range is estimated from")
  # 40 locations visited 3 times, the response missing at every row of 25 of
  # them: 45 rows with a value, at 15 locations.
  p <- fvl_layer(12, 40, 3)
  p$z[p$site <= 25] <- NA
  fvl_expect_floor(fvl_quiet(estimate_sac_range(p, "z")), 15L, 45L)
  # A predictor missing at the rows of 20 locations leaves 20 with a value.
  q <- fvl_layer(13, 40, 3)
  q$w[q$site <= 20] <- NA
  fvl_expect_floor(fvl_quiet(estimate_sac_range(q, "z", predictor_vars = "w")), 20L, 60L)
  # 40 rows at one location: the locations, not the extent, are what is short.
  set.seed(1)
  one <- sf::st_as_sf(data.frame(x = rep(5e5, 40), y = rep(5e6, 40), z = rnorm(40)),
                      coords = c("x", "y"), crs = 32632)
  fvl_expect_floor(fvl_quiet(estimate_sac_range(one, "z")), 1L, 40L)
})

test_that("the locations are counted in the n_max subsample, and the reason says where they fell short", {
  skip_if_not_installed("gstat")
  # 30 locations visited 10 times: a subsample of 30 rows misses some.
  p <- fvl_layer(14, 30, 10, R = 600)
  r <- fvl_quiet(estimate_sac_range(p, "z", n_max = 30))
  expect_true(is.na(r))
  expect_false(inherits(r, "sac_range"))
  expect_identical(names(attributes(r)), "rejected_reason")
  expect_match(attr(r, "rejected_reason"),
               paste0("^[0-9]+ distinct location\\(s\\) among the 30 point\\(s\\) of the ",
                      "n_max subsample, fewer than the 30"))
  k <- as.integer(sub(" .*", "", attr(r, "rejected_reason")))
  expect_true(k < 30L)
  # The log line says that a larger n_max keeps more of them.
  lines <- capture_spatialkit_log(estimate_sac_range(p, "z", n_max = 30), level = logger::WARN)
  expect_true(log_has(lines, "A larger n_max keeps more of the layer's locations"))
  # Without the cut every location is there.
  expect_s3_class(fvl_quiet(estimate_sac_range(p, "z", structure_alpha = 1)), "sac_range")
  # 10 locations visited 10 times: the layer itself is short, so the reason
  # gives the layer's counts and does not point at n_max.
  p10 <- fvl_layer(17, 10, 10)
  fvl_expect_floor(fvl_quiet(estimate_sac_range(p10, "z", n_max = 30)), 10L, 100L)
  lines <- capture_spatialkit_log(estimate_sac_range(p10, "z", n_max = 30), level = logger::WARN)
  expect_false(log_has(lines, "n_max"))
})

test_that("detrend = 'reml' on fewer than 30 locations fits nothing and warns of nothing", {
  skip_if_not_installed("gstat")
  skip_if_not_installed("nlme")
  p <- fvl_layer(15, 20, 3, R = 600)
  local_mocked_bindings(.reml_trend = function(...) stop(".reml_trend() was called"))
  expect_no_warning(r <- fvl_quiet(
    estimate_sac_range(p, "z", predictor_vars = "w", detrend = "reml")))
  fvl_expect_floor(r, 20L, 60L)
})

test_that("a lon/lat layer is counted after projection, with the same answer", {
  skip_if_not_installed("gstat")
  p <- fvl_layer(16, 5, 20)
  ll <- sf::st_transform(p, 4326)
  fvl_expect_floor(fvl_quiet(estimate_sac_range(ll, "z")), 5L, 100L)
})

test_that("the callers read the floor as they read the other early refusals", {
  skip_if_not_installed("gstat")
  p <- fvl_layer(101, sites = 5, visits = 20)
  why <- fvl_reason(5L, 100L)
  # make_folds(auto_range = TRUE): geometric blocks, and the warning names it.
  expect_warning(fvl_quiet(make_folds(p, k = 2, method = "block_kfold", auto_range = TRUE,
                                      response_var = "z", seed = 1)),
                 why, fixed = TRUE)
  # summarize_by_cell(deff = "variogram"): Kish, with the fallback warning,
  # where the no-structure verdict gave a design effect of 1 in every cell.
  bb <- sf::st_bbox(p); xy <- sf::st_coordinates(p)
  p$cell_id <- 1L + (xy[, 1] > mean(bb[c(1, 3)])) + 2L * (xy[, 2] > mean(bb[c(2, 4)]))
  res <- fvl_catch(fvl_quiet(summarize_by_cell(p, "z", id_col = "cell_id",
                                               deff = "variogram", estimate = "mean")))
  expect_true("spatialkit_deff_fallback" %in% res$classes)
  expect_true(any(grepl(why, res$text, fixed = TRUE)))
  expect_false("spatialkit_deff_no_structure" %in% res$classes)
  s <- res$value
  expect_identical(unique(s$deff_method), "kish (variogram fallback)")
  expect_identical(attr(s, "deff_applied")$method, "kish")
  expect_true(all(attr(s, "deff_applied")$deff > 1))
  # kriging_adequacy(): nothing to krige with, and the error says why.
  cells <- sf::st_sf(cell_id = 1:4, geometry = sf::st_make_grid(sf::st_as_sfc(bb), n = c(2, 2)))
  expect_error(fvl_quiet(suppressWarnings(kriging_adequacy(p, "z", cells, id_col = "cell_id"))),
               why, fixed = TRUE)
  # The Bayesian fit's check before fitting: the range is not estimated.
  chk <- fvl_quiet(spatialkit:::.gp_prefit_range_check(
    p, "z", character(0), resolution = c(x = 50, y = 50), gp_c = 1.5, S = c(500, 500),
    xy_scale = 1))
  expect_identical(chk$status, "not_estimated")
  expect_identical(chk$reason, why)
  # resolution_profile() on a present ladder: Cp's noise variance from the
  # finest level, said as for any variogram that carries no model.
  res <- fvl_catch(fvl_quiet(resolution_profile(p, "z", obs_per_cell = NULL, sac = NULL,
                                                estimate = NULL)))
  sac <- attr(res$value, "sac")
  expect_true(is.na(sac))
  expect_identical(attr(sac, "rejected_reason"), why)
  cpn <- res$text[res$classes == "spatialkit_cp_noise"]
  expect_length(cpn, 1L)
  expect_match(cpn, "carries no variogram model", fixed = TRUE)
  expect_false(grepl("no spatial structure", cpn, fixed = TRUE))
})

test_that("a predictor missing on every row raises no raw-response warning below 30 locations", {
  skip_if_not_installed("gstat")
  # No trend can be fitted, so the raw response would be modelled, with an R
  # warning that says so.  Below 30 locations nothing is modelled at all: the
  # NA's reason is the answer, and a warning of a fit would be untrue.
  p <- fvl_layer(101, sites = 5, visits = 20)
  p$w <- NA_real_
  for (dt in c("ols", "reml")) {
    if (dt == "reml") skip_if_not_installed("nlme")
    res <- fvl_catch(fvl_quiet(estimate_sac_range(p, "z", predictor_vars = "w", detrend = dt)))
    expect_length(res$warnings, 0L)
    fvl_expect_floor(res$value, 5L, 100L)
  }
  # Nor when the n_max subsample is what falls short.
  s <- fvl_layer(14, 30, 10, R = 600)
  s$w <- NA_real_
  res <- fvl_catch(fvl_quiet(estimate_sac_range(s, "z", predictor_vars = "w", n_max = 30)))
  expect_length(res$warnings, 0L)
  expect_match(attr(res$value, "rejected_reason"), "of the n_max subsample, fewer than the 30",
               fixed = TRUE)
  # make_folds(auto_range = TRUE) gives one warning for the one event.
  res <- fvl_catch(fvl_quiet(make_folds(p, k = 2, method = "block_kfold", auto_range = TRUE,
                                        response_var = "z", predictor_vars = "w", seed = 1)))
  expect_length(res$text, 1L)
  expect_match(res$text, fvl_reason(5L, 100L), fixed = TRUE)
  expect_false(any(grepl("RAW response", res$text, fixed = TRUE)))
  # 40 locations: the raw response is modelled, and the warning says so.
  q <- fvl_layer(102, 40, 3)
  q$w <- NA_real_
  res <- fvl_catch(fvl_quiet(estimate_sac_range(q, "z", predictor_vars = "w")))
  expect_length(res$text, 1L)
  expect_match(res$text, paste0("the OLS detrending on w failed (none of the 120 row(s) with a ",
                                "response has every predictor present and finite); the ",
                                "variogram is fitted to the RAW response instead"),
               fixed = TRUE)
  expect_s3_class(res$value, "sac_range")
  expect_false(attr(res$value, "detrended"))
  expect_identical(attr(res$value, "n_used"), 120L)
})

test_that("one location is one pair of coordinates, equal to the last bit", {
  skip_if_not_installed("gstat")
  # What the floor counts is what the structure test and the short-lag pairs
  # count (.sac_locations()): rows at one location have both coordinates
  # equal, after the projection, and nothing nearer than that is merged.  A
  # rounding, or a tolerance relative to the extent or to the coordinates,
  # would refuse these layers or give the floor another count.
  one_row_each <- function(x, y, seed) {
    set.seed(seed)
    sf::st_as_sf(data.frame(x = x, y = y, z = rnorm(length(x))),
                 coords = c("x", "y"), crs = 32632)
  }
  fitted_at_30 <- function(p) {
    r <- fvl_quiet(estimate_sac_range(p, "z", structure_alpha = 1))
    # A classed result, refused or not: the variogram was computed.
    expect_s3_class(r, "sac_range")
    expect_identical(attr(r, "n_used"), 30L)
    # No two rows share a location for the short-lag pairs either.
    expect_identical(attr(r, "zero_lag_pairs"), 0L)
    expect_false(isTRUE(grepl("distinct location", attr(r, "rejected_reason"))))
  }
  # 15 locations on a 1000 m square, each with a twin 1 mm to its east.
  set.seed(21); bx <- runif(15, 0, 1000); by <- runif(15, 0, 1000)
  fitted_at_30(one_row_each(5e5 + c(bx, bx + 0.001), 5e6 + c(by, by), seed = 31))
  # 30 locations inside a 1 m square.
  set.seed(22)
  fitted_at_30(one_row_each(5e5 + runif(30), 5e6 + runif(30), seed = 32))
  # 29 locations and one more row 1e-6 m from one of them: 30 locations.
  set.seed(23); x <- 5e5 + runif(29, 0, 1000); y <- 5e6 + runif(29, 0, 1000)
  fitted_at_30(one_row_each(c(x, x[1] + 1e-6), c(y, y[1]), seed = 33))
  # The same row on the location itself: 29.
  fvl_expect_floor(fvl_quiet(estimate_sac_range(
    one_row_each(c(x, x[1]), c(y, y[1]), seed = 33), "z")), 29L, 30L)
  # x = 0 and x = -0 are one coordinate.
  set.seed(24); x0 <- c(0, runif(28, 1, 1000)); y0 <- runif(29, 0, 1000)
  fvl_expect_floor(fvl_quiet(estimate_sac_range(
    one_row_each(c(x0, -0), c(y0, y0[1]), seed = 34), "z")), 29L, 30L)
})

test_that("the help and NEWS say what the floor counts, and when", {
  root <- testthat::test_path("..", "..")
  flat <- function(x)
    paste(unlist(lapply(x, function(el) if (is.character(el)) el else flat(el))), collapse = "")
  file <- file.path(root, "man", "estimate_sac_range.Rd")
  db <- if (file.exists(file)) tools::parse_Rd(file)
        else tools::Rd_db("spatialkit")[["estimate_sac_range.Rd"]]
  rd <- gsub("[[:space:]]+", " ", flat(db))
  says <- function(text, x) expect_true(grepl(x, text, fixed = TRUE), info = x)
  # What makes two rows one location.
  says(rd, paste0("Two rows are at one location when both their coordinates are equal to ",
                  "the last bit after the projection, so repeat visits recorded with ",
                  "slightly different coordinates count as separate locations; round the ",
                  "coordinates first to pool them."))
  # When the reason names the subsample: the layer has 30 locations or more.
  says(rd, paste0("When the layer itself has 30 locations or more and the subsample holds ",
                  "fewer, the reason of the NA says \"of the n_max subsample\""))
  says(rd, "when the layer has fewer than 30, the reason gives the layer's counts.")
  says(rd, "when the layer has 30 or more locations and the n_max subsample kept fewer")
  expect_false(grepl("when the layer has 30 locations and", rd, fixed = TRUE))
  # The subsample is part of the answer where it decides the floor.
  says(rd, paste0("On a layer whose rows crowd onto few locations the subsample can also ",
                  "decide whether 30 distinct locations remain, and so whether a range is ",
                  "estimated at all"))
  # No warning of a raw-response fit that is not made.
  says(rd, paste0("On fewer than 30 distinct locations nothing is fitted then either: the ",
                  "result is the NA, and there is no such warning."))

  news_file <- file.path(root, "NEWS.md")
  if (!file.exists(news_file)) news_file <- system.file("NEWS.md", package = "spatialkit")
  if (!nzchar(news_file) || !file.exists(news_file)) skip("NEWS.md is not there to read")
  lines <- readLines(news_file, warn = FALSE, encoding = "UTF-8")
  top <- grep("^# ", lines)
  dev <- lines[seq.int(top[1L], top[2L] - 1L)]
  entry_of <- function(section, next_section, anchor) {
    sec <- dev[seq.int(grep(paste0("^## ", section), dev), grep(paste0("^## ", next_section), dev) - 1L)]
    starts <- grep("^\\* ", sec)
    at <- grep(anchor, sec, fixed = TRUE)
    expect_length(at, 1L)
    i <- max(starts[starts <= at])
    j <- min(c(starts[starts > i], length(sec) + 1L)) - 1L
    gsub("[[:space:]]+", " ", paste(sec[i:j], collapse = " "))
  }
  # Breaking changes compares with 2.0.0, which had no kriging_adequacy()
  # and gave a one-location layer a bare NA.
  brk <- entry_of("Breaking changes", "New features",
                  "`estimate_sac_range()` needs 30 distinct locations as well as 30 rows")
  says(brk, paste0("A layer of 30 or more rows at 2 to 29 distinct locations, to which 2.0.0 ",
                   "fitted a variogram, gets the unclassed `NA`"))
  expect_false(grepl("kriging_adequacy", brk, fixed = TRUE))
  expect_false(grepl("no extent", brk, fixed = TRUE))
  fix <- entry_of("Bug fixes", "Documentation",
                  "`estimate_sac_range()` counted rows, not locations, against its floor of")
  says(fix, "When the layer has 30 or more locations and the `n_max` subsample keeps fewer")
  says(fix, "`kriging_adequacy()` stops with it")
  says(fix, "in place of \"the points have no extent\"")
  says(fix, paste0("Where a predictor is missing on every row of such a layer, there is no ",
                   "warning that the variogram is fitted to the raw response, since none is ",
                   "fitted."))
})
