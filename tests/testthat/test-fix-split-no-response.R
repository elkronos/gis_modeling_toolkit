# tests/testthat/test-fix-split-no-response.R
# ---------------------------------------------------------------------------
# resolution_profile(select_on = "split") with an estimate, on a layer whose
# rows with a response all lie in the estimation half: `loo` and the
# estimate's other columns read no row and are NA, with a warning of class
# spatialkit_no_scored_rows that gives the reason.  The check that
# `estimate` has rows to read is made on the layer, so the call is not an
# error; it used to return the NA columns with nothing said.
# ---------------------------------------------------------------------------

fsn_catch <- function(expr) {
  ws <- list()
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws[[length(ws) + 1L]] <<- w
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  list(value = value, warnings = ws,
       classes = vapply(ws, function(w) class(w)[1L], character(1)))
}
fsn_of <- function(r, class) Filter(function(w) inherits(w, class), r$warnings)
# The first warning of the class, or an empty list when there is none, so
# that the expectations on its text and fields fail rather than stop.
fsn_first <- function(r, class) {
  w <- fsn_of(r, class)
  if (length(w)) w[[1L]] else list(message = NA_character_)
}
fsn_msg <- function(w) if (inherits(w, "condition")) conditionMessage(w) else NA_character_

# 400 uniform points with a response that has structure; the split of a
# layer depends on its geometry and the seed alone.  This one splits 200/200;
# 600 points from seed 1 split into halves of different sizes.
fsn_layer <- function(n = 400, seed = 43) {
  set.seed(seed)
  x <- stats::runif(n, 0, 1e5); y <- stats::runif(n, 0, 1e5)
  sf::st_as_sf(data.frame(x = x, y = y,
                          z = sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(n, 0, 1.5)),
               coords = c("x", "y"), crs = 32617)
}
fsn_split <- function(pts)
  attr(suppressWarnings(resolution_profile(pts, "z", estimate = NULL, obs_per_cell = NULL,
                                           n_levels = 3, nstart = 2, partition = "kmeans",
                                           select_on = "split")), "split")
# The layer with no response on its selection half.
fsn_no_selection <- function(n = 400, seed = 43) {
  pts <- fsn_layer(n, seed)
  sp  <- fsn_split(pts)
  pts$z[sp$selection] <- NA
  list(pts = pts, split = sp)
}
fsn_text <- function(selection, n_response)
  sprintf(paste0("resolution_profile(): no row of the selection half has a finite value of ",
                 "'z', so `loo`, `loo_se`, `scored_per_cell`, `share_alone` and `own_weight` ",
                 "are NA at every level, and so are `rss`, `cp`, `moran_i` and `moran_z`, ",
                 "which also read the selection half. Under select_on = \"split\" the response ",
                 "is read on the selection half alone (%d rows); the %s in the estimation ",
                 "half."),
          selection,
          if (n_response == 1L) "1 row with a finite value is"
          else sprintf("%d rows with a finite value are all", n_response))
.FSN_COLS <- c("scored_per_cell", "loo", "loo_se", "share_alone", "own_weight")

test_that("a selection half with no response gives NA columns and a warning that says why", {
  skip_if_not_installed("FNN")
  L <- fsn_no_selection()
  n_sel <- length(L$split$selection)
  n_resp <- sum(is.finite(L$pts$z))
  expect_identical(n_resp, length(L$split$estimation))
  calls <- list(
    default = function() resolution_profile(L$pts, "z", select_on = "split"),
    pooled  = function() resolution_profile(L$pts, "z", estimate = "pooled",
                                            partition = "kmeans", nstart = 2,
                                            select_on = "split"),
    mean    = function() resolution_profile(L$pts, "z", estimate = "mean",
                                            select_on = "split"),
    ladder  = function() resolution_profile(L$pts, "z", estimate = "pooled",
                                            obs_per_cell = NULL, n_levels = 5,
                                            partition = "kmeans", nstart = 2,
                                            select_on = "split"))
  for (nm in names(calls)) {
    r <- fsn_catch(calls[[nm]]())
    expect_s3_class(r$value, "resolution_profile")
    expect_length(fsn_of(r, "spatialkit_no_scored_rows"), 1L)
    w <- fsn_first(r, "spatialkit_no_scored_rows")
    expect_identical(fsn_msg(w), fsn_text(n_sel, n_resp), info = nm)
    expect_identical(w$reason, "no response in the selection half", info = nm)
    expect_identical(w$n, 0L, info = nm)
    expect_identical(w$selection, n_sel, info = nm)
    expect_identical(w$n_response, n_resp, info = nm)
    expect_identical(w$column, "z", info = nm)
    p <- r$value
    for (col in .FSN_COLS) expect_true(all(is.na(p[[col]])), info = paste(nm, col))
    for (col in c("rss", "cp", "moran_i", "moran_z"))
      expect_true(all(is.na(p[[col]])), info = paste(nm, col))
    e <- attr(p, "estimate")
    expect_identical(e$n, 0L, info = nm)
    expect_identical(e$rows, "selection half", info = nm)
    expect_identical(e$loo_one_cell, NA_real_, info = nm)
    out <- utils::capture.output(print(p))
    expect_true(sprintf("  estimate    : %s; loo NA: no row of the selection half has a response",
                        if (nm == "mean") "plain means" else "pooled toward the 8 nearest cells")
                %in% out, info = nm)
    expect_false(any(grepl("fewer than two rows", out, fixed = TRUE)), info = nm)
  }
})

test_that("the warning counts the rows with a response, and a non-finite value is none", {
  skip_if_not_installed("FNN")
  L <- fsn_no_selection()
  n_sel <- length(L$split$selection)
  # One row of the layer with a response, in the estimation half: the text
  # says "1 row ... is", beside the size's own warning.
  one <- L$pts; one$z[-L$split$estimation[1L]] <- NA
  r <- fsn_catch(resolution_profile(one, "z", select_on = "split"))
  expect_identical(r$classes, c("spatialkit_not_supported", "spatialkit_no_scored_rows"))
  w <- fsn_first(r, "spatialkit_no_scored_rows")
  expect_identical(fsn_msg(w), fsn_text(n_sel, 1L))
  expect_identical(w$n_response, 1L)
  # Inf on the selection half is no response there.
  inf <- L$pts; inf$z[L$split$selection] <- Inf
  r <- fsn_catch(resolution_profile(inf, "z", select_on = "split"))
  expect_length(fsn_of(r, "spatialkit_no_scored_rows"), 1L)
  expect_identical(fsn_msg(fsn_first(r, "spatialkit_no_scored_rows")),
                   fsn_text(n_sel, length(L$split$estimation)))
  expect_identical(attr(r$value, "estimate")$n, 0L)
})

test_that("on halves of different sizes the warning gives each its own count", {
  skip_if_not_installed("FNN")
  L <- fsn_no_selection(600, 1)
  n_sel <- length(L$split$selection)
  n_est <- length(L$split$estimation)
  # The halves differ, so one half's count cannot pass for the other's.
  expect_false(n_sel == n_est)
  n_resp <- sum(is.finite(L$pts$z))
  expect_identical(n_resp, n_est)
  calls <- list(
    default = function() resolution_profile(L$pts, "z", select_on = "split"),
    mean    = function() resolution_profile(L$pts, "z", estimate = "mean",
                                            select_on = "split"))
  for (nm in names(calls)) {
    r <- fsn_catch(calls[[nm]]())
    expect_length(fsn_of(r, "spatialkit_no_scored_rows"), 1L)
    w <- fsn_first(r, "spatialkit_no_scored_rows")
    expect_identical(w$selection, n_sel, info = nm)
    expect_identical(w$n_response, n_resp, info = nm)
    expect_identical(fsn_msg(w), fsn_text(n_sel, n_resp), info = nm)
    expect_true(grepl(sprintf("the selection half alone (%d rows); the %d rows", n_sel, n_est),
                      fsn_msg(w), fixed = TRUE), info = nm)
  }
})

test_that("with a variogram given, reliability is computed and the warning does not name it", {
  skip_if_not_installed("FNN")
  skip_if_not_installed("gstat")
  L <- fsn_no_selection()
  v <- suppressWarnings(estimate_sac_range(fsn_layer(), "z"))
  r <- fsn_catch(resolution_profile(L$pts, "z", estimate = "pooled", obs_per_cell = NULL,
                                    n_levels = 4, partition = "kmeans", nstart = 2, sac = v,
                                    select_on = "split"))
  expect_length(fsn_of(r, "spatialkit_no_scored_rows"), 1L)
  msg <- fsn_msg(fsn_first(r, "spatialkit_no_scored_rows"))
  expect_identical(msg, fsn_text(length(L$split$selection), length(L$split$estimation)))
  expect_false(isTRUE(grepl("reliability", msg, fixed = TRUE)))
  p <- r$value
  expect_true(any(is.finite(p$reliability)))
  for (col in c(.FSN_COLS, "rss", "cp", "moran_i", "moran_z"))
    expect_true(all(is.na(p[[col]])), info = col)
})

test_that("the estimate adds its columns to the warning and changes nothing else", {
  skip_if_not_installed("FNN")
  L <- fsn_no_selection()
  call <- function(est) resolution_profile(L$pts, "z", estimate = est, obs_per_cell = NULL,
                                           n_levels = 5, partition = "kmeans", nstart = 2,
                                           select_on = "split")
  base <- fsn_catch(call(NULL))
  # Without an estimate the warning is raised too, and names the columns
  # that read the selection half (BUG-34).
  expect_length(fsn_of(base, "spatialkit_no_scored_rows"), 1L)
  base_w <- fsn_first(base, "spatialkit_no_scored_rows")
  for (est in c("pooled", "mean")) {
    r <- fsn_catch(call(est))
    for (col in names(base$value))
      expect_identical(r$value[[col]], base$value[[col]], info = paste(est, col))
    expect_identical(sum(r$classes == "spatialkit_no_scored_rows"), 1L)
    expect_identical(r$classes, base$classes, info = est)
    # The two warnings differ in the estimate's columns alone.
    w <- fsn_first(r, "spatialkit_no_scored_rows")
    expect_identical(fsn_msg(w), fsn_text(length(L$split$selection), length(L$split$estimation)))
    expect_identical(
      sub(paste0("`loo`, `loo_se`, `scored_per_cell`, `share_alone` and `own_weight` are NA ",
                 "at every level, and so are `rss`, `cp`, `moran_i` and `moran_z`, which also ",
                 "read the selection half"),
          paste0("`rss`, `cp`, `moran_i` and `moran_z`, which read the selection half, are NA ",
                 "at every level"),
          fsn_msg(w), fixed = TRUE),
      fsn_msg(base_w), info = est)
    for (f in c("reason", "n", "selection", "n_response", "column"))
      expect_identical(w[[f]], base_w[[f]], info = paste(est, f))
  }
})

test_that("downstream, the size still sets the count and loo refuses with the reason", {
  skip_if_not_installed("FNN")
  L <- fsn_no_selection()
  pooled <- suppressWarnings(resolution_profile(L$pts, "z", select_on = "split"))
  s <- select_resolution(pooled)
  expect_identical(s$criterion, "size")
  expect_identical(s$best, pooled$levels)
  expect_identical(s$count_rule$scored$n, 0L)
  out <- utils::capture.output(print(s$count_rule))
  expect_true("  scored      : loo NA: no row of the selection half has a response" %in% out)
  expect_false(any(grepl("in NA cells", out, fixed = TRUE)))
  msg <- paste0("select_resolution(): `loo` is NA at every level because no row of the ",
                "selection half has a response: under select_on = \"split\" the profile ",
                "scores the estimate on that half alone.")
  # The error's text, or NA when there is none.
  err <- function(expr) tryCatch({ expr; NA_character_ }, error = conditionMessage)
  expect_identical(err(select_resolution(pooled, criterion = "loo")), msg)
  mean_p <- suppressWarnings(resolution_profile(L$pts, "z", estimate = "mean",
                                                select_on = "split"))
  # "auto" reads a profile of plain means by loo.
  expect_identical(err(select_resolution(mean_p)), msg)
  expect_identical(err(select_resolution(mean_p, criterion = "loo")), msg)
})

test_that("the other criteria refuse with their own reason or the split's, and summary() gives the reason", {
  skip_if_not_installed("FNN")
  L <- fsn_no_selection()
  err <- function(expr) tryCatch({ expr; NA_character_ }, error = conditionMessage)
  pooled <- suppressWarnings(resolution_profile(L$pts, "z", select_on = "split"))
  ladder <- suppressWarnings(resolution_profile(L$pts, "z", estimate = "pooled",
                                                obs_per_cell = NULL, n_levels = 5,
                                                partition = "kmeans", nstart = 2,
                                                select_on = "split"))
  mean_p <- suppressWarnings(resolution_profile(L$pts, "z", estimate = "mean",
                                                select_on = "split"))
  # On the default profile, a ladder of one size with no variogram, cp and
  # reliability stop with the size ladder's reason.  On a ladder of counts
  # cp and moran_z name the selection half (BUG-34), and reliability, which
  # needs a variogram, keeps its own reason.  None names `loo`.
  size_text <- function(cr, what)
    sprintf(paste0("select_resolution(): `%s` is NA at every level: this profile is on a ",
                   "size ladder (obs_per_cell = 4.5), which fits no variogram unless `sac` ",
                   "is given, and `%s` takes its %s from one. Build the profile with ",
                   "`sac = NULL` to have the variogram estimated, with a `sac` of your own, ",
                   "or with obs_per_cell = NULL for the automatic ladder; criterion = ",
                   "\"size\" reads this profile."), cr, cr, what)
  expect_identical(err(select_resolution(pooled, criterion = "cp")),
                   size_text("cp", "noise variance"))
  expect_identical(err(select_resolution(pooled, criterion = "reliability")),
                   size_text("reliability", "range"))
  for (cr in c("cp", "reliability", "moran_z")) {
    e <- err(select_resolution(ladder, criterion = cr, diagnostic = cr == "moran_z"))
    expect_true(startsWith(e, sprintf("select_resolution(): `%s` is NA at every level", cr)),
                info = cr)
    expect_false(grepl("`loo` is NA", e, fixed = TRUE), info = cr)
    expect_identical(grepl("selection half", e, fixed = TRUE), cr != "reliability", info = cr)
  }
  # summary() names the selection half, not a missing response.
  split_text <- paste0("summary.resolution_profile(): no count criterion is finite at any ",
                       "level. No row of the selection half has a response, and under ",
                       "select_on = \"split\" the response is read on that half alone, so no ",
                       "criterion that reads it is finite: cp, loo and moran_z are NA at every ",
                       "level, and reliability needs a usable variogram of that half; elbow, ",
                       "read off the points alone, is NA when the WSS curve has no elbow.")
  expect_identical(err(summary(pooled)),
                   paste0(split_text, " The ladder is one size (obs_per_cell = 4.5); the count ",
                          "is set by the size, and select_resolution(criterion = \"size\") ",
                          "reads it."))
  expect_identical(err(summary(ladder)), split_text)
  e <- err(summary(mean_p))
  expect_true(startsWith(e, paste0(split_text, " The profile is on a ladder of sizes")))
  for (p in list(pooled, ladder, mean_p))
    expect_false(grepl("cp needs a response", err(summary(p)), fixed = TRUE))
  # With one row of the selection half scored, the half has a response, and
  # summary() does not say that it has none.
  one <- fsn_layer(); one$z[L$split$selection[-1L]] <- NA
  p1 <- suppressWarnings(resolution_profile(one, "z", select_on = "split"))
  expect_identical(attr(p1, "estimate")$n, 1L)
  e <- err(summary(p1))
  expect_true(startsWith(e, "summary.resolution_profile(): no count criterion is finite"))
  expect_false(grepl("No row of the selection half", e, fixed = TRUE))
})

test_that("the warning is raised where the selection half has no response, with or without an estimate, and not elsewhere", {
  skip_if_not_installed("FNN")
  pts <- fsn_layer()
  sp  <- fsn_split(pts)
  # Every row has a response.
  r <- fsn_catch(resolution_profile(pts, "z", select_on = "split"))
  expect_length(fsn_of(r, "spatialkit_no_scored_rows"), 0L)
  expect_identical(attr(r$value, "estimate")$n, length(sp$selection))
  # One row of the selection half has one: loo is NA for fewer than two rows,
  # as it was, and the print says so.
  one <- pts; one$z[sp$selection[-1L]] <- NA
  r <- fsn_catch(resolution_profile(one, "z", estimate = "pooled", partition = "kmeans",
                                    nstart = 2, select_on = "split"))
  expect_length(fsn_of(r, "spatialkit_no_scored_rows"), 0L)
  expect_identical(attr(r$value, "estimate")$n, 1L)
  expect_true(any(grepl("one cell: loo NA (fewer than two rows)",
                        utils::capture.output(print(r$value)), fixed = TRUE)))
  # select_resolution(criterion = "loo") on it keeps the reason of a level
  # with no cell of two rows.
  expect_error(select_resolution(r$value, criterion = "loo"),
               "it needs a level at which some cell holds two rows with a response",
               fixed = TRUE)
  # The count record of the size it is read by says that loo read the one
  # row, not that the half has none.
  s <- suppressWarnings(select_resolution(r$value))
  expect_identical(s$count_rule$scored$n, 1L)
  out <- utils::capture.output(print(s$count_rule))
  expect_true(paste0("  scored      : loo NA on the 1 row of the selection half in 1 cell ",
                     "(one cell: NA); 100.0% of rows alone") %in% out)
  expect_false(any(grepl("no row of the selection half has a response", out, fixed = TRUE)))
  # No response anywhere and no split: a named estimate is an error, and
  # "auto" scores none.
  none <- pts; none$z <- NA_real_
  expect_error(resolution_profile(none, "z", estimate = "pooled"),
               "`estimate` needs rows with a value of 'z', and none has a finite one.",
               fixed = TRUE)
  r <- fsn_catch(resolution_profile(none, "z"))
  expect_null(attr(r$value, "estimate"))
  expect_length(fsn_of(r, "spatialkit_no_scored_rows"), 0L)
  # The split with no estimate: nothing reads `loo`, but `rss`, `cp` and the
  # Moran columns read the selection half, and the warning says so (BUG-34).
  L <- fsn_no_selection()
  r <- fsn_catch(resolution_profile(L$pts, "z", estimate = NULL, select_on = "split"))
  expect_length(fsn_of(r, "spatialkit_no_scored_rows"), 1L)
  expect_false(grepl("`loo`", fsn_msg(fsn_first(r, "spatialkit_no_scored_rows")), fixed = TRUE))
})

test_that("a call a guard stops says nothing of the columns", {
  skip_if_not_installed("FNN")
  L <- fsn_no_selection()
  r <- local({
    withr::local_options(spatialkit.fit_work_max = 1)
    fsn_catch(resolution_profile(L$pts, "z", select_on = "split"))
  })
  expect_s3_class(r$value, "spatialkit_work_limit")
  expect_length(fsn_of(r, "spatialkit_no_scored_rows"), 0L)
})

test_that("the help says when the columns are NA and why", {
  root <- testthat::test_path("..", "..")
  flat <- function(x)
    paste(unlist(lapply(x, function(el) if (is.character(el)) el else flat(el))), collapse = "")
  file <- file.path(root, "man", "resolution_profile.Rd")
  db <- if (file.exists(file)) tools::parse_Rd(file)
        else tools::Rd_db("spatialkit")[["resolution_profile.Rd"]]
  rp <- gsub("[[:space:]]+", " ", flat(db))
  says <- function(text, x) expect_true(grepl(x, text, fixed = TRUE), info = x)
  says(rp, paste0("Under select_on = \"split\" a selection half with no finite response is ",
                  "not an error: the five columns are NA at every level, with a warning of ",
                  "class spatialkit_no_scored_rows."))
  says(rp, paste0("All five are NA at every level when no row of the selection half has a ",
                  "response (select_on = \"split\"; spatialkit_no_scored_rows)."))
  says(rp, paste0("spatialkit_no_scored_rowsUnder select_on = \"split\": no row of the ",
                  "selection half has a response"))
  says(rp, paste0("loo and the estimate's other columns read no row and are NA at every ",
                  "level too."))
  says(rp, paste0("The condition carries reason (\"no response in the selection half\"), n ",
                  "(0, the rows of the selection half scored), selection (the rows of the ",
                  "selection half), n_response (the rows of the layer with a response) and ",
                  "column."))
  # ?select_resolution says when summary() stops, and that it says why.
  file <- file.path(root, "man", "select_resolution.Rd")
  db <- if (file.exists(file)) tools::parse_Rd(file)
        else tools::Rd_db("spatialkit")[["select_resolution.Rd"]]
  sr <- gsub("[[:space:]]+", " ", flat(db))
  says(sr, paste0("When none is finite, as on a constant response, where loo is NA too, or ",
                  "when no row of the selection half has a response (select_on = \"split\"), ",
                  "summary() stops with an error that says so."))
})
