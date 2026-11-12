# tests/testthat/test-fix-split-no-response-noest.R
# ---------------------------------------------------------------------------
# resolution_profile(select_on = "split") without an estimate, on a layer
# whose rows with a response all lie in the estimation half: `rss`, `cp`,
# `moran_i` and `moran_z` read the selection half and are NA at every level.
# The warning of class spatialkit_no_scored_rows that gives the reason was
# raised only with an estimate, and summary() and select_resolution() gave
# the general reasons ("cp needs a response").  The split now records the
# rows of the selection half with a response (`n_response`), and the warning,
# summary() and select_resolution() give the reason with or without an
# estimate.
# ---------------------------------------------------------------------------

fsn2_catch <- function(expr) {
  ws <- list()
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws[[length(ws) + 1L]] <<- w
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  list(value = value, warnings = ws,
       classes = vapply(ws, function(w) class(w)[1L], character(1)))
}
fsn2_of <- function(r, class) Filter(function(w) inherits(w, class), r$warnings)
fsn2_first <- function(r, class) {
  w <- fsn2_of(r, class)
  if (length(w)) w[[1L]] else list(message = NA_character_)
}
fsn2_msg <- function(w) if (inherits(w, "condition")) conditionMessage(w) else NA_character_
fsn2_err <- function(expr) tryCatch({ expr; NA_character_ }, error = conditionMessage)

# 400 uniform points with a response that has structure, and a predictor.
fsn2_layer <- function(n = 400, seed = 43) {
  set.seed(seed)
  x <- stats::runif(n, 0, 1e5); y <- stats::runif(n, 0, 1e5)
  sf::st_as_sf(data.frame(x = x, y = y,
                          z = sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(n, 0, 1.5),
                          w = stats::rnorm(n)),
               coords = c("x", "y"), crs = 32617)
}
fsn2_split <- function(pts)
  attr(suppressWarnings(resolution_profile(pts, "z", estimate = NULL, obs_per_cell = NULL,
                                           n_levels = 3, nstart = 2, partition = "kmeans",
                                           select_on = "split")), "split")
fsn2_profile <- function(pts, ..., estimate = NULL)
  resolution_profile(pts, "z", estimate = estimate, obs_per_cell = NULL, n_levels = 5,
                     partition = "kmeans", nstart = 2, select_on = "split", ...)

# The texts, in full.
fsn2_warn_text <- function(selection, n_response, pred = FALSE)
  sprintf(paste0("resolution_profile(): no row of the selection half has a finite value of ",
                 "'z'%s, so `rss`, `cp`, `moran_i` and `moran_z`, which read the selection ",
                 "half, are NA at every level. Under select_on = \"split\" the response is ",
                 "read on the selection half alone (%d rows); the %d rows with %s are all in ",
                 "the estimation half."),
          if (pred) " and of every predictor" else "", selection, n_response,
          if (pred) "them" else "a finite value")
fsn2_est_text <- function(selection, n_response)
  sprintf(paste0("resolution_profile(): no row of the selection half has a finite value of ",
                 "'z', so `loo`, `loo_se`, `scored_per_cell`, `share_alone` and `own_weight` ",
                 "are NA at every level, and so are `rss`, `cp`, `moran_i` and `moran_z`, ",
                 "which also read the selection half. Under select_on = \"split\" the response ",
                 "is read on the selection half alone (%d rows); the %d rows with a finite ",
                 "value are all in the estimation half."),
          selection, n_response)
fsn2_select_text <- function(cr, pred = FALSE)
  sprintf(paste0("select_resolution(): `%s` is NA at every level because no row of the ",
                 "selection half has a response%s: under select_on = \"split\" the profile ",
                 "reads the response on that half alone."),
          cr, if (pred) " and every predictor" else "")
fsn2_summary_text <- function(cols = "cp and moran_z", pred = FALSE)
  sprintf(paste0("summary.resolution_profile(): no count criterion is finite at any level. ",
                 "No row of the selection half has a response%s, and under select_on = ",
                 "\"split\" the response is read on that half alone, so no criterion that ",
                 "reads it is finite: %s are NA at every level, and reliability needs a ",
                 "usable variogram of that half; elbow, read off the points alone, is NA when ",
                 "the WSS curve has no elbow."),
          if (pred) " and every predictor" else "", cols)
fsn2_general_cp <- paste0("select_resolution(): `cp` is NA at every level (it needs a response, ",
                          "and a level whose cells hold two scored rows each on average for ",
                          "its noise variance when the variogram gives no nugget).")

test_that("without an estimate, a selection half with no response gives the warning", {
  pts <- fsn2_layer()
  sp  <- fsn2_split(pts)
  pts$z[sp$selection] <- NA
  n_sel  <- length(sp$selection)
  n_resp <- sum(is.finite(pts$z))
  expect_identical(n_resp, length(sp$estimation))
  calls <- list(ladder = function() fsn2_profile(pts),
                auto   = function() resolution_profile(pts, "z", estimate = NULL,
                                                       select_on = "split"))
  for (nm in names(calls)) {
    r <- fsn2_catch(calls[[nm]]())
    expect_s3_class(r$value, "resolution_profile")
    expect_identical(sum(r$classes == "spatialkit_no_scored_rows"), 1L, info = nm)
    w <- fsn2_first(r, "spatialkit_no_scored_rows")
    expect_identical(fsn2_msg(w), fsn2_warn_text(n_sel, n_resp), info = nm)
    expect_identical(w$reason, "no response in the selection half", info = nm)
    expect_identical(w$n, 0L, info = nm)
    expect_identical(w$selection, n_sel, info = nm)
    expect_identical(w$n_response, n_resp, info = nm)
    expect_identical(w$column, "z", info = nm)
    p <- r$value
    for (col in c("rss", "cp", "moran_i", "moran_z"))
      expect_true(all(is.na(p[[col]])), info = paste(nm, col))
    expect_false("loo" %in% names(p))
    expect_null(attr(p, "estimate"))
    expect_identical(attr(p, "split")$n_response, 0L, info = nm)
    expect_identical(attr(p, "bounds")$n_response, n_resp, info = nm)
  }
})

test_that("without an estimate, summary() and select_resolution() give the split's reason", {
  pts <- fsn2_layer()
  sp  <- fsn2_split(pts)
  pts$z[sp$selection] <- NA
  p <- suppressWarnings(fsn2_profile(pts))
  expect_identical(fsn2_err(summary(p)), fsn2_summary_text())
  expect_false(grepl("cp needs a response", fsn2_err(summary(p)), fixed = TRUE))
  expect_identical(fsn2_err(select_resolution(p, criterion = "cp")), fsn2_select_text("cp"))
  expect_identical(fsn2_err(select_resolution(p, criterion = "moran_z", diagnostic = TRUE)),
                   fsn2_select_text("moran_z"))
  # "auto" without an estimate reads cp.
  expect_identical(fsn2_err(select_resolution(p)), fsn2_select_text("cp"))
  # reliability and elbow keep their own reasons: neither reads the response
  # of that half alone.
  expect_identical(fsn2_err(select_resolution(p, criterion = "reliability")),
                   paste0("select_resolution(): `reliability` is NA at every level (it needs ",
                          "a usable variogram with an identified range)."))
  expect_false(grepl("selection half", fsn2_err(select_resolution(p, criterion = "elbow")),
                     fixed = TRUE))
})

test_that("with predictors, a row counts when it has a response and every predictor", {
  # The response missing on the selection half.
  pts <- fsn2_layer()
  sp  <- fsn2_split(pts)
  a <- pts; a$z[sp$selection] <- NA
  # The response on every row, the predictor missing on the selection half.
  b <- pts; b$w[sp$selection] <- NA
  for (nm in c("response", "predictor")) {
    lay <- if (nm == "response") a else b
    r <- fsn2_catch(fsn2_profile(lay, predictor_vars = "w"))
    expect_identical(sum(r$classes == "spatialkit_no_scored_rows"), 1L, info = nm)
    w <- fsn2_first(r, "spatialkit_no_scored_rows")
    expect_identical(fsn2_msg(w),
                     fsn2_warn_text(length(sp$selection), length(sp$estimation), pred = TRUE),
                     info = nm)
    expect_identical(w$n_response, length(sp$estimation), info = nm)
    p <- r$value
    for (col in c("rss", "cp", "moran_i", "moran_z"))
      expect_true(all(is.na(p[[col]])), info = paste(nm, col))
    expect_identical(attr(p, "split")$n_response, 0L, info = nm)
    expect_identical(fsn2_err(summary(p)), fsn2_summary_text(pred = TRUE), info = nm)
    expect_identical(fsn2_err(select_resolution(p, criterion = "cp")),
                     fsn2_select_text("cp", pred = TRUE), info = nm)
    expect_identical(fsn2_err(select_resolution(p, criterion = "moran_z", diagnostic = TRUE)),
                     fsn2_select_text("moran_z", pred = TRUE), info = nm)
  }
})

test_that("with an estimate the warning is as it was, and cp and moran_z give the reason", {
  skip_if_not_installed("FNN")
  pts <- fsn2_layer()
  sp  <- fsn2_split(pts)
  pts$z[sp$selection] <- NA
  n_sel <- length(sp$selection); n_resp <- length(sp$estimation)
  for (est in c("pooled", "mean")) {
    r <- fsn2_catch(fsn2_profile(pts, estimate = est))
    expect_identical(sum(r$classes == "spatialkit_no_scored_rows"), 1L, info = est)
    w <- fsn2_first(r, "spatialkit_no_scored_rows")
    expect_identical(fsn2_msg(w), fsn2_est_text(n_sel, n_resp), info = est)
    expect_identical(w$n, 0L, info = est)
    expect_identical(w$n_response, n_resp, info = est)
    p <- r$value
    expect_identical(attr(p, "split")$n_response, 0L, info = est)
    expect_identical(attr(p, "estimate")$n, 0L, info = est)
    expect_identical(fsn2_err(summary(p)), fsn2_summary_text("cp, loo and moran_z"), info = est)
    expect_identical(fsn2_err(select_resolution(p, criterion = "cp")), fsn2_select_text("cp"),
                     info = est)
    expect_identical(fsn2_err(select_resolution(p, criterion = "moran_z", diagnostic = TRUE)),
                     fsn2_select_text("moran_z"), info = est)
  }
})

test_that("where the split is not the reason, nothing changes", {
  pts <- fsn2_layer()
  sp  <- fsn2_split(pts)
  # No row of the layer has a response: no warning, and the general texts.
  none <- pts; none$z <- NA_real_
  r <- fsn2_catch(fsn2_profile(none))
  expect_length(fsn2_of(r, "spatialkit_no_scored_rows"), 0L)
  expect_identical(attr(r$value, "split")$n_response, 0L)
  expect_identical(fsn2_err(select_resolution(r$value, criterion = "cp")), fsn2_general_cp)
  expect_true(grepl("cp needs a response", fsn2_err(summary(r$value)), fixed = TRUE))
  # Every row has a response: the split records the rows of its selection
  # half, and nothing is said.
  r <- fsn2_catch(fsn2_profile(pts))
  expect_length(fsn2_of(r, "spatialkit_no_scored_rows"), 0L)
  expect_identical(attr(r$value, "split")$n_response, length(sp$selection))
  # Geometry alone: the split records no response count.
  r <- fsn2_catch(resolution_profile(pts, estimate = NULL, obs_per_cell = NULL, n_levels = 5,
                                     partition = "kmeans", nstart = 2, select_on = "split"))
  expect_length(fsn2_of(r, "spatialkit_no_scored_rows"), 0L)
  expect_null(attr(r$value, "split")$n_response)
  # One row of the selection half has a response: the half has one, so
  # nothing is said, as with an estimate.
  one <- pts; one$z[sp$selection[-1L]] <- NA
  r <- fsn2_catch(fsn2_profile(one))
  expect_length(fsn2_of(r, "spatialkit_no_scored_rows"), 0L)
  expect_identical(attr(r$value, "split")$n_response, 1L)
  expect_false(grepl("No row of the selection half", fsn2_err(summary(r$value)), fixed = TRUE))
  expect_false(grepl("selection half", fsn2_err(select_resolution(r$value, criterion = "cp")),
                     fixed = TRUE))
  # Without the split: no `split` attribute and no warning.
  r <- fsn2_catch(resolution_profile(none, "z", estimate = NULL, obs_per_cell = NULL,
                                     n_levels = 5, partition = "kmeans", nstart = 2))
  expect_length(fsn2_of(r, "spatialkit_no_scored_rows"), 0L)
  expect_null(attr(r$value, "split"))
})

test_that("the rows of the selection half are counted on the layer, not on the subsample", {
  # Three rows of the selection half have a response, and the subsample of
  # 300 holds none of them: the half has rows with a response, so the
  # warning is not raised, and the split records three.
  set.seed(5)
  n <- 3000
  d <- data.frame(x = stats::runif(n, 0, 5000), y = stats::runif(n, 0, 5000))
  d$z <- sin(d$x / 700) + cos(d$y / 900) + stats::rnorm(n, sd = 0.3)
  big <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32615)
  sb <- attr(suppressWarnings(resolution_profile(big, "z", estimate = NULL, select_on = "split",
                                                 n_levels = 4, sample_n = 300, seed = 123)),
             "split")
  big$z[setdiff(sb$selection, sb$selection[1:3])] <- NA
  r <- fsn2_catch(resolution_profile(big, "z", estimate = NULL, select_on = "split",
                                     n_levels = 4, sample_n = 300, seed = 123))
  expect_identical(attr(r$value, "bounds")$n_scored, 0L)
  expect_identical(attr(r$value, "split")$n_response, 3L)
  expect_length(fsn2_of(r, "spatialkit_no_scored_rows"), 0L)
})

test_that("a profile whose split carries no count reads the estimate's record", {
  skip_if_not_installed("FNN")
  pts <- fsn2_layer()
  sp  <- fsn2_split(pts)
  pts$z[sp$selection] <- NA
  p <- suppressWarnings(fsn2_profile(pts, estimate = "pooled"))
  old <- attr(p, "split"); old$n_response <- NULL; attr(p, "split") <- old
  expect_identical(fsn2_err(summary(p)), fsn2_summary_text("cp, loo and moran_z"))
  expect_identical(fsn2_err(select_resolution(p, criterion = "cp")), fsn2_select_text("cp"))
})

test_that("the help says the warning is raised with or without an estimate", {
  root <- testthat::test_path("..", "..")
  flat <- function(x)
    paste(unlist(lapply(x, function(el) if (is.character(el)) el else flat(el))), collapse = "")
  rd <- function(name) {
    file <- file.path(root, "man", name)
    db <- if (file.exists(file)) tools::parse_Rd(file) else tools::Rd_db("spatialkit")[[name]]
    gsub("[[:space:]]+", " ", flat(db))
  }
  says <- function(text, x) expect_true(grepl(x, text, fixed = TRUE), info = x)
  rp <- rd("resolution_profile.Rd")
  says(rp, paste0("spatialkit_no_scored_rowsUnder select_on = \"split\": ",
                  "no row of the selection half has a response, and the rows of the layer ",
                  "with one are all in the estimation half. rss, cp, moran_i and moran_z, ",
                  "which read the selection half, are NA at every level."))
  says(rp, paste0("With predictor_vars a row counts when it has a response and every ",
                  "predictor. Not raised when no row of the layer has a response."))
  says(rp, "n (0, the rows of the selection half scored)")
  says(rp, paste0("and with a response n_response, the rows of the selection half with a ",
                  "finite response and every predictor, counted on the layer before any ",
                  "subsample."))
  # The item does not limit the warning to a profile with an estimate.
  expect_false(grepl("spatialkit_no_scored_rowsWith an estimate", rp, fixed = TRUE))
  expect_false(grepl("spatialkit_no_scored_rowsresolution_profile() with an estimate", rp,
                     fixed = TRUE))
  sr <- rd("select_resolution.Rd")
  says(sr, paste0("or when no row of the selection half has a response (select_on = ",
                  "\"split\"), summary() stops with an error that says so. In that last case ",
                  "\"cp\" and \"moran_z\", and \"loo\" with an estimate, are errors that give ",
                  "the same reason, except \"cp\" on a size ladder that fitted no variogram, ",
                  "which gives the size ladder's reason, as above."))
})

test_that("on a size ladder that fitted no variogram, cp gives the size ladder's reason", {
  # The size ladder's refusal comes first for `cp` (and "auto", which reads
  # it without an estimate); `moran_z`, which needs no variogram, gives the
  # split's reason.  The help and NEWS say so.
  pts <- fsn2_layer()
  sp  <- fsn2_split(pts)
  pts$z[sp$selection] <- NA
  r <- fsn2_catch(resolution_profile(pts, "z", estimate = NULL, obs_per_cell = c(5, 10),
                                     partition = "kmeans", nstart = 2, select_on = "split"))
  p <- r$value
  expect_length(fsn2_of(r, "spatialkit_no_scored_rows"), 1L)
  expect_identical(attr(p, "split")$n_response, 0L)
  expect_identical(attr(p, "cp_noise")$source, "none")
  for (cr in c("cp", "auto")) {
    e <- fsn2_err(select_resolution(p, criterion = cr))
    expect_true(grepl(paste0("select_resolution(): `cp` is NA at every level: this profile ",
                             "is on a size ladder"), e, fixed = TRUE), info = cr)
    expect_false(grepl("selection half", e, fixed = TRUE), info = cr)
  }
  expect_identical(fsn2_err(select_resolution(p, criterion = "moran_z", diagnostic = TRUE)),
                   fsn2_select_text("moran_z"))
})

test_that("one row of the layer with a response, in the estimation half, is the split's case", {
  # The layer has a response (one row), so summary() and select_resolution()
  # give the split's reason and not the general "needs a response".
  pts <- fsn2_layer()
  sp  <- fsn2_split(pts)
  pts$z <- NA_real_
  pts$z[sp$estimation[1L]] <- 2
  r <- fsn2_catch(fsn2_profile(pts))
  p <- r$value
  expect_length(fsn2_of(r, "spatialkit_no_scored_rows"), 1L)
  expect_identical(attr(p, "split")$n_response, 0L)
  expect_identical(attr(p, "bounds")$n_response, 1L)
  expect_identical(fsn2_err(summary(p)), fsn2_summary_text())
  expect_identical(fsn2_err(select_resolution(p, criterion = "cp")), fsn2_select_text("cp"))
  expect_identical(fsn2_err(select_resolution(p, criterion = "auto")), fsn2_select_text("cp"))
  expect_identical(fsn2_err(select_resolution(p, criterion = "moran_z", diagnostic = TRUE)),
                   fsn2_select_text("moran_z"))
})
