# ---------------------------------------------------------------------------
# The texts on Cp's noise variance in resolution_profile() and
# select_resolution() (BUG-38, BUG-162, BUG-155).
#
# Without a nugget Cp takes as its noise variance the residual mean square of
# the finest level whose cells hold at least two scored rows on average.
# Where no level's cells hold that many, Cp is NA at every level; where that
# level's residual mean square is 0, as on a layer whose locations each carry
# repeat rows of one value, Cp is NA too.  The warning, the print and
# select_resolution() gave the first reason for the second case, and the
# warnings and log lines of an unusable variogram said Cp takes its noise
# variance "from the finest level" with no condition.  The rule is unchanged;
# these tests pin the texts.
# ---------------------------------------------------------------------------

# 40 locations on a 1000 m square, each held by `reps` rows of one value.  At
# 40 cells each cell holds one location, so with three rows a location its
# residual sum of squares is exactly 0.
.cnt_layer <- function(reps = 3L, n = 40L, seed = 11) {
  set.seed(seed)
  loc <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000))
  loc$z <- sin(loc$x / 200) + cos(loc$y / 300)
  sf::st_as_sf(loc[rep(seq_len(n), each = reps), ], coords = c("x", "y"), crs = 32632)
}

# A sac made by hand, with the attributes resolution_profile() reads as
# estimate_sac_range() sets them.  `range = NA` with a `reason` is a refused
# range; `structure = FALSE` leaves the model a pure nugget.
.cnt_sac <- function(range = 250, nugget = 0, psill = 1, reason = NULL, converged = NULL,
                     at_bound = NULL, short_lag = NULL, structure = TRUE) {
  vm <- if (structure)
    data.frame(model = c("Nug", "Exp"), psill = c(nugget, psill),
               range = c(0, if (is.na(range)) 3000 else range / 3), stringsAsFactors = FALSE)
  else data.frame(model = "Nug", psill = nugget + psill, range = 0, stringsAsFactors = FALSE)
  if (!is.null(converged)) attr(vm, "converged") <- converged
  s <- structure(as.numeric(range), class = c("sac_range", "numeric"), variogram_model = vm,
                 nugget = nugget, detrended = FALSE, crs = sf::st_crs(32632))
  if (!is.null(reason)) {
    attr(s, "rejected_range") <- 9000
    attr(s, "rejected_reason") <- reason
  }
  if (!is.null(at_bound)) attr(s, "nugget_at_bound") <- at_bound
  if (!is.null(short_lag)) attr(s, "short_lag_semivariance") <- short_lag
  s
}

# The value, every R warning as a condition, and the print.
.cnt_run <- function(expr) {
  w <- list()
  val <- withCallingHandlers(expr, warning = function(c) {
    w[[length(w) + 1L]] <<- c
    invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w, text = vapply(w, conditionMessage, ""),
       print = utils::capture.output(print(val)))
}
.cnt_of <- function(r, class) Filter(function(c) inherits(c, class), r$warnings)
.cnt_profile <- function(pts, sac, levels = c(5L, 40L))
  .cnt_run(resolution_profile(pts, "z", sac = sac, levels = levels, estimate = NULL,
                              partition = "kmeans", nstart = 3L, seed = 1))
.cnt_select_error <- function(p, ...)
  tryCatch({ select_resolution(p, ...); NA_character_ }, error = conditionMessage)

# The clause of the fallback, the same in every warning and log line of an
# unusable variogram.
.CNT_CLAUSE <- paste0("`cp` takes its noise variance from the finest level whose cells hold ",
                      "at least two scored rows on average, and is NA at every level when ",
                      "none does.")
.CNT_ZERO <- paste0("the finest level whose cells hold at least two scored rows on average ",
                    "(40 cells, 3 scored rows a cell) has a residual mean square of 0, which ",
                    "as a noise variance would leave `cp` no penalty")
.CNT_PRINT_ZERO <- paste0("cp is NA at every level: the residual mean square at 40 cells (3 ",
                          "scored rows a cell), the finest level whose cells hold two on ",
                          "average, is 0")
.CNT_GENERAL <- paste0("select_resolution(): `cp` is NA at every level (it needs a response, ",
                       "and a level whose cells hold two scored rows each on average for its ",
                       "noise variance when the variogram gives no nugget).")


test_that("with the nugget at its bound, a finest level whose residual mean square is 0 is named, not called missing (BUG-38, BUG-162)", {
  r <- .cnt_profile(.cnt_layer(), .cnt_sac(at_bound = TRUE))
  p <- r$value
  expect_identical(p$levels, c(5L, 40L))
  expect_equal(p$cell_resp_mean, c(24, 3))
  # The level qualifies, and its residual sum of squares is exactly 0.
  expect_identical(p$rss[2L], 0)
  # The rule is unchanged: no noise variance, so Cp is NA at every level.
  expect_true(all(is.na(p$cp)))
  expect_true(is.na(attr(p, "one_cell")$cp))
  w <- .cnt_of(r, "spatialkit_cp_noise")
  expect_length(w, 1L)
  expect_identical(conditionMessage(w[[1L]]),
                   paste0("resolution_profile(): the variogram's nugget (0 against a sill of 1) ",
                          "is at the fit's lower bound, which is not an estimate, and ", .CNT_ZERO,
                          ", so `cp` is NA at every level."))
  expect_false(any(grepl("no level's cells hold", r$text, fixed = TRUE)))
  expect_true(is.na(w[[1L]]$source))
  expect_true(is.na(w[[1L]]$value))
  expect_identical(w[[1L]]$level, 40L)
  expect_identical(w[[1L]]$uncapped_source, "finest-level residual mean square")
  expect_identical(w[[1L]]$uncapped_value, 0)
  # The profile keeps the level that qualified.
  cn <- attr(p, "cp_noise")
  expect_true(is.list(cn))
  expect_identical(cn$source, NA_character_)
  expect_identical(cn$value, NA_real_)
  expect_identical(cn$level, 40L)
  expect_identical(cn$rows_per_cell, 3)
  expect_identical(cn$rms, 0)
  expect_match(cn$reason, .CNT_ZERO, fixed = TRUE)
  # The print says the same.
  expect_true(any(grepl(paste0("  variogram   : nugget 0 (at the fit's lower bound), partial sill 1, ",
                               "range 250; ", .CNT_PRINT_ZERO), r$print, fixed = TRUE)))
  expect_false(any(grepl("no level's cells hold", r$print, fixed = TRUE)))
  # select_resolution() names the reason the profile has, by "cp" and by
  # "auto", which reads a profile without an estimate by Cp.
  want <- paste0("select_resolution(): `cp` is NA at every level: the variogram's nugget (0 ",
                 "against a sill of 1) is at the fit's lower bound, which is not an estimate, ",
                 "and ", .CNT_ZERO, ".")
  expect_identical(.cnt_select_error(p, criterion = "cp"), want)
  expect_identical(.cnt_select_error(p), want)
  # A size ladder's record is still the only one select_resolution() reads
  # as "no variogram".
  expect_false(identical(cn$source, "none"))
})


test_that("where no level's cells hold two scored rows on average the present texts stay (BUG-38 control)", {
  r <- .cnt_profile(.cnt_layer(reps = 1L), .cnt_sac(at_bound = TRUE), levels = c(30L, 40L))
  p <- r$value
  expect_true(all(p$cell_resp_mean < 2))
  expect_true(all(is.na(p$cp)))
  expect_null(attr(p, "cp_noise"))
  w <- .cnt_of(r, "spatialkit_cp_noise")
  expect_length(w, 1L)
  expect_match(conditionMessage(w[[1L]]),
               paste0("no level's cells hold at least two scored rows on average to take a ",
                      "noise variance from, so `cp` is NA at every level."), fixed = TRUE)
  expect_true(is.na(w[[1L]]$level))
  expect_identical(w[[1L]]$uncapped_source, "variogram nugget")
  expect_true(any(grepl("cp is NA at every level: no level's cells hold two scored rows on average",
                        r$print, fixed = TRUE)))
  expect_identical(.cnt_select_error(p, criterion = "cp"), .CNT_GENERAL)
})


test_that("with no usable variogram, a finest level whose residual mean square is 0 is named after the sweep (BUG-38)", {
  pts <- .cnt_layer()
  sacs <- list(
    not_converged = .cnt_sac(range = NA, nugget = 0.5,
                             reason = "fitted range exceeds the largest lag fitted",
                             converged = FALSE),
    shortest_lag = .cnt_sac(range = NA, nugget = 0.5,
                            reason = "fitted range is below the shortest lag fitted"),
    range_alone = 250)
  why <- c(not_converged = "`sac` gives no usable variogram model",
           shortest_lag = "`sac` gives no usable variogram model",
           range_alone = "`sac` carries no variogram model and so no nugget")
  for (k in names(sacs)) {
    r <- .cnt_profile(pts, sacs[[k]])
    p <- r$value
    expect_true(all(is.na(p$cp)), info = k)
    expect_true(is.na(attr(p, "one_cell")$cp), info = k)
    # Before, nothing said why: the unusable warning, where there was one,
    # said `cp` takes its noise variance from that level.
    w <- .cnt_of(r, "spatialkit_cp_noise")
    expect_length(w, 1L)
    expect_identical(conditionMessage(w[[1L]]),
                     paste0("resolution_profile(): ", why[[k]], ", and ", .CNT_ZERO,
                            ", so `cp` is NA at every level. `reliability` needs the ",
                            "variogram's range and is NA."), info = k)
    expect_identical(w[[1L]]$level, 40L, info = k)
    expect_identical(attr(p, "cp_noise")$level, 40L, info = k)
    expect_true(any(grepl(paste0("  variogram   : none usable (reliability is NA; ",
                                 .CNT_PRINT_ZERO, ")"), r$print, fixed = TRUE)), info = k)
    expect_identical(.cnt_select_error(p, criterion = "cp"),
                     paste0("select_resolution(): `cp` is NA at every level: ", why[[k]],
                            ", and ", .CNT_ZERO, "."), info = k)
  }
})


test_that("the warnings of an unusable variogram give the fallback's condition, and cp is NA when no level meets it (BUG-155)", {
  pts <- .cnt_layer(reps = 1L)
  lv <- c(30L, 40L)
  sacs <- list(
    shortest_lag = .cnt_sac(range = NA, nugget = 0.5,
                            reason = "fitted range is below the shortest lag fitted"),
    no_structure = .cnt_sac(range = NA, nugget = 0.5, reason = "no spatial structure detected"),
    nugget_bound = .cnt_sac(range = NA, reason = "fitted range exceeds the largest lag fitted",
                            converged = TRUE, at_bound = TRUE),
    not_converged = .cnt_sac(range = NA, nugget = 0.5,
                             reason = "fitted range exceeds the largest lag fitted",
                             converged = FALSE))
  for (k in names(sacs)) {
    r <- .cnt_profile(pts, sacs[[k]], levels = lv)
    p <- r$value
    expect_true(all(p$cell_resp_mean < 2), info = k)
    expect_true(all(is.na(p$cp)), info = k)
    expect_null(attr(p, "cp_noise"))
    w <- .cnt_of(r, "spatialkit_variogram_unusable")
    expect_length(w, 1L)
    expect_true(endsWith(conditionMessage(w[[1L]]), .CNT_CLAUSE), info = k)
    expect_false(any(grepl("noise variance from the finest level.", r$text, fixed = TRUE)),
                 info = k)
  }
  # With the semivariance of nearest-neighbour pairs the hold is named inside
  # the clause.
  r <- .cnt_profile(pts, .cnt_sac(range = NA, nugget = 0.5, short_lag = 0.4,
                                  reason = "fitted range is below the shortest lag fitted"),
                    levels = lv)
  expect_true(any(endsWith(r$text, paste0(
    "`cp` takes its noise variance from the finest level whose cells hold at least two ",
    "scored rows on average, held to the semivariance of nearest-neighbour pairs (0.4), and ",
    "is NA at every level when none does."))))
  # The two log lines: a range alone, and a pure nugget of no structure.
  for (s in list(250, .cnt_sac(range = NA, reason = "no spatial structure detected",
                               structure = FALSE))) {
    log <- capture_spatialkit_log(invisible(suppressWarnings(
      resolution_profile(pts, "z", sac = s, levels = lv, estimate = NULL,
                         partition = "kmeans", nstart = 3L, seed = 1))))
    line <- grep("`reliability` is NA and `cp` takes", log, value = TRUE)
    expect_length(line, 1L)
    expect_true(endsWith(line, .CNT_CLAUSE))
  }
  # Where a level meets the condition, the fallback is unchanged.
  r <- .cnt_profile(pts, sacs$shortest_lag, levels = c(10L, 30L, 40L))
  expect_identical(attr(r$value, "cp_noise")$source, "finest-level residual mean square")
  expect_identical(attr(r$value, "cp_noise")$level, 10L)
  expect_true(all(is.finite(r$value$cp)))
})


test_that("the help gives both reasons cp has no noise variance, and the record of the second", {
  root <- testthat::test_path("..", "..")
  flat <- function(x)
    paste(unlist(lapply(x, function(el) if (is.character(el)) el else flat(el))), collapse = "")
  file <- file.path(root, "man", "resolution_profile.Rd")
  db <- if (file.exists(file)) tools::parse_Rd(file)
        else tools::Rd_db("spatialkit")[["resolution_profile.Rd"]]
  text <- gsub("[[:space:]]+", " ", flat(db))
  says <- function(x) expect_true(grepl(x, text, fixed = TRUE), info = x)
  says(paste0("When no level's cells hold that many, or the finest that does has a residual ",
              "mean square of 0 to within the rounding of the values scored (as when each of ",
              "its cells holds one location whose rows carry one value), cp is NA at every level"))
  says(paste0("or there is none and cp is NA at every level: the finest level whose cells hold ",
              "two scored rows on average has a residual mean square of 0 to within rounding, ",
              "or, with the nugget at its bound, no level's cells hold that many."))
  says(paste0("when the finest level whose cells hold two scored rows on average has a residual ",
              "mean square of 0 to within rounding, source and value NA, level, that level, ",
              "rows_per_cell, its scored rows a cell, rms, its residual mean square (0, or the ",
              "rounding it was read as), and reason, which select_resolution() gives; NULL when ",
              "cp has none for another reason"))
  says(paste0("is NA, and cp takes its noise variance from the finest level whose cells hold at ",
              "least two scored rows on average (see cp above)."))
})


test_that("with predictors, a finest level whose residual mean square is rounding error is read as 0 (BUG-38, review R-01)", {
  # The OLS residuals of rows of one location and one value differ in their
  # last bits, so the finest level's residual sum of squares is rounding
  # (1.9e-29 on macOS on values of 1) where without predictors it is exactly
  # 0.  Cp took that rounding as its noise variance, had no penalty, and
  # select_resolution() picked the finest level.  The value of the rounding
  # moves with the platform's arithmetic, so the tests pin that it is far
  # below the values and that Cp is NA, not the number.
  pts <- .cnt_layer()
  pts$x1 <- sf::st_coordinates(pts)[, 1L]
  s <- .cnt_sac(at_bound = TRUE)
  attr(s, "detrended") <- TRUE
  r <- .cnt_run(resolution_profile(pts, "z", predictor_vars = "x1", sac = s, levels = c(5L, 40L),
                                   estimate = NULL, partition = "kmeans", nstart = 3L, seed = 1))
  p <- r$value
  expect_identical(attr(p, "variable"), "residuals")
  expect_lt(p$rss[2L], 1e-20)
  expect_true(all(is.na(p$cp)))
  expect_true(is.na(attr(p, "one_cell")$cp))
  cn <- attr(p, "cp_noise")
  expect_identical(cn$source, NA_character_)
  expect_identical(cn$level, 40L)
  expect_identical(cn$rows_per_cell, 3)
  expect_lt(cn$rms, 1e-20)
  # The texts say 0, and where the residuals are not exactly equal, the
  # rounding they hold.
  zero <- paste0("\\(40 cells, 3 scored rows a cell\\) has a residual mean square of 0",
                 "( to within rounding \\([^)]+\\))?, which as a noise variance would leave ",
                 "`cp` no penalty")
  w <- .cnt_of(r, "spatialkit_cp_noise")
  expect_length(w, 1L)
  expect_match(conditionMessage(w[[1L]]),
               paste0("^resolution_profile\\(\\): the variogram's nugget \\(0 against a sill of ",
                      "1\\) is at the fit's lower bound, which is not an estimate, and the finest ",
                      "level whose cells hold at least two scored rows on average ", zero,
                      ", so `cp` is NA at every level\\.$"))
  expect_identical(w[[1L]]$level, 40L)
  expect_identical(w[[1L]]$uncapped_value, cn$rms)
  expect_true(any(grepl(paste0("cp is NA at every level: the residual mean square at 40 cells ",
                               "\\(3 scored rows a cell\\), the finest level whose cells hold two ",
                               "on average, is 0( to within rounding \\([^)]+\\))?$"), r$print)))
  expect_match(.cnt_select_error(p, criterion = "cp"),
               paste0("^select_resolution\\(\\): `cp` is NA at every level: the variogram's ",
                      "nugget .* ", zero, "\\.$"))
})


test_that("a finest level with a small residual mean square that is not rounding keeps its noise variance (BUG-38 control)", {
  # Rows of one location that differ by 1e-6: a mean square of 1e-12, far
  # above the rounding of values of 1.  The fallback takes it, with and
  # without predictors.
  pts <- .cnt_layer()
  pts$z <- pts$z + rep(c(0, 1e-6, -1e-6), 40L)
  pts$x1 <- sf::st_coordinates(pts)[, 1L]
  s <- .cnt_sac(at_bound = TRUE)
  for (pv in list(NULL, "x1")) {
    if (!is.null(pv)) attr(s, "detrended") <- TRUE
    r <- .cnt_run(resolution_profile(pts, "z", predictor_vars = pv, sac = s, levels = c(5L, 40L),
                                     estimate = NULL, partition = "kmeans", nstart = 3L, seed = 1))
    cn <- attr(r$value, "cp_noise")
    expect_identical(cn$source, "finest-level residual mean square", info = pv)
    expect_identical(cn$level, 40L, info = pv)
    expect_equal(cn$value, 1e-12, tolerance = 1e-6, info = pv)
    expect_true(all(is.finite(r$value$cp)), info = pv)
  }
})


test_that("a refused range with the nugget at its bound names the zero level too (BUG-38, review R-02)", {
  r <- .cnt_profile(.cnt_layer(), .cnt_sac(range = NA,
                                           reason = "fitted range exceeds the largest lag fitted",
                                           converged = TRUE, at_bound = TRUE))
  p <- r$value
  expect_true(all(is.na(p$cp)))
  w <- .cnt_of(r, "spatialkit_cp_noise")
  expect_length(w, 1L)
  # The variogram is kept with no correlation function, so the warning says
  # `reliability` is NA as well.
  expect_identical(conditionMessage(w[[1L]]),
                   paste0("resolution_profile(): the variogram's nugget (0 against a sill of 1) ",
                          "is at the fit's lower bound, which is not an estimate, and ", .CNT_ZERO,
                          ", so `cp` is NA at every level. `reliability` needs the variogram's ",
                          "range and is NA."))
  expect_identical(w[[1L]]$level, 40L)
  expect_identical(attr(p, "cp_noise")$level, 40L)
  expect_true(any(grepl(paste0("  variogram   : nugget 0 (at the fit's lower bound), partial sill 1, ",
                               "range unidentified; ", .CNT_PRINT_ZERO), r$print, fixed = TRUE)))
  expect_identical(.cnt_select_error(p, criterion = "cp"),
                   paste0("select_resolution(): `cp` is NA at every level: the variogram's nugget ",
                          "(0 against a sill of 1) is at the fit's lower bound, which is not an ",
                          "estimate, and ", .CNT_ZERO, "."))
})


test_that("the zero case carries the semivariance cap, and moran_z keeps its own refusal (review notes)", {
  r <- .cnt_profile(.cnt_layer(), .cnt_sac(at_bound = TRUE, short_lag = 0.4))
  w <- .cnt_of(r, "spatialkit_cp_noise")
  expect_length(w, 1L)
  expect_identical(w[[1L]]$cap, 0.4)
  expect_identical(w[[1L]]$level, 40L)
  expect_identical(attr(r$value, "cp_noise")$level, 40L)
  # The cp reason is for `cp` alone: a profile whose `moran_z` is NA at every
  # level gets moran_z's own refusal.
  p <- .cnt_profile(.cnt_layer(), .cnt_sac(at_bound = TRUE))$value
  p1 <- p[p$levels == 5L, ]
  expect_true(all(is.na(p1$moran_z)))
  e <- .cnt_select_error(p1, criterion = "moran_z", diagnostic = TRUE)
  expect_match(e, "`moran_z` is NA at every level", fixed = TRUE)
  expect_false(grepl("residual mean square of 0", e, fixed = TRUE))
})
