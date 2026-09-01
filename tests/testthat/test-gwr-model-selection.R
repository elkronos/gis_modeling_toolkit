# =============================================================================
# GWR model selection
#
# The pure helpers and the whole assembly path run without GWmodel, via the
# injectable .engine.  That is deliberate: GWmodel is in Suggests, so anything
# that can only be reached through it is untested on every CRAN check machine.
# The GWmodel-gated tests at the bottom cover only what a stub cannot -- that
# the criterion GWmodel actually returns is oriented the way we assume.
# =============================================================================

# --- .gwr_ms_criterion -------------------------------------------------------

test_that(".gwr_ms_criterion prefers a column named AICc over the positional guess", {
  # AICc deliberately in column 1, so a positional fallback would read the
  # WRONG column and this test would fail.
  m <- matrix(c(10, 20, 30,
                99, 98, 97), ncol = 2)
  colnames(m) <- c("AICc", "RSS")
  res <- .gwr_ms_criterion(m)
  expect_true(res$by_name)
  expect_identical(res$column, 1L)
  expect_equal(res$values, c(10, 20, 30))
  expect_identical(res$label, "AICc")
})

test_that(".gwr_ms_criterion matches the column name case-insensitively and ignores padding", {
  m <- matrix(c(1, 2, 5, 6), ncol = 2)
  colnames(m) <- c("RSS", "  aicc ")
  res <- .gwr_ms_criterion(m)
  expect_true(res$by_name)
  expect_identical(res$column, 2L)
  expect_equal(res$values, c(5, 6))
})

test_that(".gwr_ms_criterion reads AICc from column 3 when the table is unlabelled", {
  # GWmodel's GWR.df is documented AND built as c(bandwidth, AIC, AICc, RSS):
  # in Model.selection.r it is rbind()ed from
  #     c(bw, aic.rss[2], aic.rss[3], aic.rss[1])  with aic.rss = c(RSS, AIC, AICc)
  # over UNNAMED vectors, so it never carries column names -- which makes this
  # positional branch the path every real call takes, not a rare fallback.
  #
  # Column 2 is the UNCORRECTED AIC.  Reading it ranked models on AIC while
  # labelling the answer AICc, which selects larger models: executed against a
  # faithful GWmodel stub, the old code selected the model containing a
  # pure-noise predictor that AICc drops.
  m <- matrix(c(1, 2, 3, 40, 50, 60, 7, 8, 9, 100, 200, 300), ncol = 4)
  res <- .gwr_ms_criterion(m)
  expect_false(res$by_name)
  expect_identical(res$column, 3L)
  expect_equal(res$values, c(7, 8, 9))
  expect_true(res$shape_ok)          # the documented four columns
  expect_match(res$label, "assumed")
})

test_that(".gwr_ms_criterion flags a table that is not the documented shape", {
  # Anything other than c(bandwidth, AIC, AICc, RSS) is a layout this code has
  # never seen; it reads the same position but tells the caller the ranking is
  # unverified, so gwr_model_selection() can warn instead of asserting AICc.
  m <- matrix(c(1, 2, 3, 40, 50, 60, 7, 8, 9), ncol = 3)
  res <- .gwr_ms_criterion(m)
  expect_false(res$by_name)
  expect_identical(res$column, 3L)
  expect_false(res$shape_ok)
})

test_that(".gwr_ms_criterion falls back to column 1 for a single-column table", {
  m <- matrix(c(4, 5, 6), ncol = 1)
  res <- .gwr_ms_criterion(m)
  expect_identical(res$column, 1L)
  expect_equal(res$values, c(4, 5, 6))
  expect_false(res$shape_ok)
})

test_that(".gwr_ms_criterion accepts a bare numeric vector", {
  res <- .gwr_ms_criterion(c(3.5, 1.5, 2.5))
  expect_equal(res$values, c(3.5, 1.5, 2.5))
  expect_false(res$by_name)
})

test_that(".gwr_ms_criterion reads a data.frame column without round-tripping through character", {
  # as.matrix() on a mixed data.frame runs format() over the numeric columns,
  # which silently rounds to 7 significant digits.  Reading the column
  # directly must not lose those digits.
  df <- data.frame(model = c("m1", "m2"),
                   AICc  = c(1.123456789, 2.987654321),
                   stringsAsFactors = FALSE)
  res <- .gwr_ms_criterion(df)
  expect_true(res$by_name)
  expect_identical(res$values, c(1.123456789, 2.987654321))
})

test_that(".gwr_ms_criterion coerces non-finite values to NA", {
  m <- matrix(c(1, 2, 3, 10, Inf, NaN), ncol = 2)
  res <- .gwr_ms_criterion(m)
  expect_equal(res$values, c(10, NA_real_, NA_real_))
})

test_that(".gwr_ms_criterion rejects unusable diagnostic tables", {
  expect_error(.gwr_ms_criterion(NULL), "no diagnostic table")
  expect_error(.gwr_ms_criterion(matrix(numeric(0), nrow = 0, ncol = 0)),
               "no columns")
  expect_error(.gwr_ms_criterion(list(a = 1)), "no columns")
})


# --- .gwr_ms_varsets ---------------------------------------------------------

test_that(".gwr_ms_varsets handles GWmodel's list(DeVar, InDeVars) shape", {
  ml <- list(list("z", "a"),
             list("z", c("a", "b")),
             list("z", list("a", "b", "c")))
  out <- .gwr_ms_varsets(ml, "z")
  expect_identical(out[[1]], "a")
  expect_identical(out[[2]], c("a", "b"))
  expect_identical(out[[3]], c("a", "b", "c"))
})

test_that(".gwr_ms_varsets handles formulas and preserves term order", {
  ml <- list(z ~ b + a, z ~ a)
  out <- .gwr_ms_varsets(ml, "z")
  expect_identical(out[[1]], c("b", "a"))
  expect_identical(out[[2]], "a")
})

test_that(".gwr_ms_varsets strips the response from bare character vectors", {
  ml <- list(c("z", "a", "b"), c("a"))
  out <- .gwr_ms_varsets(ml, "z")
  expect_identical(out[[1]], c("a", "b"))
  expect_identical(out[[2]], "a")
})

test_that(".gwr_ms_varsets rejects an empty model list", {
  expect_error(.gwr_ms_varsets(list(), "z"), "empty model list")
  expect_error(.gwr_ms_varsets(NULL, "z"), "empty model list")
})


# --- .gwr_ms_table -----------------------------------------------------------

test_that(".gwr_ms_table ranks ascending when minimising", {
  vs  <- list("a", "b", c("a", "b"))
  tab <- .gwr_ms_table(vs, c(30, 10, 20))
  expect_identical(tab$rank, 1:3)
  expect_equal(tab$criterion, c(10, 20, 30))
  expect_identical(tab$variables, c("b", "a + b", "a"))
  expect_identical(tab$n_vars, c(1L, 2L, 1L))
  expect_identical(attr(tab, "varsets")[[1]], "b")
})

test_that(".gwr_ms_table breaks exact ties toward the smaller model", {
  vs  <- list(c("a", "b"), "a")
  tab <- .gwr_ms_table(vs, c(50, 50))
  expect_identical(tab$variables[1], "a")
  expect_identical(tab$n_vars[1], 1L)
  expect_identical(attr(tab, "varsets")[[1]], "a")
})

test_that(".gwr_ms_table sorts non-finite criteria last and never lets them win", {
  vs  <- list("a", "b", "c")
  tab <- .gwr_ms_table(vs, c(NA, 5, Inf))
  expect_identical(tab$variables[1], "b")
  expect_true(all(is.na(tab$criterion[2:3])))
})

test_that(".gwr_ms_table can maximise", {
  vs  <- list("a", "b")
  tab <- .gwr_ms_table(vs, c(0.2, 0.8), minimise = FALSE)
  expect_identical(tab$variables[1], "b")
})

test_that(".gwr_ms_table refuses to align mismatched lengths", {
  expect_error(.gwr_ms_table(list("a", "b"), c(1, 2, 3)),
               "cannot be aligned")
})


# --- .gwr_quietly ------------------------------------------------------------

test_that(".gwr_quietly swallows cat() output but returns the value", {
  noisy <- function() { cat("Now calibrating the model:\n"); 42 }
  out <- utils::capture.output(val <- .gwr_quietly(noisy(), quiet = TRUE))
  expect_identical(val, 42)
  expect_length(out, 0L)
})

test_that(".gwr_quietly lets output through when quiet = FALSE", {
  noisy <- function() { cat("Now calibrating the model:\n"); 42 }
  out <- utils::capture.output(val <- .gwr_quietly(noisy(), quiet = FALSE))
  expect_identical(val, 42)
  expect_match(paste(out, collapse = "\n"), "Now calibrating")
})

test_that(".gwr_quietly forces its expression exactly once", {
  calls <- 0L
  f <- function() { calls <<- calls + 1L; "x" }
  invisible(utils::capture.output(v <- .gwr_quietly(f(), quiet = TRUE)))
  expect_identical(calls, 1L)
  expect_identical(v, "x")
})

test_that(".gwr_quietly propagates errors and restores the sink", {
  before <- sink.number()
  expect_error(.gwr_quietly(stop("boom"), quiet = TRUE), "boom")
  expect_identical(sink.number(), before)
  expect_error(.gwr_quietly(stop("boom"), quiet = FALSE), "boom")
  expect_identical(sink.number(), before)
})

test_that(".gwr_quietly preserves a NULL return", {
  invisible(utils::capture.output(v <- .gwr_quietly(NULL, quiet = TRUE)))
  expect_null(v)
})


# --- assembly, via an injected engine ----------------------------------------

mk_sel_pts <- function(n = 40, seed = 11) {
  set.seed(seed)
  sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
               z = rnorm(n), a = rnorm(n), b = rnorm(n), cc = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
}

# Forward sweep over 3 candidates: 3 + 2 + 1 = 6 models.
stub_engine <- function(values = c(120, 118, 125, 110, 119, 112),
                        colnames_ = NULL, capture = NULL) {
  function(...) {
    args <- list(...)
    if (!is.null(capture)) assign("args", args, envir = capture)
    ml <- list(list("z", "a"), list("z", "b"), list("z", "cc"),
               list("z", c("a", "b")), list("z", c("a", "cc")),
               list("z", c("a", "b", "cc")))
    m <- cbind(seq_along(values), values)
    colnames(m) <- colnames_
    list(model_list = ml, gwr_df = m, bandwidth = 30,
         bandwidth_source = "supplied", used_dmat = FALSE,
         raw = list(model.list = ml, GWR.df = m))
  }
}

test_that("gwr_model_selection assembles a ranked result and picks the minimum", {
  pts <- mk_sel_pts()
  sel <- gwr_model_selection(pts, "z", c("a", "b", "cc"), bandwidth = 30,
                             .engine = stub_engine())
  expect_s3_class(sel, "gwr_model_selection")
  expect_identical(sel$best, c("a", "b"))          # criterion 110, the minimum
  expect_identical(nrow(sel$table), 6L)
  expect_identical(sel$table$rank, 1:6)
  expect_false(is.unsorted(sel$table$criterion))
  expect_identical(sel$table$criterion[1], 110)
  expect_identical(sel$bandwidth, 30)
  expect_identical(sel$n_models, 6L)
  expect_identical(sel$candidate_vars, c("a", "b", "cc"))
  expect_match(sel$criterion, "assumed")           # unlabelled stub table
})

test_that("gwr_model_selection reports a labelled criterion when GWmodel supplies one", {
  pts <- mk_sel_pts()
  sel <- gwr_model_selection(pts, "z", c("a", "b", "cc"), bandwidth = 30,
                             .engine = stub_engine(colnames_ = c("RSS", "AICc")))
  expect_identical(sel$criterion, "AICc")
})

test_that("gwr_model_selection hands the engine prepped, projected data", {
  env <- new.env()
  pts <- mk_sel_pts(n = 40)
  sel <- gwr_model_selection(pts, "z", c("a", "b", "cc"), bandwidth = 30,
                             adaptive = TRUE, kernel = "gaussian",
                             .engine = stub_engine(capture = env))
  a <- get("args", envir = env)
  expect_true(inherits(a$dat, "sf"))
  expect_identical(nrow(a$dat), nrow(pts))   # clean data: nothing dropped
  expect_false(sf::st_is_longlat(a$dat))
  expect_identical(a$response_var, "z")
  expect_identical(a$candidate_vars, c("a", "b", "cc"))
  expect_identical(a$kernel, "gaussian")
  expect_identical(a$bandwidth, 30)
  expect_true(a$adaptive)
  expect_true(a$quiet)
  expect_identical(sel$kernel, "gaussian")
})

test_that("gwr_model_selection validates its inputs", {
  pts <- mk_sel_pts()
  expect_error(gwr_model_selection(sf::st_drop_geometry(pts), "z", c("a", "b")),
               "must be an sf object")
  expect_error(gwr_model_selection(pts, "z", c("a", "nope"),
                                   .engine = stub_engine()),
               "absent from")
  expect_error(gwr_model_selection(pts, "z", "a", .engine = stub_engine()),
               "at least two candidate")
  # duplicates collapse, and then there is only one distinct candidate
  expect_error(gwr_model_selection(pts, "z", c("a", "a"),
                                   .engine = stub_engine()),
               "at least two candidate")
})

test_that("gwr_model_selection refuses a sweep larger than max_models", {
  pts <- mk_sel_pts()
  # 3 candidates -> 3*(3+1)/2 = 6 models
  expect_error(
    gwr_model_selection(pts, "z", c("a", "b", "cc"), bandwidth = 30,
                        max_models = 5L, .engine = stub_engine()),
    "above max_models = 5"
  )
  ok <- gwr_model_selection(pts, "z", c("a", "b", "cc"), bandwidth = 30,
                            max_models = 6L,
                            .engine = stub_engine(colnames_ = c("RSS", "AICc")))
  expect_identical(nrow(ok$table), 6L)
})

test_that("gwr_model_selection errors when no model produced a finite criterion", {
  pts <- mk_sel_pts()
  expect_error(
    gwr_model_selection(pts, "z", c("a", "b", "cc"), bandwidth = 30,
                        .engine = stub_engine(values = rep(NA_real_, 6),
                                              colnames_ = c("RSS", "AICc"))),
    "no model produced a finite criterion"
  )
})

test_that("gwr_model_selection keeps partially-failed sweeps but ranks failures last", {
  pts <- mk_sel_pts()
  sel <- gwr_model_selection(pts, "z", c("a", "b", "cc"), bandwidth = 30,
                             .engine = stub_engine(
                               values = c(NA, 118, NA, 110, NA, NA),
                               colnames_ = c("RSS", "AICc")))
  expect_identical(sel$best, c("a", "b"))
  expect_equal(sel$table$criterion[1:2], c(110, 118))
  expect_true(all(is.na(sel$table$criterion[3:6])))
})

test_that("print.gwr_model_selection shows the ranking and states the caveats", {
  pts <- mk_sel_pts()
  sel <- gwr_model_selection(pts, "z", c("a", "b", "cc"), bandwidth = 30,
                             .engine = stub_engine(colnames_ = c("RSS", "AICc")))
  out <- utils::capture.output(print(sel))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "forward model selection")
  expect_match(txt, "Selected: z ~ a \\+ b")
  expect_match(txt, "in-sample")
  expect_match(txt, "select_features_forward")
  # truncation notice when the table is longer than n
  out2 <- paste(utils::capture.output(print(sel, n = 2)), collapse = "\n")
  expect_match(out2, "4 more")
  invisible(utils::capture.output(vis <- withVisible(print(sel))))
  expect_false(vis$visible)
  expect_identical(vis$value, sel)
})


# --- GWmodel-gated integration ----------------------------------------------

test_that("gwr_model_selection reads GWmodel's criterion with the right orientation", {
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sp")

  set.seed(42)
  n <- 80
  dat <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
               a = rnorm(n), b = rnorm(n), noise = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
  dat$z <- 5 * dat$a - 3 * dat$b + rnorm(n, 0, 0.2)

  sel <- gwr_model_selection(dat, "z", c("a", "b", "noise"), bandwidth = 40)

  expect_s3_class(sel, "gwr_model_selection")
  expect_identical(nrow(sel$table), 6L)          # 3 + 2 + 1
  expect_identical(sel$table$rank, 1:6)
  expect_false(is.unsorted(sel$table$criterion, na.rm = TRUE))
  expect_identical(sel$bandwidth, 40L)
  expect_true(all(sel$best %in% c("a", "b", "noise")))

  # The load-bearing assertion.  If the criterion column were read positionally
  # and turned out to be something maximised (R-squared, say) rather than AICc,
  # or if the sort were reversed, the pure-noise predictor would outrank the
  # one carrying a 5:0.2 signal.  This catches that; a stub cannot.
  one <- sel$table[sel$table$n_vars == 1L, ]
  crit_a     <- one$criterion[one$variables == "a"]
  crit_noise <- one$criterion[one$variables == "noise"]
  expect_length(crit_a, 1L)
  expect_length(crit_noise, 1L)
  expect_lt(crit_a, crit_noise)
  expect_true("a" %in% sel$best)
})

test_that("gwr_model_selection selects a bandwidth when none is supplied", {
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sp")

  set.seed(7)
  n <- 60
  dat <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
               a = rnorm(n), b = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
  dat$z <- 3 * dat$a + rnorm(n, 0, 0.3)

  sel <- gwr_model_selection(dat, "z", c("a", "b"))
  expect_identical(nrow(sel$table), 3L)          # 2 + 1
  expect_match(sel$bandwidth_source, "bw\\.gwr")
  expect_true(is.finite(sel$bandwidth))
  expect_gte(sel$bandwidth, 4L)                  # n_cand + 2 floor
  expect_lte(sel$bandwidth, 60L)
  expect_true("a" %in% sel$best)
})

test_that("gwr_model_selection silences GWmodel's cat() progress unless asked", {
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sp")

  set.seed(9)
  n <- 50
  dat <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
               a = rnorm(n), b = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
  dat$z <- 2 * dat$a + rnorm(n, 0, 0.4)

  # GWmodel uses bare cat(), so suppressMessages()/suppressWarnings() do
  # nothing; only a sink can catch it.  Any residual lines here are ours
  # (.log_warn), which is why this matches on GWmodel's own strings rather
  # than asserting zero output.
  noise <- "Now calibrating|Adaptive bandwidth"

  quiet_out <- utils::capture.output(
    s1 <- gwr_model_selection(dat, "z", c("a", "b"), bandwidth = 25,
                              quiet = TRUE))
  expect_false(any(grepl(noise, quiet_out)))

  loud_out <- utils::capture.output(
    s2 <- gwr_model_selection(dat, "z", c("a", "b"), bandwidth = 25,
                              quiet = FALSE))
  expect_true(any(grepl("Now calibrating", loud_out)))

  # Silencing must not change the answer.
  expect_identical(s1$best, s2$best)
  expect_equal(s1$table$criterion, s2$table$criterion)
})
