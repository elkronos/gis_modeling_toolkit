# tests/testthat/test-fix-gwr-selection-null.R
# ---------------------------------------------------------------------------
# gwr_model_selection() scores the model without predictors and reports the
# criterion gaps.
#
# The forward sweep starts from one predictor, so by itself it always names
# at least one, whether or not any candidate helps: on a response unrelated
# to three candidates it selected one of them every time.  The table gave raw
# AICc and left the subtraction to the reader.  The sweep is also the
# package's own, on the local solve, so a candidate that is constant inside
# some windows is scored on what they estimate.
# ---------------------------------------------------------------------------

.sn_catch <- function(expr) {
  w <- character(0)
  val <- tryCatch(
    withCallingHandlers(suppressMessages(expr), warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    }),
    error = function(e) conditionMessage(e))
  list(value = val, warnings = w)
}

# The example of ?gwr_model_selection.
.sn_example <- function() {
  set.seed(1); n <- 80
  dat <- sf::st_as_sf(
    data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
               a = rnorm(n), b = rnorm(n), noise = rnorm(n)),
    coords = c("x", "y"), crs = 32632)
  dat$z <- 2 * dat$a - dat$b + rnorm(n, 0, 0.5)
  dat
}


test_that("the table carries the model without predictors and each model's gap to the best", {
  skip_if_not_installed("sp")
  dat <- .sn_example()
  sel <- suppressMessages(gwr_model_selection(dat, "z", c("a", "b", "noise"), bandwidth = 30))
  # The documented example still selects a + b.
  expect_identical(sel$best, c("a", "b"))
  tab <- sel$table
  expect_named(tab, c("rank", "n_vars", "variables", "criterion", "delta"))
  expect_identical(nrow(tab), 7L)                 # 3 + 2 + 1 and <none>
  expect_identical(sel$n_models, 7L)
  none <- which(tab$variables == "<none>")
  expect_length(none, 1L)
  expect_identical(tab$n_vars[none], 0L)
  expect_true(is.finite(tab$criterion[none]))
  # delta is the criterion minus the best: 0 in the first row.
  expect_identical(tab$delta[1L], 0)
  expect_equal(tab$delta, tab$criterion - min(tab$criterion))
  expect_false(is.unsorted(tab$delta))
  # The local mean is far behind a model that explains the response.
  expect_gt(tab$delta[none], 50)
  # It is a GWR without predictors at the sweep's bandwidth: the AICc of
  # local weighted means of the response.
  xy <- sf::st_coordinates(dat)
  cr <- .gwr_criterion(matrix(1, 80, 1), dat$z, xy, 30L, "bisquare", TRUE)
  expect_equal(tab$criterion[none], cr$value)
  w1 <- .gw_kernel_weights(sqrt((xy[, 1] - xy[1, 1])^2 + (xy[, 2] - xy[1, 2])^2),
                           30L, "bisquare", TRUE)
  sv <- .gwr_local_solve(matrix(1, 80, 1, dimnames = list(NULL, "Intercept")), dat$z,
                         xy, 30L, "bisquare", TRUE, what = "score")
  expect_equal(sv$fitted[1L], sum(w1 * dat$z) / sum(w1))
  expect_equal(sv$leverage[1L], 1 / sum(w1))
  # The models within 2 of the best, and the smallest of them.
  expect_identical(sel$within_2, tab[tab$delta <= 2, , drop = FALSE])
  expect_identical(sel$smallest_within_2, c("a", "b"))
  # raw is the forward sweep alone, in GWmodel's layout.
  expect_length(sel$raw[[1L]], 6L)
  expect_identical(dim(sel$raw[[2L]]), c(6L, 4L))
  expect_identical(sel$criterion, "AICc")
  expect_true(sel$criterion_by_name)
  expect_true(sel$criterion_verified)
})

test_that("on a response unrelated to the candidates the model without predictors is usually selected", {
  skip_if_not_installed("sp")
  # Before: one of the three candidates was selected in 8 of 8 draws, because
  # nothing else could be.
  best <- vapply(1:8, function(seed) {
    set.seed(100 + seed); m <- 150
    d <- sf::st_as_sf(data.frame(x = runif(m, 0, 1000), y = runif(m, 0, 1000),
                                 p1 = rnorm(m), p2 = rnorm(m), p3 = rnorm(m), z = rnorm(m)),
                      coords = c("x", "y"), crs = 32632)
    s <- suppressMessages(gwr_model_selection(d, "z", c("p1", "p2", "p3"), bandwidth = 60))
    expect_identical(s$table$delta[1L], 0)
    # Whatever is selected, the model without predictors is in the running.
    expect_true("<none>" %in% s$table$variables)
    if (length(s$best)) paste(s$best, collapse = "+") else "<none>"
  }, character(1))
  expect_gte(sum(best == "<none>"), 5L)
})

test_that("print() shows the gaps, the models within 2 of the best and the smallest of them", {
  # Built on a stub engine, so it runs without GWmodel: six models of a
  # sweep, and a model without predictors that the engine scored.
  set.seed(11)
  pts <- sf::st_as_sf(data.frame(x = runif(40, 0, 1000), y = runif(40, 0, 1000),
                                 z = rnorm(40), a = rnorm(40), b = rnorm(40), cc = rnorm(40)),
                      coords = c("x", "y"), crs = 32632)
  eng <- function(aicc, null = NULL) function(...) {
    ml <- list(list("z", "a"), list("z", "b"), list("z", "cc"), list("z", c("a", "b")),
               list("z", c("a", "cc")), list("z", c("a", "b", "cc")))
    m <- cbind(bandwidth = 30, AIC = aicc - 5, AICc = aicc, RSS = 1)
    out <- list(model_list = ml, gwr_df = m, bandwidth = 30, bandwidth_source = "supplied",
                used_dmat = FALSE, raw = list(ml, m))
    if (!is.null(null)) out$null_model <- c(bandwidth = 30, AIC = null - 5, AICc = null, RSS = 1)
    out
  }
  # a + b best; a + b + cc and a within 2 of it; a is the smallest of them.
  sel <- gwr_model_selection(pts, "z", c("a", "b", "cc"), bandwidth = 30,
                             .engine = eng(c(111.5, 118, 125, 110, 119, 110.8), null = 140))
  expect_identical(sel$best, c("a", "b"))
  expect_identical(sel$table$variables[1:3], c("a + b", "a + b + cc", "a"))
  expect_equal(sel$table$delta[1:3], c(0, 0.8, 1.5))
  expect_identical(sel$within_2$variables, c("a + b", "a + b + cc", "a"))
  expect_identical(sel$smallest_within_2, "a")
  expect_identical(sel$table$variables[7L], "<none>")
  txt <- paste(utils::capture.output(print(sel)), collapse = "\n")
  expect_match(txt, "delta")
  expect_match(txt, "Selected: z ~ a \\+ b")
  expect_match(txt, "Within 2 of the best \\(a difference this small does not separate models\\):")
  expect_match(txt, "z ~ a \\+ b  \\(\\+0\\.00\\)")
  expect_match(txt, "z ~ a \\+ b \\+ cc  \\(\\+0\\.80\\)")
  expect_match(txt, "z ~ a  \\(\\+1\\.50\\)")
  expect_match(txt, "Smallest of them: z ~ a\n")
  # A clear winner: the runner-up's gap is stated instead.
  clear <- gwr_model_selection(pts, "z", c("a", "b", "cc"), bandwidth = 30,
                               .engine = eng(c(120, 118, 125, 110, 119, 115), null = 140))
  txt2 <- paste(utils::capture.output(print(clear)), collapse = "\n")
  expect_match(txt2, "No other model is within 2 of the best; the runner-up is 5\\.00 away\\.")
  expect_identical(clear$smallest_within_2, c("a", "b"))
  # The model without predictors wins: best is empty, and print says so.
  nul <- gwr_model_selection(pts, "z", c("a", "b", "cc"), bandwidth = 30,
                             .engine = eng(c(120, 118, 125, 110, 119, 115), null = 109))
  expect_identical(nul$best, character(0))
  expect_identical(nul$table$variables[1L], "<none>")
  expect_identical(nul$smallest_within_2, character(0))
  txt3 <- paste(utils::capture.output(print(nul)), collapse = "\n")
  expect_match(txt3, "Selected: z ~ 1\n")
  expect_match(txt3, "The model without predictors has the lowest criterion: no candidate")
  expect_match(txt3, "z ~ 1  \\(\\+0\\.00\\)")
  expect_match(txt3, "z ~ a \\+ b  \\(\\+1\\.00\\)")
  expect_match(txt3, "Smallest of them: z ~ 1\n")
  # An engine that scores no such model gives the table it always gave, with
  # the gaps.
  old <- gwr_model_selection(pts, "z", c("a", "b", "cc"), bandwidth = 30,
                             .engine = eng(c(120, 118, 125, 110, 119, 115)))
  expect_identical(nrow(old$table), 6L)
  expect_false("<none>" %in% old$table$variables)
  expect_equal(old$table$delta, old$table$criterion - 110)
  # The table helper alone: delta, also when maximising and when nothing is
  # finite.
  expect_identical(.gwr_ms_table(list("a", character(0)), c(12, 10))$delta, c(0, 2))
  expect_identical(.gwr_ms_table(list("a", character(0)), c(12, 10))$variables, c("<none>", "a"))
  expect_equal(.gwr_ms_table(list("a", "b"), c(0.2, 0.8), minimise = FALSE)$delta, c(0, 0.6))
  expect_true(all(is.na(.gwr_ms_table(list("a", "b"), c(NA, NA))$delta)))
  expect_identical(.gwr_ms_table(list("a", "b", "c"), c(NA, 5, 7))$delta, c(0, 2, NA))
})

test_that("a well-posed sweep is the one GWmodel::gwr.model.selection() runs", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  dat <- .sn_example()
  for (cs in list(list(bw = 30L, adaptive = TRUE, kernel = "bisquare"),
                  list(bw = 300, adaptive = FALSE, kernel = "gaussian"))) {
    sel <- suppressMessages(gwr_model_selection(dat, "z", c("a", "b", "noise"), bandwidth = cs$bw,
                                                adaptive = cs$adaptive, kernel = cs$kernel))
    sp_d <- .to_sp(dat, c("z", "a", "b", "noise"))
    direct <- .gwr_quietly(GWmodel::gwr.model.selection(
      DeVar = "z", InDeVars = c("a", "b", "noise"), data = sp_d, bw = cs$bw,
      adaptive = cs$adaptive, kernel = cs$kernel,
      dMat = GWmodel::gw.dist(dp.locat = sp::coordinates(sp_d))))
    # The same models in the same order, and the same bandwidth, AIC, AICc
    # and RSS for each.
    expect_identical(sel$raw[[1L]], direct[[1L]], info = cs$kernel)
    expect_equal(unname(sel$raw[[2L]]), direct[[2L]], tolerance = 1e-9, info = cs$kernel)
    expect_identical(colnames(sel$raw[[2L]]), c("bandwidth", "AIC", "AICc", "RSS"))
  }
  # quiet = FALSE writes one line per model, and changes nothing else.
  out <- utils::capture.output(
    loud <- gwr_model_selection(dat, "z", c("a", "b", "noise"), bandwidth = 30, quiet = FALSE))
  expect_identical(sum(grepl("^Now calibrating model [1-6] of 6: z~", out)), 6L)
  expect_identical(sum(grepl("^Now calibrating the model without predictors: z~1$", out)), 1L)
  quiet <- suppressMessages(gwr_model_selection(dat, "z", c("a", "b", "noise"), bandwidth = 30))
  expect_identical(loud$table, quiet$table)
  # dmat_max_n is accepted and ignored: no distance matrix is built, so a
  # fixed bandwidth needs none.
  f0 <- suppressMessages(gwr_model_selection(dat, "z", c("a", "b"), bandwidth = 300,
                                             adaptive = FALSE, dmat_max_n = 0))
  f1 <- suppressMessages(gwr_model_selection(dat, "z", c("a", "b"), bandwidth = 300,
                                             adaptive = FALSE))
  expect_identical(f0$table, f1$table)
  expect_false(f0$used_dmat)
})

test_that("a candidate constant inside some windows is ranked on what those windows estimate", {
  skip_if_not_installed("sp")
  # z = 1 + 2 a + 3 (east of x = 500) + noise; `region` codes east and west.
  # Before: coded 0/1 the sweep stopped ("inv(): matrix is singular"); coded
  # 0.7/1.3 it returned with no warning and AICc of 1810 to 1929 for the
  # three models that hold `region`, the true driver, against 534 for `a`
  # alone.
  set.seed(42); n <- 300
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000); a <- rnorm(n); east <- x > 500
  eps <- rnorm(n, 0, 0.3); b <- rnorm(n)
  tabs <- lapply(list(c(0, 1), c(0.7, 1.3)), function(lev) {
    d <- sf::st_as_sf(data.frame(x = x, y = y, a = a, b = b, z = 1 + 2 * a + 3 * east + eps,
                                 region = ifelse(east, lev[2], lev[1])),
                      coords = c("x", "y"), crs = 32632)
    r <- .sn_catch(gwr_model_selection(d, "z", c("a", "b", "region"), bandwidth = 40))
    sel <- r$value
    expect_s3_class(sel, "gwr_model_selection")
    expect_length(r$warnings, 0L)
    expect_identical(sel$best, c("a", "region"))
    expect_identical(sel$table$variables[1:2], c("a + region", "a + region + b"))
    expect_true(all(is.finite(sel$table$criterion)))
    sel$table
  })
  # One data set in two codings: one table.
  expect_identical(tabs[[2L]]$variables, tabs[[1L]]$variables)
  expect_equal(tabs[[2L]]$criterion, tabs[[1L]]$criterion, tolerance = 1e-8)
  # The AICc fit_gwr_model() gives the same model on these rows (pinned in
  # test-fix-gwr-local-solve.R, where it was recorded from the solve).
  expect_equal(tabs[[1L]]$criterion[1L], 145.2796, tolerance = 1e-5)
})
