# tests/testthat/test-review-cv.R
# ---------------------------------------------------------------------------
# Regressions from the review of the CV runners, forward selection and
# compare_models_cv(): a parallel cv_gwr() that hung, partial fold failures
# that were only logged, and scores compared across different row sets.
# ---------------------------------------------------------------------------

.rv_gwr_pts <- function(n = 90, seed = 21) {
  set.seed(seed)
  d <- sf::st_as_sf(
    data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
               elev = rnorm(n)),
    coords = c("x", "y"), crs = 32632)
  d$price <- 10 + 2 * d$elev + rnorm(n, 0, 0.5)
  d
}


# ---------------------------------------------------------------------------
# cv_gwr(parallel = n) after a GWR fit
# ---------------------------------------------------------------------------

test_that("cv_gwr(parallel = 2) returns after a GWR has been fitted in the session", {
  # GWmodel is built with OpenMP.  Once fit_gwr_model() had run, GNU
  # libgomp's thread pool lived in the parent, every mclapply() child blocked
  # on a futex at its first OpenMP region, and cv_gwr(parallel = 2) never
  # returned.  Run in a child R process under `timeout` so a regression fails
  # here (status 124) instead of hanging the whole suite.
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sp")
  skip_if(parallel::detectCores() < 2L, "fewer than two cores: nothing forks")
  skip_if(!nzchar(Sys.which("timeout")), "coreutils `timeout` not found")

  # Load the same spatialkit this session has: the installed copy, or the
  # source tree pkgload::load_all() put in place.
  ns_path   <- getNamespaceInfo(asNamespace("spatialkit"), "path")
  installed <- file.exists(file.path(ns_path, "Meta", "package.rds"))
  if (!installed) skip_if_not_installed("pkgload")
  load_line <- if (installed)
    sprintf("suppressMessages(library(spatialkit, lib.loc = %s))",
            deparse(dirname(ns_path)))
  else
    sprintf("suppressMessages(pkgload::load_all(%s, quiet = TRUE))",
            deparse(ns_path))

  script <- tempfile(fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    load_line,
    "logger::log_threshold(logger::FATAL, namespace = 'spatialkit', index = 2)",
    "set.seed(21); n <- 90",
    "d <- sf::st_as_sf(",
    "  data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),",
    "             elev = rnorm(n)),",
    "  coords = c('x', 'y'), crs = 32632)",
    "d$price <- 10 + 2 * d$elev + rnorm(n, 0, 0.5)",
    "invisible(suppressWarnings(suppressMessages(fit_gwr_model(d, 'price', 'elev'))))",
    "seq_res <- suppressWarnings(suppressMessages(",
    "  cv_gwr(d, 'price', 'elev', k = 3, seed = 7)))",
    "w <- character(0)",
    "par_res <- withCallingHandlers(suppressMessages(",
    "  cv_gwr(d, 'price', 'elev', k = 3, seed = 7, parallel = 2L)),",
    "  warning = function(x) {",
    "    w <<- c(w, conditionMessage(x)); invokeRestart('muffleWarning') })",
    "cat('RESULT', par_res$n_folds_succeeded,",
    "    isTRUE(all.equal(par_res$overall, seq_res$overall)),",
    "    any(grepl('OpenMP', w)), '\\n')"
  ), script)

  rscript <- file.path(R.home("bin"), "Rscript")
  # `timeout` makes itself a process-group leader and signals the group, so
  # any forked worker stuck in libgomp dies with the R process.  R_TESTS is
  # cleared because R CMD check points it at a startup file relative to the
  # tests directory, which the child would fail to source.
  out <- suppressWarnings(system2(
    "timeout", c("-k", "10", "150", shQuote(rscript), shQuote(script)),
    stdout = TRUE, stderr = TRUE, env = "R_TESTS="))
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  expect_identical(as.integer(status), 0L,
                   label = paste0("child R exit status (124 = hung and killed); output:\n",
                                  paste(utils::tail(out, 15L), collapse = "\n")))
  # All three folds scored, identical to the sequential run, and a warning
  # that says why the folds did not fork.
  expect_true(any(grepl("^RESULT 3 TRUE TRUE", out)),
              label = paste(utils::tail(out, 15L), collapse = "\n"))
})


test_that("compare_models_cv(gwr_args = list(parallel = 2)) never forks GWR folds", {
  # The same hang reached through compare_models_cv().  mclapply() is
  # replaced by a stub that errors, so the old code fails here at once
  # instead of hanging on the GWR fits earlier tests have left in the session.
  skip_on_os("windows")
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sp")
  skip_if(parallel::detectCores() < 2L, "fewer than two cores: nothing forks")

  d <- .rv_gwr_pts()
  invisible(suppressWarnings(suppressMessages(fit_gwr_model(d, "price", "elev"))))
  ref <- suppressWarnings(suppressMessages(
    compare_models_cv(d, "price", "elev", models = "GWR", k = 3, seed = 7,
                      quiet = TRUE)))
  local_mocked_bindings(
    mclapply = function(...) stop("GWR folds were forked"),
    .package = "parallel")
  expect_warning(
    res <- suppressMessages(
      compare_models_cv(d, "price", "elev", models = "GWR", k = 3, seed = 7,
                        quiet = TRUE, gwr_args = list(parallel = 2L))),
    "cv_gwr\\(\\): `parallel` is ignored.*OpenMP")
  expect_identical(res$gwr_cv$n_folds_succeeded, 3L)
  expect_equal(res$overall, ref$overall)
})
