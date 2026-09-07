# ===========================================================================
# Worker-count resolution for the parallel paths.
#
# parallel::detectCores() returns NA on some platforms, and NA arithmetic
# propagates rather than erroring, so an unsanitised count reaches
# mclapply() as NA. See test-cv-parallel.R for fold-level reproducibility.
# ===========================================================================

test_that(".sanitize_core_count collapses NA/invalid input to the fallback", {
  f <- spatialkit:::.sanitize_core_count

  expect_identical(f(NA), 1L)
  expect_identical(f(NA_integer_ - 1L), 1L)   # detectCores() NA arithmetic
  expect_identical(f(0), 1L)
  expect_identical(f(-3), 1L)
  expect_identical(f(NULL), 1L)
  expect_identical(f("not a number"), 1L)
  expect_identical(f(7), 7L)
  expect_identical(f(3.9), 3L)
  expect_identical(f(NA, fallback = 2L), 2L)
})


test_that(".resolve_n_cores returns sane values for the sequential paths", {
  f <- spatialkit:::.resolve_n_cores

  expect_identical(f(parallel = FALSE), 1L)
  expect_identical(f(parallel = FALSE, n_cores = NULL), 1L)
  expect_identical(f(parallel = 1), 1L)
  expect_identical(f(parallel = NA), 1L)
  expect_identical(f(parallel = NULL), 1L)

  skip_on_os("windows")
  # An explicit request is honoured up to the machine's core count, and up to
  # two under R CMD check (_R_CHECK_LIMIT_CORES_), which is what CRAN's check
  # farm enforces; more workers than cores is not parallelism.
  cap <- function(x) {
    x <- min(x, parallel::detectCores(logical = TRUE))
    if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) x <- min(x, 2L)
    as.integer(x)
  }
  expect_identical(suppressMessages(f(parallel = 3)), cap(3L))
  expect_identical(suppressMessages(f(parallel = 2L)), cap(2L))
  expect_identical(suppressMessages(f(parallel = FALSE, n_cores = 2)), cap(2L))
  # n_cores wins over parallel, and is sanitised on the way through.
  expect_identical(suppressMessages(f(parallel = 8L, n_cores = 3L)), cap(3L))
  expect_identical(f(parallel = TRUE, n_cores = NA), 1L)
  # The machine cap says so, and never returns more than the machine has.
  n_machine <- parallel::detectCores(logical = TRUE)
  if (is.finite(n_machine) && !nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) {
    expect_message(got <- f(parallel = n_machine + 5L),
                   "workers requested on a machine with")
    expect_identical(got, as.integer(n_machine))
  }
  # Under R CMD check the cap is two whatever was asked for.
  with_check_limit <- function(expr) {
    old <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA)
    Sys.setenv("_R_CHECK_LIMIT_CORES_" = "TRUE")
    on.exit(if (is.na(old)) Sys.unsetenv("_R_CHECK_LIMIT_CORES_") else
              Sys.setenv("_R_CHECK_LIMIT_CORES_" = old))
    force(expr)
  }
  expect_lte(with_check_limit(suppressMessages(f(parallel = 8L))), 2L)
  expect_identical(with_check_limit(suppressMessages(f(parallel = 2L))),
                   min(2L, as.integer(n_machine)))
})


test_that(".resolve_n_cores auto-detection returns a usable worker count", {
  # parallel = TRUE reads detectCores(), which is machine-dependent, so the
  # assertable contract is the type and the floor -- an NA reaching mclapply()
  # is the failure this guards.
  skip_on_os("windows")
  cores <- spatialkit:::.resolve_n_cores(TRUE)
  expect_type(cores, "integer")
  expect_length(cores, 1L)
  expect_false(is.na(cores))
  expect_gte(cores, 1L)
  # Never more than the machine has; detectCores() - 1 is the documented rule.
  detected <- parallel::detectCores(logical = FALSE)
  if (!is.na(detected) && detected > 1L) expect_lte(cores, detected - 1L)
})
