# tests/testthat/test-review2-logging.R
# ---------------------------------------------------------------------------
# The package's own logger setup (R/zzz.R) and the helpers that write to it
# (R/utils.R): a user's logger configuration cannot break or tap it, a deleted
# session temp directory cannot turn a log line into an error, and a knitted
# document shows the cautions that are only logged.
# ---------------------------------------------------------------------------

# logger keeps every namespace's configuration in one internal environment.
# Tests that rebuild the "spatialkit" namespace, or change the user's global
# one, save both entries first and put them back afterwards.
.logger_namespaces <- function() get("namespaces", envir = asNamespace("logger"))

local_logger_state <- function(env = parent.frame()) {
  ns <- .logger_namespaces()
  saved_global <- get("global", envir = ns)
  saved_sk     <- get("spatialkit", envir = ns)
  withr::defer({
    assign("global", saved_global, envir = ns)
    assign("spatialkit", saved_sk, envir = ns)
  }, envir = env)
  invisible(NULL)
}

read_or_empty <- function(f) if (file.exists(f)) readLines(f) else character(0)


test_that("a global logger configuration made before loading neither breaks nor taps spatialkit's log", {
  local_logger_state()
  ns <- .logger_namespaces()
  f2 <- withr::local_tempfile()
  f3 <- withr::local_tempfile()
  # A user who configured logging BEFORE library(spatialkit): sprintf
  # formatting, a second index writing to a file, a third one to another.
  logger::log_formatter(logger::formatter_sprintf)
  logger::log_appender(logger::appender_file(f2), index = 2)
  logger::log_formatter(logger::formatter_sprintf, index = 2)
  logger::log_appender(logger::appender_file(f3), index = 3)
  logger::log_formatter(logger::formatter_glue, index = 3)
  global_before <- get("global", envir = ns)

  # As if the package were loading for the first time: logger then seeds the
  # namespace by copying all three of the user's indices.
  rm("spatialkit", envir = ns)
  spatialkit:::.onLoad(NULL, "spatialkit")
  expect_identical(get("global", envir = ns), global_before)

  # The console echo inherited formatter_sprintf, so a `%` in a message
  # aborted the caller ("too few arguments") and the promised R warning never
  # arrived; formatter_glue on index 3 did the same for a `{`.
  lines <- capture_spatialkit_log({
    expect_warning(
      spatialkit:::.warn_and_log("fold 2 skipped: %s",
                                 "object 'cov_{x' not found; 14% done"),
      "fold 2 skipped: object 'cov_{x' not found; 14% done", fixed = TRUE)
    expect_no_error(spatialkit:::.log_warn("distortion %.1f%% over {the} extent", 14))
    expect_no_error(spatialkit:::.log_info("an INFO line with %s", "50% {x}"))
  })
  expect_true(log_has(lines, "cov_\\{x' not found; 14% done"))
  expect_true(log_has(lines, "distortion 14.0% over \\{the\\} extent"))

  # And the user's own log files receive nothing from spatialkit.
  expect_identical(read_or_empty(f2), character(0))
  expect_identical(read_or_empty(f3), character(0))

  # Re-running the setup on the namespace it built is harmless.
  expect_no_error(spatialkit:::.onLoad(NULL, "spatialkit"))
  capture_spatialkit_log(
    expect_warning(spatialkit:::.warn_and_log("again %s", "5% {y}"),
                   "again 5% {y}", fixed = TRUE))
  expect_identical(read_or_empty(f3), character(0))
})


test_that("a log line that cannot be written never costs the caller its warning", {
  # A session temp directory deleted under the file trace made every call that
  # logs fail with "cannot open the connection"; .warn_and_log() logs before
  # it warns, so the R warning the manual promises became that error.
  withr::defer(logger::log_appender(spatialkit:::.sk_file_appender(),
                                    namespace = "spatialkit", index = 1))

  # The package's own trace pointed somewhere that cannot be written: the
  # line is dropped there and still reaches the console echo.
  unwritable <- file.path(withr::local_tempfile(lines = "a file"), "x", "t.log")
  logger::log_appender(spatialkit:::.sk_file_appender(function() unwritable),
                       namespace = "spatialkit", index = 1)
  lines <- capture_spatialkit_log({
    expect_warning(spatialkit:::.warn_and_log("%s(): dropped %d row(s)", "f", 3L),
                   "f(): dropped 3 row(s)", fixed = TRUE)
    expect_no_error(spatialkit:::.log_warn("still %s", "logging"))
    expect_no_error(spatialkit:::.log_info("and %s", "informing"))
  })
  expect_true(log_has(lines, "f\\(\\): dropped 3 row\\(s\\)"))
  expect_true(log_has(lines, "still logging"))

  # An appender that throws outright -- one a user installed, say -- costs
  # the log line, never the computation or the warning.
  logger::log_appender(function(lines) stop("cannot open the connection"),
                       namespace = "spatialkit", index = 1)
  expect_warning(spatialkit:::.warn_and_log("%s(): dropped %d row(s)", "g", 4L),
                 "g(): dropped 4 row(s)", fixed = TRUE)
  expect_no_error(spatialkit:::.log_warn("still %s", "running"))
  expect_no_error(spatialkit:::.log_info("and %s", "informing"))
})


test_that("the temp-file trace recreates its deleted directory and never throws", {
  dir <- withr::local_tempdir()
  f <- file.path(dir, "session", "trace.log")
  app <- spatialkit:::.sk_file_appender(function() f)
  app("first")
  expect_identical(readLines(f), "first")
  unlink(file.path(dir, "session"), recursive = TRUE)
  expect_no_error(app("second"))
  expect_identical(readLines(f), "second")
  # A path that cannot be written at all loses the line, quietly.
  blocked <- file.path(f, "not-a-directory", "trace.log")
  expect_silent(spatialkit:::.sk_file_appender(function() blocked)("third"))
  expect_false(file.exists(blocked))

  # The package's own trace follows tempdir() when a line is written rather
  # than holding the path it saw at load time, and logger's getter still
  # reports it as the call that generated it.
  expect_true(is.call(logger::log_appender(namespace = "spatialkit", index = 1)))
})


test_that("capturing a backend's stderr falls back when there is nowhere to divert it", {
  skip_if(!identical(as.integer(sink.number(type = "message")), 2L),
          "the message stream is already diverted")
  gone <- file.path(withr::local_tempdir(), "deleted", "stderr.txt")
  res <- spatialkit:::.call_capturing_stderr(function() 42, path = gone)
  expect_identical(res$value, 42)
  expect_null(res$error)
  res <- spatialkit:::.call_capturing_stderr(function() stop("boom"), path = gone)
  expect_match(conditionMessage(res$error), "boom")
  expect_identical(as.integer(sink.number(type = "message")), 2L)
})


test_that("while a document is knitted, a logged caution also reaches it as a message", {
  # knitr does not capture stderr, which is where the console echo writes, so
  # a knitted report showed the package's R warnings but none of its logged
  # cautions.
  withr::local_options(knitr.in.progress = TRUE)
  app <- spatialkit:::.sk_console_appender
  # The stderr copy is unchanged; the message is additional.
  err <- utils::capture.output(
    expect_message(app("WARN [t] a caution"), "WARN [t] a caution", fixed = TRUE),
    type = "message")
  expect_identical(err, "WARN [t] a caution")

  # Through the package's own helper and default configuration.
  utils::capture.output(
    expect_message(spatialkit:::.log_warn("%s(): 5%% of {cells} empty", "f"),
                   "f(): 5% of {cells} empty", fixed = TRUE),
    type = "message")

  # A line that is also raised as an R warning is not repeated as a message:
  # the document shows the warning already.
  utils::capture.output(
    expect_no_message(expect_warning(spatialkit:::.warn_and_log("raised once"),
                                     "raised once")),
    type = "message")
})


test_that("outside knitr the console echo is not a message", {
  withr::local_options(knitr.in.progress = NULL)
  err <- utils::capture.output(
    expect_no_message(spatialkit:::.log_warn("plain %s", "console")),
    type = "message")
  expect_true(any(grepl("plain console", err, fixed = TRUE)))
})
