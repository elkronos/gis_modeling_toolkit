# tests/testthat/helper-logging.R
# ---------------------------------------------------------------------------
# Capture spatialkit's log output.
#
# The package logs through logger::log_warn() / log_info() into the
# "spatialkit" namespace (see R/zzz.R): index 1 appends INFO+ to a temp file,
# index 2 sends WARN+ to appender_console.  Neither raises an R condition, so
# expect_warning() and expect_message() do NOT see these lines -- an assertion
# written that way passes vacuously whether or not the message is emitted.
#
# capture_spatialkit_log() temporarily installs a recording appender and
# returns the emitted lines, so tests can assert on them for real.
# ---------------------------------------------------------------------------

capture_spatialkit_log <- function(expr, level = logger::INFO) {
  captured <- character(0)

  # Register the restore FIRST.  If either mutation below (or `expr` itself)
  # throws, the recording appender must not be left installed -- it would
  # swallow every subsequent test's log output for the rest of the session and
  # silently disable appender_console.
  #
  # Restores the configuration .onLoad() installs (R/zzz.R: appender_console at
  # WARN on index 2) rather than round-tripping through the getter.
  on.exit({
    logger::log_appender(logger::appender_console,
                         namespace = "spatialkit", index = 2)
    logger::log_threshold(logger::WARN, namespace = "spatialkit", index = 2)
  }, add = TRUE)

  logger::log_appender(function(lines) captured <<- c(captured, lines),
                       namespace = "spatialkit", index = 2)
  logger::log_threshold(level, namespace = "spatialkit", index = 2)

  force(expr)
  captured
}

# TRUE when any captured line matches `pattern`.
log_has <- function(lines, pattern) {
  any(grepl(pattern, lines, fixed = FALSE))
}
