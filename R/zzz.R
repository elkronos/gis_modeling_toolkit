#' @noRd
.onLoad <- function(libname, pkgname) {
  # Set up default logging in a package-specific namespace so we never
  # overwrite the user's global logger configuration.
  # Users can reconfigure the spatialkit namespace after loading with:
  #   logger::log_appender(logger::appender_file("my_log.log"), namespace = "spatialkit")
  #   logger::log_threshold(logger::WARN, namespace = "spatialkit")

  # Index 1: full INFO+ trace to a session temp file (detailed diagnostics).
  log_path <- file.path(tempdir(), "spatialkit_model_log.log")
  logger::log_appender(logger::appender_file(log_path),
                       namespace = "spatialkit", index = 1)
  logger::log_threshold(logger::INFO, namespace = "spatialkit", index = 1)

  # Index 2: WARN+ to the console so that important problems (skipped CV
  # folds, failed predictions returning NA, extraction failures, ...) are
  # actually visible to interactive users instead of only landing in a
  # temp file nobody reads.
  logger::log_appender(logger::appender_console,
                       namespace = "spatialkit", index = 2)
  logger::log_threshold(logger::WARN, namespace = "spatialkit", index = 2)
}
