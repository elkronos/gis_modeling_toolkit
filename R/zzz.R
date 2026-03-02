#' @noRd
.onLoad <- function(libname, pkgname) {
  # Set up default logging in a package-specific namespace so we never
  # overwrite the user's global logger configuration.
  # Users can reconfigure the spatialkit namespace after loading with:
  #   logger::log_appender(logger::appender_file("my_log.log"), namespace = "spatialkit")
  #   logger::log_threshold(logger::WARN, namespace = "spatialkit")

  log_path <- file.path(tempdir(), "spatialkit_model_log.log")
  logger::log_appender(logger::appender_file(log_path), namespace = "spatialkit")
  logger::log_threshold(logger::INFO, namespace = "spatialkit")
}
