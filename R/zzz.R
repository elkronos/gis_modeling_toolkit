#' @noRd
.onLoad <- function(libname, pkgname) {
  # Set up default logging in a package-specific namespace so we never
  # overwrite the user's global logger configuration.
  # Users can reconfigure the spatialkit namespace after loading -- but note
  # that logger::log_appender() and log_threshold() BOTH default to index = 1,
  # so the two-line recipe below without an index touches only the temp-file
  # appender and leaves the console echo exactly as it was.  To redirect or
  # quieten what is printed, name index 2:
  #   logger::log_appender(logger::appender_file("my_log.log"),
  #                        namespace = "spatialkit", index = 2)
  #   logger::log_threshold(logger::FATAL, namespace = "spatialkit", index = 2)
  # spatialkit_quiet() does the second of those for you.

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


#' Quieten (or restore) spatialkit's console log
#'
#' The package writes an INFO+ trace to a session temp file (logger appender
#' index 1) and echoes WARN+ to the console (index 2).  Both
#' \code{logger::log_appender()} and \code{logger::log_threshold()} default to
#' \code{index = 1}, so the obvious two-line recipe silences the file and
#' leaves the console untouched -- which is the opposite of what anyone wants.
#' This helper names the right index.
#'
#' Note that these are log records, not R conditions: \code{suppressWarnings()}
#' and \code{tryCatch(warning = )} do not see them.  Conditions the package
#' raises as real R warnings are unaffected by this function.
#'
#' @param quiet Logical.  \code{TRUE} (default) silences the console echo;
#'   \code{FALSE} restores the WARN+ default.
#' @return Invisibly, the threshold that was in force before the change.
#' @family utilities
#' @examples
#' old <- spatialkit_quiet()      # console echo off
#' spatialkit_quiet(FALSE)        # back to WARN+
#' @export
spatialkit_quiet <- function(quiet = TRUE) {
  if (!is.logical(quiet) || length(quiet) != 1L || is.na(quiet))
    stop("spatialkit_quiet(): `quiet` must be TRUE or FALSE.", call. = FALSE)
  prev <- tryCatch(logger::log_threshold(namespace = "spatialkit", index = 2),
                   error = function(e) NULL)
  logger::log_threshold(if (quiet) logger::FATAL else logger::WARN,
                        namespace = "spatialkit", index = 2)
  invisible(prev)
}
