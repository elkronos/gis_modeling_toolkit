#' @noRd
.onLoad <- function(libname, pkgname) {
  # Set up default logging in a package-specific namespace so we never
  # overwrite the user's global logger configuration -- and, the other way
  # round, so the user's global configuration cannot break ours (see the
  # formatter note below).
  # Users can reconfigure the spatialkit namespace after loading -- but note
  # that logger::log_appender() and log_threshold() BOTH default to index = 1,
  # so the two-line recipe below without an index touches only the temp-file
  # appender and leaves the console echo exactly as it was.  To redirect or
  # quieten what is printed, name index 2:
  #   logger::log_appender(logger::appender_file("my_log.log"),
  #                        namespace = "spatialkit", index = 2)
  #   logger::log_threshold(logger::FATAL, namespace = "spatialkit", index = 2)
  # spatialkit_quiet() does the second of those for you.

  # The FORMATTER is pinned, not inherited.  logger seeds a new namespace
  # from the user's global one, so the appender/threshold lines below left
  # the formatter to be whatever the user had set -- and every helper in
  # utils.R hands logger an ALREADY-formatted string.  Under the default
  # formatter_glue a `{` in a message was re-evaluated (a fold error reading
  # "diverged at {iter=3}" logged as "diverged at 3"); under a user's
  # formatter_sprintf every message containing a literal `%` -- the CRS
  # distortion figures, the GWR collinearity percentage -- hard-errored with
  # "too few arguments", and because .warn_and_log() logs before it warns,
  # the R warning the manual promises died with it.  formatter_paste does no
  # interpolation, so the message logged is the message written.
  logger::log_formatter(logger::formatter_paste, namespace = "spatialkit")

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

#' @noRd
.onUnload <- function(libpath) {
  # logger keeps the "spatialkit" namespace alive after unloadNamespace(),
  # with index 1 still pointing at this session's temp-file path.  Nothing
  # logs into it once the package is gone, but leave nothing dangling: both
  # appenders become no-ops.  .onLoad() re-registers them on the next load.
  # Guarded so that an unload can never fail on the logger's account.
  # (A plain no-op function rather than logger::appender_void, which older
  # logger releases do not export.)
  tryCatch({
    void <- function(lines) invisible(NULL)
    logger::log_appender(void, namespace = "spatialkit", index = 1)
    logger::log_appender(void, namespace = "spatialkit", index = 2)
  }, error = function(e) NULL)
  invisible(NULL)
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
#' @param quiet \code{TRUE} (default) silences the console echo;
#'   \code{FALSE} restores the package default, WARN+.  A \pkg{logger}
#'   threshold (\code{logger::ERROR}, or the value a previous call returned)
#'   sets that level instead, which is how to put back exactly what was in
#'   force rather than the default: \code{old <- spatialkit_quiet();
#'   spatialkit_quiet(old)}.
#' @return Invisibly, the threshold that was in force before the change --
#'   a \pkg{logger} level that can be passed back as \code{quiet}.
#' @family utilities
#' @examples
#' old <- spatialkit_quiet()      # console echo off
#' spatialkit_quiet(old)          # back to whatever it was
#' spatialkit_quiet(FALSE)        # or back to the WARN+ default
#' @export
spatialkit_quiet <- function(quiet = TRUE) {
  prev <- tryCatch(logger::log_threshold(namespace = "spatialkit", index = 2),
                   error = function(e) NULL)
  level <- if (inherits(quiet, "loglevel")) {
    quiet
  } else if (is.logical(quiet) && length(quiet) == 1L && !is.na(quiet)) {
    if (quiet) logger::FATAL else logger::WARN
  } else {
    stop("spatialkit_quiet(): `quiet` must be TRUE, FALSE or a logger threshold ",
         "such as logger::ERROR (or the value an earlier call returned).",
         call. = FALSE)
  }
  logger::log_threshold(level, namespace = "spatialkit", index = 2)
  invisible(prev)
}
