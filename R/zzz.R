#' @noRd
.onLoad <- function(libname, pkgname) {
  # Set up default logging in a package-specific namespace so we never
  # overwrite the user's global logger configuration -- and, the other way
  # round, so the user's global configuration cannot break ours (see the
  # notes below).
  # Users can reconfigure the spatialkit namespace after loading -- but note
  # that logger::log_appender() and log_threshold() BOTH default to index = 1,
  # so the two-line recipe below without an index touches only the temp-file
  # appender and leaves the console echo exactly as it was.  To redirect or
  # quieten what is printed, name index 2:
  #   logger::log_appender(logger::appender_file("my_log.log"),
  #                        namespace = "spatialkit", index = 2)
  #   logger::log_threshold(logger::FATAL, namespace = "spatialkit", index = 2)
  # spatialkit_quiet() does the second of those for you.

  # EVERY setting of both indices is pinned, not inherited.  logger seeds a
  # new namespace by copying the user's whole global configuration -- every
  # index, each with its formatter, layout, appender and threshold -- and the
  # lines below overwrite only what they name.  Pinning the formatter on
  # index 1 alone left index 2 with the user's: someone who had configured two
  # global indices before loading the package got formatter_sprintf or
  # formatter_glue on the console echo, and a `%` or a `{` in a message
  # ("fold 2 skipped: object 'cov_{x' not found") aborted the function that
  # logged it, taking the R warning .warn_and_log() promises with it.  The
  # helpers in utils.R also mark every message skip_formatter(), so no
  # formatter on any index sees it; formatter_paste, which does no
  # interpolation, is the second line of defence.
  for (i in 1:2) {
    logger::log_formatter(logger::formatter_paste, namespace = "spatialkit",
                          index = i)
    logger::log_layout(logger::layout_simple, namespace = "spatialkit",
                       index = i)
  }

  # Index 1: full INFO+ trace to a session temp file (detailed diagnostics).
  # The path is resolved per line, not here; see .sk_file_appender().
  logger::log_appender(.sk_file_appender(), namespace = "spatialkit",
                       index = 1)
  logger::log_threshold(logger::INFO, namespace = "spatialkit", index = 1)

  # Index 2: WARN+ to the console so that important problems (skipped CV
  # folds, failed predictions returning NA, extraction failures, ...) are
  # actually visible to interactive users instead of only landing in a
  # temp file nobody reads.  While a document is knitted the line is also
  # sent as an R message; see .sk_console_appender().
  logger::log_appender(.sk_console_appender, namespace = "spatialkit",
                       index = 2)
  logger::log_threshold(logger::WARN, namespace = "spatialkit", index = 2)

  # Indices 3 and up can only have been copied from the user's global
  # configuration, and they kept the user's appenders: spatialkit's WARN and
  # INFO lines landed in the user's own log files.
  .sk_logger_drop_copied_indices("spatialkit")
}


# The temp-file trace (index 1).  The path is resolved when a line is written,
# not when the package loads.  A session temp directory deleted under a
# running session -- an OS cleaner, or unlink(tempdir()) -- used to turn every
# call that logs into "cannot open the connection", and because
# .warn_and_log() logs before it warns, a documented R warning became that
# error; tempdir(check = TRUE), R's own recovery, did not help, because the
# old path was fixed at load time.  The directory is recreated if it has gone
# (R's tempdir() still names it, so R still cleans it up at exit), and a line
# that still cannot be written is dropped: the trace is a diagnostic, never a
# reason for the computation to fail.  The "generator" attribute is what
# logger's getter reports for the appender, as it does for appender_file().
.sk_file_appender <- function(path = function()
                                file.path(tempdir(), "spatialkit_model_log.log")) {
  force(path)
  structure(function(lines) {
    tryCatch(suppressWarnings({
      f <- path()
      if (!dir.exists(dirname(f)))
        dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
      cat(lines, sep = "\n", file = f, append = TRUE)
    }), error = function(e) NULL)
    invisible(NULL)
  }, generator = ".sk_file_appender()")
}

# The console echo (index 2).  logger's appender_console writes to stderr with
# cat(), which knitr does not capture, so in a knitted R Markdown, Quarto or
# pkgdown document every log-only caution vanished while the R warnings next
# to it were shown.  The line still goes to stderr exactly as before, and
# while knitr is running it is ALSO sent as an R message, which the document
# shows and the chunk option `message = FALSE` hides; nothing that reached the
# console before is lost.  A line that is about to be raised as an R warning
# as well -- .warn_and_log() -- is not repeated as a message, since the
# document already shows the warning.
.sk_console_appender <- function(lines) {
  cat(lines, file = stderr(), sep = "\n")
  if (isTRUE(getOption("knitr.in.progress")) && !isTRUE(.sk_log_state$raising))
    message(paste(lines, collapse = "\n"))
  invisible(NULL)
}

# Reduce the namespace to the two indices above.  logger 0.2.2 has no public
# way to count or delete indices (later releases export delete_logger_index()),
# but its getter answers for any index past the last with the LAST index's
# settings.  Index 2's appender is this package's own function, so reading it
# back at index 3 means there is no index 3.  Without delete_logger_index() an
# extra index is switched off instead: a no-op appender, and a threshold below
# FATAL, which no log line meets.  logger 0.2.2's setters cannot address an
# index above 5, so neither can this.
.sk_logger_drop_copied_indices <- function(ns) {
  appender_at <- function(k) logger::log_appender(namespace = ns, index = k)
  has_index <- function(k) !identical(appender_at(k), appender_at(k - 1L))
  del <- tryCatch(getExportedValue("logger", "delete_logger_index"),
                  error = function(e) NULL)
  if (is.function(del)) {
    n <- 0L
    while (has_index(3L) && n < 100L) {
      del(namespace = ns, index = 3L)
      n <- n + 1L
    }
    return(invisible(NULL))
  }
  off <- structure(0L, level = "OFF", class = "loglevel")
  for (k in 3:5) {
    if (!has_index(k)) break
    logger::log_appender(.sk_log_off, namespace = ns, index = k)
    logger::log_threshold(off, namespace = ns, index = k)
  }
  invisible(NULL)
}

.sk_log_off <- function(lines) invisible(NULL)

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
#' leaves the console untouched, which is the opposite of what anyone wants.
#' This helper names the right index.
#'
#' These are log records, not R conditions: \code{suppressWarnings()}
#' and \code{tryCatch(warning = )} do not see them.  Conditions the package
#' raises as real R warnings are unaffected by this function.
#'
#' While a document is being knitted (R Markdown, Quarto, a \pkg{pkgdown}
#' article) the console echo is also sent as an R message, because
#' \pkg{knitr} does not capture what is written to the console's error stream
#' and the cautions would otherwise be missing from the output.  They appear
#' as \code{## WARN [...]} lines, and the chunk option \code{message = FALSE},
#' like \code{suppressMessages()}, keeps them out of the document.  A line the
#' package also raises as an R warning is not repeated, since the document
#' shows the warning.  Outside \pkg{knitr} nothing changes.
#'
#' @param quiet \code{TRUE} (default) silences the console echo;
#'   \code{FALSE} restores the package default, WARN+.  A \pkg{logger}
#'   threshold (\code{logger::ERROR}, or the value a previous call returned)
#'   sets that level instead, which is how to put back exactly what was in
#'   force rather than the default: \code{old <- spatialkit_quiet();
#'   spatialkit_quiet(old)}.
#' @return Invisibly, the threshold that was in force before the change, a
#'   \pkg{logger} level that can be passed back as \code{quiet}.
#' @family package options and caches
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
