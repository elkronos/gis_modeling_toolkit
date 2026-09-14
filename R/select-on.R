# =============================================================================
# Sample splitting for post-selection inference: choose on one spatially
# blocked half of the data, estimate on the other.
# =============================================================================

#' Split a point layer into two spatially blocked halves
#'
#' \code{make_folds(k = 2, method = "block_kfold")} on the layer, so the two
#' halves are made of spatial blocks rather than of interleaved points: a
#' selection made on one half does not borrow its neighbours' values from the
#' other.  Fold 1 is the selection half, fold 2 the estimation half.  Rows
#' \code{make_folds()} drops (empty or non-finite geometry) belong to neither.
#'
#' @param data_sf The layer, already reduced to points.
#' @param seed Seed for the block assignment.
#' @param caller Name for messages.
#' @return A list with \code{selection} and \code{estimation} (integer row
#'   positions in \code{data_sf}), \code{method} and \code{seed}.
#' @keywords internal
#' @noRd
.spatial_half_split <- function(data_sf, seed = 123L, caller = "select_on") {
  n <- nrow(data_sf)
  if (n < 20L)
    stop(sprintf("%s(): select_on = \"split\" needs at least 20 points to make two spatial halves; got %d.",
                 caller, n), call. = FALSE)
  f <- tryCatch(
    suppressMessages(make_folds(data_sf, k = 2L, method = "block_kfold",
                                seed = if (is.null(seed)) 123L else seed)),
    error = function(e)
      stop(sprintf("%s(): could not split the layer into two spatial halves: %s",
                   caller, conditionMessage(e)), call. = FALSE))
  a <- f$assignment
  # Row IDs are positions unless the layer carried its own `..row_id`.
  pos <- if ("..row_id" %in% names(data_sf)) match(a$row_id, data_sf$..row_id) else
    as.integer(a$row_id)
  ok  <- is.finite(pos) & pos >= 1L & pos <= n
  sel <- sort(pos[ok & a$fold == 1L])
  est <- sort(pos[ok & a$fold != 1L])
  if (length(sel) < 10L || length(est) < 10L)
    stop(sprintf(paste0("%s(): the spatial split left %d and %d points in the two ",
                        "halves; at least 10 each are needed."),
                 caller, length(sel), length(est)), call. = FALSE)
  structure(list(selection = sel, estimation = est, method = "block_kfold",
                 seed = if (is.null(seed)) 123L else seed),
            class = "spatialkit_split")
}


#' @export
print.spatialkit_split <- function(x, ...) {
  cat(sprintf("Spatial half-split (%s, seed %s): %d selection rows, %d estimation rows\n",
              x$method, format(x$seed), length(x$selection), length(x$estimation)))
  cat("  $selection and $estimation are row positions in the layer as passed.\n")
  invisible(x)
}
