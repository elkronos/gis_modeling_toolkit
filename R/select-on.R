# =============================================================================
# Sample splitting for post-selection inference: choose on one spatially
# blocked half of the data, estimate on the other.
# =============================================================================

#' Split a point layer into two spatially blocked halves
#'
#' \code{make_folds(k = 2, method = "block_kfold")} on the layer, so the two
#' halves are made of spatial blocks rather than of interleaved points: far
#' fewer selection points sit next to an estimation point than under a random
#' half.  Blocking reduces the dependence between the halves; it does not
#' remove it.  The halves share a border, and under spatial autocorrelation
#' the points on either side of it are still correlated: on the tests'
#' \code{so_field(300)} (fitted range 253), 85\% of the estimation points lie
#' within the range of a selection point (93\% within the true range of 300;
#' the median distance to the nearest one is 141, against 36 for a random
#' half).  No buffer is cut between
#' the halves.  A buffer as wide as the range would discard about 85\% of
#' that estimation half, and a narrower one would be an arbitrary fraction of
#' a range the split does not estimate; the finite-sample exposure that
#' remains is the price of splitting a contiguous layer in two.  Fold 1 is
#' the selection half, fold 2 the estimation half.  Rows \code{make_folds()}
#' drops (empty or non-finite geometry) belong to neither.
#'
#' The blocks are the default grid of \code{make_folds()} (three per half)
#' unless a block design is passed through \code{...}.  That grid is fixed by
#' the layer's extent, so the seed does not change the partition: it decides
#' only which of the two sides is the selection half (and breaks ties in the
#' packing of equal blocks).  When the default grid leaves fewer than 10
#' points in a half -- a small layer, or a small group far from the rest,
#' which always lands in a block of its own -- the split is retried on
#' finer grids (16, 36, then 100 blocks) and the first that gives both halves
#' 10 points is used, with a warning: finer blocks lengthen the border the
#' halves share.  A block design passed through \code{...} is used as given.
#'
#' @param data_sf The layer, already reduced to points.
#' @param seed Seed for the block assignment.
#' @param caller Name for messages.
#' @param ... Block design passed to \code{make_folds()}: \code{block_nx},
#'   \code{block_ny}, \code{block_size}, \code{blocks}, \code{balance_tol}.
#' @return A list of class \code{"spatialkit_split"} with \code{selection}
#'   and \code{estimation} (integer row positions in \code{data_sf}),
#'   \code{method}, \code{seed}, \code{grid} (the block grid, \code{"nx x ny"})
#'   and \code{n_blocks} (the blocks of it that hold points),
#'   \code{balance} (the larger half's size over the smaller's; the default
#'   grid accepts up to 3, and clustered layers reach 2 routinely) and
#'   \code{extent} (the bounding box of each half, in \code{data_sf}'s CRS).
#' @keywords internal
#' @noRd
.spatial_half_split <- function(data_sf, seed = 123L, caller = "select_on", ...) {
  n <- nrow(data_sf)
  if (n < 20L)
    stop(sprintf("%s(): select_on = \"split\" needs at least 20 points to make two spatial halves; got %d.",
                 caller, n), call. = FALSE)
  seed   <- if (is.null(seed)) 123L else seed
  design <- list(...)

  # One make_folds() call, with the warnings it raises held back: a split
  # that is then retried on a finer grid must not leave behind the first
  # attempt's advice, and the kept attempt's warnings are raised as before.
  attempt <- function(extra, retry = FALSE) {
    held <- list()
    f <- withCallingHandlers(
      tryCatch(
        suppressMessages(do.call(make_folds, c(list(data_sf, k = 2L,
                                                    method = "block_kfold",
                                                    seed = seed),
                                               design, extra))),
        error = function(e) {
          if (retry) return(NULL)          # a finer grid that fails is no help
          stop(sprintf("%s(): could not split the layer into two spatial halves: %s",
                       caller, conditionMessage(e)), call. = FALSE)
        }),
      warning = function(w) {
        held[[length(held) + 1L]] <<- w
        invokeRestart("muffleWarning")
      })
    if (is.null(f)) return(list(sel = integer(0), est = integer(0),
                                n_blocks = NA_integer_, grid = NA_character_,
                                warnings = list()))
    a <- f$assignment
    # Row IDs are positions unless the layer carried its own `..row_id`.
    pos <- if ("..row_id" %in% names(data_sf)) match(a$row_id, data_sf$..row_id) else
      as.integer(a$row_id)
    ok  <- is.finite(pos) & pos >= 1L & pos <= n
    list(sel = sort(pos[ok & a$fold == 1L]), est = sort(pos[ok & a$fold != 1L]),
         n_blocks = as.integer(f$params$blocks_used %||% NA_integer_),
         grid = paste(f$params$grid_nx %||% NA, f$params$grid_ny %||% NA, sep = " x "),
         warnings = held)
  }
  too_small <- function(r) length(r$sel) < 10L || length(r$est) < 10L

  r <- attempt(list())
  first <- r
  # A remote group of fewer than 10 points always sat alone in a block of the
  # six-block default grid, so the call stopped although a finer grid splits
  # the same layer 131 / 127; nothing the caller could pass would reach it.
  # Finer blocks lengthen the border the halves share, so the retries are
  # few, capped, and announced, and a design the caller chose is not second-
  # guessed.
  if (too_small(r) && !length(design)) {
    for (mult in c(8, 18, 50)) {           # about 16, 36 and 100 blocks
      r <- attempt(list(block_multiplier = mult), retry = TRUE)
      if (!too_small(r)) break
    }
    if (!too_small(r))
      .warn_and_log(paste0("%s(): the default six-block split left %d and %d points ",
                           "in the two halves (at least 10 each are needed), so ",
                           "it was made on a finer %s grid instead. The ",
                           "halves then share a longer border, and more of the ",
                           "estimation points lie close to selection points."),
                    caller, length(first$sel), length(first$est), r$grid)
  }
  if (too_small(r))
    stop(sprintf(paste0("%s(): the spatial split left %d and %d points in the two ",
                        "halves; at least 10 each are needed%s."),
                 caller, length(first$sel), length(first$est),
                 if (!length(design)) ", and grids of up to 100 blocks did no better" else ""),
         call. = FALSE)
  for (w in r$warnings) warning(w)

  sizes  <- c(length(r$sel), length(r$est))
  extent <- lapply(list(selection = r$sel, estimation = r$est), function(i)
    sf::st_bbox(sf::st_geometry(data_sf)[i]))
  structure(list(selection = r$sel, estimation = r$est, method = "block_kfold",
                 seed = seed, grid = r$grid, n_blocks = r$n_blocks,
                 balance = max(sizes) / min(sizes), extent = extent),
            class = "spatialkit_split")
}


#' @export
print.spatialkit_split <- function(x, ...) {
  cat(sprintf("Spatial half-split (%s, seed %s): %d selection rows, %d estimation rows\n",
              x$method, format(x$seed), length(x$selection), length(x$estimation)))
  cat("  $selection and $estimation are row positions in the layer as passed.\n")
  invisible(x)
}
