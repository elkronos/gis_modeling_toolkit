#' How far the held-out points actually sit from the training data
#'
#' Blocked cross-validation exists to put distance between a test point and
#' the training points that could tell you its value.  Whether it succeeded is
#' a measurement, and until now the package only offered a proxy for it: the
#' block size compared against the estimated autocorrelation range, reported
#' by \code{\link{make_folds}()} as a warning.  That comparison is about the
#' \emph{design}.  This is about the \emph{result}: for every held-out point,
#' the distance to its nearest training point, summarised per fold.
#'
#' The two can disagree, and the direction is not obvious.  Blocks wider than
#' the range still leak wherever a test point sits near a block edge with
#' training data just across it, which is most of the points in a fine block
#' grid; conversely a fold whose blocks are narrower than the range can still
#' separate well if the points inside them are clustered.  The share of
#' held-out points closer to training data than the correlation range is the
#' number that settles it, and it is the last column here.
#'
#' Nothing is estimated: the distances come from the geometry, and the range,
#' when one is shown, is the one the folds already carry (\code{make_folds(
#' auto_range = TRUE)} records it) or the one you pass as \code{sac}.
#'
#' @param folds A \code{\link{make_folds}()} result, or its \code{$folds}
#'   element (a list of \code{train}/\code{test} splits).
#' @param data_sf The layer the folds were built on.  Row identifiers are
#'   matched through \code{..row_id} when the layer carries one, and by row
#'   position otherwise, which is what \code{make_folds()} and every
#'   \code{cv_*()} do.
#' @param sac Optional: an \code{\link{estimate_sac_range}()} result or a
#'   single number, in the CRS units of \code{data_sf}.  Defaults to the range
#'   the folds carry, if any.  Supplying one adds the \code{within_range}
#'   column and the closing verdict.
#' @return A data.frame of class \code{fold_separation}, one row per fold:
#'   \code{fold}, \code{n_train}, \code{n_test}, \code{n_blocks} (\code{NA}
#'   for a scheme with no blocks), \code{min_dist} and \code{median_dist}
#'   (distance from a held-out point to its nearest training point, in CRS
#'   units), and \code{within_range} (the share of held-out points closer to
#'   training data than \code{sac}; \code{NA} without one).  Attributes:
#'   \code{method}, \code{sac_range}, \code{crs} and \code{n_unknown_ids}.
#' @family cross-validation
#' @seealso \code{\link{make_folds}()} for the fold schemes and the block
#'   sizing this measures the outcome of; \code{\link{cv_block_size_sweep}()}
#'   for choosing a block size by cross-validated error instead.
#' @examples
#' library(sf)
#' set.seed(1)
#' n <- 200
#' pts <- st_as_sf(
#'   data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000)),
#'   coords = c("x", "y"), crs = 32632
#' )
#'
#' # Random folds put a training point almost on top of every held-out one.
#' random  <- make_folds(pts, k = 4, method = "random_kfold", seed = 1)
#' print(fold_separation(random, pts, sac = 200))
#'
#' # Blocked folds hold out whole neighbourhoods, so the distances grow.
#' blocked <- make_folds(pts, k = 4, method = "block_kfold",
#'                       block_size = 250, seed = 1)
#' print(fold_separation(blocked, pts, sac = 200))
#' @export
fold_separation <- function(folds, data_sf, sac = NULL) {
  if (!inherits(data_sf, "sf"))
    stop("fold_separation(): `data_sf` must be an sf object.", call. = FALSE)
  splits <- if (is.list(folds) && !is.null(folds$folds)) folds$folds else folds
  if (!is.list(splits) || !length(splits) ||
      !all(vapply(splits, function(s) is.list(s) && !is.null(s$test), logical(1))))
    stop("fold_separation(): `folds` must be a make_folds() result or its ",
         "`$folds` element (a list of train/test splits).", call. = FALSE)

  # Same row-identity convention as make_folds() and every cv_*(): the stamped
  # column when there is one, row position otherwise.  Getting this wrong
  # would silently measure the distance between the wrong pairs of points.
  ids <- if ("..row_id" %in% names(data_sf)) data_sf[["..row_id"]] else
    seq_len(nrow(data_sf))

  pts <- data_sf
  if (!all(sf::st_geometry_type(pts, by_geometry = TRUE) == "POINT"))
    pts <- coerce_to_points(pts, "auto")
  pts <- sf::st_zm(ensure_projected(pts), drop = TRUE, what = "ZM")

  # A distance is only meaningful between finite coordinates.
  xy <- sf::st_coordinates(pts)[, 1:2, drop = FALSE]
  ok <- stats::complete.cases(xy) & is.finite(xy[, 1L]) & is.finite(xy[, 2L])
  if (!any(ok))
    stop("fold_separation(): `data_sf` has no usable coordinates.", call. = FALSE)

  sac_val <- if (!is.null(sac)) suppressWarnings(as.numeric(sac)[1L]) else
    suppressWarnings(as.numeric(folds$params$sac_range %||% NA_real_)[1L])
  if (!length(sac_val) || !is.finite(sac_val)) sac_val <- NA_real_

  fold_blocks <- folds$params$fold_blocks
  unknown <- 0L

  rows <- lapply(seq_along(splits), function(j) {
    s  <- splits[[j]]
    te <- match(s$test,  ids)
    tr <- match(s$train, ids)
    unknown <<- unknown + sum(is.na(te)) + sum(is.na(tr))
    te <- te[!is.na(te) & ok[te]]
    tr <- tr[!is.na(tr) & ok[tr]]
    d <- if (length(te) && length(tr)) .nn_dist_to(pts[te, ], pts[tr, ]) else numeric(0)
    d <- d[is.finite(d)]
    data.frame(
      fold         = j,
      n_train      = length(tr),
      n_test       = length(te),
      n_blocks     = if (is.list(fold_blocks) && length(fold_blocks) >= j)
                       length(fold_blocks[[j]]) else NA_integer_,
      min_dist     = if (length(d)) min(d) else NA_real_,
      median_dist  = if (length(d)) stats::median(d) else NA_real_,
      within_range = if (length(d) && is.finite(sac_val)) mean(d < sac_val) else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  if (unknown > 0L)
    .log_info(paste0("fold_separation(): %d fold entr(ies) name rows `data_sf` ",
                     "does not have, and were skipped. Are these the folds ",
                     "this layer was split into?"), unknown)
  structure(out, class = c("fold_separation", "data.frame"),
            method = folds$method %||% "supplied splits",
            sac_range = sac_val, n_unknown_ids = unknown,
            crs = sf::st_crs(pts)$input %||% NA_character_)
}


#' @export
print.fold_separation <- function(x, ...) {
  sac <- attr(x, "sac_range"); crs <- attr(x, "crs")
  cat(sprintf("Fold separation: %s, %d fold(s), %d held-out point(s)%s\n",
              attr(x, "method"), nrow(x), sum(x$n_test, na.rm = TRUE),
              if (is.na(crs)) "" else sprintf(" (%s)", crs)))
  if (is.finite(sac))
    cat(sprintf("  autocorrelation range: %s\n", format(sac, digits = 4)))
  df <- as.data.frame(x)
  df$min_dist    <- signif(df$min_dist, 4)
  df$median_dist <- signif(df$median_dist, 4)
  if (all(is.na(df$within_range))) df$within_range <- NULL
  else df$within_range <- sprintf("%.0f%%", 100 * df$within_range)
  if (all(is.na(df$n_blocks))) df$n_blocks <- NULL
  cat("\n"); print(df, row.names = FALSE); cat("\n")
  if (is.finite(sac)) {
    share <- stats::weighted.mean(x$within_range, x$n_test, na.rm = TRUE)
    closest <- suppressWarnings(min(x$min_dist, na.rm = TRUE))
    # The verdict is about leakage, so it is stated in terms of what leaks:
    # a held-out point closer to training data than the range is one whose
    # value the training set partly carries.
    cat(strwrap(sprintf(paste0("%.0f%% of held-out points sit closer to a ",
                               "training point than the correlation range ",
                               "(%s), and the closest is %s away. %s"),
                        100 * share, format(sac, digits = 4),
                        format(signif(closest, 3)),
                        if (share > 0.5)
                          paste("Most of the hold-out is inside the range of",
                                "its own training data, so this score is",
                                "optimistic: widen the blocks.")
                        else if (share > 0.1)
                          paste("A minority leaks, which is the usual price of",
                                "contiguous blocks at the edges.")
                        else
                          paste("Little of the hold-out is within reach of",
                                "the training data.")),
                 width = 74, prefix = "  "), sep = "\n")
  }
  invisible(x)
}
