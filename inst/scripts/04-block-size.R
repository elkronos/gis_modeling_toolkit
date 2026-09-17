# =============================================================================
# 04  Defending a block size instead of guessing one
# =============================================================================
#   source(system.file("scripts", "04-block-size.R", package = "spatialkit"))
# =============================================================================
.tour_dir <- if (nzchar(system.file("scripts", package = "spatialkit")))
  system.file("scripts", package = "spatialkit") else "."
source(file.path(.tour_dir, "_common.R"))

pts <- tour_points()

# An inverse-distance interpolator on the five nearest samples. No optional
# dependencies, and it behaves the way a spatial model is supposed to: good
# where neighbours are close, useless where they are not.
idw_fit <- function(train_sf, ...) {
  new_spatial_fit(
    subclass = "idw_fit",
    engine = list(xy = sf::st_coordinates(train_sf), z = train_sf$z, k = 5L),
    formula = z ~ 1, response_var = "z", predictor_vars = character(0),
    data_sf = train_sf)
}
predict.idw_fit <- function(object, newdata = NULL, ...) {
  e <- object$engine
  if (is.null(newdata)) return(e$z)
  q <- sf::st_coordinates(newdata)
  vapply(seq_len(nrow(q)), function(i) {
    d <- sqrt((e$xy[, 1] - q[i, 1])^2 + (e$xy[, 2] - q[i, 2])^2)
    o <- order(d)[seq_len(min(e$k, length(d)))]
    w <- 1 / pmax(d[o], 1e-9)
    sum(w * e$z[o]) / sum(w)
  }, numeric(1))
}
registerS3method("predict", "idw_fit", predict.idw_fit)

# Report where a sweep stops changing: the first size past which every further
# widening moves the score by less than `tol`.
shape <- function(x, label, tol = 0.03) {
  d   <- as.data.frame(x)
  rnd <- d$value[is.na(d$block_size)]
  blk <- d[!is.na(d$block_size), ]
  rel <- diff(blk$value) / utils::head(blk$value, -1L)
  ok  <- vapply(seq_along(rel), function(i) all(rel[i:length(rel)] < tol), logical(1))
  knee <- which(ok)[1L]
  cat(sprintf("  %-16s random %.3f | blocked %.3f at %.0f units -> %.3f at %.0f\n",
              label, rnd, blk$value[1L], blk$block_size[1L],
              blk$value[nrow(blk)], blk$block_size[nrow(blk)]))
  if (is.na(knee)) {
    cat("    never flattens: still climbing at the widest size tried\n")
  } else if (knee == 1L) {
    cat("    flat throughout: this score does not depend on block size\n")
  } else {
    cat(sprintf("    flattens from %.0f units on (%.0f%% of the %.0f-unit range)\n",
                blk$block_size[knee],
                100 * blk$block_size[knee] / attr(x, "sac_range"),
                attr(x, "sac_range")))
  }
  invisible(blk)
}

step("04.1", "Sweep the size, do not pick one")
# Script 03 used block_size = 250 with no justification. A reviewer is entitled
# to ask why. The sweep refits at several sizes and shows what the score does,
# with a random-fold run on the same data as the leaky reference to compare to.
sw <- cv_block_size_sweep(pts, "z", character(0), fit_fn = idw_fit,
                          n_sizes = 6, k = 5, seed = 1, quiet = TRUE)
print(sw)

step("04.2", "Reading the shape")
shape(sw, "interpolator")
cat("  A climb followed by a flat stretch is the shape to hope for. While the\n",
    "  blocks are narrower than the autocorrelation, test points still have\n",
    "  training neighbours nearby and the score is optimistic. Once the blocks\n",
    "  are wide enough the leak is gone, and widening them further changes\n",
    "  nothing. Report a size from the flat stretch and say the score held\n",
    "  across it. That sentence is what the sweep is for.\n", sep = "")

step("04.3", "The picture for the appendix")
if (!skip_without("ggplot2", "the sweep plot")) {
  look_for("where the curve flattens, and where that is relative to the dotted ",
           "red range line. The dashed grey line is the random-fold score, so ",
           "the vertical gap up to the flat part is what honest blocking costs.")
  show_plot(plot(sw), "04-sweep-interpolator.png")
}

step("04.4", "When the curve never flattens")
# `tour_fit` is a cubic trend surface in x and y. Widening the blocks asks it to
# extrapolate that surface further and further, so its score keeps climbing and
# there is no size at which it settles.
sw2 <- cv_block_size_sweep(pts, "z", "elev", fit_fn = tour_fit,
                           n_sizes = 6, k = 5, seed = 1, quiet = TRUE)
shape(sw2, "trend surface")
if (!skip_without("ggplot2", "the second sweep plot")) {
  look_for("no flat stretch anywhere: every widening costs more than the last.")
  show_plot(plot(sw2), "04-sweep-trend.png")
}
cat("  A curve like that is not a statement about block size. It says the model\n",
    "  degrades smoothly with distance from its training data, so there is no\n",
    "  single honest number -- the score depends entirely on how far you intend\n",
    "  to predict. Script 07 measures that distance directly.\n", sep = "")

step("04.5", "What this costs")
# Every size is a full k-fold CV: n_sizes * k fits, plus k more for the random
# baseline. `max_fits` (default 60) is the brake. Ask for more than that and the
# sweep refuses and names the budget it would have needed, so a wide sweep with
# a slow learner fails in a second instead of running for an hour.
cat(sprintf("  %d fits per sweep (%d sizes at k = %d, metric %s)\n",
            attr(sw, "n_fits"), sum(!is.na(as.data.frame(sw)$block_size)),
            attr(sw, "k"), attr(sw, "metric")))
cat("  With a slow learner start at n_sizes = 3, and widen only if the curve\n",
    "  has not flattened by the last size.\n", sep = "")

step("04.6", "The rows to distrust: sizes that leave too few blocks")
# Wide sizes can run fewer folds than you asked for (script 03.5). The `k` and
# `n_folds_succeeded` columns exist so a row resting on two folds does not sit
# unlabelled beside a row resting on five.
d <- as.data.frame(sw)
short <- d[!is.na(d$block_size) & d$k < attr(sw, "k"), ]
if (nrow(short)) {
  cat(sprintf("  %d size(s) ran fewer than k = %d folds:\n", nrow(short),
              attr(sw, "k")))
  print(short[, c("block_size", "blocks_used", "k", "n_folds_succeeded", "value")],
        row.names = FALSE, digits = 4)
  cat("  Weaker evidence, not the answer.\n")
} else {
  cat(sprintf("  every size ran the full k = %d here. Check the column anyway:\n",
              attr(sw, "k")))
  cat("  it is the first thing to go on a smaller or more clustered sample.\n")
}
