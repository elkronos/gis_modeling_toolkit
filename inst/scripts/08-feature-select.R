# =============================================================================
# 08  Choosing predictors without letting the choice flatter the score
# =============================================================================
#   source(system.file("scripts", "08-feature-select.R", package = "spatialkit"))
# =============================================================================
.tour_dir <- if (nzchar(system.file("scripts", package = "spatialkit")))
  system.file("scripts", package = "spatialkit") else "."
source(file.path(.tour_dir, "_common.R"))

pts <- tour_points()

# One candidate that drives the response, one that only looks as if it does, and
# four that are noise. `elev` is in the formula that generated `z`. `slope` is
# not, but this realisation of the field drifts west to east, and so does
# `slope`, so it predicts `z` by coincidence. `noise` and junk1 to junk3 give
# the search something to overfit to.
set.seed(11)
pts$junk1 <- stats::rnorm(nrow(pts))
pts$junk2 <- stats::rnorm(nrow(pts))
pts$junk3 <- stats::rnorm(nrow(pts))
cands <- c("elev", "slope", "noise", "junk1", "junk2", "junk3")
# Block size: a little under the range script 02 estimates on z (about 351 m).
# Blocks of 360 m make a 2 by 2 grid over the square, four blocks for five
# folds, and two blocks over the half-layer of 08.5; 300 m makes nine and
# six. Blocks this size leave some held-out points within the range of training
# data across their edges, so the scores below are the error of predicting a
# little closer to the data than range-sized blocks would hold them.
BS <- 300

xy <- sf::st_coordinates(pts)
cat(sprintf("  slope has no effect on z by construction, yet cor(slope, z) = %.2f,\n",
            stats::cor(pts$slope, pts$z)))
cat(sprintf("  because this field happens to drift west to east: cor(easting, z) = %.2f.\n",
            stats::cor(xy[, 1], pts$z)))
cat("  A variable search cannot tell that coincidence from a real effect.\n")

# select_features_forward() calls fit_fn(train_sf, predictor_vars), with two
# arguments. A learner written for cv_spatial() takes one, so it needs this
# wrapper. Without it every candidate scores exactly what the null model
# scored, and nothing is selected; the function then warns that fit_fn appears
# to be ignoring its second argument, which tells that case apart from data in
# which no candidate predicts the response.
lm_on <- function(train_sf, vars) {
  d <- sf::st_drop_geometry(train_sf)
  f <- stats::as.formula(paste("z ~",
        if (length(vars)) paste(vars, collapse = " + ") else "1"))
  new_spatial_fit(subclass = "lm_on", engine = stats::lm(f, data = d),
                  formula = f, response_var = "z", predictor_vars = vars,
                  data_sf = train_sf)
}
predict.lm_on <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(stats::fitted(object$engine))
  as.numeric(stats::predict(object$engine, newdata = sf::st_drop_geometry(newdata)))
}
registerS3method("predict", "lm_on", predict.lm_on)

step("08.1", "Forward selection, scored by blocked CV")
# Each step adds whichever remaining candidate improves cross-validated RMSE
# most, and stops when none of them improves it by more than `tol`.
sel <- select_features_forward(pts, "z", cands, fit_fn = lm_on, k = 5,
                               method = "block_kfold", block_size = BS,
                               metric = "RMSE", seed = 1, quiet = TRUE)
cat(sprintf("  selected: %s\n",
            if (length(sel$selected)) paste(sel$selected, collapse = ", ") else "(nothing)"))
cat(sprintf("  RMSE: %.3f (from %.3f with no predictors)\n",
            sel$score, sel$history$score[sel$history$step == 0L]))

step("08.2", "How far each pick beat the runner-up")
# `history` holds every candidate tried at every step. Read the gap between the
# step winner and the runner-up: a wide gap is a real predictor, a hair's
# breadth is the search choosing between noises.
h <- sel$history[sel$history$step > 0L, ]
for (st in sort(unique(h$step))) {
  hs <- h[h$step == st, ]
  hs <- hs[order(hs$score), ]
  cat(sprintf("  step %d: %-6s %.4f | runner-up %-6s %.4f | gap %.4f\n",
              st, hs$variable[1L], hs$score[1L], hs$variable[2L], hs$score[2L],
              hs$score[2L] - hs$score[1L]))
}
last <- h[h$step == max(h$step), ]
cat(sprintf("  The last step compares %d leftovers spread over %.4f of RMSE, and\n",
            nrow(last), diff(range(last$score))))
cat("  all of them are worse than the step before, so nothing was added. Had one\n",
    "  landed a thousandth the other way it would be in the model, and nothing\n",
    "  in the write-up would show how close the call had been.\n", sep = "")

if (!skip_without("ggplot2", "the selection plot")) {
  look_for("where the curve stops falling. Everything added to the right of ",
           "the flat part was chosen on a gap like the ones above.")
  show_plot(plot(sel), "08-selection.png")
}

step("08.3", "Did it take any of the decoys?")
decoys <- intersect(sel$selected, c("noise", "junk1", "junk2", "junk3"))
if (length(decoys)) {
  cat(sprintf("  yes: %s\n", paste(decoys, collapse = ", ")))
  cat("  Greedy search over six candidates on five folds picks up noise like\n",
      "  this; the margins in 08.2 show how little separated it.\n", sep = "")
} else {
  cat("  none this time, on these folds. That is one sample; 08.5 runs the\n",
      "  same search on half the data and gets a different answer.\n", sep = "")
}

step("08.4", "Tighten it with tol")
# `tol` is the improvement a candidate must exceed to be admitted, in the
# metric's own units: RMSE here, so units of `z`, not a percentage. Set it to
# the smallest improvement that would change a decision you make; it guards
# against the hairline gaps in 08.2 at no extra cost.
for (tl in c(0, 0.01, 0.05, 0.2)) {
  s <- select_features_forward(pts, "z", cands, fit_fn = lm_on, k = 5,
                               method = "block_kfold", block_size = BS,
                               metric = "RMSE", tol = tl, seed = 1, quiet = TRUE)
  cat(sprintf("  tol = %.2f -> %-20s RMSE %.3f\n", tl,
              if (length(s$selected)) paste(s$selected, collapse = ", ") else "(nothing)",
              s$score))
}
cat("  Set it from what would matter in the application, not from what leaves\n",
    "  the variables you were hoping for.\n", sep = "")

step("08.5", "The honest score for a selected model")
# The cross-validation that drove the search cannot also measure the winner.
# `select_on = "split"` runs the search on half the data, spatially split, and
# holds the other half out of the search to score with. Held out is not
# independent: the halves share a border, and rows near it are still
# correlated with the half the search saw.
# (It warns about fold imbalance: half the rows means uneven counts per block.)
sel_sp <- select_features_forward(pts, "z", cands, fit_fn = lm_on, k = 5,
                                  method = "block_kfold", block_size = BS,
                                  metric = "RMSE", seed = 1, quiet = TRUE,
                                  select_on = "split")
print(sel_sp$split)
cat(sprintf("  selected on the search half : %s\n",
            paste(sel_sp$selected, collapse = ", ")))
cat(sprintf("  score on the search half    : %.3f\n", sel_sp$score))
cat(sprintf("  score on the held-out half  : %.3f", sel_sp$score_holdout))
cat(sprintf("   (%.0f percent higher)\n",
            100 * (sel_sp$score_holdout / sel_sp$score - 1)))
cat("  Part of that gap is the selection effect. Part is that the two scores\n",
    "  come from different rows, a different part of the map and a different\n",
    "  amount of training data, so the gap is a rough estimate of the\n",
    "  selection effect.\n",
    sep = "")

step("08.6", "Two runs, two answers, one write-up")
cat(sprintf("  searched on all %d rows : %s\n", nrow(pts),
            paste(sel$selected, collapse = ", ")))
cat(sprintf("  searched on half of them: %s\n",
            paste(sel_sp$selected, collapse = ", ")))
cat("  Same layer, same learner, same seed; only the number of rows the search\n",
    "  saw changed. The variable list is not a stable property of the problem,\n",
    "  so report the procedure and its holdout score. A bare list of variables\n",
    "  presented as the ones that matter claims more than the search can\n",
    "  support.\n", sep = "")
