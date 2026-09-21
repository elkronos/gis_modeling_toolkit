# =============================================================================
# 02  Choosing a cell size, and saying why you chose it
# =============================================================================
#   source(system.file("scripts", "02-resolution.R", package = "spatialkit"))
# =============================================================================
.tour_dir <- if (nzchar(system.file("scripts", package = "spatialkit")))
  system.file("scripts", package = "spatialkit") else "."
source(file.path(.tour_dir, "_common.R"))

if (skip_without("gstat", "the whole of script 02")) {
  cat("  The range is what the ladder's floor and two of its four criteria come\n",
      "  from, and gstat fits the variogram behind it.\n", sep = "")
} else {

pts <- tour_points()

step("02.1", "How far does the field reach?")
# Every resolution question is really a range question: cells much smaller than
# the range mostly measure the same thing twice, cells much larger average the
# signal away. Start by asking the data.
rng <- suppressWarnings(estimate_sac_range(pts, "z", seed = 1))
print(rng)
if (!skip_without("ggplot2", "the variogram")) {
  look_for("the fitted curve flattening out. Where it flattens is the range; ",
           "the height it starts from is the nugget, which is noise you cannot ",
           "average away by making cells bigger.")
  show_plot(plot(rng), "02-variogram.png")
}

step("02.2", "Candidate cell counts, cheapest first")
# determine_optimal_levels() only looks at geometry: it clusters the points and
# reports the level counts where the within-cluster spread stops improving.
# Seconds, not minutes, and it is the right first move.
lv <- determine_optimal_levels(pts, max_levels = 40)
cat("  geometric candidates:", paste(lv, collapse = ", "), "\n")

# With a response and predictors it can also score levels on residual spatial
# autocorrelation. `criterion = "combined"` REQUIRES both.
lv2 <- determine_optimal_levels(pts, max_levels = 40, criterion = "combined",
                                response_var = "z", predictor_vars = "elev")
cat("  combined candidates:", paste(lv2, collapse = ", "), "\n")

step("02.3", "The full ladder")
# resolution_profile() actually builds each level and reports what it costs:
# how many points land in the smallest cell, how much of the variance survives,
# and whether the residuals are still autocorrelated.
prof <- resolution_profile(pts, response_var = "z", n_levels = 16)
print(prof)

step("02.4", "Four criteria, four answers")
# Each criterion also reports the band of levels it cannot tell apart from its
# own pick. The spread BETWEEN criteria is usually wider than any one of those
# bands, and it is the more honest measure of how much the data pins this down.
for (cr in c("cp", "reliability", "elbow", "moran_z")) {
  s <- try(select_resolution(prof, criterion = cr), silent = TRUE)
  if (inherits(s, "try-error")) {
    cat(sprintf("  %-12s not available on this profile\n", cr))
    next
  }
  # `at_floor` / `at_ceiling` mean the optimum is the end of the ladder, so the
  # ladder chose and the criterion did not. Widen n_levels and run it again.
  edge <- if (isTRUE(s$at_floor)) "  <- at the floor of the ladder"
          else if (isTRUE(s$at_ceiling)) "  <- at the ceiling of the ladder" else ""
  # The flat region is a SET: it can skip a rung, so printing its first two
  # members as "a to b" both truncated it and implied the levels between were
  # in it.  Name the others, or count them when there are too many to read.
  others <- setdiff(s$flat, s$best)
  band <- if (!length(others)) " (no tie)"
          else if (length(others) <= 4L)
            sprintf(" (tied with %s)", paste(others, collapse = ", "))
          else sprintf(" (tied with %d other levels)", length(others))
  cat(sprintf("  %-12s %2d cells%s%s\n", cr, s$best, band, edge))
}
cat("  Pick one before you look, and say which one you picked.\n")

if (!skip_without("ggplot2", "the profile plot")) {
  look_for("where the curves stop moving. A criterion whose optimum sits at ",
           "the first or last level of the ladder did not choose: the ladder ",
           "did. Widen n_levels and run it again.")
  show_plot(plot(prof), "02-profile.png", height = 7)
}

step("02.5", "The leak: choosing the resolution on data you will test on")
# Default `select_on = "all"` uses every row. That is fine when the resolution
# is a reporting decision. It is NOT fine when you will later report a test
# metric at that resolution, because the test rows helped choose it.
prof_split <- resolution_profile(pts, response_var = "z", n_levels = 16,
                                 select_on = "split")
sp <- attr(prof_split, "split")
if (!is.null(sp)) print(sp)
b_all   <- select_resolution(prof, "reliability")$best
b_split <- select_resolution(prof_split, "reliability")$best
cat(sprintf("  select_on = 'all'   -> %d cells\n", b_all))
cat(sprintf("  select_on = 'split' -> %d cells\n", b_split))
if (b_all == b_split) {
  cat("  They agree, so the choice did not depend on the rows you will test on.\n")
} else {
  cat("  They disagree, so the resolution was tuned to rows you were planning\n",
      "  to test on, and the test is no longer independent of the choice.\n",
      sep = "")
}

step("02.6", "What the chosen resolution looks like")
if (!skip_without("ggplot2", "the maps")) {
  bnd <- tour_boundary(pts)
  chosen <- select_resolution(prof, "elbow")$best
  for (k in sort(unique(c(min(prof$levels), chosen, max(prof$levels))))) {
    seeds <- get_voronoi_seeds(boundary = bnd, method = "kmeans", n = k,
                              sample_points = pts, set_seed = 1)
    tess  <- build_tessellation(seeds, boundary = bnd, method = "voronoi",
                                quiet = TRUE)
    cel   <- summarize_by_cell(assign_features_to_polygons(pts, tess$cells),
                               "z", cells_sf = tess$cells, deff = 1)
    note <- if (k == chosen)
      "the elbow pick: enough cells to show the field, few enough to fill."
    else if (k == min(prof$levels))
      "the coarsest level: every cell is well filled, but there is barely a map left."
    else
      paste("the finest level: plenty of detail, and cells thin enough that a",
            "single point moves one.")
    look_for(sprintf("%d cells -- %s", k, note))
    show_plot(plot_tessellation_map(cel, boundary = bnd, fill_col = "resp_mean_z",
                                    title = sprintf("%d cells: mean of z", k)),
              sprintf("02-cells-%02d.png", k))
  }
}

}
