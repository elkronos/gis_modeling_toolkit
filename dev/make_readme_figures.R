# ---------------------------------------------------------------------------
# Regenerate the README figures in man/figures/ from package output, on the
# North Carolina state boundary shipped with sf (the vignettes use it too).
#
# Design: the raw observations have a panel of their own (coloured points);
# each tessellation panel shows cells only, coloured by a value computed with
# assign_features_to_polygons() and summarize_by_cell().  All panels share one
# colour scale, so the reader compares panels instead of decoding overlapping
# dots and cells.
#
# Run from the package root:
#   Rscript dev/make_readme_figures.R
#
# Requires: sf, ggplot2, and pkgload (or an installed spatialkit).
#
# The figures are drawn by the working tree, loaded with pkgload::load_all(),
# as the other dev scripts are; an installed spatialkit is used only when
# pkgload is missing, and then the figures show that release, not the code
# in front of you.
#
# The script prints the cell counts and leave-one-out errors the README
# quotes.  Copy them into README.Rmd whenever the figures are rebuilt.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(spatialkit)
  }
  library(sf)
  library(ggplot2)
})
cat("spatialkit loaded from: ",
    tryCatch(getNamespaceInfo(asNamespace("spatialkit"), "path"),
             error = function(e) NA_character_), "\n", sep = "")

dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
set.seed(42)

# --- North Carolina: real state boundary + simulated observation sites -----
nc      <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
nc_proj <- ensure_projected(nc)   # a CRS chosen for the extent: Albers here
state   <- st_sf(geometry = st_union(st_geometry(nc_proj)))

pts <- st_sf(geometry = st_sample(state, 280))        # observation sites

# Smooth spatial field: high in the west (mountains), declining eastward,
# with a hotspot in the east, plus a little observation noise.
crd <- st_coordinates(pts)
xs  <- (crd[, 1] - min(crd[, 1])) / diff(range(crd[, 1]))
ys  <- (crd[, 2] - min(crd[, 2])) / diff(range(crd[, 2]))
pts$val  <- 10 - 6 * xs +
  4 * exp(-((xs - 0.80)^2 + (ys - 0.45)^2) / 0.02) +
  rnorm(nrow(pts), sd = 0.3)

# Aggregate the point values per cell with the package's own workflow.
# estimate = "mean" gives the plain mean of each cell; "pooled", the default
# of summarize_by_cell(), gives the pooled value beside it, and the panel
# then shows the pooled value.
# Voronoi cells carry their IDs in cell_id, the grids in poly_id.
cell_means <- function(cells, label, estimate = "mean") {
  id       <- if ("poly_id" %in% names(cells)) "poly_id" else "cell_id"
  assigned <- assign_features_to_polygons(pts, cells, polygon_id_col = id)
  summ     <- summarize_by_cell(assigned, response_var = "val", id_col = id,
                                cells_sf = cells, estimate = estimate)
  value    <- if (estimate == "pooled") summ$..pooled_resp_val else summ$resp_mean_val
  st_sf(method   = label,
        mean_val = value,
        geometry = st_geometry(summ))
}

# Shared plotting helper: one "Raw observations" points panel + N cell panels.
panel_plot <- function(cells_sf, lvls) {
  cells_sf$method <- factor(cells_sf$method, levels = lvls)
  pts_panel <- pts
  pts_panel$method <- factor("Raw observations", levels = lvls)
  rng <- range(c(pts$val, cells_sf$mean_val), na.rm = TRUE)

  ggplot() +
    geom_sf(data = state, fill = "grey96", color = NA) +
    geom_sf(data = cells_sf, aes(fill = mean_val),
            color = "white", linewidth = 0.25) +
    geom_sf(data = pts_panel, aes(fill = val), shape = 21, size = 1.7,
            color = "grey25", stroke = 0.15) +
    geom_sf(data = state, fill = NA, color = "grey35", linewidth = 0.3) +
    scale_fill_viridis_c(limits = rng, name = "Value",
                         na.value = "grey88") +
    facet_wrap(~method, nrow = 2) +
    guides(fill = guide_colorbar(barwidth = 12, barheight = 0.5,
                                 title.vjust = 1)) +
    theme_void(base_size = 12) +
    theme(strip.text = element_text(face = "bold",
                                    margin = margin(b = 4, t = 6)),
          legend.position = "bottom")
}

# --- Figure 1: raw field + three tessellations (2 x 2) ----------------------
# The seeds are the default kind: density-weighted k-means from 3 starts.
seeds <- voronoi_seeds_kmeans(pts, k = 30)
vor   <- create_voronoi_polygons(seeds, boundary = state, quiet = TRUE)
hexg  <- create_grid_polygons(state, target_cells = 30, type = "hex")
sqg   <- create_grid_polygons(state, target_cells = 30, type = "square")

tess_cells <- rbind(
  cell_means(vor$cells, "Voronoi (density-weighted seeds)"),
  cell_means(hexg,      "Hex grid"),
  cell_means(sqg,       "Square grid")
)
lvls1 <- c("Raw observations", "Voronoi (density-weighted seeds)",
           "Hex grid", "Square grid")

cat(sprintf("figure 1: %d sites; %d Voronoi, %d hex and %d square cells\n",
            nrow(pts), nrow(vor$cells), nrow(hexg), nrow(sqg)))

ggsave("man/figures/readme-tessellations.png",
       panel_plot(tess_cells, lvls1),
       width = 11, height = 6, dpi = 150, bg = "white")

# --- Figure 2: random folds against block folds, side by side ---------------
folds_random <- make_folds(pts, k = 5, method = "random_kfold", seed = 7)
folds_block  <- make_folds(pts, k = 5, method = "block_kfold",
                           boundary = state, seed = 7)

pts_rand <- pts
pts_rand$fold   <- factor(folds_random$assignment$fold)
pts_rand$scheme <- "Random folds"
pts_block <- pts
pts_block$fold   <- factor(folds_block$assignment$fold)
pts_block$scheme <- "Block folds (whole blocks held out)"

cv_pts <- rbind(pts_rand, pts_block)
cv_pts$scheme <- factor(cv_pts$scheme,
                        levels = c("Random folds",
                                   "Block folds (whole blocks held out)"))

p2 <- ggplot(cv_pts) +
  geom_sf(data = state, fill = "grey96", color = "grey55",
          linewidth = 0.35) +
  geom_sf(aes(fill = fold), shape = 21, size = 2,
          color = "white", stroke = 0.25) +
  scale_fill_viridis_d(option = "turbo", end = 0.92, name = "CV fold") +
  facet_wrap(~scheme, nrow = 1) +
  guides(fill = guide_legend(override.aes = list(size = 3.5), nrow = 1)) +
  theme_void(base_size = 12) +
  theme(strip.text = element_text(face = "bold",
                                  margin = margin(b = 4, t = 6)),
        legend.position = "bottom")

ggsave("man/figures/readme-spatial-cv.png", p2,
       width = 11, height = 3.4, dpi = 150, bg = "white")

# --- Figure 3: three cell counts, the default among them -------------------
# The counts come from a ladder of sizes, in observations per cell: 23, a
# coarse count; 4.5, the size resolution_profile() and select_resolution()
# use by default with the default pooled estimate; and 2.  The level of a
# size s is max(2, round(n / s)) cells.  The coarse count stays at 10 cells
# or more: below that a cell's pooled value is pulled toward the mean of all
# the other cells, not of the cells near it.
#
# Each panel shows the value summarize_by_cell() returns by default, the
# pooled value, and the profile scores each count by the leave-one-out error
# of that value; the lines printed below give the errors the README quotes.
# The panels' cells are drawn with voronoi_seeds_kmeans() at the same
# counts, the default kind of seeds (density-weighted k-means).  The profile
# fits its own centres, so the errors belong to the count, not to the exact
# cells drawn.
sizes    <- c(23, 4.5, 2)
prof_res <- resolution_profile(pts, response_var = "val", obs_per_cell = sizes)
k_def    <- select_resolution(resolution_profile(pts, response_var = "val"))$best
ks       <- as.integer(prof_res$levels)
stopifnot(length(ks) == 3L, ks[2L] == k_def)
cat(sprintf("figure 3: %d cells (%s per cell): leave-one-out error %.2f (SE %.2f)\n",
            ks, as.character(sizes), prof_res$loo, prof_res$loo_se), sep = "")

labels <- sprintf("%d cells (%s per cell%s)", ks, as.character(sizes),
                  ifelse(ks == k_def, ", the default", ""))

res_cells <- do.call(rbind, Map(function(k, lab) {
  seeds_k <- voronoi_seeds_kmeans(pts, k = k)
  vor_k   <- create_voronoi_polygons(seeds_k, boundary = state, quiet = TRUE)
  cell_means(vor_k$cells, lab, estimate = "pooled")
}, ks, labels))
lvls3 <- c("Raw observations", labels)

ggsave("man/figures/readme-resolution.png",
       panel_plot(res_cells, lvls3),
       width = 11, height = 6, dpi = 150, bg = "white")

message("Wrote man/figures/readme-tessellations.png")
message("Wrote man/figures/readme-spatial-cv.png")
message("Wrote man/figures/readme-resolution.png")
