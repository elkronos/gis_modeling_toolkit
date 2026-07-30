# ---------------------------------------------------------------------------
# Regenerate the README figures in man/figures/ using actual package output,
# on the North Carolina state boundary (the classic dataset shipped with sf,
# also used by the package vignette).
#
# Design: raw observations are shown in their own panel (coloured points);
# tessellation panels show cells only, coloured by the per-cell mean of the
# points that fall inside them (computed with assign_features_to_polygons()
# + summarize_by_cell()).  All panels share one colour scale, so the reader
# compares panels instead of decoding overlapping dots and cells.
#
# Run from the package root:
#   Rscript dev/make_readme_figures.R
#
# Requires: sf, ggplot2 (and devtools if spatialkit is not installed).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  if (requireNamespace("spatialkit", quietly = TRUE)) {
    library(spatialkit)
  } else {
    devtools::load_all(".", quiet = TRUE)
  }
  library(sf)
  library(ggplot2)
})

dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
set.seed(42)

# --- North Carolina: real state boundary + simulated observation sites -----
nc      <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
nc_proj <- ensure_projected(nc)                       # auto-picks UTM 17N
state   <- st_sf(geometry = st_union(st_geometry(nc_proj)))

pts <- st_sf(geometry = st_sample(state, 280))        # observation sites

# Smooth spatial field: high in the west (mountains), declining eastward,
# with a hotspot in the east — plus a little observation noise.
crd <- st_coordinates(pts)
xs  <- (crd[, 1] - min(crd[, 1])) / diff(range(crd[, 1]))
ys  <- (crd[, 2] - min(crd[, 2])) / diff(range(crd[, 2]))
pts$val  <- 10 - 6 * xs +
  4 * exp(-((xs - 0.80)^2 + (ys - 0.45)^2) / 0.02) +
  rnorm(nrow(pts), sd = 0.3)
pts$west <- xs   # predictor used by the model-aware resolution selection

# Aggregate the point values per cell with the package's own workflow.
cell_means <- function(cells, label) {
  assigned <- assign_features_to_polygons(pts, cells)
  summ     <- summarize_by_cell(assigned, response_var = "val",
                                cells_sf = cells)
  st_sf(method   = label,
        mean_val = summ$resp_mean_val,
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
seeds <- voronoi_seeds_kmeans(pts, k = 30)
vor   <- create_voronoi_polygons(seeds, boundary = state, quiet = TRUE)
hexg  <- create_grid_polygons(state, target_cells = 30, type = "hex")
sqg   <- create_grid_polygons(state, target_cells = 30, type = "square")

tess_cells <- rbind(
  cell_means(vor$cells, "Voronoi (k-means seeds)"),
  cell_means(hexg,      "Hex grid"),
  cell_means(sqg,       "Square grid")
)
lvls1 <- c("Raw observations", "Voronoi (k-means seeds)",
           "Hex grid", "Square grid")

ggsave("man/figures/readme-tessellations.png",
       panel_plot(tess_cells, lvls1),
       width = 11, height = 6, dpi = 150, bg = "white")

# --- Figure 2: random folds leak, spatial block folds cut geography --------
folds_random <- make_folds(pts, k = 5, method = "random_kfold", seed = 7)
folds_block  <- make_folds(pts, k = 5, method = "block_kfold",
                           boundary = state, seed = 7)

pts_rand <- pts
pts_rand$fold   <- factor(folds_random$assignment$fold)
pts_rand$scheme <- "Random k-fold (test points leak spatial signal)"
pts_block <- pts
pts_block$fold   <- factor(folds_block$assignment$fold)
pts_block$scheme <- "Spatial block k-fold (whole regions held out)"

cv_pts <- rbind(pts_rand, pts_block)
cv_pts$scheme <- factor(cv_pts$scheme,
                        levels = c("Random k-fold (test points leak spatial signal)",
                                   "Spatial block k-fold (whole regions held out)"))

p2 <- ggplot(cv_pts) +
  geom_sf(data = state, fill = "grey96", color = "grey55",
          linewidth = 0.35) +
  geom_sf(aes(fill = fold), shape = 21, size = 2,
          color = "white", stroke = 0.25) +
  scale_fill_viridis_d(option = "turbo", end = 0.92, name = "CV fold") +
  facet_wrap(~scheme, nrow = 2) +
  guides(fill = guide_legend(override.aes = list(size = 3.5), nrow = 1)) +
  theme_void(base_size = 12) +
  theme(strip.text = element_text(face = "bold",
                                  margin = margin(b = 4, t = 6)),
        legend.position = "bottom")

ggsave("man/figures/readme-spatial-cv.png", p2,
       width = 8.5, height = 7, dpi = 150, bg = "white")

# --- Figure 3: resolution sweep, with the selected k labelled ---------------
# determine_optimal_levels() combines a geometric WSS elbow with Moran's I
# on OLS residuals (val ~ west) at each candidate k; the first candidate is
# the best-ranked cell count.
k_sel <- determine_optimal_levels(pts, max_levels = 15,
                                  response_var = "val",
                                  predictor_vars = "west")[1]

ks <- sort(unique(c(max(2, k_sel - 3), k_sel, 60)))
if (length(ks) < 3) ks <- sort(unique(c(ks, 12)))
labels <- vapply(ks, function(k) {
  if (k == k_sel) sprintf("k = %d  (chosen by determine_optimal_levels)", k)
  else            sprintf("k = %d", k)
}, character(1))

res_cells <- do.call(rbind, Map(function(k, lab) {
  seeds_k <- voronoi_seeds_kmeans(pts, k = k)
  vor_k   <- create_voronoi_polygons(seeds_k, boundary = state, quiet = TRUE)
  cell_means(vor_k$cells, lab)
}, ks, labels))
lvls3 <- c("Raw observations", labels)

ggsave("man/figures/readme-resolution.png",
       panel_plot(res_cells, lvls3),
       width = 11, height = 6, dpi = 150, bg = "white")

message("Wrote man/figures/readme-tessellations.png")
message("Wrote man/figures/readme-spatial-cv.png")
message("Wrote man/figures/readme-resolution.png")
