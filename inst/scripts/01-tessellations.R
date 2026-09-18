# =============================================================================
# 01  Four ways to cut a study area into cells
# =============================================================================
#   source(system.file("scripts", "01-tessellations.R", package = "spatialkit"))
# =============================================================================
# Locate the shared setup whether the package is installed or you are sitting
# in a source checkout.
.tour_dir <- if (nzchar(system.file("scripts", package = "spatialkit")))
  system.file("scripts", package = "spatialkit") else "."
source(file.path(.tour_dir, "_common.R"))

pts <- tour_points()
bnd <- tour_boundary(pts)

step("01.1", "Voronoi cells grown from k-means seeds")
# Voronoi cells grow from SEEDS, so place the seeds first and tessellate those.
# Hand build_tessellation() the observations instead and you get one cell per
# observation, which is section 01.5.
seeds <- get_voronoi_seeds(boundary = bnd, method = "kmeans", n = 25,
                           sample_points = pts, set_seed = 1)
vor   <- build_tessellation(seeds, boundary = bnd, method = "voronoi",
                            quiet = TRUE)
stopifnot(nrow(vor$cells) == nrow(seeds))
cat(sprintf("  %d seeds -> %d cells\n", nrow(seeds), nrow(vor$cells)))

step("01.2", "Hex and square lattices (these need a boundary)")
hex <- build_tessellation(pts, boundary = bnd, method = "hex",
                          approx_n_cells = 25, quiet = TRUE)
sqr <- build_tessellation(pts, boundary = bnd, method = "square",
                          approx_n_cells = 25, quiet = TRUE)
cat(sprintf("  hex %d cells, square %d cells\n", nrow(hex$cells), nrow(sqr$cells)))

step("01.3", "Delaunay triangles")
tri <- if (skip_without("geometry", "Delaunay triangles")) NULL else
  build_tessellation(pts, boundary = bnd, method = "triangles", quiet = TRUE)
if (!is.null(tri)) cat(sprintf("  %d triangles\n", nrow(tri$cells)))

step("01.4", "Draw them")
notes <- c(
  voronoi   = paste("cells are small where samples are dense and large where",
                    "they are sparse, so the counts per cell even out."),
  hex       = paste("every cell is the same size, so a thinly sampled corner",
                    "gets a cell with almost nothing in it."),
  square    = paste("same trade as hex, plus the cells line up with the axes,",
                    "which can align with whatever laid out the sampling."),
  triangles = paste("one triangle per triple of neighbouring points, so there",
                    "are far more cells than you would choose on purpose."))
if (!skip_without("ggplot2", "the maps")) {
  for (nm in names(Filter(Negate(is.null),
                          list(voronoi = vor, hex = hex, square = sqr, triangles = tri)))) {
    t_ <- list(voronoi = vor, hex = hex, square = sqr, triangles = tri)[[nm]]
    asg <- assign_features_to_polygons(pts, t_$cells)
    cel <- summarize_by_cell(asg, "z", cells_sf = t_$cells, deff = 1)
    look_for(sprintf("%s -- %s", nm, notes[[nm]]))
    show_plot(plot_tessellation_map(cel, boundary = bnd, features_sf = pts,
                                    fill_col = "resp_mean_z",
                                    title = sprintf("%s: mean of z per cell", nm)),
              sprintf("01-%s.png", nm))
  }
}

step("01.5", "Voronoi over the observations themselves")
# `approx_n_cells` sizes the grid methods, and only those. Ask for it with
# method = "voronoi" and the package says so -- the WARN line below is the
# package talking, not this script -- before handing back one cell per
# observation: a nearest-neighbour interpolation, where every cell holds a
# single point, there is no within-cell variation, and every standard error
# comes back NA.
trap <- build_tessellation(pts, boundary = bnd, method = "voronoi",
                           approx_n_cells = 25, quiet = TRUE)
asg1 <- assign_features_to_polygons(pts, trap$cells)
cel1 <- summarize_by_cell(asg1, "z", cells_sf = trap$cells, deff = 1)
cat(sprintf("  asked for 25 cells, got %d; %d hold one point or none, %d SEs are NA\n",
            nrow(cel1), sum(cel1$n <= 1, na.rm = TRUE),
            sum(is.na(cel1$..se_resp_z))))
# The warning goes to the console and nowhere else: params$approx_n_cells is
# still NULL, so a result saved now and reopened next month carries no sign
# that 25 was ever asked for.
cat(sprintf("  the console warning is the only record: params$approx_n_cells is %s\n",
            if (is.null(trap$params$approx_n_cells)) "NULL" else
              trap$params$approx_n_cells))
cat("  Seed first (01.1), or use a grid method, which is what approx_n_cells sizes.\n")
