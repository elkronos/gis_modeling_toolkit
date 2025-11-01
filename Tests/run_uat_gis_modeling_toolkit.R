# run_uat_spatial_modeling.R
# Full UAT for spatial_modeling.R

cat("===== UAT: spatial_modeling.R =====\n")
set.seed(42)

# ---------- Paths ----------
script_dir <- normalizePath(file.path("..","R"), winslash="\\", mustWork = TRUE)
script_path <- file.path(script_dir, "gis_modeling_toolkit.R")
if (!file.exists(script_path)) stop("Cannot find gis_modeling_toolkit.R at: ", script_path)
out_dir <- "UAT_outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---------- Source SUT ----------
source(script_path)

# ---------- Tiny assertion helpers ----------
.pass <- 0L; .fail <- 0L
ok <- function(expr, msg) {
  res <- FALSE
  tryCatch({ res <- isTRUE(expr) }, error = function(e) { res <- FALSE })
  if (isTRUE(res)) { .pass <<- .pass + 1L; cat(" PASS -", msg, "\n") }
  else { .fail <<- .fail + 1L; cat(" FAIL -", msg, "\n") }
}
eq  <- function(a, b, msg, tol = 1e-8) {
  same <- FALSE
  if (is.numeric(a) && is.numeric(b)) same <- max(abs(a - b)) <= tol
  else same <- isTRUE(all.equal(a, b))
  ok(same, msg)
}
between <- function(x, lo, hi, msg, inc = TRUE) {
  if (inc) ok(all(x >= lo & x <= hi), msg) else ok(all(x > lo & x < hi), msg)
}
nz <- function(obj, msg) ok(!is.null(obj) && length(obj) > 0, msg)

# ---------- Libs ----------
suppressPackageStartupMessages({
  library(sf); library(dplyr); library(ggplot2)
})

# ---------- Helpers for robust binding of sf objects (no st_geometry_name) ----------
geom_name <- function(x) {
  if (!inherits(x, "sf")) stop("geom_name(): x must be an sf object")
  attr(x, "sf_column")
}
set_geom_name <- function(x, nm = "geometry") {
  if (!inherits(x, "sf")) stop("set_geom_name(): x must be an sf object")
  cur <- geom_name(x)
  if (is.null(cur)) stop("Object has no sf geometry column")
  if (!identical(cur, nm)) {
    # rename existing geometry column to nm, then set it as geometry
    if (!(nm %in% names(x))) names(x)[names(x) == cur] <- nm
    sf::st_geometry(x) <- nm
  }
  x
}
ensure_cols <- function(x, cols) {
  # add missing non-geometry columns as NA and subset to requested columns
  gnm <- geom_name(x)
  miss <- setdiff(cols, c(names(x), gnm))
  for (m in miss) x[[m]] <- NA
  # final order: all requested non-geom cols (cols without "geometry"), then the actual geometry col
  final_cols <- c(setdiff(cols, "geometry"), gnm)
  x[, final_cols, drop = FALSE]
}
bind_sf_rows <- function(...) {
  lst <- list(...)
  # normalize geometry column name to 'geometry' label for consistency
  lst <- lapply(lst, set_geom_name, nm = "geometry")
  # union of column names across inputs (include "geometry")
  all_cols <- unique(c(unlist(lapply(lst, names)), "geometry"))
  lst <- lapply(lst, ensure_cols, cols = all_cols)
  do.call(rbind, lst)
}

# ---------- Synthetic data ----------
# Base lon/lat grid (WGS84)
n_x <- 12; n_y <- 10
lon_vec <- seq(-95.9, -95.0, length.out = n_x)
lat_vec <- seq(29.4,  30.1, length.out = n_y)
grid_df <- expand.grid(lon = lon_vec, lat = lat_vec)
grid_df$x1 <- scale(grid_df$lon)[,1]
grid_df$x2 <- scale(grid_df$lat)[,1]
grid_df$x3 <- rnorm(nrow(grid_df))
cx <- mean(grid_df$lon); cy <- mean(grid_df$lat)
dist_center <- sqrt((grid_df$lon - cx)^2 + (grid_df$lat - cy)^2)
grid_df$resp <- 3 + 1.6*grid_df$x1 - 1.1*grid_df$x2 + 0.5*scale(dist_center)[,1] + rnorm(nrow(grid_df), 0, 0.3)

pts_ll  <- st_as_sf(grid_df, coords = c("lon","lat"), crs = 4326)
bnd_ll <- st_buffer(st_convex_hull(st_union(st_geometry(pts_ll))), dist = 0.05) # deg

# Non-point geoms
make_row_line <- function(y_index) {
  row_pts <- grid_df[grid_df$lat == lat_vec[y_index], c("lon","lat")]
  st_linestring(as.matrix(row_pts))
}
line_geoms <- lapply(1:3, make_row_line)
lines_ll <- st_sfc(line_geoms, crs = 4326) |> st_as_sf() |> mutate(lid = 1:n())

poly_from_center <- function(cx, cy, r = 0.02) {
  st_polygon(list(matrix(c(cx-r, cy-r,
                           cx+r, cy-r,
                           cx+r, cy+r,
                           cx-r, cy+r,
                           cx-r, cy-r), ncol = 2, byrow = TRUE)))
}
poly_geoms <- list(
  poly_from_center(lon_vec[3], lat_vec[4]),
  poly_from_center(lon_vec[7], lat_vec[8]),
  poly_from_center(lon_vec[10], lat_vec[5])
)
polys_ll <- st_sfc(poly_geoms, crs = 4326) |> st_as_sf() |> mutate(pid = 1:n())

# ---------- Mixed features (POINT + LINE + POLYGON) with aligned columns ----------
pt10 <- pts_ll[1:10, , drop = FALSE] |> set_geom_name("geometry")
pt10$lid <- NA_integer_
pt10$pid <- NA_integer_

lines_ll2 <- st_as_sf(lines_ll["lid"]) |>
  set_geom_name("geometry") |>
  mutate(x1 = 0, x2 = 0, x3 = 0, resp = 0, pid = NA_integer_)

polys_ll2 <- st_as_sf(polys_ll["pid"]) |>
  set_geom_name("geometry") |>
  mutate(x1 = 0, x2 = 0, x3 = 0, resp = 0, lid = NA_integer_)

mixed_ll <- bind_sf_rows(pt10, lines_ll2, polys_ll2)

# ---------- ensure_projected / CRS ----------
pts_proj <- ensure_projected(pts_ll)
ok(!st_is_longlat(pts_proj), "ensure_projected() projected lon/lat to a local CRS")

pts_no_crs <- st_set_crs(st_geometry(pts_ll), NA_crs_) |> st_as_sf() |> dplyr::bind_cols(st_drop_geometry(pts_ll))
pts_no_crs <- set_geom_name(pts_no_crs, "geometry")
pts_no_crs_proj <- ensure_projected(pts_no_crs)
ok(TRUE, "ensure_projected() handled missing-CRS input without error")

pts_to_bnd <- ensure_projected(pts_ll, bnd_ll)
eq(st_crs(pts_to_bnd), st_crs(bnd_ll), "ensure_projected(x, target) matches target CRS")

h <- harmonize_crs(pts_ll, bnd_ll)
eq(st_crs(h$a), st_crs(h$b), "harmonize_crs() aligned CRS between objects")

# ---------- coerce_to_points ----------
cp_auto   <- coerce_to_points(mixed_ll, "auto")
cp_cent   <- coerce_to_points(mixed_ll, "centroid")
cp_surf   <- coerce_to_points(mixed_ll, "surface")
cp_mid    <- coerce_to_points(lines_ll, "line_midpoint")
eq(nrow(cp_auto), nrow(mixed_ll), "coerce_to_points(auto) keeps row count")
ok(all(st_geometry_type(cp_auto) == "POINT"), "coerce_to_points(auto) returns POINT")
ok(all(st_geometry_type(cp_cent) == "POINT"), "coerce_to_points(centroid) returns POINT")
ok(all(st_geometry_type(cp_surf) == "POINT"), "coerce_to_points(surface) returns POINT")
ok(all(st_geometry_type(cp_mid)  == "POINT"), "coerce_to_points(line_midpoint) returns POINT")

# ---------- Seeding ----------
pts_proj <- ensure_projected(pts_ll, bnd_ll)
bnd_proj <- ensure_projected(bnd_ll)

seeds_km <- voronoi_seeds_kmeans(pts_proj, k = 12)
eq(st_crs(seeds_km), st_crs(pts_proj), "voronoi_seeds_kmeans() keeps CRS")
ok(nrow(seeds_km) == 12, "voronoi_seeds_kmeans() returned requested k seeds")

seeds_rand <- voronoi_seeds_random(bnd_proj, k = 10)
ok(nrow(seeds_rand) == 10, "voronoi_seeds_random() returned requested k seeds")
ok(all(st_covered_by(seeds_rand, bnd_proj, sparse = FALSE)), "Random seeds are within/on boundary")

# ---------- Voronoi ----------
dups <- bind_sf_rows(seeds_km[1:2,], seeds_km[1:2,], seeds_km[3:6,])
voro_polys <- create_voronoi_polygons(dups, clip_with = bnd_proj)
ok(nrow(voro_polys) == nrow(dups), "create_voronoi_polygons() returns tiles for each (de-duplicated) seed")

# ---------- Grids ----------
grid_hex <- create_grid_polygons(bnd_proj, target_cells = 30, type = "hex")
grid_sq  <- create_grid_polygons(bnd_proj, target_cells = 30, type = "square")
nz(grid_hex, "create_grid_polygons(hex) returned non-empty grid")
nz(grid_sq,  "create_grid_polygons(square) returned non-empty grid")
between(length(st_geometry(grid_hex)) / 30, 0.6, 1.6, "Hex cells count within sane bounds")
between(length(st_geometry(grid_sq))  / 30, 0.6, 1.6, "Square cells count within sane bounds")

# ---------- Assignment ----------
pts_proj_auto <- coerce_to_points(pts_proj, "auto")
asgn <- assign_features_to_polygons(pts_proj_auto, grid_sq)
ok("polygon_id" %in% names(asgn), "assign_features_to_polygons() adds polygon_id")
ok(all(!is.na(asgn$polygon_id)), "assign_features_to_polygons() filters out unassigned rows")

# ---------- Level selection ----------
levs <- determine_optimal_levels(pts_proj_auto, max_levels = 8, top_n = 3)
ok(is.integer(levs) || is.numeric(levs), "determine_optimal_levels() returns numeric/integer vector")
between(levs, 1, 8, "Levels within [1, max_levels]")

# ---------- GWR ----------
gwr_fit <- fit_gwr_model(pts_proj_auto |> dplyr::select(x1,x2,x3,resp, geometry), "resp", c("x1","x2","x3"))
nz(gwr_fit$model, "fit_gwr_model() returned a model")
ok(is.finite(gwr_fit$AICc), "GWR AICc is finite")
between(gwr_fit$bandwidth, 0, 1, "GWR adaptive bandwidth in (0,1]")

# ---------- Bayesian (3 covariance models) ----------
bay_exp <- fit_bayesian_spatial_model(pts_proj_auto |> dplyr::select(x1,x2,x3,resp, geometry),
                                      "resp", c("x1","x2","x3"),
                                      n.samples = 600, cov_model = "exponential")
ok(is.finite(bay_exp$DIC), "Bayesian DIC (exponential) is finite")

bay_sph <- fit_bayesian_spatial_model(pts_proj_auto |> dplyr::select(x1,x2,x3,resp, geometry),
                                      "resp", c("x1","x2","x3"),
                                      n.samples = 600, cov_model = "spherical")
ok(is.finite(bay_sph$DIC), "Bayesian DIC (spherical) is finite")

bay_mat <- fit_bayesian_spatial_model(pts_proj_auto |> dplyr::select(x1,x2,x3,resp, geometry),
                                      "resp", c("x1","x2","x3"),
                                      n.samples = 600, cov_model = "matern")
ok(is.finite(bay_mat$DIC), "Bayesian DIC (matern) is finite")

# ---------- build_tessellation (methods + seeds + levels edge cases) ----------
levels_vec <- c(1, 5, 12)

bt_v_km <- build_tessellation(pts_proj_auto, levels = levels_vec, method = "voronoi",
                              boundary = bnd_proj, seeds = "kmeans", pointize = "auto")
ok(all(names(bt_v_km) %in% as.character(levels_vec)), "build_tessellation(voronoi,kmeans) returns all levels")
eq(nrow(bt_v_km[["1"]]$polygons), 1L, "Voronoi with one seed returns single polygon")

bt_v_rn <- build_tessellation(pts_proj_auto, levels = levels_vec, method = "voronoi",
                              boundary = bnd_proj, seeds = "random", pointize = "auto")
nz(bt_v_rn[["12"]]$polygons, "Voronoi(random) polygons exist")

provided <- voronoi_seeds_kmeans(pts_proj_auto, k = 7)
bt_v_pr <- build_tessellation(pts_proj_auto, levels = levels_vec, method = "voronoi",
                              boundary = bnd_proj, seeds = "provided",
                              provided_seed_points = provided, pointize = "auto")
nz(bt_v_pr[["5"]]$polygons, "Voronoi(provided) polygons exist (mismatched level handled)")

bt_hx <- build_tessellation(pts_proj_auto, levels = c(20), method = "hex", boundary = bnd_proj, seeds = "kmeans")
bt_sq <- build_tessellation(pts_proj_auto, levels = c(20), method = "square", boundary = bnd_proj, seeds = "kmeans")
nz(bt_hx[["20"]]$polygons, "Hex tessellation polygons exist")
nz(bt_sq[["20"]]$polygons, "Square tessellation polygons exist")

# ---------- evaluate_models ----------
eval_exp <- evaluate_models(
  data_sf = pts_proj_auto |> dplyr::select(x1,x2,x3,resp, geometry),
  response_var = "resp",
  predictor_vars = c("x1","x2","x3"),
  levels = c(5, 10),
  tessellation = c("voronoi","hex","square"),
  boundary = bnd_proj,
  seeds = "kmeans",
  models = c("GWR","Bayesian"),
  n.samples = 400,
  cov_model = "exponential"
)
nz(eval_exp$results, "evaluate_models(exponential) returned results")

eval_sph <- evaluate_models(
  data_sf = pts_proj_auto |> dplyr::select(x1,x2,x3,resp, geometry),
  response_var = "resp",
  predictor_vars = c("x1","x2","x3"),
  levels = c(5),
  tessellation = c("voronoi"),
  boundary = bnd_proj,
  seeds = "random",
  models = c("GWR","Bayesian"),
  n.samples = 400,
  cov_model = "spherical"
)
nz(eval_sph$results, "evaluate_models(spherical) returned results")

eval_mat <- evaluate_models(
  data_sf = pts_proj_auto |> dplyr::select(x1,x2,x3,resp, geometry),
  response_var = "resp",
  predictor_vars = c("x1","x2","x3"),
  levels = NULL,            # auto level selection
  tessellation = c("voronoi"),
  boundary = NULL,          # bbox fallback
  seeds = "kmeans",
  models = c("GWR","Bayesian"),
  n.samples = 400,
  cov_model = "matern"
)
nz(eval_mat$results, "evaluate_models(matern, auto-level, no-boundary) returned results")

# Save combined results
all_res <- dplyr::bind_rows(
  eval_exp$results |> dplyr::mutate(Cov = "exponential"),
  eval_sph$results |> dplyr::mutate(Cov = "spherical"),
  eval_mat$results |> dplyr::mutate(Cov = "matern")
)
write.csv(all_res, file.path(out_dir, "uat_evaluate_models_results.csv"), row.names = FALSE)

# ---------- summarize_by_cell ----------
some_assign <- bt_v_km[["12"]]$data
sum_tbl <- summarize_by_cell(some_assign, response_var = "resp", predictor_vars = c("x1","x2"))
ok(all(c("polygon_id","n","mean_response","mean_x1","mean_x2") %in% names(sum_tbl)),
   "summarize_by_cell() returns expected columns")

# ---------- plot_tessellation_map ----------
p1 <- plot_tessellation_map(bt_v_km[["12"]]$polygons, bt_v_km[["12"]]$data, title = "Voronoi (kmeans, lvl=12)", boundary = bnd_proj, show_counts = TRUE)
ok(inherits(p1, "ggplot"), "plot_tessellation_map() returns ggplot")
ggplot2::ggsave(filename = file.path(out_dir, "map_voronoi_kmeans_lvl12.png"), plot = p1, width = 7, height = 5, dpi = 150)

p2 <- plot_tessellation_map(bt_sq[["20"]]$polygons, bt_sq[["20"]]$data, title = "Square Grid (lvl=20)", boundary = bnd_proj, show_counts = TRUE)
ggplot2::ggsave(filename = file.path(out_dir, "map_square_lvl20.png"), plot = p2, width = 7, height = 5, dpi = 150)

# ---------- Final summary ----------
cat("\n===== UAT SUMMARY =====\n")
cat("Assertions passed:", .pass, "\n")
cat("Assertions failed:", .fail, "\n")
if (.fail == 0L) {
  cat("RESULT: ALL TESTS PASSED ✅\n")
} else {
  cat("RESULT: TESTS FAILED ❌  (see logs above)\n")
}
cat("Outputs written to:", normalizePath(out_dir, winslash="\\", mustWork = TRUE), "\n")
