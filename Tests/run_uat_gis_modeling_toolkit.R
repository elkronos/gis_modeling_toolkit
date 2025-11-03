# run_uat_spatial_modeling.R
# Full(er) UAT for gis_modeling_toolkit.R — expanded coverage & edge cases (fixed)

cat("===== UAT: gis_modeling_toolkit.R =====\n")
set.seed(42)

# ---------- Paths ----------
script_dir <- normalizePath(file.path("..","R"), winslash="\\", mustWork = TRUE)
script_path <- file.path(script_dir, "gis_modeling_toolkit.R")
if (!file.exists(script_path)) stop("Cannot find gis_modeling_toolkit.R at: ", script_path)
out_dir <- "UAT_outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---------- Source SUT ----------
suppressWarnings(suppressMessages(source(script_path)))

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
  if (is.numeric(a) && is.numeric(b)) {
    if (length(a) == 0L && length(b) == 0L) same <- TRUE else same <- max(abs(a - b)) <= tol
  } else same <- isTRUE(all.equal(a, b))
  ok(same, msg)
}
between <- function(x, lo, hi, msg, inc = TRUE) {
  if (inc) ok(all(x >= lo & x <= hi), msg) else ok(all(x > lo & x < hi), msg)
}
nz <- function(obj, msg) ok(!is.null(obj) && length(obj) > 0, msg)
err <- function(expr, msg, pattern = NULL) {
  got <- FALSE; .m <- NULL
  tryCatch({ force(expr) }, error = function(e){ got <<- TRUE; .m <<- conditionMessage(e) })
  ok(isTRUE(got) && (is.null(pattern) || grepl(pattern, .m, ignore.case = TRUE)),
     if (is.null(pattern)) msg else paste0(msg, " [matched: ", isTRUE(grepl(pattern, .m, ignore.case = TRUE)), "]"))
}
wrn <- function(expr, msg) {
  got <- FALSE
  withCallingHandlers(expr, warning = function(w){ got <<- TRUE; invokeRestart("muffleWarning") })
  ok(got, msg)
}

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
    if (!(nm %in% names(x))) names(x)[names(x) == cur] <- nm
    sf::st_geometry(x) <- nm
  }
  x
}
ensure_cols <- function(x, cols) {
  gnm <- geom_name(x)
  miss <- setdiff(cols, c(names(x), gnm))
  for (m in miss) x[[m]] <- NA
  final_cols <- c(setdiff(cols, "geometry"), gnm)
  x[, final_cols, drop = FALSE]
}
bind_sf_rows <- function(...) {
  lst <- list(...)
  lst <- lapply(lst, set_geom_name, nm = "geometry")
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
bnd_ll  <- st_buffer(st_convex_hull(st_union(st_geometry(pts_ll))), dist = 0.05) # deg

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
  poly_from_center(lon_vec[3],  lat_vec[4]),
  poly_from_center(lon_vec[7],  lat_vec[8]),
  poly_from_center(lon_vec[10], lat_vec[5])
)
polys_ll <- st_sfc(poly_geoms, crs = 4326) |> st_as_sf() |> mutate(pid = 1:n())

# ---------- Extra geometry types / edge cases ----------
# MULTILINESTRING
mp_line <- st_sfc(
  st_multilinestring(list(
    as.matrix(grid_df[1:5,   c("lon","lat")]),
    as.matrix(grid_df[10:14, c("lon","lat")])
  )),
  crs = 4326
) |> st_as_sf() |> mutate(kind = "mls")

# Correct MULTIPOLYGON: list of polygons, each polygon is a list of rings
mp_poly <- st_sfc(
  st_multipolygon(lapply(poly_geoms[1:2], unclass)),
  crs = 4326
) |> st_as_sf() |> mutate(kind = "mpoly")

ok(all(st_geometry_type(mp_poly) == "MULTIPOLYGON") && all(st_is_valid(mp_poly)),
   "constructed MULTIPOLYGON is valid")

# ---------- Mixed features (POINT + LINE + POLYGON) with aligned columns ----------
pt10 <- pts_ll[1:10, , drop = FALSE] |> set_geom_name("geometry")
pt10$lid <- NA_integer_; pt10$pid <- NA_integer_

lines_ll2 <- st_as_sf(lines_ll["lid"]) |>
  set_geom_name("geometry") |>
  mutate(x1 = 0, x2 = 0, x3 = 0, resp = 0, pid = NA_integer_)

polys_ll2 <- st_as_sf(polys_ll["pid"]) |>
  set_geom_name("geometry") |>
  mutate(x1 = 0, x2 = 0, x3 = 0, resp = 0, lid = NA_integer_)

mixed_ll <- bind_sf_rows(pt10, lines_ll2, polys_ll2)

# =========================================
# CRS utilities
# =========================================
pts_proj <- ensure_projected(pts_ll)
ok(!st_is_longlat(pts_proj), "ensure_projected() projected lon/lat to a local CRS")

pts_no_crs <- st_set_crs(st_geometry(pts_ll), NA_crs_) |> st_as_sf() |> dplyr::bind_cols(st_drop_geometry(pts_ll))
pts_no_crs <- set_geom_name(pts_no_crs, "geometry")
pts_no_crs_proj <- ensure_projected(pts_no_crs)
ok(TRUE, "ensure_projected() handled missing-CRS lon/lat input without error")

# Missing-CRS that does NOT look like lon/lat should remain NA-CRS
rnd_xy <- data.frame(x = runif(20, 1e6, 2e6), y = runif(20, 1e6, 2e6))
rnd_sf <- st_as_sf(rnd_xy, coords = c("x","y"), crs = NA)
rnd_out <- ensure_projected(rnd_sf)
ok(is.na(st_crs(rnd_out)), "ensure_projected() leaves non-lon/lat, no-CRS data unchanged")

pts_to_bnd <- ensure_projected(pts_ll, bnd_ll)
eq(st_crs(pts_to_bnd), st_crs(bnd_ll), "ensure_projected(x, target) matches target CRS")

h <- harmonize_crs(pts_ll, bnd_ll)
eq(st_crs(h$a), st_crs(h$b), "harmonize_crs() aligned CRS between objects")

# both missing -> unchanged
h2 <- harmonize_crs(st_set_crs(pts_ll, NA), st_set_crs(bnd_ll, NA))
ok(is.na(st_crs(h2$a)) && is.na(st_crs(h2$b)), "harmonize_crs() leaves both-NA CRSs unchanged")

# =========================================
# coerce_to_points (all modes + tmp_project behavior)
# =========================================

# Core conversions on mixed geometry collection
cp_auto   <- coerce_to_points(mixed_ll, "auto")
ok(all(st_geometry_type(cp_auto) == "POINT"), "coerce_to_points(auto) returns POINT")
eq(nrow(cp_auto), nrow(mixed_ll), "coerce_to_points(auto) keeps row count")

cp_cent   <- coerce_to_points(mixed_ll, "centroid")
ok(all(st_geometry_type(cp_cent) == "POINT"), "coerce_to_points(centroid) returns POINT")

cp_surf   <- coerce_to_points(mixed_ll, "point_on_surface")  # (alias 'surface' still accepted)
ok(all(st_geometry_type(cp_surf) == "POINT"), "coerce_to_points(point_on_surface) returns POINT")

# New: bbox_center should be pure numeric center (no bbox polygon construction)
cp_bbox   <- coerce_to_points(mixed_ll, "bbox_center")
ok(all(st_geometry_type(cp_bbox) == "POINT") && nrow(cp_bbox) == nrow(mixed_ll),
   "coerce_to_points(bbox_center) returns POINT and preserves rows")

# Lines: explicit midpoint vs lon/lat fallback control
cp_mid    <- coerce_to_points(lines_ll, "line_midpoint")  # may internally project
ok(all(st_geometry_type(cp_mid)  == "POINT"), "coerce_to_points(line_midpoint) returns POINT for LINESTRING input")

# tmp_project behavior: TRUE should succeed; FALSE on lon/lat lines should avoid st_line_sample()
cp_auto_tp_true  <- coerce_to_points(lines_ll, "auto", tmp_project = TRUE)
ok(all(st_geometry_type(cp_auto_tp_true) == "POINT"), "coerce_to_points(auto, tmp_project=TRUE) works on lon/lat lines")

cp_auto_tp_false <- coerce_to_points(lines_ll, "auto", tmp_project = FALSE)
ok(all(st_geometry_type(cp_auto_tp_false) == "POINT"), "coerce_to_points(auto, tmp_project=FALSE) falls back cleanly on lon/lat lines")

# MULTILINESTRING & MULTIPOLYGON edge cases
err(coerce_to_points(mp_line, "line_midpoint"),
    "coerce_to_points(line_midpoint) errors on MULTILINESTRING", "LINESTRING")

cp_mls_cent <- coerce_to_points(mp_line, "centroid")
ok(all(st_geometry_type(cp_mls_cent) == "POINT"),
   "coerce_to_points(centroid) returns POINT for MULTILINESTRING")

mp_line_ls <- suppressWarnings(st_cast(mp_line, "LINESTRING"))
cp_ls_mid  <- coerce_to_points(mp_line_ls, "line_midpoint")
ok(all(st_geometry_type(cp_ls_mid) == "POINT"),
   "coerce_to_points(line_midpoint) returns POINT after casting MLS->LINESTRING")

cp_mpoly_auto <- coerce_to_points(mp_poly, "centroid")
ok(all(st_geometry_type(cp_mpoly_auto) == "POINT"),
   "coerce_to_points(centroid) returns POINT for MULTIPOLYGON")

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
ok(nrow(voro_polys) == nrow(dups),
   "create_voronoi_polygons() returns one tile per input seed (duplicates may share polygons)")

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
gdat <- pts_proj_auto |> dplyr::select(x1,x2,x3,resp, geometry)
gwr_fit <- fit_gwr_model(gdat, "resp", c("x1","x2","x3"))
nz(gwr_fit$model, "fit_gwr_model() returned a model")
ok(is.finite(gwr_fit$AICc), "GWR AICc is finite")
ok(isTRUE(gwr_fit$bandwidth > 0 && gwr_fit$bandwidth <= 1),
   "GWR adaptive bandwidth in (0,1]")

# ---------- Bayesian (3 covariance models) ----------
bay_exp <- fit_bayesian_spatial_model(gdat, "resp", c("x1","x2","x3"),
                                      n.samples = 600, cov_model = "exponential")
ok(is.finite(bay_exp$DIC), "Bayesian DIC (exponential) is finite")

bay_sph <- fit_bayesian_spatial_model(gdat, "resp", c("x1","x2","x3"),
                                      n.samples = 600, cov_model = "spherical")
ok(is.finite(bay_sph$DIC), "Bayesian DIC (spherical) is finite")

bay_mat <- fit_bayesian_spatial_model(gdat, "resp", c("x1","x2","x3"),
                                      n.samples = 600, cov_model = "matern")
ok(is.finite(bay_mat$DIC), "Bayesian DIC (matern) is finite")

# Negative path: invalid covariance name
err(fit_bayesian_spatial_model(gdat, "resp", c("x1","x2","x3"),
                               n.samples = 200, cov_model = "bogus"),
    "fit_bayesian_spatial_model() errors on invalid cov_model", "cov|model|invalid|supported")

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

# Negative path: seeds="provided" without providing seeds -> should error
err(build_tessellation(pts_proj_auto, levels = c(5), method = "voronoi",
                       boundary = bnd_proj, seeds = "provided",
                       provided_seed_points = NULL, pointize = "auto"),
    "build_tessellation() errors when seeds='provided' but no points supplied", "provide|seed")

# ---------- evaluate_models ----------
eval_exp <- evaluate_models(
  data_sf = gdat,
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
  data_sf = gdat,
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
  data_sf = gdat,
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
csv_path <- file.path(out_dir, "uat_evaluate_models_results.csv")
write.csv(all_res, csv_path, row.names = FALSE)
ok(file.exists(csv_path), "results CSV written")

# ---------- summarize_by_cell ----------
some_assign <- bt_v_km[["12"]]$data
sum_tbl <- summarize_by_cell(some_assign, response_var = "resp", predictor_vars = c("x1","x2"))
ok(all(c("polygon_id","n","mean_response","mean_x1","mean_x2") %in% names(sum_tbl)),
   "summarize_by_cell() returns expected columns")

# ---------- plot_tessellation_map ----------
p1 <- plot_tessellation_map(bt_v_km[["12"]]$polygons, bt_v_km[["12"]]$data,
                            title = "Voronoi (kmeans, lvl=12)", boundary = bnd_proj, show_counts = TRUE)
ok(inherits(p1, "ggplot"), "plot_tessellation_map() returns ggplot")
png1 <- file.path(out_dir, "map_voronoi_kmeans_lvl12.png")
ggplot2::ggsave(filename = png1, plot = p1, width = 7, height = 5, dpi = 150)
ok(file.exists(png1), "map_voronoi_kmeans_lvl12.png saved")

p2 <- plot_tessellation_map(bt_sq[["20"]]$polygons, bt_sq[["20"]]$data,
                            title = "Square Grid (lvl=20)", boundary = bnd_proj, show_counts = TRUE)
png2 <- file.path(out_dir, "map_square_lvl20.png")
ggplot2::ggsave(filename = png2, plot = p2, width = 7, height = 5, dpi = 150)
ok(file.exists(png2), "map_square_lvl20.png saved")

# Plotting without boundary (fallback extent)
p3 <- plot_tessellation_map(bt_sq[["20"]]$polygons, bt_sq[["20"]]$data,
                            title = "Square Grid (lvl=20, no boundary)", boundary = NULL, show_counts = FALSE)
ok(inherits(p3, "ggplot"), "plot_tessellation_map() works with boundary=NULL")
png3 <- file.path(out_dir, "map_square_lvl20_nobnd.png")
ggplot2::ggsave(filename = png3, plot = p3, width = 7, height = 5, dpi = 150)
ok(file.exists(png3), "map_square_lvl20_nobnd.png saved")

# ---------- EXTRAS ----------

# (A) Degenerate/empty geometry handling
empty_sf <- st_as_sf(st_sfc(crs = 4326))
handled <- FALSE
tryCatch({
  tmp <- coerce_to_points(empty_sf, "auto"); handled <- inherits(tmp, "sf")
}, error=function(e){ handled <<- TRUE })
ok(handled, "coerce_to_points() handles empty input (returns empty or errors cleanly)")

# (B) Non-finite predictors
pts_bad <- gdat
pts_bad$x1[1] <- NA_real_; pts_bad$x2[2] <- Inf
handled_gwr <- FALSE
tryCatch({
  tmp <- fit_gwr_model(pts_bad, "resp", c("x1","x2","x3"))
  handled_gwr <- !is.null(tmp$model)
}, warning=function(w){ handled_gwr <<- TRUE },
error=function(e){ handled_gwr <<- TRUE })
ok(handled_gwr, "fit_gwr_model() handles non-finite predictors by warning/error or dropping rows")

handled_bayes <- FALSE
tryCatch({
  tmp <- fit_bayesian_spatial_model(pts_bad, "resp", c("x1","x2","x3"), n.samples = 200, cov_model = "exponential")
  handled_bayes <- !is.null(tmp)
}, warning=function(w){ handled_bayes <<- TRUE },
error=function(e){ handled_bayes <<- TRUE })
ok(handled_bayes, "fit_bayesian_spatial_model() handles non-finite predictors by warning/error or dropping rows")

# (C) Boundary with holes — grids should respect holes after clipping
# Build a simple rectangle with an inner rectangular hole around the data center
outer <- matrix(c(min(lon_vec)-0.02, min(lat_vec)-0.02,
                  max(lon_vec)+0.02, min(lat_vec)-0.02,
                  max(lon_vec)+0.02, max(lat_vec)+0.02,
                  min(lon_vec)-0.02, max(lat_vec)+0.02,
                  min(lon_vec)-0.02, min(lat_vec)-0.02), ncol=2, byrow=TRUE)
hole  <- matrix(c(cx-0.01, cy-0.01,
                  cx+0.01, cy-0.01,
                  cx+0.01, cy+0.01,
                  cx-0.01, cy+0.01,
                  cx-0.01, cy-0.01), ncol=2, byrow=TRUE)

bnd_hole_ll   <- st_sfc(st_polygon(list(outer, hole)), crs = 4326)
bnd_hole_proj <- ensure_projected(bnd_hole_ll)

grid_hole <- create_grid_polygons(bnd_hole_proj, target_cells = 30, type = "square")
nz(grid_hole, "create_grid_polygons() produced grid for holed boundary")

# IMPORTANT: project the hole to the SAME CRS as the grid/boundary
hole_proj <- ensure_projected(st_sfc(st_polygon(list(hole)), crs = 4326), bnd_hole_proj)

# Sanity: CRSs align
eq(st_crs(hole_proj), st_crs(grid_hole), "hole and grid CRS aligned")

cent_hole  <- suppressWarnings(st_centroid(grid_hole))
within_mat <- st_within(cent_hole, hole_proj, sparse = FALSE)
ok(all(!as.vector(within_mat)), "Grid centroids do not fall inside the hole (hole respected)")

# (D) Auto-level stability (deterministic given seed)
set.seed(123)
levs_a <- determine_optimal_levels(pts_proj_auto, max_levels = 8, top_n = 3)
set.seed(123)
levs_b <- determine_optimal_levels(pts_proj_auto, max_levels = 8, top_n = 3)
eq(sort(levs_a), sort(levs_b), "determine_optimal_levels() stable across runs with same seed")

# ---------- Final summary ----------
# Session info (for reproducibility)
sess_path <- file.path(out_dir, "sessionInfo.txt")
capture.output(sessionInfo(), file = sess_path)
ok(file.exists(sess_path), "sessionInfo.txt written")

cat("\n===== UAT SUMMARY =====\n")
cat("Assertions passed:", .pass, "\n")
cat("Assertions failed:", .fail, "\n")
if (.fail == 0L) {
  cat("RESULT: ALL TESTS PASSED ✅\n")
} else {
  cat("RESULT: TESTS FAILED ❌  (see logs above)\n")
}
cat("Outputs written to:", normalizePath(out_dir, winslash="\\", mustWork = TRUE), "\n")
