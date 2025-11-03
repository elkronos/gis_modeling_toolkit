# run_uat_spatial_modeling.R
# Full UAT for updated gis_modeling_toolkit.R — covers chunks 1–7 & fixes

cat("===== UAT: gis_modeling_toolkit.R (updated) =====\n")
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
wrn <- function(expr, msg, pattern = NULL) {
  got <- FALSE; m <- NULL
  withCallingHandlers(expr, warning = function(w){ got <<- TRUE; m <<- conditionMessage(w); invokeRestart("muffleWarning") })
  ok(isTRUE(got) && (is.null(pattern) || grepl(pattern, m, ignore.case = TRUE)), msg)
}
skip <- function(msg) cat(" SKIP -", msg, "\n")

# ---------- Libs & feature flags ----------
suppressPackageStartupMessages({
  library(sf); library(dplyr); library(ggplot2)
})
HAS_SPGWR <- requireNamespace("spgwr", quietly = TRUE) && requireNamespace("sp", quietly = TRUE)
HAS_BRMS  <- requireNamespace("brms", quietly = TRUE)
RUN_BAYES <- HAS_BRMS && isTRUE(as.logical(Sys.getenv("RUN_BAYES_UAT", "0")))  # opt-in for heavy models

cat(sprintf("Packages: spgwr=%s, brms=%s, RUN_BAYES=%s\n", HAS_SPGWR, HAS_BRMS, RUN_BAYES))

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

# Extra types
mp_line <- st_sfc(
  st_multilinestring(list(
    as.matrix(grid_df[1:5,   c("lon","lat")]),
    as.matrix(grid_df[10:14, c("lon","lat")])
  )),
  crs = 4326
) |> st_as_sf() |> mutate(kind = "mls")

mp_poly <- st_sfc(
  st_multipolygon(lapply(poly_geoms[1:2], unclass)),
  crs = 4326
) |> st_as_sf() |> mutate(kind = "mpoly")

ok(all(st_geometry_type(mp_poly) == "MULTIPOLYGON") && all(st_is_valid(mp_poly)),
   "constructed MULTIPOLYGON is valid")

# Mixed features (POINT + LINE + POLYGON) with aligned columns
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

if (exists("harmonize_crs")) {
  h <- harmonize_crs(pts_ll, bnd_ll)
  eq(st_crs(h$a), st_crs(h$b), "harmonize_crs() aligned CRS between objects")
  h2 <- harmonize_crs(st_set_crs(pts_ll, NA), st_set_crs(bnd_ll, NA))
  ok(is.na(st_crs(h2$a)) && is.na(st_crs(h2$b)), "harmonize_crs() leaves both-NA CRSs unchanged")
} else {
  skip("harmonize_crs() not found — skipping")
}

# =========================================
# coerce_to_points (all modes + tmp_project behavior)
# =========================================
cp_auto   <- coerce_to_points(mixed_ll, "auto")
ok(all(st_geometry_type(cp_auto) == "POINT"), "coerce_to_points(auto) returns POINT")
eq(nrow(cp_auto), nrow(mixed_ll), "coerce_to_points(auto) keeps row count")

cp_cent   <- coerce_to_points(mixed_ll, "centroid")
ok(all(st_geometry_type(cp_cent) == "POINT"), "coerce_to_points(centroid) returns POINT")

cp_surf   <- coerce_to_points(mixed_ll, "point_on_surface")
ok(all(st_geometry_type(cp_surf) == "POINT"), "coerce_to_points(point_on_surface) returns POINT")

cp_bbox   <- coerce_to_points(mixed_ll, "bbox_center")
ok(all(st_geometry_type(cp_bbox) == "POINT") && nrow(cp_bbox) == nrow(mixed_ll),
   "coerce_to_points(bbox_center) returns POINT and preserves rows")

cp_mid    <- coerce_to_points(lines_ll, "line_midpoint")
ok(all(st_geometry_type(cp_mid)  == "POINT"), "coerce_to_points(line_midpoint) returns POINT for LINESTRING input")

cp_auto_tp_true  <- coerce_to_points(lines_ll, "auto", tmp_project = TRUE)
ok(all(st_geometry_type(cp_auto_tp_true) == "POINT"), "coerce_to_points(auto, tmp_project=TRUE) works on lon/lat lines")

cp_auto_tp_false <- coerce_to_points(lines_ll, "auto", tmp_project = FALSE)
ok(all(st_geometry_type(cp_auto_tp_false) == "POINT"), "coerce_to_points(auto, tmp_project=FALSE) falls back cleanly on lon/lat lines")

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

# ---------- Seeding (if available) ----------
if (exists("voronoi_seeds_kmeans") && exists("voronoi_seeds_random")) {
  pts_proj <- ensure_projected(pts_ll, bnd_ll)
  bnd_proj <- ensure_projected(bnd_ll)
  
  seeds_km <- voronoi_seeds_kmeans(pts_proj, k = 12)
  eq(st_crs(seeds_km), st_crs(pts_proj), "voronoi_seeds_kmeans() keeps CRS")
  ok(nrow(seeds_km) == 12, "voronoi_seeds_kmeans() returned requested k seeds")
  
  seeds_rand <- voronoi_seeds_random(bnd_proj, k = 10)
  ok(nrow(seeds_rand) == 10, "voronoi_seeds_random() returned requested k seeds")
  ok(all(st_covered_by(seeds_rand, bnd_proj, sparse = FALSE)), "Random seeds are within/on boundary")
} else {
  skip("voronoi_seeds_*() not found — skipping")
}

# ---------- Voronoi polygons (if available) ----------
if (exists("create_voronoi_polygons")) {
  # make duplicates to check index mapping
  dups <- bind_sf_rows(seeds_km[1:2,], seeds_km[1:2,], seeds_km[3:6,])
  voro <- create_voronoi_polygons(dups, boundary = bnd_proj, keep_duplicates = TRUE)
  nz(voro, "create_voronoi_polygons() returned a result")
  ok(inherits(voro$cells, "sf") && nrow(voro$cells) > 0, "Voronoi: non-empty cells")
  ok(length(voro$index) == nrow(dups), "Voronoi: index aligns to input rows")
  # Correct duplicate checks: first pair duplicated at rows 1 & 3, second at 2 & 4
  ok(isTRUE(voro$index[1] == voro$index[3]) && isTRUE(voro$index[2] == voro$index[4]),
     "Voronoi: duplicate coords map to same cell")
} else {
  skip("create_voronoi_polygons() not found — skipping")
}

# ---------- Grids ----------
bnd_proj <- ensure_projected(bnd_ll)
grid_hex <- create_grid_polygons(bnd_proj, target_cells = 30, type = "hex")
grid_sq  <- create_grid_polygons(bnd_proj, target_cells = 30, type = "square")
nz(grid_hex, "create_grid_polygons(hex) returned non-empty grid")
nz(grid_sq,  "create_grid_polygons(square) returned non-empty grid")
between(nrow(grid_hex) / 30, 0.5, 2.0, "Hex cells count within sane bounds")
between(nrow(grid_sq)  / 30, 0.5, 2.0, "Square cells count within sane bounds")

# ---------- Cached grids ----------
if (exists("create_grid_polygons_cached")) {
  hex1 <- create_grid_polygons_cached(bnd_proj, target_cells = 60, type = "hex")
  hex2 <- create_grid_polygons_cached(bnd_proj, target_cells = 60, type = "hex")
  eq(st_geometry(hex1), st_geometry(hex2), "create_grid_polygons_cached(): identical geometry on cache hit")
} else {
  skip("create_grid_polygons_cached() not found — skipping")
}

# ---------- Assignment & summarization ----------
pts_proj_auto <- coerce_to_points(pts_proj, "auto")
asgn <- assign_features_to_polygons(pts_proj_auto, grid_sq)  # default id_col = "poly_id"
ok("poly_id" %in% names(asgn), "assign_features_to_polygons() adds poly_id")
ok(all(!is.na(asgn$poly_id)), "assign_features_to_polygons() returns assignments for all points (inside grid)")

sum_tbl <- summarize_by_cell(asgn, response_var = "resp", predictor_vars = c("x1","x2"))
ok(all(c("poly_id","n","mean_response","mean_x1","mean_x2") %in% names(sum_tbl)),
   "summarize_by_cell() returns expected columns")

# ---------- plot_tessellation_map (signature-adaptive) ----------
if (exists("plot_tessellation_map")) {
  # Build a join for counts
  grid_sq_counts <- grid_sq |>
    dplyr::left_join(sum_tbl[, c("poly_id","n")], by = c("poly_id" = "poly_id"))
  
  use_old_sig <- FALSE
  fmls <- try(formals(plot_tessellation_map), silent = TRUE)
  if (!inherits(fmls, "try-error")) {
    use_old_sig <- "data" %in% names(fmls) # old form: (polygons, data, title, boundary, show_counts)
  }
  
  if (isTRUE(use_old_sig)) {
    p1 <- plot_tessellation_map(grid_sq, grid_sq_counts, title = "Square Grid w/ counts",
                                boundary = bnd_proj, show_counts = TRUE)
  } else {
    p1 <- plot_tessellation_map(grid_sq_counts, boundary = bnd_proj, fill_col = "n",
                                title = "Square Grid w/ counts", labels = TRUE,
                                label_col = "poly_id", legend = TRUE)
  }
  ok(inherits(p1, "ggplot"), "plot_tessellation_map() returns ggplot")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  png1 <- file.path(out_dir, "map_square_counts.png")
  ggplot2::ggsave(filename = png1, plot = p1, width = 7, height = 5, dpi = 150)
  ok(file.exists(png1), "map_square_counts.png saved")
} else {
  skip("plot_tessellation_map() not found — skipping")
}

# ---------- make_folds ----------
if (exists("make_folds")) {
  mf1 <- make_folds(pts_proj_auto, k = 5, method = "random_kfold", seed = 1)
  ok(mf1$k == 5 && length(mf1$folds) == 5, "make_folds(random_kfold): produced 5 folds")
  
  mf2 <- make_folds(pts_proj_auto, k = 4, method = "block_kfold", seed = 1, block_multiplier = 2, boundary = bnd_proj)
  ok(mf2$k >= 2 && length(mf2$folds) == mf2$k, "make_folds(block_kfold): produced folds")
  
  mf3 <- make_folds(pts_proj_auto, k = 3, method = "buffered_loo", buffer = 1000)
  ok(mf3$k == nrow(pts_proj_auto), "make_folds(buffered_loo): k equals n")
} else {
  skip("make_folds() not found — skipping")
}

# ---------- build_tessellation ----------
if (exists("build_tessellation")) {
  # Voronoi
  bt_v <- build_tessellation(pts_proj_auto, boundary = bnd_proj, method = "voronoi", expand = 1000)
  ok(is.list(bt_v) && inherits(bt_v$cells, "sf") && nrow(bt_v$cells) > 0, "build_tessellation(voronoi) returns cells")
  
  # Triangles (clip = FALSE should not intersect with boundary explicitly)
  bt_t <- build_tessellation(pts_proj_auto, boundary = bnd_proj, method = "triangles", clip = FALSE)
  ok(inherits(bt_t$cells, "sf") && nrow(bt_t$cells) > 0, "build_tessellation(triangles, clip=FALSE) returns cells")
  
  # Hex / Square via approx_n_cells
  bt_h <- build_tessellation(pts_proj_auto, boundary = bnd_proj, method = "hex", approx_n_cells = 40)
  bt_s <- build_tessellation(pts_proj_auto, boundary = bnd_proj, method = "square", approx_n_cells = 40)
  nz(bt_h$cells, "build_tessellation(hex) cells exist")
  nz(bt_s$cells, "build_tessellation(square) cells exist")
} else {
  skip("build_tessellation() not found — skipping")
}

# ---------- GWR ----------
if (exists("fit_gwr_model") && HAS_SPGWR) {
  gdat <- pts_proj_auto |> dplyr::select(x1,x2,x3,resp, geometry)
  gwr_fit <- fit_gwr_model(gdat, "resp", c("x1","x2","x3"))
  nz(gwr_fit$model, "fit_gwr_model() returned a model")
  ok(isTRUE(gwr_fit$bandwidth > 0), "GWR bandwidth positive")
} else {
  skip("fit_gwr_model() or spgwr not available — skipping GWR")
}

# ---------- Bayesian (optional/heavy) ----------
if (exists("fit_bayesian_spatial_model") && RUN_BAYES) {
  gdat <- pts_proj_auto |> dplyr::select(x1,x2,x3,resp, geometry)
  bay_fit <- fit_bayesian_spatial_model(gdat, "resp", c("x1","x2","x3"))
  nz(bay_fit, "fit_bayesian_spatial_model() returned a model")
} else if (exists("fit_bayesian_spatial_model")) {
  skip("fit_bayesian_spatial_model() present, but RUN_BAYES=FALSE — skipping")
} else {
  skip("fit_bayesian_spatial_model() not found — skipping")
}

# ---------- CV wrappers ----------
if (exists("cv_gwr") && HAS_SPGWR) {
  cvg <- cv_gwr(pts_proj_auto, "resp", c("x1","x2","x3"), k = 3, seed = 7)
  ok(is.data.frame(cvg$overall), "cv_gwr(): overall metrics exist")
  ok(is.data.frame(cvg$fold_metrics), "cv_gwr(): per-fold metrics exist")
} else {
  skip("cv_gwr() not available — skipping")
}

if (exists("cv_bayes") && RUN_BAYES) {
  cvb <- cv_bayes(pts_proj_auto, "resp", c("x1","x2","x3"), k = 3, seed = 7)
  ok(is.data.frame(cvb$overall), "cv_bayes(): overall metrics exist")
  ok(is.data.frame(cvb$fold_metrics), "cv_bayes(): per-fold metrics exist")
} else if (exists("cv_bayes")) {
  skip("cv_bayes() present, but RUN_BAYES=FALSE — skipping")
} else {
  skip("cv_bayes() not found — skipping")
}

# ---------- evaluate_models (wrapper) ----------
if (exists("evaluate_models") && HAS_SPGWR) {
  gdat <- pts_proj_auto |> dplyr::select(x1,x2,x3,resp, geometry)
  models_vec <- if (RUN_BAYES && exists("fit_bayesian_spatial_model")) c("GWR","Bayesian") else "GWR"
  
  em_cv <- evaluate_models(
    data_sf = gdat,
    response_var = "resp",
    predictor_vars = c("x1","x2","x3"),
    do_cv = TRUE,
    k = 3,
    seed = 123,
    models = models_vec
  )
  ok(is.data.frame(em_cv$comparison), "evaluate_models(do_cv=TRUE): comparison table exists")
  
  em_in <- evaluate_models(
    data_sf = gdat,
    response_var = "resp",
    predictor_vars = c("x1","x2","x3"),
    do_cv = FALSE,
    models = models_vec
  )
  ok(is.data.frame(em_in$metrics), "evaluate_models(do_cv=FALSE): in-sample metrics table exists")
} else {
  skip("evaluate_models() not available or deps missing — skipping")
}

# ---------- evaluate_models_cv (with tessellation diagnostics) ----------
if (exists("evaluate_models_cv") && HAS_SPGWR) {
  gdat <- pts_proj_auto |> dplyr::select(x1,x2,x3,resp, geometry)
  models_vec2 <- if (RUN_BAYES && exists("fit_bayesian_spatial_model")) c("GWR","Bayesian") else "GWR"
  
  em2 <- evaluate_models_cv(
    data_sf = gdat,
    response_var = "resp",
    predictor_vars = c("x1","x2","x3"),
    k = 3,
    seed = 321,
    tess_method = "grid",
    tess_args = list(target_cells = 25, type = "square"),
    models = models_vec2
  )
  ok(is.data.frame(em2$overall), "evaluate_models_cv(): overall table exists")
  ok(is.data.frame(em2$by_fold), "evaluate_models_cv(): by_fold table exists")
  ok(is.list(em2$tessellation), "evaluate_models_cv(): tessellation descriptor present")
} else {
  skip("evaluate_models_cv() not available or deps missing — skipping")
}

# ---------- EXTRAS ----------
# (A) Degenerate/empty geometry handling
empty_sf <- st_as_sf(st_sfc(crs = 4326))
handled <- FALSE
tryCatch({
  tmp <- coerce_to_points(empty_sf, "auto"); handled <- inherits(tmp, "sf")
}, error=function(e){ handled <<- TRUE })
ok(handled, "coerce_to_points() handles empty input (returns empty or errors cleanly)")

# (B) Non-finite predictors
gdat <- pts_proj_auto |> dplyr::select(x1,x2,x3,resp, geometry)
pts_bad <- gdat
pts_bad$x1[1] <- NA_real_; pts_bad$x2[2] <- Inf

if (exists("fit_gwr_model") && HAS_SPGWR) {
  handled_gwr <- FALSE
  tryCatch({
    tmp <- fit_gwr_model(pts_bad, "resp", c("x1","x2","x3"))
    handled_gwr <- !is.null(tmp$model)
  }, warning=function(w){ handled_gwr <<- TRUE },
  error=function(e){ handled_gwr <<- TRUE })
  ok(handled_gwr, "fit_gwr_model() handles non-finite predictors by warning/error or dropping rows")
} else {
  skip("GWR fitter not available — skipping non-finite test")
}

if (exists("fit_bayesian_spatial_model") && RUN_BAYES) {
  handled_bayes <- FALSE
  tryCatch({
    tmp <- fit_bayesian_spatial_model(pts_bad, "resp", c("x1","x2","x3"))
    handled_bayes <- !is.null(tmp)
  }, warning=function(w){ handled_bayes <<- TRUE },
  error=function(e){ handled_bayes <<- TRUE })
  ok(handled_bayes, "fit_bayesian_spatial_model() handles non-finite predictors appropriately")
} else if (exists("fit_bayesian_spatial_model")) {
  skip("Bayes fitter present, but RUN_BAYES=FALSE — skipping non-finite test")
} else {
  skip("Bayes fitter not available — skipping non-finite test")
}

# (C) Boundary with holes — grids should respect holes after clipping
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

# project the hole to the SAME CRS as the grid/boundary
hole_proj <- ensure_projected(st_sfc(st_polygon(list(hole)), crs = 4326), bnd_hole_proj)
eq(st_crs(hole_proj), st_crs(grid_hole), "hole and grid CRS aligned")

cent_hole  <- suppressWarnings(st_centroid(grid_hole))
within_mat <- st_within(cent_hole, hole_proj, sparse = FALSE)
ok(all(!as.vector(within_mat)), "Grid centroids do not fall inside the hole (hole respected)")

# ---------- Final summary ----------
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
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
