# tests/testthat/helper-core-frozen.R
# ---------------------------------------------------------------------------
# The frozen baseline of the cell pipeline, for test-core-frozen.R and
# dev/baseline-core-frozen.R.
#
# core_frozen_calls() is a fixed list of calls to resolution_profile(),
# select_resolution(), summary() of a profile, determine_optimal_levels(),
# the seeding functions, build_tessellation(), create_voronoi_polygons(),
# assign_features_to_polygons() and summarize_by_cell().
# dev/baseline-core-frozen.R runs the list on the baseline code and saves,
# for every call, what it returns, prints, raises and logs, and where it
# leaves the random-number stream (tests/testthat/fixtures/core-frozen.rds).
# The test runs the list again on the code under test and compares.
#
# A call names the arguments it is about and leaves every other one to its
# default, so the list pins what the defaults return.  When a default
# changes, the calls written for the old value have to name it.  The fixture
# stores the text of every call and the test compares it with this list, so
# a call cannot be edited without the fixture saying so.
# ---------------------------------------------------------------------------

# The layers the calls read, from the data frames the fixture stores.  The
# data are stored and not drawn again, so they are the same numbers on every
# platform.
core_frozen_layers <- function(data) {
  as_points <- function(d, crs = 32632) sf::st_as_sf(d, coords = c("x", "y"), crs = crs)
  box <- function(layer) sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(layer)))
  pts <- as_points(data$pts)
  bnd <- box(pts)
  ll  <- as_points(data$ll, 4326)
  # Regions a caller brings: the boundary cut into 4 x 4 squares, with an ID
  # of each class a caller's polygons carry.
  bb <- sf::st_bbox(pts)
  xs <- bb[["xmin"]] + (0:4) * (bb[["xmax"]] - bb[["xmin"]]) / 4
  ys <- bb[["ymin"]] + (0:4) * (bb[["ymax"]] - bb[["ymin"]]) / 4
  at <- expand.grid(i = 1:4, j = 1:4)
  squares <- sf::st_sfc(lapply(seq_len(nrow(at)), function(r) {
    i <- at$i[r]; j <- at$j[r]
    sf::st_polygon(list(rbind(c(xs[i], ys[j]), c(xs[i + 1L], ys[j]), c(xs[i + 1L], ys[j + 1L]),
                              c(xs[i], ys[j + 1L]), c(xs[i], ys[j]))))
  }), crs = 32632)
  units <- sf::st_sf(GEOID = sprintf("U%02d", 1:16), id_int = 1:16, id_dbl = (1:16) * 1e5,
                     id_fac = factor(sprintf("f%02d", 1:16)), geometry = squares)
  list(pts  = pts,                       # 400 points, a response with spatial structure
       clus = as_points(data$clus),      # 500 points, most of them in four clusters
       big  = as_points(data$big),       # 1,800 points: more than `sample_n`
       bnd  = bnd,
       ll      = ll,                     # `pts` in longitude and latitude
       bnd_ll  = box(ll),
       twice   = as_points(data$twice),  # 150 locations, each held by two rows
       lattice = as_points(data$lattice),  # a 12 x 12 lattice; a constant and an integer column
       units   = units)
}

# The calls, in the order they run.  A call reads the layers and the results
# of the calls before it by name.  The first 154 are the list the baseline
# was first taken with; the rest were added after it and are appended, so
# that the names and the order of the first are what they were.
core_frozen_calls <- function() {
  calls <- list()
  add <- function(name, expr) {
    if (!is.null(calls[[name]])) stop("core_frozen_calls(): two calls are named ", name)
    calls[[name]] <<- expr
  }
  q <- function(name, expr) add(name, substitute(expr))
  # A profile read by every criterion, and summarised.
  readings <- function(p) {
    short <- sub("^profile_", "", p)
    prof  <- as.name(p)
    for (cr in c("cp", "reliability", "elbow", "moran_z"))
      add(sprintf("select_%s_%s", short, cr),
          bquote(select_resolution(.(prof), criterion = .(cr))))
    add(sprintf("select_%s_moran_z_diagnostic", short),
        bquote(select_resolution(.(prof), criterion = "moran_z", diagnostic = TRUE)))
    add(sprintf("summary_%s", short), bquote(summary(.(prof))))
  }

  # ---- resolution_profile() ------------------------------------------------
  q("profile_default",    resolution_profile(pts, "z"))
  q("profile_geometry",   resolution_profile(clus))
  q("profile_levels",     resolution_profile(pts, "z", levels = c(4, 8, 16, 32)))
  q("profile_split",      resolution_profile(pts, "z", select_on = "split"))
  q("profile_missing",    resolution_profile(pts, "z_na"))
  q("profile_min_cell_n", resolution_profile(pts, "z", min_cell_n = 5))
  q("profile_n_levels",   resolution_profile(pts, "z", n_levels = 16))
  q("profile_density",    resolution_profile(pts, "z", partition = "density"))
  q("profile_sampled",    resolution_profile(big, "z"))
  q("profile_logical",    resolution_profile(pts, "flag"))
  q("profile_no_value",   resolution_profile(pts, "none"))
  q("profile_predictors", resolution_profile(pts, "z", predictor_vars = c("x1", "x2")))
  q("profile_noise",      resolution_profile(pts, "noise"))
  q("profile_density_clustered",
    resolution_profile(clus, "w", partition = "density"))

  # ---- select_resolution() and summary(), on every profile -----------------
  for (p in names(calls)) readings(p)
  q("select_default_unnamed", select_resolution(profile_default))
  q("select_column_subset",
    select_resolution(profile_default[, c("levels", "wss", "cp")], criterion = "cp"))
  q("select_column_subset_without",
    select_resolution(profile_default[, c("levels", "wss")], criterion = "cp"))
  q("select_row_subset", select_resolution(profile_default[3:9, ], criterion = "cp"))

  # ---- get_voronoi_seeds() ---------------------------------------------------
  q("seeds_number",    get_voronoi_seeds(method = "kmeans", n = 30, sample_points = pts))
  q("seeds_number_seeded",
    get_voronoi_seeds(method = "kmeans", n = 30, sample_points = pts, set_seed = 11))
  q("seeds_selection",
    get_voronoi_seeds(method = "kmeans", n = select_default_cp, sample_points = pts))
  q("seeds_profile",
    get_voronoi_seeds(method = "kmeans", n = profile_default, sample_points = pts))
  q("seeds_density",
    get_voronoi_seeds(method = "kmeans", n = 30, sample_points = clus,
                      partition = "density"))
  q("seeds_density_selection",
    get_voronoi_seeds(method = "kmeans", n = select_density_cp, sample_points = pts))
  q("seeds_density_profile",
    get_voronoi_seeds(method = "kmeans", n = profile_density_clustered,
                      sample_points = clus))
  q("seeds_noise_selection",
    get_voronoi_seeds(method = "kmeans", n = select_noise_cp, sample_points = pts))
  q("seeds_boundary",  get_voronoi_seeds(boundary = bnd, method = "kmeans", n = 20))
  q("seeds_random",    get_voronoi_seeds(boundary = bnd, method = "random", n = 25))
  q("seeds_provided",  get_voronoi_seeds(method = "provided", seeds = pts[1:20, ]))

  # ---- voronoi_seeds_kmeans() ------------------------------------------------
  q("kmeans_number",    voronoi_seeds_kmeans(pts, 30))
  q("kmeans_selection", voronoi_seeds_kmeans(pts, select_default_cp))
  q("kmeans_profile",   voronoi_seeds_kmeans(pts, profile_default))
  q("kmeans_density",   voronoi_seeds_kmeans(clus, 30, partition = "density"))
  q("kmeans_density_selection", voronoi_seeds_kmeans(pts, select_density_cp))
  q("kmeans_density_profile",
    voronoi_seeds_kmeans(clus, profile_density_clustered))

  # ---- build_tessellation() --------------------------------------------------
  q("cells_number",    build_tessellation(seeds_number, method = "voronoi"))
  q("cells_selection", build_tessellation(seeds_selection, method = "voronoi"))
  q("cells_profile",   build_tessellation(seeds_profile, method = "voronoi"))
  q("cells_density",   build_tessellation(seeds_density, method = "voronoi"))
  q("cells_random",    build_tessellation(seeds_random, boundary = bnd, method = "voronoi"))
  q("cells_provided",  build_tessellation(seeds_provided, method = "voronoi"))
  q("cells_kmeans_number",    build_tessellation(kmeans_number, method = "voronoi"))
  q("cells_kmeans_selection", build_tessellation(kmeans_selection, method = "voronoi"))
  q("cells_kmeans_density",   build_tessellation(kmeans_density, method = "voronoi"))
  q("cells_hex_number",
    build_tessellation(pts, boundary = bnd, method = "hex", approx_n_cells = 12))
  q("cells_hex_profile",
    build_tessellation(pts, boundary = bnd, method = "hex",
                       approx_n_cells = profile_default))

  # ---- assign_features_to_polygons() -----------------------------------------
  q("assigned_number",
    assign_features_to_polygons(pts, cells_number$cells, polygon_id_col = "cell_id"))
  q("assigned_selection",
    assign_features_to_polygons(pts, cells_selection$cells, polygon_id_col = "cell_id"))
  q("assigned_density",
    assign_features_to_polygons(clus, cells_density$cells, polygon_id_col = "cell_id"))
  q("assigned_dropped",
    assign_features_to_polygons(pts, cells_number$cells[1:20, ],
                                polygon_id_col = "cell_id"))
  q("assigned_kept",
    assign_features_to_polygons(pts, cells_number$cells[1:20, ],
                                polygon_id_col = "cell_id", keep_unassigned = TRUE))
  # One cell left with a single row and one with none.
  q("assigned_sparse", {
    id   <- assigned_number$cell_id
    ids  <- sort(unique(id))
    drop <- id == ids[2L] | (id == ids[1L] & duplicated(id))
    assigned_number[!drop, ]
  })

  # ---- summarize_by_cell() ---------------------------------------------------
  q("values_default",   summarize_by_cell(assigned_number, "z", id_col = "cell_id"))
  q("values_no_response", summarize_by_cell(assigned_number, id_col = "cell_id"))
  q("values_kish",
    summarize_by_cell(assigned_number, "z", id_col = "cell_id", deff = "kish"))
  q("sac_z", estimate_sac_range(pts, "z"))
  q("values_sac",
    summarize_by_cell(assigned_number, "z", id_col = "cell_id", deff = "variogram",
                      sac = sac_z))
  q("values_variogram",
    summarize_by_cell(assigned_number, "z", id_col = "cell_id", deff = "variogram"))
  q("values_deff_2",
    summarize_by_cell(assigned_number, "z", id_col = "cell_id", deff = 2))
  q("values_conf_level",
    summarize_by_cell(assigned_number, "z", id_col = "cell_id", conf_level = 0.90))
  q("values_cells_sf",
    summarize_by_cell(assigned_number, "z", id_col = "cell_id",
                      cells_sf = cells_number$cells))
  q("values_area",
    summarize_by_cell(assigned_number, "z", id_col = "cell_id",
                      cells_sf = cells_number$cells, area = TRUE))
  q("values_median",
    summarize_by_cell(assigned_number, "z", id_col = "cell_id", agg_funs = median))
  q("values_predictors",
    summarize_by_cell(assigned_number, "z", predictor_vars = c("x1", "x2"),
                      id_col = "cell_id", conf_level = 0.95))
  q("values_nonfinite",
    summarize_by_cell(assigned_number, "z_bad", id_col = "cell_id", conf_level = 0.95))
  q("values_sparse",
    summarize_by_cell(assigned_sparse, "z", id_col = "cell_id",
                      cells_sf = cells_number$cells, conf_level = 0.95))
  q("values_unassigned", summarize_by_cell(assigned_kept, "z", id_col = "cell_id"))
  q("values_data_frame",
    summarize_by_cell(sf::st_drop_geometry(assigned_number), "z", id_col = "cell_id"))
  q("values_selection",
    summarize_by_cell(assigned_selection, "z", id_col = "cell_id",
                      cells_sf = cells_selection$cells, deff = "kish",
                      conf_level = 0.95))
  q("values_density",
    summarize_by_cell(assigned_density, "w", id_col = "cell_id",
                      cells_sf = cells_density$cells))

  # ===========================================================================
  # Added after the first list.  Present calls on the paths the later work
  # rewrites and on inputs the first list does not read.
  # ===========================================================================

  # ---- more profiles: other present arguments, other inputs ------------------
  first <- names(calls)
  q("profile_geometry_sac",   resolution_profile(clus, sac = sac_z))
  q("profile_sac",            resolution_profile(pts, "z", sac = sac_z))
  q("profile_range_floor",    resolution_profile(pts, "z", range_floor = TRUE))
  q("profile_arguments",
    resolution_profile(pts, "z", sample_n = 200, nstart = 5, seed = 9, init = "random",
                       iter_max = 30))
  q("profile_density_arguments",
    resolution_profile(pts, "z", partition = "density", density_power = 0.5,
                       density_k = 4, nstart = 3))
  q("profile_lonlat",          resolution_profile(ll, "z"))
  q("profile_lonlat_density",  resolution_profile(ll, "z", partition = "density"))
  q("profile_twice",           resolution_profile(twice, "z"))
  q("profile_twice_density",   resolution_profile(twice, "z", partition = "density"))
  q("profile_lattice",         resolution_profile(lattice, "z"))
  q("profile_lattice_density", resolution_profile(lattice, "z", partition = "density"))
  q("profile_constant",        resolution_profile(lattice, "const"))
  q("profile_integer",         resolution_profile(lattice, "count"))
  for (p in setdiff(names(calls), first)) readings(p)
  # Layers too small for a ladder, and what is not a layer.
  q("small_five_rows",  resolution_profile(pts[1:5, ], "z"))
  q("small_five_rows_selected", select_resolution(small_five_rows))
  q("small_five_rows_count",
    get_voronoi_seeds(method = "kmeans", n = small_five_rows, sample_points = pts[1:5, ]))
  q("small_one_row",    resolution_profile(pts[1L, ], "z"))
  q("small_no_rows",    resolution_profile(pts[0L, ], "z"))
  q("small_data_frame", resolution_profile(sf::st_drop_geometry(pts), "z"))
  # Other present arguments of the readers.
  q("select_default_tol", select_resolution(profile_default, tol = 0.1, range_floor = FALSE))
  q("select_column_subset_cp",
    select_resolution(profile_default[, c("levels", "cp")], criterion = "cp"))
  q("summary_default_criteria",
    summary(profile_default, criteria = c("cp", "elbow"), tol = 0.1))
  q("print_default_digits", utils::capture.output(print(profile_default, digits = 5)))
  q("levels_default", determine_optimal_levels(pts, max_levels = 8L))
  q("levels_clustered", determine_optimal_levels(clus, max_levels = 8L))
  q("levels_morans",
    determine_optimal_levels(pts, max_levels = 8L, response_var = "z",
                             criterion = "morans_i"))
  q("levels_combined",
    determine_optimal_levels(pts, max_levels = 8L, response_var = "z",
                             criterion = "combined"))

  # ---- a count read off a profile, a selection or a vector --------------------
  q("count_seeds_no_value",
    get_voronoi_seeds(method = "kmeans", n = profile_no_value, sample_points = pts))
  q("count_kmeans_no_value", voronoi_seeds_kmeans(pts, profile_no_value))
  q("count_hex_no_value",
    build_tessellation(pts, boundary = bnd, method = "hex",
                       approx_n_cells = profile_no_value))
  q("count_seeds_geometry",
    get_voronoi_seeds(method = "kmeans", n = profile_geometry, sample_points = clus))
  q("count_kmeans_geometry", voronoi_seeds_kmeans(clus, profile_geometry))
  q("count_hex_geometry",
    build_tessellation(clus, boundary = bnd, method = "hex",
                       approx_n_cells = profile_geometry))
  q("count_seeds_no_cp",
    get_voronoi_seeds(method = "kmeans", n = profile_geometry_sac, sample_points = clus))
  q("count_kmeans_no_cp", voronoi_seeds_kmeans(clus, profile_geometry_sac))
  q("count_hex_no_cp",
    build_tessellation(clus, boundary = bnd, method = "hex",
                       approx_n_cells = profile_geometry_sac))
  q("count_seeds_noise",
    get_voronoi_seeds(method = "kmeans", n = profile_noise, sample_points = pts))
  q("count_kmeans_noise", voronoi_seeds_kmeans(pts, profile_noise))
  q("count_hex_noise",
    build_tessellation(pts, boundary = bnd, method = "hex", approx_n_cells = profile_noise))
  q("count_seeds_levels",
    get_voronoi_seeds(method = "kmeans", n = profile_levels, sample_points = pts))
  q("count_kmeans_levels_selection", voronoi_seeds_kmeans(pts, select_levels_cp))
  q("count_hex_selection",
    build_tessellation(pts, boundary = bnd, method = "hex",
                       approx_n_cells = select_default_cp))
  q("count_seeds_vector",
    get_voronoi_seeds(method = "kmeans", n = c(12, 20, 8), sample_points = pts,
                      set_seed = 3))
  q("count_kmeans_vector", voronoi_seeds_kmeans(pts, c(12, 20)))
  q("count_seeds_candidates",
    get_voronoi_seeds(method = "kmeans", n = levels_default, sample_points = pts,
                      set_seed = 3))
  q("count_seeds_other_layer",
    get_voronoi_seeds(method = "kmeans", n = select_default_cp, sample_points = clus))
  q("count_seeds_selection_density_named",
    get_voronoi_seeds(method = "kmeans", n = select_default_cp, sample_points = pts,
                      partition = "density"))
  q("count_kmeans_selection_kmeans_named",
    voronoi_seeds_kmeans(pts, select_density_cp, partition = "kmeans"))

  # ---- seeds: other present arguments, other inputs ---------------------------
  q("seeds_arguments",
    get_voronoi_seeds(method = "kmeans", n = 25, sample_points = pts, set_seed = 3,
                      kmeans_nstart = 2, kmeans_iter = 5, init = "kmeans++"))
  q("seeds_one", get_voronoi_seeds(method = "kmeans", n = 1, sample_points = pts, set_seed = 3))
  q("seeds_lonlat",
    get_voronoi_seeds(method = "kmeans", n = 25, sample_points = ll, set_seed = 3))
  q("seeds_lonlat_density",
    get_voronoi_seeds(method = "kmeans", n = 25, sample_points = ll, set_seed = 3,
                      partition = "density"))
  q("seeds_lonlat_selection",
    get_voronoi_seeds(method = "kmeans", n = select_lonlat_cp, sample_points = ll))
  q("seeds_boundary_lonlat",
    get_voronoi_seeds(boundary = bnd_ll, method = "kmeans", n = 15, set_seed = 3))
  q("seeds_random_lonlat",
    get_voronoi_seeds(boundary = bnd_ll, method = "random", n = 15, set_seed = 3))
  q("seeds_twice",
    get_voronoi_seeds(method = "kmeans", n = 20, sample_points = twice, set_seed = 3))
  q("seeds_twice_density",
    get_voronoi_seeds(method = "kmeans", n = 20, sample_points = twice, set_seed = 3,
                      partition = "density"))
  q("seeds_more_than_locations",
    get_voronoi_seeds(method = "kmeans", n = 400, sample_points = twice, set_seed = 3))
  q("seeds_lattice",
    get_voronoi_seeds(method = "kmeans", n = 16, sample_points = lattice, set_seed = 3))
  q("seeds_lattice_density",
    get_voronoi_seeds(method = "kmeans", n = 16, sample_points = lattice, set_seed = 3,
                      partition = "density"))
  q("seeds_random_density_named",
    get_voronoi_seeds(boundary = bnd, method = "random", n = 5, set_seed = 3,
                      partition = "density"))
  q("seeds_random_init_named",
    get_voronoi_seeds(boundary = bnd, method = "random", n = 5, set_seed = 3,
                      init = "kmeans++"))
  q("seeds_random_power_named",
    get_voronoi_seeds(boundary = bnd, method = "random", n = 5, set_seed = 3,
                      density_power = 0.5))
  q("seeds_provided_density_named",
    get_voronoi_seeds(method = "provided", seeds = pts[1:20, ], partition = "density"))
  q("seeds_kmeans_power_named",
    get_voronoi_seeds(method = "kmeans", n = 30, sample_points = pts, set_seed = 3,
                      density_power = 0.5))
  q("kmeans_arguments",
    voronoi_seeds_kmeans(pts, 25, set_seed = 1, nstart = 2, init = "kmeans++", iter_max = 8))
  q("kmeans_unseeded",       voronoi_seeds_kmeans(pts, 30, set_seed = NULL))
  q("kmeans_lonlat",         voronoi_seeds_kmeans(ll, 25))
  q("kmeans_lonlat_density", voronoi_seeds_kmeans(ll, 25, partition = "density"))
  q("kmeans_more_than_locations", voronoi_seeds_kmeans(twice, 200))
  q("kmeans_lattice_density", voronoi_seeds_kmeans(lattice, 16, partition = "density"))
  q("random_number", voronoi_seeds_random(bnd, 12))
  q("random_seeded", voronoi_seeds_random(bnd, 12, set_seed = 5))

  # ---- cells on those seeds ---------------------------------------------------
  q("cells_lonlat",   build_tessellation(seeds_lonlat, method = "voronoi"))
  q("cells_lonlat_selection", build_tessellation(seeds_lonlat_selection, method = "voronoi"))
  q("cells_twice",    build_tessellation(seeds_twice, method = "voronoi"))
  q("cells_lattice",  build_tessellation(seeds_lattice, method = "voronoi"))
  q("cells_one",      build_tessellation(seeds_one, method = "voronoi"))
  q("cells_boundary", build_tessellation(seeds_boundary, boundary = bnd, method = "voronoi"))
  q("cells_kmeans_profile", build_tessellation(kmeans_profile, method = "voronoi"))
  q("cells_density_selection",
    build_tessellation(seeds_density_selection, method = "voronoi"))
  q("cells_points",   build_tessellation(pts[1:40, ], method = "voronoi"))
  q("cells_quiet",    build_tessellation(seeds_number, method = "voronoi", quiet = TRUE))
  q("cells_square",
    build_tessellation(pts, boundary = bnd, method = "square", cellsize = 250))
  q("polygons_number", create_voronoi_polygons(seeds_number))

  # ---- assignments ------------------------------------------------------------
  q("assigned_default_id", assign_features_to_polygons(pts, cells_number$cells))
  q("assigned_lonlat",
    assign_features_to_polygons(ll, cells_lonlat$cells, polygon_id_col = "cell_id"))
  q("assigned_twice",
    assign_features_to_polygons(twice, cells_twice$cells, polygon_id_col = "cell_id"))
  q("assigned_lattice",
    assign_features_to_polygons(lattice, cells_lattice$cells, polygon_id_col = "cell_id"))
  q("assigned_units_character",
    assign_features_to_polygons(pts, units, polygon_id_col = "GEOID"))
  q("assigned_units_integer",
    assign_features_to_polygons(pts, units, polygon_id_col = "id_int"))
  q("assigned_units_double",
    assign_features_to_polygons(pts, units, polygon_id_col = "id_dbl"))
  q("assigned_units_factor",
    assign_features_to_polygons(pts, units, polygon_id_col = "id_fac"))
  q("assigned_units_lonlat",
    assign_features_to_polygons(ll, units, polygon_id_col = "GEOID"))
  q("assigned_lattice_units",
    assign_features_to_polygons(lattice, units, polygon_id_col = "GEOID",
                                tie_break = "first"))
  q("assigned_one_region",       { a <- assigned_number; a$cell_id <- 1L; a })
  q("assigned_one_row_regions",  assigned_number[!duplicated(assigned_number$cell_id), ])
  q("assigned_one_row",          assigned_number[1L, ])
  q("assigned_no_rows",          assigned_number[0L, ])

  # ---- summaries: other present arguments --------------------------------------
  q("values_id_found",     summarize_by_cell(assigned_number, "z"))
  q("values_not_quiet",
    summarize_by_cell(assigned_number, "z", id_col = "cell_id", quiet = FALSE))
  q("values_not_quiet_variogram",
    summarize_by_cell(assigned_number, "z", id_col = "cell_id", deff = "variogram",
                      quiet = FALSE))
  q("values_other_cells",
    summarize_by_cell(assigned_number, "z", id_col = "cell_id",
                      cells_sf = cells_number$cells[1:20, ]))
  q("values_area_missing",
    summarize_by_cell(assigned_number, "z_na", id_col = "cell_id",
                      cells_sf = cells_number$cells, area = TRUE))
  q("values_variogram_cells_reversed",
    summarize_by_cell(assigned_number, "z", id_col = "cell_id", deff = "variogram",
                      sac = sac_z, cells_sf = cells_number$cells[30:1, ]))
  q("values_variogram_max_n",
    summarize_by_cell(assigned_number, "z", id_col = "cell_id", deff = "variogram",
                      sac = sac_z, cells_sf = cells_number$cells, deff_max_n = 5))
  q("values_kish_predictors_area",
    summarize_by_cell(assigned_number, "z_na", c("x1", "x2"), id_col = "cell_id",
                      deff = "kish", cells_sf = cells_number$cells, conf_level = 0.95,
                      area = TRUE))
  q("values_predictors_only",
    summarize_by_cell(assigned_number, predictor_vars = c("x1", "x2"), id_col = "cell_id"))
  q("values_logical", summarize_by_cell(assigned_number, "flag", id_col = "cell_id"))
  q("values_functions",
    summarize_by_cell(assigned_number, "z_na", id_col = "cell_id",
                      agg_funs = list(mean = function(x) mean(x, na.rm = TRUE),
                                      top = function(x) max(x, na.rm = TRUE))))
  q("values_function_names",
    summarize_by_cell(assigned_number, "z_na", id_col = "cell_id",
                      agg_funs = c("median", "max")))
  q("values_data_frame_variogram",
    summarize_by_cell(sf::st_drop_geometry(assigned_number), "z", id_col = "cell_id",
                      deff = "variogram"))
  # The cells carry the ID as text and the rows as a number.
  q("values_id_number_text", {
    rows  <- assigned_number
    rows$cell_dbl <- rows$cell_id * 1e5
    cells <- cells_number$cells
    cells$cell_dbl <- sprintf("%.0f", cells$cell_id * 1e5)
    summarize_by_cell(rows, "z", id_col = "cell_dbl", cells_sf = cells)
  })

  # ---- summaries: other inputs --------------------------------------------------
  q("values_lonlat",
    summarize_by_cell(assigned_lonlat, "z", id_col = "cell_id",
                      cells_sf = cells_lonlat$cells))
  q("values_lonlat_variogram",
    summarize_by_cell(assigned_lonlat, "z", id_col = "cell_id",
                      cells_sf = cells_lonlat$cells, deff = "variogram"))
  q("values_lonlat_area",
    summarize_by_cell(assigned_lonlat, "z", id_col = "cell_id",
                      cells_sf = cells_lonlat$cells, area = TRUE))
  q("values_twice",
    summarize_by_cell(assigned_twice, "z", id_col = "cell_id", deff = "variogram",
                      conf_level = 0.90))
  q("values_lattice",
    summarize_by_cell(assigned_lattice, "z", id_col = "cell_id",
                      cells_sf = cells_lattice$cells, conf_level = 0.90))
  q("values_constant",
    summarize_by_cell(assigned_lattice, "const", id_col = "cell_id", deff = "kish",
                      conf_level = 0.90))
  q("values_integer", summarize_by_cell(assigned_lattice, "count", id_col = "cell_id"))
  q("values_units_character",
    summarize_by_cell(assigned_units_character, "z", id_col = "GEOID", cells_sf = units))
  q("values_units_integer",
    summarize_by_cell(assigned_units_integer, "z", id_col = "id_int", cells_sf = units,
                      area = TRUE))
  q("values_units_double",
    summarize_by_cell(assigned_units_double, "z", id_col = "id_dbl", cells_sf = units))
  q("values_units_factor",
    summarize_by_cell(assigned_units_factor, "z", id_col = "id_fac", cells_sf = units,
                      deff = "kish"))
  q("values_units_lonlat",
    summarize_by_cell(assigned_units_lonlat, "z", id_col = "GEOID", cells_sf = units,
                      conf_level = 0.90))
  q("values_one_region",
    summarize_by_cell(assigned_one_region, "z", id_col = "cell_id", deff = "kish",
                      conf_level = 0.90))
  q("values_one_row_regions",
    summarize_by_cell(assigned_one_row_regions, "z", id_col = "cell_id",
                      cells_sf = cells_number$cells, deff = "kish", conf_level = 0.90))
  q("values_one_row",
    summarize_by_cell(assigned_one_row, "z", id_col = "cell_id",
                      cells_sf = cells_number$cells))
  q("values_no_rows",
    summarize_by_cell(assigned_no_rows, "z", id_col = "cell_id",
                      cells_sf = cells_number$cells))
  calls
}

# The text of the calls, as the fixture stores it.
core_frozen_call_text <- function(calls = core_frozen_calls())
  lapply(calls, function(e) paste(deparse(e, width.cutoff = 500L), collapse = " "))

# What decides the last bits of a result: the platform, R, sf and the
# libraries behind it, the optional packages a call reads, and the linear
# algebra library.
core_frozen_platform <- function() {
  ext <- sf::sf_extSoftVersion()
  ver <- function(p) if (requireNamespace(p, quietly = TRUE))
    as.character(utils::packageVersion(p)) else NA_character_
  list(platform = R.version$platform,
       r        = R.version.string,
       sf       = ver("sf"),
       geos     = ext[["GEOS"]],
       proj     = ext[["PROJ"]],
       gdal     = ext[["GDAL"]],
       packages = vapply(c("gstat", "FNN", "Matrix", "dplyr", "s2"), ver, character(1)),
       blas     = unname(extSoftVersion()[["BLAS"]]),
       lapack   = La_library())
}

# The entries of a platform record that two runs have to share for their
# results to be held to identical().
core_frozen_platform_keys <- function()
  c("platform", "r", "sf", "geos", "proj", "packages", "blas", "lapack")

# TRUE when two platform records agree on everything that can move a last
# bit, so that results can be held to identical().
core_frozen_same_platform <- function(a, b) {
  keys <- core_frozen_platform_keys()
  identical(a[keys], b[keys])
}

# The package's functions as they run where the packages in `pkgs` are not
# installed: a copy of every function that asks requireNamespace() for one
# of them, with the question answered FALSE, as a named list for
# testthat::local_mocked_bindings().  .has_backend(), once the package has
# it, answers FALSE for them too.  The copies keep their environment, so
# every call from one package function to another meets them.
core_frozen_without <- function(pkgs = c("FNN", "Matrix")) {
  ns <- asNamespace("spatialkit")
  asks <- function(e) is.call(e) && identical(e[[1L]], quote(requireNamespace)) &&
    length(e) >= 2L && is.character(e[[2L]]) && length(e[[2L]]) == 1L && e[[2L]] %in% pkgs
  found <- FALSE
  swap <- function(e) {
    if (asks(e)) { found <<- TRUE; return(FALSE) }
    if (is.call(e) || is.pairlist(e)) {
      for (i in seq_along(e)) {
        # An argument left empty, as in x[, 1], or one that is NULL.
        if (is.null(e[[i]]) || (is.symbol(e[[i]]) && !nzchar(as.character(e[[i]])))) next
        e[[i]] <- swap(e[[i]])
      }
    }
    e
  }
  out <- list()
  for (nm in ls(ns, all.names = TRUE)) {
    f <- get(nm, envir = ns)
    if (!is.function(f) || is.primitive(f)) next
    if (!("requireNamespace" %in% all.names(body(f))) &&
        !any(grepl("requireNamespace", deparse(args(f)), fixed = TRUE))) next
    found <- FALSE
    args <- swap(formals(f))
    code <- swap(body(f))
    if (!found) next
    g <- f
    if (!is.null(args)) formals(g) <- args
    body(g) <- code
    out[[nm]] <- g
  }
  if (exists(".has_backend", envir = ns, inherits = FALSE)) {
    has <- get(".has_backend", envir = ns)
    out[[".has_backend"]] <- function(pkg, ...) if (pkg %in% pkgs) FALSE else has(pkg, ...)
  }
  out
}

# Run the calls on `data` and return one record per call:
#   value       what the call returned (NULL after an error)
#   conditions  every warning, message and error, in order: when it was
#               raised ("call" or "print"), its kind, its class, its text,
#               the function its call names (NA when it carries none) and
#               the fields it carries besides
#   output      what the call wrote to standard output
#   printed     the print of a profile, a selection or a profile summary
#   log         the lines the package logged at WARN or above, the level its
#               console echo shows, without their time stamps
#   rng         where the call left the random-number stream: `untouched`
#               when it is the state the call started from, the state itself
#               otherwise
# Everything a result can depend on besides the code is fixed for the run and
# put back afterwards: the random-number generator (seeded anew before every
# call), the options and the locale that printing and sorting read, sf's use
# of s2, the package's record of what it has logged once in this session, and
# its console log.
#
# `without` names optional packages the run is to do without, although they
# are installed (core_frozen_without()).
core_frozen_run <- function(data, calls = core_frozen_calls(), without = NULL) {
  ns <- asNamespace("spatialkit")
  # Optional packages are loaded before anything is recorded, so that what
  # loading one says is no part of a call's record.
  for (p in names(core_frozen_platform()$packages))
    suppressWarnings(suppressMessages(requireNamespace(p, quietly = TRUE)))
  if (length(without))
    testthat::local_mocked_bindings(!!!core_frozen_without(without), .package = "spatialkit")

  old_opt <- options(width = 80, digits = 7, scipen = 0, OutDec = ".",
                     useFancyQuotes = FALSE, warn = 0)
  on.exit(options(old_opt), add = TRUE)
  old_collate <- Sys.getlocale("LC_COLLATE")
  Sys.setlocale("LC_COLLATE", "C")
  on.exit(Sys.setlocale("LC_COLLATE", old_collate), add = TRUE)
  old_lang <- Sys.getenv("LANGUAGE", unset = NA)
  Sys.setenv(LANGUAGE = "en")
  on.exit(if (is.na(old_lang)) Sys.unsetenv("LANGUAGE") else Sys.setenv(LANGUAGE = old_lang),
          add = TRUE)

  old_s2 <- suppressMessages(sf::sf_use_s2(TRUE))
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)

  stream <- function() if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv())
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- stream()
  old_kind <- RNGkind()
  on.exit({
    suppressWarnings(do.call(RNGkind, as.list(old_kind)))
    if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  }, add = TRUE)

  once <- get(".warned_once", envir = ns)
  old_once <- as.list(once, all.names = TRUE)
  clear_once <- function() rm(list = ls(once, all.names = TRUE), envir = once)
  on.exit({ clear_once(); list2env(old_once, envir = once) }, add = TRUE)

  logged <- character(0)
  old_threshold <- tryCatch(logger::log_threshold(namespace = "spatialkit", index = 2),
                            error = function(e) logger::WARN)
  on.exit({
    logger::log_appender(get(".sk_console_appender", envir = ns),
                         namespace = "spatialkit", index = 2)
    logger::log_threshold(old_threshold, namespace = "spatialkit", index = 2)
  }, add = TRUE)
  logger::log_appender(function(lines) logged <<- c(logged, lines),
                       namespace = "spatialkit", index = 2)
  logger::log_threshold(logger::WARN, namespace = "spatialkit", index = 2)

  env <- list2env(core_frozen_layers(data), envir = new.env(parent = ns))
  out <- stats::setNames(vector("list", length(calls)), names(calls))
  for (nm in names(calls)) {
    clear_once()
    logged <- character(0)
    set.seed(20261012, kind = "Mersenne-Twister", normal.kind = "Inversion",
             sample.kind = "Rejection")
    start <- stream()
    rec <- core_frozen_record(calls[[nm]], env)
    end <- stream()
    rec$log <- sub("^([A-Z]+) \\[[^]]*\\] ", "\\1 ", logged)
    rec$rng <- list(untouched = identical(end, start),
                    state = if (!identical(end, start)) end)
    assign(nm, rec$value, envir = env)
    out[[nm]] <- rec
  }
  out
}

# One call: its value, conditions, output and print.
core_frozen_record <- function(expr, env) {
  conditions <- list()
  # The fields a condition carries besides its text and its call: what a
  # handler reads.  A field that is not plain data (the backtrace or the
  # parent a condition of another package may hold) is kept by its class.
  plain <- function(v) is.null(v) || is.atomic(v) ||
    (is.list(v) && !is.object(v) && all(vapply(v, plain, logical(1))))
  fields <- function(cond) {
    f <- unclass(cond)
    f <- f[setdiff(names(f), c("message", "call"))]
    lapply(f, function(v) if (plain(v)) v else list(not_recorded = class(v)))
  }
  from <- function(cond) {
    cl <- conditionCall(cond)
    if (is.null(cl)) NA_character_
    else paste(deparse(if (is.call(cl)) cl[[1L]] else cl), collapse = " ")
  }
  note <- function(when, kind, cond)
    conditions[[length(conditions) + 1L]] <<-
      list(when = when, kind = kind, class = class(cond), message = conditionMessage(cond),
           call = from(cond), fields = fields(cond))
  watch <- function(code, when) withCallingHandlers(
    tryCatch(code, error = function(e) { note(when, "error", e); NULL }),
    warning = function(w) { note(when, "warning", w); invokeRestart("muffleWarning") },
    message = function(m) { note(when, "message", m); invokeRestart("muffleMessage") })
  value <- NULL
  output <- utils::capture.output(value <- watch(eval(expr, env), "call"), type = "output")
  printed <- NULL
  if (inherits(value, c("resolution_profile", "resolution_selection", "resolution_summary")))
    printed <- utils::capture.output(invisible(watch(print(value), "print")), type = "output")
  list(value = value, conditions = conditions, output = output, printed = printed)
}

# ---- Comparing two records ---------------------------------------------------

# The coordinates of a geometry set: for each geometry, its pieces (a point,
# a line, or one ring of a polygon) as two-column matrices.
core_frozen_coords <- function(g) {
  pieces <- function(x) {
    if (is.matrix(x)) list(unname(x[, 1:2, drop = FALSE]))
    else if (is.list(x)) unlist(lapply(x, pieces), recursive = FALSE)
    else if (is.numeric(x) && length(x) >= 2L) list(matrix(as.numeric(x)[1:2], 1L, 2L))
    else list()
  }
  lapply(g, function(x) { p <- pieces(x); if (is.null(p)) list() else p })
}

# The largest gap between two pieces with the same number of vertices.  A
# ring is the same polygon wherever it starts and whichever way it runs, and
# GEOS versions differ in both, so two rings are compared at the start and
# the direction that bring them closest.  The start is not taken from one
# ring alone (its smallest vertex, say): two vertices can tie for it, and a
# last bit then decides which one starts.
core_frozen_gap <- function(a, b) {
  n <- nrow(b)
  far <- function(u, v) {
    d <- abs(u - v)
    if (all(is.na(d))) 0 else max(d, na.rm = TRUE)
  }
  ring <- function(m) nrow(m) >= 4L && !anyNA(m) && all(m[1L, ] == m[nrow(m), ])
  if (!ring(a) || !ring(b)) return(far(a, b))
  a <- a[-n, , drop = FALSE]
  b <- b[-n, , drop = FALSE]
  n <- n - 1L
  best <- Inf
  for (turned in c(FALSE, TRUE)) {
    u <- if (turned) a[n:1, , drop = FALSE] else a
    for (s in seq_len(n)) {
      best <- min(best, far(u[c(s:n, seq_len(s - 1L)), , drop = FALSE], b))
      if (best == 0) return(0)
    }
  }
  best
}

# An object in three parts.  `shape` is the object with every number and
# every geometry taken out: names, classes, attribute names, types, lengths,
# and the values that are not numbers.  `numbers` are the numeric vectors and
# `coords` the geometries' coordinates, both named by where they sit.  A
# coordinate reference system is its name alone: its text changes with PROJ.
core_frozen_parts <- function(x, root = "value") {
  numbers <- list()
  coords  <- list()
  structural <- c("names", "class", "row.names", "dim", "dimnames", "levels",
                  "sf_column", "agr")
  walk <- function(x, path) {
    if (inherits(x, "crs")) return(list(crs = x$input))
    if (inherits(x, "sfg")) x <- sf::st_sfc(x)
    if (inherits(x, "sfc")) {
      coords[[path]] <<- core_frozen_coords(x)
      return(list(class = class(x), length = length(x),
                  crs = sf::st_crs(x)$input,
                  empty = which(sf::st_is_empty(x))))
    }
    a <- attributes(x)
    # In name order: assigning to a column moves `class` to the end.
    shape <- list(type = typeof(x), length = length(x), attributes = sort(names(a)),
                  structural = a[sort(intersect(names(a), structural))])
    for (nm in sort(setdiff(names(a), structural)))
      shape[[paste0("@", nm)]] <- walk(a[[nm]], paste0(path, "@", nm))
    if (is.function(x) || is.language(x)) {
      shape$value <- deparse(x)
    } else if (is.list(x)) {
      labels <- make.unique(if (is.null(names(x))) as.character(seq_along(x)) else names(x))
      shape$elements <- stats::setNames(lapply(seq_along(x), function(i)
        walk(x[[i]], paste0(path, "$", labels[i]))), labels)
    } else if (is.numeric(x) && !is.factor(x)) {
      numbers[[path]] <<- as.vector(unclass(x))
    } else {
      shape$value <- tryCatch(as.vector(unclass(x)),
                              error = function(e) paste(class(x), collapse = "/"))
    }
    shape
  }
  shape <- walk(x, root)
  list(shape = shape, numbers = numbers, coords = coords)
}

# The differences between two objects under the rule for another platform,
# as lines of text.  `root` names the object in them.
core_frozen_compare <- function(got, want, root = "value") {
  out <- character(0)
  g <- core_frozen_parts(got, root)
  w <- core_frozen_parts(want, root)
  # Where two shapes that are not identical differ: the parts only one of
  # them has, and the parts that hold other values.
  where <- function(a, b) {
    fa <- unlist(a); fb <- unlist(b)
    both <- intersect(names(fa), names(fb))
    moved <- both[!mapply(identical, as.character(fa[both]), as.character(fb[both]))]
    few <- function(x) paste(c(utils::head(x, 6L), if (length(x) > 6L) "..."), collapse = ", ")
    d <- c(if (length(setdiff(names(fa), names(fb))))
             paste("  only now:", few(setdiff(names(fa), names(fb)))),
           if (length(setdiff(names(fb), names(fa))))
             paste("  only on the baseline:", few(setdiff(names(fb), names(fa)))),
           if (length(moved))
             paste("  other values:", few(sprintf("%s (now %s, baseline %s)", moved,
                                                  fa[moved], fb[moved]))))
    if (length(d)) d else "  identical() is FALSE, and no part by name differs"
  }
  if (!identical(g$shape, w$shape))
    return(c(sprintf("%s: names, classes, attributes or values that are not numbers differ:",
                     root),
             where(g$shape, w$shape)))
  for (nm in names(w$numbers)) {
    d <- all.equal(g$numbers[[nm]], w$numbers[[nm]], tolerance = 1e-12)
    if (!isTRUE(d)) out <- c(out, sprintf("%s: %s", nm, paste(d, collapse = "; ")))
  }
  for (nm in names(w$coords)) {
    a <- g$coords[[nm]]; b <- w$coords[[nm]]
    rows <- function(p) lapply(p, function(q) vapply(q, nrow, integer(1)))
    if (!identical(rows(a), rows(b))) {
      out <- c(out, sprintf("%s: the geometries have other numbers of vertices", nm))
      next
    }
    a <- unlist(a, recursive = FALSE); b <- unlist(b, recursive = FALSE)
    if (!identical(lapply(a, is.na), lapply(b, is.na))) {
      out <- c(out, sprintf("%s: other geometries are empty", nm))
      next
    }
    if (!length(b)) next
    all_b <- do.call(rbind, b)
    if (!any(is.finite(all_b))) next
    extent <- max(diff(range(all_b[, 1L], na.rm = TRUE)), diff(range(all_b[, 2L], na.rm = TRUE)))
    gap <- max(mapply(core_frozen_gap, a, b))
    if (!is.finite(gap) || gap > 1e-8 * max(extent, 1))
      out <- c(out, sprintf("%s: a coordinate differs by %.3g, more than 1e-8 of the extent (%.3g)",
                            nm, gap, extent))
  }
  out
}

# Text with what a platform can move taken out: every number is `#`, NaN is
# NA (R does not say which of the two a computation on missing values
# returns), and runs of blanks are one.  A changed word, a line more or a
# line fewer is left to see.
core_frozen_mask <- function(x) {
  if (is.null(x)) return(NULL)
  x <- gsub("[-+]?[0-9][0-9.,]*(?:[eE][-+]?[0-9]+)?", "#", x, perl = TRUE)
  x <- gsub("\\bNaN\\b", "NA", x, perl = TRUE)
  gsub("[[:space:]]+", " ", trimws(x))
}

# The differences between two records under the rule for another platform,
# as lines of text; none when they agree.  Identical: names, classes,
# attribute names and values that are not numbers; when each condition was
# raised, its kind, its class and the names of its fields; the number of
# lines logged, written and printed, and their text with the numbers taken
# out (core_frozen_mask()); whether the call left the random-number stream
# as it found it.  Numbers, in the value and in a condition's fields, equal
# to 1e-12 as all.equal() measures it (the mean difference over the mean
# size).  Coordinates equal to 1e-8 of the extent of the geometry set.
# Not compared: NA against NaN among the numbers, the function a condition's
# call names, and how far a call that draws moves the stream.
core_frozen_differences <- function(got, want) {
  out <- core_frozen_compare(got$value, want$value, "value")
  sig <- function(r) lapply(r$conditions, function(c) c[c("when", "kind", "class")])
  if (!identical(sig(got), sig(want))) {
    said <- function(r) paste(vapply(r$conditions, function(c)
      sprintf("%s (%s, in the %s)", c$class[1L], c$kind, c$when), character(1)), collapse = ", ")
    out <- c(out, sprintf("the conditions differ: now %s; on the baseline %s",
                          if (length(got$conditions)) said(got) else "none",
                          if (length(want$conditions)) said(want) else "none"))
  } else {
    for (i in seq_along(want$conditions)) {
      a <- got$conditions[[i]]; b <- want$conditions[[i]]
      label <- sprintf("condition %d (%s)", i, b$class[1L])
      out <- c(out, core_frozen_compare(a$fields, b$fields, paste("the fields of", label)))
      if (!identical(core_frozen_mask(a$message), core_frozen_mask(b$message)))
        out <- c(out, sprintf("the text of %s differs beyond its numbers: now \"%s\"; on the baseline \"%s\"",
                              label, a$message, b$message))
    }
  }
  texts <- c(log = "logged lines", output = "lines written to standard output",
             printed = "printed lines")
  for (part in names(texts)) {
    a <- core_frozen_mask(got[[part]]); b <- core_frozen_mask(want[[part]])
    if (identical(a, b)) next
    if (is.null(a) || is.null(b)) {
      out <- c(out, "one of the two has a print and the other none")
    } else if (length(a) != length(b)) {
      out <- c(out, sprintf("%d %s; on the baseline %d", length(a), texts[[part]], length(b)))
    } else {
      i <- which(a != b)[1L]
      out <- c(out, sprintf("the %s differ beyond their numbers, first at line %d: now \"%s\"; on the baseline \"%s\"",
                            texts[[part]], i, got[[part]][i], want[[part]][i]))
    }
  }
  if (!identical(got$rng$untouched, want$rng$untouched))
    out <- c(out, sprintf("the random-number stream: now %s; on the baseline %s",
                          if (isTRUE(got$rng$untouched)) "left as it was" else "moved",
                          if (isTRUE(want$rng$untouched)) "left as it was" else "moved"))
  out
}
