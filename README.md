# GIS Modeling Toolkit (`gis_modeling_toolkit.R`)

Data-driven spatial partitioning for R: build Voronoi/hex/square tessellations, assign features to cells, and compare models (GWR & Bayesian) to pick the best number of groups—without being locked to arbitrary admin boundaries.

---

## Why
Real-world phenomena (weather, demand, risk) don’t stop at county or ZIP borders. This toolkit lets you **derive boundaries from the data itself** and evaluate which partitioning best explains your outcomes—now with **predictive, cross-validated evaluation**.

---

## What it does
- **CRS-safe preprocessing**: auto-detects lon/lat, picks a sensible projected CRS, and harmonizes layers.  
  - If input has **no CRS** and looks like lon/lat, it’s treated as WGS84 and then projected.  
  - If input has **no CRS and does not** look like lon/lat, CRS is left **unset** (no forced projection).
- **Pointization**: converts lines/polygons to representative points for clustering/assignment. Strategies:  
  `auto`, `centroid`, `point_on_surface` (alias: `surface`), `line_midpoint`, `bbox_center`.  
  Optional `tmp_project=TRUE/FALSE` for safe midpoint sampling on lon/lat lines.
- **Tessellations**: Voronoi (seeds = `kmeans` / `random` / `provided`), hex grids, and square grids—optionally clipped to a boundary (**holes respected**).
- **Assignment**: maps features to polygons and returns per-cell groupings.
- **Level selection**: suggests good cluster counts via an elbow heuristic on within-cluster SSE.
- **Modeling**: fits GWR (`spgwr`) and Bayesian spatial models (`spBayes`) and scores them (AICc/DIC).  
  Bayesian covariance models: `exponential`, `spherical`, `matern`.
- **Predictive evaluation**: cross-validation helpers (`make_folds()`, `cv_gwr()`, `cv_bayes()`), and **fold-aware** `evaluate_models()` when you provide folds.
- **Plotting**: clean maps with polygon ID labels (or counts) over optional basemaps/boundaries.

---

## Requirements
- R (≥ 4.2 recommended)
- Packages: `logger`, `sf`, `sp`, `spgwr`, `spBayes`, `deldir`, `ggplot2`, `dplyr`, `tidyr`, `mvtnorm`
- Optional: `ggspatial` (OSM tiles), `patchwork` (multi-panel figure)

Install in R:
    
```r
install.packages(c(
  "logger","sf","sp","spgwr","spBayes","deldir",
  "ggplot2","dplyr","tidyr","mvtnorm","ggspatial","patchwork"
))
```

## Quick start

1) **Source the toolkit**
```r
source("R/gis_modeling_toolkit.R")
```

2) **Make or load points** (use your own data; here’s a simple demo)
```r
library(sf); library(dplyr)
set.seed(42)
bb  <- sf::st_as_sfc(sf::st_bbox(c(xmin = -80, ymin = 35, xmax = -79, ymax = 36), crs = 4326))
pts <- sf::st_sample(bb, size = 300, type = "random", exact = TRUE) |> sf::st_sf()
pts <- ensure_projected(pts)  # choose a sensible projected CRS
```

3) **Build tessellations** (k = 8 shown; use any of: "voronoi","hex","square")
```r
bt_v <- build_tessellation(pts, levels = 8, method = "voronoi", seeds = "kmeans")
bt_h <- build_tessellation(pts, levels = 8, method = "hex")
bt_s <- build_tessellation(pts, levels = 8, method = "square")
```

4) **Plot** (labels display polygon IDs; or show per-cell counts)
```r
p_v <- plot_tessellation_map(bt_v[["8"]]$polygons, bt_v[["8"]]$data,
                             title = "Voronoi (k=8)", label = "poly_id")
p_h <- plot_tessellation_map(bt_h[["8"]]$polygons, bt_h[["8"]]$data,
                             title = "Hex (≈8 cells)", label = "poly_id")
p_s <- plot_tessellation_map(bt_s[["8"]]$polygons, bt_s[["8"]]$data,
                             title = "Square (≈8 cells)", label = "poly_id")

ggplot2::ggsave("voronoi_k8.png", p_v, width = 8, height = 6, dpi = 200)
```

---

## Overlay on a shapefile (North Carolina example)
```r
library(sf)
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
nc_boundary <- sf::st_union(nc)

# Align CRSs
pts <- ensure_projected(pts, nc_boundary)
nc_boundary <- ensure_projected(nc_boundary)
nc <- sf::st_transform(nc, sf::st_crs(nc_boundary))

# Build tessellations clipped to the state
k <- 8
bt_v <- build_tessellation(pts, levels = k, method = "voronoi",
                           boundary = nc_boundary, seeds = "kmeans")

# Plot with counties as a basemap and boundary outline
p_v <- plot_tessellation_map(bt_v[[as.character(k)]]$polygons,
                             bt_v[[as.character(k)]]$data,
                             title = sprintf("North Carolina — Voronoi (k=%d)", k),
                             boundary = nc_boundary, basemap = nc, label = "poly_id")
```

## Pick the number of groups & best model
```r
# Suggest cluster counts from data
k_suggest <- determine_optimal_levels(pts, max_levels = 12, top_n = 3)

# Suppose you have a response + predictors on the same sf (add your columns to pts)
# pts$resp <- ...
# pts$x1 <- ...; pts$x2 <- ...
res <- evaluate_models(
  data_sf        = pts,                 # must contain columns resp & predictors
  response_var   = "resp",
  predictor_vars = c("x1","x2"),
  levels         = k_suggest,           # or set explicitly, e.g., c(6,8,10)
  tessellation   = c("voronoi","hex","square"),
  boundary       = NULL,                # or a polygon sf to clip
  seeds          = "kmeans",
  models         = c("GWR","Bayesian"),
  n.samples      = 2000,
  cov_model      = "exponential"        # also: "spherical", "matern"
)

# Results tibble: Tessellation, Level, Model, Metric (AICc for GWR, DIC for Bayesian)
dplyr::arrange(res$results, Metric)
```

---

## Make a single figure with all three maps (white background)
```r
library(patchwork)
combo <- p_v + p_h + p_s + patchwork::plot_layout(ncol = 3) &
         ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA))
ggplot2::ggsave("tessellations_side_by_side.png", combo, width = 18, height = 6, dpi = 200, bg = "white")
```

---

## Function index (brief)

- **`ensure_projected(x, target_crs=NULL)`** — ensures a projected CRS; can force to match a target CRS.  
  - No-CRS lon/lat ⇒ assume WGS84 then project.  
  - No-CRS non–lon/lat ⇒ leave CRS unset.
- **`harmonize_crs(a, b)`** — aligns CRS between two sf/sfc objects (returns list `a`, `b`).
- **`coerce_to_points(x, strategy, tmp_project=TRUE)`** — representative points from any geometry.  
  Strategies: `auto`, `centroid`, `point_on_surface` (alias: `surface`), `line_midpoint`, `bbox_center`.  
  *Note:* `line_midpoint` expects LINESTRING (will error on MULTILINESTRING; cast or use `surface`/`centroid`).  
- **`create_voronoi_polygons(points, clip_with=NULL)`** — Voronoi polygons from seed points; optional clipping.  
  Duplicate/overlapping seeds may share polygons; output preserves one tile per input seed.
- **`create_grid_polygons(boundary, target_cells, type)`** — hex/square grids clipped to a boundary; **holes respected**.
- **`voronoi_seeds_kmeans(points_sf, k)`** — k-means cluster centers as Voronoi seeds (deterministic with fixed seed).
- **`voronoi_seeds_random(boundary, k)`** — k random seeds inside a boundary (set seed for reproducibility).
- **`assign_features_to_polygons(features, polygons)`** — intersects features with tessellation; adds `poly_id`.
- **`determine_optimal_levels(points_sf, max_levels, top_n)`** — elbow heuristic over k-means WSS (stable given RNG seed).
- **`fit_gwr_model(data_sf, response, predictors)`** — fits GWR; returns model, bandwidth ∈ (0,1], AICc.
- **`fit_bayesian_spatial_model(data_sf, response, predictors, n.samples, cov_model)`** — fits Bayesian spatial model; returns samples & DIC.  
  `cov_model` ∈ `{ "exponential","spherical","matern" }` (invalid names error).
- **`build_tessellation(features_sf, levels, method, boundary=NULL, seeds="kmeans", pointize="auto", provided_seed_points=NULL, ...)`** — constructs tessellation(s) and assignments per level.  
  If `seeds="provided"`, you **must** pass `provided_seed_points`.
- **`evaluate_models(...)`** — builds tessellations, fits models across levels/methods; returns a results tibble (AICc for GWR, DIC for Bayesian).
- **`plot_tessellation_map(polygons, points, ..., label=c("poly_id","count","none"), show_counts=FALSE)`** — clean map with labels or count overlays; optional basemap/boundary.
- **`summarize_by_cell(assigned_points_sf, response, predictors)`** — per-cell counts and means.
- **`make_folds(data_sf, method=c("random_kfold","block_kfold","buffered_loo"), k=5, block_size=NULL, buffer_dist=NULL, seed=NULL)`** — creates CV folds; auto-projects for distance ops and preserves original `row_id`.
- **`cv_gwr(data_sf, response, predictors, folds)`** — per-fold GWR with adaptive bandwidth selection; returns RMSE/MAE/R² and bandwidth diagnostics.
- **`cv_bayes(data_sf, response, predictors, folds, n.samples, cov_model)`** — per-fold Bayesian GP (`spBayes::spLM`); returns RMSE/MAE/R² and optional fold DIC, with safe fallbacks if fitting fails.

---

## UAT

A complete UAT script **`run_uat_spatial_modeling.R`** is included. It exercises CRS handling, mixed-geometry pointization, tessellations (including holes in boundaries), assignments, modeling (GWR + Bayesian with `exponential`/`spherical`/`matern`), **cross-validation** (`make_folds()`, `cv_gwr()`, `cv_bayes()`), plotting, and error paths.

**Artifacts written to `UAT_outputs/`:**
- `map_voronoi_kmeans_lvl12.png`, `map_square_lvl20.png`, `map_square_lvl20_nobnd.png`
- `uat_evaluate_models_results.csv`
- `uat_cv_gwr_results.csv`, `uat_cv_bayes_results.csv`  <!-- predictive CV summaries (RMSE/MAE/R²; DIC when available) -->
- `sessionInfo.txt`

Run it with your preferred method (e.g., `Rscript run_uat_spatial_modeling.R`) and adjust `script_path` inside if needed. If `spBayes` is unavailable, the script will skip Bayesian fits gracefully and still write the other artifacts.

---

## Tips & troubleshooting
- **Missing CRS**: set or infer CRS first; `ensure_projected()` assigns WGS84 if the bbox looks like lon/lat, then projects; otherwise leaves CRS unset.
- **CRS mismatches**: align layers with `harmonize_crs(a, b)` or by transforming one to the other's CRS before intersects/joins.
- **Seeds**:
  - If `seeds = "provided"`, you **must** pass `provided_seed_points` and they should be in the **same CRS** as your features/boundary (and, if using a boundary, fall inside it).
  - If `seeds = "random"`, you **must** supply a `boundary` to sample within.
- **Line midpoints**: `line_midpoint` expects **LINESTRING**; for **MULTILINESTRING** either cast to LINESTRING or use `point_on_surface`/`centroid`.
- **Label choice**: set `label = "poly_id"` (IDs) or `label = "count"` (point counts), or use `show_counts = TRUE`.
- **OSM tiles**: set `use_osm_tiles = TRUE` (requires `ggspatial`); all layers are drawn in **EPSG:3857**.
- **Boundaries with holes**: hex/square grids respect holes after clipping; ensure boundary and grid share the same CRS.
- **Geometry validity**: invalid geometries can break spatial ops; check with `sf::st_is_valid()` and repair with `sf::st_make_valid()` before heavy intersects.
- **Non-finite predictors**: model fits will warn/error or drop rows with `NA`/`Inf`—clean inputs for best results.
- **Performance**:
  - Use projected CRSs (not lon/lat) for distance-heavy steps (GWR, buffering, Voronoi).
  - For large `n`, prefer **hex/square** grids, keep Voronoi `k` moderate, and consider subsampling in `determine_optimal_levels()`.
- **Cross-validation**: for **predictive** model selection, prefer `make_folds()` + `cv_gwr()` / `cv_bayes()` (use `block_kfold` or `buffered_loo` to reduce spatial leakage). If `spBayes` is unavailable, Bayesian CV is skipped/falls back gracefully.
- **Reproducibility**: fix RNG seeds (`set.seed(...)`) to stabilize k-means/random seeding, CV allocation, and level suggestions.

---

![NC tessellations](Examples/nc_tessellations_triptych.png)

---

## Current limitations & scaling roadmap

### Current capabilities & limits (updated)

**Where it works well today**
- Interactive → moderate workloads (thousands to low–tens of thousands of features) on a single R session.
- Robust, CRS-safe workflows: automatic local UTM pick (with 3857 fallback), safe pointization for non-point geometries, and Voronoi that correctly handles duplicate/near-duplicate seeds.
- Rapid prototyping of Voronoi/hex/square partitions, assignment, and model comparison (GWR AICc / Bayesian DIC).
- **Predictive CV available:** out-of-fold validation via `make_folds()`, `cv_gwr()`, and `cv_bayes()` with RMSE/MAE/R² summaries; supports `random_kfold`, `block_kfold`, and `buffered_loo`.
- **Hex/square grid tuning:** target a requested cell count; grids respect holes.

**Known limitations**
- **Design note:** `evaluate_models()` remains **in-sample** (AICc/DIC on the assigned training data). For predictive model/tessellation selection, use `make_folds()` + `cv_gwr()` / `cv_bayes()`. A fold-aware tessellation workflow integrated into `evaluate_models()` is on the roadmap.
- **Scaling pressure at large n (≈50k+).**
  - GWR (`spgwr`) relies on dense distance ops and slows with repeated fits (including per-fold CV).
  - Bayesian GP (`spBayes::spLM`) has \(O(n^3)\)-like behavior; memory/time grow quickly and CV multiplies cost. (`cv_bayes()` will fall back to regression-only predictions if `spBayes` is unavailable or a fold fit fails.)
- **Voronoi at high k** remains heavy; tiny deterministic jitter reduces numerical quirks but does not remove the cost.
- **Spatial joins** (`st_within`/`st_intersects`) can bottleneck without careful projection and valid/prepared geometries.
- **Geometry edge cases & projections.** Multi-zone extents, antimeridian crossers, and polar regions may need a hand-picked CRS.  
  For `coerce_to_points(mode = "line_midpoint")`, input must be **LINESTRING** (errors on MULTILINESTRING—cast or use `surface`/`centroid`).

**Cross-validation**
- **Fold generation:** `make_folds()` supports:
  - `random_kfold` (standard),
  - `block_kfold` (regular-grid spatial blocks; balanced across folds),
  - `buffered_loo` (leave-one-out with an exclusion buffer to reduce leakage).
  All modes project lon/lat automatically for distance-based steps and preserve original row indices (`row_id`).
- **GWR CV:** `cv_gwr()` fits per-fold with **adaptive bandwidth selection** (`gwr.sel`) and predicts at test coordinates; returns RMSE/MAE/R² overall and by fold, plus bandwidth diagnostics.
- **Bayesian CV:** `cv_bayes()` fits per-fold with `spBayes::spLM` and predicts with `spPredict()`; returns RMSE/MAE/R² and optional fold DIC when samples are available. Safe fallbacks handle missing packages or unstable fits.
- These CV tools measure **predictive** performance; use them to compare formulas, covariances, or tessellation designs (by reassigning features per design before CV).

**What you can do now**
- For **predictive model selection**, prefer `make_folds()` + `cv_gwr()` / `cv_bayes()` (use `block_kfold` or `buffered_loo` to reduce spatial leakage).
- For **quick, in-sample** comparisons across designs, keep using `evaluate_models()`; you can now set `cov_model` ∈ {`"exponential"`, `"spherical"`, `"matern"`}.
- For large n: start with **hex/square grids**, keep Voronoi k moderate, and use subsampling in `determine_optimal_levels()`.
- Fix RNG seeds for reproducibility (`set.seed(...)` for k-means/random seeding and CV allocation).

