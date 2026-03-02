# spatialkit

A modular R package for spatial tessellation, modeling, and cross-validation.

## Installation

```r
# Install from local source
devtools::install("path/to/spatialkit")

# Or with remotes from a repo
# remotes::install_github("your-org/spatialkit")
```

All hard dependencies (`sf`, `dplyr`, `logger`, `digest`) are installed automatically. Optional backends are declared in `Suggests` and checked at runtime:

| Feature | Requires |
|---|---|
| GWR modeling | `GWmodel`, `sp` |
| Bayesian GP modeling | `brms`, `cmdstanr` (or `rstan`) |
| Delaunay triangulation | `geometry` |
| Variogram-based SAC range | `gstat` |
| LOO cross-validation | `loo` |
| Plotting | `ggplot2` |
| Geometry repair fallback | `lwgeom` |
| Fast k-NN weights (n > 5 000) | `FNN` |

> **Note:** `FNN` is strongly recommended when using `residual_morans_i()` or
> other k-NN weight routines on datasets with more than 5 000 observations.
> Without it the dense fallback computes a full n×n distance matrix, which is
> prohibitively expensive. Install it with `install.packages("FNN")`.

## Quick Start

```r
library(spatialkit)

# Tessellation
tess <- build_tessellation(points_sf, boundary = study_area, method = "voronoi")

# Assign features to cells
assigned <- assign_features_to_polygons(features_sf, tess$cells)

# Fit a GWR model
fit <- fit_gwr_model(data_sf, response_var = "y", predictor_vars = c("x1", "x2"))

# Cross-validate
cv <- cv_gwr(data_sf, "y", c("x1", "x2"), k = 5)

# Cross-validate with parallel fold fitting (macOS / Linux)
cv <- cv_bayes(data_sf, "y", c("x1", "x2"), k = 5, parallel = TRUE)

# Compare models
comparison <- compare_models(list(GWR = gwr_fit, Bayesian = bayes_fit))
```

## Prediction

Both `predict.gwr_fit()` and `predict.bayesian_fit()` support true
out-of-sample prediction where the response variable is unknown.
`newdata` only needs to contain the predictor columns (and geometry):

```r
# newdata does NOT need the response column "y"
new_sites <- sf::st_sf(
  x1 = c(5, 6), x2 = c(7, 8),
  geometry = sf::st_sfc(sf::st_point(c(1, 2)), sf::st_point(c(3, 4)), crs = 32632)
)
preds <- predict(fit, newdata = new_sites)
```

Note that `model_metrics()` *does* require the response in `newdata`
because it needs observed values to compute error metrics.

## Parallel Cross-Validation

By default, CV folds are fitted sequentially. For expensive models —
especially Bayesian GP fits via `cv_bayes()` where each fold involves a
full MCMC run — fold-level parallelism can yield near-linear speedup.
All three CV functions (`cv_gwr()`, `cv_bayes()`, `cv_spatial()`) accept
a `parallel` argument:

```r
# Auto-detect cores (uses all physical cores minus one)
cv <- cv_bayes(data_sf, "y", c("x1", "x2"), k = 5, parallel = TRUE)

# Explicit core count
cv <- cv_gwr(data_sf, "y", c("x1", "x2"), k = 5, parallel = 4L)
```

Parallelism uses `parallel::mclapply()` (fork-based), which works on
**macOS** and **Linux**. On **Windows**, forked parallelism is not
available; the functions fall back to sequential execution with an
informative message. For Windows parallelism, consider the `future` /
`future.apply` ecosystem.

## Diagnostics

### Spatial autocorrelation and `summarize_by_cell()` standard errors

By default, the `..se_*` columns produced by `summarize_by_cell()` are
classic IID standard errors (`sd / sqrt(n)`). Under spatial
autocorrelation — the common case for spatial data — within-cell
observations are positively correlated, so the effective sample size is
smaller than `n`. **The naive SE is anticonservative**, and downstream
weighted regressions using `cell_weight` or `..se_*` will understate
uncertainty for cells with strong intra-cell correlation.

To mitigate this, pass `deff = "kish"` to apply Kish's design-effect
correction. This estimates separate intra-class correlations (ICCs) for
the response and predictor variables and deflates each cell's effective
sample size accordingly. When multiple predictor (or response) columns
are pooled into a single ICC estimate, each column is z-scored first so
that variables on different scales contribute equally:

```r
cell_summary <- summarize_by_cell(
  assigned, response_var = "y", predictor_vars = c("x1", "x2"),
  deff = "kish"   # adjust SE for within-cell spatial autocorrelation
)

# Inspect the estimated per-variable-type ICCs and per-cell design effects
attr(cell_summary, "deff_applied")
# $method   "kish"
# $icc_resp  (ICC for response)
# $icc_pred  (ICC for predictors)
# $deff      (per-cell design effects based on primary ICC)
```

You can also supply a fixed numeric design effect (e.g. `deff = 2`)
when you have an external estimate. Note that even the Kish correction
is a first-order approximation; for rigorous inference under spatial
dependence, fit an explicit spatial covariance model (see
`fit_bayesian_model()`).

### GWR collinearity checking

`fit_gwr_model()` performs both a **global** and a **local** collinearity diagnostic before fitting. The global check computes the condition number of the full predictor matrix. Because GWR fits a separate regression within each bandwidth window, global collinearity can understate the problem: spatially clustered subsets of observations may have near-zero predictor variance even when the full dataset does not. To address this, the function spot-checks the condition number of the local design sub-matrix at up to 30 randomly sampled locations and warns when a substantial fraction exhibit extreme condition numbers (> 1e6). This spot-check is not exhaustive; users with highly clustered data should consider a full local-collinearity audit as a post-fit diagnostic.

### Bandwidth fallback

When automatic bandwidth selection via `bw.gwr()` fails (e.g. due to numerical issues or degenerate data), `fit_gwr_model()` falls back to an arbitrary default bandwidth (50 nearest neighbours for adaptive mode, 1/3 of the spatial diagonal for fixed mode). Because this fallback is not optimized for the data, the resulting fit may be poor. A `warning()` is issued to alert interactive users, and the returned `gwr_fit$info` list includes `bandwidth_is_fallback = TRUE` so that downstream code (including `compare_models()`) can flag the result. Users who encounter a fallback warning should supply an explicit `bandwidth` argument based on domain knowledge or external cross-validation.

### Bayesian GP coordinate scaling (anisotropy note)

`fit_bayesian_spatial_model()` standardizes X and Y coordinates **independently** (each is centered and divided by its own standard deviation) before passing them to the GP term. This per-axis scaling stabilises the GP numerically when coordinate extents differ dramatically between axes (common in projected CRSs), but it means the isotropic squared-exponential kernel used by `brms::gp()` is **anisotropic in the original CRS** whenever `sd(X) ≠ sd(Y)`: the effective spatial correlation length differs between the easting and northing directions.

For many applied use-cases this is benign or even beneficial, but users who require strict isotropy in geographic distance should be aware of this behaviour. The stored `$info$coord_scaling` list includes a `scaling_type` element (`"anisotropic"`) so downstream code can detect the strategy used. To achieve isotropic scaling one could replace the per-axis SD with a single factor such as `max(sd(X), sd(Y))` for both axes.

## Migration Guide (from sourced scripts)

### Before (sourced scripts)

```r
source("main.R")  # dumps ~100 functions into .GlobalEnv
```

### After (package)

```r
library(spatialkit)  # namespaced, documented, dependency-managed
```

### What changed

| Aspect | Before (scripts) | After (package) |
|---|---|---|
| Loading | `source("main.R")` | `library(spatialkit)` |
| Dependencies | Manual `install.packages()` | Auto-resolved via `DESCRIPTION` |
| Namespace | Everything in `.GlobalEnv` | Exported API only; internals hidden |
| Documentation | Comments only | `?function_name` works |
| Testing | None | `testthat` infrastructure ready |
| Logging | Writes `model_log.log` in `getwd()` | Writes to `tempdir()` by default |

### API is unchanged

Every exported function keeps its original name and signature. Code like this works identically:

```r
# These calls are the same before and after migration
seeds <- get_voronoi_seeds(boundary, method = "kmeans", n = 50)
tess  <- create_voronoi_polygons(points_sf, boundary = boundary)
fit   <- fit_gwr_model(data_sf, "price", c("sqft", "beds"))
cv    <- cv_gwr(data_sf, "price", c("sqft", "beds"), k = 5)
```

### Internal functions are now truly internal

Functions prefixed with `.` (e.g., `.is_longlat`, `.safe_make_valid`, `.compute_reg_metrics`) are no longer in your global environment. If you were calling them directly, access them via `spatialkit:::` (not recommended) or use the exported wrappers.

### Logging configuration

Logging is scoped to the `"spatialkit"` namespace and will not interfere with your own logger setup. Reconfigure after loading:

```r
library(spatialkit)
logger::log_appender(logger::appender_file("my_analysis.log"), namespace = "spatialkit")
logger::log_threshold(logger::WARN, namespace = "spatialkit")  # quieter
```

## Package Structure

```
spatialkit/
├── DESCRIPTION          # Dependencies, metadata
├── NAMESPACE            # Exports and S3 method registrations
├── LICENSE
├── R/
│   ├── spatialkit-package.R   # Package-level docs
│   ├── zzz.R                  # .onLoad (logging setup)
│   ├── utils.R                # Shared helpers, metrics
│   ├── crs-geometry.R         # CRS selection, projection, coercion
│   ├── seeding.R              # Voronoi seed generation
│   ├── tessellation.R         # Voronoi, Delaunay, hex/square grids
│   ├── assignment.R           # Feature-to-polygon join
│   ├── level-selection.R      # Elbow heuristic for cluster count
│   ├── model-prep.R           # Data validation, GP length-scale bounds
│   ├── model-classes.R        # S3 class system (spatial_fit)
│   ├── model-gwr.R            # GWR fitting (GWmodel)
│   ├── model-bayesian.R       # Bayesian GP regression (brms)
│   ├── cross-validation.R     # CV infrastructure
│   ├── evaluation.R           # Model comparison and diagnostics
│   ├── plotting.R             # ggplot2 tessellation maps
│   └── stable-ids-cache.R     # Deterministic IDs, grid caching
└── tests/
    ├── testthat.R
    └── testthat/
        └── test-core.R         # Starter test suite
```

## Development

```r
# After cloning the repo:
devtools::load_all()     # Interactive development
devtools::test()         # Run tests
devtools::check()        # Full R CMD check
devtools::document()     # Regenerate NAMESPACE + man/ from roxygen2 tags
```
