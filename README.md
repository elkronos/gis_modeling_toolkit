# spatialkit

> Spatial tessellation, modeling, and cross-validation toolkit for R

[![R-CMD-check](https://github.com/elkronos/gis_modeling_toolkit/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/elkronos/gis_modeling_toolkit/actions/workflows/R-CMD-check.yaml)
[![R >= 4.1](https://img.shields.io/badge/R-%3E%3D%204.1-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://github.com/elkronos/gis_modeling_toolkit/blob/main/LICENSE.md)

## The problem this solves

You have point observations — sites, sensors, parcels, plots — and you want a
model that predicts well at locations you have not sampled. You fit something
flexible, cross-validate it, and get an R² of 0.8. Then you predict onto a grid
and the map is wrong in ways the score never hinted at.

The score was wrong because the folds were wrong. Random k-fold puts a test
point's nearest neighbours in the training set, and under spatial
autocorrelation those neighbours carry most of its signal. The model is scored
on interpolation between known points while the actual task is extrapolation
away from them. The same thing happens one level up, in variable selection: a
predictor chosen by random inner folds can be chosen for leaking rather than
for predicting.

`spatialkit` is built around fixing that. It provides folds that hold out whole
regions sized to the data's own autocorrelation range; three model backends —
geographically weighted regression (`GWmodel`), a Bayesian spatial Gaussian
process (`brms`), and random forests (`ranger`) — behind one `spatial_fit` S3
class so they can be scored on identical folds; and an area-of-applicability
estimate that marks where the resulting score actually applies. Around all of
that sits a tessellation and aggregation pipeline built on
[`sf`](https://r-spatial.github.io/sf/): CRS management, Voronoi, hex, square
and Delaunay cells with stable reproducible IDs, and cell-level summaries with
optional design-effect-corrected standard errors.

Here is the gap, on 600 synthetic points with a smooth spatial field the
predictor does not explain. Same data, same model, same settings — only the
fold construction differs:

| Fold scheme | Pooled out-of-sample R² | RMSE |
|---|---:|---:|
| `make_folds(method = "random_kfold")` (the default) | **0.824** | 1.075 |
| `make_folds(method = "block_kfold")` | **0.379** | 2.129 |

```r
library(spatialkit)
library(sf)

set.seed(42)
n     <- 600
x     <- runif(n, 0, 10000)
y     <- runif(n, 0, 10000)
field <- 3 * sin(x / 1500) * cos(y / 1500) + 2 * sin((x + y) / 4000)
site  <- st_as_sf(data.frame(elev = rnorm(n), x = x, y = y),
                  coords = c("x", "y"), crs = 32632)
site$price <- 10 + 1.5 * site$elev + field + rnorm(n, 0, 0.5)

random  <- make_folds(site, k = 5, method = "random_kfold", seed = 42)
blocked <- make_folds(site, k = 5, method = "block_kfold",  seed = 42)

# include_coords = TRUE lets the forest memorise location. That is the point:
# it is the failure mode the default (FALSE) exists to prevent, and the one
# random folds cannot see.
args <- list(include_coords = TRUE, num_trees = 300, seed = 1)
do.call(cv_rf, c(list(site, "price", "elev", folds = random),  args))$overall
do.call(cv_rf, c(list(site, "price", "elev", folds = blocked), args))$overall
```

More than half the reported skill was leakage. The gap is largest for models
that can interpolate location directly — a forest handed the coordinates, a
GWR, a Gaussian process — but it does not disappear when you leave the
coordinates out. It just stops being visible unless you measure it.

`make_folds()`'s default is `random_kfold`, deliberately: it is the baseline
you compare against, not the one to report.

## Not on CRAN

`spatialkit` is not on CRAN. `install.packages("spatialkit")` will not find it.

```r
# From GitHub
# install.packages("remotes")
remotes::install_github("elkronos/gis_modeling_toolkit")

# With the vignette (needs pandoc)
remotes::install_github("elkronos/gis_modeling_toolkit", build_vignettes = TRUE)

# From a local clone
devtools::install("path/to/spatialkit")
```

Hard dependencies (`sf`, `dplyr`, `logger`, `digest`) install automatically.
Everything else is in `Suggests` and checked at runtime:

| Feature | Install |
|---|---|
| GWR modeling | `install.packages(c("GWmodel", "sp"))` |
| Bayesian GP modeling | `install.packages("brms")` (with `rstan`; or the non-CRAN `cmdstanr` backend if you have it) |
| Random forest modeling | `install.packages("ranger")` |
| Delaunay triangulation | `install.packages("geometry")` |
| Variogram autocorrelation range | `install.packages("gstat")` |
| PSIS-LOO information criterion | `install.packages("loo")` |
| Plotting | `install.packages("ggplot2")` (`patchwork` for multi-panel layouts) |
| k-NN weights for Moran's I above n = 5,000 | `install.packages(c("FNN", "Matrix"))` |
| Tibble output from `make_folds()` / CV predictions | `install.packages("tibble")` |

`FNN` and `Matrix` are not merely a speedup. `residual_morans_i()` **errors**
above n = 5,000 without both: `FNN` supplies the kd-tree neighbour lookup and
`Matrix` keeps the weight matrix sparse, and the fallback path allocates a
dense n × n matrix. Below that threshold the fallback runs fine.

## Quick start

This first block uses nothing outside the hard dependencies. It runs on a bare
install.

```r
library(spatialkit)
library(sf)

set.seed(1)
n  <- 400
df <- data.frame(x = runif(n, 0, 5000), y = runif(n, 0, 5000))
df$elev  <- rnorm(n)
df$price <- 50 + 0.004 * df$x + 3 * df$elev + rnorm(n)
pts <- st_as_sf(df, coords = c("x", "y"), crs = 32632)

# Voronoi cells grow from SEEDS, not from the observations themselves.
# Seeding one cell per point is a nearest-neighbour interpolation, not an
# aggregation: every cell holds one observation, so every ..sd_* / ..se_*
# column comes back NA and there is no within-cell variation to estimate an
# intra-class correlation from.
seeds <- get_voronoi_seeds(sample_points = pts, method = "kmeans", n = 25)
tess  <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)

assigned <- assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id")
cells    <- summarize_by_cell(assigned, response_var = "price",
                              predictor_vars = "elev", deff = "kish")

head(as.data.frame(cells)[, c("cell_id", "n", "resp_mean_price",
                              "..sd_resp_price", "..se_resp_price")])
#>   cell_id  n resp_mean_price ..sd_resp_price ..se_resp_price
#> 1       1  7        51.44352        4.908403        4.204213
#> 2       2  4        51.60967        5.925025        5.188861
#> 3       3 11        53.47350        4.165362        3.528304

attr(cells, "deff_applied")$method
#> [1] "kish"
```

The second block fits and validates a model, on the `site` data built at the
top of this file. **It needs `ranger`**, and `gstat` if you want the fold size
derived from the data rather than chosen by hand. Swap
`fit_rf_model()`/`cv_rf()` for `fit_gwr_model()`/`cv_gwr()` (`GWmodel` + `sp`)
or `fit_bayesian_spatial_model()`/`cv_bayes()` (`brms`) — the surrounding code
is identical.

```r
# --- How big should a block be? Ask the data first --------------------------
# estimate_sac_range() fits a variogram to the OLS residuals and reports the
# distance at which spatial correlation dies out. It refuses to guess:
range <- estimate_sac_range(site, response_var = "price", predictor_vars = "elev")
range
#> NA
#> attr(,"rejected_range")   [1] 98348.32
#> attr(,"rejected_reason")  [1] "fitted range exceeds the largest lag fitted"
# ...and a logged WARN spelling out that the fitted range (98348) exceeds the
# largest lag the variogram covers (7043), so it is unidentified, not long.

# That NA is the correct answer, not a failure: this synthetic field is a smooth
# sinusoid that never levels off inside the study area, so its range is
# unidentified rather than long. `auto_range = TRUE` would ask the same question,
# log the same explanation, and fall back to geometric blocks -- so here we size
# the blocks ourselves and say so.
folds <- make_folds(site, k = 5, method = "block_kfold",
                    block_size = 2000, seed = 42)
folds$params$block_size   # what actually got used -- always worth checking
#> [1] 2000

# --- Fit, then cross-validate on those exact folds --------------------------
fit <- fit_rf_model(site, response_var = "price", predictor_vars = "elev")
cv  <- cv_rf(site, "price", "elev", folds = folds)
cv$overall        # pooled out-of-sample RMSE, MAE, R2, ...
cv$fold_metrics   # per-fold breakdown

# --- Residual spatial autocorrelation diagnostic ----------------------------
residual_morans_i(fit)
#> $observed [1] 0.7023314   $z [1] 36.3353   $p_value [1] 4.48e-289
```

A large positive residual Moran's I says the model has left spatial structure
on the table — here, because `elev` is unrelated to the spatial field.

Compare backends on the same folds. Pass `folds` explicitly: without it each
wrapper builds its own, and the comparison is between models *and* fold
schemes:

```r
comparison <- compare_models_cv(site, "price", "elev",
                                models = c("RF", "GWR", "Bayesian"),
                                folds = folds)
comparison$overall
```

`models` accepts any subset of `c("GWR", "Bayesian", "RF")` in any order.
An unrecognised name raises a `warning()` and is dropped; if nothing
recognised remains, that is an error rather than a silent fallback. A
recognised model whose backend package is not installed is dropped with a
message, so the call still returns the models that could run. Backend-specific
arguments go in `gwr_args`, `bayes_args` and `rf_args`.

## Which function do I want?

| I want to… | Use |
|---|---|
| cut my study area into cells | `build_tessellation()` (`"voronoi"`, `"hex"`, `"square"`, `"triangles"`) |
| choose how many cells | `determine_optimal_levels()` |
| put points into cells and aggregate | `assign_features_to_polygons()` → `summarize_by_cell()` |
| know how far spatial correlation reaches | `estimate_sac_range()` |
| build honest CV folds | `make_folds()` |
| see whether my folds actually separate | `plot_folds()` |
| fit a model | `fit_gwr_model()`, `fit_bayesian_spatial_model()`, `fit_rf_model()` |
| score it out of sample | `cv_gwr()`, `cv_bayes()`, `cv_rf()` |
| score *my own* learner on the same folds | `cv_spatial()` + `new_spatial_fit()` |
| pick predictors without leaking | `select_features_forward()` |
| compare backends head to head | `compare_models_cv()` |
| check for leftover spatial structure | `residual_morans_i()`, `plot(fit, type = "variogram")` |
| turn a fit into a map | `predict_surface()` |
| know where that map is extrapolation | `area_of_applicability()` |
| free memory held by the caches | `clear_grid_cache()`, `clear_fitted_cache()` |

## Function overview

| Area | Key functions |
|---|---|
| CRS & geometry | `ensure_projected()`, `harmonize_crs()`, `coerce_to_points()`, `prep_model_data()` |
| Tessellation | `build_tessellation()`, `create_voronoi_polygons()`, `create_grid_polygons()`, `create_grid_polygons_cached()`, `clip_target_for()`, `ensure_stable_poly_id()` |
| Seeding & resolution | `get_voronoi_seeds()`, `voronoi_seeds_kmeans()`, `voronoi_seeds_random()`, `determine_optimal_levels()` |
| Assignment & aggregation | `assign_features_to_polygons()`, `summarize_by_cell()` |
| Modeling | `fit_gwr_model()`, `fit_bayesian_spatial_model()`, `fit_rf_model()`, `gp_lengthscale_bounds()`, `new_spatial_fit()`; S3: `predict()`, `fitted()`, `residuals()`, `coef()`, `summary()`, `model_metrics()` |
| Cross-validation | `make_folds()` (`random_kfold`, `block_kfold`, `buffered_loo`, `leave_location_out`, `nndm`), `estimate_sac_range()`, `cv_gwr()`, `cv_bayes()`, `cv_rf()`, `cv_spatial()` |
| Variable selection | `select_features_forward()`, `gwr_model_selection()` |
| Prediction | `predict_surface()` |
| Comparison & diagnostics | `compare_models()`, `compare_models_cv()`, `evaluate_insample()`, `residual_morans_i()`, `area_of_applicability()` |
| Plotting | `plot_tessellation_map()`, `plot()` for `spatial_fit`, `plot_folds()` |
| Cache management | `clear_grid_cache()`, `clear_fitted_cache()` |

### The `spatial_fit` interface

All three backends return an object inheriting from `spatial_fit`, so
`predict()`, `fitted()`, `residuals()`, `summary()` and `model_metrics()` work
the same way on all of them. Two caveats are real and worth stating up front:

- **`coef()` is not universal.** `coef.rf_fit()` **errors by design** — a
  forest has no coefficients. Use `fit$info$importance` (permutation
  importance) instead. `coef.gwr_fit()` returns the local coefficient frame and
  `coef.bayesian_fit()` the fixed-effect posterior summary; both `stop()` when
  the backend cannot supply them rather than returning `NULL`, so
  `lapply(fits, coef)` never silently returns a short answer. Wrap in `try()`
  when sweeping a heterogeneous list.
- **`summary()` does not mean the same thing on every backend.**
  `fitted.rf_fit()` returns **out-of-bag** predictions, not in-sample ones, so
  `summary()` on an `rf_fit` reports out-of-bag metrics (and says so). A
  `gwr_fit` and a `bayesian_fit` report genuinely in-sample metrics. The two
  are not comparable. Use `compare_models_cv()`, which scores every backend the
  same way. The flag is `fit$info$fitted_are_oob`.

### Extending it: your own learner on the same folds

`cv_spatial()` is the extensibility point. Pass any `fit_fn(train_sf)` that
returns a `spatial_fit`, built with the `new_spatial_fit()` constructor, and it
plugs into the same folds, metrics and comparison tooling:

```r
lm_fit <- function(train_sf) {
  new_spatial_fit(
    subclass       = "lm_fit",
    engine         = lm(price ~ elev, sf::st_drop_geometry(train_sf)),
    formula        = price ~ elev,
    response_var   = "price",
    predictor_vars = "elev",
    data_sf        = train_sf
  )
}
predict.lm_fit <- function(object, newdata = NULL, ...)
  as.numeric(stats::predict(object$engine, sf::st_drop_geometry(newdata)))

cv <- cv_spatial(site, "price", "elev", fit_fn = lm_fit, folds = blocked)
cv$overall
cv$n_folds_attempted  # compare against n_folds_succeeded before trusting the above
```

`fit_fn` is called on the training slice only, so anything you do inside it —
scaling, tuning, an inner feature sweep — is already nested and leak-free.

## Tessellation

Three ways to cut the same geography. The first panel shows raw observations of
a smooth spatial field over North Carolina (high in the west, with an eastern
hotspot). The others show the same points aggregated into Voronoi regions grown
from k-means seeds, a hex grid, and a square grid — each cell coloured by the
mean of the points inside it (`assign_features_to_polygons()` +
`summarize_by_cell()`). All panels share one colour scale, so each tessellation
should look like a mosaic version of the raw data:

![Raw spatial field and three tessellation methods aggregating it over North Carolina](https://raw.githubusercontent.com/elkronos/gis_modeling_toolkit/main/man/figures/readme-tessellations.png)

All four methods return `cells` carrying a `cell_id` column, so downstream code
and `plot_tessellation_map(fill_col = "cell_id")` treat them alike; hex and
square grids additionally keep `poly_id`, holding the same values.

### How many cells?

Resolution is a modeling decision, not a cosmetic one: too few cells smooth the
signal away, too many leave each cell with a handful of noisy observations.
`determine_optimal_levels()` combines two criteria — a geometric elbow on the
within-cluster sum of squares from a k-means sweep, and a model-aware check
computing Moran's I on OLS residuals at each candidate resolution. The Moran's
I profile measures how much spatial structure in the response remains
*unexplained* at that scale, so the combined ranking balances parsimony against
residual spatial independence.

It returns a **vector of `top_n` candidates** (default 3), best first, with a
`"diagnostics"` attribute — not a single number:

```r
k <- determine_optimal_levels(pts, response_var = "price",
                              predictor_vars = "elev", max_levels = 15,
                              criterion = "combined")
k
#> [1] 4 5 6
k[1]                          # the top-ranked candidate
#> [1] 4
attr(k, "diagnostics")$moran_i
```

The same field cut at three resolutions, next to the raw observations — too
coarse blurs the hotspot, too fine chases noise with near-empty cells, and the
selected k preserves the trend without overfitting geography:

![Raw observations and the same spatial field tessellated at three resolutions, including the automatically selected one](https://raw.githubusercontent.com/elkronos/gis_modeling_toolkit/main/man/figures/readme-resolution.png)

## Cross-validation

`make_folds()` offers five schemes:

| Method | Holds out | Use when |
|---|---|---|
| `random_kfold` (default) | a random 1/k of rows | you want the optimistic baseline to compare against |
| `block_kfold` | contiguous spatial blocks | the general case; size blocks from `estimate_sac_range()` via `auto_range = TRUE` |
| `buffered_loo` | one point, plus everything within `buffer` | leave-one-out with an explicit exclusion radius |
| `leave_location_out` | every observation from one site (`group_var`) | repeated measurements at the same locations |
| `nndm` | one point, plus a distance-matched exclusion | you know where you will predict; pass `prediction_points` |

`nndm` implements Milà et al. (2022): rather than picking a `buffer` with
nothing to justify it, the exclusion around each held-out point is sized so the
training-to-test distance distribution reproduces the distances from your
actual prediction locations to the training data. `params$target_median` and
`params$realised_median` record how close the match came.

The difference is easy to see. On the same North Carolina sites, random folds
scatter test points among their spatially correlated training neighbours, while
block folds hold out contiguous regions:

![Random k-fold versus spatial block k-fold assignments across North Carolina](https://raw.githubusercontent.com/elkronos/gis_modeling_toolkit/main/man/figures/readme-spatial-cv.png)

`make_folds()` may return fewer folds than requested — `block_kfold` when the
grid yields fewer non-empty blocks, `leave_location_out` when there are fewer
distinct groups, and `buffered_loo`/`nndm` always return `k = n`. Read
`folds$k` rather than assuming.

## Prediction

All three backends support true out-of-sample prediction — `newdata` needs only
the predictor columns and geometry, not the response:

```r
new_sites <- st_sf(
  elev = c(0.5, -1.2),
  geometry = st_sfc(st_point(c(1000, 2000)), st_point(c(3000, 4000)), crs = 32632)
)
predict(fit, newdata = new_sites)
```

`newdata` is automatically transformed to the CRS used during fitting, rows with
missing or non-finite predictors return `NA` (output length always matches
`nrow(newdata)`), and non-point geometries are coerced to representative points.
A length mismatch between the predictions and the rows that survived cleaning is
an error, never a silent recycle. The Bayesian method additionally supports
`type = "predict"` for full posterior predictive draws and `draws = TRUE` for
the raw draw matrix; both are honoured with `newdata = NULL` as well.

`model_metrics()` *does* require the response in `newdata`, since it computes
error metrics against observed values.

### Prediction surfaces

Building `newdata` by hand is the fiddly part of producing the thing most
people actually want from a fitted spatial model — a map. `predict_surface()`
builds a regular grid over the training extent, joins covariates from the
nearest feature, clips to a boundary, predicts in chunks and returns `sf`:

```r
surf <- predict_surface(fit, n_cells = 2000, covariates = site)
plot(surf[".pred"])
```

Chunking matters for `bayesian_fit`, where the posterior draw matrix is
`n_draws x n_newdata` and a fine grid would exhaust memory long before the fit
itself would. Pass `se = TRUE` for a posterior-SD surface where the backend
exposes draws, and `boundary =` a polygon to clip to a study area.

## Parallel cross-validation

Every CV function — `cv_gwr()`, `cv_bayes()`, `cv_rf()` and `cv_spatial()` —
accepts a `parallel` argument for fold-level parallelism via
`parallel::mclapply()` (macOS/Linux; falls back to sequential on Windows with a
message). This matters most for `cv_bayes()`, where every fold is a full MCMC
run:

```r
cv <- cv_bayes(site, "price", "elev", k = 5, parallel = TRUE)  # auto-detect cores
cv <- cv_rf(site, "price", "elev", k = 5, parallel = 4L)       # explicit count
```

Results are reproducible from `seed` and identical to `parallel = FALSE`: one
RNG stream per fold is drawn in the parent process, so each fold's stream is a
function of `(seed, fold index)` alone.

## Caching

Grid construction and posterior expectations are both expensive enough to
memoise, so both are cached:

- `create_grid_polygons_cached()` memoises grids keyed on boundary geometry,
  CRS, target cell count and arguments. `clear_grid_cache()` empties it.
- `fitted()` on a `bayesian_fit` memoises `posterior_epred()` column means in
  an environment carried on the object, because `summary()`, `residuals()`,
  `model_metrics()` and `compare_models()` each call `fitted()` independently.
  `clear_fitted_cache(fit)` drops it — needed only if you mutate the engine or
  the training data by hand after fitting.

## Diagnostics & statistical notes

**Residual Moran's I.** `residual_morans_i()` computes Moran's I on model
residuals with the Cliff & Ord randomisation variance, using row-standardised
k-NN weights by default (sparse via `FNN` + `Matrix` when available) or a
user-supplied weight matrix (base or sparse `Matrix`). `compare_models()` runs
it automatically and logs a warning when residual spatial structure remains.

**Aggregation standard errors are IID unless you ask otherwise.** The `..se_*`
columns from `summarize_by_cell()` are IID standard errors at the default
`deff = 1`, which is anticonservative under within-cell spatial correlation.
Autocorrelation-aware standard errors are **opt-in**: pass `deff = "kish"` to
apply Kish's design-effect correction from estimated intra-class correlations
(separate ICCs for response and predictors), a fixed numeric design effect, or
`deff = "variogram"` to compute it from a fitted variogram instead of one
pooled correlation. Kish assumes every pair in a cell is equally correlated
regardless of separation, an assumption that degrades as cells grow; the
variogram option lets correlation decay with distance, which is what having
fitted a variogram is for. Substituting a constant off-diagonal correlation
recovers Kish exactly. Inspect what was applied via
`attr(result, "deff_applied")`.

**Area of applicability.** A fitted model returns a number for any location you
hand it, including locations whose predictor values look nothing like anything
it was trained on. Those predictions are extrapolations dressed as
interpolations, and a cross-validation score says nothing about them — the
held-out folds were drawn from the same predictor distribution as the training
data. `area_of_applicability()` implements the dissimilarity index of Meyer &
Pebesma (2021) and marks where the score applies:

```r
surf <- predict_surface(fit, n_cells = 2000, covariates = site)
aoa  <- area_of_applicability(surf, model = fit, folds = folds)

# AOA is NA wherever a predictor was NA or non-finite. `!NA` is NA, and R
# silently skips NA positions in a subscripted assignment, so test for TRUE
# explicitly and let anything else count as outside.
inside <- aoa$aoa$AOA %in% TRUE
surf$.pred[!inside] <- NA          # blank out the extrapolations
```

Pass the `make_folds()` result you actually validated with. Without it the
reference distance is each training point's nearest neighbour anywhere in the
data, which for clustered data is very close, giving a conservative area; with
it the reference distances are larger and the area is correspondingly wider.
That is not a loophole — the area of applicability is defined relative to a
performance estimate, and a spatially blocked estimate is a claim about
predicting further away.

**Random forests and location.** `fit_rf_model()` defaults to
`include_coords = FALSE`. Handing a forest the x and y coordinates lets it
reproduce the training surface almost exactly by memorising location, then fail
badly anywhere it has not seen; random cross-validation does not catch this,
because nearby points leak between folds (Meyer et al. 2019). That is the
effect measured in the table at the top of this file. Relatedly, `fitted()` on
an `rf_fit` returns **out-of-bag** predictions rather than in-sample ones —
in-sample predictions from a forest are close to memorisation and would make
`summary()` report a fictitious R². The out-of-bag error is itself a *random*
hold-out, so it is optimistic under spatial autocorrelation for the same reason
random k-fold is; use `cv_rf()` for a blocked estimate.

**Variable selection.** `select_features_forward()` scores candidates against
spatially blocked inner folds, which is the entire point of having it: random
inner folds inside blocked outer folds select variables that look predictive
only because nearby points leak between train and test, and the outer loop then
reports honest-looking numbers for a dishonestly chosen feature set. Every
candidate set is scored on the same rows, so a variable cannot win by having an
easier surviving subset. `gwr_model_selection()` is the fast in-sample
counterpart — the same forward search scored by AICc — and is worth
cross-checking against the blocked estimate when the answer matters.

**GWR collinearity.** `fit_gwr_model()` checks the global condition number of
the predictor matrix *and* spot-checks local condition numbers within bandwidth
windows at sampled locations, since spatially clustered subsets can be
collinear even when the global matrix is not.

**GWR bandwidth fallback.** If automatic bandwidth selection fails, a heuristic
fallback is used, a `warning()` is raised, and
`fit$info$bandwidth_is_fallback = TRUE` is set so downstream comparisons can
flag the result. Supply an explicit `bandwidth` if you see this.

**Bayesian GP anisotropy.** `fit_bayesian_spatial_model()` standardizes X and Y
coordinates independently before the GP term, as a conditioning step — easting
and northing often span very different ranges in a projected CRS. Because the
axes are scaled separately, a single shared length-scale would make the kernel
anisotropic in the original CRS by the arbitrary ratio `sd(X)/sd(Y)`, which
reflects the sampling layout rather than the process. The GP therefore fits one
length-scale per axis by default (`gp_iso = FALSE`), estimating directional
structure from the data; pass `gp_iso = TRUE` for a single shared length-scale.
The scaling strategy is recorded in `fit$info$coord_scaling$scaling_type`, and
a data-informed length-scale prior is derived from the inter-point distance
distribution (see `gp_lengthscale_bounds()`). `predictor_vars = character(0)`
is accepted, giving an intercept-only spatial GP — the natural null model for
asking how much of a surface is spatial structure rather than covariate effect.

## What spatialkit does not do

- **Regression only.** A non-numeric response is refused outright, and a binary
  numeric response is refused by `fit_gwr_model()` with a pointer to
  `GWmodel::ggwr.basic()`. There is no classification path.
- **Vector point data only.** No raster support. `predict_surface()` returns an
  `sf` POINT layer, not a `SpatRaster`; convert downstream if you need one.
- **No areal / lattice models.** Moran's I here is a *residual diagnostic on
  point data*. There is no CAR, SAR or spatial-lag fitter.
- **No spatio-temporal folds.** Every fold scheme is purely spatial.

## Logging

Detailed diagnostics are logged to a session temp file, and warnings are echoed
to the console. Logging is scoped to the `"spatialkit"` namespace and never
touches your global logger configuration. To customize:

```r
logger::log_appender(logger::appender_file("my_analysis.log"), namespace = "spatialkit")
logger::log_threshold(logger::WARN, namespace = "spatialkit")
```

Note that these are `logger` messages, not R conditions: `tryCatch(warning = )`
will not catch them and `suppressWarnings()` will not suppress them. Where the
documentation says a function *raises* a warning, it means a genuine R
`warning()`; where it says a function *logs* one, it means this.

## Documentation

- Every exported function is documented: see `?fit_gwr_model`, `?make_folds`,
  `?area_of_applicability`, etc.
- A worked end-to-end demo on the North Carolina boundary shipped with `sf`
  runs as a vignette: `vignette("spatialkit_nc_demo")` after installing with
  `build_vignettes = TRUE`.
- A runnable script version is installed with the package:
  `system.file("scripts", "example_nc_demo.R", package = "spatialkit")`.

## Development

```r
devtools::load_all()   # interactive development
devtools::test()       # run the test suite
devtools::document()   # regenerate NAMESPACE + man/ from roxygen2 (7.3.1)
devtools::check()      # full R CMD check (vignette build requires pandoc)
```

The checked-in `NAMESPACE` and `man/` are generated by roxygen2 7.3.1
(`RoxygenNote` in `DESCRIPTION`); `devtools::document()` reproduces them.

The README figures are generated from actual package output; regenerate them
with `Rscript dev/make_readme_figures.R`.

The test suite covers the geometry/tessellation pipeline, every exported
function, and targeted regression tests for the statistical internals (Moran's
I variance and its sparse-weights path, CV fold/row-ID alignment, CRPS,
hex-grid sizing, CRS selection at wide extents, prediction CRS and NA
alignment, and more). Tests that need an optional backend skip automatically
when it is absent. `ranger` guards the most of them, well ahead of `gstat`,
`ggplot2`, `FNN`, `sp`, `Matrix`, `geometry` and `GWmodel`; exactly one guards
on `brms`. The `backends` job in `R-CMD-check.yaml` installs all of those
except `brms`, so a green matrix means the guarded paths actually ran. `brms`
is left to the weekly `check-brms` workflow, whose header comment is explicit
about how little it currently exercises.

Contributions are welcome — please
[open an issue](https://github.com/elkronos/gis_modeling_toolkit/issues)
describing the bug or proposed change, and include a regression test with any
fix.

## Citation

```r
citation("spatialkit")
```

The package entry cites `spatialkit` itself. `citation()` also prints the
method references you should cite alongside it if you rely on
`area_of_applicability()` (Meyer & Pebesma 2021), `make_folds(method = "nndm")`
(Milà et al. 2022), the Bayesian GP basis (Riutort-Mayol et al. 2023),
permutation importance in `fit_rf_model()` (Strobl et al. 2007), coordinate
predictors and blocked validation (Meyer et al. 2019), or GWR via `GWmodel`
(Lu et al. 2014).

## Maintainer

Justin Chase <jchase.msu@gmail.com> —
[issue tracker](https://github.com/elkronos/gis_modeling_toolkit/issues)

## License

MIT © Justin Chase. See
[LICENSE.md](https://github.com/elkronos/gis_modeling_toolkit/blob/main/LICENSE.md).

## Disclaimer

This is a personal project. It is not affiliated with, endorsed by, or
connected to any organization. It uses public data sources only and was
developed independently on personal time. No confidential, proprietary, or
non-public information is included.
