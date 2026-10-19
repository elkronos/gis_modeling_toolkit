# What to save, and what to report

The other vignettes each take one decision: how many cells, which fold
scheme, which diagnostic. This one starts after those are settled. You
have regions and a validated model, and two things have to leave your
session: the regions, in a form someone else can open and join to, and
the numbers that say how far the map they are looking at can be trusted.

*Everything here runs on the hard dependencies, with two exceptions that
skip themselves with a note when absent: **ggplot2** draws the one map,
and **gstat** fits the variogram in the second half.*

## The regions as a file someone else can use

North Carolina’s counties, with the 1974 SIDS rate per thousand births,
and a hexagonal grid of regions drawn from where the counties are rather
than from the county boundaries themselves. Twelve cells are asked for;
a grid clipped to an irregular outline returns more, and the count that
comes back is the one to quote.

``` r

library(sf)
library(spatialkit)

nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
counties <- ensure_projected(nc[, c("NAME", "BIR74", "SID74")])
counties$rate <- 1000 * counties$SID74 / counties$BIR74

centres  <- coerce_to_points(counties)
boundary <- clip_target_for(centres, expand = 0.01)
regions  <- build_tessellation(centres, method = "hex", approx_n_cells = 12,
                               boundary = boundary)$cells
nrow(regions)
```

    ## [1] 18

`approx_n_cells` is approximate, as its name says: hexagons over this
bounding box came out at 18 for a request of 12, and some of them will
hold no counties at all.

### Grouping a layer you already have

[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
is usually shown on points, but it takes polygons too. `largest = TRUE`
gives each county to the region it overlaps most, which is what makes
the result a grouping: one row per county, never a county split across
two rows.

``` r

grouped <- assign_features_to_polygons(counties, regions,
                                       polygon_id_col = "cell_id",
                                       largest = TRUE)

c(counties_in = nrow(counties), rows_out = nrow(grouped),
  regions_used = length(unique(stats::na.omit(grouped$cell_id))))
```

    ##  counties_in     rows_out regions_used 
    ##          100          100           12

``` r

range(table(grouped$cell_id))
```

    ## [1]  2 13

100 counties in, 100 rows out, spread over 12 of the 18 regions, from 2
to 13 counties each. That row-count equality is the property to check
after any grouping:
[`sf::st_join()`](https://r-spatial.github.io/sf/reference/st_join.html)
would have duplicated every county that touches two regions, and the
counts per region would then be wrong by however many counties straddle
a border.

[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
reduces the grouped layer to one row per region, and passing `cells_sf`
keeps the geometry and the empty regions:

``` r

cells <- summarize_by_cell(grouped, response_var = "rate", id_col = "cell_id",
                           cells_sf = regions, conf_level = 0.95)

names(cells)
```

    ##  [1] "cell_id"                   "n"                        
    ##  [3] "resp_mean_rate"            "..sd_resp_rate"           
    ##  [5] "..se_resp_rate"            "..neff_resp_rate"         
    ##  [7] "..df_resp_rate"            "..ci_lo_resp_rate"        
    ##  [9] "..ci_hi_resp_rate"         "cell_weight"              
    ## [11] "..pooled_resp_rate"        "..pooled_w_resp_rate"     
    ## [13] "..pooled_target_resp_rate" "geometry"

``` r

sum(is.na(cells$n))
```

    ## [1] 6

The three `..pooled_` columns are each region’s pooled value, its weight
on the region’s own mean and the mean it is pulled toward;
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
adds them by default, `estimate = "mean"` leaves them out, and
“Reporting a cell value” below reads them. 6 of the 18 regions are
empty, carrying `NA` rather than being silently dropped. A region with
no data in it is a fact about the map, and one worth keeping: it is
where the model has nothing to stand on.

``` r

library(ggplot2)

ggplot() +
  geom_sf(data = cells, aes(fill = resp_mean_rate), colour = "white") +
  geom_sf(data = counties, fill = NA, colour = "grey70", linewidth = 0.2) +
  scale_fill_viridis_c(name = "SIDS rate\nper 1,000", na.value = "grey92") +
  labs(title = "Counties grouped into data-drawn regions",
       subtitle = sprintf("%d regions, %d populated; grey = no county inside",
                          nrow(cells), sum(!is.na(cells$n)))) +
  theme_void()
```

![The North Carolina counties grouped into seventeen hexagonal regions,
each shaded by the mean 1974 SIDS rate of the counties inside it. County
outlines show underneath. Five regions along the edges hold no county
and are grey.](reporting_files/figure-html/map-1.png)

### Is this location inside?

The question the regions get asked most often once they leave your
session is whether some new location falls in one. That is the same
function with `keep_unassigned = TRUE`, which returns one row per
location and `NA` for anything outside every region.

``` r

sites <- st_as_sf(
  data.frame(site = c("Raleigh", "Wilmington", "Atlantic", "Norfolk VA"),
             lon  = c(-78.64, -77.95, -74.50, -76.29),
             lat  = c( 35.78,  34.23,  35.00,  36.85)),
  coords = c("lon", "lat"), crs = 4326)

inside <- assign_features_to_polygons(st_transform(sites, st_crs(regions)),
                                      regions, polygon_id_col = "cell_id",
                                      keep_unassigned = TRUE)
st_drop_geometry(inside)
```

    ##         site cell_id
    ## 1    Raleigh      13
    ## 2 Wilmington      14
    ## 3   Atlantic      NA
    ## 4 Norfolk VA      NA

2 of the 4 sit inside a region. The default is
`keep_unassigned = FALSE`, which drops them instead; that is the right
default when you are aggregating, and the wrong one when you are
answering a membership question, because a dropped row looks like a
location you never asked about.

Two details decide whether the answer is meaningful. The locations have
to be in the regions’ CRS, which is why
[`st_transform()`](https://r-spatial.github.io/sf/reference/st_transform.html)
is there: the regions are in whatever projection
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
chose, and lon/lat degrees compared against metres will put everything
outside. And a point on a shared border falls in both cells, so
`tie_break` decides which one it is recorded in; see
[`?assign_features_to_polygons`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
for the rules.

### IDs that survive a reprojection

The `cell_id` a tessellation carries is already derived from the
geometry, so it does not change when the input points arrive in another
order. What it is not is portable: a file that has been reprojected,
re-sorted and re-read elsewhere has no memory of how it was built.
[`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md)
derives an ID from the geometry itself, in a common CRS, so the same
region gets the same ID from whoever computes it, in whatever projection
the file is in by then. The exception is fine cells stacked within about
a centimetre of the same longitude, which a reprojection can swap (see
[`?ensure_stable_poly_id`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md));
regions this size are nowhere near it.

``` r

regions$src <- seq_len(nrow(regions))

from_projected <- ensure_stable_poly_id(regions[, "src"])
from_lonlat    <- ensure_stable_poly_id(st_transform(regions[, "src"], 4326))

identical(from_projected$poly_id[order(from_projected$src)],
          from_lonlat$poly_id[order(from_lonlat$src)])
```

    ## [1] TRUE

Use that column, not the row position, as the key anything downstream
joins on.

### Writing it out

The regions are an ordinary `sf` layer, so
[`sf::st_write()`](https://r-spatial.github.io/sf/reference/st_write.html)
handles them. The format matters more than it looks:

``` r

dir <- tempfile(); dir.create(dir)

st_write(cells, file.path(dir, "regions.gpkg"), quiet = TRUE)
gpkg <- st_read(file.path(dir, "regions.gpkg"), quiet = TRUE)

suppressWarnings(st_write(cells, file.path(dir, "regions.shp"), quiet = TRUE))
shp <- st_read(file.path(dir, "regions.shp"), quiet = TRUE)

attr_names <- function(x) setdiff(names(x), attr(x, "sf_column"))
kept <- function(x) length(intersect(attr_names(x), attr_names(cells)))

c(columns = length(attr_names(cells)),
  GeoPackage = kept(gpkg), Shapefile = kept(shp))
```

    ##    columns GeoPackage  Shapefile 
    ##         13         13          2

    ##  [1] "cell_id"    "n"          "rsp_mn_"    "X__sd_r_"   "X__s_rs_"  
    ##  [6] "X__nff__"   "X__df_r_"   "X__c_l__"   "X__c_h__"   "cll_wgh"   
    ## [11] "X__pld__"   "X__pld_w__" "X__pld_t__"

A shapefile field name is capped at ten characters, so the standard
errors, confidence bounds and pooled values arrive as `X__s_rs_`,
`X__c_l__` and `X__pld__`: of the 13 columns
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
produced, 2 survive the trip. GeoPackage keeps every one of them, and
renames only the geometry column, which is a GeoPackage convention.
Write a GeoPackage unless something downstream can only read shapefiles,
and if it can only read shapefiles, rename the columns yourself first so
you choose the abbreviations.

## The numbers the run already recorded

Nothing below is recomputed for the report. Each of these is written
into the object at the time the decision was made, which is the point: a
report assembled from `$params` and `$info` describes the run that
happened rather than a reconstruction of it.

The counties carry a real response, but a weak one, so the rest of this
vignette runs on a simulated field over the same county centres, with a
known correlation length of 120 km:

``` r

set.seed(42)
xy    <- st_coordinates(centres)
D     <- as.matrix(dist(xy))
field <- as.numeric(t(chol(exp(-D / 120e3) + diag(1e-8, nrow(D)))) %*%
                      rnorm(nrow(D)))

pts <- centres
pts$elev <- rnorm(nrow(pts))
pts$z    <- 0.8 * pts$elev + 2 * field + rnorm(nrow(pts), sd = 0.3)
```

``` r

# A learner with no optional dependencies: a linear trend surface in the
# coordinates, plus the one covariate.
trend_fit <- function(train_sf, ...) {
  d  <- st_drop_geometry(train_sf)
  co <- st_coordinates(train_sf)
  d$x <- co[, 1]; d$y <- co[, 2]
  new_spatial_fit(subclass = "trend_fit",
                  engine = stats::lm(z ~ elev + x + y, data = d),
                  formula = z ~ elev + x + y, response_var = "z",
                  predictor_vars = "elev", data_sf = train_sf)
}
predict.trend_fit <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(stats::fitted(object$engine))
  d  <- st_drop_geometry(newdata)
  co <- st_coordinates(newdata)
  d$x <- co[, 1]; d$y <- co[, 2]
  as.numeric(stats::predict(object$engine, newdata = d))
}
residuals.trend_fit <- function(object, ...) stats::residuals(object$engine)
fitted.trend_fit    <- function(object, ...) stats::fitted(object$engine)
for (.g in c("predict", "residuals", "fitted"))
  registerS3method(.g, "trend_fit", get(paste0(.g, ".trend_fit")))

folds <- make_folds(pts, k = 5, method = "block_kfold", block_size = 100e3,
                    seed = 1)
cv    <- cv_spatial(pts, "z", "elev", trend_fit, folds = folds)
fit   <- trend_fit(pts)
```

## A worked report

Six items. Every one is read out of an object the run already produced.

``` r

m <- function(x) format(round(as.numeric(x)), big.mark = ",", scientific = FALSE)

report <- list(
  folds = sprintf("%s, k = %d, %d x %d blocks of %s m (%s)",
                  folds$method, folds$k, folds$params$grid_nx,
                  folds$params$grid_ny, m(folds$params$block_size),
                  folds$params$crs),
  score = sprintf("RMSE %.2f (per fold %.2f to %.2f), R2 %.3f",
                  cv$overall$RMSE, min(cv$fold_metrics$RMSE),
                  max(cv$fold_metrics$RMSE), cv$overall$R2),
  coverage = sprintf("%d of %d folds scored, %d rows dropped, %d orphan, %d unknown id",
                     cv$n_folds_succeeded, cv$n_folds_attempted, cv$n_dropped,
                     length(cv$orphan_rows), cv$n_unknown_ids)
)
```

``` r

mi <- residual_morans_i(fit)
report$residuals <- sprintf("Moran's I %.3f (z = %.1f, p = %.3g, %s null)",
                            mi$observed, mi$z, mi$p_value, mi$null)
```

``` r

surface <- predict_surface(fit, n_cells = 1500, covariates = pts,
                           boundary = counties)
aoa <- area_of_applicability(surface, model = fit, folds = folds)
report$applicability <- sprintf("%d of %d prediction cells outside (DI > %.3f)",
                                aoa$n_outside, aoa$n_new, aoa$threshold)
```

``` r

sac <- estimate_sac_range(pts, "z")
report$range <- sprintf("estimated range %s m against blocks of %s m",
                        m(sac), m(folds$params$block_size))
```

``` r

cat(paste0("- ", names(report), ": ", unlist(report)), sep = "\n")
```

    ## - folds: block_kfold, k = 5, 7 x 2 blocks of 100,000 m (+proj=aea +lat_1=34.333391 +lat_2=36.138461 +lat_0=35.559467 +lon_0=-79.400417 +datum=WGS84 +units=m +no_defs)
    ## - score: RMSE 2.25 (per fold 1.97 to 2.61), R2 -0.063
    ## - coverage: 5 of 5 folds scored, 0 rows dropped, 0 orphan, 0 unknown id
    ## - residuals: Moran's I 0.567 (z = 12.6, p = 3.62e-36, randomisation null)
    ## - applicability: 0 of 937 prediction cells outside (DI > 0.073)
    ## - range: estimated range 279,835 m against blocks of 100,000 m

Read together, those lines say the run failed its own checks, which is
what a report is for.

The blocks are 100,000 m across and the correlation reaches about
279,835 m, so training and test data sit inside one correlation length
of each other. The score is still optimistic, and the fix is blocks at
least as wide as the range, or `auto_range = TRUE` to size them from it.

The cross-validated R2 of -0.063 is below zero, which says the trend
surface predicts held-out blocks worse than their own mean does, and
Moran’s I of 0.57 on the residuals (z = 12.6) says why: most of the
structure in `z` is the field, and a plane in the coordinates cannot
represent it. That is the point at which you fit something spatial
rather than report this.

The applicability line reads 0 of 937 here because
[`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md)
joined the covariate from the nearest observation, so no grid cell holds
a combination the model never saw. On a real covariate raster it will
not be zero;
[`vignette("diagnostics")`](https://elkronos.github.io/gis_modeling_toolkit/articles/diagnostics.md)
works through what to do when it is not.

## What each line is for

| Report | From | Because |
|----|----|----|
| The fold scheme, the block size and the estimated range, in named CRS units | `folds$params` | A cross-validation score means nothing without the geometry that produced it: blocks below the range leak, and the score comes back flattering (Roberts et al. 2017) |
| The metric **with its per-fold spread**, not the pooled value alone | `cv$fold_metrics`, [`plot_cv_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_cv_metrics.md) | 0.9 in four folds and 14 in one is a model that fails in one region, and the pooled number cannot tell that from uniform mediocrity |
| Residual Moran’s I with its null | [`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md) | It is the check on whether the spatial structure the model was built for actually got modelled |
| The share of the prediction area outside the area of applicability | [`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md) | A fitted model returns a number everywhere, including ground nothing like its training data, and the CV score does not apply there (Meyer and Pebesma 2021) |
| Interval coverage, for a model with intervals | [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md), [`plot_calibration()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_calibration.md) | Good point predictions routinely come with badly calibrated intervals, so RMSE alone can prefer the model whose uncertainty is wrong (Heaton et al. 2019) |
| Rows dropped, folds that failed, orphan and unknown IDs | `cv$n_dropped`, `cv$orphan_rows`, `cv$fold_status` | A score computed over 80 percent of the data is not a score over the data, and nothing else in the output says so |

Two of the six are about what the number does **not** cover, which is
the part that gets left out. The area of applicability bounds the map in
space and the dropped-row counts bound it in sample: without them a
reader has a figure and no idea what it was computed over.

For the aggregates themselves, report the standard error that comes with
each cell mean rather than the mean alone.
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
computes it with the observation count beside it. For per-region values
like the map above, that is the default `deff = 1` standard error used
here (or the block-kriging variance, `kr_var`, from
[`kriging_adequacy()`](https://elkronos.github.io/gis_modeling_toolkit/reference/kriging_adequacy.md)).
A `deff` setting widens it for within-cell autocorrelation, which is
right only when the regional means feed a claim about the population
(here statewide) mean;
[`?summarize_by_cell`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
says which is which. A cell mean from two observations and one from
thirteen are different claims, and printed side by side they look
identical.

``` r

populated <- cells[!is.na(cells$n), ]
summary(populated$n)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   2.000   5.750   9.500   8.333  11.000  13.000

``` r

range(populated$..se_resp_rate, na.rm = TRUE)
```

    ## [1] 0.1616037 0.6601569

## Reporting a cell value

[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
returns by default, beside each cell’s mean, a pooled value: the mean
pulled toward the mean of the 8 cells around it. This section builds
cells sized for that value on a simulated layer and reads off what a
report of them has to carry. The layer is 600 points on a square of 20
km, with a smooth field plus noise:

``` r

set.seed(7)
sx <- runif(600, 0, 2e4)
sy <- runif(600, 0, 2e4)
obs <- st_as_sf(
  data.frame(x = 4e5 + sx, y = 4e6 + sy,
             yield = 50 + 6 * sin(sx / 3000) + 6 * cos(sy / 2500) +
               rnorm(600, sd = 3)),
  coords = c("x", "y"), crs = 32617)
```

The count comes from a size, 4.5 observations per cell, scored for the
pooled value, in density-weighted cells. Those are the defaults, so the
calls name none of them. The cells are the ones that were scored, and
the summary reports the same cell value:

``` r

prof  <- resolution_profile(obs, "yield")
sel   <- select_resolution(prof)
seeds <- get_voronoi_seeds(method = "kmeans", n = sel, sample_points = obs)
tess  <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)
asg   <- assign_features_to_polygons(obs, tess$cells,
                                     polygon_id_col = "cell_id")
vals  <- summarize_by_cell(asg, "yield", id_col = "cell_id",
                           cells_sf = tess$cells, obs_level = 0.90)

st_drop_geometry(vals)[1:4, c("cell_id", "n", "resp_mean_yield",
                              "..pooled_resp_yield", "..pooled_w_resp_yield",
                              "..obs_lo_resp_yield", "..obs_hi_resp_yield")]
```

    ##   cell_id  n resp_mean_yield ..pooled_resp_yield ..pooled_w_resp_yield
    ## 1       1  3        54.96304            56.37643             0.4859504
    ## 2       2  3        53.30889            55.97057             0.4859504
    ## 3       3 11        47.98111            48.73113             0.7760979
    ## 4       4  9        52.11726            52.57068             0.7393126
    ##   ..obs_lo_resp_yield ..obs_hi_resp_yield
    ## 1            50.94400            61.67927
    ## 2            50.53815            61.27341
    ## 3            43.29871            54.03397
    ## 4            47.13826            57.87352

Each earlier default is one argument away: `estimate = "mean"` in the
summary for the plain means alone,
`estimate = NULL, obs_per_cell = NULL` in the profile for the ladder of
counts, `criterion = "cp"` for Mallows’ Cp on that ladder, and
`partition = "kmeans"` in the profile and the seeding functions for
k-means cells.

Two attributes of the frame hold the wording. `attr(, "pooled")` says
what was fitted and what each interval covers, and
`attr(, "count_rule")` says how the count was set and whether these are
the cells and rows it was scored on:

``` r

attr(vals, "pooled")
```

    ## Cell values: yield, 600 rows in 133 cells
    ##   value       : ..pooled_resp_yield, pooled toward the 8 nearest cells (centroids of the rows)
    ##   variances   : inside cells 9.38e+00 (from 467 values beyond one per cell);
    ##                 between a cell and its target 2.96e+00
    ##   own weight  : median 0.56 (..pooled_w_resp_yield)
    ##   alone       : 0.5% of the rows are the only one in their cell
    ##   loo         : 1.07e+01 (leave-one-out mean squared error of the pooled value);
    ##                 one mean for all rows: 4.21e+01
    ##   interval    : 90% for ONE OBSERVATION (..obs_lo_resp_yield, ..obs_hi_resp_yield):
    ##                 value - 5.43 to value + 5.3, quantiles of 600 leave-one-out residuals
    ## Away from the observations no interval here is calibrated.

``` r

attr(vals, "count_rule")
```

    ## Count: 133 cells, set by size (4.5 rows per cell asked; 600 rows / 133 cells = 4.51)
    ##   chosen for  : the pooled value
    ##   cells       : density, 3 density-weighted k-means++ starts (power 1, 8
    ##                 neighbours; seed 123), fitted to every row;
    ##                 these are the cells the profile scored
    ##   scored      : loo 1.07e+01 on 600 rows in 133 cells (one cell: 4.21e+01); 0.5% of rows alone
    ##   checked     : the cells and rows summarised are the ones scored

A report of cell values reads four things off them.

- **The value, and which one it is.** `..pooled_resp_yield` is the
  pooled value and `resp_mean_yield` the plain mean. Both estimate the
  cell’s level: the expected value of one further observation located in
  the cell the way its observations were. Neither is the average of a
  surface over the cell’s area. Give `n` beside either.
- **The own weight.** `..pooled_w_resp_yield` is the share of the pooled
  value that is the cell’s own mean; the rest is the mean of its
  neighbours. A cell of two observations has a smaller weight than a
  cell of ten, so its value says more about its surroundings and less
  about itself. It is not a regression weight; for a regression on cell
  values use the plain mean with `cell_weight`.
- **The interval, and what it covers.** `..obs_lo_*` and `..obs_hi_*`
  are a range for one further observation in the cell, located like the
  data. The statement holds on average over observations and not cell by
  cell, and the two offsets are the same in every cell. It is not an
  interval for the cell’s level: none is returned. `..ci_*`, from
  `conf_level`, is an interval for a mean and never a range for an
  observation. The table in `attr(vals, "intervals")` gives the sentence
  for each family in the frame.
- **Where the value applies.** A cell’s value and interval describe the
  part of the cell where its observations lie, and a Voronoi cell can
  reach far past them. A cell with no observation has `NA`. Away from
  the observations no interval here is calibrated.

The lines of a report can be written from the two records, so that they
describe the run that happened:

``` r

rec <- attr(vals, "pooled")$columns$resp_yield
iv  <- attr(vals, "intervals")
cr  <- attr(vals, "count_rule")

cat(
  sprintf("- cells: %d, set by a %s of %s observations per cell; checked: %s",
          cr$cells, cr$rule, format(cr$obs_per_cell_asked), cr$check$cells),
  sprintf("- value: the cell mean pooled toward its %d nearest cells; median own weight %.2f",
          rec$neighbours, rec$own_weight),
  sprintf("- %s: %g%% for %s", iv$columns[1], 100 * iv$level[1], iv$covers[1]),
  sprintf("- leave-one-out error %.3g against %.3g for one mean; structured: %s",
          rec$loo, rec$loo_one_cell, rec$structured),
  sprintf("- %.1f%% of observations are the only one in their cell",
          100 * rec$share_alone),
  sep = "\n")
```

    ## - cells: 133, set by a size of 4.5 observations per cell; checked: match
    ## - value: the cell mean pooled toward its 8 nearest cells; median own weight 0.56
    ## - ..obs_lo_, ..obs_hi_: 90% for one further observation in the cell, located like the data; on average over observations, not cell by cell
    ## - leave-one-out error 10.7 against 42.1 for one mean; structured: TRUE
    ## - 0.5% of observations are the only one in their cell

[`cell_values_at()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cell_values_at.md)
gives the value at other locations, with each location’s distance from
the nearest observation. `dist_ratio` is that distance over the layer’s
typical spacing, and it is the number to report beside a value read at a
place that is not an observation. The third location below lies outside
the layer:

``` r

seeds$cell_id <- tess$index
elsewhere <- st_as_sf(
  data.frame(x = 4e5 + c(5000, 12000, 26000), y = 4e6 + c(5000, 9000, 10000)),
  coords = c("x", "y"), crs = 32617)
cell_values_at(elsewhere, vals, asg, "yield", id_col = "cell_id",
               seeds = seeds)
```

    ##   cell_id    value estimate   obs_lo   obs_hi filled filled_from dist_nearest
    ## 1      36 53.88710   pooled 48.45468 59.18994  FALSE          NA     257.1487
    ## 2      82 38.95122   pooled 33.51880 44.25406  FALSE          NA     650.4402
    ## 3     128 47.42130   pooled 41.98887 52.72414  FALSE          NA    6285.9837
    ##   dist_ratio
    ## 1  0.6501451
    ## 2  1.6444979
    ## 3 15.8927549

It gets a value and an interval like the others, because every location
has a nearest seed. Nothing flags it: the `dist_ratio` beside it does.
On simulated layers the interval at 0.90 held 91.5 percent of held-out
observations with a `dist_ratio` below 1 and 69.5 percent of those above
8 (`dev/sim-pooled-values.R` in the package’s repository;
[`?cell_values_at`](https://elkronos.github.io/gis_modeling_toolkit/reference/cell_values_at.md)
gives the settings and the classes between).

## What to save

Enough to rerun it, which is less than it sounds:

- the seed (`folds$params$seed`) and the fold assignment
  (`folds$assignment`),
- the regions, with
  [`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md)
  applied, as a GeoPackage,
- `folds$params` and the fit’s `$info`, which hold every decision made
  for you,
- the session’s package versions.

The regions are the one output that cannot be regenerated from the
others without care.
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
and the k-means sweep in
[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
are seeded, so they repeat, but only at the same seed and the same
package version. Save the layer.

## Next

[`vignette("getting-started")`](https://elkronos.github.io/gis_modeling_toolkit/articles/getting-started.md)
is the shortest path through the pipeline that produced all of this.
[`vignette("diagnostics")`](https://elkronos.github.io/gis_modeling_toolkit/articles/diagnostics.md)
covers the residual and applicability checks the table asks you to
report, and
[`vignette("spatial-cross-validation")`](https://elkronos.github.io/gis_modeling_toolkit/articles/spatial-cross-validation.md)
covers fold geometry and how to read a CV result down to the last row.

## References

Heaton, M. J., Datta, A., Finley, A. O., et al. (2019). A case study
competition among methods for analyzing large spatial data. *Journal of
Agricultural, Biological and Environmental Statistics* 24, 398-425.
<https://doi.org/10.1007/s13253-018-00348-w>

Meyer, H. and Pebesma, E. (2021). Predicting into unknown space?
Estimating the area of applicability of spatial prediction models.
*Methods in Ecology and Evolution* 12(9), 1620-1633.
<https://doi.org/10.1111/2041-210X.13650>

Roberts, D. R., Bahn, V., Ciuti, S., et al. (2017). Cross-validation
strategies for data with temporal, spatial, hierarchical, or
phylogenetic structure. *Ecography* 40(8), 913-929.
<https://doi.org/10.1111/ecog.02881>

    ## R version 4.6.1 (2026-06-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] ggplot2_4.0.3         spatialkit_2.0.0.9000 sf_1.1-3             
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] s2_1.1.13          sass_0.4.10        generics_0.1.4     class_7.3-23      
    ##  [5] KernSmooth_2.23-26 lattice_0.22-9     digest_0.6.39      magrittr_2.0.5    
    ##  [9] evaluate_1.0.5     grid_4.6.1         RColorBrewer_1.1-3 fastmap_1.2.0     
    ## [13] Matrix_1.7-5       jsonlite_2.0.0     e1071_1.7-17       DBI_1.3.0         
    ## [17] viridisLite_0.4.3  scales_1.4.0       textshaping_1.0.5  jquerylib_0.1.4   
    ## [21] cli_3.6.6          rlang_1.3.0        units_1.0-1        intervals_0.15.5  
    ## [25] withr_3.0.3        cachem_1.1.0       yaml_2.3.12        otel_0.2.0        
    ## [29] FNN_1.1.4.1        tools_4.6.1        dplyr_1.2.1        spacetime_1.3-4   
    ## [33] logger_0.4.3       vctrs_0.7.3        R6_2.6.1           zoo_1.9-1         
    ## [37] proxy_0.4-29       lifecycle_1.0.5    classInt_0.4-11    fs_2.1.0          
    ## [41] htmlwidgets_1.6.4  ragg_1.5.2         pkgconfig_2.0.3    desc_1.4.3        
    ## [45] pkgdown_2.2.1      pillar_1.11.1      bslib_0.12.0       gtable_0.3.6      
    ## [49] Rcpp_1.1.2         glue_1.8.1         gstat_2.1-6        systemfonts_1.3.2 
    ## [53] xfun_0.62          tibble_3.3.1       tidyselect_1.2.1   knitr_1.52        
    ## [57] farver_2.1.2       htmltools_0.5.9    labeling_0.4.3     rmarkdown_2.32    
    ## [61] xts_0.14.3         wk_0.9.5           compiler_4.6.1     S7_0.2.2          
    ## [65] sp_2.2-4
