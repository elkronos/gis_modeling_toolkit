# Select GWR predictors by forward selection on AICc

Grows a GWR model one predictor at a time and scores every intermediate
model with a corrected Akaike information criterion, the sweep
[`GWmodel::gwr.model.selection()`](https://rdrr.io/pkg/GWmodel/man/gwr.model.selection.html)
defines, and returns a ranked table with each model's distance from the
best. The model without predictors is scored too. Every model is fitted
by the package's own local solve (see "The local solve" in
[`fit_gwr_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)),
so a candidate that is constant inside some windows, a regional
indicator say, is scored on what those windows do estimate and neither
stops the sweep nor scores as rounding error.

## Usage

``` r
gwr_model_selection(
  data_sf,
  response_var,
  candidate_vars,
  bandwidth = NULL,
  adaptive = TRUE,
  kernel = c("bisquare", "gaussian", "tricube", "boxcar", "exponential"),
  bw_approach = c("AICc", "CV"),
  max_models = 200L,
  dmat_max_n = 2000L,
  quiet = TRUE,
  .engine = .gwr_ms_engine
)

# S3 method for class 'gwr_model_selection'
print(x, n = 10L, ...)
```

## Arguments

- data_sf:

  An `sf` object with response, predictors and geometry.

- response_var:

  Response column name.

- candidate_vars:

  Character vector naming at least two numeric predictors to choose
  among. Factor, character and logical candidates are refused: a factor
  is several model-matrix columns while this sweep counts it as one
  variable, so the criteria would not be comparable, and
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md),
  the documented next step, takes only numerics. Encode them as numeric
  indicators first.

- bandwidth:

  Bandwidth held fixed across all candidate models. If `NULL` (default)
  it is searched for on the model containing every candidate, by the
  grid-then-refine search of
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  (see its **Bandwidth selection**); what the search scored is returned
  in `bandwidth_search`. Integer neighbour count when `adaptive = TRUE`;
  otherwise a distance in the units of the projected CRS the sweep runs
  in, which
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
  may have chosen for you. Geographic input is projected before the
  bandwidth is used, so a value in degrees would be read as metres; a
  fixed bandwidth below a ten-thousandth of the data's extent raises a
  warning saying so, as in
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md).
  An adaptive count too small for the full model is raised, with a
  warning, to the number of candidates plus 3 for the bisquare and
  tricube kernels (which give the farthest neighbour in a window weight
  0), plus 2 for the others. That floor is enough unless several
  neighbours tie at the kernel's edge (a regular grid), which leaves a
  window fewer weighted points; then use a larger bandwidth. An adaptive
  count below 1 or above R's largest integer is refused, as in
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md).
  One above the number of observations is capped at it, with a warning.

- adaptive:

  Logical; adaptive (nearest-neighbour) bandwidth. Default `TRUE`.

- kernel:

  Weighting kernel: one of `"bisquare"` (default), `"gaussian"`,
  `"tricube"`, `"boxcar"` or `"exponential"`, as in
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md).

- bw_approach:

  Criterion for the bandwidth search on the full model: `"AICc"`
  (default, the criterion
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  uses) or `"CV"`, the leave-one-out sum of squared prediction errors.
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  has no CV option, so with `"CV"` pass the selected bandwidth to it
  explicitly. When observations share a location the search leaves out
  whole locations whichever is asked for; see **Observations that share
  a location**.

- max_models:

  Refuse the call if the sweep would exceed this many model fits.
  Default 200, which admits up to 19 candidates.

- dmat_max_n:

  Not used. The sweep takes its distances window by window and builds no
  `n` x `n` distance matrix, so a fixed bandwidth works at any `n`. The
  argument is accepted so that calls written with it run.

- quiet:

  Print no progress. Default `TRUE`, unlike
  [`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md).
  With `FALSE` one line is written (with
  [`cat()`](https://rdrr.io/r/base/cat.html)) as each model is fitted,
  to watch a long sweep progress.

- .engine:

  Internal; leave at the default.

- x:

  A `gwr_model_selection` object.

- n:

  Number of top-ranked models
  [`print()`](https://rdrr.io/r/base/print.html) shows. Default 10.

- ...:

  Not used.

## Value

An object of class `gwr_model_selection`, a list with:

- `best`:

  Character vector of the selected predictors, the model with the lowest
  criterion; `character(0)` when that is the model without predictors.

- `table`:

  Ranked data.frame of every model evaluated, with columns `rank`,
  `n_vars`, `variables`, `criterion` and `delta`, the criterion minus
  the best model's. The model without predictors is the row `<none>`.
  `criterion` is `NA`, and the model is ranked last, where it could not
  be evaluated or where AICc is undefined because the model's effective
  number of parameters \\tr(S)\\ is not below \\n - 2\\ (with a
  warning).

- `within_2`, `smallest_within_2`:

  The rows of `table` whose `delta` is at most 2, and the predictors of
  the smallest model among them.

- `criterion`:

  `"AICc"`.

- `criterion_by_name`, `criterion_column`, `criterion_verified`:

  Where the criterion was read from the diagnostic table: `TRUE`, 3 and
  `TRUE` with the built-in engine.

- `response_var`, `candidate_vars`:

  The response and the full candidate set the sweep ran over, both
  echoed by [`print()`](https://rdrr.io/r/base/print.html).

- `bandwidth`, `bandwidth_source`, `adaptive`, `kernel`:

  The smoothing held fixed across the sweep, and where it came from:
  `"supplied"`, or the criterion's name followed by
  `"search on the full model"`, the criterion being
  `"leave-one-location-out CV"` when observations share a location.

- `bandwidth_search`:

  What the search scored, as
  [`fit_gwr_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  keeps it; `NULL` when the bandwidth was supplied.

- `n_obs`, `n_models`:

  The number of observations, and the number of rows of `table`.

- `used_dmat`:

  Always `FALSE`.

- `raw`:

  The forward sweep in the layout
  [`GWmodel::gwr.model.selection()`](https://rdrr.io/pkg/GWmodel/man/gwr.model.selection.html)
  returns: the two-element list of its model list and its diagnostic
  table, without the model that has no predictors.

## Scope of the criterion

AICc penalises the effective number of parameters but is computed on the
data the model was fitted to, and under spatial autocorrelation an
in-sample criterion is optimistic in a way that a spatially blocked
estimate is not. Use this function as a fast screen and confirm with
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md),
which runs the same forward search on a spatially blocked
cross-validated score at a higher cost.

Every model is fitted at one bandwidth, as Lu et al. (2014) define the
method: comparing criteria across models needs the smoothing held fixed.
The bandwidth chosen for the full model is not optimal for smaller ones,
so refit the selected model with `bandwidth = NULL` once the set is
settled.

## The model without predictors, and the gaps

The forward sweep starts from one predictor, so by itself it always
names at least one, whether or not any candidate helps. The model
without predictors, a local mean at the sweep's bandwidth, is therefore
scored as well and ranked with the others as `<none>`. When it ranks
first, `best` is `character(0)`: no candidate improves on it.

The `delta` column of `table` is each model's criterion minus the best
one's. Differences of AICc below about 2 do not separate two models
(Burnham and Anderson 2002), so `within_2` lists the models that close
to the best and `smallest_within_2` names the one among them with the
fewest predictors, which is the one to prefer when several are that
close. `best` is still the model with the lowest criterion.

## Observations that share a location

When rows of `data_sf` have the same coordinates, the bandwidth search
leaves out whole locations, for either `bw_approach` (see
[`fit_gwr_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md),
**Observations that share a location**); `bandwidth_source` then reads
`"leave-one-location-out CV search on the full model"`. The models are
still ranked by AICc, which counts the rows at one location as
independent, so the gaps overstate what a predictor adds when a location
has an effect of its own; a warning of class
`"spatialkit_gwr_repeated_locations"` says so with the counts. A model
whose fit reproduces the mean of the observations at each location
(every such model scores the same, whatever its predictors) has
criterion `NA` and is listed last, as is one with a term that no local
window can estimate; a warning names them.

## Cost

The sweep fits `p * (p + 1) / 2` GWR models for `p` candidates (55 at p
= 10, 210 at p = 20), and the one without predictors, each over all `n`
locations. `max_models` stops the call before it runs for hours.

## Using \$raw with GWmodel directly

`raw` has the layout of GWmodel's own `list(model.list, GWR.df)`, so its
two elements have to be unpacked before GWmodel's helpers will take
them:
[`GWmodel::gwr.model.view()`](https://rdrr.io/pkg/GWmodel/man/gwr.model.view.html)
takes `(DeVar, InDeVars, model.list)`, so the call is


      GWmodel::gwr.model.view(sel$response_var, sel$candidate_vars, sel$raw[[1]])

Pass `sel$raw[[1]]`, not `sel$raw`. The diagnostic table is
`sel$raw[[2]]`, a numeric matrix whose columns are `bandwidth`, `AIC`,
`AICc`, `RSS` in that order (and are named so); the `criterion` column
of `$table` is its third column.

## Printing

[`print()`](https://rdrr.io/r/base/print.html) shows the response, the
candidate predictors, and the top-ranked models with their criterion
values and, in the `delta` column, how far each is from the best. Below
the table it lists the models within 2 of the best, which a difference
in AICc does not separate (Burnham and Anderson 2002), and names the
smallest of them. When several models are that close, treat the choice
of predictors as unsettled; when the model without predictors (`<none>`)
is among them, no candidate has been shown to help.
[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.

## References

Burnham, K. P. and Anderson, D. R. (2002). *Model Selection and
Multimodel Inference: A Practical Information-Theoretic Approach*, 2nd
edition. Springer.

Lu, B., Harris, P., Charlton, M. and Brunsdon, C. (2014). The GWmodel R
package: further topics for exploring spatial heterogeneity using
geographically weighted models. *Geo-spatial Information Science* 17(2),
85–101.
[doi:10.1080/10095020.2014.917453](https://doi.org/10.1080/10095020.2014.917453)

## See also

[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)
for the blocked cross-validated counterpart,
[`plot.gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.gwr_model_selection.md)
to draw the models, and
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
to fit the selected model.

Other predictor selection:
[`plot.feature_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.feature_selection.md),
[`plot.gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.gwr_model_selection.md),
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)

## Examples

``` r
library(sf)
set.seed(1)
n <- 80
dat <- st_as_sf(
  data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
             a = rnorm(n), b = rnorm(n), noise = rnorm(n)),
  coords = c("x", "y"), crs = 32632
)
# The slope on a rises from 1 in the west to 3 in the east.
east <- (st_coordinates(dat)[, 1] - 5e5) / 1000
dat$z <- (1 + 2 * east) * dat$a - dat$b + rnorm(n, 0, 0.3)
sel <- gwr_model_selection(dat, "z", c("a", "b", "noise"), bandwidth = 30)
sel                     # every model tried, by AICc, with the gaps
#> Geographically weighted regression - forward model selection
#> 
#>   response    : z
#>   candidates  : 3 (a, b, noise)
#>   observations: 80
#>   models      : 7
#>   bandwidth   : 30 (adaptive; supplied)
#>   kernel      : bisquare
#>   criterion   : AICc
#> 
#>  rank n_vars     variables criterion     delta
#>     1      2         a + b  70.55026   0.00000
#>     2      3 a + b + noise  81.27942  10.72916
#>     3      1             a 262.02384 191.47358
#>     4      2     a + noise 268.26330 197.71304
#>     5      1             b 371.33475 300.78449
#>     6      0        <none> 374.24430 303.69404
#>     7      1         noise 384.67142 314.12116
#> 
#> Selected: z ~ a + b
#> 
#> No other model is within 2 of the best; the runner-up is 10.73 away.
#> 
#> The criterion is in-sample and every model shares one bandwidth.
#> Confirm with select_features_forward() before relying on this.
sel$best                # the winning predictor set
#> [1] "a" "b"
# Refit with the bandwidth chosen afresh for the selected set.
fit_gwr_model(dat, "z", sel$best)
#> <GWR> spatial model fit
#>   Formula : z ~ a + b
#>   n       : 80
#>   CRS     : EPSG:32632
#>   Bandwidth: 35 neighbours (adaptive, bisquare kernel)
#>   AICc    : 70.30
#>   Eff. parameters: 22.0 of 80 observations
```
