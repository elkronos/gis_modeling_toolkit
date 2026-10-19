# Forward model selection for geographically weighted regression

Grows a GWR model one predictor at a time and scores every intermediate
model with a corrected Akaike information criterion, the sweep
[`GWmodel::gwr.model.selection()`](https://rdrr.io/pkg/GWmodel/man/gwr.model.selection.html)
defines, and returns a ranked table with each model's distance from the
best. The model without predictors is scored too. Every model is fitted
by the package's own local solve (see
[`fit_gwr_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md),
**The local solve**), so a candidate that is constant inside some
windows, a regional indicator say, is scored on what those windows do
estimate and neither stops the sweep nor scores as rounding error.

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
```

## Arguments

- data_sf:

  An `sf` object with response, predictors and geometry.

- response_var:

  Response column name.

- candidate_vars:

  Character vector naming at least two **numeric** predictors to choose
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
  otherwise a distance in the units of the **projected** CRS the sweep
  runs in, which
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

  Weighting kernel. One of `"bisquare"` (default), `"gaussian"`,
  `"tricube"`, `"boxcar"`, `"exponential"`.

- bw_approach:

  Criterion for the bandwidth search: `"AICc"` (default) or `"CV"`, the
  leave-one-out sum of squared prediction errors. Matches
  [`fit_gwr_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)'s
  default. When observations share a location the search leaves out
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

  Print no progress. Default `TRUE`. With `FALSE` one line is written
  (with [`cat()`](https://rdrr.io/r/base/cat.html)) as each model is
  fitted, to watch a long sweep progress.

- .engine:

  Internal; injectable backend used for testing.

## Value

An object of class `gwr_model_selection`, a list with: `best` (character
vector of the selected predictors, the model with the lowest criterion;
`character(0)` when that is the model without predictors); `table`
(ranked data.frame of every model evaluated, with columns `rank`,
`n_vars`, `variables`, `criterion` and `delta`, the criterion minus the
best model's; the model without predictors is the row `<none>`;
`criterion` is `NA`, and the model ranked last, where it could not be
evaluated or where AICc is undefined because the model's effective
number of parameters \\tr(S)\\ is not below \\n - 2\\, which raises a
warning); `within_2` (the rows of `table` whose `delta` is at most 2)
and `smallest_within_2` (the predictors of the smallest model among
them); `criterion` (label for the criterion actually read, noting when
it had to be located positionally); `criterion_by_name` (logical:
whether that column was found by its name rather than by the documented
position), `criterion_column` (the column it was read from) and
`criterion_verified` (logical: `FALSE` exactly when the column was read
positionally from a table that did not have the four documented columns,
which is the case the log calls unverified. Gate a script on this field
instead of on the label); `response_var` and `candidate_vars` (the
response and the full candidate set the sweep ran over, both echoed by
[`print()`](https://rdrr.io/r/base/print.html)); `bandwidth`,
`bandwidth_source`, `adaptive` and `kernel` (the smoothing held fixed
across the sweep, and where it came from: `"supplied"`, or the
criterion's name followed by `"search on the full model"`, the criterion
being `"leave-one-location-out CV"` when observations share a location);
`bandwidth_search` (what the search scored, as
[`fit_gwr_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
keeps it; `NULL` when the bandwidth was supplied); `n_obs`, `n_models`
(the rows of `table`), `used_dmat` (always `FALSE`); and `raw` (the
forward sweep in the layout
[`GWmodel::gwr.model.selection()`](https://rdrr.io/pkg/GWmodel/man/gwr.model.selection.html)
returns: the two-element list of its model list and its diagnostic
table, without the model that has no predictors).

## What this optimises, and what it does not

The criterion is **in-sample**. AICc penalises the effective number of
parameters, so it is not the same thing as maximising fit, but it is
still computed on the data the model was fitted to, and under spatial
autocorrelation an in-sample criterion is optimistic in a way that a
spatially blocked estimate is not. Treat this as fast screening.
[`select_features_forward`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)
performs the same forward search against a spatially blocked
cross-validated score; it costs far more and is the one to trust when
the answer matters. When the two disagree, the disagreement is itself
informative. It usually means a candidate is predictive only locally.

One further limitation follows from the method itself:

- **One bandwidth for every model.** Comparing criteria across models
  requires holding the smoothing fixed, but the bandwidth is itself a
  fitted quantity, and the value chosen for the full model is not
  optimal for a one-predictor model. This is how the method is defined
  (Lu et al. 2014); it is not an implementation shortcut. Refit the
  selected model with `bandwidth = NULL` to re-optimise once the
  variable set is settled.

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

[`select_features_forward`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)
for the blocked cross-validated counterpart,
[`fit_gwr_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
to fit the selected model.

Other cross-validation:
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md),
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
[`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md),
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md),
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
[`fold_separation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fold_separation.md),
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md),
[`sac_nugget()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md),
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)

## Examples

``` r
if (requireNamespace("sp", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  n <- 80
  dat <- st_as_sf(
    data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
               a = rnorm(n), b = rnorm(n), noise = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
  dat$z <- 2 * dat$a - dat$b + rnorm(n, 0, 0.5)
  sel <- gwr_model_selection(dat, "z", c("a", "b", "noise"), bandwidth = 30)
  print(sel)              # every model tried, by AICc, with the gaps
  print(sel$best)         # the winning predictor set: a and b, not noise
  fit <- fit_gwr_model(dat, "z", sel$best)
  fit
}
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
#>     1      2         a + b  136.6185   0.00000
#>     2      3 a + b + noise  150.6764  14.05793
#>     3      1             a  273.2590 136.64056
#>     4      2     a + noise  278.3346 141.71614
#>     5      1             b  360.8372 224.21866
#>     6      0        <none>  369.4529 232.83444
#>     7      1         noise  378.2502 241.63171
#> 
#> Selected: z ~ a + b
#> 
#> No other model is within 2 of the best; the runner-up is 14.06 away.
#> 
#> The criterion is in-sample and every model shares one bandwidth.
#> Confirm with select_features_forward() before relying on this.
#> [1] "a" "b"
#> <GWR> spatial model fit
#>   Formula : z ~ a + b
#>   n       : 80
#>   CRS     : EPSG:32632
#>   Bandwidth: 80 neighbours (adaptive, bisquare kernel)
#>   AICc    : 122.82
#>   Eff. parameters: 7.2 of 80 observations
```
