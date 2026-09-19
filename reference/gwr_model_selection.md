# Forward model selection for geographically weighted regression

Wraps
[`GWmodel::gwr.model.selection()`](https://rdrr.io/pkg/GWmodel/man/gwr.model.selection.html),
which grows a GWR model one predictor at a time and scores every
intermediate model with a corrected Akaike information criterion.
GWmodel returns two loosely-coupled lists; this function returns a
ranked table.

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
  among. Factor, character and logical candidates are refused: GWmodel
  fits a factor as several model-matrix columns while this sweep counts
  it as one variable, so the criteria would not be comparable, and
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md),
  the documented next step, takes only numerics. Encode them as numeric
  indicators first.

- bandwidth:

  Bandwidth held fixed across all candidate models. If `NULL` (default)
  it is selected with
  [`GWmodel::bw.gwr()`](https://rdrr.io/pkg/GWmodel/man/bw.gwr.html) on
  the model containing every candidate. Integer neighbour count when
  `adaptive = TRUE`; otherwise a distance in the units of the
  **projected** CRS the sweep runs in, which
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
  may have chosen for you. Geographic input is projected before the
  bandwidth is used, so a value in degrees would be read as metres.

- adaptive:

  Logical; adaptive (nearest-neighbour) bandwidth. Default `TRUE`.

- kernel:

  Weighting kernel. One of `"bisquare"` (default), `"gaussian"`,
  `"tricube"`, `"boxcar"`, `"exponential"`.

- bw_approach:

  Criterion for the bandwidth search: `"AICc"` (default) or `"CV"`.
  Matches
  [`fit_gwr_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)'s
  default.

- max_models:

  Refuse the call if the sweep would exceed this many model fits.
  Default 200, which admits up to 19 candidates.

- dmat_max_n:

  Precompute and reuse an `n` x `n` distance matrix when `n` is at most
  this. Default 2000 (about 32 MB). Set to 0 to disable.

- quiet:

  Discard GWmodel's progress output. Default `TRUE`, because GWmodel
  writes it with bare [`cat()`](https://rdrr.io/r/base/cat.html) that no
  [`suppressMessages()`](https://rdrr.io/r/base/message.html) can
  silence, and emits one block per candidate model, so it scales with
  the square of the candidate count. Set `FALSE` to watch a long sweep
  progress.

- .engine:

  Internal; injectable backend used for testing.

## Value

An object of class `gwr_model_selection`, a list with: `best` (character
vector of the selected predictors); `table` (ranked data.frame of every
model evaluated, with columns `rank`, `n_vars`, `variables` and
`criterion`); `criterion` (label for the criterion actually read, noting
when it had to be located positionally); `criterion_by_name` (logical:
whether that column was found by its name rather than by the documented
position), `criterion_column` (the column it was read from) and
`criterion_verified` (logical: `FALSE` exactly when the column was read
positionally from a table that did not have the four documented columns,
which is the case the log calls unverified. Gate a script on this field
instead of on the label); `response_var` and `candidate_vars` (the
response and the full candidate set the sweep ran over, both echoed by
[`print()`](https://rdrr.io/r/base/print.html)); `bandwidth`,
`bandwidth_source`, `adaptive` and `kernel` (the smoothing held fixed
across the sweep, and where it came from); `n_obs`, `n_models`,
`used_dmat`; and `raw` (GWmodel's unmodified return: the two-element
list of its model list and its diagnostic table).

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

Two further limitations follow from the method itself:

- **One bandwidth for every model.** Comparing criteria across models
  requires holding the smoothing fixed, but the bandwidth is itself a
  fitted quantity, and the value chosen for the full model is not
  optimal for a one-predictor model. This is how the method is defined
  (Lu et al. 2014); it is not an implementation shortcut. Refit the
  selected model with `bandwidth = NULL` to re-optimise once the
  variable set is settled.

- **The null model is never evaluated.** The sweep starts from one
  predictor, so the result always names at least one. It cannot tell you
  that none of the candidates help.

## Cost

The sweep fits `p * (p + 1) / 2` GWR models for `p` candidates (55 at p
= 10, 210 at p = 20), each over all `n` locations. `max_models` stops
the call before it runs for hours.

## Using \$raw with GWmodel directly

`raw` is GWmodel's own `list(model.list, GWR.df)`, so its two elements
have to be unpacked before GWmodel's own helpers will take them:
[`GWmodel::gwr.model.view()`](https://rdrr.io/pkg/GWmodel/man/gwr.model.view.html)
takes `(DeVar, InDeVars, model.list)`, so the call is


      GWmodel::gwr.model.view(sel$response_var, sel$candidate_vars, sel$raw[[1]])

Pass `sel$raw[[1]]`, not `sel$raw`. The diagnostic table is
`sel$raw[[2]]`, an unlabelled numeric matrix whose columns are
`bandwidth`, `AIC`, `AICc`, `RSS` in that order; the `criterion` column
of `$table` is its third column.

## References

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
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md),
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)

## Examples

``` r
if (requireNamespace("GWmodel", quietly = TRUE) &&
    requireNamespace("sp", quietly = TRUE)) {
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
  sel$best
  fit <- fit_gwr_model(dat, "z", sel$best)
}
```
