# Choose a cell count from the elbow of the k-means WSS curve

Fits k-means to the coordinates at \\k = 1\\ to `max_levels` cells and
returns the count at the elbow of the within-cluster sum of squares
(WSS), read on log-log axes, with its two neighbours. Points with no
cluster structure have no elbow, and the result is then `NA` with a
warning. For a count chosen from the response, use
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md).

## Usage

``` r
determine_optimal_levels(
  data_sf,
  max_levels = 12L,
  top_n = 3L,
  sample_n = 1500L,
  set_seed = 123L,
  response_var = NULL,
  predictor_vars = NULL,
  criterion = c("geometric", "morans_i", "combined"),
  select_on = c("all", "split"),
  on_no_elbow = c("na", "warn", "error")
)
```

## Arguments

- data_sf:

  An sf object. Non-point geometries are reduced to points with
  [`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md).
  Features without a location are dropped with a warning, and a layer
  that mixes XY and XYZ geometries or carries M is read in x and y (see
  [spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)).

- max_levels:

  Integer upper bound on levels. Default 12. The sweep also stops at the
  number of distinct locations (k-means cannot place more centres) and
  one short of the number of points.

- top_n:

  Integer; how many candidates to return. Default 3. Under
  `criterion = "geometric"` the candidate set is the elbow and its two
  immediate neighbours, so at most 3 values are returned however large
  `top_n` is; only the model-aware criteria can return more. Under those
  it also sets the counts scored, the elbow plus or minus
  `max(4, top_n)`: a `top_n` above 4 scores more counts, and under
  `"morans_i"` one of them can rank first. A count scores the same
  whichever `top_n` puts it in the window.

- sample_n:

  Integer; subsample size for speed. Default 1500.

- set_seed:

  Integer RNG seed for the subsample and the k-means++ restarts;
  restored afterwards. Default 123. The rows are put in coordinate order
  before either, so the answer does not depend on the order they come
  in.

- response_var:

  Name of the response column (numeric or logical), or `NULL`. Read
  under `criterion = "morans_i"`, which ranks the counts around the
  elbow by Moran's I of the cell means, and `"combined"`, which reports
  that statistic beside the elbow's ranking; under `"geometric"`, the
  default, it is not read. Logicals are read as 0/1; a factor or
  character response raises an error and is never coerced, because the
  residuals of an OLS fit to arbitrary level codes carry no meaning to
  test for autocorrelation. Rows where it, or a predictor, is missing or
  non-finite stay in the WSS sweep and the cells and are left out of
  Moran's I; a logged warning gives their number.

- predictor_vars:

  Optional predictor column names, read where `response_var` is. Must be
  numeric or logical (logicals are read as 0/1); factor/character
  columns raise an error. Without them the model-aware criteria regress
  the cell means on an intercept alone, as
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  does.

- criterion:

  One of `"geometric"` (default: the elbow and its two neighbours, read
  off the coordinates alone), `"morans_i"` (among the candidates around
  the elbow, the \\k\\ at which Moran's I of the cell means, or of their
  regression residuals, is least significant: the count at which
  neighbouring cells stop resembling each other, which need not be the
  count that represents the response best), or `"combined"` (the
  candidates around the elbow in the order of the WSS curve's log-log
  sag, the quantity the elbow is read from, with Moran's I of the cell
  means and its \\z\\ in the diagnostics; the statistic does not take
  part in the ranking). Under `"morans_i"` a warning of class
  `spatialkit_diagnostic_ranking` says that the ranking is by a
  diagnostic of the cell means; its condition carries `k` (the count
  ranked first), `criterion`, `first_computable` (the first count at
  which the statistic is finite) and `window` (the first and last count
  tried). `"morans_i"` and `"combined"` fall back to `"geometric"`, with
  a warning, when `response_var` is not given, and with a logged warning
  when no candidate clears the nine-cell floor described under
  "Model-aware criteria". A `response_var` or `predictor_vars` naming a
  column that is not in `data_sf` is an error under every criterion.
  Supplying both `response_var` and `predictor_vars` does not change the
  criterion: under `"geometric"` neither is read, and a line logged at
  the INFO level says that `"combined"` reports Moran's I.

- select_on:

  `"all"` (default) selects on every point; `"split"` reads the response
  on one spatially blocked half of the points only and returns the other
  half as the set to estimate on, so that the standard errors computed
  downstream on the chosen cells are not post-selection. The count is
  still chosen for the whole layer. See "Post-selection inference".

- on_no_elbow:

  What to do when the count is not a finding about the data: the WSS
  curve has no elbow, the ladder is too short to show one, or there are
  too few points or locations for a ladder at all. Each choice raises a
  condition of class `spatialkit_no_elbow` there, which carries `k` (the
  count the chord rule or the ladder gives), `sag` (the largest log-log
  sag of the curve; an elbow needs \\\log 1.25\\) and
  `structured = FALSE`, so a caller can catch it without matching the
  message. `"na"` (default) returns `NA_integer_` there, with a warning
  of that class, and every result under `"na"` carries the attributes
  `structured` (`TRUE` when the WSS curve has an elbow, `FALSE` with the
  `NA`) and `sag`, so one test of `attr(k, "structured")` covers both
  outcomes. `"warn"` returns the counts the chord rule or the ladder
  gives (see Value), with the same warning; its results carry neither
  attribute. `"error"` stops with an error of that class. Under the
  model-aware criteria the candidates are scored around the elbow: with
  no elbow, `"na"` and `"error"` score none, and under `"warn"` the
  ladder places that window.

## Value

An integer vector of candidate counts, best first:

- `"geometric"`:

  The elbow, then its lower and upper neighbours.

- `"morans_i"`:

  The candidates in the order of \\\|z\|\\, with a warning of class
  `spatialkit_diagnostic_ranking`.

- `"combined"`:

  The candidates around the elbow in the order of the sag.

- No elbow:

  `NA_integer_`, with a warning of class `spatialkit_no_elbow`; under
  `on_no_elbow = "warn"`, the ladder's choice and its neighbours, with
  the same warning. So it is when the ladder is too short to show an
  elbow (`max_levels` of 1 or 2), and for a layer with fewer than three
  points or every point on one location, which has no ladder (1 under
  `"warn"`). A layer with no points is an error.

`k[1]` is the top-ranked count, or `NA`, on every path, and `top_n = 1`
returns it alone. The vector can carry these attributes:

- `structured`, `sag`:

  Under the default `on_no_elbow = "na"`, on every result: `structured`
  is `TRUE` when the WSS curve has an elbow and `FALSE` with the `NA`,
  and `sag` is the curve's largest log-log sag. None under `"warn"` or
  `"error"`. The geometric path carries no other attribute.

- `diagnostics`:

  Under `"morans_i"` and `"combined"`: a list with per-\\k\\ Moran's I
  (`moran_i`) and standardised deviates (`moran_z`), the WSS curve
  (`wss`) with the relative between-restart spread at each \\k\\
  (`wss_spread`), the number of rising steps on it (`wss_bumps`), the
  restart budget (`nstart`), the elbow the scored neighbourhood was
  drawn around (`knee_k`), whether the WSS curve has an elbow
  (`structured`) and its largest log-log sag (`sag`), the \\k\\ at which
  k-means failed (`failed_k`; their `wss` entries are interpolated from
  the neighbours, not measured) and the \\k\\ the model-aware pass
  scored (`eval_ks`, the elbow's neighbourhood). Under `"combined"` it
  also carries the WSS of the re-run clustering at those \\k\\
  (`wss_eval`, the clustering Moran's I is computed on; reported, not
  ranked on), the rank of the sag of `wss` that ordered them
  (`combined_rank`, named by \\k\\; with no elbow every candidate ties
  and the order is by distance from `knee_k`) and
  `criterion = "combined"`. When `"combined"` returned the geometric
  ranking because the elbow is below ten cells, it carries
  `criterion = "geometric"` and `fallback` (the reason) in place of
  `combined_rank`. The `NA` the default returns on a model-aware call
  whose WSS curve has no elbow (or a ladder too short to show one)
  carries the attribute too, with what the sweep measured: `moran_i` and
  `moran_z` all `NA` and `eval_ks` empty, since no count was scored,
  `knee_k` the count the chord rule lands on, and `structured = FALSE`.
  Under `"combined"`, Moran's I computable at no candidate of a
  neighbourhood that reaches ten cells leaves the attribute in place,
  with `moran_i` and `moran_z` all `NA`, and a logged warning says so.
  The attribute is absent when the model-aware path falls back to the
  geometric result (no viable \\k\\ in the elbow's neighbourhood, a
  neighbourhood that ends at nine cells or fewer, or, under
  `"morans_i"`, Moran's I computable at no candidate); these fallbacks
  are logged as warnings. It is absent too on the `NA` of a layer with
  no ladder (fewer than three points, one location, or k-means failing
  at every \\k\\ from 2).

- `split`:

  With `select_on = "split"`, on every path: a list with `selection` and
  `estimation` (integer row positions in `data_sf`), `method` and
  `seed`.

A rising WSS curve is logged on every path. For a per-level table of
criteria, cell support, restart spread and the flat region, see
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md).

## How the elbow is read

Points with no cluster structure have a WSS curve close to \\c/k\\. The
classical rule, the point furthest below the chord from the first to the
last \\k\\, still finds a "knee" on such a curve on linear axes, at
about \\\sqrt{K\_{max}}\\: 4 at the default `max_levels = 12` and 13 at
160, whatever the data. On \\\log k\\ against \\\log\\ WSS that curve is
a straight line, while separated clusters fall faster than it until
there is one cell per cluster and like it after, a bend at the cluster
count. The elbow is therefore the \\k\\ whose \\\log\\ WSS sags furthest
below the straight line joining \\k = 1\\ and \\K\_{max}\\ on those
axes, and it counts as one only when the sag is at least \\\log 1.25\\
(the WSS a fifth below the power law through the ends).

An elongated extent also bends, at about its aspect ratio, because the
first cuts go across its long axis; past the threshold that bend is
reported as an elbow, and it describes the shape of the extent rather
than clusters in it. When locations repeat (stations visited many
times), the sweep can reach one cell per distinct location, where the
WSS is zero (to within \\10^{-12}\\ of the total). That \\k\\ has no
place on log axes and is left out of the line, but the fall to zero is
the sharpest bend there is: when the rest of the curve has no elbow, the
number of distinct locations is the elbow.

## When there is no elbow

With no elbow the function returns `NA`, so that a count which is not a
finding is not passed on as one. The result carries `structured = FALSE`
and the largest sag as `sag`, and a warning of class
`spatialkit_no_elbow` says why. The count the linear-axis rule gives
there is the ladder's and not the data's: about \\\sqrt{K\_{max}}\\, 3
or 4 cells at the default, close to one mean for the whole layer. The
warning carries it as `k` and the largest sag as `sag`.
`on_no_elbow = "warn"` returns that count all the same, with the same
warning, and `"error"` stops.

## Model-aware criteria

Under `criterion = "morans_i"` or `"combined"`, with a `response_var`,
the elbow is supplemented with Moran's I of the cell means at each
candidate \\k\\: the cell means are regressed on the cell-mean
predictors (on an intercept alone when no `predictor_vars` are given)
and Moran's I is computed on the residuals over each cell's 8 nearest
cells. Under `"geometric"`, the default, `response_var` and
`predictor_vars` are not read, whether one or both are given.

Moran's I of the cell means tests whether neighbouring cells are
independent. It does not measure the autocorrelation the cells leave
unexplained in the response, and it does not say how well they represent
it: on an autocorrelated response it grows as the cells get finer than
the range, because neighbouring cells then resemble each other, and it
is close to 0 at about 10 to 20 cells whatever the data, where the
neighbour graph is nearly complete.

`"morans_i"` prefers the candidate with the smallest \\\|z\|\\: the
count at which neighbouring cell means stop resembling each other, cells
about as wide as the range. Only the elbow's neighbourhood is scored, so
the response moves the answer by at most a few cells either side of the
elbow. That ranking is by a diagnostic of the cell means and not by a
count rule, and it is returned with a warning of class
`spatialkit_diagnostic_ranking` that says so and gives the count ranked
first, the first count the statistic could be computed at and the counts
scored.
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
refuses the same statistic unless asked with `diagnostic = TRUE`; here
the criterion is asked for by name, so the count is returned and the
warning carries the caveat.

`"morans_i"` ranks on the standardised deviate \\\|z\| = \|I - E\[I\]\|
/ \mathrm{sd}(I)\\, using the Cliff and Ord regression residual moments,
and not on \\\|I\|\\: both \\E\[I\]\\ and \\Var\[I\]\\ depend on the
number of cells, so \\\|I\|\\ falls as \\k\\ grows whether or not the
finer cells capture anything, and a ranking by \\\|I\|\\ would prefer
the largest candidate for arithmetic reasons alone. The cell-level
residuals are OLS residuals by construction, which is the case those
moments are derived for, but the derivation also assumes errors of equal
variance, and a cell mean over \\n_j\\ points has a variance
proportional to \\1/n_j\\, so \\z\\ is less reliable where single-point
cells sit next to cells of many points. Where structure remains,
\\\|z\|\\ mixes its size with the number of cells it is measured on,
since \\\mathrm{sd}(I)\\ shrinks as cells are added. Both quantities are
reported in the `"diagnostics"` attribute, as `moran_i` and `moran_z`.

Moran's I is computed with an 8-nearest-neighbour weight matrix, so it
carries information only once there are more than nine cells. At nine or
fewer, every cell is a neighbour of every other, the row-standardised
weight matrix is complete, and Moran's I is exactly \\-1/(k - 1)\\ for
any residual vector. The residual moments then give \\E\[I\] = I\\ and
\\\mathrm{Var}\[I\] = 0\\, so the standardised deviate is \\0/0\\: those
candidates return `NA` and are left out of the ranking of `"morans_i"`.
When no candidate in the elbow's neighbourhood clears the floor, the
call falls back to the geometric ranking and logs a warning that says
so. On points with no cluster structure the elbow sits near
\\\sqrt{K\_{max}}\\ and the neighbourhood is the elbow plus or minus
`max(4, top_n)`, so the neighbourhood reaches ten cells only from about
`max_levels = 40`. On such a layer
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
which scores Moran's z at every level of its ladder, is the model-aware
view.

`"combined"` is the elbow's ranking of the counts around it, best sag
first, with Moran's I of the cell means and its \\z\\ beside it in the
diagnostics; Moran's I does not take part in the ranking, and the result
is not an estimate of the number of clusters. The sag is read off the
sweep's WSS curve, the one the elbow is read from, so with an elbow the
first count is the elbow and the counts after it come in order of sag;
they can differ from the two neighbours `"geometric"` returns. An elbow
below ten cells is a count Moran's I cannot score: the geometric ranking
is returned, the elbow and its two neighbours, with a logged warning,
and the diagnostics record it (see Value). When the WSS curve has no
elbow (under `on_no_elbow = "warn"`) the sag is flat, every candidate
tied, and the window comes out in order of distance from the count the
chord rule lands on, which is the geometric order. The k-means re-run at
each candidate, whose cells Moran's I is computed on, can land in
another optimum than the sweep; its WSS is reported as `wss_eval` and
does not rank. Where Moran's I can be computed at no candidate of a
window that reaches ten cells (a constant response, too few cells
holding a row with one, or too few in the half that
`select_on = "split"` reads), `moran_i` and `moran_z` are `NA`, the log
gives the reason at the WARN level, and the counts are ranked as they
would be otherwise.

Use `"geometric"` when the cell count should follow the clustering of
the points, `"morans_i"` for a ranking by the diagnostic, and
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
when the response should drive the count.

## The k-means sweep

The WSS curve is read for its shape, so it is fitted to be smooth. Every
point on it is a k-means local optimum, and with a few random restarts
the curve mixes level effects with optimisation noise and can rise from
one \\k\\ to the next. Each \\k\\ is therefore fitted as the best of 25
restarts seeded by k-means++ (Arthur and Vassilvitskii 2007). A fixed
budget of restarts follows Fränti and Sieranoja (2019), and Steinley
(2003) shows why the usual handful is not enough; 25 is the package's
choice, not a published figure. A curve that still rises somewhere is
reported and not refused. A logged warning names the number of rising
steps, and the model-aware paths return it as `wss_bumps` in the
`"diagnostics"` attribute beside `wss_spread`, the relative spread of
WSS across the restarts at each \\k\\.

The sweep is plain k-means, whatever the cells are built with
afterwards. This function has no `partition` argument: the elbow is a
reading of the sum of squares that plain k-means minimises, and the
cells behind Moran's z are that partition's. A count returned here says
nothing about density-weighted cells;
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
whose cells are of that kind by default, scores counts of those.

To keep memory use and run time bounded for large `max_levels`, the
sweep records only the WSS at each \\k\\, without the cluster
assignments. Moran's I is then evaluated on demand: k-means is re-run
only for a neighbourhood around the elbow (plus or minus 4 by default,
or plus or minus `top_n` if larger). The re-run at each count draws its
starts from a stream of its own, seeded before the sweep, so for one
`set_seed` a count gets the same cells, `wss_eval` and Moran's I in any
neighbourhood that holds it, whichever `top_n` and `max_levels` set the
neighbourhood: a larger `top_n` adds counts to the neighbourhood and
leaves the scores of the others as they were.

## Post-selection inference

When the selection reads the response (here, under
`criterion = "morans_i"` with a response), everything estimated
afterwards on the chosen cells is estimated on data that already
influenced the choice, and its standard errors are post-selection ones:
descriptive, not at nominal coverage (Gao, Bien and Witten 2022; Chen
and Witten 2023 give the exact selective test for two k-means clusters,
the first link of this chain). The exposure is narrower than "the
partition was chosen on the response": every partition here is k-means
on the coordinates alone, and the response only decides which count is
ranked first. It is not zero, because the count determines every cell
the downstream standard errors are computed over. Under `"combined"` the
response reaches the diagnostics and not the ranking, which is the same
whether or not Moran's I can be computed.

`select_on = "split"` is sample splitting: the layer is cut into two
spatially blocked halves
([`make_folds`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)`(k = 2, method = "block_kfold")`),
the criteria that read the response (Moran's I on the cell means) read
the first half only, and the row positions of both halves come back in
the `"split"` attribute (`selection` and `estimation`). The WSS curve
and the k-means cells still use every point: they read coordinates
alone, and the count is for a tessellation of every point, so it is
chosen on that layer's extent and clusters rather than on half of them.
Build the tessellation on every point (cells are geometry), but
aggregate and fit on `data_sf[attr(x, "split")$estimation, ]`, whose
response the selection never saw. That keeps the selection's use of the
response out of the estimates, but only as far as the two halves are
independent: the split has no buffer, so points near the border between
them are correlated with the selection half over the autocorrelation
range, and coverage is nominal only when that range is short against the
blocks. The price is precision: half the points estimate, and García
Rasines and Young (2023) show a contiguous spatial half is less
efficient than the exchangeable split the i.i.d. theory assumes, because
the two halves are not interchangeable. Two alternatives keep the whole
sample: data thinning for count responses (Neufeld et al. 2024) and data
fission for Gaussian-like ones (Leiner et al. 2023). Neither is
implemented here; the split needs no distributional assumption. Under
`"combined"` the split changes only the rows the diagnostics read.
Selection on coordinates alone (`"geometric"`, whatever variables are
given) is not exposed in this way, and `"split"` then changes nothing
but the attribute.

The estimation rows cover only the estimation half's blocks, while the
cells cover the whole layer, so the estimation rows do not fill every
cell. A cell inside the selection half gets none and comes back `NA`; a
cell across the border between the halves is estimated from its
estimation-half points alone, and its standard error describes that
part, not the cell. Count, per cell, how many of its points are
estimation rows, and read inferential results only off cells whose
points all are.

## References

Arthur, D. and Vassilvitskii, S. (2007). k-means++: the advantages of
careful seeding. *Proceedings of the 18th Annual ACM-SIAM Symposium on
Discrete Algorithms*, 1027–1035.

Fränti, P. and Sieranoja, S. (2019). How much can k-means be improved by
using better initialization and repeats? *Pattern Recognition*, 93,
95–112.
[doi:10.1016/j.patcog.2019.04.014](https://doi.org/10.1016/j.patcog.2019.04.014)

Steinley, D. (2003). Local optima in K-means clustering: what you don't
know may hurt you. *Psychological Methods*, 8(3), 294–304.
[doi:10.1037/1082-989X.8.3.294](https://doi.org/10.1037/1082-989X.8.3.294)

Gao, L. L., Bien, J. and Witten, D. (2022). Selective inference for
hierarchical clustering. *Journal of the American Statistical
Association*, 119, 332–342.
[doi:10.1080/01621459.2022.2116331](https://doi.org/10.1080/01621459.2022.2116331)

Chen, Y. T. and Witten, D. M. (2023). Selective inference for k-means
clustering. *Journal of Machine Learning Research*, 24(152), 1–41.
<https://jmlr.org/papers/v24/22-0371.html>

García Rasines, D. and Young, G. A. (2023). Splitting strategies for
post-selection inference. *Biometrika*, 110(3), 597–614.
[doi:10.1093/biomet/asac070](https://doi.org/10.1093/biomet/asac070)

Leiner, J., Duan, B., Wasserman, L. and Ramdas, A. (2023). Data fission:
splitting a single data point. *Journal of the American Statistical
Association*, 120(549), 135–146.
[doi:10.1080/01621459.2023.2270748](https://doi.org/10.1080/01621459.2023.2270748)

Neufeld, A., Dharamshi, A., Gao, L. L. and Witten, D. (2024). Data
thinning for convolution-closed distributions. *Journal of Machine
Learning Research*, 25(57), 1–35.
<https://jmlr.org/papers/v25/23-0446.html>

## See also

[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
and
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
which accept this function's result directly as `approx_n_cells` and `n`
(the first candidate is used, and the output records that it came from
here; an `NA` from it stops them with a note on where it came from);
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
and
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
for the steps that follow.

Other resolution:
[`plot.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.resolution_profile.md),
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md),
[`summary.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.resolution_profile.md)

## Examples

``` r
library(sf)
set.seed(1)
# Two clearly separated clusters: the elbow should sit near k = 2
pts <- st_as_sf(
  data.frame(x = 5e5 + c(runif(25, 0, 10), runif(25, 90, 100)),
             y = 5e6 + c(runif(25, 0, 10), runif(25, 90, 100))),
  coords = c("x", "y"), crs = 32632
)
determine_optimal_levels(pts, max_levels = 6)   # 2 1 3: the elbow first
#> [1] 2 1 3
#> attr(,"structured")
#> [1] TRUE
#> attr(,"sag")
#> [1] 3.046636

# Points spread evenly have no elbow: NA, with structured = FALSE and a
# warning of class spatialkit_no_elbow that says why (printed here
# wrapped, in place of the warning)
even <- st_as_sf(data.frame(x = 5e5 + runif(200, 0, 100),
                            y = 5e6 + runif(200, 0, 100)),
                 coords = c("x", "y"), crs = 32632)
withCallingHandlers(
  determine_optimal_levels(even, max_levels = 6),
  spatialkit_no_elbow = function(w) {
    writeLines(strwrap(conditionMessage(w), 76))
    invokeRestart("muffleWarning")
  }
)
#> determine_optimal_levels(): the WSS curve has no elbow: on log-log axes it
#> falls in a straight line, as it does for points with no cluster structure.
#> k = 3 is where the chord rule lands on such a curve, set by the ladder (k =
#> 1 to 6, max_levels) rather than by the data. Choose the count on other
#> grounds, e.g. resolution_profile() with a response. Returning NA
#> (on_no_elbow = "na").
#> [1] NA
#> attr(,"structured")
#> [1] FALSE
#> attr(,"sag")
#> [1] 0
```
