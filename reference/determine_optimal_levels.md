# Determine an optimal number of spatial levels via an elbow heuristic

Computes a WSS curve over k=1..K_max using k-means on projected feature
coordinates and selects candidate k values around the elbow.

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
  select_on = c("all", "split")
)
```

## Arguments

- data_sf:

  An sf object. Features with empty or non-finite coordinates are
  dropped with a warning.

- max_levels:

  Integer upper bound on levels. Default 12. The sweep also stops at the
  number of distinct locations (k-means cannot place more centres) and
  one short of the number of points.

- top_n:

  Integer; how many candidates to return. Default 3. Under
  `criterion = "geometric"` the candidate set is the elbow and its two
  immediate neighbours, so at most 3 values are ever returned no matter
  how large `top_n` is; only the model-aware criteria can return more.

- sample_n:

  Integer; subsample size for speed. Default 1500.

- set_seed:

  Integer RNG seed for the subsample and the k-means++ restarts;
  restored afterwards. Default 123. The rows are put in coordinate order
  before either, so the answer does not depend on the order they come
  in.

- response_var:

  Optional response column name. When provided alongside
  `predictor_vars`, enables model-aware level selection via Moran's I on
  OLS residuals. Must be numeric or logical (logicals are read as 0/1);
  a factor or character response raises an error and is never coerced,
  because the residuals of an OLS fit to arbitrary level codes carry no
  meaning to test for autocorrelation. Rows where it, or a predictor, is
  missing or non-finite stay in the WSS sweep and the cells and are left
  out of Moran's I; a logged warning gives their number.

- predictor_vars:

  Optional predictor column names. Must be numeric or logical (logicals
  are read as 0/1); factor/character columns raise an error.

- criterion:

  One of `"geometric"` (default when no response given), `"morans_i"`
  (select the k whose residual Moran's I is least *significant*), or
  `"combined"` (rank-average of the WSS curve's log-log sag, the
  quantity the elbow is read from, and that same significance). Falls
  back to `"geometric"`, with a warning, when `response_var` or
  `predictor_vars` is not given, and with a logged warning when no
  candidate clears the nine-cell resolution floor described in
  **Details**. A `response_var` or `predictor_vars` naming a column that
  is not in `data_sf` is an error. Supplying both `response_var` and
  `predictor_vars` upgrades `"geometric"` to `"combined"`: the selection
  then depends on the response (see "Post-selection inference").

- select_on:

  `"all"` (default) selects on every point; `"split"` reads the response
  on one spatially blocked half of the points only and returns the other
  half as the set to estimate on, so that the standard errors computed
  downstream on the chosen cells are not post-selection. The count is
  still chosen for the whole layer. See "Post-selection inference".

## Value

An integer vector of candidate level counts, **best first**: under the
geometric criterion the elbow, then its lower and upper neighbours (with
a warning when the WSS curve has no elbow and the first is the ladder's
choice; see Details); under the model-aware criteria the candidates in
rank order. `k[1]` is therefore the top-ranked count on every path, and
`top_n = 1` returns it alone. When `criterion != "geometric"`, an
attribute `"diagnostics"` is attached with per-k Moran's I values
(`moran_i`) and their standardised deviates (`moran_z`), the WSS curve
(`wss`) with the relative between-restart spread at each `k`
(`wss_spread`), the number of rising steps on it (`wss_bumps`), the
restart budget (`nstart`), the geometric elbow the evaluated
neighbourhood was drawn around (`knee_k`), the `k` at which k-means
failed (`failed_k`; their `wss` entries are interpolated from the
neighbours, not measured) and the `k` the model-aware pass actually
scored (`eval_ks`, the elbow's neighbourhood). Under `"combined"` it
also carries the WSS of the re-run clustering at those `k` (`wss_eval`),
the rank average that ordered them (`combined_rank`, named by `k`) and
`criterion = "combined"`; when `"combined"` returned the geometric
ranking because the elbow is below ten cells, it carries
`criterion = "geometric"` and `fallback` (the reason) in place of
`combined_rank`. When the model-aware path itself falls back to the
geometric result (no viable k in the elbow neighbourhood, or Moran's I
could not be computed for any candidate), no diagnostics are available
and the attribute is absent. Both fallbacks are logged as warnings. The
geometric path returns a plain integer vector; a rising WSS curve is
still logged there. With `select_on = "split"` every path adds a
`"split"` attribute: a list with `selection` and `estimation` (integer
row positions in `data_sf`), `method` and `seed`. For a full per-level
table of criteria, cell support, restart spread and the flat region, see
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md).

## Details

**The elbow is read on log-log axes, and there may be none.** Points
with no cluster structure have a WSS curve close to \\c/k\\, and the
classical rule, the point furthest below the chord from the first to the
last k, still finds a "knee" on it on linear axes, at about
\\\sqrt{K\_{max}}\\: 4 at the default `max_levels = 12` and 13 at 160,
whatever the data. On \\\log k\\ against \\\log\\ WSS that curve is a
straight line, while separated clusters fall faster than it until there
is one cell per cluster and like it after, a bend at the cluster count.
The elbow is therefore the k whose \\\log\\ WSS sags furthest below the
straight line joining k = 1 and \\K\_{max}\\ on those axes, and it
counts as one only when the sag is at least \\\log 1.25\\ (the WSS a
fifth below the power law through the ends). Measured on 60 to 1500
points with ladders to 3–40, uniform layouts over squares, discs,
triangles, an L-shape, density gradients and jittered lattices sagged at
most 0.16, and two to ten separated clusters at least 0.6 once
`max_levels` passed the cluster count. With no elbow the function warns
and returns the linear-axis answer, which the ladder chose, not the
data. An elongated extent also bends, at about its aspect ratio, because
the first cuts go across its long axis (0.2–0.45 for a 4:1 rectangle);
past the threshold that bend is reported as an elbow, and it describes
the extent's shape rather than clusters in it. When locations repeat
(stations visited many times), the sweep can reach one cell per distinct
location, where the WSS is zero (to within \\10^{-12}\\ of the total).
That `k` has no place on log axes and is left out of the line, but the
fall to zero is the sharpest bend there is: when the rest of the curve
has no elbow, the number of distinct locations is the elbow.

When `response_var` and `predictor_vars` are provided, the geometric WSS
elbow is supplemented with Moran's I computed on OLS residuals at each
candidate k. The Moran's I profile measures how much spatial
autocorrelation in the response remains *unexplained* at a given
tessellation resolution. It is a direct reflection of the spatial
process being modeled instead of the mere geometric compactness of
coordinates. The combined criterion selects the k that best balances
geometric parsimony and residual spatial independence.

To keep memory use and runtime bounded for large `max_levels`, the
initial k-means sweep records only within-cluster sum-of-squares (WSS)
without retaining cluster assignments. Moran's I is then evaluated
lazily: k-means is re-run only for a focused neighbourhood around the
elbow (±4 by default, or ±`top_n` if larger), so that only the most
promising candidate k values incur the cost of the full Moran's I
computation.

**The WSS curve is read for its shape, so it is fitted to be smooth.**
Every point on it is a k-means local optimum, and with a few random
restarts the curve mixes level effects with optimisation noise: measured
on clustered layouts, a sweep over `k = 1..30` at
`stats::kmeans(nstart = 5)` *rose* at one or two steps in three of five
draws, and an earlier form of the elbow rule once selected such a bump.
Each `k` is therefore fitted as the best of 25 restarts seeded by
k-means++ (Arthur and Vassilvitskii 2007), the budget at which the gain
from further restarts saturates (Fränti and Sieranoja 2019; Steinley
2003 on why the usual handful is not enough). The same sweeps then had
no increase at all. A curve that still rises somewhere is reported but
not refused: a bumpy curve is uncertain, not unidentified. A logged
warning names the number of rising steps, and the model-aware paths
return it as `wss_bumps` in the `"diagnostics"` attribute beside
`wss_spread`, the relative spread of WSS across the restarts at each
`k`. Because the optimiser changed, a selection made by an earlier
version on a curve that had such a bump can differ from the one made
now; where the earlier curve was clean, the selection is the same.

**The model-aware criteria rank on the standardised deviate, not on
\|Moran's I\|.** Both \\E\[I\]\\ and \\Var\[I\]\\ depend on the number
of cells, so \\\|I\|\\ falls as `k` grows whether or not the finer
tessellation is capturing anything. Measured over 300 replicates of a
response with *no* spatial structure, mean \\\|I\|\\ fell monotonically
from 0.114 at `k = 10` to 0.050 at `k = 60` (\\-56\\\\), which made an
\\\|I\|\\ ranking prefer the largest candidate for arithmetic reasons
alone. Candidates are therefore ordered by \\\|z\| = \|I - E\[I\]\| /
\mathrm{sd}(I)\\ using the Cliff & Ord regression residual moments. The
cell-level residuals are OLS residuals by construction, which is the
case those moments are derived for, but the derivation also assumes
errors of equal variance, and a cell mean over \\n_j\\ points has a
variance proportional to \\1/n_j\\. Over the same runs \\z\\ had mean
\\\approx 0\\, \\\mathrm{sd} \approx 1\\ and a two-sided 5% rejection
rate of 0.040–0.057 at every `k`, and it stayed calibrated on gradient
and moderately clustered layouts; with single-point cells next to cells
of 70 or more points its mean rose to 0.2–0.34 and its rejection rate to
7–8% at 20 cells. Where structure remains, \\\|z\|\\ mixes its size with
the number of cells it is measured on, since \\\mathrm{sd}(I)\\ shrinks
as cells are added. Both quantities are reported in the `"diagnostics"`
attribute, as `moran_i` and `moran_z`.

**Resolution floor on the model-aware criteria.** Moran's I is computed
on cell-level residuals with an 8-nearest-neighbour weight matrix, so it
only carries information once there are more than nine cells. At nine or
fewer, every cell is a neighbour of every other, the row-standardised
weight matrix is complete, and Moran's I collapses to exactly \\-1/(k -
1)\\ for *any* residual vector (a function of `k` alone). The criterion
ranks on \\\|z\|\\, not on \\\|I\|\\, and at the floor the residual
moments give \\E\[I\] = I\\ and \\\mathrm{Var}\[I\] = 0\\ identically
(the algebra holds to \\10^{-16}\\), so the standardised deviate is
\\0/0\\: it carries no information about the tessellation, and whichever
way rounding noise resolves it those candidates would rank first or last
on nothing. They therefore return `NA` and are excluded from the
model-aware ranking. When no candidate in the elbow neighbourhood clears
the floor, the whole call falls back to the geometric ranking and logs a
warning that says so. That is the usual outcome well past
`max_levels = 10`: the neighbourhood is the elbow plus or minus
`max(4, top_n)`, and on points with no cluster structure the elbow sits
near \\\sqrt{K\_{max}}\\, so the neighbourhood reaches ten cells only
from about `max_levels = 40`. Measured on 1000 uniform points,
`max_levels` of 12, 20 and 30 all fell back, and 40 scored `k` = 10 and
11 alone. On such a layer
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
which scores Moran's z at every level of its ladder, is the model-aware
view. Under `criterion = "combined"`, an elbow below ten cells is a
count Moran's I cannot score, so it cannot weigh the response there: the
geometric ranking is returned, with a logged warning, and the
diagnostics record it (see Value). (Ranking the window anyway put the
smallest count Moran's I scores, ten, first whatever the response did.)
Otherwise the candidates below the floor that sit alongside candidates
above it all take the last place on the Moran's I axis, after every
candidate it scored, while still competing on the geometric axis. That
axis is the elbow's own log-log sag at each candidate; when the WSS
curve has no elbow it is flat, every candidate tied, and Moran's I alone
orders them.

`"combined"` is not an estimate of the number of clusters. Where the
elbow is at ten cells or more, Moran's I can move the pick away from it
when the response is still spatially structured at the elbow's
resolution. Use `"geometric"` when the cell count should follow the
clustering of the points, and
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
when the response should drive it.

## Post-selection inference

When the selection reads the response (here, whenever both
`response_var` and `predictor_vars` are supplied), everything estimated
afterwards on the chosen cells is estimated on data that already
influenced the choice, and its standard errors are post-selection ones:
descriptive, not at nominal coverage (Gao, Bien and Witten 2022; Chen
and Witten 2023 give the exact selective test for two k-means clusters,
the first link of this chain). The exposure is narrower than "the
partition was chosen on the response": every partition here is k-means
on the coordinates alone, and the response only decides which *count* is
ranked first. It is not zero, because the count determines every cell
the downstream standard errors are computed over.

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
Rasines and Young (2023) show a *contiguous* spatial half is less
efficient than the exchangeable split the i.i.d. theory assumes, because
the two halves are not interchangeable. Two alternatives keep the whole
sample: data thinning for count responses (Neufeld et al. 2024) and data
fission for Gaussian-like ones (Leiner et al. 2023). Neither is
implemented here; the split needs no distributional assumption, which is
why it comes first. Selection on coordinates alone (`"geometric"` with
no response) is not exposed in this way, and `"split"` then changes
nothing but the attribute.

The estimation rows cover only the estimation half's blocks, while the
cells cover the whole layer, so the estimation rows do not fill every
cell. A cell inside the selection half gets none and comes back `NA`; a
cell across the border between the halves is estimated from its
estimation-half points alone, and its standard error describes that
part, not the cell (on 400 simulated fields a nominal 95% interval
covered 0.97 in cells inside the estimation half and 0.88 in cells
across the border). Count, per cell, how many of its points are
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
here);
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
and
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
for the steps that follow.

Other aggregation:
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md),
[`kriging_adequacy()`](https://elkronos.github.io/gis_modeling_toolkit/reference/kriging_adequacy.md),
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md),
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md),
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
```
