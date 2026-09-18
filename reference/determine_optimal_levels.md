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

  An sf object.

- max_levels:

  Integer upper bound on levels. Default 12.

- top_n:

  Integer; how many candidates to return. Default 3. Under
  `criterion = "geometric"` the candidate set is the elbow and its two
  immediate neighbours, so at most 3 values are ever returned no matter
  how large `top_n` is; only the model-aware criteria can return more.

- sample_n:

  Integer; subsample size for speed. Default 1500.

- set_seed:

  Integer RNG seed. Default 123.

- response_var:

  Optional response column name. When provided alongside
  `predictor_vars`, enables model-aware level selection via Moran's I on
  OLS residuals. Must be numeric or logical (logicals are read as 0/1);
  a factor or character response raises an error rather than being
  coerced, because the residuals of an OLS fit to arbitrary level codes
  carry no meaning to test for autocorrelation.

- predictor_vars:

  Optional predictor column names. Must be numeric or logical (logicals
  are read as 0/1); factor/character columns raise an error.

- criterion:

  One of `"geometric"` (default when no response given), `"morans_i"`
  (select the k whose residual Moran's I is least *significant*), or
  `"combined"` (rank-average of WSS elbow distance and that same
  quantity). Falls back to `"geometric"` if response/predictors are
  unavailable, and also when no candidate clears the nine-cell
  resolution floor described in **Details**. Note that supplying both
  `response_var` and `predictor_vars` upgrades `"geometric"` to
  `"combined"`: the selection then depends on the response (see
  "Post-selection inference").

- select_on:

  `"all"` (default) selects on every point; `"split"` selects on one
  spatially blocked half of the points and returns the other half as the
  set to estimate on, so that the standard errors computed downstream on
  the chosen cells are not post-selection. See "Post-selection
  inference".

## Value

An integer vector of candidate level counts, **best first**: under the
geometric criterion the elbow, then its lower and upper neighbours;
under the model-aware criteria the candidates in rank order. `k[1]` is
therefore the top-ranked count on every path, and `top_n = 1` returns it
alone. When `criterion != "geometric"`, an attribute `"diagnostics"` is
attached with per-k Moran's I values (`moran_i`) and their standardised
deviates (`moran_z`), the WSS curve (`wss`) with the relative
between-restart spread at each `k` (`wss_spread`), the number of rising
steps on it (`wss_bumps`), the restart budget (`nstart`), the geometric
elbow the evaluated neighbourhood was drawn around (`knee_k`) and the
`k` at which k-means failed (`failed_k`; their `wss` entries are
interpolated from the neighbours, not measured) — except when the
model-aware path itself falls back to the geometric result (no viable k
in the elbow neighbourhood, or Moran's I could not be computed for any
candidate), in which case no diagnostics are available and the attribute
is absent. Both fallbacks are logged as warnings. The geometric path
returns a plain integer vector; a rising WSS curve is still logged
there. With `select_on = "split"` every path adds a `"split"` attribute:
a list with `selection` and `estimation` (integer row positions in
`data_sf`), `method` and `seed`. For a full per-level table — criteria,
cell support, restart spread, the flat region — see
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md).

## Details

When `response_var` and `predictor_vars` are provided, the geometric WSS
elbow is supplemented with Moran's I computed on OLS residuals at each
candidate k. The Moran's I profile measures how much spatial
autocorrelation in the response remains *unexplained* at a given
tessellation resolution — a direct reflection of the spatial process
being modeled, rather than mere geometric compactness of coordinates.
The combined criterion selects the k that best balances geometric
parsimony and residual spatial independence.

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
no increase at all. A curve that still rises somewhere is reported — a
logged warning names the number of rising steps, and the model-aware
paths return it as `wss_bumps` in the `"diagnostics"` attribute beside
`wss_spread`, the relative spread of WSS across the restarts at each `k`
— but not refused: a bumpy curve is uncertain, not unidentified. Because
the optimiser changed, a selection made by an earlier version on a curve
that had such a bump can differ from the one made now; where the earlier
curve was clean, the answer is the same.

**The model-aware criteria rank on the standardised deviate, not on
\|Moran's I\|.** Both \\E\[I\]\\ and \\Var\[I\]\\ depend on the number
of cells, so \\\|I\|\\ falls as `k` grows whether or not the finer
tessellation is capturing anything. Measured over 300 replicates of a
response with *no* spatial structure, mean \\\|I\|\\ fell monotonically
from 0.114 at `k = 10` to 0.050 at `k = 60` (\\-56\\\\), which made an
\\\|I\|\\ ranking prefer the largest candidate for arithmetic reasons
alone. Candidates are therefore ordered by \\\|z\| = \|I - E\[I\]\| /
\mathrm{sd}(I)\\ using the Cliff & Ord regression residual moments —
exact here, because the cell-level residuals are OLS residuals by
construction. Over the same runs \\z\\ had mean \\\approx 0\\,
\\\mathrm{sd} \approx 1\\ and a two-sided 5\\ 0.040–0.057 at every `k`.
Both quantities are reported in the `"diagnostics"` attribute, as
`moran_i` and `moran_z`.

**Resolution floor on the model-aware criteria.** Moran's I is computed
on cell-level residuals with an 8-nearest-neighbour weight matrix, so it
only carries information once there are more than nine cells. At nine or
fewer, every cell is a neighbour of every other, the row-standardised
weight matrix is complete, and Moran's I collapses to exactly \\-1/(k -
1)\\ for *any* residual vector — a function of `k` alone. The criterion
ranks on \\\|z\|\\, not on \\\|I\|\\, and at the floor the residual
moments give \\E\[I\] = I\\ and \\\mathrm{Var}\[I\] = 0\\ identically
(the algebra holds to \\10^{-16}\\), so the standardised deviate is
\\0/0\\: it carries no information about the tessellation, and whichever
way rounding noise resolves it those candidates would rank first or last
on nothing. They therefore return `NA` and are excluded from the
model-aware ranking. When no candidate in the elbow neighbourhood clears
the floor — which is the usual outcome for small `max_levels` — the
whole call falls back to the geometric ranking and logs a warning; raise
`max_levels` above roughly 10 if you want the model-aware criteria to
contribute. Under `criterion = "combined"`, a candidate below the floor
that sits alongside candidates above it is ranked last on the Moran's I
axis while still competing on the geometric axis.

## Post-selection inference

When the selection reads the response — here, whenever both
`response_var` and `predictor_vars` are supplied — everything estimated
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
the selection runs on the first half only, and the row positions of both
halves come back in the `"split"` attribute (`selection` and
`estimation`). Build the tessellation on every point — cells are
geometry — but aggregate and fit on
`data_sf[attr(x, "split")$estimation, ]`, which the selection never saw;
that restores nominal coverage with no new theory. The price is
precision: half the points estimate, and García Rasines and Young (2023)
show a *contiguous* spatial half is less efficient than the exchangeable
split the i.i.d. theory assumes, because the two halves are not
interchangeable. Two alternatives keep the whole sample — data thinning
for count responses (Neufeld et al. 2024) and data fission for
Gaussian-like ones (Leiner et al. 2023) — and are not implemented here;
the split needs no distributional assumption, which is why it comes
first. Selection on coordinates alone (`"geometric"` with no response)
is not exposed in this way, and `"split"` then changes nothing but the
attribute.

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
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)

## Examples

``` r
library(sf)
set.seed(1)
# Two clearly separated clusters: the elbow should sit near k = 2
pts <- st_as_sf(
  data.frame(x = c(runif(25, 0, 10), runif(25, 90, 100)),
             y = c(runif(25, 0, 10), runif(25, 90, 100))),
  coords = c("x", "y"), crs = 32632
)
determine_optimal_levels(pts, max_levels = 6)   # 2 1 3: the elbow first
#> [1] 2 1 3
```
