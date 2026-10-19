# Choosing a resolution

*This article needs two optional packages: **gstat**, which fits the
variogram behind the autocorrelation range, and **ggplot2** for the
figures. When either is missing the code is shown but not run, and a
note at the top says so.*

## The question

Before point observations can be aggregated into regions, something has
to decide how many regions. Twenty cells over a county smooth away most
of the variation you came to measure. Two thousand give you cells
holding one observation each, where the standard errors are undefined
and every cell mean is one number pretending to be an average.

Administrative boundaries settle this by fiat: you get the census tracts
that exist. Drawing the regions yourself removes the fiat and leaves the
decision with you, and the number propagates into the variance of the
cell means, the design effect, and the size of the blocks a spatial
cross-validation holds out.

Three functions answer the question, in increasing order of how much
they tell you.
[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
reads an elbow off the coordinates alone and runs on the hard
dependencies; it is the quick answer and the one
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
accepts directly.
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
scores every cell count on a ladder in five columns at once, four count
criteria and one diagnostic, and needs gstat for two of them.
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
reads one criterion off that profile, and
[`summary()`](https://rdrr.io/r/base/summary.html) on the profile reads
the count criteria side by side. If the cells will feed a model, the
profile is the one to use. If you need a count now and the data are all
you have, the elbow is defensible when the points cluster, and this
article says where it falls short. On points spread evenly there is no
elbow to read, and both functions say so: the profile leaves its `elbow`
column empty, and
[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
warns that the count it still returns was chosen by the ends of its
ladder, not by the data.

Called with the response and nothing else,
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
takes its defaults: one level at 4.5 observations per cell, scored for
the pooled value that
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
reports by default, in density-weighted cells, with no variogram fitted.
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
then reads the count off that size (“Setting a size” below). A ladder
argument such as `n_levels` gives the ladder this article reads first.
Each earlier default is one argument away:
`estimate = NULL, obs_per_cell = NULL` in the profile for the ladder
without the pooled value’s columns, `sac = NULL` for a variogram on a
ladder of sizes, `partition = "kmeans"` for k-means cells from 25
starts, `criterion = "cp"` in
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md),
and `estimate = "mean"` in
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md).

The argument that says “how many” is spelled differently by the function
it belongs to: `max_levels` bounds the elbow’s ladder, `n_levels` sets
the profile’s and `obs_per_cell` sets it in observations per cell,
`approx_n_cells` and `target_cells` are what the grid builders take, and
`n` is what
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
takes. Each of the last three accepts the object the first two return.

## A fixture with known structure

A simulated exponential field, so the structure is known before anything
is fitted. The covariance parameter `a = 100` gives an effective range
of 300 units across a 1000-unit square, and `sd = 1` adds a nugget large
enough to matter later.

``` r

library(sf)
library(spatialkit)

set.seed(42)
n  <- 500
xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
D  <- as.matrix(dist(xy))
xy$z <- as.numeric(t(chol(exp(-D / 100) + diag(1e-8, n))) %*% rnorm(n)) +
        rnorm(n, sd = 1)
pts  <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
```

## The short answer

[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
returns a few candidate cell counts, best first.

``` r

lv <- determine_optimal_levels(pts, max_levels = 40)
lv
```

    ## [1] 7 6 8

This fixture’s points are uniform, so their within-cluster sum of
squares falls like `c / k` all the way down and bends nowhere of its
own. The call warns that it found no elbow and that `max_levels` chose
`lv`: a larger `max_levels` would give a larger answer. Read an elbow
only off points that cluster; here the profile below is the better
guide.

`lv[1]` is what
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
takes as `n`, and what
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
takes as `approx_n_cells` for its `"hex"` and `"square"` methods.
Voronoi has no cell count to set: it grows one cell per point you hand
it, which is why the count is set when you place the seeds rather than
when you tessellate them. The last section shows the two calls in order.

Two things to know about this function before you rely on it. Its
`criterion = "morans_i"` and `criterion = "combined"` settings rank the
counts around the elbow by Moran’s z of the cell means, which is a
diagnostic and not a count rule; the call returns the ranking with a
warning of class `spatialkit_diagnostic_ranking` that says so. They need
a `response_var`; without one the call warns and falls back to
`"geometric"`. And its ladder runs from 1 to `max_levels` with no lower
bound from the spatial correlation of the data, so it can prefer cells
wider than the field’s own correlation range. The profile below reports
that bound, and
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
reads `reliability` from it up, which is why the two can disagree.

## The profile

[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
fits every level on a ladder and reports what each one costs, so you can
see the shape of the trade-off before picking a point on it. `n_levels`
asks for that ladder. Levels are spaced logarithmically, because cell
diameter falls as $`L^{-1/2}`$. At its defaults the profile also scores
the leave-one-out error of the pooled value, `loo`, at every level, and
draws the cells by density-weighted k-means.

``` r

prof <- resolution_profile(pts, response_var = "z", n_levels = 16)
prof
```

    ## Resolution profile: 16 levels on 500 points
    ##   ladder      : 2 to 55 cells (floor 17 from range 243, not applied, ceiling 55 from min_cell_n = 9)
    ##   variogram   : nugget 0.702, partial sill 1.42, range 243
    ##   scored on   : response; 3 density-weighted k-means++ starts (power 1, 8 neighbours) per level; WSS rises at 0 step(s)
    ##   one cell    : rss 1040, cp 2.08 (one mean for the layer; a level has to beat it)
    ##   estimate    : pooled toward the 8 nearest cells; loo on 500 rows of the layer; one cell: loo 2.08
    ##   elbow       : none; the WSS curve falls as it does with no cluster structure
    ## 
    ##  levels      wss wss_spread elbow cell_n_min cell_n_median cell_diam_median
    ##       2 52800000   0.036100    NA        228         250.0            649.0
    ##       3 32100000   0.002020    NA        150         164.0            485.0
    ##       4 21500000   0.000877    NA        110         128.0            410.0
    ##       5 17900000   0.036300    NA         71         103.0            365.0
    ##       6 14700000   0.032500    NA         71          82.0            342.0
    ##       8 10400000   0.005010    NA         54          59.5            289.0
    ##       9  8820000   0.006290    NA         42          56.0            273.0
    ##      12  6430000   0.048700    NA         31          40.5            222.0
    ##      15  5140000   0.113000    NA         24          33.0            194.0
    ##      17  4750000   0.068800    NA         20          28.0            186.0
    ##      18  4480000   0.120000    NA         15          28.0            187.0
    ##      23  3250000   0.018900    NA         15          21.0            158.0
    ##      28  2590000   0.074300    NA          9          17.0            138.0
    ##      35  2010000   0.044900    NA          9          14.0            123.0
    ##      44  1580000   0.042700    NA          4          11.0            108.0
    ##      55  1210000   0.074100    NA          4           9.0             94.4
    ##   rss   cp   moran_i moran_z reliability cells_resp cell_resp_mean
    ##  1020 2.05        NA      NA       0.760          2         250.00
    ##  1000 2.01        NA      NA       0.800          3         167.00
    ##  1010 2.03        NA      NA       0.811          4         125.00
    ##   945 1.90        NA      NA       0.814          5         100.00
    ##   929 1.88        NA      NA       0.814          6          83.30
    ##   945 1.91        NA      NA       0.811          8          62.50
    ##   911 1.85        NA      NA       0.808          9          55.60
    ##   912 1.86 -0.030200  0.9710       0.799         12          41.70
    ##   823 1.69 -0.062700  0.1200       0.790         15          33.30
    ##   824 1.69  0.000359  0.8910       0.783         17          29.40
    ##   772 1.59 -0.064300 -0.0754       0.780         18          27.80
    ##   742 1.55 -0.031700  0.1830       0.765         23          21.70
    ##   773 1.62  0.042400  1.0900       0.750         28          17.90
    ##   682 1.46  0.088000  1.7200       0.732         35          14.30
    ##   663 1.45  0.104000  2.0000       0.710         44          11.40
    ##   655 1.46  0.210000  3.8800       0.686         55           9.09
    ##  scored_per_cell  loo loo_se share_alone own_weight
    ##           250.00 2.07 0.1230           0      0.914
    ##           167.00 2.03 0.1220           0      0.921
    ##           125.00 2.05 0.1220           0      0.844
    ##           100.00 1.93 0.1160           0      0.929
    ##            83.30 1.90 0.1130           0      0.926
    ##            62.50 1.95 0.1160           0      0.863
    ##            55.60 1.89 0.1110           0      0.894
    ##            41.70 1.91 0.1140           0      0.833
    ##            33.30 1.74 0.1090           0      0.901
    ##            29.40 1.76 0.1080           0      0.879
    ##            27.80 1.65 0.0993           0      0.909
    ##            21.70 1.62 0.0975           0      0.889
    ##            17.90 1.72 0.1060           0      0.825
    ##            14.30 1.56 0.0969           0      0.831
    ##            11.40 1.57 0.0953           0      0.819
    ##             9.09 1.60 0.1000           0      0.722

One row per level. `cell_n_median` and `cell_diam_median` are usually
the first two columns worth reading: how many points a typical cell
holds, and how big it is in CRS units. `cell_diam_median` is twice the
median root-mean-square distance of a cell’s points from its centre,
about 0.8 of the side of a square cell of the same area, so compare it
against something you already know about your own data, such as the
spacing of a sampling grid or the size of a field, with that factor in
mind.

### Four criteria and a diagnostic

| criterion | measures | needs | direction |
|----|----|----|----|
| `elbow` | how far the log of the within-cluster sum of squares sags below a power law (the straight log-log line from one cell to the last level); `NA` at every level when the points have no cluster structure | coordinates only | larger |
| `cp` | Mallows’ $`C_p`$ of the piecewise-constant approximation of the response by cell means | a response, and a variogram for the nugget | smaller |
| `moran_z` | Moran’s z of the cell means (of their regression residuals when predictors are given): whether neighbouring cells resemble each other | a response | a diagnostic, not a count rule |
| `reliability` | the share of the spread in the cell means that is signal, not sampling noise | a variogram | larger |
| `loo` | the leave-one-out error of the cell value the profile scores, the pooled value by default (“Setting a size”) | a response | smaller |

`elbow` sees the coordinates and knows nothing about what you measured.
`cp` and `reliability` split the practical question in half: how well
the cells represent the field, and whether the cell values can be told
apart from noise. What `cp` estimates is the error of re-measuring the
sampled sites, each predicted by the mean of its cell, with the
variogram’s nugget as the noise of one measurement. For a new location
the penalty would also need the spread of the field inside a cell, which
the nugget leaves out. `moran_z` is a diagnostic of the cell means and
not a count rule. It grows as the cells get finer than the range,
because neighbouring cells then resemble each other, and it is near zero
at the coarsest counts it is computed at whatever the data. It does not
say how well the cells represent the response, so
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
refuses it unless asked with `diagnostic = TRUE`, and
[`summary()`](https://rdrr.io/r/base/summary.html) leaves it out unless
it is named.

``` r

plot(prof)
```

![Four stacked panels against the number of cells on a shared axis, from
2 to 55; there is no elbow panel because these uniform points have no
elbow. Mallows' Cp falls across the ladder and is lowest at 44 cells.
Reliability rises to a peak near 5 cells and declines after it. The
leave-one-out error of the pooled value falls as Cp does and is lowest
at 35 cells. The absolute Moran's z of the cell means is drawn from 12
cells on, stays near 1 or below up to 23 cells and climbs past 3 at 55.
A red dot and dotted line mark the level selected for Cp, at 44 cells,
for reliability, at 17 cells, and for the leave-one-out error, at 35
cells, with a shaded band over the levels within tolerance of each; the
Moran panel has no selected
level.](resolution_files/figure-html/plot-profile-1.png)

Two behaviours are worth knowing before reading the numbers.

**`cp` needs a nugget to have an interior optimum.** On a smooth field
the piecewise-constant approximation keeps improving as cells shrink,
the penalty is too small to stop it, and the minimum lands wherever
`min_cell_n` stops the ladder. The support ceiling is then doing the
choosing, and
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
flags it.

This fixture has a fitted nugget of 0.7 against a partial sill of 1.42,
which is why `cp` turns around inside the ladder here.

**`reliability` is highest at a coarse count**, because bigger cells
hold more points and average away more noise. Its maximum usually lies
below the range floor, where cells are wider than the range.
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
therefore reads it from the floor up, and what it returns is then the
floor: a bound on the analysis, with the range choosing and not the
curve. The next section is about spotting that.

## Choosing, and the flat region

[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
picks the level a criterion prefers and reports the band around it.

``` r

sel <- select_resolution(prof, criterion = "reliability")
sel
```

    ## Resolution by reliability: 17 cells
    ##   flat region : 17 to 18 (2 of 16 levels)
    ##   cells       : density-weighted k-means++ starts (power 1, 8 neighbours)
    ##   note        : the optimum is the range floor (area / range^2); the bound is
    ##                 choosing, not the criterion. Fewer cells would be wider than
    ##                 the range and average over more than one patch of the field.

`sel$best` is the optimum and `sel$flat` is every level within `tol` of
it. Quote the band when you write the analysis up. `tol` is relative to
the optimum’s value for `cp` and `reliability`, and relative to the
criterion’s range for `elbow` and `moran_z`, whose optima can be zero.

Two flags say when a bound chose instead of the criterion:

``` r

c(at_floor = sel$at_floor, at_ceiling = sel$at_ceiling)
```

    ##   at_floor at_ceiling 
    ##       TRUE      FALSE

`at_floor` is `TRUE` here. Reliability wanted to keep going coarser and
the correlation-range floor stopped it, so read the answer as *at most
this fine* and take the number as a bound on the analysis, not as an
optimum. The floor is a fact about the data, so the fix is a different
question, not a different `tol`.

Criteria can prefer different levels while agreeing on a region.
[`summary()`](https://rdrr.io/r/base/summary.html) on the profile reads
the count criteria at once, and closes with the levels that every band
contains:

``` r

summary(prof)
```

    ## Resolution picks: 3 criteria over 16 levels (2 to 55 cells)
    ## 
    ##   cells: density-weighted k-means++ starts (power 1, 8 neighbours)
    ## 
    ##    criterion best flat region levels in band
    ##           cp   44    35 to 55              3
    ##  reliability   17    17 to 18              2
    ##          loo   35    35 to 44              2
    ## 
    ##   reliability: the optimum is the range floor (area / range^2).
    ##   There the bound is choosing, not the criterion.
    ## 
    ##   picks span 17 to 44 cells (2.6x)
    ##   no level is in every flat region: the criteria disagree over the
    ##   whole ladder. plot() draws the curves they were read from.

A band is a set of levels, not an interval, because the criterion curves
are not monotone: a band that skips a rung prints as a comma-separated
list, while a solid run of 3 rungs prints as a range, as `cp` does with
`35 to 55`. Where the bands overlap you have a defensible set of levels,
and the last line names it; `attr(summary(prof), "common")` returns the
same levels for use in code, and `attr(summary(prof), "bands")` each
criterion’s region in full. Where they do not overlap, the criteria are
answering different questions and you have to say which one your
analysis needs. An empty intersection is a result: it says this field
has no single resolution that satisfies every way of asking.

The table is not a decision procedure, and nothing in the package picks
among its criteria for you. Cross-validating a model at each suggested
level is affordable, but it does not settle the question either: on a
simulated field the level that won on cross-validated $`R^2`$ moved with
the fold seed, and the coarsest grid won most often because its score
had the widest spread, not because it was better.

## What the ladder can support

The ladder runs from 2 cells to a ceiling the data impose. The profile
records that ceiling and, beside it, a floor from the autocorrelation
range.

``` r

str(attr(prof, "bounds"))
```

    ## List of 15
    ##  $ floor       : int 17
    ##  $ ceiling     : int 55
    ##  $ ceiling_from: chr "min_cell_n"
    ##  $ supported   : logi TRUE
    ##  $ area        : num 976371
    ##  $ range       : num 243
    ##  $ n           : int 500
    ##  $ n_sample    : int 500
    ##  $ n_distinct  : int 500
    ##  $ min_cell_n  : int 9
    ##  $ range_floor : logi FALSE
    ##  $ n_response  : int 500
    ##  $ n_scored    : int 500
    ##  $ levels_given: logi FALSE
    ##  $ predictors  : logi FALSE

The ceiling is `floor(n / min_cell_n)`, with `n` the rows that carry a
response when one is given, or the number of distinct locations when
that is smaller (one short of it when no location repeats);
`ceiling_from` says which. Past it the average cell holds too few points
to estimate anything from, or k-means has more centres to place than
distinct points. The floor is `ceiling(area / range^2)`, from the fitted
autocorrelation range: cells wider than the range average over more than
one patch of the field, mixing values the field itself keeps apart. The
range is an estimate, so by default the floor is reported and not
applied: it is a level of the ladder,
`select_resolution(criterion = "reliability")` reads from it up, and
`range_floor = TRUE` starts the ladder there.

When the floor exceeds the ceiling, the data cannot support a
tessellation that respects their own correlation structure. `supported`
is `FALSE`, a warning says so, and the ladder still runs from 2 to the
ceiling so the cost of each level stays visible. The usual causes are
too few points for the extent, or a correlation range shorter than the
spacing between observations.

## Setting a size

Of the criteria above, all but `loo` score the plain mean of a cell.
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
reports beside it, by default, a *pooled* value: the cell’s mean pulled
toward the mean of the 8 cells around it, by a weight that grows with
the number of observations the cell holds (`estimate = "pooled"`;
[`vignette("reporting")`](https://elkronos.github.io/gis_modeling_toolkit/articles/reporting.md)
shows how to read and report one). For that value no criterion sets the
count. A *size* does: the number of observations a cell is to hold on
average, `obs_per_cell`. A size `s` on `n` observations is
`max(2, round(n / s))` cells.

The default profile is one size, 4.5, scored for the pooled value, and
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
reads the count off it by that size. This is the call with nothing named
but the response:

``` r

prof_default <- resolution_profile(pts, response_var = "z")
select_resolution(prof_default)
```

    ## Resolution by size: 111 cells
    ##   size        : 4.5 rows per cell asked (obs_per_cell); 500 rows / 111 cells = 4.50
    ##   cells       : density-weighted k-means++ starts (power 1, 8 neighbours)
    ##   estimate    : pooled toward the 8 nearest cells
    ##   at this size: 0.4% of the rows are alone in their cell; median own weight 0.54
    ##   loo         : 1.54 (one cell for the layer: 2.08)
    ##   note        : the count is set by obs_per_cell = 4.5, the caller's choice;
    ##                 no criterion chose it.

`obs_per_cell` puts the profile’s ladder in those terms, one level per
size, with the centres fitted to every row. `estimate` names the cell
value the profile scores, the pooled value unless it says otherwise:

``` r

prof_size <- resolution_profile(pts, response_var = "z",
                                obs_per_cell = c(12, 9, 6, 4.5, 3))
as.data.frame(prof_size)[, c("levels", "scored_per_cell", "loo",
                             "share_alone", "own_weight")]
```

    ##   levels scored_per_cell      loo share_alone own_weight
    ## 1     42       11.904762 1.566779       0.000  0.8143550
    ## 2     56        8.928571 1.525325       0.000  0.7413916
    ## 3     83        6.024096 1.461817       0.002  0.7197274
    ## 4    111        4.504505 1.423166       0.010  0.6473814
    ## 5    167        2.994012 1.429188       0.048  0.5621471

Three of the columns say what a size leaves a cell value to stand on.

- `loo` is the leave-one-out error of the pooled value: each observation
  is predicted by the value of its cell computed without it, and the
  squared errors are averaged. It needs no variogram. `sac = "auto"`,
  the default, fits none on a ladder of sizes, so `cp` and `reliability`
  are `NA` in this profile.
- `share_alone` is the share of the observations that are the only one
  in their cell.
- `own_weight` is the median weight a pooled value gives its cell’s own
  mean. At 1 the pooled values are the plain means; at 0 each is the
  mean of the cell’s neighbours.

At the size of the default profile this ladder has other centres, and so
another `loo`, than that profile had: a profile draws the starts of its
levels in turn from one random stream under `seed`, so a level’s centres
depend on the other levels of the ladder.

`select_resolution(criterion = "size")` returns the level of a size. It
searches for no optimum, so the result has no flat region and sits on no
bound:

``` r

sel_size <- select_resolution(prof_size, criterion = "size",
                              obs_per_cell = 4.5)
sel_size
```

    ## Resolution by size: 111 cells
    ##   size        : 4.5 rows per cell asked (obs_per_cell); 500 rows / 111 cells = 4.50
    ##   cells       : density-weighted k-means++ starts (power 1, 8 neighbours)
    ##   estimate    : pooled toward the 8 nearest cells
    ##   at this size: 1.0% of the rows are alone in their cell; median own weight 0.65
    ##   loo         : 1.42 (one cell for the layer: 2.08)
    ##   note        : the count is set by obs_per_cell = 4.5, the caller's choice;
    ##                 no criterion chose it.

`loo` is in the table to compare sizes by, not to choose one. For pooled
values it tends to keep falling as the cells shrink, so its minimum
tends to be the finest size offered, and the end of the ladder would
then be choosing. On 200 simulated layers of 1,500 clustered
observations of a field plus noise, profiled at seven sizes from 18 to
2, the lowest `loo` of the pooled value was at the size 2 on 176
(`dev/sim-pooled-values.R` in the package’s repository;
[`?resolution_profile`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
gives the settings).

In the table above the lowest `loo` is at 111 cells, which is not the
finest size offered: a tendency, not a rule.

The size is therefore the caller’s choice, made on other grounds: how
many observations a cell value should stand on, which `share_alone` and
`own_weight` describe. Under a size neither the cells, which are fitted
to the coordinates, nor their number, which is a count of rows, depend
on the response, so the count is not fitted to the values it will
summarise.

A size does not ask whether the response has spatial structure the cells
can use, so the selection asks it separately. `structured` is `FALSE`
when one mean for the whole layer predicts a left-out observation as
well as the cells do. The count stands either way, since the caller set
it, and the functions that take the selection as a count then warn
(class `spatialkit_no_structure`) that the cell values will say little
that one mean does not.

``` r

data.frame(structured = sel_size$structured, loo = sel_size$loo,
           one_mean = sel_size$loo_one_cell)
```

    ##   structured      loo one_mean
    ## 1       TRUE 1.423166 2.081281

The selection goes to
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
as in “Next” below. It carries a record of how the count was set,
`sel_size$count_rule`, which the seeds, the cells and the assigned layer
hand on to
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md).
The summary then checks that it reports the cell value the count was
chosen for, on the cells and rows that were scored.

## Leave-one-out for plain means

With plain means the leave-one-out error does choose a count. Each
observation is predicted by the mean of the other observations of its
cell, or, when it is alone in its cell, by the mean of the cell of its
nearest other observation. `estimate = "mean"` scores that error in the
profile, and the default `obs_per_cell = "auto"` then offers the seven
sizes 18, 12, 9, 6, 4.5, 3 and 2:

``` r

prof_loo <- resolution_profile(pts, response_var = "z", estimate = "mean")
as.data.frame(prof_loo)[, c("levels", "scored_per_cell", "loo", "loo_se",
                            "share_alone")]
```

    ##   levels scored_per_cell      loo    loo_se share_alone
    ## 1     28       17.857143 1.640425 0.1004998       0.000
    ## 2     42       11.904762 1.601586 0.1024429       0.000
    ## 3     56        8.928571 1.632183 0.1013714       0.000
    ## 4     83        6.024096 1.645847 0.1004743       0.000
    ## 5    111        4.504505 1.664467 0.1014937       0.004
    ## 6    167        2.994012 1.600161 0.1032003       0.050
    ## 7    250        2.000000 2.008358 0.1222162       0.174

``` r

sel_loo <- select_resolution(prof_loo, criterion = "loo")
sel_loo
```

    ## Resolution by loo: 167 cells
    ##   flat region : 42, 167 (2 of 7 levels)
    ##   cells       : density-weighted k-means++ starts (power 1, 8 neighbours)
    ##   estimate    : plain means
    ##   at this size: 5.0% of the rows are alone in their cell
    ##   loo         : 1.6 (one cell for the layer: 2.08)

The lowest error is at 167 cells, inside the sizes offered, so the
criterion chose it.

`sel_loo$flat` is every level whose error is within `tol` of the lowest,
and it is the answer to quote, as for any criterion. `loo_se` is the
standard error of one level’s mean squared error. It is not the standard
error of the difference between two levels: their errors come from the
same observations.

Three things to know before relying on this criterion. The count now
depends on the values of the response, so the next section applies, and
`select_on = "split"` reads `loo` on one half of the layer. The error
describes locations among the observations: it compares counts where
data are near, and says nothing about the error of a map away from them.
And it sets the count of plain means only; for pooled values a size
does, as above. In the simulation of that section the lowest `loo` of
the plain mean was at the finest size on none of the 200 layers: it was
at 6 observations per cell on 86 of them, at 4.5 on 51 and at 9 on 49.

`criterion = "auto"` is the rule that goes with the profile’s estimate:
`"size"` for pooled values and `"loo"` for plain means. It is
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)‘s
default, and a profile built with `estimate` and passed straight to
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
as `n` is read the same way. `criterion = "cp"` reads Mallows’ Cp, which
neither profile of these two sections has, since neither fitted a
variogram.

## Choosing on one half, estimating on the other

Picking the level from the same response you then aggregate is a form of
double-dipping: the cell count was chosen to suit one realisation of the
field. `select_on = "split"` profiles part of the layer and hands back
the row positions of both parts.

``` r

prof_split <- resolution_profile(pts, response_var = "z", n_levels = 16,
                                 select_on = "split")
sp <- attr(prof_split, "split")
sp
```

    ## Spatial half-split (block_kfold, seed 123): 251 selection rows, 249 estimation rows
    ##   $selection and $estimation are row positions in the layer as passed.

The split is spatially blocked, not random, and records the method and
seed that produced it. A random half would put neighbours of every
selection point in the estimation set, and on an autocorrelated field
that leaks the structure you are choosing against. Blocking reduces that
leak without removing it: the two halves share a border, and points on
either side of it within the correlation range are still correlated, so
read the separation as a large improvement on a random half rather than
as independence.

Choose from `prof_split`, then aggregate the rows in `sp$estimation`.
The cells on the ladder are still drawn on every point, so the count it
picks is a count for the whole layer, and the ladder is as long as the
full profile’s. Only the steps that read the response (the variogram,
`cp`, `moran_z`, `loo`) use the selection half. The cost is power: those
criteria see half the data, so the profile is noisier and the band
wider. The estimation rows also fill only their own half of the layer: a
cell inside the selection half gets none of them, and a cell across the
border between the halves is estimated from the part of it on the
estimation side. Read standard errors only off cells whose points are
all estimation rows.

``` r

unlist(attr(prof_split, "bounds")[c("floor", "ceiling", "n", "supported")])
```

    ##     floor   ceiling         n supported 
    ##        21        55       500         1

The floor can differ from the full profile’s, because it comes from the
range estimated on the selection half.

## Next

Hand the selection to whichever call decides the cell count. For a
Voronoi tessellation that is
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md).
Given the
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
result itself, it returns the centres of the partition the profile
scored, and a k-means partition, weighted or not, is the Voronoi
partition of its centres, so the cells built on those seeds are the
cells the criteria judged; the bare count `sel$best` would run a new
k-means, whose cells are not. This selection sits on the range floor, so
the call warns that a bound chose the count (class
`spatialkit_bound_chosen`; warnings are not shown in this article).
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
then works on the seeds:

``` r

bnd   <- clip_target_for(pts, expand = 0.02, quiet = TRUE)
seeds <- get_voronoi_seeds(boundary = bnd, method = "kmeans", n = sel,
                           sample_points = pts)
cells <- build_tessellation(seeds, boundary = bnd, method = "voronoi",
                            quiet = TRUE)
nrow(cells$cells)
```

    ## [1] 17

For a lattice the count goes to
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
itself, and the points go in directly:

``` r

hex <- build_tessellation(pts, boundary = bnd, method = "hex",
                          approx_n_cells = sel$best, quiet = TRUE)
nrow(hex$cells)
```

    ## [1] 27

The hex count overshoots the request: the target is adjusted for packing
density and then clipped to an irregular boundary, so `approx_n_cells`
is approximate in both directions. The seeded Voronoi hits the number
exactly, because the seeds are the cells.

Passing `approx_n_cells` with `method = "voronoi"` warns and is ignored.
Without the seeding step you get one cell per observation, which is a
nearest-neighbour interpolation rather than a resolution.

[`vignette("spatialkit_nc_demo")`](https://elkronos.github.io/gis_modeling_toolkit/articles/spatialkit_nc_demo.md)
runs the whole pipeline on real boundaries.
[`?resolution_profile`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
documents how each criterion behaved on simulated fields, with the
references;
[`?determine_optimal_levels`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
covers the k-means++ restart budget and the nine-cell floor below which
the model-aware criteria carry no information.
