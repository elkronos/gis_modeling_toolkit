# dev/sim-pooled-values.R
# ---------------------------------------------------------------------------
# The simulation behind the figures the help, the vignettes and NEWS quote
# for the pooled cell values, the interval for one observation and the
# leave-one-out error (summarize_by_cell(estimate = "pooled", obs_level = ),
# cell_values_at(), resolution_profile(estimate = , obs_per_cell = )).
#
# It simulates 200 layers, takes each through the package's exported calls,
# and writes six outputs, a to f, to dev/sim-pooled-values.csv.  Every figure
# those pages quote, and every statement they make about how the new numbers
# behave, is a row of that file.  A statement the file does not bear out is
# not in the pages.
#
# Usage, from the root of the working tree:
#   Rscript --vanilla dev/sim-pooled-values.R
# The run takes about four minutes, needs the FNN package (the script
# computes nothing with it; the package's searches use it) and prints one
# line per layer.  A trial of fewer layers has to write elsewhere:
#   Rscript --vanilla dev/sim-pooled-values.R <layers> <csv file> [<rds file>]
# <rds file> keeps the per-row and per-layer tables the outputs are made of.
#
# Beside the CSV a run writes its record, dev/sim-pooled-values-run.dcf (for
# a trial, <csv file> with "-run.dcf" for ".csv"): the MD5 sums of this
# script and of the CSV, each read as lines of text, the layers, the date,
# the seconds and the versions of R, sf, FNN and the package.
# tests/testthat/test-core-docs.R compares both sums with the files, so the
# figures the pages quote are those of this script as it stands: after an
# edit here, or of the CSV, that test fails until the script is run again.
#
# SETTINGS.  Everything about layer f (1 to 200) is drawn under
# set.seed(20261012 + f), R's default generators, in the order of this list.
#
#  1. The square is 1000 by 1000.  The parents are 25 points: 25 x
#     coordinates, then 25 y coordinates, uniform on it.
#  2. A clustered draw of n locations.  Candidates are drawn 2n at a time and
#     those inside the square kept, until n are in hand; the first n are
#     used.  A candidate is, with probability 0.8, a parent drawn at random
#     with replacement plus a normal displacement of standard deviation 40 in
#     each coordinate, and otherwise uniform on the square.  One batch draws,
#     in this order: 2n uniforms (below 0.8: at a parent), 2n parents
#     (sample.int), 2n x displacements, 2n y displacements, 2n uniform x, 2n
#     uniform y.  The rows are a clustered draw of 2,100, and the further
#     locations a second clustered draw of 20,000 round the same parents.
#  3. The field is a sum of 1,000 cosine waves,
#       Z(s) = sqrt(2 / 1000) * sum(cos(s . om_j + ph_j)),
#     with the frequencies
#       om <- cbind(rnorm(1000), rnorm(1000)) / (100 * abs(rnorm(1000)))
#     (a bivariate Cauchy distribution of scale 1/100) and the phases
#       ph <- runif(1000, 0, 2 * pi).
#     Its covariance is exp(-d / 100): unit sill, range parameter 100, an
#     effective range of about 300.  It is evaluated at the rows and at the
#     further locations; that is the noise-free field.  The form gives the
#     same field at any number of locations in time proportional to their
#     number.
#  4. The response at a row is the field plus normal noise of variance 0.3.
#  5. The block: a centre uniform on the middle square, 250 to 750 in each
#     coordinate (x, then y); the 300 rows nearest it are held out, and the
#     block's radius is the distance of the furthest of them.  Held out among
#     the others: 300 of the other 1,800 rows, by sample().  The other 1,500
#     rows are the layer.  Last, 1,500 standard normal values: the white-noise
#     response of output d, on the layer's rows.
#  6. The package's calls, on the layer:
#       resolution_profile(, "y", estimate = "pooled", obs_per_cell = 4.5,
#                          partition = "density", nstart = 3, sac = "auto",
#                          seed = 123)
#       select_resolution(criterion = "size")
#       get_voronoi_seeds(method = "kmeans", n = <selection>, sample_points = )
#       build_tessellation(, method = "voronoi")
#       assign_features_to_polygons(, polygon_id_col = "cell_id")
#       summarize_by_cell(, estimate = "pooled", obs_level = 0.90,
#                         conf_level = 0.90)
#       cell_values_at(, seeds = )   for the 600 held-out rows and for the
#                                    20,000 further locations
#     For output d the same calls, without cell_values_at(), on the
#     white-noise response, which gives the same cells.  For output f a
#     profile at the seven sizes, obs_per_cell = c(18, 12, 9, 6, 4.5, 3, 2),
#     for each estimate, with the other settings of the profile above.
#
# The layer is an sf layer in EPSG:32632 with the square's corner at
# (500000, 5000000), so that its coordinates are metres in a projected CRS
# and no call transforms them.  The field is evaluated on the square's own
# coordinates, 0 to 1000.
#
# OUTPUTS, one row of the file per figure.  `value` is taken over the rows
# (or cells, or layers) of all layers together; `n` is how many; `layer_min`
# and `layer_max` are the least and the greatest of the 200 single-layer
# values, where a single layer has one.
#
#  a  The root mean squared error of the pooled value over that of the plain
#     mean, at the rows held out among the others: in all, and by the number
#     of the layer's rows in the cell the value came from (1, 2 to 4, 5 to 9,
#     10 or more).  The same in the block.
#  b  Coverage of the interval for one observation (..obs_lo_, ..obs_hi_) at
#     0.90: at the rows held out among the others, in the block, and over
#     both by dist_ratio (below 1, 1 to 2, 2 to 4, 4 to 8, above 8; a class
#     holds its upper end).
#  c  Coverage of the t interval of the plain mean (..ci_lo_, ..ci_hi_) at
#     0.90, read as a range for one observation, over the held-out rows whose
#     cell has that interval: among the others, and in the block.
#  d  On the white-noise responses: the median own weight (the median over
#     layers of a layer's median); the share of layers with tau2 = 0; the
#     share on which spatialkit_no_structure is raised, by the summary and by
#     the seeding; the correlation between a cell's pooled value and the mean
#     of the pooled values of its 8 nearest cells, beside the same for plain
#     means (the mean over layers of a layer's correlation).
#  e  The realised squared error of the pooled value over (1 - w) tau2,
#     n-weighted over cells.  Among the rows: a cell's level is the mean of
#     the noise-free field over the further locations outside the block that
#     take the cell's value (cells with at least 10).  In the block: the mean
#     over the further locations of the cell that lie within the block's
#     radius of its centre (cells with at least 10 of those).
#  f  Over the profiles at the seven sizes: the share of layers on which
#     `loo` is lowest at the finest size, for the pooled value and for the
#     plain mean (of equal errors, the coarser level), and the share at each
#     size.
#
# A held-out row or a further location takes the value of the cell
# cell_values_at() gives it: its own cell, or, when that cell has no value,
# the cell the value was filled from.
#
# `dev/` is in .Rbuildignore, so nothing here ships.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(sf))

# Load the working tree, never the installed package.
if (!file.exists("DESCRIPTION") || !requireNamespace("pkgload", quietly = TRUE))
  stop("Run from the root of the working tree, with pkgload installed.")
pkgload::load_all(".", quiet = TRUE, helpers = FALSE, export_all = FALSE)
if (!requireNamespace("FNN", quietly = TRUE))
  stop("The simulation needs the FNN package: without it the searches of the ",
       "20,000 further locations compare every pair.")
spatialkit_quiet(TRUE)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# ---- Settings ----------------------------------------------------------------------

LAYERS    <- 200L
SEED      <- 20261012L
SIDE      <- 1000
PARENTS   <- 25L
NEAR      <- 0.8
SPREAD    <- 40
N_ROWS    <- 2100L
N_FURTHER <- 20000L
WAVES     <- 1000L
RANGE     <- 100
NOISE_VAR <- 0.3
N_BLOCK   <- 300L
N_AMONG   <- 300L
SIZE      <- 4.5
SIZES     <- c(18, 12, 9, 6, 4.5, 3, 2)
LEVEL     <- 0.90
FIT_SEED  <- 123
MIN_FURTHER <- 10L
ORIGIN    <- c(5e5, 5e6)
CRS       <- 32632

OUT_DEFAULT <- file.path("dev", "sim-pooled-values.csv")
args   <- commandArgs(trailingOnly = TRUE)
layers <- if (length(args) >= 1L) suppressWarnings(as.numeric(args[[1L]])) else LAYERS
out_csv <- if (length(args) >= 2L) args[[2L]] else OUT_DEFAULT
out_rds <- if (length(args) >= 3L) args[[3L]] else NULL
if (is.na(layers) || layers != round(layers) || layers < 1 || layers > LAYERS)
  stop("<layers> is a whole number from 1 to ", LAYERS, ".")
layers <- as.integer(layers)
if (layers != LAYERS && normalizePath(out_csv, mustWork = FALSE) ==
    normalizePath(OUT_DEFAULT, mustWork = FALSE))
  stop("A trial of ", layers, " layers has to name another file than ", OUT_DEFAULT, ".")
out_run <- paste0(sub("\\.csv$", "", out_csv), "-run.dcf")

# This file, as the run was started on it, and the MD5 sum of a file read as
# lines of text, which is the same whatever the line endings of a checkout.
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(this_file) != 1L || !file.exists(this_file))
  this_file <- file.path("dev", "sim-pooled-values.R")
text_md5 <- function(file) {
  tmp <- tempfile()
  on.exit(unlink(tmp))
  con <- file(tmp, "wb")
  writeLines(readLines(file, warn = FALSE), con, sep = "\n", useBytes = TRUE)
  close(con)
  unname(tools::md5sum(tmp))
}
script_md5 <- text_md5(this_file)

# ---- The draws ---------------------------------------------------------------------

clustered_draw <- function(n, parents) {
  out <- matrix(numeric(0), 0L, 2L)
  while (nrow(out) < n) {
    m    <- 2L * n
    near <- stats::runif(m) < NEAR
    p    <- sample.int(nrow(parents), m, replace = TRUE)
    dx   <- stats::rnorm(m, 0, SPREAD)
    dy   <- stats::rnorm(m, 0, SPREAD)
    ux   <- stats::runif(m, 0, SIDE)
    uy   <- stats::runif(m, 0, SIDE)
    x    <- ifelse(near, parents[p, 1L] + dx, ux)
    y    <- ifelse(near, parents[p, 2L] + dy, uy)
    keep <- x >= 0 & x <= SIDE & y >= 0 & y <= SIDE
    out  <- rbind(out, cbind(x, y)[keep, , drop = FALSE])
  }
  unname(out[seq_len(n), , drop = FALSE])
}

wave_field <- function(xy, om, ph) {
  out <- numeric(nrow(xy))
  for (from in seq.int(1L, nrow(xy), by = 4000L)) {
    at <- from:min(from + 3999L, nrow(xy))
    out[at] <- sqrt(2 / WAVES) *
      rowSums(cos(xy[at, , drop = FALSE] %*% t(om) + rep(ph, each = length(at))))
  }
  out
}

draw_layer <- function(f) {
  set.seed(SEED + f)
  parents <- cbind(stats::runif(PARENTS, 0, SIDE), stats::runif(PARENTS, 0, SIDE))
  rows    <- clustered_draw(N_ROWS, parents)
  further <- clustered_draw(N_FURTHER, parents)
  om <- cbind(stats::rnorm(WAVES), stats::rnorm(WAVES)) / (RANGE * abs(stats::rnorm(WAVES)))
  ph <- stats::runif(WAVES, 0, 2 * pi)
  z_rows    <- wave_field(rows, om, ph)
  z_further <- wave_field(further, om, ph)
  y   <- z_rows + stats::rnorm(N_ROWS, 0, sqrt(NOISE_VAR))
  ctr <- stats::runif(2L, SIDE / 4, 3 * SIDE / 4)
  d   <- sqrt((rows[, 1L] - ctr[1L])^2 + (rows[, 2L] - ctr[2L])^2)
  block  <- rank(d, ties.method = "first") <= N_BLOCK
  radius <- max(d[block])
  among  <- seq_len(N_ROWS) %in% sample(which(!block), N_AMONG)
  layer  <- !block & !among
  white  <- stats::rnorm(sum(layer))
  list(rows = rows, further = further, y = y, z_further = z_further, ctr = ctr,
       radius = radius, block = block, among = among, layer = layer, white = white)
}

as_layer <- function(xy, ...) {
  d <- data.frame(x = ORIGIN[1L] + xy[, 1L], y_coord = ORIGIN[2L] + xy[, 2L], ...)
  sf::st_as_sf(d, coords = c("x", "y_coord"), crs = CRS)
}

# ---- The package's calls -----------------------------------------------------------

# Runs one call, keeps the class of every warning it raises and lets none
# through; messages are not shown.
seen <- new.env()
watch <- function(step, expr) {
  withCallingHandlers(
    suppressMessages(expr),
    warning = function(w) {
      cls <- class(w)[1L]
      key <- paste(step, cls, sep = ": ")
      seen[[key]] <- c(seen[[key]], conditionMessage(w))
      step_classes <<- c(step_classes, key)
      invokeRestart("muffleWarning")
    })
}
step_classes <- character(0)

cells_of <- function(lay) {
  step_classes <<- character(0)
  prof <- watch("resolution_profile", resolution_profile(
    lay, "y", estimate = "pooled", obs_per_cell = SIZE, partition = "density",
    nstart = 3, sac = "auto", seed = FIT_SEED))
  sel   <- watch("select_resolution", select_resolution(prof, criterion = "size"))
  seeds <- watch("get_voronoi_seeds",
                 get_voronoi_seeds(method = "kmeans", n = sel, sample_points = lay))
  tess  <- watch("build_tessellation", build_tessellation(seeds, method = "voronoi"))
  asg   <- watch("assign_features_to_polygons",
                 assign_features_to_polygons(lay, tess$cells, polygon_id_col = "cell_id"))
  out   <- watch("summarize_by_cell", summarize_by_cell(
    asg, "y", id_col = "cell_id", cells_sf = tess$cells, estimate = "pooled",
    obs_level = LEVEL, conf_level = LEVEL))
  seeds$cell_id <- tess$index
  list(sel = sel, seeds = seeds, tess = tess, asg = asg, out = out,
       raised = step_classes)
}

values_at <- function(new, fit) {
  at <- watch("cell_values_at",
              cell_values_at(new, fit$out, fit$asg, "y", id_col = "cell_id",
                             seeds = fit$seeds))
  # The row of the summary the value came from.
  at$from <- match(ifelse(at$filled, at$filled_from, at$cell_id), fit$out$cell_id)
  at
}

# The 8 nearest other cells of every cell, by the centroids of their rows.
neighbour_means <- function(fit) {
  xy  <- sf::st_coordinates(fit$asg)
  row <- match(fit$asg$cell_id, fit$out$cell_id)
  n   <- tabulate(row, nrow(fit$out))
  has <- which(n > 0L)
  cen <- rowsum(xy, row) / n[has]
  d2  <- as.matrix(stats::dist(cen))^2
  diag(d2) <- Inf
  k   <- min(8L, length(has) - 1L)
  nb  <- t(apply(d2, 1L, function(d) order(d)[seq_len(k)]))
  list(has = has, nb = nb, n = n[has])
}

# ---- The layers --------------------------------------------------------------------

size_class <- function(n) cut(n, c(0, 1, 4, 9, Inf),
                              labels = c("1 row", "2 to 4 rows", "5 to 9 rows",
                                         "10 rows or more"))
ratio_class <- function(r) cut(r, c(-Inf, 1, 2, 4, 8, Inf),
                               labels = c("below 1", "1 to 2", "2 to 4", "4 to 8", "above 8"))

held <- level <- white <- sizes <- info <- vector("list", layers)
t0 <- proc.time()[[3L]]
for (f in seq_len(layers)) {
  L   <- draw_layer(f)
  lay <- as_layer(L$rows[L$layer, , drop = FALSE], y = L$y[L$layer])
  fit <- cells_of(lay)
  out <- fit$out
  rec <- attr(out, "pooled")$columns$resp_y

  # a, b, c: the held-out rows
  hsel <- L$among | L$block
  at   <- values_at(as_layer(L$rows[hsel, , drop = FALSE]), fit)
  yh   <- L$y[hsel]
  stopifnot(isTRUE(all.equal(at$value, out$..pooled_resp_y[at$from])))
  held[[f]] <- data.frame(
    layer   = f,
    where   = ifelse(L$block[hsel], "block", "among"),
    n_cell  = out$n[at$from],
    e_pooled = yh - at$value,
    e_mean   = yh - out$resp_mean_y[at$from],
    in_obs  = yh >= at$obs_lo & yh <= at$obs_hi,
    in_ci   = yh >= out$..ci_lo_resp_y[at$from] & yh <= out$..ci_hi_resp_y[at$from],
    ratio   = at$dist_ratio,
    filled  = at$filled)

  # e: the level of a cell from the noise-free field at the further locations
  atf <- values_at(as_layer(L$further), fit)
  inb <- sqrt((L$further[, 1L] - L$ctr[1L])^2 + (L$further[, 2L] - L$ctr[2L])^2) <= L$radius
  g1  <- (1 - out$..pooled_w_resp_y) * rec$tau2
  one <- function(sel, where) {
    sel <- sel & !is.na(atf$from)
    cnt <- tabulate(atf$from[sel], nrow(out))
    ok  <- cnt >= MIN_FURTHER & !is.na(out$..pooled_resp_y)
    sums <- rowsum(L$z_further[sel], atf$from[sel])
    lev  <- rep(NA_real_, nrow(out))
    lev[as.integer(rownames(sums))] <- sums[, 1L] / cnt[as.integer(rownames(sums))]
    data.frame(layer = f, where = where, cells = sum(ok),
               realised = sum(out$n[ok] * (out$..pooled_resp_y[ok] - lev[ok])^2),
               stated   = sum(out$n[ok] * g1[ok]))
  }
  level[[f]] <- rbind(one(!inb, "among"), one(inb, "block"))
  # Every warning of the field's calls: the six of item 6 and cell_values_at().
  raised_field <- step_classes

  # d: white noise on the same rows, which gives the same cells
  lay_w <- lay
  lay_w$y <- L$white
  fw   <- cells_of(lay_w)
  stopifnot(identical(sf::st_coordinates(fw$seeds), sf::st_coordinates(fit$seeds)))
  ow   <- fw$out
  recw <- attr(ow, "pooled")$columns$resp_y
  nbw  <- neighbour_means(fw)
  with_nb <- function(v) stats::cor(v, rowMeans(matrix(v[nbw$nb], nrow(nbw$nb))))
  # The neighbours found here are the summary's: its target is their mean.
  mw <- ow$resp_mean_y[nbw$has]
  tg <- rowSums(matrix((nbw$n * mw)[nbw$nb], nrow(nbw$nb))) /
        rowSums(matrix(nbw$n[nbw$nb], nrow(nbw$nb)))
  stopifnot(isTRUE(all.equal(tg, ow$..pooled_target_resp_y[nbw$has])))
  white[[f]] <- data.frame(
    layer = f, own_weight = recw$own_weight, tau2_zero = recw$tau2 == 0,
    structured = recw$structured,
    raised_summary = "summarize_by_cell: spatialkit_no_structure" %in% fw$raised,
    raised_seeding = "get_voronoi_seeds: spatialkit_no_structure" %in% fw$raised,
    cor_pooled = with_nb(ow$..pooled_resp_y[nbw$has]), cor_mean = with_nb(mw))

  # f: the seven sizes, for each estimate
  lowest <- function(estimate) {
    p <- watch("resolution_profile (seven sizes)", resolution_profile(
      lay, "y", estimate = estimate, obs_per_cell = SIZES, partition = "density",
      nstart = 3, sac = "auto", seed = FIT_SEED))
    lev <- attr(p, "bounds")$size_levels
    p   <- p[order(p$levels), ]
    at  <- p$levels[which.min(p$loo)]
    as.numeric(names(lev)[match(at, lev)])
  }
  step_classes <- character(0)
  sizes[[f]] <- data.frame(layer = f, pooled = lowest("pooled"), mean = lowest("mean"))
  raised_sizes <- step_classes

  info[[f]] <- data.frame(
    layer = f, cells = nrow(fit$tess$cells), occupied = rec$cells, rows = rec$n,
    share_alone = rec$share_alone, own_weight = rec$own_weight, tau2 = rec$tau2,
    structured = rec$structured, check = attr(out, "count_rule")$check$cells,
    filled_held = sum(at$filled), filled_further = sum(atf$filled),
    warnings = paste(unique(raised_field), collapse = "; "),
    warnings_white = paste(unique(fw$raised), collapse = "; "),
    warnings_sizes = paste(unique(raised_sizes), collapse = "; "))
  # The line names every warning of the layer but the one output d counts,
  # spatialkit_no_structure on the white-noise response.
  other_white <- unique(grep("spatialkit_no_structure$", fw$raised, value = TRUE, invert = TRUE))
  said <- c(if (nzchar(info[[f]]$warnings)) paste0("raised: ", info[[f]]$warnings),
            if (length(other_white))
              paste0("on white noise: ", paste(other_white, collapse = "; ")),
            if (nzchar(info[[f]]$warnings_sizes))
              paste0("at the seven sizes: ", info[[f]]$warnings_sizes))
  cat(sprintf("layer %3d: %d cells, %d occupied, %.1f%% of rows alone, check %s, %d s%s\n",
              f, info[[f]]$cells, info[[f]]$occupied, 100 * info[[f]]$share_alone,
              info[[f]]$check, round(proc.time()[[3L]] - t0),
              if (length(said)) paste0("; ", said, collapse = "") else ""))
}
held  <- do.call(rbind, held)
level <- do.call(rbind, level)
white <- do.call(rbind, white)
sizes <- do.call(rbind, sizes)
info  <- do.call(rbind, info)

# ---- The six outputs ---------------------------------------------------------------

rows <- list()
add <- function(output, quantity, where, class, value, n, by_layer = NULL) {
  rows[[length(rows) + 1L]] <<- data.frame(
    output = output, quantity = quantity, where = where, class = class,
    value = value, n = n,
    layer_min = if (length(by_layer)) min(by_layer, na.rm = TRUE) else NA_real_,
    layer_max = if (length(by_layer)) max(by_layer, na.rm = TRUE) else NA_real_)
}
rmse_ratio <- function(x) if (nrow(x)) sqrt(sum(x$e_pooled^2) / sum(x$e_mean^2)) else NA_real_
per_layer  <- function(x, fun) vapply(split(x, x$layer), fun, numeric(1))

for (w in c("among", "block")) {
  x <- held[held$where == w, ]
  add("a", "RMSE of the pooled value over RMSE of the plain mean", w, "all rows",
      rmse_ratio(x), nrow(x), per_layer(x, rmse_ratio))
  for (cl in levels(size_class(1))) {
    xc <- x[size_class(x$n_cell) == cl, ]
    add("a", "RMSE of the pooled value over RMSE of the plain mean", w,
        paste(cl, "of the layer in the cell"), rmse_ratio(xc), nrow(xc))
  }
}
for (w in c("among", "block")) {
  x <- held[held$where == w, ]
  add("b", "coverage of the 0.90 interval for one observation", w, "all rows",
      mean(x$in_obs), nrow(x), per_layer(x, function(v) mean(v$in_obs)))
}
for (cl in levels(ratio_class(1))) {
  x <- held[ratio_class(held$ratio) == cl, ]
  add("b", "coverage of the 0.90 interval for one observation", "among and block",
      paste("dist_ratio", cl), mean(x$in_obs), nrow(x))
}
for (w in c("among", "block")) {
  x <- held[held$where == w & !is.na(held$in_ci), ]
  add("c", "coverage of the 0.90 t interval of the plain mean, read as a range for one observation",
      w, "rows whose cell has the interval", mean(x$in_ci), nrow(x),
      per_layer(x, function(v) mean(v$in_ci)))
}
add("d", "median own weight", "white noise", "layers",
    stats::median(white$own_weight), nrow(white), white$own_weight)
add("d", "share of layers with tau2 = 0", "white noise", "layers",
    mean(white$tau2_zero), nrow(white))
add("d", "share of layers on which summarize_by_cell() raises spatialkit_no_structure",
    "white noise", "layers", mean(white$raised_summary), nrow(white))
add("d", "share of layers on which get_voronoi_seeds() raises spatialkit_no_structure",
    "white noise", "layers", mean(white$raised_seeding), nrow(white))
add("d", "correlation of a cell's pooled value with the mean of its 8 nearest cells'",
    "white noise", "layers", mean(white$cor_pooled), nrow(white), white$cor_pooled)
add("d", "correlation of a cell's plain mean with the mean of its 8 nearest cells'",
    "white noise", "layers", mean(white$cor_mean), nrow(white), white$cor_mean)
for (w in c("among", "block")) {
  x <- level[level$where == w, ]
  add("e", "realised squared error of the pooled value over (1 - w) tau2", w,
      sprintf("cells with at least %d further locations", MIN_FURTHER),
      sum(x$realised) / sum(x$stated), sum(x$cells), x$realised / x$stated)
}
for (est in c("pooled", "mean")) {
  what <- if (est == "pooled") "pooled value" else "plain mean"
  add("f", "share of layers with the lowest loo at the finest size", what,
      sprintf("obs_per_cell = %s", min(SIZES)), mean(sizes[[est]] == min(SIZES)), nrow(sizes))
  for (s in SIZES)
    add("f", "share of layers with the lowest loo at this size", what,
        sprintf("obs_per_cell = %s", s), mean(sizes[[est]] == s), nrow(sizes))
}
result <- do.call(rbind, rows)
result$value     <- signif(result$value, 6)
result$layer_min <- signif(result$layer_min, 6)
result$layer_max <- signif(result$layer_max, 6)
utils::write.csv(result, out_csv, row.names = FALSE)
if (!is.null(out_rds))
  saveRDS(list(held = held, level = level, white = white, sizes = sizes, info = info,
               result = result, warnings = as.list(seen)), out_rds)

# The record of the run.  The script is read again: a file edited while the
# run was under way is not the file the figures came from.
if (!identical(text_md5(this_file), script_md5))
  stop(this_file, " changed during the run: run it again.")
seconds <- round(proc.time()[[3L]] - t0)
write.dcf(data.frame(
  script = "dev/sim-pooled-values.R", script_md5 = script_md5,
  csv_md5 = text_md5(out_csv), layers = layers, date = format(Sys.Date()),
  seconds = seconds, R = as.character(getRversion()), platform = R.version$platform,
  sf = as.character(utils::packageVersion("sf")),
  FNN = as.character(utils::packageVersion("FNN")),
  spatialkit = unname(read.dcf("DESCRIPTION", "Version")[1L, 1L]),
  stringsAsFactors = FALSE), out_run)

# ---- What the run was --------------------------------------------------------------

cat(sprintf("\n%d layers in %.0f s; R %s, sf %s, FNN %s\n", layers, seconds,
            getRversion(), utils::packageVersion("sf"), utils::packageVersion("FNN")))
cat(sprintf("cells per layer %d to %d; occupied %d to %d; rows alone in their cell %.2f%% (mean)\n",
            min(info$cells), max(info$cells), min(info$occupied), max(info$occupied),
            100 * mean(info$share_alone)))
cat(sprintf("held-out rows filled from another cell: %d; further locations filled: %d\n",
            sum(info$filled_held), sum(info$filled_further)))
cat("the summary's check of the cells, by layer:\n")
print(table(info$check, useNA = "ifany"))
cat("structured, field layers / white-noise layers:\n")
print(table(field = info$structured, useNA = "ifany"))
print(table(white = white$structured, useNA = "ifany"))
cat("warnings raised, by call and class (times raised):\n")
for (k in sort(ls(seen))) cat(sprintf("  %s: %d\n", k, length(seen[[k]])))
cat("\n")
options(width = 250)
print(result, row.names = FALSE)
cat(sprintf("\nwritten: %s and %s\n", out_csv, out_run))
