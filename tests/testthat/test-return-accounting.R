# tests/testthat/test-return-accounting.R
# ---------------------------------------------------------------------------
# Stage 10: what was already computed is returned.
#
#   10.1 / 10.6  cv_*(): fold_status (one row per fold supplied, with why
#                a fold is missing from fold_metrics), orphan_rows,
#                n_unknown_ids and n_dropped.
#   10.4         prep_model_data(): the "dropped" attribute; $info$n_dropped
#                on every fit.
#   10.5         make_folds(block_kfold): assignment$block_id,
#                params$block_sizes, params$fold_blocks, params$blocks.
#
# Everything here is additive: no existing number moves.
# ---------------------------------------------------------------------------

ra_points <- function(n = 60, seed = 1) {
  set.seed(seed)
  d <- sf::st_as_sf(data.frame(x = runif(n, 0, 100), y = runif(n, 0, 100)),
                    coords = c("x", "y"), crs = 32632)
  d$z <- rnorm(n); d$a <- d$z + rnorm(n, sd = 0.3)
  d
}
ra_fit <- function(tr, ...) lm_spatial_fit(tr, "a", "z")

test_that("fold_status names every fold supplied and says what became of it", {
  d  <- ra_points()
  fo <- make_folds(d, k = 4, method = "random_kfold", seed = 1)
  # Clean run: every fold ok, in fold order, integer fold labels that line
  # up with fold_metrics$fold and assignment$fold.
  cv <- cv_spatial(d, "a", "z", folds = fo, fit_fn = ra_fit)
  fs <- cv$fold_status
  expect_s3_class(fs, "data.frame")
  expect_named(fs, c("fold", "status", "message"))
  expect_identical(fs$fold, 1:4)
  expect_true(is.integer(fs$fold))
  expect_identical(fs$status, rep("ok", 4L))
  expect_identical(fs$message, rep("", 4L))
  expect_identical(fs$fold, cv$fold_metrics$fold)
  expect_identical(cv$orphan_rows, integer(0))
  expect_identical(cv$n_unknown_ids, 0L)
  expect_identical(cv$n_dropped, 0L)
  # The returned folds carry no bookkeeping attributes of their own.
  expect_null(attr(cv$folds, "dropped"))
  expect_null(attr(cv$folds, "orphans"))
  expect_null(attr(cv$folds, "n_unknown_ids"))

  # A fold whose fit throws: "error" with the error text, the others
  # untouched, and the gap in fold_metrics explained.
  boom <- function(tr, ...) {
    if (!(fo$folds[[2]]$test[1] %in% tr$`..row_id`)) stop("boom on fold 2")
    ra_fit(tr)
  }
  cv2 <- suppressWarnings(cv_spatial(d, "a", "z", folds = fo, fit_fn = boom))
  expect_identical(cv2$n_folds_attempted, 4L); expect_identical(cv2$n_folds_succeeded, 3L)
  fs2 <- cv2$fold_status
  expect_identical(fs2$status, c("ok", "error", "ok", "ok"))
  expect_identical(fs2$message[2], "boom on fold 2")
  expect_identical(setdiff(fs2$fold, cv2$fold_metrics$fold), 2L)

  # A fold skipped after fitting: predict() returned the wrong length.
  badlen <- function(tr, ...) {
    f <- ra_fit(tr)
    if (!(fo$folds[[3]]$test[1] %in% tr$`..row_id`)) class(f) <- c("ra_badlen", class(f))
    f
  }
  predict.ra_badlen <- function(object, newdata, ...) {
    class(object) <- setdiff(class(object), "ra_badlen")
    predict(object, newdata, ...)[1:2]
  }
  registerS3method("predict", "ra_badlen", predict.ra_badlen)
  cv3 <- suppressWarnings(cv_spatial(d, "a", "z", folds = fo, fit_fn = badlen))
  expect_identical(cv3$fold_status$status, c("ok", "ok", "skipped", "ok"))
  expect_match(cv3$fold_status$message[3], "^predict\\(\\) returned 2 value\\(s\\) for 15 test row\\(s\\)$")
})

test_that("a fold dropped before fitting, orphan rows and unknown IDs are all reported", {
  d  <- ra_points()
  fo <- make_folds(d, k = 4, method = "random_kfold", seed = 1)
  # Every test row of fold 1 has a missing response: prep_model_data() drops
  # them, the fold has an empty test set and never reaches the fitter.
  d1 <- d; d1$a[fo$folds[[1]]$test] <- NA
  cv <- suppressWarnings(cv_spatial(d1, "a", "z", folds = fo, fit_fn = ra_fit))
  expect_identical(cv$n_folds_attempted, 4L); expect_identical(cv$n_folds_succeeded, 3L)
  fs <- cv$fold_status
  expect_identical(fs$fold, 1:4)
  expect_identical(fs$status, c("dropped", "ok", "ok", "ok"))
  expect_identical(fs$message[1], "empty test set after remapping")
  expect_identical(cv$n_dropped, length(fo$folds[[1]]$test))
  # One per row the fold set names and the data no longer has -- not one per
  # mention of it, which was k times larger.
  expect_identical(cv$n_unknown_ids, length(fo$folds[[1]]$test))
  expect_identical(cv$orphan_rows, integer(0))
  # The dropped-fold record has one schema whether or not anything was
  # dropped: ifelse() on an empty logical returns logical(0), so the empty
  # frame used to carry a logical `reason` column where the populated one
  # carries character.
  rf <- spatialkit:::.remap_folds
  ids <- seq_len(20)
  ok_folds <- lapply(1:4, function(i)
    list(train = setdiff(ids, ids[((i - 1) * 5 + 1):(i * 5)]),
         test  = ids[((i - 1) * 5 + 1):(i * 5)]))
  none <- attr(rf(ok_folds, ids, 4L, 1L), "dropped")
  one  <- ok_folds; one[[2]]$test <- integer(0)
  some <- attr(suppressWarnings(rf(one, ids, 4L, 1L)), "dropped")
  expect_identical(nrow(none), 0L)
  expect_identical(vapply(none, class, character(1)), vapply(some, class, character(1)))
  expect_type(none$reason, "character")

  # Folds built on a subset of the layer: the other rows are orphans, named.
  fo40 <- make_folds(d[1:40, ], k = 3, method = "random_kfold", seed = 1)
  expect_warning(cv2 <- cv_spatial(d, "a", "z", folds = fo40, fit_fn = ra_fit),
                 "named by no fold")
  expect_identical(cv2$orphan_rows, 41:60)
  expect_identical(cv2$n_unknown_ids, 0L)
  expect_identical(cv2$fold_status$status, rep("ok", 3L))

  # The mirror image: folds naming rows the data no longer has.
  cv3 <- cv_spatial(d[1:40, ], "a", "z", folds = fo, fit_fn = ra_fit)
  expect_gt(cv3$n_unknown_ids, 0L)
  named <- unique(unlist(lapply(fo$folds, function(f) c(f$train, f$test)),
                          use.names = FALSE))
  expect_identical(cv3$n_unknown_ids, length(setdiff(named, 1:40)))
  # Every one of those rows is named by all four folds, so the old per-entry
  # count was four times this.
  expect_identical(sum(!(unlist(lapply(fo$folds, function(f) c(f$train, f$test)),
                                use.names = FALSE) %in% 1:40)),
                   cv3$n_unknown_ids * 4L)
  # And the count does not move with k: the same five absent rows are five
  # whether they are named by 3 folds or by 10.
  ids <- seq_len(20)
  absent <- 21:25
  for (kk in c(3L, 4L, 10L)) {
    grp <- split(c(ids, absent), rep_len(seq_len(kk), length(ids) + length(absent)))
    kf <- lapply(grp, function(te)
      list(train = setdiff(c(ids, absent), te), test = te))
    rm_ <- suppressWarnings(spatialkit:::.remap_folds(kf, ids, kk, 1L))
    expect_identical(attr(rm_, "n_unknown_ids"), length(absent))
  }
  expect_identical(cv3$orphan_rows, integer(0))

  # Label-vector folds and the no-folds fallback report the same shape.
  cv4 <- cv_spatial(d, "a", "z", folds = rep(1:3, length.out = 60), fit_fn = ra_fit)
  expect_identical(cv4$fold_status$fold, 1:3)
  cv5 <- suppressMessages(cv_spatial(d, "a", "z", fit_fn = ra_fit, k = 3, seed = 1))
  expect_identical(cv5$fold_status$status, rep("ok", 3L))
  expect_identical(cv5$orphan_rows, integer(0))
})

test_that("the parallel runner reports the same fold_status as the sequential one", {
  skip_on_os("windows")
  d  <- ra_points()
  fo <- make_folds(d, k = 4, method = "random_kfold", seed = 1)
  boom <- function(tr, ...) {
    if (!(fo$folds[[2]]$test[1] %in% tr$`..row_id`)) stop("boom on fold 2")
    ra_fit(tr)
  }
  seq_cv <- suppressWarnings(cv_spatial(d, "a", "z", folds = fo, fit_fn = boom))
  par_cv <- suppressWarnings(suppressMessages(
    cv_spatial(d, "a", "z", folds = fo, fit_fn = boom, parallel = 2)))
  expect_identical(par_cv$fold_status, seq_cv$fold_status)
})

test_that("every cv_* wrapper returns the accounting fields", {
  for (fn in c("cv_gwr", "cv_bayes", "cv_spatial")) {
    src <- paste(deparse(body(get(fn, asNamespace("spatialkit")))), collapse = " ")
    for (field in c("fold_status", "orphan_rows", "n_unknown_ids", "n_dropped"))
      expect_true(grepl(field, src, fixed = TRUE), label = paste(fn, "returns", field))
  }
  skip_if_not_installed("ranger")
  d <- ra_points()
  fo <- make_folds(d, k = 3, method = "random_kfold", seed = 1)
  rf <- suppressMessages(cv_rf(d, "a", "z", folds = fo))
  expect_true(all(c("fold_status", "orphan_rows", "n_unknown_ids", "n_dropped") %in% names(rf)))
  expect_identical(rf$fold_status$status, rep("ok", 3L))
})

# ---------------------------------------------------------------------------
# 10.4
# ---------------------------------------------------------------------------

test_that("prep_model_data() records what it dropped, and why, with the rows' identities", {
  d <- sf::st_as_sf(data.frame(x = 1:8, y = 8:1,
                               resp = c(1, 2, NA, 4, 5, 6, 7, NaN),
                               pred = c(1, 2, 3, 4, Inf, 6, 7, 8)),
                    coords = c("x", "y"), crs = 32632)
  sf::st_geometry(d)[6] <- sf::st_sfc(sf::st_point(), crs = 32632)
  out <- suppressWarnings(prep_model_data(d, "resp", "pred"))
  expect_equal(nrow(out), 4L)
  dr <- attr(out, "dropped")
  expect_type(dr, "list")
  expect_named(dr, c("n", "n_geometry", "which", "row_id", "reason", "n_rows"))
  expect_identical(dr$n_rows, 4L)
  expect_identical(dr$n, 4L)
  expect_identical(dr$n_geometry, 1L)
  expect_identical(dr$which, c(3L, 5L, 6L, 8L))
  expect_null(dr$row_id)
  # Precedence: geometry, then missing (NaN is missing), then non-finite.
  expect_identical(dr$reason, c("missing", "non_finite", "geometry", "missing"))
  # Row IDs travel when the layer carries them.
  d$..row_id <- 100L + seq_len(8)
  dr2 <- attr(suppressWarnings(prep_model_data(d, "resp", "pred")), "dropped")
  expect_identical(dr2$row_id, c(103L, 105L, 106L, 108L))
  # Nothing dropped: an empty but well-formed record.
  clean <- prep_model_data(d[c(1, 2, 4, 7), ], "resp", "pred")
  dr0 <- attr(clean, "dropped")
  expect_identical(dr0$n, 0L); expect_identical(dr0$which, integer(0))
  expect_identical(dr0$reason, character(0))
})

test_that("the dropped record does not survive subsetting the layer", {
  d <- sf::st_as_sf(data.frame(x = 1:8, y = 8:1,
                               resp = c(1, 2, NA, 4, 5, 6, 7, 8),
                               pred = c(1, 2, 3, 4, Inf, 6, 7, 8)),
                    coords = c("x", "y"), crs = 32632)
  out <- suppressWarnings(prep_model_data(d, "resp", "pred"))
  # After "sf", so that vctrs (dplyr::bind_rows()) sees an sf; see
  # test-review2-assignment.R.
  expect_identical(class(out), c("sf", "spatialkit_rows", "data.frame"))
  expect_identical(attr(out, "dropped")$n, 2L)
  # Every shape of `[` leaves a plain layer with no record: the positions in
  # `which` do not survive the renumbering, and `n` is not a fact about the
  # subset.
  for (sub in list(out[1:3, ], out[out$resp > 2, ], head(out, 2),
                   out[order(out$resp), ], out[, c("resp", "pred")], out[1:3, "resp"])) {
    expect_null(attr(sub, "dropped"))
    expect_identical(class(sub), c("sf", "data.frame"))
  }
  # drop = TRUE gives a column, not a layer, and is handed back untouched.
  expect_type(out[1:3, "resp", drop = TRUE], "double")
  # Assigning a column keeps the rows, so the record stays.
  keep <- out; keep$extra <- seq_len(nrow(keep))
  expect_identical(attr(keep, "dropped")$n, 2L)
  # st_transform() rebuilds the class vector with "sf" first, so `[` reaches
  # this method through sf's own method rather than by direct dispatch; the
  # record must still go.
  tr <- sf::st_transform(out, 3857)
  expect_true(inherits(tr, "spatialkit_rows"))
  expect_null(attr(tr[1:3, ], "dropped"))
  # A record whose stamp no longer matches the layer is refused, whatever
  # path put it there.
  stale <- sf::st_transform(out, 32632)
  attr(stale, "dropped") <- list(n = 99L, n_rows = 999L)
  expect_null(spatialkit:::.get_row_record(stale, "dropped"))
  # An unstamped record -- built by hand, or made by an older version -- is
  # taken at face value rather than second-guessed.
  attr(stale, "dropped") <- list(n = 7L)
  expect_identical(spatialkit:::.get_row_record(stale, "dropped")$n, 7L)
})

test_that("a layer carrying a record still satisfies S4 dispatch written for sf", {
  # The class sat ahead of "sf" (and still does on a layer saved by an
  # earlier version), and S4 looks a class up in its own table rather than
  # walking the S3 vector, so without setOldClass() every S4 method written
  # for "sf" -- methods::as(x, "Spatial") on the way into GWmodel,
  # terra::vect(), sp's coercions -- failed a prepared layer with "no method
  # or default for coercing".
  d <- sf::st_as_sf(data.frame(x = 1:8, y = 8:1,
                               resp = c(1, 2, NA, 4, 5, 6, 7, 8),
                               pred = c(1, 2, 3, 4, Inf, 6, 7, 8)),
                    coords = c("x", "y"), crs = 32632)
  out <- suppressWarnings(prep_model_data(d, "resp", "pred"))
  expect_true(methods::existsMethod("coerce", c("spatialkit_rows", "sf")) ||
              methods::is(methods::getClass("spatialkit_rows"), "classRepresentation"))
  skip_if_not_installed("sp")
  sp_obj <- methods::as(out, "Spatial")
  expect_s4_class(sp_obj, "SpatialPointsDataFrame")
  expect_equal(nrow(sp_obj@data), nrow(out))
  expect_false(methods::is(out, "Spatial"))
  # The sf generics keep working, and the geometry survives a round trip.
  expect_equal(nrow(sf::st_coordinates(out)), nrow(out))
  expect_s3_class(sf::st_transform(out, 3857), "sf")
  expect_equal(nrow(sf::st_drop_geometry(out)), nrow(out))
})

test_that("every fit carries n_dropped, and a CV run reports it too", {
  d <- ra_points(n = 80)
  d$z[c(3, 9)] <- NA; d$a[15] <- Inf
  skip_if_not_installed("ranger")
  fit <- suppressWarnings(fit_rf_model(d, "a", "z", num_trees = 50, seed = 1))
  expect_identical(fit$info$n_dropped, 3L)
  expect_equal(fit$n, 77L)
  # Prepared by the caller: the count travels on the prepared layer's
  # attribute, and is 0 when the layer carries none.
  pre <- suppressWarnings(prep_model_data(d, "a", "z"))
  fit2 <- suppressWarnings(fit_rf_model(pre, "a", "z", num_trees = 50, seed = 1,
                                        .already_prepped = TRUE))
  expect_identical(fit2$info$n_dropped, 3L)
  # A subset of the prepared layer is not the layer the record describes, so
  # the fit reports what was dropped from *it*: nothing.
  fit2b <- suppressWarnings(fit_rf_model(pre[1:20, ], "a", "z", num_trees = 50,
                                         seed = 1, .already_prepped = TRUE))
  expect_identical(fit2b$info$n_dropped, 0L)
  expect_equal(fit2b$n, 20L)
  bare <- pre; attr(bare, "dropped") <- NULL
  fit3 <- suppressWarnings(fit_rf_model(bare, "a", "z", num_trees = 50, seed = 1,
                                        .already_prepped = TRUE))
  expect_identical(fit3$info$n_dropped, 0L)
  fo <- make_folds(d, k = 3, method = "random_kfold", seed = 1)
  cv <- suppressMessages(suppressWarnings(cv_rf(d, "a", "z", folds = fo, num_trees = 50)))
  expect_identical(cv$n_dropped, 3L)
  expect_identical(cv$n_unknown_ids, 3L)
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  g <- suppressWarnings(suppressMessages(fit_gwr_model(d, "a", "z", adaptive = TRUE, bandwidth = 30)))
  expect_identical(g$info$n_dropped, 3L)
})

# ---------------------------------------------------------------------------
# 10.5
# ---------------------------------------------------------------------------

test_that("make_folds(block_kfold) returns the block design it built the folds from", {
  pts <- surf_test_points(n = 150)
  f <- make_folds(pts, k = 4, method = "block_kfold", seed = 1)
  asg <- f$assignment; pr <- f$params
  expect_true(all(c("row_id", "fold", "block_id") %in% names(asg)))
  expect_true(is.integer(asg$block_id))
  expect_false(anyNA(asg$block_id))
  # The block polygons, one per block still in the design, in the fold CRS.
  expect_s3_class(pr$blocks, "sf")
  expect_identical(pr$blocks$block_id, seq_len(nrow(pr$blocks)))
  expect_equal(nrow(pr$blocks), pr$blocks_used)
  expect_equal(sf::st_crs(pr$blocks), sf::st_crs(pts))
  expect_true(all(asg$block_id %in% pr$blocks$block_id))
  # Every point sits in the block it is labelled with.
  hit <- sf::st_intersects(pts, pr$blocks)
  expect_true(all(mapply(function(h, b) b %in% h, hit, asg$block_id)))
  # Points per block, indexed by block_id; the packing partitions the blocks;
  # a point's fold is the fold its block was packed into.
  expect_identical(pr$block_sizes, as.integer(table(factor(asg$block_id, levels = pr$blocks$block_id))))
  expect_equal(sum(pr$block_sizes), nrow(pts))
  expect_true(all(pr$block_sizes > 0L))                      # empties dropped
  expect_length(pr$fold_blocks, f$k)
  expect_identical(sort(unlist(pr$fold_blocks)), seq_len(pr$blocks_used))
  fold_of_block <- rep(seq_len(f$k), lengths(pr$fold_blocks))[order(unlist(pr$fold_blocks))]
  expect_identical(asg$fold, fold_of_block[asg$block_id])
  expect_identical(as.integer(table(asg$fold)),
                   vapply(pr$fold_blocks, function(b) sum(pr$block_sizes[b]), integer(1)))
  # Dropping the empty blocks renumbers the rest, so `source_row` is what ties
  # the returned design back to the layer the blocks came from.
  expect_true(is.integer(pr$blocks$source_row))
  expect_true(all(pr$blocks$source_row %in% seq_len(pr$n_blocks)))
  expect_false(is.unsorted(pr$blocks$source_row, strictly = TRUE))
  expect_length(unique(pr$blocks$source_row), nrow(pr$blocks))
  # Empty blocks kept: zeros in block_sizes, and the polygons include them.
  f0 <- make_folds(pts, k = 4, method = "block_kfold", seed = 1, block_nx = 8, block_ny = 8,
                   drop_empty_blocks = FALSE)
  expect_equal(nrow(f0$params$blocks), 64L)
  expect_length(f0$params$block_sizes, 64L)
  expect_true(any(f0$params$block_sizes == 0L))
  expect_equal(sum(f0$params$block_sizes), nrow(pts))
  # Nothing dropped, so source_row is the identity and the two counts agree.
  expect_identical(f0$params$blocks$source_row, 1:64)
  expect_identical(f0$params$n_blocks, 64L)
  expect_identical(f0$params$blocks_used, 64L)
  # The folds account for every block exactly once, empty ones included.
  # `blocks_used` used to be the highest block id a point fell in, so on a
  # layer whose points sit in one part of the grid it reported fewer blocks
  # than the design has, and the blocks above it were packed into no fold
  # while equally empty blocks below it were.
  expect_identical(sort(unlist(f0$params$fold_blocks)), seq_len(nrow(f0$params$blocks)))
  expect_length(unlist(f0$params$fold_blocks), nrow(f0$params$blocks))
})

test_that("the design stays consistent when the points reach only part of the grid", {
  # Points in one quadrant of an 8x8 grid, empty blocks kept: the highest
  # block id holding a point is far below the block count.
  set.seed(7)
  pts <- sf::st_as_sf(data.frame(x = runif(120, 0, 400), y = runif(120, 0, 400)),
                      coords = c("x", "y"), crs = 32632)
  pts$z <- rnorm(120)
  bnd <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(1000, 0), c(1000, 1000), c(0, 1000), c(0, 0)))), crs = 32632))
  f <- suppressWarnings(suppressMessages(
         make_folds(pts, k = 4, method = "block_kfold", seed = 1, boundary = bnd,
                    block_nx = 8, block_ny = 8, drop_empty_blocks = FALSE)))
  pr <- f$params
  expect_lt(max(f$assignment$block_id), nrow(pr$blocks))     # the case in question
  # One number for "how many blocks the design has", agreeing everywhere.
  expect_identical(pr$blocks_used, nrow(pr$blocks))
  expect_identical(pr$blocks_used, length(pr$block_sizes))
  expect_identical(pr$blocks_used, pr$n_blocks)              # nothing was dropped
  # Every block sits in exactly one fold, whether or not it holds a point.
  expect_identical(sort(unlist(pr$fold_blocks)), seq_len(nrow(pr$blocks)))
  expect_length(unlist(pr$fold_blocks), nrow(pr$blocks))
  expect_gt(sum(pr$block_sizes == 0L), 0L)
  # And the populated count is still available, without being confused for the
  # block count.
  expect_equal(sum(pr$block_sizes > 0L), length(unique(f$assignment$block_id)))
  # Packing the empty blocks moves no observation: each fold's size is the sum
  # of its blocks' sizes, and every point is in the fold its block was packed
  # into.
  expect_identical(as.integer(table(f$assignment$fold)),
                   vapply(pr$fold_blocks, function(b) sum(pr$block_sizes[b]), integer(1)))
  fold_of_block <- rep(seq_len(f$k), lengths(pr$fold_blocks))[order(unlist(pr$fold_blocks))]
  expect_identical(f$assignment$fold, fold_of_block[f$assignment$block_id])
  # Other methods carry no block column.
  fr <- make_folds(pts, k = 4, method = "random_kfold", seed = 1)
  expect_false("block_id" %in% names(fr$assignment))
  expect_null(fr$params$blocks)
})

test_that("supplied blocks come back as the design, numbered by their rows", {
  pts <- surf_test_points(n = 150)
  bnd <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(1000, 0), c(1000, 1000), c(0, 1000), c(0, 0)))), crs = 3857))
  seeds <- suppressWarnings(get_voronoi_seeds(bnd, method = "kmeans", n = 10,
                                              sample_points = pts, set_seed = 1))
  vor <- build_tessellation(seeds, boundary = bnd, method = "voronoi", quiet = TRUE)$cells
  f <- make_folds(pts, k = 3, method = "block_kfold", blocks = vor, seed = 1)
  pr <- f$params
  expect_true(pr$blocks_supplied)
  expect_equal(nrow(pr$blocks), nrow(vor))
  expect_identical(pr$blocks$block_id, seq_len(nrow(vor)))
  expect_equal(sum(pr$block_sizes), nrow(pts))
  # Each point is inside the supplied polygon its block_id names.
  hit <- sf::st_intersects(sf::st_transform(pts, sf::st_crs(pr$blocks)), pr$blocks)
  expect_true(all(mapply(function(h, b) b %in% h, hit, f$assignment$block_id)))

  # `block_id` is NOT the row of the supplied layer once empty blocks are
  # dropped -- a join by row position would mis-attribute every block after
  # the first gap -- so `source_row` carries the mapping. Only `block_id` and
  # `source_row` are kept as columns; the caller's own columns come back by
  # indexing their layer with it.
  expect_named(sf::st_drop_geometry(pr$blocks), c("block_id", "source_row"))
  expect_identical(pr$n_blocks, nrow(vor))
  # Every cell here holds a point, so nothing is dropped and the mapping is
  # the identity.
  expect_identical(pr$blocks$source_row, seq_len(nrow(vor)))
})

test_that("source_row maps the returned blocks back when the empty ones are dropped", {
  # Points in the left two fifths only, so three of the nine supplied blocks
  # hold nothing and are dropped -- which renumbers the six that remain.
  set.seed(3)
  pts <- sf::st_as_sf(data.frame(x = runif(120, 0, 400), y = runif(120, 0, 1000)),
                      coords = c("x", "y"), crs = 32632)
  pts$z <- rnorm(120)
  bnd <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(1000, 0), c(1000, 1000), c(0, 1000), c(0, 0)))), crs = 32632))
  blocks <- create_grid_polygons(bnd, target_cells = 9, type = "square")
  blocks$mylabel <- paste0("B", seq_len(nrow(blocks)))
  f  <- suppressWarnings(suppressMessages(
          make_folds(pts, k = 3, method = "block_kfold", blocks = blocks, seed = 1)))
  pr <- f$params

  expect_lt(nrow(pr$blocks), nrow(blocks))
  expect_identical(pr$n_blocks, nrow(blocks))
  expect_equal(nrow(pr$blocks), pr$blocks_used)
  # block_id is a fresh 1..n numbering, NOT the row of the supplied layer: a
  # join by row position would mis-attribute every block after the first gap.
  expect_identical(pr$blocks$block_id, seq_len(nrow(pr$blocks)))
  expect_false(identical(pr$blocks$source_row, pr$blocks$block_id))
  # source_row is what recovers them, with the caller's own columns.
  back <- blocks[pr$blocks$source_row, ]
  expect_identical(back$mylabel, paste0("B", pr$blocks$source_row))
  expect_identical(unlist(sf::st_equals(sf::st_geometry(pr$blocks),
                                        sf::st_geometry(back))),
                   seq_len(nrow(back)))
  # The blocks left out are exactly the ones no point fell in.
  left_out <- setdiff(seq_len(nrow(blocks)), pr$blocks$source_row)
  expect_gt(length(left_out), 0L)
  expect_true(all(lengths(sf::st_intersects(blocks[left_out, ], pts)) == 0L))
  # And every point's block_id still names the polygon it is inside.
  hit <- sf::st_intersects(pts, pr$blocks)
  expect_true(all(mapply(function(h, b) b %in% h, hit, f$assignment$block_id)))
})

# ---------------------------------------------------------------------------
# 10.7: the roster of one-liners
# ---------------------------------------------------------------------------

ra_bnd <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
  c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)))), crs = 32632))

test_that("assign_features_to_polygons() reports the features its tie-break decided", {
  pts <- ra_points(n = 150)
  grid <- create_grid_polygons(ra_bnd, target_cells = 16, type = "square")
  asg <- assign_features_to_polygons(pts, grid)
  ties <- attr(asg, "ties")
  expect_named(ties, c("n", "which", "rule", "n_rows"))
  expect_identical(ties$n_rows, nrow(asg))
  # The record names row positions, so it does not follow a subset.
  expect_null(attr(asg[1:10, ], "ties"))
  expect_identical(class(asg[1:10, ]), c("sf", "data.frame"))
  expect_identical(ties$n, 0L); expect_identical(ties$which, integer(0))
  expect_identical(ties$rule, "smallest_area")
  # Overlapping polygons: every feature inside a duplicated cell has a tie.
  ov <- rbind(grid, grid[1:4, ]); ov$poly_id <- seq_len(nrow(ov))
  lines <- capture_spatialkit_log(asg2 <- assign_features_to_polygons(pts, ov, tie_break = "first"))
  ties2 <- attr(asg2, "ties")
  in_dup <- which(lengths(sf::st_intersects(pts, grid[1:4, ])) > 0L)
  expect_identical(ties2$which, in_dup)
  expect_identical(ties2$n, length(in_dup))
  expect_identical(ties2$rule, "first")
  expect_true(log_has(lines, sprintf("%d of %d feature.s. fall inside more than one polygon", length(in_dup), nrow(pts))))
  # The result itself is what it always was: one row per feature.
  expect_equal(nrow(asg2), nrow(pts))
})

test_that("summarize_by_cell() records the ICCs a Kish request estimated, applied or not", {
  pts <- ra_points(n = 200)
  pts$val <- 0.05 * sf::st_coordinates(pts)[, 1] + rnorm(200, sd = 0.5)   # structured
  grid <- create_grid_polygons(ra_bnd, target_cells = 16, type = "square")
  asg <- assign_features_to_polygons(pts, grid)
  sm <- summarize_by_cell(asg, response_var = "val", predictor_vars = "z", deff = "kish")
  icc <- attr(sm, "icc")
  expect_named(icc, c("resp", "pred"))
  expect_equal(icc$resp, attr(sm, "deff_applied")$icc_resp)
  expect_true(is.finite(icc$pred))
  # Nothing applied (white-noise response, ICC clamped to 0): no deff_applied,
  # but the ICC is still on the result.
  sm0 <- summarize_by_cell(asg, response_var = "z", deff = "kish")
  expect_null(attr(sm0, "deff_applied"))
  expect_identical(attr(sm0, "icc"), list(resp = 0, pred = NA_real_))
  # Survives the cells_sf join; absent without a Kish request.
  smc <- summarize_by_cell(asg, response_var = "z", deff = "kish", cells_sf = grid)
  expect_identical(attr(smc, "icc"), list(resp = 0, pred = NA_real_))
  expect_null(attr(summarize_by_cell(asg, response_var = "z"), "icc"))
})

test_that("the variogram design effect at the cell row count is kept beside the applied one", {
  skip_if_not_installed("gstat")
  set.seed(3)
  n <- 240
  x <- runif(n, 0, 100); y <- runif(n, 0, 100)
  d <- as.matrix(stats::dist(cbind(x, y)))
  z <- as.numeric(t(chol(exp(-d / 15) + diag(0.2 + 1e-8, n))) %*% rnorm(n))
  pts <- sf::st_as_sf(data.frame(x = x, y = y, z = z, w = rnorm(n)), coords = c("x", "y"), crs = 32632)
  grid <- create_grid_polygons(ra_bnd, target_cells = 9, type = "square")
  asg <- assign_features_to_polygons(pts, grid)
  sac <- estimate_sac_range(pts, "z")
  skip_if(is.na(sac), "no identified range on this draw")
  sm <- summarize_by_cell(asg, response_var = "z", deff = "variogram", sac = sac)
  da <- attr(sm, "deff_applied")
  expect_identical(da$method, "variogram")
  expect_length(da$deff_rows, nrow(sm))
  # The response is complete, so the row-count deff is the applied one.
  expect_equal(da$deff_rows, da$deff)
  # With a missing response the two differ where rows were lost, and the
  # join to cells_sf realigns both alike (an empty cell gets NA in each).
  asg2 <- asg; asg2$z[asg2$poly_id == asg2$poly_id[1]][1:2] <- NA
  smc <- summarize_by_cell(asg2[asg2$poly_id != asg2$poly_id[nrow(asg2)], ],
                           response_var = "z", deff = "variogram", sac = sac, cells_sf = grid)
  dac <- attr(smc, "deff_applied")
  expect_length(dac$deff_rows, nrow(smc)); expect_length(dac$deff, nrow(smc))
  expect_identical(is.na(dac$deff_rows), is.na(dac$deff))
  expect_true(any(is.na(dac$deff_rows)))
  differ <- which(is.finite(dac$deff) & abs(dac$deff_rows - dac$deff) > 1e-12)
  expect_true(length(differ) >= 1L)
})

test_that("ensure_projected() attaches the projections it scored", {
  ll <- sf::st_as_sf(data.frame(lon = c(9.1, 9.2, 9.15), lat = c(48.7, 48.8, 48.75)),
                     coords = c("lon", "lat"), crs = 4326)
  p <- ensure_projected(ll)
  ch <- attr(p, "crs_choice")
  expect_s3_class(ch, "data.frame")
  expect_named(ch, c("name", "crs", "distance_error", "chosen"))
  expect_equal(nrow(ch), 1L)                       # the ordinary path: the zone alone
  expect_match(ch$name, "^UTM zone 32")
  expect_true(ch$chosen)
  expect_true(is.finite(ch$distance_error) && ch$distance_error < 0.01)
  # A continental extent: the zone and the equal-area candidates, one chosen,
  # and the chosen one is the least distorting of those measured.
  wide <- sf::st_as_sf(data.frame(lon = c(-120, -70, -95), lat = c(30, 48, 40)),
                       coords = c("lon", "lat"), crs = 4326)
  w <- suppressWarnings(ensure_projected(wide))
  cw <- attr(w, "crs_choice")
  expect_gte(nrow(cw), 2L)
  expect_identical(sum(cw$chosen), 1L)
  expect_equal(which(cw$chosen), which.min(cw$distance_error))
  expect_identical(sf::st_crs(w)$input, cw$crs[cw$chosen])
  # purpose = "area": equal-area candidates only, no zone.
  a <- ensure_projected(ll, purpose = "area")
  ca <- attr(a, "crs_choice")
  expect_false(any(grepl("^UTM", ca$name)))
  expect_identical(sum(ca$chosen), 1L)
  # A layer straddling the antimeridian gets a local equal-area projection
  # chosen for it, so it reports that choice too -- it used to be the one
  # path that picked a projection and said nothing.
  wrap <- sf::st_as_sf(data.frame(lon = c(runif(25, 179, 180), runif(25, -180, -179)),
                                  lat = runif(50, -1, 1)),
                       coords = c("lon", "lat"), crs = 4326)
  w2 <- suppressWarnings(ensure_projected(wrap))
  cwrap <- attr(w2, "crs_choice")
  expect_s3_class(cwrap, "data.frame")
  expect_equal(nrow(cwrap), 1L)
  expect_true(cwrap$chosen)
  expect_match(cwrap$name, "^Lambert azimuthal equal-area")
  expect_identical(sf::st_crs(w2)$input, cwrap$crs)
  # No local projection was chosen: nothing attached.
  expect_null(attr(ensure_projected(ll, target_crs = 3035), "crs_choice"))
  expect_null(attr(ensure_projected(sf::st_transform(ll, 32632)), "crs_choice"))
  glob <- sf::st_as_sf(data.frame(lon = c(-170, 0, 170), lat = c(0, 10, -10)),
                       coords = c("lon", "lat"), crs = 4326)
  expect_null(attr(suppressWarnings(ensure_projected(glob)), "crs_choice"))
})

test_that("the residual Moran's I result prints its diagnostics, not its weight matrix", {
  pts <- surf_test_points(n = 120)
  mi  <- residual_morans_i(lm_spatial_fit(pts, "z", "w"))
  expect_s3_class(mi, "morans_i")
  # The weight matrix is n x n, and a bare list autoprints every element: this
  # is the function's own documented example, and it emitted 1246 lines.
  out <- capture.output(print(mi))
  expect_lt(length(out), 10L)
  expect_match(out[1], "^Residual Moran's I = ")
  expect_false(any(grepl("dgCMatrix|^\\[1\\]|\\.\\.\\.", out[-5])))
  # The representation depends on the optional backends: the sparse path needs
  # FNN *and* Matrix, so a no-Suggests install gets a dense base matrix and the
  # neighbour count is deliberately not computed for it (counting non-zeros on
  # a dense n x n matrix would allocate another one).
  expect_match(out[5], "^  weights: 120 x 120 ")
  # Not retained by default, and the print line says so rather than leaving a
  # NULL `weights` unexplained.
  expect_null(mi$weights)
  expect_match(out[5], "not retained, keep_weights = TRUE to keep it", fixed = TRUE)
  expect_false(isTRUE(mi$weights_summary$kept))
  expect_identical(mi$weights_summary$n, 120L)
  if (identical(mi$weights_summary$storage, "dgCMatrix")) {
    expect_match(out[5], "8 neighbour\\(s\\) per row")
    expect_identical(mi$weights_summary$neighbours, c(8L, 8L))
  } else {
    expect_match(out[5], "^  weights: 120 x 120 matrix; not retained")
    expect_true(all(is.na(mi$weights_summary$neighbours)))
  }
  expect_match(out[3], "null: residual moments, exact for these residuals; n = 120, df = 118, design rank 2",
               fixed = TRUE)
  # print() returns its argument invisibly, and changes nothing.
  expect_identical(withVisible(print(mi))$visible, FALSE)
  expect_identical(capture.output(print(mi)), out)
  # The randomisation null reports itself as approximate and has no rank.
  outr <- capture.output(print(residual_morans_i(lm_spatial_fit(pts, "z", "w"), null = "randomisation")))
  expect_match(outr[3], "randomisation moments, approximate for these residuals", fixed = TRUE)
  expect_false(grepl("design rank", outr[3]))
  # The class is the only change to the object: every existing accessor works,
  # and `[` drops it, so subsetting a few names still yields a plain list.
  expect_true(is.list(mi))
  expect_identical(class(mi[c("observed", "z")]), "list")
  expect_equal(unlist(mi[c("observed", "z")]), c(observed = mi$observed, z = mi$z))
  # A dense weight matrix is described without being counted element-wise.
  W <- as.matrix(residual_morans_i(lm_spatial_fit(pts, "z", "w"),
                                   keep_weights = TRUE)$weights)
  dw <- spatialkit:::.morans_weights_desc(W)
  expect_match(dw, "^120 x 120 matrix$")
  expect_identical(spatialkit:::.morans_weights_desc(NULL), "not available")
})

test_that("residual_morans_i() returns its weights, the kurtosis, the design rank and exactness", {
  pts <- surf_test_points(n = 100)
  fit <- lm_spatial_fit(pts, "z", "w")
  mi <- residual_morans_i(fit, keep_weights = TRUE)
  expect_true(all(c("weights", "kurtosis", "p", "exact") %in% names(mi)))
  expect_equal(dim(mi$weights), c(100L, 100L))
  expect_true(mi$weights_summary$kept)
  # The default drops the matrix and nothing else: every other component is
  # identical, and the object is a fraction of the size.
  lite <- residual_morans_i(fit)
  expect_null(lite$weights)
  keys <- setdiff(names(mi), c("weights", "weights_summary"))
  expect_equal(lite[keys], mi[keys])
  expect_lt(as.numeric(object.size(lite)), as.numeric(object.size(mi)) / 4)
  expect_identical(lite$weights_summary[c("n", "storage", "neighbours")],
                   mi$weights_summary[c("n", "storage", "neighbours")])
  expect_true(all(abs(Matrix::rowSums(mi$weights) - 1) < 1e-8) ||
              all(abs(rowSums(as.matrix(mi$weights)) - 1) < 1e-8))
  e <- stats::residuals(fit); ec <- e - mean(e)
  expect_equal(mi$kurtosis, mean(ec^4) / mean(ec^2)^2)
  # OLS residuals on the rebuilt design: the residual null, exact, rank 2.
  expect_identical(mi$null, "residual")
  expect_true(mi$exact)
  expect_identical(mi$p, 2L)
  # Forced randomisation: approximate, no design rank.
  mr <- residual_morans_i(fit, null = "randomisation")
  expect_false(mr$exact); expect_true(is.na(mr$p))
  expect_equal(mr$kurtosis, mi$kurtosis)
  # A supplied weights matrix comes back as used (diagonal zeroed).
  W <- as.matrix(mi$weights); diag(W) <- 0.5
  ms <- suppressWarnings(residual_morans_i(fit, weights = W, keep_weights = TRUE))
  expect_equal(unname(diag(as.matrix(ms$weights))), rep(0, 100))
  # A supplied matrix is not exempt from the default either.
  expect_null(suppressWarnings(residual_morans_i(fit, weights = W))$weights)
})

test_that("fit_gwr_model() keeps the mask of non-finite local coefficients", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  pts <- surf_test_points(n = 120)
  pts$v <- rnorm(120)
  fit <- suppressWarnings(suppressMessages(fit_gwr_model(pts, "z", c("w", "v"), adaptive = TRUE, bandwidth = 40)))
  m <- fit$info$nonfinite_coef
  expect_true(is.matrix(m) && is.logical(m))
  expect_equal(dim(m), c(120L, 3L))
  expect_identical(colnames(m), c("Intercept", "w", "v"))
  expect_identical(sum(apply(m, 1L, any)), fit$info$n_local_singular)
  expect_false(any(m))
})

test_that("the Bayesian fit names the parameters that fail its convergence checks", {
  src <- paste(deparse(body(fit_bayesian_spatial_model)), collapse = " ")
  expect_true(grepl("rhat_failed <- rhat_vals[bad_rhat]", src, fixed = TRUE))
  expect_true(grepl("neff_failed <- neff_vals[low_neff]", src, fixed = TRUE))
})

test_that("determine_optimal_levels() reports the elbow and the failed k in its diagnostics", {
  set.seed(5)
  ctr <- cbind(runif(16, 0, 1000), runif(16, 0, 1000))
  g   <- sample(16, 640, TRUE)
  xy  <- ctr[g, ] + matrix(rnorm(2 * 640, sd = 20), 640)
  pts <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2], w = rnorm(640)),
                      coords = c("x", "y"), crs = 32632)
  pts$z <- pts$w + rnorm(640)
  out <- determine_optimal_levels(pts, max_levels = 30, response_var = "z",
                                  predictor_vars = "w", criterion = "combined")
  d <- attr(out, "diagnostics")
  expect_false(is.null(d))
  expect_true(is.numeric(d$knee_k) && length(d$knee_k) == 1L)
  expect_true(d$knee_k %in% seq_along(d$wss))
  expect_true(all(d$eval_ks >= d$knee_k - max(4L, 3L) & d$eval_ks <= d$knee_k + max(4L, 3L)))
  expect_identical(d$failed_k, integer(0))
  # The geometric path stays a plain integer vector, as 5.1 settled.
  expect_null(attributes(determine_optimal_levels(pts, max_levels = 30)))
})

test_that("build_tessellation() records the points it snapped to a nearest cell", {
  pts <- surf_test_points(n = 60)
  # One point a hair outside the boundary (snapped), one well outside (NA).
  extra <- sf::st_as_sf(data.frame(x = c(1000.002, 1050), y = c(500, 500), w = 0, z = 0),
                        coords = c("x", "y"), crs = 3857)
  all_pts <- rbind(pts, extra)
  bnd <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(1000, 0), c(1000, 1000), c(0, 1000), c(0, 0)))), crs = 3857))
  tess <- build_tessellation(all_pts, boundary = bnd, method = "square", approx_n_cells = 16, quiet = TRUE)
  sn <- tess$params$snapped
  expect_named(sn, c("n", "which", "distance"))
  expect_identical(sn$n, 1L)
  expect_identical(sn$which, 61L)
  expect_equal(sn$distance, 0.002, tolerance = 1e-6)
  expect_false(is.na(tess$index[61]))
  expect_true(is.na(tess$index[62]))
  # `index` itself stays a bare integer vector.
  expect_null(attributes(tess$index))
  # The other methods carry the same record.
  vor <- build_tessellation(pts, boundary = bnd, method = "voronoi", quiet = TRUE)
  expect_named(vor$params$snapped, c("n", "which", "distance"))
  expect_null(attributes(vor$index))
})

test_that("area_of_applicability() returns the scaling and the fence's outlier count", {
  pts <- surf_test_points(n = 120)
  pts$v <- rnorm(120)
  new <- surf_test_points(n = 40, seed = 9); new$v <- rnorm(40)
  res <- area_of_applicability(new, train_sf = pts, predictor_vars = c("w", "v"))
  expect_named(res$scaling, c("center", "scale"))
  expect_equal(res$scaling$center, c(w = mean(pts$w), v = mean(pts$v)))
  expect_equal(res$scaling$scale, c(w = sd(pts$w), v = sd(pts$v)))
  di <- res$train_DI[is.finite(res$train_DI)]
  fence <- stats::quantile(di, 0.75, names = FALSE) + 1.5 * stats::IQR(di)
  expect_identical(res$n_outliers, as.integer(sum(di > fence)))
  # The threshold is the largest training DI at or below that fence.
  expect_equal(res$threshold, max(di[di <= fence]))
})

test_that("get_voronoi_seeds(kmeans) returns the clustering behind the seeds", {
  pts <- ra_points(n = 120)
  seeds <- get_voronoi_seeds(ra_bnd, method = "kmeans", n = 6, sample_points = pts, set_seed = 1)
  km <- attr(seeds, "kmeans")
  expect_named(km, c("cluster", "rows", "size", "withinss", "tot_withinss", "iter", "nstart"))
  expect_length(km$cluster, 120L)
  expect_identical(km$rows, 1:120)
  expect_true(all(km$cluster %in% seeds$seed_id))
  expect_identical(as.integer(table(factor(km$cluster, levels = 1:6))), km$size)
  expect_equal(sum(km$size), 120L)
  expect_equal(sum(km$withinss), km$tot_withinss)
  expect_identical(km$nstart, 10L)
  # The seed positions are the cluster centroids.
  xy <- sf::st_coordinates(pts)
  ctr <- t(sapply(1:6, function(k) colMeans(xy[km$cluster == k, , drop = FALSE])))
  expect_equal(unname(sf::st_coordinates(seeds)), unname(ctr), tolerance = 1e-6)
  # A dropped (empty) point is absent from `rows`.
  bad <- pts; sf::st_geometry(bad)[5] <- sf::st_sfc(sf::st_point(), crs = 32632)
  seeds2 <- suppressWarnings(get_voronoi_seeds(ra_bnd, method = "kmeans", n = 6,
                                               sample_points = bad, set_seed = 1))
  expect_identical(attr(seeds2, "kmeans")$rows, setdiff(1:120, 5L))
  expect_null(attr(get_voronoi_seeds(ra_bnd, method = "random", n = 5, set_seed = 1), "kmeans"))
})

test_that("gwr_model_selection() says how its criterion column was found", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  pts <- surf_test_points(n = 100)
  pts$v <- rnorm(100)
  sel <- suppressWarnings(suppressMessages(
    gwr_model_selection(pts, "z", c("w", "v"), adaptive = TRUE, bandwidth = 40)))
  expect_true(all(c("criterion_by_name", "criterion_column", "criterion_verified") %in% names(sel)))
  expect_type(sel$criterion_by_name, "logical")
  expect_true(sel$criterion_column %in% 1:4)
  # Verified whenever the column was named or the table had its four
  # documented columns; unverified is the case the log flags.
  raw_ok <- is.matrix(sel$raw[[2]]) && ncol(sel$raw[[2]]) == 4L
  expect_identical(sel$criterion_verified, isTRUE(sel$criterion_by_name) || raw_ok)
})
