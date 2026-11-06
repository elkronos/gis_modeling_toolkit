# Blocked and buffered folds under the package's tie rule (BUG-100, BUG-101,
# BUG-202): a point on a block edge goes to the lower block, a neighbour at
# the buffer distance is excluded, and a point as near two blocks goes to the
# lower row, whatever the last bits of the coordinates say.  So the folds are
# the same in another unit or with another origin.

.ft_crs_m    <- function() sf::st_crs(26915)
.ft_crs_usft <- function() sf::st_crs("+proj=utm +zone=15 +datum=NAD83 +units=us-ft +no_defs")
.ft_tm <- function(x0, y0, units = "m") sf::st_crs(sprintf(paste0(
  "+proj=tmerc +lat_0=0 +lon_0=-93 +k=0.9996 +x_0=%s +y_0=%s +datum=NAD83 ",
  "+units=%s +no_defs"), x0, y0, units))
.ft_frames <- function(pm) list(m = pm, usft = sf::st_transform(pm, .ft_crs_usft()),
                                m01 = sf::st_transform(pm, .ft_tm(500000.1, 0.1)))
.ft_len <- function(l, metres)
  metres * (if (isTRUE(sf::st_crs(l) == .ft_crs_usft())) 3937 / 1200 else 1)

# A lattice of nx by ny points at 250 m from (450000, 4550000), with a
# response z; i and j are the column and row.
.ft_lattice <- function(nx, ny, jitter = 0) {
  g <- expand.grid(i = 0:(nx - 1), j = 0:(ny - 1))
  x <- 4.5e5 + g$i * 250; y <- 4.55e6 + g$j * 250
  if (jitter > 0) {
    x <- x + withr::with_seed(99, stats::runif(length(x), -jitter, jitter))
    y <- y + withr::with_seed(98, stats::runif(length(y), -jitter, jitter))
  }
  z <- withr::with_seed(1, sin((x - min(x)) / 1500) + cos((y - min(y)) / 1250) +
                          stats::rnorm(length(x), 0, 0.5))
  sf::st_as_sf(data.frame(x = x, y = y, z = z, i = g$i, j = g$j),
               coords = c("x", "y"), crs = .ft_crs_m())
}

.ft_quiet <- function(expr) suppressWarnings(suppressMessages(
  logger::with_log_threshold(expr, threshold = logger::FATAL, namespace = "spatialkit",
                             index = 2)))

.ft_square <- function(x0, y0, s) sf::st_polygon(list(rbind(
  c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s), c(x0, y0 + s), c(x0, y0))))

# 3 x 2 squares of 3000 m over the 37 x 25 lattice, numbered from the
# bottom left along each row, as st_make_grid() numbers a grid.
.ft_blocks <- function() {
  b <- expand.grid(ci = 0:2, cj = 0:1)
  sf::st_sf(geometry = sf::st_sfc(lapply(seq_len(nrow(b)), function(r)
    .ft_square(4.5e5 + 3000 * b$ci[r], 4.55e6 + 3000 * b$cj[r], 3000)),
    crs = .ft_crs_m()))
}

test_that("block_kfold puts a lattice point on a grid line in the lower block in every frame", {
  lat <- .ft_lattice(37, 25)
  F <- .ft_frames(lat)
  for (nb in list(c(3L, 2L), c(4L, 3L))) {
    # The lower block: a point on a vertical line goes left, on a horizontal
    # line down.  Blocks are 9000 / nx by 6000 / ny metres.
    col <- pmax(1, ceiling(lat$i * 250 / (9000 / nb[1])))
    row <- pmax(1, ceiling(lat$j * 250 / (6000 / nb[2])))
    want <- as.integer((row - 1) * nb[1] + col)
    res <- lapply(F, function(l) .ft_quiet(make_folds(
      l, k = 4, method = "block_kfold", seed = 1, block_nx = nb[1], block_ny = nb[2])))
    # Before, 469 rows (3 x 2, origin moved) and 448 / 824 rows (4 x 3, US
    # feet / origin moved) were in another fold than in metres.
    for (f in names(F)) {
      expect_identical(res[[f]]$assignment$block_id, want, info = paste(f, nb[1]))
      expect_identical(res[[f]]$assignment$fold, res$m$assignment$fold,
                       info = paste(f, nb[1]))
    }
  }
  # block_size: the same blocks in metres and in US feet.
  res <- lapply(F, function(l) .ft_quiet(make_folds(
    l, k = 4, method = "block_kfold", seed = 1, block_size = .ft_len(l, 1500))))
  expect_identical(res$usft$assignment$fold, res$m$assignment$fold)
  expect_identical(res$m01$assignment$fold, res$m$assignment$fold)
  expect_identical(res$m01$assignment$block_id, res$m$assignment$block_id)
})

test_that("block_kfold on a jittered lattice is unchanged in every frame", {
  F <- .ft_frames(.ft_lattice(37, 25, jitter = 60))
  res <- lapply(F, function(l) .ft_quiet(make_folds(
    l, k = 4, method = "block_kfold", seed = 1, block_nx = 4L, block_ny = 3L)))
  expect_identical(as.integer(table(res$m$assignment$block_id)),
                   c(76L, 79L, 74L, 80L, 79L, 70L, 69L, 82L, 78L, 83L, 78L, 77L))
  expect_identical(res$usft$assignment$fold, res$m$assignment$fold)
  expect_identical(res$m01$assignment$fold, res$m$assignment$fold)
})

test_that("the block counts of the leakage advice read an edge point as the folds do", {
  F <- .ft_frames(.ft_lattice(37, 25))
  occ <- vapply(F, function(l) {
    bb <- sf::st_bbox(l)
    c(slow = .block_occupied(l, sf::st_as_sfc(bb), 36, 24),
      fast = .block_occupied_fast(sf::st_coordinates(l), bb, 36, 24))
  }, numeric(2))
  # 36 x 24 cells of 250 m: under the lower-block rule every cell holds a
  # point.  Before, 576 / 612 in US feet and 680 / 681 with the origin moved.
  expect_true(all(occ == 864))
})

test_that("the fast block count puts a point on a grid line in the lower cell, as the polygon count does", {
  # A sparse layer on the grid lines of the whole lattice's box, which does
  # not fill the grid: the cell a point on a line goes to sets the count.
  lat <- .ft_lattice(37, 25)
  full <- .ft_frames(lat)
  for (keep in list(lat$i %in% c(0, 12, 24) & lat$j %in% c(0, 12),
                    lat$i %% 12 == 0 | lat$j %% 12 == 0,
                    lat$i %in% c(0, 18) & lat$j %in% c(0, 12, 24))) {
    F <- .ft_frames(lat[keep, ])
    for (f in names(F)) {
      bb <- sf::st_bbox(full[[f]]); xy <- sf::st_coordinates(F[[f]])
      reg <- sf::st_as_sfc(bb)
      for (nb in list(c(2, 2), c(3, 2), c(4, 4), c(6, 4))) {
        slow <- .block_occupied(F[[f]], reg, nb[1], nb[2])
        fast <- .block_occupied_fast(xy, bb, nb[1], nb[2])
        expect_identical(fast, slow, info = paste(f, sum(keep), nb[1], nb[2]))
      }
    }
  }
  # Four points at the corners of the left and middle lines of a 2 x 2 grid
  # are in the lower left cell only.
  q <- .ft_frames(lat[lat$i %in% c(0, 18) & lat$j %in% c(0, 12), ])
  for (f in names(q)) {
    bb <- sf::st_bbox(full[[f]])
    expect_identical(.block_occupied_fast(sf::st_coordinates(q[[f]]), bb, 2, 2), 1L,
                     info = f)
  }
})

test_that("a lattice point on the boundary is inside it in every frame", {
  lat <- .ft_lattice(37, 25)
  frs <- list(m = .ft_crs_m(), usft = .ft_crs_usft(), m01 = .ft_tm(500000.1, 0.1),
              km = .ft_tm(500, 0, "km"), o5 = .ft_tm(1e6 + 0.3, 1e7 + 0.9))
  run <- function(pts, bnd, crs, ...) {
    w <- character(0)
    fo <- withCallingHandlers(
      suppressMessages(logger::with_log_threshold(
        make_folds(sf::st_transform(pts, crs), k = 3, method = "block_kfold", seed = 1,
                   boundary = sf::st_transform(bnd, crs), ...),
        threshold = logger::FATAL, namespace = "spatialkit", index = 2)),
      warning = function(cnd) {
        w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
      })
    a <- sum(as.numeric(sf::st_area(sf::st_transform(fo$params$blocks, .ft_crs_m()))))
    list(block = fo$assignment$block_id, fold = fo$assignment$fold, warn = w, area = a)
  }
  # A triangle whose edges pass through 72 of the 469 points it holds, and
  # the bounding box of the whole lattice.
  tri <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(4.5e5, 4.55e6), c(4.5e5 + 9000, 4.55e6), c(4.5e5, 4.556e6), c(4.5e5, 4.55e6)))),
    crs = .ft_crs_m()))
  ins <- lat[sf::st_intersects(lat, tri, sparse = FALSE)[, 1], ]
  box <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(lat)))
  for (case in list(list(pts = ins, bnd = tri, nb = c(4L, 3L), area = 2.7e7),
                    list(pts = lat, bnd = box, nb = c(3L, 2L), area = 5.4e7))) {
    res <- lapply(frs, function(crs) run(case$pts, case$bnd, crs,
                                         block_nx = case$nb[1], block_ny = case$nb[2]))
    # Before, 8 (US feet), 49 (origin moved) and 50 (km) of the 469 points of
    # the triangle, and 71 or 72 of the 925 of the box, were outside, and the
    # block region was extended to the points' bounding box.
    for (f in names(frs)) {
      expect_false(any(grepl("outside `boundary`", res[[f]]$warn, fixed = TRUE)),
                   info = paste(f, nrow(case$pts)))
      expect_equal(res[[f]]$area, case$area, tolerance = 1e-6, info = paste(f, nrow(case$pts)))
      expect_identical(res[[f]]$block, res$m$block, info = paste(f, nrow(case$pts)))
      expect_identical(res[[f]]$fold, res$m$fold, info = paste(f, nrow(case$pts)))
    }
  }
  # A point well outside the boundary is still warned about.
  out <- rbind(ins, sf::st_as_sf(data.frame(x = 4.5e5 + 8000, y = 4.555e6, z = 0, i = NA_integer_,
                                            j = NA_integer_),
                                 coords = c("x", "y"), crs = .ft_crs_m()))
  r <- run(out, tri, .ft_crs_m(), block_nx = 4L, block_ny = 3L)
  expect_true(any(grepl("1 of 470 points fall outside `boundary`", r$warn, fixed = TRUE)))
})

test_that("the split halves of select_on = \"split\" are the same in every frame", {
  F <- .ft_frames(.ft_lattice(40, 30))
  sp <- lapply(F, function(l) sort(.ft_quiet(.spatial_half_split(l, seed = 123L))$selection))
  # Before, the selection half moved by 786 rows with the origin moved.
  expect_identical(sp$usft, sp$m)
  expect_identical(sp$m01, sp$m)
})

test_that("a lattice point on the shared edge of supplied blocks goes to the lower row in every frame", {
  lat <- .ft_lattice(37, 25)
  blk <- .ft_blocks()
  col <- pmax(1, ceiling(lat$i / 12)); row <- pmax(1, ceiling(lat$j / 12))
  want <- as.integer((row - 1) * 3 + col)
  for (f in c("m", "usft", "m01")) {
    crs <- switch(f, m = .ft_crs_m(), usft = .ft_crs_usft(), m01 = .ft_tm(500000.1, 0.1))
    fo <- .ft_quiet(make_folds(sf::st_transform(lat, crs), k = 3, method = "block_kfold",
                               seed = 1, blocks = sf::st_transform(blk, crs)))
    expect_identical(fo$assignment$block_id, want, info = f)
  }
})

test_that("a point outside the supplied blocks as near two of them goes to the lower row in every frame", {
  blk <- .ft_blocks()
  inner <- .ft_lattice(37, 25, jitter = 60)
  away <- !sf::st_is_within_distance(
    inner, sf::st_cast(sf::st_union(blk), "MULTILINESTRING"), 1, sparse = FALSE)[, 1]
  inner <- inner[sf::st_intersects(inner, sf::st_union(blk), sparse = FALSE)[, 1] & away, ]
  # 72 points 13.7 to 250 m outside the design, below, above and beside the
  # corners two blocks share: each as far from both.
  xs <- 4.5e5 + c(3000, 6000); dd <- seq(13.7, 250, length.out = 12)
  out <- rbind(data.frame(x = rep(xs, each = 12), y = 4.55e6 - dd),
               data.frame(x = rep(xs, each = 12), y = 4.55e6 + 6000 + dd),
               data.frame(x = 4.5e5 - dd, y = 4.553e6),
               data.frame(x = 4.5e5 + 9000 + dd, y = 4.553e6))
  out <- sf::st_as_sf(out, coords = c("x", "y"), crs = .ft_crs_m())
  out$z <- 0; out$i <- NA_integer_; out$j <- NA_integer_
  pts <- rbind(inner, out); io <- (nrow(inner) + 1L):nrow(pts)
  # The lower of the two blocks at each corner.
  want_out <- rep(c(1L, 2L, 4L, 5L, 1L, 3L), each = 12)
  frs <- list(m = .ft_crs_m(), usft = .ft_crs_usft(), ft = .ft_tm(500000 / 0.3048, 0, "ft"),
              km = .ft_tm(500, 0, "km"), o1 = .ft_tm(500000.1, 0.1),
              o4 = .ft_tm(0, 0), o5 = .ft_tm(1e6 + 0.3, 1e7 + 0.9))
  res <- lapply(frs, function(crs) {
    w <- character(0)
    fo <- withCallingHandlers(
      suppressMessages(logger::with_log_threshold(
        make_folds(sf::st_transform(pts, crs), k = 3, method = "block_kfold", seed = 1,
                   blocks = sf::st_transform(blk, crs)),
        threshold = logger::FATAL, namespace = "spatialkit", index = 2)),
      warning = function(cnd) {
        w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
      })
    list(block = fo$assignment$block_id, fold = fo$assignment$fold, warn = w)
  })
  # Before, 3 to 10 of the 72 points went to the other block outside metres
  # and US feet, and 3 to 627 rows to another fold.
  for (f in names(frs)) {
    expect_identical(res[[f]]$block[io], want_out, info = f)
    expect_identical(res[[f]]$block[-io], res$m$block[-io], info = f)
    expect_identical(res[[f]]$fold, res$m$fold, info = f)
    # The far-point warning and its count are as before.
    expect_true(any(grepl(sprintf("72 of %d points fall outside every supplied block",
                                  nrow(pts)), res[[f]]$warn, fixed = TRUE)), info = f)
  }
})

test_that("a point outside the supplied blocks nearer one of them takes that one", {
  blk <- .ft_blocks()
  inner <- .ft_lattice(37, 25, jitter = 60)
  inner <- inner[sf::st_intersects(inner, sf::st_union(blk), sparse = FALSE)[, 1], ]
  # 40 m below block 2, 10 m right of the corner it shares with block 1.
  out <- sf::st_as_sf(data.frame(x = 4.53e5 + 10, y = 4.55e6 - 40, z = 0,
                                 i = NA_integer_, j = NA_integer_),
                      coords = c("x", "y"), crs = .ft_crs_m())
  pts <- rbind(inner, out)
  for (crs in list(.ft_crs_m(), .ft_crs_usft(), .ft_tm(500000.1, 0.1))) {
    fo <- .ft_quiet(make_folds(sf::st_transform(pts, crs), k = 3, method = "block_kfold",
                               seed = 1, blocks = sf::st_transform(blk, crs)))
    expect_identical(fo$assignment$block_id[nrow(pts)], 2L)
  }
})

test_that("buffered_loo excludes every neighbour at the buffer distance in every frame", {
  lat <- .ft_lattice(20, 15)
  F <- .ft_frames(lat)
  for (b in c(500, 750, 600)) {
    # Excluded: squared lattice distance (in spacings) at most (b / 250)^2.
    want <- lapply(seq_len(nrow(lat)), function(r)
      which((lat$i - lat$i[r])^2 + (lat$j - lat$j[r])^2 > (b / 250)^2))
    for (f in names(F)) {
      fo <- .ft_quiet(make_folds(F[[f]], k = 1, method = "buffered_loo",
                                 buffer = .ft_len(F[[f]], b)))
      tr <- lapply(fo$folds, function(s) sort(as.integer(s$train)))
      # Before, buffer 500 m: 295 (US feet) and 253 (origin moved) of 300
      # training sets differed from the metre one; 750 m: 0 and 255.
      expect_identical(tr, want, info = paste(f, b))
    }
  }
})

test_that("buffered_loo on a jittered lattice is unchanged", {
  F <- .ft_frames(.ft_lattice(20, 15, jitter = 60))
  nt <- vapply(F, function(l) mean(vapply(.ft_quiet(make_folds(
    l, k = 1, method = "buffered_loo", buffer = .ft_len(l, 500)))$folds,
    function(s) length(s$train), 1L)), numeric(1))
  expect_equal(unname(nt), rep(289.24, 3), tolerance = 1e-12)
})

test_that(".lowest_nearest() keeps the lowest column within the tie rule", {
  d <- rbind(c(5, 5 * (1 + 1e-15), 7),
             c(5 * (1 + 1e-15), 5, 7),
             c(9, 3, 3 + 1e-6),
             c(NA, NA, NA),
             c(NA, 4, 4),
             c(Inf, Inf, 2))
  expect_identical(.lowest_nearest(d), c(1L, 1L, 2L, NA, 2L, 3L))
  # The absolute margin joins a distance that far over the least.
  expect_identical(.lowest_nearest(rbind(c(3 + 1e-6, 3)), tol = 2e-6), 1L)
  expect_identical(.lowest_nearest(rbind(c(3 + 1e-3, 3)), tol = 2e-6), 2L)
})

test_that(".edge_join() reads a point a hair over a shared edge as on it", {
  two <- function(crs) sf::st_sfc(.ft_square(0, 0, 1), .ft_square(1, 0, 1), crs = crs)
  p <- function(crs) sf::st_sfc(sf::st_point(c(1 + 1e-12, 0.5)), sf::st_point(c(1.5, 0.5)),
                                sf::st_point(c(3, 0.5)), crs = crs)
  # Planar and CRS-less: buffered by .edge_tol().
  id <- .edge_join(p(sf::NA_crs_), two(sf::NA_crs_))
  expect_identical(as.vector(id), c(1L, 2L, NA))
  expect_identical(lengths(attr(id, "hits")), c(2L, 1L, 0L))
  # In degrees with s2 on: the margin as an arc of the sphere.
  expect_identical(as.vector(.edge_join(p(4326), two(4326))), c(1L, 2L, NA))
  # No margin when every coordinate is 0.
  z <- sf::st_sfc(sf::st_point(c(0, 0)))
  expect_identical(as.vector(.edge_join(z, z)), 1L)
  # Fewer points than polygons: the points are buffered, with the same result.
  id1 <- .edge_join(p(sf::NA_crs_)[1], two(sf::NA_crs_))
  expect_identical(as.vector(id1), 1L)
  expect_identical(lengths(attr(id1, "hits")), 2L)
})

test_that(".edge_join() gives the same cells whichever layer it buffers", {
  # 925 points on 864 cells buffers the cells; 60 of them on the same cells
  # buffers the points.
  for (l in .ft_frames(.ft_lattice(37, 25))) {
    g <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(l)), n = c(36, 24))
    all <- .edge_join(l, g)
    sub <- c(1:30, 401:430)
    part <- .edge_join(l[sub, ], g, bb = c(sf::st_bbox(l), sf::st_bbox(g)))
    expect_identical(as.vector(part), as.vector(all)[sub])
    expect_identical(lapply(attr(part, "hits"), as.integer),
                     lapply(attr(all, "hits"), as.integer)[sub])
  }
})


# NNDM and kNNDM under the package's tie rule (BUG-102, BUG-103): the 20 x 15
# lattice above predicted onto a 30 x 25 lattice at 400 m that starts 3 km
# outside it, in metres, in UTM 15N US feet and with the false origin moved.
.ft_pred <- function() {
  g <- expand.grid(i = 0:29, j = 0:24)
  sf::st_as_sf(data.frame(x = 4.5e5 - 3000 + g$i * 400, y = 4.55e6 - 3000 + g$j * 400),
               coords = c("x", "y"), crs = .ft_crs_m())
}
.ft_train_sets <- function(r)
  vapply(r$folds, function(s) paste(sort(s$train), collapse = ","), "")
.ft_same_partition <- function(a, b) {
  t <- table(a, b)
  all(rowSums(t > 0) == 1) && all(colSums(t > 0) == 1)
}

test_that("nndm builds the same folds on a lattice in every frame", {
  skip_on_cran()
  F <- .ft_frames(.ft_lattice(20, 15)); P <- .ft_frames(.ft_pred())
  res <- lapply(names(F), function(f) .ft_quiet(make_folds(
    F[[f]], method = "nndm", prediction_points = P[[f]])))
  names(res) <- names(F)
  # Before, 195 (US feet) and 234 (origin moved) of the 300 training sets
  # differed from those in metres: 60670, 61078 and 60845 training rows.
  for (f in c("usft", "m01")) {
    expect_identical(.ft_train_sets(res[[f]]), .ft_train_sets(res$m), info = f)
    expect_identical(res[[f]]$params$n_removed_total, res$m$params$n_removed_total, info = f)
    expect_identical(res[[f]]$params$n_at_min_train, res$m$params$n_at_min_train, info = f)
    expect_equal(res[[f]]$params$realised_distances, .ft_len(F[[f]], res$m$params$realised_distances),
                 tolerance = 1e-9, info = f)
    expect_equal(res[[f]]$params$wasserstein, .ft_len(F[[f]], res$m$params$wasserstein),
                 tolerance = 1e-9, info = f)
  }
  expect_identical(sum(lengths(lapply(res$m$folds, `[[`, "train"))), 60670L)
})

test_that("nndm on a jittered lattice is unchanged in every frame", {
  skip_on_cran()
  F <- .ft_frames(.ft_lattice(20, 15, jitter = 60)); P <- .ft_frames(.ft_pred())
  res <- lapply(names(F), function(f) .ft_quiet(make_folds(
    F[[f]], method = "nndm", prediction_points = P[[f]])))
  for (r in res) expect_identical(sum(lengths(lapply(r$folds, `[[`, "train"))), 61090L)
  expect_identical(.ft_train_sets(res[[2]]), .ft_train_sets(res[[1]]))
  expect_identical(.ft_train_sets(res[[3]]), .ft_train_sets(res[[1]]))
})

test_that("nndm takes equidistant neighbours by the package's tie rule, with or without FNN", {
  F <- .ft_frames(.ft_lattice(12, 9))
  for (f in names(F)) {
    xy  <- sf::st_coordinates(F[[f]])
    key <- .nndm_tie_key(xy)
    a <- .nndm_neighbours(xy, 40L, key, use_fnn = TRUE)
    b <- .nndm_neighbours(xy, 40L, key, use_fnn = FALSE)
    expect_identical(a$idx, b$idx, info = f)
    expect_identical(a$d, b$d, info = f)
    # The package's own search, .knn_pick_rows() under .knn_exact_rows().
    expect_identical(a$idx, .knn_exact_rows(xy[, 1L], xy[, 2L], seq_len(nrow(xy)), 40L,
                                            tol = .TIE_TOL, num = key), info = f)
    if (f == "m") idx_m <- a$idx else expect_identical(a$idx, idx_m, info = f)
  }
  # The four nearest neighbours of an inner point go in the order of the
  # column, then the second coordinate: left, below, above, right.
  lat <- .ft_lattice(12, 9)
  xy  <- sf::st_coordinates(.ft_frames(lat)$m01)
  at  <- which(lat$i == 5 & lat$j == 4)
  nb  <- .nndm_neighbours(xy, 4L, .nndm_tie_key(xy))$idx[at, ]
  expect_identical(cbind(lat$i[nb], lat$j[nb]), cbind(c(4L, 5L, 5L, 6L), c(4L, 3L, 5L, 4L)))
  # A co-located twin is the first neighbour on either path.
  xy2 <- rbind(c(0, 0), c(0, 0), c(1, 0), c(0, 1), c(2, 2), c(3, 1))
  for (fnn in c(TRUE, FALSE)) {
    nb2 <- .nndm_neighbours(xy2, 3L, .nndm_tie_key(xy2), use_fnn = fnn)
    expect_identical(nb2$idx[1:2, 1L], c(2L, 1L))
    expect_identical(nb2$d[1:2, 1L], c(0, 0))
  }
})

test_that("nndm's tie key orders by column and row, last bits counting as equal", {
  xy <- rbind(c(1 + 1e-12, 1e-13), c(1, 1), c(0, 2), c(1 - 1e-12, 0), c(0, 1))
  # Column 1 holds rows 3 and 5 (x = 0), column 2 rows 1, 2 and 4; rows 1
  # and 4 tie on both coordinates and go by row, though row 4 has the lower
  # x and the lower y.
  expect_identical(.nndm_tie_key(xy), c(3L, 5L, 2L, 4L, 1L))
})

test_that("the ECDF deficit counts distances within the tie rule as one", {
  # Before: 1, from the step between 250 - 1e-10 and 250 + 1e-10.
  expect_identical(.ecdf_deficit(250 + 1e-10, 250 - 1e-10), 0)
  expect_identical(.ecdf_deficit(250 - 1e-10, 250 + 1e-10), 0)
  expect_identical(.ecdf_deficit(c(100, 300), c(200, 250)), 0.5)
  # Prediction points 250 or 500 m apart beside the lattice, so prediction
  # distances and held-out distances share lengths.  Before, with the same
  # nndm folds in the three frames, nndm's max_ecdf_excess was 0.42, 0.42
  # and 0.413 and its max_ecdf_deficit 0, 0 and 0.003 (250 m).
  F <- .ft_frames(.ft_lattice(20, 15))
  for (off in c(250, 500)) {
    g <- expand.grid(i = 1:20, j = 0:14)
    P <- .ft_frames(sf::st_as_sf(data.frame(x = 4.5e5 + 19 * 250 + g$i * off,
                                            y = 4.55e6 + g$j * 250),
                                 coords = c("x", "y"), crs = .ft_crs_m()))
    for (meth in c("nndm", "knndm")) {
      res <- lapply(names(F), function(f) .ft_quiet(make_folds(
        F[[f]], k = 5, method = meth, seed = 1, prediction_points = P[[f]])))
      for (nm in c("max_ecdf_excess", "max_ecdf_deficit", "loo_ecdf_deficit"))
        for (r in res[-1])
          expect_identical(r$params[[nm]], res[[1]]$params[[nm]], info = paste(meth, off, nm))
    }
  }
})

test_that("knndm builds the same folds and W on a lattice in every frame", {
  skip_on_cran()
  F <- .ft_frames(.ft_lattice(20, 15)); P <- .ft_frames(.ft_pred())
  res <- lapply(names(F), function(f) .ft_quiet(make_folds(
    F[[f]], k = 5, method = "knndm", seed = 1, prediction_points = P[[f]])))
  names(res) <- names(F)
  # Before, 170 (US feet) and 290 (origin moved) of the 300 rows were in
  # another fold, and W was 1191.21, 1223.84 and 1288.05 m.
  for (f in c("usft", "m01")) {
    expect_identical(res[[f]]$assignment$fold, res$m$assignment$fold, info = f)
    expect_identical(res[[f]]$params$chosen, res$m$params$chosen, info = f)
    expect_equal(res[[f]]$params$wasserstein, .ft_len(F[[f]], res$m$params$wasserstein),
                 tolerance = 1e-9, info = f)
    expect_equal(res[[f]]$params$candidates$W, .ft_len(F[[f]], res$m$params$candidates$W),
                 tolerance = 1e-9, info = f)
  }
})

test_that("knndm on a jittered lattice is unchanged in every frame", {
  skip_on_cran()
  F <- .ft_frames(.ft_lattice(20, 15, jitter = 60)); P <- .ft_frames(.ft_pred())
  res <- lapply(names(F), function(f) .ft_quiet(make_folds(
    F[[f]], k = 5, method = "knndm", seed = 1, prediction_points = P[[f]])))
  names(res) <- names(F)
  for (f in names(F)) {
    expect_equal(res[[f]]$params$wasserstein / .ft_len(F[[f]], 1), 1223.40, tolerance = 1e-5,
                 info = f)
    expect_identical(res[[f]]$assignment$fold, res$m$assignment$fold, info = f)
  }
})

test_that("knndm's k-means and its kmeans+ward path give the same clusters in every frame", {
  skip_on_cran()
  F <- .ft_frames(.ft_lattice(40, 30))
  cl <- lapply(F, function(l) {
    z <- .kmeans_frame(unname(sf::st_coordinates(l)))$z
    lapply(c(TRUE, FALSE), function(fnn) withr::with_seed(7, .knndm_kmeans(
      z, z[!duplicated(z), ], 150L, use_fnn = fnn)))
  })
  for (f in names(F)) {
    # With and without FNN: the same nearest centre, ties to the lower one.
    expect_identical(cl[[f]][[1]], cl[[f]][[2]], info = f)
    expect_identical(cl[[f]][[1]], cl$m[[1]], info = f)
  }
  # Ward's method continued from k-means groups.  The target is the
  # distance from a coarser lattice of prediction points.
  pr <- sf::st_coordinates(.ft_frames(.ft_pred())$m)
  res <- lapply(F, function(l) {
    xy <- unname(sf::st_coordinates(l))
    gt <- FNN::get.knnx(xy, if (isTRUE(sf::st_crs(l) == .ft_crs_m())) pr
                            else sf::st_coordinates(sf::st_transform(
                              .ft_pred(), sf::st_crs(l))), k = 1)$nn.dist[, 1]
    withr::with_seed(3, .knndm_folds(xy, gt, 5L, ward_max = 300L))
  })
  for (f in names(F)) {
    expect_identical(res[[f]]$clustering, "kmeans+ward")
    expect_identical(res[[f]]$fold, res$m$fold, info = f)
    expect_identical(res[[f]]$chosen, res$m$chosen, info = f)
  }
})

test_that("nndm compares held-out distances with phi under the tie rule", {
  skip_on_cran()
  # phi at a length of the lattice: held-out distances of 500 m or of the
  # diagonal are at phi but for their last bits.  Before, with phi = 500 m,
  # 239 of the 300 training sets in US feet differed from those in metres,
  # and n_removed_total was 2675, 1812 and 2675.
  F <- .ft_frames(.ft_lattice(20, 15)); P <- .ft_frames(.ft_pred())
  for (phi_m in c(500, 250 * sqrt(2))) {
    res <- lapply(names(F), function(f) .ft_quiet(make_folds(
      F[[f]], method = "nndm", prediction_points = P[[f]], phi = .ft_len(F[[f]], phi_m))))
    names(res) <- names(F)
    for (f in c("usft", "m01")) {
      expect_identical(.ft_train_sets(res[[f]]), .ft_train_sets(res$m), info = paste(phi_m, f))
      expect_identical(res[[f]]$params$n_removed_total, res$m$params$n_removed_total,
                       info = paste(phi_m, f))
    }
    expect_equal(res$m$params$n_removed_total, if (phi_m == 500) 2675 else 1812,
                 info = phi_m)
  }
})

test_that("knndm keeps the earlier of candidates with the same W in every frame", {
  skip_on_cran()
  # An 8 x 8 lattice predicted onto a 12 x 12 lattice at 250 / 3 m from the
  # same corner: many candidates share a W.  Before, 45 of the 64 rows in US
  # feet and with the origin moved were in another fold than in metres, with
  # 38, 63 and 63 clusters.
  g <- expand.grid(i = 0:11, j = 0:11)
  P <- .ft_frames(sf::st_as_sf(data.frame(x = 4.5e5 + g$i * 250 / 3, y = 4.55e6 + g$j * 250 / 3),
                               coords = c("x", "y"), crs = .ft_crs_m()))
  F <- .ft_frames(.ft_lattice(8, 8))
  res <- lapply(names(F), function(f) .ft_quiet(make_folds(
    F[[f]], k = 4, method = "knndm", seed = 1, prediction_points = P[[f]])))
  names(res) <- names(F)
  W <- res$m$params$candidates$W
  expect_true(any(duplicated(signif(W[is.finite(W)], 12))))
  for (f in c("usft", "m01")) {
    expect_identical(res[[f]]$assignment$fold, res$m$assignment$fold, info = f)
    expect_identical(res[[f]]$params$n_clusters, res$m$params$n_clusters, info = f)
    expect_identical(res[[f]]$params$chosen, res$m$params$chosen, info = f)
  }
})

test_that("knndm takes the earliest candidate within the tie rule of the least W", {
  # W values 1, 1 - 0.6e-7 and 1 - 1.2e-7 in turn: the second and third are
  # within one part in ten million of each other and of the least, so the
  # second is returned.  Replacing the best only when lower by more than the
  # tolerance took the third.
  w <- c(1, 1 - 0.6e-7, 1 - 1.2e-7)
  i <- 0L
  local_mocked_bindings(.fold_sep_compare = function(a, b, ...) {
    i <<- i + 1L
    list(w1 = if (i <= length(w)) w[i] else 2)
  })
  xy <- as.matrix(expand.grid(x = 0:5, y = 0:5)) * 250
  res <- withr::with_seed(1, .knndm_folds(xy, rep(300, 20), 2L, use_fnn = FALSE))
  sc  <- which(is.finite(res$candidates$W))
  expect_gte(length(sc), 3L)
  expect_identical(res$candidates$W[sc[1:3]], w)
  expect_identical(res$chosen, sc[2])
  # An exact tie goes to the earlier candidate, the random partition first.
  w <- c(1, 1, 1 - 0.6e-7); i <- 0L
  res <- withr::with_seed(1, .knndm_folds(xy, rep(300, 20), 2L, use_fnn = FALSE))
  expect_identical(res$chosen, 1L)
  expect_identical(res$fold, withr::with_seed(1, .random_fold_labels(36L, 2L)))
})

test_that("nndm's complete separation counts distances within the tie rule as equal", {
  # The largest prediction distance and the least leave-one-out distance
  # differ in their last bits only.  Before, separated was TRUE in one order
  # and FALSE in the other, so it moved with the frame.
  a <- .nndm_reach(c(250 + 1e-10, 300, 400), c(100, 250 - 1e-10))
  b <- .nndm_reach(c(250 - 1e-10, 300, 400), c(100, 250 + 1e-10))
  expect_false(a$separated)
  expect_false(b$separated)
  expect_identical(a$deficit, b$deficit)
  expect_true(.nndm_reach(c(251, 300, 400), c(100, 250))$separated)
  expect_identical(.nndm_reach(c(251, 300, 400), c(100, 250))$deficit, 1)
})
