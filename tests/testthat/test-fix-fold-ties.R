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
