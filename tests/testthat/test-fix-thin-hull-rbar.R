# tests/testthat/test-fix-thin-hull-rbar.R
# ---------------------------------------------------------------------------
# .rbar_domain() reads the domain term of resolution_profile()'s reliability,
# the mean correlation between two points of the hull, by quadrature on a
# grid.  The grid lay over the hull's bounding box at a spacing s that
# follows the hull's area, not its width.  On a strip along an axis
# narrower than s / 2, seq() started past its end and resolution_profile()
# with a fitted variogram stopped with "wrong sign in 'by' argument"; on a
# thin strip at a slant the rows crossed part of its length, or fewer than
# 10 nodes fell inside and the bounding box's value was used.  A hull
# narrower than 2 s is now laid out in the frame of its rectangle of least
# area, with its nodes spread over the width of each band; every other hull
# keeps the box grid to the bit.
#
# The layer: 300 points over a strip of 20000 x 10 m in EPSG:32615
# (seed 5), rotated about its corner.
# ---------------------------------------------------------------------------

.thr_strip_xy <- function(rot = 0, width = 10, n = 300L, len = 20000) {
  withr::local_seed(5, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
  x <- stats::runif(n, 0, len); y <- stats::runif(n, 0, width)
  z <- sin(x / 800) + stats::rnorm(n, 0, 0.3)
  th <- rot * pi / 180
  list(xy = cbind(5e5 + x * cos(th) - y * sin(th), 4e6 + x * sin(th) + y * cos(th)),
       z = z, len = diff(range(x)))
}

.thr_hull <- function(xy, crs = 32615) {
  sf::st_convex_hull(sf::st_union(sf::st_sfc(sf::st_multipoint(xy), crs = crs)))
}

# The box grid as .rbar_domain() laid it before the fix, for the layers that
# keep it.  .thr_frame_grid() below is the decided layout of a thin hull.
.thr_box_grid <- function(cor_fn, domain, area, q = 20L) {
  bb <- sf::st_bbox(domain)
  w  <- as.numeric(bb[["xmax"]] - bb[["xmin"]]); h <- as.numeric(bb[["ymax"]] - bb[["ymin"]])
  if (!is.finite(area) || area <= 0 || !is.finite(w) || !is.finite(h) || w <= 0 || h <= 0)
    return(.rbar_rect(cor_fn, w, h, q))
  s  <- max(sqrt(area) / q, sqrt(w * h / (40 * q^2)))
  gx <- seq(bb[["xmin"]] + s / 2, bb[["xmax"]], by = s)
  gy <- seq(bb[["ymin"]] + s / 2, bb[["ymax"]], by = s)
  g  <- as.matrix(expand.grid(x = gx, y = gy))
  inside <- lengths(sf::st_intersects(
    sf::st_as_sf(as.data.frame(g), coords = c("x", "y"), crs = sf::st_crs(domain)),
    domain)) > 0L
  if (sum(inside) < 10L) return(.rbar_rect(cor_fn, w, h, q))
  r <- as.numeric(cor_fn(as.numeric(stats::dist(g[inside, , drop = FALSE]))))
  mean(r[is.finite(r)])
}

# The layout of a hull narrower than 2 s, written out from the decision: in
# the frame of the rectangle of least area, spacing t = sqrt(L W) / q,
# k = max(1, floor(W / t)) equal bands across, and in each band a rank-1
# lattice of round(q^2 / (4 k)) nodes over a quarter of it with its mirror
# images in the band's two midlines.
.thr_frame_grid <- function(cor_fn, domain, q = 20L) {
  r  <- .least_area_rect(domain)
  t  <- sqrt(r$L * r$W) / q
  k  <- max(1L, as.integer(floor(r$W / t)))
  nq <- max(1L, as.integer(round(q^2 / (4 * k))))
  j  <- seq_len(nq)
  u1 <- (j - 0.5) * r$L / (2 * nq)
  v1 <- (((j - 0.5) * ((sqrt(5) - 1) / 2)) %% 1) / 2
  uq <- c(u1, r$L - u1, u1, r$L - u1)
  vq <- c(v1, v1, 1 - v1, 1 - v1)
  f  <- cbind(r$u0 + rep(uq, times = k),
              r$v0 + r$W * (rep(seq_len(k) - 1L, each = 4L * nq) + rep(vq, times = k)) / k)
  g  <- cbind(x = r$centre[1] + f[, 1] * r$a[1] + f[, 2] * r$b[1],
              y = r$centre[2] + f[, 1] * r$a[2] + f[, 2] * r$b[2])
  inside <- lengths(sf::st_intersects(
    sf::st_as_sf(as.data.frame(g), coords = c("x", "y"), crs = sf::st_crs(domain)),
    domain)) > 0L
  r <- as.numeric(cor_fn(as.numeric(stats::dist(g[inside, , drop = FALSE]))))
  list(value = mean(r[is.finite(r)]), k = k, n = nrow(g))
}

.thr_rot <- function(xy, rot) {
  th <- rot * pi / 180
  cbind(5e5 + xy[, 1] * cos(th) - xy[, 2] * sin(th), 4e6 + xy[, 1] * sin(th) + xy[, 2] * cos(th))
}

test_that(".rbar_domain() reads a thin strip over its length at every angle", {
  a  <- 1000
  cf <- function(h) exp(-h / a)
  vals <- c()
  for (rot in c(0, 0.01, 5, 30, 45)) {
    s <- .thr_strip_xy(rot)
    hull <- .thr_hull(s$xy)
    v <- .rbar_domain(cf, hull, as.numeric(sf::st_area(hull)))
    expect_true(is.finite(v))
    # Mean correlation of two points on a segment of length L.  The strip's
    # width is 1 percent of the range, so its value is the segment's.  The
    # quadrature is the mean over the off-diagonal pairs of about 400 nodes
    # along the strip, which sits below the continuous mean by about
    # (1 - rbar) / 400, 2.3 percent of it here, and the hull's tapered ends
    # drop some of the nodes spread over its width, which moves it back by 1
    # to 2 percent; 3 percent holds both on every platform.  Before the fix the values were an error
    # at 0 degrees and 64, 5, 39 and 7 percent off at the others.
    L  <- s$len
    cl <- 2 * a / L - 2 * (a / L)^2 * (1 - exp(-L / a))
    expect_equal(v, cl, tolerance = 0.03, label = sprintf("rbar at %g degrees", rot))
    vals <- c(vals, v)
  }
  # The layout is the hull's own frame, so a rotation leaves the value as it
  # is up to the rounding of the rotated coordinates and of GEOS's hull.
  expect_lt(max(vals) / min(vals) - 1, 1e-6)
})

test_that(".least_area_rect() finds the rectangle of least area over the hull's edges", {
  # Every edge direction against every vertex, in coordinates about the
  # hull's centre, as a reference for the calipers.
  brute <- function(hull) {
    xy <- hull[[1L]][[1L]][, 1:2]
    xy <- sweep(xy, 2L, colMeans(xy))
    d <- diff(xy); d <- d[rowSums(d^2) > 0, , drop = FALSE]
    e <- d / sqrt(rowSums(d^2))
    pu <- xy %*% t(e); pv <- xy %*% t(cbind(-e[, 2], e[, 1]))
    su <- apply(pu, 2L, max) - apply(pu, 2L, min)
    sv <- apply(pv, 2L, max) - apply(pv, 2L, min)
    j <- which.min(su * sv)
    c(min(su[j], sv[j]), max(su[j], sv[j]))
  }
  withr::local_seed(41, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
  for (r in seq_len(40L)) {
    n <- c(3L, 6L, 40L, 400L)[(r - 1L) %% 4L + 1L]
    th <- stats::runif(1, 0, 2 * pi)
    sx <- 10^stats::runif(1, 0, 4); sy <- 10^stats::runif(1, 0, 4)
    if (r %% 2L == 0L) {             # every point on the hull
      a <- stats::runif(n, 0, 2 * pi); x <- cos(a) * sx; y <- sin(a) * sy
    } else {
      x <- stats::runif(n) * sx; y <- stats::runif(n) * sy
    }
    xy <- cbind(1e6 + x * cos(th) - y * sin(th), 2e6 + x * sin(th) + y * cos(th))
    hull <- .thr_hull(xy)
    got <- .least_area_rect(hull)
    ref <- brute(hull)
    # The rectangle of least area is unique in area; its sides are the
    # same up to rounding in the projections, a few ulps of the hull's size.
    expect_equal(got$W * got$L, ref[1] * ref[2], tolerance = 1e-9)
    # The frame holds every vertex.
    v <- sweep(hull[[1L]][[1L]][, 1:2], 2L, got$centre)
    pu <- v %*% got$a - got$u0; pv <- v %*% got$b - got$v0
    tol <- 1e-9 * got$L
    expect_true(all(pu >= -tol & pu <= got$L + tol & pv >= -tol & pv <= got$W + tol))
  }
  # Points on a line have no rectangle with area.
  line <- .thr_hull(cbind(5e5 + 0:10 * 200, 4e6 + 0:10 * 100))
  expect_true(is.na(.least_area_rect(line)$W))
})

test_that("resolution_profile() with a fitted variogram runs on a thin strip in m and in km", {
  s <- .thr_strip_xy(0)
  pts <- sf::st_as_sf(data.frame(x = s$xy[, 1], y = s$xy[, 2], z = s$z),
                      coords = c("x", "y"), crs = 32615)
  sac <- suppressWarnings(estimate_sac_range(pts, "z"))
  rp <- suppressWarnings(resolution_profile(pts, "z", sac = sac,
                                            levels = c(2, 4, 8, 16, 32), seed = 1L))
  expect_true(all(is.finite(rp$reliability)))

  # The same strip with its coordinates in units of 1 km.
  g <- sf::st_geometry(pts) * (1 / 1000)
  km <- sf::st_sf(z = pts$z, geometry = g, crs = sf::st_crs(3857))
  sac_km <- suppressWarnings(estimate_sac_range(km, "z"))
  rp_km <- suppressWarnings(resolution_profile(km, "z", sac = sac_km,
                                               levels = c(2, 4, 8, 16, 32), seed = 1L))
  expect_true(all(is.finite(rp_km$reliability)))
})

test_that("a hull two or more grid rows wide keeps the box grid to the bit", {
  cf <- function(h) exp(-h / 300)
  mk <- function(w, h, rot, n = 600L) {
    withr::local_seed(11, .rng_kind = "Mersenne-Twister",
                      .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
    x <- stats::runif(n, 0, w); y <- stats::runif(n, 0, h)
    th <- rot * pi / 180
    .thr_hull(cbind(5e5 + x * cos(th) - y * sin(th), 4e6 + x * sin(th) + y * cos(th)))
  }
  # The layers .rbar_domain()'s help names: a square of 1000 and a strip of
  # 3000 x 120, along an axis and at 45 degrees.  Same arithmetic in one
  # process, so identical() holds on every platform.
  for (hull in list(mk(1000, 1000, 0), mk(3000, 120, 0), mk(3000, 120, 45))) {
    area <- as.numeric(sf::st_area(hull))
    expect_identical(.rbar_domain(cf, hull, area), .thr_box_grid(cf, hull, area))
  }
})

test_that(".rbar_domain() reads a thin hull that tapers over its area", {
  # A wedge, a rhombus and a lens of 20000 x 10 m.  Their width is 1 percent
  # of the range, so their mean correlation is that of two points on the
  # long axis drawn with density in proportion to the width there, w(x):
  # the reference below, by the midpoint rule on 4000 steps.  Along a centre
  # line the wedge and the rhombus read as a segment of half their length
  # (43 and 35 percent high) and the lens as one of its whole length (12
  # percent low).  The nodes spread over the band read them 0.5 to 4.0
  # percent low: the mean over the off-diagonal pairs of the n nodes inside
  # (192 to 320) sits below the continuous mean by about (1 - rbar) / (n - 1),
  # 2.6 to 3.6 percent of it, and the lattice's own error moves that by up
  # to 3 percent either way.  6 percent holds that on every platform and is
  # far from the centre line's reading.
  a <- 1000; L <- 20000; W <- 10
  cf <- function(h) exp(-h / a)
  ref <- function(wfun, n = 4000L) {
    x <- (seq_len(n) - 0.5) * L / n
    p <- wfun(x); p <- p / sum(p)
    sum(outer(p, p) * cf(abs(outer(x, x, "-"))))
  }
  tt <- seq(0, 2 * pi, length.out = 721L)[-721L]
  shapes <- list(
    wedge   = list(xy = rbind(c(0, 0), c(L, 0), c(0, W)),
                   w = function(x) 1 - x / L),
    rhombus = list(xy = rbind(c(0, W / 2), c(L / 2, 0), c(L, W / 2), c(L / 2, W)),
                   w = function(x) 1 - abs(2 * x / L - 1)),
    lens    = list(xy = cbind(L / 2 + L / 2 * cos(tt), W / 2 + W / 2 * sin(tt)),
                   w = function(x) sqrt(pmax(0, 1 - (2 * x / L - 1)^2))))
  for (nm in names(shapes)) {
    rf <- ref(shapes[[nm]]$w)
    vals <- c()
    for (rot in c(0, 20, 45)) {
      hull <- .thr_hull(.thr_rot(shapes[[nm]]$xy, rot))
      v <- .rbar_domain(cf, hull, as.numeric(sf::st_area(hull)))
      expect_equal(v, rf, tolerance = 0.06, label = sprintf("%s at %g degrees", nm, rot))
      vals <- c(vals, v)
    }
    # The mirrors make the nodes the same whichever way the frame's axes
    # point, so the angle moves the value by rounding only.
    expect_lt(max(vals) / min(vals) - 1, 1e-6, label = nm)
  }
})

test_that("the box grid is kept from a width of 2 s and the frame is used below it", {
  cf <- function(h) exp(-h / 300)
  rect <- function(L, W, rot = 0) .thr_hull(.thr_rot(rbind(c(0, 0), c(L, 0), c(L, W), c(0, W)), rot))
  # 3000 x 29 along an axis: s = sqrt(87000) / 20 = 14.75, 2 s = 29.50, so
  # the width is under 2 s and the hull takes the frame layout.  Its area
  # over the box's diagonal (29.00) does not settle it, so the rectangle is
  # computed.
  h29 <- rect(3000, 29); a29 <- as.numeric(sf::st_area(h29))
  v29 <- .rbar_domain(cf, h29, a29)
  expect_identical(v29, .thr_frame_grid(cf, h29)$value)
  expect_false(identical(v29, .thr_box_grid(cf, h29, a29)))
  # 3000 x 31: s = 15.25, 2 s = 30.50; the area over the diagonal (31.00)
  # settles the width without the rectangle, and the box grid is kept.
  h31 <- rect(3000, 31); a31 <- as.numeric(sf::st_area(h31))
  expect_identical(.rbar_domain(cf, h31, a31), .thr_box_grid(cf, h31, a31))
  # A rhombus 3000 x 40: s = 12.25, 2 s = 24.49; its area over the diagonal
  # (20.00) leaves the width open, the rectangle of least area gives 40.00,
  # at least 2 s and under 3.5 s, and the box grid is kept.
  hr <- .thr_hull(.thr_rot(rbind(c(0, 20), c(1500, 0), c(3000, 20), c(1500, 40)), 0))
  ar <- as.numeric(sf::st_area(hr))
  expect_identical(.rbar_domain(cf, hr, ar), .thr_box_grid(cf, hr, ar))
})

test_that("a hull wide enough for two bands of the frame gets two", {
  # 3000 x 32 at 45 degrees: the cap on the box's nodes raises s to 16.95,
  # 2 s = 33.90 > 32, so the frame layout; t = sqrt(3000 * 32) / 20 = 15.49
  # and k = floor(32 / 15.49) = 2 bands of round(400 / 8) * 4 = 200 nodes.
  cf <- function(h) exp(-h / 300)
  hull <- .thr_hull(.thr_rot(rbind(c(0, 0), c(3000, 0), c(3000, 32), c(0, 32)), 45))
  fr <- .thr_frame_grid(cf, hull)
  expect_identical(fr$k, 2L)
  expect_identical(fr$n, 400L)
  expect_identical(.rbar_domain(cf, hull, as.numeric(sf::st_area(hull))), fr$value)
})

test_that(".least_area_rect() returns the long side as L and an orthonormal frame", {
  withr::local_seed(43, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
  hulls <- list(
    # Four corners along an axis: GEOS starts the ring on a short side.
    .thr_hull(.thr_rot(rbind(c(0, 0), c(3000, 0), c(3000, 20), c(0, 20)), 0)),
    .thr_hull(.thr_rot(rbind(c(0, 0), c(20, 0), c(20, 3000), c(0, 3000)), 0)),
    .thr_hull(.thr_rot(rbind(c(0, 0), c(3000, 0), c(3000, 20), c(0, 20)), 30)))
  for (r in seq_len(12L)) {
    hulls[[length(hulls) + 1L]] <- .thr_hull(.thr_rot(
      cbind(stats::runif(30) * 10^stats::runif(1, 0, 4), stats::runif(30) * 10^stats::runif(1, 0, 4)),
      stats::runif(1, 0, 360)))
  }
  for (hull in hulls) {
    got <- .least_area_rect(hull)
    expect_gte(got$L, got$W)
    # Unit vectors at right angles, to the rounding of a normalised edge.
    expect_equal(sum(got$a^2), 1, tolerance = 1e-12)
    expect_equal(sum(got$b^2), 1, tolerance = 1e-12)
    expect_lt(abs(sum(got$a * got$b)), 1e-12)
  }
  r0 <- .least_area_rect(hulls[[1L]])
  expect_equal(c(r0$L, r0$W), c(3000, 20), tolerance = 1e-12)
})
