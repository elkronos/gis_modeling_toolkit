# tests/testthat/test-fix-sample-draw-size.R
# ---------------------------------------------------------------------------
# The draw sf makes to sample a boundary, bounded.
#
# sf::st_sample() draws n times the ratio of the bounding box to the area of
# points in the box and keeps those inside, at about 470 bytes per point.
# The package passed n on as given, so random and k-means seeding on a
# projected boundary of small pieces far apart, or a thin strip at a slant,
# stopped with R's "vector memory limit ... reached" (the k-means cloud of
# two 0.05 degree squares on two continents in EPSG:3857 asked for 2.8e9
# points).  The draw is now read first (.sample_draw_plan()): at most 10^6
# points (or 10 n) it is made whole, as before; above that the boundary is
# sampled piece by piece, a piece above the cap in draws of at most 10^6
# points in its box, and where the points needed in all exceed 2e7 (or
# 10 n) the call stops before any draw with an error that names the cause
# and the remedy.
#
# Every call runs with mem.maxVSize() lowered to 4096 MB, so a regression
# stops on that limit instead of taking the machine's memory, and the plan
# is checked on its own first, so a regression fails there before any draw.
# ---------------------------------------------------------------------------

quiet <- function(expr) suppressWarnings(suppressMessages(expr))
local_vsize <- function(mb = 4096, env = parent.frame()) {
  old <- mem.maxVSize()
  # R does not set the limit below the vector heap it holds, and says so in
  # a warning; after a gc() the heap is what the session uses.  Where it is
  # still larger, the test runs under the limit the session had.
  invisible(gc())
  withCallingHandlers(mem.maxVSize(mb), warning = function(w) invokeRestart("muffleWarning"))
  withr::defer(mem.maxVSize(old), envir = env)
}
sq <- function(x0, y0, d)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + d, y0), c(x0 + d, y0 + d), c(x0, y0 + d), c(x0, y0))))
layer <- function(..., crs = 4326) sf::st_sf(geometry = sf::st_sfc(list(...), crs = crs))
# A strip 50 km long and w metres wide on the diagonal in UTM 17N; its box
# is about 1.25e9 / (5e4 w) times its area.
strip_at <- function(w, o = c(5e5, 4.5e6)) {
  L <- 5e4 / sqrt(2); dw <- w / sqrt(2)
  sf::st_polygon(list(rbind(o, o + c(L, L), o + c(L, L) + c(dw, -dw), o + c(dw, -dw), o)))
}
strip <- function(w) layer(strip_at(w), crs = 32617)
inside <- function(pts, b) all(lengths(sf::st_intersects(pts, sf::st_union(b))) > 0L)
plan_of <- function(b, n) spatialkit:::.sample_draw_plan(sf::st_union(b), n)
# Every call of sf::st_sample() a seeder makes, with the points it draws in
# the box: the size asked for times the box-to-area ratio, as sf rounds it.
local_draw_log <- function(env = parent.frame()) {
  log <- new.env()
  log$box <- numeric(0)
  orig <- spatialkit:::.st_sample_quiet
  testthat::local_mocked_bindings(
    .st_sample_quiet = function(x, size, ...) {
      r <- spatialkit:::.sample_draw_plan(x, size)
      log$box <- c(log$box, round(r$draws))
      orig(x, size = size, ...)
    },
    .package = "spatialkit", .env = env)
  log
}

# Two 0.05 degree squares on two continents, projected; and three spread over
# 100 degrees of longitude.
far2_ll <- layer(sq(10, 45, 0.05), sq(120, -30, 0.05))
far2_3857  <- sf::st_transform(far2_ll, 3857)
far2_32633 <- sf::st_transform(far2_ll, 32633)
far3_3857 <- sf::st_transform(layer(sq(-50, 10, 0.1), sq(0, 20, 0.05), sq(50, 30, 0.05)), 3857)
# A diagonal strip 5 degrees long and 0.0005 wide, its box 5000 times its area.
strip_ll <- layer(sf::st_polygon(list(rbind(c(-90, 30), c(-85, 35), c(-84.9995, 34.9995),
                                            c(-89.9995, 29.9995), c(-90, 30)))))
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)

test_that("the draw plan reads sf's draw without touching the stream", {
  local_vsize()
  set.seed(1)
  before <- .Random.seed
  p <- plan_of(far2_3857, 2000)
  expect_identical(.Random.seed, before)
  expect_false(p$longlat)
  expect_gt(p$ratio, 1e6)
  expect_gt(p$draws, spatialkit:::.sample_draw_cap(2000))
  # In degrees the ratio is sf's own: sf's refusal prints the same number.
  pl <- plan_of(far2_ll, 10)
  expect_true(pl$longlat)
  expect_true(pl$sf_refuses)
  said <- tryCatch(quiet(sf::st_sample(sf::st_union(far2_ll), 10)),
                   error = function(e) conditionMessage(e))
  expect_match(said, sprintf("sampling box is %s times", format(pl$ratio)), fixed = TRUE)
  expect_false(plan_of(strip_ll, 10)$sf_refuses)
  expect_equal(plan_of(strip_ll, 10)$ratio, 5000, tolerance = 0.01)
  # The strips, and boundaries that fill their boxes.
  expect_equal(plan_of(strip(0.05), 1)$ratio, 5e5, tolerance = 1e-3)
  expect_equal(plan_of(strip(20), 1)$ratio, 1250, tolerance = 1e-3)
  nc_pr <- sf::st_transform(nc, 32617)
  for (n in c(10, 200, 10000))
    expect_lte(plan_of(nc_pr, n)$draws, spatialkit:::.sample_draw_cap(n))
  # Per piece, with `each`.
  pe <- spatialkit:::.sample_draw_plan(sf::st_geometry(far2_3857), 1, each = TRUE)
  expect_length(pe$ratio, 2L)
  expect_true(all(pe$ratio < 2))
})

test_that("projected pieces far apart are sampled piece by piece", {
  local_vsize()
  for (b in list(far2_3857, far2_32633)) {
    lab <- sf::st_crs(b)$input
    k <- quiet(get_voronoi_seeds(b, n = 10, method = "kmeans", set_seed = 1))
    expect_identical(nrow(k), 10L, label = lab)
    expect_true(inside(k, b), label = lab)
    expect_identical(quiet(get_voronoi_seeds(b, n = 10, method = "kmeans", set_seed = 1)), k,
                     label = lab)
    r <- quiet(get_voronoi_seeds(b, n = 2000, method = "random", set_seed = 1))
    expect_identical(nrow(r), 2000L, label = lab)
    expect_true(inside(r, b), label = lab)
    expect_identical(quiet(get_voronoi_seeds(b, n = 2000, method = "random", set_seed = 1)), r,
                     label = lab)
    v <- quiet(voronoi_seeds_random(b, 2000, set_seed = 1))
    expect_identical(nrow(v), 2000L, label = lab)
    expect_true(inside(v, b), label = lab)
    expect_identical(quiet(voronoi_seeds_random(b, 2000, set_seed = 1)), v, label = lab)
  }
})

test_that("projected pieces get seeds in proportion to their planar areas", {
  local_vsize()
  s <- quiet(voronoi_seeds_random(far3_3857, 3000, set_seed = 1))
  piece <- vapply(sf::st_intersects(s, far3_3857), `[`, integer(1), 1L)
  n_in <- tabulate(piece, 3L)
  a <- as.numeric(sf::st_area(far3_3857))
  expect_identical(sum(n_in), 3000L)
  # A chi-squared test of the counts against the shares of the areas, at a
  # fixed seed; the counts are random, so a p-value is the measure of fit.
  expect_gt(stats::chisq.test(n_in, p = a / sum(a))$p.value, 0.01)
})

test_that("a thin strip is drawn in bounded draws", {
  local_vsize()
  b <- strip(20)
  log <- local_draw_log()
  k <- quiet(get_voronoi_seeds(b, n = 10, method = "kmeans", set_seed = 1))
  expect_identical(nrow(k), 10L)
  expect_true(inside(k, b))
  # The cloud of 2000 points takes about 2.5e6 points in the box, in draws of
  # at most 10^6.
  expect_gte(length(log$box), 3L)
  expect_lte(max(log$box), 1e6)
  # A strip whose box is 2e6 times its area: one point asks for a fraction
  # of a point per draw, which sf rounds to 10^6 points in the box.
  thin <- strip(0.0125)
  log$box <- numeric(0)
  v <- quiet(voronoi_seeds_random(thin, 1, set_seed = 1))
  expect_identical(nrow(v), 1L)
  expect_true(inside(v, thin))
  expect_lte(max(log$box), 1e6)
})

test_that("a draw above the limit stops before it starts, with the cause and the remedy", {
  local_vsize()
  log <- local_draw_log()
  for (w in c(0.05, 1)) {
    e <- tryCatch(quiet(get_voronoi_seeds(strip(w), n = 10, method = "kmeans", set_seed = 1)),
                  error = function(e) e)
    expect_s3_class(e, "error")
    expect_match(conditionMessage(e), paste0(
      "^get_voronoi_seeds\\(\\): `boundary` cannot be sampled: the bounding box of ",
      "`boundary` is [0-9,]+ times its area"))
    expect_match(conditionMessage(e), "Pass `sample_points`", fixed = TRUE)
    expect_no_match(conditionMessage(e), "memory", fixed = TRUE)
  }
  expect_match(conditionMessage(e), "is 25,000 times its area", fixed = TRUE)
  expect_match(conditionMessage(e), "drawing 2,000 points", fixed = TRUE)
  expect_match(conditionMessage(e), "about 50,000,000 points drawn in the box", fixed = TRUE)
  # Random seeds are told to ask for fewer, in the argument's own name.
  e <- tryCatch(quiet(voronoi_seeds_random(strip(0.05), 100, set_seed = 1)),
                error = function(e) conditionMessage(e))
  expect_match(e, "^voronoi_seeds_random\\(\\): `boundary` cannot be sampled")
  expect_match(e, "Ask for fewer seeds (`k`)", fixed = TRUE)
  e <- tryCatch(quiet(get_voronoi_seeds(strip(0.05), n = 100, method = "random", set_seed = 1)),
                error = function(e) conditionMessage(e))
  expect_match(e, "Ask for fewer seeds (`n`)", fixed = TRUE)
  # Two strips far apart: the polygon is named.
  two <- layer(strip_at(0.05), strip_at(0.05, o = c(3e5, 3e6)), crs = 32617)
  e <- tryCatch(quiet(voronoi_seeds_random(two, 100, set_seed = 1)),
                error = function(e) conditionMessage(e))
  expect_match(e, "the bounding box of one of its polygons is 500,000 times that polygon's area",
               fixed = TRUE)
  # The lon/lat strip: its cloud for 200 seeds takes 5e7 points in the box.
  e <- tryCatch(quiet(get_voronoi_seeds(strip_ll, n = 200, method = "kmeans", set_seed = 1)),
                error = function(e) conditionMessage(e))
  expect_match(e, "is 5,000 times its area", fixed = TRUE)
  expect_match(e, "Pass `sample_points`", fixed = TRUE)
  # Nothing was drawn.
  expect_length(log$box, 0L)
})

test_that("`sample_points` inside a thin strip is clustered as given", {
  local_vsize()
  b <- strip(0.05)
  set.seed(2)
  t <- stats::runif(100, 0.01, 0.99)
  xy <- cbind(5e5 + t * 5e4 / sqrt(2) + 0.01, 4.5e6 + t * 5e4 / sqrt(2) - 0.005)
  sp <- sf::st_sf(geometry = sf::st_sfc(lapply(seq_len(100), function(i) sf::st_point(xy[i, ])),
                                        crs = 32617))
  expect_true(inside(sp, b))
  k <- quiet(get_voronoi_seeds(b, n = 10, method = "kmeans", sample_points = sp, set_seed = 1))
  expect_identical(nrow(k), 10L)
})

test_that("a boundary whose draw stays under the cap is drawn whole, as before", {
  local_vsize()
  testthat::local_mocked_bindings(
    .sample_by_piece = function(...) stop("the piece-by-piece path ran"),
    .package = "spatialkit")
  near_pr <- sf::st_transform(layer(sq(10, 45, 0.05), sq(11, 45, 0.05)), 3857)
  near_ll <- layer(sq(10, 45, 0.05), sq(11, 45, 0.05))
  nc_pr <- sf::st_transform(nc, 32617)
  nc_ll <- sf::st_transform(nc, 4326)
  for (b in list(near_pr, near_ll, nc_pr, nc_ll, strip(20), strip_ll)) {
    u <- spatialkit:::.with_s2(sf::st_union(b))
    for (n in c(10, 200)) {
      lab <- sprintf("%s, n = %d", sf::st_crs(b)$input, n)
      # The draw is sf's whole draw, to the bit.
      set.seed(3)
      got <- spatialkit:::.with_s2(spatialkit:::.robust_st_sample(u, n, "x", "nothing"))
      set.seed(3)
      want <- spatialkit:::.with_s2(quiet(sf::st_sample(u, size = n, type = "random",
                                                        exact = TRUE)))
      expect_identical(got, want, label = lab)
      expect_identical(nrow(quiet(voronoi_seeds_random(b, n, set_seed = 1))), as.integer(n),
                       label = lab)
    }
  }
  for (b in list(near_pr, nc_pr, nc_ll))
    for (n in c(10, 200))
      expect_identical(nrow(quiet(get_voronoi_seeds(b, n = n, method = "kmeans",
                                                    set_seed = 1))), as.integer(n))
})

# The caps and the choices made on them, read on plans built by hand, so no
# point is drawn: a wrong cap, a dropped 10 n term or a lon/lat piece sf
# refuses counted in the sum shows here.
hand_plan <- function(a1, ratio, sf_refuses = rep(FALSE, length(ratio)))
  list(longlat = FALSE, a0 = a1 * ratio, a1 = a1, ratio = ratio, draws = NA_real_,
       sf_refuses = sf_refuses)

test_that("the caps take 10 n where it is larger", {
  expect_identical(spatialkit:::.sample_draw_cap(10), 1e6)
  expect_identical(spatialkit:::.sample_draw_cap(1e5), 1e6)
  expect_identical(spatialkit:::.sample_draw_cap(5e5), 5e6)
  expect_identical(spatialkit:::.sample_total_cap(10), 2e7)
  expect_identical(spatialkit:::.sample_total_cap(2e6), 2e7)
  expect_identical(spatialkit:::.sample_total_cap(3e6), 3e7)
})

test_that("the points needed in all are refused just above the cap and not at it", {
  chk <- function(plan, n) spatialkit:::.sample_total_check(plan, n, length(plan$a1), "f",
                                                            "nothing", "Remedy.")
  # n = 1000: two pieces of equal area take 500 points each, so they need
  # 500 (r + 1) points in their boxes; the cap is 2e7.
  expect_null(chk(hand_plan(c(1, 1), c(39999, 1)), 1000))
  expect_error(chk(hand_plan(c(1, 1), c(40001, 1)), 1000),
               "about 20,000,000 points drawn in the boxes of its polygons, more than the 20,000,000",
               fixed = TRUE)
  # n = 3e6: the cap is 10 n = 3e7; the pieces need 1.5e6 (r + 1).
  expect_null(chk(hand_plan(c(1, 1), c(19, 1)), 3e6))
  expect_error(chk(hand_plan(c(1, 1), c(21, 1)), 3e6), "more than the 30,000,000 drawn at most",
               fixed = TRUE)
  # A lon/lat piece sf refuses is not counted: drawing it keeps sf's own
  # refusal (BUG-186).
  expect_null(chk(hand_plan(c(1, 1), c(1e9, 1), c(TRUE, FALSE)), 1000))
  expect_error(chk(hand_plan(c(1, 1), c(1e9, 1), c(FALSE, FALSE)), 1000), "cannot be sampled")
})

test_that("a compact boundary with a large n stays on the whole draw", {
  local_vsize()
  nc_pr <- sf::st_union(sf::st_transform(nc, 32617))
  n <- 6e5
  # Its box is about twice its area: its draw is above 10^6 points and
  # within 10 n.
  p <- spatialkit:::.sample_draw_plan(nc_pr, n)
  expect_gt(p$draws, 1e6)
  expect_lte(p$draws, 10 * n)
  # No point is drawn: sf's draw is replaced by one that records its size.
  sizes <- numeric(0)
  testthat::local_mocked_bindings(
    .st_sample_quiet = function(x, size, ...) {
      sizes <<- c(sizes, size)
      sf::st_sfc(sf::st_point(c(0, 0)), crs = sf::st_crs(x))[rep(1L, size)]
    },
    .sample_by_piece = function(...) stop("the piece-by-piece path ran"),
    .sample_bounded = function(...) stop("the bounded draw ran"),
    .package = "spatialkit")
  got <- spatialkit:::.robust_st_sample(nc_pr, n, "x", "nothing")
  expect_length(got, n)
  expect_identical(sizes, n)
})

test_that("pieces of equal area take equal shares, however much of their boxes they fill", {
  local_vsize()
  # A strip 20 m wide (box 1250 times its area) and a square of 1 km, both
  # 10^6 m^2, far apart: the whole is above the total cap, so they go piece
  # by piece, and the strip is drawn in bounded draws.
  b <- layer(strip_at(20), sq(8e5, 3e6, 1000), crs = 32617)
  a <- as.numeric(sf::st_area(b))
  expect_equal(a[1], a[2], tolerance = 1e-6)
  expect_gt(plan_of(b, 3000)$draws, spatialkit:::.sample_total_cap(3000))
  s <- quiet(voronoi_seeds_random(b, 3000, set_seed = 1))
  expect_identical(nrow(s), 3000L)
  piece <- vapply(sf::st_intersects(s, b), `[`, integer(1), 1L)
  n_in <- tabulate(piece, 2L)
  # Shared by box area, the strip would take one seed in 1251.
  expect_gt(stats::chisq.test(n_in, p = c(0.5, 0.5))$p.value, 0.01)
  expect_true(all(n_in > 1300L))
})

test_that("a bounded draw returns exactly the points asked for", {
  local_vsize()
  # sf::st_sample(exact = FALSE) returns a random number of points, so the
  # last draw often goes over.  A triangle fills half its box; a ratio of
  # 1e5 passed in makes each draw ask for 10 points.
  tri <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(100, 0), c(0, 100), c(0, 0)))),
                    crs = 32617)
  for (s in 1:20) {
    set.seed(s)
    expect_length(spatialkit:::.sample_bounded(tri, 95, 1e5, "x"), 95L)
  }
  # The strip: k = 1000 takes 1.25e6 points in its box, in two draws.
  b <- strip(20)
  for (s in 1:2) {
    v <- quiet(voronoi_seeds_random(b, 1000, set_seed = s))
    expect_identical(nrow(v), 1000L)
    expect_true(inside(v, b))
  }
  g <- quiet(get_voronoi_seeds(b, n = 1000, method = "random", set_seed = 1))
  expect_identical(nrow(g), 1000L)
})

test_that("a bounded draw refuses a ratio that would draw without end", {
  tri <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(100, 0), c(0, 100), c(0, 0)))),
                    crs = 32617)
  for (r in list(Inf, NaN, 0, 0.5))
    expect_error(spatialkit:::.sample_bounded(tri, 10, r, "f"), "^f\\(\\): internal error")
  # A ring of no area on a slant has a box with an area: its ratio is 1,
  # so sf draws n points, as sf itself does.
  flat <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 1), c(2, 2), c(0, 0)))),
                     crs = 32617)
  p <- spatialkit:::.sample_draw_plan(flat, 10)
  expect_identical(p$ratio, 1)
  expect_identical(p$draws, 10)
})

test_that("many small pieces just above the cap are drawn whole in a few bounded draws", {
  local_vsize()
  # 1000 squares of 1 km spread over 1000 km: the whole box is about 1000
  # times their area, so 1500 seeds take 1.5e6 points in the box.  Piece by
  # piece that is about a thousand calls of sf; whole, a few draws (each
  # draw keeps a random number of points, so a few more top up the last).
  set.seed(9)
  xs <- stats::runif(1000, 0, 1e6)
  ys <- stats::runif(1000, 0, 1e6)
  many <- sf::st_sf(geometry = sf::st_sfc(lapply(seq_len(1000), function(i)
    sq(xs[i], ys[i], 1000)), crs = 3857))
  p <- plan_of(many, 1500)
  expect_gt(p$draws, spatialkit:::.sample_draw_cap(1500))
  expect_lte(p$draws, 2e6)
  log <- local_draw_log()
  v <- quiet(voronoi_seeds_random(many, 1500, set_seed = 1))
  expect_identical(nrow(v), 1500L)
  expect_true(inside(v, many))
  expect_lte(length(log$box), 20L)
  expect_lte(max(log$box), 1e6)
  expect_lte(sum(log$box), 2e6)
})

test_that("the whole is drawn in bounded draws only where that costs less and fits the cap", {
  whole <- function(draws, refuses = FALSE)
    list(longlat = FALSE, ratio = draws / 1000, draws = draws, sf_refuses = refuses)
  wb <- function(w, plan) spatialkit:::.sample_whole_bounded(w, plan, 1000)
  many  <- hand_plan(rep(1, 1000), rep(1, 1000))
  three <- hand_plan(rep(1, 3), rep(1, 3))
  # 1000 pieces that fill their boxes take 1000 points in 1000 calls; the
  # whole takes 1.5e6 points in two.
  expect_true(wb(whole(1.5e6), many))
  # Three pieces take three calls: piece by piece costs less.
  expect_false(wb(whole(1.5e6), three))
  # The whole above the total cap, or refused by sf, is never drawn whole.
  expect_false(wb(whole(2.1e7), many))
  expect_false(wb(whole(1.5e6, refuses = TRUE), many))
  # Pieces whose boxes take more than the cap in all, and a whole within it.
  expect_true(wb(whole(1.5e7), hand_plan(c(1, 1), c(1e5, 1))))
})
