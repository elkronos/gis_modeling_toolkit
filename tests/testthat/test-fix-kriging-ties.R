# tests/testthat/test-fix-kriging-ties.R
# ---------------------------------------------------------------------------
# kriging_adequacy() on gridded locations.  Many locations of a lattice are
# the same distance from a cell's centre or a held-out point, so which of
# them make up the nmax of a kriging system was settled by the last bits of
# the coordinates (gstat's own search), which move with the unit and the
# origin, and by the row order of the layer.  The nmax are now chosen by the
# package's tie rule, and the rows are put in order by coordinates, cell ID
# and response first.  A location whose rows name several cells is given
# none, with a warning.  The points are read in x and y on every layer, Z
# and M dropped, as the variogram is fitted.
# ---------------------------------------------------------------------------

skip_if_not_installed("gstat")

.kt_quiet <- function(expr) suppressWarnings(suppressMessages(
  logger::with_log_threshold(expr, threshold = logger::FATAL, namespace = "spatialkit",
                             index = 2)))

.kt_catch <- function(expr) {
  got <- list()
  value <- logger::with_log_threshold(
    withCallingHandlers(expr, warning = function(w) {
      got[[length(got) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }), threshold = logger::FATAL, namespace = "spatialkit", index = 2)
  list(value = value,
       classes = vapply(got, function(w) class(w)[1L], character(1)),
       texts = vapply(got, conditionMessage, character(1)))
}

.kt_cols <- c("n", "mean", "se", "plain_var", "kr_pred", "kr_var", "kr_ratio",
              "kr_gain", "kr_shift", "kr_n_used")
.kt_table <- function(ka) sf::st_drop_geometry(ka)[, .kt_cols]
.kt_cv <- function(ka) unlist(attr(ka, "cv")[c("zscore_var", "zscore_mean", "rmse", "n_pred")])

.kt_store <- new.env(parent = emptyenv())

# A 20 x 16 lattice at 250 m in EPSG:26915 with a short-range response, and
# the same lattice jittered by up to 60 m (no tied distances); cells of 2 x 2
# locations; the variogram fitted in metres.  The cross-validation runs on
# four folds of 5 x 4 blocks of the lattice, labelled by row, so that every
# frame has the same folds and only the kriging can differ.
.kt_frames_layer <- function(jitter) {
  key <- paste0("frames", jitter)
  if (is.null(.kt_store[[key]])) {
    g <- expand.grid(i = 0:19, j = 0:15)
    jx <- if (jitter) withr::with_seed(99, stats::runif(nrow(g), -60, 60)) else 0
    jy <- if (jitter) withr::with_seed(98, stats::runif(nrow(g), -60, 60)) else 0
    x <- 4.5e5 + g$i * 250 + jx; y <- 4.55e6 + g$j * 250 + jy
    z <- withr::with_seed(3, sin((x - min(x)) / 300 + 0.4) * cos((y - min(y)) / 260) +
                            stats::rnorm(length(x), 0, 0.4))
    fold <- (g$i %/% 5L) %% 2L * 2L + (g$j %/% 4L) %% 2L + 1L
    pts <- sf::st_as_sf(data.frame(x = x, y = y, z = z, fold = fold), coords = c("x", "y"),
                        crs = 26915)
    box <- sf::st_bbox(c(xmin = 4.5e5 - 125, ymin = 4.55e6 - 125,
                         xmax = 4.5e5 + 4875, ymax = 4.55e6 + 3875), crs = sf::st_crs(26915))
    cells <- sf::st_sf(cell_id = 1:80,
                       geometry = sf::st_make_grid(sf::st_as_sfc(box), n = c(10, 8)))
    asg <- assign_features_to_polygons(pts, cells, polygon_id_col = "cell_id", quiet = TRUE)
    sac <- .kt_quiet(estimate_sac_range(pts, "z"))
    .kt_store[[key]] <- list(asg = asg, cells = cells, sac = sac)
  }
  .kt_store[[key]]
}

# The layer and its variogram in another frame: US feet, or the false
# origin moved by 0.1 m.  The variogram carries that frame's CRS and unit,
# so the points stay in it.
.kt_frame <- function(L, frame) {
  crs <- switch(frame,
    m    = sf::st_crs(26915),
    usft = sf::st_crs("+proj=utm +zone=15 +datum=NAD83 +units=us-ft +no_defs"),
    m01  = sf::st_crs(paste0("+proj=tmerc +lat_0=0 +lon_0=-93 +k=0.9996 +x_0=500000.1 ",
                             "+y_0=0.1 +datum=NAD83 +units=m +no_defs")))
  if (frame == "m") return(L)
  s <- L$sac
  attr(s, "crs") <- crs
  if (frame == "usft") {
    f <- 3937 / 1200
    s[] <- as.numeric(s) * f
    vm <- attr(s, "variogram_model"); vm$range <- vm$range * f
    attr(s, "variogram_model") <- vm
  }
  list(asg = sf::st_transform(L$asg, crs), cells = sf::st_transform(L$cells, crs), sac = s)
}

.kt_krige_frame <- function(L, nmax)
  .kt_quiet(kriging_adequacy(L$asg, "z", L$cells, id_col = "cell_id", sac = L$sac,
                             folds = L$asg$fold, nmax = nmax))

test_that("kriging_adequacy() on a lattice gives the same cells and cross-validation in US feet and with the origin moved", {
  # nmax = 12 keeps every cell's 4 locations among the 12 nearest its centre;
  # nmax = 3 enlarges every cell.  Before, gstat's search settled the ties:
  # kr_pred moved by up to 0.03 in about half the cells.
  L <- .kt_frames_layer(FALSE)
  for (nmax in c(12L, 3L)) {
    ref <- .kt_krige_frame(L, nmax)
    if (nmax == 3L) expect_true(all(ref$kr_n_used > 3L))
    for (frame in c("usft", "m01")) {
      got <- .kt_krige_frame(.kt_frame(L, frame), nmax)
      expect_identical(got$kr_n_used, ref$kr_n_used)
      expect_equal(.kt_table(got), .kt_table(ref), tolerance = 1e-9)
      expect_equal(.kt_cv(got), .kt_cv(ref), tolerance = 1e-9)
    }
  }
})

# An 18 x 8 lattice at 250 m in EPSG:26915, cells of 3 x 1 locations, so a
# cell's two end locations are as far from its centre as the locations
# above and below its middle one: with nmax = 4 whether a cell holds a
# location its nmax nearest leave out rests on a tie at the furthest own
# location.  Folds of 6 x 8 blocks, labelled by row.
.kt_strip_layer <- function() {
  if (is.null(.kt_store$strip)) {
    g <- expand.grid(i = 0:17, j = 0:7)
    x <- 4.5e5 + g$i * 250; y <- 4.55e6 + g$j * 250
    z <- withr::with_seed(3, sin((x - min(x)) / 300 + 0.4) * cos((y - min(y)) / 260) +
                            stats::rnorm(length(x), 0, 0.4))
    pts <- sf::st_as_sf(data.frame(x = x, y = y, z = z, fold = (g$i %/% 6L) %% 3L + 1L),
                        coords = c("x", "y"), crs = 26915)
    box <- sf::st_bbox(c(xmin = 4.5e5 - 125, ymin = 4.55e6 - 125,
                         xmax = 4.5e5 + 4375, ymax = 4.55e6 + 1875), crs = sf::st_crs(26915))
    cells <- sf::st_sf(cell_id = 1:48,
                       geometry = sf::st_make_grid(sf::st_as_sfc(box), n = c(6, 8)))
    asg <- assign_features_to_polygons(pts, cells, polygon_id_col = "cell_id", quiet = TRUE)
    sac <- .kt_quiet(estimate_sac_range(pts, "z"))
    .kt_store$strip <- list(asg = asg, cells = cells, sac = sac)
  }
  .kt_store$strip
}

test_that("whether a cell is enlarged is counted with the tie rule's tolerance, the same in US feet and with the origin moved", {
  # Before, and with the count made without the tolerance, the end
  # locations' last bits decided whether the cell was enlarged: kr_n_used
  # differed between metres and US feet in most of the 48 cells.
  L <- .kt_strip_layer()
  expect_identical(as.vector(table(L$asg$cell_id)), rep(3L, 48))
  ref <- .kt_krige_frame(L, 4L)
  for (frame in c("usft", "m01")) {
    got <- .kt_krige_frame(.kt_frame(L, frame), 4L)
    expect_identical(got$kr_n_used, ref$kr_n_used, label = frame)
    expect_equal(.kt_table(got), .kt_table(ref), tolerance = 1e-9, label = frame)
    expect_equal(.kt_cv(got), .kt_cv(ref), tolerance = 1e-9, label = frame)
  }
})

test_that("on the jittered lattice the frames agree and the result is what it was", {
  # No distances tie, so gstat's own search and the rule pick the same sets.
  # The figures were taken before the change of the neighbour choice.
  L <- .kt_frames_layer(TRUE)
  ref <- .kt_krige_frame(L, 12L)
  for (frame in c("usft", "m01")) {
    got <- .kt_krige_frame(.kt_frame(L, frame), 12L)
    expect_equal(.kt_table(got), .kt_table(ref), tolerance = 1e-9)
    expect_equal(.kt_cv(got), .kt_cv(ref), tolerance = 1e-9)
  }
  expect_equal(sum(ref$kr_pred), 1.8915855223, tolerance = 1e-9)
  expect_equal(sum(ref$kr_var), 1.6506119569, tolerance = 1e-9)
  expect_equal(unname(.kt_cv(ref)[1:3]), c(1.0642127840, -0.0087527581, 0.6619738610),
               tolerance = 1e-9)
})

test_that("kriging each cell and held-out location on its own from the rule's pick gives what gstat's own search gives where the nmax are one set", {
  # The two routes a kriging system can take: gstat's nmax search in one
  # call, and a call of its own on the locations .tie_pick() chose with
  # nmax = Inf.  On the jittered lattice every system has one answer, so
  # forcing every one onto the second route changes nothing but rounding.
  L <- .kt_frames_layer(TRUE)
  ref <- .kt_krige_frame(L, 12L)
  pick <- .tie_pick
  local_mocked_bindings(.tie_pick = function(d2, k, x, y) {
    p <- pick(d2, k, x, y)
    p$tied <- TRUE
    p
  })
  got <- .kt_krige_frame(L, 12L)
  expect_identical(got$kr_n_used, ref$kr_n_used)
  expect_equal(.kt_table(got), .kt_table(ref), tolerance = 1e-9)
  expect_equal(.kt_cv(got), .kt_cv(ref), tolerance = 1e-9)
})

# A 24 x 24 lattice at 100 m in EPSG:32615, 3 x 3 cells, nmax = 20; and the
# same jittered by up to 20 m.
.kt_order_layer <- function(jitter) {
  key <- paste0("order", jitter)
  if (is.null(.kt_store[[key]])) {
    g <- expand.grid(x = seq(50, 2350, by = 100), y = seq(50, 2350, by = 100))
    withr::with_seed(3, {
      if (jitter) {
        g$x <- g$x + stats::runif(nrow(g), -20, 20)
        g$y <- g$y + stats::runif(nrow(g), -20, 20)
      }
      g$z <- sin(g$x / 700) + cos(g$y / 900) + stats::rnorm(nrow(g), 0, 0.3)
    })
    pts <- sf::st_as_sf(g, coords = c("x", "y"), crs = 32615)
    box <- sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 2400, ymax = 2400), crs = sf::st_crs(32615))
    cells <- sf::st_sf(poly_id = 1:9, geometry = sf::st_make_grid(sf::st_as_sfc(box), n = c(3, 3)))
    asg <- assign_features_to_polygons(pts, cells, quiet = TRUE)
    sac <- .kt_quiet(estimate_sac_range(pts, "z"))
    .kt_store[[key]] <- list(asg = asg, cells = cells, sac = sac)
  }
  .kt_store[[key]]
}
.kt_rows <- function(x, o) x[o, , drop = FALSE]

test_that("kriging_adequacy() gives the same result for a lattice in any row order", {
  # Before, the rows reversed moved kr_pred in 5 of 9 cells and the
  # cross-validation in the fourth digit.
  for (jitter in c(FALSE, TRUE)) {
    L <- .kt_order_layer(jitter)
    n <- nrow(L$asg)
    run <- function(o) .kt_catch(kriging_adequacy(.kt_rows(L$asg, o), "z", L$cells,
                                                  sac = L$sac, nmax = 20))
    ref <- run(seq_len(n))
    expect_identical(run(rev(seq_len(n))), ref)
    expect_identical(run(withr::with_seed(5, sample(n))), ref)
  }
})

test_that("a location whose rows name two cells is given neither, with a warning, in any row order", {
  # One location of cell 5 given a second row labelled cell 6.  Before, the
  # location was cell 5's or cell 6's by which row came first, and kr_pred of
  # both cells moved with the order, with nothing said.
  L <- .kt_order_layer(TRUE)
  extra <- L$asg[L$asg$poly_id == 5, ][1, ]
  extra$poly_id <- 6L
  extra$z <- extra$z + 0.1
  b <- rbind(L$asg, extra)
  run <- function(x) .kt_catch(kriging_adequacy(x, "z", L$cells, sac = L$sac, nmax = 20))
  ref <- run(b)
  rev_b <- .kt_rows(b, rev(seq_len(nrow(b))))
  expect_identical(run(rev_b), ref)
  hit <- grepl("location(s) carry rows with more than one cell ID", ref$texts, fixed = TRUE)
  expect_equal(sum(hit), 1L)
  expect_match(ref$texts[hit], "^kriging_adequacy\\(\\): 1 location\\(s\\)")
  # Each row still counts in the cell it names.
  ka <- ref$value
  expect_identical(ka$n[ka$poly_id == 6] - sum(L$asg$poly_id == 6), 1L)
  expect_identical(ka$n[ka$poly_id == 5], sum(L$asg$poly_id == 5))
  # Without the second row, no warning of the kind.
  plain <- run(L$asg)
  expect_false(any(grepl("more than one cell ID", plain$texts, fixed = TRUE)))
  # The location is in neither cell's own set.  nmax = 20 enlarges every
  # cell to its own locations plus 20 outside it, so cell 5 is kriged from
  # one location fewer than without the second row, and cell 6 from as many.
  pv <- plain$value
  expect_true(all(pv$kr_n_used == as.vector(table(L$asg$poly_id)) + 20L))
  expect_identical(ka$kr_n_used[ka$poly_id == 5], pv$kr_n_used[pv$poly_id == 5] - 1L)
  expect_identical(ka$kr_n_used[ka$poly_id == 6], pv$kr_n_used[pv$poly_id == 6])
})

test_that("a vector of fold labels is still read by the rows' positions as given", {
  # The rows are put in coordinate order before the folds are read, so the
  # labels have to follow them: a vector of labels in the layer's order and
  # the same layer and labels reversed give the same folds.
  L <- .kt_order_layer(TRUE)
  n <- nrow(L$asg)
  lab <- rep(1:4, length.out = n)
  ref <- .kt_quiet(kriging_adequacy(L$asg, "z", L$cells, sac = L$sac, folds = lab))
  got <- .kt_quiet(kriging_adequacy(.kt_rows(L$asg, rev(seq_len(n))), "z", L$cells,
                                    sac = L$sac, folds = rev(lab)))
  expect_identical(attr(got, "cv"), attr(ref, "cv"))
  expect_identical(attr(ref, "cv")$k, 4L)
  # Labels that follow the coordinates rather than the rows would give other folds.
  wrong <- .kt_quiet(kriging_adequacy(L$asg, "z", L$cells, sac = L$sac,
                                      folds = lab[order(sf::st_coordinates(L$asg)[, 1])]))
  expect_false(identical(attr(wrong, "cv")$rmse, attr(ref, "cv")$rmse))
})

test_that("a location seen twice is held out under the fold label of its first row as given", {
  # The jittered 20 x 16 lattice, every fourth location given a second row
  # with a lower response and the next fold's label.  At such a location the
  # second row comes first in coordinate, cell ID and response order, but
  # the label read is that of the row the caller gave first.
  L <- .kt_frames_layer(TRUE)
  a <- L$asg
  r <- a[seq(1, nrow(a), by = 4), ]
  r$z <- r$z - 0.5
  r$fold <- r$fold %% 4L + 1L
  cv <- function(x) attr(.kt_quiet(kriging_adequacy(x, "z", L$cells, id_col = "cell_id",
                                                     sac = L$sac, folds = x$fold,
                                                     nmax = 12)), "cv")
  ref <- cv(rbind(a, r))
  # The second rows' labels play no part.
  r_same <- r; r_same$fold <- a$fold[seq(1, nrow(a), by = 4)]
  expect_identical(cv(rbind(a, r_same)), ref)
  # Nor does the order of the rows, while each location's first row stays.
  expect_identical(cv(rbind(.kt_rows(a, withr::with_seed(4, sample(nrow(a)))),
                            .kt_rows(r, rev(seq_len(nrow(r)))))), ref)
  # With the second rows given first, their labels are the ones read.
  first_r <- cv(rbind(r, a))
  a_r <- a; a_r$fold[seq(1, nrow(a), by = 4)] <- r$fold
  expect_identical(cv(rbind(r, a_r)), first_r)
  expect_false(isTRUE(all.equal(first_r$rmse, ref$rmse)))
})

test_that(".canonical_order() does not depend on the order given or on the last bits", {
  x <- rep(c(0, 100, 200), each = 3); y <- rep(c(0, 100, 200), 3)
  id <- c("b", "a", "B", "c", "a", "A", "b", "b", "a")
  o <- .canonical_order(x, y, id)
  expect_identical(o, 1:9)
  p <- c(9L, 3L, 5L, 1L, 7L, 2L, 8L, 4L, 6L)
  expect_identical(p[.canonical_order(x[p], y[p], id[p])], o)
  # Last bits: x of one column moved by a few units in the last place.
  x2 <- x + c(0, 0, 0, 1e-12, -1e-12, 2e-12, 0, 0, 0)
  expect_identical(.canonical_order(x2, rev(y)), .canonical_order(x, rev(y)))
  # Rows at one place go by the further keys, text byte by byte.
  xs <- rep(5, 4); ys <- rep(7, 4)
  expect_identical(.canonical_order(xs, ys, c("b", "a", "B", NA), c(1, 2, 3, 4)),
                   c(3L, 2L, 1L, 4L))
})

# 300 points in EPSG:32617 on a smooth field, 16 grid cells.
.kt_zm_layer <- function() {
  if (is.null(.kt_store$zm)) {
    d <- withr::with_seed(7, {
      n <- 300
      x <- stats::runif(n, 0, 1e5); y <- stats::runif(n, 0, 1e5)
      cx <- stats::runif(30, 0, 1e5); cy <- stats::runif(30, 0, 1e5)
      z <- rowSums(sapply(1:30, function(k)
        exp(-((x - cx[k])^2 + (y - cy[k])^2) / (2 * 8e3^2)))) + stats::rnorm(n, 0, 0.1)
      data.frame(x = 5e5 + x, y = 5e6 + y, z = z, h = stats::runif(n, 0, 1000))
    })
    pts <- sf::st_as_sf(d[, c("x", "y", "z")], coords = c("x", "y"), crs = 32617)
    cells <- sf::st_sf(cell_id = 1:16,
                       geometry = sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(pts)), n = c(4, 4)))
    asg <- assign_features_to_polygons(pts, cells, polygon_id_col = "cell_id", quiet = TRUE)
    sac <- .kt_quiet(estimate_sac_range(pts, "z", cutoff = 1))
    sac_other <- .kt_quiet(estimate_sac_range(sf::st_transform(pts, 32616), "z", cutoff = 1))
    .kt_store$zm <- list(asg = asg, cells = cells, sac = sac, sac_other = sac_other, h = d$h)
  }
  .kt_store$zm
}
# The layer's points given a third coordinate: `third` per row, as Z or M.
.kt_with <- function(asg, third, dim) {
  xy <- sf::st_coordinates(asg)
  g <- lapply(seq_len(nrow(asg)), function(i) sf::st_point(c(xy[i, 1:2], third[i]), dim = dim))
  sf::st_set_geometry(asg, sf::st_sfc(g, crs = sf::st_crs(asg)))
}

test_that("kriging_adequacy() kriges points XYZ or XYM throughout as their XY copy", {
  # gstat::krige() kriged points XYZ throughout in three dimensions with a
  # variogram fitted in x and y: heights of 0 to 1000 m moved kr_var by
  # about a tenth.
  L <- .kt_zm_layer()
  variants <- list(z_random = .kt_with(L$asg, L$h, "XYZ"),
                   z_const = .kt_with(L$asg, rep(3, nrow(L$asg)), "XYZ"),
                   m_random = .kt_with(L$asg, L$h, "XYM"))
  for (sac_name in c("sac", "sac_other", "none")) {
    sac <- if (sac_name == "none") NULL else L[[sac_name]]
    run <- function(a) .kt_catch(kriging_adequacy(a, "z", L$cells, id_col = "cell_id",
                                                  sac = sac, seed = 1))
    ref <- run(L$asg)
    for (v in names(variants)) {
      got <- run(variants[[v]])
      expect_identical(got, ref, label = paste(v, sac_name))
      expect_false(any(got$classes == "spatialkit_zm_dropped"))
    }
  }
  # Reading an XY layer in x and y leaves it as it is.
  expect_identical(sf::st_zm(L$asg, drop = TRUE, what = "ZM"), L$asg)
})

test_that("?kriging_adequacy states the tie rule, its bound, the row order and the location named by two cells", {
  file <- test_path("..", "..", "man", "kriging_adequacy.Rd")
  rd <- if (file.exists(file)) tools::parse_Rd(file)
        else tools::Rd_db("spatialkit")[["kriging_adequacy.Rd"]]
  expect_s3_class(rd, "Rd")
  txt <- gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  for (s in c(
    "So the nmax are chosen by the package's tie rule: squared distances within one part in ten million of each other count as equal",
    "as long as no coordinate moves by more than about a two-hundred-millionth of the distance between neighbouring locations",
    "so the result does not depend on the order of the rows either",
    "A location whose rows name more than one cell (rows at one place carrying different IDs) is no cell's own",
    "Any Z or M dimension is dropped, on every layer"))
    expect_match(txt, s, fixed = TRUE)
})
