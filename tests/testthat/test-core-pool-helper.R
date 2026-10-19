# tests/testthat/test-core-pool-helper.R
# ---------------------------------------------------------------------------
# The reference functions of helper-pool-reference.R are the text they were
# copied as, and they run on the layouts the engine is checked on.
#
# The functions are a fixed point for the tests of R/pool.R, so the first
# test holds their text still.  The second holds the layouts still: their
# settings, and the coordinates, response and regions each gives.  The others
# run the functions on every layout and check what they return against the
# same quantities worked out row by row: a copy that had lost a line would
# not pass.
# ---------------------------------------------------------------------------

test_that("the copied block of the reference helper has not been edited", {
  lines <- sub("\r$", "", readLines(test_path("helper-pool-reference.R"), warn = FALSE))
  from <- grep("^# ---- core-factorial/20_estimates.R", lines)
  to   <- grep("^# ---- Layouts", lines)
  expect_length(from, 1L)
  expect_length(to, 1L)
  block <- paste(lines[from:(to - 1L)], collapse = "\n")
  # The block as it was copied from the experiments' files.  A change here is
  # a change to the reference, not a correction of it.
  expect_identical(digest::digest(block, algo = "sha256", serialize = FALSE),
                   "d3cc783800fa057f4762a5ea57f78c183833c69e816b766ef1b3c869b926529c")
  expect_identical(EST, c("P", "R", "L", "PAR"))
  expect_identical(NEIGHBOURS, 8L)
  for (f in c("oneway", "reml_oneway", "weight", "reldiff", "parent_of", "region_estimates",
              "predict_partition", "loo_plain", "pick_rung", "rung_k", "wfun", "reml_means",
              "fit_regions", "loo_int"))
    expect_true(is.function(get(f)), info = f)
})

test_that("the layouts are the thirteen of the engine's tests, and are the same on every call", {
  expect_length(POOL_LAYOUTS, 13L)
  # The settings, written out: rows, rows per region, clustered or not, the
  # number of regions where it is set, and the seed each is drawn under.
  stated <- data.frame(
    N         = c(600, 600, 600, 3000, 3000, 3000, 3000, 40, 12, 9, 60, 60, 90),
    per       = c(18, 4.5, 2, 4.5, 9, 3, 2, 4.5, 4.5, 4.5, NA, NA, NA),
    clustered = c(FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE,
                  FALSE, FALSE, TRUE),
    K         = c(rep(NA, 10), 2, 3, 4),
    seed      = c(826, 731.5, 714, 3131.5, 3163, 3121, 3114, 171.5, 143.5, 140.5, 5, 6, 7))
  pick <- function(l, nm) if (is.null(l[[nm]])) NA else l[[nm]]
  seed_of <- function(l) if (is.null(l$seed))
    eval(formals(pool_layout)$seed, list(N = l$N, per = l$per)) else l$seed
  expect_identical(vapply(POOL_LAYOUTS, function(l) as.numeric(l$N), numeric(1)), stated$N)
  expect_identical(vapply(POOL_LAYOUTS, function(l) as.numeric(pick(l, "per")), numeric(1)),
                   stated$per)
  expect_identical(vapply(POOL_LAYOUTS, function(l) l$clustered, logical(1)), stated$clustered)
  expect_identical(vapply(POOL_LAYOUTS, function(l) as.numeric(pick(l, "K")), numeric(1)),
                   stated$K)
  expect_identical(vapply(POOL_LAYOUTS, function(l) as.numeric(seed_of(l)), numeric(1)),
                   stated$seed)
  # What each layout is, from one run: the regions it has, the sums of its
  # coordinates and of its response, the sum of squares of the coordinates
  # about the regions' centroids, its first row, its largest region and its
  # regions of one row.  A change to the draws, to their order, to the
  # response or to the k-means fit moves one of them.
  pinned <- list(
    list(K = 33L, x = 29874706.2895214, y = 28792422.7186013, z = 140.246634888009,
         wss = 26439897402.9235, first = c(45574.0485340357, 11603.8269363344, 1.76399348173672),
         largest = 37L, alone = 0L,
         regions = "fd062acf5a72dade85aea1c49438c1130a36c56bf834399e75578b321b23328b"),
    list(K = 133L, x = 29715273.1257956, y = 29639554.288378, z = 83.4625253875932,
         wss = 4598880412.62402, first = c(81091.2917600945, 12437.1200334281, -1.09252907117771),
         largest = 10L, alone = 2L,
         regions = "36c86cd123197947c275c8ad86b8247b684094167e4e4cf4c949d4d2b9c9f270"),
    list(K = 300L, x = 39526653.8854346, y = 40511077.5110761, z = -60.6694540007996,
         wss = 316262114.417961, first = c(83434.6807526603, 74346.5297039839, -0.288675094614457),
         largest = 5L, alone = 107L,
         regions = "a024f5e5195dfebef544634998f11229b811bcb22d55f8f81a88144a26803d5f"),
    list(K = 667L, x = 136283944.421336, y = 197566173.908574, z = 1155.78586429603,
         wss = 1467946939.7624, first = c(68475.5925537745, 86202.6189963955, 0.67636526770625),
         largest = 13L, alone = 39L,
         regions = "ccaa86f8fe1e76b33827478ac0b87c174e2904954bb01260e04b1885992d2dee"),
    list(K = 333L, x = 148742602.818901, y = 149702122.419886, z = 576.068167900529,
         wss = 11185811292.6901, first = c(85614.8979626596, 17524.1708289832, -0.773295764742362),
         largest = 19L, alone = 0L,
         regions = "c0a4fef54697d3cf27179d46a6c787a7bb3f6834a1e70b3fc9b069e75c15d659"),
    list(K = 1000L, x = 135971717.360894, y = 130139974.012379, z = 913.383691126545,
         wss = 713700440.715928, first = c(51671.8402377044, 21982.1038269469, -0.70027314610776),
         largest = 9L, alone = 134L,
         regions = "cd404c3466e24e4902acb8cda89c10a5212ff7c0bc766fa3f5b6b470c738e6a8"),
    list(K = 1500L, x = 147678263.059002, y = 151589194.867504, z = 582.11146843136,
         wss = 1042374972.84667, first = c(69767.6845593378, 53949.3277901784, -0.939162754885976),
         largest = 7L, alone = 536L,
         regions = "a6f35a014d2e5c495a3644f910936e166a572b087511cf83268d1a00e43e4f54"),
    list(K = 9L, x = 1998256.63326774, y = 1874450.39760787, z = 9.79517962950725,
         wss = 4964049010.7189, first = c(16371.9046860933, 59489.6582653746, 0.456385099586902),
         largest = 8L, alone = 0L,
         regions = "bb67f9970cb629dda9732b1ef532d7caedf290d9ed8882f3a7978299779a6b4c"),
    list(K = 3L, x = 593708.877358586, y = 732662.02150844, z = 4.77098301529861,
         wss = 5621047760.91425, first = c(94322.5271301344, 95603.4328555688, -0.400924170733543),
         largest = 7L, alone = 0L,
         regions = "5e918ee4e0970ac5b38003809833bafa5565e92efbd74aab888f51fd71e1734d"),
    list(K = 2L, x = 655398.420360455, y = 363021.137519904, z = -4.56240929685216,
         wss = 5220797540.42634, first = c(91986.9286115127, 27714.2144500128, -1.06983221268626),
         largest = 7L, alone = 0L,
         regions = "954e7ac196f2aafcfdf9f38642edd8aee913bbd73b35bf11ced2c33693e2ad15"),
    list(K = 2L, x = 2969596.72642406, y = 3207225.74295942, z = 15.5424211974207,
         wss = 56285508381.0287, first = c(20021.4452575892, 62396.6502491385, -0.364812325134488),
         largest = 31L, alone = 0L,
         regions = "d3a490a96f0c08a035ed84c064030965e7ab5d0378ec488a95567b0f3bc96328"),
    list(K = 3L, x = 3396383.72983318, y = 2914499.51124378, z = -5.18568712226186,
         wss = 38488574836.6697, first = c(60626.8297648057, 25096.3960308582, 0.136392097837547),
         largest = 25L, alone = 0L,
         regions = "cde65ba054d2d4f3f6321e5e1069b7cb60663b509d1f63af208dcbb95b448381"),
    list(K = 4L, x = 4324757.12586078, y = 4356933.26675007, z = 49.7154318092928,
         wss = 37068536222.5694, first = c(79334.3344396412, 8184.14742803582, 0.0718807282317518),
         largest = 31L, alone = 0L,
         regions = "94a7db8adadb8a204add1b98ba92670fec824e8db4a538573b740d170411dc6e"))
  # Which row is in which region is held where the platform is the one the
  # numbers were taken on (the frozen fixture's): elsewhere a k-means fit
  # differs in its last bits.
  here <- core_frozen_platform()
  same_platform <- core_frozen_same_platform(
    here, readRDS(test_path("fixtures", "core-frozen.rds"))$info)
  for (i in seq_along(POOL_LAYOUTS)) {
    L   <- do.call(pool_layout, POOL_LAYOUTS[[i]])
    P   <- pinned[[i]]
    tag <- sprintf("layout %d", i)
    n   <- tabulate(L$gi, L$K)
    cen <- rowsum(L$xy, L$gi) / n
    expect_identical(L$K, P$K, info = tag)
    expect_equal(sum(L$xy[, 1]), P$x, tolerance = 1e-9, info = tag)
    expect_equal(sum(L$xy[, 2]), P$y, tolerance = 1e-9, info = tag)
    expect_equal(sum(L$y), P$z, tolerance = 1e-9, info = tag)
    expect_equal(sum((L$xy - cen[L$gi, ])^2), P$wss, tolerance = 1e-9, info = tag)
    expect_equal(c(L$xy[1, ], L$y[1]), P$first, tolerance = 1e-9, info = tag)
    expect_identical(max(n), P$largest, info = tag)
    expect_identical(sum(n == 1L), P$alone, info = tag)
    if (same_platform)
      expect_identical(digest::digest(paste(L$gi, collapse = ","), algo = "sha256",
                                      serialize = FALSE), P$regions, info = tag)
  }
  shape <- vapply(POOL_LAYOUTS, function(l) {
    L <- do.call(pool_layout, l)
    c(N = nrow(L$xy), K = L$K, rows = length(L$y), groups = length(unique(L$gi)),
      top = max(L$gi))
  }, numeric(5))
  expect_identical(unname(shape["N", ]),
                   c(600, 600, 600, 3000, 3000, 3000, 3000, 40, 12, 9, 60, 60, 90))
  # rung_k() of the rows and the size, less the centres that took no row.
  wanted <- vapply(POOL_LAYOUTS, function(l)
    if (is.null(l$K)) rung_k(l$N, l$per) else l$K, numeric(1))
  expect_true(all(shape["K", ] <= wanted))
  expect_true(all(shape["K", ] >= 0.9 * wanted))
  expect_identical(unname(shape["K", 11:13]), c(2, 3, 4))
  expect_identical(shape["groups", ], shape["K", ])
  expect_identical(shape["top", ], shape["K", ])
  # The regions are those of the nearest centre: every row is nearer the
  # centroid of its own region than that of any other.
  for (l in POOL_LAYOUTS) {
    L   <- do.call(pool_layout, l)
    cen <- rowsum(L$xy, L$gi) / tabulate(L$gi)
    d2  <- outer(L$xy[, 1], cen[, 1], "-")^2 + outer(L$xy[, 2], cen[, 2], "-")^2
    expect_gte(mean(max.col(-d2, ties.method = "first") == L$gi), 0.98)
  }
  # A layout does not depend on the session's random numbers, and leaves them.
  set.seed(1); a <- do.call(pool_layout, POOL_LAYOUTS[[3L]]); after <- stats::runif(1)
  set.seed(2); b <- do.call(pool_layout, POOL_LAYOUTS[[3L]])
  expect_identical(a, b)
  set.seed(1)
  expect_identical(stats::runif(1), after)
  # Nor on the kind of generator the session has set, which it leaves too.
  old_kind <- RNGkind()
  old_seed <- get(".Random.seed", envir = globalenv())
  withr::defer({
    suppressWarnings(do.call(RNGkind, as.list(old_kind)))
    assign(".Random.seed", old_seed, envir = globalenv())
  })
  other_kind <- c("L'Ecuyer-CMRG", "Box-Muller", "Rounding")
  suppressWarnings(do.call(RNGkind, as.list(other_kind)))
  set.seed(1)
  drawn <- stats::runif(1)
  set.seed(1)
  expect_identical(do.call(pool_layout, POOL_LAYOUTS[[3L]]), a)
  # R-devel reports a fourth kind (the binomial generator); the three set here
  # are the first three on every version.
  expect_identical(RNGkind()[seq_along(other_kind)], other_kind)
  expect_identical(stats::runif(1), drawn)
  suppressWarnings(do.call(RNGkind, as.list(old_kind)))
  # Clustered and uniform layouts differ in how even the regions are.
  cv <- function(l) { n <- tabulate(do.call(pool_layout, l)$gi); stats::sd(n) / mean(n) }
  expect_gt(cv(POOL_LAYOUTS[[4L]]), cv(POOL_LAYOUTS[[2L]]))
})

test_that("the reference functions run on every layout and agree with a row-by-row count", {
  skip_if_not_installed("FNN")
  for (l in POOL_LAYOUTS) {
    L   <- do.call(pool_layout, l)
    y   <- L$y; gi <- L$gi; K <- L$K; xy <- L$xy; N <- length(y)
    tag <- sprintf("N = %d, K = %d", N, K)
    nn1 <- pool_nearest_other(xy)
    # Another row, and the nearest one.
    expect_false(any(nn1 == seq_len(N)), info = tag)
    if (N <= 600L) {
      d <- as.matrix(stats::dist(xy))
      diag(d) <- Inf
      expect_identical(nn1, unname(apply(d, 1L, which.min)), info = tag)
    }
    R   <- region_estimates(y, gi, K, xy, want = c("P", "L"))
    Fr  <- fit_regions(gi, gi[1:3], y, xy, xy[1:3, , drop = FALSE], nn1)

    # The plain mean and the one-way parts.
    n <- tabulate(gi, K)
    m <- vapply(seq_len(K), function(j) mean(y[gi == j]), numeric(1))
    expect_equal(unname(R$est[, "P"]), m, tolerance = 1e-12, info = tag)
    expect_identical(R$a$n, n, info = tag)
    expect_equal(R$a$sigma2, sum((y - m[gi])^2) / (N - K), tolerance = 1e-12, info = tag)

    # The pooled value: its neighbours, target, weight and the value itself.
    k_nb <- min(8L, K - 1L)
    expect_identical(dim(Fr$nb), c(K, k_nb), info = tag)
    expect_false(any(Fr$nb == seq_len(K)), info = tag)
    cen <- t(vapply(seq_len(K), function(j) colMeans(xy[gi == j, , drop = FALSE]), numeric(2)))
    for (j in unique(c(1L, K, round(seq(1, K, length.out = 6))))) {
      d <- sqrt((cen[, 1] - cen[j, 1])^2 + (cen[, 2] - cen[j, 2])^2)
      expect_setequal(Fr$nb[j, ], order(d)[2:(k_nb + 1L)])
    }
    target <- vapply(seq_len(K), function(j) mean(y[gi %in% Fr$nb[j, ]]), numeric(1))
    expect_equal(Fr$t, target, tolerance = 1e-10, info = tag)
    tau2 <- max(0, (sum(n * (m - target)^2) - K * R$a$sigma2) / N)
    expect_equal(Fr$tau2, tau2, tolerance = 1e-9, info = tag)
    expect_true(all(Fr$w >= 0 & Fr$w <= 1), info = tag)
    expect_identical(unname(Fr$est[, "L"]), unname(R$est[, "L"]), info = tag)
    expect_equal(unname(R$est[, "L"]), Fr$w * m + (1 - Fr$w) * target, tolerance = 1e-10,
                 info = tag)

    # Leave one row out, for a sample of rows: the region's count and mean
    # again without the row, the variances and the targets held.
    rows <- if (N <= 600L) seq_len(N) else round(seq(1, N, length.out = 300))
    e_mean <- e_pool <- numeric(length(rows))
    for (s in seq_along(rows)) {
      i <- rows[s]; j <- gi[i]
      rest <- y[gi == j & seq_len(N) != i]
      if (length(rest)) {
        w <- if (tau2 + R$a$sigma2 / length(rest) > 0)
          tau2 / (tau2 + R$a$sigma2 / length(rest)) else 1
        e_mean[s] <- y[i] - mean(rest)
        e_pool[s] <- y[i] - (w * mean(rest) + (1 - w) * target[j])
      } else {
        e_mean[s] <- y[i] - m[gi[nn1[i]]]
        e_pool[s] <- y[i] - target[j]
      }
    }
    expect_equal(Fr$eP[rows], e_mean, tolerance = 1e-9, info = tag)
    expect_equal(Fr$eL[rows], e_pool, tolerance = 1e-9, info = tag)
    expect_identical(Fr$loo_P, loo_plain(y, gi, R$a, nn1), info = tag)
    expect_identical(Fr$loo_P, mean(Fr$eP^2), info = tag)
    expect_identical(Fr$loo_L, mean(Fr$eL^2), info = tag)
    expect_identical(Fr$loo_L_se, stats::sd(Fr$eL^2) / sqrt(N), info = tag)

    # The interval from the residuals.
    q <- stats::quantile(Fr$eL, c(0.05, 0.95), names = FALSE, type = 7)
    expect_identical(loo_int(Fr$est[, "L"], Fr$eL, 1, 1),
                     cbind(Fr$est[, "L"] + q[1], Fr$est[, "L"] + q[2]), info = tag)
  }
})

test_that("the rest of the reference functions do what their names say", {
  skip_if_not_installed("FNN")
  # The count at a size: n over m, rounded to the even number at a half, and
  # never below two.
  expect_identical(rung_k(600, 4.5), 133L)
  expect_identical(rung_k(9, 2), 4L)        # 4.5 rounds to 4
  expect_identical(rung_k(27, 2), 14L)      # 13.5 rounds to 14
  expect_identical(rung_k(12, 18), 2L)
  expect_identical(pick_rung(c(3, 1, 1, 2)), 2L)
  expect_identical(weight(0, 0, 3), 1)
  expect_identical(weight(1, 2, 2), 0.5)
  expect_identical(wfun(0, 0), 1)
  expect_identical(wfun(3, 1), 0.75)
  expect_identical(reldiff(c(1, 2, -4), c(1, 4, -2)), c(0, 0.5, 0.5))

  L <- do.call(pool_layout, POOL_LAYOUTS[[2L]])
  # A response with no structure, drawn so that the moment estimate of the
  # between-region variance is negative: it is set to 0, every weight is 0
  # and the pooled value is the target.
  noise <- withr::with_seed(903, stats::rnorm(length(L$y)), .rng_kind = "Mersenne-Twister",
                            .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
  Fn <- fit_regions(L$gi, L$gi[1:3], noise, L$xy, L$xy[1:3, , drop = FALSE],
                    pool_nearest_other(L$xy))
  expect_lt(Fn$Q - Fn$K * Fn$sigma2, 0)
  expect_identical(Fn$tau2, 0)
  expect_identical(Fn$w, rep(0, Fn$K))
  expect_identical(unname(Fn$est[, "L"]), Fn$t)
  expect_identical(unname(region_estimates(noise, L$gi, L$K, L$xy, want = "L")$est[, "L"]),
                   Fn$t)

  # The nearest other row is not defined by the reference where a location
  # repeats: FNN returns a row as its own neighbour there.
  expect_error(pool_nearest_other(rbind(L$xy[1:50, ], L$xy[1:50, ])),
               "two rows share a location")

  a <- oneway(L$y, L$gi, L$K)
  r <- reml_oneway(a)
  expect_true(is.finite(r$sigma2) && r$sigma2 > 0 && is.finite(r$tau2) && r$tau2 >= 0)
  expect_error(oneway(L$y[1:5], 1:5, 5L), "no group holds two")
  rm0 <- reml_means(c(1, 1, 1), c(1, 1, 1))
  expect_identical(rm0$tau2, 0)
  expect_gt(reml_means(c(-5, 0, 5, 10), rep(0.01, 4))$tau2, 1)

  # A held-out row in a region with no training row takes the value of the
  # region of its nearest training row.
  hold <- L$gi %in% c(1L, 2L) | seq_along(L$y) %% 7 == 0
  P <- predict_partition(L$gi[!hold], L$gi[hold], L$y[!hold], L$xy[!hold, ], L$xy[hold, ],
                         want = c("P", "R", "L"))
  expect_identical(P$n_empty, sum(L$gi[hold] %in% c(1L, 2L)))
  expect_identical(dim(P$pred), c(sum(hold), 3L))
  expect_true(all(is.finite(P$pred)))
  tr_mean <- vapply(split(L$y[!hold], L$gi[!hold]), mean, numeric(1))
  te <- which(hold)
  kept <- !(L$gi[te] %in% c(1L, 2L))
  expect_equal(unname(P$pred[kept, "P"]), unname(tr_mean[as.character(L$gi[te][kept])]),
               tolerance = 1e-12)
  lost <- te[!kept][1L]
  tr <- which(!hold)
  near <- tr[which.min((L$xy[tr, 1] - L$xy[lost, 1])^2 + (L$xy[tr, 2] - L$xy[lost, 2])^2)]
  expect_equal(unname(P$pred[which(te == lost), "P"]),
               unname(tr_mean[as.character(L$gi[near])]), tolerance = 1e-12)
  expect_identical(P$n_own[which(te == lost)], 0L)
})
