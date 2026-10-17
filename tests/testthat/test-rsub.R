# Racing subordinators (Math/RDMSWTN_Derivation.tex).  Dense references are
# the inverse-Gaussian closed forms written directly with pnorm(); packet
# references are brute-force Skellam sums over Poisson pmfs.  Uniform start
# points are integrated by quadrature over the distance.

rsub_dense_cdf <- function(u, v, d, s) {
  c <- s^2; z <- sqrt(c * d)
  pnorm((v * u - d) / z) - exp(2 * v * u / c + pnorm(-(v * u + d) / z, log.p = TRUE))
}
rsub_dense_pdf <- function(u, v, d, s) {
  c <- s^2; z <- sqrt(c * d)
  2 * v / z * dnorm((v * u - d) / z) -
    2 * v / c * exp(2 * v * u / c + pnorm(-(v * u + d) / z, log.p = TRUE))
}
# Packets: N ~ Pois(2 v u / c) events have arrived, K ~ Pois(2 d / c) packets
# beyond the first are needed.  log P(N - K >= 1) / log P(N - K <= 0) and the
# rate times P(N = K).
# The joint mass of N = K peaks near sqrt(a b), so sum past both Poisson tails.
rsub_pk_k <- function(a, b) 0:(qpois(1e-17, max(a, b), lower.tail = FALSE) + 30)
rsub_lse <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
rsub_pk_lF <- function(u, v, d, s, upper = TRUE) {
  c <- s^2; a <- 2 * v * u / c; b <- 2 * d / c; k <- rsub_pk_k(a, b)
  rsub_lse(dpois(k, b, log = TRUE) + ppois(k, a, lower.tail = !upper, log.p = TRUE))
}
rsub_pk_cdf <- function(u, v, d, s)
  vapply(seq_along(d), function(i) exp(rsub_pk_lF(u, v, d[i], s)), 0)
rsub_pk_pdf <- function(u, v, d, s) {
  c <- s^2; a <- 2 * v * u / c
  vapply(d, function(di) {
    b <- 2 * di / c; k <- rsub_pk_k(a, b)
    2 * v / c * exp(rsub_lse(dpois(k, a, log = TRUE) + dpois(k, b, log = TRUE)))
  }, 0)
}
rsub_ref <- function(packets) {
  if (packets) list(cdf = rsub_pk_cdf, pdf = rsub_pk_pdf)
  else list(cdf = rsub_dense_cdf, pdf = rsub_dense_pdf)
}
# Uniform start point by quadrature over the distance.
rsub_ref_avg <- function(fn, u, v, B, A, s) {
  if (A == 0) return(fn(u, v, B, s))
  integrate(function(d) fn(u, v, d, s), B, B + A, rel.tol = 1e-12)$value / A
}

rsub_cases <- list(c(1.3, 1.7, 0, .8), c(1.3, .9, 1.1, .8), c(.3, 1, .5, .4),
                   c(2, 4, 3, 1), c(2, .05, 1e-3, 3), c(.5, 0, 1, 1),
                   c(1, 2, 1e-6, .3), c(4, .3, .2, 2))

test_that("density and CDF match the closed forms and start-point quadrature", {
  for (packets in c(FALSE, TRUE)) {
    ref <- rsub_ref(packets)
    for (cs in rsub_cases) {
      v <- cs[1]; B <- cs[2]; A <- cs[3]; s <- cs[4]
      med <- (B + A / 2) / v
      info <- paste(c(packets, cs), collapse = "/")
      for (u in med * c(.05, .3, 1, 2, 4)) {
        expect_equal(EMC2:::prsub(u, v, B, A, s, packets),
                     rsub_ref_avg(ref$cdf, u, v, B, A, s),
                     tolerance = 1e-9, info = info)
        expect_equal(EMC2:::drsub(u, v, B, A, s, packets),
                     rsub_ref_avg(ref$pdf, u, v, B, A, s),
                     tolerance = 1e-9, info = info)
        expect_equal(EMC2:::prsub(u, v, B, A, s, packets) +
                       EMC2:::prsub(u, v, B, A, s, packets, lower_tail = FALSE), 1,
                     tolerance = 1e-14, info = info)
      }
    }
  }
})

test_that("the dense point CDF is the inverse-Gaussian survival of the evidence", {
  pig <- function(x, mu, lam) pnorm(sqrt(lam / x) * (x / mu - 1)) +
    exp(2 * lam / mu) * pnorm(-sqrt(lam / x) * (x / mu + 1))
  v <- 1.3; d <- 1.7; s <- .8
  for (u in c(.2, .9, 2.5)) {
    m <- v * u
    expect_equal(EMC2:::prsub(u, v, d, 0, s), 1 - pig(d, m, m^2 / s^2),
                 tolerance = 1e-12)
  }
})

test_that("the packet density is the Bessel closed form", {
  v <- 1.3; d <- 1.7; s <- .8; c <- s^2
  for (u in c(.2, .9, 2.5)) {
    a <- 2 * v * u / c; b <- 2 * d / c
    f <- 2 * v / c * exp(-(sqrt(a) - sqrt(b))^2) *
      besselI(2 * sqrt(a * b), 0, expon.scaled = TRUE)
    expect_equal(EMC2:::drsub(u, v, d, 0, s, TRUE), f, tolerance = 1e-12)
  }
})

test_that("the CDF is the integral of the density", {
  for (packets in c(FALSE, TRUE)) for (cs in rsub_cases[c(1, 2, 4, 6, 8)]) {
    v <- cs[1]; B <- cs[2]; A <- cs[3]; s <- cs[4]
    for (u in c(.1, .7, 2)) {
      Fi <- integrate(function(x) EMC2:::drsub(x, v, B, A, s, packets), 0, u,
                      rel.tol = 1e-12, subdivisions = 1000)$value
      expect_equal(EMC2:::prsub(u, v, B, A, s, packets), Fi, tolerance = 1e-9,
                   info = paste(c(packets, cs, u), collapse = "/"))
    }
  }
})

test_that("the density at zero is positive and matches its limit", {
  v <- .7; d <- 1.2; s <- .9; r <- d / s^2
  f0 <- 2 * v / sqrt(s^2 * d) * dnorm(sqrt(r)) - 2 * v / s^2 * pnorm(-sqrt(r))
  expect_gt(f0, 0)
  expect_equal(EMC2:::drsub(1e-9, v, d, 0, s), f0, tolerance = 1e-6)
  # Packets: the first event arrives at rate 2 v / c and suffices when no
  # further packet is needed.
  expect_equal(EMC2:::drsub(1e-12, v, d, 0, s, TRUE),
               2 * v / s^2 * exp(-2 * d / s^2), tolerance = 1e-9)
})

test_that("packets near zero distance finish at the first event", {
  v <- .7; s <- .9; u <- c(.1, .8, 3)
  expect_equal(EMC2:::prsub(u, v, 1e-12, 0, s, TRUE), pexp(u, 2 * v / s^2),
               tolerance = 1e-10)
  expect_equal(EMC2:::drsub(u, v, 1e-12, 0, s, TRUE), dexp(u, 2 * v / s^2),
               tolerance = 1e-10)
})

test_that("log survivor and density keep their digits deep in the upper tail", {
  # The start-point average is dominated by the far end, where the log
  # integrand falls at slope kappa; integrate over the last 60 / kappa.
  tail_ref <- function(g, B, A) {
    d2 <- B + A; h <- 1e-6 * d2
    kap <- (g(d2) - g(d2 - h)) / h
    w <- if (A == 0) 0 else min(A, 60 / kap)
    l0 <- g(d2)
    if (w == 0) return(l0)
    l0 + log(integrate(function(d) exp(g(d) - l0), d2 - w, d2,
                       rel.tol = 1e-13, subdivisions = 5000)$value / A)
  }
  cases <- list(c(2, 0, 3, .05), c(2, 1, .5, .4), c(.3, 4, 3, 1), c(1, 2, 0, .2))
  for (packets in c(FALSE, TRUE)) for (cs in cases) {
    v <- cs[1]; B <- cs[2]; A <- cs[3]; s <- cs[4]
    for (u in (B + A) / v * c(3, 8)) {
      gS <- function(d) EMC2:::prsub(rep(u, length(d)), v, d, 0, s, packets,
                                     FALSE, TRUE)
      gf <- function(d) EMC2:::drsub(rep(u, length(d)), v, d, 0, s, packets, TRUE)
      info <- paste(c(packets, cs, u), collapse = "/")
      expect_equal(EMC2:::prsub(u, v, B, A, s, packets, FALSE, TRUE),
                   tail_ref(gS, B, A), tolerance = 1e-10, info = info)
      expect_equal(EMC2:::drsub(u, v, B, A, s, packets, TRUE),
                   tail_ref(gf, B, A), tolerance = 1e-10, info = info)
      if (A == 0) {
        ref <- if (packets) rsub_pk_lF(u, v, B, s, upper = FALSE) else {
          z <- sqrt(s^2 * B); a <- (v * u - B) / z; b <- (v * u + B) / z
          lm <- function(x) pnorm(-x, log.p = TRUE) - dnorm(x, log = TRUE)
          dnorm(a, log = TRUE) + lm(a) + log1p(exp(lm(b) - lm(a)))
        }
        expect_equal(EMC2:::prsub(u, v, B, 0, s, packets, FALSE, TRUE), ref,
                     tolerance = 1e-10, info = info)
      }
    }
  }
})

test_that("the log CDF keeps its digits at the leading edge", {
  for (packets in c(FALSE, TRUE))
    for (cs in list(c(1, 1, 0, .3), c(1, 1, .5, .3), c(2, 0, 1, .5))) {
      v <- cs[1]; B <- cs[2]; A <- cs[3]; s <- cs[4]
      for (u in c(.05, .1)) {
        ref <- integrate(function(x) EMC2:::drsub(x, v, B, A, s, packets), 0, u,
                         rel.tol = 1e-12, abs.tol = 0)$value
        expect_equal(EMC2:::prsub(u, v, B, A, s, packets, log_out = TRUE),
                     log(ref), tolerance = 1e-9,
                     info = paste(c(packets, cs, u), collapse = "/"))
      }
    }
})

test_that("the evidence-scale invariance holds", {
  u <- c(.1, .5, 1.5)
  k <- 3.7
  for (packets in c(FALSE, TRUE)) {
    expect_equal(EMC2:::drsub(u, 1.2, .8, .4, .6, packets),
                 EMC2:::drsub(u, k * 1.2, k * .8, k * .4, sqrt(k) * .6, packets),
                 tolerance = 1e-12)
    expect_equal(EMC2:::prsub(u, 1.2, .8, .4, .6, packets),
                 EMC2:::prsub(u, k * 1.2, k * .8, k * .4, sqrt(k) * .6, packets),
                 tolerance = 1e-12)
  }
})

test_that("the two evidence laws converge as packets shrink", {
  # With c -> 0 at fixed v, both first-passage laws approach the ballistic
  # time d / v; relative to their common spread they share the Gaussian limit.
  u <- seq(.85, 1.15, by = .05); v <- 1; d <- 1; s <- .02
  expect_lt(max(abs(EMC2:::prsub(u, v, d, 0, s, TRUE) -
                      EMC2:::prsub(u, v, d, 0, s))), 0.01)
})

test_that("the simulator reproduces the CDF", {
  skip_on_cran()
  local_test_rng(kind = "Mersenne-Twister")
  n <- 20000
  probe <- c(.4, .7, 1, 1.5, 2.5)
  pars <- cbind(v = 1.3, B = .9, A = 1.1, t0 = .2, s = .8)[rep(1, n), ]
  for (packets in c(FALSE, TRUE)) for (cpp in c(TRUE, FALSE)) {
    sim <- withr::with_options(list(emc2.cpp_rfun = cpp),
                               EMC2:::.rfun_RSUB(factor(rep("r", n)), pars, packets))
    th <- EMC2:::prsub(probe - .2, 1.3, .9, 1.1, .8, packets)
    emp <- vapply(probe, function(x) mean(sim$rt <= x), numeric(1))
    expect_lt(max(abs(emp - th)), 0.015)
    expect_true(all(sim$rt > .2))
  }
  # Race: the response probability is the integral of the winner's subdensity.
  lR <- factor(rep(c("a", "b"), n), levels = c("a", "b"))
  p2 <- cbind(v = rep(c(2, .5), n), B = 1, A = .5, t0 = 0, s = 1)
  for (packets in c(FALSE, TRUE)) {
    sim <- EMC2:::.rfun_RSUB(lR, p2, packets)
    pa <- integrate(function(u) EMC2:::drsub(u, 2, 1, .5, 1, packets) *
                      EMC2:::prsub(u, .5, 1, .5, 1, packets, lower_tail = FALSE),
                    0, Inf, rel.tol = 1e-10)$value
    expect_lt(abs(mean(sim$R == "a") - pa), 0.015)
  }
})

rsub_fixture <- function(packets, n = 60) {
  set.seed(20261017)
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  dat <- dat[seq_len(n), ]
  model <- if (packets) function() RSUB("packets") else RSUB
  des <- suppressMessages(design(data = dat, model = model, matchfun = matchfun,
                                 formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1,
                                                s ~ lM)))
  invisible(utils::capture.output(p <- sampled_pars(des)))
  p[] <- 0
  n0 <- names(p)
  p[c("v", "v_lMTRUE", "B", "A", "t0", "s", "s_lMTRUE")] <-
    c(log(1.5), 0.9, log(.8), log(.3), log(.15), log(.7), 0.3)
  stopifnot(identical(names(p), n0))
  list(dat = dat, des = des, p = p)
}

rsub_ll <- function(emc, p) {
  model <- emc[[1]]$model()
  p_types <- names(model$p_types)
  dadm <- emc[[1]]$data[[1]]
  designs <- list()
  for (nm in p_types) {
    designs[[nm]] <- attr(dadm, "designs")[[nm]][
      attr(attr(dadm, "designs")[[nm]], "expand"), , drop = FALSE]
  }
  EMC2:::calc_ll_oo(matrix(p, nrow = 1, dimnames = list(NULL, names(p))), dadm,
                    constants = attr(dadm, "constants"), designs = designs,
                    type = model$c_name, bounds = model$bound,
                    transforms = model$transform,
                    pretransforms = model$pre_transform, p_types = p_types,
                    min_ll = log(1e-10), trend = model$trend)
}

# Independent race likelihood from the reference forms, with an optional
# [LT, UT] truncation normaliser.
rsub_ref_ll <- function(dadm, pars, packets, LT = 0, UT = Inf) {
  ref <- rsub_ref(packets)
  S <- function(t, i) {
    u <- t - pars[i, "t0"]
    if (u <= 0) return(1)
    if (is.infinite(u)) return(0)
    1 - rsub_ref_avg(ref$cdf, u, pars[i, "v"], pars[i, "B"], pars[i, "A"],
                     pars[i, "s"])
  }
  n_lR <- nlevels(dadm$lR)
  trial <- rep(seq_len(nrow(dadm) / n_lR), each = n_lR)
  win <- as.logical(dadm$winner)
  ll <- vapply(split(seq_len(nrow(dadm)), trial), function(ix) {
    w <- ix[win[ix]]
    num <- log(rsub_ref_avg(ref$pdf, dadm$rt[w] - pars[w, "t0"], pars[w, "v"],
                            pars[w, "B"], pars[w, "A"], pars[w, "s"])) +
      sum(log(vapply(ix[!win[ix]], function(i) S(dadm$rt[i], i), 0)))
    Z <- prod(vapply(ix, function(i) S(LT, i), 0)) -
      prod(vapply(ix, function(i) S(UT, i), 0))
    num - log(Z)
  }, numeric(1))
  sum(pmax(ll, log(1e-10)))
}

test_that("the compiled race likelihood matches the reference", {
  skip_on_cran()
  for (packets in c(FALSE, TRUE)) {
    fx <- rsub_fixture(packets)
    emc <- suppressMessages(make_emc(fx$dat, fx$des, type = "single", n_chains = 1,
                                     compress = FALSE, rt_resolution = NULL))
    dadm <- emc[[1]]$data[[1]]
    pars <- EMC2:::get_pars_matrix_oo(fx$p, dadm, emc[[1]]$model())
    expect_equal(rsub_ll(emc, fx$p), rsub_ref_ll(dadm, pars, packets),
                 tolerance = 1e-8, info = packets)
  }
})

test_that("truncation normalises with the race survivor", {
  skip_on_cran()
  LT <- 0.2; UT <- 1.2
  for (packets in c(FALSE, TRUE)) {
    fx <- rsub_fixture(packets)
    d2 <- fx$dat[fx$dat$rt >= LT & fx$dat$rt <= UT, ]
    d2$LT <- LT; d2$UT <- UT
    emc <- suppressMessages(make_emc(d2, fx$des, type = "single", n_chains = 1,
                                     compress = FALSE, rt_resolution = NULL))
    dadm <- emc[[1]]$data[[1]]
    pars <- EMC2:::get_pars_matrix_oo(fx$p, dadm, emc[[1]]$model())
    expect_equal(rsub_ll(emc, fx$p), rsub_ref_ll(dadm, pars, packets, LT, UT),
                 tolerance = 1e-8, info = packets)
  }
})

test_that("narrow start ranges at B near zero keep the CDF digits", {
  # The average of the point CDF over d, in log d, resolves the transition
  # near d ~ m^2 / s^2 that a rule over the whole range misses.
  ref <- function(u, v, B, A, s, pk) {
    f <- function(y) exp(y) * EMC2:::prsub(u, v, exp(y), 0, s, packets = pk)
    lo <- log(max(B, 1e-40)); hi <- log(B + A); m <- .11 * u
    br <- c(lo, log(m^2 / s^2) + seq(-12, 12, 2), log(m) + seq(-4, 4, 2), hi)
    br <- sort(unique(br[br >= lo & br <= hi]))
    I <- sum(vapply(seq_along(br[-1]), function(k)
      integrate(f, br[k], br[k + 1], rel.tol = 1e-14)$value, 0))
    if (B == 0) I <- I + exp(lo)  # F_d -> 1 below exp(lo)
    log(I / A)
  }
  for (pk in c(FALSE, TRUE)) for (B in c(0, 1e-9)) for (A in c(1.5e-5, 1e-4))
    for (u in c(4.2e-7, 1e-5, 2e-4)) {
      lF <- EMC2:::prsub(u, .11, B, A, 2.1, packets = pk, log_out = TRUE)
      expect_equal(lF, ref(u, .11, B, A, 2.1, pk), tolerance = 1e-10,
                   info = paste(pk, B, A, u))
      lS <- EMC2:::prsub(u, .11, B, A, 2.1, packets = pk, lower_tail = FALSE,
                         log_out = TRUE)
      expect_equal(exp(lF) + exp(lS), 1, tolerance = 1e-14)
    }
})

test_that("v = 0 never finishes and invalid rows simulate as omissions", {
  for (pk in c(FALSE, TRUE)) {
    expect_identical(EMC2:::drsub(c(.1, 1), 0, 1, .5, 1, packets = pk), c(0, 0))
    expect_identical(EMC2:::prsub(c(.1, 1), 0, 1, .5, 1, packets = pk), c(0, 0))
    expect_identical(EMC2:::prsub(c(.1, 1), 0, 1, .5, 1, packets = pk,
                                  lower_tail = FALSE), c(1, 1))
    lR <- factor(c("a", "b"), levels = c("a", "b"))
    pars <- cbind(v = c(0, 1), B = c(1, 0), A = c(.5, 0), t0 = .2, s = 1)
    sim_c <- EMC2:::rrsub_cpp(pars[c(1, 1), ], levels(lR), c(TRUE, TRUE), pk)
    expect_true(all(!is.finite(sim_c$rt)))
    # B = A = 0 is not a valid row in either simulator.
    sim_r <- EMC2:::rRSUB(lR, pars[c(2, 2), ], packets = pk)
    sim_c <- EMC2:::rrsub_cpp(pars[c(2, 2), ], levels(lR), c(TRUE, TRUE), pk)
    expect_true(all(!is.finite(sim_r$rt)))
    expect_true(all(!is.finite(sim_c$rt)))
  }
})

test_that("RSUB exposes the race-model contract", {
  m <- RSUB()
  expect_identical(m$type, "RACE")
  expect_identical(m$c_name, "RSUB")
  expect_identical(names(m$p_types)[1:5], m$p_types_canonical)
  expect_identical(RSUB("packets")$c_name, "RSUB_PK")
  expect_identical(names(RSUB("packets")$p_types), names(m$p_types))
})
