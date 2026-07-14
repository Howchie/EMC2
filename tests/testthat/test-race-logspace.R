# Numerical-stability tests for the race-model log-space fallbacks
# (RDMSWTN / RDMGBM / SWTN / truncated ex-Gaussian + the shared normal-tail
# helpers).  Throughout, log_out must remain a pure output-scale choice:
# natural-scale results are unchanged in ordinary regimes, and log outputs
# stay finite past the natural under/overflow points.

test_that("RDMSWTN k0 CDF: natural and log outputs agree in the safe region", {
  for (tt in c(0.3, 0.8, 2)) {
    for (A in c(0, 0.4)) {
      p <- EMC2:::pwald(tt, 1.2, 1.5, A, 1, 0)
      lp <- EMC2:::pwald(tt, 1.2, 1.5, A, 1, 0, log_out = TRUE)
      expect_equal(lp, log(p), tolerance = 1e-10)
    }
  }
})

test_that("RDMSWTN k0 log CDF stays finite past natural underflow", {
  # Early time with a distant threshold: the natural CDF underflows to zero,
  # the log CDF must stay finite (point start and uniform start).
  for (A in c(0, 0.3)) {
    p <- EMC2:::pwald(0.05, 1, 15, A, 1, 0)
    lp <- EMC2:::pwald(0.05, 1, 15, A, 1, 0, log_out = TRUE)
    expect_equal(p, 0)
    expect_true(is.finite(lp))
    expect_lt(lp, log(.Machine$double.xmin))
  }
})

test_that("RDMSWTN A>0 log-CDF fallback matches an adaptive reference in deep tails", {
  # The closed-form fallback splits the start-point average into its two
  # positive components; it should agree with an independent adaptive
  # log-space quadrature to near machine precision, including mu <= 0.
  ref_adaptive <- function(t, mu, b, A) {
    dlo <- max(0, b - A); dhi <- b
    lF <- function(d) {
      a1 <- pnorm((mu * t - d) / sqrt(t), log.p = TRUE)
      a2 <- 2 * mu * d + pnorm((-mu * t - d) / sqrt(t), log.p = TRUE)
      m <- pmax(a1, a2)
      m + log1p(exp(pmin(a1, a2) - m))
    }
    m <- max(lF(seq(dlo, dhi, length.out = 101)))
    v <- integrate(function(d) exp(lF(d) - m), dlo, dhi, rel.tol = 1e-12)$value
    m + log(v) - log(A)
  }
  for (cs in list(c(0.05, 1, 15, 0.3), c(0.02, 2, 8, 0.5),
                  c(0.1, -3, 6, 0.4), c(0.03, 0, 10, 0.5))) {
    got <- EMC2:::pwald(cs[1], cs[2], cs[3], cs[4], 1, 0,
                        log_out = TRUE, posdrift = FALSE)
    expect_equal(got, ref_adaptive(cs[1], cs[2], cs[3], cs[4]),
                 tolerance = 1e-8)
  }
})

test_that("RDMSWTN k0 density: log output finite in the early-time tail", {
  for (tt in c(0.3, 1, 3)) {
    for (A in c(0, 0.4)) {
      d <- EMC2:::dwald(tt, 1.2, 1.5, A, 1, 0)
      ld <- EMC2:::dwald(tt, 1.2, 1.5, A, 1, 0, log_out = TRUE)
      expect_equal(ld, log(d), tolerance = 1e-10)
    }
  }
  for (A in c(0, 0.3)) {
    ld <- EMC2:::dwald(0.05, 1, 15, A, 1, 0, log_out = TRUE)
    expect_true(is.finite(ld))
    expect_lt(ld, log(.Machine$double.xmin))
  }
})

test_that("killed Wald eventual-hit mass survives exp underflow (Erlang-2, t = Inf)", {
  # Large kill rates push exp(eta1 * d) below the natural floor; the log
  # output must stay finite instead of collapsing to -Inf.
  p_nat <- EMC2:::pwald(Inf, 1, 2, 0.5, 1, 0, 0, 1e6, FALSE, 2, FALSE, FALSE)
  lp <- EMC2:::pwald(Inf, 1, 2, 0.5, 1, 0, 0, 1e6, TRUE, 2, FALSE, FALSE)
  expect_equal(p_nat, 0)
  expect_true(is.finite(lp))
  expect_lt(lp, log(.Machine$double.xmin))
  # Ordinary kill rates keep their natural/log agreement.
  p <- EMC2:::pwald(Inf, 1, 1.5, 0.4, 1, 0, 0, 0.8, FALSE, 2, FALSE, FALSE)
  lp2 <- EMC2:::pwald(Inf, 1, 1.5, 0.4, 1, 0, 0, 0.8, TRUE, 2, FALSE, FALSE)
  expect_equal(lp2, log(p), tolerance = 1e-10)
})

test_that("SWTN positive-drift densities keep finite log tails", {
  # Central agreement (sv > 0, posdrift).
  for (tt in c(0.4, 1, 2.5)) {
    d <- EMC2:::dswtn(tt, 1, 1.5, 1, 0, 0.5)
    ld <- EMC2:::dswtn(tt, 1, 1.5, 1, 0, 0.5, log_out = TRUE)
    expect_equal(ld, log(d), tolerance = 1e-10)
    p <- EMC2:::pswtn(tt, 1, 1.5, 1, 0, 0.5)
    lp <- EMC2:::pswtn(tt, 1, 1.5, 1, 0, 0.5, log_out = TRUE)
    expect_equal(lp, log(p), tolerance = 1e-10)
  }
  # Early-time tail: natural zero, finite log.
  ld <- EMC2:::dswtn(0.05, 1, 15, 1, 0, 0.5, log_out = TRUE)
  expect_true(is.finite(ld))
  expect_lt(ld, log(.Machine$double.xmin))
  lp <- EMC2:::pswtn(0.05, 1, 15, 1, 0, 0.5, log_out = TRUE)
  expect_true(is.finite(lp))
  expect_lt(lp, log(.Machine$double.xmin))
})

test_that("RDMSWTN joint A/sv density falls back to the log quadrature", {
  # Central agreement first.
  for (tt in c(0.4, 1, 2.5)) {
    d <- EMC2:::drdmswtn(tt, 1, 1.5, 0.4, 1, 0, 0.5)
    ld <- EMC2:::drdmswtn(tt, 1, 1.5, 0.4, 1, 0, 0.5, log_out = TRUE)
    expect_equal(ld, log(d), tolerance = 1e-10)
    p <- EMC2:::prdmswtn(tt, 1, 1.5, 0.4, 1, 0, 0.5)
    lp <- EMC2:::prdmswtn(tt, 1, 1.5, 0.4, 1, 0, 0.5, log_out = TRUE)
    expect_equal(lp, log(p), tolerance = 1e-8)
  }
  # Early-time tail: the closed-form bivariate expressions cancel to zero on
  # the natural scale; the log fallback must produce a finite value.
  ld <- EMC2:::drdmswtn(0.05, 1, 15, 0.4, 1, 0, 0.5, log_out = TRUE)
  expect_true(is.finite(ld))
  expect_lt(ld, log(.Machine$double.xmin))
  lp <- EMC2:::prdmswtn(0.05, 1, 15, 0.4, 1, 0, 0.5, log_out = TRUE)
  expect_true(is.finite(lp))
  expect_lt(lp, log(.Machine$double.xmin))
})

test_that("RDMGBM SPV density: natural fast path and log fallback agree", {
  # Central agreement in the SPV branch (A > 0).
  for (tt in c(0.3, 0.8, 2)) {
    d <- EMC2:::dgbm(tt, 1, 3, 0.5, 1)
    ld <- EMC2:::dgbm(tt, 1, 3, 0.5, 1, log_out = TRUE)
    expect_equal(ld, log(d), tolerance = 1e-10)
    p <- EMC2:::pgbm(tt, 1, 3, 0.5, 1, 0, 0, 0, FALSE, 1L, FALSE, 1)
    lp <- EMC2:::pgbm(tt, 1, 3, 0.5, 1, 0, 0, 0, TRUE, 1L, FALSE, 1)
    expect_equal(lp, log(p), tolerance = 1e-10)
  }
  # Early time, distant boundary: exp(exp_factor) * integral underflows
  # naturally; the signed-log evaluation must stay finite.
  d0 <- EMC2:::dgbm(0.005, 1, 100, 0.5, 1)
  ld0 <- EMC2:::dgbm(0.005, 1, 100, 0.5, 1, log_out = TRUE)
  expect_equal(d0, 0)
  expect_true(is.finite(ld0))
  expect_lt(ld0, log(.Machine$double.xmin))
})

test_that("RDMGBM eventual-hit mass survives prefactor underflow (Erlang-2, t = Inf)", {
  lp <- EMC2:::pgbm(Inf, 0.1, 50, 0.5, 1, 0, 0, 1e6, TRUE, 2L, FALSE, 1)
  expect_true(is.finite(lp))
  expect_lt(lp, log(.Machine$double.xmin))
  # Ordinary rates: natural/log agreement.
  p <- EMC2:::pgbm(Inf, 0.5, 3, 0.5, 1, 0, 0, 0.8, FALSE, 2L, FALSE, 1)
  lp2 <- EMC2:::pgbm(Inf, 0.5, 3, 0.5, 1, 0, 0, 0.8, TRUE, 2L, FALSE, 1)
  expect_equal(lp2, log(p), tolerance = 1e-10)
})

test_that("ex-Gaussian upper tail has a direct log implementation", {
  mu <- 0.4; sigma <- 0.05; tau <- 0.1
  # Central agreement.
  for (q in c(0.3, 0.5, 0.8)) {
    s <- EMC2:::pTEXG_vec(q, mu, sigma, tau, lb = -Inf,
                          lower_tail = FALSE, log_p = FALSE)
    ls <- EMC2:::pTEXG_vec(q, mu, sigma, tau, lb = -Inf,
                           lower_tail = FALSE, log_p = TRUE)
    expect_equal(ls, log(s), tolerance = 1e-8)
  }
  # Saturated lower tail: previously -Inf / 0, now the exponential tail.
  ls5 <- EMC2:::pTEXG_vec(5, mu, sigma, tau, lb = -Inf,
                          lower_tail = FALSE, log_p = TRUE)
  expect_true(is.finite(ls5))
  # Reference: S(q) ~ exp((mu - q)/tau + sigma^2/(2 tau^2)) for q >> mu.
  ref <- (mu - 5) / tau + sigma^2 / (2 * tau^2)
  expect_equal(ls5, ref, tolerance = 1e-6)
  s5 <- EMC2:::pTEXG_vec(5, mu, sigma, tau, lb = -Inf,
                         lower_tail = FALSE, log_p = FALSE)
  expect_equal(s5, exp(ls5), tolerance = 1e-12)
})

test_that("truncated ex-Gaussian survives a far-upper-tail truncation bound", {
  mu <- 0.4; sigma <- 0.05; tau <- 0.1; lb <- 5
  # The truncation normaliser 1 - F(5) is ~1e-20: the natural subtraction
  # collapses, the survivor-log normaliser must not.  The renormalised
  # density must integrate to one over the retained region.
  f <- function(x) EMC2:::dTEXG_vec(x, mu, sigma, tau, lb)
  expect_true(all(is.finite(log(f(c(5.05, 5.2, 5.5))))))
  expect_equal(integrate(f, lb, 9)$value, 1, tolerance = 1e-4)
  # Truncated upper-tail probabilities stay consistent with the density.
  lsurv <- EMC2:::pTEXG_vec(5.2, mu, sigma, tau, lb,
                            lower_tail = FALSE, log_p = TRUE)
  expect_true(is.finite(lsurv))
  expect_equal(exp(lsurv), integrate(f, 5.2, 9)$value, tolerance = 1e-4)
})

test_that("SWTNspv vector wrappers keep finite log tails with sv > 0", {
  ld <- EMC2:::dSWTNspv(t = c(0.5, 0.05), v = c(1, 1), b = c(1.5, 15),
                        A = c(0.4, 0.4), s = 1, t0 = 0, sv = 0.5,
                        log_out = TRUE)
  d <- EMC2:::dSWTNspv(t = c(0.5, 0.05), v = c(1, 1), b = c(1.5, 15),
                       A = c(0.4, 0.4), s = 1, t0 = 0, sv = 0.5,
                       log_out = FALSE)
  expect_equal(ld[1], log(d[1]), tolerance = 1e-10)
  expect_true(all(is.finite(ld)))
  expect_lt(ld[2], log(.Machine$double.xmin))
})

test_that("REXG raw survivor keeps its tail after CDF saturation", {
  # Likelihood-level: a two-choice REXG race where the loser has nearly
  # always finished by the observed RT.  The raw survivor kernel previously
  # collapsed the trial to the min_ll floor via 1 - CDF; it must now yield a
  # finite log-likelihood without NaNs in either regime.
  lI_fun <- function(d) factor(rep(1, nrow(d)), levels = 1)
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  design_rexg <- design(
    data = dat, model = REXG, matchfun = matchfun,
    formula = list(mu ~ lM, sigma ~ 1, tau ~ 1),
    functions = list(lI = lI_fun)
  )
  emc <- make_emc(dat, design_rexg, type = "single", n_chains = 1, compress = FALSE)
  dadm <- emc[[1]]$data[[1]]
  model <- emc[[1]]$model()
  p_types <- names(model$p_types)
  designs <- list()
  for (nm in p_types) {
    designs[[nm]] <- attr(dadm, "designs")[[nm]][attr(attr(dadm, "designs")[[nm]], "expand"), , drop = FALSE]
  }
  ll_for <- function(p) {
    EMC2:::calc_ll_oo(
      matrix(p, nrow = 1, dimnames = list(NULL, p_types)),
      dadm, constants = attr(dadm, "constants"), designs = designs,
      type = model$c_name, bounds = model$bound,
      transforms = model$transform, pretransforms = model$pre_transform,
      p_types = p_types, min_ll = log(1e-10), trend = model$trend
    )
  }
  # Ordinary parameters and an extreme mu contrast (loser far in its tail).
  for (mu_d in c(0.2, 5)) {
    p <- c(0.4, mu_d, log(0.05), log(0.1))
    ll <- ll_for(p)
    expect_true(is.finite(ll))
  }
})
