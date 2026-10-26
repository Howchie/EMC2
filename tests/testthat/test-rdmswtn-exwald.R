# RDMSWTN(exwald = TRUE): Wald accumulators plus an Exp(mean tau) hesitation,
# start point uniform on [0, A], sv = 0 (src/model_RDMSWTN_exwald.cpp).

exw_d <- function(t, v, B, A, s, tau, log = FALSE) {
  n <- max(length(t), length(v), length(B), length(A), length(s), length(tau))
  r <- function(x) rep(x, length.out = n)
  EMC2:::dRDMSWTN_exwald_cpp(r(t), r(v), r(B), r(A), r(0), r(s), r(tau), log_out = log)
}
exw_p <- function(t, v, B, A, s, tau, lower = TRUE, log = FALSE) {
  n <- max(length(t), length(v), length(B), length(A), length(s), length(tau))
  r <- function(x) rep(x, length.out = n)
  EMC2:::pRDMSWTN_exwald_cpp(r(t), r(v), r(B), r(A), r(0), r(s), r(tau),
                             log_out = log, lower_tail = lower)
}

# Convolution of the Wald density with the exponential, in log space about
# its peak, averaged over the distance by adaptive quadrature.
exw_ref <- function(t, v, B, A, s, tau) {
  lam <- 1 / tau
  point <- function(a) {
    lg <- function(u) log(a) - 0.5 * log(2 * pi * s^2 * u^3) -
      (a - v * u)^2 / (2 * s^2 * u) + log(lam) - lam * (t - u)
    m <- optimize(lg, c(1e-9, t), maximum = TRUE)$objective
    m + log(integrate(function(u) exp(lg(u) - m), 0, t, rel.tol = 1e-13,
                      subdivisions = 5000L)$value)
  }
  if (A == 0) return(point(B))
  l0 <- point(B + A / 2)
  l0 + log(integrate(function(a) exp(vapply(a, point, numeric(1)) - l0), B, B + A,
                     rel.tol = 1e-12)$value / A)
}

test_that("the ex-Wald density matches the convolution on both kappa branches", {
  grid <- expand.grid(t = c(0.08, 0.4, 1.5), v = c(0.6, 3), B = c(0.3, 1),
                      A = c(0, 0.7), s = c(1, 0.6), tau = c(0.02, 0.3, 2))
  # kappa^2 = v^2 - 2 s^2 / tau: both signs occur in the grid
  k2 <- grid$v^2 - 2 * grid$s^2 / grid$tau
  expect_true(any(k2 > 0) && any(k2 < 0))
  got <- exw_d(grid$t, grid$v, grid$B, grid$A, grid$s, grid$tau, log = TRUE)
  ref <- mapply(exw_ref, grid$t, grid$v, grid$B, grid$A, grid$s, grid$tau)
  expect_equal(got, ref, tolerance = 1e-8)
})

test_that("CDF integrates the density and S = 1 - F", {
  for (pp in list(c(0.6, 1, 0.5, 1, 0.3), c(3, 0.4, 0.8, 0.7, 2), c(1, 1, 0, 1, 0.05))) {
    for (t in c(0.15, 0.7, 2)) {
      F <- integrate(function(x) exw_d(x, pp[1], pp[2], pp[3], pp[4], pp[5]), 0, t,
                     rel.tol = 1e-11)$value
      expect_equal(exw_p(t, pp[1], pp[2], pp[3], pp[4], pp[5]), F, tolerance = 1e-8)
      expect_equal(exw_p(t, pp[1], pp[2], pp[3], pp[4], pp[5]) +
                   exw_p(t, pp[1], pp[2], pp[3], pp[4], pp[5], lower = FALSE), 1,
                   tolerance = 1e-12)
    }
  }
})

test_that("tau -> 0 recovers the Wald and huge tau keeps the CDF's digits", {
  wald <- exw_d(0.5, 2, 0.5, 0.5, 1, 0)
  expect_equal(exw_d(0.5, 2, 0.5, 0.5, 1, 1e-8), wald, tolerance = 1e-6)
  expect_equal(exw_p(0.5, 2, 0.5, 0.5, 1, 1e-8), exw_p(0.5, 2, 0.5, 0.5, 1, 0),
               tolerance = 1e-6)
  # F_X ~ lambda int_0^t (t - u) f_D(u) du as lambda t -> 0
  approx <- 1e-7 * integrate(function(u) (0.3 - u) * exw_d(u, 2, 0.5, 0.5, 1, 0),
                             0, 0.3, rel.tol = 1e-12)$value
  expect_equal(exw_p(0.3, 2, 0.5, 0.5, 1, 1e7), approx, tolerance = 1e-5)
})

test_that("log density and log survivor stay finite in the tails", {
  ld <- exw_d(c(1e-3, 30), 2, 0.5, 0.5, 1, 0.2, log = TRUE)
  expect_true(all(is.finite(ld)))
  # lambda > v^2 / 2: the Wald tail (rate v^2 / 2) outlasts the hesitation
  ls <- exw_p(c(100, 200), 2, 0.5, 0.5, 1, 0.2, lower = FALSE, log = TRUE)
  expect_equal(diff(ls) / 100, -2, tolerance = 1e-2)
})

test_that("the ex-Wald toggle refuses unsupported options", {
  expect_error(RDMSWTN(exwald = TRUE, drift_dist = "gamma"), "drift variability")
  expect_error(RDMSWTN(exwald = TRUE, erlang_type = "local_kill"), "clocks")
  expect_error(RDMSWTN(exwald = TRUE, erlang_shape = 2L), "clocks")
  expect_error(RDMSWTN(exwald = TRUE, posdrift = FALSE), "posdrift")
  expect_error(RDMSWTN(exwald = TRUE, correlated = TRUE), "correlated")
  m <- RDMSWTN(exwald = TRUE)
  expect_identical(m$c_name, "RDMSWTN_EXWALD")
  expect_true("tau" %in% names(m$p_types))
  expect_false(any(c("sv", "mG", "mK") %in% names(m$p_types)))
})

exw_fixture <- function(n = 80, seed = 20261026) {
  set.seed(seed)
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  dat <- dat[seq_len(n), ]
  des <- suppressMessages(design(
    data = dat, model = RDMSWTN(exwald = TRUE), matchfun = matchfun,
    formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, tau ~ 1),
    constants = c(s = log(1)), report_p_vector = FALSE))
  p <- sampled_pars(des)
  p[c("v", "v_lMTRUE", "B", "A", "t0", "tau")] <-
    c(log(2.5), 0.6, log(1.0), log(0.6), log(0.15), log(0.25))
  emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1,
                                   compress = FALSE, rt_resolution = NULL))
  list(dat = dat, des = des, p = p, emc = emc)
}

test_that("the compiled race likelihood matches the R reference", {
  skip_on_cran()
  fx <- exw_fixture()
  model <- fx$emc[[1]]$model()
  dadm <- fx$emc[[1]]$data[[1]]
  designs <- list()
  for (nm in names(model$p_types)) {
    designs[[nm]] <- attr(dadm, "designs")[[nm]][
      attr(attr(dadm, "designs")[[nm]], "expand"), , drop = FALSE]
  }
  got <- EMC2:::calc_ll_oo(matrix(fx$p, nrow = 1, dimnames = list(NULL, names(fx$p))),
                           dadm, constants = attr(dadm, "constants"), designs = designs,
                           type = model$c_name, bounds = model$bound,
                           transforms = model$transform,
                           pretransforms = model$pre_transform,
                           p_types = names(model$p_types),
                           min_ll = log(1e-10), trend = model$trend)
  expect_true(is.finite(got))
  pars <- EMC2:::get_pars_matrix_oo(fx$p, dadm, model)
  d <- model$dfun(dadm$rt, pars)
  S <- 1 - model$pfun(dadm$rt, pars)
  n_lR <- nlevels(dadm$lR)
  trial <- rep(seq_len(nrow(dadm) / n_lR), each = n_lR)
  win <- as.logical(dadm$winner)
  ref <- sum(vapply(split(seq_len(nrow(dadm)), trial), function(ix)
    max(log(d[ix][win[ix]]) + sum(log(S[ix][!win[ix]])), log(1e-10)), numeric(1)))
  expect_equal(got, ref, tolerance = 1e-8)
})

test_that("simulated finishing times follow the ex-Wald CDF", {
  skip_on_cran()
  set.seed(11)
  lR <- factor(rep("a", 4000), levels = "a")
  probs <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  for (A in c(0, 1.2)) {
    pars <- cbind(v = 1.8, B = 0.6, A = A, t0 = 0.2, s = 1, tau = 0.4)[rep(1, 4000), ]
    sim <- EMC2:::rRDMSWTN_exwald(lR, pars)
    qs <- quantile(sim$rt, probs)
    expect_equal(exw_p(unname(qs) - 0.2, 1.8, 0.6, A, 1, 0.4), probs, tolerance = 0.03)
  }
})

test_that("the ex-Wald stays exact when v / s is very large", {
  # s -> 0 makes D uniform on [B / v, (B + A) / v], so X = D + Exp(tau) has a
  # closed form; tau -> 0 must then reproduce the Wald.
  v <- 3.2e7; B <- 1.5e6; A <- 1.04e7; t <- c(0.17, 0.3, 0.8)
  lo <- B / v; hi <- (B + A) / v; w <- hi - lo
  for (s in c(1e-3, 1.8e-16)) {
    for (tau in c(0.1, 1e8)) {
      lam <- 1 / tau
      f_ref <- (exp(-lam * pmax(t - hi, 0)) - exp(-lam * pmax(t - lo, 0))) / w
      Fe <- function(u) ifelse(u > 0, u + expm1(-lam * u) / lam, 0)
      F_ref <- (Fe(t - lo) - Fe(t - hi)) / w
      expect_equal(exw_d(t, v, B, A, s, tau), f_ref, tolerance = 1e-6)
      expect_equal(exw_p(t, v, B, A, s, tau), F_ref, tolerance = 1e-6)
    }
    for (tau in c(1e-23, 1e-17)) {
      expect_equal(exw_d(t, v, B, A, s, tau), exw_d(t, v, B, A, s, 0), tolerance = 1e-8)
      expect_equal(exw_p(t, v, B, A, s, tau), exw_p(t, v, B, A, s, 0), tolerance = 1e-8)
    }
  }
})

test_that("tiny and huge tau stay on their limits for point and uniform starts", {
  t <- c(0.2, 0.5, 1.2)
  for (A in c(0, 0.5)) {
    # the hesitation shifts the Wald by about tau
    for (tau in c(1e-20, 1e-12, 1e-6)) {
      expect_equal(exw_d(t, 2, 0.5, A, 1, tau), exw_d(t, 2, 0.5, A, 1, 0), tolerance = 50 * tau + 1e-10)
      expect_equal(exw_p(t, 2, 0.5, A, 1, tau), exw_p(t, 2, 0.5, A, 1, 0), tolerance = 50 * tau + 1e-10)
    }
    for (tau in c(1e4, 1e8)) {
      ref <- vapply(t, exw_ref, numeric(1), v = 2, B = 0.5, A = A, s = 1, tau = tau)
      expect_equal(exw_d(t, 2, 0.5, A, 1, tau, log = TRUE), ref, tolerance = 1e-7)
    }
  }
})

test_that("near-deterministic accumulators with tiny tau give the Wald limit", {
  pars <- cbind(v = c(32478964.758, 26158550.348, 32478461.050, 32459885.039),
                B = c(1626093.468, 1596038.667, 1493455.266, 1193807.148),
                A = 10430218.649, s = 1.8032268788e-16,
                tau = c(2.165e-23, 3.699e-17, 1.037e-21, 4.964e-20))
  t <- c(1.5916667, 1.375, 0.9583333, 0.8583333) - 0.788365056
  d <- exw_d(t, pars[, "v"], pars[, "B"], pars[, "A"], pars[, "s"], pars[, "tau"])
  p <- exw_p(t, pars[, "v"], pars[, "B"], pars[, "A"], pars[, "s"], pars[, "tau"])
  expect_equal(d, exw_d(t, pars[, "v"], pars[, "B"], pars[, "A"], pars[, "s"], 0))
  expect_equal(p, exw_p(t, pars[, "v"], pars[, "B"], pars[, "A"], pars[, "s"], 0))
  expect_true(all(is.finite(d) & d >= 0 & p >= 0 & p <= 1))
})
