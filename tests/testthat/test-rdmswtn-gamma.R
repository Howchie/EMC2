# RDMSWTN(drift_dist = "gamma"): Wald accumulators with Gamma(k, rate tau)
# drift, k = (v/sv)^2, tau = v/sv^2 (src/rdmswtn_gamma.h).

gdrift_wald_d <- function(t, v, b, s) {
  b / (s * sqrt(2 * pi * t^3)) * exp(-(b - v * t)^2 / (2 * s^2 * t))
}
gdrift_wald_p <- function(t, v, b, s) {
  pnorm((v * t - b) / (s * sqrt(t))) +
    exp(2 * v * b / s^2 + pnorm(-(b + v * t) / (s * sqrt(t)), log.p = TRUE))
}

# Reference by quadrature over u = log V across the whole drift law.
gdrift_ref <- function(t, v, b, s, sv, cdf = FALSE) {
  k <- (v / sv)^2; tau <- v / sv^2
  g <- function(u) {
    x <- exp(u)
    w <- exp(dgamma(x, k, tau, log = TRUE) + u)
    w * if (cdf) gdrift_wald_p(t, x, b, s) else gdrift_wald_d(t, x, b, s)
  }
  centre <- log(k / tau); half <- 12 / sqrt(k) + 8 / sqrt(1 + k)
  # P(V < e^lo) <= (tau e^lo)^k / Gamma(k + 1) < 1e-16
  lo <- min(centre - half, (lgamma(k + 1) - 37) / k - log(tau))
  hi <- centre + max(half, log(1 + 40 * b / (s^2 * t)))
  integrate(g, lo, hi, rel.tol = 1e-11, abs.tol = 0, subdivisions = 5000L)$value
}

gdrift_d <- function(t, v, b, s, sv, log = FALSE) {
  n <- max(length(t), length(v), length(b), length(s), length(sv)); r <- function(x) rep(x, length.out = n)
  EMC2:::dRDMSWTN_gdrift_cpp(r(t), r(v), r(b), r(0), r(s), r(sv), log_out = log)
}
gdrift_p <- function(t, v, b, s, sv, lower = TRUE) {
  n <- max(length(t), length(v), length(b), length(s), length(sv)); r <- function(x) rep(x, length.out = n)
  EMC2:::pRDMSWTN_gdrift_cpp(r(t), r(v), r(b), r(0), r(s), r(sv), lower_tail = lower)
}

test_that("gamma-drift density and CDF match quadrature over the drift law", {
  grid <- expand.grid(t = c(0.15, 0.4, 1, 3), v = c(0.7, 2.5), b = c(0.8, 2),
                      s = c(1, 0.7), sv = c(0.2, 0.8, 2))
  d <- gdrift_d(grid$t, grid$v, grid$b, grid$s, grid$sv)
  p <- gdrift_p(grid$t, grid$v, grid$b, grid$s, grid$sv)
  rd <- mapply(gdrift_ref, grid$t, grid$v, grid$b, grid$s, grid$sv)
  rp <- mapply(gdrift_ref, grid$t, grid$v, grid$b, grid$s, grid$sv,
               MoreArgs = list(cdf = TRUE))
  expect_lt(max(abs(d / rd - 1)), 1e-8)
  expect_lt(max(abs(p / rp - 1)), 1e-8)
})

test_that("log CDF and log survivor are complementary and keep tail digits", {
  t <- c(0.02, 0.1, 0.5, 2, 10, 40)
  for (sv in c(0.3, 1.5)) {
    lF <- log(gdrift_p(t, 2, 1.5, 1, sv))
    lS <- EMC2:::pRDMSWTN_gdrift_cpp(t, rep(2, 6), rep(1.5, 6), rep(0, 6),
                                      rep(1, 6), rep(sv, 6), log_out = TRUE,
                                      lower_tail = FALSE)
    expect_true(all(is.finite(lS)))
    expect_lt(max(abs(exp(lF) + exp(lS) - 1)), 1e-14)
  }
  # far right tail: the survivor is resolved, not rounded to zero
  lS <- EMC2:::pRDMSWTN_gdrift_cpp(60, 3, 1, 0, 1, 0.3, log_out = TRUE, lower_tail = FALSE)
  expect_true(is.finite(lS) && lS < -50)
})

test_that("broad drift laws at large distances match an EMG-kernel reference", {
  # S = int G_Y(y) K(y) dy, F = F_W(0) + int Q_Y(y) K(y) dy with
  # K = -dS_W/dy smooth around delta; Y = V sqrt(t) / s ~ Gamma(k, p).  Logs of
  # f and of the smaller tail, the integrands scaled by g_Y(delta).
  ref <- function(t, v, b, sv, s) {
    k <- (v / sv)^2; p <- v / sv^2 * s / sqrt(t); d <- b / (s * sqrt(t))
    lg <- dgamma(d, k, p, log = TRUE)
    lK <- function(y) log(2 * d) + dnorm(y - d, log = TRUE) +
      pnorm(-(y + d), log.p = TRUE) - dnorm(y + d, log = TRUE)
    I <- function(g) integrate(g, max(0, d - 40), d + 40, rel.tol = 1e-12, abs.tol = 0,
                               subdivisions = 4000L, stop.on.error = FALSE)$value
    lf <- log(d / t) + lg + log(I(function(y) exp(dgamma(y, k, p, log = TRUE) - lg + dnorm(y - d, log = TRUE))))
    lower <- pgamma(d, k, p, log.p = TRUE) < pgamma(d, k, p, lower.tail = FALSE, log.p = TRUE)
    lt <- lg + log(I(function(y) exp(pgamma(y, k, p, lower.tail = lower, log.p = TRUE) - lg + lK(y))))
    c(lf = lf, lt = lt, side = if (lower) 3 else 2)
  }
  # the last two rows: small s, a steep drift density at delta (tilted series)
  grid <- data.frame(t = c(2e-4, 1e-3, 3e-3, 5e-4, 0.32, 0.09), v = c(0.3, 1, 2, 0.05, 10, 24.5),
                     b = c(0.5, 1.5, 2, 1, 2.77, 1.65), sv = c(3, 2, 6, 4, 0.09, 0.23),
                     s = c(1, 1, 1, 1, 0.0046, 0.0099))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    r <- ref(g$t, g$v, g$b, g$sv, g$s)
    got <- c(gdrift_d(g$t, g$v, g$b, g$s, g$sv, log = TRUE),
             log(gdrift_p(g$t, g$v, g$b, g$s, g$sv)),
             EMC2:::pRDMSWTN_gdrift_cpp(g$t, g$v, g$b, 0, g$s, g$sv, log_out = TRUE, lower_tail = FALSE))
    expect_lt(abs(got[1] - r[["lf"]]), 1e-9 * max(1, abs(r[["lf"]])))
    expect_lt(abs(got[r[["side"]]] - r[["lt"]]), 1e-9 * max(1, abs(r[["lt"]])))
  }
  # a burn-in corner (tiny mean, huge sv and b) stays finite and fast
  el <- system.time(for (i in 1:50) x <- gdrift_d(0.5266, 0.0010006, 42178.6, 1, 3146.95, log = TRUE))
  expect_true(is.finite(x))
  expect_lt(el[["elapsed"]], 0.5)
})

test_that("sv -> 0 recovers the Wald law smoothly", {
  for (t in c(0.1, 0.6, 3)) {
    wd <- gdrift_wald_d(t, 2, 1.2, 1); ws <- 1 - gdrift_wald_p(t, 2, 1.2, 1)
    e <- sapply(c(1e-3, 1e-4, 1e-5, 0), function(sv)
      c(gdrift_d(t, 2, 1.2, 1, sv) / wd - 1, gdrift_p(t, 2, 1.2, 1, sv, FALSE) / ws - 1))
    expect_lt(max(abs(e[, 3:4])), 1e-8)
    # O(sv^2): a tenfold smaller sv shrinks the gap a hundredfold
    expect_equal(e[1, 1] / e[1, 2], 100, tolerance = 1e-2)
  }
})

test_that("the gamma chart refuses unsupported options", {
  expect_error(RDMSWTN(drift_dist = "gamma", erlang_type = "local_kill"), "clocks")
  expect_error(RDMSWTN(drift_dist = "gamma", erlang_shape = 2L), "clocks")
  expect_error(RDMSWTN(drift_dist = "gamma", posdrift = FALSE), "posdrift")
  expect_error(RDMSWTN(drift_dist = "gamma", correlated = TRUE), "correlated")
  m <- RDMSWTN(drift_dist = "gamma")
  expect_identical(m$c_name, "RDMSWTN_GDRIFT")
  expect_false(any(c("A", "mG", "mK") %in% names(m$p_types)))
})

gdrift_fixture <- function(n = 80, seed = 20261019) {
  set.seed(seed)
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  dat <- dat[seq_len(n), ]
  des <- suppressMessages(design(
    data = dat, model = RDMSWTN(drift_dist = "gamma"), matchfun = matchfun,
    formula = list(v ~ lM, B ~ 1, t0 ~ 1, sv ~ 1),
    constants = c(s = log(1)), report_p_vector = FALSE))
  p <- sampled_pars(des)
  p[] <- c(log(2.5), 0.6, log(1.4), log(0.15), log(0.9))
  emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1,
                                   compress = FALSE, rt_resolution = NULL))
  list(dat = dat, des = des, p = p, emc = emc)
}

gdrift_ll <- function(emc, p) {
  model <- emc[[1]]$model()
  dadm <- emc[[1]]$data[[1]]
  designs <- list()
  for (nm in names(model$p_types)) {
    designs[[nm]] <- attr(dadm, "designs")[[nm]][
      attr(attr(dadm, "designs")[[nm]], "expand"), , drop = FALSE]
  }
  EMC2:::calc_ll_oo(matrix(p, nrow = 1, dimnames = list(NULL, names(p))), dadm,
                    constants = attr(dadm, "constants"), designs = designs,
                    type = model$c_name, bounds = model$bound,
                    transforms = model$transform,
                    pretransforms = model$pre_transform,
                    p_types = names(model$p_types),
                    min_ll = log(1e-10), trend = model$trend)
}

test_that("the compiled race likelihood matches the R reference", {
  skip_on_cran()
  fx <- gdrift_fixture()
  got <- gdrift_ll(fx$emc, fx$p)
  expect_true(is.finite(got))
  dadm <- fx$emc[[1]]$data[[1]]
  pars <- EMC2:::get_pars_matrix_oo(fx$p, dadm, fx$emc[[1]]$model())
  model <- fx$emc[[1]]$model()
  d <- model$dfun(dadm$rt, pars)
  S <- 1 - model$pfun(dadm$rt, pars)
  n_lR <- nlevels(dadm$lR)
  trial <- rep(seq_len(nrow(dadm) / n_lR), each = n_lR)
  win <- as.logical(dadm$winner)
  ref <- sum(vapply(split(seq_len(nrow(dadm)), trial), function(ix)
    max(log(d[ix][win[ix]]) + sum(log(S[ix][!win[ix]])), log(1e-10)), numeric(1)))
  expect_equal(got, ref, tolerance = 1e-8)
})

test_that("simulated finishing times follow the gamma-drift CDF", {
  skip_on_cran()
  set.seed(7)
  lR <- factor(rep("a", 4000), levels = "a")
  pars <- cbind(v = 1.8, B = 1.3, t0 = 0.2, s = 1, sv = 1.1)[rep(1, 4000), ]
  sim <- EMC2:::rRDMSWTN_gdrift(lR, pars)
  qs <- quantile(sim$rt, c(0.1, 0.3, 0.5, 0.7, 0.9))
  expect_equal(unname(gdrift_p(qs - 0.2, 1.8, 1.3, 1, 1.1)),
               c(0.1, 0.3, 0.5, 0.7, 0.9), tolerance = 0.03)
})
