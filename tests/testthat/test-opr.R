# Opponent Poisson race (Math/Opponent_Poisson_Race.md).  References are the
# note's closed forms evaluated with besselI() and pchisq(), and quadrature of
# the density.

opr_ref_pdf <- function(t, K, a, b) {
  if (b == 0) return(dgamma(t, K, a))
  z <- 2 * t * sqrt(a * b)
  exp(log(K) - log(t) + 0.5 * K * (log(a) - log(b)) - (sqrt(a) - sqrt(b))^2 * t +
        log(besselI(z, nu = K, expon.scaled = TRUE)))
}
opr_ref_cdf <- function(t, K, a, b) {
  if (b == 0) return(pgamma(t, K, a))
  pchisq(2 * a * t, df = 2 * K, ncp = 2 * b * t) +
    exp(K * (log(a) - log(b)) +
          pchisq(2 * b * t, df = 2 * (K + 1), ncp = 2 * a * t, log.p = TRUE))
}

opr_cases <- list(c(3, 30, 10), c(1, 15, 5), c(2.5, 8, 8), c(7.3, 40, 35),
                  c(4, 10, 20), c(0.6, 5, 1), c(12.2, 100, 50), c(3, 30, 0),
                  c(3, 30, 1e-6), c(1.5, 3, 0.2))

test_that("density and CDF match the Bessel / noncentral chi-squared forms", {
  u <- c(.01, .05, .1, .3, .5, 1, 2, 4)
  for (cs in opr_cases) {
    K <- cs[1]; a <- cs[2]; b <- cs[3]
    d <- EMC2:::dopr(u, K, a, b)
    p <- EMC2:::popr(u, K, a, b)
    s <- EMC2:::popr(u, K, a, b, lower_tail = FALSE)
    info <- paste(cs, collapse = "/")
    expect_equal(d, sapply(u, opr_ref_pdf, K = K, a = a, b = b),
                 tolerance = 1e-11, info = info)
    expect_equal(p, sapply(u, opr_ref_cdf, K = K, a = a, b = b),
                 tolerance = 1e-11, info = info)
    expect_equal(p + s, rep(1, length(u)), tolerance = 1e-14, info = info)
    h <- if (a >= b) 1 else (a / b)^K
    expect_equal(EMC2:::popr(Inf, K, a, b), h, tolerance = 1e-14, info = info)
    expect_equal(EMC2:::popr(Inf, K, a, b, lower_tail = FALSE), 1 - h,
                 tolerance = 1e-14, info = info)
  }
})

test_that("the CDF is the integral of the density", {
  for (cs in opr_cases[c(1, 3, 5, 6, 7)]) {
    K <- cs[1]; a <- cs[2]; b <- cs[3]
    for (u in c(.1, .7, 2)) {
      Fi <- integrate(function(s) EMC2:::dopr(s, K, a, b), 0, u,
                      rel.tol = 1e-12, subdivisions = 1000)$value
      expect_equal(EMC2:::popr(u, K, a, b), Fi, tolerance = 1e-10)
    }
  }
})

test_that("the log survivor keeps its digits deep in the upper tail", {
  tail_ref <- function(u, K, a, b) {
    lf0 <- EMC2:::dopr(u, K, a, b, log_out = TRUE)
    lf0 + log(integrate(function(s) exp(EMC2:::dopr(s, K, a, b, log_out = TRUE) - lf0),
                        u, Inf, rel.tol = 1e-13, subdivisions = 2000)$value)
  }
  for (cs in list(c(3, 30, 10), c(12.2, 100, 50), c(0.6, 5, 1), c(1, 1000, 999))) {
    K <- cs[1]; a <- cs[2]; b <- cs[3]
    for (u in K / (a - b) * c(3, 30, 100)) {
      ls <- EMC2:::popr(u, K, a, b, lower_tail = FALSE, log_out = TRUE)
      expect_lt(abs(expm1(ls - tail_ref(u, K, a, b))), 1e-9)
    }
  }
})

test_that("the log CDF keeps its digits at the leading edge", {
  for (cs in list(c(2.5, 3, 1), c(1, 15, 5), c(4, 10, 20), c(1.3, 10, 8))) {
    K <- cs[1]; a <- cs[2]; b <- cs[3]
    for (u in c(1e-3, 1e-6, 1e-9)) {
      ref <- integrate(function(x) EMC2:::dopr(x, K, a, b), 0, u,
                       rel.tol = 1e-12, abs.tol = 0)$value
      expect_equal(EMC2:::popr(u, K, a, b, log_out = TRUE), log(ref),
                   tolerance = 1e-9, info = paste(c(cs, u), collapse = "/"))
    }
  }
})

test_that("defective accumulators satisfy f(u; a, b) = h f(u; b, a)", {
  u <- c(.05, .3, 1, 3)
  K <- 2.7; a <- 6; b <- 11
  h <- (a / b)^K
  expect_equal(EMC2:::dopr(u, K, a, b), h * EMC2:::dopr(u, K, b, a),
               tolerance = 1e-13)
  expect_equal(EMC2:::popr(u, K, a, b), h * EMC2:::popr(u, K, b, a),
               tolerance = 1e-13)
})

test_that("moments match K / (a - b) and K (a + b) / (a - b)^3", {
  K <- 3.4; a <- 20; b <- 8
  S <- function(u) EMC2:::popr(u, K, a, b, lower_tail = FALSE)
  m1 <- integrate(S, 0, Inf, rel.tol = 1e-12)$value
  m2 <- 2 * integrate(function(u) u * S(u), 0, Inf, rel.tol = 1e-12)$value
  expect_equal(m1, K / (a - b), tolerance = 1e-9)
  expect_equal(m2 - m1^2, K * (a + b) / (a - b)^3, tolerance = 1e-8)
})

test_that("the simulator reproduces the signed Poisson walk and the CDF", {
  skip_on_cran()
  local_test_rng(kind = "Mersenne-Twister")
  n <- 4000
  # Integer criterion: first passage of N+ - N- to K by event simulation.
  walk <- function(K, a, b) {
    x <- 0; t <- 0
    while (x < K) {
      t <- t + rexp(1, a + b)
      x <- x + if (runif(1) < a / (a + b)) 1 else -1
    }
    t
  }
  lR <- factor(rep("r", n))
  pars <- cbind(K = 3, lambda_p = 12, lambda_m = 4, t0 = 0)[rep(1, n), ]
  sim <- EMC2:::ropr_cpp(pars, levels(lR), rep(TRUE, n), 1L)$rt
  ref <- replicate(n, walk(3, 12, 4))
  expect_gt(suppressWarnings(ks.test(sim, ref))$p.value, 1e-3)

  # Noninteger criterion and both charts, against the analytic CDF.
  probe <- c(.1, .2, .4, .8)
  pars_v <- cbind(K = 2.6, v = 8, rho = 0.7, t0 = 0.1)[rep(1, 2 * n), ]
  lR2 <- factor(rep(c("a", "b"), n), levels = c("a", "b"))
  for (cpp in c(TRUE, FALSE)) {
    sim2 <- withr::with_options(list(emc2.cpp_rfun = cpp),
                                EMC2:::.rfun_OPR(lR2, pars_v, 0L))
    th <- 1 - EMC2:::popr(probe - .1, 2.6, 8 * 1.7, 8 * .7, lower_tail = FALSE)^2
    emp <- vapply(probe, function(x) mean(sim2$rt <= x), numeric(1))
    expect_lt(max(abs(emp - th)), 0.03)
    expect_true(all(sim2$rt > .1))
  }

  # Background chart: lambda_p = v + lambda0, lambda_m = lambda0.
  pars_b <- cbind(K = 2.6, v = 8, lambda0 = 5.6, t0 = 0.1)[rep(1, 2 * n), ]
  simb <- EMC2:::.rfun_OPR(lR2, pars_b, 2L)
  emp <- vapply(probe, function(x) mean(simb$rt <= x), numeric(1))
  expect_lt(max(abs(emp - th)), 0.03)

  # Defective accumulator: omission rate 1 - h.
  pars_d <- cbind(K = 2, lambda_p = 5, lambda_m = 8, t0 = 0)[rep(1, n), ]
  simd <- EMC2:::ropr_cpp(pars_d, "r", rep(TRUE, n), 1L)
  expect_lt(abs(mean(!is.finite(simd$rt)) - (1 - (5 / 8)^2)), 0.03)
  fin <- is.finite(simd$rt)
  expect_lt(abs(mean(simd$rt[fin] <= .3) -
                  EMC2:::popr(.3, 2, 8, 5)), 0.03)
})

opr_fixture <- function(parameterization, n = 60) {
  set.seed(20261017)
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  dat <- dat[seq_len(n), ]
  model <- function() OPR(parameterization = parameterization)
  r <- switch(parameterization, rates = c("lambda_p", "lambda_m"),
              background = c("v", "lambda0"), c("v", "rho"))
  formula <- list(K ~ 1, stats::as.formula(paste(r[1], "~ lM")),
                  stats::as.formula(paste(r[2], "~ lM")), t0 ~ 1)
  des <- suppressMessages(design(data = dat, model = model, matchfun = matchfun,
                                 formula = formula))
  invisible(utils::capture.output(p <- sampled_pars(des)))
  p[] <- 0
  n0 <- names(p)
  if (parameterization == "rates") {
    p[c("K", "lambda_p", "lambda_p_lMTRUE", "lambda_m", "lambda_m_lMTRUE", "t0")] <-
      c(log(1.5), log(4), 0.9, log(9), -1.2, log(.15))
  } else if (parameterization == "background") {
    p[c("K", "v", "v_lMTRUE", "lambda0", "lambda0_lMTRUE", "t0")] <-
      c(log(1.4), log(4), 0.9, log(7), -0.5, log(.15))
  } else {
    p[c("K", "v", "v_lMTRUE", "rho", "rho_lMTRUE", "t0")] <-
      c(log(1.2), log(5), 0.8, log(.6), -0.4, log(.15))
  }
  stopifnot(identical(names(p), n0))
  list(dat = dat, des = des, p = p)
}

opr_ll <- function(emc, p) {
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
# [LT, UT] truncation normaliser that retains the never-finish atom.
opr_ref_ll <- function(dadm, pars, LT = 0, UT = Inf) {
  r <- EMC2:::.opr_rates(pars)
  K <- pars[, "K"]; t0 <- pars[, "t0"]
  S <- function(t, i) if (t <= t0[i]) 1 else
    if (is.infinite(t)) 1 - (if (r$a[i] >= r$b[i]) 1 else (r$a[i] / r$b[i])^K[i]) else
      1 - opr_ref_cdf(t - t0[i], K[i], r$a[i], r$b[i])
  n_lR <- nlevels(dadm$lR)
  trial <- rep(seq_len(nrow(dadm) / n_lR), each = n_lR)
  win <- as.logical(dadm$winner)
  ll <- vapply(split(seq_len(nrow(dadm)), trial), function(ix) {
    w <- ix[win[ix]]
    num <- log(opr_ref_pdf(dadm$rt[w] - t0[w], K[w], r$a[w], r$b[w])) +
      sum(log(vapply(ix[!win[ix]], function(i) S(dadm$rt[i], i), 0)))
    Z <- prod(vapply(ix, function(i) S(LT, i), 0)) -
      prod(vapply(ix, function(i) S(UT, i), 0)) +
      if (is.finite(UT)) prod(vapply(ix, function(i) S(Inf, i), 0)) else 0
    num - log(Z)
  }, numeric(1))
  sum(pmax(ll, log(1e-10)))
}

test_that("the compiled race likelihood matches the reference on both charts", {
  skip_on_cran()
  for (ch in c("ratio", "background", "rates")) {
    fx <- opr_fixture(ch)
    emc <- suppressMessages(make_emc(fx$dat, fx$des, type = "single", n_chains = 1,
                                     compress = FALSE, rt_resolution = NULL))
    dadm <- emc[[1]]$data[[1]]
    pars <- EMC2:::get_pars_matrix_oo(fx$p, dadm, emc[[1]]$model())
    expect_equal(opr_ll(emc, fx$p), opr_ref_ll(dadm, pars), tolerance = 1e-8,
                 info = ch)
  }
})

test_that("truncation keeps the defective atom in the normaliser", {
  skip_on_cran()
  fx <- opr_fixture("rates")
  LT <- 0.2; UT <- 1.2
  d2 <- fx$dat[fx$dat$rt >= LT & fx$dat$rt <= UT, ]
  d2$LT <- LT; d2$UT <- UT
  emc <- suppressMessages(make_emc(d2, fx$des, type = "single", n_chains = 1,
                                   compress = FALSE, rt_resolution = NULL))
  dadm <- emc[[1]]$data[[1]]
  pars <- EMC2:::get_pars_matrix_oo(fx$p, dadm, emc[[1]]$model())
  r <- EMC2:::.opr_rates(pars)
  expect_true(any(r$b > r$a))
  expect_equal(opr_ll(emc, fx$p), opr_ref_ll(dadm, pars, LT, UT),
               tolerance = 1e-8)
})

test_that("OPR exposes the race-model contract", {
  for (m in list(OPR(), OPR("background"), OPR("rates"))) {
    expect_identical(m$type, "RACE")
    expect_true("t0" %in% m$p_types_canonical)
    expect_identical(names(m$p_types)[1:4], m$p_types_canonical)
  }
  expect_identical(OPR()$c_name, "OPR")
  expect_identical(OPR("rates")$c_name, "OPR_RATES")
  expect_identical(OPR("background")$c_name, "OPR_BG")
})
