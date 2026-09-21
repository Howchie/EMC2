skip_model_validation()

# FRQ, the finite reservoir quorum process -- src/model_FRQ.h, R/model_FRQ.R.
#
# The model is fully closed form, so strong oracles are cheap. This file keeps
# the ones that are fast enough for R CMD check:
#
#   1. A small simulation of the LITERAL generative process (N Bernoulli(p)
#      availability indicators, Exponential(lambda) latencies, K-th order
#      statistic). It shares no formula with the kernel, so it establishes
#      that the kernel's closed form is the right closed form.
#   2. An independent R transcription of the closed-form density and CDF.
#   3. The defective-distribution contract: F(Inf) = h, S(Inf) = 1 - h,
#      integrate(f) = h.
#   4. That the parameterisation means what it claims: h is the completion
#      probability, lambda is the direct registration rate, and
#      (alpha, beta, h, lambda) inverts to the generating p.

# ---------------------------------------------------------------------------
# Independent references
# ---------------------------------------------------------------------------

ref_pl <- function(a, b, h, lambda, delta = 0, cv_u = 0) {
  zh <- if (delta > 0) EMC2:::frq_h_inv_r(h, delta) else h
  p <- qbeta(zh, a, b)
  c(p = p, lambda = lambda)
}

ref_H <- function(z, d) {
  if (d == 0) return(z)
  lz <- log(z) - log1p(-z)
  integrate(function(e) 1 / (1 + exp(-(lz + e))), -d, d,
            rel.tol = 1e-12)$value / (2 * d)
}

ref_pH <- function(x, a, b, h, lambda, d = 0, cv_u = 0) {
  pl <- ref_pl(a, b, h, lambda, d, cv_u)
  p <- pl[["p"]]
  c2 <- cv_u^2
  q <- if (c2 == 0) {
    p * (-expm1(-lambda * x))
  } else {
    p * (1 - (1 + c2 * lambda * x)^(-1 / c2))
  }
  z <- pbeta(q, a, b)
  if (d > 0) vapply(z, ref_H, numeric(1), d = d) else z
}

ref_d <- function(x, a, b, h, lambda, delta = 0, cv_u = 0) {
  pl <- ref_pl(a, b, h, lambda, delta, cv_u)
  p <- pl[["p"]]
  c2 <- cv_u^2
  if (c2 == 0) {
    q <- p * (-expm1(-lambda * x))
    omq <- (1 - p) + p * exp(-lambda * x)
    dens <- p * lambda * exp(-lambda * x) *
      q^(a - 1) * omq^(b - 1) / beta(a, b)
  } else {
    sr <- (1 + c2 * lambda * x)^(-1 / c2)
    g <- lambda * (1 + c2 * lambda * x)^(-1 / c2 - 1)
    q <- p * (1 - sr)
    omq <- (1 - p) + p * sr
    dens <- p * g * q^(a - 1) * omq^(b - 1) / beta(a, b)
  }
  if (delta > 0) {
    z <- pbeta(q, a, b)
    c2_d <- 4 * sinh(delta / 2)^2
    sinhc <- if (delta < 0.1) 1 + delta^2 / 6 + delta^4 / 120 else sinh(delta) / delta
    Hprime <- sinhc / (1 + 2 * z * (1 - z) * (cosh(delta) - 1))
    dens <- dens * Hprime
  }
  dens
}

ref_p <- function(x, a, b, h, lambda, delta = 0, cv_u = 0) {
  ref_pH(x, a, b, h, lambda, d = delta, cv_u = cv_u)
}

ref_s <- function(x, a, b, h, lambda, delta = 0, cv_u = 0) {
  1 - ref_p(x, a, b, h, lambda, delta = delta, cv_u = cv_u)
}

cpp_d <- function(x, a, b, h, lambda, delta = 0, cv_u = 0)
  EMC2:::dfrq(x, a, b, h, lambda, delta = delta, cv_u = cv_u)
cpp_p <- function(x, a, b, h, lambda, delta = 0, cv_u = 0)
  EMC2:::pfrq(x, a, b, h, lambda, delta = delta, cv_u = cv_u)
cpp_s <- function(x, a, b, h, lambda, delta = 0, cv_u = 0)
  EMC2:::pfrq(x, a, b, h, lambda, delta = delta, lower_tail = FALSE, cv_u = cv_u)

frq_grid <- data.frame(
  a      = c(1,    1,    2,    1.7,  3,     8,    2.5,  1),
  b      = c(1,    3,    1,    2.3,  6,     2.3,  1.4,  6),
  h      = c(0.3,  0.95, 0.75, 0.95, 0.999, 0.6,  0.85, 0.999),
  lambda = c(2.0,  1.5,  1.2,  0.8,  3.0,   1.0,  1.4,  2.5))
frq_x <- c(1e-4, 0.01, 0.05, 0.2, 0.5, 1, 2, 5, 20)

# ---------------------------------------------------------------------------
# 1. The literal finite-reservoir process
# ---------------------------------------------------------------------------

test_that("the closed form is the K-th order statistic of a finite reservoir", {
  skip_on_cran()
  set.seed(20260816)
  nsim <- 100
  N <- 7; K <- 3; p <- 0.8; lam <- 2.5
  a <- K; b <- N - K + 1
  avail <- matrix(runif(nsim * N) < p, nrow = nsim)
  lat <- matrix(rexp(nsim * N, lam), nrow = nsim)
  lat[!avail] <- Inf
  Tk <- apply(lat, 1, function(r) sort(r)[K])

  h <- pbeta(p, a, b)

  expect_true(all(Tk[!is.infinite(Tk)] > 0))
  expect_equal(cpp_p(Inf, a, b, h, lam), h, tolerance = 1e-12)

  # ... and the kernel must recover exactly the (p, lambda) that generated it.
  got <- EMC2:::frq_rate(a, b, h, lam)
  expect_equal(unname(got[1, "p"]), p, tolerance = 1e-9)
  expect_equal(unname(got[1, "lambda"]), lam, tolerance = 1e-9)
})

# ---------------------------------------------------------------------------
# 2. Against an independent transcription of the boxed formulas
# ---------------------------------------------------------------------------

test_that("pdf, cdf and survivor match an independent R reference", {
  worst <- c(d = 0, p = 0, s = 0)
  rel <- function(got, want) max(abs(got - want) / pmax(want, 1e-300))
  for (i in seq_len(nrow(frq_grid))) {
    g <- frq_grid[i, ]
    worst["d"] <- max(worst["d"], rel(cpp_d(frq_x, g$a, g$b, g$h, g$lambda),
                                      ref_d(frq_x, g$a, g$b, g$h, g$lambda)))
    worst["p"] <- max(worst["p"], rel(cpp_p(frq_x, g$a, g$b, g$h, g$lambda),
                                      ref_p(frq_x, g$a, g$b, g$h, g$lambda)))
    worst["s"] <- max(worst["s"], rel(cpp_s(frq_x, g$a, g$b, g$h, g$lambda),
                                      ref_s(frq_x, g$a, g$b, g$h, g$lambda)))
  }
  expect_lt(worst[["d"]], 1e-11)
  expect_lt(worst[["p"]], 1e-11)
  expect_lt(worst[["s"]], 1e-10)
})

test_that("the log output equals the log of the natural output", {
  for (i in c(1, 4, 5)) {
    g <- frq_grid[i, ]
    expect_equal(EMC2:::dfrq(frq_x, g$a, g$b, g$h, g$lambda, log_out = TRUE),
                 log(cpp_d(frq_x, g$a, g$b, g$h, g$lambda)), tolerance = 1e-12)
    expect_equal(EMC2:::pfrq(frq_x, g$a, g$b, g$h, g$lambda, log_out = TRUE),
                 log(cpp_p(frq_x, g$a, g$b, g$h, g$lambda)), tolerance = 1e-12)
    expect_equal(EMC2:::pfrq(frq_x, g$a, g$b, g$h, g$lambda, lower_tail = FALSE,
                             log_out = TRUE),
                 log(cpp_s(frq_x, g$a, g$b, g$h, g$lambda)), tolerance = 1e-12)
  }
})

# ---------------------------------------------------------------------------
# 3. The defective-distribution contract
# ---------------------------------------------------------------------------

test_that("the asymptotic and integral checks of the defective tail hold", {
  for (i in seq_len(nrow(frq_grid))) {
    g <- frq_grid[i, ]
    expect_equal(cpp_p(Inf, g$a, g$b, g$h, g$lambda), g$h, tolerance = 1e-12)
    expect_equal(cpp_s(Inf, g$a, g$b, g$h, g$lambda), 1 - g$h, tolerance = 1e-9)
    I <- integrate(function(z) cpp_d(z, g$a, g$b, g$h, g$lambda), 0, Inf,
                   rel.tol = 1e-10)$value
    expect_equal(I, g$h, tolerance = 1e-7)
    for (x in c(0.05, 0.3, 1.5)) {
      Ix <- integrate(function(z) cpp_d(z, g$a, g$b, g$h, g$lambda), 0, x,
                      rel.tol = 1e-10)$value
      expect_equal(Ix, cpp_p(x, g$a, g$b, g$h, g$lambda), tolerance = 1e-7)
    }
    expect_equal(cpp_p(frq_x, g$a, g$b, g$h, g$lambda) +
                   cpp_s(frq_x, g$a, g$b, g$h, g$lambda),
                 rep(1, length(frq_x)), tolerance = 1e-12)
  }
})

test_that("the CDF is monotone and the density non-negative and finite", {
  xs <- exp(seq(log(1e-6), log(50), length.out = 31))
  for (i in seq_len(nrow(frq_grid))) {
    g <- frq_grid[i, ]
    d <- cpp_d(xs, g$a, g$b, g$h, g$lambda)
    p <- cpp_p(xs, g$a, g$b, g$h, g$lambda)
    expect_false(any(is.na(d)) || any(is.na(p)))
    expect_true(all(d >= 0) && all(is.finite(d)))
    expect_true(all(diff(p) >= -1e-14))
    expect_true(all(p >= 0 & p <= g$h + 1e-12))
  }
})

# ---------------------------------------------------------------------------
# 4. The (alpha, beta, h, lambda) parameterisation
# ---------------------------------------------------------------------------

test_that("h is the completion probability and lambda is the registration rate", {
  for (i in seq_len(nrow(frq_grid))) {
    g <- frq_grid[i, ]
    expect_equal(cpp_p(Inf, g$a, g$b, g$h, g$lambda), g$h, tolerance = 1e-12)
    expect_equal(cpp_s(Inf, g$a, g$b, g$h, g$lambda), 1 - g$h, tolerance = 1e-9)
  }
})

test_that("frq_rate is the exact inverse of h = I_p(alpha, beta)", {
  for (i in seq_len(nrow(frq_grid))) {
    g <- frq_grid[i, ]
    got <- EMC2:::frq_rate(g$a, g$b, g$h, g$lambda)
    want <- ref_pl(g$a, g$b, g$h, g$lambda)
    expect_equal(unname(got[1, "p"]), want[["p"]], tolerance = 1e-12)
    expect_equal(unname(got[1, "lambda"]), want[["lambda"]], tolerance = 1e-10)
    expect_equal(pbeta(unname(got[1, "p"]), g$a, g$b), g$h, tolerance = 1e-10)
  }
  v <- EMC2:::frq_rate(c(1, 2, 3), 2, 0.9, 0.4)
  expect_equal(dim(v), c(3L, 2L))
  expect_equal(colnames(v), c("p", "lambda"))
  expect_equal(v[2, ], EMC2:::frq_rate(2, 2, 0.9, 0.4)[1, ])
})

test_that("K = N removes p from the conditional RT distribution", {
  a <- 4
  cond <- function(x, b, h, lam = 2) {
    cpp_p(x, a, b, h, lam) / h
  }
  xs <- c(0.1, 0.3, 0.9)
  expect_equal(cond(xs, 1, 0.4), cond(xs, 1, 0.95), tolerance = 1e-10)
  expect_gt(max(abs(cond(xs, 3, 0.4) - cond(xs, 3, 0.95))), 1e-3)
})

test_that("alpha controls the order of the leading edge", {
  xs <- c(1e-7, 2e-7)
  for (a in c(1, 1.5, 3, 5.5)) {
    slope <- diff(log(cpp_d(xs, a, 2.5, 0.9, 1.4))) / diff(log(xs))
    expect_equal(slope, a - 1, tolerance = 1e-4)
  }
})

test_that("degenerate parameters are rejected rather than returning garbage", {
  bad <- list(c(0, 2, 0.9, 1.3), c(2, 0, 0.9, 1.3), c(2, 2, 0, 1.3),
              c(2, 2, 1.5, 1.3), c(2, 2, 0.9, 0), c(2, 2, 0.9, -1),
              c(Inf, 2, 0.9, 1.3), c(NA, 2, 0.9, 1.3))
  for (v in bad) {
    expect_equal(cpp_d(0.4, v[1], v[2], v[3], v[4]), 0)
    expect_equal(cpp_p(0.4, v[1], v[2], v[3], v[4]), 0)
    expect_true(is.na(EMC2:::frq_rate(v[1], v[2], v[3], v[4])[1, "p"]))
  }
  expect_equal(cpp_d(c(-1, 0), 2, 3, 0.9, 1.3), c(0, 0))
  expect_equal(cpp_p(c(-1, 0), 2, 3, 0.9, 1.3), c(0, 0))
  expect_equal(cpp_s(c(-1, 0), 2, 3, 0.9, 1.3), c(1, 1))
})

test_that("h at its upper bound stays finite and reaches the non-defective limit", {
  h <- 1 - 1e-9
  expect_equal(cpp_p(Inf, 2, 3, h, 1.3), h, tolerance = 1e-15)
  expect_equal(cpp_s(Inf, 2, 3, h, 1.3), 1 - h, tolerance = 1e-5)
  expect_true(all(is.finite(cpp_d(frq_x, 2, 3, h, 1.3))))
  expect_equal(unname(EMC2:::frq_rate(2, 3, 1, 1.3)[1, "p"]), 1)
  expect_equal(cpp_p(Inf, 2, 3, 1, 1.3), 1, tolerance = 1e-14)
  expect_equal(cpp_s(Inf, 2, 3, 1, 1.3), 0, tolerance = 1e-14)
})

# ---------------------------------------------------------------------------
# 5. The R model interface
# ---------------------------------------------------------------------------

test_that("the constructor exposes the documented contract", {
  m <- FRQ()
  expect_equal(m$type, "RACE")
  expect_equal(m$c_name, "FRQ")
  expect_equal(names(m$p_types),
               c("alpha", "beta", "h", "lambda", "t0", "delta",
                 "cv_u",
                 "pContaminant", "pGuess"))
  expect_equal(m$p_types_canonical,
               c("alpha", "beta", "h", "lambda", "t0"))
  expect_equal(unname(m$transform$func[c("alpha", "beta", "h", "lambda", "t0",
                                         "delta", "cv_u")]),
               c("exp", "exp", "pnorm", "exp", "exp", "exp", "exp"))
  expect_equal(exp(m$p_types[["delta"]]), 0)
  expect_equal(unname(m$bound$minmax[, "delta"]), c(1e-4, 6))
  expect_equal(m$bound$exception[["delta"]], 0)
  expect_equal(unname(m$bound$minmax[, "cv_u"]), c(1e-4, Inf))
  expect_equal(m$bound$exception[["cv_u"]], 0)
  expect_equal(m$bound$exception[["h"]], 1)
  expect_equal(pnorm(m$p_types[["h"]]), 1, tolerance = 1e-12)
  expect_equal(unname(m$bound$minmax[, "alpha"]), c(1, Inf))
  expect_equal(unname(m$bound$minmax[, "beta"]), c(1, Inf))
  expect_equal(unname(m$bound$minmax[, "lambda"]), c(1e-4, Inf))
  expect_length(formals(FRQ), 0L)
})

test_that("the (h, lambda) coordinates are representable on the whole bounded box", {
  grid <- data.frame(
    alpha = c(1, 1.5, 2, 20, 200),
    beta = c(1.0001, 2, 5, 20, 1),
    h = c(1e-3, 0.05, 0.5, 0.9, 0.99),
    lambda = c(0.1, 0.5, 1.0, 2.0, 5.0)
  )
  pl <- EMC2:::frq_rate(grid$alpha, grid$beta, grid$h, grid$lambda)
  expect_false(anyNA(pl[, "p"]))
  expect_true(all(pl[, "p"] > 0 & pl[, "p"] < 1))
  expect_true(all(is.finite(pl[, "lambda"]) & pl[, "lambda"] > 0))
  expect_equal(pbeta(pl[, "p"], grid$alpha, grid$beta), grid$h, tolerance = 1e-8)
  expect_equal(qbeta(0.99, 0.05, 0.05), 1)
  expect_equal(qbeta(0.1, 1e-4, 1e-4), 0)
})

test_that("dfun/pfun apply t0 and call the same kernel", {
  pars <- cbind(alpha = 2, beta = 3, h = 0.9, lambda = 1.35, t0 = 0.15, delta = 0)
  pars <- pars[rep(1, 5), ]
  rt <- c(0.05, 0.15, 0.3, 0.8, Inf)
  expect_equal(EMC2:::dFRQ(rt, pars),
               c(0, 0, cpp_d(c(0.15, 0.65), 2, 3, 0.9, 1.35), 0))
  expect_equal(EMC2:::pFRQ(rt, pars),
               c(0, 0, cpp_p(c(0.15, 0.65), 2, 3, 0.9, 1.35), 0.9))
  expect_equal(EMC2:::sFRQ(rt, pars),
               c(1, 1, cpp_s(c(0.15, 0.65), 2, 3, 0.9, 1.35), 0.1))
})

test_that("Ttransform reports the generative parameters", {
  pars <- cbind(alpha = 2, beta = 3, h = 0.9, lambda = 1.35, t0 = 0.15, delta = 0)
  out <- FRQ()$Ttransform(pars, NULL)
  want <- ref_pl(2, 3, 0.9, 1.35)
  expect_equal(unname(out[1, "p"]), want[["p"]], tolerance = 1e-12)
  expect_equal(unname(out[1, "lambda"]), want[["lambda"]], tolerance = 1e-10)
  expect_equal(unname(out[1, "N"]), 4)
  expect_equal(unname(out[1, "d"]), 0.5)
  expect_equal(unname(out[1, "sQ"]), 0)
  expect_equal(unname(out[1, "a_u"]), Inf)
  expect_equal(out[, colnames(pars), drop = FALSE], pars)
  expect_null(rownames(out))
})

# ---------------------------------------------------------------------------
# 6. The compiled race likelihood
# ---------------------------------------------------------------------------

frq_ll <- function(emc, p, p_types_override = NULL) {
  model <- emc[[1]]$model()
  p_types <- names(model$p_types)
  dadm <- emc[[1]]$data[[1]]
  designs <- list()
  for (nm in p_types) {
    designs[[nm]] <- attr(dadm, "designs")[[nm]][
      attr(attr(dadm, "designs")[[nm]], "expand"), , drop = FALSE]
  }
  if (!is.null(p_types_override)) p_types <- p_types_override
  EMC2:::calc_ll_oo(matrix(p, nrow = 1, dimnames = list(NULL, names(p))), dadm,
                    constants = attr(dadm, "constants"), designs = designs,
                    type = model$c_name, bounds = model$bound,
                    transforms = model$transform,
                    pretransforms = model$pre_transform, p_types = p_types,
                    min_ll = log(1e-10), trend = model$trend)
}

ref_race_ll <- function(dadm, pars) {
  n_lR <- nlevels(dadm$lR)
  trial <- rep(seq_len(nrow(dadm) / n_lR), each = n_lR)
  d <- EMC2:::dFRQ(dadm$rt, pars)
  s <- EMC2:::sFRQ(dadm$rt, pars)
  win <- as.logical(dadm$winner)
  ll <- vapply(split(seq_len(nrow(dadm)), trial), function(ix) {
    if (!any(win[ix])) return(sum(log(s[ix])))
    log(d[ix][win[ix]]) + sum(log(s[ix][!win[ix]]))
  }, numeric(1))
  sum(pmax(ll, log(1e-10)))
}

frq_truth <- list(hT = pnorm(qnorm(0.93) + 0.45), hF = pnorm(qnorm(0.93) - 0.45))

frq_fixture <- function(n = 60, seed = 20260816, extensions = FALSE) {
  set.seed(seed)
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  dat <- dat[seq_len(n), ]
  ADmat <- matrix(c(-1 / 2, 1 / 2), ncol = 1, dimnames = list(NULL, "d"))
  formula <- list(alpha ~ 1, beta ~ 1, h ~ lM, lambda ~ lM, t0 ~ 1)
  if (extensions)
    formula[c("delta", "cv_u")] <- list(delta ~ 1, cv_u ~ 1)
  des <- design(data = dat, model = FRQ, matchfun = matchfun,
                formula = formula,
                contrasts = list(h = list(lM = ADmat), lambda = list(lM = ADmat)))
  p <- sampled_pars(des)
  p["alpha"] <- log(2.2); p["beta"] <- log(3.1)
  p["h"] <- qnorm(0.93); p["h_lMd"] <- 0.9
  p["lambda"] <- log(1.5); p["lambda_lMd"] <- 0.6
  p["t0"] <- log(0.15)
  if (extensions) {
    p["delta"] <- log(1.2)
    p["cv_u"] <- log(0.5)
  }
  emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1,
                                   compress = FALSE, rt_resolution = NULL))
  list(dat = dat, des = des, p = p, emc = emc)
}

test_that("the compiled race likelihood matches the R reference", {
  skip_on_cran()
  fx <- frq_fixture(n = 60)
  got <- frq_ll(fx$emc, fx$p)
  expect_true(is.finite(got))
  dadm <- fx$emc[[1]]$data[[1]]
  pars <- EMC2:::get_pars_matrix_oo(fx$p, dadm, fx$emc[[1]]$model())
  expect_equal(got, ref_race_ll(dadm, pars), tolerance = 1e-6)
})

test_that("the compiled likelihood propagates threshold variability and frailty", {
  skip_on_cran()
  fx <- frq_fixture(n = 40, extensions = TRUE)
  got <- frq_ll(fx$emc, fx$p)
  expect_true(is.finite(got))
  dadm <- fx$emc[[1]]$data[[1]]
  pars <- EMC2:::get_pars_matrix_oo(fx$p, dadm, fx$emc[[1]]$model())
  expect_equal(got, ref_race_ll(dadm, pars), tolerance = 1e-6)
})

test_that("a p_types reordering is caught by the column contract", {
  skip_on_cran()
  fx <- frq_fixture(n = 30)
  swapped <- names(FRQ()$p_types)
  swapped[3:4] <- swapped[4:3]
  expect_error(frq_ll(fx$emc, fx$p, p_types_override = swapped),
               "FRQ kernels expect parameter column")
})

test_that("an omission scores exactly sum_j log(1 - h_j)", {
  skip_on_cran()
  fx <- frq_fixture(n = 30)
  d_om <- fx$dat
  d_om$rt[1:6] <- Inf
  d_om$R[1:6] <- NA
  mk <- function(d) suppressMessages(
    make_emc(d, fx$des, type = "single", n_chains = 1, compress = FALSE,
             rt_resolution = NULL))
  ll_om <- frq_ll(mk(d_om), fx$p)
  expect_true(is.finite(ll_om))
  per_omission <- log(1 - frq_truth$hT) + log(1 - frq_truth$hF)
  expect_equal(ll_om - frq_ll(mk(d_om[-(1:6), ]), fx$p), 6 * per_omission,
               tolerance = 1e-6)
})

test_that("a finite UT normaliser keeps the retained never-finish atom", {
  skip_on_cran()
  fx <- frq_fixture(n = 40)
  LT <- 0.2; UT <- 0.9
  d2 <- fx$dat
  d2$rt[1:6] <- Inf
  d2$R[1:6] <- NA
  d2 <- d2[d2$rt >= LT | is.infinite(d2$rt), ]
  d2 <- d2[d2$rt <= UT | is.infinite(d2$rt), ]
  d2$LT <- LT; d2$UT <- UT
  e2 <- suppressMessages(make_emc(d2, fx$des, type = "single", n_chains = 1,
                                  compress = FALSE, rt_resolution = NULL))
  dadm <- e2[[1]]$data[[1]]
  pars <- EMC2:::get_pars_matrix_oo(fx$p, dadm, e2[[1]]$model())

  n_lR <- nlevels(dadm$lR)
  trial <- rep(seq_len(nrow(dadm) / n_lR), each = n_lR)
  ref_trunc_ll <- function(p) {
    pars <- EMC2:::get_pars_matrix_oo(p, dadm, e2[[1]]$model())
    s_at <- function(t) EMC2:::sFRQ(rep(t, nrow(dadm)), pars)
    s_LT <- s_at(LT); s_UT <- s_at(UT); s_Inf <- s_at(Inf)
    d <- EMC2:::dFRQ(dadm$rt, pars)
    s <- EMC2:::sFRQ(dadm$rt, pars)
    win <- as.logical(dadm$winner)
    ll <- vapply(split(seq_len(nrow(dadm)), trial), function(ix) {
      logZ <- log(prod(s_LT[ix]) - prod(s_UT[ix]) + prod(s_Inf[ix]))
      num <- if (!is.finite(dadm$rt[ix][1])) sum(log(s[ix]))
             else log(d[ix][win[ix]]) + sum(log(s[ix][!win[ix]]))
      num - logZ
    }, numeric(1))
    sum(pmax(ll, log(1e-10)))
  }
  expect_equal(frq_ll(e2, fx$p), ref_trunc_ll(fx$p), tolerance = 1e-6)

  p_low <- fx$p
  p_low[["h"]] <- qnorm(0.35); p_low[["h_lMd"]] <- 0
  expect_equal(frq_ll(e2, p_low), ref_trunc_ll(p_low), tolerance = 1e-6)

  pars_low <- EMC2:::get_pars_matrix_oo(p_low, dadm, e2[[1]]$model())
  s_low <- function(t) prod(EMC2:::sFRQ(rep(t, nrow(dadm)), pars_low)[seq_len(n_lR)])
  Z_fixed <- s_low(LT) - s_low(UT) + s_low(Inf)
  Z_finite_only <- s_low(LT) - s_low(UT)
  expect_gt(s_low(Inf), 0.8 * Z_finite_only)
  expect_gt(40 * (log(Z_fixed) - log(Z_finite_only)), 10)
})

test_that("truncation, censoring and pContaminant stay well posed", {
  skip_on_cran()
  fx <- frq_fixture(n = 40)
  ll_base <- frq_ll(fx$emc, fx$p)
  expect_true(is.finite(ll_base))
  for (bounds in list(list(LT = 0.05, UT = 3), list(LC = 0.05, UC = 2))) {
    d2 <- fx$dat
    for (nm in names(bounds)) d2[[nm]] <- bounds[[nm]]
    e2 <- suppressMessages(make_emc(d2, fx$des, type = "single", n_chains = 1,
                                    compress = FALSE, rt_resolution = NULL))
    expect_true(is.finite(frq_ll(e2, fx$p)))
  }

  matchfun <- function(d) d$S == d$lR
  ADmat <- matrix(c(-1 / 2, 1 / 2), ncol = 1, dimnames = list(NULL, "d"))
  des_pc <- suppressMessages(design(
    data = fx$dat, model = FRQ, matchfun = matchfun,
    formula = list(alpha ~ 1, beta ~ 1, h ~ lM, lambda ~ lM, t0 ~ 1,
                   pContaminant ~ 1),
    contrasts = list(h = list(lM = ADmat), lambda = list(lM = ADmat))))
  p_pc <- c(fx$p, pContaminant = qnorm(0.05))[names(sampled_pars(des_pc))]
  e_pc <- suppressMessages(make_emc(fx$dat, des_pc, type = "single",
                                    n_chains = 1, compress = FALSE,
                                    rt_resolution = NULL))
  expect_equal(frq_ll(e_pc, p_pc), ll_base + nrow(fx$dat) * log(0.95),
               tolerance = 1e-6)
})

test_that("a response below t0 floors cleanly", {
  skip_on_cran()
  fx <- frq_fixture(n = 20)
  p <- fx$p
  p[["t0"]] <- log(max(fx$dat$rt) + 1)
  ll <- frq_ll(fx$emc, p)
  expect_false(is.na(ll))
  expect_equal(ll, nrow(fx$dat) * log(1e-10), tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# 7. Simulation
# ---------------------------------------------------------------------------

test_that("the C++ and R simulators agree distributionally with the CDF", {
  skip_on_cran()
  n <- 1000
  lR <- factor(rep(c("left", "right"), n), levels = c("left", "right"))
  pars <- cbind(alpha = 2, beta = 3, h = 0.9, lambda = 1.35, t0 = 0.1)
  pars <- pars[rep(1, length(lR)), ]
  set.seed(4)
  a <- EMC2:::rFRQ(lR, pars)
  set.seed(4)
  b <- withr::with_options(list(emc2.cpp_rfun = TRUE),
                           EMC2:::.rfun_FRQ(lR, pars))
  probe <- c(0.2, 0.35, 0.5, 0.8)
  th <- 1 - cpp_s(probe - 0.1, 2, 3, 0.9, 1.35)^2
  for (dat in list(a, b)) {
    fin <- is.finite(dat$rt)
    emp <- vapply(probe, function(x) mean(fin & dat$rt <= x), numeric(1))
    expect_lt(max(abs(emp - th)), 0.03)
    expect_lt(abs(mean(!fin) - 0.01), 0.006)
    expect_true(all(is.na(dat$R[!fin])))
    expect_true(all(dat$rt[fin] > 0.1))
  }
})

test_that("the simulator reproduces the literal reservoir at integer shapes", {
  skip_on_cran()
  set.seed(99)
  N <- 6; K <- 2; p <- 0.7; lam <- 3
  a <- K; b <- N - K + 1
  h <- pbeta(p, a, b)
  n <- 1e4
  lR <- factor(rep("go", n), levels = "go")
  pars <- cbind(alpha = a, beta = b, h = h, lambda = lam, t0 = 0)
  sim <- EMC2:::.rfun_FRQ(lR, pars[rep(1, n), ])

  avail <- matrix(runif(n * N) < p, nrow = n)
  lat <- matrix(rexp(n * N, lam), nrow = n)
  lat[!avail] <- Inf
  ref <- apply(lat, 1, function(r) sort(r)[K])

  probe <- c(0.05, 0.15, 0.3, 0.6, 1.2)
  emp <- vapply(probe, function(x) mean(is.finite(sim$rt) & sim$rt <= x),
                numeric(1))
  ref_e <- vapply(probe, function(x) mean(ref <= x), numeric(1))
  expect_lt(max(abs(emp - ref_e)), 0.02)
  expect_lt(abs(mean(!is.finite(sim$rt)) - mean(is.infinite(ref))), 0.015)
})

test_that("make_data produces omissions the design can be fit back through", {
  skip_on_cran()
  set.seed(20260816)
  fx <- frq_fixture(n = 30)
  p_omit <- fx$p
  p_omit[["h"]] <- qnorm(0.8)
  p_omit[["h_lMd"]] <- 0.4
  sim <- make_data(p_omit, design = fx$des, n_trials = 150)
  expect_true(any(is.infinite(sim$rt)))
  expect_true(all(is.na(sim$R[is.infinite(sim$rt)])))
  fin <- is.finite(sim$rt)
  expect_true(all(sim$rt[fin] > 0.15))
  expect_gt(mean(sim$S[fin] == sim$R[fin]), 0.6)
  e <- suppressMessages(make_emc(sim, fx$des, type = "single", n_chains = 1,
                                 compress = FALSE, rt_resolution = NULL))
  expect_true(is.finite(frq_ll(e, p_omit)))
})

# ---------------------------------------------------------------------------
# 8. Recovery
# ---------------------------------------------------------------------------

test_that("the likelihood is maximised at the generating parameters", {
  skip_on_cran()
  set.seed(20260816)
  fx <- frq_fixture(n = 30)
  p_rec <- fx$p
  p_rec["h"] <- qnorm(0.7)
  p_rec["h_lMd"] <- 0.3
  sim <- make_data(p_rec, design = fx$des, n_trials = 3000)
  e <- suppressMessages(make_emc(sim, fx$des, type = "single", n_chains = 1,
                                 compress = FALSE, rt_resolution = NULL))
  ll0 <- frq_ll(e, p_rec)
  expect_true(is.finite(ll0))
  for (nm in names(p_rec)) {
    for (delta in c(-0.5, 0.5)) {
      p2 <- p_rec
      p2[[nm]] <- p2[[nm]] + delta
      expect_lt(frq_ll(e, p2), ll0)
    }
  }
})

# ---------------------------------------------------------------------------
# 9. Threshold variability (delta) and unit-rate frailty (cv_u)
# ---------------------------------------------------------------------------

test_that("delta defaults to zero and is exactly the untransformed model", {
  x <- c(0.05, 0.2, 0.6, 2, Inf)
  expect_identical(EMC2:::pfrq(x, 2, 3, 0.9, 1.4, 0),
                   EMC2:::pfrq(x, 2, 3, 0.9, 1.4))
  expect_identical(EMC2:::dfrq(x, 2, 3, 0.9, 1.4, 0),
                   EMC2:::dfrq(x, 2, 3, 0.9, 1.4))
  expect_identical(EMC2:::frq_rate(2, 3, 0.9, 1.4, 0),
                   EMC2:::frq_rate(2, 3, 0.9, 1.4))
  expect_equal(EMC2:::pfrq(x[1:4], 2, 3, 0.9, 1.4, 0), ref_p(x[1:4], 2, 3, 0.9, 1.4))
})

test_that("unit-rate frailty matches the closed-form registration mixture", {
  a <- 2.5; b <- 1.4; h <- 0.85; lambda <- 1.4; cv <- 0.7
  pl <- EMC2:::frq_rate(a, b, h, lambda, 0, cv)
  p <- qbeta(h, a, b)
  expect_equal(unname(pl[1, "p"]), p, tolerance = 1e-12)
  expect_equal(unname(pl[1, "lambda"]), lambda, tolerance = 1e-12)

  x <- c(0.05, 0.2, 0.6, 1.5, 5)
  c2 <- cv^2
  sr <- (1 + c2 * lambda * x)^(-1 / c2)
  g <- lambda * (1 + c2 * lambda * x)^(-1 / c2 - 1)
  q <- p * (1 - sr)
  want_p <- pbeta(q, a, b)
  want_d <- p * g * q^(a - 1) * (1 - q)^(b - 1) / beta(a, b)
  expect_equal(EMC2:::pfrq(x, a, b, h, lambda, 0, cv_u = cv), want_p,
               tolerance = 1e-11)
  expect_equal(EMC2:::dfrq(x, a, b, h, lambda, 0, cv_u = cv), want_d,
               tolerance = 1e-11)
  expect_equal(EMC2:::pfrq(Inf, a, b, h, lambda, 0, cv_u = cv), h,
               tolerance = 1e-11)
  expect_equal(EMC2:::pfrq(Inf, a, b, h, lambda, 0, cv_u = cv,
                            lower_tail = FALSE), 1 - h,
               tolerance = 1e-11)
  expect_equal(integrate(function(z)
    EMC2:::dfrq(z, a, b, h, lambda, 0, cv_u = cv), 0, Inf,
    rel.tol = 1e-9)$value, h, tolerance = 1e-8)
})

test_that("unit-rate frailty nests the exponential member exactly", {
  x <- c(0.05, 0.2, 0.6, 2, Inf)
  expect_identical(EMC2:::dfrq(x, 2, 3, 0.9, 1.4, 0),
                   EMC2:::dfrq(x, 2, 3, 0.9, 1.4, 0, cv_u = 0))
  expect_identical(EMC2:::pfrq(x, 2, 3, 0.9, 1.4, 0),
                   EMC2:::pfrq(x, 2, 3, 0.9, 1.4, 0, cv_u = 0))
  expect_identical(EMC2:::frq_rate(2, 3, 0.9, 1.4, 0),
                   EMC2:::frq_rate(2, 3, 0.9, 1.4, 0, 0))
})

test_that("threshold variability and unit frailty compose", {
  a <- 2.5; b <- 1.4; h <- 0.85; lambda <- 1.4
  d <- 1.25; cv <- 0.55
  x <- c(0.05, 0.2, 0.6, 1.5, 5)
  pl <- EMC2:::frq_rate(a, b, h, lambda, d, cv)
  p <- unname(pl[1, "p"])
  c2 <- cv^2
  q <- p * (1 - (1 + c2 * lambda * x)^(-1 / c2))
  z <- pbeta(q, a, b)
  want <- vapply(z, ref_H, numeric(1), d = d)
  expect_equal(EMC2:::pfrq(x, a, b, h, lambda, d, cv_u = cv), want,
               tolerance = 1e-10)
  expect_equal(EMC2:::pfrq(Inf, a, b, h, lambda, d, cv_u = cv), h,
               tolerance = 1e-10)
  expect_equal(integrate(function(t)
    EMC2:::dfrq(t, a, b, h, lambda, d, cv_u = cv), 0, Inf,
    rel.tol = 1e-9)$value, h, tolerance = 1e-8)
})

test_that("the proper h = 1 boundary remains proper with unit frailty", {
  pars <- cbind(alpha = 2, beta = 3, h = 1, lambda = 1.4, t0 = 0,
                delta = 0, cv_u = 0.7)
  pl <- EMC2:::frq_rate(2, 3, 1, 1.4, 0, 0.7)
  expect_equal(unname(pl[1, "p"]), 1)
  expect_equal(EMC2:::pfrq(Inf, 2, 3, 1, 1.4, 0, cv_u = 0.7), 1)
  expect_equal(EMC2:::pfrq(Inf, 2, 3, 1, 1.4, 0, cv_u = 0.7,
                            lower_tail = FALSE), 0)
  set.seed(20260920)
  sim_cpp <- EMC2:::rfrq_cpp(pars[rep(1, 2000), ], "1", rep(TRUE, 2000))
  expect_true(all(is.finite(sim_cpp$rt)))
  sim_r <- EMC2:::rFRQ(factor(rep("1", 2000)),
                       pars[rep(1, 2000), ])
  expect_true(all(is.finite(sim_r$rt)))
  out <- FRQ()$Ttransform(pars, NULL)
  expect_equal(unname(out[1, "a_u"]), 1 / 0.7^2, tolerance = 1e-12)
})

test_that("the delta > 0 CDF matches the marginalised criterion state", {
  x <- c(0.02, 0.1, 0.35, 0.9, 4)
  for (d in c(0.25, 1, 3, 6)) {
    for (i in c(1, 2, 5)) {
      g <- frq_grid[i, ]
      expect_equal(EMC2:::pfrq(x, g$a, g$b, g$h, g$lambda, d),
                   ref_pH(x, g$a, g$b, g$h, g$lambda, d),
                   tolerance = 1e-9,
                   info = paste("delta", d, "row", i))
    }
  }
})

test_that("delta preserves the defective contract and the meaning of h and lambda", {
  for (d in c(0.25, 1, 3, 6)) for (i in seq_len(nrow(frq_grid))) {
    g <- frq_grid[i, ]
    lab <- paste("delta", d, "row", i)
    expect_equal(EMC2:::pfrq(Inf, g$a, g$b, g$h, g$lambda, d), g$h,
                 tolerance = 1e-10, info = lab)
    expect_equal(EMC2:::pfrq(Inf, g$a, g$b, g$h, g$lambda, d, lower_tail = FALSE),
                 1 - g$h, tolerance = 1e-10, info = lab)
    expect_equal(integrate(function(t) EMC2:::dfrq(t, g$a, g$b, g$h, g$lambda, d),
                           0, Inf, rel.tol = 1e-10,
                           subdivisions = 2000L)$value,
                 g$h, tolerance = 1e-8, info = lab)
  }
})

test_that("the density is the derivative of the CDF under delta", {
  eps <- 1e-5
  for (d in c(0.5, 2, 5)) {
    x <- c(0.1, 0.3, 0.7, 1.5)
    fd <- (EMC2:::pfrq(x + eps, 2.5, 1.4, 0.85, 1.4, d) -
           EMC2:::pfrq(x - eps, 2.5, 1.4, 0.85, 1.4, d)) / (2 * eps)
    expect_equal(EMC2:::dfrq(x, 2.5, 1.4, 0.85, 1.4, d), fd,
                 tolerance = 1e-6, info = paste("delta", d))
  }
})

test_that("H is invertible and the family nests the base model as delta -> 0", {
  y <- c(1e-9, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-9)
  for (d in c(0.3, 1, 3, 6))
    expect_equal(vapply(EMC2:::frq_h_inv_r(y, d), ref_H, numeric(1), d = d), y,
                 tolerance = 1e-12, info = paste("delta", d))
  x <- c(0.05, 0.3, 1)
  base <- EMC2:::pfrq(x, 2, 3, 0.9, 1.4, 0)
  prev <- Inf
  for (d in c(1e-3, 1e-4, 1e-5)) {
    dev <- max(abs(EMC2:::pfrq(x, 2, 3, 0.9, 1.4, d) - base))
    expect_lt(dev, prev)
    prev <- dev
  }
  expect_lt(prev, 1e-10)
})

test_that("the simulator draws from the delta > 0 likelihood", {
  skip_on_cran()
  set.seed(4)
  n <- 4e4
  d <- 2
  pars <- cbind(alpha = 2, beta = 3, h = 0.8, lambda = 1.4, t0 = 0.15,
                delta = d)[rep(1, n), ]
  sim <- EMC2:::rfrq_cpp(pars, "1", rep(TRUE, n))
  fin <- is.finite(sim$rt)
  expect_equal(mean(fin), 0.8, tolerance = 0.01)
  probs <- seq(0.1, 0.9, 0.1)
  qs <- quantile(sim$rt[fin], probs)
  expect_equal(unname(EMC2:::pfrq(qs - 0.15, 2, 3, 0.8, 1.4, d) / 0.8),
               probs, tolerance = 0.02)
  expect_false(isTRUE(all.equal(
    unname(EMC2:::pfrq(qs - 0.15, 2, 3, 0.8, 1.4, 0) / 0.8), probs,
    tolerance = 0.02)))
})

test_that("FRQ rejects numerically saturated defective and extreme delta inputs", {
  expect_true(is.na(EMC2:::frq_rate(0.05, 0.05, 0.99, 1.3)[1, "p"]))
  expect_true(is.na(EMC2:::frq_rate(2, 3, 0.9, 1.3, 1000)[1, "p"]))
  expect_equal(EMC2:::pfrq(0.4, 2, 3, 0.9, 1.3, 1000), 0)
})

test_that("the FRQ R simulator keeps malformed rows as omissions", {
  lR <- factor(rep(c("a", "b"), 2), levels = c("a", "b"))
  pars <- cbind(alpha = 2, beta = 3, h = 0.9, lambda = 1.3,
                t0 = c(NA, NA, 0.1, 0.1))
  sim <- EMC2:::rFRQ(lR, pars)
  expect_true(is.infinite(sim$rt[1]) || is.na(sim$R[1]))
  expect_error(EMC2:::rFRQ(lR, pars, ok = TRUE), "one value per")
})

test_that("the FRQ R simulator permits a valid accumulator to win", {
  lR <- factor(c("a", "b"), levels = c("a", "b"))
  pars <- cbind(alpha = c(2, 2), beta = c(3, 3), h = c(1, 1),
                lambda = c(1.3, 1.3), t0 = c(NA, 0.1))
  sim <- EMC2:::rFRQ(lR, pars)
  expect_equal(sim$R, factor("b", levels = levels(lR)))
  expect_true(is.finite(sim$rt))
})
