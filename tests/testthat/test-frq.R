# FRQ, the finite reservoir quorum process -- src/model_FRQ.h, R/model_FRQ.R.
#
# The model is fully closed form, so strong oracles are cheap.  This file keeps
# the ones that are fast enough for R CMD check:
#
#   1. A small simulation of the LITERAL generative process (N Bernoulli(p)
#      availability indicators, Exponential(lambda) latencies, K-th order
#      statistic).  It shares no formula with the kernel, so it is what
#      establishes that Math/FRQ.tex's closed form is the right closed form.
#   2. An independent R transcription of the boxed formulas.
#   3. The doc's implementation contract (Sec. 20): F(Inf) = h, S(Inf) = 1 - h,
#      integrate(f) = h.
#   4. That the parameterisation means what it claims: tau is the conditional
#      median, and (alpha, beta, h, tau) inverts to the generating (p, lambda).
#
# The heavier validation -- large-sample Monte Carlo over a full parameter
# grid, and a sampler-based recovery study -- lives in
# WorkingTests/test_frq_validation.R, which is not run by R CMD check.

# ---------------------------------------------------------------------------
# Independent references
# ---------------------------------------------------------------------------

# (alpha, beta, h, tau) -> (p, lambda), transcribed from Math/FRQ.tex Sec. 7.
ref_pl <- function(a, b, h, tau) {
  p <- qbeta(h, a, b)
  u <- qbeta(0.5 * h, a, b)
  c(p = p, lambda = -log1p(-u / p) / tau)
}

ref_d <- function(x, a, b, h, tau) {
  pl <- ref_pl(a, b, h, tau)
  q <- pl[["p"]] * (-expm1(-pl[["lambda"]] * x))
  pl[["p"]] * pl[["lambda"]] * exp(-pl[["lambda"]] * x) *
    q^(a - 1) * (1 - q)^(b - 1) / beta(a, b)
}

ref_p <- function(x, a, b, h, tau) {
  pl <- ref_pl(a, b, h, tau)
  pbeta(pl[["p"]] * (-expm1(-pl[["lambda"]] * x)), a, b)
}

ref_s <- function(x, a, b, h, tau) {
  pl <- ref_pl(a, b, h, tau)
  pbeta(pl[["p"]] * (-expm1(-pl[["lambda"]] * x)), a, b, lower.tail = FALSE)
}

cpp_d <- function(x, a, b, h, tau) EMC2:::dfrq(x, a, b, h, tau)
cpp_p <- function(x, a, b, h, tau) EMC2:::pfrq(x, a, b, h, tau)
cpp_s <- function(x, a, b, h, tau) EMC2:::pfrq(x, a, b, h, tau,
                                               lower_tail = FALSE)

# Eight cells rather than a full crossing: enough to hit both xlogy branches
# (alpha == 1, beta == 1), an interior point, a near-non-defective h, and
# non-integer shapes.
frq_grid <- data.frame(
  a   = c(1,    1,    2,    1.7,  3,    8,    2.5,  1),
  b   = c(1,    3,    1,    2.3,  6,    2.3,  1.4,  6),
  h   = c(0.3,  0.95, 0.75, 0.95, 0.999, 0.6, 0.85, 0.999),
  tau = c(0.15, 0.4,  0.4,  1.2,  0.25, 0.4, 0.35, 0.15))
frq_x <- c(1e-4, 0.01, 0.05, 0.2, 0.5, 1, 2, 5, 20)

# ---------------------------------------------------------------------------
# 1. The literal finite-reservoir process
# ---------------------------------------------------------------------------

test_that("the closed form is the K-th order statistic of a finite reservoir", {
  # THE anchor test.  Nothing here knows the incomplete beta function exists:
  # it builds N potential evidence units, keeps each with probability p, gives
  # the survivors Exponential(lambda) latencies, and takes the K-th smallest.
  # If F(x) = I_{q(x)}(K, N-K+1) were wrong, this is what would catch it.
  skip_on_cran()
  set.seed(20260816)
  nsim <- 2e4
  N <- 7; K <- 3; p <- 0.8; lam <- 2.5
  a <- K; b <- N - K + 1
  avail <- matrix(runif(nsim * N) < p, nrow = nsim)
  lat <- matrix(rexp(nsim * N, lam), nrow = nsim)
  lat[!avail] <- Inf
  Tk <- apply(lat, 1, function(r) sort(r)[K])

  # Reach the kernel through its OWN coordinates: convert (p, lambda) to
  # (h, tau) as the doc says, then hand the kernel (h, tau).  This exercises
  # the parameterisation and the CDF in one step.
  h <- pbeta(p, a, b)
  tau <- -log1p(-qbeta(0.5 * h, a, b) / p) / lam

  probe <- c(0.05, 0.1, 0.2, 0.4, 0.8, 1.5, 3)
  emp <- vapply(probe, function(x) mean(Tk <= x), numeric(1))
  expect_lt(max(abs(emp - cpp_p(probe, a, b, h, tau))), 0.015)

  # The defective mass is a property of the process, not a fitted constant:
  # fewer than K available units means no response, ever.
  expect_lt(abs(mean(is.infinite(Tk)) - (1 - h)), 0.01)
  expect_equal(cpp_p(Inf, a, b, h, tau), h, tolerance = 1e-12)

  # ... and the kernel must recover exactly the (p, lambda) that generated it.
  got <- EMC2:::frq_rate(a, b, h, tau)
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
    worst["d"] <- max(worst["d"], rel(cpp_d(frq_x, g$a, g$b, g$h, g$tau),
                                      ref_d(frq_x, g$a, g$b, g$h, g$tau)))
    worst["p"] <- max(worst["p"], rel(cpp_p(frq_x, g$a, g$b, g$h, g$tau),
                                      ref_p(frq_x, g$a, g$b, g$h, g$tau)))
    worst["s"] <- max(worst["s"], rel(cpp_s(frq_x, g$a, g$b, g$h, g$tau),
                                      ref_s(frq_x, g$a, g$b, g$h, g$tau)))
  }
  expect_lt(worst[["d"]], 1e-11)
  expect_lt(worst[["p"]], 1e-11)
  expect_lt(worst[["s"]], 1e-10)
})

test_that("the log output equals the log of the natural output", {
  for (i in c(1, 4, 5)) {
    g <- frq_grid[i, ]
    expect_equal(EMC2:::dfrq(frq_x, g$a, g$b, g$h, g$tau, log_out = TRUE),
                 log(cpp_d(frq_x, g$a, g$b, g$h, g$tau)), tolerance = 1e-12)
    expect_equal(EMC2:::pfrq(frq_x, g$a, g$b, g$h, g$tau, log_out = TRUE),
                 log(cpp_p(frq_x, g$a, g$b, g$h, g$tau)), tolerance = 1e-12)
    expect_equal(EMC2:::pfrq(frq_x, g$a, g$b, g$h, g$tau, lower_tail = FALSE,
                             log_out = TRUE),
                 log(cpp_s(frq_x, g$a, g$b, g$h, g$tau)), tolerance = 1e-12)
  }
})

# ---------------------------------------------------------------------------
# 3. The doc's implementation contract (Sec. 20)
# ---------------------------------------------------------------------------

test_that("the asymptotic and integral checks of Sec. 20 hold", {
  for (i in seq_len(nrow(frq_grid))) {
    g <- frq_grid[i, ]
    # cdf(+Inf) = h and survival(+Inf) = 1 - h.  These break if the tail is
    # clamped to one, and they are what makes an omission score log(1 - h).
    expect_equal(cpp_p(Inf, g$a, g$b, g$h, g$tau), g$h, tolerance = 1e-12)
    expect_equal(cpp_s(Inf, g$a, g$b, g$h, g$tau), 1 - g$h, tolerance = 1e-9)
    # integrate(f) = h: a separate claim, the density agreeing with the CDF
    # rather than the CDF with itself.
    I <- integrate(function(z) cpp_d(z, g$a, g$b, g$h, g$tau), 0, Inf,
                   rel.tol = 1e-10)$value
    expect_equal(I, g$h, tolerance = 1e-7)
    # ... and at finite t as well, not only in total.
    for (x in c(0.05, 0.3, 1.5)) {
      Ix <- integrate(function(z) cpp_d(z, g$a, g$b, g$h, g$tau), 0, x,
                      rel.tol = 1e-10)$value
      expect_equal(Ix, cpp_p(x, g$a, g$b, g$h, g$tau), tolerance = 1e-7)
    }
    # CDF and survivor are exact complements.
    expect_equal(cpp_p(frq_x, g$a, g$b, g$h, g$tau) +
                   cpp_s(frq_x, g$a, g$b, g$h, g$tau),
                 rep(1, length(frq_x)), tolerance = 1e-12)
  }
})

test_that("the CDF is monotone and the density non-negative and finite", {
  xs <- exp(seq(log(1e-6), log(50), length.out = 300))
  for (i in seq_len(nrow(frq_grid))) {
    g <- frq_grid[i, ]
    d <- cpp_d(xs, g$a, g$b, g$h, g$tau)
    p <- cpp_p(xs, g$a, g$b, g$h, g$tau)
    expect_false(any(is.na(d)) || any(is.na(p)))
    expect_true(all(d >= 0) && all(is.finite(d)))
    expect_true(all(diff(p) >= -1e-14))
    expect_true(all(p >= 0 & p <= g$h + 1e-12))
  }
})

# ---------------------------------------------------------------------------
# 4. The (alpha, beta, h, tau) parameterisation
# ---------------------------------------------------------------------------

test_that("tau is the conditional median decision time", {
  # The whole point of the parameterisation: tau is not a rate, a scale, or a
  # marginal median.  P(T <= tau | T < Inf) = 1/2 means F(tau) = h/2 exactly.
  expect_equal(EMC2:::frq_quantile_level(), 0.5)
  for (i in seq_len(nrow(frq_grid))) {
    g <- frq_grid[i, ]
    expect_equal(cpp_p(g$tau, g$a, g$b, g$h, g$tau), g$h / 2, tolerance = 1e-12)
  }
})

test_that("frq_rate is the exact inverse of h = I_p(alpha, beta)", {
  for (i in seq_len(nrow(frq_grid))) {
    g <- frq_grid[i, ]
    got <- EMC2:::frq_rate(g$a, g$b, g$h, g$tau)
    want <- ref_pl(g$a, g$b, g$h, g$tau)
    expect_equal(unname(got[1, "p"]), want[["p"]], tolerance = 1e-12)
    expect_equal(unname(got[1, "lambda"]), want[["lambda"]], tolerance = 1e-10)
    # h = I_p(alpha, beta) is the defining relation; check the round trip.
    expect_equal(pbeta(unname(got[1, "p"]), g$a, g$b), g$h, tolerance = 1e-10)
  }
  # Vectorises and recycles -- the Ttransform calls it a column at a time.
  v <- EMC2:::frq_rate(c(1, 2, 3), 2, 0.9, 0.4)
  expect_equal(dim(v), c(3L, 2L))
  expect_equal(colnames(v), c("p", "lambda"))
  expect_equal(v[2, ], EMC2:::frq_rate(2, 2, 0.9, 0.4)[1, ])
})

test_that("K = N removes p from the conditional RT distribution", {
  # Math/FRQ.tex Sec. 4's exceptional case.  With beta = 1 (K = N),
  # F(x) = [p G(x)]^N, so F(x)/h = G(x)^N no longer depends on p: availability
  # becomes a pure omission mechanism once you condition on responding.  For
  # beta > 1 it emphatically does not, and that is the claim distinguishing
  # FRQ from an omission-contaminated race.
  a <- 4
  # Hold lambda fixed while h varies, so p is the only thing that moves.
  cond <- function(x, b, h, lam = 2) {
    tau <- -log1p(-qbeta(0.5 * h, a, b) / qbeta(h, a, b)) / lam
    cpp_p(x, a, b, h, tau) / h
  }
  xs <- c(0.1, 0.3, 0.9)
  expect_equal(cond(xs, 1, 0.4), cond(xs, 1, 0.95), tolerance = 1e-10)
  expect_gt(max(abs(cond(xs, 3, 0.4) - cond(xs, 3, 0.95))), 1e-3)
})

test_that("alpha controls the order of the leading edge", {
  # Sec. 5: G(x) ~ lambda x near zero, so f(x) is proportional to x^(alpha-1)
  # and the log-log slope of the density at the leading edge is alpha - 1.
  xs <- c(1e-7, 2e-7)
  for (a in c(1, 1.5, 3, 5.5)) {
    slope <- diff(log(cpp_d(xs, a, 2.5, 0.9, 0.4))) / diff(log(xs))
    expect_equal(slope, a - 1, tolerance = 1e-4)
  }
})

test_that("degenerate parameters are rejected rather than returning garbage", {
  bad <- list(c(0, 2, 0.9, 0.3), c(2, 0, 0.9, 0.3), c(2, 2, 0, 0.3),
              c(2, 2, 1.5, 0.3), c(2, 2, 0.9, 0), c(2, 2, 0.9, -1),
              c(Inf, 2, 0.9, 0.3), c(NA, 2, 0.9, 0.3))
  for (v in bad) {
    expect_equal(cpp_d(0.4, v[1], v[2], v[3], v[4]), 0)
    expect_equal(cpp_p(0.4, v[1], v[2], v[3], v[4]), 0)
    expect_true(is.na(EMC2:::frq_rate(v[1], v[2], v[3], v[4])[1, "p"]))
  }
  # Non-positive decision times: zero density, survivor one.
  expect_equal(cpp_d(c(-1, 0), 2, 3, 0.9, 0.3), c(0, 0))
  expect_equal(cpp_p(c(-1, 0), 2, 3, 0.9, 0.3), c(0, 0))
  expect_equal(cpp_s(c(-1, 0), 2, 3, 0.9, 0.3), c(1, 1))
})

test_that("h at its upper bound stays finite and reaches the non-defective limit", {
  # 1 - 1e-9 is the model's upper bound for h; below it the defective mass is
  # still resolved to a sensible relative accuracy.  (Pushing h closer to one
  # than that loses 1 - p to rounding inside qbeta, which is why the bound is
  # where it is rather than at one.)
  h <- 1 - 1e-9
  expect_equal(cpp_p(Inf, 2, 3, h, 0.3), h, tolerance = 1e-15)
  expect_equal(cpp_s(Inf, 2, 3, h, 0.3), 1 - h, tolerance = 1e-5)
  expect_true(all(is.finite(cpp_d(frq_x, 2, 3, h, 0.3))))
  # h = 1 exactly is the non-defective boundary: p = 1 and S(Inf) = 0.
  expect_equal(unname(EMC2:::frq_rate(2, 3, 1, 0.3)[1, "p"]), 1)
  expect_equal(cpp_p(Inf, 2, 3, 1, 0.3), 1, tolerance = 1e-14)
  expect_equal(cpp_s(Inf, 2, 3, 1, 0.3), 0, tolerance = 1e-14)
})

# ---------------------------------------------------------------------------
# 5. The R model interface
# ---------------------------------------------------------------------------

test_that("the constructor exposes the documented contract", {
  m <- FRQ()
  expect_equal(m$type, "RACE")
  expect_equal(m$c_name, "FRQ")
  # p_types ORDER is the kernel column order (src/col_registry.h); a silent
  # reordering here is a silent wrong answer, so it is pinned literally.
  expect_equal(names(m$p_types),
               c("alpha", "beta", "h", "tau", "t0", "pContaminant", "pGuess"))
  expect_equal(m$p_types_canonical, c("alpha", "beta", "h", "tau", "t0"))
  expect_equal(unname(m$transform$func[c("alpha", "beta", "h", "tau", "t0")]),
               c("exp", "exp", "pnorm", "exp", "exp"))
  # The default must not declare an implausible omission rate for a parameter
  # the user left out of the formula.
  expect_equal(pnorm(m$p_types[["h"]]), 0.95, tolerance = 1e-12)
  # relax only moves the shape bounds; the kernel is shared.
  expect_equal(unname(m$bound$minmax[, "alpha"]), c(1, Inf))
  expect_equal(unname(m$bound$minmax[, "beta"]), c(1, Inf))
  expect_equal(unname(FRQ(relax = TRUE)$bound$minmax[, "alpha"]),
               c(1e-4, Inf))
  expect_equal(FRQ(relax = TRUE)$c_name, "FRQ")
})

test_that("dfun/pfun apply t0 and call the same kernel", {
  pars <- cbind(alpha = 2, beta = 3, h = 0.9, tau = 0.35, t0 = 0.15)
  pars <- pars[rep(1, 5), ]
  rt <- c(0.05, 0.15, 0.3, 0.8, Inf)
  # Below t0 the density and CDF are zero; at rt = Inf the CDF is h, not one.
  expect_equal(EMC2:::dFRQ(rt, pars),
               c(0, 0, cpp_d(c(0.15, 0.65), 2, 3, 0.9, 0.35), 0))
  expect_equal(EMC2:::pFRQ(rt, pars),
               c(0, 0, cpp_p(c(0.15, 0.65), 2, 3, 0.9, 0.35), 0.9))
  expect_equal(EMC2:::sFRQ(rt, pars),
               c(1, 1, cpp_s(c(0.15, 0.65), 2, 3, 0.9, 0.35), 0.1))
})

test_that("Ttransform reports the generative parameters", {
  pars <- cbind(alpha = 2, beta = 3, h = 0.9, tau = 0.35, t0 = 0.15)
  out <- FRQ()$Ttransform(pars, NULL)
  want <- ref_pl(2, 3, 0.9, 0.35)
  expect_equal(unname(out[1, "p"]), want[["p"]], tolerance = 1e-12)
  expect_equal(unname(out[1, "lambda"]), want[["lambda"]], tolerance = 1e-10)
  expect_equal(unname(out[1, "N"]), 4)          # alpha + beta - 1
  expect_equal(unname(out[1, "d"]), 0.5)        # alpha / N
  # Ttransform must not disturb the leading columns the kernel indexes, and
  # must not invent row names that then travel with the parameter matrix.
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

# The independent race likelihood of Sec. 14/20:
#   log L_i(t) = log f_i(t) + sum_{j != i} log S_j(t),
#   log L(omit) = sum_j log S_j(Inf) = sum_j log(1 - h_j).
ref_race_ll <- function(dadm, pars) {
  n_lR <- nlevels(dadm$lR)
  trial <- rep(seq_len(nrow(dadm) / n_lR), each = n_lR)
  d <- EMC2:::dFRQ(dadm$rt, pars)
  s <- EMC2:::sFRQ(dadm$rt, pars)
  win <- as.logical(dadm$winner)
  ll <- vapply(split(seq_len(nrow(dadm)), trial), function(ix) {
    if (!any(win[ix])) return(sum(log(s[ix])))   # omission: rt is Inf
    log(d[ix][win[ix]]) + sum(log(s[ix][!win[ix]]))
  }, numeric(1))
  sum(pmax(ll, log(1e-10)))
}

# h and tau per accumulator under the fixture's ADmat contrast on lM.
frq_truth <- list(hT = pnorm(qnorm(0.93) + 0.45), hF = pnorm(qnorm(0.93) - 0.45))

frq_fixture <- function(n = 60, seed = 20260816) {
  set.seed(seed)
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  dat <- dat[seq_len(n), ]
  ADmat <- matrix(c(-1 / 2, 1 / 2), ncol = 1, dimnames = list(NULL, "d"))
  des <- suppressMessages(design(
    data = dat, model = FRQ, matchfun = matchfun,
    formula = list(alpha ~ 1, beta ~ 1, h ~ lM, tau ~ lM, t0 ~ 1),
    contrasts = list(h = list(lM = ADmat), tau = list(lM = ADmat))))
  p <- c(alpha = log(2), beta = log(3), h = qnorm(0.93), h_lMd = 0.9,
         tau = log(0.35), tau_lMd = -0.5, t0 = log(0.15))
  list(dat = dat, des = des, p = p[names(sampled_pars(des))],
       emc = suppressMessages(make_emc(dat, des, type = "single", n_chains = 1,
                                       compress = FALSE, rt_resolution = NULL)))
}

test_that("the compiled race likelihood matches the R reference", {
  skip_on_cran()
  fx <- frq_fixture()
  got <- frq_ll(fx$emc, fx$p)
  expect_true(is.finite(got))
  dadm <- fx$emc[[1]]$data[[1]]
  pars <- EMC2:::get_pars_matrix_oo(fx$p, dadm, fx$emc[[1]]$model())
  expect_equal(got, ref_race_ll(dadm, pars), tolerance = 1e-6)
})

test_that("a p_types reordering is caught by the column contract", {
  # The kernel indexes alpha/beta/h/tau/t0 positionally, so a silent
  # reordering would make it read tau out of h's slot.  validate_col_prefix
  # must stop rather than return a plausible number.
  skip_on_cran()
  fx <- frq_fixture(n = 30)
  swapped <- names(FRQ()$p_types)
  swapped[3:4] <- swapped[4:3]
  expect_error(frq_ll(fx$emc, fx$p, p_types_override = swapped),
               "FRQ kernels expect parameter column")
})

test_that("an omission scores exactly sum_j log(1 - h_j)", {
  # Sec. 15 and Sec. 20's no_response contract.  With an infinite observation
  # window an omission's likelihood is the product of the accumulators'
  # never-finish masses -- true only if the survivor saturates at 1 - h
  # instead of falling to zero.
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
  # Dropping the six omissions must remove exactly 6 * per_omission.
  expect_equal(ll_om - frq_ll(mk(d_om[-(1:6), ]), fx$p), 6 * per_omission,
               tolerance = 1e-6)
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

  # pContaminant is a Bernoulli omission rate multiplying the whole
  # likelihood, and FRQ's intrinsic 1 - h defect is a separate mechanism.  On
  # data without omissions the mixture must shift the total by exactly
  # n * log(1 - pc) -- if the two were confounded it would not.
  matchfun <- function(d) d$S == d$lR
  ADmat <- matrix(c(-1 / 2, 1 / 2), ncol = 1, dimnames = list(NULL, "d"))
  des_pc <- suppressMessages(design(
    data = fx$dat, model = FRQ, matchfun = matchfun,
    formula = list(alpha ~ 1, beta ~ 1, h ~ lM, tau ~ lM, t0 ~ 1,
                   pContaminant ~ 1),
    contrasts = list(h = list(lM = ADmat), tau = list(lM = ADmat))))
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
  p[["t0"]] <- log(max(fx$dat$rt) + 1)   # every trial now has density zero
  ll <- frq_ll(fx$emc, p)
  expect_false(is.na(ll))
  expect_equal(ll, nrow(fx$dat) * log(1e-10), tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# 7. Simulation
# ---------------------------------------------------------------------------

test_that("the C++ and R simulators agree distributionally with the CDF", {
  skip_on_cran()
  n <- 4000
  lR <- factor(rep(c("left", "right"), n), levels = c("left", "right"))
  pars <- cbind(alpha = 2, beta = 3, h = 0.9, tau = 0.35, t0 = 0.1)
  pars <- pars[rep(1, length(lR)), ]
  set.seed(4)
  a <- EMC2:::rFRQ(lR, pars)
  set.seed(4)
  b <- withr::with_options(list(emc2.cpp_rfun = TRUE),
                           EMC2:::.rfun_FRQ(lR, pars))
  probe <- c(0.2, 0.35, 0.5, 0.8)
  # Race CDF for two identical accumulators: 1 - S(x)^2, with S saturating at
  # 1 - h.  Distributional, not stream, equivalence -- as for every rfun port.
  th <- 1 - cpp_s(probe - 0.1, 2, 3, 0.9, 0.35)^2
  for (dat in list(a, b)) {
    fin <- is.finite(dat$rt)
    emp <- vapply(probe, function(x) mean(fin & dat$rt <= x), numeric(1))
    expect_lt(max(abs(emp - th)), 0.03)
    # Race-level omission rate is prod(1 - h_i) = (1 - h)^2 (Sec. 15).
    expect_lt(abs(mean(!fin) - 0.01), 0.006)
    expect_true(all(is.na(dat$R[!fin])))
    expect_true(all(dat$rt[fin] > 0.1))
  }
})

test_that("the simulator reproduces the literal reservoir at integer shapes", {
  # A second route to the Sec. 6 claim, now through the package's simulator
  # rather than its density: for integer alpha = K and beta = N - K + 1 the
  # Beta-quorum draw must be indistinguishable from building N cues.
  skip_on_cran()
  set.seed(99)
  N <- 6; K <- 2; p <- 0.7; lam <- 3
  a <- K; b <- N - K + 1
  h <- pbeta(p, a, b)
  tau <- -log1p(-qbeta(0.5 * h, a, b) / p) / lam
  n <- 1e4
  lR <- factor(rep("go", n), levels = "go")
  pars <- cbind(alpha = a, beta = b, h = h, tau = tau, t0 = 0)
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
  sim <- make_data(fx$p, design = fx$des, n_trials = 150)
  expect_true(any(is.infinite(sim$rt)))              # intrinsic omissions
  expect_true(all(is.na(sim$R[is.infinite(sim$rt)])))
  fin <- is.finite(sim$rt)
  expect_true(all(sim$rt[fin] > 0.15))
  # h_lMd and tau_lMd both favour the matching accumulator.
  expect_gt(mean(sim$S[fin] == sim$R[fin]), 0.6)
  e <- suppressMessages(make_emc(sim, fx$des, type = "single", n_chains = 1,
                                 compress = FALSE, rt_resolution = NULL))
  expect_true(is.finite(frq_ll(e, fx$p)))
})

# ---------------------------------------------------------------------------
# 8. Recovery
# ---------------------------------------------------------------------------

test_that("the likelihood is maximised at the generating parameters", {
  # A cheap, deterministic stand-in for a recovery study: simulate, then check
  # the profile likelihood in each estimated coordinate peaks at the truth.  A
  # sign error, a swapped column, or an h/tau confusion all show up here
  # without paying for a sampler run.  (The full sampler-based recovery lives
  # in WorkingTests/test_frq_validation.R.)
  skip_on_cran()
  set.seed(20260816)
  fx <- frq_fixture(n = 30)
  # 3000 trials, not fewer: beta is by far the weakest coordinate (a +-0.3
  # perturbation costs ~1.4 nats against ~40-2000 for the others, which is the
  # identifiability note in ?FRQ showing up quantitatively), so at a few
  # hundred trials sampling noise flips its profile.
  sim <- make_data(fx$p, design = fx$des, n_trials = 3000)
  e <- suppressMessages(make_emc(sim, fx$des, type = "single", n_chains = 1,
                                 compress = FALSE, rt_resolution = NULL))
  ll0 <- frq_ll(e, fx$p)
  expect_true(is.finite(ll0))
  for (nm in names(fx$p)) {
    for (delta in c(-0.3, 0.3)) {
      p2 <- fx$p
      p2[[nm]] <- p2[[nm]] + delta
      expect_lt(frq_ll(e, p2), ll0)
    }
  }
})
