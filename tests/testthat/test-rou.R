# Racing Ornstein-Uhlenbeck (leaky accumulator) -- src/fpe_race.h,
# src/model_ROU.h, R/model_ROU.R.
#
# The organising principle of this file is the k = 0 oracle.  At zero leak the
# ROU is the racing diffusion model, but it gets there through the
# Fokker-Planck solver rather than by dispatching to the Wald formulae, so every
# comparison against RDM below exercises the whole path -- design, Ttransform,
# column order, cache, kernels, race combination -- against an independently
# trusted implementation.  See race_ou_integration_plan.md section 4.3.

wald_pdf <- function(t, mu, b) b / sqrt(2 * pi * t^3) * exp(-(b - mu * t)^2 / (2 * t))
wald_cdf <- function(t, mu, b) {
  pnorm((mu * t - b) / sqrt(t)) + exp(2 * mu * b) * pnorm(-(b + mu * t) / sqrt(t))
}

rou_pars <- function(t, v, k, B, A, t0 = 0, s = 1) {
  n <- length(t)
  cbind(v = rep(v, n), k = rep(k, n), B = rep(B, n), A = rep(A, n),
        t0 = rep(t0, n), s = rep(s, n))
}

tq <- seq(0.05, 2.0, by = 0.05)

test_that("k = 0 reproduces the Wald through the PDE, not around it", {
  p <- rou_pars(tq, v = 1.5, k = 0, B = 1, A = 0)
  expect_lt(max(abs(EMC2:::dROU(tq, p) - wald_pdf(tq, 1.5, 1))), 3e-3)
  expect_lt(max(abs(EMC2:::pROU(tq, p) - wald_cdf(tq, 1.5, 1))), 3e-4)

  # with start-point variability, against the package's own RDM
  p2 <- rou_pars(tq, v = 1.5, k = 0, B = 1, A = 0.5)
  rdm <- cbind(v = 1.5, B = 1, A = 0.5, t0 = 0, s = 1)[rep(1, length(tq)), ]
  expect_lt(max(abs(EMC2:::dROU(tq, p2) - EMC2:::dRDM(tq, rdm))), 3e-3)
  expect_lt(max(abs(EMC2:::pROU(tq, p2) - EMC2:::pRDM(tq, rdm))), 3e-4)
})

test_that("s is scaled out exactly, as in the RDM", {
  # (v, B, A) -> (2v, 2B, 2A) with s -> 2s is the same process observed in a
  # rescaled state, so the first-passage density is unchanged.  k has units of
  # 1/time and must NOT be rescaled; if it were, this would fail at k > 0.
  a <- EMC2:::dROU(tq, rou_pars(tq, v = 1.5, k = 2, B = 1, A = 0.5, s = 1))
  b <- EMC2:::dROU(tq, rou_pars(tq, v = 3.0, k = 2, B = 2, A = 1.0, s = 2))
  expect_identical(a, b)
})

test_that("the leak slows the race and creates a long finite tail", {
  # theta = v/k is the asymptote.  With v = 1.5 and b = 1, k = 3 puts the
  # asymptote at 0.5, BELOW threshold, so finishing depends on noise alone.
  cdf2 <- vapply(c(0, 1, 3), function(kk) {
    EMC2:::pROU(2.0, rou_pars(2.0, v = 1.5, k = kk, B = 1, A = 0))
  }, numeric(1))
  expect_true(all(diff(cdf2) < 0))     # more leak, less mass by 2 s
  expect_gt(cdf2[1], 0.95)
  expect_lt(cdf2[3], 0.85)
  expect_identical(EMC2:::pROU(Inf, rou_pars(Inf, v = 1.5, k = 3, B = 1, A = 0)), 1)

  # and the density must stay a density
  d <- EMC2:::dROU(tq, rou_pars(tq, v = 1.5, k = 3, B = 1, A = 0.5))
  expect_true(all(d >= 0))
  expect_false(anyNA(d))
})

test_that("the reference simulator agrees with the solved cdf", {
  # The simulator steps the EXACT OU transition (k -> 0 taken analytically) with
  # a Brownian-bridge crossing correction, and shares nothing with the solver
  # but the parameterisation.  MC error at 2e5 draws is ~1e-3.
  skip_on_cran()
  set.seed(20260729)
  probe <- c(0.2, 0.4, 0.7, 1.0, 1.5, 2.0)
  N <- 2e5
  for (kk in c(0, 2)) {
    ht <- EMC2:::rou_hit_times_vec(rep(1.5, N), rep(kk, N), rep(1, N),
                                   rep(0.5, N), rep(1, N), 5e-4, 10)
    emp <- vapply(probe, function(x) mean(ht <= x), numeric(1))
    th <- EMC2:::pROU(probe, rou_pars(probe, v = 1.5, k = kk, B = 1, A = 0.5))
    expect_lt(max(abs(emp - th)), 5e-3)
  }
})

test_that("one solve serves every row sharing a parameter tuple", {
  # This is the whole reason a PDE model is affordable here: cost tracks the
  # number of distinct parameter rows, not the number of trials.
  n <- 500
  r <- EMC2:::rou_pdf_cdf_vec(seq(0.2, 2, length.out = n), rep(1.5, n), rep(1, n),
                              rep(1, n), rep(0.5, n), rep(0, n), rep(1, n))
  expect_identical(r$n_solves, 1L)

  v <- rep(c(1.0, 1.5, 2.0), length.out = n)
  r2 <- EMC2:::rou_pdf_cdf_vec(seq(0.2, 2, length.out = n), v, rep(1, n),
                               rep(1, n), rep(0.5, n), rep(0, n), rep(1, n))
  expect_identical(r2$n_solves, 3L)
})

test_that("batched fixed-boundary solves reproduce scalar solves", {
  # Eight distinct tuples exercise one full AVX-512 batch when available and
  # two four-lane batches otherwise. Mixed horizons and starts verify lane-local
  # schedules rather than merely checking the equal-clock case.
  rt <- c(0.05, 0.2, 0.4, 0.65, 0.9, 1.2, 1.55, 2)
  v <- seq(1, 2.4, length.out = 8)
  k <- seq(0.5, 1.9, length.out = 8)
  A <- rep(c(0, 0.5), 4)
  pars <- cbind(v = v, k = k, B = 1, A = A, t0 = 0, s = 1)

  batched_pdf <- EMC2:::dROU(rt, pars)
  batched_cdf <- EMC2:::pROU(rt, pars)
  scalar_pdf <- vapply(seq_along(rt), function(i) {
    EMC2:::dROU(rt[i], pars[i, , drop = FALSE])
  }, numeric(1))
  scalar_cdf <- vapply(seq_along(rt), function(i) {
    EMC2:::pROU(rt[i], pars[i, , drop = FALSE])
  }, numeric(1))

  expect_equal(batched_pdf, scalar_pdf, tolerance = 1e-13)
  expect_equal(batched_cdf, scalar_cdf, tolerance = 1e-13)
})

test_that("sparse finite-time output reproduces the full solver grid", {
  rt <- c(0.07, 0.18, 0.31, 0.52, 0.78, 1.05, 1.42, 1.9)
  pars <- cbind(
    v = seq(1.0, 2.4, length.out = 8),
    k = seq(0.4, 1.8, length.out = 8),
    B = rep(c(0.9, 1.1), 4),
    A = rep(c(0, 0.35), 4),
    t0 = rep(c(0, 0.04), 4),
    s = 1
  )

  withr::local_options(emc2.rou_sparse_output = FALSE)
  full_pdf <- EMC2:::dROU(rt, pars)
  full_cdf <- EMC2:::pROU(rt, pars)
  options(emc2.rou_sparse_output = TRUE)
  sparse_pdf <- EMC2:::dROU(rt, pars)
  sparse_cdf <- EMC2:::pROU(rt, pars)

  expect_equal(sparse_pdf, full_pdf, tolerance = 1e-13)
  expect_equal(sparse_cdf, full_cdf, tolerance = 1e-13)
})

test_that(".rfun_ROU respects emc2.cpp_rfun option and falls back to R simulator", {
  lR <- factor(c("left", "right"))
  pars <- matrix(c(v = 1.5, k = 1, B = 1, A = 0.5, t0 = 0.2, s = 1,
                   v = 1.0, k = 1, B = 1, A = 0.5, t0 = 0.2, s = 1),
                 nrow = 2, byrow = TRUE,
                 dimnames = list(NULL, c("v", "k", "B", "A", "t0", "s")))

  set.seed(42)
  res_cpp <- EMC2:::.rfun_ROU(lR, pars, ok = rep(TRUE, 2))
  expect_s3_class(res_cpp, "data.frame")
  expect_named(res_cpp, c("R", "rt"))

  op <- options(emc2.cpp_rfun = FALSE)
  on.exit(options(op), add = TRUE)

  set.seed(42)
  res_r <- EMC2:::.rfun_ROU(lR, pars, ok = rep(TRUE, 2))
  expect_s3_class(res_r, "data.frame")
  expect_named(res_r, c("R", "rt"))
})

test_that(".rfun_ROU applies timed winners on the C++ path", {
  withr::local_options(emc2.cpp_rfun = TRUE)
  lR <- factor(c("left", "time"), levels = c("left", "time"))
  pars <- matrix(
    c(0, 0, 100, 0, 0, 1,
      0, 0,   0, 0, 0, 1),
    nrow = 2, byrow = TRUE,
    dimnames = list(NULL, c("v", "k", "B", "A", "t0", "s"))
  )

  out <- EMC2:::.rfun_ROU(lR, pars)
  expect_identical(as.character(out$R), "left")
  expect_identical(out$isTime, TRUE)
  expect_identical(out$rt, 0)
})

test_that(".rfun_ROU forwards simulator resolution and horizon options", {
  withr::local_options(list(
    emc2.cpp_rfun = TRUE,
    emc2.rou_sim_dt = 1e-5,
    emc2.rou_sim_tmax = 1e-4
  ))
  lR <- factor("left")
  pars <- matrix(
    c(v = 1, k = 0, B = 0.01, A = 0, t0 = 0, s = 1e-12),
    nrow = 1, dimnames = list(NULL, c("v", "k", "B", "A", "t0", "s"))
  )

  out <- EMC2:::.rfun_ROU(lR, pars)
  expect_true(is.na(out$R))
  expect_identical(out$rt, Inf)
})


# ---------------------------------------------------------------------------
# Package-level: design -> make_emc -> C++ likelihood.
# ---------------------------------------------------------------------------

rou_ll_fixture <- function() {
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  ADmat <- matrix(c(-1 / 2, 1 / 2), ncol = 1, dimnames = list(NULL, "d"))

  mk <- function(model, extra_constants) {
    des <- design(data = dat, model = model, matchfun = matchfun,
                  formula = if (identical(model, RDM)) {
                    list(v ~ lM, B ~ E + lR, A ~ 1, t0 ~ 1)
                  } else {
                    list(v ~ lM, k ~ 1, B ~ E + lR, A ~ 1, t0 ~ 1)
                  },
                  contrasts = list(v = list(lM = ADmat)),
                  constants = extra_constants)
    # rt_resolution MUST be pinned to NULL on both sides.  ROU declares
    # compress_ok = FALSE, so make_emc() forces rt_resolution = NULL for it
    # regardless; leaving the RDM on the 1/60 default would floor its rt and
    # nothing else, and the ROU-vs-RDM comparisons below would then be measuring
    # the flooring bias (~0.15 per trial) rather than the solver's
    # discretisation error -- which refining the mesh makes worse, not better.
    list(des = des,
         emc = make_emc(dat, des, type = "single", n_chains = 1, compress = TRUE,
                        rt_resolution = NULL))
  }
  list(dat = dat,
       rdm = mk(RDM, c(s = log(1))),
       rou0 = mk(ROU, c(s = log(1), k = log(0))),
       rouk = mk(ROU, c(s = log(1))))
}

rou_ll <- function(fx, p) {
  model <- fx$emc[[1]]$model()
  p_types <- names(model$p_types)
  dadm <- fx$emc[[1]]$data[[1]]
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

test_that("ROU at k = 0 reproduces the RDM log-likelihood on forstmann", {
  skip_on_cran()
  fx <- suppressMessages(rou_ll_fixture())
  pv <- c(v = log(1.6), v_lMd = 0.9, B = log(1.2), B_Eneutral = 0.1,
          B_Eaccuracy = 0.2, B_lRright = 0.05, A = log(0.4), t0 = log(0.2))

  l_rdm <- rou_ll(fx$rdm, pv[names(sampled_pars(fx$rdm$des))])
  l_rou <- rou_ll(fx$rou0, pv[names(sampled_pars(fx$rou0$des))])
  expect_true(is.finite(l_rdm) && is.finite(l_rou))
  # Discretisation error, per trial, at the shipped default resolution.
  expect_lt(abs(l_rou - l_rdm) / nrow(fx$dat), 0.01)
})

test_that("the ROU likelihood converges on the RDM as the grid is refined", {
  # The residual above must be DISCRETISATION error, which means it has to
  # shrink when the mesh and time step are refined.  If it did not, the two
  # models would simply be different -- an accuracy claim that is not tested by
  # refinement is not an accuracy claim at all.
  skip_on_cran()
  fx <- suppressMessages(rou_ll_fixture())
  pv <- c(v = log(1.6), v_lMd = 0.9, B = log(1.2), B_Eneutral = 0.1,
          B_Eaccuracy = 0.2, B_lRright = 0.05, A = log(0.4), t0 = log(0.2))
  l_rdm <- rou_ll(fx$rdm, pv[names(sampled_pars(fx$rdm$des))])
  p_rou <- pv[names(sampled_pars(fx$rou0$des))]

  old <- options(emc2.fpe_nx = 128L, emc2.fpe_dt = 8e-3)
  on.exit(options(old), add = TRUE)
  coarse <- abs(rou_ll(fx$rou0, p_rou) - l_rdm)
  options(emc2.fpe_nx = 512L, emc2.fpe_dt = 5e-4)
  fine <- abs(rou_ll(fx$rou0, p_rou) - l_rdm)

  expect_lt(fine, coarse)
  expect_lt(fine, 0.25 * coarse)
})

test_that("a free leak moves the likelihood and stays finite", {
  skip_on_cran()
  fx <- suppressMessages(rou_ll_fixture())
  pv <- c(v = log(1.6), v_lMd = 0.9, B = log(1.2), B_Eneutral = 0.1,
          B_Eaccuracy = 0.2, B_lRright = 0.05, A = log(0.4), t0 = log(0.2))
  nm <- names(sampled_pars(fx$rouk$des))
  lls <- vapply(c(-20, log(0.5), log(2)), function(kv) {
    rou_ll(fx$rouk, c(pv, k = kv)[nm])
  }, numeric(1))
  expect_true(all(is.finite(lls)))
  expect_true(all(diff(lls) < 0))          # forstmann is not a leaky data set
  # k -> 0 must agree with the model that has k fixed at zero
  expect_equal(lls[1], rou_ll(fx$rou0, pv[names(sampled_pars(fx$rou0$des))]),
               tolerance = 1e-8)
})

test_that("ROU handles censored and truncated designs", {
  # Unknown-winner censoring needs only scalar CDF calls at the censoring and
  # truncation bounds (log_surv_cm, particle_ll.cpp), which the solve cache
  # serves by interpolation.  Check against RDM at k = 0 so the comparison is
  # against a trusted answer rather than against itself.
  skip_on_cran()
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  # Censor the slowest responses: rt above the cut is known only to be above it.
  cut <- as.numeric(quantile(dat$rt, 0.9))
  upper <- max(dat$rt) + 0.5
  dat$rt[dat$rt > cut] <- Inf
  dat$LT <- 0.1
  dat$UC <- cut
  dat$UT <- upper
  ADmat <- matrix(c(-1 / 2, 1 / 2), ncol = 1, dimnames = list(NULL, "d"))

  mk <- function(model, form, consts) {
    des <- design(data = dat, model = model, matchfun = matchfun, formula = form,
                  contrasts = list(v = list(lM = ADmat)), constants = consts)
    list(des = des, emc = suppressMessages(
      # NULL on both sides: ROU forces it (compress_ok = FALSE), so asking for
      # 0.02 would bin the RDM's rt and not ROU's, and the comparison below
      # would be measuring that difference instead of the solver's accuracy.
      make_emc(dat, des, type = "single", n_chains = 1, compress = TRUE,
               rt_resolution = NULL)))
  }
  fr <- suppressMessages(mk(RDM, list(v ~ lM, B ~ E + lR, A ~ 1, t0 ~ 1),
                            c(s = log(1))))
  fo <- suppressMessages(mk(ROU, list(v ~ lM, k ~ 1, B ~ E + lR, A ~ 1, t0 ~ 1),
                            c(s = log(1), k = log(0))))

  pv <- c(v = log(1.6), v_lMd = 0.9, B = log(1.2), B_Eneutral = 0.1,
          B_Eaccuracy = 0.2, B_lRright = 0.05, A = log(0.4), t0 = log(0.2))
  l_rdm <- rou_ll(fr, pv[names(sampled_pars(fr$des))])
  l_rou <- rou_ll(fo, pv[names(sampled_pars(fo$des))])
  expect_true(is.finite(l_rdm) && is.finite(l_rou))
  expect_lt(abs(l_rou - l_rdm) / nrow(dat), 0.01)

  withr::local_options(emc2.rou_sparse_output = FALSE)
  full_grid <- rou_ll(fo, pv[names(sampled_pars(fo$des))])
  options(emc2.rou_sparse_output = TRUE)
  sparse <- rou_ll(fo, pv[names(sampled_pars(fo$des))])
  expect_equal(sparse, full_grid, tolerance = 1e-12)
})

test_that("rROU simulates data the model can then fit back", {
  skip_on_cran()
  set.seed(20260729)
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  ADmat <- matrix(c(-1 / 2, 1 / 2), ncol = 1, dimnames = list(NULL, "d"))
  des <- suppressMessages(design(
    data = dat, model = ROU, matchfun = matchfun,
    formula = list(v ~ lM, k ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
    contrasts = list(v = list(lM = ADmat)), constants = c(s = log(1))))
  p <- c(v = log(2), v_lMd = 1, k = log(1), B = log(1), A = log(0.3),
         t0 = log(0.2))[names(sampled_pars(des))]
  sim <- make_data(p, design = des, n_trials = 40)

  expect_true(all(is.finite(sim$rt)))
  expect_true(all(sim$rt > 0.2))
  expect_false(anyNA(sim$R))
  # Accuracy must exceed chance: v_lMd = 1 favours the matching accumulator.
  expect_gt(mean(sim$S == sim$R), 0.6)
})

# ============================================================================
# Collapsing boundaries (step 12).  The whole point of a PDE likelihood is that
# it can do what the analytic race models cannot, so these are not an extension
# to be checked loosely.
# ============================================================================

rou_bnd_pars <- function(n, v = 1.5, k = 0.5, B = 1, A = 0.4, Binf = 0.3,
                         tau = 0.4, pw = NULL) {
  p <- cbind(v = rep(v, n), k = rep(k, n), B = rep(B, n), A = rep(A, n),
             t0 = rep(0, n), s = rep(1, n), Binf = rep(Binf, n),
             tau = rep(tau, n))
  if (!is.null(pw)) p <- cbind(p, pw = rep(pw, n))
  p
}

test_that("a collapse to Binf == b0 is exactly the fixed-bound model", {
  # Not a nicety: the degeneracy test in FPE_Boundary::set_kind is what restores
  # the one-time operator factorisation, so if this drifts the model silently
  # pays 3x for a boundary that does not move.
  rt <- seq(0.05, 2, length.out = 120)
  n <- length(rt)
  fixed <- cbind(v = rep(1.5, n), k = rep(0.5, n), B = rep(1, n),
                 A = rep(0.4, n), t0 = rep(0, n), s = rep(1, n))
  for (kind in c("exponential", "linear_additive", "linear_multiplicative", "weibull")) {
    pw <- if (kind == "weibull") 1.5 else NULL
    coll <- rou_bnd_pars(n, Binf = 1.4, pw = pw)   # Binf == b0 = B + A
    expect_identical(EMC2:::dROU(rt, coll, kind = kind), EMC2:::dROU(rt, fixed))
    expect_identical(EMC2:::pROU(rt, coll, kind = kind), EMC2:::pROU(rt, fixed))
  }
})

test_that("a collapsing bound finishes strictly sooner, at every t", {
  rt <- seq(0.05, 3, length.out = 150)
  n <- length(rt)
  fixed <- cbind(v = rep(1.5, n), k = rep(0.5, n), B = rep(1, n),
                 A = rep(0.4, n), t0 = rep(0, n), s = rep(1, n))
  for (kind in c("exponential", "linear_additive", "linear_multiplicative", "weibull")) {
    pw <- if (kind == "weibull") 1.5 else NULL
    coll <- rou_bnd_pars(n, pw = pw)
    cf <- EMC2:::pROU(rt, fixed); cc <- EMC2:::pROU(rt, coll, kind = kind)
    expect_true(all(cc >= cf - 1e-9), info = kind)
    expect_gt(max(cc - cf), 0.05)
    # A slow accumulator gains finite-time response mass when the bound moves.
    expect_gt(max(cc), max(cf))
  }
})

test_that("a very slow collapse converges back to the fixed bound", {
  rt <- seq(0.05, 1.5, length.out = 80)
  n <- length(rt)
  fixed <- cbind(v = rep(1.5, n), k = rep(0.5, n), B = rep(1, n),
                 A = rep(0.4, n), t0 = rep(0, n), s = rep(1, n))
  cf <- EMC2:::pROU(rt, fixed)
  err <- sapply(c(5, 50, 500), function(tau)
    max(abs(EMC2:::pROU(rt, rou_bnd_pars(n, tau = tau), kind = "exponential") - cf)))
  expect_true(all(diff(err) < 0))
  expect_lt(err[3], 5e-3)
})

test_that("the collapsing solver agrees with the reference simulator", {
  skip_on_cran()
  set.seed(20260730)
  N <- 1e5
  for (kind in c("exponential", "linear_additive", "linear_multiplicative", "weibull")) {
    pw <- if (kind == "weibull") 1.5 else NULL
    ps <- rou_bnd_pars(N, pw = pw)
    bk <- c(exponential = 2L, linear_additive = 3L, linear_multiplicative = 4L, weibull = 1L)[[kind]]
    ht <- EMC2:::rou_hit_times_vec(ps[, "v"], ps[, "k"], ps[, "B"], ps[, "A"],
                            ps[, "s"], 2e-4, 30, bk, ps[, "Binf"],
                            ps[, "tau"], if (is.null(pw)) numeric(0) else ps[, "pw"])
    tt <- seq(0.05, 2.5, length.out = 40)
    emp <- sapply(tt, function(x) mean(ht <= x))
    sol <- EMC2:::pROU(tt, rou_bnd_pars(length(tt), pw = pw), kind = kind)
    # MC standard error is ~1/sqrt(N) = 3.2e-3; allow twice that.
    expect_lt(max(abs(emp - sol)), 7e-3, label = kind)
  }
})

test_that("a collapsing-bound design fits through the sampled likelihood", {
  skip_on_cran()
  set.seed(20260731)
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  des <- suppressMessages(design(
    data = dat, model = function() ROU(boundary_collapse = "exponential"),
    matchfun = matchfun,
    formula = list(v ~ 1, k ~ 1, B ~ 1, A ~ 1, t0 ~ 1, Binf ~ 1, tau ~ 1),
    constants = c(s = log(1))))
  nm <- names(sampled_pars(des))
  expect_true(all(c("Binf", "tau") %in% nm))

  invisible(capture.output(emc <- make_emc(dat, des, type = "single",
                                           n_chains = 1, compress = TRUE)))
  model <- emc[[1]]$model()
  expect_identical(model$c_name, "ROU_BEXP")
  dadm <- emc[[1]]$data[[1]]
  p_types <- names(model$p_types)
  designs <- lapply(p_types, function(x)
    attr(dadm, "designs")[[x]][attr(attr(dadm, "designs")[[x]], "expand"), , drop = FALSE])
  names(designs) <- p_types
  ll_of <- function(p) EMC2:::calc_ll_oo(
    matrix(p, nrow = 1, dimnames = list(NULL, nm)), dadm,
    constants = attr(dadm, "constants"), designs = designs,
    type = model$c_name, bounds = model$bound, transforms = model$transform,
    pretransforms = model$pre_transform, p_types = p_types,
    min_ll = log(1e-10), trend = model$trend)

  p <- c(v = log(2), k = log(0.5), B = log(1.2), A = log(0.4), t0 = log(0.2),
         Binf = log(0.5), tau = log(0.5))[nm]
  expect_true(is.finite(ll_of(p)))
  # Binf == B is the fixed bound, and must reproduce the fixed-bound model's
  # likelihood exactly -- through a different code path, with the shape columns
  # present but degenerate.
  # b0 = B + A, and Binf is measured from zero, so the degenerate value is B + A.
  p_deg <- p; p_deg[["Binf"]] <- log(exp(p[["B"]]) + exp(p[["A"]]))
  des_f <- suppressMessages(design(
    data = dat, model = ROU, matchfun = matchfun,
    formula = list(v ~ 1, k ~ 1, B ~ 1, A ~ 1, t0 ~ 1), constants = c(s = log(1))))
  invisible(capture.output(emc_f <- make_emc(dat, des_f, type = "single",
                                             n_chains = 1, compress = TRUE)))
  mf <- emc_f[[1]]$model(); pt_f <- names(mf$p_types); df <- emc_f[[1]]$data[[1]]
  dg <- lapply(pt_f, function(x)
    attr(df, "designs")[[x]][attr(attr(df, "designs")[[x]], "expand"), , drop = FALSE])
  names(dg) <- pt_f
  nmf <- names(sampled_pars(des_f))
  ll_f <- EMC2:::calc_ll_oo(matrix(p[nmf], nrow = 1, dimnames = list(NULL, nmf)),
    df, constants = attr(df, "constants"), designs = dg, type = mf$c_name,
    bounds = mf$bound, transforms = mf$transform,
    pretransforms = mf$pre_transform, p_types = pt_f, min_ll = log(1e-10),
    trend = mf$trend)
  expect_equal(ll_of(p_deg), ll_f)
  # And a real collapse must actually change the likelihood.
  expect_false(isTRUE(all.equal(ll_of(p), ll_f)))

  sim <- make_data(p, design = des, n_trials = 40)
  expect_true(all(is.finite(sim$rt)))
  expect_false(anyNA(sim$R))
})

test_that("the collapse asymptote is measured from zero, and can reach it", {
  # The anchor: with a wide start range and a fast collapse to Binf = 0.1, the
  # bound sweeps DOWN THROUGH [0, A] and nearly everything finishes early.  Were
  # the asymptote Binf + A = 0.9 instead, the bound would still sit above most
  # start points and almost nothing would.
  rt <- c(0.05, 0.1, 0.2, 0.4, 1.0); n <- length(rt)
  cdf_at <- function(A, Binf, tau = 0.05, bk = 2L, pw = numeric(0))
    EMC2:::rou_pdf_cdf_vec(rt, rep(0.05, n), rep(3, n), rep(1, n), rep(A, n),
                           rep(0, n), rep(1, n), 384L, 2e-3, 8, 32, bk,
                           rep(Binf, n), rep(tau, n), pw)$cdf
  expect_gt(cdf_at(0.8, 0.1)[2], 0.4)
  expect_lt(cdf_at(0.8, 0.9)[2], 0.1)

  # Binf = 0 is an ordinary interior value: the domain is [x_lo, b(t)] with
  # x_lo < 0, so a bound reaching the start point does not degenerate it.
  for (bk in c(1L, 2L, 3L, 4L)) {
    pw <- if (bk == 1L) rep(1.5, n) else numeric(0)
    c0 <- cdf_at(0.4, 0, tau = 0.4, bk = bk, pw = pw)
    expect_true(all(is.finite(c0)))
    expect_true(all(diff(c0) >= -1e-12))
    # and it is the continuous limit of Binf -> 0, not a special case
    expect_lt(max(abs(cdf_at(0.4, 1e-6, tau = 0.4, bk = bk, pw = pw) - c0)), 1e-6)
  }
  expect_identical(unname(ROU(boundary_collapse = "exponential")$bound$minmax[, "Binf"]),
                   c(0, Inf))
})

test_that("the likelihood is invariant to the scale of (v, s, b, A)", {
  # dX = (v - kX)dt + s dW absorbed at b is unchanged by X -> cX, which sends
  # (v, s, b, A) -> (cv, cs, cb, cA) and leaves k and t0 alone.  This is why
  # only one of v, s and b can be identified as a scale, and why fixing B while
  # freeing s recovers exactly the model with s fixed and B free.
  rt <- c(0.25, 0.4, 0.7, 1.2, 2)
  n  <- length(rt)
  mk <- function(cc) cbind(v = rep(2.5 * cc, n), k = rep(1.5, n),
                           B = rep(1.0 * cc, n), A = rep(0.3 * cc, n),
                           t0 = rep(0.2, n), s = rep(1.0 * cc, n))
  d1 <- EMC2:::dROU(rt, mk(1)); p1 <- EMC2:::pROU(rt, mk(1))
  expect_true(all(is.finite(d1)) && any(d1 > 0))
  for (cc in c(0.5, 1.7, 3)) {
    expect_equal(EMC2:::dROU(rt, mk(cc)), d1, tolerance = 1e-10)
    expect_equal(EMC2:::pROU(rt, mk(cc)), p1, tolerance = 1e-10)
  }
})

# ---------------------------------------------------------------------------
# Alternative parameterisations -- ROU(parameterization = ).
#
# These are restricted charts on (v, k, s), all reaching the same solver and
# cache.  The organising test is therefore "is it the rate likelihood at the
# mapped point?"; each chart's defining map is checked directly. See
# src/fpe_race.h, rou_map_to_rate().
# ---------------------------------------------------------------------------

rou_map <- function(par_kind, p1, p2, p3, B, A = 0) {
  n <- max(length(p1), length(p2), length(p3), length(B), length(A))
  rep_n <- function(x) rep_len(x, n)
  EMC2:::rou_to_rate_vec(par_kind, rep_n(p1), rep_n(p2), rep_n(p3),
                         rep_n(B), rep_n(A))
}

test_that("the curvature map fixes the reference crossing time", {
  b <- 1.4; tstar <- 0.6; ss <- 0.35
  ks <- c(0, 1e-8, 0.05, 0.5, 1, 3, 5)
  m <- rou_map(1L, tstar, ks, ss, b)

  EX <- ifelse(m$k > 0, (m$v / m$k) * -expm1(-m$k * tstar), m$v * tstar)
  expect_equal(EX, rep(b, length(ks)), tolerance = 1e-12)
  expect_equal(m$k, ks, tolerance = 1e-14)
  expect_equal(m$s, rep(ss, length(ks)), tolerance = 1e-14)
  expect_identical(m$k[1], 0)
  expect_equal(m$v[1], b / tstar, tolerance = 1e-14)
})

test_that("the equilibrium map is centered on the mean start", {
  b <- 1.4; A <- 0.6; d <- b + A / 2; tk <- 0.5; chi <- 0.3
  qs <- c(0.3, 0.8, 1, 1.6)
  m <- rou_map(2L, tk, qs, chi, B = b, A = A)
  expect_equal((m$v / m$k - A / 2) / d, qs, tolerance = 1e-14)
  expect_equal(m$k, rep(1 / tk, length(qs)), tolerance = 1e-14)
  expect_equal(m$s * sqrt(tk) / d, rep(chi, length(qs)), tolerance = 1e-14)
})

test_that("the rate parameterisation is an exact pass-through", {
  m <- rou_map(0L, c(2, 3), c(0.5, 1.5), c(1, 2), c(1, 4), c(0.1, 0.2))
  expect_identical(m$v, c(2, 3))
  expect_identical(m$k, c(0.5, 1.5))
  expect_identical(m$s, c(1, 2))
  # ... and the rate model itself is untouched by the new argument.
  expect_identical(ROU()$c_name, "ROU")
  expect_identical(names(ROU()$p_types),
                   c("v", "k", "B", "A", "t0", "s", "pContaminant"))
  a <- ROU()
  b <- ROU(parameterization = "rate")
  expect_identical(a$c_name, b$c_name)
  expect_identical(a$p_types, b$p_types)
  expect_identical(a$transform, b$transform)
  expect_identical(a$bound, b$bound)
})

test_that("the curvature map uses b = B + A and equilibrium uses B + A/2", {
  m1 <- rou_map(1L, 0.6, 0.8, 0.35, B = 1.0, A = 0.4)
  m2 <- rou_map(1L, 0.6, 0.8, 0.35, B = 1.4, A = 0.0)
  expect_equal(unlist(m1), unlist(m2), tolerance = 1e-14)
  me <- rou_map(2L, 0.5, 0.8, 0.35, B = 1.0, A = 0.4)
  expect_equal((me$v / me$k - 0.2) / 1.2, 0.8, tolerance = 1e-14)
})

test_that("unusable rows fall out of the map rather than being substituted", {
  bad <- rou_map(1L, c(0, -1, NA, 0.6), 0.8, 0.35, 1.4)
  expect_true(all(is.na(bad$v[1:3])))       # tstar <= 0, tstar < 0, NA
  expect_true(is.finite(bad$v[4]))
  expect_true(is.na(rou_map(1L, 0.6, -0.1, 0.35, 1.4)$v))   # k < 0
  expect_true(is.na(rou_map(2L, 0, 0.8, 0.3, 1.4)$v))       # tk <= 0
})

test_that("the alternative densities are the rate density at the mapped point", {
  rt <- seq(0.25, 2.5, by = 0.05)
  rep_rows <- function(p) p[rep(1, length(rt)), , drop = FALSE]

  as_rate <- function(m, B, A, t0) {
    cbind(v = m$v, k = m$k, B = B, A = A, t0 = t0, s = m$s)
  }

  pc <- rep_rows(cbind(tstar = 0.6, k = 0.8, s = 0.35, B = 1.2, A = 0.3,
                       t0 = 0.2))
  mc <- rou_map(1L, pc[, "tstar"], pc[, "k"], pc[, "s"], pc[, "B"], pc[, "A"])
  pr <- as_rate(mc, pc[, "B"], pc[, "A"], pc[, "t0"])
  expect_identical(EMC2:::dROU(rt, pc, par = "curvature"), EMC2:::dROU(rt, pr))
  expect_identical(EMC2:::pROU(rt, pc, par = "curvature"), EMC2:::pROU(rt, pr))

  pe <- rep_rows(cbind(tk = 0.5, theta = 0.8, chi = 0.5, B = 1.2, A = 0.3,
                       t0 = 0.2))
  me <- rou_map(2L, pe[, "tk"], pe[, "theta"], pe[, "chi"], pe[, "B"], pe[, "A"])
  pr2 <- as_rate(me, pe[, "B"], pe[, "A"], pe[, "t0"])
  expect_identical(EMC2:::dROU(rt, pe, par = "equilibrium"), EMC2:::dROU(rt, pr2))
  expect_identical(EMC2:::pROU(rt, pe, par = "equilibrium"), EMC2:::pROU(rt, pr2))
})

test_that("theta < 1 leaves the accumulator subthreshold and escaping by noise", {
  # The regime the equilibrium parameterisation exists to reach: the mean
  # relaxes BELOW the bound, so responses are noise-driven escapes and a large
  # survivor mass is still there at a plausible deadline.
  rt <- rep(1.5, 3)
  qs <- c(0.6, 1.0, 1.6)
  cdf <- vapply(qs, function(qq) {
    p <- cbind(tk = 0.4, theta = qq, chi = 0.35, B = 1, A = 0, t0 = 0.2)
    EMC2:::pROU(1.5, p, par = "equilibrium")
  }, numeric(1))
  expect_true(all(diff(cdf) > 0))          # higher equilibrium, more finishers
  expect_lt(cdf[1], 0.5)                   # subthreshold: most trials survive
  expect_gt(cdf[3], 0.9)

  # The late hazard flattens to a constant -- the escape signature, and the
  # thing that identifies theta within a single deadline.
  tt <- seq(2, 6, by = 0.5)
  p <- cbind(tk = 0.4, theta = 0.6, chi = 0.35, B = 1, A = 0,
             t0 = 0)[rep(1, length(tt)), ]
  h <- EMC2:::dROU(tt, p, par = "equilibrium") /
       (1 - EMC2:::pROU(tt, p, par = "equilibrium"))
  expect_true(all(is.finite(h)) && all(h > 0))
  expect_lt(max(abs(diff(h[-1]))) / mean(h[-1]), 0.05)
})

test_that("the alternative simulators are the rate simulator at the same point", {
  lR <- factor(rep(c("left", "right"), 40), levels = c("left", "right"))
  pc <- cbind(tstar = 0.6, k = 0.8, s = 0.35, B = 1.2, A = 0.3,
              t0 = 0.2)[rep(1, length(lR)), ]
  mc <- rou_map(1L, pc[, "tstar"], pc[, "k"], pc[, "s"], pc[, "B"], pc[, "A"])
  pr <- cbind(v = mc$v, k = mc$k, B = pc[, "B"], A = pc[, "A"],
              t0 = pc[, "t0"], s = mc$s)

  for (use_cpp in c(TRUE, FALSE)) {
    withr::local_options(emc2.cpp_rfun = use_cpp)
    set.seed(4242)
    a <- EMC2:::.rfun_ROU(lR, pc, par = "curvature")
    set.seed(4242)
    b <- EMC2:::.rfun_ROU(lR, pr)
    expect_identical(a$R, b$R)
    expect_equal(a$rt, b$rt, tolerance = 1e-12)
  }
})

# ---------------------------------------------------------------------------
# Package level: the SAMPLED likelihood must agree too, which is the test that
# exercises the C++ column layouts, the adapter dispatch and the cache key.
# ---------------------------------------------------------------------------

rou_par_fixture <- function(model, formula, constants) {
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  des <- design(data = dat, model = model,
                matchfun = function(d) d$S == d$lR,
                formula = formula, constants = constants)
  list(dat = dat, des = des,
       emc = make_emc(dat, des, type = "single", n_chains = 1, compress = TRUE,
                      rt_resolution = NULL))
}

test_that("the curvature likelihood equals the rate likelihood at the map", {
  skip_on_cran()
  b <- 1.2; tstar <- 0.6; kk <- 0.8; ss <- 0.35
  m <- rou_map(1L, tstar, kk, ss, b, 0)

  fc <- suppressMessages(rou_par_fixture(
    function() ROU(parameterization = "curvature"),
    list(tstar ~ 1, k ~ 1, s ~ 1, B ~ 1, t0 ~ 1), c(A = log(0))))
  fr <- suppressMessages(rou_par_fixture(
    ROU, list(v ~ 1, k ~ 1, B ~ 1, t0 ~ 1),
    c(A = log(0), s = log(m$s))))

  l_c <- rou_ll(fc, c(tstar = log(tstar), k = log(kk), s = log(ss),
                      B = log(b), t0 = log(0.2)))
  l_r <- rou_ll(fr, c(v = log(m$v), k = log(m$k), B = log(b), t0 = log(0.2)))
  expect_true(is.finite(l_c))
  expect_equal(l_c, l_r, tolerance = 1e-10)
})

test_that("the equilibrium likelihood equals the rate likelihood at the map", {
  skip_on_cran()
  b <- 1.2; tk <- 0.5; qq <- 0.9; chi <- 0.45
  m <- rou_map(2L, tk, qq, chi, b, 0)

  fe <- suppressMessages(rou_par_fixture(
    function() ROU(parameterization = "equilibrium"),
    list(tk ~ 1, theta ~ 1, chi ~ 1, B ~ 1, t0 ~ 1), c(A = log(0))))
  fr <- suppressMessages(rou_par_fixture(
    ROU, list(v ~ 1, k ~ 1, B ~ 1, t0 ~ 1),
    c(A = log(0), s = log(m$s))))

  l_e <- rou_ll(fe, c(tk = log(tk), theta = log(qq), chi = log(chi),
                      B = log(b), t0 = log(0.2)))
  l_r <- rou_ll(fr, c(v = log(m$v), k = log(m$k), B = log(b), t0 = log(0.2)))
  expect_true(is.finite(l_e))
  expect_equal(l_e, l_r, tolerance = 1e-10)
})

test_that("k = 0 through the sampled likelihood is the Wiener race", {
  skip_on_cran()
  b <- 1.2; tstar <- 0.6; ss <- b / sqrt(tstar)
  fc <- suppressMessages(rou_par_fixture(
    function() ROU(parameterization = "curvature"),
    list(tstar ~ 1, s ~ 1, B ~ 1, t0 ~ 1),
    c(A = log(0), k = log(0))))
  fr <- suppressMessages(rou_par_fixture(
    ROU, list(v ~ 1, B ~ 1, t0 ~ 1),
    c(A = log(0), k = log(0), s = log(ss))))

  l_c <- rou_ll(fc, c(tstar = log(tstar), s = log(ss), B = log(b),
                      t0 = log(0.2)))
  l_r <- rou_ll(fr, c(v = log(b / tstar), B = log(b), t0 = log(0.2)))
  expect_equal(l_c, l_r, tolerance = 1e-10)
})

test_that("a subthreshold equilibrium design fits a censored data set", {
  # The case Option B is for: an UC deadline, omissions coded rt = Inf with
  # R = NA, and a theta that the sampler is free to move across 1.
  skip_on_cran()
  set.seed(20260801)
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  deadline <- 1.0
  des <- design(data = dat, model = function() ROU(parameterization = "equilibrium"),
                matchfun = matchfun,
                formula = list(tk ~ 1, theta ~ lM, chi ~ 1, B ~ 1, t0 ~ 1),
                contrasts = list(theta = list(lM = matrix(c(-1/2, 1/2), ncol = 1,
                                                      dimnames = list(NULL, "d")))),
                constants = c(A = log(0)))
  sim <- suppressMessages(make_data(
    c(tk = log(0.4), theta = log(0.9), theta_lMd = 0.5, chi = log(0.4), B = log(1),
      t0 = log(0.2)), design = des, n_trials = 40))

  # Impose the deadline: past UC we know only that no response had happened.
  sim$UC <- deadline
  sim$UT <- Inf
  late <- !is.finite(sim$rt) | sim$rt > deadline
  sim$rt[late] <- Inf
  sim$R[late] <- NA
  expect_gt(mean(late), 0.05)               # the regime must actually censor

  emc <- suppressMessages(make_emc(sim, des, type = "single", n_chains = 1,
                                   compress = TRUE, rt_resolution = NULL))
  fx <- list(emc = emc)
  ll <- rou_ll(fx, c(tk = log(0.4), theta = log(0.9), theta_lMd = 0.5, chi = log(0.4),
                     B = log(1), t0 = log(0.2)))
  expect_true(is.finite(ll))

  # A theta far from the generating one must be worse -- the deadline plus the
  # omission rate has to carry information about where the process settles.
  worse <- rou_ll(fx, c(tk = log(0.4), theta = log(0.4), theta_lMd = 0.5,
                        chi = log(0.4), B = log(1), t0 = log(0.2)))
  expect_lt(worse, ll)
})
