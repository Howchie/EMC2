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

test_that("the leak slows the race and makes the upper tail defective", {
  # theta = v/k is the asymptote.  With v = 1.5 and b = 1, k = 3 puts the
  # asymptote at 0.5, BELOW threshold, so finishing depends on noise alone.
  cdf2 <- vapply(c(0, 1, 3), function(kk) {
    EMC2:::pROU(2.0, rou_pars(2.0, v = 1.5, k = kk, B = 1, A = 0))
  }, numeric(1))
  expect_true(all(diff(cdf2) < 0))     # more leak, less mass by 2 s
  expect_gt(cdf2[1], 0.95)
  expect_lt(cdf2[3], 0.85)

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
    list(des = des,
         emc = make_emc(dat, des, type = "single", n_chains = 1, compress = TRUE))
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
  dat$rt[dat$rt > cut] <- Inf
  ADmat <- matrix(c(-1 / 2, 1 / 2), ncol = 1, dimnames = list(NULL, "d"))

  mk <- function(model, form, consts) {
    des <- design(data = dat, model = model, matchfun = matchfun, formula = form,
                  contrasts = list(v = list(lM = ADmat)), constants = consts)
    list(des = des, emc = suppressMessages(
      make_emc(dat, des, type = "single", n_chains = 1, compress = TRUE,
               rt_resolution = 0.02)))
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
  for (kind in c("exponential", "linear", "weibull")) {
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
  for (kind in c("exponential", "linear", "weibull")) {
    pw <- if (kind == "weibull") 1.5 else NULL
    coll <- rou_bnd_pars(n, pw = pw)
    cf <- EMC2:::pROU(rt, fixed); cc <- EMC2:::pROU(rt, coll, kind = kind)
    expect_true(all(cc >= cf - 1e-9), info = kind)
    expect_gt(max(cc - cf), 0.05)
    # A leaky accumulator that could never finish against a fixed bound can
    # finish against a collapsing one, so the defect must shrink too.
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
  for (kind in c("exponential", "linear", "weibull")) {
    pw <- if (kind == "weibull") 1.5 else NULL
    ps <- rou_bnd_pars(N, pw = pw)
    bk <- c(exponential = 2L, linear = 3L, weibull = 1L)[[kind]]
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
  for (bk in c(1L, 2L, 3L)) {
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
