# Focused observable contracts for the ballistic operational-time warp.
# The slow simulation block is gated by skip_model_validation(); algebra and
# likelihood checks are intentionally small and run in the ordinary test pass.

tw_models <- list(
  LBA = function() LBA(),
  BAwL = function() BAwL(),
  BAwD = function() BAwD(),
  BAwDp = function() BAwDp(),
  BAwF = function() BAwF(),
  BAwR = function() BAwR(),
  BTAwL_T = function() BTAwLTransient(),
  BTAwL_S = function() BTAwLSustained(),
  BTAwL = function() BTAwL()
)

# Natural-scale matrices for the default launch parameterisation of each
# constructor.  The extra `b` column is what the single-accumulator wrappers
# consume; `B` is retained for Ttransform endpoint checks.
tw_nat_pars <- function(name, eta = 0, n = 1L) {
  vals <- switch(name,
    LBA = c(v = 2, sv = .35, B = 1, A = .2, b = 1.2, t0 = .1),
    BAwL = c(v = 2, sv = .35, B = 1, A = .2, b = 1.2, t0 = .1,
             k = .45, lambda_g = 0, lambda_k = 0),
    BAwD = c(mean = .25, cv = .35, B = 1, A = .2, b = 1.2, t0 = .1,
             k = .55, ell = .75),
    BAwDp = c(mean = .25, cv = .35, B = 1, A = .2, b = 1.2, t0 = .1,
              k = .8, lambda = .35),
    BAwF = c(mean = .25, cv = .35, B = 1, A = .2, b = 1.2, t0 = .1,
             k = .55),
    BAwR = c(mean = .25, cv = .35, B = 1, A = .2, b = 1.2, t0 = .1,
             kappa = .25, p = 1),
    BTAwL_T = c(mu = .25, sigma = .35, B = 1, A = .2, b = 1.2, t0 = .1,
                k = .55, tau = 1.1),
    BTAwL_S = c(mu = .25, sigma = .35, B = 1, A = .2, b = 1.2, t0 = .1,
                k = .55, tau_s = 1.1),
    BTAwL = c(mu = .25, sigma = .35, B = 1, A = .2, b = 1.2, t0 = .1,
              k = .55, tau_s = 1.1, tau_t = 1.1, pi = .5)
  )
  vals <- c(vals, eta = eta)
  out <- matrix(rep(unname(vals), length.out = n * length(vals)),
                nrow = n, byrow = TRUE,
                dimnames = list(NULL, names(vals)))
  out
}

# ---------------------------------------------------------------------------
# Shared warp algebra and guards
# ---------------------------------------------------------------------------

test_that("time-warp helpers are identity, inverse, and limit safe", {
  u <- c(1e-6, .01, .1, .5, 1, 3, 10)
  eta <- c(-3, -1, -.3, 0, .3, 1, 3)
  for (e in eta) {
    expect_equal(EMC2:::.tw_inv(EMC2:::.tw_fwd(u, e), e), u,
                 tolerance = 1e-10, info = paste("eta", e))
  }
  expect_identical(EMC2:::.tw_fwd(u, 0), u)
  expect_identical(EMC2:::.tw_inv(u, 0), u)
  expect_identical(EMC2:::.tw_log_jac(u, 0), numeric(length(u)))

  # Explicit omega -> 0 limits avoid an unstable expm1/log1p quotient.
  expect_equal(EMC2:::.tw_fwd(u, -40), log1p(u), tolerance = 0)
  expect_equal(EMC2:::.tw_inv(u, -40), expm1(u), tolerance = 1e-12)

  expect_true(all(is.infinite(EMC2:::.tw_fwd(Inf, c(-3, 0, 3)))))
  expect_equal(EMC2:::.tw_fwd(1, 1000), Inf)
  expect_equal(EMC2:::.tw_inv(1, 1000), 0)
  expect_equal(EMC2:::.tw_log_jac(1, 1000), Inf)

  guarded <- c(-1, 0, NA_real_, Inf)
  expect_identical(EMC2:::.tw_fwd(guarded, .3), guarded)
  expect_identical(EMC2:::.tw_inv(guarded, .3), guarded)
  expect_identical(EMC2:::.tw_log_jac(guarded, .3), rep(0, length(guarded)))
  # NaN is an invalid active parameter, not the eta = 0 identity.
  expect_true(is.nan(EMC2:::.tw_fwd(1, NaN)))
  expect_true(is.nan(EMC2:::.tw_log_jac(1, NaN)))
  expect_true(is.nan(EMC2:::.tw_inv(1, NaN)))
})


# ---------------------------------------------------------------------------
# R wrapper composition: every supported ballistic constructor has the same
# outer map, with no Jacobian on the CDF and one Jacobian on the density.
# ---------------------------------------------------------------------------

test_that("all ballistic R wrappers compose density and CDF with the warp", {
  u <- c(.18, .35, .65)
  for (nm in names(tw_models)) {
    model <- tw_models[[nm]]()
    for (e in c(-.6, .6)) {
      warped <- tw_nat_pars(nm, e, length(u))
      parent <- tw_nat_pars(nm, 0, length(u))
      rt <- warped[, "t0"] + u
      s <- EMC2:::.tw_fwd(u, e)
      d_ref <- model$dfun(parent[, "t0"] + s, parent) *
        EMC2:::.tw_jac(u, e)
      p_ref <- model$pfun(parent[, "t0"] + s, parent)
      expect_equal(model$dfun(rt, warped), d_ref, tolerance = 2e-8,
                   info = paste(nm, "density", e))
      expect_equal(model$pfun(rt, warped), p_ref, tolerance = 2e-8,
                   info = paste(nm, "CDF", e))

      # Centered differences check the observable CDF/PDF contract without
      # requiring numDeriv (and stay away from hard endpoints).
      h <- pmin(1e-5, u / 10)
      fd <- (model$pfun(rt + h, warped) - model$pfun(rt - h, warped)) / (2 * h)
      expect_equal(fd, model$dfun(rt, warped), tolerance = 2e-4,
                   info = paste(nm, "finite difference", e))
    }
  }
})

# Intrinsic never-finish probability is on the operational endpoint, hence it
# cannot change under a one-to-one finite-time reparameterisation.
test_that("the intrinsic defect is invariant under eta", {
  defective <- c("BAwD", "BAwDp", "BAwF", "BAwR", "BTAwL_T")
  for (nm in defective) {
    model <- tw_models[[nm]]()
    p0 <- tw_nat_pars(nm, 0)
    base <- as.numeric(model$pfun(Inf, p0))
    for (e in c(-2, -.5, .5, 2)) {
      got <- as.numeric(model$pfun(Inf, tw_nat_pars(nm, e)))
      expect_equal(got, base, tolerance = 1e-12, info = paste(nm, e))
    }
  }
})

# ---------------------------------------------------------------------------
# Endpoint reporting and race ordering
# ---------------------------------------------------------------------------

test_that("finite endpoints are inverse-warped and add t0 on the physical clock", {
  cases <- list(
    BAwD = list(model = BAwD(), op = function(p)
      EMC2:::bawd_tmax_vec(p[, "A"], p[, "b"], p[, "k"], p[, "ell"])),
    BAwDp = list(model = BAwDp(), op = function(p)
      -log(p[, "lambda"]) / p[, "k"]),
    BAwF = list(model = BAwF(), op = function(p)
      EMC2:::bawf_tmax_vec(p[, "A"], p[, "b"], p[, "k"], Inf)),
    BAwR = list(model = BAwR(), op = function(p)
      EMC2:::bawr_tmax_vec(p[, "A"], p[, "b"], p[, "kappa"], p[, "p"])),
    BTAwL_T = list(model = BTAwLTransient(), op = function(p)
      EMC2:::btawl_tmax_vec(p[, "k"], p[, "tau"]))
  )
  for (nm in names(cases)) {
    p <- tw_nat_pars(nm, .55)
    endpoint <- cases[[nm]]$model$Ttransform(p, NULL)
    op <- cases[[nm]]$op(p)
    expect_equal(endpoint[, "Tmax"], unname(EMC2:::.tw_inv(op, .55)),
                 tolerance = 1e-12, info = paste(nm, "inverse endpoint"))
    expect_identical(as.numeric(endpoint[, "rt_max"]),
                     as.numeric(endpoint[, "t0"] + endpoint[, "Tmax"]))

    eps <- .01 * max(1, endpoint[, "Tmax"])
    below <- endpoint[, "rt_max"] - eps
    above <- endpoint[, "rt_max"] + eps
    expect_true(is.finite(cases[[nm]]$model$dfun(below, p)) &&
                  cases[[nm]]$model$dfun(below, p) > 0,
                info = paste(nm, "below endpoint"))
    expect_equal(cases[[nm]]$model$dfun(above, p), 0, tolerance = 1e-12,
                 info = paste(nm, "endpoint density"))
    expect_equal(cases[[nm]]$model$pfun(above, p),
                 cases[[nm]]$model$pfun(endpoint[, "rt_max"], p),
                 tolerance = 1e-12, info = paste(nm, "endpoint CDF"))
  }
})

test_that("a common warp preserves the race winner but changes physical time", {
  old <- options(emc2.cpp_rfun = TRUE)
  on.exit(options(old), add = TRUE)
  model <- LBA()
  levels <- c("left", "right")
  lR <- factor(rep(levels, 100), levels = levels)
  p0 <- tw_nat_pars("LBA", 0, length(lR))
  p1 <- tw_nat_pars("LBA", 1, length(lR))
  attr(p0, "ok") <- attr(p1, "ok") <- rep(TRUE, length(lR))
  set.seed(90210)
  x0 <- model$rfun(list(lR = lR), p0)
  set.seed(90210)
  x1 <- model$rfun(list(lR = lR), p1)
  expect_identical(x0$R, x1$R)
  expect_true(any(abs(x0$rt - x1$rt) > 1e-10))
})

# ---------------------------------------------------------------------------
# A compiled adapter comparison against an independent exported-kernel
# reference.  This also exercises name-based eta resolution in calc_ll_oo.
# ---------------------------------------------------------------------------

tw_bawr_context <- function(dat, include_eta = TRUE) {
  form <- list(mean ~ 1, cv ~ 1, B ~ 1, A ~ 1, t0 ~ 1,
               kappa ~ 1, p ~ 1)
  if (include_eta) form <- c(form, list(eta ~ 1))
  des <- design(data = dat, Rlevels = levels(dat$R), model = BAwR,
                matchfun = function(d) as.character(d$S) == as.character(d$lR),
                formula = form, report_p_vector = FALSE)
  emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1,
                                   compress = FALSE, rt_resolution = NULL))
  list(design = des, emc = emc)
}

tw_calc_ll <- function(ctx, p) {
  dadm <- ctx$emc[[1]]$data[[1]]
  model <- ctx$emc[[1]]$model()
  dz <- attr(dadm, "designs")
  designs <- lapply(names(model$p_types), function(nm)
    dz[[nm]][attr(dz[[nm]], "expand"), , drop = FALSE])
  names(designs) <- names(model$p_types)
  constants <- attr(dadm, "constants")
  if (is.null(constants)) constants <- NA
  as.numeric(EMC2:::calc_ll_oo(
    matrix(p, nrow = 1, dimnames = list(NULL, names(p))), dadm,
    constants = constants, designs = designs, type = model$c_name,
    bounds = model$bound, transforms = model$transform,
    pretransforms = model$pre_transform, p_types = names(model$p_types),
    min_ll = log(1e-12), trend = model$trend))
}

tw_ref_bawr_ll <- function(dadm, pars, min_ll = log(1e-12)) {
  nr <- length(levels(dadm$lR))
  ntrial <- nrow(dadm) / nr
  one <- numeric(ntrial)
  for (j in seq_len(ntrial)) {
    idx <- ((j - 1L) * nr + 1L):(j * nr)
    rt <- dadm$rt[idx[1L]]
    if (is.finite(rt)) {
      w <- idx[which(dadm$winner[idx])[1L]]
      u <- rt - pars[w, "t0"]
      s <- EMC2:::.tw_fwd(u, pars[w, "eta"])
      term <- EMC2:::dbawr(s, pars[w, "A"], pars[w, "b"], pars[w, "mean"],
                           pars[w, "cv"], pars[w, "kappa"], pars[w, "p"],
                           launch = 1L) * EMC2:::.tw_jac(u, pars[w, "eta"])
      for (i in setdiff(idx, w)) {
        ui <- rt - pars[i, "t0"]
        si <- EMC2:::.tw_fwd(ui, pars[i, "eta"])
        term <- term * (1 - EMC2:::pbawr(si, pars[i, "A"], pars[i, "b"],
                                         pars[i, "mean"], pars[i, "cv"],
                                         pars[i, "kappa"], pars[i, "p"],
                                         launch = 1L))
      }
    } else {
      term <- prod(vapply(idx, function(i) {
        ui <- Inf - pars[i, "t0"]
        si <- EMC2:::.tw_fwd(ui, pars[i, "eta"])
        1 - EMC2:::pbawr(si, pars[i, "A"], pars[i, "b"], pars[i, "mean"],
                         pars[i, "cv"], pars[i, "kappa"], pars[i, "p"],
                         launch = 1L)
      }, numeric(1)))
    }
    one[j] <- max(log(max(term, 0)), min_ll)
  }
  sum(one[attr(dadm, "expand")])
}

test_that("compiled BAwR likelihood equals its warped exported-kernel reference", {
  dat <- data.frame(
    subjects = factor("s1"), S = factor("stim"),
    R = factor("left", levels = c("left", "right")), rt = .72,
    LT = 0, UT = Inf, LC = 0, UC = Inf
  )
  ctx <- tw_bawr_context(dat, include_eta = TRUE)
  p <- sampled_pars(ctx$design, doMap = FALSE)
  p[c("mean", "cv", "B", "A", "t0", "kappa", "p", "eta")] <-
    c(.25, log(.35), log(1), log(.2), log(.1), log(.25), log(1), .45)
  dadm <- ctx$emc[[1]]$data[[1]]
  pars <- EMC2:::get_pars_matrix_oo(p, dadm, ctx$emc[[1]]$model())
  expect_equal(tw_calc_ll(ctx, p), tw_ref_bawr_ll(dadm, pars), tolerance = 1e-9)
})

# ---------------------------------------------------------------------------
# Parent recovery, censoring, and truncation.  The make_data check deliberately
# compares the same RNG stream, not merely summary statistics.
# ---------------------------------------------------------------------------

test_that("eta = 0 recovers the parent design and simulator exactly", {
  dat <- data.frame(
    subjects = factor(rep("s1", 4)), S = factor(rep("stim", 4)),
    R = factor(rep(c("left", "right"), 2), levels = c("left", "right")),
    rt = c(.5, .6, .7, .8), LT = 0, UT = Inf, LC = 0, UC = Inf
  )
  no_eta <- tw_bawr_context(dat, include_eta = FALSE)
  with_eta <- tw_bawr_context(dat, include_eta = TRUE)
  base <- c(mean = .25, cv = log(.35), B = log(1), A = log(.2),
            t0 = log(.1), kappa = log(.25), p = log(1))
  p0 <- sampled_pars(no_eta$design, doMap = FALSE)
  p1 <- sampled_pars(with_eta$design, doMap = FALSE)
  p0[names(base)] <- base
  p1[names(base)] <- base
  p1["eta"] <- 0
  expect_identical(tw_calc_ll(no_eta, p0), tw_calc_ll(with_eta, p1))
  set.seed(711)
  d0 <- make_data(p0, no_eta$design, n_trials = 20)
  set.seed(711)
  d1 <- make_data(p1, with_eta$design, n_trials = 20)
  expect_identical(d0$R, d1$R)
  expect_identical(d0$rt, d1$rt)
})

test_that("upper censoring uses the warped survivor and truncation normalises it", {
  # Two independent racers: a +Inf response censored at UC has survivor equal
  # to the product of the two parent survivors at c_eta(UC - t0).
  dat_uc <- data.frame(
    subjects = factor("s1"), S = factor("stim"),
    R = factor(NA_character_, levels = c("left", "right")), rt = Inf,
    LT = 0, UT = Inf, LC = 0, UC = .85
  )
  uc <- tw_bawr_context(dat_uc, include_eta = TRUE)
  p <- sampled_pars(uc$design, doMap = FALSE)
  p[c("mean", "cv", "B", "A", "t0", "kappa", "p", "eta")] <-
    c(.25, log(.35), log(1), log(.2), log(.1), log(.25), log(1), -.5)
  dadm <- uc$emc[[1]]$data[[1]]
  pars <- EMC2:::get_pars_matrix_oo(p, dadm, uc$emc[[1]]$model())
  u <- dat_uc$UC - pars[, "t0"]
  s <- EMC2:::.tw_fwd(u, pars[, "eta"])
  ref_uc <- sum(log1p(-vapply(seq_len(nrow(pars)), function(i)
    EMC2:::pbawr(s[i], pars[i, "A"], pars[i, "b"], pars[i, "mean"],
                 pars[i, "cv"], pars[i, "kappa"], pars[i, "p"], launch = 1L),
    numeric(1))))
  expect_equal(tw_calc_ll(uc, p), ref_uc, tolerance = 1e-9)

  # One racer gives a compact independent reference for a finite truncation
  # window: f(rt) / (F(UT) - F(LT)), all on the warped physical clock.
  dat_tr <- data.frame(
    subjects = factor("s1"), S = factor("stim"),
    R = factor("left", levels = "left"), rt = .62,
    LT = .25, UT = .95, LC = 0, UC = Inf
  )
  tr <- tw_bawr_context(dat_tr, include_eta = TRUE)
  p_tr <- sampled_pars(tr$design, doMap = FALSE)
  p_tr[c("mean", "cv", "B", "A", "t0", "kappa", "p", "eta")] <-
    c(.25, log(.35), log(1), log(.2), log(.1), log(.25), log(1), .45)
  pars_tr <- EMC2:::get_pars_matrix_oo(
    p_tr, tr$emc[[1]]$data[[1]], tr$emc[[1]]$model())[1, ]
  tt <- c(dat_tr$rt, dat_tr$LT, dat_tr$UT)
  ss <- EMC2:::.tw_fwd(tt - pars_tr["t0"], pars_tr["eta"])
  f <- EMC2:::dbawr(ss[1], pars_tr["A"], pars_tr["b"], pars_tr["mean"],
                    pars_tr["cv"], pars_tr["kappa"], pars_tr["p"],
                    launch = 1L) * EMC2:::.tw_jac(
                      tt[1] - pars_tr["t0"], pars_tr["eta"])
  z <- EMC2:::pbawr(ss[3], pars_tr["A"], pars_tr["b"], pars_tr["mean"],
                    pars_tr["cv"], pars_tr["kappa"], pars_tr["p"],
                    launch = 1L) -
    EMC2:::pbawr(ss[2], pars_tr["A"], pars_tr["b"], pars_tr["mean"],
                 pars_tr["cv"], pars_tr["kappa"], pars_tr["p"],
                 launch = 1L)
  # BAwR is defective: finite upper truncation retains the intrinsic +Inf
  # atom, so it belongs in the conditional normaliser alongside [LT, UT].
  defect <- 1 - EMC2:::pbawr(Inf, pars_tr["A"], pars_tr["b"], pars_tr["mean"],
                             pars_tr["cv"], pars_tr["kappa"], pars_tr["p"],
                             launch = 1L)
  expect_equal(tw_calc_ll(tr, p_tr), log(f) - log(z + defect),
               tolerance = 1e-09)
})

# Refused variants fail at design construction rather than silently ignoring
# an eta formula.
test_that("clock and correlated BAwL variants refuse eta", {
  dat <- data.frame(subjects = factor("s1"), S = factor("stim"),
                    R = factor("left", levels = c("left", "right")), rt = .7)
  form <- list(v ~ 1, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1, eta ~ 1)
  expect_error(
    design(data = dat, Rlevels = levels(dat$R),
           model = function() BAwL(erlang_type = "local_kill"),
           matchfun = function(d) as.character(d$S) == as.character(d$lR), formula = form),
    "eta"
  )
  expect_error(
    design(data = dat, Rlevels = levels(dat$R),
           model = function() BAwL(correlated = TRUE),
           matchfun = function(d) as.character(d$S) == as.character(d$lR), formula = form),
    "eta"
  )
})

test_that("logical-rule and non-ballistic models refuse eta", {
  dat <- data.frame(subjects = factor("s1"), S = factor("stim"),
                    R = factor("A", levels = c("A", "B")), rt = .7)
  expect_error(
    design(data = dat, Rlevels = levels(dat$R), model = LogicalRulesLBA,
           matchfun = function(d) TRUE, formula = list(eta ~ 1)),
    "eta"
  )
  expect_error(
    design(data = dat, Rlevels = levels(dat$R), model = DDM,
           formula = list(eta ~ 1)),
    "eta"
  )
})

# ---------------------------------------------------------------------------
# Slow Monte-Carlo standard: both the compiled and pure-R simulator should
# sample F_eta(u) = F_0(c_eta(u)).
# ---------------------------------------------------------------------------

test_that("simulators follow the warped accumulation-time CDF", {
  skip_model_validation()
  old <- options(emc2.cpp_rfun = TRUE)
  on.exit(options(old), add = TRUE)
  model <- LBA()
  n <- 2000L
  lR <- factor(rep("one", n), levels = "one")
  grid <- c(.2, .4, .8, 1.2)
  for (cpp in c(TRUE, FALSE)) {
    options(emc2.cpp_rfun = cpp)
    for (e in c(-.7, .7)) {
      p <- tw_nat_pars("LBA", e, n)
      attr(p, "ok") <- rep(TRUE, n)
      set.seed(if (cpp) 418 else 419)
      out <- model$rfun(list(lR = lR), p)
      u <- out$rt - .1
      empirical <- vapply(grid, function(x) mean(u <= x & is.finite(u)), numeric(1))
      reference <- EMC2:::plba(EMC2:::.tw_fwd(grid, e), .2, 1.2, 2, .35,
                               posdrift = TRUE)
      expect_lt(max(abs(empirical - reference)), .07,
                label = paste("cpp", cpp, "eta", e))
    }
  }
})
