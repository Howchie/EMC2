# BAwD fading-clearance (fixed gamma in {0, 1/2, 2/3, 3/4, 1}) -- integration tests.
#
# This file covers approved plan tests 8-10: the R and C++ simulator paths,
# constructor/adapter wiring (c_name -> ContextForRaceModels -> compiled
# likelihood), and the make_data -> calc_ll_oo round trip. Kernel-level
# accuracy (density/CDF closed forms, the frozen-mass integral, the negative-
# moment survivor primitive) is covered elsewhere; tests/testthat/test-bawd.R is
# left untouched.
#
# The reference likelihood below is deliberately NOT a port of test-bawd.R's
# closed-form H()/Jf() algebra (which would risk silently reproducing the
# same gamma-generalisation mistakes as the kernel). Instead it is built from
# the model definition directly: for a fixed launch strength V, the
# accumulator's trajectory X(u) = z + V q(u) - ell c_gamma(u) rises to a
# single peak and then falls, so "hit by time u" is exactly
# "X(min(u, u_peak(V))) >= b" -- a single deterministic threshold on z. The
# CDF is then a single numerical integral over the launch-strength
# distribution (no change of variables, no H-function). This was checked
# against test-bawd.R's trusted gamma = 0 closed forms during development
# (max abs difference ~3e-11 over the same parameter sets test-bawd.R uses)
# and is used here unmodified at gamma > 0, where the same threshold argument
# holds verbatim. The density is the numerical derivative of that CDF.

# ---------------------------------------------------------------------------
# Independent reference: CDF via a single integral over the launch strength
# ---------------------------------------------------------------------------

gi_ref_q <- function(u, k) if (k < 1e-12) u else -expm1(-k * u) / k

gi_ref_cgamma <- function(u, k, gamma) {
  if (k < 1e-12) return(u)
  if (gamma <= 1e-12) return(u)
  if (gamma >= 1 - 1e-12) return(gi_ref_q(u, k))
  -expm1(-gamma * k * u) / (gamma * k)
}

# Time at which Xdot(u) = V exp(-ku) - ell exp(-gamma k u) changes sign; Inf if
# it never does (V <= ell: peak at u = 0; the co-decay gamma = 1 case with
# V > ell: peak at u = Inf).
gi_ref_upeak <- function(V, ell, k, gamma) {
  if (!isTRUE(V > ell)) return(0)
  if (k <= 1e-12) return(Inf)
  if (gamma >= 1 - 1e-12) return(Inf)
  log(V / ell) / ((1 - gamma) * k)
}

gi_ref_X <- function(u, V, k, ell, gamma) V * gi_ref_q(u, k) - ell * gi_ref_cgamma(u, k, gamma)

# Required start point z such that the trajectory is exactly tangent to b by
# time u (given the trajectory rises to a single peak and then falls, "hit by
# u" for a static start z is z >= this threshold).
gi_ref_zeff <- function(u, V, b, k, ell, gamma) {
  ueff <- min(u, gi_ref_upeak(V, ell, k, gamma))
  b - gi_ref_X(ueff, V, k, ell, gamma)
}
gi_ref_zeff_vec <- function(u, V, b, k, ell, gamma)
  vapply(V, function(vv) gi_ref_zeff(u, vv, b, k, ell, gamma), numeric(1))

gi_ref_F_normal <- function(u, v, sv, b, A, k, ell, gamma, posdrift = TRUE) {
  denom <- if (posdrift) pnorm(v / sv) else 1
  lower_v <- if (posdrift) 0 else -Inf
  if (A <= 0) {
    Vc <- tryCatch(
      uniroot(function(V) gi_ref_zeff(u, V, b, k, ell, gamma),
              interval = c(max(ell, 0) + 1e-9, 1e10),
              extendInt = "downX", tol = 1e-13)$root,
      error = function(e) NA)
    if (is.na(Vc)) return(0)
    return((1 - pnorm(Vc, v, sv)) / denom)
  }
  f <- function(V) dnorm(V, v, sv) *
    pmin(pmax((A - gi_ref_zeff_vec(u, V, b, k, ell, gamma)) / A, 0), 1)
  integrate(f, lower_v, Inf, rel.tol = 1e-8)$value / denom
}

gi_ref_F_logn <- function(u, mu, sigma, b, A, k, ell, gamma) {
  if (A <= 0) {
    Vc <- tryCatch(
      uniroot(function(V) gi_ref_zeff(u, V, b, k, ell, gamma),
              interval = c(max(ell, 0) + 1e-9, 1e10),
              extendInt = "downX", tol = 1e-13)$root,
      error = function(e) NA)
    if (is.na(Vc)) return(0)
    return(plnorm(Vc, mu, sigma, lower.tail = FALSE))
  }
  f <- function(V) dlnorm(V, mu, sigma) *
    pmin(pmax((A - gi_ref_zeff_vec(u, V, b, k, ell, gamma)) / A, 0), 1)
  integrate(f, 0, Inf, rel.tol = 1e-8)$value
}

# Central-difference derivative of the CDF, Richardson-extrapolated. The
# integral defining F() is only accurate to ~1e-9 in absolute terms (the
# integrand has a kink where a start point crosses into/out of saturation),
# so an overly small h amplifies that noise; h on the order of 1e-3 with
# Richardson extrapolation was checked against the compiled kernel at
# gamma = 0 during development and matches to 7+ significant digits.
gi_fd <- function(Ffun, u, h = 1e-3) {
  h <- min(h, u / 4)
  d1 <- (Ffun(u + h) - Ffun(u - h)) / (2 * h)
  d2 <- (Ffun(u + h / 2) - Ffun(u - h / 2)) / h
  (4 * d2 - d1) / 3
}
gi_ref_f_normal <- function(u, v, sv, b, A, k, ell, gamma, posdrift = TRUE)
  gi_fd(function(uu) gi_ref_F_normal(uu, v, sv, b, A, k, ell, gamma, posdrift), u)
gi_ref_f_logn <- function(u, mu, sigma, b, A, k, ell, gamma)
  gi_fd(function(uu) gi_ref_F_logn(uu, mu, sigma, b, A, k, ell, gamma), u)

# Independent first-passage oracle (root-finding, not Newton) for the
# per-draw simulator check: root of X(u) = d on the rising limb.
gi_ref_fp <- function(V, z, b, k, ell, gamma) {
  d <- b - z
  if (k < 1e-12) return(if (V > ell) d / (V - ell) else Inf)
  if (ell <= 1e-12) {
    x <- 1 - k * d / V
    return(if (x > 0) -log(x) / k else Inf)
  }
  if (gamma >= 1 - 1e-12) {
    if (!isTRUE(V > ell)) return(Inf)
    x <- 1 - k * d / (V - ell)
    return(if (x > 0) -log(x) / k else Inf)
  }
  if (!isTRUE(V > ell)) return(Inf)
  up <- gi_ref_upeak(V, ell, k, gamma)
  if (gi_ref_X(up, V, k, ell, gamma) < d) return(Inf)
  uniroot(function(u) gi_ref_X(u, V, k, ell, gamma) - d, c(0, up), tol = 1e-13)$root
}

cpp_p_g <- function(u, p1, p2, b, A, k, ell, launch, gamma, posdrift = TRUE) {
  EMC2:::pbawd(t = u, A = A, b = b, p1 = p1, p2 = p2, k = k, ell = ell,
               launch = as.integer(launch), posdrift = posdrift, gamma = gamma)
}

# ---------------------------------------------------------------------------
# Constructor / c_name wiring
# ---------------------------------------------------------------------------

test_that("BAwD() encodes fixed gamma regimes into c_name", {
  expect_equal(BAwD(gamma = 0.5)$c_name, "BAwD_LOGN_GAM12")
  expect_equal(BAwD("normal", gamma = 1)$c_name, "BAwD_GAM100")
  expect_equal(BAwD("lognormal", gamma = 2 / 3)$c_name,
               "BAwD_LOGN_GAM23")
  expect_equal(BAwD("lognormal", gamma = 0.75)$c_name,
               "BAwD_LOGN_GAM34")
  expect_equal(BAwD("normal")$c_name, "BAwD")
  expect_equal(BAwD("normal", posdrift = FALSE, gamma = 0.75)$c_name,
               "BAwDIO_GAM34")
  expect_error(BAwD(gamma = 0.25), "one of 0, 1/2, 2/3, 3/4, 1")
  expect_error(BAwD(gamma = c(0, 0.5)), "one of 0, 1/2, 2/3, 3/4, 1")
  expect_false("gamma" %in% names(BAwD(gamma = 0.5)$p_types))
  expect_false("gamma" %in% names(BAwD(gamma = 2 / 3)$p_types))
  expect_equal(BAwD(gamma = 0.5)$gamma, 0.5)
  expect_equal(BAwD(gamma = 2 / 3)$gamma, 2 / 3)
})

# ---------------------------------------------------------------------------
# Simulation: both rfun paths against the compiled CDF, and per-draw hit
# times against an independent root-finding oracle
# ---------------------------------------------------------------------------

test_that("the R and C++ simulators agree distributionally with pbawd across gamma", {
  skip_on_cran()
  lR <- factor(rep(c("left", "right"), 6000), levels = c("left", "right"))
  for (gamma in c(0.5, 2 / 3, 0.75, 1)) {
    for (launch in c(0L, 1L)) {
      nm <- if (launch == 1L) c("mu", "sigma") else c("v", "sv")
      pars <- cbind(0.7, 0.6, 0.8, 0.3, 0.1, 1.2, 0.5)
      colnames(pars) <- c(nm, "B", "A", "t0", "k", "ell")
      pars <- pars[rep(1, length(lR)), ]
      pars <- cbind(pars, b = pars[, "B"] + pars[, "A"])

      set.seed(11)
      a <- EMC2:::rBAwD(lR, pars, launch = launch, gamma = gamma)
      set.seed(11)
      b <- withr::with_options(list(emc2.cpp_rfun = TRUE),
        EMC2:::.rfun_BAwD(lR, pars, launch = launch, gamma = gamma))

      for (dat in list(a, b)) {
        fin <- is.finite(dat$rt)
        probe <- c(0.2, 0.35, 0.5, 0.8)
        th <- vapply(probe, function(x) {
          Fx <- cpp_p_g(x - 0.1, 0.7, 0.6, 1.1, 0.3, 1.2, 0.5, launch, gamma)
          1 - (1 - Fx)^2
        }, numeric(1))
        emp <- vapply(probe, function(x) mean(fin & dat$rt <= x), numeric(1))
        expect_lt(max(abs(emp - th)), 2.5e-2)
        expect_true(all(is.na(dat$R[!is.finite(dat$rt)])))
      }
      # gamma < 1 has a finite T_max, so omissions must occur; gamma = 1 has
      # none, but the drive asymptote (V - ell)/k can still fall short.
      expect_true(any(!is.finite(a$rt)))
    }
  }
})

test_that(".bawd_hit_time matches an independent root-finding oracle", {
  skip_on_cran()
  set.seed(99)
  n <- 200
  k <- 1.1; ell <- 0.7; b <- 1.3
  for (gamma in c(0, 0.5, 2 / 3, 0.75, 1)) {
    V <- rlnorm(n, 0.4, 0.6)
    z <- runif(n, 0, 0.4)
    got <- mapply(EMC2:::.bawd_hit_time, V, b - z, k, ell, MoreArgs = list(gamma = gamma))
    ref <- mapply(gi_ref_fp, V, z, MoreArgs = list(b = b, k = k, ell = ell, gamma = gamma))
    finite_both <- is.finite(got) & is.finite(ref)
    expect_gt(sum(finite_both), 10L)
    expect_equal(got[finite_both], ref[finite_both], tolerance = 1e-9)
    expect_equal(is.finite(got), is.finite(ref))
  }
})

# ---------------------------------------------------------------------------
# Compiled adapter: c_name -> ContextForRaceModels -> race likelihood
# ---------------------------------------------------------------------------

ref_race_ll_gamma <- function(dadm, pars, launch, gamma, min_ll = log(1e-10)) {
  nm <- if (launch == 1L) c("mu", "sigma") else c("v", "sv")
  Fu <- function(i, u) {
    if (!isTRUE(u > 0)) return(0)
    if (launch == 1L)
      gi_ref_F_logn(u, pars[i, nm[1]], pars[i, nm[2]], pars[i, "b"], pars[i, "A"],
                    pars[i, "k"], pars[i, "ell"], gamma)
    else
      gi_ref_F_normal(u, pars[i, nm[1]], pars[i, nm[2]], pars[i, "b"], pars[i, "A"],
                      pars[i, "k"], pars[i, "ell"], gamma)
  }
  fu <- function(i, u) {
    if (!isTRUE(u > 0)) return(0)
    if (launch == 1L)
      gi_ref_f_logn(u, pars[i, nm[1]], pars[i, nm[2]], pars[i, "b"], pars[i, "A"],
                    pars[i, "k"], pars[i, "ell"], gamma)
    else
      gi_ref_f_normal(u, pars[i, nm[1]], pars[i, nm[2]], pars[i, "b"], pars[i, "A"],
                      pars[i, "k"], pars[i, "ell"], gamma)
  }
  n_lR <- length(levels(dadm$lR))
  n_tr <- nrow(dadm) / n_lR
  lt <- numeric(n_tr)
  for (j in seq_len(n_tr)) {
    idx <- ((j - 1) * n_lR + 1):(j * n_lR)
    rt <- dadm$rt[idx[1]]
    w <- idx[which(dadm$winner[idx])]
    term <- fu(w, rt - pars[w, "t0"])
    for (i in setdiff(idx, w))
      term <- term * (1 - Fu(i, rt - pars[i, "t0"]))
    lt[j] <- max(log(max(term, 0)), min_ll)
  }
  sum(lt[attr(dadm, "expand")])
}

# Direct kernel assembly complements the independent numerical oracle above:
# it isolates c_name -> adapter routing from finite-difference error at seams.
ref_race_ll_kernel <- function(dadm, pars, launch, gamma,
                               min_ll = log(1e-10)) {
  nm <- if (launch == 1L) c("mu", "sigma") else c("v", "sv")
  Fu <- function(i, u) {
    if (!isTRUE(u > 0)) return(0)
    EMC2:::pbawd(t = u, A = pars[i, "A"], b = pars[i, "b"],
                 p1 = pars[i, nm[1]], p2 = pars[i, nm[2]],
                 k = pars[i, "k"], ell = pars[i, "ell"],
                 launch = launch, gamma = gamma)
  }
  fu <- function(i, u) {
    if (!isTRUE(u > 0)) return(0)
    EMC2:::dbawd(t = u, A = pars[i, "A"], b = pars[i, "b"],
                 p1 = pars[i, nm[1]], p2 = pars[i, nm[2]],
                 k = pars[i, "k"], ell = pars[i, "ell"],
                 launch = launch, gamma = gamma)
  }
  n_lR <- length(levels(dadm$lR))
  n_tr <- nrow(dadm) / n_lR
  lt <- numeric(n_tr)
  for (j in seq_len(n_tr)) {
    idx <- ((j - 1) * n_lR + 1):(j * n_lR)
    rt <- dadm$rt[idx[1]]
    w <- idx[which(dadm$winner[idx])]
    term <- fu(w, rt - pars[w, "t0"])
    for (i in setdiff(idx, w))
      term <- term * (1 - Fu(i, rt - pars[i, "t0"]))
    lt[j] <- max(log(max(term, 0)), min_ll)
  }
  sum(lt[attr(dadm, "expand")])
}

bawd_gamma_ll_fixture <- function(model, form, consts, dat_n = 24) {
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  dat <- dat[seq_len(dat_n), ]
  des <- suppressMessages(design(data = dat, model = model, matchfun = matchfun,
                                 formula = form, constants = consts))
  emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1,
                                   compress = FALSE, rt_resolution = NULL))
  list(dat = dat, des = des, emc = emc)
}

bawd_gamma_ll <- function(fx, p) {
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

test_that("the compiled adapter routes fixed gamma regimes correctly", {
  skip_on_cran()

  # Lognormal launch, gamma = 1/2.
  fx_ln <- bawd_gamma_ll_fixture(function() BAwD(gamma = 0.5),
    list(mu ~ 1, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1, ell ~ 1), NULL)
  p_ln <- c(mu = 0.9, sigma = log(0.6), B = log(0.8), A = log(0.3),
            t0 = log(0.15), k = log(0.8),
            ell = log(0.5))[names(sampled_pars(fx_ln$des))]
  ll_gam <- bawd_gamma_ll(fx_ln, p_ln)
  expect_true(is.finite(ll_gam))

  dadm_ln <- fx_ln$emc[[1]]$data[[1]]
  pars_ln <- EMC2:::get_pars_matrix_oo(p_ln, dadm_ln, fx_ln$emc[[1]]$model())
  ref_ll <- ref_race_ll_gamma(dadm_ln, pars_ln, launch = 1L, gamma = 0.5)
  expect_equal(ll_gam, ref_ll, tolerance = 1e-4)

  # Same design/data/parameters at gamma = 0 must route to a DIFFERENT
  # likelihood -- this is what pins the c_name -> ctx.bawd_gamma hop; equal
  # values would mean the suffix never reached the adapter.
  fx0 <- bawd_gamma_ll_fixture(BAwD,
    list(mu ~ 1, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1, ell ~ 1), NULL)
  ll0 <- bawd_gamma_ll(fx0, p_ln)
  expect_true(is.finite(ll0))
  expect_gt(abs(ll_gam - ll0), 1e-3)
  expect_equal(fx0$emc[[1]]$model()$c_name, "BAwD_LOGN")

  # Lognormal launch, gamma = 2/3 (the newly added interior regime).
  fx_23 <- bawd_gamma_ll_fixture(function() BAwD(gamma = 2 / 3),
    list(mu ~ 1, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1, ell ~ 1), NULL)
  p_23 <- p_ln
  ll_23 <- bawd_gamma_ll(fx_23, p_23)
  expect_true(is.finite(ll_23))
  dadm_23 <- fx_23$emc[[1]]$data[[1]]
  pars_23 <- EMC2:::get_pars_matrix_oo(p_23, dadm_23, fx_23$emc[[1]]$model())
  ref_23 <- ref_race_ll_gamma(dadm_23, pars_23, launch = 1L, gamma = 2 / 3)
  expect_equal(ll_23, ref_23, tolerance = 1e-4)
  expect_gt(abs(ll_23 - ll_gam), 1e-3)
  expect_equal(fx_23$emc[[1]]$model()$c_name, "BAwD_LOGN_GAM23")

  # Normal launch, gamma = 3/4.
  fx_no <- bawd_gamma_ll_fixture(
    function() BAwD("normal", gamma = 0.75),
    list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1, ell ~ 1), c(sv = log(1)))
  # Use the ell = 0.5 regime this test was written for; a shorter window
  # floors most trials at min_ll and the loose numeric reference below then
  # has nothing left to compare against.
  p_no <- c(v = 3, B = log(0.8), A = log(0.3), t0 = log(0.15), k = log(0.8),
            ell = log(0.5))[names(sampled_pars(fx_no$des))]
  ll_no <- bawd_gamma_ll(fx_no, p_no)
  expect_true(is.finite(ll_no))
  dadm_no <- fx_no$emc[[1]]$data[[1]]
  pars_no <- EMC2:::get_pars_matrix_oo(p_no, dadm_no, fx_no$emc[[1]]$model())
  ref_no <- ref_race_ll_kernel(dadm_no, pars_no, launch = 0L, gamma = 0.75)
  expect_equal(ll_no, ref_no, tolerance = 1e-8)
  numeric_no <- ref_race_ll_gamma(dadm_no, pars_no, launch = 0L, gamma = 0.75)
  expect_equal(ll_no, numeric_no, tolerance = 1e-2)
  expect_equal(fx_no$emc[[1]]$model()$c_name, "BAwD_GAM34")
})

# ---------------------------------------------------------------------------
# make_data round trip
# ---------------------------------------------------------------------------

test_that("make_data with gamma one-half produces omissions and finite fits", {
  skip_on_cran()
  set.seed(20260818)
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  ADmat <- matrix(c(-1 / 2, 1 / 2), ncol = 1, dimnames = list(NULL, "d"))
  des <- suppressMessages(design(
    data = dat, model = function() BAwD(gamma = 0.5), matchfun = matchfun,
    formula = list(mu ~ lM, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1,
                   ell ~ 1),
    contrasts = list(mu = list(lM = ADmat))))
  p <- c(mu = 1.2, mu_lMd = 0.8, sigma = log(0.5), B = log(0.7),
         A = log(0.3), t0 = log(0.15), k = log(0.7),
         ell = log(1.0))[names(sampled_pars(des))]

  sim <- make_data(p, design = des, n_trials = 60)
  expect_true(any(is.infinite(sim$rt)))
  expect_true(all(is.na(sim$R[is.infinite(sim$rt)])))
  fin <- is.finite(sim$rt)
  expect_true(all(sim$rt[fin] > 0.15))
  expect_gt(mean(sim$S[fin] == sim$R[fin]), 0.55)

  e <- suppressMessages(make_emc(sim, des, type = "single", n_chains = 1))
  model <- e[[1]]$model()
  expect_equal(model$c_name, "BAwD_LOGN_GAM12")
  dadm <- e[[1]]$data[[1]]
  designs <- list()
  for (nm in names(model$p_types)) {
    designs[[nm]] <- attr(dadm, "designs")[[nm]][
      attr(attr(dadm, "designs")[[nm]], "expand"), , drop = FALSE]
  }
  ll <- EMC2:::calc_ll_oo(matrix(p, nrow = 1, dimnames = list(NULL, names(p))), dadm,
                          constants = attr(dadm, "constants"), designs = designs,
                          type = model$c_name, bounds = model$bound,
                          transforms = model$transform,
                          pretransforms = model$pre_transform,
                          p_types = names(model$p_types),
                          min_ll = log(1e-10), trend = model$trend)
  expect_true(is.finite(ll))

  # Co-decay: still finite likelihood and still capable of omissions, but no
  # finite T_max.
  expect_equal(EMC2:::bawd_tmax(0.3, 1.1, 0.8, 1, 1), Inf)
  des1 <- suppressMessages(design(
    data = dat, model = function() BAwD(gamma = 1), matchfun = matchfun,
    formula = list(mu ~ lM, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1),
    contrasts = list(mu = list(lM = ADmat)), constants = c(ell = log(1))))
  p1 <- p[names(sampled_pars(des1))]
  sim1 <- make_data(p1, design = des1, n_trials = 60)
  expect_true(any(is.infinite(sim1$rt)))
  e1 <- suppressMessages(make_emc(sim1, des1, type = "single", n_chains = 1))
  m1 <- e1[[1]]$model()
  expect_equal(m1$c_name, "BAwD_LOGN_GAM100")
  dadm1 <- e1[[1]]$data[[1]]
  designs1 <- list()
  for (nm in names(m1$p_types)) {
    designs1[[nm]] <- attr(dadm1, "designs")[[nm]][
      attr(attr(dadm1, "designs")[[nm]], "expand"), , drop = FALSE]
  }
  ll1 <- EMC2:::calc_ll_oo(matrix(p1, nrow = 1, dimnames = list(NULL, names(p1))), dadm1,
                           constants = attr(dadm1, "constants"), designs = designs1,
                           type = m1$c_name, bounds = m1$bound,
                           transforms = m1$transform,
                           pretransforms = m1$pre_transform,
                           p_types = names(m1$p_types),
                           min_ll = log(1e-10), trend = m1$trend)
  expect_true(is.finite(ll1))
})
