# Independent behavioural checks for the finite-rho BAwD family.
# The reference below uses the defining trajectory and root finding; it does not
# duplicate the closed forms used by the compiled implementation.

rho_grid <- c(1, 2, 4, Inf)
gamma_grid <- c(0, .5, 2 / 3, .75, 1)

rho_h <- function(x, rho) if (is.infinite(rho)) exp(-x) else (1 + x / rho)^(-rho)
rho_distance <- function(u, V, z, b, k, ell, gamma, rho) {
  q <- if (is.infinite(rho)) (1 - exp(-k * u)) / k else {
    L <- log1p(k * u / rho)
    rho * (if (rho == 1) L else
      if (abs((1 - rho) * L) < 1e-8) L else
        expm1((1 - rho) * L) / (1 - rho)) / k
  }
  r <- if (gamma == 0) u else if (is.infinite(rho))
    (1 - exp(-gamma * k * u)) / (gamma * k) else {
      L <- log1p(k * u / rho); c <- 1 - rho * gamma
      rho * (if (rho * gamma == 1) L else
        if (abs(c * L) < 1e-8) L else expm1(c * L) / c) / k
    }
  V * q - ell * r - (b - z)
}
rho_vcrit <- function(u, z, b, k, ell, gamma, rho) {
  f <- function(v) rho_distance(u, v, z, b, k, ell, gamma, rho)
  if (f(0) >= 0) return(0)
  hi <- max(ell, 1)
  while (f(hi) < 0 && hi < 1e8) hi <- hi * 2
  if (f(hi) < 0) return(Inf)
  uniroot(f, c(0, hi), tol = 1e-11)$root
}
rho_ref_hit <- function(V, d, k, ell, gamma, rho) {
  if (V <= ell) return(Inf)
  f <- function(u) rho_distance(u, V, 0, d, k, ell, gamma, rho)
  hi <- 1
  while (f(hi) < 0 && hi < 1e7) hi <- hi * 2
  if (f(hi) < 0) return(Inf)
  uniroot(f, c(0, hi), tol = 1e-11)$root
}
rho_ref_cdf <- function(u, p1, p2, b, A, k, ell, launch, gamma, rho) {
  if (u <= 0) return(0)
  zcdf <- function(z) vapply(z, function(zz) {
    vc <- rho_ref_first_vcrit(u, zz, b, k, ell, gamma, rho)
    if (!is.finite(vc)) return(0)
    if (launch == 0L) pnorm(vc, p1, p2, lower.tail = FALSE) / pnorm(p1 / p2)
    else pnorm((log(vc) - p1) / p2, lower.tail = FALSE)
  }, numeric(1))
  if (A <= 0) return(zcdf(0))
  integrate(zcdf, 0, A, rel.tol = 2e-7)$value / A
}
rho_frozen_psi <- function(x, gamma, rho) {
  if (is.infinite(rho)) {
    q <- 1 - exp(-x)
    r <- if (gamma == 0) x else (1 - exp(-gamma * x)) / gamma
    return(exp((1 - gamma) * x) * q - r)
  }
  L <- log1p(x / rho)
  q <- if (rho == 1) L else expm1((1 - rho) * L) / (1 - rho)
  r <- if (gamma == 0) expm1(L) else if (rho * gamma == 1) L else
    expm1((1 - rho * gamma) * L) / (1 - rho * gamma)
  rho * (exp(rho * (1 - gamma) * L) * q - r)
}
rho_frozen_peak <- function(z, b, k, ell, gamma, rho) {
  target <- k * (b - z) / ell
  f <- function(x) rho_frozen_psi(x, gamma, rho) - target
  hi <- 1
  while (f(hi) < 0) hi <- hi * 2
  x <- uniroot(f, c(0, hi), tol = 1e-12)$root
  v <- if (is.infinite(rho)) ell * exp((1 - gamma) * x) else
    ell * exp(rho * (1 - gamma) * log1p(x / rho))
  c(x = x, v = v)
}
rho_frozen_vcrit <- function(z, b, k, ell, gamma, rho)
  rho_frozen_peak(z, b, k, ell, gamma, rho)[["v"]]

rho_ref_first_vcrit <- function(u, z, b, k, ell, gamma, rho) {
  if (gamma < 1) {
    peak <- rho_frozen_peak(z, b, k, ell, gamma, rho)
    if (is.infinite(u) || k * u >= peak[["x"]]) return(peak[["v"]])
  }
  if (is.infinite(u) && gamma == 1 && rho == 1) return(ell)
  rho_vcrit(u, z, b, k, ell, gamma, rho)
}

rho_ref_frozen <- function(p1, p2, b, A, k, ell, launch, gamma, rho) {
  integrand <- function(z) vapply(z, function(zz) {
    v <- rho_frozen_vcrit(zz, b, k, ell, gamma, rho)
    if (launch == 0L) pnorm(v, p1, p2, lower.tail = FALSE) / pnorm(p1 / p2)
    else pnorm((log(v) - p1) / p2, lower.tail = FALSE)
  }, numeric(1))
  integrate(integrand, 0, A, rel.tol = 1e-9)$value / A
}
rho_cpp_p <- function(u, p1, p2, b, A, k, ell, launch, gamma, rho) {
  EMC2:::pbawd(t = u, A = A, b = b, p1 = p1, p2 = p2, k = k, ell = ell,
               launch = launch, gamma = gamma, rho = rho)
}
rho_cpp_d <- function(u, p1, p2, b, A, k, ell, launch, gamma, rho) {
  EMC2:::dbawd(t = u, A = A, b = b, p1 = p1, p2 = p2, k = k, ell = ell,
               launch = launch, gamma = gamma, rho = rho)
}

# Headline identities, including the two values quoted in the model definition.
test_that("rho=2 gamma=0 has the square-root critical corner", {
  expect_equal(EMC2:::bawd_tmax(0, 1, 2, 1, gamma = 0, rho = 2),
               1, tolerance = 2e-12)
  mu <- 1.2; sg <- .9
  p2 <- rho_cpp_p(Inf, mu, sg, 1, 0, 2, 1, 1L, 0, 2)
  expect_equal(p2, plnorm(4, mu, sg, lower.tail = FALSE), tolerance = 2e-9)
  kexp <- 1.25643120862617
  expect_equal(EMC2:::bawd_tmax(0, 1, kexp, 1, gamma = 0, rho = Inf),
               1, tolerance = 2e-12)
  pinf <- rho_cpp_p(Inf, mu, sg, 1, 0, kexp, 1, 1L, 0, Inf)
  expect_equal(pinf, plnorm(3.51286241725235, mu, sg, lower.tail = FALSE),
               tolerance = 2e-9)
  expect_lt(p2, pinf)
})

test_that("finite rho grid is finite, monotone, and density integrates to CDF", {
  for (launch in 0:1) for (rho in rho_grid) for (gamma in gamma_grid) {
    a <- list(p1 = if (launch == 0) 1.4 else .7,
              p2 = if (launch == 0) .8 else .65, b = 1.2, A = .35,
              k = .9, ell = .8, launch = launch, gamma = gamma, rho = rho)
    tt <- c(.08, .3, .8, 1.7, 4)
    pp <- do.call(rho_cpp_p, c(list(u = tt), a))
    dd <- do.call(rho_cpp_d, c(list(u = tt), a))
    expect_true(all(is.finite(pp)) && all(pp >= -1e-10 & pp <= 1 + 1e-10))
    expect_true(all(diff(pp) >= -2e-8))
    expect_true(all(is.finite(dd)) && all(dd >= -1e-9))
    for (u in tt[c(2, 4)]) {
      mass <- integrate(function(x) do.call(rho_cpp_d, c(list(u = x), a)),
                        1e-7, u, rel.tol = 2e-5)$value
      expect_equal(mass, do.call(rho_cpp_p, c(list(u = u), a)), tolerance = 3e-4)
    }
    # An independently root-found trajectory oracle on representative points.
    if (gamma %in% c(0, .5, 1) && rho %in% rho_grid) {
      expect_equal(do.call(rho_cpp_p, c(list(u = .8), a)),
                   do.call(rho_ref_cdf, c(list(u = .8), a)), tolerance = 3e-4)
    }
  }
})

test_that("rho approaching infinity is continuous", {
  for (launch in 0:1) for (gamma in gamma_grid) {
    z <- list(p1 = if (launch == 0) 1.3 else .8, p2 = .7, b = 1.1, A = .2,
              k = .75, ell = .9, launch = launch, gamma = gamma)
    for (u in c(.2, 1, 3)) {
      expect_equal(do.call(rho_cpp_p, c(list(u = u, rho = 1e7), z)),
                   do.call(rho_cpp_p, c(list(u = u, rho = Inf), z)), tolerance = 2e-6)
    }
  }
})

test_that("gamma one corners have the stated limiting survivor", {
  for (rho in c(1, 2, 4)) {
    expect_true(is.infinite(EMC2:::bawd_tmax(.3, 1.1, .8, 1, 1, rho = rho)))
    got <- rho_cpp_p(Inf, .9, .6, 1.1, .3, .8, 1, 1L, 1, rho)
    want <- if (rho == 1) pnorm((log(1) - .9) / .6, lower.tail = FALSE) else
      integrate(function(z) pnorm((log(1 + (rho - 1) / rho * .8 * (1.1 - z)) - .9) / .6,
                                   lower.tail = FALSE), 0, .3)$value / .3
    expect_equal(got, want, tolerance = 2e-5)
    expect_lt(rho_cpp_p(20, .9, .6, 1.1, .3, .8, 1, 1L, 1, rho), want)
  }
})

test_that("negative-m and logratio lognormal primitives agree with integration", {
  for (m in c(-.5, -.25, 0)) for (v in c(.7, 2)) {
    got <- exp(EMC2:::lognormal_power_stoploss_log(v, .2, .8, m))
    want <- integrate(function(w) w^(-(m + 1)) * pnorm((.2 - log(w)) / .8), v, Inf,
                      rel.tol = 2e-9)$value
    expect_equal(got, want, tolerance = 2e-6)
  }
  v <- 1.7
  got <- exp(EMC2:::lognormal_logratio_stoploss_log(v, .2, .8, log(.6)))
  want <- integrate(function(w) log(w / .6) * pnorm((.2 - log(w)) / .8), v, Inf,
                    rel.tol = 2e-9)$value
  expect_equal(got, want, tolerance = 2e-7)
  v_deep <- exp(8); mu_deep <- 0; sigma_deep <- .7; log_ell_deep <- 7.5
  got_deep <- exp(EMC2:::lognormal_logratio_stoploss_log(
    v_deep, mu_deep, sigma_deep, log_ell_deep))
  x_deep <- (log(v_deep) - mu_deep) / sigma_deep
  want_deep <- integrate(function(y)
    (mu_deep + sigma_deep * y - log_ell_deep) *
      pnorm(y, lower.tail = FALSE) * sigma_deep *
      exp(mu_deep + sigma_deep * y),
    x_deep, Inf, rel.tol = 1e-8)$value
  expect_equal(got_deep, want_deep, tolerance = 2e-6)
  expect_true(is.finite(EMC2:::lognormal_power_stoploss_log(exp(8), 0, .7, -.5)))
})

test_that("frozen fallback, boundary hygiene, simulator, and wiring", {
  # A narrow frozen interval exercises the ill-conditioned closed-form fallback.
  for (rho in c(1, 2, 4, Inf)) {
    x <- rho_cpp_p(1.1, .8, .7, 4, 1e-7, 2.5, .9, 1L, .75, rho)
    expect_true(is.finite(x) && x >= 0 && x <= 1)
  }
  # Independent first-passage integration over the frozen start-point mass.
  for (rho in c(1, 2, 4)) {
    got <- rho_cpp_p(Inf, .8, .7, 4, 1e-7, 2.5, .9, 1L, .75, rho)
    want <- rho_ref_frozen(.8, .7, 4, 1e-7, 2.5, .9, 1L, .75, rho)
    expect_equal(got, want, tolerance = 2e-7)
  }
  # The reference must remain first-passage based after the trajectory peaks.
  expect_equal(rho_ref_first_vcrit(Inf, 0, 1, 2, 1, 0, 2), 4,
               tolerance = 2e-12)
  expect_equal(
    rho_cpp_p(2, .8, .7, 1, 0, 2, 1, 1L, 0, 2),
    rho_cpp_p(Inf, .8, .7, 1, 0, 2, 1, 1L, 0, 2),
    tolerance = 2e-12
  )
  for (bad in list(-1, .5, NaN, -Inf))
    expect_error(EMC2::BAwD(gamma = .5, rho = bad), "rho|1, 2, 4|Inf")
  for (rho in c(1, 2, 4, Inf)) for (g in gamma_grid) {
    x <- rho_cpp_p(c(1e-8, 1, 1e5), .8, .7, 1, .2, 1e-7, 1e-8, 1L, g, rho)
    expect_true(all(is.finite(x)) && all(x >= 0 & x <= 1))
  }
  expect_equal(EMC2::BAwD(gamma = .5, rho = 2)$c_name, "BAwD_LOGN_GAM12_RHO2")
  expect_equal(EMC2::BAwD()$c_name, "BAwD_LOGN")
  expect_error(EMC2::BAwD(rho = 3), "rho|1, 2, 4|Inf")
  for (rho in c(1, 2, 4)) for (g in c(0, .5, 1)) {
    for (v in c(1.4, 2.5, 5)) {
      got <- EMC2:::.bawd_hit_time(v, 1.1, .9, 1, gamma = g, rho = rho)
      want <- rho_ref_hit(v, 1.1, .9, 1, g, rho)
      expect_equal(got, want, tolerance = 2e-8)
    }
  }
  for (rho in c(1, 2, 4)) for (g in c(0, .5, 1)) {
    got <- EMC2:::.bawd_hit_time(3, 1.1, .9, 0, gamma = g, rho = rho)
    want <- if (rho == 1) expm1(.9 * 1.1 / 3) / .9 else {
      s <- .9 * 1.1 * (rho - 1) / (3 * rho)
      if (s >= 1) Inf else rho * ((1 - s)^(-1 / (rho - 1)) - 1) / .9
    }
    expect_equal(got, want, tolerance = 2e-12)
  }
  lR_many <- factor(rep(c("L", "R"), 2000), levels = c("L", "R"))
  pars_many <- cbind(mu = rep(.8, length(lR_many)),
                     sigma = rep(.7, length(lR_many)),
                     b = rep(1.1, length(lR_many)),
                     A = rep(.2, length(lR_many)),
                     t0 = rep(.1, length(lR_many)),
                     k = rep(.9, length(lR_many)),
                     ell = rep(1, length(lR_many)))
  set.seed(20260819)
  r_ref <- EMC2:::rBAwD(lR_many, pars_many, launch = 1L, gamma = .5, rho = 2)
  set.seed(20260820)
  r_cpp <- withr::with_options(
    list(emc2.cpp_rfun = TRUE),
    EMC2:::.rfun_BAwD(lR_many, pars_many, launch = 1L, gamma = .5, rho = 2))
  for (dat in list(r_ref, r_cpp)) {
    for (u in c(.3, .7, 1.2)) {
      f <- rho_cpp_p(u - .1, .8, .7, 1.1, .2, .9, 1, 1L, .5, 2)
      expect_lt(abs(mean(is.finite(dat$rt) & dat$rt <= u) - (1 - (1 - f)^2)),
                4e-2)
    }
  }
  expect_true(all(is.na(r_ref$R[!is.finite(r_ref$rt)])))
  pars_k0 <- cbind(mu = rep(2, length(lR_many)),
                   sigma = rep(.01, length(lR_many)),
                   b = rep(.8, length(lR_many)),
                   A = rep(0, length(lR_many)),
                   t0 = rep(.1, length(lR_many)),
                   k = rep(0, length(lR_many)),
                   ell = rep(.5, length(lR_many)))
  set.seed(20260822)
  k0_cpp <- withr::with_options(
    list(emc2.cpp_rfun = TRUE),
    EMC2:::.rfun_BAwD(lR_many, pars_k0, launch = 1L, gamma = .5, rho = 2))
  expect_true(all(is.finite(k0_cpp$rt)))
  pars <- data.frame(mu = rep(.8, 2), sigma = rep(.7, 2), b = rep(1.1, 2),
                     A = rep(.2, 2), t0 = rep(.1, 2), k = rep(.9, 2),
                     ell = rep(1, 2))
  lR <- factor(c("R1", "R2"), levels = c("R1", "R2"))
  set.seed(42); r1 <- EMC2:::rBAwD(lR, pars, launch = 1L, gamma = .5, rho = 2)
  expect_true(all(c("R", "rt") %in% names(r1)))
})

test_that("finite rho reaches the sampled likelihood adapter", {
  skip_on_cran()
  matchfun <- function(d) d$S == d$lR
  dat <- EMC2::forstmann[EMC2::forstmann$subjects ==
                           unique(EMC2::forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  dat <- dat[seq_len(min(24, nrow(dat))), ]
  des <- suppressMessages(EMC2::design(
    data = dat, model = function() EMC2::BAwD(gamma = .5, rho = 2),
    matchfun = matchfun,
    formula = list(mu ~ 1, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1),
    constants = c(ell = log(1))))
  p <- c(mu = .9, sigma = log(.6), B = log(.8), A = log(.3),
         t0 = log(.15), k = log(.8))[names(EMC2::sampled_pars(des))]
  set.seed(20260821)
  sim <- EMC2::make_data(p, design = des, n_trials = 12)
  e <- suppressMessages(EMC2::make_emc(sim, des, type = "single", n_chains = 1))
  model <- e[[1]]$model()
  expect_equal(model$c_name, "BAwD_LOGN_GAM12_RHO2")
  dadm <- e[[1]]$data[[1]]
  designs <- list()
  for (nm in names(model$p_types)) {
    designs[[nm]] <- attr(dadm, "designs")[[nm]][
      attr(attr(dadm, "designs")[[nm]], "expand"), , drop = FALSE]
  }
  ll <- EMC2:::calc_ll_oo(
    matrix(p, nrow = 1, dimnames = list(NULL, names(p))), dadm,
    constants = attr(dadm, "constants"), designs = designs,
    type = model$c_name, bounds = model$bound, transforms = model$transform,
    pretransforms = model$pre_transform, p_types = names(model$p_types),
    min_ll = log(1e-10), trend = model$trend)

  pars <- EMC2:::get_pars_matrix_oo(p, dadm, model)
  n_lR <- length(levels(dadm$lR))
  n_tr <- nrow(dadm) / n_lR
  lt <- numeric(n_tr)
  for (j in seq_len(n_tr)) {
    idx <- ((j - 1) * n_lR + 1):(j * n_lR)
    rt <- dadm$rt[idx[1]]
    if (!is.finite(rt) || !any(dadm$winner[idx])) {
      term <- 1
      for (i in idx) {
        fi <- rho_cpp_p(Inf, pars[i, "mu"], pars[i, "sigma"],
                        pars[i, "b"], pars[i, "A"], pars[i, "k"], pars[i, "ell"],
                        1L, .5, 2)
        term <- term * (1 - fi)
      }
    } else {
      w <- idx[which(dadm$winner[idx])]
      f <- rho_cpp_d(rt - pars[w, "t0"], pars[w, "mu"], pars[w, "sigma"],
                     pars[w, "b"], pars[w, "A"], pars[w, "k"], pars[w, "ell"],
                     1L, .5, 2)
      term <- f
      for (i in setdiff(idx, w)) {
        fi <- rho_cpp_p(rt - pars[i, "t0"], pars[i, "mu"], pars[i, "sigma"],
                        pars[i, "b"], pars[i, "A"], pars[i, "k"], pars[i, "ell"],
                        1L, .5, 2)
        term <- term * (1 - fi)
      }
    }
    lt[j] <- max(log(max(term, 0)), log(1e-10))
  }
  expect_equal(ll, sum(lt[attr(dadm, "expand")]), tolerance = 1e-8)
})
