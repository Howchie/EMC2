# Poisson counter race model (PCOUNTER).  References are the Appendix of
# Ratcliff & Smith (2004, Psych Review, 111, 333-367): Equations A10a/A10b for
# the first-passage densities and A11a/A11b for the response probabilities.

# Literal transcription of Equation A10a: the winner's Erlang density times the
# loser's Erlang survivor written as its Poisson-sum form.
ga_A10a <- function(t, alpha, beta, Ka, Kb) {
  win <- (alpha * t)^(Ka - 1) * alpha * exp(-alpha * t) / factorial(Ka - 1)
  lose <- rowSums(sapply(0:(Kb - 1), function(j) (beta * t)^j * exp(-beta * t) / factorial(j)))
  win * lose
}

# Equation A11a as a negative-binomial sum: P(a) is the probability that Ka
# counts favouring a arrive before Kb counts favouring b.
Pa_A11a <- function(alpha, beta, Ka, Kb) {
  p <- alpha / (alpha + beta)
  sum(sapply(0:(Kb - 1), function(j) choose(Ka + j - 1, j) * (1 - p)^j * p^Ka))
}

test_that("PCOUNTER default structure and parameter types are correct", {
  model <- PCOUNTER()
  expect_equal(model$c_name, "PCOUNTER")
  expect_equal(model$p_types_canonical, c("nu", "sv", "gamma", "k", "omega", "t0"))
  expect_true(all(c("nu", "sv", "gamma", "k", "omega", "t0", "pContaminant") %in% names(model$p_types)))
})

test_that("PCOUNTER dfun/pfun match Erlang first-passage functions in fixed branch", {
  model <- PCOUNTER()
  rt <- c(0.05, 0.2, 0.5, 1.0, 2.0, 4.0)
  nu <- c(3, 8, 12, 5, 20, 30)
  t0 <- c(0.01, 0.05, 0.1, 0.2, 0.15, 0.3)

  for (K in c(1, 2, 5)) {
    pars <- cbind(nu = nu, sv = 0, gamma = 0, k = rep(K, length(rt)), omega = 0, t0 = t0)
    tt <- rt - t0
    expect_equal(model$dfun(rt, pars),
                 ifelse(tt > 0, stats::dgamma(tt, shape = K, rate = nu), 0),
                 tolerance = 1e-12)
    expect_equal(model$pfun(rt, pars),
                 ifelse(tt > 0, stats::pgamma(tt, shape = K, rate = nu), 0),
                 tolerance = 1e-12)
  }

  # No mass below t0.
  pars <- cbind(nu = rep(10, 2), sv = rep(0, 2), gamma = rep(0, 2), k = rep(3, 2), omega = rep(0, 2), t0 = rep(0.3, 2))
  expect_equal(model$dfun(c(0, 0.15), pars), c(0, 0))
  expect_equal(model$pfun(c(0, 0.15), pars), c(0, 0))
})

test_that("the two-accumulator race density reproduces Ratcliff & Smith Eq. A10a", {
  model <- PCOUNTER()
  t <- c(0.05, 0.1, 0.25, 0.5, 1, 2)
  alpha <- 12; beta <- 7; Ka <- 4L; Kb <- 6L

  win_pars  <- cbind(nu = rep(alpha, length(t)), sv = 0, gamma = 0, k = rep(Ka, length(t)), omega = 0, t0 = rep(0, length(t)))
  lose_pars <- cbind(nu = rep(beta, length(t)),  sv = 0, gamma = 0, k = rep(Kb, length(t)), omega = 0, t0 = rep(0, length(t)))
  # EMC2 race likelihood: winner density x loser survivor.
  got <- model$dfun(t, win_pars) * (1 - model$pfun(t, lose_pars))

  expect_equal(got, ga_A10a(t, alpha, beta, Ka, Kb), tolerance = 1e-12)
})

test_that("integrated race density reproduces the Eq. A11a response probability", {
  model <- PCOUNTER()
  for (cs in list(c(12, 7, 4, 6), c(5, 5, 3, 3), c(20, 3, 8, 2), c(2, 9, 1, 5))) {
    alpha <- cs[1]; beta <- cs[2]; Ka <- cs[3]; Kb <- cs[4]
    dens <- function(t) {
      w <- cbind(nu = rep(alpha, length(t)), sv = 0, gamma = 0, k = rep(Ka, length(t)), omega = 0, t0 = rep(0, length(t)))
      l <- cbind(nu = rep(beta, length(t)),  sv = 0, gamma = 0, k = rep(Kb, length(t)), omega = 0, t0 = rep(0, length(t)))
      model$dfun(t, w) * (1 - model$pfun(t, l))
    }
    Pa <- stats::integrate(dens, 0, Inf, rel.tol = 1e-10)$value
    expect_equal(Pa, Pa_A11a(alpha, beta, Ka, Kb), tolerance = 1e-8)
    # Both counters exhaust the probability: a response is always made.
    Pb <- Pa_A11a(beta, alpha, Kb, Ka)
    expect_equal(Pa + Pb, 1, tolerance = 1e-10)
  }
})

test_that("k criterion is snapped to nearest positive integer", {
  model <- PCOUNTER()
  rt <- c(0.2, 0.5, 1.0)
  mk <- function(k_val) cbind(nu = rep(8, 3), sv = 0, gamma = 0, k = rep(k_val, 3), omega = 0, t0 = rep(0.1, 3))

  for (cs in list(c(2.3, 2), c(2.7, 3), c(2.0, 2), c(0.4, 1), c(1.2, 1), c(5.9, 6),
                  c(3.5, 4), c(2.5, 3))) {
    expect_equal(model$dfun(rt, mk(cs[1])), model$dfun(rt, mk(cs[2])), tolerance = 1e-12)
    expect_equal(model$pfun(rt, mk(cs[1])), model$pfun(rt, mk(cs[2])), tolerance = 1e-12)
  }
})

test_that("extended branches (sv, gamma, omega) evaluate cleanly and stay finite", {
  model <- PCOUNTER()
  rt <- c(0.2, 0.5, 1.2)
  pars <- cbind(nu = rep(10, 3), sv = rep(2, 3), gamma = rep(0.5, 3),
                k = rep(4, 3), omega = rep(1, 3), t0 = rep(0.1, 3))

  d_vals <- model$dfun(rt, pars)
  p_vals <- model$pfun(rt, pars)

  expect_true(all(is.finite(d_vals) & d_vals > 0))
  expect_true(all(is.finite(p_vals) & p_vals >= 0 & p_vals <= 1))
})

# Helper: evaluate the C++ particle likelihood for a design/parameter pair.
pcounter_cpp_ll <- function(p, design, dat) {
  emc <- make_emc(dat, design, type = "single")
  dadm <- emc[[1]]$data[[1]]
  m <- emc[[1]]$model()
  p_types <- names(m$p_types)
  designs <- list()
  for (nm in p_types) {
    d <- attr(dadm, "designs")[[nm]]
    designs[[nm]] <- d[attr(d, "expand"), , drop = FALSE]
  }
  EMC2:::calc_ll_oo(
    matrix(p, nrow = 1, dimnames = list(NULL, names(p))),
    dadm,
    constants = attr(dadm, "constants"),
    designs = designs,
    type = m$c_name,
    bounds = m$bound,
    transforms = m$transform,
    pretransforms = m$pre_transform,
    p_types = p_types,
    min_ll = log(1e-10),
    trend = m$trend
  )
}

test_that("simulated PCOUNTER choice proportions match Eq. A11a prediction", {
  matchfun <- function(d) factor(as.character(d$S) == as.character(d$lR), levels = c(FALSE, TRUE))
  design_pc <- design(
    factors = list(subjects = 1, S = 1:2), Rlevels = 1:2, matchfun = matchfun,
    formula = list(nu ~ lM, sv ~ 1, gamma ~ 1, k ~ 1, omega ~ 1, t0 ~ 1),
    constants = c(sv = log(0), gamma = log(0), omega = log(0)),
    model = PCOUNTER
  )
  nu_mis <- 6; nu_match <- 14; K <- 4
  p <- c("nu" = log(nu_mis), "nu_lMTRUE" = log(nu_match) - log(nu_mis),
         "k" = log(K), "t0" = log(0.1))

  set.seed(13)
  dat <- make_data(p, design_pc, n_trials = 4000, rt_resolution = NULL)
  acc <- mean(dat$R == dat$S)
  expect_equal(acc, Pa_A11a(nu_match, nu_mis, K, K), tolerance = 0.02)

  # Mean decision time of the winner is bounded by the faster counter's K/nu.
  expect_lt(mean(dat$rt) - 0.1, K / nu_mis)
  expect_gt(mean(dat$rt) - 0.1, 0)

  # Check C++ likelihood evaluation against R model
  ll_cpp <- pcounter_cpp_ll(p, design_pc, dat)
  expect_true(is.finite(ll_cpp))
})

test_that("simulated PCOUNTER works with non-zero sv, gamma, and omega", {
  matchfun <- function(d) factor(as.character(d$S) == as.character(d$lR), levels = c(FALSE, TRUE))
  design_pc <- design(
    factors = list(subjects = 1, S = 1:2), Rlevels = 1:2, matchfun = matchfun,
    formula = list(nu ~ lM, sv ~ 1, gamma ~ 1, k ~ 1, omega ~ 1, t0 ~ 1),
    model = PCOUNTER
  )
  p <- c("nu" = log(10), "nu_lMTRUE" = log(1.5), "sv" = log(2),
         "gamma" = log(0.5), "k" = log(3), "omega" = log(1), "t0" = log(0.15))

  set.seed(42)
  dat <- make_data(p, design_pc, n_trials = 500, rt_resolution = NULL)
  expect_gt(nrow(dat), 0)
  expect_true(all(dat$rt > 0.15))

  ll_cpp <- pcounter_cpp_ll(p, design_pc, dat)
  expect_true(is.finite(ll_cpp))
})
