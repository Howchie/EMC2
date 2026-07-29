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

test_that("PCOUNTER dfun/pfun are the Erlang first-passage functions with t0 shift", {
  model <- PCOUNTER()
  rt <- c(0.05, 0.2, 0.5, 1.0, 2.0, 4.0)
  alpha <- c(3, 8, 12, 5, 20, 30)
  t0 <- c(0.01, 0.05, 0.1, 0.2, 0.15, 0.3)

  for (K in c(1, 2, 5)) {
    pars <- cbind(alpha = alpha, K = rep(K, length(rt)), t0 = t0)
    tt <- rt - t0
    expect_equal(model$dfun(rt, pars),
                 ifelse(tt > 0, stats::dgamma(tt, shape = K, rate = alpha), 0),
                 tolerance = 1e-12)
    expect_equal(model$pfun(rt, pars),
                 ifelse(tt > 0, stats::pgamma(tt, shape = K, rate = alpha), 0),
                 tolerance = 1e-12)
  }

  # No mass below t0.
  pars <- cbind(alpha = rep(10, 2), K = rep(3, 2), t0 = rep(0.3, 2))
  expect_equal(model$dfun(c(0, 0.15), pars), c(0, 0))
  expect_equal(model$pfun(c(0, 0.15), pars), c(0, 0))
})

test_that("the two-accumulator race density reproduces Ratcliff & Smith Eq. A10a", {
  model <- PCOUNTER()
  t <- c(0.05, 0.1, 0.25, 0.5, 1, 2)
  alpha <- 12; beta <- 7; Ka <- 4L; Kb <- 6L

  win_pars  <- cbind(alpha = rep(alpha, length(t)), K = rep(Ka, length(t)), t0 = rep(0, length(t)))
  lose_pars <- cbind(alpha = rep(beta, length(t)),  K = rep(Kb, length(t)), t0 = rep(0, length(t)))
  # EMC2 race likelihood: winner density x loser survivor.
  got <- model$dfun(t, win_pars) * (1 - model$pfun(t, lose_pars))

  expect_equal(got, ga_A10a(t, alpha, beta, Ka, Kb), tolerance = 1e-12)
})

test_that("integrated race density reproduces the Eq. A11a response probability", {
  model <- PCOUNTER()
  for (cs in list(c(12, 7, 4, 6), c(5, 5, 3, 3), c(20, 3, 8, 2), c(2, 9, 1, 5))) {
    alpha <- cs[1]; beta <- cs[2]; Ka <- cs[3]; Kb <- cs[4]
    dens <- function(t) {
      w <- cbind(alpha = rep(alpha, length(t)), K = rep(Ka, length(t)), t0 = rep(0, length(t)))
      l <- cbind(alpha = rep(beta, length(t)),  K = rep(Kb, length(t)), t0 = rep(0, length(t)))
      model$dfun(t, w) * (1 - model$pfun(t, l))
    }
    Pa <- stats::integrate(dens, 0, Inf, rel.tol = 1e-10)$value
    expect_equal(Pa, Pa_A11a(alpha, beta, Ka, Kb), tolerance = 1e-8)
    # Both counters exhaust the probability: a response is always made.
    Pb <- Pa_A11a(beta, alpha, Kb, Ka)
    expect_equal(Pa + Pb, 1, tolerance = 1e-10)
  }
})

test_that("integer_K snaps the criterion to the nearest count", {
  cont <- PCOUNTER()
  int  <- PCOUNTER(integer_K = TRUE)
  rt <- c(0.2, 0.5, 1.0)
  mk <- function(K) cbind(alpha = rep(8, 3), K = rep(K, 3), t0 = rep(0.1, 3))

  # Rounding, not ceiling, so each integer sits at the centre of its interval.
  for (cs in list(c(2.3, 2), c(2.7, 3), c(2.0, 2), c(0.4, 1), c(1.2, 1), c(5.9, 6),
                  c(3.5, 4), c(2.5, 3))) {
    expect_equal(int$dfun(rt, mk(cs[1])), cont$dfun(rt, mk(cs[2])), tolerance = 1e-12)
    expect_equal(int$pfun(rt, mk(cs[1])), cont$pfun(rt, mk(cs[2])), tolerance = 1e-12)
  }
  # Continuous mode leaves K alone.
  expect_false(isTRUE(all.equal(cont$dfun(rt, mk(2.3)), cont$dfun(rt, mk(2)))))
  expect_equal(int$c_name, "PCOUNTER_INTK")
  expect_equal(cont$c_name, "PCOUNTER")
})

# Helper: evaluate the C++ particle likelihood for a design/parameter pair.
pcounter_cpp_ll <- function(p, design, dat) {
  emc <- make_emc(dat, design, type = "single")
  dadm <- emc[[1]]$data[[1]]
  model <- emc[[1]]$model()
  p_types <- names(model$p_types)
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
    type = model$c_name,
    bounds = model$bound,
    transforms = model$transform,
    pretransforms = model$pre_transform,
    p_types = p_types,
    min_ll = log(1e-10),
    trend = model$trend
  )
}

test_that("PCOUNTER C++ likelihood equals the equivalent RGAMMA likelihood", {
  matchfun <- function(d) d$S == d$lR
  design_pc <- design(
    factors = list(subjects = 1, S = 1:2), Rlevels = 1:2, matchfun = matchfun,
    formula = list(alpha ~ lM, K ~ 1, t0 ~ 1), model = PCOUNTER
  )
  design_rg <- design(
    factors = list(subjects = 1, S = 1:2), Rlevels = 1:2, matchfun = matchfun,
    formula = list(lambda ~ lM, shape ~ 1, shift ~ 1), model = RGAMMA
  )

  p_pc <- c("alpha" = log(8), "alpha_lMTRUE" = 0.7, "K" = log(4), "t0" = log(0.15))
  p_rg <- c("lambda" = log(8), "lambda_lMTRUE" = 0.7, "shape" = log(4), "shift" = log(0.15))
  expect_equal(names(sampled_pars(design_pc)),
               sub("^shift", "t0", sub("^shape", "K", sub("^lambda", "alpha", names(sampled_pars(design_rg))))))

  set.seed(11)
  dat <- make_data(p_pc, design_pc, n_trials = 60)
  expect_true(all(is.finite(dat$rt)))
  expect_gt(min(dat$rt), 0.15)

  ll_pc <- pcounter_cpp_ll(p_pc, design_pc, dat)
  ll_rg <- pcounter_cpp_ll(p_rg, design_rg, dat)
  expect_true(is.finite(ll_pc))
  expect_equal(ll_pc, ll_rg, tolerance = 1e-10)
})

test_that("PCOUNTER_INTK C++ likelihood equals RGAMMA at the snapped criterion", {
  design_pc <- design(
    factors = list(subjects = 1, S = 1), Rlevels = 1:2,
    formula = list(alpha ~ 1, K ~ 1, t0 ~ 1), model = PCOUNTER(integer_K = TRUE)
  )
  design_rg <- design(
    factors = list(subjects = 1, S = 1), Rlevels = 1:2,
    formula = list(lambda ~ 1, shape ~ 1, shift ~ 1), model = RGAMMA
  )

  K_raw <- 3.4
  p_pc <- c("alpha" = log(9), "K" = log(K_raw), "t0" = log(0.12))
  p_rg <- c("lambda" = log(9), "shape" = log(floor(K_raw + 0.5)), "shift" = log(0.12))

  set.seed(12)
  dat <- make_data(p_pc, design_pc, n_trials = 60)

  expect_equal(pcounter_cpp_ll(p_pc, design_pc, dat),
               pcounter_cpp_ll(p_rg, design_rg, dat), tolerance = 1e-10)

  # The snap is a step: any K in [2.5, 3.5) gives the same likelihood, and a K
  # in the next interval up does not.
  p_same <- p_pc; p_same["K"] <- log(2.6)
  p_next <- p_pc; p_next["K"] <- log(3.6)
  expect_equal(pcounter_cpp_ll(p_same, design_pc, dat),
               pcounter_cpp_ll(p_pc, design_pc, dat), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(pcounter_cpp_ll(p_next, design_pc, dat),
                                pcounter_cpp_ll(p_pc, design_pc, dat))))
})

test_that("integer_K exposes the snapped criterion as a derived parameter", {
  d <- design(
    factors = list(subjects = 1, S = 1), Rlevels = 1:2,
    formula = list(alpha ~ 1, K ~ 1, t0 ~ 1), model = PCOUNTER(integer_K = TRUE)
  )
  mp <- mapped_pars(d, c(alpha = log(9), K = log(3.4), t0 = log(0.12)))
  expect_true("K_int" %in% colnames(mp))
  expect_true(all(mp$K_int == 3))
  expect_true(all(mp$K == 3.4))

  # Continuous mode has no snapped column: K is already the criterion.
  d_cont <- design(
    factors = list(subjects = 1, S = 1), Rlevels = 1:2,
    formula = list(alpha ~ 1, K ~ 1, t0 ~ 1), model = PCOUNTER
  )
  expect_false("K_int" %in% colnames(mapped_pars(d_cont, c(alpha = log(9), K = log(3.4), t0 = log(0.12)))))
})

test_that("the snapped criterion is what get_pars()/recovery() summarise", {
  # recovery() reaches get_pars() through its dots, so map/add_recalculated
  # must surface K_int for BOTH the posterior samples and the true values --
  # otherwise recovery is scored on a raw K that is only identified up to its
  # interval.
  d <- design(
    factors = list(subjects = 1, S = 1), Rlevels = 1:2,
    formula = list(alpha ~ 1, K ~ 1, t0 ~ 1), model = PCOUNTER(integer_K = TRUE)
  )
  p <- c(alpha = log(9), K = log(3), t0 = log(0.12))
  set.seed(21)
  dat <- make_data(p, d, n_trials = 150)
  emc <- make_emc(dat, d, type = "single", n_chains = 1)
  emc <- fit(emc, stage = "preburn", iter = 60, verbose = FALSE)

  # remove_constants = FALSE matters here: a perfectly identified criterion is
  # constant across draws, and the default would drop K_int for exactly the
  # fits where it is best determined.
  est <- get_pars(emc, selection = "alpha", map = TRUE, add_recalculated = TRUE,
                  merge_chains = TRUE, by_subject = TRUE, remove_constants = FALSE)
  truth <- get_pars(emc, selection = "alpha", map = TRUE, add_recalculated = TRUE,
                    merge_chains = TRUE, by_subject = TRUE, remove_constants = FALSE,
                    true_pars = p)
  expect_true("K_int" %in% colnames(est[[1]][[1]]))
  expect_true("K_int" %in% colnames(truth[[1]][[1]]))
  # The true value is snapped the same way, so both sides are on the count scale.
  expect_true(all(truth[[1]][[1]][, "K_int"] == 3))
  expect_true(all(est[[1]][[1]][, "K_int"] %% 1 == 0))
})

test_that("simulated PCOUNTER choice proportions match the Eq. A11a prediction", {
  matchfun <- function(d) d$S == d$lR
  design_pc <- design(
    factors = list(subjects = 1, S = 1:2), Rlevels = 1:2, matchfun = matchfun,
    formula = list(alpha ~ lM, K ~ 1, t0 ~ 1), model = PCOUNTER
  )
  alpha_mis <- 6; alpha_match <- 14; K <- 4
  p <- c("alpha" = log(alpha_mis), "alpha_lMTRUE" = log(alpha_match) - log(alpha_mis),
         "K" = log(K), "t0" = log(0.1))

  set.seed(13)
  dat <- make_data(p, design_pc, n_trials = 4000, rt_resolution = NULL)
  acc <- mean(dat$R == dat$S)
  expect_equal(acc, Pa_A11a(alpha_match, alpha_mis, K, K), tolerance = 0.02)

  # Mean decision time of the winner is bounded by the faster counter's K/alpha.
  expect_lt(mean(dat$rt) - 0.1, K / alpha_mis)
  expect_gt(mean(dat$rt) - 0.1, 0)
})
