test_that("REXG dfun/pfun use shifted zero-truncated ex-Gaussian", {
  rt <- c(-0.2, 0.1, 0.2, 0.5, 1.0, 2.0, Inf)
  mu <- c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40)
  sigma <- c(0.15, 0.20, 0.12, 0.18, 0.25, 0.30, 0.22)
  tau <- c(0.20, 0.25, 0.30, 0.18, 0.22, 0.35, 0.28)
  t0 <- rep(0.1, length(rt))
  pars <- cbind(mu = mu, sigma = sigma, tau = tau, t0 = t0)
  trunc_pars <- cbind(mu = mu, sigma = sigma, tau = tau, exg_lb = 0)
  model <- REXG()

  d <- model$dfun(rt, pars)
  p <- model$pfun(rt, pars)
  at_or_below <- rt <= t0
  above <- rt > t0 & is.finite(rt)

  expect_equal(d[at_or_below], rep(0, sum(at_or_below)))
  expect_equal(p[at_or_below], rep(0, sum(at_or_below)))
  expect_equal(1 - p[at_or_below], rep(1, sum(at_or_below)))
  expect_equal(d[above], EMC2:::dtexGaussian(rt[above] - t0[above],
                                              trunc_pars[above, , drop = FALSE]),
               tolerance = 1e-12)
  expect_equal(p[above], EMC2:::ptexGaussian(rt[above] - t0[above],
                                              trunc_pars[above, , drop = FALSE]),
               tolerance = 1e-12)
  expect_equal(d[length(d)], 0)
  expect_equal(p[length(p)], 1)
  invalid_t0 <- pars[1:3, , drop = FALSE]
  invalid_t0[, "t0"] <- c(NA, Inf, -1)
  expect_equal(model$pfun(rep(Inf, 3), invalid_t0), numeric(3))
})

test_that("REXG uses the common 0.05-second t0 lower bound", {
  model <- REXG()
  expect_equal(unname(model$bound$minmax[1, "t0"]), 0.05)
  # The 0.05 floor is the point; t0 = 0 stays an exception, as in every other
  # race model.  Without it REXG's own p_types default (log(0)) is outside its
  # bound, so a design that leaves t0 out of the formula fixes it at a value
  # do_bound() rejects and the likelihood floors for every proposal.
  expect_equal(unname(model$bound$exception[["t0"]]), 0)
})

test_that("REXG simulation draws positive processes and races shifted finishes", {
  model <- REXG()
  lR <- factor(rep(c("a", "b"), 40), levels = c("a", "b"))
  pars <- cbind(
    mu = rep(0.35, length(lR)),
    sigma = rep(0.12, length(lR)),
    tau = rep(0.18, length(lR)),
    t0 = rep(c(0.6, 0.05), 40)
  )
  attr(pars, "ok") <- rep(TRUE, nrow(pars))

  set.seed(11)
  out <- model$rfun(data = list(lR = lR), pars = pars)

  set.seed(11)
  dt <- EMC2:::rtexG(
    length(lR), mu = pars[, "mu"], sigma = pars[, "sigma"],
    tau = pars[, "tau"], lb = rep(0, length(lR))
  )
  finish <- matrix(dt, nrow = 2) + matrix(pars[, "t0"], nrow = 2)
  winner <- max.col(-t(finish), ties.method = "first")
  process_winner <- max.col(-t(matrix(dt, nrow = 2)), ties.method = "first")
  pick <- cbind(winner, seq_len(ncol(finish)))
  expected_R <- factor(levels(lR)[winner], levels = levels(lR))
  winning_t0 <- matrix(pars[, "t0"], nrow = 2)[pick]

  expect_true(all(dt > 0))
  expect_true(any(winner != process_winner))
  expect_identical(out$R, expected_R)
  expect_equal(out$rt, finish[pick], tolerance = 1e-12)
  expect_true(all(out$rt > winning_t0))
})

test_that("REXG compiled likelihood includes zero-truncation normalizer", {
  lI_fun <- function(d) factor(rep(1, nrow(d)), levels = 1)
  design_rexg <- design(
    factors = list(subjects = 1, S = 1),
    Rlevels = 1,
    formula = list(mu ~ 1, sigma ~ 1, tau ~ 1, t0 ~ 1),
    functions = list(lI = lI_fun),
    model = plain_link(REXG)
  )
  dat <- data.frame(
    subjects = factor(1), S = factor(1), R = factor(1, levels = 1), rt = 0.75
  )
  dadm <- EMC2:::design_model(dat, design_rexg, compress = FALSE, rt_resolution = NULL)
  model <- attr(dadm, "model")()
  p <- c(mu = log(0.35), sigma = log(0.12), tau = log(0.18), t0 = log(0.1))
  p_types <- names(model$p_types)
  designs <- list()
  for (nm in p_types) {
    designs[[nm]] <- attr(dadm, "designs")[[nm]][
      attr(attr(dadm, "designs")[[nm]], "expand"), , drop = FALSE
    ]
  }

  ll <- EMC2:::calc_ll_oo(
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
  trunc_pars <- cbind(mu = 0.35, sigma = 0.12, tau = 0.18, exg_lb = 0)
  expected <- log(EMC2:::dtexGaussian(0.75 - 0.1, trunc_pars))

  expect_equal(as.numeric(ll), expected, tolerance = 1e-10)
})

test_that("REXG race truncation evaluates shared endpoint tuples", {
  data("forstmann", package = "EMC2")
  lI_fun <- function(d) factor(rep(1, nrow(d)), levels = 1)
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ,
                   drop = FALSE]
  dat$subjects <- droplevels(dat$subjects)
  des <- design(
    data = dat, model = REXG, matchfun = matchfun,
    formula = list(mu ~ 1, sigma ~ 1, tau ~ 1, t0 ~ 1),
    functions = list(lI = lI_fun), UT = 2
  )
  emc <- make_emc(dat, des, type = "single", n_chains = 1, compress = FALSE)
  dadm <- emc[[1L]]$data[[1L]]
  model <- emc[[1L]]$model()
  p_types <- names(model$p_types)
  designs <- lapply(p_types, function(nm) {
    d <- attr(dadm, "designs")[[nm]]
    d[attr(d, "expand"), , drop = FALSE]
  })
  names(designs) <- p_types
  p <- stats::setNames(rep(-5, length(p_types)), p_types)
  sampled <- rownames(emc[[1L]]$samples$alpha)
  p[sampled] <- as.numeric(emc[[1L]]$samples$alpha[, 1L, 1L])
  particles <- matrix(rep(p, each = 4L), nrow = 4L,
                      dimnames = list(NULL, p_types))

  ll <- EMC2:::calc_ll_oo(
    particles, dadm, constants = attr(dadm, "constants"), designs = designs,
    type = model$c_name, bounds = model$bound, transforms = model$transform,
    pretransforms = model$pre_transform, p_types = p_types,
    min_ll = log(1e-10), trend = model$trend
  )
  expect_true(all(is.finite(ll)))
  expect_identical(ll, rep(ll[[1L]], length(ll)))
})

test_that("REXG compiled log density stays accurate through both tails", {
  lI_fun <- function(d) factor(rep(1, nrow(d)), levels = 1)
  model <- plain_link(REXG)()
  des <- design(
    factors = list(subjects = 1, S = 1), Rlevels = 1,
    formula = list(mu ~ 1, sigma ~ 1, tau ~ 1, t0 ~ 1),
    functions = list(lI = lI_fun), model = plain_link(REXG)
  )
  mu <- 0.4; sigma <- 0.03; tau <- 0.1; t0 <- 0.05
  # Includes the old z = (x - mu) / sigma - sigma / tau approximation
  # boundary and points deep enough into the left tail to expose cancellation.
  x <- c(1e-4, 0.01, 0.03, 0.169, 0.4, 0.55, 1.3, 5)
  log_sum_exp_R <- function(a, b) {
    h <- max(a, b)
    h + log(exp(a - h) + exp(b - h))
  }
  log_survival_zero <- log_sum_exp_R(
    pnorm(-mu / sigma, lower.tail = FALSE, log.p = TRUE),
    mu / tau + 0.5 * (sigma / tau)^2 +
      pnorm(-mu / sigma - sigma / tau, log.p = TRUE)
  )
  expected <- vapply(x, function(xx) {
    y <- (xx - mu) / sigma
    -log(tau) + (mu - xx) / tau + 0.5 * (sigma / tau)^2 +
      pnorm(y - sigma / tau, log.p = TRUE) - log_survival_zero
  }, numeric(1))

  p <- c(mu = log(mu), sigma = log(sigma), tau = log(tau), t0 = log(t0))
  p_types <- names(model$p_types)
  for (i in seq_along(x)) {
    dat <- data.frame(subjects = factor(1), S = factor(1),
                      R = factor(1, levels = 1), rt = t0 + x[i])
    invisible(capture.output({
      dadm <- EMC2:::design_model(dat, des, compress = FALSE, rt_resolution = NULL)
      model_i <- attr(dadm, "model")()
      designs <- lapply(p_types, function(nm) {
        d <- attr(dadm, "designs")[[nm]]
        d[attr(d, "expand"), , drop = FALSE]
      })
      names(designs) <- p_types
      ll <- EMC2:::calc_ll_oo(
        matrix(p, nrow = 1, dimnames = list(NULL, names(p))), dadm,
        constants = attr(dadm, "constants"), designs = designs,
        type = model_i$c_name, bounds = model_i$bound,
        transforms = model_i$transform, pretransforms = model_i$pre_transform,
        p_types = p_types, min_ll = -1000, trend = model_i$trend
      )
    }))
    expect_equal(as.numeric(ll), expected[i], tolerance = 2e-10)
  }
})

test_that("REXG compiled density handles small scales and extreme ratios", {
  lI_fun <- function(d) factor(rep(1, nrow(d)), levels = 1)
  model_factory <- plain_link(REXG)
  des <- design(
    factors = list(subjects = 1, S = 1), Rlevels = 1,
    formula = list(mu ~ 1, sigma ~ 1, tau ~ 1, t0 ~ 1),
    functions = list(lI = lI_fun), model = model_factory
  )
  cases <- list(
    list(mu = 0.4, sigma = 0.03, tau = 0.1,
         x = c(0.0001, 0.169, 0.38877, 0.4, 1, 5)),
    list(mu = 0.4, sigma = 0.1, tau = 0.001,
         x = c(0.001, 0.1, 0.4, 0.6, 1)),
    list(mu = 0.4, sigma = 1e-4, tau = 2e-4,
         x = c(0.0001, 0.3995, 0.4, 0.4005, 0.5)),
    list(mu = 0.4, sigma = 1e-6, tau = 1e-6,
         x = c(0.399995, 0.4, 0.400005))
  )
  log_sum_exp_R <- function(a, b) {
    h <- max(a, b)
    h + log(exp(a - h) + exp(b - h))
  }

  for (case in cases) {
    mu <- case$mu; sigma <- case$sigma; tau <- case$tau; t0 <- 0.05
    log_survival_zero <- log_sum_exp_R(
      pnorm(-mu / sigma, lower.tail = FALSE, log.p = TRUE),
      mu / tau + 0.5 * (sigma / tau)^2 +
        pnorm(-mu / sigma - sigma / tau, log.p = TRUE)
    )
    expected <- vapply(case$x, function(xx) {
      y <- (xx - mu) / sigma
      -log(tau) + (mu - xx) / tau + 0.5 * (sigma / tau)^2 +
        pnorm(y - sigma / tau, log.p = TRUE) - log_survival_zero
    }, numeric(1))
    p <- c(mu = log(mu), sigma = log(sigma), tau = log(tau), t0 = log(t0))
    p_types <- names(model_factory()$p_types)

    for (i in seq_along(case$x)) {
      dat <- data.frame(subjects = factor(1), S = factor(1),
                        R = factor(1, levels = 1), rt = t0 + case$x[i])
      invisible(capture.output({
        dadm <- EMC2:::design_model(dat, des, compress = FALSE,
                                     rt_resolution = NULL)
        model_i <- attr(dadm, "model")()
        designs <- lapply(p_types, function(nm) {
          d <- attr(dadm, "designs")[[nm]]
          d[attr(d, "expand"), , drop = FALSE]
        })
        names(designs) <- p_types
        ll <- EMC2:::calc_ll_oo(
          matrix(p, nrow = 1, dimnames = list(NULL, names(p))), dadm,
          constants = attr(dadm, "constants"), designs = designs,
          type = model_i$c_name, bounds = model_i$bound,
          transforms = model_i$transform, pretransforms = model_i$pre_transform,
          p_types = p_types, min_ll = -1e20, trend = model_i$trend
        )
      }))
      expect_equal(as.numeric(ll), expected[i], tolerance = 2e-10)
    }
  }
})
