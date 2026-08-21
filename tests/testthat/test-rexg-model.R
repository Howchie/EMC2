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

test_that("REXG simulation draws positive processes and races shifted finishes", {
  model <- REXG()
  lR <- factor(rep(c("a", "b"), 40), levels = c("a", "b"))
  pars <- cbind(
    mu = rep(0.35, length(lR)),
    sigma = rep(0.12, length(lR)),
    tau = rep(0.18, length(lR)),
    t0 = rep(c(0.6, 0), 40)
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
    model = REXG
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
