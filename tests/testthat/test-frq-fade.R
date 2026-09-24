fade_pars <- function(n, alpha = 2, beta = 3, lambda = 1.3, kappa = 0.7,
                      t0 = 0.12, delta = 0, cv_u = 0) {
  out <- cbind(alpha = rep(alpha, n), beta = rep(beta, n),
    lambda = rep(lambda, n), kappa = rep(kappa, n), t0 = rep(t0, n),
    delta = rep(delta, n), cv_u = rep(cv_u, n))
  out
}

fade_calc_ll <- function(dat, TC = list()) {
  matchfun <- function(d) as.character(d$S) == as.character(d$lR)
  des <- design(data = dat, model = FRQfade, matchfun = matchfun,
    formula = list(alpha ~ 1, beta ~ 1, lambda ~ 1, kappa ~ 1, t0 ~ 1),
    TC = TC)
  emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1,
    compress = FALSE, rt_resolution = NULL))
  dadm <- emc[[1]]$data[[1]]
  model <- emc[[1]]$model()
  p_types <- names(model$p_types)
  designs <- lapply(p_types, function(nm) {
    d <- attr(dadm, "designs")[[nm]]
    d[attr(d, "expand"), , drop = FALSE]
  })
  names(designs) <- p_types
  p <- model$p_types
  p["alpha"] <- 0
  p["beta"] <- log(2)
  p["lambda"] <- log(1.3)
  p["kappa"] <- log(0.8)
  p["t0"] <- log(0.1)
  ll <- EMC2:::calc_ll_oo(
    matrix(p, nrow = 1, dimnames = list(NULL, p_types)), dadm,
    constants = attr(dadm, "constants"), designs = designs,
    type = model$c_name, bounds = model$bound,
    transforms = model$transform, pretransforms = model$pre_transform,
    p_types = p_types, min_ll = log(1e-10), trend = model$trend)
  list(ll = as.numeric(ll), model = model, dadm = dadm, p = p)
}

test_that("zero fading rate is the proper FRQ boundary", {
  rt <- c(0.05, 0.15, 0.3, 0.8, 2, Inf)
  pars <- fade_pars(length(rt), alpha = 2.4, beta = 1.7,
    lambda = 0.95, kappa = 0, t0 = 0.1, delta = 0.6, cv_u = 0.35)
  old_pars <- cbind(alpha = pars[, "alpha"], beta = pars[, "beta"], h = 1,
    lambda = pars[, "lambda"], t0 = pars[, "t0"], delta = pars[, "delta"],
    cv_u = pars[, "cv_u"])

  expect_equal(dFRQfade(rt, pars), EMC2:::dFRQ(rt, old_pars), tolerance = 1e-12)
  expect_equal(pFRQfade(rt, pars), EMC2:::pFRQ(rt, old_pars), tolerance = 1e-12)
  expect_equal(sFRQfade(rt, pars), EMC2:::sFRQ(rt, old_pars), tolerance = 1e-12)
  summary <- EMC2:::frq_fade_summary(2.4, 1.7, 0.95, 0, 0.6, 0.35)
  expect_equal(unname(summary[1, "q_inf"]), 1)
  expect_equal(unname(summary[1, "h"]), 1)
})

test_that("fading opportunity gives its derived defective endpoint", {
  pars <- fade_pars(1, alpha = 1.6, beta = 2.8, lambda = 1.1,
    kappa = 0.75, t0 = 0, delta = 0.45, cv_u = 0.3)
  summary <- EMC2:::frq_fade_summary(1.6, 2.8, 1.1, 0.75, 0.45, 0.3)
  h <- unname(summary[1, "h"])

  expect_equal(pFRQfade(Inf, pars), h, tolerance = 1e-13)
  expect_equal(sFRQfade(Inf, pars), 1 - h, tolerance = 1e-13)
  mass <- integrate(function(u) EMC2:::dfrqfade(
    u, 1.6, 2.8, 1.1, 0.75, 0.45, 0.3), 0, Inf,
    rel.tol = 1e-8)$value
  expect_equal(mass, h, tolerance = 2e-7)

  expect_false(isTRUE(all.equal(
    unname(EMC2:::frq_fade_summary(1, 1, 1.1, 0.75, 0, 0)[1, "h"]),
    unname(EMC2:::frq_fade_summary(1, 1, 1.1, 0.75, 0.8, 0)[1, "h"]))))
  expect_false(isTRUE(all.equal(
    unname(EMC2:::frq_fade_summary(1, 1, 1.1, 0.75, 0, 0)[1, "h"]),
    unname(EMC2:::frq_fade_summary(1, 1, 1.1, 0.75, 0, 0.5)[1, "h"]))))

  tiny_kappa <- fade_pars(1, alpha = 1.6, beta = 2.8, lambda = 1.1,
    kappa = 1e-12, t0 = 0, delta = 0.45, cv_u = 0.3)
  zero_kappa <- tiny_kappa
  zero_kappa[, "kappa"] <- 0
  probe <- c(0.1, 0.5, 2)
  tiny_kappa <- tiny_kappa[rep(1, length(probe)), , drop = FALSE]
  zero_kappa <- zero_kappa[rep(1, length(probe)), , drop = FALSE]
  expect_equal(dFRQfade(probe, tiny_kappa), dFRQfade(probe, zero_kappa),
               tolerance = 2e-11)
  expect_equal(pFRQfade(probe, tiny_kappa), pFRQfade(probe, zero_kappa),
               tolerance = 2e-11)
})

test_that("zero frailty and zero threshold variation use the closed form", {
  u <- c(0.01, 0.2, 0.8, 2)
  lambda <- 1.25
  kappa <- 0.65
  A <- -expm1(-kappa * u) / kappa
  q <- -expm1(-lambda * A)
  density <- lambda * exp(-kappa * u - lambda * A)

  expect_equal(EMC2:::pfrqfade(u, 1, 1, lambda, kappa, 0, 0), q,
               tolerance = 1e-13)
  expect_equal(EMC2:::dfrqfade(u, 1, 1, lambda, kappa, 0, 0), density,
               tolerance = 1e-13)
})

test_that("FRQfade enforces fitted shape and rate support", {
  model <- FRQfade()
  expect_identical(model$c_name, "FRQ_FADE")
  expect_identical(model$p_types_canonical,
    c("alpha", "beta", "lambda", "kappa", "t0"))
  expect_false("h" %in% names(model$p_types))
  expect_equal(unname(model$bound$minmax[1, "alpha"]), 1)
  expect_equal(unname(model$bound$exception[["kappa"]]), 0)
  expect_equal(unname(EMC2:::frq_fade_summary(0.9, 2, 1, 0.5, 0, 0)[1, "h"]),
               NA_real_)
  pars <- cbind(alpha = 2, beta = 3, lambda = 1.3, kappa = 0.8,
                t0 = 0, delta = 0, cv_u = 0)
  mapped <- model$Ttransform(pars, NULL)
  expect_true(all(c("h", "q_inf", "N", "d", "sQ", "a_u") %in%
                  colnames(mapped)))
  expect_equal(unname(mapped[1, "h"]),
    unname(EMC2:::frq_fade_summary(2, 3, 1.3, 0.8, 0, 0)[1, "h"]))
})

test_that("the compiled default simulator follows completion and RT CDFs", {
  skip_on_cran()
  withr::local_options(emc2.cpp_rfun = TRUE)
  n <- 24000
  lR <- factor(rep("go", n), levels = "go")
  pars <- fade_pars(n, alpha = 1, beta = 1, lambda = 1.2,
    kappa = 0.85, t0 = 0.1, delta = 0.65, cv_u = 0.25)
  set.seed(317)
  sim <- FRQfade()$rfun(list(lR = lR), pars)
  h <- EMC2:::frq_fade_summary(1, 1, 1.2, 0.85, 0.65, 0.25)[1, "h"]
  probes <- c(0.3, 0.6, 1)
  analytic <- pfrqfade(probes - 0.1, 1, 1, 1.2, 0.85, 0.65, 0.25)
  empirical <- vapply(probes, function(x) mean(is.finite(sim$rt) & sim$rt <= x),
                      numeric(1))

  expect_lt(abs(mean(is.finite(sim$rt)) - h), 0.012)
  expect_lt(max(abs(empirical - analytic)), 0.015)
  expect_true(all(is.na(sim$R[!is.finite(sim$rt)])))
  expect_true(all(sim$rt[is.finite(sim$rt)] > 0.1))

  reference <- withr::with_options(list(emc2.cpp_rfun = FALSE), {
    set.seed(318)
    FRQfade()$rfun(list(lR = lR), pars)
  })
  ref_empirical <- vapply(probes, function(x)
    mean(is.finite(reference$rt) & reference$rt <= x), numeric(1))
  expect_lt(abs(mean(is.finite(reference$rt)) - h), 0.02)
  expect_lt(max(abs(ref_empirical - analytic)), 0.025)
})

test_that("FRQfade race likelihood preserves omission mass and finite truncation", {
  omission <- data.frame(subjects = factor("s1"),
    S = factor("left", levels = c("left", "right")),
    R = factor(NA_character_, levels = c("left", "right")), rt = Inf)
  got <- fade_calc_ll(omission)
  h <- EMC2:::frq_fade_summary(2, 3, 1.3, 0.8, 0, 0)[1, "h"]
  expect_equal(got$ll, 2 * log1p(-unname(h)), tolerance = 1e-8)

  # Both accumulators share (alpha, beta, lambda, kappa, t0) = (2, 3, 1.3,
  # 0.8, 0.1): the winner contributes f*S, and a finite UT divides by Z.
  finite <- data.frame(subjects = factor("s1"),
    S = factor(c("left", "right"), levels = c("left", "right")),
    R = factor(c("left", "right"), levels = c("left", "right")),
    rt = c(0.55, 0.9))
  UT <- 1.25
  f <- function(t) EMC2:::dfrqfade(t - 0.1, 2, 3, 1.3, 0.8, 0, 0)
  S <- function(t) EMC2:::pfrqfade(t - 0.1, 2, 3, 1.3, 0.8, 0, 0,
                                   lower_tail = FALSE)
  raw <- sum(log(f(finite$rt)) + log(S(finite$rt)))
  # Default: omissions are observable, so the S(Inf) atom stays in Z.
  truncated <- fade_calc_ll(finite, TC = list(UT = UT))
  expect_equal(truncated$ll, raw - 2 * log(1 - S(UT)^2 + S(Inf)^2),
               tolerance = 1e-8)
  # filter_defective drops the atom: Z is the probability of a response by UT.
  filtered <- fade_calc_ll(finite, TC = list(UT = UT, filter_defective = TRUE))
  expect_equal(filtered$ll, raw - 2 * log(1 - S(UT)^2), tolerance = 1e-8)
})

test_that("the compiled simulator defaults absent delta and cv_u to zero", {
  withr::local_options(emc2.cpp_rfun = TRUE)
  n <- 50
  lR <- factor(rep("go", n), levels = "go")
  full <- fade_pars(n)
  bare <- full[, c("alpha", "beta", "lambda", "kappa", "t0")]
  set.seed(11)
  a <- EMC2:::.rfun_FRQfade(lR, full)
  set.seed(11)
  b <- EMC2:::.rfun_FRQfade(lR, bare)
  expect_identical(a, b)
})
