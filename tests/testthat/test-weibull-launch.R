test_that("Weibull launch constructors expose a common positive parameterisation", {
  skip_model_validation()

  models <- list(
    BAwL = EMC2::BAwL(drift_distribution = "weibull"),
    BAwD = EMC2::BAwD(drift_distribution = "weibull"),
    BAwDD = EMC2::BAwDD(drift_distribution = "weibull"),
    BAwDp = EMC2::BAwDp(drift_distribution = "weibull"),
    BAwF = EMC2::BAwF(drift_distribution = "weibull"),
    BAwR = EMC2::BAwR(drift_distribution = "weibull"),
    BTAwL = EMC2::BTAwL(drift_distribution = "weibull"),
    BTAwLTransient = EMC2::BTAwLTransient(drift_distribution = "weibull"),
    BTAwLSustained = EMC2::BTAwLSustained(drift_distribution = "weibull")
  )
  expect_true(all(vapply(models, function(m) grepl("_WEIB", m$c_name), logical(1))))
  expect_true(all(vapply(models, function(m)
    identical(unname(m$p_types_canonical[1:2]), c("shape", "mean")), logical(1))))
  expect_true(all(vapply(models, function(m)
    identical(unname(m$transform$func[1:2]), c("exp", "exp")), logical(1))))
  expect_error(EMC2::BAwD(drift_distribution = "weibull", posdrift = FALSE), "posdrift")
  expect_error(EMC2::BAwF(drift_distribution = "weibull", posdrift = FALSE), "posdrift")
  expect_error(EMC2::BAwR(drift_distribution = "weibull", posdrift = FALSE), "posdrift")
  expect_error(EMC2::BAwL(drift_distribution = "weibull", posdrift = FALSE), "posdrift")
  expect_error(EMC2::BTAwL(drift_distribution = "weibull", posdrift = FALSE), "posdrift")

  # Existing defaults remain unchanged while Weibull is opt-in.
  expect_identical(EMC2::BAwL()$drift_distribution, "normal")
  expect_identical(EMC2::BAwD()$drift_distribution, "lognormal")
  expect_identical(EMC2::BAwF()$drift_distribution, "lognormal")
  expect_identical(EMC2::BAwR()$drift_distribution, "lognormal")
  expect_identical(EMC2::BTAwL()$drift_distribution, "lognormal")
})

test_that("Weibull launch mean is the arithmetic mean", {
  shape <- 1.7
  target_mean <- 2.4
  scale <- target_mean / gamma(1 + 1 / shape)
  expect_equal(EMC2:::.weibull_scale_from_mean(shape, target_mean), scale,
               tolerance = 1e-14)
  set.seed(20260828)
  draws <- EMC2:::.rweibull_mean(100000L, shape, target_mean)
  expect_equal(base::mean(draws), target_mean, tolerance = 0.02)
})

test_that("all Weibull ballistic simulators use the compiled default path", {
  skip_model_validation()
  n_trial <- 40L
  lR <- factor(rep(c("left", "right"), n_trial), levels = c("left", "right"))
  ok <- rep(TRUE, length(lR))
  make_common <- function(n = length(lR))
    cbind(shape = rep(2, n), mean = rep(3, n), b = 1.2, A = .2,
          t0 = .1, k = .2)

  cases <- list(
    BAwL = list(fun = EMC2:::.rfun_BAwL,
                pars = cbind(make_common(), lambda_g = .1, lambda_k = .1),
                args = list(launch = 3L)),
    BAwD = list(fun = EMC2:::.rfun_BAwD,
                pars = cbind(make_common(), ell = .2),
                args = list(launch = 3L, gamma = 0, rho = Inf)),
    BAwDp = list(fun = EMC2:::.rfun_BAwDp,
                 pars = cbind(make_common(), lambda = .2),
                 args = list(launch = 3L)),
    BAwF = list(fun = EMC2:::.rfun_BAwF,
                pars = make_common(),
                args = list(launch = 3L, rho = Inf)),
    BAwR = list(fun = EMC2:::.rfun_BAwR,
                pars = cbind(make_common(), kappa = .2, p = 1.4),
                args = list(launch = 3L)),
    BTAwLTransient = list(fun = EMC2:::.rfun_BTAwL,
                          pars = cbind(make_common(), tau = 2),
                          args = list(launch = 3L, mode = "transient")),
    BTAwLSustained = list(fun = EMC2:::.rfun_BTAwL,
                          pars = cbind(make_common(), tau_s = .9),
                          args = list(launch = 3L, mode = "sustained")),
    BTAwL = list(fun = EMC2:::.rfun_BTAwL,
                 pars = cbind(make_common(), tau_s = .9, tau_t = .8, pi = .4),
                 args = list(launch = 3L, mode = "full"))
  )

  withr::with_options(list(emc2.cpp_rfun = TRUE), {
    for (case in cases) {
      out <- do.call(case$fun, c(list(lR = lR, pars = case$pars, ok = ok), case$args))
      expect_equal(nrow(out), n_trial, info = names(case))
      expect_true(all(is.na(out$R) == is.infinite(out$rt)), info = names(case))
      expect_true(any(is.finite(out$rt)), info = names(case))
    }
  })
})

test_that("Weibull BAwL closed form agrees with an independent start-point integral", {
  skip_model_validation()
  shape <- 1.7; mean <- 2.4; A <- .35; b <- 1.4; k <- .3
  scale <- mean / gamma(1 + 1 / shape)
  cdf_ref <- function(t) {
    integrate(function(z) {
      e <- exp(-k * t)
      w <- k * (b - z * e) / (1 - e)
      ifelse(w > 0, pweibull(w, shape, scale, lower.tail = FALSE), 1)
    }, 0, A, rel.tol = 1e-11)$value / A
  }
  tt <- c(.35, .6, .9)
  got <- EMC2:::pleakyba(tt, A, b, shape, mean, k, TRUE, 3L)
  expect_equal(got, vapply(tt, cdf_ref, numeric(1)), tolerance = 2e-7)
  expect_true(all(diff(got) >= -1e-10))
  expect_true(all(EMC2:::dleakyba(tt, A, b, shape, mean, k, TRUE, 3L) >= 0))
})

test_that("Weibull CDF/PDF kernels are finite and monotone across ballistic families", {
  skip_model_validation()
  tt <- 10^seq(-3, 3, length.out = 80)
  p <- list(
    BAwD = EMC2:::pbawd(tt, .5, 2, 2, 3, .2, .7, 3L, TRUE, FALSE, 0, Inf),
    BAwDp = EMC2:::pbawdp(tt, .5, 2, 2, 3, .2, .3, 3L, TRUE, FALSE),
    BAwF = EMC2:::pbawf(tt, .5, 2, 2, 3, .2, 3L, TRUE, FALSE, Inf),
    BAwR = EMC2:::pbawr(tt, .5, 2, 2, 3, .2, 1.4, 3L, TRUE, FALSE),
    BTAwLTransient = EMC2:::pBTAwLTransient(
      tt, data.frame(shape = rep(2, length(tt)), mean = 3, b = 2, A = .5,
                     t0 = 0, k = .2, tau = .9), launch = 3L),
    BTAwLSustained = EMC2:::pBTAwLSustained(
      tt, data.frame(shape = rep(2, length(tt)), mean = 3, b = 2, A = .5,
                     t0 = 0, k = .2, tau_s = .9), launch = 3L)
  )
  for (x in p) {
    expect_true(all(is.finite(x)))
    expect_true(all(x >= 0 & x <= 1))
    expect_true(all(diff(x) >= -1e-8))
  }
})

test_that("Weibull kernels handle extreme time limits gracefully", {
  skip_model_validation()
  tt_ext <- c(1e-12, 1e12, Inf)
  p_bawl <- EMC2:::pleakyba(tt_ext, .5, 2, 2, 3, .2, TRUE, 3L)
  expect_equal(p_bawl[1], 0)
  expect_true(is.finite(p_bawl[2]) && p_bawl[2] > 0 && p_bawl[2] <= 1)
  expect_equal(p_bawl[2], p_bawl[3], tolerance = 1e-10)

  p_bawd <- EMC2:::pbawd(tt_ext, .5, 2, 2, 3, .2, .7, 3L, TRUE, FALSE, 0, Inf)
  expect_equal(p_bawd[1], 0)
  expect_true(is.finite(p_bawd[2]) && p_bawd[2] > 0 && p_bawd[2] <= 1)
  expect_equal(p_bawd[2], p_bawd[3], tolerance = 1e-10)
})
