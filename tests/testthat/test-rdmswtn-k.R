test_that("pswtn reduces to pwald when sv is zero", {
  t <- 0.8
  b <- 1.4
  v <- 1.1
  s <- 1.0
  lambda <- 0.35

  expect_equal(
    EMC2:::pswtn(t, b, v, sv = 0, s = s, lambda = lambda),
    EMC2:::pwald(t, b, v, sigma = s, A = 0, k = lambda),
    tolerance = 1e-12
  )
})

test_that("killed wald and rdmswtn infinite-time masses are defective", {
  p_wald_inf <- EMC2:::pwald(Inf, b = 2, mu = 3, sigma = 1, A = 0.4, k = 1, log_out = FALSE)
  p_rdmswtn_sv0_inf <- EMC2:::prdmswtn(Inf, b = 2, mu_drift = 3, A = 0.4, sv = 0,
                                s = 1, c = 1, lambda = 1, n_gauss_nodes = 20, log_out = FALSE)
  p_rdmswtn_sv1_inf <- EMC2:::prdmswtn(Inf, b = 2, mu_drift = 3, A = 0.4, sv = 1,
                                s = 1, c = 1, lambda = 1, n_gauss_nodes = 20, log_out = FALSE)

  expect_gt(p_wald_inf, 0)
  expect_lte(p_wald_inf, 1)
  expect_gt(p_rdmswtn_sv0_inf, 0)
  expect_lte(p_rdmswtn_sv0_inf, 1)
  expect_gt(p_rdmswtn_sv1_inf, 0)
  expect_lt(p_rdmswtn_sv1_inf, 1)
})

test_that("killed swtn cdf is locally consistent with the pdf", {
  t <- 0.9
  h <- 1e-6
  b <- 1.3
  v <- 1.0
  sv <- 0.45
  s <- 1.0
  lambda <- 0.6

  pdf <- EMC2:::dswtn(t, b, v, sv = sv, s = s, lambda = lambda)
  slope <- (EMC2:::pswtn(t + h, b, v, sv = sv, s = s, lambda = lambda) -
            EMC2:::pswtn(t - h, b, v, sv = sv, s = s, lambda = lambda)) / (2 * h)

  expect_equal(slope, pdf, tolerance = 2e-4)
})

test_that("rRDMSWTN returns omissions when killing wins the race", {
  set.seed(1)

  n_trials <- 200
  lR <- factor(rep(c("r1", "r2"), n_trials), levels = c("r1", "r2"))
  pars <- matrix(
    rep(c(v = 1.0, b = 1.1, A = 0.1, t0 = 0.2, sv = 0.25, lambda = 1e6), 2 * n_trials),
    ncol = 6,
    byrow = TRUE,
    dimnames = list(NULL, c("v", "b", "A", "t0", "sv", "lambda"))
  )

  sim <- EMC2:::rRDMSWTN(lR, pars)

  expect_true(all(is.na(sim$R)))
  expect_true(all(is.infinite(sim$rt)))
})

test_that("rRDMSWTN nests to the no-kill model when lambda is zero", {
  set.seed(1)

  designRDMSWTN <- design(
    factors = list(S = "Target", subjects = 1, L = c("L", "M", "H")),
    Rlevels = c("Go"),
    formula = list(v ~ L, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1, lambda ~ 1),
    constants = c(s = log(1)),
    model = RDMSWTN(),
    UC = 3,
    report_p_vector = FALSE
  )

  p_vec <- sampled_pars(designRDMSWTN)
  p_vec[1:3] <- log(c(3, 2, 1.5))
  p_vec["B"] <- log(1.1)
  p_vec["A"] <- log(0.4)
  p_vec["t0"] <- log(0.2)
  p_vec["sv"] <- log(0.25)
  p_vec["lambda"] <- log(0)

  expect_silent(dat <- make_data(p_vec, design = designRDMSWTN, n_trials = 200))
  expect_true(is.data.frame(dat))
  expect_true(all(is.finite(dat$rt)))
  expect_false(anyNA(dat$R))
})

test_that("RDMSWTN requires lambda in model-level parameter matrices", {
  rt <- c(0.5, 0.8)
  pars_lambda <- cbind(v = 1.1, b = 1.4, A = 0.2, t0 = 0.1, sv = 0.3, s = 1.0, lambda = 0.4)

  d1 <- EMC2:::dRDMSWTN(rt, pars_lambda)
  p1 <- EMC2:::pRDMSWTN(rt, pars_lambda)

  expect_true(all(is.finite(d1)))
  expect_true(all(is.finite(p1)))
  expect_error(EMC2:::dRDMSWTN(rt, cbind(v = 1.1, b = 1.4, A = 0.2, t0 = 0.1, sv = 0.3, s = 1.0, k = 0.4)), "requires parameter column 'lambda'")
  expect_error(EMC2:::pRDMSWTN(rt, cbind(v = 1.1, b = 1.4, A = 0.2, t0 = 0.1, sv = 0.3, s = 1.0, k = 0.4)), "requires parameter column 'lambda'")
})
