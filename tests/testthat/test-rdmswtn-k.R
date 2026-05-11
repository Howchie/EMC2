test_that("pswtn reduces to pwald when sv is zero", {
  t <- 0.8
  b <- 1.4
  v <- 1.1
  s <- 1.0
  lambda <- 0.35

  expect_equal(
    EMC2:::pswtn(t, v, b, sv = 0, s = s, lambda_k = lambda),
    EMC2:::pwald(t, v, b, sigma = s, A = 0, lambda_k = lambda),
    tolerance = 1e-12
  )
})

test_that("prdmswtn forwards posdrift in the sv=0 Wald path", {
  t <- 0.8
  b <- 1.4
  A <- 0.25
  v <- -0.6
  s <- 1.0

  expect_equal(
    EMC2:::prdmswtn(t, v, b, A, sv = 0, s = s, posdrift = FALSE),
    EMC2:::pwald(t, v, b, A = A, sigma = s, posdrift = FALSE),
    tolerance = 1e-12
  )
})

test_that("killed wald and rdmswtn infinite-time masses are defective", {
  p_wald_inf <- EMC2:::pwald(Inf, mu = 3, b = 2, sigma = 1, A = 0.4, lambda_k = 1,
                              log_out = FALSE, posdrift = FALSE)
  p_rdmswtn_sv0_inf <- EMC2:::prdmswtn(Inf, mu_drift = 3, b = 2, A = 0.4, sv = 0,
                                        s = 1, lambda_k = 1, n_gauss_nodes = 20,
                                        log_out = FALSE, posdrift = FALSE)
  p_rdmswtn_sv1_inf <- EMC2:::prdmswtn(Inf, mu_drift = 3, b = 2, A = 0.4, sv = 1,
                                        s = 1, lambda_k = 1, n_gauss_nodes = 20,
                                        log_out = FALSE, posdrift = FALSE)

  expect_gt(p_wald_inf, 0)
  expect_lte(p_wald_inf, 1)
  expect_gt(p_rdmswtn_sv0_inf, 0)
  expect_lte(p_rdmswtn_sv0_inf, 1)
  expect_gt(p_rdmswtn_sv1_inf, 0)
  expect_lt(p_rdmswtn_sv1_inf, 1)
})

test_that("posdrift=TRUE truncates SWTN drift at zero; posdrift=FALSE gives hit mass", {
  b <- 1.5; A <- 0.4; v <- -0.3; sv <- 0.5; s <- 1.0

  # posdrift=TRUE: positive-drift truncated distribution integrates to 1
  expect_equal(EMC2:::prdmswtn(Inf, v, b, A, sv = sv, s = s, posdrift = TRUE),  1.0, tolerance = 1e-12)
  expect_equal(EMC2:::pswtn(Inf, v, b,    sv = sv, s = s, posdrift = TRUE),     1.0, tolerance = 1e-12)

  # posdrift=FALSE: defective mass, strictly less than 1 for negative mean drift
  p_rdm <- EMC2:::prdmswtn(Inf, v, b, A, sv = sv, s = s, posdrift = FALSE)
  p_swtn <- EMC2:::pswtn(Inf, v, b,   sv = sv, s = s, posdrift = FALSE)
  expect_gt(p_rdm,  0); expect_lt(p_rdm,  1)
  expect_gt(p_swtn, 0); expect_lt(p_swtn, 1)

  # Literal positive-drift truncation is not conditional-hit normalization.
  t <- 1.2
  p_def  <- EMC2:::prdmswtn(t, v, b, A, sv = sv, s = s, posdrift = FALSE)
  p_norm <- EMC2:::prdmswtn(t, v, b, A, sv = sv, s = s, posdrift = TRUE)
  expect_true(p_norm > 0 && p_norm < 1)
  expect_false(isTRUE(all.equal(p_norm, p_def / p_rdm, tolerance = 1e-4)))
})

test_that("dswtn/drdmswtn posdrift=TRUE integrates to 1", {
  b <- 1.2; A <- 0.3; v <- -0.2; sv <- 0.4; s <- 1.0

  # Use adaptive quadrature: negative-drift Wald has heavy tails that defeat fixed grids.
  int_swtn <- integrate(
    function(t) vapply(t, function(ti) EMC2:::dswtn(ti, v, b, sv = sv, s = s, posdrift = TRUE), numeric(1)),
    lower = 1e-6, upper = Inf, subdivisions = 500L, rel.tol = 1e-4
  )
  expect_equal(int_swtn$value, 1.0, tolerance = 0.01)

  int_rdm <- integrate(
    function(t) vapply(t, function(ti) EMC2:::drdmswtn(ti, v, b, A, sv = sv, s = s, posdrift = TRUE), numeric(1)),
    lower = 1e-6, upper = Inf, subdivisions = 300L, rel.tol = 1e-3
  )
  expect_equal(int_rdm$value, 1.0, tolerance = 0.02)
})

test_that("dswtn posdrift=FALSE integrates to hit mass", {
  tt <- seq(0.02, 100, by = 0.05)
  v <- -0.4; b <- 1.0; sv <- 0.5; s <- 1.0

  d_def  <- sapply(tt, function(t) EMC2:::dswtn(t, v, b, sv = sv, s = s, posdrift = FALSE))
  p_inf  <- EMC2:::pswtn(Inf, v, b, sv = sv, s = s, posdrift = FALSE)
  expect_equal(sum(d_def) * 0.05, p_inf, tolerance = 0.02)
})

test_that("killed swtn cdf is locally consistent with the pdf", {
  t <- 0.9
  h <- 1e-6
  b <- 1.3
  v <- 1.0
  sv <- 0.45
  s <- 1.0
  lambda <- 0.6

  pdf <- EMC2:::dswtn(t, v, b, sv = sv, s = s, lambda_k = lambda)
  slope <- (EMC2:::pswtn(t + h, v, b, sv = sv, s = s, lambda_k = lambda) -
            EMC2:::pswtn(t - h, v, b, sv = sv, s = s, lambda_k = lambda)) / (2 * h)

  expect_equal(slope, pdf, tolerance = 2e-4)
})

test_that("wald small-k limit is numerically stable (q near 0 regime)", {
  t <- 0.95
  b <- 1.35
  mu <- 1.05
  sigma <- 1.0
  A <- 0.25
  t0 <- 0.2
  k_small <- 1e-10

  cdf0 <- EMC2:::pwald(t, mu, b, sigma = sigma, A = A, t0 = t0, lambda_k = 0.0, guess = FALSE)
  cdfk <- EMC2:::pwald(t, mu, b, sigma = sigma, A = A, t0 = t0, lambda_k = k_small, guess = FALSE)
  pdf0 <- EMC2:::dwald(t, mu, b, sigma = sigma, A = A, t0 = t0, lambda_k = 0.0, guess = FALSE)
  pdfk <- EMC2:::dwald(t, mu, b, sigma = sigma, A = A, t0 = t0, lambda_k = k_small, guess = FALSE)

  expect_equal(cdfk, cdf0, tolerance = 1e-8)
  expect_equal(pdfk, pdf0, tolerance = 1e-7)
})

test_that("joint A and sv RDMSWTN density matches explicit GL over thresholds", {
  skip_if_not_installed("statmod")
  gl <- statmod::gauss.quad(80, kind = "legendre")
  ref_density <- function(rt, v, b, A, s, t0, sv) {
    center <- b - 0.5 * A
    half_width <- 0.5 * A
    thresholds <- center + half_width * gl$nodes
    0.5 * sum(gl$weights * vapply(thresholds, function(threshold) {
      EMC2:::dswtn(
        t = rt, mu_drift = v, threshold = threshold, s = s, t0 = t0,
        sv = sv, lambda_g = 0, lambda_k = 0, log_out = FALSE, posdrift = FALSE
      )
    }, numeric(1)))
  }

  cases <- data.frame(
    rt = c(0.55, 0.95, 1.7, 3.0),
    v  = c(1.1,  0.8,  1.6, 0.35),
    b  = c(1.4,  1.2,  1.8, 1.1),
    A  = c(0.25, 0.45, 0.65, 0.2),
    s  = c(1.0,  0.8,  1.2, 0.9),
    t0 = c(0.15, 0.2,  0.1, 0.25),
    sv = c(0.35, 0.6,  0.25, 0.9)
  )

  for (i in seq_len(nrow(cases))) {
    pars <- cases[i, ]
    closed <- EMC2:::drdmswtn(
      t = pars$rt, mu_drift = pars$v, b = pars$b, A = pars$A, s = pars$s,
      t0 = pars$t0, sv = pars$sv, lambda_g = 0, lambda_k = 0,
      n_gauss_nodes = 20, log_out = FALSE, posdrift = FALSE
    )
    gl_ref <- ref_density(pars$rt, pars$v, pars$b, pars$A, pars$s,
                          pars$t0, pars$sv)
    expect_equal(closed, gl_ref, tolerance = 1e-6,
                 label = paste("row", i))
  }
})

test_that("joint A and sv RDMSWTN log density is consistent with density scale", {
  args <- list(t = 1.1, mu_drift = 0.9, b = 1.35, A = 0.3, s = 0.95,
               t0 = 0.2, sv = 0.5, lambda_g = 0, lambda_k = 0,
               n_gauss_nodes = 20)
  dens     <- do.call(EMC2:::drdmswtn, c(args, list(log_out = FALSE)))
  log_dens <- do.call(EMC2:::drdmswtn, c(args, list(log_out = TRUE)))

  expect_equal(log_dens, log(dens), tolerance = 1e-12)
})

test_that("combined local guess+kill exponential cdf matches independent pdf integration", {
  rt <- 0.95
  t0 <- 0.2
  lambda_g <- 0.45
  lambda_k <- 0.8

  swtn_cdf <- EMC2:::pSWTNspv(
    t = rt, v = 1.1, b = 1.35, A = 0.25, t0 = t0, sv = 0.0, s = 1.0,
    lambda_g = lambda_g, lambda_k = lambda_k, kill_shape = 1L
  )
  swtn_int <- integrate(function(x) {
    EMC2:::dSWTNspv(
      t = x, v = 1.1, b = 1.35, A = 0.25, t0 = t0, sv = 0.0, s = 1.0,
      lambda_g = lambda_g, lambda_k = lambda_k, kill_shape = 1L
    )
  }, lower = 0, upper = rt, subdivisions = 200L, rel.tol = 1e-10)$value

  gbm_cdf <- EMC2:::pGBMspv(
    t = rt, v = 1.0, b = 1.55, A = 0.25, t0 = t0, s = 0.9,
    lambda_g = lambda_g, lambda_k = lambda_k, kill_shape = 1L
  )
  gbm_int <- integrate(function(x) {
    EMC2:::dGBMspv(
      t = x, v = 1.0, b = 1.55, A = 0.25, t0 = t0, s = 0.9,
      lambda_g = lambda_g, lambda_k = lambda_k, kill_shape = 1L
    )
  }, lower = 0, upper = rt, subdivisions = 200L, rel.tol = 1e-10)$value

  expect_equal(swtn_cdf, swtn_int, tolerance = 1e-7)
  expect_equal(gbm_cdf, gbm_int, tolerance = 1e-7)
})

test_that("combined local guess+kill exponential cdf preserves pre-t0 guess domain", {
  rt <- 0.08
  t0 <- 0.3
  lambda_g <- 0.6
  lambda_k <- 0.9
  rate <- lambda_g + lambda_k
  expected <- lambda_g / rate * (1 - exp(-rate * rt))

  expect_equal(
    EMC2:::pSWTNspv(
      t = rt, v = 1.1, b = 1.35, A = 0.2, t0 = t0, sv = 0.1, s = 1.0,
      lambda_g = lambda_g, lambda_k = lambda_k, kill_shape = 1L
    ),
    expected,
    tolerance = 1e-12
  )
  expect_equal(
    EMC2:::pGBMspv(
      t = rt, v = 1.0, b = 1.5, A = 0.2, t0 = t0, s = 1.0,
      lambda_g = lambda_g, lambda_k = lambda_k, kill_shape = 1L
    ),
    expected,
    tolerance = 1e-12
  )
})

test_that("combined local guess+kill Erlang-2 point-start cdf is closed form for Wald and GBM", {
  rt <- 0.95
  t0 <- 0.2
  lambda_g <- 0.45
  lambda_k <- 0.8

  wald_cdf <- EMC2:::pSWTNspv(
    t = rt, v = 1.1, b = 1.35, A = 0, t0 = t0, sv = 0, s = 1.0,
    lambda_g = lambda_g, lambda_k = lambda_k, kill_shape = 2L
  )
  wald_int <- integrate(function(x) {
    EMC2:::dSWTNspv(
      t = x, v = 1.1, b = 1.35, A = 0, t0 = t0, sv = 0, s = 1.0,
      lambda_g = lambda_g, lambda_k = lambda_k, kill_shape = 2L
    )
  }, lower = 0, upper = rt, subdivisions = 500L, rel.tol = 1e-10)$value

  gbm_cdf <- EMC2:::pGBMspv(
    t = rt, v = 1.0, b = 1.55, A = 0, t0 = t0, s = 0.9,
    lambda_g = lambda_g, lambda_k = lambda_k, kill_shape = 2L
  )
  gbm_int <- integrate(function(x) {
    EMC2:::dGBMspv(
      t = x, v = 1.0, b = 1.55, A = 0, t0 = t0, s = 0.9,
      lambda_g = lambda_g, lambda_k = lambda_k, kill_shape = 2L
    )
  }, lower = 0, upper = rt, subdivisions = 500L, rel.tol = 1e-10)$value

  expect_equal(wald_cdf, wald_int, tolerance = 1e-10)
  expect_equal(gbm_cdf, gbm_int, tolerance = 1e-10)
})

test_that("combined local guess+kill Erlang-2 SPV cdf is closed form for Wald and GBM", {
  rt <- 0.95
  t0 <- 0.2
  lambda_g <- 0.45
  lambda_k <- 0.8

  wald_cdf <- EMC2:::pSWTNspv(
    t = rt, v = 1.1, b = 1.35, A = 0.25, t0 = t0, sv = 0, s = 1.0,
    lambda_g = lambda_g, lambda_k = lambda_k, kill_shape = 2L
  )
  wald_int <- integrate(function(x) {
    EMC2:::dSWTNspv(
      t = x, v = 1.1, b = 1.35, A = 0.25, t0 = t0, sv = 0, s = 1.0,
      lambda_g = lambda_g, lambda_k = lambda_k, kill_shape = 2L
    )
  }, lower = 0, upper = rt, subdivisions = 800L, rel.tol = 1e-10)$value

  gbm_cdf <- EMC2:::pGBMspv(
    t = rt, v = 1.0, b = 1.55, A = 0.25, t0 = t0, s = 0.9,
    lambda_g = lambda_g, lambda_k = lambda_k, kill_shape = 2L
  )
  gbm_int <- integrate(function(x) {
    EMC2:::dGBMspv(
      t = x, v = 1.0, b = 1.55, A = 0.25, t0 = t0, s = 0.9,
      lambda_g = lambda_g, lambda_k = lambda_k, kill_shape = 2L
    )
  }, lower = 0, upper = rt, subdivisions = 800L, rel.tol = 1e-10)$value

  expect_equal(wald_cdf, wald_int, tolerance = 1e-10)
  expect_equal(gbm_cdf, gbm_int, tolerance = 1e-10)
})

test_that("combined local guess+kill SWTN Erlang-2 cdf is locally consistent with the pdf", {
  t <- 0.85
  h <- 1e-5
  b <- 1.25
  v <- 1.05
  sv <- 0.35
  s <- 1.0
  lambda_g <- 0.4
  lambda_k <- 0.7

  pdf <- EMC2:::dSWTNspv(
    t = t, v = v, b = b, A = 0, t0 = 0, sv = sv, s = s,
    lambda_g = lambda_g, lambda_k = lambda_k, kill_shape = 2L
  )
  slope <- (
    EMC2:::pSWTNspv(
      t = t + h, v = v, b = b, A = 0, t0 = 0, sv = sv, s = s,
      lambda_g = lambda_g, lambda_k = lambda_k, kill_shape = 2L
    ) -
    EMC2:::pSWTNspv(
      t = t - h, v = v, b = b, A = 0, t0 = 0, sv = sv, s = s,
      lambda_g = lambda_g, lambda_k = lambda_k, kill_shape = 2L
    )
  ) / (2 * h)

  expect_equal(slope, pdf, tolerance = 5e-4)
})

test_that("combined local guess+kill Erlang-2 respects posdrift in the sv=0 analytic CDF path", {
  rt <- 0.9
  t0 <- 0.15
  lambda_g <- 0.45
  lambda_k <- 0.7

  cdf <- EMC2:::pSWTNspv(
    t = rt, v = -0.2, b = 1.3, A = 0.25, t0 = t0, sv = 0, s = 1.0,
    lambda_g = lambda_g, lambda_k = lambda_k, kill_shape = 2L,
    posdrift = TRUE
  )
  pdf_int <- integrate(function(x) {
    EMC2:::dSWTNspv(
      t = x, v = -0.2, b = 1.3, A = 0.25, t0 = t0, sv = 0, s = 1.0,
      lambda_g = lambda_g, lambda_k = lambda_k, kill_shape = 2L,
      posdrift = TRUE
    )
  }, lower = 0, upper = rt, subdivisions = 800L, rel.tol = 1e-10)$value

  expect_equal(cdf, pdf_int, tolerance = 1e-8)
})

test_that("rRDMSWTN local kill+guess yields both omissions and responses", {
  set.seed(11)
  n_trials <- 1000
  lR <- factor(rep(c("r1", "r2"), n_trials), levels = c("r1", "r2"))
  trial_rows <- rbind(
    c(v = 1e-6, b = 1.1, A = 0.1, t0 = 0.2, sv = 1e-8, lambda_g = 0.35, lambda_k = 1.0),
    c(v = 1e-6, b = 1.1, A = 0.1, t0 = 0.2, sv = 1e-8, lambda_g = 0.0,  lambda_k = 1.0)
  )
  pars <- trial_rows[rep(seq_len(nrow(trial_rows)), n_trials), , drop = FALSE]

  sim <- EMC2:::rRDMSWTN(lR, pars, erlang_shape = 2L, erlang_type = "local_kill_guess")

  expect_gt(mean(is.na(sim$R)), 0.05)
  expect_gt(mean(!is.na(sim$R)), 0.05)
})

test_that("rRDMSWTN returns omissions when killing wins the race", {
  set.seed(1)

  n_trials <- 200
  lR <- factor(rep(c("r1", "r2"), n_trials), levels = c("r1", "r2"))
  pars <- matrix(
    rep(c(v = 1.0, b = 1.1, A = 0.1, t0 = 0.2, sv = 0.25, lambda_g = 0, lambda_k = 1e6), 2 * n_trials),
    ncol = 7,
    byrow = TRUE,
    dimnames = list(NULL, c("v", "b", "A", "t0", "sv", "lambda_g", "lambda_k"))
  )

  sim <- EMC2:::rRDMSWTN(lR, pars, erlang_type = "local_kill")

  expect_true(all(is.na(sim$R)))
  expect_true(all(is.infinite(sim$rt)))
})

test_that("rRDMSWTN nests to the no-kill model when lambda_k is zero", {
  set.seed(1)

  designRDMSWTN <- design(
    factors = list(S = "Target", subjects = 1, L = c("L", "M", "H")),
    Rlevels = c("Go"),
    formula = list(v ~ L, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1, mG ~ 1, mK ~ 1),
    constants = c(s = log(1)),
    model = RDMSWTN(erlang_type = "local_kill"),
    report_p_vector = FALSE
  )

  p_vec <- sampled_pars(designRDMSWTN)
  p_vec[1:3] <- log(c(3, 2, 1.5))
  p_vec["B"] <- log(1.1)
  p_vec["A"] <- log(0.4)
  p_vec["t0"] <- log(0.2)
  p_vec["sv"] <- log(0.25)
  p_vec["mG"] <- log(1e10)
  p_vec["mK"] <- log(1e10)

  expect_silent(dat <- make_data(p_vec, design = designRDMSWTN, n_trials = 200))
  expect_true(is.data.frame(dat))
  expect_true(all(is.finite(dat$rt)))
  expect_false(anyNA(dat$R))
})

test_that("RDMSWTN requires lambda_g/lambda_k in model-level parameter matrices", {
  rt <- c(0.5, 0.8)
  pars_lambda <- cbind(v = 1.1, b = 1.4, A = 0.2, t0 = 0.1, sv = 0.3, s = 1.0,
                       lambda_g = 0.0, lambda_k = 0.4)

  d1 <- EMC2:::dRDMSWTN(rt, pars_lambda)
  p1 <- EMC2:::pRDMSWTN(rt, pars_lambda)

  expect_true(all(is.finite(d1)))
  expect_true(all(is.finite(p1)))
  expect_error(EMC2:::dRDMSWTN(rt, cbind(v = 1.1, b = 1.4, A = 0.2, t0 = 0.1, sv = 0.3, s = 1.0, k = 0.4)), "lambda_g")
  expect_error(EMC2:::pRDMSWTN(rt, cbind(v = 1.1, b = 1.4, A = 0.2, t0 = 0.1, sv = 0.3, s = 1.0, k = 0.4)), "lambda_g")
})

test_that("RDMSWTN and RDMGBM wrappers preserve rt<t0 Erlang guess mass", {
  rt <- c(0.05, 0.12, 0.20)
  t0 <- 0.35
  ks <- 2L
  lg <- 0.8

  ref_pdf <- dgamma(rt, shape = ks, rate = lg)
  ref_cdf <- pgamma(rt, shape = ks, rate = lg)

  pars_swtn <- cbind(v = 1.0, b = 1.2, A = 0.1, t0 = t0, sv = 0.2, s = 1.0,
                     lambda_g = lg, lambda_k = 0.0)
  pars_gbm  <- cbind(v = 1.0, b = 1.2, A = 0.1, t0 = t0, s = 1.0,
                     lambda_g = lg, lambda_k = 0.0)

  expect_equal(EMC2:::dRDMSWTN(rt, pars_swtn, erlang = ks), ref_pdf, tolerance = 5e-6)
  expect_equal(EMC2:::pRDMSWTN(rt, pars_swtn, erlang = ks), ref_cdf, tolerance = 5e-6)
  expect_equal(EMC2:::dRDMGBM(rt, pars_gbm, erlang = ks), ref_pdf, tolerance = 5e-6)
  expect_equal(EMC2:::pRDMGBM(rt, pars_gbm, erlang = ks), ref_cdf, tolerance = 5e-6)
})

test_that("mixed RDMSWTN local Erlang mode weights Erlang-1 and Erlang-2 kernels", {
  rt <- 0.8
  args <- list(
    t = rt, v = 1.1, b = 1.35, A = 0.25, s = 1.0, t0 = 0.2, sv = 0.1,
    lambda_g = 0.45, lambda_k = 0.8, posdrift = TRUE
  )
  omega <- 0.35

  d1 <- do.call(EMC2:::dSWTNspv, c(args, list(kill_shape = 1L)))
  args_e2 <- args
  args_e2$lambda_g <- 2 * args_e2$lambda_g
  args_e2$lambda_k <- 2 * args_e2$lambda_k
  d2 <- do.call(EMC2:::dSWTNspv, c(args_e2, list(kill_shape = 2L)))
  p1 <- do.call(EMC2:::pSWTNspv, c(args, list(kill_shape = 1L)))
  p2 <- do.call(EMC2:::pSWTNspv, c(args_e2, list(kill_shape = 2L)))

  expect_equal(
    do.call(EMC2:::dSWTNspv, c(args, list(kill_shape = 3L, erlang_omega = omega))),
    omega * d1 + (1 - omega) * d2,
    tolerance = 1e-12
  )
  expect_equal(
    do.call(EMC2:::pSWTNspv, c(args, list(kill_shape = 3L, erlang_omega = omega))),
    omega * p1 + (1 - omega) * p2,
    tolerance = 1e-12
  )
})

test_that("RDMSWTN timer means are transformed to likelihood rates by Erlang shape", {
  pars_base <- cbind(v = 1, B = 1, A = 0.2, t0 = 0.2, s = 1, sv = 0.1,
                     mG = 2, mK = 4)
  pars_mixed <- cbind(pars_base, omega = 0.25)

  fixed_e1 <- RDMSWTN(erlang_shape = 1L, erlang_type = "local_kill_guess")$Ttransform(pars_base, NULL)
  fixed_e2 <- RDMSWTN(erlang_shape = 2L, erlang_type = "local_kill_guess")$Ttransform(pars_base, NULL)
  mixed    <- RDMSWTN(erlang_shape = "mixed", erlang_type = "local_kill_guess")$Ttransform(pars_mixed, NULL)

  expect_equal(unname(fixed_e1[, "lambda_g"]), unname(1 / pars_base[, "mG"]), tolerance = 1e-12)
  expect_equal(unname(fixed_e1[, "lambda_k"]), unname(1 / pars_base[, "mK"]), tolerance = 1e-12)
  expect_equal(unname(fixed_e2[, "lambda_g"]), unname(2 / pars_base[, "mG"]), tolerance = 1e-12)
  expect_equal(unname(fixed_e2[, "lambda_k"]), unname(2 / pars_base[, "mK"]), tolerance = 1e-12)
  expect_equal(unname(mixed[, "lambda_g"]), unname(1 / pars_base[, "mG"]), tolerance = 1e-12)
  expect_equal(unname(mixed[, "lambda_k"]), unname(1 / pars_base[, "mK"]), tolerance = 1e-12)
})

test_that("mixed RDMSWTN rfun supports local Erlang simulation only", {
  model <- RDMSWTN(erlang_shape = "mixed", erlang_type = "local_kill_guess")
  lR <- factor(rep(c("left", "right"), 20), levels = c("left", "right"))
  pars <- matrix(
    rep(c(v = 1.1, B = 1.0, A = 0.2, t0 = 0.1, s = 1.0, sv = 0.2,
          lambda_g = 0.8, lambda_k = 0.6, omega = 0.35, b = 1.2), 20),
    ncol = 10, byrow = TRUE,
    dimnames = list(NULL, c("v", "B", "A", "t0", "s", "sv", "lambda_g", "lambda_k", "omega", "b"))
  )

  expect_silent(sim <- model$rfun(list(lR = lR), pars))
  expect_s3_class(sim, "data.frame")
  expect_equal(nrow(sim), length(lR) / length(levels(lR)))
  expect_error(
    RDMSWTN(erlang_shape = "mixed", erlang_type = "global_kill"),
    "local Erlang processes only"
  )
})
