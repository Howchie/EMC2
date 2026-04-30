test_that("BAwL density shrinks with higher lambda for finite RT", {
  rt <- c(0.45, 0.70, 1.10)
  t0 <- 0.15
  pars0 <- cbind(v = 1.2, sv = 0.4, A = 0.2, b = 1.3, t0 = t0, k = 0.25, lambda = 0.0)
  pars1 <- cbind(v = 1.2, sv = 0.4, A = 0.2, b = 1.3, t0 = t0, k = 0.25, lambda = 0.8)
  pars0 <- pars0[rep(1, length(rt)), , drop = FALSE]
  pars1 <- pars1[rep(1, length(rt)), , drop = FALSE]

  d0 <- EMC2:::dBAwL(rt, pars0)
  d1 <- EMC2:::dBAwL(rt, pars1)

  expect_true(all(d1 < d0))
})

test_that("BAwL k->0 branch stays finite and monotone as a CDF", {
  rt <- seq(0.25, 1.25, length.out = 10)
  pars <- cbind(v = 1.1, sv = 0.5, A = 0.25, b = 1.35, t0 = 0.1, k = 1e-12, lambda = 0.2)
  pars <- pars[rep(1, length(rt)), , drop = FALSE]
  p <- EMC2:::pBAwL(rt, pars)

  expect_true(all(is.finite(p)))
  expect_true(all(p >= 0))
  expect_true(all(p <= 1))
  expect_true(all(diff(p) >= -1e-8))
})

test_that("BAwL simulator yields more omissions for larger lambda", {
  set.seed(123)
  n_trials <- 600
  lR <- factor(rep(c("r1", "r2"), n_trials), levels = c("r1", "r2"))
  base <- c(v = 1.1, sv = 0.35, b = 1.2, A = 0.2, t0 = 0.2, k = 0.2)
  pars_lo <- matrix(rep(c(base, lambda = 0.05), 2 * n_trials), ncol = 7, byrow = TRUE,
                    dimnames = list(NULL, c("v", "sv", "b", "A", "t0", "k", "lambda")))
  pars_hi <- matrix(rep(c(base, lambda = 2.0), 2 * n_trials), ncol = 7, byrow = TRUE,
                    dimnames = list(NULL, c("v", "sv", "b", "A", "t0", "k", "lambda")))

  sim_lo <- EMC2:::rBAwL(lR, pars_lo)
  sim_hi <- EMC2:::rBAwL(lR, pars_hi)

  omit_lo <- mean(is.na(sim_lo$R))
  omit_hi <- mean(is.na(sim_hi$R))
  expect_gt(omit_hi, omit_lo)
})

test_that("BAwL Erlang-2 CDF matches direct numerical integration", {
  rt <- c(0.55, 0.9, 1.25)
  pars <- cbind(v = 1.15, sv = 0.45, A = 0.22, b = 1.28, t0 = 0.15, k = 0.3, lambda = 0.6)
  pars <- pars[rep(1, length(rt)), , drop = FALSE]

  approx <- EMC2:::pBAwL(rt, pars, erlang = 2L)
  direct <- vapply(seq_along(rt), function(i) {
    dt <- rt[i] - pars[i, "t0"]
    integrate(
      function(x) {
        EMC2:::dleakyba(
          t = x, A = pars[i, "A"], b = pars[i, "b"], v = pars[i, "v"],
          sv = pars[i, "sv"], k = pars[i, "k"], posdrift = TRUE
        ) * exp(-pars[i, "lambda"] * x) * (1 + pars[i, "lambda"] * x)
      },
      lower = 0, upper = dt, subdivisions = 400L, rel.tol = 1e-8
    )$value
  }, numeric(1))

  expect_lt(max(abs(approx - direct)), 2e-3)
})
