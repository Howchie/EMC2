test_that("Gaussian screening ratios include the shear Jacobian and reverse exactly", {
  withr::local_seed(34)
  p <- 3L
  n <- 2L
  H <- crossprod(matrix(rnorm(p^2), p))
  fits <- list(list(prec = H, lin = rnorm(p)), NULL)
  a <- matrix(rnorm(p * n), p)
  b <- matrix(rnorm(p * n), p)
  log_lik <- function(x) sum(x[, 1] * fits[[1]]$lin) -
    drop(crossprod(x[, 1], H %*% x[, 1])) / 2
  delta <- EMC2:::.emc_shear_lik_delta(a, b, fits)
  expect_equal(delta, log_lik(b) - log_lik(a), tolerance = 1e-12)
  expect_equal(EMC2:::.emc_shear_lik_delta(b, a, fits), -delta, tolerance = 1e-12)
  expect_equal(EMC2:::.emc_shear_lik_delta(a, b, rep(list(NULL), n)), 0)
  fits[[1]]$lin[1] <- NA_real_
  expect_true(is.na(EMC2:::.emc_shear_lik_delta(a, b, fits)))
})

test_that("delayed acceptance screens likelihood calls and applies the exact correction", {
  cur <- list(mu = 0, V = matrix(1), lp = 0)
  prop <- list(state = list(mu = 1, V = matrix(1), lp = 0.2),
               alpha = matrix(1), log_jac = -0.4)
  x <- list(cur = cur, cur_hier = -0.3, alpha = matrix(0), ll = 0, rounds = 0L, wpool = NULL)
  sampler <- list(par_names = "x", lik_laplace = list(list(prec = matrix(2), lin = 0,
                 check = list(pareto_k = 0, outside = 0, laplace_status = "ok"))))
  calls <- 0L
  local_mocked_bindings(
    .emc_shear_propose = function(...) prop,
    .emc_shear_log_hier = function(...) -0.7,
    .emc_wpool_subject_ll = function(pool, rows, ...) {
      calls <<- calls + 1L
      list(ll = -2, pool = pool)
    }, .package = "EMC2")
  # Surrogate: -1; hierarchy/prior/Jacobian: -0.6; exact likelihood: -2.
  local_mocked_bindings(runif = function(...) 0.9, .package = "stats")
  ctx <- list(fits = sampler$lik_laplace, spec = NULL, delayed = TRUE)
  rejected <- EMC2:::.emc_shear_mh(x, c(1, 0), sampler, ctx)
  expect_equal(calls, 0L)
  expect_false(rejected$accepted)
  expect_identical(rejected$alpha, x$alpha)
  local_mocked_bindings(runif = function(...) 0, .package = "stats")
  accepted <- EMC2:::.emc_shear_mh(x, c(1, 0), sampler, ctx)
  expect_equal(calls, 1L)
  expect_true(accepted$accepted)
  expect_equal(accepted$accept_prob, exp(-1.6) * exp(-1), tolerance = 1e-12)
  expect_identical(accepted$ll, -2)
})

test_that("delayed acceptance obeys detailed balance even with a wrong Gaussian", {
  withr::local_seed(5)
  for (i in seq_len(1000L)) {
    target_ratio <- rnorm(1, 0, 5)
    surrogate_ratio <- rnorm(1, 0, 10)
    log_accept <- min(0, surrogate_ratio) + min(0, target_ratio - surrogate_ratio)
    log_reverse <- min(0, -surrogate_ratio) + min(0, -target_ratio + surrogate_ratio)
    expect_equal(log_accept - log_reverse, target_ratio, tolerance = 1e-12)
  }
})

test_that("a poor or unchecked surrogate cannot reject an exact-acceptable move", {
  x <- list(cur = list(mu = 0, V = matrix(1), lp = 0), cur_hier = 0,
            alpha = matrix(0), ll = 0, rounds = 0L, wpool = NULL)
  prop <- list(state = list(mu = 0, V = matrix(1), lp = 0),
               alpha = matrix(10), log_jac = 0)
  calls <- 0L
  local_mocked_bindings(
    .emc_shear_propose = function(...) prop,
    .emc_shear_log_hier = function(...) 0,
    .emc_wpool_subject_ll = function(pool, ...) {
      calls <<- calls + 1L
      list(ll = 0, pool = pool)
    }, .package = "EMC2")
  for (check in list(NULL, list(pareto_k = Inf, outside = 0, laplace_status = "ok"))) {
    sampler <- list(par_names = "x", n_pars = 1L, n_subjects = 1L, type = "diagonal-gamma",
                    prior = list(shape = 2, rate = 0.3),
                    lik_laplace = list(list(prec = matrix(2), lin = 0, check = check)))
    ctx <- EMC2:::.emc_shear_context(sampler)
    expect_false(ctx$delayed)
    out <- EMC2:::.emc_shear_mh(x, c(0, 0), sampler, ctx)
    expect_true(out$accepted)
    expect_identical(out$alpha, prop$alpha)
  }
  expect_equal(calls, 2L)
})
