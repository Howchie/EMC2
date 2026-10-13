test_that("marginal ensemble scores use the selected subjects and current group state", {
  alpha <- matrix(c(1, -0.2, 0.5, 2, 0.1, -1), 3, 2)
  loadings <- matrix(c(1, 0.7, -0.3, 0.2, 1, 0.4), 3, 2)
  mu <- c(0.1, -0.4, 0.2)
  precision <- c(2, 0.5, 3)
  sampler <- list(n_pars = 3L, n_subjects = 2L, n_factors = 2L)
  drawn_covariance <- NULL
  local_mocked_bindings(rmvnorm = function(n, sigma, ...) {
    drawn_covariance <<- sigma
    matrix(0, n, nrow(sigma))
  }, .package = "EMC2")
  for (type in c("factor", "infnt_factor")) {
    sampler$type <- type
    sampler$n_factors <- if (type == "factor") 2L else NULL
    pars <- list(tmu = mu, lambda_untransf = loadings, lambda = loadings * 2,
                 sig_err_inv = precision, epsilon_inv = precision,
                 psi_inv = c(0.7, 1.4), eta = matrix(99, 2, 2))
    L <- if (type == "factor") loadings else loadings * 2
    P <- if (type == "factor") diag(pars$psi_inv) else diag(2)
    V <- solve(P + t(L) %*% diag(precision) %*% L)
    expected <- t(V %*% t(L) %*% diag(precision) %*% (alpha - mu))
    out <- EMC2:::.emc_ensemble_refresh_scores(pars, sampler, alpha)
    expect_equal(out$eta, expected, tolerance = 1e-12)
    expect_equal(drawn_covariance, V, tolerance = 1e-12)
    expect_identical(out[names(out) != "eta"], pars[names(pars) != "eta"])
  }
})

test_that("SEM score restoration conditions on covariates and correlated disturbances", {
  sampler <- list(type = "SEM", n_pars = 3L, n_subjects = 2L, n_factors = 2L,
                  sem_settings = list(covariates = matrix(c(-0.5, 1), 2, 1)))
  pars <- list(tmu = c(0.1, -0.4, 0.2),
               lambda = matrix(c(1, 0.7, -0.3, 0.2, 1, 0.4), 3, 2),
               K = matrix(c(0.3, -0.4, 0.7), 3, 1), G = matrix(c(0.8, -0.2), 2, 1),
               B = matrix(c(0, 0.4, 0, 0), 2),
               delta_inv = matrix(c(1.5, 0.4, 0.4, 0.7), 2),
               epsilon_inv = c(2, 0.5, 3), eta = matrix(99, 2, 2))
  alpha <- matrix(c(1, -0.2, 0.5, 2, 0.1, -1), 3, 2)
  E <- diag(pars$epsilon_inv)
  B0 <- diag(2) - pars$B
  P <- t(B0) %*% pars$delta_inv %*% B0
  x <- t(sampler$sem_settings$covariates)
  m0 <- solve(B0, pars$G %*% x)
  V <- solve(P + t(pars$lambda) %*% E %*% pars$lambda)
  expected <- t(V %*% (P %*% m0 + t(pars$lambda) %*% E %*%
                         (alpha - pars$tmu - pars$K %*% x)))
  local_mocked_bindings(rmvnorm = function(n, sigma, ...) matrix(0, n, nrow(sigma)),
                       .package = "EMC2")
  out <- EMC2:::.emc_ensemble_refresh_scores(pars, sampler, alpha)
  expect_equal(out$eta, expected, tolerance = 1e-12)
})

test_that("score restoration preserves the full conditional covariance", {
  withr::local_seed(52)
  n <- 6000L
  sampler <- list(type = "factor", n_pars = 2L, n_subjects = n, n_factors = 2L)
  pars <- list(tmu = c(0.1, -0.2), lambda_untransf = matrix(c(1, 0.5, 0.4, 1.2), 2),
               sig_err_inv = c(2, 0.7), psi_inv = c(0.8, 1.4), eta = matrix(0, n, 2))
  alpha <- matrix(c(0.4, -0.1), 2, n)
  L <- pars$lambda_untransf
  V <- solve(diag(pars$psi_inv) + t(L) %*% diag(pars$sig_err_inv) %*% L)
  m <- drop(V %*% t(L) %*% diag(pars$sig_err_inv) %*% (alpha[, 1] - pars$tmu))
  scores <- EMC2:::.emc_ensemble_refresh_scores(pars, sampler, alpha)$eta
  expect_true(all(abs(colMeans(scores) - m) < 6 * sqrt(diag(V) / n)))
  variance_mcse <- sqrt((V^2 + outer(diag(V), diag(V))) / (n - 1))
  expect_true(all(abs(cov(scores) - V) < 6 * variance_mcse))
})
