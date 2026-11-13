mvtnorm_logd <- function(x, m, S) {
  R <- chol(S)
  u <- backsolve(R, t(x) - m, transpose = TRUE)
  -0.5 * colSums(u^2) - sum(log(diag(R))) - 0.5 * ncol(x) * log(2 * pi)
}

test_that("Gaussian KL matches the closed form and is zero at equality", {
  m0 <- c(0, 1); S0 <- diag(c(1, 4))
  m1 <- c(0.5, 1); S1 <- diag(c(2, 4))
  expected <- 0.5 * (2 / 1 + 0.25 / 1 - 1 + log(1 / 2))
  expect_equal(EMC2:::.emc_gauss_kl(m1, S1, m0, S0), expected)
  expect_equal(EMC2:::.emc_gauss_kl(m0, S0, m0, S0), 0)
  expect_true(is.na(EMC2:::.emc_gauss_kl(m0, -S0, m0, S0)))
})

check_key <- function(mode, cov, config = "c", subject = 1L)
  list(config = config, subject = subject, pooled = NULL,
       mode = mode, covariance = cov)

test_that("pooling needs the same configuration and a near Gaussian", {
  S <- diag(2)
  draws <- list(components = list(list(mean = c(0, 0), cov = S, n = 4L)),
                x = matrix(0, 4, 2), ll = rep(0, 4))
  entries <- list(list(key = check_key(c(0, 0), S), check = list(),
                       draws = draws))
  near <- EMC2:::.emc_laplace_check_pool(entries, check_key(c(0.01, 0), S))
  expect_true(near$stable)
  expect_equal(nrow(near$pool$x), 4L)
  far <- EMC2:::.emc_laplace_check_pool(entries, check_key(c(2, 0), S))
  expect_false(far$stable)
  expect_null(far$pool)
  other <- EMC2:::.emc_laplace_check_pool(entries,
                                          check_key(c(0, 0), S, config = "d"))
  expect_null(other$pool)
  expect_null(EMC2:::.emc_laplace_check_pool(list(), check_key(c(0, 0), S))$pool)
})

test_that("pooled check weights recover the target and cap the pool", {
  set.seed(1)
  p <- 2L
  prior <- list(mu = c(0, 0), var = diag(c(100, 100)))
  post_mu <- c(0.3, -0.2)
  post_S <- matrix(c(0.04, 0.01, 0.01, 0.02), 2)
  Pi <- solve(post_S)
  # Likelihood chosen so prior x likelihood is exactly N(post_mu, post_S).
  lik_prec <- Pi - solve(prior$var)
  lik_mu <- drop(solve(lik_prec, Pi %*% post_mu))
  ll_batch <- function(x) {
    d <- t(x) - lik_mu
    -0.5 * colSums(d * (lik_prec %*% d))
  }
  first <- EMC2:::.emc_laplace_is_start(post_mu, post_S, prior, ll_batch,
                                        keep_draws = TRUE)
  expect_lt(first$pareto_k, 0.5)
  expect_equal(nrow(first$draws$x), 256L)
  pool <- first$draws
  for (i in 1:6) {
    mode <- post_mu + c(0.01, 0) * i
    res <- EMC2:::.emc_laplace_is_start(mode, post_S, prior, ll_batch,
                                        pool = pool)
    expect_equal(res$status, "ok")
    expect_lt(res$pareto_k, 0.5)
    pool <- res$draws
    expect_lte(nrow(pool$x), EMC2:::.EMC_LAPLACE_POOL_MAX)
    expect_equal(sum(vapply(pool$components, `[[`, integer(1), "n")),
                 nrow(pool$x))
  }
  # Self-normalised weights against the mixture recover the target mean.
  n <- nrow(pool$x)
  log_q <- sapply(pool$components, function(cm)
    log(cm$n / n) + mvtnorm_logd(pool$x, cm$mean, cm$cov))
  log_q <- apply(log_q, 1L, function(r) max(r) + log(sum(exp(r - max(r)))))
  log_t <- mvtnorm_logd(pool$x, post_mu, post_S)
  w <- exp(log_t - log_q); w <- w / sum(w)
  expect_equal(drop(crossprod(w, pool$x)), post_mu, tolerance = 0.05)
})
