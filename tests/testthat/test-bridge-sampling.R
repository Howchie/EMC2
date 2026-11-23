# Log inverse-Wishart density from log determinants and an explicit inverse.
iwish_reference <- function(W, v, S) {
  p <- nrow(W)
  0.5 * v * as.numeric(determinant(S, logarithm = TRUE)$modulus) -
    0.5 * (v + p + 1) * as.numeric(determinant(W, logarithm = TRUE)$modulus) -
    0.5 * sum(diag(S %*% solve(W))) -
    sum(lgamma((v + 1 - seq_len(p)) / 2)) - 0.5 * v * p * log(2) -
    0.25 * p * (p - 1) * log(pi)
}

log_chol <- function(W) {
  L <- t(chol(W))
  diag(L) <- log(diag(L))
  L[lower.tri(L, diag = TRUE)]
}

test_that("inverse-Wishart log density is exact on well- and ill-conditioned matrices", {
  local_test_rng(1)
  for (p in c(1L, 2L, 5L, 9L, 13L)) {
    for (v in c(p - 0.5, p + 1, p + 10)) {
      for (condition in c(1, 1e6, 1e9, 1e12)) {
        correlation <- cov2cor(crossprod(matrix(rnorm(p * (p + 3)), p + 3)))
        sds <- exp(seq(0, log(condition) / 2, length.out = p))[sample(p)]
        W <- outer(sds, sds) * correlation
        scale <- exp(rnorm(p, 0, log(condition + 1) / 4))
        expected <- iwish_reference(W, v, diag(scale, p))
        expect_equal(EMC2:::log_diwish(W, v, diag(scale, p)), expected, tolerance = 1e-8)
        expect_equal(EMC2:::log_diwish_log_chol(log_chol(W), v, log(scale)), expected,
                     tolerance = 1e-8)
        expect_equal(EMC2:::log_diwish_log_chol(log_chol(solve(W)), v, log(scale),
                                                covariance = FALSE),
                     expected, tolerance = 1e-6)
      }
    }
  }
  # A full scale matrix.
  W <- matrix(c(2, 0.5, 0.5, 3), 2)
  S <- matrix(c(3, -1, -1, 2), 2)
  expect_equal(EMC2:::log_diwish(W, 4, S), iwish_reference(W, 4, S), tolerance = 1e-12)
})

test_that("inverse-Wishart log density is not repaired off the positive-definite cone", {
  expect_identical(EMC2:::log_diwish(matrix(c(1, 2, 2, 1), 2), 3, diag(2)), -Inf)
  expect_identical(EMC2:::log_diwish(diag(c(1, 0)), 3, diag(2)), -Inf)
  expect_identical(EMC2:::log_diwish(diag(2), 3, diag(c(1, -1))), -Inf)
  expect_identical(EMC2:::log_diwish_log_chol(c(0, NaN, 0), 3, c(0, 0)), -Inf)
  expect_identical(EMC2:::log_diwish_log_chol(c(-800, 0.5, 0), 3, c(0, 0)), -Inf)
  # Variances this far apart stay on the log scale.
  W <- diag(c(1, 1e-8))
  S <- diag(c(3, 1e-8))
  expect_equal(EMC2:::log_diwish(W, 4, S), iwish_reference(W, 4, S), tolerance = 1e-12)
  expect_equal(EMC2:::log_diwish_log_chol(c(300, 0, -300), 3, c(600, -600)),
               -1 - (lgamma(1.5) + lgamma(1) + 3 * log(2) + 0.5 * log(pi)),
               tolerance = 1e-12)
  expect_error(EMC2:::log_diwish(diag(2), 1, diag(2)), "degrees of freedom")
})

# Draws from the standard hierarchy's prior in the layout bridge sampling
# uses: subject parameters, then the group block.
standard_prior_draws <- function(n, par_group, n_subjects, v, A) {
  p <- length(par_group)
  theta_mu <- matrix(rnorm(p * n), p, n)
  a_half <- matrix(1 / rgamma(p * n, shape = 0.5, rate = 1 / A^2), p, n)
  theta_var <- array(0, c(p, p, n))
  alpha <- array(0, c(p, n_subjects, n))
  for (i in seq_len(n)) {
    for (group in unique(par_group)) {
      idx <- which(par_group == group)
      d <- length(idx)
      theta_var[idx, idx, i] <- if (d > 1) {
        solve(rWishart(1, v + d - 1, diag(a_half[idx, i] / (2 * v), d))[, , 1])
      } else {
        1 / rgamma(1, shape = v / 2, rate = v / a_half[idx, i])
      }
    }
    alpha[, , i] <- theta_mu[, i] +
      t(chol(theta_var[, , i])) %*% matrix(rnorm(p * n_subjects), p)
  }
  samples <- list(
    subjects = seq_len(n_subjects), par_names = paste0("p", seq_len(p)), n_pars = p,
    n_subjects = n_subjects, par_group = par_group, group_designs = NULL,
    is_blocked = duplicated(par_group) | duplicated(par_group, fromLast = TRUE),
    prior = list(theta_mu_mean = rep(0, p), theta_mu_var = diag(p), v = v, A = rep(A, p)),
    samples = list(theta_mu = theta_mu, a_half = a_half, theta_var = theta_var)
  )
  group <- EMC2:::bridge_add_group(matrix(numeric(), n, 0L), samples, seq_len(n), "standard")
  list(draws = cbind(t(matrix(alpha, p * n_subjects, n)), group),
       info = EMC2:::make_info(samples, "standard"),
       data = vector("list", n_subjects))
}

test_that("prior draws under a constant likelihood have log marginal likelihood zero", {
  local_test_rng(11)
  local_mocked_bindings(
    calc_ll_manager = function(proposals, dadm, model, ...) numeric(nrow(proposals)),
    .package = "EMC2"
  )
  # One covariance block with a wide spread of scales, then blocks of three
  # and two parameters beside an uncorrelated one.
  cases <- list(list(par_group = c(1, 1, 1, 1), v = 4, A = 0.3),
                list(par_group = c(1, 1, 1, 2, 2, 3), v = 8, A = 1))
  for (settings in cases) {
    case <- standard_prior_draws(30000L, settings$par_group, n_subjects = 5L,
                                 v = settings$v, A = settings$A)
    n <- nrow(case$draws)
    for (split in list(seq_len(n / 2), seq.int(n / 2 + 1, n))) {
      fit <- EMC2:::bridge_fit(case$draws, split)
      log_ml <- EMC2:::bridge_sampling(case$draws, case$data, case$info, n / 2, split, fit,
                                       1L, 1L, 5000L, 1e-5, 1e-10, 1e-6, FALSE)
      expect_lt(abs(log_ml), 0.5)
    }
  }
})

test_that("bridge updates retain mass when both sums underflow on the natural scale", {
  expect_equal(EMC2:::bridge_update(c(1000, 1001), c(-1002, -1003), log(0.5), log(0.5)),
               -2, tolerance = 1e-12)
  expect_equal(EMC2:::bridge_log_add(c(-Inf, -1000, NA_real_), c(-Inf, -1000, -1000)),
               c(-Inf, -1000 + log(2), NA_real_))
})

test_that("Bayes factor comparison reproduces the seeded LNR marginal deviance", {
  local_test_rng(12, kind = "Mersenne-Twister")
  result <- compare(list(m = samples_LNR), BayesFactor = TRUE, WAIC = FALSE,
                    both_splits = FALSE, cores_for_props = 1L, print_summary = FALSE)
  expect_equal(result$MD, -571.823715913, tolerance = 1e-8)
  expect_equal(result$wMD, 1)
})

test_that("proposal blocks preserve draws and densities at singleton boundaries", {
  local_test_rng(5)
  samples <- EMC2:::merge_chains(samples_LNR)
  info <- EMC2:::make_info(samples, samples$type)
  draws <- EMC2:::bridge_draws(samples_LNR, "sample", info)
  fit <- EMC2:::bridge_fit(draws, seq_len(70))
  # Only the last block contains a single draw.
  idx <- 71:77
  gen <- matrix(rnorm(length(idx) * ncol(draws)), length(idx))
  evaluate <- function(block_size) {
    EMC2:::eval.unnormalized.posterior(draws, idx, gen, samples$data, fit$mean,
                                     fit$chol, info, 1L, 1L, FALSE, block_size)
  }
  expect_equal(evaluate(3L), evaluate(length(idx)), tolerance = 1e-10)
})

test_that("two bridge estimates swap complementary halves and use all fitting draws", {
  local_test_rng()
  for (odd in c(FALSE, TRUE)) {
    emc <- samples_LNR
    if (odd) emc <- subset(emc, length.out = 49L)
    samples <- emc[[1L]]
    info <- EMC2:::make_info(samples, samples$type)
    draws <- EMC2:::bridge_draws(emc, "sample", info)
    n_eff <- round(ess_summary(emc, selection = "alpha", stat = "median", stat_only = TRUE) / 2)
    halves <- split(seq_len(nrow(draws)), seq_len(nrow(draws)) > round(nrow(draws) / 2))

    set.seed(12)
    expected <- vapply(halves, function(idx) {
      fitted <- draws[idx, , drop = FALSE]
      fit <- list(mean = colMeans(fitted),
                  chol = chol(Matrix::nearPD(var(fitted), base.matrix = TRUE)$mat))
      EMC2:::bridge_sampling(draws, samples$data, info, n_eff, idx, fit,
                             1L, 1L, 5000L, 1e-5, 1e-10, 1e-6, FALSE)
    }, numeric(1))

    set.seed(12)
    actual <- suppressWarnings(run_bridge_sampling(emc, both_splits = TRUE, cores_for_props = 1L))
    expect_equal(actual, mean(expected), tolerance = 1e-10)
  }
})
