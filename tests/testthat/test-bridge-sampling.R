test_that("native inverse-Wishart retains covariance repair for ill-conditioned blocks", {
  reference <- function(W, v, S, repair = FALSE) {
    chol_W <- chol(if (repair) Matrix::nearPD(W, base.matrix = TRUE)$mat else W)
    chol_S <- chol(if (repair) Matrix::nearPD(S, base.matrix = TRUE)$mat else S)
    p <- nrow(S)
    log_denom <- sum(lgamma((v + 1 - seq_len(p)) / 2)) +
      0.5 * v * p * log(2) + 0.25 * p * (p - 1) * log(pi)
    exp(v * sum(log(diag(chol_S))) - (v + p + 1) * sum(log(diag(chol_W))) -
          0.5 * sum(S * chol2inv(chol_W)) - log_denom)
  }

  W <- matrix(c(2, 0.5, 0.5, 3), 2)
  S <- diag(c(3, 2))
  expect_equal(EMC2:::bridge_diwish(W, 4, S), reference(W, 4, S), tolerance = 1e-12)

  W <- diag(c(1, 1e-8))
  S <- diag(c(3, 1e-8))
  expect_true(is.na(EMC2:::bridge_diwish(W, 4, S)))
  expect_equal(EMC2:::robust_diwish(W, 4, S), reference(W, 4, S, repair = TRUE),
               tolerance = 1e-12)
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
