test_that("PSIS weights and Pareto k match loo::psis", {
  skip_if_not_installed("loo")
  set.seed(812)
  cases <- list(
    rnorm(256),
    rnorm(256, sd = 4),
    log(rexp(256)) * 3,
    c(rnorm(250), 40 + rnorm(6)),
    rnorm(37),
    seq(-2, 2, length.out = 100))
  for (lr in cases) {
    ours <- EMC2:::.emc_psis(lr)
    ref <- suppressWarnings(loo::psis(lr, r_eff = 1))
    expect_equal(ours$pareto_k, as.numeric(ref$diagnostics$pareto_k),
                 tolerance = 1e-12)
    expect_equal(ours$weights,
                 as.numeric(stats::weights(ref, log = FALSE, normalize = TRUE)),
                 tolerance = 1e-12)
  }
})

test_that("PSIS leaves the RNG alone and handles short or flat tails", {
  set.seed(9)
  before <- .Random.seed
  short <- EMC2:::.emc_psis(c(0.1, -0.3, 0.2))
  expect_identical(.Random.seed, before)
  expect_identical(short$pareto_k, Inf)
  expect_equal(sum(short$weights), 1)
  flat <- EMC2:::.emc_psis(c(rnorm(200), rep(5, 56)))
  expect_equal(sum(flat$weights), 1)
})
