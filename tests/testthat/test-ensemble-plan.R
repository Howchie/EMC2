# The between-block plan for the ensemble update (R/ensemble_plan.R): pool
# size and rounds from a block's record, on constructed records with known
# answers.

plan_sampler <- function(p, n, n_iter, M, R, J, e, c_round, c_member, C = 1) {
  set.seed(2)
  mu <- matrix(rnorm(p * n_iter), p, n_iter)
  tv <- array(0, c(p, p, n_iter))
  for (t in seq_len(n_iter)) tv[, , t] <- diag(exp(2 * rnorm(p)), p)
  blk <- EMC2:::.emc_ensemble_block_new(p, n, M, R)
  blk$n <- n_iter
  v <- c(apply(mu, 1, var), apply(0.5 * log(apply(tv, 3, diag)), 1, var))
  # Per-round displacement sums that standardise to J on every coordinate.
  blk$jump <- outer(v * n_iter, J)
  blk$eff <- rep(e, n) * n_iter
  blk$moved <- rep(0.5, n) * n_iter
  blk$iter_wall <- C * n_iter
  blk$rounds_wall <- c_round * n_iter * R
  blk$pool_wall <- c_member * n_iter * (M - 1)
  list(samples = list(theta_mu = mu, theta_var = tv, idx = n_iter, stage = rep("sample", n_iter)),
       ensemble_stats = list(block = blk))
}

test_that("rounds follow the increments the reselections bring, priced against a round", {
  # iid draws: the chain's own lag-1 displacement is 2 per coordinate, so a
  # unit of wall time must buy 2 (C = 1).
  s <- plan_sampler(2, 4, 400, 8L, 4L, J = c(1, 1.5, 1.6, 1.62), e = 0.5,
                    c_round = 0.01, c_member = 0.01)
  out <- EMC2:::.emc_ensemble_plan_update(s)
  rec <- out$ensemble_plan_history[[1]]
  expect_equal(rec$J_tot, 2, tolerance = 0.2)
  expect_equal(rec$G, 0.62, tolerance = 1e-6)
  # Every increment (0.5, 0.1, 0.02) per 0.01 s beats 2 per s: the rounds grow.
  expect_identical(out$ensemble_plan$rounds, 8L)
  expect_null(out$ensemble_stats$block)
  # The last increment below the rate: stop after the last that qualifies.
  s <- plan_sampler(2, 4, 400, 8L, 4L, J = c(1, 1.5, 1.6, 1.61), e = 0.5,
                    c_round = 0.01, c_member = 0.01)
  expect_identical(EMC2:::.emc_ensemble_plan_update(s)$ensemble_plan$rounds, 3L)
  # Nothing beyond the first Gibbs draw: the floor of two rounds.
  s <- plan_sampler(2, 4, 400, 8L, 4L, J = c(1, 1.001, 1.002, 1.003), e = 0.5,
                    c_round = 0.01, c_member = 0.01)
  expect_identical(EMC2:::.emc_ensemble_plan_update(s)$ensemble_plan$rounds, 2L)
})

test_that("the pool follows the refreshed share's marginal value, priced against a member", {
  J <- c(1, 1.5, 1.6, 1.62)
  # D = G / m(7) with e = 0.5: 0.62 / 0.778; D m'(7) = 0.0197 per member.
  cheap <- plan_sampler(2, 4, 400, 8L, 4L, J, e = 0.5, c_round = 0.01, c_member = 0.001)
  expect_identical(EMC2:::.emc_ensemble_plan_update(cheap)$ensemble_plan$pool, 15L)
  mid <- plan_sampler(2, 4, 400, 8L, 4L, J, e = 0.5, c_round = 0.01, c_member = 0.02)
  expect_identical(EMC2:::.emc_ensemble_plan_update(mid)$ensemble_plan$pool, 5L)
  dear <- plan_sampler(2, 4, 400, 8L, 4L, J, e = 0.5, c_round = 0.01, c_member = 0.1)
  expect_identical(EMC2:::.emc_ensemble_plan_update(dear)$ensemble_plan$pool, 2L)
  # No displacement from the reselections: one fresh member.
  flat <- plan_sampler(2, 4, 400, 8L, 4L, c(1, 1, 1, 1), e = 0.5, c_round = 0.01, c_member = 0.001)
  expect_identical(EMC2:::.emc_ensemble_plan_update(flat)$ensemble_plan$pool, 2L)
  # A record too short to price leaves the plan alone.
  short <- plan_sampler(2, 4, 1, 8L, 4L, J, e = 0.5, c_round = 0.01, c_member = 0.001)
  expect_null(EMC2:::.emc_ensemble_plan_update(short)$ensemble_plan)
})

test_that("the plan sets the options the move reads, and bounds them", {
  withr::local_options(list(emc2.ensemble_pool = 32L, emc2.ensemble_rounds = 10L))
  cfg <- EMC2:::.emc_ensemble_options(list(ensemble_plan = list(pool = 6L, rounds = 3L)))
  expect_identical(c(cfg$pool, cfg$rounds), c(6L, 3L))
  cfg <- EMC2:::.emc_ensemble_options(list(ensemble_plan = list(pool = 1L, rounds = 0L)))
  expect_identical(c(cfg$pool, cfg$rounds), c(2L, 1L))
  cfg <- EMC2:::.emc_ensemble_options(NULL)
  expect_identical(c(cfg$pool, cfg$rounds), c(32L, 10L))
})
