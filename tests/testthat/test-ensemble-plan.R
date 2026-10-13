# The between-block plan for the ensemble update (R/ensemble_plan.R): pool
# size and rounds from a block's record, on constructed records with known
# answers.

# Draws are iid (lag-1 displacement 2 per coordinate) except the first group
# mean, an AR(1) with correlation `rho1`. `J` is the per-round displacement on
# every coordinate, or a matrix with one row per coordinate (means, then log
# sds).
plan_sampler <- function(p, n, n_iter, M, R, J, e, c_round, c_member, C = 1, rho1 = 0) {
  set.seed(2)
  mu <- matrix(rnorm(p * n_iter), p, n_iter, dimnames = list(paste0("m", seq_len(p)), NULL))
  if (rho1 > 0) mu[1, ] <- as.numeric(stats::filter(sqrt(1 - rho1^2) * mu[1, ], rho1,
                                                    method = "recursive"))
  tv <- array(0, c(p, p, n_iter))
  for (t in seq_len(n_iter)) tv[, , t] <- diag(exp(2 * rnorm(p)), p)
  blk <- EMC2:::.emc_ensemble_block_new(p, n, M, R)
  blk$n <- n_iter
  v <- c(apply(mu, 1, var), apply(0.5 * log(apply(tv, 3, diag)), 1, var))
  if (!is.matrix(J)) J <- matrix(J, 2 * p, R, byrow = TRUE)
  # Per-round displacement sums that standardise to J.
  blk$jump <- v * n_iter * J
  blk$eff <- rep(e, n) * n_iter
  blk$moved <- rep(0.5, n) * n_iter
  blk$iter_wall <- C * n_iter
  blk$rounds_wall <- c_round * n_iter * R
  blk$pool_wall <- c_member * n_iter * (M - 1)
  list(samples = list(theta_mu = mu, theta_var = tv, idx = n_iter, stage = rep("sample", n_iter)),
       ensemble_stats = list(block = blk))
}

plan_of <- function(...) EMC2:::.emc_ensemble_plan_update(plan_sampler(...))$ensemble_plan

test_that("rounds follow the displacement the reselections bring, priced against a round", {
  s <- plan_sampler(2, 4, 400, 8L, 4L, J = c(1, 1.5, 1.8, 2.1), e = 0.5,
                    c_round = 0.01, c_member = 0.01)
  out <- EMC2:::.emc_ensemble_plan_update(s)
  rec <- out$ensemble_plan_history[[1]]
  expect_equal(rec$J_tot, 2, tolerance = 0.2)
  expect_equal(rec$G, 1.1, tolerance = 1e-6)
  # Each round adds 0.3 for 0.01 s against 2 per s: the rounds grow.
  expect_identical(out$ensemble_plan$rounds, 8L)
  expect_null(out$ensemble_stats$block)
  # The last round adds nothing: drop it.
  expect_identical(plan_of(2, 4, 400, 8L, 4L, J = c(1, 1.5, 1.6, 1.6), e = 0.5,
                           c_round = 0.01, c_member = 0.01)$rounds, 3L)
  # Nothing beyond the first Gibbs draw: two rounds and one fresh member.
  flat <- plan_of(2, 4, 400, 8L, 4L, J = c(1, 1, 1, 1), e = 0.5, c_round = 0.01, c_member = 0.01)
  expect_identical(c(flat$pool, flat$rounds), c(2L, 2L))
})

test_that("the pool follows the refreshed share's value, priced against a member", {
  J <- c(1, 1.5, 2, 2.5)
  cheap <- plan_of(2, 4, 400, 8L, 4L, J, e = 0.5, c_round = 0.001, c_member = 1e-4)$pool
  mid <- plan_of(2, 4, 400, 8L, 4L, J, e = 0.5, c_round = 0.001, c_member = 0.02)$pool
  dear <- plan_of(2, 4, 400, 8L, 4L, J, e = 0.5, c_round = 0.001, c_member = 0.1)$pool
  # Cheap members: the pool doubles; dearer ones: fewer.
  expect_identical(cheap, 15L)
  expect_true(cheap > mid && mid > dear)
  expect_lte(dear, 3L)
  # Without a measured price for the pool rows the pool is kept.
  expect_identical(plan_of(2, 4, 400, 8L, 4L, J, e = 0.5, c_round = 0.001, c_member = 0)$pool, 8L)
  # A record too short to price leaves the plan alone.
  expect_null(plan_of(2, 4, 1, 8L, 4L, J, e = 0.5, c_round = 0.01, c_member = 0.001))
})

test_that("the slowest coordinate decides the pool", {
  # The first group mean mixes slowly (lag-1 displacement about 0.2, against
  # about 2 elsewhere).
  slow_gain <- rbind(c(0, 0.05, 0.1, 0.15), matrix(0, 3, 4))
  fast_gain <- rbind(0, matrix(c(0, 0.3, 0.4, 0.5), 3, 4, byrow = TRUE))
  s <- plan_sampler(2, 4, 400, 8L, 4L, 1 + slow_gain, e = 0.5, c_round = 0.001,
                    c_member = 0.01, rho1 = 0.9)
  out <- EMC2:::.emc_ensemble_plan_update(s)
  rec <- out$ensemble_plan_history[[1]]
  expect_lt(rec$T_min, 0.4)
  expect_identical(rec$worst, "m1")
  # The reselections move the slow coordinate: the pool and the rounds stay
  # or grow.
  expect_gte(out$ensemble_plan$pool, 8L)
  expect_identical(out$ensemble_plan$rounds, 8L)
  # They move only coordinates that are already fast: the pool is cut,
  # although the gain averaged over coordinates is larger.
  fast <- plan_of(2, 4, 400, 8L, 4L, 1 + fast_gain, e = 0.5, c_round = 0.001,
                  c_member = 0.01, rho1 = 0.9)
  expect_identical(c(fast$pool, fast$rounds), c(2L, 2L))
})

test_that("pool rows are priced from their worker CPU over the chain's cores", {
  s <- plan_sampler(2, 4, 400, 8L, 4L, c(1, 1, 1, 1), e = 0.5, c_round = 0.01, c_member = 0)
  s <- EMC2:::.emc_ensemble_record_rows(s, list(cpu = 0.3, cores = 3))
  expect_equal(s$ensemble_stats$block$pool_wall, 0.1)
  s <- EMC2:::.emc_ensemble_record_rows(s, list(cpu = NULL, cores = 3))
  expect_equal(s$ensemble_stats$block$pool_wall, 0.1)
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
