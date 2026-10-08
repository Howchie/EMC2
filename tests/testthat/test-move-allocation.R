# Allocation of shear rounds and subject HMC (R/move_allocation.R).
local_rng_guard()

# rho of a unit under allocation x, by the composition the solver assumes.
alloc_unit_rho <- function(u, x) {
  switch(u$kind,
    fixed = max(u$rho),
    shear = u$rho * EMC2:::.emc_alloc_shear_factor(u$g, x),
    hmc = max(u$rho * (1 - x + x * u$h)))
}

test_that("the allocation is as good as a dense grid over rounds and rates", {
  local_test_rng(7)
  tau <- function(r) (1 + r) / (1 - r)
  for (trial in 1:30) {
    units <- list(
      list(kind = "shear", rho = stats::runif(1, 0, 0.99), g = stats::runif(1, 0.05, 0.95),
           cost = stats::runif(1, 0.02, 0.3)),
      list(kind = "hmc", rho = stats::runif(2, 0, 0.995), h = stats::runif(2, 0, 0.9),
           cost = stats::runif(1, 0.02, 0.5)),
      list(kind = "fixed", rho = stats::runif(1, 0, 0.9), cost = 0))
    C <- 1
    xs <- seq(0, 4, by = 0.02); rs <- seq(0, 1, by = 0.01)
    best <- Inf
    for (x in xs) for (r in rs) {
      rho <- max(alloc_unit_rho(units[[1]], x), alloc_unit_rho(units[[2]], r),
                 alloc_unit_rho(units[[3]], 0))
      best <- min(best, tau(rho) * (C + x * units[[1]]$cost + r * units[[2]]$cost))
    }
    sol <- EMC2:::.emc_alloc_solve(units, C)
    rho_sol <- max(vapply(seq_along(units), function(k)
      alloc_unit_rho(units[[k]], sol$x[k]), numeric(1)))
    # The solver's own accounting is honest.
    expect_lte(rho_sol, sol$rho + 1e-9)
    expect_equal(sol$cost, C + sum(sol$x * vapply(units, `[[`, 0, "cost")))
    if (best < 0.88 * sol$objective0) {
      expect_lte(sol$objective, best * 1.03)
    } else {
      expect_lte(sol$objective, sol$objective0 * (1 + 1e-12))
    }
    expect_true(sol$x[1] >= 0 && sol$x[1] <= EMC2:::.EMC_ALLOC_MAX_ROUNDS)
    expect_true(sol$x[2] >= 0 && sol$x[2] <= 1)
  }
})

test_that("nothing is spent below a floor the moves cannot reach", {
  solve <- EMC2:::.emc_alloc_solve
  # Every coordinate fixed: nothing to buy.
  expect_equal(solve(list(list(kind = "fixed", rho = 0.9, cost = 0)), 1)$x, 0)
  # A slow group coordinate no shear can reach bounds the subjects' HMC.
  sol <- solve(list(list(kind = "fixed", rho = 0.99, cost = 0),
                    list(kind = "hmc", rho = c(0.9, 0.5), h = c(0.1, 0.1), cost = 0.01)), 1)
  expect_equal(sol$x, c(0, 0))
  # A stuck subject with an effective, cheap HMC gets a rate.
  sol <- solve(list(list(kind = "fixed", rho = 0.3, cost = 0),
                    list(kind = "hmc", rho = c(0.995, 0.2), h = c(0.2, 0.5), cost = 0.05)), 1)
  expect_gt(sol$x[2], 0)
  # A shear round that does not move its coordinate is never bought.
  sol <- solve(list(list(kind = "shear", rho = 0.99, g = 1, cost = 0.01)), 1)
  expect_equal(sol$x, 0)
  # One that does is bought, within the per-coordinate cap.
  sol <- solve(list(list(kind = "shear", rho = 0.99, g = 0.4, cost = 0.05),
                    list(kind = "fixed", rho = 0.2, cost = 0)), 1)
  expect_gt(sol$x[1], 0)
  expect_lte(sol$x[1], EMC2:::.EMC_ALLOC_MAX_ROUNDS)
  # A small predicted gain is not worth moves.
  sol <- solve(list(list(kind = "shear", rho = 0.25, g = 0.3, cost = 0.2),
                    list(kind = "fixed", rho = 0.2, cost = 0)), 1)
  expect_equal(sol$x, c(0, 0))
})

test_that("expected rounds compose the measured round factor", {
  factor <- EMC2:::.emc_alloc_shear_factor
  expect_equal(factor(c(0.1, 0.5, 0.9), 1), c(0.1, 0.5, 0.9))
  expect_equal(factor(0.9, 0), 1)
  expect_equal(factor(0.2, 2.5), 0.2^2 * (0.5 + 0.5 * 0.2))
})

test_that("the rest of the iteration is recorded against the previous final state", {
  p <- 2; n <- 3
  s <- list(n_pars = p, n_subjects = n,
            samples = list(idx = 2L, theta_mu = cbind(c(0, 0), c(1, 2)),
                           theta_var = array(c(diag(2), diag(c(4, 9))), c(p, p, 2)),
                           alpha = array(c(rep(0, 6), 1:6), c(p, n, 2))))
  s <- EMC2:::.emc_alloc_record_group(s, list(tmu = c(3, 4), tvar = diag(c(1, 16))))
  d <- c(2, 2, -log(2), log(4 / 3))
  expect_equal(s$alloc$window$G, tcrossprod(d))
  s <- EMC2:::.emc_alloc_record_subjects(s, matrix(1:6, p) + 1)
  expect_equal(s$alloc$window$A[[2]], matrix(1, p, p))
  s <- EMC2:::.emc_alloc_record_cost(s, 0.5, 0.1, 2L, list(NULL, list(cpu = 0.2), NULL))
  w <- s$alloc$window
  expect_identical(c(w$n, w$nG, w$nA, w$rounds), c(1L, 1L, 1L, 2L))
  expect_equal(c(w$cpu_rest, w$cpu_shear), c(0.5, 0.1))
  expect_equal(w$cpu_hmc, c(0, 0.2, 0))
  expect_identical(w$hmc_att, c(0L, 1L, 0L))
})

test_that("a chain held in one region measures slow against the pooled spread", {
  local_test_rng(3)
  p <- 1; n <- 1; m <- 200
  chain <- function(centre, sd) {
    a <- centre + sd * stats::rnorm(m)
    list(n_pars = p, n_subjects = n, type = "diagonal-gamma", par_names = "x",
         lik_laplace = list(NULL),
         samples = list(idx = m, theta_mu = matrix(a, 1),
                        theta_var = array(1, c(1, 1, m)), alpha = array(a, c(1, 1, m))),
         alloc = list(window = list(n = m - 1L, nG = m - 1L, nA = m - 1L,
                                    G = diag(c(sum(diff(a)^2), 0)),
                                    A = list(matrix(sum(diff(a)^2))),
                                    cpu_rest = 1, cpu_shear = 0, rounds = 0L,
                                    cpu_hmc = 0, hmc_att = 0L)))
  }
  # Two chains mixing over the same spread look fast; with one of them held
  # at a distance, the held one's small jumps are slow against the pooled sd.
  mixed <- list(chain(0, 1), chain(0, 1))
  held <- list(chain(0, 1), chain(4, 0.05))
  rho <- function(emc, j) EMC2:::.emc_alloc_measure(emc, j, NULL, diag(1))$s_rho[1, 1]
  expect_lt(rho(mixed, 1), 0.2)
  expect_gt(rho(held, 2), 0.99)
  # The lower bound only removes noise: a stuck coordinate stays stuck.
  expect_equal(EMC2:::.emc_alloc_rho_lo(0, 1, 10), 1)
  expect_lt(EMC2:::.emc_alloc_rho_lo(2, 1, 100), 0.05)
})

# rho of a subject unit under particle rate q and HMC rate r.
alloc_subject_rho <- function(u, q, r) {
  base <- if (is.null(u$p)) u$rho else u$rho * (1 - q * (1 - u$p))
  base * if (is.null(u$h)) 1 else (1 - r * (1 - u$h))
}

test_that("particle and HMC rates are as good as a dense grid over both", {
  local_test_rng(11)
  tau <- function(r) (1 + r) / (1 - r)
  qs <- EMC2:::.EMC_ALLOC_PARTICLE_RATES
  for (trial in 1:20) {
    subj <- list(kind = "subject", rho = stats::runif(2, 0.5, 1),
                 p = stats::runif(2, 0.1, 0.9), c_p = stats::runif(1, 0.05, 0.5),
                 h = stats::runif(2, 0, 0.9), c_h = stats::runif(1, 0.02, 0.5))
    units <- list(subj, list(kind = "fixed", rho = stats::runif(1, 0, 0.9), cost = 0))
    C <- 1
    level <- subj$rho * subj$p
    best <- Inf
    for (q in qs) for (r in seq(0, 1, by = 0.005)) {
      rho_s <- alloc_subject_rho(subj, q, r)
      if (any(rho_s > level + 1e-12)) next
      rho <- max(rho_s, units[[2]]$rho)
      best <- min(best, tau(rho) * (C + q * subj$c_p + r * subj$c_h))
    }
    sol <- EMC2:::.emc_alloc_solve(units, C)
    q <- sol$q[1]; r <- sol$x[1]
    rho_s <- alloc_subject_rho(subj, q, r)
    # No coordinate ends slower than under its particle step alone, and the
    # solver's accounting is honest.
    expect_true(all(rho_s <= level + 1e-9))
    expect_lte(max(rho_s, units[[2]]$rho), sol$rho + 1e-9)
    expect_equal(sol$cost, C + q * subj$c_p + r * subj$c_h)
    expect_true(q %in% qs)
    if (best < 0.88 * sol$objective0) {
      expect_lte(sol$objective, best * 1.03)
    } else {
      expect_lte(sol$objective, sol$objective0 * (1 + 1e-12))
    }
  }
})

test_that("the particle step stands down only for a cheaper, better HMC rate", {
  solve <- EMC2:::.emc_alloc_solve
  fixed <- list(kind = "fixed", rho = 0.2, cost = 0)
  subj <- list(kind = "subject", rho = c(0.95, 0.9), p = c(0.6, 0.5), c_p = 1,
               h = c(0.1, 0.1), c_h = 0.2)
  # Its zero allocation runs every particle step and prices it.
  sol <- solve(list(fixed, replace(subj, "h", list(NULL))), 1)
  expect_identical(c(sol$q[2], sol$x[2]), c(1, 0))
  expect_equal(sol$objective0, EMC2:::.emc_alloc_tau(0.95 * 0.6) * 2)
  # HMC that moves the subject more for less replaces it.
  sol <- solve(list(fixed, subj), 1)
  expect_lt(sol$q[2], 1)
  expect_gt(sol$x[2], 0)
  # An HMC step that moves it less never replaces it.
  sol <- solve(list(fixed, replace(subj, "h", list(c(0.9, 0.9)))), 1)
  expect_identical(sol$q[2], 1)
  # Without a measured particle factor the step always runs, as before.
  sol <- solve(list(fixed, replace(subj, "p", list(NULL))), 1)
  expect_identical(sol$q[2], 1)
})

test_that("particle steps are recorded with their share of the cost", {
  p <- 2; n <- 3
  s <- list(n_pars = p, n_subjects = n)
  before <- matrix(0, p, n)
  after <- matrix(c(1, 1, 0, 0, 2, 0), p)
  s <- EMC2:::.emc_alloc_record_particle(s, before, after, c(TRUE, FALSE, TRUE), 3,
                                         c(1, 5, 2))
  w <- s$alloc$window
  expect_identical(w$nP, c(1L, 0L, 1L))
  expect_equal(w$P[[1]], matrix(1, p, p))
  expect_equal(w$P[[3]], diag(c(4, 0)))
  expect_equal(w$cpu_part, c(1, 0, 2))
  # Without worker times the cost is shared equally.
  s <- EMC2:::.emc_alloc_record_particle(s, before, after, c(TRUE, TRUE, FALSE), 2, NULL)
  expect_equal(s$alloc$window$cpu_part, c(2, 1, 2))
})
