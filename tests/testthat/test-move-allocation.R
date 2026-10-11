# Allocation of shear rounds and subject HMC (R/move_allocation.R).
local_rng_guard()

# rho of a unit under allocation x, by the composition the solver assumes.
alloc_unit_rho <- function(u, x) {
  switch(u$kind,
    fixed = max(u$rho),
    shear = u$rho * EMC2:::.emc_alloc_factor(u$g, x, if (is.null(u$a)) u$g else u$a),
    hmc = max(u$rho * (1 - x + x * u$h)))
}

test_that("the allocation is as good as a dense grid over rounds and rates", {
  local_test_rng(7)
  tau <- function(r) (1 + r) / (1 - r)
  for (trial in 1:30) {
    g <- stats::runif(1, 0.05, 0.95)
    # Half the trials saturate: a below g.
    a <- if (trial %% 2 == 0) stats::runif(1, 0, g) else g
    units <- list(
      list(kind = "shear", rho = stats::runif(1, 0, 0.99), g = g, a = a,
           cost = stats::runif(1, 0.02, 0.3)),
      list(kind = "hmc", rho = stats::runif(2, 0, 0.995), h = stats::runif(2, 0, 0.9),
           cost = stats::runif(1, 0.02, 0.5)),
      list(kind = "fixed", rho = stats::runif(1, 0, 0.9), cost = 0))
    C <- 1
    xs <- seq(0, 4, by = 0.02); rs <- seq(0, 1, by = 0.01)
    best <- Inf
    for (x in xs) for (r in rs) {
      # The fixed coordinate counts for nothing: no move lowers it.
      rho <- max(alloc_unit_rho(units[[1]], x), alloc_unit_rho(units[[2]], r))
      best <- min(best, tau(rho) * (C + x * units[[1]]$cost + r * units[[2]]$cost))
    }
    sol <- EMC2:::.emc_alloc_solve(units, C)
    rho_sol <- max(vapply(1:2, function(k) alloc_unit_rho(units[[k]], sol$x[k]), numeric(1)))
    expect_equal(sol$tau0_all, tau(max(vapply(units, alloc_unit_rho, 0, x = 0))))
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

test_that("a coordinate no move can lower does not bound the others", {
  solve <- EMC2:::.emc_alloc_solve
  # Every coordinate fixed: nothing to buy.
  expect_equal(solve(list(list(kind = "fixed", rho = 0.9, cost = 0)), 1)$x, 0)
  # A slow group coordinate no shear can reach does not stop the subjects'
  # HMC: the allocation is over what the moves can reach.
  sol <- solve(list(list(kind = "fixed", rho = 0.99, cost = 0),
                    list(kind = "hmc", rho = c(0.9, 0.5), h = c(0.1, 0.1), cost = 0.01)), 1)
  expect_equal(sol$x[1], 0)
  expect_gt(sol$x[2], 0)
  expect_equal(sol$tau0, EMC2:::.emc_alloc_tau(0.9))
  # A subject whose particle step and HMC do not move it (a drifting chain
  # reads every jump as nothing against the pooled variance) does not veto
  # the rounds that move the group.
  stuck <- list(kind = "subject", rho = c(0.9999, 0.999), p = c(1, 1), c_p = 0.1,
                h = c(1, 1), c_h = 0.1)
  sol <- solve(list(list(kind = "shear", rho = 0.95, g = 0.4, cost = 0.05), stuck), 1)
  expect_gt(sol$x[1], 0)
  expect_identical(c(sol$q[2], sol$x[2]), c(1, 0))
  expect_equal(sol$cost, 1 + sol$x[1] * 0.05 + 0.1)
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
  sol <- solve(list(list(kind = "shear", rho = 0.25, g = 0.8, cost = 0.2),
                    list(kind = "fixed", rho = 0.2, cost = 0)), 1)
  expect_equal(sol$x, c(0, 0))
})

test_that("expected rounds compose the measured round factor", {
  factor <- EMC2:::.emc_alloc_factor
  expect_equal(factor(c(0.1, 0.5, 0.9), 1), c(0.1, 0.5, 0.9))
  expect_equal(factor(0.9, 0), 1)
  expect_equal(factor(0.2, 2.5), 0.2^2 * (0.5 + 0.5 * 0.2))
  # Saturating moves: one move is g, many approach 1 - c, c = (1 - g) / (1 - a).
  expect_equal(factor(0.8, 1, 0.5), 0.8)
  expect_equal(factor(0.8, 2, 0.5), 1 - 0.4 * (1 - 0.25))
  expect_equal(factor(0.8, 200, 0.5), 0.6)
  expect_equal(factor(0.8, 3, 0), 0.8)
  # The inverse gives the least expected moves reaching a target, and Inf
  # past the floor.
  inv <- EMC2:::.emc_alloc_factor_inverse
  local_test_rng(4)
  for (trial in 1:200) {
    g <- stats::runif(1, 0.05, 0.98)
    a <- if (trial %% 2 == 0) stats::runif(1, 0, g) else g
    t <- stats::runif(1, 0.01, 0.99)
    x <- inv(g, a, t)
    floor_t <- if (a >= g) 0 else 1 - (1 - g) / (1 - a)
    if (t <= floor_t) {
      expect_identical(x, Inf)
    } else {
      expect_equal(factor(g, x, a), t, tolerance = 1e-9)
      if (x > 1e-6) expect_gt(factor(g, x - 1e-6, a), t)
    }
  }
})

test_that("the moves' correlation within their conditional is recovered from net jumps", {
  local_test_rng(9)
  # Moves as an autoregression with correlation a on a share c of a unit
  # variance; iterations run two or three of them.
  sim <- function(a, c, ks, n = 20000) {
    m1 <- 0; nm <- 0; net <- 0; K <- matrix(0L, 1, max(ks))
    for (i in seq_len(n)) {
      k <- ks[(i %% length(ks)) + 1L]
      x0 <- x <- stats::rnorm(1, 0, sqrt(c))
      for (j in seq_len(k)) {
        y <- a * x + sqrt(c * (1 - a^2)) * stats::rnorm(1)
        m1 <- m1 + (y - x)^2; nm <- nm + 1; x <- y
      }
      net <- net + (x - x0)^2
      K[1, k] <- K[1, k] + 1L
    }
    g <- 1 - (m1 / nm) / 2
    list(g = g, a = EMC2:::.emc_alloc_saturation(g, (net / n) / 2, K))
  }
  for (case in list(c(a = 0.3, c = 0.4), c(a = 0.7, c = 0.5), c(a = 0, c = 0.2))) {
    est <- sim(case[["a"]], case[["c"]], 2:3)
    expect_equal(est$g, 1 - case[["c"]] * (1 - case[["a"]]), tolerance = 0.02)
    expect_lt(abs(est$a - case[["a"]]), 0.06)
  }
  # The whole variance (c = 1): geometric.
  est <- sim(0.6, 1, 2:3)
  expect_lt(abs(est$a - est$g), 0.06)
  # One move per iteration does not identify a: geometric.
  expect_identical(EMC2:::.emc_alloc_saturation(0.7, 0.2, matrix(c(50L, 0L), 1)), 0.7)
  # Per-coordinate rows, and a coordinate without iterations.
  K <- rbind(c(0L, 10L), c(0L, 0L))
  expect_identical(EMC2:::.emc_alloc_saturation(c(0.8, 0.8), c(0.2, 0.2), K)[2], 0.8)
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
      best <- min(best, tau(max(rho_s)) * (C + q * subj$c_p + r * subj$c_h))
    }
    sol <- EMC2:::.emc_alloc_solve(units, C)
    q <- sol$q[1]; r <- sol$x[1]
    rho_s <- alloc_subject_rho(subj, q, r)
    # No coordinate ends slower than under its particle step alone, and the
    # solver's accounting is honest.
    expect_true(all(rho_s <= level + 1e-9))
    expect_lte(max(rho_s), sol$rho + 1e-9)
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
  expect_equal(sol$cost, 2)
  # Without HMC nothing can be lowered: the target is empty, the slowest
  # coordinate still reported.
  expect_equal(sol$tau0_all, EMC2:::.emc_alloc_tau(0.95 * 0.6))
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

test_that("under the refresh floor the cheaper kernel refreshes every subject", {
  local_test_rng(12)
  solve <- EMC2:::.emc_alloc_solve
  # A group coordinate the rounds barely move: a floor near 0.8.
  fixed <- list(kind = "shear", rho = 0.8, g = 0.999, cost = 0)
  # HMC moves one coordinate less than the particle step but costs half.
  subj <- list(kind = "subject", rho = c(0.9, 0.8), p = c(0.5, 0.3), c_p = 1,
               h = c(0.4, 0.5), c_h = 0.5)
  sol <- solve(list(fixed, subj), 1)
  expect_identical(c(sol$q[2], sol$x[2]), c(1, 0))
  withr::local_options(emc2.alloc_subject_floor = "refresh")
  sol <- solve(list(fixed, subj), 1)
  expect_identical(c(sol$q[2], sol$x[2]), c(0, 1))
  expect_equal(sol$cost, 1 + subj$c_h)
  # Every subject is refreshed (q + r >= 1) and meets the target; HMC that
  # costs more than the particle step never replaces it.
  dear <- replace(subj, "c_h", list(1.5))
  sol <- solve(list(fixed, dear), 1)
  expect_identical(c(sol$q[2], sol$x[2]), c(1, 0))
  for (trial in 1:20) {
    s <- list(kind = "subject", rho = stats::runif(2, 0.5, 1),
              p = stats::runif(2, 0.1, 0.9), c_p = stats::runif(1, 0.05, 0.5),
              h = stats::runif(2, 0, 0.9), c_h = stats::runif(1, 0.02, 0.5))
    units <- list(s, list(kind = "fixed", rho = stats::runif(1, 0, 0.9), cost = 0))
    sol <- solve(units, 1)
    q <- sol$q[1]; r <- sol$x[1]
    expect_gte(q + r, 1 - 1e-9)
    expect_true(all(alloc_subject_rho(s, q, r) <= sol$rho + 1e-9))
    expect_equal(sol$cost, 1 + q * s$c_p + r * s$c_h)
  }
})

test_that("particle steps are recorded at their marginal (worker CPU) cost", {
  p <- 2; n <- 3
  s <- list(n_pars = p, n_subjects = n)
  before <- matrix(0, p, n)
  after <- matrix(c(1, 1, 0, 0, 2, 0), p)
  # The round's wall cost (6) beyond the workers' CPU stays in the rest of
  # the iteration.
  s <- EMC2:::.emc_alloc_record_particle(s, before, after, c(TRUE, FALSE, TRUE), 6,
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

test_that("HMC transitions cost their worker CPU and the round a fixed cost", {
  s <- list(n_pars = 2, n_subjects = 3)
  tasks <- list(list(cpu = 1.5, work = 1), NULL, list(cpu = 2.5, work = 2))
  s <- EMC2:::.emc_alloc_record_cost(s, 0.5, 0, 0L, tasks, hmc_round = 7)
  s <- EMC2:::.emc_alloc_record_cost(s, 0.5, 0, 0L, list(NULL, NULL, NULL), hmc_round = 0)
  w <- s$alloc$window
  expect_equal(w$cpu_hmc, c(1, 0, 2))
  expect_identical(w$hmc_att, c(1L, 0L, 1L))
  expect_identical(w$hmc_rounds, 1L)
  expect_equal(c(w$hmc_round_cost, w$hmc_round_work), c(7, 3))
})

test_that("a plan pays the HMC round's fixed cost only if it uses HMC", {
  fixed <- list(kind = "fixed", rho = 0.2, cost = 0)
  subj <- list(kind = "subject", rho = c(0.95, 0.9), p = c(0.6, 0.5), c_p = 1,
               h = c(0.5, 0.5), c_h = 0.2)
  u <- list(units = list(fixed, subj, subj), C_rest = 1, joint = NULL)
  # Cheap and better per transition: HMC replaces the particle steps.
  sol <- EMC2:::.emc_alloc_solve_units(c(u, hmc_fixed = 0.1))
  expect_true(all(sol$x[2:3] > 0))
  expect_true(all(sol$q[2:3] < 1))
  # A round costing more than the particle steps it would save is not run.
  sol <- EMC2:::.emc_alloc_solve_units(c(u, hmc_fixed = 50))
  expect_equal(sol$x[2:3], c(0, 0))
  expect_equal(sol$q[2:3], c(1, 1))
  # Either way the reference is the iteration without moves or HMC.
  expect_equal(sol$objective0, EMC2:::.emc_alloc_tau(0.95 * 0.6) * 3)
})

test_that("production measurements use every production draw and pool the stage's jumps", {
  local_test_rng(5)
  p <- 1; n <- 1; burn <- 100; m <- 400; w_n <- 99L
  chain <- function(centre) {
    a <- c(stats::rnorm(burn, centre + 5), centre + stats::rnorm(m))
    k <- burn + m
    list(n_pars = p, n_subjects = n, type = "diagonal-gamma", par_names = "x",
         lik_laplace = list(NULL),
         samples = list(idx = k, stage = rep(c("burn", "sample"), c(burn, m)),
                        theta_mu = matrix(a, 1), theta_var = array(1, c(1, 1, k)),
                        alpha = array(a, c(1, 1, k))),
         alloc = list(stage = "sample", pool = list(G = diag(c(30, 0)), nG = 50L),
                      window = list(n = w_n, nG = w_n, nA = w_n,
                                    G = diag(c(20, 0)), A = list(matrix(20)),
                                    cpu_rest = 1, cpu_shear = 0, rounds = 0L,
                                    cpu_hmc = 0, hmc_att = 0L)))
  }
  emc <- list(chain(0), chain(0.5))
  prod <- EMC2:::.emc_alloc_measure(emc, 1, NULL, diag(1), "sample")
  x <- c(emc[[1]]$samples$theta_mu[burn + seq_len(m)],
         emc[[2]]$samples$theta_mu[burn + seq_len(m)])
  expect_equal(prod$g_var[1], stats::var(x))
  expect_identical(prod$pairs[["group"]], w_n + 50L)
  expect_equal(prod$pool$G, diag(c(50, 0)))
  expect_equal(prod$g_rho[1], EMC2:::.emc_alloc_rho_lo(50 / (w_n + 50), stats::var(x), w_n + 50))
  # The first production rebuild measures the window alone and starts the pool.
  first <- emc
  first[[1]]$alloc$stage <- "adapt"
  m1 <- EMC2:::.emc_alloc_measure(first, 1, NULL, diag(1), "sample")
  expect_identical(m1$pairs[["group"]], w_n)
  expect_null(m1$pool)
  # Outside production: the window's draws, the window's jumps.
  b <- EMC2:::.emc_alloc_measure(emc, 1, NULL, diag(1), "burn")
  k <- EMC2:::.emc_alloc_draws(emc[[1]], w_n)
  expect_equal(b$g_var[1], stats::var(c(emc[[1]]$samples$theta_mu[k],
                                        emc[[2]]$samples$theta_mu[k])))
  expect_identical(b$pairs[["group"]], w_n)
})

test_that("an unmeasured shear coordinate counts at its kind's measured factor", {
  p <- 2; n <- 1; k <- 50
  a <- matrix(stats::rnorm(p * k), p)
  ch <- list(n_pars = p, n_subjects = n, type = "diagonal-gamma", par_names = c("x", "y"),
             lik_laplace = list(NULL),
             samples = list(idx = k, theta_mu = a, theta_var = array(diag(2), c(p, p, k)),
                            alpha = array(a, c(p, 1, k))),
             shear = list(plan = list(dirs = diag(p),
                                      last = list(E = c(10, 0, 0, 0), n = c(10L, 0L, 0L, 0L),
                                                  moved = integer(4)))),
             alloc = list(window = NULL))
  withr::local_options(list(emc2.shear_move = TRUE))
  local_mocked_bindings(.emc_shear_active = function(sampler) TRUE, .package = "EMC2")
  ms <- EMC2:::.emc_alloc_measure(list(ch), 1, NULL, diag(p))
  g1 <- 1 - 1 / (2 * ms$g_var[1])
  expect_equal(ms$g_fac[1], g1)
  # The other mean direction borrows the measured one; the log sds, with no
  # measured round of their own, borrow it too.
  expect_equal(ms$g_fac[2:4], rep(g1, 3))
  expect_identical(ms$g_measured, c(TRUE, FALSE, FALSE, FALSE))
  # Without any measured round, the prior factor.
  ch$shear$plan$last$n[] <- 0L
  ms0 <- EMC2:::.emc_alloc_measure(list(ch), 1, NULL, diag(p))
  expect_equal(ms0$g_fac, rep(EMC2:::.EMC_ALLOC_SHEAR_PRIOR_G, 4))
})

test_that("an unmeasured shear coordinate counts at the prior factor, uncapped", {
  ms <- list(s_rho = matrix(0.5, 1, 1), g_rho = c(0.976, 0.6), g_fac = c(0.5, 0.9),
             g_measured = c(FALSE, TRUE), C_rest = 1, c_round = 0.05, coll = TRUE,
             c_fac = c(0.965, 0.5), c_measured = TRUE, c_coll = 0.1, c_rate = 1,
             W = list(diag(1)), s_var = matrix(1, 1, 1), c_hmc = NA_real_)
  u <- EMC2:::.emc_alloc_units(ms, NULL, TRUE, FALSE)
  expect_null(u$units[[1]]$cap)
  expect_identical(u$units[[1]]$g, 0.5)
  # The slowest coordinate's shear is unmeasured: the plan leans on its
  # prior factor rather than conceding to a floor at the guess.
  sol <- EMC2:::.emc_alloc_solve_units(u)
  expect_gt(sol$x[1], 0.25)
  expect_lt(sol$tau, EMC2:::.emc_alloc_tau(0.976 * EMC2:::.emc_alloc_factor(0.5, 0.25)))
})

test_that("collapsed round allocations remain specific to each chain", {
  local_mocked_bindings(
    .emc_alloc_wanted = function(sampler) TRUE,
    .emc_prepare_chain_rng = function(emc) emc,
    .emc_hmc_options = function() list(enabled = FALSE),
    .emc_alloc_draws = function(ch, n) 1L,
    .emc_alloc_measure = function(emc, j, ...) {
      list(chain_id = j, shear = FALSE, g_rho = numeric(), g_fac = numeric(),
           g_sat = numeric(), coll = TRUE, c_fac = 0.5, c_sat = 0.5,
           c_coll = 0.1, c_rate = 1, s_rho = matrix(0.5, 1L, 1L),
           s_var = matrix(1, 1L, 1L), p_fac = matrix(1, 1L, 1L),
           c_part = 0.1, c_hmc = NA_real_, stored = 0L, pairs = 0L,
           C_rest = 1, c_round = 0.1, pool = NULL, W = list(diag(1)))
    },
    .emc_alloc_units = function(ms, ...) {
      list(chain_id = ms$chain_id, units = list(list(h = NULL)), C_rest = 1,
           c_round = 0.1, joint = list(cost = 0.1), hmc_fixed = 0)
    },
    .emc_alloc_solve_units = function(u) {
      list(x = 0, joint = c(1.125, 0)[u$chain_id], tau = 1, tau0 = 1,
           cost = 1, objective = 1, objective0 = 1)
    },
    .package = "EMC2")

  chain <- function() list(
    n_pars = 1L, n_subjects = 1L,
    samples = list(theta_var = array(1, c(1L, 1L, 1L))),
    collapsed = list(plan = list(rounds = 0, eps = 0.4)),
    alloc = list(window = NULL, history = list()))
  out <- EMC2:::.emc_alloc_rebuild(list(chain(), chain()), "sample", 1L, 10L)
  expect_equal(vapply(out, function(ch) ch$collapsed$plan$rounds, numeric(1)),
               c(1.125, 0))
})
