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
  s <- EMC2:::.emc_alloc_record_cost(s, 0.5, list(NULL, list(cpu = 0.2), NULL),
                                     moves = list(shear = c(cpu = 0.1, own = 0.05, n = 2)))
  w <- s$alloc$window
  expect_identical(c(w$n, w$nG, w$nA, w$fam$shear$moves), c(1L, 1L, 1L, 2L))
  expect_equal(c(w$cpu_rest, w$fam$shear$cpu, w$fam$shear$own), c(0.5, 0.1, 0.05))
  expect_equal(w$cpu_hmc, c(0, 0.2, 0))
  expect_identical(w$hmc_att, c(0L, 1L, 0L))
})

test_that("the allocation preserves the between-chain rho floor", {
  rho <- c(0.2, 0.8)
  b <- c(0.5, 0.5)
  expect_equal(EMC2:::.emc_alloc_remove_factor(c(0.3, 0.6, 0.4), c(0.5, 0.5, 0)), c(0.6, 1, 1))
  unit <- list(kind = "shear", rho = rho, b_share = b, g = c(0.1, 0.2),
               cost = 0, cap = 100)
  predicted <- EMC2:::.emc_alloc_predict(list(unit),
    list(x = 100, q = 1, joint = NULL), max_rounds = 100)[[1]]
  expect_true(all(predicted >= b - 1e-12))
  expect_equal(predicted, b + (1 - b) * rho * unit$g^100, tolerance = 1e-12)

  hmc <- list(kind = "subject", rho = c(0.2, 0.8), b_share = b,
              p = c(1, 1), h = c(0, 0), c_p = 0, c_h = 0)
  need <- EMC2:::.emc_alloc_unit_need(hmc, rs = c(0.49, 0.5, 0.7), max_rounds = 4)
  expect_true(all(is.infinite(need$cost[1])))
  expect_equal(need$x[2:3], c(1, 0.5))
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
  # A stuck coordinate reads one; a jump of twice the variance reads zero.
  expect_equal(EMC2:::.emc_alloc_rho(0, 1), 1)
  expect_equal(EMC2:::.emc_alloc_rho(2, 1), 0)
  expect_equal(EMC2:::.emc_alloc_rho(0.5, 1), 0.75)
})

# rho of a subject unit under particle rate q and HMC rate r.
# A subject's rho under particle rate q and HMC rate r, the two kernels
# selected jointly as R/sampling.R runs them (one uniform per subject).
alloc_subject_rho <- function(u, q, r) {
  P <- if (is.null(u$p)) rep(1, length(u$rho)) else u$p
  H <- if (is.null(u$h)) rep(1, length(u$rho)) else u$h
  if (is.null(u$p)) q <- 1
  if (is.null(u$h)) r <- 0
  u$rho * vapply(seq_along(u$rho), function(i)
    EMC2:::.emc_alloc_dev_factor(q, r, P[i], H[i]), numeric(1))
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

test_that("a cheaper HMC step that moves a coordinate less keeps the particle step", {
  solve <- EMC2:::.emc_alloc_solve
  # A group coordinate the rounds barely move: a floor near 0.8.
  fixed <- list(kind = "shear", rho = 0.8, g = 0.999, cost = 0)
  # HMC moves one coordinate less than the particle step but costs half.
  subj <- list(kind = "subject", rho = c(0.9, 0.8), p = c(0.5, 0.3), c_p = 1,
               h = c(0.4, 0.5), c_h = 0.5)
  sol <- solve(list(fixed, subj), 1)
  expect_identical(c(sol$q[2], sol$x[2]), c(1, 0))
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
  s <- EMC2:::.emc_alloc_record_cost(s, 0.5, tasks, moves = list(hmc = 7))
  s <- EMC2:::.emc_alloc_record_cost(s, 0.5, list(NULL, NULL, NULL), moves = list(hmc = 0))
  w <- s$alloc$window
  expect_equal(w$cpu_hmc, c(1, 0, 2))
  expect_identical(w$hmc_att, c(1L, 0L, 1L))
  # Tasks without rows (frame failures) count none.
  expect_equal(w$rows_hmc, c(0, 0, 0))
  s2 <- EMC2:::.emc_alloc_record_cost(list(n_pars = 2, n_subjects = 3), 0.5,
    list(list(cpu = 1.5, work = 1, rows = 12), NULL, list(cpu = 2.5, work = 2, rows = 30)))
  expect_equal(s2$alloc$window$rows_hmc, c(12, 0, 30))
  expect_identical(w$fam$hmc$iters, 1L)
  expect_identical(w$fam$hmc$moves, 2L)
  expect_equal(c(w$fam$hmc$cpu, w$fam$hmc$own), c(7, 3))
  # Every family splits the same way: marginal per move, fixed per iteration
  # with moves (carried while the window is short).
  for (i in 1:9) s <- EMC2:::.emc_alloc_record_cost(s, 0.5, tasks, moves = list(hmc = 7))
  m <- EMC2:::.emc_alloc_marginal(s$alloc$window, "hmc")
  expect_equal(m$fixed, 4 / 1)
  expect_equal(m$c, 3 / 2)
  short <- EMC2:::.emc_alloc_marginal(w, "hmc", last = list(fixed = list(hmc = 0.25)))
  expect_equal(short$fixed, 0.25)
  expect_true(is.na(short$c))
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

test_that("particle rates held below one keep the plan with HMC", {
  fixed <- list(kind = "fixed", rho = 0.2, cost = 0)
  subj <- list(kind = "subject", rho = c(0.95, 0.9), p = c(0.6, 0.5), c_p = 1,
               h = c(0.5, 0.5), c_h = 0.2)
  u <- list(units = list(fixed, subj, subj), C_rest = 1, joint = NULL, hmc_fixed = 50)
  # Unheld, the dear round is not run and every particle step is kept.
  expect_identical(EMC2:::.emc_alloc_solve_units(u)$choice, "none")
  # Held at the rates shared across chains, the plan without HMC cannot run
  # them, so HMC makes up the subjects' movement.
  u$units[2:3] <- lapply(u$units[2:3], function(x) c(x, q_fixed = 0.5))
  sol <- EMC2:::.emc_alloc_solve_units(u)
  expect_false(identical(sol$choice, "none"))
  expect_identical(sol$options[["none"]], Inf)
  expect_equal(sol$q[2:3], c(0.5, 0.5))
  expect_true(all(sol$x[2:3] > 0))
})

test_that("every chain runs the HMC choice with the least summed objective", {
  opts <- list(c(none = 10, "0" = 8, "1" = 9.5),
               c(none = 10, "0" = 8.5, "1" = 8),
               c(none = 10, "0" = 9.8, "1" = 9.9))
  local_mocked_bindings(
    .emc_alloc_wanted = function(sampler) TRUE,
    .emc_prepare_chain_rng = function(emc) emc,
    .emc_hmc_support = function(sampler) list(ok = FALSE),
    .emc_alloc_draws = function(ch, n) 1L,
    .emc_alloc_measure = function(emc, j, ...) {
      list(chain_id = j, shear = FALSE, g_rho = numeric(), coll = FALSE,
           s_rho = matrix(0.5, 1L, 1L), p_fac = matrix(1, 1L, 1L), c_part = 0.1,
           c_hmc = NA_real_, stored = 0L, pairs = 0L, C_rest = 1, pool = NULL)
    },
    .emc_alloc_units = function(ms, ...) {
      list(chain_id = ms$chain_id, units = list(list(h = NULL)), C_rest = 1, hmc_fixed = 0)
    },
    .emc_alloc_solve_units = function(u) {
      o <- opts[[u$chain_id]]
      choice <- if (is.null(u$choice_fixed)) names(o)[which.min(o)] else u$choice_fixed
      list(x = 0, joint = 0, tau = 1, tau0 = 1, cost = 1, objective = o[[choice]],
           objective0 = 10, options = o, choice = choice)
    },
    .package = "EMC2")
  chain <- function() list(
    n_pars = 1L, n_subjects = 1L,
    samples = list(theta_var = array(1, c(1L, 1L, 1L))),
    alloc = list(window = NULL, history = list()))
  out <- EMC2:::.emc_alloc_rebuild(list(chain(), chain(), chain()), "sample", 1L, 10L)
  # Chain 1's own choice (0) wins after summing the options. Chain 3's
  # objective does not clear the local gain threshold, but it is finite and
  # the aggregate HMC objective clears the threshold.
  expect_identical(out[[1]]$alloc$history[[1]]$choice, "0")
  expect_identical(vapply(out, function(ch) ch$alloc$last$choice, ""), rep("0", 3))
  expect_true(out[[1]]$alloc$last$choice_fixed)
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
  within_g <- mean(c(stats::var(emc[[1]]$samples$theta_mu[burn + seq_len(m)]),
                     stats::var(emc[[2]]$samples$theta_mu[burn + seq_len(m)])))
  expect_equal(prod$g_var[1], within_g)
  expect_identical(prod$pairs[["group"]], w_n + 50L)
  expect_equal(prod$pool$G, diag(c(50, 0)))
  expect_equal(prod$g_rho[1], min(1, EMC2:::.emc_alloc_rho(50 / (w_n + 50), within_g) /
                                       (1 - prod$g_white[1])))
  expect_equal(prod$g_b_share[1], max(0, stats::var(x) - within_g) / stats::var(x))
  # The first production rebuild measures the window alone and starts the pool.
  first <- emc
  first[[1]]$alloc$stage <- "adapt"
  m1 <- EMC2:::.emc_alloc_measure(first, 1, NULL, diag(1), "sample")
  expect_identical(m1$pairs[["group"]], w_n)
  expect_null(m1$pool)
  # Outside production: the window's draws, the window's jumps.
  b <- EMC2:::.emc_alloc_measure(emc, 1, NULL, diag(1), "burn")
  k <- EMC2:::.emc_alloc_draws(emc[[1]], w_n)
  within_b <- mean(c(stats::var(emc[[1]]$samples$theta_mu[k]),
                     stats::var(emc[[2]]$samples$theta_mu[k])))
  expect_equal(b$g_var[1], within_b)
  expect_identical(b$pairs[["group"]], w_n)
})

test_that("a group coordinate is priced on its slow part, the group's draw given the subjects being white", {
  share <- EMC2:::.emc_alloc_white_share
  cap <- EMC2:::.EMC_ALLOC_WHITE_CAP
  S <- matrix(c(2, 0.5, 0.5, 1), 2)
  Dinv <- matrix(c(1, 1, 0, 2), 2)
  var <- c(1, 4, 0.05, 0.5)
  # Mean coordinates: d' S d / N of their variance; log sds: 1 / (2 N).
  w <- share(S, Dinv, 10, var, rho = rep(0, 4))
  expect_equal(w[1:2], rowSums((Dinv %*% S) * Dinv) / 10 / var[1:2])
  expect_equal(w[4], 1 / 20 / 0.5)
  # No more than the cap's share of what decorrelated at lag one; nothing
  # where nothing did, and nothing without a group covariance.
  expect_equal(w[3], cap)
  expect_equal(share(S, Dinv, 10, var, rho = c(0.9, 0.9, 1, 0.99))[3:4], c(0, cap * 0.01))
  expect_equal(share(NULL, Dinv, 10, var, rho = rep(0, 4))[1:2], c(0, 0))

  # An iteration of white noise of share s on an AR(1) part of coefficient r
  # has lag-one autocorrelation s r; the unit's rho is r again, and a move's
  # jump counts against the slow part's variance.
  set.seed(9)
  p <- 1; N <- 25; k <- 4000; s_true <- 0.2; r <- 0.8
  slow <- as.numeric(stats::filter(rnorm(k), r, "recursive")) * sqrt(1 - r^2) * sqrt(1 - s_true)
  mu <- slow + rnorm(k, sd = sqrt(s_true))
  Sig <- s_true * N
  ch <- list(n_pars = p, n_subjects = N, type = "diagonal-gamma", par_names = "x",
             lik_laplace = rep(list(NULL), N),
             samples = list(idx = k, theta_mu = matrix(mu, 1), theta_var = array(Sig, c(1, 1, k)),
                            alpha = array(0, c(1, N, k))),
             shear = list(plan = list(dirs = diag(1),
                                      last = list(E = c(0.2, 0), n = c(10L, 0L), moved = integer(2)))),
             alloc = list(window = modifyList(EMC2:::.emc_alloc_empty_window(N, p),
                                              list(n = k - 1L, nG = k - 1L,
                                                   G = diag(c(sum(diff(mu)^2), 0))))))
  local_mocked_bindings(.emc_shear_active = function(sampler) TRUE,
                        .emc_alloc_draws = function(ch, m) seq_len(ch$samples$idx),
                        .package = "EMC2")
  ms <- EMC2:::.emc_alloc_measure(list(ch), 1, NULL, diag(p))
  v <- ms$g_var[1]
  lag1 <- EMC2:::.emc_alloc_rho(sum(diff(mu)^2) / (k - 1), v)
  expect_equal(ms$g_white[1], min(s_true * N / N / v, cap * (1 - lag1)))
  expect_equal(ms$g_rho[1], lag1 / (1 - ms$g_white[1]))
  expect_equal(ms$g_rho[1], r, tolerance = 0.03)
  g_obs <- 1 - (0.2 / 10) / (2 * v * (1 - ms$g_white[1]))
  expect_equal(ms$g_fac[1], g_obs)
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
  local_mocked_bindings(.emc_shear_active = function(sampler) TRUE, .package = "EMC2")
  ms <- EMC2:::.emc_alloc_measure(list(ch), 1, NULL, diag(p))
  g1 <- max(0, 1 - 1 / (2 * ms$g_var[1] * (1 - ms$g_white[1])))
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

test_that("measured zero movement is a factor of exactly one, and buys nothing", {
  p <- 2; n <- 1; k <- 50
  a <- matrix(stats::rnorm(p * k), p)
  ch <- list(n_pars = p, n_subjects = n, type = "diagonal-gamma", par_names = c("x", "y"),
             lik_laplace = list(NULL),
             samples = list(idx = k, theta_mu = a, theta_var = array(diag(2), c(p, p, k)),
                            alpha = array(a, c(p, 1, k))),
             shear = list(plan = list(dirs = diag(p),
                                      last = list(E = c(0, 0, 0, 0), n = c(10L, 0L, 0L, 0L),
                                                  moved = integer(4)))),
             alloc = list(window = NULL))
  local_mocked_bindings(.emc_shear_active = function(sampler) TRUE, .package = "EMC2")
  ms <- EMC2:::.emc_alloc_measure(list(ch), 1, NULL, diag(p))
  expect_identical(ms$g_fac[1], 1)
  # A coordinate at rho 1 whose rounds do not move it is not a kernel the
  # solver can use: no rounds, no fabricated gain.
  sol <- EMC2:::.emc_alloc_solve(list(list(kind = "shear", rho = 1, g = 1, cost = 0.05)), 1)
  expect_identical(sol$x, 0)
  expect_equal(sol$objective, sol$objective0)
})

test_that("a measured collapsed move that never moves earns no credit", {
  p <- 2; n <- 1; k <- 50
  a <- matrix(stats::rnorm(p * k), p)
  v <- array(0, c(p, p, k))
  for (i in seq_len(k)) v[, , i] <- diag(exp(stats::rnorm(p)))
  ch <- list(n_pars = p, n_subjects = n, type = "diagonal-gamma", par_names = c("x", "y"),
             lik_laplace = list(NULL),
             samples = list(idx = k, theta_mu = a, theta_var = v,
                            alpha = array(a, c(p, 1, k))),
             collapsed = list(plan = list(C = diag(2 * p), rounds = 1,
                                          last = list(E = matrix(0, 2 * p, 2 * p), n = 600,
                                                      acc = 0, moved = 0L))),
             alloc = list(window = NULL))
  local_mocked_bindings(.emc_collapsed_active = function(sampler) TRUE, .package = "EMC2")
  ms <- EMC2:::.emc_alloc_measure(list(ch), 1, NULL, diag(p))
  expect_true(ms$c_measured)
  expect_identical(ms$c_fac, rep(1, 2 * p))
  # Movement on one coordinate is shaded towards the move's own mean factor.
  ch$collapsed$plan$last$E[1, 1] <- 60
  ms <- EMC2:::.emc_alloc_measure(list(ch), 1, NULL, diag(p))
  obs <- c(max(0, 1 - 0.1 / (2 * ms$g_var[1] * (1 - ms$g_white[1]))), rep(1, 3))
  expect_equal(ms$c_fac, (600 * obs + 20 * mean(obs)) / 620)
})

test_that("subject HMC factors weigh their transitions against the subjects' pooled factor", {
  n0 <- EMC2:::.EMC_ALLOC_ROUND_EVIDENCE
  p <- 1
  ev <- function(n, M) list(n = n, M = matrix(M), cpu = n * 0.01)
  ms <- list(s_rho = matrix(0.9, 2, 1), g_rho = 0.5, g_fac = 0.5, g_measured = FALSE,
             C_rest = 1, c_round = 0.1, coll = FALSE, W = list(diag(1), diag(1)),
             s_var = matrix(1, 2, 1), c_hmc = c(0.1, 0.1), c_part = c(NA, NA), hmc_fixed = 0)
  st <- list(supported = TRUE,
             subjects = list(list(eligible = TRUE), list(eligible = TRUE)),
             evidence = list(ev(10, 16), ev(100, 0)))  # subject 1 moves (h 0.2), subject 2 does not (h 1)
  u <- EMC2:::.emc_alloc_units(ms, st, TRUE, FALSE)
  h1 <- 1 - 16 / 10 / 2; h2 <- 1
  ref <- (10 * h1 + 100 * h2) / 110
  expect_equal(u$units[[2]]$h, (10 * h1 + n0 * ref) / (10 + n0))
  expect_equal(u$units[[3]]$h, (100 * h2 + n0 * ref) / (100 + n0))
  # Decay halves a subject's evidence at a production rebuild.
  d <- EMC2:::.emc_hmc_decay_evidence(ev(10, 16))
  expect_equal(c(d$n, d$M[1], d$cpu), EMC2:::.EMC_ALLOC_EVIDENCE_DECAY * c(10, 16, 0.1))
})

test_that("a subject without HMC in the block is priced by the rows its transitions evaluate", {
  ms <- list(s_rho = matrix(0.9, 2, 1), g_rho = 0.5, g_fac = 0.5, g_measured = FALSE,
             C_rest = 1, c_round = 0.1, coll = FALSE, W = list(diag(1), diag(1)),
             s_var = matrix(1, 2, 1), c_hmc = c(0.1, NA), c_part = c(NA, NA), hmc_fixed = 0,
             data_rows = c(100, 300), hmc_row_cpu = 2e-6)
  ev <- function(cpu) list(n = 10, M = matrix(8), cpu = cpu, arows = c(20, 40))
  st <- list(supported = TRUE, subjects = list(list(eligible = TRUE), list(eligible = TRUE)),
             evidence = list(ev(1), ev(7)))
  u <- EMC2:::.emc_alloc_units(ms, st, TRUE, FALSE)
  # Measured in the block: its own cost. Otherwise mean rows per transition
  # x data rows x the chain's time per row, whatever its pilots' CPU was.
  expect_equal(u$units[[2]]$c_h, 0.1)
  expect_equal(u$units[[3]]$c_h, 30 * 300 * 2e-6)
  st$evidence[[2]]$cpu <- 99
  expect_equal(EMC2:::.emc_alloc_units(ms, st, TRUE, FALSE)$units[[3]]$c_h, 30 * 300 * 2e-6)
  # Without a time per row, the evidence's own CPU.
  ms$hmc_row_cpu <- NA_real_
  expect_equal(EMC2:::.emc_alloc_units(ms, st, TRUE, FALSE)$units[[3]]$c_h, 99 / 10)
  expect_equal(EMC2:::.emc_alloc_data_rows(list(data.frame(x = 1:3), data.frame(x = 1:2))), 5)
})

test_that("a move earns its evidence at an exploration rate and runs no faster than its evidence allows", {
  n0 <- EMC2:::.EMC_ALLOC_ROUND_EVIDENCE
  B <- 100
  ms <- list(s_rho = matrix(0.5, 1, 1), g_rho = c(0.99, 0.6), g_fac = c(0.5, 0.9),
             g_measured = c(FALSE, TRUE), g_n = c(0, 400), C_rest = 1, c_round = 0.05,
             coll = TRUE, c_fac = rep(0.5, 2), c_measured = FALSE, c_n = 0, c_coll = 0.1,
             c_rate = 0, W = list(diag(1)), s_var = matrix(1, 1, 1), c_hmc = NA_real_)
  u <- EMC2:::.emc_alloc_units(ms, NULL, TRUE, FALSE, B = B)
  # No rounds yet: the floor collects n0 rounds by the next rebuild; a shear
  # coordinate may still run once per iteration at its pooled factor, while
  # the collapsed move is capped at its floor. 400 rounds: no floor, the
  # usual cap.
  expect_equal(c(u$units[[1]]$floor, u$units[[1]]$cap), c(n0 / B, 1))
  expect_equal(c(u$units[[2]]$floor, u$units[[2]]$cap), c(0, EMC2:::.EMC_ALLOC_MAX_ROUNDS))
  expect_equal(c(u$joint$floor, u$joint$cap), c(n0 / B, n0 / B))
  sol <- EMC2:::.emc_alloc_solve_units(u)
  # The unmeasured coordinate is the slowest, so it is bought past its floor.
  expect_gt(sol$x[1], n0 / B)
  expect_lte(sol$x[1], 1 + 1e-9)
  expect_equal(sol$joint, n0 / B)
  # Half the evidence in: the floor collects the rest, and a measured move
  # may run at one move per iteration.
  half <- ms; half$g_n <- c(n0 / 2, 400); half$c_n <- n0 / 2
  uh <- EMC2:::.emc_alloc_units(half, NULL, TRUE, FALSE, B = B)
  expect_equal(c(uh$units[[1]]$floor, uh$units[[1]]$cap), c(n0 / 2 / B, 1))
  expect_equal(uh$joint$floor, n0 / 2 / B)
  # Without B (no rebuild horizon) nothing is floored or capped by evidence.
  u0 <- EMC2:::.emc_alloc_units(ms, NULL, TRUE, FALSE)
  expect_equal(u0$units[[1]]$floor, 0)
  expect_equal(u0$joint$floor, 0)
  # The plan pays for the exploration it runs.
  expect_gte(sol$cost, 1 + n0 / B * (0.05 + 0.1) - 1e-9)
})

test_that("the evidence cap is twice the rate over the iterations the evidence spans", {
  # A move run at a steady rate over production blocks that double, its
  # evidence and the iterations behind it discounted alike: the cap is twice
  # that rate whatever the next block's length.
  d <- EMC2:::.EMC_ALLOC_EVIDENCE_DECAY
  x <- 1.5
  iters <- 200 + d * 100
  ms <- list(s_rho = matrix(0.5, 1, 1), g_rho = c(0.99, 0.6), g_fac = c(0.5, 0.9),
             g_measured = c(TRUE, TRUE), g_n = c(x * iters, 400), g_iters = c(iters, 0),
             C_rest = 1, c_round = 0.05, coll = TRUE, c_fac = rep(0.5, 2), c_measured = TRUE,
             c_n = x * iters, c_iters = iters, c_coll = 0.1, c_rate = 0, W = list(diag(1)),
             s_var = matrix(1, 1, 1), c_hmc = NA_real_)
  u <- EMC2:::.emc_alloc_units(ms, NULL, TRUE, FALSE, B = 400)
  expect_equal(u$units[[1]]$cap, 2 * x)
  expect_equal(u$joint$cap, 2 * x)
  # Without the span the evidence is read against the next block.
  expect_equal(u$units[[2]]$cap, 2 * 400 / 400)
  ms$g_iters <- NULL
  expect_equal(EMC2:::.emc_alloc_units(ms, NULL, TRUE, FALSE, B = 400)$units[[1]]$cap,
               2 * x * iters / 400)
  # A measured move that ran little may still run at one move per iteration:
  # a plan that ran none of it does not cap it down for good.
  low <- ms; low$g_iters <- c(iters, 0); low$g_n <- c(0.01 * iters, 400)
  low$c_n <- 0.01 * iters
  ul <- EMC2:::.emc_alloc_units(low, NULL, TRUE, FALSE, B = 400)
  expect_equal(c(ul$units[[1]]$cap, ul$joint$cap), c(1, 1))
  # The span of a closing block is the allocation window's iteration count.
  expect_equal(EMC2:::.emc_alloc_block_iters(list(alloc = list(window = list(n = 37L)))), 37)
  expect_equal(EMC2:::.emc_alloc_block_iters(list()), 0)
})

test_that("a chain held to shared particle rates is solved at those rates", {
  fixed <- list(kind = "fixed", rho = 0.2, cost = 0)
  subj <- list(kind = "subject", rho = c(0.95, 0.9), p = c(0.6, 0.5), c_p = 1,
               h = c(0.3, 0.3), c_h = 0.2)
  own <- EMC2:::.emc_alloc_solve(list(fixed, subj), 1)
  expect_lt(own$q[2], 1)
  held <- EMC2:::.emc_alloc_solve(list(fixed, replace(subj, "q_fixed", list(1))), 1)
  expect_identical(held$q[2], 1)
  expect_equal(held$cost, 1 + 1 * subj$c_p + held$x[2] * subj$c_h)
  expect_gte(held$cost, own$cost - 1e-9)
})

test_that("a few rounds shade a coordinate's factor towards its kind's, and log-sd rounds pool across chains", {
  p <- 2; n <- 1; k <- 50
  a <- matrix(stats::rnorm(p * k), p)
  tv <- array(0, c(p, p, k)); for (t in seq_len(k)) tv[, , t] <- diag(exp(stats::rnorm(p)))
  chain <- function(last, stage = "sample") list(
    n_pars = p, n_subjects = n, type = "diagonal-gamma", par_names = c("x", "y"),
    lik_laplace = list(NULL),
    samples = list(idx = k, stage = rep("sample", k), theta_mu = a,
                   theta_var = tv, alpha = array(a, c(p, 1, k))),
    shear = list(plan = list(dirs = diag(p), stage = stage, last = last)),
    alloc = list(window = NULL))
  local_mocked_bindings(.emc_shear_active = function(sampler) TRUE, .package = "EMC2")
  n0 <- EMC2:::.EMC_ALLOC_ROUND_EVIDENCE
  # Mean direction 1: 50 rounds moving well; direction 2: 5 rounds that did
  # not move at all. The five shade direction 2 towards direction 1's factor
  # instead of declaring it immovable.
  ch <- chain(list(E = c(50, 0, 0, 0), n = c(50, 5, 0, 0), moved = integer(4)))
  ms <- EMC2:::.emc_alloc_measure(list(ch), 1, NULL, diag(p))
  g1 <- max(0, 1 - 1 / (2 * ms$g_var[1] * (1 - ms$g_white[1])))
  # Two measured coordinates leave the spread at its floor 2 / n0: log
  # movements shrink towards their precision-weighted mean, n0 rounds' worth.
  t2 <- 2 / n0
  l <- log(c(1 - g1, EMC2:::.EMC_ALLOC_MIN_MOVE)); v <- 2 / c(50, 5)
  ref <- sum(l / (t2 + v)) / sum(1 / (t2 + v))
  lk <- (l / v + ref / t2) / (1 / v + 1 / t2)
  expect_equal(ms$g_fac[1:2], 1 - exp(lk))
  expect_lt(ms$g_fac[2], 1)
  # With evidence the measured zero movement wins: 500 rounds of nothing.
  ch2 <- chain(list(E = c(50, 0, 0, 0), n = c(50, 500, 0, 0), moved = integer(4)))
  ms2 <- EMC2:::.emc_alloc_measure(list(ch2), 1, NULL, diag(p))
  expect_gt(ms2$g_fac[2], 0.98)
  # A second chain's log-sd rounds count for this chain's log sds, not for
  # its mean directions.
  other <- chain(list(E = c(0, 0, 1, 0), n = c(0, 0, 10, 0), moved = integer(4)))
  mine <- chain(list(E = c(0, 0, 0, 0), n = c(0, 0, 0, 0), moved = integer(4)))
  msp <- EMC2:::.emc_alloc_measure(list(mine, other), 1, NULL, diag(p), stage = "sample")
  g3 <- max(0, 1 - 0.1 / (2 * msp$g_var[3] * (1 - msp$g_white[3])))
  expect_true(msp$g_measured[3]); expect_false(msp$g_measured[1])
  expect_equal(msp$g_fac[3], (10 * g3 + n0 * g3) / (10 + n0))
  # Outside production nothing pools.
  msb <- EMC2:::.emc_alloc_measure(list(mine, other), 1, NULL, diag(p), stage = "burn")
  expect_false(msb$g_measured[3])
})

test_that("a well-measured slow coordinate keeps its own factor among quick ones", {
  n0 <- EMC2:::.EMC_ALLOC_ROUND_EVIDENCE
  # Four quick mean directions and one slow ridge, all well measured.
  g_obs <- c(0.6, 0.55, 0.65, 0.6, 0.94)
  n <- c(50, 50, 50, 50, 60)
  g <- EMC2:::.emc_alloc_shrink_g(g_obs, n, rep(1L, 5))
  # The spread between them far exceeds their noise, so the ridge keeps its
  # factor within its noise, where a fixed n0-round weight credits it with
  # the quick directions' movement.
  fixed <- (60 * 0.94 + n0 * sum(n[1:4] * g_obs[1:4]) / 200) / (60 + n0)
  expect_lt(abs(log(1 - g[5]) - log(0.06)), 2 * sqrt(2 / 60))
  expect_gt(g[5], fixed + 0.02)
  # Equal coordinates still pool: noise alone is shrunk away.
  same <- EMC2:::.emc_alloc_shrink_g(c(0.7, 0.75, 0.72), c(5, 5, 5), rep(1L, 3))
  expect_lt(diff(range(same)), diff(range(c(0.7, 0.75, 0.72))))
  # An unmeasured coordinate takes its kind's centre; with nothing measured,
  # the prior.
  un <- EMC2:::.emc_alloc_shrink_g(c(0.6, 0.6, 0.6, NA), c(30, 30, 30, 0), rep(1L, 4))
  expect_equal(un[4], 0.6)
  expect_equal(EMC2:::.emc_alloc_shrink_g(rep(NA, 2), c(0, 0), 1:2),
               rep(EMC2:::.EMC_ALLOC_SHEAR_PRIOR_G, 2))
})

test_that("rounds are priced at their own time, the frames as a fixed cost paid once", {
  s <- list(n_pars = 2, n_subjects = 3)
  for (i in 1:12) s <- EMC2:::.emc_alloc_record_cost(s, 1, list(NULL, NULL, NULL),
    moves = list(shear = c(cpu = 2.5, own = 0.5, n = 5), collapsed = c(cpu = 1, own = 0.8, n = 1)))
  s <- EMC2:::.emc_alloc_record_cost(s, 1, list(NULL, NULL, NULL),
    moves = list(shear = c(cpu = 0, own = 0, n = 0)))
  w <- s$alloc$window
  expect_equal(w$fam$shear$own, 6)
  expect_identical(w$fam$shear$iters, 12L)
  expect_identical(w$fam$shear$moves, 60L)
  ms <- EMC2:::.emc_alloc_marginal(w, "shear"); expect_equal(c(ms$c, ms$fixed), c(0.1, 2))
  mc <- EMC2:::.emc_alloc_marginal(w, "collapsed"); expect_equal(c(mc$c, mc$fixed), c(0.8, 0.2))
  # The collapsed move pays its fixed cost once when any move runs.
  fixed <- list(kind = "fixed", rho = 0.95, cost = 0)
  joint <- list(idx = 1L, g = 0.3, cost = 0.05, cap = 4, fixed = 0.4)
  sol <- EMC2:::.emc_alloc_solve(list(fixed), 1, joint = joint)
  expect_gt(sol$joint, 0)
  expect_equal(sol$cost, 1 + 0.4 + sol$joint * 0.05)
  expect_identical(EMC2:::.emc_alloc_solve(list(fixed), 1, joint = replace(joint, "fixed", list(80)))$joint, 0)
  # A plan with rounds pays the fixed cost; one whose rounds cannot repay it
  # runs none.
  fixed <- list(kind = "fixed", rho = 0.2, cost = 0)
  sh <- list(kind = "shear", rho = 0.95, g = 0.4, cost = 0.02)
  u <- list(units = list(sh, fixed), C_rest = 1, joint = NULL, hmc_fixed = 0, shear_fixed = 0.3)
  sol <- EMC2:::.emc_alloc_solve_units(u)
  expect_gt(sol$x[1], 0)
  expect_equal(sol$cost, 1 + 0.3 + sol$x[1] * 0.02)
  expect_equal(sol$tau0, EMC2:::.emc_alloc_tau(0.95))
  dear <- EMC2:::.emc_alloc_solve_units(replace(u, "shear_fixed", list(50)))
  expect_identical(dear$x[1], 0)
  expect_equal(dear$cost, 1)
})

test_that("an unmeasured shear coordinate counts at the prior factor, uncapped", {
  ms <- list(s_rho = matrix(0.5, 1, 1), g_rho = c(0.976, 0.6), g_fac = c(0.5, 0.9),
             g_measured = c(FALSE, TRUE), C_rest = 1, c_round = 0.05, coll = TRUE,
             c_fac = c(0.965, 0.5), c_measured = TRUE, c_coll = 0.1, c_rate = 1,
             W = list(diag(1)), s_var = matrix(1, 1, 1), c_hmc = NA_real_)
  u <- EMC2:::.emc_alloc_units(ms, NULL, TRUE, FALSE)
  expect_equal(u$units[[1]]$cap, EMC2:::.EMC_ALLOC_MAX_ROUNDS)
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
    .emc_hmc_support = function(sampler) list(ok = FALSE),
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

test_that("a family's fixed cost is charged at the plan's activation probability", {
  sh <- list(kind = "shear", rho = 0.9, g = 0.5, cost = 0.05)
  subj <- list(kind = "subject", rho = 0.9, p = NULL, c_p = 0, h = 0.5, c_h = 0.1)
  act <- EMC2:::.emc_alloc_activation
  # independent fractional rounds; a whole round always activates
  expect_equal(act(list(sh, sh), cbind(c(0.1, 0.2), c(1.5, 0)), c(0, 0), 10, 0),
               10 * c(1 - 0.9 * 0.8, 1))
  # HMC subjects at their total rates, including a held floor
  expect_equal(act(list(subj, subj, sh), cbind(c(0.5, 0.5, 3)), c(0, 0.5, 0), 0, 4),
               4 * (1 - 0.5 * 0))
  expect_equal(act(list(subj, sh), cbind(c(0, 1)), c(0, 0), 0, 4), 0)
  # The compiled product is R's long-double prod() exactly, non-finite rates
  # counting as one.
  set.seed(3)
  X <- matrix(runif(400, -0.1, 1.2), 80)
  X[c(3, 90, 200)] <- c(NA, Inf, -Inf)
  add <- runif(80) * (runif(80) < 0.5)
  R <- pmin(X + add, 1)
  R[!is.finite(R)] <- 1
  expect_identical(EMC2:::emc_alloc_act_term_cpp(X, add), 1 - apply(1 - R, 2, prod))
  # A solved plan pays the fixed cost only on the iterations that run rounds.
  sol <- EMC2:::.emc_alloc_solve(list(sh), 1, shear_fixed = 0.5)
  expect_gt(sol$x, 0)
  expect_equal(sol$cost, 1 + sol$x * 0.05 + 0.5 * min(1, sol$x))
})

test_that("a subject's rates are priced as the kernels are selected, jointly", {
  need <- EMC2:::.emc_alloc_unit_need
  u <- list(kind = "subject", rho = 0.9, p = 0.6, c_p = 1, h = 0.2, c_h = 1, q_fixed = 0.6)
  rs <- c(0.5, 0.3, 0.1)
  got <- need(u, rs, 10)
  # One target met while the kernels exclude each other (r <= 1 - q), one
  # only once they overlap, one out of reach.
  expect_lt(got$x[1], 1 - 0.6)
  expect_gt(got$x[2], 1 - 0.6)
  for (k in 1:2)
    expect_equal(alloc_subject_rho(u, 0.6, got$x[k]), rs[k], tolerance = 1e-10)
  expect_identical(got$cost[3], Inf)
  # A held HMC rate stands where it meets the target and nowhere else.
  held <- need(replace(u, "r_fixed", list(0.7)), rs, 10)
  expect_identical(held$x[1:2], c(0.7, 0.7))
  expect_identical(held$cost[3], Inf)
  held <- need(replace(u, "r_fixed", list(0.5)), rs, 10)
  expect_identical(held$x[1], 0.5)
  expect_identical(held$cost[2:3], c(Inf, Inf))
})

# A subject unit with every optional field drawn at random.
alloc_rand_subject <- function(k) {
  u <- list(kind = "subject", rho = stats::runif(k, 0.3, 0.999),
            p = stats::runif(k, 0.05, 0.95), c_p = stats::runif(1, 0.05, 0.5),
            h = stats::runif(k, 0, 1.1), c_h = stats::runif(1, 0.02, 0.5))
  if (stats::runif(1) < 0.15) u$p <- NULL
  if (stats::runif(1) < 0.15) u$h <- NULL
  if (stats::runif(1) < 0.2) u$q_fixed <- sample(EMC2:::.EMC_ALLOC_PARTICLE_RATES, 1)
  if (stats::runif(1) < 0.2) u$r_fixed <- stats::runif(1)
  if (stats::runif(1) < 0.3) u$r_min <- sample(c(0.5, 1), 1)
  u
}

test_that("subjects priced together match each subject priced alone", {
  local_test_rng(21)
  for (trial in 1:20) {
    us <- c(list(list(kind = "shear", rho = 0.9, g = 0.6, cost = 0.05),
                 list(kind = "fixed", rho = 0.5, cost = 0)),
            lapply(stats::rpois(12, 2) + 1L, alloc_rand_subject))
    movable <- c(TRUE, FALSE, stats::runif(12) < 0.85)
    rs <- sort(c(stats::runif(25, 0.05, 1), 1, 0.999999))
    got <- EMC2:::.emc_alloc_units_need(us, movable, rs, 4)
    for (i in seq_along(us)) {
      ref <- if (movable[i]) EMC2:::.emc_alloc_unit_need(us[[i]], rs, 4)
        else EMC2:::.emc_alloc_unit_still(us[[i]], length(rs))
      expect_equal(got[[i]]$cost, ref$cost, tolerance = 1e-12)
      ok <- is.finite(ref$cost)
      expect_equal(got[[i]]$x[ok], ref$x[ok], tolerance = 1e-12)
      if (!is.null(ref$q)) expect_identical(got[[i]]$q[ok], ref$q[ok])
    }
  }
})

test_that("shear coordinates priced together match each priced alone, bit for bit", {
  local_test_rng(23)
  rand_shear <- function() {
    g <- sample(c(stats::runif(1, 0.05, 0.98), 0, 1, 1 - 1e-10), 1, prob = c(0.85, 0.05, 0.05, 0.05))
    u <- list(kind = "shear", rho = stats::runif(1, 0.05, 0.999), g = g,
              cost = stats::runif(1, 0.01, 0.2))
    r <- stats::runif(1)
    if (r < 0.4) u$a <- g * stats::runif(1)
    else if (r < 0.5) u$a <- -0.1
    else if (r < 0.7) u$a <- g
    if (stats::runif(1) < 0.5) u$cap <- stats::runif(1, 0, 4)
    if (stats::runif(1) < 0.4) u$floor <- stats::runif(1, 0, 0.5)
    u
  }
  for (trial in 1:20) {
    us <- replicate(15, rand_shear(), simplify = FALSE)
    movable <- stats::runif(15) < 0.9
    rs <- sort(c(stats::runif(40, 0, 1), 0, 1))
    got <- EMC2:::.emc_alloc_units_need(us, movable, rs, 4)
    for (i in seq_along(us)) {
      ref <- if (movable[i]) EMC2:::.emc_alloc_unit_need(us[[i]], rs, 4)
        else EMC2:::.emc_alloc_unit_still(us[[i]], length(rs))
      expect_identical(got[[i]], ref)
    }
  }
})

test_that("a plan's predicted autocorrelations are those it was solved for", {
  local_test_rng(22)
  for (trial in 1:15) {
    units <- c(list(list(kind = "shear", rho = stats::runif(1, 0.5, 0.99),
                         g = stats::runif(1, 0.2, 0.95), cost = stats::runif(1, 0.01, 0.2)),
                    list(kind = "fixed", rho = stats::runif(1, 0, 0.5), cost = 0)),
               lapply(1:3, function(i)
                 list(kind = "subject", rho = stats::runif(2, 0.3, 0.999),
                      p = stats::runif(2, 0.05, 0.95), c_p = stats::runif(1, 0.05, 0.5),
                      h = stats::runif(2, 0, 0.95), c_h = stats::runif(1, 0.02, 0.5))))
    sol <- EMC2:::.emc_alloc_solve(units, 1)
    pred <- sol$pred
    expect_length(pred, length(units))
    expect_equal(pred[[1]], units[[1]]$rho *
                   EMC2:::.emc_alloc_factor(units[[1]]$g, sol$x[1], units[[1]]$g))
    expect_equal(pred[[2]], units[[2]]$rho)
    for (i in 3:5) expect_equal(pred[[i]], alloc_subject_rho(units[[i]], sol$q[i], sol$x[i]))
    # The solver's target bounds every coordinate it ranges over.
    expect_lte(max(unlist(pred[-2])), sol$rho + 1e-9)
  }
})

test_that("the solver's objective is the model's objective at the plan it returns", {
  local_test_rng(23)
  subj <- function() {
    list(kind = "subject", rho = stats::runif(1, 0.3, 0.98), p = stats::runif(1, 0.1, 0.9),
         c_p = stats::runif(1, 0.05, 0.4), h = stats::runif(1, 0, 0.9),
         c_h = stats::runif(1, 0.02, 0.4), d_rho = stats::runif(1, 0.1, 0.9),
         p_dev = stats::runif(1, 0.1, 0.9), h_dev = stats::runif(1, 0.2, 1),
         q_old = 1)
  }
  tau <- EMC2:::.emc_alloc_tau
  for (trial in 1:25) {
    rho_g <- stats::runif(1, 0.5, 0.99)
    u <- list(units = list(list(kind = "shear", rho = rho_g, g = stats::runif(1, 0.2, 0.95),
                                cost = stats::runif(1, 0.01, 0.2)), subj(), subj()),
              C_rest = 1, joint = NULL, hmc_fixed = stats::runif(1, 0, 0.2),
              shear_fixed = 0)
    # The subjects' deviations are coordinates of their units.
    u$units[2:3] <- lapply(u$units[2:3], EMC2:::.emc_alloc_dev_rows)
    sol <- EMC2:::.emc_alloc_solve_units(u)
    rates <- EMC2:::.emc_alloc_plan_rates(u$units, sol)
    q <- rates$q; r <- rates$r
    g <- u$units[[1]]
    rho <- g$rho * EMC2:::.emc_alloc_factor(g$g, sol$x[1], g$g)
    for (i in 2:3) {
      s <- u$units[[i]]
      rho_s <- alloc_subject_rho(s, q[i], r[i])
      # No subject coordinate, its deviations included, slower than under its
      # particle step alone.
      expect_true(all(rho_s <= s$rho * s$p + 1e-9))
      rho <- max(rho, rho_s)
    }
    cost <- 1 + sol$x[1] * u$units[[1]]$cost +
      sum(vapply(2:3, function(i) q[i] * u$units[[i]]$c_p + r[i] * u$units[[i]]$c_h, 0)) +
      u$hmc_fixed * (1 - prod(1 - r[2:3]))
    # The solver reports its target, which the rates meet: never below the
    # model's objective, and above it only by the target's resolution.
    expect_gte(sol$objective, tau(rho) * cost * (1 - 1e-9))
    expect_lte(sol$objective, tau(rho) * cost * (1 + 1e-3))
    expect_lte(sol$objective, sol$objective0 * (1 + 1e-12))
  }
})

test_that("measurement variances and log sds match their per-draw definitions", {
  local_test_rng(24)
  x <- matrix(stats::rnorm(5 * 40, sd = 3), 5)
  x[2, 7] <- NA
  expect_equal(EMC2:::.emc_alloc_row_var(x), apply(x, 1L, stats::var), tolerance = 1e-12)
  expect_identical(EMC2:::.emc_alloc_row_var(x[, 1, drop = FALSE]), rep(NA_real_, 5))
  p <- 3
  V <- array(0, c(p, p, 6), dimnames = list(c("a", "b", "c"), c("a", "b", "c"), NULL))
  for (t in 1:6) V[, , t] <- crossprod(matrix(stats::rnorm(p * p), p))
  V[2, 2, 4] <- 0
  k <- c(2L, 4L, 5L)
  ref <- vapply(k, function(t) EMC2:::.emc_alloc_lsd(V[, , t]), numeric(p))
  expect_identical(EMC2:::.emc_alloc_lsd_draws(V, k), ref)
})

test_that("chains measured from shared pooled covariances match direct variances", {
  local_test_rng(25)
  p <- 3; n <- 4; k <- 40
  chain <- function() {
    a <- array(stats::rnorm(p * n * k), c(p, n, k))
    mu <- matrix(stats::rnorm(p * k), p)
    tv <- array(0, c(p, p, k)); for (t in seq_len(k)) tv[, , t] <- diag(exp(stats::rnorm(p)))
    list(n_pars = p, n_subjects = n, type = "diagonal-gamma", par_names = c("x", "y", "z"),
         lik_laplace = vector("list", n),
         samples = list(idx = k, stage = rep("sample", k), theta_mu = mu, theta_var = tv, alpha = a),
         shear = list(plan = list(dirs = matrix(stats::rnorm(p * p), p) + diag(2, p),
                                  stage = "sample",
                                  last = list(E = numeric(2 * p), n = numeric(2 * p),
                                              moved = integer(2 * p)))),
         alloc = list(window = NULL))
  }
  emc <- list(chain(), chain())
  emc[[2]]$shear$plan$dirs <- emc[[1]]$shear$plan$dirs
  hst <- list(subjects = lapply(seq_len(n), function(s) list(W = matrix(stats::rnorm(p * p), p))))
  local_mocked_bindings(.emc_shear_active = function(sampler) TRUE, .package = "EMC2")
  cache <- new.env()
  for (j in 1:2) {
    shared <- EMC2:::.emc_alloc_measure(emc, j, hst, diag(p), "sample", cache = cache)
    alone <- EMC2:::.emc_alloc_measure(emc, j, hst, diag(p), "sample")
    expect_equal(shared, alone, tolerance = 1e-12)
  }
  # Against each chain's transformed draws: the solver uses within-chain
  # variance and carries the between-chain share as a separate rho floor.
  mu_by_chain <- lapply(emc, function(ch) ch$samples$theta_mu)
  lsd_by_chain <- lapply(emc, function(ch)
    vapply(seq_len(k), function(t) EMC2:::.emc_alloc_lsd(ch$samples$theta_var[, , t]), numeric(p)))
  Dinv <- solve(emc[[1]]$shear$plan$dirs)
  group_by_chain <- lapply(seq_along(emc), function(i)
    rbind(Dinv %*% mu_by_chain[[i]], lsd_by_chain[[i]]))
  group_within <- Reduce(`+`, lapply(group_by_chain, function(x)
    apply(x, 1L, stats::var))) / length(emc)
  expect_equal(unname(shared$g_var), group_within, tolerance = 1e-12)
  pooled_group <- do.call(cbind, group_by_chain)
  group_total <- apply(pooled_group, 1L, stats::var)
  expect_equal(shared$g_b_share, pmin(1, pmax(0, (group_total - group_within) / group_total)),
               tolerance = 1e-12)
  for (s in seq_len(n)) {
    xs <- lapply(emc, function(ch) hst$subjects[[s]]$W %*% ch$samples$alpha[, s, ])
    within <- Reduce(`+`, lapply(xs, function(x) apply(x, 1L, stats::var))) / length(xs)
    expect_equal(shared$s_var[s, ], within, tolerance = 1e-12)
  }
  # The deviations' covariance, computed once the window asks for it.
  st <- EMC2:::.emc_alloc_draw_stats(emc, lapply(emc, function(ch) seq_len(k)), TRUE)
  mu <- do.call(cbind, mu_by_chain)
  x <- do.call(cbind, lapply(emc, function(ch) ch$samples$alpha[, 2, ]))
  expect_equal(st$d_cov[[2]], stats::cov(t(x - mu)), tolerance = 1e-12)
})

test_that("a subject's deviations join its unit as coordinates the particle step must keep refreshed", {
  u <- list(kind = "subject", rho = c(0.5, 0.6), p = c(0.8, 0.9), c_p = 0.1, h = c(0.4, 0.5), c_h = 0.05,
            d_rho = c(0.7, 0.3), p_dev = c(0.5, 0.95), h_dev = c(0.9, 0.2), q_old = 0.5)
  d <- EMC2:::.emc_alloc_dev_rows(u)
  expect_length(d$rho, 4L)
  expect_equal(d$p, c(0.8, 0.9, 0.5, 0.95))
  expect_equal(d$h, c(0.4, 0.5, 0.9, 0.2))
  # The measured deviation rho has the particle step's factor at its rate divided out.
  expect_equal(d$rho[3:4], c(0.7, 0.3) / (1 - 0.5 * (1 - c(0.5, 0.95))))
  uf <- u
  uf$b_share <- c(0.2, 0.4)
  uf$d_b_share <- c(0.3, 0.6)
  df <- EMC2:::.emc_alloc_dev_rows(uf)
  den <- 1 - uf$q_old * (1 - uf$p_dev)
  # The measured rho is within-chain: only the particle factor comes out, and
  # the between-chain share rides along to be added back when priced.
  expect_equal(df$rho[3:4], pmin(1, uf$d_rho / den))
  expect_equal(df$b_share, c(uf$b_share, uf$d_b_share))
  # Without deviation factors for a kernel the unit has, nothing is added.
  expect_identical(EMC2:::.emc_alloc_dev_rows(u[setdiff(names(u), "h_dev")]), u[setdiff(names(u), "h_dev")])
  expect_identical(EMC2:::.emc_alloc_dev_rows(u[setdiff(names(u), "d_rho")]), u[setdiff(names(u), "d_rho")])
  # A particle step that refreshes a deviation HMC barely moves cannot be stood down.
  rs <- seq(0.3, 0.95, by = 0.05)
  nd <- EMC2:::.emc_alloc_unit_need(d, rs, 4)
  expect_true(all(nd$q[is.finite(nd$cost)] == 1))
  nd0 <- EMC2:::.emc_alloc_unit_need(u, rs, 4)
  expect_true(any(nd0$q[is.finite(nd0$cost)] < 1))
})

test_that("a cheaper HMC option within the minimum gain of the best is chosen", {
  o <- c(none = 10, "0" = 9.5, "1" = 9.2)
  expect_identical(EMC2:::.emc_alloc_choose(o, c(5, 4, 6)), "0")
  expect_identical(EMC2:::.emc_alloc_choose(o, c(3, 4, 6)), "none")
  expect_identical(EMC2:::.emc_alloc_choose(o, c(6, 6, 6)), "1")
  expect_identical(EMC2:::.emc_alloc_choose(c(none = 10, "0" = 8, "1" = 7.5), c(3, 4, 6)), "0")
  expect_identical(EMC2:::.emc_alloc_choose(o), "1")
  expect_identical(EMC2:::.emc_alloc_choose(c(none = Inf, "0" = 9, "1" = 9.5), c(Inf, 4, 6)), "0")
})

test_that("a curved shear coordinate runs a round in every iteration", {
  ms <- list(s_rho = matrix(0.5, 1, 1), g_rho = c(0.3, 0.98), g_fac = c(0.5, 0.95),
             g_n = c(200, 200), g_iters = c(100, 100), g_measured = c(TRUE, TRUE),
             C_rest = 1, c_round = 0.2, coll = FALSE, W = list(diag(1)),
             s_var = matrix(1, 1, 1), c_hmc = NA_real_)
  plain <- EMC2:::.emc_alloc_units(ms, NULL, TRUE, FALSE, B = 100)
  expect_identical(plain$units[[1]]$floor, 0)
  curved <- EMC2:::.emc_alloc_units(c(ms, list(g_curved = 1L)), NULL, TRUE, FALSE, B = 100)
  expect_identical(curved$units[[1]]$floor, 1)
  expect_gte(curved$units[[1]]$cap, 1)
  expect_identical(curved$units[[2]]$floor, 0)
  # The plan holds it there although the other coordinate sets the target.
  sol <- EMC2:::.emc_alloc_solve_units(curved)
  expect_gte(sol$x[1], 1)
  expect_lt(EMC2:::.emc_alloc_solve_units(plain)$x[1], 1)
})

# A reversible rest-of-iteration kernel on q coordinates with covariance cov:
# its jump covariance, slow directions decaying at rho.
joint_system <- function(q, rho, seed) {
  set.seed(seed)
  cov <- crossprod(matrix(rnorm(q * q), q)) / q + diag(q)
  e <- eigen(cov)
  h <- e$vectors %*% (sqrt(e$values) * t(e$vectors))
  o <- qr.Q(qr(matrix(rnorm(q * q), q)))
  list(cov = cov, G = 2 * h %*% (o %*% ((1 - rho) * t(o))) %*% h)
}

test_that("group coordinates priced jointly reduce to their own pricing when uncoupled", {
  v <- c(0.5, 1, 2, 0.7, 1.3, 0.9)
  rho <- c(0.6, 0.9, 0.97, 0.8, 0.95, 0.7)
  cp <- EMC2:::.emc_alloc_coupling(diag(v), v, diag(2 * v * (1 - rho)), matrix(0, 3, 3), diag(3),
                                   1e9, 1e6)
  expect_equal(diag(cp$Ar), rho, tolerance = 1e-6)
  g <- c(0.4, 0.8, 0.9, 0.5, 0.7, NA)
  a <- c(0.4, 0.5, 0.9, 0.2, 0.7, NA)
  x <- c(0.5, 2, 3.25, 0, 1, 0)
  fac <- c(1, 0.9, 0.8, 1, 0.95, 0.7)
  own <- diag(cp$Ar) * fac * c(mapply(EMC2:::.emc_alloc_factor, g[1:5], x[1:5], a[1:5]), 1)
  expect_equal(EMC2:::.emc_alloc_group_rho(cp, g, a, x, fac)[, 1], own, tolerance = 1e-10)
  # The cheapest rounds to a target are each coordinate's own, to the path's step.
  gf <- list(idx = 1:6, g = g, a = a, floor = rep(0, 6), cap = c(rep(4, 5), 0),
             cost = rep(0.05, 6), b = rep(0, 6))
  rs <- seq(0.72, 0.97, length.out = 12)
  nd <- EMC2:::.emc_alloc_group_need(cp, gf, rep(1, 6), c(rep(TRUE, 5), FALSE), rs)
  sep <- EMC2:::emc_alloc_shears_need_cpp(diag(cp$Ar)[1:5], rep(0, 5), g[1:5], a[1:5],
                                          gf$cost[1:5], gf$cap[1:5], gf$floor[1:5], rs)
  xj <- do.call(rbind, lapply(nd[1:5], `[[`, "x"))
  xs <- do.call(rbind, lapply(sep, `[[`, "x"))
  fin <- is.finite(colSums(xs))
  expect_identical(is.finite(Reduce(`+`, lapply(nd, `[[`, "cost"))), fin)
  expect_true(all(xj[, fin] >= xs[, fin] - 1e-9))
  expect_true(all(xj[, fin] <= xs[, fin] + 0.5))
})

test_that("a coupled group system prices each coordinate's rounds within its conditional", {
  q <- 8
  s <- joint_system(q, seq(0.5, 0.98, length.out = q), 11)
  cp <- EMC2:::.emc_alloc_coupling(s$cov, diag(s$cov), s$G, diag(4) * 0.01, diag(4), 50, 1e6)
  expect_identical(cp$shrink, EMC2:::.EMC_ALLOC_COUPLING_SHRINK)
  g <- c(0.5, 0.9, 0.7, 0.95, 0.6, 0.8, 0.85, NA)
  a <- g * c(1, 0.5, 0.8, 1, 0.6, 0.9, 0.7, NA)
  x <- c(1, 0.5, 2.5, 0, 3, 0.25, 1.75, 0)
  fac <- seq(1, 0.8, length.out = q)
  A <- diag(fac) %*% cp$Ar
  for (k in seq_len(q)) {
    th <- if (is.finite(g[k]) && x[k] > 0) 1 - EMC2:::.emc_alloc_factor(g[k], x[k], a[k]) else 0
    A <- (diag(q) - th * outer(replace(numeric(q), k, 1), cp$Q[k, ] / cp$Q[k, k])) %*% A
  }
  tau <- 2 * diag(solve(diag(q) - A) %*% cp$C) / diag(cp$C) - 1
  expect_equal(EMC2:::.emc_alloc_group_rho(cp, g, a, x, fac)[, 1], (tau - 1) / (tau + 1),
               tolerance = 1e-10)
  # A round that fully redraws one coordinate within its conditional leaves
  # the autocorrelation the other coordinates carry into it.
  none <- EMC2:::.emc_alloc_group_rho(cp, rep(NA_real_, q), rep(0, q), numeric(q))[, 1]
  k <- which.max(none)
  full <- EMC2:::.emc_alloc_group_rho(cp, replace(rep(NA_real_, q), k, 0), rep(0, q),
                                      replace(numeric(q), k, 1))[k, 1]
  expect_gt(full, 0)
  expect_lt(full, none[k])
})

test_that("the joint rounds for a target meet it within floors and caps", {
  q <- 8
  s <- joint_system(q, seq(0.6, 0.985, length.out = q), 5)
  cp <- EMC2:::.emc_alloc_coupling(s$cov, diag(s$cov), s$G, diag(4) * 0.01, diag(4), 50, 1e6)
  g <- c(0.5, 0.9, 0.7, 0.95, 0.6, 0.8, 0.85, NA)
  a <- g * c(1, 0.5, 0.8, 1, 0.6, 0.9, 0.7, NA)
  gf <- list(idx = seq_len(q), g = g, a = a, floor = c(0.25, rep(0, q - 1)),
             cap = c(rep(3, q - 1), 0), cost = seq(0.03, 0.06, length.out = q),
             b = rep(0.02, q))
  movable <- c(rep(TRUE, q - 1), FALSE)
  parts <- EMC2:::.emc_alloc_group_parts(cp, gf, rep(1, q))
  rs <- seq(max(parts$best[movable]) - 0.05, max(parts$levels[movable]) + 0.005, length.out = 30)
  nd <- EMC2:::.emc_alloc_group_need(cp, gf, rep(1, q), movable, rs)
  X <- do.call(rbind, lapply(nd, `[[`, "x"))
  cost <- Reduce(`+`, lapply(nd, `[[`, "cost"))
  expect_true(any(is.finite(cost)) && any(!is.finite(cost)))
  expect_equal(cost[length(rs)], sum(gf$floor * gf$cost))
  expect_true(all(diff(cost[is.finite(cost)]) <= 1e-12))
  for (t in which(is.finite(cost))) {
    r <- gf$b + (1 - gf$b) * EMC2:::.emc_alloc_group_rho(cp, g, a, X[, t])[, 1]
    expect_true(all(r[movable] <= rs[t] + 1e-9))
    expect_true(all(X[, t] <= gf$cap + 1e-9 & X[, t] >= gf$floor - 1e-9))
  }
  # Steps updated in place give the path a full rebuild at every step does.
  path <- function(rebuild) EMC2:::emc_alloc_group_need_cpp(
    cp$C, cp$Q, cp$Ar, rep(1, q), g, a, gf$floor, gf$cap, gf$cost, gf$b, numeric(q), movable, rs,
    rebuild)
  expect_equal(path(1000L)$x, path(1L)$x, tolerance = 1e-9)
  expect_identical(path(1000L)$reached, path(1L)$reached)
})

test_that("coupling is shrunk away where the window cannot support it", {
  q <- 6
  s <- joint_system(q, rep(0.9, q), 3)
  few <- EMC2:::.emc_alloc_coupling(s$cov, diag(s$cov), s$G, diag(3) * 0.01, diag(3), 50, q)
  expect_identical(few$shrink, 1)
  expect_equal(few$C[upper.tri(few$C)], numeric(15))
  expect_equal(few$Ar[upper.tri(few$Ar)], numeric(15))
  expect_null(EMC2:::.emc_alloc_coupling(s$cov, diag(s$cov), s$G * NA, diag(3), diag(3), 50, 1e3))
  expect_null(EMC2:::.emc_alloc_coupling(s$cov, diag(s$cov), s$G, diag(3), diag(3), 50, 1))
  expect_null(EMC2:::.emc_alloc_coupling(s$cov, diag(s$cov), -s$G, diag(3), diag(3), 50, 1e3))
})

test_that("the solver prices a coupled group block on the joint system", {
  q <- 4
  s <- joint_system(q, c(0.7, 0.9, 0.95, 0.98), 7)
  cp <- EMC2:::.emc_alloc_coupling(s$cov, diag(s$cov), s$G, diag(2) * 0.001, diag(2), 200, 1e6)
  g <- c(0.5, 0.6, 0.7, 0.6)
  units <- c(lapply(seq_len(q), function(k)
    list(kind = "shear", rho = cp$Ar[k, k], b_share = 0, g = g[k], a = g[k], cost = 0.05,
         floor = 0, cap = 4)),
    list(list(kind = "subject", rho = 0.5, b_share = 0, p = NULL, c_p = 0, h = NULL, c_h = 0)))
  sol <- EMC2:::.emc_alloc_solve(units, 1, coupling = cp[c("C", "Q", "Ar")])
  xg <- sol$x[seq_len(q)]
  expect_true(all(is.finite(xg)))
  expect_gt(sum(xg), 0)
  expect_equal(unlist(sol$pred[seq_len(q)]), EMC2:::.emc_alloc_group_rho(cp, g, g, xg)[, 1],
               tolerance = 1e-9)
  # Without the coupling the same units are priced one by one.
  alone <- EMC2:::.emc_alloc_solve(units, 1)
  expect_false(isTRUE(all.equal(alone$x[seq_len(q)], xg)))
})

test_that("round factors level off at what the rounds can refresh", {
  cap <- EMC2:::.emc_alloc_cap_rounds
  fac <- EMC2:::.emc_alloc_factor
  # A fifth of the variance is refreshable: no number of rounds does more.
  out <- cap(0.95, 0.95, 0.2)
  expect_identical(out$g, 0.95)
  expect_equal(fac(out$g, 1, out$a), 0.95)
  expect_equal(fac(out$g, 200, out$a), 0.8, tolerance = 1e-6)
  expect_gt(fac(out$g, 2, out$a), 0.95^2)
  # One round cannot refresh more than the share either.
  out <- cap(0.5, 0.5, 0.2)
  expect_equal(c(out$g, out$a), c(0.8, 0))
  expect_equal(fac(out$g, 3, out$a), 0.8)
  # A measured saturation below the cap stands; an unmeasured share changes nothing.
  expect_equal(cap(0.9, 0.1, 0.5)$a, 0.1)
  expect_identical(cap(c(0.9, NA, 0.7), c(0.6, NA, NA), c(NA, 0.3, NA)),
                   list(g = c(0.9, NA, 0.7), a = c(0.6, NA, NA)))
  expect_equal(cap(0.9, 0.9, 1)$a, 0.9)
})

test_that("a block's autocorrelation time allows for what its common mean removes", {
  series_tau <- EMC2:::.emc_alloc_series_tau
  centred <- EMC2:::.emc_alloc_centred_acov
  # The expected autocovariance of an estimate about the common mean: a
  # white series loses its mean's variance at lag zero, and an exponential
  # one more at every lag the slower it is.
  expect_equal(centred(NA_real_, 100L, 3L, 3L)[1], 1 - 1 / 300)
  expect_equal(centred(NA_real_, 100L, 3L, 3L)[3], -98 / (3 * 100^2))
  expect_equal(centred(0, 100L, 2L, 3L)[1], 1 - 1 / 300)
  slow <- centred(0.98, 200L, 50L, 3L)
  expect_lt(slow[1], 0.9)
  expect_true(all(slow < 0.98^(0:50)))
  set.seed(8)
  mc <- rowMeans(replicate(400, {
    x <- vapply(1:3, function(j) as.numeric(arima.sim(list(ar = 0.9), 60, n.start = 500)),
                numeric(60)) * sqrt(1 - 0.81)
    x <- x - mean(x)
    rowMeans(apply(x, 2, function(v) stats::acf(v, 5, type = "covariance", plot = FALSE,
                                                demean = FALSE)$acf))
  }))
  expect_equal(mc, centred(0.9, 60L, 5L, 3L), tolerance = 0.05)
  # A white part on a slow one: tau = w + (1 - w) tau_slow, which autocovariances
  # about each chain's own mean fall well short of.
  sim <- function(B, w, ts) {
    phi <- (ts - 1) / (ts + 1)
    lapply(1:3, function(j) sqrt(w) * rnorm(B) + sqrt((1 - w) * (1 - phi^2)) *
             as.numeric(arima.sim(list(ar = phi), B, n.start = 2000)))
  }
  est <- t(replicate(60, {
    xs <- sim(300L, 0.2, 60)
    c(series_tau(xs)$tau, EMC2:::.emc_shear_window_iact(xs))
  }))
  truth <- 0.2 + 0.8 * 60
  expect_gt(median(est[, 1]), 0.75 * truth)
  expect_lt(median(est[, 1]), 1.3 * truth)
  expect_lt(median(est[, 2]), 0.65 * truth)
  one <- series_tau(sim(300L, 0.2, 60))
  expect_true(one$white >= 0 && one$white < 1 && one$slow > one$tau)
  # A coordinate the block resolves takes the shape-free estimate.
  fast <- series_tau(sim(400L, 0, 3))
  expect_true(is.na(fast$slow))
  expect_gt(fast$tau, 2)
  expect_lt(fast$tau, 4.5)
  expect_true(is.na(series_tau(list(rnorm(10), rnorm(10)))$tau))
  expect_true(is.na(series_tau(list(c(rnorm(50), NA), rnorm(51)))$tau))
})

test_that("the group's excess raises the joint system's time", {
  q <- 8
  s <- joint_system(q, seq(0.6, 0.985, length.out = q), 5)
  cp <- EMC2:::.emc_alloc_coupling(s$cov, diag(s$cov), s$G, diag(4) * 0.01, diag(4), 50, 1e6)
  g <- c(0.5, 0.9, 0.7, 0.95, 0.6, 0.8, 0.85, NA)
  gf <- list(idx = seq_len(q), g = g, a = g, floor = numeric(q), cap = c(rep(3, q - 1), 0),
             cost = rep(0.05, q), b = rep(0.02, q))
  movable <- c(rep(TRUE, q - 1), FALSE)
  tau <- function(r) (1 + r) / (1 - r)
  plain <- EMC2:::.emc_alloc_group_parts(cp, gf, rep(1, q))
  expect_equal(EMC2:::.emc_alloc_group_parts(cp, gf, rep(1, q), numeric(q)), plain,
               tolerance = 1e-12)
  ex <- c(2, 0, -0.3, 0.5, 0, 0, 0, 0)
  with <- EMC2:::.emc_alloc_group_parts(cp, gf, rep(1, q), ex)
  expect_equal(tau(with$levels), pmax(1, tau(plain$levels) * (1 + ex)), tolerance = 1e-8)
  expect_equal(tau(with$best), pmax(1, tau(plain$best) * (1 + ex)), tolerance = 1e-8)
  rs <- seq(max(with$best[movable]) + 0.002, max(with$levels[movable]), length.out = 12)
  nd <- EMC2:::.emc_alloc_group_need(cp, gf, rep(1, q), movable, rs, ex)
  X <- do.call(rbind, lapply(nd, `[[`, "x"))
  cost <- Reduce(`+`, lapply(nd, `[[`, "cost"))
  expect_true(any(is.finite(cost)))
  for (t in which(is.finite(cost))) {
    r <- gf$b + (1 - gf$b) * EMC2:::.emc_alloc_group_rho(cp, g, g, X[, t])[, 1]
    expect_true(all((pmax(1, tau(r) * (1 + ex)))[movable] <= tau(rs[t]) * (1 + 1e-8)))
  }
})

test_that("a measured time is shrunk towards the system's by the blocks behind it", {
  set.seed(3)
  ar <- function(n, rho) as.numeric(stats::filter(rnorm(n), rho, "recursive"))
  ys <- function(n) lapply(1:3, function(i) rbind(ar(n, 0.95), ar(n, 0.5)))
  # A time of 39 against a system that says 10, and one of 3 against 3.
  long <- EMC2:::.emc_alloc_excess(ys(4000), c(10, 3))
  expect_gt(long[1], 1.5)
  expect_lt(long[1], 3.5)
  expect_lt(abs(long[2]), 0.85)
  short <- EMC2:::.emc_alloc_excess(ys(150), c(10, 3))
  expect_lt(abs(short[1]), long[1] / 2)
  expect_identical(EMC2:::.emc_alloc_excess(lapply(1:3, function(i) matrix(rnorm(8), 2)),
                                            c(10, 3)), c(0, 0))
  expect_identical(EMC2:::.emc_alloc_excess(ys(500), c(Inf, NA)), c(0, 0))
})

test_that("holds on the subjects' deviations scale the group's excess", {
  subject <- function(h, base = 8) {
    u <- list(kind = "subject", rho = c(0.5, 0.6, 0.8, 0.85), b_share = rep(0, 4),
              p = rep(1, 4), c_p = 0, c_h = 1, n_dev = 2L, tau_base = base)
    u["h"] <- list(h)
    u
  }
  group <- list(kind = "shear", rho = 0.9, g = 0.5, cost = 1)
  units <- list(group, subject(rep(0.3, 4)), subject(NULL), subject(rep(0.3, 4), NA))
  movable <- rep(TRUE, 4)
  expect_identical(EMC2:::.emc_alloc_excess_holds(units, NULL, movable),
                   list(list(units = units, ex = NULL, hold = Inf)))
  holds <- EMC2:::.emc_alloc_excess_holds(units, c(30, -2), movable)
  tau <- function(r) (1 + r) / (1 - r)
  # Nothing held: a subject with a measured base at its particle step's time
  # over the base, one without at one.
  level <- mean(tau(c(0.8, 0.85))) / 8
  expect_equal(holds[[1]]$ex, c(30 * (2 * level + 1) / 3, -2))
  expect_identical(holds[[1]]$units, units)
  expect_gt(length(holds), 1L)
  for (hd in holds[-1]) {
    # Only the subject that can reach the hold is limited, on its deviations.
    expect_null(hd$units[[3]]$lim)
    expect_null(hd$units[[4]]$lim)
    lim <- hd$units[[2]]$lim
    expect_identical(is.na(lim), c(TRUE, TRUE, FALSE, FALSE))
    expect_equal(tau(lim[3:4]), rep(hd$hold * 8, 2))
    expect_equal(hd$ex, c(30 * (min(hd$hold, level) + level + 1) / 3, -2))
    # The limited subject's need reaches the limit at every target.
    rs <- c(0.5, 0.9, 0.99)
    nd <- EMC2:::.emc_alloc_units_need(hd$units[2], TRUE, rs, 4)[[1]]
    expect_true(all(is.finite(nd$cost)))
    expect_gt(nd$x[3], 0)
    expect_lte(max(tau(c(0.8, 0.85) * (1 - nd$x[3] * 0.7))), hd$hold * 8 * (1 + 1e-8))
    expect_equal(EMC2:::.emc_alloc_unit_need(hd$units[[2]], rs, 4)$x, nd$x)
  }
  expect_true(all(diff(vapply(holds, `[[`, numeric(1), "hold")) < 0))
})

test_that("the pooled blocks' plans accumulate only where the plan recorded them", {
  ch <- list(n_pars = 2L, n_subjects = 3L,
             alloc = list(last = list(rounds = c(1, 0, 0.5, NA), collapsed_rounds = 0.2,
                                      dev_tau = c(4, NA, 6))))
  first <- EMC2:::.emc_alloc_pool_plans(NULL, ch, list(n = 100L))
  expect_equal(first, list(n = 100, x = c(100, 0, 50, 0), xc = 20, dev = c(400, 0, 600),
                           dw = c(100, 0, 100)))
  again <- EMC2:::.emc_alloc_pool_plans(c(list(G = 1, nG = 1), first), ch, list(n = 50L))
  expect_equal(again$n, 150)
  expect_equal(again$dev / pmax(again$dw, 1), c(4, 0, 6))
  ch$alloc$last$dev_tau <- NULL
  expect_identical(EMC2:::.emc_alloc_pool_plans(first, ch, list(n = 100L)), first)
})
