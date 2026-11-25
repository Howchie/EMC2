# Allocation of shear rounds and collapsed moves (R/move_allocation.R).
local_rng_guard()

# rho of a unit under allocation x, by the composition the solver assumes.
alloc_unit_rho <- function(u, x) {
  switch(u$kind,
    fixed = max(u$rho),
    shear = u$rho * EMC2:::.emc_alloc_factor(u$g, x, if (is.null(u$a)) u$g else u$a))
}

test_that("the allocation is as good as a dense grid over rounds", {
  local_test_rng(7)
  tau <- function(r) (1 + r) / (1 - r)
  for (trial in 1:30) {
    shear <- function() {
      g <- stats::runif(1, 0.05, 0.95)
      # Half the trials saturate: a below g.
      a <- if (trial %% 2 == 0) stats::runif(1, 0, g) else g
      list(kind = "shear", rho = stats::runif(1, 0, 0.99), g = g, a = a,
           cost = stats::runif(1, 0.02, 0.3))
    }
    units <- list(shear(), shear(), list(kind = "fixed", rho = stats::runif(1, 0, 0.9), cost = 0))
    C <- 1
    xs <- seq(0, 4, by = 0.04)
    best <- Inf
    for (x1 in xs) for (x2 in xs) {
      # The fixed coordinate counts for nothing: no move lowers it.
      rho <- max(alloc_unit_rho(units[[1]], x1), alloc_unit_rho(units[[2]], x2))
      best <- min(best, tau(rho) * (C + x1 * units[[1]]$cost + x2 * units[[2]]$cost))
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
    expect_true(all(sol$x[1:2] >= 0 & sol$x[1:2] <= EMC2:::.EMC_ALLOC_MAX_ROUNDS))
  }
})

test_that("a coordinate no move can lower does not bound the others", {
  solve <- EMC2:::.emc_alloc_solve
  # Every coordinate fixed: nothing to buy.
  expect_equal(solve(list(list(kind = "fixed", rho = 0.9, cost = 0)), 1)$x, 0)
  # A slow coordinate no shear can reach does not stop the rounds on another:
  # the allocation is over what the moves can reach.
  sol <- solve(list(list(kind = "fixed", rho = 0.99, cost = 0),
                    list(kind = "shear", rho = 0.9, g = 0.4, cost = 0.01)), 1)
  expect_equal(sol$x[1], 0)
  expect_gt(sol$x[2], 0)
  expect_equal(sol$tau0, EMC2:::.emc_alloc_tau(0.9))
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
  s <- EMC2:::.emc_alloc_record_cost(s, 0.5, moves = list(shear = c(cpu = 0.1, own = 0.05, n = 2)))
  w <- s$alloc$window
  expect_identical(c(w$n, w$nG, w$fam$shear$moves), c(1L, 1L, 2L))
  expect_equal(c(w$cpu_rest, w$fam$shear$cpu, w$fam$shear$own), c(0.5, 0.1, 0.05))
})

test_that("the allocation preserves the between-chain rho floor", {
  rho <- c(0.2, 0.8)
  b <- c(0.5, 0.5)
  unit <- list(kind = "shear", rho = rho, b_share = b, g = c(0.1, 0.2),
               cost = 0, cap = 100)
  predicted <- EMC2:::.emc_alloc_predict(list(unit), list(x = 100, joint = NULL),
                                         max_rounds = 100)[[1]]
  expect_true(all(predicted >= b - 1e-12))
  expect_equal(predicted, b + (1 - b) * rho * unit$g^100, tolerance = 1e-12)
})

test_that("a chain held in one region measures slow against the pooled spread", {
  local_test_rng(3)
  p <- 1; n <- 1; m <- 200
  chain <- function(centre, sd) {
    a <- centre + sd * stats::rnorm(m)
    list(n_pars = p, n_subjects = n, type = "diagonal-gamma", par_names = "x",
         samples = list(idx = m, theta_mu = matrix(a, 1),
                        theta_var = array(1, c(1, 1, m)), alpha = array(a, c(1, 1, m))),
         alloc = list(window = modifyList(EMC2:::.emc_alloc_empty_window(p),
                                          list(n = m - 1L, nG = m - 1L, cpu_rest = 1,
                                               G = diag(c(sum(diff(a)^2), 0))))))
  }
  # Two chains mixing over the same spread look fast; with one of them held
  # at a distance, the held one's small jumps are slow against the pooled sd.
  mixed <- list(chain(0, 1), chain(0, 1))
  held <- list(chain(0, 1), chain(4, 0.05))
  rho <- function(emc, j) EMC2:::.emc_alloc_measure(emc, j)$g_rho[1]
  expect_lt(rho(mixed, 1), 0.2)
  expect_gt(rho(held, 2), 0.99)
  # A stuck coordinate reads one; a jump of twice the variance reads zero.
  expect_equal(EMC2:::.emc_alloc_rho(0, 1), 1)
  expect_equal(EMC2:::.emc_alloc_rho(2, 1), 0)
  expect_equal(EMC2:::.emc_alloc_rho(0.5, 1), 0.75)
})

test_that("production measurements use every production draw and pool the stage's jumps", {
  local_test_rng(5)
  p <- 1; n <- 1; burn <- 100; m <- 400; w_n <- 99L
  chain <- function(centre) {
    a <- c(stats::rnorm(burn, centre + 5), centre + stats::rnorm(m))
    k <- burn + m
    list(n_pars = p, n_subjects = n, type = "diagonal-gamma", par_names = "x",
         samples = list(idx = k, stage = rep(c("burn", "sample"), c(burn, m)),
                        theta_mu = matrix(a, 1), theta_var = array(1, c(1, 1, k)),
                        alpha = array(a, c(1, 1, k))),
         alloc = list(stage = "sample", pool = list(G = diag(c(30, 0)), nG = 50L),
                      window = modifyList(EMC2:::.emc_alloc_empty_window(p),
                                          list(n = w_n, nG = w_n, cpu_rest = 1,
                                               G = diag(c(20, 0))))))
  }
  emc <- list(chain(0), chain(0.5))
  prod <- EMC2:::.emc_alloc_measure(emc, 1, "sample")
  x <- c(emc[[1]]$samples$theta_mu[burn + seq_len(m)],
         emc[[2]]$samples$theta_mu[burn + seq_len(m)])
  within_g <- mean(c(stats::var(emc[[1]]$samples$theta_mu[burn + seq_len(m)]),
                     stats::var(emc[[2]]$samples$theta_mu[burn + seq_len(m)])))
  expect_equal(prod$g_var[1], within_g)
  expect_identical(prod$pairs, w_n + 50L)
  expect_equal(prod$pool$G, diag(c(50, 0)))
  expect_equal(prod$g_rho[1], min(1, EMC2:::.emc_alloc_rho(50 / (w_n + 50), within_g) /
                                       (1 - prod$g_white[1])))
  expect_equal(prod$g_b_share[1], max(0, stats::var(x) - within_g) / stats::var(x))
  # The first production rebuild measures the window alone and starts the pool.
  first <- emc
  first[[1]]$alloc$stage <- "adapt"
  m1 <- EMC2:::.emc_alloc_measure(first, 1, "sample")
  expect_identical(m1$pairs, w_n)
  expect_null(m1$pool)
  # Outside production: the window's draws, the window's jumps.
  b <- EMC2:::.emc_alloc_measure(emc, 1, "burn")
  k <- EMC2:::.emc_alloc_draws(emc[[1]], w_n)
  within_b <- mean(c(stats::var(emc[[1]]$samples$theta_mu[k]),
                     stats::var(emc[[2]]$samples$theta_mu[k])))
  expect_equal(b$g_var[1], within_b)
  expect_identical(b$pairs, w_n)
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
             samples = list(idx = k, theta_mu = matrix(mu, 1), theta_var = array(Sig, c(1, 1, k)),
                            alpha = array(0, c(1, N, k))),
             shear = list(plan = list(dirs = diag(1),
                                      last = list(E = c(0.2, 0), n = c(10L, 0L), moved = integer(2)))),
             alloc = list(window = modifyList(EMC2:::.emc_alloc_empty_window(p),
                                              list(n = k - 1L, nG = k - 1L,
                                                   G = diag(c(sum(diff(mu)^2), 0))))))
  local_mocked_bindings(.emc_shear_active = function(sampler) TRUE,
                        .emc_alloc_draws = function(ch, m) seq_len(ch$samples$idx),
                        .package = "EMC2")
  ms <- EMC2:::.emc_alloc_measure(list(ch), 1)
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
             samples = list(idx = k, theta_mu = a, theta_var = array(diag(2), c(p, p, k)),
                            alpha = array(a, c(p, 1, k))),
             shear = list(plan = list(dirs = diag(p),
                                      last = list(E = c(10, 0, 0, 0), n = c(10L, 0L, 0L, 0L),
                                                  moved = integer(4)))),
             alloc = list(window = NULL))
  local_mocked_bindings(.emc_shear_active = function(sampler) TRUE, .package = "EMC2")
  ms <- EMC2:::.emc_alloc_measure(list(ch), 1)
  g1 <- max(0, 1 - 1 / (2 * ms$g_var[1] * (1 - ms$g_white[1])))
  expect_equal(ms$g_fac[1], g1)
  # The other mean direction borrows the measured one; the log sds, with no
  # measured round of their own, borrow it too.
  expect_equal(ms$g_fac[2:4], rep(g1, 3))
  expect_identical(ms$g_measured, c(TRUE, FALSE, FALSE, FALSE))
  # Without any measured round, the prior factor.
  ch$shear$plan$last$n[] <- 0L
  ms0 <- EMC2:::.emc_alloc_measure(list(ch), 1)
  expect_equal(ms0$g_fac, rep(EMC2:::.EMC_ALLOC_SHEAR_PRIOR_G, 4))
})

test_that("measured zero movement is a factor of exactly one, and buys nothing", {
  p <- 2; n <- 1; k <- 50
  a <- matrix(stats::rnorm(p * k), p)
  ch <- list(n_pars = p, n_subjects = n, type = "diagonal-gamma", par_names = c("x", "y"),
             samples = list(idx = k, theta_mu = a, theta_var = array(diag(2), c(p, p, k)),
                            alpha = array(a, c(p, 1, k))),
             shear = list(plan = list(dirs = diag(p),
                                      last = list(E = c(0, 0, 0, 0), n = c(10L, 0L, 0L, 0L),
                                                  moved = integer(4)))),
             alloc = list(window = NULL))
  local_mocked_bindings(.emc_shear_active = function(sampler) TRUE, .package = "EMC2")
  ms <- EMC2:::.emc_alloc_measure(list(ch), 1)
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
             samples = list(idx = k, theta_mu = a, theta_var = v,
                            alpha = array(a, c(p, 1, k))),
             collapsed = list(plan = list(C = diag(2 * p), rounds = 1,
                                          last = list(E = matrix(0, 2 * p, 2 * p), n = 600,
                                                      acc = 0, moved = 0L))),
             alloc = list(window = NULL))
  local_mocked_bindings(.emc_collapsed_active = function(sampler) TRUE, .package = "EMC2")
  ms <- EMC2:::.emc_alloc_measure(list(ch), 1)
  expect_true(ms$c_measured)
  expect_identical(ms$c_fac, rep(1, 2 * p))
  # Movement on one coordinate is shaded towards the move's own mean factor.
  ch$collapsed$plan$last$E[1, 1] <- 60
  ms <- EMC2:::.emc_alloc_measure(list(ch), 1)
  obs <- c(max(0, 1 - 0.1 / (2 * ms$g_var[1] * (1 - ms$g_white[1]))), rep(1, 3))
  expect_equal(ms$c_fac, (600 * obs + 20 * mean(obs)) / 620)
})

test_that("a move earns its evidence at an exploration rate and runs no faster than its evidence allows", {
  n0 <- EMC2:::.EMC_ALLOC_ROUND_EVIDENCE
  B <- 100
  ms <- list(g_rho = c(0.99, 0.6), g_fac = c(0.5, 0.9),
             g_measured = c(FALSE, TRUE), g_n = c(0, 400), C_rest = 1, c_round = 0.05,
             coll = TRUE, c_fac = rep(0.5, 2), c_measured = FALSE, c_n = 0, c_coll = 0.1,
             c_rate = 0)
  u <- EMC2:::.emc_alloc_units(ms, TRUE, B = B)
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
  uh <- EMC2:::.emc_alloc_units(half, TRUE, B = B)
  expect_equal(c(uh$units[[1]]$floor, uh$units[[1]]$cap), c(n0 / 2 / B, 1))
  expect_equal(uh$joint$floor, n0 / 2 / B)
  # Without B (no rebuild horizon) nothing is floored or capped by evidence.
  u0 <- EMC2:::.emc_alloc_units(ms, TRUE)
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
  ms <- list(g_rho = c(0.99, 0.6), g_fac = c(0.5, 0.9),
             g_measured = c(TRUE, TRUE), g_n = c(x * iters, 400), g_iters = c(iters, 0),
             C_rest = 1, c_round = 0.05, coll = TRUE, c_fac = rep(0.5, 2), c_measured = TRUE,
             c_n = x * iters, c_iters = iters, c_coll = 0.1, c_rate = 0)
  u <- EMC2:::.emc_alloc_units(ms, TRUE, B = 400)
  expect_equal(u$units[[1]]$cap, 2 * x)
  expect_equal(u$joint$cap, 2 * x)
  # Without the span the evidence is read against the next block.
  expect_equal(u$units[[2]]$cap, 2 * 400 / 400)
  ms$g_iters <- NULL
  expect_equal(EMC2:::.emc_alloc_units(ms, TRUE, B = 400)$units[[1]]$cap,
               2 * x * iters / 400)
  # A measured move that ran little may still run at one move per iteration:
  # a plan that ran none of it does not cap it down for good.
  low <- ms; low$g_iters <- c(iters, 0); low$g_n <- c(0.01 * iters, 400)
  low$c_n <- 0.01 * iters
  ul <- EMC2:::.emc_alloc_units(low, TRUE, B = 400)
  expect_equal(c(ul$units[[1]]$cap, ul$joint$cap), c(1, 1))
  # The span of a closing block is the allocation window's iteration count.
  expect_equal(EMC2:::.emc_alloc_block_iters(list(alloc = list(window = list(n = 37L)))), 37)
  expect_equal(EMC2:::.emc_alloc_block_iters(list()), 0)
})

test_that("a few rounds shade a coordinate's factor towards its kind's, and log-sd rounds pool across chains", {
  p <- 2; n <- 1; k <- 50
  a <- matrix(stats::rnorm(p * k), p)
  tv <- array(0, c(p, p, k)); for (t in seq_len(k)) tv[, , t] <- diag(exp(stats::rnorm(p)))
  chain <- function(last, stage = "sample") list(
    n_pars = p, n_subjects = n, type = "diagonal-gamma", par_names = c("x", "y"),
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
  ms <- EMC2:::.emc_alloc_measure(list(ch), 1)
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
  ms2 <- EMC2:::.emc_alloc_measure(list(ch2), 1)
  expect_gt(ms2$g_fac[2], 0.98)
  # A second chain's log-sd rounds count for this chain's log sds, not for
  # its mean directions.
  other <- chain(list(E = c(0, 0, 1, 0), n = c(0, 0, 10, 0), moved = integer(4)))
  mine <- chain(list(E = c(0, 0, 0, 0), n = c(0, 0, 0, 0), moved = integer(4)))
  msp <- EMC2:::.emc_alloc_measure(list(mine, other), 1, stage = "sample")
  g3 <- max(0, 1 - 0.1 / (2 * msp$g_var[3] * (1 - msp$g_white[3])))
  expect_true(msp$g_measured[3]); expect_false(msp$g_measured[1])
  expect_equal(msp$g_fac[3], (10 * g3 + n0 * g3) / (10 + n0))
  # Outside production nothing pools.
  msb <- EMC2:::.emc_alloc_measure(list(mine, other), 1, stage = "burn")
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
  for (i in 1:12) s <- EMC2:::.emc_alloc_record_cost(s, 1,
    moves = list(shear = c(cpu = 2.5, own = 0.5, n = 5), collapsed = c(cpu = 1, own = 0.8, n = 1)))
  s <- EMC2:::.emc_alloc_record_cost(s, 1, moves = list(shear = c(cpu = 0, own = 0, n = 0)))
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
  u <- list(units = list(sh, fixed), C_rest = 1, joint = NULL, shear_fixed = 0.3)
  sol <- EMC2:::.emc_alloc_solve_units(u)
  expect_gt(sol$x[1], 0)
  expect_equal(sol$cost, 1 + 0.3 + sol$x[1] * 0.02)
  expect_equal(sol$tau0, EMC2:::.emc_alloc_tau(0.95))
  dear <- EMC2:::.emc_alloc_solve_units(replace(u, "shear_fixed", list(50)))
  expect_identical(dear$x[1], 0)
  expect_equal(dear$cost, 1)
})

test_that("an unmeasured shear coordinate counts at the prior factor, uncapped", {
  ms <- list(g_rho = c(0.976, 0.6), g_fac = c(0.5, 0.9),
             g_measured = c(FALSE, TRUE), C_rest = 1, c_round = 0.05, coll = TRUE,
             c_fac = c(0.965, 0.5), c_measured = TRUE, c_coll = 0.1, c_rate = 1)
  u <- EMC2:::.emc_alloc_units(ms, TRUE)
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
           c_coll = 0.1, c_rate = 1, stored = 0L, pairs = 0L,
           C_rest = 1, c_round = 0.1, pool = NULL)
    },
    .emc_alloc_units = function(ms, ...) {
      list(chain_id = ms$chain_id, units = list(), C_rest = 1,
           c_round = 0.1, joint = list(cost = 0.1))
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

test_that("the shear's fixed cost is charged at the plan's activation probability", {
  sh <- list(kind = "shear", rho = 0.9, g = 0.5, cost = 0.05)
  fixed <- list(kind = "fixed", rho = 0.9, cost = 0)
  act <- EMC2:::.emc_alloc_activation
  # independent fractional rounds; a whole round always activates
  expect_equal(act(list(sh, sh), cbind(c(0.1, 0.2), c(1.5, 0)), 10), 10 * c(1 - 0.9 * 0.8, 1))
  # only shear units activate it
  expect_equal(act(list(fixed, sh), cbind(c(1, 0.5)), 4), 4 * 0.5)
  expect_equal(act(list(fixed), cbind(1), 4), 0)
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
  shear <- function() list(kind = "shear", rho = stats::runif(1, 0.5, 0.99),
                           g = stats::runif(1, 0.2, 0.95), cost = stats::runif(1, 0.01, 0.2))
  for (trial in 1:15) {
    units <- list(shear(), list(kind = "fixed", rho = stats::runif(1, 0, 0.5), cost = 0), shear())
    sol <- EMC2:::.emc_alloc_solve(units, 1)
    pred <- sol$pred
    expect_length(pred, length(units))
    for (i in c(1, 3))
      expect_equal(pred[[i]], units[[i]]$rho *
                     EMC2:::.emc_alloc_factor(units[[i]]$g, sol$x[i], units[[i]]$g))
    expect_equal(pred[[2]], units[[2]]$rho)
    # The solver's target bounds every coordinate it ranges over.
    expect_lte(max(unlist(pred[-2])), sol$rho + 1e-9)
  }
})

test_that("the solver's objective is the model's objective at the plan it returns", {
  local_test_rng(23)
  shear <- function() list(kind = "shear", rho = stats::runif(1, 0.5, 0.99),
                           g = stats::runif(1, 0.2, 0.95), cost = stats::runif(1, 0.01, 0.2))
  tau <- EMC2:::.emc_alloc_tau
  for (trial in 1:25) {
    u <- list(units = list(shear(), shear()), C_rest = 1, joint = NULL,
              shear_fixed = stats::runif(1, 0, 0.2))
    sol <- EMC2:::.emc_alloc_solve_units(u)
    rho <- max(vapply(1:2, function(i) {
      g <- u$units[[i]]
      g$rho * EMC2:::.emc_alloc_factor(g$g, sol$x[i], g$g)
    }, numeric(1)))
    cost <- 1 + sum(sol$x * vapply(u$units, `[[`, 0, "cost")) +
      u$shear_fixed * (1 - prod(1 - pmin(sol$x, 1)))
    # The solver reports its target, which the rounds meet: never below the
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
         samples = list(idx = k, stage = rep("sample", k), theta_mu = mu, theta_var = tv, alpha = a),
         shear = list(plan = list(dirs = matrix(stats::rnorm(p * p), p) + diag(2, p),
                                  stage = "sample",
                                  last = list(E = numeric(2 * p), n = numeric(2 * p),
                                              moved = integer(2 * p)))),
         alloc = list(window = NULL))
  }
  emc <- list(chain(), chain())
  emc[[2]]$shear$plan$dirs <- emc[[1]]$shear$plan$dirs
  local_mocked_bindings(.emc_shear_active = function(sampler) TRUE, .package = "EMC2")
  cache <- new.env()
  for (j in 1:2) {
    shared <- EMC2:::.emc_alloc_measure(emc, j, "sample", cache = cache)
    alone <- EMC2:::.emc_alloc_measure(emc, j, "sample")
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
})

test_that("a curved shear coordinate runs a round in every iteration", {
  ms <- list(g_rho = c(0.3, 0.98), g_fac = c(0.5, 0.95),
             g_n = c(200, 200), g_iters = c(100, 100), g_measured = c(TRUE, TRUE),
             C_rest = 1, c_round = 0.2, coll = FALSE)
  plain <- EMC2:::.emc_alloc_units(ms, TRUE, B = 100)
  expect_identical(plain$units[[1]]$floor, 0)
  curved <- EMC2:::.emc_alloc_units(c(ms, list(g_curved = 1L)), TRUE, B = 100)
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
    cp$C, cp$Q, cp$Ar, rep(1, q), g, a, gf$floor, gf$cap, gf$cost, gf$b, movable, rs, rebuild)
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
  units <- lapply(seq_len(q), function(k)
    list(kind = "shear", rho = cp$Ar[k, k], b_share = 0, g = g[k], a = g[k], cost = 0.05,
         floor = 0, cap = 4))
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

