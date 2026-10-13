local_rng_guard()  # see helper-rng.R: keep this file's RNG changes inside it

# Each subject's particle count is chosen by the move allocation
# (R/move_allocation.R) from every step's own clouds: their expected squared
# jump at fractions of their size (new_particle) and the step's cost.

# A cloud of n iid draws from proposal N(0, s^2 I) for target N(0, I), with
# its reference drawn from the target: log weights and squared jumps.
gauss_cloud <- function(n, p, s, x0 = stats::rnorm(p)) {
  x <- matrix(stats::rnorm(n * p, sd = s), n)
  lw <- function(z) rowSums(stats::dnorm(z, log = TRUE)) -
    rowSums(stats::dnorm(z, sd = s, log = TRUE))
  l <- c(lw(matrix(x0, 1L)), lw(x))
  w <- exp(l - max(l))
  list(w0 = w[1L], w = w[-1L], d = rowSums(sweep(x, 2L, x0)^2))
}

test_that("a cloud's leading draws measure the jump of a smaller cloud", {
  local_test_rng(1)
  frac <- EMC2:::.EMC_COUNT_FRACTIONS
  reps <- 3000L
  # Two components of 40 and 24 draws; at fraction f each keeps its first
  # round(f * n_k), so the subset is a cloud of the smaller size.
  sizes <- c(40L, 24L)
  prefix <- full <- matrix(0, reps, length(frac))
  for (r in seq_len(reps)) {
    x0 <- stats::rnorm(6)
    cl <- gauss_cloud(sum(sizes), 6, 1.6, x0)
    prefix[r, ] <- EMC2:::.emc_count_prefix_curve(cl$w0, cl$w, cl$d, sizes)
    for (k in seq_along(frac)) {
      small <- gauss_cloud(sum(pmax(1L, round(frac[k] * sizes))), 6, 1.6, x0)
      full[r, k] <- sum(small$w * small$d) / (small$w0 + sum(small$w))
    }
  }
  se <- sqrt(apply(prefix, 2L, stats::var) / reps + apply(full, 2L, stats::var) / reps)
  expect_true(all(abs(colMeans(prefix) - colMeans(full)) < 4 * se))
  # The jump grows with the cloud here, so the measured curve rises.
  expect_true(all(diff(colMeans(prefix)) > 0))
})

test_that("measuring the curve leaves the particle step's draws unchanged", {
  local_mocked_bindings(calc_ll_pooled = function(proposals, ...)
    -0.5 * rowSums(proposals^2), .package = "EMC2")
  p <- 3L
  alpha <- stats::setNames(c(0.1, -0.2, 0.3), letters[1:3])
  pm <- list(list(epsilon = rep(0.5, 3L), mix = EMC2:::get_default_mix("adapt"),
                  n_particles = 60L, iter = 0L, proposal_counts = integer(4L),
                  acc_counts = integer(4L), gd_good = FALSE))
  run <- function(evidence) {
    set.seed(5L)
    out <- EMC2:::new_particle(
      1L, NULL, pm, chains_mu = alpha, chains_var = diag(p), eff_mu = alpha,
      eff_var = diag(p), prev_ll = 0, parameters = NULL, stage = "sample",
      type = "standard",
      tune = list(components = rep(1L, p), shared_ll_idx = rep(1L, p),
                  count_evidence = evidence),
      current_alpha = alpha, population_mu = alpha * 0, population_var = diag(p))
    list(out = out, seed = .Random.seed)
  }
  with_ev <- run(TRUE)
  without <- run(FALSE)
  expect_identical(with_ev$out$proposal, without$out$proposal)
  expect_identical(with_ev$seed, without$seed)
  ev <- with_ev$out$pm_settings[[1L]]$count_ev
  expect_length(ev$curve, length(EMC2:::.EMC_COUNT_FRACTIONS))
  expect_true(all(ev$curve >= 0))
  # One likelihood call per conditional-IS cloud.
  expect_length(ev$rows, 2L)
  expect_null(without$out$pm_settings[[1L]]$count_ev)
})

test_that("kappa interpolates the measured curve and extends it at bounded slopes", {
  frac <- EMC2:::.EMC_COUNT_FRACTIONS
  curve <- 10 * frac^0.4
  kap <- function(f) EMC2:::.emc_alloc_kappa(curve, f)
  expect_equal(vapply(frac, kap, 0), frac^0.4)
  expect_equal(kap(0.6), 0.6^0.4)
  expect_equal(kap(sqrt(2)), sqrt(2)^0.4)
  expect_equal(kap(0.25), 0.25^0.4)
  # A curve rising faster than the count is extended at slope one.
  expect_equal(EMC2:::.emc_alloc_kappa(10 * frac^3, 2), 2)
  # A falling curve is extended flat.
  expect_equal(EMC2:::.emc_alloc_kappa(c(3, 2, 1), 2), 1)
})

count_state <- function(curve, m_meas = 1, base = 100, row_cpu = NA_real_) {
  list(base = list(base), m_now = m_meas, m_meas = m_meas,
       curve = matrix(curve, 1L), row_cpu = row_cpu)
}

test_that("candidate counts carry the scaled factors, their cost and the floor", {
  frac <- EMC2:::.EMC_COUNT_FRACTIONS
  u <- list(kind = "subject", rho = c(0.9, 0.8), p = c(0.4, 0.6), c_p = 0.1,
            p_dev = c(0.5, 0.7))
  # Jump proportional to the square root of the count; each data row of a
  # likelihood row costs 1e-5 seconds over a fixed cost.
  st <- count_state(10 * sqrt(frac), m_meas = 0.5, base = 200, row_cpu = 1e-5)
  v <- EMC2:::.emc_alloc_count_candidates(u, st, 1L, data_rows = 50)
  expect_equal(v$cands$m[1L], 0.5)
  expect_equal(sort(v$cands$m), sort(c(0.5 * c(1, sqrt(0.5), 0.5, sqrt(2)), 1)))
  for (k in seq_along(v$cands$m)) {
    kap <- sqrt(v$cands$m[k] / 0.5)
    expect_equal(v$cands$p[[k]], 1 - kap * (1 - u$p))
    expect_equal(v$cands$p_dev[[k]], 1 - kap * (1 - u$p_dev))
    rows <- max(25, round(200 * v$cands$m[k]))
    expect_equal(v$cands$c[k], 0.1 + 1e-5 * 50 * (rows - 100))
  }
  # The floor is the step at its given count (multiplier one): the jump
  # there extends the measured curve beyond its count.
  expect_equal(v$p_level, 1 - sqrt(2) * (1 - u$p))
  # Counts stay within [min(25, base), 2 base].
  st <- count_state(10 * sqrt(frac), m_meas = 0.2, base = 100)
  v <- EMC2:::.emc_alloc_count_candidates(u, st, 1L, data_rows = 50)
  expect_gte(min(v$cands$m), 0.25 - 1e-12)
  st <- count_state(10 * sqrt(frac), m_meas = 2, base = 100)
  expect_lte(max(EMC2:::.emc_alloc_count_candidates(u, st, 1L, 50)$cands$m), 2)
  # Without a measured curve there are no candidates.
  st <- count_state(rep(NA_real_, 3))
  expect_null(EMC2:::.emc_alloc_count_candidates(u, st, 1L, 50)$cands)
})

test_that("a subject cut far below its given count can still meet its floor", {
  frac <- EMC2:::.EMC_COUNT_FRACTIONS
  u <- list(kind = "subject", rho = c(0.9, 0.8), p = c(0.7, 0.8), c_p = 0.1,
            p_dev = c(0.8, 0.9))
  st <- count_state(10 * sqrt(frac), m_meas = 0.25, base = 400)
  v <- EMC2:::.emc_alloc_count_candidates(u, st, 1L, data_rows = 50)
  expect_true(any(abs(v$cands$m - 1) < 1e-12))
  # Without HMC only the given count's step reaches the floor, and it does.
  floor <- max(EMC2:::.emc_alloc_subject_floor(v))
  rs <- sort(c(floor, 0.95, 0.99))
  need <- EMC2:::.emc_alloc_units_need(list(v), TRUE, rs, 4)[[1L]]
  expect_true(all(is.finite(need$cost)))
  expect_equal(need$m[1L], 1)
})

test_that("each subject takes its cheapest count at every target", {
  local_test_rng(7)
  rs <- sort(c(stats::runif(30, 0.05, 1), 1))
  for (trial in 1:20) {
    k <- sample(1:3, 1)
    u <- list(kind = "subject", rho = stats::runif(k, 0.6, 0.999),
              p = stats::runif(k, 0.2, 0.9), c_p = stats::runif(1, 0.05, 1),
              h = if (stats::runif(1) < 0.6) stats::runif(k, 0.05, 0.9),
              c_h = stats::runif(1, 0.05, 1))
    frac <- EMC2:::.EMC_COUNT_FRACTIONS
    st <- count_state(10 * frac^stats::runif(1, 0.05, 0.9), m_meas = 1,
                      row_cpu = stats::runif(1, 1e-6, 1e-4))
    u <- EMC2:::.emc_alloc_count_candidates(u, st, 1L, data_rows = 100)
    got <- EMC2:::.emc_alloc_units_need(list(u), TRUE, rs, 4)[[1L]]
    # Each candidate priced alone as its own step, against the same floor
    # (the step at its given count).
    cost <- do.call(rbind, lapply(seq_along(u$cands$m), function(c) {
      v <- u
      v$p <- u$cands$p[[c]]
      v$c_p <- u$cands$c[c]
      v$cands <- NULL
      EMC2:::.emc_alloc_unit_need(v, rs, 4)$cost
    }))
    expect_equal(got$cost, apply(cost, 2L, min), tolerance = 1e-12)
    ok <- is.finite(got$cost)
    chosen <- match(got$m[ok], u$cands$m)
    expect_equal(cost[cbind(chosen, which(ok))], got$cost[ok], tolerance = 1e-12)
    expect_false(anyNA(chosen))
  }
})

test_that("a cheaper count that still meets the target replaces a dearer one", {
  frac <- EMC2:::.EMC_COUNT_FRACTIONS
  fixed <- list(kind = "fixed", rho = 0.2, cost = 0)
  # The step barely moves less at half the count, and costs half as much.
  u <- list(kind = "subject", rho = c(0.95, 0.9), p = c(0.3, 0.3), c_p = 1)
  u <- EMC2:::.emc_alloc_count_candidates(
    u, count_state(c(0.97, 0.99, 1), m_meas = 1, base = 200, row_cpu = NA_real_),
    1L, data_rows = 10)
  # The floor (the step at its given count) is held by a little HMC, too
  # weak to replace the particle step but cheap enough to make up the gap.
  with_h <- c(u, list(h = c(0.9, 0.9), c_h = 0.05))
  sol <- EMC2:::.emc_alloc_solve(list(fixed, with_h), 1)
  expect_identical(sol$q[2], 1)
  expect_lt(sol$m[2], 1)
  # Without HMC the floor rules out a weaker step: the count holds.
  sol <- EMC2:::.emc_alloc_solve(list(fixed, u), 1)
  m <- if (is.null(sol$m)) NA_real_ else sol$m[2]
  expect_true(is.na(m) || m >= 1)
})

test_that("applied counts follow the multiplier within bounds and stand the ESS rule down", {
  ch <- list(samples = structure(list(), pm_settings = list(
    list(list(n_particles = 100), list(n_particles = 40)),
    list(list(n_particles = 100), list(n_particles = 40)))),
    alloc = list(count = list(base = list(c(100, 40), c(100, 40)), m = c(0.5, 3))))
  out <- EMC2:::.emc_alloc_apply_counts(ch)
  pm <- attr(out$samples, "pm_settings")
  expect_equal(vapply(pm[[1]], `[[`, 0, "n_particles"), c(50, 25))
  expect_equal(vapply(pm[[2]], `[[`, 0, "n_particles"), c(200, 80))
  expect_true(all(vapply(unlist(pm, recursive = FALSE), function(b)
    isTRUE(b$count_managed), logical(1))))
  # The weight-ESS rule leaves a managed count alone.
  pmb <- list(iter = 1000L, proposal_counts = rep(10, 4), acc_counts = rep(3, 4),
              mix = rep(0.25, 4), epsilon = rep(1, 3), n_particles = 100,
              gd_good = TRUE)
  tune <- list(n0 = 10L, p_accept = rep(0.7, 3), mix_adapt = 0.1,
               target_ESS = 30, ESS_scale = 0.5, max_particles = 500)
  w <- c(0.5, rep(0.01, 40))
  upd <- function(b) EMC2:::update_pm_settings(b, rep(1, 4), rep(10, 4),
                                               sum(w)^2 / sum(w^2), tune, n_pars = 5L)
  expect_false(upd(pmb)$n_particles == 100)
  expect_identical(upd(c(pmb, count_managed = TRUE))$n_particles, 100)
})
