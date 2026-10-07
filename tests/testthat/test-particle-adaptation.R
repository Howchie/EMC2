# C17: the quantity the particle count adapts from, recorded so it can be
# compared with the quantity that matters.
#
# `update_pm_settings` sets the particle count from
# `sum(w)^2 / sum(w^2)` -- the effective size of ONE proposal cloud. That is not
# the effective sample size of the saved chain, which is about autocorrelation
# across iterations and is what a posterior is worth. The audit's experiment is
# whether the first tracks the second; this makes both observable. The rule
# itself is unchanged, because changing it needs the PMwG invariant argument and
# frozen-versus-adaptive validation.

# Shapes as the sampler has them: one mixture weight per proposal, the first
# being the group-level proposal, and one epsilon and target acceptance for each
# of the others.
pm_stub <- function(n_mix = 4L, gd_good = TRUE, n_particles = 100) {
  list(iter = 1000L, proposal_counts = rep(10, n_mix), acc_counts = rep(3, n_mix),
       mix = rep(1 / n_mix, n_mix), epsilon = rep(1, n_mix - 1L),
       n_particles = n_particles, gd_good = gd_good)
}
tune_stub <- function() {
  list(n0 = 10L, p_accept = rep(0.7, 3), mix_adapt = 0.1,
       target_ESS = 30, ESS_scale = 0.5, max_particles = 500)
}
# `update_pm_settings` takes, per component, how many of its particles
# outweighed the reference particle, plus the weight ESS. For a single
# conditional-IS step over weights `w` (reference first, then each
# component's particles in turn) those are:
ups <- function(pm, w, particle_numbers, tune, n_pars) {
  ends <- 1L + cumsum(particle_numbers)
  starts <- ends - particle_numbers + 1L
  better <- mapply(function(a, b) if (b >= a) sum(w[a:b] > w[1]) else 0,
                   starts, ends)
  EMC2:::update_pm_settings(pm, better, particle_numbers,
                            sum(w)^2 / sum(w^2), tune, n_pars = n_pars)
}

test_that("the weight ESS is recorded whether or not the rule fires", {
  w <- c(0.5, runif(40))
  fired <- ups(pm_stub(gd_good = TRUE), w, rep(10, 4L),
                                     tune_stub(), n_pars = 5L)
  idle <- ups(pm_stub(gd_good = FALSE), w, rep(10, 4L),
                                    tune_stub(), n_pars = 5L)
  expect_false(is.null(fired$weight_ess))
  expect_false(is.null(idle$weight_ess))
  # It is the quantity the rule uses, exactly.
  expect_equal(fired$weight_ess, sum(w)^2 / sum(w^2))
  expect_equal(idle$weight_ess, sum(w)^2 / sum(w^2))
})

test_that("the weight ESS follows the current cloud, not the first one", {
  # The settings carry over from one iteration to the next, so the value must
  # be replaced on every call -- including calls where the rule does not fire.
  tune <- tune_stub()
  st <- ups(pm_stub(gd_good = FALSE),
                                  c(1, rep(1e-12, 40)), rep(10, 4L), tune,
                                  n_pars = 5L)
  expect_lt(st$weight_ess, 1.001)
  st <- ups(st, rep(1, 41), rep(10, 4L), tune,
                                  n_pars = 5L)
  expect_equal(st$weight_ess, 41)
  st$gd_good <- TRUE
  st <- ups(st, rep(1, 21), rep(5, 4L), tune,
                                  n_pars = 5L)
  expect_equal(st$weight_ess, 21)
})

test_that("recording it does not change the adaptation", {
  # The particle count must follow the same rule it always did.
  w <- c(0.5, runif(60))
  tune <- tune_stub()
  got <- ups(pm_stub(), w, rep(15, 4L), tune,
                                   n_pars = 5L)
  ess <- sum(w)^2 / sum(w^2)
  want <- max(25, min(tune$max_particles,
                      round(100 * (tune$target_ESS / ess)^tune$ESS_scale)))
  expect_identical(got$n_particles, want)
  # And the floor and ceiling still bind.
  tiny <- ups(pm_stub(n_particles = 26),
                                    c(1, rep(1e-12, 50)), rep(12, 4L), tune,
                                    n_pars = 5L)
  expect_gte(tiny$n_particles, 25)
  expect_lte(tiny$n_particles, tune$max_particles)
})

test_that("weight ESS and chain ESS are different quantities", {
  # The point of recording it. A cloud concentrated on one draw has a weight ESS
  # near 1 however well the chain is mixing; a flat cloud has a weight ESS near
  # its own size however stuck the chain is. Neither statement is about
  # autocorrelation, which is what ess_bulk measures.
  flat <- rep(1, 200)
  spike <- c(1, rep(1e-12, 199))
  expect_equal(sum(flat)^2 / sum(flat^2), 200)
  expect_lt(sum(spike)^2 / sum(spike^2), 1.001)
  # A heavily autocorrelated chain has a chain ESS far below its length, and the
  # weight ESS above says nothing about it either way.
  set.seed(31)
  slow <- matrix(c(cumsum(rnorm(200)), cumsum(rnorm(200))), ncol = 2L)
  expect_lt(EMC2:::ess_bulk(slow), 200)
  fast <- matrix(rnorm(400), ncol = 2L)
  expect_gt(EMC2:::ess_bulk(fast), EMC2:::ess_bulk(slow))
  # And a chain that never moved has no chain ESS at all, while its weight ESS
  # would still be a perfectly finite number -- which is the confusion.
  expect_true(is.na(EMC2:::ess_bulk(matrix(rep(1.5, 400), ncol = 2L))))
})

test_that("the schema carries the weight ESS beside the chain's", {
  sch <- EMC2:::.emc_profile_schema
  expect_true("particle_weight_ess" %in% names(sch))
  expect_match(sch[["particle_weight_ess"]]$desc, "proposal cloud")
  row <- EMC2:::.emc_profile_row(iteration = 1L, particle_weight_ess = 42)
  expect_identical(row$particle_weight_ess, 42)
  # Unmeasured stays NA rather than 0: a cloud of size zero is not a claim
  # anyone should be able to read out of an ordinary run.
  expect_true(is.na(EMC2:::.emc_profile_row(iteration = 1L)$particle_weight_ess))
})

test_that("the report names it for what it is", {
  rows <- lapply(1:3, function(i) {
    EMC2:::.emc_profile_row(iteration = i, stage = "sample", total = 0.1,
                            subject_cpu_max = 0.01, worker_cpu_max = 0.02,
                            subject_cpu_sum = 0.03, particle_weight_ess = 27)
  })
  out <- capture.output(
    EMC2:::.emc_profile_report(EMC2:::.emc_profile_bind(rows), drop_first = FALSE))
  txt <- paste(out, collapse = " ")
  expect_match(txt, "weight_ess")
  # Named so it cannot be read as the chain's ESS, which is the confusion the
  # whole item is about.
  expect_match(txt, "NOT chain ESS")
})

test_that("a tune without search_width still adapts epsilon at every stage", {
  # `run_stages()` defaults `search_width` to NULL. It used to reach
  # `set_p_accept()` as `1/NULL`, give an empty target, and empty every epsilon
  # at the first adaptation -- after which half of every proposal cloud was NaN
  # and every update repeated the previous state.
  n_pars <- 5L
  for (stage in c("preburn", "burn", "adapt", "sample")) {
    unset <- EMC2:::check_tune_settings(list(), n_pars, stage, 100)
    given <- EMC2:::check_tune_settings(list(search_width = 1), n_pars, stage, 100)
    expect_identical(unset$p_accept, given$p_accept, info = stage)
    pm <- EMC2:::check_sampling_settings(list(list()), stage, n_pars, 100)[[1L]]
    n_mix <- length(pm$mix)
    pm$iter <- unset$n0 + 1
    w <- c(0.5, runif(n_mix * 10))
    got <- ups(pm, w, rep(10, n_mix), unset,
                                     n_pars = n_pars)
    # One epsilon per non-group proposal, all finite -- the shape
    # `new_particle()` indexes when it scales proposals 2..n_mix.
    expect_length(got$epsilon, n_mix - 1L)
    expect_true(length(got$epsilon) > 0L && all(is.finite(got$epsilon)),
                info = stage)
  }
})

test_that("the local step settles at the random-walk scale without overshoot", {
  # Share of particles that beat the reference for a Gaussian conditional with a
  # matched covariance; the controller should settle at 2.38 / sqrt(d).
  settle <- function(n_pars, eps0) {
    tune <- EMC2:::check_tune_settings(list(), n_pars, "adapt", 100)
    tune$n0 <- 0L
    pm <- pm_stub(); pm$iter <- 0L; pm$epsilon <- rep(eps0, 3L)
    path <- numeric(600)
    for (i in seq_along(path)) {
      share <- pnorm(-pm$epsilon[1L] * sqrt(n_pars / 2))
      pm <- EMC2:::update_pm_settings(
        pm, c(0, round(1e4 * share), 0, 0), c(0, 1e4, 1e4, 1e4), 1, tune,
        n_pars = n_pars, layout = c("chains", "", ""))
      path[i] <- pm$epsilon[1L] * sqrt(n_pars)
    }
    path
  }
  for (n_pars in c(14L, 71L)) {
    # The first call records the layout and starts at the random-walk scale.
    from_small <- settle(n_pars, 0.05)
    expect_equal(from_small[1L], 2.38)
    expect_lt(abs(from_small[600L] - 2.38), 0.05)
    expect_lt(max(from_small), 2.38 * 1.05)
  }
})

test_that("a local step restarts when its covariance changes, and has a floor", {
  tune <- EMC2:::check_tune_settings(list(), 24L, "burn", 100)
  tune$n0 <- 0L
  pm <- pm_stub(n_mix = 3L); pm$iter <- 0L
  never <- function(pm, layout, n = 1L) {
    for (i in seq_len(n)) pm <- EMC2:::update_pm_settings(
      pm, rep(0, 3L), rep(10, 3L), 1, tune, n_pars = 24L, layout = layout)
    pm
  }
  # The first call starts both random walks at the random-walk scale.
  pm <- never(pm, c("group", "chains"))
  expect_equal(pm$epsilon, rep(2.38 / sqrt(24), 2L))
  # No particle ever beats the reference: the step shrinks, but not below
  # a quarter of that scale.
  shrunk <- never(pm, c("group", "chains"), 50L)
  expect_true(all(shrunk$epsilon < pm$epsilon))
  pm$epsilon <- rep(1e-3, 2L)
  pm <- never(pm, c("group", "chains"))
  expect_equal(pm$epsilon, rep(2.38 / sqrt(24) / 4, 2L))
  # The efficient-proposal mixture reuses slot 1 for the chains random walk.
  pm$mix <- rep(0.25, 4L); pm$proposal_counts <- pm$acc_counts <- rep(0, 4L)
  pm <- EMC2:::update_pm_settings(pm, rep(0, 4L), rep(10, 4L), 1, tune,
                                  n_pars = 24L, layout = c("chains", "", ""))
  expect_equal(pm$epsilon[1L], 2.38 / sqrt(24))
  # Independence components keep their (unused) multiplier.
  expect_equal(pm$epsilon[2L], 2.38 / sqrt(24) / 4)
})

test_that("mixture weights follow squared jump per particle, without targets", {
  tune <- EMC2:::check_tune_settings(list(), 10L, "adapt", 100)
  tune$n0 <- 0L
  pm <- pm_stub(gd_good = FALSE); pm$iter <- 0L
  pm$proposal_counts <- pm$acc_counts <- rep(0, 4L)
  pm$jump_counts <- numeric(4L)
  # Same particle counts; component 3 delivers three times the jump of 4, and
  # the others none. The weights converge to those shares (with the floor).
  for (i in 1:2000) pm <- EMC2:::update_pm_settings(
    pm, rep(0, 4L), rep(10, 4L), 1, tune, n_pars = 10L,
    layout = c("chains", "", ""), jump = c(0, 0, 3, 1))
  expect_equal(pm$mix, c(0.02, 0.02, 0.72, 0.24) / 1, tolerance = 0.02)
  # Without jumps nothing moves.
  pm2 <- pm_stub(gd_good = FALSE); pm2$jump_counts <- numeric(4L)
  pm2 <- EMC2:::update_pm_settings(pm2, rep(0, 4L), rep(10, 4L), 1, tune,
                                   n_pars = 10L, jump = numeric(4L))
  expect_equal(pm2$mix, pm_stub()$mix)
})

test_that("the efficient proposal mean follows the current group mean", {
  m <- c(1, 2)
  attr(m, "B") <- matrix(c(0.5, 0, 0, 0.25), 2)
  attr(m, "given") <- c(0, 0)
  expect_equal(EMC2:::.emc_eff_mean(m, c(2, 4)), c(2, 3))
  expect_equal(EMC2:::.emc_eff_mean(m, c(0, 0)), c(1, 2))
  expect_equal(EMC2:::.emc_eff_mean(c(1, 2), c(2, 4)), c(1, 2))
})

test_that("window covariances keep the posterior correlations", {
  # The chain proposal is an unscaled independence proposal: its covariance
  # must be the window's own, correlations included. Independent chains use
  # their own window; otherwise the window pools every chain.
  emc <- samples_LNR
  for (j in seq_along(emc)) emc[[j]]$samples$stage[] <- "adapt"
  own <- lapply(emc, function(ch) ch$samples$alpha[, 1, 1:50])
  for (independent in c(TRUE, FALSE)) {
    withr::local_options(list(emc2.independent_chains = independent))
    out <- EMC2:::create_chain_proposals(emc, samples_idx = 1:50, stage = "adapt",
                                          do_block = FALSE)
    a <- if (independent) own[[1]] else do.call(cbind, own)
    expect_equal(out[[1]]$chains_var[[1]], stats::cov(t(a)), ignore_attr = TRUE)
    expect_equal(out[[1]]$ensemble_q$var[[1]], stats::cov(t(a)), ignore_attr = TRUE)
  }
})
