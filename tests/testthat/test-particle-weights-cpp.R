local_rng_guard()

legacy_particle_weights <- function(proposals, log_likelihood,
                                    group_log_density, mix, components,
                                    means, rootis, log_constants,
                                    auxiliary_components, reference_ll) {
  terms <- matrix(0, nrow(proposals), length(components))
  for (j in seq_along(components)) {
    k <- components[[j]]
    density <- if (k == 1L) group_log_density else
      EMC2:::fast_dmvnorm_rooti(
        proposals, means[[k]], rootis[[k]], log_constants[[k]])
    terms[, j] <- log(mix[[k]]) + density
  }
  largest <- terms[, 1L]
  if (ncol(terms) > 1L)
    for (j in 2:ncol(terms)) largest <- pmax(largest, terms[, j])
  log_mix <- largest + log(rowSums(exp(terms - largest)))
  bad <- is.infinite(log_mix) | is.na(log_mix)
  if (any(bad)) {
    finite <- log_mix[!bad]
    log_mix[bad] <- if (length(finite)) min(finite) else -1e10
  }
  auxiliary <- numeric(nrow(proposals))
  for (k in auxiliary_components)
    auxiliary <- auxiliary + EMC2:::fast_dmvnorm_rooti(
      proposals, means[[k]], rootis[[k]], log_constants[[k]])
  total <- log_likelihood + reference_ll - log_likelihood[[1L]]
  log_weight <- total + group_log_density + auxiliary - log_mix
  exp(log_weight - max(log_weight))
}

test_that("C++ particle weights match the R reference to floating precision", {
  proposals <- matrix(c(
    -1.1, .3, .8, 1.2, -.4, .2, .7, -.9, .1, 1.4, .5, -.8,
    -.3, .9, -1.2, .6, .2, .4, -1.4, .1, .5, -.7, .3, 1.1,
    .2, -.5, 1.3, -.8, .6, -.2), ncol = 3L)
  covs <- list(diag(c(.7, 1.1, 1.4)),
               matrix(c(1.2, .1, 0, .1, .8, .05, 0, .05, 1.5), 3L),
               diag(c(1.4, .9, .6)),
               matrix(c(.9, .1, .05, .1, 1.3, 0, .05, 0, .8), 3L))
  factors <- lapply(covs, EMC2:::.chol_factor, p = 3L)
  rootis <- lapply(factors, `[[`, "rooti")
  log_constants <- vapply(factors, `[[`, numeric(1), "sum_log_diag") -
    .5 * 3L * log(2 * pi)
  means <- list(c(0, 0, 0), c(.1, -.2, .3), c(-.2, .1, .1), c(.3, .2, -.1))
  mix <- c(.25, .3, .15, .3)
  ll <- c(-4, -3, -5, -2, -4, -6, -2, -3, -5, -4)
  group_ll <- EMC2:::fast_dmvnorm_rooti(
    proposals, means[[1L]], rootis[[1L]], log_constants[[1L]])
  cases <- list(
    list(components = c(1L, 3L, 4L), auxiliary = integer(0)),
    list(components = 2L, auxiliary = 2L),
    list(components = c(1L, 2L, 3L, 4L), auxiliary = 2L))

  for (case in cases) {
    expected <- legacy_particle_weights(
      proposals, ll, group_ll, mix, case$components, means, rootis,
      log_constants, case$auxiliary, reference_ll = -.7)
    actual <- EMC2:::emc2_particle_log_weights(
      proposals, ll, group_ll, mix, case$components, means, rootis,
      log_constants, case$auxiliary, reference_log_likelihood = -.7)
    expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 2e-15)
  }
})

test_that("the C++ weight path preserves an EMC2 particle transition", {
  local_mocked_bindings(
    calc_ll_pooled = function(proposals, ...) -0.5 * rowSums(proposals^2),
    .package = "EMC2")
  p <- 5L
  alpha <- stats::setNames(rep(0, p), paste0("p", seq_len(p)))
  settings_for <- function(config) list(list(
    epsilon = config$epsilon, mix = config$mix,
    n_particles = config$n_particles, iter = 0L,
    proposal_counts = integer(length(config$mix)),
    acc_counts = integer(length(config$mix)), gd_good = FALSE))
  tune <- list(components = rep(1L, p), shared_ll_idx = rep(1L, p), n0 = Inf)
  args <- list(
    s = 1L, data = NULL, pm_settings = NULL,
    eff_mu = alpha + .2, eff_var = diag(p) * .7,
    chains_mu = alpha - .1, chains_var = diag(p) * .8,
    prev_ll = 0, parameters = NULL, model = NULL,
    stage = "sample", type = "single", tune = tune,
    current_alpha = alpha, population_mu = alpha,
    population_var = diag(p))
  configs <- list(
    list(stage = "preburn", mix = c(.4, .6), epsilon = c(.7, .9),
         n_particles = 20L),
    list(stage = "burn", mix = c(.15, .35, .5), epsilon = c(.7, .9),
         n_particles = 30L),
    list(stage = "adapt", mix = c(.05, .3, .3, .35), epsilon = c(.7, .8, .9),
         n_particles = 30L),
    list(stage = "sample", mix = c(.05, .3, .3, .35), epsilon = c(.7, .8, .9),
         n_particles = 30L))

  old_option <- getOption("emc2.particle_weights_cpp")
  on.exit(options(emc2.particle_weights_cpp = old_option), add = TRUE)
  for (config in configs) {
    args$stage <- config$stage
    args$pm_settings <- settings_for(config)
    for (seed in c(481L, 1234L, 9001L)) {
      options(emc2.particle_weights_cpp = FALSE)
      set.seed(seed)
      reference <- do.call(EMC2:::new_particle, args)
      reference_rng <- .Random.seed

      options(emc2.particle_weights_cpp = TRUE)
      set.seed(seed)
      candidate <- do.call(EMC2:::new_particle, args)
      candidate_rng <- .Random.seed

      expect_identical(candidate, reference,
                       info = paste(config$stage, "seed", seed))
      expect_identical(candidate_rng, reference_rng,
                       info = paste(config$stage, "RNG", seed))
    }
  }
})
