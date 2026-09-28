# The conditional importance-sampling particle update must leave the subject
# posterior invariant. Two past defects broke this while every chain still
# looked converged: (1) numbers_from_proportion() forced at least one particle
# per mixture component, so rare components were over-sampled relative to the
# mixture density used in the weights; (2) proposal components centred on the
# current state made the reference particle's proposal density depend on
# itself. Both inflated subject-level spread, and so the group variance, by an
# amount that depended on each fit's adapted proposal settings.

local_rng_guard()

test_that("particle allocation is plain multinomial (no forced minimum)", {
  set.seed(1)
  counts <- replicate(200, EMC2:::numbers_from_proportion(c(.98, .01, .01), 5))
  expect_true(any(counts[2, ] == 0))
  expect_true(all(colSums(counts) == 5))
})

test_that("new_particle leaves a conjugate Gaussian posterior invariant", {
  skip_model_validation()
  ybar <- c(1.5, -0.5); Li <- diag(1 / c(0.3, 0.5))
  mu <- c(0, 0); Sig <- matrix(c(1, .3, .3, 1), 2)
  post_var <- solve(solve(Sig) + Li)
  post_mu <- as.vector(post_var %*% (solve(Sig) %*% mu + Li %*% ybar))
  local_mocked_bindings(
    calc_ll_pooled = function(proposals, ...) {
      d <- sweep(proposals, 2, ybar)
      -0.5 * rowSums((d %*% Li) * d)
    },
    .package = "EMC2")
  bm_se <- function(x, nb = 40) {
    b <- split(x, cut(seq_along(x), nb)); sd(vapply(b, mean, 0)) / sqrt(nb)
  }
  run <- function(stage, mix, eps, iters = 20000) {
    set.seed(3)
    pm <- list(list(epsilon = eps, mix = mix, n_particles = 5, iter = 0,
                    proposal_counts = rep(0, length(mix)),
                    acc_counts = rep(0, length(mix)), gd_good = FALSE))
    tune <- list(components = c(1, 1), shared_ll_idx = c(1, 1), n0 = Inf)
    a <- c(a = 0, b = 0)
    ll <- -0.5 * sum((a - ybar)^2 * diag(Li))
    out <- matrix(NA_real_, iters, 2)
    for (t in seq_len(iters)) {
      r <- EMC2:::new_particle(
        1, NULL, pm, eff_mu = c(.8, -.2), eff_var = diag(2) * .3,
        chains_mu = c(.8, -.2), chains_var = diag(2) * .2, prev_ll = ll,
        parameters = NULL, stage = stage, type = "standard", tune = tune,
        current_alpha = a, population_mu = mu, population_var = Sig)
      a <- r$proposal; ll <- r$ll; out[t, ] <- a
    }
    out[-(1:1000), 1]
  }
  # Mixes with a state-centred component and a rare component: the old
  # kernel missed the mean here by z = -18 (sample) and z = -22 (adapt).
  for (cfg in list(list("sample", c(.2, .4, .2, .2), c(1, 1, 1)),
                   list("adapt", c(.3, .4, .3), c(1, 1)),
                   list("sample", c(.3, .02, .3, .38), c(1, 1, 1)))) {
    x <- run(cfg[[1]], cfg[[2]], cfg[[3]])
    expect_lt(abs(mean(x) - post_mu[1]) / bm_se(x), 4.5)
    expect_lt(abs(var(x) - post_var[1, 1]) / bm_se((x - mean(x))^2), 4.5)
  }
})
