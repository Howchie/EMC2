# Regression tests for design(marginalise = "t0"): the shared subject-level t0
# is integrated out of the MH proposal (fixed-eta Stage 1) and reconstructed
# into the stored alpha. These guard the two defects found in review:
#   (1) a broken marginalised index that froze the whole chain, and
#   (2) accidental integration leaking into downstream (predict/make_data/IC).

matchfun <- function(d) d$S == d$lR

# Two subjects keeps the sampler exercise cheap.
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2], ]
dat$subjects <- droplevels(dat$subjects)

pmean <- c(v = 1, B = log(.5), B_Eneutral = log(1.5), B_Eaccuracy = log(2),
           B_lRright = 0, A = log(0.25), t0 = log(.2))
psd <- c(v = 1, B = 0.3, B_Eneutral = 0.3, B_Eaccuracy = 0.3,
         B_lRright = 0.3, A = 0.4, t0 = .5)

make_marginal_emc <- function() {
  des <- design(data = dat, model = LBA,
                formula = list(v ~ 1, sv ~ 1, B ~ E + lR, A ~ 1, t0 ~ 1),
                constants = c(sv = log(1)), marginalise = "t0",
                report_p_vector = FALSE)
  pr <- prior(des, type = "standard", pmean = pmean, psd = psd)
  list(des = des,
       emc = make_emc(dat, des, type = "standard", prior = pr,
                      compress = TRUE, n_chains = 1))
}

test_that("marginalise is carried on the design but never tagged on the dadm", {
  mm <- make_marginal_emc()
  expect_identical(attr(mm$des, "marginalise"), "t0")
  # The sampler resolves the fixed eta from the existing t0 prior.
  expect_identical(mm$emc[[1]]$marginalise$param, "t0")
  expect_true(is.finite(mm$emc[[1]]$marginalise$mu))
  expect_true(mm$emc[[1]]$marginalise$sigma > 0)
  # Exactly t0 is flagged as marginalised in the sampled vector.
  idx <- mm$emc[[1]]$marginalised_idx
  expect_identical(mm$emc[[1]]$par_names[idx], "t0")
  # Downstream consumers must see an ordinary (untagged) dadm.
  expect_null(attr(mm$emc[[1]]$data[[1]], "marginalise"))
})

test_that("a plain likelihood on the dadm does NOT integrate t0", {
  mm <- make_marginal_emc()
  s1 <- mm$emc[[1]]
  d1 <- s1$data[[1]]
  P <- matrix(pmean[s1$par_names], nrow = 1L,
              dimnames = list(NULL, s1$par_names))
  ll_plain <- EMC2:::calc_ll_manager(P, dadm = d1, s1$model)
  grid <- EMC2:::compute_marginal_grid(P, d1, s1$model, s1$marginalise)
  ll_marg <- EMC2:::marginal_ll_from_grid(grid$log_terms)
  # The two must differ substantially: if they matched, the plain path would be
  # silently integrating (the fallback bug this feature had).
  expect_gt(abs(ll_plain - ll_marg), 1)
})

test_that("each proposal's quadrature is independent of the rest of the batch", {
  # The rule used to be centred on mass pooled over the whole particle batch, so
  # a scattered batch (preburn/burn) gave per-particle quadrature errors of many
  # nats that did NOT cancel in the MH weights.  Rules are now per particle: a
  # particle's marginal ll must not depend on its companions.
  mm <- make_marginal_emc()
  s1 <- mm$emc[[1]]
  d1 <- s1$data[[1]]
  set.seed(11)
  np <- 8L
  P <- matrix(rep(pmean[s1$par_names], each = np), nrow = np,
              dimnames = list(NULL, s1$par_names))
  P[-1, ] <- P[-1, ] + matrix(rnorm((np - 1L) * ncol(P), 0, 0.15), np - 1L)
  batch <- EMC2:::marginal_ll_from_grid(
    EMC2:::compute_marginal_grid(P, d1, s1$model, s1$marginalise)$log_terms)
  solo <- vapply(seq_len(np), function(i) {
    EMC2:::marginal_ll_from_grid(
      EMC2:::compute_marginal_grid(P[i, , drop = FALSE], d1, s1$model,
                                   s1$marginalise)$log_terms)[1]
  }, numeric(1))
  expect_true(all(is.finite(batch)))
  expect_equal(batch, solo, tolerance = 1e-6)
  # Every particle gets its own node row, and the rule is centred where that
  # particle's own conditional posterior actually sits.
  grid <- EMC2:::compute_marginal_grid(P, d1, s1$model, s1$marginalise)
  expect_identical(dim(grid$nodes), dim(grid$log_terms))
  expect_identical(nrow(grid$nodes), np)
  modes <- grid$nodes[cbind(seq_len(np), max.col(grid$log_terms))]
  expect_gt(diff(range(modes)), 0)
})

test_that("a warm start reproduces the cold rule, and a stale one is discarded", {
  # The accepted particle's Laplace fit is carried to the next iteration to
  # replace the pilot scan. It may only ever save work: the answer must not
  # depend on whether a hint was supplied, and a hint that no longer brackets
  # the batch must be thrown away rather than quietly narrowing the rule.
  mm <- make_marginal_emc()
  s1 <- mm$emc[[1]]
  d1 <- s1$data[[1]]
  set.seed(13)
  np <- 6L
  P <- matrix(rep(pmean[s1$par_names], each = np), nrow = np,
              dimnames = list(NULL, s1$par_names))
  P[-1, ] <- P[-1, ] + matrix(rnorm((np - 1L) * ncol(P), 0, 1e-3), np - 1L)
  cold_grid <- EMC2:::compute_marginal_grid(P, d1, s1$model, s1$marginalise)
  cold <- EMC2:::marginal_ll_from_grid(cold_grid$log_terms)
  warm <- EMC2:::marginal_warm_state(cold_grid, 1L)
  expect_true(all(is.finite(warm)))

  warm_grid <- EMC2:::compute_marginal_grid(P, d1, s1$model, s1$marginalise,
                                            warm = warm)
  expect_identical(warm_grid$warm_used, 1L)
  expect_equal(EMC2:::marginal_ll_from_grid(warm_grid$log_terms), cold,
               tolerance = 1e-4)

  # A hint pointing somewhere else entirely: the probe must detect that and fall
  # back to the pilot, giving exactly the cold answer.
  stale <- c(mode = warm[["mode"]] + 1.5, sd = warm[["sd"]] * 50)
  stale_grid <- EMC2:::compute_marginal_grid(P, d1, s1$model, s1$marginalise,
                                             warm = stale)
  expect_identical(stale_grid$warm_used, 0L)
  expect_identical(EMC2:::marginal_ll_from_grid(stale_grid$log_terms), cold)

  # Backoff: a rejected hint is retried with geometrically growing delay, so a
  # batch whose spread the hint can never cover pays at most one probe in 16.
  st <- NULL
  waits <- integer(0)
  for (k in 1:6) {
    attempted <- is.null(st) || st$backoff <= 0L
    st <- EMC2:::marginal_warm_backoff(st, attempted = attempted, used = FALSE)
    if (attempted) waits <- c(waits, st$backoff)
  }
  expect_true(all(diff(waits) > 0))
  expect_equal(EMC2:::marginal_warm_backoff(st, attempted = TRUE, used = TRUE)$backoff, 0L)
})

test_that("the marginalise chain mixes (no freeze) and pins t0's group level", {
  skip_on_cran()
  skip_on_os("windows")
  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
  mm <- make_marginal_emc()
  # Preburn-only keeps this fast and deterministic (no gd auto-extension).
  emc <- EMC2:::run_emc(mm$emc, stage = "preburn",
                        stop_criteria = list(iter = 12),
                        cores_for_chains = 1, particle_factor = 15,
                        verbose = FALSE)
  samples <- emc[[1]]$samples
  alpha <- samples$alpha
  idx <- samples$idx

  # (1) The freeze bug produced EXACTLY zero movement in every parameter for
  # every subject. Every non-marginalised coordinate must actually move.
  non_t0 <- rownames(alpha) != "t0"
  movement <- apply(alpha[non_t0, , , drop = FALSE], 1L,
                    function(x) diff(range(x)))
  expect_true(all(movement > 0))

  # (2) t0 is reconstructed into alpha (kept for predict/make_data) and, being a
  # posterior draw per iteration, varies and stays finite.
  expect_true(all(is.finite(alpha["t0", , ])))
  expect_gt(diff(range(alpha["t0", , ])), 0)

  # (3) Stage-1 holds t0's group level at the fixed prior eta: gibbs pins the
  # group mean/variance deterministically each iteration.
  expect_equal(unname(samples$theta_mu["t0", idx]),
               emc[[1]]$marginalise$mu, tolerance = 1e-8)
  expect_equal(unname(samples$theta_var["t0", "t0", idx]),
               emc[[1]]$marginalise$sigma^2, tolerance = 1e-8)
})

test_that("marginalise supports single-subject alpha-only samplers", {
  skip_on_cran()
  skip_on_os("windows")

  dat_single <- data.frame(
    subjects = factor(rep("s1", 8)),
    S = factor(c("left", "right", "left", "right", "left", "right", "left", "right")),
    R = factor(c("left", "right", "nogo", "left", "right", "nogo", "left", "nogo"),
               levels = c("left", "right", "nogo")),
    rt = c(.45, .50, Inf, .40, .48, Inf, .42, Inf)
  )
  matchfun_single <- function(d) as.character(d$S) == as.character(d$lR)
  des_single <- design(
    data = dat_single, model = LBA, matchfun = matchfun_single,
    formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
    constants = c(sv = log(1)), marginalise = "t0",
    report_p_vector = FALSE
  )

  set.seed(123)
  emc <- make_emc(dat_single, des_single, type = "single",
                  compress = FALSE, n_chains = 1)
  sampler <- emc[[1]]
  expect_identical(sampler$type, "single")
  expect_null(sampler$samples$theta_mu)
  expect_identical(sampler$marginalise$param, "t0")

  emc <- EMC2:::run_emc(
    emc, stage = "preburn", stop_criteria = list(iter = 3),
    cores_for_chains = 1, particle_factor = 5, verbose = FALSE
  )
  sampler <- emc[[1]]
  samples <- sampler$samples
  idx <- samples$idx
  expect_true(all(is.finite(samples$alpha["t0", 1, seq_len(idx)])))
  expect_gt(diff(range(samples$alpha["t0", 1, seq_len(idx)])), 0)

  # The quadrature path is active for this alpha-only sampler.  The integrated
  # and ordinary likelihoods should differ; the stored subject likelihood is
  # computed at the accepted particle's proposal grid, whose adaptive nodes
  # can differ from a fresh one-row grid here.
  proposal <- matrix(samples$alpha[, 1, idx], nrow = 1L,
                     dimnames = list(NULL, sampler$par_names))
  grid <- EMC2:::compute_marginal_grid(
    proposal, sampler$data[[1]], sampler$model, sampler$marginalise
  )
  ll_marg <- EMC2:::marginal_ll_from_grid(grid$log_terms)[1]
  ll_plain <- EMC2:::calc_ll_manager(proposal, sampler$data[[1]], sampler$model)
  expect_true(is.finite(samples$subj_ll[1, idx]))
  expect_true(is.finite(ll_marg))
  expect_gt(abs(ll_plain - ll_marg), 0.1)
})
