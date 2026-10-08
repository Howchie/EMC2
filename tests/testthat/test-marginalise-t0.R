local_rng_guard()  # see helper-rng.R: keep this file's RNG changes inside it
# Regression tests for design(marginalise = "t0"): the shared subject-level t0
# is integrated out of the particle step against its conditional group prior
# and redrawn into the stored alpha from its conditional posterior. These guard the two defects found in review:
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
  # Evaluations outside the sampler integrate against the t0 prior.
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

test_that("the starting point only saves work, and screening only skips negligible particles", {
  # Each rule starts at a predicted peak (the subject's last batch, or its
  # posterior moments). The start may only change how many kernel passes a
  # particle takes, not where its rule ends up beyond quadrature error.
  mm <- make_marginal_emc()
  s1 <- mm$emc[[1]]
  d1 <- s1$data[[1]]
  pn <- s1$par_names
  set.seed(13)
  np <- 6L
  P <- matrix(rep(pmean[pn], each = np), nrow = np, dimnames = list(NULL, pn))
  P[-1, ] <- P[-1, ] + matrix(rnorm((np - 1L) * ncol(P), 0, 0.05), np - 1L)
  cold <- EMC2:::compute_marginal_grid(P, d1, s1$model, s1$marginalise)
  expect_true(all(is.finite(cold$ll)))

  warm_spec <- s1$marginalise
  warm_spec$centre <- cold$mode
  warm_spec$probe_sd <- cold$sd
  warm <- EMC2:::compute_marginal_grid(P, d1, s1$model, warm_spec)
  expect_lte(sum(warm$rounds), sum(cold$rounds))
  expect_lt(max(abs(warm$ll - cold$ll)), 0.01)

  stale_spec <- warm_spec
  stale_spec$centre <- cold$mode + 10 * cold$sd
  stale <- EMC2:::compute_marginal_grid(P, d1, s1$model, stale_spec)
  expect_gt(sum(stale$rounds), sum(warm$rounds))
  expect_lt(max(abs(stale$ll - cold$ll)), 0.01)

  # The regression state from one batch predicts the next batch's peaks.
  V <- diag(psd[pn]^2)
  dimnames(V) <- list(pn, pn)
  cond <- EMC2:::.marginal_conditional(pmean[pn], V,
                                       EMC2:::.marginal_par_idx(pn, s1$marginalise))
  big <- matrix(rep(pmean[pn], each = 40L), nrow = 40L, dimnames = list(NULL, pn))
  big <- big + matrix(rnorm(40L * ncol(big), 0, 0.03), 40L)
  gb <- EMC2:::compute_marginal_grid(big, d1, s1$model, s1$marginalise)
  state <- EMC2:::.marginal_state(gb, big, cond)
  expect_false(is.null(state))
  spec <- EMC2:::.marginal_batch_spec(s1$marginalise, cond, P, state = state)
  expect_lt(max(abs(spec$centre - cold$mode) / cold$sd), 1)

  # A particle whose weight is far below the best keeps its first rule;
  # the ones that can be selected are untouched.
  scr_spec <- stale_spec
  scr_spec$log_offset <- c(0, -1e4, rep(0, np - 2L))
  scr <- EMC2:::compute_marginal_grid(P, d1, s1$model, scr_spec)
  expect_gte(scr$screened, 1L)
  expect_equal(scr$rounds[2], 1L)
  live <- stale$ll >= max(stale$ll) - 10
  expect_identical(scr$ll[live], stale$ll[live])
})

test_that("the marginalise chain mixes (no freeze) and estimates t0's group level", {
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

  # (3) t0's group mean, variance and covariances are sampled like any other
  # parameter's.
  it <- 2:idx
  expect_gt(diff(range(samples$theta_mu["t0", it])), 0)
  expect_gt(diff(range(samples$theta_var["t0", "t0", it])), 0)
  expect_true(any(samples$theta_var["t0", "v", it] != 0))
  expect_true(all(is.finite(samples$subj_ll[, it])))
})

test_that("the t0 draw follows its conditional posterior", {
  # Given the other parameters, t0 has density proportional to
  # L(y | x, t0) N(t0 | conditional group prior). Compare the draws with that
  # density on a fine grid.
  mm <- make_marginal_emc()
  s1 <- mm$emc[[1]]
  d1 <- s1$data[[1]]
  pn <- s1$par_names
  V <- diag(psd[pn]^2)
  dimnames(V) <- list(pn, pn)
  V["t0", "v"] <- V["v", "t0"] <- 0.3 * psd[["t0"]] * psd[["v"]]
  m_idx <- EMC2:::.marginal_par_idx(pn, s1$marginalise)
  cond <- EMC2:::.marginal_conditional(pmean[pn], V, m_idx)
  x <- pmean[pn]
  x[["v"]] <- x[["v"]] + 0.2
  spec <- EMC2:::.marginal_batch_spec(s1$marginalise, cond,
                                      matrix(x, 1, dimnames = list(NULL, pn)))
  grid <- EMC2:::compute_marginal_grid(matrix(x, 1, dimnames = list(NULL, pn)),
                                       d1, s1$model, spec)
  set.seed(7)
  draws <- vapply(1:1500, function(i) {
    EMC2:::.marginal_draw(x, cond, as.numeric(grid$nodes[1, ]),
                          as.numeric(grid$log_terms[1, ]), grid$mode[1],
                          grid$sd[1], d1, s1$model)$x[["t0"]]
  }, numeric(1))
  g <- seq(grid$mode[1] - 8 * grid$sd[1], grid$mode[1] + 8 * grid$sd[1],
           length.out = 801)
  P <- matrix(x, length(g), length(x), byrow = TRUE, dimnames = list(NULL, pn))
  P[, "t0"] <- g
  lp <- EMC2:::calc_ll_manager(P, d1, s1$model) +
    dnorm(g, spec$mu, spec$sigma, log = TRUE)
  w <- exp(lp - max(lp)); w <- w / sum(w)
  m <- sum(w * g); sdv <- sqrt(sum(w * (g - m)^2))
  expect_lt(abs(mean(draws) - m), 4 * sdv / sqrt(length(draws)))
  expect_equal(sd(draws), sdv, tolerance = 0.1)
  # Continuous, not confined to the quadrature nodes.
  expect_gt(length(unique(draws)), 100)
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

test_that("flexible marginalise n_nodes options are parsed and propagated correctly", {
  # 1. Vector form: c("t0", n_nodes = 10)
  des1 <- design(data = dat, model = LBA,
                 formula = list(v ~ 1, sv ~ 1, B ~ E + lR, A ~ 1, t0 ~ 1),
                 constants = c(sv = log(1)), marginalise = c("t0", n_nodes = 10),
                 report_p_vector = FALSE)
  m_attr1 <- attr(des1, "marginalise")
  expect_equal(m_attr1$param, "t0")
  expect_equal(m_attr1$n_nodes, 10L)

  pr1 <- prior(des1, type = "standard", pmean = pmean, psd = psd)
  emc1 <- make_emc(dat, des1, type = "standard", prior = pr1, compress = TRUE, n_chains = 1)
  expect_equal(emc1[[1]]$marginalise$n_nodes, 10L)
  expect_equal(emc1[[1]]$marginalise$param, "t0")

  # 2. Named vector form: c(param = "t0", n_nodes = 8)
  des2 <- design(data = dat, model = LBA,
                 formula = list(v ~ 1, sv ~ 1, B ~ E + lR, A ~ 1, t0 ~ 1),
                 constants = c(sv = log(1)), marginalise = c(param = "t0", n_nodes = 8),
                 report_p_vector = FALSE)
  m_attr2 <- attr(des2, "marginalise")
  expect_equal(m_attr2$param, "t0")
  expect_equal(m_attr2$n_nodes, 8L)

  # 3. List form: list("t0", n_nodes = 10)
  des3 <- design(data = dat, model = LBA,
                 formula = list(v ~ 1, sv ~ 1, B ~ E + lR, A ~ 1, t0 ~ 1),
                 constants = c(sv = log(1)), marginalise = list("t0", n_nodes = 10),
                 report_p_vector = FALSE)
  m_attr3 <- attr(des3, "marginalise")
  expect_equal(m_attr3$param, "t0")
  expect_equal(m_attr3$n_nodes, 10L)

  # 4. Standard string default: marginalise = "t0" defaults n_nodes to 5
  des4 <- design(data = dat, model = LBA,
                 formula = list(v ~ 1, sv ~ 1, B ~ E + lR, A ~ 1, t0 ~ 1),
                 constants = c(sv = log(1)), marginalise = "t0",
                 report_p_vector = FALSE)
  pr4 <- prior(des4, type = "standard", pmean = pmean, psd = psd)
  emc4 <- make_emc(dat, des4, type = "standard", prior = pr4, compress = TRUE, n_chains = 1)
  expect_equal(emc4[[1]]$marginalise$n_nodes, 5L)
})
