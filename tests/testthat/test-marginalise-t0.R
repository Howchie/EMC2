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
