# emc2.independent_chains: each chain's proposals AND sub-block partition come
# from its own draws.

test_that("independent chains block on their own covariances", {
  local_mocked_bindings(create_chain_proposals = function(emc, ...) emc,
                        .package = "EMC2")
  pn <- c("a", "b", "c", "d")
  cor_pairs <- function(pairs) {
    R <- diag(4); dimnames(R) <- list(pn, pn)
    for (pr in pairs) R[pr[1], pr[2]] <- R[pr[2], pr[1]] <- 0.95
    R
  }
  chain <- function(R) {
    d <- list(list(), list())
    attr(d, "shared_ll_idx") <- rep(1L, 4)
    list(chains_var = list(R, R), data = d,
         samples = list(stage = rep("burn", 3), idx = 3L))
  }
  emc <- list(chain(cor_pairs(list(c("a", "b"), c("c", "d")))),
              chain(cor_pairs(list(c("a", "c"), c("b", "d")))))
  withr::local_options(list(emc2.independent_chains = TRUE))
  out <- EMC2:::add_proposals(emc, "burn", 1, n_blocks = 2)
  c1 <- attr(out[[1]]$data, "components")
  c2 <- attr(out[[2]]$data, "components")
  expect_identical(unname(c1[["a"]] == c1[["b"]]), TRUE)
  expect_identical(unname(c2[["a"]] == c2[["c"]]), TRUE)
  expect_false(identical(unname(c1), unname(c2)))
  withr::local_options(list(emc2.independent_chains = FALSE))
  shared <- EMC2:::add_proposals(emc, "burn", 1, n_blocks = 2)
  expect_identical(attr(shared[[1]]$data, "components"),
                   attr(shared[[2]]$data, "components"))
})

test_that("the ensemble scratch store keeps static matrices whole", {
  n <- 3L
  samples <- list(stage = rep("sample", n), idx = 2L,
                  theta_mu = matrix(1:6, 2, n),
                  last_theta_var_inv = diag(3))
  st <- EMC2:::.emc_ensemble_store(samples)
  expect_identical(st$last_theta_var_inv, diag(3))
  expect_identical(st$theta_mu, matrix(3:4, 2, 1))
})

test_that("chains build their own proposal windows by default", {
  withr::local_options(list(emc2.independent_chains = NULL, emc2.group_move = "ensemble"))
  expect_true(EMC2:::.emc_independent_chains())
  emc <- EMC2:::restore_duplicates(samples_LNR)
  for (j in seq_along(emc)) {
    emc[[j]]$samples$stage[] <- "sample"
    emc[[j]]$burn_mode_var <- NULL
    emc[[j]]$burn_mode_mu <- NULL
  }
  before <- EMC2:::create_chain_proposals(emc, samples_idx = 1:50, stage = "sample")
  other <- emc
  other[[2]]$samples$alpha <- other[[2]]$samples$alpha + 100
  after <- EMC2:::create_chain_proposals(other, samples_idx = 1:50, stage = "sample")
  expect_identical(before[[1]]$chains_mu, after[[1]]$chains_mu)
  expect_identical(before[[1]]$chains_var, after[[1]]$chains_var)
  expect_identical(before[[1]]$ensemble_q, after[[1]]$ensemble_q)
  expect_false(identical(before[[2]]$chains_mu, after[[2]]$chains_mu))
})

test_that("efficient proposals condition on the same chain's draws", {
  withr::local_options(list(emc2.independent_chains = NULL))
  withr::local_seed(sum(utf8ToInt("independent-conditional-proposals")))
  emc <- EMC2:::restore_duplicates(samples_LNR)
  set.seed(sum(utf8ToInt("independent-conditional-proposals")))
  before <- EMC2:::create_eff_proposals(emc, n_cores = 1)
  other <- emc
  other[[2]]$samples$alpha <- other[[2]]$samples$alpha + 100
  set.seed(sum(utf8ToInt("independent-conditional-proposals")))
  after <- EMC2:::create_eff_proposals(other, n_cores = 1)
  expect_identical(before[[1]]$eff_mu, after[[1]]$eff_mu)
  expect_identical(before[[1]]$eff_var, after[[1]]$eff_var)
  expect_false(identical(before[[2]]$eff_mu, after[[2]]$eff_mu))
})

test_that("incomplete objects recover configuration without borrowing adaptive state", {
  first <- samples_LNR[[1]]
  first$shear <- list(log_scale = 3)
  second <- samples_LNR[[2]]
  second$rng <- list(gibbs = 1:7)
  prior <- second$prior
  samples <- second$samples
  prop_var <- attr(second, "prop_var")
  expect_warning(out <- EMC2:::restore_duplicates(list(first, second), warn = TRUE),
                 "lacks chain-local adaptive state")
  # Summaries of such an object do not warn; only resuming sampling does.
  expect_silent(EMC2:::restore_duplicates(list(first, second)))
  expect_identical(out[[2]]$type, first$type)
  expect_identical(out[[2]]$data, first$data)
  expect_identical(out[[2]]$prior, prior)
  expect_identical(out[[2]]$samples, samples)
  expect_identical(out[[2]]$rng, second$rng)
  expect_identical(attr(out[[2]], "prop_var"), prop_var)
  expect_null(out[[2]]$shear)
  expect_null(out[[2]]$chains_var)
  expect_true(out[[2]]$init)
})
