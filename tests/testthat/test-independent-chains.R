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
