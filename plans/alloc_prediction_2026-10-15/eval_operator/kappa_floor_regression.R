# Exercise the preserved implementation independently of later working-tree edits.
suppressPackageStartupMessages(library(EMC2, lib.loc = ".Rlib_alloc_review"))
e <- new.env(parent = asNamespace("EMC2"))
sys.source("plans/alloc_prediction_2026-10-15/eval_operator/evaluated_allocation.R", e)
testthat::test_that("corrected kappa pricing rejects targets below its floor", {
  rs <- c(0.2, 0.6, 0.7, 0.95)
  reachable <- rs >= e$.emc_alloc_eff(0, 4)
  testthat::expect_equal(e$.emc_alloc_eff_inv(rs[!reachable], 4), -1)
  testthat::expect_equal(e$.emc_alloc_eff(e$.emc_alloc_eff_inv(rs[reachable], 4), 4), rs[reachable])
  u <- list(kind = "shear", rho = 0, g = 0.5, a = 0.5, cost = 1, floor = 0, cap = 4, kappa = 4)
  one <- e$.emc_alloc_unit_need(u, rs, 4)
  batch <- e$.emc_alloc_units_need(list(u), TRUE, rs, 4)[[1]]
  testthat::expect_equal(one$cost, c(Inf, 0, 0, 0))
  testthat::expect_equal(batch$cost, one$cost)
})
testthat::test_that("the stage tau estimator matches a long AR reference", {
  set.seed(11, kind = "L'Ecuyer-CMRG")
  ar <- function(n, rho) {
    x <- numeric(n)
    for (t in 2:n) x[t] <- rho * x[t - 1L] + rnorm(1, sd = sqrt(1 - rho^2))
    x
  }
  n <- 20000L
  mixed <- lapply(1:3, function(i) cbind(ar(n, 0.5), ar(n, 0.9)))
  tau <- e$.emc_alloc_geyer_tau(mixed)
  testthat::expect_equal(tau, c(3, 19), tolerance = 0.25)
  a <- simplify2array(mixed)
  ref <- vapply(1:2, function(k) n * 3 / posterior::ess_basic(a[, k, ]), numeric(1))
  testthat::expect_equal(tau, ref, tolerance = 0.05)
})
