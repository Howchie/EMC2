# par_groups labels need not be 1..k: c(2, 2, 5, 5) is two covariance blocks.
# The sampler used to compare labels with table() positions, so the second
# block was silently treated as diagonal.

test_that("non-contiguous par_groups labels block like contiguous ones", {
  base <- list(nuisance = rep(FALSE, 4), par_names = paste0("p", 1:4))
  info <- function(pg) EMC2:::add_info_standard(base, par_groups = pg)$is_blocked
  expect_identical(info(c(2, 2, 5, 5)), rep(TRUE, 4))
  expect_identical(info(c(2, 2, 5, 7)), c(TRUE, TRUE, FALSE, FALSE))
  expect_identical(info(c(1, 1, 2, 3)), c(TRUE, TRUE, FALSE, FALSE))
})

test_that("the prior sampler zeroes only cross-block covariances", {
  set.seed(1)
  pr <- EMC2:::get_prior_standard(n_pars = 4, sample = TRUE, N = 50,
                                  selection = "sigma2", par_groups = c(2, 2, 5, 5),
                                  design = NULL)
  v <- pr$theta_var
  expect_true(all(v[1, 2, ] != 0))
  expect_true(all(v[3, 4, ] != 0))
  expect_true(all(v[1, 3, ] == 0))
})
