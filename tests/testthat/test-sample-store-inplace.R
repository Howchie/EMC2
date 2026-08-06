test_that("native sample writers fill the public array representation", {
  alpha <- array(NA_real_, c(2, 3, 4))
  ll <- matrix(NA_real_, 3, 4)
  proposals <- rbind(matrix(1:6, 2, 3), c(-1, -2, -3))

  EMC2:::emc_set_particle_slice(alpha, ll, proposals, 3L, 2L)

  expect_equal(alpha[, , 3], proposals[1:2, ])
  expect_equal(ll[, 3], proposals[3, ])
  expect_true(all(is.na(alpha[, , -3])))

  theta <- array(NA_real_, c(2, 2, 4))
  EMC2:::emc_set_last_slice(theta, 2L, as.numeric(diag(2)))
  expect_equal(theta[, , 2], diag(2))
  expect_true(all(is.na(theta[, , -2])))
})

test_that("initial sample stores detach replicated chain references", {
  shared <- list(alpha = array(NA_real_, c(2, 2, 1)),
                 subj_ll = matrix(NA_real_, 2, 1), idx = 0L)
  chain_1 <- shared
  chain_2 <- shared
  chain_1 <- EMC2:::emc_clone_sample_store(chain_1)
  proposals <- rbind(matrix(1:4, 2, 2), c(-1, -2))

  EMC2:::emc_set_particle_slice(chain_1$alpha, chain_1$subj_ll,
                                proposals, 1L, 2L)

  expect_equal(chain_1$alpha[, , 1], proposals[1:2, ])
  expect_true(all(is.na(chain_2$alpha)))
  expect_true(all(is.na(chain_2$subj_ll)))
})

test_that("block extension detaches history before native writes", {
  original <- array(seq_len(12), c(2, 3, 2))
  old_binding <- original
  extended <- EMC2:::extend_obj(original, 2L)

  EMC2:::emc_set_last_slice(extended, 3L, rep(99, 6))

  expect_identical(original, old_binding)
  expect_equal(extended[, , 3], matrix(99, 2, 3))
  expect_equal(extended[, , 1:2], original)
})

test_that("block extension copies the old contiguous prefix", {
  original <- array(seq_len(24), c(2, 3, 4),
                    dimnames = list(c("a", "b"), letters[1:3], NULL))
  extended <- EMC2:::extend_obj(original, 3L)

  expect_equal(dim(extended), c(2, 3, 7))
  expect_equal(extended[, , 1:4], original)
  expect_true(all(is.na(extended[, , 5:7])))
  expect_identical(dimnames(extended)[1:2], dimnames(original)[1:2])
})
