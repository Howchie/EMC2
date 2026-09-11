# C14: joining blocks of sample history.
#
# `abind()` was doing this for every element of the sample store. For the big
# numeric arrays the join is a memcpy of two contiguous prefixes, because
# extending the LAST dimension leaves both operands contiguous. The replacement
# has to be indistinguishable from abind for every shape the store holds.

test_that("the join matches abind for every shape the sample store holds", {
  set.seed(11)
  cases <- list(
    alpha    = list(array(rnorm(3 * 4 * 5), c(3, 4, 5)),
                    array(rnorm(3 * 4 * 2), c(3, 4, 2))),
    subj_ll  = list(matrix(rnorm(4 * 5), 4, 5), matrix(rnorm(4 * 3), 4, 3)),
    theta_mu = list(matrix(rnorm(3 * 5), 3, 5), matrix(rnorm(3 * 1), 3, 1)),
    theta_var = list(array(rnorm(3 * 3 * 5), c(3, 3, 5)),
                     array(rnorm(3 * 3 * 4), c(3, 3, 4))),
    # An empty block: a stage can be extended by nothing.
    empty    = list(array(rnorm(2 * 2 * 3), c(2, 2, 3)),
                    array(numeric(0), c(2, 2, 0)))
  )
  for (nm in names(cases)) {
    a <- cases[[nm]][[1L]]
    b <- cases[[nm]][[2L]]
    expect_identical(EMC2:::.emc_join_iteration(a, b), abind::abind(a, b),
                     info = nm)
  }
})

test_that("non-array elements keep abind's behaviour exactly", {
  # The stage vector, scalar indices and anything else in the store are not
  # where the time is, and the rest of the file depends on abind's semantics
  # for them.
  for (pair in list(list("sample", "sample"),
                    list(c("burn", "burn"), c("adapt")),
                    list(5, 3),
                    list(1:3, 4:5))) {
    expect_identical(EMC2:::.emc_join_iteration(pair[[1L]], pair[[2L]]),
                     abind::abind(pair[[1L]], pair[[2L]]))
  }
})

test_that("dimnames survive the join the way abind leaves them", {
  a <- array(rnorm(2 * 3 * 2), c(2, 3, 2),
             dimnames = list(c("v", "B"), paste0("s", 1:3), c("i1", "i2")))
  b <- array(rnorm(2 * 3 * 1), c(2, 3, 1),
             dimnames = list(c("v", "B"), paste0("s", 1:3), "i3"))
  got <- EMC2:::.emc_join_iteration(a, b)
  want <- abind::abind(a, b)
  expect_identical(got, want)
  expect_identical(dimnames(got)[[1L]], c("v", "B"))
  expect_identical(dimnames(got)[[3L]], c("i1", "i2", "i3"))
  # Parameter names are what `get_pars` keys on, so losing them would be a
  # silent renaming of the posterior.
  a2 <- a; dimnames(a2) <- NULL
  expect_identical(EMC2:::.emc_join_iteration(a2, b), abind::abind(a2, b))
})

test_that("mismatched leading dimensions fall back rather than corrupt", {
  # A shape the fast path cannot handle must go to abind and fail the same way,
  # not silently produce an array of the wrong size.
  a <- array(rnorm(2 * 3 * 2), c(2, 3, 2))
  b <- array(rnorm(3 * 3 * 2), c(3, 3, 2))
  expect_error(EMC2:::.emc_join_iteration(a, b))
})

test_that("the join preserves values exactly, not approximately", {
  set.seed(12)
  a <- array(rnorm(5 * 7 * 11), c(5, 7, 11))
  b <- array(rnorm(5 * 7 * 4), c(5, 7, 4))
  got <- EMC2:::.emc_join_iteration(a, b)
  expect_bit_identical(as.numeric(got[, , 1:11]), as.numeric(a), "first block")
  expect_bit_identical(as.numeric(got[, , 12:15]), as.numeric(b), "second block")
  # Including the awkward values a sample store really carries.
  a[1, 1, 1] <- NA_real_; a[2, 2, 2] <- NaN; a[3, 3, 3] <- Inf; a[4, 4, 4] <- -Inf
  expect_identical(EMC2:::.emc_join_iteration(a, b), abind::abind(a, b))
})

test_that("the join takes any number of blocks, as merge_chains needs", {
  # `merge_chains` folds one array per CHAIN, not two, and a two-argument join
  # silently broke it -- caught by a multi-block gate rather than by a unit
  # test, which is why this one exists.
  set.seed(13)
  parts3 <- list(array(rnorm(2 * 3 * 4), c(2, 3, 4)),
                 array(rnorm(2 * 3 * 2), c(2, 3, 2)),
                 array(rnorm(2 * 3 * 5), c(2, 3, 5)))
  expect_identical(do.call(EMC2:::.emc_join_iteration, parts3),
                   do.call(abind::abind, parts3))
  parts5 <- replicate(5L, matrix(rnorm(3 * 2), 3, 2), simplify = FALSE)
  expect_identical(do.call(EMC2:::.emc_join_iteration, parts5),
                   do.call(abind::abind, parts5))
  # One operand, and none.
  expect_identical(EMC2:::.emc_join_iteration(parts3[[1L]]), parts3[[1L]])
  expect_null(EMC2:::.emc_join_iteration())
  # Mixed naming across three, which is where the axis-by-axis rule shows.
  a <- array(rnorm(8), c(2, 2, 2), dimnames = list(c("v", "B"), NULL, NULL))
  b <- array(rnorm(4), c(2, 2, 1))
  cc <- array(rnorm(4), c(2, 2, 1), dimnames = list(NULL, c("x", "y"), "i9"))
  expect_identical(EMC2:::.emc_join_iteration(a, b, cc), abind::abind(a, b, cc))
})
