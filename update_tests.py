import re

with open('tests/testthat/test-splitlognormal.R', 'r') as f:
    text = f.read()

# Add a new test for BAwL and extreme deltas.
test_to_add = """
testthat::test_that("extreme delta produces valid splitlognormal parameters without overflow", {
  skip_if_not(exists(".bawd_split_params", envir = asNamespace("EMC2"), inherits = FALSE),
             "split launch helpers are unavailable")
  split_params <- getFromNamespace(".bawd_split_params", "EMC2")
  
  mu <- 0
  sigma <- 1
  delta <- -1500
  
  sp <- split_params(mu, sigma, delta)
  expect_true(is.finite(sp$c))
  expect_true(is.finite(sp$sigma_L))
  expect_true(is.finite(sp$sigma_R))
  expect_true(is.finite(sp$a))
})

testthat::test_that("delta zero preserves BAwL lognormal path", {
  skip_if_not(exists("dbawl", envir = asNamespace("EMC2"), inherits = FALSE),
              "dbawl unavailable")
  expect_equal(EMC2:::dbawl(1, 0, 1, 0, 1, 0.5, launch = 1),
               EMC2:::dbawl(1, 0, 1, 0, 1, 0.5, launch = 2, delta = 0), tolerance = 1e-12)
  expect_equal(EMC2:::pbawl(1, 0, 1, 0, 1, 0.5, launch = 1),
               EMC2:::pbawl(1, 0, 1, 0, 1, 0.5, launch = 2, delta = 0), tolerance = 1e-12)
})
"""

text += test_to_add

with open('tests/testthat/test-splitlognormal.R', 'w') as f:
    f.write(text)
print("Tests appended")
