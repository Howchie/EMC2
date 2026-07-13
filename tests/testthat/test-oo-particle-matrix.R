test_that("particle matrices restore dropped zero columns and constants", {
  design_mat <- matrix(c(1, 1, 0, 0), nrow = 2, ncol = 2,
                       dimnames = list(NULL, c("active", "dropped")))
  dadm <- structure(data.frame(row = 1:2),
                    sampled_p_names = "active",
                    p_names = c("active", "dropped", "fixed"),
                    constants = c(fixed = 2),
                    designs = list(main = design_mat))

  particles <- matrix(c(3, 4), nrow = 2, ncol = 1,
                      dimnames = list(NULL, "active"))
  out <- EMC2:::.oo_particle_matrix(particles, dadm,
                                    keep_all_columns = TRUE)

  expect_identical(colnames(out), c("active", "dropped", "fixed"))
  expect_equal(out[, "active"], c(3, 4))
  expect_equal(out[, "dropped"], c(0, 0))
  expect_equal(out[, "fixed"], c(2, 2))
})

test_that("particle matrices still reject a missing active design column", {
  dadm <- structure(data.frame(row = 1:2),
                    sampled_p_names = "active",
                    p_names = c("active", "missing"),
                    designs = list(main = matrix(c(1, 1), ncol = 1,
                                                 dimnames = list(NULL, "active"))))
  particles <- matrix(c(3, 4), nrow = 2, ncol = 1,
                      dimnames = list(NULL, "active"))

  expect_error(
    EMC2:::.oo_particle_matrix(particles, dadm, keep_all_columns = TRUE),
    "missing"
  )
})

test_that("the OO mapper consumes compressed design rows", {
  design_mat <- matrix(c(1, 2), nrow = 2, ncol = 1,
                       dimnames = list(NULL, "p"))
  attr(design_mat, "expand") <- c(1L, 2L, 1L, 2L)
  dadm <- structure(data.frame(row = 1:4),
                    sampled_p_names = "p",
                    p_names = "p",
                    designs = list(x = design_mat))
  identity <- list(func = c(p = "identity"),
                   lower = c(p = -Inf), upper = c(p = Inf))
  model <- list(transform = identity, pre_transform = identity,
                trend = NULL)

  out <- EMC2:::get_pars_oo(c(p = 3), dadm, model)
  expect_equal(as.numeric(out[, "x"]), c(3, 6, 3, 6))
})
