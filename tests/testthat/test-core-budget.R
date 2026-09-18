# Exercise the bounded allocator's explicit spare-capacity option. The sampler's
# default route leaves automatic donation disabled until a combined executor is
# available.
budget <- function(...) {
  EMC2:::.particle_core_budget(..., blas_threads = 1L, allow_spare = TRUE)
}

test_that("spare cores are used when subjects are the binding constraint", {
  # The audit's own example.
  b <- budget(8L, n_cores = 32L, r_cores = 1L, total_cores = 32L)
  expect_identical(b$subject, 8L)
  expect_gt(b$likelihood, 1L)
  expect_identical(b$subject * b$likelihood, 32L)
})

test_that("the budget is never exceeded", {
  for (n_subj in c(1L, 2L, 8L, 32L, 140L)) {
    for (total in c(1L, 2L, 4L, 8L, 16L, 32L, 64L)) {
      for (r in c(1L, 2L, 4L)) {
        b <- budget(n_subj, n_cores = total, r_cores = r, total_cores = total)
        lbl <- sprintf("subjects %d, total %d, r %d", n_subj, total, r)
        expect_gte(b$subject, 1L, label = lbl)
        expect_gte(b$likelihood, 1L, label = lbl)
        expect_lte(b$subject * b$likelihood, total, label = lbl)
        # Never more subject workers than there are subjects to give them.
        expect_lte(b$subject, max(1L, n_subj), label = lbl)
      }
    }
  }
})

test_that("the budget is counted in threads, not processes", {
  # Eight workers each running a four-threaded BLAS is thirty-two threads. A
  # budget that counted only processes would oversubscribe by that factor.
  single <- budget(32L, n_cores = 32L, r_cores = 1L, total_cores = 32L)
  four <- EMC2:::.particle_core_budget(32L, n_cores = 32L, r_cores = 1L,
                                       total_cores = 32L, blas_threads = 4L,
                                       allow_spare = TRUE)
  expect_identical(single$subject * single$likelihood, 32L)
  expect_lte(four$subject * four$likelihood * 4L, 32L)
  expect_lt(four$subject, single$subject)
  # The BLAS width is an input to the budget, not part of its answer: callers
  # have always received exactly two elements and still do.
  expect_identical(names(four), c("subject", "likelihood"))
})

test_that("a single subject still gets its particles split", {
  b <- budget(1L, n_cores = 16L, r_cores = 1L, total_cores = 16L)
  expect_identical(b$subject, 1L)
  expect_gt(b$likelihood, 1L)
  expect_lte(b$likelihood, 16L)
})

test_that("more subjects than cores gives every core a subject", {
  b <- budget(140L, n_cores = 8L, r_cores = 1L, total_cores = 8L)
  expect_identical(b$subject, 8L)
  expect_identical(b$likelihood, 1L)
})

test_that("without a global budget the old contract is unchanged", {
  # The no-total_cores path is what callers outside the sampler use, and it
  # must keep returning r_cores as the inner width.
  b <- EMC2:::.particle_core_budget(8L, n_cores = 4L, r_cores = 2L)
  expect_identical(b$subject, 4L)
  expect_identical(b$likelihood, 2L)
  b1 <- EMC2:::.particle_core_budget(1L, n_cores = 4L, r_cores = 1L)
  expect_identical(b1$subject, 1L)
  expect_identical(b1$likelihood, 4L)
})

test_that("the thread policy is set where the library supports it", {
  # Not by exporting an environment variable after the library has initialised,
  # which is read once at load and silently ignored afterwards.
  skip_if_not_installed("RhpcBLASctl")
  previous <- EMC2:::.emc_set_blas_threads(2L)
  skip_if(is.null(previous), "no supported way to set thread policy here")
  on.exit(EMC2:::.emc_set_blas_threads(previous), add = TRUE)
  expect_identical(RhpcBLASctl::blas_get_num_procs(), 2L)
  EMC2:::.emc_set_blas_threads(previous)
  expect_identical(RhpcBLASctl::blas_get_num_procs(), previous)
})

test_that("an absent thread library is assumed to be single-threaded", {
  # A wrong guess here costs throughput, never correctness, so the absence of
  # the optional package must not be an error.
  n <- with_mocked_bindings(EMC2:::.emc_blas_threads(),
                            requireNamespace = function(...) FALSE,
                            .package = "base")
  expect_gte(n, 1L)
})

test_that("the mcmapply fallback keeps exactly the inner width it asked for", {
  # Each of that route's children draws its own proposals, and an inner
  # mclapply inside the likelihood advances the child's L'Ecuyer stream before
  # it draws again to choose a particle. Widening the inner budget there would
  # move the sampler's draws, so a fit would stop being reproducible from a
  # fixed seed for a reason unrelated to the model.
  spare <- budget(8L, n_cores = 32L, r_cores = 1L, total_cores = 32L)
  strict <- EMC2:::.particle_core_budget(8L, n_cores = 32L, r_cores = 1L,
                                         total_cores = 32L, blas_threads = 1L)
  expect_gt(spare$likelihood, 1L)
  expect_identical(strict$likelihood, 1L)
  expect_identical(strict$subject, spare$subject)
  # An explicit r_cores is still honoured on that route -- the caller asked.
  asked <- EMC2:::.particle_core_budget(8L, n_cores = 32L, r_cores = 4L,
                                        total_cores = 32L, blas_threads = 1L)
  expect_gte(asked$likelihood, 4L)
})

test_that("the default keeps the requested inner width", {
  b <- EMC2:::.particle_core_budget(3L, n_cores = 8L, r_cores = 1L,
                                    total_cores = 8L, blas_threads = 1L)
  expect_identical(b$likelihood, 1L)
  expect_identical(b$subject, 3L)
  # And the return shape callers have always had.
  expect_identical(names(b), c("subject", "likelihood"))
})
