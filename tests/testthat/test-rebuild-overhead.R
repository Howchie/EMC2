local_rng_guard()

test_that("rebuild claims preserve order, NULL values and subject streams", {
  skip_on_os("windows")
  set.seed(71)
  seeds <- .emc_subject_streams(31L)
  X <- stats::setNames(seq_along(seeds), paste0("subject", seq_along(seeds)))
  work <- function(i) {
    assign(".Random.seed", seeds[[i]], envir = globalenv())
    if (i == 7L) return(NULL)
    list(draw = rnorm(4), seed = .Random.seed)
  }
  expected <- .emc_with_preserved_rng(lapply(X, work))
  before <- .Random.seed
  for (i in seq_len(12L)) {
    got <- .emc_with_preserved_rng(.emc_rebuild_lapply(X, work, 2L, rev(seq_along(X))))
    expect_identical(got, expected)
    expect_identical(.Random.seed, before)
    expect_false(.emc_rebuild_queue$busy)
  }
})

test_that("rebuild queues propagate worker errors and can be reused", {
  skip_on_os("windows")
  expect_error(suppressWarnings(.emc_rebuild_lapply(1:12, function(i) {
    if (i == 3L) stop("rebuild failed")
    i
  }, 2L)), "rebuild failed")
  expect_false(.emc_rebuild_queue$busy)
  expect_identical(.emc_rebuild_lapply(1:12, identity, 2L), as.list(1:12))
})

test_that("single-core and short rebuilds keep the serial path", {
  expect_identical(.emc_rebuild_lapply(integer(), identity, 1L), list())
  expect_identical(.emc_rebuild_lapply(1:3, identity, 1L), as.list(1:3))
})


test_that("direct array slices match the split index implementation", {
  previous <- function(a, n) setNames(lapply(
    split(a, arrayInd(seq_along(a), dim(a))[, n]),
    array, dim = dim(a)[-n], dimnames(a)[-n]), dimnames(a)[[n]])
  for (dims in list(c(3L, 4L), c(2L, 3L, 5L), c(1L, 4L, 2L), c(2L, 1L, 3L, 4L))) {
    a <- array(seq_len(prod(dims)), dims)
    for (named in c(FALSE, TRUE)) {
      if (named) dimnames(a) <- lapply(dims, function(n) paste0("x", seq_len(n)))
      for (n in seq_along(dims)) expect_identical(split_along_dim(a, n), previous(a, n))
    }
  }
  a <- array(c(NA_real_, NaN, Inf, -Inf, -0, 1), c(1L, 2L, 3L))
  for (n in 1:3) expect_identical(split_along_dim(a, n), previous(a, n))
})
