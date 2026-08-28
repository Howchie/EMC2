test_that("pointwise likelihood workers are one-shot and match serial results", {
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("mkfifo")))

  data(forstmann)
  ids <- levels(forstmann$subjects)[1:3]
  dat <- forstmann[forstmann$subjects %in% ids, ]
  dat$subjects <- droplevels(dat$subjects)
  des <- design(data = dat, model = LBA,
                formula = list(v ~ 1, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1))
  emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1))

  n_iter <- 4L
  n_subj <- length(emc[[1]]$data)
  emc[[1]]$samples$alpha <- array(
    rep(c(1, 1, .5, .5, .2), n_subj * n_iter),
    dim = c(5L, n_subj, n_iter),
    dimnames = list(c("v", "sv", "B", "A", "t0"),
                    names(emc[[1]]$data), NULL))
  emc[[1]]$samples$stage <- rep("sample", n_iter)
  emc[[1]]$samples$idx <- n_iter
  emc[[1]]$samples$subj_ll <- matrix(0, n_subj, n_iter)

  serial <- EMC2:::.ll_matrix_pooled(emc, cores = 1)
  parallel <- EMC2:::.ll_matrix_pooled(emc, cores = 2)
  expect_equal(parallel, serial, tolerance = 1e-12)
})

test_that("LOO runs in an isolated wrapper when multiple cores are requested", {
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("mkfifo")))

  set.seed(11)
  ll <- matrix(rnorm(2000), nrow = 100L, ncol = 20L)
  serial <- suppressWarnings(loo::loo(ll, cores = 1L)$estimates[
    "looic", "Estimate"])
  parallel <- suppressWarnings(EMC2:::loo_from_ll(ll, cores = 2L))
  expect_equal(parallel, serial, tolerance = 1e-12)
})
