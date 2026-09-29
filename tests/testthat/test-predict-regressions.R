local_rng_guard()

test_that("predict keeps a free parameter that is constant across posterior draws", {
  # get_pars() drops rows without variation across draws by default. A free
  # parameter can be constant (e.g. integer sampling); the subject path then
  # lost it and the hyper path paired a short mu with a full Sigma.
  set.seed(5)
  emc <- EMC2:::restore_duplicates(samples_LNR)
  for (j in seq_along(emc)) emc[[j]]$samples$alpha["s", , ] <- -0.5
  sims <- suppressMessages(predict(emc, n_post = 3))
  expect_s3_class(sims, "data.frame")
  expect_gt(nrow(sims), 0L)

  hyper <- EMC2:::restore_duplicates(samples_LNR)
  for (j in seq_along(hyper)) hyper[[j]]$samples$theta_mu["s", ] <- -0.5
  sims <- suppressMessages(predict(hyper, n_post = 3, hyper = TRUE))
  expect_s3_class(sims, "data.frame")
})

test_that("design constants never reach the matrix predict() draws from", {
  # predict() takes sampled parameters (map = FALSE); design constants are
  # added later by the mapper. So keeping constant rows changes nothing for a
  # model with no constant draws. (Checked 2026-09-29: predictions identical
  # before and after on RDM with constants s, A and LBA with constant sv.)
  for (selection in c("alpha", "mu")) {
    keep <- get_pars(samples_LNR, selection = selection, map = FALSE,
                     return_mcmc = FALSE, merge_chains = TRUE,
                     remove_constants = FALSE)
    default <- get_pars(samples_LNR, selection = selection, map = FALSE,
                        return_mcmc = FALSE, merge_chains = TRUE)
    expect_identical(keep, default)
  }
})

test_that("joint models return their data and predictions", {
  # get_data() named the joint data via sampled_pars(emc), which itself calls
  # get_data(): a C stack overflow for every joint emc.
  design_in <- get_design(samples_LNR)[[1]]
  data_in <- get_data(samples_LNR)
  joint <- suppressMessages(make_emc(
    list(data_in, data_in), list(a = design_in, b = design_in),
    prior_list = prior(list(a = design_in, b = design_in)), n_chains = 1))
  dat <- get_data(joint)
  expect_named(dat, c("a", "b"))
  expect_identical(nrow(dat$a), nrow(data_in))
  expect_true(all(grepl("^(a|b)[|]", names(sampled_pars(joint)))))
})
