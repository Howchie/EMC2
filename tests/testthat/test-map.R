RNGkind("L'Ecuyer-CMRG")
set.seed(123)

dat <- forstmann
dat$covariate <- 1:nrow(forstmann)

des <- design(data = dat, formula = list(v ~ covariate*E, B ~ E, t0 ~ S),
              model = LBA)

test_that("mapped_pars", {
  expect_snapshot(mapped_pars(des))
  expect_snapshot(mapped_pars(des, p_vector= rnorm(length(sampled_pars(des)))))
  expect_snapshot(mapped_pars(prior(des, mu_mean = c('v_covariate'  = 1))))
  expect_snapshot(mapped_pars(samples_LNR))
  expect_snapshot(mapped_pars(get_prior(samples_LNR)))
  expect_snapshot(mapped_pars(get_design(samples_LNR)))
})

test_that("map", {
  expect_snapshot(credint(samples_LNR, selection = "mu", map = "E"))
  expect_snapshot(credint(samples_LNR, selection = "mu", map = list(~ E*S)))
  expect_snapshot(credint(samples_LNR, selection = "mu", map = TRUE))
})

test_that("batched mapped draws preserve mapped summaries", {
  old_batch_size <- getOption("EMC2.map_batch_size")
  on.exit(options(EMC2.map_batch_size = old_batch_size), add = TRUE)

  options(EMC2.map_batch_size = 1L)
  one_draw <- get_pars(samples_LNR, selection = "mu", map = TRUE,
                       length.out = 8, merge_chains = TRUE)
  options(EMC2.map_batch_size = 4L)
  batched <- get_pars(samples_LNR, selection = "mu", map = TRUE,
                      length.out = 8, merge_chains = TRUE)

  expect_identical(batched, one_draw)
})

test_that("credible-interval margins accommodate rotated labels", {
  short <- .credint_label_mar("m_short", 90)
  long <- .credint_label_mar("m_a_much_longer_mapped_parameter_label", 90)

  expect_gte(short[1L], 4.1)
  expect_gt(long[1L], short[1L])
})
