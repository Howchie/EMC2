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

test_that("mapped_pars can restrict mappings to observed factor combinations", {
  dat <- data.frame(
    subjects = factor(rep("s1", 3)),
    block = factor(c("b1", "b1", "b2"), levels = c("b1", "b2")),
    E = factor(c("easy", "hard", "hard"), levels = c("easy", "hard")),
    R = factor(c("left", "right", "left"), levels = c("left", "right")),
    rt = c(0.5, 0.6, 0.7)
  )
  des <- design(
    data = dat,
    model = DDM,
    formula = list(v ~ 0 + block:E, a ~ 1, t0 ~ 1, s ~ 1,
                   Z ~ 1, sv ~ 1, SZ ~ 1),
    report_p_vector = FALSE
  )
  p_vec <- sampled_pars(des)

  all_combinations <- mapped_pars(des, p_vec)
  observed_combinations <- mapped_pars(des, p_vec, data = dat)

  expect_equal(nrow(all_combinations), 4)
  expect_equal(nrow(observed_combinations), 3)
  expect_false(any(observed_combinations$block == "b2" &
                     observed_combinations$E == "easy"))

  emc <- make_emc(dat, des, type = "single", compress = FALSE, n_chains = 1)
  expect_equal(mapped_pars(emc, p_vec), observed_combinations)
  expect_equal(nrow(mapped_pars(emc, p_vec, use_data = FALSE)),
               nrow(all_combinations))

  expected_sampled <- sampled_pars(des, data = dat)
  expect_identical(names(sampled_pars(emc)), names(expected_sampled))
  expect_false(any(grepl("blockb2:Eeasy", names(sampled_pars(emc)), fixed = TRUE)))
  expect_true(any(grepl("blockb2:Eeasy",
                         names(sampled_pars(emc, use_data = FALSE)), fixed = TRUE)))

  # A saved EMC object may predate the data attribute on its design.  The
  # emc method must still use the data stored in the object by default.
  legacy_emc <- emc
  legacy_design <- attr(legacy_emc[[1]]$prior, "design")
  legacy_design <- lapply(legacy_design, function(x) {
    attr(x, "data") <- NULL
    x
  })
  attr(legacy_emc[[1]]$prior, "design") <- legacy_design
  expect_identical(names(sampled_pars(legacy_emc)), names(expected_sampled))

  emc_sampled <- attr(emc[[1]]$data[[1]], "sampled_p_names")
  expect_false(any(grepl("blockb2:Eeasy", emc_sampled, fixed = TRUE)))

  des_with_unused_constant <- des
  des_with_unused_constant$constants <- c(
    des_with_unused_constant$constants,
    v_blockb2.Eeasy = 0
  )
  names(des_with_unused_constant$constants)[
    names(des_with_unused_constant$constants) == "v_blockb2.Eeasy"
  ] <- "v_blockb2:Eeasy"
  emc_with_unused_constant <- make_emc(
    dat, des_with_unused_constant, type = "single", compress = FALSE,
    n_chains = 1
  )
  expect_false(any(grepl(
    "blockb2:Eeasy",
    attr(emc_with_unused_constant[[1]]$data[[1]], "p_names"),
    fixed = TRUE
  )))

  emc_full <- make_emc(dat, des, type = "single", compress = FALSE,
                       n_chains = 1, use_data = FALSE)
  emc_full_p <- attr(emc_full[[1]]$data[[1]], "p_names")
  expect_true(any(grepl("blockb2:Eeasy", emc_full_p, fixed = TRUE)))
})

test_that("data-aware make_emc keeps a shared design for subjects with missing cells", {
  dat <- data.frame(
    subjects = factor(c("s1", "s1", "s2")),
    block = factor(c("b1", "b1", "b2"), levels = c("b1", "b2")),
    E = factor(c("easy", "hard", "hard"), levels = c("easy", "hard")),
    R = factor(c("left", "right", "left"), levels = c("left", "right")),
    rt = c(0.5, 0.6, 0.7)
  )
  des <- design(
    data = dat,
    model = DDM,
    formula = list(v ~ 0 + block:E, a ~ 1, t0 ~ 1, s ~ 1,
                   Z ~ 1, sv ~ 1, SZ ~ 1),
    report_p_vector = FALSE
  )

  emc <- make_emc(dat, des, type = "standard", compress = FALSE,
                  n_chains = 1)
  s2 <- emc[[1]]$data[["s2"]]
  v_dm <- attr(s2, "designs")$v
  v_used <- v_dm[attr(v_dm, "expand"), , drop = FALSE]

  expect_true("v_blockb1:Eeasy" %in% attr(emc[[1]]$data[[1]], "sampled_p_names"))
  expect_false(any(v_used[, "v_blockb1:Eeasy"] != 0))
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
