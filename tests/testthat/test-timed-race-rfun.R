test_that("race rfuns convert a winning time accumulator into a timed guess", {
  set.seed(123)
  n_trials <- 200
  lR <- factor(rep(c("left", "right", "time"), n_trials), levels = c("left", "right", "time"))

  lba_pars <- rbind(
    c(v = 1e-6, sv = 1e-8, b = 1.1, A = 0.0, t0 = 0.2),
    c(v = 1e-6, sv = 1e-8, b = 1.1, A = 0.0, t0 = 0.2),
    c(v = 8.0,  sv = 1e-8, b = 1.1, A = 0.0, t0 = 0.2)
  )
  lba_pars <- lba_pars[rep(seq_len(nrow(lba_pars)), n_trials), , drop = FALSE]
  sim_lba <- EMC2:::rLBA(lR, lba_pars)
  expect_false(any(as.character(sim_lba$R) == "time"))
  expect_true(all(as.character(sim_lba$R) %in% c("left", "right")))

  lnr_pars <- rbind(
    c(m = 2.0,  s = 1e-8, t0 = 0.2),
    c(m = 2.0,  s = 1e-8, t0 = 0.2),
    c(m = -8.0, s = 1e-8, t0 = 0.2)
  )
  lnr_pars <- lnr_pars[rep(seq_len(nrow(lnr_pars)), n_trials), , drop = FALSE]
  sim_lnr <- EMC2:::rLNR(lR, lnr_pars)
  expect_false(any(as.character(sim_lnr$R) == "time"))
  expect_true(all(as.character(sim_lnr$R) %in% c("left", "right")))

  rdm_pars <- rbind(
    c(v = 1e-6, B = 1.1, A = 0.0, t0 = 0.2),
    c(v = 1e-6, B = 1.1, A = 0.0, t0 = 0.2),
    c(v = 8.0,  B = 1.1, A = 0.0, t0 = 0.2)
  )
  rdm_pars <- rdm_pars[rep(seq_len(nrow(rdm_pars)), n_trials), , drop = FALSE]
  sim_rdm <- EMC2:::rRDM(lR, rdm_pars)
  expect_false(any(as.character(sim_rdm$R) == "time"))
  expect_true(all(as.character(sim_rdm$R) %in% c("left", "right")))
})

test_that("predict does not emit the synthetic time response level", {
  set.seed(456)

  timed_design <- design(
    factors = list(subjects = 1, S = "stim"),
    Rlevels = c("left", "right", "time"),
    matchfun = function(d) as.character(d$S) == as.character(d$lR),
    model = LBA,
    formula = list(v ~ lR, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1)
  )

  p_vec <- sampled_pars(timed_design, doMap = FALSE)
  p_vec["v_lRleft"] <- 1e-6
  p_vec["v_lRright"] <- 1e-6
  p_vec["v_lRtime"] <- 8.0
  p_vec["sv"] <- log(0.1)
  p_vec["B"] <- log(1.1)
  p_vec["A"] <- log(0.1)
  p_vec["t0"] <- log(0.2)
  p_sd <- setNames(rep(1e-8, length(p_vec)), names(p_vec))

  pr <- prior(timed_design, mu_mean = p_vec, mu_sd = p_sd)
  template <- data.frame(
    subjects = factor(rep(1, 40)),
    S = factor(rep("stim", 40), levels = "stim"),
    R = factor(rep(NA_character_, 40), levels = c("left", "right", "time")),
    rt = NA_real_
  )
  pp <- predict(pr, data = template, n_post = 2, n_cores = 1)

  expect_false(any(as.character(pp$R) == "time"))
  expect_true(all(as.character(stats::na.omit(pp$R)) %in% c("left", "right")))
})
