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

test_that("RDMSWTN scalar-rt wrappers vectorize over accumulator rows", {
  pars <- cbind(
    v = c(3, 1, 2),
    B = 1,
    A = 0,
    t0 = c(0.2, 0.2, 0.05),
    s = 1,
    sv = 0,
    lambda_g = 0,
    lambda_k = 0,
    pContaminant = 0,
    b = 1
  )

  expected_p <- EMC2:::pSWTNspv(
    rep(0.55, 3),
    v = pars[, "v"], b = pars[, "b"], A = pars[, "A"],
    s = pars[, "s"], t0 = pars[, "t0"], sv = pars[, "sv"],
    lambda_g = 0, lambda_k = 0
  )
  expected_d <- EMC2:::dSWTNspv(
    rep(0.55, 3),
    v = pars[, "v"], b = pars[, "b"], A = pars[, "A"],
    s = pars[, "s"], t0 = pars[, "t0"], sv = pars[, "sv"],
    lambda_g = 0, lambda_k = 0
  )

  expect_equal(EMC2:::pRDMSWTN(0.55, pars), expected_p)
  expect_equal(EMC2:::dRDMSWTN(0.55, pars), expected_d)
})

test_that("rSWTN samples from the same start-point convention as pSWTNspv", {
  set.seed(789)
  n <- 40000
  probs <- c(0.2, 0.5, 0.8)
  sim <- EMC2:::rSWTN(n, b = 1, v = 2, A = 0.3, sv = 0, k = 0)
  empirical <- unname(stats::quantile(sim, probs))
  expected <- vapply(probs, function(p) {
    stats::uniroot(
      function(t) EMC2:::pSWTNspv(t, v = 2, b = 1, A = 0.3, s = 1, t0 = 0,
                                  sv = 0, lambda_g = 0, lambda_k = 0) - p,
      c(1e-6, 5)
    )$root
  }, numeric(1))

  expect_equal(empirical, expected, tolerance = 0.02)
})

test_that("timed RDMSWTN likelihood matches manual mixture with truncation", {
  timed_design <- design(
    factors = list(subjects = 1, S = "stim"),
    Rlevels = c("left", "right", "time"),
    matchfun = function(d) as.character(d$S) == as.character(d$lR),
    model = RDMSWTN(),
    formula = list(v ~ 0 + lR, B ~ 1, A ~ 1, t0 ~ 0 + lR, s ~ 1, sv ~ 1),
    constants = c(lambda_g = log(0), lambda_k = log(0)),
    report_p_vector = FALSE
  )

  dat <- data.frame(
    subjects = factor(1),
    S = factor("stim"),
    R = factor("left", levels = c("left", "right", "time")),
    rt = 0.55,
    LT = 0.25,
    UT = 1.5
  )
  emc <- make_emc(dat, timed_design, type = "single", compress = FALSE, n_chains = 1)
  model_obj <- emc[[1]]$model()
  dadm <- emc[[1]]$data[[1]]

  p_vec <- sampled_pars(timed_design, doMap = FALSE)
  p_vec[] <- 0
  p_vec["v_lRleft"] <- log(3)
  p_vec["v_lRright"] <- log(1.2)
  p_vec["v_lRtime"] <- log(2)
  p_vec["B"] <- log(1)
  p_vec["A"] <- log(0)
  p_vec["t0_lRleft"] <- log(0.2)
  p_vec["t0_lRright"] <- log(0.2)
  p_vec["t0_lRtime"] <- log(0.05)
  p_vec["s"] <- log(1)
  p_vec["sv"] <- log(0)

  p_mat <- matrix(p_vec, nrow = 1, dimnames = list(NULL, names(p_vec)))
  designs <- lapply(names(model_obj$p_types), function(p) {
    attr(dadm, "designs")[[p]][attr(attr(dadm, "designs")[[p]], "expand"), , drop = FALSE]
  })
  names(designs) <- names(model_obj$p_types)
  constants <- attr(dadm, "constants")
  if (is.null(constants)) constants <- NA

  ll_cpp <- EMC2:::calc_ll_oo(
    p_mat, dadm,
    constants = constants,
    designs = designs,
    type = model_obj$c_name,
    bounds = model_obj$bound,
    transforms = model_obj$transform,
    pretransforms = model_obj$pre_transform,
    p_types = names(model_obj$p_types),
    min_ll = log(1e-10),
    trend = model_obj$trend
  )

  pars <- as.matrix(mapped_pars(timed_design, p_vec, data = dat)[, names(model_obj$p_types), drop = FALSE])
  pars <- cbind(pars, b = pars[, "B"] + pars[, "A"])
  rt <- 0.55
  f <- EMC2:::dSWTNspv(rep(rt, 3), pars[, "v"], pars[, "b"], pars[, "A"],
                       pars[, "s"], pars[, "t0"], pars[, "sv"], lambda_g = 0, lambda_k = 0)
  F_rt <- EMC2:::pSWTNspv(rep(rt, 3), pars[, "v"], pars[, "b"], pars[, "A"],
                          pars[, "s"], pars[, "t0"], pars[, "sv"], lambda_g = 0, lambda_k = 0)
  F_lt <- EMC2:::pSWTNspv(rep(0.25, 3), pars[, "v"], pars[, "b"], pars[, "A"],
                          pars[, "s"], pars[, "t0"], pars[, "sv"], lambda_g = 0, lambda_k = 0)
  F_ut <- EMC2:::pSWTNspv(rep(1.5, 3), pars[, "v"], pars[, "b"], pars[, "A"],
                          pars[, "s"], pars[, "t0"], pars[, "sv"], lambda_g = 0, lambda_k = 0)
  S_rt <- 1 - F_rt
  manual_density <- f[1] * S_rt[2] * S_rt[3] + f[3] * S_rt[1] * S_rt[2] / 2
  trunc_prob <- prod(1 - F_lt) - prod(1 - F_ut)
  expect_equal(as.numeric(ll_cpp), log(manual_density) - log(trunc_prob), tolerance = 5e-3)
})

test_that("RDMSWTN IO C++ likelihood keeps defective Wald normalization with truncation", {
  timed_design <- design(
    factors = list(subjects = 1, S = "stim"),
    Rlevels = c("left", "right"),
    matchfun = function(d) as.character(d$S) == as.character(d$lR),
    model = RDMSWTN(posdrift = FALSE),
    formula = list(v ~ 0 + lR, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1),
    constants = c(lambda_g = log(0), lambda_k = log(0), sv = log(0)),
    report_p_vector = FALSE
  )

  dat <- data.frame(
    subjects = factor(1),
    S = factor("stim"),
    R = factor("left", levels = c("left", "right")),
    rt = 0.7,
    LT = 0.2,
    UT = 2.0
  )
  emc <- make_emc(dat, timed_design, type = "single", compress = FALSE, n_chains = 1)
  model_obj <- emc[[1]]$model()
  dadm <- emc[[1]]$data[[1]]

  p_vec <- sampled_pars(timed_design, doMap = FALSE)
  p_vec[] <- 0
  p_vec["v_lRleft"] <- -0.7
  p_vec["v_lRright"] <- -0.4
  p_vec["B"] <- log(1.0)
  p_vec["A"] <- log(0.2)
  p_vec["t0"] <- log(0.1)
  p_vec["s"] <- log(1)

  p_mat <- matrix(p_vec, nrow = 1, dimnames = list(NULL, names(p_vec)))
  ll_cpp <- EMC2:::calc_ll_oo(
    p_mat, dadm,
    constants = attr(dadm, "constants"),
    designs = EMC2:::.oo_expanded_designs(dadm),
    type = model_obj$c_name,
    bounds = model_obj$bound,
    transforms = model_obj$transform,
    pretransforms = model_obj$pre_transform,
    p_types = names(model_obj$p_types),
    min_ll = log(1e-10),
    trend = model_obj$trend
  )

  pars <- as.matrix(mapped_pars(timed_design, p_vec, data = dat)[, names(model_obj$p_types), drop = FALSE])
  pars <- cbind(pars, b = pars[, "B"] + pars[, "A"])
  f <- EMC2:::dSWTNspv(rep(dat$rt, 2), pars[, "v"], pars[, "b"], pars[, "A"],
                       pars[, "s"], pars[, "t0"], pars[, "sv"],
                       lambda_g = 0, lambda_k = 0, posdrift = FALSE)
  F_rt <- EMC2:::pSWTNspv(rep(dat$rt, 2), pars[, "v"], pars[, "b"], pars[, "A"],
                          pars[, "s"], pars[, "t0"], pars[, "sv"],
                          lambda_g = 0, lambda_k = 0, posdrift = FALSE)
  F_lt <- EMC2:::pSWTNspv(rep(dat$LT, 2), pars[, "v"], pars[, "b"], pars[, "A"],
                          pars[, "s"], pars[, "t0"], pars[, "sv"],
                          lambda_g = 0, lambda_k = 0, posdrift = FALSE)
  F_ut <- EMC2:::pSWTNspv(rep(dat$UT, 2), pars[, "v"], pars[, "b"], pars[, "A"],
                          pars[, "s"], pars[, "t0"], pars[, "sv"],
                          lambda_g = 0, lambda_k = 0, posdrift = FALSE)
  manual_density <- f[1] * (1 - F_rt[2])
  trunc_prob <- prod(1 - F_lt) - prod(1 - F_ut)

  expect_equal(as.numeric(ll_cpp), log(manual_density) - log(trunc_prob), tolerance = 1e-8)
})

test_that("TRDM split transforms are preserved for C++ likelihood mapping", {
  skip_if_not_installed("dplyr")

  dat <- expand.grid(
    subjects = factor(1),
    sa = factor(c("speed", "neutral", "accuracy"), levels = c("speed", "neutral", "accuracy")),
    S = factor(c("left", "right"), levels = c("left", "right")),
    rep = 1:2
  )
  dat$R <- factor(dat$S, levels = c("left", "right", "time"))
  dat$rt <- seq(0.35, 0.9, length.out = nrow(dat))
  dat$LT <- 0.25
  dat$UT <- 1.5
  matchfun <- function(d) as.numeric(d$lR) == as.numeric(d$S)

  timed_design <- design(
    data = dat,
    factors = list(Rlevels = c("left", "right", "time")),
    model = RDMSWTN(),
    transform = list(func = c(v = "exp")),
    matchfun = matchfun,
    functions = list(
      match = function(d) {
        dplyr::case_when(
          d$lM == TRUE & !(d$lR == "time") ~ 0.5,
          d$lM == FALSE & !(d$lR == "time") ~ -0.5,
          d$lR == "time" ~ 0
        )
      },
      E = function(d) ifelse(d$lR == "time", 0, 1),
      Time = function(d) ifelse(d$lR == "time", 1, 0),
      Resp = function(d) {
        dplyr::case_when(
          d$lR == "left" ~ 0.5,
          d$lR == "right" ~ -0.5,
          d$lR == "time" ~ 0
        )
      }
    ),
    formula = list(
      v ~ 0 + sa:E + match + sa:Time,
      B ~ 0 + sa:E + Resp + Time,
      t0 ~ 0 + (E:sa + Time),
      s ~ 0 + (E + Time),
      A ~ 1,
      sv ~ 1
    ),
    constants = c(
      A = log(0), sv = log(0), lambda_g = log(0), lambda_k = log(0),
      "B_saspeed:E" = log(1), "B_Time" = log(1), t0_Time = log(0.05)
    ),
    pre_transform_terms = list(
      B = c("B_saspeed:E", "B_saneutral:E", "B_saaccuracy:E", "B_Time"),
      v = c("v_saspeed:E", "v_saneutral:E", "v_saaccuracy:E",
            "v_saspeed:Time", "v_saneutral:Time", "v_saaccuracy:Time")
    ),
    bound = list(minmax = cbind(v = c(0, Inf))),
    report_p_vector = FALSE
  )

  p_vec <- sampled_pars(timed_design, doMap = FALSE)
  vals <- c(c(3.99), log(c(2.2, 1.89, 1.59, 1.66, 1.12, 0.86)),
            -0.11, log(c(1.05, 1.08, 0.24, 0.29, 0.29, 1.33, 0.51)))
  names(vals) <- names(p_vec)
  p_vec[names(vals)] <- vals

  emc <- make_emc(dat, timed_design, type = "single", compress = FALSE, n_chains = 1)
  dadm <- emc[[1]]$data[[1]]
  model_obj <- emc[[1]]$model()
  designs <- EMC2:::.oo_expanded_designs(dadm)
  constants <- attr(dadm, "constants")
  if (is.null(constants)) constants <- NA
  p_mat <- matrix(p_vec, nrow = 1, dimnames = list(NULL, names(p_vec)))

  cpp <- EMC2:::get_pars_c_wrapper_oo(
    p_mat, dadm, constants, designs, model_obj$bound, model_obj$transform,
    model_obj$pre_transform, model_obj$trend, FALSE, FALSE, 1L
  )
  r_map <- mapped_pars(timed_design, p_vec, data = dat, digits = 10)
  key <- paste(dadm$sa, dadm$S, dadm$lR, sep = "::")
  expect_equal(cpp[key == "speed::left::left", "v"][1], r_map[key == "speed::left::left", "v"][1], tolerance = 1e-6)
  expect_equal(cpp[key == "speed::left::right", "v"][1], r_map[key == "speed::left::right", "v"][1], tolerance = 1e-6)
  expect_equal(cpp[key == "neutral::left::left", "v"][1], r_map[key == "neutral::left::left", "v"][1], tolerance = 1e-6)
  expect_equal(cpp[key == "accuracy::left::right", "v"][1], r_map[key == "accuracy::left::right", "v"][1], tolerance = 1e-6)
  expect_equal(cpp[key == "speed::left::left", "B"][1], 0.945, tolerance = 1e-10)
  expect_equal(cpp[key == "speed::left::right", "B"][1], 1.055, tolerance = 1e-10)
  expect_lt(max(cpp[, "v"]), 5)
})

test_that("calc_ll_oo does not depend on designs emc2_transform attribute", {
  dat <- data.frame(
    subjects = factor(c(1, 1, 1, 1)),
    S = factor(c("left", "left", "right", "right")),
    E = factor(c("easy", "hard", "easy", "hard"), levels = c("easy", "hard")),
    R = factor(c("left", "left", "right", "right"), levels = c("left", "right")),
    rt = c(0.45, 0.55, 0.5, 0.6)
  )
  ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
  des <- design(
    data = dat,
    model = RDM,
    transform = list(func = c(v = "exp")),
    matchfun = function(d) d$S == d$lR,
    formula = list(v ~ E + lM, B ~ 1, A ~ 1, t0 ~ 1),
    contrasts = list(v = list(lM = ADmat)),
    pre_transform_terms = list(v = c("v", "v_Ehard")),
    bound = list(minmax = cbind(v = c(0, Inf))),
    report_p_vector = FALSE
  )

  emc <- make_emc(dat, des, type = "single", compress = TRUE, n_chains = 1)
  dadm <- emc[[1]]$data[[1]]
  model_obj <- emc[[1]]$model()
  p_vec <- sampled_pars(des, doMap = FALSE)
  p_mat <- matrix(p_vec, nrow = 1, dimnames = list(NULL, names(p_vec)))
  constants <- attr(dadm, "constants")
  if (is.null(constants)) constants <- NA

  designs_keep <- EMC2:::.oo_expanded_designs(dadm)
  ll_keep <- EMC2:::calc_ll_oo(
    p_mat, dadm, constants, designs_keep, model_obj$c_name, model_obj$bound,
    model_obj$transform, model_obj$pre_transform, names(model_obj$p_types), log(1e-10), model_obj$trend
  )

  designs_drop <- designs_keep
  attr(designs_drop, "emc2_transform") <- NULL
  ll_drop <- EMC2:::calc_ll_oo(
    p_mat, dadm, constants, designs_drop, model_obj$c_name, model_obj$bound,
    model_obj$transform, model_obj$pre_transform, names(model_obj$p_types), log(1e-10), model_obj$trend
  )

  expect_equal(ll_keep, ll_drop, tolerance = 1e-12)
})
