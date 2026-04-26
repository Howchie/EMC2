build_ll_ctx <- function(data, design) {
  emc <- make_emc(data, design, type = "single", compress = FALSE, n_chains = 1)
  model <- emc[[1]]$model()
  dadm <- emc[[1]]$data[[1]]
  p_types <- names(model$p_types)
  designs <- setNames(vector("list", length(p_types)), p_types)
  for (p in p_types) {
    dm <- attr(dadm, "designs")[[p]]
    designs[[p]] <- dm[attr(dm, "expand"), , drop = FALSE]
  }
  constants <- attr(dadm, "constants")
  if (is.null(constants)) constants <- NA
  list(dadm = dadm, model = model, p_types = p_types, designs = designs, constants = constants)
}

calc_ctx_ll <- function(ctx, p_mat) {
  EMC2:::calc_ll_oo(
    p_mat,
    ctx$dadm,
    constants = ctx$constants,
    designs = ctx$designs,
    type = ctx$model$c_name,
    bounds = ctx$model$bound,
    transforms = ctx$model$transform,
    pretransforms = ctx$model$pre_transform,
    p_types = ctx$p_types,
    min_ll = log(1e-10),
    trend = ctx$model$trend
  )
}

set_detection_params <- function(p_vec, values) {
  for (nm in names(values)) {
    if (nm %in% names(p_vec)) p_vec[[nm]] <- values[[nm]]
  }
  p_vec
}

detection_funcs <- list(
  GoA = function(d) ifelse(d$lR == "A", 1, 0),
  GoB = function(d) ifelse(d$lR == "B", 1, 0),
  NoGo = function(d) ifelse(d$lR == "nogo", 1, 0)
)

logical_rule_funcs <- c(
  detection_funcs,
  list(
    NegA = function(d) ifelse(d$lR == "n_A", 1, 0),
    NegB = function(d) ifelse(d$lR == "n_B", 1, 0)
  )
)

test_that("mixed 2AFC and GNG pooled race data dispatches by active nogo accumulator", {
  pooled <- data.frame(
    subjects = factor(rep("s1", 6)),
    Condition = factor(c("2afc", "2afc", "2afc", "gng", "gng", "gng"),
                       levels = c("2afc", "gng")),
    RACE = factor(c("2", "2", "2", "3", "3", "3"), levels = c("2", "3")),
    R = factor(c("left", NA, NA, "right", NA, NA), levels = c("left", "right", "nogo")),
    rt = c(0.35, Inf, -Inf, 0.42, Inf, -Inf),
    LC = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2),
    UC = c(0.8, 0.8, 0.8, 0.8, 0.8, 0.8)
  )

  pooled_design <- design(
    data = pooled,
    model = LBA,
    formula = list(v ~ 1, B ~ 1, t0 ~ 1, A ~ 1)
  )
  pooled_ctx <- build_ll_ctx(pooled, pooled_design)

  afc_data <- droplevels(subset(pooled, Condition == "2afc", select = c(subjects, R, rt, LC, UC)))
  afc_data$R <- factor(afc_data$R, levels = c("left", "right"))
  afc_design <- design(
    data = afc_data,
    model = LBA,
    formula = list(v ~ 1, B ~ 1, t0 ~ 1, A ~ 1)
  )
  afc_ctx <- build_ll_ctx(afc_data, afc_design)

  gng_data <- droplevels(subset(pooled, Condition == "gng", select = c(subjects, R, rt, LC, UC)))
  gng_data$R <- factor(gng_data$R, levels = c("left", "right", "nogo"))
  gng_design <- design(
    data = gng_data,
    model = LBA,
    formula = list(v ~ 1, B ~ 1, t0 ~ 1, A ~ 1)
  )
  gng_ctx <- build_ll_ctx(gng_data, gng_design)

  expect_identical(pooled_ctx$p_types, afc_ctx$p_types)
  expect_identical(pooled_ctx$p_types, gng_ctx$p_types)

  set.seed(42)
  p_mat <- matrix(rnorm(16 * length(pooled_ctx$p_types), sd = 0.2), nrow = 16)
  colnames(p_mat) <- pooled_ctx$p_types

  pooled_ll <- calc_ctx_ll(pooled_ctx, p_mat)
  split_ll <- calc_ctx_ll(afc_ctx, p_mat) + calc_ctx_ll(gng_ctx, p_mat)

  expect_true(all(is.finite(pooled_ll)))
  expect_true(all(is.finite(split_ll)))
  expect_equal(pooled_ll, split_ll, tolerance = 1e-8)
})

test_that("make_missing only disables truncation on rows with active nogo", {
  dat <- data.frame(
    subjects = factor(c("s1", "s1")),
    RACE = factor(c("2", "3"), levels = c("2", "3")),
    R = factor(c(NA, NA), levels = c("left", "right", "nogo")),
    rt = c(0.5, 0.5)
  )

  out <- make_missing(dat, LT = 0.2, UT = 0.9, rt_resolution = NULL)

  expect_equal(out$LT[1], 0.2)
  expect_equal(out$UT[1], 0.9)
  expect_equal(out$LT[2], 0)
  expect_true(is.infinite(out$UT[2]))
})

test_that("make_data supports pooled mixed 2AFC and GNG race templates", {
  template <- data.frame(
    subjects = factor(rep("s1", 4)),
    Condition = factor(c("2afc", "2afc", "gng", "gng"), levels = c("2afc", "gng")),
    S = factor(c("left", "right", "left", "stop"), levels = c("left", "right", "stop")),
    RACE = factor(c("2", "2", "3", "3"), levels = c("2", "3")),
    R = factor(rep(NA_character_, 4), levels = c("left", "right", "nogo"))
  )

  matchfun <- function(d) {
    (d$Condition == "2afc" & as.character(d$S) == as.character(d$lR)) |
      (d$Condition == "gng" &
         ((d$S == "left" & d$lR == "left") | (d$S == "stop" & d$lR == "nogo")))
  }

  design_mixed <- design(
    data = template,
    Rlevels = c("left", "right", "nogo"),
    matchfun = matchfun,
    model = LBA,
    formula = list(v ~ 0 + lM, B ~ 0 + lR, t0 ~ 1, A ~ 1),
    constants = c(sv = log(1))
  )

  p_vector <- sampled_pars(design_mixed, doMap = FALSE)
  if ("v_lMFALSE" %in% names(p_vector)) p_vector[["v_lMFALSE"]] <- 0.3
  if ("v_lMTRUE" %in% names(p_vector)) p_vector[["v_lMTRUE"]] <- 1.6
  if ("B_lRleft" %in% names(p_vector)) p_vector[["B_lRleft"]] <- log(0.8)
  if ("B_lRright" %in% names(p_vector)) p_vector[["B_lRright"]] <- log(0.8)
  if ("B_lRnogo" %in% names(p_vector)) p_vector[["B_lRnogo"]] <- log(0.6)
  if ("t0" %in% names(p_vector)) p_vector[["t0"]] <- log(0.2)
  if ("A" %in% names(p_vector)) p_vector[["A"]] <- log(0.3)

  set.seed(7)
  dat <- make_data(
    p_vector,
    design_mixed,
    data = template,
    expand = 40,
    TC = list(LT = 0.2, UT = 0.9, UC = 0.8)
  )

  expect_true("RACE" %in% names(dat))
  expect_true(all(dat$LT[dat$Condition == "gng"] == 0))
  expect_true(all(is.infinite(dat$UT[dat$Condition == "gng"])))
  expect_true(all(dat$LT[dat$Condition == "2afc"] == 0.2))
  expect_true(all(dat$UT[dat$Condition == "2afc"] == 0.9))
})

test_that("OR_DETECTION_ANALYTIC matches RedundantTarget analytic likelihood", {
  rt_template <- data.frame(
    subjects = factor(rep("s1", 3)),
    S = factor(c("A", "B", "AB"), levels = c("A", "B", "AB")),
    R = factor(rep(NA_character_, 3), levels = c("yes"))
  )
  rt_matchfun <- function(d) {
    (d$S == "A" & d$lR == "A") | (d$S == "B" & d$lR == "B") | (d$S == "AB")
  }
  rt_design <- design(
    data = rt_template,
    Rlevels = c("yes"),
    fixed_accumulator_roles = factor(c("A", "B"), levels = c("A", "B")),
    matchfun = rt_matchfun,
    model = RedundantTargetLBA,
    formula = list(v ~ 0 + GoA + GoB, B ~ 0 + GoA + GoB, t0 ~ 1, A ~ 1),
    constants = c(sv = log(1), pContaminant = qnorm(0)),
    functions = detection_funcs
  )
  p_rt <- sampled_pars(rt_design, doMap = FALSE)
  p_rt <- set_detection_params(
    p_rt,
    c(v_GoA = 1.4, v_GoB = 1.1, B_GoA = log(0.8), B_GoB = log(0.75), t0 = log(0.2), A = log(0.3))
  )
  dat_rt <- make_data(p_rt, rt_design, data = rt_template, expand = 80)
  rt_ctx <- build_ll_ctx(dat_rt, rt_design)

  lr_data <- within(dat_rt, {
    LogicalRule <- factor("OR_DETECTION_ANALYTIC", levels = c("OR_DETECTION_ANALYTIC", "OR"))
  })
  lr_design <- design(
    data = lr_data,
    Rlevels = c("yes"),
    fixed_accumulator_roles = factor(c("A", "B", "n_A", "n_B", "nogo"),
                                     levels = c("A", "B", "n_A", "n_B", "nogo")),
    matchfun = function(d) d$lR %in% c("A", "B"),
    model = LogicalRulesLBA,
    formula = list(v ~ 0 + GoA + GoB + NoGo, B ~ 0 + GoA + GoB + NoGo, t0 ~ 1, A ~ 1),
    constants = c(sv = log(1)),
    functions = detection_funcs
  )
  p_lr <- sampled_pars(lr_design, doMap = FALSE)
  p_lr <- set_detection_params(
    p_lr,
    c(v_GoA = 1.4, v_GoB = 1.1, B_GoA = log(0.8), B_GoB = log(0.75),
      v_NoGo = 0.2, B_NoGo = log(0.9), t0 = log(0.2), A = log(0.3))
  )
  lr_ctx <- build_ll_ctx(lr_data, lr_design)

  ll_rt <- calc_ctx_ll(rt_ctx, matrix(p_rt, nrow = 1, dimnames = list(NULL, names(p_rt))))
  ll_lr <- calc_ctx_ll(lr_ctx, matrix(p_lr, nrow = 1, dimnames = list(NULL, names(p_lr))))
  expect_equal(ll_lr, ll_rt, tolerance = 1e-8)
})

test_that("OR_DETECTION_GNG matches RedundantTarget nogo likelihood", {
  rt_template <- data.frame(
    subjects = factor(rep("s1", 3)),
    S = factor(c("A", "B", "AB"), levels = c("A", "B", "AB")),
    R = factor(rep(NA_character_, 3), levels = c("yes", "nogo"))
  )
  rt_matchfun <- function(d) {
    (d$S == "A" & d$lR == "A") |
      (d$S == "B" & d$lR == "B") |
      (d$S == "AB" & d$lR %in% c("A", "B")) |
      d$lR == "nogo"
  }
  rt_design <- design(
    data = rt_template,
    Rlevels = c("yes", "nogo"),
    fixed_accumulator_roles = factor(c("A", "B", "nogo"), levels = c("A", "B", "nogo")),
    matchfun = rt_matchfun,
    model = RedundantTargetLBA,
    formula = list(v ~ 0 + GoA + GoB + NoGo, B ~ 0 + GoA + GoB + NoGo, t0 ~ 1, A ~ 1),
    constants = c(sv = log(1), pContaminant = qnorm(0)),
    UC = 0.8,
    functions = detection_funcs
  )
  p_rt <- sampled_pars(rt_design, doMap = FALSE)
  p_rt <- set_detection_params(
    p_rt,
    c(v_GoA = 1.5, v_GoB = 1.2, v_NoGo = 1.8,
      B_GoA = log(0.8), B_GoB = log(0.75), B_NoGo = log(0.65),
      t0 = log(0.2), A = log(0.3))
  )
  dat_rt <- make_data(p_rt, rt_design, data = rt_template, expand = 80, TC = list(UC = 0.8))
  rt_ctx <- build_ll_ctx(dat_rt, rt_design)

  lr_data <- within(dat_rt, {
    LogicalRule <- factor("OR_DETECTION_GNG", levels = c("OR_DETECTION_GNG", "OR"))
  })
  lr_design <- design(
    data = lr_data,
    Rlevels = c("yes", "nogo"),
    fixed_accumulator_roles = factor(c("A", "B", "n_A", "n_B", "nogo"),
                                     levels = c("A", "B", "n_A", "n_B", "nogo")),
    matchfun = function(d) d$lR %in% c("A", "B", "nogo"),
    model = LogicalRulesLBA,
    formula = list(v ~ 0 + GoA + GoB + NoGo, B ~ 0 + GoA + GoB + NoGo, t0 ~ 1, A ~ 1),
    constants = c(sv = log(1)),
    UC = 0.8,
    functions = detection_funcs
  )
  p_lr <- sampled_pars(lr_design, doMap = FALSE)
  p_lr <- set_detection_params(
    p_lr,
    c(v_GoA = 1.5, v_GoB = 1.2, v_NoGo = 1.8,
      B_GoA = log(0.8), B_GoB = log(0.75), B_NoGo = log(0.65),
      t0 = log(0.2), A = log(0.3))
  )
  lr_ctx <- build_ll_ctx(lr_data, lr_design)

  ll_rt <- calc_ctx_ll(rt_ctx, matrix(p_rt, nrow = 1, dimnames = list(NULL, names(p_rt))))
  ll_lr <- calc_ctx_ll(lr_ctx, matrix(p_lr, nrow = 1, dimnames = list(NULL, names(p_lr))))
  expect_equal(ll_lr, ll_rt, tolerance = 1e-8)
})

test_that("pooled analytic and GNG detection rules recover rowwise likelihoods", {
  pooled_template <- data.frame(
    subjects = factor(rep("s1", 6)),
    S = factor(c("A", "B", "AB", "A", "B", "AB"), levels = c("A", "B", "AB")),
    LogicalRule = factor(
      c(rep("OR_DETECTION_ANALYTIC", 3), rep("OR_DETECTION_GNG", 3)),
      levels = c("OR_DETECTION_ANALYTIC", "OR_DETECTION_GNG", "OR")
    ),
    R = factor(rep(NA_character_, 6), levels = c("yes", "nogo"))
  )

  pooled_matchfun <- function(d) {
    is_analytic <- d$LogicalRule == "OR_DETECTION_ANALYTIC"
    is_gng <- d$LogicalRule == "OR_DETECTION_GNG"
    (is_analytic & d$lR %in% c("A", "B")) |
      (is_gng & d$lR %in% c("A", "B", "nogo"))
  }

  pooled_design <- design(
    data = pooled_template,
    Rlevels = c("yes", "nogo"),
    fixed_accumulator_roles = factor(c("A", "B", "n_A", "n_B", "nogo"),
                                     levels = c("A", "B", "n_A", "n_B", "nogo")),
    matchfun = pooled_matchfun,
    model = LogicalRulesLBA,
    formula = list(v ~ 0 + GoA + GoB + NoGo, B ~ 0 + GoA + GoB + NoGo, t0 ~ 1, A ~ 1),
    constants = c(sv = log(1)),
    UC = 0.8,
    functions = detection_funcs
  )
  p_pool <- sampled_pars(pooled_design, doMap = FALSE)
  p_pool <- set_detection_params(
    p_pool,
    c(v_GoA = 1.35, v_GoB = 1.1, v_NoGo = 2.2,
      B_GoA = log(0.8), B_GoB = log(0.75), B_NoGo = log(0.55),
      t0 = log(0.2), A = log(0.3))
  )

  set.seed(11)
  pooled_dat <- make_data(p_pool, pooled_design, data = pooled_template, expand = 200, TC = list(UC = 0.8))
  pooled_ctx <- build_ll_ctx(pooled_dat, pooled_design)

  analytic_dat <- droplevels(subset(pooled_dat, LogicalRule == "OR_DETECTION_ANALYTIC"))
  analytic_dat$R <- factor(analytic_dat$R, levels = c("yes"))
  analytic_design <- design(
    data = analytic_dat,
    Rlevels = c("yes"),
    fixed_accumulator_roles = factor(c("A", "B"), levels = c("A", "B")),
    matchfun = function(d) {
      (d$S == "A" & d$lR == "A") | (d$S == "B" & d$lR == "B") | (d$S == "AB")
    },
    model = RedundantTargetLBA,
    formula = list(v ~ 0 + GoA + GoB, B ~ 0 + GoA + GoB, t0 ~ 1, A ~ 1),
    constants = c(sv = log(1), pContaminant = qnorm(0)),
    functions = detection_funcs
  )
  p_analytic <- sampled_pars(analytic_design, doMap = FALSE)
  p_analytic <- set_detection_params(
    p_analytic,
    c(v_GoA = 1.35, v_GoB = 1.1, B_GoA = log(0.8), B_GoB = log(0.75), t0 = log(0.2), A = log(0.3))
  )
  analytic_ctx <- build_ll_ctx(analytic_dat, analytic_design)

  gng_dat <- droplevels(subset(pooled_dat, LogicalRule == "OR_DETECTION_GNG"))
  gng_dat$R <- factor(gng_dat$R, levels = c("yes", "nogo"))
  gng_design <- design(
    data = gng_dat,
    Rlevels = c("yes", "nogo"),
    fixed_accumulator_roles = factor(c("A", "B", "nogo"), levels = c("A", "B", "nogo")),
    matchfun = function(d) {
      (d$S == "A" & d$lR == "A") |
        (d$S == "B" & d$lR == "B") |
        (d$S == "AB" & d$lR %in% c("A", "B")) |
        d$lR == "nogo"
    },
    model = RedundantTargetLBA,
    formula = list(v ~ 0 + GoA + GoB + NoGo, B ~ 0 + GoA + GoB + NoGo, t0 ~ 1, A ~ 1),
    constants = c(sv = log(1), pContaminant = qnorm(0)),
    UC = 0.8,
    functions = detection_funcs
  )
  p_gng <- sampled_pars(gng_design, doMap = FALSE)
  p_gng <- set_detection_params(
    p_gng,
    c(v_GoA = 1.35, v_GoB = 1.1, v_NoGo = 2.2,
      B_GoA = log(0.8), B_GoB = log(0.75), B_NoGo = log(0.55),
      t0 = log(0.2), A = log(0.3))
  )
  gng_ctx <- build_ll_ctx(gng_dat, gng_design)

  expect_true(any(is.na(gng_dat$R)))

  ll_pool <- calc_ctx_ll(pooled_ctx, matrix(p_pool, nrow = 1, dimnames = list(NULL, names(p_pool))))
  ll_split <- calc_ctx_ll(analytic_ctx, matrix(p_analytic, nrow = 1, dimnames = list(NULL, names(p_analytic)))) +
    calc_ctx_ll(gng_ctx, matrix(p_gng, nrow = 1, dimnames = list(NULL, names(p_gng))))

  expect_equal(ll_pool, ll_split, tolerance = 1e-8)
})

test_that("full four-accumulator OR rules are unchanged by adding a dormant nogo role", {
  or_template <- data.frame(
    subjects = factor(rep("s1", 4)),
    S = factor(c("NN", "AN", "NB", "AB"), levels = c("NN", "AN", "NB", "AB")),
    LogicalRule = factor(rep("OR", 4), levels = c("OR")),
    R = factor(rep(NA_character_, 4), levels = c("yes", "no"))
  )

  or_matchfun <- function(d) dplyr::case_when(
    d$S == "NN" & d$lR == "n_A" ~ TRUE,
    d$S == "NN" & d$lR == "n_B" ~ TRUE,
    d$S == "AN" & d$lR == "A" ~ TRUE,
    d$S == "AN" & d$lR == "n_B" ~ TRUE,
    d$S == "NB" & d$lR == "n_A" ~ TRUE,
    d$S == "NB" & d$lR == "B" ~ TRUE,
    d$S == "AB" & d$lR == "A" ~ TRUE,
    d$S == "AB" & d$lR == "B" ~ TRUE,
    TRUE ~ FALSE
  )

  base_design <- design(
    data = or_template,
    Rlevels = c("yes", "no"),
    fixed_accumulator_roles = factor(c("A", "B", "n_A", "n_B"), levels = c("A", "B", "n_A", "n_B")),
    matchfun = or_matchfun,
    model = LogicalRulesLBA,
    formula = list(v ~ 0 + GoA + GoB + NegA + NegB,
                   B ~ 0 + GoA + GoB + NegA + NegB,
                   t0 ~ 1, A ~ 1),
    constants = c(sv = log(1)),
    functions = logical_rule_funcs
  )
  p_base <- sampled_pars(base_design, doMap = FALSE)
  p_base <- set_detection_params(
    p_base,
    c(v_GoA = 1.45, v_GoB = 1.15, v_NegA = 0.75, v_NegB = 0.7,
      B_GoA = log(0.8), B_GoB = log(0.78), B_NegA = log(0.85), B_NegB = log(0.88),
      t0 = log(0.2), A = log(0.3))
  )

  set.seed(21)
  or_dat <- make_data(p_base, base_design, data = or_template, expand = 200)
  base_ctx <- build_ll_ctx(or_dat, base_design)

  pooled_design <- design(
    data = transform(or_dat, R = factor(R, levels = c("yes", "no", "nogo"))),
    Rlevels = c("yes", "no", "nogo"),
    fixed_accumulator_roles = factor(c("A", "B", "n_A", "n_B", "nogo"),
                                     levels = c("A", "B", "n_A", "n_B", "nogo")),
    matchfun = or_matchfun,
    model = LogicalRulesLBA,
    formula = list(v ~ 0 + GoA + GoB + NegA + NegB + NoGo,
                   B ~ 0 + GoA + GoB + NegA + NegB + NoGo,
                   t0 ~ 1, A ~ 1),
    constants = c(sv = log(1)),
    functions = logical_rule_funcs
  )
  p_pool <- sampled_pars(pooled_design, doMap = FALSE)
  p_pool <- set_detection_params(
    p_pool,
    c(v_GoA = 1.45, v_GoB = 1.15, v_NegA = 0.75, v_NegB = 0.7, v_NoGo = 2.1,
      B_GoA = log(0.8), B_GoB = log(0.78), B_NegA = log(0.85), B_NegB = log(0.88), B_NoGo = log(0.58),
      t0 = log(0.2), A = log(0.3))
  )
  pooled_ctx <- build_ll_ctx(transform(or_dat, R = factor(R, levels = c("yes", "no", "nogo"))), pooled_design)

  ll_base <- calc_ctx_ll(base_ctx, matrix(p_base, nrow = 1, dimnames = list(NULL, names(p_base))))
  ll_pool <- calc_ctx_ll(pooled_ctx, matrix(p_pool, nrow = 1, dimnames = list(NULL, names(p_pool))))

  expect_equal(ll_pool, ll_base, tolerance = 1e-8)
})
