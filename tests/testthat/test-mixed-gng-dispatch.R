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

test_that("logical-rule GL path handles unequal t0 within target/nontarget pairs", {
  dat <- data.frame(
    subjects = factor("s1"),
    S = factor("AB", levels = "AB"),
    LogicalRule = factor("OR", levels = "OR"),
    R = factor("yes", levels = c("yes", "no")),
    rt = 0.42
  )

  des <- design(
    data = dat,
    Rlevels = c("yes", "no"),
    fixed_accumulator_roles = factor(c("A", "B", "n_A", "n_B"),
                                     levels = c("A", "B", "n_A", "n_B")),
    matchfun = function(d) d$lR %in% c("A", "B", "n_A", "n_B"),
    model = LogicalRulesLBA,
    formula = list(v ~ 0 + GoA + GoB + NegA + NegB,
                   B ~ 0 + GoA + GoB + NegA + NegB,
                   t0 ~ 0 + GoA + GoB + NegA + NegB,
                   A ~ 1),
    constants = c(sv = log(1)),
    functions = logical_rule_funcs
  )

  p <- sampled_pars(des, doMap = FALSE)
  p[] <- 0
  p[c("v_GoA", "v_GoB", "v_NegA", "v_NegB")] <- c(1.20, 1.10, 1.40, 1.00)
  p[c("B_GoA", "B_GoB", "B_NegA", "B_NegB")] <- log(c(0.80, 0.85, 0.70, 0.90))
  p[c("t0_GoA", "t0_GoB", "t0_NegA", "t0_NegB")] <- log(c(0.25, 0.20, 0.05, 0.20))
  p["A"] <- log(0.20)

  ctx <- build_ll_ctx(dat, des)
  p_mat <- matrix(p, nrow = 1, dimnames = list(NULL, names(p)))
  ll_cpp <- calc_ctx_ll(ctx, p_mat)

  pars <- EMC2:::get_pars_c_wrapper_oo(
    p_mat, ctx$dadm, ctx$constants, ctx$designs,
    ctx$model$bound, ctx$model$transform, ctx$model$pre_transform,
    ctx$model$trend, FALSE, TRUE
  )
  role_par <- function(role) {
    row <- match(role, as.character(ctx$dadm$lR))
    out <- pars[row, c("v", "sv", "B", "A", "t0"), drop = FALSE]
    out <- cbind(out, b = out[, "B"] + out[, "A"])
    out
  }
  d1 <- function(u, p_row) {
    p_eval <- p_row[rep(1, length(u)), , drop = FALSE]
    EMC2:::dLBA(u, p_eval)
  }
  p1 <- function(u, p_row) {
    p_eval <- p_row[rep(1, length(u)), , drop = FALSE]
    EMC2:::pLBA(u, p_eval)
  }
  win_before <- function(t, winner, loser) {
    if (!(t > 0)) return(0)
    stats::integrate(
      function(u) d1(u, winner) * (1 - p1(u, loser)),
      lower = 0, upper = t, rel.tol = 1e-10, abs.tol = 1e-12
    )$value
  }

  t <- ctx$dadm$rt[1]
  A <- role_par("A")
  B <- role_par("B")
  nA <- role_par("n_A")
  nB <- role_par("n_B")

  FA <- p1(t, A); FB <- p1(t, B)
  FnA <- p1(t, nA); FnB <- p1(t, nB)
  gA_yes <- d1(t, A) * (1 - FnA)
  gB_yes <- d1(t, B) * (1 - FnB)
  S_dec_A <- (1 - FA) * (1 - FnA)
  S_dec_B <- (1 - FB) * (1 - FnB)
  GA_no <- win_before(t, nA, A)
  GB_no <- win_before(t, nB, B)
  p_ref <- gA_yes * (GB_no + S_dec_B) + gB_yes * (GA_no + S_dec_A)

  expect_equal(as.numeric(ll_cpp), log(p_ref), tolerance = 2e-4)
})
