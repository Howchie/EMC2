# Monte-Carlo validation of the logical-rules censoring branches in
# c_log_likelihood_logicalrules (particle_ll.cpp).
#
# A lower-censored trial (rt == -Inf) codes an overt rule response in
# [LT, LC]; an upper-censored trial (rt == +Inf) codes no overt response by
# UC. The C++ masses are compared against empirical outcome proportions from
# LogicalRules_rfun simulation (via make_data), which defines the reference
# outcome semantics (apply_logical_rules / the detection branch).

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

calc_ctx_ll <- function(ctx, p_vec) {
  EMC2:::calc_ll_oo(
    matrix(p_vec, nrow = 1, dimnames = list(NULL, names(p_vec))),
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

lr_funcs <- list(
  GoA = function(d) ifelse(d$lR == "A", 1, 0),
  GoB = function(d) ifelse(d$lR == "B", 1, 0),
  NegA = function(d) ifelse(d$lR == "n_A", 1, 0),
  NegB = function(d) ifelse(d$lR == "n_B", 1, 0),
  NoGo = function(d) ifelse(d$lR == "nogo", 1, 0)
)

# Unequal channels (A vs B) so the general (non-channels_equal) path runs.
lr_pvec <- function(design) {
  p <- sampled_pars(design, doMap = FALSE)
  p[] <- 0
  vals <- c(v_GoA = 1.60, v_GoB = 1.10, v_NegA = 0.95, v_NegB = 1.25,
            B_GoA = log(0.75), B_GoB = log(0.95), B_NegA = log(0.85), B_NegB = log(0.70),
            v_NoGo = 1.30, B_NoGo = log(0.80),
            t0 = log(0.2), A = log(0.3))
  for (nm in names(vals)) if (nm %in% names(p)) p[[nm]] <- vals[[nm]]
  p
}

make_rule_design <- function(template, Rlevels) {
  design(
    data = template,
    Rlevels = Rlevels,
    fixed_accumulator_roles = factor(c("A", "B", "n_A", "n_B"),
                                     levels = c("A", "B", "n_A", "n_B")),
    matchfun = function(d) d$lR %in% c("A", "B", "n_A", "n_B"),
    model = LogicalRulesLBA,
    formula = list(v ~ 0 + GoA + GoB + NegA + NegB,
                   B ~ 0 + GoA + GoB + NegA + NegB,
                   t0 ~ 1, A ~ 1),
    constants = c(sv = log(1)),
    functions = lr_funcs
  )
}

censor_case <- function(template, rule, Rlevels, R, rt, LC = NA, UC = NA) {
  d <- data.frame(
    subjects = factor("s1"),
    S = factor("AB", levels = "AB"),
    LogicalRule = factor(rule, levels = rule),
    R = factor(R, levels = Rlevels),
    rt = rt
  )
  if (!is.na(LC)) d$LC <- LC
  if (!is.na(UC)) d$UC <- UC
  d
}

mc_tol <- function(p_hat, n) 4 * sqrt(pmax(p_hat, 1e-4) * (1 - pmax(p_hat, 1e-4)) / n) + 2e-4

test_that("OR/AND/XOR lower- and upper-censor masses match simulation", {
  skip_on_cran()
  n_sim <- 1e5
  LC <- 0.55
  UC <- 1.4

  for (rule in c("OR", "AND", "XOR")) {
    template <- data.frame(
      subjects = factor("s1"),
      S = factor("AB", levels = "AB"),
      LogicalRule = factor(rule, levels = rule),
      R = factor(NA_character_, levels = c("yes", "no"))
    )
    des <- make_rule_design(template, Rlevels = c("yes", "no"))
    p <- lr_pvec(des)

    set.seed(101)
    sim <- make_data(p, des, data = template, expand = n_sim)

    for (r in c("yes", "no")) {
      p_mc <- mean(sim$R == r & sim$rt < LC, na.rm = TRUE)
      case <- censor_case(template, rule, c("yes", "no"), R = r, rt = -Inf, LC = LC)
      ctx <- build_ll_ctx(case, des)
      p_ll <- exp(calc_ctx_ll(ctx, p))
      expect_lt(abs(p_ll - p_mc), mc_tol(p_mc, n_sim),
                label = sprintf("%s lower-censor R=%s |%g - %g|", rule, r, p_ll, p_mc))
    }

    # Upper censor: no overt response by UC (resp-agnostic mass; R must still
    # be recorded for rules <= 4, but does not enter the omission probability).
    p_mc_uc <- mean(sim$rt > UC, na.rm = TRUE)
    case_uc <- censor_case(template, rule, c("yes", "no"), R = "yes", rt = Inf, UC = UC)
    ctx_uc <- build_ll_ctx(case_uc, des)
    p_ll_uc <- exp(calc_ctx_ll(ctx_uc, p))
    expect_lt(abs(p_ll_uc - p_mc_uc), mc_tol(p_mc_uc, n_sim),
              label = sprintf("%s upper-censor |%g - %g|", rule, p_ll_uc, p_mc_uc))
  }
})

test_that("ID lower-censor masses match simulation for every response identity", {
  skip_on_cran()
  n_sim <- 1e5
  LC <- 0.75
  id_levels <- c("NN", "AN", "NB", "AB")

  template <- data.frame(
    subjects = factor("s1"),
    S = factor("AB", levels = "AB"),
    LogicalRule = factor("ID", levels = "ID"),
    R = factor(NA_character_, levels = id_levels)
  )
  des <- make_rule_design(template, Rlevels = id_levels)
  p <- lr_pvec(des)

  set.seed(102)
  sim <- make_data(p, des, data = template, expand = n_sim)

  for (r in id_levels) {
    p_mc <- mean(sim$R == r & sim$rt < LC, na.rm = TRUE)
    case <- censor_case(template, "ID", id_levels, R = r, rt = -Inf, LC = LC)
    ctx <- build_ll_ctx(case, des)
    p_ll <- exp(calc_ctx_ll(ctx, p))
    expect_lt(abs(p_ll - p_mc), mc_tol(p_mc, n_sim),
              label = sprintf("ID lower-censor R=%s |%g - %g|", r, p_ll, p_mc))
  }
})

test_that("OR_DETECTION_GNG censor masses exclude nogo wins and match simulation", {
  skip_on_cran()
  n_sim <- 5e4
  LC <- 0.5
  UC <- 1.2
  stim_levels <- c("NN", "AN", "NB", "AB")

  template <- data.frame(
    subjects = factor(rep("s1", 4)),
    S = factor(stim_levels, levels = stim_levels),
    LogicalRule = factor(rep("OR_DETECTION_GNG", 4), levels = "OR_DETECTION_GNG"),
    R = factor(rep(NA_character_, 4), levels = c("yes", "no"))
  )
  des <- design(
    data = template,
    Rlevels = c("yes", "no"),
    fixed_accumulator_roles = factor(c("A", "B", "nogo"), levels = c("A", "B", "nogo")),
    matchfun = function(d) d$lR %in% c("A", "B", "nogo"),
    model = LogicalRulesLBA,
    formula = list(v ~ 0 + GoA + GoB + NoGo,
                   B ~ 0 + GoA + GoB + NoGo,
                   t0 ~ 1, A ~ 1),
    constants = c(sv = log(1)),
    functions = lr_funcs
  )
  p <- lr_pvec(des)

  set.seed(103)
  sim <- make_data(p, des, data = template, expand = n_sim)

  gng_case <- function(S, R, rt, LC = NA, UC = NA) {
    d <- data.frame(
      subjects = factor("s1"),
      S = factor(S, levels = stim_levels),
      LogicalRule = factor("OR_DETECTION_GNG", levels = "OR_DETECTION_GNG"),
      R = factor(R, levels = c("yes", "no")),
      rt = rt
    )
    if (!is.na(LC)) d$LC <- LC
    if (!is.na(UC)) d$UC <- UC
    d
  }

  # Lower censor on go-stimulus conditions: mass = P(go detector beats nogo
  # AND rt < LC), with R recorded ("yes") or missing (same single overt type).
  for (S in c("AN", "NB", "AB")) {
    p_mc <- mean(sim$S == S & !is.na(sim$R) & sim$rt < LC, na.rm = TRUE) /
      mean(sim$S == S)
    for (r in c("yes", NA_character_)) {
      ctx <- build_ll_ctx(gng_case(S, r, -Inf, LC = LC), des)
      p_ll <- exp(calc_ctx_ll(ctx, p))
      expect_lt(abs(p_ll - p_mc), mc_tol(p_mc, n_sim / 4),
                label = sprintf("GNG lower-censor S=%s R=%s |%g - %g|", S, r, p_ll, p_mc))
    }
  }

  # NN condition: no go detector races, so an overt response before LC is
  # impossible — a nogo finish is a withheld response, not an observation.
  ctx_nn <- build_ll_ctx(gng_case("NN", NA_character_, -Inf, LC = LC), des)
  expect_equal(calc_ctx_ll(ctx_nn, p), log(1e-10))

  # Upper censor (regression check on the already-correct branch): no overt
  # response by UC = nogo won, or nothing finished by UC.
  for (S in c("NN", "AN", "AB")) {
    denom <- mean(sim$S == S)
    p_mc <- mean(sim$S == S & (is.na(sim$R) | sim$rt > UC), na.rm = TRUE) / denom
    ctx <- build_ll_ctx(gng_case(S, NA_character_, Inf, UC = UC), des)
    p_ll <- exp(calc_ctx_ll(ctx, p))
    expect_lt(abs(p_ll - p_mc), mc_tol(p_mc, n_sim / 4),
              label = sprintf("GNG upper-censor S=%s |%g - %g|", S, p_ll, p_mc))
  }
})
