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
  # R contains overt rule responses in the simulation data, but the design
  # scaffold must expose the latent accumulator roles used by the formulas.
  design_template <- template
  design_template$R <- factor(NA_character_, levels = c("A", "B", "n_A", "n_B"))
  design(
    data = design_template,
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

# OR_DETECTION_GNG is the four-horse OR task (A, n_A, B, n_B) whose "no"
# outcome is a withheld response (rt = Inf). Drifts are driven by the stimulus
# S (present targets fast, absent targets slow), so — unlike the old
# detector-gating model — a false alarm ("yes") remains possible on NN trials.
gng_stim_levels <- c("NN", "AN", "NB", "AB")

gng_design <- function() {
  template <- data.frame(
    subjects = factor(rep("s1", 4)),
    S = factor(gng_stim_levels, levels = gng_stim_levels),
    LogicalRule = factor(rep("OR_DETECTION_GNG", 4), levels = "OR_DETECTION_GNG"),
    R = factor(rep(NA_character_, 4), levels = c("A", "B", "n_A", "n_B"))
  )
  design(
    data = template,
    Rlevels = c("yes", "no"),
    fixed_accumulator_roles = factor(c("A", "B", "n_A", "n_B"),
                                     levels = c("A", "B", "n_A", "n_B")),
    matchfun = function(d) d$lR %in% c("A", "B", "n_A", "n_B"),
    model = LogicalRulesLBA,
    formula = list(v ~ 0 + GoA:S + GoB:S + NegA:S + NegB:S,
                   B ~ 0 + GoA + GoB + NegA + NegB,
                   t0 ~ 1, A ~ 1),
    constants = c(sv = log(1)),
    functions = lr_funcs
  )
}

gng_pvec <- function(des) {
  p <- sampled_pars(des, doMap = FALSE)
  p[] <- 0
  setv <- function(nm, val) if (nm %in% names(p)) p[[nm]] <<- val
  for (S in gng_stim_levels) {
    Apres <- S %in% c("AN", "AB"); Bpres <- S %in% c("NB", "AB")
    setv(paste0("v_GoA:S", S),  if (Apres) 1.8 else 0.5)
    setv(paste0("v_GoB:S", S),  if (Bpres) 1.5 else 0.4)
    setv(paste0("v_NegA:S", S), if (Apres) 0.5 else 1.6)
    setv(paste0("v_NegB:S", S), if (Bpres) 0.6 else 1.5)
  }
  setv("B_GoA", log(0.75)); setv("B_GoB", log(0.95))
  setv("B_NegA", log(0.85)); setv("B_NegB", log(0.70))
  setv("t0", log(0.2)); setv("A", log(0.3))
  p
}

test_that("OR_DETECTION_GNG censor masses match four-horse simulation", {
  skip_on_cran()
  n_sim <- 5e4
  LC <- 0.5
  UC <- 1.2

  des <- gng_design()
  p <- gng_pvec(des)

  template <- data.frame(
    subjects = factor(rep("s1", 4)),
    S = factor(gng_stim_levels, levels = gng_stim_levels),
    LogicalRule = factor(rep("OR_DETECTION_GNG", 4), levels = "OR_DETECTION_GNG"),
    R = factor(rep(NA_character_, 4), levels = c("yes", "no"))
  )
  set.seed(103)
  sim <- make_data(p, des, data = template, expand = n_sim)

  gng_case <- function(S, R, rt, LC = NA, UC = NA) {
    d <- data.frame(
      subjects = factor("s1"),
      S = factor(S, levels = gng_stim_levels),
      LogicalRule = factor("OR_DETECTION_GNG", levels = "OR_DETECTION_GNG"),
      R = factor(R, levels = c("yes", "no")),
      rt = rt
    )
    if (!is.na(LC)) d$LC <- LC
    if (!is.na(UC)) d$UC <- UC
    d
  }

  # Lower censor: an overt go response ("yes") fired in [0, LC]. This is now a
  # four-horse subrace outcome, possible even on NN (a false alarm). The mass is
  # response-type agnostic (single overt type), so R recorded or missing match.
  for (S in gng_stim_levels) {
    p_mc <- mean(sim$S == S & !is.na(sim$R) & sim$rt < LC, na.rm = TRUE) /
      mean(sim$S == S)
    for (r in c("yes", NA_character_)) {
      ctx <- build_ll_ctx(gng_case(S, r, -Inf, LC = LC), des)
      p_ll <- exp(calc_ctx_ll(ctx, p))
      expect_lt(abs(p_ll - p_mc), mc_tol(p_mc, n_sim / 4),
                label = sprintf("GNG lower-censor S=%s R=%s |%g - %g|", S, r, p_ll, p_mc))
    }
  }

  # NN false alarms are now possible: the lower-censor go mass is strictly
  # positive (the old detector-gating model floored this to min_ll).
  ll_nn_low <- calc_ctx_ll(build_ll_ctx(gng_case("NN", "yes", -Inf, LC = LC), des), p)
  expect_gt(ll_nn_low, log(1e-6))

  # Upper censor (rt = +Inf): no overt go response by UC = both subraces
  # resolved "no" by UC, or neither channel had said "yes" yet.
  for (S in gng_stim_levels) {
    denom <- mean(sim$S == S)
    p_mc <- mean(sim$S == S & (is.na(sim$R) | sim$rt > UC), na.rm = TRUE) / denom
    ctx <- build_ll_ctx(gng_case(S, NA_character_, Inf, UC = UC), des)
    p_ll <- exp(calc_ctx_ll(ctx, p))
    expect_lt(abs(p_ll - p_mc), mc_tol(p_mc, n_sim / 4),
              label = sprintf("GNG upper-censor S=%s |%g - %g|", S, p_ll, p_mc))
  }
})
