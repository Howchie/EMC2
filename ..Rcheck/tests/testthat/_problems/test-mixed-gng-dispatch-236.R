# Extracted from test-mixed-gng-dispatch.R:236

# prequel ----------------------------------------------------------------------
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

# test -------------------------------------------------------------------------
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
