# Extracted from test-guess-dispatch.R:172

# prequel ----------------------------------------------------------------------
guess_dispatch_context <- function(model, formula, par_values, rt = 0.1) {
  design_obj <- design(
    factors = list(subjects = 1, S = "stim"),
    Rlevels = c("left", "right"),
    formula = formula,
    model = model,
    report_p_vector = FALSE
  )

  dat <- data.frame(
    subjects = factor(1),
    S = factor("stim"),
    R = factor("left", levels = c("left", "right")),
    rt = rt
  )

  emc <- make_emc(dat, design_obj, type = "single", compress = FALSE, n_chains = 1)
  model_obj <- emc[[1]]$model()
  dadm <- emc[[1]]$data[[1]]

  p_vec <- sampled_pars(design_obj, doMap = FALSE)
  p_vec[names(par_values)] <- unname(par_values)
  p_mat <- matrix(p_vec, nrow = 1, dimnames = list(NULL, names(p_vec)))

  designs <- list()
  for (p in names(model_obj$p_types)) {
    designs[[p]] <- attr(dadm, "designs")[[p]][attr(attr(dadm, "designs")[[p]], "expand"), , drop = FALSE]
  }
  constants <- attr(dadm, "constants")
  if (is.null(constants)) constants <- NA

  list(
    model = model_obj,
    dadm = dadm,
    p_mat = p_mat,
    p_types = names(model_obj$p_types),
    designs = designs,
    constants = constants
  )
}
calc_ll_guess <- function(ctx) {
  EMC2:::calc_ll_oo(
    ctx$p_mat, ctx$dadm,
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
guess_dispatch_cases <- list(
  list(
    label = "RDMGBM",
    formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, lambda_g ~ 1, lambda_k ~ 1),
    pars = c(v = log(1.2), B = log(1.1), A = log(0.2), t0 = log(0.3), s = log(1.0),
             lambda_g = log(0.8), lambda_k = log(0.6)),
    local_guess = EMC2::RDMGBM(erlang_type = "local_guess"),
    local_kill_guess = EMC2::RDMGBM(erlang_type = "local_kill_guess")
  ),
  list(
    label = "RDMSWTN",
    formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1, mG ~ 1, mK ~ 1),
    pars = c(v = 1.2, B = log(1.1), A = log(0.2), t0 = log(0.3), s = log(1.0),
             sv = log(0.25), mG = log(1 / 0.8), mK = log(1 / 0.6)),
    local_guess = EMC2::RDMSWTN(erlang_type = "local_guess"),
    local_kill_guess = EMC2::RDMSWTN(erlang_type = "local_kill_guess")
  ),
  list(
    label = "BAwL",
    formula = list(v ~ 1, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1, lambda_g ~ 1, lambda_k ~ 1),
    pars = c(v = 1.2, sv = log(1.0), B = log(1.1), A = log(0.2), t0 = log(0.3),
             k = log(0.15), lambda_g = log(0.8), lambda_k = log(0.6)),
    local_guess = EMC2::BAwL(erlang_type = "local_guess"),
    local_kill_guess = EMC2::BAwL(erlang_type = "local_kill_guess")
  )
)

# test -------------------------------------------------------------------------
for (case in guess_dispatch_cases[c(1, 2)]) {
    pars_kill_only <- case$pars
    if (case$label == "RDMSWTN") pars_kill_only["mG"] <- log(1e14) else pars_kill_only["lambda_g"] <- log(0)

    ctx_local_kill <- guess_dispatch_context(
      if (case$label == "RDMGBM") EMC2::RDMGBM(erlang_type = "local_kill") else EMC2::RDMSWTN(erlang_type = "local_kill"),
      case$formula,
      pars_kill_only,
      rt = 0.8
    )
    ctx_local_kill_guess <- guess_dispatch_context(
      case$local_kill_guess,
      case$formula,
      pars_kill_only,
      rt = 0.8
    )

    expect_equal(
      calc_ll_guess(ctx_local_kill_guess),
      calc_ll_guess(ctx_local_kill),
      tolerance = 1e-12,
      info = paste(case$label, "aggregate likelihood should keep the kill branch")
    )

    model_fun_kill <- function() ctx_local_kill$model
    model_fun_kill_guess <- function() ctx_local_kill_guess$model

    expect_equal(
      calc_ll_pw(ctx_local_kill_guess$p_mat, ctx_local_kill_guess$dadm, model_fun_kill_guess),
      calc_ll_pw(ctx_local_kill$p_mat, ctx_local_kill$dadm, model_fun_kill),
      tolerance = 1e-12,
      info = paste(case$label, "particlewise likelihood should keep the kill branch")
    )
  }
