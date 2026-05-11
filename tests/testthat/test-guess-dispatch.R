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

test_that("guess-model calc_ll_oo dispatch handles rt<t0 correctly", {
  min_ll <- log(1e-10)

  for (case in guess_dispatch_cases) {
    ctx_guess <- guess_dispatch_context(case$local_guess, case$formula, case$pars)
    pars_no_guess <- case$pars
    if (case$label == "RDMSWTN") pars_no_guess["mG"] <- log(1e14) else pars_no_guess["lambda_g"] <- log(0)
    ctx_guess_off <- guess_dispatch_context(case$local_guess, case$formula, pars_no_guess)
    ctx_kill_guess <- guess_dispatch_context(case$local_kill_guess, case$formula, case$pars)
    ctx_kill_guess_off <- guess_dispatch_context(case$local_kill_guess, case$formula, pars_no_guess)

    ll_guess <- calc_ll_guess(ctx_guess)
    ll_guess_off <- calc_ll_guess(ctx_guess_off)
    ll_kill_guess <- calc_ll_guess(ctx_kill_guess)
    ll_kill_guess_off <- calc_ll_guess(ctx_kill_guess_off)

    expect_equal(
      ll_guess_off,
      min_ll,
      tolerance = 1e-12,
      info = paste(case$label, "local_guess should collapse to min_ll when lambda_g is zero")
    )
    expect_true(
      is.finite(ll_guess) && ll_guess > ll_guess_off,
      info = paste(case$label, "local_guess should retain pre-t0 guess likelihood")
    )
    expect_equal(
      ll_kill_guess_off,
      min_ll,
      tolerance = 1e-12,
      info = paste(case$label, "local_kill_guess should collapse to min_ll when lambda_g is zero")
    )
    expect_true(
      is.finite(ll_kill_guess) && ll_kill_guess > ll_kill_guess_off,
      info = paste(case$label, "local_kill_guess should retain pre-t0 guess likelihood")
    )
  }
})

test_that("guess-model calc_ll_oo_pw matches calc_ll_oo", {
  for (case in guess_dispatch_cases) {
    for (variant in c("local_guess", "local_kill_guess")) {
      ctx <- guess_dispatch_context(case[[variant]], case$formula, case$pars)
      ll <- calc_ll_guess(ctx)
      model_fun <- function() ctx$model
      ll_pw <- calc_ll_pw(ctx$p_mat, ctx$dadm, model_fun)

      expect_equal(nrow(ll_pw), 1L, info = paste(case$label, variant))
      expect_equal(
        sum(ll_pw[1, ]),
        ll,
        tolerance = 1e-12,
        info = paste(case$label, variant, "particlewise likelihood should sum to aggregate likelihood")
      )
      expect_true(
        all(is.finite(ll_pw[1, ])) && all(ll_pw[1, ] > log(1e-10)),
        info = paste(case$label, variant, "particlewise likelihood should stay finite for pre-t0 guess trials")
      )
    }
  }
})

test_that("mixed RDMSWTN guess dispatch carries omega through particle likelihood", {
  formula <- list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1, mG ~ 1, mK ~ 1, omega ~ 1)
  pars <- c(v = 1.2, B = log(1.1), A = log(0.2), t0 = log(0.3), s = log(1.0),
            sv = log(0.25), mG = log(1 / 0.8), mK = log(1 / 0.6), omega = qnorm(0.35))
  model <- EMC2::RDMSWTN(erlang_shape = "mixed", erlang_type = "local_kill_guess")
  ctx <- guess_dispatch_context(model, formula, pars)

  ll <- calc_ll_guess(ctx)
  model_fun <- function() ctx$model
  ll_pw <- calc_ll_pw(ctx$p_mat, ctx$dadm, model_fun)

  expect_match(ctx$model$c_name, "RDMSWTN_EMIX")
  expect_true(is.finite(ll) && ll > log(1e-10))
  expect_equal(sum(ll_pw[1, ]), ll, tolerance = 1e-12)
})

test_that("local_kill_guess falls back to local_kill when lambda_g is zero after t0", {
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
})
