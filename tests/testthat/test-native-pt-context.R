make_native_pt_context_fixture <- function() {
  des <- design(
    factors = list(subjects = 1, S = 1), Rlevels = 1, model = REXG,
    formula = list(mu ~ 1, sigma ~ 1, tau ~ 1))
  centre <- log(c(mu = .4, sigma = .1, tau = .2))
  set.seed(17)
  dat <- make_data(centre, des, n_trials = 120, rt_resolution = .001)
  emc <- suppressMessages(make_emc(
    dat, des, type = "single", n_chains = 1, rt_resolution = .001,
    compress = FALSE, verbose = FALSE))
  dadm <- emc[[1L]]$data[[1L]]
  model_fn <- emc[[1L]]$model
  model <- model_fn()
  proposals <- matrix(rep(centre, 3L), nrow = 3L, byrow = TRUE,
                      dimnames = list(NULL, names(centre)))
  proposals[2L, 1L] <- proposals[2L, 1L] + .03
  proposals[3L, 2L] <- proposals[3L, 2L] - .02
  constants <- attr(dadm, "constants")
  if (is.null(constants)) constants <- NA_real_
  list(dadm = dadm, model_fn = model_fn, model = model,
       proposals = proposals, constants = constants,
       designs = EMC2:::.oo_expanded_designs(dadm, expand = FALSE),
       p_types = names(model$p_types))
}

calc_context_ll <- function(fx, proposals, context = NULL, constants = fx$constants) {
  EMC2:::calc_ll_oo(
    proposals, fx$dadm, constants, fx$designs, fx$model$c_name,
    fx$model$bound, fx$model$transform, fx$model$pre_transform,
    fx$p_types, log(1e-10), fx$model$trend,
    pt_mapper_context = context)
}

test_that("native parameter contexts preserve likelihoods as particles change", {
  fx <- make_native_pt_context_fixture()
  context <- EMC2:::emc2_create_pt_mapper_context(
    fx$proposals[1L, , drop = FALSE], fx$dadm, fx$constants, fx$designs,
    fx$model$bound, fx$model$transform, fx$model$pre_transform, fx$p_types)

  expect_identical(calc_context_ll(fx, fx$proposals),
                   calc_context_ll(fx, fx$proposals, context))
  changed <- fx$proposals
  changed[, 1L] <- changed[, 1L] + c(.01, -.02, .04)
  expect_identical(calc_context_ll(fx, changed),
                   calc_context_ll(fx, changed, context))

  stale_columns <- fx$proposals
  colnames(stale_columns) <- paste0("stale", seq_len(ncol(stale_columns)))
  expect_error(calc_context_ll(fx, stale_columns, context),
               "stale parameter-mapping context columns")
  shorter_data <- fx$dadm[-nrow(fx$dadm), , drop = FALSE]
  expect_error(EMC2:::calc_ll_oo(
    fx$proposals, shorter_data, fx$constants, fx$designs, fx$model$c_name,
    fx$model$bound, fx$model$transform, fx$model$pre_transform,
    fx$p_types, log(1e-10), fx$model$trend,
    pt_mapper_context = context),
    "stale parameter-mapping context dimensions")
})

test_that("manager retains a context only for validated internal calls", {
  fx <- make_native_pt_context_fixture()
  state <- EMC2:::.emc_native_pt_context_state
  saved <- list(pid = state$pid, key = state$key, ptr = state$ptr)
  old_option <- getOption("emc2.native_pt_context")
  on.exit({
    state$pid <- saved$pid
    state$key <- saved$key
    state$ptr <- saved$ptr
    options(emc2.native_pt_context = old_option)
  }, add = TRUE)

  state$pid <- NA_integer_
  state$key <- NULL
  state$ptr <- NULL
  options(emc2.native_pt_context = TRUE)
  public <- EMC2:::calc_ll_manager(fx$proposals, fx$dadm, fx$model_fn)
  expect_null(state$ptr)

  internal <- EMC2:::calc_ll_manager(
    fx$proposals, fx$dadm, fx$model_fn, ll_cache_validated = TRUE)
  expect_identical(internal, public)
  expect_true(inherits(state$ptr, "EMC2PtMapperContext"))
  first_context <- state$ptr
  expect_identical(EMC2:::.emc_native_pt_context_get(
    fx$proposals, fx$dadm, fx$model, fx$constants, fx$designs,
    fx$p_types), first_context)

  changed_model <- fx$model
  changed_model$c_name <- paste0(changed_model$c_name, "_variant")
  rekeyed <- EMC2:::.emc_native_pt_context_get(
    fx$proposals, fx$dadm, changed_model, fx$constants, fx$designs,
    fx$p_types)
  expect_false(identical(rekeyed, first_context))

  constants_one <- c(context_probe = 1)
  constants_two <- c(context_probe = 2)
  first_constant_context <- EMC2:::.emc_native_pt_context_get(
    fx$proposals, fx$dadm, fx$model, constants_one, fx$designs, fx$p_types)
  changed_constants <- EMC2:::.emc_native_pt_context_get(
    fx$proposals, fx$dadm, fx$model, constants_two, fx$designs, fx$p_types)
  expect_false(identical(first_constant_context, changed_constants))

  state$ptr <- NULL
  state$key <- NULL
  options(emc2.native_pt_context = FALSE)
  expect_null(EMC2:::.emc_native_pt_context_get(
    fx$proposals, fx$dadm, fx$model, fx$constants, fx$designs,
    fx$p_types))
  expect_null(state$ptr)
})
