test_that("LogicalRulesRDMSWTN uses the shared logical-rules race path", {
  roles <- c("A", "B", "n_A", "n_B")
  lr_functions <- list(
    GoA = function(d) ifelse(d$lR == "A", 1, 0),
    GoB = function(d) ifelse(d$lR == "B", 1, 0),
    NegA = function(d) ifelse(d$lR == "n_A", 1, 0),
    NegB = function(d) ifelse(d$lR == "n_B", 1, 0)
  )
  template <- data.frame(
    subjects = factor("s1"),
    S = factor("AB", levels = "AB"),
    LogicalRule = factor("OR", levels = "OR"),
    R = factor(NA_character_, levels = c("yes", "no"))
  )
  design_template <- template
  design_template$R <- factor(NA_character_, levels = roles)
  des <- design(
    data = design_template,
    Rlevels = c("yes", "no"),
    fixed_accumulator_roles = factor(roles, levels = roles),
    matchfun = function(d) d$lR %in% roles,
    model = LogicalRulesRDMSWTN,
    formula = list(
      v ~ 0 + GoA + GoB + NegA + NegB,
      B ~ 0 + GoA + GoB + NegA + NegB,
      t0 ~ 1, A ~ 1, s ~ 1, sv ~ 1
    ),
    functions = lr_functions
  )

  model <- des$model()
  expect_match(model$c_name, "RDMSWTN")
  expect_match(model$c_name, "LogicalRules")
  expect_equal(names(model$p_types)[seq_len(6)], c("v", "B", "A", "t0", "s", "sv"))

  expect_match(
    LogicalRulesRDMSWTN(erlang_shape = 2L, erlang_type = "local_kill")$c_name,
    "RDMSWTN_E2_LOCAL_KILL_LogicalRules"
  )
  expect_match(LogicalRulesRDMSWTN(posdrift = FALSE)$c_name, "RDMSWTN_IO_LogicalRules")

  p <- sampled_pars(des, doMap = FALSE)
  p[] <- 0
  set_value <- function(name, value) {
    if (name %in% names(p)) p[[name]] <<- value
  }
  for (role in roles) {
    set_value(paste0("v_", role), log(if (role %in% c("A", "B")) 1.6 else 1.1))
    set_value(paste0("B_", role), log(0.8))
  }
  set_value("t0", log(0.2))
  set_value("A", log(0.25))
  set_value("s", log(1))
  set_value("sv", log(0.15))

  set.seed(11)
  simulated <- make_data(p, des, data = template, expand = 40)
  expect_true(all(is.finite(simulated$rt)))
  expect_true(all(as.character(simulated$R) %in% c("yes", "no")))

  emc <- make_emc(simulated, des, type = "single", compress = FALSE, n_chains = 1)
  model <- emc[[1]]$model()
  dadm <- emc[[1]]$data[[1]]
  p_types <- names(model$p_types)
  designs <- setNames(vector("list", length(p_types)), p_types)
  for (param in p_types) {
    dm <- attr(dadm, "designs")[[param]]
    designs[[param]] <- dm[attr(dm, "expand"), , drop = FALSE]
  }
  constants <- attr(dadm, "constants")
  if (is.null(constants)) constants <- NA
  ll <- EMC2:::calc_ll_oo(
    matrix(p, nrow = 1, dimnames = list(NULL, names(p))), dadm,
    constants = constants, designs = designs, type = model$c_name,
    bounds = model$bound, transforms = model$transform,
    pretransforms = model$pre_transform, p_types = p_types,
    min_ll = log(1e-10), trend = model$trend
  )
  expect_true(is.finite(ll))
})
