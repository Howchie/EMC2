# Regression checks for the cross-subrace LogicalRules capacity route.

lr_capacity_functions <- list(
  GoA = function(d) ifelse(d$lR == "A", 1, 0),
  GoB = function(d) ifelse(d$lR == "B", 1, 0),
  NegA = function(d) ifelse(d$lR == "n_A", 1, 0),
  NegB = function(d) ifelse(d$lR == "n_B", 1, 0)
)

lr_capacity_data <- function() {
  data.frame(
    subjects = factor(rep("s1", 4)),
    S = factor(c("AN", "NB", "AB", "NN"),
               levels = c("NN", "AN", "NB", "AB")),
    LogicalRule = factor(rep("OR", 4), levels = "OR"),
    R = factor(rep("yes", 4), levels = c("yes", "no")),
    rt = rep(0.6, 4)
  )
}

lr_capacity_design <- function(data, capacity) {
  design_data <- data
  design_data$R <- factor(NA_character_,
                          levels = c("A", "B", "n_A", "n_B"))
  design(
    data = design_data,
    Rlevels = c("yes", "no"),
    fixed_accumulator_roles = factor(c("A", "B", "n_A", "n_B"),
                                     levels = c("A", "B", "n_A", "n_B")),
    matchfun = function(d) TRUE,
    model = if (capacity) {
      function() LogicalRulesLBA(capacity = TRUE)
    } else {
      LogicalRulesLBA
    },
    formula = if (capacity) {
      list(v ~ 0 + GoA + GoB + NegA + NegB,
           B ~ 0 + GoA + GoB + NegA + NegB,
           t0 ~ 1, A ~ 1, kappa ~ 1, tau ~ 1)
    } else {
      list(v ~ 0 + GoA + GoB + NegA + NegB,
           B ~ 0 + GoA + GoB + NegA + NegB,
           t0 ~ 1, A ~ 1)
    },
    constants = c(sv = log(1)),
    functions = lr_capacity_functions
  )
}

lr_capacity_parameters <- function(design, tau = 0) {
  p <- sampled_pars(design, doMap = FALSE)
  p[] <- 0
  values <- c(
    v_GoA = 1.6, v_GoB = 1.1, v_NegA = 0.95, v_NegB = 1.25,
    B_GoA = log(0.75), B_GoB = log(0.95),
    B_NegA = log(0.85), B_NegB = log(0.70),
    t0 = log(0.2), A = log(0.3)
  )
  for (nm in names(values)) {
    if (nm %in% names(p)) p[nm] <- values[[nm]]
  }
  if (all(c("kappa", "tau") %in% names(p))) {
    p["kappa"] <- 0                 # log(1)
    p["tau"] <- if (tau == 0) log(0) else log(tau)
  }
  p
}

lr_capacity_detection_design <- function(rule, data, roles) {
  design_data <- data
  design_data$R <- factor(NA_character_, levels = roles)
  design(
    data = design_data,
    Rlevels = c("yes", "no"),
    fixed_accumulator_roles = factor(roles, levels = roles),
    matchfun = function(d) TRUE,
    model = function() LogicalRulesLBA(capacity = TRUE),
    formula = list(v ~ 0 + GoA + GoB, B ~ 0 + GoA + GoB,
                   t0 ~ 1, A ~ 1, kappa ~ 1, tau ~ 1),
    constants = c(sv = log(1)),
    functions = lr_capacity_functions
  )
}

lr_capacity_ll_pw <- function(data, des, p) {
  emc <- make_emc(data, des, type = "single", compress = FALSE, n_chains = 1)
  model <- emc[[1]]$model()
  dadm <- emc[[1]]$data[[1]]
  p_types <- names(model$p_types)
  designs <- setNames(lapply(p_types, function(nm) {
    dm <- attr(dadm, "designs")[[nm]]
    dm[attr(dm, "expand"), , drop = FALSE]
  }), p_types)
  do.call(EMC2:::calc_ll_oo_pw, list(
    matrix(p, nrow = 1, dimnames = list(NULL, names(p))), dadm,
    constants = attr(dadm, "constants"), designs = designs,
    type = model$c_name, bounds = model$bound,
    transforms = model$transform, pretransforms = model$pre_transform,
    p_types = p_types, min_ll = log(1e-10), trend = model$trend
  ))[1, ]
}

test_that("capacity baseline is exactly the ordinary LogicalRules route", {
  data <- lr_capacity_data()
  ordinary <- lr_capacity_design(data, FALSE)
  capacity <- lr_capacity_design(data, TRUE)
  p0 <- lr_capacity_parameters(ordinary)
  p1 <- lr_capacity_parameters(capacity, tau = 0)
  p1[names(p0)] <- p0

  ll0 <- lr_capacity_ll_pw(data, ordinary, p0)
  ll1 <- lr_capacity_ll_pw(data, capacity, p1)
  expect_equal(ll1, ll0, tolerance = 1e-12)
})

test_that("only the AB trial enters the shared-capacity route", {
  data <- lr_capacity_data()
  capacity <- lr_capacity_design(data, TRUE)
  p0 <- lr_capacity_parameters(capacity, tau = 0)
  p1 <- lr_capacity_parameters(capacity, tau = 0.4)

  ll0 <- lr_capacity_ll_pw(data, capacity, p0)
  ll1 <- lr_capacity_ll_pw(data, capacity, p1)
  expect_equal(ll1[c(1, 2, 4)], ll0[c(1, 2, 4)], tolerance = 1e-12)
  expect_gt(abs(ll1[3] - ll0[3]), 1e-7)
})

test_that("capacity simulator accepts one shared factor per AB trial", {
  data <- lr_capacity_data()
  capacity <- lr_capacity_design(data, TRUE)
  p <- lr_capacity_parameters(capacity, tau = 0.4)
  withr::local_options(emc2.cpp_rfun = TRUE)
  set.seed(1)
  simulated <- make_data(p, capacity, data = data, expand = 100)
  expect_equal(nrow(simulated), 400)
  expect_true(all(is.finite(simulated$rt)))
})

test_that("capacity simulator loads only A and B on AB trials", {
  finish_cpp <- getFromNamespace("logicalrules_capacity_finish_cpp", "EMC2")
  pars <- matrix(
    rep(c(1, 0, 1, 0, 0, 2, 0), 8),
    ncol = 7, byrow = TRUE,
    dimnames = list(NULL, c("v", "sv", "b", "A", "t0", "kappa", "tau"))
  )
  finish <- finish_cpp(
    pars,
    c("A", "B", "n_A", "n_B"),
    c("AB", "AN"),
    posdrift = FALSE
  )

  expect_equal(unname(finish[1, c("A", "B")]), c(0.5, 0.5), tolerance = 1e-12)
  expect_equal(unname(finish[1, c("n_A", "n_B")]), c(1, 1), tolerance = 1e-12)
  expect_equal(unname(finish[2, ]), rep(1, 4), tolerance = 1e-12)
})

test_that("capacity detection routes evaluate analytic and GNG AB trials", {
  data <- lr_capacity_data()[3, , drop = FALSE]
  data$LogicalRule <- factor("OR_DETECTION_ANALYTIC",
                             levels = "OR_DETECTION_ANALYTIC")
  analytic <- lr_capacity_detection_design(
    "OR_DETECTION_ANALYTIC", data, c("A", "B"))

  # GNG is the full four-horse OR task (a withheld response in place of the
  # overt "no"), so it uses the ordinary A/n_A/B/n_B capacity design.
  gng_data <- data
  gng_data$LogicalRule <- factor("OR_DETECTION_GNG", levels = "OR_DETECTION_GNG")
  gng <- lr_capacity_design(gng_data, TRUE)

  expect_true(is.finite(lr_capacity_ll_pw(data, analytic,
                                         lr_capacity_parameters(analytic, 0.4))))
  # A finite-RT "yes" (go) and a withheld (rt = Inf) AB trial both evaluate.
  expect_true(is.finite(lr_capacity_ll_pw(gng_data, gng,
                                         lr_capacity_parameters(gng, 0.4))))
  gng_withheld <- gng_data
  gng_withheld$R <- factor(NA_character_, levels = c("yes", "no"))
  gng_withheld$rt <- Inf
  expect_true(is.finite(lr_capacity_ll_pw(gng_withheld, gng,
                                         lr_capacity_parameters(gng, 0.4))))
})

test_that("detection ignores an invalid inactive detector row", {
  data <- data.frame(
    subjects = factor(c("s1", "s1")),
    S = factor(c("A", "B"), levels = c("A", "B", "AB")),
    LogicalRule = factor(rep("OR_DETECTION_ANALYTIC", 2),
                         levels = "OR_DETECTION_ANALYTIC"),
    R = factor(rep("yes", 2), levels = c("yes", "no")),
    rt = c(0.6, 0.6)
  )
  analytic <- lr_capacity_detection_design(
    "OR_DETECTION_ANALYTIC", data, c("A", "B"))
  p <- lr_capacity_parameters(analytic, tau = 0)
  p["v_GoB"] <- NA_real_

  ll <- lr_capacity_ll_pw(data, analytic, p)
  expect_true(is.finite(ll[1]))
  expect_equal(ll[2], log(1e-10))
})

test_that("capacity censoring rejects impossible positive quadrature masses", {
  data <- lr_capacity_data()[3, , drop = FALSE]
  data$LogicalRule <- factor("OR_DETECTION_ANALYTIC",
                             levels = "OR_DETECTION_ANALYTIC")
  data$rt <- Inf
  data$R <- factor(NA_character_, levels = c("yes", "no"))
  data$LT <- 0.1
  data$UC <- 2.5
  data$UT <- Inf
  analytic <- lr_capacity_detection_design(
    "OR_DETECTION_ANALYTIC", data, c("A", "B"))
  p <- lr_capacity_parameters(analytic, tau = 0.4)

  # These are deliberately extreme sampler-scale values representative of
  # the runaway alpha draws.  The censored event and truncation masses remain
  # probabilities, so their conditional log likelihood cannot be positive.
  p[c("v_GoA", "v_GoB")] <- c(426.73, 101.42)
  p[c("B_GoA", "B_GoB")] <- 327.61
  p["A"] <- 342.65
  p["t0"] <- -2.74
  p["kappa"] <- 292.43
  p["tau"] <- 347.41

  withr::local_envvar(c(
    EMC2_LRCAP_SCAN_N = "12",
    EMC2_LRCAP_FINE_N = "12"
  ))
  ll <- lr_capacity_ll_pw(data, analytic, p)
  expect_true(all(is.finite(ll)))
  expect_true(all(ll <= 1e-8))
})

test_that("capacity route counters distinguish ordinary and AB trials", {
  skip_if_not(is.loaded("_EMC2_lr_capacity_counter_values", PACKAGE = "EMC2"))
  data <- lr_capacity_data()
  capacity <- lr_capacity_design(data, TRUE)
  p <- lr_capacity_parameters(capacity, tau = 0.4)
  withr::local_envvar(EMC2_LRCAP_COUNTERS = "1")
  .Call("_EMC2_lr_capacity_counters_reset", PACKAGE = "EMC2")
  invisible(lr_capacity_ll_pw(data, capacity, p))
  counters <- .Call("_EMC2_lr_capacity_counter_values", PACKAGE = "EMC2")
  expect_equal(counters$ordinary_trials, 3)
  expect_equal(counters$capacity_choice_trials, 1)
  expect_gt(counters$factor_node_evaluations, 0)
})

test_that("truncation normalisers are reused across RTs", {
  skip_if_not(is.loaded("_EMC2_lr_capacity_counter_values", PACKAGE = "EMC2"))
  data <- lr_capacity_data()
  data <- rbind(data[1:3, , drop = FALSE], data[3, , drop = FALSE],
                data[4, , drop = FALSE])
  data$rt <- c(0.5, 0.6, 0.7, 0.8, 0.9)
  data$LT <- rep(0.1, nrow(data))
  data$UT <- rep(1.2, nrow(data))
  capacity <- lr_capacity_design(data, TRUE)
  p <- lr_capacity_parameters(capacity, tau = 0.4)

  withr::local_envvar(c(
    EMC2_LRCAP_COUNTERS = "1",
    EMC2_LRCAP_SCAN_N = "12",
    EMC2_LRCAP_FINE_N = "12"
  ))
  invisible(lr_capacity_ll_pw(data, capacity, p))
  counters <- .Call("_EMC2_lr_capacity_counter_values", PACKAGE = "EMC2")

  # Two AB trials have different RTs but the same parameter cell and window:
  # two numerator passes plus one shared denominator pass.
  expect_equal(counters$capacity_choice_trials, 2)
  expect_equal(counters$factor_node_evaluations, 3 * (12 + 12))
})
