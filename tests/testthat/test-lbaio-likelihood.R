test_that("LBAIO likelihood matches between R and C++ for finite RTs", {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(42)

  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

  design_lbaio <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = LBAIO,
    formula = list(B ~ 1, v ~ 0 + lM, A ~ 1, t0 ~ 1),
    constants = c(sv = log(1))
  )

  p_vector <- sampled_pars(design_lbaio, doMap = FALSE)
  p_vector[["B"]] <- log(1.0)
  p_vector[["A"]] <- log(0.4)
  p_vector[["t0"]] <- log(0.25)
  p_vector[["v_lMFALSE"]] <- -0.2
  p_vector[["v_lMTRUE"]] <- 1.0

  dat <- make_data(
    p_vector, design_lbaio, n_trials = 200,
    UC = Inf,
    UCresponse = FALSE,
    rtContaminantNA = FALSE
  )

  dat <- dat[is.finite(dat$rt) & !is.na(dat$R), , drop = FALSE]
  expect_gt(nrow(dat), 20)

  dadm <- EMC2:::design_model(dat, design_lbaio, verbose = FALSE)

  ll_r <- EMC2:::calc_ll_R(p_vector, LBAIO(), dadm)

  p_matrix <- matrix(p_vector, nrow = 1)
  colnames(p_matrix) <- names(p_vector)
  model <- attr(dadm, "model")()
  p_types <- names(model$p_types)

  designs <- list()
  for (p in p_types) {
    designs[[p]] <- attr(dadm, "designs")[[p]][attr(attr(dadm, "designs")[[p]], "expand"), , drop = FALSE]
  }
  constants <- attr(dadm, "constants")
  if (is.null(constants)) constants <- NA

  ll_c <- EMC2:::calc_ll(
    p_matrix, dadm, constants, designs, model$c_name,
    model$bound, model$transform, model$pre_transform,
    p_types, log(1e-10), model$trend
  )

  expect_equal(as.numeric(ll_c), as.numeric(ll_r), tolerance = 1e-6)
})

test_that("LBAIO likelihood matches between R and C++ with deadline censoring", {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(43)

  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

  design_lbaio <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = LBAIO,
    formula = list(B ~ 1, v ~ 0 + lM, A ~ 1, t0 ~ 1),
    constants = c(sv = log(1))
  )

  p_vector <- sampled_pars(design_lbaio, doMap = FALSE)
  p_vector[["B"]] <- log(1.0)
  p_vector[["A"]] <- log(0.4)
  p_vector[["t0"]] <- log(0.25)
  p_vector[["v_lMFALSE"]] <- -0.2
  p_vector[["v_lMTRUE"]] <- 1.0

  dat <- make_data(
    p_vector, design_lbaio, n_trials = 200,
    UC = 0.6,
    UCresponse = FALSE,
    rtContaminantNA = FALSE
  )

  expect_true(any(is.infinite(dat$rt)))

  dadm <- EMC2:::design_model(dat, design_lbaio, verbose = FALSE)

  ll_r <- EMC2:::calc_ll_R(p_vector, LBAIO(), dadm)

  p_matrix <- matrix(p_vector, nrow = 1)
  colnames(p_matrix) <- names(p_vector)
  model <- attr(dadm, "model")()
  p_types <- names(model$p_types)

  designs <- list()
  for (p in p_types) {
    designs[[p]] <- attr(dadm, "designs")[[p]][attr(attr(dadm, "designs")[[p]], "expand"), , drop = FALSE]
  }
  constants <- attr(dadm, "constants")
  if (is.null(constants)) constants <- NA

  ll_c <- EMC2:::calc_ll(
    p_matrix, dadm, constants, designs, model$c_name,
    model$bound, model$transform, model$pre_transform,
    p_types, log(1e-10), model$trend
  )

  expect_equal(as.numeric(ll_c), as.numeric(ll_r), tolerance = 1e-6)
})
