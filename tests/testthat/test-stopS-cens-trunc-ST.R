test_that("SSexG truncation normalisation supports stop-triggered accumulators (R vs C++)", {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(44)

  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

  lIfun_ST <- function(d) {
    # Mark one accumulator as stop-triggered (ST == 1); the others are GO (== 2).
    is_st <- as.character(d$lR) %in% c("STOP", "stop")
    factor(ifelse(is_st, 1, 2), levels = 1:2)
  }

  mySSD_function <- function(d) SSD_function(d, SSD = 0.20, pSSD = 0.6)

  design_ss <- design(
    model = SSexG,
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right", "STOP"),
    matchfun = matchfun,
    functions = list(lI = lIfun_ST, SSD = mySSD_function),
    formula = list(mu ~ 0 + lR, sigma ~ 1, tau ~ 1, muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
  )

  p_vector <- sampled_pars(design_ss, doMap = FALSE)
  mu_names <- names(p_vector)[grepl("^mu_", names(p_vector))]
  stop_mu <- mu_names[grepl("STOP|stop", mu_names)]
  p_vector[mu_names] <- log(0.60)
  p_vector[stop_mu] <- log(0.30)
  p_vector[["sigma"]] <- log(0.10)
  p_vector[["tau"]] <- log(0.15)
  p_vector[["muS"]] <- log(0.35)
  p_vector[["sigmaS"]] <- log(0.06)
  p_vector[["tauS"]] <- log(0.10)
  p_vector[["gf"]] <- qnorm(0.15)
  p_vector[["tf"]] <- qnorm(0.20)

  dat <- make_data(p_vector, design_ss, n_trials = 250, LT = 0.05, UT = 1.20, UC = Inf)
  expect_gt(nrow(dat), 40)

  dadm <- EMC2:::design_model(dat, design_ss, verbose = FALSE)
  expect_true("LT" %in% names(dadm))
  expect_true("UT" %in% names(dadm))
  expect_gt(sum(dadm$lI == levels(dadm$lI)[1]), 0)
  expect_true(any(dadm$winner & dadm$lI == levels(dadm$lI)[1] & is.finite(dadm$SSD)))

  model <- attr(dadm, "model")()
  pars <- EMC2:::get_pars_matrix(p_vector, dadm, model)
  ll_r <- model$log_likelihood_cens_trunc(pars, dadm, model, min_ll = log(1e-10), normalise_trunc = TRUE)

  p_matrix <- matrix(p_vector, nrow = 1)
  colnames(p_matrix) <- names(p_vector)
  p_types <- names(model$p_types)

  designs <- list()
  for (p in p_types) {
    designs[[p]] <- attr(dadm, "designs")[[p]][attr(attr(dadm, "designs")[[p]], "expand"), , drop = FALSE]
  }
  constants <- attr(dadm, "constants")
  if (is.null(constants)) constants <- NA
  pre_transform <- model$pre_transform
  if (is.null(pre_transform)) pre_transform <- list()
  trend <- model$trend
  if (is.null(trend)) trend <- list()

  ll_c <- EMC2:::calc_ll(
    p_matrix, dadm, constants, designs, "SSexG",
    model$bound, model$transform, pre_transform,
    p_types, log(1e-10), trend
  )

  expect_equal(as.numeric(ll_c), as.numeric(ll_r), tolerance = 1e-6)
})

