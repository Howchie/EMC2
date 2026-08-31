omx_test_values <- function(model) {
  m <- if (is.function(model)) model() else model
  p <- m$p_types
  out <- p
  for (nm in names(out)) {
    out[[nm]] <- switch(nm,
      mu = log(2), sigma = log(.4),
      v = .8, sv = log(.3),
      B = log(1), A = log(.2), t0 = log(.1),
      k = log(.5), ell = log(.1), kappa = log(.5), p = log(1),
      lambda = qnorm(.5), tau = log(1), tau_s = log(1), tau_t = log(1),
      pi = qnorm(.5), mG = log(1), mK = log(1),
      shape = log(1), mean = log(1), delta = 0,
      eta = .2, s = log(1), alpha = qnorm(.7), beta = log(1),
      h = qnorm(.95),
      pContaminant = qnorm(.01), pGuess = qnorm(0),
      out[[nm]])
  }
  out
}

omx_test_emc <- function(model, tc = list()) {
  data <- data.frame(
    R = factor(rep(c("a", "b"), each = 4)),
    rt = rep(.4, 8), subjects = factor(rep(1, 8))
  )
  m <- if (is.function(model)) model() else model
  formula <- lapply(names(m$p_types), function(nm)
    stats::as.formula(paste0(nm, " ~ 1")))
  des <- design(data = data, model = model, formula = formula, TC = tc)
  emc <- make_emc(data, des, type = "single", compress = FALSE, n_chains = 1)
  pnames <- names(sampled_pars(des))
  values <- omx_test_values(model)[pnames]
  alpha <- array(as.numeric(values), dim = c(length(values), 1L, 1L),
                 dimnames = list(pnames, "1", NULL))
  emc[[1]]$samples$alpha <- alpha
  emc[[1]]$samples$stage <- "sample"
  emc[[1]]$samples$idx <- 1L
  emc[[1]]$init <- TRUE
  emc[[1]]$par_names <- pnames
  emc
}

test_that("decomposition retains newer BAwD options and finite UT", {
  model <- BAwD(drift_distribution = "lognormal", gamma = .5, rho = 2)
  emc <- omx_test_emc(model, list(UC = .8, UT = 1.5, UCresponse = FALSE))
  out <- decompose_omissions(emc, stat = "mean")

  expect_match(out$c_name, "GAM12_RHO2")
  expect_false(any(out$summary$component == "unexplained"))
  expect_equal(out$draws$contaminant, .01, tolerance = 1e-8)
  expect_true(all(abs(out$draws$omission_total + out$draws$responded - 1) < 1e-8))
  expect_true(any(out$draws$censor_slow > 0))
  expect_true("asymptotic_subthreshold" %in% out$summary$component)
  audit <- omission_mechanisms(emc)
  expect_true("censor_slow" %in% audit$mechanism)
  expect_false("live" %in% names(audit))
})

test_that("BAwD splits dead launches from the turnaround residual", {
  emc <- omx_test_emc(function() BAwD("normal", gamma = .5, rho = 2),
                       list(UC = Inf, UT = Inf))
  out <- decompose_omissions(emc, stat = "mean")
  expect_true(all(c("dead_launch", "asymptotic_subthreshold") %in%
                  as.character(out$summary$component)))
  expect_gt(out$draws$dead_launch, 0)
  expect_gt(out$draws$asymptotic_subthreshold, 0)

  prep <- EMC2:::.omx_prep(emc)
  mech <- EMC2:::.omx_mechanisms(prep$model_list, prep$pars_ref)
  split <- EMC2:::.omx_split_inf(prep$model_list, mech, prep$pars_ref,
                                 nrow(prep$trial), prep$n_acc)
  expect_equal(split$dead_launch + split$asymptotic_subthreshold,
               split$residual, tolerance = 1e-12)
})

test_that("IO clock and BAwL asymptotic mass are separated", {
  model <- BAwL(posdrift = FALSE, erlang_type = "local_kill")
  emc <- omx_test_emc(model)
  audit <- omission_mechanisms(emc)
  expect_true(all(c("kill", "negative_drift",
                    "asymptotic_subthreshold") %in% audit$mechanism))
  expect_false("leak" %in% audit$mechanism)
  expect_match(audit$description[
    audit$mechanism == "asymptotic_subthreshold"
  ], "V < k b", fixed = TRUE)
  out <- decompose_omissions(emc, stat = "mean")
  expect_true(all(c("kill", "negative_drift", "asymptotic_subthreshold") %in%
                    as.character(out$summary$component)))
  expect_false("leak" %in% as.character(out$summary$component))
  expect_gt(out$draws$asymptotic_subthreshold, 0)
  expect_true(all(out$draws$omission_total >= 0 & out$draws$omission_total <= 1))
  printed <- capture.output(print(out))
  expect_false(any(grepl("censor_slow", printed, fixed = TRUE)))
  expect_false(any(grepl("live", printed, fixed = TRUE)))

  lba <- decompose_omissions(
    omx_test_emc(function() LBA(posdrift = FALSE)), stat = "mean"
  )
  expect_true("negative_drift" %in% as.character(lba$summary$component))
})

test_that("supported deterministic race constructors are accepted", {
  models <- list(
    BAwD, BAwDp, BAwF, BAwR, BTAwL,
    BTAwLTransient, BTAwLSustained, LBA, FRQ
  )
  for (model in models) {
    emc <- omx_test_emc(model)
    out <- decompose_omissions(emc, stat = "mean")
    expect_s3_class(out, "omission_decomposition")
    expect_true(all(out$draws$omission_total >= -1e-10 &
                    out$draws$omission_total <= 1 + 1e-10),
                info = model()$c_name)
  }
})

test_that("shared kill clocks and finite-reservoir defects are decomposed", {
  emc_kill <- omx_test_emc(function()
    BAwL(erlang_type = "global_kill"))
  kill <- decompose_omissions(emc_kill, stat = "mean")
  expect_true("kill" %in% as.character(kill$summary$component))
  expect_true(all(abs(kill$draws$omission_total + kill$draws$responded - 1) < 1e-8))

  emc_frq <- omx_test_emc(FRQ)
  frq <- decompose_omissions(emc_frq, stat = "mean")
  expect_true(any(frq$draws$asymptotic_subthreshold > 0))
  expect_true(all(abs(frq$draws$omission_total + frq$draws$responded - 1) < 1e-8))
})
