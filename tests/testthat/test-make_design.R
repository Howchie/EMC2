design_data <- design(data = data.frame(forstmann, CO = 1:nrow(forstmann)),
            model=LBA,matchfun=function(d)d$S==d$lR,
            formula=list(v~lM,sv~lM,B~E+lR,t0~E2 + CO),
            contrasts=list(v = list(lM=matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d")))),
            constants=c(sv=log(1)),
            functions = list(
              E2 = function(d) factor(d$E!="speed",labels=c("speed","nonspeed"))
            ))


design_custom <- design(factors = list(S = c("left", "right"),
                                         subjects = 1:3),
                          Rlevels = c("left", "right"), model = LNR,
                          formula =list(m~0+S,s~1, t0~1),
                          constants=c(s=log(1)))

test_that("design", {
  expect_snapshot(
    str(design_data, give.attr = FALSE)
  )
  expect_snapshot(
    str(design_custom, give.attr = FALSE)
  )
})

test_that("auto covariate detection ignores unused numeric columns", {
  dat <- data.frame(forstmann, CO = seq_len(nrow(forstmann)),
                    UNUSED_NUM = seq_len(nrow(forstmann)))
  des <- design(
    data = dat,
    model = LBA,
    matchfun = function(d) d$S == d$lR,
    formula = list(v ~ lM, sv ~ lM, B ~ E + lR, t0 ~ E2 + CO),
    contrasts = list(v = list(lM = matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d")))),
    constants = c(sv = log(1)),
    functions = list(
      E2 = function(d) factor(d$E != "speed", labels = c("speed", "nonspeed"))
    )
  )

  expect_true("CO" %in% des$Fcovariates)
  expect_false("UNUSED_NUM" %in% des$Fcovariates)
})

test_that("pre_transform_terms applies selected race-model contrasts on the natural scale", {
  dat <- data.frame(
    subjects = factor(c(1, 1, 1, 1)),
    S = factor(c("left", "left", "right", "right")),
    E = factor(c("easy", "hard", "easy", "hard"), levels = c("easy", "hard")),
    R = factor(c("left", "left", "right", "right"), levels = c("left", "right")),
    rt = c(0.5, 0.6, 0.55, 0.65)
  )
  ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
  matchfun <- function(d) d$S == d$lR

  design_default <- design(
    data = dat,
    model = RDM,
    matchfun = matchfun,
    formula = list(v ~ E + lM, B ~ 1, A ~ 1, t0 ~ 1),
    contrasts = list(v = list(lM = ADmat)),
    report_p_vector = FALSE
  )

  design_natural <- design(
    data = dat,
    model = RDM,
    matchfun = matchfun,
    formula = list(v ~ E + lM, B ~ 1, A ~ 1, t0 ~ 1),
    contrasts = list(v = list(lM = ADmat)),
    pre_transform_terms = list(v = c("v", "v_Ehard")),
    report_p_vector = FALSE
  )

  p_vec <- c(v = log(2), v_Ehard = log(2), v_lMd = 0.5, B = log(1.5), A = 0, t0 = log(0.2))

  default_map <- mapped_pars(design_default, p_vec, digits = 10)
  natural_map <- mapped_pars(design_natural, p_vec, digits = 10)

  default_vals <- with(default_map, setNames(v, paste(E, lM, sep = "::")))
  natural_vals <- with(natural_map, setNames(v, paste(E, lM, sep = "::")))

  expect_equal(natural_vals[["easy::FALSE"]], 2 - 0.25, tolerance = 1e-8)
  expect_equal(natural_vals[["easy::TRUE"]], 2 + 0.25, tolerance = 1e-8)
  expect_equal(natural_vals[["hard::FALSE"]], 4 - 0.25, tolerance = 1e-8)
  expect_equal(natural_vals[["hard::TRUE"]], 4 + 0.25, tolerance = 1e-8)

  expect_equal(default_vals[["easy::FALSE"]], exp(log(2) - 0.25), tolerance = 1e-8)
  expect_equal(default_vals[["hard::TRUE"]], exp(log(2) + log(2) + 0.25), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(unname(default_vals), unname(natural_vals), tolerance = 1e-8)))
})

test_that("pre_transform_terms only rejects natural-scale negatives when bounds require positivity", {
  dat <- data.frame(
    subjects = factor(c(1, 1, 1, 1)),
    S = factor(c("left", "left", "right", "right")),
    E = factor(c("easy", "hard", "easy", "hard"), levels = c("easy", "hard")),
    R = factor(c("left", "left", "right", "right"), levels = c("left", "right")),
    rt = c(0.5, 0.6, 0.55, 0.65)
  )
  ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
  matchfun <- function(d) d$S == d$lR

  design_pos <- design(
    data = dat,
    model = RDM,
    matchfun = matchfun,
    formula = list(v ~ 1, B ~ E + lM, A ~ 1, t0 ~ 1),
    contrasts = list(B = list(lM = ADmat)),
    pre_transform_terms = list(B = c("B", "B_Ehard")),
    report_p_vector = FALSE
  )

  design_free <- design(
    data = dat,
    model = RDM,
    matchfun = matchfun,
    formula = list(v ~ 1, B ~ E + lM, A ~ 1, t0 ~ 1),
    contrasts = list(B = list(lM = ADmat)),
    pre_transform_terms = list(B = c("B", "B_Ehard")),
    bound = list(minmax = cbind(B = c(-Inf, Inf))),
    report_p_vector = FALSE
  )

  p_vec <- c(v = log(1), B = log(0.2), B_Ehard = log(2), B_lMd = 1, A = 0, t0 = log(0.2))

  dadm_pos <- EMC2:::design_model(dat, design_pos, compress = FALSE, verbose = FALSE)
  dadm_free <- EMC2:::design_model(dat, design_free, compress = FALSE, verbose = FALSE)

  pars_pos <- EMC2:::get_pars_matrix_oo(p_vec, dadm_pos, design_pos$model())
  pars_free <- EMC2:::get_pars_matrix_oo(p_vec, dadm_free, design_free$model())

  map_free <- mapped_pars(design_free, p_vec, digits = 10)
  free_vals <- with(map_free, setNames(B, paste(E, lM, sep = "::")))
  expect_equal(free_vals[["easy::FALSE"]], 0.2 - 0.5, tolerance = 1e-8)
  expect_lt(free_vals[["easy::FALSE"]], 0)

  ok_pos <- attr(pars_pos, "ok")
  ok_free <- attr(pars_free, "ok")
  pos_key <- paste(dadm_pos$E, dadm_pos$lM, sep = "::")
  free_key <- paste(dadm_free$E, dadm_free$lM, sep = "::")

  expect_true(any(pos_key == "easy::FALSE"))
  expect_false(any(ok_pos[pos_key == "easy::FALSE"]))
  expect_true(all(ok_free[free_key == "easy::FALSE"]))
})

test_that("make_data returns trialwise parameters that respect pre_transform_terms", {
  dat <- data.frame(
    subjects = factor(c(1, 1, 1, 1)),
    S = factor(c("left", "left", "right", "right")),
    E = factor(c("easy", "hard", "easy", "hard"), levels = c("easy", "hard")),
    R = factor(c("left", "left", "right", "right"), levels = c("left", "right")),
    rt = c(0.5, 0.6, 0.55, 0.65)
  )
  ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
  matchfun <- function(d) d$S == d$lR

  design_nat <- design(
    data = dat,
    model = RDM,
    matchfun = matchfun,
    formula = list(v ~ E + lM, B ~ 1, A ~ 1, t0 ~ 1),
    contrasts = list(v = list(lM = ADmat)),
    pre_transform_terms = list(v = c("v", "v_Ehard")),
    report_p_vector = FALSE
  )

  p_vec <- c(v = log(2), v_Ehard = log(2), v_lMd = 0.5, B = log(1.5), A = 0, t0 = log(0.2))
  sim <- make_data(
    p_vec,
    design = design_nat,
    data = dat,
    return_trialwise_parameters = TRUE
  )

  twp <- attr(sim, "trialwise_parameters")
  dadm <- EMC2:::design_model(
    EMC2:::add_accumulators(dat, design_nat$matchfun, simulate = TRUE, type = design_nat$model()$type),
    design_nat, design_nat$model, add_acc = FALSE, compress = FALSE, verbose = FALSE, rt_check = FALSE
  )
  key <- paste(dadm$E, dadm$lM, sep = "::")
  vals <- setNames(twp[, "v"], key)

  expect_equal(vals[["easy::FALSE"]], 2 - 0.25, tolerance = 1e-8)
  expect_equal(vals[["easy::TRUE"]], 2 + 0.25, tolerance = 1e-8)
  expect_equal(vals[["hard::FALSE"]], 4 - 0.25, tolerance = 1e-8)
  expect_equal(vals[["hard::TRUE"]], 4 + 0.25, tolerance = 1e-8)
})

test_that("predict.emc propagates pre_transform_terms through make_data", {
  dat <- data.frame(
    subjects = factor(c(1, 1, 1, 1)),
    S = factor(c("left", "left", "right", "right")),
    E = factor(c("easy", "hard", "easy", "hard"), levels = c("easy", "hard")),
    R = factor(c("left", "left", "right", "right"), levels = c("left", "right")),
    rt = c(0.5, 0.6, 0.55, 0.65)
  )
  ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
  matchfun <- function(d) d$S == d$lR

  design_nat <- design(
    data = dat,
    model = RDM,
    matchfun = matchfun,
    formula = list(v ~ E + lM, B ~ 1, A ~ 1, t0 ~ 1),
    contrasts = list(v = list(lM = ADmat)),
    pre_transform_terms = list(v = c("v", "v_Ehard")),
    report_p_vector = FALSE
  )

  p_vec <- c(v = log(2), v_Ehard = log(2), v_lMd = 0.5, B = log(1.5), A = 0, t0 = log(0.2), s = log(1))
  make_data_called <- FALSE

  emc_obj <- samples_LNR

  pred <- testthat::with_mocked_bindings(
    {
      EMC2:::predict.emc(
        emc_obj,
        n_post = 1,
        stat = "mean",
        n_cores = 1,
        return_trialwise_parameters = TRUE
      )
    },
    get_design = function(x) list(design_nat),
    make_data = function(parameters, design, data, ...) {
      make_data_called <<- TRUE
      expect_identical(design$model()$transform$pre_sum_terms$v, c("v", "v_Ehard"))
      out <- data
      out$rt <- 0.5
      out$R <- factor(data$R, levels = levels(data$R))
      out
    },
    .env = asNamespace("EMC2")
  )
  expect_true(make_data_called)
  expect_true(is.data.frame(pred))
})

test_that("prior map=TRUE path respects pre_transform_terms in mapped values", {
  dat <- data.frame(
    subjects = factor(c(1, 1, 1, 1)),
    S = factor(c("left", "left", "right", "right")),
    E = factor(c("easy", "hard", "easy", "hard"), levels = c("easy", "hard")),
    R = factor(c("left", "left", "right", "right"), levels = c("left", "right")),
    rt = c(0.5, 0.6, 0.55, 0.65)
  )
  ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
  matchfun <- function(d) d$S == d$lR

  design_nat <- design(
    data = dat,
    model = RDM,
    matchfun = matchfun,
    formula = list(v ~ E + lM, B ~ 1, A ~ 1, t0 ~ 1),
    contrasts = list(v = list(lM = ADmat)),
    pre_transform_terms = list(v = c("v", "v_Ehard")),
    report_p_vector = FALSE
  )

  mu_mean <- c(v = log(2), v_Ehard = log(2), v_lMd = 0.5, B = log(1.5), A = 0, t0 = log(0.2), s = log(1))
  mu_sd <- setNames(rep(1e-6, length(mu_mean)), names(mu_mean))
  pr <- prior(design_nat, type = "single", mu_mean = mu_mean, mu_sd = mu_sd)

  mapped <- parameters(pr, selection = "alpha", N = 8, map = TRUE)

  expect_equal(mean(mapped$v_Eeasy_lMFALSE), 2 - 0.25, tolerance = 1e-3)
  expect_equal(mean(mapped$v_Eeasy_lMTRUE), 2 + 0.25, tolerance = 1e-3)
  expect_equal(mean(mapped$v_Ehard_lMFALSE), 4 - 0.25, tolerance = 1e-3)
  expect_equal(mean(mapped$v_Ehard_lMTRUE), 4 + 0.25, tolerance = 1e-3)
})
