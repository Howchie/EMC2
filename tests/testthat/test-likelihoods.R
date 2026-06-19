RNGkind("L'Ecuyer-CMRG")
set.seed(123)

calc_lls <- function(emc, n_particles=1e3) {
  calc_ll_context(emc, n_particles)$new
}

calc_ll_context <- function(emc, n_particles=1e3, particle_sd=1) {
  model <- emc[[1]]$model()
  p_types <- names(model$p_types)
  dadm <- emc[[1]]$data[[1]]

  designs <- list()
  for(p in p_types){
    designs[[p]] <- attr(dadm,"designs")[[p]][attr(attr(dadm,"designs")[[p]],"expand"),,drop=FALSE]
  }
  constants <- attr(dadm, "constants")
  if(is.null(constants)) constants <- NA

  # make p_mat

  p_mat <- matrix(rnorm(n_particles*length(p_types), sd = particle_sd), ncol=length(p_types))
  colnames(p_mat) <- p_types

  lls_old <- EMC2:::calc_ll(p_mat, dadm, constants = constants, designs = designs, type = model$c_name,
                            model$bound, model$transform, model$pre_transform, p_types = p_types,
                            min_ll = log(1e-10), model$trend)
  lls_new <- EMC2:::calc_ll_oo(p_mat, dadm, constants = constants, designs = designs, type = model$c_name,
                               model$bound, model$transform, model$pre_transform, p_types = p_types,
                               min_ll = log(1e-10), model$trend)
  list(old = lls_old, new = lls_new, dadm = dadm, p_mat = p_mat,
       constants = constants, designs = designs, model = model, p_types = p_types)
}

time_min <- function(expr, reps=5) {
  expr <- substitute(expr)
  envir <- parent.frame()
  min(vapply(seq_len(reps), function(i) system.time(eval(expr, envir))[["elapsed"]], numeric(1)))
}

censored_rt_count <- function(emc) {
  sum(!is.finite(emc[[1]]$data[[1]]$rt))
}

expect_uc_path <- function(base, mild, moderate, label, n_particles=5000,
                           mild_input_censored, moderate_input_censored,
                           mild_limit=8, moderate_limit=25) {
  set.seed(12345)
  base_ctx <- calc_ll_context(base, n_particles, particle_sd = 0.25)
  set.seed(12345)
  mild_ctx <- calc_ll_context(mild, n_particles, particle_sd = 0.25)
  set.seed(12345)
  moderate_ctx <- calc_ll_context(moderate, n_particles, particle_sd = 0.25)

  expect_true(all(is.finite(mild_ctx$new)))
  expect_true(all(is.finite(moderate_ctx$new)))
  mild_ok <- is.finite(mild_ctx$old) & is.finite(mild_ctx$new)
  moderate_ok <- is.finite(moderate_ctx$old) & is.finite(moderate_ctx$new)
  expect_true(any(mild_ok))
  expect_true(any(moderate_ok))
  expect_lt(max(abs(mild_ctx$old[mild_ok] - mild_ctx$new[mild_ok])), 1e-8)
  expect_lt(max(abs(moderate_ctx$old[moderate_ok] - moderate_ctx$new[moderate_ok])), 1e-8)
  expect_gt(censored_rt_count(mild), 0)
  expect_gte(censored_rt_count(moderate), censored_rt_count(mild))
  expect_gt(mild_input_censored, 0)
  expect_gt(moderate_input_censored, mild_input_censored)

  gc()
  base_time <- time_min(EMC2:::calc_ll_oo(base_ctx$p_mat, base_ctx$dadm,
                                          constants = base_ctx$constants,
                                          designs = base_ctx$designs,
                                          type = base_ctx$model$c_name,
                                          bounds = base_ctx$model$bound,
                                          transforms = base_ctx$model$transform,
                                          pretransforms = base_ctx$model$pre_transform,
                                          p_types = base_ctx$p_types,
                                          min_ll = log(1e-10),
                                          trend = base_ctx$model$trend))
  gc()
  mild_time <- time_min(EMC2:::calc_ll_oo(mild_ctx$p_mat, mild_ctx$dadm,
                                          constants = mild_ctx$constants,
                                          designs = mild_ctx$designs,
                                          type = mild_ctx$model$c_name,
                                          bounds = mild_ctx$model$bound,
                                          transforms = mild_ctx$model$transform,
                                          pretransforms = mild_ctx$model$pre_transform,
                                          p_types = mild_ctx$p_types,
                                          min_ll = log(1e-10),
                                          trend = mild_ctx$model$trend))
  gc()
  moderate_time <- time_min(EMC2:::calc_ll_oo(moderate_ctx$p_mat, moderate_ctx$dadm,
                                              constants = moderate_ctx$constants,
                                              designs = moderate_ctx$designs,
                                              type = moderate_ctx$model$c_name,
                                              bounds = moderate_ctx$model$bound,
                                              transforms = moderate_ctx$model$transform,
                                              pretransforms = moderate_ctx$model$pre_transform,
                                              p_types = moderate_ctx$p_types,
                                              min_ll = log(1e-10),
                                              trend = moderate_ctx$model$trend))

  expect_lt(mild_time / base_time, mild_limit)
  expect_lt(moderate_time / base_time, moderate_limit)
}

# Simplest design, no trend -----------------------------------------------
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
matchfun=function(d)d$S==d$lR
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2],]
dat$subjects <- droplevels(dat$subjects)

# LNR ---------------------------------------------------------------------
design_LNR <- design(data = dat,model=LNR,matchfun=matchfun,
                     formula=list(m~lM,s~1,t0~1),
                     contrasts=list(m=list(lM=ADmat)))
LNR_s1 <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 2, compress=FALSE)

test_that("LNR", {
  expect_snapshot(calc_lls(LNR_s1))
})

LNR_s2 <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 2, compress=TRUE)
test_that("LNR compressed", {
  expect_snapshot(calc_lls(LNR_s2))
})



# RDM ----------------------------------------------------------
design_RDM <- design(data = dat,model=RDM,matchfun=matchfun,
                     formula=list(v~lM,s~1,t0~1,A~1,B~1),
                     contrasts=list(v=list(lM=ADmat)))
RDM_s <- make_emc(dat, design_RDM, rt_resolution = 0.05, n_chains = 2, compress=FALSE)
test_that("RDM", {
  expect_snapshot(calc_lls(RDM_s))
})
RDM_s2 <- make_emc(dat, design_RDM, rt_resolution = 0.05, n_chains = 2, compress=TRUE)
test_that("RDM compressed", {
  expect_snapshot(calc_lls(RDM_s2))
})


# LBA ----------------------------------------------------------
design_LBA <- design(data = dat,model=LBA,matchfun=matchfun,
                     formula=list(v~lM,sv~1,t0~1,A~1,B~1),
                     contrasts=list(v=list(lM=ADmat)))
LBA_s <- make_emc(dat, design_LBA, rt_resolution = 0.05, n_chains = 2, compress=FALSE)
test_that("LBA", {
  expect_snapshot(calc_lls(LBA_s))
})
LBA_s2 <- make_emc(dat, design_LBA, rt_resolution = 0.05, n_chains = 2, compress=TRUE)
test_that("LBA compressed", {
  expect_snapshot(calc_lls(LBA_s2))
})

test_that("LBAIO keeps intrinsic omission mass at infinity", {
  design_lbaio <- design(
    factors = list(subjects = "s1"),
    Rlevels = c("yes"),
    model = LBA(posdrift = FALSE),
    formula = list(v ~ 1, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
    constants = c(sv = log(1), B = log(1), A = log(0.5), t0 = log(0.2))
  )

  pars <- matrix(
    c(v = 0.5, sv = 1, B = 1, A = 0.5, t0 = 0.2, b = 1.5),
    nrow = 1,
    dimnames = list(NULL, c("v", "sv", "B", "A", "t0", "b"))
  )
  expect_equal(
    EMC2:::pLBA(Inf, pars, posdrift = FALSE),
    pnorm(0, mean = 0.5, sd = 1, lower.tail = FALSE),
    tolerance = 1e-12
  )

  dat <- data.frame(
    subjects = factor("s1"),
    R = factor(NA, levels = "yes"),
    rt = Inf,
    UC = Inf
  )
  emc <- make_emc(dat, design_lbaio, type = "single", n_chains = 1, compress = FALSE, verbose = FALSE)
  dadm <- emc[[1]]$data[[1]]
  model <- emc[[1]]$model()
  p_mat <- matrix(0.5, nrow = 1, dimnames = list(NULL, "v"))
  designs <- list(v = attr(dadm, "designs")$v[attr(attr(dadm, "designs")$v, "expand"), , drop = FALSE])
  ll <- EMC2:::calc_ll_oo(
    p_mat, dadm,
    constants = attr(dadm, "constants"),
    designs = designs,
    type = model$c_name,
    bounds = model$bound,
    transforms = model$transform,
    pretransforms = model$pre_transform,
    p_types = names(model$p_types),
    min_ll = log(1e-10),
    trend = model$trend
  )

  expect_equal(as.numeric(ll), pnorm(0, 0.5, 1, log.p = TRUE), tolerance = 1e-8)

  dat_known <- data.frame(
    subjects = factor("s1"),
    R = factor("yes", levels = "yes"),
    rt = Inf,
    UC = Inf
  )
  emc_known <- make_emc(dat_known, design_lbaio, type = "single", n_chains = 1,
                        compress = FALSE, verbose = FALSE)
  dadm_known <- emc_known[[1]]$data[[1]]
  ll_known <- EMC2:::calc_ll_oo(
    p_mat, dadm_known,
    constants = attr(dadm_known, "constants"),
    designs = list(v = attr(dadm_known, "designs")$v[
      attr(attr(dadm_known, "designs")$v, "expand"), , drop = FALSE
    ]),
    type = model$c_name,
    bounds = model$bound,
    transforms = model$transform,
    pretransforms = model$pre_transform,
    p_types = names(model$p_types),
    min_ll = log(1e-10),
    trend = model$trend
  )

  expect_equal(as.numeric(ll_known), pnorm(0, 0.5, 1, log.p = TRUE), tolerance = 1e-8)

  set.seed(1001)
  sim <- make_data(p_mat[1, ], design_lbaio, n_trials = 200, TC = list(UC = Inf))
  expect_true(any(is.infinite(sim$rt)))
  expect_true(all(is.na(sim$R[is.infinite(sim$rt)])))

  dat_trunc <- data.frame(
    subjects = factor("s1"),
    R = factor("yes", levels = "yes"),
    rt = 0.6,
    LT = 0.3,
    UC = Inf
  )
  emc_trunc <- make_emc(dat_trunc, design_lbaio, type = "single", n_chains = 1,
                        compress = FALSE, verbose = FALSE)
  dat_notrunc <- dat_trunc
  dat_notrunc$LT <- 0
  emc_notrunc <- make_emc(dat_notrunc, design_lbaio, type = "single", n_chains = 1,
                          compress = FALSE, verbose = FALSE)
  ll_trunc <- EMC2:::calc_ll_manager(p_mat, emc_trunc[[1]]$data[[1]], emc_trunc[[1]]$model)
  ll_notrunc <- EMC2:::calc_ll_manager(p_mat, emc_notrunc[[1]]$data[[1]], emc_notrunc[[1]]$model)
  S_LT <- 1 - EMC2:::pLBA(0.3, pars, posdrift = FALSE)
  S_Inf <- pnorm(0, mean = 0.5, sd = 1)
  expect_equal(as.numeric(ll_trunc - ll_notrunc), -log(S_LT - S_Inf), tolerance = 1e-8)
})


# WDM ----------------------------------------------------------
design_WDM <- design(data = dat, model=DDM,
                     formula =list(v~1,a~1, t0~1, s~1, Z~1, sv~1, SZ~1),
                     constants=c(s=1, sv=log(0), SZ=qnorm(0)))
WDM_s <- make_emc(dat, design_WDM, rt_resolution = 0.05, n_chains = 2)
test_that("WDM", {
  expect_snapshot(calc_lls(WDM_s))
})
WDM_s2 <- make_emc(dat, design_WDM, rt_resolution = 0.05, n_chains = 2, compress=TRUE)
test_that("WDM compressed", {
  expect_snapshot(calc_lls(WDM_s2))
})


# DDM ---
design_DDM <- design(data = dat,model=DDM,
                     formula =list(v~1,a~1, t0~1, s~1, Z~1, sv~1, SZ~1, st0~1),
                     constants=c(s=1))
DDM_s <- make_emc(dat, design_DDM, rt_resolution = 0.05, n_chains = 2)
test_that("DDM", {
  expect_snapshot(calc_lls(DDM_s))
})

test_that("DDM finite RT path remains fast with mild censoring", {
  skip_on_cran()

  uc_rt_resolution <- 1/60
  DDM_uc_base <- make_emc(dat, design_DDM, rt_resolution = uc_rt_resolution, n_chains = 2)
  DDM_mild_dat <- make_missing(dat, UC = 1.05, verbose = FALSE, rt_resolution = uc_rt_resolution)
  DDM_moderate_dat <- make_missing(dat, UC = 0.7, verbose = FALSE, rt_resolution = uc_rt_resolution)
  DDM_mild_uc <- make_emc(DDM_mild_dat,
                          design_DDM, rt_resolution = uc_rt_resolution, n_chains = 2)
  DDM_moderate_uc <- make_emc(DDM_moderate_dat,
                              design_DDM, rt_resolution = uc_rt_resolution, n_chains = 2)

  expect_uc_path(DDM_uc_base, DDM_mild_uc, DDM_moderate_uc, "DDM",
                 mild_input_censored = sum(!is.finite(DDM_mild_dat$rt)),
                 moderate_input_censored = sum(!is.finite(DDM_moderate_dat$rt)),
                 mild_limit = 8, moderate_limit = 25)
})

test_that("LNR finite RT path remains fast with mild censoring", {
  skip_on_cran()

  uc_rt_resolution <- 1/60
  LNR_uc_base <- make_emc(dat, design_LNR, rt_resolution = uc_rt_resolution, n_chains = 2, compress = TRUE)
  LNR_mild_dat <- make_missing(dat, UC = 1.05, verbose = FALSE, rt_resolution = uc_rt_resolution)
  LNR_moderate_dat <- make_missing(dat, UC = 0.7, verbose = FALSE, rt_resolution = uc_rt_resolution)
  LNR_mild_uc <- make_emc(LNR_mild_dat,
                          design_LNR, rt_resolution = uc_rt_resolution, n_chains = 2, compress = TRUE)
  LNR_moderate_uc <- make_emc(LNR_moderate_dat,
                              design_LNR, rt_resolution = uc_rt_resolution, n_chains = 2, compress = TRUE)

  expect_uc_path(LNR_uc_base, LNR_mild_uc, LNR_moderate_uc, "LNR",
                 mild_input_censored = sum(!is.finite(LNR_mild_dat$rt)),
                 moderate_input_censored = sum(!is.finite(LNR_moderate_dat$rt)),
                 mild_limit = 8, moderate_limit = 25)
})
