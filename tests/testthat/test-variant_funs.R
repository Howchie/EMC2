local_rng_guard()  # see helper-rng.R: keep this file's RNG changes inside it
# When working with lM it is useful to design  an "average and difference"
# contrast matrix, which for binary responses has a simple canonical from:
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
# We also define a match function for lM
matchfun=function(d)d$S==d$lR

# Drop most subjects
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2],]
dat$subjects <- droplevels(dat$subjects)

design_LNR <- design(data = dat,model=LNR,matchfun=matchfun,
                          formula=list(m~lM,s~1,t0~1),
                          contrasts=list(m=list(lM=ADmat)))



# Set seed for mclapply
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
# Run the models for the different variants
# Only preburn, so bit risky since adapt could fail
N <- 25
idx <- N + 1

test_that("run_diag", {
  skip_on_os("windows")

  LNR_diag <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 2, type = "diagonal")
  LNR_diag <- init_chains(LNR_diag, cores_for_chains = 1, particles = 5)
  LNR_diag <- run_emc(LNR_diag, cores_for_chains = 1, stop_criteria = list(iter = N), stage = "preburn", particle_factor = 10)


  expect_snapshot(
    LNR_diag[[1]]$samples$alpha[,,idx], variant = Sys.info()[1]
  )
  expect_snapshot(
    LNR_diag[[1]]$samples$theta_mu[,idx], variant = Sys.info()[1]
  )
  expect_snapshot(
    LNR_diag[[1]]$samples$theta_var[,,idx], variant = Sys.info()[1]
  )
  expect_snapshot(
    compare(list(diag = LNR_diag),
            stage = "preburn", cores_for_props = 1),
    variant = Sys.info()[1]
  )
})

test_that("run_blocked", {
  skip_on_os("windows")

  LNR_blocked <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 2, type = "blocked",
                          par_group = c(1,2,3,3))
  LNR_blocked <- init_chains(LNR_blocked, cores_for_chains = 1, particles = 5)
  LNR_blocked <- run_emc(LNR_blocked, cores_for_chains = 1, stop_criteria = list(iter = N), stage = "preburn", particle_factor = 10)

  expect_snapshot(
    LNR_blocked[[1]]$samples$alpha[,,idx], variant = Sys.info()[1]
  )
  expect_snapshot(
    LNR_blocked[[1]]$samples$theta_mu[,idx], variant = Sys.info()[1]
  )
  expect_snapshot(
    LNR_blocked[[1]]$samples$theta_var[,,idx], variant = Sys.info()[1]
  )
  # Information criteria only: bridge sampling on a 25-iteration preburn fails
  # to converge for some seeds, and blocked bridge sampling has its own exact
  # tests in test-group-bridge.R.
  expect_snapshot(
    compare(list(blocked = LNR_blocked), BayesFactor = FALSE,
            stage = "preburn", cores_for_props = 1),
    variant = Sys.info()[1]
  )
})

test_that("run_single", {
  skip_on_os("windows")

  LNR_single <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 2, type = "single")
  LNR_single <- init_chains(LNR_single, cores_for_chains = 1, particles = 5)
  LNR_single <- run_emc(LNR_single, cores_for_chains = 1, stop_criteria = list(iter = N), stage = "preburn", particle_factor = 10)

  expect_snapshot(
    LNR_single[[1]]$samples$alpha[,,idx], variant = Sys.info()[1]
  )
  expect_snapshot(
    compare(list(single = LNR_single),
            stage = "preburn", cores_for_props = 1),
    variant = Sys.info()[1]
  )
})

test_that("a one-parameter prior keeps its supplied mean variance", {
  des <- suppressMessages(design(data = dat, model = LNR, matchfun = matchfun,
                                 formula = list(m ~ 1, s ~ 1, t0 ~ 1),
                                 constants = c(s = log(0.5), t0 = log(0.2))))
  for (type in c("diagonal-gamma", "standard")) {
    pr <- prior(des, type = type, theta_mu_var = matrix(2.25, 1, 1))
    expect_equal(dim(pr$theta_mu_var), c(1L, 1L))
    expect_equal(pr$theta_mu_var[1, 1], 2.25)
  }
})

test_that("a one-parameter hierarchy runs its group update", {
  des <- suppressMessages(design(data = dat, model = LNR, matchfun = matchfun,
                                 formula = list(m ~ 1, s ~ 1, t0 ~ 1),
                                 constants = c(s = log(0.5), t0 = log(0.2))))
  for (type in c("diagonal-gamma", "standard")) {
    set.seed(11)
    emc <- suppressMessages(make_emc(dat, des, n_chains = 1L, type = type))
    emc <- suppressMessages(run_emc(emc, stage = "preburn", cores_for_chains = 1,
      cores_per_chain = 1, step_size = 5, max_tries = 1, verbose = FALSE,
      verboseProgress = FALSE,
      stop_criteria = list(iter = 5, max_gd = Inf, min_unique = 0, min_es = 0)))
    s <- emc[[1]]$samples
    expect_equal(s$idx, 6)
    expect_true(all(is.finite(s$theta_mu[1, seq_len(s$idx)])))
    expect_gt(length(unique(s$theta_mu[1, seq_len(s$idx)])), 1L)
    for (sel in c("mu", "sigma2")) {
      draws <- get_pars(emc, selection = sel, stage = "preburn")
      expect_identical(ncol(as.matrix(draws[[1]])), 1L)
      utils::capture.output(sm <- summary(emc, selection = sel, stage = "preburn"))
      expect_equal(dim(sm[[1]]), c(1L, 5L))
    }
  }
})

test_that("last_sample keeps a one-parameter group covariance as a 1 x 1 matrix", {
  store <- list(idx = 2L, theta_mu = matrix(c(0.1, 0.2), 1, 2, dimnames = list("v", NULL)),
                theta_var = array(c(1, 2), c(1, 1, 2), dimnames = list("v", "v", NULL)),
                a_half = matrix(c(1, 1), 1, 2), last_theta_var_inv = matrix(0.5, 1, 1))
  for (f in c("last_sample_standard", "last_sample_diag_gamma")) {
    last <- getFromNamespace(f, "EMC2")(store)
    expect_true(is.matrix(last$tvar))
    expect_equal(last$tvar, matrix(2, 1, 1, dimnames = list("v", "v")))
  }
})

test_that("the native factor group draw is the R form's draw", {
  skip_on_os("windows")
  for (k in 1:2) {
    fac <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 1, type = "factor",
                    n_factors = k)
    fac <- init_chains(fac, cores_for_chains = 1, particles = 5)
    sampler <- fac[[1]]
    alpha <- sampler$samples$alpha[, , sampler$samples$idx]
    set.seed(11)
    seed <- .Random.seed
    r_form <- EMC2:::gibbs_step_factor_r(sampler, alpha)
    after <- .Random.seed
    assign(".Random.seed", seed, envir = globalenv())
    native <- EMC2:::.emc_gibbs_factor_native(sampler, alpha)
    expect_false(is.null(native))
    expect_identical(native, r_form)
    expect_identical(.Random.seed, after)
  }
})

test_that("the native infinite-factor group draw is the R form's draw", {
  skip_on_os("windows")
  for (k in 1:2) {
    fac <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 1,
                    type = "infnt_factor", n_factors = k)
    fac <- init_chains(fac, cores_for_chains = 1, particles = 5)
    sampler <- fac[[1]]
    alpha <- sampler$samples$alpha[, , sampler$samples$idx]
    set.seed(11)
    seed <- .Random.seed
    r_form <- EMC2:::gibbs_step_infnt_factor_r(sampler, alpha)
    after <- .Random.seed
    assign(".Random.seed", seed, envir = globalenv())
    native <- EMC2:::.emc_gibbs_infnt_factor_native(sampler, alpha)
    expect_false(is.null(native))
    expect_identical(native, r_form)
    expect_identical(.Random.seed, after)
  }
})

test_that("the native multivariate normal draw is mvtnorm's", {
  set.seed(4)
  for (i in 1:40) {
    k <- sample(1:5, 1)
    n <- sample(c(1, 3, 12), 1)
    sigma <- crossprod(matrix(rnorm(k * (k + 2)), k + 2, k))
    mean <- if (i %% 2) rnorm(k)
    seed <- .Random.seed
    theirs <- if (is.null(mean)) mvtnorm::rmvnorm(n, sigma = sigma) else
      mvtnorm::rmvnorm(n, mean, sigma)
    after <- .Random.seed
    assign(".Random.seed", seed, envir = globalenv())
    native <- EMC2:::emc_rmvnorm_cpp(n, mean, sigma)
    expect_identical(native, theirs)
    expect_identical(.Random.seed, after)
  }
  # A covariance mvtnorm refuses is left to it.
  expect_null(EMC2:::emc_rmvnorm_cpp(1L, NULL, matrix(c(1, 0.5, 0, 1), 2)))
})

test_that("an SEM without covariates takes a group draw", {
  skip_on_os("windows")
  L <- matrix(c(1, Inf, Inf, Inf), 4, 1)
  sem <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 1, type = "SEM",
                  sem_settings = list(Lambda_mat = L))
  sem <- init_chains(sem, cores_for_chains = 1, particles = 5)
  sampler <- sem[[1]]
  set.seed(5)
  pars <- EMC2:::gibbs_step_SEM(sampler, sampler$samples$alpha[, , sampler$samples$idx])
  expect_true(all(is.finite(pars$tmu)))
  expect_true(all(is.finite(pars$tvar)))
})

test_that("the native latent ensemble rounds are the R rounds", {
  skip_on_os("windows")
  for (type in c("factor", "infnt_factor")) {
    fac <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 1, type = type,
                    n_factors = 2)
    fac <- init_chains(fac, cores_for_chains = 1, particles = 5)
    sampler <- fac[[1]]
    n <- sampler$n_subjects
    p <- sampler$n_pars
    M <- 4L
    set.seed(21)
    pars <- EMC2:::gibbs_step(sampler, sampler$samples$alpha[, , sampler$samples$idx], type)
    stack <- matrix(rnorm(n * M * p, sd = 0.5), n * M, p,
                    dimnames = list(NULL, sampler$par_names))
    base_w <- rnorm(n * M)
    k0 <- sample.int(M, n, replace = TRUE)
    cfg <- list(rounds = 3L)
    seed <- .Random.seed
    r_form <- EMC2:::.emc_ensemble_rounds_r_form(sampler, pars, stack, base_w, k0, M, cfg)
    after <- .Random.seed
    assign(".Random.seed", seed, envir = globalenv())
    native <- EMC2:::.emc_ensemble_rounds_latent(sampler, pars, stack, base_w, k0, M, cfg)
    expect_false(is.null(native))
    expect_identical(native, r_form)
    expect_identical(.Random.seed, after)
  }
})
