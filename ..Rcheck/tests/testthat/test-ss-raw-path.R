# The stop-signal branch of calc_ll_oo evaluates likelihoods directly from
# ParamTable column pointers (src/ss_raw.h) and only falls back to
# materializing a per-trial NumericMatrix when a required column is missing
# or the EMC2_SS_FORCE_MATERIALIZE env var is set. These tests pin the two
# paths to each other: any divergence means the raw kernels no longer match
# the materialized c_log_likelihood_ss implementation.

ss_ll_both_paths <- function(p_vector, design_obj, n_trials, seed,
                             TC = NULL, make_data_args = list(),
                             n_particles = 25) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  args <- c(list(p_vector, design_obj, n_trials = n_trials), make_data_args)
  if (!is.null(TC)) args$TC <- TC
  dat <- do.call(make_data, args)
  emc <- make_emc(dat, design_obj, type = "single", verbose = FALSE)

  set.seed(seed + 1)
  p_mat <- matrix(rep(as.numeric(p_vector), each = n_particles),
                  nrow = n_particles) +
    matrix(rnorm(n_particles * length(p_vector), sd = 0.25), nrow = n_particles)
  colnames(p_mat) <- names(p_vector)

  ll_raw <- EMC2:::calc_ll_manager(p_mat, emc[[1]]$data[[1]], emc[[1]]$model)
  Sys.setenv(EMC2_SS_FORCE_MATERIALIZE = "1")
  on.exit(Sys.unsetenv("EMC2_SS_FORCE_MATERIALIZE"), add = TRUE)
  ll_mat <- EMC2:::calc_ll_manager(p_mat, emc[[1]]$data[[1]], emc[[1]]$model)
  list(raw = as.numeric(ll_raw), mat = as.numeric(ll_mat))
}

test_that("SSEXG raw path matches materialized fallback (staircase, UC+LC)", {
  lIfun_go <- function(d) factor(rep(2, nrow(d)), levels = 1:2)
  designSSE <- design(
    model = SSEXG, factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = function(d) as.character(d$S) == as.character(d$lR),
    functions = list(lI = lIfun_go),
    covariates = "SSD",
    formula = list(mu ~ 0 + lM, sigma ~ 1, tau ~ 1,
                   muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1))
  pv <- sampled_pars(designSSE, doMap = FALSE)
  pv[c("mu_lMFALSE","mu_lMTRUE","sigma","tau","muS","sigmaS","tauS","gf","tf")] <-
    c(log(0.65), log(0.50), log(0.05), log(0.15), log(0.20), log(0.04), log(0.08),
      qnorm(0.08), qnorm(0.08))
  res <- ss_ll_both_paths(pv, designSSE, n_trials = 200, seed = 401,
                          TC = list(UC = 1.2, LC = 0.3),
                          make_data_args = list(functions = list(SSD = make_ssd())))
  expect_true(all(is.finite(res$raw)))
  expect_equal(res$raw, res$mat, tolerance = 1e-12)
})

test_that("SSEXG raw path matches materialized fallback (ST accumulator)", {
  lIfun_st <- function(d) factor(ifelse(as.character(d$lR) == "st", 1, 2), levels = 1:2)
  AccTypeFun <- function(d) {
    is_st <- as.character(d$lR) == "st"
    is_match <- !is_st & (as.character(d$S) == as.character(d$lR))
    factor(ifelse(is_st, "st", ifelse(is_match, "match", "mismatch")),
           levels = c("mismatch", "match", "st"))
  }
  mySSD <- function(d) SSD_function(d, SSD = c(0.20, 0.35), pSSD = c(0.25, 0.25))
  designSSEst <- design(
    model = SSEXG, factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right", "st"),
    matchfun = function(d) as.character(d$S) == as.character(d$lR),
    functions = list(lI = lIfun_st, SSD = mySSD, AccType = AccTypeFun),
    formula = list(mu ~ 0 + AccType, sigma ~ 1, tau ~ 1,
                   muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1))
  pv <- sampled_pars(designSSEst, doMap = FALSE)
  pv[c("mu_AccTypemismatch","mu_AccTypematch","mu_AccTypest","sigma","tau",
       "muS","sigmaS","tauS","gf","tf")] <-
    c(log(0.65), log(0.50), log(0.38), log(0.05), log(0.15),
      log(0.20), log(0.04), log(0.08), qnorm(0.05), qnorm(0.01))
  res <- ss_ll_both_paths(pv, designSSEst, n_trials = 200, seed = 402,
                          TC = list(UC = 1.2))
  expect_true(all(is.finite(res$raw)))
  expect_equal(res$raw, res$mat, tolerance = 1e-12)
})

test_that("SSRDEX raw path matches materialized fallback", {
  lIfun_go <- function(d) factor(rep(2, nrow(d)), levels = 1:2)
  mySSD <- function(d) SSD_function(d, SSD = c(0.20, 0.35), pSSD = c(0.25, 0.25))
  designSSR <- design(
    model = SSRDEX, factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = function(d) as.character(d$S) == as.character(d$lR),
    functions = list(lI = lIfun_go, SSD = mySSD),
    formula = list(v ~ 0 + lM, B ~ 1, A ~ 1, t0 ~ 1,
                   muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1))
  pv <- sampled_pars(designSSR, doMap = FALSE)
  pv[c("v_lMFALSE","v_lMTRUE","B","A","t0","muS","sigmaS","tauS","gf","tf")] <-
    c(log(0.2), log(1.5), log(1), log(0.4), log(0.15),
      log(0.2), log(0.08), log(0.1), qnorm(0.08), qnorm(0.08))
  res <- ss_ll_both_paths(pv, designSSR, n_trials = 200, seed = 403,
                          TC = list(UC = 1.5))
  expect_true(all(is.finite(res$raw)))
  expect_equal(res$raw, res$mat, tolerance = 1e-12)
})
