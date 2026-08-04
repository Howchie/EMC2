# Extracted from test-trend.R:78

# prequel ----------------------------------------------------------------------
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
matchfun=function(d)d$S==d$lR
n_trials <- 10
covariate1 <- rnorm(n_trials*2)
covariate2 <- rnorm(n_trials*2)
covariate2[1:5] <- NA
trend <- make_trend(par_names = "m",
                    cov_names = list(c("covariate1", "covariate2")),
                    kernels = "exp_incr")
design_base <- design(factors = list(subjects = 1, S = 1:2),
                      Rlevels = 1:2,
                      covariates = c('covariate1', 'covariate2'),
                      matchfun = matchfun,
                      trend = trend,
                      formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                      contrasts = list(lM = ADmat),
                      model = LNR)
p_vector <- sampled_pars(design_base, doMap = FALSE)
p_vector[1:6] <- c(-1, 1.5, log(1), log(.2), log(.2), log(.2))
dat <- make_data(p_vector, design_base, n_trials = n_trials,
                 covariates = data.frame(covariate1 = covariate1, covariate2 = covariate2))
LNR2cov <- make_emc(dat, design_base, compress = F, n_chains = 1, type = "single")
trend_2types <- make_trend(par_names = c("m", "m_lMd"),
                           cov_names = list(c("covariate1", "covariate2"), "covariate1"),
                           kernels = c("exp_incr", "pow_decr"))
design_base_shared <- design(data = dat,
                             trend = trend_2types,
                             formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                             contrasts = list(lM = ADmat),
                             matchfun = matchfun,
                             model = LNR)
LNR2cov_shared <- make_emc(dat, design_base_shared, compress = FALSE, n_chains = 1, type = "single")
trend_premap <- make_trend(
  par_names = c("m", "m_lMd"),
  cov_names = list("covariate1", "covariate2"),
  kernels = c("exp_incr", "poly2"),
  phase = "premap"
)
design_premap <- design(
  data = dat,
  trend = trend_premap,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1, m_lMd.d1 ~ lR),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR
)
LNR_premap <- make_emc(dat, design_premap, compress = FALSE, n_chains = 1, type = "single")

# test -------------------------------------------------------------------------
expect_snapshot(init_chains(LNR_premap, particles = 3, cores_per_chain = 1)[[1]]$samples)
