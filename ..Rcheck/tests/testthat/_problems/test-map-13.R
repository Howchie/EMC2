# Extracted from test-map.R:13

# prequel ----------------------------------------------------------------------
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
dat <- forstmann
dat$covariate <- 1:nrow(forstmann)
des <- design(data = dat, formula = list(v ~ covariate*E, B ~ E, t0 ~ S),
              model = LBA)

# test -------------------------------------------------------------------------
expect_snapshot(mapped_pars(des))
expect_snapshot(mapped_pars(des, p_vector= rnorm(length(sampled_pars(des)))))
expect_snapshot(mapped_pars(prior(des, mu_mean = c('v_covariate'  = 1))))
