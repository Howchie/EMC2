# Extracted from test-map.R:20

# prequel ----------------------------------------------------------------------
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
dat <- forstmann
dat$covariate <- 1:nrow(forstmann)
des <- design(data = dat, formula = list(v ~ covariate*E, B ~ E, t0 ~ S),
              model = LBA)

# test -------------------------------------------------------------------------
expect_snapshot(credint(samples_LNR, selection = "mu", map = "E"))
