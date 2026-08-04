# Extracted from test-map.R:148

# prequel ----------------------------------------------------------------------
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
dat <- forstmann
dat$covariate <- 1:nrow(forstmann)
des <- design(data = dat, formula = list(v ~ covariate*E, B ~ E, t0 ~ S),
              model = LBA)

# test -------------------------------------------------------------------------
short <- .credint_label_mar("m_short", 90)
