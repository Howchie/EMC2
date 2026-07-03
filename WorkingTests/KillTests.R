rm(list = ls())
library(EMC2)
library(dplyr)

designRDMSWTN <- design(
    factors = list(
        S = "Target", subjects = 1, L = c("L", "M", "H")
    ),
    Rlevels = c("Go"),
    formula = list(v ~ L, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1, mK ~ 1),
    constants = c(s = log(1)),
    model = RDMSWTN(erlang_type = "local_kill", erlang_shape=1), UC = 3
)

designRDM <- design(
  factors = list(
    S = "Target", subjects = 1, L = c("L", "M", "H")
  ),
  Rlevels = c("Go"),
  formula = list(v ~ L, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1, mK ~ 1),
  constants = c(s = log(1), sv=log(0)),
  model = RDMSWTN(erlang_type = "local_kill", erlang_shape=1), UC = 3
)

designBAwL <- design(
  factors = list(
    S = "Target", subjects = 1, L = c("L", "M", "H")
  ),
  Rlevels = c("Go"),
  formula = list(v ~ L, B ~ 1, A ~ 1, t0 ~ 1, sv ~ L, k~1, mK ~ 1),
  constants = c(sv = log(1), k=log(0)),
  model = BAwL(erlang_type = "local_kill", erlang_shape=1), UC = 3
)

p_vec <- sampled_pars(designRDMSWTN)
p_vec[1:3] = log(c(3, 2, 1.5))
p_vec["B"] <- log(1.1)
p_vec["A"] <- log(0.4)
p_vec["t0"] <- log(0.2)
p_vec["sv"] <- log(1)
p_vec["mK"] <- log(2)

dat <- make_data(p_vec, design = designRDMSWTN, n_trials = 500)
mean(is.finite(dat$rt))
mean(dat$rt[is.finite(dat$rt)])

emc1 = make_emc(dat, designRDMSWTN, type="single")
emc1 = fit(emc1)

emc2 = make_emc(dat, designBAwL, type="single")
emc2 = fit(emc2)

emc3 = make_emc(dat, designRDM, type="single")
emc3 = fit(emc3)
