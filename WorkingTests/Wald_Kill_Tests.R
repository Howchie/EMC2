rm(list = ls())
library(EMC2)
library(dplyr)

designRDMSWTN <- design(
    factors = list(
        S = "Target", subjects = 1, L = c("L", "M", "H")
    ),
    Rlevels = c("Go"),
    formula = list(v ~ L, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1, lambda ~ 1),
    constants = c(s = log(1)),
    model = RDMSWTN(erlang_type = "global_kill"), UC = 3
)

p_vec <- sampled_pars(designRDMSWTN)
p_vec[1:3] = log(c(3, 2, 1.5))
p_vec["B"] <- log(1.1)
p_vec["A"] <- log(0.4)
p_vec["t0"] <- log(0.2)
p_vec["sv"] <- log(0)
p_vec["lambda"] <- log(0)

dat <- make_data(p_vec, design = designRDMSWTN, n_trials = 500)
mean(is.finite(dat$rt))
mean(dat$rt[is.finite(dat$rt)])
likelihood(p_vec, data = NULL, model = model)
