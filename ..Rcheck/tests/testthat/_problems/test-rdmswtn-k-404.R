# Extracted from test-rdmswtn-k.R:404

# test -------------------------------------------------------------------------
set.seed(1)
designRDMSWTN <- design(
    factors = list(S = "Target", subjects = 1, L = c("L", "M", "H")),
    Rlevels = c("Go"),
    formula = list(v ~ L, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1, mG ~ 1, mK ~ 1),
    constants = c(s = log(1)),
    model = RDMSWTN(erlang_type = "local_kill"),
    report_p_vector = FALSE
  )
p_vec <- sampled_pars(designRDMSWTN)
p_vec[1:3] <- log(c(3, 2, 1.5))
p_vec["B"] <- log(1.1)
p_vec["A"] <- log(0.4)
p_vec["t0"] <- log(0.2)
p_vec["sv"] <- log(0.25)
p_vec["mG"] <- log(1e10)
p_vec["mK"] <- log(1e10)
expect_silent(dat <- make_data(p_vec, design = designRDMSWTN, n_trials = 200))
