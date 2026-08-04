# Extracted from test-timed-race-rfun.R:65

# test -------------------------------------------------------------------------
set.seed(456)
timed_design <- design(
    factors = list(subjects = 1, S = "stim"),
    Rlevels = c("left", "right", "time"),
    matchfun = function(d) as.character(d$S) == as.character(d$lR),
    model = LBA,
    formula = list(v ~ lR, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1)
  )
p_vec <- sampled_pars(timed_design, doMap = FALSE)
p_vec["v_lRleft"] <- 1e-6
p_vec["v_lRright"] <- 1e-6
p_vec["v_lRtime"] <- 8.0
p_vec["sv"] <- log(0.1)
p_vec["B"] <- log(1.1)
p_vec["A"] <- log(0.1)
p_vec["t0"] <- log(0.2)
p_sd <- setNames(rep(1e-8, length(p_vec)), names(p_vec))
pr <- prior(timed_design, mu_mean = p_vec, mu_sd = p_sd)
template <- data.frame(
    subjects = factor(rep(1, 40)),
    S = factor(rep("stim", 40), levels = "stim"),
    R = factor(rep(NA_character_, 40), levels = c("left", "right", "time")),
    rt = NA_real_
  )
pp <- predict(pr, data = template, n_post = 2, n_cores = 1)
