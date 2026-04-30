rm(list=ls())
source("WorkingTests/LogicalRulesContrasts.R")
source("WorkingTests/SFT_Functions.R")
library(EMC2)
library(dplyr)

designLR <- design(
  fixed_accumulator_roles = factor(c("A","n_A","B","n_B"),levels=c("A","n_A","B","n_B")),
  matchfun=matchfun,
  model=LogicalRulesLBA,constants=c('sv'=log(1)),
  factors=list(subjects=1,S=factor(rep(c("AB","AN","NB","NN"),2),levels=c("AB","AN","NB","NN")),LogicalRule="OR"),
  Rlevels=c("yes","no"),
  formula=list(v~0+mismatch+MatchY+MatchN+DT,B~0 + Yes + No,t0~1,A~1,sv~1),
  functions=funcs, UT=Inf
  )

p_vec = sampled_pars(designLR)
p_vec["v_mismatch"] = 1; p_vec["v_MatchY"] = 2.5; p_vec["v_MatchN"] = 2.8; p_vec["v_DT"] = -0.25
p_vec["B_Yes"] = log(1); p_vec["B_No"] = log(1.2)
p_vec["t0"] = log(.2); p_vec["A"] = log(.4)
mapped_pars(designLR,p_vec)

dat = make_data(p_vec,designLR, n_trials = 500) %>%
  .annotate_trials(.) %>%
  mutate(Correct=correctfun(.))
sum = dat %>%
  group_by(S) %>%
  summarise(Acc=mean(Correct), RT=mean(rt))
print(sum)

cap = .compute_capacity(dat,"OR")
cap_df = build_ct_df(cap$OR_Ct$times,list(list(fn=cap$OR_Ct$Ct,rt=dat$rt,label="Ct")))
plot(cap_df$Time,cap_df$Ct, type='l',col='blue')

## Test Likelihood
source("WorkingTests/test_likelihood_plotfuns_ah.R")
timing = system.time(profile_plot_test(dat, designLR, p_vec, n_cores = 3, range=1,
                  layout = c(3,3), use_c = TRUE,  figure_title = "LogicalRules", natural=TRUE))
emc = make_emc(dat,designLR,type="single")
emc = fit(emc,
          stop_criteria = list(
            sample = list(
              iter = 1000,
              max_gd = 1.10,
              max_flat_loc = 0.5,
              flat_selection = c("alpha", "subj_ll"),
              flat_p1 = 1/3,
              flat_p2 = 1/3,
              max_sample_iter = 5000
            ),
            cores_per_chain = 3,
            cores_for_chains = 3
          ),
          max_tries = 30)
recovery(emc, true_pars = p_vec)
post_predict <- predict(emc, n_post = 50)
plot_pars(emc, post_predict = post_predict, true_pars = p_vec)
save.image("TestLogicalRules.RData")
