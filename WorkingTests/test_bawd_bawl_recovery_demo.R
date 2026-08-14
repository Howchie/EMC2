## RLF Race Model Simulation Recovery Test
## Testing low alpha (1.1) and high alpha (1.8) single-subject single-condition recovery

rm(list = ls())
library(EMC2)

n_trials = 1000
label = "bawd_test"
cores_for_chains = 3
cores_per_chain = 10
matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

designBAwD <- design(
  factors = list(subjects = 1, S = c("left", "right")),
  Rlevels = c("left", "right"),
  matchfun = matchfun,
  model = BAwD(),
  formula = list(mu ~ S*lM, B ~ 1, A ~ 1, t0 ~ 1, sigma ~ 1, ell~1, k ~ 1),
  constants = c(ell=log(1),`mu_Sright`=0),
  TC = list(UC=2.5)
)

designBAwL <- design(
  factors = list(subjects = 1, S = c("left", "right")),
  Rlevels = c("left", "right"),
  matchfun = matchfun,
  model = BAwL(drift_distribution="lognormal"),
  formula = list(mu ~ S*lM, B ~ 1, A ~ 1, t0 ~ 1, sigma ~ 1, k ~ 1),
  constants = c(B=log(1),`mu_Sright`=0),
  TC = list(UC=2.5)
)

p_vector <- sampled_pars(designBAwD, doMap = FALSE)
p_vector["sigma"] <- log(.3)
p_vector["B"] <- log(.8)
p_vector["A"] <- log(.4)
p_vector["t0"] <- log(0.2)
p_vector["mu"] <- log(2)
p_vector["mu_lMTRUE"] <- log(2.5)
p_vector["mu_Sright:lMTRUE"] <- -.2
# For alpha: mapped_alpha = 1 + (2-1)*pnorm(alpha_sampled)
# So alpha_sampled = qnorm(alpha_mapped - 1)
p_vector["k"] <- log(1.2)

RNGkind("L'Ecuyer-CMRG")
set.seed(42)

cat("Generating simulated dataset...\n")
t_sim <- system.time({
  dat <- make_data(p_vector, designBAwD, n_trials = n_trials)
})
cat(sprintf("Data generation time: %.3f sec\n", t_sim["elapsed"]))

cat("Data summary:\n")
cat("Counts by response:\n")
print(table(dat$R, useNA = "ifany"))
cat("Mean RT by stimulus:\n")
print(tapply(dat$rt, dat$S, function(x) mean(x[is.finite(x)], na.rm = TRUE)))

emc <- make_emc(dat, designBAwD, type = "single")
emc2 <- make_emc(dat, designBAwL, type = "single")
cat("\nFitting BAwD model (MCMC)...\n")
t_fit <- system.time({
  fit_res <- fit(
    emc,
    fileName = "samples.RData",
    cores_for_chains = cores_for_chains,
    cores_per_chain = cores_per_chain,
    stop_criteria = list(
      sample = list(
        iter = 1000,
        max_gd = 1.10,
        max_flat_loc = 0.5,
        flat_selection = c("alpha", "subj_ll"),
        flat_p1 = 1/3,
        flat_p2 = 1/3,
        max_sample_iter = 5000
      )
    ),
    max_tries = 30
  )
})
t_fit2 <- system.time({
  fit_res2 <- fit(
    emc2,
    fileName = "samples.RData"),
    cores_for_chains = cores_for_chains,
    cores_per_chain = cores_per_chain,
    stop_criteria = list(
      sample = list(
        iter = 1000,
        max_gd = 1.10,
        max_flat_loc = 0.5,
        flat_selection = c("alpha", "subj_ll"),
        flat_p1 = 1/3,
        flat_p2 = 1/3,
        max_sample_iter = 5000
      )
    ),
    max_tries = 30
  )
})
cat(sprintf("\nFit Completed! System time for fit (%s):\n", label))
print(t_fit)
print(t_fit2)
modelList = list("BAwD" = fit_res,"BAwL" = fit_res2)
comparison = compare(modelList,cores_for_props = 4,cores_for_loo = 8)
cat("\n--- Parameter Recovery (Sampled Scale) ---\n")
rec_stats <- recovery(fit_res, true_pars = p_vector)
print(rec_stats)


predData = predict(fit_res)
predData2 = predict(fit_res2)
pdf("fit-BAwD.pdf")
plot_cdf(dat, post_predict=predData,
         functions=list(Correct=matchfun), defective_factor = "R", factors = "S")
plot_cdf(dat, post_predict=predData2,
         functions=list(Correct=matchfun), defective_factor = "R", factors = "S")
dev.off()

invisible(list(
  data = dat,
  pred = predData,
  design = designRLF,
  true_pars = p_vector,
  fit = fit_res,
  system_time = t_fit,
  mapped_summary = mapped_summary
))

