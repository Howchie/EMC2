repo <- "/data/work/EMC2_dev_oo"
pkgload::load_all(repo, quiet = TRUE)
load(file.path(repo, "plans/convergence/burn_pool_200.RData"))
emc <- run_emc(emc, stage = "burn",
               stop_criteria = list(iter = 300L, mean_gd = 1.1,
                                    selection = c("alpha", "mu"),
                                    omit_mpsrf = TRUE),
               cores_for_chains = 3L, cores_per_chain = 4L,
               step_size = 100L, max_tries = 3L, trim = FALSE,
               verbose = TRUE,
               fileName = file.path(repo, "plans/convergence/burn_wide_eff_trial.RData"))
save(emc, file = file.path(repo, "plans/convergence/burn_wide_eff_trial.RData"))
