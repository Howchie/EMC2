# Summarise burn_arms.R results by case, build label and stopping rule.
#   Rscript benchmarks/mode-start/summarise_burn_arms.R
out <- Sys.getenv("MODE_START_OUT", "benchmarks/mode-start/results")
f <- list.files(out, pattern = "^burnarm_.*\\.rds$", full.names = TRUE)
res <- do.call(rbind, lapply(f, readRDS))
options(width = 200)
cols <- c("case", "arm", "stop", "seed", "burn_iter", "n_pairs_rhat_1.1", "n_pairs",
          "n_subj_rhat_1.1", "n_subj_stuck", "max_rhat_alpha", "max_rhat_mu", "eps_median")
print(res[order(res$case, res$stop, res$arm, res$seed), cols], row.names = FALSE, digits = 3)
