suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv('EMC_LIB')))
custom_ll <- function(pars, dadm, ...) {
  idx <- dadm$S == 'left'
  sum(log(c(dnorm(dadm$rt[idx], pars[1], exp(pars[3])),
            dnorm(dadm$rt[!idx], pars[2], exp(pars[4])))))
}
des <- design(model = custom_ll, custom_p_vector = c('muL', 'muR', 'sdL', 'sdR'))
set.seed(101)
emc <- make_emc(get_data(samples_LNR), des, type = 'standard', nuisance_non_hyper = 3:4,
                n_chains = 3)
options(emc2.fit_host = FALSE)
start <- proc.time()
out <- init_chains(emc, cores_for_chains = 2, particles = 10)
timing <- proc.time() - start
saveRDS(list(fit = out, seed = .Random.seed, timing = timing, session = sessionInfo()), Sys.getenv('OUT'))
