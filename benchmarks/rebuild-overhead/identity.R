# Compare baseline and candidate outputs from separate runs with a saved state.
# HOLD_ALLOCATION=1 fixes move counts, as in the sampler's invariance tests,
# so a wall-clock measurement cannot change the compared transition.
lib <- Sys.getenv("EMC_LIB")
suppressPackageStartupMessages(library(EMC2, lib.loc = lib))
cat("EMC2 from", find.package("EMC2"), "\n")
ns <- asNamespace("EMC2")
if (Sys.getenv("HOLD_ALLOCATION") == "1") {
  assignInNamespace(".emc_alloc_solve_units", function(u) {
    hmc <- vapply(u$units, function(unit) identical(unit$kind, "subject") &&
      !is.null(unit$h), logical(1))
    x <- vapply(u$units, function(unit) if (identical(unit$kind, "shear")) 1 else 0,
      numeric(1))
    x[hmc] <- 0.5
    list(x = x, q = ifelse(hmc, 0.5, 1), joint = if (!is.null(u$joint)) 1,
      rho = NA_real_, tau = NA_real_, cost = u$C_rest,
      objective = NA_real_, objective0 = NA_real_, tau0 = NA_real_)
  }, ns)
}
e <- readRDS(Sys.getenv("EMC_IN"))
set.seed(713)
stage <- Sys.getenv("STAGE", "burn")
n <- as.integer(Sys.getenv("ITER", "4"))
x <- run_emc(e, stage, list(iter = n, selection = c("mu", "sigma2", "alpha")),
  step_size = max(1L, n %/% 2L), max_tries = 2L,
  cores_for_chains = 3L, cores_per_chain = 3L, trim = FALSE, verbose = TRUE)
x <- ns$restore_duplicates(x)
keep <- lapply(x, function(ch) {
  s <- ch$samples
  list(samples = s[c("alpha", "theta_mu", "theta_var", "subj_ll", "idx", "stage")],
    rng = ch$rng, chain_mu = ch$chains_mu, chain_var = ch$chains_var,
    eff_mu = ch$eff_mu, eff_var = ch$eff_var,
    shear_rounds = ch$shear$plan$rounds, shear_steps = ch$shear$plan$log_step,
    collapsed_eps = ch$collapsed$plan$eps, collapsed_rounds = ch$collapsed$plan$rounds,
    particle_rates = ch$alloc$particle$rate, hmc_rates = ch$subject_hmc$rates)
})
saveRDS(keep, Sys.getenv("IDENTITY_OUT"))
