pkgload::load_all(".", quiet = TRUE, compile = FALSE)
options(emc2.failure_policy = "strict")
fit_dir <- "/tmp/collapsed_transport_followup"
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
chain <- readRDS(file.path(fit_dir, "safe_fit_ess_start.rds"))$start[[2L]]
par_name <- "v_LoadHigh"
coordinate <- match(par_name, chain$par_names)
stopifnot(!is.na(coordinate))

ctx <- EMC2:::.emc_collapsed_context(chain)
idx <- chain$samples$idx
cur <- EMC2:::.emc_shear_state(chain$samples$theta_mu[, idx],
  as.matrix(chain$samples$theta_var[, , idx]), ctx$fits, chain$prior, ctx$spec)
stopifnot(!is.null(cur), is.finite(cur$lp))
phi0 <- EMC2:::.emc_collapsed_to_phi(cur$mu, cur$V, ctx$chart)
base_plan <- chain$collapsed$plan
base_C <- base_plan$C
base <- list(cur = cur, alpha = matrix(chain$samples$alpha[, , idx],
              chain$n_pars), ll = chain$samples$subj_ll[, idx], rounds = 0L,
              wpool_part = NULL)

n_pilot <- 64L
factors <- c(1, 2, 4, 8)
posterior_var <- stats::var(chain$samples$theta_mu[coordinate,
  EMC2:::.emc_alloc_draws(chain, 200L)])
seed_streams <- list(checkpoint = chain$rng$gibbs, independent = NULL)
set.seed(20261010)
seed_streams$independent <- get(".Random.seed", envir = globalenv())

result <- list()
for (replicate in seq_along(seed_streams)) {
  seed0 <- seed_streams[[replicate]]
  for (factor in factors) {
    plan <- base_plan
    plan$C[coordinate, ] <- plan$C[coordinate, ] * factor
    C_candidate <- EMC2:::.emc_collapsed_preconditioner(plan$C, plan$scale,
      plan$t0_ratio, plan$t0_idx, plan$t0_hyper_ratio, plan$t0_mean_idx)
    assign(".Random.seed", seed0, envir = globalenv())
    plan$eps <- EMC2:::.emc_collapsed_tune_eps(phi0, C_candidate, plan$eps,
      ctx, plan$length)
    seed_pilots <- get(".Random.seed", envir = globalenv())
    acceptance <- movement <- movement_c <- numeric(n_pilot)
    elapsed <- system.time(for (i in seq_len(n_pilot)) {
      assign(".Random.seed", seed_pilots, envir = globalenv())
      x <- base
      x$wpool <- NULL
      x$moved <- FALSE
      x <- EMC2:::.emc_collapsed_round(x, plan, chain, ctx)
      seed_pilots <- get(".Random.seed", envir = globalenv())
      acceptance[i] <- x$acc
      movement[i] <- x$acc * x$jump[coordinate]^2
      normalized <- forwardsolve(base_C, x$phi_jump)
      movement_c[i] <- x$acc * sum(normalized^2)
    })[["elapsed"]]
    result[[length(result) + 1L]] <- data.frame(
      parameter = par_name, replicate = names(seed_streams)[replicate],
      row_factor = factor, eps = plan$eps,
      acceptance = mean(acceptance), target_esjd = mean(movement),
      target_esjd_per_posterior_variance = mean(movement) / posterior_var,
      total_C_metric_esjd = mean(movement_c), seconds = elapsed,
      proposals = n_pilot, proposals_per_second = n_pilot / elapsed)
    cat("replicate", names(seed_streams)[replicate], "factor", factor,
        "eps", plan$eps, "acceptance", mean(acceptance),
        "target ESJD", mean(movement), "normalized total ESJD",
        mean(movement_c), "seconds", elapsed, "\n")
    flush.console()
  }
}
result <- do.call(rbind, result)
write.csv(result,
  file.path(artifact_dir, "directional_scale_v_LoadHigh_extended.csv"),
  row.names = FALSE)
saveRDS(result,
  file.path(artifact_dir, "directional_scale_v_LoadHigh_extended.rds"))
