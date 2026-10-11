pkgload::load_all(".", quiet = TRUE, compile = FALSE)

artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
state <- readRDS("/tmp/lba_sampling_audit/shared_allocation_current.rds")
chain <- state$chains[[2L]]
captured <- state$captured[[2L]]$units
p <- chain$n_pars
plan <- chain$collapsed$plan
ctx <- EMC2:::.emc_collapsed_context(chain)
idx <- chain$samples$idx
cur <- EMC2:::.emc_shear_state(chain$samples$theta_mu[, idx],
  as.matrix(chain$samples$theta_var[, , idx]), ctx$fits, chain$prior, ctx$spec)
stopifnot(!is.null(cur), is.finite(cur$lp))
phi0 <- EMC2:::.emc_collapsed_to_phi(cur$mu, cur$V, ctx$chart)
alpha0 <- matrix(chain$samples$alpha[, , idx], p)
ll0 <- as.numeric(chain$samples$subj_ll[, idx])
base <- list(cur = cur, alpha = alpha0, ll = ll0, rounds = 0L,
             wpool = NULL, wpool_part = NULL)

# The allocator's first group coordinate is the current bottleneck.
k <- which.max(vapply(captured$units[seq_along(captured$joint$idx)],
                      function(u) max(u$rho), numeric(1)))
D <- chain$shear$plan$dirs
Dinv <- solve(D)
d <- D[, k]
d[startsWith(chain$par_names, "t0")] <- 0
stopifnot(sum(d^2) > 0)
u <- d / sqrt(sum(d^2))

t0_chart <- EMC2:::.emc_collapsed_t0_chart_idx(chain$par_names, ctx$chart)
t0_mean <- EMC2:::.emc_collapsed_t0_mean_chart_idx(chain$par_names)
C_base <- EMC2:::.emc_collapsed_preconditioner(plan$C, plan$scale,
  plan$t0_ratio, t0_chart, plan$t0_hyper_ratio, t0_mean,
  mean_idx = integer(), mean_ratio = 1)
C_current <- EMC2:::.emc_collapsed_preconditioner(plan$C, plan$scale,
  plan$t0_ratio, t0_chart, plan$t0_hyper_ratio, t0_mean,
  plan$mean_idx, plan$mean_ratio)

production_draws <- lapply(state$chains, function(ch) {
  which(ch$samples$stage[seq_len(ch$samples$idx)] == "sample")
})
coordinate_draws <- do.call(cbind, lapply(seq_along(state$chains), function(j) {
  ch <- state$chains[[j]]
  Dinv %*% ch$samples$theta_mu[, production_draws[[j]], drop = FALSE]
}))
coordinate_var <- apply(coordinate_draws, 1L, stats::var)
target_var <- coordinate_var[k]
logsd_draws <- do.call(cbind, lapply(seq_along(state$chains), function(j) {
  ch <- state$chains[[j]]
  vapply(production_draws[[j]], function(t)
    EMC2:::.emc_alloc_lsd(ch$samples$theta_var[, , t]), numeric(p))
}))
logsd_var <- apply(logsd_draws, 1L, stats::var)

n_probe <- 64L
set.seed(20261010L)
trial_seeds <- replicate(n_probe, sample.int(.Machine$integer.max, 1L))
base_rng <- chain$rng$gibbs

make_C <- function(factor, base = C_base) {
  T <- diag(ctx$chart$dim)
  T[seq_len(p), seq_len(p)] <- diag(p) + (factor - 1) * tcrossprod(u)
  Sigma <- T %*% base %*% t(base) %*% t(T)
  t(chol((Sigma + t(Sigma)) / 2))
}

evaluate <- function(label, C, actual_plan = FALSE) {
  candidate <- plan
  if (!actual_plan) {
    candidate$C <- C
    candidate$scale <- 1
    candidate$t0_ratio <- 1
    candidate$t0_hyper_ratio <- 1
    candidate$mean_idx <- integer()
    candidate$mean_ratio <- 1
  }
  assign(".Random.seed", base_rng, envir = globalenv())
  candidate$eps <- EMC2:::.emc_collapsed_tune_eps(phi0, C,
    if (actual_plan) plan$eps else plan$eps, ctx, plan$length)
  assign(".Random.seed", base_rng, envir = globalenv())
  acceptance <- target <- group_mu <- group_logsd <- numeric(n_probe)
  started <- proc.time()[["elapsed"]]
  for (i in seq_len(n_probe)) {
    set.seed(trial_seeds[i])
    x <- base
    x$moved <- FALSE
    x <- EMC2:::.emc_collapsed_round(x, candidate, chain, ctx)
    acceptance[i] <- x$acc
    jump_mu <- drop(Dinv %*% x$jump[seq_len(p)])
    jump_logsd <- x$jump[p + seq_len(p)]
    target[i] <- x$acc * jump_mu[k]^2 / target_var
    group_mu[i] <- x$acc * sum(jump_mu^2 / coordinate_var)
    group_logsd[i] <- x$acc * sum(jump_logsd^2 / logsd_var)
  }
  elapsed <- proc.time()[["elapsed"]] - started
  data.frame(candidate = label, eps = candidate$eps,
    exact_acceptance = mean(acceptance), target_factor =
      max(0, 1 - mean(target) / 2), target_esjd_per_variance = mean(target),
    all_group_mu_esjd = mean(group_mu), all_group_logsd_esjd = mean(group_logsd),
    seconds = elapsed, seconds_per_proposal = elapsed / n_probe)
}

results <- list(evaluate("current_raw_row", C_current, actual_plan = TRUE))
for (factor in c(1.5, 2, 3, 4, 6))
  results[[length(results) + 1L]] <- evaluate(
    paste0("add_shear_direction_x", factor), make_C(factor, C_current))
results <- do.call(rbind, results)
dominant_raw <- which.max(abs(D[, k]))
cat("chain", 2L, "slow shear coordinate", k,
    "dominant raw parameter", chain$par_names[dominant_raw],
    "baseline rho", max(captured$units[[k]]$rho),
    "t0 direction share", sum(D[startsWith(chain$par_names, "t0"), k]^2) /
      sum(D[, k]^2), "\n")
print(results, row.names = FALSE)
write.csv(results, file.path(artifact_dir,
  "allocator_direction_additive_paired.csv"), row.names = FALSE)
saveRDS(list(chain = 2L, coordinate = k, target_var = target_var,
  n_probe = n_probe, results = results), file.path(artifact_dir,
  "allocator_direction_additive_paired.rds"))
