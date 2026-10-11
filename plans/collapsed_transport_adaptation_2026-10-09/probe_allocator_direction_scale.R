pkgload::load_all(".", quiet = TRUE, compile = FALSE)

fit_dir <- "/tmp/collapsed_transport_followup"
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
probe <- readRDS(file.path(fit_dir, "capped_allocator_probe.rds"))
chain <- probe$start[[2L]]
options(emc2.subject_hmc = FALSE)

p <- chain$n_pars
plan <- chain$collapsed$plan
ctx <- EMC2:::.emc_collapsed_context(chain)
idx <- chain$samples$idx
cur <- EMC2:::.emc_shear_state(chain$samples$theta_mu[, idx],
  as.matrix(chain$samples$theta_var[, , idx]), ctx$fits, chain$prior, ctx$spec)
stopifnot(!is.null(cur), is.finite(cur$lp))
phi0 <- EMC2:::.emc_collapsed_to_phi(cur$mu, cur$V, ctx$chart)
alpha <- matrix(chain$samples$alpha[, , idx], p)
base <- list(cur = cur, alpha = alpha, ll = chain$samples$subj_ll[, idx],
             rounds = 0L, wpool = NULL, wpool_part = NULL)

t0_chart <- EMC2:::.emc_collapsed_t0_chart_idx(chain$par_names, ctx$chart)
t0_mean <- EMC2:::.emc_collapsed_t0_mean_chart_idx(chain$par_names)
C_no_raw_target <- EMC2:::.emc_collapsed_preconditioner(
  plan$C, plan$scale, plan$t0_ratio, t0_chart, plan$t0_hyper_ratio,
  t0_mean, mean_idx = integer())
S <- C_no_raw_target %*% t(C_no_raw_target)

# Allocator group coordinates are D^{-1} mu, so its slowest mean direction
# in raw parameter space is the corresponding column of D. Keep t0 group
# means fixed in this probe, matching the separate support-safe scale policy.
D <- chain$shear$plan$dirs
rho <- vapply(probe$captured[[2L]]$units$units[seq_len(p)], `[[`, 0, "rho")
k <- which.max(rho)
d <- D[, k]
d[which(startsWith(chain$par_names, "t0"))] <- 0
stopifnot(sum(d^2) > 0)
u <- d / sqrt(sum(d^2))

posterior_draws <- EMC2:::.emc_alloc_draws(chain, 200L)
Dinv <- solve(D)
coordinate_draws <- drop(Dinv[k, , drop = FALSE] %*%
  chain$samples$theta_mu[, posterior_draws, drop = FALSE])
target_var <- stats::var(coordinate_draws)

set.seed(20261010L)
trial_seeds <- replicate(24L, sample.int(.Machine$integer.max, 1L))
base_rng <- chain$rng$gibbs

make_candidate <- function(label, direction_factor = NULL, legacy = FALSE) {
  if (legacy) {
    C <- EMC2:::.emc_collapsed_preconditioner(plan$C, plan$scale,
      plan$t0_ratio, t0_chart, plan$t0_hyper_ratio, t0_mean,
      plan$mean_idx, plan$mean_ratio)
  } else if (is.null(direction_factor)) {
    C <- C_no_raw_target
  } else {
    T <- diag(ctx$chart$dim)
    T[seq_len(p), seq_len(p)] <- diag(p) +
      (direction_factor - 1) * tcrossprod(u)
    Sigma <- T %*% S %*% t(T)
    C <- t(chol((Sigma + t(Sigma)) / 2))
  }
  candidate <- plan
  candidate$C <- C
  candidate$scale <- 1
  candidate$t0_ratio <- 1
  candidate$t0_hyper_ratio <- 1
  candidate$mean_idx <- integer()
  candidate$mean_ratio <- 1
  assign(".Random.seed", base_rng, envir = globalenv())
  candidate$eps <- EMC2:::.emc_collapsed_tune_eps(phi0, C, plan$eps,
    ctx, plan$length)
  started <- proc.time()[["elapsed"]]
  acceptance <- jump2 <- numeric(length(trial_seeds))
  for (i in seq_along(trial_seeds)) {
    set.seed(trial_seeds[i])
    x <- base
    x$moved <- FALSE
    x <- EMC2:::.emc_collapsed_round(x, candidate, chain, ctx)
    acceptance[i] <- x$acc
    jump <- sum(Dinv[k, ] * x$jump[seq_len(p)])
    jump2[i] <- x$acc * jump^2
  }
  elapsed <- proc.time()[["elapsed"]] - started
  data.frame(candidate = label, eps = candidate$eps,
    acceptance = mean(acceptance), target_esjd = mean(jump2),
    target_esjd_per_variance = mean(jump2) / target_var,
    seconds = elapsed, proposals = length(trial_seeds))
}

results <- list(make_candidate("current_raw_row", legacy = TRUE),
                make_candidate("no_extra_mean_scale"))
for (f in c(2, 4, 8, 16))
  results[[length(results) + 1L]] <- make_candidate(
    paste0("shear_direction_", f), direction_factor = f)
results <- do.call(rbind, results)
t0_rows <- which(startsWith(chain$par_names, "t0"))
cat("selected shear coordinate", k, "rho", rho[k],
    "t0 direction share", sum(D[t0_rows, k]^2) / sum(D[, k]^2), "\n")
print(results, row.names = FALSE)
write.csv(results, file.path(artifact_dir,
  "allocator_direction_scale_probe.csv"), row.names = FALSE)
saveRDS(list(direction = u, shear_coordinate = k, shear_rho = rho[k],
  t0_direction_share = sum(D[t0_rows, k]^2) / sum(D[, k]^2), results = results),
  file.path(artifact_dir, "allocator_direction_scale_probe.rds"))
