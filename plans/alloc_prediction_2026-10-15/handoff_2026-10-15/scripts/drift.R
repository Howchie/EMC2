# Persistent transport along the chain-separating ridge direction: variogram slope vs one-step jump, per arm and per CPU.
SP <- Sys.getenv("SP"); A <- "/data/work/EMC2_dev_oo/plans/alloc_prediction_2026-10-15/arms"
load("/data/work/PM/NirvanaHons_Nback/Samples/samples_control_Exp1_frqk7r_blk_gL_sepA_sepL_sepV_sv_t0T_mt.RData")
p <- emc[[1]]$n_pars; pn <- emc[[1]]$par_names; idx <- 952:1051
mu <- lapply(emc, function(e) t(e$samples$theta_mu[, idx]))
wsd <- sqrt(colMeans(t(sapply(mu, function(x) apply(x, 2, var)))))
Z <- sweep(t(sapply(mu, colMeans)), 2, colMeans(t(sapply(mu, colMeans)))) %*% diag(1 / wsd)
u <- svd(Z)$v[, 1] / wsd                                    # leading chain-separating direction, raw mu units
saveRDS(u, file.path(SP, "ridge_u.rds")); cat("direction loadings (scaled):", paste(sprintf("%s %+.2f", pn, svd(Z)$v[, 1])[order(-abs(svd(Z)$v[, 1]))][1:6], collapse = ", "), "\n")
cat(sprintf("checkpoint window: between-chain sd along it %.2f within-sd units\n", sd(Z %*% svd(Z)$v[, 1])))
vg <- function(xs, K = 40) sapply(1:K, function(k) mean(sapply(xs, function(x) mean((x[-(1:k)] - x[1:(length(x) - k)])^2))))
fs <- c(old = file.path(A, "old_1.rds"), new = file.path(A, "new_1.rds"), slow = file.path(A, "slow_1.rds"),
        slow2 = file.path(A, "slow2_1.rds"), point = file.path(SP, "arm_point_1.rds"), kappa = file.path(SP, "arm_kappa_1.rds"))
ck <- lapply(emc, function(e) drop(t(e$samples$theta_mu[, 152:1051]) %*% u))
out <- NULL
for (n in names(fs)) { r <- readRDS(fs[n]); xs <- lapply(r$draws, function(g) drop(g[, 1:p] %*% u))
  g <- vg(xs); slope <- coef(lm(g[10:40] ~ I(10:40)))[2]; bounce <- coef(lm(g[10:40] ~ I(10:40)))[1]
  fin <- sapply(xs, function(x) mean(tail(x, 100))); start <- sapply(xs, function(x) mean(head(x, 50)))
  out <- rbind(out, data.frame(arm = n, cpu_it = r$cpu / 400, jump1 = g[1], bounce = bounce, drift_per_it = slope,
    drift_per_kcpu = 1000 * slope / (r$cpu / 400), spread_start = sd(start), spread_end = sd(fin))) }
print(format(out, digits = 3), row.names = FALSE)
cat("(units: squared within-sd along the ridge direction; spread = sd of the 3 chain means)\n")
