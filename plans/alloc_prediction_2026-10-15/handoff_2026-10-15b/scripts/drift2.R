# Ridge movement per arm: variogram slope (persistent transport) along the chain-separating direction u, gap between chain means, per CPU.
A <- "/data/work/EMC2_dev_oo/plans/alloc_prediction_2026-10-15"
load("/data/work/PM/NirvanaHons_Nback/Samples/samples_control_Exp1_frqk7r_blk_gL_sepA_sepL_sepV_sv_t0T_mt.RData")
p <- emc[[1]]$n_pars; pn <- emc[[1]]$par_names; idx <- 952:1051
mu <- lapply(emc, function(e) t(e$samples$theta_mu[, idx]))
wsd <- sqrt(colMeans(t(sapply(mu, function(x) apply(x, 2, var)))))
Z <- sweep(t(sapply(mu, colMeans)), 2, colMeans(t(sapply(mu, colMeans)))) %*% diag(1 / wsd)
sv <- svd(Z); u <- sv$v[, 1] / wsd
cat("ridge direction (scaled loadings):", paste(sprintf("%s %+.2f", pn, sv$v[, 1])[order(-abs(sv$v[, 1]))][1:6], collapse = ", "), "\n")
cat(sprintf("checkpoint: chain means along u (within-sd units): %s\n", paste(round(Z %*% sv$v[, 1], 2), collapse = " ")))
vg <- function(xs, K = 60) sapply(1:K, function(k) mean(sapply(xs, function(x) mean((x[-(1:k)] - x[1:(length(x) - k)])^2))))
files <- c(old_1 = file.path(A, "arms/old_1.rds"), new_1 = file.path(A, "arms/new_1.rds"),
  new_2 = file.path(A, "handoff_2026-10-15/arms/pair_new_2.rds"), new_3 = file.path(A, "handoff_2026-10-15/arms/pair_new_3.rds"),
  point_1 = file.path(A, "handoff_2026-10-15/arms/arm_point_1.rds"), point_2 = file.path(A, "handoff_2026-10-15/arms/pair_point_2.rds"),
  point_3 = file.path(A, "handoff_2026-10-15/arms/pair_point_3.rds"), kappa_1 = file.path(A, "handoff_2026-10-15/arms/arm_kappa_1.rds"),
  transport_1 = file.path(A, "eval_transport/arm_transport_1.rds"), slow_1 = file.path(A, "arms/slow_1.rds"), slow2_1 = file.path(A, "arms/slow2_1.rds"),
  path_2 = file.path(A, "handoff_2026-10-15/arms/arm_path_2.rds"))
out <- NULL
for (n in names(files)) { r <- readRDS(files[n]); xs <- lapply(r$draws, function(g) drop(g[, 1:p] %*% u))
  g <- vg(xs); fit <- coef(lm(g[10:60] ~ I(10:60))); slope <- fit[2]; bounce <- fit[1]
  m <- sapply(xs, mean); fin <- sapply(xs, function(x) mean(tail(x, 100))); start <- sapply(xs, function(x) mean(head(x, 100)))
  # per-chain net displacement over the run (end-start block means), in within-sd units
  disp <- fin - start
  rounds <- sum(r$moves["shear", ]); hmc <- sum(r$moves["hmc", ])
  out <- rbind(out, data.frame(arm = n, cpu = round(r$cpu), rounds = rounds, hmc = hmc, jump1 = g[1], bounce = bounce, slope = slope,
    slope_per_kcpu = 1000 * slope / (r$cpu / 400), spread_start = sd(start), spread_end = sd(fin),
    disp = paste(sprintf("%+.2f", disp), collapse = " "), disp_rms = sqrt(mean(disp^2)))) }
print(format(out, digits = 3), row.names = FALSE)
cat("(squared within-sd units along u; slope = persistent transport per iteration; disp = chain mean shift first 100 -> last 100)\n")
