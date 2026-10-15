# PM screening arms side by side: CPU, ESS per kCPU, ridge drift, lsd ESS.
suppressPackageStartupMessages(library(posterior))
SP <- Sys.getenv("SP"); A <- "/data/work/EMC2_dev_oo/plans/alloc_prediction_2026-10-15"
load("/data/work/PM/NirvanaHons_Nback/Samples/samples_control_Exp1_frqk7r_blk_gL_sepA_sepL_sepV_sv_t0T_mt.RData")
p <- emc[[1]]$n_pars; pn <- emc[[1]]$par_names
u <- readRDS(file.path(A, "handoff_2026-10-15/scripts/ridge_u.rds"))
files <- c(new_1 = file.path(A, "arms/new_1.rds"), new_2 = file.path(A, "handoff_2026-10-15/arms/pair_new_2.rds"), new_3 = file.path(A, "handoff_2026-10-15/arms/pair_new_3.rds"),
  point_1 = file.path(A, "handoff_2026-10-15/arms/arm_point_1.rds"), point_2 = file.path(A, "handoff_2026-10-15/arms/pair_point_2.rds"), point_3 = file.path(A, "handoff_2026-10-15/arms/pair_point_3.rds"),
  cand_1 = file.path(SP, "arms/cand_1.rds"), cand_2 = file.path(SP, "arms/cand_2.rds"), cand_3 = file.path(SP, "arms/cand_3.rds"))
files <- files[file.exists(files)]
vg <- function(xs, K = 60) sapply(1:K, function(k) mean(sapply(xs, function(x) mean((x[-(1:k)] - x[1:(length(x) - k)])^2))))
out <- NULL
for (n in names(files)) { r <- readRDS(files[n]); G <- r$draws
  xs <- lapply(G, function(g) drop(g[, 1:p] %*% u)); g <- vg(xs); slope <- coef(lm(g[10:60] ~ I(10:60)))[2]
  ridge <- paste0("mu.", c("v_LoadHigh", "v_LoadLow", "K", "t0_LoadHigh", "t0_LoadLow", "sv")); lsd <- grep("^lsd", names(r$ess))
  out <- rbind(out, data.frame(arm = n, cpu = round(r$cpu), wall = round(r$elapsed), hmc = sum(r$moves["hmc", ]), rounds = sum(r$moves["shear", ]),
    ridge_ess_kcpu = 1000 * min(r$ess[ridge]) / r$cpu, med_ess_kcpu = 1000 * median(r$ess) / r$cpu, lsd_min_ess = min(r$ess[lsd]), lsd_med_ess = median(r$ess[lsd]),
    max_rhat = max(r$rhat), ridge_rhat = max(r$rhat[ridge]), drift_kcpu = 1000 * slope / (r$cpu / 400))) }
print(format(out, digits = 3), row.names = FALSE)
