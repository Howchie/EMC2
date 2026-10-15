# Assemble (block b draws, plan from rebuild at end of b, block b+1 draws) for the PM burn and the arms.
suppressMessages(library(EMC2, lib.loc = "/data/work/EMC2_dev_oo/.lib"))
SP <- "SCRATCH"
PL <- "/data/work/EMC2_dev_oo/plans/alloc_prediction_2026-10-15"
load("/data/work/PM/NirvanaHons_Nback/Samples/samples_control_Exp1_frqk7r_blk_gL_sepA_sepL_sepV_sv_t0T_mt.RData")
p <- emc[[1]]$n_pars; D <- solve(emc[[1]]$shear$plan$dirs)
lsd <- function(V) 0.5 * log(pmax(apply(V, 3, diag), .Machine$double.xmin))
X <- lapply(emc, function(ch) { s <- ch$samples; cbind(t(D %*% s$theta_mu[, 1:s$idx]), t(lsd(s$theta_var[, , 1:s$idx]))) })
# the plan in force after each rebuild: the chain's last record at that iteration
plans_of <- function(h) { it <- vapply(h, `[[`, numeric(1), "iteration"); h[!duplicated(it, fromLast = TRUE)] }
pmplans <- lapply(emc, function(ch) plans_of(ch$alloc$history))
burn <- list(X = X, plans = pmplans)
arms <- lapply(c("old_1", "new_1", "slow_1", "slow2_1"), function(f) { r <- readRDS(file.path(PL, "arms", paste0(f, ".rds")))
  n0 <- 1051
  list(X = lapply(r$draws, function(g) rbind(matrix(NA, n0, 2 * p), cbind(g[, 1:p] %*% t(D), g[, p + 1:p]))),
       plans = lapply(r$plans, plans_of)) })
names(arms) <- c("old_1", "new_1", "slow_1", "slow2_1")
saveRDS(list(pm_burn = burn, arms = arms, p = p), file.path(SP, "pairs.rds"))
for (nm in names(arms)) print(sapply(arms[[nm]]$plans, function(h) sapply(h, function(x) paste(x$iteration, x$pairs[1], x$stage))))
