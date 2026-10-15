suppressMessages(library(EMC2, lib.loc = paste0(Sys.getenv("SP"), "/scr/libnew")))
load("Samples/samples_control_Exp1_frqk7r_blk_gL_sepA_sepL_sepV_sv_t0T_mt.RData")
p <- emc[[1]]$n_pars
for (f in c("old_1", "new_1")) {
  r <- readRDS(file.path(Sys.getenv("SP"), "scr", paste0(f, ".rds")))
  for (j in 1:3) {
    D <- solve(emc[[j]]$shear$plan$dirs)
    X <- lapply(1:3, function(i) { g <- r$draws[[i]][1:100, ]; cbind(g[, 1:p] %*% t(D), g[, p + 1:p]) })
    pv <- apply(do.call(rbind, X), 2, var); xj <- X[[j]]
    esjd <- colMeans(diff(xj)^2)
    pooled <- 1 - esjd / (2 * pv); within <- 1 - esjd / (2 * apply(xj, 2, var))
    m <- esjd / (2 * pv); lo <- pmax(0, 1 - m * (1 + 2 * sqrt(2 / 99)))
    pred <- unlist(r$plans[[j]][[2]]$pred_group)
    cat(sprintf("%s ch%d  mu: pred-pooled %+.2f pred-within %+.2f pooledvar/within %.2f | lsd: pred-pooled %+.2f | UCB shift %+.2f\n", f, j,
      median(pred[1:p] - pooled[1:p]), median(pred[1:p] - within[1:p]), median(pv[1:p] / apply(xj, 2, var)[1:p]),
      median(pred[p + 1:p] - pooled[p + 1:p]), median(lo - pooled)))
  }
}
