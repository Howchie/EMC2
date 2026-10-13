# Exactly `iter` sample iterations from a checkpoint under the between-block
# rule (max_tries = 1), rebuilding every 100; prints each block's plan record.
args <- commandArgs(TRUE)
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
cat("EMC2 from", find.package("EMC2"), "\n")
ns <- asNamespace("EMC2")
fx <- args[1]; iter <- as.integer(args[2]); out <- args[3]
emc <- readRDS(fx)
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(emc), quiet = TRUE)
n0 <- emc[[1]]$sample_transitions
if (is.null(n0)) n0 <- sum(emc[[1]]$samples$stage == "sample")
for (j in seq_along(emc)) {
  emc[[j]]$production_options$emc2.sample_rebuilds <- as.integer(seq(n0 + 100L, n0 + iter + 100L, by = 100L))
  emc[[j]]$sample_proposals_next <- n0 + 100L
  emc[[j]]$ensemble_stats$block <- NULL
  emc[[j]]$ensemble_plan_history <- NULL
  emc[[j]]$ensemble_plan <- NULL
}
set.seed(1)
idx0 <- emc[[1]]$samples$idx
t0 <- proc.time()
emc <- suppressWarnings(suppressMessages(fit(emc, stage = "sample", iter = iter, max_tries = 1,
                                             cores_per_chain = 3, cores_for_chains = 3, verbose = FALSE)))
pt <- proc.time() - t0
emc <- ns$restore_duplicates(emc)
cat(sprintf("iterations %d, wall %.0f s (%.2f s/iter), tree cpu %.0f s\n", emc[[1]]$samples$idx - idx0,
            pt[["elapsed"]], pt[["elapsed"]] / max(1, emc[[1]]$samples$idx - idx0),
            sum(pt[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE)))
tab <- do.call(rbind, lapply(seq_along(emc), function(j) do.call(rbind, lapply(emc[[j]]$ensemble_plan_history, function(r)
  data.frame(chain = j, n = r$n, M = r$M, R = r$R, C = r$C, c_round = r$c_round, c_member = r$c_member,
             G = r$G, J_tot = r$J_tot, T_min = if (is.null(r$T_min)) NA else r$T_min,
             worst = if (is.null(r$worst)) NA else r$worst, eff = r$eff, pool = r$pool, rounds = r$rounds)))))
print(tab, digits = 3, row.names = FALSE)
print(ns$ensemble_diagnostics(emc))
saveRDS(list(tab = tab, wall = pt, iters = emc[[1]]$samples$idx - idx0), out)
