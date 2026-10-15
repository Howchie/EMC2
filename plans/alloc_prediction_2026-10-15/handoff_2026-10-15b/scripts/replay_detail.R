# LIB=<lib> Rscript replay_detail.R <label> <emc.rds>: replay the next allocator rebuild from a saved fit and
# show, per chain, the options, the plan, and the units that bind the target.
args <- commandArgs(TRUE); label <- args[1]
suppressMessages(library(EMC2, lib.loc = Sys.getenv("LIB"))); ns <- asNamespace("EMC2")
options(emc2.alloc_transport = FALSE)
if (grepl("RData$", args[2])) load(args[2]) else emc <- readRDS(args[2])
for (j in 2:length(emc)) for (nm in setdiff(names(emc[[1]]), names(emc[[j]]))) emc[[j]][[nm]] <- emc[[1]][[nm]]
emc <- ns$restore_custom_kernel_pointers(emc, quiet = TRUE); class(emc) <- "emc"
p <- emc[[1]]$n_pars; pn <- emc[[1]]$par_names; nsub <- emc[[1]]$n_subjects
u <- tryCatch(readRDS("/data/work/EMC2_dev_oo/plans/alloc_prediction_2026-10-15/handoff_2026-10-15/scripts/ridge_u.rds"), error = function(e) NULL)
base <- local({ set.seed(1, kind = "L'Ecuyer-CMRG"); get(".Random.seed", envir = globalenv()) })
for (j in seq_along(emc)) { base <- parallel::nextRNGStream(base); emc[[j]]$rng <- ns$.emc_chain_rng(base, emc[[j]]$n_subjects) }
h0 <- lengths(lapply(emc, function(ch) ch$alloc$history))
stage <- emc[[1]]$samples$stage[emc[[1]]$samples$idx]
t0 <- proc.time()[["elapsed"]]
out <- ns$.emc_alloc_rebuild(emc, stage, 27, 100)
cat(sprintf("== %s: stage %s idx %d, window n %s, rebuild %.1f s\n", label, stage, emc[[1]]$samples$idx,
            paste(sapply(emc, function(ch) ch$alloc$window$n), collapse = "/"), proc.time()[["elapsed"]] - t0))
tau1 <- function(r) { r <- pmin(pmax(r, 0), 1 - 1e-9); (1 + r) / (1 - r) }
tot_cost <- 0
for (j in seq_along(out)) {
  h <- out[[j]]$alloc$history; r <- h[[length(h)]]
  dirs <- emc[[j]]$shear$plan$dirs
  k_ridge <- if (is.null(u) || is.null(dirs)) NA else which.max(abs(drop(t(dirs) %*% u)) / sqrt(colSums(dirs^2)))
  gname <- c(paste0("dir", seq_len(p)), paste0("lsd.", pn))
  kap <- if (is.null(r$group_kappa)) rep(1, 2 * p) else r$group_kappa
  eff <- function(rho, k) { if (is.null(k)) rho else { t <- k * tau1(rho); (t - 1) / (t + 1) } }
  g_eff0 <- eff(r$group_rho, kap); g_eff <- eff(r$pred_group, kap)
  s_pred <- apply(r$pred_subject, 1, max); s0 <- apply(r$pred0_subject, 1, max)
  cat(sprintf("ch%d choice %-4s options [%s] | tau %.1f (tau0 %.1f) cost %.2f (C_rest %.2f part %.2f rounds %.2fx%.2f hmc %.2f coll %.2f) obj %.1f\n",
    j, r$choice, paste(names(r$options), round(r$options, 1), if (is.null(r$option_cost)) "" else paste0("/", round(r$option_cost, 1)), sep = "", collapse = " "), r$tau, r$tau0, r$cost, r$C_rest,
    sum(r$particle_shared * r$c_part, na.rm = TRUE), sum(r$rounds), r$c_round, sum(r$rate * r$c_hmc, na.rm = TRUE) + (if (any(r$rate > 0)) r$hmc_fixed else 0),
    r$collapsed_rounds * r$collapsed_cost, r$objective))
  cat(sprintf("     ridge dir %s: rho0 %.3f kappa %.1f tau0 %.0f -> plan rounds %.2f (g %.3f a %.3f) pred tau %.0f | q mean %.2f (q<1 %d) hmc mean %.2f (>0 %d)\n",
    k_ridge, r$group_rho[k_ridge], kap[k_ridge], kap[k_ridge] * tau1(r$group_rho[k_ridge]), r$rounds[k_ridge], r$group_g[k_ridge], r$group_a[k_ridge],
    kap[k_ridge] * tau1(r$pred_group[k_ridge]), mean(r$particle_shared), sum(r$particle_shared < 1), mean(r$rate), sum(r$rate > 0)))
  # slowest predicted units under the plan, and before it
  allp <- c(setNames(kap * tau1(r$pred_group), gname), setNames(tau1(s_pred), paste0("s", seq_len(nsub))))
  all0 <- c(setNames(kap * tau1(r$group_rho), gname), setNames(tau1(s0), paste0("s", seq_len(nsub))))
  o <- order(-allp)[1:8]
  cat(sprintf("     slowest under plan: %s\n     slowest before:     %s\n",
    paste(sprintf("%s %.0f", names(allp)[o], allp[o]), collapse = " "),
    paste(sprintf("%s %.0f", names(all0)[order(-all0)[1:8]], all0[order(-all0)[1:8]]), collapse = " ")))
  cat(sprintf("     kappa: median %.1f, mu dirs %s | lsd %s\n", median(kap), paste(round(kap[1:p], 1), collapse = " "), paste(round(kap[p + 1:p], 1), collapse = " ")))
  tot_cost <- tot_cost + r$cost
}
cat(sprintf("total planned cost per iteration over chains: %.2f\n", tot_cost))
saveRDS(lapply(out, function(ch) ch$alloc$history[[length(ch$alloc$history)]]), file.path(Sys.getenv("SP"), paste0("replaydetail_", label, ".rds")))
