# LIB=<lib> Rscript replay_q.R <label> [rule]: replay one allocator rebuild at the PM checkpoint.
args <- commandArgs(TRUE); label <- args[1]
suppressMessages(library(EMC2, lib.loc = Sys.getenv("LIB"))); ns <- asNamespace("EMC2")
load("/data/work/PM/NirvanaHons_Nback/Samples/samples_control_Exp1_frqk7r_blk_gL_sepA_sepL_sepV_sv_t0T_mt.RData")
for (j in 2:length(emc)) for (nm in setdiff(names(emc[[1]]), names(emc[[j]]))) emc[[j]][[nm]] <- emc[[1]][[nm]]
emc <- ns$restore_custom_kernel_pointers(emc, quiet = TRUE); class(emc) <- "emc"
base <- local({ set.seed(1, kind = "L'Ecuyer-CMRG"); get(".Random.seed", envir = globalenv()) })
for (j in seq_along(emc)) { base <- parallel::nextRNGStream(base); emc[[j]]$rng <- ns$.emc_chain_rng(base, emc[[j]]$n_subjects) }
h0 <- lengths(lapply(emc, function(ch) ch$alloc$history))
t0 <- proc.time()[["elapsed"]]
out <- ns$.emc_alloc_rebuild(emc, "sample", 27, 100)
el <- proc.time()[["elapsed"]] - t0
cat(sprintf("== %s (rebuild %.1f s)\n", label, el))
for (j in seq_along(out)) {
  h <- out[[j]]$alloc$history; new <- h[(h0[j] + 1):length(h)]
  for (k in seq_along(new)) { r <- new[[k]]
    cat(sprintf(" ch%d call%d choice %-4s obj %7.1f cost %5.2f tau %5.1f | q mean %.2f (q<1: %3d) | hmc mean %.2f | rounds %.2f coll %.2f\n",
      j, k, r$choice, r$objective, r$cost, r$tau, mean(r$particle_rate), sum(r$particle_rate < 1), mean(r$rate),
      sum(r$rounds), r$collapsed_rounds)) }
  q <- out[[j]]$alloc$particle$rate
  cat(sprintf(" ch%d shared q mean %.2f  final hmc mean %.2f\n", j, mean(q), mean(out[[j]]$alloc$last$rate)))
}
saveRDS(lapply(out, function(ch) ch$alloc), file.path(Sys.getenv("SP"), paste0("replay_", label, ".rds")))
