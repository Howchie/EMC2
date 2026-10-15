# LIB=<lib> Rscript replay_state.R <label> <emc.rds>: replay the next allocator rebuild from a saved fit.
args <- commandArgs(TRUE); label <- args[1]
suppressMessages(library(EMC2, lib.loc = Sys.getenv("LIB"))); ns <- asNamespace("EMC2")
options(emc2.alloc_transport = FALSE)
emc <- readRDS(args[2])
for (j in 2:length(emc)) for (nm in setdiff(names(emc[[1]]), names(emc[[j]]))) emc[[j]][[nm]] <- emc[[1]][[nm]]
emc <- ns$restore_custom_kernel_pointers(emc, quiet = TRUE); class(emc) <- "emc"
cat(sprintf("== %s: idx %d, window n %s\n", label, emc[[1]]$samples$idx, paste(sapply(emc, function(ch) ch$alloc$window$n), collapse = "/")))
base <- local({ set.seed(1, kind = "L'Ecuyer-CMRG"); get(".Random.seed", envir = globalenv()) })
for (j in seq_along(emc)) { base <- parallel::nextRNGStream(base); emc[[j]]$rng <- ns$.emc_chain_rng(base, emc[[j]]$n_subjects) }
h0 <- lengths(lapply(emc, function(ch) ch$alloc$history))
out <- ns$.emc_alloc_rebuild(emc, "sample", 27, 100)
for (j in seq_along(out)) {
  h <- out[[j]]$alloc$history; new <- h[(h0[j] + 1):length(h)]
  for (k in seq_along(new)) { r <- new[[k]]
    cat(sprintf(" ch%d call%d choice %-4s obj %6.1f cost %5.2f tau %5.1f tau0 %6.1f | q %.2f hmc %.2f | rounds %.2f | options %s\n",
      j, k, r$choice, r$objective, r$cost, r$tau, r$tau0, mean(r$particle_rate), mean(r$rate), sum(r$rounds),
      paste(names(r$options), round(r$options, 1), sep = "=", collapse = " "))) }
}
saveRDS(lapply(out, function(ch) ch$alloc$history[(h0[1] + 1):length(ch$alloc$history)]), file.path(Sys.getenv("SP"), paste0("replaystate_", label, ".rds")))
