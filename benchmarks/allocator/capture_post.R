# The allocator's exact input at a rebuild: add_proposals() on a saved
# pre-rebuild state with the rebuild captured instead of run.
#   EMC_LIB=<lib> Rscript benchmarks/allocator/capture_post.R <pre_k.rds> <post_k.rds> [n_cores]
args <- commandArgs(TRUE)
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
cat("EMC2 from", find.package("EMC2"), "\n")
ns <- asNamespace("EMC2")
pre <- readRDS(args[1])
emc <- ns$restore_custom_kernel_pointers(pre$emc, quiet = TRUE)
for (j in seq_along(emc)) emc[[j]]$ws1_trace <- NULL
n_cores <- if (length(args) >= 3) as.integer(args[3]) else pre$n_cores
captured <- NULL
cap <- function(emc, stage, n_cores, B) { captured <<- list(emc = emc, stage = stage, n_cores = n_cores, B = B); emc }
environment(cap) <- environment()
orig <- ns$.emc_alloc_rebuild
utils::assignInNamespace(".emc_alloc_rebuild", cap, ns)
if (!is.null(pre$seed)) assign(".Random.seed", pre$seed, envir = globalenv())
invisible(ns$add_proposals(emc, pre$stage, n_cores, pre$n_blocks, pre$block_iter))
captured$seed <- get(".Random.seed", envir = globalenv())
saveRDS(captured, args[2], compress = FALSE)
cat("saved; B =", captured$B, "\n")
