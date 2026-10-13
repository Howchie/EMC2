# Profile one allocator rebuild on a captured input (capture_post.R).
#   EMC_LIB=<lib> Rscript benchmarks/allocator/prof_rebuild.R <post.rds> [reps]
args <- commandArgs(TRUE)
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
cat("EMC2 from", find.package("EMC2"), "\n")
ns <- asNamespace("EMC2")
post <- readRDS(args[1])
options(post$options)
reps <- if (length(args) >= 2) as.integer(args[2]) else 2L
run <- function() { assign(".Random.seed", post$seed, envir = globalenv())
  ns$.emc_alloc_rebuild(post$emc, post$stage, post$n_cores, post$B) }
invisible(run())
el <- replicate(reps, system.time(run())[["elapsed"]])
cat("rebuild elapsed:", round(el, 2), "\n")
Rprof(tf <- tempfile(), interval = 0.01)
out <- run()
Rprof(NULL)
s <- summaryRprof(tf)
print(head(s$by.total[, 1:2], 30))
saveRDS(lapply(out, function(ch) ch$alloc$last), sub("\\.rds$", paste0("_plan_", basename(dirname(dirname(Sys.getenv("EMC_LIB")))), basename(Sys.getenv("EMC_LIB")), ".rds"), args[1]))
