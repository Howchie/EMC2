# Stop-success integration fixture for the configurable-bound EXG normalization
# cache. Each call uses two distinct go parameter tuples and a nonzero lower
# bound for both go and stop distributions.
suppressPackageStartupMessages(
  library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
ns <- asNamespace("EMC2")
texg <- get("ss_texg_stop_success_value", ns)
rdex <- get("ss_rdex_stop_success_value", ns)

texg_pars <- rbind(
  c(.50, .05, .08, .20, .03, .05, 0, 0, .05, .04),
  c(.45, .06, .12, .20, .03, .05, 0, 0, .08, .04))
rdex_pars <- rbind(
  c(2.5, .6, .2, .15, 1, .20, .03, .05, 0, 0, .04),
  c(2.1, .7, .25, .12, 1, .20, .03, .05, 0, 0, .04))
SSD <- .15
n_iter <- as.integer(Sys.getenv("N_ITER", "500"))
n_rep <- as.integer(Sys.getenv("N_REP", "9"))

call_texg <- function() texg(SSD, texg_pars, method = "gl", n_nodes = 128L)
call_rdex <- function() rdex(SSD, rdex_pars, method = "gl", n_nodes = 128L)
value_texg <- call_texg()
value_rdex <- call_rdex()

time_calls <- function(f) {
  invisible(f())
  replicate(n_rep, {
    start <- proc.time()[["elapsed"]]
    for (i in seq_len(n_iter)) f()
    proc.time()[["elapsed"]] - start
  })
}

texg_times <- time_calls(call_texg)
rdex_times <- time_calls(call_rdex)
out <- Sys.getenv("OUT")
if (nzchar(out)) saveRDS(list(
  texg = list(value = value_texg, times = texg_times),
  rdex = list(value = value_rdex, times = rdex_times),
  n_iter = n_iter, n_rep = n_rep, session = sessionInfo()), out)
cat("EXG stop-success", value_texg, "median", median(texg_times),
    "seconds /", n_iter, "calls\n")
cat("RDEX stop-success", value_rdex, "median", median(rdex_times),
    "seconds /", n_iter, "calls\n")
