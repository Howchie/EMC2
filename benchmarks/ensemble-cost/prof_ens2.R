# Per-iteration cost of the ensemble update and its parts on a saved sample-stage state.
#   EMC_LIB=<lib> EMC_IN=<fit.rds> ITER=200 Rscript prof_ens2.R [out.rds]
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
cat("EMC2 from", find.package("EMC2"), "\n")
ns <- asNamespace("EMC2")
emc <- readRDS(Sys.getenv("EMC_IN"))
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(emc), quiet = TRUE)
chain <- as.integer(Sys.getenv("CHAIN", "1"))
emc <- emc[chain]; class(emc) <- "emc"
iter <- as.integer(Sys.getenv("ITER", "200"))
s <- emc[[1]]
cat(sprintf("type %s p %d subjects %d pool %s rounds %s\n", s$type, s$n_pars, s$n_subjects,
            s$production_options$emc2.ensemble_pool, s$production_options$emc2.ensemble_rounds))
fns <- c(".emc_ensemble_move", ".emc_ensemble_draw_pools", ".emc_ensemble_draw_q",
         ".emc_ensemble_log_q", ".emc_ensemble_log_q_m2", ".emc_ensemble_pool_ll",
         ".emc_ensemble_store", ".emc_ensemble_subject_prior", ".emc_ensemble_log_g",
         ".emc_ensemble_select_all", ".emc_ensemble_refresh_scores",
         "gibbs_step", "gibbs_step_standard", "gibbs_step_diag_gamma", "riwish", "rwish",
         "fill_samples", "fill_samples_standard", "fill_samples_diag_gamma")
tm <- new.env()
for (f in fns) {
  if (!exists(f, envir = ns, inherits = FALSE)) next
  g <- get(f, envir = ns)
  tm[[f]] <- c(t = 0, n = 0)
  e <- new.env(parent = environment(g))
  e$f <- f; e$g <- g; e$tm <- tm
  w <- function(...) {
    t0 <- proc.time()[["elapsed"]]
    on.exit({ v <- tm[[f]]; v[["t"]] <- v[["t"]] + proc.time()[["elapsed"]] - t0
              v[["n"]] <- v[["n"]] + 1; tm[[f]] <- v })
    g(...)
  }
  environment(w) <- e
  assignInNamespace(f, w, ns)
}
set.seed(1)
prof_file <- Sys.getenv("PROF_FILE", "")
if (nzchar(prof_file)) Rprof(prof_file, interval = 0.002)
t0 <- proc.time()[["elapsed"]]
emc <- suppressMessages(fit(emc, stage = "sample", iter = iter, cores_per_chain = as.integer(Sys.getenv("CORES", "9")),
                            cores_for_chains = 1, verbose = FALSE))
wall <- proc.time()[["elapsed"]] - t0
if (nzchar(prof_file)) Rprof(NULL)
n_it <- tm[[".emc_ensemble_move"]][["n"]]
if (n_it < 1) n_it <- iter
cat(sprintf("sample stage: %d iterations, wall %.1f s, %.1f ms/iter\n", n_it, wall, 1000 * wall / n_it))
tab <- do.call(rbind, lapply(fns[fns %in% ls(tm, all.names = TRUE)], function(f) {
  v <- tm[[f]]
  data.frame(fn = f, calls = v[["n"]], calls_per_iter = v[["n"]] / n_it,
             ms_per_iter = 1000 * v[["t"]] / n_it, ms_per_call = if (v[["n"]] > 0) 1000 * v[["t"]] / v[["n"]] else NA)
}))
print(tab, row.names = FALSE, digits = 3)
st <- emc[[1]]$ensemble_stats
if (!is.null(st)) cat(sprintf("ensemble stats: sweeps %d moved %.3f wess %.3f\n", st$sweeps, st$moved / st$sweeps, st$wess / st$sweeps))
if (nzchar(prof_file)) { sm <- summaryRprof(prof_file); cat("---- by.self\n"); print(head(sm$by.self, 30)) }
out <- commandArgs(TRUE)
if (length(out)) saveRDS(list(tab = tab, wall = wall, iter = n_it, type = s$type, p = s$n_pars,
                              n = s$n_subjects, lib = find.package("EMC2"), stats = st), out[1])
