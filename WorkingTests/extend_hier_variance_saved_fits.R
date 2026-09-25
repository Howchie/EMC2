## Continue the two saved applied fits and compare old versus new SD draws.
## Runtime is several minutes and the saved project fits are required.
## REPRO_FITS_DIR selects the source directory; REPRO_OUT_DIR selects RDS output.
suppressMessages(library(EMC2))

fits_dir <- Sys.getenv("REPRO_FITS_DIR", "/data/work/PM/NirvanaHons_Nback/Fits")
out_dir <- Sys.getenv("REPRO_OUT_DIR", "/tmp")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
files <- c(
  "Match_Control_Exp1_rdm7s_mass_noA_t0R_vs2_mt_d42.RData",
  "Match_Control_Exp1_rdm7s_mass_noA_t0R_UT0_vs2_mt_d42.RData"
)

sd_medians <- function(emc, draws) {
  by_chain <- lapply(emc, function(chain) {
    v <- chain$samples$theta_var
    vapply(draws, function(t) sqrt(diag(v[, , t])), numeric(dim(v)[1L]))
  })
  apply(do.call(cbind, by_chain), 1L, median)
}

results <- vector("list", length(files))
for (i in seq_along(files)) {
  path <- file.path(fits_dir, files[[i]])
  if (!file.exists(path)) stop("Missing saved fit: ", path)
  box <- new.env(parent = emptyenv())
  load(path, envir = box)
  emc <- box$emc
  stopifnot(length(emc) == 3L,
            dim(emc[[1L]]$samples$theta_var)[3L] == 1000L)

  ## Stored chain RNG states govern the chains. This fixes the parent RNG used
  ## while proposals are constructed, for an independently rerunnable test.
  set.seed(20260925L + i)
  emc <- run_emc(emc, stage = "sample",
    stop_criteria = list(iter = 500L, max_gd = 1.1,
      selection = c("alpha", "mu", "sigma2"), omit_mpsrf = TRUE),
    cores_for_chains = 3L, cores_per_chain = 9L,
    step_size = 100L, max_tries = 6L, verbose = TRUE)

  count <- dim(emc[[1L]]$samples$theta_var)[3L]
  result <- list(
    file = files[[i]],
    count = count,
    old = sd_medians(emc, seq_len(1000L)),
    new = sd_medians(emc, 1001L:count),
    sigma2_rhat = max(unlist(gd_summary(emc, selection = "sigma2")), na.rm = TRUE)
  )
  names(result$old) <- names(result$new) <- emc[[1L]]$par_names
  results[[i]] <- result
  saveRDS(result, file.path(out_dir, paste0("hier_variance_extend_", i, ".rds")))
  cat("Fit", i, "sample draws", count, "max sigma2 Rhat", result$sigma2_rhat,
      "v.d_match_LoadHigh old/new SD medians",
      result$old[["v.d_match_LoadHigh"]],
      result$new[["v.d_match_LoadHigh"]], "\n")
}
cat("Maximum between-fit SD-median gap: old",
    max(abs(results[[1L]]$old - results[[2L]]$old)), "new",
    max(abs(results[[1L]]$new - results[[2L]]$new)), "\n")
