# Chain spread and cores_per_chain scaling on a real fit:
# /data/work/PM/AnnikaHons_PM1_3/PMDC_RDMSWTN_control.R (RDMSWTN, control condition).
#
#   Rscript WorkingTests/bench_annika_real.R
#
# EMC_CPC  cores_per_chain grid   EMC_ITER  iterations per timed block

lib <- Sys.getenv("EMC_LIB", .libPaths()[1])
self <- "WorkingTests/bench_annika_real.R"
root <- "/data/work/PM/AnnikaHons_PM1_3"
setup_file <- Sys.getenv("EMC_SETUP", file.path(tempdir(), "annika_emc.rds"))
iter <- as.integer(Sys.getenv("EMC_ITER", "25"))

# FunctionsControlRDM.R builds its designs at source time and refers to `data`
# by name, so the whole setup has to happen in the global environment -- inside
# a function `data` resolves to utils::data instead.
build <- function() {
  suppressMessages({library(dplyr); library(tidyr); library(forcats)})
  load(file.path(root, "CleanData.RData"))
  assign("data", envir = globalenv(), annika_dat %>%
    select(subjects, trialType, blockIndex, trialIndex, item, itemCategory, Nback,
           task, rt, R, S, Cond, UC, LT) %>%
    filter(Cond == "Control") %>%
    mutate(trialIndex = trialIndex + 1,
           trialType = fct_recode(trialType, "Novel" = "F", "Lure" = "L", "Target" = "N"),
           R = factor(R, levels = c("non_target", "target")),
           S = factor(S, levels = c("non_target", "target")),
           repeatTrial = factor(dplyr::case_when(
             trialType == "Target" | trialType == "Lure" ~ "Repeat", TRUE ~ "Novel"),
             levels = c("Novel", "Repeat"))))
  owd <- setwd(root); on.exit(setwd(owd))
  sys.source("FunctionsControlRDM.R", envir = globalenv())
  data <- get("data", envir = globalenv())
  subDesign <- get("models_rdm", envir = globalenv())[["WeigardModelRDM"]]
  pg <- sampled_pars(subDesign, data = data) + 1
  pg[grepl("^B", names(pg))] <- 2; pg[grepl("^A", names(pg))] <- 2
  pg[grepl("^t0", names(pg))] <- 3; pg[grepl("^pContaminant", names(pg))] <- 4
  make_emc(data = data, design = subDesign, par_groups = pg,
           prior = priorRDM(subDesign, data = data), rt_resolution = 1/60)
}

if (identical(Sys.getenv("EMC_CHILD"), "")) {
  .libPaths(c(lib, .libPaths())); suppressMessages(library(EMC2))
  RNGkind("L'Ecuyer-CMRG"); set.seed(2024)
  emc <- suppressMessages(build())
  # Burn in once, serially, so every timed configuration starts from the same
  # state and measures steady-state cost rather than start-up transients.
  emc <- suppressMessages(run_emc(emc, stage = "preburn",
    stop_criteria = list(iter = 15, max_gd = Inf, min_unique = 0, min_es = 0),
    cores_for_chains = 3, cores_per_chain = 8, verbose = FALSE, step_size = 15,
    max_tries = 1))
  saveRDS(emc, setup_file)
  d <- attr(emc[[1]]$data, "data_list")
  cat(sprintf("\nRDMSWTN control: %d subjects, %d pars, %d chains\n",
              emc[[1]]$n_subjects, emc[[1]]$n_pars, length(emc)))
  for (cpc in as.integer(strsplit(Sys.getenv("EMC_CPC", "1,2,4,8,16"), ",")[[1]]))
    for (flat in strsplit(Sys.getenv("EMC_FLAT", "FALSE"), ",")[[1]])
      system2("Rscript", c(self, cpc, flat),
              env = c("EMC_CHILD=1", paste0("EMC_LIB=", lib),
                      paste0("EMC_SETUP=", setup_file), paste0("EMC_ITER=", iter)))
  quit(save = "no")
}

# --- child --------------------------------------------------------------------
.libPaths(c(lib, .libPaths())); suppressMessages(library(EMC2))
a <- commandArgs(TRUE)
cpc <- as.integer(a[1])
flat <- if (length(a) > 1) as.logical(a[2]) else FALSE
dyn <- as.logical(Sys.getenv("EMC_DYN", "TRUE"))
options(emc2.flat_parallel = flat, emc2.dynamic_cores = dyn)
emc <- readRDS(setup_file)
stamp <- file.path(tempdir(), "annika_times"); dir.create(stamp, showWarnings = FALSE)
unlink(list.files(stamp, full.names = TRUE))
options(emc2.chain_timer_dir = stamp)
RNGkind("L'Ecuyer-CMRG"); set.seed(11)
t <- system.time(emc <- suppressMessages(run_emc(emc, stage = "burn",
  stop_criteria = list(iter = iter, max_gd = Inf, min_unique = 0, min_es = 0),
  cores_for_chains = 3, cores_per_chain = cpc, verbose = FALSE,
  step_size = iter, max_tries = 1)))
cpu <- sum(t[c("user.self", "sys.self", "user.child", "sys.child")])
ct <- sort(vapply(list.files(stamp, full.names = TRUE), readRDS, numeric(1)))
# The flat path has no per-chain processes, so it writes no timings.
if (!length(ct)) ct <- t[["elapsed"]]
cat(sprintf("flat=%-5s dyn=%-5s cpc=%2d total=%2d elapsed=%6.1fs cpu=%7.1fs util=%.2f chains=[%s] spread=%4.1f%% idle_core_s=%5.1f master=%4.1fs\n",
            flat, dyn, cpc, 3 * cpc, t[["elapsed"]], cpu, cpu / (t[["elapsed"]] * 3 * cpc),
            paste(sprintf("%.1f", ct), collapse = " "),
            100 * (max(ct) - min(ct)) / max(ct), cpc * sum(max(ct) - ct),
            t[["elapsed"]] - max(ct)))
