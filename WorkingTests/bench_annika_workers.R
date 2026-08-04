# Does the nested path stop spawning subject workers after the first batch?
#
# Runs the real RDMSWTN control fit across several blocks while a background
# poller counts live R processes, so worker spawning can be seen over time
# rather than inferred from a single elapsed time.
#
#   Rscript WorkingTests/bench_annika_workers.R
#
# Writes /tmp/worker_counts.txt : "<seconds since start> <n R processes>".

lib <- Sys.getenv("EMC_LIB", .libPaths()[1])
root <- "/data/work/PM/AnnikaHons_PM1_3"
cpc <- as.integer(Sys.getenv("EMC_CPC", "8"))
iter <- as.integer(Sys.getenv("EMC_ITER", "40"))
step <- as.integer(Sys.getenv("EMC_STEP", "10"))

.libPaths(c(lib, .libPaths())); suppressMessages(library(EMC2))
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
owd <- setwd(root)
sys.source("FunctionsControlRDM.R", envir = globalenv())
setwd(owd)
subDesign <- get("models_rdm", envir = globalenv())[["WeigardModelRDM"]]
data <- get("data", envir = globalenv())
pg <- sampled_pars(subDesign, data = data) + 1
pg[grepl("^B", names(pg))] <- 2; pg[grepl("^A", names(pg))] <- 2
pg[grepl("^t0", names(pg))] <- 3; pg[grepl("^pContaminant", names(pg))] <- 4

RNGkind("L'Ecuyer-CMRG"); set.seed(2024)
emc <- suppressMessages(make_emc(data = data, design = subDesign, par_groups = pg,
                                 prior = priorRDM(subDesign, data = data),
                                 rt_resolution = 1/60))

out <- "/tmp/worker_counts.txt"
poll <- sprintf(
  "start=$(date +%%s.%%N); while true; do n=$(pgrep -c -u $(id -u) -x R); echo \"$(echo \"$(date +%%s.%%N) - $start\" | bc) $n\"; sleep 0.25; done > %s",
  out)
poller <- system2("bash", c("-c", shQuote(poll)), wait = FALSE)
on.exit(system("pkill -f 'pgrep -c -u' ; pkill -f 'while true; do n=' ", ignore.stderr = TRUE))

cat("--- preburn (the initial all-core burst) ---\n")
emc <- suppressMessages(run_emc(emc, stage = "preburn",
  stop_criteria = list(iter = 15, max_gd = Inf, min_unique = 0, min_es = 0),
  cores_for_chains = 3, cores_per_chain = cpc, verbose = FALSE, step_size = 15,
  max_tries = 1))
cat("--- burn, several blocks (watch for respawning) ---\n")
emc <- suppressMessages(run_emc(emc, stage = "burn",
  stop_criteria = list(iter = iter, max_gd = Inf, min_unique = 0, min_es = 0),
  cores_for_chains = 3, cores_per_chain = cpc, verbose = FALSE, step_size = step,
  max_tries = 1))
cat("done\n")
