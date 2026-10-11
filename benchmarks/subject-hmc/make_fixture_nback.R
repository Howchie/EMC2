# N-back fixtures from a fitted emc object (real data, kept out of the repo:
# written to FIXTURE_DIR).
#   EMC_LIB=<lib> Rscript benchmarks/subject-hmc/make_fixture_nback.R <emc.rds> <fixture_dir>
# Writes fixture_nback.rds (all subjects) and fixture_nback31.rds (every
# third subject, a cheaper proxy with the same design and prior). Each fixture
# holds data, design, prior and the sampler type make_emc() needs.
args <- commandArgs(TRUE)
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
emc <- readRDS(args[1]); dir <- args[2]
dir.create(dir, showWarnings = FALSE, recursive = TRUE)
data <- get_data(emc)
design <- get_design(emc)
prior <- emc[[1]]$prior
type <- emc[[1]]$type
subjects <- as.character(unique(data$subjects))
write_fx <- function(keep, name) {
  d <- data[data$subjects %in% keep, ]
  d$subjects <- factor(as.character(d$subjects), levels = keep)
  fx <- list(data = d, design = design, prior = prior, type = type,
             settings = list(source = basename(args[1]), subjects = keep))
  saveRDS(fx, file.path(dir, paste0("fixture_", name, ".rds")))
  cat(sprintf("%s: %d subjects, %d rows\n", name, length(keep), nrow(d)))
}
write_fx(subjects, "nback")
write_fx(subjects[seq(1, length(subjects), by = 3)], "nback31")
