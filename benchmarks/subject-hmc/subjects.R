# Per-subject chain diagnostics of saved fit_arm.R fits (SAVE_EMC=1).
#   EMC_LIB=<lib> Rscript benchmarks/subject-hmc/subjects.R <fit_emc.rds> [par ...]
# Per subject: max rank-normalised R-hat and min bulk ESS over its
# parameters, and for each named parameter the narrowest chain's sd relative
# to the pooled sd (a chain held in one region shows a small ratio).
args <- commandArgs(TRUE)
suppressPackageStartupMessages({ library(EMC2, lib.loc = Sys.getenv("EMC_LIB")); library(posterior) })
emc <- readRDS(args[1])
pars <- if (length(args) > 1) args[-1] else c("A", "B", "t0")
pn <- emc[[1]]$par_names
pars <- intersect(pars, pn)
keep <- lapply(emc, function(ch) which(ch$samples$stage == "sample" & seq_along(ch$samples$stage) <= ch$samples$idx))
n <- min(lengths(keep))
tab <- do.call(rbind, lapply(seq_len(emc[[1]]$n_subjects), function(s) {
  a <- simplify2array(lapply(seq_along(emc), function(j)
    t(emc[[j]]$samples$alpha[, s, utils::tail(keep[[j]], n), drop = TRUE])))   # iter x par x chain
  a <- aperm(a, c(1, 3, 2))
  rh <- apply(a, 3, rhat); es <- apply(a, 3, ess_bulk)
  ratio <- vapply(pars, function(k) {
    x <- a[, , match(k, pn)]
    min(apply(x, 2, stats::sd)) / stats::sd(as.vector(x))
  }, numeric(1))
  data.frame(subject = s, max_rhat = round(max(rh), 3), min_ess = round(min(es)),
             worst_par = pn[which.max(rh)], t(round(ratio, 2)))
}))
print(tab, row.names = FALSE)
cat(sprintf("subjects with max R-hat > 1.1: %d; > 1.01: %d; min ESS < 100: %d\n",
            sum(tab$max_rhat > 1.1), sum(tab$max_rhat > 1.01), sum(tab$min_ess < 100)))
