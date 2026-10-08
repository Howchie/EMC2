# Shared setup: subject 11, chain 1 of the saved lba71sim production checkpoint,
# group state fixed at the chain's last draw.
.libPaths(c(Sys.getenv("LIB"), .libPaths()))
suppressMessages({ library(EMC2); library(posterior) }); ns <- asNamespace("EMC2")
source("R/subject_hmc.R")
CKPT <- "benchmarks/production-adaptation/results_shear/shear_v4/emc_lba71sim_doubling_shearK4_2.rds"
e <- ns$restore_duplicates(readRDS(CKPT))
s <- 11L; ch <- e[[1]]; pn <- ch$par_names; p <- length(pn); dat <- ch$data[[s]]; idx <- ch$samples$idx
mu <- ch$samples$theta_mu[, idx]; V <- ch$samples$theta_var[, , idx]
rows_used <- 0
ll_fn <- function(X) { colnames(X) <- pn; rows_used <<- rows_used + nrow(X)
  as.numeric(ns$calc_ll_manager(X, dadm = dat, model = ch$model, r_cores = 1L)) }
floor_s <- .emc_hmc_ll_floor(dat)
H <- ch$lik_laplace[[s]]$prec
fr <- .emc_hmc_frame(H, mu, V)
a_stuck <- ch$samples$alpha[, s, idx]; ll_stuck <- ch$samples$subj_ll[s, idx]
iA <- match("A", pn); it0 <- match("t0", pn); iB <- match("B", pn)
