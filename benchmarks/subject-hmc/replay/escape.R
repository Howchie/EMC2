source("benchmarks/subject-hmc/replay/setup11.R")
library(parallel)
N <- as.integer(Sys.getenv("N", "300"))
# Count the particle step's likelihood rows through the namespace function.
pcount <- new.env(); pcount$rows <- 0
orig <- ns$calc_ll_manager
counting <- function(proposals, ...) { pcount$rows <- pcount$rows + (if (is.matrix(proposals)) nrow(proposals) else 1L); orig(proposals, ...) }
environment(counting) <- environment()
unlockBinding("calc_ll_manager", ns); assign("calc_ll_manager", counting, envir = ns)
lap <- ns$.emc_laplace_conditional(structure(numeric(p), lik = ch$lik_laplace[[s]]), list(mu = mu, var = V))
eff_mu <- ns$.emc_attach_lik(ch$eff_mu, ch$lik_laplace)[[s]]
pm0 <- attr(ch$samples, "pm_settings")[[s]]
tune <- list(components = rep(1, p), shared_ll_idx = rep(1, p))
alpha_all <- ch$samples$alpha[, , idx]
particle <- function(a, ll, pm) {
  pars <- list(tmu = mu, tvar = V, subj_mu = matrix(mu, p, ch$n_subjects), alpha = alpha_all); pars$alpha[, s] <- a
  before <- pcount$rows
  r <- ns$new_particle(s, dat, pm, eff_mu = eff_mu, eff_var = ch$eff_var[[s]], chains_mu = ch$chains_mu[[s]],
                       chains_var = ch$chains_var[[s]], prev_ll = ll, parameters = pars, model = ch$model,
                       stage = "sample", type = ch$type, tune = tune)
  list(a = r$proposal, ll = r$ll, pm = r$pm_settings, rows = pcount$rows - before)
}
EPS <- as.numeric(Sys.getenv("EPS", "0.2"))
run <- function(start, arm, seed) {
  set.seed(seed); a <- start; ll <- ll_fn(matrix(a, 1)); pm <- pm0
  A <- numeric(N); rows <- numeric(N); acc <- rep(NA, N)
  rate <- switch(arm, P = 0, H = 1, PH = 1, PH10 = 0.1)
  for (t in 1:N) {
    r <- 0
    if (arm != "H") { pr <- particle(a, ll, pm); a <- pr$a; ll <- pr$ll; pm <- pr$pm; r <- r + pr$rows }
    if (rate > 0 && runif(1) < rate) { rows_used <<- 0
      tr <- .emc_hmc_transition(a, ll, fr, ll_fn, EPS, 1e-4, floor_s); a <- tr$alpha; ll <- tr$ll; r <- r + tr$rows; acc[t] <- tr$acc_prob }
    A[t] <- a[iA]; rows[t] <- r
  }
  list(A = A, rows = rows, acc = acc)
}
# Reference conditional: long HMC chain from the Laplace mean.
ref <- local({ set.seed(9); a <- lap$mu; ll <- ll_fn(matrix(a, 1)); A <- numeric(1500)
  for (t in 1:1500) { tr <- .emc_hmc_transition(a, ll, fr, ll_fn, EPS, 1e-4, floor_s); a <- tr$alpha; ll <- tr$ll; A[t] <- a[iA] }; A[-(1:100)] })
cat(sprintf("reference (HMC, 1400 its): A median %.3f mean %.3f sd %.3f ESS %.0f\n", median(ref), mean(ref), sd(ref), ess_bulk(ref)))
set.seed(4); Ls <- t(chol(lap$var))
starts <- list(stuck_c1 = a_stuck, chain2 = e[[2]]$samples$alpha[, s, e[[2]]$samples$idx], chain3 = e[[3]]$samples$alpha[, s, e[[3]]$samples$idx])
for (k in 1:3) starts[[paste0("disp", k)]] <- lap$mu + 2 * drop(Ls %*% rnorm(p))
jobs <- expand.grid(start = names(starts), arm = c("P", "H", "PH", "PH10"), stringsAsFactors = FALSE)
res <- mclapply(seq_len(nrow(jobs)), function(k) run(starts[[jobs$start[k]]], jobs$arm[k], 1000 + k), mc.cores = min(24, nrow(jobs)))
med <- median(ref)
for (k in seq_len(nrow(jobs))) {
  r <- res[[k]]; if (!is.list(r)) { cat(jobs$start[k], jobs$arm[k], "FAILED\n"); next }
  A0 <- starts[[jobs$start[k]]][iA]; side <- sign(A0 - med)
  esc <- if (abs(A0 - med) < 0.5 * sd(ref)) NA else which(sign(r$A - med) != side)[1]
  keep <- -(1:(N %/% 3))
  cat(sprintf("%-8s %-4s A0 %6.3f | escape it %4s | last2/3: A mean %6.3f sd %.3f ESS %5.1f | rows/it %6.0f | ESS/krow %.4f | acc %.2f\n",
    jobs$start[k], jobs$arm[k], A0, ifelse(is.na(esc), "-", esc), mean(r$A[keep]), sd(r$A[keep]), ess_bulk(r$A[keep]),
    mean(r$rows), 1000 * ess_bulk(r$A[keep]) / sum(r$rows[keep]), mean(r$acc, na.rm = TRUE)))
}
saveRDS(list(jobs = jobs, res = res, ref = ref, starts = starts), Sys.getenv("OUT", "escape.rds"))
