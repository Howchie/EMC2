# The subject-marginal posterior along the ridge, log p(t | data) up to a constant, from a saved chain state by
# importance sampling with common random numbers: per subject K particles from a mixture of its conditional at
# each grid t (inflated), reused for every t; repeated over seeds for the correlated-pseudo-marginal variance.
# LIB=<lib> Rscript marginal.R <chain> <K> <seeds>
args <- commandArgs(TRUE); j <- as.integer(args[1]); K <- as.integer(args[2]); nseed <- as.integer(args[3])
suppressMessages(library(EMC2, lib.loc = Sys.getenv("LIB"))); ns <- asNamespace("EMC2")
load("/data/work/PM/NirvanaHons_Nback/Samples/samples_control_Exp1_frqk7r_blk_gL_sepA_sepL_sepV_sv_t0T_mt.RData")
for (i in 2:length(emc)) for (nm in setdiff(names(emc[[1]]), names(emc[[i]]))) emc[[i]][[nm]] <- emc[[1]][[nm]]
emc <- ns$restore_custom_kernel_pointers(emc, quiet = TRUE)
ch <- emc[[j]]; p <- ch$n_pars; pn <- ch$par_names; n <- ch$n_subjects; idx <- ch$samples$idx
mu <- ch$samples$theta_mu[, idx]; V <- ch$samples$theta_var[, , idx]; alpha <- matrix(ch$samples$alpha[, , idx], p)
Vinv <- chol2inv(chol(V))
u <- readRDS("/data/work/EMC2_dev_oo/plans/alloc_prediction_2026-10-15/handoff_2026-10-15/scripts/ridge_u.rds")
d <- u / sum(u^2)                                  # a step of t moves the group mean t within-sd along the ridge
pos <- function(m) sum(u * m)
others <- sapply(seq_along(emc), function(i) pos(emc[[i]]$samples$theta_mu[, idx]))
cat(sprintf("chain %d at ridge position %.2f; chains at %s\n", j, pos(mu), paste(round(others, 2), collapse = " ")))
tg <- c(-3, -2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2, 3, 4, 5, 6)
# per-subject conditional covariance at the current point and its shift per unit t
Sig <- Del <- vector("list", n)
for (s in seq_len(n)) {
  Q <- ch$lik_laplace[[s]]$prec; if (!is.matrix(Q) || any(!is.finite(Q))) Q <- matrix(0, p, p)
  A <- (Q + t(Q)) / 2 + Vinv; e <- eigen(A, symmetric = TRUE); e$values <- pmax(e$values, 1e-6 * max(e$values))
  Sig[[s]] <- e$vectors %*% diag(1 / e$values) %*% t(e$vectors); Del[[s]] <- drop(Sig[[s]] %*% Vinv %*% d)
}
ldmvn <- function(X, m, S) { L <- chol(S); z <- backsolve(L, t(X) - m, transpose = TRUE); -0.5 * colSums(z^2) - sum(log(diag(L))) - 0.5 * p * log(2 * pi) }
lse <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
prior_mu <- function(t) { m <- mu + t * d; -0.5 * drop(t(m - ch$prior$theta_mu_mean) %*% ch$prior$theta_mu_invar %*% (m - ch$prior$theta_mu_mean)) }
one_seed <- function(seed) {
  set.seed(seed)
  res <- parallel::mclapply(seq_len(n), function(s) {
    S2 <- 2 * Sig[[s]]; cen <- lapply(tg, function(t) alpha[, s] + Del[[s]] * t)
    comp <- sample.int(length(tg), K, replace = TRUE)
    X <- t(sapply(seq_len(K), function(k) cen[[comp[k]]] + drop(chol(S2) %*% rnorm(p)))); colnames(X) <- pn
    lq <- apply(sapply(cen, function(c) ldmvn(X, c, S2)), 1, lse) - log(length(tg))
    ll <- ns$.emc_ensemble_ll_checked(X, ch$data[[s]], ch$model)
    ll[!is.finite(ll)] <- -1e10
    sapply(tg, function(t) { lw <- ll + ldmvn(X, mu + t * d, V) - lq; c(lse(lw) - log(K), exp(2 * lse(lw) - lse(2 * lw))) })
  }, mc.cores = 24)
  logZ <- sapply(res, function(r) r[1, ]); ess <- sapply(res, function(r) r[2, ])
  list(logZ = logZ, ess = ess)
}
runs <- lapply(seq_len(nseed), one_seed)
lp <- sapply(runs, function(r) rowSums(r$logZ) + sapply(tg, prior_mu))
lp <- sweep(lp, 2, lp[tg == 0, ])
ess_med <- apply(sapply(runs, function(r) apply(r$ess, 1, median)), 1, mean)
ess_min <- apply(sapply(runs, function(r) apply(r$ess, 1, min)), 1, mean)
out <- data.frame(t = tg, pos = pos(mu) + tg, log_marg = rowMeans(lp), sd_over_seeds = apply(lp, 1, sd), ess_median = ess_med, ess_min = ess_min)
cat(sprintf("K = %d particles per subject, %d seeds; log p(t|data) - log p(0|data), with the correlated-estimate sd over seeds\n", K, nseed))
print(format(out, digits = 3), row.names = FALSE)
saveRDS(list(out = out, lp = lp, runs = runs), file.path(Sys.getenv("SP"), sprintf("marginal_ch%d_K%d.rds", j, K)))
