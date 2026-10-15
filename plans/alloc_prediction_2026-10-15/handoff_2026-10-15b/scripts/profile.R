# Two profiles of the posterior along the ridge from a saved chain state, t in within-sd units of the group mean:
#  fixed  : the joint log density along the shear's fixed-residual line (each subject shifted by its Laplace
#           conditional shift Delta_s t), relative to t = 0;
#  marginal: the Laplace (profile) approximation to log p(t | data): each subject re-optimised in its
#           conditional at mu + t d, relative to t = 0.
# LIB=<lib> Rscript profile.R <chain>
args <- commandArgs(TRUE); j <- as.integer(args[1])
suppressMessages(library(EMC2, lib.loc = Sys.getenv("LIB"))); ns <- asNamespace("EMC2")
load("/data/work/PM/NirvanaHons_Nback/Samples/samples_control_Exp1_frqk7r_blk_gL_sepA_sepL_sepV_sv_t0T_mt.RData")
for (i in 2:length(emc)) for (nm in setdiff(names(emc[[1]]), names(emc[[i]]))) emc[[i]][[nm]] <- emc[[1]][[nm]]
emc <- ns$restore_custom_kernel_pointers(emc, quiet = TRUE)
ch <- emc[[j]]; p <- ch$n_pars; pn <- ch$par_names; n <- ch$n_subjects; idx <- ch$samples$idx
mu <- ch$samples$theta_mu[, idx]; V <- ch$samples$theta_var[, , idx]; alpha <- matrix(ch$samples$alpha[, , idx], p)
Vinv <- chol2inv(chol(V)); Lv <- chol(V)
u <- readRDS("/data/work/EMC2_dev_oo/plans/alloc_prediction_2026-10-15/handoff_2026-10-15/scripts/ridge_u.rds")
d <- u / sum(u^2); pos <- function(m) sum(u * m)
others <- sapply(seq_along(emc), function(i) pos(emc[[i]]$samples$theta_mu[, idx]))
cat(sprintf("chain %d at ridge position %.2f; chains at %s\n", j, pos(mu), paste(round(others, 2), collapse = " ")))
tg <- c(-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4, 5, 6)
lprior <- function(a, m) -0.5 * sum(backsolve(Lv, a - m, transpose = TRUE)^2)
ll1 <- function(s, a) ns$.emc_ensemble_ll_checked(matrix(a, 1, dimnames = list(NULL, pn)), ch$data[[s]], ch$model)
llm <- function(s, A) ns$.emc_ensemble_ll_checked(A, ch$data[[s]], ch$model)
res <- parallel::mclapply(seq_len(n), function(s) {
  Q <- ch$lik_laplace[[s]]$prec; Q <- (Q + t(Q)) / 2
  A <- Q + Vinv; e <- eigen(A, symmetric = TRUE); e$values <- pmax(e$values, 1e-6 * max(e$values))
  Sig <- e$vectors %*% diag(1 / e$values) %*% t(e$vectors); Del <- drop(Sig %*% Vinv %*% d)
  a0 <- alpha[, s]; l0 <- ll1(s, a0) + lprior(a0, mu)
  fixed <- marg <- numeric(length(tg)); a_prev <- a0
  for (k in seq_along(tg)) {
    t <- tg[k]; m <- mu + t * d
    af <- a0 + Del * t; fixed[k] <- ll1(s, af) + lprior(af, m) - l0
    f <- function(a) { v <- ll1(s, a) + lprior(a, m); if (!is.finite(v)) 1e10 else -v }
    gr <- function(a) { h <- 1e-4; X <- rbind(t(a + h * diag(p)), t(a - h * diag(p))); colnames(X) <- pn
      v <- llm(s, X) + apply(X, 1, lprior, m = m); v[!is.finite(v)] <- -1e10; -(v[1:p] - v[p + 1:p]) / (2 * h) }
    starts <- list(af, a_prev); best <- NULL
    for (st in starts) { o <- tryCatch(optim(st, f, gr, method = "BFGS", control = list(maxit = 60)), error = function(e) NULL)
      if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o }
    marg[k] <- -best$value - l0; a_prev <- best$par
  }
  rbind(fixed = fixed, marg = marg)
}, mc.cores = 24)
fixed <- sapply(res, function(r) r["fixed", ]); marg <- sapply(res, function(r) r["marg", ])
pm <- sapply(tg, function(t) { m <- mu + t * d; -0.5 * drop(t(m - ch$prior$theta_mu_mean) %*% ch$prior$theta_mu_invar %*% (m - ch$prior$theta_mu_mean)) })
out <- data.frame(t = tg, pos = pos(mu) + tg, fixed_line = rowSums(fixed) + pm - pm[tg == 0], marginal = rowSums(marg) + pm - pm[tg == 0],
  marg_subject_q10 = apply(marg, 1, quantile, 0.1), marg_subject_med = apply(marg, 1, median), marg_subject_q90 = apply(marg, 1, quantile, 0.9))
cat("log density relative to t = 0 (summed over subjects, group prior included); per-subject marginal quantiles\n")
print(format(out, digits = 3), row.names = FALSE)
saveRDS(list(out = out, fixed = fixed, marg = marg, tg = tg, pos0 = pos(mu)), file.path(Sys.getenv("SP"), sprintf("profile_ch%d.rds", j)))
