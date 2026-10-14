# Offline dragging with per-subject bridge lengths chosen exactly.
# Each subject's bridge length is K_long with probability p(|W1|), else K_short,
# where W1 is the subject's one-step work (endpoint log-density difference) at
# the bridge's start. The reverse move shares the midpoint chart, so its start
# statistic is the same difference at the bridge's final point; the acceptance
# carries q(K | reverse) / q(K | forward) per subject.
# Usage: Rscript variable_drag.R fixture.rds output.rds [draws] [reps] [designs]
#   designs: comma-separated Kshort:Klong:tau (tau in nats; tau=0 means uniform Kshort)
suppressPackageStartupMessages(library(EMC2))
ns <- asNamespace("EMC2"); args <- commandArgs(TRUE)
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(args[1])), quiet = TRUE)
draws <- if (length(args) >= 3) as.integer(strsplit(args[3], ",")[[1]]) else c(1001, 1000, 980, 950, 900, 800, 700, 600)
REPS <- if (length(args) >= 4) as.integer(args[4]) else 4L
designs <- if (length(args) >= 5) strsplit(args[5], ",")[[1]] else c("1:0:0", "4:0:0", "16:0:0", "1:64:2", "4:64:2", "4:32:2", "4:64:5")
designs <- lapply(designs, function(x) { v <- as.numeric(strsplit(x, ":")[[1]]); list(Ks = v[1], Kl = v[2], tau = v[3], name = x) })
EPS <- 0.02
ch <- emc[[1]]; p <- ch$n_pars; n <- ch$n_subjects; l <- match("B_Espeed", ch$par_names)
idx <- ns$.emc_chain_proposal_window(emc, seq_len(ch$samples$idx), "sample")
ctx <- ns$.emc_shear_context(ch); spec <- ctx$spec; walls <- ctx$walls; fit <- ctx$packed
# Measured lbaE likelihood call cost (ms) by rows; a shear round is one 8-row call per subject.
call_ms <- function(rows) approx(c(1, 2, 8, 16), c(0.230, 0.298, 0.638, 1.098), rows, rule = 2)$y
shear_round_ms <- n * call_ms(8)
RNGkind("L'Ecuyer-CMRG"); set.seed(sum(utf8ToInt("variable-drag")))
p_long <- function(w, tau) if (tau <= 0) 0 else EPS + (1 - 2 * EPS) * w^2 / (w^2 + tau^2)
log_q <- function(long, w, tau) { if (tau <= 0) return(0); pl <- p_long(w, tau); if (long) log(pl) else log1p(-pl) }
results <- list()
for (d in draws) for (jump in c(-1, -0.5, 0.5, 1)) {
  mu <- ch$samples$theta_mu[, d]; V <- ch$samples$theta_var[, , d]
  alpha <- matrix(ch$samples$alpha[, , d], p)
  path <- ns$.emc_shear_coord_path(ch, idx, p + l, mu, V, c(0, jump / 2, jump))
  st <- lapply(path, function(g) ns$.emc_shear_state(g$mu, g$V, fit, ch$prior, spec))
  gp <- function(g) ns$.emc_shear_log_prior_v(g$V, spec) + sum((spec$size + 1) * g$u_sd)
  group_delta <- gp(path[[3]]) - gp(path[[1]])
  subj <- lapply(seq_len(n), function(s) {
    f <- lapply(st, function(g) list(m = g$fm$m[, s], R = g$fm$R[, , s], logdet = g$fm$logdet[s]))
    haswall <- !is.null(walls) && walls$has[s]
    transport <- function(a, f0, f1) {
      if (haswall) ns$.emc_shear_wall_map(a, f0, f1, walls$U[, s], walls$b[s], walls$Uperp[[s]]) else
        list(a = f1$m + drop(backsolve(f1$R, f0$R %*% (a - f0$m))), log_jac = f0$logdet - f1$logdet)
    }
    am <- transport(alpha[, s], f[[1]], f[[2]])$a
    z0 <- drop(f[[2]]$R %*% (am - f[[2]]$m))
    eval_z <- function(z) {
      a <- f[[2]]$m + drop(backsolve(f[[2]]$R, z))
      tr <- lapply(c(1, 3), function(g) transport(a, f[[2]], f[[g]]))
      X <- t(vapply(tr, `[[`, numeric(p), "a"))
      ll <- ns$.emc_burn_ll_batch(X, ch$data[[s]], ch$model, ch$par_names)
      vapply(1:2, function(i) {
        g <- c(1, 3)[i]; zz <- forwardsolve(t(st[[g]]$chol_V), X[i, ] - mu)
        ll[i] - 0.5 * sum(zz^2) - sum(log(diag(st[[g]]$chol_V))) + tr[[i]]$log_jac - f[[2]]$logdet
      }, 1)
    }
    list(z0 = z0, eval_z = eval_z, start = eval_z(z0))
  })
  bridge <- function(sb, K) {
    z <- sb$z0; F <- sb$start; work <- F[2] - F[1]
    for (i in seq_len(K - 1L)) {
      beta <- i / K
      eps <- sample(c(0.4, 0.8, 1), 1)
      zp <- sqrt(1 - eps^2) * z + eps * rnorm(p)
      FF <- sb$eval_z(zp)
      lr <- sum(c(1 - beta, beta) * (FF - F)) + 0.5 * (sum(zp^2) - sum(z^2))
      if (log(runif(1)) < lr) { z <- zp; F <- FF }
      work <- work + F[2] - F[1]
    }
    list(work = work / K, w_end = abs(F[2] - F[1]), calls = K)
  }
  for (ds in designs) {
    lr <- numeric(REPS); calls <- matrix(0L, REPS, n); nlong <- numeric(REPS); qpen <- numeric(REPS)
    for (r in seq_len(REPS)) {
      tot <- group_delta
      for (s in seq_len(n)) {
        w0 <- abs(subj[[s]]$start[2] - subj[[s]]$start[1])
        long <- runif(1) < p_long(w0, ds$tau)
        b <- bridge(subj[[s]], if (long) ds$Kl else ds$Ks)
        q <- log_q(long, b$w_end, ds$tau) - log_q(long, w0, ds$tau)
        tot <- tot + b$work + q
        qpen[r] <- qpen[r] + q; calls[r, s] <- b$calls; nlong[r] <- nlong[r] + long
      }
      lr[r] <- tot
    }
    acc <- mean(pmin(1, exp(lr)))
    ms_serial <- mean(rowSums(call_ms(2) * calls))
    row <- data.frame(draw = d, jump = jump, design = ds$name, acceptance = acc,
                      mean_log_ratio = mean(lr), mean_q_penalty = mean(qpen), n_long = mean(nlong),
                      rows_per_subject = 2 * mean(calls), ms_total = ms_serial,
                      ms_critical = mean(apply(calls, 1, max)) * call_ms(2),
                      esjd = jump^2 * acc)
    print(row, digits = 3); results[[length(results) + 1]] <- row
  }
  saveRDS(list(results = do.call(rbind, results), shear_round_ms = shear_round_ms,
               head = system("git rev-parse HEAD", intern = TRUE)), args[2])
}
res <- do.call(rbind, results)
agg <- aggregate(cbind(acceptance, esjd, mean_q_penalty, n_long, rows_per_subject, ms_total) ~ design + abs(jump), res, mean)
agg$esjd_per_ms <- agg$esjd / agg$ms_total
print(agg, digits = 3)
cat("shear round: ~", round(shear_round_ms), "ms total likelihood; measured E per round 0.10 (window sd^2) at its tuned scale\n")
