# Offline reversible conditional bridges in a symmetric midpoint chart.
# Usage: Rscript dragging.R fixture.rds output.rds [draws] [K] [reps] [cores]
suppressPackageStartupMessages(library(EMC2))
ns <- asNamespace("EMC2"); args <- commandArgs(TRUE)
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(args[1])), quiet = TRUE)
draws <- if (length(args) >= 3) as.integer(strsplit(args[3], ",")[[1]]) else c(1001, 980, 900)
Ks <- if (length(args) >= 4) as.integer(strsplit(args[4], ",")[[1]]) else c(1, 4, 8)
REPS <- if (length(args) >= 5) as.integer(args[5]) else 4L
NC <- if (length(args) >= 6) as.integer(args[6]) else 3L
ch <- emc[[1]]; p <- ch$n_pars; n <- ch$n_subjects; l <- match("B_Espeed", ch$par_names)
idx <- ns$.emc_chain_proposal_window(emc, seq_len(ch$samples$idx), "sample")
ctx <- ns$.emc_shear_context(ch); spec <- ctx$spec; walls <- ctx$walls
fit <- ctx$packed
RNGkind("L'Ecuyer-CMRG")
set.seed(sum(utf8ToInt("conditional-bridge-replay")))
results <- list()
profile <- rep(1, n)
allocated <- Sys.getenv("BRIDGE_ALLOCATION") == "work"
if (allocated) {
  pilot_draws <- setdiff(c(600, 700, 800, 900), draws)
  stopifnot(length(pilot_draws) >= 3)
  residual <- matrix(0, n, 0)
  for (dd in pilot_draws) {
    mu <- ch$samples$theta_mu[, dd]; V <- ch$samples$theta_var[, , dd]
    aa <- matrix(ch$samples$alpha[, , dd], p)
    gg <- c(-1, -0.5, 0, 0.5, 1); g0 <- which(gg == 0)
    path <- ns$.emc_shear_coord_path(ch, idx, p + l, mu, V, gg)
    states <- lapply(path, function(g) ns$.emc_shear_state(g$mu, g$V, fit, ch$prior, spec))
    tr <- lapply(states, function(g) ns$.emc_shear_transport(aa, states[[g0]], g, walls)$alpha)
    rr <- vapply(seq_len(n), function(s) {
      X <- t(vapply(tr, function(a) a[, s], numeric(p)))
      ll <- ns$.emc_burn_ll_batch(X, ch$data[[s]], ch$model, ch$par_names)
      r <- ll - drop(X %*% fit$lin[, s]) + 0.5 * rowSums((X %*% fit$cube[, , s]) * X)
      (r - r[g0])[-g0]
    }, numeric(length(gg) - 1L))
    residual <- cbind(residual, t(rr))
  }
  profile <- sqrt(rowMeans(residual^2))
  cat("frozen work allocation: pilot draws", pilot_draws, "\n")
}
for (d in draws) for (jump in c(-1, -0.5, 0.5, 1)) {
  mu <- ch$samples$theta_mu[, d]; V <- ch$samples$theta_var[, , d]
  alpha <- matrix(ch$samples$alpha[, , d], p)
  path <- ns$.emc_shear_coord_path(ch, idx, p + l, mu, V, c(0, jump / 2, jump))
  st <- lapply(path, function(g) ns$.emc_shear_state(g$mu, g$V, fit, ch$prior, spec))
  gp <- function(g) ns$.emc_shear_log_prior_v(g$V, spec) + sum((spec$size + 1) * g$u_sd)
  group_delta <- gp(path[[3]]) - gp(path[[1]])
  for (K in Ks) {
    lengths <- rep(K, n)
    if (allocated && K > 1) {
      lengths <- rep(1L, n)
      # Greedily minimize sum_s work_s^2 / K_s at a fixed bridge budget.
      for (extra in seq_len(n * (K - 1))) {
        gain <- profile^2 / (lengths * (lengths + 1))
        gain[lengths >= 64L] <- -Inf
        ss <- which.max(gain); lengths[ss] <- lengths[ss] + 1L
      }
    }
    begin <- proc.time()
    jobs <- parallel::mclapply(seq_len(n), function(s) {
      K_local <- lengths[s]
      f <- lapply(st, function(g) list(m = g$fm$m[, s], R = g$fm$R[, , s], logdet = g$fm$logdet[s]))
      haswall <- !is.null(walls) && walls$has[s]
      transport <- function(a, f0, f1) {
        if (haswall) ns$.emc_shear_wall_map(a, f0, f1, walls$U[, s], walls$b[s], walls$Uperp[[s]]) else
          list(a = f1$m + drop(backsolve(f1$R, f0$R %*% (a - f0$m))), log_jac = f0$logdet - f1$logdet)
      }
      am <- transport(alpha[, s], f[[1]], f[[2]])$a
      z0 <- drop(f[[2]]$R %*% (am - f[[2]]$m))
      calls <- 0L
      eval_z <- function(z) {
        a <- f[[2]]$m + drop(backsolve(f[[2]]$R, z))
        tr <- lapply(c(1, 3), function(g) transport(a, f[[2]], f[[g]]))
        X <- t(vapply(tr, `[[`, numeric(p), "a"))
        ll <- ns$.emc_burn_ll_batch(X, ch$data[[s]], ch$model, ch$par_names)
        calls <<- calls + 1L
        vapply(1:2, function(i) {
          g <- c(1, 3)[i]; zz <- forwardsolve(t(st[[g]]$chol_V), X[i, ] - mu)
          ll[i] - 0.5 * sum(zz^2) - sum(log(diag(st[[g]]$chol_V))) + tr[[i]]$log_jac - f[[2]]$logdet
        }, 1)
      }
      start <- eval_z(z0)
      repl <- t(vapply(seq_len(REPS), function(r) {
        z <- z0; F <- start; work <- F[2] - F[1]; moved <- 0
        for (i in seq_len(K_local - 1L)) {
          beta <- i / K_local
          # A fixed mixture of pCN scales is reversible with the bridge density.
          eps <- sample(c(0.4, 0.8, 1), 1)
          zp <- sqrt(1 - eps^2) * z + eps * rnorm(p)
          FF <- eval_z(zp)
          lr <- sum(c(1 - beta, beta) * (FF - F)) + 0.5 * (sum(zp^2) - sum(z^2))
          if (log(runif(1)) < lr) { z <- zp; F <- FF; moved <- moved + 1 }
          work <- work + F[2] - F[1]
        }
        c(work = work / K_local, moved = moved, displacement = sum((z - z0)^2))
      }, numeric(3)))
      list(repl = repl, calls = calls, cpu = sum(proc.time()[1:2]))
    }, mc.cores = NC, mc.set.seed = TRUE)
    work <- Reduce(`+`, lapply(jobs, function(j) j$repl[, "work"]))
    lr <- work + group_delta
    row <- data.frame(draw = d, jump = jump, K = K, reps = REPS,
                      acceptance = mean(pmin(1, exp(lr))), mean_log_ratio = mean(lr), sd_log_ratio = sd(lr),
                      bridge_accept = if (K > 1) sum(vapply(jobs, function(j) sum(j$repl[, "moved"]), 1)) / (REPS * sum(lengths - 1)) else NA,
                      max_K = max(lengths), active_subjects = sum(lengths > 1),
                      elapsed = (proc.time() - begin)[[3]],
                      likelihood_calls = sum(vapply(jobs, `[[`, 1L, "calls")),
                      likelihood_rows = 2 * sum(vapply(jobs, `[[`, 1L, "calls")))
    print(row)
    results[[length(results) + 1]] <- row
    saveRDS(list(results = do.call(rbind, results), head = system("git rev-parse HEAD", intern = TRUE),
                 idx = idx, seed = "conditional-bridge-replay", profile = profile,
                 allocated = allocated), args[2])
  }
}
