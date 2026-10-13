# Decision stability: block-resample one block's iterations into the
# allocator's inputs and solve again, holding everything else fixed.
#   EMC_LIB=<lib> Rscript benchmarks/allocator/bootstrap.R <run_dir> <k> [n_rep] [variants]
# k: a rebuild whose input state pre_<k>.rds carries a full block trace (k >= 2).
# variants (comma-separated, default "all,rest,moves,costs"): which inputs
# are resampled:
#   rest   the rest-of-iteration jumps (group G; subjects A and D);
#   moves  the moves' own evidence (shear and collapsed block records, HMC
#          transitions, particle steps);
#   costs  every recorded cost;
#   all    all of them, from the same resampled iterations.
# Also replays the rebuild once without resampling and checks that it
# reproduces the plan the run made (replay determinism).
# Writes <run_dir>/boot_<k>.rds.
args <- commandArgs(TRUE)
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
cat("EMC2 from", find.package("EMC2"), "\n")
ns <- asNamespace("EMC2")
run_dir <- args[1]
k <- as.integer(args[2])
n_rep <- if (length(args) >= 3) as.integer(args[3]) else 20L
variants <- if (length(args) >= 4) strsplit(args[4], ",")[[1]] else c("all", "rest", "moves", "costs")
L_block <- 10L

pre <- readRDS(file.path(run_dir, sprintf("pre_%02d.rds", k)))
actual <- readRDS(file.path(run_dir, sprintf("alloc_%02d.rds", k)))
emc0 <- ns$restore_custom_kernel_pointers(pre$emc, quiet = TRUE)
traces <- lapply(emc0, `[[`, "ws1_trace")
for (j in seq_along(emc0)) emc0[[j]]$ws1_trace <- NULL
# A resumed production fit runs under the options its chains recorded.
options(ns$.emc_record_production_options(emc0)$options)

# The state the allocator receives: add_proposals() as the run made it, with
# the rebuild itself captured instead of run.
captured <- NULL
orig_rebuild <- ns$.emc_alloc_rebuild
capture <- function(emc, stage, n_cores, B) { captured <<- list(emc = emc, stage = stage, n_cores = n_cores, B = B); emc }
environment(capture) <- environment()
utils::assignInNamespace(".emc_alloc_rebuild", capture, ns)
if (!is.null(pre$seed)) assign(".Random.seed", pre$seed, envir = globalenv())
invisible(ns$add_proposals(emc0, pre$stage, pre$n_cores, pre$n_blocks, pre$block_iter))
utils::assignInNamespace(".emc_alloc_rebuild", orig_rebuild, ns)
post <- captured
seed_post <- get(".Random.seed", envir = globalenv())

solve_plan <- function(emc) {
  assign(".Random.seed", seed_post, envir = globalenv())
  out <- ns$.emc_alloc_rebuild(emc, post$stage, post$n_cores, post$B)
  lapply(out, function(ch) ch$alloc$last)
}
summ <- function(recs) do.call(rbind, lapply(seq_along(recs), function(j) {
  r <- recs[[j]]
  p2 <- length(r$rounds); p <- p2 / 2
  data.frame(chain = j, tau = r$tau, tau0 = r$tau0, regime = r$regime,
             rounds_mu = sum(r$rounds[seq_len(p)]), rounds_lsd = sum(r$rounds[p + seq_len(p)]),
             coll = r$collapsed_rounds, hmc_n = sum(r$rate > 0), hmc_mean = mean(r$rate),
             q_mean = mean(r$particle_rate), q_n = sum(r$particle_rate < 1),
             bind = which.max(c(r$pred_group, apply(r$pred_subject, 1, max))))
}))

# Replay determinism: the same input, the same plan.
replay <- solve_plan(post$emc)
same <- vapply(seq_along(replay), function(j) {
  a <- actual$records[[j]]; b <- replay[[j]]
  isTRUE(all.equal(a$rounds, b$rounds)) && isTRUE(all.equal(a$rate, b$rate)) &&
    isTRUE(all.equal(a$collapsed_rounds, b$collapsed_rounds)) &&
    isTRUE(all.equal(a$particle_rate, b$particle_rate)) && isTRUE(all.equal(a$tau, b$tau))
}, logical(1))
cat("replay reproduces the run's plan per chain:", same, "\n")
for (j in which(!same)) {
  a <- actual$records[[j]]; b <- replay[[j]]
  diffs <- names(a)[!vapply(names(a), function(f) isTRUE(all.equal(a[[f]], b[[f]])), logical(1))]
  cat(" chain", j, "differs in:", diffs, "\n")
  for (f in intersect(diffs, c("tau", "rounds", "hmc_h", "calib_cpu", "collapsed_rounds", "rate")))
    cat("   ", f, "run:", signif(a[[f]], 4), "\n    replay:", signif(b[[f]], 4), "\n")
}

# --- reconstruction of the block's inputs from resampled iterations ----------
pad <- function(M, w) if (ncol(M) < w) cbind(M, matrix(0L, nrow(M), w - ncol(M))) else M
diff_rec <- function(a, b) {  # b - a for block records (K padded)
  out <- b
  for (f in setdiff(names(b), "K")) out[[f]] <- b[[f]] - (if (is.null(a)) 0 else a[[f]])
  if (!is.null(b$K)) {
    w <- max(ncol(b$K), if (is.null(a)) 1L else ncol(a$K))
    out$K <- pad(b$K, w) - (if (is.null(a)) 0L else pad(a$K, w))
  }
  out
}
sum_rec <- function(incs, idx) {
  if (!length(incs)) return(NULL)
  w <- max(vapply(incs, function(x) if (is.null(x$K)) 1L else ncol(x$K), integer(1)))
  out <- incs[[1]]
  for (f in names(out)) out[[f]] <- 0 * (if (f == "K") pad(out$K, w) else out[[f]])
  for (i in idx) for (f in names(out)) out[[f]] <- out[[f]] + (if (f == "K") pad(incs[[i]]$K, w) else incs[[i]][[f]])
  out
}

chain_inputs <- function(tr, idx, p, n) {
  w <- ns$.emc_alloc_empty_window(n, p)
  w$D <- replicate(n, matrix(0, p, p), simplify = FALSE)
  for (i in idx) {
    d <- tr$group[[i]]
    if (all(is.finite(d))) { w$G <- w$G + tcrossprod(d); w$nG <- w$nG + 1L }
    sj <- tr$subjects[[i]]
    if (all(is.finite(sj$d))) {
      dd <- sj$d - sj$dmu
      for (s in seq_len(n)) {
        w$A[[s]] <- w$A[[s]] + tcrossprod(sj$d[, s])
        w$D[[s]] <- w$D[[s]] + tcrossprod(dd[, s])
      }
      w$nA <- w$nA + 1L
    }
    pt <- tr$particle[[i]]
    ran <- which(pt$ran)
    if (length(ran)) {
      t_ran <- if (is.null(pt$times)) rep(0, length(ran)) else pmax(0, pt$times[ran])
      share <- if (sum(t_ran) > 0) t_ran / sum(t_ran) else rep(1 / length(ran), length(ran))
      for (m in seq_along(ran)) {
        s <- ran[m]; dv <- pt$d[, s]
        if (!all(is.finite(dv))) next
        w$P[[s]] <- w$P[[s]] + tcrossprod(dv)
        w$nP[s] <- w$nP[s] + 1L
        w$cpu_part[s] <- w$cpu_part[s] + if (sum(t_ran) > 0) t_ran[m] else pt$cost * share[m]
      }
    }
    cs <- tr$cost[[i]]
    w$n <- w$n + 1L
    w$cpu_rest <- w$cpu_rest + cs$rest
    moves <- cs$moves
    work <- 0; att <- 0L
    for (s in seq_along(cs$hmc)) {
      tk <- cs$hmc[[s]]
      if (is.null(tk)) next
      c_s <- if (!is.null(tk$work)) tk$work else if (!is.null(tk$cost)) tk$cost else tk$cpu
      w$cpu_hmc[s] <- w$cpu_hmc[s] + c_s
      w$hmc_att[s] <- w$hmc_att[s] + 1L
      work <- work + c_s; att <- att + 1L
    }
    if (att > 0L) moves$hmc <- c(cpu = as.numeric(moves$hmc)[1], own = work, n = att)
    for (fam in intersect(names(moves), ns$.EMC_ALLOC_FAMILIES)) {
      m <- moves[[fam]]
      if (is.null(m) || length(m) < 3L || !is.finite(m[[3]]) || m[[3]] <= 0) next
      f <- w$fam[[fam]]
      f$cpu <- f$cpu + as.numeric(m[[1]]); f$own <- f$own + as.numeric(m[[2]])
      f$iters <- f$iters + 1L; f$moves <- f$moves + as.integer(round(m[[3]]))
      w$fam[[fam]] <- f
    }
  }
  w
}

hmc_incs <- function(tr, idx, n, p) {
  out <- lapply(seq_len(n), function(s) list(n = 0, M = matrix(0, p, p), cpu = 0))
  for (i in idx) for (s in seq_along(tr$cost[[i]]$hmc)) {
    tk <- tr$cost[[i]]$hmc[[s]]
    if (is.null(tk) || is.null(tk$disp)) next
    out[[s]]$n <- out[[s]]$n + 1
    out[[s]]$M <- out[[s]]$M + tk$acc_prob * tcrossprod(tk$disp)
    out[[s]]$cpu <- out[[s]]$cpu + (if (is.null(tk$cpu)) 0 else tk$cpu)
  }
  out
}

# Direction matching of the shear plan rebuild (R/shear_move.R), so a block
# record maps onto the new plan's coordinates as the rebuild mapped it.
shear_carry <- function(old_dirs, new_dirs, rec, p) {
  rec <- ns$.emc_shear_record(rec, p)
  out <- ns$.emc_shear_record(NULL, p)
  out$K <- matrix(0L, 2L * p, ncol(rec$K))
  scale <- c(rep(NA_real_, p), rep(1, p)); match <- c(rep(NA_integer_, p), p + seq_len(p))
  for (kk in seq_len(p)) {
    d <- new_dirs[, kk]
    cs <- abs(drop(crossprod(old_dirs, d))) / (sqrt(colSums(old_dirs^2)) * sqrt(sum(d^2)))
    m <- which.max(cs)
    if (cs[m] >= 0.9) { match[kk] <- m; scale[kk] <- sqrt(sum(old_dirs[, m]^2) / sum(d^2)) }
  }
  got <- which(!is.na(match))
  out$E[got] <- rec$E[match[got]] * scale[got]^2
  out$n[got] <- rec$n[match[got]]
  out$moved[got] <- rec$moved[match[got]]
  out$N[got] <- rec$N[match[got]] * scale[got]^2
  out$K[got, ] <- rec$K[match[got], ]
  out
}
sub_rec <- function(a, b, sign = 1) {  # a + sign * b, K padded
  w <- max(ncol(a$K), ncol(b$K))
  out <- a
  for (f in setdiff(names(a), "K")) if (!is.null(b[[f]])) out[[f]] <- a[[f]] + sign * b[[f]]
  out$K <- pad(a$K, w) + sign * pad(b$K, w)
  out
}

boot_state <- function(idx_by_chain, what) {
  emc <- post$emc
  for (j in seq_along(emc)) {
    tr <- traces[[j]]
    ch <- emc[[j]]
    p <- ch$n_pars; n <- ch$n_subjects
    m <- length(tr$cost)
    act <- seq_len(m)
    ib <- idx_by_chain[[j]]
    pick <- function(field) if (field %in% what || "all" %in% what) ib else act
    w_rest <- chain_inputs(tr, pick("rest"), p, n)
    w_mov <- chain_inputs(tr, pick("moves"), p, n)
    w_cost <- chain_inputs(tr, pick("costs"), p, n)
    w <- w_rest
    w$P <- w_mov$P; w$nP <- w_mov$nP
    w$cpu_part <- w_cost$cpu_part * ifelse(w_cost$nP > 0, w_mov$nP / pmax(1, w_cost$nP), 0)
    w$cpu_rest <- w_cost$cpu_rest * w_rest$n / max(1, w_cost$n)
    w$fam <- w_cost$fam
    w$cpu_hmc <- w_cost$cpu_hmc; w$hmc_att <- w_cost$hmc_att
    ch$alloc$window <- w
    # HMC evidence: the block's transitions replaced by the resampled ones.
    if (!is.null(ch$subject_hmc$evidence)) {
      a <- hmc_incs(tr, act, n, p); b <- hmc_incs(tr, pick("moves"), n, p)
      for (s in seq_len(n)) {
        ev <- ch$subject_hmc$evidence[[s]]
        if (is.null(ev)) next
        ev$n <- max(0, ev$n - a[[s]]$n + b[[s]]$n)
        ev$M <- ev$M - a[[s]]$M + b[[s]]$M
        ev$cpu <- max(0, ev$cpu - a[[s]]$cpu + b[[s]]$cpu)
        ch$subject_hmc$evidence[[s]] <- ev
      }
    }
    # Shear and collapsed records: the block's contribution replaced.
    incs_s <- Map(function(a, b) diff_rec(a, b), c(list(NULL), lapply(tr$cost[-m], `[[`, "shear_block")),
                  lapply(tr$cost, `[[`, "shear_block"))
    if (is.list(ch$shear$plan) && length(incs_s) && !is.null(incs_s[[1]]) &&
        identical(ch$shear$plan$stage, "sample")) {
      old_dirs <- emc0[[j]]$shear$plan$dirs
      blk_a <- shear_carry(old_dirs, ch$shear$plan$dirs, sum_rec(incs_s, act), p)
      blk_b <- shear_carry(old_dirs, ch$shear$plan$dirs, sum_rec(incs_s, pick("moves")), p)
      ch$shear$plan$last <- sub_rec(sub_rec(ns$.emc_shear_record(ch$shear$plan$last, p), blk_a, -1), blk_b, 1)
    }
    incs_c <- Map(function(a, b) diff_rec(a, b), c(list(NULL), lapply(tr$cost[-m], `[[`, "coll_block")),
                  lapply(tr$cost, `[[`, "coll_block"))
    if (is.list(ch$collapsed$plan) && length(incs_c) && !is.null(incs_c[[1]]) &&
        identical(ch$collapsed$plan$stage, "sample") && is.list(ch$collapsed$plan$last) &&
        ch$collapsed$plan$last$n >= 5) {
      blk_a <- sum_rec(incs_c, act); blk_b <- sum_rec(incs_c, pick("moves"))
      lst <- ns$.emc_collapsed_record(ch$collapsed$plan$last, p)
      ch$collapsed$plan$last <- sub_rec(sub_rec(lst, blk_a, -1), blk_b, 1)
    }
    emc[[j]] <- ch
  }
  emc
}

# Identity resampling reproduces the actual inputs.
ident <- lapply(traces, function(tr) seq_along(tr$cost))
chk <- solve_plan(boot_state(ident, "all"))
cat("identity reconstruction reproduces the plan per chain:",
    vapply(seq_along(chk), function(j) isTRUE(all.equal(chk[[j]]$rounds, replay[[j]]$rounds)) &&
             isTRUE(all.equal(chk[[j]]$rate, replay[[j]]$rate)) &&
             isTRUE(all.equal(chk[[j]]$tau, replay[[j]]$tau, tolerance = 1e-6)), logical(1)), "\n")

mbb <- function(m) {  # moving-block bootstrap indices
  starts <- sample.int(max(1L, m - L_block + 1L), ceiling(m / L_block), replace = TRUE)
  unlist(lapply(starts, function(s) s:min(m, s + L_block - 1L)))[seq_len(m)]
}
# The resampling draws from its own stream: solve_plan resets the global
# seed to the rebuild's, which would otherwise repeat every resample.
set.seed(1000 + k)
boot_seed <- get(".Random.seed", envir = globalenv())
res <- list(base = summ(replay))
for (v in variants) {
  reps <- lapply(seq_len(n_rep), function(b) {
    assign(".Random.seed", boot_seed, envir = globalenv())
    idx <- lapply(traces, function(tr) mbb(length(tr$cost)))
    boot_seed <<- get(".Random.seed", envir = globalenv())
    s <- summ(solve_plan(boot_state(idx, v)))
    s$rep <- b
    s
  })
  res[[v]] <- do.call(rbind, reps)
  cat("variant", v, "done\n")
}
res$same <- same
saveRDS(res, file.path(run_dir, sprintf("boot_%02d.rds", k)))
base <- res$base
for (v in variants) {
  r <- res[[v]]
  cat("\n==", v, "\n")
  print(do.call(rbind, lapply(split(r, r$chain), function(x) {
    b <- base[base$chain == x$chain[1], ]
    data.frame(chain = b$chain, tau = b$tau, tau_sd = sd(x$tau),
               mu = b$rounds_mu, mu_q = paste(round(quantile(x$rounds_mu, c(.1, .9)), 2), collapse = "-"),
               lsd = b$rounds_lsd, lsd_q = paste(round(quantile(x$rounds_lsd, c(.1, .9)), 2), collapse = "-"),
               coll = b$coll, coll_q = paste(round(quantile(x$coll, c(.1, .9)), 2), collapse = "-"),
               hmc_n = b$hmc_n, hmc_q = paste(round(quantile(x$hmc_n, c(.1, .9))), collapse = "-"),
               regime_flip = mean(x$regime != b$regime), q_n = b$q_n,
               q_q = paste(round(quantile(x$q_n, c(.1, .9))), collapse = "-"),
               bind_same = mean(x$bind == b$bind))
  })), digits = 3)
}
