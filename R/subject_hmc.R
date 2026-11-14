# Subject-conditional Hamiltonian Monte Carlo.
#
# The last update of an iteration, after the subject, ensemble and shear
# steps: a subject may take one HMC transition on its conditional posterior
# with the final group state (mu, V) held fixed:
#   log pi_s(a) = ll_s(a) - (a - mu)' V^-1 (a - mu) / 2 + const,
# in the sampler's own alpha coordinates. The particle step moves locally or
# jumps between fitted Gaussians; neither follows a curved ridge of a
# subject's likelihood. A trajectory does. Coming last, its rest-of-iteration
# measure (R/move_allocation.R) is exactly everything else the iteration does.
#
# Kernel. The likelihood fit lik_laplace (precision H_s, positive
# semi-definite, see .emc_lik_canonical) and the current group prior give
# Q = H_s + V^-1, positive definite whenever V is. With Q = R'R, the whitened
# coordinates are z = R (a - mu), so a = mu + S z with S = R^-1 and
# S S' = Q^-1. Leapfrog runs in z with unit mass. The likelihood gradient is
# a forward difference along the columns of S at a fixed scale h, the p
# probes and the new centre in one likelihood call (a central difference
# would double the rows for an error leapfrog already exceeds). A deterministic force
# keeps leapfrog volume preserving and reversible, so the exact Hamiltonian
# acceptance ratio corrects it (Neal 2011, section 5.5.3): the finite
# differences change efficiency, not the target. Every position the
# trajectory evaluates, the start included, must have finite likelihoods above
# the bound floor at the centre and all probes; otherwise the whole trajectory
# is rejected. That rule depends on positions only, and a path and its
# reverse visit the same positions, so it keeps the kernel reversible. The
# step size eps_s is jittered uniformly by +-10% independently of the state,
# and L = ceiling((pi/2) / eps) steps make a quarter period of the whitened
# Gaussian approximation.
#
# Random numbers come from the subject's own stream, after its particle step,
# in this order: the selection uniform (drawn on the master, only when the
# rate is in (0, 1)), the step-size jitter, p momenta, the acceptance uniform.
# A subject with rate 0 draws nothing.
#
# Adaptation. Rates come from the move allocation (R/move_allocation.R).
# Rates, step sizes and the likelihood metric change only at the proposal
# rebuilds (every block in burn and adapt, at 100, 400, 1600, ... production
# transitions) and are fixed in between: each block runs a fixed kernel that
# leaves the subject conditional invariant. The exception is a subject's
# first burn block, whose draws are discarded: every eligible subject then
# moves at an even rate while its step size adapts after each transition
# (Robbins-Monro on log eps towards an acceptance target), until a fixed
# number of transitions, so the first allocation reads measured moves rather
# than a short, transient window.
# Later, step sizes move by acceptance-weighted squared jump per likelihood
# row, compared between the two halves of the jitter. A subject first
# eligible outside burn has its step size bracketed with one-step trials
# (Hoffman and Gelman 2014, algorithm 4) and pilot transitions on private
# copies of its state; they inform the allocation and never replace a draw.

# Bounds on the step size in whitened units: above 2 leapfrog is unstable on
# the Gaussian approximation, and the lower bound caps a trajectory at 79
# steps. Every adapted quantity (eps, rates) stays in a fixed compact set.
.EMC_HMC_EPS_RANGE <- c(0.02, 2)
# A subject's first burn block runs at an even rate while its step size is
# found by stochastic approximation towards this acceptance rate.
.EMC_HMC_WARMUP_RATE <- 0.5
.EMC_HMC_WARMUP_ACCEPT <- 0.75
# Transitions after which a subject's warmup stops for the rest of the
# block: its step size has settled, and the allocation has the evidence it
# needs.
.EMC_HMC_WARMUP_MOVES <- 20L
# Finite-difference step of the likelihood gradient, pilot trajectories run
# for a subject that lacks evidence, and the moves that count as evidence.
.EMC_HMC_H <- 1e-4
.EMC_HMC_PILOTS <- 20L
.EMC_HMC_MIN_EVIDENCE <- 10L

.emc_hmc_clip_eps <- function(eps) {
  if (!is.finite(eps)) return(eps)
  min(.EMC_HMC_EPS_RANGE[2], max(.EMC_HMC_EPS_RANGE[1], eps))
}

# Hierarchies whose subject conditional is the Gaussian prior times the
# subject likelihood, with one likelihood per subject.
.emc_hmc_support <- function(sampler) {
  why <- if (!isTRUE(sampler$type %in% c("standard", "diagonal-gamma"))) {
    paste("hierarchy type", sampler$type)
  } else if (any(sampler$nuisance)) {
    "nuisance parameters"
  } else if (!is.null(sampler$gd) ||
             (!is.null(sampler$group_designs) && length(sampler$group_designs))) {
    "group design"
  } else if (!is.data.frame(sampler$data[[1L]])) {
    "joint likelihood"
  } else if (length(unique(attr(sampler$data, "shared_ll_idx"))) > 1L) {
    "joint likelihood"
  } else if (is.null(sampler$par_names) || length(sampler$par_names) != sampler$n_pars) {
    "parameter names"
  } else NULL
  list(ok = is.null(why), reason = if (is.null(why)) "" else why)
}

# A likelihood value at or below the bound floor (every trial at min_ll) marks
# a position outside the parameter bounds.
.emc_hmc_ll_floor <- function(dadm) {
  n <- length(attr(dadm, "expand"))
  if (!n) return(-Inf)
  f <- n * log(1e-10)
  f + 1e-9 * abs(f)
}

# --- kernel ---------------------------------------------------------------

.emc_hmc_frame <- function(prec, mu, V) {
  p <- length(mu)
  if (!is.matrix(prec) || nrow(prec) != p || ncol(prec) != p ||
      !all(is.finite(prec)) || !is.matrix(V) || nrow(V) != p) return(NULL)
  Vc <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(Vc)) return(NULL)
  Vi <- chol2inv(Vc)
  Q <- prec + Vi
  R <- tryCatch(chol((Q + t(Q)) / 2), error = function(e) NULL)
  if (is.null(R) || !all(is.finite(R))) return(NULL)
  list(mu = as.numeric(mu), Vi = Vi, R = R, S = backsolve(R, diag(p)), p = p)
}

# Log target, whitened gradient and likelihood at z. ll_centre, when known,
# is the centre's likelihood and saves its evaluation.
.emc_hmc_eval <- function(z, fr, ll_fn, h, floor, ll_centre = NULL) {
  p <- fr$p
  a <- fr$mu + drop(fr$S %*% z)
  X <- rbind(if (is.null(ll_centre)) a, matrix(a, p, p, byrow = TRUE) + t(fr$S) * h)
  ll <- ll_fn(X)
  rows <- nrow(X)
  if (length(ll) != rows || anyNA(ll) || !all(is.finite(ll)) || any(ll <= floor) ||
      (!is.null(ll_centre) && !(is.finite(ll_centre) && ll_centre > floor)))
    return(list(ok = FALSE, rows = rows))
  off <- if (is.null(ll_centre)) 1L else 0L
  lc <- if (is.null(ll_centre)) ll[1L] else ll_centre
  d <- a - fr$mu
  Vid <- drop(fr$Vi %*% d)
  g <- (ll[off + seq_len(p)] - lc) / h - drop(crossprod(fr$S, Vid))
  list(ok = TRUE, a = a, z = z, ll = lc, lp = lc - 0.5 * sum(d * Vid), g = g,
       rows = rows)
}

# One HMC transition from a0 (likelihood ll0). Draws, in order: the jitter,
# p momenta, the acceptance uniform; always all of them.
.emc_hmc_transition <- function(a0, ll0, fr, ll_fn, eps_base, h, floor) {
  p <- fr$p
  U <- stats::runif(1L, 0.9, 1.1)
  v0 <- stats::rnorm(p)
  u_acc <- stats::runif(1L)
  eps <- eps_base * U
  L <- max(1L, as.integer(ceiling((pi / 2) / eps)))
  out <- list(alpha = a0, ll = ll0, accepted = FALSE, acc_prob = 0,
              disp = numeric(p), U = U, eps = eps, L = L, steps = 0L,
              rows = 0L, status = "abort", energy = NA_real_)
  cur <- .emc_hmc_eval(drop(fr$R %*% (a0 - fr$mu)), fr, ll_fn, h, floor,
                       ll_centre = ll0)
  out$rows <- cur$rows
  if (!cur$ok) return(out)
  v <- v0 + 0.5 * eps * cur$g
  z <- cur$z
  for (l in seq_len(L)) {
    z <- z + eps * v
    new <- .emc_hmc_eval(z, fr, ll_fn, h, floor)
    out$rows <- out$rows + new$rows
    out$steps <- l
    if (!new$ok) return(out)
    v <- v + (if (l < L) eps else 0.5 * eps) * new$g
  }
  energy <- (-new$lp + 0.5 * sum(v^2)) - (-cur$lp + 0.5 * sum(v0^2))
  out$energy <- energy
  out$status <- "ok"
  out$acc_prob <- if (is.finite(energy)) min(1, exp(-energy)) else 0
  out$disp <- new$a - a0
  if (u_acc < out$acc_prob) {
    out$alpha <- new$a
    out$ll <- new$ll
    out$accepted <- TRUE
  }
  out
}

# Step size where a single leapfrog step's acceptance crosses one half,
# doubling or halving from eps0 (Hoffman and Gelman 2014, algorithm 4).
.emc_hmc_bracket <- function(a0, ll0, fr, ll_fn, h, floor, eps0, max_trials = 12L) {
  rows <- 0L
  cur <- .emc_hmc_eval(drop(fr$R %*% (a0 - fr$mu)), fr, ll_fn, h, floor,
                       ll_centre = ll0)
  rows <- rows + cur$rows
  if (!cur$ok) return(list(eps = NA_real_, rows = rows, trials = 0L))
  one <- function(eps) {
    v0 <- stats::rnorm(fr$p)
    v <- v0 + 0.5 * eps * cur$g
    new <- .emc_hmc_eval(cur$z + eps * v, fr, ll_fn, h, floor)
    rows <<- rows + new$rows
    if (!new$ok) return(0)
    v <- v + 0.5 * eps * new$g
    e <- (-new$lp + 0.5 * sum(v^2)) - (-cur$lp + 0.5 * sum(v0^2))
    if (is.finite(e)) min(1, exp(-e)) else 0
  }
  eps <- eps0
  acc <- one(eps)
  dir <- if (acc > 0.5) 1 else -1
  trials <- 1L
  while (trials < max_trials) {
    nxt <- eps * 2^dir
    acc <- one(nxt)
    trials <- trials + 1L
    if ((dir > 0 && acc <= 0.5) || (dir < 0 && acc > 0.5)) {
      eps <- if (dir > 0) eps else nxt
      break
    }
    eps <- nxt
  }
  list(eps = eps, rows = rows, trials = trials)
}

# --- measurement ----------------------------------------------------------

# Whitening matrix of the reference frame: rows are the eigenvectors of
# Sigma_ref = (H + Vbar^-1)^-1 divided by their sd, each signed so its largest
# entry is positive.
.emc_hmc_ref_whiten <- function(prec, Vbar) {
  p <- nrow(Vbar)
  Vi <- tryCatch(chol2inv(chol(Vbar)), error = function(e) NULL)
  if (is.null(Vi) || !is.matrix(prec) || !all(is.finite(prec))) return(NULL)
  Q <- prec + Vi
  e <- tryCatch(eigen((Q + t(Q)) / 2, symmetric = TRUE), error = function(e) NULL)
  if (is.null(e) || any(!is.finite(e$values)) || min(e$values) <= 0) return(NULL)
  E <- e$vectors
  sgn <- apply(E, 2L, function(v) { k <- which.max(abs(v)); if (v[k] < 0) -1 else 1 })
  E <- sweep(E, 2L, sgn, "*")
  # Sigma_ref = E diag(1 / values) E', so its whitening is diag(sqrt(values)) E'.
  sqrt(e$values) * t(E)
}

# --- the round ------------------------------------------------------------

# One transition of subject s from (a0, ll0) given the group state.
.emc_hmc_run <- function(a0, ll0, data, model, group_mu, group_var, plan_s,
                         r_cores = 1L, par_names,
                         ll_cache_validated = FALSE) {
  started <- proc.time()
  task <- list(attempted = TRUE, rows = 0L, status = "frame")
  fr <- .emc_hmc_frame(plan_s$prec, group_mu, as.matrix(group_var))
  if (!is.null(fr)) {
    ll_fn <- function(X) {
      colnames(X) <- par_names
      tryCatch(as.numeric(calc_ll_manager(X, dadm = data, model = model,
                                          r_cores = r_cores,
                                          ll_cache_validated = ll_cache_validated)),
        error = function(e) {
          cls <- .emc_reject_record(e, "particle")
          if (identical(.emc_failure_action(cls), "abort"))
            .emc_failure_abort(cls, "particle", e)
          NULL
        })
    }
    tr <- .emc_hmc_transition(a0, ll0, fr, ll_fn, plan_s$eps, plan_s$h, plan_s$floor)
    task[c("alpha", "ll", "accepted", "acc_prob", "disp", "U", "eps", "L", "steps",
           "rows", "status", "energy")] <-
      tr[c("alpha", "ll", "accepted", "acc_prob", "disp", "U", "eps", "L", "steps",
           "rows", "status", "energy")]
    task$score <- if (is.null(plan_s$W)) NA_real_ else
      tr$acc_prob * mean(drop(plan_s$W %*% tr$disp)^2)
  }
  spent <- proc.time() - started
  task$cpu <- unname(spent[["user.self"]] + spent[["sys.self"]])
  task
}

# Worker side of a round: each listed subject from its own stream.
.emc_hmc_compute <- function(msg, ctx) .emc_with_preserved_rng({
  tasks <- vector("list", length(msg$subs))
  seeds <- vector("list", length(msg$subs))
  native <- .emc_native_hmc(msg, ctx)
  for (k in seq_along(msg$subs)) {
    if (isTRUE(native$done[k])) {
      tasks[[k]] <- native$tasks[[k]]
      seeds[[k]] <- native$seeds[[k]]
      next
    }
    s <- msg$subs[k]
    assign(".Random.seed", msg$seeds[[k]], envir = globalenv())
    plan_s <- ctx$hmc[[s]]
    plan_s$eps <- msg$eps[k]
    data_s <- .emc_pooled_set_one(ctx$data[[s]], msg$pooled)
    tasks[[k]] <- .emc_hmc_run(msg$alpha[, k], msg$ll[k], data_s, ctx$model,
                               msg$mu[, k], msg$var, plan_s, ctx$r_cores,
                               ctx$par_names,
                               ll_cache_validated = isTRUE(ctx$ll_cache_validated))
    seeds[[k]] <- get(".Random.seed", envir = globalenv())
  }
  list(tasks = tasks, seeds = seeds)
})

# The transitions of a request in one native call (src/native_particle.cpp),
# for subjects whose likelihood arguments are already prepared. NULL, or
# `done` FALSE for a subject, leaves the R path above to run it from the
# stream the request carries.
.emc_native_hmc <- function(msg, ctx) {
  if (!isTRUE(getOption("emc2.native_particle_block", TRUE)) ||
      !isTRUE(ctx$ll_cache_validated) || !is.null(msg$pooled) ||
      !is.null(ctx$r_cores) && ctx$r_cores > 1L || !length(msg$subs) ||
      !is.character(ctx$par_names) || !is.data.frame(ctx$data[[msg$subs[1L]]]))
    return(NULL)
  tryCatch(
    emc2_native_hmc_subjects(msg$alpha, msg$ll, msg$mu, msg$var, msg$seeds,
                             msg$eps, ctx$hmc[msg$subs], ctx$data[msg$subs],
                             ctx$model, ctx$par_names),
    error = function(e) NULL)
}

# The round, from the final state of the other updates. Draws each
# subject's selection uniform from its own stream on the master, runs the
# selected transitions on the workers (or in-process), and returns the
# updated proposals and streams with one task per subject (NULL when no
# transition ran). pick: per subject, whether the iteration's joint draw
# with the particle step (R/sampling.R) selected HMC (NA where the round
# draws its own selection).
.emc_hmc_round <- function(plan, proposals, pars, type, streams, data, model,
                           wpool = NULL, n_cores = 1L, r_cores = 1L, par_names,
                           pick = NULL, ll_cache_validated = FALSE) {
  n <- length(streams)
  p <- nrow(proposals) - 1L
  tasks <- vector("list", n)
  go <- logical(n)
  .emc_with_preserved_rng(for (s in seq_len(n)) {
    ps <- plan[[s]]
    if (is.null(ps) || !(ps$rate > 0)) next
    if (ps$rate >= 1) {
      go[s] <- TRUE
      next
    }
    if (length(pick) >= s && !is.na(pick[s])) {
      go[s] <- pick[s]
      next
    }
    assign(".Random.seed", streams[[s]], envir = globalenv())
    go[s] <- stats::runif(1L) < ps$rate
    streams[[s]] <- get(".Random.seed", envir = globalenv())
  })
  subs <- which(go)
  if (!length(subs))
    return(list(proposals = proposals, streams = streams, tasks = tasks, pool = wpool,
                plan = plan))
  gl <- lapply(subs, function(s) get_group_level(pars, s, type))
  # Workers hold the data of their template, so the current pooled values
  # travel with each request.
  pooled <- attr(data[[1L]], "pooled_values")
  msg_for <- function(ss) {
    k <- match(ss, subs)
    list(kind = "hmc", subs = ss, notify = FALSE, pooled = pooled,
         alpha = proposals[seq_len(p), ss, drop = FALSE], ll = proposals[p + 1L, ss],
         mu = matrix(vapply(gl[k], function(g) as.numeric(g$mu), numeric(p)), p),
         var = as.matrix(gl[[1L]]$var), seeds = streams[ss],
         eps = vapply(ss, function(s) plan[[s]]$eps, numeric(1)))
  }
  results <- vector("list", 0L)
  if (!is.null(wpool) && isTRUE(wpool$alive) && wpool$n > 1L) {
    part <- .emc_lpt_partition(.emc_subject_cost(data)[subs], min(wpool$n, length(subs)))
    part <- lapply(part, function(i) subs[i])
    sent <- logical(length(part))
    for (w in seq_along(part)) {
      if (!length(part[[w]])) next
      sent[w] <- .emc_wpool_request(wpool, w, msg_for(part[[w]]))
      if (!sent[w]) {
        .emc_wpool_degraded(.emc_wpool_transport_error("send to a worker failed"))
        wpool$alive <- FALSE
        break
      }
    }
    deadline <- .emc_wpool_deadline()
    for (w in seq_along(part)) {
      if (!length(part[[w]])) next
      res <- if (sent[w]) .emc_wpool_reply(wpool, w, deadline) else NULL
      if (is.null(res) || !is.null(res$failed) || length(res$tasks) != length(part[[w]])) {
        if (sent[w]) {
          .emc_wpool_degraded("subject HMC worker gave no valid reply")
          wpool$alive <- FALSE
        }
        res <- .emc_hmc_compute(msg_for(part[[w]]), list(
          data = data, model = model, hmc = plan, r_cores = r_cores,
          par_names = par_names, ll_cache_validated = ll_cache_validated))
      }
      results[[length(results) + 1L]] <- list(subs = part[[w]], res = res)
    }
  } else {
    ctx <- list(data = data, model = model, hmc = plan, r_cores = r_cores,
                par_names = par_names, ll_cache_validated = ll_cache_validated)
    hmc_cores <- max(1L, min(n_cores, length(subs)))
    if (hmc_cores <= 1L) {
      # Each subject runs from its own stream, so one request serves all.
      results[[1L]] <- list(subs = subs, res = .emc_hmc_compute(msg_for(subs), ctx))
    } else {
      out <- auto_mclapply(subs, function(s) .emc_hmc_compute(msg_for(s), ctx),
                           mc.cores = hmc_cores, mc.set.seed = FALSE)
      for (k in seq_along(subs))
        results[[k]] <- list(subs = subs[k], res = out[[k]])
    }
  }
  for (r in results) for (k in seq_along(r$subs)) {
    s <- r$subs[k]
    tk <- r$res$tasks[[k]]
    tasks[[s]] <- tk
    streams[[s]] <- r$res$seeds[[k]]
    if (isTRUE(tk$accepted)) {
      proposals[seq_len(p), s] <- tk$alpha
      proposals[p + 1L, s] <- tk$ll
    }
    if (isTRUE(plan[[s]]$warmup)) {
      plan[[s]] <- .emc_hmc_warmup_step(plan[[s]], tk)
      tasks[[s]]$eps_next <- plan[[s]]$eps
    }
  }
  list(proposals = proposals, streams = streams, tasks = tasks, pool = wpool, plan = plan)
}

# One warmup update of a subject's step size from a transition's acceptance
# probability (an aborted trajectory counts as zero).
.emc_hmc_warmup_step <- function(ps, tk) {
  acc <- if (is.null(tk$acc_prob) || !is.finite(tk$acc_prob)) 0 else tk$acc_prob
  ps$n_warm <- ps$n_warm + 1L
  ps$eps <- .emc_hmc_clip_eps(ps$eps * exp((acc - .EMC_HMC_WARMUP_ACCEPT) / sqrt(ps$n_warm)))
  if (ps$n_warm >= .EMC_HMC_WARMUP_MOVES) {
    ps$warmup <- FALSE
    ps$rate <- 0
  }
  ps
}

# --- plan state -----------------------------------------------------------

.emc_hmc_active <- function(sampler) {
  !is.null(sampler$subject_hmc) && isTRUE(sampler$subject_hmc$supported)
}

# Each subject's HMC rate under a worker plan (zero without one).
.emc_hmc_plan_rates <- function(plan, n) {
  if (is.null(plan)) return(numeric(n))
  vapply(seq_len(n), function(s)
    if (length(plan) >= s && !is.null(plan[[s]])) plan[[s]]$rate else 0, numeric(1))
}

# What a transition needs from the plan (sent to the workers once per block);
# NULL when no subject has a positive rate.
.emc_hmc_worker_plan <- function(sampler) {
  if (!.emc_hmc_active(sampler)) return(NULL)
  st <- sampler$subject_hmc
  plan <- lapply(st$subjects, function(x) {
    if (!isTRUE(x$eligible) || !(x$rate > 0)) return(NULL)
    x[c("prec", "eps", "rate", "h", "floor", "W", "warmup", "n_warm")]
  })
  if (all(vapply(plan, is.null, logical(1)))) NULL else plan
}

.emc_hmc_new_evidence <- function(p, key) {
  list(key = key, n = 0L, M = matrix(0, p, p), accepted = 0L, aborts = 0L,
       rows = 0, cpu = 0, U = numeric(0), score = numeric(0), arows = numeric(0),
       accp = numeric(0))
}

.emc_hmc_decay_evidence <- function(ev) {
  if (is.null(ev)) return(ev)
  d <- .EMC_ALLOC_EVIDENCE_DECAY
  ev$n <- ev$n * d
  ev$M <- ev$M * d
  ev$cpu <- ev$cpu * d
  ev
}

.emc_hmc_add_evidence <- function(ev, task) {
  ev$n <- ev$n + 1L
  ev$M <- ev$M + task$acc_prob * tcrossprod(task$disp)
  ev$accepted <- ev$accepted + as.integer(isTRUE(task$accepted))
  ev$aborts <- ev$aborts + as.integer(!identical(task$status, "ok"))
  ev$rows <- ev$rows + task$rows
  ev$cpu <- ev$cpu + if (is.null(task$cpu)) 0 else task$cpu
  ev$U <- c(ev$U, task$U)
  ev$score <- c(ev$score, task$score)
  ev$arows <- c(ev$arows, task$rows)
  ev$accp <- c(ev$accp, task$acc_prob)
  ev
}

# Add one round's transitions to each subject's evidence and the totals.
.emc_hmc_record_tasks <- function(st, tasks) {
  if (is.null(st) || !isTRUE(st$supported)) return(st)
  if (is.list(tasks)) {
    r <- emc_hmc_record_tasks_cpp(st, tasks)
    if (!is.null(r)) return(r)
  }
  for (s in seq_along(tasks)) {
    tk <- tasks[[s]]
    if (is.null(tk)) next
    st$totals$attempts <- st$totals$attempts + 1L
    st$totals$rows <- st$totals$rows + tk$rows
    st$totals$cpu <- st$totals$cpu + tk$cpu
    if (is.null(tk$disp)) next
    st$totals$accepted <- st$totals$accepted + as.integer(isTRUE(tk$accepted))
    st$totals$aborts <- st$totals$aborts + as.integer(!identical(tk$status, "ok"))
    ev <- st$evidence[[s]]
    if (is.null(ev)) ev <- .emc_hmc_new_evidence(
      st$p, .emc_hmc_key(st$subjects[[s]]$eps))
    ev <- .emc_hmc_add_evidence(ev, tk)
    # Warmup moves are kept as the evidence of the step size they reach.
    if (!is.null(tk$eps_next)) {
      st$subjects[[s]]$eps <- tk$eps_next
      ev$key <- .emc_hmc_key(tk$eps_next)
    }
    st$evidence[[s]] <- ev
  }
  st
}

# HMC measurements are comparable while the step size is unchanged; the
# jumps are kept in alpha coordinates, so a new likelihood metric only
# changes the frame they are read in.
.emc_hmc_key <- function(eps) signif(eps, 12)

# A chain's HMC state at a rebuild: per subject, eligibility (a usable
# likelihood precision and frame), the reference whitening W (the frame its
# moves are measured in), the bound floor, and a step size moved by the
# last evidence; rates are set by the move allocation.
.emc_hmc_plan_chain <- function(sampler, stage, Vbar) {
  sup <- .emc_hmc_support(sampler)
  old <- sampler$subject_hmc
  n <- sampler$n_subjects
  p <- sampler$n_pars
  st <- list(supported = sup$ok, reason = sup$reason,
             version = if (is.null(old$version)) 1L else old$version + 1L,
             stage = stage, p = p, n = n, h = .EMC_HMC_H)
  if (!sup$ok) return(st)
  lik <- sampler$lik_laplace
  evidence <- if (length(old$evidence) == n) old$evidence else vector("list", n)
  # Production evidence decays at each rebuild like every move family's, so
  # a subject whose HMC stopped is measured afresh when it runs again.
  if (identical(stage, "sample") && identical(old$stage, "sample"))
    evidence <- lapply(evidence, .emc_hmc_decay_evidence)
  subjects <- vector("list", n)
  for (s in seq_len(n)) {
    fit <- if (length(lik) >= s) lik[[s]] else NULL
    prec <- if (is.list(fit)) fit$prec else NULL
    ok_prec <- is.matrix(prec) && nrow(prec) == p && ncol(prec) == p &&
      all(is.finite(prec)) && max(abs(prec - t(prec))) <= 1e-8 * max(1, abs(prec))
    prior <- tryCatch(.emc_burn_subject_prior(sampler, s), error = function(e) NULL)
    fr <- if (ok_prec && !is.null(prior)) .emc_hmc_frame(prec, prior$mu, prior$var)
    W <- if (ok_prec) .emc_hmc_ref_whiten(prec, Vbar)
    reason <- if (!is.list(fit)) "no likelihood fit" else if (!ok_prec)
      "invalid likelihood precision" else if (is.null(fr) || is.null(W))
        "factorisation failed" else ""
    old_s <- if (length(old$subjects) == n) old$subjects[[s]] else NULL
    sj <- list(eligible = !nzchar(reason), reason = reason,
               prec = if (ok_prec) prec, W = W, h = .EMC_HMC_H,
               floor = .emc_hmc_ll_floor(sampler$data[[s]]), rate = 0,
               eps = if (!is.null(old_s$eps)) old_s$eps else NA_real_,
               n_tuned = if (is.null(old_s$n_tuned)) 0L else old_s$n_tuned,
               calib_cpu = 0, warmup = FALSE, n_warm = 0L)
    if (identical(stage, "burn") && sj$eligible && is.na(sj$eps)) {
      sj$warmup <- TRUE
      sj$eps <- .emc_hmc_clip_eps(p^(-1 / 4))
    } else {
      sj <- .emc_hmc_update_eps(sj, evidence[[s]])
    }
    ev <- evidence[[s]]
    if (!is.null(ev) && (is.na(sj$eps) || !identical(ev$key, .emc_hmc_key(sj$eps))))
      evidence[s] <- list(NULL)
    subjects[[s]] <- sj
  }
  st$subjects <- subjects
  st$evidence <- evidence
  st$totals <- if (is.null(old$totals)) list(attempts = 0L, accepted = 0L,
    aborts = 0L, rows = 0, cpu = 0, calib_cpu = 0) else old$totals
  st$history <- old$history
  st
}

# Bounded log-step update from the two halves of the jitter: a higher score
# (acceptance-weighted squared jump per likelihood row) on one side moves eps
# towards it. A window that cannot resolve the slope keeps eps; one where
# trajectories are mostly rejected halves it.
.emc_hmc_update_eps <- function(sj, ev) {
  if (is.null(ev) || is.na(sj$eps) || ev$n < 4L) return(sj)
  if (mean(ev$accp) < 0.2) {
    sj$eps <- .emc_hmc_clip_eps(sj$eps / 2)
    return(sj)
  }
  eff <- ev$score / pmax(1, ev$arows)
  hi <- ev$U > 1
  if (sum(hi) < 3L || sum(!hi) < 3L || anyNA(eff)) return(sj)
  d <- mean(eff[hi]) - mean(eff[!hi])
  se <- sqrt(stats::var(eff[hi]) / sum(hi) + stats::var(eff[!hi]) / sum(!hi))
  if (!is.finite(se) || se <= 0 || abs(d) < 2 * se) return(sj)
  sj$n_tuned <- sj$n_tuned + 1L
  sj$eps <- .emc_hmc_clip_eps(sj$eps * exp(sign(d) * min(log(2), log(1.5) / sqrt(sj$n_tuned))))
  sj
}

# Calibration of one subject on a private copy of its state with the group
# fixed: a bracketed step size, then pilot transitions. Returns the evidence,
# the CPU spent and the advanced RNG stream.
.emc_hmc_calibrate <- function(sampler, st, s, bracket = TRUE, pilots = 0L) {
  started <- proc.time()
  assign(".Random.seed", sampler$rng$subjects[[s]], envir = globalenv())
  sj <- st$subjects[[s]]
  prior <- .emc_burn_subject_prior(sampler, s)
  data <- sampler$data[[s]]
  pn <- sampler$par_names
  # Failures are tallied and returned for the master to count: this runs in a
  # forked worker, whose counters and warnings end with it.
  failed <- list(counts = stats::setNames(integer(length(.EMC_REJECT_CLASSES)),
                                          .EMC_REJECT_CLASSES), msgs = list())
  ll_fn <- function(X) {
    colnames(X) <- pn
    tryCatch(as.numeric(calc_ll_manager(X, dadm = data, model = sampler$model,
                                        r_cores = 1L)), error = function(e) {
      cls <- .emc_classify_failure(e)
      if (identical(.emc_failure_action(cls), "abort"))
        .emc_failure_abort(cls, "pilot", e)
      failed$counts[[cls]] <<- failed$counts[[cls]] + 1L
      if (is.null(failed$msgs[[cls]]))
        failed$msgs[[cls]] <<- tryCatch(conditionMessage(e), error = function(e2) "")
      NULL
    })
  }
  fr <- .emc_hmc_frame(sj$prec, prior$mu, prior$var)
  a <- as.numeric(sampler$samples$alpha[, s, sampler$samples$idx])
  ll <- sampler$samples$subj_ll[s, sampler$samples$idx]
  rows <- 0
  eps <- sj$eps
  if (!is.null(fr) && bracket) {
    b <- .emc_hmc_bracket(a, ll, fr, ll_fn, sj$h, sj$floor,
                          if (is.finite(eps)) eps else st$p^(-1 / 4))
    rows <- rows + b$rows
    if (is.finite(b$eps)) eps <- .emc_hmc_clip_eps(b$eps)
  }
  ev <- NULL
  if (!is.null(fr) && is.finite(eps) && pilots > 0L) {
    ev <- .emc_hmc_new_evidence(st$p, .emc_hmc_key(eps))
    for (k in seq_len(pilots)) {
      t0 <- proc.time()
      tr <- .emc_hmc_transition(a, ll, fr, ll_fn, eps, sj$h, sj$floor)
      t1 <- proc.time() - t0
      tr$cpu <- unname(t1[["user.self"]] + t1[["sys.self"]])
      rows <- rows + tr$rows
      tr$score <- tr$acc_prob * mean(drop(sj$W %*% tr$disp)^2)
      ev <- .emc_hmc_add_evidence(ev, tr)
      a <- tr$alpha
      ll <- tr$ll
    }
  }
  spent <- proc.time() - started
  list(eps = if (is.null(fr)) NA_real_ else eps, evidence = ev, rows = rows,
       cpu = unname(spent[["user.self"]] + spent[["sys.self"]]),
       seed = get(".Random.seed", envir = globalenv()), rejects = failed)
}

.emc_hmc_apply_calibration <- function(st, s, r) {
  sj <- st$subjects[[s]]
  sj$calib_cpu <- r$cpu
  st$totals$calib_cpu <- st$totals$calib_cpu + r$cpu
  if (!is.finite(r$eps)) {
    sj$eligible <- FALSE
    sj$reason <- "no usable step size"
  } else {
    if (!identical(sj$eps, r$eps)) {
      sj$n_tuned <- 0L
      if (!is.null(st$evidence[[s]]) &&
          !identical(st$evidence[[s]]$key, .emc_hmc_key(r$eps)))
        st$evidence[s] <- list(NULL)
    }
    sj$eps <- r$eps
    if (!is.null(r$evidence)) {
      ev <- st$evidence[[s]]
      st$evidence[[s]] <- if (is.null(ev)) r$evidence else .emc_hmc_merge_evidence(ev, r$evidence)
    }
  }
  st$subjects[[s]] <- sj
  st
}

.emc_hmc_merge_evidence <- function(a, b) {
  for (k in c("n", "M", "accepted", "aborts", "rows", "cpu"))
    a[[k]] <- a[[k]] + b[[k]]
  for (k in c("U", "score", "arows", "accp")) a[[k]] <- c(a[[k]], b[[k]])
  a
}
