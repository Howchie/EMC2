pkgload::load_all(".", quiet = TRUE)

archive_root <- Sys.getenv(
  "EMC2_FULLC_ROOT",
  ".Rtmp/claude-1000/-data-work-EMC2-dev-oo/3baace33-74ca-4333-ba0b-800d175b2295/scratchpad/runs/final")
lba_path <- Sys.getenv("EMC2_LBA_RDS", "/tmp/lba_sampling_audit/snapshot.rds")
artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
n_proposals <- as.integer(Sys.getenv("EMC2_MASK_PROPOSALS", "64"))
cases <- list(
  T30 = list(path = file.path(archive_root, "T30_fullc_1_emc.rds"), chain = 1L),
  corr15 = list(path = file.path(archive_root, "corr15_fullc_1_emc.rds"), chain = 1L),
  Forstmann = list(path = file.path(archive_root, "forstmann_fullc_2_emc.rds"), chain = 2L),
  LBA = list(path = lba_path, chain = 1L))

with_binding <- function(name, value, code) {
  ns <- asNamespace("EMC2")
  old <- get(name, envir = ns, inherits = FALSE)
  unlockBinding(name, ns)
  assign(name, value(old), envir = ns)
  lockBinding(name, ns)
  on.exit({
    unlockBinding(name, ns)
    assign(name, old, envir = ns)
    lockBinding(name, ns)
  }, add = TRUE)
  force(code)
}

load_chain <- function(case) {
  x <- EMC2:::restore_custom_kernel_pointers(
    EMC2:::restore_duplicates(readRDS(case$path)), quiet = TRUE)
  x[[case$chain]]
}

# Reconstruct the selector used before the package narrowed its production
# mask to t0. Retained only to reproduce the historical all-bound arm.
legacy_all_hard_fixed_idx <- function(sampler) {
  model <- EMC2:::.oo_model_list(sampler$model)
  bounds <- model$bound$minmax
  transform <- model$transform
  constrained <- character()
  if (is.matrix(bounds) && nrow(bounds) == 2L && !is.null(colnames(bounds))) {
    pars <- colnames(bounds)
    func <- stats::setNames(rep("identity", length(pars)), pars)
    if (length(transform$func)) {
      use <- intersect(pars, names(transform$func))
      func[use] <- transform$func[use]
    }
    lower <- stats::setNames(rep(-Inf, length(pars)), pars)
    upper <- stats::setNames(rep(Inf, length(pars)), pars)
    lower[func %in% c("exp", "pnorm")] <- 0
    upper[func == "pnorm"] <- 1
    if (length(transform$lower)) {
      use <- intersect(pars, names(transform$lower))
      lower[use] <- transform$lower[use]
    }
    if (length(transform$upper)) {
      use <- intersect(pars, names(transform$upper))
      upper[use] <- transform$upper[use]
    }
    constrained <- pars[bounds[1L, pars] > lower | bounds[2L, pars] < upper]
  }
  has_t0 <- "t0" %in% names(model$p_types) ||
    "t0" %in% names(transform$func)
  fixed_pars <- unique(c(constrained, if (has_t0) "t0"))
  if (!length(fixed_pars)) return(integer())
  which(EMC2:::get_p_types(sampler$par_names) %in% fixed_pars)
}

measure_arm <- function(chain, fixed_idx, arm, seed) {
  n <- chain$n_subjects
  i <- chain$samples$idx
  alpha <- matrix(chain$samples$alpha[, , i, drop = FALSE], chain$n_pars)
  ctx_seconds <- state_seconds <- numeric(6L)
  for (r in seq_along(ctx_seconds)) {
    ctx_started <- proc.time()[["elapsed"]]
    next_ctx <- EMC2:::.emc_collapsed_context(chain, fixed_idx,
      alpha[fixed_idx, , drop = FALSE])
    ctx_seconds[r] <- proc.time()[["elapsed"]] - ctx_started
    state_started <- proc.time()[["elapsed"]]
    next_cur <- EMC2:::.emc_shear_state(chain$samples$theta_mu[, i],
      as.matrix(chain$samples$theta_var[, , i]), next_ctx$fits,
      chain$prior, next_ctx$spec)
    state_seconds[r] <- proc.time()[["elapsed"]] - state_started
    if (r == 1L) {
      ctx <- next_ctx
      cur <- next_cur
    }
  }
  ctx_seconds <- median(ctx_seconds)
  state_seconds <- median(state_seconds)
  stopifnot(!is.null(cur), is.finite(cur$lp))
  x0 <- list(cur = cur, alpha = alpha,
    ll = chain$samples$subj_ll[, i], rounds = 0L, wpool = NULL,
    wpool_part = NULL, moved = FALSE)
  plan <- chain$collapsed$plan
  phi_jumps <- vector("list", n_proposals)
  proposed_rows <- vector("list", n_proposals)
  inner_acceptance <- numeric(n_proposals)
  exact_called <- logical(n_proposals)
  proposal_i <- 0L

  with_binding("collapsed_trajectories_cpp", function(old) {
    function(phi0, C, eps, n_traj, len, model) {
      out <- old(phi0, C, eps, n_traj, len, model)
      proposal_i <<- proposal_i + 1L
      phi_jumps[[proposal_i]] <<- as.numeric(out$phi - phi0)
      inner_acceptance[proposal_i] <<- mean(out$acc)
      out
    }
  }, {
    with_binding(".emc_wpool_subject_ll", function(old) {
      function(pool, rows, sampler, part = NULL, source) {
        if (identical(source, "collapsed")) {
          exact_called[proposal_i] <<- TRUE
          proposed_rows[[proposal_i]] <<- rows
        }
        old(pool, rows, sampler, part, source)
      }
    }, {
      started <- proc.time()[["elapsed"]]
      acc <- jump <- numeric(n_proposals)
      moved <- logical(n_proposals)
      for (k in seq_len(n_proposals)) {
        set.seed(seed + k - 1L)
        x <- x0
        x$moved <- FALSE
        x <- EMC2:::.emc_collapsed_round(x, plan, chain, ctx)
        acc[k] <- x$acc
        moved[k] <- isTRUE(x$moved)
        jump[k] <- sum(forwardsolve(plan$C, phi_jumps[[k]])^2)
      }
      elapsed <- proc.time()[["elapsed"]] - started
      support <- matrix(NA_integer_, n_proposals, 2L,
        dimnames = list(NULL, c("global_bound", "t0_rt")))
      model <- EMC2:::.oo_model_list(chain$model)
      for (k in which(exact_called)) {
        rows <- proposed_rows[[k]]
        bad_bound <- bad_t0 <- FALSE
        for (s in seq_along(rows)) {
          pars <- EMC2:::get_pars_oo(rows[[s]], chain$data[[s]], model)
          if (length(model$bound$minmax))
            bad_bound <- bad_bound || any(!EMC2:::do_bound(pars, model$bound))
          if ("t0" %in% colnames(pars) && "rt" %in% names(chain$data[[s]]))
            bad_t0 <- bad_t0 || any(pars[, "t0"] >= chain$data[[s]]$rt,
              na.rm = TRUE)
        }
        support[k, ] <- c(as.integer(bad_bound), as.integer(bad_t0))
      }
      movement <- sum(acc * jump)
      data.frame(arm = arm, subjects = n, group_parameters = chain$n_pars,
        fixed_coefficients = length(fixed_idx), proposals = n_proposals,
        expected_acceptance = mean(acc), moved_rate = mean(moved),
        mean_inner_acceptance = mean(inner_acceptance),
        inner_movement_rate = mean(vapply(phi_jumps, function(z)
          any(z != 0), logical(1))),
        exact_eval_rate = mean(exact_called),
        chart_esjd_per_second = movement /
          (elapsed + n_proposals * (ctx_seconds + state_seconds)),
        round_seconds = elapsed / n_proposals,
        context_and_state_seconds = ctx_seconds + state_seconds,
        bound_hit_rate = mean(support[, "global_bound"], na.rm = TRUE),
        t0_hit_rate = mean(support[, "t0_rt"], na.rm = TRUE))
    })
  })
}

result <- list()
row <- 0L
for (case_name in names(cases)) {
  chain <- load_chain(cases[[case_name]])
  type <- EMC2:::get_p_types(chain$par_names)
  t0_idx <- which(type == "t0")
  auto_idx <- legacy_all_hard_fixed_idx(chain)
  masks <- list(unconditioned = integer(), t0_only = t0_idx,
    all_hard_bounds = auto_idx)
  for (arm in names(masks)) {
    row <- row + 1L
    result[[row]] <- measure_arm(chain, masks[[arm]], arm,
      seed = 20261010L + 100L * match(case_name, names(cases)))
    result[[row]]$case <- case_name
    cat(case_name, arm, "acceptance", result[[row]]$expected_acceptance,
      "ESJD/s", result[[row]]$chart_esjd_per_second,
      "round ms", 1000 * result[[row]]$round_seconds, "\n")
    flush.console()
  }
}
result <- do.call(rbind, result)
write.csv(result, file.path(artifact_dir,
  "conditional_mask_efficiency_historical.csv"), row.names = FALSE)
print(result, row.names = FALSE)
