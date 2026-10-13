# Instrumentation for the allocator audit (see README.md). Sourced by
# calibrate.R before fit(); every wrapper calls the package's own function
# and only records what passed through, so the sampler runs unchanged.
#
# Chain processes: each iteration's raw measurement increments are appended
# to sampler$ws1_trace (the jumps and costs the allocation window sums, the
# moves' cumulative block records, the HMC tasks), so a rebuild's inputs can
# be block-resampled later.
# Master: every production rebuild saves the state add_proposals() received
# (`pre_<k>.rds`, the block's trace included), and the allocator's per-chain
# measurements, units and records (`alloc_<k>.rds`).

ws1_install_trace <- function(out_dir, save_states = TRUE) {
  ns <- asNamespace("EMC2")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  wrap <- function(name, make) {
    orig <- get(name, envir = ns)
    f <- make(orig)
    environment(f) <- list2env(list(orig = orig, out_dir = out_dir, save_states = save_states),
                               parent = environment(orig))
    utils::assignInNamespace(name, f, ns)
  }
  push <- function(sampler, field, x) {
    tr <- sampler$ws1_trace
    if (is.null(tr)) tr <- list()
    tr[[field]] <- c(tr[[field]], list(x))
    sampler$ws1_trace <- tr
    sampler
  }
  environment(push) <- ns
  assign("ws1_push", push, envir = globalenv())

  wrap(".emc_alloc_record_group", function(orig) function(sampler, pars) {
    idx <- sampler$samples$idx
    d <- c(as.numeric(pars$tmu), .emc_alloc_lsd(pars$tvar)) -
      c(sampler$samples$theta_mu[, idx], .emc_alloc_lsd(sampler$samples$theta_var[, , idx]))
    ws1_push(orig(sampler, pars), "group", d)
  })
  wrap(".emc_alloc_record_subjects", function(orig) function(sampler, alpha, mu = NULL) {
    idx <- sampler$samples$idx
    old <- matrix(sampler$samples$alpha[, , idx], nrow(alpha))
    dmu <- if (is.null(mu)) numeric(nrow(alpha)) else
      as.numeric(mu) - as.numeric(sampler$samples$theta_mu[, idx])
    ws1_push(orig(sampler, alpha, mu), "subjects", list(d = alpha - old, dmu = dmu))
  })
  wrap(".emc_alloc_record_particle", function(orig) function(sampler, before, after, ran, cost, times) {
    ws1_push(orig(sampler, before, after, ran, cost, times), "particle",
             list(d = after - before, ran = ran, cost = cost, times = times))
  })
  wrap(".emc_alloc_record_cost", function(orig) function(sampler, rest, tasks, moves = list()) {
    out <- orig(sampler, rest, tasks, moves)
    hmc <- lapply(tasks, function(tk) if (is.null(tk)) NULL else
      list(disp = tk$disp, acc_prob = tk$acc_prob, cpu = tk$cpu, work = tk$work,
           cost = tk$cost))
    ws1_push(out, "cost", list(rest = rest, moves = moves, hmc = hmc,
      shear_block = out$shear$plan$block, coll_block = out$collapsed$plan$block))
  })

  assign("ws1_k", 0L, envir = globalenv())
  assign("ws1_depth", 0L, envir = globalenv())
  wrap("add_proposals", function(orig) function(emc, stage, n_cores, n_blocks, block_iter = 100L) {
    due <- stage != "preburn" && (stage != "sample" || .emc_sample_rebuild_due(emc) ||
                                    is.null(emc[[1L]]$sample_proposals_next))
    if (due) {
      k_rebuild <- get("ws1_k", globalenv()) + 1L
      assign("ws1_k", k_rebuild, envir = globalenv())
      seed <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
      if (isTRUE(save_states))
        saveRDS(list(emc = emc, stage = stage, n_cores = n_cores, n_blocks = n_blocks,
                     block_iter = block_iter, seed = seed),
                file.path(out_dir, sprintf("pre_%02d.rds", k_rebuild)), compress = FALSE)
      for (j in seq_along(emc)) emc[[j]]$ws1_trace <- NULL
      assign("ws1_capture", list(), envir = globalenv())
    }
    out <- orig(emc, stage, n_cores, n_blocks, block_iter)
    if (due) {
      k_rebuild <- get("ws1_k", globalenv())
      cap <- get("ws1_capture", envir = globalenv())
      cap$records <- lapply(out, function(ch) ch$alloc$last)
      cap$stage <- stage
      cap$transitions <- out[[1L]]$sample_transitions
      saveRDS(cap, file.path(out_dir, sprintf("alloc_%02d.rds", k_rebuild)))
    }
    out
  })
  capture <- function(field, x) {
    if (!exists("ws1_capture", globalenv())) return(invisible())
    cap <- get("ws1_capture", envir = globalenv())
    cap[[field]] <- c(cap[[field]], list(x))
    assign("ws1_capture", cap, envir = globalenv())
  }
  assign("ws1_capture_add", capture, envir = globalenv())
  wrap(".emc_alloc_rebuild", function(orig) function(emc, stage, n_cores, B) {
    ws1_capture_add("rebuild_input", list(stage = stage, n_cores = n_cores, B = B))
    orig(emc, stage, n_cores, B)
  })
  wrap(".emc_alloc_measure", function(orig) function(emc, j, hmc_state, Vbar, stage = NULL, ...) {
    ms <- orig(emc, j, hmc_state, Vbar, stage, ...)
    ws1_capture_add("measure", list(j = j, ms = ms, dirs = emc[[j]]$shear$plan$dirs,
                                    pool_G = emc[[j]]$alloc$pool))
    ms
  })
  wrap(".emc_alloc_units", function(orig) function(ms, hmc_state, shear_on, optimistic, B = NULL) {
    u <- orig(ms, hmc_state, shear_on, optimistic, B)
    ws1_capture_add("units", list(optimistic = optimistic, u = u))
    u
  })
  wrap(".emc_alloc_solve_units", function(orig) function(u) {
    depth <- get("ws1_depth", globalenv()) + 1L
    assign("ws1_depth", depth, envir = globalenv())
    on.exit(assign("ws1_depth", depth - 1L, envir = globalenv()))
    sol <- orig(u)
    if (depth == 1L) ws1_capture_add("solve", list(u = u, sol = sol))
    sol
  })
  invisible(TRUE)
}
