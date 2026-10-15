# Offline allocator decisions from identical measurements; no transitions run.
args <- commandArgs(TRUE)
root <- normalizePath(args[1])
scratch <- file.path(root, ".Rtmp/claude-1000/-data-work-EMC2-dev-oo/19c4c6df-25fb-4858-8ef0-fdc4d9b6acea/scratchpad")
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB", file.path(root, ".Rlib_alloc_review"))))
ns <- asNamespace("EMC2")
outdir <- file.path(root, "plans/alloc_prediction_2026-10-15/eval_operator")
env <- new.env(parent = ns)
sys.source(file.path(outdir, "evaluated_allocation.R"), env)
old <- new.env(parent = ns)
baseline <- system2("git", c("-C", shQuote(root), "show", "bf677e14:R/move_allocation.R"), stdout = TRUE)
eval(parse(text = baseline), envir = old)
variants <- strsplit(Sys.getenv("VARIANTS", "baseline,point,point_no_credit,kappa,kappa_no_credit"), ",", fixed = TRUE)[[1]]
source(file.path(outdir, "subject_objective.R"))
source(file.path(outdir, "conditional_subjects.R"))
files <- if (length(args) > 1L) args[-1L] else c("pair_new_2", "pair_new_3", "pair_point_2", "pair_point_3", "arm_pointsave_1")
records <- list()
for (name in files) {
  if (name == "checkpoint") {
    box <- new.env()
    load("/data/work/PM/NirvanaHons_Nback/Samples/samples_control_Exp1_frqk7r_blk_gL_sepA_sepL_sepV_sv_t0T_mt.RData", envir = box)
    emc <- ns$restore_duplicates(box$emc)
    stage <- "sample"
  } else {
    emc <- ns$restore_duplicates(readRDS(file.path(scratch, paste0(name, "_emc.rds"))))
    stage <- emc[[1]]$samples$stage[emc[[1]]$samples$idx]
  }
  Vbar <- lapply(emc, function(ch) {
    k <- env$.emc_alloc_draws(ch, max(2L, ch$alloc$window$n))
    apply(ch$samples$theta_var[, , k, drop = FALSE], 1:2, mean)
  })
  hst <- lapply(seq_along(emc), function(j) ns$.emc_hmc_plan_chain(emc[[j]], stage, Vbar[[j]]))
  cache <- new.env()
  ms <- lapply(seq_along(emc), function(j) env$.emc_alloc_measure(emc, j, hst[[j]], Vbar[[j]], stage, cache))
  # Match the original lower-bound factors as well as its base rhos.
  cache_old <- new.env()
  ms_old <- lapply(seq_along(emc), function(j) old$.emc_alloc_measure(emc, j, hst[[j]], Vbar[[j]], stage, cache_old))
  ms_cond <- if ("conditional" %in% variants) conditional_measure(emc, ms, hst, env)
  for (variant in variants) {
    ev <- new.env(parent = ns)
    sys.source(file.path(outdir, if (variant == "latest") "latest_allocation.R" else "evaluated_allocation.R"), ev)
    if (grepl("no_credit", variant) || variant == "conditional") ev$.emc_alloc_dev_ratio <- function(...) 1
    inputs <- lapply(seq_along(emc), function(j) {
      m <- if (variant == "baseline") ms_old[[j]] else ms[[j]]
      if (variant == "conditional") m <- ms_cond[[j]]
      if (!grepl("kappa", variant)) m$g_kappa <- NULL
      ev$.emc_alloc_units(m, hst[[j]], m$shear, FALSE, 100L)
    })
    policy_solve <- if (variant == "baseline") old$.emc_alloc_solve_units else ev$.emc_alloc_solve_units
    if (variant %in% c("own_subject", "baseline_subject")) {
      mode <- if (variant == "own_subject") "own" else "baseline"
      policy_solve <- function(u) subject_objective(u, ev, mode)
    }
    solutions <- lapply(inputs, policy_solve)
    keys <- Reduce(intersect, lapply(solutions, function(s) names(s$options)))
    choice <- if (length(keys)) keys[which.min(Reduce(`+`, lapply(solutions, function(s) s$options[keys])))] else NULL
    if (variant == "latest" && length(keys)) choice <- ev$.emc_alloc_choose(
      Reduce(`+`, lapply(solutions, function(s) s$options[keys])),
      Reduce(`+`, lapply(solutions, function(s) s$option_cost[keys])))
    if (!is.null(choice)) for (j in seq_along(inputs)) {
      inputs[[j]]$regime_fixed <- choice
      solutions[[j]] <- policy_solve(inputs[[j]])
    }
    ng <- length(ms[[1]]$g_rho)
    subjects <- ng + seq_len(emc[[1]]$n_subjects)
    q <- do.call(pmax, lapply(solutions, function(s) s$q[subjects]))
    for (j in seq_along(inputs)) {
      if (any(solutions[[j]]$q[subjects] < q - 1e-9)) {
        for (s in seq_along(subjects)) {
          i <- subjects[s]
          if (!is.null(inputs[[j]]$units[[i]]$p) && !is.null(inputs[[j]]$units[[i]]$h))
            inputs[[j]]$units[[i]]$q_fixed <- q[s]
        }
        solutions[[j]] <- policy_solve(inputs[[j]])
      }
      sol <- solutions[[j]]
      rates <- ev$.emc_alloc_plan_rates(inputs[[j]]$units, sol)
      rates$r[subjects] <- pmin(1, sol$x[subjects] + sol$regime)
      row <- data.frame(state = name, variant = variant, chain = j, choice = sol$choice,
                        tau = sol$tau, cost = sol$cost, objective = sol$objective,
                        q = mean(q), hmc = mean(rates$r[subjects]), shear = sum(sol$x[seq_len(ng)]),
                        collapsed = sol$joint)
      records[[length(records) + 1L]] <- row
      print(row, row.names = FALSE)
    }
    saveRDS(list(inputs = inputs, solutions = solutions, measurements = ms),
            file.path(outdir, paste0("replay_", name, "_", variant, ".rds")))
    write.csv(do.call(rbind, records), file.path(outdir, Sys.getenv("REPLAY_CSV", "replay.csv")), row.names = FALSE)
  }
}
