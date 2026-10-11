pkgload::load_all(".", quiet = TRUE, compile = FALSE)
x <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS("/tmp/lba_sampling_audit/snapshot.rds")),
  quiet = TRUE)
options(emc2.subject_hmc = FALSE)
for (j in seq_along(x)) {
  x[[j]]$production_options$emc2.subject_hmc <- FALSE
  x[[j]]$sample_proposals_next <- x[[j]]$sample_transitions
}

# Capture the independent per-chain optimizer inputs before the allocator
# synchronizes all chains to the largest selected collapsed rate.
.captured_collapsed_allocations <- list()
trace(".emc_alloc_solve_units", where = asNamespace("EMC2"), print = FALSE,
  exit = quote({
    .GlobalEnv$.captured_collapsed_allocations[[
      length(.GlobalEnv$.captured_collapsed_allocations) + 1L]] <-
      list(units = u, solution = returnValue())
  }))
x <- EMC2:::add_proposals(x, "sample", n_cores = 12L, n_blocks = 1L,
                          block_iter = 200L)
untrace(".emc_alloc_solve_units", where = asNamespace("EMC2"))
captured <- .captured_collapsed_allocations
stopifnot(length(captured) == length(x))

allocation_diagnostic <- do.call(rbind, lapply(seq_along(captured), function(j) {
  u <- captured[[j]]$units
  sol <- captured[[j]]$solution
  group_units <- u$units[seq_along(u$joint$idx)]
  slow <- which.max(vapply(group_units, function(z) max(z$rho), numeric(1)))
  data.frame(chain = j, selected_rounds = sol$joint,
    rest_cost = u$C_rest, calibration_cost = u$calibration_cost,
    move_cost = u$joint$cost, move_fixed_cost = u$joint$fixed_cost,
    slow_group_index = slow,
    slow_rho = max(group_units[[slow]]$rho),
    move_factor = u$joint$g[slow],
    move_saturation = u$joint$a[slow],
    objective = sol$objective,
    objective_without_move = sol$objective_without_joint)
}))
print(allocation_diagnostic, row.names = FALSE)
write.csv(allocation_diagnostic,
  "plans/collapsed_transport_adaptation_2026-10-09/allocator_current_diagnostic.csv",
  row.names = FALSE)
saveRDS(list(chains = x, captured = captured),
  "/tmp/lba_sampling_audit/shared_allocation_current.rds")

candidate_rates <- seq(0, 4, by = 0.125)
objective <- t(vapply(candidate_rates, function(rate) {
  vapply(captured, function(item) {
    u <- item$units
    units <- u$units
    if (rate > 0 && !is.null(u$joint)) {
      joint <- u$joint
      factor <- EMC2:::.emc_alloc_factor(joint$g, rate,
        if (is.null(joint$a)) joint$g else joint$a)
      for (i in seq_along(joint$idx))
        units[[joint$idx[i]]]$rho <- units[[joint$idx[i]]]$rho * factor[i]
    } else {
      rate <- 0
    }
    EMC2:::.emc_alloc_solve(units, u$C_rest + rate *
      if (rate > 0 && !is.null(u$joint)) u$joint$cost else 0)$objective
  }, numeric(1))
}, numeric(length(captured))))

cat("chain-local selected rates:", paste(vapply(captured, function(z)
  z$solution$joint, numeric(1)), collapse = ","), "\n")
cat("rate; chain objectives; worst-chain objective\n")
for (i in seq_along(candidate_rates))
  cat(candidate_rates[i], paste(signif(objective[i, ], 6), collapse = ","),
      max(objective[i, ]), "\n")
best <- which.min(apply(objective, 1L, max))
cat("minimax common rate:", candidate_rates[best], "\n")
cat("synchronized max-of-local rates:", max(vapply(captured, function(z)
  z$solution$joint, numeric(1))), "\n")
cat("actual plans after allocator synchronization:",
    paste(vapply(x, function(ch) ch$collapsed$plan$rounds, numeric(1)),
          collapse = ","), "\n")
