library(EMC2)

snapshot <- "/tmp/lba_sampling_audit/snapshot.rds"
fit <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS(snapshot)), quiet = TRUE)

records <- list()
for (j in seq_along(fit)) {
  ch <- fit[[j]]
  obj <- readRDS(sprintf("/tmp/lba_sampling_audit/proposals_chain%d.rds", j))
  bounds <- ch$model()$bound$minmax
  exceptions <- ch$model()$bound$exception
  for (k in seq_along(obj$props)) {
    proposal_bad <- list()
    subject_bad <- list()
    support_bad <- integer()
    for (s in seq_len(ch$n_subjects)) {
      alpha <- matrix(obj$props[[k]][, s], 1L,
                      dimnames = list(NULL, ch$par_names))
      pars <- EMC2:::get_pars_oo(alpha, ch$data[[s]], ch$model)
      names_here <- intersect(colnames(bounds), colnames(pars))
      bad <- sweep(pars[, names_here, drop = FALSE], 2L,
                   bounds[1L, names_here], `<`) |
        sweep(pars[, names_here, drop = FALSE], 2L,
              bounds[2L, names_here], `>`)
      for (name in intersect(names(exceptions), names_here))
        bad[, name] <- bad[, name] & pars[, name] != exceptions[name]
      for (name in names_here) {
        if (is.null(proposal_bad[[name]])) proposal_bad[[name]] <- FALSE
        if (is.null(subject_bad[[name]])) subject_bad[[name]] <- 0L
        subject_crossing <- any(bad[, name])
        proposal_bad[[name]] <- proposal_bad[[name]] || subject_crossing
        subject_bad[[name]] <- subject_bad[[name]] + as.integer(subject_crossing)
      }
      rt <- ch$data[[s]]$rt
      support_bad[s] <- any(is.finite(rt) & pars[, "t0"] >= rt)
    }
    for (name in names(proposal_bad)) {
      records[[length(records) + 1L]] <- data.frame(
        chain = j, proposal = k, scale = obj$scales[k],
        parameter = name, any_invalid = proposal_bad[[name]],
        subjects_invalid = subject_bad[[name]],
        t0_rt_support = if (name == "t0") sum(support_bad) else 0L)
    }
  }
}

detail <- do.call(rbind, records)
summary <- aggregate(cbind(any_invalid, subjects_invalid, t0_rt_support) ~
                       scale + parameter, detail, mean)
write.csv(detail,
          "plans/lba_sampling_audit_2026-10-09/bounds_by_parameter_detail.csv",
          row.names = FALSE)
write.csv(summary,
          "plans/lba_sampling_audit_2026-10-09/bounds_by_parameter_summary.csv",
          row.names = FALSE)
print(summary[order(summary$scale, -summary$any_invalid), ], row.names = FALSE)
