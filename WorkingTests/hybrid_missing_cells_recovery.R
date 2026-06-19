# Hybrid within/between-subjects recovery stress test
#
# This script is for manual use from the repo root, e.g.
#   source("WorkingTests/hybrid_missing_cells_recovery.R")
#   obj <- build_hybrid_missing_test()
#   fit <- run_hybrid_missing_fit(obj, iter = 300)
#
# The core question it targets is whether EMC2's hierarchical sampler can cope
# with subject-level parameters that correspond to design cells not observed for
# some subjects. EMC2 will still sample a full alpha vector for every subject;
# the unobserved cell parameters are then informed only through the hierarchical
# prior and, if allowed, cross-parameter covariance.

suppressPackageStartupMessages({
  library(EMC2)
})


make_hybrid_subject_table <- function(n_subj = 24, seed = 1) {
  set.seed(seed)
  stopifnot(n_subj %% 2 == 0)
  subjects <- factor(seq_len(n_subj))
  between_group <- rep(c("med", "high"), each = n_subj / 2)
  data.frame(
    subjects = subjects,
    between_group = factor(between_group, levels = c("med", "high"))
  )
}


build_hybrid_design <- function(subjects) {
  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)
  design(
    factors = list(
      subjects = subjects,
      S = c("left", "right"),
      Cond3 = c("low", "med", "high")
    ),
    Rlevels = c("left", "right"),
    functions = list(Correct = function(d) ifelse(d$lM==TRUE,1,0),
                     Error = function(d) ifelse(d$lM==FALSE,1,0)),
    matchfun = matchfun,
    model = RDM,
    formula = list(
      v ~ 0+Error+Correct:Cond3,
      B ~ Cond3,
      A ~ 1,
      t0 ~ 1,
      s ~ 1
    ),
    constants = c(
      s = log(1), A=log(0)
    ),
    report_p_vector = TRUE
  )
}


make_true_group_means <- function(design_obj) {
  p <- sampled_pars(design_obj, doMap = FALSE)
  p[] <- c(
    v_Error = log(.6),
    "v_Correct:Cond3low" = log(1.6),
    "v_Correct:Cond3med" = log(1.2),
    "v_Correct:Cond3high" = log(0.75),
    B = log(0.9),
    B_Cond3med = log(1.05) - log(0.9),
    B_Cond3high = log(1.2) - log(0.9),
    #A = log(0.35),
    t0 = log(0.28)
  )
  p
}


make_true_covariance <- function(design_obj) {
  par_names <- names(sampled_pars(design_obj, doMap = FALSE))
  sds <- c(
    v_Error = 0.12,
    "v_Correct:Cond3low" = 0.20,
    "v_Correct:Cond3med" = 0.20,
    "v_Correct:Cond3high" = 0.20,
    B = 0.10,
    B_Cond3med = 0.10,
    B_Cond3high = 0.10,
    A = 0.08,
    t0 = 0.08,
    s = 0.05
  )
  Sigma <- diag(sds[par_names]^2)
  dimnames(Sigma) <- list(par_names, par_names)

  Sigma["B_Cond3med", "B_Cond3high"] <- 0.006
  Sigma["B_Cond3high", "B_Cond3med"] <- 0.006
  Sigma["B", "t0"] <- -0.003
  Sigma["t0", "B"] <- -0.003
  Sigma["v_Correct:Cond3low", "B_Cond3med"] <- 0.004
  Sigma["v_Correct:Cond3med", "B_Cond3high"] <- 0.004
  Sigma["v_Correct:Cond3high", "B_Cond3med"] <- 0.004
  Sigma["v_Correct:Cond3low", "B_Cond3high"] <- 0.004
  Sigma["B_Cond3med", "v_Correct:Cond3low"] <- 0.004
  Sigma["B_Cond3med", "v_Correct:Cond3high"] <- 0.004
  Sigma["B_Cond3high", "v_Correct:Cond3low"] <- 0.004
  Sigma["B_Cond3high", "v_Correct:Cond3med"] <- 0.004
  Sigma["v_Correct:Cond3low", "v_Correct:Cond3med"] <- 0.010
  Sigma["v_Correct:Cond3low", "v_Correct:Cond3high"] <- 0.010
  Sigma["v_Correct:Cond3med", "v_Correct:Cond3low"] <- 0.010
  Sigma["v_Correct:Cond3high", "v_Correct:Cond3med"] <- 0.010
  Sigma["v_Correct:Cond3high", "v_Correct:Cond3low"] <- 0.010
  Sigma["v_Correct:Cond3med", "v_Correct:Cond3high"] <- 0.010
  Sigma
}


apply_hybrid_missingness <- function(data, subject_table) {
  merged <- merge(data, subject_table, by = "subjects", all.x = TRUE, sort = FALSE)
  keep <- merged$Cond3 == "low" |
    (merged$between_group == "med" & merged$Cond3 == "med") |
    (merged$between_group == "high" & merged$Cond3 == "high")
  out <- merged[keep, , drop = FALSE]
  out$Cond3 <- factor(out$Cond3, levels = c("low", "med", "high"))
  out$subjects <- factor(out$subjects, levels = levels(subject_table$subjects))
  out[order(out$subjects, out$trials), , drop = FALSE]
}


apply_blockwise_exclusions <- function(data,
                                       drop_low_subjects = integer(0),
                                       drop_between_subjects = integer(0)) {
  out <- data
  if (length(drop_low_subjects) > 0) {
    out <- out[!(out$subjects %in% drop_low_subjects & out$Cond3 == "low"), , drop = FALSE]
  }
  if (length(drop_between_subjects) > 0) {
    is_between_cell <- (out$between_group == "med" & out$Cond3 == "med") |
      (out$between_group == "high" & out$Cond3 == "high")
    out <- out[!(out$subjects %in% drop_between_subjects & is_between_cell), , drop = FALSE]
  }
  out$Cond3 <- factor(out$Cond3, levels = c("low", "med", "high"))
  out$subjects <- factor(out$subjects, levels = levels(data$subjects))
  out[order(out$subjects, out$trials), , drop = FALSE]
}


summarise_subject_cells <- function(data) {
  tabs <- xtabs(~ subjects + Cond3, data = data)
  present <- (tabs > 0) * 1L
  counts <- rowSums(present)
  list(
    trial_counts = tabs,
    cell_presence = present,
    n_cells_per_subject = counts
  )
}


inspect_subject_designs <- function(data, design_obj) {
  dadm <- EMC2:::design_model(
    data = data,
    design = design_obj,
    compress = FALSE,
    verbose = FALSE
  )
  by_subj <- EMC2:::dm_list(dadm)
  out <- lapply(by_subj, function(x) {
    b_dm_full <- attr(x, "designs")[["B"]]
    b_dm <- b_dm_full[attr(b_dm_full, "expand"), , drop = FALSE]
    data.frame(
      parameter = colnames(b_dm),
      column_sum = colSums(b_dm),
      column_nonzero = colSums(abs(b_dm) > 0)
    )
  })
  out
}


build_hybrid_missing_test <- function(n_subj = 24,
                                      n_trials = 160,
                                      seed = 1,
                                      exclusion_seed = 2,
                                      n_drop_low = 0,
                                      n_drop_between = 0,
                                      compress = TRUE) {
  subject_table <- make_hybrid_subject_table(n_subj = n_subj, seed = seed)
  design_obj <- build_hybrid_design(subject_table$subjects)
  true_mu <- make_true_group_means(design_obj)
  true_Sigma <- make_true_covariance(design_obj)

  set.seed(seed)
  true_alpha <- make_random_effects(
    design = design_obj,
    group_means = true_mu,
    n_subj = n_subj,
    covariances = true_Sigma
  )
  rownames(true_alpha) <- as.character(subject_table$subjects)

  full_data <- make_data(
    parameters = true_alpha,
    design = design_obj,
    n_trials = n_trials
  )
  full_data$Cond3 <- factor(full_data$Cond3, levels = c("low", "med", "high"))
  full_data$subjects <- factor(full_data$subjects, levels = levels(subject_table$subjects))

  hybrid_data <- apply_hybrid_missingness(full_data, subject_table)

  set.seed(exclusion_seed)
  drop_low_subjects <- if (n_drop_low > 0) {
    sample(levels(subject_table$subjects), n_drop_low)
  } else {
    character(0)
  }
  remaining <- setdiff(levels(subject_table$subjects), drop_low_subjects)
  drop_between_subjects <- if (n_drop_between > 0) {
    sample(remaining, n_drop_between)
  } else {
    character(0)
  }

  stressed_data <- apply_blockwise_exclusions(
    hybrid_data,
    drop_low_subjects = drop_low_subjects,
    drop_between_subjects = drop_between_subjects
  )

  prior_obj <- prior(design_obj, type = "standard")
  emc_obj <- make_emc(
    data = stressed_data,
    design = design_obj,
    type = "standard",
    prior_list = prior_obj,
    compress = compress,
    n_chains = 3
  )

  list(
    design = design_obj,
    prior = prior_obj,
    emc = emc_obj,
    subject_table = subject_table,
    true_mu = true_mu,
    true_Sigma = true_Sigma,
    true_alpha = true_alpha,
    full_data = full_data,
    hybrid_data = hybrid_data,
    stressed_data = stressed_data,
    dropped = list(
      low = drop_low_subjects,
      between = drop_between_subjects
    ),
    cell_summary = summarise_subject_cells(stressed_data),
    subject_designs = inspect_subject_designs(stressed_data, design_obj)
  )
}


run_hybrid_missing_fit <- function(obj,
                                   iter = 400,
                                   particle_factor = 20,
                                   cores_per_chain = 1,
                                   cores_for_chains = 3,
                                   step_size = 100,
                                   verbose = TRUE) {
  fit(
    obj$emc,
    iter = iter,
    particle_factor = particle_factor,
    cores_per_chain = cores_per_chain,
    cores_for_chains = cores_for_chains,
    step_size = step_size,
    verbose = verbose
  )
}


make_hybrid_group_subject_table <- function(n_subj = 24, seed = 11) {
  set.seed(seed)
  subjects <- factor(seq_len(n_subj))
  load_group <- sample(c("Low", "Med"), size = n_subj, replace = TRUE)
  data.frame(
    subjects = subjects,
    LoadGroup = factor(load_group, levels = c("Low", "Med"))
  )
}


make_hybrid_group_skeleton <- function(subject_table) {
  rows <- lapply(seq_len(nrow(subject_table)), function(i) {
    cur_sub <- subject_table$subjects[i]
    cur_group <- as.character(subject_table$LoadGroup[i])
    expand.grid(
      subjects = cur_sub,
      S = c("left", "right"),
      Load = factor(c(cur_group, "High"), levels = c("Low", "Med", "High"))
    )
  })
  out <- do.call(rbind, rows)
  out$subjects <- factor(out$subjects, levels = levels(subject_table$subjects))
  out$S <- factor(out$S, levels = c("left", "right"))
  out$Load <- factor(out$Load, levels = c("Low", "Med", "High"))
  out$R <- factor("left", levels = c("left", "right"))
  out[order(out$subjects), , drop = FALSE]
}


build_hybrid_group_subject_design <- function(subjects) {
  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)
  design(
    factors = list(
      subjects = 1,
      S = c("left", "right"),
      Load = c("Low", "Med", "High")
    ),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    functions = list(
      HighTag = function(d) factor(
        ifelse(as.character(d$Load) == "High", "High", "Base"),
        levels = c("Base", "High")
      ),
      Correct = function(d) ifelse(d$lM==TRUE,1,0),
      Error = function(d) ifelse(d$lM==FALSE,1,0)
    ),
    model = RDM,
    formula = list(
      v ~ 0+Error+(Correct:HighTag),
      B ~ HighTag,
      A ~ 1,
      t0 ~ 1,
      s ~ 1
    ),
    constants = c(
      s = log(1)
    ),
    report_p_vector = TRUE
  )
}


build_hybrid_group_group_design <- function(subject_data, subject_design) {
  group_design(
    formula = list(
      B ~ LoadGroup,
      `v_Correct:HighTagBase` ~ LoadGroup
    ),
    data = subject_data,
    subject_design = subject_design
  )
}


make_hybrid_group_beta <- function(subject_design, group_design_obj) {
  beta <- sampled_pars(subject_design, group_design = group_design_obj, doMap = FALSE)
  beta[] <- c(
    v_Error = log(.6),
    "v_Correct:HighTagBase" = log(1.25),
    "v_Correct:HighTagBase_LoadGroupMed" = log(1.05) - log(1.25),
    "v_Correct:HighTagHigh" = log(.6),
    B = log(0.88),
    B_LoadGroupMed = log(1.00) - log(0.88),
    B_HighTagHigh = log(1.18) - log(0.88),
    A = log(0.35),
    t0 = log(0.28)
  )
  beta
}


make_hybrid_group_covariance <- function(subject_design) {
  par_names <- names(sampled_pars(subject_design, doMap = FALSE))
  sds <- c(
    v_Error = 0.12,
    "v_Correct:HighTagBase" = 0.20,
    "v_Correct:HighTagHigh" = 0.20,
    B = 0.10,
    B_HighTagHigh = 0.10,
    A = 0.08,
    t0 = 0.08,
    s = 0.05
  )
  Sigma <- diag(sds[par_names]^2)
  dimnames(Sigma) <- list(par_names, par_names)
  Sigma["B", "t0"] <- -0.003
  Sigma["t0", "B"] <- -0.003
  Sigma["v_Correct:HighTagBase", "B_HighTagHigh"] <- 0.004
  Sigma["v_Correct:HighTagHigh", "B_HighTagHigh"] <- 0.004
  Sigma["B_HighTagHigh", "v_Correct:HighTagBase"] <- 0.004
  Sigma["B_HighTagHigh", "v_Correct:HighTagHigh"] <- 0.004
  Sigma
}


draw_hybrid_group_alpha <- function(beta, Sigma, subject_design, group_design_obj, subjects) {
  beta_mat <- matrix(beta, ncol = 1, dimnames = list(names(beta), NULL))
  Sigma_cube <- array(
    Sigma,
    dim = c(nrow(Sigma), ncol(Sigma), 1),
    dimnames = list(rownames(Sigma), colnames(Sigma), NULL)
  )
  alpha_arr <- EMC2:::get_alphas(
    mu = beta_mat,
    var = Sigma_cube,
    sub_names = as.character(subjects),
    group_design = group_design_obj
  )
  alpha <- t(alpha_arr[, , 1, drop = FALSE][, , 1])
  alpha
}


build_hybrid_group_design_test <- function(n_subj = 24,
                                           n_trials = 160,
                                           seed = 11,
                                           compress = TRUE) {
  subject_table <- make_hybrid_group_subject_table(n_subj = n_subj, seed = seed)
  skeleton <- make_hybrid_group_skeleton(subject_table)
  subject_design <- build_hybrid_group_subject_design(subject_table$subjects)
  group_design_obj <- build_hybrid_group_group_design(subject_table, subject_design)
  true_beta <- make_hybrid_group_beta(subject_design, group_design_obj)
  true_Sigma <- make_hybrid_group_covariance(subject_design)

  set.seed(seed)
  true_alpha <- draw_hybrid_group_alpha(
    beta = true_beta,
    Sigma = true_Sigma,
    subject_design = subject_design,
    group_design_obj = group_design_obj,
    subjects = subject_table$subjects
  )

  sim_data <- make_data(
    parameters = true_alpha,
    design = subject_design,
    data = skeleton,
    expand = n_trials,
    conditional_on_data = FALSE
  )
  sim_data$subjects <- factor(sim_data$subjects, levels = levels(subject_table$subjects))
  sim_data$Load <- factor(sim_data$Load, levels = c("Low", "Med", "High"))
  sim_data$LoadGroup <- subject_table$LoadGroup[match(sim_data$subjects, subject_table$subjects)]

  prior_obj <- prior(subject_design, type = "standard", group_design = group_design_obj)
  emc_obj <- make_emc(
    data = sim_data,
    design = subject_design,
    type = "standard",
    prior_list = prior_obj,
    group_design = group_design_obj,
    compress = compress,
    n_chains = 3
  )

  list(
    design = subject_design,
    group_design = group_design_obj,
    prior = prior_obj,
    emc = emc_obj,
    subject_table = subject_table,
    skeleton = skeleton,
    true_beta = true_beta,
    true_Sigma = true_Sigma,
    true_alpha = true_alpha,
    data = sim_data,
    cell_summary = summarise_subject_cells(transform(sim_data, Cond3 = Load))
  )
}


print_hybrid_group_design_summary <- function(obj) {
  cat("\nLoadGroup by subject:\n")
  print(obj$subject_table)
  cat("\nObserved Load cells per subject:\n")
  print(obj$cell_summary$cell_presence)
  cat("\nTrue beta coefficients:\n")
  print(obj$true_beta)
  cat("\nGroup design summary:\n")
  print(summary(obj$group_design))
  invisible(obj)
}


compare_to_truth <- function(fit_obj, obj) {
  out <- list(alpha = recovery(fit_obj, true_pars = obj$true_alpha, selection = "alpha"))
  if (!is.null(obj$true_mu)) {
    out$mu <- recovery(fit_obj, true_pars = obj$true_mu, selection = "mu")
  }
  if (!is.null(obj$true_beta)) {
    out$beta <- recovery(fit_obj, true_pars = obj$true_beta, selection = "beta")
  }
  out
}


print_hybrid_missing_summary <- function(obj) {
  cat("\nDropped low-cell subjects:\n")
  print(obj$dropped$low)
  cat("\nDropped between-cell subjects:\n")
  print(obj$dropped$between)

  cat("\nObserved cell presence by subject (1 = present):\n")
  print(obj$cell_summary$cell_presence)

  cat("\nNumber of observed Cond3 cells per subject:\n")
  print(obj$cell_summary$n_cells_per_subject)

  cat("\nFirst subject design summaries for parameter B:\n")
  print(obj$subject_designs[[1]])

  invisible(obj)
}


# Suggested manual runs:
#
# 1) Pure hybrid design: everyone has low + exactly one of med/high.
obj_hybrid <- build_hybrid_missing_test(
  n_subj = 24,
  n_trials = 160,
  n_drop_low = 0,
  n_drop_between = 0
)
print_hybrid_missing_summary(obj_hybrid)
fit_hybrid <- run_hybrid_missing_fit(obj_hybrid, iter = 1000)
#
# 2) Stress case: some subjects lose one more block and are left with one cell.
obj_stress <- build_hybrid_missing_test(
  n_subj = 24,
  n_trials = 160,
  n_drop_low = 8,
  n_drop_between = 7
)
print_hybrid_missing_summary(obj_stress)
fit_stress <- run_hybrid_missing_fit(obj_stress, iter = 1000)
pdf("test.pdf")
compare_to_truth(fit_stress, obj_stress)
dev.off()
#
# 3) True hybrid within/between example:
# everybody gets High, then each subject also gets either Low or Med.
# The shared High effect is handled in the subject-level design and the
# Low-vs-Med difference is handled through group_design() for both threshold
# and the non-High correct drift.
obj_group <- build_hybrid_group_design_test(
  n_subj = 24,
  n_trials = 160
)
print_hybrid_group_design_summary(obj_group)
fit_group <- run_hybrid_missing_fit(obj_group, iter = 1000)
pdf("test.pdf")
compare_to_truth(fit_group, obj_group)
dev.off()
