# Minimal fit templates for two hybrid race-model designs.
#
# This file keeps only the pieces needed to fit real data that already exist.
# It does not simulate data or run recovery checks.
#
# Usage:
#   source("WorkingTests/hybrid_missing_cells_fit_templates.R")
#   des1 <- build_missing_cells_rdm_design(my_data)
#   emc1 <- make_missing_cells_emc(my_data, des1)
#   fit1 <- fit_missing_cells_emc(emc1, iter = 1000)
#
#   source("WorkingTests/hybrid_missing_cells_fit_templates.R")
#   des2 <- build_group_design_rdm_subject_design()
#   gd2  <- build_group_design_rdm_group_design(subject_data, des2)
#   emc2 <- make_group_design_rdm_emc(trial_data, subject_data, des2, gd2)
#   fit2 <- fit_group_design_rdm_emc(emc2, iter = 1000)

suppressPackageStartupMessages({
  library(EMC2)
})


# -------------------------------------------------------------------------
# 1) Missing-cells version
# -------------------------------------------------------------------------
#
# Intended for the case where the design is conceptually within-subject, but
# some cells are missing after filtering/exclusion.
#
# Required trial-level columns in `data`:
#   - subjects
#   - S
#   - R
#   - rt
#   - Cond3   with levels low / med / high
#
# The design uses:
#   - drift:    v ~ 0 + Error + Correct:Cond3
#   - threshold B ~ Cond3

build_missing_cells_rdm_design <- function(data = NULL) {
  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)

  design(
    data = data,
    factors = if (is.null(data)) {
      list(subjects = 1, S = c("left", "right"), Cond3 = c("low", "med", "high"))
    } else {
      NULL
    },
    Rlevels = if (is.null(data)) c("left", "right") else NULL,
    functions = list(
      Correct = function(d) ifelse(d$lM == TRUE, 1, 0),
      Error = function(d) ifelse(d$lM == FALSE, 1, 0)
    ),
    matchfun = matchfun,
    model = RDM,
    formula = list(
      v ~ 0 + Error + (Correct:Cond3),
      B ~ Cond3,
      A ~ 1,
      t0 ~ 1,
      s ~ 1
    ),
    constants = c(s = log(1)),
    report_p_vector = FALSE
  )
}


make_missing_cells_emc <- function(data,
                                   design,
                                   compress = TRUE,
                                   n_chains = 3,
                                   prior_list = NULL) {
  if (is.null(prior_list)) {
    prior_list <- prior(design, type = "standard")
  }

  make_emc(
    data = data,
    design = design,
    type = "standard",
    prior_list = prior_list,
    compress = compress,
    n_chains = n_chains
  )
}


fit_missing_cells_emc <- function(emc,
                                  iter = 1000,
                                  particle_factor = 20,
                                  cores_per_chain = 1,
                                  cores_for_chains = 3,
                                  step_size = 100,
                                  verbose = TRUE) {
  fit(
    emc,
    iter = iter,
    particle_factor = particle_factor,
    cores_per_chain = cores_per_chain,
    cores_for_chains = cores_for_chains,
    step_size = step_size,
    verbose = verbose
  )
}


# -------------------------------------------------------------------------
# 2) True hybrid within/between version using group_design()
# -------------------------------------------------------------------------
#
# Intended for the case where everyone has High, and each subject also has
# exactly one of Low or Med. The shared High effect is modeled at the
# subject level; Low-vs-Med differences are modeled through group_design().
#
# Required trial-level columns in `trial_data`:
#   - subjects
#   - S
#   - R
#   - rt
#   - Load      with levels Low / Med / High
#
# Required subject-level columns in `subject_data`:
#   - subjects
#   - LoadGroup with levels Low / Med
#
# The design uses:
#   - drift:    v ~ 0 + Error + Correct:HighTag
#   - threshold B ~ HighTag
# and group-level effects on:
#   - B
#   - v_Correct:HighTagBase

build_group_design_rdm_subject_design <- function() {
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
      Correct = function(d) ifelse(d$lM == TRUE, 1, 0),
      Error = function(d) ifelse(d$lM == FALSE, 1, 0)
    ),
    model = RDM,
    formula = list(
      v ~ 0 + Error + (Correct:HighTag),
      B ~ HighTag,
      A ~ 1,
      t0 ~ 1,
      s ~ 1
    ),
    constants = c(s = log(1)),
    report_p_vector = FALSE
  )
}


build_group_design_rdm_group_design <- function(subject_data, subject_design) {
  group_design(
    formula = list(
      B ~ LoadGroup,
      `v_Correct:HighTagBase` ~ LoadGroup
    ),
    data = subject_data,
    subject_design = subject_design
  )
}


make_group_design_rdm_emc <- function(trial_data,
                                      subject_data,
                                      subject_design,
                                      group_design,
                                      compress = TRUE,
                                      n_chains = 3,
                                      prior_list = NULL) {
  if (is.null(prior_list)) {
    prior_list <- prior(
      subject_design,
      type = "standard",
      group_design = group_design
    )
  }

  make_emc(
    data = trial_data,
    design = subject_design,
    type = "standard",
    prior_list = prior_list,
    group_design = group_design,
    compress = compress,
    n_chains = n_chains
  )
}


fit_group_design_rdm_emc <- function(emc,
                                     iter = 1000,
                                     particle_factor = 20,
                                     cores_per_chain = 1,
                                     cores_for_chains = 3,
                                     step_size = 100,
                                     verbose = TRUE) {
  fit(
    emc,
    iter = iter,
    particle_factor = particle_factor,
    cores_per_chain = cores_per_chain,
    cores_for_chains = cores_for_chains,
    step_size = step_size,
    verbose = verbose
  )
}
