rm(list = ls())
library(EMC2)
library(dplyr)

# =============================================================================
# HOW TRIALS ARE SPECIFIED FOR LogicalRulesLBA
#
# Every row of the template / data frame must carry a `LogicalRule` column.
# The six supported values and what they require:
#
#   "OR"                    – 4-accumulator competing-channels gate.
#                             Accumulators: A, n_A, B, n_B (required).
#                             S levels:     NN, AN, NB, AB.
#                             Responses:    "yes" (first YES fires) / "no"
#                                           (last NO fires).
#
#   "AND"                   – same 4 accumulators; YES needs both channels YES.
#
#   "XOR"                   – same; YES iff exactly one channel YES.
#
#   "ID"                    – same; response identifies winning channel directly
#                             (R ∈ {NN, AN, NB, AB}).
#
#   "OR_DETECTION_ANALYTIC" – 2-accumulator detection (no nogo).
#                             Accumulators: A, B (required).
#                             S levels:     AN, NB, AB  (NN → p_j = 0 → min_ll).
#                             Responses:    "yes" (first detector fires) / omit.
#
#   "OR_DETECTION_GNG"      – 3-accumulator detection with nogo deadline.
#                             Accumulators: A, B, nogo (all required).
#                             S levels:     NN, AN, NB, AB.
#                             Responses:    "yes" / "nogo" (nogo wins before UC).
#                             Requires:     UC (upper criterion) in design or TC.
#
# Rows with different LogicalRule values can coexist in one data frame / design.
# Extra accumulator roles specified in fixed_accumulator_roles are silently
# dormant for rules that don't use them (e.g. nogo is ignored for OR).
#
# make_data() dispatches via LogicalRules_rfun when the data has a LogicalRule
# column and the model is LogicalRulesLBA.  All six rule codes are handled.
# =============================================================================

# --- Shared setup -------------------------------------------------------------

UC_val <- 0.8   # response deadline for OR_DETECTION_GNG

# All five accumulator roles.  Extra roles are dormant for rules that don't
# need them; no special action required.
acc_roles <- factor(c("A","n_A","B","n_B","nogo"),
                    levels = c("A","n_A","B","n_B","nogo"))

# Rule-aware matchfun.
# For OR:        matched = accumulator activated by stimulus (standard LR logic).
# For detection: matched = go detectors A, B, and nogo (they compete directly).
# n_A, n_B are dormant on detection trials; their lM value doesn't affect the
# likelihood but must produce valid (finite) drift rates.
joint_matchfun <- function(d) {
  dplyr::case_when(
    d$LogicalRule == "OR" & d$S == "NN" & d$lR %in% c("n_A","n_B") ~ TRUE,
    d$LogicalRule == "OR" & d$S == "AN" & d$lR %in% c("A","n_B")   ~ TRUE,
    d$LogicalRule == "OR" & d$S == "NB" & d$lR %in% c("n_A","B")   ~ TRUE,
    d$LogicalRule == "OR" & d$S == "AB" & d$lR %in% c("A","B")     ~ TRUE,
    d$LogicalRule %in% c("OR_DETECTION_ANALYTIC","OR_DETECTION_GNG") &
      d$lR %in% c("A","B","nogo")                                   ~ TRUE,
    TRUE ~ FALSE
  )
}

# Contrast functions used in the formula.
# `match`  = 1 for non-nogo accumulators that are activated by the stimulus.
# `nogo`   = 1 for the nogo accumulator (exclusive; no overlap with match).
# Mismatched non-nogo accumulators get 0 on both → v = 0 (effectively dormant).
joint_funcs <- list(
  match = function(d) ifelse(as.logical(d$lM) & d$lR != "nogo", 1, 0),
  nogo  = function(d) ifelse(d$lR == "nogo", 1, 0)
)

# --- Trial template -----------------------------------------------------------
# OR:               all four S levels (NN elicits "no" under OR gate).
# OR_DETECTION_ANALYTIC: no NN (p_j = 0 → min_ll for NN × analytic detection).
# OR_DETECTION_GNG: all four S levels (NN → nogo wins → nogo response).

trial_template <- data.frame(
  subjects    = factor("s1"),
  S           = factor(c("NN","AN","NB","AB",
                          "AN","NB","AB",
                          "NN","AN","NB","AB"),
                        levels = c("NN","AN","NB","AB")),
  LogicalRule = factor(c(rep("OR",                    4),
                          rep("OR_DETECTION_ANALYTIC", 3),
                          rep("OR_DETECTION_GNG",      4)),
                        levels = c("OR","OR_DETECTION_ANALYTIC","OR_DETECTION_GNG")),
  R           = factor(rep(NA_character_, 11),
                        levels = c("yes","no","nogo"))
)

joint_design <- design(
  data                    = trial_template,
  Rlevels                 = c("yes","no","nogo"),
  fixed_accumulator_roles = acc_roles,
  matchfun                = joint_matchfun,
  model                   = LogicalRulesLBA,
  constants               = c(sv = log(1)),
  functions               = joint_funcs,
  formula                 = list(v ~ 0 + match + nogo,
                                  B  ~ 1,
                                  t0 ~ 1,
                                  A  ~ 1,
                                  sv ~ 1),
  UC = UC_val
)

# Inspect parameter names
p_vec <- sampled_pars(joint_design, doMap = FALSE)
cat("Parameter names:", paste(names(p_vec), collapse = ", "), "\n")

# True parameters (natural scale shown in comments)
p_vec["v_match"] <-  2.0          # drift for activated accumulator
p_vec["v_nogo"]  <-  1.5          # drift for nogo accumulator
p_vec["B"]       <-  log(0.8)     # threshold ~ 0.8
p_vec["t0"]      <-  log(0.2)     # non-decision time ~ 0.2 s
p_vec["A"]       <-  log(0.4)     # start-point noise ~ 0.4

cat("\nMapped (natural-scale) parameters:\n")
print(round(mapped_pars(joint_design, p_vec), 3))

# --- Simulate data ------------------------------------------------------------
set.seed(123)
dat <- make_data(p_vec, joint_design,
                 data   = trial_template,
                 expand = 400,
                 TC     = list(UC = UC_val))

cat("\nSimulated data summary (n rows =", nrow(dat), "):\n")
dat %>%
  group_by(LogicalRule, S) %>%
  summarise(
    n       = n(),
    p_yes   = round(mean(R == "yes",  na.rm = TRUE), 2),
    p_no    = round(mean(R == "no",   na.rm = TRUE), 2),
    p_nogo  = round(mean(R == "nogo", na.rm = TRUE), 2),
    RT_mean = round(mean(rt[is.finite(rt)], na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  print(n = Inf)

# --- Likelihood profile -------------------------------------------------------
# Each parameter is swept ± range around its true value; peak should sit at
# the true value.
source("WorkingTests/test_likelihood_plotfuns_ah.R")
timing <- system.time(
  profile_plot_test(
    dat, joint_design, p_vec,
    n_cores      = 3,
    range        = 1,
    layout       = c(2, 3),
    use_c        = TRUE,
    figure_title = "OR variants joint design",
    natural      = TRUE
  )
)
cat(sprintf("\nLikelihood profile elapsed: %.1f s\n", timing["elapsed"]))

# --- MCMC recovery ------------------------------------------------------------
emc <- make_emc(dat, joint_design, type = "single")

emc <- fit(
  emc,
  stop_criteria = list(
    sample = list(
      iter            = 1000,
      max_gd          = 1.10,
      max_flat_loc    = 0.5,
      flat_selection  = c("alpha","subj_ll"),
      flat_p1         = 1/3,
      flat_p2         = 1/3,
      max_sample_iter = 5000
    ),
    cores_per_chain = 3,
    cores_for_chains = 3
  ),
  max_tries = 30
)

recovery(emc, true_pars = p_vec)
post_predict <- predict(emc, n_post = 50)
plot_pars(emc, post_predict = post_predict, true_pars = p_vec)
save.image("TestORVariants_Recovery.RData")
