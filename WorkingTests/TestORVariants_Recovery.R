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
#   "OR_DETECTION_GNG"      – four-accumulator OR task with a withheld no.
#                             Accumulators: A, n_A, B, n_B (all required).
#                             S levels:     NN, AN, NB, AB.
#                             Responses:    "yes" / withheld (R = NA, rt = Inf).
#
# Rows with different LogicalRule values can coexist in one data frame / design.
# Extra accumulator roles specified in fixed_accumulator_roles are silently
# dormant for rules that don't use them.
#
# make_data() dispatches via LogicalRules_rfun when the data has a LogicalRule
# column and the model is LogicalRulesLBA.  All six rule codes are handled.
# =============================================================================

# --- Shared setup -------------------------------------------------------------

UC_val <- 0.8   # optional upper-censoring deadline for GNG omissions

# All four LogicalRules GNG roles.  Analytic detection uses only A and B;
# n_A and n_B are dormant there and active in the GNG/ordinary OR rules.
acc_roles <- factor(c("A","n_A","B","n_B"),
                    levels = c("A","n_A","B","n_B"))

# Rule-aware matchfun.
# For OR:        matched = accumulator activated by stimulus (standard LR logic).
# For analytic detection: matched = active target detectors A and B.
# For GNG/OR: matched = the target/nontarget pair for each stimulus condition.
joint_matchfun <- function(d) {
  is_analytic <- d$LogicalRule == "OR_DETECTION_ANALYTIC"
  has_a <- d$S %in% c("AN", "AB")
  has_b <- d$S %in% c("NB", "AB")
  ifelse(is_analytic,
    (d$lR == "A" & has_a) | (d$lR == "B" & has_b),
    (d$lR == "A" & has_a) | (d$lR == "n_A" & !has_a) |
      (d$lR == "B" & has_b) | (d$lR == "n_B" & !has_b)
  )
}

# Role indicators used in the formula.  The likelihood applies the logical
# rule; these functions only provide distinct parameter cells for the roles.
joint_funcs <- list(
  GoA  = function(d) ifelse(d$lR == "A", 1, 0),
  GoB  = function(d) ifelse(d$lR == "B", 1, 0),
  NegA = function(d) ifelse(d$lR == "n_A", 1, 0),
  NegB = function(d) ifelse(d$lR == "n_B", 1, 0)
)

# --- Trial template -----------------------------------------------------------
# OR:               all four S levels (NN elicits "no" under OR gate).
# OR_DETECTION_ANALYTIC: no NN (p_j = 0 → min_ll for NN × analytic detection).
# OR_DETECTION_GNG: all four S levels (NN permits false alarms; no is withheld).

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
                        levels = c("yes","no"))
)

joint_design <- design(
  data                    = trial_template,
  Rlevels                 = c("yes","no"),
  fixed_accumulator_roles = acc_roles,
  matchfun                = joint_matchfun,
  model                   = LogicalRulesLBA,
  constants               = c(sv = log(1)),
  functions               = joint_funcs,
  formula                 = list(v ~ 0 + GoA + GoB + NegA + NegB,
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
p_vec["v_GoA"]   <-  2.0          # drift for target A
p_vec["v_GoB"]   <-  1.8          # drift for target B
p_vec["v_NegA"]  <-  1.2          # drift for nontarget A
p_vec["v_NegB"]  <-  1.4          # drift for nontarget B
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
    p_yes   = round(mean(R == "yes", na.rm = TRUE), 2),
    p_no    = round(mean(R == "no",  na.rm = TRUE), 2),
    p_withheld = round(mean(is.na(R)), 2),
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
