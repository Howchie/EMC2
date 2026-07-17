source("R/EMC2-package.R")
devtools::load_all()

design <- design(
  factors = list(subjects=1, LogicalRule=c("OR", "AND", "XOR", "ID")),
  Rlevels = c("A", "n_A", "B", "n_B", "yes", "no"),
  model = LogicalRulesLBA(),
  formula = list(v~LogicalRule, B~LogicalRule, A~1, t0~1, pContaminant~1)
)
pars <- c(v = 1, v_LogicalRuleAND = 1, v_LogicalRuleXOR = 1, v_LogicalRuleID = 1,
          B = 1, B_LogicalRuleAND = 1, B_LogicalRuleXOR = 1, B_LogicalRuleID = 1,
          A = 1, t0 = 1, pContaminant = qnorm(0.5))
data <- make_data(pars, design, n_trials=100)
print(sum(is.na(data$R))/nrow(data))
