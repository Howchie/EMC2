test_that("LogicalRules detection treats missing stimulus as NN", {
  template <- data.frame(
    subjects = factor(c("s1", "s1")),
    S = factor(c("A", NA_character_), levels = c("A", "B", "AB")),
    LogicalRule = factor(rep("OR_DETECTION_ANALYTIC", 2),
                         levels = "OR_DETECTION_ANALYTIC"),
    R = factor(NA_character_, levels = c("yes", "no")),
    rt = NA_real_
  )

  design_template <- template
  design_template$R <- factor(NA_character_, levels = c("A", "B"))
  des <- design(
    data = design_template,
    Rlevels = c("yes", "no"),
    fixed_accumulator_roles = factor(c("A", "B"), levels = c("A", "B")),
    matchfun = function(d) {
      (d$lR == "A" & d$S %in% c("A", "AB")) |
        (d$lR == "B" & d$S %in% c("B", "AB"))
    },
    model = LogicalRulesLBA,
    formula = list(v ~ 0 + lR, B ~ 1, A ~ 1, t0 ~ 1, sv ~ 1)
  )

  p <- sampled_pars(des, doMap = FALSE)
  set.seed(1)
  simulated <- make_data(p, des, data = template)

  expect_true(is.na(simulated$R[2]))
  expect_true(is.infinite(simulated$rt[2]))
})
