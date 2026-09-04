RNGkind("L'Ecuyer-CMRG")
set.seed(123)


# Create subject-level design
subj_design <- design(data = forstmann, model = DDM,
                      formula = list(v ~ S, a ~ E, t0 ~ 1),
                      contrasts = list(S = contr.helmert))
# Add some age covariate and roughly demean
# Demeaning is important to ensure that the interpretation of the group-level intercept
# is the mean of the group (i.e., 'mu' still represents the group-level mean)
forstmann$age <- (as.numeric(forstmann$subjects) -mean(as.numeric(forstmann$subjects)))/
  max(as.numeric(forstmann$subjects))
# Create fake group column
forstmann$group <- ifelse(forstmann$subjects %in%
                            unique(forstmann$subjects)[seq(1, 19, 2)], "A", "B")

# Create group-level design matrices
group_des <- group_design(
  formula = list(v_S1 ~ age + group, a ~ age),
  data = forstmann,
  subject_design = subj_design,
  contrasts = list(group = contr.bayes)
)

pri <- prior(subj_design, group_design = group_des)
emc <- make_emc(forstmann, subj_design, group_design = group_des, prior_list = pri, compress = F, n_chains = 1)

test_that("groups", {
  expect_snapshot(print(group_des))
  expect_snapshot(summary(group_des))
  expect_snapshot(credint(pri, map = TRUE))
  expect_snapshot(credint(pri, selection = "beta"))
  expect_snapshot(init_chains(emc, particles = 10, cores_per_chain = 1)[[1]]$samples)
})

test_that("mapped_pars evaluates group-level predictors", {
  p <- sampled_pars(subj_design, group_design = group_des)
  p[] <- 0
  p[c("v", "v_S1", "v_S1_age", "v_S1_group1")] <- c(1, 2, 0, 1)

  mapped <- mapped_pars(pri, p)
  expect_true("group" %in% names(mapped))
  expect_setequal(unique(as.character(mapped$group)), c("A", "B"))

  # Priors carry the group design, so the default path should include its
  # predictors even when no coefficient vector is supplied explicitly.
  expect_true("group" %in% names(mapped_pars(pri)))

  # The group coefficient changes v_S1, which changes v differently for the
  # two groups even though `group` is absent from the subject-level formulae.
  left <- mapped[mapped$S == "left", , drop = FALSE]
  expect_equal(length(unique(left$v)), 2L)

  mapped_A <- mapped_pars(pri, p,
                          data = forstmann[forstmann$group == "A", , drop = FALSE])
  expect_setequal(unique(as.character(mapped_A$group)), "A")
})
