rm(list = ls())

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Need pkgload to load local EMC2 sources.")
}
pkgload::load_all(".")

RNGkind("L'Ecuyer-CMRG")
set.seed(20260502)

matchfun <- function(d) d$S == d$lR
ADmat <- matrix(c(-1 / 2, 1 / 2), ncol = 1, dimnames = list(NULL, "d"))

template <- data.frame(
  subjects = factor(c("s1", "s1")),
  S = factor(c("left", "right"), levels = c("left", "right")),
  R = factor(c("left", "right"), levels = c("left", "right"))
)

make_design_rgamma <- function(shape) {
  design(
    data = template,
    matchfun = matchfun,
    model = RGAMMA,
    formula = list(lambda ~ lM),
    contrasts = list(lambda = list(lM = ADmat)),
    constants = c(shape = log(shape), shift = log(0))
  )
}

fit_one <- function(dat, design_obj, label) {
  cat("\n--- Fitting:", label, "---\n")
  emc <- make_emc(dat, design_obj, type = "single", n_chains = 3, compress = FALSE)
  fit(
    emc,
    stop_criteria = list(
      sample = list(
        iter = 140,
        max_gd = 1.2,
        max_flat_loc = 0.9,
        flat_selection = c("alpha", "subj_ll"),
        flat_p1 = 1 / 3,
        flat_p2 = 1 / 3,
        max_sample_iter = 500
      ),
      cores_per_chain = 1,
      cores_for_chains = 3
    ),
    max_tries = 6
  )
}

summ_recovery <- function(fit_obj, p_true, label) {
  cat("\nRecovery:", label, "\n")
  rec <- recovery(fit_obj, true_pars = p_true)
  print(rec)
  rec
}

design_s1 <- make_design_rgamma(1L)
design_s2 <- make_design_rgamma(2L)

# Keep generating parameters on the same sampling scale for both datasets.
p_true <- sampled_pars(design_s1, doMap = FALSE)
p_true["lambda"] <- log(2.0)
p_true["lambda_lMd"] <- 0.6

cat("True parameter vector (sampling scale):\n")
print(round(p_true, 4))

cat("\nTrue mapped parameters in shape-1 design:\n")
mp <- mapped_pars(design_s1, p_true)
num_cols <- vapply(mp, is.numeric, logical(1))
mp[num_cols] <- lapply(mp[num_cols], function(x) round(x, 4))
print(mp)

set.seed(20260502)
dat_s1 <- make_data(p_true, design_s1, n_trials = 600)
set.seed(20260503)
dat_s2 <- make_data(p_true, design_s2, n_trials = 600)

cat("\nData summaries:\n")
cat("shape-1 data rt mean =", round(mean(dat_s1$rt), 4), "; n =", nrow(dat_s1), "\n")
cat("shape-2 data rt mean =", round(mean(dat_s2$rt), 4), "; n =", nrow(dat_s2), "\n")

# Fit both model types to both datasets (matched + mismatched)
fit_s1_on_s1 <- fit_one(dat_s1, design_s1, "model shape=1 on data shape=1")
fit_s2_on_s1 <- fit_one(dat_s1, design_s2, "model shape=2 on data shape=1")
fit_s1_on_s2 <- fit_one(dat_s2, design_s1, "model shape=1 on data shape=2")
fit_s2_on_s2 <- fit_one(dat_s2, design_s2, "model shape=2 on data shape=2")

rec_s1_on_s1 <- summ_recovery(fit_s1_on_s1, p_true, "model shape=1 on data shape=1")
rec_s2_on_s2 <- summ_recovery(fit_s2_on_s2, p_true, "model shape=2 on data shape=2")

cat("\n--- Model comparison on data generated from shape=1 ---\n")
cmp_on_s1 <- compare(
  list(model_shape1 = fit_s1_on_s1, model_shape2 = fit_s2_on_s1),
  cores_for_props = 1
)
print(cmp_on_s1)

cat("\n--- Model comparison on data generated from shape=2 ---\n")
cmp_on_s2 <- compare(
  list(model_shape1 = fit_s1_on_s2, model_shape2 = fit_s2_on_s2),
  cores_for_props = 1
)
print(cmp_on_s2)

saveRDS(
  list(
    true_pars = p_true,
    data = list(shape1 = dat_s1, shape2 = dat_s2),
    fits = list(
      s1_on_s1 = fit_s1_on_s1,
      s2_on_s1 = fit_s2_on_s1,
      s1_on_s2 = fit_s1_on_s2,
      s2_on_s2 = fit_s2_on_s2
    ),
    recovery = list(s1_on_s1 = rec_s1_on_s1, s2_on_s2 = rec_s2_on_s2),
    comparison = list(on_shape1 = cmp_on_s1, on_shape2 = cmp_on_s2)
  ),
  file = "WorkingTests/sim_recover_rgamma.rds"
)

cat("\nSaved to WorkingTests/sim_recover_rgamma.rds\n")
