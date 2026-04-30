rm(list = ls())
library(EMC2)

RNGkind("L'Ecuyer-CMRG")
set.seed(20260424)

# Binary match function used by both tasks
matchfun <- function(d) d$S == d$lR

extract_mean_vector <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.list(x)) {
    vecs <- lapply(x, extract_mean_vector)
    vecs <- Filter(Negate(is.null), vecs)
    if (length(vecs) == 0) return(NULL)
    all_names <- unique(unlist(lapply(vecs, names)))
    mat <- matrix(NA_real_, nrow = length(all_names), ncol = length(vecs),
                  dimnames = list(all_names, NULL))
    for (i in seq_along(vecs)) mat[names(vecs[[i]]), i] <- vecs[[i]]
    return(rowMeans(mat, na.rm = TRUE))
  }
  if (is.matrix(x)) {
    v <- rowMeans(x)
    names(v) <- rownames(x)
    return(v)
  }
  if (is.numeric(x) && !is.null(names(x))) return(x)
  NULL
}

extract_mean_matrix <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.list(x)) {
    mats <- lapply(x, extract_mean_matrix)
    mats <- Filter(Negate(is.null), mats)
    if (length(mats) == 0) return(NULL)
    arr <- simplify2array(mats)
    return(apply(arr, c(1, 2), mean, na.rm = TRUE))
  }
  if (is.array(x) && length(dim(x)) == 3) {
    return(apply(x, c(1, 2), mean, na.rm = TRUE))
  }
  if (is.matrix(x)) return(x)
  NULL
}

n_subj <- 8

# Use average/difference coding for lM (same style as package tests)
ADmat <- matrix(c(-1 / 2, 1 / 2), ncol = 1, dimnames = list(NULL, "d"))

# -------------------------------
# 1) Build two task designs
# -------------------------------
design_lba <- design(
  factors = list(subjects = as.character(seq_len(n_subj)), S = c("left", "right")),
  Rlevels = c("left", "right"),
  matchfun = matchfun,
  model = LBA,
  formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1),
  contrasts = list(v = list(lM = ADmat)),
  constants = c(sv = log(1))
)

design_rdm <- design(
  factors = list(subjects = as.character(seq_len(n_subj)), S = c("left", "right")),
  Rlevels = c("left", "right"),
  matchfun = matchfun,
  model = RDM,
  formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1),
  contrasts = list(v = list(lM = ADmat)),
  constants = c(s = log(1))
)

# ----------------------------------------
# 2) Specify known group means (true mu)
# ----------------------------------------
mu_lba <- sampled_pars(design_lba, doMap = FALSE)
mu_lba[] <- c(
  v = 1.20,
  v_lMd = 0.90,
  B = log(1.8),
  A = log(0.5),
  t0 = log(0.25)
)

mu_rdm <- sampled_pars(design_rdm, doMap = FALSE)
mu_rdm[] <- c(
  v = log(1.10),
  v_lMd = 0.55,
  B = log(1.7),
  A = log(0.45),
  t0 = log(0.22)
)

# ------------------------------------------------------------
# 3) Simulate subject-level effects with joint cross-task Sigma
# ------------------------------------------------------------
joint_par_names <- names(sampled_pars(list(lba = design_lba, rdm = design_rdm), doMap = FALSE))
mu_joint <- c(
  setNames(mu_lba, paste0("lba|", names(mu_lba))),
  setNames(mu_rdm, paste0("rdm|", names(mu_rdm)))
)
mu_joint <- mu_joint[joint_par_names]

p <- length(mu_joint)
Sigma_joint <- diag(0.04, p)
colnames(Sigma_joint) <- rownames(Sigma_joint) <- joint_par_names

# Add stronger within-task and cross-task covariance to make recovery visible.
Sigma_joint["lba|t0", "lba|B"] <- 0.018
Sigma_joint["lba|B", "lba|t0"] <- 0.018
Sigma_joint["rdm|t0", "rdm|B"] <- 0.016
Sigma_joint["rdm|B", "rdm|t0"] <- 0.016
Sigma_joint["lba|t0", "rdm|t0"] <- 0.022
Sigma_joint["rdm|t0", "lba|t0"] <- 0.022
Sigma_joint["lba|B", "rdm|B"] <- 0.020
Sigma_joint["rdm|B", "lba|B"] <- 0.020
Sigma_joint["lba|v", "rdm|v"] <- 0.015
Sigma_joint["rdm|v", "lba|v"] <- 0.015

subj_joint <- mvtnorm::rmvnorm(n_subj, mean = mu_joint, sigma = Sigma_joint)
subj_joint <- as.matrix(subj_joint)
colnames(subj_joint) <- joint_par_names
rownames(subj_joint) <- as.character(seq_len(n_subj))

subj_lba <- subj_joint[, grepl("^lba\\|", colnames(subj_joint)), drop = FALSE]
subj_rdm <- subj_joint[, grepl("^rdm\\|", colnames(subj_joint)), drop = FALSE]
colnames(subj_lba) <- sub("^lba\\|", "", colnames(subj_lba))
colnames(subj_rdm) <- sub("^rdm\\|", "", colnames(subj_rdm))

# -----------------------------------
# 4) Simulate data for each component
# -----------------------------------
dat_lba <- make_data(subj_lba, design_lba, n_trials = 120)
dat_rdm <- make_data(subj_rdm, design_rdm, n_trials = 120)

cat("Rows per task:\n")
print(c(lba = nrow(dat_lba), rdm = nrow(dat_rdm)))

# -----------------------------------------
# 5) Build and fit a joint hierarchical EMC
# -----------------------------------------
joint_design <- list(lba = design_lba, rdm = design_rdm)
joint_data <- list(lba = dat_lba, rdm = dat_rdm)

prior_joint <- prior(joint_design, type = "standard")

emc_joint <- make_emc(
  data = joint_data,
  design = joint_design,
  type = "standard",
  prior_list = prior_joint,
  n_chains = 2,
  compress = TRUE
)

# Keep this deliberately short for a working example script.
# Increase iter/max_sample_iter for real recovery runs.
fit_joint <- fit(
  emc_joint,
  stop_criteria = list(
    sample = list(
      iter = 150,
      max_gd = 1.20,
      max_flat_loc = 0.80,
      flat_selection = c("alpha", "subj_ll"),
      flat_p1 = 1 / 3,
      flat_p2 = 1 / 3,
      max_sample_iter = 400
    ),
    cores_per_chain = 1,
    cores_for_chains = 2
  ),
  max_tries = 10
)

# ---------------------------------
# 6) Quick summary against truth
# ---------------------------------
mu_post <- extract_mean_vector(get_pars(
  fit_joint,
  selection = "mu",
  map = FALSE,
  return_mcmc = FALSE,
  merge_chains = TRUE
))
cor_post <- extract_mean_matrix(get_pars(
  fit_joint,
  selection = "correlation",
  map = FALSE,
  return_mcmc = FALSE,
  merge_chains = TRUE
))
true_cor <- cov2cor(Sigma_joint)

cat("\nPosterior mu (posterior mean):\n")
if (is.null(mu_post)) {
  cat("Could not extract posterior mean vector from fit object.\n")
} else {
  print(round(mu_post, 3))
}

cat("\nTrue generating values (joint names):\n")
truth_joint <- mu_joint
print(round(truth_joint, 3))

cov_pairs <- rbind(
  c("lba|t0", "rdm|t0"),
  c("lba|B", "rdm|B"),
  c("lba|v", "rdm|v"),
  c("lba|t0", "lba|B"),
  c("rdm|t0", "rdm|B")
)
cat("\nTrue vs recovered correlations (posterior mean):\n")
for (i in seq_len(nrow(cov_pairs))) {
  p1 <- cov_pairs[i, 1]
  p2 <- cov_pairs[i, 2]
  rec <- if (!is.null(cor_post)) cor_post[p1, p2] else NA_real_
  cat(sprintf("%s <-> %s: true=%.3f, recovered=%.3f\n", p1, p2, true_cor[p1, p2], rec))
}

saveRDS(
  list(
    true = list(mu = mu_joint, Sigma = Sigma_joint, cor = true_cor),
    subj_pars = list(lba = subj_lba, rdm = subj_rdm),
    data = joint_data,
    fit = fit_joint
  ),
  file = "WorkingTests/joint_lba_rdm_sim_fit.rds"
)

cat("\nSaved results to WorkingTests/joint_lba_rdm_sim_fit.rds\n")
