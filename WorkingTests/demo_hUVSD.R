# Demo script for Hierarchical Unequal-Variance Signal Detection (hUVSD)
# Verifies R vs C++ likelihoods and demonstrates recovery fit.

rm(list=ls())
library(EMC2)

# 1. Define Design
# Following Lages (2024) hUVSD model parameterization:
# d  : Sensitivity (distance between means)
# c  : Bias (threshold deviation from midpoint)
# sd : Signal SD (Noise SD is fixed at 1.0)
design_huvsd <- design(
  Rlevels = c("No", "Yes"),
  factors = list(subjects = 1:3, S = c("Noise", "Signal")),
  formula = list(d ~ 1, c ~ 1, sd ~ 1),
  model = hUVSD
)

# 2. Simulate Data
# True parameters (on the scale they are sampled: sd is log-scale)
# d = 1.5, c = 0.2, sd = log(1.5) -> Signal SD = 1.5
p_vector <- c(d = 1.5, c = 0.2, sd = log(1.5))
data <- make_data(p_vector, design = design_huvsd, n_trials = 500)

# 3. Setup DADM (Data Augmented Design Matrix)
dadm <- make_dadm(data, design_huvsd)

# 4. Compare R vs C++ Likelihoods
cat("--- Likelihood Verification ---\n")
model_huvsd <- hUVSD()

# Use subject 1 for comparison
dadm1 <- dadm[[1]]
p_mat <- matrix(p_vector, nrow = 1)
colnames(p_mat) <- names(p_vector)

# R Likelihood
pars_mapped <- mapped_pars(dadm1, p_vector)
ll_r <- model_huvsd$log_likelihood(pars_mapped, dadm1, model_huvsd)

# C++ Likelihood (via calc_ll_oo)
# Note: calc_ll_oo is an internal function used by samplers
ll_c <- EMC2:::calc_ll_oo(
  particle_matrix = p_mat,
  data = dadm1,
  constants = attr(dadm1, "constants"),
  designs = attr(dadm1, "design_list"),
  type = model_huvsd$c_name,
  bounds = model_huvsd$bound,
  transforms = model_huvsd$transform,
  pretransforms = list(),
  p_types = names(model_huvsd$p_types),
  min_ll = log(1e-10)
)

cat("R Log-Likelihood Total: ", sum(ll_r), "\n")
cat("C++ Log-Likelihood Total:", ll_c, "\n")
cat("Difference:              ", sum(ll_r) - ll_c, "\n\n")

# 5. Recovery Fit Setup
cat("--- Recovery Fit Setup ---\n")
# Define priors
prior_huvsd <- prior(design_huvsd, 
                     mu_mean = c(d=1, c=0, sd=0), 
                     mu_sd = c(d=2, c=1, sd=1))

# Initialize samplers
samplers <- make_emc(data, design_huvsd, prior_huvsd)

# The following lines would run the fit:
# samplers <- run_emc(samplers, stage = "init", iter = 100)
# samplers <- run_emc(samplers, stage = "sample", iter = 500)
# summary(samplers)

cat("hUVSD demo completed successfully.\n")
