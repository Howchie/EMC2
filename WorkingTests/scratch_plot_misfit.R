rm(list = ls())
library(EMC2)
library(ggplot2)
library(dplyr)
library(tidyr)

matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)
designROU <- design(
  factors = list(subjects = 1, S = c("left", "right")),
  Rlevels = c("left", "right"),
  matchfun = matchfun,
  model = ROU(parameterization="equilibrium"),
  formula = list(theta ~ lM, B ~ 1, A ~ 1, t0 ~ 1, chi ~ 1, tk ~ 1),
  constants = c(B = log(1),A=log(0))
)
p_vector <- sampled_pars(designROU, doMap = FALSE)
p_vector["t0"] <- log(0.15)
p_vector["chi"] <- log(1)
p_vector["theta"] <- log(1)
p_vector["theta_lMTRUE"] <- log(2.5)
p_vector["tk"] <- log(.5)

set.seed(42)
dat <- make_data(p_vector, designROU, n_trials = 10000)
dat$Correct <- as.numeric(dat$S==dat$R)

designRDM <- design(
  factors = list(subjects = 1, S = c("left", "right")),
  Rlevels = c("left", "right"),
  matchfun = matchfun,
  model = RDM(),
  formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1),
  constants = c(s = log(1),A=log(0))
)

emc2 <- make_emc(dat, designRDM, type = "single", rt_resolution = 1/60)

cat("Fitting RDM...\n")
fit_res2 <- fit(
  emc2,
  cores_for_chains = 3,
  cores_per_chain = 4,
  stop_criteria = list(
    sample = list(
      iter = 200,
      max_gd = 1.10,
      max_flat_loc = 2,
      flat_selection = c("alpha", "subj_ll"),
      flat_p1 = 1/3,
      flat_p2 = 1/3,
      max_sample_iter = 500
    )
  ),
  max_tries = 1
)

cat("Predicting...\n")
pred <- predict(fit_res2)
pred$Correct <- as.numeric(pred$S == pred$R)

dat$Type <- "True Data (ROU)"
pred$Type <- "Best RDM Fit"

combined <- bind_rows(
  dat %>% select(rt, Correct, Type),
  pred %>% select(rt, Correct, Type)
) %>%
  mutate(Response = ifelse(Correct == 1, "Correct", "Error"))

# Plot density
p1 <- ggplot(combined, aes(x = rt, color = Type, linetype = Type)) +
  geom_density(linewidth = 1) +
  facet_wrap(~Response, scales = "free_y") +
  theme_minimal() +
  labs(title = "RT Densities: Data vs RDM Fit",
       x = "Response Time (s)", y = "Density") +
  coord_cartesian(xlim=c(0, 1.5))

# Plot CDF
# Calculate empirical CDFs
cdf_data <- combined %>%
  group_by(Type, Response) %>%
  arrange(rt) %>%
  mutate(cdf = row_number() / n())

p2 <- ggplot(cdf_data, aes(x = rt, y = cdf, color = Type, linetype = Type)) +
  geom_line(linewidth = 1) +
  facet_wrap(~Response) +
  theme_minimal() +
  labs(title = "RT CDFs: Data vs RDM Fit",
       x = "Response Time (s)", y = "Cumulative Probability") +
  coord_cartesian(xlim=c(0, 1.5))

library(gridExtra)
p <- arrangeGrob(p1, p2, ncol = 1)

ggsave("/home/ubuntu/.gemini/antigravity-cli/brain/d1585385-429c-4651-9eea-3db9cab0703a/scratch/rdm_misfit.png", p, width = 8, height = 8)
cat("Plot saved.\n")
