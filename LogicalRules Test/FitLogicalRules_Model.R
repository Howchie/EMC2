## Script to simulate from, and test recovery of, RDMSWTN
rm(list=ls())
Packages <- c("EMC2","dplyr","sft")
lapply(Packages, library, character.only = TRUE) # to install change 'library' to 'install.packages'

## Read in Eidels 2015 data
load("LogicalRules Test/Eidels2015.RData")
data = all_data_acc
data = all_data_acc[all_data_acc$subjects=="AW", ]

source("LogicalRules Test/Functions.R")
source("LogicalRules Test/Contrasts.R")
designBushmakin <- design(
  fixed_accumulator_roles = factor(c("A","n_A","B","n_B"),levels=c("A","n_A","B","n_B")),
  matchfun=correctfun,
  model=function(){LogicalRulesLBA(posdrift=TRUE)},constants=c('sv_lMFALSE:LogicalRuleAND'=log(1),'sv_lMFALSE:LogicalRuleOR'=log(1),p_LogicalRuleAND=qnorm(1),q_LogicalRuleAND=qnorm(0.5)),
  formula=list(v~0+Stim,B~(0+Resp:LogicalRule),t0~0+LogicalRule,A~0+LogicalRule,sv~0+lM:LogicalRule,p~0+LogicalRule,q~0+LogicalRule),
  functions=list(Stim=RateYNMatchMismatch_Contrast,Resp=YN_Contrast,nDots=nTargets_Contrast),
  factors=list(subjects=factor(rep(1,8)),S=factor(rep(c("DT","UT","LT","NT"),2),levels=c("DT","UT","LT","NT")),LogicalRule=factor(rep(c("OR","AND"),each=4),levels=c("AND","OR"))),
  Rlevels=c("yes","no")
)

BushmakinModel = make_emc(data,designBushmakin,type="single",fileName = 'samples.RData')
BushmakinModel=fit(BushmakinModel,max_tries = 50, cores_per_chain=1, cores_for_chains=3, 
    stop_criteria = list(
    preburn = list(iter = 10), burn = list(mean_gd = 2.5), adapt = list(min_unique = 20),
    sample = list(iter = 25, max_gd = 2), verbose = FALSE, particle_factor = 30, step_size = 25))
save.image("BushmakinModel.RData")
# Pred Data
p_vec=EMC2::credint(BushmakinModel, probs = .5)[[1]]
pred=make_data(p_vec,designBushmakin,n_trials=50)
## Plotting
quantiles <- c(.1,.5,.9)#seq(from=0.05,to=0.95,by=0.05)
design_cells=data.frame(label=factor(levels(interaction(data$S,data$LogicalRule))))
pred$Correct=correctfun(pred);data$Correct=correctfun(data)
pred=pred[pred$rt<4,]
# Extract Quantiles by Design Cell x PM_trial x Correct
plotdata <- purrr::map_dfr(design_cells$label, function(lbl) {
  parts <- strsplit(as.character(lbl), "\\.")[[1]]
  names(parts) <- c("S", "LogicalRule")
  compute_quantile_points_by_subject(
    data_list = list("Recovered" = pred),
    factors = list(S = parts["S"], LogicalRule = parts["LogicalRule"]),
    design_label = lbl,
    correct_col = "Correct",
    quantiles=quantiles
  )
}) %>%
  group_by(Quantile, Source, Subset, DesignCell) %>%
  summarise(
    RT_mean = mean(RT, na.rm = TRUE),
    RT_se = sd(RT, na.rm = TRUE) / sqrt(n()),
    RT_lower = quantile(RT, 0.025, na.rm = TRUE),
    RT_upper = quantile(RT, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  separate(DesignCell, into = c("S", "LogicalRule"), sep = "\\.",remove = FALSE) %>%
  mutate(S=factor(S,levels=c("DT","UT","LT","NT")),LogicalRule=factor(LogicalRule,levels=c("AND","OR"))) %>%
  arrange(LogicalRule,S) %>%
  as.data.frame()

# Add any needed columns for grouping/factors
group_vars <- c("S", "LogicalRule", "Correct")

# Calculate quantiles per draw per group
pred_quantiles <- pred %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(
    quantile_vals = list(quantile(rt, probs = quantiles, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  unnest_wider(quantile_vals, names_sep = "_") %>%
  pivot_longer(
    starts_with("quantile_vals_"),
    names_to = "Quantile",
    names_prefix = "quantile_vals_",
    values_to = "RT"
  ) %>%
  mutate(Quantile = as.numeric(gsub("%", "", Quantile)) / 100) %>%
  group_by(S, LogicalRule, Correct, Quantile) %>%
  summarise(
    RT_mean = mean(RT, na.rm = TRUE),
    RT_lower = quantile(RT, 0.025, na.rm = TRUE),
    RT_upper = quantile(RT, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Source="Recovered") %>%
  filter(Correct==TRUE) %>%
  bind_rows(data %>%
              group_by(across(all_of(group_vars))) %>%
              summarise(
                quantile_vals = list(quantile(rt, probs = quantiles, na.rm = TRUE)),
                .groups = "drop"
              ) %>%
              unnest_wider(quantile_vals, names_sep = "_") %>%
              pivot_longer(
                starts_with("quantile_vals_"),
                names_to = "Quantile",
                names_prefix = "quantile_vals_",
                values_to = "RT"
              ) %>%
              mutate(Quantile = as.numeric(gsub("%", "", Quantile)) / 100) %>%
              group_by(S, LogicalRule, Correct, Quantile) %>%
              summarise(
                RT_mean = mean(RT, na.rm = TRUE),
                .groups = "drop"
              ) %>%
              mutate(Source="Data") %>%
              filter(Correct==TRUE))

# Extract Quantiles by Design Cell x PM_trial x Correct
accuracy_data <- purrr::map_dfr(design_cells$label, function(lbl) {
  parts <- strsplit(as.character(lbl), "\\.")[[1]]
  names(parts) <- c("S", "LogicalRule")
  
  compute_accuracy_by_subject(
    data_list = list("Recovered" = pred),
    factors = list(S = parts["S"], LogicalRule = parts["LogicalRule"]),
    design_label = lbl,
    correct_col = "Correct"
  )
}) %>%
  group_by(Source, Subset, DesignCell) %>%
  summarise(
    Accuracy = mean(Accuracy),
    .groups = "drop"
  ) %>%
  separate(DesignCell, into = c("S", "LogicalRule"), sep = "\\.",remove = FALSE) %>%
  mutate(S=factor(S,levels=c("DT","UT","LT","NT")),LogicalRule=factor(LogicalRule,levels=c("AND","OR"))) %>%
  arrange(LogicalRule,S) %>%
  as.data.frame()

# Filtered subsets for clarity
plotdata_correct <- plotdata[plotdata$Subset == "correct", ]
# Lines: only for "Recovered"
lines_data <- subset(plotdata_correct, Source == "Recovered")
# Points: only for "Observed"
points_data <- subset(plotdata_correct, Source == "Data")

ggplot() +
  # Lines for "Recovered"
  geom_line(data = lines_data, 
            aes(x = S, y = RT_mean, group = interaction(Quantile,LogicalRule)), 
            size = 1) +
  # Points for "Observed"
  geom_point(data = plotdata_correct, 
             aes(x = S, y = RT_mean, shape = interaction(LogicalRule,Source)), 
             size = 3, fill = "white") +
  theme_bw(base_size = 14) +
  labs(
    title = "RTs",
    y = "RT (s)",
    x = "Quantile"
  ) +
  facet_wrap(~LogicalRule) +
  scale_shape_manual(values = c(21, 24, 22, 23)) +
  scale_color_brewer(palette = "Set1") +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    panel.grid.minor = element_blank()
  ) + ylim(0.3,2)

## Accuracy Data
# Lines: only for "Recovered"
lines_data <- subset(accuracy_data, Source == "Recovered")
# Points: only for "Observed"
points_data <- subset(accuracy_data, Source == "Data")

ggplot() +
  # Lines for "Recovered"
  geom_line(data = lines_data, 
            aes(x = S, y = Accuracy, linetype = LogicalRule, group = interaction(Source,LogicalRule)), 
            size = 1) +
  # Points for "Observed"
  geom_point(data = points_data, 
             aes(x = S, y = Accuracy, shape = LogicalRule), 
             size = 2, fill = "white") +
  theme_bw(base_size = 14) +
  labs(
    title = "Accuracy",
    y = "Accuracy",
    x = "S"
  ) +
  scale_shape_manual(values = c(21, 24, 22, 23)) +
  scale_color_brewer(palette = "Set1") +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    panel.grid.minor = element_blank()
  ) + ylim(.5,1.00)

# 
# tmp1 <- msm::rtnorm(1000, mean = 0.904 , sd=0.259,lower=0)
# tmp2 <- msm::rtnorm(1000, mean = -4.201, sd=1,lower=0)
# 
# # Plot density of first
# plot(density(tmp1), col = "blue", lwd = 2, 
#      main = "Density comparison", xlab = "Value", ylim = c(0, max(density(tmp1)$y, density(tmp2)$y)))
# # Overlay density of second
# lines(density(tmp2), col = "red", lwd = 2)
# 
# legend("topright", legend = c("tmp1", "tmp2"), col = c("blue", "red"), lwd = 2)
# 
truncnorm_mean <- function(mu, sigma, a, b) {
  alpha <- (a - mu) / sigma
  beta <- (b - mu) / sigma
  Z <- pnorm(beta) - pnorm(alpha)
  mu + (dnorm(alpha) - dnorm(beta)) / Z * sigma
}
# Example: mean 0, sd 1, truncation [-1, 2]
truncnorm_mean(, 0.320011, -1, 2)