## Script to simulate from, and test recovery of, RDMSWTN
rm(list=ls())
Packages <- c("EMC2","dplyr","ggplot2","ggforce")
lapply(Packages, library, character.only = TRUE) # to install change 'library' to 'install.packages'

## Read in Eidels 2015 data
load("LogicalRules Test/Eidels2015.RData"); data=all_data_acc
data = all_data_acc[all_data_acc$subjects=="AW", ]

source("LogicalRules Test/Functions.R")
source("LogicalRules Test/Contrasts.R")
designBushmakin <- design(data=data,
  fixed_accumulator_roles = factor(c("A","n_A","B","n_B"),levels=c("A","n_A","B","n_B")),
  matchfun=matchfun,
  model=function(){LogicalRulesLBA(posdrift=TRUE)},constants=c('sv_lMFALSE'=log(1)),
  formula=list(v~0+Stim,B~(0+Resp),t0~0+nDots,A~1,sv~0+lM),
  functions=list(Stim=RateYNMatchMismatch_ContrastShared,Resp=YN_Contrast,nDots=nTargets_Contrast),
  Rlevels=factor(c("yes","no"),levels=c("yes","no")),UC=4.01
)

p_vector <- sampled_pars(designBushmakin,doMap=FALSE)
p_vector[grepl("^B", names(p_vector))] <- log(1)
p_vector[grepl("^t0", names(p_vector))] <- log(0.3)
p_vector[grepl("^A", names(p_vector))] <- log(1)
p_vector[grepl("^sv", names(p_vector))] <- log(1)
p_vector[regexpr("^v",names(p_vector))==1] <- 1
p_vector[regexpr("^v.*Match",names(p_vector))==1] <- 2
p_vector[regexpr("^v.*Mismatch",names(p_vector))==1] <- 1

# variances
s_vector <- sampled_pars(designBushmakin,doMap=FALSE)
s_vector[grepl("^B", names(s_vector))] <- 1
s_vector[grepl("^t0", names(s_vector))] <- .2
s_vector[grepl("^A", names(s_vector))] <- 1
s_vector[regexpr("^sv",names(s_vector))==1] <- 1
s_vector[regexpr("^v.*Mismatch",names(s_vector))==1] <- 1
s_vector[regexpr("^v.*Match",names(s_vector))==1] <- 1

priorBushmakinModel= prior(designBushmakin,mu_mean=p_vector,mu_sd=s_vector,type="single")

BushmakinModel = make_emc(data,designBushmakin,type="single",fileName = 'samples.RData',prior=priorBushmakinModel)
BushmakinModel=fit(BushmakinModel,max_tries = 50, cores_per_chain=4, cores_for_chains=3, 
    stop_criteria = list(
    preburn = list(iter = 10), burn = list(mean_gd = 2.5), adapt = list(min_unique = 20),
    sample = list(iter = 25, max_gd = 2)), verbose = TRUE, verbose_progress=TRUE, fileName = 'samples.RData')
save.image("BushmakinModel.RData")

# Pred Data
pred = predict(BushmakinModel)
## Plotting
quantiles <- c(.1,.5,.9)#seq(from=0.05,to=0.95,by=0.05)
design_cells=data.frame(label=factor(levels(interaction(data$S,data$LogicalRule))))
pred$Correct=correctfun(pred)
data$Correct=correctfun(data)
pred=pred[pred$rt<4,]

## Densities
density_data_all <- purrr::map_dfr(design_cells$label, function(lbl) {
  parts <- strsplit(as.character(lbl), "\\.")[[1]]
  names(parts) <- c("S", "LogicalRule")
  
  prepare_density_data(
    data_list = list("Data"=data, "Recovered"=pred),
    correct_col = "Correct",
    factors = list(S = parts["S"], LogicalRule = parts["LogicalRule"]),
    design_label = lbl
  )
})  %>%
  separate(DesignCell, into = c("S", "LogicalRule"), sep = "\\.",remove = FALSE) %>%
  mutate(S=factor(S,levels=c("DT","UT","LT","NT")),LogicalRule=factor(LogicalRule,levels=c("AND","OR"))) %>%
  arrange(LogicalRule,S) %>%
  as.data.frame()

# Densities
ggplot(density_data_all, aes(x = RT, color = Source, fill = Source)) +
  geom_density(alpha = 0.3) +
  facet_wrap_paginate(Subset ~ DesignCell,ncol=4,nrow=4,page=1) +
  scale_color_manual(values = c("Data" = "red", "Recovered" = "blue")) +
  scale_fill_manual(values = c("Data" = "red", "Recovered" = "blue")) +
  theme_minimal() +
  labs(title = "RT Density: Data vs Recovered",
       x = "RT (seconds)", y = "Density")
# 
# Extract Quantiles by Design Cell x PM_trial x Correct
plotdata <- purrr::map_dfr(design_cells$label, function(lbl) {
  parts <- strsplit(as.character(lbl), "\\.")[[1]]
  names(parts) <- c("S", "LogicalRule")
  
  compute_quantile_points_by_subject(
    data_list = list("Data"=data, "Recovered"=pred),
    factors = list(S = parts["S"], LogicalRule = parts["LogicalRule"]),
    design_label = lbl,
    correct_col = "Correct",
    quantiles=quantiles,
    group_level = FALSE
  )
}) %>%
  group_by(Quantile, Source, Subset, DesignCell) %>%
  summarise(
    RT_mean = mean(RT, na.rm = TRUE),
    RT_se = sd(RT, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  separate(DesignCell, into = c("S", "LogicalRule"), sep = "\\.",remove = FALSE) %>%
  mutate(S=factor(S,levels=c("DT","UT","LT","NT")),LogicalRule=factor(LogicalRule,levels=c("AND","OR"))) %>%
  arrange(LogicalRule,S) %>%
  as.data.frame()
# Extract Quantiles by Design Cell x PM_trial x Correct
accuracy_data <- purrr::map_dfr(design_cells$label, function(lbl) {
  parts <- strsplit(as.character(lbl), "\\.")[[1]]
  names(parts) <- c("S", "LogicalRule")
  
  compute_accuracy_by_subject(
    data_list = list("Data"=data, "Recovered"=pred),
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
  ) + ylim(0.3,1.2)

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