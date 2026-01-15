matchfun <- function(d) dplyr::case_when(d$S =="NT" & d$lR=="n_A" ~ TRUE,
                                         d$S =="NT" & d$lR=="n_B" ~ TRUE,
                                         d$S =="UT" & d$lR=="A" ~ TRUE,
                                         d$lR=="n_A_flip" ~ TRUE,
                                         d$S =="UT" & d$lR=="n_B" ~ TRUE,
                                         d$S =="LT" & d$lR=="n_A" ~ TRUE,
                                         d$S =="LT" & d$lR=="B" ~ TRUE,
                                         d$lR=="n_B_flip" ~ TRUE,
                                         d$S =="DT" & d$lR=="A" ~ TRUE,
                                         d$S =="DT" & d$lR=="B" ~ TRUE,
                                         TRUE ~ FALSE)
correctfun <- function(d) dplyr::case_when(d$LogicalRule=="OR" & d$S =="NT" & d$R=="no" ~ TRUE,
                                           d$LogicalRule=="OR" & !(d$S =="NT") & d$R=="yes" ~ TRUE,
                                           d$LogicalRule=="AND" & d$S =="DT" & d$R=="yes" ~ TRUE,
                                           d$LogicalRule=="AND" & !(d$S =="DT") & d$R=="no" ~ TRUE,
                                           d$LogicalRule=="ID" & d$S==d$R ~ TRUE,
                                           TRUE ~ FALSE)
stimLM_Contrast <- function(d) factor(dplyr::case_when(d$lM =="TRUE" & (d$lR=="A"|d$lR=="B") & d$S=="DT"~ "Y_DT",
                                                       d$lM =="TRUE" & (d$lR=="A"|d$lR=="B") & (d$S=="UT"|d$S=="LT")~ "Y_ST",
                                                       d$lM =="TRUE" & (d$lR=="n_A"|d$lR=="n_B")~ "N",
                                                       TRUE ~ "Mismatch"),levels=c("Mismatch","N","Y_ST","Y_DT"))
RateYNMatchMismatch_Contrast <- function(d) factor(dplyr::case_when(
  d$LogicalRule=="OR" & d$lM =="TRUE" & (d$lR=="A"|d$lR=="B")~ "MatchYP_OR",
  d$LogicalRule=="OR" & d$lM =="FALSE" & (d$lR=="A"|d$lR=="B")~ "MismatchYA_OR",
  d$LogicalRule=="OR" & d$lM =="TRUE" & (d$lR=="n_A"|d$lR=="n_B") ~ "MatchNA_OR",
  d$LogicalRule=="OR" & (d$lR=="n_A_flip"|d$lR=="n_B_flip") ~ "MatchNA_OR",
  d$LogicalRule=="OR" & d$lM =="FALSE" & (d$lR=="n_A"|d$lR=="n_B")~ "MismatchNP_OR",
  d$LogicalRule=="AND" & d$lM =="TRUE" & (d$lR=="A"|d$lR=="B")~ "MatchYP_AND",
  d$LogicalRule=="AND" & d$lM =="FALSE" & (d$lR=="A"|d$lR=="B")~ "MismatchYA_AND",
  d$LogicalRule=="AND" & d$lM =="TRUE" & (d$lR=="n_A"|d$lR=="n_B")~ "MatchNA_AND",
  d$LogicalRule=="AND" & (d$lR=="n_A_flip"|d$lR=="n_B_flip") ~ "MatchNA_AND",
  d$LogicalRule=="AND" & d$lM =="FALSE" & (d$lR=="n_A"|d$lR=="n_B")~ "MismatchNP_AND",
  d$lR=="guess" ~"Guess"))
RateYNMatchMismatch_ContrastShared <- function(d) factor(dplyr::case_when(
  d$lM =="TRUE" & (d$lR=="A"|d$lR=="B")~ "MatchYP",
  d$lM =="FALSE" & (d$lR=="A"|d$lR=="B")~ "MismatchYA",
  d$lM =="TRUE" & (d$lR=="n_A"|d$lR=="n_B")~ "MatchNA",
  d$lM =="FALSE" & (d$lR=="n_A"|d$lR=="n_B")~ "MismatchNP",
  d$lR=="guess" ~"Guess"))
Bias_Contrast <- function(d) factor(dplyr::case_when(
  d$LogicalRule=="OR" & d$lR=="A"~ "Bias",
  d$LogicalRule=="OR" & d$lR=="B"~ "Bias",
  d$LogicalRule=="OR" & d$lR=="n_A"~ "Other",
  d$LogicalRule=="OR" & d$lR=="n_B"~ "Other",
  d$LogicalRule=="AND" & d$lR=="A"~ "Other",
  d$LogicalRule=="AND" & d$lR=="B"~ "Other",
  d$LogicalRule=="AND" & d$lR=="n_A"~ "Bias",
  d$LogicalRule=="AND" & d$lR=="n_B"~ "Bias"))
Threshold_Location_Contrast <- function(d) factor(dplyr::case_when(
  d$lR=="A"~ "UT",
  d$lR=="B"~ "LT",
  d$lR=="n_A"~ "UT",
  d$lR=="n_B"~ "LT",
  d$lR=="n_A_flip"~ "UT",
  d$lR=="n_B_flip"~ "LT",
  d$lR=="guess"~"guess"
))
Capacity_Contrast <- function(d) factor(dplyr::case_when(d$lM =="TRUE" & (d$lR=="A"|d$lR=="B") & d$S=="DT"~ "DT",
                                                         TRUE ~ "Base"))
nTargets_Contrast <- function(d) factor(dplyr::case_when(d$S=="DT"~ "2T",
                                                         (d$S=="UT"|d$S=="LT")~ "1T",
                                                         d$S =="NT"~ "0T"),levels=c("0T","1T","2T"))
YN_Contrast <- function(d) factor(dplyr::case_when((d$lR=="A"|d$lR=="B")~ "Y",
                                                   d$lR=="n_A"|d$lR=="n_B"|d$lR=="n_A_flip"|d$lR=="n_B_flip"~ "N",
                                                   d$lR=="guess" ~"Guess"))


Rlevels <- c("A", "B", "n_A", "n_B", "guess","n_A_flip","n_B_flip")
Conditions <- c("AND", "OR")
Slevels = c("DT", "UT", "LT","NT")

drift_df <- expand.grid(
  lR = Rlevels,
  S = Slevels,
  LogicalRule = Conditions,
  stringsAsFactors = FALSE
)

threshold_df <- expand.grid(
  lR = Rlevels,
  S = Slevels,
  LogicalRule = Conditions,
  stringsAsFactors = FALSE
)

sv_df <- expand.grid(
  lR = Rlevels,
  S = Slevels,
  LogicalRule = Conditions,
  stringsAsFactors = FALSE
)

A_df <- expand.grid(
  lR = Rlevels,
  S = Slevels,
  LogicalRule = Conditions,
  stringsAsFactors = FALSE
)

# Combine to create consistent row labels
drift_df$cell_label <- apply(drift_df, 1, function(x) paste(x, collapse = "."))
drift_cells <- drift_df$cell_label
threshold_df$cell_label <- apply(threshold_df, 1, function(x) paste(x, collapse = "."))
threshold_cells <- threshold_df$cell_label
sv_df$cell_label <- apply(sv_df, 1, function(x) paste(x, collapse = "."))
sv_cells <- sv_df$cell_label
A_df$cell_label <- apply(A_df, 1, function(x) paste(x, collapse = "."))
A_cells <- A_df$cell_label

drift_df_reduced = drift_df
threshold_df_reduced = threshold_df

sv_names  =  c("matchY_UT", "matchN_UT", "mismatchY_UT", "mismatchN_UT",
               "matchY_LT", "matchN_LT", "mismatchY_LT", "mismatchN_LT","guess")
rate_names  =  c("OR_quality_P.UT", "OR_quantity_P.UT", "OR_quality_P.LT", "OR_quantity_P.LT", 
                 "OR_quality_A.UT", "OR_quantity_A.UT", "OR_quality_A.LT", "OR_quantity_A.LT",
                 "AND_quality_P.UT", "AND_quantity_P.UT", "AND_quality_P.LT", "AND_quantity_P.LT", 
                 "AND_quality_A.UT", "AND_quantity_A.UT", "AND_quality_A.LT", "AND_quantity_A.LT","guess")
rate_names_reduced  =  c("match_P.UT", "mismatch_P.UT", "match_P.LT", "mismatch_P.LT", 
                         "match_A.UT", "mismatch_A.UT", "match_A.LT", "mismatch_A.LT",
                         "guess")
threshold_names <- c("OR_Y_UT","OR_N_UT","OR_Y_LT","OR_N_LT","AND_Y_UT","AND_N_UT","AND_Y_LT","AND_N_LT","guess")
threshold_names_reduced <- c("Y_UT","N_UT","Y_LT","N_LT","guess")
# A_names = c("UT_bias","UT_other","LT_bias","LT_other","guess")
A_names = c("UT","LT","guess")

# Add columns to the cell_df -- this lets us use df operations to set the contrasts
for (p in rate_names) {
  drift_df[[p]] <- 0
}
for (p in rate_names_reduced) {
  drift_df_reduced[[p]] <- 0
}
for (p in sv_names) {
  sv_df[[p]] <- 0
}
for (p in threshold_names) {
  threshold_df[[p]] <- 0
}
for (p in threshold_names_reduced) {
  threshold_df_reduced[[p]] <- 0
}
for (p in A_names) {
  A_df[[p]] <- 0
}

drift_df = drift_df %>%
  mutate(lM=matchfun(drift_df),
         OR_quality_P.UT=case_when(LogicalRule=="OR" &  (S=="DT"|S=="UT") & (lR=="A"|lR=="n_A_flip") & lM==TRUE ~ 0.5,LogicalRule=="OR" & (S=="DT"|S=="UT") & lR=="n_A" & lM==FALSE ~ -0.5,TRUE ~ 0),
         OR_quantity_P.UT=case_when(LogicalRule=="OR" &  (S=="DT"|S=="UT") & (lR=="n_A"|lR=="n_A_flip"| lR=="A") ~ 0.5, TRUE ~ 0),
         OR_quality_P.LT=case_when(LogicalRule=="OR" &  (S=="DT"|S=="LT") & (lR=="B"|lR=="n_B_flip") & lM==TRUE ~ 0.5,LogicalRule=="OR" &(S=="DT"|S=="LT") & lR=="n_B" & lM==FALSE ~ -0.5,TRUE ~ 0),
         OR_quantity_P.LT=case_when(LogicalRule=="OR" &  (S=="DT"|S=="LT") & (lR=="n_B"|lR=="n_B_flip"| lR=="B") ~ 0.5, TRUE ~ 0),
         OR_quality_A.UT=case_when(LogicalRule=="OR" &  (S=="LT"|S=="NT") & (lR=="n_A"|lR=="n_A_flip") & lM==TRUE ~ 0.5,LogicalRule=="OR" &(S=="LT"|S=="NT") & lR=="A" & lM==FALSE ~ -0.5, TRUE ~ 0),
         OR_quantity_A.UT=case_when(LogicalRule=="OR" &  (S=="LT"|S=="NT") & (lR=="n_A"|lR=="n_A_flip"| lR=="A") ~ 0.5,TRUE ~ 0),
         OR_quality_A.LT=case_when(LogicalRule=="OR" &  (S=="UT"|S=="NT") & (lR=="n_B"|lR=="n_B_flip") & lM==TRUE ~ 0.5,LogicalRule=="OR" & (S=="UT"|S=="NT")& lR=="B" & lM==FALSE ~ -0.5, TRUE ~ 0),
         OR_quantity_A.LT=case_when(LogicalRule=="OR" &  (S=="UT"|S=="NT") & (lR=="n_B"|lR=="n_B_flip"| lR=="B") ~ 0.5, TRUE ~ 0),
         AND_quality_P.UT=case_when(LogicalRule=="AND" &  (S=="DT"|S=="UT") & (lR=="A"|lR=="n_A_flip") & lM==TRUE ~ 0.5,LogicalRule=="AND" & (S=="DT"|S=="UT") & lR=="n_A" & lM==FALSE ~ -0.5,TRUE ~ 0),
         AND_quantity_P.UT=case_when(LogicalRule=="AND" & (S=="DT"|S=="UT") & (lR=="n_A"|lR=="n_A_flip"| lR=="A") ~ 0.5, TRUE ~ 0),
         AND_quality_P.LT=case_when(LogicalRule=="AND" & (S=="DT"|S=="LT") & (lR=="B"|lR=="n_B_flip") & lM==TRUE ~ 0.5,LogicalRule=="AND" & (S=="DT"|S=="LT") & lR=="n_B" & lM==FALSE ~ -0.5,TRUE ~ 0),
         AND_quantity_P.LT=case_when(LogicalRule=="AND" & (S=="DT"|S=="LT") & (lR=="n_B"|lR=="n_B_flip"| lR=="B") ~ 0.5, TRUE ~ 0),
         AND_quality_A.UT=case_when(LogicalRule=="AND" & (S=="LT"|S=="NT") & (lR=="n_A"|lR=="n_A_flip") & lM==TRUE ~ 0.5,LogicalRule=="AND" & (S=="LT"|S=="NT") & lR=="A" & lM==FALSE ~ -0.5, TRUE ~ 0),
         AND_quantity_A.UT=case_when(LogicalRule=="AND" & (S=="LT"|S=="NT") & (lR=="n_A"|lR=="n_A_flip"| lR=="A") ~ 0.5,TRUE ~ 0),
         AND_quality_A.LT=case_when(LogicalRule=="AND" & (S=="UT"|S=="NT") & (lR=="n_B"|lR=="n_B_flip") & lM==TRUE ~ 0.5,LogicalRule=="AND" & (S=="UT"|S=="NT")& lR=="B" & lM==FALSE ~ -0.5, TRUE ~ 0),
         AND_quantity_A.LT=case_when(LogicalRule=="AND" & (S=="UT"|S=="NT") & (lR=="n_B"|lR=="n_B_flip"| lR=="B") ~ 0.5, TRUE ~ 0),
         guess = case_when(lR=="guess"~ 1,TRUE ~ 0)
  )

drift_df_reduced = drift_df_reduced %>%
  mutate(lM=matchfun(drift_df),
         match_P.UT=case_when((S=="DT"|S=="UT") & (lR=="A"|lR=="n_A_flip") & lM==TRUE ~ 1,TRUE ~ 0),
         mismatch_P.UT=case_when((S=="DT"|S=="UT") & (lR=="n_A") ~ 1, TRUE ~ 0),
         match_P.LT=case_when((S=="DT"|S=="LT") & (lR=="B"|lR=="n_B_flip") & lM==TRUE ~ 1,TRUE ~ 0),
         mismatch_P.LT=case_when((S=="DT"|S=="LT") & (lR=="n_B") ~ 1, TRUE ~ 0),
         match_A.UT=case_when((S=="LT"|S=="NT") & (lR=="n_A"|lR=="n_A_flip") & lM==TRUE ~ 1, TRUE ~ 0),
         mismatch_A.UT=case_when((S=="LT"|S=="NT") & (lR=="A") ~ 1,TRUE ~ 0),
         match_A.LT=case_when((S=="UT"|S=="NT") & (lR=="n_B"|lR=="n_B_flip") & lM==TRUE ~ 1, TRUE ~ 0),
         mismatch_A.LT=case_when((S=="UT"|S=="NT") & (lR=="B") ~ 1, TRUE ~ 0),
         guess = case_when(lR=="guess"~ 1,TRUE ~ 0)
  )

sv_df = sv_df %>%
  mutate(lM=matchfun(sv_df),
         matchY_UT=case_when((lR=="A")  & lM==TRUE ~ 1, TRUE ~ 0),
         matchN_UT=case_when((lR=="n_A"|lR=="n_A_flip") & lM==TRUE ~ 1, TRUE ~ 0),
         mismatchY_UT=case_when((lR=="A") & lM==FALSE ~ 1, TRUE ~ 0),
         mismatchN_UT=case_when((lR=="n_A") & lM==FALSE ~ 1,TRUE ~ 0),
         matchY_LT=case_when((lR=="B")  & lM==TRUE ~ 1, TRUE ~ 0),
         matchN_LT=case_when((lR=="n_B"|lR=="n_B_flip") & lM==TRUE ~ 1, TRUE ~ 0),
         mismatchY_LT=case_when((lR=="B") & lM==FALSE ~ 1, TRUE ~ 0),
         mismatchN_LT=case_when((lR=="n_B") & lM==FALSE ~ 1,TRUE ~ 0),
         guess = case_when(lR=="guess"~ 1,TRUE ~ 0)
  )

# A_df = A_df %>%
#   mutate(UT_bias=case_when((lR=="A" & LogicalRule=="OR") | ((lR=="n_A"|lR=="n_A_flip") & LogicalRule=="AND") ~ 1, TRUE ~ 0),
#          UT_other=case_when(((lR=="n_A"|lR=="n_A_flip") & LogicalRule=="OR") | (lR=="A" & LogicalRule=="AND") ~ 1, TRUE ~ 0),
#          LT_bias=case_when((lR=="B" & LogicalRule=="OR") | ((lR=="n_B"|lR=="n_B_flip") & LogicalRule=="AND") ~ 1, TRUE ~ 0),
#          LT_other=case_when(((lR=="n_B"|lR=="n_B_flip") & LogicalRule=="OR") | (lR=="B" & LogicalRule=="AND") ~ 1,TRUE ~ 0),
#          guess = case_when(lR=="guess"~ 1,TRUE ~ 0)
#   )

A_df = A_df %>%
  mutate(UT=case_when((lR=="A"|lR=="n_A"|lR=="n_A_flip") ~ 1, TRUE ~ 0),
         LT=case_when((lR=="B"|lR=="n_B"|lR=="n_B_flip") ~ 1,TRUE ~ 0),
         guess = case_when(lR=="guess"~ 1,TRUE ~ 0)
  )


threshold_df = threshold_df %>%
  mutate(OR_Y_UT=case_when(LogicalRule=="OR" & lR=="A" ~ 1, TRUE ~ 0),
         OR_N_UT=case_when(LogicalRule=="OR" & (lR=="n_A"|lR=="n_A_flip") ~ 1, TRUE ~ 0),
         OR_Y_LT=case_when(LogicalRule=="OR" & lR=="B" ~ 1, TRUE ~ 0),
         OR_N_LT=case_when(LogicalRule=="OR" & (lR=="n_B"|lR=="n_B_flip") ~ 1,TRUE ~ 0),
         AND_Y_UT=case_when(LogicalRule=="AND" & lR=="A" ~ 1, TRUE ~ 0),
         AND_N_UT=case_when(LogicalRule=="AND" & (lR=="n_A"|lR=="n_A_flip") ~ 1, TRUE ~ 0),
         AND_Y_LT=case_when(LogicalRule=="AND" & lR=="B" ~ 1, TRUE ~ 0),
         AND_N_LT=case_when(LogicalRule=="AND" & (lR=="n_B"|lR=="n_B_flip") ~ 1,TRUE ~ 0),
         guess = case_when(lR=="guess"~ 1,TRUE ~ 0)
  )

threshold_df_reduced = threshold_df_reduced %>%
  mutate(Y_UT=case_when(lR=="A" ~ 1, TRUE ~ 0),
         N_UT=case_when((lR=="n_A"|lR=="n_A_flip") ~ 1, TRUE ~ 0),
         Y_LT=case_when(lR=="B" ~ 1, TRUE ~ 0),
         N_LT=case_when((lR=="n_B"|lR=="n_B_flip") ~ 1,TRUE ~ 0),
         guess = case_when(lR=="guess"~ 1,TRUE ~ 0)
  )



# select just the drift contrast columns, used for EMC2 design
rate_design = drift_df %>%
  dplyr::select(all_of(rate_names)) %>%
  as.matrix()
rownames(rate_design) <- drift_df$cell_label

rate_design_reduced = drift_df_reduced %>%
  dplyr::select(all_of(rate_names_reduced)) %>%
  as.matrix()
rownames(rate_design_reduced) <- drift_df_reduced$cell_label

matchfun <- function(d) dplyr::case_when(d$S =="NT" & d$lR=="n_A" ~ TRUE,
                                         d$S =="NT" & d$lR=="n_B" ~ TRUE,
                                         d$S =="UT" & d$lR=="A" ~ TRUE,
                                         d$lR=="n_A_flip" ~ TRUE,
                                         d$S =="UT" & d$lR=="n_B" ~ TRUE,
                                         d$S =="LT" & d$lR=="n_A" ~ TRUE,
                                         d$S =="LT" & d$lR=="B" ~ TRUE,
                                         d$lR=="n_B_flip" ~ TRUE,
                                         d$S =="DT" & d$lR=="A" ~ TRUE,
                                         d$S =="DT" & d$lR=="B" ~ TRUE,
                                         TRUE ~ FALSE)
correctfun <- function(d) dplyr::case_when(d$LogicalRule=="OR" & d$S =="NT" & d$R=="no" ~ TRUE,
                                           d$LogicalRule=="OR" & !d$S =="NT" & d$R=="yes" ~ TRUE,
                                           d$LogicalRule=="AND" & d$S =="DT" & d$R=="yes" ~ TRUE,
                                           d$LogicalRule=="AND" & !d$S =="DT" & d$R=="no" ~ TRUE,
                                           TRUE ~ FALSE)
stimLM_Contrast <- function(d) factor(dplyr::case_when(d$lM =="TRUE" & (d$lR=="A"|d$lR=="B") & d$S=="DT"~ "Y_DT",
                                                       d$lM =="TRUE" & (d$lR=="A"|d$lR=="B") & (d$S=="UT"|d$S=="LT")~ "Y_ST",
                                                       d$lM =="TRUE" & (d$lR=="n_A"|d$lR=="n_B")~ "N",
                                                       TRUE ~ "Mismatch"),levels=c("Mismatch","N","Y_ST","Y_DT"))
RateYNMatchMismatch_Contrast <- function(d) factor(dplyr::case_when(
  d$LogicalRule=="OR" & d$lM =="TRUE" & (d$lR=="A"|d$lR=="B")~ "MatchYP_OR",
  d$LogicalRule=="OR" & d$lM =="FALSE" & (d$lR=="A"|d$lR=="B")~ "MismatchYA_OR",
  d$LogicalRule=="OR" & d$lM =="TRUE" & (d$lR=="n_A"|d$lR=="n_B") ~ "MatchNA_OR",
  d$LogicalRule=="OR" & (d$lR=="n_A_flip"|d$lR=="n_B_flip") ~ "MatchNA_OR",
  d$LogicalRule=="OR" & d$lM =="FALSE" & (d$lR=="n_A"|d$lR=="n_B")~ "MismatchNP_OR",
  d$LogicalRule=="AND" & d$lM =="TRUE" & (d$lR=="A"|d$lR=="B")~ "MatchYP_AND",
  d$LogicalRule=="AND" & d$lM =="FALSE" & (d$lR=="A"|d$lR=="B")~ "MismatchYA_AND",
  d$LogicalRule=="AND" & d$lM =="TRUE" & (d$lR=="n_A"|d$lR=="n_B")~ "MatchNA_AND",
  d$LogicalRule=="AND" & (d$lR=="n_A_flip"|d$lR=="n_B_flip") ~ "MatchNA_AND",
  d$LogicalRule=="AND" & d$lM =="FALSE" & (d$lR=="n_A"|d$lR=="n_B")~ "MismatchNP_AND",
  d$lR=="guess" ~"Guess"))
RateYNMatchMismatch_ContrastShared <- function(d) factor(dplyr::case_when(
  d$lM =="TRUE" & (d$lR=="A"|d$lR=="B")~ "MatchYP",
  d$lM =="FALSE" & (d$lR=="A"|d$lR=="B")~ "MismatchYA",
  d$lM =="TRUE" & (d$lR=="n_A"|d$lR=="n_B")~ "MatchNA",
  d$lM =="FALSE" & (d$lR=="n_A"|d$lR=="n_B")~ "MismatchNP",
  d$lR=="guess" ~"Guess"))
Bias_Contrast <- function(d) factor(dplyr::case_when(
  d$LogicalRule=="OR" & d$lR=="A"~ "Bias",
  d$LogicalRule=="OR" & d$lR=="B"~ "Bias",
  d$LogicalRule=="OR" & d$lR=="n_A"~ "Other",
  d$LogicalRule=="OR" & d$lR=="n_B"~ "Other",
  d$LogicalRule=="OR" & d$lR=="n_A_flip"~ "Other",
  d$LogicalRule=="OR" & d$lR=="n_B_flip"~ "Other",
  d$LogicalRule=="AND" & d$lR=="A"~ "Other",
  d$LogicalRule=="AND" & d$lR=="B"~ "Other",
  d$LogicalRule=="AND" & d$lR=="n_A"~ "Bias",
  d$LogicalRule=="AND" & d$lR=="n_B"~ "Bias",
  d$LogicalRule=="AND" & d$lR=="n_A_flip"~ "Bias",
  d$LogicalRule=="AND" & d$lR=="n_B_flip"~ "Bias",
  d$lR=="guess"~"guess"))
Location_Contrast <- function(d) factor(dplyr::case_when(
  d$lR=="A"~ "UT",
  d$lR=="B"~ "LT",
  d$lR=="n_A"~ "UT",
  d$lR=="n_B"~ "LT",
  d$lR=="n_A_flip"~ "UT",
  d$lR=="n_B_flip"~ "LT",
  d$lR=="guess"~"guess"
))
Capacity_Contrast <- function(d) factor(dplyr::case_when(d$lM =="TRUE" & (d$lR=="A"|d$lR=="B") & d$S=="DT"~ "DT",
                                                         TRUE ~ "Base"))
nTargets_Contrast <- function(d) factor(dplyr::case_when(d$S=="DT"~ "2T",
                                                         (d$S=="UT"|d$S=="LT")~ "1T",
                                                         d$S =="NT"~ "0T"),levels=c("0T","1T","2T"))
YN_Contrast <- function(d) factor(dplyr::case_when((d$lR=="A"|d$lR=="B")~ "Y",
                                                   d$lR=="n_A"|d$lR=="n_B"|d$lR=="n_A_flip"|d$lR=="n_B_flip"~ "N",
                                                   d$lR=="guess" ~"Guess"))


Rlevels <- c("A", "B", "n_A", "n_B", "guess","n_A_flip","n_B_flip")
Conditions <- c("AND", "OR")
Slevels = c("DT", "UT", "LT","NT")

drift_df <- expand.grid(
  lR = Rlevels,
  S = Slevels,
  LogicalRule = Conditions,
  stringsAsFactors = FALSE
)

threshold_df <- expand.grid(
  lR = Rlevels,
  S = Slevels,
  LogicalRule = Conditions,
  stringsAsFactors = FALSE
)

sv_df <- expand.grid(
  lR = Rlevels,
  S = Slevels,
  LogicalRule = Conditions,
  stringsAsFactors = FALSE
)

A_df <- expand.grid(
  lR = Rlevels,
  S = Slevels,
  LogicalRule = Conditions,
  stringsAsFactors = FALSE
)

# Combine to create consistent row labels
drift_df$cell_label <- apply(drift_df, 1, function(x) paste(x, collapse = "."))
drift_cells <- drift_df$cell_label
threshold_df$cell_label <- apply(threshold_df, 1, function(x) paste(x, collapse = "."))
threshold_cells <- threshold_df$cell_label
sv_df$cell_label <- apply(sv_df, 1, function(x) paste(x, collapse = "."))
sv_cells <- sv_df$cell_label
A_df$cell_label <- apply(A_df, 1, function(x) paste(x, collapse = "."))
A_cells <- A_df$cell_label

drift_df_reduced = drift_df
threshold_df_reduced = threshold_df

# sv_names  =  c("matchY_UT", "matchN_UT", "mismatchY_UT", "mismatchN_UT",
#                "matchY_LT", "matchN_LT", "mismatchY_LT", "mismatchN_LT","guess")
sv_names  =  c("matchY_UT", "matchN_UT", 
               "matchY_LT", "matchN_LT", "mismatch", "guess")
rate_names  =  c("OR_quality_P.UT", "OR_quantity_P.UT", "OR_quality_P.LT", "OR_quantity_P.LT", 
                 "OR_quality_A.UT", "OR_quantity_A.UT", "OR_quality_A.LT", "OR_quantity_A.LT",
                 "AND_quality_P.UT", "AND_quantity_P.UT", "AND_quality_P.LT", "AND_quantity_P.LT", 
                 "AND_quality_A.UT", "AND_quantity_A.UT", "AND_quality_A.LT", "AND_quantity_A.LT","guess")
rate_names_reduced  =  c("match_P.UT", "mismatch_P.UT", "match_P.LT", "mismatch_P.LT", 
                         "match_A.UT", "mismatch_A.UT", "match_A.LT", "mismatch_A.LT",
                         "guess")
threshold_names <- c("OR_Y_UT","OR_N_UT","OR_Y_LT","OR_N_LT","AND_Y_UT","AND_N_UT","AND_Y_LT","AND_N_LT","guess")
threshold_names_reduced <- c("Y_UT","N_UT","Y_LT","N_LT","guess")
# A_names = c("UT_bias","UT_other","LT_bias","LT_other","guess")
A_names = c("UT","LT","guess")

# Add columns to the cell_df -- this lets us use df operations to set the contrasts
for (p in rate_names) {
  drift_df[[p]] <- 0
}
for (p in rate_names_reduced) {
  drift_df_reduced[[p]] <- 0
}
for (p in sv_names) {
  sv_df[[p]] <- 0
}
for (p in threshold_names) {
  threshold_df[[p]] <- 0
}
for (p in threshold_names_reduced) {
  threshold_df_reduced[[p]] <- 0
}
for (p in A_names) {
  A_df[[p]] <- 0
}

drift_df = drift_df %>%
  mutate(lM=matchfun(drift_df),
         OR_quality_P.UT=case_when(LogicalRule=="OR" &  (S=="DT"|S=="UT") & (lR=="A"|lR=="n_A_flip") & lM==TRUE ~ 0.5,LogicalRule=="OR" & (S=="DT"|S=="UT") & lR=="n_A" & lM==FALSE ~ -0.5,TRUE ~ 0),
         OR_quantity_P.UT=case_when(LogicalRule=="OR" &  (S=="DT"|S=="UT") & (lR=="n_A"|lR=="n_A_flip"| lR=="A") ~ 0.5, TRUE ~ 0),
         OR_quality_P.LT=case_when(LogicalRule=="OR" &  (S=="DT"|S=="LT") & (lR=="B"|lR=="n_B_flip") & lM==TRUE ~ 0.5,LogicalRule=="OR" &(S=="DT"|S=="LT") & lR=="n_B" & lM==FALSE ~ -0.5,TRUE ~ 0),
         OR_quantity_P.LT=case_when(LogicalRule=="OR" &  (S=="DT"|S=="LT") & (lR=="n_B"|lR=="n_B_flip"| lR=="B") ~ 0.5, TRUE ~ 0),
         OR_quality_A.UT=case_when(LogicalRule=="OR" &  (S=="LT"|S=="NT") & (lR=="n_A"|lR=="n_A_flip") & lM==TRUE ~ 0.5,LogicalRule=="OR" &(S=="LT"|S=="NT") & lR=="A" & lM==FALSE ~ -0.5, TRUE ~ 0),
         OR_quantity_A.UT=case_when(LogicalRule=="OR" &  (S=="LT"|S=="NT") & (lR=="n_A"|lR=="n_A_flip"| lR=="A") ~ 0.5,TRUE ~ 0),
         OR_quality_A.LT=case_when(LogicalRule=="OR" &  (S=="UT"|S=="NT") & (lR=="n_B"|lR=="n_B_flip") & lM==TRUE ~ 0.5,LogicalRule=="OR" & (S=="UT"|S=="NT")& lR=="B" & lM==FALSE ~ -0.5, TRUE ~ 0),
         OR_quantity_A.LT=case_when(LogicalRule=="OR" &  (S=="UT"|S=="NT") & (lR=="n_B"|lR=="n_B_flip"| lR=="B") ~ 0.5, TRUE ~ 0),
         AND_quality_P.UT=case_when(LogicalRule=="AND" &  (S=="DT"|S=="UT") & (lR=="A"|lR=="n_A_flip") & lM==TRUE ~ 0.5,LogicalRule=="AND" & (S=="DT"|S=="UT") & lR=="n_A" & lM==FALSE ~ -0.5,TRUE ~ 0),
         AND_quantity_P.UT=case_when(LogicalRule=="AND" & (S=="DT"|S=="UT") & (lR=="n_A"|lR=="n_A_flip"| lR=="A") ~ 0.5, TRUE ~ 0),
         AND_quality_P.LT=case_when(LogicalRule=="AND" & (S=="DT"|S=="LT") & (lR=="B"|lR=="n_B_flip") & lM==TRUE ~ 0.5,LogicalRule=="AND" & (S=="DT"|S=="LT") & lR=="n_B" & lM==FALSE ~ -0.5,TRUE ~ 0),
         AND_quantity_P.LT=case_when(LogicalRule=="AND" & (S=="DT"|S=="LT") & (lR=="n_B"|lR=="n_B_flip"| lR=="B") ~ 0.5, TRUE ~ 0),
         AND_quality_A.UT=case_when(LogicalRule=="AND" & (S=="LT"|S=="NT") & (lR=="n_A"|lR=="n_A_flip") & lM==TRUE ~ 0.5,LogicalRule=="AND" & (S=="LT"|S=="NT") & lR=="A" & lM==FALSE ~ -0.5, TRUE ~ 0),
         AND_quantity_A.UT=case_when(LogicalRule=="AND" & (S=="LT"|S=="NT") & (lR=="n_A"|lR=="n_A_flip"| lR=="A") ~ 0.5,TRUE ~ 0),
         AND_quality_A.LT=case_when(LogicalRule=="AND" & (S=="UT"|S=="NT") & (lR=="n_B"|lR=="n_B_flip") & lM==TRUE ~ 0.5,LogicalRule=="AND" & (S=="UT"|S=="NT")& lR=="B" & lM==FALSE ~ -0.5, TRUE ~ 0),
         AND_quantity_A.LT=case_when(LogicalRule=="AND" & (S=="UT"|S=="NT") & (lR=="n_B"|lR=="n_B_flip"| lR=="B") ~ 0.5, TRUE ~ 0),
         guess = case_when(lR=="guess"~ 1,TRUE ~ 0)
  )

drift_df_reduced = drift_df_reduced %>%
  mutate(lM=matchfun(drift_df),
         match_P.UT=case_when((S=="DT"|S=="UT") & (lR=="A"|lR=="n_A_flip") & lM==TRUE ~ 1,TRUE ~ 0),
         mismatch_P.UT=case_when((S=="DT"|S=="UT") & (lR=="n_A") ~ 1, TRUE ~ 0),
         match_P.LT=case_when((S=="DT"|S=="LT") & (lR=="B"|lR=="n_B_flip") & lM==TRUE ~ 1,TRUE ~ 0),
         mismatch_P.LT=case_when((S=="DT"|S=="LT") & (lR=="n_B") ~ 1, TRUE ~ 0),
         match_A.UT=case_when((S=="LT"|S=="NT") & (lR=="n_A"|lR=="n_A_flip") & lM==TRUE ~ 1, TRUE ~ 0),
         mismatch_A.UT=case_when((S=="LT"|S=="NT") & (lR=="A") ~ 1,TRUE ~ 0),
         match_A.LT=case_when((S=="UT"|S=="NT") & (lR=="n_B"|lR=="n_B_flip") & lM==TRUE ~ 1, TRUE ~ 0),
         mismatch_A.LT=case_when((S=="UT"|S=="NT") & (lR=="B") ~ 1, TRUE ~ 0),
         guess = case_when(lR=="guess"~ 1,TRUE ~ 0)
  )

# sv_df = sv_df %>%
#   mutate(lM=matchfun(sv_df),
#          matchY_UT=case_when((lR=="A")  & lM==TRUE ~ 1, TRUE ~ 0),
#          matchN_UT=case_when((lR=="n_A"|lR=="n_A_flip") & lM==TRUE ~ 1, TRUE ~ 0),
#          mismatchY_UT=case_when((lR=="A") & lM==FALSE ~ 1, TRUE ~ 0),
#          mismatchN_UT=case_when((lR=="n_A") & lM==FALSE ~ 1,TRUE ~ 0),
#          matchY_LT=case_when((lR=="B")  & lM==TRUE ~ 1, TRUE ~ 0),
#          matchN_LT=case_when((lR=="n_B"|lR=="n_B_flip") & lM==TRUE ~ 1, TRUE ~ 0),
#          mismatchY_LT=case_when((lR=="B") & lM==FALSE ~ 1, TRUE ~ 0),
#          mismatchN_LT=case_when((lR=="n_B") & lM==FALSE ~ 1,TRUE ~ 0),
#          guess = case_when(lR=="guess"~ 1,TRUE ~ 0)
#   )
sv_df = sv_df %>%
  mutate(lM=matchfun(sv_df),
         matchY_UT=case_when((lR=="A")  & lM==TRUE ~ 1, TRUE ~ 0),
         matchN_UT=case_when((lR=="n_A"|lR=="n_A_flip") & lM==TRUE ~ 1, TRUE ~ 0),
         matchY_LT=case_when((lR=="B")  & lM==TRUE ~ 1, TRUE ~ 0),
         matchN_LT=case_when((lR=="n_B"|lR=="n_B_flip") & lM==TRUE ~ 1, TRUE ~ 0),
         mismatch=case_when(lM==FALSE ~ 1, TRUE ~ 0),
         guess = case_when(lR=="guess"~ 1,TRUE ~ 0)
  )

# A_df = A_df %>%
#   mutate(UT_bias=case_when((lR=="A" & LogicalRule=="OR") | ((lR=="n_A"|lR=="n_A_flip") & LogicalRule=="AND") ~ 1, TRUE ~ 0),
#          UT_other=case_when(((lR=="n_A"|lR=="n_A_flip") & LogicalRule=="OR") | (lR=="A" & LogicalRule=="AND") ~ 1, TRUE ~ 0),
#          LT_bias=case_when((lR=="B" & LogicalRule=="OR") | ((lR=="n_B"|lR=="n_B_flip") & LogicalRule=="AND") ~ 1, TRUE ~ 0),
#          LT_other=case_when(((lR=="n_B"|lR=="n_B_flip") & LogicalRule=="OR") | (lR=="B" & LogicalRule=="AND") ~ 1,TRUE ~ 0),
#          guess = case_when(lR=="guess"~ 1,TRUE ~ 0)
#   )

A_df = A_df %>%
  mutate(UT=case_when((lR=="A"|lR=="n_A"|lR=="n_A_flip") ~ 1, TRUE ~ 0),
         LT=case_when((lR=="B"|lR=="n_B"|lR=="n_B_flip") ~ 1,TRUE ~ 0),
         guess = case_when(lR=="guess"~ 1,TRUE ~ 0)
  )


threshold_df = threshold_df %>%
  mutate(OR_Y_UT=case_when(LogicalRule=="OR" & lR=="A" ~ 1, TRUE ~ 0),
         OR_N_UT=case_when(LogicalRule=="OR" & (lR=="n_A"|lR=="n_A_flip") ~ 1, TRUE ~ 0),
         OR_Y_LT=case_when(LogicalRule=="OR" & lR=="B" ~ 1, TRUE ~ 0),
         OR_N_LT=case_when(LogicalRule=="OR" & (lR=="n_B"|lR=="n_B_flip") ~ 1,TRUE ~ 0),
         AND_Y_UT=case_when(LogicalRule=="AND" & lR=="A" ~ 1, TRUE ~ 0),
         AND_N_UT=case_when(LogicalRule=="AND" & (lR=="n_A"|lR=="n_A_flip") ~ 1, TRUE ~ 0),
         AND_Y_LT=case_when(LogicalRule=="AND" & lR=="B" ~ 1, TRUE ~ 0),
         AND_N_LT=case_when(LogicalRule=="AND" & (lR=="n_B"|lR=="n_B_flip") ~ 1,TRUE ~ 0),
         guess = case_when(lR=="guess"~ 1,TRUE ~ 0)
  )

threshold_df_reduced = threshold_df_reduced %>%
  mutate(Y_UT=case_when(lR=="A" ~ 1, TRUE ~ 0),
         N_UT=case_when((lR=="n_A"|lR=="n_A_flip") ~ 1, TRUE ~ 0),
         Y_LT=case_when(lR=="B" ~ 1, TRUE ~ 0),
         N_LT=case_when((lR=="n_B"|lR=="n_B_flip") ~ 1,TRUE ~ 0),
         guess = case_when(lR=="guess"~ 1,TRUE ~ 0)
  )



# select just the drift contrast columns, used for EMC2 design
rate_design = drift_df %>%
  dplyr::select(all_of(rate_names)) %>%
  as.matrix()
rownames(rate_design) <- drift_df$cell_label

rate_design_reduced = drift_df_reduced %>%
  dplyr::select(all_of(rate_names_reduced)) %>%
  as.matrix()
rownames(rate_design_reduced) <- drift_df_reduced$cell_label


sv_design = sv_df %>%
  dplyr::select(all_of(sv_names)) %>%
  as.matrix()
rownames(sv_design) <- sv_df$cell_label

A_design = A_df %>%
  dplyr::select(all_of(A_names)) %>%
  as.matrix()
rownames(A_design) <- A_df$cell_label

threshold_design = threshold_df %>%
  dplyr::select(all_of(threshold_names)) %>%
  as.matrix()
rownames(threshold_design) <- threshold_df$cell_label

threshold_design_reduced = threshold_df_reduced %>%
  dplyr::select(all_of(threshold_names_reduced)) %>%
  as.matrix()
rownames(threshold_design_reduced) <- threshold_df_reduced$cell_label
sv_design = sv_df %>%
  dplyr::select(all_of(sv_names)) %>%
  as.matrix()
rownames(sv_design) <- sv_df$cell_label

A_design = A_df %>%
  dplyr::select(all_of(A_names)) %>%
  as.matrix()
rownames(A_design) <- A_df$cell_label

threshold_design = threshold_df %>%
  dplyr::select(all_of(threshold_names)) %>%
  as.matrix()
rownames(threshold_design) <- threshold_df$cell_label

threshold_design_reduced = threshold_df_reduced %>%
  dplyr::select(all_of(threshold_names_reduced)) %>%
  as.matrix()
rownames(threshold_design_reduced) <- threshold_df_reduced$cell_label