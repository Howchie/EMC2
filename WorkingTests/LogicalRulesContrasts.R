funcs = list(mismatch = function(d) ifelse(d$lM==TRUE, 0, 1),
             match = function(d) ifelse(d$lM==TRUE, 1, 0),
             LocationA=function(d) ifelse(d$lR=="A"|d$lR=="n_A",1,0),
             LocationB=function(d) ifelse(d$lR=="B"|d$lR=="n_B",1,0),
             Yes = function(d) ifelse(d$lR=="A"|d$lR=="B",1,0),
             No = function(d) ifelse(d$lR=="n_A"|d$lR=="n_B",1,0),
             DT = function(d) dplyr::case_when(d$lM =="TRUE" & (d$lR=="A"|d$lR=="B") & d$S=="AB"~ 1, TRUE ~ 0),
             BaseY = function(d) dplyr::case_when(d$lM =="TRUE" & (d$lR=="A"|d$lR=="B")~ 1, TRUE ~ 0),
             BaseN = function(d) dplyr::case_when(d$lM =="TRUE" & (d$lR=="n_A"|d$lR=="n_B")~ 1, TRUE ~ 0),
             MatchY = function(d) dplyr::case_when(d$match & d$Yes ~ 1, TRUE ~ 0),
             MatchN = function(d) dplyr::case_when(d$match & d$No ~ 1, TRUE ~ 0)
)
matchfun <- function(d) dplyr::case_when(d$S =="NN" & d$lR=="n_A" ~ TRUE,
                                         d$S =="NN" & d$lR=="n_B" ~ TRUE,
                                         d$S =="AN" & d$lR=="A" ~ TRUE,
                                         d$S =="AN" & d$lR=="n_B" ~ TRUE,
                                         d$S =="NB" & d$lR=="n_A" ~ TRUE,
                                         d$S =="NB" & d$lR=="B" ~ TRUE,
                                         d$S =="AB" & d$lR=="A" ~ TRUE,
                                         d$S =="AB" & d$lR=="B" ~ TRUE,
                                         TRUE ~ FALSE)

correctfun <- function(d) {
  cond <- d$LogicalRule
  s <- d$S
  r <- d$R
  
  res <- logical(length(cond))
  
  if (any(idx_or <- cond == "OR")) {
    s_or <- s[idx_or]
    r_or <- r[idx_or]
    res[idx_or] <- (s_or == "NN" & r_or == "no") | (s_or != "NN" & r_or == "yes")
  }
  
  if (any(idx_and <- cond == "AND")) {
    s_and <- s[idx_and]
    r_and <- r[idx_and]
    res[idx_and] <- (s_and == "AB" & r_and == "yes") | (s_and != "AB" & r_and == "no")
  }
  
  if (any(idx_xor <- cond == "XOR")) {
    s_xor <- s[idx_xor]
    r_xor <- r[idx_xor]
    res[idx_xor] <- (s_xor %in% c("AN", "NB") & r_xor == "yes") | (s_xor %in% c("AB", "NN") & r_xor == "no")
  }
  
  if (any(idx_id <- cond == "ID")) {
    res[idx_id] <- (s[idx_id] == r[idx_id])
  }
  
  as.numeric(res)
}
