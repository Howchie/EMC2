
library(EMC2)
source("WorkingTests/test_likelihood_plotfuns_ah.R")

run_lba_test <- function(n_trials = 10000, UC = Inf, label = "test", seed = 123) {
  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)
  designLBA <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = LBA(posdrift=TRUE), 
    formula = list(B ~ 1, v ~ 0+lM, A ~ 1, t0 ~ 1, sv ~ 0+lM),
    constants = c(A = log(0.4))
  )
  p_vector <- sampled_pars(designLBA, doMap = FALSE)
  p_vector[["B"]] <- log(1.2)
  p_vector[["A"]] <- log(.4)
  p_vector[["t0"]] <- log(0.15)
  p_vector[["v_lMFALSE"]] <- .4
  p_vector[["v_lMTRUE"]] <- 1.2
  p_vector[["sv_lMFALSE"]] <- log(1.2)
  p_vector[["sv_lMTRUE"]] <- log(1)

  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  dat <- make_data(p_vector, designLBA, n_trials = n_trials, TC=list(UC = UC))
  cat("Min RT in dat: ", min(dat$rt, na.rm=TRUE), "\n")
  cat("Smallest 10 RTs:\n")
  print(head(sort(dat$rt), 10))
  cat("True t0 (natural): ", exp(p_vector["t0"]), "\n")
  
  emc <- make_emc(dat, designLBA, type = "single", rt_resolution = 1/60)
  
  emc <- fit(emc, stop_criteria = list(
    sample = list(iter = 1000, max_gd = 1.1)
  ), cores_for_chains = 3, cores_per_chain = 1, max_tries=1)
  
  cat("\nRecovery for ", label, ":\n")
  samples <- emc[[1]]$samples$alpha[,,emc[[1]]$samples$stage=="sample"]
  post_means <- apply(samples, 1, mean)
  print(rbind(True=p_vector, PostMean=post_means))
  
  dadm <- emc[[1]]$data[[1]]
  ll <- EMC2:::calc_ll_manager(t(as.matrix(p_vector)), dadm, designLBA$model)
  cat("Log-likelihood for true pars: ", ll, "\n")
  
  return(list(emc=emc, ll=ll))
}

cat("--- Running UC = Inf with Seed 123 ---\n")
res_123 <- run_lba_test(UC = Inf, label = "Seed_123", seed = 123)

cat("\n--- Running UC = Inf with Seed 456 ---\n")
res_456_inf <- run_lba_test(UC = Inf, label = "Seed_456_Inf", seed = 456)

cat("\n--- Running Seed 123 with rt_resolution = NULL ---\n")
run_lba_test_null <- function(n_trials = 10000, seed = 123) {
  # (re-copy logic to ensure it's clean)
  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)
  designLBA <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = matchfun,
    model = LBA(posdrift=TRUE), 
    formula = list(B ~ 1, v ~ 0+lM, A ~ 1, t0 ~ 1, sv ~ 0+lM),
    constants = c(A = log(0.4))
  )
  p_vector <- sampled_pars(designLBA, doMap = FALSE)
  p_vector[["B"]] <- log(1.2)
  p_vector[["A"]] <- log(.4)
  p_vector[["t0"]] <- log(0.15)
  p_vector[["v_lMFALSE"]] <- .4
  p_vector[["v_lMTRUE"]] <- 1.2
  p_vector[["sv_lMFALSE"]] <- log(1.2)
  p_vector[["sv_lMTRUE"]] <- log(1)

  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  dat <- make_data(p_vector, designLBA, n_trials = n_trials, TC=list(UC = Inf))
  
  emc <- make_emc(dat, designLBA, type = "single", rt_resolution = NULL)
  
  emc <- fit(emc, stop_criteria = list(
    sample = list(iter = 1000, max_gd = 1.1)
  ), cores_for_chains = 3, cores_per_chain = 1, max_tries=1)
  
  cat("\nRecovery for Seed 123 res NULL:\n")
  samples <- emc[[1]]$samples$alpha[,,emc[[1]]$samples$stage=="sample"]
  post_means <- apply(samples, 1, mean)
  print(rbind(True=p_vector, PostMean=post_means))
}

run_lba_test_null(seed = 123)
