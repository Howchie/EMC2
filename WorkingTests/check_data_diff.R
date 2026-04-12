
library(EMC2)

run_lba_data_check <- function(n_trials = 10000, seed = 123) {
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
  dat_inf <- make_data(p_vector, designLBA, n_trials = n_trials, TC=list(UC = Inf))
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  dat_100 <- make_data(p_vector, designLBA, n_trials = n_trials, TC=list(UC = 100))
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  dat_300 <- make_data(p_vector, designLBA, n_trials = n_trials, TC=list(UC = 300))

  cat("--- Data Comparison ---\n")
  cat("Max RT in dat_inf: ", max(dat_inf$rt, na.rm=TRUE), "\n")
  cat("Number of Inf RTs in dat_inf: ", sum(is.infinite(dat_inf$rt)), "\n")
  cat("Number of Inf RTs in dat_100: ", sum(is.infinite(dat_100$rt)), "\n")
  cat("Number of Inf RTs in dat_300: ", sum(is.infinite(dat_300$rt)), "\n")
  
  cat("Min RT in dat_inf: ", min(dat_inf$rt, na.rm=TRUE), "\n")
  
  cat("--- Seed Comparison ---\n")
  for (s in c(123, 456, 789)) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(s)
    dat_raw <- make_data(p_vector, designLBA, n_trials = 10000, TC=list(UC = Inf))
    cat("Seed ", s, ": Max RT=", max(dat_raw$rt, na.rm=TRUE), " Min RT=", min(dat_raw$rt, na.rm=TRUE), "\n")
  }

  # Check LL for true pars
  dadm_inf <- EMC2:::design_model(dat_inf, designLBA, rt_resolution = 1/60)
  dadm_100 <- EMC2:::design_model(dat_100, designLBA, rt_resolution = 1/60)
  dadm_300 <- EMC2:::design_model(dat_300, designLBA, rt_resolution = 1/60)
  
  dadm_inf <- EMC2:::.cache_ll_data_attrs(dadm_inf)
  dadm_100 <- EMC2:::.cache_ll_data_attrs(dadm_100)
  dadm_300 <- EMC2:::.cache_ll_data_attrs(dadm_300)
  
  ll_inf <- EMC2:::calc_ll_manager(t(as.matrix(p_vector)), dadm_inf, designLBA$model)
  ll_100 <- EMC2:::calc_ll_manager(t(as.matrix(p_vector)), dadm_100, designLBA$model)
  ll_300 <- EMC2:::calc_ll_manager(t(as.matrix(p_vector)), dadm_300, designLBA$model)
  
  cat("\n--- LL Comparison ---\n")
  cat("LL Inf: ", ll_inf, "\n")
  cat("LL 100: ", ll_100, "\n")
  cat("LL 300: ", ll_300, "\n")
  cat("Diff LL (100 - Inf): ", ll_100 - ll_inf, "\n")
  cat("Diff LL (300 - Inf): ", ll_300 - ll_inf, "\n")

  cat("\n--- Attribute Check ---\n")
  cat("emc2_all_finite_trials Inf: ", isTRUE(attr(dadm_inf, "emc2_all_finite_trials")), "\n")
  cat("emc2_all_finite_trials 100: ", isTRUE(attr(dadm_100, "emc2_all_finite_trials")), "\n")
  cat("emc2_all_finite_trials 300: ", isTRUE(attr(dadm_300, "emc2_all_finite_trials")), "\n")
  
  # Also check if they have partition attributes
  cat("has finite_rt_mask Inf: ", !is.null(attr(dadm_inf, "finite_rt_mask")), "\n")
  cat("has finite_rt_mask 100: ", !is.null(attr(dadm_100, "finite_rt_mask")), "\n")
  cat("has finite_rt_mask 300: ", !is.null(attr(dadm_300, "finite_rt_mask")), "\n")
}

run_lba_data_check()
