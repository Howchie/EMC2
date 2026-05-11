library(testthat)
library(EMC2)

test_that("pointwise log-likelihood sums to total log-likelihood", {
  # Setup a simple model
  data(forstmann)
  dat <- forstmann[forstmann$subjects == unique(forstmann$subjects)[1],]
  dat$subjects <- droplevels(dat$subjects)
  design_LNR <- design(data = dat, model=LNR, 
                       formula=list(m~1, s~1, t0~1))
  emc <- make_emc(dat, design_LNR, type = "single")
  
  # Get some parameters
  model <- emc[[1]]$model
  p_types <- names(model()$p_types)
  p_vector <- setNames(rep(0, length(p_types)), p_types)
  p_vector["t0"] <- -1
  
  p_mat <- matrix(p_vector, nrow=1)
  colnames(p_mat) <- p_types
  
  # Pick a subject-level dadm
  dadm <- emc[[1]]$data[[1]]
  
  # Calculate total LL
  total_ll <- EMC2:::calc_ll_manager(p_mat, dadm, model)
  
  # Calculate pointwise LL
  pw_ll <- calc_ll_pw(p_mat, dadm, model)
  
  expect_equal(sum(pw_ll[1,]), as.numeric(total_ll), tolerance = 1e-10)
})
