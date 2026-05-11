library(testthat)
library(EMC2)

RNGkind("L'Ecuyer-CMRG")
set.seed(123)

test_that("pw_ll extraction with race models (LBA/LNR)", {
  data(forstmann)
  # Take a small subset for speed
  dat <- forstmann[forstmann$subjects == unique(forstmann$subjects)[1],]
  dat$subjects <- droplevels(dat$subjects)
  
  # LBA model has n_lR = 2 for this data
  design_LBA <- design(data = dat, model=LBA, 
                      formula=list(v~1, sv~1, B~1, A~1, t0~1))
  
  # make_emc with compression (default)
  emc <- make_emc(dat, design_LBA, type = "single", n_chains = 1)
  
  # n_pars for LBA is 5
  n_pars <- 5
  emc[[1]]$samples$alpha <- array(rnorm(n_pars * 1 * 10), dim = c(n_pars, 1, 10))
  dimnames(emc[[1]]$samples$alpha) <- list(c("v", "sv", "B", "A", "t0"), "as1t", NULL)
  emc[[1]]$samples$stage <- rep("sample", 10)
  emc[[1]]$samples$idx <- 10
  
  # 1. Test add_pw_ll
  emc <- add_pw_ll(emc)
  expect_true(!is.null(emc[[1]]$samples$pw_ll))
  # dat has 810 trials. LBA DADM has 1620 rows. 
  # pw_ll should have 810 rows (one per trial).
  expect_equal(nrow(emc[[1]]$samples$pw_ll), nrow(dat))
  
  # 2. Test extract_pw_ll
  pw_ll_df <- extract_pw_ll(emc, stage = "sample")
  # Should have one row per trial
  expect_equal(nrow(pw_ll_df), nrow(dat))
  expect_true("ll" %in% names(pw_ll_df))
})

test_that("pw_ll extraction with DDM (already trial-wise)", {
  data(forstmann)
  dat <- forstmann[forstmann$subjects == unique(forstmann$subjects)[1],]
  dat$subjects <- droplevels(dat$subjects)
  
  design_DDM <- design(data = dat, model=DDM, 
                      formula=list(v~1, a~1, sv~1, t0~1, Z~1, SZ~1, st0~1))
  
  emc <- make_emc(dat, design_DDM, type = "single", n_chains = 1)
  
  n_pars <- 7
  emc[[1]]$samples$alpha <- array(rnorm(n_pars * 1 * 10), dim = c(n_pars, 1, 10))
  dimnames(emc[[1]]$samples$alpha) <- list(c("v", "a", "sv", "t0", "st0", "Z", "SZ"), "as1t", NULL)
  emc[[1]]$samples$stage <- rep("sample", 10)
  emc[[1]]$samples$idx <- 10
  
  emc <- add_pw_ll(emc)
  pw_ll_df <- extract_pw_ll(emc, stage = "sample")
  
  expect_equal(nrow(pw_ll_df), nrow(dat))
})
