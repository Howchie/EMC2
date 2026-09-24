local_rng_guard()  # see helper-rng.R: keep this file's RNG changes inside it
test_that("plot_pars", {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
  vdiffr::expect_doppelganger("density_plots", plot_pars(samples_LNR, selection = "sigma2"))
  vdiffr::expect_doppelganger("density_alpha", plot_pars(samples_LNR, all_subjects = T))
})
