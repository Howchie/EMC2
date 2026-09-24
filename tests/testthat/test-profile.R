local_rng_guard()  # see helper-rng.R: keep this file's RNG changes inside it
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
test_that("profile_plot", {
  vdiffr::expect_doppelganger("profile_plot",
    profile_plot(get_data(samples_LNR), get_design(samples_LNR)[[1]], credint(samples_LNR, probs = .5)[[1]]))
})
