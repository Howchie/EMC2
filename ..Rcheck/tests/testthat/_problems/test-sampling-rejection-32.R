# Extracted from test-sampling-rejection.R:32

# test -------------------------------------------------------------------------
pm_settings <- list(list(epsilon = 1))
parameters <- list(alpha = matrix(c(1, 2), ncol = 1))
out <- safe_new_particle(
    s = 1,
    data = NULL,
    pm_settings = pm_settings,
    prev_ll = -7,
    parameters = parameters,
    model = NULL,
    stage = "preburn",
    type = "standard",
    tune = list(components = c(1, 1), shared_ll_idx = c(1, 1)),
    r_cores = 1
  )
