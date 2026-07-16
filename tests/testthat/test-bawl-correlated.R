test_that("BAwLcorr supports correlated racers and independent rho=0 racers", {
  skip_on_cran()
  n <- 2500
  lR <- factor(rep(c("correct", "error", "pm"), n),
               levels = c("correct", "error", "pm"))
  pars <- cbind(
    v = rep(c(1.0, 0.4, 0.2), n), sv = 1, b = 1.5, A = 0.3,
    t0 = 0.2, k = 0.2, lambda_g = 0, lambda_k = 0,
    rho = rep(c(0.6, -0.6, 0), n)
  )
  ok <- rep(TRUE, nrow(pars))

  set.seed(101)
  cpp <- EMC2:::rbawl_corr_cpp(pars, levels(lR), ok, TRUE, 1L, FALSE, FALSE)
  set.seed(102)
  ref <- EMC2:::rBAwL_corr(lR, pars, ok = ok, posdrift = TRUE,
                           erlang = 1L, guess = FALSE, global = FALSE)

  expect_length(cpp$R, n)
  expect_equal(nrow(ref), n)
  expect_lt(abs(mean(cpp$R == 1, na.rm = TRUE) -
                mean(ref$R == "correct", na.rm = TRUE)), 0.05)
  expect_lt(abs(mean(cpp$R == 3, na.rm = TRUE) -
                mean(ref$R == "pm", na.rm = TRUE)), 0.05)
  expect_equal(BAwLcorr()$c_name, "BAwL_CORR")
  expect_true("rho" %in% names(BAwLcorr()$p_types))
})

test_that("BAwLcorr likelihood uses the shared race path and supports truncation", {
  skip_on_cran()
  des <- design(
    factors = list(subjects = 1, S = "stim"),
    Rlevels = c("correct", "error", "pm"),
    formula = list(v ~ 0 + lR, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1,
                   k ~ 1, mG ~ 1, mK ~ 1, rho ~ 0 + lR),
    constants = c(rho_lRpm = 0),
    model = BAwLcorr(), report_p_vector = FALSE
  )
  dat <- data.frame(
    subjects = factor(1), S = factor("stim"),
    R = factor("correct", levels = c("correct", "error", "pm")),
    rt = 0.7, LT = 0.2, UT = 1.5
  )
  emc <- make_emc(dat, des, type = "single", n_chains = 1,
                  compress = FALSE, verbose = FALSE)
  model <- emc[[1]]$model()
  dadm <- emc[[1]]$data[[1]]
  p <- sampled_pars(des, doMap = FALSE)
  p[grep("^v_lR", names(p))] <- c(1.0, 0.4, 0.2)
  p[c("sv", "B", "A", "t0", "k", "mG", "mK")] <-
    log(c(1, 1.5, 0.3, 0.2, 0.2, 1, 1))
  p[grep("^rho_lR", names(p))] <- c(0.6, -0.6)

  designs <- lapply(names(model$p_types), function(nm) {
    dm <- attr(dadm, "designs")[[nm]]
    dm[attr(dm, "expand"), , drop = FALSE]
  })
  names(designs) <- names(model$p_types)
  args <- list(
    particle_matrix = matrix(p, nrow = 1, dimnames = list(NULL, names(p))),
    data = dadm, constants = attr(dadm, "constants"), designs = designs,
    type = model$c_name, bounds = model$bound, transforms = model$transform,
    pretransforms = model$pre_transform, p_types = names(model$p_types),
    min_ll = log(1e-10), trend = model$trend
  )
  ll <- do.call(EMC2:::calc_ll_oo, args)
  ll_pw <- do.call(EMC2:::calc_ll_oo_pw, args)

  expect_true(is.finite(ll))
  expect_equal(as.numeric(ll), as.numeric(ll_pw), tolerance = 1e-8)
})
