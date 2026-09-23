library(EMC2)

# design(TC = list(filter_defective = TRUE)): at a finite UT the never-finish
# atom of a defective race is truncated along with the finite tail, so the
# normaliser is P(LT <= T <= UT).  The default keeps the atom in the retained
# sample space, P(LT <= T <= UT) + P(T = Inf).  Every check below compares the
# engine against the finite-window mass integrated from UNTRUNCATED densities,
# and against P(T = Inf) read from an untruncated omission trial.

fd_context <- function(dat, model, formula, constants = NULL, fd = FALSE) {
  des <- design(data = dat, Rlevels = levels(dat$R),
                matchfun = function(d) as.character(d$S) == as.character(d$lR),
                formula = formula, constants = constants, model = model,
                report_p_vector = FALSE,
                TC = if (fd) list(filter_defective = TRUE) else NULL)
  emc <- make_emc(dat, des, type = "single", n_chains = 1, compress = FALSE,
                  verbose = FALSE, rt_resolution = NULL)
  list(emc = emc[[1]], design = des)
}

fd_trial_ll <- function(ctx, p) {
  dadm <- ctx$emc$data[[1]]
  model <- ctx$emc$model()
  designs <- lapply(names(model$p_types), function(nm) {
    dm <- attr(dadm, "designs")[[nm]]
    dm[attr(dm, "expand"), , drop = FALSE]
  })
  names(designs) <- names(model$p_types)
  as.numeric(EMC2:::calc_ll_oo_pw(
    particle_matrix = matrix(p, nrow = 1, dimnames = list(NULL, names(p))),
    data = dadm, constants = attr(dadm, "constants"), designs = designs,
    type = model$c_name, bounds = model$bound, transforms = model$transform,
    pretransforms = model$pre_transform, p_types = names(model$p_types),
    min_ll = log(1e-300), trend = model$trend))
}

# Finite response, the other response, a pure omission and an upper-censored
# trial with unknown response (UC < UT).
fd_data <- function(LT, UT, UC, Rlev) {
  data.frame(subjects = factor(1), S = factor(Rlev[1], levels = Rlev),
             R = factor(c(Rlev[1], Rlev[2], NA, NA), levels = Rlev),
             rt = c(.8, 1.1, Inf, Inf), LT = LT, UT = UT, LC = 0,
             UC = c(Inf, Inf, Inf, UC))
}

fd_window_mass <- function(model, formula, constants, p, a, b, Rlev,
                           n = 2000) {
  x <- seq(a, b, length.out = n + 1)
  d <- data.frame(subjects = factor(1), S = factor(Rlev[1], levels = Rlev),
                  R = factor(rep(Rlev, each = length(x)), levels = Rlev),
                  rt = rep(x, length(Rlev)))
  y <- rowSums(matrix(exp(fd_trial_ll(fd_context(d, model, formula, constants), p)),
                      ncol = length(Rlev)))
  h <- (b - a) / n
  h / 3 * (y[1] + y[n + 1] + 4 * sum(y[seq(2, n, 2)]) + 2 * sum(y[seq(3, n - 1, 2)]))
}

expect_atom_semantics <- function(model, formula, constants, fill,
                                  Rlev = c("left", "right"), tol = 1e-6) {
  LT <- .3; UT <- 2.5; UC <- 1.2
  truncated <- fd_data(LT, UT, UC, Rlev)
  keep_ctx <- fd_context(truncated, model, formula, constants)
  p <- fill(sampled_pars(keep_ctx$design, doMap = FALSE))
  plain <- fd_trial_ll(fd_context(fd_data(0, Inf, UC, Rlev), model, formula,
                                  constants), p)
  S_inf <- exp(plain[3])
  expect_gt(S_inf, 1e-3)                     # the atom is material
  Z_drop <- fd_window_mass(model, formula, constants, p, LT, UT, Rlev)
  W_cens <- fd_window_mass(model, formula, constants, p, UC, UT, Rlev)
  Z_keep <- Z_drop + S_inf

  keep <- fd_trial_ll(keep_ctx, p)
  expect_equal(keep, c(plain[1:2] - log(Z_keep), log(S_inf / Z_keep),
                       log((W_cens + S_inf) / Z_keep)), tolerance = tol)

  drop <- fd_trial_ll(fd_context(truncated[-3, ], model, formula, constants,
                                 fd = TRUE), p)
  expect_equal(drop, c(plain[1:2] - log(Z_drop), log(W_cens / Z_drop)),
               tolerance = tol)
}

test_that("ordinary race normalises with and without the defective atom", {
  expect_atom_semantics(
    LBA(posdrift = FALSE), list(v ~ lM, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1), NULL,
    function(p) {
      p[] <- 0; p["v"] <- -.4; p["v_lMTRUE"] <- 1.6; p["sv"] <- log(.6)
      p["B"] <- log(.8); p["A"] <- log(.4); p["t0"] <- log(.2); p
    })
  expect_atom_semantics(
    RDMSWTN(posdrift = FALSE), list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1),
    c(sv = log(0)),
    function(p) {
      p[] <- 0; p["v"] <- -.5; p["v_lMTRUE"] <- .45; p["A"] <- log(.3)
      p["t0"] <- log(.2); p
    })
})

bawl_corr_fill <- function(p) {
  p[] <- 0
  p["v"] <- -.4; p[grep("^v_", names(p))] <- 1.6
  p["sv"] <- log(.6); p["B"] <- log(.8); p["A"] <- log(.4); p["t0"] <- log(.2)
  p["k"] <- log(0); p["rho"] <- qnorm(.75)          # rho = .5
  p
}
bawl_corr_formula <- list(v ~ lM, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1,
                          mG ~ 1, mK ~ 1, rho ~ 1)
rdm_corr_fill <- function(p) {
  p[] <- 0; p["v"] <- -.5; p[grep("^v_", names(p))] <- .45
  p["A"] <- log(.3); p["t0"] <- log(.2)
  if ("sv" %in% names(p)) p["sv"] <- log(.5)
  p["rho"] <- qnorm(.75)
  p
}

test_that("correlated routes honour the same retained sample space", {
  skip_model_validation()
  # Exact BVN pair and the finishing-time copula are deterministic.
  expect_atom_semantics(BAwLcorr(posdrift = FALSE), bawl_corr_formula, NULL,
                        bawl_corr_fill)
  expect_atom_semantics(
    RDMSWTNcorr(posdrift = FALSE, correlate = "times"),
    list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, rho ~ 1), c(sv = log(0)),
    rdm_corr_fill)
  # Three loaded racers take the Gauss-Hermite factor route, whose default
  # node count carries ~1e-5 quadrature error of its own.
  expect_atom_semantics(BAwLcorr(posdrift = FALSE), bawl_corr_formula, NULL,
                        bawl_corr_fill, Rlev = c("a", "b", "c"), tol = 1e-4)
})

test_that("filter_defective validates the data it is declared for", {
  dat <- data.frame(subjects = factor(1), S = factor("left", levels = c("left", "right")),
                    R = factor(c("left", NA), levels = c("left", "right")),
                    rt = c(.8, Inf), UT = 2.5, UC = Inf)
  f <- list(v ~ lM, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1)
  mk <- function(formula, TC) {
    des <- design(data = dat, Rlevels = levels(dat$R),
                  matchfun = function(d) as.character(d$S) == as.character(d$lR),
                  formula = formula, model = LBA(posdrift = FALSE),
                  report_p_vector = FALSE, TC = TC)
    make_emc(dat, des, type = "single", n_chains = 1, compress = FALSE,
             verbose = FALSE, rt_resolution = NULL)
  }
  expect_error(mk(f, list(filter_defective = NA)), "single TRUE or FALSE")
  expect_error(mk(f, list(filter_defective = TRUE)), "omissions")
  # Retained model omissions are the default, and contaminant omissions are
  # applied after truncation, so either explains the omitted trial.
  emc <- mk(f, NULL)
  expect_null(attr(emc[[1]]$data[[1]], "emc2_filter_defective"))
  emc <- mk(c(f, pContaminant ~ 1), list(filter_defective = TRUE))
  expect_true(attr(emc[[1]]$data[[1]], "emc2_filter_defective"))
})
