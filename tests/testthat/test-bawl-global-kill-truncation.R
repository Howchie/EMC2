skip_model_validation()

# The global kill applies its shared clock as a per-trial factor at the
# likelihood-assembly sites rather than inside cdf1 (apply_lk_to_racers is
# false), so the truncation normaliser has to switch the kill back on for the
# racer.  Before this was fixed, Z was P(no racer crossed by LT) instead of
# P(no response by LT) = 1 - int_0^LT f_race S_K, and global_kill silently
# disagreed with the equivalent local_kill fit.

bawl_kill_ctx <- function(erlang_type, dat, Rlevels, v_formula) {
  des <- design(
    factors = list(subjects = 1, S = "stim"),
    Rlevels = Rlevels,
    formula = c(v_formula, list(sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1,
                                mG ~ 1, mK ~ 1)),
    model = function() BAwL(erlang_type = erlang_type),
    report_p_vector = FALSE
  )
  emc <- make_emc(dat, des, type = "single", compress = FALSE, n_chains = 1,
                  rt_resolution = NULL)
  list(design = des, emc = emc)
}

bawl_kill_ll <- function(ctx, p) {
  model <- ctx$emc[[1]]$model()
  dadm <- ctx$emc[[1]]$data[[1]]
  p_mat <- matrix(p, nrow = 1, dimnames = list(NULL, names(p)))
  dz <- attr(dadm, "designs")
  designs <- lapply(names(model$p_types), function(q)
    dz[[q]][attr(dz[[q]], "expand"), , drop = FALSE])
  names(designs) <- names(model$p_types)
  constants <- attr(dadm, "constants")
  if (is.null(constants)) constants <- NA
  as.numeric(EMC2:::calc_ll_oo(
    p_mat, dadm, constants = constants, designs = designs,
    type = model$c_name, bounds = model$bound, transforms = model$transform,
    pretransforms = model$pre_transform, p_types = names(model$p_types),
    min_ll = log(1e-300), trend = model$trend))
}

# Killed sub-density / sub-CDF for one accumulator, from a mapped parameter row.
bawl_kill_d <- function(t, r) EMC2:::dkilledleakyba(
  t, r["v"], r["b"], r["A"], r["sv"], r["t0"], r["k"], 0, r["lambda_k"],
  TRUE, TRUE, 1L, FALSE, 1, 0L)
bawl_kill_p <- function(t, r) EMC2:::pkilledleakyba(
  t, r["v"], r["b"], r["A"], r["sv"], r["t0"], r["k"], 0, r["lambda_k"],
  TRUE, FALSE, 1L, FALSE, 1, 0L)

bawl_kill_pars <- function(ctx, p)
  EMC2:::get_pars_matrix_oo(p, ctx$emc[[1]]$data[[1]], ctx$emc[[1]]$model())

test_that("global_kill keeps the kill in the lower-truncation normaliser", {
  skip_on_cran()
  LT <- 0.45
  dat <- data.frame(subjects = factor(1), S = factor("stim"),
                    R = factor("left", levels = "left"), rt = 0.75, LT = LT)
  loc <- bawl_kill_ctx("local_kill", dat, "left", list(v ~ 1))
  glo <- bawl_kill_ctx("global_kill", dat, "left", list(v ~ 1))
  p <- sampled_pars(loc$design, doMap = FALSE)
  p[c("v", "sv", "B", "A", "t0", "k", "mG", "mK")] <-
    c(1.2, log(.6), log(.9), log(.3), log(.08), log(.5), log(1), log(.35))

  r <- bawl_kill_pars(loc, p)[1, ]
  # Z = P(no response before LT) = 1 - int_0^LT f_R S_K, not S_R(LT) S_K(LT).
  reference <- bawl_kill_d(dat$rt, r) - log(1 - bawl_kill_p(LT, r))

  expect_equal(bawl_kill_ll(loc, p), reference, tolerance = 1e-10)
  expect_equal(bawl_kill_ll(glo, p), reference, tolerance = 1e-10)
})

test_that("local_kill factorises correctly over independent accumulators", {
  skip_on_cran()
  LT <- 0.45
  dat <- data.frame(subjects = factor(1), S = factor("stim"),
                    R = factor("left", levels = c("left", "right")),
                    rt = 0.75, LT = LT)
  ctx <- bawl_kill_ctx("local_kill", dat, c("left", "right"), list(v ~ 0 + lR))
  p <- sampled_pars(ctx$design, doMap = FALSE)
  p[grep("^v_", names(p))] <- c(1.2, .7)
  p[c("sv", "B", "A", "t0", "k", "mG", "mK")] <-
    c(log(.6), log(.9), log(.3), log(.08), log(.5), log(1), log(.35))
  pars <- bawl_kill_pars(ctx, p)

  # Independent per-accumulator clocks DO factorise into per-racer killed
  # sub-CDFs, so the product form is exact here.
  numerator <- bawl_kill_d(dat$rt, pars[1, ]) +
    log1p(-bawl_kill_p(dat$rt, pars[2, ]))
  log_Z <- sum(log1p(-c(bawl_kill_p(LT, pars[1, ]),
                        bawl_kill_p(LT, pars[2, ]))))
  expect_equal(bawl_kill_ll(ctx, p), numerator - log_Z, tolerance = 1e-10)
})

test_that("global_kill refuses truncation it cannot normalise", {
  skip_on_cran()
  dat <- data.frame(subjects = factor(1), S = factor("stim"),
                    R = factor("left", levels = c("left", "right")),
                    rt = 0.75, LT = 0.45)
  ctx <- bawl_kill_ctx("global_kill", dat, c("left", "right"), list(v ~ 0 + lR))
  p <- sampled_pars(ctx$design, doMap = FALSE)
  p[grep("^v_", names(p))] <- c(1.2, .7)
  p[c("sv", "B", "A", "t0", "k", "mG", "mK")] <-
    c(log(.6), log(.9), log(.3), log(.08), log(.5), log(1), log(.35))

  # A shared clock does not factorise into per-accumulator survivors, so there
  # is no product form to correct: refuse rather than return a wrong Z.
  expect_error(bawl_kill_ll(ctx, p), "does not support truncation")

  # Without truncation the same model is fine.
  dat0 <- dat
  dat0$LT <- 0
  ctx0 <- bawl_kill_ctx("global_kill", dat0, c("left", "right"),
                        list(v ~ 0 + lR))
  expect_true(is.finite(bawl_kill_ll(ctx0, p)))
})
