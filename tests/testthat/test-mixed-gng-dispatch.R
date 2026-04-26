build_ll_ctx <- function(data, design) {
  emc <- make_emc(data, design, type = "single", compress = FALSE, n_chains = 1)
  model <- emc[[1]]$model()
  dadm <- emc[[1]]$data[[1]]
  p_types <- names(model$p_types)
  designs <- setNames(vector("list", length(p_types)), p_types)
  for (p in p_types) {
    dm <- attr(dadm, "designs")[[p]]
    designs[[p]] <- dm[attr(dm, "expand"), , drop = FALSE]
  }
  constants <- attr(dadm, "constants")
  if (is.null(constants)) constants <- NA
  list(dadm = dadm, model = model, p_types = p_types, designs = designs, constants = constants)
}

calc_ctx_ll <- function(ctx, p_mat) {
  EMC2:::calc_ll_oo(
    p_mat,
    ctx$dadm,
    constants = ctx$constants,
    designs = ctx$designs,
    type = ctx$model$c_name,
    bounds = ctx$model$bound,
    transforms = ctx$model$transform,
    pretransforms = ctx$model$pre_transform,
    p_types = ctx$p_types,
    min_ll = log(1e-10),
    trend = ctx$model$trend
  )
}

test_that("mixed 2AFC and GNG pooled race data dispatches by active nogo accumulator", {
  pooled <- data.frame(
    subjects = factor(rep("s1", 6)),
    Condition = factor(c("2afc", "2afc", "2afc", "gng", "gng", "gng"),
                       levels = c("2afc", "gng")),
    RACE = factor(c("2", "2", "2", "3", "3", "3"), levels = c("2", "3")),
    R = factor(c("left", NA, NA, "right", NA, NA), levels = c("left", "right", "nogo")),
    rt = c(0.35, Inf, -Inf, 0.42, Inf, -Inf),
    LC = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2),
    UC = c(0.8, 0.8, 0.8, 0.8, 0.8, 0.8)
  )

  pooled_design <- design(
    data = pooled,
    model = LBA,
    formula = list(v ~ 1, B ~ 1, t0 ~ 1, A ~ 1)
  )
  pooled_ctx <- build_ll_ctx(pooled, pooled_design)

  afc_data <- droplevels(subset(pooled, Condition == "2afc", select = c(subjects, R, rt, LC, UC)))
  afc_data$R <- factor(afc_data$R, levels = c("left", "right"))
  afc_design <- design(
    data = afc_data,
    model = LBA,
    formula = list(v ~ 1, B ~ 1, t0 ~ 1, A ~ 1)
  )
  afc_ctx <- build_ll_ctx(afc_data, afc_design)

  gng_data <- droplevels(subset(pooled, Condition == "gng", select = c(subjects, R, rt, LC, UC)))
  gng_data$R <- factor(gng_data$R, levels = c("left", "right", "nogo"))
  gng_design <- design(
    data = gng_data,
    model = LBA,
    formula = list(v ~ 1, B ~ 1, t0 ~ 1, A ~ 1)
  )
  gng_ctx <- build_ll_ctx(gng_data, gng_design)

  expect_identical(pooled_ctx$p_types, afc_ctx$p_types)
  expect_identical(pooled_ctx$p_types, gng_ctx$p_types)

  set.seed(42)
  p_mat <- matrix(rnorm(16 * length(pooled_ctx$p_types), sd = 0.2), nrow = 16)
  colnames(p_mat) <- pooled_ctx$p_types

  pooled_ll <- calc_ctx_ll(pooled_ctx, p_mat)
  split_ll <- calc_ctx_ll(afc_ctx, p_mat) + calc_ctx_ll(gng_ctx, p_mat)

  expect_true(all(is.finite(pooled_ll)))
  expect_true(all(is.finite(split_ll)))
  expect_equal(pooled_ll, split_ll, tolerance = 1e-8)
})

test_that("make_missing only disables truncation on rows with active nogo", {
  dat <- data.frame(
    subjects = factor(c("s1", "s1")),
    RACE = factor(c("2", "3"), levels = c("2", "3")),
    R = factor(c(NA, NA), levels = c("left", "right", "nogo")),
    rt = c(0.5, 0.5)
  )

  out <- make_missing(dat, LT = 0.2, UT = 0.9, rt_resolution = NULL)

  expect_equal(out$LT[1], 0.2)
  expect_equal(out$UT[1], 0.9)
  expect_equal(out$LT[2], 0)
  expect_true(is.infinite(out$UT[2]))
})

test_that("make_data supports pooled mixed 2AFC and GNG race templates", {
  template <- data.frame(
    subjects = factor(rep("s1", 4)),
    Condition = factor(c("2afc", "2afc", "gng", "gng"), levels = c("2afc", "gng")),
    S = factor(c("left", "right", "left", "stop"), levels = c("left", "right", "stop")),
    RACE = factor(c("2", "2", "3", "3"), levels = c("2", "3")),
    R = factor(rep(NA_character_, 4), levels = c("left", "right", "nogo"))
  )

  matchfun <- function(d) {
    (d$Condition == "2afc" & as.character(d$S) == as.character(d$lR)) |
      (d$Condition == "gng" &
         ((d$S == "left" & d$lR == "left") | (d$S == "stop" & d$lR == "nogo")))
  }

  design_mixed <- design(
    data = template,
    Rlevels = c("left", "right", "nogo"),
    matchfun = matchfun,
    model = LBA,
    formula = list(v ~ 0 + lM, B ~ 0 + lR, t0 ~ 1, A ~ 1),
    constants = c(sv = log(1))
  )

  p_vector <- sampled_pars(design_mixed, doMap = FALSE)
  if ("v_lMFALSE" %in% names(p_vector)) p_vector[["v_lMFALSE"]] <- 0.3
  if ("v_lMTRUE" %in% names(p_vector)) p_vector[["v_lMTRUE"]] <- 1.6
  if ("B_lRleft" %in% names(p_vector)) p_vector[["B_lRleft"]] <- log(0.8)
  if ("B_lRright" %in% names(p_vector)) p_vector[["B_lRright"]] <- log(0.8)
  if ("B_lRnogo" %in% names(p_vector)) p_vector[["B_lRnogo"]] <- log(0.6)
  if ("t0" %in% names(p_vector)) p_vector[["t0"]] <- log(0.2)
  if ("A" %in% names(p_vector)) p_vector[["A"]] <- log(0.3)

  set.seed(7)
  dat <- make_data(
    p_vector,
    design_mixed,
    data = template,
    expand = 40,
    TC = list(LT = 0.2, UT = 0.9, UC = 0.8)
  )

  expect_true("RACE" %in% names(dat))
  expect_true(all(dat$LT[dat$Condition == "gng"] == 0))
  expect_true(all(is.infinite(dat$UT[dat$Condition == "gng"])))
  expect_true(all(dat$LT[dat$Condition == "2afc"] == 0.2))
  expect_true(all(dat$UT[dat$Condition == "2afc"] == 0.9))
})
