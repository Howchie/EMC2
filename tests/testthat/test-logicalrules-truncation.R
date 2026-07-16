# Validation of the logical-rules truncation normaliser:
#   Z = P(overt rule response in [LT, UT]) = N(LT) - N(UT)
# with N(t) the rule-correct no-response probability (UT = Inf: Z = N(LT)).
#
# Strategy: (1) self-normalization — the truncated outcome space must sum to
# one, integrating per-trial densities from calc_ll_oo_pw over the window and
# adding censor masses (this pits the normaliser's channel-state math against
# the independently-computed GL-pass densities); (2) direct MC check — the
# log-likelihood shift from adding LT/UT equals -log P(RT in window) estimated
# by simulation.

build_ll_ctx <- function(data, design) {
  # rt_resolution = NULL: the integration grids below must not be floored to
  # the default 1/60 s resolution (500 grid points would collapse to ~48).
  emc <- make_emc(data, design, type = "single", compress = FALSE, n_chains = 1,
                  rt_resolution = NULL)
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

ctx_args <- function(ctx, p_vec) {
  list(
    matrix(p_vec, nrow = 1, dimnames = list(NULL, names(p_vec))),
    ctx$dadm, constants = ctx$constants, designs = ctx$designs,
    type = ctx$model$c_name, bounds = ctx$model$bound,
    transforms = ctx$model$transform, pretransforms = ctx$model$pre_transform,
    p_types = ctx$p_types, min_ll = log(1e-10), trend = ctx$model$trend
  )
}

calc_ctx_ll <- function(ctx, p_vec) do.call(EMC2:::calc_ll_oo, ctx_args(ctx, p_vec))
calc_ctx_ll_pw <- function(ctx, p_vec) drop(do.call(EMC2:::calc_ll_oo_pw, ctx_args(ctx, p_vec)))

lr_funcs <- list(
  GoA = function(d) ifelse(d$lR == "A", 1, 0),
  GoB = function(d) ifelse(d$lR == "B", 1, 0),
  NegA = function(d) ifelse(d$lR == "n_A", 1, 0),
  NegB = function(d) ifelse(d$lR == "n_B", 1, 0),
  NoGo = function(d) ifelse(d$lR == "nogo", 1, 0)
)

lr_pvec <- function(design) {
  p <- sampled_pars(design, doMap = FALSE)
  p[] <- 0
  vals <- c(v_GoA = 1.60, v_GoB = 1.10, v_NegA = 0.95, v_NegB = 1.25,
            B_GoA = log(0.75), B_GoB = log(0.95), B_NegA = log(0.85), B_NegB = log(0.70),
            v_NoGo = 1.30, B_NoGo = log(0.80),
            t0 = log(0.2), A = log(0.3))
  for (nm in names(vals)) if (nm %in% names(p)) p[[nm]] <- vals[[nm]]
  p
}

make_rule_design <- function(template, Rlevels) {
  # R contains overt rule responses in the simulation data, but the design
  # scaffold must expose the latent accumulator roles used by the formulas.
  design_template <- template
  design_template$R <- factor(NA_character_, levels = c("A", "B", "n_A", "n_B"))
  design(
    data = design_template,
    Rlevels = Rlevels,
    fixed_accumulator_roles = factor(c("A", "B", "n_A", "n_B"),
                                     levels = c("A", "B", "n_A", "n_B")),
    matchfun = function(d) d$lR %in% c("A", "B", "n_A", "n_B"),
    model = LogicalRulesLBA,
    formula = list(v ~ 0 + GoA + GoB + NegA + NegB,
                   B ~ 0 + GoA + GoB + NegA + NegB,
                   t0 ~ 1, A ~ 1),
    constants = c(sv = log(1)),
    functions = lr_funcs
  )
}

lr_rows <- function(rule, Rlevels, R, rt, ...) {
  extras <- list(...)
  d <- data.frame(
    subjects = factor(rep("s1", length(rt))),
    S = factor(rep("AB", length(rt)), levels = "AB"),
    LogicalRule = factor(rep(rule, length(rt)), levels = rule),
    R = factor(R, levels = Rlevels),
    rt = rt
  )
  for (nm in names(extras)) d[[nm]] <- extras[[nm]]
  d
}

# Trapezoid-integrate exp(ll) over the trials of ctx$dadm, grouped by response.
integrate_pw <- function(ctx, p_vec, n_acc = 4) {
  lls <- calc_ctx_ll_pw(ctx, p_vec)
  first_rows <- seq(1, nrow(ctx$dadm), by = n_acc)
  rt <- ctx$dadm$rt[first_rows]
  R <- ctx$dadm$R[first_rows]
  total <- 0
  for (r in levels(droplevels(R))) {
    pick <- which(R == r)
    o <- pick[order(rt[pick])]
    x <- rt[o]
    y <- exp(lls[o])
    total <- total + sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  }
  total
}

test_that("OR and ID truncated outcome space integrates to one", {
  skip_on_cran()
  LT <- 0.3
  UT <- 1.1
  n_grid <- 500

  for (spec in list(list(rule = "OR", Rlevels = c("yes", "no")),
                    list(rule = "ID", Rlevels = c("NN", "AN", "NB", "AB")))) {
    template <- lr_rows(spec$rule, spec$Rlevels, R = NA_character_, rt = NA_real_)
    template$rt <- NULL
    des <- make_rule_design(template, spec$Rlevels)
    p <- lr_pvec(des)

    grid <- seq(LT + 1e-4, UT - 1e-4, length.out = n_grid)
    dat <- do.call(rbind, lapply(spec$Rlevels, function(r) {
      lr_rows(spec$rule, spec$Rlevels, R = rep(r, n_grid), rt = grid, LT = LT, UT = UT)
    }))
    ctx <- build_ll_ctx(dat, des)
    total <- integrate_pw(ctx, p)
    expect_equal(total, 1.0, tolerance = 5e-3,
                 label = sprintf("%s truncated total mass (%g)", spec$rule, total))
  }
})

test_that("truncation shift equals -log P(window) from simulation (OR)", {
  skip_on_cran()
  n_sim <- 1e5
  LT <- 0.3
  UT <- 1.1
  Rlevels <- c("yes", "no")

  template <- lr_rows("OR", Rlevels, R = NA_character_, rt = NA_real_)
  template$rt <- NULL
  des <- make_rule_design(template, Rlevels)
  p <- lr_pvec(des)

  set.seed(104)
  sim <- make_data(p, des, data = template, expand = n_sim)
  p_window <- mean(sim$rt >= LT & sim$rt <= UT, na.rm = TRUE)

  rt0 <- 0.6
  ctx_plain <- build_ll_ctx(lr_rows("OR", Rlevels, "yes", rt0), des)
  ctx_trunc <- build_ll_ctx(lr_rows("OR", Rlevels, "yes", rt0, LT = LT, UT = UT), des)
  logZ <- calc_ctx_ll(ctx_plain, p) - calc_ctx_ll(ctx_trunc, p)
  expect_lt(abs(exp(logZ) - p_window),
            4 * sqrt(p_window * (1 - p_window) / n_sim) + 2e-4,
            label = sprintf("OR Z (%g) vs MC window mass (%g)", exp(logZ), p_window))

  # LT-only truncation (UT = Inf): Z = N(LT)
  p_lt <- mean(sim$rt >= LT, na.rm = TRUE)
  ctx_lt <- build_ll_ctx(lr_rows("OR", Rlevels, "yes", rt0, LT = LT), des)
  logZ_lt <- calc_ctx_ll(ctx_plain, p) - calc_ctx_ll(ctx_lt, p)
  expect_lt(abs(exp(logZ_lt) - p_lt),
            4 * sqrt(p_lt * (1 - p_lt) / n_sim) + 2e-4,
            label = sprintf("OR Z_LT (%g) vs MC (%g)", exp(logZ_lt), p_lt))
})

test_that("truncation + censoring outcome space integrates to one (OR)", {
  skip_on_cran()
  LT <- 0.25
  LC <- 0.55
  UC <- 0.95
  UT <- 1.5
  Rlevels <- c("yes", "no")
  n_grid <- 400

  template <- lr_rows("OR", Rlevels, R = NA_character_, rt = NA_real_)
  template$rt <- NULL
  des <- make_rule_design(template, Rlevels)
  p <- lr_pvec(des)

  # Observable finite window [LC, UC]
  grid <- seq(LC + 1e-4, UC - 1e-4, length.out = n_grid)
  dat_fin <- do.call(rbind, lapply(Rlevels, function(r) {
    lr_rows("OR", Rlevels, R = rep(r, n_grid), rt = grid, LT = LT, UT = UT)
  }))
  finite_mass <- integrate_pw(build_ll_ctx(dat_fin, des), p)

  # Lower-censored masses per response
  low_mass <- sum(vapply(Rlevels, function(r) {
    exp(calc_ctx_ll(build_ll_ctx(
      lr_rows("OR", Rlevels, r, -Inf, LT = LT, UT = UT, LC = LC), des), p))
  }, numeric(1)))

  # Upper-censored mass (response in (UC, UT]; resp-agnostic, recorded R required)
  up_mass <- exp(calc_ctx_ll(build_ll_ctx(
    lr_rows("OR", Rlevels, "yes", Inf, LT = LT, UT = UT, UC = UC), des), p))

  total <- finite_mass + low_mass + up_mass
  expect_equal(total, 1.0, tolerance = 5e-3,
               label = sprintf("OR censored+truncated total (%g = %g + %g + %g)",
                               total, finite_mass, low_mass, up_mass))
})

test_that("OR_DETECTION_GNG truncated outcome space integrates to one", {
  skip_on_cran()
  LT <- 0.25
  LC <- 0.45
  UC <- 0.9
  UT <- 1.4
  stim_levels <- c("NN", "AN", "NB", "AB")
  n_grid <- 400

  template <- data.frame(
    subjects = factor(rep("s1", 4)),
    S = factor(stim_levels, levels = stim_levels),
    LogicalRule = factor(rep("OR_DETECTION_GNG", 4), levels = "OR_DETECTION_GNG"),
    R = factor(rep(NA_character_, 4), levels = c("yes", "no"))
  )
  design_template <- template
  design_template$R <- factor(NA_character_, levels = c("A", "B", "nogo"))
  des <- design(
    data = design_template,
    Rlevels = c("yes", "no"),
    fixed_accumulator_roles = factor(c("A", "B", "nogo"), levels = c("A", "B", "nogo")),
    matchfun = function(d) d$lR %in% c("A", "B", "nogo"),
    model = LogicalRulesLBA,
    formula = list(v ~ 0 + GoA + GoB + NoGo,
                   B ~ 0 + GoA + GoB + NoGo,
                   t0 ~ 1, A ~ 1),
    constants = c(sv = log(1)),
    functions = lr_funcs
  )
  p <- lr_pvec(des)

  gng_rows <- function(S, R, rt, ...) {
    extras <- list(...)
    d <- data.frame(
      subjects = factor(rep("s1", length(rt))),
      S = factor(rep(S, length(rt)), levels = stim_levels),
      LogicalRule = factor(rep("OR_DETECTION_GNG", length(rt)),
                           levels = "OR_DETECTION_GNG"),
      R = factor(R, levels = c("yes", "no")),
      rt = rt
    )
    for (nm in names(extras)) d[[nm]] <- extras[[nm]]
    d
  }

  for (S in c("AN", "AB")) {
    grid <- seq(LC + 1e-4, UC - 1e-4, length.out = n_grid)
    finite_mass <- integrate_pw(build_ll_ctx(
      gng_rows(S, rep("yes", n_grid), grid, LT = LT, UT = UT), des), p, n_acc = 3)
    low_mass <- exp(calc_ctx_ll(build_ll_ctx(
      gng_rows(S, "yes", -Inf, LT = LT, UT = UT, LC = LC), des), p))
    up_mass <- exp(calc_ctx_ll(build_ll_ctx(
      gng_rows(S, NA_character_, Inf, LT = LT, UT = UT, UC = UC), des), p))
    total <- finite_mass + low_mass + up_mass
    expect_equal(total, 1.0, tolerance = 5e-3,
                 label = sprintf("GNG %s censored+truncated total (%g = %g + %g + %g)",
                                 S, total, finite_mass, low_mass, up_mass))
  }
})
