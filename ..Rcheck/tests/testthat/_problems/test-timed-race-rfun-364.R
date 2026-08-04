# Extracted from test-timed-race-rfun.R:364

# test -------------------------------------------------------------------------
skip_if_not_installed("dplyr")
dat <- expand.grid(
    subjects = factor(1),
    sa = factor(c("speed", "neutral", "accuracy"), levels = c("speed", "neutral", "accuracy")),
    S = factor(c("left", "right"), levels = c("left", "right")),
    rep = 1:2
  )
dat$R <- factor(dat$S, levels = c("left", "right", "time"))
dat$rt <- seq(0.35, 0.9, length.out = nrow(dat))
dat$LT <- 0.25
dat$UT <- 1.5
matchfun <- function(d) as.numeric(d$lR) == as.numeric(d$S)
timed_design <- design(
    data = dat,
    factors = list(Rlevels = c("left", "right", "time")),
    model = RDMSWTN(),
    transform = list(func = c(v = "exp")),
    matchfun = matchfun,
    functions = list(
      match = function(d) {
        dplyr::case_when(
          d$lM == TRUE & !(d$lR == "time") ~ 0.5,
          d$lM == FALSE & !(d$lR == "time") ~ -0.5,
          d$lR == "time" ~ 0
        )
      },
      E = function(d) ifelse(d$lR == "time", 0, 1),
      Time = function(d) ifelse(d$lR == "time", 1, 0),
      Resp = function(d) {
        dplyr::case_when(
          d$lR == "left" ~ 0.5,
          d$lR == "right" ~ -0.5,
          d$lR == "time" ~ 0
        )
      }
    ),
    formula = list(
      v ~ 0 + sa:E + match + sa:Time,
      B ~ 0 + sa:E + Resp + Time,
      t0 ~ 0 + (E:sa + Time),
      s ~ 0 + (E + Time),
      A ~ 1,
      sv ~ 1
    ),
    constants = c(
      A = log(0), sv = log(0),
      "B_saspeed:E" = log(1), "B_Time" = log(1), t0_Time = log(0.05)
    ),
    pre_transform_terms = list(
      B = c("B_saspeed:E", "B_saneutral:E", "B_saaccuracy:E", "B_Time"),
      v = c("v_saspeed:E", "v_saneutral:E", "v_saaccuracy:E",
            "v_saspeed:Time", "v_saneutral:Time", "v_saaccuracy:Time")
    ),
    bound = list(minmax = cbind(v = c(0, Inf))),
    report_p_vector = FALSE
  )
p_vec <- sampled_pars(timed_design, doMap = FALSE)
vals <- c(c(3.99), log(c(2.2, 1.89, 1.59, 1.66, 1.12, 0.86)),
            -0.11, log(c(1.05, 1.08, 0.24, 0.29, 0.29, 1.33, 0.51)))
names(vals) <- names(p_vec)
p_vec[names(vals)] <- vals
emc <- make_emc(dat, timed_design, type = "single", compress = FALSE, n_chains = 1)
dadm <- emc[[1]]$data[[1]]
model_obj <- emc[[1]]$model()
designs <- EMC2:::.oo_expanded_designs(dadm)
constants <- attr(dadm, "constants")
if (is.null(constants)) constants <- NA
p_mat <- matrix(p_vec, nrow = 1, dimnames = list(NULL, names(p_vec)))
cpp <- EMC2:::get_pars_c_wrapper_oo(
    p_mat, dadm, constants, designs, model_obj$bound, model_obj$transform,
    model_obj$pre_transform, model_obj$trend, FALSE, FALSE, 1L
  )
r_map <- mapped_pars(timed_design, p_vec, data = dat, digits = 10)
