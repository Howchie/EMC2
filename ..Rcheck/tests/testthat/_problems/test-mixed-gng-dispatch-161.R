# Extracted from test-mixed-gng-dispatch.R:161

# prequel ----------------------------------------------------------------------
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
set_detection_params <- function(p_vec, values) {
  for (nm in names(values)) {
    if (nm %in% names(p_vec)) p_vec[[nm]] <- values[[nm]]
  }
  p_vec
}
detection_funcs <- list(
  GoA = function(d) ifelse(d$lR == "A", 1, 0),
  GoB = function(d) ifelse(d$lR == "B", 1, 0),
  NoGo = function(d) ifelse(d$lR == "nogo", 1, 0)
)
logical_rule_funcs <- c(
  detection_funcs,
  list(
    NegA = function(d) ifelse(d$lR == "n_A", 1, 0),
    NegB = function(d) ifelse(d$lR == "n_B", 1, 0)
  )
)

# test -------------------------------------------------------------------------
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
