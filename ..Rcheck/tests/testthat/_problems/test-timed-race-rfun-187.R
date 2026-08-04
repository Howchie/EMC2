# Extracted from test-timed-race-rfun.R:187

# test -------------------------------------------------------------------------
timed_design <- design(
    factors = list(subjects = 1, S = "stim"),
    Rlevels = c("left", "right", "time"),
    matchfun = function(d) as.character(d$S) == as.character(d$lR),
    model = RDMSWTN(),
    formula = list(v ~ 0 + lR, B ~ 1, A ~ 1, t0 ~ 0 + lR, s ~ 1, sv ~ 1),
    report_p_vector = FALSE
  )
dat <- data.frame(
    subjects = factor(1),
    S = factor("stim"),
    R = factor("left", levels = c("left", "right", "time")),
    rt = 0.55,
    LT = 0.25,
    UT = 1.5
  )
emc <- make_emc(dat, timed_design, type = "single", compress = FALSE, n_chains = 1)
model_obj <- emc[[1]]$model()
dadm <- emc[[1]]$data[[1]]
p_vec <- sampled_pars(timed_design, doMap = FALSE)
p_vec[] <- 0
p_vec["v_lRleft"] <- log(3)
p_vec["v_lRright"] <- log(1.2)
p_vec["v_lRtime"] <- log(2)
p_vec["B"] <- log(1)
p_vec["A"] <- log(0)
p_vec["t0_lRleft"] <- log(0.2)
p_vec["t0_lRright"] <- log(0.2)
p_vec["t0_lRtime"] <- log(0.05)
p_vec["s"] <- log(1)
p_vec["sv"] <- log(0)
p_mat <- matrix(p_vec, nrow = 1, dimnames = list(NULL, names(p_vec)))
designs <- lapply(names(model_obj$p_types), function(p) {
    attr(dadm, "designs")[[p]][attr(attr(dadm, "designs")[[p]], "expand"), , drop = FALSE]
  })
names(designs) <- names(model_obj$p_types)
constants <- attr(dadm, "constants")
if (is.null(constants)) constants <- NA
ll_cpp <- EMC2:::calc_ll_oo(
    p_mat, dadm,
    constants = constants,
    designs = designs,
    type = model_obj$c_name,
    bounds = model_obj$bound,
    transforms = model_obj$transform,
    pretransforms = model_obj$pre_transform,
    p_types = names(model_obj$p_types),
    min_ll = log(1e-10),
    trend = model_obj$trend
  )
ll_pw <- EMC2:::calc_ll_oo_pw(
    p_mat, dadm,
    constants = constants,
    designs = designs,
    type = model_obj$c_name,
    bounds = model_obj$bound,
    transforms = model_obj$transform,
    pretransforms = model_obj$pre_transform,
    p_types = names(model_obj$p_types),
    min_ll = log(1e-10),
    trend = model_obj$trend
  )
pars <- as.matrix(mapped_pars(timed_design, p_vec, data = dat)[, names(model_obj$p_types), drop = FALSE])
