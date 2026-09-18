custom_trend_test_data <- data.frame(
  subjects = factor(rep("s1", 8)),
  S = factor(rep(c("left", "right"), 4)),
  R = factor(rep(c("left", "right"), 4), levels = c("left", "right")),
  rt = rep(c(0.55, 0.65), 4),
  distance = rep(c(0, 1), 4)
)

custom_trend_match <- function(d) d$S == d$lR
custom_trend_contrast <- matrix(
  c(-0.5, 0.5), ncol = 1,
  dimnames = list(NULL, "difference")
)

test_that("design functions can be declared static for prediction", {
  trend <- make_trend(
    par_names = "m", cov_names = "race_sign", kernels = "lin_incr",
    at = NULL
  )
  funs <- list(race_sign = function(d) ifelse(d$lR == "right", 1, -1))

  conservative <- design(
    data = custom_trend_test_data,
    model = LNR,
    matchfun = custom_trend_match,
    formula = list(m ~ lM, s ~ 1, t0 ~ 1),
    contrasts = list(lM = custom_trend_contrast),
    functions = funs,
    trend = trend
  )
  static <- design(
    data = custom_trend_test_data,
    model = LNR,
    matchfun = custom_trend_match,
    formula = list(m ~ lM, s ~ 1, t0 ~ 1),
    contrasts = list(lM = custom_trend_contrast),
    functions = funs,
    behavioral_functions = character(0),
    trend = trend
  )

  expect_true(EMC2:::has_conditional_covariates(conservative))
  expect_false(EMC2:::has_conditional_covariates(static))
  expect_identical(static$Fbehavioral, character(0))
})

test_that("behavioral_functions validates function output names", {
  expect_error(
    design(
      data = custom_trend_test_data,
      model = LNR,
      matchfun = custom_trend_match,
      formula = list(m ~ lM, s ~ 1, t0 ~ 1),
      contrasts = list(lM = custom_trend_contrast),
      functions = list(race_sign = function(d) 1),
      behavioral_functions = "missing"
    ),
    "must name entries in functions"
  )
})

test_that("pointer restoration is a no-op for ordinary joint models", {
  plain_model <- function() list(trend = NULL)
  object <- list(list(model = list(plain_model, plain_model)))
  class(object) <- "emc"
  expect_identical(restore_custom_kernel_pointers(object), object)
})

test_that("serialized custom trends restore their external pointers", {
  source_file <- testthat::test_path(
    "fixtures", "serialized_custom_trend.cpp"
  )
  registered <- register_trend(
    trend_parameters = "slope",
    file = source_file,
    transforms = c(slope = "identity"),
    base = "add"
  )
  trend <- make_trend(
    par_names = "m", cov_names = "distance", kernels = "custom",
    custom_trend = registered, at = NULL
  )
  registration <- attr(trend[[1]], "custom_registration")
  expect_identical(registration$file, normalizePath(source_file))
  expect_identical(registration$trend_parameters, "slope")

  des <- design(
    data = custom_trend_test_data,
    model = LNR,
    matchfun = custom_trend_match,
    formula = list(m ~ lM, s ~ 1, t0 ~ 1),
    contrasts = list(lM = custom_trend_contrast),
    trend = trend
  )
  object <- make_emc(
    custom_trend_test_data, des, compress = FALSE,
    n_chains = 1, type = "single"
  )
  object <- unserialize(serialize(object, NULL))
  dead_ptr <- attr(object[[1]]$model()$trend[[1]], "custom_ptr")
  expect_false(EMC2:::.custom_kernel_pointer_valid(dead_ptr))

  restored <- restore_custom_kernel_pointers(object, quiet = TRUE)
  live_ptr <- attr(restored[[1]]$model()$trend[[1]], "custom_ptr")
  design_ptr <- attr(get_design(restored)[[1]]$model()$trend[[1]], "custom_ptr")
  expect_true(EMC2:::.custom_kernel_pointer_valid(live_ptr))
  expect_true(EMC2:::.custom_kernel_pointer_valid(design_ptr))

  mapped_p <- sampled_pars(get_design(object)[[1]], doMap = FALSE)
  mapped_p[] <- 0
  mapped_p["s"] <- log(0.3)
  mapped_p["t0"] <- log(0.2)
  expect_no_error(mapped_pars(object, p_vector = mapped_p))

  restored_design <- get_design(restored)[[1]]
  p <- sampled_pars(restored_design, doMap = FALSE)
  p[] <- 0
  p["s"] <- log(0.3)
  p["t0"] <- log(0.2)
  simulated <- make_data(
    p, design = restored_design, data = custom_trend_test_data,
    conditional_on_data = TRUE
  )
  expect_s3_class(simulated, "data.frame")
  expect_equal(nrow(simulated), nrow(custom_trend_test_data))
})

test_that("serialized joint custom trends restore in every model", {
  source_file <- testthat::test_path(
    "fixtures", "serialized_custom_trend.cpp"
  )
  registered <- register_trend(
    trend_parameters = "slope", file = source_file,
    transforms = c(slope = "identity"), base = "add"
  )
  trend <- make_trend(
    par_names = "m", cov_names = "distance", kernels = "custom",
    custom_trend = registered, at = NULL
  )
  model <- design(
    data = custom_trend_test_data, model = LNR,
    matchfun = custom_trend_match,
    formula = list(m ~ lM, s ~ 1, t0 ~ 1),
    contrasts = list(lM = custom_trend_contrast), trend = trend
  )$model
  object <- list(list(model = list(model, model)))
  class(object) <- "emc"
  object <- unserialize(serialize(object, NULL))
  restored <- restore_custom_kernel_pointers(object, quiet = TRUE)
  ptrs <- get_custom_kernel_pointers(restored)
  expect_length(ptrs, 2)
  expect_true(all(vapply(ptrs, function(x)
    EMC2:::.custom_kernel_pointer_valid(x[[1]]), logical(1))))
})

test_that("unconditional custom trends refresh behavioral design functions", {
  source_file <- testthat::test_path(
    "fixtures", "serialized_custom_trend.cpp"
  )
  registered <- register_trend(
    trend_parameters = "slope", file = source_file,
    transforms = c(slope = "identity"), base = "add"
  )
  trend <- make_trend(
    par_names = "m", cov_names = "feedback", kernels = "custom",
    custom_trend = registered, at = NULL
  )
  data <- custom_trend_test_data[, c("subjects", "S", "R", "rt")]
  des <- design(
    data = data, model = LNR, matchfun = custom_trend_match,
    formula = list(m ~ lM, s ~ 1, t0 ~ 1),
    contrasts = list(lM = custom_trend_contrast),
    functions = list(feedback = function(d) as.numeric(d$R == "left")),
    trend = trend
  )
  p <- sampled_pars(des, doMap = FALSE)
  p[] <- 0
  p[c("s", "t0", "m.slope")] <- c(log(0.3), log(0.2), 0.8)
  set.seed(42)
  scalar <- make_data(p, des, data = data, conditional_on_data = FALSE,
                      return_functions = TRUE)
  set.seed(42)
  automatic <- make_data(p, des, data = data, return_functions = TRUE)
  set.seed(42)
  vectorised <- make_data(p, des, data = data, conditional_on_data = FALSE,
                          use_vectorised = TRUE, return_functions = TRUE)
  expect_true(all(scalar$feedback == as.numeric(scalar$R == "left")))
  expect_identical(scalar, automatic)
  expect_identical(scalar, vectorised)

  static <- design(
    data = data, model = LNR, matchfun = custom_trend_match,
    formula = list(m ~ lM, s ~ 1, t0 ~ 1),
    contrasts = list(lM = custom_trend_contrast),
    functions = list(feedback = function(d) as.numeric(d$R == "left")),
    behavioral_functions = character(0), trend = trend
  )
  set.seed(42)
  static_automatic <- make_data(p, static, data = data, return_functions = TRUE)
  set.seed(42)
  static_conditional <- make_data(p, static, data = data,
                                  conditional_on_data = TRUE,
                                  return_functions = TRUE)
  expect_identical(static_automatic, static_conditional)

  multi_data <- data
  multi_data$subjects <- factor(rep(c("s1", "s2"), each = nrow(data) / 2))
  multi_des <- design(
    data = multi_data, model = LNR, matchfun = custom_trend_match,
    formula = list(m ~ lM, s ~ 1, t0 ~ 1),
    contrasts = list(lM = custom_trend_contrast),
    functions = list(feedback = function(d) as.numeric(d$R == "left")),
    trend = trend
  )
  multi_p <- sampled_pars(multi_des, doMap = FALSE)
  multi_p <- rbind(s1 = multi_p, s2 = multi_p)
  multi_p[] <- 0
  multi_p[, c("s", "t0", "m.slope")] <-
    matrix(rep(c(log(0.3), log(0.2), 0.8), each = nrow(multi_p)),
           nrow = nrow(multi_p))
  set.seed(42)
  multi_vectorised <- make_data(
    multi_p, multi_des, data = multi_data, conditional_on_data = FALSE,
    use_vectorised = TRUE, return_functions = TRUE
  )
  expect_equal(nrow(multi_vectorised), nrow(multi_data))
  expect_true(all(multi_vectorised$feedback ==
                  as.numeric(multi_vectorised$R == "left")))
})

test_that("legacy objects accept a registered custom trend as pointer source", {
  source_file <- testthat::test_path(
    "fixtures", "serialized_custom_trend.cpp"
  )
  registered <- register_trend(
    trend_parameters = "slope", file = source_file,
    transforms = c(slope = "identity"), base = "add"
  )
  trend <- make_trend(
    par_names = "m", cov_names = "distance", kernels = "custom",
    custom_trend = registered, at = NULL
  )
  des <- design(
    data = custom_trend_test_data,
    model = LNR,
    matchfun = custom_trend_match,
    formula = list(m ~ lM, s ~ 1, t0 ~ 1),
    contrasts = list(lM = custom_trend_contrast),
    trend = trend
  )
  object <- make_emc(
    custom_trend_test_data, des, compress = FALSE,
    n_chains = 1, type = "single"
  )
  object <- unserialize(serialize(object, NULL))
  object[[1]]$model <- local({
    model_list <- object[[1]]$model()
    attr(model_list$trend[[1]], "custom_registration") <- NULL
    function() model_list
  })

  restored <- fix_custom_kernel_pointers(object, registered)
  ptr <- attr(restored[[1]]$model()$trend[[1]], "custom_ptr")
  expect_true(EMC2:::.custom_kernel_pointer_valid(ptr))
})
