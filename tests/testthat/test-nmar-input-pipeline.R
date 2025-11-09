test_that("parse_nmar_spec captures core pieces", {
  df <- data.frame(Y = c(1, NA, 3), X = 1:3, Z = 4:6)
  spec <- NMAR:::parse_nmar_spec(Y ~ X | Z + Z, df)
  expect_s3_class(spec, "nmar_input_spec")
  expect_equal(spec$outcome, "Y")
  expect_equal(spec$auxiliary_vars, "X")
  expect_equal(spec$response_predictors, "Z")
  expect_false(spec$is_survey)
  expect_true(is.data.frame(spec$data))
})

test_that("parse_nmar_spec works for survey designs", {
  skip_if_not_installed("survey")
  df <- data.frame(Y = c(1, NA, 3), X = rnorm(3), w = c(1, 1.5, 0.5))
  design <- survey::svydesign(ids = ~1, data = df, weights = ~w)
  spec <- NMAR:::parse_nmar_spec(Y ~ X, design)
  expect_true(spec$is_survey)
  expect_s3_class(spec$original_data, "survey.design")
  expect_equal(nrow(spec$data), nrow(df))
})

test_that("validate_nmar_args enforces defaults", {
  df <- data.frame(Y = c(1, NA, 3, 4), X = rnorm(4), Z = rnorm(4))
  spec <- NMAR:::parse_nmar_spec(Y ~ X | X, df)
  expect_error(NMAR:::validate_nmar_args(spec, list()), "mutually exclusive", fixed = FALSE)
})

test_that("response-side dot expansion enforces overlap rules", {
  df <- data.frame(Y = c(1, NA, 2), X = rnorm(3), Z = rnorm(3))
  spec <- NMAR:::parse_nmar_spec(Y ~ X + Z | ., df)
  expect_error(
    NMAR:::validate_nmar_args(spec, list()),
    "mutually exclusive",
    fixed = FALSE
  )
})

test_that("validate_nmar_args can be relaxed for EL", {
  df <- data.frame(Y = c(1, NA, 3, 4), X = rnorm(4), Z = rnorm(4))
  eng_obj <- el_engine(auxiliary_means = c(X = 0), variance_method = "none")
  spec <- NMAR:::parse_nmar_spec(Y ~ X | Y, df, traits = NMAR:::engine_traits(eng_obj))
  traits <- NMAR:::engine_traits(eng_obj)
  expect_silent(NMAR:::validate_nmar_args(spec, traits))
  expect_equal(spec$response_predictors_raw, "Y")
})

test_that("default traits block outcome on response model", {
  df <- data.frame(Y = c(1, NA, 2), X = rnorm(3))
  spec <- NMAR:::parse_nmar_spec(Y ~ 1 | Y, df)
  expect_error(
    NMAR:::validate_nmar_args(spec, list()),
    "Outcome variable cannot appear",
    fixed = FALSE
  )
})

test_that("nonparam engine accepts multi-outcome formulas", {
  df <- data.frame(
    Voted_A = c(10, 12),
    Voted_B = c(5, 8),
    Other = c(2, 1),
    Gender = c(0, 1),
    Refusal = c(3, 4)
  )
  traits_np <- NMAR:::engine_traits(exptilt_nonparam_engine(refusal_col = "Refusal"))
  spec <- NMAR:::parse_nmar_spec(Voted_A + Voted_B + Other ~ Gender, df, traits = traits_np)
  expect_silent(NMAR:::validate_nmar_args(spec, traits_np))
  traits_default <- NMAR:::engine_traits(exptilt_engine())
  expect_error(NMAR:::validate_nmar_args(spec, traits_default))
})

test_that("single-outcome traits reject multi-outcome formulas", {
  df <- data.frame(
    Y1 = c(1, NA, 2),
    Y2 = c(3, 4, NA),
    X = rnorm(3)
  )
  spec <- NMAR:::parse_nmar_spec(cbind(Y1, Y2) ~ X, df, traits = NMAR:::NMAR_DEFAULT_TRAITS)
  traits <- NMAR:::NMAR_DEFAULT_TRAITS
  expect_error(
    NMAR:::validate_nmar_args(spec, traits),
    "exactly one outcome",
    fixed = FALSE
  )
})

test_that("parse_nmar_spec resolves dots and factors", {
  df <- data.frame(
    Y = c(1, NA, 3),
    X = 1:3,
    Z = rnorm(3),
    F = factor(c("a", "b", "a"))
  )
  spec <- NMAR:::parse_nmar_spec(Y ~ . | F, df)
  expect_setequal(spec$auxiliary_vars, c("X", "Z", "F"))
  expect_equal(spec$response_predictors, "F")
  traits <- NMAR:::engine_traits(el_engine(auxiliary_means = c(X = 0)))
  expect_silent(NMAR:::validate_nmar_args(spec, traits))
})

test_that("partition rebuild preserves formula environments", {
  local_fn <- function(x) x^2
  env <- environment()
  df <- data.frame(
    Y = c(1, 2, NA, 4),
    X = rnorm(4)
  )
  formula <- Y ~ local_fn(X)
  environment(formula) <- env
  spec <- NMAR:::parse_nmar_spec(formula, df)
  traits <- NMAR:::engine_traits(el_engine(auxiliary_means = NULL, variance_method = "none"))
  NMAR:::validate_nmar_args(spec, traits)
  task <- NMAR:::new_nmar_task(spec, traits)
  design <- NMAR:::prepare_nmar_design(task)
  rebuilt <- design$user_formula
  expect_identical(attr(rebuilt, ".Environment"), env)
  mm <- stats::model.matrix(rebuilt, df)
  expect_true(any(grepl("local_fn", colnames(mm))))
})

test_that("parse_nmar_spec resolves scalars from the formula environment", {
  df <- data.frame(Y = c(1, NA, 3), X = rnorm(3))
  degree <- 2
  fml <- Y ~ poly(X, degree)
  environment(fml) <- environment()
  spec <- NMAR:::parse_nmar_spec(fml, df)
  expect_s3_class(spec, "nmar_input_spec")
  traits <- NMAR:::engine_traits(el_engine(auxiliary_means = NULL, variance_method = "none"))
  expect_silent(NMAR:::validate_nmar_args(spec, traits))
  expect_equal(spec$auxiliary_vars, "X")
})

test_that("new_nmar_task preserves raw predictor metadata", {
  df <- data.frame(Y = c(1, NA, 3), X = 1:3, Z = rnorm(3))
  eng_obj2 <- el_engine(auxiliary_means = c(X = 0), variance_method = "none")
  spec <- NMAR:::parse_nmar_spec(Y ~ X | Z + Y, df, traits = NMAR:::engine_traits(eng_obj2))
  traits <- NMAR:::engine_traits(eng_obj2)
  task <- NMAR:::new_nmar_task(spec, traits)
  expect_equal(task$response_predictors_raw, c("Z", "Y"))
  expect_equal(task$auxiliary_vars_raw, "X")
})

test_that("outcome transformations must be numeric vectors", {
  df <- data.frame(Y = c(1, 2, NA, 4), X = rnorm(4))
  spec <- NMAR:::parse_nmar_spec(I(as.character(Y)) ~ X, df)
  traits <- NMAR:::engine_traits(el_engine(auxiliary_means = NULL, variance_method = "none"))
  NMAR:::validate_nmar_args(spec, traits)
  task <- NMAR:::new_nmar_task(spec, traits)
  expect_error(
    NMAR:::prepare_nmar_design(task),
    "numeric vector",
    fixed = FALSE
  )
})

test_that("outcome transformations respect helper-first symbols", {
  df <- data.frame(
    Y = c(1, NA, 3),
    denom = c(2, 2, 2),
    X = rnorm(3)
  )
  spec <- NMAR:::parse_nmar_spec(I(denom - Y) ~ X, df)
  traits <- NMAR:::engine_traits(el_engine(variance_method = "none"))
  expect_silent(NMAR:::validate_nmar_args(spec, traits))
})

test_that("multi-column outcome transformations handle helper variables", {
  scale_val <- 2
  df <- data.frame(
    Y = c(1, NA, 3),
    X = rnorm(3)
  )
  fml <- cbind(Y, scale_val - Y) ~ X
  environment(fml) <- environment()
  traits_multi <- utils::modifyList(
    NMAR:::NMAR_DEFAULT_TRAITS,
    list(requires_single_outcome = FALSE)
  )
  spec <- NMAR:::parse_nmar_spec(fml, df, traits = traits_multi)
  expect_silent(NMAR:::validate_nmar_args(spec, traits_multi))
  task <- NMAR:::new_nmar_task(spec, traits_multi)
  design <- NMAR:::prepare_nmar_design(task)
  expect_equal(length(design$outcome_columns), 2)
  expect_true(all(design$outcome_columns %in% names(design$data)))
})

test_that("multi-column outcomes receive deterministic, valid column names", {
  helper_matrix <- function(v) {
    mat <- cbind(v, v * 2)
    colnames(mat) <- c("", NA_character_)
    mat
  }
  env <- environment()
  df <- data.frame(
    Y = c(1, NA, 3),
    X = rnorm(3)
  )
  fml <- cbind(Y, helper_matrix(Y)) ~ X
  environment(fml) <- env
  traits_multi <- utils::modifyList(
    NMAR:::NMAR_DEFAULT_TRAITS,
    list(requires_single_outcome = FALSE)
  )
  spec <- NMAR:::parse_nmar_spec(fml, df, traits = traits_multi)
  NMAR:::validate_nmar_args(spec, traits_multi)
  task <- NMAR:::new_nmar_task(spec, traits_multi)
  design <- NMAR:::prepare_nmar_design(task)
  expect_true(all(nchar(design$outcome_columns) > 0))
  expect_equal(anyDuplicated(design$outcome_columns), 0L)
  expect_true(any(grepl("helper_matrix", design$outcome_columns, fixed = TRUE)))
  expect_true(all(design$outcome_columns %in% names(design$data)))
})

test_that("engine formulas reference the sanitized outcome column names", {
  helper_tbl <- function(v) {
    cbind(named_helper = v + 1)
  }
  env <- environment()
  df <- data.frame(
    Y_main = c(1, NA, 4, NA),
    X = rnorm(4)
  )
  fml <- cbind(Y_main, helper_tbl(Y_main), Y_main * 0 + 3) ~ X
  environment(fml) <- env
  traits_multi <- utils::modifyList(
    NMAR:::NMAR_DEFAULT_TRAITS,
    list(requires_single_outcome = FALSE)
  )
  spec <- NMAR:::parse_nmar_spec(fml, df, traits = traits_multi)
  NMAR:::validate_nmar_args(spec, traits_multi)
  task <- NMAR:::new_nmar_task(spec, traits_multi)
  design <- NMAR:::prepare_nmar_design(task)
  lhs_vars <- all.vars(design$engine_formula[[2]])
  expect_setequal(lhs_vars, design$outcome_columns)
  expect_true(length(lhs_vars) == length(design$outcome_columns))
})

test_that("survey designs pick up derived outcome columns", {
  skip_if_not_installed("survey")
  df <- data.frame(
    Y = c(1, 2, NA, 4),
    X = rnorm(4),
    w = c(1, 1.5, 0.5, 2)
  )
  design <- survey::svydesign(ids = ~1, data = df, weights = ~w)
  spec <- NMAR:::parse_nmar_spec(log(Y + 10) ~ X, design)
  traits <- NMAR:::engine_traits(el_engine(auxiliary_means = NULL, variance_method = "none"))
  NMAR:::validate_nmar_args(spec, traits)
  task <- NMAR:::new_nmar_task(spec, traits)
  info <- NMAR:::prepare_nmar_design(task)
  expect_true(info$is_survey)
  expect_true(info$outcome_column %in% names(info$survey_design$variables))
  expect_equal(info$survey_design$variables[[info$outcome_column]], info$data[[info$outcome_column]])
})

test_that("respondents-only survey designs require auxiliary metadata", {
  skip_if_not_installed("survey")
  df <- data.frame(
    Y = rnorm(5, mean = 5),
    X = rnorm(5),
    w = c(1, 1.5, 0.7, 1.2, 0.9)
  )
  design <- survey::svydesign(ids = ~1, data = df, weights = ~w)
  engine <- el_engine(
    auxiliary_means = c(X = 0),
    n_total = sum(df$w),
    variance_method = "none"
  )
  traits <- NMAR:::engine_traits(engine)
  spec <- NMAR:::parse_nmar_spec(log(Y) ~ X, design, traits = traits)
  expect_silent(NMAR:::validate_nmar_args(spec, traits))
  task <- NMAR:::new_nmar_task(spec, traits)
  info <- NMAR:::prepare_nmar_design(task, auxiliary_means = engine$auxiliary_means)
  expect_true(info$is_survey)
  expect_false(anyNA(info$data[[info$outcome_column]]))
  expect_equal(
    unname(info$weights),
    as.numeric(stats::weights(info$survey_design))
  )
  expect_equal(
    info$survey_design$variables[[info$outcome_column]],
    info$data[[info$outcome_column]]
  )
})

test_that("multi-outcome validation enforces trait restrictions", {
  df <- data.frame(
    Y1 = c(1, 2, 3),
    Y2 = c(4, 5, 6),
    X = c(0.1, 0.2, 0.3),
    Z = c(1, 1, 0)
  )
  traits_np <- NMAR:::engine_traits(exptilt_nonparam_engine(refusal_col = "Z"))
  spec <- NMAR:::parse_nmar_spec(Y1 + Y2 ~ X | Z, df, traits = traits_np)
  expect_silent(NMAR:::validate_nmar_args(spec, traits_np))
  traits_block <- NMAR:::engine_traits(exptilt_engine())
  expect_error(NMAR:::validate_nmar_args(spec, traits_block))
})

test_that("multi-outcome respondents-only inputs honor trait flags", {
  df <- data.frame(
    Y1 = rnorm(4),
    Y2 = rnorm(4),
    X = rnorm(4)
  )
  traits_block <- utils::modifyList(
    NMAR:::NMAR_DEFAULT_TRAITS,
    list(requires_single_outcome = FALSE, allow_respondents_only = FALSE)
  )
  spec <- NMAR:::parse_nmar_spec(Y1 + Y2 ~ X, df, traits = traits_block)
  expect_error(
    NMAR:::validate_nmar_args(spec, traits_block),
    "must contain NA",
    fixed = FALSE
  )

  traits_allow <- utils::modifyList(
    traits_block,
    list(allow_respondents_only = TRUE)
  )
  spec2 <- NMAR:::parse_nmar_spec(Y1 + Y2 ~ X, df, traits = traits_allow)
  expect_silent(NMAR:::validate_nmar_args(spec2, traits_allow))
})

test_that("respondents-only inputs require allow_respondents_only trait", {
  df <- data.frame(Y = rnorm(5), X = rnorm(5))
  spec <- NMAR:::parse_nmar_spec(Y ~ X, df)
  traits_default <- NMAR:::engine_traits(exptilt_engine())
  expect_error(NMAR:::validate_nmar_args(spec, traits_default), "must contain NA", fixed = FALSE)
  traits_el <- NMAR:::engine_traits(el_engine(auxiliary_means = c(X = 0), n_total = 100, variance_method = "none"))
  expect_silent(NMAR:::validate_nmar_args(spec, traits_el))
})

test_that("environment-only auxiliary predictors are accepted", {
  df <- data.frame(Y = c(1, NA, 3))
  offset_vec <- rnorm(3)
  env <- environment()
  fml <- Y ~ offset_vec
  environment(fml) <- env
  spec <- NMAR:::parse_nmar_spec(fml, df)
  traits <- NMAR:::engine_traits(el_engine(variance_method = "none"))
  expect_silent(NMAR:::validate_nmar_args(spec, traits))
  task <- NMAR:::new_nmar_task(spec, traits)
  design <- NMAR:::prepare_nmar_design(task)
  expect_true("offset_vec" %in% colnames(design$design_matrices$auxiliary))
})

test_that("environment-only response predictors are accepted", {
  df <- data.frame(Y = c(1, NA, 3), X = rnorm(3))
  offset_resp <- rnorm(3)
  env <- environment()
  fml <- Y ~ X | offset_resp
  environment(fml) <- env
  spec <- NMAR:::parse_nmar_spec(fml, df)
  traits <- NMAR:::engine_traits(el_engine(variance_method = "none"))
  expect_silent(NMAR:::validate_nmar_args(spec, traits))
  task <- NMAR:::new_nmar_task(spec, traits)
  design <- NMAR:::prepare_nmar_design(task)
  resp_mm <- design$design_matrices$response
  expect_true(is.matrix(resp_mm))
  expect_true("offset_resp" %in% colnames(resp_mm))
})

test_that("outcome selection favors columns with missingness", {
  df <- data.frame(
    Y_obs = c(1, NA, 3),
    helper = 5:7,
    X = rnorm(3)
  )
  fml <- I(helper - Y_obs) ~ X
  spec <- NMAR:::parse_nmar_spec(fml, df)
  traits <- NMAR:::engine_traits(el_engine(variance_method = "none"))
  expect_silent(NMAR:::validate_nmar_args(spec, traits))
  expect_identical(spec$outcome_primary, "Y_obs")
  expect_identical(spec$outcome_helpers, "helper")
})

test_that("survey designs keep environment-only auxiliaries", {
  skip_if_not_installed("survey")
  df <- data.frame(
    Y = c(1, 2, NA, 4),
    w = c(1, 2, 1, 2)
  )
  offset_vec <- rnorm(nrow(df))
  env <- environment()
  fml <- Y ~ offset_vec
  environment(fml) <- env
  design <- survey::svydesign(ids = ~1, data = df, weights = ~w)
  spec <- NMAR:::parse_nmar_spec(fml, design)
  traits <- NMAR:::engine_traits(el_engine(variance_method = "none"))
  expect_silent(NMAR:::validate_nmar_args(spec, traits))
  task <- NMAR:::new_nmar_task(spec, traits)
  info <- NMAR:::prepare_nmar_design(task)
  aux_mm <- info$design_matrices$auxiliary
  expect_true("offset_vec" %in% colnames(aux_mm))
  expect_equal(
    info$survey_design$variables[[info$outcome_column]],
    info$data[[info$outcome_column]]
  )
})

test_that("outcome cannot enter response RHS via helper-first transforms", {
  df <- data.frame(
    Y = c(1, NA, 3),
    denom = c(2, 2, 2),
    X = rnorm(3)
  )
  fml <- I(denom - Y) ~ X | Y
  traits <- NMAR:::engine_traits(exptilt_engine())
  spec <- NMAR:::parse_nmar_spec(fml, df, traits = traits)
  expect_error(
    NMAR:::validate_nmar_args(spec, traits),
    "Outcome variable cannot appear",
    fixed = FALSE
  )
})

test_that("Y ~ 1 | X implies no auxiliaries and response includes outcome", {
  df <- data.frame(Y = c(1, NA, 3), X = rnorm(3))
  spec <- NMAR:::parse_nmar_spec(Y ~ 1 | X, df)
# No auxiliaries when RHS aux is `1`
  expect_length(spec$auxiliary_vars, 0)
  expect_equal(spec$response_predictors, "X")
# Internal response formula includes outcome implicitly
  split <- NMAR:::nmar_split_partitioned_formula(NMAR:::nmar_rebuild_partitioned_formula(spec$formula, spec$response_rhs_lang, spec$aux_rhs_lang, env = spec$environment))
  forms <- NMAR:::nmar_build_internal_formulas(delta_name = "..d..", outcome_var = split$outcome_var, aux_rhs_lang = split$aux_rhs_lang, response_rhs_lang = split$response_rhs_lang, env = split$env)
  resp_str <- paste(deparse(forms$response), collapse = " ")
  expect_true(grepl("..d.. ~ Y + X", resp_str, fixed = TRUE))
})

test_that("No bar: Y ~ X has empty explicit response predictors; response is delta ~ Y", {
  df <- data.frame(Y = c(1, NA, 3), X = rnorm(3))
  spec <- NMAR:::parse_nmar_spec(Y ~ X, df)
  expect_length(spec$response_predictors_raw, 0)
  split <- NMAR:::nmar_split_partitioned_formula(spec$formula)
  forms <- NMAR:::nmar_build_internal_formulas(delta_name = "..d..", outcome_var = split$outcome_var, aux_rhs_lang = split$aux_rhs_lang, response_rhs_lang = split$response_rhs_lang, env = split$env)
  resp_str <- paste(deparse(forms$response), collapse = " ")
  expect_equal(resp_str, "..d.. ~ Y")
})

test_that("Blueprint source_variables reflect RHS only", {
  df <- data.frame(Y = c(1, NA, 3), X = rnorm(3), Z = rnorm(3))
  bp <- NMAR:::nmar_build_formula_blueprint(outcome_vars = "Y", aux_expr = quote(X), response_expr = quote(Y + Z), data = df, env = parent.frame())
  expect_equal(bp$aux$source_variables, "X")
# Explicit RHS may include Y; ensure only RHS names appear
  expect_setequal(bp$response$source_variables, c("Y", "Z"))
})

test_that("Auxiliary intercepts follow trait flag", {
  df <- data.frame(Y = c(1, NA, 3), X = rnorm(3))
  traits_drop <- NMAR:::engine_traits(el_engine(variance_method = "none"))
  spec_drop <- NMAR:::parse_nmar_spec(Y ~ 1 + X, df, traits = traits_drop)
  task_drop <- NMAR:::new_nmar_task(spec_drop, traits_drop)
  design_drop <- NMAR:::prepare_nmar_design(task_drop)
  aux_mm <- design_drop$design_matrices$auxiliary
  expect_equal(colnames(aux_mm), "X")

  traits_keep <- NMAR:::NMAR_DEFAULT_TRAITS
  traits_keep$drop_auxiliary_intercept <- FALSE
  spec_keep <- NMAR:::parse_nmar_spec(Y ~ 1 + X, df, traits = traits_keep)
  task_keep <- NMAR:::new_nmar_task(spec_keep, traits_keep)
  design_keep <- NMAR:::prepare_nmar_design(task_keep)
  aux_keep <- design_keep$design_matrices$auxiliary
  expect_true("(Intercept)" %in% colnames(aux_keep))
  expect_true("X" %in% colnames(aux_keep))
})

test_that("LHS helpers are allowed in outcome transformations", {
  df <- data.frame(Y = c(1, NA, 3), denom = 2, X = rnorm(3))
  spec <- NMAR:::parse_nmar_spec(I(Y / denom) ~ X, df)
  traits <- NMAR:::engine_traits(el_engine(variance_method = "none"))
  expect_silent(NMAR:::validate_nmar_args(spec, traits))
  task <- NMAR:::new_nmar_task(spec, traits)
  design <- NMAR:::prepare_nmar_design(task)
  expect_equal(spec$outcome_primary, "Y")
  expect_equal(spec$outcome_helpers, "denom")
  derived <- design$data[[design$outcome_column]]
  expect_true(anyNA(derived))
})

test_that("blueprint caches survive user-side data mutation", {
  df <- data.frame(
    Y = c(1, NA, 3, NA),
    X = rnorm(4),
    Z = rnorm(4)
  )
  spec <- NMAR:::parse_nmar_spec(Y ~ X + Z, df)
  traits <- NMAR:::engine_traits(el_engine(auxiliary_means = c(X = 0)))
  NMAR:::validate_nmar_args(spec, traits)

# Mutate the original data frame after parsing; spec$data should retain the old snapshot
  df$X <- NULL
  task <- NMAR:::new_nmar_task(spec, traits)
  design <- NMAR:::prepare_nmar_design(task)
  expect_true("X" %in% colnames(design$design_matrices$auxiliary))

# If spec$data itself is tampered with, the shared pipeline should fail loudly
  spec$data$X <- NULL
  task_broken <- NMAR:::new_nmar_task(spec, traits)
  expect_error(
    NMAR:::prepare_nmar_design(task_broken),
    "not found",
    fixed = FALSE
  )
})
