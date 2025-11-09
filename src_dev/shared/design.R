#' Input parsing and design preparation
#'
#' `parse_nmar_spec()` turns a user formula and data object into a normalized
#' `nmar_input_spec`. RHS partitions are separated, the original environments are
#' preserved, and a blueprint of `terms` objects is cached so that
#' `model.matrix()` can be re-run later without re-parsing the formula. Outcome
#' symbols are tagged by provenance (data column vs environment helper), and
#' single-outcome workflows pick the first outcome column that actually contains
#' missingness so that helper-first expressions still resolve to the intended
#' estimand.
#'
#' `prepare_nmar_design()` uses the spec plus engine traits to produce a
#' method-agnostic `design_info` bundle containing
#' * Derived outcome columns (LHS transforms are evaluated once and stored)
#' * Consistent weights for data frames or survey designs
#' * Design matrices for auxiliaries/response predictors with intercept handling
#'   controlled by engine traits
#' * User-facing and engine-ready formulas so diagnostics can report the input
#'   verbatim while engines consume the normalized version
#'
#' Non-finite values produced by LHS transformations are allowed only on rows
#' where the input variables are missing (nonrespondents); for observed rows we
#' fail fast to avoid silent NA propagation.
#'
#' @keywords internal

parse_nmar_spec <- function(formula, data, env = parent.frame(), traits = NMAR_DEFAULT_TRAITS) {
  validator$assert_formula_two_sided(formula, name = "formula")
  if (!is.environment(env)) {
    stop("`env` must be an environment for formula evaluation.", call. = FALSE)
  }
  if (is.null(environment(formula))) environment(formula) <- env
  traits <- utils::modifyList(NMAR_DEFAULT_TRAITS, traits %||% list())

  outcome_expr <- formula[[2L]]
  outcome_vars_all <- unique(all.vars(outcome_expr))
  outcome_is_multi <- nmar_is_multi_outcome_expr(outcome_expr)
  if (length(outcome_vars_all) == 0L) {
    stop("The formula must specify at least one outcome variable on the left-hand side.", call. = FALSE)
  }
  outcome_label <- paste(deparse(outcome_expr, width.cutoff = 500L), collapse = " ")

  rhs_parts <- nmar_partition_rhs(formula[[3L]])
  aux_expr_user <- rhs_parts$aux_expr
  resp_expr <- rhs_parts$response_expr

  validator$assert_data_frame_or_survey(data, name = "data")

  is_survey <- inherits(data, "survey.design")
  data_df <- if (is_survey) data$variables else data
  if (!is.data.frame(data_df)) {
    stop("Unable to access variables from the supplied data object.", call. = FALSE)
  }
  outcome_vars_data <- outcome_vars_all[outcome_vars_all %in% names(data_df)]
  outcome_vars_data <- unique(outcome_vars_data)
  outcome_vars_env <- setdiff(outcome_vars_all, outcome_vars_data)
# Only data-backed outcome symbols participate in validation; environment-only
# helpers are tracked so we can allow constructs such as poly(X, degree).
  if (!length(outcome_vars_data)) {
    stop(
      "Outcome variable",
      if (length(outcome_vars_all) > 1) "s" else "",
      " not found in data: ",
      paste(outcome_vars_all, collapse = ", "),
      call. = FALSE
    )
  }

  outcome_na_flags <- if (length(outcome_vars_data)) {
    vapply(outcome_vars_data, function(var) anyNA(data_df[[var]]), logical(1L))
  } else {
    logical(0)
  }
# Prefer the first outcome column that exhibits missingness so helper-first
# expressions (e.g., I(helper - Y)) still treat Y as the estimand. When all
# outcome columns are fully observed we fall back to the first data symbol.
  outcome_primary <- if (any(outcome_na_flags, na.rm = TRUE)) {
    outcome_vars_data[which(outcome_na_flags)[1L]]
  } else {
    outcome_vars_data[[1L]]
  }
  outcome_helpers <- setdiff(outcome_vars_data, outcome_primary)
  outcome_vars <- if (isTRUE(outcome_is_multi)) outcome_vars_data else outcome_primary

# Normalize auxiliary RHS once traits are known: drop explicit intercepts if
# requested and make sure outcome variables do not leak onto the auxiliary
# side when the user relied on '.' expansion.
  aux_expr <- nmar_normalize_auxiliary_expr(
    aux_expr = aux_expr_user,
    drop_intercept = isTRUE(traits$drop_auxiliary_intercept),
    data = data_df,
    env = environment(formula),
    exclude = outcome_vars_all
  )

  normalized_formula <- formula
  if (!is.null(resp_expr)) {
    normalized_formula[[3L]] <- aux_expr_user
  }

  aux_vars_raw <- unique(setdiff(all.vars(aux_expr_user), "."))
  response_vars_raw <- if (is.null(resp_expr)) {
    character()
  } else {
    unique(setdiff(all.vars(resp_expr), "."))
  }

  blueprint <- nmar_build_formula_blueprint(
    outcome_vars = outcome_vars,
    aux_expr = aux_expr,
    response_expr = resp_expr,
    data = data_df,
    env = environment(formula)
  )

  outcome_names_all <- outcome_vars_data
# Blueprint source variables include whatever symbols model.matrix() kept
# after dot expansion. Remove outcome variables so auxiliaries remain strictly
# covariates.
  aux_source_data <- blueprint$aux$source_variables_data %||% blueprint$aux$source_variables %||% character()
  aux_source_env <- blueprint$aux$source_variables_env %||% character()
  response_source_data <- blueprint$response$source_variables_data %||% blueprint$response$source_variables %||% character()
  response_source_env <- blueprint$response$source_variables_env %||% character()
  auxiliary_vars <- unique(setdiff(aux_source_data, outcome_names_all))
  auxiliary_vars_env <- unique(setdiff(aux_source_env, outcome_names_all))
  response_predictors <- unique(response_source_data)
  response_predictors_env <- unique(response_source_env)

  structure(
    list(
      formula = normalized_formula,
      formula_original = formula,
      outcome = outcome_vars,
      outcome_primary = outcome_primary,
      outcome_helpers = outcome_helpers,
      outcome_vars_data = outcome_vars_data,
      outcome_vars_env = outcome_vars_env,
      outcome_expr = outcome_expr,
      outcome_label = outcome_label,
      outcome_is_multi = outcome_is_multi,
      auxiliary_vars = auxiliary_vars,
      auxiliary_vars_env = auxiliary_vars_env,
      auxiliary_vars_raw = aux_vars_raw,
      response_predictors = response_predictors,
      response_predictors_env = response_predictors_env,
      response_predictors_raw = response_vars_raw,
      aux_rhs_lang = aux_expr,
      aux_rhs_lang_user = aux_expr_user,
      response_rhs_lang = resp_expr,
      response_rhs_lang_user = resp_expr,
      data = data_df,
      original_data = data,
      is_survey = is_survey,
      environment = environment(formula),
      blueprint = blueprint
    ),
    class = "nmar_input_spec"
  )
}

#' Create a normalized NMAR task object
#' @keywords internal
new_nmar_task <- function(spec, traits, trace_level = 1) {
  if (!inherits(spec, "nmar_input_spec")) {
    stop("`spec` must be created by `parse_nmar_spec()`.", call. = FALSE)
  }
  structure(
    list(
      formula = spec$formula,
      formula_original = spec$formula_original,
      outcome = spec$outcome,
      outcome_primary = spec$outcome_primary,
      outcome_helpers = spec$outcome_helpers,
      outcome_vars_data = spec$outcome_vars_data,
      outcome_vars_env = spec$outcome_vars_env,
      outcome_expr = spec$outcome_expr,
      outcome_label = spec$outcome_label,
      outcome_is_multi = spec$outcome_is_multi,
      auxiliary_vars = spec$auxiliary_vars,
      auxiliary_vars_env = spec$auxiliary_vars_env,
      auxiliary_vars_raw = spec$auxiliary_vars_raw,
      response_predictors = spec$response_predictors,
      response_predictors_env = spec$response_predictors_env,
      response_predictors_raw = spec$response_predictors_raw,
      aux_rhs_lang = spec$aux_rhs_lang,
      aux_rhs_lang_user = spec$aux_rhs_lang_user,
      response_rhs_lang = spec$response_rhs_lang,
      response_rhs_lang_user = spec$response_rhs_lang_user,
      data = spec$data,
      original_data = spec$original_data,
      is_survey = spec$is_survey,
      environment = spec$environment,
      blueprint = spec$blueprint,
      requires_single_outcome = isTRUE(traits$requires_single_outcome),
      trace_level = trace_level
    ),
    class = "nmar_task"
  )
}

#' Prepare common NMAR design components
#' @keywords internal
prepare_nmar_design <- function(task,
                                standardize = TRUE,
                                auxiliary_means = NULL,
    include_response = TRUE,
    include_auxiliary = TRUE,
    data = task$data,
                                design_weights = NULL) {
  if (!inherits(task, "nmar_task")) {
    stop("`task` must be created by `new_nmar_task()`.", call. = FALSE)
  }
  data_df <- data
  if (!is.data.frame(data_df)) {
    stop("`data` must be a data.frame when preparing the NMAR design.", call. = FALSE)
  }

  outcome_info <- nmar_prepare_outcome_column(
    data = data_df,
    expr = task$outcome_expr,
    env = task$environment
  )
  data_df <- outcome_info$data
  outcome_columns <- outcome_info$columns
  outcome_column <- outcome_columns[[1L]]
  if (isTRUE(task$requires_single_outcome) && length(outcome_columns) != 1L) {
    stop("The formula must have exactly one outcome variable on the left-hand side.", call. = FALSE)
  }

# Survey designs keep their own copy of the data inside design$variables. When
# LHS transformations create derived columns we copy them back so that
# survey::weights() and friends can find the variables referenced in the
# formula. Only the new columns are inserted; user data is left untouched.
  if (isTRUE(task$is_survey) && !is.null(task$original_data) && length(outcome_info$added_columns)) {
    for (col in outcome_info$added_columns) {
      task$original_data$variables[[col]] <- data_df[[col]]
    }
  }

  weights <- design_weights
  if (is.null(weights)) {
    if (isTRUE(task$is_survey)) {
      weights <- stats::weights(task$original_data)
      if (is.null(weights)) {
        stop("Unable to retrieve design weights from the supplied survey design.", call. = FALSE)
      }
    } else {
      weights <- rep(1, nrow(data_df))
    }
  }
  weights <- as.numeric(weights)
  if (length(weights) != nrow(data_df)) {
    stop("`design_weights` must have the same length as the supplied data.", call. = FALSE)
  }

  engine_formula <- nmar_replace_formula_lhs(task$formula, outcome_columns)
  user_base_formula <- task$formula_original %||% task$formula
  user_formula <- nmar_rebuild_partitioned_formula(
    base_formula = user_base_formula,
    response_rhs_lang = task$response_rhs_lang_user %||% task$response_rhs_lang,
    aux_rhs_lang = task$aux_rhs_lang_user %||% task$aux_rhs_lang,
    env = task$environment
  )

  design_matrices <- nmar_materialize_design_matrices(
    blueprint = task$blueprint,
    data = data_df,
    include_auxiliary = include_auxiliary,
    include_response = include_response
  )
  list(
    data = data_df,
    outcome = task$outcome,
    auxiliary_vars = if (isTRUE(include_auxiliary)) task$auxiliary_vars else character(),
    response_predictors = if (isTRUE(include_response)) task$response_predictors else character(),
    weights = weights,
    is_survey = isTRUE(task$is_survey),
    survey_design = if (isTRUE(task$is_survey)) task$original_data else NULL,
    formula = task$formula,
    user_formula = user_formula,
    engine_formula = engine_formula,
    standardize = standardize,
    auxiliary_means = auxiliary_means,
    aux_rhs_lang = if (isTRUE(include_auxiliary)) task$aux_rhs_lang else NULL,
    response_rhs_lang = if (isTRUE(include_response)) task$response_rhs_lang else NULL,
    blueprint = task$blueprint,
    outcome_column = outcome_column,
    outcome_columns = outcome_columns,
    outcome_label = task$outcome_label,
    design_matrices = design_matrices,
    environment = task$environment
  )
}

#' @keywords internal
nmar_materialize_design_matrices <- function(blueprint,
                                             data,
                                             include_auxiliary = TRUE,
                                             include_response = TRUE) {
  aux_matrix <- NULL
  response_matrix <- NULL
  if (isTRUE(include_auxiliary) && !is.null(blueprint$aux$terms)) {
    aux_matrix <- nmar_model_matrix_from_terms(blueprint$aux, data)
  }
  if (isTRUE(include_response) && !is.null(blueprint$response$terms)) {
    response_matrix <- nmar_model_matrix_from_terms(blueprint$response, data)
  }
  list(auxiliary = aux_matrix, response = response_matrix)
}

#' @keywords internal
nmar_prepare_outcome_column <- function(data, expr, env) {
  if (!is.call(expr) && is.symbol(expr)) {
    column <- as.character(expr)
    if (!column %in% names(data)) {
      stop("Outcome variable '", column, "' not found in data.", call. = FALSE)
    }
    return(list(
      data = data,
      columns = column,
      added_columns = character()
    ))
  }

  outcome_formula <- stats::as.formula(call("~", expr))
  attr(outcome_formula, ".Environment") <- env
  outcome_frame <- stats::model.frame(
    outcome_formula,
    data = data,
    na.action = stats::na.pass,
    drop.unused.levels = FALSE
  )
  outcome_values <- outcome_frame[[1]]
# Track rows where any input variable in the LHS expression is missing; we
# do not penalize non-finite results that arise solely from missing inputs.
  input_missing <- Reduce(
    "|",
    lapply(
      all.vars(expr),
      function(var) if (var %in% names(data)) is.na(data[[var]]) else rep(FALSE, nrow(data))
    ),
    init = rep(FALSE, nrow(data))
  )

  coerce_numeric_column <- function(values, label) {
    if (!is.numeric(values)) {
      stop("Outcome transformation must evaluate to numeric vector(s). Offending component: ", label, call. = FALSE)
    }
    bad_idx <- which(!input_missing & !is.finite(values))
    if (length(bad_idx)) {
      stop(
        "Outcome transformation produced non-finite values in component '", label, "'. ",
        "First offending row: ", bad_idx[1],
        call. = FALSE
      )
    }
    values
  }

  append_column <- function(col_values, base_name, data_env) {
    if (!is.null(base_name) &&
        base_name %in% names(data_env) &&
        identical(col_values, data_env[[base_name]])) {
      return(list(name = base_name, data = data_env, added = FALSE))
    }
    new_name <- base_name %||% "..nmar_outcome.."
    if (new_name %in% names(data_env)) {
      new_name <- nmar_make_unique_colname(new_name, names(data_env))
    }
    data_env[[new_name]] <- col_values
    list(name = new_name, data = data_env, added = TRUE)
  }

  if (!is.data.frame(outcome_values) && !(is.matrix(outcome_values) && ncol(outcome_values) > 1L)) {
    if (is.matrix(outcome_values)) {
      outcome_values <- outcome_values[, 1L, drop = TRUE]
    }
    outcome_values <- coerce_numeric_column(outcome_values, deparse(expr, width.cutoff = 60))
    new_column <- nmar_make_unique_colname("..nmar_outcome..", names(data))
    data[[new_column]] <- outcome_values
    return(list(
      data = data,
      columns = new_column,
      added_columns = new_column
    ))
  }

  if (is.matrix(outcome_values)) {
    outcome_df <- as.data.frame(outcome_values, stringsAsFactors = FALSE)
  } else {
    outcome_df <- outcome_values
  }

  if (!ncol(outcome_df)) {
    stop("Outcome transformation must produce at least one column.", call. = FALSE)
  }

  new_columns <- character(ncol(outcome_df))
  added_columns <- logical(ncol(outcome_df))
  for (idx in seq_len(ncol(outcome_df))) {
    col_label <- colnames(outcome_df)[idx] %||% paste0("..nmar_outcome..", idx)
    col_values <- coerce_numeric_column(outcome_df[[idx]], col_label)
    insertion <- append_column(col_values, col_label, data)
    data <- insertion$data
    new_columns[[idx]] <- insertion$name
    added_columns[[idx]] <- isTRUE(insertion$added)
  }

  list(
    data = data,
    columns = new_columns,
    added_columns = new_columns[added_columns]
  )
}

#' @keywords internal
nmar_replace_formula_lhs <- function(formula, outcome_columns) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` must be a two-sided formula.", call. = FALSE)
  }
  new_formula <- formula
  if (!length(outcome_columns)) {
    stop("`outcome_columns` must contain at least one column name.", call. = FALSE)
  }
  lhs_expr <- if (length(outcome_columns) == 1L) {
    as.name(outcome_columns[[1L]])
  } else {
    as.call(c(as.name("cbind"), lapply(outcome_columns, as.name)))
  }
  new_formula[[2L]] <- lhs_expr
  attr(new_formula, ".Environment") <- environment(formula)
  new_formula
}
