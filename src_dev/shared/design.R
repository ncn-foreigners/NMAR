#' Input parsing and design preparation
#'
#' parse_nmar_spec() turns a user formula and data into a normalized
#' nmar_input_spec with resolved RHS partitions, preserved environments,
#' and a blueprint of terms for fast model.matrix materialization.
#'
#' prepare_nmar_design() uses the spec and traits to construct a
#' method-agnostic design_info bundle with:
#' - data with derived outcome column(s) (LHS transforms are allowed)
#' - consistent weights (from survey.design or ones for data.frames)
#' - design matrices for auxiliaries and explicit response RHS (no intercepts)
#' - user and engine formulas (partitioned and LHS-rewritten)
#'
#' Non-finite values produced by LHS transformations are allowed only on rows
#' where the input variables are missing (nonrespondents); for observed rows we
#' error early with a clear message. This prevents subtle NA/Inf propagation.
#'
#' @keywords internal

parse_nmar_spec <- function(formula, data, env = parent.frame()) {
  validator$assert_formula_two_sided(formula, name = "formula")
  if (!is.environment(env)) {
    stop("`env` must be an environment for formula evaluation.", call. = FALSE)
  }
  if (is.null(environment(formula))) environment(formula) <- env

  outcome_expr <- formula[[2L]]
  outcome_vars <- unique(all.vars(outcome_expr))
  if (length(outcome_vars) == 0L) {
    stop("The formula must specify at least one outcome variable on the left-hand side.", call. = FALSE)
  }
  outcome_label <- paste(deparse(outcome_expr, width.cutoff = 500L), collapse = " ")

  rhs_parts <- nmar_partition_rhs(formula[[3L]])
  aux_expr <- rhs_parts$aux_expr
  resp_expr <- rhs_parts$response_expr

  validator$assert_data_frame_or_survey(data, name = "data")

  is_survey <- inherits(data, "survey.design")
  data_df <- if (is_survey) data$variables else data
  if (!is.data.frame(data_df)) {
    stop("Unable to access variables from the supplied data object.", call. = FALSE)
  }

  normalized_formula <- formula
  if (!is.null(resp_expr)) {
    normalized_formula[[3L]] <- aux_expr
  }

  aux_vars_raw <- unique(setdiff(all.vars(aux_expr), "."))
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

  auxiliary_vars <- blueprint$aux$source_variables
  response_predictors <- blueprint$response$source_variables

  structure(
    list(
      formula = normalized_formula,
      outcome = outcome_vars,
      outcome_expr = outcome_expr,
      outcome_label = outcome_label,
      auxiliary_vars = auxiliary_vars,
      auxiliary_vars_raw = aux_vars_raw,
      response_predictors = response_predictors,
      response_predictors_raw = response_vars_raw,
      aux_rhs_lang = aux_expr,
      response_rhs_lang = resp_expr,
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
      outcome = spec$outcome,
      outcome_expr = spec$outcome_expr,
      outcome_label = spec$outcome_label,
      auxiliary_vars = spec$auxiliary_vars,
      auxiliary_vars_raw = spec$auxiliary_vars_raw,
      response_predictors = spec$response_predictors,
      response_predictors_raw = spec$response_predictors_raw,
      aux_rhs_lang = spec$aux_rhs_lang,
      response_rhs_lang = spec$response_rhs_lang,
      data = spec$data,
      original_data = spec$original_data,
      is_survey = spec$is_survey,
      environment = spec$environment,
      traits = traits,
      blueprint = spec$blueprint,
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
  user_formula <- nmar_rebuild_partitioned_formula(
    base_formula = task$formula,
    response_rhs_lang = task$response_rhs_lang,
    aux_rhs_lang = task$aux_rhs_lang,
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
    aux_terms = if (isTRUE(include_auxiliary)) task$blueprint$aux$terms else NULL,
    response_terms = if (isTRUE(include_response)) task$blueprint$response$terms else NULL,
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
