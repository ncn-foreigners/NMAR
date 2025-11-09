#' NMAR input parsing and validation helpers
#'
#' These developer-facing utilities centralize the parsing of user supplied
#' formulas/data and allow engines to declare small trait lists that control
#' validation strictness. Keeping all argument checks in one place ensures
#' consistent error messages and makes it easier to extend the package with
#' additional engines.
#'
#' @keywords internal

#' Default engine trait values
#'
#' Canonical default values for engine validation traits. Engines can override
#' these via their engine_traits methods. Defined as a constant to ensure
#' consistency across the package.
#'
#' @keywords internal
NMAR_DEFAULT_TRAITS <- list(
  allow_outcome_in_missingness = FALSE,
  allow_covariate_overlap = FALSE,
  requires_single_outcome = TRUE,
  allow_respondents_only = FALSE
)

#' Extract auxiliary and response expressions from partitioned RHS
#'
#' Splits a formula RHS containing `|` into auxiliary (left of |) and response
#' (right of |) components. If no `|` is present, returns the entire RHS as
#' auxiliary with NULL response.
#'
#' @param rhs_expr Right-hand side language object (may contain `|`)
#' @return Named list with `aux_expr` and `response_expr` (NULL if no `|`)
#' @keywords internal
nmar_partition_rhs <- function(rhs_expr) {
  aux_expr <- rhs_expr
  response_expr <- NULL
  if (is.call(rhs_expr) && identical(rhs_expr[[1L]], as.name("|"))) {
    aux_expr <- rhs_expr[[2L]]
    response_expr <- rhs_expr[[3L]]
  }
  list(aux_expr = aux_expr, response_expr = response_expr)
}

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

# Support partitioned RHS using `|`: y ~ aux | response
  rhs_parts <- nmar_partition_rhs(formula[[3L]])
  aux_expr <- rhs_parts$aux_expr
  resp_expr <- rhs_parts$response_expr

  validator$assert_data_frame_or_survey(data, name = "data")

  is_survey <- inherits(data, "survey.design")
  data_df <- if (is_survey) data$variables else data
  if (!is.data.frame(data_df)) {
    stop("Unable to access variables from the supplied data object.", call. = FALSE)
  }

# Rebuild a normalized formula whose RHS is just auxiliaries (the response side is
# tracked separately as `response_predictors`). This ensures downstream code that
# prints or stores the user formula keeps the original environment.
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

#' Engine trait declarations
#'
#' Public S3 generic returning a small list of declarative flags that describe
#' how an engine expects input validation to behave. This function is also used
#' internally by the input pipeline. Users can call it on an engine object to
#' introspect behaviour such as whether outcome variables may appear on the
#' response RHS explicitly or whether multiple outcomes are supported.
#'
#' Important: The NMAR response model always includes the outcome implicitly by
#' design. The `allow_outcome_in_missingness` trait governs whether the outcome
#' may be specified explicitly on the response RHS (e.g., `Y ~ aux | Y + Z`).
#'
#' @param engine An object inheriting from class `nmar_engine`.
#' @return A named list of trait flags. The current fields are:
#'   - `allow_outcome_in_missingness`: logical.
#'   - `allow_covariate_overlap`: logical.
#'   - `requires_single_outcome`: logical.
#'   Engines may add traits over time; callers should use `$` with care
#'   and rely on presence checks when needed.
#' @keywords engine_view
#' @export
engine_traits <- function(engine) {
  UseMethod("engine_traits")
}
#' @keywords engine_view
#' @export
engine_traits.default <- function(engine) {
  NMAR_DEFAULT_TRAITS
}
#' @keywords engine_view
#' @export
engine_traits.nmar_engine <- function(engine) {
# Parent class falls back to the baseline defaults; specific engines
# override in their class-specific methods.
  engine_traits.default(engine)
}

#' @keywords engine_view
#' @export
engine_traits.nmar_engine_el <- function(engine) {
# Activate respondents-only mode only if n_total is provided by the engine
  allow_resp_only <- !is.null(engine$n_total)

  utils::modifyList(
    engine_traits.default(engine),
    list(
      allow_outcome_in_missingness = TRUE,
      allow_covariate_overlap = TRUE,
      allow_respondents_only = allow_resp_only
    )
  )
}

#' @export
engine_traits.nmar_engine_exptilt <- function(engine) {
  engine_traits.default(engine)
}

#' @export
engine_traits.nmar_engine_exptilt_nonparam <- function(engine) {
  utils::modifyList(
    engine_traits.default(engine),
    list(requires_single_outcome = FALSE)
  )
}

#' Validate parsed NMAR inputs
#'
#' @param spec Object produced by [parse_nmar_spec()].
#' @param traits List of engine traits produced by [engine_traits()].
#'
#' @keywords internal
validate_nmar_args <- function(spec, traits = list()) {
  if (!inherits(spec, "nmar_input_spec")) {
    stop("`spec` must be created by `parse_nmar_spec()`.", call. = FALSE)
  }
  traits <- utils::modifyList(NMAR_DEFAULT_TRAITS, traits)

  if (traits$requires_single_outcome && length(spec$outcome) != 1L) {
    stop("The formula must have exactly one outcome variable on the left-hand side.", call. = FALSE)
  }

  aux_vars_check <- unique(c(spec$auxiliary_vars, spec$auxiliary_vars_raw))
# Enforce engine policy using only the explicit response RHS variables.
# The outcome is always included implicitly in the response model by design;
# the allow_outcome_in_missingness trait governs explicit usage only.
  response_vars_check <- unique(spec$response_predictors_raw)

  validate_predictor_relationships(
    outcomes = spec$outcome,
    auxiliary_vars = aux_vars_check,
    response_vars = response_vars_check,
    allow_outcome_in_missingness = traits$allow_outcome_in_missingness,
    allow_covariate_overlap = traits$allow_covariate_overlap
  )

  if (length(spec$outcome) == 1L) {
    validate_data(
      data = spec$original_data,
      outcome_variable = spec$outcome[[1]],
      covariates_for_outcome = aux_vars_check,
      covariates_for_missingness = response_vars_check,
      allow_outcome_in_missingness = traits$allow_outcome_in_missingness,
      allow_covariate_overlap = traits$allow_covariate_overlap,
      allow_respondents_only = isTRUE(traits$allow_respondents_only)
    )
  } else {
    validate_multi_outcome_data(spec$data, spec$outcome)
    nmar_validate_covariates(spec$data, aux_vars_check, block_label = "auxiliary")
    nmar_validate_covariates(spec$data, response_vars_check, block_label = "response")
  }

  invisible(spec)
}

#' Create a normalized NMAR task object
#'
#' @param spec Parsed NMAR specification produced by [parse_nmar_spec()].
#' @param traits Engine traits returned by [engine_traits()].
#' @param trace_level Integer 0-3 controlling verbosity during estimation.
#'
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
#'
#' @param task An object created by [new_nmar_task()].
#' @param standardize Logical flag forwarded to engines.
#' @param auxiliary_means Optional named vector of auxiliary means.
#' @param include_response Logical; include response predictors.
#' @param include_auxiliary Logical; include auxiliary predictors.
#' @param data Optional override of the data-frame backing the task.
#' @param design_weights Optional numeric vector of design weights.
#'
#' @return A list containing the trimmed data, predictor sets, weights, survey
#'   design (if applicable), formula, and standardization settings.
#'
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

#' Rebuild a partitioned formula y ~ aux | response
#'
#' Given a base formula whose RHS contains auxiliaries only (the normalized
#' `task$formula`) and language objects for the auxiliary and response partitions,
#' construct a new partitioned formula `y ~ aux | response` without string
#' manipulation. If `response_rhs_lang` is `NULL`, returns an unpartitioned
#' `y ~ aux` formula. The formula environment is preserved.
#'
#' @param base_formula A two-sided formula (y ~ aux-only) whose environment is set.
#' @param response_rhs_lang An R language object for the response-only predictors (right of `|`), or `NULL`.
#' @param aux_rhs_lang An R language object for the auxiliary predictors (left of `|`), or `NULL` to reuse the base RHS.
#' @param env Optional formula environment override; defaults to the base formula environment.
#' @keywords internal
nmar_rebuild_partitioned_formula <- function(base_formula,
                                            response_rhs_lang = NULL,
                                            aux_rhs_lang = NULL,
                                            env = NULL) {
  if (!inherits(base_formula, "formula") || length(base_formula) != 3L) {
    stop("`base_formula` must be a two-sided formula.", call. = FALSE)
  }
  lhs <- base_formula[[2L]]
  rhs_aux <- aux_rhs_lang %||% base_formula[[3L]]
# If no response language provided, return base or aux-adjusted formula
  if (is.null(response_rhs_lang)) {
    f <- call("~", lhs, rhs_aux)
  } else {
    rhs_bar <- call("|", rhs_aux, response_rhs_lang)
    f <- call("~", lhs, rhs_bar)
  }
  if (is.null(env)) env <- environment(base_formula)
  if (is.null(env)) env <- parent.frame()
  class(f) <- "formula"
  attr(f, ".Environment") <- env
  f
}

# Split a partitioned formula y ~ aux | response into language parts
# Returns list(outcome_var, aux_rhs_lang, response_rhs_lang, env)
nmar_split_partitioned_formula <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` must be a two-sided formula.", call. = FALSE)
  }
  env <- environment(formula)
  if (is.null(env)) env <- parent.frame()
  outcome_var <- all.vars(formula[[2L]])
  if (length(outcome_var) != 1L) stop("LHS must be a single outcome variable.", call. = FALSE)

  rhs_parts <- nmar_partition_rhs(formula[[3L]])
  list(
    outcome_var = outcome_var[[1L]],
    aux_rhs_lang = rhs_parts$aux_expr,
    response_rhs_lang = rhs_parts$response_expr,
    env = env
  )
}

# Build internal EL formulas (outcome, response, auxiliary) using language objects
# delta_name must already be unique in the provided data
nmar_build_internal_formulas <- function(delta_name, outcome_var, aux_rhs_lang, response_rhs_lang, env) {
# outcome: y ~ 1
  outcome_fml <- stats::as.formula(call("~", as.name(outcome_var), 1L))
  attr(outcome_fml, ".Environment") <- env
# response: delta ~ outcome + response_rhs_lang (or ~ outcome if NULL)
  rhs_resp <- if (is.null(response_rhs_lang)) {
    as.name(outcome_var)
  } else {
    call("+", as.name(outcome_var), response_rhs_lang)
  }
  response_fml <- stats::as.formula(call("~", as.name(delta_name), rhs_resp))
  attr(response_fml, ".Environment") <- env
# auxiliary: ~ 0 + aux_rhs_lang (or NULL)
  auxiliary_fml <- NULL
  if (!is.null(aux_rhs_lang)) {
    rhs_aux <- call("+", 0, aux_rhs_lang) # 0 + x removes intercept
    auxiliary_fml <- stats::as.formula(call("~", rhs_aux))
    attr(auxiliary_fml, ".Environment") <- env
  }
  list(outcome = outcome_fml, response = response_fml, auxiliary = auxiliary_fml)
}

# Generate a unique column name that does not collide with names(data)
nmar_make_unique_colname <- function(base, data_names) {
  nm <- base
  if (!(nm %in% data_names)) return(nm)
  i <- 1L
  repeat {
    cand <- paste0(base, i)
    if (!(cand %in% data_names)) return(cand)
    i <- i + 1L
  }
}

# Build formula terms/metadata so engines can reuse resolved model matrices
nmar_build_formula_blueprint <- function(outcome_vars,
                                         aux_expr,
                                         response_expr,
                                         data,
                                         env) {
  aux_formula <- nmar_make_aux_formula(outcome_vars, aux_expr, env)
  aux_terms <- nmar_build_terms_info(aux_formula, data, drop_response = TRUE)

  response_formula <- nmar_make_response_formula(outcome_vars, response_expr, env)
  response_info <- nmar_build_terms_info(response_formula, data, drop_response = TRUE)

  list(
    aux = aux_terms,
    response = response_info,
    outcome = outcome_vars
  )
}

nmar_make_aux_formula <- function(outcome_vars, aux_expr, env) {
  if (is.null(aux_expr)) return(NULL)
  if (length(all.vars(aux_expr)) == 0L) return(NULL)
  lhs <- nmar_make_lhs_expr(outcome_vars)
  rhs <- call("+", 0, aux_expr)
  if (is.null(lhs)) {
    f <- call("~", rhs)
  } else {
    f <- call("~", lhs, rhs)
  }
  formula <- stats::as.formula(f, env = env)
  formula
}

nmar_make_response_formula <- function(outcome_vars, response_expr, env) {
  if (is.null(response_expr)) return(NULL)
  if (length(all.vars(response_expr)) == 0L) return(NULL)
  lhs <- nmar_make_lhs_expr(outcome_vars)
  rhs <- call("+", 0, response_expr)
  if (is.null(lhs)) {
    f <- call("~", rhs)
  } else {
    f <- call("~", lhs, rhs)
  }
  stats::as.formula(f, env = env)
}

nmar_build_terms_info <- function(formula, data, drop_response = FALSE) {
  if (is.null(formula)) {
    return(list(
      formula = NULL,
      terms = NULL,
      column_names = character(),
      source_variables = character(),
      xlevels = NULL,
      contrasts = NULL
    ))
  }
  vars_needed <- unique(setdiff(all.vars(formula), "."))
  missing_vars <- setdiff(vars_needed, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass, drop.unused.levels = FALSE)
  tr <- attr(mf, "terms")
  tr_use <- tr
  if (isTRUE(drop_response) && attr(tr, "response") > 0) {
    tr_use <- stats::delete.response(tr)
  }
  mm <- stats::model.matrix(tr_use, mf)
# Derive source variable names from the RHS terms only (exclude response)
# by converting the terms object back to a response-less formula.
  rhs_vars <- unique(setdiff(all.vars(stats::formula(tr_use)), "."))
  list(
    formula = formula,
    terms = tr_use,
    column_names = colnames(mm),
    source_variables = rhs_vars,
    xlevels = attr(mf, "xlevels"),
    contrasts = attr(mm, "contrasts")
  )
}

nmar_make_lhs_expr <- function(outcome_vars) {
  if (length(outcome_vars) == 0L) return(NULL)
  lhs <- as.name(outcome_vars[[1L]])
  if (length(outcome_vars) == 1L) return(lhs)
  for (var in outcome_vars[-1L]) {
    lhs <- call("+", lhs, as.name(var))
  }
  lhs
}

validate_multi_outcome_data <- function(data, outcome_vars) {
  missing_outcomes <- setdiff(outcome_vars, names(data))
  if (length(missing_outcomes) > 0) {
    stop("Outcome variables not found in data: ", paste(missing_outcomes, collapse = ", "))
  }
  for (outcome_var in outcome_vars) {
    col <- data[[outcome_var]]
    if (!is.numeric(col)) {
      bad_val <- col[which(!is.numeric(col))[1]]
      stop(
        "Outcome variable '", outcome_var, "' must be numeric.\n",
        "First invalid value: '", bad_val, "' at row ", which(!is.numeric(col))[1]
      )
    }
    if (all(is.na(col))) {
      stop(
        "Outcome variable '", outcome_var, "' cannot be entirely NA.",
        call. = FALSE
      )
    }
  }
}

nmar_model_matrix_from_terms <- function(terms_info, data) {
  if (is.null(terms_info$terms)) {
    return(NULL)
  }
  mf <- stats::model.frame(
    terms_info$terms,
    data = data,
    na.action = stats::na.pass,
    drop.unused.levels = FALSE,
    xlev = terms_info$xlevels
  )
  stats::model.matrix(terms_info$terms, mf, contrasts.arg = terms_info$contrasts)
}

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

# Handle vector output
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

# Handle multi-column output (matrix/data.frame)
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

validate_predictor_relationships <- function(outcomes,
                                             auxiliary_vars,
                                             response_vars,
                                             allow_outcome_in_missingness,
                                             allow_covariate_overlap) {
  auxiliary_vars <- unique(auxiliary_vars %||% character())
  response_vars <- unique(response_vars %||% character())
  if (!allow_outcome_in_missingness) {
    outcome_on_rhs <- intersect(outcomes, response_vars)
    if (length(outcome_on_rhs) > 0) {
      stop(
        "Outcome variable cannot appear on the response-model side unless the engine allows it.\n",
        "Variables: ", paste(outcome_on_rhs, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (!allow_covariate_overlap) {
    overlap <- intersect(auxiliary_vars, response_vars)
    if (length(overlap) > 0) {
      stop(
        "Covariate sets must be mutually exclusive. Overlapping variables: ",
        paste(overlap, collapse = ", "),
        call. = FALSE
      )
    }
  }
}

nmar_validate_covariates <- function(data, vars, block_label) {
  vars <- unique(vars)
  if (!length(vars)) return(invisible(NULL))
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0) {
    stop(
      block_label, " covariates not found in data: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }
  for (var in vars) {
    column <- data[[var]]
    valid_type <- is.numeric(column) || is.logical(column) || is.factor(column) || is.character(column)
    if (!valid_type) {
      stop(
        "Covariate '", var, "' in ", block_label, " predictors must be numeric, logical, factor, or character.\n",
        "Detected type: ", paste(class(column), collapse = "/"),
        call. = FALSE
      )
    }
    if (is.numeric(column)) {
      bad_indices <- which(!is.na(column) & !is.finite(column))
      if (length(bad_indices) > 0) {
        first_bad_idx <- bad_indices[1]
        bad_val <- column[first_bad_idx]
        stop(
          block_label, " covariate '", var, "' contains non-finite values (Inf, -Inf, or NaN).\n",
          "First non-finite value: ", bad_val, " at row ", first_bad_idx,
          call. = FALSE
        )
      }
    }
    if (anyNA(column)) {
      stop(
        block_label, " covariate '", var, "' contains NA values.\n",
        "First NA at row ", which(is.na(column))[1],
        call. = FALSE
      )
    }
  }
}
