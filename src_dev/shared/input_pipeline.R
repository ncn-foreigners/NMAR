#' NMAR input parsing and validation helpers
#'
#' These developer-facing utilities centralize the parsing of user supplied
#' formulas/data and allow engines to declare small trait lists that control
#' validation strictness. Keeping all argument checks in one place ensures
#' consistent error messages and makes it easier to extend the package with
#' additional engines.
#'
#' @keywords internal
parse_nmar_spec <- function(formula, data, response_predictors = NULL, env = parent.frame()) {
  if (!inherits(formula, "formula") || length(formula) != 3) {
    stop("`formula` must be a two-sided formula like y ~ x1 + x2.", call. = FALSE)
  }
  if (!is.environment(env)) env <- parent.frame()
  if (is.null(environment(formula))) environment(formula) <- env

  outcome_vars <- unique(all.vars(formula[[2]]))
  if (length(outcome_vars) == 0L) {
    stop("The formula must specify at least one outcome variable on the left-hand side.", call. = FALSE)
  }
  auxiliary_vars <- unique(all.vars(formula[[3]]))

  if (is.null(response_predictors)) {
    response_predictors <- character()
  }
  if (!is.character(response_predictors)) {
    stop("`response_predictors` must be a character vector of variable names or NULL.", call. = FALSE)
  }
  if (length(response_predictors) && anyNA(response_predictors)) {
    stop("`response_predictors` cannot contain NA values.", call. = FALSE)
  }
  response_predictors <- unique(response_predictors)

  if (!inherits(data, c("data.frame", "survey.design"))) {
    stop("`data` must be a data.frame or survey.design object.", call. = FALSE)
  }

  is_survey <- inherits(data, "survey.design")
  data_df <- if (is_survey) data$variables else data
  if (!is.data.frame(data_df)) {
    stop("Unable to access variables from the supplied data object.", call. = FALSE)
  }

  structure(
    list(
      formula = formula,
      outcome = outcome_vars,
      auxiliary_vars = auxiliary_vars,
      response_predictors = response_predictors,
      data = data_df,
      original_data = data,
      is_survey = is_survey,
      environment = environment(formula)
    ),
    class = "nmar_input_spec"
  )
}

#' Engine trait declarations
#'
#' Public S3 generic returning a small list of declarative flags that describe
#' how an engine expects input validation to behave. This function is also used
#' internally by the input pipeline. Users can call it on an engine object to
#' introspect behaviour such as whether outcome variables may appear in the
#' missingness model or whether multiple outcomes are supported.
#'
#' @param engine An object inheriting from class `nmar_engine`.
#' @return A named list of trait flags. The current fields are:
#'   - `allow_outcome_in_missingness`: logical.
#'   - `allow_covariate_overlap`: logical.
#'   - `requires_single_outcome`: logical.
#'   Engines may add traits over time; callers should use `$` with care
#'   and rely on presence checks when needed.
#' @export
engine_traits <- function(engine) {
  UseMethod("engine_traits")
}

#' @export
engine_traits.default <- function(engine) {
  list(
    allow_outcome_in_missingness = FALSE,
    allow_covariate_overlap = FALSE,
    requires_single_outcome = TRUE
  )
}

#' @export
engine_traits.nmar_engine <- function(engine) {
# Parent class falls back to the baseline defaults; specific engines
# override in their class-specific methods.
  engine_traits.default(engine)
}

#' @export
engine_traits.nmar_engine_el <- function(engine) {
  utils::modifyList(
    engine_traits.default(engine),
    list(
      allow_outcome_in_missingness = TRUE,
      allow_covariate_overlap = TRUE
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
  default_traits <- engine_traits(NULL)
  traits <- utils::modifyList(default_traits, traits)

  if (traits$requires_single_outcome && length(spec$outcome) != 1L) {
    stop("The formula must have exactly one outcome variable on the left-hand side.", call. = FALSE)
  }

  if (length(spec$outcome) == 1L) {
    validate_data(
      data = spec$original_data,
      outcome_variable = spec$outcome[[1]],
      covariates_for_outcome = spec$auxiliary_vars,
      covariates_for_missingness = spec$response_predictors,
      allow_outcome_in_missingness = traits$allow_outcome_in_missingness,
      allow_covariate_overlap = traits$allow_covariate_overlap
    )
  } else {
    missing_outcomes <- setdiff(spec$outcome, names(spec$data))
    if (length(missing_outcomes) > 0) {
      stop("Outcome variables not found in data: ", paste(missing_outcomes, collapse = ", "))
    }
    for (outcome_var in spec$outcome) {
      col <- spec$data[[outcome_var]]
      if (!is.numeric(col)) {
        bad_val <- col[which(!is.numeric(col))[1]]
        stop(
          "Outcome variable '", outcome_var, "' must be numeric.\n",
          "First invalid value: '", bad_val, "' at row ", which(!is.numeric(col))[1]
        )
      }
      if (anyNA(col)) {
        stop(
          "Outcome variable '", outcome_var, "' contains NA values.\n",
          "First NA at row ", which(is.na(col))[1]
        )
      }
    }
  }

  invisible(spec)
}

#' Create a normalized NMAR task object
#'
#' @param spec Parsed NMAR specification produced by [parse_nmar_spec()].
#' @param traits Engine traits returned by [engine_traits()].
#'
#' @keywords internal
new_nmar_task <- function(spec, traits) {
  if (!inherits(spec, "nmar_input_spec")) {
    stop("`spec` must be created by `parse_nmar_spec()`.", call. = FALSE)
  }
  structure(
    list(
      formula = spec$formula,
      outcome = spec$outcome,
      auxiliary_vars = spec$auxiliary_vars,
      response_predictors = spec$response_predictors,
      data = spec$data,
      original_data = spec$original_data,
      is_survey = spec$is_survey,
      environment = spec$environment,
      traits = traits
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
  list(
    data = data_df,
    outcome = task$outcome,
    auxiliary_vars = if (isTRUE(include_auxiliary)) task$auxiliary_vars else character(),
    response_predictors = if (isTRUE(include_response)) task$response_predictors else character(),
    weights = as.numeric(weights),
    is_survey = isTRUE(task$is_survey),
    survey_design = if (isTRUE(task$is_survey)) task$original_data else NULL,
    formula = task$formula,
    standardize = standardize,
    auxiliary_means = auxiliary_means
  )
}
