# R/validation_general.R (lub gdziekolwiek chcesz umieścić tę funkcję)

# #' @importFrom jsonlite read_json
# #' @importFrom utils modifyList

#' Validate and complete settings for a specified NMAR imputation method.
#'
#' This function dynamically loads the schema for a given method from a JSON file,
#' applies default values, and performs general validations (type, allowed values, min/max).
#' It then calls method-specific validation if defined.
#'
#' @param method_name A character string specifying the name of the imputation method
#'   (e.g., "exptilt"). This name must correspond to an entry in the
#'   'methods' object of your `inst/extdata/method_params.json` file.
#' @param settings A list of user-provided settings for the specified method.
#'   Any unspecified parameters will be filled with their default values from the schema.
#' @return A validated and complete list of settings for the specified method.
#' @keywords internal
validate_method_settings <- function(method_name, settings = list()) {
browser()
  # --- 1. Load the method schemas directly from the JSON file ---
  json_path <- system.file("extdata", "method_params.json", package = "nmar")

  if (!file.exists(json_path)) {
    stop("Method parameters JSON file (method_params.json) not found in 'inst/extdata/'. Package might be corrupted or not installed correctly.")
  }

  all_schemas <- jsonlite::read_json(json_path, simplifyVector = TRUE)

  # Check if the requested method exists in the schema
  if (!method_name %in% names(all_schemas$methods)) {
    stop(paste0("Unknown method '", method_name, "'. Available methods in schema: ", paste(names(all_schemas$methods), collapse = ", "), "."))
  }

  method_schema <- all_schemas$methods[[method_name]]

  # --- 2. Initialize final_settings with default values from the schema ---
  final_settings <- list()
  for (param_def in method_schema$parameters) {
    if (!is.null(param_def$default)) {
      final_settings[[param_def$name]] <- param_def$default
      cat(paste0("Setting '", param_def$name, "' initialized with default value: ", deparse(param_def$default), "\n"))
    }
  }

  # --- 3. Merge user-provided settings, overwriting defaults ---
  final_settings <- modifyList(final_settings, settings)

  # --- 4. Perform general dynamic validation based on the schema ---
  for (param_def in method_schema$parameters) {
    param_name <- param_def$name
    param_value <- final_settings[[param_name]]

    # # Skip validation for truly absent optional parameters (no default, not provided)
    # if (is.null(param_value) && !isTRUE(param_def$required) && is.null(param_def$default)) {
    #   next
    # }

    # Validate 'required' status
    if (isTRUE(param_def$required) && is.null(param_value)) {
      stop(paste0("Setting '", param_name, "' is required for '", method_name, "' method but is missing."))
    }

    # Validate 'type'
    if (!is.null(param_def$type) && !is.null(param_value)) {
      is_type_valid <- switch(
        param_def$type,
        "string" = is.character(param_value) && length(param_value) == 1,
        "number" = is.numeric(param_value) && length(param_value) == 1,
        "integer" = is.numeric(param_value) && param_value == as.integer(param_value) && length(param_value) == 1,
        "boolean" = is.logical(param_value) && length(param_value) == 1,
        FALSE # Default for unknown types
      )
      if (!is_type_valid) {
        stop(paste0("Setting '", param_name, "' must be a single ", param_def$type, " for '", method_name, "' method (received: ", deparse(param_value), ", type: ", class(param_value)[1], ")."))
      }
    }

    # Validate 'allowed_values'
    if (!is.null(param_def$allowed_values) && !is.null(param_value)) {
      if (!(param_value %in% param_def$allowed_values)) {
        stop(paste0("Setting '", param_name, "' has an invalid value '", param_value, "' for '", method_name, "' method. Allowed values are: ", paste(param_def$allowed_values, collapse = ", "), "."))
      }
    }

    # Validate 'min'
    if (!is.null(param_def$min) && !is.null(param_value) && is.numeric(param_value)) {
      if (param_value < param_def$min) {
        stop(paste0("Setting '", param_name, "' for '", method_name, "' method must be greater than or equal to ", param_def$min, " (received: ", param_value, ")."))
      }
      # You can add 'max' validation here if you add it to your schema.
    }
  }

  # --- 5. Call method-specific validations (placeholder) ---
  # These are validations that are unique to a specific method
  # and cannot be described solely by the general schema rules (e.g., param1 < param2)
  if (method_name == "exptilt") {
    # Example: Cross-parameter validation specific to 'exptilt'
    if (final_settings$min_iter >= final_settings$max_iter) {
      stop("Setting 'min_iter' must be strictly less than 'max_iter' for 'exptilt' method.")
    }
    # Add other 'exptilt' specific validations here.
    # For very complex specific validations, you could encapsulate them in another
    # internal function, e.g.: final_settings <- .validate_exptilt_specific_rules(final_settings)
  }
  # Add 'else if (method_name == "another_method") { ... }' for future methods' specific validations.

  return(final_settings)
}
