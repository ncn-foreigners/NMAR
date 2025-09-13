validate_method_settings <- function(method_name, arguments, schemas) {
  # Pobierz parametry dla wybranej metody
  method_params <- schemas$methods[[method_name]]$parameters

  # Sprawdź czy metoda istnieje
  if (is.null(method_params)) {
    stop("Method '", method_name, "' not found in schema")
  }

  # Przekształć na data.frame dla łatwiejszej manipulacji
  params_df <- do.call(rbind, lapply(method_params, as.data.frame))

  # Walidacja każdego przekazanego argumentu
  for (arg_name in names(arguments)) {
    # Znajdź specyfikację parametru
    param_spec <- params_df[params_df$name == arg_name, ]

    # Sprawdź czy parametr istnieje w schemacie
    if (nrow(param_spec) == 0) {
      stop("Parameter '", arg_name, "' is not valid for method '", method_name, "'")
    }

    # Pobierz wartość argumentu
    arg_value <- arguments[[arg_name]]

    # Walidacja wartości dozwolonych (dla typów string)
    if (param_spec$type == "string" && !is.null(param_spec$allowed_values)) {
      if (!arg_value %in% param_spec$allowed_values[[1]]) {
        stop("Invalid value '", arg_value, "' for parameter '", arg_name,
             "'. Allowed values: ", paste(param_spec$allowed_values[[1]], collapse = ", "))
      }
    }

    # Walidacja przedziałów (dla typów number/integer)
    if (!is.null(param_spec$min)) {
      if (arg_value < param_spec$min) {
        stop("Parameter '", arg_name, "' must be greater than or equal to ", param_spec$min)
      }
    }

    if (!is.null(param_spec$max)) {
      if (arg_value > param_spec$max) {
        stop("Parameter '", arg_name, "' must be less than or equal to ", param_spec$max)
      }
    }

    # Walidacja typu danych
    expected_type <- param_spec$type
    if (expected_type == "integer" && !is.integer(arg_value)) {
      if (is.numeric(arg_value) && arg_value == round(arg_value)) {
        arguments[[arg_name]] <- as.integer(arg_value)
      } else {
        stop("Parameter '", arg_name, "' must be an integer")
      }
    } else if (expected_type == "number" && !is.numeric(arg_value)) {
      stop("Parameter '", arg_name, "' must be a number")
    } else if (expected_type == "string" && !is.character(arg_value)) {
      stop("Parameter '", arg_name, "' must be a string")
    }
  }

  # Uzupełnij brakujące argumenty wartościami domyślnymi
  missing_args <- setdiff(params_df$name, names(arguments))
  for (arg_name in missing_args) {
    param_spec <- params_df[params_df$name == arg_name, ]
    if (!param_spec$required) {
      arguments[[arg_name]] <- param_spec$default[[1]]
    } else {
      stop("Required parameter '", arg_name, "' is missing")
    }
  }

  return(arguments)
}
