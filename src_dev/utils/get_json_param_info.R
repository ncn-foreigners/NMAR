get_json_param_info <- function(schema, method_name, param_name) {
  params <- schema$methods[[method_name]]$parameters
  param_info <- params[params$name == param_name, ]

  if (nrow(param_info) == 0) {
    stop("Parameter '", param_name, "' not found in method '", method_name, "'")
  }

  as.list(param_info)
}

