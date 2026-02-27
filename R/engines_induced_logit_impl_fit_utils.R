#' Induced-logit fitting utilities
#'
#' @keywords internal
#' @noRd
il_assert_scalar_logical <- function(x, name) {
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("Internal error: `name` must be a nonempty string.", call. = FALSE)
  }
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("Internal error: `", name, "` must be TRUE/FALSE.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
il_quote_name_for_formula <- function(x) {
  x <- as.character(x)
  x <- gsub("`", "''", x, fixed = TRUE)
  paste0("`", x, "`")
}

#' @keywords internal
#' @noRd
il_reformulate <- function(response, term_labels, intercept = TRUE) {
  if (!is.character(response) || length(response) != 1L || is.na(response) || !nzchar(response)) {
    stop("Internal error: `response` must be a nonempty string.", call. = FALSE)
  }
  if (is.null(term_labels)) term_labels <- character(0)
  if (!is.character(term_labels)) stop("Internal error: `term_labels` must be character().", call. = FALSE)
  stats::reformulate(termlabels = term_labels, response = response, intercept = isTRUE(intercept))
}

#' @keywords internal
#' @noRd
il_glm_control <- function(control) {
  if (is.null(control)) control <- list()
  validator_assert_list(control, name = "control")
  out <- tryCatch(do.call(stats::glm.control, control), error = function(e) e)
  if (inherits(out, "error")) {
    stop("Invalid `control` for glm.control(): ", out$message, call. = FALSE)
  }
  out
}

#' @keywords internal
#' @noRd
il_scale_matrix <- function(mat, weights = NULL, weight_mask = NULL) {
  if (!is.matrix(mat)) stop("Internal error: expected matrix.", call. = FALSE)
  if (ncol(mat) == 0L) return(list(mat = mat, recipe = NULL))
  recipe <- create_nmar_scaling_recipe(
    mat,
    intercept_col = "(Intercept)",
    weights = weights,
    weight_mask = weight_mask,
    warn_on_constant = FALSE
  )
  list(mat = apply_nmar_scaling(mat, recipe), recipe = recipe)
}

#' Scale a mu-model matrix as a pure reparameterization
#'
#' Centering predictors is a reparameterization only when the mu-model includes an
#' intercept. When the intercept is omitted, centering changes the column space
#' (equivalent to implicitly adding a constant column) and therefore changes the
#' fitted values. For no-intercept mu models we therefore scale (divide by sd)
#' but do not center.
#'
#' @keywords internal
#' @noRd
il_scale_mu_matrix <- function(mat, has_intercept, weights = NULL, weight_mask = NULL) {
  if (!is.matrix(mat)) stop("Internal error: expected matrix.", call. = FALSE)
  if (ncol(mat) == 0L) return(list(mat = mat, recipe = NULL))

  recipe <- create_nmar_scaling_recipe(
    mat,
    intercept_col = "(Intercept)",
    weights = weights,
    weight_mask = weight_mask,
    warn_on_constant = FALSE
  )
  if (!isTRUE(has_intercept)) {
    for (nm in names(recipe)) recipe[[nm]]$mean <- 0
  }
  recipe <- validate_nmar_scaling_recipe(recipe)
  list(mat = apply_nmar_scaling(mat, recipe), recipe = recipe)
}

#' @keywords internal
#' @noRd
il_response_model_diagnostics_matrix <- function(mm) {
  if (!is.matrix(mm)) stop("Internal error: response-model matrix must be a matrix.", call. = FALSE)
  if (ncol(mm) == 0L) {
    return(list(rank = 0L, ncol = 0L, condition_number = NA_real_))
  }

  qr_mm <- qr(mm)
  rank_mm <- as.integer(qr_mm$rank)
  p <- as.integer(ncol(mm))

  kappa_mm <- tryCatch(
    as.numeric(kappa(mm, exact = FALSE)),
    error = function(e) NA_real_
  )
  if (!is.finite(kappa_mm)) {
    d <- tryCatch(svd(mm, nu = 0, nv = 0)$d, error = function(e) numeric(0))
    if (length(d) > 0L && all(is.finite(d)) && max(d) > 0) {
      dmin <- suppressWarnings(min(d[d > 0]))
      if (is.finite(dmin) && dmin > 0) {
        kappa_mm <- as.numeric(max(d) / dmin)
      }
    }
  }

  list(rank = rank_mm, ncol = p, condition_number = as.numeric(kappa_mm))
}

#' @keywords internal
#' @noRd
il_enforce_response_model_identifiability <- function(diag,
                                                      condition_warn = 1e8,
                                                      condition_fail = 1e12) {
  if (!is.list(diag)) stop("Internal error: `diag` must be a list.", call. = FALSE)
  rank_mm <- as.integer(diag$rank %||% NA_integer_)
  p <- as.integer(diag$ncol %||% NA_integer_)
  kappa_mm <- as.numeric(diag$condition_number %||% NA_real_)

  if (!is.finite(rank_mm) || !is.finite(p) || rank_mm < p) {
    stop(
      "Induced-logit response model is not identifiable: design matrix is rank deficient.\n  ",
      sprintf("rank = %s, ncol = %s", as.character(rank_mm), as.character(p)),
      call. = FALSE
    )
  }

  if (!is.finite(kappa_mm) || kappa_mm > condition_fail) {
    stop(
      "Induced-logit response model is nearly non-identifiable (extreme conditioning).\n  ",
      sprintf("condition number = %s (threshold = %g).", format(kappa_mm, digits = 6), condition_fail),
      call. = FALSE
    )
  }

  warn_msg <- NULL
  if (kappa_mm > condition_warn) {
    warn_msg <- paste0(
      "Induced-logit response model is weakly identified: high condition number (",
      format(kappa_mm, digits = 6), ")."
    )
  }

  list(diag = diag, warning = warn_msg)
}
