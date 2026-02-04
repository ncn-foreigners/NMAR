#' NMAR numeric settings
#'
#' @details
#' Centralized access to numeric thresholds used across the package.
#'
#' - `nmar.eta_cap`: scalar > 0. Caps the response-model linear predictor
#' to avoid extreme link values in Newton updates. Default 50.
#' - `nmar.grad_eps`: finite-difference step size epsilon for numeric
#' gradients of smooth functionals. Default 1e-6.
#' - `nmar.grad_d`: relative step adjustment for numeric gradients.
#' Default 1e-3.
#'
#' @return A named list with entries `eta_cap`, `grad_eps`, and `grad_d`.
#' @keywords internal
nmar_get_numeric_settings <- function() {
  eta_cap <- getOption("nmar.eta_cap", default = 50)
  if (!is.numeric(eta_cap) || length(eta_cap) != 1L || !is.finite(eta_cap) || eta_cap <= 0) eta_cap <- 50
  grad_eps <- getOption("nmar.grad_eps", default = 1e-6)
  if (!is.numeric(grad_eps) || length(grad_eps) != 1L || !is.finite(grad_eps) || grad_eps <= 0) grad_eps <- 1e-6
  grad_d <- getOption("nmar.grad_d", default = 1e-3)
  if (!is.numeric(grad_d) || length(grad_d) != 1L || !is.finite(grad_d) || grad_d <= 0) grad_d <- 1e-3
  list(eta_cap = eta_cap, grad_eps = grad_eps, grad_d = grad_d)
}

#' @keywords internal
get_eta_cap <- function() {
  nmar_get_numeric_settings()$eta_cap
}

#' EL denominator floor
#'
#' Returns the small positive floor \eqn{\delta} used to guard the empirical
#' likelihood denominator \eqn{D_i(\theta)} away from zero.
#'
#' @keywords internal
nmar_get_el_denom_floor <- function() {
  val <- getOption("nmar.el_denom_floor", 1e-8)
  if (!is.numeric(val) || length(val) != 1L || !is.finite(val) || val <= 0) 1e-8 else val
}

#' Weighted linear algebra
#'
#' Compute X' diag(w) X efficiently. If w >= 0, use SPD crossprod(X*sqrt(w)).
#' Otherwise, fall back to X' (diag(w) X) via crossprod(X, X*w).
#'
#' @keywords internal
shared_weighted_gram <- function(X, w) {
  w <- as.numeric(w)
  if (length(w) != nrow(X)) stop("shared_weighted_gram: length(w) must equal nrow(X)", call. = FALSE)
  if (all(is.finite(w)) && all(w >= 0)) {
    crossprod(X * sqrt(w))
  } else {
    crossprod(X, X * w)
  }
}

#' @keywords internal
shared_weighted_Xty <- function(X, w, y) {
# Compute X' (w * y) with elementwise row weights
  w <- as.numeric(w)
  y <- as.numeric(y)
  if (length(w) != nrow(X) || length(y) != nrow(X)) stop("shared_weighted_Xty: lengths must match nrow(X)", call. = FALSE)
  crossprod(X, w * y)
}

#' @keywords internal
shared_weighted_XtY <- function(X, w, Y) {
# Compute X' (diag(w) Y) efficiently via row-wise scaling of Y
  w <- as.numeric(w)
  if (length(w) != nrow(X) || nrow(Y) != nrow(X)) stop("shared_weighted_XtY: nrow mismatch", call. = FALSE)
  crossprod(X, sweep(Y, 1, w, "*"))
}



#' Parse nleqslv control list for compatibility
#' @keywords internal
sanitize_nleqslv_control <- function(ctrl) {
  if (is.null(ctrl) || !is.list(ctrl)) return(list())

  allowed <- c("xtol", "ftol", "btol", "maxit", "trace", "stepmax", "delta", "allowSing")

  ctrl <- ctrl[setdiff(names(ctrl), c("global", "xscalm", "method"))]
  unknown <- setdiff(names(ctrl), allowed)
  if (length(unknown) > 0) {
    warning(sprintf("Ignoring unknown nleqslv control fields: %s", paste(unknown, collapse = ", ")), call. = FALSE)
  }
  out <- ctrl[names(ctrl) %in% allowed]

  num_pos <- function(x, nm) {
    if (!is.null(x) && (!is.finite(x) || x <= 0)) {
      warning(sprintf("Coercing control$%s to a positive finite value; using default.", nm), call. = FALSE)
      return(NULL)
    }
    x
  }
  out$xtol <- num_pos(out$xtol, "xtol")
  out$ftol <- num_pos(out$ftol, "ftol")
  out$btol <- num_pos(out$btol, "btol")
  out$maxit <- if (!is.null(out$maxit) && (!is.finite(out$maxit) || out$maxit <= 0)) {
    warning("Coercing control$maxit to a positive integer; using default.", call. = FALSE)
    NULL
  } else out$maxit
  if (!is.null(out$trace) && !is.logical(out$trace)) {
    warning("Coercing control$trace to logical; using default.", call. = FALSE)
    out$trace <- NULL
  }
  out$stepmax <- num_pos(out$stepmax, "stepmax")
  out$delta <- num_pos(out$delta, "delta")
  out
}

#' Extract top-level nleqslv arguments from a control-like list
#' @keywords internal
extract_nleqslv_top <- function(ctrl) {
  if (is.null(ctrl) || !is.list(ctrl)) return(list())
  out <- list()
  if (!is.null(ctrl$global)) out$global <- ctrl$global
  if (!is.null(ctrl$xscalm)) out$xscalm <- ctrl$xscalm
  if (!is.null(ctrl$method)) out$method <- ctrl$method
  out
}

#' Validate top-level nleqslv arguments and coerce invalid to defaults
#' @keywords internal
validate_nleqslv_top <- function(top) {
  if (is.null(top) || !is.list(top)) return(list())
  out <- list()
  if (!is.null(top$global)) {
    allowed_g <- c("dbldog", "pwldog", "qline", "cline", "gline", "hook", "none")
    if (is.character(top$global) && top$global[1] %in% allowed_g) {
      out$global <- top$global[1]
    } else {
      warning("Unknown nleqslv 'global' value in control; using default 'qline'.", call. = FALSE)
      out$global <- "qline"
    }
  }
  if (!is.null(top$xscalm)) {
    allowed_x <- c("auto", "fixed")
    if (is.character(top$xscalm) && top$xscalm[1] %in% allowed_x) {
      out$xscalm <- top$xscalm[1]
    } else {
      warning("Unknown nleqslv 'xscalm' value in control; using default 'auto'.", call. = FALSE)
      out$xscalm <- "auto"
    }
  }

  out
}

#' Prefer explicit solver_args over control-provided top-level args
#' @keywords internal
merge_nleqslv_top <- function(solver_args, control_top) {
  res <- control_top %||% list()
  if (!is.list(res)) res <- list()
  if (is.list(solver_args)) {
    if (!is.null(solver_args$global)) res$global <- solver_args$global
    if (!is.null(solver_args$xscalm)) res$xscalm <- solver_args$xscalm
  }
  res
}
