#' Numeric thresholds and options
"_PACKAGE"

#' NMAR numeric settings
#'
#' @details
#' Centralized access to numeric thresholds used across the package.
#' These settings can be overridden via options() for advanced users:
#' - `nmar.eta_cap`: scalar > 0. Caps the response-model linear predictor
#'   to avoid extreme link values in Newton updates. Default 50.
#' - `nmar.grad_eps`: finite-difference step size epsilon for numeric
#'   gradients of smooth functionals. Default 1e-6.
#' - `nmar.grad_d`: relative step adjustment for numeric gradients.
#'   Default 1e-3.
#'
#' The defaults are chosen to be conservative and stable across typical
#' NMAR problems. Engines should retrieve values via this helper rather
#' than hard-coding numbers, so documentation stays consistent.
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

#' Numeric gradient helper
#' @description Thin wrapper around numDeriv::grad for internal use across engines.
#' @param x Numeric vector at which to evaluate the gradient.
#' @param func Function mapping numeric vector to scalar; receives `x`.
#'
#' @keywords internal
grad_numeric <- function(x, func) {
  s <- nmar_get_numeric_settings()
  numDeriv::grad(func = func, x = x, method.args = list(d = s$grad_d, eps = s$grad_eps))
}

#' Weighted linear algebra helpers
#'
#' @keywords internal
shared_weighted_gram <- function(X, w) {
# Compute X' diag(w) X efficiently. If w >= 0, use SPD crossprod(X*sqrt(w))
# Otherwise, fall back to X' (diag(w) X) via crossprod(X, X*w)
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
# Compute X' (w ∘ y) with elementwise row weights
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



#' Sanitize nleqslv control list for compatibility
#' @keywords internal
sanitize_nleqslv_control <- function(ctrl) {
  if (is.null(ctrl) || !is.list(ctrl)) return(list())
# Keep a conservative whitelist to avoid unknown-name errors on older versions
  allowed <- c("xtol", "ftol", "btol", "maxit", "trace", "stepmax", "delta", "allowSing")
  ctrl[names(ctrl) %in% allowed]
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

#' Validate top-level nleqslv arguments (coerce invalid to defaults)
#' @keywords internal
validate_nleqslv_top <- function(top) {
  if (is.null(top) || !is.list(top)) return(list())
  out <- list()
  if (!is.null(top$global)) {
    allowed_g <- c("dbldog", "pwldog", "qline", "none")
    out$global <- if (is.character(top$global) && top$global[1] %in% allowed_g) top$global[1] else "dbldog"
  }
  if (!is.null(top$xscalm)) {
    allowed_x <- c("auto", "fixed")
    out$xscalm <- if (is.character(top$xscalm) && top$xscalm[1] %in% allowed_x) top$xscalm[1] else "auto"
  }
# We ignore method from user-facing API to keep a single solver path
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
