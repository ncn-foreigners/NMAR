#' Numerical helpers
#' @keywords internal
get_eta_cap <- function() {
  cap <- getOption("nmar.eta_cap", default = 50)
  if (!is.numeric(cap) || length(cap) != 1L || !is.finite(cap) || cap <= 0) cap <- 50
  cap
}

#' Numeric gradient helper
#' @description Thin wrapper around numDeriv::grad for internal use across engines.
#' @param x Numeric vector at which to evaluate the gradient.
#' @param func Function mapping numeric vector to scalar; receives `x`.
#' @keywords internal
#' @noRd
grad_numeric <- function(x, func) {
# Allow users to tune finite-difference settings via options; fall back to stable defaults
  eps <- getOption("nmar.grad_eps", default = 1e-6)
  d <- getOption("nmar.grad_d", default = 1e-3)
  numDeriv::grad(func = func, x = x, method.args = list(d = d, eps = eps))
}

#' Invert a Jacobian matrix with numerically robust fallbacks.
#'
#' Attempts a plain inverse when the matrix is well conditioned, otherwise
#' applies ridge regularisation or a pseudo-inverse as requested.
#'
#' @param A_matrix Numeric square matrix to invert.
#' @param variance_ridge Logical or positive numeric: apply ridge stabilisation
#'   (`TRUE` selects an adaptive ridge, numeric values supply the ridge size).
#' @param variance_pseudoinverse Logical; if `TRUE`, compute an SVD-based
#'   pseudo-inverse.
#' @param kappa_threshold Numeric tolerance controlling when the plain inverse is
#'   attempted.
#' @param ridge_scale Optional numeric multiplier for the adaptive ridge.
#' @param svd_tol Optional numeric tolerance below which singular values are
#'   truncated when forming the pseudo-inverse.
#' @keywords internal
invert_jacobian <- function(A_matrix,
                            variance_ridge = FALSE,
                            variance_pseudoinverse = FALSE,
                            kappa_threshold = 1e8,
                            ridge_scale = NULL,
                            svd_tol = NULL) {
# Compute condition number (may error on singular matrices)
  kappa_val <- tryCatch(kappa(A_matrix), error = function(e) Inf)

# 1) Try plain inverse if seemingly well-conditioned.
  if (is.finite(kappa_val) && kappa_val <= kappa_threshold) {
    plain <- tryCatch(
      {
        list(inv = -solve(A_matrix), invert_rule = "plain", used_pinv = FALSE, used_ridge = FALSE, kappa = kappa_val)
      },
      error = function(e1) NULL
    )
    if (!is.null(plain)) {
      return(plain)
    }
  }

# 2) If ridge requested, attempt adaptive ridge.
# variance_ridge may be logical or a positive numeric epsilon.
  if (!identical(variance_ridge, FALSE)) {
    eps <- NA_real_
    if (isTRUE(variance_ridge)) {
# Adaptive epsilon: scale by spectral norm
      smax <- tryCatch(max(svd(A_matrix, nu = 0, nv = 0)$d), error = function(e) NA_real_)
      if (!is.finite(smax) || smax <= 0) smax <- 1
# Make the base epsilon more aggressive; allow override via option
      base_opt <- getOption("nmar.var_ridge_base", default = 1e-6)
      base <- if (is.numeric(ridge_scale) && is.finite(ridge_scale) && ridge_scale > 0) ridge_scale else base_opt
      eps <- base * smax
    } else if (is.numeric(variance_ridge) && is.finite(variance_ridge) && variance_ridge > 0) {
      eps <- variance_ridge
    }
    if (is.finite(eps) && eps > 0) {
      A_ridge <- A_matrix + diag(eps, nrow(A_matrix))
      ridged <- tryCatch(
        {
          list(inv = -solve(A_ridge), invert_rule = "ridge", used_pinv = FALSE, used_ridge = TRUE, kappa = kappa_val, ridge_epsilon = eps)
        },
        error = function(e2) NULL
      )
      if (!is.null(ridged)) {
        return(ridged)
      }
    }
  }

# 3) If pseudo-inverse requested, use an SVD-based construction.
  if (isTRUE(variance_pseudoinverse)) {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Jacobian matrix is singular; pseudo-inverse requested but 'MASS' is not installed.", call. = FALSE)
    }
    if (is.null(svd_tol)) {
      pinv <- MASS::ginv(A_matrix)
    } else {
      s <- svd(A_matrix)
      d <- s$d
      d_inv <- if (length(d) == 0) numeric(0) else ifelse(d > svd_tol, 1 / d, 0)
      if (length(d_inv) == 0) {
        pinv <- matrix(0, nrow = ncol(A_matrix), ncol = nrow(A_matrix))
      } else {
        pinv <- s$v %*% (diag(d_inv, nrow = length(d_inv), ncol = length(d_inv)) %*% t(s$u))
      }
    }
    return(list(inv = -pinv, invert_rule = "pinv", used_pinv = TRUE, used_ridge = FALSE, kappa = kappa_val))
  }

# 4) Final attempt: plain inverse (even if ill-conditioned), else error.
  plain2 <- tryCatch(
    {
      list(inv = -solve(A_matrix), invert_rule = "plain", used_pinv = FALSE, used_ridge = FALSE, kappa = kappa_val)
    },
    error = function(e3) NULL
  )
  if (!is.null(plain2)) {
    return(plain2)
  }

  stop("Jacobian is singular or ill-conditioned; try variance_pseudoinverse=TRUE or variance_ridge=TRUE, or use variance_method='bootstrap'.", call. = FALSE)
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
