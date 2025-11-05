#' Inner Wu multiplier solve for EL (block method)
#'
#' Solves g(lambda) = sum_i d_i * u_i / (1 + u_i^T lambda) = 0 for
#' lambda = (lambda_x, lambda_W) at fixed (beta, W), using a modified Newton
#' method with backtracking that enforces positivity of denominators and
#' monotone decrease of the merit ||g(lambda)||_∞ (Wu's modified NR).
#'
#' @param U n x p matrix with rows u_i = c(Xc_i, w_i(beta) - W).
#' @param d numeric length-n vector of nonnegative base weights (respondents).
#' @param denom_floor small positive floor for denominators (global guard).
#' @param kappa feasibility buffer multiplier (>= 1); require D_i >= kappa*denom_floor.
#' @param tol gradient infinity-norm tolerance.
#' @param maxit maximum iterations for Newton updates.
#' @return list(lambda, D, ok, iters, gnorm, Lval)
#' @keywords internal
el_inner_wu_solve <- function(U, d, denom_floor, kappa = 10, tol = 1e-8, maxit = 50) {
  n <- nrow(U)
  p <- ncol(U)
  if (length(d) != n) stop("el_inner_wu_solve: length(d) must equal nrow(U)", call. = FALSE)
  d <- as.numeric(d)
  if (any(d < 0)) stop("el_inner_wu_solve: weights d must be nonnegative", call. = FALSE)

# Start at zero multipliers (feasible: D_i = 1)
  lambda <- rep(0, p)
  names(lambda) <- colnames(U)

# Helper closures
  compute_state <- function(lam) {
    D <- as.numeric(1 + U %*% lam)
    D <- pmax(D, .Machine$double.xmin) # temporary to avoid NaNs in log
    g <- as.numeric(t(U) %*% (d / D)) # gradient of L(lambda) = sum d log D
    L <- sum(d * log(D)) # CONVEX objective to minimize
    list(D = D, g = g, L = L)
  }

  st <- compute_state(lambda)
  gnorm <- max(abs(st$g))
  it <- 0L
  ok <- TRUE

# Newton iterations on root g(lambda) = 0 with line search on ||g||
  while (gnorm > tol && it < maxit) {
    it <- it + 1L
    D <- st$D
# Build SPD Hessian H = sum d_i * (u_i u_i^T) / D_i^2
    Wdiag <- (d / (D^2))
# Compute crossprod(U * sqrt(Wdiag)) to avoid forming large diag
    Uw <- sweep(U, 1, sqrt(Wdiag), "*")
    H <- crossprod(Uw)
# Small ridge to handle near singularity
    diag(H) <- diag(H) + 1e-12
# Newton step for root-finding: solve J Δ = -g where J = -H ⇒ H Δ = g
    step <- tryCatch(solve(H, st$g), error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) {
      ok <- FALSE
      break
    }
# Backtracking on the merit ||g|| (sufficient decrease) with positivity guard
    alpha <- 1.0
# Precompute feasibility threshold
    thr <- kappa * denom_floor
    improved <- FALSE
    g_old <- st$g
    for (ls in 1:30) {
      lam_new <- lambda + alpha * step
      D_new <- as.numeric(1 + U %*% lam_new)
      if (min(D_new) >= thr) {
        g_new <- as.numeric(t(U) %*% (d / D_new))
# sufficient reduction in infinity norm of gradient
        if (max(abs(g_new)) <= (1 - 1e-4 * alpha) * gnorm) {
          improved <- TRUE
          break
        }
      }
      alpha <- alpha / 2
    }
    if (!improved) {
      ok <- FALSE
      break
    }
# Accept
    lambda <- lam_new
    g_new <- as.numeric(t(U) %*% (d / D_new))
    st <- list(D = D_new, g = g_new, L = sum(d * log(D_new)))
    gnorm <- max(abs(st$g))
  }

  list(lambda = lambda, D = st$D, ok = ok && (gnorm <= tol), iters = it, gnorm = gnorm, Lval = st$L)
}

#' Inner solve for lambda_x only with fixed lambda_W
#'
#' Solves g_x(lambda_x) = sum_i d_i Xc_i / (offset_i + Xc_i^T lambda_x) = 0
#' where offset_i = 1 + lambda_W * (w_i - W) is fixed (QLS eq. 10).
#' This reduces the inner dimension and improves stability.
#'
#' @param Xc n x p matrix of centered auxiliaries.
#' @param offset length-n vector: 1 + lambda_W * (w - W).
#' @inheritParams el_inner_wu_solve
#' @return list(lambda_x, D, ok, iters, gnorm)
#' @keywords internal
el_inner_wu_solve_xonly <- function(Xc, offset, d, denom_floor, kappa = 10, tol = 1e-8, maxit = 50) {
  n <- nrow(Xc)
  p <- ncol(Xc)
  if (length(d) != n) stop("el_inner_wu_solve_xonly: length(d) must equal nrow(Xc)", call. = FALSE)
  d <- as.numeric(d)
  if (any(d < 0)) stop("el_inner_wu_solve_xonly: weights d must be nonnegative", call. = FALSE)
  if (length(offset) != n) stop("el_inner_wu_solve_xonly: offset length mismatch", call. = FALSE)

  lambda <- rep(0, p)
  names(lambda) <- colnames(Xc)

  compute_state <- function(lam) {
    D <- as.numeric(offset + Xc %*% lam)
    D <- pmax(D, .Machine$double.xmin)
    g <- as.numeric(t(Xc) %*% (d / D))
    list(D = D, g = g)
  }

  st <- compute_state(lambda)
  gnorm <- max(abs(st$g))
  it <- 0L
  ok <- TRUE
  thr <- kappa * denom_floor

  while (gnorm > tol && it < maxit) {
    it <- it + 1L
    D <- st$D
    Wdiag <- (d / (D^2))
    Xw <- sweep(Xc, 1, sqrt(Wdiag), "*")
    H <- crossprod(Xw)
    diag(H) <- diag(H) + 1e-12
    step <- tryCatch(solve(H, st$g), error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) { ok <- FALSE; break }
    alpha <- 1.0
    improved <- FALSE
    for (ls in 1:30) {
      lam_new <- lambda + alpha * step
      D_new <- as.numeric(offset + Xc %*% lam_new)
      if (min(D_new) >= thr) {
        g_new <- as.numeric(t(Xc) %*% (d / D_new))
        if (max(abs(g_new)) <= (1 - 1e-4 * alpha) * gnorm) {
          improved <- TRUE
          break
        }
      }
      alpha <- alpha / 2
    }
    if (!improved) { ok <- FALSE; break }
    lambda <- lam_new
    st <- list(D = D_new, g = g_new)
    gnorm <- max(abs(st$g))
  }
  list(lambda_x = lambda, D = st$D, ok = ok && (gnorm <= tol), iters = it, gnorm = gnorm)
}
