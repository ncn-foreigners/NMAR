#' Concentrated log-likelihood objective (block method)
#'
#' Builds closures to evaluate the concentrated log-likelihood and its
#' gradient with respect to (beta, z) where W = logistic(z), at fixed
#' design matrices (scaled) and respondent weights. The inner multiplier
#' (lambda_x, lambda_W) is obtained via the Wu solver at each evaluation.
#'
#' @param family response model family (logit/probit) with linkinv and mu.eta.
#' @param Z n x K_beta response-model matrix (scaled).
#' @param Xc n x K_aux centered auxiliary matrix (scaled; may have zero columns).
#' @param d length-n respondent weights (nonnegative).
#' @param N_pop scalar population total (on same scale as d).
#' @param denom_floor small positive floor for denominators.
#' @return list with functions: value(theta), gradient(theta), and last_inner()
#'   that returns a list with lambda, D, ok from the most recent evaluation.
#' @keywords internal
el_build_concentrated_obj <- function(family, Z, Xc, d, N_pop, denom_floor,
                                      kappa = 10, inner_tol = 1e-10, inner_maxit = 80,
                                      eta_cap = get_eta_cap(),
                                      fix_lambda_W = FALSE,
                                      numeric_grad = FALSE) {
  n <- nrow(Z)
  K_beta <- ncol(Z)
  K_aux <- if (is.null(Xc)) 0L else ncol(Xc)
  if (length(d) != n) stop("el_build_concentrated_obj: length(d) must equal nrow(Z)", call. = FALSE)

# Reusable work buffers
  last <- new.env(parent = emptyenv())
  last$lambda <- NULL
  last$D <- NULL
  last$ok <- FALSE

# Helper to run the inner solve given beta and W
  run_inner <- function(beta, W) {
# Compute response probabilities and derivative
    eta <- as.vector(Z %*% beta)
    eta <- pmax(pmin(eta, eta_cap), -eta_cap)
    w <- family$linkinv(eta)
    w <- pmin(pmax(w, 1e-12), 1 - 1e-12)
    if (isTRUE(fix_lambda_W)) {
# QLS (10) with weighted generalization: lambda_W = (N/sum d - 1)/(1 - W)
      lambda_W <- ((N_pop / sum(d)) - 1) / (1 - W)
      if (K_aux > 0) {
        offset <- as.numeric(1 + lambda_W * (w - W))
        out <- el_inner_wu_solve_xonly(Xc = Xc, offset = offset, d = d,
                                       denom_floor = denom_floor, kappa = kappa,
                                       tol = inner_tol, maxit = inner_maxit)
        lam_full <- c(out$lambda_x, `(w-W)` = lambda_W)
        last$lambda <- lam_full
        last$D <- out$D
        last$ok <- isTRUE(out$ok)
        list(lambda = lam_full, D = out$D, ok = out$ok, w = w)
      } else {
        D <- as.numeric(1 + lambda_W * (w - W))
        last$lambda <- c(`(w-W)` = lambda_W)
        last$D <- D
        last$ok <- all(D > denom_floor)
        list(lambda = last$lambda, D = D, ok = last$ok, w = w)
      }
    } else {
# Full inner solve (lambda_x, lambda_W)
      if (K_aux > 0) {
        U <- cbind(Xc, w - W)
        colnames(U) <- c(colnames(Xc), "(w-W)")
      } else {
        U <- matrix(w - W, ncol = 1)
        colnames(U) <- "(w-W)"
      }
      out <- el_inner_wu_solve(U = U, d = d, denom_floor = denom_floor,
                               kappa = kappa, tol = inner_tol, maxit = inner_maxit)
      last$lambda <- out$lambda
      last$D <- out$D
      last$ok <- isTRUE(out$ok)
      list(lambda = out$lambda, D = out$D, ok = out$ok, w = w)
    }
  }

# Value function: returns negative concentrated loglik for minimizers
  value_fn <- function(theta) {
    beta <- theta[seq_len(K_beta)]
    z <- theta[K_beta + 1L]
    W <- plogis(z)
    W <- min(max(W, 1e-12), 1 - 1e-12)
    inner <- run_inner(beta, W)
    if (!inner$ok || any(!is.finite(inner$D))) {
      return(1e50) # penalty
    }
# ℓ(β, W) = sum d log w + (N_pop - sum d) log(1 - W) - sum d log D
    val <- sum(d * log(inner$w)) + (N_pop - sum(d)) * log(1 - W) - sum(d * log(inner$D))
    return(-val)
  }

  gradient_fn_core <- function(theta) {
    beta <- theta[seq_len(K_beta)]
    z <- theta[K_beta + 1L]
    W <- plogis(z)
    W <- min(max(W, 1e-12), 1 - 1e-12)
    inner <- run_inner(beta, W)
    if (!inner$ok || any(!is.finite(inner$D))) {
# Return a large gradient to push optimizer away from infeasible points
      return(rep(1e6, K_beta + 1L))
    }
    D <- inner$D
    w <- inner$w
    eta <- as.vector(Z %*% beta)
    eta <- pmax(pmin(eta, eta_cap), -eta_cap)
    m <- family$mu.eta(eta)
# Extract lambda_W (last component)
    lambda_W <- as.numeric(last$lambda[length(last$lambda)])
# ∂ℓ/∂β = Σ d m Z (1/w - lambda_W / D)
    term <- d * (m * (1 / w - lambda_W / D))
    grad_beta <- as.numeric(crossprod(Z, term))
# ∂ℓ/∂W = −(N−Σd)/(1−W) + lambda_W Σ d / D
    grad_W <- (-(N_pop - sum(d)) / (1 - W)) + lambda_W * sum(d / D)
# Chain for z: dW/dz = W(1 − W)
    grad_z <- grad_W * (W * (1 - W))
    return(-c(grad_beta, grad_z)) # negative because we minimize -ℓ
  }

  gradient_fn <- function(theta) {
    if (!isTRUE(numeric_grad)) return(gradient_fn_core(theta))
    grad_numeric(theta, value_fn)
  }

  last_inner <- function() {
    list(lambda = last$lambda, D = last$D, ok = last$ok)
  }

  list(value = value_fn, gradient = gradient_fn, last_inner = last_inner)
}
