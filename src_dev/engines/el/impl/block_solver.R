#' Outer optimizers for concentrated EL (block method)
#'
#' Minimizes the negative concentrated log-likelihood using a user-selected
#' optimizer. Requires closures from `el_build_concentrated_obj()` providing
#' value and gradient. The inner Wu solver is invoked inside those closures.
#'
#' Supported methods: "BFGS" (default), "L-BFGS-B", "nlminb".
#'
#' @keywords internal
el_run_block_outer <- function(value_fn, grad_fn, init, method = "BFGS", control = list()) {
  safe_value <- function(x) {
    v <- value_fn(x)
    if (!is.finite(v)) 1e50 else v
  }
  safe_grad <- function(x) {
    g <- grad_fn(x)
    g[!is.finite(g)] <- 0
    g
  }
  method <- toupper(method %||% "BFGS")
  if (method == "NLMINB") {
    op <- tryCatch(
      stats::nlminb(start = init, objective = safe_value, gradient = safe_grad, control = control),
      error = function(e) NULL
    )
    if (is.null(op)) return(list(theta = init, convergence = 1L, message = "nlminb failed", value = NA_real_, counts = list()))
    return(list(theta = op$par, convergence = if (is.null(op$convergence)) 0L else op$convergence, message = op$message %||% "", value = op$objective, counts = list(evaluations = op$evaluations)))
  }
  if (method == "L-BFGS-B") {
# Optional bound on z (last coordinate): keep W = plogis(z) in a stable interior
    lower <- rep(-Inf, length(init)); upper <- rep(Inf, length(init))
    if (length(init) >= 1) {
# z bounds ~ logit([1e-6, 1-1e-6])
      lower[length(init)] <- qlogis(1e-6); upper[length(init)] <- qlogis(1 - 1e-6)
    }
    op <- tryCatch(
      optim(par = init, fn = safe_value, gr = safe_grad, method = "L-BFGS-B", lower = lower, upper = upper, control = control),
      error = function(e) NULL
    )
  } else {
    op <- tryCatch(
      optim(par = init, fn = safe_value, gr = safe_grad, method = "BFGS", control = control),
      error = function(e) NULL
    )
  }
  if (is.null(op)) return(list(theta = init, convergence = 1L, message = "optim failed", value = NA_real_, counts = list()))
  list(theta = op$par, convergence = op$convergence, message = op$message %||% "", value = op$value, counts = op$counts)
}
