il_force_mixed_01 <- function(r) {
  r <- as.integer(r)
  if (length(r) == 0L) return(r)
  if (all(r == 1L)) r[1L] <- 0L
  if (all(r == 0L)) r[1L] <- 1L
  r
}

il_sim_iid_df <- function(n,
                          mu_fn,
                          p_fn,
                          x1 = rnorm(n),
                          x2 = rnorm(n),
                          eps = rnorm(n)) {
  stopifnot(is.numeric(n), length(n) == 1L, n >= 2L)
  stopifnot(is.function(mu_fn), is.function(p_fn))
  if (length(x1) != n) stop("x1 length mismatch")
  if (length(x2) != n) stop("x2 length mismatch")
  if (length(eps) != n) stop("eps length mismatch")

  mu <- as.numeric(mu_fn(x1 = x1, x2 = x2))
  if (length(mu) != n) stop("mu_fn must return length n")
  y <- mu + eps

  p <- as.numeric(p_fn(x1 = x1, x2 = x2, y = y, mu = mu))
  if (length(p) != n) stop("p_fn must return length n")
  r <- il_force_mixed_01(rbinom(n, 1, p))

  data.frame(
    y_miss = ifelse(r == 1L, y, NA_real_),
    x1 = x1,
    x2 = x2,
    y = y,
    mu = mu,
    r = r,
    check.names = FALSE
  )
}
