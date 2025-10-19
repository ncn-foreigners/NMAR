test_that("Mean and beta variance identities hold in a benign aux cell", {
  set.seed(9101)
# Benign cell: auxiliaries present, normal errors, logit
  df <- make_iid_nmar(n = 600, alpha = 0.4, seed = 9101)
  resp_idx <- which(!is.na(df$Y_miss))
  respondent_data <- df[resp_idx, , drop = FALSE]
  respondent_weights <- rep(1, nrow(respondent_data))
  N_pop <- nrow(df)

# Build response model matrix with (Intercept) + Y_miss (and optional Z if needed)
  response_fml <- stats::as.formula("..delta.. ~ Y_miss")
  delta_name <- "..delta.."; respondent_data[[delta_name]] <- 1L
  Z_un <- model.matrix(update(response_fml, NULL ~ .), data = respondent_data)
  X_un <- model.matrix(~ X - 1, data = respondent_data)
  mu_x_un <- colMeans(model.matrix(~ X - 1, data = df))

  sc <- NMAR:::validate_and_apply_nmar_scaling(
    standardize = TRUE, has_aux = TRUE,
    response_model_matrix_unscaled = Z_un,
    auxiliary_matrix_unscaled = X_un,
    mu_x_unscaled = mu_x_un,
    weights = respondent_weights
  )
  Xbeta <- sc$response_model_matrix_scaled
  Xaux <- sc$auxiliary_matrix_scaled
  mu_x <- sc$mu_x_scaled
  K_beta <- ncol(Xbeta); K_aux <- ncol(Xaux)
  fam <- NMAR:::logit_family()
  eq <- NMAR:::el_build_equation_system(fam, Xbeta, Xaux, respondent_weights, N_pop, sum(respondent_weights), mu_x)
  Aj <- NMAR:::build_el_jacobian(fam, Xbeta, Xaux, respondent_weights, N_pop, sum(respondent_weights), mu_x)

# Solve for theta
  init <- c(rep(0, K_beta), 0, rep(0, K_aux))
  sol <- nleqslv::nleqslv(x = init, fn = eq, jac = Aj, method = "Newton",
                          control = list(ftol = 1e-8, xtol = 1e-8, maxit = 100), global = "dbldog")
  expect_lte(sol$termcd, 2)
  theta <- sol$x
  A <- Aj(theta)
  eta <- as.vector(Xbeta %*% theta[1:K_beta])
  w <- fam$linkinv(eta)
  W <- stats::plogis(theta[K_beta + 1])
  lamW <- ((N_pop / sum(respondent_weights)) - 1) / (1 - W)
  denom <- 1 + lamW * (w - W) + as.vector(sweep(Xaux, 2, mu_x, "-") %*% theta[(K_beta + 2):length(theta)])
  denom <- pmax(as.numeric(denom), 1e-8)

# Score contributions and B (respondent totals, centered)
  U <- NMAR:::el_compute_score_contrib(fam, Xbeta, Xaux, mu_x, eta, w, W, denom, lamW, respondent_weights, df)
  Uc <- scale(U, center = TRUE, scale = FALSE)
  Bm <- crossprod(Uc)

# Σ via two linear solves with correct right solve
  Xsolve <- solve(A, Bm)
  Sig <- Xsolve %*% solve(t(A))
  Sig <- 0.5 * (Sig + t(Sig))

# Mean identity: x' B x vs g' Σ g
  Xc <- sweep(Xaux, 2, mu_x, "-")
  g <- NMAR:::el_grad_g_analytic(fam, X_beta = Xbeta, Xc = Xc, mu_x_scaled = mu_x,
                                 respondent_weights = respondent_weights, eta_i_hat = eta, w_i_hat = w,
                                 W_hat = W, denominator_hat = denom, lambda_W_hat = lamW,
                                 outcome_vec = respondent_data$Y_miss, n_resp_weighted = sum(respondent_weights), N_pop = N_pop)
  xvec <- solve(t(A), g)
  var_two <- as.numeric(t(xvec) %*% Bm %*% xvec)
  var_sig <- as.numeric(t(g) %*% Sig %*% g)
  expect_true(is.finite(var_two) && is.finite(var_sig))
  expect_lt(abs(var_two - var_sig) / max(1e-8, abs(var_sig)), 1e-8)

# Beta identity: Vb_two vs Σ[beta,beta]
  Ebeta <- matrix(0, nrow = ncol(A), ncol = K_beta)
  Ebeta[seq_len(K_beta), ] <- diag(K_beta)
  Xbeta_solve <- solve(t(A), Ebeta)
  Vb_two <- crossprod(Xbeta_solve, Bm %*% Xbeta_solve)
  Vb_two <- 0.5 * (Vb_two + t(Vb_two))
  Vb_sig <- Sig[1:K_beta, 1:K_beta, drop = FALSE]
  Vb_sig <- 0.5 * (Vb_sig + t(Vb_sig))
  denomM <- max(1e-8, sqrt(sum(Vb_sig^2)))
  expect_lt(sqrt(sum((Vb_two - Vb_sig)^2)) / denomM, 1e-8)

# PSD checks
  min_eig <- function(M) suppressWarnings(min(eigen(M, symmetric = TRUE, only.values = TRUE)$values))
  expect_gte(min_eig(0.5 * (Bm + t(Bm))), -1e-10)
  expect_gte(min_eig(Vb_two), -1e-10)
  expect_gte(min_eig(Sig), -1e-10)
})
