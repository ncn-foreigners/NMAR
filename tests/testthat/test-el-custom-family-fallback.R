test_that("EL runs with a custom family lacking second derivatives", {
  df <- make_iid_nmar(n = 250, alpha = 0.6, seed = 4401)

  custom_logit <- list(
    name = "logit",
    linkinv = stats::plogis,
    mu.eta = function(eta) {
      p <- stats::plogis(eta)
      p * (1 - p)
    },
    score_eta = function(eta, y) {
      p <- stats::plogis(eta)
      ifelse(y == 1, 1 - p, -p)
    }
  )

  eng <- el_engine(
    auxiliary_means = c(X = 0),
    variance_method = "none",
    standardize = TRUE,
    family = custom_logit
  )

  fit <- nmar(Y_miss ~ X, data = df, engine = eng, trace_level = 0)
  expect_s3_class(fit, "nmar_result_el")
  expect_true(isTRUE(fit$converged))
  expect_true(is.finite(as.numeric(fit$y_hat)))
})

test_that("analytic Jacobian is disabled when family lacks d2mu.deta2", {
  df <- make_iid_nmar(n = 120, alpha = 0.6, seed = 4402)
  spec <- NMAR:::el_prepare_inputs(
    formula = Y_miss ~ X,
    data = df,
    weights = NULL,
    n_total = NULL
  )

  resp_idx <- which(spec$respondent_mask)
  Z_un <- spec$missingness_design
  X_un <- spec$aux_design_full[spec$respondent_mask, , drop = FALSE]
  aux_means <- c(X = 0)
  sc <- NMAR:::validate_and_apply_nmar_scaling(TRUE, ncol(X_un) > 0, Z_un, X_un, aux_means)

  custom_logit <- list(
    name = "logit",
    linkinv = stats::plogis,
    mu.eta = function(eta) {
      p <- stats::plogis(eta)
      p * (1 - p)
    },
    score_eta = function(eta, y) {
      p <- stats::plogis(eta)
      ifelse(y == 1, 1 - p, -p)
    }
  )

  jac_fun <- NMAR:::el_build_jacobian(
    family = custom_logit,
    missingness_model_matrix = sc$response_model_matrix_scaled,
    auxiliary_matrix = sc$auxiliary_matrix_scaled,
    respondent_weights = rep(1, length(resp_idx)),
    N_pop = nrow(spec$analysis_data),
    n_resp_weighted = length(resp_idx),
    mu_x_scaled = sc$mu_x_scaled
  )
  expect_null(jac_fun)
})
