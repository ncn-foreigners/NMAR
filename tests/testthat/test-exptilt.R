generate_test_data <- function(n_rows=400,n_cols=2,case=2,x_var=0.5,eps_var=0.5,a=0.8,b=-0.2) {
  X <- as.data.frame(replicate(n_cols, rnorm(n_rows, 0, x_var)^1))
  colnames(X) <- paste0("x", 1:n_cols)

  coefs <- runif(n_cols, 0.95, 1.05)
  eps <- rnorm(n_rows, 0, eps_var)
  if (case==1){
    X$Y <- as.vector(-1 + as.matrix(X) %*% coefs + eps)
  }
  else if(case==2){
    X$Y <- -2 + 0.5 * exp(as.matrix(X) %*% coefs) + eps
  }
  else if(case==3){
    X$Y <- -1 + sin(2* as.matrix(X) %*% coefs) + eps
  }
  else if(case==4){
    X$Y <- -1 + 0.4 * as.matrix(X)^3 %*% coefs + eps
  }

  X <- X[order(X$Y), ]
  Y_original = X$Y

  pi_obs <- 1 / (1 + exp(-(a + b * X$Y)))

  mask <- runif(nrow(X)) > pi_obs
  mask[1]=TRUE
  mask[2]=FALSE
  X$Y[mask] <- NA
  return(list(X = X, Y_original = Y_original))
}

data_OK <- generate_test_data(n_rows = 5, n_cols = 2, case = 1)$X
data_OK_full <- generate_test_data(n_rows = 100, n_cols = 2, case = 1)$X

formula_faultys <- c(
  list(
    outcome = ~ Y,
    covariates_outcome = ~ x1 + x2,
    covariates_missingness = ~ x1 + x2
  ),
  list(
    outcome = ~ not_existing_col,
    covariates_outcome = ~ x1 + x2,
    covariates_missingness = ~ x1 + x2
  ),
  list(
    outcome = ~ Y,
    covariates_outcome = ~ x1 + not_existing_col,
    covariates_missingness = ~x2
  ),
  list(
    outcome = ~ Y,
    covariates_outcome = ~ x1 + x2,
    covariates_missingness = ~ not_existing_col
  ),
  list(
    outcome = ~ NULL,
    covariates_outcome = ~ NULL,
    covariates_missingness = ~ NULL
  ),
  list(
    outcome = ~ Y,
    covariates_outcome = ~ Y,
    covariates_missingness = ~ NULL
  )
)


test_that("exptilt returns OK data for correct input", {
  formula_ok <- list(
    outcome = ~ Y,
    covariates_outcome = ~ x1 + x2,
    covariates_missingness = ~ NULL
  )

  exptilt_config_ok <- exptilt(
    prob_model_type = 'probit',
    y_dens = 'normal',
    tol_value = 0.01,
    min_iter = 1,
    max_iter = 5
  )

  res <- nmar(formula = formula_ok, data = data_OK_full, engine = exptilt_config_ok)

  expect_s3_class(res, "nmar_result")
  expect_type(res, "list")

  expect_type(res$est_mean, "double")
  expect_type(res$est_var, "double")

  expect_length(res$est_mean, 1)
  expect_length(res$est_var, 1)

  expect_true(is.finite(res$est_mean))
  expect_true(is.finite(res$est_var))

  expect_gt(res$est_var, 0)
})


test_that("exptilt returns error for faulty data and OK formula", {
  formula_ok <- list(
    outcome = ~ Y,
    covariates_outcome = ~ x1 + x2,
    covariates_missingness = ~ NULL
  )

  exptilt_config_ok <- exptilt(
    prob_model_type = 'probit',
    y_dens = 'normal',
    tol_value = 0.01,
    min_iter = 1,
    max_iter = 5
  )

  # 1. No NA in Y
  data_constant_y <- data_OK
  data_constant_y$Y <- 42

  expect_error(nmar(formula = formula_ok, data = data_constant_y, engine = exptilt_config_ok))

  # 2. NaN in X
  data_nan <- data_OK
  data_nan$x1[1] <- NaN

  expect_error(nmar(formula = formula_ok, data = data_nan, engine = exptilt_config_ok))

  # 3. Inf in X
  data_inf <- data_OK
  data_inf$x1[1] <- Inf

  expect_error(nmar(formula = formula_ok, data = data_inf, engine = exptilt_config_ok))

  # 4. Non-numeric in X
  data_char <- data_OK
  data_char$x1[1] <- "invalid_value"

  expect_error(nmar(formula = formula_ok, data = data_char, engine = exptilt_config_ok))
})

# Tests for OK data but faulty formulas
test_that("exptilt returns error for faulty formulas with OK data", {
  exptilt_config_ok <- exptilt(
    prob_model_type = 'probit',
    y_dens = 'normal',
    tol_value = 0.01,
    min_iter = 1,
    max_iter = 5
  )


  for (i in seq_along(formula_faultys)) {
    expect_error(nmar(formula = formula_faultys[[i]], data = data_OK, engine = exptilt_config_ok))
  }
})

# Tests for bad config
test_that("exptilt returns error for faulty configurations", {
  formula_ok <- list(
    outcome = ~ Y,
    covariates_outcome = ~ x1 + x2,
    covariates_missingness = ~ NULL
  )

  # 1
  exptilt_faulty1 <- exptilt(
    prob_model_type = 'invalid_type',
    y_dens = 'normal',
    tol_value = 0.01
  )

  expect_error(nmar(formula = formula_ok, data = data_OK, engine = exptilt_faulty1))

  # 2
  exptilt_faulty2 <- exptilt(
    prob_model_type = 'probit',
    y_dens = 'invalid_distribution',
    tol_value = 0.01
  )

  expect_error(nmar(formula = formula_ok, data = data_OK, engine = exptilt_faulty2))

  # 3.
  exptilt_faulty3 <- exptilt(
    prob_model_type = 'probit',
    y_dens = 'normal',
    tol_value = -0.01
  )

  expect_error(nmar(formula = formula_ok, data = data_OK, engine = exptilt_faulty3))

  # 4.
  exptilt_faulty4 <- exptilt(
    prob_model_type = 'probit',
    y_dens = 'normal',
    min_iter = 10,
    max_iter = 5
  )

  expect_error(nmar(formula = formula_ok, data = data_OK, engine = exptilt_faulty4))
})

# Emptyness
test_that("exptilt returns error for empty data", {
  formula_ok <- list(
    outcome = ~ Y,
    covariates_outcome = ~ x1 + x2,
    covariates_missingness = ~ NULL
  )

  exptilt_config_ok <- exptilt(
    prob_model_type = 'probit',
    y_dens = 'normal',
    tol_value = 0.01
  )


  data_empty <- data_OK[0, ]

  expect_error(nmar(formula = formula_ok, data = data_empty, engine = exptilt_config_ok))
})


test_that("exptilt returns error for data with missing values in covariates", {
  formula_ok <- list(
    outcome = ~ Y,
    covariates_outcome = ~ x1 + x2,
    covariates_missingness = ~ NULL
  )

  exptilt_config_ok <- exptilt(
    prob_model_type = 'probit',
    y_dens = 'normal',
    tol_value = 0.01
  )


  data_na_covariates <- data_OK
  data_na_covariates$x1[1] <- NA

  expect_error(nmar(formula = formula_ok, data = data_na_covariates, engine = exptilt_config_ok))
})
