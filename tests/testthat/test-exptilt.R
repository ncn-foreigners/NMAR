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

test_that("exptilt returns OK data for correct input", {
  formula_ok <- Y ~ x1 + x2
  exptilt_config_ok <- exptilt_engine(
    family = 'probit',
    y_dens = 'normal',
    tol_value = 0.01,
    min_iter = 1,
    max_iter = 5
  )

  res <- nmar(formula = formula_ok, data = data_OK_full, engine = exptilt_config_ok)

  expect_s3_class(res, "nmar_result")
  expect_type(res, "list")
  expect_type(res$y_hat, "double")
  expect_type(res$se, "double")
  expect_length(res$y_hat, 1)
  expect_length(res$se, 1)
  expect_true(is.finite(res$y_hat))
  expect_true(is.finite(res$se))
  expect_gt(res$se, 0)
})

test_that("exptilt returns error for faulty data", {
  formula_ok <- Y ~ x1 + x2
  exptilt_config_ok <- exptilt_engine(
    family = 'probit',
    y_dens = 'normal',
    tol_value = 0.01,
    min_iter = 1,
    max_iter = 5
  )

  data_constant_y <- data_OK
  data_constant_y$Y <- 42
  expect_error(nmar(formula = formula_ok, data = data_constant_y, engine = exptilt_config_ok))

  data_nan <- data_OK
  data_nan$x1[1] <- NaN
  expect_error(nmar(formula = formula_ok, data = data_nan, engine = exptilt_config_ok))

  data_inf <- data_OK
  data_inf$x1[1] <- Inf
  expect_error(nmar(formula = formula_ok, data = data_inf, engine = exptilt_config_ok))

  data_char <- data_OK
  data_char$x1[1] <- "invalid_value"
  expect_error(nmar(formula = formula_ok, data = data_char, engine = exptilt_config_ok))
})

test_that("exptilt returns error for faulty formulas", {
  exptilt_config_ok <- exptilt_engine(
    family = 'probit',
    y_dens = 'normal',
    tol_value = 0.01,
    min_iter = 1,
    max_iter = 5
  )

  expect_error(nmar(formula = not_existing_col ~ x1 + x2, data = data_OK, engine = exptilt_config_ok))
  expect_error(nmar(formula = Y ~ x1 + not_existing_col, data = data_OK, engine = exptilt_config_ok))
  expect_error(nmar(formula = Y ~ x1 + x2, data = data_OK, engine = exptilt_config_ok, response_predictors = "not_existing_col"))
})

test_that("exptilt returns error for bad config", {
  formula_ok <- Y ~ x1 + x2

  expect_error(nmar(
    formula = formula_ok,
    data = data_OK,
    engine = exptilt_engine(family = 'invalid_type', y_dens = 'normal')
  ))

  expect_error(nmar(
    formula = formula_ok,
    data = data_OK,
    engine = exptilt_engine(family = 'probit', y_dens = 'invalid_distribution')
  ))

  expect_error(nmar(
    formula = formula_ok,
    data = data_OK,
    engine = exptilt_engine(family = 'probit', y_dens = 'normal', tol_value = -0.01)
  ))

  expect_error(nmar(
    formula = formula_ok,
    data = data_OK,
    engine = exptilt_engine(family = 'probit', y_dens = 'normal', min_iter = 10, max_iter = 5)
  ))
})

test_that("exptilt returns error for empty data", {
  formula_ok <- Y ~ x1 + x2
  exptilt_config_ok <- exptilt_engine(
    family = 'probit',
    y_dens = 'normal',
    tol_value = 0.01
  )

  data_empty <- data_OK[0, ]
  expect_error(nmar(formula = formula_ok, data = data_empty, engine = exptilt_config_ok))
})

test_that("exptilt returns error for missing values in covariates", {
  formula_ok <- Y ~ x1 + x2
  exptilt_config_ok <- exptilt_engine(
    family = 'probit',
    y_dens = 'normal',
    tol_value = 0.01
  )

  data_na_covariates <- data_OK
  data_na_covariates$x1[1] <- NA
  expect_error(nmar(formula = formula_ok, data = data_na_covariates, engine = exptilt_config_ok))
})
