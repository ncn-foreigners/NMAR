test_that("el_prepare_inputs handles non-default contrasts for factor auxiliaries", {
  set.seed(456)
  n <- 60
  F <- factor(sample(letters[1:4], n, replace = TRUE))
  Y <- 1 + rnorm(n)
  R <- runif(n) < 0.7
  df <- data.frame(Y_miss = Y, F = F)
  df$Y_miss[!R] <- NA_real_

  contrasts(df$F) <- stats::contr.sum(nlevels(df$F))

  spec <- el_prepare_inputs(Y_miss ~ F, df)
  aux <- spec$aux_design_full
  expect_true(is.matrix(aux))
  expect_equal(ncol(aux), nlevels(df$F) - 1)

  eng <- el_engine(auxiliary_means = NULL, variance_method = "none")
  fit <- nmar(Y_miss ~ F, data = df, engine = eng)
  expect_s3_class(fit, "nmar_result_el")
})

test_that("auxiliary factors use k-1 coding with sum contrasts and explicit +0", {
  set.seed(789)
  n <- 80
  F <- factor(sample(letters[1:3], n, replace = TRUE))
  Y <- rnorm(n)
  R <- runif(n) < 0.7
  df <- data.frame(Y_miss = Y, F = F)
  df$Y_miss[!R] <- NA_real_

# Sum contrasts
  contrasts(df$F) <- stats::contr.sum(nlevels(df$F))

  spec1 <- el_prepare_inputs(Y_miss ~ F, df)
  spec2 <- el_prepare_inputs(Y_miss ~ 0 + F, df)

  aux1 <- spec1$aux_design_full
  aux2 <- spec2$aux_design_full

  expect_true(is.matrix(aux1))
  expect_true(is.matrix(aux2))
  expect_equal(ncol(aux1), nlevels(df$F) - 1)
  expect_equal(ncol(aux2), nlevels(df$F) - 1)
  expect_identical(colnames(aux1), colnames(aux2))
})
