test_that("EL engine runs for survey.design and returns CI (skip if survey missing)", {
  skip_if_not_installed("survey")
  library(survey)
  data(api)
  set.seed(3103)
  apiclus1_df <- apiclus1
  apiclus1_df$api00_miss <- apiclus1_df$api00
  y_std <- scale(apiclus1_df$api00)[, 1]
  prob <- plogis(-0.5 + 0.4 * y_std + 0.2 * scale(apiclus1_df$ell)[, 1])
  miss_idx <- runif(nrow(apiclus1_df)) > prob
  apiclus1_df$api00_miss[miss_idx] <- NA

  dclus1 <- svydesign(id = ~dnum, weights = ~pw, data = apiclus1_df, fpc = ~fpc)
  pop_mean_ell <- mean(apiclus1_df$ell)

  eng <- el_engine(variance_method = "delta", auxiliary_means = c(ell = pop_mean_ell))
  fml <- api00_miss ~ ell
  res <- nmar(formula = fml, data = dclus1, engine = eng)

  expect_s3_class(res, "nmar_result_el")
  expect_true(isTRUE(res$converged))
  expect_true(is.numeric(res$se) && res$se >= 0)
  ci <- nmar:::confint.nmar_result_el(res)
  expect_true(is.matrix(ci) && nrow(ci) == 1)
  expect_equal(rownames(ci), "api00_miss")
})
