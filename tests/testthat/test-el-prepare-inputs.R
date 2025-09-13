test_that("prepare_el_inputs creates unique delta var name when colliding", {
  df <- data.frame(`..nmar_delta..` = 1:5, Y = c(1, 2, NA, 4, 5), X = rnorm(5))
  res <- nmar:::prepare_el_inputs(Y ~ X, df, NULL)
  fml <- res$formula_list$response
  lhs <- all.vars(fml)[1]
  expect_true(lhs != "..nmar_delta..")
  expect_true(lhs %in% names(res$data))
})
