skip_if_not_installed("survey")

test_that("survey designs reuse EL prep workflow", {
  set.seed(100)
  df <- data.frame(
    Y_miss = c(1, 2, NA, 4),
    X = rnorm(4),
    w = c(1, 2, 1, 1)
  )
  design <- survey::svydesign(ids = ~1, data = df, weights = ~w)

  prep <- el_prepare_inputs(Y_miss ~ X, df)
  design_svy <- el_prepare_inputs(Y_miss ~ X, design$variables)

  expect_equal(prep$missingness_design, design_svy$missingness_design)
  expect_equal(prep$aux_design_full[prep$respondent_mask, , drop = FALSE],
               design_svy$aux_design_full[design_svy$respondent_mask, , drop = FALSE])
  expect_equal(prep$respondent_mask, design_svy$respondent_mask)
})

test_that("survey strata augmentation appends dummies and implied means", {
  skip_if_not_installed("survey")
  set.seed(123)
  df <- data.frame(
    Y_miss = c(1, NA, 2, NA, 3, NA),
    X = rnorm(6),
    strata = factor(c("A", "A", "B", "B", "B", "B"))
  )
  w <- c(2, 2, 3, 3, 3, 3)
  design <- survey::svydesign(ids = ~1, strata = ~strata, data = df, weights = ~w)

  eng <- el_engine(auxiliary_means = c(X = 0), variance_method = "none", strata_augmentation = TRUE)
  fit <- nmar(Y_miss ~ X, data = design, engine = eng)

# augmented auxiliary matrix should contain strata dummies
  aux_mat <- fit$diagnostics$auxiliary_matrix
  expect_true(any(grepl("^strata_", colnames(aux_mat))))

# implied strata means should match N_h / N_pop for non-reference levels
  N_pop <- fit$sample$n_total
  strata_fac <- design$variables$strata
  W_h <- tapply(w, strata_fac, sum) / N_pop
  strata_means_fit <- fit$diagnostics$auxiliary_means[grepl("^strata_", names(fit$diagnostics$auxiliary_means))]
  expect_equal(unname(strata_means_fit), unname(W_h[setdiff(levels(strata_fac), levels(strata_fac)[1])]), ignore_attr = TRUE)
})

test_that("survey strata augmentation works with sample-implied auxiliary means", {
  skip_if_not_installed("survey")
  set.seed(456)

  df <- make_iid_nmar(n = 200, alpha = 0.4, include_z = FALSE, seed = 2L)
  df$strata <- factor(rep(letters[1:4], length.out = nrow(df)))
  df$w <- runif(nrow(df), 0.5, 2)

  design <- survey::svydesign(ids = ~1, strata = ~strata, weights = ~w, data = df)
  eng <- el_engine(variance_method = "none", strata_augmentation = TRUE)

  fit <- nmar(Y_miss ~ X, data = design, engine = eng)

  expect_s3_class(fit, "nmar_result_el")
  expect_true(isTRUE(fit$converged))

  aux_means <- fit$diagnostics$auxiliary_means
  expect_true("X" %in% names(aux_means))
  expect_true(any(grepl("^strata_", names(aux_means))))
})

test_that("survey prep stores delta column and uses rescaled weights", {
  skip_if_not_installed("survey")
  set.seed(200)
  df <- data.frame(
    Y_miss = c(rnorm(5), NA),
    X = rnorm(6),
    w = c(1, 2, 3, 4, 5, 6)
  )
  des <- survey::svydesign(ids = ~1, weights = ~w, data = df)

  n_total <- 2 * sum(weights(des))
  res <- suppressWarnings(el.survey.design(
    data = des,
    formula = Y_miss ~ X,
    auxiliary_means = c(X = 0),
    n_total = n_total,
    variance_method = "none"
  ))
  des_after <- res$sample$design
  expect_true("..nmar_delta.." %in% names(des_after$variables))
  expect_equal(sum(stats::weights(res, scale = "population")), n_total)
})

test_that("respondents-only survey with strata augmentation warns on shares", {
  skip_if_not_installed("survey")
  df <- data.frame(
    Y_miss = c(1, 2, 3),
    strata = factor(c("A", "B", "B")),
    w = c(1, 1, 1)
  )
  des <- survey::svydesign(ids = ~1, strata = ~strata, data = df, weights = ~w)
  eng <- el_engine(variance_method = "none", strata_augmentation = TRUE, auxiliary_means = c(), n_total = sum(weights(des)))
  expect_warning(
    nmar(Y_miss ~ 1, data = des, engine = eng),
    regexp = "strata augmentation with respondents-only",
    fixed = FALSE
  )
})

test_that("el_prepare_inputs carries survey metadata and totals", {
  skip_if_not_installed("survey")
  set.seed(321)
  df <- data.frame(
    Y_miss = c(rnorm(4), NA, 1),
    X = rnorm(6),
    w = c(1, 2, 1, 1, 3, 4)
  )
  des <- survey::svydesign(ids = ~1, weights = ~w, data = df)
  n_total <- sum(weights(des)) * 1.5
  spec <- el_prepare_inputs(
    formula = Y_miss ~ X,
    data = des$variables,
    weights = as.numeric(weights(des)),
    n_total = n_total,
    design_object = des
  )
  expect_s3_class(spec$analysis_data, "survey.design")
  expect_true("..nmar_delta.." %in% names(spec$analysis_data$variables))
  expect_equal(spec$N_pop, n_total)
  expect_equal(sum(spec$respondent_weights), sum(spec$respondent_mask * weights(des)))
})

test_that("survey strata extraction supports interactions", {
  skip_if_not_installed("survey")
  set.seed(404)
  df <- data.frame(
    Y_miss = c(1, 2, 2, NA, 3),
    A = factor(c("a", "a", "b", "b", "b")),
    B = factor(c("x", "y", "x", "y", "y")),
    w = c(1, 2, 3, 4, 5)
  )
  des <- survey::svydesign(ids = ~1, strata = ~ A + B, data = df, weights = ~w)
  eng <- el_engine(variance_method = "none")
  fit <- nmar(Y_miss ~ 1, data = des, engine = eng)
# Interaction of A and B should be used for augmentation; expect strata_ levels
  aug_means <- fit$diagnostics$auxiliary_means
  strata_terms <- names(aug_means)[grepl("^strata_", names(aug_means))]
  expect_true(length(strata_terms) > 0)
})
