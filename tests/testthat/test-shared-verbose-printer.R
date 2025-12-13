test_that("create_verboser validates trace_level and is silent at 0", {
  expect_error(NMAR:::create_verboser(99), "trace_level", fixed = FALSE)

  v0 <- NMAR:::create_verboser(0)
  out0 <- capture.output(v0("hello", level = 1, type = "info"))
  expect_length(out0, 0)
})

test_that("create_verboser prints only at or below trace_level", {
  v1 <- NMAR:::create_verboser(1)

  out1 <- capture.output(v1("hello", level = 1, type = "info"))
  expect_true(any(grepl("\\[INFO\\]", out1)))
  expect_true(any(grepl("hello", out1)))

  out2 <- capture.output(v1("nope", level = 2, type = "info"))
  expect_length(out2, 0)
})

test_that("create_verboser prints object summaries across types", {
  v3 <- NMAR:::create_verboser(3)

  out_num <- capture.output(v3("x", level = 1, type = "detail", obj = 1.234))
  expect_true(any(grepl("Value:", out_num)))

  out_short <- capture.output(v3("x", level = 1, type = "detail", obj = 1:3, max_print = 10))
  expect_true(any(grepl("Values:", out_short)))

  out_long <- capture.output(v3("x", level = 1, type = "detail", obj = 1:100, max_print = 10))
  expect_true(any(grepl("Length:", out_long)))
  expect_true(any(grepl("Range:", out_long)))
  expect_true(any(grepl("Mean:", out_long)))

  out_mat_small <- capture.output(v3("x", level = 1, type = "detail", obj = matrix(1:4, nrow = 2), max_print = 10))
  expect_true(any(grepl("Dimensions:", out_mat_small)))

  out_mat_big <- capture.output(v3("x", level = 1, type = "detail", obj = matrix(1:400, nrow = 20), max_print = 10))
  expect_true(any(grepl("Dimensions:", out_mat_big)))
  expect_true(any(grepl("Range:", out_mat_big)))

  out_other <- capture.output(
    v3(
      "x",
      level = 1,
      type = "detail",
      obj = structure(list(a = 1), class = "nmar_test_obj")
    )
  )
  expect_true(any(grepl("\\$a", out_other)))
})
