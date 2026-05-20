test_that("prior constructors return correct S3 classes", {
  expect_s3_class(flat(),        c("prior_flat",     "bpm_prior"))
  expect_s3_class(jeffreys(),    c("prior_jeffreys", "bpm_prior"))
  expect_s3_class(log_f(),       c("prior_logf",     "bpm_prior"))
  expect_s3_class(log_f(m = 1), c("prior_logf",     "bpm_prior"))
  expect_s3_class(bridge(),      c("prior_bridge",   "bpm_prior"))
})

test_that("log_f validates m", {
  expect_error(log_f(m = 0),   "`m`")
  expect_error(log_f(m = -1),  "`m`")
  expect_error(log_f(m = "a"), "`m`")
  expect_no_error(log_f(m = 0.5))
  expect_no_error(log_f(m = 10))
})

test_that("log_f stores m on the object", {
  expect_equal(log_f(m = 3)$m, 3)
  expect_equal(log_f()$m, 2)
})

test_that("print.bpm_prior runs without error", {
  expect_output(print(flat()),     "flat")
  expect_output(print(jeffreys()), "Jeffreys")
  expect_output(print(log_f(2)),   "log-F")
  expect_output(print(bridge()),   "Ridge")
})
