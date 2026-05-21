test_that("project_pm() returns a bpmproj_pm object", {
  fit  <- bpm(y ~ x1 + x2, data = toy, prior = flat())
  proj <- project_pm(fit)
  expect_s3_class(proj, "bpmproj_pm")
})

test_that("project_pm() self-projection has correct coefficient names", {
  fit  <- bpm(y ~ x1 + x2, data = toy, prior = flat())
  proj <- project_pm(fit)
  expect_named(coef(proj), c("(Intercept)", "x1", "x2"))
})

test_that("project_pm() coefficients close to main coefs with large n, flat prior", {
  fit  <- bpm(y ~ x1 + x2, data = toy, prior = flat())
  proj <- project_pm(fit)
  expect_equal(coef(proj), coef(fit), tolerance = 0.15)
})

test_that("project_pm() requires data when model = FALSE", {
  fit <- bpm(y ~ x1 + x2, data = toy, prior = flat(), model = FALSE)
  expect_error(project_pm(fit), "Supply `data`|model = TRUE")
  proj <- project_pm(fit, data = toy)
  expect_s3_class(proj, "bpmproj_pm")
})

test_that("predict.bpmproj_pm type='response' returns values in (0, 1)", {
  proj <- project_pm(bpm(y ~ x1 + x2, data = toy, prior = flat()))
  p    <- predict(proj, new1)
  expect_true(p > 0 && p < 1)
})

test_that("predict.bpmproj_pm type='link' returns a finite numeric", {
  proj <- project_pm(bpm(y ~ x1 + x2, data = toy, prior = flat()))
  eta  <- predict(proj, new1, type = "link")
  expect_true(is.numeric(eta) && is.finite(eta))
})

test_that("predict.bpmproj_pm is close to predict(fit, method='pm') with large n", {
  fit    <- bpm(y ~ x1 + x2, data = toy, prior = flat())
  proj   <- project_pm(fit)
  p_pm   <- predict(fit, new1, method = "pm")
  p_proj <- predict(proj, new1)
  expect_equal(p_proj, p_pm, tolerance = 0.03)
})

test_that("predict.bpmproj_pm handles multiple rows", {
  proj   <- project_pm(bpm(y ~ x1 + x2, data = toy, prior = flat()))
  newdf  <- data.frame(x1 = c(-1, 0, 1), x2 = c(0, 0, 0))
  p      <- predict(proj, newdf)
  expect_length(p, 3L)
  expect_true(all(p > 0 & p < 1))
})

test_that("custom projection onto fewer predictors works", {
  fit  <- bpm(y ~ x1 + x2, data = toy, prior = flat())
  proj <- project_pm(fit, formula = ~ x1, data = toy)
  expect_s3_class(proj, "bpmproj_pm")
  expect_named(coef(proj), c("(Intercept)", "x1"))
  p <- predict(proj, new1)
  expect_true(p > 0 && p < 1)
})

test_that("print.bpmproj_pm runs without error", {
  proj <- project_pm(bpm(y ~ x1 + x2, data = toy, prior = flat()))
  expect_output(print(proj), "PM-Projected")
})
