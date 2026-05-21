test_that("bpmproject() returns a bpmproj_pm object for type='pm'", {
  fit  <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  proj <- bpmproject(fit)
  expect_s3_class(proj, "bpmproj_pm")
})

test_that("bpmproject() self-projection has correct coefficient names", {
  fit  <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  proj <- bpmproject(fit)
  expect_named(coef(proj), c("(Intercept)", "x1", "x2"))
})

test_that("bpmproject() coefficients close to main coefs with large n, flat prior", {
  fit  <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  proj <- bpmproject(fit)
  expect_equal(coef(proj), coef(fit), tolerance = 0.15)
})

test_that("bpmproject() requires data when model = FALSE", {
  fit <- bpmfit(y ~ x1 + x2, data = toy, prior = flat(), model = FALSE)
  expect_error(bpmproject(fit), "Supply `data`|model = TRUE")
  proj <- bpmproject(fit, data = toy)
  expect_s3_class(proj, "bpmproj_pm")
})

test_that("bpmproject(type='full') errors with not-yet-implemented message", {
  fit <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  expect_error(bpmproject(fit, type = "full"), "not yet implemented")
})

test_that("predict on bpmproject result (type='response') returns values in (0, 1)", {
  proj <- bpmproject(bpmfit(y ~ x1 + x2, data = toy, prior = flat()))
  p    <- predict(proj, new1)
  expect_true(p > 0 && p < 1)
})

test_that("predict on bpmproject result (type='link') returns a finite numeric", {
  proj <- bpmproject(bpmfit(y ~ x1 + x2, data = toy, prior = flat()))
  eta  <- predict(proj, new1, type = "link")
  expect_true(is.numeric(eta) && is.finite(eta))
})

test_that("predict on bpmproject result is close to predict(fit, method='pm') with large n", {
  fit    <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  proj   <- bpmproject(fit)
  p_pm   <- predict(fit, new1, method = "pm")
  p_proj <- predict(proj, new1)
  expect_equal(p_proj, p_pm, tolerance = 0.03)
})

test_that("predict on bpmproject result handles multiple rows", {
  proj  <- bpmproject(bpmfit(y ~ x1 + x2, data = toy, prior = flat()))
  newdf <- data.frame(x1 = c(-1, 0, 1), x2 = c(0, 0, 0))
  p     <- predict(proj, newdf)
  expect_length(p, 3L)
  expect_true(all(p > 0 & p < 1))
})

test_that("custom projection onto fewer predictors works", {
  fit  <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  proj <- bpmproject(fit, formula = ~ x1, data = toy)
  expect_s3_class(proj, "bpmproj_pm")
  expect_named(coef(proj), c("(Intercept)", "x1"))
  p <- predict(proj, new1)
  expect_true(p > 0 && p < 1)
})

test_that("print on bpmproject result runs without error", {
  proj <- bpmproject(bpmfit(y ~ x1 + x2, data = toy, prior = flat()))
  expect_output(print(proj), "PM-Projected")
})
