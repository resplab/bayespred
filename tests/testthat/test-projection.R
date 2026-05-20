test_that("projpred = TRUE stores projection with correct names", {
  fit <- bpm(y ~ x1 + x2, data = toy, prior = flat(), projpred = TRUE)
  expect_false(is.null(fit$projection))
  expect_named(fit$projection$coefficients, c("(Intercept)", "x1", "x2"))
  expect_equal(fit$projection$fit_method, "self_projection")
})

test_that("coef(type = 'projection') returns projected coefs", {
  fit <- bpm(y ~ x1 + x2, data = toy, prior = flat(), projpred = TRUE)
  expect_equal(coef(fit, type = "projection"), fit$projection$coefficients)
})

test_that("coef(type = 'main') returns main coefs", {
  fit <- bpm(y ~ x1 + x2, data = toy, prior = flat(), projpred = TRUE)
  expect_equal(coef(fit, type = "main"), fit$coefficients)
})

test_that("coef(type = 'projection') errors with helpful message when no projection", {
  fit <- bpm(y ~ x1, data = toy, prior = flat())
  expect_error(coef(fit, type = "projection"), "add_projection")
})

test_that("projected coefs are close to main coefs with large n, flat prior", {
  # With large n and no shrinkage, PM ~= PE, so self-projection ~= main coefs.
  fit <- bpm(y ~ x1 + x2, data = toy, prior = flat(), projpred = TRUE)
  expect_equal(coef(fit, type = "projection"), coef(fit), tolerance = 0.15)
})

test_that("pm_proj is close to pm with large n", {
  fit    <- bpm(y ~ x1 + x2, data = toy, prior = flat(), projpred = TRUE)
  p_pm   <- predict(fit, new1, method = "pm")
  p_proj <- predict(fit, new1, method = "pm_proj")
  expect_equal(p_proj, p_pm, tolerance = 0.03)
})
