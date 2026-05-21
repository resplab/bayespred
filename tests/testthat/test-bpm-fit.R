test_that("bpmfit() returns a bpmfit object with expected top-level fields", {
  fit <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  expect_s3_class(fit, "bpmfit")
  expected <- c("posterior", "prior", "formula", "call", "model", "fit_method")
  expect_true(all(expected %in% names(fit)))
})

test_that("bpmfit()$posterior is a bpm object with expected fields", {
  fit  <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  post <- fit$posterior
  expect_s3_class(post, "bpm")
  expected_post <- c("coefficients", "vcov", "family", "terms", "contrasts", "xlevels")
  expect_true(all(expected_post %in% names(post)))
})

test_that("posterior() returns the bpm object from bpmfit", {
  fit  <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  post <- posterior(fit)
  expect_s3_class(post, "bpm")
  expect_identical(post, fit$posterior)
})

test_that("flat prior: point estimates match glm()", {
  fit_bpm <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  fit_glm <- glm(y ~ x1 + x2, data = toy, family = binomial())
  expect_equal(coef(fit_bpm), coef(fit_glm), tolerance = 1e-6)
})

test_that("flat prior: vcov matches glm()", {
  fit_bpm <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  fit_glm <- glm(y ~ x1 + x2, data = toy, family = binomial())
  expect_equal(vcov(fit_bpm), vcov(fit_glm), tolerance = 1e-5)
})

test_that("jeffreys prior: converges and returns named coefs", {
  fit <- bpmfit(y ~ x1 + x2, data = toy, prior = jeffreys())
  expect_s3_class(fit, "bpmfit")
  expect_false(anyNA(coef(fit)))
  expect_named(coef(fit), c("(Intercept)", "x1", "x2"))
})

test_that("log_f prior: converges and returns named coefs", {
  fit <- bpmfit(y ~ x1 + x2, data = toy, prior = log_f(m = 2))
  expect_s3_class(fit, "bpmfit")
  expect_false(anyNA(coef(fit)))
})

test_that("log_f prior: larger m shrinks slopes more toward 0", {
  fit_m1 <- bpmfit(y ~ x1 + x2, data = toy, prior = log_f(m = 1))
  fit_m5 <- bpmfit(y ~ x1 + x2, data = toy, prior = log_f(m = 5))
  expect_true(sum(coef(fit_m5)[-1]^2) <= sum(coef(fit_m1)[-1]^2))
})

test_that("bridge prior: converges and returns named coefs", {
  fit <- bpmfit(y ~ x1 + x2, data = toy, prior = bridge())
  expect_s3_class(fit, "bpmfit")
  expect_false(anyNA(coef(fit)))
  expect_named(coef(fit), c("(Intercept)", "x1", "x2"))
})

test_that("vcov is symmetric and positive definite for all priors", {
  priors <- list(flat(), jeffreys(), log_f(m = 2), bridge())
  for (pr in priors) {
    fit <- bpmfit(y ~ x1 + x2, data = toy, prior = pr)
    V   <- vcov(fit)
    expect_equal(V, t(V), tolerance = 1e-10)
    expect_true(all(eigen(V, only.values = TRUE)$values > 0))
  }
})

test_that("bpmfit() errors on unsupported family", {
  expect_error(bpmfit(y ~ x1, data = toy, family = gaussian()),       "logit")
  expect_error(bpmfit(y ~ x1, data = toy, family = binomial("probit")), "logit")
})

test_that("bpmfit() errors when prior is not a bpm_prior", {
  expect_error(bpmfit(y ~ x1, data = toy, prior = list(type = "flat")), "`prior`")
})

test_that("model = FALSE omits the model frame", {
  fit <- bpmfit(y ~ x1, data = toy, prior = flat(), model = FALSE)
  expect_null(fit$model)
})

test_that("fit_method tag is set correctly for each prior", {
  expect_equal(bpmfit(y ~ x1, data = toy, prior = flat())$fit_method,     "glm")
  expect_equal(bpmfit(y ~ x1, data = toy, prior = jeffreys())$fit_method, "brglmFit")
  expect_equal(bpmfit(y ~ x1, data = toy, prior = log_f())$fit_method,    "data_augmentation")
  expect_equal(bpmfit(y ~ x1, data = toy, prior = bridge())$fit_method,   "mgcv_ridge")
})

test_that("coef() and vcov() on bpm object match those on bpmfit object", {
  fit  <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  post <- posterior(fit)
  expect_equal(coef(post),  coef(fit))
  expect_equal(vcov(post),  vcov(fit))
})

test_that("predict.bpm() gives same result as predict.bpmfit() for new data", {
  fit  <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  post <- posterior(fit)
  expect_equal(predict(post, new1), predict(fit, new1))
})

test_that("predict.bpm() errors when newdata is missing", {
  post <- posterior(bpmfit(y ~ x1, data = toy, prior = flat()))
  expect_error(predict(post), "newdata")
})
