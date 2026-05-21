# ---- glm compatibility -------------------------------------------------------

test_that("type='response', method='pe' matches predict.glm(type='response')", {
  fit_bpm <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  fit_glm <- glm(y ~ x1 + x2, data = toy, family = binomial())
  p_bpm   <- predict(fit_bpm, new1, type = "response", method = "pe")
  p_glm   <- predict(fit_glm, new1, type = "response")
  expect_equal(p_bpm, unname(p_glm), tolerance = 1e-6)
})

test_that("type='link' matches predict.glm(type='link')", {
  fit_bpm <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  fit_glm <- glm(y ~ x1 + x2, data = toy, family = binomial())
  p_bpm   <- predict(fit_bpm, new1, type = "link")
  p_glm   <- predict(fit_glm, new1, type = "link")
  expect_equal(p_bpm, unname(p_glm), tolerance = 1e-6)
})

test_that("type='link', se.fit=TRUE matches predict.glm structure", {
  fit_bpm <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  fit_glm <- glm(y ~ x1 + x2, data = toy, family = binomial())
  out_bpm <- predict(fit_bpm, new1, type = "link", se.fit = TRUE)
  out_glm <- predict(fit_glm, new1, type = "link", se.fit = TRUE)
  expect_named(out_bpm, c("fit", "se.fit", "residual.scale"))
  expect_equal(out_bpm$fit,    unname(out_glm$fit),    tolerance = 1e-6)
  expect_equal(out_bpm$se.fit, unname(out_glm$se.fit), tolerance = 1e-5)
  expect_equal(out_bpm$residual.scale, 1)
})

test_that("type='response', se.fit=TRUE returns correct structure", {
  fit <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  out <- predict(fit, new1, type = "response", se.fit = TRUE)
  expect_named(out, c("fit", "se.fit", "residual.scale"))
  expect_true(out$fit > 0 && out$fit < 1)
  expect_true(out$se.fit > 0)
  expect_equal(out$residual.scale, 1)
})

test_that("type='response', se.fit=TRUE: se is delta-method approx (PE case)", {
  fit    <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  out    <- predict(fit, new1, type = "response", method = "pe", se.fit = TRUE)
  mu     <- out$fit
  # delta-method: se_response = mu*(1-mu) * se_link
  se_lnk <- predict(fit, new1, type = "link", se.fit = TRUE)$se.fit
  expect_equal(out$se.fit, mu * (1 - mu) * se_lnk, tolerance = 1e-8)
})

test_that("type='terms' raises a clear error", {
  fit <- bpmfit(y ~ x1, data = toy, prior = flat())
  expect_error(predict(fit, new1, type = "terms"), "not implemented")
})

test_that("se.fit and interval cannot be combined", {
  fit <- bpmfit(y ~ x1, data = toy, prior = flat())
  expect_error(predict(fit, new1, se.fit = TRUE, interval = 0.95), "cannot be used together")
})

# ---- type = "link" with interval ---------------------------------------------

test_that("type='link' with interval returns symmetric data frame with se.link", {
  fit <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  out <- predict(fit, new1, type = "link", interval = 0.95)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("fit", "lwr", "upr", "se.link"))
  # Symmetric around fit on the link scale
  expect_equal(out$fit - out$lwr, out$upr - out$fit, tolerance = 1e-10)
  expect_true(out$lwr < out$fit && out$fit < out$upr)
  # se.link is positive
  expect_true(out$se.link > 0)
  # lwr/upr are consistent with se.link
  z <- qnorm(0.975)
  expect_equal(out$lwr, out$fit - z * out$se.link, tolerance = 1e-10)
})

test_that("type='response' with interval has asymmetric bounds (logistic transform)", {
  fit  <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  out  <- predict(fit, new1, type = "response", interval = 0.95)
  # Asymmetric after plogis transform
  lower_gap <- out$fit - out$lwr
  upper_gap <- out$upr - out$fit
  expect_false(isTRUE(all.equal(lower_gap, upper_gap)))
})

# ---- Bayesian methods --------------------------------------------------------

test_that("all methods return values strictly in (0, 1)", {
  fit <- bpmfit(y ~ x1 + x2, data = toy, prior = log_f(m = 2))
  for (meth in c("pe", "pm", "pm_mackay")) {
    p <- predict(fit, new1, method = meth)
    expect_true(p > 0 && p < 1, label = paste("method =", meth))
  }
})

test_that("pm and pm_mackay agree closely for moderate uncertainty", {
  fit  <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  p_pm <- predict(fit, new1, method = "pm")
  p_mk <- predict(fit, new1, method = "pm_mackay")
  expect_equal(p_pm, p_mk, tolerance = 0.01)
})

test_that("method is ignored when type = 'link'", {
  fit <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  # All methods should give identical link predictions
  vals <- sapply(c("pe", "pm", "pm_mackay"),
                 function(m) predict(fit, new1, type = "link", method = m))
  expect_equal(unname(diff(vals)), c(0, 0), tolerance = 1e-10)
})

# ---- interval ----------------------------------------------------------------

test_that("response interval: correct columns including se.link, lwr < fit < upr", {
  fit <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  out <- predict(fit, new1, interval = 0.95)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("fit", "lwr", "upr", "se.link"))
  expect_true(out$lwr < out$fit && out$fit < out$upr)
  expect_true(out$se.link > 0)
})

test_that("wider interval is wider than narrower interval", {
  fit  <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  ci95 <- predict(fit, new1, interval = 0.95)
  ci80 <- predict(fit, new1, interval = 0.80)
  expect_true(ci80$lwr > ci95$lwr)
  expect_true(ci80$upr < ci95$upr)
})

test_that("predict errors when interval is out of range", {
  fit <- bpmfit(y ~ x1, data = toy, prior = flat())
  expect_error(predict(fit, new1, interval = 1.2),  "interval")
  expect_error(predict(fit, new1, interval = -0.1), "interval")
})

test_that("interval = 0 returns fit and se.link only (type = 'response')", {
  fit <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  out <- predict(fit, new1, interval = 0)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("fit", "se.link"))
  expect_true(out$fit > 0 && out$fit < 1)
  expect_true(out$se.link > 0)
})

test_that("interval = 0 returns fit and se.link only (type = 'link')", {
  fit <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  out <- predict(fit, new1, type = "link", interval = 0)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("fit", "se.link"))
  expect_true(is.finite(out$fit))
  expect_true(out$se.link > 0)
})

test_that("interval = 0 se.link matches the se.link from a full interval", {
  fit  <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  se0  <- predict(fit, new1, interval = 0)$se.link
  se95 <- predict(fit, new1, interval = 0.95)$se.link
  expect_equal(se0, se95)
})

# ---- multi-row and edge cases ------------------------------------------------

test_that("predict on multiple rows returns correct length / nrow", {
  fit   <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  newdf <- data.frame(x1 = c(-1, 0, 1), x2 = c(0, 0, 0))
  expect_length(predict(fit, newdf), 3L)
  out_r <- predict(fit, newdf, interval = 0.95)
  expect_equal(nrow(out_r), 3L)
  expect_named(out_r, c("fit", "lwr", "upr", "se.link"))
  out_l <- predict(fit, newdf, type = "link", interval = 0.95)
  expect_equal(nrow(out_l), 3L)
  expect_named(out_l, c("fit", "lwr", "upr", "se.link"))
})

test_that("predict without newdata works when model = TRUE", {
  fit <- bpmfit(y ~ x1, data = toy, prior = flat(), model = TRUE)
  expect_length(predict(fit), nrow(toy))
})

test_that("predict without newdata errors when model = FALSE", {
  fit <- bpmfit(y ~ x1, data = toy, prior = flat(), model = FALSE)
  expect_error(predict(fit), "newdata")
})
