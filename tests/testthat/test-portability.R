# Critical: bpmfit objects must survive saveRDS/readRDS and predict correctly
# using only base R (no mgcv or brglm2 needed after loading).

test_that("flat-prior bpmfit survives save/reload and predicts identically", {
  fit         <- bpmfit(y ~ x1 + x2, data = toy, prior = flat())
  pred_before <- predict(fit, new1)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(fit, tmp)
  pred_after <- predict(readRDS(tmp), new1)
  expect_equal(pred_before, pred_after)
})

test_that("bpmproject survives save/reload and predicts identically", {
  proj        <- bpmproject(bpmfit(y ~ x1 + x2, data = toy, prior = log_f(m = 2)))
  pred_before <- predict(proj, new1)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(proj, tmp)
  pred_after <- predict(readRDS(tmp), new1)
  expect_equal(pred_before, pred_after)
})

test_that("bridge-prior bpmfit survives save/reload and predicts identically", {
  fit        <- bpmfit(y ~ x1 + x2, data = toy, prior = bridge())
  pred_before <- predict(fit, new1)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(fit, tmp)
  pred_after <- predict(readRDS(tmp), new1)
  expect_equal(pred_before, pred_after)
})

test_that("bpmfit object stores no mgcv/brglm2 class objects", {
  fit  <- bpmfit(y ~ x1 + x2, data = toy, prior = bridge())
  # Check top-level and posterior-level for any gam/brglm class objects
  top_classes  <- unlist(lapply(fit, function(x) class(x)))
  post_classes <- unlist(lapply(fit$posterior, function(x) class(x)))
  expect_length(grep("gam|brglm", c(top_classes, post_classes), value = TRUE), 0L)
})

test_that("predict after reload works with factor predictor", {
  fit <- bpmfit(y ~ x1 + grp, data = toy, prior = flat())
  nd  <- data.frame(x1 = 0.2, grp = factor("B", levels = c("A", "B")))
  pred_before <- predict(fit, nd)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(fit, tmp)
  pred_after <- predict(readRDS(tmp), nd)
  expect_equal(pred_before, pred_after)
})

test_that("jeffreys-prior bpmfit survives save/reload", {
  fit        <- bpmfit(y ~ x1 + x2, data = toy, prior = jeffreys())
  pred_before <- predict(fit, new1)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(fit, tmp)
  pred_after <- predict(readRDS(tmp), new1)
  expect_equal(pred_before, pred_after)
})

test_that("bpm posterior object survives save/reload and predicts identically", {
  post        <- posterior(bpmfit(y ~ x1 + x2, data = toy, prior = log_f(m = 2)))
  pred_before <- predict(post, new1)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(post, tmp)
  pred_after <- predict(readRDS(tmp), new1)
  expect_equal(pred_before, pred_after)
})
