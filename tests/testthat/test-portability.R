# Critical: bpm objects must survive saveRDS/readRDS and predict correctly
# using only base R (no mgcv or brglm2 needed after loading).

test_that("flat-prior bpm survives save/reload and predicts identically", {
  fit         <- bpm(y ~ x1 + x2, data = toy, prior = flat())
  pred_before <- predict(fit, new1)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(fit, tmp)
  pred_after <- predict(readRDS(tmp), new1)
  expect_equal(pred_before, pred_after)
})

test_that("bpmproj_pm survives save/reload and predicts identically", {
  proj        <- project_pm(bpm(y ~ x1 + x2, data = toy, prior = log_f(m = 2)))
  pred_before <- predict(proj, new1)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(proj, tmp)
  pred_after <- predict(readRDS(tmp), new1)
  expect_equal(pred_before, pred_after)
})

test_that("bridge-prior bpm survives save/reload and predicts identically", {
  fit        <- bpm(y ~ x1 + x2, data = toy, prior = bridge())
  pred_before <- predict(fit, new1)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(fit, tmp)
  pred_after <- predict(readRDS(tmp), new1)
  expect_equal(pred_before, pred_after)
})

test_that("bpm object stores no mgcv/brglm2 class objects", {
  fit  <- bpm(y ~ x1 + x2, data = toy, prior = bridge())
  objs <- unlist(lapply(fit, function(x) class(x)))
  expect_length(grep("gam|brglm", objs, value = TRUE), 0L)
})

test_that("predict after reload works with factor predictor", {
  fit <- bpm(y ~ x1 + grp, data = toy, prior = flat())
  nd  <- data.frame(x1 = 0.2, grp = factor("B", levels = c("A", "B")))
  pred_before <- predict(fit, nd)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(fit, tmp)
  pred_after <- predict(readRDS(tmp), nd)
  expect_equal(pred_before, pred_after)
})

test_that("jeffreys-prior bpm survives save/reload", {
  fit        <- bpm(y ~ x1 + x2, data = toy, prior = jeffreys())
  pred_before <- predict(fit, new1)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(fit, tmp)
  pred_after <- predict(readRDS(tmp), new1)
  expect_equal(pred_before, pred_after)
})
