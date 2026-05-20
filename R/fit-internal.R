# Internal per-prior fit functions.
# Each accepts the design matrix X (with intercept column), binary response y,
# and family object; returns list(coefficients, vcov, fit_method, converged).
# vcov is on the original (unstandardised) predictor scale.

.fit_flat <- function(X, y, family) {
  fit <- glm.fit(X, y, family = family)
  list(
    coefficients = fit$coefficients,
    vcov         = .logistic_vcov(X, fit$fitted.values),
    fit_method   = "glm",
    converged    = isTRUE(fit$converged)
  )
}

.fit_jeffreys <- function(X, y, family) {
  fit <- brglm2::brglm_fit(X, y, family = family,
                           control = brglm2::brglmControl(type = "AS_mean"))
  list(
    coefficients = fit$coefficients,
    vcov         = .logistic_vcov(X, fit$fitted.values),
    fit_method   = "brglmFit",
    converged    = isTRUE(fit$converged)
  )
}

.fit_logf <- function(X, y, m, family) {
  p             <- ncol(X)
  intercept_col <- which(colnames(X) == "(Intercept)")
  slope_cols    <- setdiff(seq_len(p), intercept_col)
  n_slope       <- length(slope_cols)
  n_pseudo      <- 2L * n_slope

  X_pseudo <- matrix(0, n_pseudo, p, dimnames = list(NULL, colnames(X)))
  y_pseudo <- numeric(n_pseudo)
  w_pseudo <- rep(m / 2, n_pseudo)

  for (k in seq_len(n_slope)) {
    j <- slope_cols[k]
    X_pseudo[2L * k - 1L, j] <- 1; y_pseudo[2L * k - 1L] <- 1
    X_pseudo[2L * k,      j] <- 1; y_pseudo[2L * k]      <- 0
  }

  X_aug <- rbind(X, X_pseudo)
  y_aug <- c(y, y_pseudo)
  w_aug <- c(rep(1, length(y)), w_pseudo)

  # Non-integer weights (m/2 when m is odd) are intentional — suppress the
  # binomial "non-integer #successes" warning that glm.fit would otherwise emit.
  fit <- suppressWarnings(
    glm.fit(X_aug, y_aug, weights = w_aug, family = family)
  )
  list(
    coefficients = fit$coefficients,
    vcov         = .logistic_vcov(X_aug, fit$fitted.values, w_aug),
    fit_method   = "data_augmentation",
    converged    = isTRUE(fit$converged)
  )
}

.fit_bridge <- function(X, y, family) {
  Xs_raw <- X[, -1L, drop = FALSE]

  if (ncol(Xs_raw) == 0L)
    return(.fit_flat(X, y, family))   # intercept-only: ridge collapses to flat

  Xs  <- scale(Xs_raw)
  ctr <- attr(Xs, "scaled:center")
  scl <- attr(Xs, "scaled:scale")

  if (any(scl == 0))
    stop("Bayesian Ridge requires all predictors to have non-zero variance.")

  mus <- c(0, ctr)
  sds <- c(1, scl)

  p        <- ncol(Xs)
  S        <- diag(p)
  pen_list <- list(Xs = list(S))

  dat <- data.frame(y = y, Xs = I(Xs))
  fit <- mgcv::gam(
    y ~ Xs,
    data    = dat,
    family  = family,
    method  = "REML",
    paraPen = pen_list
  )

  beta_z <- coef(fit)
  V_z    <- vcov(fit, unconditional = TRUE)

  # Transform from standardised to original scale via delta method.
  # beta[1]   = beta_z[1] - sum_j (center_j / scale_j) * beta_z[j+1]
  # beta[j+1] = beta_z[j+1] / scale_j
  D          <- diag(1 / sds)
  D[1L, -1L] <- -mus[-1L] / sds[-1L]

  beta <- as.numeric(D %*% beta_z)
  V    <- D %*% V_z %*% t(D)

  nms         <- colnames(X)
  names(beta) <- nms
  rownames(V) <- nms
  colnames(V) <- nms

  list(
    coefficients = beta,
    vcov         = V,
    fit_method   = "mgcv_ridge",
    converged    = isTRUE(fit$converged)
  )
}

# Compute (X^T W X)^{-1} where W_i = mu_i*(1-mu_i)*weight_i.
# Falls back to SVD pseudoinverse for near-singular matrices.
.logistic_vcov <- function(X, mu, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, length(mu))
  W   <- mu * (1 - mu) * weights
  XW  <- X * sqrt(W)
  M   <- crossprod(XW)
  nms <- colnames(X)
  V   <- tryCatch(
    solve(M),
    error = function(e) {
      warning("Near-singular information matrix; using SVD pseudoinverse.",
              call. = FALSE)
      sv  <- svd(M)
      tol <- max(sv$d) * .Machine$double.eps * max(dim(M))
      inv <- ifelse(sv$d > tol, 1 / sv$d, 0)
      sv$v %*% diag(inv, nrow = length(inv)) %*% t(sv$u)
    }
  )
  rownames(V) <- nms
  colnames(V) <- nms
  V
}
