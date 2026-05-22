# Internal per-prior fit functions.
# Each accepts the design matrix X (with intercept column), binary response y,
# and family object; returns list(coefficients, vcov, fit_method, converged).
# vcov is on the original (unstandardised) predictor scale.

# Extract (X'WX)^{-1} from a glm.fit / brglm_fit result, exactly as
# summary.glm() does it: chol2inv of the upper-triangular R factor stored in
# the QR decomposition at the final IWLS iteration.
.vcov_from_fit <- function(fit, nms) {
  p1 <- seq_len(fit$rank)
  V  <- chol2inv(fit$qr$qr[p1, p1, drop = FALSE])
  dimnames(V) <- list(nms, nms)
  V
}

.fit_flat <- function(X, y, family) {
  fit <- glm.fit(X, y, family = family)
  list(
    coefficients = fit$coefficients,
    vcov         = .vcov_from_fit(fit, colnames(X)),
    fit_method   = "glm",
    converged    = isTRUE(fit$converged)
  )
}

.fit_jeffreys <- function(X, y, family) {
  fit <- brglm2::brglm_fit(X, y, family = family,
                           control = brglm2::brglmControl(type = "AS_mean"))
  list(
    coefficients = fit$coefficients,
    vcov         = .vcov_from_fit(fit, colnames(X)),
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
    vcov         = .vcov_from_fit(fit, colnames(X)),
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
  dimnames(V) <- list(nms, nms)

  list(
    coefficients = beta,
    vcov         = V,
    fit_method   = "mgcv_ridge",
    converged    = isTRUE(fit$converged)
  )
}
