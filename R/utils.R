# Gauss-Hermite quadrature nodes and weights (physicist convention).
# Uses eigendecomposition of the symmetric tridiagonal Jacobi matrix.
# Returns list(nodes, weights) satisfying:
#   integral f(x) exp(-x^2) dx ~= sum_k weights[k] * f(nodes[k])
# with sum(weights) = sqrt(pi).
.gh_nodes_weights <- function(K) {
  K <- as.integer(K)
  if (K == 1L) return(list(nodes = 0, weights = sqrt(pi)))
  e   <- sqrt(seq_len(K - 1L) / 2)
  J   <- matrix(0, K, K)
  idx <- seq_len(K - 1L)
  J[cbind(idx, idx + 1L)] <- e
  J[cbind(idx + 1L, idx)] <- e
  eig <- eigen(J, symmetric = TRUE)
  list(
    nodes   = eig$values,
    weights = sqrt(pi) * eig$vectors[1L, ]^2
  )
}

# Posterior mean of plogis(eta) when eta ~ N(mu, sigma^2).
# Uses K-point Gauss-Hermite quadrature. Vectorised over mu and sigma.
.pm_quadrature <- function(mu, sigma, K = 30L) {
  gh  <- .gh_nodes_weights(K)
  xk  <- sqrt(2) * gh$nodes       # scaled nodes
  wk  <- gh$weights / sqrt(pi)    # normalised weights (sum to 1)
  # t_mat[i, k] = mu[i] + sigma[i] * xk[k]
  t_mat <- sweep(outer(sigma, xk, `*`), 1L, mu, `+`)
  as.numeric(plogis(t_mat) %*% wk)
}

# MacKay (1992) closed-form approximation to the logistic-normal posterior mean.
.pm_mackay <- function(mu, sigma) {
  as.numeric(plogis(mu / sqrt(1 + pi * sigma^2 / 8)))
}

# Apply cached self-projection coefficients.
.pm_projection <- function(X_new, beta_proj) {
  as.numeric(plogis(X_new %*% beta_proj))
}

# Extract factor levels from a model frame (mirrors stats:::.getXlevels).
.get_xlevels <- function(terms_obj, mf) {
  vars <- vapply(attr(terms_obj, "variables"), deparse, "")[-1L]
  xlev <- lapply(vars, function(v) {
    col <- mf[[v]]
    if (!is.null(levels(col))) levels(col) else NULL
  })
  names(xlev) <- vars
  xlev[!vapply(xlev, is.null, NA)]
}
