# Shared toy dataset used across test files.
set.seed(42)
n_toy <- 300
toy <- data.frame(
  x1  = rnorm(n_toy),
  x2  = rnorm(n_toy),
  grp = factor(sample(c("A", "B"), n_toy, replace = TRUE)),
  y   = rbinom(n_toy, 1L,
               prob = plogis(0.4 + 0.6 * rnorm(n_toy) - 0.4 * rnorm(n_toy)))
)
new1 <- data.frame(x1 = 0.5, x2 = -0.3, grp = factor("A", levels = c("A", "B")))
