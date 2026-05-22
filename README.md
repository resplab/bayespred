
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bayespred

<!-- badges: start -->

<!-- badges: end -->

`bayespred` implements the pragmatic Bayesian workflow for clinical
prediction modelling described in Sadatsafavi & Riley (2026). It fits
logistic regression under four shrinkage priors and produces
**self-contained** model objects whose `predict()` method requires only
base R at deployment time.

| Prior              | Constructor    | Backend                       |
|--------------------|----------------|-------------------------------|
| Flat (unpenalised) | `flat()`       | `glm.fit`                     |
| Jeffreys           | `jeffreys()`   | `brglm2::brglm_fit`           |
| log-F(m, m)        | `log_f(m = 2)` | data augmentation + `glm.fit` |
| Bayesian Ridge     | `bridge()`     | `mgcv::gam` + REML            |

Three prediction methods on a `bpmfit` object:

| Method | Argument | Description |
|----|----|----|
| Plug-in | `"pe"` | `plogis(X β̂)` |
| Posterior mean (quadrature) | `"pm"` *(default)* | 30-point Gauss-Hermite |
| MacKay approximation | `"pm_mackay"` | Closed-form PM approximation |

## Installation

``` r
# Install dependencies
install.packages(c("mgcv", "brglm2"))

# Install from GitHub
# install.packages("pak")
pak::pak("resplab/bayespred")
```

## Quick start

``` r
library(bayespred)

set.seed(1)
d <- data.frame(x1 = rnorm(500), x2 = rnorm(500),
                y  = rbinom(500, 1, 0.25))

# Fit under log-F(2) prior (default)
fit <- bpmfit(y ~ x1 + x2, data = d, prior = log_f(m = 2))

# New patient — plain values, no need to specify factor levels
new_pt <- data.frame(x1 = 0.5, x2 = -1.2)

# Posterior mean (recommended)
predict(fit, new_pt)
#> [1] 0.2214141

# With 95% credible interval
predict(fit, new_pt, interval = 0.95)
#>         fit       lwr       upr     se.fit
#> 1 0.2214141 0.1695492 0.2806577 0.02836097

# Compare all four priors
priors <- list(flat = flat(), jeffreys = jeffreys(),
               logf = log_f(m = 2), ridge = bridge())
sapply(priors, function(p) {
  predict(bpmfit(y ~ x1 + x2, data = d, prior = p), new_pt)
})
#>      flat  jeffreys      logf     ridge 
#> 0.2212171 0.2229100 0.2214141 0.2415931
```

## PM projection

`bpmproj_pm()` regresses the posterior-mean soft labels onto a design
matrix, producing a standalone `bpmproj_pm` object — a simplified
deployable model with no covariance matrix. All three aspects default to
the main fit but can be overridden: `formula` (predictor set), `family`
(link function), and `data` (can be the development sample, local site
data, or any external dataset).

``` r
# Self-projection: same predictors as the main model
proj <- bpmproj_pm(fit)
coef(proj)
#> (Intercept)          x1          x2 
#>  -1.0187703  -0.1398351   0.1387906
predict(proj, new_pt)
#> [1] 0.221795

# Custom projection: reduce to a single predictor
proj_simple <- bpmproj_pm(fit, formula = ~ x1, data = d)
coef(proj_simple)
#> (Intercept)          x1 
#>  -1.0199888  -0.1451582
predict(proj_simple, new_pt)
#> [1] 0.251135
```

The `bpmproj_pm` object is fully self-contained and requires only base R
at deployment time.

## Portability

``` r
tmp <- tempfile(fileext = ".rds")
saveRDS(fit, tmp)

# In any R session — base R only needed at predict time:
fit2 <- readRDS(tmp)
predict(fit2, new_pt)
#> [1] 0.2214141
unlink(tmp)
```

## Likelihood and posterior

For multi-centre or federated settings, `likelihood()` returns the
unpenalised MLE and observed Fisher information — the data’s
contribution independent of the prior. `posterior()` returns the MAP
estimate and Laplace-approximated covariance used for prediction.

``` r
lik  <- likelihood(fit)   # unpenalised MLE (computes on-the-fly)
post <- posterior(fit)    # MAP + posterior covariance
```

## Reference

Sadatsafavi M, Riley RD (2026). Progression to the mean: A practical
Bayesian workflow for the development and deployment of clinical
prediction models *\[<https://arxiv.org/abs/2605.19163>\]*.
