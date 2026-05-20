
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BayesCPM

<!-- badges: start -->

<!-- badges: end -->

`BayesCPM` implements the pragmatic Bayesian workflow for clinical
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

Four prediction methods:

| Method | Argument | Description |
|----|----|----|
| Plug-in | `"pe"` | `plogis(X β̂)` |
| Posterior mean (quadrature) | `"pm"` *(default)* | 30-point Gauss-Hermite |
| MacKay approximation | `"pm_mackay"` | Closed-form PM approximation |
| Self-projection | `"pm_proj"` | Simplified linear predictor (opt-in) |

## Installation

``` r
# Install dependencies
install.packages(c("mgcv", "brglm2"))

# Install from GitHub
# install.packages("pak")
pak::pak("resplab/BayesCPM")
```

## Quick start

``` r
library(BayesCPM)

set.seed(1)
d <- data.frame(x1 = rnorm(500), x2 = rnorm(500),
                y  = rbinom(500, 1, 0.25))

# Fit under log-F(2) prior (default)
fit <- bpm(y ~ x1 + x2, data = d, prior = log_f(m = 2))

# New patient
new_pt <- data.frame(x1 = 0.5, x2 = -1.2)

# Posterior mean (recommended)
predict(fit, new_pt)
#> [1] 0.2214141

# With 95% credible interval
predict(fit, new_pt, interval = 0.95)
#>         fit       lwr       upr   se.link
#> 1 0.2214141 0.1695492 0.2806577 0.1652133

# Compare all four priors
priors <- list(flat = flat(), jeffreys = jeffreys(),
               logf = log_f(m = 2), ridge = bridge())
sapply(priors, function(p) {
  predict(bpm(y ~ x1 + x2, data = d, prior = p), new_pt)
})
#>      flat  jeffreys      logf     ridge 
#> 0.2212171 0.2229100 0.2214141 0.2415931
```

## Self-projection

``` r
fit_sp <- bpm(y ~ x1 + x2, data = d, prior = log_f(m = 2), projpred = TRUE)
coef(fit_sp, type = "projection")
#> (Intercept)          x1          x2 
#>  -1.0187703  -0.1398351   0.1387906
predict(fit_sp, new_pt, method = "pm_proj")
#> [1] 0.221795
```

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

## Reference

Sadatsafavi M, Riley RD (2026). A practical Bayesian workflow for
clinical prediction modelling.
