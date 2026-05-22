
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

# Low birth weight study (Hosmer & Lemeshow 1989), shipped with R via MASS
data(birthwt, package = "MASS")
birthwt$race <- factor(birthwt$race, levels = 1:3,
                       labels = c("white", "black", "other"))

# Fit under log-F(2) prior (default)
fit <- bpmfit(low ~ age + lwt + race + smoke + ht + ui,
              data = birthwt, prior = log_f(m = 2))

# New patient
new_pt <- data.frame(age = 26, lwt = 130, race = "white",
                     smoke = 1, ht = 0, ui = 0)

# Posterior mean (recommended)
predict(fit, new_pt)
#> [1] 0.2609625

# With 95% credible interval
predict(fit, new_pt, interval = 0.95)
#>         fit       lwr       upr     se.fit
#> 1 0.2609625 0.1589038 0.3869804 0.05872089

# Compare all four priors
priors <- list(flat = flat(), jeffreys = jeffreys(),
               logf = log_f(m = 2), ridge = bridge())
sapply(priors, function(p) {
  predict(bpmfit(low ~ age + lwt + race + smoke + ht + ui,
                 data = birthwt, prior = p), new_pt)
})
#>      flat  jeffreys      logf     ridge 
#> 0.2493921 0.2581651 0.2609625 0.2691777
```

## PM projection

`bpmproject()` regresses the posterior-mean soft labels onto a design
matrix, producing a standalone deployable model with no covariance
matrix. All three aspects default to the main fit but can be overridden:
`formula` (predictor set), `family` (link function), and `data` (can be
the development sample, local site data, or any external dataset).

``` r
# Self-projection: same predictors as the main model
proj <- bpmproject(fit)
coef(proj)
#> (Intercept)         age         lwt   raceblack   raceother       smoke 
#>  0.44166751 -0.01964943 -0.01408265  1.04757057  0.73220901  0.87796708 
#>          ht          ui 
#>  1.46541520  0.77642062
predict(proj, new_pt)
#> [1] 0.2646391

# Custom projection: reduce to key predictors
proj_simple <- bpmproject(fit, formula = ~ lwt + smoke + ht, data = birthwt)
coef(proj_simple)
#> (Intercept)         lwt       smoke          ht 
#>  0.89934199 -0.01592838  0.61634750  1.45267675
predict(proj_simple, new_pt)
#> [1] 0.3647052

# Linear probability model — useful for nomograms
proj_linear <- bpmproject(fit, family = gaussian(link = "identity"))
coef(proj_linear)
#>  (Intercept)          age          lwt    raceblack    raceother        smoke 
#>  0.525533333 -0.002728915 -0.002519029  0.194795303  0.130626195  0.164639888 
#>           ht           ui 
#>  0.300494513  0.168657481
```

The object returned by `bpmproject()` is fully self-contained and
requires only base R at deployment time.

## Likelihood and posterior

For multi-centre or federated settings, `likelihood()` returns the
unpenalised MLE and observed Fisher information — the data’s
contribution independent of the prior. `posterior()` returns the `bpm`
object: the MAP estimate, Laplace-approximated covariance, and link
function. A `bpm` object is self-contained and can be used with
`predict()` directly.

``` r
lik  <- likelihood(fit)    # unpenalised MLE (computes on-the-fly)
post <- posterior(fit)     # bpm object: MAP + posterior covariance + link

# predict() works directly on the bpm object (newdata always required)
predict(post, new_pt)
#> [1] 0.2609625
```

## Reference

Sadatsafavi M, Riley RD (2026). Progression to the mean: A practical
Bayesian workflow for the development and deployment of clinical
prediction models *\[<https://arxiv.org/abs/2605.19163>\]*.
