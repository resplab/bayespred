# bayespred 0.1.0

* Initial CRAN submission.
* Fits logistic regression under four shrinkage priors: flat, Jeffreys,
  log-F(m), and Bayesian Ridge.
* Posterior mean prediction via 30-point Gauss-Hermite quadrature,
  MacKay approximation, and plug-in estimate.
* `bpmproj_pm()` for PM projection onto a simplified deployable model.
* `likelihood()` and `posterior()` for federated / multi-centre inference.
