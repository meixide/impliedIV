# impliedIV

Longitudinal **implied interventions**: what happens to an outcome if a
binary time-varying instrument (a prescription, a recommendation, a
monitoring visit) fires `alpha` times as often as it does now, when the
treatment it encourages (taking the drug, adherence) may be confounded with
the outcome by unmeasured factors.

- **Estimand.** `psi(alpha) = E[Y]` under the stochastic intervention
  `P(R(j) = 1 | history) -> alpha * P(R(j) = 1 | history)` at every visit. It
  never asks for the instrument where it does not occur, so there are no
  positivity problems. Contrasts are with `alpha = 1`, the current world.
- **Implied intervention on the treatment.**
  `g_j(alpha) = pi0 + alpha * lambda_j * (pi1 - pi0)`: only the
  instrument-driven part of treatment moves.
- **Estimator.** Longitudinal TMLE (sequential regression with targeting),
  monotone censoring, influence-curve and bootstrap inference.
- **Nuisances.** Estimated with Highly Adaptive Principal Components
  ([hapc](https://github.com/meixide/hapc)) by default, or with logistic
  regression for a quick pass.
- **Tools.**
  - physician-preference instruments, built leave-one-patient-out from
    strictly prior orders, with a trend term (`preference_instrument()`);
  - a test of "only the current instrument acts on the current treatment"
    (`lag_sufficiency_test()`);
  - calibration of alpha to "x% more treatment" (`calibrate_alpha()`).

Methodology: C. García Meixide and M. J. van der Laan, *Longitudinal implied
interventions* (in preparation).

## Install

```r
install.packages(c("remotes", "glmnet"))
remotes::install_github("meixide/hapc")        # nuisance learner (C++; needs a compiler)
remotes::install_github("meixide/impliedIV", build_vignettes = TRUE)
```

## Quick start (simulated data)

```r
library(impliedIV)
sim <- simulate_lii(n = 1000, K = 2, seed = 1)          # prescription -> adherence, hidden confounder
d   <- lii_data(sim$data, id = "id", time = "visit", instrument = "R", treatment = "A",
                outcome = "Y", covariates = "L", baseline = "W")
fit <- lii_tmle(d, alphas = c(0.5, 0.75, 1, 1.25, 1.5), learner = learner_glm)  # or learner_hapc
fit
plot(fit, diff = TRUE)
lii_diagnostics(fit)
```

**On real data, follow [AGENTS.md](AGENTS.md)**: data mapping, the order of
the analyses, diagnostics, what may and may not be shared, and an
interpretation template. Validation status: [VALIDATION.md](VALIDATION.md).
