# impliedIV 0.2.0 (2026-11-10)

* `lii_tmle(fixed = TRUE)` (default): estimator and standard errors for the
  parameter defined with the true instrument intensity, from the full
  efficient influence curve (new closed-form term for the estimation of the
  intensity). Columns `estimate_fix`, `se_fix`, ... in `estimates`.
* `lii_nuisance(pooled = FALSE)`: one model per visit instead of a pooled
  model (much faster with hapc).
* `lii_naive()`: comparator that ignores the instrument.
* `simulate_lii(confounder = "fresh")` and `truth_treatment()`.
* `inst/validation/eic_check.R`: exact numerical check of the efficient
  influence curve, the second-order expansion and the polynomial
  dose-response in a discrete model.

# impliedIV 0.1.0 (2026-11-09)

* `lii_project()`: KL projection of a desired treatment intervention (a
  factor, one factor per visit, or a rule) onto alpha-scalings of the
  instrument (common or per visit), with the LTMLE of its effect.
* First version: `lii_data()`, `lii_nuisance()`, `lii_tmle()` (LTMLE for
  alpha-scaling of a binary instrument, monotone censoring), `lii_implied()`,
  `lii_treatment_outcome()` + `calibrate_alpha()`, `lii_bootstrap()`,
  `lii_diagnostics()`, `preference_instrument()`, `binarize_instrument()`,
  `lag_sufficiency_test()`, `simulate_lii()`, learners `learner_hapc()`
  (default) and `learner_glm()`.
