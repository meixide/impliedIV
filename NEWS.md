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
