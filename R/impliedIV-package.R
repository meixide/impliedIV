#' impliedIV: longitudinal implied interventions by scaling an instrument
#'
#' Workflow: [lii_data()] -> [lii_describe()] -> [lii_project()] (desired
#' treatment intervention -> closest instrument intervention and its effect)
#' and/or [lii_tmle()] (a grid of instrument scalings) ->
#' [lii_diagnostics()] -> optionally [lii_bootstrap()],
#' [lii_treatment_outcome()] + [calibrate_alpha()]. Tools for building
#' instruments: [preference_instrument()], [binarize_instrument()],
#' [lag_sufficiency_test()]. Simulated data: [simulate_lii()].
#'
#' Read `vignette("impliedIV")` and, for an automated agent, the file
#' `AGENTS.md` at the top of the source package.
#' @keywords internal
"_PACKAGE"
