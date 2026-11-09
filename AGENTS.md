# impliedIV — manual for an AI agent (or a human) running it on real data

This file is self-contained. Read all of it before touching the data. The
statistical background is in `vignette("impliedIV")` and in the help pages
(`?lii_tmle`, `?lii_data`, `?preference_instrument`).

## 1. What the package estimates

There is a binary **instrument** `R(j)` at each visit `j = 0..K` (e.g. "the
drug was prescribed at this visit", "the patient was monitored"). It
encourages a binary **treatment** `A(j)` ("the patient took the drug",
adherence). The treatment may be confounded with the outcome by things that
were not measured. The outcome `Y` is measured after visit `K`.

The package estimates `psi(alpha)`: the mean outcome, had nobody been censored,
if the instrument fired `alpha` times as often as it does now, for every type
of patient. Formally, its probability given the covariate-treatment history,
`lambda_j(F_j)`, becomes `alpha * lambda_j(F_j)`.
- `alpha = 1`: the current world.
- `alpha = 1.2`: 20% more prescriptions, proportionally to current practice.
- `alpha = 0.5`: half as many.

The usual question is about the treatment ("what if patients took the drug
20% more often?"), which is not identified when the treatment is confounded.
`lii_project()` turns it into the closest identifiable question. It finds the
instrument scaling `alpha` (one per visit if wanted) whose implied treatment
mechanism is closest, in Kullback–Leibler divergence, to the desired one.
Then it estimates the effect of that projected intervention.
- The instrument can move the treatment only between `pi0` (it never fires)
  and `pi1` (it always fires).
- A desired intervention outside that range is not reachable, and the result
  is then the effect of the closest reachable intervention.
- Always report `target_p_treat` against `implied_p_treat`.

It also reports the **implied intervention** on the treatment,
`g_j(alpha) = pi0 + alpha * lambda_j * (pi1 - pi0)`, where `pi_r = P(A(j)=1 | R(j)=r, history)`.
`pi0` is untouched; only the instrument-driven part of treatment moves.

Assumptions (say them in any report):
1. **Sequential randomisation of the instrument.** Given the observed past,
   the instrument is as good as random. The covariates must contain what the
   prescriber reacts to.
2. **Positivity of the instrument.** `0 < P(R(j)=1 | past) < 1`.
3. **Censoring explained by the observed past.**
4. For reading `psi(alpha)` as the effect of the implied intervention on the
   *treatment*:
   - the instrument affects the outcome only through the treatment
     (exclusion);
   - only the current instrument acts on the current treatment (test with
     `lag_sufficiency_test()`).

No assumption of "no unmeasured confounding of the treatment" is needed.

## 2. Non-negotiable rules

- **Patient-level data never leave the analyst's machine.**
  - Do not paste rows, do not upload data, and do not send any file
    containing individual rows.
  - `fit$ic` (one value per subject) and saved fit objects (`saveRDS(fit)`)
    count as individual-level data.
  - Share only the aggregated tables listed in section 6.
- **Do not edit the package code.** If something fails or looks wrong, write
  down the exact call, the error and `sessionInfo()`, and report it to the
  package author (Carlos García Meixide, carlos.garcia@icmat.es).
- **The instrument and the treatment must be different variables.**
  "Prescribed" is an instrument; "took it / filled it" is a treatment.
  Setting `R(j)` = "the doctor prescribed" and `A(j)` = "the doctor prescribed"
  (or the same thing measured twice) makes the analysis meaningless. If the
  data only have prescriptions, see section 4.
- **No hazard ratios, no Cox models.** The estimands are means and risk
  differences.
- **Do not read meaning into estimates without support.** Check the
  diagnostics first (section 5).

## 3. Installation and check

```r
install.packages(c("remotes", "glmnet", "testthat"))
remotes::install_github("meixide/hapc")     # nuisance learner (C++; needs a compiler)
remotes::install_github("meixide/impliedIV", build_vignettes = TRUE)
# alternative: git clone https://github.com/meixide/impliedIV && devtools::install("impliedIV")
library(impliedIV)
```

Run the tests from a clone of the repository: `devtools::test()` (≈ 5 s, all must
pass). If hapc does not install, the whole pipeline still runs with
`learner = learner_glm`, but report that hapc was not used.

## 4. Mapping the data

`lii_data()` needs a **long** table, one row per (subject, visit):

| Argument | Meaning | Rules |
|---|---|---|
| `id` | subject | no NA |
| `time` | orders the visits within subject | visits are re-indexed 0, 1, 2, … in this order |
| `covariates` | `L(j)`, measured **at or before** the visit, before the instrument | numeric or factor, no NA (impute or add missingness indicators first) |
| `baseline` | constant within subject (age at entry, sex, …) | optional |
| `instrument` | `R(j)`, 0/1 | no NA at observed visits |
| `treatment` | `A(j)`, 0/1, after the instrument | no NA at observed visits |
| `outcome` | `Y`, after the last visit; taken from the subject's last row | binary or continuous; NA = censored |

**Choosing visits and `K`.**
- A "visit" is a decision point at which the instrument can fire: an
  encounter, or a fixed period such as a month.
- `K` is the last visit index used. Pick it so that `lii_describe(d)` shows
  enough subjects still observed at visit `K` (a few hundred at least).
- Subjects with fewer visits are treated as censored; censoring is handled by
  inverse weighting.
- Keep `K` small at first (2–5).

**Recommended mapping for antibiotic prescribing (prescription → adherence).**
- `R(j)` = antibiotic prescribed at encounter `j`.
- `A(j)` = prescription filled or taken (pharmacy dispensing, refill or
  self-report).
- `L(j)` = what the clinician saw: symptoms, vital signs, labs, diagnosis
  codes, prior antibiotic use.
- `Y` = clinical outcome after the last encounter (resolution, complication,
  return visit, resistant infection, …).

**If the data contain only prescriptions (no adherence)** and the candidate
instrument is the prescriber's preference:
1. Build it with `preference_instrument(orders, provider, patient, time, treated, window, min_orders, trend = TRUE)`.
   - It uses only orders strictly **before** the encounter.
   - It **excludes the patient's own orders**, which is necessary: otherwise
     the instrument contains the patient's own decision.
   - It also returns the trend from the previous window.
2. Test the "current instrument only" assumption with
   `lag_sufficiency_test(data, id, time, instrument = "pref", treatment, covariates, lags = 2)`.
   It works with a continuous instrument. A small p-value means older
   windows still predict today's prescription beyond the current window.
   Try a shorter window, or the (current, trend) pair. AUROC/AUPRC stability
   across windows does **not** test this.
3. `lii_tmle()` needs a binary instrument. `binarize_instrument(pref)` gives
   "seen by a high-preference prescriber". This is only a sensitivity
   analysis: compare several cut-offs and flag it as such. Then `A(j)` =
   prescribed.

## 5. Workflow

```r
library(impliedIV)
d <- lii_data(df, id = "patient_id", time = "enc_date", instrument = "rx",
              treatment = "filled", outcome = "y", covariates = c(...), baseline = c(...), K = 3)
d
lii_describe(d)                       # n per visit, P(R), P(A), P(A | R = 1) - P(A | R = 0)

## 1) pipeline check, fast
lii_reset_log()
fit_glm <- lii_tmle(d, alphas = c(0.5, 0.75, 1, 1.25, 1.5), learner = learner_glm)

## 2) main analysis, nuisances with hapc
lii_reset_log()
nu  <- lii_nuisance(d, learner = learner_hapc)          # fits that do not depend on alpha

## 2a) the question about the TREATMENT: project it and estimate its effect
pr <- lii_project(d, target = 1.2, nuisance = nu, learner = learner_hapc)   # "20% more treatment"
pr   # alpha, target vs implied treatment probability per visit, reachable range, effect + CI
## one factor per visit: target = c(1.1, 1.2, 1.3); any rule: target = function(visit, history) ...

## 2b) the whole curve in alpha (instrument scalings)
fit <- lii_tmle(d, alphas = c(0.5, 0.75, 1, 1.25, 1.5), learner = learner_hapc, nuisance = nu)
fit; plot(fit); plot(fit, diff = TRUE)
dg <- lii_diagnostics(fit)

## 3) interpretable x-axis: how much more treatment does each alpha mean?
ft <- lii_tmle(lii_treatment_outcome(d), alphas = fit$settings$alphas, learner = learner_hapc)
calibrate_alpha(ft, ratio = c(1.1, 1.2))

## 4) optional: bootstrap for the parameter defined with the true lambda (slow: B full fits)
bt <- lii_bootstrap(d, alphas = c(0.5, 1, 1.5), learner = learner_hapc, B = 200)
```

**Runtime.**
- With `learner_glm`, seconds.
- With `learner_hapc`, a fit's cost grows roughly like (rows)^3. In testing,
  2700 rows took about 2 min and 1000 rows a few seconds; n = 1000, K = 2
  and 3 alphas took 7 min in total.
- Pooled models use `n x (K+1)` rows, subsampled to `max_n = 4000` rows per
  fit (`learner_args = list(max_n = ...)`), so expect ≈ 5–6 min per pooled
  fit on a laptop.
- There are 5 pooled fits, plus 2 outcome fits per visit and per alpha.
  Start with few alphas.
- The bootstrap multiplies everything by `B`; run it last, on a server if
  possible.

**Read the diagnostics before the estimates.**

| Check | Where | Problem if |
|---|---|---|
| Subjects per visit | `lii_describe(d)` | < 200 at visit `K` → lower `K` |
| First-stage strength `pi1 - pi0` | `dg$first_stage` | median < 0.05 → weak instrument; a flat curve in alpha says nothing about the treatment |
| Capping | `fit$estimates$p_capped`, `dg$alpha_max_uncapped` | > 5% → alphas above `1 / max(lambda)` change the intervention; report and prefer smaller alphas |
| Weights | `fit$estimates$min_n_eff`, `max_weight` | `min_n_eff` < 10% of n, or huge weights → unstable; narrow the alpha range |
| Censoring | `dg$censoring` | probabilities near 0 |
| Learner fits | `dg$fits` (`lii_fit_log()`) | many `upper_edge`/`null_model` (learner ignored the regressors), `hapc_failed_glm_used`, `subsampled` |
| Assumption: current instrument only | `lag_sufficiency_test()` | p < 0.05 |
| hapc vs glm | `fit` vs `fit_glm` | large disagreement → investigate before reporting |
| Projection reaches the target | `pr$projection`: `implied_p_treat` vs `target_p_treat`, `reachable_min/max` | target outside the reachable range → the effect is that of the closest reachable intervention; say so |

**Inference.**
- `fit$estimates` SEs come from the influence curve. They are valid for the
  *data-adaptive* parameter, which treats the estimated instrument intensity
  as known. They can be slightly short for the parameter defined with the
  true intensity: in simulation, coverage was 0.92–0.95 with glm.
- `lii_bootstrap()` targets the latter.
- Report which one you used.

## 6. What to report back (aggregated only)

1. Data summary.
   - `print(d)` and `lii_describe(d)`.
   - Definitions of `R`, `A`, `Y`, the visits, `K`, the covariates.
   - Cohort selection and number of subjects excluded, with reasons.
2. `pr$alpha`, `pr$projection`, `pr$effect`; `fit$estimates` and `fit_glm$estimates`.
3. `subset(fit$implied, visit == 0)`, or the whole table.
4. `dg$first_stage`, `dg$alpha_max_uncapped`, `dg$weights`, `dg$censoring`,
   `dg$fits`.
5. `lag_sufficiency_test(...)` output (coefficients and p-value), if
   applicable.
6. `calibrate_alpha(...)` and the bootstrap table, if run.
7. `sessionInfo()` and the package version (`packageVersion("impliedIV")`).
8. Anything unexpected, verbatim.

Never: rows of data, `fit$ic`, `nu` or `fit$nuisance` (they hold
predictions for each subject), saved `.rds` fits.

## 7. Interpretation template

> "If [prescriptions] were [20% more frequent] at every encounter,
> proportionally to current practice, the [risk of Y] would be [estimate]
> (95% CI [..]), compared with [estimate at alpha = 1] now (difference [..]).
> This corresponds to [x%] more [treatment] in that world
> (`calibrate_alpha`). It relies on [assumptions 1–3]. Reading it as the
> effect of the implied change in [treatment] also relies on assumption 4
> (lag test p = [..])."

Keep "estimate" and "interpretation" separate. A causal claim inherits every
assumption; say which ones are plausible in these data, and why.
