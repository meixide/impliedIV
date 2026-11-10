# impliedIV — validation status

All checks use `simulate_lii()`: prescription → adherence, K = 2, and an
unmeasured confounder of adherence and outcome. The prescription law depends
on the past only through `F_j`. Censoring after each visit depends on `L(j)`;
in this DGP about 23% of subjects are censored before the outcome. Truth:
Monte Carlo with 10^6 draws of the intervened system, without censoring.
Scripts are in `inst/validation/`; results are in
`inst/validation/results/`.

## 1. Unit tests (`devtools::test()`)

52 expectations, all passing (≈ 10 s), including a check that learner_hapc calls cv.hapc with family = "binomial", norm = "1" (needed for the soft labels) and works with hard and soft labels, and the projection: a reachable target is reached, KL decreases, its effect equals lii_tmle at the same alpha, per-visit targets, and targets beyond the reachable range.
- Data preparation and censoring.
- Rejection of non-binary instruments and of NA covariates.
- Rescaling of continuous outcomes.
- At α = 1 without censoring, the estimate is within 0.5 SE of the sample mean.
- The influence curve has mean 0.
- Agreement with the truth at n = 3000.
- The implied intervention is affine in α and equals π0 at α = 0.
- The treatment-count outcome and calibration.
- The bootstrap.
- `preference_instrument()` against a brute-force computation.
- `lag_sufficiency_test()`.

`R CMD check --no-manual`: no errors. The only warning came from the
sandbox's lack of network access to CRAN (unstated test dependencies `glue`
and `waldo`, now declared in Suggests).

## 2. Influence-curve CIs, `learner_glm` (`coverage_glm.R 200 1000`)

200 replications, n = 1000.

| α | truth | bias | sd | mean SE | coverage | contrast with α = 1: bias | sd | mean SE | coverage |
|---|---|---|---|---|---|---|---|---|---|
| 0.50 | 0.3172 | −0.0004 | 0.0222 | 0.0211 | 0.925 | 0.0004 | 0.0134 | 0.0128 | 0.920 |
| 0.75 | 0.3028 | −0.0011 | 0.0181 | 0.0174 | 0.935 | −0.0002 | 0.0061 | 0.0058 | 0.920 |
| 1.00 | 0.2882 | −0.0009 | 0.0166 | 0.0160 | 0.950 | — | — | — | — |
| 1.25 | 0.2756 | −0.0020 | 0.0174 | 0.0165 | 0.935 | −0.0011 | 0.0059 | 0.0055 | 0.915 |
| 1.50 | 0.2631 | −0.0024 | 0.0202 | 0.0189 | 0.925 | −0.0015 | 0.0121 | 0.0112 | 0.930 |

- Bias is negligible.
- The influence-curve SE is about 95% of the Monte Carlo sd, so coverage is
  0.92–0.95. The Monte Carlo error of a coverage with 200 replications is
  about ±0.03.
- At first the shortfall was attributed to the variability of λ̂, which the
  influence-curve SE ignores. **Section 3 rules that out in this DGP**: the
  bootstrap, which refits λ̂, gives the same SEs. The cause is still open:
  finite-sample bias of the SE estimator, glm misspecification of the
  sequential regressions, or Monte Carlo error.

## 3. Bootstrap CIs, `learner_glm` (`coverage_bootstrap_glm.R 100 100 1000`)

100 replications, B = 100 resamples, n = 1000, percentile CIs.

| α | coverage | contrast coverage | mean bootstrap SE | mean contrast bootstrap SE |
|---|---|---|---|---|
| 0.5 | 0.91 | 0.92 | 0.0213 | 0.0128 |
| 1.0 | 0.94 | — | 0.0163 | — |
| 1.5 | 0.96 | 0.91 | 0.0193 | 0.0114 |

- The bootstrap SEs are practically equal to the influence-curve SEs of
  section 2 (0.0211, 0.0160 and 0.0189; contrasts 0.0128 and 0.0112).
- Coverage is in the same range, with Monte Carlo error ±0.045 for 100
  replications.
- In this DGP, the extra influence-curve term from estimating λ is
  negligible, and the bootstrap brings no improvement over the
  influence-curve CIs. It remains the principled route for the fixed
  parameter, and its cost is B full fits.

## 4. End-to-end run with `learner_hapc` (`run_hapc.R 1000`)

n = 1000, K = 2, α ∈ {0.5, 1, 1.5}; hapc 2.7.1 (`cv.hapc`, binomial,
`norm = "1"`, `max_degree = 2`, `n^0.7` PCs, 5 folds).

| α | estimate | 95% CI | truth |
|---|---|---|---|
| 0.5 | 0.344 | [0.300, 0.389] | 0.317 |
| 1.0 | 0.291 | [0.259, 0.323] | 0.289 |
| 1.5 | 0.253 | [0.215, 0.291] | 0.262 |

- Every lambda was interior to its grid; no fit failed.
- First stage: π1 − π0 had median 0.45.
- **Time: 403 s**, almost all of it in the five pooled nuisance fits (≈ 2 min
  each with 2735 rows). The outcome regressions take about 1 s each.
- hapc cost grows roughly like n³ in the rows of a fit, so the default cap
  `max_n = 4000` rows per fit means ≈ 5–6 min per pooled fit on a laptop.
- This is a single run: coverage with hapc has not been studied.

## 5. Theory checked exactly in a discrete model (`eic_check.R`)

Two visits, seven binary variables, 128 cells; every functional is computed
exactly.
- Numerical pathwise derivatives of psi(alpha) along random paths agree with
  the inner product of the efficient influence curve and the score to
  3e-11. Without the term for the estimation of the intensity they differ
  by up to 0.02.
- The two sides of the exact second-order expansion agree to 1e-17.
- Third differences of alpha -> psi(alpha) vanish (polynomial of degree 2).

## 6. Simulation study of the paper (`sim_paper.R`, `make_tables.R`)

Three visits, no censoring, per-visit nuisance fits, two scenarios for the
unmeasured confounder (persistent, renewed at each visit). Estimators:
`lii_tmle(fixed = TRUE)` and the comparator `lii_naive()`. Summary in
`inst/validation/results/paper/summary.csv`; the paper's Section 7 has the
tables. In short, for the contrast psi(alpha) - psi(1):
- logistic regression, n = 500, 1000, 2000 (200 to 345 replicates): bias
  below 0.001, coverage 0.92-0.97, with or without the intensity term;
- hapc, n = 500, persistent scenario (about 100 replicates): bias about a
  tenth of the contrast, coverage about 0.9; more replicates and larger n
  need a cluster;
- the comparator that ignores the instrument recovers one half to three
  quarters of the contrast, and its coverage goes to zero with n;
- `lag_sufficiency_test()` rejects at about the nominal rate in the renewed
  scenario and increasingly often in the persistent one.

## Not validated yet

- Coverage with hapc, and of the bootstrap with hapc.
- K > 2, rare instruments, strong censoring, continuous outcomes in
  simulation.
- `binarize_instrument()` (by design only a sensitivity analysis).
- Interpretation of `lag_sufficiency_test()` with overlapping windows: its
  power, and the size under heavy overlap.
