# TITE-AIDE rebuild

This folder is an independent, modular reconstruction of the AIDE/IPDE Phase I-II TITE simulator specified in `Presentation 8-10-2026/AIDE_IPDE_Phase_I_II_TITE_Rebuild_Specification.pdf`.

The implementation has one event clock for both allocation strategies. A closed cohort is re-opened only after an assignable new-patient or structurally eligible retreat event; every opened cohort has one frozen dose. The current `n_eval` rule gates both stay and escalation decisions; de-escalation is immediately permitted. A blocked new trigger is dropped before entering the queue; a blocked retreat trigger stays queued.

## Layout

- `R/` - configuration, state, endpoint, model, decision, event, final-analysis, and OC modules.
- `inst/jags/` - the four delayed-outcome efficacy likelihoods specified for E-RD, E-RS, E-AD, and E-AS.  The R engine uses summary-only native adapters by default; these files support a JAGS backend without changing the public contract.
- `tests/testthat/` - base-R executable regression checks for the finalised event-loop rules.
- `TITE-AIDE.R` - source this one file to load the public API.

## Minimal use

```r
source("TITE-AIDE.R")
config <- aide_phase12_config(allocation = "one_stage", Nmax = 30, cohort_size = 3)
scenario <- aide_phase12_scenario(
  p_true = c(.05, .12, .24, .38, .55),
  e_true = c(.10, .22, .38, .45, .42),
  toxicity_ipde_dgm = list(alpha_true = .10),
  efficacy_ipde_dgm = list(alpha_true = .15)
)
trial <- simulate_AIDE_phase_I_II_event(config, scenario, seed = 1)
```

Toxicity uses the bundled JAGS skeleton TITE-CRM models: a random discount-`r` model or an additive previous-dose-`alpha` model. It uses the TITE likelihood `1 - w * theta` for pending non-DLTs. Efficacy uses one of the four bundled delayed-outcome JAGS models. `toxicity$beta_prior_mean` and `toxicity$beta_prior_sd` describe the CRM curve prior; `toxicity$carryover_prior` is the Beta prior for `r` or additive `alpha`. `toxicity$prior` is intentionally unsupported because it is not a CRM prior. Both `rjags` and `coda` must be available in the active R library.

The result retains only posterior summaries and audit tables; neither posterior draws nor model diagnostics are stored.

## Operating-characteristic runs

The runner and extractor mirror the project-level Phase I/II OC setup, including its 37-scenario source file and `IDX` seed convention:

```r
Rscript run_oc_TITE_AIDE_phase_I_II.R 1 1 1
Rscript extract_result_TITE_AIDE_phase_I_II.R
```

They write JAGS TITE-AIDE results under `oc_results_cluster_TITE_AIDE_phase_I_II_*` and extraction files under `OC_summary_TITE_AIDE_phase_I_II/`.

## Optional individual recycle screens

Set these in `aide_tite_oc_settings()` to apply a screen to a structurally eligible retreat candidate at the frozen next dose of an open cohort:

```r
apply_individual_toxicity_risk = TRUE,
toxicity_ipde_overdose_cutoff = .95,
apply_individual_efficacy_benefit = TRUE,
efficacy_ipde_min_increment = .05
```

The toxicity screen blocks a retreat when the posterior probability that its IPDE toxicity exceeds the target is above the cutoff. The efficacy screen blocks it when `P(IPDE efficacy) - P(regular efficacy)` is less than or equal to the requested minimum increment. A blocked retreat stays in the queue for reconsideration under a later cohort dose.

## Deliberate compatibility boundary

`simulate_AIDE_phase_I_II()` is retained as a lightweight compatibility wrapper for `config`/`scenario` calls and basic legacy scalar inputs.  Legacy recycle-destination controls (`ipde_design`, `flexible_ipde`, `random_crm_*`) are intentionally rejected because the rebuild has one design-selected cohort dose.
