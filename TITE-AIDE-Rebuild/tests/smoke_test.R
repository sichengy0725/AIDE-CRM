source("TITE-AIDE.R")
if (!requireNamespace("rjags", quietly = TRUE) || !requireNamespace("coda", quietly = TRUE))
  stop("smoke_test.R requires rjags and coda in the active R library.")
config <- aide_phase12_config(Nmax = 12L, cohort_size = 3L, n_eval = 1L)
scenario <- aide_phase12_scenario(c(.05, .15, .30, .45), c(.10, .28, .45, .40),
                                  toxicity_ipde_dgm = list(alpha_true = .10),
                                  efficacy_ipde_dgm = list(alpha_true = .15))
trial <- simulate_AIDE_phase_I_II_event(config, scenario, seed = 11L)
stopifnot(nrow(trial$admin) > 0L, all(trial$admin$dose == trial$admin$decision_next_dose))
print(list(administrations = nrow(trial$admin), decisions = nrow(trial$decision_log),
           MTD = trial$final$MTD, OBD = trial$final$OBD, stop_reason = trial$stop_reason))
