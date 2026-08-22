source("TITE-AIDE.R")
if (!requireNamespace("rjags", quietly = TRUE) || !requireNamespace("coda", quietly = TRUE)) {
  message("model-matrix tests skipped: rjags and coda are required")
  quit(save = "no", status = 0L)
}
scenario <- aide_phase12_scenario(c(.05, .14, .28), c(.12, .30, .42),
                                  list(alpha_true = .10), list(alpha_true = .10))
for (toxicity_model in aide_phase12_model_choices()$toxicity) {
  for (efficacy_model in aide_phase12_model_choices()$efficacy) {
    config <- aide_phase12_config(Nmax = 9L, cohort_size = 3L, n_eval = 1L,
                                  toxicity = list(
                                    model = toxicity_model,
                                    skeleton = c(.05, .12, .25),
                                    time = list(arrival_rate = 0.018, t0 = 0, T_assess = 28),
                                    target = .30,
                                    cutoff = .95,
                                    beta_prior_mean = 0,
                                    beta_prior_sd = sqrt(2),
                                    carryover_prior = c(.15, .85)
                                  ),
      efficacy = list(model = efficacy_model, prior = c(0.5, 0.5), carryover_prior = c(0.15, 0.85), threshold = .20,
                      futility_cutoff = .90, min_eff_n_for_futility = 3L))
    trial <- simulate_AIDE_phase_I_II_event(config, scenario, seed = 100L)
    stopifnot(nrow(trial$admin) > 0L, is.null(trial$final$toxicity$draws), is.null(trial$final$efficacy$draws))
  }
}
# two_stage <- aide_phase12_config(allocation = "one_stage", Nmax = 12L, cohort_size = 3L, s1_Max = 3L, n_eval = 1L)
# trial <- simulate_AIDE_phase_I_II_event(two_stage, scenario, seed = 4L)
# stopifnot(nrow(trial$admin) > 0L, trial$stage1$n >= 0L, trial$stage2$n >= 0L)
# oc <- get_oc_sim_AIDE_phase_I_II(scenario, two_stage, ntrial = 2L, seed = 80L)
stopifnot(nrow(oc$trial_summary) == 2L)
message("model-matrix tests passed")
