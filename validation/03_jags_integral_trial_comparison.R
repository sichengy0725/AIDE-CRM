args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
validation_dir <- if (length(file_arg)) {
  dirname(normalizePath(file_arg[1L], winslash = "/"))
} else {
  source_files <- vapply(sys.frames(), function(frame) {
    value <- frame$ofile
    if (is.null(value)) "" else as.character(value)[1L]
  }, character(1))
  source_file <- tail(source_files[nzchar(source_files)], 1L)
  if (length(source_file)) dirname(normalizePath(source_file, winslash = "/")) else normalizePath(getwd(), winslash = "/")
}

project_dir <- dirname(validation_dir)
source(file.path(project_dir, "TITE-AIDE-Rebuild", "TITE-AIDE.R"))

scenario_table <- utils::read.csv(
  file.path(project_dir, "Set_5dose_adaptive_r_37_true_MTD_OBD_summary_lambda1.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
scenario18 <- scenario_table[scenario_table$Scenario == 18L, , drop = FALSE]
scenario <- aide_phase12_scenario(
  p_true = as.numeric(scenario18[paste0("Tox_Dose", 1:5)]),
  e_true = as.numeric(scenario18[paste0("Eff_Dose", 1:5)]),
  toxicity_ipde_dgm = list(alpha_true = 0),
  efficacy_ipde_dgm = list(alpha_true = 0)
)

make_comparison_config <- function(toxicity_model, efficacy_model) {
  aide_phase12_config(
    allocation = "one_stage",
    cohort_size = 3L,
    Nmax = 6L,
    s1_Max = 6L,
    N_s2 = 6L,
    n_eval = 3L,
    m_U = 6L,
    cycle_max = 1L,
    start_dose = 1L,
    time = list(arrival_rate = 1 / 14, t0 = 0, T_assess = 28,
                dlt_dist = "uniform", efficacy_dist = "uniform"),
    toxicity = list(
      model = toxicity_model,
      target = as.numeric(scenario18$Target_Toxicity),
      cutoff = 0.95,
      skeleton = c(0.15, 0.20, 0.30, 0.35, 0.45),
      beta_prior_mean = 0,
      beta_prior_sd = sqrt(2),
      carryover_prior = c(0.15, 0.85),
      n_chains = 2L,
      n_adapt = 500L,
      n_burnin = 500L,
      n_iter = 2000L,
      thin = 1L
    ),
    efficacy = list(
      model = efficacy_model,
      prior = c(0.5, 0.5),
      carryover_prior = c(0.15, 0.85),
      threshold = 0.20,
      futility_cutoff = 0.85,
      min_eff_n_for_futility = 3L,
      n_chains = 2L,
      n_adapt = 500L,
      n_burnin = 500L,
      n_iter = 2000L,
      thin = 1L
    ),
    integration = list(r_grid_n = 9L, beta_lower = -8, beta_upper = 8, rel_tol = 1e-5),
    utility = list(type = 1L, lambda_T = 1, scores = c(0, 1, -1, 0), tie_break = "lower_dose")
  )
}

jags_config <- make_comparison_config("discount_r", "shared_carryover")
integral_config <- make_comparison_config("integral_discount_r", "integral_shared_carryover")
trial_seed <- 20260818L

jags_time <- system.time({
  jags_trial <- simulate_AIDE_phase_I_II_event(jags_config, scenario, seed = trial_seed)
})

integral_time <- system.time({
  integral_trial <- simulate_AIDE_phase_I_II_event(integral_config, scenario, seed = trial_seed)
})

trial_time_comparison <- data.frame(
  scenario = 18L,
  version = c("JAGS", "integral"),
  toxicity_model = c("discount_r", "integral_discount_r"),
  efficacy_model = c("shared_carryover", "integral_shared_carryover"),
  user_seconds = c(unname(jags_time["user.self"]), unname(integral_time["user.self"])),
  system_seconds = c(unname(jags_time["sys.self"]), unname(integral_time["sys.self"])),
  elapsed_seconds = c(unname(jags_time["elapsed"]), unname(integral_time["elapsed"])),
  n_administrations = c(nrow(jags_trial$admin), nrow(integral_trial$admin)),
  trial_end_time = c(jags_trial$final$final_time, integral_trial$final$final_time),
  MTD = c(jags_trial$final$MTD, integral_trial$final$MTD),
  OBD = c(jags_trial$final$OBD, integral_trial$final$OBD),
  stringsAsFactors = FALSE
)

results_dir <- file.path(validation_dir, "results")
dir.create(results_dir, showWarnings = FALSE)
write.csv(scenario18, file.path(results_dir, "scenario18_truth.csv"), row.names = FALSE)
write.csv(trial_time_comparison, file.path(results_dir, "scenario18_jags_integral_trial_time.csv"), row.names = FALSE)
saveRDS(
  list(jags_trial = jags_trial, integral_trial = integral_trial, trial_time_comparison = trial_time_comparison),
  file.path(results_dir, "scenario18_jags_integral_trials.rds")
)

print(trial_time_comparison)
