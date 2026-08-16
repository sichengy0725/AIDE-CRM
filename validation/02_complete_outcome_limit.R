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

source(file.path(dirname(validation_dir), "TITE-AIDE-Rebuild", "TITE-AIDE.R"))
source(file.path(validation_dir, "validation_models.R"))

T_assess <- 28
complete_admin <- data.frame(
  admin_id = 1:12,
  patient_id = 1:12,
  cycle = 1L,
  t_start = seq(0, 22, by = 2),
  dose = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L),
  previous_dose = c(1L, 1L, 1L, 1L, 2L, 1L, 2L, 2L, 3L, 2L, 3L, 2L),
  assignment_type = c("new", "new", "retreat", "new", "new", "retreat", "new", "new", "new", "retreat", "new", "retreat"),
  dlt_final = c(0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 0L),
  eff_final = c(0L, 1L, 1L, 0L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 1L),
  stringsAsFactors = FALSE
)
complete_admin$assessment_end <- complete_admin$t_start + T_assess
complete_admin$t_dlt <- ifelse(complete_admin$dlt_final == 1L, complete_admin$t_start + 9, Inf)
complete_admin$t_response <- ifelse(complete_admin$eff_final == 1L, complete_admin$t_start + 12, Inf)

t_final <- max(complete_admin$assessment_end)
tite_toxicity_data <- aide_phase12_tite_toxicity_data(complete_admin, t_final, T_assess)
tite_efficacy_data <- aide_phase12_tite_efficacy_data(complete_admin, t_final, T_assess)
complete_toxicity_data <- transform(tite_toxicity_data, y = as.integer(y), weight = NULL, ascertained = NULL)
complete_efficacy_data <- transform(tite_efficacy_data, y = as.integer(y), weight = NULL, ascertained = NULL)

toxicity_models <- aide_phase12_model_choices()$toxicity
efficacy_models <- aide_phase12_model_choices()$efficacy
ndose <- 3L
current_dose <- 3L

make_validation_config <- function(toxicity_model, efficacy_model) {
  aide_phase12_config(
    allocation = "one_stage",
    Nmax = nrow(complete_admin),
    cohort_size = 3L,
    n_eval = 1L,
    m_U = 6L,
    start_dose = current_dose,
    time = list(arrival_rate = 0.2, t0 = 0, T_assess = T_assess),
    toxicity = list(
      model = toxicity_model,
      target = 0.30,
      cutoff = 0.95,
      skeleton = c(0.05, 0.15, 0.30),
      beta_prior_mean = 0,
      beta_prior_sd = sqrt(2),
      carryover_prior = c(1, 9),
      n_chains = 2L,
      n_adapt = 1000L,
      n_burnin = 1000L,
      n_iter = 8000L,
      thin = 1L
    ),
    efficacy = list(
      model = efficacy_model,
      prior = c(1, 1),
      carryover_prior = c(1, 9),
      threshold = 0.20,
      futility_cutoff = 0.90,
      min_eff_n_for_futility = 3L,
      n_chains = 2L,
      n_adapt = 1000L,
      n_burnin = 1000L,
      n_iter = 8000L,
      thin = 1L
    ),
    utility = list(type = 1L, lambda_T = 1, scores = c(0, 1, -1, 0), tie_break = "lower_dose")
  )
}

toxicity_fits <- lapply(seq_along(toxicity_models), function(i) {
  toxicity_model <- toxicity_models[i]
  config <- make_validation_config(toxicity_model, efficacy_models[1L])
  set.seed(20260820L + i)
  complete <- validation_fit_complete_toxicity(complete_toxicity_data, config, ndose)
  set.seed(20260920L + i)
  tite <- aide_fit_toxicity(tite_toxicity_data, config, ndose)
  list(model = toxicity_model, complete = complete, tite = tite)
})
names(toxicity_fits) <- toxicity_models

efficacy_fits <- lapply(seq_along(efficacy_models), function(i) {
  efficacy_model <- efficacy_models[i]
  config <- make_validation_config(toxicity_models[1L], efficacy_model)
  set.seed(20261020L + i)
  complete <- validation_fit_complete_efficacy(complete_efficacy_data, config, ndose)
  set.seed(20261120L + i)
  tite <- aide_fit_efficacy(tite_efficacy_data, config, ndose)
  list(model = efficacy_model, complete = complete, tite = tite)
})
names(efficacy_fits) <- efficacy_models

posterior_comparisons <- list()
decision_comparisons <- list()
comparison_id <- 0L

for (toxicity_model in toxicity_models) {
  for (efficacy_model in efficacy_models) {
    comparison_id <- comparison_id + 1L
    config <- make_validation_config(toxicity_model, efficacy_model)
    complete_toxicity_fit <- toxicity_fits[[toxicity_model]]$complete
    tite_toxicity_fit <- toxicity_fits[[toxicity_model]]$tite
    complete_efficacy_fit <- efficacy_fits[[efficacy_model]]$complete
    tite_efficacy_fit <- efficacy_fits[[efficacy_model]]$tite

    complete_state <- validation_decision_state(complete_admin, ndose, current_dose)
    tite_state <- validation_decision_state(complete_admin, ndose, current_dose)
    scenario <- list(ndose = ndose)
    complete_decision <- aide_make_design_decision(complete_state, complete_toxicity_fit, complete_efficacy_fit, config, scenario)
    tite_decision <- aide_make_design_decision(tite_state, tite_toxicity_fit, tite_efficacy_fit, config, scenario)
    complete_final <- validation_final_selection(complete_admin, complete_toxicity_fit, complete_efficacy_fit, config)
    tite_final <- validation_final_selection(complete_admin, tite_toxicity_fit, tite_efficacy_fit, config)

    posterior_comparisons[[comparison_id]] <- rbind(
      data.frame(
        toxicity_model = toxicity_model,
        efficacy_model = efficacy_model,
        endpoint = "toxicity",
        dose = seq_len(ndose),
        complete_p_regular_mean = complete_toxicity_fit$p_regular_mean,
        tite_p_regular_mean = tite_toxicity_fit$p_regular_mean,
        absolute_difference = abs(complete_toxicity_fit$p_regular_mean - tite_toxicity_fit$p_regular_mean),
        complete_p_ipde_mean = complete_toxicity_fit$p_ipde_mean,
        tite_p_ipde_mean = tite_toxicity_fit$p_ipde_mean,
        complete_carryover_mean = complete_toxicity_fit$carryover_mean,
        tite_carryover_mean = tite_toxicity_fit$carryover_mean
      ),
      data.frame(
        toxicity_model = toxicity_model,
        efficacy_model = efficacy_model,
        endpoint = "efficacy",
        dose = seq_len(ndose),
        complete_p_regular_mean = complete_efficacy_fit$p_regular_mean,
        tite_p_regular_mean = tite_efficacy_fit$p_regular_mean,
        absolute_difference = abs(complete_efficacy_fit$p_regular_mean - tite_efficacy_fit$p_regular_mean),
        complete_p_ipde_mean = complete_efficacy_fit$p_ipde_mean,
        tite_p_ipde_mean = tite_efficacy_fit$p_ipde_mean,
        complete_carryover_mean = complete_efficacy_fit$carryover_mean,
        tite_carryover_mean = tite_efficacy_fit$carryover_mean
      )
    )

    decision_comparisons[[comparison_id]] <- data.frame(
      toxicity_model = toxicity_model,
      efficacy_model = efficacy_model,
      complete_MTD = complete_final$MTD,
      tite_MTD = tite_final$MTD,
      complete_OBD = complete_final$OBD,
      tite_OBD = tite_final$OBD,
      complete_toxicity_recommendation = complete_decision$toxicity_recommended_dose,
      tite_toxicity_recommendation = tite_decision$toxicity_recommended_dose,
      complete_next_dose = complete_decision$next_dose,
      tite_next_dose = tite_decision$next_dose,
      complete_action = complete_decision$action,
      tite_action = tite_decision$action,
      complete_futile_doses = paste(which(complete_decision$efficacy_futile), collapse = ","),
      tite_futile_doses = paste(which(tite_decision$efficacy_futile), collapse = ","),
      stringsAsFactors = FALSE
    )
  }
}

posterior_comparisons <- do.call(rbind, posterior_comparisons)
decision_comparisons <- do.call(rbind, decision_comparisons)
utility_comparisons <- do.call(rbind, lapply(seq_len(nrow(decision_comparisons)), function(i) {
  toxicity_model <- decision_comparisons$toxicity_model[i]
  efficacy_model <- decision_comparisons$efficacy_model[i]
  data.frame(
    toxicity_model = toxicity_model,
    efficacy_model = efficacy_model,
    dose = seq_len(ndose),
    complete_utility = validation_final_selection(complete_admin, toxicity_fits[[toxicity_model]]$complete,
                                                   efficacy_fits[[efficacy_model]]$complete,
                                                   make_validation_config(toxicity_model, efficacy_model))$utilities,
    tite_utility = validation_final_selection(complete_admin, toxicity_fits[[toxicity_model]]$tite,
                                               efficacy_fits[[efficacy_model]]$tite,
                                               make_validation_config(toxicity_model, efficacy_model))$utilities
  )
}))
utility_comparisons$absolute_difference <- abs(utility_comparisons$complete_utility - utility_comparisons$tite_utility)

results_dir <- file.path(validation_dir, "results")
dir.create(results_dir, showWarnings = FALSE)
write.csv(complete_admin, file.path(results_dir, "complete_outcome_admin.csv"), row.names = FALSE)
write.csv(posterior_comparisons, file.path(results_dir, "complete_outcome_posterior_comparisons.csv"), row.names = FALSE)
write.csv(utility_comparisons, file.path(results_dir, "complete_outcome_utility_comparisons.csv"), row.names = FALSE)
write.csv(decision_comparisons, file.path(results_dir, "complete_outcome_decision_comparisons.csv"), row.names = FALSE)
saveRDS(list(toxicity_fits = toxicity_fits, efficacy_fits = efficacy_fits), file.path(results_dir, "complete_outcome_fits.rds"))

print(posterior_comparisons)
print(utility_comparisons)
print(decision_comparisons)
