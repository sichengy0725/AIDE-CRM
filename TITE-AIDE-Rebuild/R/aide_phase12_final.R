aide_final_analysis <- function(state, config, scenario) {
  if (!nrow(state$admin)) return(list(MTD = NA_integer_, OBD = NA_integer_, no_selection_reason = "no_administrations",
                                       elimination = state$eliminated, futility = state$futile, utilities = rep(NA_real_, scenario$ndose)))
  t_final <- max(state$admin$assessment_end)
  interim <- aide_build_interim_data(state, t_final, config)
  tox <- aide_fit_toxicity(interim$toxicity, config, scenario$ndose)
  eff <- aide_fit_efficacy(interim$efficacy, config, scenario$ndose)
  eliminated <- state$eliminated | tox$eliminated; futile <- state$futile | eff$futile
  tried <- sort(unique(state$admin$dose)); tox_candidates <- tried[!eliminated[tried]]
  if (!length(tox_candidates)) return(list(MTD = NA_integer_, OBD = NA_integer_, no_selection_reason = "no_safe_tried_dose",
    elimination = eliminated, futility = futile, utilities = aide_compute_utility(eff$p_regular_mean, tox$p_regular_mean, config$utility),
    toxicity = tox, efficacy = eff))
  MTD <- tox_candidates[which.min(abs(tox$p_regular_mean[tox_candidates] - config$toxicity$target))]
  obd_candidates <- tox_candidates[!futile[tox_candidates] & tox_candidates <= MTD]
  OBD <- if (length(obd_candidates)) aide_select_lowest_tied(eff$p_regular_mean, obd_candidates) else NA_integer_
  list(MTD = MTD, OBD = OBD, no_selection_reason = if (is.na(OBD)) "all_safe_doses_futile" else NA_character_,
       elimination = eliminated, futility = futile, utilities = aide_compute_utility(eff$p_regular_mean, tox$p_regular_mean, config$utility),
       toxicity = tox, efficacy = eff, final_time = t_final)
}

aide_summarize_trial <- function(state, config, scenario) {
  final <- aide_final_analysis(state, config, scenario)
  list(config = config, scenario = scenario, admin = state$admin, event_log = state$logs$event_log,
       queue_log = state$queue, cohort_log = state$logs$cohort_log, n_eval_log = state$logs$n_eval_log,
       decision_log = state$logs$decision_log, retreat_log = state$logs$retreat_log,
       stage1 = list(n = sum(state$admin$stage == "stage1"), transitioned = state$stage$transitioned),
       stage2 = list(n = sum(state$admin$stage == "stage2")), final = final,
       events_remaining = state$events, queue_remaining = state$queue[state$queue$status == "waiting", , drop = FALSE],
       stop_reason = state$stop_reason)
}

simulate_AIDE_phase_I_II <- function(config = NULL, scenario = NULL, p_true = NULL, e_true = NULL, ...) {
  dots <- list(...)
  deprecated <- intersect(names(dots), c("ipde_design", "flexible_ipde", "random_crm_recycle_efficacy_delta", "enrollment_scheme"))
  if (length(deprecated)) stop("Deprecated rebuild controls are not supported: ", paste(deprecated, collapse = ", "), ".")
  if (is.null(config)) config <- do.call(aide_phase12_config, dots)
  if (is.null(scenario)) scenario <- aide_phase12_scenario(p_true, e_true)
  simulate_AIDE_phase_I_II_event(config, scenario, seed = dots$seed %||% 1L)
}
