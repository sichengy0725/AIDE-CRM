aide_empty_df <- function(cols) as.data.frame(setNames(replicate(length(cols), logical(0), simplify = FALSE), cols), stringsAsFactors = FALSE)
aide_add_row <- function(x, row) { x[nrow(x) + 1L, names(row)] <- row; x }

aide_phase12_initialize_state <- function(config, scenario, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  aide_phase12_validate_config(config, scenario)
  state <- list(
    t_now = config$time$t0, active = TRUE, stop_reason = NA_character_,
    events = aide_empty_df(c("event_id", "time", "seq", "event_type", "patient_id", "admin_id")),
    queue = aide_empty_df(c("queue_id", "time", "seq", "type", "patient_id", "source_admin_id", "status", "last_reason")),
    admin = aide_empty_df(c("admin_id", "patient_id", "cohort_id", "decision_id", "stage", "assignment_type", "dose", "decision_next_dose", "previous_dose", "previous_admin_id", "cycle", "t_arrival", "t_start", "t_dlt", "t_response", "assessment_end", "dlt_final", "eff_final", "true_p_tox", "true_p_eff", "true_toxicity_state", "true_efficacy_state")),
    cohort = list(open = TRUE, cohort_id = 1L, opened_time = config$time$t0,
                  next_dose = config$design$start_dose, decision_action = "initial",
                  stage = if (config$design$allocation == "two_stage") "stage1" else "one_stage",
                  capacity = config$design$cohort_size, filled = 0L, decision_id = 0L,
                  allocation_doses = config$design$start_dose,
                  allocation_probabilities = setNames(1, paste0("D", config$design$start_dose)),
                  individual_randomization = FALSE, randomization_draws = 0L,
                  risk_context = NULL),
    current_dose = config$design$start_dose,
    stage = list(active = if (config$design$allocation == "two_stage") "stage1" else "one_stage", transitioned = FALSE),
    eliminated = rep(FALSE, scenario$ndose), futile = rep(FALSE, scenario$ndose),
    logs = list(event_log = aide_empty_df(c("event_id", "time", "seq", "event_type", "patient_id", "admin_id", "action")),
    decision_log = aide_empty_df(c("decision_id", "time", "trigger_event_id", "current_dose", "stage", "action", "next_dose", "MTD", "provisional_OBD", "stage2_allocation", "allocation_probabilities", "n_eval_blocked", "trigger_type", "trigger_disposition", "n_eval_required", "n_eval_observed", "stop_trial", "stop_reason")),
                cohort_log = aide_empty_df(c("cohort_id", "time", "decision_id", "stage", "action", "next_dose", "filled", "closed_time")),
                n_eval_log = aide_empty_df(c("decision_id", "time", "current_dose", "n_eval_required", "n_eval_observed", "evaluated_admin_ids", "blocked", "trigger_type", "trigger_disposition")),
                retreat_log = aide_empty_df(c("queue_id", "time", "patient_id", "source_admin_id", "eligible", "reason", "status", "cohort_id", "assigned_dose"))),
    counters = list(event = 0L, seq = 0L, queue = 0L, admin = 0L, patient = 0L, decision = 0L, cohort = 1L, arrivals = 0L)
  )
  state
}

aide_phase12_validate_state <- function(state, config) {
  if (isTRUE(state$cohort$open) && state$cohort$filled >= state$cohort$capacity) stop("Open cohort cannot be full.")
  if (nrow(state$admin) && any(state$admin$dose != state$admin$decision_next_dose)) stop("Cohort dose mismatch.")
  if (isTRUE(state$cohort$open) && nrow(state$admin) &&
      !isTRUE(state$cohort$individual_randomization)) {
    id <- state$admin$cohort_id == state$cohort$cohort_id
    if (any(state$admin$dose[id] != state$cohort$next_dose)) stop("Open cohort is not dose homogeneous.")
  }
  if (isTRUE(state$cohort$individual_randomization) && nrow(state$admin)) {
    id <- state$admin$cohort_id == state$cohort$cohort_id
    if (any(!state$admin$dose[id] %in% state$cohort$allocation_doses))
      stop("Randomized cohort contains a dose outside its frozen candidate set.")
  }
  invisible(TRUE)
}
