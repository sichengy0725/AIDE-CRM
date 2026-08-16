aide_event_schedule <- function(state, time, event_type, patient_id = NA_integer_, admin_id = NA_integer_) {
  state$counters$event <- state$counters$event + 1L; state$counters$seq <- state$counters$seq + 1L
  state$events <- aide_add_row(state$events, list(event_id = state$counters$event, time = as.numeric(time), seq = state$counters$seq,
                                                   event_type = event_type, patient_id = patient_id, admin_id = admin_id))
  state
}

aide_event_pop_next <- function(state) {
  if (!nrow(state$events)) return(list(event = NULL, state = state))
  i <- order(state$events$time, state$events$seq)[1L]; event <- state$events[i, , drop = FALSE]
  state$events <- state$events[-i, , drop = FALSE]
  list(event = event, state = state)
}

aide_schedule_new_arrival <- function(state, config) {
  if (state$counters$arrivals >= config$design$Nmax * 10L) return(state)
  state$counters$arrivals <- state$counters$arrivals + 1L; state$counters$patient <- state$counters$patient + 1L
  arrival_time <- if (state$counters$arrivals == 1L) config$time$t0 else state$t_now + stats::rexp(1L, config$time$arrival_rate)
  aide_event_schedule(state, arrival_time, "new_arrival", state$counters$patient)
}

aide_queue_add <- function(state, time, seq, type, patient_id, source_admin_id = NA_integer_, reason = "waiting") {
  duplicate <- nrow(state$queue) && any(state$queue$type == type & state$queue$patient_id == patient_id &
                                        (if (type == "new") TRUE else state$queue$source_admin_id == source_admin_id) &
                                        state$queue$status %in% c("waiting", "assigned"))
  if (duplicate) return(state)
  state$counters$queue <- state$counters$queue + 1L
  state$queue <- aide_add_row(state$queue, list(queue_id = state$counters$queue, time = time, seq = seq, type = type,
                                                patient_id = patient_id, source_admin_id = source_admin_id,
                                                status = "waiting", last_reason = reason))
  state
}

aide_retreat_structural_eligibility <- function(admin_row, t_now, config) {
  if (nrow(admin_row) != 1L) stop("Retreat eligibility requires one administration.")
  if (admin_row$assessment_end > t_now) return(list(eligible = FALSE, reason = "assessment_pending"))
  if (admin_row$dlt_final == 1L) return(list(eligible = FALSE, reason = "dlt"))
  if (admin_row$eff_final == 1L) return(list(eligible = FALSE, reason = "efficacy_response"))
  if (admin_row$cycle >= config$design$cycle_max) return(list(eligible = FALSE, reason = "cycle_max"))
  list(eligible = TRUE, reason = "eligible")
}

aide_queue_retreat_opportunity <- function(state, admin_id, config) {
  i <- match(admin_id, state$admin$admin_id); if (is.na(i)) return(list(state = state, queue_id = NA_integer_, eligible = FALSE))
  row <- state$admin[i, , drop = FALSE]; check <- aide_retreat_structural_eligibility(row, state$t_now, config)
  qid <- NA_integer_
  if (check$eligible) {
    before <- nrow(state$queue)
    state <- aide_queue_add(state, state$t_now, state$counters$seq, "retreat", row$patient_id, row$admin_id, check$reason)
    if (nrow(state$queue) > before) qid <- state$queue$queue_id[nrow(state$queue)] else qid <- state$queue$queue_id[which(state$queue$source_admin_id == row$admin_id)[1L]]
  }
  state$logs$retreat_log <- aide_add_row(state$logs$retreat_log, list(queue_id = qid, time = state$t_now, patient_id = row$patient_id,
                                                                      source_admin_id = row$admin_id, eligible = check$eligible,
                                                                      reason = check$reason, status = if (check$eligible) "waiting" else "removed_ineligible",
                                                                      cohort_id = NA_integer_, assigned_dose = NA_integer_))
  list(state = state, queue_id = qid, eligible = check$eligible)
}

aide_queue_next_new <- function(state) { i <- which(state$queue$type == "new" & state$queue$status == "waiting"); if (!length(i)) integer(0) else i[order(state$queue$time[i], state$queue$seq[i])][1L] }

aide_cohort_allocation <- function(state) {
  doses <- as.integer(state$cohort$allocation_doses %||% state$cohort$next_dose)
  doses <- unique(doses[is.finite(doses)])
  probabilities <- as.numeric(state$cohort$allocation_probabilities %||% numeric(0))
  if (length(probabilities) != length(doses) || any(!is.finite(probabilities)) ||
      any(probabilities < 0) || sum(probabilities) <= 0) {
    probabilities <- rep(1 / length(doses), length(doses))
  } else {
    probabilities <- probabilities / sum(probabilities)
  }
  list(doses = doses, probabilities = probabilities)
}

aide_retreat_individual_risk_screen <- function(state, queue_index, config,
                                                 candidate_doses = NULL) {
  if (is.null(candidate_doses)) candidate_doses <- aide_cohort_allocation(state)$doses
  candidate_doses <- as.integer(candidate_doses)
  if (!isTRUE(state$cohort$open) || is.null(state$cohort$risk_context)) {
    return(list(allowed = length(candidate_doses) > 0L, reason = "eligible",
                eligible_doses = candidate_doses))
  }
  risk <- state$cohort$risk_context
  eligible <- candidate_doses
  if (isTRUE(config$recycle$apply_individual_toxicity_risk)) {
    prob <- risk$prob_toxicity_ipde_overdose[eligible]
    eligible <- eligible[is.finite(prob) & prob <= config$recycle$toxicity_ipde_overdose_cutoff]
    if (!length(eligible)) {
      return(list(allowed = FALSE, reason = "blocked_individual_toxicity_risk",
                  eligible_doses = integer(0), toxicity_overdose_probability = prob))
    }
  }
  if (isTRUE(config$recycle$apply_individual_efficacy_benefit)) {
    benefit <- risk$p_efficacy_ipde[eligible] - risk$p_efficacy[eligible]
    eligible <- eligible[is.finite(benefit) & benefit > config$recycle$efficacy_ipde_min_increment]
    if (!length(eligible)) {
      return(list(allowed = FALSE, reason = "blocked_individual_efficacy_benefit",
                  eligible_doses = integer(0), efficacy_increment = benefit))
    }
  }
  list(allowed = TRUE, reason = "eligible", eligible_doses = eligible)
}

aide_queue_next_retreat <- function(state, config) {
  i <- which(state$queue$type == "retreat" & state$queue$status == "waiting")
  if (!length(i)) return(list(state = state, index = integer(0)))
  for (k in i[order(state$queue$time[i], state$queue$seq[i])]) {
    a <- state$admin[match(state$queue$source_admin_id[k], state$admin$admin_id), , drop = FALSE]
    e <- aide_retreat_structural_eligibility(a, state$t_now, config)
    if (e$eligible) {
      risk <- aide_retreat_individual_risk_screen(state, k, config)
      if (risk$allowed) return(list(state = state, index = k))
      if (state$queue$last_reason[k] != risk$reason) {
        state$queue$last_reason[k] <- risk$reason
        state$logs$retreat_log <- aide_add_row(state$logs$retreat_log, list(queue_id = state$queue$queue_id[k], time = state$t_now,
          patient_id = state$queue$patient_id[k], source_admin_id = state$queue$source_admin_id[k], eligible = FALSE,
          reason = risk$reason, status = "waiting", cohort_id = state$cohort$cohort_id, assigned_dose = state$cohort$next_dose))
      }
      next
    }
    state$queue$status[k] <- "removed_ineligible"; state$queue$last_reason[k] <- e$reason
  }
  list(state = state, index = integer(0))
}

aide_select_assignment_opportunity <- function(state, config) {
  i <- aide_queue_next_new(state); if (length(i)) return(list(state = state, index = i))
  aide_queue_next_retreat(state, config)
}

aide_assignment_dose <- function(state, queue_index, config) {
  allocation <- aide_cohort_allocation(state)
  doses <- allocation$doses
  probabilities <- allocation$probabilities
  if (!length(doses)) stop("Open cohort has no allocation dose.")
  if (state$queue$type[queue_index] == "retreat") {
    risk <- aide_retreat_individual_risk_screen(state, queue_index, config, doses)
    if (!risk$allowed) stop("Attempted to assign a retreat patient without an eligible cohort dose.")
    keep <- match(risk$eligible_doses, doses)
    doses <- risk$eligible_doses
    probabilities <- probabilities[keep]
    probabilities <- probabilities / sum(probabilities)
  }
  if (!isTRUE(state$cohort$individual_randomization)) return(doses[1L])
  sample(doses, size = 1L, prob = probabilities)
}

aide_append_administration <- function(state, queue_index, config, scenario) {
  q <- state$queue[queue_index, , drop = FALSE]; is_retreat <- q$type == "retreat"
  previous <- NA_integer_; cycle <- 1L; arrival <- q$time
  if (is_retreat) { prev <- state$admin[match(q$source_admin_id, state$admin$admin_id), , drop = FALSE]; previous <- prev$dose; cycle <- prev$cycle + 1L; arrival <- prev$t_arrival }
  ## For Stage II top-two allocation, this draw is made for this individual
  ## only. The candidate set and probabilities remain frozen for the cohort.
  dose <- aide_assignment_dose(state, queue_index, config)
  if (isTRUE(state$cohort$individual_randomization)) {
    state$cohort$randomization_draws <- state$cohort$randomization_draws + 1L
  }
  p_tox <- if (is_retreat) aide_phase12_ipde_probability(scenario$p_true, scenario$toxicity_ipde_dgm$alpha_true, dose) else scenario$p_true[dose]
  p_eff <- if (is_retreat) aide_phase12_ipde_probability(scenario$e_true, scenario$efficacy_ipde_dgm$alpha_true, dose) else scenario$e_true[dose]
  tox <- aide_draw_binary_event_time(p_tox, state$t_now, config$time$T_assess, config$time$dlt_dist)
  eff <- aide_draw_binary_event_time(p_eff, state$t_now, config$time$T_assess, config$time$efficacy_dist)
  state$counters$admin <- state$counters$admin + 1L; id <- state$counters$admin
  state$admin <- aide_add_row(state$admin, list(admin_id = id, patient_id = q$patient_id, cohort_id = state$cohort$cohort_id,
    decision_id = state$cohort$decision_id, stage = state$cohort$stage, assignment_type = q$type, dose = dose,
    decision_next_dose = dose, previous_dose = previous, cycle = cycle, t_arrival = arrival, t_start = state$t_now,
    t_dlt = tox$event_time, t_response = eff$event_time, assessment_end = tox$window_end, dlt_final = tox$outcome,
    eff_final = eff$outcome, true_p_tox = p_tox, true_p_eff = p_eff))
  if (is.finite(tox$event_time)) state <- aide_event_schedule(state, tox$event_time, "dlt_observed", q$patient_id, id)
  if (is.finite(eff$event_time)) state <- aide_event_schedule(state, eff$event_time, "efficacy_response", q$patient_id, id)
  state <- aide_event_schedule(state, tox$window_end, "assessment_complete", q$patient_id, id)
  state$queue$status[queue_index] <- "assigned"; state$queue$last_reason[queue_index] <- "assigned"
  state$cohort$filled <- state$cohort$filled + 1L
  state
}

aide_fill_open_cohort <- function(state, config, scenario) {
  while (state$cohort$open && state$cohort$filled < state$cohort$capacity && nrow(state$admin) < config$design$Nmax) {
    selected <- aide_select_assignment_opportunity(state, config); state <- selected$state; i <- selected$index; if (!length(i)) break
    state <- aide_append_administration(state, i, config, scenario)
  }
  if (state$cohort$filled >= state$cohort$capacity) {
    state$cohort$open <- FALSE
    state$logs$cohort_log <- aide_add_row(state$logs$cohort_log, list(cohort_id = state$cohort$cohort_id, time = state$t_now,
      decision_id = state$cohort$decision_id, stage = state$cohort$stage, action = state$cohort$decision_action,
      next_dose = state$cohort$next_dose, filled = state$cohort$filled, closed_time = state$t_now))
  }
  stage2_dose_cap <- state$cohort$stage == "stage2" &&
    any(tabulate(state$admin$dose[state$admin$stage == "stage2"],
                 nbins = length(state$eliminated)) >= config$design$N_s2)
  if (stage2_dose_cap) {
    state$active <- FALSE; state$stop_reason <- "stage2_dose_cap"
  } else if (nrow(state$admin) >= config$design$Nmax) {
    state$active <- FALSE; state$stop_reason <- "administration_cap"
  }
  state
}

aide_open_cohort <- function(state, decision) {
  state$counters$cohort <- state$counters$cohort + 1L
  allocation_doses <- as.integer(decision$allocation_doses %||% decision$next_dose)
  allocation_probabilities <- as.numeric(decision$allocation_probabilities %||% 1)
  if (length(allocation_probabilities) != length(allocation_doses) ||
      any(!is.finite(allocation_probabilities)) || any(allocation_probabilities < 0) ||
      sum(allocation_probabilities) <= 0) {
    allocation_probabilities <- rep(1 / length(allocation_doses), length(allocation_doses))
  } else {
    allocation_probabilities <- allocation_probabilities / sum(allocation_probabilities)
  }
  state$cohort <- list(open = TRUE, cohort_id = state$counters$cohort, opened_time = state$t_now, next_dose = decision$next_dose,
    decision_action = decision$action, stage = decision$stage, capacity = state$cohort$capacity, filled = 0L,
    decision_id = state$counters$decision, allocation_doses = allocation_doses,
    allocation_probabilities = allocation_probabilities,
    individual_randomization = isTRUE(decision$individual_randomization), randomization_draws = 0L)
  state$cohort$risk_context <- list(p_efficacy = decision$p_efficacy, p_efficacy_ipde = decision$p_efficacy_ipde,
    prob_toxicity_ipde_overdose = decision$prob_toxicity_ipde_overdose)
  state$current_dose <- decision$next_dose; state$stage$active <- decision$stage; state$stage$transitioned <- state$stage$transitioned || isTRUE(decision$stage_transition)
  state
}

aide_record_decision <- function(state, decision, event, gate, config) {
  allocation_probabilities <- decision$allocation_probabilities %||% numeric(0)
  allocation_text <- if (length(allocation_probabilities)) {
    paste(paste0(names(allocation_probabilities), "=", formatC(allocation_probabilities,
          format = "f", digits = 3L)), collapse = ";")
  } else {
    NA_character_
  }
  state$logs$decision_log <- aide_add_row(state$logs$decision_log, list(decision_id = state$counters$decision, time = state$t_now,
    trigger_event_id = event$event_id, current_dose = decision$current_dose, stage = decision$stage, action = decision$action,
    next_dose = decision$next_dose, MTD = decision$MTD %||% NA_integer_,
    provisional_OBD = decision$provisional_OBD %||% NA_integer_,
    stage2_allocation = decision$stage2_allocation %||% NA_character_,
    allocation_probabilities = allocation_text, n_eval_blocked = gate$blocked, trigger_type = ifelse(event$event_type == "new_arrival", "new", "retreat"),
    trigger_disposition = gate$trigger_disposition, n_eval_required = config$design$n_eval, n_eval_observed = gate$n_eval_observed,
    stop_trial = decision$stop_trial, stop_reason = decision$stop_reason))
  state
}

simulate_AIDE_phase_I_II_event <- function(config, scenario, seed = 1L) {
  state <- aide_phase12_initialize_state(config, scenario, seed)
  state <- aide_schedule_new_arrival(state, config)
  while (state$active && nrow(state$events)) {
    popped <- aide_event_pop_next(state); event <- popped$event; state <- popped$state; state$t_now <- event$time
    state$logs$event_log <- aide_add_row(state$logs$event_log, list(event_id = event$event_id, time = event$time, seq = event$seq,
      event_type = event$event_type, patient_id = event$patient_id, admin_id = event$admin_id, action = "processed"))
    trigger_type <- NA_character_; trigger_queue_id <- NA_integer_
    if (event$event_type == "new_arrival") {
      state <- aide_schedule_new_arrival(state, config)
      trigger_type <- "new"
    } else if (event$event_type == "assessment_complete") {
      retreat <- aide_queue_retreat_opportunity(state, event$admin_id, config); state <- retreat$state
      if (retreat$eligible) { trigger_type <- "retreat"; trigger_queue_id <- retreat$queue_id }
    }
    # Only new arrivals and newly eligible retreat opportunities are assignable events.
    if (is.na(trigger_type)) next
    if (state$cohort$open) {
      if (trigger_type == "new") state <- aide_queue_add(state, state$t_now, event$seq, "new", event$patient_id)
      state <- aide_fill_open_cohort(state, config, scenario)
      aide_phase12_validate_state(state, config)
      next
    }
    interim <- aide_build_interim_data(state, state$t_now, config)
    tox_fit <- aide_fit_toxicity(interim$toxicity, config, scenario$ndose)
    eff_fit <- aide_fit_efficacy(interim$efficacy, config, scenario$ndose)
    decision <- aide_make_design_decision(state, tox_fit, eff_fit, config, scenario)
    state$eliminated <- decision$toxicity_eliminated; state$futile <- decision$efficacy_futile
    state$counters$decision <- state$counters$decision + 1L
    gate <- aide_apply_n_eval_gate(state, decision, config, trigger_type, trigger_queue_id); state <- gate$state
    state <- aide_record_decision(state, decision, event, gate, config)
    if (decision$action %in% c("stay", "escalate")) state$logs$n_eval_log <- aide_add_row(state$logs$n_eval_log, list(
      decision_id = state$counters$decision, time = state$t_now, current_dose = decision$current_dose,
      n_eval_required = config$design$n_eval, n_eval_observed = gate$n_eval_observed,
      evaluated_admin_ids = paste(gate$qualifying_admin_ids, collapse = ";"), blocked = gate$blocked,
      trigger_type = trigger_type, trigger_disposition = gate$trigger_disposition))
    if (decision$stop_trial) { state$active <- FALSE; state$stop_reason <- decision$stop_reason; break }
    if (is.na(decision$next_dose)) next
    if (!gate$allowed) next
    if (trigger_type == "new") state <- aide_queue_add(state, state$t_now, event$seq, "new", event$patient_id)
    state <- aide_open_cohort(state, decision)
    state <- aide_fill_open_cohort(state, config, scenario)
    aide_phase12_validate_state(state, config)
  }
  if (is.na(state$stop_reason)) state$stop_reason <- if (nrow(state$events)) "stopped" else "no_assignable_events"
  aide_summarize_trial(state, config, scenario)
}
