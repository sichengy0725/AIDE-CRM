aide_select_lowest_tied <- function(score, candidates) {
  if (!length(candidates)) return(NA_integer_)
  candidates[which.max(replace(score[candidates], is.na(score[candidates]), -Inf))][1L]
}

aide_build_admissible_set <- function(state, toxicity_fit, efficacy_fit, config, scenario) {
  ndose <- scenario$ndose; doses <- seq_len(ndose)
  tried <- unique(state$admin$dose)
  excluded <- state$eliminated | state$futile
  if (config$monitoring$restrict_to_tried && length(tried)) excluded <- excluded | !(doses %in% tried)
  if (config$monitoring$no_skipping && length(tried)) excluded <- excluded | doses > max(c(state$current_dose, min(ndose, max(tried) + 1L)))
  list(doses = doses[!excluded], exclusion_reasons = which(excluded))
}

aide_action_from_next <- function(next_dose, current_dose) {
  if (is.na(next_dose)) "stop" else if (next_dose > current_dose) "escalate" else if (next_dose < current_dose) "de_escalate" else "stay"
}

aide_one_stage_decision <- function(state, toxicity_recommendation, admissible, toxicity_fit, efficacy_fit, config) {
  d <- state$current_dose
  if (toxicity_recommendation$stop) return(list(stop_trial = TRUE, stop_reason = "dose_1_overtoxicity", next_dose = NA_integer_))
  if (toxicity_recommendation$action == "de_escalate") return(list(stop_trial = FALSE, stop_reason = NA_character_, next_dose = toxicity_recommendation$recommended_dose, candidate_doses = toxicity_recommendation$recommended_dose))
  current_n <- sum(state$admin$dose == d)
  local <- if (toxicity_recommendation$action == "stay" && current_n >= config$design$m_U) c(d - 1L, d) else c(d - 1L, d, d + 1L)
  local <- intersect(admissible$doses, local[local >= 1L & local <= length(toxicity_fit$p_regular_mean)])
  if (!length(local)) return(list(stop_trial = TRUE, stop_reason = "no_admissible_dose", next_dose = NA_integer_))
  utility <- aide_compute_utility(efficacy_fit$p_regular_mean, toxicity_fit$p_regular_mean, config$utility)
  next_dose <- aide_select_lowest_tied(utility, local)
  list(stop_trial = FALSE, stop_reason = NA_character_, next_dose = next_dose, candidate_doses = local, utility = utility)
}

aide_two_stage_decision <- function(state, toxicity_recommendation, admissible, toxicity_fit, efficacy_fit, config) {
  stage1_n <- sum(state$admin$stage == "stage1")
  active <- state$stage$active
  transition <- FALSE
  if (active == "stage1" && stage1_n >= config$design$s1_Max && toxicity_recommendation$action == "stay") { active <- "stage2"; transition <- TRUE }
  if (toxicity_recommendation$stop) return(list(stop_trial = TRUE, stop_reason = "dose_1_overtoxicity", next_dose = NA_integer_, stage = active, stage_transition = transition))
  if (active == "stage1") return(list(stop_trial = FALSE, stop_reason = NA_character_, next_dose = toxicity_recommendation$recommended_dose, candidate_doses = toxicity_recommendation$recommended_dose, stage = active, stage_transition = transition))
  candidates <- admissible$doses
  if (config$monitoring$stage2_mtd_ceiling) candidates <- candidates[candidates <= state$current_dose]
  if (!length(candidates)) return(list(stop_trial = TRUE, stop_reason = "no_stage2_candidate", next_dose = NA_integer_, stage = active, stage_transition = transition))
  next_dose <- aide_select_lowest_tied(efficacy_fit$p_regular_mean, candidates)
  list(stop_trial = FALSE, stop_reason = NA_character_, next_dose = next_dose, candidate_doses = candidates, stage = active, stage_transition = transition)
}

aide_make_design_decision <- function(state, toxicity_fit, efficacy_fit, config, scenario) {
  state$eliminated <- state$eliminated | toxicity_fit$eliminated
  futility <- aide_update_futility(efficacy_fit, state$futile); state$futile <- futility$futile
  adm <- aide_build_admissible_set(state, toxicity_fit, efficacy_fit, config, scenario)
  tox_rec <- aide_toxicity_recommendation(state$current_dose, toxicity_fit, config, state$eliminated)
  raw <- if (config$design$allocation == "one_stage") aide_one_stage_decision(state, tox_rec, adm, toxicity_fit, efficacy_fit, config) else aide_two_stage_decision(state, tox_rec, adm, toxicity_fit, efficacy_fit, config)
  next_dose <- raw$next_dose %||% NA_integer_
  action <- aide_action_from_next(next_dose, state$current_dose)
  list(stop_trial = isTRUE(raw$stop_trial), stop_reason = raw$stop_reason %||% NA_character_,
       stage = raw$stage %||% state$stage$active, stage_transition = raw$stage_transition %||% FALSE,
       current_dose = state$current_dose, action = action, next_dose = next_dose,
       candidate_doses = raw$candidate_doses %||% integer(0), exclusion_reasons = adm$exclusion_reasons,
       p_toxicity = toxicity_fit$p_regular_mean, p_toxicity_ipde = toxicity_fit$p_ipde_mean,
       prob_toxicity_ipde_overdose = toxicity_fit$prob_ipde_overtox_by_dose,
       p_efficacy = efficacy_fit$p_regular_mean, p_efficacy_ipde = efficacy_fit$p_ipde_mean,
       utility = raw$utility %||% aide_compute_utility(efficacy_fit$p_regular_mean, toxicity_fit$p_regular_mean, config$utility),
       toxicity_action = tox_rec$action, toxicity_recommended_dose = tox_rec$recommended_dose,
       toxicity_eliminated = state$eliminated, efficacy_futile = state$futile,
       futility = futility)
}

aide_apply_n_eval_gate <- function(state, decision, config, trigger_type, trigger_queue_id = NA_integer_) {
  out <- aide_count_evaluated_current_dose(state$admin, state$current_dose, state$t_now, config$time$T_assess)
  # Latest design override: n_eval gates both a stay and an escalation.
  # De-escalation remains immediately permissible regardless of n_eval.
  blocked <- decision$action %in% c("stay", "escalate") && out$n < config$design$n_eval
  disposition <- "allowed"; new_dropped <- FALSE; retreat_retained <- FALSE
  if (blocked && trigger_type == "new") { disposition <- "new_dropped_n_eval"; new_dropped <- TRUE }
  if (blocked && trigger_type == "retreat") {
    disposition <- "retreat_retained_n_eval"; retreat_retained <- TRUE
    if (!is.na(trigger_queue_id)) { i <- match(trigger_queue_id, state$queue$queue_id); if (!is.na(i)) { state$queue$status[i] <- "waiting"; state$queue$last_reason[i] <- disposition } }
  }
  list(state = state, allowed = !blocked, blocked = blocked, n_eval_observed = out$n,
       qualifying_admin_ids = out$qualifying_admin_ids, trigger_disposition = disposition,
       new_patient_dropped = new_dropped, retreat_opportunity_retained = retreat_retained)
}

# Compatibility alias retained for callers of the prior rebuild draft.
aide_apply_escalation_gate <- aide_apply_n_eval_gate
