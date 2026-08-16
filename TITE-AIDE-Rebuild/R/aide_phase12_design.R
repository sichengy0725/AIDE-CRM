aide_select_lowest_tied <- function(score, candidates) {
  if (!length(candidates)) return(NA_integer_)
  candidates[which.max(replace(score[candidates], is.na(score[candidates]), -Inf))][1L]
}

aide_current_mtd <- function(toxicity_fit, eliminated, target) {
  candidates <- which(!eliminated)
  if (!length(candidates)) return(NA_integer_)
  aide_select_lowest_tied(
    -abs(toxicity_fit$p_regular_mean - target), candidates
  )
}

aide_observed_efficacy_response <- function(state) {
  ndose <- length(state$eliminated)
  if (!nrow(state$admin)) return(rep(FALSE, ndose))
  tabulate(
    as.integer(state$admin$dose[
      state$admin$eff_final == 1L &
        is.finite(state$admin$t_response) &
        state$admin$t_response <= state$t_now
    ]),
    nbins = ndose
  ) > 0L
}

aide_allocation_probabilities <- function(utility, doses) {
  weights <- utility[doses]
  if (all(is.finite(weights)) && all(weights >= 0) && sum(weights) > 0) {
    return(weights / sum(weights))
  }
  rep(1 / length(doses), length(doses))
}

aide_build_admissible_set <- function(state, toxicity_fit, efficacy_fit, config, scenario) {
  ndose <- scenario$ndose; doses <- seq_len(ndose)
  excluded <- state$eliminated | state$futile
  ## Interim allocation never imposes a tried-dose restriction. No-skipping
  ## governs only one-stage upward movement, not the admissible set.
  list(
    doses = doses[!excluded],
    mtd = aide_current_mtd(
      toxicity_fit, state$eliminated, config$toxicity$target
    ),
    exclusion_reasons = which(excluded)
  )
}

aide_action_from_next <- function(next_dose, current_dose) {
  if (is.na(next_dose)) "stop" else if (next_dose > current_dose) "escalate" else if (next_dose < current_dose) "de_escalate" else "stay"
}

aide_one_stage_decision <- function(state, toxicity_recommendation, admissible, toxicity_fit, efficacy_fit, config) {
  d <- state$current_dose
  if (toxicity_recommendation$stop) {
    return(list(stop_trial = TRUE, stop_reason = "dose_1_overtoxicity", next_dose = NA_integer_))
  }
  candidates <- admissible$doses[admissible$doses <= admissible$mtd]
  if (!length(candidates)) {
    return(list(stop_trial = TRUE, stop_reason = "no_admissible_dose", next_dose = NA_integer_))
  }
  utility <- aide_compute_utility(efficacy_fit$p_regular_mean, toxicity_fit$p_regular_mean, config$utility)
  provisional_obd <- aide_select_lowest_tied(utility, candidates)

  ## Safety takes precedence over the utility rule. When the current dose is
  ## no longer safe, move to the highest globally admissible lower dose.
  if (toxicity_recommendation$action == "de_escalate") {
    lower <- candidates[candidates < d]
    if (!length(lower)) {
      return(list(stop_trial = TRUE, stop_reason = "no_safe_lower_one_stage_dose", next_dose = NA_integer_))
    }
    return(list(stop_trial = FALSE, stop_reason = NA_character_, next_dose = max(lower),
                candidate_doses = candidates, utility = utility,
                provisional_obd = provisional_obd))
  }

  if (provisional_obd > d) {
    highest_tried <- if (nrow(state$admin)) max(state$admin$dose) else d
    if (provisional_obd <= highest_tried) {
      next_dose <- provisional_obd
    } else {
      upward <- candidates[candidates > highest_tried & candidates <= provisional_obd]
      if (!length(upward)) {
        return(list(stop_trial = TRUE, stop_reason = "no_admissible_no_skip_escalation", next_dose = NA_integer_))
      }
      next_dose <- min(upward)
    }
  } else if (provisional_obd == d) {
    next_dose <- d
  } else {
    response_observed <- aide_observed_efficacy_response(state)
    ## For a utility target below the current dose, use the lowest
    ## response-qualified dose in the closed interval from the target to the
    ## current dose. Responses below the target cannot redirect allocation.
    response_qualified <- candidates[
      candidates >= provisional_obd &
        candidates <= d &
        response_observed[candidates]
    ]
    if (length(response_qualified)) {
      next_dose <- min(response_qualified)
    } else if (d %in% candidates) {
      next_dose <- d
    } else {
      lower <- candidates[candidates < d]
      if (!length(lower)) {
        return(list(stop_trial = TRUE, stop_reason = "no_safe_lower_one_stage_dose", next_dose = NA_integer_))
      }
      next_dose <- max(lower)
    }
  }
  list(stop_trial = FALSE, stop_reason = NA_character_, next_dose = next_dose,
       candidate_doses = candidates, utility = utility,
       provisional_obd = provisional_obd)
}

aide_two_stage_decision <- function(state, toxicity_recommendation, admissible, toxicity_fit, efficacy_fit, config) {
  stage1_n <- sum(state$admin$stage == "stage1")
  active <- state$stage$active
  transition <- FALSE
  if (active == "stage1" && stage1_n >= config$design$s1_Max && toxicity_recommendation$action == "stay") { active <- "stage2"; transition <- TRUE }
  if (toxicity_recommendation$stop) return(list(stop_trial = TRUE, stop_reason = "dose_1_overtoxicity", next_dose = NA_integer_, stage = active, stage_transition = transition))
  if (active == "stage1") {
    ## Stage I remains toxicity-directed, but a dose declared futile cannot be
    ## treated again.  This mirrors the non-TITE Stage-I admissible-dose rule:
    ## retain the toxicity recommendation when admissible; otherwise use the
    ## closest lower admissible dose on that recommendation path.
    current <- state$current_dose
    recommended <- as.integer(toxicity_recommendation$recommended_dose)
    if (recommended %in% admissible$doses) {
      return(list(stop_trial = FALSE, stop_reason = NA_character_, next_dose = recommended,
                  candidate_doses = recommended, stage = active, stage_transition = transition))
    }
    path <- if (recommended > current) {
      seq.int(current, min(recommended, length(state$eliminated)))
    } else {
      seq_len(max(1L, min(recommended, length(state$eliminated))))
    }
    candidates <- path[path %in% admissible$doses]
    if (!length(candidates)) {
      return(list(stop_trial = TRUE, stop_reason = "no_stage1_admissible_dose",
                  next_dose = NA_integer_, candidate_doses = integer(0),
                  stage = active, stage_transition = transition))
    }
    return(list(stop_trial = FALSE, stop_reason = NA_character_,
                next_dose = max(candidates), candidate_doses = candidates,
                stage = active, stage_transition = transition))
  }
  candidates <- admissible$doses[admissible$doses <= admissible$mtd]
  if (!length(candidates)) return(list(stop_trial = TRUE, stop_reason = "no_stage2_candidate", next_dose = NA_integer_, stage = active, stage_transition = transition))
  utility <- aide_compute_utility(efficacy_fit$p_regular_mean, toxicity_fit$p_regular_mean, config$utility)
  ranked <- candidates[order(-utility[candidates], candidates)]
  provisional_obd <- aide_select_lowest_tied(utility, candidates)
  response_observed <- aide_observed_efficacy_response(state)
  mtd_dose <- if (admissible$mtd %in% candidates) admissible$mtd else max(candidates)
  if (identical(config$monitoring$stage2_allocation, "top2_randomized")) {
    top_two <- ranked[seq_len(min(2L, length(ranked)))]
    if (length(top_two) == 1L) {
      allocation_doses <- top_two
    } else {
      responders_in_top_two <- top_two[response_observed[top_two]]
      allocation_doses <- if (length(responders_in_top_two) == 2L) {
        top_two
      } else if (length(responders_in_top_two) == 1L) {
        c(responders_in_top_two, mtd_dose)
      } else {
        lower_than_mtd <- candidates[candidates < mtd_dose]
        c(mtd_dose, if (length(lower_than_mtd)) max(lower_than_mtd) else mtd_dose)
      }
      allocation_doses <- unique(allocation_doses)
      if (length(allocation_doses) == 1L && length(candidates) > 1L) {
        alternatives <- candidates[
          candidates < mtd_dose & !(candidates %in% allocation_doses)
        ]
        if (!length(alternatives)) alternatives <- setdiff(candidates, allocation_doses)
        if (length(alternatives)) allocation_doses <- c(allocation_doses, max(alternatives))
      }
    }
    probabilities <- aide_allocation_probabilities(utility, allocation_doses)
    ## Freeze the candidate set and its probabilities when the triggering
    ## patient arrives. The event layer then randomizes each cohort position
    ## independently from this fixed set; it does not randomize one dose for
    ## the whole cohort. The MTD is retained as the Stage II reference dose
    ## for the decision log, safety recommendation, and n_eval gate.
    next_dose <- mtd_dose
    allocation_probabilities <- setNames(probabilities, paste0("D", allocation_doses))
    individual_randomization <- length(allocation_doses) > 1L
  } else {
    response_doses <- candidates[response_observed[candidates]]
    next_dose <- if (length(response_doses)) {
      aide_select_lowest_tied(utility, response_doses)
    } else {
      mtd_dose
    }
    allocation_probabilities <- setNames(1, paste0("D", next_dose))
    allocation_doses <- next_dose
    individual_randomization <- FALSE
  }
  list(stop_trial = FALSE, stop_reason = NA_character_, next_dose = next_dose,
       candidate_doses = candidates, utility = utility, stage = active,
       stage_transition = transition,
       provisional_obd = provisional_obd,
       stage2_allocation = config$monitoring$stage2_allocation,
       allocation_probabilities = allocation_probabilities,
       allocation_doses = allocation_doses,
       individual_randomization = individual_randomization)
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
       MTD = adm$mtd,
       provisional_OBD = raw$provisional_obd %||% NA_integer_,
       stage2_allocation = raw$stage2_allocation %||% "not_applicable",
       allocation_probabilities = raw$allocation_probabilities %||% numeric(0),
       allocation_doses = raw$allocation_doses %||% next_dose,
       individual_randomization = raw$individual_randomization %||% FALSE,
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
