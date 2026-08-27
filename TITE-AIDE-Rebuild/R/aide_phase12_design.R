aide_select_lowest_tied <- function(score, candidates) {
  if (!length(candidates)) return(NA_integer_)
  candidates[which.max(replace(score[candidates], is.na(score[candidates]), -Inf))][1L]
}

aide_current_mtd <- function(toxicity_fit, eliminated, target) {
  ## A dose-1 toxicity finding stops the trial rather than creating a
  ## dose-specific toxicity admissibility set.
  if (isTRUE(eliminated[1L])) return(NA_integer_)
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

aide_one_stage_rule <- function(state, admissible, toxicity_fit, efficacy_fit, config) {
  d <- state$current_dose
  candidates <- admissible$doses[admissible$doses <= admissible$mtd]
  if (!length(candidates)) {
    return(list(stop_trial = TRUE, stop_reason = "no_admissible_dose",
                next_dose = NA_integer_))
  }
  utility <- aide_compute_utility(
    efficacy_fit$p_regular_mean, toxicity_fit$p_regular_mean, config$utility
  )
  response_observed <- aide_observed_efficacy_response(state)
  response_doses <- candidates[response_observed[candidates]]
  if (length(response_doses)) {
    response_floor <- min(response_doses)
    response_candidates <- candidates[candidates >= response_floor]
    target <- aide_select_lowest_tied(utility, response_candidates)
    target_action <- "response_floor_utility_target"
  } else {
    response_floor <- NA_integer_
    response_candidates <- candidates
    ## The MTD is the no-response target. If a standing elimination removes
    ## it, retain the highest remaining admissible dose below the MTD ceiling.
    target <- if (admissible$mtd %in% candidates) {
      admissible$mtd
    } else {
      max(candidates)
    }
    target_action <- "mtd_no_response_fallback"
  }

  allocation_mode <- config$design$allocation_mode
  if (allocation_mode == "upward_no_skipping" && target <= d) {
    next_dose <- target
    action <- target_action
  } else if (allocation_mode == "one_level_toward_obd") {
    if (target == d) {
      next_dose <- target
      action <- target_action
    } else {
      path <- if (target > d) {
        seq.int(d + 1L, target)
      } else {
        seq.int(d - 1L, target, by = -1L)
      }
      eligible_path <- path[path %in% admissible$doses]
      if (!length(eligible_path)) {
        return(list(stop_trial = TRUE,
                    stop_reason = "no_admissible_path_to_one_stage_target",
                    next_dose = NA_integer_))
      }
      next_dose <- eligible_path[1L]
      bypassed <- next_dose != path[1L]
      action <- if (target > d) {
        if (bypassed) "one_level_escalate_bypass_excluded" else "one_level_escalate"
      } else {
        if (bypassed) "one_level_deescalate_bypass_excluded" else "one_level_deescalate"
      }
    }
  } else {
    ## Inspect every dose on the upward path.  Using the highest tried dose
    ## here would incorrectly allow a previously tried high dose to bypass an
    ## intervening dose that has never been tried.
    path <- seq.int(d + 1L, target)
    tried <- unique(as.integer(state$admin$dose))
    first_untried <- path[path %in% admissible$doses & !(path %in% tried)]
    next_dose <- if (length(first_untried)) min(first_untried) else target
    action <- if (next_dose == target) target_action else "no_skip_untried_dose"
  }
  list(
    stop_trial = FALSE,
    stop_reason = NA_character_,
    next_dose = next_dose,
    target_dose = target,
    candidate_doses = response_candidates,
    utility = utility,
    provisional_obd = target,
    response_floor = response_floor,
    allocation_mode = allocation_mode,
    allocation_action = action
  )
}

aide_one_stage_decision <- function(state, toxicity_recommendation, admissible, toxicity_fit, efficacy_fit, config) {
  if (toxicity_recommendation$stop) {
    return(list(stop_trial = TRUE, stop_reason = "dose_1_overtoxicity", next_dose = NA_integer_))
  }
  aide_one_stage_rule(state, admissible, toxicity_fit, efficacy_fit, config)
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

  ## Stage II is deliberately not a separate allocation design.  After the
  ## toxicity-driven Stage I transition, it uses the one-stage rule above,
  ## including the response floor, MTD fallback, and selected allocation mode.
  ## No Stage-II-specific randomization or utility mode is used.
  stage2 <- aide_one_stage_decision(
    state, toxicity_recommendation, admissible, toxicity_fit, efficacy_fit,
    config
  )
  stage2$stage <- active
  stage2$stage_transition <- transition
  stage2$stage2_allocation <- "one_stage"
  stage2$allocation_doses <- stage2$next_dose
  stage2$allocation_probabilities <- if (is.na(stage2$next_dose)) {
    numeric(0)
  } else {
    setNames(1, paste0("D", stage2$next_dose))
  }
  stage2$individual_randomization <- FALSE
  return(stage2)

}

aide_make_design_decision <- function(state, toxicity_fit, efficacy_fit, config, scenario) {
  ## toxicity_fit$eliminated can mark dose 1 only. Efficacy futility remains
  ## the sole dose-by-dose elimination mechanism.
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
        futility = futility,
        toxicity_fit = toxicity_fit,
        efficacy_fit = efficacy_fit)
}

aide_apply_n_eval_gate <- function(state, decision, config, trigger_type, trigger_queue_id = NA_integer_) {
  out <- aide_count_evaluated_current_dose(state$admin, state$current_dose, state$t_now, config$time$T_assess)
  ## With no open cohort, the revised event clock checks this condition before
  ## fitting either model.  It gates the decision itself (not merely upward
  ## movement), and all waiting events remain in the queue when it is unmet.
  blocked <- out$n < config$design$n_eval
  disposition <- if (blocked) "decision_unavailable_n_eval" else "allowed"
  list(state = state, allowed = !blocked, blocked = blocked, n_eval_observed = out$n,
       qualifying_admin_ids = out$qualifying_admin_ids, trigger_disposition = disposition,
       new_patient_dropped = FALSE, retreat_opportunity_retained = blocked && trigger_type == "retreat")
}

# Compatibility alias retained for callers of the prior rebuild draft.
aide_apply_escalation_gate <- aide_apply_n_eval_gate
