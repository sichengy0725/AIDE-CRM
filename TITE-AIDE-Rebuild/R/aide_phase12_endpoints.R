aide_draw_binary_event_time <- function(probability, start_time, window, dist = "uniform") {
  outcome <- stats::rbinom(1L, 1L, probability)
  if (!outcome) return(list(outcome = 0L, event_time = Inf, window_end = start_time + window))
  delay <- if (identical(dist, "front_loaded")) stats::rbeta(1L, 1, 3) * window else stats::runif(1L, 0, window)
  list(outcome = 1L, event_time = start_time + delay, window_end = start_time + window)
}

aide_phase12_ipde_probability <- function(p, alpha, dose) {
  lower <- if (dose <= 1L) 0 else p[dose - 1L]
  min(1, p[dose] + alpha * lower)
}

## Alternative truth mechanisms from the August 31 specification.  The
## settings are shared by toxicity and efficacy; only their regular scenario
## probability vectors differ.
aide_phase12_tite_true_generation_settings <- function(
    model,
    ndose,
    dose_specific_alpha = NULL,
    random_effect_eta = 1,
    effective_dose_values = NULL,
    effective_dose_alpha = 0) {
  choices <- c(
    "legacy",
    "shared_multicycle",
    "dose_specific_geometric",
    "shared_patient_logistic",
    "effective_dose_geometric"
  )
  model <- match.arg(model, choices)
  if (is.null(dose_specific_alpha)) {
    ## Gives (0.2, 0.3, ..., 0.6) for the five-dose example in the PDF.
    dose_specific_alpha <- seq(0.2, 0.6, length.out = ndose)
  }
  dose_specific_alpha <- as.numeric(dose_specific_alpha)
  if (length(dose_specific_alpha) != ndose ||
      any(!is.finite(dose_specific_alpha)) || any(dose_specific_alpha < 0)) {
    stop("true_dose_specific_alpha must be non-negative with one value per dose.")
  }
  if (length(random_effect_eta) != 1L || !is.finite(random_effect_eta) ||
      random_effect_eta <= 0) {
    stop("true_random_effect_eta must be one positive finite value.")
  }
  if (is.null(effective_dose_values)) {
    effective_dose_values <- if (ndose == 5L) {
      c(15, 20, 30, 35, 45)
    } else {
      seq_len(ndose)
    }
  }
  effective_dose_values <- as.numeric(effective_dose_values)
  if (length(effective_dose_values) != ndose ||
      any(!is.finite(effective_dose_values)) ||
      any(diff(effective_dose_values) <= 0)) {
    stop("true_effective_dose_values must be strictly increasing with one value per dose.")
  }
  if (length(effective_dose_alpha) != 1L ||
      !is.finite(effective_dose_alpha) || effective_dose_alpha < 0) {
    stop("true_effective_dose_alpha must be one non-negative finite value.")
  }
  list(
    model = model,
    dose_specific_alpha = dose_specific_alpha,
    random_effect_eta = as.numeric(random_effect_eta),
    effective_dose_values = effective_dose_values,
    effective_dose_alpha = as.numeric(effective_dose_alpha)
  )
}

aide_phase12_tite_true_endpoint_probability <- function(
    settings,
    regular_probability,
    current_dose,
    history_doses = integer(0),
    patient_random_effect = 0,
    legacy_alpha = 0,
    legacy_multicycle = FALSE,
    previous_state = 0) {
  current_dose <- as.integer(current_dose)
  history_doses <- as.integer(history_doses)
  if (settings$model == "legacy") {
    state <- regular_probability[current_dose] + legacy_alpha * previous_state
    probability <- if (!length(history_doses)) {
      regular_probability[current_dose]
    } else if (legacy_multicycle) {
      min(1, state)
    } else {
      aide_phase12_ipde_probability(
        regular_probability, legacy_alpha, current_dose
      )
    }
    return(list(state = state, probability = probability,
                effective_dose = NA_real_))
  }
  if (settings$model == "shared_multicycle") {
    state <- regular_probability[current_dose] + legacy_alpha * previous_state
    return(list(state = state, probability = min(1, state),
                effective_dose = NA_real_))
  }
  if (settings$model == "dose_specific_geometric") {
    lag <- rev(seq_along(history_doses))
    state <- regular_probability[current_dose] + if (length(history_doses)) {
      sum(settings$dose_specific_alpha[history_doses] ^ lag *
            regular_probability[history_doses])
    } else {
      0
    }
    return(list(state = state, probability = min(1, state),
                effective_dose = NA_real_))
  }
  if (settings$model == "shared_patient_logistic") {
    probability <- stats::plogis(
      stats::qlogis(regular_probability[current_dose]) + patient_random_effect
    )
    return(list(state = probability, probability = probability,
                effective_dose = NA_real_))
  }
  lag <- rev(seq_along(history_doses))
  effective_dose <- settings$effective_dose_values[current_dose] +
    if (length(history_doses)) {
      sum(settings$effective_dose_alpha ^ lag *
            settings$effective_dose_values[history_doses])
    } else {
      0
  }
  dose_values <- settings$effective_dose_values
  ndose <- length(dose_values)
  if (ndose == 1L || effective_dose <= dose_values[1L]) {
    probability <- regular_probability[1L]
  } else if (effective_dose >= dose_values[ndose]) {
    k <- ndose - 1L
    probability <- regular_probability[k] +
      (effective_dose - dose_values[k]) /
        (dose_values[ndose] - dose_values[k]) *
        (regular_probability[ndose] - regular_probability[k])
  } else {
    k <- max(which(dose_values <= effective_dose))
    probability <- regular_probability[k] +
      (effective_dose - dose_values[k]) /
        (dose_values[k + 1L] - dose_values[k]) *
        (regular_probability[k + 1L] - regular_probability[k])
  }
  list(
    state = effective_dose,
    probability = min(1, max(0, probability)),
    effective_dose = effective_dose
  )
}

aide_update_endpoint_status <- function(admin, t_now, endpoint = c("toxicity", "efficacy")) {
  endpoint <- match.arg(endpoint)
  event_time <- if (endpoint == "toxicity") admin$t_dlt else admin$t_response
  final <- if (endpoint == "toxicity") admin$dlt_final else admin$eff_final
  observed <- ifelse(final == 1L & event_time <= t_now, 1L,
                     ifelse(admin$assessment_end <= t_now, 0L, NA_integer_))
  weight <- ifelse(is.na(observed), pmin(1, pmax(0, (t_now - admin$t_start) / (admin$assessment_end - admin$t_start))), 1)
  list(y = observed, weight = weight, ascertained = !is.na(observed))
}

aide_phase12_tite_toxicity_data <- function(admin, t_now, T_assess) {
  if (!nrow(admin)) return(data.frame())
  s <- aide_update_endpoint_status(admin, t_now, "toxicity")
  previous_admin_id <- if ("previous_admin_id" %in% names(admin)) {
    admin$previous_admin_id
  } else {
    rep.int(0L, nrow(admin))
  }
  data.frame(admin_id = admin$admin_id, patient_id = admin$patient_id, cycle = admin$cycle, t_start = admin$t_start,
             dose = admin$dose, previous_dose = admin$previous_dose,
             previous_admin_id = previous_admin_id,
             ipde = as.integer(admin$assignment_type == "retreat"), y = s$y, weight = s$weight,
             ascertained = s$ascertained)
}

aide_phase12_tite_efficacy_data <- function(admin, t_now, T_assess) {
  if (!nrow(admin)) return(data.frame())
  s <- aide_update_endpoint_status(admin, t_now, "efficacy")
  previous_admin_id <- if ("previous_admin_id" %in% names(admin)) {
    admin$previous_admin_id
  } else {
    rep.int(0L, nrow(admin))
  }
  data.frame(admin_id = admin$admin_id, patient_id = admin$patient_id, cycle = admin$cycle, t_start = admin$t_start,
             dose = admin$dose, previous_dose = admin$previous_dose,
             previous_admin_id = previous_admin_id,
             ipde = as.integer(admin$assignment_type == "retreat"), y = s$y, weight = s$weight,
             ascertained = s$ascertained)
}

aide_build_interim_data <- function(state, t_now, config) list(
  toxicity = aide_phase12_tite_toxicity_data(state$admin, t_now, config$time$T_assess),
  efficacy = aide_phase12_tite_efficacy_data(state$admin, t_now, config$time$T_assess)
)

aide_count_evaluated_current_dose <- function(admin, current_dose, t_now, T_assess) {
  if (!nrow(admin)) return(list(n = 0L, qualifying_admin_ids = integer(0)))
  at_dose <- admin$dose == current_dose
  qualified <- at_dose & ((admin$dlt_final == 1L & admin$t_dlt <= t_now) |
                          (admin$assessment_end <= t_now & admin$dlt_final == 0L))
  ids <- admin$admin_id[qualified]
  list(n = length(ids), qualifying_admin_ids = ids)
}
