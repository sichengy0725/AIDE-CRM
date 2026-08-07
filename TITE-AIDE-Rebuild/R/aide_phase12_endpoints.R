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
  data.frame(admin_id = admin$admin_id, patient_id = admin$patient_id, cycle = admin$cycle, t_start = admin$t_start,
             dose = admin$dose, previous_dose = admin$previous_dose,
             ipde = as.integer(admin$assignment_type == "retreat"), y = s$y, weight = s$weight,
             ascertained = s$ascertained)
}

aide_phase12_tite_efficacy_data <- function(admin, t_now, T_assess) {
  if (!nrow(admin)) return(data.frame())
  s <- aide_update_endpoint_status(admin, t_now, "efficacy")
  data.frame(admin_id = admin$admin_id, patient_id = admin$patient_id, cycle = admin$cycle, t_start = admin$t_start,
             dose = admin$dose, previous_dose = admin$previous_dose,
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
