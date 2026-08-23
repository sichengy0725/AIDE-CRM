.this_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1L])
source(file.path(dirname(dirname(dirname(normalizePath(.this_file)))), "TITE-AIDE.R"))

aide_test_multicycle_fit <- function(interim, config, ndose) {
  if (ndose != 3L) stop("This recycle-risk fixture requires three doses.")
  history <- aide_multicycle_history_data(interim)
  regular_draws <- rbind(
    c(.05, .10, .10),
    c(.10, .25, .30)
  )
  carryover_draws <- c(.50, .50)
  list(
    model = "multicycle_additive",
    p_regular_mean = colMeans(regular_draws),
    history_data = history,
    regular_draws = regular_draws,
    carryover_draws = carryover_draws,
    state_draws = aide_multicycle_state_draws(
      regular_draws, carryover_draws, history
    )
  )
}

aide_test_append_admin <- function(state, admin_id, patient_id, cycle, dose,
                                   t_start, assessment_end, previous_admin_id = 0L,
                                   previous_dose = NA_integer_, assignment_type = "new") {
  state$admin <- aide_add_row(state$admin, list(
    admin_id = as.integer(admin_id), patient_id = as.integer(patient_id),
    cohort_id = 1L, decision_id = 17L, stage = "one_stage",
    assignment_type = assignment_type, dose = as.integer(dose),
    decision_next_dose = as.integer(dose), previous_dose = previous_dose,
    previous_admin_id = as.integer(previous_admin_id), cycle = as.integer(cycle),
    t_arrival = as.numeric(t_start), t_start = as.numeric(t_start),
    t_dlt = Inf, t_response = Inf, assessment_end = as.numeric(assessment_end),
    dlt_final = 0L, eff_final = 0L, true_p_tox = 0, true_p_eff = 0,
    true_toxicity_state = 0, true_efficacy_state = 0
  ))
  state
}

aide_test_recycle_fixture <- function(config, scenario) {
  state <- aide_phase12_initialize_state(config, scenario, seed = 22L)
  state$t_now <- 10
  state$counters$decision <- 17L
  state$cohort$cohort_id <- 9L
  state$cohort$next_dose <- 2L
  state$cohort$decision_id <- 17L
  state$cohort$capacity <- 1L
  state$cohort$filled <- 0L
  state$cohort$allocation_doses <- 2L
  state$cohort$allocation_probabilities <- setNames(1, "D2")

  ## Administrations 1:2 were available when the cohort opened.  The later
  ## cycle-3 administration is the source of the newly waiting recycle.
  state <- aide_test_append_admin(state, 1L, 1L, 1L, 1L, 0, 2)
  state <- aide_test_append_admin(state, 2L, 1L, 2L, 2L, 2, 4,
                                  previous_admin_id = 1L, previous_dose = 1L,
                                  assignment_type = "retreat")
  state <- aide_test_append_admin(state, 3L, 1L, 3L, 3L, 4, 6,
                                  previous_admin_id = 2L, previous_dose = 2L,
                                  assignment_type = "retreat")
  ## This pending administration confirms that the refreshed fit receives
  ## normal TITE partial-follow-up information as well.
  state <- aide_test_append_admin(state, 4L, 2L, 1L, 1L, 9, 19)
  state$counters$admin <- 4L

  stale_interim <- aide_phase12_tite_toxicity_data(
    state$admin[seq_len(2L), , drop = FALSE], state$t_now, config$time$T_assess
  )
  state$cohort$risk_context <- list(
    toxicity_fit = aide_test_multicycle_fit(
      stale_interim, config, ndose = length(state$eliminated)
    ),
    toxicity_fit_time = 4,
    toxicity_fit_n_admin = 2L
  )
  state <- aide_queue_add(state, state$t_now, 1L, "retreat", 1L, 3L)
  state
}

run_recycle_risk_refresh_tests <- function() {
  cfg <- aide_phase12_config(
    Nmax = 12L, cohort_size = 1L, n_eval = 1L, cycle_max = 4L,
    toxicity = list(model = "multicycle_additive"),
    recycle = list(
      apply_individual_toxicity_risk = TRUE,
      toxicity_ipde_overdose_cutoff = .75
    )
  )
  sce <- aide_phase12_scenario(
    c(.05, .15, .30), c(.10, .30, .45),
    list(alpha_true = .10), list(alpha_true = .10)
  )

  original_fit <- aide_fit_toxicity
  fit_calls <- 0L
  captured_interim <- NULL
  assign("aide_fit_toxicity", function(interim, config, ndose) {
    fit_calls <<- fit_calls + 1L
    captured_interim <<- interim
    aide_test_multicycle_fit(interim, config, ndose)
  }, envir = .GlobalEnv)
  on.exit(assign("aide_fit_toxicity", original_fit, envir = .GlobalEnv), add = TRUE)

  state <- aide_test_recycle_fixture(cfg, sce)
  source_admin_id <- state$queue$source_admin_id[1L]
  stopifnot(!(source_admin_id %in%
    state$cohort$risk_context$toxicity_fit$history_data$admin_id))
  old_next_dose <- state$cohort$next_dose
  old_decision_id <- state$cohort$decision_id
  old_decision_count <- state$counters$decision
  old_decision_log_rows <- nrow(state$logs$decision_log)

  state <- aide_refresh_open_cohort_recycle_risk(state, cfg)
  refreshed_fit <- state$cohort$risk_context$toxicity_fit
  stopifnot(
    fit_calls == 1L,
    source_admin_id %in% refreshed_fit$history_data$admin_id,
    all(state$queue$source_admin_id[
      state$queue$type == "retreat" & state$queue$status == "waiting"
    ] %in% refreshed_fit$history_data$admin_id),
    identical(state$cohort$next_dose, old_next_dose),
    identical(state$cohort$decision_id, old_decision_id),
    identical(state$counters$decision, old_decision_count),
    nrow(state$logs$decision_log) == old_decision_log_rows
  )
  pending_row <- match(4L, captured_interim$admin_id)
  stopifnot(is.na(captured_interim$y[pending_row]),
            captured_interim$weight[pending_row] > 0,
            captured_interim$weight[pending_row] < 1)

  ## A second screen at the same event state reuses the just-refreshed fit.
  state <- aide_refresh_open_cohort_recycle_risk(state, cfg)
  stopifnot(fit_calls == 1L)
  screen <- aide_retreat_individual_risk_screen(state, 1L, cfg)
  stopifnot(screen$allowed, identical(screen$eligible_doses, old_next_dose))

  ## The posterior recursion uses the full state from all prior cycles, not
  ## merely the regular probability at the immediate prior dose.
  next_draws <- aide_multicycle_next_probability_draws(
    refreshed_fit, source_admin_id, old_next_dose
  )
  regular <- refreshed_fit$regular_draws
  alpha <- refreshed_fit$carryover_draws
  expected_full_state <- regular[, 2L] + alpha * (
    regular[, 3L] + alpha * (regular[, 2L] + alpha * regular[, 1L])
  )
  immediate_previous_dose_only <- regular[, 2L] + alpha * regular[, 3L]
  stopifnot(isTRUE(all.equal(next_draws, expected_full_state)),
            any(abs(next_draws - immediate_previous_dose_only) > 1e-8))

  ## The recycle threshold is strict: an overdose probability equal to the
  ## cutoff is not eligible.  The fixture deliberately gives probability .5.
  strict_cfg <- cfg
  strict_cfg$recycle$toxicity_ipde_overdose_cutoff <- .5
  strict_screen <- aide_retreat_individual_risk_screen(state, 1L, strict_cfg)
  stopifnot(!strict_screen$allowed,
            identical(strict_screen$reason, "blocked_individual_toxicity_risk"),
            identical(as.numeric(strict_screen$toxicity_overdose_probability), .5))

  ## A non-multicycle toxicity model does not refresh an individual risk fit.
  nonmulticycle_cfg <- aide_phase12_config(
    Nmax = 12L, cohort_size = 1L, n_eval = 1L, cycle_max = 4L,
    toxicity = list(model = "discount_r"),
    recycle = list(
      apply_individual_toxicity_risk = TRUE,
      toxicity_ipde_overdose_cutoff = .75
    )
  )
  nonmulticycle_state <- aide_test_recycle_fixture(nonmulticycle_cfg, sce)
  nonmulticycle_risk <- nonmulticycle_state$cohort$risk_context
  calls_before_nonmulticycle <- fit_calls
  nonmulticycle_state <- aide_refresh_open_cohort_recycle_risk(
    nonmulticycle_state, nonmulticycle_cfg
  )
  stopifnot(
    fit_calls == calls_before_nonmulticycle,
    identical(nonmulticycle_state$cohort$risk_context, nonmulticycle_risk)
  )

  ## The fill path refreshes once before selecting the waiting retreat and
  ## assigns it only to the already frozen cohort dose.
  fill_state <- aide_test_recycle_fixture(cfg, sce)
  fit_calls <- 0L
  fill_state <- aide_fill_open_cohort(fill_state, cfg, sce)
  assigned <- fill_state$admin[fill_state$admin$admin_id == 5L, , drop = FALSE]
  stopifnot(
    fit_calls == 1L,
    nrow(assigned) == 1L,
    identical(assigned$dose, old_next_dose),
    identical(fill_state$cohort$next_dose, old_next_dose),
    identical(fill_state$cohort$decision_id, old_decision_id),
    identical(fill_state$counters$decision, old_decision_count)
  )
}

run_recycle_risk_refresh_tests()
message("recycle-risk refresh tests passed")
