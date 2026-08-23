.this_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1L])
source(file.path(dirname(dirname(dirname(normalizePath(.this_file)))), "TITE-AIDE.R"))

make_allocation_state <- function(config, stage, current_dose, doses,
                                  response_doses = integer(0)) {
  scenario <- aide_phase12_scenario(rep(.2, 5), rep(.4, 5))
  state <- aide_phase12_initialize_state(config, scenario, seed = 1L)
  state$t_now <- 10
  state$current_dose <- current_dose
  state$stage$active <- stage
  state$cohort$open <- FALSE
  state$admin <- data.frame(
    admin_id = seq_along(doses), patient_id = seq_along(doses),
    cohort_id = seq_along(doses), decision_id = seq_along(doses),
    stage = rep(stage, length(doses)), assignment_type = rep("new", length(doses)),
    dose = doses, decision_next_dose = doses,
    previous_dose = rep(NA_integer_, length(doses)), cycle = rep(1L, length(doses)),
    t_arrival = 0, t_start = 0, t_dlt = Inf,
    t_response = ifelse(doses %in% response_doses, 1, Inf),
    assessment_end = 28, dlt_final = 0L,
    eff_final = as.integer(doses %in% response_doses),
    true_p_tox = .2, true_p_eff = .4
  )
  list(state = state, scenario = scenario)
}

make_fits <- function(p_toxicity, p_efficacy) {
  list(
    toxicity = list(p_regular_mean = p_toxicity, eliminated = rep(FALSE, length(p_toxicity))),
    efficacy = list(p_regular_mean = p_efficacy, futile = rep(FALSE, length(p_efficacy)))
  )
}

toxicity_stay <- list(action = "stay", recommended_dose = 3L, stop = FALSE)
one_stage_config <- aide_phase12_config(
  allocation = "one_stage", Nmax = 18L, cohort_size = 3L
)
fits <- make_fits(c(.1, .2, .3, .35, .5), c(.2, .4, .9, .8, .1))

## Toxicity never removes doses one by one. Only an excessive posterior risk
## at dose 1 marks toxicity elimination and triggers trial termination.
stopifnot(identical(
  aide_toxicity_elimination(c(.20, .99, .98, .97, .96), cutoff = .95),
  rep(FALSE, 5L)
))
dose1_overtoxic <- aide_toxicity_elimination(
  c(.96, .99, .98, .97, .96), cutoff = .95
)
stopifnot(identical(dose1_overtoxic, c(TRUE, FALSE, FALSE, FALSE, FALSE)))

## In particular, a high posterior overdose probability at the current
## non-dose-1 level no longer forces an extra de-escalation.
current_dose_fit <- list(
  p_regular_mean = c(.05, .25, .45, .60, .75),
  prob_overtox_by_dose = c(.10, .99, .99, .99, .99)
)
current_dose_rec <- aide_toxicity_recommendation(
  current_dose = 2L, toxicity_fit = current_dose_fit,
  config = one_stage_config, eliminated = rep(FALSE, 5L)
)
stopifnot(
  identical(current_dose_rec$action, "stay"),
  identical(current_dose_rec$recommended_dose, 2L),
  !current_dose_rec$stop,
  is.na(aide_current_mtd(current_dose_fit, dose1_overtoxic, .30))
)
dose1_stop <- aide_one_stage_decision(
  list(),
  aide_toxicity_recommendation(
    1L, current_dose_fit, one_stage_config, dose1_overtoxic
  ),
  list(), list(), list(), one_stage_config
)
stopifnot(
  dose1_stop$stop_trial,
  identical(dose1_stop$stop_reason, "dose_1_overtoxicity")
)

## Interim MTD is model based over all non-eliminated doses; an untried
## closest dose can guide allocation, while final selection filters to tried.
model_based_fit <- list(p_regular_mean = c(.10, .30, .28, .55, .70))
stopifnot(identical(
  aide_current_mtd(
    model_based_fit, rep(FALSE, 5L), target = .30
  ),
  2L
))

## No OBD is always represented by the shared integer sentinel, including a
## trial that terminates before any administration.
empty_scenario <- aide_phase12_scenario(rep(.2, 5L), rep(.4, 5L))
empty_state <- aide_phase12_initialize_state(one_stage_config, empty_scenario, seed = 1L)
empty_final <- aide_final_analysis(empty_state, one_stage_config, empty_scenario)
stopifnot(identical(empty_final$OBD, aide_phase12_no_obd),
          identical(aide_phase12_no_obd, 99L))

## A response at dose 2 creates the interval 2:MTD. Utility can therefore
## select dose 3, while utility at dose 1 is outside the candidate interval.
case <- make_allocation_state(one_stage_config, "one_stage", 4L, 1:4, 2L)
admissible <- aide_build_admissible_set(
  case$state, fits$toxicity, fits$efficacy, one_stage_config, case$scenario
)
admissible$mtd <- 4L
decision <- aide_one_stage_decision(
  case$state, toxicity_stay, admissible, fits$toxicity, fits$efficacy,
  one_stage_config
)
stopifnot(decision$provisional_obd == 3L, decision$next_dose == 3L,
          decision$response_floor == 2L)

## With no observed response, the MTD is the target. The no-skipping rule
## uses the path from the current dose: a previously tried dose 4 cannot let
## allocation jump over untried dose 3.
case <- make_allocation_state(one_stage_config, "one_stage", 2L, c(1L, 2L, 4L))
admissible <- aide_build_admissible_set(
  case$state, fits$toxicity, fits$efficacy, one_stage_config, case$scenario
)
admissible$mtd <- 4L
decision <- aide_one_stage_decision(
  case$state, toxicity_stay, admissible, fits$toxicity, fits$efficacy,
  one_stage_config
)
stopifnot(decision$provisional_obd == 4L, decision$next_dose == 3L,
          decision$candidate_doses[1L] == 1L)

## After Stage I, two-stage allocation is exactly the one-stage rule. It has
## one allocation dose and cannot top-two randomize.
two_stage_config <- aide_phase12_config(
  allocation = "two_stage", Nmax = 18L, cohort_size = 3L, s1_Max = 3L
)
case <- make_allocation_state(two_stage_config, "stage2", 2L, c(1L, 2L, 4L))
admissible <- aide_build_admissible_set(
  case$state, fits$toxicity, fits$efficacy, two_stage_config, case$scenario
)
admissible$mtd <- 4L
decision <- aide_two_stage_decision(
  case$state, toxicity_stay, admissible, fits$toxicity, fits$efficacy,
  two_stage_config
)
stopifnot(decision$next_dose == 3L,
          identical(decision$allocation_doses, 3L),
          identical(decision$stage2_allocation, "one_stage"),
          !decision$individual_randomization)

## Final OBD candidates must satisfy all four restrictions: tried, observed
## response, non-futile, and no higher than the final toxicity-model MTD.
stopifnot(identical(
  aide_response_qualified_obd_candidates(
    tried_doses = c(1L, 2L, 3L, 4L),
    eliminated = c(FALSE, FALSE, FALSE, FALSE, FALSE),
    futile = c(FALSE, FALSE, TRUE, FALSE, FALSE),
    mtd = 4L,
    response_observed = c(FALSE, TRUE, TRUE, TRUE, TRUE)
  ),
  c(2L, 4L)
))

message("allocation-rule tests passed")
