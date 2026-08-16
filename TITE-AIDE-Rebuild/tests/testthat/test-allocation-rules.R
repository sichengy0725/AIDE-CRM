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
    stage = rep(stage, length(doses)),
    assignment_type = rep("new", length(doses)),
    dose = doses, decision_next_dose = doses,
    previous_dose = rep(NA_integer_, length(doses)),
    cycle = rep(1L, length(doses)),
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
    toxicity = list(
      p_regular_mean = p_toxicity,
      eliminated = rep(FALSE, length(p_toxicity))
    ),
    efficacy = list(
      p_regular_mean = p_efficacy,
      futile = rep(FALSE, length(p_efficacy))
    )
  )
}

toxicity_stay <- list(action = "stay", recommended_dose = 3L, stop = FALSE)
two_stage_config <- aide_phase12_config(
  allocation = "two_stage", Nmax = 18L, cohort_size = 3L, s1_Max = 3L,
  monitoring = list(stage2_allocation = "highest_utility")
)
fits <- make_fits(c(.1, .2, .3, .4, .5), c(.9, .4, .2, .1, .1))

## Highest-utility Stage II selects a response-producing dose even when it is
## not the ordinary utility leader; if there are no observed responses, it
## falls back to the current MTD.
case <- make_allocation_state(two_stage_config, "stage2", 3L, 1:3, 2L)
admissible <- aide_build_admissible_set(
  case$state, fits$toxicity, fits$efficacy, two_stage_config, case$scenario
)
stopifnot(aide_two_stage_decision(
  case$state, toxicity_stay, admissible, fits$toxicity, fits$efficacy,
  two_stage_config
)$next_dose == 2L)

case <- make_allocation_state(two_stage_config, "stage2", 3L, 1:3)
admissible <- aide_build_admissible_set(
  case$state, fits$toxicity, fits$efficacy, two_stage_config, case$scenario
)
stopifnot(aide_two_stage_decision(
  case$state, toxicity_stay, admissible, fits$toxicity, fits$efficacy,
  two_stage_config
)$next_dose == 3L)

## Top-two randomization pairs the only responding ordinary leader with the
## MTD; with no responses it uses the MTD and highest admissible lower dose.
two_stage_config$monitoring$stage2_allocation <- "top2_randomized"
case <- make_allocation_state(two_stage_config, "stage2", 3L, 1:3, 1L)
admissible <- aide_build_admissible_set(
  case$state, fits$toxicity, fits$efficacy, two_stage_config, case$scenario
)
decision <- aide_two_stage_decision(
  case$state, toxicity_stay, admissible, fits$toxicity, fits$efficacy,
  two_stage_config
)
stopifnot(identical(sort(names(decision$allocation_probabilities)), c("D1", "D3")))

case <- make_allocation_state(two_stage_config, "stage2", 3L, 1:3)
admissible <- aide_build_admissible_set(
  case$state, fits$toxicity, fits$efficacy, two_stage_config, case$scenario
)
decision <- aide_two_stage_decision(
  case$state, toxicity_stay, admissible, fits$toxicity, fits$efficacy,
  two_stage_config
)
stopifnot(identical(sort(names(decision$allocation_probabilities)), c("D2", "D3")))

one_stage_config <- aide_phase12_config(
  allocation = "one_stage", Nmax = 18L, cohort_size = 3L
)
one_stage_fits <- make_fits(c(.1, .2, .25, .3, .45), c(.9, .8, .7, .1, .1))

## With a lower utility target, one-stage allocation selects the highest
## response-qualified destination; absent a response, it stays at the current
## admissible dose.
case <- make_allocation_state(one_stage_config, "one_stage", 4L, 1:4, 3L)
admissible <- aide_build_admissible_set(
  case$state, one_stage_fits$toxicity, one_stage_fits$efficacy,
  one_stage_config, case$scenario
)
stopifnot(aide_one_stage_decision(
  case$state, toxicity_stay, admissible, one_stage_fits$toxicity,
  one_stage_fits$efficacy, one_stage_config
)$next_dose == 3L)

case <- make_allocation_state(one_stage_config, "one_stage", 4L, 1:4)
admissible <- aide_build_admissible_set(
  case$state, one_stage_fits$toxicity, one_stage_fits$efficacy,
  one_stage_config, case$scenario
)
stopifnot(aide_one_stage_decision(
  case$state, toxicity_stay, admissible, one_stage_fits$toxicity,
  one_stage_fits$efficacy, one_stage_config
)$next_dose == 4L)

## Upward movement remains no-skipping when an untried dose is the utility
## target.
upward_fits <- make_fits(c(.1, .15, .2, .25, .3), c(.1, .2, .3, .4, .9))
case <- make_allocation_state(one_stage_config, "one_stage", 3L, 1:3)
admissible <- aide_build_admissible_set(
  case$state, upward_fits$toxicity, upward_fits$efficacy,
  one_stage_config, case$scenario
)
stopifnot(aide_one_stage_decision(
  case$state, toxicity_stay, admissible, upward_fits$toxicity,
  upward_fits$efficacy, one_stage_config
)$next_dose == 4L)

message("allocation-rule tests passed")
