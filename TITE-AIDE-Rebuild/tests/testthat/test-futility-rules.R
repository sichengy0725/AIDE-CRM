.this_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1L])
source(file.path(dirname(dirname(dirname(normalizePath(.this_file)))), "TITE-AIDE.R"))

## The futility calculation has a fixed Beta(1,1) prior and ignores every
## pending efficacy record, irrespective of its TITE follow-up weight.
config <- aide_phase12_config(
  efficacy = list(
    threshold = .20,
    futility_cutoff = .85,
    min_eff_n_for_futility = 3L
  )
)
interim <- data.frame(
  dose = c(rep(1L, 10L), rep(1L, 5L), 2L),
  y = c(rep(0L, 10L), rep(NA_integer_, 5L), 1L),
  ascertained = c(rep(TRUE, 10L), rep(FALSE, 5L), TRUE),
  weight = c(rep(1, 10L), seq(.1, .5, length.out = 5L), 1)
)
futility <- aide_beta_binomial_futility(interim, config$efficacy, ndose = 2L)
stopifnot(
  identical(futility$prior, c(1, 1)),
  identical(futility$n_fully_ascertained_by_dose, c(10L, 1L)),
  identical(futility$y_fully_ascertained_by_dose, c(0L, 1L)),
  isTRUE(all.equal(futility$prob_below_threshold[1L], stats::pbeta(.20, 1, 11))),
  isTRUE(futility$futile[1L]),
  !futility$futile[2L]
)

## Stage I stays toxicity-directed, but it must redirect when the toxicity
## recommendation is no longer efficacy-admissible.
scenario <- aide_phase12_scenario(rep(.2, 3L), rep(.4, 3L))
stage1_config <- aide_phase12_config(
  allocation = "two_stage", Nmax = 12L, cohort_size = 3L, s1_Max = 12L
)
state <- aide_phase12_initialize_state(stage1_config, scenario, seed = 1L)
state$cohort$open <- FALSE
state$current_dose <- 3L
state$admin <- data.frame(
  admin_id = 1:3, patient_id = 1:3, cohort_id = 1:3, decision_id = 1:3,
  stage = rep("stage1", 3L), assignment_type = rep("new", 3L),
  dose = 1:3, decision_next_dose = 1:3, previous_dose = rep(NA_integer_, 3L),
  cycle = rep(1L, 3L), t_arrival = 0, t_start = 0, t_dlt = Inf,
  t_response = Inf, assessment_end = 28, dlt_final = 0L, eff_final = 0L,
  true_p_tox = .2, true_p_eff = .4
)
state$futile[3L] <- TRUE
toxicity_fit <- list(p_regular_mean = c(.1, .2, .3), eliminated = rep(FALSE, 3L))
efficacy_fit <- list(p_regular_mean = c(.3, .3, .3), futile = state$futile)
admissible <- aide_build_admissible_set(
  state, toxicity_fit, efficacy_fit, stage1_config, scenario
)
rm(scenario)
decision <- aide_two_stage_decision(
  state,
  list(action = "stay", recommended_dose = 3L, stop = FALSE),
  admissible, toxicity_fit, efficacy_fit, stage1_config
)
stopifnot(!decision$stop_trial, decision$stage == "stage1", decision$next_dose == 2L)

message("futility-rule tests passed")
