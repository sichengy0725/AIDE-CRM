## Regression checks for the non-TITE AIDE Phase I/II designs.
## Run with: Rscript test_AIDE_phase_I_II.R

source("AIDE_phase_I_II.R")

## Utility calculations and the exposed PDF-design controls do not require
## JAGS, so verify them even in lightweight R installations.
utility_e <- c(0.40, 0.70)
utility_t <- c(0.10, 0.30)
stopifnot(identical(
  aide_phase12_compute_utility(utility_e, utility_t, utility_type = 1L),
  utility_e
))
stopifnot(isTRUE(all.equal(
  aide_phase12_compute_utility(
    utility_e, utility_t, utility_type = 2L, lambda_T = 2
  ),
  c(0.20, 0.10)
)))
stopifnot(isTRUE(all.equal(
  aide_phase12_compute_utility(
    utility_e, utility_t, utility_type = 3L,
    utility_scores = c(u00 = 0, u01 = 2, u10 = -1, u11 = 3)
  ),
  c(0.78, 1.52)
)))
phase12_formals <- formals(simulate_AIDE_phase_I_II)
stopifnot(
  identical(phase12_formals$m_U, 6L),
  all(c("continuous", "ipde_first") %in%
      eval(phase12_formals$enrollment_scheme)),
  all(c("previous_dose") %in% eval(phase12_formals$crm_r_model)),
  all(c(
    "shared_carryover", "dose_specific_carryover",
    "previous_dose_additive", "dose_specific_previous_dose_additive"
  ) %in% eval(phase12_formals$efficacy_model))
)
stopifnot(
  all(c("apply_ipde_toxicity_rule", "ipde_toxicity_cutoff") %in%
        names(phase12_formals)),
  !any(c(
    "apply_random_crm_recycle_toxicity_rule",
    "random_crm_recycle_toxicity_cutoff",
    "apply_random_crm_recycle_efficacy_rule",
    "random_crm_recycle_efficacy_delta",
    "ipde_design",
    "flexible_ipde"
  ) %in% names(phase12_formals))
)

## A--F: arithmetic checks for the multicycle patient-level toxicity gate.
## The posterior has one deterministic draw to make the target quantities
## directly auditable.
multicycle_post <- cbind(
  "p[1]" = 0.10,
  "p[2]" = 0.20,
  "p[3]" = 0.30,
  alpha = 0.50
)
three_cycle_history <- crm_multicycle_history_data(
  data.frame(
    id = rep(1L, 3L), ncycle = 1:3, dose = 1:3,
    y = c(0L, 0L, 0L), type = c("new", "retreat", "retreat")
  ),
  ndose = 3L
)
three_cycle_states <- crm_multicycle_toxicity_posterior_states(
  multicycle_post, three_cycle_history
)
## A. 0.10, 0.25, and 0.425 are the full uncapped states.
stopifnot(isTRUE(all.equal(
  as.numeric(three_cycle_states$state_draws[1L, 2:4]),
  c(0.10, 0.25, 0.425)
)))

## B. No intermediate state cap: 0.80 + 0.50 * 0.80 = 1.20.
uncapped_post <- cbind("p[1]" = 0.80, alpha = 0.50)
uncapped_history <- crm_multicycle_history_data(
  data.frame(
    id = c(1L, 1L), ncycle = c(1L, 2L), dose = c(1L, 1L),
    y = c(0L, 0L), type = c("new", "retreat")
  ),
  ndose = 1L
)
uncapped_state <- crm_multicycle_toxicity_posterior_states(
  uncapped_post, uncapped_history
)
uncapped_gate <- crm_multicycle_recycle_toxicity_gate(
  uncapped_post, uncapped_history, 1L, 1L, phi = 0.90, cutoff = 0.95
)
stopifnot(
  isTRUE(all.equal(uncapped_state$state_draws[1L, 3L], 1.20)),
  isTRUE(all.equal(uncapped_gate$q_next_draws, 1.00))
)

## C. A D1 -> D2 history and proposed D3 administration use
## p[D3] + alpha * S[D2] = 0.425.
two_cycle_history <- crm_multicycle_history_data(
  data.frame(
    id = c(1L, 1L), ncycle = c(1L, 2L), dose = c(1L, 2L),
    y = c(0L, 0L), type = c("new", "retreat")
  ),
  ndose = 3L
)
two_cycle_gate <- crm_multicycle_recycle_toxicity_gate(
  fit = multicycle_post, history = two_cycle_history, patient_id = 1L,
  next_dose = 3L, phi = 0.50, cutoff = 0.95
)
stopifnot(
  isTRUE(all.equal(two_cycle_gate$theta_posterior_mean, 0.425)),
  isTRUE(all.equal(two_cycle_gate$q_next_draws, 0.425)),
  isTRUE(all.equal(two_cycle_gate$q_next_draws, 0.30 + 0.50 * 0.25))
)

## E. With only a D1 history, the gate reduces to the two-cycle formula
## p[D2] + alpha * p[D1] = 0.25.
one_cycle_history <- crm_multicycle_history_data(
  data.frame(id = 1L, ncycle = 1L, dose = 1L, y = 0L, type = "new"),
  ndose = 3L
)
one_cycle_gate <- crm_multicycle_recycle_toxicity_gate(
  multicycle_post, one_cycle_history, 1L, 2L, phi = 0.50, cutoff = 0.95
)
stopifnot(isTRUE(all.equal(one_cycle_gate$q_next_draws, 0.20 + 0.50 * 0.10)))

## D. Patients at current D2 can have different full histories and gates.
two_patient_history <- crm_multicycle_history_data(
  data.frame(
    id = c(1L, 1L, 2L, 2L),
    ncycle = c(1L, 2L, 1L, 2L),
    dose = c(1L, 2L, 3L, 2L),
    y = c(0L, 0L, 0L, 0L),
    type = c("new", "retreat", "new", "retreat")
  ),
  ndose = 3L
)
patient_one_gate <- crm_multicycle_recycle_toxicity_gate(
  multicycle_post, two_patient_history, 1L, 3L, phi = 0.45, cutoff = 0.95
)
patient_two_gate <- crm_multicycle_recycle_toxicity_gate(
  multicycle_post, two_patient_history, 2L, 3L, phi = 0.45, cutoff = 0.95
)
stopifnot(
  patient_one_gate$current_dose == 2L,
  patient_two_gate$current_dose == 2L,
  isTRUE(all.equal(patient_one_gate$theta_posterior_mean, 0.425)),
  isTRUE(all.equal(patient_two_gate$theta_posterior_mean, 0.475)),
  patient_one_gate$probability_over_phi == 0,
  patient_two_gate$probability_over_phi == 1
)

## F. alpha = 0 removes all history from the proposed-cycle probability.
zero_alpha_post <- multicycle_post
zero_alpha_post[, "alpha"] <- 0
zero_alpha_gate <- crm_multicycle_recycle_toxicity_gate(
  zero_alpha_post, two_patient_history, 2L, 3L, phi = 0.50, cutoff = 0.95
)
stopifnot(isTRUE(all.equal(zero_alpha_gate$q_next_draws, 0.30)))

## G. The individual toxicity gate is never active for non-multicycle CRM,
## BOIN, or disabled multicycle configurations.
stopifnot(
  aide_phase12_multicycle_toxicity_gate_is_active("CRM", "multicycle_additive", TRUE),
  !aide_phase12_multicycle_toxicity_gate_is_active("CRM", "random", TRUE),
  !aide_phase12_multicycle_toxicity_gate_is_active("BOIN", "multicycle_additive", TRUE),
  !aide_phase12_multicycle_toxicity_gate_is_active("CRM", "multicycle_additive", FALSE)
)
cat("Multicycle toxicity gate checks A-G passed.\n")

## The common one-stage rule first uses the lowest response-producing dose as
## the candidate-interval floor, and otherwise uses the MTD as its target.
## No-skipping is evaluated from the current dose rather than the highest dose
## ever tried.
response_rule <- aide_phase12_one_stage_allocation(
  current_dose = 4L,
  mtd = 4L,
  admissible_doses = 1:4,
  utility = c(.95, .40, .80, .70),
  response_observed = c(FALSE, TRUE, FALSE, FALSE),
  tried_doses = 1:4
)
stopifnot(response_rule$target == 3L, response_rule$dose == 3L,
          response_rule$response_floor == 2L)
no_response_rule <- aide_phase12_one_stage_allocation(
  current_dose = 2L,
  mtd = 4L,
  admissible_doses = 1:4,
  utility = c(.20, .40, .90, .80),
  response_observed = rep(FALSE, 4L),
  tried_doses = c(1L, 2L, 4L)
)
stopifnot(no_response_rule$target == 4L, no_response_rule$dose == 3L,
          identical(no_response_rule$action, "no_skip_untried_dose"))

## Option B moves one level toward the interim OBD.  A standing futility
## exclusion at the adjacent dose is bypassed, but the move cannot go beyond
## the selected OBD.
one_level_rule <- aide_phase12_one_stage_allocation(
  current_dose = 2L,
  mtd = 5L,
  admissible_doses = c(1L, 2L, 4L, 5L),
  utility = c(.1, .2, .3, .9, .8),
  response_observed = rep(FALSE, 5L),
  tried_doses = 1:5,
  allocation_mode = "one_level_toward_obd"
)
stopifnot(
  one_level_rule$target == 5L,
  one_level_rule$dose == 4L,
  identical(one_level_rule$action, "one_level_escalate_bypass_excluded")
)
one_level_down_rule <- aide_phase12_one_stage_allocation(
  current_dose = 5L,
  mtd = 5L,
  admissible_doses = c(1L, 2L, 3L, 5L),
  utility = c(.1, .9, .4, .2, .3),
  response_observed = c(TRUE, FALSE, FALSE, FALSE, FALSE),
  tried_doses = 1:5,
  allocation_mode = "one_level_toward_obd"
)
stopifnot(
  one_level_down_rule$target == 2L,
  one_level_down_rule$dose == 3L,
  identical(one_level_down_rule$action, "one_level_deescalate_bypass_excluded")
)

## Alternative true generators in the August 31 specification.
truth_probability <- c(.10, .20, .30, .40, .50)
dose_specific_truth <- aide_phase12_true_generation_settings(
  model = "dose_specific_geometric", ndose = 5L,
  dose_specific_alpha = c(.2, .3, .4, .5, .6)
)
dose_specific_cycle3 <- aide_phase12_true_endpoint_probability(
  settings = dose_specific_truth,
  regular_probability = truth_probability,
  current_dose = 3L,
  history_doses = c(1L, 2L)
)
stopifnot(isTRUE(all.equal(
  dose_specific_cycle3$probability,
  .30 + .3 * .20 + .2^2 * .10
)))

logistic_truth <- aide_phase12_true_generation_settings(
  model = "shared_patient_logistic", ndose = 5L
)
logistic_probability <- aide_phase12_true_endpoint_probability(
  settings = logistic_truth,
  regular_probability = truth_probability,
  current_dose = 2L,
  patient_random_effect = .5
)
stopifnot(isTRUE(all.equal(
  logistic_probability$probability, stats::plogis(stats::qlogis(.20) + .5)
)))

effective_truth <- aide_phase12_true_generation_settings(
  model = "effective_dose_geometric", ndose = 5L,
  effective_dose_values = c(15, 20, 30, 35, 45),
  effective_dose_alpha = .5
)
effective_cycle3 <- aide_phase12_true_endpoint_probability(
  settings = effective_truth,
  regular_probability = truth_probability,
  current_dose = 3L,
  history_doses = c(1L, 2L)
)
stopifnot(
  isTRUE(all.equal(effective_cycle3$effective_dose, 43.75)),
  isTRUE(all.equal(effective_cycle3$probability, .4875))
)
carryover_map <- lapply(
  c(
    "discount_shared", "discount_dose_specific",
    "additive_shared", "additive_dose_specific"
  ),
  aide_phase12_carryover_model_map
)
stopifnot(
  identical(vapply(carryover_map, `[[`, character(1), "toxicity_model"), c(
    "discount_shared", "discount_shared", "additive_shared", "additive_shared"
  )),
  identical(vapply(carryover_map, `[[`, character(1), "crm_r_model"), c(
    "random", "random", "previous_dose", "previous_dose"
  )),
  identical(vapply(carryover_map, `[[`, character(1), "efficacy_model"), c(
    "shared_carryover", "dose_specific_carryover",
    "previous_dose_additive", "dose_specific_previous_dose_additive"
  ))
)
stopifnot(
  !aide_phase12_dose_threshold_reached(c(3L, 3L, 3L), 6L),
  aide_phase12_dose_threshold_reached(c(3L, 6L, 3L), 6L)
)
## The additive previous-dose models must pair a recycled administration with
## that same patient's immediately preceding dose, not the last global dose.
previous_dose_admin <- data.frame(
  id = c(1L, 2L, 1L, 2L),
  ncycle = c(1L, 1L, 2L, 2L),
  dose = c(1L, 2L, 3L, 1L),
  y = c(0L, 0L, 0L, 1L),
  eff = c(0L, 1L, 1L, 0L),
  type = c("new", "new", "retreat", "retreat")
)
tox_previous_dose <- crm_previous_dose_index(
  crm_prepare_dat(previous_dose_admin, ndose = 3L)
)
stopifnot(identical(tox_previous_dose, c(1L, 1L, 1L, 2L)))
eff_previous_dose <- aide_phase12_previous_dose_efficacy_data(
  previous_dose_admin, ndose = 3L
)
stopifnot(
  identical(eff_previous_dose$dose, c(1L, 3L, 2L, 1L)),
  identical(eff_previous_dose$ipde, c(0L, 1L, 0L, 1L)),
  identical(eff_previous_dose$previous_dose, c(1L, 1L, 1L, 2L))
)
stopifnot(file.exists("previous_dose_additive_CRM.bug"))
stopifnot(file.exists("previous_dose_additive_beta_binomial_efficacy.jags"))
stopifnot(file.exists("previous_dose_additive_dose_specific_beta_binomial_efficacy.jags"))
stopifnot(file.exists("beta_binomial_shared_carryover.jags"))

## Additive toxicity and efficacy models implement all intended Beta priors
## through ratios of independent, equal-rate Gamma variables. The toxicity
## runner retains its public a_r/b_r interface and sends those names to JAGS.
toxicity_model_text <- paste(readLines("previous_dose_additive_CRM.bug"), collapse = "\n")
efficacy_model_text <- paste(readLines("previous_dose_additive_beta_binomial_efficacy.jags"), collapse = "\n")
efficacy_dose_specific_additive_model_text <- paste(
  readLines("previous_dose_additive_dose_specific_beta_binomial_efficacy.jags"),
  collapse = "\n"
)
crm_helper_text <- paste(readLines("AIDE_CRM_helper_modified.R"), collapse = "\n")
stopifnot(
  !grepl("dbeta", toxicity_model_text, fixed = TRUE),
  !grepl("dbeta", efficacy_model_text, fixed = TRUE),
  grepl("g_alpha_1 ~ dgamma(a_r, 1)", toxicity_model_text, fixed = TRUE),
  grepl("g_alpha_2 ~ dgamma(b_r, 1)", toxicity_model_text, fixed = TRUE),
  grepl("g_regular_1[j] ~ dgamma(a_regular[j], 1)", efficacy_model_text, fixed = TRUE),
  grepl("g_regular_2[j] ~ dgamma(b_regular[j], 1)", efficacy_model_text, fixed = TRUE),
  grepl("g_alpha_1[j] ~ dgamma(a_alpha[j], 1)",
        efficacy_dose_specific_additive_model_text, fixed = TRUE),
  grepl("alpha[dose[i]] * p_regular[previous_dose[i]]",
        efficacy_dose_specific_additive_model_text, fixed = TRUE),
  grepl("jags_data$a_r <- a_r", crm_helper_text, fixed = TRUE),
  grepl("jags_data$b_r <- b_r", crm_helper_text, fixed = TRUE)
)

## The preceding discount-factor r models use the same Gamma-ratio form for
## both toxicity and dose-specific IPDE efficacy carryover parameters.
random_crm_model_text <- paste(readLines("random_CRM.bug"), collapse = "\n")
discount_efficacy_model_text <- paste(readLines("beta_binomial_dose_specific.jags"), collapse = "\n")
shared_discount_efficacy_model_text <- paste(
  readLines("beta_binomial_shared_carryover.jags"), collapse = "\n"
)
stopifnot(
  !grepl("dbeta", random_crm_model_text, fixed = TRUE),
  !grepl("dbeta", discount_efficacy_model_text, fixed = TRUE),
  !grepl("dbeta", shared_discount_efficacy_model_text, fixed = TRUE),
  grepl("g_r_1 ~ dgamma(a_r, 1)", random_crm_model_text, fixed = TRUE),
  grepl("g_r_2 ~ dgamma(b_r, 1)", random_crm_model_text, fixed = TRUE),
  grepl("g_regular_1[j] ~ dgamma(a_r[j], 1)", discount_efficacy_model_text, fixed = TRUE),
  grepl("g_carry_1[j] ~ dgamma(a_carry[j], 1)", discount_efficacy_model_text, fixed = TRUE),
  grepl("g_carry_2[j] ~ dgamma(b_carry[j], 1)", discount_efficacy_model_text, fixed = TRUE),
  grepl("g_carry_1 ~ dgamma(a_carry, 1)", shared_discount_efficacy_model_text, fixed = TRUE),
  grepl("r_e <- g_carry_1 / (g_carry_1 + g_carry_2)",
        shared_discount_efficacy_model_text, fixed = TRUE)
)

## The efficacy-futility rule uses the observed efficacy count at the dose:
## Pr(p_E < eta | y, n) = pbeta(eta, y + 1, n - y + 1).
beta_eff_admin <- data.frame(
  dose = rep(1L, 20L),
  eff = rep(0L, 20L)
)
beta_futility <- aide_phase12_beta_binomial_futility(
  admin = beta_eff_admin,
  ndose = 2L,
  efficacy_threshold = 0.20,
  futility_cutoff = 0.95,
  min_eff_n_for_futility = 3L
)
stopifnot(
  beta_futility$futility_eliminated[1L] == 1L,
  isTRUE(all.equal(
    beta_futility$prob_below_threshold[1L],
    pbeta(0.20, 1L, 21L)
  ))
)
beta_success_admin <- data.frame(
  dose = rep(1L, 20L),
  eff = rep(1L, 20L)
)
stopifnot(!aide_phase12_beta_binomial_futility(
  admin = beta_success_admin,
  ndose = 2L,
  efficacy_threshold = 0.20,
  futility_cutoff = 0.95,
  min_eff_n_for_futility = 3L
)$futility_eliminated[1L])

if (!requireNamespace("rjags", quietly = TRUE)) {
  cat("AIDE Phase I/II regression checks skipped: rjags/JAGS is not installed.\n")
  quit(save = "no", status = 0L)
}

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

report_trial_decisions <- function(label, fit) {
  cat("\n========================================\n")
  cat(label, "-- cohort allocation decisions\n")
  cat("========================================\n")
  if (nrow(fit$decision_log) == 0L) {
    cat("No interim allocation decisions were made.\n")
  } else {
    print(fit$decision_log, row.names = FALSE)
  }
  cat("\n", label, " -- IPDE/recycling decisions and gate quantities\n", sep = "")
  if (nrow(fit$recycling_decision_log) == 0L) {
    cat("No patient was eligible for recycling consideration.\n")
  } else {
    print(fit$recycling_decision_log, row.names = FALSE)
  }
}

## Direct beta-binomial futility calculation: 0 responses from 20 patients
## must be declared futile for a 20% efficacy threshold.
eff_admin <- data.frame(
  dose = rep(1L, 20L),
  eff = rep(0L, 20L),
  y = rep(0L, 20L)
)
futility <- aide_phase12_efficacy_summary(
  admin = eff_admin,
  ndose = 2L,
  efficacy_prior = c(1, 1),
  efficacy_threshold = 0.20,
  futility_cutoff = 0.95,
  min_eff_n_for_futility = 3L
)
assert_true(
  futility$futility_eliminated[1L] == 1L,
  "Beta-binomial futility rule did not eliminate an ineffective dose."
)

## Recycled efficacy is generated as p2 + alpha * p1, capped at one.
assert_true(
  isTRUE(all.equal(
    aide_phase12_ipde_efficacy_probability(
      e_regular = c(0.10, 0.30),
      e_ipde_base = c(0.10, 0.30),
      previous_dose = 1L,
      current_dose = 2L,
      alpha = 0.50
    ),
    0.35
  )),
  "IPDE efficacy was not generated as p2 + alpha * p1."
)
assert_true(
  isTRUE(all.equal(
    aide_phase12_ipde_toxicity_probability(
      p_regular = c(0.10, 0.30),
      p_ipde_base = c(0.10, 0.30),
      previous_dose = 1L,
      current_dose = 2L,
      alpha = 0.50
    ),
    0.35
  )),
  "IPDE toxicity was not generated as p2 + alpha * p1."
)

## The dose-specific beta-binomial efficacy gate treats the recycled efficacy
## at d2 as r_d2 + (1-r_d2)p_d2 and compares it with regular efficacy at d1.
recycle_eff_admin <- rbind(
  data.frame(dose = rep(1L, 20L), eff = rep(0L, 20L), type = "new"),
  data.frame(dose = rep(2L, 20L), eff = rep(0L, 20L), type = "new"),
  data.frame(dose = rep(2L, 20L), eff = rep(1L, 20L), type = "retreat")
)
efficacy_carryover_fit <- aide_phase12_efficacy_posterior(
  recycle_eff_admin,
  ndose = 2L,
  efficacy_prior = c(1, 1),
  efficacy_carryover_prior = c(1, 9)
)
assert_true(
  efficacy_carryover_fit$carryover_posterior_mean[2L] > 0.2 &&
    efficacy_carryover_fit$ipde_posterior_mean[2L] >
      efficacy_carryover_fit$regular_posterior_mean[2L],
  "Dose-specific efficacy carryover parameter r was not fitted from IPDE data."
)
eligible_eff_gate <- aide_phase12_dose_specific_efficacy_recycle_gate(
  recycle_eff_admin,
  current_dose = 1L,
  next_dose = 2L,
  ndose = 2L,
  efficacy_prior = c(1, 1),
  delta = 0.20
)
assert_true(
  eligible_eff_gate$allowed && eligible_eff_gate$posterior_mean_increment > 0.20,
  "Dose-specific beta-binomial efficacy gate did not allow a sufficient efficacy improvement."
)
assert_true(
  !aide_phase12_dose_specific_efficacy_recycle_gate(
    recycle_eff_admin,
    current_dose = 1L,
    next_dose = 2L,
    ndose = 2L,
    efficacy_prior = c(1, 1),
    delta = 0.95
  )$allowed,
  "Dose-specific beta-binomial efficacy gate did not enforce the user delta threshold."
)
negative_increment_gate <- aide_phase12_dose_specific_efficacy_recycle_gate(
  admin = data.frame(
    dose = rep(1L, 3L),
    eff = c(1L, 1L, 0L),
    type = rep("new", 3L)
  ),
  current_dose = 1L,
  next_dose = 2L,
  ndose = 2L,
  efficacy_prior = c(1, 1),
  delta = 0.20
)
assert_true(
  !negative_increment_gate$allowed &&
    negative_increment_gate$posterior_mean_increment < 0,
  "A negative IPDE efficacy increment was incorrectly allowed."
)
common_args <- list(
  p_true = c(0, 0, 0),
  e_true = c(0.20, 0.60, 0.40),
  N_s1 = 6L,
  N_s2 = 12L,
  Nmax = 12L,
  C = 3L,
  T_assess = 1,
  cycle_max = 1L,
  efficacy_threshold = 0.20,
  seed = 114L
)

two_stage <- do.call(
  simulate_AIDE_phase_I_II,
  c(common_args, list(allocation = "two_stage"))
)
assert_true(two_stage$final$n_admin == 12L, "Two-stage design did not reach Nmax.")
assert_true(two_stage$stage1$MTD %in% 1:3, "Two-stage design did not select a Stage-I MTD.")
assert_true(two_stage$final$OBD %in% 1:3, "Two-stage design did not select an OBD.")
assert_true(
  !two_stage$final$recycling_rules$active,
  "Random-CRM recycling rules were incorrectly applied to the BOIN design."
)
assert_true(
  isTRUE(all.equal(
    two_stage$final$utility$toxicity,
    two_stage$final$toxicity_utility_fit$phat
  )),
  "Utility did not use the selected toxicity model's estimate."
)
assert_true(
  two_stage$stage1$threshold_reached || two_stage$final$n_admin == common_args$Nmax,
  "Stage I ended before a dose reached N_s1."
)
stage1_counts <- tabulate(
  two_stage$admin$dose[two_stage$admin$stage == "stage1"],
  nbins = length(common_args$p_true)
)
stage2_counts <- tabulate(
  two_stage$admin$dose[two_stage$admin$stage == "stage2"],
  nbins = length(common_args$p_true)
)
assert_true(
  nrow(two_stage$admin[two_stage$admin$stage == "stage2", , drop = FALSE]) == 0L ||
    isTRUE(two_stage$stage1$threshold_reached),
  "Stage II began before a dose reached N_s1 with a stay decision."
)
if (isTRUE(two_stage$stage1$threshold_reached)) {
  transition_rows <- subset(
    two_stage$decision_log,
    stage == "stage1" &
      current_dose == two_stage$stage1$transition_dose &
      n_current >= common_args$N_s1 &
      toxicity_action == "stay" &
      toxicity_next_dose == current_dose &
      allocated_dose == current_dose
  )
  assert_true(
    nrow(transition_rows) > 0L,
    "Stage I transition was not preceded by a stay decision at N_s1."
  )
}
assert_true(
  all(c("t_eff", "t_tox", "t_eval") %in% names(two_stage$admin)) &&
    all(two_stage$admin$t_eval == two_stage$admin$t_start +
          pmax(two_stage$admin$t_tox, two_stage$admin$t_eff)),
  "Cohort decisions were not delayed until both gen.tite endpoint outcomes were observed."
)
assert_true(
  all(!two_stage$final$admissible | two_stage$final$tried_dose),
  "Final OBD candidates were not restricted to tried doses."
)

one_stage_args <- common_args
one_stage_args$seed <- 115L
one_stage <- do.call(
  simulate_AIDE_phase_I_II,
  c(one_stage_args, list(allocation = "one_stage"))
)
assert_true(one_stage$final$n_admin == 12L, "One-stage design did not reach Nmax.")
assert_true(one_stage$final$MTD %in% 1:3, "One-stage design did not select an MTD.")
assert_true(one_stage$final$OBD %in% 1:3, "One-stage design did not select an OBD.")
assert_true(
  all(one_stage$admin$stage == "one_stage"),
  "One-stage design recorded an incorrect administration stage."
)
assert_true(
  all(one_stage$decision_log$allocated_dose %in% seq_along(one_stage$final$tried_dose)),
  "One-stage allocation returned an invalid dose."
)

## Under continuous accrual, new patients accumulate in the waiting queue
## while cohort 1 is evaluated. The priority option determines cohort 2.
priority_args <- list(
  p_true = c(0, 0),
  e_true = c(0.20, 0.40),
  allocation = "two_stage",
  N_s1 = 6L,
  Nmax = 6L,
  C = 3L,
  cycle_max = 2L,
  arrival_rate = 10,
  T_assess = 1,
  seed = 77L
)
recycle_first_trial <- do.call(
  simulate_AIDE_phase_I_II,
  c(priority_args, list(enrollment_scheme = "ipde_first"))
)
new_first_trial <- do.call(
  simulate_AIDE_phase_I_II,
  c(priority_args, list(enrollment_scheme = "continuous"))
)
assert_true(
  all(recycle_first_trial$admin$type[recycle_first_trial$admin$cohort == 2L] == "retreat") &&
    all(new_first_trial$admin$type[new_first_trial$admin$cohort == 2L] == "new"),
  "Enrollment priority did not switch between recycle-first and new-first."
)
assert_true(
  nrow(recycle_first_trial$waiting_queue) == 0L &&
    nrow(new_first_trial$waiting_queue) == 0L,
  "IPDE-first recruitment retained a continuous new-patient waiting queue."
)
assert_true(
  nrow(recycle_first_trial$arrival_schedule) ==
    sum(recycle_first_trial$admin$type == "new") &&
    nrow(new_first_trial$arrival_schedule) == priority_args$Nmax,
  "IPDE-first recruitment pre-enrolled new patients instead of recruiting only open cohort slots."
)

## A recycling-enabled BOIN run exercises the decision log. Recycling-gate
## quantities are NA because those gates are intentionally inactive outside
## multicycle CRM; efficacy itself is nevertheless fitted with JAGS.
recycling_report_trial <- simulate_AIDE_phase_I_II(
  p_true = c(0, 0, 0),
  e_true = c(0.20, 0.40, 0.50),
  allocation = "two_stage",
  N_s1 = 6L,
  Nmax = 12L,
  C = 3L,
  cycle_max = 2L,
  T_assess = 1,
  seed = 124L
)
assert_true(
  nrow(recycling_report_trial$recycling_decision_log) > 0L &&
    any(recycling_report_trial$recycling_decision_log$selected_for_recycling),
  "Recycling-enabled test did not record its IPDE decisions."
)

oc <- get_oc_sim_AIDE_phase_I_II(
  p_true = common_args$p_true,
  e_true = common_args$e_true,
  ntrial = 2L,
  N_s1 = common_args$N_s1,
  N_s2 = common_args$N_s2,
  Nmax = common_args$Nmax,
  C = common_args$C,
  T_assess = common_args$T_assess,
  allocation = "two_stage"
)
assert_true(
  length(oc$OBD_by_trial) == 2L && length(oc$OBD_selection_percent) == 3L,
  "Operating-characteristic wrapper did not return dose-level OBD results."
)

report_trial_decisions("Two-stage test", two_stage)
report_trial_decisions("One-stage test", one_stage)
report_trial_decisions("Recycling-log test", recycling_report_trial)

cat("AIDE Phase I/II regression checks passed.\n")
