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
      eval(phase12_formals$enrollment_scheme))
)

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

## Random-CRM recycling toxicity gate: reject a proposed recycled dose when
## its posterior Pr{r + (1-r)p_d > phi} reaches the supplied cutoff.
safe_post <- cbind(
  r = rep(0.10, 20L),
  "p[1]" = rep(0.10, 20L),
  "p[2]" = rep(0.20, 20L)
)
unsafe_post <- cbind(
  r = rep(0.90, 20L),
  "p[1]" = rep(0.10, 20L),
  "p[2]" = rep(0.20, 20L)
)
assert_true(
  (safe_tox_gate <- aide_phase12_random_crm_recycle_toxicity_gate(
    safe_post, next_dose = 2L, phi = 0.30, cutoff = 0.95
  ))$allowed,
  "Random-CRM recycling toxicity gate rejected a safe proposed dose."
)
assert_true(
  !(unsafe_tox_gate <- aide_phase12_random_crm_recycle_toxicity_gate(
    unsafe_post, next_dose = 2L, phi = 0.30, cutoff = 0.95
  ))$allowed,
  "Random-CRM recycling toxicity gate did not reject an unsafe proposed dose."
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
cat("\nRandom-CRM IPDE gate calculation report\n")
print(data.frame(
  case = c("toxicity_safe", "toxicity_unsafe", "efficacy_increment"),
  toxicity_probability_over_phi = c(
    safe_tox_gate$probability_over_phi,
    unsafe_tox_gate$probability_over_phi,
    NA_real_
  ),
  toxicity_cutoff = c(safe_tox_gate$cutoff, unsafe_tox_gate$cutoff, NA_real_),
  efficacy_p_regular_current = c(
    NA_real_, NA_real_, eligible_eff_gate$p_regular_current
  ),
  efficacy_theta_ipde_next = c(
    NA_real_, NA_real_, eligible_eff_gate$theta_ipde_next
  ),
  efficacy_posterior_mean_increment = c(
    NA_real_, NA_real_, eligible_eff_gate$posterior_mean_increment
  ),
  efficacy_delta = c(NA_real_, NA_real_, eligible_eff_gate$delta),
  allowed = c(
    safe_tox_gate$allowed,
    unsafe_tox_gate$allowed,
    eligible_eff_gate$allowed
  )
), row.names = FALSE)

common_args <- list(
  p_true = c(0, 0, 0),
  e_true = c(0.20, 0.60, 0.40),
  N_s1 = 6L,
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
  all(two_stage$admin$stage[1:6] == "stage1") &&
    all(two_stage$admin$stage[7:12] == "stage2"),
  "Two-stage administrations were not partitioned at N_s1."
)
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
  all(abs(one_stage$decision_log$allocated_dose -
            one_stage$decision_log$current_dose) <= 1L),
  "One-stage allocation left the required local candidate set."
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
## random CRM; efficacy itself is nevertheless fitted with JAGS.
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

## flexible_ipde permits a patient to recycle downward, including to dose 1.
## In this fixed-seed scenario, patients treated at dose 2 are reused at dose 1.
flexible_ipde_trial <- simulate_AIDE_phase_I_II(
  p_true = c(0, 0.65, 0.65),
  e_true = c(0.20, 0.40, 0.50),
  allocation = "two_stage",
  N_s1 = 9L,
  Nmax = 12L,
  C = 3L,
  cycle_max = 2L,
  ipde_design = 2L,
  flexible_ipde = TRUE,
  toxicity_ipde_alpha = 0.10,
  efficacy_ipde_alpha = 0.50,
  T_assess = 1,
  seed = 18L
)
downward_recycle <- with(
  flexible_ipde_trial$recycling_decision_log,
  selected_for_recycling & current_dose > 1L & next_dose == 1L
)
assert_true(
  any(downward_recycle),
  "flexible_ipde did not permit a patient to recycle to dose 1."
)
assert_true(
  isTRUE(flexible_ipde_trial$final$recycling_rules$flexible_ipde),
  "The flexible IPDE option was not retained in the trial output."
)
recycled_admin <- flexible_ipde_trial$admin[
  flexible_ipde_trial$admin$type == "retreat",
  , drop = FALSE
]
expected_recycled_efficacy <- vapply(seq_len(nrow(recycled_admin)), function(i) {
  one <- recycled_admin[i, , drop = FALSE]
  previous_rows <- flexible_ipde_trial$admin[
    flexible_ipde_trial$admin$id == one$id &
      flexible_ipde_trial$admin$row_id < one$row_id,
    , drop = FALSE
  ]
  previous <- previous_rows[which.max(previous_rows$row_id), , drop = FALSE]
  min(1, c(0.20, 0.40, 0.50)[one$dose] +
        0.50 * c(0.20, 0.40, 0.50)[previous$dose])
}, numeric(1))
expected_recycled_toxicity <- vapply(seq_len(nrow(recycled_admin)), function(i) {
  one <- recycled_admin[i, , drop = FALSE]
  previous_rows <- flexible_ipde_trial$admin[
    flexible_ipde_trial$admin$id == one$id &
      flexible_ipde_trial$admin$row_id < one$row_id,
    , drop = FALSE
  ]
  previous <- previous_rows[which.max(previous_rows$row_id), , drop = FALSE]
  min(1, c(0.00, 0.65, 0.65)[one$dose] +
        0.10 * c(0.00, 0.65, 0.65)[previous$dose])
}, numeric(1))
assert_true(
  isTRUE(all.equal(
    recycled_admin$true_efficacy_probability,
    expected_recycled_efficacy
  )),
  "Simulated recycled efficacy did not use p2 + alpha * p1."
)
assert_true(
  isTRUE(all.equal(
    recycled_admin$true_toxicity_probability,
    expected_recycled_toxicity
  )),
  "Simulated recycled toxicity did not use p2 + alpha * p1."
)

oc <- get_oc_sim_AIDE_phase_I_II(
  p_true = common_args$p_true,
  e_true = common_args$e_true,
  ntrial = 2L,
  N_s1 = common_args$N_s1,
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
report_trial_decisions("Flexible-IPDE test", flexible_ipde_trial)

cat("AIDE Phase I/II regression checks passed.\n")
