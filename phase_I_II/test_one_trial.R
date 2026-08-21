## ============================================================
## One-scenario IPDE test for the non-TITE AIDE Phase I/II design
##
## Run from the project directory with:
##   Rscript phase_I_II/test_one_trial.R
##
## This fixed-seed scenario exercises flexible IPDE recycling, including a
## patient recycled from dose 3 to dose 2. It reports every allocation and
## every IPDE eligibility decision.
## ============================================================

source("phase_I_II/AIDE_phase_I_II.R")

## Change this to "one_stage" to run the one-stage design.
allocation_to_run <- "two_stage"

p_true <- c(0.05, 0.10, 0.20, 0.30, 0.50)
e_true <- c(0.20, 0.40, 0.50, 0.5, 0.5)
toxicity_ipde_alpha <- 0
efficacy_ipde_alpha <- 0

fit <- simulate_AIDE_phase_I_II(
  p_true = p_true,
  e_true = e_true,
  toxicity_ipde_alpha = toxicity_ipde_alpha,
  efficacy_ipde_alpha = efficacy_ipde_alpha,
  allocation = allocation_to_run,
  model = "CRM",
  crm_r_model = "previous_dose",
  crm_skeleton = c(0.1,0.2,0.3,0.4,0.5),
  target = 0.30,
  utility_type = 3L,
  utility_scores = c(u00 = 0, u01 = 40, u10 = 60, u11 = 100),
  N_s1 = 6L,
  N_s2 = 30L,
  Nmax = 30L,
  C = 3L,
  cycle_max = 2L,
  arrival_rate = 1/28,
  ## Continuous enrollment gives waiting new patients priority, then fills
  ## the remaining cohort positions with eligible recycled patients.
  enrollment_scheme = "continuous",
  apply_ipde_toxicity_rule = TRUE,
  apply_ipde_efficacy_rule = FALSE,
  efficacy_threshold = 0.05,
  T_assess = 28,
  seed = 113L
)

cat("\n========================================\n")
cat("AIDE Phase I/II flexible-IPDE scenario\n")
cat("========================================\n")
cat("Final MTD:", fit$final$MTD, "\n")
cat("Final OBD:", fit$final$OBD, "\n")
cat("IPDE toxicity generation: min(1, p[d2] + alpha * p[d1]), alpha =",
    toxicity_ipde_alpha, "\n")
cat("IPDE efficacy generation: min(1, e[d2] + alpha * e[d1]), alpha =",
    efficacy_ipde_alpha, "\n\n")

cat("Dose-level true probabilities\n")
print(data.frame(
  dose = seq_along(p_true),
  regular_toxicity = p_true,
  regular_efficacy = e_true
), row.names = FALSE)

cat("\nAll cohort allocation decisions\n")
print(fit$decision_log, row.names = FALSE)

cat("\nAll administrations (including recycled patients)\n")
print(
  fit$admin[, c(
    "row_id", "cohort", "stage", "id", "ncycle", "type", "dose", "y", "eff",
    "true_toxicity_probability", "true_efficacy_probability"
  )],
  row.names = FALSE
)

cat("\nAll IPDE eligibility decisions and selected recycled patients\n")
if (nrow(fit$recycling_decision_log) == 0L) {
  cat("No patients were considered for recycling.\n")
} else {
  print(fit$recycling_decision_log, row.names = FALSE)
}

cat("\nWaiting new-patient queue at trial end\n")
if (nrow(fit$waiting_queue) == 0L) {
  cat("No arrived new patients remain unassigned.\n")
} else {
  print(fit$waiting_queue, row.names = FALSE)
}
