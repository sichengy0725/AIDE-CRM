## ============================================================
## One-scenario IPDE test for the non-TITE AIDE Phase I/II design
##
## Run from the project directory with:
##   Rscript run_AIDE_phase_I_II_example.R
##
## This fixed-seed scenario exercises flexible IPDE recycling, including a
## patient recycled from dose 3 to dose 2. It reports every allocation and
## every IPDE eligibility decision.
## ============================================================

source("AIDE_phase_I_II.R")

p_true <- c(0.05, 0.10, 0.20, 0.30, 0.50)
e_true <- c(0.20, 0.40, 0.50, 0.5, 0.5)
p_ipde <- p_true
toxicity_ipde_alpha <- 0.3
e_ipde <- e_true
efficacy_ipde_alpha <- 0.3

fit <- simulate_AIDE_phase_I_II(
  p_true = p_true,
  e_true = e_true,
  p_ipde = p_ipde,
  toxicity_ipde_alpha = toxicity_ipde_alpha,
  e_ipde = e_ipde,
  efficacy_ipde_alpha = efficacy_ipde_alpha,
  allocation = "two_stage",
  model = "CRM",
  crm_r_model = "random",
  crm_skeleton = c(0.1,0.2,0.3,0.4,0.5),
  target = 0.30,
  N_s1 = 30L,
  N_s2 = 60L,
  Nmax = 60L,
  C = 3L,
  cycle_max = 2L,
  ipde_design = 2L,
  flexible_ipde = TRUE,
  ## IPDE-first recruitment fills cohort positions with eligible recycled
  ## patients before recruiting the remaining number of new patients.
  enrollment_scheme = "ipde_first",
  efficacy_threshold = 0.05,
  T_assess = 1,
  seed = 113L
)

cat("\n========================================\n")
cat("AIDE Phase I/II flexible-IPDE scenario\n")
cat("========================================\n")
cat("Final MTD:", fit$final$MTD, "\n")
cat("Final OBD:", fit$final$OBD, "\n")
cat("IPDE toxicity generation: min(1, p2 + alpha * p1), alpha =",
    toxicity_ipde_alpha, "\n")
cat("IPDE efficacy generation: min(1, p2 + alpha * p1), alpha =",
    efficacy_ipde_alpha, "\n\n")

cat("Dose-level true probabilities\n")
print(data.frame(
  dose = seq_along(p_true),
  regular_toxicity = p_true,
  recycled_toxicity_base_p2 = p_ipde,
  regular_efficacy = e_true,
  recycled_efficacy_base_p2 = e_ipde
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
