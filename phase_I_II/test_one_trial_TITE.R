## ============================================================
## One-scenario JAGS TITE AIDE Phase I/II example
##
## Run from the project directory with:
##   Rscript phase_I_II/test_one_trial_TITE.R
## ============================================================

source("phase_I_II/AIDE_phase_I_II_TITE.R")

allocation_to_run <- "two_stage"
p_true <- c(.05, .10, .20, .30, .50)
e_true <- c(.20, .40, .50, .50, .50)

fit <- simulate_AIDE_phase_I_II_TITE(
  p_true = p_true,
  e_true = e_true,
  toxicity_ipde_alpha = 0,
  efficacy_ipde_alpha = 0,
  allocation = allocation_to_run,
  Nmax = 30L,
  N_s1 = 6L,
  N_s2 = 30L,
  C = 3L,
  cycle_max = 1L,
  n_eval = 3L,
  m_U = 6L,
  enrollment_scheme = "continuous",
  arrival_rate = 1 / 28,
  T_assess = 28,
  dlt_dist = "uniform",
  efficacy_dist = "uniform",
  carryover_model = "additive_shared",
  crm_skeleton = c(.10, .20, .30, .40, .50),
  target = .30,
  utility_type = 3L,
  utility_scores = c(u00 = 0, u01 = 40, u10 = 60, u11 = 100),
  efficacy_threshold = .20,
  futility_cutoff = .85,
  min_eff_n_for_futility = 3L,
  apply_ipde_toxicity_rule = TRUE,
  apply_ipde_efficacy_rule = FALSE,
  seed = 113L
)

cat("\n========================================\n")
cat("JAGS TITE AIDE Phase I/II trial\n")
cat("========================================\n")
cat("Final MTD:", fit$final$MTD, "\n")
cat("Final OBD:", fit$final$OBD, "\n")
cat("Stop reason:", fit$final$stop_reason, "\n\n")

cat("Dose-level true probabilities\n")
print(data.frame(dose = seq_along(p_true), regular_toxicity = p_true,
                 regular_efficacy = e_true), row.names = FALSE)

cat("\nEvent-time cohort decisions\n")
print(fit$decision_log, row.names = FALSE)

cat("\nAll administrations (including recycled patients)\n")
print(fit$admin_phase12, row.names = FALSE)

cat("\nTITE n_eval decisions\n")
print(fit$n_eval_log, row.names = FALSE)

cat("\nIPDE eligibility decisions\n")
print(fit$recycling_decision_log, row.names = FALSE)
