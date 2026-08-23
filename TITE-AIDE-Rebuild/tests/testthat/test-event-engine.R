.this_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1L])
source(file.path(dirname(dirname(dirname(normalizePath(.this_file)))), "TITE-AIDE.R"))

cfg <- aide_phase12_config(Nmax = 9L, cohort_size = 3L, n_eval = 3L)
sce <- aide_phase12_scenario(c(.05, .15, .30), c(.10, .30, .45),
                             list(alpha_true = .10), list(alpha_true = .10))

# Starting cohort is exempt from n_eval and stays dose homogeneous.
if (requireNamespace("rjags", quietly = TRUE) && requireNamespace("coda", quietly = TRUE)) {
  trial <- simulate_AIDE_phase_I_II_event(cfg, sce, seed = 41L)
  stopifnot(nrow(trial$admin) > 0L, all(trial$admin$dose == trial$admin$decision_next_dose))
}

# With no open cohort, a model-based decision is unavailable until the
# required number of fully evaluated DLT outcomes is present. The triggering
# new patient remains eligible to wait; it is not dropped.
state <- aide_phase12_initialize_state(cfg, sce, seed = 2L)
state$cohort$open <- FALSE
decision <- list(action = "escalate")
gate <- aide_apply_escalation_gate(state, decision, cfg, "new")
stopifnot(gate$blocked, !gate$new_patient_dropped,
          gate$trigger_disposition == "decision_unavailable_n_eval")

# A retreat trigger is also retained, and the gate applies before every
# proposed model decision, including de-escalation.
state <- aide_queue_add(state, 0, 1, "retreat", 99L, 1L)
gate <- aide_apply_escalation_gate(state, decision, cfg, "retreat", 1L)
stopifnot(gate$blocked, gate$retreat_opportunity_retained,
          gate$trigger_disposition == "decision_unavailable_n_eval")

# The same closed-cohort gate applies to stay and de-escalation decisions.
stopifnot(aide_apply_escalation_gate(state, list(action = "stay"), cfg, "new")$blocked)
stopifnot(aide_apply_escalation_gate(state, list(action = "de_escalate"), cfg, "new")$blocked)

# Non-multicycle TITE models have no individual toxicity recycle screen.
risk_cfg <- aide_phase12_config(recycle = list(
  apply_individual_toxicity_risk = TRUE, toxicity_ipde_overdose_cutoff = .90
))
risk_state <- aide_phase12_initialize_state(risk_cfg, sce, seed = 3L)
risk_state$cohort$next_dose <- 1L
risk_state$cohort$risk_context <- list(
  prob_toxicity_ipde_overdose = c(.96, .2, .1),
  p_efficacy = c(.20, .30, .40), p_efficacy_ipde = c(.22, .32, .50)
)
screen <- aide_retreat_individual_risk_screen(risk_state, 1L, risk_cfg)
stopifnot(screen$allowed, screen$reason == "eligible")

message("event-engine tests passed")
