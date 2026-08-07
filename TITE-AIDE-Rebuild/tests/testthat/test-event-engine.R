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

# A proposed escalation with no evaluated administrations blocks a new trigger.
state <- aide_phase12_initialize_state(cfg, sce, seed = 2L)
state$cohort$open <- FALSE
decision <- list(action = "escalate")
gate <- aide_apply_escalation_gate(state, decision, cfg, "new")
stopifnot(gate$blocked, gate$new_patient_dropped, gate$trigger_disposition == "new_dropped_n_eval")

# The corresponding retreat trigger is retained rather than dropped.
state <- aide_queue_add(state, 0, 1, "retreat", 99L, 1L)
gate <- aide_apply_escalation_gate(state, decision, cfg, "retreat", 1L)
stopifnot(gate$blocked, gate$retreat_opportunity_retained, gate$trigger_disposition == "retreat_retained_n_eval")

# Stay is now gated by n_eval; de-escalation remains immediately permissible.
stopifnot(aide_apply_escalation_gate(state, list(action = "stay"), cfg, "new")$blocked)
stopifnot(!aide_apply_escalation_gate(state, list(action = "de_escalate"), cfg, "new")$blocked)

# Optional individual-risk recycle screens use the open cohort's frozen dose.
risk_cfg <- aide_phase12_config(recycle = list(
  apply_individual_toxicity_risk = TRUE, toxicity_ipde_overdose_cutoff = .90,
  apply_individual_efficacy_benefit = TRUE, efficacy_ipde_min_increment = .05
))
risk_state <- aide_phase12_initialize_state(risk_cfg, sce, seed = 3L)
risk_state$cohort$next_dose <- 1L
risk_state$cohort$risk_context <- list(
  prob_toxicity_ipde_overdose = c(.96, .2, .1),
  p_efficacy = c(.20, .30, .40), p_efficacy_ipde = c(.22, .32, .50)
)
screen <- aide_retreat_individual_risk_screen(risk_state, 1L, risk_cfg)
stopifnot(!screen$allowed, screen$reason == "blocked_individual_toxicity_risk")
risk_state$cohort$risk_context$prob_toxicity_ipde_overdose[1L] <- .10
screen <- aide_retreat_individual_risk_screen(risk_state, 1L, risk_cfg)
stopifnot(!screen$allowed, screen$reason == "blocked_individual_efficacy_benefit")
risk_state$cohort$risk_context$p_efficacy_ipde[1L] <- .30
stopifnot(aide_retreat_individual_risk_screen(risk_state, 1L, risk_cfg)$allowed)

message("event-engine tests passed")
