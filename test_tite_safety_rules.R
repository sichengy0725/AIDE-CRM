## Regression checks for the opt-in TITE suspension and random-CRM IPDE rules.
## Run with: Rscript test_tite_safety_rules.R

source("AIDE_CRM_helper_modified.R")
source("AIDE_CRM_helper_TITE.R")

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

## The probability must be calculated from paired posterior draws, not from
## the posterior means of r and p[2].
post <- cbind(
  "r" = c(0.10, 0.90),
  "p[2]" = c(0.10, 0.10)
)
safety <- crm_random_ipde_p2_safety(post, target = 0.50)
assert_true(
  isTRUE(all.equal(safety$theta_ipde_p2_hat, 0.55)),
  "Incorrect posterior mean for r + (1-r) p[2]."
)
assert_true(
  isTRUE(all.equal(safety$prob_ipde_p2_over_target, 0.50)),
  "Incorrect posterior probability for r + (1-r) p[2] > target."
)

## Replace the JAGS-dependent decision and final-selection functions with
## deterministic versions so the event-clock rules can be tested with base R.
crm_move <- function(current_dose, ndose, ...) {
  list(
    next_dose = min(current_dose + 1L, ndose),
    action = "escalate",
    eliminated = rep(0L, ndose),
    prob_ipde_p2_over_target = 0.96
  )
}

select.mtd.crm <- function(target, dat, ndose, ...) {
  list(
    MTD = 1L,
    phat = rep(target, ndose),
    pj_iso = rep(target, ndose),
    p_hat = rep(target, ndose),
    theta_ipde_hat = rep(NA_real_, ndose),
    r_hat = NA_real_,
    eliminated = rep(0L, ndose),
    earlystop = 0L,
    stop = 0L
  )
}

aide_crm_tite_draw_dlt <- function(pi_val, t_start, window, ...) {
  list(y = 0L, t_tox = Inf, t_eval = t_start + window)
}

sim_args <- list(
  p_true = c(0, 0),
  target = 0.30,
  N_pat = 10L,
  Nmax_eff = 10L,
  cohortsize = 3L,
  T_assess = 1,
  cycle_max = 2L,
  arrival_rate = 100,
  n_eval_escalate = 3L,
  crm_r_model = "random",
  crm_skeleton = c(0.20, 0.30),
  seed = 71L
)

## With the suspension rule off, arrivals are retained in the backlog.
queued <- do.call(
  simulate_AIDE_CRM_TITE_design,
  c(sim_args, list(drop_new_patients_when_suspended = FALSE))
)

## With it on, arrivals after the trial enters suspension are dropped.
dropped <- do.call(
  simulate_AIDE_CRM_TITE_design,
  c(sim_args, list(drop_new_patients_when_suspended = TRUE))
)
assert_true(queued$dropped_new_patients == 0L, "Queued suspension run dropped a new patient.")
assert_true(dropped$dropped_new_patients > 0L, "Suspended new arrivals were not dropped.")
assert_true(
  nrow(dropped$admin) < nrow(queued$admin),
  "Dropping suspended arrivals did not change the treated-patient count."
)

## The random-CRM gate rejects all queued and future IPDE events when its
## posterior probability exceeds the configured cutoff.
gated <- do.call(
  simulate_AIDE_CRM_TITE_design,
  c(
    sim_args,
    list(
      drop_new_patients_when_suspended = FALSE,
      random_crm_ipde_safety = TRUE,
      random_crm_ipde_prob_cutoff = 0.95
    )
  )
)
assert_true(gated$dropped_ipde_events > 0L, "Random-CRM gate did not reject IPDE events.")
assert_true(
  !any(gated$admin$type == "retreat"),
  "An IPDE assignment was made despite an unsafe random-CRM posterior."
)

cat("TITE safety-rule regression checks passed.\n")
