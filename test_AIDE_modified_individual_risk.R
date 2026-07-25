## Regression checks for the optional random-CRM IPDE individual-risk rule.
## Run with: Rscript test_AIDE_modified_individual_risk.R

source("AIDE_BOIN_helper.R")
source("AIDE_modified.R")

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

## The criterion must use paired posterior draws of r and p[next_dose].
post <- cbind(
  "r" = c(0.10, 0.90),
  "p[2]" = c(0.10, 0.10)
)
paired_gate <- aide_random_crm_recycle_toxicity_gate(
  post = post,
  next_dose = 2L,
  phi = 0.50,
  cutoff = 0.75
)
assert_true(
  isTRUE(all.equal(paired_gate$theta_posterior_mean, 0.55)),
  "Individual-risk gate did not use paired posterior draws."
)
assert_true(
  isTRUE(all.equal(paired_gate$probability_over_phi, 0.50)),
  "Individual-risk gate calculated the wrong posterior exceedance probability."
)

## Replace the JAGS-dependent CRM functions with deterministic versions.  The
## first post-cohort decision moves to dose 2, where the IPDE rule is tested.
crm_move <- function(current_dose, ndose, ...) {
  list(
    next_dose = if (current_dose < ndose) current_dose + 1L else current_dose,
    action = "escalate",
    mu_hat = 0.20,
    n_eff = 3L,
    stop_trial = FALSE,
    earlystop = 0L,
    eliminated = rep(0L, ndose)
  )
}

unsafe_post <- cbind(
  "r" = rep(0.90, 100L),
  "p[1]" = rep(0.10, 100L),
  "p[2]" = rep(0.10, 100L)
)
crm_fit <- function(...) list(post = unsafe_post)

select.mtd.crm <- function(target, dat, ndose, ...) {
  list(
    MTD = 1L,
    phat = rep(target, ndose),
    pj_iso = rep(target, ndose),
    eliminated = rep(0L, ndose),
    r_hat = NA_real_,
    earlystop = 0L,
    stop = 0L
  )
}

sim_args <- list(
  N_pat = 9L,
  Nmax_eff = 9L,
  C = 3L,
  T_assess = 1,
  cycle_max = 2L,
  arrival_rate = 100,
  new_pat_first = 2L,
  p_true = c(0, 0),
  model = "CRM",
  crm_r_model = "random",
  crm_skeleton = c(0.15, 0.30),
  seed = 41L
)

gated <- do.call(
  simulate_AIDE_design,
  c(sim_args, list(apply_random_crm_recycle_toxicity_rule = TRUE))
)
assert_true(
  !any(gated$admin$type == "retreat"),
  "An unsafe IPDE assignment was made while the individual-risk rule was enabled."
)
assert_true(
  nrow(gated$ipde_risk_log) > 0L &&
    all(!gated$ipde_risk_log$allowed) &&
    all(gated$ipde_risk_log$probability_over_target > 0.95),
  "The unsafe random-CRM posterior was not recorded as an IPDE prohibition."
)

ungated <- do.call(
  simulate_AIDE_design,
  c(sim_args, list(apply_random_crm_recycle_toxicity_rule = FALSE))
)
assert_true(
  any(ungated$admin$type == "retreat"),
  "Disabling the individual-risk rule did not restore eligible IPDE assignments."
)
assert_true(
  nrow(ungated$ipde_risk_log) > 0L &&
    all(!ungated$ipde_risk_log$rule_enabled) &&
    all(ungated$ipde_risk_log$allowed),
  "The IPDE risk log did not record the disabled rule correctly."
)

cat("AIDE individual-risk regression checks passed.\n")
