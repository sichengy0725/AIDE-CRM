## Regression checks for the multicycle CRM posterior contract.
## Run with: Rscript test_crm_multicycle_posterior_contract.R

source("AIDE_CRM_helper_modified.R")

dat <- data.frame(
  id = c(1L, 2L),
  ncycle = c(1L, 1L),
  dose = c(1L, 2L),
  y = c(0L, 0L),
  type = c("new", "new")
)
posterior <- cbind(
  "p[1]" = c(0.10, 0.15),
  "p[2]" = c(0.20, 0.25),
  alpha = c(0.20, 0.30)
)
regular_fit <- list(
  p_hat = c(0.10, 0.20),
  theta_ipde_hat = c(0.10, 0.20),
  r_hat = NA_real_,
  prob_overtox = c(0.01, 0.02),
  prob_ipde_p2_over_target = 0.05,
  eliminated = c(0L, 0L),
  stop = 0L,
  earlystop = 0L,
  post = posterior
)
stopped_fit <- regular_fit
stopped_fit$eliminated <- c(1L, 1L)
stopped_fit$stop <- 1L
stopped_fit$earlystop <- 1L

## The MTD selector must retain its backend fit even when dose 1 stops.
stopped_selection <- select.mtd.crm(
  target = 0.30,
  dat = dat,
  ndose = 2L,
  skeleton = c(0.10, 0.20),
  r_model = "multicycle_additive",
  fit = stopped_fit
)
stopifnot(
  identical(stopped_selection$MTD, 99L),
  identical(stopped_selection$fit, stopped_fit)
)

## Every crm_move branch reached after a supplied fit keeps that same fit.
stopped_move <- crm_move(
  current_dose = 1L,
  ndose = 2L,
  dat = dat,
  skeleton = c(0.10, 0.20),
  r_model = "multicycle_additive",
  fit = stopped_fit
)
all_eliminated_move <- crm_move(
  current_dose = 1L,
  ndose = 2L,
  dat = dat,
  skeleton = c(0.10, 0.20),
  r_model = "multicycle_additive",
  elimi = c(1L, 1L),
  fit = regular_fit
)
regular_move <- crm_move(
  current_dose = 1L,
  ndose = 2L,
  dat = dat,
  skeleton = c(0.10, 0.20),
  r_model = "multicycle_additive",
  fit = regular_fit
)
stopifnot(
  identical(stopped_move$fit, stopped_fit),
  identical(all_eliminated_move$fit, regular_fit),
  identical(regular_move$fit, regular_fit)
)

## A normal selector result continues to supply draw-level p_j and alpha to
## the patient-specific gate; it is not reduced to posterior means.
regular_selection <- select.mtd.crm(
  target = 0.30,
  dat = dat,
  ndose = 2L,
  skeleton = c(0.10, 0.20),
  r_model = "multicycle_additive",
  fit = regular_fit
)
history <- crm_multicycle_history_data(dat, ndose = 2L)
gate <- crm_multicycle_recycle_toxicity_gate(
  fit = regular_selection,
  history = history,
  patient_id = 1L,
  next_dose = 2L,
  phi = 0.30,
  cutoff = 0.95
)
stopifnot(
  identical(regular_selection$fit, regular_fit),
  length(gate$q_next_draws) == nrow(posterior),
  isTRUE(all.equal(gate$q_next_draws, posterior[, "p[2]"] +
    posterior[, "alpha"] * posterior[, "p[1]"]))
)

## The non-TITE simulator should reuse this backend posterior at an unchanged
## administration history for both the CRM move and recycle screening.
source("AIDE_phase_I_II.R")
check_simulator_cache <- function() {
  originals <- mget(
    c(
      "crm_move", "select.mtd.crm", "aide_phase12_efficacy_summary",
      "aide_phase12_utility", "aide_phase12_load_efficacy_jags"
    ),
    envir = .GlobalEnv
  )
  on.exit(
    for (name in names(originals)) assign(name, originals[[name]], envir = .GlobalEnv),
    add = TRUE
  )

  posterior_refits <- 0L
  move_received_cached_fit <- logical(0)
  backend_fit <- regular_fit

  assign("select.mtd.crm", function(target, dat, ndose, fit = NULL, ...) {
    if (is.null(fit)) {
      posterior_refits <<- posterior_refits + 1L
      fit <- backend_fit
    }
    list(
      MTD = as.integer(ndose),
      phat = fit$p_hat,
      pj_iso = fit$p_hat,
      p_hat = fit$p_hat,
      theta_ipde_hat = fit$theta_ipde_hat,
      r_hat = fit$r_hat,
      prob_overtox = fit$prob_overtox,
      prob_ipde_p2_over_target = fit$prob_ipde_p2_over_target,
      eliminated = fit$eliminated,
      stop = fit$stop,
      earlystop = fit$earlystop,
      fit = fit
    )
  }, envir = .GlobalEnv)
  assign("crm_move", function(current_dose, ndose, fit = NULL, ...) {
    move_received_cached_fit <<- c(move_received_cached_fit, !is.null(fit))
    list(
      next_dose = as.integer(current_dose),
      action = "stay",
      p_hat = fit$p_hat,
      theta_ipde_hat = fit$theta_ipde_hat,
      r_hat = fit$r_hat,
      prob_overtox = fit$prob_overtox,
      prob_ipde_p2_over_target = fit$prob_ipde_p2_over_target,
      eliminated = fit$eliminated,
      stop_trial = FALSE,
      stop = 0L,
      earlystop = 0L,
      fit = fit
    )
  }, envir = .GlobalEnv)
  assign("aide_phase12_efficacy_summary", function(ndose, ...) {
    list(futility_eliminated = rep.int(0L, ndose))
  }, envir = .GlobalEnv)
  assign("aide_phase12_load_efficacy_jags", function(...) invisible(TRUE),
         envir = .GlobalEnv)
  assign("aide_phase12_utility", function(ndose, toxicity_estimate, ...) {
    list(
      toxicity = as.numeric(toxicity_estimate),
      efficacy = rep.int(0.50, ndose),
      utility = rev(seq_len(ndose))
    )
  }, envir = .GlobalEnv)

  trial <- simulate_AIDE_phase_I_II(
    p_true = c(0, 0),
    e_true = c(0, 0),
    p_ipde = c(0, 0),
    e_ipde = c(0, 0),
    model = "CRM",
    crm_r_model = "multicycle_additive",
    crm_skeleton = c(0.10, 0.20),
    allocation = "one_stage",
    N_s1 = 3L,
    Nmax = 6L,
    C = 3L,
    cycle_max = 2L,
    enrollment_scheme = "ipde_first",
    seed = 7L
  )
  stopifnot(
    identical(posterior_refits, 2L),
    length(move_received_cached_fit) > 0L,
    all(move_received_cached_fit),
    any(trial$admin$type == "retreat")
  )
}
check_simulator_cache()

## A terminal selector result must return no recycling candidates before the
## draw-by-draw gate is evaluated.
check_terminal_recycle_guard <- function() {
  originals <- mget(
    c(
      "crm_move", "select.mtd.crm", "aide_phase12_efficacy_summary",
      "aide_phase12_utility", "aide_phase12_load_efficacy_jags",
      "crm_multicycle_recycle_toxicity_gate"
    ),
    envir = .GlobalEnv
  )
  on.exit(
    for (name in names(originals)) assign(name, originals[[name]], envir = .GlobalEnv),
    add = TRUE
  )

  selector_calls <- 0L
  assign("select.mtd.crm", function(target, dat, ndose, fit = NULL, ...) {
    selector_calls <<- selector_calls + 1L
    if (is.null(fit)) fit <- regular_fit
    is_terminal <- selector_calls == 2L
    list(
      MTD = if (is_terminal) 99L else as.integer(ndose),
      phat = fit$p_hat,
      pj_iso = fit$p_hat,
      p_hat = fit$p_hat,
      theta_ipde_hat = fit$theta_ipde_hat,
      r_hat = fit$r_hat,
      prob_overtox = fit$prob_overtox,
      prob_ipde_p2_over_target = fit$prob_ipde_p2_over_target,
      ## Keep this test's trial running long enough to exercise eligibility;
      ## the terminal state itself is what blocks recycling.
      eliminated = c(0L, 0L),
      stop = as.integer(is_terminal),
      earlystop = as.integer(is_terminal),
      fit = fit
    )
  }, envir = .GlobalEnv)
  assign("crm_move", function(current_dose, ndose, fit = NULL, ...) {
    list(
      next_dose = as.integer(current_dose),
      action = "stay",
      p_hat = fit$p_hat,
      theta_ipde_hat = fit$theta_ipde_hat,
      r_hat = fit$r_hat,
      prob_overtox = fit$prob_overtox,
      prob_ipde_p2_over_target = fit$prob_ipde_p2_over_target,
      eliminated = fit$eliminated,
      stop_trial = FALSE,
      stop = 0L,
      earlystop = 0L,
      fit = fit
    )
  }, envir = .GlobalEnv)
  assign("aide_phase12_efficacy_summary", function(ndose, ...) {
    list(futility_eliminated = rep.int(0L, ndose))
  }, envir = .GlobalEnv)
  assign("aide_phase12_utility", function(ndose, toxicity_estimate, ...) {
    list(
      toxicity = as.numeric(toxicity_estimate),
      efficacy = rep.int(0.50, ndose),
      utility = rev(seq_len(ndose))
    )
  }, envir = .GlobalEnv)
  assign("aide_phase12_load_efficacy_jags", function(...) invisible(TRUE),
         envir = .GlobalEnv)
  assign("crm_multicycle_recycle_toxicity_gate", function(...) {
    stop("The recycle gate must not run after terminal toxicity selection.")
  }, envir = .GlobalEnv)

  trial <- simulate_AIDE_phase_I_II(
    p_true = c(0, 0),
    e_true = c(0, 0),
    p_ipde = c(0, 0),
    e_ipde = c(0, 0),
    model = "CRM",
    crm_r_model = "multicycle_additive",
    crm_skeleton = c(0.10, 0.20),
    allocation = "one_stage",
    N_s1 = 3L,
    Nmax = 6L,
    C = 3L,
    cycle_max = 2L,
    enrollment_scheme = "ipde_first",
    seed = 8L
  )
  stopifnot(!any(trial$admin$type == "retreat"))
}
check_terminal_recycle_guard()

cat("Multicycle CRM posterior-contract checks passed.\n")
