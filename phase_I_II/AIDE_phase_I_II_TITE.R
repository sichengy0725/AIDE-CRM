## ============================================================
## AIDE Phase I/II TITE design (JAGS)
##
## This is the time-to-event counterpart to AIDE_phase_I_II.R.  It keeps the
## Phase I/II public interface, result fields, OC summaries, and one-/two-
## stage choices while delegating the event queue and delayed-outcome JAGS
## models to TITE-AIDE-Rebuild.  No integral approximation is available here.
##
## TITE rules inherited from the rebuild:
##   * toxicity and efficacy use partial follow-up at each decision time;
##   * a closed cohort reopens only for a new arrival or eligible retreat;
##   * new arrivals have priority over recycled patients;
##   * n_eval blocks stay/escalate, but not de-escalation; and
##   * every opened cohort has one frozen dose.
## ============================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

.aide_phase12_tite_source_files <- vapply(sys.frames(), function(frame) {
  value <- frame$ofile %||% ""
  if (length(value)) as.character(value)[1L] else ""
}, character(1))
.aide_phase12_tite_this_file <- tail(
  .aide_phase12_tite_source_files[
    basename(.aide_phase12_tite_source_files) == "AIDE_phase_I_II_TITE.R"
  ],
  1L
)
if (!length(.aide_phase12_tite_this_file) || !nzchar(.aide_phase12_tite_this_file)) {
  .aide_phase12_tite_candidates <- c(
    file.path(getwd(), "phase_I_II", "AIDE_phase_I_II_TITE.R"),
    file.path(getwd(), "AIDE_phase_I_II_TITE.R")
  )
  .aide_phase12_tite_this_file <- .aide_phase12_tite_candidates[
    file.exists(.aide_phase12_tite_candidates)
  ][1L]
}
if (is.na(.aide_phase12_tite_this_file) || !file.exists(.aide_phase12_tite_this_file)) {
  stop("Cannot locate AIDE_phase_I_II_TITE.R.")
}
.aide_phase12_tite_dir <- dirname(normalizePath(.aide_phase12_tite_this_file,
                                                 winslash = "/", mustWork = TRUE))
.aide_phase12_tite_rebuild <- file.path(dirname(.aide_phase12_tite_dir),
                                         "TITE-AIDE-Rebuild", "TITE-AIDE.R")
if (!file.exists(.aide_phase12_tite_rebuild)) {
  stop("The JAGS TITE implementation is missing: ", .aide_phase12_tite_rebuild)
}
source(.aide_phase12_tite_rebuild)

aide_phase12_tite_require_jags <- function() {
  if (!requireNamespace("rjags", quietly = TRUE) ||
      !requireNamespace("coda", quietly = TRUE)) {
    stop("The TITE Phase I/II design requires rjags, coda, and a working JAGS installation.")
  }
  invisible(TRUE)
}

aide_phase12_tite_model_map <- function(carryover_model = NULL,
                                         toxicity_model = NULL,
                                         crm_r_model = NULL,
                                         efficacy_model = NULL) {
  if (!is.null(carryover_model)) {
    carryover_model <- match.arg(
      carryover_model,
      c("discount_shared", "discount_dose_specific",
        "additive_shared", "additive_dose_specific")
    )
    mapped <- switch(
      carryover_model,
      discount_shared = list(toxicity = "discount_r", efficacy = "shared_carryover"),
      discount_dose_specific = list(toxicity = "discount_r", efficacy = "dose_specific_carryover"),
      additive_shared = list(toxicity = "additive_alpha", efficacy = "previous_dose_additive"),
      additive_dose_specific = list(toxicity = "additive_alpha", efficacy = "dose_specific_previous_dose_additive")
    )
  } else {
    if (is.null(toxicity_model)) {
      crm_r_model <- crm_r_model %||% "previous_dose"
      toxicity_model <- switch(
        crm_r_model,
        random = "discount_r",
        discount_r = "discount_r",
        previous_dose = "additive_alpha",
        previous_dose_additive = "additive_alpha",
        additive_alpha = "additive_alpha",
        stop("TITE-AIDE supports crm_r_model = 'random' or 'previous_dose'.")
      )
    }
    toxicity_model <- match.arg(toxicity_model, c("discount_r", "additive_alpha"))
    mapped <- list(
      toxicity = toxicity_model,
      efficacy = if (toxicity_model == "discount_r") {
        "shared_carryover"
      } else {
        "previous_dose_additive"
      }
    )
  }
  if (!is.null(efficacy_model)) mapped$efficacy <- match.arg(
    efficacy_model,
    c("dose_specific_carryover", "shared_carryover",
      "dose_specific_previous_dose_additive", "previous_dose_additive")
  )
  mapped
}

aide_phase12_tite_admin_view <- function(admin) {
  data.frame(
    row_id = admin$admin_id,
    id = admin$patient_id,
    t_arrival = admin$t_arrival,
    t_start = admin$t_start,
    t_tox = admin$t_dlt,
    t_eff = admin$t_response,
    t_eval = admin$assessment_end,
    dose = admin$dose,
    y = admin$dlt_final,
    eff = admin$eff_final,
    true_toxicity_probability = admin$true_p_tox,
    true_efficacy_probability = admin$true_p_eff,
    ncycle = admin$cycle,
    cohort = admin$cohort_id,
    stage = admin$stage,
    type = admin$assignment_type,
    stringsAsFactors = FALSE
  )
}

aide_phase12_tite_attach_phase12_fields <- function(trial) {
  admin_view <- aide_phase12_tite_admin_view(trial$admin)
  final <- trial$final
  n_admin <- nrow(trial$admin)
  final$utility <- list(
    toxicity = final$toxicity$p_regular_mean,
    efficacy = final$efficacy$p_regular_mean,
    utility = final$utilities,
    utility_type = trial$config$utility$type,
    lambda_T = trial$config$utility$lambda_T,
    utility_scores = trial$config$utility$scores
  )
  final$earlystop <- as.integer(!trial$stop_reason %in%
    c("administration_cap", "stage2_dose_cap"))
  final$stop_reason <- trial$stop_reason
  final$n_admin <- n_admin
  final$n_unique_patients <- length(unique(trial$admin$patient_id))
  final$trial_time <- if (n_admin) {
    max(trial$admin$assessment_end) - min(trial$admin$t_arrival)
  } else {
    NA_real_
  }
  trial$admin_phase12 <- admin_view
  trial$recycling_decision_log <- trial$retreat_log
  trial$final <- final
  trial
}

simulate_AIDE_phase_I_II_TITE <- function(
    p_true,
    e_true,
    toxicity_ipde_alpha = 0,
    efficacy_ipde_alpha = 0,
    allocation = c("two_stage", "one_stage"),
    Nmax = 30L,
    N_s1 = 15L,
    N_s2 = Nmax,
    C = 3L,
    cycle_max = 2L,
    n_eval = C,
    m_U = 6L,
    stage2_allocation = NULL,
    enrollment_scheme = "continuous",
    arrival_rate = 1 / 28,
    t0 = 0,
    T_assess = 28,
    dlt_dist = c("uniform", "front_loaded"),
    efficacy_dist = dlt_dist,
    target = .30,
    cutoff = .95,
    carryover_model = NULL,
    toxicity_model = NULL,
    crm_r_model = "previous_dose",
    crm_skeleton = NULL,
    crm_alpha_sd = sqrt(2),
    crm_a_r = 1,
    crm_b_r = 9,
    crm_n_chains = 2L,
    crm_n_adapt = 500L,
    crm_n_burnin = 500L,
    crm_n_iter = 2000L,
    crm_thin = 1L,
    efficacy_prior = c(1, 1),
    efficacy_carryover_prior = c(1, 9),
    efficacy_model = NULL,
    efficacy_n_chains = 3L,
    efficacy_n_adapt = 1000L,
    efficacy_n_burnin = 1000L,
    efficacy_n_iter = 4000L,
    efficacy_thin = 2L,
    efficacy_threshold = .20,
    futility_cutoff = .95,
    min_eff_n_for_futility = 0L,
    utility_type = 3L,
    lambda_T = 1,
    utility_scores = c(u00 = 0, u01 = 40, u10 = 60, u11 = 100),
    apply_ipde_toxicity_rule = TRUE,
    ipde_toxicity_cutoff = cutoff,
    apply_ipde_efficacy_rule = FALSE,
    ipde_efficacy_delta = .05,
    seed = NULL,
    verbose = FALSE) {
  aide_phase12_tite_require_jags()
  allocation <- match.arg(allocation)
  dlt_dist <- match.arg(dlt_dist)
  efficacy_dist <- match.arg(efficacy_dist, c("uniform", "front_loaded"))
  if (is.null(stage2_allocation)) {
    stage2_allocation <- if (allocation == "two_stage") "top2_randomized" else "highest_utility"
  }
  stage2_allocation <- match.arg(stage2_allocation,
                                 c("highest_utility", "top2_randomized"))
  if (!identical(enrollment_scheme, "continuous")) {
    stop("The TITE rebuild implements continuous new-patient-first enrollment only.")
  }
  mapped <- aide_phase12_tite_model_map(
    carryover_model = carryover_model,
    toxicity_model = toxicity_model,
    crm_r_model = crm_r_model,
    efficacy_model = efficacy_model
  )
  config <- aide_phase12_config(
    allocation = allocation,
    cohort_size = as.integer(C),
    Nmax = as.integer(Nmax),
    s1_Max = as.integer(N_s1),
    N_s2 = as.integer(N_s2),
    n_eval = as.integer(n_eval),
    m_U = as.integer(m_U),
    cycle_max = as.integer(cycle_max),
    time = list(arrival_rate = arrival_rate, t0 = t0, T_assess = T_assess,
                dlt_dist = dlt_dist, efficacy_dist = efficacy_dist),
    toxicity = list(
      model = mapped$toxicity, target = target, cutoff = cutoff,
      skeleton = crm_skeleton, beta_prior_mean = 0,
      beta_prior_sd = crm_alpha_sd, carryover_prior = c(crm_a_r, crm_b_r),
      n_chains = as.integer(crm_n_chains), n_adapt = as.integer(crm_n_adapt),
      n_burnin = as.integer(crm_n_burnin), n_iter = as.integer(crm_n_iter),
      thin = as.integer(crm_thin)
    ),
    efficacy = list(
      model = mapped$efficacy, prior = efficacy_prior,
      carryover_prior = efficacy_carryover_prior, threshold = efficacy_threshold,
      futility_cutoff = futility_cutoff,
      min_eff_n_for_futility = as.integer(min_eff_n_for_futility),
      n_chains = as.integer(efficacy_n_chains), n_adapt = as.integer(efficacy_n_adapt),
      n_burnin = as.integer(efficacy_n_burnin), n_iter = as.integer(efficacy_n_iter),
      thin = as.integer(efficacy_thin)
    ),
    utility = list(type = as.integer(utility_type), lambda_T = lambda_T,
                   scores = utility_scores),
    monitoring = list(stage2_allocation = stage2_allocation),
    recycle = list(
      priority = "new_first",
      apply_individual_toxicity_risk = isTRUE(apply_ipde_toxicity_rule),
      toxicity_ipde_overdose_cutoff = ipde_toxicity_cutoff,
      apply_individual_efficacy_benefit = isTRUE(apply_ipde_efficacy_rule),
      efficacy_ipde_min_increment = ipde_efficacy_delta
    ),
    reporting = list(verbose = isTRUE(verbose))
  )
  scenario <- aide_phase12_scenario(
    p_true = p_true,
    e_true = e_true,
    toxicity_ipde_dgm = list(alpha_true = toxicity_ipde_alpha),
    efficacy_ipde_dgm = list(alpha_true = efficacy_ipde_alpha)
  )
  aide_phase12_tite_attach_phase12_fields(
    simulate_AIDE_phase_I_II_event(config, scenario, seed = seed %||% 1L)
  )
}

get_oc_sim_AIDE_phase_I_II_TITE <- function(p_true, e_true, ntrial = 1000L,
                                             seed = 1L, store_raw = FALSE, ...) {
  ntrial <- as.integer(ntrial)
  if (length(ntrial) != 1L || is.na(ntrial) || ntrial < 1L)
    stop("ntrial must be a positive integer.")
  ndose <- length(p_true)
  trials <- lapply(seq_len(ntrial), function(i) {
    simulate_AIDE_phase_I_II_TITE(p_true = p_true, e_true = e_true,
                                  seed = seed + i - 1L, ...)
  })
  mtd <- vapply(trials, function(x) x$final$MTD, integer(1))
  obd <- vapply(trials, function(x) x$final$OBD, integer(1))
  stopped <- vapply(trials, function(x) x$final$earlystop, integer(1))
  dose_counts <- function(x, type = NULL) {
    dat <- x$admin
    if (!is.null(type)) dat <- dat[dat$assignment_type == type, , drop = FALSE]
    tabulate(dat$dose, nbins = ndose)
  }
  unique_counts <- function(x) vapply(seq_len(ndose), function(dose) {
    length(unique(x$admin$patient_id[x$admin$dose == dose]))
  }, numeric(1))
  extract_final <- function(x, name) {
    value <- x[[name]]
    if (is.null(value)) rep(NA_real_, ndose) else as.numeric(value)
  }
  tox <- do.call(rbind, lapply(trials, function(x) extract_final(x$final$toxicity, "p_regular_mean")))
  eff <- do.call(rbind, lapply(trials, function(x) extract_final(x$final$efficacy, "p_regular_mean")))
  utility <- do.call(rbind, lapply(trials, function(x) as.numeric(x$final$utilities)))
  carry <- vapply(trials, function(x) x$final$toxicity$posterior_carryover_mean, numeric(1))
  duration <- vapply(trials, function(x) x$final$trial_time, numeric(1))
  administrations <- vapply(trials, function(x) x$final$n_admin, integer(1))
  unique_patients <- vapply(trials, function(x) x$final$n_unique_patients, integer(1))
  out <- list(
    p_true = as.numeric(p_true), e_true = as.numeric(e_true), ntrial = ntrial,
    MTD_by_trial = mtd, OBD_by_trial = obd, early_stop_by_trial = stopped,
    n_administrations_by_trial = administrations,
    n_unique_patients_by_trial = unique_patients,
    duration_by_trial = duration,
    n_by_dose_by_trial = do.call(rbind, lapply(trials, dose_counts)),
    unique_n_by_dose_by_trial = do.call(rbind, lapply(trials, unique_counts)),
    nipde_by_dose_by_trial = do.call(rbind, lapply(trials, dose_counts, type = "retreat")),
    crm_pj_by_trial = tox, efficacy_pj_by_trial = eff,
    utility_by_trial = utility,
    r_hat_by_trial = matrix(rep(carry, ndose), ncol = ndose),
    MTD_selection_percent = 100 * tabulate(mtd[!is.na(mtd) & mtd >= 1L & mtd <= ndose], nbins = ndose) / ntrial,
    OBD_selection_percent = 100 * tabulate(obd[!is.na(obd) & obd >= 1L & obd <= ndose], nbins = ndose) / ntrial,
    early_stop_percent = 100 * mean(stopped),
    mean_administrations = mean(administrations),
    mean_administrations_by_dose = colMeans(do.call(rbind, lapply(trials, dose_counts))),
    mean_unique_patients = mean(unique_patients),
    mean_unique_patients_by_dose = colMeans(do.call(rbind, lapply(trials, unique_counts))),
    mean_ipde_doses_by_dose = colMeans(do.call(rbind, lapply(trials, dose_counts, type = "retreat"))),
    mean_duration = mean(duration, na.rm = TRUE),
    trial_summary = data.frame(
      trial = seq_len(ntrial), MTD = mtd, OBD = obd,
      administrations = administrations, unique_patients = unique_patients,
      duration = duration, early_stop = stopped
    )
  )
  if (isTRUE(store_raw)) out$raw <- trials
  out
}
