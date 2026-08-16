aide_phase12_root <- function() {
  root <- getOption("tite_aide_root")
  if (is.null(root)) normalizePath(getwd(), winslash = "/", mustWork = FALSE) else root
}

aide_phase12_model_choices <- function() list(
  toxicity = c("discount_r", "additive_alpha"),
  efficacy = c("dose_specific_carryover", "shared_carryover",
               "dose_specific_previous_dose_additive", "previous_dose_additive")
)

aide_phase12_merge <- function(x, defaults) modifyList(defaults, x %||% list())
`%||%` <- function(x, y) if (is.null(x)) y else x

aide_phase12_config <- function(
    allocation = c("one_stage", "two_stage"), cohort_size = 3L, Nmax = 30L,
    s1_Max = 15L, N_s2 = Nmax, n_eval = cohort_size, m_U = 6L,
    cycle_max = 2L, start_dose = 1L,
    time = list(arrival_rate = .2, t0 = 0, T_assess = 28,
                dlt_dist = "uniform", efficacy_dist = "uniform"),
    toxicity = list(model = "discount_r", target = .30, cutoff = .95,
                    carryover_prior = c(1, 9), skeleton = NULL,
                    beta_prior_mean = 0, beta_prior_sd = sqrt(2),
                    n_chains = 2L, n_adapt = 500L, n_burnin = 500L,
                    n_iter = 2000L, thin = 1L),
    efficacy = list(model = "shared_carryover", prior = c(1, 1),
                    carryover_prior = c(1, 9), threshold = .20,
                    futility_cutoff = .90, min_eff_n_for_futility = 3L,
                    n_chains = 3L, n_adapt = 1000L, n_burnin = 1000L,
                    n_iter = 4000L, thin = 2L),
    monitoring = list(stage2_allocation = "highest_utility"),
    utility = list(type = 1L, lambda_T = 1, scores = c(0, 1, -1, 0),
                   tie_break = "lower_dose"),
    recycle = list(priority = "new_first", apply_individual_toxicity_risk = FALSE,
                   toxicity_ipde_overdose_cutoff = .95,
                   apply_individual_efficacy_benefit = FALSE,
                   efficacy_ipde_min_increment = 0),
    reporting = list(verbose = FALSE, store_raw_tables = TRUE), ...) {
  allocation <- match.arg(allocation)
  config <- list(
    design = list(allocation = allocation, cohort_size = as.integer(cohort_size),
                  Nmax = as.integer(Nmax), s1_Max = as.integer(s1_Max),
                  N_s2 = as.integer(N_s2), n_eval = as.integer(n_eval),
                  m_U = as.integer(m_U), cycle_max = as.integer(cycle_max),
                  start_dose = as.integer(start_dose)),
    time = aide_phase12_merge(time, list(arrival_rate = .2, t0 = 0, T_assess = 28,
                                         dlt_dist = "uniform", efficacy_dist = "uniform")),
    toxicity = aide_phase12_merge(toxicity, list(model = "discount_r", target = .30,
                                                 cutoff = .95, skeleton = NULL,
                                                 beta_prior_mean = 0, beta_prior_sd = sqrt(2),
                                                 n_chains = 2L, n_adapt = 500L, n_burnin = 500L,
                                                 n_iter = 2000L, thin = 1L,
                                                 carryover_prior = c(1, 9), backend = "jags")),
    efficacy = aide_phase12_merge(efficacy, list(model = "shared_carryover", prior = c(1, 1),
                                                 carryover_prior = c(1, 9), threshold = .20,
                                                 futility_cutoff = .90,
                                                 min_eff_n_for_futility = 3L,
                                                 n_chains = 3L, n_adapt = 1000L, n_burnin = 1000L,
                                                 n_iter = 4000L, thin = 2L,
                                                 backend = "jags")),
    monitoring = aide_phase12_merge(monitoring, list(
      stage2_allocation = "highest_utility"
    )),
    utility = aide_phase12_merge(utility, list(type = 1L, lambda_T = 1,
                                               scores = c(0, 1, -1, 0),
                                               tie_break = "lower_dose")),
    recycle = aide_phase12_merge(recycle, list(priority = "new_first",
                                               drop_new_trigger_on_n_eval = TRUE,
                                               retain_retreat_trigger_on_n_eval = TRUE,
                                               apply_individual_toxicity_risk = FALSE,
                                               toxicity_ipde_overdose_cutoff = .95,
                                               apply_individual_efficacy_benefit = FALSE,
                                               efficacy_ipde_min_increment = 0)),
    reporting = aide_phase12_merge(reporting, list(verbose = FALSE, store_raw_tables = TRUE))
  )
  aide_phase12_validate_config(config)
  config
}

aide_phase12_scenario <- function(p_true, e_true,
                                  toxicity_ipde_dgm = list(alpha_true = 0),
                                  efficacy_ipde_dgm = list(alpha_true = 0),
                                  endpoint_time_dgm = list()) {
  p_true <- as.numeric(p_true); e_true <- as.numeric(e_true)
  if (!length(p_true) || length(p_true) != length(e_true) ||
      any(!is.finite(c(p_true, e_true))) || any(c(p_true, e_true) < 0 | c(p_true, e_true) > 1))
    stop("p_true and e_true must have the same positive length and values in [0, 1].")
  tx <- aide_phase12_merge(toxicity_ipde_dgm, list(alpha_true = 0))
  ef <- aide_phase12_merge(efficacy_ipde_dgm, list(alpha_true = 0))
  if (any(!is.finite(c(tx$alpha_true, ef$alpha_true))) || any(c(tx$alpha_true, ef$alpha_true) < 0))
    stop("IPDE DGM alpha_true values must be finite and non-negative.")
  list(p_true = p_true, e_true = e_true, toxicity_ipde_dgm = tx,
       efficacy_ipde_dgm = ef, endpoint_time_dgm = endpoint_time_dgm, ndose = length(p_true))
}

aide_phase12_validate_config <- function(config, scenario = NULL) {
  d <- config$design; choices <- aide_phase12_model_choices()
  numeric_design <- unlist(d[c("cohort_size", "Nmax", "s1_Max", "N_s2", "n_eval", "m_U", "cycle_max", "start_dose")])
  if (!d$allocation %in% c("one_stage", "two_stage") ||
      any(!is.finite(numeric_design)) || any(unlist(d[c("cohort_size", "Nmax", "n_eval", "cycle_max", "start_dose")]) < 1))
    stop("Invalid design configuration.")
  if (d$s1_Max < 1L || d$s1_Max > d$Nmax || d$N_s2 < 1L || d$N_s2 > d$Nmax)
    stop("s1_Max and N_s2 must be positive administration thresholds no greater than Nmax.")
  if (!is.finite(config$time$T_assess) || config$time$T_assess <= 0 ||
      !is.finite(config$time$arrival_rate) || config$time$arrival_rate <= 0)
    stop("time requires one positive common T_assess and positive arrival_rate.")
  if (!config$toxicity$model %in% choices$toxicity || !config$efficacy$model %in% choices$efficacy)
    stop("Selected toxicity or efficacy TITE backend is unavailable.")
  if (!identical(config$recycle$priority, "new_first")) stop("The rebuilt TITE engine requires recycle$priority = 'new_first'.")
  if (!is.finite(config$recycle$toxicity_ipde_overdose_cutoff) || config$recycle$toxicity_ipde_overdose_cutoff <= 0 || config$recycle$toxicity_ipde_overdose_cutoff >= 1 ||
      !is.finite(config$recycle$efficacy_ipde_min_increment) || config$recycle$efficacy_ipde_min_increment < 0)
    stop("Individual recycle-risk thresholds are invalid.")
  if (!is.null(config$toxicity$prior))
    stop("toxicity$prior is not a CRM prior. Use beta_prior_mean/beta_prior_sd for the CRM curve and carryover_prior for r or alpha.")
  if (length(config$toxicity$carryover_prior) != 2L || length(config$efficacy$prior) != 2L ||
       any(c(config$toxicity$carryover_prior, config$efficacy$prior) <= 0)) stop("Carryover and efficacy priors must be positive beta parameters.")
  if (!is.finite(config$efficacy$threshold) || config$efficacy$threshold <= 0 || config$efficacy$threshold >= 1 ||
      !is.finite(config$efficacy$futility_cutoff) || config$efficacy$futility_cutoff <= 0 || config$efficacy$futility_cutoff >= 1 ||
      !is.finite(config$efficacy$min_eff_n_for_futility) || config$efficacy$min_eff_n_for_futility < 0 ||
      config$efficacy$min_eff_n_for_futility != as.integer(config$efficacy$min_eff_n_for_futility)) {
    stop("Efficacy futility settings are invalid.")
  }
  if (!is.null(config$toxicity$skeleton) &&
      (any(!is.finite(config$toxicity$skeleton)) || any(config$toxicity$skeleton <= 0 | config$toxicity$skeleton >= 1)))
    stop("toxicity$skeleton must contain probabilities strictly between 0 and 1.")
  if (!config$monitoring$stage2_allocation %in% c("highest_utility", "top2_randomized"))
    stop("monitoring$stage2_allocation must be 'highest_utility' or 'top2_randomized'.")
  if (!is.finite(config$toxicity$beta_prior_mean) || !is.finite(config$toxicity$beta_prior_sd) ||
      config$toxicity$beta_prior_sd <= 0 || any(as.integer(unlist(config$toxicity[c("n_chains", "n_adapt", "n_burnin", "n_iter", "thin")])) < c(1L, 0L, 0L, 1L, 1L)))
    stop("The CRM beta prior or JAGS MCMC controls are invalid.")
  if (any(as.integer(unlist(config$efficacy[c("n_chains", "n_adapt", "n_burnin", "n_iter", "thin")])) < c(1L, 0L, 0L, 1L, 1L)))
    stop("The efficacy JAGS MCMC controls are invalid.")
  if (!is.null(scenario) && d$start_dose > scenario$ndose) stop("start_dose exceeds the number of doses.")
  invisible(TRUE)
}
