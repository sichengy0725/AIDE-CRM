get_oc_sim_AIDE_phase_I_II <- function(scenario, config, ntrial = 100L, seed = 1L, store_raw_tables = FALSE) {
  trials <- lapply(seq_len(ntrial), function(i) simulate_AIDE_phase_I_II_event(config, scenario, seed = seed + i - 1L))
  mtd <- vapply(trials, function(x) x$final$MTD %||% NA_integer_, integer(1))
  obd <- vapply(trials, function(x) x$final$OBD %||% aide_phase12_no_obd, integer(1))
  ## Normalize legacy trial objects that recorded no OBD as NA.
  obd[is.na(obd)] <- aide_phase12_no_obd
  ndose <- scenario$ndose
  dose_counts <- function(x, type = NULL) {
    dat <- x$admin
    if (!is.null(type)) dat <- dat[dat$assignment_type == type, , drop = FALSE]
    tabulate(dat$dose, nbins = ndose)
  }
  unique_counts <- function(x) vapply(seq_len(ndose), function(j) length(unique(x$admin$patient_id[x$admin$dose == j])), numeric(1))
  duration <- vapply(trials, function(x) if (nrow(x$admin)) max(x$admin$assessment_end) else 0, numeric(1))
  tox_p <- do.call(rbind, lapply(trials, function(x) x$final$toxicity$p_regular_mean))
  eff_p <- do.call(rbind, lapply(trials, function(x) x$final$efficacy$p_regular_mean))
  utility <- do.call(rbind, lapply(trials, function(x) x$final$utilities))
  r_hat <- do.call(rbind, lapply(trials, function(x) rep(x$final$toxicity$posterior_carryover_mean, ndose)))
  summaries <- data.frame(trial = seq_len(ntrial), MTD = mtd, OBD = obd,
    administrations = vapply(trials, function(x) nrow(x$admin), integer(1)),
    duration = duration,
    n_eval_blocks = vapply(trials, function(x) sum(x$n_eval_log$blocked), integer(1)))
  out <- list(ntrial = as.integer(ntrial), trial_summary = summaries,
              MTD_by_trial = mtd, OBD_by_trial = obd,
              early_stop_by_trial = vapply(trials, function(x) x$stop_reason != "administration_cap", logical(1)),
              n_administrations_by_trial = summaries$administrations,
              n_unique_patients_by_trial = vapply(trials, function(x) length(unique(x$admin$patient_id)), integer(1)),
              duration_by_trial = duration,
              n_by_dose_by_trial = do.call(rbind, lapply(trials, dose_counts)),
              unique_n_by_dose_by_trial = do.call(rbind, lapply(trials, unique_counts)),
              nipde_by_dose_by_trial = do.call(rbind, lapply(trials, dose_counts, type = "retreat")),
              crm_pj_by_trial = tox_p, efficacy_pj_by_trial = eff_p,
              utility_by_trial = utility, r_hat_by_trial = r_hat,
              n_eval_blocks_by_trial = summaries$n_eval_blocks,
              MTD_selection = prop.table(table(factor(mtd, levels = seq_len(ndose)), useNA = "ifany")),
              OBD_selection = prop.table(table(factor(obd, levels = seq_len(ndose)), useNA = "ifany")),
              No_OBD_selection_percent = 100 * mean(obd == aide_phase12_no_obd),
              mean_administrations = mean(summaries$administrations), mean_duration = mean(summaries$duration))
  if (store_raw_tables) out$trials <- trials
  out
}
