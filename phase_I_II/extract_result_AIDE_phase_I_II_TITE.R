## ============================================================
## Extract JAGS TITE Phase I/II operating characteristics.
##
## The settings below intentionally mirror run_oc_AIDE_phase_I_II_TITE.R,
## matching the existing non-TITE Phase I/II extractor structure.
## ============================================================

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
if (length(script_file) == 1L && file.exists(script_file)) {
  phase12_dir <- dirname(normalizePath(script_file, winslash = "/"))
} else {
  phase12_dir <- normalizePath("phase_I_II", winslash = "/", mustWork = TRUE)
}
project_root <- dirname(phase12_dir)
setwd(project_root)
source(file.path(phase12_dir, "AIDE_phase_I_II_TITE.R"))

## ============================================================
## User settings: copy these from run_oc_AIDE_phase_I_II_TITE.R
## ============================================================

scenario_file <- "Set_5dose_adaptive_r_37_true_MTD_OBD_summary_lambda1.csv"
scenario_set_name <- tools::file_path_sans_ext(basename(scenario_file))
results_root <- file.path(project_root, paste0("oc_results_cluster_phase_I_II_TITE_", scenario_set_name))
out_dir <- file.path(project_root, "OC_summary_phase_I_II_TITE")
jobs_expected <- 1:1000
ntrial_per_job <- 1L
scenario_id_list <- 1:38

two_stage_sizes <- data.frame(
  allocation = "two_stage", Nmax = 30L, N_s1 = 6L, N_s2 = 30L,
  stringsAsFactors = FALSE
)
one_stage_sizes <- data.frame(
  allocation = "one_stage", Nmax = 30L, N_s1 = 30L, N_s2 = 30L,
  stringsAsFactors = FALSE
)
design_size_grid <- rbind(two_stage_sizes, one_stage_sizes)
allocation_mode_grid <- data.frame(
  allocation_mode_id = c(1L, 2L),
  allocation_mode = c("upward_no_skipping", "one_level_toward_obd"),
  stringsAsFactors = FALSE
)
model_grid <- data.frame(
  model_id = "multicycle_additive",
  carryover_model = "multicycle_additive",
  crm_r_model = "multicycle_additive",
  efficacy_model = "multicycle_additive",
  stringsAsFactors = FALSE
)
crm_prior_grid <- data.frame(crm_a_r = .15, crm_b_r = .85, stringsAsFactors = FALSE)
efficacy_prior_grid <- data.frame(
  efficacy_a = .5, efficacy_b = .5, carry_a = .15, carry_b = .85,
  stringsAsFactors = FALSE
)
utility_grid <- rbind(
  data.frame(utility_type = 2L, lambda_T = .3),
  data.frame(utility_type = 3L, lambda_T = 1)
)
arrival_grid <- data.frame(arrival_rate = 1 / 56)
fmt_truth_setting <- function(x) {
  out <- format(round(as.numeric(x), 4L), scientific = FALSE, trim = TRUE)
  out <- sub("(\\.[0-9]*?)0+$", "\\1", out)
  out <- sub("\\.$", "", out)
  gsub("\\.", "p", gsub("-", "m", out))
}
alpha_grid <- data.frame(
  toxicity_ipde_alpha = c(0, .3, .6, .9),
  efficacy_ipde_alpha = c(0, .3, .6, .9)
)

## Must match the TITE runner's model-specific truth settings exactly.
true_dose_specific_alpha <- c(0.2, 0.3, 0.4, 0.5, 0.6)
true_random_effect_eta <- 1
true_effective_dose_values <- c(15, 20, 30, 35, 45)
true_effective_dose_alpha_grid <- c(0.3, 0.6, 0.9)
true_generation_grid <- rbind(
  data.frame(
    true_generation_id = 1L, true_generation = "legacy", alpha_grid,
    true_random_effect_eta = true_random_effect_eta,
    true_random_effect_eta_tag = "na",
    true_effective_dose_alpha = 0,
    true_effective_dose_alpha_tag = "na",
    stringsAsFactors = FALSE
  ),
  data.frame(
    true_generation_id = 2L, true_generation = "dose_specific_geometric",
    toxicity_ipde_alpha = 0, efficacy_ipde_alpha = 0,
    true_random_effect_eta = true_random_effect_eta,
    true_random_effect_eta_tag = "na",
    true_effective_dose_alpha = 0,
    true_effective_dose_alpha_tag = "na",
    stringsAsFactors = FALSE
  ),
  data.frame(
    true_generation_id = 3L, true_generation = "shared_patient_logistic",
    toxicity_ipde_alpha = 0, efficacy_ipde_alpha = 0,
    true_random_effect_eta = true_random_effect_eta,
    true_random_effect_eta_tag = fmt_truth_setting(true_random_effect_eta),
    true_effective_dose_alpha = 0,
    true_effective_dose_alpha_tag = "na",
    stringsAsFactors = FALSE
  ),
  data.frame(
    true_generation_id = 4L, true_generation = "effective_dose_geometric",
    toxicity_ipde_alpha = 0, efficacy_ipde_alpha = 0,
    true_random_effect_eta = true_random_effect_eta,
    true_random_effect_eta_tag = "na",
    true_effective_dose_alpha = true_effective_dose_alpha_grid,
    true_effective_dose_alpha_tag = fmt_truth_setting(true_effective_dose_alpha_grid),
    stringsAsFactors = FALSE
  )
)

utility_scores <- c(u00 = 0, u01 = 40, u10 = 60, u11 = 100)
cohort_size <- 3L
cycle_max <- 2L
T_assess <- 28
n_eval <- 3L
m_U <- 6L
dlt_dist <- "uniform"
efficacy_dist <- "uniform"
crm_skeleton <- c(.15, .20, .30, .35, .45)
toxicity_elimination_cutoff <- .95
efficacy_threshold <- .20
futility_cutoff <- .85
min_eff_n_for_futility <- 0L
apply_ipde_toxicity_rule <- TRUE
ipde_toxicity_cutoff <- .95

## ============================================================
## Task reconstruction helpers
## ============================================================

fmt_short <- function(x, digits = 4L) {
  out <- format(round(as.numeric(x), digits), scientific = FALSE, trim = TRUE)
  out <- sub("(\\.[0-9]*?)0+$", "\\1", out)
  out <- sub("\\.$", "", out)
  gsub("\\.", "p", gsub("-", "m", out))
}
cross_join <- function(x, y) merge(x, y, by = NULL, sort = FALSE)
read_phase12_truth <- function(file) {
  if (!file.exists(file)) stop("Scenario summary does not exist: ", file)
  truth <- utils::read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("Scenario", "Source_Scenario", "Scenario_Group", "Attempt",
                "True_MTD_Level", "Target_Toxicity", paste0("Tox_Dose", 1:5),
                paste0("Eff_Dose", 1:5))
  missing <- setdiff(required, names(truth))
  if (length(missing)) stop("Scenario summary is missing: ", paste(missing, collapse = ", "))
  truth$Scenario <- as.integer(truth$Scenario)
  if (anyNA(truth$Scenario) || anyDuplicated(truth$Scenario))
    stop("Scenario must be unique, non-missing integers.")
  truth
}
calculate_true_obd <- function(p_true, e_true, true_mtd, utility_type,
                               lambda_T, utility_scores) {
  if (is.na(true_mtd) || true_mtd < 1L || true_mtd > length(p_true)) return(NA_integer_)
  score <- aide_compute_utility(e_true, p_true,
    list(type = utility_type, lambda_T = lambda_T, scores = utility_scores))
  candidates <- seq_len(as.integer(true_mtd))
  candidates[which.max(score[candidates])]
}
make_phase12_tite_config_tag <- function(task) {
  paste0(
    "P12TITE-a", if (task$allocation == "two_stage") "2s" else "1s",
    "-om", task$allocation_mode_id,
    "-N", task$Nmax, "-s1", task$N_s1,
    "-u", task$utility_type, "-l", fmt_short(task$lambda_T),
    "-cm", task$carryover_model, "-tm", task$crm_r_model,
    "-em", task$efficacy_model,
    "-rp", fmt_short(task$crm_a_r), "x", fmt_short(task$crm_b_r),
    "-ep", fmt_short(task$efficacy_a), "x", fmt_short(task$efficacy_b),
    "-cp", fmt_short(task$carry_a), "x", fmt_short(task$carry_b),
    "-eth", fmt_short(task$efficacy_threshold),
    "-fut", fmt_short(task$futility_cutoff),
    "-w", fmt_short(task$T_assess), "-c", task$C,
    "-cyc", task$cycle_max, "-ne", task$n_eval,
    "-rate", fmt_short(task$arrival_rate),
    "-ta", fmt_short(task$toxicity_ipde_alpha),
    "-ea", fmt_short(task$efficacy_ipde_alpha),
    "-gm", task$true_generation_id,
    "-ge", task$true_random_effect_eta_tag,
    "-ga", task$true_effective_dose_alpha_tag
  )
}

if (!dir.exists(results_root)) stop("Results directory does not exist: ", results_root)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
truth <- read_phase12_truth(file.path(project_root, scenario_file))
if (!is.null(scenario_id_list)) {
  scenario_id_list <- as.integer(scenario_id_list)
  missing_scenarios <- setdiff(scenario_id_list, truth$Scenario)
  if (length(missing_scenarios))
    stop("scenario_id_list is not in scenario_file: ", paste(missing_scenarios, collapse = ", "))
  truth <- truth[truth$Scenario %in% scenario_id_list, , drop = FALSE]
}
setting_grid <- Reduce(cross_join, list(
  design_size_grid, allocation_mode_grid, model_grid, crm_prior_grid,
  efficacy_prior_grid, utility_grid, arrival_grid, true_generation_grid
))
setting_grid$setting_id <- seq_len(nrow(setting_grid))

tasks <- list()
task_id <- 1L
for (setting_row in seq_len(nrow(setting_grid))) {
  setting <- setting_grid[setting_row, , drop = FALSE]
  for (scenario_row in seq_len(nrow(truth))) {
    p_true <- as.numeric(unlist(truth[scenario_row, paste0("Tox_Dose", 1:5), drop = FALSE], use.names = FALSE))
    e_true <- as.numeric(unlist(truth[scenario_row, paste0("Eff_Dose", 1:5), drop = FALSE], use.names = FALSE))
    tasks[[task_id]] <- list(
      task_id = task_id, setting_id = as.integer(setting$setting_id),
      Scenario = as.integer(truth$Scenario[scenario_row]),
      Source_Scenario = as.integer(truth$Source_Scenario[scenario_row]),
      Scenario_Group = as.character(truth$Scenario_Group[scenario_row]),
      Attempt = as.integer(truth$Attempt[scenario_row]),
      true_mtd = as.integer(truth$True_MTD_Level[scenario_row]),
      target = as.numeric(truth$Target_Toxicity[scenario_row]),
      p_true = p_true, e_true = e_true,
      allocation = as.character(setting$allocation),
      allocation_mode_id = as.integer(setting$allocation_mode_id),
      allocation_mode = as.character(setting$allocation_mode),
      Nmax = as.integer(setting$Nmax),
      N_s1 = as.integer(setting$N_s1), N_s2 = as.integer(setting$N_s2),
      model_id = as.character(setting$model_id),
      carryover_model = as.character(setting$carryover_model),
      crm_r_model = as.character(setting$crm_r_model),
      efficacy_model = as.character(setting$efficacy_model),
      crm_a_r = as.numeric(setting$crm_a_r), crm_b_r = as.numeric(setting$crm_b_r),
      efficacy_a = as.numeric(setting$efficacy_a), efficacy_b = as.numeric(setting$efficacy_b),
      carry_a = as.numeric(setting$carry_a), carry_b = as.numeric(setting$carry_b),
      utility_type = as.integer(setting$utility_type), lambda_T = as.numeric(setting$lambda_T),
      arrival_rate = as.numeric(setting$arrival_rate),
      toxicity_ipde_alpha = as.numeric(setting$toxicity_ipde_alpha),
      efficacy_ipde_alpha = as.numeric(setting$efficacy_ipde_alpha),
      true_generation_id = as.integer(setting$true_generation_id),
      true_generation = as.character(setting$true_generation),
      true_dose_specific_alpha = true_dose_specific_alpha,
      true_random_effect_eta = as.numeric(setting$true_random_effect_eta),
      true_random_effect_eta_tag = as.character(setting$true_random_effect_eta_tag),
      true_effective_dose_values = true_effective_dose_values,
      true_effective_dose_alpha = as.numeric(setting$true_effective_dose_alpha),
      true_effective_dose_alpha_tag = as.character(setting$true_effective_dose_alpha_tag),
      C = cohort_size, cycle_max = cycle_max, T_assess = T_assess,
      n_eval = n_eval, m_U = m_U, dlt_dist = dlt_dist, efficacy_dist = efficacy_dist,
      crm_skeleton = crm_skeleton, utility_scores = utility_scores,
      toxicity_elimination_cutoff = toxicity_elimination_cutoff,
      efficacy_threshold = efficacy_threshold, futility_cutoff = futility_cutoff,
      min_eff_n_for_futility = min_eff_n_for_futility,
      apply_ipde_toxicity_rule = apply_ipde_toxicity_rule,
      ipde_toxicity_cutoff = ipde_toxicity_cutoff
    )
    tasks[[task_id]]$true_obd <- calculate_true_obd(
      p_true, e_true, tasks[[task_id]]$true_mtd,
      tasks[[task_id]]$utility_type, tasks[[task_id]]$lambda_T, utility_scores
    )
    task_id <- task_id + 1L
  }
}

## ============================================================
## Extract and summarize available result files
## ============================================================

jobs_expected <- sort(unique(as.integer(jobs_expected)))
if (!length(jobs_expected) || anyNA(jobs_expected) || any(jobs_expected < 1L))
  stop("jobs_expected must contain positive integer IDX values.")
task_file <- function(task, job_i) {
  file.path(results_root, make_phase12_tite_config_tag(task),
            sprintf("P12TITE-SC%02d-j%04d.rds", task$Scenario, as.integer(job_i)))
}
trial_matrix <- function(results, field, ndose) {
  value <- do.call(rbind, lapply(results, `[[`, field))
  if (ncol(value) != ndose) stop("Unexpected number of dose columns in ", field, ".")
  value
}
safe_mean <- function(x) if (!length(x) || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
summarize_task <- function(results, task) {
  ndose <- length(task$p_true)
  ntrial <- sum(vapply(results, function(x) as.integer(x$ntrial), integer(1)))
  mtd <- unlist(lapply(results, `[[`, "MTD_by_trial"), use.names = FALSE)
  obd <- unlist(lapply(results, `[[`, "OBD_by_trial"), use.names = FALSE)
  ## Normalize legacy result files that encoded no OBD as NA.
  obd[is.na(obd)] <- aide_phase12_no_obd
  trial_summary <- do.call(rbind, lapply(results, `[[`, "trial_summary"))
  dose_n <- colMeans(trial_matrix(results, "n_by_dose_by_trial", ndose))
  unique_n <- colMeans(trial_matrix(results, "unique_n_by_dose_by_trial", ndose))
  ipde_n <- colMeans(trial_matrix(results, "nipde_by_dose_by_trial", ndose))
  tox <- colMeans(trial_matrix(results, "crm_pj_by_trial", ndose))
  eff <- colMeans(trial_matrix(results, "efficacy_pj_by_trial", ndose))
  utility <- colMeans(trial_matrix(results, "utility_by_trial", ndose))
  carry <- colMeans(trial_matrix(results, "r_hat_by_trial", ndose))
  mtd_pct <- 100 * tabulate(mtd[!is.na(mtd) & mtd >= 1L & mtd <= ndose], nbins = ndose) / ntrial
  obd_pct <- 100 * tabulate(obd[!is.na(obd) & obd >= 1L & obd <= ndose], nbins = ndose) / ntrial
  no_obd_pct <- 100 * mean(obd == aide_phase12_no_obd)
  data.frame(
    Task_ID = task$task_id, Scenario = task$Scenario,
    Source_Scenario = task$Source_Scenario, Scenario_Group = task$Scenario_Group,
    Attempt = task$Attempt, Allocation = task$allocation,
    Allocation_Mode_ID = task$allocation_mode_id,
    Allocation_Mode = task$allocation_mode,
    Stage2_Allocation = "one_stage", Model_ID = task$model_id,
    Carryover_Model = task$carryover_model, Efficacy_Model = task$efficacy_model,
    Nmax = task$Nmax, N_s1 = task$N_s1, N_s2 = task$N_s2,
    n_eval = task$n_eval, Cycle_Max = task$cycle_max,
    Arrival_Rate = task$arrival_rate, T_assess = task$T_assess,
    True_MTD = task$true_mtd, True_OBD = task$true_obd,
    Target_Toxicity = task$target, Utility_Type = task$utility_type,
    Lambda_T = task$lambda_T, Dose = seq_len(ndose),
    Toxicity_IPDE_Alpha = task$toxicity_ipde_alpha,
    Efficacy_IPDE_Alpha = task$efficacy_ipde_alpha,
    True_Generation_ID = task$true_generation_id,
    True_Generation = task$true_generation,
    True_Dose_Specific_Alpha = paste(task$true_dose_specific_alpha, collapse = ","),
    True_Random_Effect_Eta = task$true_random_effect_eta,
    True_Effective_Dose_Values = paste(task$true_effective_dose_values, collapse = ","),
    True_Effective_Dose_Alpha = task$true_effective_dose_alpha,
    True_DLT_rate = task$p_true, True_Efficacy_rate = task$e_true,
    Estimated_CRM_pj = tox, Estimated_Efficacy = eff,
    Estimated_Utility = utility, Carryover_Mean = carry,
    MTD_Selection_pct = mtd_pct, OBD_Selection_pct = obd_pct,
    No_OBD_Selection_pct = no_obd_pct,
    Pts_Treated = dose_n, Unique_Pts_by_Dose = unique_n, IPDE_Doses = ipde_n,
    Total_Administrations = safe_mean(trial_summary$administrations),
    Total_Unique_Patients = safe_mean(trial_summary$unique_patients),
    Duration = safe_mean(trial_summary$duration),
    Early_Stop_pct = 100 * safe_mean(trial_summary$early_stop),
    ntrial = ntrial, stringsAsFactors = FALSE
  )
}

summaries <- list()
missing_jobs <- list()
for (task in tasks) {
  paths <- vapply(jobs_expected, task_file, character(1), task = task)
  found <- paths[file.exists(paths)]
  missing <- jobs_expected[!file.exists(paths)]
  missing_jobs[[length(missing_jobs) + 1L]] <- data.frame(
    Task_ID = task$task_id, Scenario = task$Scenario, Model_ID = task$model_id,
    Folder = dirname(paths[1L]), n_found = length(found), n_expected = length(paths),
    missing_IDX = paste(missing, collapse = ","), stringsAsFactors = FALSE
  )
  if (length(found)) summaries[[length(summaries) + 1L]] <- summarize_task(lapply(found, readRDS), task)
}
if (!length(summaries)) stop("No TITE result files were found.")

dose_summary <- do.call(rbind, summaries)
missing_summary <- do.call(rbind, missing_jobs)
model3_alpha_tags <- unique(
  true_generation_grid$true_effective_dose_alpha_tag[
    true_generation_grid$true_generation_id == 4L
  ]
)
model2_eta_tags <- unique(
  true_generation_grid$true_random_effect_eta_tag[
    true_generation_grid$true_generation_id == 3L
  ]
)
tag <- paste0(
  "TITE_AIDE_phase_I_II",
  "_allocmode", paste(unique(allocation_mode_grid$allocation_mode_id), collapse = "_"),
  "_generation_model", paste(unique(true_generation_grid$true_generation_id), collapse = "_"),
  "_model3alpha", paste(model3_alpha_tags, collapse = "_"),
  "_model2eta", paste(model2_eta_tags, collapse = "_"),
  "_IDX_", sprintf("%04d", min(jobs_expected)),
              "_to_", sprintf("%04d", max(jobs_expected)))
utils::write.csv(dose_summary, file.path(out_dir, paste0(tag, "_dose_summary.csv")), row.names = FALSE)
utils::write.csv(missing_summary, file.path(out_dir, paste0(tag, "_missing_jobs.csv")), row.names = FALSE)
saveRDS(dose_summary, file.path(out_dir, paste0(tag, "_dose_summary.rds")))
saveRDS(missing_summary, file.path(out_dir, paste0(tag, "_missing_jobs.rds")))
cat("Saved TITE Phase I/II summaries to", out_dir, "\n")
