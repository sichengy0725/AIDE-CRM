## ============================================================
## Extract TITE-AIDE Phase I/II operating characteristics.
## Keep scenario and grid settings synchronized with the runner.
## ============================================================

source("TITE-AIDE.R")
settings <- aide_tite_oc_settings()
scenario_set <- tools::file_path_sans_ext(basename(settings$scenario_file))
results_root <- paste0("oc_results_cluster_TITE_AIDE_phase_I_II_", scenario_set)
out_dir <- "OC_summary_TITE_AIDE_phase_I_II"
jobs_expected <- 1:1000
ntrial_per_job <- 1L
if (!dir.exists(results_root)) stop("Results directory does not exist: ", results_root)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

truth <- aide_tite_read_truth(settings$scenario_file, settings$scenario_ids)
tasks <- aide_tite_make_tasks(settings, truth, 1L, ntrial_per_job, getwd())

task_file <- function(task, job) {
  task$job_i <- as.integer(job); task$ntrial <- ntrial_per_job; task$seed <- settings$seed_base + (job - 1L) * ntrial_per_job
  file.path(results_root, aide_tite_tag(task), sprintf("tite-sc%02d-j%04d.rds", task$Scenario, task$job_i))
}
trial_matrix <- function(results, field, ndose) do.call(rbind, lapply(results, function(x) as.matrix(x[[field]])))
safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
summarize_task <- function(results, task) {
  ndose <- length(task$p_true); ntrial <- sum(vapply(results, function(x) x$ntrial, integer(1)))
  mtd <- unlist(lapply(results, `[[`, "MTD_by_trial")); obd <- unlist(lapply(results, `[[`, "OBD_by_trial"))
  dose_n <- colMeans(trial_matrix(results, "n_by_dose_by_trial", ndose))
  unique_n <- colMeans(trial_matrix(results, "unique_n_by_dose_by_trial", ndose))
  ipde_n <- colMeans(trial_matrix(results, "nipde_by_dose_by_trial", ndose))
  tox <- colMeans(trial_matrix(results, "crm_pj_by_trial", ndose)); eff <- colMeans(trial_matrix(results, "efficacy_pj_by_trial", ndose))
  utility <- colMeans(trial_matrix(results, "utility_by_trial", ndose)); carry <- colMeans(trial_matrix(results, "r_hat_by_trial", ndose))
  summaries <- do.call(rbind, lapply(results, `[[`, "trial_summary"))
  mtd_pct <- 100 * tabulate(mtd[!is.na(mtd) & mtd >= 1L & mtd <= ndose], nbins = ndose) / ntrial
  obd_pct <- 100 * tabulate(obd[!is.na(obd) & obd >= 1L & obd <= ndose], nbins = ndose) / ntrial
  data.frame(Task_ID = task$task_id, Scenario = task$Scenario, Source_Scenario = task$Source_Scenario,
    Scenario_Group = task$Scenario_Group, Attempt = task$Attempt, Allocation = task$allocation, Model_ID = task$model_id,
    Toxicity_Model = task$toxicity_model, Efficacy_Model = task$efficacy_model, Nmax = task$Nmax, s1_Max = task$s1_Max,
    N_s2 = task$N_s2, n_eval = task$n_eval, Cycle_Max = task$cycle_max, Arrival_Rate = task$arrival_rate,
    True_MTD = task$true_mtd, Target_Toxicity = task$target, Dose = seq_len(ndose), True_DLT_rate = task$p_true,
    True_Efficacy_rate = task$e_true, Estimated_CRM_pj = tox, Estimated_Efficacy = eff, Estimated_Utility = utility,
    Carryover_Mean = carry, MTD_Selection_pct = mtd_pct,
    OBD_Selection_pct = obd_pct,
    Pts_Treated = dose_n, Unique_Pts_by_Dose = unique_n, IPDE_Doses = ipde_n,
    Total_Administrations = safe_mean(summaries$administrations), Total_Unique_Patients = safe_mean(vapply(results, function(x) mean(x$n_unique_patients_by_trial), numeric(1))),
    Duration = safe_mean(summaries$duration), n_eval_blocks = safe_mean(summaries$n_eval_blocks), ntrial = ntrial)
}

all_summaries <- list(); missing <- list(); k <- 0L
for (task in tasks) {
  paths <- vapply(jobs_expected, function(job) task_file(task, job), character(1)); found <- paths[file.exists(paths)]
  k <- k + 1L; missing[[k]] <- data.frame(Task_ID = task$task_id, Scenario = task$Scenario, Model_ID = task$model_id,
    Folder = dirname(paths[1L]), n_found = length(found), n_expected = length(paths), missing_IDX = paste(jobs_expected[!file.exists(paths)], collapse = ","))
  if (!length(found)) next
  all_summaries[[length(all_summaries) + 1L]] <- summarize_task(lapply(found, readRDS), task)
}
if (!length(all_summaries)) stop("No result files were found.")
dose_summary <- do.call(rbind, all_summaries); missing_summary <- do.call(rbind, missing)
tag <- paste0("TITE_AIDE_phase_I_II_IDX_", sprintf("%04d", min(jobs_expected)), "_to_", sprintf("%04d", max(jobs_expected)))
utils::write.csv(dose_summary, file.path(out_dir, paste0(tag, "_dose_summary.csv")), row.names = FALSE)
utils::write.csv(missing_summary, file.path(out_dir, paste0(tag, "_missing_jobs.csv")), row.names = FALSE)
saveRDS(dose_summary, file.path(out_dir, paste0(tag, "_dose_summary.rds")))
saveRDS(missing_summary, file.path(out_dir, paste0(tag, "_missing_jobs.rds")))
cat("Saved summaries to", out_dir, "\n")
