## ============================================================
## Extract operating characteristics from run_oc_AIDE_phase_I_II.R
##
## Usage:
##   Rscript extract_result_AIDE_phase_I_II.R [first_IDX] [last_IDX]
##
## Examples:
##   Rscript extract_result_AIDE_phase_I_II.R 1 1000
##   Rscript extract_result_AIDE_phase_I_II.R 24 24
##
## By default this reads oc_results_AIDE_phase_I_II/ and writes a dose-level
## summary, a table-style summary, and a missing-job log to
## OC_summary_AIDE_phase_I_II/.  Set AIDE_PHASE12_RESULTS_DIR or
## AIDE_PHASE12_SUMMARY_DIR to override either directory.
## ============================================================

args <- commandArgs(trailingOnly = TRUE)
script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
if (length(script_file) == 1L && file.exists(script_file)) {
  setwd(dirname(normalizePath(script_file)))
}

parse_positive_integer <- function(x, name) {
  value <- suppressWarnings(as.integer(x))
  if (length(value) != 1L || is.na(value) || value < 1L) {
    stop(name, " must be a positive integer.")
  }
  value
}

first_idx <- if (length(args) >= 1L) {
  parse_positive_integer(args[1L], "first_IDX")
} else {
  1L
}
last_idx <- if (length(args) >= 2L) {
  parse_positive_integer(args[2L], "last_IDX")
} else {
  1000L
}
if (last_idx < first_idx) stop("last_IDX must be at least first_IDX.")

results_root <- Sys.getenv(
  "AIDE_PHASE12_RESULTS_DIR",
  unset = "oc_results_AIDE_phase_I_II"
)
out_dir <- Sys.getenv(
  "AIDE_PHASE12_SUMMARY_DIR",
  unset = "OC_summary_AIDE_phase_I_II"
)
if (!dir.exists(results_root)) stop("Results directory does not exist: ", results_root)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

task_id_from_file <- function(path) {
  as.integer(sub("^task_([0-9]+)_.*$", "\\1", basename(path)))
}

job_id_from_file <- function(path) {
  as.integer(sub("^.*_job_([0-9]+)\\.rds$", "\\1", basename(path)))
}

task_value <- function(task, name, default = NA) {
  value <- task[[name]]
  if (is.null(value) || length(value) == 0L) return(default)
  value[[1L]]
}

result_truth <- function(result) {
  truth <- result$truth
  if (is.null(truth)) truth <- list()
  list(
    p_true = if (!is.null(truth$p_true)) truth$p_true else result$p_true,
    e_true = if (!is.null(truth$e_true)) truth$e_true else result$e_true,
    target = task_value(truth, "target"),
    true_mtd = task_value(truth, "true_mtd"),
    true_obd = task_value(truth, "true_obd")
  )
}

safe_col_mean <- function(mat) {
  out <- colMeans(mat, na.rm = TRUE)
  out[is.nan(out)] <- NA_real_
  out
}

as_trial_matrix <- function(value, ndose, ntrial, label) {
  if (is.null(dim(value))) {
    value <- as.numeric(value)
    if (length(value) == ndose) {
      value <- matrix(value, nrow = 1L, ncol = ndose)
    } else if (length(value) == ntrial * ndose) {
      value <- matrix(value, nrow = ntrial, ncol = ndose, byrow = TRUE)
    } else {
      stop(label, " has an unexpected length.")
    }
  } else {
    value <- as.matrix(value)
  }
  if (ncol(value) != ndose || nrow(value) != ntrial) {
    stop(label, " has dimensions inconsistent with ntrial and the number of doses.")
  }
  storage.mode(value) <- "double"
  value
}

extract_trial_matrix <- function(results, by_trial_field, mean_field, ndose) {
  rows <- lapply(results, function(result) {
    ntrial <- as.integer(result$ntrial)
    value <- result[[by_trial_field]]
    if (!is.null(value)) {
      return(as_trial_matrix(value, ndose, ntrial, by_trial_field))
    }
    value <- result[[mean_field]]
    if (!is.null(value) && length(value) == ndose) {
      return(matrix(rep(as.numeric(value), ntrial),
                    nrow = ntrial, ncol = ndose, byrow = TRUE))
    }
    matrix(NA_real_, nrow = ntrial, ncol = ndose)
  })
  do.call(rbind, rows)
}

extract_trial_vector <- function(results,
                                 by_trial_field,
                                 mean_field = NULL,
                                 mean_multiplier = 1,
                                 required = FALSE) {
  out <- lapply(results, function(result) {
    ntrial <- as.integer(result$ntrial)
    value <- result[[by_trial_field]]
    if (!is.null(value)) {
      value <- as.numeric(value)
      if (length(value) != ntrial) {
        stop(by_trial_field, " has length inconsistent with ntrial.")
      }
      return(value)
    }
    if (!is.null(mean_field) && !is.null(result[[mean_field]])) {
      return(rep(as.numeric(result[[mean_field]]) * mean_multiplier, ntrial))
    }
    if (required) stop("Missing required result field: ", by_trial_field)
    rep(NA_real_, ntrial)
  })
  unlist(out, use.names = FALSE)
}

selection_rate <- function(selection, ndose, denominator) {
  keep <- !is.na(selection) & selection >= 1L & selection <= ndose
  100 * tabulate(as.integer(selection[keep]), nbins = ndose) / denominator
}

make_metadata <- function(result, task_id) {
  task <- result$task
  if (is.null(task)) task <- list()
  truth <- result_truth(result)
  runner <- result$runner_settings
  if (is.null(runner)) runner <- list()

  data.frame(
    Task_ID = as.integer(task_id),
    Scenario = as.integer(task_value(task, "Scenario")),
    Source_Scenario = as.integer(task_value(task, "Source_Scenario")),
    Scenario_Group = as.character(task_value(task, "Scenario_Group")),
    Scenario_Attempt = as.integer(task_value(task, "Attempt")),
    True_MTD = as.integer(truth$true_mtd),
    True_OBD = as.integer(truth$true_obd),
    Target_Toxicity = as.numeric(truth$target),
    Allocation = as.character(task_value(task, "allocation")),
    Nmax = as.integer(task_value(task, "Nmax")),
    N_s1 = as.integer(task_value(task, "N_s1")),
    Model = "CRM",
    CRM_r_Model = "random",
    CRM_Prior_a = as.numeric(task_value(task, "crm_a_r")),
    CRM_Prior_b = as.numeric(task_value(task, "crm_b_r")),
    Efficacy_Prior_a = as.numeric(task_value(task, "efficacy_a")),
    Efficacy_Prior_b = as.numeric(task_value(task, "efficacy_b")),
    Efficacy_Carryover_Prior_a = as.numeric(task_value(task, "carry_a")),
    Efficacy_Carryover_Prior_b = as.numeric(task_value(task, "carry_b")),
    Utility_Type = as.integer(task_value(task, "utility_type")),
    Lambda_T = as.numeric(task_value(task, "lambda_T")),
    m_U = as.integer(task_value(runner, "m_U", 6L)),
    Enrollment_Scheme = as.character(task_value(task, "enrollment_scheme")),
    Arrival_Rate = as.numeric(task_value(task, "arrival_rate")),
    T_assess = as.numeric(task_value(runner, "T_assess", 28)),
    Cycle_Max = as.integer(task_value(runner, "cycle_max", 1L)),
    IPDE_Design = as.integer(task_value(task, "ipde_design")),
    Flexible_IPDE = as.integer(isTRUE(task_value(task, "flexible_ipde", FALSE))),
    Toxicity_IPDE_Alpha = as.numeric(task_value(task, "toxicity_ipde_alpha")),
    Efficacy_IPDE_Alpha = as.numeric(task_value(task, "efficacy_ipde_alpha")),
    stringsAsFactors = FALSE
  )
}

summarize_task <- function(results, task_id) {
  first <- results[[1L]]
  truth <- result_truth(first)
  p_true <- as.numeric(truth$p_true)
  e_true <- as.numeric(truth$e_true)
  if (length(p_true) < 1L || length(p_true) != length(e_true)) {
    stop("Task ", task_id, " does not contain matching toxicity and efficacy truth vectors.")
  }
  ndose <- length(p_true)
  ntrial <- sum(vapply(results, function(x) as.integer(x$ntrial), integer(1)))

  mtd <- extract_trial_vector(results, "MTD_by_trial", required = TRUE)
  obd <- extract_trial_vector(results, "OBD_by_trial", required = TRUE)
  design_early_stop <- extract_trial_vector(
    results, "early_stop_by_trial", "early_stop_percent", mean_multiplier = 0.01
  )
  n_admin <- extract_trial_vector(
    results, "n_administrations_by_trial", "mean_administrations"
  )
  n_unique <- extract_trial_vector(
    results, "n_unique_patients_by_trial", "mean_unique_patients"
  )
  duration <- extract_trial_vector(results, "duration_by_trial", "mean_duration")

  n_by_dose <- safe_col_mean(extract_trial_matrix(
    results, "n_by_dose_by_trial", "mean_administrations_by_dose", ndose
  ))
  unique_n_by_dose <- safe_col_mean(extract_trial_matrix(
    results, "unique_n_by_dose_by_trial", "mean_unique_patients_by_dose", ndose
  ))
  nipde_by_dose <- safe_col_mean(extract_trial_matrix(
    results, "nipde_by_dose_by_trial", "mean_ipde_doses_by_dose", ndose
  ))
  crm_pj <- safe_col_mean(extract_trial_matrix(
    results, "crm_pj_by_trial", "crm_pj_mean", ndose
  ))
  efficacy_pj <- safe_col_mean(extract_trial_matrix(
    results, "efficacy_pj_by_trial", "efficacy_pj_mean", ndose
  ))
  utility <- safe_col_mean(extract_trial_matrix(
    results, "utility_by_trial", "utility_mean", ndose
  ))
  r_hat <- safe_col_mean(extract_trial_matrix(
    results, "r_hat_by_trial", "r_hat_mean", ndose
  ))

  total_admin <- mean(n_admin, na.rm = TRUE)
  if (is.nan(total_admin)) total_admin <- NA_real_
  total_unique <- mean(n_unique, na.rm = TRUE)
  if (is.nan(total_unique)) total_unique <- NA_real_
  mean_duration <- mean(duration, na.rm = TRUE)
  if (is.nan(mean_duration)) mean_duration <- NA_real_

  mtd_selection <- selection_rate(mtd, ndose, ntrial)
  obd_selection <- selection_rate(obd, ndose, ntrial)
  early_stopping <- 100 * mean(mtd == 99L, na.rm = TRUE)
  if (is.nan(early_stopping)) early_stopping <- NA_real_
  design_early_stopping <- 100 * mean(design_early_stop, na.rm = TRUE)
  if (is.nan(design_early_stopping)) design_early_stopping <- NA_real_
  no_obd <- 100 * mean(is.na(obd) | obd == 99L)

  metadata <- make_metadata(first, task_id)
  dose_summary <- metadata[rep(1L, ndose), , drop = FALSE]
  dose_summary$Dose <- seq_len(ndose)
  dose_summary$True_DLT_rate <- p_true
  dose_summary$True_Efficacy_rate <- e_true
  dose_summary$Estimated_CRM_pj <- crm_pj
  dose_summary$Estimated_Efficacy <- efficacy_pj
  dose_summary$Estimated_Utility <- utility
  dose_summary$r_hat <- r_hat
  dose_summary$r_cap <- NA_real_
  dose_summary$r_use <- NA_real_
  dose_summary$MTD_Selection_pct <- mtd_selection
  dose_summary$OBD_Selection_pct <- obd_selection
  dose_summary$Pts_Treated <- n_by_dose
  dose_summary$Unique_Pts_by_Dose <- unique_n_by_dose
  dose_summary$IPDE_Doses <- nipde_by_dose
  dose_summary$Total_Administrations <- total_admin
  dose_summary$Total_Unique_Patients <- total_unique
  dose_summary$Early_Stopping_pct <- early_stopping
  dose_summary$Design_Early_Stopping_pct <- design_early_stopping
  dose_summary$No_OBD_Selection_pct <- no_obd
  dose_summary$Duration <- mean_duration
  dose_summary$n_valid <- ntrial
  dose_summary$ntrial_from_files <- ntrial

  metrics <- c(
    "True DLT rate",
    "True efficacy rate",
    "Estimated CRM pj",
    "Estimated efficacy",
    "Estimated utility",
    "r_hat",
    "r_cap",
    "r_use",
    "MTD Selection %",
    "OBD Selection %",
    "# Pts Treated",
    "# Unique Pts by Dose",
    "# IPDE Doses",
    "# Total Administrations",
    "# Total Unique Patients",
    "% Early Stopping",
    "% Design Early Stopping",
    "% No OBD Selection",
    "Duration"
  )
  table_summary <- metadata[rep(1L, length(metrics)), , drop = FALSE]
  table_summary$Metric <- metrics
  dose_cols <- paste0("D", seq_len(ndose))
  for (dose_col in dose_cols) table_summary[[dose_col]] <- NA_real_
  table_summary$Total <- NA_real_
  table_summary$Duration <- NA_real_
  table_summary$n_valid <- ntrial
  table_summary$ntrial_from_files <- ntrial

  table_summary[1L, dose_cols] <- p_true
  table_summary[2L, dose_cols] <- e_true
  table_summary[3L, dose_cols] <- crm_pj
  table_summary[4L, dose_cols] <- efficacy_pj
  table_summary[5L, dose_cols] <- utility
  table_summary[6L, dose_cols] <- r_hat
  table_summary[9L, dose_cols] <- mtd_selection
  table_summary[10L, dose_cols] <- obd_selection
  table_summary[11L, dose_cols] <- n_by_dose
  table_summary[12L, dose_cols] <- unique_n_by_dose
  table_summary[13L, dose_cols] <- nipde_by_dose
  table_summary[14L, "Total"] <- total_admin
  table_summary[15L, "Total"] <- total_unique
  table_summary[16L, "Total"] <- early_stopping
  table_summary[17L, "Total"] <- design_early_stopping
  table_summary[18L, "Total"] <- no_obd
  table_summary[19L, "Duration"] <- mean_duration

  list(dose_summary = dose_summary, table_summary = table_summary, ntrial = ntrial)
}

raw_files <- list.files(
  results_root,
  pattern = "^task_[0-9]+_scenario_[0-9]+_.*_job_[0-9]+\\.rds$",
  full.names = TRUE
)
if (length(raw_files) == 0L) stop("No Phase I/II task RDS files found in: ", results_root)

job_ids <- vapply(raw_files, job_id_from_file, integer(1))
task_ids <- vapply(raw_files, task_id_from_file, integer(1))
keep <- job_ids >= first_idx & job_ids <= last_idx
raw_files <- raw_files[keep]
job_ids <- job_ids[keep]
task_ids <- task_ids[keep]
if (length(raw_files) == 0L) {
  stop("No RDS files found for IDX ", first_idx, " through ", last_idx, ".")
}

manifest_file <- file.path(results_root, "AIDE_phase_I_II_task_manifest.csv")
manifest_task_ids <- integer(0)
if (file.exists(manifest_file)) {
  manifest <- utils::read.csv(manifest_file, stringsAsFactors = FALSE)
  if ("task_id" %in% names(manifest)) manifest_task_ids <- as.integer(manifest$task_id)
}
all_task_ids <- sort(unique(c(task_ids, manifest_task_ids)))

dose_summaries <- list()
table_summaries <- list()
missing_log <- list()
expected_jobs <- seq.int(first_idx, last_idx)

for (task_id in all_task_ids) {
  these <- which(task_ids == task_id)
  files_for_task <- raw_files[these]
  jobs_for_task <- job_ids[these]
  missing_jobs <- setdiff(expected_jobs, jobs_for_task)
  ntrial_from_files <- 0L
  if (length(files_for_task) > 0L) {
    results <- lapply(files_for_task, readRDS)
    ntrial_from_files <- sum(vapply(results, function(x) as.integer(x$ntrial), integer(1)))
    one <- summarize_task(results, task_id)
    dose_summaries[[as.character(task_id)]] <- one$dose_summary
    table_summaries[[as.character(task_id)]] <- one$table_summary
  }
  missing_log[[as.character(task_id)]] <- data.frame(
    Task_ID = task_id,
    n_found = length(files_for_task),
    n_expected = length(expected_jobs),
    n_missing = length(missing_jobs),
    ntrial_from_files = ntrial_from_files,
    missing_IDX = paste(missing_jobs, collapse = ","),
    stringsAsFactors = FALSE
  )
}

if (length(dose_summaries) == 0L) stop("No readable Phase I/II result files were found.")

round_numeric_df <- function(dat, digits = 4L) {
  is_num <- vapply(dat, is.numeric, logical(1))
  dat[is_num] <- lapply(dat[is_num], round, digits = digits)
  dat
}

dose_summary <- do.call(rbind, dose_summaries)
table_summary <- do.call(rbind, table_summaries)
missing_summary <- do.call(rbind, missing_log)
tag <- paste0("AIDE_phase_I_II_IDX_", sprintf("%04d", first_idx), "_to_", sprintf("%04d", last_idx))

dose_csv <- file.path(out_dir, paste0(tag, "_dose_summary.csv"))
table_csv <- file.path(out_dir, paste0(tag, "_table_summary.csv"))
missing_csv <- file.path(out_dir, paste0(tag, "_missing_jobs.csv"))
dose_rds <- file.path(out_dir, paste0(tag, "_dose_summary.rds"))
table_rds <- file.path(out_dir, paste0(tag, "_table_summary.rds"))
missing_rds <- file.path(out_dir, paste0(tag, "_missing_jobs.rds"))

utils::write.csv(round_numeric_df(dose_summary), dose_csv, row.names = FALSE)
utils::write.csv(round_numeric_df(table_summary), table_csv, row.names = FALSE)
utils::write.csv(missing_summary, missing_csv, row.names = FALSE)
saveRDS(dose_summary, dose_rds)
saveRDS(table_summary, table_rds)
saveRDS(missing_summary, missing_rds)

cat("Saved dose-level summary:", dose_csv, "\n")
cat("Saved table-style summary:", table_csv, "\n")
cat("Saved missing-job log:", missing_csv, "\n")
