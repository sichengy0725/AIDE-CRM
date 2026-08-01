## ============================================================
## Extract operating characteristics from run_oc_AIDE_phase_I_II.R
##
## Edit the settings below, then run this file interactively (for example,
## with RStudio's Source button).  Keep the scenario and design settings in
## sync with run_oc_AIDE_phase_I_II.R.  The extractor rebuilds the same task
## grid so that it can report missing result files for every requested task.
## ============================================================

## ------------------------------------------------------------
## User settings: match these to run_oc_AIDE_phase_I_II.R.
## ------------------------------------------------------------

scenario_file <- "Set_5dose_adaptive_r_37_true_MTD_OBD_summary_lambda1.csv"
scenario_set_name <- tools::file_path_sans_ext(basename(scenario_file))
results_root <- paste0("oc_results_cluster_AIDE_phase_I_II_", scenario_set_name)
out_dir <- "OC_summary_AIDE_phase_I_II"

## Match the one-trial-per-IDX LSF array that produced the current files.
## The file names shown are job-1001 through job-2000, giving 1,000 trials
## per task. Change this range if a different IDX array was submitted.
jobs_expected <- 1001:2000
ntrial_per_job_expected <- 1L
seed_base_expected <- 1L
block_id_expected <- 1L

## The runner uses all 37 scenarios in the supplied truth file.  Keep this
## list synchronized if you run a subset of scenarios.
scenario_id_list <- 1:37

## For two-stage allocation, N_s1 and N_s2 are per-dose administration
## thresholds: Stage I ends at N_s1 only after a stay decision at the current
## dose, and Stage II ends at N_s2. Nmax remains the total administration
## limit. One-stage does not use either threshold.
two_stage_sizes <- data.frame(
  allocation = "two_stage",
  Nmax = c(30L, 60L),
  N_s1 = c(15L, 30L),
  N_s2 = c(30L, 60L),
  stringsAsFactors = FALSE
)
one_stage_sizes <- data.frame(
  allocation = "one_stage",
  Nmax = c(30L, 60L),
  N_s1 = c(30L, 60L),
  N_s2 = c(30L, 60L),
  stringsAsFactors = FALSE
)
design_size_grid <- rbind(two_stage_sizes, one_stage_sizes)

## Keep this paired model grid synchronized with the OC runner.
model_grid <- data.frame(
  model_id = c("random_carryover", "previous_dose_additive"),
  crm_r_model = c("random", "previous_dose"),
  efficacy_model = c("dose_specific_carryover", "previous_dose_additive"),
  stringsAsFactors = FALSE
)

## Beta-equivalent prior shapes for random CRM r or the additive toxicity
## alpha. The JAGS models implement these through equal-rate Gamma ratios.
crm_prior_grid <- data.frame(
  crm_prior_id = "r_beta_0p15_0p85",
  crm_a_r = 0.15,
  crm_b_r = 0.85,
  stringsAsFactors = FALSE
)

## Beta-equivalent prior shapes for regular efficacy, dose-specific carryover,
## and the shared additive efficacy alpha. The JAGS models implement these
## through equal-rate Gamma ratios.
efficacy_prior_grid <- data.frame(
  efficacy_prior_id = "regular_beta_0p15_0p85_carry_beta_0p15_0p85",
  efficacy_a = 0.15,
  efficacy_b = 0.85,
  carry_a = 0.15,
  carry_b = 0.85,
  efficacy_additive_alpha_a = 0.15,
  efficacy_additive_alpha_b = 0.85,
  stringsAsFactors = FALSE
)

enrollment_schemes <- data.frame(
  enrollment_scheme = c("continuous", "ipde_first"),
  stringsAsFactors = FALSE
)
lambda_T_grid <- 0.3
utility_grid <- rbind(
  # data.frame(utility_type = 1L, lambda_T = 1, stringsAsFactors = FALSE),
  data.frame(utility_type = 2L, lambda_T = lambda_T_grid, stringsAsFactors = FALSE),
  data.frame(utility_type = 3L, lambda_T = 1, stringsAsFactors = FALSE)
)
arrival_grid <- data.frame(arrival_rate = 1 / 56, stringsAsFactors = FALSE)
ipde_grid <- data.frame(
  ipde_design = 1L,
  flexible_ipde = TRUE,
  stringsAsFactors = FALSE
)
alpha_grid <- data.frame(
  toxicity_ipde_alpha = 0,
  efficacy_ipde_alpha = 0,
  stringsAsFactors = FALSE
)

## Fixed design settings, retained here as a compact record of the matching
## run configuration.  They do not change the task IDs, but should also match
## the corresponding values in run_oc_AIDE_phase_I_II.R.
cohort_size <- 3L
cycle_max <- 1L
T_assess <- 28
m_U <- 6L
crm_skeleton <- c(0.15, 0.20, 0.30, 0.35, 0.45)
crm_alpha_sd <- sqrt(2)
efficacy_threshold <- 0.20
futility_cutoff <- 0.95
min_eff_n_for_futility <- 0L
apply_random_crm_recycle_toxicity_rule <- FALSE
apply_random_crm_recycle_efficacy_rule <- FALSE

if (!dir.exists(results_root)) stop("Results directory does not exist: ", results_root)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cross_join <- function(x, y) {
  merge(x, y, by = NULL, sort = FALSE)
}

read_phase12_truth <- function(file) {
  if (!file.exists(file)) stop("Scenario summary does not exist: ", file)
  raw <- utils::read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)

  tox_cols <- paste0("Tox_Dose", 1:5)
  eff_cols <- paste0("Eff_Dose", 1:5)
  required <- c(
    "Scenario", "Source_Scenario", "Scenario_Group", "True_MTD_Level",
    "Target_Toxicity", tox_cols, eff_cols
  )
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0L) {
    stop("Scenario summary is missing: ", paste(missing, collapse = ", "))
  }

  truth_columns <- c(
    "Scenario", "Source_Scenario", "Scenario_Group", "Attempt",
    "True_MTD_Level", "Target_Toxicity", tox_cols, eff_cols
  )
  truth <- raw[, truth_columns[truth_columns %in% names(raw)], drop = FALSE]
  truth$Scenario <- as.integer(truth$Scenario)
  if (anyNA(truth$Scenario) || anyDuplicated(truth$Scenario)) {
    stop("Scenario must contain unique non-missing integer values.")
  }
  truth
}

truth <- read_phase12_truth(scenario_file)
if (!is.null(scenario_id_list)) {
  scenario_id_list <- as.integer(scenario_id_list)
  missing_scenarios <- setdiff(scenario_id_list, truth$Scenario)
  if (length(missing_scenarios) > 0L) {
    stop("scenario_id_list is not in scenario_file: ",
         paste(missing_scenarios, collapse = ", "))
  }
  truth <- truth[truth$Scenario %in% scenario_id_list, , drop = FALSE]
}

setting_grid <- Reduce(
  cross_join,
  list(
    design_size_grid,
    model_grid,
    crm_prior_grid,
    efficacy_prior_grid,
    enrollment_schemes,
    utility_grid,
    arrival_grid,
    ipde_grid,
    alpha_grid
  )
)
setting_grid$setting_id <- seq_len(nrow(setting_grid))
expected_tasks <- cross_join(
  truth[, c("Scenario", "Source_Scenario", "Scenario_Group", "Attempt"), drop = FALSE],
  setting_grid
)
expected_tasks$task_id <- seq_len(nrow(expected_tasks))
expected_tasks <- expected_tasks[, c("task_id", setdiff(names(expected_tasks), "task_id")), drop = FALSE]

jobs_expected <- sort(unique(as.integer(jobs_expected)))
if (length(jobs_expected) == 0L || anyNA(jobs_expected) || any(jobs_expected < 1L)) {
  stop("jobs_expected must contain one or more positive integer IDX values.")
}
ntrial_per_job_expected <- as.integer(ntrial_per_job_expected)
if (length(ntrial_per_job_expected) != 1L || is.na(ntrial_per_job_expected) ||
    ntrial_per_job_expected < 1L) {
  stop("ntrial_per_job_expected must be a positive integer.")
}
ntrial_expected <- length(jobs_expected) * ntrial_per_job_expected

fmt_short <- function(x, digits = 2L) {
  x <- round(as.numeric(x), digits)
  out <- format(x, scientific = FALSE, trim = TRUE)
  out <- sub("(\\.[0-9]*?)0+$", "\\1", out)
  out <- sub("\\.$", "", out)
  out[out == "-0"] <- "0"
  out[out == ""] <- "0"
  out <- gsub("-", "m", out)
  gsub("\\.", "p", out)
}

fmt_param <- function(x, digits = 4L) {
  fmt_short(x, digits = digits)
}

make_phase12_config_tag <- function(task) {
  allocation_tag <- if (identical(task$allocation, "two_stage")) "2s" else "1s"
  enrollment_tag <- if (identical(task$enrollment_scheme, "continuous")) {
    "cont"
  } else {
    "ipdefirst"
  }
  paste0(
    "P12",
    "-a", allocation_tag,
    "-N", fmt_short(task$Nmax),
    "-s1", fmt_short(task$N_s1),
    "-s2", fmt_short(task$N_s2),
    "-u", task$utility_type,
    "-l", fmt_short(task$lambda_T),
    "-en", enrollment_tag,
    "-tm", task$crm_r_model,
    "-em", task$efficacy_model,
    "-rp", fmt_param(task$crm_a_r), "x", fmt_param(task$crm_b_r),
    "-ep", fmt_param(task$efficacy_a), "x", fmt_param(task$efficacy_b),
    "-cp", fmt_param(task$carry_a), "x", fmt_param(task$carry_b),
    "-ap", fmt_param(task$efficacy_additive_alpha_a), "x",
    fmt_param(task$efficacy_additive_alpha_b),
    "-w", fmt_short(T_assess),
    "-c", fmt_short(cohort_size),
    "-cyc", fmt_short(cycle_max),
    "-rate", fmt_short(task$arrival_rate),
    "-ip", task$ipde_design,
    "-fx", as.integer(isTRUE(task$flexible_ipde)),
    "-ta", fmt_short(task$toxicity_ipde_alpha),
    "-ea", fmt_short(task$efficacy_ipde_alpha),
    "-tg", as.integer(isTRUE(apply_random_crm_recycle_toxicity_rule)),
    "-eg", as.integer(isTRUE(apply_random_crm_recycle_efficacy_rule))
  )
}

make_phase12_folder <- function(task) make_phase12_config_tag(task)

make_phase12_raw_rds_filenames <- function(task, jobs = jobs_expected) {
  jobs <- as.integer(jobs)
  seed_block <- as.integer(seed_base_expected) +
    (jobs - 1L) * as.integer(ntrial_per_job_expected)
  scenario_name <- paste0("P12-SC", as.integer(task$Scenario))
  paste0(
    scenario_name,
    "-CRM-phase12_", task$model_id,
    "-j", jobs,
    "-b", as.integer(block_id_expected),
    "-s", seed_block,
    "-n", as.integer(ntrial_per_job_expected),
    ".rds"
  )
}

make_phase12_combined_rds_filenames <- function(task, jobs = jobs_expected) {
  paste0(
    make_phase12_group_key(task),
    "-job-", as.integer(jobs),
    "-combined.rds"
  )
}

find_phase12_result_files <- function(task, jobs = jobs_expected) {
  folderpath <- file.path(results_root, make_phase12_folder(task))
  raw_paths <- file.path(
    folderpath,
    make_phase12_raw_rds_filenames(task, jobs = jobs)
  )
  combined_names <- make_phase12_combined_rds_filenames(task, jobs = jobs)

  ## Prefer the raw task file.  The combined file contains the same summary
  ## fields and is accepted as a fallback, including the older flat layout.
  combined_in_folder <- file.path(folderpath, combined_names)
  combined_in_root <- file.path(results_root, combined_names)
  raw_exists <- file.exists(raw_paths)
  combined_folder_exists <- file.exists(combined_in_folder)
  combined_root_exists <- file.exists(combined_in_root)

  selected_paths <- raw_paths
  selected_type <- rep("raw", length(jobs))
  use_combined_folder <- !raw_exists & combined_folder_exists
  use_combined_root <- !raw_exists & !combined_folder_exists & combined_root_exists
  selected_paths[use_combined_folder] <- combined_in_folder[use_combined_folder]
  selected_type[use_combined_folder] <- "combined (setting folder)"
  selected_paths[use_combined_root] <- combined_in_root[use_combined_root]
  selected_type[use_combined_root] <- "combined (results root)"

  found <- raw_exists | combined_folder_exists | combined_root_exists
  list(
    folderpath = folderpath,
    raw_paths = raw_paths,
    result_paths = selected_paths,
    result_type = selected_type,
    found = found
  )
}

make_phase12_group_key <- function(task) {
  paste0("P12-SC", task$Scenario, "_", make_phase12_config_tag(task))
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
    Sample_Size = as.integer(task_value(task, "Nmax")),
    N_s1 = as.integer(task_value(task, "N_s1")),
    N_s2 = as.integer(task_value(task, "N_s2", task_value(task, "Nmax"))),
    Model = "CRM",
    Model_ID = as.character(task_value(
      task, "model_id", task_value(runner, "model_id")
    )),
    CRM_r_Model = as.character(task_value(task, "crm_r_model", "random")),
    CRM_Prior_a = as.numeric(task_value(task, "crm_a_r")),
    CRM_Prior_b = as.numeric(task_value(task, "crm_b_r")),
    Efficacy_Model = as.character(task_value(
      task, "efficacy_model", "dose_specific_carryover"
    )),
    Efficacy_Prior_a = as.numeric(task_value(task, "efficacy_a")),
    Efficacy_Prior_b = as.numeric(task_value(task, "efficacy_b")),
    Efficacy_Carryover_Prior_a = as.numeric(task_value(task, "carry_a")),
    Efficacy_Carryover_Prior_b = as.numeric(task_value(task, "carry_b")),
    Efficacy_Additive_Alpha_Prior_a = as.numeric(task_value(
      task, "efficacy_additive_alpha_a"
    )),
    Efficacy_Additive_Alpha_Prior_b = as.numeric(task_value(
      task, "efficacy_additive_alpha_b"
    )),
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

## ------------------------------------------------------------
## Main extraction: explicit expected-file discovery, matching
## extract_result_AIDE.R. Each task accepts its raw RDS file, or its
## per-job combined RDS file when only combined outputs were retained.
## ------------------------------------------------------------

all.dose.summary <- list()
all.table.summary <- list()
missing.log <- list()
miss.idx <- 1L

for (task_row in seq_len(nrow(expected_tasks))) {
  task <- expected_tasks[task_row, , drop = FALSE]
  discovery <- find_phase12_result_files(task, jobs = jobs_expected)
  folderpath <- discovery$folderpath
  files.use <- discovery$result_paths[discovery$found]
  missing.jobs <- jobs_expected[!discovery$found]
  example_file <- if (length(files.use) > 0L) {
    files.use[1L]
  } else {
    discovery$raw_paths[1L]
  }

  cat("\n====================================\n")
  cat("Task:", task$task_id, "of", nrow(expected_tasks), "\n")
  cat("Scenario:", task$Scenario, "\n")
  cat("Source scenario:", task$Source_Scenario, "\n")
  cat("Scenario group:", task$Scenario_Group, "\n")
  cat("Scenario attempt:", task$Attempt, "\n")
  cat("Model: CRM\n")
  cat("CRM r model:", task$crm_r_model, "\n")
  cat("Efficacy model:", task$efficacy_model, "\n")
  cat("Allocation:", task$allocation, "\n")
  cat("Nmax:", task$Nmax, "\n")
  cat("N_s1:", task$N_s1, "\n")
  if ("N_s2" %in% names(task)) cat("N_s2:", task$N_s2, "\n")
  cat("CRM prior:", task$crm_prior_id, "\n")
  cat("Efficacy prior:", task$efficacy_prior_id, "\n")
  cat("Utility type:", task$utility_type, "\n")
  cat("Lambda T:", task$lambda_T, "\n")
  cat("Enrollment scheme:", task$enrollment_scheme, "\n")
  cat("Arrival rate:", task$arrival_rate, "\n")
  cat("IPDE design:", task$ipde_design, "\n")
  cat("Flexible IPDE:", task$flexible_ipde, "\n")
  cat("Folder:", folderpath, "\n")
  cat("Result file (raw preferred):", example_file, "\n")
  cat("Found", length(files.use), "files; missing", length(missing.jobs), "jobs.\n")
  cat("====================================\n")
  flush.console()

  missing.log[[miss.idx]] <- data.frame(
    task,
    Folder = folderpath,
    Example_file = example_file,
    n_found = length(files.use),
    n_expected = length(jobs_expected),
    n_missing = length(missing.jobs),
    ntrial_expected = ntrial_expected,
    ntrial_from_files = 0L,
    missing_IDX = paste(missing.jobs, collapse = ","),
    stringsAsFactors = FALSE
  )
  miss.idx <- miss.idx + 1L

  if (length(files.use) == 0L) {
    warning("No files found for task ", task$task_id, ": ", example_file)
    next
  }

  t.read <- Sys.time()
  results <- lapply(files.use, readRDS)
  ntrial_from_files <- sum(vapply(results, function(x) as.integer(x$ntrial), integer(1)))
  one <- summarize_task(results, task_id = task$task_id)

  missing.log[[miss.idx - 1L]]$ntrial_from_files <- ntrial_from_files
  cat(
    "read/summarize time:",
    round(difftime(Sys.time(), t.read, units = "secs"), 2),
    "sec\n"
  )

  if (one$ntrial != ntrial_expected) {
    cat(
      "Warning: expected", ntrial_expected,
      "trials but found", one$ntrial, "\n"
    )
  }

  key <- make_phase12_group_key(task)
  all.dose.summary[[key]] <- one$dose_summary
  all.table.summary[[key]] <- one$table_summary
}

if (length(all.dose.summary) == 0L) {
  stop("No files were found. Check results_root, filenames, and working directory.")
}

round_numeric_df <- function(dat, digits = 4L) {
  is_num <- vapply(dat, is.numeric, logical(1))
  dat[is_num] <- lapply(dat[is_num], round, digits = digits)
  dat
}

dose_summary <- do.call(rbind, all.dose.summary)
table_summary <- do.call(rbind, all.table.summary)
missing_summary <- do.call(rbind, missing.log)
job_tag <- if (length(jobs_expected) > 1L &&
               identical(jobs_expected, seq.int(min(jobs_expected), max(jobs_expected)))) {
  paste0(sprintf("%04d", min(jobs_expected)), "_to_", sprintf("%04d", max(jobs_expected)))
} else if (length(jobs_expected) == 1L) {
  sprintf("%04d", jobs_expected)
} else {
  paste0("selected_", length(jobs_expected), "_jobs")
}
tag <- paste0("AIDE_phase_I_II_IDX_", job_tag)

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
