## ============================================================
## run_oc_AIDE_phase_I_II_TITE.R
## Cluster-style JAGS TITE runner for the Phase I/II AIDE design.
##
## This deliberately follows the layout of run_oc_AIDE_phase_I_II.R:
## all editable settings and the complete task grid live in this file.
##
## Usage: Rscript phase_I_II/run_oc_AIDE_phase_I_II_TITE.R [workers] [ntrial] [IDX]
## ============================================================

args <- commandArgs(trailingOnly = TRUE)
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

fmt_short <- function(x, digits = 4L) {
  out <- format(round(as.numeric(x), digits), scientific = FALSE, trim = TRUE)
  out <- sub("(\\.[0-9]*?)0+$", "\\1", out)
  out <- sub("\\.$", "", out)
  gsub("\\.", "p", gsub("-", "m", out))
}

parse_positive_integer <- function(x, name) {
  value <- suppressWarnings(as.integer(x))
  if (length(value) != 1L || is.na(value) || value < 1L)
    stop(name, " must be a positive integer.")
  value
}

get_job_index <- function(arguments) {
  if (length(arguments) >= 3L && !is.na(suppressWarnings(as.integer(arguments[3L]))))
    return(parse_positive_integer(arguments[3L], "IDX"))
  for (name in c("JOB_I", "LSB_JOBINDEX")) {
    value <- Sys.getenv(name, unset = "")
    if (nzchar(value)) return(parse_positive_integer(value, name))
  }
  1L
}

cross_join <- function(x, y) merge(x, y, by = NULL, sort = FALSE)

read_phase12_truth <- function(file) {
  if (!file.exists(file)) stop("Scenario summary does not exist: ", file)
  truth <- utils::read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  tox_cols <- paste0("Tox_Dose", 1:5)
  eff_cols <- paste0("Eff_Dose", 1:5)
  required <- c("Scenario", "Source_Scenario", "Scenario_Group", "Attempt",
                "True_MTD_Level", "Target_Toxicity", tox_cols, eff_cols)
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
  score <- aide_compute_utility(
    e_true, p_true,
    list(type = utility_type, lambda_T = lambda_T, scores = utility_scores)
  )
  candidates <- seq_len(as.integer(true_mtd))
  candidates[which.max(score[candidates])]
}

make_phase12_tite_config_tag <- function(task) {
  paste0(
    "P12TITE-a", if (task$allocation == "two_stage") "2s" else "1s",
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
    "-ea", fmt_short(task$efficacy_ipde_alpha)
  )
}

## ============================================================
## User settings
## ============================================================

scenario_file <- "Set_5dose_adaptive_r_37_true_MTD_OBD_summary_lambda1.csv"
scenario_set_name <- tools::file_path_sans_ext(basename(scenario_file))
ntrial.total <- 1L
seed_base <- 1L
scenario_id_list <- 1:38

## Stage II is not a separate allocation mode: the two-stage design simply
## transitions to the same one-stage rule used by the one-stage design.
two_stage_sizes <- data.frame(
  allocation = "two_stage", Nmax = 30L, N_s1 = 6L, N_s2 = 30L,
  stringsAsFactors = FALSE
)
one_stage_sizes <- data.frame(
  allocation = "one_stage", Nmax = 30L, N_s1 = 30L, N_s2 = 30L,
  stringsAsFactors = FALSE
)
design_size_grid <- rbind(two_stage_sizes, one_stage_sizes)

model_grid <- data.frame(
  model_id = "multicycle_additive",
  carryover_model = "multicycle_additive",
  crm_r_model = "multicycle_additive",
  efficacy_model = "multicycle_additive",
  stringsAsFactors = FALSE
)
crm_prior_grid <- data.frame(
  crm_a_r = .15, crm_b_r = .85, stringsAsFactors = FALSE
)
efficacy_prior_grid <- data.frame(
  efficacy_a = .5, efficacy_b = .5, carry_a = .15, carry_b = .85,
  stringsAsFactors = FALSE
)
utility_grid <- rbind(
  data.frame(utility_type = 2L, lambda_T = .3),
  data.frame(utility_type = 3L, lambda_T = 1)
)
arrival_grid <- data.frame(arrival_rate = 1 / 56)
alpha_grid <- data.frame(
  toxicity_ipde_alpha = c(0, .3, .6, .9),
  efficacy_ipde_alpha = c(0, .3, .6, .9)
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
apply_ipde_efficacy_rule <- FALSE
ipde_efficacy_delta <- .05

## ============================================================
## Command-line arguments
## ============================================================

workers <- if (length(args) >= 1L) parse_positive_integer(args[1L], "workers") else 1L
if (length(args) >= 2L) ntrial.total <- parse_positive_integer(args[2L], "ntrial")
job_i <- get_job_index(args)
if (workers != 1L)
  warning("workers is retained for run_oc_AIDE_phase_I_II compatibility; execution is sequential.")
aide_phase12_tite_require_jags()

truth <- read_phase12_truth(file.path(project_root, scenario_file))
if (!is.null(scenario_id_list)) {
  scenario_id_list <- as.integer(scenario_id_list)
  missing_scenarios <- setdiff(scenario_id_list, truth$Scenario)
  if (length(missing_scenarios))
    stop("scenario_id_list is not in scenario_file: ", paste(missing_scenarios, collapse = ", "))
  truth <- truth[truth$Scenario %in% scenario_id_list, , drop = FALSE]
}
setting_grid <- Reduce(cross_join, list(
  design_size_grid, model_grid, crm_prior_grid, efficacy_prior_grid,
  utility_grid, arrival_grid, alpha_grid
))
setting_grid$setting_id <- seq_len(nrow(setting_grid))

results_root <- file.path(
  project_root, paste0("oc_results_cluster_phase_I_II_TITE_", scenario_set_name)
)
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
job_seed <- seed_base + (job_i - 1L) * ntrial.total

cat("LSF job/block index:", job_i, "\n")
cat("Trials per setting in this job:", ntrial.total, "\n")
cat("TITE Phase I/II setting combinations:", nrow(setting_grid), "\n")

## ============================================================
## Build task list and run each task
## ============================================================

tasks <- list()
task_id <- 1L
for (setting_row in seq_len(nrow(setting_grid))) {
  setting <- setting_grid[setting_row, , drop = FALSE]
  for (scenario_row in seq_len(nrow(truth))) {
    p_true <- as.numeric(unlist(truth[scenario_row, paste0("Tox_Dose", 1:5), drop = FALSE], use.names = FALSE))
    e_true <- as.numeric(unlist(truth[scenario_row, paste0("Eff_Dose", 1:5), drop = FALSE], use.names = FALSE))
    tasks[[task_id]] <- list(
      task_id = task_id, setting_id = as.integer(setting$setting_id),
      job_i = job_i, ntrial = ntrial.total, seed = job_seed,
      Scenario = as.integer(truth$Scenario[scenario_row]),
      Source_Scenario = as.integer(truth$Source_Scenario[scenario_row]),
      Scenario_Group = as.character(truth$Scenario_Group[scenario_row]),
      Attempt = as.integer(truth$Attempt[scenario_row]),
      true_mtd = as.integer(truth$True_MTD_Level[scenario_row]),
      target = as.numeric(truth$Target_Toxicity[scenario_row]),
      p_true = p_true, e_true = e_true,
      allocation = as.character(setting$allocation), Nmax = as.integer(setting$Nmax),
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
      C = cohort_size, cycle_max = cycle_max, T_assess = T_assess,
      n_eval = n_eval, m_U = m_U, dlt_dist = dlt_dist, efficacy_dist = efficacy_dist,
      crm_skeleton = crm_skeleton, utility_scores = utility_scores,
      toxicity_elimination_cutoff = toxicity_elimination_cutoff,
      efficacy_threshold = efficacy_threshold, futility_cutoff = futility_cutoff,
      min_eff_n_for_futility = min_eff_n_for_futility,
      apply_ipde_toxicity_rule = apply_ipde_toxicity_rule,
      ipde_toxicity_cutoff = ipde_toxicity_cutoff,
      apply_ipde_efficacy_rule = apply_ipde_efficacy_rule,
      ipde_efficacy_delta = ipde_efficacy_delta
    )
    tasks[[task_id]]$true_obd <- calculate_true_obd(
      p_true, e_true, tasks[[task_id]]$true_mtd,
      tasks[[task_id]]$utility_type, tasks[[task_id]]$lambda_T, utility_scores
    )
    task_id <- task_id + 1L
  }
}

run_one_phase12_tite_task <- function(task) {
  tag <- make_phase12_tite_config_tag(task)
  outdir <- file.path(results_root, tag)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  outfile <- file.path(outdir, sprintf("P12TITE-SC%02d-j%04d.rds", task$Scenario, task$job_i))
  cat("Task", task$task_id, "/", length(tasks),
      "scenario", task$Scenario, "job", task$job_i, "\n")
  result <- get_oc_sim_AIDE_phase_I_II_TITE(
    p_true = task$p_true, e_true = task$e_true, ntrial = task$ntrial,
    seed = task$seed, store_raw = FALSE,
    allocation = task$allocation, Nmax = task$Nmax, N_s1 = task$N_s1,
    N_s2 = task$N_s2, C = task$C, cycle_max = task$cycle_max,
    n_eval = task$n_eval, m_U = task$m_U,
    enrollment_scheme = "continuous",
    arrival_rate = task$arrival_rate, T_assess = task$T_assess,
    dlt_dist = task$dlt_dist, efficacy_dist = task$efficacy_dist,
    target = task$target, cutoff = task$toxicity_elimination_cutoff,
    carryover_model = task$carryover_model, crm_r_model = task$crm_r_model,
    crm_skeleton = task$crm_skeleton, crm_a_r = task$crm_a_r, crm_b_r = task$crm_b_r,
    efficacy_prior = c(task$efficacy_a, task$efficacy_b),
    efficacy_carryover_prior = c(task$carry_a, task$carry_b),
    efficacy_model = task$efficacy_model,
    efficacy_threshold = task$efficacy_threshold, futility_cutoff = task$futility_cutoff,
    min_eff_n_for_futility = task$min_eff_n_for_futility,
    utility_type = task$utility_type, lambda_T = task$lambda_T,
    utility_scores = task$utility_scores,
    toxicity_ipde_alpha = task$toxicity_ipde_alpha,
    efficacy_ipde_alpha = task$efficacy_ipde_alpha,
    apply_ipde_toxicity_rule = task$apply_ipde_toxicity_rule,
    ipde_toxicity_cutoff = task$ipde_toxicity_cutoff,
    apply_ipde_efficacy_rule = task$apply_ipde_efficacy_rule,
    ipde_efficacy_delta = task$ipde_efficacy_delta
  )
  result$task <- task
  result$truth <- list(p_true = task$p_true, e_true = task$e_true,
                       target = task$target, true_mtd = task$true_mtd,
                       true_obd = task$true_obd)
  result$runner_settings <- list(method = "JAGS TITE-AIDE",
    new_patient_priority = "new_first", n_eval_rule = "stay_and_escalate_blocked",
    seed = task$seed, ntrial = task$ntrial, T_assess = task$T_assess,
    cycle_max = task$cycle_max)
  saveRDS(result, outfile)
  cat("Saved:", outfile, "\n")
  outfile
}

files <- vapply(tasks, run_one_phase12_tite_task, character(1))
cat("Completed", length(files), "JAGS TITE Phase I/II tasks.\n")
