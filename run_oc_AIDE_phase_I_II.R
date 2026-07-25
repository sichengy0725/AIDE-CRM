## ============================================================
## run_oc_AIDE_phase_I_II.R
##
## Array-friendly operating-characteristic runner for the proposed
## efficacy-enabled AIDE Phase I/II design.  Truth is extracted directly
## from Set_5dose_adaptive_r_37_true_MTD_OBD_summary_lambda1.csv.
##
## Usage (matching run_oc_AIDE.R):
##   Rscript run_oc_AIDE_phase_I_II.R [workers] [ntrial] [IDX]
##
## Examples:
##   Rscript run_oc_AIDE_phase_I_II.R 1 100 1
##   Rscript run_oc_AIDE_phase_I_II.R 1 100 2
##
## Each IDX runs every scenario/setting task with a non-overlapping seed
## block. Thus IDX = 1,...,10 with ntrial = 100 yields 1,000 trials per task.
## ============================================================

args <- commandArgs(trailingOnly = TRUE)

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
if (length(script_file) == 1L && file.exists(script_file)) {
  setwd(dirname(normalizePath(script_file)))
}

source("AIDE_phase_I_II.R")

parse_positive_integer <- function(x, name) {
  value <- suppressWarnings(as.integer(x))
  if (length(value) != 1L || is.na(value) || value < 1L) {
    stop(name, " must be a positive integer.")
  }
  value
}

get_job_index <- function(args) {
  if (length(args) >= 3L) {
    out <- suppressWarnings(as.integer(args[3L]))
    if (!is.na(out) && out >= 1L) return(out)
  }

  for (name in c("JOB_I", "LSB_JOBINDEX")) {
    out <- suppressWarnings(as.integer(Sys.getenv(name, unset = "")))
    if (!is.na(out) && out >= 1L) return(out)
  }

  1L
}

fmt <- function(x, digits = 3L) {
  out <- format(round(as.numeric(x), digits), trim = TRUE, scientific = FALSE)
  out <- sub("(\\.[0-9]*?)0+$", "\\1", out)
  out <- sub("\\.$", "", out)
  gsub("\\.", "p", out)
}

cross_join <- function(x, y) {
  merge(x, y, by = NULL, sort = FALSE)
}

calculate_true_obd <- function(p_true, e_true, true_mtd, utility_type,
                               lambda_T, utility_scores) {
  if (length(p_true) != length(e_true) || true_mtd < 1L ||
      true_mtd > length(p_true)) {
    return(NA_integer_)
  }
  candidates <- seq_len(as.integer(true_mtd))
  utility <- aide_phase12_compute_utility(
    efficacy = e_true,
    toxicity = p_true,
    utility_type = utility_type,
    lambda_T = lambda_T,
    utility_scores = utility_scores
  )
  candidates[which.max(utility[candidates])]
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

  truth <- raw[, unique(c(
    "Scenario", "Source_Scenario", "Scenario_Group", "Attempt",
    "True_MTD_Level", "Target_Toxicity", tox_cols, eff_cols,
    "OBD_Level_Utility1", "OBD_Level_Utility2", "OBD_Level_Utility3"
  ))[unique(c(
    "Scenario", "Source_Scenario", "Scenario_Group", "Attempt",
    "True_MTD_Level", "Target_Toxicity", tox_cols, eff_cols,
    "OBD_Level_Utility1", "OBD_Level_Utility2", "OBD_Level_Utility3"
  )) %in% names(raw)], drop = FALSE]

  truth$Scenario <- as.integer(truth$Scenario)
  if (anyNA(truth$Scenario) || anyDuplicated(truth$Scenario)) {
    stop("Scenario must contain unique non-missing integer values.")
  }
  probability_columns <- c(tox_cols, eff_cols, "Target_Toxicity")
  truth[probability_columns] <- lapply(truth[probability_columns], as.numeric)
  if (any(!is.finite(as.matrix(truth[probability_columns]))) ||
      any(as.matrix(truth[probability_columns]) < 0 | as.matrix(truth[probability_columns]) > 1)) {
    stop("All extracted toxicity and efficacy truth values must lie in [0, 1].")
  }
  truth
}

## ------------------------------------------------------------
## User settings: edit these compact grids for a smaller/larger run.
## ------------------------------------------------------------

scenario_file <- "Set_5dose_adaptive_r_37_true_MTD_OBD_summary_lambda1.csv"
output_dir <- "oc_results_AIDE_phase_I_II"
requested_workers <- if (length(args) >= 1L) {
  parse_positive_integer(args[1L], "workers")
} else {
  1L
}
ntrial <- if (length(args) >= 2L) parse_positive_integer(args[2L], "ntrial") else 100L
job_i <- get_job_index(args)
seed_base <- 1L
job_seed <- seed_base + (job_i - 1L) * ntrial
if (job_seed > .Machine$integer.max) {
  stop("The IDX and ntrial combination produces a seed outside R's integer range.")
}

## Two-stage uses paired (Nmax, N_s1) settings.  One-stage uses only Nmax;
## N_s1 is set to Nmax because it is unused by that allocation rule.
two_stage_sizes <- data.frame(
  allocation = "two_stage",
  Nmax = c(30L, 45L, 60L),
  N_s1 = c(15L, 24L, 30L),
  stringsAsFactors = FALSE
)
one_stage_sizes <- data.frame(
  allocation = "one_stage",
  Nmax = c(30L, 45L, 60L),
  N_s1 = c(30L, 45L, 60L),
  stringsAsFactors = FALSE
)
design_size_grid <- rbind(two_stage_sizes, one_stage_sizes)

## Priors for the random-carryover CRM.
crm_prior_grid <- data.frame(
  crm_prior_id = "r_beta_0p15_0p85",
  crm_a_r = 0.15,
  crm_b_r = 0.85,
  stringsAsFactors = FALSE
)

## Independent Beta priors for regular efficacy and dose-specific efficacy
## carryover.  Each row supplies Beta(a, b) for every dose.
efficacy_prior_grid <- data.frame(
  efficacy_prior_id = "regular_beta_0p15_0p85_carry_beta_0p15_0p85",
  efficacy_a = 0.15,
  efficacy_b = 0.85,
  carry_a = 0.15,
  carry_b = 0.85,
  stringsAsFactors = FALSE
)

enrollment_schemes <- data.frame(
  enrollment_scheme = c("continuous", "ipde_first"),
  stringsAsFactors = FALSE
)
## lambda_T is the toxicity penalty in Utility 2 only. Utility 1 and Utility
## 3 do not depend on lambda_T, so they are each run once.
lambda_T_grid <- 0.3
utility_grid <- rbind(
  data.frame(utility_type = 1L, lambda_T = 1, stringsAsFactors = FALSE),
  data.frame(utility_type = 2L, lambda_T = lambda_T_grid, stringsAsFactors = FALSE),
  data.frame(utility_type = 3L, lambda_T = 1, stringsAsFactors = FALSE)
)
## One new-patient arrival per 56 days, expressed as a per-day rate.
arrival_grid <- data.frame(arrival_rate = 1 / 56, stringsAsFactors = FALSE)
ipde_grid <- data.frame(
  ipde_design = c(1L, 2L),
  flexible_ipde = c(TRUE, FALSE),
  stringsAsFactors = FALSE
)
alpha_grid <- data.frame(
  toxicity_ipde_alpha = c(0, 0.3, 0.6, 0.9),
  efficacy_ipde_alpha = c(0, 0.3, 0.6, 0.9),
  stringsAsFactors = FALSE
)

## Utility 3 scores reproduce the values used in the supplied lambda-1
## summary: u00 = 40, u01 = 100, u10 = 0, u11 = 60.
utility_scores <- c(u00 = 40, u01 = 100, u10 = 0, u11 = 60)

## Fixed design settings kept outside the grid.
cohort_size <- 3L
cycle_max <- 2L
T_assess <- 28
m_U <- 6L
crm_skeleton <- c(0.15, 0.20, 0.30, 0.35, 0.45)
crm_alpha_sd <- sqrt(2)
efficacy_threshold <- 0.20
futility_cutoff <- 0.95
min_eff_n_for_futility <- 0L

## No optional random-CRM IPDE safety or efficacy gates are applied here.
## The Phase I/II simulator's intrinsic toxicity and futility monitoring is
## left intact, because it is part of the design being evaluated.
apply_random_crm_recycle_toxicity_rule <- FALSE
apply_random_crm_recycle_efficacy_rule <- FALSE

## ------------------------------------------------------------
## Extract truth and build one task per scenario/setting combination.
## ------------------------------------------------------------

truth <- read_phase12_truth(scenario_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

truth_export <- truth
names(truth_export)[match(paste0("Tox_Dose", 1:5), names(truth_export))] <-
  paste0("true_toxicity_dose", 1:5)
names(truth_export)[match(paste0("Eff_Dose", 1:5), names(truth_export))] <-
  paste0("true_efficacy_dose", 1:5)
truth_export_file <- file.path(
  output_dir, "Set_5dose_adaptive_r_37_phase12_truth.csv"
)
if (job_i == 1L || !file.exists(truth_export_file)) {
  utils::write.csv(truth_export, truth_export_file, row.names = FALSE)
}

setting_grid <- Reduce(
  cross_join,
  list(
    design_size_grid,
    crm_prior_grid,
    efficacy_prior_grid,
    enrollment_schemes,
    utility_grid,
    arrival_grid,
    ipde_grid,
    alpha_grid
  )
)
tasks <- cross_join(
  truth[, c("Scenario", "Source_Scenario", "Scenario_Group", "Attempt"), drop = FALSE],
  setting_grid
)
tasks$task_id <- seq_len(nrow(tasks))
tasks <- tasks[, c("task_id", setdiff(names(tasks), "task_id")), drop = FALSE]

manifest_file <- file.path(output_dir, "AIDE_phase_I_II_task_manifest.csv")
if (job_i == 1L || !file.exists(manifest_file)) {
  utils::write.csv(tasks, manifest_file, row.names = FALSE)
}

cat("Extracted scenarios:", nrow(truth), "\n")
cat("Phase I/II tasks:", nrow(tasks), "\n")
cat("Trials per task in this job:", ntrial, "\n")
cat("Job IDX:", job_i, "\n")
cat("First seed for this job:", job_seed, "\n")
if (requested_workers != 1L) {
  warning("workers is retained for run_oc_AIDE compatibility; this runner is sequential.")
}

if (!requireNamespace("rjags", quietly = TRUE)) {
  stop("This runner requires rjags and a working JAGS installation.")
}
if (!requireNamespace("coda", quietly = TRUE)) {
  stop("This runner requires the coda package.")
}

for (task_row in seq_len(nrow(tasks))) {
  task <- tasks[task_row, , drop = FALSE]
  scenario_row <- match(task$Scenario, truth$Scenario)
  p_true <- as.numeric(unlist(
    truth[scenario_row, paste0("Tox_Dose", 1:5), drop = FALSE],
    use.names = FALSE
  ))
  e_true <- as.numeric(unlist(
    truth[scenario_row, paste0("Eff_Dose", 1:5), drop = FALSE],
    use.names = FALSE
  ))
  true_mtd <- as.integer(truth$True_MTD_Level[scenario_row])
  true_obd <- calculate_true_obd(
    p_true = p_true,
    e_true = e_true,
    true_mtd = true_mtd,
    utility_type = as.integer(task$utility_type),
    lambda_T = task$lambda_T,
    utility_scores = utility_scores
  )

  cat(
    "Running task", task$task_id, "of", nrow(tasks),
    "for scenario", task$Scenario, "\n"
  )
  result <- get_oc_sim_AIDE_phase_I_II(
    p_true = p_true,
    e_true = e_true,
    ntrial = ntrial,
    seed = job_seed,
    store_raw = FALSE,

    allocation = task$allocation,
    Nmax = as.integer(task$Nmax),
    N_s1 = as.integer(task$N_s1),
    N_pat = as.integer(task$Nmax),
    C = cohort_size,
    cycle_max = cycle_max,
    m_U = m_U,

    model = "CRM",
    target = as.numeric(truth$Target_Toxicity[scenario_row]),
    crm_r_model = "random",
    crm_skeleton = crm_skeleton,
    crm_alpha_sd = crm_alpha_sd,
    crm_a_r = task$crm_a_r,
    crm_b_r = task$crm_b_r,

    efficacy_prior = c(task$efficacy_a, task$efficacy_b),
    efficacy_carryover_prior = c(task$carry_a, task$carry_b),
    utility_type = as.integer(task$utility_type),
    lambda_T = task$lambda_T,
    utility_scores = utility_scores,

    enrollment_scheme = task$enrollment_scheme,
    arrival_rate = task$arrival_rate,
    T_assess = T_assess,

    ipde_design = as.integer(task$ipde_design),
    flexible_ipde = as.logical(task$flexible_ipde),
    toxicity_ipde_alpha = task$toxicity_ipde_alpha,
    efficacy_ipde_alpha = task$efficacy_ipde_alpha,

    efficacy_threshold = efficacy_threshold,
    futility_cutoff = futility_cutoff,
    min_eff_n_for_futility = min_eff_n_for_futility,
    apply_random_crm_recycle_toxicity_rule = apply_random_crm_recycle_toxicity_rule,
    apply_random_crm_recycle_efficacy_rule = apply_random_crm_recycle_efficacy_rule
  )

  result$task <- as.list(task)
  result$truth <- list(
    p_true = p_true,
    e_true = e_true,
    target = truth$Target_Toxicity[scenario_row],
    true_mtd = true_mtd,
    true_obd = true_obd,
    utility_type = as.integer(task$utility_type)
  )
  result$runner_settings <- list(
    utility_scores = utility_scores,
    lambda_T_grid = lambda_T_grid,
    m_U = m_U,
    optional_random_crm_ipde_gates = FALSE,
    job_i = job_i,
    ntrial = ntrial,
    seed = job_seed
  )

  outfile <- file.path(
    output_dir,
    paste0(
      "task_", sprintf("%06d", task$task_id),
      "_scenario_", task$Scenario,
      "_", task$allocation,
      "_Nmax", task$Nmax,
      "_utility", task$utility_type,
      "_", task$enrollment_scheme,
      "_job_", sprintf("%04d", job_i),
      ".rds"
    )
  )
  saveRDS(result, outfile)
  cat("Saved:", outfile, "\n")
}
