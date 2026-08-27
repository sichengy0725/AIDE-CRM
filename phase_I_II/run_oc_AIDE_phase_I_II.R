## ============================================================
## run_oc_AIDE_phase_I_II.R
## Cluster-style runner for the efficacy-enabled AIDE Phase I/II design.
##
## This follows run_oc_AIDE.R's construction: build an explicit task list
## for one IDX, run each task sequentially, save raw RDS files in a
## setting-specific folder, and save combined files. Progress and errors are
## written only to the scheduler's standard output/error files.
##
## Usage:
##   Rscript run_oc_AIDE_phase_I_II.R [workers] [ntrial] [IDX]
## ============================================================

lib_path <- "/rsrch8/home/biostatistics/syang10/R/x86_64-pc-linux-gnu-library/4.4"
if (dir.exists(lib_path)) {
  .libPaths(c(lib_path, .libPaths()))
}
library(rjags)
library(coda)
args <- commandArgs(trailingOnly = TRUE)

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
if (length(script_file) == 1L && file.exists(script_file)) {
  phase12_dir <- dirname(normalizePath(script_file))
} else {
  phase12_dir <- normalizePath("phase_I_II")
}
setwd(dirname(phase12_dir))

## Use the maintained Phase I/II implementation at the project root.  The
## phase_I_II copy is retained for historical runs but does not include the
## arbitrary-cycle carryover backend.
source("AIDE_phase_I_II.R")

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

  lsb_name <- Sys.getenv("LSB_JOBNAME", unset = "")
  if (nzchar(lsb_name)) {
    match <- regexpr("[0-9]+$", lsb_name)
    if (match > 0L) {
      out <- suppressWarnings(as.integer(regmatches(lsb_name, match)))
      if (!is.na(out) && out >= 1L) return(out)
    }
  }

  1L
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
    "Scenario", "Source_Scenario", "Scenario_Group", "Attempt",
    "True_MTD_Level", "Target_Toxicity", tox_cols, eff_cols
  )
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0L) {
    stop("Scenario summary is missing: ", paste(missing, collapse = ", "))
  }

  truth <- raw[, c(
    "Scenario", "Source_Scenario", "Scenario_Group", "Attempt",
    "True_MTD_Level", "Target_Toxicity", tox_cols, eff_cols,
    intersect(c("OBD_Level_Utility1", "OBD_Level_Utility2", "OBD_Level_Utility3"), names(raw))
  ), drop = FALSE]
  truth$Scenario <- as.integer(truth$Scenario)
  if (anyNA(truth$Scenario) || anyDuplicated(truth$Scenario)) {
    stop("Scenario must contain unique non-missing integer values.")
  }

  probability_columns <- c(tox_cols, eff_cols, "Target_Toxicity")
  truth[probability_columns] <- lapply(truth[probability_columns], as.numeric)
  probabilities <- as.matrix(truth[probability_columns])
  if (any(!is.finite(probabilities)) || any(probabilities < 0 | probabilities > 1)) {
    stop("All extracted toxicity and efficacy truth values must lie in [0, 1].")
  }
  truth
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
    "-om", task$allocation_mode_id,
    "-N", fmt_short(task$Nmax),
    "-s1", fmt_short(task$N_s1),
    "-u", task$utility_type,
    "-l", fmt_short(task$lambda_T),
    "-en", enrollment_tag,
    "-cm", task$carryover_model,
    "-tm", task$crm_r_model,
    "-em", task$efficacy_model,
    "-rp", fmt_param(task$crm_a_r), "x", fmt_param(task$crm_b_r),
    "-ep", fmt_param(task$efficacy_a), "x", fmt_param(task$efficacy_b),
    "-cp", fmt_param(task$carry_a), "x", fmt_param(task$carry_b),
    "-ap", fmt_param(task$efficacy_additive_alpha_a), "x",
    fmt_param(task$efficacy_additive_alpha_b),
    "-eth", fmt_short(task$efficacy_threshold),
    "-fut", fmt_short(task$futility_cutoff),
    "-w", fmt_short(task$T_assess),
    "-c", fmt_short(task$C),
    "-cyc", fmt_short(task$cycle_max),
    "-rate", fmt_short(task$arrival_rate),
    "-ta", fmt_short(task$toxicity_ipde_alpha),
    "-ea", fmt_short(task$efficacy_ipde_alpha),
    ## Generation-model codes: 1 = default, 2--4 = PDF Models 1--3.
    "-gm", task$true_generation_id,
    "-ge", task$true_random_effect_eta_tag,
    "-ga", task$true_effective_dose_alpha_tag,
    "-tg", as.integer(isTRUE(task$apply_ipde_toxicity_rule)),
    "-tc", fmt_short(task$ipde_toxicity_cutoff)
  )
}

make_phase12_folder <- function(task) make_phase12_config_tag(task)

make_phase12_group_key <- function(task) {
  paste0("P12-SC", task$Scenario, "_", make_phase12_config_tag(task))
}

combine_phase12_results <- function(files) {
  results <- lapply(files, readRDS)
  if (length(results) == 1L) return(results[[1L]])

  out <- results[[1L]]
  out$ntrial <- sum(vapply(results, function(x) as.integer(x$ntrial), integer(1)))

  vector_fields <- c(
    "MTD_by_trial", "OBD_by_trial", "early_stop_by_trial",
    "n_administrations_by_trial", "n_unique_patients_by_trial", "duration_by_trial"
  )
  matrix_fields <- c(
    "n_by_dose_by_trial", "unique_n_by_dose_by_trial", "nipde_by_dose_by_trial",
    "crm_pj_by_trial", "efficacy_pj_by_trial", "utility_by_trial", "r_hat_by_trial"
  )

  for (field in vector_fields) {
    if (all(vapply(results, function(x) !is.null(x[[field]]), logical(1)))) {
      out[[field]] <- unlist(lapply(results, `[[`, field), use.names = FALSE)
    }
  }
  for (field in matrix_fields) {
    if (all(vapply(results, function(x) !is.null(x[[field]]), logical(1)))) {
      out[[field]] <- do.call(rbind, lapply(results, `[[`, field))
    }
  }

  ndose <- length(out$p_true)
  mtd <- out$MTD_by_trial
  obd <- out$OBD_by_trial
  out$MTD_selection_percent <- 100 * tabulate(
    mtd[!is.na(mtd) & mtd >= 1L & mtd <= ndose], nbins = ndose
  ) / out$ntrial
  out$OBD_selection_percent <- 100 * tabulate(
    obd[!is.na(obd) & obd >= 1L & obd <= ndose], nbins = ndose
  ) / out$ntrial
  out$early_stop_percent <- 100 * mean(out$early_stop_by_trial)
  out$mean_administrations <- mean(out$n_administrations_by_trial)
  out$mean_administrations_by_dose <- colMeans(out$n_by_dose_by_trial)
  out$mean_unique_patients <- mean(out$n_unique_patients_by_trial)
  out$mean_unique_patients_by_dose <- colMeans(out$unique_n_by_dose_by_trial)
  out$mean_ipde_doses_by_dose <- colMeans(out$nipde_by_dose_by_trial)
  out$mean_duration <- mean(out$duration_by_trial, na.rm = TRUE)
  out$combined_files <- files
  out
}

## ============================================================
## User settings
## ============================================================

scenario_file <- "Set_5dose_adaptive_r_37_true_MTD_OBD_summary_lambda1.csv"
scenario_set_name <- tools::file_path_sans_ext(basename(scenario_file))

## Match run_oc_AIDE.R: one trial per array IDX by default.
ntrial.total <- 1L
seed_base <- 1L

## Use all scenario IDs in the truth file by default.
scenario_id_list <- 1:38

## Stage II is not a separate allocation mode: after the cumulative Stage I
## threshold and a stay decision, it uses the same one-stage rule to Nmax.
two_stage_sizes <- data.frame(
  allocation = "two_stage",
  Nmax = c(30L),
  N_s1 = c(6L),
  N_s2 = c(30L),
  stringsAsFactors = FALSE
)
one_stage_sizes <- data.frame(
  allocation = "one_stage",
  Nmax = c(30L),
  N_s1 = c(30L),
  N_s2 = c(30L),
  stringsAsFactors = FALSE
)
design_size_grid <- rbind(two_stage_sizes, one_stage_sizes)

## One-stage allocation options: option 1 is the historical upward
## no-skipping rule; option 2 moves one non-futile level toward the OBD.
allocation_mode_grid <- data.frame(
  allocation_mode_id = c(1L, 2L),
  allocation_mode = c("upward_no_skipping", "one_level_toward_obd"),
  stringsAsFactors = FALSE
)

## The arbitrary-cycle shared additive model propagates each patient's
## uncapped endpoint state across all administrations.
model_grid <- data.frame(
  model_id = c("multicycle_additive"),
  carryover_model = c("multicycle_additive"),
  crm_r_model = c("multicycle_additive"),
  efficacy_model = c(
    "multicycle_additive"
  ),
  stringsAsFactors = FALSE
)
## Beta(a, b) prior for random CRM r or the additive toxicity alpha,
## depending on crm_r_model.
crm_prior_grid <- data.frame(
  crm_prior_id = "r_beta_0p15_0p85",
  crm_a_r = 0.15,
  crm_b_r = 0.85,
  stringsAsFactors = FALSE
)
## The carryover prior applies to the discount efficacy models; the
## additive-alpha prior applies to the additive efficacy models.
efficacy_prior_grid <- data.frame(
  efficacy_prior_id = "regular_beta_0p5_0p5_carry_beta_1_1",
  efficacy_a = c(0.5),
  efficacy_b = c(0.5),
  carry_a = 0.15,
  carry_b = 0.85,
  efficacy_additive_alpha_a = c(0.15),
  efficacy_additive_alpha_b = c(0.85),
  stringsAsFactors = FALSE
)
enrollment_schemes <- data.frame(
  enrollment_scheme = c("continuous"),
  stringsAsFactors = FALSE
)
lambda_T_grid <- 0.3
utility_grid <- rbind(
  # data.frame(utility_type = 1L, lambda_T = 1, stringsAsFactors = FALSE),
  data.frame(utility_type = 2L, lambda_T = lambda_T_grid, stringsAsFactors = FALSE),
  data.frame(utility_type = 3L, lambda_T = 1, stringsAsFactors = FALSE)
)
arrival_grid <- data.frame(arrival_rate = c(1 / 56), stringsAsFactors = FALSE)
alpha_grid <- data.frame(
  toxicity_ipde_alpha = c(0, 0.3, 0.6, 0.9),
  efficacy_ipde_alpha = c(0, 0.3, 0.6, 0.9),
  stringsAsFactors = FALSE
)

## Alternative-truth parameters. Model 1 uses its fixed source-dose alpha
## vector; Model 2 uses editable eta; Model 3 alpha is user-controlled.
true_dose_specific_alpha <- c(0.2, 0.3, 0.4, 0.5, 0.6)
true_random_effect_eta <- 1
true_effective_dose_values <- c(15, 20, 30, 35, 45)
true_effective_dose_alpha_grid <- c(0.3, 0.6, 0.9)

## The legacy generator is evaluated over its endpoint-alpha grid. Each
## alternative DGM has only its relevant settings, avoiding redundant runs.
true_generation_grid <- rbind(
  data.frame(
    true_generation_id = 1L,
    true_generation = "legacy",
    alpha_grid,
    true_random_effect_eta = true_random_effect_eta,
    true_random_effect_eta_tag = "na",
    true_effective_dose_alpha = 0,
    true_effective_dose_alpha_tag = "na",
    stringsAsFactors = FALSE
  ),
  data.frame(
    true_generation_id = 2L,
    true_generation = "dose_specific_geometric",
    toxicity_ipde_alpha = 0,
    efficacy_ipde_alpha = 0,
    true_random_effect_eta = true_random_effect_eta,
    true_random_effect_eta_tag = "na",
    true_effective_dose_alpha = 0,
    true_effective_dose_alpha_tag = "na",
    stringsAsFactors = FALSE
  ),
  data.frame(
    true_generation_id = 3L,
    true_generation = "shared_patient_logistic",
    toxicity_ipde_alpha = 0,
    efficacy_ipde_alpha = 0,
    true_random_effect_eta = true_random_effect_eta,
    true_random_effect_eta_tag = fmt_short(true_random_effect_eta),
    true_effective_dose_alpha = 0,
    true_effective_dose_alpha_tag = "na",
    stringsAsFactors = FALSE
  ),
  data.frame(
    true_generation_id = 4L,
    true_generation = "effective_dose_geometric",
    toxicity_ipde_alpha = 0,
    efficacy_ipde_alpha = 0,
    true_random_effect_eta = true_random_effect_eta,
    true_random_effect_eta_tag = "na",
    true_effective_dose_alpha = true_effective_dose_alpha_grid,
    true_effective_dose_alpha_tag = fmt_short(true_effective_dose_alpha_grid),
    stringsAsFactors = FALSE
  )
)

utility_scores <- c(u00 = 0, u01 = 40, u10 = 60, u11 = 100)
cohort_size <- 3L
cycle_max <- 2L
T_assess <- 28
crm_skeleton <- c(0.15, 0.20, 0.30, 0.35, 0.45)
crm_alpha_sd <- sqrt(2)
efficacy_threshold <- 0.20
futility_cutoff <- 0.85
min_eff_n_for_futility <- 0L
apply_ipde_toxicity_rule <- TRUE
ipde_toxicity_cutoff <- 0.95

## ============================================================
## Command-line arguments, matching run_oc_AIDE.R
## ============================================================

requested_workers <- if (length(args) >= 1L) {
  parse_positive_integer(args[1L], "workers")
} else {
  1L
}
if (length(args) >= 2L) ntrial.total <- parse_positive_integer(args[2L], "ntrial")
job_i <- get_job_index(args)
if (requested_workers != 1L) {
  warning("workers is retained for run_oc_AIDE compatibility; this runner is sequential.")
}
job_seed <- seed_base + (job_i - 1L) * ntrial.total
if (job_seed > .Machine$integer.max) {
  stop("The IDX and ntrial combination produces a seed outside R's integer range.")
}

if (!requireNamespace("rjags", quietly = TRUE)) {
  stop("This runner requires rjags and a working JAGS installation.")
}
if (!requireNamespace("coda", quietly = TRUE)) {
  stop("This runner requires the coda package.")
}

truth <- read_phase12_truth(scenario_file)
scenario_id_list <- as.integer(scenario_id_list)
missing_scenarios <- setdiff(scenario_id_list, truth$Scenario)
if (length(missing_scenarios) > 0L) {
  stop("scenario_id_list is not in scenario_file: ", paste(missing_scenarios, collapse = ", "))
}
truth <- truth[truth$Scenario %in% scenario_id_list, , drop = FALSE]

setting_grid <- Reduce(
  cross_join,
  list(
    design_size_grid,
    allocation_mode_grid,
    model_grid,
    crm_prior_grid,
    efficacy_prior_grid,
    enrollment_schemes,
    utility_grid,
    arrival_grid,
    true_generation_grid
  )
)
setting_grid$setting_id <- seq_len(nrow(setting_grid))

workdir <- getwd()
outdir_root <- paste0("oc_results_cluster_phase_I_II_", scenario_set_name)
dir.create(outdir_root, recursive = TRUE, showWarnings = FALSE)

cat("LSF job/block index:", job_i, "\n")
cat("Working directory:", workdir, "\n")
cat("Execution mode: sequential\n")
cat("Requested workers argument:", requested_workers, "\n")
cat("Trials per setting in this job:", ntrial.total, "\n")
cat("Scenario set:", scenario_set_name, "\n")
cat("Scenario IDs:", paste(scenario_id_list, collapse = ", "), "\n")
cat("Phase I/II setting combinations:", nrow(setting_grid), "\n")
cat("First seed for this job:", job_seed, "\n")

## ============================================================
## Worker function
## ============================================================

run_one_phase12_task <- function(task) {
  setwd(task$workdir)
  source("AIDE_phase_I_II.R")
  library(rjags)
  library(coda)

  foldername <- make_phase12_folder(task)
  outdir <- file.path(task$outdir_root, foldername)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  scenario_name <- paste0("P12-SC", task$Scenario)
  method_tag <- paste0("phase12_", task$model_id)
  filename <- paste0(
    scenario_name,
    "-CRM-", method_tag,
    "-j", task$job_i,
    "-b", task$block_id,
    "-s", task$seed.block,
    "-n", task$ntrial.block,
    ".rds"
  )
  outfile <- file.path(outdir, filename)
  cat("Task", task$task_id, "/", task$setting_id,
      "scenario", task$Scenario, "job", task$job_i, "\n")

  result <- get_oc_sim_AIDE_phase_I_II(
      p_true = task$p_true,
      e_true = task$e_true,
      ntrial = task$ntrial.block,
      seed = task$seed.block,
      store_raw = FALSE,

      allocation = task$allocation,
      allocation_mode = task$allocation_mode,
      Nmax = task$Nmax,
      N_s1 = task$N_s1,
      N_s2 = task$N_s2,
      N_pat = task$Nmax,
      C = task$C,
      cycle_max = task$cycle_max,

      model = "CRM",
      target = task$target,
      carryover_model = task$carryover_model,
      crm_r_model = task$crm_r_model,
      crm_skeleton = task$crm_skeleton,
      crm_alpha_sd = task$crm_alpha_sd,
      crm_a_r = task$crm_a_r,
      crm_b_r = task$crm_b_r,

      efficacy_prior = c(task$efficacy_a, task$efficacy_b),
      efficacy_carryover_prior = c(task$carry_a, task$carry_b),
      efficacy_model = task$efficacy_model,
      efficacy_additive_alpha_prior = c(
        task$efficacy_additive_alpha_a, task$efficacy_additive_alpha_b
      ),
      utility_type = task$utility_type,
      lambda_T = task$lambda_T,
      utility_scores = task$utility_scores,

      enrollment_scheme = task$enrollment_scheme,
      arrival_rate = task$arrival_rate,
      T_assess = task$T_assess,

      toxicity_ipde_alpha = task$toxicity_ipde_alpha,
      efficacy_ipde_alpha = task$efficacy_ipde_alpha,
      true_generation = task$true_generation,
      true_dose_specific_alpha = task$true_dose_specific_alpha,
      true_random_effect_eta = task$true_random_effect_eta,
      true_effective_dose_values = task$true_effective_dose_values,
      true_effective_dose_alpha = task$true_effective_dose_alpha,

      efficacy_threshold = task$efficacy_threshold,
      futility_cutoff = task$futility_cutoff,
      min_eff_n_for_futility = task$min_eff_n_for_futility,
      apply_ipde_toxicity_rule = task$apply_ipde_toxicity_rule,
      ipde_toxicity_cutoff = task$ipde_toxicity_cutoff
  )

  result$task <- task
  result$truth <- list(
    p_true = task$p_true,
    e_true = task$e_true,
    target = task$target,
    true_mtd = task$true_mtd,
    true_obd = task$true_obd,
    utility_type = task$utility_type
  )
  result$runner_settings <- list(
    model_id = task$model_id,
    carryover_model = task$carryover_model,
    crm_r_model = task$crm_r_model,
    efficacy_model = task$efficacy_model,
    allocation_mode_id = task$allocation_mode_id,
    allocation_mode = task$allocation_mode,
    true_generation_id = task$true_generation_id,
    true_generation = task$true_generation,
    true_random_effect_eta_tag = task$true_random_effect_eta_tag,
    true_effective_dose_alpha_tag = task$true_effective_dose_alpha_tag,
    true_generation_parameters = list(
      dose_specific_alpha = task$true_dose_specific_alpha,
      random_effect_eta = task$true_random_effect_eta,
      effective_dose_values = task$true_effective_dose_values,
      effective_dose_alpha = task$true_effective_dose_alpha
    ),
    efficacy_additive_alpha_prior = c(
      task$efficacy_additive_alpha_a, task$efficacy_additive_alpha_b
    ),
    utility_scores = task$utility_scores,
    stage2_allocation = "one_stage",
    cycle_max = task$cycle_max,
    T_assess = task$T_assess,
    multicycle_ipde_toxicity_gate = list(
      enabled = task$apply_ipde_toxicity_rule,
      cutoff = task$ipde_toxicity_cutoff
    ),
    job_i = task$job_i,
    block_id = task$block_id,
    ntrial = task$ntrial.block,
    seed = task$seed.block
  )
  saveRDS(result, outfile)
  cat("Saved:", outfile, "\n")

  list(
    file = outfile,
    folder = outdir,
    task = task
  )
}

## ============================================================
## Build task list
## ============================================================

tasks <- list()
task_id <- 1L

for (setting_row in seq_len(nrow(setting_grid))) {
  setting <- setting_grid[setting_row, , drop = FALSE]
  for (scenario_row in seq_len(nrow(truth))) {
    p_true <- as.numeric(unlist(
      truth[scenario_row, paste0("Tox_Dose", 1:5), drop = FALSE],
      use.names = FALSE
    ))
    e_true <- as.numeric(unlist(
      truth[scenario_row, paste0("Eff_Dose", 1:5), drop = FALSE],
      use.names = FALSE
    ))
    true_mtd <- as.integer(truth$True_MTD_Level[scenario_row])
    utility_type <- as.integer(setting$utility_type)
    lambda_T <- as.numeric(setting$lambda_T)

    tasks[[task_id]] <- list(
      task_id = task_id,
      setting_id = as.integer(setting$setting_id),
      workdir = workdir,
      outdir_root = outdir_root,
      scenario_set = scenario_set_name,
      job_i = job_i,
      block_id = 1L,
      ntrial.block = ntrial.total,
      seed.block = job_seed,

      Scenario = as.integer(truth$Scenario[scenario_row]),
      Source_Scenario = as.integer(truth$Source_Scenario[scenario_row]),
      Scenario_Group = as.character(truth$Scenario_Group[scenario_row]),
      Attempt = as.integer(truth$Attempt[scenario_row]),
      true_mtd = true_mtd,
      true_obd = calculate_true_obd(
        p_true, e_true, true_mtd, utility_type, lambda_T, utility_scores
      ),
      target = as.numeric(truth$Target_Toxicity[scenario_row]),
      p_true = p_true,
      e_true = e_true,

      allocation = as.character(setting$allocation),
      allocation_mode_id = as.integer(setting$allocation_mode_id),
      allocation_mode = as.character(setting$allocation_mode),
      Nmax = as.integer(setting$Nmax),
      N_s1 = as.integer(setting$N_s1),
      N_s2 = as.integer(setting$N_s2),
      model_id = as.character(setting$model_id),
      carryover_model = as.character(setting$carryover_model),
      crm_r_model = as.character(setting$crm_r_model),
      crm_prior_id = as.character(setting$crm_prior_id),
      crm_a_r = as.numeric(setting$crm_a_r),
      crm_b_r = as.numeric(setting$crm_b_r),
      efficacy_prior_id = as.character(setting$efficacy_prior_id),
      efficacy_a = as.numeric(setting$efficacy_a),
      efficacy_b = as.numeric(setting$efficacy_b),
      carry_a = as.numeric(setting$carry_a),
      carry_b = as.numeric(setting$carry_b),
      efficacy_model = as.character(setting$efficacy_model),
      efficacy_additive_alpha_a = as.numeric(setting$efficacy_additive_alpha_a),
      efficacy_additive_alpha_b = as.numeric(setting$efficacy_additive_alpha_b),
      utility_type = utility_type,
      lambda_T = lambda_T,
      enrollment_scheme = as.character(setting$enrollment_scheme),
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

      C = cohort_size,
      cycle_max = cycle_max,
      T_assess = T_assess,
      crm_skeleton = crm_skeleton,
      crm_alpha_sd = crm_alpha_sd,
      utility_scores = utility_scores,
      efficacy_threshold = efficacy_threshold,
      futility_cutoff = futility_cutoff,
      min_eff_n_for_futility = min_eff_n_for_futility,
      apply_ipde_toxicity_rule = apply_ipde_toxicity_rule,
      ipde_toxicity_cutoff = ipde_toxicity_cutoff
    )
    task_id <- task_id + 1L
  }
}

cat("Number of tasks:", length(tasks), "\n")

## ============================================================
## Run jobs sequentially
## ============================================================

task_results <- lapply(tasks, run_one_phase12_task)
cat("\nAll sequential Phase I/II tasks completed.\n")

## ============================================================
## Combine chunk results by setting within this LSF job
## ============================================================

group_key <- vapply(task_results, function(x) make_phase12_group_key(x$task), character(1))
groups <- split(task_results, group_key)

for (group_name in names(groups)) {
  files <- vapply(groups[[group_name]], function(x) x$file, character(1))
  folder <- groups[[group_name]][[1L]]$folder
  combined <- combine_phase12_results(files)
  combined_file <- file.path(folder, paste0(group_name, "-job-", job_i, "-combined.rds"))
  saveRDS(combined, combined_file)
  cat("Combined:", combined_file, "\n")
  cat("  ntrial.total =", combined$ntrial, "\n")
}

cat("\nAll AIDE Phase I/II cluster runs completed.\n")
