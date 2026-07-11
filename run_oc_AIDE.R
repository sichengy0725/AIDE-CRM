## ============================================================
## run_oc_AIDE.R
## Cluster runner for AIDE-BOIN / AIDE-CRM / AIDE-CFO
## Uses the prior settings from methods_prior.R
## Compatible with job_01 to job_02000.lsf
## ============================================================

rm(list = ls())

## -------------------------------
## Working directory / library path
## -------------------------------

# setwd("/rsrch8/home/biostatistics/syang10/AIDE")

lib_path <- "/rsrch8/home/biostatistics/syang10/R/x86_64-pc-linux-gnu-library/4.4"
if (dir.exists(lib_path)) {
  .libPaths(c(lib_path, .libPaths()))
}

## -------------------------------
## Source project files
## -------------------------------

## Newest pooled-r_mle BOIN code.
## Required files in this working directory:
##   AIDE_CRM_helper_final.R
##   AIDE_BOIN_helper.R
##   AIDE_modified.R
source("AIDE_BOIN_helper.R")
source("AIDE_CRM_helper_final.R")
source("AIDE_modified.R")

## If CRM or PRIDE-backed CFO is used, load rjags/coda when available.
## Loading here is safe even for BOIN if installed.
if (requireNamespace("rjags", quietly = TRUE)) library(rjags)
if (requireNamespace("coda", quietly = TRUE)) library(coda)

## ============================================================
## Helper functions
## ============================================================

fmt_param <- function(x, digits = 4) {
  x <- round(as.numeric(x), digits)
  gsub("\\.", "p", format(x, scientific = FALSE, trim = TRUE))
}

fmt_num <- function(x) {
  gsub("\\.", "p", as.character(x))
}

## Short numeric formatter for folder/file names.
## Rounds to 2 digits by default and removes unnecessary trailing zeros.
fmt_short <- function(x, digits = 2) {
  x <- round(as.numeric(x), digits)
  out <- format(x, scientific = FALSE, trim = TRUE)
  ## Remove trailing zeros only after a decimal point; preserve exact zero as "0".
  out <- sub("(\\.[0-9]*?)0+$", "\\1", out)
  out <- sub("\\.$", "", out)
  out[out == "-0"] <- "0"
  out[out == ""] <- "0"
  out <- gsub("-", "m", out)
  out <- gsub("\\.", "p", out)
  out
}

make_random_scenario_filename <- function(
    ndose,
    nscenario,
    target,
    target_diff_below,
    target_diff_above,
    out_dir = "scenario_sets"
) {
  file.path(
    out_dir,
    paste0(
      "random_scenarios",
      "_ndose", as.integer(ndose),
      "_nscenario", as.integer(nscenario),
      "_target", fmt_short(target, digits = 4),
      "_tdiffbelow", fmt_short(target_diff_below, digits = 4),
      "_tdiffabove", fmt_short(target_diff_above, digits = 4),
      ".csv"
    )
  )
}

with_sink <- function(file, expr) {
  con <- file(file, open = "wt")
  sink(con)
  sink(con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(con)
  }, add = TRUE)
  force(expr)
}

get_job_index <- function(args) {
  if (length(args) >= 3) {
    out <- suppressWarnings(as.integer(args[3]))
    if (!is.na(out) && out >= 1L) return(out)
  }

  env_job_i <- Sys.getenv("JOB_I", unset = "")
  if (nzchar(env_job_i)) {
    out <- suppressWarnings(as.integer(env_job_i))
    if (!is.na(out) && out >= 1L) return(out)
  }

  env_lsb_index <- Sys.getenv("LSB_JOBINDEX", unset = "")
  if (nzchar(env_lsb_index)) {
    out <- suppressWarnings(as.integer(env_lsb_index))
    if (!is.na(out) && out >= 1L) return(out)
  }

  env_lsb_name <- Sys.getenv("LSB_JOBNAME", unset = "")
  if (nzchar(env_lsb_name)) {
    m <- regexpr("[0-9]+$", env_lsb_name)
    if (m > 0) {
      out <- suppressWarnings(as.integer(regmatches(env_lsb_name, m)))
      if (!is.na(out) && out >= 1L) return(out)
    }
  }

  1L
}

get_scenario_batch_index <- function(args, default = 1L) {
  if (length(args) >= 4) {
    out <- suppressWarnings(as.integer(args[4]))
    if (!is.na(out) && out >= 1L) return(out)
    stop("Scenario batch index in args[4] must be a positive integer.")
  }

  env_batch_i <- Sys.getenv("SCENARIO_BATCH_I", unset = "")
  if (nzchar(env_batch_i)) {
    out <- suppressWarnings(as.integer(env_batch_i))
    if (!is.na(out) && out >= 1L) return(out)
    stop("SCENARIO_BATCH_I must be a positive integer.")
  }

  if (length(args) >= 3) {
    out <- suppressWarnings(as.integer(args[3]))
    if (!is.na(out) && out >= 1L) return(out)
  }

  default
}

make_p_ipde <- function(p_base, alpha_true) {
  if (length(alpha_true) != 1L || !is.finite(alpha_true) || alpha_true < 0) {
    stop("alpha_true must be a single finite nonnegative value.")
  }

  p_base <- as.numeric(p_base)

  ## IPDE truth: dose j has current-dose toxicity plus carryover from
  ## the immediately previous dose. Dose 1 has no previous-dose term.
  prev_dose_p <- c(0, p_base[-length(p_base)])

  pmin(p_base + alpha_true * prev_dose_p, 1)
}

load_scenario_file <- function(file) {
  if (!is.character(file) || length(file) != 1L || !nzchar(file)) {
    stop("file must be a single non-empty character string.")
  }

  if (!file.exists(file)) {
    stop("Scenario file does not exist: ", file)
  }

  file_ext <- tolower(sub("^.*\\.([^.]*)$", "\\1", basename(file)))

  scenario_meta <- if (file_ext == "csv") {
    utils::read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (file_ext == "rds") {
    readRDS(file)
  } else {
    stop("Scenario file must be a .csv or .rds file: ", file)
  }

  if (!is.data.frame(scenario_meta)) {
    stop("Scenario file must contain a data.frame.")
  }

  if (!"Scenario" %in% names(scenario_meta)) {
    scenario_meta$Scenario <- seq_len(nrow(scenario_meta))
  }

  if (!"True_MTD" %in% names(scenario_meta) && "MTD" %in% names(scenario_meta)) {
    scenario_meta$True_MTD <- scenario_meta$MTD
  }

  if (!"Source_Scenario" %in% names(scenario_meta)) {
    scenario_meta$Source_Scenario <- scenario_meta$Scenario
  }

  if (!"Scenario_Group" %in% names(scenario_meta)) {
    scenario_meta$Scenario_Group <- "scenario_file"
  }

  if (!"Attempt" %in% names(scenario_meta)) {
    scenario_meta$Attempt <- NA_integer_
  }

  required_meta_cols <- c("Scenario", "Source_Scenario", "Scenario_Group", "True_MTD", "Attempt")
  missing_meta_cols <- setdiff(required_meta_cols, names(scenario_meta))
  if (length(missing_meta_cols) > 0L) {
    stop(
      "Scenario file is missing required columns: ",
      paste(missing_meta_cols, collapse = ", ")
    )
  }

  integer_cols <- intersect(
    c("Scenario", "Source_Scenario", "True_MTD", "Attempt"),
    names(scenario_meta)
  )
  scenario_meta[integer_cols] <- lapply(scenario_meta[integer_cols], as.integer)

  if (anyNA(scenario_meta$Scenario) || anyDuplicated(scenario_meta$Scenario)) {
    stop("Scenario file must have non-missing unique Scenario values.")
  }

  scenario_meta
}

make_aide_folder <- function(task) {
  method_tag <- if (task$model == "BOIN") {
    paste0(task$decision_method, "-", task$r_estimator)
  } else if (task$model == "CRM") {
    paste0("crm_", task$crm_r_model)
  } else {
    paste0("cfo_", task$cfo_method)
  }

  ## Keep folder names short to avoid filesystem path limits on the cluster.
  ## Per-setting details are still recorded in each task log and saved RDS object.
  scenario_set <- if (!is.null(task$scenario_set)) task$scenario_set else "default"

  paste0(
    scenario_set,
    "-model-", task$model,
    "-opt-", method_tag,
    "-w-", fmt_short(task$T_assess),
    "-c-", fmt_short(task$C),
    "-cyc-", fmt_short(task$cycle_max),
    "-rate-", fmt_short(task$arrival_rate),
    "-Nmax-", fmt_short(task$Nmax_eff),
    "-tried-", as.integer(isTRUE(task$restrict_to_tried)),
    if (isTRUE(task$restrict_to_target)) "-ptarget-1" else ""
  )
}

combine_oc_AIDE_results <- function(files) {
  x <- lapply(files, readRDS)

  ntrial_total <- sum(vapply(x, function(z) z$ntrial, numeric(1)))
  ndose <- x[[1]]$ndose

  final_MTD <- unlist(lapply(x, `[[`, "final_MTD"), use.names = FALSE)

  bind_trial_matrix <- function(field) {
    mats <- lapply(x, function(z) z[[field]])
    mats <- mats[!vapply(mats, is.null, logical(1))]

    if (length(mats) == 0L) {
      return(NULL)
    }

    do.call(rbind, mats)
  }

  mean_by_col <- function(mat) {
    if (is.null(mat)) {
      return(rep(NA_real_, ndose))
    }

    apply(mat, 2, function(z) {
      if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE)
    })
  }

  pj_iso_by_trial <- bind_trial_matrix("pj_iso_by_trial")
  r_hat_by_trial <- bind_trial_matrix("r_hat_by_trial")
  r_use_by_trial <- bind_trial_matrix("r_use_by_trial")
  r_cap_by_trial <- bind_trial_matrix("r_cap_by_trial")

  total_admin <- unlist(lapply(x, `[[`, "total_admin"), use.names = FALSE)
  total_unique <- unlist(lapply(x, `[[`, "total_unique"), use.names = FALSE)
  duration <- unlist(lapply(x, `[[`, "duration"), use.names = FALSE)

  earlystop <- unlist(lapply(x, function(z) {
    if (is.null(z$earlystop)) {
      rep(NA_integer_, z$ntrial)
    } else {
      z$earlystop
    }
  }), use.names = FALSE)

  earlystop_count <- sum(earlystop, na.rm = TRUE)

  sel_count <- tabulate(
    final_MTD[!is.na(final_MTD) & final_MTD >= 1L & final_MTD <= ndose],
    nbins = ndose
  )

  stop_count <- sum(final_MTD == 99L, na.rm = TRUE)

  na_count <- sum(
    is.na(final_MTD) |
      (!is.na(final_MTD) & final_MTD != 99L &
         (final_MTD < 1L | final_MTD > ndose))
  )

  weighted_mean_vec <- function(field) {
    Reduce(
      "+",
      lapply(x, function(z) z[[field]] * z$ntrial)
    ) / ntrial_total
  }

  out <- x[[1]]

  out$ntrial <- ntrial_total
  out$sel_count <- sel_count
  out$stop_count <- stop_count
  out$na_count <- na_count

  out$selection_pct <- 100 * sel_count / ntrial_total
  out$early_stop_pct <- 100 * earlystop_count / ntrial_total
  out$na_pct <- 100 * na_count / ntrial_total

  out$earlystop <- earlystop
  out$earlystop_count <- earlystop_count

  out$final_MTD <- final_MTD

  out$pj_iso_by_trial <- pj_iso_by_trial
  out$pj_iso_mean <- mean_by_col(pj_iso_by_trial)

  out$r_hat_by_trial <- r_hat_by_trial
  out$r_hat_mean <- mean_by_col(r_hat_by_trial)

  out$r_use_by_trial <- r_use_by_trial
  out$r_use_mean <- mean_by_col(r_use_by_trial)

  out$r_cap_by_trial <- r_cap_by_trial
  out$r_cap_mean <- mean_by_col(r_cap_by_trial)

  out$n_by_dose <- weighted_mean_vec("n_by_dose")
  out$unique_n_by_dose <- weighted_mean_vec("unique_n_by_dose")
  out$nipde_by_dose <- weighted_mean_vec("nipde_by_dose")

  out$total_admin <- total_admin
  out$total_unique <- total_unique
  out$duration <- duration

  out$total_admin_mean <- mean(total_admin, na.rm = TRUE)
  out$total_unique_mean <- mean(total_unique, na.rm = TRUE)
  out$duration_mean <- mean(duration, na.rm = TRUE)

  out$raw_trials <- NULL

  class(out) <- "oc_AIDE"
  out
}

## ============================================================
## User settings
## ============================================================

target_BOIN <- 0.30
cutoff_equiv <- 0.95

## Switch this between 5L and 8L to change both scenarios and priors.
scenario_dose_count <- 5L
if (!scenario_dose_count %in% c(5L, 8L)) {
  stop("scenario_dose_count must be 5L or 8L.")
}

## Set TRUE to use a scenario CSV generated by random_scenarios.R.
## The file name is built from these parameters and must match the
## random_scenarios.R settings used to generate the CSV.
use_random_scenario_file <- TRUE
random_nscenario <- 10000L
random_target <- target_BOIN
random_target_diff_below <- 0.05
random_target_diff_above <- 0.05
random_scenario_batch_size <- 10L
random_scenario_batch_i <- 1L

args <- commandArgs(trailingOnly = TRUE)

## This is trials per LSF job, not total across all 2000 jobs.
ntrial.total <- 1L

seed_base <- 1L

T_assess_equiv <- 28
C_equiv <- 3L

## Maximum number of cycles/dosings per patient to evaluate.
## Set to c(1L, 2L, 3L) to compare no-IPDE vs 2-dose vs 3-dose IPDE.
cycle_max_list_equiv <- c(1L, 2L, 3L)

Nmax_eff_list_equiv <- c(30L)
ipde_design_equiv <- 2L

t0 <- 0
day_obs <- 0

## Use continuous enrollment with a backlog queue, matching the AIDE
## illustration more closely than just-in-time recruitment.
continuous_enrollment_equiv <- TRUE
dlt_dist_equiv <- 2
dlt_alpha_equiv <- 0.5

d_cap_aide <- 100
dose_cap_aide <- 3L

store_raw <- FALSE
verbose <- FALSE

## Choose AIDE models here.
## Backfill standard BOIN/CRM for selected scenarios.
model_list_aide <- c("BOIN", "CRM")
## model_list_aide <- c("BOIN")
## model_list_aide <- c("BOIN", "CRM", "CFO")
## model_list_aide <- c("CRM")
## model_list_aide <- c("CFO")

## BOIN methods and r estimator.
method_list_boin <- c("approx1", "approx2")
## method_list_boin <- c("boin", "approx1", "approx2")

## BOIN r options:
##   r_fixed    : use fixed r_carry
##   r_mle      : local smoothed MLE truncated by adjacent regular toxicity
##   r_adaptive : global posterior estimator of r, integrating p_k out
r_estimator_list_boin <- c("r_fixed")
## r_estimator_list_boin <- c("r_adaptive")
## r_estimator_list_boin <- c("r_fixed", "r_mle", "r_adaptive")

## CRM versions.
## Supported AIDE-CRM methods:
##   fixed/r_fixed : discount CRM with fixed r
##   random        : discount CRM with one random r
##   level         : level-specific random r, available only for 5-dose runs
##   alpha_crm     : alpha-CRM effective-dose model
##   cumu_crm      : logistic cumulative-dose CRM/IPCRM
crm_r_model_list <- if (scenario_dose_count == 5L) {
  c("fixed", "random", "level", "alpha_crm", "cumu_crm")
} else {
  c("fixed", "random", "alpha_crm", "cumu_crm")
}
## For quick tests, use for example:
## crm_r_model_list <- c("alpha_crm")
## crm_r_model_list <- c("cumu_crm")
## crm_r_model_list <- c("fixed")

## Final MTD selection gate.
## TRUE: select among tried/non-eliminated doses.
## FALSE: allow all non-eliminated doses to enter final selection.
restrict_to_tried_list_aide <- c(TRUE)

## Final MTD target-safety gate.
## TRUE: final selected MTD must have estimated p_j <= target.
## FALSE: select the admissible dose closest to target, as before.
restrict_to_target_list_aide <- c(FALSE)

## CRM settings from methods_prior.R.
## Power CRM / alpha-CRM prior: theta ~ N(0, 2).
theta_mean <- 0
theta_sd <- sqrt(2)

if (scenario_dose_count == 5L) {
  ## 5-dose prior settings.
  q_skeleton <- c(0.15, 0.20, 0.30, 0.35, 0.45)
  dose_alpha_mg <- c(15, 20, 30, 35, 45)
  dose_ipcrm_raw <- c(15, 20, 30, 35, 45)
  fixed_intercept <- -2.8
  beta1_shape <- 2.5
  cfo_skeleton_default <- c(0.005, 0.01, 0.05, 0.10, 0.30)
} else if (scenario_dose_count == 8L) {
  ## 8-dose prior settings.
  q_skeleton <- c(0.12, 0.16, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45)
  dose_alpha_mg <- c(10, 20, 30, 40, 50, 60, 70, 80)
  dose_ipcrm_raw <- c(10, 20, 30, 40, 50, 60, 70, 80)
  fixed_intercept <- -2.2
  beta1_shape <- 2.2
  cfo_skeleton_default <- c(0.002, 0.008, 0.012, 0.04, 0.08, 0.10, 0.20, 0.35)
}

crm_skeleton_default <- q_skeleton
crm_alpha_sd_default <- theta_sd
crm_a_r_default <- target_BOIN / 2
crm_b_r_default <- 1 - crm_a_r_default

## New CRM helper uses one uniform model_file argument.
## Keep NULL to use defaults:
##   fixed  -> fix_CRM.bug
##   random -> random_CRM.bug
##   level  -> random_CRM_level.bug
crm_model_file_default <- NULL
crm_fixed_model_file_default <- "fix_CRM.bug"
crm_random_model_file_default <- "random_CRM.bug"
crm_level_model_file_default <- "random_CRM_level.bug"

## alpha-CRM settings.
## Baseline model: p_j = S(d_j)^exp(theta), with S(d_j)=skeleton_j.
## For IPDE observations, effective dose uses actual dose amounts and
## calendar-time gaps based on crm_time_col_default.
crm_dose_values_alpha_default <- dose_alpha_mg
crm_time_col_default <- "t_start"
crm_alpha_grid_default <- seq(0.01, 0.99, length.out = 61)
crm_alpha_T_default <- T_assess_equiv
crm_theta_prior_mean_default <- theta_mean
crm_theta_prior_sd_default <- theta_sd
crm_alpha_L_default <- 8
crm_alpha_rel_tol_default <- 1e-8
crm_alpha_eps_default <- 1e-12
crm_alpha_n_draw_prior_default <- 5000

## Logistic cumulative CRM / IPCRM settings.
## beta0 ~ t(fixed_intercept, precision = beta0_prec, df = beta0_df)
## beta1 ~ Gamma(beta1_shape, beta1_rate)
## beta2 ~ Exp(1), but beta2 drops out at baseline because cumu.d = 0.
beta0_prec <- 2
beta0_df <- 1
beta1_rate <- 1.6

## IPCRM / cumulative CRM current dose scores.
## These values are passed directly to the cumulative CRM helper.
dose_ipcrm <- dose_ipcrm_raw / (2 * stats::sd(dose_ipcrm_raw))

crm_dose_scores_cumu_default <- dose_ipcrm
crm_cumu_model_file_default <- NULL
crm_cumu_beta0_mean_default <- fixed_intercept
crm_cumu_beta0_prec_default <- beta0_prec
crm_cumu_beta0_df_default <- beta0_df
crm_cumu_beta1_shape_default <- beta1_shape
crm_cumu_beta1_rate_default <- beta1_rate
crm_cumu_beta2_rate_default <- 1
crm_cumu_include_current_default <- FALSE

## CFO / PRIDE settings from methods_prior.R.
cfo_method_list_aide <- c("empirical", "pride")
cfo_model_file_default <- "PRIDE.bug"
cfo_sigma2_beta_default <- 30
cfo_eta_default <- 1
cfo_pk_method_default <- "approx"
cfo_n_mc_w_default <- 200
cfo_m_use_default <- 1000
cfo_use_monotone_pair_default <- FALSE

crm_n_chains <- 2
crm_n_adapt <- 500
crm_n_burnin <- 500
crm_n_iter <- 2000
crm_thin <- 1

cfo_n_chains <- 3
cfo_n_adapt <- 1000
cfo_n_burnin <- 2000
cfo_n_iter <- 5000
cfo_thin <- 2

alpha_true_list <- c(0, 0.3, 0.6, 0.9)
accrual_list <- c(1 / 14)
## For BOIN with r_estimator = "r_mle", r_carry is not used.
## Keep one placeholder value to avoid duplicate identical BOIN runs.
r_carry_list_boin <- c(0)

## For CRM or BOIN with r_estimator = "r_fixed", r_carry is used.
r_carry_list_crm <- c(0)
r_carry_list_fixed_boin <- c(0)

## Adaptive/global r estimator settings.  These match the carryover
## estimator demo and use posterior mean r as the plug-in value.
r_adaptive_prior_default <- c(target_BOIN / 2, 1 - target_BOIN / 2)
p_adaptive_prior_default <- c(target_BOIN / 2, 1 - target_BOIN / 2)
r_adaptive_max_default <- 0.99
r_adaptive_plug_in_default <- "mean"
r_adaptive_rel_tol_default <- 1e-6

## Scenario set.
scenario_dir <- "scenario_sets"

if (isTRUE(use_random_scenario_file)) {
  scenario_file <- make_random_scenario_filename(
    ndose = scenario_dose_count,
    nscenario = random_nscenario,
    target = random_target,
    target_diff_below = random_target_diff_below,
    target_diff_above = random_target_diff_above,
    out_dir = scenario_dir
  )
  scenario_set_name <- tools::file_path_sans_ext(basename(scenario_file))

  if (!identical(as.numeric(random_target), as.numeric(target_BOIN))) {
    warning("random_target differs from target_BOIN; run_oc_AIDE will still use target_BOIN for design decisions.")
  }

  random_scenario_batch_i <- get_scenario_batch_index(
    args,
    default = random_scenario_batch_i
  )

  if (is.na(random_scenario_batch_size) || random_scenario_batch_size < 1L) {
    stop("random_scenario_batch_size must be a positive integer.")
  }

  random_scenario_start <- (random_scenario_batch_i - 1L) * random_scenario_batch_size + 1L
  random_scenario_end <- min(
    random_scenario_batch_i * random_scenario_batch_size,
    as.integer(random_nscenario)
  )

  if (random_scenario_start > as.integer(random_nscenario)) {
    stop(
      "Random scenario batch ",
      random_scenario_batch_i,
      " starts after random_nscenario = ",
      random_nscenario,
      "."
    )
  }

  scenarios <- random_scenario_start:random_scenario_end
} else if (scenario_dose_count == 5L) {
  scenario_set_name <- "Set_5dose_adaptive_r_37"
  scenario_file <- file.path(scenario_dir, paste0(scenario_set_name, ".csv"))

  ## Default 5-dose run: newly added high-gap scenarios.
  scenarios <- 26:37
} else if (scenario_dose_count == 8L) {
  scenario_set_name <- "Set_8dose_adaptive_r_40"
  scenario_file <- file.path(scenario_dir, paste0(scenario_set_name, ".csv"))

  ## Default 8-dose run: newly added high-gap scenarios.
  scenarios <- 21:40
}

scenario_meta <- load_scenario_file(scenario_file)

dose_col_names <- paste0("Dose", seq_along(crm_skeleton_default))
if (!all(dose_col_names %in% names(scenario_meta))) {
  stop("scenario_meta is missing required dose columns: ",
       paste(setdiff(dose_col_names, names(scenario_meta)), collapse = ", "))
}

scenario_matrix <- as.matrix(scenario_meta[dose_col_names])
rownames(scenario_matrix) <- as.character(scenario_meta$Scenario)

## Choose which scenario indexes to run by overriding scenarios above.
## scenarios <- seq_len(nrow(scenario_matrix))

scenario_id_list <- as.integer(scenarios)
if (anyNA(scenario_id_list) || any(!scenario_id_list %in% scenario_meta$Scenario)) {
  stop("All entries in scenarios must match scenario_meta$Scenario.")
}

if (ncol(scenario_matrix) != length(crm_skeleton_default)) {
  stop("Scenarios and CRM skeleton have different lengths.")
}

if (length(crm_dose_values_alpha_default) != length(crm_skeleton_default)) {
  stop("alpha-CRM dose values and CRM skeleton have different lengths.")
}

if (length(crm_dose_scores_cumu_default) != length(crm_skeleton_default)) {
  stop("IPCRM dose scores and CRM skeleton have different lengths.")
}

if (length(cfo_skeleton_default) != length(crm_skeleton_default)) {
  stop("CFO skeleton and CRM skeleton have different lengths.")
}

## ============================================================
## Command-line arguments
##
## args[1] = worker count, retained for backward compatibility and ignored
## args[2] = number of trials per LSF job
## args[3] = job index, 1,...,2000
## args[4] = random scenario batch index, if using random scenario file
## ============================================================

if (length(args) >= 1) {
  requested_workers <- suppressWarnings(as.integer(args[1]))
} else {
  requested_workers <- 1L
}

if (length(args) >= 2) {
  ntrial.total <- as.integer(args[2])
}

job_i <- get_job_index(args)

if (is.na(job_i) || job_i < 1L) {
  stop("job_i must be a positive integer.")
}

if (is.na(requested_workers) || requested_workers < 1L) {
  stop("args[1] worker count must be a positive integer when supplied.")
}

if (is.na(ntrial.total) || ntrial.total < 1L) {
  stop("ntrial.total must be a positive integer.")
}

if (requested_workers != 1L) {
  warning("Worker-count argument ignored; running sequentially with one worker.")
}

cat("LSF job/block index:", job_i, "\n")
cat("Working directory:", getwd(), "\n")
cat("Execution mode: sequential\n")
cat("Requested workers argument:", requested_workers, "\n")
cat("Trials per setting in this job:", ntrial.total, "\n")
cat("Models:", paste(model_list_aide, collapse = ", "), "\n")
cat("BOIN methods:", paste(method_list_boin, collapse = ", "), "\n")
cat("BOIN r estimators:", paste(r_estimator_list_boin, collapse = ", "), "\n")
cat("BOIN adaptive r prior:", paste(r_adaptive_prior_default, collapse = ", "), "\n")
cat("BOIN adaptive p prior:", paste(p_adaptive_prior_default, collapse = ", "), "\n")
cat("BOIN adaptive r max:", r_adaptive_max_default, "\n")
cat("BOIN adaptive r plug-in:", r_adaptive_plug_in_default, "\n")
cat("CRM r models:", paste(crm_r_model_list, collapse = ", "), "\n")
cat("CFO methods:", paste(cfo_method_list_aide, collapse = ", "), "\n")
cat("Nmax_eff list:", paste(Nmax_eff_list_equiv, collapse = ", "), "\n")
cat("restrict_to_tried list:", paste(restrict_to_tried_list_aide, collapse = ", "), "\n")
cat("restrict_to_target list:", paste(restrict_to_target_list_aide, collapse = ", "), "\n")
cat("Continuous enrollment:", continuous_enrollment_equiv, "\n")
cat("Cycle max list:", paste(cycle_max_list_equiv, collapse = ", "), "\n")
cat("Scenario dose count:", scenario_dose_count, "\n")
cat("CRM skeleton:", paste(crm_skeleton_default, collapse = ", "), "\n")
cat("alpha-CRM dose values:", paste(crm_dose_values_alpha_default, collapse = ", "), "\n")
cat("alpha-CRM theta prior mean/sd:", crm_theta_prior_mean_default, crm_theta_prior_sd_default, "\n")
cat("cumu CRM dose scores:", paste(round(crm_dose_scores_cumu_default, 6), collapse = ", "), "\n")
cat("cumu CRM beta0 mean/prec/df:", crm_cumu_beta0_mean_default, crm_cumu_beta0_prec_default, crm_cumu_beta0_df_default, "\n")
cat("cumu CRM beta1 shape/rate:", crm_cumu_beta1_shape_default, crm_cumu_beta1_rate_default, "\n")
cat("CFO skeleton:", paste(cfo_skeleton_default, collapse = ", "), "\n")
cat("CFO sigma2_beta / eta:", cfo_sigma2_beta_default, cfo_eta_default, "\n")
cat("CFO model file:", cfo_model_file_default, "\n")
cat("Scenario set:", scenario_set_name, "\n")
cat("Scenario file:", scenario_file, "\n")
if (isTRUE(use_random_scenario_file)) {
  cat("Random scenario batch index:", random_scenario_batch_i, "\n")
  cat("Random scenario batch size:", random_scenario_batch_size, "\n")
  cat("Random scenario range:", random_scenario_start, "to", random_scenario_end, "\n")
}
cat("Scenario IDs:", paste(scenario_id_list, collapse = ", "), "\n")

trial_blocks <- ntrial.total
block_start <- 1L

cat("Trial blocks:", paste(trial_blocks, collapse = ", "), "\n")
cat("Block starts:", paste(block_start, collapse = ", "), "\n")

## ============================================================
## Worker function
## ============================================================

run_one_aide_task <- function(task) {
  setwd(task$workdir)

  if (!is.null(task$lib_path) && dir.exists(task$lib_path)) {
    .libPaths(c(task$lib_path, .libPaths()))
  }

  source("AIDE_BOIN_helper.R")
  source("AIDE_CRM_helper_final.R")
  source("AIDE_modified.R")

  if (requireNamespace("rjags", quietly = TRUE)) library(rjags)
  if (requireNamespace("coda", quietly = TRUE)) library(coda)

  foldername <- make_aide_folder(task)

  outdir <- file.path(task$outdir_root, foldername)
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  method_tag <- if (task$model == "BOIN") {
    paste0(task$decision_method, "-", task$r_estimator)
  } else if (task$model == "CRM") {
    paste0("crm_", task$crm_r_model)
  } else {
    paste0("cfo_", task$cfo_method)
  }

  scenario_name <- paste0(task$scenario_set, "_SC", task$scenario_id)

  filename <- paste0(
    scenario_name,
    "-", task$model,
    "-", method_tag,
    "-a", fmt_short(task$alpha_true),
    "-r", fmt_short(task$r_carry),
    "-cyc", fmt_short(task$cycle_max),
    "-tried", as.integer(isTRUE(task$restrict_to_tried)),
    if (isTRUE(task$restrict_to_target)) "-ptarget1" else "",
    "-j", task$job_i,
    "-b", task$block_id,
    "-s", task$seed.block,
    "-n", task$ntrial.block,
    ".rds"
  )

  outfile <- file.path(outdir, filename)

  logfile <- file.path(
    outdir,
    paste0(
      scenario_name,
      "-", task$model,
      "-", method_tag,
      "-cyc", fmt_short(task$cycle_max),
      "-tried", as.integer(isTRUE(task$restrict_to_tried)),
      if (isTRUE(task$restrict_to_target)) "-ptarget1" else "",
      "-j", task$job_i,
      "-b", task$block_id,
      ".log"
    )
  )

  with_sink(logfile, {
    cat("====================================\n")
    cat("AIDE cluster task\n")
    cat("====================================\n")
    cat("Scenario set:", task$scenario_set, "\n")
    cat("Scenario:", task$scenario_id, "\n")
    cat("Scenario group:", task$scenario_group, "\n")
    cat("Model:", task$model, "\n")
    cat("Method:", method_tag, "\n")
    cat("Job:", task$job_i, "\n")
    cat("Block:", task$block_id, "\n")
    cat("ntrial.block:", task$ntrial.block, "\n")
    cat("Seed:", task$seed.block, "\n")
    cat("Target:", task$target, "\n")
    cat("alpha_true:", task$alpha_true, "\n")
    cat("r_carry:", task$r_carry, "\n")
    cat("Accrual:", task$arrival_rate, "\n")
    cat("T_assess:", task$T_assess, "\n")
    cat("C:", task$C, "\n")
    cat("cycle_max:", task$cycle_max, "\n")
    cat("Nmax_eff:", task$Nmax_eff, "\n")
    cat("dose_cap:", task$dose_cap, "\n")
    cat("restrict_to_tried:", task$restrict_to_tried, "\n")
    cat("restrict_to_target:", task$restrict_to_target, "\n")
    cat("continuous_enrollment:", task$continuous_enrollment, "\n")
    cat("p.true:", paste(task$p.true, collapse = ", "), "\n")
    cat("p.true_ipde:", paste(task$p.true_ipde, collapse = ", "), "\n")

    if (task$model == "BOIN") {
      cat("BOIN r estimator:", task$r_estimator, "\n")
      if (identical(task$r_estimator, "r_adaptive")) {
        cat("BOIN adaptive r prior:", paste(task$r_adaptive_prior, collapse = ", "), "\n")
        cat("BOIN adaptive p prior:", paste(task$p_adaptive_prior, collapse = ", "), "\n")
        cat("BOIN adaptive r max:", task$r_adaptive_max, "\n")
        cat("BOIN adaptive r plug-in:", task$r_adaptive_plug_in, "\n")
        cat("BOIN adaptive rel.tol:", task$r_adaptive_rel_tol, "\n")
      }
    }

    if (task$model == "CRM") {
      cat("CRM r model:", task$crm_r_model, "\n")
      cat("CRM skeleton:", paste(task$crm_skeleton, collapse = ", "), "\n")
      cat("CRM alpha sd:", task$crm_alpha_sd, "\n")
      cat("CRM a_r:", task$crm_a_r, "\n")
      cat("CRM b_r:", task$crm_b_r, "\n")
      cat("CRM model file:", task$crm_model_file, "\n")
      cat("CRM fixed model file:", task$crm_fixed_model_file, "\n")
      cat("CRM random model file:", task$crm_random_model_file, "\n")
      cat("CRM level model file:", task$crm_level_model_file, "\n")
      if (task$crm_r_model == "alpha_crm") {
        cat("alpha-CRM dose values:", paste(task$crm_dose_values, collapse = ", "), "\n")
        cat("alpha-CRM time col:", task$crm_time_col, "\n")
        cat("alpha-CRM alpha_T:", task$crm_alpha_T, "\n")
        cat("alpha-CRM theta prior mean/sd:", task$crm_theta_prior_mean, task$crm_theta_prior_sd, "\n")
        cat("alpha-CRM alpha grid length:", length(task$crm_alpha_grid), "\n")
      }
      if (task$crm_r_model == "cumu_crm") {
        cat("cumu CRM dose scores:", paste(task$crm_dose_scores, collapse = ", "), "\n")
        cat("cumu CRM beta0 mean/prec/df:", task$crm_cumu_beta0_mean, task$crm_cumu_beta0_prec, task$crm_cumu_beta0_df, "\n")
        cat("cumu CRM beta1 shape/rate:", task$crm_cumu_beta1_shape, task$crm_cumu_beta1_rate, "\n")
        cat("cumu CRM beta2 rate:", task$crm_cumu_beta2_rate, "\n")
        cat("cumu CRM include current:", task$crm_cumu_include_current, "\n")
      }
    }

    if (task$model == "CFO") {
      cat("CFO method:", task$cfo_method, "\n")
      cat("CFO skeleton:", paste(task$cfo_skeleton, collapse = ", "), "\n")
      cat("CFO model file:", task$cfo_model_file, "\n")
      cat("CFO sigma2_beta:", task$cfo_sigma2_beta, "\n")
      cat("CFO eta:", task$cfo_eta, "\n")
      cat("CFO pk method:", task$cfo_pk_method, "\n")
      cat("CFO n_mc_w:", task$cfo_n_mc_w, "\n")
      cat("CFO m_use:", task$cfo_m_use, "\n")
      cat("CFO use monotone pair:", task$cfo_use_monotone_pair, "\n")
    }

    cat("====================================\n\n")

    res <- get_oc_sim_AIDE(
      target = task$target,
      p.true = task$p.true,
      p.true_ipde = task$p.true_ipde,

      ntrial = task$ntrial.block,
      seed = task$seed.block,

      model = task$model,
      ipde_design = task$ipde_design,

      N_pat = task$N_pat,
      Nmax_eff = task$Nmax_eff,
      C = task$C,
      T_assess = task$T_assess,
      cycle_max = task$cycle_max,

      arrival_rate = task$arrival_rate,
      t0 = task$t0,
      continuous_enrollment = task$continuous_enrollment,

      cutoff = task$cutoff,

      d.cap = task$d.cap,
      dose_cap = task$dose_cap,
      day_obs = task$day_obs,

      dlt_dist = task$dlt_dist,
      dlt_alpha = task$dlt_alpha,

      decision_method = task$decision_method,
      mtd_method = task$mtd_method,
      restrict_to_tried = task$restrict_to_tried,
      restrict_to_target = task$restrict_to_target,
      r_carry = task$r_carry,
      r_estimator = task$r_estimator,
      r_adaptive_prior = task$r_adaptive_prior,
      p_adaptive_prior = task$p_adaptive_prior,
      r_adaptive_max = task$r_adaptive_max,
      r_adaptive_plug_in = task$r_adaptive_plug_in,
      r_adaptive_rel_tol = task$r_adaptive_rel_tol,

      crm_r_model = task$crm_r_model,
      crm_skeleton = task$crm_skeleton,
      crm_alpha_sd = task$crm_alpha_sd,
      crm_a_r = task$crm_a_r,
      crm_b_r = task$crm_b_r,
      crm_model_file = task$crm_model_file,
      crm_fixed_model_file = task$crm_fixed_model_file,
      crm_random_model_file = task$crm_random_model_file,
      crm_level_model_file = task$crm_level_model_file,

      crm_dose_values = task$crm_dose_values,
      crm_time_col = task$crm_time_col,
      crm_alpha_grid = task$crm_alpha_grid,
      crm_alpha_T = task$crm_alpha_T,
      crm_theta_prior_mean = task$crm_theta_prior_mean,
      crm_theta_prior_sd = task$crm_theta_prior_sd,
      crm_alpha_L = task$crm_alpha_L,
      crm_alpha_rel_tol = task$crm_alpha_rel_tol,
      crm_alpha_eps = task$crm_alpha_eps,
      crm_alpha_n_draw_prior = task$crm_alpha_n_draw_prior,

      crm_dose_scores = task$crm_dose_scores,
      crm_cumu_model_file = task$crm_cumu_model_file,
      crm_cumu_beta0_mean = task$crm_cumu_beta0_mean,
      crm_cumu_beta0_prec = task$crm_cumu_beta0_prec,
      crm_cumu_beta0_df = task$crm_cumu_beta0_df,
      crm_cumu_beta1_shape = task$crm_cumu_beta1_shape,
      crm_cumu_beta1_rate = task$crm_cumu_beta1_rate,
      crm_cumu_beta2_rate = task$crm_cumu_beta2_rate,
      crm_cumu_include_current = task$crm_cumu_include_current,

      crm_n_chains = task$crm_n_chains,
      crm_n_adapt = task$crm_n_adapt,
      crm_n_burnin = task$crm_n_burnin,
      crm_n_iter = task$crm_n_iter,
      crm_thin = task$crm_thin,

      cfo_method = task$cfo_method,
      cfo_skeleton = task$cfo_skeleton,
      cfo_model_file = task$cfo_model_file,
      cfo_sigma2_beta = task$cfo_sigma2_beta,
      cfo_eta = task$cfo_eta,
      cfo_pk_method = task$cfo_pk_method,
      cfo_n_mc_w = task$cfo_n_mc_w,
      cfo_m_use = task$cfo_m_use,
      cfo_use_monotone_pair = task$cfo_use_monotone_pair,
      cfo_restrict_to_tried = task$cfo_restrict_to_tried,
      cfo_n_chains = task$cfo_n_chains,
      cfo_n_adapt = task$cfo_n_adapt,
      cfo_n_burnin = task$cfo_n_burnin,
      cfo_n_iter = task$cfo_n_iter,
      cfo_thin = task$cfo_thin,

      store_raw = task$store_raw,
      verbose = task$verbose
    )

    res$scenario_set <- task$scenario_set
    res$scenario_name <- scenario_name
    res$scenario_id <- task$scenario_id
    res$source_scenario <- task$source_scenario
    res$true_mtd <- task$true_mtd
    res$scenario_attempt <- task$scenario_attempt
    res$Nmax_eff <- task$Nmax_eff

    saveRDS(res, outfile)

    cat("\nSaved:", outfile, "\n")
  })

  list(
    file = outfile,
    folder = outdir,
    scenario_set = task$scenario_set,
    scenario_id = task$scenario_id,
    scenario_name = scenario_name,
    source_scenario = task$source_scenario,
    scenario_group = task$scenario_group,
    true_mtd = task$true_mtd,
    scenario_attempt = task$scenario_attempt,
    Nmax_eff = task$Nmax_eff,
    model = task$model,
    decision_method = task$decision_method,
    mtd_method = task$mtd_method,
    r_estimator = task$r_estimator,
    crm_r_model = task$crm_r_model,
    cfo_method = task$cfo_method,
    restrict_to_tried = task$restrict_to_tried,
    restrict_to_target = task$restrict_to_target,
    alpha_true = task$alpha_true,
    r_carry = task$r_carry,
    arrival_rate = task$arrival_rate,
    cycle_max = task$cycle_max,
    continuous_enrollment = task$continuous_enrollment,
    block_id = task$block_id,
    ntrial.block = task$ntrial.block,
    seed.block = task$seed.block
  )
}

## ============================================================
## Build task list
## ============================================================

workdir <- getwd()

## Save runs under a scenario-specific output root so old files are not overwritten.
outdir_root <- paste0("oc_results_cluster_AIDE_", scenario_set_name)
if (!dir.exists(outdir_root)) {
  dir.create(outdir_root, recursive = TRUE)
}

tasks <- list()
task_id <- 1L

## Non-overlapping seed range across job_i = 1,...,2000.
## If ntrial.total = 1, job 1 uses seed 1, job 2 uses seed 2, etc.
job_seed_offset <- seed_base + (job_i - 1L) * ntrial.total

for (Nmax_eff_aide in Nmax_eff_list_equiv) {
  for (restrict_to_tried_aide in restrict_to_tried_list_aide) {
    for (restrict_to_target_aide in restrict_to_target_list_aide) {
      for (model_aide in model_list_aide) {
        for (scenario_id in scenario_id_list) {
          p_base <- as.numeric(scenario_matrix[as.character(scenario_id), ])
          scenario_row <- match(scenario_id, scenario_meta$Scenario)

          for (alpha_true in alpha_true_list) {
            p_ipde <- make_p_ipde(p_base, alpha_true)

            for (accrual in accrual_list) {
              for (cycle_max_aide in cycle_max_list_equiv) {

          if (model_aide == "BOIN") {
            method_loop <- method_list_boin
            r_estimator_loop <- r_estimator_list_boin
            crm_loop <- NA_character_
            cfo_loop <- NA_character_
          } else if (model_aide == "CRM") {
            method_loop <- "crm"
            r_estimator_loop <- NA_character_
            crm_loop <- crm_r_model_list
            cfo_loop <- NA_character_
          } else {
            method_loop <- "cfo"
            r_estimator_loop <- NA_character_
            crm_loop <- NA_character_
            cfo_loop <- cfo_method_list_aide
          }

          for (decision_method_aide in method_loop) {
            for (r_estimator_aide in r_estimator_loop) {
              for (crm_r_model_aide in crm_loop) {
                for (cfo_method_aide in cfo_loop) {

                  if (model_aide == "BOIN") {
                    mtd_method_aide <- decision_method_aide
                    r_carry_loop <- if (identical(r_estimator_aide, "r_fixed")) {
                      r_carry_list_fixed_boin
                    } else {
                      r_carry_list_boin
                    }
                  } else if (model_aide == "CRM") {
                    mtd_method_aide <- paste0("crm_", crm_r_model_aide)
                    r_carry_loop <- r_carry_list_crm
                  } else {
                    mtd_method_aide <- paste0("cfo_", cfo_method_aide)
                    r_carry_loop <- c(0)
                  }

                  for (r_carry_aide in r_carry_loop) {
                    for (block_id in seq_along(trial_blocks)) {
                      ntrial.block <- trial_blocks[block_id]

                      ## Non-overlapping seed for each block and each job.
                      seed.block <- job_seed_offset + block_start[block_id] - 1L

                      tasks[[task_id]] <- list(
                        workdir = workdir,
                        lib_path = lib_path,
                        outdir_root = outdir_root,
                        scenario_set = scenario_set_name,

                        job_i = job_i,
                        block_id = block_id,
                        ntrial.block = ntrial.block,
                        seed.block = seed.block,

                        scenario_id = scenario_id,
                        source_scenario = scenario_meta$Source_Scenario[scenario_row],
                        scenario_group = scenario_meta$Scenario_Group[scenario_row],
                        p.true = p_base,
                        p.true_ipde = p_ipde,

                        target = target_BOIN,
                        cutoff = cutoff_equiv,

                        model = model_aide,
                        ipde_design = ipde_design_equiv,

                        true_mtd = scenario_meta$True_MTD[scenario_row],
                        scenario_attempt = scenario_meta$Attempt[scenario_row],

                        N_pat = Nmax_eff_aide,
                        Nmax_eff = Nmax_eff_aide,
                        C = C_equiv,
                        T_assess = T_assess_equiv,
                        cycle_max = cycle_max_aide,

                        arrival_rate = accrual,
                        t0 = t0,
                        continuous_enrollment = continuous_enrollment_equiv,

                        d.cap = d_cap_aide,
                        dose_cap = dose_cap_aide,
                        day_obs = day_obs,

                        dlt_dist = dlt_dist_equiv,
                        dlt_alpha = dlt_alpha_equiv,

                        decision_method = if (model_aide == "BOIN") decision_method_aide else "boin",
                        mtd_method = if (model_aide == "BOIN") mtd_method_aide else NULL,
                        restrict_to_tried = restrict_to_tried_aide,
                        restrict_to_target = restrict_to_target_aide,
                        r_carry = r_carry_aide,
                        r_estimator = if (model_aide == "BOIN") r_estimator_aide else "r_fixed",
                        r_adaptive_prior = r_adaptive_prior_default,
                        p_adaptive_prior = p_adaptive_prior_default,
                        r_adaptive_max = r_adaptive_max_default,
                        r_adaptive_plug_in = r_adaptive_plug_in_default,
                        r_adaptive_rel_tol = r_adaptive_rel_tol_default,

                        crm_r_model = if (model_aide == "CRM") crm_r_model_aide else "fixed",
                        crm_skeleton = crm_skeleton_default,
                        crm_alpha_sd = crm_alpha_sd_default,
                        crm_a_r = crm_a_r_default,
                        crm_b_r = crm_b_r_default,
                        crm_model_file = crm_model_file_default,
                        crm_fixed_model_file = crm_fixed_model_file_default,
                        crm_random_model_file = crm_random_model_file_default,
                        crm_level_model_file = crm_level_model_file_default,

                        ## alpha-CRM inputs
                        crm_dose_values = crm_dose_values_alpha_default,
                        crm_time_col = crm_time_col_default,
                        crm_alpha_grid = crm_alpha_grid_default,
                        crm_alpha_T = crm_alpha_T_default,
                        crm_theta_prior_mean = crm_theta_prior_mean_default,
                        crm_theta_prior_sd = crm_theta_prior_sd_default,
                        crm_alpha_L = crm_alpha_L_default,
                        crm_alpha_rel_tol = crm_alpha_rel_tol_default,
                        crm_alpha_eps = crm_alpha_eps_default,
                        crm_alpha_n_draw_prior = crm_alpha_n_draw_prior_default,

                        ## cumulative CRM / IPCRM inputs
                        crm_dose_scores = crm_dose_scores_cumu_default,
                        crm_cumu_model_file = crm_cumu_model_file_default,
                        crm_cumu_beta0_mean = crm_cumu_beta0_mean_default,
                        crm_cumu_beta0_prec = crm_cumu_beta0_prec_default,
                        crm_cumu_beta0_df = crm_cumu_beta0_df_default,
                        crm_cumu_beta1_shape = crm_cumu_beta1_shape_default,
                        crm_cumu_beta1_rate = crm_cumu_beta1_rate_default,
                        crm_cumu_beta2_rate = crm_cumu_beta2_rate_default,
                        crm_cumu_include_current = crm_cumu_include_current_default,

                        crm_n_chains = crm_n_chains,
                        crm_n_adapt = crm_n_adapt,
                        crm_n_burnin = crm_n_burnin,
                        crm_n_iter = crm_n_iter,
                        crm_thin = crm_thin,

                        ## CFO / PRIDE inputs
                        cfo_method = if (model_aide == "CFO") cfo_method_aide else "empirical",
                        cfo_skeleton = cfo_skeleton_default,
                        cfo_model_file = cfo_model_file_default,
                        cfo_sigma2_beta = cfo_sigma2_beta_default,
                        cfo_eta = cfo_eta_default,
                        cfo_pk_method = cfo_pk_method_default,
                        cfo_n_mc_w = cfo_n_mc_w_default,
                        cfo_m_use = cfo_m_use_default,
                        cfo_use_monotone_pair = cfo_use_monotone_pair_default,
                        cfo_restrict_to_tried = restrict_to_tried_aide,
                        cfo_n_chains = cfo_n_chains,
                        cfo_n_adapt = cfo_n_adapt,
                        cfo_n_burnin = cfo_n_burnin,
                        cfo_n_iter = cfo_n_iter,
                        cfo_thin = cfo_thin,

                        alpha_true = alpha_true,

                        store_raw = store_raw,
                        verbose = verbose
                      )

                      task_id <- task_id + 1L
              }
            }
          }
        }
      }
    }
  }
}
      }
    }
      }
    }
  }
}

cat("Number of tasks:", length(tasks), "\n")

## ============================================================
## Run jobs sequentially
## ============================================================

task_results <- lapply(tasks, run_one_aide_task)

cat("\nAll sequential tasks completed.\n")

## ============================================================
## Combine chunk results by setting within this LSF job
## ============================================================

group_key <- vapply(
  task_results,
  function(z) {
    method_tag <- if (z$model == "BOIN") {
      paste0(z$decision_method, "-", z$r_estimator)
    } else if (z$model == "CRM") {
      paste0("crm_", z$crm_r_model)
    } else {
      paste0("cfo_", z$cfo_method)
    }

    key_parts <- c(
      paste0(z$scenario_set, "_SC", z$scenario_id),
      z$model,
      method_tag,
      paste0("a", fmt_short(z$alpha_true)),
      paste0("r", fmt_short(z$r_carry)),
      paste0("rate", fmt_short(z$arrival_rate)),
      paste0("cyc", fmt_short(z$cycle_max)),
      paste0("Nmax", fmt_short(z$Nmax_eff)),
      paste0("cont", as.integer(isTRUE(z$continuous_enrollment))),
      paste0("tried", as.integer(isTRUE(z$restrict_to_tried)))
    )

    if (isTRUE(z$restrict_to_target)) {
      key_parts <- c(key_parts, "ptarget1")
    }

    paste(key_parts, collapse = "_")
  },
  character(1)
)

groups <- split(task_results, group_key)

for (g in names(groups)) {
  files <- vapply(groups[[g]], function(z) z$file, character(1))
  folder <- groups[[g]][[1]]$folder

  combined <- combine_oc_AIDE_results(files)

  combined_file <- file.path(
    folder,
    paste0(g, "-job-", job_i, "-combined.rds")
  )

  saveRDS(combined, combined_file)

  cat("Combined:", combined_file, "\n")
  cat("  ntrial.total =", combined$ntrial, "\n")
}

cat("\nAll AIDE cluster runs completed.\n")
