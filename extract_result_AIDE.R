## ============================================================
## Extract AIDE OC from cluster jobs
## Matched to run_oc_AIDE_cluster_SetB_8dose.R
##
## Newest runner naming rules:
##
## Folder:
##   oc_results_cluster_AIDE_SetB_8dose/
##   SetB_8dose-model-{model}-opt-{method_tag}-w-{T_assess}-c-{C}-cyc-{cycle_max}-
##   rate-{arrival_rate}-Nmax-{Nmax_eff}/
##
## Combined file:
##   SetB_8dose_SC{sc}_{model}_{method_tag}_a{alpha_true}_r{r_carry}_rate{arrival_rate}_
##   cyc{cycle_max}_cont{0/1}-job-{job}-combined.rds
##
## Supports:
##   BOIN: boin / approx1 / approx2 with r_fixed or r_mle
##   CRM : fixed / random / level / alpha_crm / cumu_crm
## ============================================================

rm(list = ls())

## -------------------------------
## Settings matched to newest runner
## -------------------------------

## setwd("/rsrch8/home/biostatistics/syang10/AIDE")

scenario_set_name <- "SetB_8dose"
results_root <- paste0("oc_results_cluster_AIDE_", scenario_set_name)

## Set B: eight-dose true DLT scenarios
p_true_setB <- rbind(
  `1` = c(0.07, 0.10, 0.12, 0.14, 0.17, 0.19, 0.22, 0.26),
  `2` = c(0.05, 0.08, 0.10, 0.15, 0.20, 0.22, 0.30, 0.45),
  `3` = c(0.05, 0.10, 0.12, 0.15, 0.20, 0.30, 0.50, 0.60),
  `4` = c(0.02, 0.10, 0.15, 0.18, 0.30, 0.44, 0.52, 0.60),
  `5` = c(0.05, 0.10, 0.22, 0.30, 0.46, 0.53, 0.59, 0.66),
  `6` = c(0.20, 0.30, 0.47, 0.51, 0.56, 0.60, 0.64, 0.69),
  `7` = c(0.30, 0.35, 0.50, 0.55, 0.60, 0.65, 0.70, 0.80),
  `8` = c(0.40, 0.41, 0.43, 0.45, 0.47, 0.49, 0.52, 0.56)
)
scenario_id_list <- as.integer(rownames(p_true_setB))
ndose_expected <- ncol(p_true_setB)

jobs.expected <- 1:2000

## If each LSF job used ntrial.total = 2, expected total is 4000.
## Change this if you pass a different trials-per-job argument to the runner.
ntrial_per_job_expected <- 1L
ntrial.expected <- length(jobs.expected) * ntrial_per_job_expected

target <- 0.30

T_assess <- 28
C <- 3L
cycle_max_list <- c(1L, 2L, 3L)
Nmax_eff <- 30L
dose_cap <- 3L
continuous_enrollment <- TRUE

alpha_true_list <- c(0, 0.3, 0.6, 0.9)
arrival_rate_list <- c(1 / 14)

## BOIN settings matched to runner
boin_method_list <- c("approx1")
boin_r_estimator_list <- c("r_fixed")
boin_r_carry_list_r_mle <- c(0)
boin_r_carry_list_fixed <- c(0)

## CRM settings matched to the actual SetB_8dose runner outputs.
## Your current example folder/file is for alpha_crm with r0:
##   SetB_8dose-model-CRM-opt-crm_alpha_crm-...
##   SetB_8dose_SC7_CRM_crm_alpha_crm_a0_r0_rate0p07_cyc1_cont1-job-1630-combined.rds
##
## Edit this vector to extract other CRM backends. The r values are selected
## below by crm_r_carry_values.
# crm_r_model_list <- c("alpha_crm")
crm_r_model_list <- c("fixed", "level", "random", "alpha_crm", "cumu_crm")

## r values by CRM backend, matching the names produced by the SetB runner.
## fixed CRM uses r_carry; alpha_crm/cumu_crm do not use fixed r, so the
## runner convention is usually r0 in the filename.
crm_r_carry_values <- list(
  fixed     = c(0),
  random    = c(0),
  level     = c(0),
  alpha_crm = c(0),
  cumu_crm  = c(0)
)

get_crm_r_loop <- function(crm_r_model) {
  out <- crm_r_carry_values[[crm_r_model]]
  if (is.null(out)) {
    stop("No r_carry values specified for CRM model: ", crm_r_model)
  }
  out
}

## skeleton for power CRM and alpha-CRM
crm_skeleton <- c(0.04, 0.08, 0.12, 0.16, 0.30, 0.47, 0.54, 0.60)
q_skeleton <- crm_skeleton

## Power CRM / alpha-CRM prior: theta ~ N(0, 2)
theta_mean <- 0
theta_sd <- sqrt(2)
crm_alpha_sd <- theta_sd
crm_theta_prior_mean <- theta_mean
crm_theta_prior_sd <- theta_sd
crm_a_r <- target / 2
crm_b_r <- 1 - crm_a_r

## alpha-CRM actual dose amounts, mg
crm_dose_values_alpha <- c(10, 20, 30, 40, 50, 60, 70, 80)
dose_alpha_mg <- crm_dose_values_alpha
crm_alpha_grid_length <- 61L

## IPCRM / cumulative CRM current dose scores and priors
## beta0 ~ t(fixed_intercept, precision = beta0_prec, df = beta0_df)
## beta1 ~ Gamma(beta1_shape, beta1_rate)
## beta2 ~ Exp(1), but beta2 drops out at baseline when cumu.d = 0
crm_dose_scores_raw <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5)
dose_ipcrm <- crm_dose_scores_raw
## Runner passes these scores directly; do not standardize them here.
crm_dose_scores_cumu <- dose_ipcrm
fixed_intercept <- -3
beta0_prec <- 2
beta0_df <- 1
beta1_shape <- 2.83
beta1_rate <- 1.21
beta2_rate <- 1
crm_cumu_beta0_mean <- fixed_intercept
crm_cumu_beta0_prec <- beta0_prec
crm_cumu_beta0_df <- beta0_df
crm_cumu_beta1_shape <- beta1_shape
crm_cumu_beta1_rate <- beta1_rate
crm_cumu_beta2_rate <- beta2_rate
crm_cumu_include_current <- FALSE

## Extract these models. For newest CRM-only runs, leave as c("CRM").
model_list <- c("CRM")
## model_list <- c("BOIN", "CRM")

out.dir <- paste0("OC_summary_from_parallel_AIDE_jobs_", scenario_set_name)
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

## -------------------------------
## Helpers
## -------------------------------

fmt_num <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return("NA")
  gsub("\\.", "p", as.character(x))
}

## Must match runner fmt_short exactly.
fmt_short <- function(x, digits = 2) {
  x <- round(as.numeric(x), digits)
  out <- format(x, scientific = FALSE, trim = TRUE)
  out <- sub("(\\.[0-9]*?)0+$", "\\1", out)
  out <- sub("\\.$", "", out)
  out[out == "-0"] <- "0"
  out[out == ""] <- "0"
  out <- gsub("-", "m", out)
  out <- gsub("\\.", "p", out)
  out
}

calc_select_rate_pct <- function(sel, ndose, denom) {
  out <- rep(0, ndose)
  if (denom <= 0) return(out)
  for (k in seq_len(ndose)) {
    out[k] <- 100 * sum(sel == k, na.rm = TRUE) / denom
  }
  out
}

get_field <- function(x, candidates, required = TRUE) {
  for (nm in candidates) {
    if (!is.null(x[[nm]])) return(x[[nm]])
  }
  if (required) {
    stop("None of these fields found: ", paste(candidates, collapse = ", "))
  }
  NULL
}

make_method_tag <- function(model,
                            boin_method = NULL,
                            boin_r_estimator = NULL,
                            crm_r_model = NULL) {
  if (model == "BOIN") {
    return(paste0(boin_method, "-", boin_r_estimator))
  }
  paste0("crm_", crm_r_model)
}

make_foldername <- function(model,
                            method_tag,
                            cycle_max,
                            arrival_rate,
                            scenario_set = scenario_set_name) {
  paste0(
    scenario_set,
    "-model-", model,
    "-opt-", method_tag,
    "-w-", fmt_short(T_assess),
    "-c-", fmt_short(C),
    "-cyc-", fmt_short(cycle_max),
    "-rate-", fmt_short(arrival_rate),
    "-Nmax-", fmt_short(Nmax_eff)
  )
}

make_group_key <- function(sc,
                           model,
                           method_tag,
                           alpha_true,
                           r_carry,
                           arrival_rate,
                           cycle_max,
                           continuous_enrollment,
                           scenario_set = scenario_set_name) {
  paste(
    paste0(scenario_set, "_SC", sc),
    model,
    method_tag,
    paste0("a", fmt_short(alpha_true)),
    paste0("r", fmt_short(r_carry)),
    paste0("rate", fmt_short(arrival_rate)),
    paste0("cyc", fmt_short(cycle_max)),
    paste0("cont", as.integer(isTRUE(continuous_enrollment))),
    sep = "_"
  )
}

make_filename <- function(sc,
                          model,
                          method_tag,
                          alpha_true,
                          r_carry,
                          arrival_rate,
                          cycle_max,
                          continuous_enrollment,
                          job) {
  paste0(
    make_group_key(
      sc = sc,
      model = model,
      method_tag = method_tag,
      alpha_true = alpha_true,
      r_carry = r_carry,
      arrival_rate = arrival_rate,
      cycle_max = cycle_max,
      continuous_enrollment = continuous_enrollment
    ),
    "-job-", job,
    "-combined.rds"
  )
}

bind_matrix_field <- function(res.list, field_candidates, required = FALSE) {
  mat.list <- lapply(res.list, function(x) {
    z <- get_field(x, field_candidates, required = FALSE)
    if (is.null(z)) return(NULL)
    as.matrix(z)
  })
  
  keep <- !vapply(mat.list, is.null, logical(1))
  mat.list <- mat.list[keep]
  
  if (length(mat.list) == 0L) {
    if (required) {
      stop("None of these matrix fields found: ", paste(field_candidates, collapse = ", "))
    }
    return(NULL)
  }
  
  do.call(rbind, mat.list)
}

mean_from_matrix <- function(mat, ndose) {
  if (is.null(mat)) return(rep(NA_real_, ndose))
  apply(mat, 2, function(z) {
    if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE)
  })
}

extract_mean_vector <- function(res.list,
                                by_trial_candidates = NULL,
                                mean_candidates = NULL,
                                ndose = NULL,
                                required = FALSE) {
  if (is.null(ndose)) {
    ndose <- get_field(res.list[[1]], c("ndose"))
  }
  
  if (!is.null(by_trial_candidates)) {
    mat <- bind_matrix_field(
      res.list,
      field_candidates = by_trial_candidates,
      required = FALSE
    )
    if (!is.null(mat)) return(mean_from_matrix(mat, ndose = ndose))
  }
  
  if (!is.null(mean_candidates)) {
    vals <- lapply(res.list, get_field, candidates = mean_candidates, required = FALSE)
    keep <- !vapply(vals, is.null, logical(1))
    vals <- vals[keep]
    
    if (length(vals) > 0L) {
      ntrial <- vapply(res.list[keep], function(x) x$ntrial, numeric(1))
      ntrial.total <- sum(ntrial)
      out <- rep(0, ndose)
      for (i in seq_along(vals)) {
        out <- out + as.numeric(vals[[i]]) * ntrial[i]
      }
      return(out / ntrial.total)
    }
  }
  
  if (required) {
    stop(
      "Cannot extract field from candidates: ",
      paste(c(by_trial_candidates, mean_candidates), collapse = ", ")
    )
  }
  
  rep(NA_real_, ndose)
}

extract_mean_metric <- function(res.list, field_candidates, ndose = NULL, required = TRUE) {
  vals <- lapply(res.list, get_field, candidates = field_candidates, required = FALSE)
  keep <- !vapply(vals, is.null, logical(1))
  vals <- vals[keep]
  
  if (length(vals) == 0L) {
    if (required) stop("No metric field found: ", paste(field_candidates, collapse = ", "))
    if (is.null(ndose)) ndose <- get_field(res.list[[1]], c("ndose"))
    return(rep(NA_real_, ndose))
  }
  
  ntrial <- vapply(res.list[keep], function(x) x$ntrial, numeric(1))
  ntrial.total <- sum(ntrial)
  
  if (is.null(ndose)) ndose <- length(vals[[1]])
  
  out <- rep(0, ndose)
  for (i in seq_along(vals)) {
    out <- out + as.numeric(vals[[i]]) * ntrial[i]
  }
  out / ntrial.total
}

extract_est_pj <- function(res.list, model, ndose) {
  if (model == "CRM") {
    return(extract_mean_vector(
      res.list,
      by_trial_candidates = c(
        "crm_pj_by_trial",
        "p_hat_by_trial",
        "phat_by_trial",
        "pj_by_trial",
        "pj_iso_by_trial"
      ),
      mean_candidates = c(
        "crm_pj_mean",
        "p_hat_mean",
        "phat_mean",
        "pj_mean",
        "pj_iso_mean"
      ),
      ndose = ndose,
      required = TRUE
    ))
  }
  
  extract_mean_vector(
    res.list,
    by_trial_candidates = c("pj_iso_by_trial"),
    mean_candidates = c("pj_iso_mean"),
    ndose = ndose,
    required = TRUE
  )
}

safe_unlist_field <- function(res.list, candidates, required = TRUE) {
  vals <- lapply(res.list, get_field, candidates = candidates, required = required)
  unlist(vals, use.names = FALSE)
}

summarize_aide_files <- function(files.use,
                                 sc,
                                 model,
                                 method_tag,
                                 boin_method,
                                 boin_r_estimator,
                                 crm_r_model,
                                 alpha_true,
                                 r_carry,
                                 arrival_rate,
                                 cycle_max,
                                 continuous_enrollment) {
  res.list <- lapply(files.use, readRDS)
  
  ndose <- get_field(res.list[[1]], c("ndose"))
  ntrial.total <- sum(vapply(res.list, function(x) x$ntrial, numeric(1)))
  
  p.true <- as.numeric(get_field(res.list[[1]], c("p.true")))
  p.true_ipde <- as.numeric(get_field(res.list[[1]], c("p.true_ipde")))
  p.true_setB_expected <- if (as.character(sc) %in% rownames(p_true_setB)) {
    as.numeric(p_true_setB[as.character(sc), ])
  } else {
    rep(NA_real_, length(p.true))
  }
  if (!all(is.na(p.true_setB_expected))) {
    if (length(p.true) != length(p.true_setB_expected) ||
        any(abs(p.true - p.true_setB_expected) > 1e-10)) {
      warning(
        "Scenario ", sc,
        ": p.true stored in the RDS does not match Set B in this extractor. ",
        "The output keeps the RDS p.true values; check that the runner used Set B."
      )
    }
  }
  
  final_MTD <- safe_unlist_field(res.list, c("final_MTD"))
  
  est_pj <- extract_est_pj(res.list, model = model, ndose = ndose)
  
  selection_pct <- calc_select_rate_pct(
    sel = final_MTD,
    ndose = ndose,
    denom = length(final_MTD)
  )
  
  early_stop_pct <- 100 * sum(final_MTD == 99L, na.rm = TRUE) / length(final_MTD)
  
  n_by_dose <- extract_mean_metric(res.list, c("n_by_dose"), ndose = ndose)
  unique_n_by_dose <- extract_mean_metric(res.list, c("unique_n_by_dose"), ndose = ndose)
  nipde_by_dose <- extract_mean_metric(res.list, c("nipde_by_dose"), ndose = ndose)
  
  total_admin <- safe_unlist_field(res.list, c("total_admin"))
  total_unique <- safe_unlist_field(res.list, c("total_unique"))
  duration <- safe_unlist_field(res.list, c("duration"))
  
  total_admin_mean <- mean(total_admin, na.rm = TRUE)
  total_unique_mean <- mean(total_unique, na.rm = TRUE)
  duration_mean <- mean(duration, na.rm = TRUE)
  
  ## Available mostly for BOIN; for CRM fixed/random/level this will be used if stored.
  r_hat <- extract_mean_vector(
    res.list,
    by_trial_candidates = c("r_hat_by_trial"),
    mean_candidates = c("r_hat_mean"),
    ndose = ndose,
    required = FALSE
  )
  
  r_cap <- extract_mean_vector(
    res.list,
    by_trial_candidates = c("r_cap_by_trial"),
    mean_candidates = c("r_cap_mean"),
    ndose = ndose,
    required = FALSE
  )
  
  r_use <- extract_mean_vector(
    res.list,
    by_trial_candidates = c("r_use_by_trial"),
    mean_candidates = c("r_use_mean"),
    ndose = ndose,
    required = FALSE
  )
  
  dose_summary <- data.frame(
    Scenario_Set = scenario_set_name,
    Scenario_Name = paste0(scenario_set_name, "_SC", sc),
    Scenario = sc,
    Model = model,
    Method = method_tag,
    BOIN_Method = ifelse(model == "BOIN", boin_method, NA_character_),
    BOIN_r_estimator = ifelse(model == "BOIN", boin_r_estimator, NA_character_),
    CRM_r_model = ifelse(model == "CRM", crm_r_model, NA_character_),
    Alpha_true = alpha_true,
    r_carry = r_carry,
    Accrual = arrival_rate,
    T_assess = T_assess,
    Cycle_Max = cycle_max,
    Continuous_Enrollment = as.integer(isTRUE(continuous_enrollment)),
    Dose = seq_len(ndose),
    
    True_DLT_rate = p.true,
    True_IPDE_DLT_rate = p.true_ipde,
    Estimated_pj = est_pj,
    
    r_hat = r_hat,
    r_cap = r_cap,
    r_use = r_use,
    
    MTD_Selection_pct = selection_pct,
    Pts_Treated = n_by_dose,
    Unique_Pts_by_Dose = unique_n_by_dose,
    IPDE_Doses = nipde_by_dose,
    
    Total_Administrations = total_admin_mean,
    Total_Unique_Patients = total_unique_mean,
    Early_Stopping_pct = early_stop_pct,
    Duration = duration_mean,
    
    n_valid = length(final_MTD),
    ntrial_from_files = ntrial.total,
    CRM_Skeleton = paste(crm_skeleton, collapse = ","),
    CRM_Theta_Prior_Mean = crm_theta_prior_mean,
    CRM_Theta_Prior_SD = crm_theta_prior_sd,
    CRM_r_Prior_a = crm_a_r,
    CRM_r_Prior_b = crm_b_r,
    stringsAsFactors = FALSE
  )
  
  if (model == "CRM" && crm_r_model == "alpha_crm") {
    dose_summary$CRM_Dose_Values <- paste(crm_dose_values_alpha, collapse = ",")
    dose_summary$CRM_Theta_Prior_Mean <- crm_theta_prior_mean
    dose_summary$CRM_Theta_Prior_SD <- crm_theta_prior_sd
    dose_summary$CRM_Alpha_Grid_Length <- crm_alpha_grid_length
  } else {
    dose_summary$CRM_Dose_Values <- NA_character_
    ## Keep theta prior fields from the shared metadata columns above.
    dose_summary$CRM_Alpha_Grid_Length <- NA_integer_
  }
  
  if (model == "CRM" && crm_r_model == "cumu_crm") {
    dose_summary$CRM_Dose_Scores_Raw <- paste(crm_dose_scores_raw, collapse = ",")
    dose_summary$CRM_Dose_Scores <- paste(round(crm_dose_scores_cumu, 8), collapse = ",")
    dose_summary$CRM_Cumu_Beta0_Mean <- crm_cumu_beta0_mean
    dose_summary$CRM_Cumu_Beta0_Prec <- crm_cumu_beta0_prec
    dose_summary$CRM_Cumu_Beta0_DF <- crm_cumu_beta0_df
    dose_summary$CRM_Cumu_Beta1_Shape <- crm_cumu_beta1_shape
    dose_summary$CRM_Cumu_Beta1_Rate <- crm_cumu_beta1_rate
    dose_summary$CRM_Cumu_Beta2_Rate <- crm_cumu_beta2_rate
    dose_summary$CRM_Cumu_Include_Current <- as.integer(isTRUE(crm_cumu_include_current))
  } else {
    dose_summary$CRM_Dose_Scores_Raw <- NA_character_
    dose_summary$CRM_Dose_Scores <- NA_character_
    dose_summary$CRM_Cumu_Beta0_Mean <- NA_real_
    dose_summary$CRM_Cumu_Beta0_Prec <- NA_real_
    dose_summary$CRM_Cumu_Beta0_DF <- NA_real_
    dose_summary$CRM_Cumu_Beta1_Shape <- NA_real_
    dose_summary$CRM_Cumu_Beta1_Rate <- NA_real_
    dose_summary$CRM_Cumu_Beta2_Rate <- NA_real_
    dose_summary$CRM_Cumu_Include_Current <- NA_integer_
  }
  
  dose_cols <- paste0("D", seq_len(ndose))
  
  metrics <- c(
    "True DLT rate",
    "True IPDE DLT rate",
    ifelse(model == "CRM", "Estimated CRM pj", "Estimated pj_iso"),
    "r_hat",
    "r_cap",
    "r_use",
    "MTD Selection %",
    "# Pts Treated",
    "# Unique Pts by Dose",
    "# IPDE Doses",
    "# Total Administrations",
    "# Total Unique Patients",
    "% Early Stopping",
    "Duration"
  )
  
  table_summary <- data.frame(
    Scenario_Set = scenario_set_name,
    Scenario_Name = paste0(scenario_set_name, "_SC", sc),
    Scenario = sc,
    Model = model,
    Method = method_tag,
    BOIN_Method = ifelse(model == "BOIN", boin_method, NA_character_),
    BOIN_r_estimator = ifelse(model == "BOIN", boin_r_estimator, NA_character_),
    CRM_r_model = ifelse(model == "CRM", crm_r_model, NA_character_),
    Alpha_true = alpha_true,
    r_carry = r_carry,
    Accrual = arrival_rate,
    T_assess = T_assess,
    Cycle_Max = cycle_max,
    Continuous_Enrollment = as.integer(isTRUE(continuous_enrollment)),
    Metric = metrics,
    stringsAsFactors = FALSE
  )
  
  for (dc in dose_cols) table_summary[[dc]] <- NA_real_
  table_summary$Total <- NA_real_
  table_summary$Duration <- NA_real_
  table_summary$n_valid <- length(final_MTD)
  table_summary$ntrial_from_files <- ntrial.total
  
  table_summary[1, dose_cols] <- p.true
  table_summary[2, dose_cols] <- p.true_ipde
  table_summary[3, dose_cols] <- est_pj
  table_summary[4, dose_cols] <- r_hat
  table_summary[5, dose_cols] <- r_cap
  table_summary[6, dose_cols] <- r_use
  table_summary[7, dose_cols] <- selection_pct
  table_summary[8, dose_cols] <- n_by_dose
  table_summary[9, dose_cols] <- unique_n_by_dose
  table_summary[10, dose_cols] <- nipde_by_dose
  
  table_summary[11, "Total"] <- total_admin_mean
  table_summary[12, "Total"] <- total_unique_mean
  table_summary[13, "Total"] <- early_stop_pct
  table_summary[14, "Duration"] <- duration_mean
  
  list(
    dose_summary = dose_summary,
    table_summary = table_summary,
    n_valid = length(final_MTD),
    ntrial_from_files = ntrial.total
  )
}

## -------------------------------
## Main extraction
## -------------------------------

all.dose.summary <- list()
all.table.summary <- list()
missing.log <- list()

miss.idx <- 1L

for (model in model_list) {
  if (model == "BOIN") {
    method_loop <- boin_method_list
    boin_r_estimator_loop <- boin_r_estimator_list
    crm_loop <- NA_character_
  } else {
    method_loop <- NA_character_
    boin_r_estimator_loop <- NA_character_
    crm_loop <- crm_r_model_list
  }
  
  for (sc in scenario_id_list) {
    for (alpha_true in alpha_true_list) {
      for (arrival_rate in arrival_rate_list) {
        for (cycle_max in cycle_max_list) {
          for (method0 in method_loop) {
            for (boin_r_estimator in boin_r_estimator_loop) {
              for (crm_r_model in crm_loop) {
                
                if (model == "BOIN") {
                  r_loop <- if (identical(boin_r_estimator, "r_fixed")) {
                    boin_r_carry_list_fixed
                  } else {
                    boin_r_carry_list_r_mle
                  }
                } else {
                  r_loop <- get_crm_r_loop(crm_r_model)
                }
                
                for (r_carry in r_loop) {
                  method_tag <- make_method_tag(
                    model = model,
                    boin_method = if (model == "BOIN") method0 else NULL,
                    boin_r_estimator = if (model == "BOIN") boin_r_estimator else NULL,
                    crm_r_model = if (model == "CRM") crm_r_model else NULL
                  )
                  
                  foldername <- make_foldername(
                    model = model,
                    method_tag = method_tag,
                    cycle_max = cycle_max,
                    arrival_rate = arrival_rate
                  )
                  
                  folderpath <- file.path(results_root, foldername)
                  
                  files <- file.path(
                    folderpath,
                    vapply(
                      jobs.expected,
                      function(j) {
                        make_filename(
                          sc = sc,
                          model = model,
                          method_tag = method_tag,
                          alpha_true = alpha_true,
                          r_carry = r_carry,
                          arrival_rate = arrival_rate,
                          cycle_max = cycle_max,
                          continuous_enrollment = continuous_enrollment,
                          job = j
                        )
                      },
                      character(1)
                    )
                  )
                  
                  exists.vec <- file.exists(files)
                  files.use <- files[exists.vec]
                  missing.jobs <- jobs.expected[!exists.vec]
                  
                  ## Fallback: discover matching combined files in the folder.
                  ## This is useful if a run was started at a non-1 job index,
                  ## or if jobs.expected does not exactly match the submitted array.
                  if (length(files.use) == 0L && dir.exists(folderpath)) {
                    strict_key <- make_group_key(
                      sc = sc,
                      model = model,
                      method_tag = method_tag,
                      alpha_true = alpha_true,
                      r_carry = r_carry,
                      arrival_rate = arrival_rate,
                      cycle_max = cycle_max,
                      continuous_enrollment = continuous_enrollment
                    )
                    pattern <- paste0("^", strict_key, "-job-[0-9]+-combined\\.rds$")
                    files.use <- list.files(
                      folderpath,
                      pattern = pattern,
                      full.names = TRUE
                    )
                    if (length(files.use) > 0L) {
                      found_jobs <- as.integer(sub(
                        paste0(".*-job-([0-9]+)-combined\\.rds$"),
                        "\\1",
                        basename(files.use)
                      ))
                      missing.jobs <- setdiff(jobs.expected, found_jobs)
                    }
                  }
                  
                  cat("\n====================================\n")
                  cat("Scenario set:", scenario_set_name, "\n")
                  cat("Scenario:", sc, "\n")
                  cat("Model:", model, "\n")
                  cat("Method tag:", method_tag, "\n")
                  cat("BOIN method:", ifelse(model == "BOIN", method0, NA), "\n")
                  cat("BOIN r estimator:", ifelse(model == "BOIN", boin_r_estimator, NA), "\n")
                  cat("CRM r model:", ifelse(model == "CRM", crm_r_model, NA), "\n")
                  cat("alpha_true:", alpha_true, "\n")
                  cat("r_carry:", r_carry, "\n")
                  cat("Accrual:", arrival_rate, "\n")
                  cat("Cycle max:", cycle_max, "\n")
                  cat("Continuous:", continuous_enrollment, "\n")
                  cat("Folder:", folderpath, "\n")
                  cat("Example file:", files[1], "\n")
                  cat("Found", length(files.use), "files; missing", length(missing.jobs), "jobs.\n")
                  cat("====================================\n")
                  
                  missing.log[[miss.idx]] <- data.frame(
                    Scenario_Set = scenario_set_name,
                    Scenario_Name = paste0(scenario_set_name, "_SC", sc),
                    Scenario = sc,
                    Model = model,
                    Method = method_tag,
                    BOIN_Method = ifelse(model == "BOIN", method0, NA_character_),
                    BOIN_r_estimator = ifelse(model == "BOIN", boin_r_estimator, NA_character_),
                    CRM_r_model = ifelse(model == "CRM", crm_r_model, NA_character_),
                    Alpha_true = alpha_true,
                    r_carry = r_carry,
                    Accrual = arrival_rate,
                    T_assess = T_assess,
                    Cycle_Max = cycle_max,
                    Continuous_Enrollment = as.integer(isTRUE(continuous_enrollment)),
                    Folder = folderpath,
                    Example_file = files[1],
                    n_found = length(files.use),
                    n_missing = length(missing.jobs),
                    missing_jobs = paste(missing.jobs, collapse = ","),
                    stringsAsFactors = FALSE
                  )
                  miss.idx <- miss.idx + 1L
                  
                  if (length(files.use) == 0L) {
                    warning("No files found for folder: ", folderpath)
                    next
                  }
                  
                  t.read <- Sys.time()
                  
                  one <- summarize_aide_files(
                    files.use = files.use,
                    sc = sc,
                    model = model,
                    method_tag = method_tag,
                    boin_method = method0,
                    boin_r_estimator = boin_r_estimator,
                    crm_r_model = crm_r_model,
                    alpha_true = alpha_true,
                    r_carry = r_carry,
                    arrival_rate = arrival_rate,
                    cycle_max = cycle_max,
                    continuous_enrollment = continuous_enrollment
                  )
                  
                  cat(
                    "read/summarize time:",
                    round(difftime(Sys.time(), t.read, units = "secs"), 2),
                    "sec\n"
                  )
                  
                  if (one$n_valid != ntrial.expected) {
                    cat(
                      "Warning: expected", ntrial.expected,
                      "trials but found", one$n_valid, "\n"
                    )
                  }
                  
                  key <- make_group_key(
                    sc = sc,
                    model = model,
                    method_tag = method_tag,
                    alpha_true = alpha_true,
                    r_carry = r_carry,
                    arrival_rate = arrival_rate,
                    cycle_max = cycle_max,
                    continuous_enrollment = continuous_enrollment
                  )
                  
                  all.dose.summary[[key]] <- one$dose_summary
                  all.table.summary[[key]] <- one$table_summary
                }
              }
            }
          }
        }
      }
    }
  }
}

if (length(all.dose.summary) == 0L) {
  stop("No files were found. Check results_root, folder names, filenames, and working directory.")
}

all.dose.summary.df <- do.call(rbind, all.dose.summary)
all.table.summary.df <- do.call(rbind, all.table.summary)
missing.log.df <- do.call(rbind, missing.log)

round_numeric_df <- function(dat, digits = 4) {
  num_cols <- vapply(dat, is.numeric, logical(1))
  dat[num_cols] <- lapply(dat[num_cols], round, digits)
  dat
}

all.dose.summary.out <- round_numeric_df(all.dose.summary.df, digits = 4)
all.table.summary.out <- round_numeric_df(all.table.summary.df, digits = 4)

out.tag <- paste0(
  paste0("All_AIDE_OC_", scenario_set_name),
  "_models", paste(model_list, collapse = "_"),
  if ("BOIN" %in% model_list) paste0("_boin", paste(boin_method_list, collapse = "_")) else "",
  if ("BOIN" %in% model_list) paste0("_rest", paste(boin_r_estimator_list, collapse = "_")) else "",
  if ("CRM" %in% model_list) paste0("_crm", paste(crm_r_model_list, collapse = "_")) else "",
  "_target", fmt_num(target),
  "_w", fmt_short(T_assess),
  "_c", fmt_short(C),
  "_cyc", paste(fmt_short(cycle_max_list), collapse = "_"),
  "_rate", paste(fmt_short(arrival_rate_list), collapse = "_"),
  "_Nmax", fmt_short(Nmax_eff),
  "_dosecap", fmt_short(dose_cap),
  "_cont", as.integer(isTRUE(continuous_enrollment)),
  "_jobs", min(jobs.expected), "to", max(jobs.expected)
)

out.dose.csv <- file.path(out.dir, paste0(out.tag, "_dose_summary.csv"))
out.table.csv <- file.path(out.dir, paste0(out.tag, "_table_summary.csv"))
out.missing.csv <- file.path(out.dir, paste0(out.tag, "_missing_jobs.csv"))

out.dose.rds <- file.path(out.dir, paste0(out.tag, "_dose_summary.rds"))
out.table.rds <- file.path(out.dir, paste0(out.tag, "_table_summary.rds"))
out.missing.rds <- file.path(out.dir, paste0(out.tag, "_missing_jobs.rds"))

write.csv(all.dose.summary.out, out.dose.csv, row.names = FALSE)
write.csv(all.table.summary.out, out.table.csv, row.names = FALSE)
write.csv(missing.log.df, out.missing.csv, row.names = FALSE)

saveRDS(all.dose.summary.df, out.dose.rds)
saveRDS(all.table.summary.df, out.table.rds)
saveRDS(missing.log.df, out.missing.rds)

cat("\nSaved dose-level summary:", out.dose.csv, "\n")
cat("Saved table-style summary:", out.table.csv, "\n")
cat("Saved missing-job log:", out.missing.csv, "\n")
