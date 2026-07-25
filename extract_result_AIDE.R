## ============================================================
## Extract AIDE OC from cluster jobs
## Matched to run_oc_AIDE.R
##
## Newest runner naming rules:
##
## Folder:
##   oc_results_cluster_AIDE_{scenario_set_name}/
##   {scenario_set_name}-model-{model}-opt-{method_tag}-w-{T_assess}-c-{C}-
##   cyc-{cycle_max}-rate-{arrival_rate}-Nmax-{Nmax_eff}-ipde-{1/2}-
##   newfirst-{1/2}-neval-{n_eval_escalate}-tried-{0/1}/
##
## Raw RDS file written by run_oc_AIDE.R:
##   {scenario_set_name}_SC{sc}-{model}-{method_tag}-a{alpha_true}-r{r_carry}-
##   cyc{cycle_max}-tried{0/1}-j{job}-b{block}-s{seed}-n{ntrial}.rds
##
## Supports:
##   BOIN: boin / approx1 / approx2 with r_fixed, r_mle, or r_adaptive
##   CRM : fixed / random / level / alpha_crm / cumu_crm.  Current
##         run_oc_AIDE.R writes fixed CRM as crm_r_fixed; older result folders
##         may use the r_fixed alias and are discovered automatically.
##   CFO : empirical / pride
## ============================================================

rm(list = ls())

## -------------------------------
## Settings matched to newest runner
## -------------------------------

## setwd("/rsrch8/home/biostatistics/syang10/AIDE")

## Current non-random 5-dose runner branch: scenarios 26:37.
scenario_file <- file.path("scenario_sets", "Set_5dose_adaptive_r_37.csv")
scenario_set_name <- tools::file_path_sans_ext(basename(scenario_file))
results_root <- paste0("oc_results_cluster_AIDE_", scenario_set_name)

scenario_meta <- utils::read.csv(scenario_file, check.names = FALSE, stringsAsFactors = FALSE)
dose_col_names <- grep("^Dose[0-9]+$", names(scenario_meta), value = TRUE)
p_true_scenarios <- as.matrix(scenario_meta[dose_col_names])
rownames(p_true_scenarios) <- as.character(scenario_meta$Scenario)
scenario_id_list <- 1:6
ndose_expected <- ncol(p_true_scenarios)

jobs.expected <- 1:2000

## If each LSF job used ntrial.total = 1, expected total is 2000.
## Change this if you pass a different trials-per-job argument to the runner.
ntrial_per_job_expected <- 1L
seed_base_expected <- 1L
block_id_expected <- 1L
ntrial.expected <- length(jobs.expected) * ntrial_per_job_expected

target <- 0.30

T_assess <- 28
C <- 3L
cycle_max_list <- c(2L)
Nmax_eff_list <- c(30L)
## 1 = recycle from any lower dose; 2 = adjacent lower dose only.
ipde_design_list <- c(1L)
## 1 = prioritize new patients; 2 = prioritize eligible IPDE patients.
new_pat_first_list <- c(1L)
## Escalation is allowed once this many patients are evaluated at the dose.
n_eval_escalate_list <- c(3L)
## Runner settings: d.cap is the per-dose enrollment cap; dose_cap limits doses.
d_cap <- 100L
dose_cap <- 3L
continuous_enrollment <- TRUE
restrict_to_tried_list <- c(TRUE)

alpha_true_list <- c(0, 0.3, 0.6, 0.9)
arrival_rate_list <- c(1 / 28, 1 / 56)

## BOIN settings matched to runner
boin_method_list <- c("approx1", "approx2")
boin_r_estimator_list <- c("fixed")
boin_r_carry_list_r_mle <- c(0)
boin_r_carry_list_fixed <- c(0)

## CRM settings matched to the current 5-dose random-scenario runner.
crm_r_model_list <- c("crm_r_fixed", "random")

## r values by CRM backend, matching the names produced by run_oc_AIDE.R.
## fixed CRM uses r_carry; alpha_crm/cumu_crm do not use fixed r, so the
## runner convention is usually r0 in the filename.
crm_r_carry_values <- list(
  crm_r_fixed = c(0),
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

## 5-dose skeleton and dose scores from run_oc_AIDE.R.
crm_skeleton <- c(0.15, 0.20, 0.30, 0.35, 0.45)
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
crm_dose_values_alpha <- c(15, 20, 30, 35, 45)
dose_alpha_mg <- crm_dose_values_alpha
crm_alpha_grid_length <- 61L

## IPCRM / cumulative CRM current dose scores and priors
## beta0 ~ t(fixed_intercept, precision = beta0_prec, df = beta0_df)
## beta1 ~ Gamma(beta1_shape, beta1_rate)
## beta2 ~ Exp(1), but beta2 drops out at baseline when cumu.d = 0
crm_dose_scores_raw <- c(15, 20, 30, 35, 45)
dose_ipcrm <- crm_dose_scores_raw / (2 * stats::sd(crm_dose_scores_raw))
crm_dose_scores_cumu <- dose_ipcrm
fixed_intercept <- -2.8
beta0_prec <- 2
beta0_df <- 1
beta1_shape <- 2.5
beta1_rate <- 1.6
beta2_rate <- 1
crm_cumu_beta0_mean <- fixed_intercept
crm_cumu_beta0_prec <- beta0_prec
crm_cumu_beta0_df <- beta0_df
crm_cumu_beta1_shape <- beta1_shape
crm_cumu_beta1_rate <- beta1_rate
crm_cumu_beta2_rate <- beta2_rate
crm_cumu_include_current <- FALSE

## CFO / PRIDE settings matched to methods_prior.R and run_oc_AIDE.R.
cfo_method_list <- c("empirical", "pride")
cfo_skeleton <- c(0.002, 0.008, 0.012, 0.04, 0.08, 0.10, 0.20, 0.35)
cfo_model_file <- "PRIDE.bug"
cfo_sigma2_beta <- 30
cfo_eta <- 1
cfo_pk_method <- "approx"
cfo_n_mc_w <- 200
cfo_m_use <- 1000
cfo_use_monotone_pair <- FALSE
cfo_n_chains <- 3
cfo_n_adapt <- 1000
cfo_n_burnin <- 2000
cfo_n_iter <- 5000
cfo_thin <- 2

model_list <- c("CRM")

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
                            crm_r_model = NULL,
                            cfo_method = NULL) {
  if (model == "BOIN") {
    return(paste0(boin_method, "-", boin_r_estimator))
  }
  if (model == "CRM") {
    if (identical(crm_r_model, "fixed")) return("crm_r_fixed")
    return(crm_r_model)
  }
  paste0("cfo_", cfo_method)
}

make_foldername <- function(model,
                            method_tag,
                            cycle_max,
                            arrival_rate,
                            Nmax_eff,
                            ipde_design,
                            new_pat_first,
                            n_eval_escalate,
                            restrict_to_tried,
                            scenario_set = scenario_set_name) {
  paste0(
    scenario_set,
    "-model-", model,
    "-opt-", method_tag,
    "-w-", fmt_short(T_assess),
    "-c-", fmt_short(C),
    "-cyc-", fmt_short(cycle_max),
    "-rate-", fmt_short(arrival_rate),
    "-Nmax-", fmt_short(Nmax_eff),
    "-ipde-", ipde_design,
    "-newfirst-", new_pat_first,
    "-neval-", n_eval_escalate,
    "-tried-", as.integer(isTRUE(restrict_to_tried))
  )
}

resolve_folderpath <- function(foldername) {
  nested <- file.path(results_root, foldername)
  if (dir.exists(nested) || !dir.exists(foldername)) {
    return(nested)
  }
  foldername
}

method_tag_candidates <- function(model, method_tag) {
  if (identical(model, "CRM") && identical(method_tag, "crm_r_fixed")) {
    ## crm_r_fixed is the tag emitted by the current runner.  Keep r_fixed as
    ## a legacy read-only alias for previously submitted result folders.
    return(c("crm_r_fixed", "r_fixed"))
  }
  method_tag
}

resolve_result_folder <- function(model,
                                  method_tag,
                                  cycle_max,
                                  arrival_rate,
                                  Nmax_eff,
                                  ipde_design,
                                  new_pat_first,
                                  n_eval_escalate,
                                  restrict_to_tried) {
  tags <- method_tag_candidates(model, method_tag)
  attempts <- lapply(tags, function(tag) {
    foldername <- make_foldername(
      model = model,
      method_tag = tag,
      cycle_max = cycle_max,
      arrival_rate = arrival_rate,
      Nmax_eff = Nmax_eff,
      ipde_design = ipde_design,
      new_pat_first = new_pat_first,
      n_eval_escalate = n_eval_escalate,
      restrict_to_tried = restrict_to_tried
    )
    list(
      method_tag = tag,
      foldername = foldername,
      folderpath = resolve_folderpath(foldername)
    )
  })

  found <- vapply(attempts, function(x) dir.exists(x$folderpath), logical(1))
  attempts[[if (any(found)) which(found)[1L] else 1L]]
}

make_group_key <- function(sc,
                           model,
                           method_tag,
                           alpha_true,
                           r_carry,
                           arrival_rate,
                           cycle_max,
                           Nmax_eff,
                           ipde_design,
                           continuous_enrollment,
                           new_pat_first,
                           n_eval_escalate,
                           restrict_to_tried,
                           scenario_set = scenario_set_name) {
  paste(
    paste0(scenario_set, "_SC", sc),
    model,
    method_tag,
    paste0("a", fmt_short(alpha_true)),
    paste0("r", fmt_short(r_carry)),
    paste0("rate", fmt_short(arrival_rate)),
    paste0("cyc", fmt_short(cycle_max)),
    paste0("Nmax", fmt_short(Nmax_eff)),
    paste0("ipde", ipde_design),
    paste0("cont", as.integer(isTRUE(continuous_enrollment))),
    paste0("newfirst", new_pat_first),
    paste0("neval", n_eval_escalate),
    paste0("tried", as.integer(isTRUE(restrict_to_tried))),
    sep = "_"
  )
}

## Each raw task RDS has the fields used by summarize_aide_files(), so exact
## raw files written by run_oc_AIDE.R can be summarized directly.
make_raw_rds_filenames <- function(sc,
                                   model,
                                   method_tag,
                                   alpha_true,
                                   r_carry,
                                   cycle_max,
                                   restrict_to_tried,
                                   jobs = jobs.expected,
                                   block_id = block_id_expected,
                                   seed_base = seed_base_expected,
                                   ntrial_per_job = ntrial_per_job_expected,
                                   scenario_set = scenario_set_name) {
  jobs <- as.integer(jobs)
  seed_block <- as.integer(seed_base) + (jobs - 1L) * as.integer(ntrial_per_job)

  paste0(
    scenario_set, "_SC", sc,
    "-", model,
    "-", method_tag,
    "-a", fmt_short(alpha_true),
    "-r", fmt_short(r_carry),
    "-cyc", fmt_short(cycle_max),
    "-tried", as.integer(isTRUE(restrict_to_tried)),
    "-j", jobs,
    "-b", block_id,
    "-s", seed_block,
    "-n", as.integer(ntrial_per_job),
    ".rds"
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
                                 cfo_method,
                                 alpha_true,
                                 r_carry,
                                 arrival_rate,
                                 Nmax_eff,
                                 cycle_max,
                                 continuous_enrollment,
                                 ipde_design,
                                 new_pat_first,
                                 n_eval_escalate,
                                 restrict_to_tried) {
  res.list <- lapply(files.use, readRDS)

  ndose <- get_field(res.list[[1]], c("ndose"))
  ntrial.total <- sum(vapply(res.list, function(x) x$ntrial, numeric(1)))

  p.true <- as.numeric(get_field(res.list[[1]], c("p.true")))
  p.true_ipde <- as.numeric(get_field(res.list[[1]], c("p.true_ipde")))
  p.true_expected <- if (as.character(sc) %in% rownames(p_true_scenarios)) {
    as.numeric(p_true_scenarios[as.character(sc), ])
  } else {
    rep(NA_real_, length(p.true))
  }
  if (!all(is.na(p.true_expected))) {
    if (length(p.true) != length(p.true_expected) ||
        any(abs(p.true - p.true_expected) > 1e-10)) {
      warning(
        "Scenario ", sc,
        ": p.true stored in the RDS does not match this extractor's scenario set. ",
        "The output keeps the RDS p.true values; check that the runner and extractor match."
      )
    }
  }

  source_scenario <- get_field(res.list[[1]], c("source_scenario"), required = FALSE)
  if (is.null(source_scenario) && sc %in% scenario_meta$Scenario) {
    source_scenario <- scenario_meta$Source_Scenario[match(sc, scenario_meta$Scenario)]
  }
  if (is.null(source_scenario)) source_scenario <- NA_integer_

  true_mtd <- get_field(res.list[[1]], c("true_mtd"), required = FALSE)
  if (is.null(true_mtd) && sc %in% scenario_meta$Scenario) {
    true_mtd <- scenario_meta$True_MTD[match(sc, scenario_meta$Scenario)]
  }
  if (is.null(true_mtd)) true_mtd <- NA_integer_

  scenario_attempt <- get_field(res.list[[1]], c("scenario_attempt"), required = FALSE)
  if (is.null(scenario_attempt) && sc %in% scenario_meta$Scenario) {
    scenario_attempt <- scenario_meta$Attempt[match(sc, scenario_meta$Scenario)]
  }
  if (is.null(scenario_attempt)) scenario_attempt <- NA_integer_

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
    Source_Scenario = source_scenario,
    True_MTD = true_mtd,
    Scenario_Attempt = scenario_attempt,
    Model = model,
    Method = method_tag,
    BOIN_Method = ifelse(model == "BOIN", boin_method, NA_character_),
    BOIN_r_estimator = ifelse(model == "BOIN", boin_r_estimator, NA_character_),
    CRM_r_model = ifelse(model == "CRM", crm_r_model, NA_character_),
    CFO_Method = ifelse(model == "CFO", cfo_method, NA_character_),
    Alpha_true = alpha_true,
    r_carry = r_carry,
    Accrual = arrival_rate,
    T_assess = T_assess,
    Nmax_eff = Nmax_eff,
    Cycle_Max = cycle_max,
    IPDE_Design = ipde_design,
    New_Patient_First = new_pat_first,
    N_Eval_Escalate = n_eval_escalate,
    D_Cap = d_cap,
    Dose_Cap = dose_cap,
    Continuous_Enrollment = as.integer(isTRUE(continuous_enrollment)),
    Restrict_To_Tried = as.integer(isTRUE(restrict_to_tried)),
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
    CFO_Skeleton = ifelse(model == "CFO", paste(cfo_skeleton, collapse = ","), NA_character_),
    CFO_Model_File = ifelse(model == "CFO" && cfo_method == "pride", cfo_model_file, NA_character_),
    CFO_Sigma2_Beta = ifelse(model == "CFO", cfo_sigma2_beta, NA_real_),
    CFO_Eta = ifelse(model == "CFO", cfo_eta, NA_real_),
    CFO_PK_Method = ifelse(model == "CFO", cfo_pk_method, NA_character_),
    CFO_N_MC_W = ifelse(model == "CFO", cfo_n_mc_w, NA_integer_),
    CFO_M_Use = ifelse(model == "CFO", cfo_m_use, NA_integer_),
    CFO_Use_Monotone_Pair = ifelse(model == "CFO", as.integer(isTRUE(cfo_use_monotone_pair)), NA_integer_),
    CFO_N_Chains = ifelse(model == "CFO", cfo_n_chains, NA_integer_),
    CFO_N_Adapt = ifelse(model == "CFO", cfo_n_adapt, NA_integer_),
    CFO_N_Burnin = ifelse(model == "CFO", cfo_n_burnin, NA_integer_),
    CFO_N_Iter = ifelse(model == "CFO", cfo_n_iter, NA_integer_),
    CFO_Thin = ifelse(model == "CFO", cfo_thin, NA_integer_),
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
    ifelse(model == "CRM", "Estimated CRM pj", ifelse(model == "CFO", "Estimated CFO pj_iso", "Estimated pj_iso")),
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
    Source_Scenario = source_scenario,
    True_MTD = true_mtd,
    Scenario_Attempt = scenario_attempt,
    Model = model,
    Method = method_tag,
    BOIN_Method = ifelse(model == "BOIN", boin_method, NA_character_),
    BOIN_r_estimator = ifelse(model == "BOIN", boin_r_estimator, NA_character_),
    CRM_r_model = ifelse(model == "CRM", crm_r_model, NA_character_),
    CFO_Method = ifelse(model == "CFO", cfo_method, NA_character_),
    Alpha_true = alpha_true,
    r_carry = r_carry,
    Accrual = arrival_rate,
    T_assess = T_assess,
    Nmax_eff = Nmax_eff,
    Cycle_Max = cycle_max,
    IPDE_Design = ipde_design,
    New_Patient_First = new_pat_first,
    N_Eval_Escalate = n_eval_escalate,
    D_Cap = d_cap,
    Dose_Cap = dose_cap,
    Continuous_Enrollment = as.integer(isTRUE(continuous_enrollment)),
    Restrict_To_Tried = as.integer(isTRUE(restrict_to_tried)),
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

setting_grid <- expand.grid(
  Nmax_eff = Nmax_eff_list,
  ipde_design = ipde_design_list,
  new_pat_first = new_pat_first_list,
  n_eval_escalate = n_eval_escalate_list,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

for (setting_i in seq_len(nrow(setting_grid))) {
  Nmax_eff <- setting_grid$Nmax_eff[setting_i]
  ipde_design <- setting_grid$ipde_design[setting_i]
  new_pat_first <- setting_grid$new_pat_first[setting_i]
  n_eval_escalate <- setting_grid$n_eval_escalate[setting_i]
  for (restrict_to_tried in restrict_to_tried_list) {
    for (model in model_list) {
      if (model == "BOIN") {
        method_loop <- boin_method_list
        boin_r_estimator_loop <- boin_r_estimator_list
        crm_loop <- NA_character_
        cfo_loop <- NA_character_
      } else if (model == "CRM") {
        method_loop <- NA_character_
        boin_r_estimator_loop <- NA_character_
        crm_loop <- crm_r_model_list
        cfo_loop <- NA_character_
      } else {
        method_loop <- NA_character_
        boin_r_estimator_loop <- NA_character_
        crm_loop <- NA_character_
        cfo_loop <- cfo_method_list
      }

  for (sc in scenario_id_list) {
    for (alpha_true in alpha_true_list) {
      for (arrival_rate in arrival_rate_list) {
        for (cycle_max in cycle_max_list) {
          for (method0 in method_loop) {
            for (boin_r_estimator in boin_r_estimator_loop) {
              for (crm_r_model in crm_loop) {
                for (cfo_method in cfo_loop) {

                  if (model == "BOIN") {
                    r_loop <- if (identical(boin_r_estimator, "r_fixed")) {
                      boin_r_carry_list_fixed
                    } else {
                      boin_r_carry_list_r_mle
                    }
                  } else if (model == "CRM") {
                    r_loop <- get_crm_r_loop(crm_r_model)
                  } else {
                    r_loop <- c(0)
                  }

                  for (r_carry in r_loop) {
                    method_tag <- make_method_tag(
                      model = model,
                      boin_method = if (model == "BOIN") method0 else NULL,
                      boin_r_estimator = if (model == "BOIN") boin_r_estimator else NULL,
                      crm_r_model = if (model == "CRM") crm_r_model else NULL,
                      cfo_method = if (model == "CFO") cfo_method else NULL
                    )

                    folder_info <- resolve_result_folder(
                      model = model,
                      method_tag = method_tag,
                      cycle_max = cycle_max,
                      arrival_rate = arrival_rate,
                      Nmax_eff = Nmax_eff,
                      ipde_design = ipde_design,
                      new_pat_first = new_pat_first,
                      n_eval_escalate = n_eval_escalate,
                      restrict_to_tried = restrict_to_tried
                    )
                    folderpath <- folder_info$folderpath
                    file_method_tag <- folder_info$method_tag

                    raw_names <- make_raw_rds_filenames(
                      sc = sc,
                      model = model,
                      method_tag = file_method_tag,
                      alpha_true = alpha_true,
                      r_carry = r_carry,
                      cycle_max = cycle_max,
                      restrict_to_tried = restrict_to_tried,
                      jobs = jobs.expected
                    )
                    raw_paths <- file.path(folderpath, raw_names)
                    existing_raw <- file.exists(raw_paths)
                    files.use <- raw_paths[existing_raw]
                    missing.jobs <- jobs.expected[!existing_raw]
                    example_file <- if (length(files.use) > 0L) {
                      files.use[1L]
                    } else {
                      raw_paths[1L]
                    }

                    cat("\n====================================\n")
                    cat("Scenario set:", scenario_set_name, "\n")
                    cat("Scenario:", sc, "\n")
                    cat("Model:", model, "\n")
                    cat("Method tag:", method_tag, "\n")
                    cat("BOIN method:", ifelse(model == "BOIN", method0, NA), "\n")
                    cat("BOIN r estimator:", ifelse(model == "BOIN", boin_r_estimator, NA), "\n")
                    cat("CRM r model:", ifelse(model == "CRM", crm_r_model, NA), "\n")
                    cat("CFO method:", ifelse(model == "CFO", cfo_method, NA), "\n")
                    cat("alpha_true:", alpha_true, "\n")
                    cat("r_carry:", r_carry, "\n")
                    cat("Accrual:", arrival_rate, "\n")
                    cat("Cycle max:", cycle_max, "\n")
                    cat("Nmax eff:", Nmax_eff, "\n")
                    cat("IPDE design:", ipde_design, "\n")
                    cat("New patient first:", new_pat_first, "\n")
                    cat("N evaluated for escalation:", n_eval_escalate, "\n")
                    cat("Continuous:", continuous_enrollment, "\n")
                    cat("Restrict to tried:", restrict_to_tried, "\n")
                    cat("Folder:", folderpath, "\n")
                    cat("Raw-RDS file:", example_file, "\n")
                    cat("Found", length(files.use), "files; missing", length(missing.jobs), "jobs.\n")
                    cat("====================================\n")

                    missing.log[[miss.idx]] <- data.frame(
                      Scenario_Set = scenario_set_name,
                      Scenario_Name = paste0(scenario_set_name, "_SC", sc),
                      Scenario = sc,
                      Source_Scenario = scenario_meta$Source_Scenario[match(sc, scenario_meta$Scenario)],
                      True_MTD = scenario_meta$True_MTD[match(sc, scenario_meta$Scenario)],
                      Scenario_Attempt = scenario_meta$Attempt[match(sc, scenario_meta$Scenario)],
                      Model = model,
                      Method = method_tag,
                      BOIN_Method = ifelse(model == "BOIN", method0, NA_character_),
                      BOIN_r_estimator = ifelse(model == "BOIN", boin_r_estimator, NA_character_),
                      CRM_r_model = ifelse(model == "CRM", crm_r_model, NA_character_),
                      CFO_Method = ifelse(model == "CFO", cfo_method, NA_character_),
                      Alpha_true = alpha_true,
                      r_carry = r_carry,
                      Accrual = arrival_rate,
                      T_assess = T_assess,
                      Nmax_eff = Nmax_eff,
                      Cycle_Max = cycle_max,
                      IPDE_Design = ipde_design,
                      New_Patient_First = new_pat_first,
                      N_Eval_Escalate = n_eval_escalate,
                      D_Cap = d_cap,
                      Dose_Cap = dose_cap,
                      Continuous_Enrollment = as.integer(isTRUE(continuous_enrollment)),
                      Restrict_To_Tried = as.integer(isTRUE(restrict_to_tried)),
                      Folder = folderpath,
                      Example_file = example_file,
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
                      cfo_method = cfo_method,
                      alpha_true = alpha_true,
                      r_carry = r_carry,
                      arrival_rate = arrival_rate,
                      Nmax_eff = Nmax_eff,
                      cycle_max = cycle_max,
                      continuous_enrollment = continuous_enrollment,
                      ipde_design = ipde_design,
                      new_pat_first = new_pat_first,
                      n_eval_escalate = n_eval_escalate,
                      restrict_to_tried = restrict_to_tried
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
                      Nmax_eff = Nmax_eff,
                      ipde_design = ipde_design,
                      continuous_enrollment = continuous_enrollment,
                      new_pat_first = new_pat_first,
                      n_eval_escalate = n_eval_escalate,
                      restrict_to_tried = restrict_to_tried
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
  if ("CFO" %in% model_list) paste0("_cfo", paste(cfo_method_list, collapse = "_")) else "",
  "_target", fmt_num(target),
  "_w", fmt_short(T_assess),
  "_c", fmt_short(C),
  "_cyc", paste(fmt_short(cycle_max_list), collapse = "_"),
  "_rate", paste(fmt_short(arrival_rate_list), collapse = "_"),
  "_Nmax", paste(fmt_short(Nmax_eff_list), collapse = "_"),
  "_ipde", paste(ipde_design_list, collapse = "_"),
  "_newfirst", paste(new_pat_first_list, collapse = "_"),
  "_neval", paste(n_eval_escalate_list, collapse = "_"),
  "_dcap", fmt_short(d_cap),
  "_dosecap", fmt_short(dose_cap),
  "_cont", as.integer(isTRUE(continuous_enrollment)),
  "_tried", paste(as.integer(restrict_to_tried_list), collapse = "_"),
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
