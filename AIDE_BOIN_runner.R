## ============================================================
## AIDE-BOIN / AIDE-CFO runner and OC table printer
## Source this after defining scenarios / global simulation settings,
## or edit the fallback settings below.
## ============================================================
source("AIDE_BOIN_helper.R")
source("AIDE_CRM_helper_final.R")
source("AIDE_modified.R")

## ------------------------------------------------------------
## Helper functions
## ------------------------------------------------------------
fmt_param <- function(x, digits = 4) {
  x <- round(as.numeric(x), digits)
  gsub("\\.", "p", format(x, scientific = FALSE, trim = TRUE))
}

with_sink <- function(file, expr) {
  con <- file(file, open = "wt")
  sink(con)
  on.exit({
    sink()
    close(con)
  }, add = TRUE)
  force(expr)
}

print_oc_table_aide <- function(
    x,
    scenario,
    design,
    method = "AIDE-BOIN",
    true_ipde_alpha = NA_real_,
    discount_alpha = NA_real_,
    allow_suspend = NA,
    accrual = NA_real_,
    digits = 2
) {
  K <- x$ndose

  if (is.null(K) || length(K) != 1L || is.na(K)) {
    K <- length(x$p.true)
  }

  K <- as.integer(K)
  dose_cols <- paste0("D", seq_len(K))

  metrics <- c(
    "True DLT rate",
    "True IPDE DLT rate",
    "Estimated pj_iso",
    "MTD Selection %",
    "# Pts Treated",
    "# Unique Pts by Dose",
    "# IPDE Doses",
    "# Total Administrations",
    "# Total Unique Patients",
    "% Early Stopping",
    "Duration"
  )

  tab <- data.frame(
    Metric = metrics,
    matrix(NA_real_, nrow = length(metrics), ncol = K)
  )

  names(tab)[2:(K + 1)] <- dose_cols
  tab$Total <- NA_real_
  tab$Duration <- NA_real_

  tab[1, dose_cols] <- x$p.true
  tab[2, dose_cols] <- x$p.true_ipde

  if (!is.null(x$pj_iso_mean)) {
    tab[3, dose_cols] <- x$pj_iso_mean
  }

  tab[4, dose_cols] <- x$selection_pct
  tab[5, dose_cols] <- x$n_by_dose
  tab[6, dose_cols] <- x$unique_n_by_dose
  tab[7, dose_cols] <- x$nipde_by_dose

  tab[8, "Total"] <- x$total_admin_mean
  tab[9, "Total"] <- x$total_unique_mean
  tab[10, "Total"] <- x$early_stop_pct
  tab[11, "Duration"] <- x$duration_mean

  num_cols <- setdiff(names(tab), "Metric")
  tab[num_cols] <- lapply(tab[num_cols], function(z) round(z, digits))

  cat("\n")
  cat("------------------------------------------------------------\n")
  cat("Scenario:", scenario, "\n")
  cat("Method:", method, "\n")
  cat("AIDE design:", design, "\n")
  cat("Decision method:", x$decision_method, "\n")
  cat("Final MTD method:", x$mtd_method, "\n")
  cat("Restrict to tried:", x$restrict_to_tried, "\n")
  cat("True IPDE alpha:", true_ipde_alpha, "\n")
  cat("Discount alpha:", discount_alpha, "\n")
  cat("r_carry:", x$r_carry, "\n")
  cat("Allow suspend label:", allow_suspend, "\n")
  cat("Accrual rate:", accrual, "\n")
  cat("ntrial:", x$ntrial, "\n")
  cat("------------------------------------------------------------\n")

  print(tab, row.names = FALSE)

  invisible(tab)
}

## ------------------------------------------------------------
## Fallback settings
## If these objects already exist in your workspace, they will not be overwritten.
## ------------------------------------------------------------

if (!exists("target_BOIN")) target_BOIN <- 0.3
if (!exists("cutoff_equiv")) cutoff_equiv <- 0.95
if (!exists("ntrial")) ntrial <- 3000
if (!exists("seed0")) seed0 <- 1

if (!exists("T_assess_equiv")) T_assess_equiv <- 3
if (!exists("C_equiv")) C_equiv <- 3L
if (!exists("cycle_max_equiv")) cycle_max_equiv <- 2L
if (!exists("Nmax_eff_equiv")) Nmax_eff_equiv <- 30L
if (!exists("ipde_design_equiv")) ipde_design_equiv <- 2L

if (!exists("t0")) t0 <- 0
if (!exists("day_obs")) day_obs <- 0
if (!exists("dlt_dist_equiv")) dlt_dist_equiv <- 2
if (!exists("dlt_alpha_equiv")) dlt_alpha_equiv <- 0.5
if (!exists("arrival_type_label")) arrival_type_label <- "on-demand exponential"
if (!exists("allow_suspend_equiv")) allow_suspend_equiv <- NA
if (!exists("store_raw")) store_raw <- FALSE
if (!exists("verbose")) verbose <- FALSE

if (!exists("alpha_true_list")) alpha_true_list <- c(0)
if (!exists("accrual_list")) accrual_list <- c(1 / 2)

## Use your real scenarios object if already defined.
## Fallback is only for a quick syntax/run check.
if (!exists("scenarios")) {
  scenarios <- rbind(
    c(0.07, 0.12, 0.17, 0.22, 0.30),
    c(0.05, 0.10, 0.18, 0.30, 0.40),
    c(0.15, 0.20, 0.30, 0.35, 0.45),
    c(0.15, 0.30, 0.38, 0.45, 0.55),
    c(0.30, 0.35, 0.40, 0.45, 0.50),
    c(0.50, 0.55, 0.60, 0.65, 0.70)
  )
}

## ------------------------------------------------------------
## Runner
## ------------------------------------------------------------

outdir_aide <- "oc_results_equiv_AIDE"
if (!dir.exists(outdir_aide)) dir.create(outdir_aide, recursive = TRUE)

all_res_aide <- list()

## Prior settings aligned with methods_prior.R.
q_skeleton_prior <- c(0.15, 0.20, 0.30, 0.35, 0.45)
dose_alpha_mg_prior <- c(15, 20, 30, 35, 45)
dose_ipcrm_prior <- c(15, 20, 30, 35, 45)
dose_ipcrm_prior <- dose_ipcrm_prior / (2 * stats::sd(dose_ipcrm_prior))

crm_skeleton_aide <- q_skeleton_prior
crm_dose_values_aide <- dose_alpha_mg_prior
crm_dose_scores_aide <- dose_ipcrm_prior
crm_theta_prior_mean_aide <- 0
crm_theta_prior_sd_aide <- sqrt(2)
crm_cumu_beta0_mean_aide <- -2.8
crm_cumu_beta0_prec_aide <- 2
crm_cumu_beta0_df_aide <- 1
crm_cumu_beta1_shape_aide <- 2.5
crm_cumu_beta1_rate_aide <- 1.6
crm_cumu_beta2_rate_aide <- 1

cfo_skeleton_aide <- c(0.005, 0.01, 0.05, 0.1, 0.3)
cfo_sigma2_beta_aide <- 30
cfo_eta_aide <- 1
cfo_model_file_aide <- "PRIDE.bug"
cfo_pk_method_aide <- "approx"
cfo_n_mc_w_aide <- 200
cfo_m_use_aide <- 1000

model_list_aide <- c("BOIN", "CFO")

## old cap used in your original AIDE code
d_cap_aide <- 100

## new escalation gate: no escalation until current dose has >= 3 treated
dose_cap_aide <- 3L

## Final MTD selection gate.
## TRUE: select among tried/non-eliminated doses.
## FALSE: allow all non-eliminated doses to enter final selection.
restrict_to_tried_aide <- TRUE

## If you only want one method, edit these vectors.
## BOIN options: c("boin", "approx1", "approx2")
method_list_boin_aide <- c("approx1", "approx2")
## CFO options: c("empirical", "pride")
cfo_method_list_aide <- c("empirical", "pride")
alpha_true_list <- c(0)
## Discount value used by approx1 / approx2.
## Set to 0.1 for your usual working setting.
r_carry_aide <- 0
ntrial = 3000
for (model_aide in model_list_aide) {
  method_loop_aide <- if (model_aide == "BOIN") method_list_boin_aide else cfo_method_list_aide

  for (method_aide in method_loop_aide) {
    for (accrual in accrual_list) {
      for (alpha_true in alpha_true_list) {
        if (model_aide == "BOIN") {
          decision_method_aide <- method_aide
          mtd_method_aide <- method_aide
          cfo_method_aide <- "empirical"
          method_tag_aide <- method_aide
          method_family_aide <- "AIDE-BOIN"
        } else {
          decision_method_aide <- "boin"
          mtd_method_aide <- NULL
          cfo_method_aide <- method_aide
          method_tag_aide <- paste0("cfo_", cfo_method_aide)
          method_family_aide <- "AIDE-CFO"
        }

      method_label_aide <- paste0(
        "AIDE-equiv",
        "-", method_tag_aide,
        "-w", T_assess_equiv,
        "-c", C_equiv,
        "-cyc", cycle_max_equiv,
        "-rate", accrual,
        "-dosecap", dose_cap_aide,
        "-tried", as.integer(restrict_to_tried_aide),
        "-rcarry", r_carry_aide
      )

      outfile <- file.path(
        outdir_aide,
        paste0(
          "aide_equiv",
          "_model", model_aide,
          "_method", method_tag_aide,
          "_a", fmt_param(alpha_true),
          "_w", T_assess_equiv,
          "_c", C_equiv,
          "_cyc", cycle_max_equiv,
          "_rate", fmt_param(accrual),
          "_dosecap", dose_cap_aide,
          "_tried", as.integer(restrict_to_tried_aide),
          "_rcarry", fmt_param(r_carry_aide),
          "_n", ntrial,
          ".txt"
        )
      )

      with_sink(outfile, {

        cat("========================================\n")
        cat(method_family_aide, "equivalence check\n")
        cat("========================================\n")
        cat("model =", model_aide, "\n")
        cat("method =", method_label_aide, "\n")
        cat("decision_method =", ifelse(model_aide == "BOIN", decision_method_aide, method_tag_aide), "\n")
        cat("mtd_method =", ifelse(is.null(mtd_method_aide), method_tag_aide, mtd_method_aide), "\n")
        cat("target =", target_BOIN, "\n")
        cat("alpha_true =", alpha_true, "\n")
        cat("r_carry =", r_carry_aide, "\n")
        cat("ntrial =", ntrial, "\n")
        cat("seed0 =", seed0, "\n")
        cat("settings:",
            "cohort=", C_equiv,
            "Nmax=", Nmax_eff_equiv,
            "cycle_max=", cycle_max_equiv,
            "window=", T_assess_equiv,
            "arrival=", arrival_type_label,
            "rate=", accrual,
            "d.cap=", d_cap_aide,
            "dose_cap=", dose_cap_aide,
            "restrict_to_tried=", restrict_to_tried_aide,
            "design=", ipde_design_equiv,
            "\n")
        if (model_aide == "CFO") {
          cat("CFO skeleton:", paste(cfo_skeleton_aide, collapse = ", "), "\n")
          cat("CFO sigma2_beta:", cfo_sigma2_beta_aide, "\n")
          cat("CFO eta:", cfo_eta_aide, "\n")
          cat("CFO model file:", cfo_model_file_aide, "\n")
          cat("CFO pk method:", cfo_pk_method_aide, "\n")
          cat("CFO m_use:", cfo_m_use_aide, "\n")
        }
        cat("outfile =", outfile, "\n")
        cat("========================================\n\n")

        for (i in seq_len(nrow(scenarios))) {

          p_base <- as.numeric(scenarios[i, ])

          p_ipde <- p_base
          if (length(p_ipde) >= 2) {
            p_ipde[-1] <- pmin(
              p_base[-1] + alpha_true * p_base[-length(p_base)],
              1
            )
          }

          res <- get_oc_sim_AIDE(
            target = target_BOIN,
            p.true = p_base,
            p.true_ipde = p_ipde,

            ntrial = ntrial,
            seed = seed0,

            model = model_aide,
            ipde_design = ipde_design_equiv,

            N_pat = Nmax_eff_equiv,
            Nmax_eff = Nmax_eff_equiv,
            C = C_equiv,
            T_assess = T_assess_equiv,
            cycle_max = cycle_max_equiv,

            arrival_rate = accrual,
            t0 = t0,

            cutoff = cutoff_equiv,

            ## old stopping cap
            d.cap = d_cap_aide,

            ## new escalation gate
            dose_cap = dose_cap_aide,

            day_obs = day_obs,

            dlt_dist = dlt_dist_equiv,
            dlt_alpha = dlt_alpha_equiv,

            decision_method = decision_method_aide,
            mtd_method = mtd_method_aide,
            restrict_to_tried = restrict_to_tried_aide,
            r_carry = r_carry_aide,
            store_raw = store_raw,
            verbose = verbose,
            crm_skeleton = crm_skeleton_aide,
            crm_dose_values = crm_dose_values_aide,
            crm_dose_scores = crm_dose_scores_aide,
            crm_theta_prior_mean = crm_theta_prior_mean_aide,
            crm_theta_prior_sd = crm_theta_prior_sd_aide,
            crm_cumu_beta0_mean = crm_cumu_beta0_mean_aide,
            crm_cumu_beta0_prec = crm_cumu_beta0_prec_aide,
            crm_cumu_beta0_df = crm_cumu_beta0_df_aide,
            crm_cumu_beta1_shape = crm_cumu_beta1_shape_aide,
            crm_cumu_beta1_rate = crm_cumu_beta1_rate_aide,
            crm_cumu_beta2_rate = crm_cumu_beta2_rate_aide,
            cfo_method = cfo_method_aide,
            cfo_skeleton = cfo_skeleton_aide,
            cfo_model_file = cfo_model_file_aide,
            cfo_sigma2_beta = cfo_sigma2_beta_aide,
            cfo_eta = cfo_eta_aide,
            cfo_pk_method = cfo_pk_method_aide,
            cfo_n_mc_w = cfo_n_mc_w_aide,
            cfo_m_use = cfo_m_use_aide
          )

          key <- paste(
            model_aide,
            method_tag_aide,
            "alpha", alpha_true,
            "accrual", accrual,
            "scenario", i,
            sep = "_"
          )
          all_res_aide[[key]] <- res

          print_oc_table_aide(
            x = res,
            scenario = i,
            design = ipde_design_equiv,
            method = method_label_aide,
            true_ipde_alpha = alpha_true,
            discount_alpha = ifelse(model_aide == "BOIN" && decision_method_aide != "boin", r_carry_aide, 0),
            allow_suspend = allow_suspend_equiv,
            accrual = accrual,
            digits = 2
          )
        }
      })

      cat("Finished AIDE:", outfile, "\n")
    }
  }
}
}

saveRDS(
  all_res_aide,
  file = file.path(outdir_aide, "all_res_aide.rds")
)
