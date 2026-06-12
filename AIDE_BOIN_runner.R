## ============================================================
## AIDE-BOIN runner and OC table printer
## Source this after defining scenarios / global simulation settings,
## or edit the fallback settings below.
## ============================================================
source("AIDE.R")

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

model_aide <- "BOIN"

## old cap used in your original AIDE code
d_cap_aide <- 100

## new escalation gate: no escalation until current dose has >= 3 treated
dose_cap_aide <- 3L

## Run all three versions.
## If you only want one method, set method_list_aide <- c("boin"), c("approx1"), or c("approx2").
method_list_aide <- c("approx1", "approx2")
alpha_true_list <- c(0)
## Discount value used by approx1 / approx2.
## Set to 0.1 for your usual working setting.
r_carry_aide <- 0
ntrial = 3000
for (decision_method_aide in method_list_aide) {
  for (accrual in accrual_list) {
    for (alpha_true in alpha_true_list) {
      mtd_method_aide <- decision_method_aide

      method_label_aide <- paste0(
        "AIDE-equiv",
        "-", decision_method_aide,
        "-w", T_assess_equiv,
        "-c", C_equiv,
        "-cyc", cycle_max_equiv,
        "-rate", accrual,
        "-dosecap", dose_cap_aide,
        "-rcarry", r_carry_aide
      )

      outfile <- file.path(
        outdir_aide,
        paste0(
          "aide_equiv",
          "_method", decision_method_aide,
          "_a", fmt_param(alpha_true),
          "_w", T_assess_equiv,
          "_c", C_equiv,
          "_cyc", cycle_max_equiv,
          "_rate", fmt_param(accrual),
          "_dosecap", dose_cap_aide,
          "_rcarry", fmt_param(r_carry_aide),
          "_n", ntrial,
          ".txt"
        )
      )

      with_sink(outfile, {

        cat("========================================\n")
        cat("AIDE-BOIN equivalence check\n")
        cat("========================================\n")
        cat("method =", method_label_aide, "\n")
        cat("decision_method =", decision_method_aide, "\n")
        cat("mtd_method =", mtd_method_aide, "\n")
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
            "design=", ipde_design_equiv,
            "\n")
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
            r_carry = r_carry_aide,

            store_raw = store_raw,
            verbose = verbose
          )

          key <- paste(
            decision_method_aide,
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
            discount_alpha = ifelse(decision_method_aide == "boin", 0, r_carry_aide),
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

saveRDS(
  all_res_aide,
  file = file.path(outdir_aide, "all_res_aide.rds")
)
