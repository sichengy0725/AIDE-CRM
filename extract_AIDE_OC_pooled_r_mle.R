## ============================================================
## Extract AIDE OC from cluster jobs using fixed folder/file names
## No recursive file searching
##
## Works with run_oc_AIDE_cluster.R output:
##   folder:
##     SC1-model-BOIN-method-approx1-...
##   file:
##     SC_1_modelBOIN_methodapprox1_alphaTrue0_rcarry0_rate0p5-job-1-combined.rds
##
## For CRM, Estimated pj uses CRM posterior estimate of p_j.
## ============================================================

rm(list = ls())

## -------------------------------
## Settings
## -------------------------------

## setwd("/rsrch8/home/biostatistics/syang10/AIDE")

scenario_id_list <- c(1, 2, 3, 4, 5, 6)

jobs.expected <- 1:2000
ntrial.expected <- 2000

target <- 0.30

T_assess <- 3
C <- 3L
cycle_max <- 2L
Nmax_eff <- 30L
dose_cap <- 3L

alpha_true_list <- c(0)
arrival_rate_list <- c(1 / 2)

## BOIN settings
boin_method_list <- c("boin", "approx1", "approx2")
boin_r_carry_list <- c(0, 0.1, 0.2, 0.4)

## CRM settings
crm_r_model_list <- c("random")
crm_fixed_r_list <- c(0, 0.1, 0.2)

crm_skeleton <- c(0.15, 0.20, 0.30, 0.35, 0.45)
crm_alpha_sd <- 2
crm_a_r <- 1
crm_b_r <- 9

## Choose models to extract
# model_list <- c("BOIN", "CRM")
model_list <- c("CRM")
## model_list <- c("BOIN")

out.dir <- "OC_summary_from_parallel_AIDE_jobs"
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

## -------------------------------
## Helpers
## -------------------------------

fmt_num <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return("NA")
  gsub("\\.", "p", as.character(x))
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

make_p_ipde <- function(p_base, alpha_true) {
  p_ipde <- p_base
  
  if (length(p_ipde) >= 2L) {
    p_ipde[-1] <- pmin(
      p_base[-1] + alpha_true * p_base[-length(p_base)],
      1
    )
  }
  
  p_ipde
}

make_crm_tag <- function(model, crm_r_model) {
  if (model == "BOIN") {
    return("none")
  }
  
  paste0(
    "skeleton-", paste(fmt_num(crm_skeleton), collapse = "_"),
    "-alphasd-", fmt_num(crm_alpha_sd),
    "-ar-", fmt_num(crm_a_r),
    "-br-", fmt_num(crm_b_r)
  )
}

make_method_tag <- function(model, boin_method = NULL, crm_r_model = NULL) {
  if (model == "BOIN") {
    return(boin_method)
  }
  
  paste0("crm_", crm_r_model)
}

make_foldername <- function(sc,
                            model,
                            method_tag,
                            crm_tag,
                            alpha_true,
                            r_carry,
                            arrival_rate) {
  
  paste0(
    "SC", sc,
    "-model-", model,
    "-method-", method_tag,
    "-crm-", crm_tag,
    "-target-", fmt_num(target),
    "-alphaTrue-", fmt_num(alpha_true),
    "-rcarry-", fmt_num(r_carry),
    "-w-", T_assess,
    "-c-", C,
    "-cyc-", cycle_max,
    "-rate-", fmt_num(arrival_rate),
    "-Nmax-", Nmax_eff,
    "-dosecap-", dose_cap
  )
}

make_filename <- function(sc,
                          model,
                          method_tag,
                          alpha_true,
                          r_carry,
                          arrival_rate,
                          job) {
  
  group_key <- paste(
    "SC", sc,
    paste0("model", model),
    paste0("method", method_tag),
    paste0("alphaTrue", fmt_num(alpha_true)),
    paste0("rcarry", fmt_num(r_carry)),
    paste0("rate", fmt_num(arrival_rate)),
    sep = "_"
  )
  
  paste0(group_key, "-job-", job, "-combined.rds")
}

## Estimate pj extraction.
## For BOIN, this is pj_iso_mean.
## For CRM, this should be the posterior mean estimate of current-dose p_j.
## In the modified get_oc_sim_AIDE(), CRM estimates are stored in pj_iso_mean /
## pj_iso_by_trial because the wrapper saves fit$final$pj_iso or fit$final$phat.
extract_est_pj_mat <- function(res.list, model) {
  
  mat.list <- lapply(res.list, function(x) {
    
    if (model == "CRM") {
      ## Preferred names if you later store CRM p_j under a clearer name.
      z <- get_field(
        x,
        candidates = c(
          "crm_pj_by_trial",
          "p_hat_by_trial",
          "phat_by_trial",
          "pj_by_trial",
          "pj_iso_by_trial"
        ),
        required = FALSE
      )
      
      if (!is.null(z)) return(as.matrix(z))
      
      zmean <- get_field(
        x,
        candidates = c(
          "crm_pj_mean",
          "p_hat_mean",
          "phat_mean",
          "pj_mean",
          "pj_iso_mean"
        ),
        required = FALSE
      )
      
      if (!is.null(zmean)) {
        return(matrix(as.numeric(zmean), nrow = 1L))
      }
      
      stop("Cannot find CRM p_j estimates in one result object.")
    }
    
    ## BOIN / approx1 / approx2
    z <- get_field(
      x,
      candidates = c("pj_iso_by_trial"),
      required = FALSE
    )
    
    if (!is.null(z)) return(as.matrix(z))
    
    zmean <- get_field(
      x,
      candidates = c("pj_iso_mean"),
      required = TRUE
    )
    
    matrix(as.numeric(zmean), nrow = 1L)
  })
  
  do.call(rbind, mat.list)
}

extract_mean_metric <- function(res.list, field_candidates, ndose = NULL) {
  vals <- lapply(res.list, get_field, candidates = field_candidates)
  
  ntrial <- vapply(res.list, function(x) x$ntrial, numeric(1))
  ntrial.total <- sum(ntrial)
  
  if (is.null(ndose)) {
    ndose <- length(vals[[1]])
  }
  
  out <- rep(0, ndose)
  
  for (i in seq_along(vals)) {
    out <- out + as.numeric(vals[[i]]) * ntrial[i]
  }
  
  out / ntrial.total
}

## Build one OC summary from all job-level combined RDS files
summarize_aide_files <- function(files.use,
                                 sc,
                                 model,
                                 method_tag,
                                 crm_r_model,
                                 alpha_true,
                                 r_carry,
                                 arrival_rate) {
  
  res.list <- lapply(files.use, readRDS)
  
  ndose <- get_field(res.list[[1]], c("ndose"))
  ntrial.total <- sum(vapply(res.list, function(x) x$ntrial, numeric(1)))
  
  p.true <- as.numeric(get_field(res.list[[1]], c("p.true")))
  p.true_ipde <- as.numeric(get_field(res.list[[1]], c("p.true_ipde")))
  
  final_MTD <- unlist(
    lapply(res.list, get_field, candidates = c("final_MTD")),
    use.names = FALSE
  )
  
  est_pj_mat <- extract_est_pj_mat(res.list, model = model)
  
  est_pj <- apply(est_pj_mat, 2, function(z) {
    if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE)
  })
  
  selection_pct <- calc_select_rate_pct(
    sel = final_MTD,
    ndose = ndose,
    denom = length(final_MTD)
  )
  
  early_stop_pct <- 100 * sum(final_MTD == 99L, na.rm = TRUE) / length(final_MTD)
  
  n_by_dose <- extract_mean_metric(
    res.list,
    field_candidates = c("n_by_dose"),
    ndose = ndose
  )
  
  unique_n_by_dose <- extract_mean_metric(
    res.list,
    field_candidates = c("unique_n_by_dose"),
    ndose = ndose
  )
  
  nipde_by_dose <- extract_mean_metric(
    res.list,
    field_candidates = c("nipde_by_dose"),
    ndose = ndose
  )
  
  total_admin <- unlist(
    lapply(res.list, get_field, candidates = c("total_admin")),
    use.names = FALSE
  )
  
  total_unique <- unlist(
    lapply(res.list, get_field, candidates = c("total_unique")),
    use.names = FALSE
  )
  
  duration <- unlist(
    lapply(res.list, get_field, candidates = c("duration")),
    use.names = FALSE
  )
  
  total_admin_mean <- mean(total_admin, na.rm = TRUE)
  total_unique_mean <- mean(total_unique, na.rm = TRUE)
  duration_mean <- mean(duration, na.rm = TRUE)
  
  ## Long dose-level summary, easier for plotting / merging.
  dose_summary <- data.frame(
    Scenario = sc,
    Model = model,
    Method = method_tag,
    CRM_r_model = ifelse(model == "CRM", crm_r_model, NA_character_),
    Alpha_true = alpha_true,
    r_carry = r_carry,
    Accrual = arrival_rate,
    Dose = seq_len(ndose),
    
    True_DLT_rate = p.true,
    True_IPDE_DLT_rate = p.true_ipde,
    
    ## For BOIN this is pj_iso; for CRM this is CRM posterior estimate of p_j.
    Estimated_pj = est_pj,
    
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
    stringsAsFactors = FALSE
  )
  
  ## Printed-table-style summary, matching your txt output.
  dose_cols <- paste0("D", seq_len(ndose))
  
  table_summary <- data.frame(
    Scenario = sc,
    Model = model,
    Method = method_tag,
    CRM_r_model = ifelse(model == "CRM", crm_r_model, NA_character_),
    Alpha_true = alpha_true,
    r_carry = r_carry,
    Accrual = arrival_rate,
    Metric = c(
      "True DLT rate",
      "True IPDE DLT rate",
      ifelse(model == "CRM", "Estimated CRM pj", "Estimated pj_iso"),
      "MTD Selection %",
      "# Pts Treated",
      "# Unique Pts by Dose",
      "# IPDE Doses",
      "# Total Administrations",
      "# Total Unique Patients",
      "% Early Stopping",
      "Duration"
    ),
    matrix(NA_real_, nrow = 11L, ncol = ndose),
    Total = NA_real_,
    Duration = NA_real_,
    n_valid = length(final_MTD),
    ntrial_from_files = ntrial.total,
    stringsAsFactors = FALSE
  )
  
  names(table_summary)[8:(7 + ndose)] <- dose_cols
  
  table_summary[1, dose_cols] <- p.true
  table_summary[2, dose_cols] <- p.true_ipde
  table_summary[3, dose_cols] <- est_pj
  table_summary[4, dose_cols] <- selection_pct
  table_summary[5, dose_cols] <- n_by_dose
  table_summary[6, dose_cols] <- unique_n_by_dose
  table_summary[7, dose_cols] <- nipde_by_dose
  
  table_summary[8, "Total"] <- total_admin_mean
  table_summary[9, "Total"] <- total_unique_mean
  table_summary[10, "Total"] <- early_stop_pct
  table_summary[11, "Duration"] <- duration_mean
  
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

idx <- 1L
miss.idx <- 1L

for (model in model_list) {
  
  if (model == "BOIN") {
    method_loop <- boin_method_list
    crm_loop <- NA_character_
  } else {
    method_loop <- "crm"
    crm_loop <- crm_r_model_list
  }
  
  for (sc in scenario_id_list) {
    for (alpha_true in alpha_true_list) {
      for (arrival_rate in arrival_rate_list) {
        
        for (crm_r_model in crm_loop) {
          
          if (model == "BOIN") {
            r_loop <- boin_r_carry_list
          } else if (model == "CRM" && crm_r_model == "fixed") {
            r_loop <- crm_fixed_r_list
          } else if (model == "CRM" && crm_r_model == "random") {
            ## r_carry is only a label for random CRM in the current runner.
            ## Keep this consistent with the r_carry_list used when running.
            r_loop <- c(0)
          } else {
            stop("Unknown model / crm_r_model combination.")
          }
          
          for (method0 in method_loop) {
            for (r_carry in r_loop) {
              
              method_tag <- make_method_tag(
                model = model,
                boin_method = if (model == "BOIN") method0 else NULL,
                crm_r_model = if (model == "CRM") crm_r_model else NULL
              )
              
              crm_tag <- make_crm_tag(
                model = model,
                crm_r_model = crm_r_model
              )
              
              foldername <- make_foldername(
                sc = sc,
                model = model,
                method_tag = method_tag,
                crm_tag = crm_tag,
                alpha_true = alpha_true,
                r_carry = r_carry,
                arrival_rate = arrival_rate
              )
              
              files <- file.path(
                foldername,
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
                      job = j
                    )
                  },
                  character(1)
                )
              )
              
              exists.vec <- file.exists(files)
              files.use <- files[exists.vec]
              missing.jobs <- jobs.expected[!exists.vec]
              
              cat("\n====================================\n")
              cat("Scenario:", sc, "\n")
              cat("Model:", model, "\n")
              cat("Method:", method_tag, "\n")
              cat("CRM r model:", ifelse(model == "CRM", crm_r_model, NA), "\n")
              cat("r_carry:", r_carry, "\n")
              cat("Accrual:", arrival_rate, "\n")
              cat("Folder:", foldername, "\n")
              cat("Found", length(files.use), "files; missing", length(missing.jobs), "jobs.\n")
              cat("====================================\n")
              
              missing.log[[miss.idx]] <- data.frame(
                Scenario = sc,
                Model = model,
                Method = method_tag,
                CRM_r_model = ifelse(model == "CRM", crm_r_model, NA_character_),
                Alpha_true = alpha_true,
                r_carry = r_carry,
                Accrual = arrival_rate,
                Folder = foldername,
                n_found = length(files.use),
                n_missing = length(missing.jobs),
                missing_jobs = paste(missing.jobs, collapse = ","),
                stringsAsFactors = FALSE
              )
              miss.idx <- miss.idx + 1L
              
              if (length(files.use) == 0L) {
                warning("No files found for folder: ", foldername)
                next
              }
              
              t.read <- Sys.time()
              
              one <- summarize_aide_files(
                files.use = files.use,
                sc = sc,
                model = model,
                method_tag = method_tag,
                crm_r_model = crm_r_model,
                alpha_true = alpha_true,
                r_carry = r_carry,
                arrival_rate = arrival_rate
              )
              
              cat("read/summarize time:",
                  round(difftime(Sys.time(), t.read, units = "secs"), 2),
                  "sec\n")
              
              if (one$n_valid != ntrial.expected) {
                cat(
                  "Warning: expected", ntrial.expected,
                  "trials but found", one$n_valid, "\n"
                )
              }
              
              key <- paste(
                "SC", sc,
                model,
                method_tag,
                paste0("alpha", fmt_num(alpha_true)),
                paste0("r", fmt_num(r_carry)),
                paste0("rate", fmt_num(arrival_rate)),
                sep = "_"
              )
              
              all.dose.summary[[key]] <- one$dose_summary
              all.table.summary[[key]] <- one$table_summary
              
              idx <- idx + 1L
            }
          }
        }
      }
    }
  }
}

if (length(all.dose.summary) == 0L) {
  stop("No files were found. Check folder names, filenames, and working directory.")
}

all.dose.summary.df <- do.call(rbind, all.dose.summary)
all.table.summary.df <- do.call(rbind, all.table.summary)
missing.log.df <- do.call(rbind, missing.log)

## Round numeric columns for readable CSV output.
round_numeric_df <- function(dat, digits = 4) {
  num_cols <- vapply(dat, is.numeric, logical(1))
  dat[num_cols] <- lapply(dat[num_cols], round, digits)
  dat
}

all.dose.summary.out <- round_numeric_df(all.dose.summary.df, digits = 4)
all.table.summary.out <- round_numeric_df(all.table.summary.df, digits = 4)

out.tag <- paste0(
  "All_AIDE_OC",
  "_models", paste(model_list, collapse = "_"),
  "_target", fmt_num(target),
  "_w", T_assess,
  "_c", C,
  "_cyc", cycle_max,
  "_Nmax", Nmax_eff,
  "_dosecap", dose_cap,
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