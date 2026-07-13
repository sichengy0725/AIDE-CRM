## ============================================================
## Extract TITE AIDE-CRM operating characteristics from cluster jobs
## Matched to run_oc_AIDE_TITE.R
##
## This extractor discovers all combined TITE result files under the
## configured results root, groups them by scenario and design setting, and
## aggregates across LSF job indices.  It writes the same core metrics as the
## standard AIDE extractor: dose selection, allocation, unique patients,
## IPDE administrations, stopping, and duration.
## ============================================================

rm(list = ls())

## -------------------------------
## User settings
## -------------------------------

## Run from the AIDE project directory, or set this before running the file.
## setwd("/rsrch8/home/biostatistics/syang10/AIDE")

## Defaults match the fixed 5-dose TITE runner.  Environment variables make
## it easy to use a different scenario set or results directory on the cluster.
scenario_set_name <- Sys.getenv(
  "AIDE_TITE_SCENARIO_SET",
  unset = "Set_5dose_adaptive_r_37"
)

default_results_root <- paste0("oc_results_cluster_AIDE_TITE_", scenario_set_name)
results_root <- Sys.getenv(
  "AIDE_TITE_RESULTS_ROOT",
  unset = default_results_root
)

## Expected LSF job indices.  Used only for the missing-job report.
## Accepted examples: "1:2000", "1,2,5:10", or "" to skip the check.
expected_jobs_spec <- Sys.getenv("AIDE_TITE_EXPECTED_JOBS", unset = "1:2000")

default_out_dir <- paste0("OC_summary_from_parallel_AIDE_TITE_jobs_", scenario_set_name)
out.dir <- Sys.getenv("AIDE_TITE_SUMMARY_DIR", unset = default_out_dir)

## -------------------------------
## Helpers
## -------------------------------

parse_job_spec <- function(x) {
  x <- trimws(x)
  if (!nzchar(x)) return(integer(0))

  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  out <- unlist(lapply(parts, function(part) {
    if (grepl("^[0-9]+:[0-9]+$", part)) {
      ends <- as.integer(strsplit(part, ":", fixed = TRUE)[[1]])
      return(seq.int(min(ends), max(ends)))
    }
    if (grepl("^[0-9]+$", part)) return(as.integer(part))
    stop("Invalid AIDE_TITE_EXPECTED_JOBS component: ", part)
  }), use.names = FALSE)

  sort(unique(as.integer(out)))
}

fmt_short <- function(x, digits = 2) {
  x <- round(as.numeric(x), digits)
  out <- format(x, scientific = FALSE, trim = TRUE)
  out <- sub("(\\.[0-9]*?)0+$", "\\1", out)
  out <- sub("\\.$", "", out)
  out[out == "-0"] <- "0"
  out[out == ""] <- "0"
  out <- gsub("-", "m", out)
  gsub("\\.", "p", out)
}

decode_number <- function(x) {
  x <- gsub("p", ".", x, fixed = TRUE)
  x <- gsub("m", "-", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

collapse_integer_ranges <- function(x) {
  x <- sort(unique(as.integer(x[!is.na(x)])))
  if (length(x) == 0L) return("")

  starts <- x[c(TRUE, diff(x) != 1L)]
  ends <- x[c(diff(x) != 1L, TRUE)]
  paste(ifelse(starts == ends, starts, paste0(starts, "-", ends)), collapse = ",")
}

get_field <- function(x, candidates, required = TRUE) {
  for (nm in candidates) {
    if (!is.null(x[[nm]])) return(x[[nm]])
  }
  if (required) {
    stop("Missing required result field: ", paste(candidates, collapse = ", "))
  }
  NULL
}

get_scalar_field <- function(x, candidates, default = NA) {
  z <- get_field(x, candidates, required = FALSE)
  if (is.null(z) || length(z) == 0L) default else z[[1L]]
}

parse_combined_file <- function(path) {
  file <- basename(path)
  pattern <- paste0(
    "^(.*)_SC([0-9]+)_([^_]+)_(.*)_a([^_]+)_r([^_]+)_rate([^_]+)",
    "_cyc([^_]+)_Nmax([^_]+)_cont([01])_tried([01])(_ptarget1)?",
    "-job-([0-9]+)-combined\\.rds$"
  )
  hit <- regexec(pattern, file, perl = TRUE)
  vals <- regmatches(file, hit)[[1L]]
  if (length(vals) == 0L) return(NULL)

  folder <- basename(dirname(path))
  folder_hit <- regexec("-w-([^-]+)-c-([^-]+)-cyc-", folder, perl = TRUE)
  folder_vals <- regmatches(folder, folder_hit)[[1L]]

  data.frame(
    file = normalizePath(path, winslash = "/", mustWork = FALSE),
    folder = folder,
    ## T_assess and cohort size are encoded in the parent folder rather than
    ## the filename, so include that folder when forming an aggregation key.
    group_key = paste(
      normalizePath(dirname(path), winslash = "/", mustWork = FALSE),
      sub("-job-[0-9]+-combined\\.rds$", "", file),
      sep = "::"
    ),
    Scenario_Set = vals[2L],
    Scenario = as.integer(vals[3L]),
    Model = vals[4L],
    Method = vals[5L],
    Alpha_true = decode_number(vals[6L]),
    r_carry = decode_number(vals[7L]),
    Accrual = decode_number(vals[8L]),
    Cycle_Max = as.integer(decode_number(vals[9L])),
    Nmax_eff = as.integer(decode_number(vals[10L])),
    Continuous_Enrollment = as.integer(vals[11L]),
    Restrict_To_Tried = as.integer(vals[12L]),
    Restrict_To_Target = as.integer(!is.na(vals[13L]) && nzchar(vals[13L])),
    Job = as.integer(vals[14L]),
    T_assess = if (length(folder_vals) > 0L) decode_number(folder_vals[2L]) else NA_real_,
    Cohort_Size = if (length(folder_vals) > 0L) as.integer(decode_number(folder_vals[3L])) else NA_integer_,
    stringsAsFactors = FALSE
  )
}

trial_vector <- function(res, candidates, ntrial) {
  z <- get_field(res, candidates, required = FALSE)
  if (is.null(z)) return(rep(NA_real_, ntrial))

  z <- as.numeric(z)
  if (length(z) == ntrial) return(z)
  if (length(z) == 1L) return(rep(z, ntrial))

  stop(
    "Field ", candidates[1L], " has length ", length(z),
    "; expected 1 or ntrial = ", ntrial, "."
  )
}

weighted_dose_mean <- function(results, candidates, ndose) {
  numerator <- rep(0, ndose)
  denominator <- 0

  for (res in results) {
    z <- get_field(res, candidates, required = FALSE)
    if (is.null(z)) next
    z <- as.numeric(z)
    if (length(z) != ndose) {
      stop(
        "Field ", candidates[1L], " has length ", length(z),
        "; expected ndose = ", ndose, "."
      )
    }
    ntrial <- as.numeric(res$ntrial)
    numerator <- numerator + z * ntrial
    denominator <- denominator + ntrial
  }

  if (denominator == 0) rep(NA_real_, ndose) else numerator / denominator
}

weighted_estimated_p <- function(results, ndose) {
  candidates <- c(
    "p_hat_mean", "phat_mean", "pj_mean", "pj_iso_mean",
    "p_hat", "phat", "pj_iso"
  )
  weighted_dose_mean(results, candidates, ndose)
}

scenario_metadata <- function(res, scenario_id, scenario_set) {
  list(
    Scenario_Set = scenario_set,
    Scenario_Name = as.character(get_scalar_field(
      res, c("scenario_name"), paste0(scenario_set, "_SC", scenario_id)
    )),
    Source_Scenario = get_scalar_field(res, c("source_scenario"), NA_integer_),
    True_MTD = get_scalar_field(res, c("true_mtd"), NA_integer_),
    Scenario_Attempt = get_scalar_field(res, c("scenario_attempt"), NA_integer_),
    p.true = as.numeric(get_field(res, c("p.true"))),
    p.true_ipde = as.numeric(get_field(res, c("p.true_ipde")))
  )
}

summarize_tite_group <- function(meta, results) {
  first <- results[[1L]]
  ndose <- as.integer(first$ndose)
  ntrial.total <- sum(vapply(results, function(x) as.numeric(x$ntrial), numeric(1)))

  final_mtd <- unlist(lapply(results, function(res) {
    trial_vector(res, c("final_MTD", "final_mtd_by_trial"), as.integer(res$ntrial))
  }), use.names = FALSE)

  if (length(final_mtd) != ntrial.total) {
    stop("final_MTD length does not match total ntrial for ", meta$group_key[1L])
  }

  scenario <- scenario_metadata(
    res = first,
    scenario_id = meta$Scenario[1L],
    scenario_set = meta$Scenario_Set[1L]
  )
  if (length(scenario$p.true) != ndose || length(scenario$p.true_ipde) != ndose) {
    stop("Stored true toxicity vectors do not match ndose for ", meta$group_key[1L])
  }

  selection_pct <- vapply(
    seq_len(ndose),
    function(dose) 100 * sum(final_mtd == dose, na.rm = TRUE) / ntrial.total,
    numeric(1)
  )

  n_by_dose <- weighted_dose_mean(results, c("n_by_dose", "npatients"), ndose)
  unique_n_by_dose <- weighted_dose_mean(
    results, c("unique_n_by_dose", "nuniquepatients"), ndose
  )
  nipde_by_dose <- weighted_dose_mean(results, c("nipde_by_dose", "nipde"), ndose)
  ntox_by_dose <- weighted_dose_mean(results, c("ntox"), ndose)
  estimated_pj <- weighted_estimated_p(results, ndose)

  total_admin <- unlist(lapply(results, function(res) {
    trial_vector(res, c("total_admin", "total_admin_by_trial", "totaln"), as.integer(res$ntrial))
  }), use.names = FALSE)
  total_unique <- unlist(lapply(results, function(res) {
    trial_vector(res, c("total_unique", "total_unique_by_trial"), as.integer(res$ntrial))
  }), use.names = FALSE)
  duration <- unlist(lapply(results, function(res) {
    trial_vector(res, c("duration", "duration_by_trial"), as.integer(res$ntrial))
  }), use.names = FALSE)

  early_stop_pct <- 100 * sum(final_mtd == 99L, na.rm = TRUE) / ntrial.total
  na_pct <- 100 * sum(
    is.na(final_mtd) |
      (!is.na(final_mtd) & final_mtd != 99L &
         (final_mtd < 1L | final_mtd > ndose))
  ) / ntrial.total

  common <- data.frame(
    Scenario_Set = scenario$Scenario_Set,
    Scenario_Name = scenario$Scenario_Name,
    Scenario = meta$Scenario[1L],
    Source_Scenario = scenario$Source_Scenario,
    True_MTD = scenario$True_MTD,
    Scenario_Attempt = scenario$Scenario_Attempt,
    Model = meta$Model[1L],
    Method = meta$Method[1L],
    CRM_r_model = if (identical(meta$Model[1L], "CRM")) meta$Method[1L] else NA_character_,
    Alpha_true = meta$Alpha_true[1L],
    r_carry = meta$r_carry[1L],
    Accrual = meta$Accrual[1L],
    T_assess = meta$T_assess[1L],
    Cohort_Size = meta$Cohort_Size[1L],
    Nmax_eff = meta$Nmax_eff[1L],
    Cycle_Max = meta$Cycle_Max[1L],
    Continuous_Enrollment = meta$Continuous_Enrollment[1L],
    Restrict_To_Tried = meta$Restrict_To_Tried[1L],
    Restrict_To_Target = meta$Restrict_To_Target[1L],
    n_jobs_found = nrow(meta),
    ntrial_from_files = ntrial.total,
    stringsAsFactors = FALSE
  )

  dose_summary <- common[rep(1L, ndose), , drop = FALSE]
  dose_summary$Dose <- seq_len(ndose)
  dose_summary$True_DLT_rate <- scenario$p.true
  dose_summary$True_IPDE_DLT_rate <- scenario$p.true_ipde
  dose_summary$Estimated_pj <- estimated_pj
  dose_summary$MTD_Selection_pct <- selection_pct
  dose_summary$Pts_Treated <- n_by_dose
  dose_summary$Unique_Pts_by_Dose <- unique_n_by_dose
  dose_summary$IPDE_Doses <- nipde_by_dose
  dose_summary$DLTs <- ntox_by_dose
  dose_summary$Total_Administrations <- mean(total_admin, na.rm = TRUE)
  dose_summary$Total_Unique_Patients <- mean(total_unique, na.rm = TRUE)
  dose_summary$Early_Stopping_pct <- early_stop_pct
  dose_summary$No_Selection_pct <- na_pct
  dose_summary$Duration <- mean(duration, na.rm = TRUE)
  dose_summary$Estimated_pj_available <- as.integer(any(is.finite(estimated_pj)))

  metrics <- c(
    "True DLT rate",
    "True IPDE DLT rate",
    "Estimated toxicity probability",
    "MTD Selection %",
    "# Pts Treated",
    "# Unique Pts by Dose",
    "# IPDE Doses",
    "# DLTs",
    "# Total Administrations",
    "# Total Unique Patients",
    "% Early Stopping",
    "% No Selection",
    "Duration"
  )

  table_summary <- common[rep(1L, length(metrics)), , drop = FALSE]
  table_summary$Metric <- metrics
  dose_cols <- paste0("D", seq_len(ndose))
  for (dose_col in dose_cols) table_summary[[dose_col]] <- NA_real_
  table_summary$Total <- NA_real_
  table_summary$Duration <- NA_real_

  table_summary[1L, dose_cols] <- scenario$p.true
  table_summary[2L, dose_cols] <- scenario$p.true_ipde
  table_summary[3L, dose_cols] <- estimated_pj
  table_summary[4L, dose_cols] <- selection_pct
  table_summary[5L, dose_cols] <- n_by_dose
  table_summary[6L, dose_cols] <- unique_n_by_dose
  table_summary[7L, dose_cols] <- nipde_by_dose
  table_summary[8L, dose_cols] <- ntox_by_dose
  table_summary[9L, "Total"] <- mean(total_admin, na.rm = TRUE)
  table_summary[10L, "Total"] <- mean(total_unique, na.rm = TRUE)
  table_summary[11L, "Total"] <- early_stop_pct
  table_summary[12L, "Total"] <- na_pct
  table_summary[13L, "Duration"] <- mean(duration, na.rm = TRUE)

  list(dose_summary = dose_summary, table_summary = table_summary)
}

round_numeric_df <- function(dat, digits = 4) {
  is_numeric <- vapply(dat, is.numeric, logical(1))
  dat[is_numeric] <- lapply(dat[is_numeric], round, digits = digits)
  dat
}

## -------------------------------
## Discover and summarize results
## -------------------------------

jobs.expected <- parse_job_spec(expected_jobs_spec)

if (!dir.exists(results_root)) {
  stop("TITE results root does not exist: ", results_root)
}
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

files <- list.files(
  results_root,
  pattern = "-job-[0-9]+-combined\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)
if (length(files) == 0L) {
  stop("No combined TITE result files found under: ", results_root)
}

metadata <- lapply(files, parse_combined_file)
keep <- !vapply(metadata, is.null, logical(1))
if (any(!keep)) {
  warning("Skipping ", sum(!keep), " combined file(s) whose names do not match the TITE runner convention.")
}
metadata <- do.call(rbind, metadata[keep])
if (nrow(metadata) == 0L) {
  stop("No combined TITE result files could be parsed under: ", results_root)
}

groups <- split(metadata, metadata$group_key)
dose_outputs <- vector("list", length(groups))
table_outputs <- vector("list", length(groups))
missing_outputs <- vector("list", length(groups))

for (i in seq_along(groups)) {
  meta <- groups[[i]]
  meta <- meta[order(meta$Job), , drop = FALSE]

  cat(
    "Summarizing scenario ", meta$Scenario[1L],
    ", model ", meta$Model[1L],
    ", method ", meta$Method[1L],
    ", ", nrow(meta), " job file(s).\n",
    sep = ""
  )

  results <- lapply(meta$file, readRDS)
  one <- summarize_tite_group(meta, results)
  dose_outputs[[i]] <- one$dose_summary
  table_outputs[[i]] <- one$table_summary

  missing_jobs <- if (length(jobs.expected) == 0L) integer(0) else {
    setdiff(jobs.expected, meta$Job)
  }
  missing_outputs[[i]] <- data.frame(
    Scenario_Set = meta$Scenario_Set[1L],
    Scenario = meta$Scenario[1L],
    Model = meta$Model[1L],
    Method = meta$Method[1L],
    Alpha_true = meta$Alpha_true[1L],
    r_carry = meta$r_carry[1L],
    Accrual = meta$Accrual[1L],
    Cycle_Max = meta$Cycle_Max[1L],
    Nmax_eff = meta$Nmax_eff[1L],
    Restrict_To_Tried = meta$Restrict_To_Tried[1L],
    Restrict_To_Target = meta$Restrict_To_Target[1L],
    n_jobs_found = nrow(meta),
    n_jobs_expected = length(jobs.expected),
    n_jobs_missing = length(missing_jobs),
    missing_jobs = collapse_integer_ranges(missing_jobs),
    stringsAsFactors = FALSE
  )
}

dose_summary <- do.call(rbind, dose_outputs)
table_summary <- do.call(rbind, table_outputs)
missing_log <- do.call(rbind, missing_outputs)

dose_summary_out <- round_numeric_df(dose_summary, digits = 4)
table_summary_out <- round_numeric_df(table_summary, digits = 4)

out.tag <- paste0("All_AIDE_TITE_OC_", scenario_set_name)
out.dose.csv <- file.path(out.dir, paste0(out.tag, "_dose_summary.csv"))
out.table.csv <- file.path(out.dir, paste0(out.tag, "_table_summary.csv"))
out.missing.csv <- file.path(out.dir, paste0(out.tag, "_missing_jobs.csv"))
out.dose.rds <- file.path(out.dir, paste0(out.tag, "_dose_summary.rds"))
out.table.rds <- file.path(out.dir, paste0(out.tag, "_table_summary.rds"))
out.missing.rds <- file.path(out.dir, paste0(out.tag, "_missing_jobs.rds"))

write.csv(dose_summary_out, out.dose.csv, row.names = FALSE)
write.csv(table_summary_out, out.table.csv, row.names = FALSE)
write.csv(missing_log, out.missing.csv, row.names = FALSE)
saveRDS(dose_summary, out.dose.rds)
saveRDS(table_summary, out.table.rds)
saveRDS(missing_log, out.missing.rds)

cat("\nSaved dose-level summary: ", out.dose.csv, "\n", sep = "")
cat("Saved table-style summary: ", out.table.csv, "\n", sep = "")
cat("Saved missing-job log: ", out.missing.csv, "\n", sep = "")
