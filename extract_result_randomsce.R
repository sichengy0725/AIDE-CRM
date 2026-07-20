## ============================================================
## Extract AIDE OC results for random scenarios
## Matched to run_oc_AIDE.R random scenario mode
##
## Output averages scenario-level metrics over the 10,000 random
## scenarios used by run_oc_AIDE.R:
##   1. average selection of each dose level
##   2. average allocation of each dose level
##   3. average number of total unique patients
##   4. average trial duration
##   5. average MTD selection
##   6. average MTD allocation
##   7. average overdose selection
##   8. average overdose allocation
##   9. average underdose allocation
##  10. MCSE of the average MTD selection rate
##  11. paired MTD-selection difference from the oracle (mean and SD)
##
## Selection metrics are percentages. Allocation metrics use n_by_dose,
## the AIDE dose-administration allocation stored by get_oc_sim_AIDE().
## To use patient-level dose allocation instead, set:
##   allocation_field <- "unique_n_by_dose"
## ============================================================

rm(list = ls())

## -------------------------------
## Settings matched to run_oc_AIDE.R
## -------------------------------

## setwd("/rsrch8/home/biostatistics/syang10/AIDE")

target <- 0.30

scenario_dose_count <- 5L
random_nscenario <- 10000L
random_target <- target
random_target_diff_below <- 0.05
random_target_diff_above <- 0.05
scenario_dir <- "scenario_sets"

if (!isTRUE(all.equal(random_target_diff_below, random_target_diff_above))) {
  stop("random_target_diff_below and random_target_diff_above must be equal.")
}
target_gap <- random_target_diff_below
if (!target_gap %in% c(0.05, 0.10, 0.15)) {
  stop("target_gap must be one of: 0.05, 0.10, 0.15.")
}

## run_oc_AIDE.R defaults
T_assess <- 28
C <- 3L
cycle_max_list <- c(1L, 2L, 3L)
Nmax_eff_list <- c(30L)
dose_cap <- 3L
continuous_enrollment <- TRUE
restrict_to_tried_list <- c(TRUE)
restrict_to_target_list <- c(FALSE)

## CFO results are available only for alpha_true = 0.6 and 0.9; retain the
## full alpha grid for BOIN and CRM.
alpha_true_list_non_cfo <- c(0, 0.3, 0.6, 0.9)
alpha_true_list_cfo <- c(0.6, 0.9)
arrival_rate_list <- c(1 / 14)

model_list <- c("BOIN", "CRM", "CFO")

boin_method_list <- c("approx1", "approx2")
boin_r_estimator_list <- c("r_fixed")
boin_r_carry_list_r_mle <- c(0)
boin_r_carry_list_fixed <- c(0)

crm_r_model_list <- if (scenario_dose_count == 5L) {
  c("r_fixed", "random", "level", "alpha_crm", "cumu_crm")
} else {
  c("r_fixed", "random", "alpha_crm", "cumu_crm")
}

crm_r_carry_values <- list(
  r_fixed   = c(0),
  random    = c(0),
  level     = c(0),
  alpha_crm = c(0),
  cumu_crm  = c(0)
)

## The reference in the paired comparison is the no-carryover fixed CRM.
## It is taken from the alpha_true = 0 setting, regardless of the alpha_true
## value for the method being compared.
oracle_model <- "CRM"
oracle_crm_r_model <- "r_fixed"
oracle_alpha_true <- 0
oracle_r_carry <- 0

cfo_method_list <- c("empirical", "pride")

## Allocation metric requested here. This matches the existing extractor's
## "# Pts Treated" / n_by_dose metric, which is dose administrations.
allocation_field <- "n_by_dose"

## Progress print frequency while reading files for one setting.
progress_every_files <- 1000L

## -------------------------------
## Helpers copied to match runner names
## -------------------------------

fmt_num <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return("NA")
  gsub("\\.", "p", as.character(x))
}

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

fmt_gap <- function(x) {
  gsub("\\.", "p", formatC(x, format = "f", digits = 2))
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

load_scenario_file <- function(file) {
  if (!file.exists(file)) {
    stop("Scenario file does not exist: ", file)
  }

  scenario_meta <- utils::read.csv(
    file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

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
    scenario_meta$Scenario_Group <- "random_scenario"
  }
  if (!"Attempt" %in% names(scenario_meta)) {
    scenario_meta$Attempt <- NA_integer_
  }

  required <- c("Scenario", "Source_Scenario", "Scenario_Group", "True_MTD", "Attempt")
  missing <- setdiff(required, names(scenario_meta))
  if (length(missing) > 0L) {
    stop("Scenario file is missing required columns: ", paste(missing, collapse = ", "))
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

get_crm_r_loop <- function(crm_r_model) {
  out <- crm_r_carry_values[[crm_r_model]]
  if (is.null(out)) {
    stop("No r_carry values specified for CRM model: ", crm_r_model)
  }
  out
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
    ## The runner writes the fixed-r CRM as crm_r_fixed.  Accept the older
    ## aliases here as well so the extracted label is always unambiguous.
    if (crm_r_model %in% c("fixed", "r_fixed", "crm_fixed", "crm_r_fixed")) {
      return("crm_r_fixed")
    }
    return(paste0("crm_", crm_r_model))
  }
  paste0("cfo_", cfo_method)
}

make_foldername <- function(model,
                            method_tag,
                            cycle_max,
                            arrival_rate,
                            Nmax_eff,
                            restrict_to_tried,
                            restrict_to_target,
                            scenario_set) {
  paste0(
    scenario_set,
    "-model-", model,
    "-opt-", method_tag,
    "-w-", fmt_short(T_assess),
    "-c-", fmt_short(C),
    "-cyc-", fmt_short(cycle_max),
    "-rate-", fmt_short(arrival_rate),
    "-Nmax-", fmt_short(Nmax_eff),
    "-tried-", as.integer(isTRUE(restrict_to_tried)),
    if (isTRUE(restrict_to_target)) "-ptarget-1" else ""
  )
}

resolve_folderpath <- function(results_root, foldername) {
  nested <- file.path(results_root, foldername)
  if (dir.exists(nested) || !dir.exists(foldername)) {
    return(nested)
  }
  foldername
}

regex_escape <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x, perl = TRUE)
}

make_file_pattern <- function(scenario_set,
                              model,
                              method_tag,
                              alpha_true,
                              r_carry,
                              arrival_rate,
                              cycle_max,
                              Nmax_eff,
                              continuous_enrollment,
                              restrict_to_tried,
                              restrict_to_target) {
  paste0(
    "^",
    regex_escape(scenario_set),
    "_SC[0-9]+_",
    regex_escape(model),
    "_",
    regex_escape(method_tag),
    "_a",
    regex_escape(fmt_short(alpha_true)),
    "_r",
    regex_escape(fmt_short(r_carry)),
    "_rate",
    regex_escape(fmt_short(arrival_rate)),
    "_cyc",
    regex_escape(fmt_short(cycle_max)),
    "_Nmax",
    regex_escape(fmt_short(Nmax_eff)),
    "_cont",
    as.integer(isTRUE(continuous_enrollment)),
    "_tried",
    as.integer(isTRUE(restrict_to_tried)),
    if (isTRUE(restrict_to_target)) "_ptarget1" else "",
    "-job-[0-9]+-combined\\.rds$"
  )
}

extract_scenario_id_from_file <- function(files, scenario_set) {
  rx <- paste0("^", regex_escape(scenario_set), "_SC([0-9]+)_.*$")
  as.integer(sub(rx, "\\1", basename(files), perl = TRUE))
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

pad_numeric <- function(x, len, fill = 0) {
  out <- rep(fill, len)
  if (is.null(x) || length(x) == 0L) {
    return(out)
  }
  z <- as.numeric(x)
  L <- min(length(z), len)
  out[seq_len(L)] <- z[seq_len(L)]
  out
}

sum_count_finite <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(c(sum = 0, count = 0))
  }
  z <- as.numeric(x)
  keep <- is.finite(z)
  c(sum = sum(z[keep]), count = sum(keep))
}

calc_mtd_selection_mcse_pct <- function(mtd_selection_pct, ntrial) {
  ## The reported MTD-selection percentage is the equally weighted mean of
  ## scenario-level rates.  Conditional on the sampled scenarios, its Monte
  ## Carlo variance is sum_s p_s (1 - p_s) / n_s divided by S^2.
  keep <- is.finite(mtd_selection_pct) & is.finite(ntrial) & ntrial > 0
  if (!any(keep)) return(NA_real_)

  p_hat <- pmin(pmax(mtd_selection_pct[keep] / 100, 0), 1)
  n_s <- ntrial[keep]
  100 * sqrt(sum(p_hat * (1 - p_hat) / n_s) / length(p_hat)^2)
}

make_data_key <- function(dat, columns) {
  if (nrow(dat) == 0L) return(character(0))

  pieces <- lapply(columns, function(nm) {
    x <- dat[[nm]]
    if (is.numeric(x)) {
      return(formatC(x, digits = 15, format = "fg", flag = "#"))
    }
    x <- as.character(x)
    x[is.na(x)] <- "<NA>"
    x
  })
  do.call(paste, c(pieces, sep = "\r"))
}

collapse_integer_ranges <- function(x) {
  x <- sort(unique(as.integer(x)))
  x <- x[!is.na(x)]
  if (length(x) == 0L) return("")

  starts <- x[c(TRUE, diff(x) != 1L)]
  ends <- x[c(diff(x) != 1L, TRUE)]
  paste(
    ifelse(starts == ends, as.character(starts), paste0(starts, "-", ends)),
    collapse = ","
  )
}

new_accumulator <- function(scenario_id_list, ndose) {
  nscenario <- length(scenario_id_list)
  acc <- new.env(parent = emptyenv())
  acc$scenario_id_list <- scenario_id_list
  acc$scenario_pos <- stats::setNames(seq_along(scenario_id_list), as.character(scenario_id_list))
  acc$ndose <- ndose
  acc$file_count <- integer(nscenario)
  acc$ntrial <- numeric(nscenario)
  acc$selection_count <- matrix(0, nrow = nscenario, ncol = ndose)
  acc$allocation_sum <- matrix(0, nrow = nscenario, ncol = ndose)
  acc$total_unique_sum <- numeric(nscenario)
  acc$total_unique_count <- numeric(nscenario)
  acc$duration_sum <- numeric(nscenario)
  acc$duration_count <- numeric(nscenario)
  acc
}

get_result_ntrial <- function(res) {
  ntrial <- get_field(res, c("ntrial"), required = FALSE)
  if (!is.null(ntrial) && length(ntrial) > 0L && is.finite(as.numeric(ntrial[1]))) {
    return(as.numeric(ntrial[1]))
  }

  final_mtd <- get_field(res, c("final_MTD"), required = FALSE)
  if (!is.null(final_mtd)) {
    return(length(final_mtd))
  }

  1
}

get_selection_count <- function(res, ndose, ntrial) {
  sel_count <- get_field(res, c("sel_count"), required = FALSE)
  if (!is.null(sel_count)) {
    return(pad_numeric(sel_count, ndose, fill = 0))
  }

  selection_pct <- get_field(res, c("selection_pct"), required = FALSE)
  if (!is.null(selection_pct)) {
    return(pad_numeric(selection_pct, ndose, fill = 0) * ntrial / 100)
  }

  final_mtd <- get_field(res, c("final_MTD"), required = FALSE)
  if (!is.null(final_mtd)) {
    final_mtd <- as.integer(final_mtd)
    return(tabulate(final_mtd[!is.na(final_mtd) & final_mtd >= 1L & final_mtd <= ndose], nbins = ndose))
  }

  rep(0, ndose)
}

accumulate_result <- function(acc, res, scenario_id, allocation_field) {
  pos <- acc$scenario_pos[[as.character(scenario_id)]]
  if (is.null(pos) || is.na(pos)) {
    warning("Skipping result for scenario not in scenario file: ", scenario_id)
    return(invisible(FALSE))
  }

  ndose <- acc$ndose
  ntrial <- get_result_ntrial(res)

  acc$file_count[pos] <- acc$file_count[pos] + 1L
  acc$ntrial[pos] <- acc$ntrial[pos] + ntrial
  acc$selection_count[pos, ] <- acc$selection_count[pos, ] +
    get_selection_count(res, ndose = ndose, ntrial = ntrial)

  allocation <- get_field(res, c(allocation_field), required = TRUE)
  acc$allocation_sum[pos, ] <- acc$allocation_sum[pos, ] +
    pad_numeric(allocation, ndose, fill = 0) * ntrial

  total_unique <- get_field(res, c("total_unique"), required = FALSE)
  if (!is.null(total_unique)) {
    sc <- sum_count_finite(total_unique)
  } else {
    total_unique_mean <- get_field(res, c("total_unique_mean"), required = TRUE)
    sc <- c(sum = as.numeric(total_unique_mean[1]) * ntrial, count = ntrial)
  }
  acc$total_unique_sum[pos] <- acc$total_unique_sum[pos] + sc[["sum"]]
  acc$total_unique_count[pos] <- acc$total_unique_count[pos] + sc[["count"]]

  duration <- get_field(res, c("duration"), required = FALSE)
  if (!is.null(duration)) {
    sc <- sum_count_finite(duration)
  } else {
    duration_mean <- get_field(res, c("duration_mean"), required = TRUE)
    sc <- c(sum = as.numeric(duration_mean[1]) * ntrial, count = ntrial)
  }
  acc$duration_sum[pos] <- acc$duration_sum[pos] + sc[["sum"]]
  acc$duration_count[pos] <- acc$duration_count[pos] + sc[["count"]]

  invisible(TRUE)
}

make_setting_meta <- function(model,
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
                              restrict_to_tried,
                              restrict_to_target,
                              scenario_set_name) {
  crm_r_model_label <- if (
    model == "CRM" && crm_r_model %in% c("fixed", "r_fixed", "crm_fixed", "crm_r_fixed")
  ) {
    "r_fixed"
  } else {
    crm_r_model
  }

  data.frame(
    Scenario_Set = scenario_set_name,
    Model = model,
    Method = method_tag,
    BOIN_Method = if (model == "BOIN") boin_method else NA_character_,
    BOIN_r_estimator = if (model == "BOIN") boin_r_estimator else NA_character_,
    CRM_r_model = if (model == "CRM") crm_r_model_label else NA_character_,
    CFO_Method = if (model == "CFO") cfo_method else NA_character_,
    Alpha_true = alpha_true,
    r_carry = r_carry,
    Accrual = arrival_rate,
    T_assess = T_assess,
    Nmax_eff = Nmax_eff,
    Dose_Cap = dose_cap,
    Cycle_Max = cycle_max,
    Continuous_Enrollment = as.integer(isTRUE(continuous_enrollment)),
    Restrict_To_Tried = as.integer(isTRUE(restrict_to_tried)),
    Restrict_To_Target = as.integer(isTRUE(restrict_to_target)),
    stringsAsFactors = FALSE
  )
}

summarize_accumulator <- function(acc, scenario_meta, setting_meta, n_files_found) {
  ndose <- acc$ndose
  scenario_id_list <- acc$scenario_id_list
  scenario_rows <- match(scenario_id_list, scenario_meta$Scenario)
  true_mtd <- scenario_meta$True_MTD[scenario_rows]

  found <- acc$ntrial > 0
  found_idx <- which(found)
  n_found <- length(found_idx)
  n_expected <- length(scenario_id_list)
  missing_ids <- scenario_id_list[!found]

  avg_selection_by_dose <- rep(NA_real_, ndose)
  avg_allocation_by_dose <- rep(NA_real_, ndose)
  avg_total_unique <- NA_real_
  avg_duration <- NA_real_
  avg_mtd_selection <- NA_real_
  avg_mtd_allocation <- NA_real_
  avg_overdose_selection <- NA_real_
  avg_overdose_allocation <- NA_real_
  avg_underdose_allocation <- NA_real_
  mcse_mtd_selection <- NA_real_
  scenario_mtd_summary <- cbind(
    setting_meta[0, , drop = FALSE],
    data.frame(
      Scenario = integer(0),
      ntrial = numeric(0),
      MTD_Selection_pct = numeric(0),
      stringsAsFactors = FALSE
    )
  )

  if (n_found > 0L) {
    ntrial_found <- acc$ntrial[found_idx]

    selection_pct <- sweep(
      acc$selection_count[found_idx, , drop = FALSE],
      1,
      ntrial_found,
      "/"
    ) * 100

    allocation <- sweep(
      acc$allocation_sum[found_idx, , drop = FALSE],
      1,
      ntrial_found,
      "/"
    )

    total_unique_mean <- acc$total_unique_sum[found_idx] /
      ifelse(acc$total_unique_count[found_idx] > 0, acc$total_unique_count[found_idx], NA_real_)

    duration_mean <- acc$duration_sum[found_idx] /
      ifelse(acc$duration_count[found_idx] > 0, acc$duration_count[found_idx], NA_real_)

    avg_selection_by_dose <- colMeans(selection_pct, na.rm = TRUE)
    avg_allocation_by_dose <- colMeans(allocation, na.rm = TRUE)
    avg_total_unique <- mean(total_unique_mean, na.rm = TRUE)
    avg_duration <- mean(duration_mean, na.rm = TRUE)

    mtd_selection <- rep(NA_real_, n_found)
    mtd_allocation <- rep(NA_real_, n_found)
    overdose_selection <- rep(NA_real_, n_found)
    overdose_allocation <- rep(NA_real_, n_found)
    underdose_allocation <- rep(NA_real_, n_found)

    for (i in seq_along(found_idx)) {
      tm <- true_mtd[found_idx[i]]
      if (is.na(tm) || tm < 1L || tm > ndose) {
        next
      }

      mtd_selection[i] <- selection_pct[i, tm]
      mtd_allocation[i] <- allocation[i, tm]

      overdose_doses <- if (tm < ndose) seq.int(tm + 1L, ndose) else integer(0)
      underdose_doses <- if (tm > 1L) seq_len(tm - 1L) else integer(0)

      overdose_selection[i] <- if (length(overdose_doses) > 0L) {
        sum(selection_pct[i, overdose_doses], na.rm = TRUE)
      } else {
        0
      }

      overdose_allocation[i] <- if (length(overdose_doses) > 0L) {
        sum(allocation[i, overdose_doses], na.rm = TRUE)
      } else {
        0
      }

      underdose_allocation[i] <- if (length(underdose_doses) > 0L) {
        sum(allocation[i, underdose_doses], na.rm = TRUE)
      } else {
        0
      }
    }

    avg_mtd_selection <- mean(mtd_selection, na.rm = TRUE)
    avg_mtd_allocation <- mean(mtd_allocation, na.rm = TRUE)
    avg_overdose_selection <- mean(overdose_selection, na.rm = TRUE)
    avg_overdose_allocation <- mean(overdose_allocation, na.rm = TRUE)
    avg_underdose_allocation <- mean(underdose_allocation, na.rm = TRUE)
    mcse_mtd_selection <- calc_mtd_selection_mcse_pct(
      mtd_selection_pct = mtd_selection,
      ntrial = ntrial_found
    )

    scenario_mtd_summary <- cbind(
      setting_meta[rep(1L, n_found), , drop = FALSE],
      data.frame(
        Scenario = scenario_id_list[found_idx],
        ntrial = ntrial_found,
        MTD_Selection_pct = mtd_selection,
        stringsAsFactors = FALSE
      )
    )
  }

  wide_summary <- cbind(
    setting_meta,
    data.frame(
      n_scenarios_expected = n_expected,
      n_scenarios_found = n_found,
      n_scenarios_missing = length(missing_ids),
      n_files_found = n_files_found,
      ntrial_total_from_files = sum(acc$ntrial),
      avg_ntrial_per_found_scenario = if (n_found > 0L) mean(acc$ntrial[found_idx]) else NA_real_,
      stringsAsFactors = FALSE
    )
  )

  for (j in seq_len(ndose)) {
    wide_summary[[paste0("Avg_Selection_pct_D", j)]] <- avg_selection_by_dose[j]
  }
  for (j in seq_len(ndose)) {
    wide_summary[[paste0("Avg_Allocation_D", j)]] <- avg_allocation_by_dose[j]
  }

  wide_summary$Avg_Total_Unique_Patients <- avg_total_unique
  wide_summary$Avg_Trial_Duration <- avg_duration
  wide_summary$Avg_MTD_Selection_pct <- avg_mtd_selection
  wide_summary$MCSE_MTD_Selection_pct <- mcse_mtd_selection
  wide_summary$Avg_MTD_Allocation <- avg_mtd_allocation
  wide_summary$Avg_Overdose_Selection_pct <- avg_overdose_selection
  wide_summary$Avg_Overdose_Allocation <- avg_overdose_allocation
  wide_summary$Avg_Underdose_Allocation <- avg_underdose_allocation

  dose_cols <- paste0("D", seq_len(ndose))
  metrics <- c(
    "Average dose selection %",
    "Average dose allocation",
    "Average total unique patients",
    "Average trial duration",
    "Average MTD selection %",
    "MCSE of MTD selection %",
    "Average MTD allocation",
    "Average overdose selection %",
    "Average overdose allocation",
    "Average underdose allocation",
    "Mean paired MTD selection difference vs oracle %",
    "SD paired MTD selection difference vs oracle %"
  )

  table_summary <- setting_meta[rep(1L, length(metrics)), , drop = FALSE]
  rownames(table_summary) <- NULL
  table_summary$Metric <- metrics

  for (dc in dose_cols) table_summary[[dc]] <- NA_real_
  table_summary$Total <- NA_real_
  table_summary$Duration <- NA_real_
  table_summary$n_scenarios_expected <- n_expected
  table_summary$n_scenarios_found <- n_found
  table_summary$n_scenarios_missing <- length(missing_ids)
  table_summary$n_files_found <- n_files_found
  table_summary$ntrial_total_from_files <- sum(acc$ntrial)
  table_summary$Oracle_Reference <- NA_character_
  table_summary$n_scenarios_paired_with_oracle <- NA_integer_

  table_summary[1, dose_cols] <- avg_selection_by_dose
  table_summary[2, dose_cols] <- avg_allocation_by_dose
  table_summary[3, "Total"] <- avg_total_unique
  table_summary[4, "Duration"] <- avg_duration
  table_summary[5, "Total"] <- avg_mtd_selection
  table_summary[6, "Total"] <- mcse_mtd_selection
  table_summary[7, "Total"] <- avg_mtd_allocation
  table_summary[8, "Total"] <- avg_overdose_selection
  table_summary[9, "Total"] <- avg_overdose_allocation
  table_summary[10, "Total"] <- avg_underdose_allocation

  missing_summary <- cbind(
    setting_meta,
    data.frame(
      n_scenarios_expected = n_expected,
      n_scenarios_found = n_found,
      n_scenarios_missing = length(missing_ids),
      n_files_found = n_files_found,
      missing_scenarios = collapse_integer_ranges(missing_ids),
      stringsAsFactors = FALSE
    )
  )

  list(
    wide_summary = wide_summary,
    table_summary = table_summary,
    missing_summary = missing_summary,
    scenario_mtd_summary = scenario_mtd_summary
  )
}

round_numeric_df <- function(dat, digits = 4) {
  num_cols <- vapply(dat, is.numeric, logical(1))
  dat[num_cols] <- lapply(dat[num_cols], round, digits)
  dat
}

## -------------------------------
## Scenario metadata and outputs
## -------------------------------

scenario_file <- make_random_scenario_filename(
  ndose = scenario_dose_count,
  nscenario = random_nscenario,
  target = random_target,
  target_diff_below = random_target_diff_below,
  target_diff_above = random_target_diff_above,
  out_dir = scenario_dir
)

scenario_set_name <- tools::file_path_sans_ext(basename(scenario_file))
results_root <- paste0("oc_results_cluster_AIDE_", scenario_set_name)

scenario_meta <- load_scenario_file(scenario_file)
dose_col_names <- paste0("Dose", seq_len(scenario_dose_count))
if (!all(dose_col_names %in% names(scenario_meta))) {
  stop(
    "scenario_meta is missing required dose columns: ",
    paste(setdiff(dose_col_names, names(scenario_meta)), collapse = ", ")
  )
}

scenario_id_list <- scenario_meta$Scenario
ndose_expected <- length(dose_col_names)

out.dir <- paste0("OC_summary_randomsce_targetgap", fmt_gap(target_gap))
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

cat("Scenario set:", scenario_set_name, "\n")
cat("Scenario file:", scenario_file, "\n")
cat("Results root:", results_root, "\n")
cat("Number of scenarios to average:", length(scenario_id_list), "\n")
cat("Allocation field:", allocation_field, "\n")

## -------------------------------
## Main extraction
## -------------------------------

all.wide.summary <- list()
all.table.summary <- list()
all.missing.summary <- list()
all.scenario.mtd.summary <- list()

summary_idx <- 1L
missing_idx <- 1L

for (Nmax_eff in Nmax_eff_list) {
  for (restrict_to_tried in restrict_to_tried_list) {
    for (restrict_to_target in restrict_to_target_list) {
      for (model in model_list) {
        alpha_true_loop <- if (model == "CFO") {
          alpha_true_list_cfo
        } else {
          alpha_true_list_non_cfo
        }

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

        for (alpha_true in alpha_true_loop) {
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

                        foldername <- make_foldername(
                          model = model,
                          method_tag = method_tag,
                          cycle_max = cycle_max,
                          arrival_rate = arrival_rate,
                          Nmax_eff = Nmax_eff,
                          restrict_to_tried = restrict_to_tried,
                          restrict_to_target = restrict_to_target,
                          scenario_set = scenario_set_name
                        )

                        folderpath <- resolve_folderpath(results_root, foldername)

                        setting_meta <- make_setting_meta(
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
                          restrict_to_tried = restrict_to_tried,
                          restrict_to_target = restrict_to_target,
                          scenario_set_name = scenario_set_name
                        )

                        cat("\n====================================\n")
                        cat("Model:", model, "\n")
                        cat("Method tag:", method_tag, "\n")
                        cat("alpha_true:", alpha_true, "\n")
                        cat("Cycle max:", cycle_max, "\n")
                        cat("Nmax eff:", Nmax_eff, "\n")
                        cat("Folder:", folderpath, "\n")

                        if (!dir.exists(folderpath)) {
                          cat("Folder not found; skipping.\n")
                          acc <- new_accumulator(scenario_id_list, ndose_expected)
                          one <- summarize_accumulator(
                            acc = acc,
                            scenario_meta = scenario_meta,
                            setting_meta = setting_meta,
                            n_files_found = 0L
                          )
                          all.missing.summary[[missing_idx]] <- one$missing_summary
                          missing_idx <- missing_idx + 1L
                          next
                        }

                        pattern <- make_file_pattern(
                          scenario_set = scenario_set_name,
                          model = model,
                          method_tag = method_tag,
                          alpha_true = alpha_true,
                          r_carry = r_carry,
                          arrival_rate = arrival_rate,
                          cycle_max = cycle_max,
                          Nmax_eff = Nmax_eff,
                          continuous_enrollment = continuous_enrollment,
                          restrict_to_tried = restrict_to_tried,
                          restrict_to_target = restrict_to_target
                        )

                        files.use <- list.files(
                          folderpath,
                          pattern = pattern,
                          full.names = TRUE
                        )

                        cat("Found", length(files.use), "combined files.\n")

                        acc <- new_accumulator(scenario_id_list, ndose_expected)

                        if (length(files.use) > 0L) {
                          file_scenario_ids <- extract_scenario_id_from_file(
                            files.use,
                            scenario_set = scenario_set_name
                          )

                          t.read <- Sys.time()

                          for (i in seq_along(files.use)) {
                            res <- readRDS(files.use[i])
                            accumulate_result(
                              acc = acc,
                              res = res,
                              scenario_id = file_scenario_ids[i],
                              allocation_field = allocation_field
                            )

                            if (progress_every_files > 0L &&
                                i %% progress_every_files == 0L) {
                              cat("  read", i, "of", length(files.use), "files\n")
                            }
                          }

                          cat(
                            "Read/summarize time:",
                            round(difftime(Sys.time(), t.read, units = "secs"), 2),
                            "sec\n"
                          )
                        }

                        one <- summarize_accumulator(
                          acc = acc,
                          scenario_meta = scenario_meta,
                          setting_meta = setting_meta,
                          n_files_found = length(files.use)
                        )

                        all.missing.summary[[missing_idx]] <- one$missing_summary
                        missing_idx <- missing_idx + 1L

                        if (length(files.use) == 0L) {
                          next
                        }

                        all.wide.summary[[summary_idx]] <- one$wide_summary
                        all.table.summary[[summary_idx]] <- one$table_summary
                        all.scenario.mtd.summary[[summary_idx]] <- one$scenario_mtd_summary
                        summary_idx <- summary_idx + 1L

                        cat(
                          "Scenarios found:",
                          one$wide_summary$n_scenarios_found,
                          "of",
                          one$wide_summary$n_scenarios_expected,
                          "\n"
                        )
                        cat("====================================\n")
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

if (length(all.wide.summary) == 0L) {
  stop("No combined RDS files were found. Check results_root, folder names, and run settings.")
}

wide.summary.df <- do.call(rbind, all.wide.summary)
table.summary.df <- do.call(rbind, all.table.summary)
missing.summary.df <- do.call(rbind, all.missing.summary)
scenario.mtd.summary.df <- do.call(rbind, all.scenario.mtd.summary)

## ------------------------------------------------------------
## Pair every method with the alpha_true = 0 fixed-CRM oracle
## ------------------------------------------------------------

setting_key_cols <- c(
  "Scenario_Set", "Model", "Method", "BOIN_Method", "BOIN_r_estimator",
  "CRM_r_model", "CFO_Method", "Alpha_true", "r_carry", "Accrual",
  "T_assess", "Nmax_eff", "Dose_Cap", "Cycle_Max",
  "Continuous_Enrollment", "Restrict_To_Tried", "Restrict_To_Target"
)
oracle_match_key_cols <- c(
  "Scenario_Set", "Accrual", "T_assess", "Nmax_eff", "Dose_Cap",
  "Cycle_Max", "Continuous_Enrollment", "Restrict_To_Tried",
  "Restrict_To_Target"
)

scenario.mtd.summary.df$.setting_key <- make_data_key(
  scenario.mtd.summary.df,
  setting_key_cols
)
scenario.mtd.summary.df$.oracle_match_key <- make_data_key(
  scenario.mtd.summary.df,
  oracle_match_key_cols
)
wide.summary.df$.setting_key <- make_data_key(wide.summary.df, setting_key_cols)
table.summary.df$.setting_key <- make_data_key(table.summary.df, setting_key_cols)

is_oracle <-
  scenario.mtd.summary.df$Model == oracle_model &
  scenario.mtd.summary.df$CRM_r_model == oracle_crm_r_model &
  abs(scenario.mtd.summary.df$Alpha_true - oracle_alpha_true) < 1e-12 &
  abs(scenario.mtd.summary.df$r_carry - oracle_r_carry) < 1e-12

oracle_setting_map <- unique(scenario.mtd.summary.df[
  is_oracle,
  c(".setting_key", ".oracle_match_key"),
  drop = FALSE
])
oracle_key_count <- table(oracle_setting_map$.oracle_match_key)
ambiguous_oracle_keys <- names(oracle_key_count)[oracle_key_count > 1L]
if (length(ambiguous_oracle_keys) > 0L) {
  warning(
    "Multiple alpha_true = 0 fixed-CRM oracle settings were found for ",
    length(ambiguous_oracle_keys),
    " design setting(s); oracle differences are reported as NA for those settings."
  )
}

setting_keys <- unique(scenario.mtd.summary.df$.setting_key)
oracle_comparison <- lapply(setting_keys, function(setting_key) {
  method_data <- scenario.mtd.summary.df[
    scenario.mtd.summary.df$.setting_key == setting_key,
    c("Scenario", "MTD_Selection_pct", ".oracle_match_key"),
    drop = FALSE
  ]
  match_key <- method_data$.oracle_match_key[1L]
  candidate_oracles <- oracle_setting_map$.setting_key[
    oracle_setting_map$.oracle_match_key == match_key
  ]

  out <- data.frame(
    .setting_key = setting_key,
    Oracle_Reference = NA_character_,
    n_scenarios_paired_with_oracle = NA_integer_,
    Mean_Paired_MTD_Selection_Diff_vs_Oracle_pct = NA_real_,
    SD_Paired_MTD_Selection_Diff_vs_Oracle_pct = NA_real_,
    stringsAsFactors = FALSE
  )

  if (length(candidate_oracles) != 1L) {
    return(out)
  }

  oracle_data <- scenario.mtd.summary.df[
    scenario.mtd.summary.df$.setting_key == candidate_oracles,
    c("Scenario", "MTD_Selection_pct"),
    drop = FALSE
  ]
  names(oracle_data)[2L] <- "Oracle_MTD_Selection_pct"
  names(method_data)[2L] <- "Method_MTD_Selection_pct"

  paired <- merge(method_data, oracle_data, by = "Scenario", all = FALSE)
  paired <- paired[
    is.finite(paired$Method_MTD_Selection_pct) &
      is.finite(paired$Oracle_MTD_Selection_pct),
    , drop = FALSE
  ]

  if (nrow(paired) == 0L) {
    return(out)
  }

  ## Positive values mean the method has higher MTD selection than the
  ## alpha_true = 0 fixed-CRM oracle for the same random scenario.
  delta <- paired$Method_MTD_Selection_pct - paired$Oracle_MTD_Selection_pct
  out$Oracle_Reference <- "CRM r_fixed (alpha_true = 0, r_carry = 0)"
  out$n_scenarios_paired_with_oracle <- nrow(paired)
  out$Mean_Paired_MTD_Selection_Diff_vs_Oracle_pct <- mean(delta)
  out$SD_Paired_MTD_Selection_Diff_vs_Oracle_pct <- if (length(delta) > 1L) {
    stats::sd(delta)
  } else {
    NA_real_
  }
  out
})
oracle_comparison.df <- do.call(rbind, oracle_comparison)

comparison_index <- match(
  wide.summary.df$.setting_key,
  oracle_comparison.df$.setting_key
)
wide.summary.df$Oracle_Reference <- oracle_comparison.df$Oracle_Reference[comparison_index]
wide.summary.df$n_scenarios_paired_with_oracle <-
  oracle_comparison.df$n_scenarios_paired_with_oracle[comparison_index]
wide.summary.df$Mean_Paired_MTD_Selection_Diff_vs_Oracle_pct <-
  oracle_comparison.df$Mean_Paired_MTD_Selection_Diff_vs_Oracle_pct[comparison_index]
wide.summary.df$SD_Paired_MTD_Selection_Diff_vs_Oracle_pct <-
  oracle_comparison.df$SD_Paired_MTD_Selection_Diff_vs_Oracle_pct[comparison_index]

comparison_index <- match(
  table.summary.df$.setting_key,
  oracle_comparison.df$.setting_key
)
table.summary.df$Oracle_Reference <- oracle_comparison.df$Oracle_Reference[comparison_index]
table.summary.df$n_scenarios_paired_with_oracle <-
  oracle_comparison.df$n_scenarios_paired_with_oracle[comparison_index]

mean_diff_rows <- table.summary.df$Metric ==
  "Mean paired MTD selection difference vs oracle %"
sd_diff_rows <- table.summary.df$Metric ==
  "SD paired MTD selection difference vs oracle %"
table.summary.df$Total[mean_diff_rows] <-
  oracle_comparison.df$Mean_Paired_MTD_Selection_Diff_vs_Oracle_pct[
    comparison_index[mean_diff_rows]
  ]
table.summary.df$Total[sd_diff_rows] <-
  oracle_comparison.df$SD_Paired_MTD_Selection_Diff_vs_Oracle_pct[
    comparison_index[sd_diff_rows]
  ]

wide.summary.df$.setting_key <- NULL
table.summary.df$.setting_key <- NULL

wide.summary.out <- round_numeric_df(wide.summary.df, digits = 4)
table.summary.out <- round_numeric_df(table.summary.df, digits = 4)

out.tag <- paste0("randomsce_targetgap", fmt_gap(target_gap))

out.wide.csv <- file.path(out.dir, paste0(out.tag, "_wide_summary.csv"))
out.table.csv <- file.path(out.dir, paste0(out.tag, "_table_summary.csv"))
out.missing.csv <- file.path(out.dir, paste0(out.tag, "_missing_scenarios.csv"))

out.wide.rds <- file.path(out.dir, paste0(out.tag, "_wide_summary.rds"))
out.table.rds <- file.path(out.dir, paste0(out.tag, "_table_summary.rds"))
out.missing.rds <- file.path(out.dir, paste0(out.tag, "_missing_scenarios.rds"))

write.csv(wide.summary.out, out.wide.csv, row.names = FALSE)
write.csv(table.summary.out, out.table.csv, row.names = FALSE)
write.csv(missing.summary.df, out.missing.csv, row.names = FALSE)

saveRDS(wide.summary.df, out.wide.rds)
saveRDS(table.summary.df, out.table.rds)
saveRDS(missing.summary.df, out.missing.rds)

cat("\nSaved random-scenario wide summary:", out.wide.csv, "\n")
cat("Saved random-scenario table summary:", out.table.csv, "\n")
cat("Saved random-scenario missing-scenario log:", out.missing.csv, "\n")
