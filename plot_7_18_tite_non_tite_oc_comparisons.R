## ============================================================
## Presentation 7-18-2026: TITE versus non-TITE CRM OC plots
##
## The source summaries contain r_fixed and random CRM only.
## Produces four baseline-scenario figures (alpha = 0, 0.3, 0.6, 0.9)
## using the eight OC metrics from the Presentation 7-6 layout.
## ============================================================

target <- 0.30
ncycle_use <- 2L
nmax_use <- 30L
restrict_to_target_use <- 0L
accrual_rate_use <- 0.0179
accrual_days_label <- 56L
alpha_values <- c(0, 0.3, 0.6, 0.9)

input_dir <- "Presentation 7-18-2026"
non_tite_file <- file.path(
  input_dir,
  "All_AIDE_OC_Set_5dose_adaptive_r_37_modelsCRM_crmcrm_r_fixed_random_target0p3_w28_c3_cyc2_rate0p02_Nmax30_dose_summary.csv"
)
tite_file <- file.path(
  input_dir,
  "All_AIDE_TITE_OC_Set_5dose_adaptive_r_37_modelsCRM_crmr_fixed_random_target0p3_w28_c3_cyc2_rate0p02_Nmax30_ipde1_newfirst1_neval3_dcap100_dosecap3_dose_summary.csv"
)
out_dir <- file.path(input_dir, "Plots", "OC Ncycle2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_cols <- c(
  "Scenario", "True_MTD", "Model", "CRM_r_model", "Alpha_true",
  "Accrual", "Nmax_eff", "Cycle_Max", "Dose", "True_DLT_rate",
  "MTD_Selection_pct", "Pts_Treated", "IPDE_Doses",
  "Total_Unique_Patients", "Early_Stopping_pct", "Duration"
)

read_summary <- function(filename, source_name) {
  if (!file.exists(filename)) {
    stop("Cannot find ", source_name, " file: ", filename)
  }
  data <- read.csv(filename, stringsAsFactors = FALSE, check.names = FALSE)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop(source_name, " is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!"Restrict_To_Target" %in% names(data)) {
    data$Restrict_To_Target <- restrict_to_target_use
  }
  data
}

non_tite <- read_summary(non_tite_file, "Non-TITE")
tite <- read_summary(tite_file, "TITE")

scenario_set <- list(
  label = "Baseline Scenarios",
  scenarios = 1:6
)

design_specs <- data.frame(
  key = c("r_fixed_non_tite", "r_fixed_tite", "random_non_tite", "random_tite"),
  label = c("r_fixed CRM", "r_fixed CRM TITE", "random CRM", "random CRM TITE"),
  color = c("#1f77b4", "#1f77b4", "#e377c2", "#e377c2"),
  pch = c(16, 1, 18, 5),
  lty = c(1, 2, 1, 2),
  lwd = rep(2, 4),
  stringsAsFactors = FALSE
)

first_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA_real_ else x[[1]]
}

alpha_tag <- function(alpha_value) {
  value <- format(alpha_value, trim = TRUE, scientific = FALSE)
  paste0("alpha", gsub("\\.", "p", value))
}

filter_rows <- function(data, alpha_value) {
  data[
    data$Model == "CRM" &
      data$Scenario %in% scenario_set$scenarios &
      data$Nmax_eff == nmax_use &
      data$Cycle_Max == ncycle_use &
      data$Restrict_To_Target == restrict_to_target_use &
      abs(data$Alpha_true - alpha_value) < 1e-12 &
      abs(data$Accrual - accrual_rate_use) < 1e-10,
    ,
    drop = FALSE
  ]
}

add_design <- function(data, key) {
  data$Design <- key
  data
}

select_plot_rows <- function(alpha_value) {
  non_tite_current <- filter_rows(non_tite, alpha_value)
  tite_current <- filter_rows(tite, alpha_value)
  pieces <- list(
    add_design(
      non_tite_current[non_tite_current$CRM_r_model == "crm_r_fixed", , drop = FALSE],
      "r_fixed_non_tite"
    ),
    add_design(
      tite_current[tite_current$CRM_r_model == "r_fixed", , drop = FALSE],
      "r_fixed_tite"
    ),
    add_design(
      non_tite_current[non_tite_current$CRM_r_model == "random", , drop = FALSE],
      "random_non_tite"
    ),
    add_design(
      tite_current[tite_current$CRM_r_model == "random", , drop = FALSE],
      "random_tite"
    )
  )
  common_cols <- Reduce(intersect, lapply(pieces, names))
  do.call(rbind, lapply(pieces, function(data) data[, common_cols, drop = FALSE]))
}

validate_rows <- function(rows) {
  observed <- table(
    factor(rows$Design, levels = design_specs$key),
    factor(rows$Scenario, levels = scenario_set$scenarios)
  )
  if (any(observed != 5L)) {
    bad <- which(observed != 5L, arr.ind = TRUE)
    details <- apply(
      bad,
      1,
      function(index) {
        paste0(
          rownames(observed)[index[[1]]], "/scenario ", colnames(observed)[index[[2]]],
          " (", observed[index[[1]], index[[2]]], " rows)"
        )
      }
    )
    stop(
      "Each design/scenario combination must have five dose rows. Problems: ",
      paste(details, collapse = "; ")
    )
  }
  invisible(rows)
}

summarize_group <- function(data) {
  data <- data[order(data$Dose), , drop = FALSE]
  true_rates <- data$True_DLT_rate
  total_treated <- sum(data$Pts_Treated, na.rm = TRUE)
  if (!is.finite(total_treated) || total_treated <= 0L) {
    total_treated <- NA_real_
  }

  true_mtd <- first_value(data$True_MTD)
  selected_toxic <- sum(data$MTD_Selection_pct[true_rates > target + 1e-12], na.rm = TRUE)
  treated_toxic <- sum(data$Pts_Treated[true_rates > target + 1e-12], na.rm = TRUE)
  treated_subtherapeutic <- sum(data$Pts_Treated[true_rates < target - 1e-12], na.rm = TRUE)

  if (is.na(true_mtd) || true_mtd < 1L || true_mtd > nrow(data)) {
    correct_selection <- first_value(data$Early_Stopping_pct)
    treated_mtd <- 0
  } else {
    correct_selection <- data$MTD_Selection_pct[data$Dose == true_mtd][1]
    treated_mtd <- data$Pts_Treated[data$Dose == true_mtd][1]
  }

  data.frame(
    Scenario = first_value(data$Scenario),
    Design = first_value(data$Design),
    Correct_MTD_Selection = correct_selection,
    Overdose_Selection = selected_toxic,
    Patients_Allocated_MTD = 100 * treated_mtd / total_treated,
    Patients_Overdosed = 100 * treated_toxic / total_treated,
    Patients_Underdosed = 100 * treated_subtherapeutic / total_treated,
    Total_IPDE = sum(data$IPDE_Doses, na.rm = TRUE),
    Sample_Size = first_value(data$Total_Unique_Patients),
    Duration = first_value(data$Duration),
    stringsAsFactors = FALSE
  )
}

summarize_plot_data <- function(rows) {
  validate_rows(rows)
  groups <- split(rows, paste(rows$Design, rows$Scenario, sep = "__"))
  output <- do.call(rbind, lapply(groups, summarize_group))
  output$Scenario <- factor(output$Scenario, levels = scenario_set$scenarios)
  output$Scenario_Index <- as.integer(output$Scenario)
  output$Design <- factor(output$Design, levels = design_specs$key)
  output <- output[order(output$Design, output$Scenario_Index), , drop = FALSE]
  row.names(output) <- NULL
  output
}

metric_limits <- function(plot_data, metric) {
  values <- plot_data[[metric]]
  values <- values[is.finite(values)]
  if (length(values) == 0L) return(c(0, 1))
  if (metric %in% c("Correct_MTD_Selection", "Patients_Overdosed", "Patients_Underdosed")) {
    return(c(0, 100))
  }
  if (metric == "Overdose_Selection") {
    return(c(0, max(70, ceiling(max(values) / 10) * 10)))
  }
  if (metric == "Patients_Allocated_MTD") {
    return(c(0, max(80, ceiling(max(values) / 10) * 10)))
  }
  if (metric == "Total_IPDE") {
    return(c(0, max(10, ceiling(max(values) + 1))))
  }
  padding <- max(1, diff(range(values)) * 0.08)
  c(max(0, floor(min(values) - padding)), ceiling(max(values) + padding))
}

plot_metric <- function(plot_data, metric, panel_title, ylab, ylim) {
  xvals <- seq_along(scenario_set$scenarios)
  plot(
    range(xvals), ylim,
    type = "n", axes = FALSE,
    xlab = "Scenario", ylab = ylab, main = panel_title,
    font.main = 2, cex.main = 1.12, cex.lab = 1.02
  )
  axis(1, at = xvals, labels = scenario_set$scenarios)
  axis(2)
  box()
  abline(v = xvals, col = "#d9d9d9", lwd = 0.8)
  abline(h = pretty(ylim), col = "#e6e6e6", lwd = 0.8)
  for (ii in seq_len(nrow(design_specs))) {
    spec <- design_specs[ii, , drop = FALSE]
    data <- plot_data[plot_data$Design == spec$key, , drop = FALSE]
    data <- data[order(data$Scenario_Index), , drop = FALSE]
    keep <- is.finite(data[[metric]])
    if (sum(keep) == 0L) next
    lines(
      data$Scenario_Index[keep], data[[metric]][keep],
      type = "b", col = spec$color, pch = spec$pch,
      lty = spec$lty, lwd = spec$lwd, cex = 1
    )
  }
}

draw_plot <- function(plot_data, alpha_value) {
  metrics <- list(
    list("Correct_MTD_Selection", "A. Percentage of correct MTD selection", "Percentage"),
    list("Overdose_Selection", "B. Percentage of overdosing selection", "Percentage"),
    list("Patients_Allocated_MTD", "C. Percentage of patients allocated to the MTD", "Percentage"),
    list("Patients_Overdosed", "D. Percentage of patients overdosed", "Percentage"),
    list("Patients_Underdosed", "E. Percentage of patients underdosed", "Percentage"),
    list("Total_IPDE", "F. Average total number of IPDE", "Average"),
    list("Sample_Size", "G. Average sample size", "Average"),
    list("Duration", "H. Average trial duration", "Days")
  )

  layout(
    rbind(c(1, 2), c(3, 4), c(5, 6), c(7, 8), c(9, 9), c(10, 10)),
    heights = c(1, 1, 1, 1, 0.34, 0.24)
  )
  par(oma = c(0, 0, 3.0, 0), mar = c(3.4, 4.2, 2.1, 1.4), mgp = c(2.2, 0.7, 0))
  for (metric in metrics) {
    plot_metric(
      plot_data, metric[[1]], metric[[2]], metric[[3]],
      metric_limits(plot_data, metric[[1]])
    )
  }

  par(mar = c(0, 0, 0, 0), xpd = NA)
  plot.new()
  legend(
    "center", legend = design_specs$label,
    col = design_specs$color, pch = design_specs$pch,
    lty = design_specs$lty, lwd = design_specs$lwd,
    ncol = nrow(design_specs), bty = "n", cex = 0.96,
    x.intersp = 0.80, y.intersp = 1.15
  )

  plot.new()
  text(
    0.5, 0.72,
    paste0("Scenarios shown: ", paste(scenario_set$scenarios, collapse = ", "), ". True MTD labels come from the dose-summary CSV."),
    cex = 0.82
  )
  text(
    0.5, 0.42,
    "TITE and non-TITE summaries use accrual rate 0.0179 (56-day arrival setting), ncycle = 2, and Nmax = 30.",
    cex = 0.82
  )
  text(
    0.5, 0.14,
    "Patient-overdose percentages use treated counts above target divided by total treated counts.",
    cex = 0.82
  )

  mtext(
    paste0(
      "Baseline Scenarios: TITE vs non-TITE CRM (alpha = ",
      format(alpha_value, trim = TRUE, scientific = FALSE),
      ", ncycle = ", ncycle_use, ", Nmax = ", nmax_use,
      ", restrict_to_target = ", restrict_to_target_use, ")"
    ),
    outer = TRUE, side = 3, line = 0.8, cex = 1.32, font = 2
  )
}

out_files <- character(0)
for (alpha_value in alpha_values) {
  plot_data <- summarize_plot_data(select_plot_rows(alpha_value))
  output_file <- file.path(
    out_dir,
    paste0(
      "Set_5dose_baseline_TITE_vs_nonTITE_OC_", alpha_tag(alpha_value),
      "_ncycle2_Nmax30_rate0p0179.pdf"
    )
  )
  pdf(output_file, width = 13.5, height = 10.5, onefile = FALSE, useDingbats = FALSE)
  draw_plot(plot_data, alpha_value)
  dev.off()
  out_files <- c(out_files, output_file)
}

cat("Wrote ", length(out_files), " PDF files:\n", sep = "")
cat(paste0("  ", out_files), sep = "\n")
cat("\n")
