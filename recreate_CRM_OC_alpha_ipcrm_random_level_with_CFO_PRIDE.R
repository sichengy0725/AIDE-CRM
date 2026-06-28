target <- 0.30
alpha_values <- c(0, 0.3, 0.6, 0.9)
cycle_values <- 1:3

crm_file <- file.path(
  "Results",
  "All_AIDE_OC_Set_5dose_methods_prior_modelsCRM_crmfixed_level_random_alpha_crm_cumu_crm_target0p3_w28_c3_cyc1_2_3_rate0p07_Nmax30_dosecap3_cont1_tried1_jobs1to2000_dose_summary.csv"
)

cfo_file <- file.path(
  "Results",
  "All_AIDE_OC_Set_5dose_methods_prior_modelsBOIN_CFO_boinapprox1_approx2_restr_fixed_cfoempirical_pride_target0p3_w28_c3_cyc1_2_3_rate0p07_Nmax30_dosecap3_cont1_tried1_dose_summary.csv"
)

out_dir <- "Results"

crm <- read.csv(crm_file, stringsAsFactors = FALSE)
cfo <- read.csv(cfo_file, stringsAsFactors = FALSE)

required_cols <- c(
  "Scenario", "Method", "CRM_r_model", "CFO_Method", "Model", "Alpha_true",
  "Cycle_Max", "Dose", "True_DLT_rate", "MTD_Selection_pct", "Pts_Treated",
  "IPDE_Doses", "Total_Unique_Patients", "Early_Stopping_pct", "Duration"
)

missing_crm <- setdiff(required_cols, names(crm))
missing_cfo <- setdiff(required_cols, names(cfo))
if (length(missing_crm) > 0) {
  stop("CRM summary is missing columns: ", paste(missing_crm, collapse = ", "))
}
if (length(missing_cfo) > 0) {
  stop("BOIN/CFO summary is missing columns: ", paste(missing_cfo, collapse = ", "))
}

design_specs <- data.frame(
  key = c(
    "crm_fixed", "crm_alpha", "cfo_pride", "crm_ip", "crm_random",
    "oracle"
  ),
  label = c(
    "AIDE-CRM fixed, r = 0", "AIDE-alpha-CRM", "CFO PRIDE",
    "AIDE-IP-CRM", "AIDE-CRM random", "Oracle CRM"
  ),
  color = c("#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e", "#e377c2", "#17becf"),
  pch = c(16, 17, 8, 15, 18, 4),
  lty = c(1, 1, 4, 1, 1, 3),
  lwd = c(2, 2, 2, 2, 2, 2),
  stringsAsFactors = FALSE
)

alpha_tag <- function(alpha_value) {
  tag <- format(alpha_value, trim = TRUE, scientific = FALSE)
  paste0("alpha", gsub("\\.", "p", tag))
}

plot_scenario_label <- function(scenario) {
  ifelse(scenario <= 5, 6 - scenario, scenario)
}

first_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  x[[1]]
}

label_crm_design <- function(r_model) {
  out <- rep(NA_character_, length(r_model))
  out[r_model == "fixed"] <- "crm_fixed"
  out[r_model == "alpha_crm"] <- "crm_alpha"
  out[r_model == "cumu_crm"] <- "crm_ip"
  out[r_model == "random"] <- "crm_random"
  out
}

make_plot_rows <- function(alpha_value, ncycle) {
  crm_use <- subset(
    crm,
    abs(Alpha_true - alpha_value) < 1e-12 &
      Cycle_Max == ncycle &
      CRM_r_model %in% c("fixed", "alpha_crm", "cumu_crm", "random")
  )
  crm_use$Design <- label_crm_design(crm_use$CRM_r_model)
  crm_use$Is_Oracle <- FALSE

  cfo_pride <- subset(
    cfo,
    Model == "CFO" &
      Method == "cfo_pride" &
      CFO_Method == "pride" &
      abs(Alpha_true - alpha_value) < 1e-12 &
      Cycle_Max == ncycle
  )
  cfo_pride$Design <- "cfo_pride"
  cfo_pride$Is_Oracle <- FALSE

  oracle <- subset(
    crm,
    abs(Alpha_true - alpha_value) < 1e-12 &
      Cycle_Max == 1 &
      CRM_r_model == "fixed"
  )
  oracle$Design <- "oracle"
  oracle$Is_Oracle <- TRUE

  common_cols <- Reduce(intersect, list(names(crm_use), names(cfo_pride), names(oracle)))
  combined <- rbind(crm_use[common_cols], cfo_pride[common_cols], oracle[common_cols])
  combined <- combined[!is.na(combined$Design), ]
  combined
}

summarize_group <- function(d) {
  scenario <- first_value(d$Scenario)
  design <- first_value(d$Design)
  is_oracle <- isTRUE(first_value(d$Is_Oracle))

  dose_order <- order(d$Dose)
  d <- d[dose_order, ]

  true_rates <- d$True_DLT_rate
  safe_doses <- which(true_rates <= target + 1e-12)
  mtd_idx <- if (length(safe_doses) > 0) safe_doses[which.max(true_rates[safe_doses])] else NA_integer_

  total_treated <- sum(d$Pts_Treated, na.rm = TRUE)
  if (!is.finite(total_treated) || total_treated <= 0) {
    total_treated <- NA_real_
  }

  selected_toxic <- sum(d$MTD_Selection_pct[true_rates > target + 1e-12], na.rm = TRUE)
  treated_toxic <- sum(d$Pts_Treated[true_rates > target + 1e-12], na.rm = TRUE)
  treated_subtherapeutic <- sum(d$Pts_Treated[true_rates < target - 1e-12], na.rm = TRUE)

  if (is.na(mtd_idx)) {
    correct_selection <- first_value(d$Early_Stopping_pct)
    treated_mtd <- 0
  } else {
    correct_selection <- d$MTD_Selection_pct[mtd_idx]
    treated_mtd <- d$Pts_Treated[mtd_idx]
  }

  sample_size <- if (is_oracle) total_treated else first_value(d$Total_Unique_Patients)
  duration <- if (is_oracle) NA_real_ else first_value(d$Duration)

  data.frame(
    Scenario = scenario,
    Design = design,
    Correct_MTD_Selection = correct_selection,
    Overdose_Selection = selected_toxic,
    Patients_Allocated_MTD = 100 * treated_mtd / total_treated,
    Patients_Overdosed = 100 * treated_toxic / total_treated,
    Patients_Underdosed = 100 * treated_subtherapeutic / total_treated,
    Total_IPDE = sum(d$IPDE_Doses, na.rm = TRUE),
    Sample_Size = sample_size,
    Duration = duration,
    stringsAsFactors = FALSE
  )
}

summarize_plot_data <- function(rows) {
  split_rows <- split(rows, paste(rows$Design, rows$Scenario, sep = "__"))
  out <- do.call(rbind, lapply(split_rows, summarize_group))
  out$Original_Scenario <- out$Scenario
  out$Scenario <- plot_scenario_label(out$Original_Scenario)
  out$Design <- factor(out$Design, levels = design_specs$key)
  out <- out[order(out$Design, out$Scenario), ]
  row.names(out) <- NULL
  out
}

metric_limits <- function(plot_data, metric) {
  if (metric %in% c(
    "Correct_MTD_Selection", "Patients_Overdosed", "Patients_Underdosed"
  )) {
    return(c(0, 100))
  }
  if (metric == "Overdose_Selection") {
    return(c(0, max(70, ceiling(max(plot_data[[metric]], na.rm = TRUE) / 10) * 10)))
  }
  if (metric == "Patients_Allocated_MTD") {
    return(c(0, max(80, ceiling(max(plot_data[[metric]], na.rm = TRUE) / 10) * 10)))
  }
  if (metric == "Total_IPDE") {
    max_val <- max(plot_data[[metric]], na.rm = TRUE)
    return(c(0, max(10, ceiling(max_val + 1))))
  }

  vals <- plot_data[[metric]]
  vals <- vals[is.finite(vals)]
  pad <- max(1, diff(range(vals)) * 0.08)
  c(max(0, floor(min(vals) - pad)), ceiling(max(vals) + pad))
}

plot_metric <- function(plot_data, metric, panel_title, ylab, ylim) {
  scenarios <- sort(unique(plot_data$Scenario))
  plot(
    range(scenarios),
    ylim,
    type = "n",
    axes = FALSE,
    xlab = "Scenario",
    ylab = ylab,
    main = panel_title,
    font.main = 2,
    cex.main = 1.15,
    cex.lab = 1.05
  )
  axis(1, at = scenarios, labels = scenarios)
  axis(2)
  box()
  abline(v = scenarios, col = "#d9d9d9", lwd = 0.8)
  abline(h = pretty(ylim), col = "#e6e6e6", lwd = 0.8)

  for (ii in seq_len(nrow(design_specs))) {
    spec <- design_specs[ii, ]
    d <- plot_data[plot_data$Design == spec$key, ]
    d <- d[order(d$Scenario), ]
    y <- d[[metric]]
    x <- d$Scenario
    keep <- is.finite(y)
    if (sum(keep) == 0) {
      next
    }
    lines(
      x[keep],
      y[keep],
      type = "b",
      col = spec$color,
      pch = spec$pch,
      lty = spec$lty,
      lwd = spec$lwd,
      cex = 1
    )
  }
}

plot_one_cycle <- function(alpha_value, ncycle) {
  plot_data <- summarize_plot_data(make_plot_rows(alpha_value, ncycle))

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

  out_file <- file.path(
    out_dir,
    paste0("CRM_OC_alpha_ipcrm_random_cfo_pride_", alpha_tag(alpha_value), "_ncycle", ncycle, ".pdf")
  )
  pdf(out_file, width = 13.5, height = 10.5, onefile = FALSE, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)

  layout(
    rbind(c(1, 2), c(3, 4), c(5, 6), c(7, 8), c(9, 9), c(10, 10)),
    heights = c(1, 1, 1, 1, 0.34, 0.2)
  )
  par(oma = c(0, 0, 2.8, 0), mar = c(3.4, 4.2, 2.1, 1.4), mgp = c(2.2, 0.7, 0))

  for (metric in metrics) {
    metric_name <- metric[[1]]
    plot_metric(
      plot_data,
      metric_name,
      metric[[2]],
      metric[[3]],
      metric_limits(plot_data, metric_name)
    )
  }

  par(mar = c(0, 0, 0, 0), xpd = NA)
  plot.new()
  legend(
    "center",
    legend = design_specs$label,
    col = design_specs$color,
    pch = design_specs$pch,
    lty = design_specs$lty,
    lwd = design_specs$lwd,
    ncol = 3,
    bty = "n",
    cex = 1.05,
    x.intersp = 0.85,
    y.intersp = 1.15
  )

  plot.new()
  text(
    0.5,
    0.82,
    "Scenarios 1-5 are plotted in reverse order: displayed scenarios 1-5 correspond to original scenarios 5-1; scenario 6 is unchanged.",
    cex = 0.83
  )
  text(
    0.5,
    0.52,
    "Scenario 6 note: early stop is treated as correct MTD selection; all non-stop selections are treated as overdosing.",
    cex = 0.83
  )
  text(
    0.5,
    0.22,
    "Patient-overdose percentages use dose-level treated counts above the true MTD divided by total treated counts. Oracle duration is not plotted; oracle sample size is sum of dose-level number treated.",
    cex = 0.83
  )

  mtext(
    paste0(
      "Operating Characteristics Across Six Scenarios (ncycle = ",
      ncycle,
      ", Carryover alpha = ",
      format(alpha_value, trim = TRUE, scientific = FALSE),
      ")"
    ),
    outer = TRUE,
    side = 3,
    line = 0.8,
    cex = 1.55,
    font = 2
  )

  invisible(out_file)
}

out_files <- unlist(
  lapply(alpha_values, function(alpha_value) {
    vapply(cycle_values, function(ncycle) plot_one_cycle(alpha_value, ncycle), character(1))
  }),
  use.names = FALSE
)
cat("Wrote:\n")
cat(paste0("  ", out_files, collapse = "\n"), "\n")
