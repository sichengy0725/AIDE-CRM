## ============================================================
## Presentation 7-6-2026: alpha_crm and cumu_crm Nmax plots
## ============================================================

target <- 0.30
ncycle_use <- 2L
alpha_values <- c(0, 0.3, 0.6, 0.9)
nmax_values <- c(30L, 45L, 60L)
restrict_to_target_use <- 0L

input_file <- file.path(
  "Presentation 7-5",
  "All_AIDE_OC_Set_5dose_adaptive_r_37_modelsCRM_BOIN_CFO_boinapprox1_approx2_restr_fixed_crmr_fixed_random_alpha_crm_cumu_crm_cfoempirical_pride_dose_summary.csv"
)

out_dir <- file.path(
  "Presentation 7-6-2026",
  "Plots",
  "Sensitivity",
  "Sample Size"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c(
  "Scenario", "Scenario_Group", "True_MTD", "Model", "CRM_r_model",
  "Alpha_true", "Nmax_eff", "Cycle_Max", "Restrict_To_Target",
  "Dose", "True_DLT_rate", "MTD_Selection_pct", "Pts_Treated",
  "IPDE_Doses", "Total_Unique_Patients", "Early_Stopping_pct", "Duration"
)
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0L) {
  stop("Input summary is missing columns: ", paste(missing_cols, collapse = ", "))
}

scenario_sets <- list(
  baseline_scenarios = list(
    label = "Baseline Scenarios",
    scenarios = 1:6
  ),
  larger_gap_below_MTD_s26_31 = list(
    label = "Larger Gap Below MTD",
    scenarios = 26:31
  ),
  larger_gap_above_MTD_s32_37 = list(
    label = "Larger Gap Above MTD",
    scenarios = 32:37
  )
)

method_specs <- data.frame(
  key = c("alpha_crm", "cumu_crm"),
  label = c("AIDE-alpha-CRM", "AIDE-IP-CRM"),
  stringsAsFactors = FALSE
)

nmax_specs <- data.frame(
  key = paste0("nmax_", nmax_values),
  label = paste0("Nmax = ", nmax_values),
  value = nmax_values,
  color = c("#1f77b4", "#ff7f0e", "#2ca02c"),
  pch = c(16, 15, 17),
  lty = c(1, 1, 1),
  lwd = c(2, 2, 2),
  stringsAsFactors = FALSE
)

alpha_tag <- function(alpha_value) {
  tag <- format(alpha_value, trim = TRUE, scientific = FALSE)
  paste0("alpha", gsub("\\.", "p", tag))
}

first_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  x[[1]]
}

method_label <- function(method_key) {
  out <- method_specs$label[method_specs$key == method_key]
  if (length(out) != 1L) {
    stop("Unknown method_key: ", method_key)
  }
  out
}

make_rows <- function(method_key,
                      alpha_value,
                      scenario_ids) {
  base_filter <- dat$Scenario %in% scenario_ids &
    dat$Cycle_Max == ncycle_use &
    dat$Restrict_To_Target == restrict_to_target_use &
    abs(dat$Alpha_true - alpha_value) < 1e-12 &
    dat$Model == "CRM" &
    dat$CRM_r_model == method_key

  out <- lapply(seq_len(nrow(nmax_specs)), function(ii) {
    spec <- nmax_specs[ii, ]
    rows <- dat[base_filter & dat$Nmax_eff == spec$value, , drop = FALSE]
    rows$Line_Key <- spec$key
    rows
  })

  do.call(rbind, out)
}

summarize_group <- function(d) {
  d <- d[order(d$Dose), , drop = FALSE]

  scenario <- first_value(d$Scenario)
  scenario_group <- first_value(d$Scenario_Group)
  true_mtd <- first_value(d$True_MTD)
  line_key <- first_value(d$Line_Key)

  true_rates <- d$True_DLT_rate
  total_treated <- sum(d$Pts_Treated, na.rm = TRUE)
  if (!is.finite(total_treated) || total_treated <= 0) {
    total_treated <- NA_real_
  }

  selected_toxic <- sum(d$MTD_Selection_pct[true_rates > target + 1e-12], na.rm = TRUE)
  treated_toxic <- sum(d$Pts_Treated[true_rates > target + 1e-12], na.rm = TRUE)
  treated_subtherapeutic <- sum(d$Pts_Treated[true_rates < target - 1e-12], na.rm = TRUE)

  if (is.na(true_mtd) || true_mtd < 1L || true_mtd > nrow(d)) {
    correct_selection <- first_value(d$Early_Stopping_pct)
    treated_mtd <- 0
  } else {
    correct_selection <- d$MTD_Selection_pct[d$Dose == true_mtd][1]
    treated_mtd <- d$Pts_Treated[d$Dose == true_mtd][1]
  }

  data.frame(
    Scenario = scenario,
    Scenario_Group = scenario_group,
    True_MTD = true_mtd,
    Line_Key = line_key,
    Correct_MTD_Selection = correct_selection,
    Overdose_Selection = selected_toxic,
    Patients_Allocated_MTD = 100 * treated_mtd / total_treated,
    Patients_Overdosed = 100 * treated_toxic / total_treated,
    Patients_Underdosed = 100 * treated_subtherapeutic / total_treated,
    Total_IPDE = sum(d$IPDE_Doses, na.rm = TRUE),
    Sample_Size = first_value(d$Total_Unique_Patients),
    Duration = first_value(d$Duration),
    stringsAsFactors = FALSE
  )
}

summarize_plot_data <- function(rows, scenario_ids) {
  if (nrow(rows) == 0L) {
    stop("No rows available after filtering.")
  }

  split_rows <- split(rows, paste(rows$Line_Key, rows$Scenario, sep = "__"))
  out <- do.call(rbind, lapply(split_rows, summarize_group))
  out$Scenario <- factor(out$Scenario, levels = scenario_ids)
  out$Scenario_Index <- as.integer(out$Scenario)
  out$Line_Key <- factor(out$Line_Key, levels = nmax_specs$key)
  out <- out[order(out$Line_Key, out$Scenario_Index), ]
  row.names(out) <- NULL
  out
}

metric_limits <- function(plot_data, metric) {
  vals <- plot_data[[metric]]
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) {
    return(c(0, 1))
  }

  if (metric %in% c("Correct_MTD_Selection", "Patients_Overdosed", "Patients_Underdosed")) {
    return(c(0, 100))
  }
  if (metric == "Overdose_Selection") {
    return(c(0, max(70, ceiling(max(vals) / 10) * 10)))
  }
  if (metric == "Patients_Allocated_MTD") {
    return(c(0, max(80, ceiling(max(vals) / 10) * 10)))
  }
  if (metric == "Total_IPDE") {
    return(c(0, max(10, ceiling(max(vals) + 1))))
  }

  pad <- max(1, diff(range(vals)) * 0.08)
  c(max(0, floor(min(vals) - pad)), ceiling(max(vals) + pad))
}

plot_metric <- function(plot_data, metric, panel_title, ylab, ylim, scenario_ids) {
  xvals <- seq_along(scenario_ids)

  plot(
    range(xvals),
    ylim,
    type = "n",
    axes = FALSE,
    xlab = "Scenario",
    ylab = ylab,
    main = panel_title,
    font.main = 2,
    cex.main = 1.12,
    cex.lab = 1.02
  )
  axis(1, at = xvals, labels = scenario_ids)
  axis(2)
  box()
  abline(v = xvals, col = "#d9d9d9", lwd = 0.8)
  abline(h = pretty(ylim), col = "#e6e6e6", lwd = 0.8)

  for (ii in seq_len(nrow(nmax_specs))) {
    spec <- nmax_specs[ii, ]
    d <- plot_data[plot_data$Line_Key == spec$key, ]
    d <- d[order(d$Scenario_Index), ]
    y <- d[[metric]]
    x <- d$Scenario_Index
    keep <- is.finite(y)
    if (sum(keep) == 0L) {
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

plot_one_setting <- function(set_key,
                             set_info,
                             method_key,
                             alpha_value) {
  scenario_ids <- set_info$scenarios
  rows <- make_rows(method_key, alpha_value, scenario_ids)
  plot_data <- summarize_plot_data(rows, scenario_ids)

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

  file_stub <- paste0(
    "Nmax_5dose_", set_key, "_", method_key, "_OC_",
    alpha_tag(alpha_value),
    "_ncycle", ncycle_use,
    "_ptarget", restrict_to_target_use
  )
  out_file <- file.path(out_dir, paste0(file_stub, ".pdf"))

  pdf(out_file, width = 13.5, height = 10.5, onefile = FALSE, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)

  layout(
    rbind(c(1, 2), c(3, 4), c(5, 6), c(7, 8), c(9, 9), c(10, 10)),
    heights = c(1, 1, 1, 1, 0.34, 0.24)
  )
  par(oma = c(0, 0, 3.0, 0), mar = c(3.4, 4.2, 2.1, 1.4), mgp = c(2.2, 0.7, 0))

  for (metric in metrics) {
    metric_name <- metric[[1]]
    plot_metric(
      plot_data = plot_data,
      metric = metric_name,
      panel_title = metric[[2]],
      ylab = metric[[3]],
      ylim = metric_limits(plot_data, metric_name),
      scenario_ids = scenario_ids
    )
  }

  par(mar = c(0, 0, 0, 0), xpd = NA)
  plot.new()
  legend(
    "center",
    legend = nmax_specs$label,
    col = nmax_specs$color,
    pch = nmax_specs$pch,
    lty = nmax_specs$lty,
    lwd = nmax_specs$lwd,
    ncol = nrow(nmax_specs),
    bty = "n",
    cex = 1.02,
    x.intersp = 0.85,
    y.intersp = 1.15
  )

  plot.new()
  text(
    0.5,
    0.82,
    paste0(
      "Scenarios shown: ",
      paste(scenario_ids, collapse = ", "),
      ". True MTD labels come from the dose-summary CSV."
    ),
    cex = 0.82
  )
  text(
    0.5,
    0.54,
    "Line colors compare Nmax values within the same method, alpha, ncycle, and target gate.",
    cex = 0.82
  )
  text(
    0.5,
    0.26,
    "Patient-overdose percentages use dose-level treated counts above target divided by total treated counts.",
    cex = 0.82
  )

  mtext(
    paste0(
      set_info$label,
      " - ", method_label(method_key),
      " Nmax Comparison",
      " (alpha = ", format(alpha_value, trim = TRUE, scientific = FALSE),
      ", ncycle = ", ncycle_use,
      ", restrict_to_target = ", restrict_to_target_use,
      ")"
    ),
    outer = TRUE,
    side = 3,
    line = 0.8,
    cex = 1.42,
    font = 2
  )

  invisible(out_file)
}

out_files <- character(0)
for (set_key in names(scenario_sets)) {
  for (method_key in method_specs$key) {
    for (alpha_value in alpha_values) {
      out_files <- c(
        out_files,
        plot_one_setting(
          set_key = set_key,
          set_info = scenario_sets[[set_key]],
          method_key = method_key,
          alpha_value = alpha_value
        )
      )
    }
  }
}

cat("Wrote:\n")
cat(paste0("  ", out_files, collapse = "\n"), "\n")
