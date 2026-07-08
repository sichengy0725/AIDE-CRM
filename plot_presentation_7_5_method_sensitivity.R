## ============================================================
## Presentation 7-5: method-specific sensitivity OC plots
## ============================================================

target <- 0.30
ncycle_use <- 2L
alpha_values <- c(0, 0.3, 0.6, 0.9)
nmax_values <- c(30L, 45L, 60L)
restrict_to_target_values <- c(0L, 1L)

main_file <- file.path(
  "Presentation 7-5",
  "All_AIDE_OC_Set_5dose_adaptive_r_37_modelsCRM_BOIN_CFO_boinapprox1_approx2_restr_fixed_crmr_fixed_random_alpha_crm_cumu_crm_cfoempirical_pride_dose_summary.csv"
)
boin_adaptive_file <- file.path(
  "Presentation 7-5",
  "All_AIDE_OC_Set_5dose_adaptive_r_37_modelsBOIN_boinapprox1_approx2_restr_adaptive_target0p3_w28_c3_cyc2_rate0p07_Nmax30_45_60_dosecap3_cont1_tried1_ptarget1_0_jobs1to2000_dose_summary.csv"
)

out_root <- file.path("Presentation 7-5", "Scenarios", "method_sensitivity_OC_ncycle2")
out_dir_nmax <- file.path(out_root, "compare_Nmax")
out_dir_target <- file.path(out_root, "compare_restrict_to_target")
dir.create(out_dir_nmax, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_target, recursive = TRUE, showWarnings = FALSE)

dat_main <- read.csv(main_file, stringsAsFactors = FALSE, check.names = FALSE)
dat_boin <- read.csv(boin_adaptive_file, stringsAsFactors = FALSE, check.names = FALSE)

missing_from_boin <- setdiff(names(dat_main), names(dat_boin))
missing_from_main <- setdiff(names(dat_boin), names(dat_main))
if (length(missing_from_boin) > 0L || length(missing_from_main) > 0L) {
  stop(
    "Main and BOIN adaptive summaries have different columns. ",
    "Missing from BOIN adaptive: ", paste(missing_from_boin, collapse = ", "),
    ". Missing from main: ", paste(missing_from_main, collapse = ", ")
  )
}

dat <- rbind(dat_main, dat_boin[, names(dat_main), drop = FALSE])

required_cols <- c(
  "Scenario", "Scenario_Group", "True_MTD", "Model", "Method", "CRM_r_model",
  "BOIN_Method", "BOIN_r_estimator", "Alpha_true", "Nmax_eff", "Cycle_Max",
  "Restrict_To_Target", "Dose", "True_DLT_rate", "MTD_Selection_pct",
  "Pts_Treated", "IPDE_Doses", "Total_Unique_Patients",
  "Early_Stopping_pct", "Duration"
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
  key = c("crm_random", "boin_approx1_r_adaptive", "boin_approx2_r_adaptive"),
  label = c(
    "CRM random",
    "BOIN approx1 r-adaptive",
    "BOIN approx2 r-adaptive"
  ),
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

target_specs <- data.frame(
  key = paste0("ptarget_", restrict_to_target_values),
  label = paste0("restrict_to_target = ", restrict_to_target_values),
  value = restrict_to_target_values,
  color = c("#1f77b4", "#d62728"),
  pch = c(16, 0),
  lty = c(1, 2),
  lwd = c(2, 2),
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

method_filter <- function(method_key) {
  if (identical(method_key, "crm_random")) {
    return(dat$Model == "CRM" & dat$CRM_r_model == "random")
  }
  if (identical(method_key, "boin_approx1_r_adaptive")) {
    return(
      dat$Model == "BOIN" &
        dat$Method == "approx1-r_adaptive" &
        dat$BOIN_Method == "approx1" &
        dat$BOIN_r_estimator == "r_adaptive"
    )
  }
  if (identical(method_key, "boin_approx2_r_adaptive")) {
    return(
      dat$Model == "BOIN" &
        dat$Method == "approx2-r_adaptive" &
        dat$BOIN_Method == "approx2" &
        dat$BOIN_r_estimator == "r_adaptive"
    )
  }
  stop("Unknown method_key: ", method_key)
}

method_label <- function(method_key) {
  out <- method_specs$label[method_specs$key == method_key]
  if (length(out) != 1L) {
    stop("Unknown method_key: ", method_key)
  }
  out
}

make_nmax_rows <- function(method_key,
                           alpha_value,
                           restrict_to_target,
                           scenario_ids) {
  base_filter <- dat$Scenario %in% scenario_ids &
    dat$Cycle_Max == ncycle_use &
    dat$Restrict_To_Target == restrict_to_target &
    abs(dat$Alpha_true - alpha_value) < 1e-12 &
    method_filter(method_key)

  out <- lapply(seq_len(nrow(nmax_specs)), function(ii) {
    spec <- nmax_specs[ii, ]
    rows <- dat[base_filter & dat$Nmax_eff == spec$value, , drop = FALSE]
    rows$Line_Key <- spec$key
    rows
  })

  out <- do.call(rbind, out)
  out
}

make_target_rows <- function(method_key,
                             alpha_value,
                             nmax_eff,
                             scenario_ids) {
  base_filter <- dat$Scenario %in% scenario_ids &
    dat$Cycle_Max == ncycle_use &
    dat$Nmax_eff == nmax_eff &
    abs(dat$Alpha_true - alpha_value) < 1e-12 &
    method_filter(method_key)

  out <- lapply(seq_len(nrow(target_specs)), function(ii) {
    spec <- target_specs[ii, ]
    rows <- dat[base_filter & dat$Restrict_To_Target == spec$value, , drop = FALSE]
    rows$Line_Key <- spec$key
    rows
  })

  out <- do.call(rbind, out)
  out
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

summarize_plot_data <- function(rows, scenario_ids, line_specs) {
  if (nrow(rows) == 0L) {
    stop("No rows available after filtering.")
  }

  split_rows <- split(rows, paste(rows$Line_Key, rows$Scenario, sep = "__"))
  out <- do.call(rbind, lapply(split_rows, summarize_group))
  out$Scenario <- factor(out$Scenario, levels = scenario_ids)
  out$Scenario_Index <- as.integer(out$Scenario)
  out$Line_Key <- factor(out$Line_Key, levels = line_specs$key)
  out <- out[order(out$Line_Key, out$Scenario_Index), ]
  row.names(out) <- NULL

  expected <- line_specs$key
  found <- sort(unique(as.character(out$Line_Key)))
  missing <- setdiff(expected, found)
  if (length(missing) > 0L) {
    warning("Missing comparison lines: ", paste(missing, collapse = ", "), call. = FALSE)
  }

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

plot_metric <- function(plot_data, metric, panel_title, ylab, ylim, scenario_ids, line_specs) {
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

  for (ii in seq_len(nrow(line_specs))) {
    spec <- line_specs[ii, ]
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

plot_with_lines <- function(plot_data,
                            line_specs,
                            scenario_ids,
                            title,
                            footer_line,
                            out_file) {
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
      scenario_ids = scenario_ids,
      line_specs = line_specs
    )
  }

  par(mar = c(0, 0, 0, 0), xpd = NA)
  plot.new()
  legend(
    "center",
    legend = line_specs$label,
    col = line_specs$color,
    pch = line_specs$pch,
    lty = line_specs$lty,
    lwd = line_specs$lwd,
    ncol = nrow(line_specs),
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
    footer_line,
    cex = 0.82
  )
  text(
    0.5,
    0.26,
    "Patient-overdose percentages use dose-level treated counts above target divided by total treated counts.",
    cex = 0.82
  )

  mtext(title, outer = TRUE, side = 3, line = 0.8, cex = 1.42, font = 2)

  invisible(out_file)
}

plot_nmax_setting <- function(set_key,
                              set_info,
                              method_key,
                              alpha_value,
                              restrict_to_target) {
  scenario_ids <- set_info$scenarios
  rows <- make_nmax_rows(
    method_key = method_key,
    alpha_value = alpha_value,
    restrict_to_target = restrict_to_target,
    scenario_ids = scenario_ids
  )
  plot_data <- summarize_plot_data(rows, scenario_ids, nmax_specs)

  file_stub <- paste0(
    "Nmax_5dose_", set_key, "_", method_key, "_OC_",
    alpha_tag(alpha_value),
    "_ncycle", ncycle_use,
    "_ptarget", restrict_to_target
  )
  out_file <- file.path(out_dir_nmax, paste0(file_stub, ".pdf"))

  title <- paste0(
    set_info$label,
    " - ", method_label(method_key),
    " Nmax Comparison",
    " (alpha = ", format(alpha_value, trim = TRUE, scientific = FALSE),
    ", ncycle = ", ncycle_use,
    ", restrict_to_target = ", restrict_to_target,
    ")"
  )
  footer_line <- "Line colors compare Nmax values within the same method, alpha, ncycle, and target gate."

  plot_with_lines(plot_data, nmax_specs, scenario_ids, title, footer_line, out_file)
}

plot_target_setting <- function(set_key,
                                set_info,
                                method_key,
                                alpha_value,
                                nmax_eff) {
  scenario_ids <- set_info$scenarios
  rows <- make_target_rows(
    method_key = method_key,
    alpha_value = alpha_value,
    nmax_eff = nmax_eff,
    scenario_ids = scenario_ids
  )
  plot_data <- summarize_plot_data(rows, scenario_ids, target_specs)

  file_stub <- paste0(
    "RestrictTarget_5dose_", set_key, "_", method_key, "_OC_",
    alpha_tag(alpha_value),
    "_ncycle", ncycle_use,
    "_Nmax", nmax_eff
  )
  out_file <- file.path(out_dir_target, paste0(file_stub, ".pdf"))

  title <- paste0(
    set_info$label,
    " - ", method_label(method_key),
    " Target Gate Comparison",
    " (alpha = ", format(alpha_value, trim = TRUE, scientific = FALSE),
    ", ncycle = ", ncycle_use,
    ", Nmax = ", nmax_eff,
    ")"
  )
  footer_line <- "Line colors compare restrict_to_target values within the same method, alpha, ncycle, and Nmax."

  plot_with_lines(plot_data, target_specs, scenario_ids, title, footer_line, out_file)
}

out_files <- character(0)

for (set_key in names(scenario_sets)) {
  for (method_key in method_specs$key) {
    for (restrict_to_target in restrict_to_target_values) {
      for (alpha_value in alpha_values) {
        out_files <- c(
          out_files,
          plot_nmax_setting(
            set_key = set_key,
            set_info = scenario_sets[[set_key]],
            method_key = method_key,
            alpha_value = alpha_value,
            restrict_to_target = restrict_to_target
          )
        )
      }
    }

    for (nmax_eff in nmax_values) {
      for (alpha_value in alpha_values) {
        out_files <- c(
          out_files,
          plot_target_setting(
            set_key = set_key,
            set_info = scenario_sets[[set_key]],
            method_key = method_key,
            alpha_value = alpha_value,
            nmax_eff = nmax_eff
          )
        )
      }
    }
  }
}

cat("Wrote:\n")
cat(paste0("  ", out_files, collapse = "\n"), "\n")
