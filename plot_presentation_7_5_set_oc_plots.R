## ============================================================
## Presentation 7-5: set-specific OC plots for 5-dose scenarios
## ============================================================

target <- 0.30
ncycle_use <- 2L
alpha_values <- c(0, 0.3, 0.6, 0.9)

input_file <- file.path(
  "Presentation 7-5",
  "All_AIDE_OC_Set_5dose_adaptive_r_37_modelsCRM_BOIN_CFO_boinapprox1_approx2_restr_fixed_crmr_fixed_random_alpha_crm_cumu_crm_cfoempirical_pride_dose_summary.csv"
)
boin_adaptive_file <- file.path(
  "Presentation 7-5",
  "All_AIDE_OC_Set_5dose_adaptive_r_37_modelsBOIN_boinapprox1_approx2_restr_adaptive_target0p3_w28_c3_cyc2_rate0p07_Nmax30_45_60_dosecap3_cont1_tried1_ptarget1_0_jobs1to2000_dose_summary.csv"
)

out_dir <- file.path("Presentation 7-5", "Scenarios", "set_specific_OC_ncycle2")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

dat_main <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
dat_boin_adaptive <- read.csv(boin_adaptive_file, stringsAsFactors = FALSE, check.names = FALSE)

missing_from_boin <- setdiff(names(dat_main), names(dat_boin_adaptive))
missing_from_main <- setdiff(names(dat_boin_adaptive), names(dat_main))
if (length(missing_from_boin) > 0L || length(missing_from_main) > 0L) {
  stop(
    "Main and BOIN adaptive summaries have different columns. ",
    "Missing from BOIN adaptive: ", paste(missing_from_boin, collapse = ", "),
    ". Missing from main: ", paste(missing_from_main, collapse = ", ")
  )
}

dat <- rbind(
  dat_main,
  dat_boin_adaptive[, names(dat_main), drop = FALSE]
)

required_cols <- c(
  "Scenario", "Scenario_Group", "True_MTD", "Model", "Method", "CRM_r_model",
  "CFO_Method", "BOIN_Method", "BOIN_r_estimator", "Alpha_true", "r_carry",
  "Nmax_eff", "Cycle_Max", "Restrict_To_Target", "Dose", "True_DLT_rate",
  "MTD_Selection_pct",
  "Pts_Treated", "IPDE_Doses", "Total_Unique_Patients",
  "Early_Stopping_pct", "Duration"
)
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0L) {
  stop("Input summary is missing columns: ", paste(missing_cols, collapse = ", "))
}

scenario_sets <- list(
  baseline = list(
    label = "Baseline Scenarios",
    scenarios = 1:6
  ),
  larger_gaps_below_MTD_s26_31 = list(
    label = "Larger Gaps Below MTD",
    scenarios = 26:31
  ),
  larger_gaps_above_MTD_s32_37 = list(
    label = "Larger Gaps Above MTD",
    scenarios = 32:37
  )
)

design_specs <- data.frame(
  key = c(
    "cfo_pride", "boin_approx1_r_adaptive", "boin_approx2_r_adaptive",
    "crm_alpha", "crm_ip", "crm_fixed", "crm_random", "oracle"
  ),
  label = c(
    "CFO PRIDE",
    "BOIN approx1 r-adaptive",
    "BOIN approx2 r-adaptive",
    "AIDE-alpha-CRM",
    "AIDE-IP-CRM",
    "AIDE-CRM fixed, r = 0",
    "AIDE-CRM random",
    "Oracle CRM"
  ),
  color = c("#9467bd", "#d62728", "#8c564b", "#2ca02c", "#ff7f0e", "#1f77b4", "#e377c2", "#17becf"),
  pch = c(8, 0, 2, 17, 15, 16, 18, 4),
  lty = c(4, 2, 2, 1, 1, 1, 1, 3),
  lwd = c(2, 2, 2, 2, 2, 2, 2, 2),
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

label_crm_design <- function(r_model) {
  out <- rep(NA_character_, length(r_model))
  out[r_model == "r_fixed"] <- "crm_fixed"
  out[r_model == "alpha_crm"] <- "crm_alpha"
  out[r_model == "cumu_crm"] <- "crm_ip"
  out[r_model == "random"] <- "crm_random"
  out
}

make_plot_rows <- function(alpha_value,
                           nmax_eff,
                           restrict_to_target,
                           scenario_ids) {
  base_filter <- dat$Scenario %in% scenario_ids &
    dat$Cycle_Max == ncycle_use &
    dat$Nmax_eff == nmax_eff &
    dat$Restrict_To_Target == as.integer(isTRUE(restrict_to_target))
  
  crm_use <- dat[
    base_filter &
      dat$Model == "CRM" &
      abs(dat$Alpha_true - alpha_value) < 1e-12 &
      dat$CRM_r_model %in% c("r_fixed", "alpha_crm", "cumu_crm", "random") &
      (dat$CRM_r_model != "r_fixed" | abs(dat$r_carry) < 1e-12),
    ,
    drop = FALSE
  ]
  crm_use$Design <- label_crm_design(crm_use$CRM_r_model)
  crm_use$Is_Oracle <- FALSE
  
  cfo_pride <- dat[
    base_filter &
      dat$Model == "CFO" &
      dat$Method == "cfo_pride" &
      dat$CFO_Method == "pride" &
      abs(dat$Alpha_true - alpha_value) < 1e-12,
    ,
    drop = FALSE
  ]
  cfo_pride$Design <- "cfo_pride"
  cfo_pride$Is_Oracle <- FALSE

  boin_approx1_r_adaptive <- dat[
    base_filter &
      dat$Model == "BOIN" &
      dat$Method == "approx1-r_adaptive" &
      dat$BOIN_Method == "approx1" &
      dat$BOIN_r_estimator == "r_adaptive" &
      abs(dat$Alpha_true - alpha_value) < 1e-12,
    ,
    drop = FALSE
  ]
  boin_approx1_r_adaptive$Design <- "boin_approx1_r_adaptive"
  boin_approx1_r_adaptive$Is_Oracle <- FALSE

  boin_approx2_r_adaptive <- dat[
    base_filter &
      dat$Model == "BOIN" &
      dat$Method == "approx2-r_adaptive" &
      dat$BOIN_Method == "approx2" &
      dat$BOIN_r_estimator == "r_adaptive" &
      abs(dat$Alpha_true - alpha_value) < 1e-12,
    ,
    drop = FALSE
  ]
  boin_approx2_r_adaptive$Design <- "boin_approx2_r_adaptive"
  boin_approx2_r_adaptive$Is_Oracle <- FALSE
  
  oracle <- dat[
    base_filter &
      dat$Model == "CRM" &
      dat$CRM_r_model == "r_fixed" &
      abs(dat$r_carry) < 1e-12 &
      abs(dat$Alpha_true) < 1e-12,
    ,
    drop = FALSE
  ]
  oracle$Design <- "oracle"
  oracle$Is_Oracle <- TRUE
  
  common_cols <- Reduce(
    intersect,
    list(
      names(crm_use),
      names(cfo_pride),
      names(boin_approx1_r_adaptive),
      names(boin_approx2_r_adaptive),
      names(oracle)
    )
  )
  combined <- rbind(
    crm_use[common_cols],
    cfo_pride[common_cols],
    boin_approx1_r_adaptive[common_cols],
    boin_approx2_r_adaptive[common_cols],
    oracle[common_cols]
  )
  combined <- combined[!is.na(combined$Design), , drop = FALSE]
  
  expected_designs <- design_specs$key
  found_designs <- sort(unique(combined$Design))
  missing_designs <- setdiff(expected_designs, found_designs)
  if (length(missing_designs) > 0L) {
    warning(
      "Missing designs for alpha=", alpha_value,
      ", Nmax=", nmax_eff,
      ", restrict_to_target=", restrict_to_target,
      ", scenarios=", paste(scenario_ids, collapse = ","),
      ": ", paste(missing_designs, collapse = ", "),
      call. = FALSE
    )
  }
  
  combined
}

summarize_group <- function(d) {
  d <- d[order(d$Dose), , drop = FALSE]
  
  scenario <- first_value(d$Scenario)
  scenario_group <- first_value(d$Scenario_Group)
  true_mtd <- first_value(d$True_MTD)
  design <- first_value(d$Design)
  is_oracle <- isTRUE(first_value(d$Is_Oracle))
  
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
  
  sample_size <- if (is_oracle) total_treated else first_value(d$Total_Unique_Patients)
  duration <- if (is_oracle) NA_real_ else first_value(d$Duration)
  
  data.frame(
    Scenario = scenario,
    Scenario_Group = scenario_group,
    True_MTD = true_mtd,
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

summarize_plot_data <- function(rows, scenario_ids) {
  if (nrow(rows) == 0L) {
    stop("No rows available after filtering.")
  }
  
  split_rows <- split(rows, paste(rows$Design, rows$Scenario, sep = "__"))
  out <- do.call(rbind, lapply(split_rows, summarize_group))
  out$Scenario <- factor(out$Scenario, levels = scenario_ids)
  out$Scenario_Index <- as.integer(out$Scenario)
  out$Design <- factor(out$Design, levels = design_specs$key)
  out <- out[order(out$Design, out$Scenario_Index), ]
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
  
  for (ii in seq_len(nrow(design_specs))) {
    spec <- design_specs[ii, ]
    d <- plot_data[plot_data$Design == spec$key, ]
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
                             alpha_value,
                             nmax_eff,
                             restrict_to_target) {
  scenario_ids <- set_info$scenarios
  plot_rows <- make_plot_rows(
    alpha_value = alpha_value,
    nmax_eff = nmax_eff,
    restrict_to_target = restrict_to_target,
    scenario_ids = scenario_ids
  )
  plot_data <- summarize_plot_data(plot_rows, scenario_ids)
  
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
    "Set_5dose_", set_key, "_OC_",
    alpha_tag(alpha_value),
    "_ncycle", ncycle_use,
    "_Nmax", nmax_eff,
    "_ptarget", as.integer(isTRUE(restrict_to_target))
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
    legend = design_specs$label,
    col = design_specs$color,
    pch = design_specs$pch,
    lty = design_specs$lty,
    lwd = design_specs$lwd,
    ncol = 4,
    bty = "n",
    cex = 0.98,
    x.intersp = 0.8,
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
    "Oracle CRM uses CRM r_fixed with r = 0 and alpha = 0. BOIN uses approx1/approx2 with r_adaptive from the adaptive summary.",
    cex = 0.82
  )
  text(
    0.5,
    0.26,
    "Patient-overdose percentages use dose-level treated counts above target divided by total treated counts. Oracle duration is not plotted.",
    cex = 0.82
  )
  
  mtext(
    paste0(
      set_info$label,
      " (alpha = ",
      format(alpha_value, trim = TRUE, scientific = FALSE),
      ", ncycle = ", ncycle_use,
      ", Nmax = ", nmax_eff,
      ", restrict_to_target = ", as.integer(isTRUE(restrict_to_target)),
      ")"
    ),
    outer = TRUE,
    side = 3,
    line = 0.8,
    cex = 1.45,
    font = 2
  )
  
  invisible(out_file)
}

nmax_values <- sort(unique(dat$Nmax_eff))
restrict_to_target_values <- sort(unique(dat$Restrict_To_Target))

out_files <- character(0)
for (set_key in names(scenario_sets)) {
  for (nmax_eff in nmax_values) {
    for (restrict_to_target in restrict_to_target_values) {
      for (alpha_value in alpha_values) {
        out_files <- c(
          out_files,
          plot_one_setting(
            set_key = set_key,
            set_info = scenario_sets[[set_key]],
            alpha_value = alpha_value,
            nmax_eff = nmax_eff,
            restrict_to_target = as.logical(restrict_to_target)
          )
        )
      }
    }
  }
}

cat("Wrote:\n")
cat(paste0("  ", out_files, collapse = "\n"), "\n")
